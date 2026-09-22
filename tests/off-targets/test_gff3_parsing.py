"""Tests for GFF3 support in the annotation parser.

README advertised "GTF/GFF3 or BED", but the attribute regex only matched GTF's
`key "value";` syntax. GFF3's `key=value;` yielded null gene_id/transcript_id and
exp_value=0.0 with no error — and that zero flows into
``W_i = Expression_i * exp(-dG/RT)``, so a GFF3 run produced silently wrong
off-target weights rather than failing.

Two GFF3 dialects matter in practice and are covered here:
  * GENCODE  — carries literal ``gene_id=`` / ``transcript_id=`` attributes
  * Ensembl  — carries ``Parent=transcript:ENST...`` and needs prefix stripping
"""

from pathlib import Path

import polars as pl
import pytest

from sioff.services.annotation_parser import AnnotationParser


@pytest.fixture
def parser() -> AnnotationParser:
    return AnnotationParser()


def _write(tmp_path: Path, name: str, body: str) -> Path:
    path = tmp_path / name
    path.write_text(body)
    return path


GENCODE_GFF3 = """\
##gff-version 3
chr1\tHAVANA\texon\t1200\t1350\t.\t+\t.\tID=exon:ENST00000456328.2:1;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;RPKM=1000
chr1\tHAVANA\texon\t1400\t1550\t.\t+\t.\tID=exon:ENST00000456328.2:2;Parent=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;RPKM=1000
chr1\tHAVANA\ttranscript\t1000\t11000\t.\t+\t.\tID=ENST00000456328.2;gene_id=ENSG00000290825.1;transcript_id=ENST00000456328.2;RPKM=1000
chr2\tHAVANA\texon\t2200\t2350\t.\t-\t.\tID=exon:ENST00000999999.1:1;Parent=ENST00000999999.1;gene_id=ENSG00000111111.1;transcript_id=ENST00000999999.1;RPKM=25.5
"""

ENSEMBL_GFF3 = """\
##gff-version 3
chr1\tensembl\texon\t1200\t1350\t.\t+\t.\tParent=transcript:ENST00000456328;Name=ENSE00002234944;exon_id=ENSE00002234944;RPKM=750
chr1\tensembl\texon\t1400\t1550\t.\t+\t.\tParent=transcript:ENST00000456328;Name=ENSE00003582793;exon_id=ENSE00003582793;RPKM=750
"""


class TestGff3Gencode:
    def test_extracts_literal_gene_and_transcript_ids(self, parser, tmp_path):
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        df = parser.load_gtf(path)

        row = df.row(0, named=True)
        assert row["gene_id"] == "ENSG00000290825.1"
        assert row["transcript_id"] == "ENST00000456328.2"

    def test_extracts_the_score_attribute(self, parser, tmp_path):
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        df = parser.load_gtf(path)

        assert df.row(0, named=True)["exp_value"] == 1000.0
        assert (
            df.filter(pl.col("chrom") == "chr2").row(0, named=True)["exp_value"] == 25.5
        )

    def test_filters_by_feature_like_gtf(self, parser, tmp_path):
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        df = parser.load_gtf(path, feature="exon")

        assert df.height == 3  # 3 exons, the transcript row excluded

    def test_coordinates_and_strand_survive(self, parser, tmp_path):
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        row = parser.load_gtf(path).row(0, named=True)

        assert (row["chrom"], row["start"], row["end"], row["strand"]) == (
            "chr1",
            1200,
            1350,
            "+",
        )


class TestGff3Ensembl:
    def test_falls_back_to_parent_with_type_prefix_stripped(self, parser, tmp_path):
        path = _write(tmp_path, "ensembl.gff3", ENSEMBL_GFF3)

        df = parser.load_gtf(path)

        assert df.row(0, named=True)["transcript_id"] == "ENST00000456328"

    def test_gene_id_falls_back_to_the_transcript_when_absent(self, parser, tmp_path):
        """Ensembl exon lines carry no gene_id; the row must still be usable.

        Falling back keeps every row keyed by something real rather than dropping
        it, which matters because these rows carry the expression values.
        """
        path = _write(tmp_path, "ensembl.gff3", ENSEMBL_GFF3)

        df = parser.load_gtf(path)

        assert df.row(0, named=True)["gene_id"] == "ENST00000456328"


class TestGff3Detection:
    def test_detected_from_attribute_syntax_not_just_the_extension(
        self, parser, tmp_path
    ):
        """A misnamed .gtf file holding GFF3 attributes must still parse.

        `format: auto` is the default, and annotation files are routinely
        misnamed, so sniffing the attribute column is what makes auto honest.
        """
        path = _write(tmp_path, "actually_gff3.gtf", GENCODE_GFF3)

        df = parser.load_gtf(path)

        assert df.row(0, named=True)["gene_id"] == "ENSG00000290825.1"

    def test_explicit_gff3_format_is_accepted(self, parser, tmp_path):
        path = _write(tmp_path, "ann.gff3", GENCODE_GFF3)

        df = parser.load_gtf(path, format="gff3")

        assert df.row(0, named=True)["transcript_id"] == "ENST00000456328.2"

    def test_gtf_still_parses_unchanged(self, parser):
        """Regression guard: GFF3 support must not disturb the GTF path."""
        gtf = Path(__file__).parent / "data" / "expression_data.gtf"

        row = parser.load_gtf(gtf).row(0, named=True)

        assert row["gene_id"] == "gene_1"
        assert row["transcript_id"] == "transcript_21"
        assert row["exp_value"] == 1000.0


class TestGff3UrlDecoding:
    def test_percent_escapes_are_decoded(self, parser, tmp_path):
        """The GFF3 spec requires ; = % , and tab to be URL-escaped in values."""
        body = (
            "##gff-version 3\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\t"
            "gene_id=gene%3Bwith%3Dchars;transcript_id=tx1;RPKM=5\n"
        )
        path = _write(tmp_path, "escaped.gff3", body)

        df = parser.load_gtf(path)

        assert df.row(0, named=True)["gene_id"] == "gene;with=chars"


class TestExtractionFailureIsLoud:
    def test_raises_when_no_row_yields_an_id(self, parser, tmp_path):
        body = (
            "##gff-version 3\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tsomething_else=x;RPKM=5\n"
        )
        path = _write(tmp_path, "noids.gff3", body)

        with pytest.raises(ValueError, match="transcript_id"):
            parser.load_gtf(path)

    def test_raises_when_the_score_column_matches_nothing(self, parser, tmp_path):
        """The original silent-wrong-answer path: exp_value=0.0 for every row."""
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        with pytest.raises(ValueError, match="NONEXISTENT"):
            parser.load_gtf(path, score_col="NONEXISTENT")

    def test_the_error_names_the_file_and_the_detected_format(self, parser, tmp_path):
        path = _write(tmp_path, "gencode.gff3", GENCODE_GFF3)

        with pytest.raises(ValueError) as excinfo:
            parser.load_gtf(path, score_col="NONEXISTENT")

        message = str(excinfo.value)
        assert "gencode.gff3" in message
        assert "gff3" in message

    def test_a_partially_populated_score_column_does_not_raise(self, parser, tmp_path):
        """Only *total* failure is an error; sparse annotation stays allowed."""
        body = (
            "##gff-version 3\n"
            "chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id=g1;transcript_id=t1;RPKM=5\n"
            "chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id=g2;transcript_id=t2\n"
        )
        path = _write(tmp_path, "sparse.gff3", body)

        df = parser.load_gtf(path)

        assert df["exp_value"].to_list() == [5.0, 0.0]
