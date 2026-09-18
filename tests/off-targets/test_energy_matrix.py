"""Energy-parameter-set (DSM) selection must cover everything risearch ships.

risearch bundles five nearest-neighbour scoring models, verified against the
installed bindings:

    t04          Turner 2004            (RNA-RNA, default)
    t99          Turner 1999            (RNA-RNA)
    slh04        SantaLucia-Hicks 2004  (DNA-DNA)
    s95-rna-dna  Sugimoto 1995          (RNA query / DNA target)
    s95-dna-rna  Sugimoto 1995          (DNA query / RNA target)

RIOT previously hard-rejected everything but t04 and t99, which put the two
Sugimoto hybrid tables out of reach — the ones that apply to DNA-based oligos.
The choice is not cosmetic: on the shipped fixtures t04 and t99 return different
hit counts, so it changes which off-targets are found and their energies.
"""

import pytest

from riot.services.risearch_service import VALID_DSM_IDS, RIsearchError, RIsearchService


class TestValidDsmIds:
    def test_every_model_risearch_ships_is_allowed(self):
        assert VALID_DSM_IDS == frozenset(
            {"t04", "t99", "slh04", "s95-rna-dna", "s95-dna-rna"}
        )

    @pytest.mark.parametrize("matrix", sorted(VALID_DSM_IDS))
    def test_each_one_passes_validation(self, tmp_path, matrix):
        """Validation runs before risearch is imported, so this needs no extra."""
        index = tmp_path / "g.idx"
        index.touch()
        (tmp_path / "q.fa").write_text(">q\nACGU\n")

        service = RIsearchService()
        service._target_registry[str(index)] = tmp_path / "g.fa"

        # Reaching the import means validation accepted the matrix. Without the
        # optional dependency installed that surfaces as a RIsearchError naming
        # the missing module, never one complaining about the matrix.
        try:
            service.run_search(
                query_path=tmp_path / "q.fa",
                index_path=index,
                target_fasta=tmp_path / "g.fa",
                matrix=matrix,
            )
        except RIsearchError as exc:
            assert "matrix" not in str(exc), f"{matrix} was rejected: {exc}"
        except Exception:
            pass  # any other failure is downstream of validation

    def test_an_unknown_id_is_rejected_with_the_valid_set_named(self, tmp_path):
        index = tmp_path / "g.idx"
        index.touch()
        (tmp_path / "q.fa").write_text(">q\nACGU\n")
        service = RIsearchService()
        service._target_registry[str(index)] = tmp_path / "g.fa"

        with pytest.raises(RIsearchError) as excinfo:
            service.run_search(
                query_path=tmp_path / "q.fa",
                index_path=index,
                target_fasta=tmp_path / "g.fa",
                matrix="bogus",
            )

        message = str(excinfo.value)
        assert "bogus" in message
        for valid in VALID_DSM_IDS:
            assert valid in message, f"error should list {valid}"

    def test_the_bare_directory_name_s95_is_not_a_valid_id(self, tmp_path):
        """`s95` is the data directory; the ids are s95-rna-dna / s95-dna-rna.

        risearch itself rejects it with "unknown DSM id 's95'", so accepting it
        here would only defer the failure.
        """
        index = tmp_path / "g.idx"
        index.touch()
        (tmp_path / "q.fa").write_text(">q\nACGU\n")
        service = RIsearchService()
        service._target_registry[str(index)] = tmp_path / "g.fa"

        with pytest.raises(RIsearchError, match="s95"):
            service.run_search(
                query_path=tmp_path / "q.fa",
                index_path=index,
                target_fasta=tmp_path / "g.fa",
                matrix="s95",
            )


class TestPublicApiExposesTheSearchKnobs:
    """`riot.search` must offer what the CLI offers.

    `--matrix`/`-z` and the seed-geometry flags existed on the CLI while the
    Python API silently pinned every library user to the defaults.
    """

    @pytest.mark.parametrize(
        "parameter", ["matrix", "seed_start", "seed_end", "seed_wobble"]
    )
    def test_parameter_is_present(self, parameter):
        import inspect

        import riot

        assert parameter in inspect.signature(riot.search).parameters

    def test_defaults_match_the_core_layer(self):
        import inspect

        import riot
        from riot.core.risearch import run_search as core_search

        api = inspect.signature(riot.search).parameters
        core = inspect.signature(core_search).parameters

        for name in ("matrix", "seed_start", "seed_end", "seed_wobble"):
            assert api[name].default == core[name].default, (
                f"{name} default drifted between the API and core layers"
            )
