"""
Tests for IntersectionService overlap semantics.

Focuses on the specific behaviours that were fixed:
  1. Any-overlap (not containment-only)
  2. Boundary-touching intervals do NOT overlap (BED half-open)
"""

import polars as pl

from riot.services.intersection_service import IntersectionService


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def pred(start, end, chrom="chr1", strand="+", sid="si1", energy=-15.0):
    return pl.DataFrame(
        {
            "sirna_id": [sid],
            "chrom": [chrom],
            "start": [start],
            "end": [end],
            "strand": [strand],
            "energy": [energy],
        }
    )


def trans(start, end, chrom="chr1", strand="+", gene="g1", tid="t1", exp=10.0):
    return pl.DataFrame(
        {
            "chrom": [chrom],
            "start": [start],
            "end": [end],
            "strand": [strand],
            "gene_id": [gene],
            "transcript_id": [tid],
            "exp_value": [exp],
        }
    )


def hits(svc, p, t):
    return svc.intersect(p, t).height


# ---------------------------------------------------------------------------
# Any-overlap tests
# ---------------------------------------------------------------------------


class TestAnyOverlap:
    """Partial overlap in either direction must be detected."""

    def test_pred_overlaps_right_end_of_transcript(self):
        # pred [10,20], trans [5,15]: pred_start=10 < trans_end=15 ✓  pred_end=20 > trans_start=5 ✓
        svc = IntersectionService()
        assert hits(svc, pred(10, 20), trans(5, 15)) == 1

    def test_pred_overlaps_left_end_of_transcript(self):
        # pred [10,20], trans [15,30]: pred_start=10 < trans_end=30 ✓  pred_end=20 > trans_start=15 ✓
        svc = IntersectionService()
        assert hits(svc, pred(10, 20), trans(15, 30)) == 1

    def test_pred_contained_inside_transcript(self):
        svc = IntersectionService()
        assert hits(svc, pred(10, 20), trans(5, 30)) == 1

    def test_transcript_contained_inside_pred(self):
        svc = IntersectionService()
        assert hits(svc, pred(5, 30), trans(10, 20)) == 1

    def test_exact_same_coordinates(self):
        svc = IntersectionService()
        assert hits(svc, pred(10, 20), trans(10, 20)) == 1


# ---------------------------------------------------------------------------
# Boundary-touching tests (the NCLS +1 fix)
# ---------------------------------------------------------------------------


class TestBoundaryTouching:
    """
    BED intervals are half-open [start, end).
    pred_end == trans_start means the intervals share one coordinate but do NOT overlap.
    These must return zero hits.
    """

    def test_pred_end_equals_trans_start_no_hit(self):
        # pred [10, 20], trans [20, 30]: pred_end = trans_start = 20 → no overlap
        svc = IntersectionService()
        assert hits(svc, pred(10, 20), trans(20, 30)) == 0

    def test_pred_start_equals_trans_end_no_hit(self):
        # pred [20, 30], trans [10, 20]: pred_start = trans_end = 20 → no overlap
        svc = IntersectionService()
        assert hits(svc, pred(20, 30), trans(10, 20)) == 0

    def test_adjacent_intervals_no_hit(self):
        # pred [0, 5], trans [5, 10] → touching, no overlap
        svc = IntersectionService()
        assert hits(svc, pred(0, 5), trans(5, 10)) == 0

    def test_one_base_overlap_is_a_hit(self):
        # pred [9, 11], trans [10, 20]: pred_start=9 < trans_end=20, pred_end=11 > trans_start=10 → overlap
        svc = IntersectionService()
        assert hits(svc, pred(9, 11), trans(10, 20)) == 1


# ---------------------------------------------------------------------------
# Strand specificity tests
# ---------------------------------------------------------------------------


class TestStrandSpecificity:
    def test_same_strand_matches(self):
        svc = IntersectionService()
        assert hits(svc, pred(10, 20, strand="+"), trans(5, 30, strand="+")) == 1
        assert hits(svc, pred(10, 20, strand="-"), trans(5, 30, strand="-")) == 1

    def test_opposite_strand_no_match(self):
        svc = IntersectionService()
        assert hits(svc, pred(10, 20, strand="+"), trans(5, 30, strand="-")) == 0
        assert hits(svc, pred(10, 20, strand="-"), trans(5, 30, strand="+")) == 0
