"""A run shows only the tracks it asked for. No model is loaded here.

Head ids and tissues are neutral placeholders, not a reused example locus.
"""

import unittest

from ntv3_tools.browser import _ucsc_lane
from ntv3_tools.errors import RequestError
from ntv3_tools.registry import head_indices, model_head_for
from ntv3_tools.resolved import gene_and_clinvar_lanes

_HEADS = ["ENCSR100AAA", "ENCSR200BBB", "CNhs99999_P"]


class HeadIndexTests(unittest.TestCase):
    def test_exact_heads_map_to_their_columns(self) -> None:
        self.assertEqual(
            head_indices(names=_HEADS, heads=["CNhs99999_P", "ENCSR100AAA"]), [2, 0]
        )

    def test_a_missing_id_raises_and_is_named(self) -> None:
        with self.assertRaisesRegex(RequestError, r"\['ENCSR300CCC'\]") as caught:
            head_indices(names=_HEADS, heads=["ENCSR100AAA", "ENCSR300CCC"])
        self.assertIn("No other head is substituted", str(caught.exception))

    def test_a_near_match_is_never_substituted(self) -> None:
        # "ENCSR100" is a substring of exactly one head; it must not stand in for it.
        with self.assertRaisesRegex(RequestError, r"\['ENCSR100'\]"):
            head_indices(names=_HEADS, heads=["ENCSR100"])


class ModelHeadTests(unittest.TestCase):
    def test_a_row_with_a_model_head_keeps_it(self) -> None:
        row = {
            "query": "liver DNase",
            "gt_accession": "ENCSR100AAA",
            "model_head": "CNhs99999_P",
        }
        self.assertEqual(model_head_for(row=row, available=_HEADS), "CNhs99999_P")

    def test_an_encode_accession_is_not_used_as_the_model_head(self) -> None:
        row = {"query": "liver DNase", "gt_accession": "ENCSR200BBB", "model_head": ""}
        with self.assertRaisesRegex(
            RequestError,
            "ENCODE accession ENCSR200BBB would have been substituted as model head "
            "ENCSR200BBB",
        ):
            model_head_for(row=row, available=_HEADS)

    def test_no_head_and_no_matching_accession_also_raises(self) -> None:
        row = {"query": "liver DNase", "gt_accession": "ENCSR900ZZZ"}
        with self.assertRaisesRegex(RequestError, "no validated model head"):
            model_head_for(row=row, available=_HEADS)


class LaneQuestionTests(unittest.TestCase):
    def test_the_human_lanes_offered_are_genes_and_clinvar(self) -> None:
        self.assertEqual(
            [choice.track for choice in gene_and_clinvar_lanes(assembly="hg38")],
            ["knownGene", "clinvarMain"],
        )

    def test_the_mouse_gene_lane_offered_is_ncbi_refseq(self) -> None:
        self.assertEqual(
            [choice.track for choice in gene_and_clinvar_lanes(assembly="mm10")],
            ["ncbiRefSeq"],
        )

    def test_a_mouse_gene_row_shows_name2_only_when_it_has_one(self) -> None:
        # Saved ncbiRefSeq payload, genePred columns; no network.
        row = {
            "name": "NM_000000002.1",
            "name2": "Gene4",
            "chrom": "chr5",
            "strand": "-",
            "txStart": 500,
            "txEnd": 900,
            "exonStarts": "500,800,",
            "exonEnds": "600,900,",
        }
        payload = {"ncbiRefSeq": [row, {**row, "name": "NM_000000003.1", "name2": ""}]}

        lane = _ucsc_lane(track="ncbiRefSeq", label="NCBI RefSeq", payload=payload)

        self.assertEqual(lane["kind"], "genes")
        self.assertEqual(
            [item["n"] for item in lane["items"]], ["Gene4", "NM_000000003.1"]
        )
        self.assertEqual(lane["items"][0]["b"], [[500, 600], [800, 900]])


if __name__ == "__main__":
    unittest.main()
