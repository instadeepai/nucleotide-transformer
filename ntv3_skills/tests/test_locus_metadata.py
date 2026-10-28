"""Sourced UCSC knownGene locus-metadata tests."""

import json
import tempfile
import unittest
from pathlib import Path

from ntv3_tools.locus import known_gene_metadata
from ntv3_tools.lookups import LookupRequest, run_lookup


class KnownGeneMetadataTests(unittest.TestCase):
    def test_hbe1_minus_strand_uses_zero_based_last_base_as_tss(self) -> None:
        payload = {
            "_request_url": (
                "https://api.genome.ucsc.edu/getData/track?"
                "genome=hg38;track=knownGene;chrom=chr11;start=5260000;end=5276384"
            ),
            "_label": "knownGene",
            "knownGene": [
                {
                    "name": "ENST00000380237.6",
                    "geneName": "HBE1",
                    "chrom": "chr11",
                    "strand": "-",
                    "chromStart": 5_281_649,
                    "chromEnd": 5_283_028,
                    "rank": 2,
                    "tag": "GENCODE Primary",
                },
                {
                    "name": "ENST00000499845.2",
                    "geneName": "HBE1",
                    "chrom": "chr11",
                    "strand": "-",
                    "chromStart": 5_281_600,
                    "chromEnd": 5_283_200,
                    "rank": 1,
                },
                {
                    "name": "ENST00000380236.2",
                    "geneName": "HBG2",
                    "chrom": "chr11",
                    "strand": "-",
                    "chromStart": 5_276_000,
                    "chromEnd": 5_278_000,
                    "rank": 1,
                    "tag": "MANE Select",
                },
            ],
        }

        metadata = known_gene_metadata(
            gene="HBE1",
            chrom="chr11",
            start=5_260_000,
            end=5_290_000,
            payload=payload,
        )

        self.assertEqual(metadata["strand"], "-")
        self.assertEqual(metadata["tss"], 5_283_027)
        self.assertEqual(metadata["transcript"], "ENST00000380237.6")
        self.assertEqual(metadata["source"]["track"], "knownGene")
        self.assertIn("HBE1 · gene", metadata["elements"])
        self.assertIn("HBG2 · nearby gene", metadata["elements"])

    def test_explicit_coordinates_keep_null_focus_and_named_nearby_gene(self) -> None:
        payload = {
            "_request_url": "https://api.genome.ucsc.edu/getData/track?track=knownGene",
            "_label": "knownGene",
            "knownGene": [
                {
                    "name": "ENST_NEARBY",
                    "geneName": "HBG1",
                    "chrom": "chr11",
                    "strand": "-",
                    "chromStart": 5_270_000,
                    "chromEnd": 5_271_000,
                    "rank": 1,
                }
            ],
        }

        metadata = known_gene_metadata(
            gene="NOT_A_GENE",
            chrom="chr11",
            start=5_260_000,
            end=5_276_384,
            payload=payload,
        )

        self.assertIsNone(metadata["tss"])
        self.assertIsNone(metadata["strand"])
        self.assertEqual(metadata["elements"], ["HBG1 · nearby gene"])

    def test_lookup_fetches_only_the_approved_lanes(self) -> None:
        for lanes, expected, has_metadata in (
            (["clinvarMain"], ["clinvarMain"], False),
            ([], [], False),
            (["knownGene"], ["knownGene"], True),
        ):
            with self.subTest(lanes=lanes):
                seen, metadata = self._lookup(ucsc_tracks=lanes)
                self.assertEqual(seen, expected)
                self.assertEqual(metadata is not None, has_metadata)
                if metadata is not None:
                    self.assertEqual(metadata["source"]["track"], "knownGene")

    def _lookup(self, ucsc_tracks: list[str]) -> tuple[list[str], dict | None]:
        seen_tracks: list[str] = []

        def encode_get(
            _url: str,
            _params: dict[str, str] | None = None,
        ) -> dict[str, object]:
            return {
                "@graph": [
                    {
                        "accession": "ENCFF079PRE",
                        "assembly": "GRCh38",
                        "output_type": "signal",
                        "href": (
                            "https://www.encodeproject.org/files/ENCFF079PRE/"
                            "@@download/ENCFF079PRE.bigWig"
                        ),
                    }
                ]
            }

        def fetch_ucsc(**kwargs: object) -> dict[str, object]:
            track = str(kwargs["track"])
            seen_tracks.append(track)
            return {
                "_label": track,
                "_request_url": f"https://api.genome.ucsc.edu/{track}",
                track: [],
            }

        with tempfile.TemporaryDirectory() as tmp:
            run_lookup(
                request=LookupRequest(
                    species="human",
                    tissue="K562",
                    names=["K562 CAGE"],
                    out_dir=tmp,
                    chrom="chr11",
                    start=5_260_000,
                    end=5_276_384,
                    ucsc_tracks=ucsc_tracks,
                    gene="HBE1",
                ),
                get_json=encode_get,
                fetch_ucsc=fetch_ucsc,
            )
            metadata_path = Path(tmp, "locus_metadata.json")
            metadata = (
                json.loads(metadata_path.read_text()) if metadata_path.exists() else None
            )
        return seen_tracks, metadata


if __name__ == "__main__":
    unittest.main()
