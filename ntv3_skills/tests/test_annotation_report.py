"""A finished GFF3 in a run becomes the NTv3 annotation lane of report.json."""

import gzip
import json
import tempfile
import unittest
from pathlib import Path

from ntv3_tools.annotation_report import annotation_payload, read_gff3
from ntv3_tools.errors import ArtifactError, RequestError
from ntv3_tools.report import write_viewer_bundle

from tests.run_fixtures import (
    ANNOTATION_GFF3,
    ANNOTATION_PLAN,
    predicted_run,
    write_annotation,
)


class ReadGff3Tests(unittest.TestCase):
    def test_genes_on_both_strands_become_half_open_lanes_with_cds_exons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "annotation.gff3")
            path.write_text(ANNOTATION_GFF3)
            gff3 = read_gff3(path)

        self.assertEqual(
            gff3.genes,
            [
                {
                    "chrom": "chr11",
                    "name": "gene_1_p",
                    "start": 5_266_000,
                    "end": 5_267_000,
                    "strand": "+",
                    "exons": [
                        {"start": 5_266_000, "end": 5_266_200},
                        {"start": 5_266_700, "end": 5_267_000},
                    ],
                },
                {
                    "chrom": "chr11",
                    "name": "gene_2_m",
                    "start": 5_268_000,
                    "end": 5_269_500,
                    "strand": "-",
                    "exons": [
                        {"start": 5_268_000, "end": 5_268_300},
                        {"start": 5_269_200, "end": 5_269_500},
                    ],
                },
            ],
        )
        self.assertEqual(gff3.regions, [{"chrom": "chr11", "length": 135_086_622}])
        self.assertEqual(gff3.provenance["decoder"], "hmm_constrained")
        self.assertEqual(gff3.provenance["min-cds-length"], "201")

    def test_gzipped_gff3_reads_the_same(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "annotation.gff3.gz")
            path.write_bytes(gzip.compress(ANNOTATION_GFF3.encode()))
            gff3 = read_gff3(path)

        self.assertEqual(
            [gene["name"] for gene in gff3.genes], ["gene_1_p", "gene_2_m"]
        )

    def test_cds_without_its_gene_is_refused(self) -> None:
        orphan = ANNOTATION_GFF3.replace("Parent=gene_2_m", "Parent=gene_9_m")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "annotation.gff3")
            path.write_text(orphan)
            with self.assertRaisesRegex(ArtifactError, "gene_9_m"):
                read_gff3(path)


class AnnotationPayloadTests(unittest.TestCase):
    def test_explain_values_come_from_the_gff3_header_and_the_plan(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            payload = annotation_payload(root=Path(tmp))

        assert payload is not None
        self.assertEqual(payload["label"], "NTv3 annotation")
        self.assertEqual(payload["model"], "InstaDeepAI/NTv3_annot_mammals")
        self.assertEqual(payload["modelName"], "new_ntv3_700M_annotate_mammals_1mb")
        self.assertEqual(payload["clade"], "mammals")
        self.assertEqual(payload["decoder"], "hmm_constrained")
        self.assertEqual(payload["windowBp"], 1_048_576)
        self.assertEqual(payload["keepFraction"], 0.75)
        self.assertEqual(payload["minCdsLength"], 201)
        self.assertTrue(payload["inFrameStopFilter"])
        self.assertEqual(payload["gff3"], "annotation.gff3")
        self.assertEqual(
            payload["sources"], ["annotation_plan.json", "annotation.gff3"]
        )
        self.assertEqual(len(payload["genes"]), 2)
        self.assertNotIn("probabilities", payload)

    def test_the_gff3_header_wins_over_the_plan(self) -> None:
        header = ANNOTATION_GFF3.replace(
            "##ntv3-decoder hmm_constrained", "##ntv3-decoder argmax"
        )
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp), gff3=header)
            payload = annotation_payload(root=Path(tmp))

        assert payload is not None
        self.assertEqual(payload["decoder"], "argmax")

    def test_unrecorded_values_stay_unrecorded(self) -> None:
        bare = "\n".join(
            line for line in ANNOTATION_GFF3.splitlines() if "##ntv3-" not in line
        )
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "genes.gff3").write_text(bare + "\n")
            payload = annotation_payload(root=Path(tmp))

        assert payload is not None
        for key in (
            "model",
            "clade",
            "decoder",
            "windowBp",
            "keepFraction",
            "minCdsLength",
            "inFrameStopFilter",
        ):
            self.assertIsNone(payload[key], key)
        self.assertEqual(payload["sources"], ["genes.gff3"])

    def test_no_gff3_means_no_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(annotation_payload(root=Path(tmp)))

    def test_two_gff3_files_are_refused_rather_than_guessed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "a.gff3").write_text(ANNOTATION_GFF3)
            Path(tmp, "b.gff3").write_text(ANNOTATION_GFF3)
            with self.assertRaisesRegex(ArtifactError, "a.gff3"):
                annotation_payload(root=Path(tmp))


class ReportAnnotationTests(unittest.TestCase):
    def test_predict_run_with_a_gff3_carries_the_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            write_annotation(Path(tmp))
            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())

        self.assertEqual(report["tracks"]["assays"][0]["modelHead"], "CNhs11250_P")
        self.assertEqual(
            [gene["name"] for gene in report["annotation"]["genes"]],
            ["gene_1_p", "gene_2_m"],
        )

    def test_predict_run_without_a_gff3_has_no_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())

        self.assertIsNone(report["annotation"])

    def test_annotate_run_is_the_same_report_without_tracks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())
            written = sorted(path.name for path in Path(tmp).iterdir())

        self.assertEqual(report["tracks"]["assays"], [])
        self.assertIsNone(report["explain"])
        self.assertEqual(report["locus"]["chrom"], "chr11")
        self.assertEqual(report["locus"]["start"], 0)
        self.assertEqual(report["locus"]["end"], 135_086_622)
        self.assertEqual(report["modelShort"], ANNOTATION_PLAN["model"])
        self.assertEqual(report["annotation"]["clade"], "mammals")
        self.assertEqual(
            written,
            [
                "annotation.gff3",
                "annotation_plan.json",
                "browser.html",
                "critic",
                "ntv3_report.html",
                "report.json",
                "ro-crate-metadata.json",
            ],
        )

    def test_annotate_run_without_a_gff3_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "annotation_plan.json").write_text(json.dumps(ANNOTATION_PLAN))
            with self.assertRaisesRegex(ArtifactError, "no GFF3"):
                write_viewer_bundle(run_dir=tmp)
            self.assertFalse(Path(tmp, "report.json").exists())


_TWO_CONTIGS = (
    "##gff-version 3\n"
    "##sequence-region chrA 1 5000\n"
    "##sequence-region chrB 1 9000\n"
    "chrA\tNTv3\tgene\t101\t900\t.\t+\t.\tID=a1\n"
    "chrA\tNTv3\tCDS\t101\t900\t.\t+\t0\tParent=a1\n"
    "chrB\tNTv3\tgene\t101\t900\t.\t+\t.\tID=b1\n"
    "chrB\tNTv3\tCDS\t101\t900\t.\t+\t0\tParent=b1\n"
    "chrB\tNTv3\tgene\t2001\t2900\t.\t-\t.\tID=b2\n"
    "chrB\tNTv3\tCDS\t2001\t2900\t.\t-\t0\tParent=b2\n"
)


class AnnotateContigTests(unittest.TestCase):
    def test_the_page_shows_the_contig_with_the_most_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp), gff3=_TWO_CONTIGS)
            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())

        self.assertEqual(
            (report["locus"]["chrom"], report["locus"]["end"]), ("chrB", 9000)
        )
        self.assertEqual(
            {
                k: report["runParameters"][k]
                for k in ("contig", "genes_on_contig", "contigs", "genes")
            },
            {"contig": "chrB", "genes_on_contig": 2, "contigs": 2, "genes": 3},
        )

    def test_contig_chooses_another_and_an_unknown_one_is_refused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp), gff3=_TWO_CONTIGS)
            report = json.loads(
                write_viewer_bundle(run_dir=tmp, contig="chrA").read_text()
            )
            with self.assertRaisesRegex(RequestError, "Which contig\\? chrZ"):
                write_viewer_bundle(run_dir=tmp, refresh=True, contig="chrZ")

        self.assertEqual(
            (report["locus"]["chrom"], report["locus"]["end"]), ("chrA", 5000)
        )

    def test_contig_is_refused_for_a_predict_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            predicted_run(Path(tmp))
            with self.assertRaisesRegex(RequestError, "--contig"):
                write_viewer_bundle(run_dir=tmp, contig="chr11")


class AnnotateCriticTests(unittest.TestCase):
    def test_the_gene_claim_is_audited_with_its_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            Path(tmp, "annotation_manifest.json").write_text('{"genes": 2}')
            report = json.loads(write_viewer_bundle(run_dir=tmp).read_text())
            log = json.loads(Path(tmp, "critic", "audit_log.json").read_text())

        self.assertEqual(
            (report["critic"]["kept"], report["critic"]["suppressed"]), (1, 0)
        )
        self.assertEqual([row["claim_id"] for row in log], ["annotate"])
        self.assertIn(
            "predicted 2 genes with the hmm_constrained decoder on chr11",
            log[0]["text"],
        )

    def test_no_manifest_suppresses_the_claim_and_off_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            suppressed = json.loads(write_viewer_bundle(run_dir=tmp).read_text())
        with tempfile.TemporaryDirectory() as tmp:
            write_annotation(Path(tmp))
            plan = json.loads(Path(tmp, "annotation_plan.json").read_text())
            plan["want_critic"] = False
            Path(tmp, "annotation_plan.json").write_text(json.dumps(plan))
            off = json.loads(write_viewer_bundle(run_dir=tmp).read_text())

        self.assertEqual(suppressed["critic"]["suppressed"], 1)
        self.assertEqual(off["critic"]["enabled"], False)

    def test_a_locus_annotation_added_on_refresh_is_audited(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            predicted_run(root)
            plan = root / "plan.yaml"
            plan.write_text(
                plan.read_text().replace("want_critic: false", "want_critic: true")
            )
            write_viewer_bundle(run_dir=tmp)
            before = json.loads((root / "critic" / "audit_log.json").read_text())
            write_annotation(root)
            (root / "annotation_manifest.json").write_text('{"genes": 2}')
            write_viewer_bundle(run_dir=tmp, refresh=True)
            after = json.loads((root / "critic" / "audit_log.json").read_text())

        self.assertNotIn("annotate", [row["claim_id"] for row in before])
        self.assertIn("annotate", [row["claim_id"] for row in after])


if __name__ == "__main__":
    unittest.main()
