"""Unit tests for deterministic literature record packaging."""

import json
import tempfile
import unittest
from pathlib import Path

from ntv3_tools.critic import Claim, audit_claims
from ntv3_tools.errors import RequestError
from ntv3_tools.literature import package_literature
from ntv3_tools.report import _trace_from_tree


def valid_record(**overrides: object) -> dict[str, object]:
    """Create one valid primary-source literature record."""
    record: dict[str, object] = {
        "citekey": "Avsec2021Enformer",
        "title": "Effective gene expression prediction from sequence",
        "authors": ["Žiga Avsec", "Vikram Agarwal"],
        "year": 2021,
        "source_url": "https://www.nature.com/articles/s41592-021-01252-x",
        "doi": "10.1038/s41592-021-01252-x",
        "retrieval_date": "2026-09-08",
        "search_query": "sequence model gene expression primary paper",
        "relevance": "Supports the claim that long-range sequence context helps.",
    }
    record.update(overrides)
    return record


def write_records(run_dir: Path, records: list[dict[str, object]]) -> Path:
    """Write test input records under a run directory."""
    path = run_dir / "literature_records.json"
    path.write_text(json.dumps(records))
    return path


class LiteratureValidationTests(unittest.TestCase):
    def test_rejects_missing_blank_and_duplicate_fields(self) -> None:
        invalid_records = [
            [valid_record(title=" ")],
            [valid_record(), valid_record()],
            [valid_record(retrieval_date="09/08/2026")],
            [valid_record(source_url="ftp://nature.com/paper")],
        ]
        for records in invalid_records:
            with (
                self.subTest(records=records),
                tempfile.TemporaryDirectory() as tmp,
            ):
                run_dir = Path(tmp)
                (run_dir / "plan.yaml").touch()
                records_path = write_records(run_dir=run_dir, records=records)
                with self.assertRaises(RequestError):
                    package_literature(
                        run_dir=run_dir,
                        records_path=records_path,
                    )

    def test_rejects_placeholder_and_non_primary_sources(self) -> None:
        invalid_records = [
            [valid_record(citekey="TODO2026")],
            [
                valid_record(
                    source_url="https://example.com/paper",
                    doi=None,
                )
            ],
            [
                valid_record(
                    source_url="https://scholar.google.com/scholar?q=enformer",
                    doi=None,
                )
            ],
            [valid_record(doi="doi:placeholder")],
        ]
        for records in invalid_records:
            with (
                self.subTest(records=records),
                tempfile.TemporaryDirectory() as tmp,
            ):
                run_dir = Path(tmp)
                (run_dir / "plan.yaml").touch()
                records_path = write_records(run_dir=run_dir, records=records)
                with self.assertRaises(RequestError):
                    package_literature(
                        run_dir=run_dir,
                        records_path=records_path,
                    )

    def test_rejects_records_file_outside_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run_dir = root / "run"
            run_dir.mkdir()
            (run_dir / "plan.yaml").touch()
            records_path = write_records(run_dir=root, records=[valid_record()])

            with self.assertRaises(RequestError):
                package_literature(run_dir=run_dir, records_path=records_path)


class LiteratureOutputTests(unittest.TestCase):
    def test_report_trace_includes_packaged_literature(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "plan.yaml").touch()
            records_path = write_records(run_dir=run_dir, records=[valid_record()])
            package_literature(run_dir=run_dir, records_path=records_path)

            literature_trace = next(
                row
                for row in _trace_from_tree(root=run_dir)
                if row["name"] == "literature"
            )

            self.assertEqual(literature_trace["status"], "ok")
            self.assertEqual(
                literature_trace["artifact"],
                "literature/search_manifest.json",
            )

    def test_output_is_deterministic_and_escaped(self) -> None:
        with (
            tempfile.TemporaryDirectory() as first_tmp,
            tempfile.TemporaryDirectory() as second_tmp,
        ):
            outputs: list[tuple[str, str]] = []
            for tmp in (first_tmp, second_tmp):
                run_dir = Path(tmp)
                (run_dir / "plan.yaml").touch()
                records_path = write_records(
                    run_dir=run_dir,
                    records=[
                        valid_record(
                            title="Sequence & expression: 100% reproducible",
                        )
                    ],
                )
                bib_path, manifest_path = package_literature(
                    run_dir=run_dir,
                    records_path=records_path,
                )
                outputs.append((bib_path.read_text(), manifest_path.read_text()))

            self.assertEqual(outputs[0], outputs[1])
            self.assertIn(r"Sequence \& expression: 100\% reproducible", outputs[0][0])
            manifest = json.loads(outputs[0][1])
            self.assertEqual(
                manifest["records"][0]["search_query"],
                "sequence model gene expression primary paper",
            )
            self.assertEqual(
                manifest["records"][0]["relevance"],
                "Supports the claim that long-range sequence context helps.",
            )

    def test_bibtex_is_compatible_with_critic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "plan.yaml").touch()
            records_path = write_records(run_dir=run_dir, records=[valid_record()])
            bib_path, _ = package_literature(
                run_dir=run_dir,
                records_path=records_path,
            )

            decisions = audit_claims(
                claims=[
                    Claim(
                        claim_id="literature",
                        text="Long-range context supports prediction.",
                        evidence_paths=["literature/citations.bib"],
                        citekey="Avsec2021Enformer",
                        kind="literature",
                    )
                ],
                run_dir=str(run_dir),
                bib_path=str(bib_path),
            )

            self.assertEqual(decisions[0].verdict, "KEEP")

    def test_existing_output_refuses_without_modification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "plan.yaml").touch()
            records_path = write_records(run_dir=run_dir, records=[valid_record()])
            literature_dir = run_dir / "literature"
            literature_dir.mkdir()
            bib_path = literature_dir / "citations.bib"
            bib_path.write_text("keep me")

            with self.assertRaises(RequestError):
                package_literature(run_dir=run_dir, records_path=records_path)

            self.assertEqual(bib_path.read_text(), "keep me")
            self.assertFalse((literature_dir / "search_manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
