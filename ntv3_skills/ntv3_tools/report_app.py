"""Write one run as the NTv3 report app: one HTML file a Claude artifact renders.

``report_app.html`` is the report app (the report sections and GenomeBrowser)
built into a single file by ``report_app_build.mjs``, with one empty
``ntv3-report`` data tag. ``write_report_app`` puts the run's ``report.json``
in that tag and writes ``results/<run>/ntv3_report.html``. The page loads no
script or stylesheet from the network, and draws the UCSC rows the run saved
without calling UCSC. Report time needs neither Node nor the app source.

The app source is ``ntv3_skills/report/``. Rebuild the shell after it changes
(after ``pnpm install --frozen-lockfile`` there):

    python -m ntv3_tools.report_app --app ntv3_skills/report
"""

import argparse
import json
import os
import subprocess
from pathlib import Path
from typing import Any

from ntv3_tools.browser import json_for_script
from ntv3_tools.errors import ArtifactError

SHELL_PATH = Path(__file__).with_name("report_app.html")
BUILD_SCRIPT = Path(__file__).with_name("report_app_build.mjs")
REPORT_APP_HTML = "ntv3_report.html"
DATA_TAG = '<script type="application/json" id="ntv3-report">null</script>'


def render_report_app(report: dict[str, Any]) -> str:
    """Return the report app page with one run's report embedded.

    Args:
        report: Parsed ``report.json``.

    Returns:
        The HTML document.

    Raises:
        ArtifactError: The shell is missing or does not hold exactly one empty
            data tag.
    """
    if not SHELL_PATH.is_file():
        raise ArtifactError(
            f"{SHELL_PATH.name} is missing. Build it with "
            "`python -m ntv3_tools.report_app --app ntv3_skills/report`."
        )
    shell = SHELL_PATH.read_text()
    if shell.count(DATA_TAG) != 1:
        raise ArtifactError(f"{SHELL_PATH.name} must hold exactly one empty data tag.")
    return shell.replace(
        DATA_TAG, DATA_TAG.replace(">null<", f">{json_for_script(report)}<")
    )


def write_report_app(run_dir: str | Path) -> Path:
    """Write ``ntv3_report.html`` from the run's ``report.json``.

    Args:
        run_dir: Run directory holding ``report.json``.

    Returns:
        Path to ``ntv3_report.html``, the file to publish as the artifact.

    Raises:
        ArtifactError: ``report.json`` or the shell is missing.
    """
    root = Path(run_dir)
    report_path = root / "report.json"
    if not report_path.is_file():
        raise ArtifactError(f"{REPORT_APP_HTML} needs {report_path}. Run report first.")
    path = root / REPORT_APP_HTML
    path.write_text(render_report_app(report=json.loads(report_path.read_text())))
    return path


def build_shell(app_dir: str | Path, node: str = "node") -> Path:
    """Rebuild ``report_app.html`` from the report app with its own Vite.

    Args:
        app_dir: Report app directory with ``index.html``, ``src/`` and
            ``node_modules/``.
        node: Node executable.

    Returns:
        Path to the rebuilt shell.
    """
    subprocess.run(
        [node, str(BUILD_SCRIPT), str(Path(app_dir).resolve()), str(SHELL_PATH)],
        check=True,
    )
    return SHELL_PATH


def main() -> None:
    """Rebuild the shell: ``python -m ntv3_tools.report_app --app <dir>``."""
    parser = argparse.ArgumentParser(prog="ntv3-report-app")
    parser.add_argument(
        "--app",
        default=os.environ.get("NTV3_REPORT_APP"),
        help="Report app directory (default $NTV3_REPORT_APP)",
    )
    parser.add_argument("--node", default="node")
    args = parser.parse_args()
    if not args.app:
        parser.error("Pass --app or set NTV3_REPORT_APP to the report app directory.")
    build_shell(app_dir=args.app, node=args.node)


if __name__ == "__main__":
    main()
