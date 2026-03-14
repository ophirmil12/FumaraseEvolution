"""
chimerax_run.py - Generate ChimeraX .cxc scripts with correct paths and run them.

Reads the template .cxc files, replaces STRUCTURE_DIR with the actual path
from config.py, writes resolved copies to results/structure/, then runs each
via ChimeraX headless mode.

Usage:
    python scripts/structure/chimerax_run.py

    # Single script:
    python scripts/structure/chimerax_run.py --script class1_AF

    # Dry run (write resolved .cxc files without running ChimeraX):
    python scripts/structure/chimerax_run.py --dry-run

Requires ChimeraX to be installed and accessible as `chimerax` on PATH,
or set CHIMERAX_BIN environment variable:
    export CHIMERAX_BIN="/path/to/ChimeraX"
"""

import sys
import os
import subprocess
import argparse
import logging
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from config import PATHS

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

STRUCTURE_DIR = Path(PATHS["results"]) / "structure"
SCRIPTS_DIR   = Path(__file__).parent

# Template .cxc files (relative to this script)
SCRIPTS = {
    "class1_AF":       "chimerax_class1_AF.cxc",
    "class2_AF":       "chimerax_class2_AF.cxc",
    "class2_tetramer": "chimerax_class2_tetramer.cxc",
}


def resolve_cxc(template_path: Path, structure_dir: Path) -> Path:
    """
    Read template .cxc, replace STRUCTURE_DIR placeholder with actual path,
    write resolved file to structure_dir. Returns path to resolved file.
    """
    content = template_path.read_text()
    resolved = content.replace("STRUCTURE_DIR", str(structure_dir).replace("\\", "/"))

    out_path = structure_dir / template_path.name
    out_path.write_text(resolved)
    log.info(f"  Resolved .cxc -> {out_path}")
    return out_path


def run_chimerax(cxc_path: Path) -> None:
    """Run a .cxc script headlessly via ChimeraX."""
    chimerax_bin = os.environ.get("CHIMERAX_BIN", "chimerax")

    cmd = [chimerax_bin, "--nogui", "--script", str(cxc_path)]
    log.info(f"  Running: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout:
        log.info(result.stdout.strip())
    if result.stderr:
        log.warning(result.stderr.strip())
    if result.returncode != 0:
        raise RuntimeError(
            f"ChimeraX exited with code {result.returncode} for {cxc_path.name}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Run ChimeraX structure visualization scripts."
    )
    parser.add_argument(
        "--script",
        choices=list(SCRIPTS.keys()) + ["all"],
        default="all",
        help="Which script to run (default: all)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write resolved .cxc files without running ChimeraX",
    )
    args = parser.parse_args()

    STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)

    to_run = SCRIPTS if args.script == "all" else {args.script: SCRIPTS[args.script]}

    for label, cxc_filename in to_run.items():
        log.info(f"\n{'='*50}")
        log.info(f" {label}: {cxc_filename}")
        log.info(f"{'='*50}")

        template_path = SCRIPTS_DIR / cxc_filename
        if not template_path.exists():
            log.error(f"  Template not found: {template_path}")
            continue

        resolved_path = resolve_cxc(template_path, STRUCTURE_DIR)

        if args.dry_run:
            log.info("  [dry-run] Skipping ChimeraX execution")
            continue

        try:
            run_chimerax(resolved_path)
            log.info(f"  Done — outputs in {STRUCTURE_DIR}")
        except RuntimeError as e:
            log.error(str(e))

    log.info("\nAll done.")
    log.info(f"Outputs: {STRUCTURE_DIR}")


if __name__ == "__main__":
    main()
