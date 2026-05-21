"""Minimal command-line smoke test for loading a PDB input file."""

from __future__ import annotations

import argparse
from pathlib import Path

from loxdockaudit.pdb_parser import get_chain_ids, load_structure


def main() -> None:
    """Parse a PDB file and print a compact structure summary."""
    parser = argparse.ArgumentParser(description="Load a PDB file with LOXDockAudit.")
    parser.add_argument("--input", required=True, type=Path, help="Path to an input PDB file")
    args = parser.parse_args()

    structure = load_structure(str(args.input))
    chains = get_chain_ids(structure)
    print(f"Loaded {args.input}")
    print(f"Chains: {', '.join(chains) if chains else 'none'}")


if __name__ == "__main__":
    main()
