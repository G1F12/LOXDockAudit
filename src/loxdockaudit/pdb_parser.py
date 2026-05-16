"""PDB parsing utilities for LOXDockAudit."""

from __future__ import annotations

import hashlib
from io import StringIO
from pathlib import Path

from Bio.PDB import PDBParser, Polypeptide
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure

from loxdockaudit.preprocessing import clean_pdb_lines


def load_structure(pdb_path: str, structure_id: str = "model") -> Structure:
    """Load a PDB file and return a Biopython Structure object."""
    path = Path(pdb_path)
    if not path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    parser = PDBParser(QUIET=True)
    raw_lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    cleaned_lines = clean_pdb_lines(raw_lines)
    if cleaned_lines:
        return parser.get_structure(structure_id, StringIO("".join(cleaned_lines)))
    return parser.get_structure(structure_id, str(path))


def get_atoms_by_selection(
    structure: Structure,
    chains: list[str],
    resi_range: tuple[int, int] | None = None,
    resi_list: list[int] | None = None,
    atom_names: list[str] | None = None,
) -> list[Atom]:
    """
    Return a list of Atom objects filtered by chain, residue, and atom name.

    Hydrogens are always excluded.
    """
    selected_atoms: list[Atom] = []
    chain_set = set(chains)
    resi_set = set(resi_list) if resi_list is not None else None
    atom_name_set = set(atom_names) if atom_names is not None else None

    for model in structure:
        for chain in model:
            if chain.id not in chain_set:
                continue

            for residue in chain:
                residue_number = residue.id[1]
                if resi_range is not None:
                    start, end = resi_range
                    if not start <= residue_number <= end:
                        continue
                if resi_set is not None and residue_number not in resi_set:
                    continue

                for atom in residue:
                    if _is_hydrogen(atom):
                        continue
                    if atom_name_set is not None and atom.name not in atom_name_set:
                        continue
                    selected_atoms.append(atom)

    return selected_atoms


def get_sequence_from_chain(structure: Structure, chain_id: str) -> str:
    """
    Return one-letter amino acid sequence for a chain.

    HETATM residues are skipped. Returns an empty string if the chain is not
    found.
    """
    model = next(structure.get_models(), None)
    if model is None or chain_id not in model:
        return ""

    residues: list[str] = []
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue

        one_letter = _three_to_one(residue.resname)
        if one_letter is not None:
            residues.append(one_letter)

    return "".join(residues)


def get_chain_ids(structure: Structure) -> list[str]:
    """Return list of all chain IDs present in first model."""
    model = next(structure.get_models(), None)
    if model is None:
        return []
    return [chain.id for chain in model]


def compute_sha256(pdb_path: str) -> str:
    """Return hex SHA256 hash of file contents."""
    path = Path(pdb_path)
    if not path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_hydrogen(atom: Atom) -> bool:
    return atom.element.strip().upper() == "H"


def _three_to_one(residue_name: str) -> str | None:
    converter = getattr(Polypeptide, "three_to_one", None)
    normalized_name = residue_name.strip().upper()

    try:
        if converter is not None:
            return converter(normalized_name)
        return Polypeptide.index_to_one(Polypeptide.three_to_index(normalized_name))
    except (KeyError, ValueError):
        return None
