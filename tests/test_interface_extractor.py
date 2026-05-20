from __future__ import annotations

from pathlib import Path

from loxdockaudit.interface_extractor import (
    contact_pairs_between_chains,
    get_active_site_contacts,
    get_interface_residues,
)


def test_get_interface_residues_detects_chain_contacts(tmp_path: Path) -> None:
    pdb_path = tmp_path / "interface.pdb"
    pdb_path.write_text(_interface_pdb_text(), encoding="utf-8")

    interface = get_interface_residues(pdb_path, "A", "D", cutoff_angstrom=5.0)

    assert interface.chain_a_residues == frozenset({11})
    assert interface.chain_b_residues == frozenset({100})


def test_contact_pairs_between_chains_returns_residue_pairs(tmp_path: Path) -> None:
    pdb_path = tmp_path / "interface.pdb"
    pdb_path.write_text(_interface_pdb_text(), encoding="utf-8")

    pairs = contact_pairs_between_chains(pdb_path, "A", "D", cutoff_angstrom=5.0)

    assert pairs == frozenset({(11, 100)})


def test_get_active_site_contacts_returns_closest_substrate_residue(tmp_path: Path) -> None:
    pdb_path = tmp_path / "interface.pdb"
    pdb_path.write_text(_interface_pdb_text(), encoding="utf-8")

    contacts = get_active_site_contacts(pdb_path, [100], substrate_chain="A")

    assert len(contacts) == 1
    assert contacts[0].receptor_resnum == 100
    assert contacts[0].ligand_resnum == 11
    assert contacts[0].distance == 3.0


def _interface_pdb_text() -> str:
    return "\n".join(
        [
            "ATOM      1  CB  LYS A  11       3.000   0.000   0.000  1.00 80.00           C  ",
            "ATOM      2  NZ  LYS A  11       4.000   0.000   0.000  1.00 80.00           N  ",
            "ATOM      3  CA  ALA A  12      30.000   0.000   0.000  1.00 80.00           C  ",
            "ATOM      4  CA  HIS D 100       0.000   0.000   0.000  1.00 80.00           C  ",
            "ATOM      5  CA  HIS D 101      20.000   0.000   0.000  1.00 80.00           C  ",
            "TER",
            "END",
        ]
    )
