from Bio.PDB import PDBParser

path = "examples/data/r5_real_clean2/model_1.pdb"

parser = PDBParser(QUIET=False, PERMISSIVE=True)

try:
    structure = parser.get_structure("x", path)
    print("PARSED OK")
    for chain in structure.get_chains():
        residues = list(chain.get_residues())
        atoms = list(chain.get_atoms())
        print("chain", chain.id, "residues", len(residues), "atoms", len(atoms))
except Exception as e:
    print("PARSE FAILED")
    print(type(e).__name__)
    print(e)