# PDB Preprocessing

LOXDockAudit v0.2.1 includes an official preprocessing path for HDOCK-style PDB
files. The parser uses this path before handing structures to Biopython.

## Rules

- Strip `MODEL` records.
- Strip `ENDMDL` records.
- Strip `REMARK` records.
- Preserve `TER` records.
- Ignore duplicate or upstream `END` records and emit one final `END` when
  coordinates are retained.
- Drop malformed `ATOM` or `HETATM` rows whose coordinate columns cannot be
  parsed.
- Normalize blank chain IDs to `X`.
- Detect HDOCK `HEADER lig*.pdb` sections and map ligand coordinate records to
  chain `D`, matching the packaged LOXDockAudit config convention.

## Python API

```python
from loxdockaudit.preprocessing import clean_pdb_file, clean_pdb_directory

clean_pdb_file("raw/model_1.pdb", "cleaned/model_1.pdb")
clean_pdb_directory("raw_models", "cleaned_models")
```

`pdb_parser.load_structure()` automatically applies the same preprocessing
rules in memory, so normal CLI runs benefit from the hardening without requiring
a separate cleaning command.
