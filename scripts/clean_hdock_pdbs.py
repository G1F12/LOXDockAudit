from pathlib import Path

SRC = Path("examples/data/r5_real")
DST = Path("examples/data/r5_real_clean")
DST.mkdir(parents=True, exist_ok=True)

for i in range(1, 11):
    inp = SRC / f"model_{i}.pdb"
    out = DST / f"model_{i}.pdb"

    lines_out = []
    section = None
    atom_serial = 1

    with inp.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")

            if line.startswith("HEADER rec"):
                section = "rec"
                continue

            if line.startswith("HEADER lig"):
                section = "lig"
                continue

            if not line.startswith(("ATOM", "HETATM")):
                continue

            # receptor оставляем как A/B/C если есть
            # ligand принудительно делаем chain D
            if section == "lig":
                chain = "D"
            else:
                chain = line[21] if len(line) > 21 and line[21].strip() else "A"

            # пересобираем строку PDB с новым serial и chain
            new_line = (
                f"{line[:6]}"
                f"{atom_serial:5d}"
                f"{line[11:21]}"
                f"{chain}"
                f"{line[22:]}"
            )

            lines_out.append(new_line + "\n")
            atom_serial += 1

    lines_out.append("END\n")

    with out.open("w", encoding="utf-8", newline="\n") as f:
        f.writelines(lines_out)

    print(f"saved {out} atoms={atom_serial - 1}")