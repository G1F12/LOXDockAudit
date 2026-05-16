from pathlib import Path

SRC = Path("examples/data/r5_real")
DST = Path("examples/data/r5_real_clean2")
DST.mkdir(parents=True, exist_ok=True)

def normalize_atom_name(atom_name: str) -> str:
    return atom_name[:4].ljust(4)

for i in range(1, 11):
    inp = SRC / f"model_{i}.pdb"
    out = DST / f"model_{i}.pdb"

    section = None
    atom_serial = 1
    residue_seen = {}
    lines_out = []

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

            record = "ATOM"
            atom_name = normalize_atom_name(line[12:16].strip())
            altloc = " "
            resname = line[17:20].strip()[:3].rjust(3)

            if section == "lig":
                chain = "D"
            else:
                chain = line[21] if len(line) > 21 and line[21].strip() else "A"

            try:
                resseq = int(line[22:26])
            except ValueError:
                continue

            icode = " "

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            occ = 1.00
            bfac = 0.00
            element = line[76:78].strip() if len(line) >= 78 else atom_name.strip()[0]
            element = element[:2].rjust(2)

            # Skip duplicate atom in same residue; Biopython can choke on these
            key = (chain, resseq, icode, resname, atom_name)
            if key in residue_seen:
                continue
            residue_seen[key] = True

            new_line = (
                f"{record:<6}{atom_serial:5d} "
                f"{atom_name:<4}{altloc}"
                f"{resname:>3} "
                f"{chain}{resseq:4d}{icode}   "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{occ:6.2f}{bfac:6.2f}          "
                f"{element:>2}"
                "\n"
            )
            lines_out.append(new_line)
            atom_serial += 1

    lines_out.append("END\n")

    with out.open("w", encoding="ascii", newline="\n") as f:
        f.writelines(lines_out)

    print(f"saved {out} atoms={atom_serial - 1}")