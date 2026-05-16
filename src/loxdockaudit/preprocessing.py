"""PDB preprocessing utilities for HDOCK-style model files."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, TypedDict


SKIPPED_RECORDS = {"MODEL", "ENDMDL", "REMARK"}
COORDINATE_RECORDS = {"ATOM", "HETATM"}


class _ParsedCoordinate(TypedDict):
    record: str
    serial: int
    atom_name: str
    altloc: str
    resname: str
    chain_id: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str


def clean_pdb_lines(lines: Iterable[str]) -> list[str]:
    """
    Return parser-friendly PDB lines for docking output.

    The cleaner strips HDOCK wrapper records, preserves coordinate and TER
    records, drops malformed coordinate rows, normalizes blank chain IDs to X,
    and emits a single END record when coordinates were retained.
    """
    cleaned: list[str] = []
    saw_coordinates = False
    hdock_section: str | None = None

    for raw_line in lines:
        line = raw_line.rstrip("\r\n")
        record = line[:6].strip().upper()

        if record == "HEADER":
            header_text = line.lower()
            if "lig" in header_text:
                hdock_section = "ligand"
            elif "rec" in header_text:
                hdock_section = "receptor"
            continue
        if record in SKIPPED_RECORDS:
            continue
        if record == "END":
            continue
        if record == "TER":
            cleaned.append(f"{line}\n")
            continue
        if record not in COORDINATE_RECORDS:
            continue
        normalized_line = _normalize_coordinate_line(
            line,
            chain_override="D" if hdock_section == "ligand" else None,
        )
        if normalized_line is None:
            continue

        saw_coordinates = True
        cleaned.append(f"{normalized_line}\n")

    if saw_coordinates:
        cleaned.append("END\n")
    return cleaned


def clean_pdb_file(input_path: str, output_path: str) -> dict[str, Any]:
    """
    Clean one PDB file and write the parser-friendly result.

    Returns a small summary dictionary with input/output paths and retained line
    counts so callers can log preprocessing provenance.
    """
    source = Path(input_path)
    destination = Path(output_path)
    if not source.is_file():
        raise FileNotFoundError(f"PDB file not found: {input_path}")

    original_lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    cleaned_lines = clean_pdb_lines(original_lines)

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(cleaned_lines), encoding="utf-8")

    return {
        "input_path": str(source),
        "output_path": str(destination),
        "input_lines": len(original_lines),
        "output_lines": len(cleaned_lines),
        "coordinate_lines": sum(
            1 for line in cleaned_lines if line[:6].strip().upper() in COORDINATE_RECORDS
        ),
    }


def clean_pdb_directory(input_dir: str, output_dir: str) -> list[dict[str, Any]]:
    """
    Clean all PDB files in a directory using natural model-number ordering.
    """
    from loxdockaudit.sorting import natural_sort_key

    source_dir = Path(input_dir)
    destination_dir = Path(output_dir)
    if not source_dir.is_dir():
        raise FileNotFoundError(f"PDB directory not found: {input_dir}")

    summaries: list[dict[str, Any]] = []
    for pdb_path in sorted(source_dir.glob("*.pdb"), key=lambda path: natural_sort_key(str(path))):
        summaries.append(
            clean_pdb_file(
                str(pdb_path),
                str(destination_dir / pdb_path.name),
            )
        )
    return summaries


def _normalize_coordinate_line(line: str, chain_override: str | None = None) -> str | None:
    if len(line) >= 54:
        parsed = _parse_fixed_width_coordinate_line(line)
    else:
        parsed = _parse_split_coordinate_line(line)
    if parsed is None:
        return None

    element = parsed["element"] or _infer_element(str(parsed["atom_name"]))
    chain_id = chain_override or str(parsed["chain_id"])
    return (
        f"{parsed['record']:<6}{int(parsed['serial']):>5} "
        f"{str(parsed['atom_name']):<4}"
        f"{str(parsed['altloc'])[:1]:1}"
        f"{str(parsed['resname']):>3} "
        f"{chain_id[:1]:1}"
        f"{int(parsed['resseq']):>4}"
        f"{str(parsed['icode'])[:1]:1}"
        f"   {float(parsed['x']):>8.3f}"
        f"{float(parsed['y']):>8.3f}"
        f"{float(parsed['z']):>8.3f}"
        f"{float(parsed['occupancy']):>6.2f}"
        f"{float(parsed['bfactor']):>6.2f}"
        f"          {element:>2}"
    )


def _parse_fixed_width_coordinate_line(line: str) -> _ParsedCoordinate | None:
    try:
        return {
            "record": line[:6].strip().upper(),
            "serial": int(line[6:11]),
            "atom_name": line[12:16].strip(),
            "altloc": line[16:17].strip(),
            "resname": line[17:20].strip(),
            "chain_id": line[21:22].strip() or "X",
            "resseq": int(line[22:26]),
            "icode": line[26:27].strip(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "occupancy": _safe_float(line[54:60], 1.0),
            "bfactor": _safe_float(line[60:66], 0.0),
            "element": line[76:78].strip() if len(line) >= 78 else "",
        }
    except (ValueError, IndexError):
        return _parse_split_coordinate_line(line)


def _parse_split_coordinate_line(line: str) -> _ParsedCoordinate | None:
    parts = line.split()
    if len(parts) < 9:
        return None
    try:
        return {
            "record": parts[0].upper(),
            "serial": int(parts[1]),
            "atom_name": parts[2],
            "altloc": "",
            "resname": parts[3],
            "chain_id": parts[4] if len(parts[4]) == 1 else "X",
            "resseq": int(parts[5] if len(parts[4]) == 1 else parts[4]),
            "icode": "",
            "x": float(parts[6] if len(parts[4]) == 1 else parts[5]),
            "y": float(parts[7] if len(parts[4]) == 1 else parts[6]),
            "z": float(parts[8] if len(parts[4]) == 1 else parts[7]),
            "occupancy": 1.0,
            "bfactor": 0.0,
            "element": "",
        }
    except (ValueError, IndexError):
        return None


def _safe_float(value: str, default: float) -> float:
    try:
        return float(value)
    except ValueError:
        return default


def _infer_element(atom_name: str) -> str:
    letters = "".join(character for character in atom_name if character.isalpha())
    if not letters:
        return ""
    if len(letters) >= 2 and letters[:2].title() in {"FE", "ZN", "CU", "MG", "CA", "MN"}:
        return letters[:2].title()
    return letters[0].upper()
