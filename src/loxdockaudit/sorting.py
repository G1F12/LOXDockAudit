"""Sorting helpers for LOXDockAudit file handling."""

from __future__ import annotations

import re
from pathlib import Path


def natural_sort_key(path: str) -> list[int | str]:
    """
    Sort key that orders model_2 before model_10.

    Splits filename into text and integer chunks.
    """
    return [
        int(chunk) if chunk.isdigit() else chunk.lower()
        for chunk in re.split(r"(\d+)", Path(path).name)
    ]
