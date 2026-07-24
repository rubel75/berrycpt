"""
vasp_symmetry.py

Module for parsing crystal symmetry operations from a VASP OUTCAR file.

Each symmetry operation includes:
    - R : 3×3 integer rotation matrix (direct-space)
    - tau_g : 3-element translation vector ("grid translation")
    - tau_p : 3-element translation vector ("non-primitive translation")

Only proper/unimodular rotations (det = ±1) are kept.

Usage:
    from vasp_symmetry import parse_outcar_symmetry
    sym_ops = parse_outcar_symmetry("OUTCAR")

    for op in sym_ops:
        print("R =")
        print(op.R)
        print("tau_p =", op.tau_p)
        print("tau_g =", op.tau_g)
"""

import re
import numpy as np
from typing import List, Optional
from dataclasses import dataclass

# --- Data structure ---------------------------------------------------------

@dataclass
class SymOp:
    """Container for one symmetry operation."""
    R: np.ndarray        # 3×3 integer rotation matrix (direct-space)
    tau_g: np.ndarray    # grid translation vector (3 floats)
    tau_p: np.ndarray    # non-primitive translation vector (3 floats)


# --- Internal helpers -------------------------------------------------------

# Regular expressions for parsing
_INT3_LINE = re.compile(r"(-?\d+)\s+(-?\d+)\s+(-?\d+)")
_FLOAT_TOKEN = re.compile(r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+")


def _parse_three_ints(line: str) -> Optional[List[int]]:
    """Extract three integers from a line."""
    m = _INT3_LINE.search(line)
    if not m:
        return None
    return [int(m.group(1)), int(m.group(2)), int(m.group(3))]


def _parse_three_floats_on_line(line: str) -> Optional[List[float]]:
    """Extract the last three floats on a line (VASP sometimes prints extra text)."""
    toks = _FLOAT_TOKEN.findall(line)
    if len(toks) < 3:
        return None
    vals = [float(t.replace('D', 'E').replace('d', 'e')) for t in toks[-3:]]
    return vals


# --- Main parser ------------------------------------------------------------

def parse_outcar_symmetry(path: str) -> List[SymOp]:
    """
    Parse all symmetry operations printed in blocks starting with
    'irot :' and 'isymop:' in a VASP OUTCAR file.

    Parameters
    ----------
    path : str
        Path to the OUTCAR file.

    Returns
    -------
    List[SymOp]
        List of symmetry operations (each with R, tau_g, tau_p).
    """
    ops: List[SymOp] = []
    number_of_operations = None

    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    i, n = 0, len(lines)
    while i < n:
        line = lines[i]

        # --- Find the line reporting total number of symmetry operations
        if "Subroutine INISYM returns:" in line:
            words = line.split()
            try:
                number_of_operations = int(words[4])
            except (IndexError, ValueError):
                raise ValueError(
                    f"Cannot extract number of symmetry operations from line:\n{line}"
                )
            print(f"[info] Expected {number_of_operations} symmetry operations")

        # --- Find each symmetry block
        if re.search(r'^\s*irot\s*:\s*\d+', line, flags=re.IGNORECASE):
            words = line.split()
            try:
                indx_of_operation = int(words[2])
            except (IndexError, ValueError):
                raise ValueError(f"Cannot parse symmetry index from line:\n{line}")

            print(f"\n[info] Extracting symmetry operation {indx_of_operation}")

            # Locate the 'isymop:' header
            j = i + 1
            while j < n and 'isymop' not in lines[j].lower():
                j += 1
            if j >= n:
                raise ValueError('Premature end of file while searching for isymop block.')

            # --- Read rotation matrix (3×3 integers)
            R_rows: List[List[int]] = []

            # First line: may have "isymop:" prefix
            parts = lines[j].split(':', 1)
            row = _parse_three_ints(parts[1] if len(parts) > 1 else lines[j])
            if row:
                R_rows.append(row)
            j += 1

            # Remaining two lines
            while j < n and len(R_rows) < 3:
                row = _parse_three_ints(lines[j])
                if row is not None:
                    R_rows.append(row)
                j += 1

            if len(R_rows) != 3:
                raise ValueError(f"isymop block near line {i+1} did not have 3 rows.")

            R = np.array(R_rows, dtype=int)

            # --- Optional translations
            tau_g = np.zeros(3)
            tau_p = np.zeros(3)

            # gtrans
            while j < n and 'gtrans' not in lines[j].lower():
                j += 1
            if j < n:
                gvals = _parse_three_floats_on_line(lines[j])
                if gvals:
                    tau_g = np.array(gvals, dtype=float)
                j += 1

            # ptrans
            while j < n and 'ptrans' not in lines[j].lower():
                j += 1
            if j < n:
                pvals = _parse_three_floats_on_line(lines[j])
                if pvals:
                    tau_p = np.array(pvals, dtype=float)
                j += 1

            det = round(np.linalg.det(R))
            if det not in (+1, -1):
                raise ValueError(
                    f"Rotation matrix determinant {det} invalid for operation {indx_of_operation}"
                )

            ops.append(SymOp(R=R, tau_g=tau_g, tau_p=tau_p))

            print(f"  R =\n{R}")
            print(f"  non-primitive translation tau_p = {tau_p}")
            print(f"  grid translation tau_g         = {tau_g}")

            i = j
        else:
            i += 1

    print(f"\n[info] Parsed {len(ops)} symmetry operations from OUTCAR")

    if number_of_operations is not None and len(ops) != number_of_operations:
        raise ValueError(
            f"Mismatch: expected {number_of_operations}, extracted {len(ops)} operations."
        )

    return ops


# --- Example direct run -----------------------------------------------------

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python vasp_symmetry.py OUTCAR")
        sys.exit(1)

    sym_ops = parse_outcar_symmetry(sys.argv[1])
    print(f"\nTotal symmetry operations: {len(sym_ops)}")
