"""
vasp_lattice.py

Module for reading real-space (direct) and reciprocal lattice vectors
from a VASP OUTCAR file. Only the last occurrence of the lattice section
is retained.

Usage:
    from vasp_lattice import read_lattice

    lat = read_lattice("OUTCAR")
    print("Direct lattice (A):")
    print(lat.A)
    print("Reciprocal lattice (B):")
    print(lat.B)
"""

import numpy as np
import re

class Lattice:
    """Container for lattice vectors."""
    def __init__(self, A=None, B=None):
        # Direct lattice vectors in real space (3x3 numpy array)
        self.A = A
        # Reciprocal lattice vectors (3x3 numpy array)
        # Note: these DO NOT include a factor of 2π
        self.B = B


def read_lattice(filename):
    """
    Read the last occurrence of real and reciprocal lattice vectors
    from a VASP OUTCAR file.

    Parameters
    ----------
    filename : str
        Path to OUTCAR file.

    Returns
    -------
    Lattice
        Object containing direct (A) and reciprocal (B) lattice vectors.

    Raises
    ------
    ValueError
        If the required section or numerical data cannot be found.
    """

    # Read the file lines
    with open(filename, "r") as f:
        lines = f.readlines()

    # Find all sections that start with the lattice header
    header = "direct lattice vectors"
    indices = [i for i, line in enumerate(lines) if header in line]

    if not indices:
        raise ValueError(f"No '{header}' section found in {filename}")

    # Use the last occurrence
    start = indices[-1]
    # The next three lines contain the vectors
    section = lines[start + 1 : start + 4]

    if len(section) < 3:
        raise ValueError("Incomplete lattice vector section found.")

    A = np.zeros((3, 3), dtype=float)
    B = np.zeros((3, 3), dtype=float)

    # Regular expression: six float numbers per line
    float_pattern = re.compile(r"[-+]?\d*\.\d+(?:[Ee][-+]?\d+)?")

    for i, line in enumerate(section):
        nums = float_pattern.findall(line)
        if len(nums) != 6:
            raise ValueError(f"Invalid lattice line format at line {start + i + 1}.")
        vals = [float(x) for x in nums]
        A[i, :] = vals[:3]
        B[i, :] = vals[3:]

    return Lattice(A=A, B=B)


# Example usage if run directly
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python vasp_lattice.py OUTCAR")
        sys.exit(1)

    lat = read_lattice(sys.argv[1])
    print("Direct lattice vectors (A):")
    print(lat.A)
    print("\nReciprocal lattice vectors (B) [without 2π factor]:")
    print(lat.B)
