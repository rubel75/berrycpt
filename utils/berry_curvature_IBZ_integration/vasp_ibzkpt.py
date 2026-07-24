"""Read VASP IBZKPT files.

The parser reads the irreducible k-point coordinates and their integer
multiplicities.  It intentionally reads exactly the number of k-point rows
specified in the header, so optional tetrahedron information following the
k-point list does not cause a parsing failure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import math


class IBZKPTError(ValueError):
    """Raised when an IBZKPT file is incomplete or inconsistent."""


@dataclass(frozen=True)
class KPoint:
    """One irreducible k-point and its full-BZ multiplicity."""

    kfrac: tuple[float, float, float]
    mult: int


@dataclass
class IBZ:
    """Irreducible k-point set read from an IBZKPT file."""

    comment: str = ""
    kpoints: list[KPoint] = field(default_factory=list)
    nktot: int = 0
    multktot: int = 0

    def __getitem__(self, index: int) -> KPoint:
        return self.kpoints[index]

    def __len__(self) -> int:
        return len(self.kpoints)

    def __iter__(self) -> Iterator[KPoint]:
        return iter(self.kpoints)


def _next_nonempty(lines: list[str], start: int) -> tuple[int, str]:
    for index in range(start, len(lines)):
        text = lines[index].strip()
        if text:
            return index, text
    raise IBZKPTError("unexpected end of file")


def _parse_positive_integer_like(token: str, context: str) -> int:
    try:
        value = float(token)
    except ValueError as exc:
        raise IBZKPTError(f"{context}: invalid numeric value {token!r}") from exc

    if not math.isfinite(value):
        raise IBZKPTError(f"{context}: value is not finite")

    integer = int(round(value))
    if integer <= 0 or abs(value - integer) > 1.0e-8:
        raise IBZKPTError(
            f"{context}: expected a positive integer-like value, got {token!r}"
        )
    return integer


def read_ibzkpt(filename: str | Path) -> IBZ:
    """Read irreducible k-points and multiplicities from ``IBZKPT``.

    Blank lines are ignored.  Exactly ``nktot`` k-point rows are consumed.
    Additional rows, such as VASP tetrahedron data, are left uninterpreted.
    """

    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"IBZKPT file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise IBZKPTError("IBZKPT file is empty")

    comment_index, comment = _next_nonempty(lines, 0)
    count_index, count_text = _next_nonempty(lines, comment_index + 1)
    mode_index, mode_text = _next_nonempty(lines, count_index + 1)

    nktot = _parse_positive_integer_like(count_text.split()[0], "k-point count")

    if not mode_text.lower().startswith("reciprocal"):
        raise IBZKPTError(
            "IBZKPT must use reciprocal fractional coordinates, "
            f"found header {mode_text!r}"
        )

    kpoints: list[KPoint] = []
    cursor = mode_index + 1
    while len(kpoints) < nktot:
        try:
            line_index, text = _next_nonempty(lines, cursor)
        except IBZKPTError as exc:
            raise IBZKPTError(
                f"expected {nktot} k-point rows, found only {len(kpoints)}"
            ) from exc
        cursor = line_index + 1

        fields = text.split()
        if len(fields) < 4:
            raise IBZKPTError(
                f"line {line_index + 1}: expected three coordinates and a weight"
            )

        try:
            coordinates = tuple(float(token) for token in fields[:3])
        except ValueError as exc:
            raise IBZKPTError(
                f"line {line_index + 1}: invalid reciprocal coordinate"
            ) from exc

        if not all(math.isfinite(value) for value in coordinates):
            raise IBZKPTError(
                f"line {line_index + 1}: reciprocal coordinates must be finite"
            )

        multiplicity = _parse_positive_integer_like(
            fields[3], f"line {line_index + 1} multiplicity"
        )
        kpoints.append(
            KPoint(
                kfrac=(coordinates[0], coordinates[1], coordinates[2]),
                mult=multiplicity,
            )
        )

    return IBZ(
        comment=comment,
        kpoints=kpoints,
        nktot=nktot,
        multktot=sum(point.mult for point in kpoints),
    )


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Read VASP IBZKPT coordinates and multiplicities."
    )
    parser.add_argument("ibzkpt", help="path to the VASP IBZKPT file")
    args = parser.parse_args()

    ibz = read_ibzkpt(args.ibzkpt)
    print(f"Comment: {ibz.comment}")
    print(f"Irreducible k-points: {ibz.nktot}")
    print(f"Total multiplicity: {ibz.multktot}")
    print(f"First point: {ibz[0].kfrac}, multiplicity={ibz[0].mult}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
