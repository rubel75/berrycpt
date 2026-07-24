"""Construct full-Brillouin-zone k-point stars from VASP IBZ data.

This module is the second layer of the Berry-curvature reconstruction code.
It uses only reciprocal fractional coordinates and already validated
``SymmetryGeometry`` operations.  It does not transform Berry curvature yet.

For one irreducible representative k, every symmetry operation generates

    k' = K @ k

where ``K = W**(-T)``.  Coordinates differing by an integer reciprocal-lattice
vector are treated as identical.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from symmetry_geometry import SymmetryGeometry, load_vasp_symmetry_geometry
from vasp_ibzkpt import IBZ, KPoint, read_ibzkpt


FloatVector = NDArray[np.float64]
KPointKey = tuple[int, int, int]


class KPointStarError(ValueError):
    """Raised when symmetry stars do not reproduce VASP multiplicities."""


def _as_kpoint(value: ArrayLike, name: str = "k-point") -> FloatVector:
    coordinate = np.asarray(value, dtype=np.float64)
    if coordinate.shape != (3,):
        raise KPointStarError(f"{name} must have shape (3,), got {coordinate.shape}")
    if not np.all(np.isfinite(coordinate)):
        raise KPointStarError(f"{name} contains non-finite values")
    return coordinate


def canonicalize_kpoint(coordinate: ArrayLike, tolerance: float = 1.0e-8) -> FloatVector:
    """Return a reciprocal coordinate in the half-open cell [-1/2, 1/2).

    Values within ``tolerance`` of zero are set to exactly zero.  Values within
    ``tolerance`` of +1/2 or -1/2 are represented as exactly -1/2 so the two
    periodic boundary representations have one canonical form.
    """

    if tolerance <= 0.0:
        raise KPointStarError("k-point tolerance must be positive")

    kpoint = _as_kpoint(coordinate)
    wrapped = kpoint - np.floor(kpoint + 0.5)

    wrapped[np.abs(wrapped) <= tolerance] = 0.0
    boundary = (np.abs(wrapped - 0.5) <= tolerance) | (
        np.abs(wrapped + 0.5) <= tolerance
    )
    wrapped[boundary] = -0.5

    if np.any(wrapped < -0.5 - tolerance) or np.any(wrapped >= 0.5 + tolerance):
        raise KPointStarError(
            f"failed to canonicalize reciprocal coordinate {coordinate!r}"
        )
    return wrapped


def periodic_difference(first: ArrayLike, second: ArrayLike) -> FloatVector:
    """Return ``first - second`` reduced modulo reciprocal-lattice integers."""

    difference = _as_kpoint(first, "first k-point") - _as_kpoint(
        second, "second k-point"
    )
    return difference - np.rint(difference)


def periodic_equivalent(
    first: ArrayLike,
    second: ArrayLike,
    tolerance: float = 1.0e-8,
) -> bool:
    """Return whether two fractional reciprocal coordinates are equivalent."""

    if tolerance <= 0.0:
        raise KPointStarError("k-point tolerance must be positive")
    return bool(np.max(np.abs(periodic_difference(first, second))) <= tolerance)


def kpoint_key(coordinate: ArrayLike, tolerance: float = 1.0e-8) -> KPointKey:
    """Return a hashable periodic key for a reciprocal coordinate.

    The key is used only after canonicalization.  Any key collision is checked
    again with ``periodic_equivalent`` before points are merged.
    """

    canonical = canonicalize_kpoint(coordinate, tolerance)
    scaled = np.rint(canonical / tolerance)
    limit = np.iinfo(np.int64).max
    if np.any(np.abs(scaled) > limit):
        raise KPointStarError("k-point tolerance is too small for integer keying")
    values = scaled.astype(np.int64)
    return int(values[0]), int(values[1]), int(values[2])


@dataclass(frozen=True)
class StarMember:
    """One unique full-BZ point generated from an irreducible representative."""

    kfrac: FloatVector
    operation_indices: tuple[int, ...]

    @property
    def operation_count(self) -> int:
        return len(self.operation_indices)


@dataclass(frozen=True)
class KPointStar:
    """The complete symmetry star of one irreducible k-point."""

    ibz_index: int
    representative: FloatVector
    expected_multiplicity: int
    members: tuple[StarMember, ...]
    group_order: int

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def stabilizer_order(self) -> int:
        if self.size == 0 or self.group_order % self.size != 0:
            raise KPointStarError(
                f"IBZ point {self.ibz_index}: group order {self.group_order} "
                f"is not divisible by star size {self.size}"
            )
        return self.group_order // self.size


@dataclass(frozen=True)
class IBZStarSet:
    """Validated symmetry stars for all representatives in an IBZKPT file."""

    stars: tuple[KPointStar, ...]
    operation_count: int
    total_multiplicity: int
    unique_full_kpoints: int

    @property
    def star_size_distribution(self) -> dict[int, int]:
        return dict(sorted(Counter(star.size for star in self.stars).items()))

    @property
    def stabilizer_distribution(self) -> dict[int, int]:
        return dict(
            sorted(Counter(star.stabilizer_order for star in self.stars).items())
        )

    def full_kpoints(self) -> FloatVector:
        """Return all reconstructed full-BZ points in star order."""

        points = [member.kfrac for star in self.stars for member in star.members]
        if not points:
            return np.empty((0, 3), dtype=np.float64)
        return np.stack(points, axis=0)


def _validate_operations(operations: Sequence[SymmetryGeometry]) -> None:
    if not operations:
        raise KPointStarError("no symmetry operations were supplied")

    indices = [operation.index for operation in operations]
    if len(indices) != len(set(indices)):
        raise KPointStarError("symmetry operation indices must be unique")


def build_kpoint_star(
    ibz_index: int,
    representative: ArrayLike,
    expected_multiplicity: int,
    operations: Sequence[SymmetryGeometry],
    tolerance: float = 1.0e-8,
) -> KPointStar:
    """Generate and validate one irreducible k-point star."""

    if ibz_index < 1:
        raise KPointStarError(f"IBZ index must be positive, got {ibz_index}")
    if expected_multiplicity < 1:
        raise KPointStarError(
            f"IBZ point {ibz_index}: multiplicity must be positive"
        )
    _validate_operations(operations)

    representative_array = canonicalize_kpoint(representative, tolerance)
    grouped: dict[KPointKey, tuple[FloatVector, list[int]]] = {}

    for operation in operations:
        mapped = canonicalize_kpoint(
            operation.map_k_fractional(representative_array), tolerance
        )
        key = kpoint_key(mapped, tolerance)

        if key not in grouped:
            grouped[key] = (mapped, [operation.index])
            continue

        stored_coordinate, operation_indices = grouped[key]
        if not periodic_equivalent(stored_coordinate, mapped, tolerance):
            raise KPointStarError(
                f"IBZ point {ibz_index}: k-point key collision between "
                f"{stored_coordinate} and {mapped}. Reduce the tolerance."
            )
        operation_indices.append(operation.index)

    members = tuple(
        StarMember(
            kfrac=coordinate.copy(),
            operation_indices=tuple(sorted(operation_indices)),
        )
        for coordinate, operation_indices in sorted(
            grouped.values(), key=lambda item: tuple(item[0])
        )
    )

    star_size = len(members)
    if star_size != expected_multiplicity:
        raise KPointStarError(
            f"IBZ point {ibz_index} at {representative_array}: generated star "
            f"size {star_size}, but IBZKPT multiplicity is {expected_multiplicity}"
        )

    group_order = len(operations)
    if group_order % star_size != 0:
        raise KPointStarError(
            f"IBZ point {ibz_index}: group order {group_order} is not divisible "
            f"by star size {star_size}"
        )

    expected_stabilizer = group_order // star_size
    operation_counts = [member.operation_count for member in members]
    if any(count != expected_stabilizer for count in operation_counts):
        raise KPointStarError(
            f"IBZ point {ibz_index}: operations are not distributed uniformly "
            f"over the star. Expected {expected_stabilizer} operations per member, "
            f"found {operation_counts}"
        )

    if not any(
        periodic_equivalent(member.kfrac, representative_array, tolerance)
        for member in members
    ):
        raise KPointStarError(
            f"IBZ point {ibz_index}: representative is absent from its own star. "
            "The supplied operations may not contain the identity."
        )

    return KPointStar(
        ibz_index=ibz_index,
        representative=representative_array,
        expected_multiplicity=expected_multiplicity,
        members=members,
        group_order=group_order,
    )


def build_ibz_stars(
    ibz: IBZ,
    operations: Sequence[SymmetryGeometry],
    tolerance: float = 1.0e-8,
) -> IBZStarSet:
    """Build all stars and verify that they partition the reconstructed FBZ."""

    _validate_operations(operations)
    if ibz.nktot != len(ibz.kpoints):
        raise KPointStarError(
            f"IBZ object reports {ibz.nktot} points but contains "
            f"{len(ibz.kpoints)}"
        )
    if ibz.nktot < 1:
        raise KPointStarError("IBZ contains no k-points")

    stars: list[KPointStar] = []
    global_points: dict[KPointKey, tuple[int, FloatVector]] = {}

    for ibz_index, point in enumerate(ibz.kpoints, start=1):
        star = build_kpoint_star(
            ibz_index=ibz_index,
            representative=point.kfrac,
            expected_multiplicity=point.mult,
            operations=operations,
            tolerance=tolerance,
        )
        stars.append(star)

        for member in star.members:
            key = kpoint_key(member.kfrac, tolerance)
            if key not in global_points:
                global_points[key] = (ibz_index, member.kfrac)
                continue

            previous_index, previous_coordinate = global_points[key]
            if not periodic_equivalent(previous_coordinate, member.kfrac, tolerance):
                raise KPointStarError(
                    "global k-point key collision between non-equivalent points. "
                    "Reduce the tolerance."
                )
            raise KPointStarError(
                f"IBZ points {previous_index} and {ibz_index} generate the same "
                f"full-BZ point {member.kfrac}. The IBZ representatives are not "
                "symmetry-inequivalent under the supplied group."
            )

    total_multiplicity = sum(point.mult for point in ibz.kpoints)
    if total_multiplicity != ibz.multktot:
        raise KPointStarError(
            f"IBZ object reports total multiplicity {ibz.multktot}, but the "
            f"point-wise sum is {total_multiplicity}"
        )

    unique_full_kpoints = len(global_points)
    if unique_full_kpoints != total_multiplicity:
        raise KPointStarError(
            f"reconstructed {unique_full_kpoints} unique full-BZ points, but "
            f"the total IBZKPT multiplicity is {total_multiplicity}"
        )

    return IBZStarSet(
        stars=tuple(stars),
        operation_count=len(operations),
        total_multiplicity=total_multiplicity,
        unique_full_kpoints=unique_full_kpoints,
    )


def load_vasp_ibz_stars(
    outcar: str | Path,
    ibzkpt: str | Path,
    geometry_tolerance: float = 1.0e-7,
    kpoint_tolerance: float = 1.0e-8,
) -> tuple[IBZ, list[SymmetryGeometry], IBZStarSet]:
    """Read OUTCAR and IBZKPT, then construct and validate every star."""

    _, _, operations = load_vasp_symmetry_geometry(
        outcar, tolerance=geometry_tolerance
    )
    ibz = read_ibzkpt(ibzkpt)
    stars = build_ibz_stars(ibz, operations, tolerance=kpoint_tolerance)
    return ibz, operations, stars


def _format_distribution(distribution: dict[int, int]) -> str:
    return ", ".join(
        f"{key}: {count}" for key, count in distribution.items()
    )


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct and validate full-BZ k-point stars from VASP OUTCAR "
            "symmetry operations and IBZKPT multiplicities."
        )
    )
    parser.add_argument("outcar", help="OUTCAR containing NWRITE=3 symmetry output")
    parser.add_argument("ibzkpt", help="matching VASP IBZKPT file")
    parser.add_argument(
        "--geometry-tolerance",
        type=float,
        default=1.0e-7,
        help="OUTCAR lattice and Cartesian-matrix tolerance, default: 1e-7",
    )
    parser.add_argument(
        "--kpoint-tolerance",
        type=float,
        default=1.0e-8,
        help="periodic fractional-coordinate tolerance, default: 1e-8",
    )
    parser.add_argument(
        "--show-stars",
        action="store_true",
        help="print every representative, star size, and stabilizer order",
    )
    args = parser.parse_args()

    ibz, operations, star_set = load_vasp_ibz_stars(
        outcar=args.outcar,
        ibzkpt=args.ibzkpt,
        geometry_tolerance=args.geometry_tolerance,
        kpoint_tolerance=args.kpoint_tolerance,
    )

    print(f"Validated symmetry operations: {len(operations)}")
    print(f"Irreducible k-points: {ibz.nktot}")
    print(f"Total IBZKPT multiplicity: {star_set.total_multiplicity}")
    print(f"Unique reconstructed FBZ points: {star_set.unique_full_kpoints}")
    print(
        "Star-size distribution (size: number of IBZ points): "
        + _format_distribution(star_set.star_size_distribution)
    )
    print(
        "Stabilizer-order distribution (order: number of IBZ points): "
        + _format_distribution(star_set.stabilizer_distribution)
    )

    if args.show_stars:
        for star in star_set.stars:
            representative = " ".join(f"{value: .12f}" for value in star.representative)
            print(
                f"IBZ {star.ibz_index:4d}: k=({representative}) "
                f"star={star.size:2d} stabilizer={star.stabilizer_order:2d}"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
