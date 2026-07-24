"""Symmetry-reconstruct Berry curvature from irreducible VASP k-points.

This is the fourth layer of the reconstruction workflow.  It combines

* validated VASP symmetry geometry,
* validated IBZKPT stars, and
* the occupation-weighted band-0 BerryCPT vector at each IBZ representative.

For every unique full-BZ star member, all symmetry operations that map the
representative onto that member are applied independently.  Their transformed
Berry-curvature vectors are compared before their arithmetic mean is used as
the symmetry-projected value for that member.

If several operations map the representative k point to the same target k
point, they differ by an element of the little co-group.  Exact symmetry
therefore requires all transformed Berry-curvature vectors to agree.  A
failure of this test is retained as a diagnostic rather than hidden by taking
only the first operation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from berrycpt_bcurv import (
    BerryCurvatureData,
    IBZBerryCurvature,
    pair_with_ibzkpt,
    read_bcurv_ij,
)
from ibz_stars import (
    IBZStarSet,
    KPointStar,
    StarMember,
    build_ibz_stars,
    load_vasp_ibz_stars,
    periodic_equivalent,
)
from symmetry_geometry import SymmetryGeometry
from vasp_ibzkpt import IBZ, read_ibzkpt


FloatVector = NDArray[np.float64]
FloatMatrix = NDArray[np.float64]
ConsistencyPolicy = Literal["warn", "error", "ignore"]


class BerryReconstructionError(ValueError):
    """Raised when Berry-curvature symmetry reconstruction is inconsistent."""


def _readonly_vector(value: ArrayLike, name: str) -> FloatVector:
    vector = np.asarray(value, dtype=np.float64)
    if vector.shape != (3,):
        raise BerryReconstructionError(
            f"{name} must have shape (3,), got {vector.shape}"
        )
    if not np.all(np.isfinite(vector)):
        raise BerryReconstructionError(f"{name} contains non-finite values")
    vector = vector.copy()
    vector.setflags(write=False)
    return vector


def _readonly_matrix(value: ArrayLike, name: str) -> FloatMatrix:
    matrix = np.asarray(value, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[1] != 3:
        raise BerryReconstructionError(
            f"{name} must have shape (n, 3), got {matrix.shape}"
        )
    if not np.all(np.isfinite(matrix)):
        raise BerryReconstructionError(f"{name} contains non-finite values")
    matrix = matrix.copy()
    matrix.setflags(write=False)
    return matrix


def _validate_tolerances(atol: float, rtol: float) -> None:
    if not np.isfinite(atol) or atol < 0.0:
        raise BerryReconstructionError(
            f"absolute consistency tolerance must be finite and nonnegative, got {atol}"
        )
    if not np.isfinite(rtol) or rtol < 0.0:
        raise BerryReconstructionError(
            f"relative consistency tolerance must be finite and nonnegative, got {rtol}"
        )


def _validate_policy(policy: str) -> ConsistencyPolicy:
    if policy not in {"warn", "error", "ignore"}:
        raise BerryReconstructionError(
            "consistency policy must be one of 'warn', 'error', or 'ignore'"
        )
    return policy  # type: ignore[return-value]


def _pairwise_spread(vectors: FloatMatrix) -> tuple[float, float, FloatVector]:
    """Return maximum pairwise norm, scale, and component ranges."""

    count = vectors.shape[0]
    if count < 1:
        raise BerryReconstructionError("cannot compare an empty vector set")

    component_range = np.ptp(vectors, axis=0)
    scale = float(np.max(np.linalg.norm(vectors, axis=1)))

    maximum = 0.0
    for first in range(count):
        for second in range(first + 1, count):
            distance = float(np.linalg.norm(vectors[first] - vectors[second]))
            maximum = max(maximum, distance)

    return maximum, scale, _readonly_vector(component_range, "component range")


@dataclass(frozen=True)
class ReconstructedStarMember:
    """Berry curvature reconstructed at one unique full-BZ k-point."""

    kfrac: FloatVector
    operation_indices: tuple[int, ...]
    candidate_omegas_xyz: FloatMatrix
    omega_xyz: FloatVector
    max_pairwise_difference: float
    relative_spread: float
    component_range_xyz: FloatVector
    consistency_threshold: float
    consistent: bool

    @property
    def operation_count(self) -> int:
        return len(self.operation_indices)


@dataclass(frozen=True)
class ReconstructedKPointStar:
    """One IBZ representative and all reconstructed full-BZ vectors."""

    ibz_index: int
    representative: FloatVector
    source_omega_xyz: FloatVector
    members: tuple[ReconstructedStarMember, ...]
    star_sum_xyz: FloatVector

    @property
    def size(self) -> int:
        return len(self.members)

    @property
    def inconsistent_member_count(self) -> int:
        return sum(not member.consistent for member in self.members)

    @property
    def duplicate_mapping_member_count(self) -> int:
        return sum(member.operation_count > 1 for member in self.members)

    @property
    def max_pairwise_difference(self) -> float:
        return max(
            (member.max_pairwise_difference for member in self.members),
            default=0.0,
        )


@dataclass(frozen=True)
class BerryCurvatureReconstruction:
    """Complete symmetry-reconstructed full-BZ Berry-curvature data."""

    stars: tuple[ReconstructedKPointStar, ...]
    total_multiplicity: int
    full_sum_xyz: FloatVector
    full_average_xyz: FloatVector
    consistency_atol: float
    consistency_rtol: float
    consistency_policy: ConsistencyPolicy

    @property
    def full_kpoint_count(self) -> int:
        return sum(star.size for star in self.stars)

    @property
    def duplicate_mapping_member_count(self) -> int:
        return sum(star.duplicate_mapping_member_count for star in self.stars)

    @property
    def inconsistent_member_count(self) -> int:
        return sum(star.inconsistent_member_count for star in self.stars)

    @property
    def consistent_member_count(self) -> int:
        return self.full_kpoint_count - self.inconsistent_member_count

    @property
    def max_pairwise_difference(self) -> float:
        return max(
            (star.max_pairwise_difference for star in self.stars),
            default=0.0,
        )

    @property
    def max_relative_spread(self) -> float:
        return max(
            (
                member.relative_spread
                for star in self.stars
                for member in star.members
            ),
            default=0.0,
        )

    def inconsistent_members(
        self,
    ) -> list[tuple[ReconstructedKPointStar, ReconstructedStarMember]]:
        issues = [
            (star, member)
            for star in self.stars
            for member in star.members
            if not member.consistent
        ]
        return sorted(
            issues,
            key=lambda item: item[1].max_pairwise_difference,
            reverse=True,
        )


def _operation_lookup(
    operations: Sequence[SymmetryGeometry],
) -> dict[int, SymmetryGeometry]:
    if not operations:
        raise BerryReconstructionError("no symmetry operations were supplied")

    lookup: dict[int, SymmetryGeometry] = {}
    for operation in operations:
        if operation.index in lookup:
            raise BerryReconstructionError(
                f"duplicate symmetry operation index {operation.index}"
            )
        lookup[operation.index] = operation
    return lookup


def _reconstruct_member(
    member: StarMember,
    source_omega_xyz: ArrayLike,
    operation_lookup: dict[int, SymmetryGeometry],
    consistency_atol: float,
    consistency_rtol: float,
) -> ReconstructedStarMember:
    if not member.operation_indices:
        raise BerryReconstructionError(
            f"star member {member.kfrac} has no generating symmetry operations"
        )

    source = _readonly_vector(source_omega_xyz, "source Berry curvature")
    candidates: list[FloatVector] = []
    for operation_index in member.operation_indices:
        try:
            operation = operation_lookup[operation_index]
        except KeyError as exc:
            raise BerryReconstructionError(
                f"star member references missing operation {operation_index}"
            ) from exc
        candidates.append(operation.map_berry_curvature(source))

    candidate_matrix = _readonly_matrix(
        np.stack(candidates, axis=0), "candidate Berry curvatures"
    )
    mean = _readonly_vector(
        np.mean(candidate_matrix, axis=0), "mean Berry curvature"
    )
    maximum, scale, component_range = _pairwise_spread(candidate_matrix)
    threshold = float(consistency_atol + consistency_rtol * scale)
    consistent = bool(maximum <= threshold)

    if scale == 0.0:
        relative_spread = 0.0 if maximum == 0.0 else float("inf")
    else:
        relative_spread = maximum / scale

    return ReconstructedStarMember(
        kfrac=_readonly_vector(member.kfrac, "star-member k-point"),
        operation_indices=tuple(member.operation_indices),
        candidate_omegas_xyz=candidate_matrix,
        omega_xyz=mean,
        max_pairwise_difference=maximum,
        relative_spread=relative_spread,
        component_range_xyz=component_range,
        consistency_threshold=threshold,
        consistent=consistent,
    )


def _format_inconsistency(
    star: ReconstructedKPointStar,
    member: ReconstructedStarMember,
) -> str:
    krep = " ".join(f"{value:.12g}" for value in star.representative)
    ktarget = " ".join(f"{value:.12g}" for value in member.kfrac)
    operations = ",".join(str(index) for index in member.operation_indices)
    return (
        f"IBZ {star.ibz_index}, representative ({krep}), target ({ktarget}), "
        f"operations [{operations}]: max pairwise difference "
        f"{member.max_pairwise_difference:.6e} exceeds threshold "
        f"{member.consistency_threshold:.6e}"
    )


def _validate_inputs(
    paired: IBZBerryCurvature,
    star_set: IBZStarSet,
    operations: Sequence[SymmetryGeometry],
    kpoint_tolerance: float,
) -> dict[int, SymmetryGeometry]:
    if kpoint_tolerance <= 0.0:
        raise BerryReconstructionError("k-point tolerance must be positive")

    operation_lookup = _operation_lookup(operations)
    if len(operation_lookup) != star_set.operation_count:
        raise BerryReconstructionError(
            f"star set was built with {star_set.operation_count} operations, but "
            f"{len(operation_lookup)} operations were supplied"
        )

    if len(star_set.stars) != len(paired.bcurv.records):
        raise BerryReconstructionError(
            f"star set contains {len(star_set.stars)} representatives, but "
            f"BerryCPT contains {len(paired.bcurv.records)} records"
        )
    if len(star_set.stars) != paired.ibz.nktot:
        raise BerryReconstructionError(
            f"star set contains {len(star_set.stars)} representatives, but "
            f"IBZKPT reports {paired.ibz.nktot}"
        )

    for expected_index, (point, star, record) in enumerate(
        zip(paired.ibz.kpoints, star_set.stars, paired.bcurv.records), start=1
    ):
        if star.ibz_index != expected_index:
            raise BerryReconstructionError(
                f"star position {expected_index} has IBZ index {star.ibz_index}"
            )
        if record.kp_index != expected_index:
            raise BerryReconstructionError(
                f"BerryCPT position {expected_index} has KP index {record.kp_index}"
            )
        if star.expected_multiplicity != point.mult:
            raise BerryReconstructionError(
                f"IBZ {expected_index}: star multiplicity "
                f"{star.expected_multiplicity} differs from IBZKPT weight {point.mult}"
            )
        if not periodic_equivalent(
            star.representative, point.kfrac, tolerance=kpoint_tolerance
        ):
            raise BerryReconstructionError(
                f"IBZ {expected_index}: star representative {star.representative} "
                f"does not match IBZKPT coordinate {point.kfrac}"
            )

    if star_set.total_multiplicity != paired.ibz.multktot:
        raise BerryReconstructionError(
            f"star-set multiplicity {star_set.total_multiplicity} differs from "
            f"IBZKPT total {paired.ibz.multktot}"
        )

    return operation_lookup


def reconstruct_berry_curvature(
    paired: IBZBerryCurvature,
    star_set: IBZStarSet,
    operations: Sequence[SymmetryGeometry],
    *,
    consistency_atol: float = 1.0e-8,
    consistency_rtol: float = 1.0e-6,
    consistency_policy: ConsistencyPolicy = "warn",
    kpoint_tolerance: float = 1.0e-8,
) -> BerryCurvatureReconstruction:
    """Reconstruct occupation-weighted Berry curvature throughout the FBZ.

    Every operation listed for a star member is tried independently.  The
    candidate vectors are compared with the criterion

        max_pairwise_norm <= atol + rtol * max_candidate_norm.

    Their arithmetic mean is then used as the symmetry-projected vector at the
    unique target k-point.  ``consistency_policy='error'`` aborts at the first
    failed comparison.  ``'warn'`` and ``'ignore'`` retain all diagnostics in
    the returned object, with warning output handled by the command-line layer.
    """

    _validate_tolerances(consistency_atol, consistency_rtol)
    policy = _validate_policy(consistency_policy)
    operation_lookup = _validate_inputs(
        paired, star_set, operations, kpoint_tolerance
    )

    reconstructed_stars: list[ReconstructedKPointStar] = []
    full_sum = np.zeros(3, dtype=np.float64)

    for star, record in zip(star_set.stars, paired.bcurv.records):
        source = _readonly_vector(
            record.total_omega_xyz,
            f"IBZ {star.ibz_index} band-0 Berry curvature",
        )
        members = tuple(
            _reconstruct_member(
                member=member,
                source_omega_xyz=source,
                operation_lookup=operation_lookup,
                consistency_atol=consistency_atol,
                consistency_rtol=consistency_rtol,
            )
            for member in star.members
        )
        star_sum = _readonly_vector(
            np.sum(np.stack([member.omega_xyz for member in members]), axis=0),
            f"IBZ {star.ibz_index} star sum",
        )
        reconstructed_star = ReconstructedKPointStar(
            ibz_index=star.ibz_index,
            representative=_readonly_vector(
                star.representative, "IBZ representative"
            ),
            source_omega_xyz=source,
            members=members,
            star_sum_xyz=star_sum,
        )

        if policy == "error":
            for member in members:
                if not member.consistent:
                    raise BerryReconstructionError(
                        _format_inconsistency(reconstructed_star, member)
                    )

        reconstructed_stars.append(reconstructed_star)
        full_sum += star_sum

    total_multiplicity = star_set.total_multiplicity
    full_average = full_sum / float(total_multiplicity)

    return BerryCurvatureReconstruction(
        stars=tuple(reconstructed_stars),
        total_multiplicity=total_multiplicity,
        full_sum_xyz=_readonly_vector(full_sum, "full-BZ Berry-curvature sum"),
        full_average_xyz=_readonly_vector(
            full_average, "full-BZ Berry-curvature average"
        ),
        consistency_atol=float(consistency_atol),
        consistency_rtol=float(consistency_rtol),
        consistency_policy=policy,
    )


def write_reconstructed_fbz(
    reconstruction: BerryCurvatureReconstruction,
    filename: str | Path,
) -> None:
    """Write one line per unique reconstructed full-BZ k-point."""

    path = Path(filename)
    with path.open("w", encoding="utf-8") as stream:
        stream.write(
            "# fbz_index ibz_index k1 k2 k3 Omega_x Omega_y Omega_z "
            "n_operations max_pairwise_norm relative_spread consistent operations\n"
        )
        fbz_index = 0
        for star in reconstruction.stars:
            for member in star.members:
                fbz_index += 1
                operations = ",".join(str(i) for i in member.operation_indices)
                stream.write(
                    f"{fbz_index:6d} {star.ibz_index:6d} "
                    f"{member.kfrac[0]: .15e} {member.kfrac[1]: .15e} "
                    f"{member.kfrac[2]: .15e} "
                    f"{member.omega_xyz[0]: .15e} "
                    f"{member.omega_xyz[1]: .15e} "
                    f"{member.omega_xyz[2]: .15e} "
                    f"{member.operation_count:3d} "
                    f"{member.max_pairwise_difference: .8e} "
                    f"{member.relative_spread: .8e} "
                    f"{int(member.consistent):1d} {operations}\n"
                )


def load_and_reconstruct(
    outcar: str | Path,
    ibzkpt: str | Path,
    bcurv: str | Path,
    *,
    geometry_tolerance: float = 1.0e-7,
    kpoint_tolerance: float = 1.0e-8,
    consistency_atol: float = 1.0e-8,
    consistency_rtol: float = 1.0e-6,
    consistency_policy: ConsistencyPolicy = "warn",
) -> tuple[
    IBZ,
    list[SymmetryGeometry],
    IBZStarSet,
    BerryCurvatureData,
    BerryCurvatureReconstruction,
]:
    """Read all matching files and perform the complete stage-4 reconstruction."""

    ibz_data, operations, star_set = load_vasp_ibz_stars(
        outcar=outcar,
        ibzkpt=ibzkpt,
        geometry_tolerance=geometry_tolerance,
        kpoint_tolerance=kpoint_tolerance,
    )
    bcurv_data = read_bcurv_ij(bcurv)
    paired = pair_with_ibzkpt(bcurv_data, ibz_data)
    reconstruction = reconstruct_berry_curvature(
        paired=paired,
        star_set=star_set,
        operations=operations,
        consistency_atol=consistency_atol,
        consistency_rtol=consistency_rtol,
        consistency_policy=consistency_policy,
        kpoint_tolerance=kpoint_tolerance,
    )
    return ibz_data, operations, star_set, bcurv_data, reconstruction


def _format_vector(vector: ArrayLike) -> str:
    values = np.asarray(vector, dtype=np.float64)
    return " ".join(f"{value: .12e}" for value in values)


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct occupation-weighted Berry curvature from IBZ points, "
            "testing every symmetry operation that reaches each unique star member."
        )
    )
    parser.add_argument("outcar", help="matching OUTCAR with NWRITE=3 symmetry data")
    parser.add_argument("ibzkpt", help="matching IBZKPT")
    parser.add_argument("bcurv", help="matching BerryCPT bcurv_ij.dat")
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
        "--consistency-atol",
        type=float,
        default=1.0e-8,
        help="absolute tolerance for duplicate-operation vectors, default: 1e-8",
    )
    parser.add_argument(
        "--consistency-rtol",
        type=float,
        default=1.0e-6,
        help="relative tolerance for duplicate-operation vectors, default: 1e-6",
    )
    parser.add_argument(
        "--consistency-policy",
        choices=("warn", "error", "ignore"),
        default="warn",
        help="action when duplicate-operation vectors disagree, default: warn",
    )
    parser.add_argument(
        "--show-inconsistent",
        type=int,
        default=10,
        metavar="N",
        help="print the N largest duplicate-operation disagreements, default: 10",
    )
    parser.add_argument(
        "--output",
        help="optional output file containing all reconstructed full-BZ points",
    )
    args = parser.parse_args()

    if args.show_inconsistent < 0:
        parser.error("--show-inconsistent must be nonnegative")

    ibz, operations, star_set, bcurv, result = load_and_reconstruct(
        outcar=args.outcar,
        ibzkpt=args.ibzkpt,
        bcurv=args.bcurv,
        geometry_tolerance=args.geometry_tolerance,
        kpoint_tolerance=args.kpoint_tolerance,
        consistency_atol=args.consistency_atol,
        consistency_rtol=args.consistency_rtol,
        consistency_policy=args.consistency_policy,
    )

    print(f"Validated symmetry operations: {len(operations)}")
    print(f"IBZ representatives: {ibz.nktot}")
    print(f"BerryCPT band-0 records: {len(bcurv.records)}")
    print(f"Reconstructed full-BZ points: {result.full_kpoint_count}")
    print(f"Total IBZKPT multiplicity: {result.total_multiplicity}")
    print(
        "Members reached by multiple operations: "
        f"{result.duplicate_mapping_member_count}"
    )
    print(
        "Members failing duplicate-operation consistency: "
        f"{result.inconsistent_member_count}"
    )
    print(
        "Maximum pairwise vector difference: "
        f"{result.max_pairwise_difference:.12e}"
    )
    print(f"Maximum relative spread: {result.max_relative_spread:.12e}")
    print("Symmetry-reconstructed full-BZ sum:")
    print(_format_vector(result.full_sum_xyz))
    print("Symmetry-reconstructed full-BZ average:")
    print(_format_vector(result.full_average_xyz))

    if args.consistency_policy == "warn" and result.inconsistent_member_count:
        print(
            "WARNING: duplicate-operation vectors disagree. Their arithmetic "
            "mean was used as the symmetry-projected value."
        )
        for rank, (star, member) in enumerate(
            result.inconsistent_members()[: args.show_inconsistent], start=1
        ):
            print(f"  {rank:3d}. {_format_inconsistency(star, member)}")
            for operation_index, candidate in zip(
                member.operation_indices, member.candidate_omegas_xyz
            ):
                print(
                    f"       operation {operation_index:3d}: "
                    f"{_format_vector(candidate)}"
                )
            print(f"       mean:          {_format_vector(member.omega_xyz)}")

    if args.output:
        write_reconstructed_fbz(result, args.output)
        print(f"Wrote reconstructed full-BZ data to {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
