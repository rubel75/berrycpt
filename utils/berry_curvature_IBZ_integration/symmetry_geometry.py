"""Symmetry geometry for VASP k-point and Berry-curvature reconstruction.

Matrix conventions
------------------
VASP lattice vectors are stored as rows of ``A`` and reciprocal lattice
vectors without the factor 2*pi are stored as rows of ``B``. Therefore

    A @ B.T = I

Fractional direct coordinates and fractional reciprocal coordinates are
column vectors. A VASP direct-space symmetry matrix ``W`` acts as

    x' = W @ x + tau

The corresponding reciprocal-space transformation is

    k' = W**(-T) @ k

For row-wise lattice matrices, the Cartesian polar-vector rotation is

    Q = A.T @ W @ B

and the spatial transformation of an axial vector is

    Q_axial = det(Q) * Q

Berry curvature is an axial vector that is odd under time reversal. For an
antiunitary operation, both the reciprocal coordinate and Berry curvature
receive one additional minus sign.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray


IntMatrix = NDArray[np.int64]
FloatMatrix = NDArray[np.float64]
FloatVector = NDArray[np.float64]


class SymmetryGeometryError(ValueError):
    """Raised when a symmetry or lattice consistency check fails."""


def _as_float_matrix(value: ArrayLike, name: str) -> FloatMatrix:
    matrix = np.asarray(value, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise SymmetryGeometryError(f"{name} must have shape (3, 3), got {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise SymmetryGeometryError(f"{name} contains non-finite values")
    return matrix


def _as_float_vector(value: ArrayLike, name: str) -> FloatVector:
    vector = np.asarray(value, dtype=np.float64)
    if vector.shape != (3,):
        raise SymmetryGeometryError(f"{name} must have shape (3,), got {vector.shape}")
    if not np.all(np.isfinite(vector)):
        raise SymmetryGeometryError(f"{name} contains non-finite values")
    return vector


def _as_integer_matrix(value: ArrayLike, name: str = "matrix") -> IntMatrix:
    matrix_float = np.asarray(value, dtype=np.float64)
    if matrix_float.shape != (3, 3):
        raise SymmetryGeometryError(
            f"{name} must have shape (3, 3), got {matrix_float.shape}"
        )
    if not np.all(np.isfinite(matrix_float)):
        raise SymmetryGeometryError(f"{name} contains non-finite values")

    matrix_int = np.rint(matrix_float).astype(np.int64)
    if not np.array_equal(matrix_float, matrix_int.astype(np.float64)):
        raise SymmetryGeometryError(f"{name} must contain exact integers")
    return matrix_int


def determinant_3x3_integer(matrix: ArrayLike) -> int:
    """Return the exact integer determinant of a 3 by 3 integer matrix."""
    m = _as_integer_matrix(matrix)
    a, b, c = (int(x) for x in m[0])
    d, e, f = (int(x) for x in m[1])
    g, h, i = (int(x) for x in m[2])
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def inverse_unimodular_integer(matrix: ArrayLike) -> IntMatrix:
    """Return the exact integer inverse of a unimodular 3 by 3 matrix.

    A crystallographic fractional-coordinate rotation must have determinant
    +1 or -1. Its inverse is therefore integer-valued.
    """
    m = _as_integer_matrix(matrix, "direct-space rotation")
    det = determinant_3x3_integer(m)
    if det not in (-1, 1):
        raise SymmetryGeometryError(
            f"direct-space rotation must have determinant +1 or -1, got {det}"
        )

    a, b, c = (int(x) for x in m[0])
    d, e, f = (int(x) for x in m[1])
    g, h, i = (int(x) for x in m[2])

    cofactor = np.array(
        [
            [e * i - f * h, -(d * i - f * g), d * h - e * g],
            [-(b * i - c * h), a * i - c * g, -(a * h - b * g)],
            [b * f - c * e, -(a * f - c * d), a * e - b * d],
        ],
        dtype=np.int64,
    )
    inverse = det * cofactor.T

    identity = np.eye(3, dtype=np.int64)
    if not np.array_equal(m @ inverse, identity):
        raise SymmetryGeometryError("failed exact right-inverse validation")
    if not np.array_equal(inverse @ m, identity):
        raise SymmetryGeometryError("failed exact left-inverse validation")
    return inverse


def validate_lattice(
    direct_lattice: ArrayLike,
    reciprocal_lattice: ArrayLike,
    tolerance: float = 1.0e-8,
) -> tuple[FloatMatrix, FloatMatrix]:
    """Validate row-wise direct and reciprocal lattice matrices."""
    if tolerance <= 0.0:
        raise SymmetryGeometryError("tolerance must be positive")

    a = _as_float_matrix(direct_lattice, "direct lattice A")
    b = _as_float_matrix(reciprocal_lattice, "reciprocal lattice B")

    det_a = float(np.linalg.det(a))
    if abs(det_a) <= tolerance:
        raise SymmetryGeometryError("direct lattice is singular or nearly singular")

    identity = np.eye(3, dtype=np.float64)
    product_error = float(np.max(np.abs(a @ b.T - identity)))
    if product_error > tolerance:
        raise SymmetryGeometryError(
            "direct and reciprocal lattices are inconsistent: "
            f"max|A @ B.T - I| = {product_error:.3e}"
        )

    expected_b = np.linalg.inv(a).T
    reciprocal_error = float(np.max(np.abs(b - expected_b)))
    if reciprocal_error > tolerance:
        raise SymmetryGeometryError(
            "reciprocal lattice does not equal inv(A).T: "
            f"maximum error = {reciprocal_error:.3e}"
        )

    return a, b


@dataclass
class SymmetryGeometry:
    """One validated VASP symmetry operation in several representations."""

    index: int
    direct_rotation: IntMatrix
    reciprocal_rotation: IntMatrix
    cartesian_rotation: FloatMatrix
    axial_rotation: FloatMatrix
    grid_translation: FloatVector
    primitive_translation: FloatVector
    antiunitary: bool = False

    @property
    def determinant(self) -> int:
        return determinant_3x3_integer(self.direct_rotation)

    def map_direct_fractional(self, coordinate: ArrayLike) -> FloatVector:
        """Map a direct fractional coordinate, including VASP translations."""
        x = _as_float_vector(coordinate, "direct fractional coordinate")
        return (
            self.direct_rotation @ x
            + self.grid_translation
            + self.primitive_translation
        )

    def map_k_fractional(self, coordinate: ArrayLike) -> FloatVector:
        """Map a reciprocal fractional coordinate.

        Antiunitary operations include the k -> -k action of time reversal.
        """
        k = _as_float_vector(coordinate, "reciprocal fractional coordinate")
        sign = -1.0 if self.antiunitary else 1.0
        return sign * (self.reciprocal_rotation @ k)

    def map_polar_vector(self, vector: ArrayLike) -> FloatVector:
        """Map a Cartesian polar vector under the spatial operation."""
        v = _as_float_vector(vector, "Cartesian polar vector")
        return self.cartesian_rotation @ v

    def map_axial_vector(self, vector: ArrayLike) -> FloatVector:
        """Map a Cartesian axial vector under the spatial operation only."""
        v = _as_float_vector(vector, "Cartesian axial vector")
        return self.axial_rotation @ v

    def map_berry_curvature(self, vector: ArrayLike) -> FloatVector:
        """Map Cartesian Berry curvature, including time reversal if present."""
        omega = _as_float_vector(vector, "Berry-curvature vector")
        sign = -1.0 if self.antiunitary else 1.0
        return sign * (self.axial_rotation @ omega)


def build_symmetry_geometry(
    index: int,
    direct_rotation: ArrayLike,
    direct_lattice: ArrayLike,
    reciprocal_lattice: ArrayLike,
    grid_translation: ArrayLike = (0.0, 0.0, 0.0),
    primitive_translation: ArrayLike = (0.0, 0.0, 0.0),
    antiunitary: bool = False,
    tolerance: float = 1.0e-8,
) -> SymmetryGeometry:
    """Construct and validate one symmetry operation."""
    if index < 1:
        raise SymmetryGeometryError(f"symmetry index must be positive, got {index}")

    a, b = validate_lattice(direct_lattice, reciprocal_lattice, tolerance)
    w = _as_integer_matrix(direct_rotation, "direct-space rotation W")
    det_w = determinant_3x3_integer(w)
    if det_w not in (-1, 1):
        raise SymmetryGeometryError(
            f"operation {index}: det(W) must be +1 or -1, got {det_w}"
        )

    w_inverse = inverse_unimodular_integer(w)
    reciprocal_rotation = w_inverse.T.copy()

    cartesian_rotation = a.T @ w @ b
    det_q = float(np.linalg.det(cartesian_rotation))
    orthogonality_error = float(
        np.max(np.abs(cartesian_rotation.T @ cartesian_rotation - np.eye(3)))
    )
    determinant_error = abs(det_q - float(det_w))

    if orthogonality_error > tolerance:
        raise SymmetryGeometryError(
            f"operation {index}: Cartesian rotation is not orthogonal, "
            f"max|Q.T @ Q - I| = {orthogonality_error:.3e}"
        )
    if determinant_error > tolerance:
        raise SymmetryGeometryError(
            f"operation {index}: det(Q) = {det_q:.12g} disagrees with "
            f"det(W) = {det_w}"
        )

    direct_consistency_error = float(
        np.max(np.abs(a.T @ w - cartesian_rotation @ a.T))
    )
    if direct_consistency_error > tolerance:
        raise SymmetryGeometryError(
            f"operation {index}: direct Cartesian transformation is inconsistent, "
            f"maximum error = {direct_consistency_error:.3e}"
        )

    reciprocal_consistency_error = float(
        np.max(
            np.abs(
                b.T @ reciprocal_rotation
                - cartesian_rotation @ b.T
            )
        )
    )
    if reciprocal_consistency_error > tolerance:
        raise SymmetryGeometryError(
            f"operation {index}: reciprocal Cartesian transformation is inconsistent, "
            f"maximum error = {reciprocal_consistency_error:.3e}"
        )

    axial_rotation = float(det_w) * cartesian_rotation

    return SymmetryGeometry(
        index=index,
        direct_rotation=w.copy(),
        reciprocal_rotation=reciprocal_rotation,
        cartesian_rotation=cartesian_rotation,
        axial_rotation=axial_rotation,
        grid_translation=_as_float_vector(grid_translation, "grid translation"),
        primitive_translation=_as_float_vector(
            primitive_translation, "primitive translation"
        ),
        antiunitary=bool(antiunitary),
    )


def build_symmetry_geometries(
    parsed_operations: Iterable[object],
    direct_lattice: ArrayLike,
    reciprocal_lattice: ArrayLike,
    antiunitary_flags: Sequence[bool] | None = None,
    tolerance: float = 1.0e-8,
) -> list[SymmetryGeometry]:
    """Convert operations returned by ``vasp_symmetry.py`` into validated objects."""
    operations = list(parsed_operations)
    if not operations:
        raise SymmetryGeometryError("no symmetry operations were supplied")

    if antiunitary_flags is None:
        flags = [False] * len(operations)
    else:
        flags = list(antiunitary_flags)
        if len(flags) != len(operations):
            raise SymmetryGeometryError(
                "number of antiunitary flags does not match number of operations"
            )

    geometries: list[SymmetryGeometry] = []
    for position, (operation, antiunitary) in enumerate(zip(operations, flags), start=1):
        if not hasattr(operation, "R"):
            raise SymmetryGeometryError(
                f"operation {position} has no direct-space rotation attribute R"
            )

        index = int(getattr(operation, "index", position))
        tau_g = getattr(operation, "tau_g", np.zeros(3))
        tau_p = getattr(operation, "tau_p", np.zeros(3))

        geometries.append(
            build_symmetry_geometry(
                index=index,
                direct_rotation=operation.R,
                direct_lattice=direct_lattice,
                reciprocal_lattice=reciprocal_lattice,
                grid_translation=tau_g,
                primitive_translation=tau_p,
                antiunitary=antiunitary,
                tolerance=tolerance,
            )
        )

    return geometries


def load_vasp_symmetry_geometry(
    outcar: str | Path,
    antiunitary_flags: Sequence[bool] | None = None,
    tolerance: float = 1.0e-8,
) -> tuple[FloatMatrix, FloatMatrix, list[SymmetryGeometry]]:
    """Read one OUTCAR using the existing parsers and validate its symmetry geometry."""
    from vasp_lattice import read_lattice
    from vasp_symmetry import parse_outcar_symmetry

    outcar_path = Path(outcar)
    if not outcar_path.is_file():
        raise FileNotFoundError(f"OUTCAR not found: {outcar_path}")

    lattice = read_lattice(str(outcar_path))
    parsed_operations = parse_outcar_symmetry(str(outcar_path))
    a, b = validate_lattice(lattice.A, lattice.B, tolerance)
    geometries = build_symmetry_geometries(
        parsed_operations,
        direct_lattice=a,
        reciprocal_lattice=b,
        antiunitary_flags=antiunitary_flags,
        tolerance=tolerance,
    )
    return a, b, geometries


def axial_group_projector(operations: Sequence[SymmetryGeometry]) -> FloatMatrix:
    """Return the average Berry-curvature transformation over a group."""
    if not operations:
        raise SymmetryGeometryError("cannot construct a projector from an empty group")
    matrices = []
    for operation in operations:
        sign = -1.0 if operation.antiunitary else 1.0
        matrices.append(sign * operation.axial_rotation)
    return np.mean(np.stack(matrices, axis=0), axis=0)


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate VASP direct, reciprocal, Cartesian, and axial symmetry matrices."
    )
    parser.add_argument("outcar", help="VASP OUTCAR containing NWRITE=3 symmetry output")
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-8,
        help="absolute matrix-validation tolerance, default: 1e-8",
    )
    parser.add_argument(
        "--show-matrices",
        action="store_true",
        help="print W, K, Q, and det(Q) Q for every operation",
    )
    args = parser.parse_args()

    a, b, operations = load_vasp_symmetry_geometry(
        args.outcar, tolerance=args.tolerance
    )

    print(f"Validated {len(operations)} symmetry operations")
    print(f"max|A @ B.T - I| = {np.max(np.abs(a @ b.T - np.eye(3))):.3e}")

    for operation in operations:
        q = operation.cartesian_rotation
        error = np.max(np.abs(q.T @ q - np.eye(3)))
        print(
            f"operation {operation.index:3d}: det={operation.determinant:+d} "
            f"orthogonality_error={error:.3e}"
        )
        if args.show_matrices:
            print("W =")
            print(operation.direct_rotation)
            print("K = W^(-T) =")
            print(operation.reciprocal_rotation)
            print("Q =")
            print(operation.cartesian_rotation)
            print("det(Q) Q =")
            print(operation.axial_rotation)

    print("Berry-curvature group projector =")
    print(axial_group_projector(operations))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
