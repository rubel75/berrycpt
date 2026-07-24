
# Berry-curvature IBZ integration

This utility reconstructs the occupation-weighted Berry curvature calculated at VASP irreducible Brillouin-zone k-points onto the full Brillouin-zone mesh.

It uses the symmetry operations from `OUTCAR`, the irreducible k-point mesh from `IBZKPT`, and the band-0 occupation-weighted Berry curvature written by BerryCPT.

## Requirements

* Python 3.10 or newer
* NumPy

```bash
python -m pip install numpy
```

## Input files

The following files must come from the same calculation:

* `OUTCAR` containing the detailed VASP symmetry information generated with `NWRITE = 3`
* `IBZKPT` containing the irreducible k-points and their multiplicities
* BerryCPT `bcurv_ij.dat`

Both current and legacy BerryCPT output formats are supported.

Current format:

```text
# KP: <index> NBANDS_OUT: <number> NBANDS_TRANS: <number>
<band> <degenerate_block> <Omega_yz> <Omega_zx> <Omega_xy>
```

Legacy format:

```text
# KP: <index> NVB: <number> NEMAX: <number>
<band> <Omega_yz> <Omega_zx> <Omega_xy>
```

The reconstruction uses the final band-0 row, which contains the occupation-weighted Berry-curvature total.

The Cartesian axial-vector convention is

```text
Omega_x = Omega_yz
Omega_y = Omega_zx
Omega_z = Omega_xy
```

## Usage

Run from this directory:

```bash
python berrycpt_reconstruct.py OUTCAR IBZKPT bcurv_ij.dat \
    --output reconstructed_bcurv_fbz.dat
```

The program:

1. reads and validates the lattice and symmetry operations
2. reconstructs the full k-point star of every IBZ representative
3. transforms the Berry-curvature axial vector under each symmetry operation
4. compares symmetry operations that reach the same full-BZ point
5. averages equivalent symmetry paths
6. writes the reconstructed full-BZ Berry-curvature field

Use

```bash
python berrycpt_reconstruct.py --help
```

to list the available tolerances and consistency options.

## Output

The output file contains one record for each reconstructed full-BZ k-point:

```text
# fbz_index ibz_index k1 k2 k3 Omega_x Omega_y Omega_z n_operations max_pairwise_norm relative_spread consistent operations
```

The program also prints the unnormalised full-BZ sum and the simple k-point average of the reconstructed Berry curvature.

It does not convert the result to anomalous Hall conductivity.

## Important assumptions

* `OUTCAR`, `IBZKPT`, and `bcurv_ij.dat` must use the same structure, magnetic configuration, symmetry setting, and k-point mesh.
* BerryCPT records are paired with `IBZKPT` by k-point index.
* The current `OUTCAR` parser does not automatically identify antiunitary or primed magnetic symmetry operations.
