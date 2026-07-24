# Tests

This directory contains regression tests for reconstructing Berry curvature from the irreducible Brillouin zone to the full Brillouin zone.

## Input

The `input/` directory contains:

* `OUTCAR`
* `IBZKPT`
* `bcurv_ij.dat`, using the current BerryCPT output format
* `bcurv_ij_legacy.dat`, using the legacy BerryCPT output format

## Expected output

The `expected_output/` directory contains the reconstructed full-BZ data:

```text
reconstructed_bcurv_fbz.dat
reconstructed_bcurv_fbz_legacy.dat
```

It also contains the corresponding combined standard-output and standard-error logs:

```text
reconstruction.log
reconstruction_legacy.log
```

The logs include the parsed symmetry operations, reconstruction statistics, full-BZ sums and averages, consistency warnings, and output-file confirmation.

## Reproducing the reference results

Run from the parent utility directory:

```bash
python berrycpt_reconstruct.py \
    tests/input/OUTCAR \
    tests/input/IBZKPT \
    tests/input/bcurv_ij.dat \
    --output tests/expected_output/reconstructed_bcurv_fbz.dat \
    > tests/expected_output/reconstruction.log 2>&1
```

For the legacy BerryCPT format:

```bash
python berrycpt_reconstruct.py \
    tests/input/OUTCAR \
    tests/input/IBZKPT \
    tests/input/bcurv_ij_legacy.dat \
    --output tests/expected_output/reconstructed_bcurv_fbz_legacy.dat \
    > tests/expected_output/reconstruction_legacy.log 2>&1
```

## Unit tests

Run:

```bash
python -m unittest tests/test_reconstruction.py -v
```

The unit tests execute the reconstruction for both BerryCPT formats and compare:

* the generated full-BZ data with the corresponding reference DAT file
* the combined standard-output and standard-error stream with the corresponding reference log

The output-file path in the final log line is normalized because the unit test writes the generated DAT file to a temporary directory.
