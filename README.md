# BerryCPT

BerryCPT calculates band-resolved Berry curvature, orbital angular momentum (OAM), and generalized orbital-angular-momentum matrices from density-functional-theory matrix elements.

Supported electronic-structure codes:

- [WIEN2k](https://www.wien2k.at)
- [VASP](https://www.vasp.at)

BerryCPT is written in Fortran and is intended for Linux-based systems.

## Current status

The current development branch contains a major refactoring and extension of the code after the `v1.0.0` release.

Validated in the current branch:

- ordinary Berry curvature for MoS2 using WIEN2k
- ordinary OAM for MoS2 using WIEN2k
- explicit Fermi-energy occupations using `--efermiry`
- ordinary Berry curvature for MoS2 using VASP
- ordinary OAM for MoS2 using VASP
- explicit Fermi-energy occupations using `--efermiev`

Implemented and still undergoing broader validation:

- spin-resolved WIEN2k Berry-curvature quantities
- spin-up-projected OAM
- generalized OAM matrices
- additional occupation and spin configurations

Release `v1.0.0` preserves the pre-refactor implementation. That release contains a known pre-existing bug in the VASP `EIGENVAL` reader that limits imported eigenvalues to two decimal places. The current development version reads the full precision present in `EIGENVAL`.

## Installation

Clone the repository:

```bash
git clone https://github.com/rubel75/berrycpt
cd berrycpt
```

The supplied `Makefile` uses Intel MKL and requires the environment variable `MKLROOT` to be defined. Initialize the Intel oneAPI or MKL environment, then verify:

```bash
echo $MKLROOT
```

The default compiler configuration uses Intel `ifort` with debugging checks enabled. Alternative GNU Fortran settings are included as commented lines in the `Makefile`.

Compile with:

```bash
make
```

Useful build targets are:

```bash
make clean
make veryclean
```

`make clean` removes object and module files. `make veryclean` also removes the executable and editor backup files.

## Input preparation

BerryCPT requires matrix elements, eigenvalues, and either an explicit
Fermi energy or occupations calculated by the DFT code.

Instructions for generating the required matrix-element files are available
on the Wiki:

- [Generating WIEN2k `case.mommat2` files](https://github.com/rubel75/mstar/wiki/Generate-case.mommat2-file-in-WIEN2k)
- [Generating WIEN2k `case.mommat2` and VASP `WAVEDER` files](https://github.com/rubel75/mstar/wiki)

The matrix-element and eigenvalue files used in one BerryCPT calculation
must originate from the same DFT calculation.

The number of unoccupied bands must be convergence-tested. Contributions
from remote bands decrease as \(1/(\Delta E)^2\) for Berry curvature, but
only as \(1/\Delta E\) for OAM and generalized OAM. Consequently, OAM
normally requires a larger energy window than Berry curvature.

## Command-line syntax

General form:

```bash
berrycpt MATRIX_FILE [MATRIX_FILE_UP MATRIX_FILE_DN] \
    --enefile ENERGY_FILE \
    OCCUPATION_MODE
```

Specify either one matrix-element file or three matrix-element files. Three matrix-element files are supported only for WIEN2k calculations with total, spin-up, and spin-down projected matrix elements.

All calculations require `--enefile` and exactly one occupation mode.

For VASP, `OCCUPATION_MODE` is exactly one of:

```text
--efermiev EF
--efermiry EF
--dftocc
```

For WIEN2k, `OCCUPATION_MODE` is exactly one of:

```text
--efermiev EF
--efermiry EF
--dftocc --occfile OCCUPATION_FILE
```

The options have the following meanings:

- `--efermiev EF` constructs step-function occupations using a Fermi energy in eV
- `--efermiry EF` constructs step-function occupations using a Fermi energy in Ry
- `--dftocc` uses occupations supplied by the DFT calculation
- `--occfile FILE` specifies the WIEN2k occupation file and is required with `--dftocc`

Both explicit-Fermi-energy options are accepted for either input format. The examples below use eV for VASP and Ry for WIEN2k as convenient conventions, not as restrictions.

For VASP, DFT occupations are read directly from `EIGENVAL`, and `--occfile` must not be supplied. For WIEN2k, `--dftocc` is incomplete without `--occfile`.

Display the built-in help with:

```bash
berrycpt --help
```

## Execution examples

Each configuration below shows both occupation modes: an explicit Fermi energy and the occupations supplied by the DFT calculation. Choose one mode for each run.

### VASP

Using an explicit Fermi energy:

```bash
/path/to/berrycpt/berrycpt WAVEDER \
    --enefile EIGENVAL \
    --efermiev EF
```

Using DFT occupations from `EIGENVAL`:

```bash
/path/to/berrycpt/berrycpt WAVEDER \
    --enefile EIGENVAL \
    --dftocc
```

The primary VASP input file must have the basename `WAVEDER` or `WAVEDERF`. When `WAVEDERF` is supplied as the identifying name, BerryCPT reads the companion binary `WAVEDER` file from the same directory.

### WIEN2k without spin polarization or SOC

Using an explicit Fermi energy:

```bash
/path/to/berrycpt/berrycpt case.mommat2 \
    --enefile case.energy \
    --efermiry EF
```

Using DFT occupations:

```bash
/path/to/berrycpt/berrycpt case.mommat2 \
    --enefile case.energy \
    --dftocc \
    --occfile case.weight
```

### Spin-polarized WIEN2k without SOC

Run the spin channels separately using matching matrix-element, energy, and occupation files.

Using an explicit Fermi energy:

```bash
/path/to/berrycpt/berrycpt case.mommat2up \
    --enefile case.energyup \
    --efermiry EF

/path/to/berrycpt/berrycpt case.mommat2dn \
    --enefile case.energydn \
    --efermiry EF
```

Using DFT occupations:

```bash
/path/to/berrycpt/berrycpt case.mommat2up \
    --enefile case.energyup \
    --dftocc \
    --occfile case.weightup

/path/to/berrycpt/berrycpt case.mommat2dn \
    --enefile case.energydn \
    --dftocc \
    --occfile case.weightdn
```

### WIEN2k with SOC and without spin polarization

#### Ordinary quantities from the total matrix elements

Using an explicit Fermi energy:

```bash
/path/to/berrycpt/berrycpt case.mommat2 \
    --enefile case.energyso \
    --efermiry EF
```

Using DFT occupations:

```bash
/path/to/berrycpt/berrycpt case.mommat2 \
    --enefile case.energyso \
    --dftocc \
    --occfile case.weight
```

#### Ordinary and spin-resolved quantities

Supply the total, spin-up, and spin-down projected matrix elements in that order.

Using an explicit Fermi energy:

```bash
/path/to/berrycpt/berrycpt \
    case.mommat2 case.mommat2up case.mommat2dn \
    --enefile case.energyso \
    --efermiry EF
```

Using DFT occupations:

```bash
/path/to/berrycpt/berrycpt \
    case.mommat2 case.mommat2up case.mommat2dn \
    --enefile case.energyso \
    --dftocc \
    --occfile case.weight
```

A validated example using DFT occupations is:

```bash
../../berrycpt \
    CoSi-berrycpt.mommat2 \
    CoSi-berrycpt.mommat2up \
    CoSi-berrycpt.mommat2dn \
    --enefile CoSi-berrycpt.energyso \
    --dftocc \
    --occfile CoSi-berrycpt.weight
```

A spin-polarized SOC workflow has not yet been included among the validated cases documented here. Use matching spin-labelled WIEN2k files and verify the results independently before relying on that configuration.

## Degenerate bands

Bands are assigned to contiguous degenerate groups at every k-point and spin channel. Adjacent bands are grouped when their energy separation does not exceed `1.0E-5 Ha`.

For an isolated band, BerryCPT evaluates the usual band-resolved perturbation-theory expression. For a degenerate group, it constructs a Hermitian effective matrix and reports its eigenvalues. The reported values inside such a block therefore correspond to eigenvalues of the effective operator, not expectation values attached to the original unrotated DFT eigenvectors.

The target-band cutoff is adjusted so that it never divides a degenerate group. For WIEN2k, the initial target range is up to twice the largest occupied-band index, limited by the number of available bands. For VASP, the initial target range is read from `NBANDS_CDER` in `WAVEDER`.

## Output conventions

The pseudovector component order is:

```text
yz, zx, xy
```

These correspond to Voigt indices 4, 5, and 6.

Units:

- Berry curvature: `bohr^2`
- OAM: `hbar`
- generalized OAM: `hbar`

Existing output files are replaced. Every output file contains a detailed header describing its columns, k-point organization, band ranges, degenerate-block labels, and units. At normal termination, BerryCPT prints a summary of the output files that were created and closed successfully.

### Ordinary Berry curvature

Output file:

```text
bcurv_ij[-up/-dn].dat
```

The optional suffix identifies an ordinary spin channel, for example a collinear VASP spin channel or a separately processed WIEN2k `up` or `dn` input.

For an isolated band `n`, the generalized operator form is

```text
Omega_n^(A,B)(alpha,beta)
  = -2 Im sum_(l != n)
      A_alpha(n,l) B_beta(l,n) / (E_n - E_l)^2 .
```

Ordinary Berry curvature uses the total velocity or momentum matrix elements in both operator positions, `A = v` and `B = v`.

Each k-point block contains:

```text
# KP: k  NBANDS_OUT: Nout  NBANDS_TRANS: Ntrans
band  block  Omega_yz  Omega_zx  Omega_xy
```

The final record for each k-point has `band = 0` and `block = 0`. It contains the occupation-weighted total

```text
Omega_total = sum_n f_n Omega_n .
```

### Spin-resolved Berry-curvature quantities

These files are written only when three WIEN2k matrix-element files are supplied:

```text
bcurv_ij-up.dat
bcurv_ij-dn.dat
bcurv_ij-up-dn.dat
```

The implemented operator combinations are:

```text
bcurv_ij-up.dat:      A = v_up,          B = v
bcurv_ij-dn.dat:      A = v_dn,          B = v
bcurv_ij-up-dn.dat:   A = v_up - v_dn,   B = v
```

The third quantity is labelled the sigma-z-normalized spin Berry curvature.

For non-degenerate bands, linear combinations of these quantities follow directly from the operator definitions. Inside a degenerate block, the three effective matrices are diagonalized separately. Their ordered eigenvalues must not be combined band by band. Linear relations remain valid at the matrix level and for block traces.

### Ordinary OAM

Output file:

```text
oam_ij[-up/-dn].dat
```

For an isolated band `n`,

```text
L_n^(A,B)(alpha,beta)
  = 2 Im sum_(l != n)
      A_alpha(n,l) B_beta(l,n) / (E_n - E_l) .
```

Ordinary OAM uses `A = v` and `B = v`.

Each record contains:

```text
band  block  L_yz  L_zx  L_xy
```

No occupation-weighted OAM total is written.

### Spin-up-projected OAM

Output file:

```text
oam_ij-sigma_z-up.dat
```

This file is written only when the three spin-resolved WIEN2k matrix-element files are supplied. The current implementation uses the spin-up projected matrix elements in both operator positions:

```text
A = v_up
B = v_up
```

The output has the same band, block, and component structure as ordinary OAM.

### Generalized OAM

Output file:

```text
goam_ij_nm[-up/-dn].dat
```

For each component `c`, BerryCPT writes the complete complex Hermitian matrix

```text
L_c(i,j) = <u_i|L_c|u_j> .
```

The first matrix index is the bra-band index and the second is the ket-band index. The matrix remains in the original DFT eigenstate basis and is not diagonalized within degenerate blocks.

For Cartesian components `alpha` and `beta`, define

```text
P_ijl = A_alpha(i,l) B_beta(l,j)
      - B_beta(i,l) A_alpha(l,j) .
```

BerryCPT evaluates

```text
L_ij^(alpha,beta)
  = (i/2) sum_l P_ijl
      [ Q_i(l)/(E_l-E_i) + Q_j(l)/(E_l-E_j) ] .
```

`Q_i(l)` is zero when the intermediate state `l` belongs to the same degenerate block as external state `i`, and one otherwise. `Q_j(l)` is defined analogously. These factors exclude singular same-block intermediate-state couplings.

Only the upper triangle is calculated explicitly. The lower triangle is reconstructed by Hermitian conjugation. The full matrix is written using complex `(Re,Im)` fields.

The formulation was developed from Eq. (2) of Faria Junior et al., “Generalized many-body exciton g-factors: magnetic hybridization and non-monotonic Rydberg series in monolayer WSe2”:

- [arXiv:2505.18468](https://arxiv.org/abs/2505.18468)

## Numerical convergence

Berry curvature, OAM, and generalized OAM all depend on intermediate-state sums. For a remote intermediate state separated by an energy `Delta E`, the energy denominator suppresses its contribution as:

```text
Berry curvature:  1/(Delta E)^2
OAM:              1/(Delta E)
Generalized OAM:  1/(Delta E)
```

Berry curvature therefore generally converges faster with the number and energy range of unoccupied bands. OAM and generalized OAM normally require a larger intermediate-state window. This comparison concerns the denominator dependence. The matrix elements, near-degenerate states, and quality of the high-energy eigenstates also affect the observed convergence.

Convergence must be checked with respect to:

- the number of unoccupied bands
- the upper energy range used to generate matrix elements
- the quality of high-energy eigenstates
- k-point sampling
- the DFT basis and related numerical parameters

Agreement for one material or one band range does not establish convergence for another system.

## Notes

- All output files are overwritten if they already exist.
- `WAVEDER` and `EIGENVAL` must originate from the same VASP calculation.
- WIEN2k matrix-element, energy, and occupation files must represent the same calculation and compatible band ranges.
- The built-in `--help` text should be kept synchronized with this README whenever the command-line interface changes.
