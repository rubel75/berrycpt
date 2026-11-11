# BerryCPT

BerryCPT calculates band-resolved Berry curvature, orbital angular momentum (OAM), and generalized orbital-angular-momentum matrices from density-functional-theory matrix elements.

Supported electronic-structure codes:

- [WIEN2k](https://www.wien2k.at)
- [VASP](https://www.vasp.at)

BerryCPT is written in Fortran and is intended for Linux-based systems.

## Status

The `main` branch is the reference version of BerryCPT. Every output file records the Git commit of the build in its first header line, and the same commit is printed at the end of each run. Results should be reported together with that commit. Older tagged releases predate the current implementation and should not be used.

The following functionality has been validated, as described in the accompanying paper (in preparation):

- Berry curvature and OAM from WIEN2k and VASP (MoS2 monolayer, compared between the two codes and with Wannier interpolation)
- sigma_z spin Berry curvature from WIEN2k with SOC (Bi2Se3 in the R-3m and Pnma phases, Sb2TeSe2)
- OAM compared with circular-dichroism ARPES (WSe2 monolayer, CoSi, CoSn)
- spin-projected OAM from WIEN2k with SOC (CoSn)
- both occupation modes, explicit Fermi energy and DFT occupations

Generalized OAM matrices are implemented but are not part of that validation. Spin-polarized calculations with SOC have not been validated.

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

The default `Makefile` configuration is the production build with Intel `ifort` (`-O2`). A debugging build and GNU Fortran settings are included as commented lines in the `Makefile`.

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
- [Generating VASP `WAVEDER` files](https://github.com/rubel75/mstar/wiki/Generate-WAVEDER-file-in-VASP)

The matrix-element and eigenvalue files used in one BerryCPT calculation
must originate from the same DFT calculation.

The number of unoccupied bands must be convergence-tested. Contributions
from remote bands decrease as $1/(\Delta E)^2$ for Berry curvature, but
only as $1/\Delta E$ for OAM and generalized OAM. Consequently, OAM
normally requires a larger energy window than Berry curvature.

### Matrix elements and limitations

All quantities are built from velocity matrix elements
$v_{\alpha,nn'} = \langle u_{n\mathbf{k}} | \hat{v}_\alpha | u_{n'\mathbf{k}} \rangle$
between cell-periodic Bloch states at the same $\mathbf{k}$.

WIEN2k (`OPTIC`) writes momentum matrix elements
$\mathbf{p}_{nn'} = \langle u_{n\mathbf{k}} | -\mathrm{i}\hbar\nabla | u_{n'\mathbf{k}} \rangle$
and, with SOC, their spin-diagonal parts $\mathbf{p}^{\uparrow\uparrow}_{nn'}$ and
$\mathbf{p}^{\downarrow\downarrow}_{nn'}$, with
$\mathbf{p}_{nn'} = \mathbf{p}^{\uparrow\uparrow}_{nn'} + \mathbf{p}^{\downarrow\downarrow}_{nn'}$
(WIEN2k 24.1, `SRC_optic/sph-UP.F` and `outmat.f`). BerryCPT converts them as
$\mathbf{v}_{nn'} = \mathbf{p}_{nn'}/m_0$. This has the following consequences:

- The conversion holds only for a local potential. A nonlocal potential
  $\hat{V}_\mathrm{NL}$ adds $(\mathrm{i}/\hbar)[\hat{V}_\mathrm{NL},\hat{\mathbf{r}}]$
  to the velocity operator, which WIEN2k does not provide. WIEN2k calculations
  with BerryCPT are therefore restricted to semilocal exchange-correlation
  functionals and exclude hybrid functionals.
- With SOC, the matrix elements are taken between spinor states, so SOC enters
  through the wave functions. The spin-orbit correction to the velocity
  operator itself, $(\hbar/4m_0^2c_0^2)\,\hat{\boldsymbol{\sigma}}\times\nabla V$,
  is not included.

VASP (`WAVEDER`) provides the derivative overlaps in the form
$\langle u_{n\mathbf{k}} | -\mathrm{i}\nabla_{\mathbf{k}} u_{n'\mathbf{k}} \rangle = -\mathbf{r}_{nn'}$,
where $\mathbf{r}_{nn'} = \mathrm{i}\langle u_{n\mathbf{k}} | \nabla_{\mathbf{k}} u_{n'\mathbf{k}} \rangle$
are the dipole matrix elements. They include the nonlocal projector augmented-wave
terms, so the resulting matrix elements are velocity matrix elements. This
convention, including the factor $-\mathrm{i}$, is the one defined in the VASP 6.6.1
source (`linear_optics.F`). The VASP wiki page on `WAVEDER` (revision of
15 April 2026, read 11 November 2026) gives the overlap without the factor
$-\mathrm{i}$. BerryCPT follows the source code. BerryCPT uses

```math
\mathbf{v}_{nn'} = \frac{\mathrm{i}}{\hbar}\left(E_{n} - E_{n'}\right)\mathbf{r}_{nn'} ,
\qquad n \neq n' .
```

VASP does not write spin-diagonal matrix elements, so spin-resolved quantities
are available only from WIEN2k.

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

A spin-polarized SOC workflow has not been validated. Use matching spin-labelled WIEN2k files and verify the results independently before relying on that configuration.

## Tests

The folder `test_files` contains one test for each DFT code, with the input files, the DFT input and output in `etc`, and the reference results in `expected_output`. Run each test from inside its folder with the command given in its `README`:

```bash
cd test_files/WIEN2k
../../berrycpt CoSi-berrycpt.mommat2 CoSi-berrycpt.mommat2up CoSi-berrycpt.mommat2dn \
    --enefile CoSi-berrycpt.energyso --efermiry 0.8782005017

cd ../VASP
../../berrycpt WAVEDER --enefile EIGENVAL --efermiev 11.97
```

Both tests use CoSi. The WIEN2k test (SOC with spin-resolved matrix elements) produces every output file type. The VASP test produces the ordinary quantities. Apart from the first header line, which records the Git commit, the results should agree with `expected_output` up to rounding in the last digits, which depends on the compiler and the LAPACK library.

## Degenerate bands

Bands are assigned to contiguous degenerate groups at every k-point and spin channel. Adjacent bands are grouped when their energy separation does not exceed $10^{-5}$ Ha. This tolerance is fixed in the code (`energy_tolerance` in `berrycpt.f90`).

For an isolated band, BerryCPT evaluates the usual band-resolved perturbation-theory expression (see [Output conventions](#output-conventions)). For a degenerate group $D$, it constructs a Hermitian effective matrix within the group from the couplings to all states $n' \notin D$ and reports its eigenvalues. For the Berry curvature,

```math
M_{nm}^{\Omega} = \mathrm{i}\hbar^2 \sum_{n' \notin D}
\frac{v_{\alpha,nn'}\, v_{\beta,n'm} - v_{\beta,nn'}\, v_{\alpha,n'm}}
{(E_{n'} - E_{n})(E_{n'} - E_{m})} ,
\qquad n, m \in D ,
```

and for the OAM,

```math
M_{nm}^{L} = \frac{\mathrm{i} m_0 \hbar}{2} \sum_{n' \notin D}
\left( v_{\alpha,nn'}\, v_{\beta,n'm} - v_{\beta,nn'}\, v_{\alpha,n'm} \right)
\left[ \frac{1}{E_{n'} - E_{n}} + \frac{1}{E_{n'} - E_{m}} \right] ,
\qquad n, m \in D .
```

For a single band ($n = m$) both reduce to the non-degenerate expressions. The reported values inside a degenerate block are eigenvalues of the effective matrix, not expectation values attached to the original unrotated DFT eigenvectors.

The target-band cutoff is adjusted so that it never divides a degenerate group. For WIEN2k, the initial target range is up to twice the largest occupied-band index, limited by the number of available bands. For VASP, the initial target range is read from `NBANDS_CDER` in `WAVEDER`.

## Output conventions

The pseudovector component order is:

```text
yz, zx, xy
```

These correspond to Voigt indices 4, 5, and 6, so that, for example, $\Omega_{xy,n}$ is the $z$ component of the pseudovector $\boldsymbol{\Omega}_n$.

Units and conventions:

- Berry curvature and spin Berry curvature: bohr$^2$
- OAM, spin-projected OAM, and generalized OAM: $\hbar$
- Wave vectors include the factor $2\pi$. The Berry curvature therefore has the dimension of length squared, and its integral over a closed two-dimensional Brillouin zone equals $2\pi$ times the Chern number, $\int_\mathrm{BZ} \Omega_{xy}(\mathbf{k})\, d^2k = 2\pi C$.

The units are also stated in the header of every output file. Existing output files are replaced. Every output file contains a detailed header describing its columns, k-point organization, band ranges, degenerate-block labels, and units. At normal termination, BerryCPT prints a summary of the output files that were created and closed successfully.

The expressions below are written in the same form as in the accompanying paper. The $\mathbf{k}$ label is omitted, and all matrix elements and energies refer to the same $\mathbf{k}$ point. Internally BerryCPT works with velocity matrix elements in Hartree atomic units ($\hbar = m_0 = 1$), where a WIEN2k momentum matrix element $\mathbf{p}_{nn'}$ and the velocity $\mathbf{p}_{nn'}/m_0$ have the same numerical value.

### Ordinary Berry curvature

Output file:

```text
bcurv_ij[-up/-dn].dat
```

The optional suffix identifies an ordinary spin channel, for example a collinear VASP spin channel or a separately processed WIEN2k `up` or `dn` input.

For an isolated band $n$,

```math
\Omega_{\alpha\beta,n} = -2\hbar^2 \,\mathrm{Im} \sum_{n' \neq n}
\frac{v_{\alpha,nn'}\, v_{\beta,n'n}}{(E_{n} - E_{n'})^2} .
```

Each k-point block contains:

```text
# KP: k  NBANDS_OUT: Nout  NBANDS_TRANS: Ntrans
band  block  Omega_yz  Omega_zx  Omega_xy
```

The final record for each k-point has `band = 0` and `block = 0`. It contains the occupation-weighted total

```math
\Omega_{\alpha\beta} = \sum_{n} f_{n}\, \Omega_{\alpha\beta,n} ,
```

where $f_n$ are the occupations (step-function or from the DFT calculation).

### Spin-resolved Berry-curvature quantities

These files are written only when three WIEN2k matrix-element files are supplied:

```text
bcurv_ij-up.dat
bcurv_ij-dn.dat
bcurv_ij-up-dn.dat
```

The file `bcurv_ij-up-dn.dat` contains the $\sigma_z$ spin Berry curvature

```math
\Omega^{(\sigma_z)}_{\alpha\beta,n} = -2\hbar^2 \,\mathrm{Im} \sum_{n' \neq n}
\frac{J^{(\sigma_z)}_{\alpha,nn'}\, v_{\beta,n'n}}{(E_{n} - E_{n'})^2} ,
\qquad
J^{(\sigma_z)}_{\alpha,nn'} = v^{\uparrow\uparrow}_{\alpha,nn'} - v^{\downarrow\downarrow}_{\alpha,nn'} ,
```

where $J^{(\sigma_z)}_{\alpha,nn'}$ is the matrix element of the $\sigma_z$-normalized spin current operator $\hat{J}^{(\sigma_z)}_\alpha = \frac{1}{2}\{\hat{v}_\alpha, \hat{\sigma}_z\}$. The physical spin current, associated with $\hat{s}_z = (\hbar/2)\hat{\sigma}_z$, is larger by the factor $\hbar/2$. The spin Berry curvature depends on the choice of the spin quantization axis set in the WIEN2k SOC calculation.

The files `bcurv_ij-up.dat` and `bcurv_ij-dn.dat` contain the spin-up and spin-down contributions, obtained by replacing $J^{(\sigma_z)}_{\alpha,nn'}$ with $v^{\uparrow\uparrow}_{\alpha,nn'}$ or $v^{\downarrow\downarrow}_{\alpha,nn'}$. Their difference is the $\sigma_z$ spin Berry curvature.

For non-degenerate bands, linear combinations of these quantities follow directly from the operator definitions. Inside a degenerate block, the three effective matrices are diagonalized separately. Their ordered eigenvalues must not be combined band by band. Linear relations remain valid at the matrix level and for block traces.

### Ordinary OAM

Output file:

```text
oam_ij[-up/-dn].dat
```

For an isolated band $n$,

```math
L_{\alpha\beta,n} = 2 m_0 \hbar \,\mathrm{Im} \sum_{n' \neq n}
\frac{v_{\alpha,nn'}\, v_{\beta,n'n}}{E_{n} - E_{n'}} .
```

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

This file is written only when the three spin-resolved WIEN2k matrix-element files are supplied. It contains the OAM projected onto the spin-up subspace,

```math
\mathbf{L}^{(\sigma_z\uparrow)}_{n} = \langle u_{n} | \hat{P}_{\sigma_z\uparrow}\, \hat{\mathbf{L}}\, \hat{P}_{\sigma_z\uparrow} | u_{n} \rangle ,
\qquad
\hat{P}_{\sigma_z\uparrow} = \tfrac{1}{2}\left(I_2 + \hat{\sigma}_z\right) ,
```

evaluated from the spin-diagonal velocity matrix elements as

```math
L^{(\sigma_z\uparrow)}_{\alpha\beta,n} \approx 2 m_0 \hbar \,\mathrm{Im} \sum_{n' \neq n}
\frac{v^{\uparrow\uparrow}_{\alpha,nn'}\, v^{\uparrow\uparrow}_{\beta,n'n}}{E_{n} - E_{n'}} .
```

This is an approximation. For spinors that are not spin-pure, the spin-up and spin-down projections evaluated in this way do not add up to the ordinary OAM, $L^{(\sigma_z\uparrow)} + L^{(\sigma_z\downarrow)} \neq L$. The missing terms are the products $v^{\uparrow\uparrow}_{\alpha,nn'} v^{\downarrow\downarrow}_{\beta,n'n}$ and $v^{\downarrow\downarrow}_{\alpha,nn'} v^{\uparrow\uparrow}_{\beta,n'n}$, which enter $L$ but neither spin-projected sum.

The output has the same band, block, and component structure as ordinary OAM.

### Generalized OAM

Output file:

```text
goam_ij_nm[-up/-dn].dat
```

For each component, BerryCPT writes the complete complex Hermitian matrix $L_{\alpha\beta,nm} = \langle u_{n} | \hat{L}_{\alpha\beta} | u_{m} \rangle$. The first matrix index is the bra-band index and the second is the ket-band index. The matrix remains in the original DFT eigenstate basis and is not diagonalized within degenerate blocks. Its off-diagonal elements therefore depend on the phases of the DFT eigenvectors, which each DFT code fixes arbitrarily at every $\mathbf{k}$ point. Only gauge-invariant combinations are physical: the diagonal elements of non-degenerate bands, $|L_{\alpha\beta,nm}|^2$ between non-degenerate bands, and the eigenvalues or traces of the sub-blocks belonging to degenerate subspaces.

BerryCPT evaluates

```math
L_{\alpha\beta,nm} = \frac{\mathrm{i} m_0 \hbar}{2} \sum_{n'}
\left( v_{\alpha,nn'}\, v_{\beta,n'm} - v_{\beta,nn'}\, v_{\alpha,n'm} \right)
\left[ \frac{Q_{n}(n')}{E_{n'} - E_{n}} + \frac{Q_{m}(n')}{E_{n'} - E_{m}} \right] .
```

$Q_n(n')$ is zero when the intermediate state $n'$ belongs to the same degenerate block as the external state $n$, and one otherwise. $Q_m(n')$ is defined analogously. These factors exclude singular same-block intermediate-state couplings. The diagonal elements of a non-degenerate band reduce to the ordinary OAM.

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

## How to cite

Please cite BerryCPT in publications and other scholarly work that use this software or results generated with it, giving the Git commit from which the results originated:

```text
Rubel, O. (2024). BerryCPT: Calculates Berry curvature using perturbation
theory [Computer software]. https://github.com/rubel75/berrycpt,
commit <commit hash>.
```

The commit is printed at the end of every run and in the first header line of every output file, in the form `Git revision: <date>-<commit hash>`.
