# BerryCPT
BerryCPT is a code for calculating Berry curvature and orbital angular momentum (OAM) tensors using perturbation theory and DFT outputs. Supported codes:
* [WIEN2k](http://www.wien2k.at)
* [VASP](https://www.vasp.at)

It is written in Fortran and intended for Linux-based systems.

### Current status

**Near-final version:** All planned features have been implemented. The code is currently being tested and reviewed before release.


### Installation:
First, clone the GitHub repository

`$ git clone https://github.com/rubel75/berrycpt`

**WIEN2k compatibility note**:
The format of the `case.mommat2` file has slightly changed in v20.1 to enable calculations for the number of bands greater than 9999. This `berrycpt` code is compatible with these changes (a new non-fixed format). If you would like to use this code in conjunction with WIEN2k _prior_ to v20.1, you need to comment/uncomment two lines in the `berrycpt/read_mommat_pij.f90` file making the following changes before compiling the code:

```
    READ(cline,'(3X,2I4,6E13.6,F13.8)',ERR=20) bii ,bjj, & !...
        p1_Re, p1_Im, p2_Re, p2_Im, p3_Re, p3_Im, dEij(bii,bjj)
    
    ! WIEN2k after May 2020 (the case.mommat2 file is not a fixed format)
    !READ(cline,*,ERR=30) bii ,bjj, & !...
    !    p1_Re, p1_Im, p2_Re, p2_Im, p3_Re, p3_Im, dEij(bii,bjj)
```

The `makefile` is set up for the Intel Fortran compiler `ifort` and Intel MKL. (To compile with `gfortran`, you need to uncomment corresponding `FC`, `FCFLAGS`, and `FLFLAGS` variables in the `makefile`). To compile, simply execute

`$ cd berrycpt; make`


### Execution
First, you need to perform a standard SCF (self-consistent field) calculation and generate a file that contains optical matrix elements (`case.mommat2[up/dn]` in WIEN2k or `WAVEDER` in VASP). Tips on how to do this can be found at this [Wiki page](https://github.com/rubel75/mstar/wiki). Once the file is ready, execute

`/path/to/berrycpt/berrycpt case.mommat2[up/dn] -nvb XX # if you use this GitHub version and WIEN2k (see the compatibility note above)`

`/path/to/berrycpt/berrycpt WAVEDER -nvb XX # VASP`

`/path/to/berrycpt/berrycpt WAVEDER -efermiev 3.123 # VASP`

Options:

  * `[-up/-dn]` tells `berrycpt` to read `case.mommat2[up/dn]` files (needed for spin-polarized calculations or calculations with spin-orbit coupling).

  * `-nvb XX` sets the number of bands used in the calculation. For total Berry curvature, this corresponds to the number of occupied bands. For band-resolved Berry curvature and OAM, it defines the range of bands to include (not necessarily all occupied).

  * `-efermiev 3.123` sets Fermi energy (in eV) to determine occupied bands (VASP only).


### Output

- `bcurv_ij.dat`  
  Total Berry curvature tensor computed from the full spinor wavefunctions.

- `bcurv_ij-up.dat`, `bcurv_ij-dn.dat`  
  Spin-projected Berry curvature components. These are computed by replacing one of the momentum matrix elements of the standard expression:

      Im[ ⟨uₙ | v_a | uₘ⟩ × ⟨uₘ | v_b | uₙ⟩ ]

  with a spin-projected form:

      ⟨uₘ^↑ | v_a | uₙ^↑⟩   or   ⟨uₘ^↓ | v_a | uₙ^↓⟩

  This is **not** equivalent to evaluating the Berry curvature with the **spin current operator**, which is defined as:

      J_a^σ_i = (1/2) { v_a, σ_i }

  Hence, `bcurv_ij-up.dat` and `bcurv_ij-dn.dat` reflect diagonal components of the spin current operator J_a^σ_z.


- `bcurv_ij-up-dn.dat`  
  Spin-z Berry curvature tensor components. This file contains the difference between the spin-up and spin-down projected Berry curvature:

      Ω_ab^{s_z} = Ω_ab^↑ − Ω_ab^↓

  which is equivalent to computing the Berry curvature using the **spin current operator**:

      J_b^{s_z} = (1/2) { v_a, σ_z }

- `oam_ij.dat`  
  Orbital angular momentum (OAM) tensor for each k-point and band.

- `oam_ij-sigma_z-up.dat`  
  Spin-up projected orbital angular momentum, computed using second-order perturbation theory based on velocity (or momentum) matrix elements between spinor Bloch states:

      L_{ab,n}^{σ_z,↑}(k) = sum_{m ≠ n} Im[ v_a(n,m)^↑↑ * v_b(m,n)^↑↑ ] / (ε_m(k) - ε_n(k))
  
  where:
  - `v_a(n,m)^↑↑` = ⟨ u_{n,k,↑} | v_a | u_{m,k,↑} ⟩ is the velocity matrix element between spin-up components along direction `a`,
  - `ε_n(k)` and `ε_m(k)` are the band energies at the k-point,
  - `a` and `b` are Cartesian directions (`x`, `y`, or `z`),
  - the sum runs over all bands `m ≠ n`.
 
  This expression corresponds to the spin-up-z projection of the orbital angular momentum tensor in the form:

      L_{ab}^{σ_z,↑}(k) = ⟨ u_{n,k} | L_{ab} * P_{↑_z} | u_{n,k} ⟩

  with `P_{↑_z} = 1/2 * (I + σ_z) = [[1, 0], [0, 0]]` being the projection operator onto spin-up states. Here `I = [[1, 0], [0, 1]]` is the identity matrix.

> **Note**: All files will be overwritten if they exist from a previous run.
