# BerryCPT
Berry curvature calculation with DFT using perturbation theory. Currently supported codes:
* [WIEN2k](http://www.wien2k.at)
* [VASP](https://www.vasp.at)

It is written in Fortran and intended for Linux OS

### Current Version

It is a work in progress (incomplete). Please wait until a fully functional version is released.


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

`x berrycpt [-up/-dn] [-settol 1.0e-5] # if you use the version built into WIEN2k starting with v20.1`

`/path/to/berrycpt/berrycpt case.mommat2[up/dn] [1e-5] # if you use this GitHub version and WIEN2k (see the compatibility note above)`

`/path/to/berrycpt/berrycpt WAVEDER [1e-5] # VASP`

Options:

  * `[-up/-dn]` tells `berrycpt` to read `case.mommat2[up/dn]` files (needed for spin-polarized calculations or calculations with spin-orbit coupling).

  * `[-settol 1.0e-5]` or `[1e-5]` is (optional) degeneracy energy tolerance [Ha], which is max dE for 2 states to be considered as degenerate (default value is 1.0e-6 Ha). Eigenvalues of states that are physically expected to be degenerate (for example light and heavy holes in GaAs at Gamma point) are not exactly the same for numerical reasons. If we treat them as non-degenerate states, a large error will occur due to the term p^2/dE when dE -> 0. These states need to be "bundled" together and treated as degenerate. It is recommended to inspect eigenvalues (in the case of GaAs it will be the top valence bands) and evaluate the energy difference dE between states that should be treated as degenerate. Then set the tolerance slightly greater than this value.


### Output

`bcurv_ij.dat` - contains elements of the Berry curvature tensor for each k point and band.

The file headers explain the content.