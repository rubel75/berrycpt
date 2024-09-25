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

`/path/to/berrycpt/berrycpt case.mommat2[up/dn] -nvb XX # if you use this GitHub version and WIEN2k (see the compatibility note above)`

`/path/to/berrycpt/berrycpt WAVEDER -nvb XX # VASP`

`/path/to/berrycpt/berrycpt WAVEDER -efermiev 3.123 # VASP`

Options:

  * `[-up/-dn]` tells `berrycpt` to read `case.mommat2[up/dn]` files (needed for spin-polarized calculations or calculations with spin-orbit coupling).

  * `-nvb XX` sets the number of occupied bands.

  * `-efermiev 3.123` sets the Fermi energy in (eV) to determine which bands are occupied (available in VASP only).


### Output

`bcurv_ij.dat` - contains elements of the Berry curvature tensor for each k point and band.

The file headers explain the content.