PROGRAM berrycpt
! Includes subprograms:
!     degenbc.f90       -- Handles perturbation theory for degenerate states when 
!                          calculating Berry curvature
!     degenoam.f90      -- Handles perturbation theory for degenerate states when 
!                          calculating orbital angular momentum (OAM)
!     eigvz.f90         -- Solves a complex eigenvalue problem for a Hermitian matrix
!     read_mommat_nb.f90 -- Determines the number of bands in a k-point block of 
!                          the mommat file
!     read_mommat_pij.f90 -- Reads momentum matrix elements <i|p_a|j> (a = x, y, z) 
!                          and energy differences E_i - E_j in a k-point block of the mommat file
!     read_numlines.f90 -- Reads the number of lines in a file
!
! (c) Oleg Rubel, Sep 2024
!
! Execution:
!   $ ./berrycpt arg1 -nvb arg3 [-so]
!     arg1  - input mommat file name from WIEN2k or VASP WAVEDER file
!     arg3  - number of occupied bands
!     -so   - optional switch for WIEN2k with spin-orbit coupling (SOC) 
!             (including spin-resolved Berry curvature)
!
!   $ ./berrycpt WAVEDER -efermiev arg3  # VASP only!
!     arg3  - Fermi energy (in eV) to determine occupied bands
!
! Example (WIEN2k, spin-polarized, SOC):
!   ./berrycpt case.mommat2 -nvb 56 -so
!
! Example (VASP with WAVEDER, using Fermi energy):
!   ./berrycpt WAVEDER -efermiev 6.52
!
! Output:
!     bcurv_ij.dat         - contains elements of the Berry curvature tensor
!     bcurv_ij-up.dat      - spin-up Berry curvature tensor elements
!     bcurv_ij-dn.dat      - spin-down Berry curvature tensor elements
!     bcurv_ij-up-dn.dat   - spin_z Berry curvature elements
!     oam_ij.dat           - contains elements of the orbital angular momentum tensor
!     (Note: Existing files will be overwritten)
!
! Tips:
! (1) Writing of the mommat file is _not_ enabled by default in WIEN2k.
!     To enable writing, edit the case.inop file and change
!     OFF to ON in the following line:
!         ON           ON/OFF   WRITEs MME to unit 4
!         -^
!
! (2) Ensure _plenty_ of empty bands are included during SCF:
!     (a) Increase "de" in case.in1(c) above 5 Ry:
!         K-VECTORS FROM UNIT:4   -9.0      10.0    10   emin / de (emax = Ef + de) / nband
!                                     ------------------^
!
!     (b) For SOC calculations, extend "Emax" in case.inso up to 5 Ry:
!         -10 5.0                Emin, Emax
!         -----^
!
!     (c) The default "Emax" in case.inop is 3 Ry (usually OK),
!         but consider testing convergence and increasing up to 5 Ry:
!         -5.0 3.5 9999         Emin, Emax for matrix elements, NBvalMAX
!         ------^
!
! (3) For VASP calculations:
!     - Set LOPTICS = .TRUE.
!     - Increase NBANDS (at least x3)
!     - Disable finite-difference derivatives: LPEAD = .FALSE.

!! Variables

USE OMP_LIB
USE precision_mod, ONLY: sp, dp
USE read_mommat_nb_mod, ONLY: read_mommat_nb
USE read_numlines_mod, ONLY: read_numlines
USE read_numlines_vasp_mod, ONLY: read_numlines_vasp
USE read_mommat_pij_mod, ONLY: read_mommat_pij
USE read_mommat_pij_vasp_mod, ONLY: read_mommat_pij_vasp
USE finddegenblocks_mod, ONLY: finddegenblocks
USE degenbc_mod, ONLY: degenbc
USE degenoam_mod, ONLY: degenoam

IMPLICIT NONE

!! Variables

CHARACTER(len=256) :: &
    arg1, arg2, arg3, arg4, & ! command line input arguments
    fnameinp, & ! input file with momentum or dipole matrix elements
    fnameinpUP, fnameinpDN, & ! spin-resolved mom. matr. elements (WIEN2k only)
    fnameout2,  & ! output file names
    fnameout21, fnameout22, fnameout23, & ! ... spin UP, DN, UP-DN Berry curvature
    fnameout3, fnameout31,  & ! output file names for OAM
    wformat2, & ! format for writing/reading data
    charspin ! spin component for spin-polarized calculations
INTEGER :: &
    nltot, & ! total number of lines in mommat file
    nltotUP, nltotDN, & ! total number of lines in mommat2[up/dn] files (WIEN2k)
    nktot, & ! total number of k-points in WAVEDER file
    nstot, & ! total number of spins in WAVEDER file
    iline, & ! current line number during reading of mommat file
    ilineUP, ilineDN, & ! ... for up/dn components (WIEN2k)
    ispin, & ! current spin (1/2)
    nb, & ! number of bands for a current k-point
    nbUP, nbDN, & ! ... for up/dn components (WIEN2k unly)
    nbb, & ! number of band-to-band transition 
    nbcder, & ! the max number of bands for which the Berry curvature
              ! is calculated
              ! In WIEN2k we set it = to the total number of bands
              ! In VASP it is set in WAVEDER file as NBANDS_CDER
    ivb, i, j, ikpt, m, & ! counters
    idg1, idg2, & ! indices of the first and last degenerate states in a block
    n, & ! band indices
    nvb, & ! number of occupied bands [1:nvb]
    nvbinpt=0, &    ! input number of occupied bands [1:nvb]
                    ! (init. 0 as a flag of the variable being set)
    ivoigt, & ! Voigt index (1..6)
    alpha, beta, & ! Cartesian directions 1,2,3 = x,y,z
    ierr, & ! error code
    ngroups, & ! number of groups of degenerate energies
    ig, & ! group index
    nmg ! number of bands as members in a degenerate group
INTEGER, ALLOCATABLE :: &
    nbocck(:,:), & ! k- and spin-resolved number of occupied bands
    dg_group(:), & ! assignment of each band to a group of degenerate bands
    members(:) ! list of bands that are members of a degenerate group
REAL(kind=sp) :: &
    efermi ! Fermi energy in [eV]
REAL(kind=sp), ALLOCATABLE :: &
    dEij(:,:), & ! energy differences E_i-E_j [Ha]
    dEijdg(:,:), & ! energy differences E_i-E_j [Ha] for degenerate bands
    dEijUP(:,:), dEijDN(:,:), & ! ... for up/dn components (WIEN2k)
    dEijks(:,:,:,:) ! k- and spin-dependent energy differences E_i-E_j [Ha]
REAL(kind=dp), ALLOCATABLE :: &
    bcurv(:,:), & ! array to store Berry curvature
    bcurvdg(:), & ! array to store Berry curvature of a block of degenerate bands
    oam(:,:), & ! array to store orbital angular momentum (OAM)
    oamdg(:) ! array to store OAM of a block of degenerate bands
COMPLEX(kind=sp), ALLOCATABLE :: &
    pij(:,:,:), & ! momentum matrix elements [at.u.]
    pijA(:,:), pijB(:,:), & ! subset of momentum matrix elements [at.u.] for degenerate bands
    pijUP(:,:,:), pijDN(:,:,:), & ! ... for up/dn components (WIEN2k)
    pijks(:,:,:,:,:) ! k- and spin-dependent momentum matrix elements [at.u.]
LOGICAL :: &
    fmommatend, & ! end of mommat file
    file_exists, &
    wien2k, & ! true if this is a WIEN2k calculation
    spinor ! the wave function is a spinor

!! Get command line input arguments

DO i = 1, iargc() ! loop through all command line arguments to search for help
    CALL GETARG(i,arg1)
    IF ( TRIM(arg1)=='-h' .or.  TRIM(arg1)=='--h' .or. & !...
            TRIM(arg1)=='-help' .or. TRIM(arg1)=='--help') THEN
        GOTO 911 ! print help and STOP
    END IF
END DO
WRITE(*,'(A,I0)') ' Detected input arguments = ', iargc()
IF ((iargc() < 3) .OR. (iargc() > 4)) THEN ! check number of input arguments
    WRITE(*,'(A)') ' Expected 3 or 4 input arguments'
    GOTO 912 ! print error and STOP
END IF
! 1st argument
CALL GETARG(1,arg1) ! mommat file name
WRITE (*,*) 'Input mommat file = ' // TRIM(arg1)
fnameinp = TRIM(arg1)
! 2nd argument
CALL GETARG(2,arg2) ! switch for the number of valence bands [-nvb]
IF (TRIM(arg2)=='-nvb') THEN
    ! 3rd argument can be the number of occupied bands
    CALL GETARG(3,arg3)
    read(arg3,*,IOSTAT=ierr) nvbinpt
    IF (ierr /= 0) THEN ! input error for 3rd argument
        WRITE(*,*) '3rd argument is "' // TRIM(arg3) // '"'
        WRITE(*,*) 'Error detected for the 3rd input argument (must be a number of occupied bands)'
        GOTO 912 ! print error and STOP
    END IF
    WRITE (*,'(A,I0)') & !...
        ' Number of occupied bands = ', nvbinpt
ELSE IF (TRIM(arg2)=='-efermiev') THEN
    ! 3rd argument can be the Fermi energy, which will be used to determine
    ! spin and k-specific number of occupied bands
    CALL GETARG(3,arg3)
    read(arg3,*,IOSTAT=ierr) efermi
    IF (ierr /= 0) THEN ! input error for 3rd argument
        WRITE(*,*) '3rd argument is "' // TRIM(arg3) // '"'
        WRITE(*,*) 'Error detected for the 3rd input argument (must be a Fermi energy)'
        GOTO 912 ! print error and STOP
    END IF
    WRITE (*,'(A,F10.6)') & !...
        ' Fermi energy (eV) = ', efermi
ELSE ! impossible
    WRITE (*,*) 'The second argument is "' // TRIM(arg2) // &
        '", while expected "-nvb" or "-efermiev"'
    GOTO 912 ! print error and STOP
END IF
! 4th argument
spinor = .false. ! default is _not_ a spinor wave function (matters for WIEN2k only)
IF ( iargc()==4 ) THEN
    CALL GETARG(4,arg4) ! switch for the spinor wave function [-so]
    IF ( TRIM(arg4)=='-so' ) THEN
        spinor = .true.
        WRITE(*,*) 'The wave function is a spinor'
    ELSE ! impossible
        WRITE (*,*) 'The 4th argument is "' // TRIM(arg4) // &
            '", while expected "-so"'
        GOTO 912 ! print error and STOP
    END IF
END IF

!! WIEN2k or VASP?

IF ( TRIM(fnameinp) == 'WAVEDER' .OR. TRIM(fnameinp) == 'WAVEDERF' ) THEN ! VASP
    wien2k = .false.
    WRITE (*,*) 'Assume VASP calculation'
    ! change to a binary file in case it is pointed at the formatted 
    ! file WAVEDER
    fnameinp = 'WAVEDER'
    WRITE(*,*) 'Assumed VASP calculation based on the input file name.'
ELSE ! WIEN2k (default)
    wien2k = .true.
    WRITE(*,*) 'Assumed WIEN2k calculation based on the input file name.'
    WRITE(*,*) '(If you would like to read VASP file, the input file should'
    WRITE(*,*) 'have the exact name WAVEDER.)'
    IF ( spinor ) THEN
        ! take last two letters of the input file name
        charspin = fnameinp( LEN(TRIM(fnameinp))-1 : LEN(TRIM(fnameinp)) )
        IF ( (charspin=='up') .OR. (charspin=='dn')) THEN
            WRITE(*,'(A)') ' The input file name ends with ' // TRIM(charspin)
            WRITE(*,*) 'and the spinor option is selected. This is inconsistent.'
            WRITE(*,*) 'Instead you should point to the momentum matrix elements'
            WRITE(*,*) 'file with spin up+dn, which is case.mommat2'
            WRITE(*,*) 'Check execution options with "berrycpt -h"'
            ERROR STOP
        END IF
        fnameinpUP = TRIM(fnameinp)//'up'
        fnameinpDN = TRIM(fnameinp)//'dn'
        WRITE(*,*) 'It will be assumed that WIEN2k spin-resolved momentum matrix'
        WRITE(*,'(A)') ' elements are stored in the files: ' // &
                TRIM(fnameinpUP) // ' and ' // TRIM(fnameinpDN)
    END IF
END IF

!! Determine number of lines in mommat file & read VASP WAVEDER header

! check if the file exists
INQUIRE(FILE=fnameinp, EXIST=file_exists)
IF ( file_exists ) THEN
    WRITE(*,*) 'The input file ' // TRIM(fnameinp) // ' was found.'
ELSE IF ( .not.(file_exists) ) THEN
    WRITE(*,*) 'The input file ' // TRIM(fnameinp) // ' does not exist. Exiting'
    ERROR STOP
END IF
IF ( wien2k .AND. spinor ) THEN ! check spin-resolved momentum matrix elements
    INQUIRE(FILE=fnameinpUP, EXIST=file_exists)
    IF ( file_exists ) THEN
        WRITE(*,*) 'The input file ' // TRIM(fnameinpUP) // ' was found.'
    ELSE IF ( .not.(file_exists) ) THEN
        WRITE(*,*) 'The input file ' // TRIM(fnameinpUP) // ' does not exist. Exiting'
        ERROR STOP
    END IF
    INQUIRE(FILE=fnameinpDN, EXIST=file_exists)
    IF ( file_exists ) THEN
        WRITE(*,*) 'The input file ' // TRIM(fnameinpDN) // ' was found.'
    ELSE IF ( .not.(file_exists) ) THEN
        WRITE(*,*) 'The input file ' // TRIM(fnameinpDN) // ' does not exist. Exiting'
        ERROR STOP
    END IF
END IF
IF ( file_exists .AND. (.not.(wien2k)) ) THEN ! check EIGENVAL for VASP
    INQUIRE(FILE='EIGENVAL', EXIST=file_exists)
    IF ( .not.(file_exists) ) THEN
        WRITE(*,*) 'The file EIGENVAL is also required, but it ', &
            'does not exist. Exiting'
        ERROR STOP
    ELSE
        WRITE(*,*) 'The file EIGENVAL is also required. ', &
            'It is found.'
    END IF
END IF

!! Get the length of input files

IF ( wien2k ) THEN
    CALL read_numlines(fnameinp, 1, & ! <- args in 
            nltot) ! -> args out
    WRITE (*,'(A,I0)') '  number of lines in mommat file = ', nltot
    IF (nltot < 10) THEN
        WRITE(*,*) 'The file ' // TRIM(fnameinp) // ' is too short ' // & !...
            'and most likely useless. Stopping'
        ERROR STOP
    END IF
    nstot = 1 
    IF ( spinor ) THEN
        ! get the length of spin-resolved momentum matr. el. files
        CALL read_numlines(fnameinpUP, 1, & ! <- args in 
                nltotUP) ! -> args out
        CALL read_numlines(fnameinpDN, 1, & ! <- args in 
                nltotDN) ! -> args out
        IF ( (nltot.NE.nltotUP) .OR. (nltotUP.NE.nltotDN) ) THEN
            WRITE(*,*) 'Error: inconsistency in the length of input files' 
            WRITE(*,'(3A,I0)') '  the length of ', &
                    TRIM(fnameinp), ' = ', nltot
            WRITE(*,'(3A,I0)') '  the length of ', &
                    TRIM(fnameinpUP), ' = ', nltotUP
            WRITE(*,'(3A,I0)') '  the length of ', &
                    TRIM(fnameinpDN), ' = ', nltotDN
            WRITE(*,*) 'Stopping' 
            ERROR STOP
        END IF
    END IF
ELSE ! VASP
    CALL read_numlines_vasp(fnameinp, 1, & ! <- args in 
        nstot, nktot, nbcder, nb) ! -> args out
    WRITE (*,'(A,I0)') '  number of spins in WAVEDER file = ', nstot
    WRITE (*,'(A,I0)') '  number of k-points in WAVEDER file = ', nktot
    WRITE (*,'(A,I0)') '  smaller number of bands in WAVEDER file = ', nbcder
    WRITE (*,'(A,I0)') '  number of bands per k-point in EIGENVAL file = ', &
        nb
    ! Memory estimate: 
    ! 4 bytes * 2 (complex) * 3 (3D x,y,z) * 10e9 (Bytes -> GB)
    WRITE (*,'(A,1X,F5.1,1X,A)') ' Memory required to store matrix elements'//&
        ' from WAVEDER file', nstot*nktot*(nb**2.)*4*2*3/10.0**9., 'GB'
    WRITE (*,'(A,1X,F5.1,1X,A)') ' Memory required to store dE_ij'//&
        ' from EIGENVAL file', nstot*nktot*(nb**2.)*4*1/10.0**9., 'GB'
    WRITE (*,'(A,1X,F5.1,1X,A)') ' It will take additional', &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB overhead'//&
        ' to run the calculation'
    WRITE (*,'(A,1X,F5.1,1X,A)') ' Overall, you will need at least', &
        nstot*nktot*(nb**2.)*4*2*3/10.0**9. + &
        nstot*nktot*(nb**2.)*4*1/10.0**9. + &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB (+ 20% incidental) of RAM'//&
        ' to run the calculation'
END IF
! open input mommat/WAVEDER file(s) for reading (file ID = 1, [11, 12])
IF ( wien2k ) THEN
    OPEN (1, file = TRIM(fnameinp), status = 'old', iostat=ierr) ! case.mommat2
    IF (ierr /= 0) THEN
        WRITE(*,*) 'Error opening output file ', TRIM(fnameinp)
        ERROR STOP
    END IF
    IF ( spinor ) THEN
        OPEN (11, file = TRIM(fnameinpUP), status = 'old', &
            iostat=ierr) ! case.mommat2up
        IF (ierr /= 0) THEN
            WRITE(*,*) 'Error opening output file ', TRIM(fnameinpUP)
            ERROR STOP
        END IF
        OPEN (12, file = TRIM(fnameinpDN), status = 'old', &
            iostat=ierr) ! case.mommat2dn
        IF (ierr /= 0) THEN
            WRITE(*,*) 'Error opening output file ', TRIM(fnameinpDN)
            ERROR STOP
        END IF
    END IF
ELSE ! VASP
    OPEN (1, file = TRIM(fnameinp), form = 'unformatted', &
        status = 'old', iostat=ierr)
    IF (ierr /= 0) THEN
        WRITE(*,*) 'Error opening output file ', TRIM(fnameinp)
        ERROR STOP
    END IF
END IF

!! Read VASP WAVEDER and EIGENVAL files to determine matrix elements

IF ( .NOT.(wien2k) ) THEN
    CALL read_mommat_pij_vasp (1, nstot, nktot, nbcder, nb, efermi, & ! <- args in
            pijks, dEijks, nbocck) ! -> args out (allocated inside)
END IF

! Loop over spins (for VASP only)
! In WIEN2k nstot=1 since case.momat2up and case.momat2dn files should be
! read one after the other any ways
WRITE(*,*) 'Entering the main loop...'
DO ispin = 1, nstot

    !! Spin suffix
    
    IF ( wien2k ) THEN
        ! take last two letters of the input file name
        charspin = fnameinp( LEN(TRIM(fnameinp))-1 : LEN(TRIM(fnameinp)) )
        IF ( TRIM(charspin) == 'up' ) THEN
            charspin = '-up'
        ELSE IF ( TRIM(charspin) == 'dn' ) THEN
            charspin = '-dn'
        ELSE
            charspin = '' ! no spin identity
        END IF
    ELSE ! VASP
        IF (ispin==1 .AND. nstot==2) THEN
            charspin = '-up'
        ELSE IF (ispin==2 .AND. nstot==2) THEN
            charspin = '-dn'
        ELSE
            charspin = '' ! no spin identity
        END IF
    END IF

    !! Prepare output files

    fnameout2 = 'bcurv_ij'//TRIM(charspin)//'.dat'
    OPEN (2, file = TRIM(fnameout2), status = 'UNKNOWN') ! output file 1
    WRITE(2,'(A)') '# This file is generated by berrycpt'
    WRITE(2,'(A)') '# the output contains components of the Berry curvature tensor Omega_ij'
    WRITE(2,'(A)') '# that are grouped by k-point index and then by the band index.'
    WRITE(2,'(A)') '# The last entry for each k-point with the band index "0" is the'
    WRITE(2,'(A)') '# total Berry curvature of all occupied bands.'
    WRITE(2,'(A)') '# Columns correspond to values of Omega_ij'
    WRITE(2,'(A)') '# band  4=yz;      5=xz;      6=xy'
    IF ( wien2k .AND. spinor ) THEN
        fnameout21 = 'bcurv_ij-up.dat'
        OPEN (21, file = TRIM(fnameout21), status = 'UNKNOWN') ! output file 1
        WRITE(21,'(A)') '# This file is generated by berrycpt'
        WRITE(21,'(A)') '# the output contains components of the spin UP Berry curvature tensor Omega_ij'
        WRITE(21,'(A)') '# that are grouped by k-point index and then by the band index.'
        WRITE(21,'(A)') '# The last entry for each k-point with the band index "0" is the'
        WRITE(21,'(A)') '# total Berry curvature of all occupied bands.'
        WRITE(21,'(A)') '# Columns correspond to values of Omega_ij'
        WRITE(21,'(A)') '# band  4=yz;     5=xz;      6=xy'
        fnameout22 = 'bcurv_ij-dn.dat'
        OPEN (22, file = TRIM(fnameout22), status = 'UNKNOWN') ! output file 2
        WRITE(22,'(A)') '# This file is generated by berrycpt'
        WRITE(22,'(A)') '# the output contains components of the spin DN Berry curvature tensor Omega_ij'
        WRITE(22,'(A)') '# that are grouped by k-point index and then by the band index.'
        WRITE(22,'(A)') '# The last entry for each k-point with the band index "0" is the'
        WRITE(22,'(A)') '# total Berry curvature of all occupied bands.'
        WRITE(22,'(A)') '# Columns correspond to values of Omega_ij'
        WRITE(22,'(A)') '# band  4=yz;      5=xz;      6=xy'
        fnameout23 = 'bcurv_ij-up-dn.dat'
        OPEN (23, file = TRIM(fnameout23), status = 'UNKNOWN') ! output file 3
        WRITE(23,'(A)') '# This file is generated by berrycpt'
        WRITE(23,'(A)') '# the output contains components of the spin UP-DN Berry curvature tensor Omega_ij'
        WRITE(23,'(A)') '# that are grouped by k-point index and then by the band index.'
        WRITE(23,'(A)') '# The last entry for each k-point with the band index "0" is the'
        WRITE(23,'(A)') '# total Berry curvature of all occupied bands.'
        WRITE(23,'(A)') '# Columns correspond to values of Omega_ij'
        WRITE(23,'(A)') '# band  4=yz;      5=xz;      6=xy'
    END IF

    fnameout3 = 'oam_ij'//TRIM(charspin)//'.dat'
    OPEN (3, file = TRIM(fnameout3), status = 'UNKNOWN') ! output file 1
    WRITE(3,'(A)') '# This file is generated by berrycpt'
    WRITE(3,'(A)') '# the output contains components of the orbital angular momentum tensor L_ij'
    WRITE(3,'(A)') '# that are grouped by k-point index and then by the band index.'
    WRITE(3,'(A)') '# Columns correspond to values of L_ij'
    WRITE(3,'(A)') '# band  4=yz;      5=xz;      6=xy'
    IF ( wien2k .AND. spinor ) THEN
        fnameout31 = 'oam_ij-sigma_z-up.dat'
        OPEN (31, file = TRIM(fnameout31), status = 'UNKNOWN') ! output file 1
        WRITE(31,'(A)') '# This file is generated by berrycpt'
        WRITE(31,'(A)') '# the output contains components of the spin-projected'
        WRITE(31,'(A)') '# orbital angular momentum tensor expressed as'
        WRITE(31,'(A)') '# L^{sigma_z, up}_{ij,n} = <u|L_{ij}*1/2*(I+sigma_z)|u>'
        WRITE(31,'(A)') '# that are grouped by k-point index and then by the band index.'
        WRITE(31,'(A)') '# Columns correspond to values of L^{sigma_z, up}_{ij,n}'
        WRITE(31,'(A)') '# band  4=yz;     5=xz;      6=xy'
    END IF

    !! Main part nested loops
    
    fmommatend = .false. ! end of mommat file is not reached
    ikpt = 0 ! initialize the counter for k-points
    nbb = 0 ! initialize the number of band-to-band transitions
    iline = 0 ! initialize the number of lines to skip

    DO WHILE (.not.(fmommatend)) ! loop over k-points until the file ends
        ikpt = ikpt + 1 ! count number of k-points

        !! Determine number of bands in mommat file

        IF (wien2k) THEN ! do for WIEN2k only (nb was read above in case of VASP)
            CALL read_mommat_nb(1, & ! <- args in
                    iline, & ! <-> args in-out
                    nb) ! -> args out
            ! the max number of bands for which the Berry curvature is calculated
            ! here (WIEN2k) we set it = to the total number of bands
            nbcder = nb
            IF ( spinor ) THEN
                ! same for case.mommat2[up/dn]
                CALL read_mommat_nb(11, & ! <- args in
                        ilineUP, & ! <-> args in-out
                        nbUP) ! -> args out
                CALL read_mommat_nb(12, & ! <- args in
                        ilineDN, & ! <-> args in-out
                        nbDN) ! -> args out
                IF ( (nb.NE.nbUP) .OR. (nb.NE.nbDN) .OR. &
                        (nbUP.NE.nbDN) ) THEN
                    WRITE(*,'(A,I0)') ' K-point ', ikpt
                    WRITE(*,'(A,I0,A)') '  file ' // TRIM(fnameinp) // ' has ', &
                            nb, ' bands'
                    WRITE(*,'(A,I0,A)') '  file ' // TRIM(fnameinpUP) // ' has ',&
                            nbUP, ' bands'
                    WRITE(*,'(A,I0,A)') '  file ' // TRIM(fnameinpDN) // ' has ', &
                            nbDN, ' bands'
                    ERROR STOP 'Error: inconsistency in the number of bands for this k-point' 
                END IF
            END IF
        END IF

        !! Read <b_i|p_a|b_j> from mommat file, a=1,2,3 (x,y,z)

        IF ( wien2k ) THEN
            nbb = (nb+nb**2)/2 ! number of band-to-band transitions
            ! spin up+dn momentum matrix elements
            CALL read_mommat_pij (1, nb, nbb, & ! <- args in
                    iline, & ! <-> args in-out
                    pij, dEij) ! -> args out (allocated inside)
            IF ( spinor ) THEN
                ! spin UP momentum matrix elements
                CALL read_mommat_pij (11, nbUP, nbb, & ! <- args in
                        ilineUP, & ! <-> args in-out
                        pijUP, dEijUP) ! -> args out (allocated inside)
                ! spin DN momentum matrix elements
                CALL read_mommat_pij (12, nbDN, nbb, & ! <- args in
                        ilineDN, & ! <-> args in-out
                        pijDN, dEijDN) ! -> args out (allocated inside)
            END IF
        ELSE ! VASP
            IF (.NOT. ALLOCATED(pij)) THEN
                ALLOCATE(pij(3, nb, nb))
            END IF
            IF (.NOT. ALLOCATED(dEij)) THEN
                ALLOCATE(dEij(nb, nb))
            END IF
            pij = pijks(:,:,:,ikpt,ispin)
            dEij = dEijks(:,:,ikpt,ispin)
        END IF

        !! VASP k- and spin-specific number of occupied bands

        IF (wien2k) THEN
            nvb = nvbinpt ! use input value for the last occupied band
        ELSE ! VASP
            IF ( nvbinpt > 0 ) THEN
                nvb = nvbinpt
            ELSE
                nvb = nbocck(ikpt,ispin)
            END IF
        END IF

        !! Identify band degeneracies

        CALL finddegenblocks(nb, dEij, 1.0e-5, & ! <- args in 
            dg_group, ngroups) ! -> args out (dg_group is allocated inside)

        !! Handle situation when the degenerate block extends past 'nvb'

        DO n = nvb+1, nb
            IF (dg_group(nvb) .NE. dg_group(n)) THEN ! end of deceneracy group
                IF (n-1 > nvb) THEN
                    WRITE(*,'(A,I0,A,I0,A,I0,A)') 'WARNING: While processing ' &
                        //'k-point ', ikpt, &
                        ', the last occupied band was updated from ', &
                        nvbinpt, ' to ', n-1, ' due to degeneracy.'
                END IF
                nvb = n-1
                EXIT ! DO loop
            END IF
        END DO

        !! Write information about the current k-point

        WRITE(2,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
        WRITE(3,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
        IF ( wien2k .AND. spinor ) THEN
            WRITE(21,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                    '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
            WRITE(22,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                    '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
            WRITE(23,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                    '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
            WRITE(31,'(A,I0,1X,A,I0,1X,A,I0)') & !...
                    '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
        END IF
        
        !! Prepare output formats for Berry curvatures
        
        WRITE(wformat2,'(I0)') nbcder ! make a character of the length 'nbcder'
        ! format line to WRITE Berry curvatures
        WRITE(wformat2,'(I0)') LEN(TRIM(wformat2))
        wformat2 = '(I' // TRIM(wformat2) // ',1X,5(es10.3,1X),es10.3)'
            
        !! Loop through blocks of occupied bands

        ! array to store non-spin-resolved Berry curvatures and OAM
        ! Size 3 is because of only off-diogonal components of the Berry
        ! curvature and OAM tensor (4=yz; 5=xz; 6=xy) are not = 0
        ALLOCATE( bcurv(nvb,3) )
        ALLOCATE( oam(nvb,3) )
        bcurv = 0.0_dp
        oam = 0.0_dp
        ! Loop over groups of degenerate bands in the range of [1:nvb].
        ! Each band is treated as a degenerate group of size 1 if not part 
        ! of a larger degenerate block.
        DO ig = 1, dg_group(nvb)
            nmg = COUNT(dg_group == ig) ! number of group members
            ALLOCATE( members(nmg) )
            ! get bands that are members of the group
            m = 0
            DO n = 1, nb
                IF (dg_group(n) == ig) THEN
                    m = m + 1
                    members(m) = n
                ELSE IF (m > 0) THEN
                    ! break the loop once first condition is not true any more,
                    ! which marks the end of a degenerate block
                    EXIT
                END IF
            END DO
            ! loop over Voigt indices
            DO ivoigt = 4, 6
                ! handle Voigt notations
                SELECT CASE (ivoigt)
                    CASE (4); alpha = 2; beta = 3 ! 4=yz
                    CASE (5); alpha = 1; beta = 3 ! 5=xz
                    CASE (6); alpha = 1; beta = 2 ! 6=xy
                END SELECT
                ! degenerate group 1st and last band
                idg1 = members(1)
                idg2 = members(nmg)
                IF (idg1 .GT. idg2) THEN
                    WRITE(*,'(2(A,I0))') 'idg1=', idg1, ' idg2=', idg2
                    WRITE(*,*) 'Degenerate group misidentified: idg1 > idg2'
                    ERROR STOP
                END IF
                ! allocate group-specific arrays
                ALLOCATE( bcurvdg(nmg), pijA(nmg,nb), &
                    pijB(nb,nmg), dEijdg(nmg,nb) )
                ! get group-specific matrix elements and energies
                DO m = 1, nmg
                    pijA(m,:) = pij(alpha, members(m), :)
                    pijB(:,m) = pij(beta, :, members(m))
                    dEijdg(m,:) = dEij(members(m), :)
                END DO
                ! solve degenerate PT problem
                CALL degenbc(nb, idg1, idg2, & ! <- args in
                    pijA, pijB, dEijdg, & ! <- args in 
                    bcurvdg) ! -> args out
                CALL degenoam(nb, idg1, idg2, & ! <- args in
                    pijA, pijB, dEijdg, & ! <- args in 
                    oamdg) ! -> args out (oamdg allocated inside)
                ! store group result for the Berry curvature and OAM in array
                DO m = 1, nmg
                    bcurv(members(m), ivoigt-3) = bcurvdg(m)
                    oam(members(m), ivoigt-3) = oamdg(m)
                END DO
                DEALLOCATE(pijA, pijB, dEijdg, bcurvdg, oamdg)
            END DO ! loop over 'ivoigt'
            DEALLOCATE( members)
        END DO ! loop over 'ig'
        ! store the Berry curvature and OAM band-by-band for all occupied states
        DO ivb = 1, nvb
            WRITE(2,TRIM(wformat2)) ivb, (bcurv(ivb,j), j=1,3)
            WRITE(3,TRIM(wformat2)) ivb, (oam(ivb,j), j=1,3)
        END DO
        ! store the total Berry curvature for all occupied states
        WRITE(2,TRIM(wformat2)) 0, (SUM(bcurv(:,j)), j=1,3)
        DEALLOCATE( bcurv, oam )

        ! calculate spin-resolved and spin Berry curvature or OAM
        IF ( wien2k .AND. spinor ) THEN
            ! spin-resolved OAM
            ! OAM^{sigma_z, up}_{ab,n} = <u|L_{ab}*1/2*(I+sigma_z)|u>
            ! here I -- identity matrix,
            !      sigma_z -- Pauli matrix
            ! 1/2*(I+sigma_z) = (1 0; 0 0), leading to the PT sum over states
            ! OAM^{sigma_z, up}_{ab,n} = 
            !    ... sum_m 2*Im[p_{nm}^{up,up}*p_{mn}^{up,up}]/(E_n - E_m)
            ! in atomic units
            ALLOCATE( oam(nvb,3) )
            oam = 0.0_dp
            ! Loop over groups of degenerate bands in the range of [1:nvb].
            ! Even if a band is not degenerate, it is still considered as
            ! degeneracy of the size 1.
            DO ig = 1, dg_group(nvb)
                nmg = COUNT(dg_group == ig) ! number of group members
                ALLOCATE( members(nmg) )
                ! get bands that are members of the group
                m = 0
                DO n = 1, nb
                    IF (dg_group(n) == ig) THEN
                        m = m + 1
                        members(m) = n
                    ELSE IF (m > 0) THEN
                        ! break the loop once first condition is not true any more,
                        ! which marks the end of a degenerate block
                        EXIT
                    END IF
                END DO
                ! loop over Voigt indices
                DO ivoigt = 4, 6
                    ! handle Voigt notations
                    SELECT CASE (ivoigt)
                        CASE (4); alpha = 2; beta = 3 ! 4=yz
                        CASE (5); alpha = 1; beta = 3 ! 5=xz
                        CASE (6); alpha = 1; beta = 2 ! 6=xy
                    END SELECT
                    ! degenerate group 1st and last band
                    idg1 = members(1)
                    idg2 = members(nmg)
                    IF (idg1 .GT. idg2) THEN
                        WRITE(*,'(2(A,I0))') 'idg1=', idg1, ' idg2=', idg2
                        ERROR STOP 'Wrong boundaries of the degenerate block'
                    END IF
                    ! allocate group-specific arrays
                    ALLOCATE( oamdg(nmg), pijA(nmg,nb), pijB(nb,nmg), &
                        dEijdg(nmg,nb))
                    ! get group-specific matrix elements and energies
                    DO m = 1, nmg
                        pijA(m,:) = pijUP(alpha, members(m), :)
                        pijB(:,m) = pijUP(beta, :, members(m))
                        dEijdg(m,:) = dEij(members(m), :)
                    END DO
                    ! solve degenerate PT problem
                    CALL degenoam(nb, idg1, idg2, & ! <- args in
                        pijA, pijB, dEijdg, & ! <- args in 
                        oamdg) ! -> args out
                    ! store group result for OAM in array
                    DO m = 1, nmg
                        oam(members(m), ivoigt-3) = oamdg(m)
                    END DO
                    DEALLOCATE(pijA, pijB, dEijdg, oamdg)
                END DO ! loop over 'ivoigt'
                DEALLOCATE( members)
            END DO ! loop over 'ig'
            ! store OAM band-by-band for all occupied states
            DO ivb = 1, nvb
                WRITE(31,TRIM(wformat2)) ivb, (oam(ivb,j), j=1,3)
            END DO
            DEALLOCATE( oam )

            ! spin UP Berry curvature
            ALLOCATE( bcurv(nvb,3) )
            bcurv = 0.0_dp
            ! Loop over groups of degenerate bands in the range of [1:nvb].
            ! Even if a band is not degenerate, it is still considered as
            ! degeneracy of the size 1.
            DO ig = 1, dg_group(nvb)
                nmg = COUNT(dg_group == ig) ! number of group members
                ALLOCATE( members(nmg) )
                ! get bands that are members of the group
                m = 0
                DO n = 1, nb
                    IF (dg_group(n) == ig) THEN
                        m = m + 1
                        members(m) = n
                    ELSE IF (m > 0) THEN
                        ! break the loop once first condition is not true any more,
                        ! which marks the end of a degenerate block
                        EXIT
                    END IF
                END DO
                ! loop over Voigt indices
                DO ivoigt = 4, 6
                    ! handle Voigt notations
                    SELECT CASE (ivoigt)
                        CASE (4); alpha = 2; beta = 3 ! 4=yz
                        CASE (5); alpha = 1; beta = 3 ! 5=xz
                        CASE (6); alpha = 1; beta = 2 ! 6=xy
                    END SELECT
                    ! degenerate group 1st and last band
                    idg1 = members(1)
                    idg2 = members(nmg)
                    IF (idg1 .GT. idg2) THEN
                        WRITE(*,'(2(A,I0))') 'idg1=', idg1, ' idg2=', idg2
                        ERROR STOP 'Wrong boundaries of the degenerate block'
                    END IF
                    ! allocate group-specific arrays
                    ALLOCATE( bcurvdg(nmg), pijA(nmg,nb), pijB(nb,nmg), &
                        dEijdg(nmg,nb))
                    ! get group-specific matrix elements and energies
                    DO m = 1, nmg
                        pijA(m,:) = pijUP(alpha, members(m), :)
                        pijB(:,m) = pij(beta, :, members(m))
                        dEijdg(m,:) = dEij(members(m), :)
                    END DO
                    ! solve degenerate PT problem
                    CALL degenbc(nb, idg1, idg2, & ! <- args in
                        pijA, pijB, dEijdg, & ! <- args in 
                        bcurvdg) ! -> args out
                    ! store group result for the Berry curvature in array
                    DO m = 1, nmg
                        bcurv(members(m), ivoigt-3) = bcurvdg(m)
                    END DO
                    DEALLOCATE(pijA, pijB, dEijdg, bcurvdg)
                END DO ! loop over 'ivoigt'
                DEALLOCATE( members)
            END DO ! loop over 'ig'
            ! store the Berry curvature band-by-band for all occupied states
            DO ivb = 1, nvb
                WRITE(21,TRIM(wformat2)) ivb, (bcurv(ivb,j), j=1,3)
            END DO
            ! store the total Berry curvature for all occupied states
            WRITE(21,TRIM(wformat2)) 0, (SUM(bcurv(:,j)), j=1,3)
            DEALLOCATE( bcurv )

            ! spin DN Berry curvature
            ALLOCATE( bcurv(nvb,3) )
            bcurv = 0.0_dp
            ! Loop over groups of degenerate bands in the range of [1:nvb].
            ! Even if a band is not degenerate, it is still considered as
            ! degeneracy of the size 1.
            DO ig = 1, dg_group(nvb)
                nmg = COUNT(dg_group == ig) ! number of group members
                ALLOCATE( members(nmg) )
                ! get bands that are members of the group
                m = 0
                DO n = 1, nb
                    IF (dg_group(n) == ig) THEN
                        m = m + 1
                        members(m) = n
                    ELSE IF (m > 0) THEN
                        ! break the loop once first condition is not true any more,
                        ! which marks the end of a degenerate block
                        EXIT
                    END IF
                END DO
                ! loop over Voigt indices
                DO ivoigt = 4, 6
                    ! handle Voigt notations
                    SELECT CASE (ivoigt)
                        CASE (4); alpha = 2; beta = 3 ! 4=yz
                        CASE (5); alpha = 1; beta = 3 ! 5=xz
                        CASE (6); alpha = 1; beta = 2 ! 6=xy
                    END SELECT
                    ! degenerate group 1st and last band
                    idg1 = members(1)
                    idg2 = members(nmg)
                    IF (idg1 .GT. idg2) THEN
                        WRITE(*,'(2(A,I0))') 'idg1=', idg1, ' idg2=', idg2
                        ERROR STOP 'Wrong boundaries of the degenerate block'
                    END IF
                    ! allocate group-specific arrays
                    ALLOCATE( bcurvdg(nmg), pijA(nmg,nb), pijB(nb,nmg), &
                        dEijdg(nmg,nb))
                    ! get group-specific matrix elements and energies
                    DO m = 1, nmg
                        pijA(m,:) = pijDN(alpha, members(m), :)
                        pijB(:,m) = pij(beta, :, members(m))
                        dEijdg(m,:) = dEij(members(m), :)
                    END DO
                    ! solve degenerate PT problem
                    CALL degenbc(nb, idg1, idg2, & ! <- args in
                        pijA, pijB, dEijdg, & ! <- args in 
                        bcurvdg) ! -> args out
                    ! store group result for the Berry curvature in array
                    DO m = 1, nmg
                        bcurv(members(m), ivoigt-3) = bcurvdg(m)
                    END DO
                    DEALLOCATE(pijA, pijB, dEijdg, bcurvdg)
                END DO ! loop over 'ivoigt'
                DEALLOCATE( members)
            END DO ! loop over 'ig'
            ! store the Berry curvature band-by-band for all occupied states
            DO ivb = 1, nvb
                WRITE(22,TRIM(wformat2)) ivb, (bcurv(ivb,j), j=1,3)
            END DO
            ! store the total Berry curvature for all occupied states
            WRITE(22,TRIM(wformat2)) 0, (SUM(bcurv(:,j)), j=1,3)
            DEALLOCATE( bcurv )

            ! spin UP-DN Berry curvature (Omega^z)
            ALLOCATE( bcurv(nvb,3) )
            bcurv = 0.0_dp
            ! Loop over groups of degenerate bands in the range of [1:nvb].
            ! Even if a band is not degenerate, it is still considered as
            ! degeneracy of the size 1.
            DO ig = 1, dg_group(nvb)
                nmg = COUNT(dg_group == ig) ! number of group members
                ALLOCATE( members(nmg) )
                ! get bands that are members of the group
                m = 0
                DO n = 1, nb
                    IF (dg_group(n) == ig) THEN
                        m = m + 1
                        members(m) = n
                    ELSE IF (m > 0) THEN
                        ! break the loop once first condition is not true any more,
                        ! which marks the end of a degenerate block
                        EXIT
                    END IF
                END DO
                ! loop over Voigt indices
                DO ivoigt = 4, 6
                    ! handle Voigt notations
                    SELECT CASE (ivoigt)
                        CASE (4); alpha = 2; beta = 3 ! 4=yz
                        CASE (5); alpha = 1; beta = 3 ! 5=xz
                        CASE (6); alpha = 1; beta = 2 ! 6=xy
                    END SELECT
                    ! degenerate group 1st and last band
                    idg1 = members(1)
                    idg2 = members(nmg)
                    IF (idg1 .GT. idg2) THEN
                        WRITE(*,'(2(A,I0))') 'idg1=', idg1, ' idg2=', idg2
                        ERROR STOP 'Wrong boundaries of the degenerate block'
                    END IF
                    ! allocate group-specific arrays
                    ALLOCATE( bcurvdg(nmg), pijA(nmg,nb), pijB(nb,nmg), &
                        dEijdg(nmg,nb))
                    ! get group-specific matrix elements and energies
                    DO m = 1, nmg
                        pijA(m,:) = pijUP(alpha, members(m), :) - &!...
                            pijDN(alpha, members(m), :)
                        pijB(:,m) = pij(beta, :, members(m))
                        dEijdg(m,:) = dEij(members(m), :)
                    END DO
                    ! solve degenerate PT problem
                    CALL degenbc(nb, idg1, idg2, & ! <- args in
                        pijA, pijB, dEijdg, & ! <- args in 
                        bcurvdg) ! -> args out
                    ! store group result for the Berry curvature in array
                    DO m = 1, nmg
                        bcurv(members(m), ivoigt-3) = bcurvdg(m)
                    END DO
                    DEALLOCATE(pijA, pijB, dEijdg, bcurvdg)
                END DO ! loop over 'ivoigt'
                DEALLOCATE( members)
            END DO ! loop over 'ig'
            ! store the Berry curvature band-by-band for all occupied states
            DO ivb = 1, nvb
                WRITE(23,TRIM(wformat2)) ivb, (bcurv(ivb,j), j=1,3)
            END DO
            ! store the total Berry curvature for all occupied states
            WRITE(23,TRIM(wformat2)) 0, (SUM(bcurv(:,j)), j=1,3)
            DEALLOCATE( bcurv )
        END IF ! WIEN2k spin-resolved Berry curvature

        DEALLOCATE( pij, dEij, dg_group ) ! all k-point specific variables
        IF ( wien2k .AND. spinor ) THEN
            DEALLOCATE( pijUP, dEijUP, pijDN, dEijDN ) ! up/dn components
        END IF
        
        !! Output progress to screen
        
        IF (wien2k) THEN
            WRITE(*,'(A,I0,A,I0,A,I3,A)') & !...
                ' KP: ', ikpt, ' bands: ', nb, & !...
                ' progress: ', INT(100*iline/nltot), '%'
        ELSE ! VASP
            IF (ispin==1 .AND. nstot==2) THEN
                WRITE(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2), '%'
            ELSE IF (ispin==2 .AND. nstot==2) THEN
                WRITE(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2 + 50), '%'
            ELSE
                WRITE(*,'(A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot), '%'
            END IF
        END IF
        
        !! Check if the end of file mommat is reached
        
        IF ( wien2k ) THEN
            IF (iline == nltot) THEN ! end of mommat file, exit WHILE loop
                fmommatend = .true.
            END IF
        ELSE ! VASP
            IF (ikpt == nktot) THEN ! end of WAVEDER file, exit WHILE loop
                fmommatend = .true.
            END IF
        END IF

    END DO ! loop over k-points
    
    !! Close output files
    
    CLOSE (2) ! output file with Berry curvature
    IF ( wien2k .AND. spinor ) THEN
        CLOSE (21) ! spin-resolved UP Berry curvature
        CLOSE (22) ! spin-resolved DN Berry curvature
        CLOSE (23) ! spin-resolved UP-DN Berry curvature
    END IF
    CLOSE (3) ! output file with OAM
    IF ( wien2k .AND. spinor ) THEN
        CLOSE (31) ! spin-resolved UP OAM
    END IF

END DO ! loop over spins

!! Close input files

CLOSE (1) ! input  file 1 (case.mommat2 or WAVEDER)
IF ( wien2k .AND. spinor ) THEN
    CLOSE (11) ! case.mommat2up
    CLOSE (12) ! case.mommat2dn
END IF

WRITE(*,*) 'Summary of the output:'
WRITE(*,*) '(1) Components of the Berry curvature tensor are stored in the file'
WRITE(*,*) '    ' // TRIM(fnameout2)
WRITE(*,*) '    See the file header for the description'
WRITE(*,*) '(2) Components of the orbital angular momentum tensor are stored in the file'
WRITE(*,*) '    ' // TRIM(fnameout3)
IF ( wien2k .AND. spinor ) THEN
    WRITE(*,*) '(3) Spin Berry curvature (WIEN2k only) is stored in'
    WRITE(*,*) '    ' // TRIM(fnameout21) // ', ' // &
        TRIM(fnameout22) // ', ' // TRIM(fnameout23)
    WRITE(*,*) '(4) Spin OAM sigma_z, UP (WIEN2k only) is stored in'
        WRITE(*,*) '    ' // TRIM(fnameout31)
END IF
WRITE(*,*) ''
WRITE(*,*) 'Suggested reference:'
WRITE(*,'(A,/A,/,A)') '[1] Rubel, O. (2024). BerryCPT: Berry curvature and orbital ', &
    '    angular momentum from perturbation theory [Computer software]. ', &
    '    GitHub. https://github.com/rubel75/berrycpt'
STOP ! end of the main code



!! Help section

911 & ! label for GOTO statement
WRITE(*,*) 'Suggested execution (WIEN2k, no SP, no SOC):'
WRITE(*,*) '$ ./berrycpt mass.mommat2 -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (WIEN2k, with SP, no SOC):'
WRITE(*,*) '$ ./berrycpt mass.mommat2up -nvb 26'
WRITE(*,*) '$ ./berrycpt mass.mommat2dn -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (WIEN2k, with SOC to produce spin-resolved output):'
WRITE(*,*) '$ ./berrycpt mass.mommat2 -nvb 26 -so'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (VASP) with the number of occupied bands:'
WRITE(*,*) '$ ./berrycpt WAVEDER -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (VASP) with the Fermi energy (eV):'
WRITE(*,*) '$ ./berrycpt WAVEDER -efermiev 3.123'
WRITE(*,*) ''
WRITE(*,*) 'Output:'
WRITE(*,*) 'bcurv_ij.dat - contains elements of the Berry curvature tensor.'
WRITE(*,*) 'bcurv_ij-[up,dn].dat - contains elements of the spin-resolved'
WRITE(*,*) '     Berry curvature tensor.'
WRITE(*,*) 'If the files already exist from a previous run, they will be overwritten'
WRITE(*,*) ''
WRITE(*,*) 'Tips:'
WRITE(*,*) '(1): Writing of the mommat file is _not_ default in WIEN2k.'
WRITE(*,*) '     To enable writing, edit the case.inop file and change'
WRITE(*,*) '     OFF to ON in the following line:'
WRITE(*,*) '     ON           ON/OFF   WRITEs MME to unit 4'
WRITE(*,*) '     -^'
WRITE(*,*) '(2): Make sure to get _plenty_ empty bands during SCF.'
WRITE(*,*) '     This requires modification in several input files:'
WRITE(*,*) '     (a) extend "de" in case.in1(c) above 5 Ry'
WRITE(*,'(A)') '          K-VECTORS FROM UNIT:4   -9.0      10.0'//&!...
    '    10   emin / de (emax=Ef+de) / nband'
WRITE(*,*) '         -----------------------------------^'
WRITE(*,*) '     (b) if you do SOC calculation, extend "Emax" in'
WRITE(*,*) '         case.inso up to 5 Ry'
WRITE(*,*) '         -10 5.0                Emin, Emax'
WRITE(*,*) '         -----^'
WRITE(*,*) '     (c) default "Emax" in case.inop is 3 Ry, which should'
WRITE(*,*) '         be OK, but it can be good to test the convergence'
WRITE(*,*) '         and push this parameter up to 5 Ry'
WRITE(*,*) '         -5.0 3.5 9999 Emin, Emax for matrix elements, NBvalMAX'
WRITE(*,*) '         ------^'
WRITE(*,*) '(3): In VASP calculations use LOPTICS = .TRUE., increase'
WRITE(*,*) '     (at least x3) the number NBANDS = XXXX, and disable a finite'
WRITE(*,*) '     differences derivative of the cell-periodic part of'
WRITE(*,*) '     the orbitals LPEAD =.FALSE.'
STOP

!! Error section

912 & ! label for GOTO statement
WRITE(*,*) 'Error detected for the number of input arguments.'
WRITE(*,*) 'There should be 3 or 4 arguments'
WRITE(*,*) 'Suggested execution (WIEN2k, no SP, no SOC):'
WRITE(*,*) '$ ./berrycpt mass.mommat2 -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (WIEN2k, with SP, no SOC):'
WRITE(*,*) '$ ./berrycpt mass.mommat2up -nvb 26'
WRITE(*,*) '$ ./berrycpt mass.mommat2dn -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (WIEN2k, with SOC):'
WRITE(*,*) '$ ./berrycpt mass.mommat2 -nvb 26 -so'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (VASP) with the number of occupied bands:'
WRITE(*,*) '$ ./berrycpt WAVEDER -nvb 26'
WRITE(*,*) ' '
WRITE(*,*) 'Suggested execution (VASP) with the Fermi energy (eV):'
WRITE(*,*) '$ ./berrycpt WAVEDER -efermiev 3.123'
ERROR STOP

END PROGRAM berrycpt
