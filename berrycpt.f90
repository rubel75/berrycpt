PROGRAM berrycpt
! includes subprograms:
!     dgenen.f90 -- Group degenerate states
!     eigvz.f90 -- Solve a complex eigenvalue problem for a Hermitian matrix
!     read_mommat_nb.f90 -- Determine number of bands in a k-point block of 
!       mommat file
!     read_mommat_pij.f90 -- Read momentum matrix elements <i|p_a|j> (a=x,y,z) 
!       and energy differences E_i - E_j in a k-point block of mommat file
!     read_numlines.f90 -- Read number of lines in a file
!
! (c) Oleg Rubel, Mar 2021
!
! Execution:
!   $ ./berrycpt arg1 [arg2]
! arg1 - input mommat file name
! arg2 - (optional) degeneracy energy tolerance [Ha]
!        max dE for 2 states to be considered as degenerate
!        default value is 1.0e-6 Ha
!
! Output:
! bcurv_ij.dat - contains elements of the Berry curvature tensor.
! If the files exist from a previous run, they will be removed
!
! Tips:
! (1): Writing of the mommat file is _not_ default in WIEN2k.
!      To enable writing, edit the case.inop file and change
!      OFF to ON in the following line:
!      ON           ON/OFF   writes MME to unit 4
!      -^
! (2): Make sure to get _plenty_ empty bands during SCF.
!      This requires modification in several input files:
!      (a) extend "de" in case.in1(c) above 5 Ry
!          K-VECTORS FROM UNIT:4   -9.0      10.0    10   emin / de (emax=Ef+de) / nband
!          -----------------------------------^
!      (b) if you do SOC calculation, extend "Emax" in
!          case.inso up to 5 Ry
!          -10 5.0                Emin, Emax
!          -----^
!      (c) default "Emax" in case.inop is 3 Ry, which should
!          be OK, but it can be good to test the convergence
!          and push this parameter up to 5 Ry
!          -5.0 3.5 9999 Emin, Emax for matrix elements, NBvalMAX
!          ------^
! (3): In VASP calculations use LOPTICS = .TRUE., increase (at least x3) 
!      the number NBANDS = XXXX, and disable a finite differences derivative 
!      of the cell-periodic part of the orbitals LPEAD =.FALSE. Tests show 
!      that m* values calculated with LPEAD =.TRUE. make no sense.

!! Variables

USE OMP_LIB

implicit none
CHARACTER(len=256) :: &
    arg1, arg2, arg3, & ! command line input arguments
    fnameinp, fnameout2, & ! input/output file names
    wformat2, & ! format for writing/reading data
    charspin ! spin component for spin-polarized calculations
INTEGER :: &
    nltot, & ! total number of lines in mommat file
    nktot, & ! total number of k-points in WAVEDER file
    nstot, & ! total number of spins in WAVEDER file
    iline, & ! current line number during reading of mommat file
    ispin, & ! current spin (1/2)
    nb, & ! number of bands for a current k-point
    nbb, & ! number of band-to-band transition 
    nbcder, & ! the max number of bands for which m* is calculated
              ! In WIEN2k we set it = to the total number of bands
              ! In VASP it is set in WAVEDER file as NBANDS_CDER
    ivb, i, j, ikpt, & ! counters
    n, & ! band indices
    nvb, & ! number of occupied bands [1:nvb]
    ivoigt, & ! Voigt index (1..6)
    alpha, beta, & ! Cartesian directions 1,2,3 = x,y,z
    ierr ! error code
INTEGER, ALLOCATABLE :: &
    nbocck(:,:) ! k- and spin-resolved number of occupied bands
REAL(kind=4) :: &
    dE, & ! energy difference [Ha]
    efermi ! Fermi energy in [eV]
REAL(kind=4), ALLOCATABLE :: &
    dEij(:,:), & ! energy differences E_i-E_j [Ha]
    dEijks(:,:,:,:), & ! k- and spin-dependent energy differences E_i-E_j [Ha]
    bcurv(:,:), & ! array to store Berry curvature
    omega, & ! Berry curvature (intermediate)
    p2 ! product of momentum matrix elements
COMPLEX(kind=4), ALLOCATABLE :: &
    pij(:,:,:), & ! momentum matrix elements [at.u.]
    pijks(:,:,:,:,:) ! k- and spin-dependent momentum matrix elements [at.u.]
LOGICAL :: &
    fmommatend, & ! end of mommat file
    file_exists, &
    wien2k ! true if this is a WIEN2k calculation

!! Get command line input arguments

DO i = 1, iargc() ! loop through all command line arguments to search for help
    CALL GETARG(i,arg1)
    IF ( TRIM(arg1)=='-h' .or.  TRIM(arg1)=='--h' .or. & !...
            TRIM(arg1)=='-help' .or. TRIM(arg1)=='--help') THEN
        GOTO 911 ! print help and STOP
    END IF
END DO
write(*,'(A,I0)') ' Detected input arguments = ', iargc()
IF (iargc() .ne. 3) THEN ! check number of input arguments
    write(*,'(A)') ' Expected 3 input arguments'
    GOTO 912 ! print error and STOP
END IF
! 1st argument
CALL GETARG(1,arg1) ! mommat file name
write (*,*) 'Input mommat file = ', TRIM(arg1)
fnameinp = TRIM(arg1)
! 2nd argument
CALL GETARG(2,arg2) ! switch for the number of valence bands [-nvb]
IF (TRIM(arg2)=='-nvb') THEN
    ! 3rd argument can be the number of occupied bands
    CALL GETARG(3,arg3)
    read(arg3,*,IOSTAT=ierr) nvb
    IF (ierr /= 0) THEN ! input error for 3rd argument
        write(*,*) '3rd argument is "', TRIM(arg3), '"'
        write(*,*) 'Error detected for the 3rd input argument (must be a number of occupied bands)'
        GOTO 912 ! print error and STOP
    END IF
    write (*,'(A,I0)') & !...
        ' Number of occupied bands = ', nvb
ELSE IF (TRIM(arg2)=='-efermiev') THEN
    ! 3rd argument can be the Fermi energy, which will be used to determine
    ! spin and k-specific number of occupied bands
    CALL GETARG(3,arg3)
    read(arg3,*,IOSTAT=ierr) efermi
    IF (ierr /= 0) THEN ! input error for 3rd argument
        write(*,*) '3rd argument is "', TRIM(arg3), '"'
        write(*,*) 'Error detected for the 3rd input argument (must be a Fermi energy)'
        GOTO 912 ! print error and STOP
    END IF
    write (*,'(A,F10.6)') & !...
        ' Fermi energy (eV) = ', nvb
ELSE ! impossible
    write (*,*) 'The second argument is "', TRIM(arg2), &
        '", while expected "-nvb" or "-efermiev"'
    GOTO 912 ! print error and STOP
END IF

!! WIEN2k or VASP?

IF ( TRIM(fnameinp) == 'WAVEDER' .OR. TRIM(fnameinp) == 'WAVEDERF' ) THEN ! VASP
    wien2k = .false.
    write (*,*) 'Assume VASP calculation'
    ! change to a binary file in case it is pointed at the formatted 
    ! file WAVEDER
    fnameinp = 'WAVEDER'
    write(*,*) 'Assumed VASP calculation based on the input file name.'
ELSE ! WIEN2k (default)
    wien2k = .true.
    write(*,*) 'Assumed WIEN2k calculation based on the input file name.'
    write(*,*) '(If you would like to read VASP file, the input file should'
    write(*,*) 'have the exact name WAVEDER.)'
END IF

!! Determine number of lines in mommat file & read VASP WAVEDER header

! check if the file exists
INQUIRE(FILE=fnameinp, EXIST=file_exists)
IF ( file_exists ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' was found.'
ELSE IF ( .not.(file_exists) ) THEN
    write(*,*) 'The input file ', TRIM(fnameinp), ' does not exist. Exiting'
    STOP
END IF
IF ( file_exists .AND. (.not.(wien2k)) ) THEN ! check EIGENVAL for VASP
    INQUIRE(FILE='EIGENVAL', EXIST=file_exists)
    IF ( .not.(file_exists) ) THEN
        write(*,*) 'The file EIVENVAL is also required, but it ', &
            'does not exist. Exiting'
        STOP
    ELSE
        write(*,*) 'The file EIVENVAL is also required. ', &
            'It is found.'
    END IF
END IF
IF (wien2k) THEN
    CALL read_numlines(fnameinp, 1, & ! <- args in 
                nltot) ! -> args out
    write (*,'(A,I0)') '  number of lines in mommat file = ', nltot
    IF (nltot < 10) THEN
        write(*,*) 'The file ', TRIM(fnameinp), ' is too short ', & !...
            'and most likely useless. Stopping'
        STOP
    END IF
    ! assume 1 spin since case.momat2up and case.momat2dn files should be
    ! read one after the other any ways
    nstot = 1 
ELSE ! VASP
    CALL read_numlines_vasp(fnameinp, 1, & ! <- args in 
        nstot, nktot, nbcder, nb) ! -> args out
    write (*,'(A,I0)') '  number of spins in WAVEDER file = ', nstot
    write (*,'(A,I0)') '  number of k-points in WAVEDER file = ', nktot
    write (*,'(A,I0)') '  smaller number of bands in WAVEDER file = ', nbcder
    write (*,'(A,I0)') '  number of bands per k-point in EIGENVAL file = ', &
        nb
    ! Memory estimate: 
    ! 4 bytes * 2 (complex) * 3 (3D x,y,z) * 10e9 (Bytes -> GB)
    write (*,'(A,1X,F5.1,1X,A)') ' Memory required to store matrix elements'//&
        ' from WAVEDER file', nstot*nktot*(nb**2.)*4*2*3/10.0**9., 'GB'
    write (*,'(A,1X,F5.1,1X,A)') ' Memory required to store dE_ij'//&
        ' from EIGENVAL file', nstot*nktot*(nb**2.)*4*1/10.0**9., 'GB'
    write (*,'(A,1X,F5.1,1X,A)') ' It will take additional', &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB overhead'//&
        ' to run the calculation'
    write (*,'(A,1X,F5.1,1X,A)') ' Overall, you will need at least', &
        nstot*nktot*(nb**2.)*4*2*3/10.0**9. + &
        nstot*nktot*(nb**2.)*4*1/10.0**9. + &
        1*1*(nb**2.)*4*(2*3+1)/10.0**9., 'GB (+ 20% incidental) of RAM'//&
        ' to run the calculation'
END IF
! open input mommat/WAVEDER file for reading (file ID=1)
IF (wien2k) THEN
    OPEN (1, file = TRIM(fnameinp), status = 'old')
ELSE ! VASP
    OPEN (1, file = TRIM(fnameinp), form = 'unformatted', status = 'old')
END IF

!! Read VASP WAVEDER and EIGENVAL files to determine matrix elements

IF (.not.(wien2k)) THEN
    ALLOCATE( pijks(3,nb,nb,nktot,nstot), dEijks(nb,nb,nktot,nstot) )
    ALLOCATE( nbocck(nktot,nstot) )
    CALL read_mommat_pij_vasp (1, nstot, nktot, nbcder, nb, efermi, & ! <- args in
            pijks, dEijks, nbocck) ! -> args out
END IF

! Loop over spins (for VASP only)
! In WIEN2k nstot=1 since case.momat2up and case.momat2dn files should be
! read one after the other any ways
write(*,*) 'Entering the main loop...'
DO ispin = 1, nstot

    !! Spin suffix
    
    IF (wien2k) THEN
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
    write(2,'(A)') '# This file is generated by berrycpt'
    write(2,'(A)') '# the output contains components of the Berry curvature tensor that'
    write(2,'(A)') '# are grouped by k-point index and then by the band index.'
    write(2,'(A)') '# The last entry for each k-point with the band index "0" is the'
    write(2,'(A)') '# total Berry curvature of all occupied bands.'
    write(2,'(A)') '# Columns correspond to Cartesian directions for Omega_ij'
    write(2,'(A)') '# band 1=xx;     2=yy;      3=zz;     4=yz;'//&!...
        '       5=xz;      6=xy'

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
            ! the max number of bands for which m* is calculated
            ! here (WIEN2k) we set it = to the total number of bands
            nbcder = nb
        END IF

        !! Read <b_i|p_a|b_j> from mommat file, a=1,2,3 (x,y,z)

        ALLOCATE( pij(3,nb,nb), dEij(nb,nb) )
        IF (wien2k) THEN
            nbb = (nb+nb**2)/2 ! number of band-to-band transitions
            CALL read_mommat_pij (1, nb, nbb, & ! <- args in
                    iline, & ! <-> args in-out
                    pij, dEij) ! -> args out
        ELSE ! VASP
            pij = pijks(:,:,:,ikpt,ispin)
            dEij = dEijks(:,:,ikpt,ispin)
        END IF

        !! VASP k- and spin-specific number of occupied bands

        IF (.not.(wien2k)) THEN
            nvb = nbocck(ikpt,ispin)
        END IF

        !! Write information about the current k-point

        write(2,'(A,I0,1X,A,I0,1X,A,I0)') & !...
            '# KP: ', ikpt, 'NVB: ', nvb, 'NEMAX: ', nb
        
        !! preapare output formats for Berry curvatures
        
        write(wformat2,'(I0)') nbcder ! make a character of the length 'nbcder'
        ! format line to write Berry curvatures
        write(wformat2,'(I0)') LEN(TRIM(wformat2))
        wformat2 = '(I' // TRIM(wformat2) // ',1X,5(es10.3,1X),es10.3)'
            
        !! Loop through blocks of occupied bands

        ! array to store Berry curvatures; size 6 is because of 6 Voigt indices
        ALLOCATE( bcurv(nvb,6) )
        bcurv = 0.0
        DO ivb = 1,nvb
            ! loop over Voigt indecies: 1=xx; 2=yy; 3=zz; 4=yz; 5=xz; 6=xy
            DO ivoigt = 1,6
                ! handle Voigt notations
                SELECT CASE (ivoigt)
                    CASE (1) ! 1=xx
                        alpha = 1
                        beta = 1
                    CASE (2) ! 2=yy
                        alpha = 2
                        beta = 2
                    CASE (3) ! 3=zz
                        alpha = 3
                        beta = 3
                    CASE (4) ! 4=yz
                        alpha = 2
                        beta = 3
                    CASE (5) ! 5=xz
                        alpha = 1
                        beta = 3
                    CASE (6) ! 6=xy
                        alpha = 1
                        beta = 2
                    CASE DEFAULT
                END SELECT
                DO n = nvb+1, nb
                    ! extract complex part of the product
                    ! take into account that i*5i = -5
                    p2 = (-1.0) * AIMAG( pij(alpha,ivb,n)*pij(beta,n,ivb) - & !...
                        pij(beta,ivb,n)*pij(alpha,n,ivb) )
                    dE = dEij(ivb,n)
                    omega = p2/(dE*dE)
                    ! make sure omega is finite (not NaN and not Inf)
                    IF (omega /= omega .or. abs(omega) > HUGE(abs(omega))) THEN
                        write(*,*) 'ikpt =', ikpt
                        write(*,*) 'ivb =', ivb
                        write(*,*) 'n =', n
                        write(*,*) 'alpha =', alpha
                        write(*,*) 'beta =', beta
                        write(*,*) 'pij(alpha,ivb,n) =', &
                            pij(alpha,ivb,n)
                        write(*,*) 'pij(beta,n,ivb) =', &
                            pij(beta,n,ivb)
                        write(*,*) 'pij(beta,ivb,n) =', &
                            pij(beta,ivb,n)
                        write(*,*) 'pij(alpha,n,ivb) =', &
                            pij(alpha,n,ivb)
                        write(*,*) 'dEij(ivb,n) =', dEij(ivb,n)
                        write(*,*) 'omega = ', omega
                        write(*,*) 'dE = ', dE
                        write(*,*) 'p2 = ', p2
                        STOP 'Error: omega is not finite'
                    END IF
                    ! update Mnm array
                    bcurv(ivb,ivoigt) = bcurv(ivb,ivoigt) + omega
                END DO ! loop over 'n'
            END DO ! loop over 'ivoigt'
        END DO ! loop over 'ivb'
        ! store the Berry curvature band-by-band for all occupied states
        DO ivb = 1, nvb
            write(2,TRIM(wformat2)) ivb, (bcurv(ivb,j), j=1,6)
        END DO
        ! store the total Berry curvature for all occupied states
        write(2,TRIM(wformat2)) 0, (SUM(bcurv(:,j)), j=1,6)
        DEALLOCATE( bcurv )
        DEALLOCATE( pij, dEij ) ! all k-point specific variables
        
        !! Output progress to screen
        
        IF (wien2k) THEN
            write(*,'(A,I0,A,I0,A,I3,A)') & !...
                ' KP: ', ikpt, ' bands: ', nb, & !...
                ' progress: ', INT(100*iline/nltot), '%'
        ELSE ! VASP
            IF (ispin==1 .AND. nstot==2) THEN
                write(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2), '%'
            ELSE IF (ispin==2 .AND. nstot==2) THEN
                write(*,'(A,I0,A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' spin: ', ispin, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot/2 + 50), '%'
            ELSE
                write(*,'(A,I0,A,I0,A,I3,A)') & !...
                    ' KP: ', ikpt, ' bands: ', nb, & !...
                    ' progress: ', INT(100*ikpt/nktot), '%'
            END IF
        END IF
        
        !! Check if the end of file mommat is reached
        
        IF (wien2k) THEN
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
    
    CLOSE (2) ! output file 1
    
END DO ! loop over spins

CLOSE (1) ! input  file 1

write(*,*) 'Summary of the output:'
write(*,*) '(1) Components of the Berry curvature tensor are stored in the file'
write(*,*) '    ', TRIM(fnameout2)
write(*,*) '    See the file header for the description'
write(*,*) ''
write(*,*) 'Suggested reference:'
write(*,*) '[1] O. Rubel, F. Tran, X. Rocquefelte, and P. Blaha "Perturbation'
write(*,*) '    approach to ab initio effective mass calculations"'
write(*,*) '    Comp. Phys. Commun. 261, 107648 (2021).'
write(*,*) '    https://doi.org/10.1016/j.cpc.2020.107648'
STOP ! end of the main code



!! Help section

911 & ! label for GOTO statement
write(*,*) 'Execution:'
write(*,*) '  $ ./berrycpt arg1 [arg2]'
write(*,*) 'arg1 - input case.mommat2 (WIEN2k) or WAVEDER (VASP) file name'
write(*,*) 'arg2 - (optional) degeneracy energy tolerance [Ha]'
write(*,*) '       max dE for 2 states to be considered as degenerate'
write(*,*) '       default value is 1.0e-6 Ha'
write(*,*) ''
write(*,*) 'Output:'
write(*,*) 'bcurv_ij.dat - contains elements of the Berry curvature tensor.'
write(*,*) 'If the files exist from a previous run, they will be removed'
write(*,*) ''
write(*,*) 'Tips:'
write(*,*) '(1): Writing of the mommat file is _not_ default in WIEN2k.'
write(*,*) '     To enable writing, edit the case.inop file and change'
write(*,*) '     OFF to ON in the following line:'
write(*,*) '     ON           ON/OFF   writes MME to unit 4'
write(*,*) '     -^'
write(*,*) '(2): Make sure to get _plenty_ empty bands during SCF.'
write(*,*) '     This requires modification in several input files:'
write(*,*) '     (a) extend "de" in case.in1(c) above 5 Ry'
write(*,'(A)') '          K-VECTORS FROM UNIT:4   -9.0      10.0'//&!...
    '    10   emin / de (emax=Ef+de) / nband'
write(*,*) '         -----------------------------------^'
write(*,*) '     (b) if you do SOC calculation, extend "Emax" in'
write(*,*) '         case.inso up to 5 Ry'
write(*,*) '         -10 5.0                Emin, Emax'
write(*,*) '         -----^'
write(*,*) '     (c) default "Emax" in case.inop is 3 Ry, which should'
write(*,*) '         be OK, but it can be good to test the convergence'
write(*,*) '         and push this parameter up to 5 Ry'
write(*,*) '         -5.0 3.5 9999 Emin, Emax for matrix elements, NBvalMAX'
write(*,*) '         ------^'
write(*,*) '(3): In VASP calculations use LOPTICS = .TRUE., increase'
write(*,*) '     (at least x3) the number NBANDS = XXXX, and disable a finite'
write(*,*) '     differences derivative of the cell-periodic part of'
write(*,*) '     the orbitals LPEAD =.FALSE. Tests show that m* values'
write(*,*) '     calculated with LPEAD =.TRUE. make no sense.'
STOP 0

!! Error section

912 & ! label for GOTO statement
write(*,*) 'Error detected for the number of input arguments.'
write(*,*) 'There should be 3 arguments'
write(*,*) 'Suggested execution (WIEN2k):'
write(*,*) '$ ./berrycpt mass.mommat2 -nvb 26'
write(*,*) ' '
write(*,*) 'Suggested execution (VASP) with the number of occupied bands:'
write(*,*) '$ ./berrycpt WAVEDER -nvb 26'
write(*,*) ' '
write(*,*) 'Suggested execution (VASP) with the Fermi energy (eV):'
write(*,*) '$ ./berrycpt WAVEDER -efermiev 3.123'
STOP 1

END PROGRAM berrycpt
