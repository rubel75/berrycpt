PROGRAM berrycpt
! BerryCPT calculates band-resolved Berry curvature, orbital angular
! momentum, and generalized orbital-angular-momentum matrices from
! momentum or wave-function-derivative matrix elements.
!
! Supported input:
!
!   WIEN2k
!     - one case.mommat2 file for ordinary calculations
!     - three matrix-element files, total, spin-up, and spin-down,
!       for spin-resolved calculations with spin-orbit coupling
!     - case.energy, case.energyup/dn, case.energyso, or case.energysoup/dn
!     - case.weight or case.weightup/dn when DFT occupations
!       are requested
!
!   VASP
!     - WAVEDER or WAVEDERF
!     - EIGENVAL from the same calculation
!
! Exactly one occupation mode must be selected:
!
!   --efermiev EF
!       Step-function occupations using the Fermi energy EF in eV.
!
!   --efermiry EF
!       Step-function occupations using the Fermi energy EF in Ry.
!
!   --dftocc
!       Occupations supplied by the DFT calculation. WIEN2k also
!       requires --occfile. VASP occupations are read from EIGENVAL.
!
! Use
!
!   ./berrycpt --help
!
! for the complete command-line syntax.
!
! Examples:
!
!   VASP:
!
!     ./berrycpt WAVEDER --enefile EIGENVAL --efermiev 6.52
!
!     ./berrycpt WAVEDER --enefile EIGENVAL --dftocc
!
!   WIEN2k without spin-resolved matrix elements:
!
!     ./berrycpt case.mommat2 --enefile case.energy --efermiry 0.73
!
!   WIEN2k with spin-resolved matrix elements and SOC:
!
!     ./berrycpt case.mommat2 case.mommat2up case.mommat2dn \
!         --enefile case.energyso --dftocc --occfile case.weight
!
! Calculated quantities:
!
!   bcurv_ij[-up/-dn].dat
!       Ordinary Berry curvature resolved by band and k-point.
!
!   bcurv_ij-up.dat
!   bcurv_ij-dn.dat
!       Spin-up and spin-down contributions to the Berry curvature.
!       These files are produced when spin-resolved WIEN2k matrix
!       elements are supplied.
!
!   bcurv_ij-up-dn.dat
!       Sigma_z-normalized spin Berry curvature.
!
!   oam_ij[-up/-dn].dat
!       Orbital angular momentum resolved by band and k-point.
!
!   oam_ij-sigma_z-up.dat
!       Spin-up-projected orbital angular momentum.
!
!   goam_ij_nm[-up/-dn].dat
!       Complex generalized orbital-angular-momentum matrices in the
!       original DFT eigenstate basis.
!
! The pseudovector component order is
!
!   yz, zx, xy
!
! corresponding to Voigt indices 4, 5, and 6. Berry curvature is
! reported in bohr^2. Orbital angular momentum is reported in hbar.
!
! Existing output files are replaced. Each output file contains a
! detailed header describing its columns, k-point organization,
! degenerate-block labels, and matrix layout. At normal termination,
! berrycpt prints a summary of all output files that were created.
!
! Tips for generating matrix elements
!
! (1) WIEN2k does not write the mommat file by default.
!
!     In case.inop, enable writing of momentum matrix elements by
!     changing OFF to ON:
!
!         ON           ON/OFF   WRITEs MME to unit 4
!         -^
!
! (2) Include enough empty bands and a sufficiently large energy
!     range. Berry curvature and orbital angular momentum contain
!     sums over intermediate states, so their convergence must be
!     tested with respect to the upper band and energy limits.
!
!     In case.in1 or case.in1c, increase and convergence-test "de":
!
!         K-VECTORS FROM UNIT:4   -9.0      10.0    10   emin / de / nband
!                                  ----------^
!
!     For calculations with SOC, ensure that Emax in case.inso covers
!     all states required in the intermediate-state sums:
!
!         -10 5.0                Emin, Emax
!           ---^
!
!     Ensure that the matrix-element range in case.inop is at least
!     as large as the corresponding eigenvalue and SOC ranges:
!
!         -5.0 3.5 9999         Emin, Emax for matrix elements, NBvalMAX
!                ---^
!
!     The numerical values above are examples, not universal
!     convergence settings.
!
! (3) For VASP, generate WAVEDER using:
!
!         LOPTICS = .TRUE.
!         LPEAD   = .FALSE.
!
!     LPEAD = .FALSE. is required because berrycpt uses the full
!     interband derivative matrix stored in WAVEDER.
!
!     Include enough empty bands:
!
!         NBANDS = XXXX
!
!     Approximately three times the number of occupied bands can be
!     used as an initial estimate, but NBANDS must be convergence-tested.
!     WAVEDER and EIGENVAL must come from the same VASP calculation.
!
!     GW PAW potentials _MUST BE USED_ for VASP calculations because
!     they are optimized for accurately representing unoccupied states
!     far above the Fermi level and provide a more suitable PAW
!     projector basis for the high-energy empty states entering the
!     intermediate-state sums.
!
! Copyright (c) 2024-2026 Oleg Rubel

!! Variables

USE write_progress_mod, ONLY: initialize_progress, write_progress
USE precision_mod, ONLY: sp, dp
USE command_line_args_mod, ONLY: command_line_args
USE validate_input_files_mod, ONLY: validate_input_files
USE read_waveder_header_mod, ONLY: read_waveder_header
USE read_eigenvalues_wien2k_mod, ONLY: read_eigenvalues_wien2k
USE read_eigenvalues_vasp_mod, ONLY: read_eigenvalues_vasp
USE read_occupations_wien2k_mod, ONLY: read_occupations_wien2k
USE construct_occupations_mod, ONLY: construct_occupations
USE filename_contains_energyso_mod, ONLY: filename_contains_energyso
USE find_degenerate_groups_mod, ONLY: find_degenerate_groups
USE set_spin_suffix_mod, ONLY: set_spin_suffix
USE open_mommat_files_mod, ONLY: open_mommat_files, close_mommat_files
USE read_mommat_pij_vasp_mod, ONLY: read_mommat_pij_vasp
USE open_output_files_mod, ONLY: &
    open_output_files, close_output_files, write_output_summary, &
    output_units_type
USE read_mommat_pij_wien2k_kpoint_mod, ONLY: read_mommat_pij_wien2k_kpoint
USE calculate_bcurv_kpoint_mod, ONLY: calculate_bcurv_kpoint
USE write_bcurv_kpoint_mod, ONLY: write_bcurv_kpoint
USE calculate_oam_kpoint_mod, ONLY: calculate_oam_kpoint
USE write_oam_kpoint_mod, ONLY: write_oam_kpoint
USE calculate_goam_kpoint_mod, ONLY: calculate_goam_kpoint
USE write_goam_kpoint_mod, ONLY: write_goam_kpoint
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan

IMPLICIT NONE

! Output-file units for the current spin channel
TYPE(output_units_type) :: output_units

! Input and derived file names
CHARACTER(LEN=:), ALLOCATABLE :: &
    fnameinp, &   ! primary matrix-element file
    fnameinpUP, & ! spin-up matrix-element file (WIEN2k only)
    fnameinpDN, & ! spin-down matrix-element file (WIEN2k only)
    fnameene, &   ! eigenvalue file
    fnameocc, &   ! occupation file used with --dftocc (WIEN2k only)
    fnamebase     ! basename of the primary matrix-element file

! Suffix appended to spin-dependent output file names
CHARACTER(LEN=3) :: charspin

! Global calculation dimensions
INTEGER :: &
    nstot, &      ! total number of spin channels
    nktot, &      ! total number of k-points
    nb, &         ! maximum number of available bands
    nbcder = 0    ! WAVEDER target-band limit; ignored for WIEN2k

! Main-loop indices and current k-point dimensions
INTEGER :: &
    ispin, &      ! current spin-channel index
    ikpt, &       ! current k-point index
    i, &          ! target-band index used to construct pijSPIN
    j, &          ! transition-band index used to construct pijSPIN
    alpha, &      ! Cartesian component used to construct pijSPIN
    nbki, &       ! number of available bands at the current k-point and spin
    nbcderki      ! number of target bands at the current k-point and spin

! File-name parsing
INTEGER :: &
    ipath         ! position of the final directory separator

! Units of already-open WIEN2k matrix-element files
INTEGER :: &
    unitinp = -1, &
    unitinpUP = -1, &
    unitinpDN = -1

! Band counts and degenerate-group metadata
INTEGER, ALLOCATABLE :: &
    nbk(:,:), &         ! available bands at each k-point and spin
    dg_groupk(:,:,:), & ! degenerate-group label of every available band
    ngroupsk(:,:), &    ! number of degenerate groups at each k-point and spin
    nbcderk(:,:)         ! target-band limit after degeneracy adjustment

! Fermi-level and occupation metadata
REAL(KIND=sp) :: &
    efermi, &          ! Fermi energy [Ha]
    full_occupation    ! occupation assigned to a fully occupied band

! Eigenvalues, occupations, and k-point weights
REAL(KIND=sp), ALLOCATABLE :: &
    ene(:,:,:), &      ! eigenvalues (band, k-point, spin) [Ha]
    occup(:,:,:), &    ! occupations (band, k-point, spin)
    wk(:)              ! relative WIEN2k k-point weights

! Degeneracy and occupation thresholds
REAL(KIND=sp), PARAMETER :: &
    energy_tolerance = 1.0E-5_sp, &    ! adjacent-band tolerance [Ha]
    occupation_tolerance = 1.0E-5_sp   ! threshold for nonzero occupation

! Energy differences at the current k-point
REAL(KIND=sp), ALLOCATABLE :: &
    dEij(:,:)          ! E_transition-E_target (target, transition) [Ha]

! Complete VASP energy-difference data
REAL(KIND=sp), ALLOCATABLE :: &
    dEijks(:,:,:,:)    ! E_j-E_i (target, transition, k-point, spin) [Ha]

! Matrix elements at the current k-point
COMPLEX(KIND=sp), ALLOCATABLE :: &
    pij(:,:,:), &      ! total matrix elements (component, target, transition) [a.u.]
    pijUP(:,:,:), &    ! spin-up-projected matrix elements (WIEN2k only) [a.u.]
    pijDN(:,:,:), &    ! spin-down-projected matrix elements (WIEN2k only) [a.u.]
    pijSPIN(:,:,:)     ! pijUP-pijDN used for spin Berry curvature [a.u.]

! Complete VASP matrix-element data
COMPLEX(KIND=sp), ALLOCATABLE :: &
    pijks(:,:,:,:,:)   ! (component, target, transition, k-point, spin) [a.u.]

! Band-resolved Berry curvature and orbital angular momentum
REAL(KIND=dp), ALLOCATABLE :: &
    bcurv(:,:), &      ! Omega_yz, Omega_zx, Omega_xy by target band [bohr^2]
    oam(:,:)           ! L_yz, L_zx, L_xy by target band [hbar]

! Generalized orbital-angular-momentum matrices
COMPLEX(KIND=dp), ALLOCATABLE :: &
    goam(:,:,:)        ! (bra band, ket band, yz/zx/xy) [hbar]

! Calculation mode
LOGICAL :: &
    wien2k, &            ! true for WIEN2k input and false for VASP input
    soc_sp_resolv_pij, & ! separate WIEN2k spin-up/down matrices are supplied
    dftocc, &             ! use occupations read from the DFT calculation
    soc                   ! WIEN2k calculation includes spin-orbit coupling

!! Get command line input arguments

CALL command_line_args( &
    fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc, dftocc, efermi) ! out
WRITE(*,'(A)') ' Command line arguments were initialized successfully.'

! Three matrix-element files indicate spin-resolved WIEN2k SOC input
soc_sp_resolv_pij = LEN_TRIM(fnameinpUP) > 0 .AND. LEN_TRIM(fnameinpDN) > 0
! Distinguish from SOC calculation without spin-resolved matrix elements
soc = soc_sp_resolv_pij .OR. filename_contains_energyso(fnameene)

!! WIEN2k or VASP?

! Extract the file name without the directory path
ipath = SCAN(TRIM(fnameinp), '/', BACK=.TRUE.)
IF (ipath > 0) THEN
    fnamebase = fnameinp(ipath+1:LEN_TRIM(fnameinp))
ELSE
    fnamebase = TRIM(fnameinp)
END IF
IF (fnamebase == 'WAVEDER' .OR. fnamebase == 'WAVEDERF') THEN ! VASP
    wien2k = .FALSE.
    IF (soc_sp_resolv_pij) THEN
        WRITE(*,*) 'Three matrix-element files can be used only with WIEN2k.'
        WRITE(*,*) 'Check execution options with "berrycpt -h"'
        ERROR STOP
    END IF
    ! Change to the binary file while preserving its directory path
    IF (fnamebase == 'WAVEDERF') THEN
        IF (ipath > 0) THEN
            fnameinp = fnameinp(:ipath)//'WAVEDER'
        ELSE
            fnameinp = 'WAVEDER'
        END IF
    END IF
    WRITE(*,*) 'Assumed VASP calculation based on the input file name.'
ELSE ! WIEN2k
    wien2k = .TRUE.
    WRITE(*,*) 'Assumed WIEN2k calculation based on the input file name.'
    WRITE(*,*) '(If you would like to read a VASP file, its file name should'
    WRITE(*,*) 'be WAVEDER or WAVEDERF.)'
    IF (soc_sp_resolv_pij) THEN
        WRITE(*,*) 'WIEN2k spin-resolved momentum matrix elements will be read'
        WRITE(*,'(A)') ' from the files: '//TRIM(fnameinpUP)//' and '// &
            TRIM(fnameinpDN)
    END IF
END IF

!! Validate input files

CALL validate_input_files( &
    fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc, & ! in
    wien2k, soc_sp_resolv_pij, dftocc)                                 ! in

!! Read input metadata

IF (.NOT. wien2k) THEN
    CALL read_waveder_header( &
        fnameinp, &                  ! in
        nstot, nktot, nbcder, nb)    ! out, scalar WAVEDER dimensions
END IF

!! Read eigenvalues from WIEN2k or VASP. Return them in [Ha]

IF (wien2k) THEN
    CALL read_eigenvalues_wien2k( &
        fnameene, &                  ! in
        nktot, nb, &                 ! out, scalar dimensions
        nbk, &                       ! out, allocated as nbk(nktot,1)
        ene, &                       ! out, allocated as ene(nb,nktot,1)
        wk)                          ! out, allocated as wk(nktot)

    nstot = 1
ELSE
    CALL read_eigenvalues_vasp( &
        fnameene, &                  ! in
        nstot, nktot, nb, &          ! in, dimensions read from WAVEDER
        nbk, &                       ! out, allocated as nbk(nktot,nstot)
        ene, &                       ! out, allocated as ene(nb,nktot,nstot)
        occup)                       ! out, allocated as occup(nb,nktot,nstot)
END IF

!! Read WIEN2k band occupations from `case.weight` only if requested
!! Should be generated by:
!!   x lapw2 [-so]

IF (wien2k .AND. dftocc) THEN
    CALL read_occupations_wien2k( &
        fnameocc, nktot, nb, &       ! in
        nbk, ene, wk, soc, &         ! in
        occup)                       ! out, allocated as occup(nb,nktot,1)
END IF

!! Construct occupation if Fermi energy is provided

IF ( .NOT. dftocc ) THEN
    IF (IEEE_IS_NAN(efermi)) THEN
        WRITE(*,'(/,A)') ' ERROR: Fermi energy was not initialized'
        ERROR STOP 1
    END IF
    CALL construct_occupations( &
        wien2k, soc, fnameinp, fnameene, &    ! in
        nstot, nktot, nb, nbk, ene, efermi, & ! in
        occup, &                              ! in/out, allocated or resized for
                                              ! WIEN2k as occup(nb,nktot,nstot);
                                              ! existing VASP array is overwritten
        full_occupation)                      ! out, full-band occupation
END IF

!! Find groups of degenerate bands

CALL find_degenerate_groups( &
    wien2k, nstot, nktot, nb, nbcder, nbk, ene, occup, & ! in
    energy_tolerance, occupation_tolerance, &            ! in
    dg_groupk, &  ! out, allocated as dg_groupk(nb,nktot,nstot)
                  !      = degenerate-group label (band,k-point,spin)
    ngroupsk, &   ! out, allocated as ngroupsk(nktot,nstot)
                  !      = number of groups at each k-point and spin
    nbcderk)      ! out, allocated as nbcderk(nktot,nstot)
                  !      = target-band limit after degeneracy adjustment

!! Read the complete VASP WAVEDER record once

IF (.NOT. wien2k) THEN
    CALL read_mommat_pij_vasp( &
        fnameinp, nstot, nktot, nbcder, nb, ene, & ! in
        pijks, &  ! out, allocated as
                ! pijks(3,nbcder,nb,nktot,nstot)
                ! = (Cartesian component, target band, transition band,
                !    k-point, spin)
        dEijks) ! out, allocated as
                ! dEijks(nbcder,nb,nktot,nstot)
                ! = E_transition-E_target [Ha]
END IF

!! Open WIEN2k matrix-element files

IF (wien2k) THEN
    CALL open_mommat_files( &
        fnameinp, fnameinpUP, fnameinpDN, & ! in
        soc_sp_resolv_pij, &                ! in
        unitinp, &                          ! out, connected unit for the total file
        unitinpUP, &                        ! out, connected spin-up unit or -1
        unitinpDN)                          ! out, connected spin-down unit or -1
END IF

!! Main processing loops

WRITE(*,'(A)') ' Entering the main spin and k-point loops.'
CALL initialize_progress( &
    nstot, nktot) ! in

DO ispin = 1, nstot

    CALL set_spin_suffix( &
        wien2k, nstot, ispin, fnameinp, & ! in
        charspin)                         ! out

    CALL open_output_files( &
        wien2k, soc_sp_resolv_pij, charspin, & ! in
        output_units)           ! out, scalar derived type;
                                ! inactive unit fields remain -1

    DO ikpt = 1, nktot

        nbki = nbk(ikpt,ispin)
        nbcderki = nbcderk(ikpt,ispin)

        !! Obtain rectangular matrices for the current k-point and spin

        IF (wien2k) THEN
            ! Read one k-point from the already-open sequential files.
            CALL read_mommat_pij_wien2k_kpoint( &
                unitinp, unitinpUP, unitinpDN, &           ! in
                soc_sp_resolv_pij, ikpt, nbki, nbcderki, & ! in
                pij, &  ! out, allocated as pij(3,nbcderki,nbki)
                        !      = (Cartesian component, target band, transition band)
                dEij, & ! out, allocated as dEij(nbcderki,nbki)
                        !      = E_transition-E_target [Ha]
                pijUP, &! out, allocated as pijUP(3,nbcderki,nbki)
                        !      only when soc_sp_resolv_pij is true
                pijDN)  ! out, allocated as pijDN(3,nbcderki,nbki)
                        !      only when soc_sp_resolv_pij is true
        ELSE
            ! Extract the current VASP matrix elements and energy differences.
            ALLOCATE( &
                pij(3,nbcderki,nbki), &
                dEij(nbcderki,nbki))

            pij = pijks(:,1:nbcderki,1:nbki,ikpt,ispin)
            dEij = dEijks(1:nbcderki,1:nbki,ikpt,ispin)
        END IF

        !! Ordinary Berry curvature

        CALL calculate_bcurv_kpoint( &
            nbki, nbcderki, &                       ! in
            dg_groupk(1:nbki,ikpt,ispin), &         ! in
            pij, pij, dEij, &                       ! in
            bcurv)                                  ! out, allocated as
                                                    ! bcurv(nbcderki,3)
                                                    ! = (target band, yz/zx/xy)
                                                    ! in bohr^2

        CALL write_bcurv_kpoint( &
            output_units%bcurv, ikpt, nbki, nbcderki, & ! in
            dg_groupk(1:nbki,ikpt,ispin), &             ! in
            occup(1:nbcderki,ikpt,ispin), bcurv)        ! in

        DEALLOCATE(bcurv)

        !! Spin-resolved Berry-curvature quantities

        IF (wien2k .AND. soc_sp_resolv_pij) THEN

            !! Spin-up contribution to the Berry curvature
            !
            ! First operator:  pijA = pijUP
            ! Second operator: pijB = pij

            CALL calculate_bcurv_kpoint( &
                nbki, nbcderki, &                         ! in
                dg_groupk(1:nbki,ikpt,ispin), &           ! in
                pijUP, pij, dEij, &                       ! in
                bcurv)                                    ! out

            CALL write_bcurv_kpoint( &
                output_units%bcurv_up, ikpt, nbki, nbcderki, & ! in
                dg_groupk(1:nbki,ikpt,ispin), &                ! in
                occup(1:nbcderki,ikpt,ispin), bcurv)           ! in

            DEALLOCATE(bcurv)

            !! Spin-down contribution to the Berry curvature
            !
            ! First operator:  pijA = pijDN
            ! Second operator: pijB = pij

            CALL calculate_bcurv_kpoint( &
                nbki, nbcderki, &                         ! in
                dg_groupk(1:nbki,ikpt,ispin), &           ! in
                pijDN, pij, dEij, &                       ! in
                bcurv)                                    ! out

            CALL write_bcurv_kpoint( &
                output_units%bcurv_dn, ikpt, nbki, nbcderki, & ! in
                dg_groupk(1:nbki,ikpt,ispin), &                ! in
                occup(1:nbcderki,ikpt,ispin), bcurv)           ! in

            DEALLOCATE(bcurv)

            !! Sigma_z-normalized spin Berry curvature
            !
            ! First operator:  pijA = pijUP - pijDN
            ! Second operator: pijB = pij
            !
            ! Construct the first-operator matrix explicitly. Scalar loops are
            ! used so that the subtraction does not require an array-expression
            ! temporary.

            ALLOCATE(pijSPIN(3,nbcderki,nbki))

            DO j = 1, nbki
                DO i = 1, nbcderki
                    DO alpha = 1, 3
                        pijSPIN(alpha,i,j) = &
                            pijUP(alpha,i,j) - pijDN(alpha,i,j)
                    END DO
                END DO
            END DO

            CALL calculate_bcurv_kpoint( &
                nbki, nbcderki, &                         ! in
                dg_groupk(1:nbki,ikpt,ispin), &           ! in
                pijSPIN, pij, dEij, &                     ! in
                bcurv)                                    ! out

            CALL write_bcurv_kpoint( &
                output_units%bcurv_spin, ikpt, nbki, nbcderki, & ! in
                dg_groupk(1:nbki,ikpt,ispin), &                  ! in
                occup(1:nbcderki,ikpt,ispin), bcurv)             ! in

            DEALLOCATE(bcurv)
            DEALLOCATE(pijSPIN)

        END IF

        !! Ordinary orbital angular momentum

        ! The ordinary OAM uses the total momentum matrix elements in both
        ! operator positions:
        !
        !     pijA = pij
        !     pijB = pij

        CALL calculate_oam_kpoint( &
            nbki, nbcderki, &                       ! in
            dg_groupk(1:nbki,ikpt,ispin), &         ! in
            pij, pij, dEij, &                       ! in
            oam)                                    ! out, allocated as
                                                    ! oam(nbcderki,3)
                                                    ! = (target band, yz/zx/xy)
                                                    ! in hbar

        CALL write_oam_kpoint( &
            output_units%oam, ikpt, nbki, nbcderki, & ! in
            dg_groupk(1:nbki,ikpt,ispin), oam)        ! in

        DEALLOCATE(oam)

        !! Spin orbital angular momentum

        IF (wien2k .AND. soc_sp_resolv_pij) THEN

            ! The spin OAM implemented here is the spin-up-projected OAM.
            ! Unlike the spin-related Berry-curvature quantities, it uses the
            ! spin-up momentum matrix elements in both operator positions:
            !
            !     pijA = pijUP
            !     pijB = pijUP
            !
            ! It is therefore not calculated using pijUP-pijDN, nor does it
            ! use the total momentum matrix pij in the second operator
            ! position.

            CALL calculate_oam_kpoint( &
                nbki, nbcderki, &                         ! in
                dg_groupk(1:nbki,ikpt,ispin), &           ! in
                pijUP, pijUP, dEij, &                     ! in
                oam)                                      ! out

            CALL write_oam_kpoint( &
                output_units%oam_up, &                    ! in
                ikpt, nbki, nbcderki, &                   ! in
                dg_groupk(1:nbki,ikpt,ispin), oam)        ! in

            DEALLOCATE(oam)

        END IF

        !! Generalized orbital angular momentum

        CALL calculate_goam_kpoint( &
            nbki, nbcderki, &                   ! in
            dg_groupk(1:nbki,ikpt,ispin), &     ! in
            pij, pij, dEij, &                   ! in
            goam)                               ! out, allocated as
                                                ! goam(nbcderki,nbcderki,3)
                                                ! = <u_i|L_c|u_j>, with
                                                !   i = bra band
                                                !   j = ket band
                                                !   c = yz, zx, xy
                                                ! units: hbar

        CALL write_goam_kpoint( &
            output_units%goam, ikpt, nbki, nbcderki, & ! in
            dg_groupk(1:nbcderki,ikpt,ispin), goam)    ! in

        DEALLOCATE(goam)

        !! Release current k-point matrices

        IF (ALLOCATED(pij)) DEALLOCATE(pij)
        IF (ALLOCATED(dEij)) DEALLOCATE(dEij)
        IF (ALLOCATED(pijUP)) DEALLOCATE(pijUP)
        IF (ALLOCATED(pijDN)) DEALLOCATE(pijDN)

        CALL write_progress( &
            ispin, nstot, ikpt, nktot) ! in

    END DO

    CALL close_output_files( &
        output_units) ! in/out

END DO

!! Close WIEN2k matrix-element files

IF (wien2k) THEN
    CALL close_mommat_files( &
        unitinp, unitinpUP, unitinpDN) ! in/out
END IF

!! Release complete VASP matrices

IF (ALLOCATED(pijks)) DEALLOCATE(pijks)
IF (ALLOCATED(dEijks)) DEALLOCATE(dEijks)

WRITE(*,'(A)') ' Completed the main spin and k-point loops.'

CALL write_output_summary()

WRITE(*,'(/,A)') &
    ' Please cite BerryCPT in publications and other scholarly work'
WRITE(*,'(A)') &
    ' that use this software or results generated with it.'
WRITE(*,'(A)') &
    ' Suggested citation:'
    WRITE(*,'(A)') &
    ' Rubel, O. (2024). BerryCPT: Berry curvature and orbital angular'
WRITE(*,'(A)') &
    ' momentum from DFT calculations [Computer software].'
WRITE(*,'(A)') &
    ' https://github.com/rubel75/berrycpt'

END PROGRAM berrycpt
