MODULE construct_occupations_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: construct_occupations


CONTAINS

SUBROUTINE construct_occupations( &
        wien2k, soc, fnameinp, fnameene, &
        nstot, nktot, nb, nbk, ene, efermi, &
        occup, full_occupation)

! Calculation type
LOGICAL, INTENT(IN) :: &
    wien2k, & ! Whether this is a WIEN2k calculation
    soc       ! Whether the WIEN2k calculation includes SOC

! Names of the matrix-element and eigenvalue files
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, fnameene

! Eigenvalue-array dimensions
INTEGER, INTENT(IN) :: nstot, nktot, nb

! Number of eigenvalues at each k-point and spin
INTEGER, INTENT(IN) :: nbk(:,:)

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=sp), INTENT(IN) :: ene(:,:,:)

! Fermi energy [Ha]
REAL(KIND=sp), INTENT(IN) :: efermi

! Band-, k-point-, and spin-dependent occupations
REAL(KIND=sp), ALLOCATABLE, INTENT(INOUT) :: occup(:,:,:)

! Occupation assigned to a fully occupied band
REAL(KIND=sp), INTENT(OUT) :: full_occupation

! Band, k-point, and spin counters
INTEGER :: iband, ikpt, ispin

CALL validate_dimensions( &
    nstot, nktot, nb, nbk, ene) ! in

IF (wien2k) THEN
    CALL determine_wien2k_full_occupation( &
        soc, fnameinp, fnameene, & ! in
        full_occupation)           ! out

    IF (ALLOCATED(occup)) THEN
        IF (ANY(SHAPE(occup) /= [nb,nktot,nstot])) THEN
            DEALLOCATE(occup)
            ALLOCATE(occup(nb,nktot,nstot))
        END IF
    ELSE
        ALLOCATE(occup(nb,nktot,nstot))
    END IF
ELSE
    ! VASP occupations have already been read from EIGENVAL. Preserve
    ! their spin multiplicity by inspecting the occupation associated
    ! with the lowest eigenvalue in each spin channel before replacing
    ! the occupation array with a step function.
    IF (.NOT. ALLOCATED(occup)) THEN
        CALL occupation_construction_error( &
            'VASP occupations must be read before they are constructed') ! in
    END IF

    IF (ANY(SHAPE(occup) /= [nb,nktot,nstot])) THEN
        CALL occupation_construction_error( &
            'The VASP occupation array has inconsistent dimensions') ! in
    END IF

    CALL determine_vasp_full_occupation( &
        nstot, nktot, nbk, ene, occup, & ! in
        full_occupation)                 ! out
END IF

occup = 0.0_sp

DO ispin = 1, nstot
    DO ikpt = 1, nktot
        DO iband = 1, nbk(ikpt,ispin)
            IF (ene(iband,ikpt,ispin) < efermi) THEN
                occup(iband,ikpt,ispin) = full_occupation
            END IF
        END DO
    END DO
END DO

WRITE(*,'(A)') &
    ' Band occupations were constructed from the Fermi energy successfully.'
WRITE(*,'(A,ES16.8)') ' Fermi energy [Ha] = ', efermi
WRITE(*,'(A,F4.1)') ' Full band occupation = ', full_occupation

END SUBROUTINE construct_occupations


SUBROUTINE determine_wien2k_full_occupation( &
        soc, fnameinp, fnameene, &
        full_occupation)

! Whether the WIEN2k calculation includes SOC
LOGICAL, INTENT(IN) :: soc

! Names of the WIEN2k matrix-element and energy files
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, fnameene

! Occupation of a fully occupied WIEN2k band
REAL(KIND=sp), INTENT(OUT) :: full_occupation

! Spin-channel suffixes inferred from the input basenames
CHARACTER(LEN=2) :: matrix_spin, energy_spin

! An SOC eigenstate has occupation one.
IF (soc) THEN
    full_occupation = 1.0_sp
    RETURN
END IF

matrix_spin = filename_spin_suffix(fnameinp)
energy_spin = filename_spin_suffix(fnameene)

IF (LEN_TRIM(matrix_spin) == 0 .AND. LEN_TRIM(energy_spin) == 0) THEN
    ! A scalar-relativistic non-spin-polarized WIEN2k band represents
    ! two spin states.
    full_occupation = 2.0_sp
ELSE IF (LEN_TRIM(matrix_spin) == 0 .OR. LEN_TRIM(energy_spin) == 0) THEN
    WRITE(*,'(/,A)') &
        ' ERROR: Inconsistent WIEN2k spin-channel filename suffixes'
    WRITE(*,'(A)') ' Matrix-element file: '//TRIM(fnameinp)
    WRITE(*,'(A)') ' Energy file: '//TRIM(fnameene)
    WRITE(*,'(A)') &
        ' Both basenames must end in "up" or "dn" for a separate spin channel.'
    ERROR STOP 1
ELSE IF (matrix_spin /= energy_spin) THEN
    WRITE(*,'(/,A)') &
        ' ERROR: Conflicting WIEN2k spin-channel filename suffixes'
    WRITE(*,'(A)') ' Matrix-element file: '//TRIM(fnameinp)
    WRITE(*,'(A)') ' Energy file: '//TRIM(fnameene)
    WRITE(*,'(A,A)') ' Matrix-element spin channel = ', matrix_spin
    WRITE(*,'(A,A)') ' Energy-file spin channel = ', energy_spin
    ERROR STOP 1
ELSE
    ! A separate non-SOC WIEN2k spin channel has occupation one.
    full_occupation = 1.0_sp
END IF

END SUBROUTINE determine_wien2k_full_occupation


PURE FUNCTION filename_spin_suffix(fname) RESULT(spin_suffix)

! Input filename
CHARACTER(LEN=*), INTENT(IN) :: fname

! Conventional WIEN2k spin-channel suffix
CHARACTER(LEN=2) :: spin_suffix

! Lowercase filename
CHARACTER(LEN=LEN(fname)) :: fname_lower

! Character positions in the filename
INTEGER :: i, basename_start, basename_end

spin_suffix = ''

IF (LEN_TRIM(fname) < 2) RETURN

fname_lower = lowercase(fname)
basename_start = 1
basename_end = LEN_TRIM(fname_lower)

! Locate the beginning of the basename.
DO i = basename_end, 1, -1
    IF (fname_lower(i:i) == '/' .OR. &
        IACHAR(fname_lower(i:i)) == 92) THEN

        basename_start = i + 1
        EXIT
    END IF
END DO

IF (basename_end - basename_start + 1 < 2) RETURN

SELECT CASE (fname_lower(basename_end-1:basename_end))
    CASE ('up', 'dn')
        spin_suffix = fname_lower(basename_end-1:basename_end)
END SELECT

END FUNCTION filename_spin_suffix


PURE FUNCTION lowercase(string) RESULT(string_lower)

! Input character string
CHARACTER(LEN=*), INTENT(IN) :: string

! Lowercase character string
CHARACTER(LEN=LEN(string)) :: string_lower

! Character position and ASCII code
INTEGER :: i, character_code

string_lower = string

DO i = 1, LEN(string)
    character_code = IACHAR(string(i:i))

    IF (character_code >= IACHAR('A') .AND. &
        character_code <= IACHAR('Z')) THEN

        string_lower(i:i) = ACHAR( &
            character_code + IACHAR('a') - IACHAR('A'))
    END IF
END DO

END FUNCTION lowercase


SUBROUTINE determine_vasp_full_occupation( &
        nstot, nktot, nbk, ene, occup, &
        full_occupation)

! Eigenvalue-array dimensions
INTEGER, INTENT(IN) :: nstot, nktot

! Number of eigenvalues at each k-point and spin
INTEGER, INTENT(IN) :: nbk(:,:)

! VASP eigenvalues [Ha] and occupations read from EIGENVAL
REAL(KIND=sp), INTENT(IN) :: &
    ene(:,:,:), &
    occup(:,:,:)

! Occupation of a fully occupied VASP band
REAL(KIND=sp), INTENT(OUT) :: full_occupation

! Band, k-point, and spin counters
INTEGER :: iband, ikpt, ispin

! Lowest eigenvalue and its occupation in the current spin channel
REAL(KIND=sp) :: &
    lowest_energy, &
    lowest_energy_occupation, &
    spin_full_occupation

! Numerical tolerance for identifying an exact spin multiplicity
REAL(KIND=sp), PARAMETER :: occupation_tolerance = 1.0E-5_sp

full_occupation = 0.0_sp

DO ispin = 1, nstot
    lowest_energy = HUGE(0.0_sp)
    lowest_energy_occupation = -HUGE(0.0_sp)

    DO ikpt = 1, nktot
        DO iband = 1, nbk(ikpt,ispin)
            IF (ene(iband,ikpt,ispin) < lowest_energy) THEN
                lowest_energy = ene(iband,ikpt,ispin)
                lowest_energy_occupation = occup(iband,ikpt,ispin)
            END IF
        END DO
    END DO

    IF (ABS(lowest_energy_occupation - 1.0_sp) <= &
        occupation_tolerance) THEN

        spin_full_occupation = 1.0_sp
    ELSE IF (ABS(lowest_energy_occupation - 2.0_sp) <= &
        occupation_tolerance) THEN

        spin_full_occupation = 2.0_sp
    ELSE
        WRITE(*,'(/,A)') &
            ' ERROR: Unable to determine the full VASP band occupation'
        WRITE(*,'(A,I0)') ' Spin channel = ', ispin
        WRITE(*,'(A,ES16.8)') ' Lowest eigenvalue [Ha] = ', lowest_energy
        WRITE(*,'(A,ES16.8)') &
            ' Occupation of the lowest eigenvalue = ', &
            lowest_energy_occupation
        ERROR STOP 1
    END IF

    IF (ispin == 1) THEN
        full_occupation = spin_full_occupation
    ELSE IF (ABS(spin_full_occupation - full_occupation) > &
        occupation_tolerance) THEN

        WRITE(*,'(/,A)') &
            ' ERROR: Inconsistent full occupations in VASP spin channels'
        WRITE(*,'(A,F4.1)') &
            ' Full occupation in spin channel 1 = ', full_occupation
        WRITE(*,'(A,I0,A,F4.1)') &
            ' Full occupation in spin channel ', ispin, ' = ', &
            spin_full_occupation
        ERROR STOP 1
    END IF
END DO

END SUBROUTINE determine_vasp_full_occupation


SUBROUTINE validate_dimensions( &
        nstot, nktot, nb, nbk, ene)

! Eigenvalue-array dimensions
INTEGER, INTENT(IN) :: nstot, nktot, nb

! Number of eigenvalues at each k-point and spin
INTEGER, INTENT(IN) :: nbk(:,:)

! Band-, k-point-, and spin-dependent eigenvalues
REAL(KIND=sp), INTENT(IN) :: ene(:,:,:)

IF (nstot <= 0 .OR. nktot <= 0 .OR. nb <= 0) THEN
    CALL occupation_construction_error( &
        'Invalid eigenvalue-array dimensions were supplied') ! in
END IF

IF (SIZE(nbk,1) /= nktot .OR. SIZE(nbk,2) /= nstot) THEN
    CALL occupation_construction_error( &
        'The nbk array has inconsistent dimensions') ! in
END IF

IF (ANY(nbk <= 0) .OR. ANY(nbk > nb)) THEN
    CALL occupation_construction_error( &
        'The nbk array contains an invalid number of bands') ! in
END IF

IF (SIZE(ene,1) /= nb .OR. SIZE(ene,2) /= nktot .OR. &
    SIZE(ene,3) /= nstot) THEN

    CALL occupation_construction_error( &
        'The eigenvalue array has inconsistent dimensions') ! in
END IF

END SUBROUTINE validate_dimensions


SUBROUTINE occupation_construction_error( &
        message)

! Description of the occupation-construction error
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
ERROR STOP 1

END SUBROUTINE occupation_construction_error

END MODULE construct_occupations_mod
