MODULE read_occupations_wien2k_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_occupations_wien2k


CONTAINS

SUBROUTINE read_occupations_wien2k( &
    fnameocc, nktot, nb, nbk, ene, wk, soc, &
    occup)

! Name of the WIEN2k case.weight file
CHARACTER(LEN=*), INTENT(IN) :: fnameocc

! Eigenvalue-array dimensions
INTEGER, INTENT(IN) :: nktot, nb

! Number of eigenvalues at each k-point and spin
INTEGER, INTENT(IN) :: nbk(:,:)

! Eigenvalues [Ha] and relative k-point weights
REAL(KIND=sp), INTENT(IN) :: &
    ene(:,:,:), &
    wk(:)

! Whether the WIEN2k calculation includes spin-orbit coupling
LOGICAL, INTENT(IN) :: soc

! Band-, k-point-, and spin-dependent occupations
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: &
    occup(:,:,:)

! Current record and integration method
CHARACTER(LEN=1024) :: cline
CHARACTER(LEN=16) :: integration_method, ctext

! Input unit, file metadata, counters, and I/O status
INTEGER :: &
    unitocc, &
    nktot_file, nb_file, &
    ikpt, ikpt_file, ikpt_record, &
    iband, iband_record, &
    ierr

! Energies read from case.weight and the corresponding reference energy
REAL(KIND=sp) :: &
    efermi_ry, energy_min_ry, &
    energy_ry, energy_reference_ry

! k-point-weight normalization and occupation conversion
REAL(KIND=sp) :: &
    band_weight, &
    wk_sum, wk_normalized, &
    weight_correction, &
    occupation, occupation_max

! Status of the k-point header
LOGICAL :: kpoint_header

! I/O error message
CHARACTER(LEN=256) :: iomsg

! Numerical tolerances
REAL(KIND=sp), PARAMETER :: &
    energy_tolerance = 1.0E-5_sp, &
    occupation_tolerance = 1.0E-5_sp

IF (nktot <= 0 .OR. nb <= 0) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Invalid eigenvalue-array dimensions were supplied') ! in
END IF

IF (SIZE(nbk,1) /= nktot .OR. SIZE(nbk,2) < 1) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'The nbk array is inconsistent with nktot') ! in
END IF

IF (SIZE(ene,1) /= nb .OR. SIZE(ene,2) /= nktot .OR. &
    SIZE(ene,3) < 1) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'The eigenvalue array has inconsistent dimensions') ! in
END IF

IF (SIZE(wk) /= nktot) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'The k-point-weight array is inconsistent with nktot') ! in
END IF

wk_sum = SUM(wk)

IF (wk_sum <= 0.0_sp .OR. ANY(wk <= 0.0_sp)) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Invalid WIEN2k k-point weights were supplied') ! in
END IF

! Without SOC, case.weight may contain occupations up to two.
! For SOC, WIEN2k writes an additional factor of two in the
! integration weights. Remove this factor so that a fully
! occupied spinor band has occupation one.
weight_correction = 1.0_sp
occupation_max = 2.0_sp

IF (soc) THEN
    weight_correction = 2.0_sp
    occupation_max = 1.0_sp
END IF

CALL open_occupation_file( &
    fnameocc, & ! in
    unitocc)    ! out

! The first record contains the Fermi energy, lower energy
! limit, and integration method
READ(unitocc,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

IF (ierr /= 0) THEN
    CLOSE(unitocc)
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to read the first case.weight record: '// &
        TRIM(iomsg)) ! in
END IF

READ(cline,*,IOSTAT=ierr) &
    efermi_ry, energy_min_ry, integration_method

IF (ierr /= 0) THEN
    CLOSE(unitocc)
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to interpret the first case.weight record') ! in
END IF

! The second record contains the number of k-points
READ(unitocc,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

IF (ierr /= 0) THEN
    CLOSE(unitocc)
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to read the case.weight k-point count: '// &
        TRIM(iomsg)) ! in
END IF

READ(cline,*,IOSTAT=ierr) nktot_file, ctext

IF (ierr /= 0 .OR. TRIM(ADJUSTL(ctext)) /= 'K-POINTS') THEN
    CLOSE(unitocc)
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to interpret the case.weight k-point count') ! in
END IF

IF (nktot_file /= nktot) THEN
    CLOSE(unitocc)
    WRITE(*,'(/,A)') ' ERROR: Inconsistent WIEN2k k-point counts'
    WRITE(*,'(A,I0,A,I0)') &
        ' case.energy = ', nktot, ', case.weight = ', nktot_file
    ERROR STOP 1
END IF

ALLOCATE(occup(nb,nktot,1))
occup = 0.0_sp

DO ikpt = 1, nktot
    CALL read_next_nonempty_line( &
        unitocc, fnameocc, & ! in
        cline)               ! out

    CALL parse_weight_kpoint_header( &
        cline, &                           ! in
        kpoint_header, ikpt_file, nb_file) ! out

    IF (.NOT. kpoint_header) THEN
        CLOSE(unitocc)
        WRITE(*,'(/,A)') ' ERROR: Invalid case.weight k-point header'
        WRITE(*,'(A)') ' '//TRIM(cline)
        ERROR STOP 1
    END IF

    IF (ikpt_file /= ikpt) THEN
        CLOSE(unitocc)
        WRITE(*,'(/,A)') ' ERROR: Inconsistent case.weight k-point index'
        WRITE(*,'(A,I0,A,I0)') &
            ' Expected k-point ', ikpt, ', found ', ikpt_file
        ERROR STOP 1
    END IF

    IF (nb_file /= nbk(ikpt,1)) THEN
        CLOSE(unitocc)
        WRITE(*,'(/,A)') ' ERROR: Inconsistent number of WIEN2k bands'
        WRITE(*,'(A,I0,A,I0,A,I0)') &
            ' K-point ', ikpt, ': case.energy = ', nbk(ikpt,1), &
            ', case.weight = ', nb_file
        ERROR STOP 1
    END IF

    wk_normalized = wk(ikpt) / wk_sum

    IF (wk_normalized <= 0.0_sp) THEN
        CLOSE(unitocc)
        WRITE(*,'(/,A)') ' ERROR: Invalid normalized WIEN2k k-point weight'
        WRITE(*,'(A,I0)') ' K-point ', ikpt
        WRITE(*,'(A,ES16.8)') &
            ' Normalized k-point weight = ', wk_normalized
        ERROR STOP 1
    END IF

    DO iband = 1, nb_file
        READ(unitocc,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

        IF (ierr /= 0) THEN
            CLOSE(unitocc)
            CALL occupation_file_error( &
                fnameocc, & ! in
                'Incomplete case.weight band block: '// &
                TRIM(iomsg)) ! in
        END IF

        READ(cline,*,IOSTAT=ierr) &
            ikpt_record, iband_record, energy_ry, band_weight

        IF (ierr /= 0) THEN
            CLOSE(unitocc)
            WRITE(*,'(/,A)') ' ERROR: Invalid case.weight band record'
            WRITE(*,'(A)') ' '//TRIM(cline)
            ERROR STOP 1
        END IF

        IF (ikpt_record /= ikpt .OR. iband_record /= iband) THEN
            CLOSE(unitocc)
            WRITE(*,'(/,A)') ' ERROR: Inconsistent case.weight band indices'
            WRITE(*,'(A,I0,A,I0)') &
                ' Expected k-point ', ikpt, ', band ', iband
            WRITE(*,'(A,I0,A,I0)') &
                ' Found k-point ', ikpt_record, ', band ', iband_record
            ERROR STOP 1
        END IF

        energy_reference_ry = 2.0_sp * ene(iband,ikpt,1)

        IF (ABS(energy_ry - energy_reference_ry) > &
            energy_tolerance * MAX(1.0_sp, ABS(energy_ry), &
            ABS(energy_reference_ry))) THEN

            CLOSE(unitocc)
            WRITE(*,'(/,A)') &
                ' ERROR: Inconsistent case.energy and case.weight eigenvalues'
            WRITE(*,'(A,I0,A,I0)') &
                ' K-point ', ikpt, ', band ', iband
            WRITE(*,'(A,ES16.8)') &
                ' case.energy eigenvalue [Ry] = ', energy_reference_ry
            WRITE(*,'(A,ES16.8)') &
                ' case.weight eigenvalue [Ry] = ', energy_ry
            ERROR STOP 1
        END IF

        occupation = band_weight / &
            (weight_correction * wk_normalized)

        ! Remove small numerical deviations from exact occupations
        IF (ABS(occupation) <= occupation_tolerance) THEN
            occupation = 0.0_sp
        ELSE IF (ABS(occupation - 1.0_sp) <= &
            occupation_tolerance) THEN

            occupation = 1.0_sp
        ELSE IF (.NOT. soc .AND. ABS(occupation - 2.0_sp) <= &
            occupation_tolerance) THEN

            occupation = 2.0_sp
        END IF

        IF (occupation < 0.0_sp .OR. &
            occupation > occupation_max) THEN

            CLOSE(unitocc)
            WRITE(*,'(/,A)') ' ERROR: Invalid WIEN2k band occupation'
            WRITE(*,'(A,I0,A,I0)') &
                ' K-point ', ikpt, ', band ', iband
            WRITE(*,'(A,ES16.8)') &
                ' Normalized occupation = ', occupation
            WRITE(*,'(A,ES16.8)') &
                ' Maximum allowed occupation = ', occupation_max
            ERROR STOP 1
        END IF

        occup(iband,ikpt,1) = occupation
    END DO
END DO

CLOSE(unitocc,IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to close the occupation file: '//TRIM(iomsg)) ! in
END IF

WRITE(*,'(A)') ' WIEN2k occupations were read successfully.'
WRITE(*,'(A,A)') &
    ' Integration method = ', TRIM(integration_method)

END SUBROUTINE read_occupations_wien2k


SUBROUTINE parse_weight_kpoint_header( &
    cline, &
    kpoint_header, ikpt, nb)

! Current case.weight record
CHARACTER(LEN=*), INTENT(IN) :: cline

! Parsed k-point metadata
LOGICAL, INTENT(OUT) :: kpoint_header
INTEGER, INTENT(OUT) :: ikpt, nb

! Record positions and internal read status
INTEGER :: &
    ipos_kpoint, ipos_colon, ipos_bands, &
    ierr

kpoint_header = .FALSE.
ikpt = 0
nb = 0

ipos_kpoint = INDEX(cline,'K-PNT')
ipos_colon = INDEX(cline,':')
ipos_bands = INDEX(cline,'BANDS')

IF (ipos_kpoint <= 0 .OR. ipos_colon <= ipos_kpoint + 5 .OR. &
    ipos_bands <= ipos_colon + 1) RETURN

READ(cline(ipos_kpoint+5:ipos_colon-1),*,IOSTAT=ierr) ikpt
IF (ierr /= 0) RETURN

READ(cline(ipos_colon+1:ipos_bands-1),*,IOSTAT=ierr) nb
IF (ierr /= 0) RETURN

IF (ikpt <= 0 .OR. nb <= 0) RETURN

kpoint_header = .TRUE.

END SUBROUTINE parse_weight_kpoint_header


SUBROUTINE read_next_nonempty_line( &
    unitocc, fnameocc, &
    cline)

! Input unit and file name
INTEGER, INTENT(IN) :: unitocc
CHARACTER(LEN=*), INTENT(IN) :: fnameocc

! First nonempty record
CHARACTER(LEN=*), INTENT(OUT) :: cline

! I/O status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

DO
    READ(unitocc,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr /= 0) THEN
        CALL occupation_file_error( &
            fnameocc, & ! in
            'Unable to read the next case.weight record: '// &
            TRIM(iomsg)) ! in
    END IF

    IF (LEN_TRIM(cline) > 0) EXIT
END DO

END SUBROUTINE read_next_nonempty_line


SUBROUTINE open_occupation_file( &
    fnameocc, &
    unitocc)

! Name of the case.weight file
CHARACTER(LEN=*), INTENT(IN) :: fnameocc

! Input unit
INTEGER, INTENT(OUT) :: unitocc

! I/O status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

OPEN(NEWUNIT=unitocc,FILE=TRIM(fnameocc),FORM='FORMATTED', &
    STATUS='OLD',ACTION='READ',IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL occupation_file_error( &
        fnameocc, & ! in
        'Unable to open the occupation file: '//TRIM(iomsg)) ! in
END IF

END SUBROUTINE open_occupation_file


SUBROUTINE occupation_file_error( &
    fnameocc, message)

! Name of the occupation file and description of the error
CHARACTER(LEN=*), INTENT(IN) :: fnameocc, message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
WRITE(*,'(A)') ' Occupation file: '//TRIM(fnameocc)
ERROR STOP 1

END SUBROUTINE occupation_file_error

END MODULE read_occupations_wien2k_mod