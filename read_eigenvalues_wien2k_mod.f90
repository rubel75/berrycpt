MODULE read_eigenvalues_wien2k_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_eigenvalues_wien2k


CONTAINS

SUBROUTINE read_eigenvalues_wien2k( &
        fnameene, &
        nktot, nb, nbk, ene, wk)

! Name of the WIEN2k case.energy file
CHARACTER(LEN=*), INTENT(IN) :: fnameene

! Number of k-points and maximum number of eigenvalues
INTEGER, INTENT(OUT) :: nktot, nb

! Number of eigenvalues at each k-point and spin
INTEGER, ALLOCATABLE, INTENT(OUT) :: &
    nbk(:,:)

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: &
    ene(:,:,:)

! Relative k-point weights from the WIEN2k energy file
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: &
    wk(:)

! Text from the current record and optional k-point label
CHARACTER(LEN=1024) :: cline
CHARACTER(LEN=10) :: kname

! Input unit, band counts, counters, and I/O status
INTEGER :: &
    unitene, &
    nb_current, &
    ikpt, iband, iband_file, ierr

! K-point metadata and eigenvalue in Ry
REAL(KIND=sp) :: &
    kvec(3), &
    kweight, &
    energy_ry

! Status of the current record
LOGICAL :: kpoint_header

! I/O error message
CHARACTER(LEN=256) :: iomsg

CALL open_energy_file( &
    fnameene, & ! in
    unitene)    ! out

! First pass: determine the number of k-points and the maximum
! number of eigenvalues at any k-point
nktot = 0
nb = 0

DO
    READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr < 0) EXIT

    IF (ierr > 0) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'Unable to scan the WIEN2k energy file: '//TRIM(iomsg)) ! in
    END IF

    CALL parse_wien2k_kpoint_header( &
        cline, &                             ! in
        kpoint_header, kvec, kname, &        ! out
        nb_current, kweight)                 ! out

    IF (.NOT. kpoint_header) CYCLE

    nktot = nktot + 1
    nb = MAX(nb, nb_current)

    DO iband = 1, nb_current
        READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

        IF (ierr /= 0) THEN
            CLOSE(unitene)
            CALL eigenvalue_file_error( &
                fnameene, & ! in
                'Incomplete WIEN2k eigenvalue block: '//TRIM(iomsg)) ! in
        END IF

        READ(cline,*,IOSTAT=ierr) iband_file, energy_ry

        IF (ierr /= 0 .OR. iband_file /= iband) THEN
            CLOSE(unitene)
            WRITE(*,'(/,A)') ' ERROR: Invalid WIEN2k eigenvalue record'
            WRITE(*,'(A)') ' '//TRIM(cline)
            ERROR STOP 1
        END IF
    END DO
END DO

IF (nktot <= 0 .OR. nb <= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'No WIEN2k k-point records were found') ! in
END IF

REWIND(unitene,IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to rewind the energy file: '//TRIM(iomsg)) ! in
END IF

ALLOCATE(nbk(nktot,1))
ALLOCATE(ene(nb,nktot,1))
ALLOCATE(wk(nktot))

nbk = 0
ene = HUGE(0.0_sp)
wk = 0.0_sp
ikpt = 0

! Second pass: read eigenvalues and k-point weights
DO
    READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

    IF (ierr < 0) EXIT

    IF (ierr > 0) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'Unable to read the WIEN2k energy file: '//TRIM(iomsg)) ! in
    END IF

    CALL parse_wien2k_kpoint_header( &
        cline, &                             ! in
        kpoint_header, kvec, kname, &        ! out
        nb_current, kweight)                 ! out

    IF (.NOT. kpoint_header) CYCLE

    ikpt = ikpt + 1

    IF (ikpt > nktot .OR. nb_current > nb) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'WIEN2k metadata changed between read passes') ! in
    END IF

    IF (kweight <= 0.0_sp) THEN
        CLOSE(unitene)
        CALL eigenvalue_file_error( &
            fnameene, & ! in
            'A non-positive WIEN2k k-point weight was found') ! in
    END IF

    nbk(ikpt,1) = nb_current
    wk(ikpt) = kweight

    DO iband = 1, nb_current
        READ(unitene,'(A)',IOSTAT=ierr,IOMSG=iomsg) cline

        IF (ierr /= 0) THEN
            CLOSE(unitene)
            CALL eigenvalue_file_error( &
                fnameene, & ! in
                'Incomplete WIEN2k eigenvalue block: '//TRIM(iomsg)) ! in
        END IF

        READ(cline,*,IOSTAT=ierr) iband_file, energy_ry

        IF (ierr /= 0 .OR. iband_file /= iband) THEN
            CLOSE(unitene)
            WRITE(*,'(/,A)') ' ERROR: Invalid WIEN2k eigenvalue record'
            WRITE(*,'(A)') ' '//TRIM(cline)
            ERROR STOP 1
        END IF

        ! Convert the eigenvalue from Ry to Ha
        ene(iband,ikpt,1) = energy_ry / 2.0_sp
    END DO
END DO

IF (ikpt /= nktot) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'The number of WIEN2k k-points changed between read passes') ! in
END IF

IF (ANY(nbk(:,1) <= 0) .OR. MAXVAL(nbk(:,1)) /= nb) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Invalid WIEN2k band-count metadata were read') ! in
END IF

IF (ANY(wk <= 0.0_sp)) THEN
    CLOSE(unitene)
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Invalid WIEN2k k-point weights were read') ! in
END IF

CLOSE(unitene,IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to close the energy file: '//TRIM(iomsg)) ! in
END IF

WRITE(*,'(A)') ' WIEN2k eigenvalues were read successfully.'
WRITE(*,'(A,I0)') ' Number of k-points = ', nktot
WRITE(*,'(A,I0)') ' Maximum number of eigenvalues = ', nb

END SUBROUTINE read_eigenvalues_wien2k


SUBROUTINE parse_wien2k_kpoint_header( &
        cline, &
        kpoint_header, kvec, kname, nb, kweight)

! Current WIEN2k energy-file record
CHARACTER(LEN=*), INTENT(IN) :: cline

! Parsed k-point metadata
LOGICAL, INTENT(OUT) :: kpoint_header
REAL(KIND=sp), INTENT(OUT) :: kvec(3), kweight
CHARACTER(LEN=10), INTENT(OUT) :: kname
INTEGER, INTENT(OUT) :: nb

! Unused matrix-size field and internal read status
INTEGER :: idummy, ierr

kpoint_header = .FALSE.
kvec = 0.0_sp
kname = ''
nb = 0
kweight = 0.0_sp

! WIEN2k case.energy k-point record:
! kx, ky, kz, optional label, unused matrix size,
! number of eigenvalues, weight
READ(cline,'(3E19.12,A10,2I6,F5.1)',IOSTAT=ierr) &
    kvec, kname, idummy, nb, kweight

IF (ierr /= 0) RETURN
IF (nb <= 0) RETURN
IF (kweight <= 0.0_sp) RETURN

kpoint_header = .TRUE.

END SUBROUTINE parse_wien2k_kpoint_header


SUBROUTINE open_energy_file( &
        fnameene, &
        unitene)

! Name of the energy file
CHARACTER(LEN=*), INTENT(IN) :: fnameene

! Input unit
INTEGER, INTENT(OUT) :: unitene

! I/O status and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

OPEN(NEWUNIT=unitene,FILE=TRIM(fnameene),FORM='FORMATTED', &
    STATUS='OLD',ACTION='READ',IOSTAT=ierr,IOMSG=iomsg)

IF (ierr /= 0) THEN
    CALL eigenvalue_file_error( &
        fnameene, & ! in
        'Unable to open the energy file: '//TRIM(iomsg)) ! in
END IF

END SUBROUTINE open_energy_file


SUBROUTINE eigenvalue_file_error( &
        fnameene, message)

! Name of the energy file and description of the error
CHARACTER(LEN=*), INTENT(IN) :: fnameene, message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
WRITE(*,'(A)') ' Energy file: '//TRIM(fnameene)
ERROR STOP 1

END SUBROUTINE eigenvalue_file_error

END MODULE read_eigenvalues_wien2k_mod