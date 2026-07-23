MODULE open_mommat_files_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: open_mommat_files, close_mommat_files

CONTAINS

SUBROUTINE open_mommat_files( &
        fnameinp, fnameinpUP, fnameinpDN, soc_sp_resolv_pij, &
        unitinp, unitinpUP, unitinpDN)

! Names of the matrix-element files
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, fnameinpUP, fnameinpDN

! Read spin-resolved WIEN2k matrix-element files
LOGICAL, INTENT(IN) :: soc_sp_resolv_pij

! File units for the matrix-element files
INTEGER, INTENT(OUT) :: unitinp, unitinpUP, unitinpDN

unitinp = -1
unitinpUP = -1
unitinpDN = -1

CALL open_single_mommat_file( &
    fnameinp, & ! in
    unitinp)    ! out

IF (soc_sp_resolv_pij) THEN
    CALL open_single_mommat_file( &
        fnameinpUP, & ! in
        unitinpUP)    ! out

    CALL open_single_mommat_file( &
        fnameinpDN, & ! in
        unitinpDN)    ! out
END IF

WRITE(*,*) 'Input matrix-element files were opened successfully.'

END SUBROUTINE open_mommat_files

SUBROUTINE open_single_mommat_file( &
        fnameinp, &
        unitinp)

! Name of the matrix-element file
CHARACTER(LEN=*), INTENT(IN) :: fnameinp

! File unit for the matrix-element file
INTEGER, INTENT(OUT) :: unitinp

! Error code and message
INTEGER :: ierr
CHARACTER(LEN=256) :: iomsg

OPEN( &
    NEWUNIT=unitinp, &
    FILE=TRIM(fnameinp), &
    STATUS='old', &
    ACTION='read', &
    FORM='formatted', &
    IOSTAT=ierr, &
    IOMSG=iomsg)

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR while opening the matrix-element file:'
    WRITE(*,'(A)') ' '//TRIM(fnameinp)
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP
END IF

END SUBROUTINE open_single_mommat_file



SUBROUTINE close_mommat_files( &
    unitinp, unitinpUP, unitinpDN) ! in/out

! Units of matrix-element files opened by open_mommat_files
INTEGER, INTENT(INOUT) :: &
    unitinp, &
    unitinpUP, &
    unitinpDN

CALL close_single_mommat_file( &
    unitinp) ! in/out
CALL close_single_mommat_file( &
    unitinpUP) ! in/out
CALL close_single_mommat_file( &
    unitinpDN) ! in/out

WRITE(*,'(A)') ' Input matrix-element files were closed successfully.'

END SUBROUTINE close_mommat_files


SUBROUTINE close_single_mommat_file( &
    unitinp) ! in/out

! Connected input-file unit, reset after successful closing
INTEGER, INTENT(INOUT) :: unitinp

! Actual file name associated with the Fortran unit
CHARACTER(LEN=512) :: filename

! Input/output status information
INTEGER :: ierr
LOGICAL :: opened
CHARACTER(LEN=512) :: iomsg

IF (unitinp == -1) RETURN

filename = ''
INQUIRE( &
    UNIT=unitinp, &
    OPENED=opened, &
    NAME=filename, &
    IOSTAT=ierr, &
    IOMSG=iomsg)

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to inquire about a matrix-element file'
    WRITE(*,'(A,I0)') ' Unit: ', unitinp
    WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
    ERROR STOP 1
END IF

IF (.NOT. opened) THEN
    WRITE(*,'(/,A)') ' ERROR: Matrix-element file unit is not connected'
    WRITE(*,'(A,I0)') ' Unit: ', unitinp
    ERROR STOP 1
END IF

CLOSE( &
    UNIT=unitinp, &
    STATUS='KEEP', &
    IOSTAT=ierr, &
    IOMSG=iomsg)

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to close a matrix-element file'
    WRITE(*,'(A)') ' File: ' // TRIM(filename)
    WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
    ERROR STOP 1
END IF

unitinp = -1

END SUBROUTINE close_single_mommat_file

END MODULE open_mommat_files_mod