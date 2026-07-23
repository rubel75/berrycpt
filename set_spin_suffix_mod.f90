MODULE set_spin_suffix_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: set_spin_suffix

CONTAINS

SUBROUTINE set_spin_suffix( &
        wien2k, nstot, ispin, fnameinp, &
        charspin)

! Calculation type
LOGICAL, INTENT(IN) :: wien2k

! Total number of spin channels and current spin channel
INTEGER, INTENT(IN) :: nstot, ispin

! Name of the matrix-element input file
CHARACTER(LEN=*), INTENT(IN) :: fnameinp

! Spin suffix appended to output filenames
CHARACTER(LEN=*), INTENT(OUT) :: charspin

! Length of the matrix-element input filename
INTEGER :: filename_length

IF (LEN(charspin) < 3) THEN
    WRITE(*,'(/,A)') ' ERROR: The character variable charspin is too short'
    WRITE(*,'(A,I0)') ' Its length is ', LEN(charspin)
    WRITE(*,'(A)') ' A length of at least 3 characters is required.'
    ERROR STOP 1
END IF

IF (nstot < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: The number of spin channels must be positive'
    WRITE(*,'(A,I0)') ' nstot = ', nstot
    ERROR STOP 1
END IF

IF (ispin < 1 .OR. ispin > nstot) THEN
    WRITE(*,'(/,A)') ' ERROR: The spin-channel index is outside its valid range'
    WRITE(*,'(A,I0)') ' ispin = ', ispin
    WRITE(*,'(A,I0)') ' nstot = ', nstot
    ERROR STOP 1
END IF

charspin = ''

IF (wien2k) THEN
    filename_length = LEN_TRIM(fnameinp)

    IF (filename_length >= 2) THEN
        SELECT CASE (fnameinp(filename_length-1:filename_length))
        CASE ('up')
            charspin = '-up'
        CASE ('dn')
            charspin = '-dn'
        END SELECT
    END IF
ELSE
    IF (nstot == 2) THEN
        SELECT CASE (ispin)
        CASE (1)
            charspin = '-up'
        CASE (2)
            charspin = '-dn'
        END SELECT
    END IF
END IF

END SUBROUTINE set_spin_suffix

END MODULE set_spin_suffix_mod
