MODULE filename_contains_energyso_mod

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: filename_contains_energyso


    CONTAINS

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


    LOGICAL FUNCTION filename_contains_energyso(fname) RESULT(found)

    ! Name of the WIEN2k energy file
    CHARACTER(LEN=*), INTENT(IN) :: fname

    ! Lowercase file name
    CHARACTER(LEN=LEN(fname)) :: fname_lower

    ! Character position and beginning of the basename
    INTEGER :: i, basename_start

    found = .FALSE.

    IF (LEN_TRIM(fname) == 0) RETURN

    fname_lower = lowercase(fname)
    basename_start = 1

    ! Locate the beginning of the basename
    DO i = LEN_TRIM(fname_lower), 1, -1
        IF (fname_lower(i:i) == '/' .OR. &
            IACHAR(fname_lower(i:i)) == 92) THEN

            basename_start = i + 1
            EXIT
        END IF
    END DO

    found = INDEX( &
        fname_lower(basename_start:LEN_TRIM(fname_lower)), &
        'energyso') > 0

    END FUNCTION filename_contains_energyso

    END MODULE filename_contains_energyso_mod