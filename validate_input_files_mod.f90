MODULE validate_input_files_mod

IMPLICIT NONE

PRIVATE
PUBLIC :: validate_input_files

CONTAINS

SUBROUTINE validate_input_files( &
        fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc, &
        wien2k, soc_sp_resolv_pij, dftocc)

! Names of the matrix-element, energy, and occupation files
CHARACTER(LEN=*), INTENT(IN) :: &
    fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc

! Input-file configuration
LOGICAL, INTENT(IN) :: wien2k, soc_sp_resolv_pij, dftocc

IF (.NOT. wien2k .AND. soc_sp_resolv_pij) THEN
    CALL input_file_error( &
        'Three matrix-element files can be used only with WIEN2k')
END IF

CALL require_input_file( &
    fnameinp, 'matrix-element file') ! in

IF (soc_sp_resolv_pij) THEN
    CALL require_input_file( &
        fnameinpUP, 'spin-up matrix-element file') ! in
    CALL require_input_file( &
        fnameinpDN, 'spin-down matrix-element file') ! in
END IF

! Eigenvalue file is required for all calculations
IF (LEN_TRIM(fnameene) == 0) THEN
    CALL input_file_error( &
        'All calculations require "--enefile"')
END IF

IF (wien2k) THEN
    IF (dftocc .AND. LEN_TRIM(fnameocc) == 0) THEN
        CALL input_file_error( &
            'WIEN2k calculations with DFT occupations require "--occfile"')
    END IF
ELSE
    IF (LEN_TRIM(fnameocc) > 0) THEN
        CALL input_file_error( &
            'VASP occupations are read from EIGENVAL; "--occfile" is not used')
    END IF
END IF

IF (LEN_TRIM(fnameene) > 0) THEN
    CALL require_input_file( &
        fnameene, 'energy file') ! in
END IF

IF (LEN_TRIM(fnameocc) > 0) THEN
    CALL require_input_file( &
        fnameocc, 'occupation file') ! in
END IF

WRITE(*,'(A)') ' Input files were validated successfully.'

END SUBROUTINE validate_input_files

SUBROUTINE require_input_file( &
        filename, description)

! Name and description of the input file
CHARACTER(LEN=*), INTENT(IN) :: filename, description

! Status of the input file
LOGICAL :: file_exists

IF (LEN_TRIM(filename) == 0) THEN
    CALL input_file_error( &
        'No '//TRIM(description)//' was specified')
END IF

INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
IF (.NOT. file_exists) THEN
    CALL input_file_error( &
        'The '//TRIM(description)//' '//TRIM(filename)//' does not exist')
END IF

WRITE(*,'(A)') ' The input file '//TRIM(filename)//' was found.'

END SUBROUTINE require_input_file

SUBROUTINE input_file_error( &
        message)

! Description of the input-file error
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
ERROR STOP 1

END SUBROUTINE input_file_error

END MODULE validate_input_files_mod
