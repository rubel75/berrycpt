MODULE command_line_args_mod

USE precision_mod, ONLY: sp
USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_value, ieee_quiet_nan

IMPLICIT NONE

PRIVATE
PUBLIC :: command_line_args

CONTAINS

SUBROUTINE command_line_args( &
        fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc, dftocc, efermi)

! Names of the matrix-element, energy, and occupation files
CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: &
    fnameinp, fnameinpUP, fnameinpDN, fnameene, fnameocc

! Fermi energy in [Ha]
REAL(KIND=sp), INTENT(OUT) :: efermi

! Use DFT occupations instead of the step function
LOGICAL, INTENT(OUT) :: dftocc

! Command-line argument and the value following an option
CHARACTER(LEN=:), ALLOCATABLE :: argument, value

! Number and index of command-line arguments and input status
INTEGER :: narg, iarg, ierr

! Number of matrix-element files
INTEGER :: ninp

! Status of the command-line arguments
LOGICAL :: fnameene_set, fnameocc_set, efermi_set

fnameinp = ''
fnameinpUP = ''
fnameinpDN = ''
fnameene = ''
fnameocc = ''
efermi = IEEE_VALUE(efermi, ieee_quiet_nan)
dftocc = .FALSE.

ninp = 0
fnameene_set = .FALSE.
fnameocc_set = .FALSE.
efermi_set = .FALSE.

narg = COMMAND_ARGUMENT_COUNT()
WRITE(*,'(A,I0)') ' Detected input arguments = ', narg

IF (narg == 0) THEN
    CALL command_line_error('No input arguments were provided')
END IF

iarg = 1
DO WHILE (iarg <= narg)
    CALL get_argument( &
        iarg, &      ! in
        argument)    ! out
    SELECT CASE (TRIM(argument))
        CASE ('-h', '--h', '-help', '--help')
            CALL print_help()
            STOP
        CASE ('--enefile')
            IF (fnameene_set) THEN
                CALL command_line_error('The option "--enefile" was specified more than once')
            END IF
            iarg = iarg + 1
            CALL get_option_value( &
                iarg, narg, '--enefile', & ! in
                value)                     ! out
            CALL copy_argument( &
                value, &     ! in
                fnameene)    ! out
            fnameene_set = .TRUE.
        CASE ('--occfile')
            IF (fnameocc_set) THEN
                CALL command_line_error('The option "--occfile" was specified more than once')
            END IF
            iarg = iarg + 1
            CALL get_option_value( &
                iarg, narg, '--occfile', & ! in
                value)                     ! out
            CALL copy_argument( &
                value, &     ! in
                fnameocc)    ! out
            fnameocc_set = .TRUE.
        CASE ('--efermiev')
            IF (efermi_set) THEN
                CALL command_line_error('The Fermi energy was specified more than once')
            END IF
            iarg = iarg + 1
            CALL get_option_value( &
                iarg, narg, '--efermiev', & ! in
                value)                      ! out
            READ(value,*,IOSTAT=ierr) efermi
            IF (ierr /= 0) THEN
                CALL command_line_error( &
                    'The value following "--efermiev" must be a Fermi energy in eV')
            END IF
            efermi = efermi / 27.21139_sp ! [eV] -> [Ha]
            efermi_set = .TRUE.
        CASE ('--efermiry')
            IF (efermi_set) THEN
                CALL command_line_error('The Fermi energy was specified more than once')
            END IF
            iarg = iarg + 1
            CALL get_option_value( &
                iarg, narg, '--efermiry', & ! in
                value)                      ! out
            READ(value,*,IOSTAT=ierr) efermi ! in [Ry]
            IF (ierr /= 0) THEN
                CALL command_line_error( &
                    'The value following "--efermiry" must be a Fermi energy in Ry')
            END IF
            efermi = efermi / 2.0_sp ! [Ry] -> [Ha]
            efermi_set = .TRUE.
        CASE ('--dftocc')
            IF (dftocc) THEN
                CALL command_line_error('The option "--dftocc" was specified more than once')
            END IF
            dftocc = .TRUE.
        CASE DEFAULT
            IF (LEN_TRIM(argument) == 0) THEN
                CALL command_line_error('An empty command-line argument was detected')
            ELSE IF (argument(1:1) == '-') THEN
                CALL command_line_error('Unknown option: '//TRIM(argument))
            END IF
            ninp = ninp + 1
            SELECT CASE (ninp)
                CASE (1)
                    CALL copy_argument( &
                        argument, & ! in
                        fnameinp)   ! out
                CASE (2)
                    CALL copy_argument( &
                        argument, &  ! in
                        fnameinpUP)  ! out
                CASE (3)
                    CALL copy_argument( &
                        argument, &  ! in
                        fnameinpDN)  ! out
                CASE DEFAULT
                    CALL command_line_error( &
                        'Specify either one or three matrix-element files')
            END SELECT
    END SELECT
    iarg = iarg + 1
END DO

IF (ninp /= 1 .AND. ninp /= 3) THEN
    CALL command_line_error( &
        'Specify either one matrix-element file or three files: total, up, and down')
END IF

IF (efermi_set .EQV. dftocc) THEN
    CALL command_line_error( &
        'Specify exactly one occupation mode: "--efermiev X.YY", "--efermiry X.YY", or "--dftocc"')
END IF

IF (fnameocc_set .AND. .NOT. dftocc) THEN
    CALL command_line_error( &
        'The option "--occfile" can be used only together with "--dftocc"')
END IF

WRITE(*,'(A)') ' Input matrix-element file = '//TRIM(fnameinp)
IF (ninp == 3) THEN
    WRITE(*,'(A)') ' Input spin-up matrix-element file = '//TRIM(fnameinpUP)
    WRITE(*,'(A)') ' Input spin-down matrix-element file = '//TRIM(fnameinpDN)
END IF

IF (fnameene_set) THEN
    WRITE(*,'(A)') ' Input energy file = '//TRIM(fnameene)
END IF

IF (dftocc) THEN
    WRITE(*,'(A)') ' Occupation mode = DFT occupations'
    IF (fnameocc_set) THEN
        WRITE(*,'(A)') ' Input occupation file = '//TRIM(fnameocc)
    END IF
ELSE
    WRITE(*,'(A,F10.6)') ' Fermi energy (Ha) = ', efermi
END IF

END SUBROUTINE command_line_args

SUBROUTINE get_argument( &
        iarg, &
        argument)

! Position of the command-line argument
INTEGER, INTENT(IN) :: iarg

! Command-line argument
CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: argument

! Length of the argument and command-line status
INTEGER :: larg, istat

CALL GET_COMMAND_ARGUMENT( &
    iarg, &              ! in
    LENGTH=larg, &       ! out
    STATUS=istat)        ! out
IF (istat /= 0) THEN
    CALL command_line_error('Unable to determine the command-line argument length')
END IF

ALLOCATE(CHARACTER(LEN=MAX(1,larg)) :: argument)

CALL GET_COMMAND_ARGUMENT( &
    iarg, &              ! in
    VALUE=argument, &    ! out
    STATUS=istat)        ! out
IF (istat /= 0) THEN
    CALL command_line_error('Unable to read a command-line argument')
END IF

END SUBROUTINE get_argument

SUBROUTINE get_option_value( &
        iarg, narg, option, &
        value)

! Position and total number of command-line arguments
INTEGER, INTENT(IN) :: iarg, narg

! Command-line option
CHARACTER(LEN=*), INTENT(IN) :: option

! Value following the command-line option
CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: value

IF (iarg > narg) THEN
    CALL command_line_error('The option "'//TRIM(option)//'" requires a value')
END IF

CALL get_argument( &
    iarg, &  ! in
    value)   ! out

IF (LEN_TRIM(value) == 0) THEN
    CALL command_line_error('The option "'//TRIM(option)//'" requires a non-empty value')
END IF

END SUBROUTINE get_option_value

SUBROUTINE copy_argument( &
        source, & ! in
        destination) ! out

! Command-line argument
CHARACTER(LEN=*), INTENT(IN) :: source

! Destination variable for the command-line argument
CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: destination

destination = TRIM(source)

END SUBROUTINE copy_argument

SUBROUTINE command_line_error( &
        message)

! Description of the command-line error
CHARACTER(LEN=*), INTENT(IN) :: message

WRITE(*,'(/,A)') ' ERROR: '//TRIM(message)
CALL print_help()
ERROR STOP 1

END SUBROUTINE command_line_error

SUBROUTINE print_help()

    WRITE(*,'(/,A)') ' Usage VASP:'
    WRITE(*,'(A)') '   berrycpt WAVEDER --enefile EIGENVAL --efermiev EF' // &
                    ' # Step-function occupations at E_F'
    WRITE(*,'(A)') '   berrycpt WAVEDER --enefile EIGENVAL --dftocc' // &
                    ' # DFT occupations from EIGENVAL'

    WRITE(*,'(/,A)') ' Usage WIEN2k without spin polarization or SOC:'
    WRITE(*,'(A)') '   berrycpt case.mommat2 --enefile case.energy' // &
                    ' --efermiry EF # Step-function occupations at E_F'
    WRITE(*,'(A)') '   berrycpt case.mommat2 --enefile case.energy' // &
                    ' --dftocc --occfile case.weight # DFT occupations'

    WRITE(*,'(/,A)') ' Usage spin-polarized WIEN2k without SOC:'
    WRITE(*,'(A)') '   berrycpt case.mommat2up --enefile case.energyup' // &
                    ' --efermiry EF # Spin-up channel'
    WRITE(*,'(A)') '   berrycpt case.mommat2up --enefile case.energyup' // &
                    ' --dftocc --occfile case.weightup # Spin-up channel'
    WRITE(*,'(A)') '   Replace the "up" suffix by "dn" for the spin-down channel.'

    WRITE(*,'(/,A)') ' Usage WIEN2k with SOC and without spin polarization:'
    WRITE(*,'(A)') '   berrycpt case.mommat2 --enefile case.energyso' // &
                    ' --efermiry EF # Ordinary quantities'
    WRITE(*,'(A)') '   berrycpt case.mommat2 --enefile case.energyso' // &
                    ' --dftocc --occfile case.weight # Ordinary quantities'
    WRITE(*,'(A)') '   berrycpt case.mommat2 case.mommat2up case.mommat2dn' // &
                    ' --enefile case.energyso --efermiry EF'
    WRITE(*,'(A)') '   berrycpt case.mommat2 case.mommat2up case.mommat2dn' // &
                    ' --enefile case.energyso --dftocc --occfile case.weight'

    WRITE(*,'(/,A)') ' Matrix-element files:'
    WRITE(*,'(A)') '   Specify one matrix-element file or three files:' // &
                    ' total, spin-up, and spin-down.'
    WRITE(*,'(A)') '   Three files are supported only for WIEN2k calculations' // &
                    ' with SOC.'

    WRITE(*,'(/,A)') ' Options:'
    WRITE(*,'(A)') '   --enefile FILE   Band-energy file required for all calculations'
    WRITE(*,'(A)') '   --efermiev EF    Step-function occupations using EF in eV'
    WRITE(*,'(A)') '   --efermiry EF    Step-function occupations using EF in Ry'
    WRITE(*,'(A)') '   --dftocc         Use occupations calculated by the DFT code'
    WRITE(*,'(A)') '   --occfile FILE   WIEN2k occupation file used with --dftocc'
    WRITE(*,'(A)') '   -h, --help       Print this help message'

END SUBROUTINE print_help

END MODULE command_line_args_mod