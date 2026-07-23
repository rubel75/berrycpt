MODULE write_progress_mod

USE precision_mod, ONLY: dp
USE, INTRINSIC :: iso_fortran_env, ONLY: int64, OUTPUT_UNIT

IMPLICIT NONE

PRIVATE
PUBLIC :: initialize_progress, write_progress

! Number of progress intervals requested over the complete calculation.
! A value of 20 gives updates approximately every 5%.
INTEGER, PARAMETER :: target_number_of_updates = 20

! Progress state retained between calls
INTEGER, SAVE :: &
    expected_nstot = 0, & ! expected number of spin channels
    expected_nktot = 0, & ! expected number of k-points per spin channel
    total_blocks = 0, & ! total number of spin/k-point blocks
    report_interval = 1, & ! number of completed blocks between reports
    next_report = 1 ! next completed-block count that triggers a report
INTEGER(KIND=int64), SAVE :: &
    start_count = 0_int64, & ! system-clock count at initialization
    clock_rate = 0_int64, & ! system-clock counts per second
    clock_max = 0_int64 ! largest system-clock count
LOGICAL, SAVE :: &
    progress_initialized = .FALSE. ! progress state has been initialized

CONTAINS

SUBROUTINE initialize_progress( &
    nstot, nktot) ! in

IMPLICIT NONE

! Calculation dimensions
INTEGER, INTENT(IN) :: &
    nstot, & ! number of spin channels
    nktot ! number of k-points per spin channel

IF (nstot < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: Number of spin channels must be positive'
    ERROR STOP 1
END IF

IF (nktot < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: Number of k-points must be positive'
    ERROR STOP 1
END IF

expected_nstot = nstot
expected_nktot = nktot
total_blocks = nstot*nktot

! Use approximately target_number_of_updates reports over the complete
! calculation. Small calculations are reported after every block.
report_interval = MAX( &
    1, &
    CEILING(REAL(total_blocks,KIND=dp) / &
        REAL(target_number_of_updates,KIND=dp)))
next_report = report_interval

CALL SYSTEM_CLOCK( &
    count=start_count, &
    count_rate=clock_rate, &
    count_max=clock_max)

IF (clock_rate <= 0_int64) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid system-clock rate'
    ERROR STOP 1
END IF

progress_initialized = .TRUE.

WRITE(*,'(A,I0,A,I0,A)') &
    ' Progress: 0.0% (0/', total_blocks, &
    ' spin/k-point blocks, approximately ', &
    target_number_of_updates, ' updates planned)'
FLUSH(OUTPUT_UNIT)

END SUBROUTINE initialize_progress


SUBROUTINE write_progress( &
    ispin, nstot, ikpt, nktot) ! in

IMPLICIT NONE

! Current loop indices and dimensions
INTEGER, INTENT(IN) :: &
    ispin, & ! current spin-channel index
    nstot, & ! total number of spin channels
    ikpt, & ! current k-point index
    nktot ! total number of k-points per spin channel

! Progress counters
INTEGER :: &
    completed_blocks ! completed spin/k-point blocks
INTEGER(KIND=int64) :: &
    current_count ! current system-clock count

! Progress timing and completion
REAL(KIND=dp) :: &
    elapsed_time, & ! elapsed wall-clock time [s]
    percent_complete ! completed fraction [%]

IF (.NOT. progress_initialized) THEN
    WRITE(*,'(/,A)') ' ERROR: Progress reporting was not initialized'
    ERROR STOP 1
END IF

IF (nstot /= expected_nstot .OR. nktot /= expected_nktot) THEN
    WRITE(*,'(/,A)') ' ERROR: Progress dimensions differ from initialization'
    ERROR STOP 1
END IF

IF (ispin < 1 .OR. ispin > nstot) THEN
    WRITE(*,'(/,A)') ' ERROR: Spin index is outside the valid progress range'
    ERROR STOP 1
END IF

IF (ikpt < 1 .OR. ikpt > nktot) THEN
    WRITE(*,'(/,A)') ' ERROR: K-point index is outside the valid progress range'
    ERROR STOP 1
END IF

completed_blocks = (ispin-1)*nktot + ikpt

! Always report the first and final blocks. Between them, report only when
! the next approximately 5% threshold has been reached.
IF (completed_blocks /= 1 .AND. &
    completed_blocks < next_report .AND. &
    completed_blocks /= total_blocks) RETURN

percent_complete = 100.0_dp * &
    REAL(completed_blocks,KIND=dp) / REAL(total_blocks,KIND=dp)

CALL SYSTEM_CLOCK(count=current_count)

IF (current_count >= start_count) THEN
    elapsed_time = REAL(current_count-start_count,KIND=dp) / &
        REAL(clock_rate,KIND=dp)
ELSE
    elapsed_time = ( &
        REAL(clock_max-start_count,KIND=dp) + &
        REAL(current_count,KIND=dp) + 1.0_dp) / &
        REAL(clock_rate,KIND=dp)
END IF
elapsed_time = MAX(0.0_dp,elapsed_time)

WRITE(*,'(A,F6.1,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F12.1,A)') &
    ' Progress: ', percent_complete, '% (', &
    completed_blocks, '/', total_blocks, '), spin ', &
    ispin, '/', nstot, ', k-point ', ikpt, '/', nktot, &
    ', elapsed ', elapsed_time, ' s'
FLUSH(OUTPUT_UNIT)

DO WHILE (next_report <= completed_blocks)
    next_report = next_report + report_interval
END DO

END SUBROUTINE write_progress

END MODULE write_progress_mod
