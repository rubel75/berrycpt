MODULE write_bcurv_kpoint_mod

USE precision_mod, ONLY: sp, dp

IMPLICIT NONE

PRIVATE
PUBLIC :: write_bcurv_kpoint

CONTAINS

SUBROUTINE write_bcurv_kpoint( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, occup, bcurv)        ! in

! Write band-resolved Berry curvature and its occupation-weighted total for
! one k-point. The routine is independent of the physical operator used to
! construct bcurv. It can therefore be reused for ordinary, spin-resolved,
! and spin Berry-curvature output by passing the corresponding file unit.
!
! Output structure
! ----------------
! The k-point record is
!
!   # KP: <ikpt> NBANDS_OUT: <nbcderki> NBANDS_TRANS: <nbki>
!
! Each band record is
!
!   <band> <degenerate_block> <Omega_yz> <Omega_zx> <Omega_xy>
!
! The final record for the k-point is
!
!   0 0 <total_yz> <total_zx> <total_xy>
!
! The widths of the band and degenerate-block integer fields are determined
! separately for every k-point from the largest value written in each column.
! Fixed-width integer descriptors are then constructed dynamically. This keeps
! all Berry-curvature columns aligned when an index gains another digit.
!
! where the two zero indices identify the occupation-weighted total rather
! than a physical band or degenerate block. Degenerate-block labels are local
! to the current k-point and spin channel. Bands carrying the same label were
! treated together by degenerate perturbation theory. Within such a block,
! the reported Berry curvatures are eigenvalues of the effective Hermitian
! Berry-curvature matrix and should not be interpreted as expectation values
! associated with the original unrotated DFT eigenvectors.
!
! Total Berry curvature
! ---------------------
! For each pseudovector component c, the total is
!
!   Omega_total(c) = SUM_n f_n Omega_n(c),
!
! where f_n is the DFT or constructed occupation. Occupations enter as sp and
! are explicitly converted to dp before multiplication by the dp Berry
! curvature. The sum is accumulated in dp using Kahan compensated summation.
!
! Pseudovector convention
! -----------------------
! bcurv(:,1) = Omega_yz
! bcurv(:,2) = Omega_zx
! bcurv(:,3) = Omega_xy

! Connected formatted output unit and current k-point metadata
INTEGER, INTENT(IN) :: &
    unit, &      ! connected output-file unit
    ikpt, &      ! current k-point index
    nbki, &      ! number of transition-space bands at this k-point
    nbcderki     ! number of target bands written to the output

! Degenerate-group labels and band occupations
INTEGER, INTENT(IN) :: &
    dg_group(:)  ! degenerate-block label for every transition-space band
REAL(KIND=sp), INTENT(IN) :: &
    occup(:)     ! occupation of every target band

! Band-resolved Berry-curvature pseudovector components
REAL(KIND=dp), INTENT(IN) :: &
    bcurv(:,:)   ! (target band, yz/zx/xy) [bohr^2]

! Occupation-weighted total Berry curvature
REAL(KIND=dp) :: &
    bcurv_total(3)

! Kahan compensated-summation variables
REAL(KIND=dp) :: &
    sum_value, &  ! running compensated sum
    correction, & ! compensation for lost low-order bits
    term, &       ! current occupation-weighted contribution
    corrected, &  ! contribution after applying the compensation
    updated       ! updated running sum

! Loop indices
INTEGER :: &
    iband, & ! target-band index
    icomp    ! pseudovector-component index

! Dynamically determined integer-field widths
INTEGER :: &
    band_width, &  ! number of digits required for the largest band index
    block_width, & ! number of digits required for the largest block index
    max_block      ! largest degenerate-block index written at this k-point

! Dynamically constructed output format
CHARACTER(LEN=128) :: &
    output_format
CHARACTER(LEN=32) :: &
    integer_text

! Input/output status information
INTEGER :: &
    ierr
CHARACTER(LEN=512) :: &
    iomsg

CALL validate_bcurv_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, occup, bcurv)        ! in

! Determine a fixed width for each integer column from the largest value that
! will be written for the current k-point. I0 is used only in these internal
! writes to count digits and to construct the final fixed-width format. It is
! not used for the band-resolved output records.
!
! Berry-curvature components are passed to formatted WRITE as three scalar
! items. A row slice such as bcurv(iband,1:3) is noncontiguous in Fortran
! column-major storage and causes Intel Fortran to create an array temporary.
WRITE(integer_text,'(I0)') nbcderki
band_width = LEN_TRIM(integer_text)

max_block = MAXVAL(dg_group(1:nbcderki))
WRITE(integer_text,'(I0)') max_block
block_width = LEN_TRIM(integer_text)

WRITE( &
    output_format, &
    '("(I",I0,",1X,I",I0,",3(1X,ES16.8E3))")') &
    band_width, block_width

! Form the occupation-weighted total in dp. Kahan compensation is applied
! independently to each pseudovector component.
bcurv_total = 0.0_dp
DO icomp = 1, 3
    sum_value = 0.0_dp
    correction = 0.0_dp

    DO iband = 1, nbcderki
        term = REAL(occup(iband),KIND=dp) * bcurv(iband,icomp)
        corrected = term - correction
        updated = sum_value + corrected
        correction = (updated - sum_value) - corrected
        sum_value = updated
    END DO

    bcurv_total(icomp) = sum_value
END DO

WRITE( &
    unit, '(A,I0,1X,A,I0,1X,A,I0)', &
    IOSTAT=ierr, IOMSG=iomsg) &
    '# KP: ', ikpt, &
    'NBANDS_OUT: ', nbcderki, &
    'NBANDS_TRANS: ', nbki
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to write Berry-curvature k-point header'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

DO iband = 1, nbcderki
    WRITE( &
        unit, TRIM(output_format), &
        IOSTAT=ierr, IOMSG=iomsg) &
        iband, dg_group(iband), &
        bcurv(iband,1), bcurv(iband,2), bcurv(iband,3)
    IF (ierr /= 0) THEN
        WRITE(*,'(/,A)') ' ERROR: Unable to write band Berry curvature'
        WRITE(*,'(A,I0)') ' K-point = ', ikpt
        WRITE(*,'(A,I0)') ' Band = ', iband
        WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
        ERROR STOP 1
    END IF
END DO

WRITE( &
    unit, TRIM(output_format), &
    IOSTAT=ierr, IOMSG=iomsg) &
    0, 0, &
    bcurv_total(1), bcurv_total(2), bcurv_total(3)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to write total Berry curvature'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

FLUSH(unit, IOSTAT=ierr, IOMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to flush Berry-curvature output'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF


END SUBROUTINE write_bcurv_kpoint


SUBROUTINE validate_bcurv_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, occup, bcurv)        ! in

! Validate the dimensions and metadata required to write one k-point.

! Connected formatted output unit and current k-point metadata
INTEGER, INTENT(IN) :: &
    unit, &
    ikpt, &
    nbki, &
    nbcderki

! Degenerate-group labels and band occupations
INTEGER, INTENT(IN) :: &
    dg_group(:)
REAL(KIND=sp), INTENT(IN) :: &
    occup(:)

! Band-resolved Berry-curvature pseudovector components
REAL(KIND=dp), INTENT(IN) :: &
    bcurv(:,:)

! Validation loop index and input/output status
INTEGER :: &
    iband, &
    ierr
LOGICAL :: &
    opened
CHARACTER(LEN=512) :: &
    iomsg

INQUIRE( &
    UNIT=unit, &
    OPENED=opened, &
    IOSTAT=ierr, &
    IOMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to inquire about Berry-curvature output unit'
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

IF (.NOT. opened) THEN
    WRITE(*,'(/,A)') ' ERROR: Berry-curvature output unit is not open'
    WRITE(*,'(A,I0)') ' Unit = ', unit
    ERROR STOP 1
END IF

IF (ikpt < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: K-point index must be positive'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    ERROR STOP 1
END IF

IF (nbki < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: Number of transition-space bands must be positive'
    WRITE(*,'(A,I0)') ' nbki = ', nbki
    ERROR STOP 1
END IF

IF (nbcderki < 1 .OR. nbcderki > nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid number of Berry-curvature output bands'
    WRITE(*,'(A,I0)') ' nbcderki = ', nbcderki
    WRITE(*,'(A,I0)') ' nbki = ', nbki
    ERROR STOP 1
END IF

IF (SIZE(dg_group) /= nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected degenerate-group array size'
    WRITE(*,'(A,I0)') ' Actual size = ', SIZE(dg_group)
    WRITE(*,'(A,I0)') ' Expected size = ', nbki
    ERROR STOP 1
END IF

IF (SIZE(occup) /= nbcderki) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected occupation array size'
    WRITE(*,'(A,I0)') ' Actual size = ', SIZE(occup)
    WRITE(*,'(A,I0)') ' Expected size = ', nbcderki
    ERROR STOP 1
END IF

IF (ANY(SHAPE(bcurv) /= [nbcderki,3])) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected Berry-curvature array shape'
    WRITE(*,'(A,2(I0,1X))') ' Actual shape = ', SHAPE(bcurv)
    WRITE(*,'(A,2(I0,1X))') ' Expected shape = ', nbcderki, 3
    ERROR STOP 1
END IF

IF (dg_group(1) /= 1) THEN
    WRITE(*,'(/,A)') ' ERROR: The first degenerate-block label must be one'
    WRITE(*,'(A,I0)') ' First label = ', dg_group(1)
    ERROR STOP 1
END IF

DO iband = 2, nbki
    IF (dg_group(iband) < dg_group(iband-1) .OR. &
            dg_group(iband) > dg_group(iband-1)+1) THEN
        WRITE(*,'(/,A)') ' ERROR: Degenerate-block labels are not contiguous'
        WRITE(*,'(A,I0)') ' Band = ', iband
        WRITE(*,'(A,I0)') ' Previous block = ', dg_group(iband-1)
        WRITE(*,'(A,I0)') ' Current block = ', dg_group(iband)
        ERROR STOP 1
    END IF
END DO

END SUBROUTINE validate_bcurv_output_inputs

END MODULE write_bcurv_kpoint_mod
