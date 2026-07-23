MODULE write_oam_kpoint_mod

USE precision_mod, ONLY: dp

IMPLICIT NONE

PRIVATE
PUBLIC :: write_oam_kpoint

CONTAINS

SUBROUTINE write_oam_kpoint( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, oam)                 ! in

! Write band-resolved orbital angular momentum for one k-point. The routine is
! independent of the physical operators used to construct oam. It can therefore
! be reused for ordinary and spin-projected OAM output by passing the
! corresponding output-file unit.
!
! Output structure
! ----------------
! The k-point record is
!
!   # KP: <ikpt> NBANDS_OUT: <nbcderki> NBANDS_TRANS: <nbki>
!
! Each band record is
!
!   <band> <degenerate_block> <L_yz> <L_zx> <L_xy>
!
! No occupation-weighted total is written. This preserves the behavior of the
! original OAM output, which reported only band-resolved values.
!
! The widths of the band and degenerate-block integer fields are determined
! separately for every k-point from the largest value written in each column.
! Fixed-width integer descriptors are then constructed dynamically. This keeps
! all OAM columns aligned when an index gains another digit.
!
! Degenerate-block labels are local to the current k-point and spin channel.
! Bands carrying the same positive block ID were treated together by degenerate
! perturbation theory. Within such a block, the reported OAM values are the
! eigenvalues of the effective Hermitian OAM matrix and should not be
! interpreted as expectation values associated with the original unrotated DFT
! eigenvectors.
!
! Pseudovector convention
! -----------------------
! oam(:,1) = L_yz
! oam(:,2) = L_zx
! oam(:,3) = L_xy

! Connected formatted output unit and current k-point metadata
INTEGER, INTENT(IN) :: &
    unit, &      ! connected output-file unit
    ikpt, &      ! current k-point index
    nbki, &      ! number of transition-space bands at this k-point
    nbcderki     ! number of target bands written to the output

! Degenerate-group labels
INTEGER, INTENT(IN) :: &
    dg_group(:)  ! degenerate-block label for every transition-space band

! Band-resolved orbital-angular-momentum pseudovector components
REAL(KIND=dp), INTENT(IN) :: &
    oam(:,:)     ! (target band, yz/zx/xy)

! Loop index
INTEGER :: &
    iband

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

CALL validate_oam_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, oam)                 ! in

! Determine a fixed width for each integer column from the largest value that
! will be written for the current k-point. I0 is used only in these internal
! writes to count digits and to construct the final fixed-width format. It is
! not used for the band-resolved output records.
!
! OAM components are passed to formatted WRITE as three scalar items. A row
! slice such as oam(iband,1:3) is noncontiguous in Fortran column-major storage
! and causes Intel Fortran to create an array temporary.
WRITE(integer_text,'(I0)') nbcderki
band_width = LEN_TRIM(integer_text)

max_block = MAXVAL(dg_group(1:nbcderki))
WRITE(integer_text,'(I0)') max_block
block_width = LEN_TRIM(integer_text)

WRITE( &
    output_format, &
    '("(I",I0,",1X,I",I0,",3(1X,ES16.8E3))")') &
    band_width, block_width

WRITE( &
    unit, '(A,I0,1X,A,I0,1X,A,I0)', &
    IOSTAT=ierr, IOMSG=iomsg) &
    '# KP: ', ikpt, &
    'NBANDS_OUT: ', nbcderki, &
    'NBANDS_TRANS: ', nbki
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to write OAM k-point header'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

DO iband = 1, nbcderki
    WRITE( &
        unit, TRIM(output_format), &
        IOSTAT=ierr, IOMSG=iomsg) &
        iband, dg_group(iband), &
        oam(iband,1), oam(iband,2), oam(iband,3)
    IF (ierr /= 0) THEN
        WRITE(*,'(/,A)') ' ERROR: Unable to write band OAM'
        WRITE(*,'(A,I0)') ' K-point = ', ikpt
        WRITE(*,'(A,I0)') ' Band = ', iband
        WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
        ERROR STOP 1
    END IF
END DO

FLUSH(unit, IOSTAT=ierr, IOMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to flush OAM output'
    WRITE(*,'(A,I0)') ' K-point = ', ikpt
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

END SUBROUTINE write_oam_kpoint


SUBROUTINE validate_oam_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, oam)                 ! in

! Validate the dimensions and metadata required to write one k-point.

! Connected formatted output unit and current k-point metadata
INTEGER, INTENT(IN) :: &
    unit, &
    ikpt, &
    nbki, &
    nbcderki

! Degenerate-group labels
INTEGER, INTENT(IN) :: &
    dg_group(:)

! Band-resolved orbital-angular-momentum pseudovector components
REAL(KIND=dp), INTENT(IN) :: &
    oam(:,:)

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
    WRITE(*,'(/,A)') ' ERROR: Unable to inquire about OAM output unit'
    WRITE(*,'(A)') ' Reason: '//TRIM(iomsg)
    ERROR STOP 1
END IF

IF (.NOT. opened) THEN
    WRITE(*,'(/,A)') ' ERROR: OAM output unit is not open'
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
    WRITE(*,'(/,A)') ' ERROR: Invalid number of OAM output bands'
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

IF (ANY(SHAPE(oam) /= [nbcderki,3])) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected OAM array shape'
    WRITE(*,'(A,2(I0,1X))') ' Actual shape = ', SHAPE(oam)
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

END SUBROUTINE validate_oam_output_inputs

END MODULE write_oam_kpoint_mod
