MODULE write_goam_kpoint_mod

USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
USE precision_mod, ONLY: dp

IMPLICIT NONE

PRIVATE
PUBLIC :: write_goam_kpoint

CONTAINS

SUBROUTINE write_goam_kpoint( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, goam)                ! in

! Write the generalized orbital-angular-momentum matrices for one k-point.
!
! MATRIX CONVENTION
! -----------------
!
! The input array follows
!
!     goam(i,j,icomp) = <u_i|L_icomp|u_j>,
!
! where i is the bra-band index, j is the ket-band index, and
!
!     icomp = 1  ->  L_yz  (Voigt 4),
!     icomp = 2  ->  L_zx  (Voigt 5),
!     icomp = 3  ->  L_xy  (Voigt 6).
!
! Both external matrix indices span the smaller target-band range
! 1:nbcderki. The larger transition-band range 1:nbki has already been used
! by calculate_goam_kpoint when summing over intermediate states. It is
! reported in the k-point header so that the calculation range is explicit.
!
! OUTPUT STRUCTURE
! ----------------
!
! The k-point record is
!
!   # KP: <ikpt> NBANDS_OUT: <nbcderki> NBANDS_TRANS: <nbki>
!
! Each component begins with, for example,
!
!   # COMPONENT: yz VOIGT: 4
!   # COLUMNS: ket bands 1:<nbcderki>
!
! Each following matrix row is
!
!   <bra band> <degenerate block> (Re L_i1,Im L_i1) ... (Re L_iN,Im L_iN)
!
! The band and block fields have separately determined fixed widths for the
! current k-point. The matrix columns therefore remain aligned when either
! integer index gains another digit. The complex matrix elements are written
! from scalar expressions inside an implied DO loop, avoiding a noncontiguous
! array section in the input/output list.
!
! The degenerate-block labels are local to the current k-point and spin
! channel. They identify the degenerate subspaces used when excluding
! singular intermediate-state couplings in calculate_goam_kpoint. The matrix
! itself remains in the original DFT eigenstate basis.
!
! HERMITICITY
! -----------
!
! calculate_goam_kpoint constructs a Hermitian matrix,
!
!     goam(j,i,icomp) = CONJG(goam(i,j,icomp)).
!
! The writer outputs the complete matrix rather than only one triangle. This
! makes the file directly usable without requiring a reconstruction step.

IMPLICIT NONE

! Connected output unit and current k-point dimensions
INTEGER, INTENT(IN) :: &
    unit, &      ! connected formatted output unit
    ikpt, &      ! current k-point index
    nbki, &      ! number of transition bands used in intermediate-state sums
    nbcderki     ! number of external bands written as matrix rows and columns

! Degenerate-block assignment for the external bands
INTEGER, INTENT(IN) :: &
    dg_group(:)  ! local degenerate-block labels

! Generalized orbital-angular-momentum matrices
COMPLEX(KIND=dp), INTENT(IN) :: &
    goam(:,:,:)  ! (bra band, ket band, yz/zx/xy) [hbar]

! Component labels and corresponding Voigt indices
CHARACTER(LEN=2), PARAMETER :: &
    component_label(3) = ['yz', 'zx', 'xy']
INTEGER, PARAMETER :: &
    voigt_index(3) = [4, 5, 6]

! Dynamically constructed row format and integer-width text
CHARACTER(LEN=160) :: &
    output_format
CHARACTER(LEN=32) :: &
    band_width_text, &
    block_width_text, &
    matrix_size_text

! Input/output status information
INTEGER :: &
    ierr
CHARACTER(LEN=512) :: &
    iomsg

! Loop indices and integer-field widths
INTEGER :: &
    band_width, &
    block_width, &
    i, &
    icomp, &
    j, &
    max_block

CALL validate_goam_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, goam)                ! in

! Determine separately the fixed widths needed by the band and block fields.
WRITE(band_width_text,'(I0)') nbcderki
band_width = LEN_TRIM(band_width_text)

max_block = MAXVAL(dg_group(1:nbcderki))
WRITE(block_width_text,'(I0)') max_block
block_width = LEN_TRIM(block_width_text)

WRITE(band_width_text,'(I0)') band_width
WRITE(block_width_text,'(I0)') block_width
WRITE(matrix_size_text,'(I0)') nbcderki

output_format = &
    '(I' // TRIM(band_width_text) // ',1X,I' // &
    TRIM(block_width_text) // ',' // TRIM(matrix_size_text) // &
    '(1X,"(",ES16.8E3,",",ES16.8E3,")"))'

WRITE( &
    unit, '(A,I0,1X,A,I0,1X,A,I0)', &
    IOSTAT=ierr, IOMSG=iomsg) &
    '# KP: ', ikpt, 'NBANDS_OUT: ', nbcderki, 'NBANDS_TRANS: ', nbki
CALL check_write_status( &
    ierr, iomsg, ikpt, 0, 0) ! in

DO icomp = 1, 3
    WRITE( &
        unit, '(A,A,1X,A,I0)', &
        IOSTAT=ierr, IOMSG=iomsg) &
        '# COMPONENT: ', component_label(icomp), &
        'VOIGT: ', voigt_index(icomp)
    CALL check_write_status( &
        ierr, iomsg, ikpt, icomp, 0) ! in

    WRITE( &
        unit, '(A,I0)', &
        IOSTAT=ierr, IOMSG=iomsg) &
        '# COLUMNS: ket bands 1:', nbcderki
    CALL check_write_status( &
        ierr, iomsg, ikpt, icomp, 0) ! in

    DO i = 1, nbcderki
        WRITE( &
            unit, TRIM(output_format), &
            IOSTAT=ierr, IOMSG=iomsg) &
            i, dg_group(i), &
            (REAL(goam(i,j,icomp),KIND=dp), &
             AIMAG(goam(i,j,icomp)), j=1,nbcderki)

        CALL check_write_status( &
            ierr, iomsg, ikpt, icomp, i) ! in
    END DO
END DO

END SUBROUTINE write_goam_kpoint


SUBROUTINE validate_goam_output_inputs( &
    unit, ikpt, nbki, nbcderki, & ! in
    dg_group, goam)                ! in

! Validate dimensions, labels, values, and the output-unit connection.

IMPLICIT NONE

! Connected output unit and current k-point dimensions
INTEGER, INTENT(IN) :: &
    unit, &
    ikpt, &
    nbki, &
    nbcderki

! Degenerate-block assignment
INTEGER, INTENT(IN) :: &
    dg_group(:)

! Generalized orbital-angular-momentum matrices
COMPLEX(KIND=dp), INTENT(IN) :: &
    goam(:,:,:)

! Output-unit connection status
LOGICAL :: &
    unit_is_open

! Loop indices for finite-value validation
INTEGER :: &
    i, &
    icomp, &
    j

INQUIRE(UNIT=unit,OPENED=unit_is_open)
IF (.NOT. unit_is_open) THEN
    WRITE(*,'(/,A)') ' ERROR: Generalized-OAM output unit is not open'
    WRITE(*,'(A,I0)') ' Unit = ', unit
    ERROR STOP 1
END IF

IF (ikpt < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid k-point index for generalized-OAM output'
    WRITE(*,'(A,I0)') ' ikpt = ', ikpt
    ERROR STOP 1
END IF

IF (nbki < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid transition-band count for generalized-OAM output'
    WRITE(*,'(A,I0)') ' nbki = ', nbki
    ERROR STOP 1
END IF

IF (nbcderki < 1 .OR. nbcderki > nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid external-band count for generalized-OAM output'
    WRITE(*,'(A,I0)') ' nbcderki = ', nbcderki
    WRITE(*,'(A,I0)') ' nbki = ', nbki
    ERROR STOP 1
END IF

IF (SIZE(dg_group) < nbcderki) THEN
    WRITE(*,'(/,A)') ' ERROR: Degenerate-group array is too small for generalized-OAM output'
    WRITE(*,'(A,I0)') ' SIZE(dg_group) = ', SIZE(dg_group)
    WRITE(*,'(A,I0)') ' Required size = ', nbcderki
    ERROR STOP 1
END IF

IF (ANY(dg_group(1:nbcderki) < 1)) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid degenerate-block label in generalized-OAM output'
    ERROR STOP 1
END IF

IF (SIZE(goam,1) /= nbcderki .OR. &
    SIZE(goam,2) /= nbcderki .OR. &
    SIZE(goam,3) /= 3) THEN
    WRITE(*,'(/,A)') ' ERROR: Generalized-OAM array has an unexpected shape'
    WRITE(*,'(A,3(1X,I0))') &
        ' Actual shape = ', SIZE(goam,1), SIZE(goam,2), SIZE(goam,3)
    WRITE(*,'(A,3(1X,I0))') &
        ' Expected shape = ', nbcderki, nbcderki, 3
    ERROR STOP 1
END IF

DO icomp = 1, 3
    DO j = 1, nbcderki
        DO i = 1, nbcderki
            IF (.NOT. IEEE_IS_FINITE(REAL(goam(i,j,icomp),KIND=dp)) .OR. &
                .NOT. IEEE_IS_FINITE(AIMAG(goam(i,j,icomp)))) THEN
                WRITE(*,'(/,A)') ' ERROR: Non-finite generalized-OAM output value'
                WRITE(*,'(A,I0)') ' Bra-band index = ', i
                WRITE(*,'(A,I0)') ' Ket-band index = ', j
                WRITE(*,'(A,I0)') ' Component index = ', icomp
                ERROR STOP 1
            END IF
        END DO
    END DO
END DO

END SUBROUTINE validate_goam_output_inputs


SUBROUTINE check_write_status( &
    ierr, iomsg, ikpt, icomp, iband) ! in

! Report a formatted-output failure with the current matrix location.

IMPLICIT NONE

! Input/output status and current output location
INTEGER, INTENT(IN) :: &
    ierr, &
    ikpt, &
    icomp, &
    iband
CHARACTER(LEN=*), INTENT(IN) :: &
    iomsg

IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to write generalized orbital angular momentum'
    WRITE(*,'(A,I0)') ' K-point index = ', ikpt
    IF (icomp > 0) WRITE(*,'(A,I0)') ' Component index = ', icomp
    IF (iband > 0) WRITE(*,'(A,I0)') ' Bra-band index = ', iband
    WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
    ERROR STOP 1
END IF

END SUBROUTINE check_write_status

END MODULE write_goam_kpoint_mod
