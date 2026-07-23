MODULE calculate_goam_kpoint_mod

USE precision_mod, ONLY: sp, dp
USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_IS_FINITE

IMPLICIT NONE

PRIVATE
PUBLIC :: calculate_goam_kpoint

CONTAINS

SUBROUTINE calculate_goam_kpoint( &
    nbki, nbcderki, &               ! in
    dg_group, &                     ! in
    pijA, pijB, dEij, &             ! in
    goam)                           ! out

! Calculate the generalized orbital-angular-momentum matrix at one k-point.
!
! The returned matrix is restricted to the target-band range:
!
!     i,j = 1,...,nbcderki,
!
! while the intermediate-state sum retains the complete transition range:
!
!     l = 1,...,nbki.
!
! Therefore, bands above nbcderki do not appear as external row or column
! indices, but their virtual couplings are included in every matrix element.
!
! MATRIX-ELEMENT CONVENTION
! -------------------------
!
! The rectangular input arrays store
!
!     pijA(alpha,i,l) = <u_i|A_alpha|u_l>,
!     pijB(beta ,i,l) = <u_i|B_beta |u_l>,
!
! for target state i <= nbcderki and transition state l <= nbki. The reverse
! orientation required below is reconstructed from Hermiticity:
!
!     <u_l|B_beta |u_j> = CONJG(pijB(beta,j,l)),
!     <u_l|A_alpha|u_j> = CONJG(pijA(alpha,j,l)).
!
! The output convention is
!
!     goam(i,j,icomp) = <u_i|L_icomp|u_j>,
!
! with components
!
!     icomp = 1: L_yz,
!     icomp = 2: L_zx,
!     icomp = 3: L_xy.
!
! GENERALIZED OAM FORMULA
! -----------------------
!
! Define
!
!     P_ijl^(alpha,beta) =
!         <u_i|A_alpha|u_l><u_l|B_beta|u_j>
!       - <u_i|B_beta |u_l><u_l|A_alpha|u_j>.
!
! With
!
!     dEij(i,l) = E_l-E_i,
!
! the generalized OAM matrix is evaluated as
!
!     L_ij^(alpha,beta) = (i/2) SUM_l P_ijl^(alpha,beta)
!         [ Q_i(l)/(E_l-E_i) + Q_j(l)/(E_l-E_j) ],
!
! where Q_i(l) is zero when l belongs to the same degenerate block as i and
! one otherwise. Q_j(l) is defined analogously. These projectors exclude
! singular couplings within either external state's degenerate subspace.
!
! ISOLATED-BAND DIAGONAL
! ----------------------
!
! For an isolated state i, setting j=i gives
!
!     L_ii = i SUM_(l /= i) P_iil/(E_l-E_i)
!
!          = 2 SUM_(l /= i)
!              Im[<u_i|A_alpha|u_l><u_l|B_beta|u_i>]/(E_i-E_l),
!
! which is exactly the validated non-degenerate formula used by the regular
! OAM calculation. The diagonal is implemented as a separate branch so that
! it follows the same arithmetic sequence as the regular OAM path rather than
! evaluating two nominally identical half-terms.
!
! HERMITICITY
! -----------
!
! Only the upper triangle is calculated explicitly. The lower triangle is
! reconstructed as
!
!     goam(j,i,icomp) = CONJG(goam(i,j,icomp)).
!
! The diagonal is explicitly reduced to its real part after accumulation.
!
! PRECISION AND SUMMATION
! -----------------------
!
! The momentum matrix elements and energy differences enter in sp, matching
! the readers and the regular OAM implementation. To reproduce the validated
! regular-OAM diagonal as closely as possible:
!
!   1. each antisymmetrized momentum product is formed in sp;
!   2. each stored energy difference is read in sp;
!   3. numerator and denominator are promoted separately to dp;
!   4. generalized-OAM terms are evaluated in dp;
!   5. complex Kahan compensated summation is performed in dp;
!   6. the returned generalized-OAM matrix remains COMPLEX(KIND=dp).
!
! No dp result is converted back to sp.

! Current k-point dimensions and degenerate-group assignment
INTEGER, INTENT(IN) :: &
    nbki, &      ! number of transition bands at the current k-point and spin
    nbcderki     ! number of external target bands at the current k-point
INTEGER, INTENT(IN) :: &
    dg_group(:)  ! contiguous degenerate-block label for every transition band

! Rectangular operator matrices and energy differences
COMPLEX(KIND=sp), INTENT(IN) :: &
    pijA(:,:,:), & ! <u_i|A_alpha|u_l> [a.u.]
    pijB(:,:,:)    ! <u_i|B_beta|u_l> [a.u.]
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:)      ! E_l-E_i for target i and transition state l [Ha]

! Generalized orbital-angular-momentum matrices
COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT) :: &
    goam(:,:,:)    ! (bra band, ket band, yz/zx/xy) [hbar]

! Kahan summation values for one matrix element
COMPLEX(KIND=dp) :: &
    goam_sum, &       ! accumulated generalized-OAM matrix element
    correction, & ! Kahan compensation
    term          ! generalized-OAM contribution from one intermediate band

! Antisymmetrized momentum product
COMPLEX(KIND=sp) :: &
    p2                ! P_ijl formed in single precision

! Double-precision reciprocal denominator factor
REAL(KIND=dp) :: &
    inverse_dE_sum    ! Q_i/dE_i + Q_j/dE_j

! Loop indices and Cartesian component directions
INTEGER :: &
    icomp, & ! pseudovector component: 1=yz, 2=zx, 3=xy
    alpha, & ! first Cartesian operator direction
    beta, &  ! second Cartesian operator direction
    i, &     ! external bra-band index
    j, &     ! external ket-band index
    l, &     ! intermediate transition-band index
    ierr     ! allocation status

! Inclusion flag for one off-diagonal intermediate-state contribution
LOGICAL :: &
    has_contribution

! Allocation diagnostic
CHARACTER(LEN=512) :: &
    iomsg

! Cartesian pairs corresponding to {L_yz, L_zx, L_xy}
INTEGER, PARAMETER :: &
    alpha_component(3) = [2, 3, 1], &
    beta_component(3)  = [3, 1, 2]

CALL validate_goam_inputs( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij)           ! in

ALLOCATE( &
    goam(nbcderki,nbcderki,3), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') &
        ' ERROR: Unable to allocate generalized orbital-angular-momentum array'
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF
goam = (0.0_dp,0.0_dp)

DO icomp = 1, 3
    alpha = alpha_component(icomp)
    beta = beta_component(icomp)

    DO i = 1, nbcderki
        DO j = i, nbcderki

            goam_sum = (0.0_dp,0.0_dp)
            correction = (0.0_dp,0.0_dp)

            IF (i == j) THEN

                ! Follow the validated regular-OAM diagonal arithmetic exactly:
                !
                !     L_ii = i SUM_l P_iil/dEij(i,l).
                !
                ! The same-block test also excludes l=i.
                DO l = 1, nbki
                    IF (dg_group(l) /= dg_group(i)) THEN
                        CALL check_nonzero_denominator( &
                            i, l, dEij(i,l)) ! in

                        p2 = &
                            pijA(alpha,i,l) * CONJG(pijB(beta,i,l)) - &
                            pijB(beta,i,l) * CONJG(pijA(alpha,i,l))

                        term = (0.0_dp,1.0_dp) * &
                            CMPLX(p2,KIND=dp) / REAL(dEij(i,l),KIND=dp)

                        CALL check_finite_term( &
                            i, j, l, icomp, term) ! in

                        CALL add_kahan_complex( &
                            term, &                ! in
                            goam_sum, correction) ! in/out
                    END IF
                END DO

                ! A Hermitian operator has a real diagonal. Remove only the
                ! residual imaginary roundoff after the complete dp sum.
                goam(i,i,icomp) = &
                    CMPLX(REAL(goam_sum,KIND=dp),0.0_dp,KIND=dp)

            ELSE

                ! General off-diagonal matrix element:
                !
                !   L_ij = (i/2) SUM_l P_ijl
                !       [Q_i(l)/dEij(i,l) + Q_j(l)/dEij(j,l)].
                DO l = 1, nbki
                    inverse_dE_sum = 0.0_dp
                    has_contribution = .FALSE.

                    IF (dg_group(l) /= dg_group(i)) THEN
                        CALL check_nonzero_denominator( &
                            i, l, dEij(i,l)) ! in
                        inverse_dE_sum = inverse_dE_sum + &
                            1.0_dp / REAL(dEij(i,l),KIND=dp)
                        has_contribution = .TRUE.
                    END IF

                    IF (dg_group(l) /= dg_group(j)) THEN
                        CALL check_nonzero_denominator( &
                            j, l, dEij(j,l)) ! in
                        inverse_dE_sum = inverse_dE_sum + &
                            1.0_dp / REAL(dEij(j,l),KIND=dp)
                        has_contribution = .TRUE.
                    END IF

                    IF (has_contribution) THEN
                        p2 = &
                            pijA(alpha,i,l) * CONJG(pijB(beta,j,l)) - &
                            pijB(beta,i,l) * CONJG(pijA(alpha,j,l))

                        term = (0.0_dp,0.5_dp) * &
                            CMPLX(p2,KIND=dp) * inverse_dE_sum

                        CALL check_finite_term( &
                            i, j, l, icomp, term) ! in

                        CALL add_kahan_complex( &
                            term, &                ! in
                            goam_sum, correction) ! in/out
                    END IF
                END DO

                goam(i,j,icomp) = goam_sum
                goam(j,i,icomp) = CONJG(goam_sum)

            END IF
        END DO
    END DO
END DO

END SUBROUTINE calculate_goam_kpoint


SUBROUTINE add_kahan_complex( &
    term, &                  ! in
    accumulated, correction) ! in/out

! Add one complex dp term using Kahan compensated summation.

! Term to add
COMPLEX(KIND=dp), INTENT(IN) :: &
    term

! Accumulator and Kahan correction
COMPLEX(KIND=dp), INTENT(INOUT) :: &
    accumulated, &
    correction

! Kahan working values
COMPLEX(KIND=dp) :: &
    corrected_term, &
    temporary

corrected_term = term - correction
temporary = accumulated + corrected_term
correction = (temporary - accumulated) - corrected_term
accumulated = temporary

END SUBROUTINE add_kahan_complex


SUBROUTINE check_nonzero_denominator( &
    iband, lband, denominator) ! in

! Guard against a zero denominator outside the identified degenerate block.

! Band indices
INTEGER, INTENT(IN) :: &
    iband, &
    lband

! Stored energy difference
REAL(KIND=sp), INTENT(IN) :: &
    denominator

IF (ABS(denominator) <= TINY(denominator)) THEN
    WRITE(*,'(/,A)') &
        ' ERROR: Zero generalized-OAM denominator outside a degenerate block'
    WRITE(*,'(A,I0)') ' External band = ', iband
    WRITE(*,'(A,I0)') ' Intermediate band = ', lband
    ERROR STOP 1
END IF

END SUBROUTINE check_nonzero_denominator


SUBROUTINE check_finite_term( &
    i, j, l, icomp, term) ! in

! Verify that one generalized-OAM contribution is finite.

! Matrix and component indices
INTEGER, INTENT(IN) :: &
    i, &
    j, &
    l, &
    icomp

! Generalized-OAM contribution
COMPLEX(KIND=dp), INTENT(IN) :: &
    term

IF (.NOT. IEEE_IS_FINITE(REAL(term,KIND=dp)) .OR. &
        .NOT. IEEE_IS_FINITE(AIMAG(term))) THEN
    WRITE(*,'(/,A)') ' ERROR: Non-finite generalized-OAM contribution'
    WRITE(*,'(A,I0)') ' Bra-band index = ', i
    WRITE(*,'(A,I0)') ' Ket-band index = ', j
    WRITE(*,'(A,I0)') ' Intermediate-band index = ', l
    WRITE(*,'(A,I0)') ' Component index = ', icomp
    WRITE(*,'(A,2ES16.8E3)') &
        ' Contribution = ', REAL(term,KIND=dp), AIMAG(term)
    ERROR STOP 1
END IF

END SUBROUTINE check_finite_term


SUBROUTINE validate_goam_inputs( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij)           ! in

! Validate the interface and degenerate-block invariants required here.

! Current k-point dimensions and group labels
INTEGER, INTENT(IN) :: &
    nbki, &
    nbcderki
INTEGER, INTENT(IN) :: &
    dg_group(:)

! Rectangular operator matrices and energy differences
COMPLEX(KIND=sp), INTENT(IN) :: &
    pijA(:,:,:), &
    pijB(:,:,:)
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:)

! Validation loop index
INTEGER :: &
    iband

IF (nbki < 1) THEN
    WRITE(*,'(/,A)') ' ERROR: nbki must be positive'
    ERROR STOP 1
END IF

IF (nbcderki < 1 .OR. nbcderki > nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid generalized-OAM target-band limit'
    WRITE(*,'(A,I0)') ' nbcderki = ', nbcderki
    WRITE(*,'(A,I0)') ' nbki = ', nbki
    ERROR STOP 1
END IF

IF (SIZE(dg_group) /= nbki) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected dg_group size'
    WRITE(*,'(A,I0)') ' Actual size = ', SIZE(dg_group)
    WRITE(*,'(A,I0)') ' Expected size = ', nbki
    ERROR STOP 1
END IF

IF (ANY(SHAPE(pijA) /= [3,nbcderki,nbki])) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected shape of pijA'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(pijA)
    WRITE(*,'(A,3(I0,1X))') &
        ' Expected shape = ', 3, nbcderki, nbki
    ERROR STOP 1
END IF

IF (ANY(SHAPE(pijB) /= [3,nbcderki,nbki])) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected shape of pijB'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(pijB)
    WRITE(*,'(A,3(I0,1X))') &
        ' Expected shape = ', 3, nbcderki, nbki
    ERROR STOP 1
END IF

IF (ANY(SHAPE(dEij) /= [nbcderki,nbki])) THEN
    WRITE(*,'(/,A)') ' ERROR: Unexpected shape of dEij'
    WRITE(*,'(A,2(I0,1X))') ' Actual shape = ', SHAPE(dEij)
    WRITE(*,'(A,2(I0,1X))') &
        ' Expected shape = ', nbcderki, nbki
    ERROR STOP 1
END IF

IF (ANY(dg_group < 1)) THEN
    WRITE(*,'(/,A)') ' ERROR: Degenerate-group labels must be positive'
    ERROR STOP 1
END IF

DO iband = 2, nbki
    IF (dg_group(iband) < dg_group(iband-1) .OR. &
            dg_group(iband) > dg_group(iband-1)+1) THEN
        WRITE(*,'(/,A)') &
            ' ERROR: Degenerate-group labels are not contiguous'
        WRITE(*,'(A,I0)') ' Band index = ', iband
        WRITE(*,'(A,I0)') &
            ' Previous group = ', dg_group(iband-1)
        WRITE(*,'(A,I0)') ' Current group = ', dg_group(iband)
        ERROR STOP 1
    END IF
END DO

IF (nbcderki < nbki) THEN
    IF (dg_group(nbcderki) == dg_group(nbcderki+1)) THEN
        WRITE(*,'(/,A)') &
            ' ERROR: nbcderki splits a degenerate group'
        WRITE(*,'(A,I0)') ' nbcderki = ', nbcderki
        WRITE(*,'(A,I0)') &
            ' Split group = ', dg_group(nbcderki)
        ERROR STOP 1
    END IF
END IF

END SUBROUTINE validate_goam_inputs

END MODULE calculate_goam_kpoint_mod
