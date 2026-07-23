MODULE degenoam_mod
CONTAINS

SUBROUTINE degenoam(nb, idg1, idg2, pijA, pijB, dEij, & ! <- args in 
    oam) ! -> args out

! Orbital angular momentum (OAM) for a block of degenerate bands

USE precision_mod, ONLY: sp, dp
USE eigvz_mod, ONLY: eigvz
USE is_hermitian_mod, ONLY: is_hermitian
use, intrinsic :: ieee_arithmetic ! needed for IEEE_IS_FINITE
IMPLICIT NONE

!! Variables in/out

INTEGER, intent(in) :: &
    nb, & ! number of bands for a current k-point
    idg1, idg2 ! intext of the first and last degenerate state in the block
COMPLEX(kind=sp), intent(in) :: &
    pijA(:,:), & ! momentum matrix elements [at.u.]
    pijB(:,:)
REAL(kind=sp), intent(in) :: &
    dEij(:,:) ! energy differences E_n-E_i [Ha]
REAL(kind=dp), ALLOCATABLE, intent(out) :: &
    oam(:)  ! Berry curvature for a block of degenerate bands 
            ! (allocated inside CALL eigvd(...) )

!! Variables internal

COMPLEX(kind=dp), ALLOCATABLE :: &
    M(:,:), & ! Matrix similar to Eq. (6) in mstar paper (https://doi.org/10.1016/j.cpc.2020.107648)
    Mcorr(:,:) ! intermediates for Kahan summation
COMPLEX(kind=dp) :: &
    Lnln, & ! component leading to M(m,n)
    temp, corrected_term ! intermediates for Kahan summation
COMPLEX(kind=sp) :: &
    p2 ! antisymmetrized product of momentum matrix elements
REAL(kind=dp) :: &
    dEi, dEj, & ! E_n-E_i and E_n-E_j [Ha]
    inverse_gap_sym ! symmetrized reciprocal energy denominator [Ha^-1]
INTEGER :: &
    n, & ! band indices
    i, j, & ! counter
    ndg ! number of degenerate bands

!! Parameters

CHARACTER(len=256), PARAMETER :: &
    wformat_m = '(*( :, "(", ES12.4E3, ",", ES12.4E3, ")", : , 1X ))'

ndg = 1+idg2-idg1 ! number of degenerate bands
ALLOCATE( M(ndg,ndg), Mcorr(ndg,ndg) )

!! Check dimensions of INTENT(IN) arrays

! Check dimensions of pijA
IF (SIZE(pijA, 1) /= ndg .OR. SIZE(pijA, 2) /= nb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: pijA shape is (', SIZE(pijA,1), ',', SIZE(pijA,2), &
        ') but expected (', ndg, ',', nb, ')'
    ERROR STOP 'Inconsistent pijA dimensions'
END IF
! Check dimensions of pijB
IF (SIZE(pijB, 1) /= nb .OR. SIZE(pijB, 2) /= ndg) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: pijB shape is (', SIZE(pijB,1), ',', SIZE(pijB,2), &
        ') but expected (', nb, ',', ndg, ')'
    ERROR STOP 'Inconsistent pijB dimensions'
END IF
! Check dimensions of dEij
IF (SIZE(dEij, 1) /= ndg .OR. SIZE(dEij, 2) /= nb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: dEij shape is (', SIZE(dEij,1), ',', SIZE(dEij,2), &
        ') but expected (', ndg, ',', nb, ')'
    ERROR STOP 'Inconsistent dEij dimensions'
END IF

!! Construct the effective Hermitian OAM matrix
!
! For a degenerate block, the matrix element is
!
!   M_ij = (i/2) SUM_n P_ijn *
!          [1/(E_n-E_i) + 1/(E_n-E_j)],
!
! where
!
!   P_ijn = A_i,n B_n,j - B_i,n A_n,j.
!
! Intermediate states n inside the current degenerate block are excluded.
! The input convention is dEij(i,n) = E_n-E_i.
!
! For i = j, the expression reduces to
!
!   M_ii = i SUM_n P_iin/(E_n-E_i)
!        = 2 SUM_n Im[A_i,n B_n,i]/(E_i-E_n),
!
! which is exactly the validated non-degenerate OAM formula. For E_i=E_j,
! the expression is also identical to the previous degenerate formula. The
! difference appears only for numerically near-degenerate external states.

DO i = 1, ndg
    DO j = i, ndg
        M(i,j) = (0.0_dp, 0.0_dp) ! initialize
        Mcorr(i,j) = (0.0_dp, 0.0_dp)
        DO n = 1, nb
            IF (n < idg1 .OR. n > idg2) THEN ! ignore degenerate bands
                ! note that pijB(i,n) = CONJG(pijB(n,i))
                ! Products are formed in sp, then promoted to dp before
                ! division and Kahan compensated summation.
                p2 = pijA(i,n)*pijB(n,j) - &
                    CONJG(pijB(n,i))*CONJG(pijA(j,n))

                dEi = REAL(dEij(i,n), KIND=dp)
                dEj = REAL(dEij(j,n), KIND=dp)

                ! Preserve the validated isolated-band operation directly.
                ! For i /= j, use the Hermitian average of reciprocal gaps.
                IF (i == j) THEN
                    inverse_gap_sym = 1.0_dp/dEi
                ELSE
                    inverse_gap_sym = 0.5_dp*(1.0_dp/dEi + 1.0_dp/dEj)
                END IF

                Lnln = (0.0_dp, 1.0_dp) * CMPLX(p2, KIND=dp) * &
                    CMPLX(inverse_gap_sym, 0.0_dp, KIND=dp)

                ! make sure Lnln is finite (not NaN and not Inf)
                IF (.NOT. IEEE_IS_FINITE(AIMAG(Lnln)) &
                    .OR. .NOT. IEEE_IS_FINITE(REAL(Lnln,dp))) THEN
                    WRITE(*,*) 'i =', i
                    WRITE(*,*) 'j =', j
                    WRITE(*,*) 'n =', n
                    WRITE(*,*) 'pijA(i,n) =', pijA(i,n)
                    WRITE(*,*) 'pijB(n,j) =', pijB(n,j)
                    WRITE(*,*) 'dEij(i,n) =', dEij(i,n)
                    WRITE(*,*) 'dEij(j,n) =', dEij(j,n)
                    WRITE(*,*) 'inverse_gap_sym = ', inverse_gap_sym
                    WRITE(*,*) 'Lnln = ', Lnln
                    WRITE(*,*) 'p2 = ', p2
                    ERROR STOP 'Error: Lnln is not finite'
                END IF
                ! Kahan summation into M(i,j)
                corrected_term = Lnln - Mcorr(i,j)
                temp = M(i,j) + corrected_term
                Mcorr(i,j) = (temp - M(i,j)) - corrected_term
                M(i,j) = temp
            END IF
        END DO ! n
        IF (i /= j) M(j,i) = CONJG(M(i,j)) ! Hermitian M
    END DO ! j
END DO ! i

!! OAM is eigenvalues of the M matrix

IF (ndg > 1) THEN ! more than 1 degenerate band, M is a matrix
    ! M is Hermitian and can be diagonalized with eigvz.
    IF (.not. is_hermitian(M)) THEN
        WRITE(*,'(A)') 'M='
        DO i = 1, ndg
            WRITE(*,wformat_m) (M(i, j), j=1,ndg)
        END DO
        ERROR STOP 'in degenoam: Matrix "M" is not Hermitian'
    END IF
    CALL eigvz(ndg, M, & ! <- args in 
        oam) ! -> args out (allocated inside)
ELSE ! band is not degenerate, M is a complex number (not matrix)
    ALLOCATE( oam(ndg) )
    oam(1) = REAL(M(1,1))
END IF

DEALLOCATE( M, Mcorr )
RETURN
END SUBROUTINE degenoam

END MODULE degenoam_mod
