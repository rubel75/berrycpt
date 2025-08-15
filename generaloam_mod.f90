MODULE generaloam_mod
CONTAINS

SUBROUTINE generaloam(nb, nvb, idg1, idg2, pijA, pijB, dEij, & ! <- args in 
    oamgnrlzdg) ! -> args out

! Orbital angular momentum (OAM) for a block of degenerate bands

USE precision_mod, ONLY: sp, dp
USE eigvz_mod, ONLY: eigvz
use, intrinsic :: ieee_arithmetic ! needed for IEEE_IS_FINITE
IMPLICIT NONE

!! Variables in/out

INTEGER, intent(in) :: &
    nb, & ! number of bands for a current k-point
    nvb, & ! number of valence bands for a current k-point
    idg1, idg2 ! intext of the first and last degenerate state in the block
COMPLEX(kind=sp), intent(in) :: &
    pijA(:,:), & ! momentum matrix elements [at.u.]
    pijB(:,:)
REAL(kind=sp), intent(in) :: &
    dEij(:,:) ! energy differences E_i-E_j [Ha]
COMPLEX(kind=dp), ALLOCATABLE, intent(out) :: &
    oamgnrlzdg(:,:)  ! Berry curvature for a block of degenerate bands 
            ! (allocated inside CALL eigvd(...) )

!! Variables internal

COMPLEX(kind=dp), ALLOCATABLE :: &
    M(:,:), & ! Matrix similar to Eq. (6) in mstar paper (https://doi.org/10.1016/j.cpc.2020.107648)
    Mcorr(:,:) ! intermediates for Kahan summation
COMPLEX(kind=dp) :: &
    Lnlm, & ! component leading to M(m,n)
    temp, corrected_term, & ! intermediates for Kahan summation    
    p2 ! product of momentum matrix elements
INTEGER :: &
    n, & ! band indices
    i, j, & ! counter
    ndg ! number of degenerate bands

ndg = 1+idg2-idg1 ! number of degenerate bands
ALLOCATE( M(nvb,nvb), Mcorr(nvb,nvb), oamgnrlzdg(nvb,nvb) )

!! Check dimensions of INTENT(IN) arrays

! Check dimensions of pijA
IF (SIZE(pijA, 1) /= nvb .OR. SIZE(pijA, 2) /= nb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: pijA shape is (', SIZE(pijA,1), ',', SIZE(pijA,2), &
        ') but expected (', nvb, ',', nb, ')'
    ERROR STOP 'Inconsistent pijA dimensions'
END IF
! Check dimensions of pijB
IF (SIZE(pijB, 1) /= nb .OR. SIZE(pijB, 2) /= nvb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: pijB shape is (', SIZE(pijB,1), ',', SIZE(pijB,2), &
        ') but expected (', nb, ',', nvb, ')'
    ERROR STOP 'Inconsistent pijB dimensions'
END IF
! Check dimensions of dEij
IF (SIZE(dEij, 1) /= nb .OR. SIZE(dEij, 2) /= nb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A)') &
        'Error: dEij shape is (', SIZE(dEij,1), ',', SIZE(dEij,2), &
        ') but expected (', nb, ',', nb, ')'
    ERROR STOP 'Inconsistent dEij dimensions'
END IF

!! construct M matrix
    
DO i = idg1, idg2
    DO j = i, nvb
        M(i,j) = (0.0_dp, 0.0_dp)! initialize
        Mcorr(i,j) = (0.0_dp, 0.0_dp)
        DO n = 1, nb
            ! The summation done in 2 terms
            ! Term 1:
            ! L1_{i,j} = (i*hbar/2*m0)*SUM_{n /= i, with i \in D} (px_{i,n}*py_{n,j} - py_{i,n}*px_{n,j})/(E_{i} - E_{n})
            ! here D is a subset of degenerate bands in the range [idg1, idg2]
            IF ( n < idg1 .or. n > idg2 ) THEN ! ignore degenerate bands (n={D})
                ! note that pijB(i,n) = CONJG(pijB(n,i))
                p2 = pijA(i,n)*pijB(n,j) - CONJG(pijB(n,i))*CONJG(pijA(j,n)) ! single precision
                ! double precision
                Lnlm = CMPLX(p2, kind=dp)/REAL(dEij(i,n), kind=dp)
                ! make sure Lnlm is finite (not NaN and not Inf)
                IF (.not. IEEE_IS_FINITE(AIMAG(Lnlm)) &
                    .or. .not. IEEE_IS_FINITE(REAL(Lnlm,dp))) THEN
                    WRITE(*,*) 'i =', i
                    WRITE(*,*) 'j =', j
                    WRITE(*,*) 'n =', n
                    WRITE(*,*) 'pijA(i,n) =', pijA(i,n)
                    WRITE(*,*) 'pijB(n,j) =', pijB(n,j)
                    WRITE(*,*) 'dEij(i,n) =', dEij(i,n)
                    WRITE(*,*) 'Lnlm = ', Lnlm
                    WRITE(*,*) 'p2 = ', p2
                    ERROR STOP 'Error: Lnlm is not finite'
                END IF
                ! Kahan summation into M(i,j)
                corrected_term = Lnlm - Mcorr(i,j)
                temp = M(i,j) + corrected_term
                Mcorr(i,j) = (temp - M(i,j)) - corrected_term
                M(i,j) = temp
            END IF ! term 1
            ! Term 2:
            ! L2_{i,j} = (i*hbar/2*m0)*SUM_{n /= j, with i \in D} (px_{i,n}*py_{n,j} - py_{i,n}*px_{n,j})/(E_{i} - E_{n})
            ! here D is a subset of degenerate bands in the range [idg1, idg2]
            IF ( n /= j ) THEN ! ignore degenerate bands (n=j)
                ! ignore degenerate bands (j={D} and n={D} at the same time)
                IF ( (j >= idg1 .or. j <= idg2) .and. (n >= idg1 .or. n <= idg2) ) THEN 
                    CYCLE ! go to next 'n' and skip the term
                END IF
                ! note that pijB(i,n) = CONJG(pijB(n,i))
                p2 = pijA(i,n)*pijB(n,j) - CONJG(pijB(n,i))*CONJG(pijA(j,n)) ! single precision
                ! double precision
                Lnlm = CMPLX(p2, kind=dp)/REAL(dEij(j,n), kind=dp)
                ! make sure Lnlm is finite (not NaN and not Inf)
                IF (.not. IEEE_IS_FINITE(AIMAG(Lnlm)) &
                    .or. .not. IEEE_IS_FINITE(REAL(Lnlm,dp))) THEN
                    WRITE(*,*) 'i =', i
                    WRITE(*,*) 'j =', j
                    WRITE(*,*) 'n =', n
                    WRITE(*,*) 'pijA(i,n) =', pijA(i,n)
                    WRITE(*,*) 'pijB(n,j) =', pijB(n,j)
                    WRITE(*,*) 'dEij(i,n) =', dEij(i,n)
                    WRITE(*,*) 'Lnlm = ', Lnlm
                    WRITE(*,*) 'p2 = ', p2
                    ERROR STOP 'Error: Lnlm is not finite'
                END IF
                ! Kahan summation into M(i,j)
                corrected_term = Lnlm - Mcorr(i,j)
                temp = M(i,j) + corrected_term
                Mcorr(i,j) = (temp - M(i,j)) - corrected_term
                M(i,j) = temp
            END IF ! term 2
        END DO ! n
        IF (i /= j) M(j,i) = CONJG(M(i,j)) ! symmetrize Hermitian M
    END DO ! j
END DO ! i

! Preparing the final result to be returned

oamgnrlzdg = ((0.0_dp, 1.0_dp)/2.0_dp) * M ! L = (i*hbar/2*m0)*M
! Because we work in a.u. hbar = 1, m0 = 1

DEALLOCATE( M, Mcorr )
RETURN
END SUBROUTINE generaloam

END MODULE generaloam_mod
