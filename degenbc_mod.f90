MODULE degenbc_mod
CONTAINS

SUBROUTINE degenbc(nb, idg1, idg2, pijA, pijB, dEij, & ! <- args in 
    bcurv) ! -> args out

! Berry curvature for a block of degenerate bands

USE precision_mod, ONLY: sp, dp
USE eigvd_mod, ONLY: eigvd
IMPLICIT NONE

!! Variables in/out

INTEGER, intent(in) :: &
    nb, & ! number of bands for a current k-point
    idg1, idg2 ! intext of the first and last degenerate state in the block
COMPLEX(kind=sp), intent(in) :: &
    pijA(:,:), & ! momentum matrix elements [at.u.]
    pijB(:,:)
REAL(kind=sp), intent(in) :: &
    dEij(:,:) ! energy differences E_i-E_j [Ha]
REAL(kind=dp), ALLOCATABLE, intent(out) :: &
    bcurv(:)    ! Berry curvature for a block of degenerate bands 
                ! (allocated inside CALL eigvd(...) )

!! Variables internal

REAL(kind=dp) :: &
    omega, & ! component leading to M(m,n)
    temp, corrected_term ! intermediates for Kahan summation
REAL(kind=dp), ALLOCATABLE :: &
    M(:,:), & ! Matrix similar to Eq. (6) in mstar paper (https://doi.org/10.1016/j.cpc.2020.107648)
    Mcorr(:,:) ! intermediates for Kahan summation
REAL(kind=sp) :: &
    p2, & ! product of momentum matrix elements
    dE ! energy difference [Ha]
INTEGER :: &
    n, & ! band index
    i, j, & ! counter
    ndg ! number of degenerate bands

!! construct M matrix

ndg = 1+idg2-idg1 ! number of degenerate bands
ALLOCATE( M(ndg,ndg), Mcorr(ndg,ndg) )
DO i = 1, ndg
    DO j = 1, ndg
        M(i,j) = 0.0_dp ! initialize
        Mcorr(i,j) = 0.0_dp
        DO n = 1, nb
            ! extract complex part of the product
            ! take into account that i*5i = -5
            IF (n < idg1 .or. n > idg2) THEN ! ignore degenerate bands
                p2 = (-2.0_sp) * AIMAG( pijA(i,n)*pijB(n,j)) ! single precision
                dE = (dEij(i,n) + dEij(j,n))/2.0_sp ! single precision
                ! double precision
                omega = REAL(p2, dp)/REAL(dE*dE, dp)
                ! make sure omega is finite (not NaN and not Inf)
                IF (omega /= omega .or. abs(omega) > HUGE(abs(omega))) THEN
                    WRITE(*,*) 'i =', i
                    WRITE(*,*) 'j =', j
                    WRITE(*,*) 'n =', n
                    WRITE(*,*) 'pijA(i,n) =', pijA(i,n)
                    WRITE(*,*) 'pijB(n,j) =', pijB(n,j)
                    WRITE(*,*) 'dEij(i,n) =', dEij(i,n)
                    WRITE(*,*) 'dEij(j,n) =', dEij(j,n)
                    WRITE(*,*) 'omega = ', omega
                    WRITE(*,*) 'dE = ', dE
                    WRITE(*,*) 'p2 = ', p2
                    ERROR STOP 'Error: omega is not finite'
                END IF
                ! Kahan summation into M(i,j)
                corrected_term = omega - Mcorr(i,j)
                temp = M(i,j) + corrected_term
                Mcorr(i,j) = (temp - M(i,j)) - corrected_term
                M(i,j) = temp
            END IF
        END DO ! n
    END DO ! j
END DO ! i

!! Berry curvature is given by the eigenvalues of the M matrix

! convert "M" to single precision to match eigvS
CALL eigvd(ndg, M, & ! <- args in 
    bcurv) ! -> args out (allocated inside)

DEALLOCATE( M, Mcorr )
RETURN
END SUBROUTINE degenbc

END MODULE degenbc_mod