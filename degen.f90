SUBROUTINE degen(nb, idg1, idg2, pijA, pijB, dEij, & ! <- args in 
    bcurv) ! -> args out

! Berry curvature for a block of degenerate bands

!! Variables in/out

implicit none
INTEGER, intent(in) :: &
    nb, & ! number of bands for a current k-point
    idg1, idg2 ! intext of the first and last degenerate state in the block
COMPLEX(kind=4), intent(in) :: &
    pijA(1+idg2-idg1,nb), & ! momentum matrix elements [at.u.]
    pijB(nb,1+idg2-idg1)
REAL(kind=4), intent(in) :: &
    dEij(1+idg2-idg1,nb) ! energy differences E_i-E_j [Ha]
REAL(kind=4), intent(out) :: &
    bcurv(1+idg2-idg1) ! Berry curvature for a block of degenerate bands

!! Variables internal

REAL(kind=4) :: &
    M(1+idg2-idg1,1+idg2-idg1), & ! Matrix similar to Eq. (6) in mstar paper (https://doi.org/10.1016/j.cpc.2020.107648)
    p2, & ! product of momentum matrix elements
    dE, & ! energy difference [Ha]
    omega ! component leading to M(m,n)
INTEGER :: &
    n, idg, & ! band indices
    i, j, & ! counter
    ndg ! number of degenerate bands


!! construct M matrix

ndg = 1+idg2-idg1 ! number of degenerate bands
DO i = 1, ndg
    DO j = 1, ndg
        idg = idg1 + i - 1
        M(i,j) = 0.0 ! initialize
        DO n = 1, nb
            ! extract complex part of the product
            ! take into account that i*5i = -5
            IF (n < idg1 .or. n > idg2) THEN ! ignore degenerate bands
                p2 = (-2.0) * AIMAG( pijA(i,n)*pijB(n,j))
                dE = (dEij(i,n) + dEij(j,n))/2.0
                omega = p2/(dE*dE)
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
                M(i,j) = M(i,j) + omega
            END IF
        END DO ! n
    END DO ! j
END DO ! i

!! Berry corvature is eigenvalues of the M matrix

CALL eigvs(ndg, M, & ! <- args in 
    bcurv) ! -> args out

RETURN
END SUBROUTINE degen
