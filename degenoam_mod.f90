MODULE degenoam_mod
CONTAINS

SUBROUTINE degenoam(nb, idg1, idg2, pijA, pijB, dEij, & ! <- args in 
    oam) ! -> args out

! Orbital angular momentum (OAM) for a block of degenerate bands

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
    oam(:)  ! Berry curvature for a block of degenerate bands 
            ! (allocated inside CALL eigvd(...) )

!! Variables internal

REAL(kind=dp) :: &
    Lnln, & ! component leading to M(m,n)
    temp, corrected_term ! intermediates for Kahan summation
REAL(kind=dp), ALLOCATABLE :: &
    M(:,:), & ! Matrix similar to Eq. (6) in mstar paper (https://doi.org/10.1016/j.cpc.2020.107648)
    Mcorr(:,:) ! intermediates for Kahan summation
REAL(kind=sp) :: &
    p2, & ! product of momentum matrix elements
    dE ! energy difference [Ha]
INTEGER :: &
    n, idg, & ! band indices
    i, j, & ! counter
    ndg ! number of degenerate bands

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

!! construct M matrix
    
DO i = 1, ndg
    DO j = i, ndg
        idg = idg1 + i - 1
        M(i,j) = 0.0 ! initialize
        DO n = 1, nb
            ! extract complex part of the product
            ! take into account that i*5i = -5
            IF (n < idg1 .or. n > idg2) THEN ! ignore degenerate bands
                p2 = (-2.0_sp) * AIMAG( pijA(i,n)*pijB(n,j)) ! single precision
                dE = (dEij(i,n) + dEij(j,n))/2.0_sp ! single precision
                ! double precision
                Lnln = REAL(-p2, dp)/REAL(dE, dp)
                ! make sure Lnln is finite (not NaN and not Inf)
                IF (Lnln /= Lnln .or. abs(Lnln) > HUGE(abs(Lnln))) THEN
                    WRITE(*,*) 'i =', i
                    WRITE(*,*) 'j =', j
                    WRITE(*,*) 'n =', n
                    WRITE(*,*) 'pijA(i,n) =', pijA(i,n)
                    WRITE(*,*) 'pijB(n,j) =', pijB(n,j)
                    WRITE(*,*) 'dEij(i,n) =', dEij(i,n)
                    WRITE(*,*) 'dEij(j,n) =', dEij(j,n)
                    WRITE(*,*) 'Lnln = ', Lnln
                    WRITE(*,*) 'dE = ', dE
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
        IF (i /= j) M(j,i) = M(i,j) ! symmetrize M
    END DO ! j
END DO ! i

!! OAM is eigenvalues of the M matrix

CALL eigvd(ndg, M, & ! <- args in 
    oam) ! -> args out (allocated inside)

DEALLOCATE( M, Mcorr )
RETURN
END SUBROUTINE degenoam

END MODULE degenoam_mod
