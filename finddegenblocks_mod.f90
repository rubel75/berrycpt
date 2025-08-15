MODULE finddegenblocks_mod
CONTAINS

SUBROUTINE finddegenblocks(nb, dEij, etol, dg_group, ngroups)
! Identify degenerate groups of bands based on energy differences.
!
! Given a symmetric matrix dEij(i,j) of energy differences between bands,
! this subroutine groups together bands whose pairwise energy differences
! are less than the specified tolerance etol, indicating degeneracy.
! Each group is assigned a unique integer ID.
!
! Inputs:
!   - dEij(nb, nb) : matrix of band energy differences (ε_i - ε_j)
!   - nb           : number of bands
!   - etol         : tolerance for considering two bands degenerate
!
! Outputs:
!   - dg_group(nb) : array assigning each band a group ID (1 to ngroups)
!   - ngroups      : total number of degenerate groups found

USE precision_mod, ONLY: sp
IMPLICIT NONE
INTEGER, INTENT(IN) :: nb
REAL(KIND=sp), INTENT(IN)  :: dEij(:, :) ! size nb x nb
REAL(KIND=sp), INTENT(IN)  :: etol
INTEGER, allocatable, INTENT(OUT):: dg_group(:)
INTEGER, INTENT(OUT):: ngroups

LOGICAL, allocatable :: visited(:)
INTEGER, allocatable :: stack(:)
INTEGER :: top, i, j, k, gid

allocate( dg_group(nb) ) ! allocate INTENT(OUT) array, but not deallocate here
allocate ( visited(nb), stack(nb) )

!! Check dimensions of INTENT(IN) arrays

IF (SIZE(dEij,1) /= nb .OR. SIZE(dEij,2) /= nb) THEN
    WRITE(*,'(A,I0,A,I0,A,I0,A)') &
        'Error: dEij shape is (', SIZE(dEij,1), ',', SIZE(dEij,2), &
        ') but expected square (', nb, ')'
    ERROR STOP 'Inconsistent dEij dimensions'
END IF

!! Create groups

visited = .FALSE.
dg_group = 0
gid = 0
DO i = 1, nb
    IF (.NOT. visited(i)) THEN
        gid = gid + 1
        top = 1
        stack(top) = i
        visited(i) = .TRUE.
        dg_group(i) = gid
        DO WHILE (top > 0)
            j = stack(top)
            top = top - 1
            DO k = 1, nb
                IF (k /= j .AND. ABS(dEij(j,k)) < etol .AND. .NOT. visited(k)) THEN
                    top = top + 1
                    stack(top) = k
                    visited(k) = .TRUE.
                    dg_group(k) = gid
                END IF
            END DO
        END DO
    END IF
END DO

ngroups = gid

deallocate ( visited, stack )
RETURN
END SUBROUTINE finddegenblocks

END MODULE finddegenblocks_mod