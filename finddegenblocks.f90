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

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nb
    REAL*4, INTENT(IN)  :: dEij(nb, nb)
    REAL*4, INTENT(IN)  :: etol
    INTEGER, INTENT(OUT):: dg_group(nb)
    INTEGER, INTENT(OUT):: ngroups

    LOGICAL :: visited(nb)
    INTEGER :: stack(nb)
    INTEGER :: top, i, j, k, gid

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
END SUBROUTINE finddegenblocks