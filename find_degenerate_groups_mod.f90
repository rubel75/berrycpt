MODULE find_degenerate_groups_mod

USE precision_mod, ONLY: sp

IMPLICIT NONE

PRIVATE
PUBLIC :: find_degenerate_groups


CONTAINS

SUBROUTINE find_degenerate_groups( &
        wien2k, nstot, nktot, nb, nbcder, nbk, ene, occup, &
        energy_tolerance, occupation_tolerance, &
        dg_groupk, ngroupsk, nbcderk)

! Calculation type and eigenvalue metadata
LOGICAL, INTENT(IN) :: wien2k
INTEGER, INTENT(IN) :: nstot, nktot, nb

! VASP derivative-band limit read from WAVEDER
! This value is ignored for WIEN2k.
INTEGER, INTENT(IN) :: nbcder

! Number of available eigenvalues at each k-point and spin
INTEGER, INTENT(IN) :: nbk(:,:)

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=sp), INTENT(IN) :: ene(:,:,:)

! Band-, k-point-, and spin-dependent occupations
REAL(KIND=sp), INTENT(IN) :: occup(:,:,:)

! Numerical tolerances for degeneracies [Ha] and nonzero occupations
REAL(KIND=sp), INTENT(IN) :: &
    energy_tolerance, occupation_tolerance

! Degenerate-group assignment for each band
! A value of zero marks a band that is absent at a given k-point and spin.
INTEGER, ALLOCATABLE, INTENT(OUT) :: dg_groupk(:,:,:)

! Number of degenerate groups at each k-point and spin
INTEGER, ALLOCATABLE, INTENT(OUT) :: ngroupsk(:,:)

! Derivative-band limit after k-point-dependent degeneracy adjustment
INTEGER, ALLOCATABLE, INTENT(OUT) :: nbcderk(:,:)

! Last occupied band at each k-point and spin
INTEGER :: nbocck(nktot,nstot)

! Common derivative-band limit and largest occupied-band index
INTEGER :: nbcder_common, nbocc_max

! Band, k-point, spin, and group counters
INTEGER :: &
    iband, ikpt, ispin, igroup, &
    candidate, first_band, cut_group

! Adjacent-band energy difference [Ha]
REAL(KIND=sp) :: energy_gap

IF (nstot <= 0 .OR. nktot <= 0 .OR. nb <= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid eigenvalue metadata were supplied'
    WRITE(*,'(3(A,I0))') &
        ' Number of spin channels = ', nstot, &
        ', number of k-points = ', nktot, &
        ', maximum number of bands = ', nb
    ERROR STOP 1
END IF

IF (energy_tolerance <= 0.0_sp) THEN
    WRITE(*,'(/,A)') ' ERROR: Degeneracy tolerance must be positive'
    WRITE(*,'(A,ES14.6E3)') &
        ' Degeneracy tolerance [Ha] = ', energy_tolerance
    ERROR STOP 1
END IF

IF (occupation_tolerance < 0.0_sp) THEN
    WRITE(*,'(/,A)') ' ERROR: Occupation tolerance must not be negative'
    WRITE(*,'(A,ES14.6E3)') &
        ' Occupation tolerance = ', occupation_tolerance
    ERROR STOP 1
END IF

IF (ANY(SHAPE(nbk) /= [nktot,nstot])) THEN
    WRITE(*,'(/,A)') ' ERROR: Inconsistent dimensions of nbk'
    WRITE(*,'(A,2(I0,1X))') ' Actual shape = ', SHAPE(nbk)
    WRITE(*,'(A,2(I0,1X))') ' Expected shape = ', nktot, nstot
    ERROR STOP 1
END IF

IF (ANY(SHAPE(ene) /= [nb,nktot,nstot])) THEN
    WRITE(*,'(/,A)') ' ERROR: Inconsistent dimensions of ene'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(ene)
    WRITE(*,'(A,3(I0,1X))') ' Expected shape = ', nb, nktot, nstot
    ERROR STOP 1
END IF

IF (ANY(SHAPE(occup) /= [nb,nktot,nstot])) THEN
    WRITE(*,'(/,A)') ' ERROR: Inconsistent dimensions of occup'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(occup)
    WRITE(*,'(A,3(I0,1X))') ' Expected shape = ', nb, nktot, nstot
    ERROR STOP 1
END IF

IF (ANY(nbk < 1) .OR. ANY(nbk > nb)) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid number of bands at a k-point and spin'
    WRITE(*,'(A,I0)') ' Minimum value in nbk = ', MINVAL(nbk)
    WRITE(*,'(A,I0)') ' Maximum value in nbk = ', MAXVAL(nbk)
    WRITE(*,'(A,I0)') ' Allowed maximum = ', nb
    ERROR STOP 1
END IF

ALLOCATE(dg_groupk(nb,nktot,nstot))
ALLOCATE(ngroupsk(nktot,nstot))
ALLOCATE(nbcderk(nktot,nstot))

dg_groupk = 0
ngroupsk = 0
nbocck = 0
nbcderk = 0

! Form chain-connected groups from adjacent ordered eigenvalues and
! determine the last band with a nonzero occupation at each k-point.
DO ispin = 1, nstot
    DO ikpt = 1, nktot
        igroup = 1
        dg_groupk(1,ikpt,ispin) = igroup

        IF (occup(1,ikpt,ispin) > occupation_tolerance) THEN
            nbocck(ikpt,ispin) = 1
        END IF

        DO iband = 2, nbk(ikpt,ispin)
            energy_gap = ene(iband,ikpt,ispin) - &
                ene(iband-1,ikpt,ispin)

            IF (energy_gap < -energy_tolerance) THEN
                WRITE(*,'(/,A)') &
                    ' ERROR: Eigenvalues are not ordered by band index'
                WRITE(*,'(A,I0)') ' Spin channel = ', ispin
                WRITE(*,'(A,I0)') ' K-point = ', ikpt
                WRITE(*,'(A,I0)') ' Lower band index = ', iband-1
                WRITE(*,'(A,ES14.6E3)') &
                    ' Lower-band energy [Ha] = ', &
                    ene(iband-1,ikpt,ispin)
                WRITE(*,'(A,I0)') ' Upper band index = ', iband
                WRITE(*,'(A,ES14.6E3)') &
                    ' Upper-band energy [Ha] = ', &
                    ene(iband,ikpt,ispin)
                ERROR STOP 1
            END IF

            IF (energy_gap > energy_tolerance) THEN
                igroup = igroup + 1
            END IF

            dg_groupk(iband,ikpt,ispin) = igroup

            IF (occup(iband,ikpt,ispin) > occupation_tolerance) THEN
                nbocck(ikpt,ispin) = iband
            END IF
        END DO

        ngroupsk(ikpt,ispin) = igroup

        IF (nbocck(ikpt,ispin) == 0) THEN
            WRITE(*,'(/,A)') ' ERROR: No occupied band was identified'
            WRITE(*,'(A,I0)') ' Spin channel = ', ispin
            WRITE(*,'(A,I0)') ' K-point = ', ikpt
            WRITE(*,'(A,ES14.6E3)') &
                ' Occupation tolerance = ', occupation_tolerance
            ERROR STOP 1
        END IF
    END DO
END DO

! Use one global occupied-band reference for WIEN2k so that the common
! derivative-band limit does not vary because of k-dependent occupations.
nbocc_max = MAXVAL(nbocck)

IF (wien2k) THEN
    nbcder_common = MIN(2 * nbocc_max, nb)
ELSE
    IF (nbcder <= 0 .OR. nbcder > nb) THEN
        WRITE(*,'(/,A)') &
            ' ERROR: Invalid NBANDS_CDER value was read from WAVEDER'
        WRITE(*,'(A,I0)') ' NBANDS_CDER = ', nbcder
        WRITE(*,'(A,I0)') ' Total number of bands = ', nb
        ERROR STOP 1
    END IF

    nbcder_common = nbcder
END IF

IF (nbocc_max > nbcder_common) THEN
    WRITE(*,'(/,A)') &
        ' ERROR: Derivative matrix elements do not include all occupied bands'
    WRITE(*,'(A,I0)') ' Highest occupied-band index = ', nbocc_max
    WRITE(*,'(A,I0)') &
        ' Common derivative-band limit = ', nbcder_common
    ERROR STOP 1
END IF

! Start every k-point and spin from the same common limit. Restrict it only
! by the number of available eigenvalues. If the boundary cuts a degenerate
! group, stop at the end of the preceding group.
DO ispin = 1, nstot
    DO ikpt = 1, nktot
        candidate = MIN(nbcder_common, nbk(ikpt,ispin))
        nbcderk(ikpt,ispin) = candidate

        IF (candidate < nbk(ikpt,ispin)) THEN
            IF (dg_groupk(candidate,ikpt,ispin) == &
                    dg_groupk(candidate+1,ikpt,ispin)) THEN

                cut_group = dg_groupk(candidate,ikpt,ispin)
                first_band = candidate

                DO WHILE (first_band > 1)
                    IF (dg_groupk(first_band-1,ikpt,ispin) /= &
                            cut_group) EXIT
                    first_band = first_band - 1
                END DO

                nbcderk(ikpt,ispin) = first_band - 1
            END IF
        END IF

        IF (nbcderk(ikpt,ispin) < nbocck(ikpt,ispin)) THEN
            WRITE(*,'(/,A)') &
                ' ERROR: Degeneracy adjustment excludes an occupied band'
            WRITE(*,'(A,I0)') ' Spin channel = ', ispin
            WRITE(*,'(A,I0)') ' K-point = ', ikpt
            WRITE(*,'(A,I0)') &
                ' Last occupied band = ', nbocck(ikpt,ispin)
            WRITE(*,'(A,I0)') &
                ' Common derivative-band limit = ', candidate
            WRITE(*,'(A,I0)') &
                ' Adjusted derivative-band limit = ', &
                nbcderk(ikpt,ispin)
            ERROR STOP 1
        END IF
    END DO
END DO

END SUBROUTINE find_degenerate_groups

END MODULE find_degenerate_groups_mod
