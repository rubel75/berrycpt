MODULE calculate_bcurv_kpoint_mod

USE precision_mod, ONLY: sp, dp
USE degenbc_mod, ONLY: degenbc

IMPLICIT NONE
PRIVATE

PUBLIC :: calculate_bcurv_kpoint

CONTAINS

SUBROUTINE calculate_bcurv_kpoint( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij, &         ! in
    bcurv)                      ! out

! Calculate the three pseudovector components of the Berry curvature for
! every target band at one k-point and one spin channel.
!
! Array and energy-difference conventions
! ---------------------------------------
! pijA(alpha,i,j) = <i|A_alpha|j> in atomic units
! pijB(beta,i,j)  = <i|B_beta |j> in atomic units
! dEij(i,j)       = E_j - E_i in Hartree
!
! Both operator arrays use rectangular storage:
!     (Cartesian component, target band, all bands)
! where the target-band index runs only through nbcderki and the transition
! band index runs through nbki. The reverse matrix elements required by the
! perturbation-theory expression are reconstructed using Hermiticity:
!     <l|B_beta|j> = CONJG(<j|B_beta|l>).
!
! Pseudovector convention
! -----------------------
! bcurv(:,1) = Omega_yz
! bcurv(:,2) = Omega_zx
! bcurv(:,3) = Omega_xy
!
! Exact non-degenerate expression used
! ------------------------------------
! In atomic units, where hbar = 1, the result for an isolated target band n
! is
!
!   Omega_n^(A,B)(alpha,beta)
!       = -2 Im SUM_{l /= n}
!           A_alpha(n,l) B_beta(l,n) / (E_n - E_l)^2.
!
! Since dEij(n,l) = E_l - E_n, the squared denominator is exactly
! dEij(n,l)^2. For ordinary Berry curvature A = B = velocity. For the
! sigma_z-normalized spin Berry curvature, A is the spin-current operator
! and B is the velocity operator.
!
! Exact degenerate expression used
! --------------------------------
! For a contiguous degenerate block D = {idg1,...,idg2}, degenbc constructs
! the non-Abelian Hermitian Berry-curvature matrix
!
!   M_ij^(A,B) = i SUM_{l not in D}
!       [ A_alpha(i,l) B_beta(l,j)
!       - B_beta(i,l) A_alpha(l,j) ]
!       / [ (E_l-E_i)(E_l-E_j) ].
!
! In the stored convention, the denominator is
!
!   dEij(i,l) * dEij(j,l).
!
! For i = j, this reduces exactly to the validated squared denominator of the
! non-degenerate formula. For E_i = E_j, it also reduces to the earlier
! degenerate expression based on the square of one common gap. The difference
! appears only when states grouped as degenerate have slightly different
! numerical energies. The Berry curvatures assigned to the block are the
! eigenvalues of M. The established degenbc -> eigvz -> ZHEEVD path, including
! all MKL-facing declarations and workspace handling, is unchanged.
!
! Precision and summation
! -----------------------
! The input momentum/velocity matrices and energy differences are sp, matching
! the matrix-element readers. Group-specific arrays remain sp. The two products
! in the antisymmetrized numerator are formed in sp. Each energy difference is
! then promoted separately to dp, and their product is evaluated in dp. Each
! effective-matrix term, complex Kahan compensated summation, the effective
! matrix, its eigenvalues, and the final bcurv array are dp. No dp result is
! converted back to sp.

! Current k-point dimensions and degenerate-group assignment
INTEGER, INTENT(IN) :: &
    nbki, &      ! number of bands at the current k-point and spin
    nbcderki     ! number of target bands at the current k-point and spin
INTEGER, INTENT(IN) :: &
    dg_group(:)  ! contiguous degenerate-group label for each band

! Rectangular operator matrices and energy differences
COMPLEX(KIND=sp), INTENT(IN) :: &
    pijA(:,:,:), & ! first operator matrix elements [a.u.]
    pijB(:,:,:)    ! second operator matrix elements [a.u.]
REAL(KIND=sp), INTENT(IN) :: &
    dEij(:,:)      ! E_j-E_i for target band i and transition band j [Ha]

! Berry-curvature pseudovector components
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: &
    bcurv(:,:)     ! (target band, yz/zx/xy) [bohr^2]

! Group-specific matrices passed to the established degenerate-PT routine
COMPLEX(KIND=sp), ALLOCATABLE :: &
    pijAdg(:,:), & ! <i|A_alpha|l> for target states i in the current group
    pijBdg(:,:)    ! <l|B_beta|j> for target states j in the current group
REAL(KIND=sp), ALLOCATABLE :: &
    dEijdg(:,:)    ! E_l-E_i for target states i in the current group [Ha]
REAL(KIND=dp), ALLOCATABLE :: &
    bcurvdg(:)     ! eigenvalues of the effective matrix for one group

! Loop indices and current-group bounds
INTEGER :: &
    icomp, & ! pseudovector component: 1=yz, 2=zx, 3=xy
    alpha, & ! first Cartesian operator direction
    beta, &  ! second Cartesian operator direction
    idg1, &  ! first band in the current degenerate group
    idg2, &  ! last band in the current degenerate group
    iband, & ! absolute target-band index
    m, &     ! local target-band index inside the current group
    ndg, &   ! number of bands in the current degenerate group
    ierr     ! allocation status

! Allocation diagnostic
CHARACTER(LEN=512) :: &
    iomsg

! Cartesian pairs corresponding to {Omega_yz, Omega_zx, Omega_xy}
INTEGER, PARAMETER :: &
    alpha_component(3) = [2, 3, 1], &
    beta_component(3)  = [3, 1, 2]

CALL validate_bcurv_inputs( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij)           ! in

ALLOCATE( &
    bcurv(nbcderki,3), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to allocate Berry-curvature array'
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF
bcurv = 0.0_dp

! Degenerate groups are contiguous and numbered consecutively. The cutoff
! nbcderki is required to lie at a group boundary, so every group encountered
! here is complete.
idg1 = 1
DO WHILE (idg1 <= nbcderki)

    idg2 = idg1
    DO WHILE (idg2 < nbcderki)
        IF (dg_group(idg2+1) /= dg_group(idg1)) EXIT
        idg2 = idg2 + 1
    END DO
    ndg = idg2 - idg1 + 1

    ALLOCATE( &
        pijAdg(ndg,nbki), &
        pijBdg(nbki,ndg), &
        dEijdg(ndg,nbki), &
        STAT=ierr, ERRMSG=iomsg)
    IF (ierr /= 0) THEN
        WRITE(*,'(/,A)') &
            ' ERROR: Unable to allocate degenerate-group Berry arrays'
        WRITE(*,'(A,I0,A,I0)') &
            ' Group band range = ', idg1, ':', idg2
        WRITE(*,'(A)') ' '//TRIM(iomsg)
        ERROR STOP 1
    END IF

    DO icomp = 1, 3
        alpha = alpha_component(icomp)
        beta = beta_component(icomp)

        ! Form the two operator orientations required by degenbc. The first
        ! operator is available directly as <i|A_alpha|l>. Rectangular storage
        ! contains <j|B_beta|l>, so Hermiticity gives
        ! <l|B_beta|j> = CONJG(<j|B_beta|l>).
        DO m = 1, ndg
            iband = idg1 + m - 1
            pijAdg(m,:) = pijA(alpha,iband,1:nbki)
            pijBdg(:,m) = CONJG(pijB(beta,iband,1:nbki))
            dEijdg(m,:) = dEij(iband,1:nbki)
        END DO

        CALL degenbc( &
            nbki, idg1, idg2, &          ! in
            pijAdg, pijBdg, dEijdg, &    ! in
            bcurvdg)                     ! out

        IF (.NOT. ALLOCATED(bcurvdg)) THEN
            WRITE(*,'(/,A)') &
                ' ERROR: degenbc did not allocate its output array'
            ERROR STOP 1
        END IF

        IF (SIZE(bcurvdg) /= ndg) THEN
            WRITE(*,'(/,A)') &
                ' ERROR: Unexpected Berry-curvature block size from degenbc'
            WRITE(*,'(A,I0)') ' Returned size = ', SIZE(bcurvdg)
            WRITE(*,'(A,I0)') ' Expected size = ', ndg
            ERROR STOP 1
        END IF

        bcurv(idg1:idg2,icomp) = bcurvdg
        DEALLOCATE(bcurvdg)
    END DO

    DEALLOCATE(pijAdg, pijBdg, dEijdg)
    idg1 = idg2 + 1
END DO

END SUBROUTINE calculate_bcurv_kpoint


SUBROUTINE validate_bcurv_inputs( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij)           ! in

! Validate only the interface invariants required by this module.

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
    WRITE(*,'(/,A)') ' ERROR: Invalid current derivative-band limit'
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

END SUBROUTINE validate_bcurv_inputs

END MODULE calculate_bcurv_kpoint_mod
