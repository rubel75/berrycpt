MODULE calculate_oam_kpoint_mod

USE precision_mod, ONLY: sp, dp
USE degenoam_mod, ONLY: degenoam

IMPLICIT NONE
PRIVATE

PUBLIC :: calculate_oam_kpoint

CONTAINS

SUBROUTINE calculate_oam_kpoint( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij, &         ! in
    oam)                        ! out

! Calculate the three pseudovector components of the orbital angular momentum
! for every target band at one k-point and one spin channel.
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
! oam(:,1) = L_yz
! oam(:,2) = L_zx
! oam(:,3) = L_xy
!
! In Cartesian-vector notation, these are the x, y, and z components,
! respectively. The result is expressed in atomic units, where hbar = 1,
! and can be interpreted as orbital angular momentum in units of hbar.
!
! Exact non-degenerate expression used
! ------------------------------------
! For an isolated target band n, the established degenoam routine evaluates
!
!   L_n^(A,B)(alpha,beta)
!       = 2 Im SUM_{l /= n}
!           A_alpha(n,l) B_beta(l,n) / (E_n - E_l).
!
! Since dEij(n,l) = E_l - E_n, the identical expression used by the code is
!
!   L_n^(A,B)(alpha,beta)
!       = -2 Im SUM_{l /= n}
!           A_alpha(n,l) B_beta(l,n) / dEij(n,l).
!
! This follows directly from the diagonal element of the effective matrix
! constructed below. For ordinary orbital angular momentum, A = B = p.
!
! The spin-projected OAM used elsewhere in berrycpt is not constructed like
! the spin Berry curvature. Its spin-up contribution uses
!
!   A = p^up,  B = p^up,
!
! and a possible spin-down contribution analogously uses
!
!   A = p^dn,  B = p^dn.
!
! The mixed choices (p^up,p), (p^dn,p), and (p^up-p^dn,p) belong to the
! spin-resolved Berry-curvature definitions and are not the spin-projected
! OAM definition implemented by the existing code.
!
! Exact degenerate expression used
! --------------------------------
! For a contiguous degenerate block D = {idg1,...,idg2}, degenoam constructs
! the Hermitian effective OAM matrix
!
!   M_ij^(A,B) = (i/2) SUM_{l not in D}
!       [ A_alpha(i,l) B_beta(l,j)
!       - B_beta(i,l) A_alpha(l,j) ]
!       * [ 1/(E_l-E_i) + 1/(E_l-E_j) ].
!
! In the stored convention, the reciprocal-gap factor is
!
!   0.5 * [ 1/dEij(i,l) + 1/dEij(j,l) ].
!
! For i = j, this reduces exactly to the validated non-degenerate OAM formula.
! For E_i = E_j, it also reduces to the earlier degenerate expression based on
! one common gap. The difference appears only when states grouped as degenerate
! have slightly different numerical energies. The OAM values assigned to the
! block are the eigenvalues of M. The established degenoam -> eigvz -> ZHEEVD
! path, including all MKL-facing declarations and workspace handling, is
! unchanged.
!
! The numerator in degenoam is evaluated as
!
!   A_alpha(i,l) B_beta(l,j)
!       - CONJG(B_beta(l,i)) CONJG(A_alpha(j,l)),
!
! which is algebraically the same expression when A and B are Hermitian.
!
! Precision and compensated summation
! -----------------------------------
! The input momentum matrices and energy differences are sp, matching the
! matrix-element readers. Group-specific arrays remain sp. Each product in the
! antisymmetrized numerator is formed in sp. The two energy differences are
! promoted separately to dp before their reciprocals and average are evaluated.
! Each effective-matrix term, complex Kahan compensated summation, the
! Hermitian matrix, its eigenvalues, and the final oam array are dp. No dp
! result is converted back to sp.

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

! Orbital-angular-momentum pseudovector components
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: &
    oam(:,:)       ! (target band, yz/zx/xy) [hbar]

! Group-specific matrices passed to the established degenerate-PT routine
COMPLEX(KIND=sp), ALLOCATABLE :: &
    pijAdg(:,:), & ! <i|A_alpha|l> for target states i in the current group
    pijBdg(:,:)    ! <l|B_beta|j> for target states j in the current group
REAL(KIND=sp), ALLOCATABLE :: &
    dEijdg(:,:)    ! E_l-E_i for target states i in the current group [Ha]
REAL(KIND=dp), ALLOCATABLE :: &
    oamdg(:)       ! eigenvalues of the effective matrix for one group

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

! Cartesian pairs corresponding to {L_yz, L_zx, L_xy}
INTEGER, PARAMETER :: &
    alpha_component(3) = [2, 3, 1], &
    beta_component(3)  = [3, 1, 2]

CALL validate_oam_inputs( &
    nbki, nbcderki, dg_group, & ! in
    pijA, pijB, dEij)           ! in

ALLOCATE( &
    oam(nbcderki,3), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    WRITE(*,'(/,A)') ' ERROR: Unable to allocate orbital-angular-momentum array'
    WRITE(*,'(A)') ' '//TRIM(iomsg)
    ERROR STOP 1
END IF
oam = 0.0_dp

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
            ' ERROR: Unable to allocate degenerate-group OAM arrays'
        WRITE(*,'(A,I0,A,I0)') &
            ' Group band range = ', idg1, ':', idg2
        WRITE(*,'(A)') ' '//TRIM(iomsg)
        ERROR STOP 1
    END IF

    DO icomp = 1, 3
        alpha = alpha_component(icomp)
        beta = beta_component(icomp)

        ! Form the two operator orientations required by degenoam. The first
        ! operator is available directly as <i|A_alpha|l>. Rectangular storage
        ! contains <j|B_beta|l>, so Hermiticity gives
        ! <l|B_beta|j> = CONJG(<j|B_beta|l>).
        DO m = 1, ndg
            iband = idg1 + m - 1
            pijAdg(m,:) = pijA(alpha,iband,1:nbki)
            pijBdg(:,m) = CONJG(pijB(beta,iband,1:nbki))
            dEijdg(m,:) = dEij(iband,1:nbki)
        END DO

        CALL degenoam( &
            nbki, idg1, idg2, &          ! in
            pijAdg, pijBdg, dEijdg, &    ! in
            oamdg)                       ! out

        IF (.NOT. ALLOCATED(oamdg)) THEN
            WRITE(*,'(/,A)') &
                ' ERROR: degenoam did not allocate its output array'
            ERROR STOP 1
        END IF

        IF (SIZE(oamdg) /= ndg) THEN
            WRITE(*,'(/,A)') &
                ' ERROR: Unexpected OAM block size from degenoam'
            WRITE(*,'(A,I0)') ' Returned size = ', SIZE(oamdg)
            WRITE(*,'(A,I0)') ' Expected size = ', ndg
            ERROR STOP 1
        END IF

        oam(idg1:idg2,icomp) = oamdg
        DEALLOCATE(oamdg)
    END DO

    DEALLOCATE(pijAdg, pijBdg, dEijdg)
    idg1 = idg2 + 1
END DO

END SUBROUTINE calculate_oam_kpoint


SUBROUTINE validate_oam_inputs( &
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

END SUBROUTINE validate_oam_inputs

END MODULE calculate_oam_kpoint_mod
