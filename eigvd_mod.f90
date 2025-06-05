MODULE eigvd_mod
CONTAINS

SUBROUTINE eigvd(n, H, & ! <- args in 
            EIGV) ! -> args out

! Solve a _real_ matrix eigenvalue problem double precision

USE precision_mod, ONLY: sp, dp
implicit none

!! Variables in-out

INTEGER, intent(in) :: &
    n ! size of H matrixs
REAL(kind=dp), intent(in) :: &
    H(:,:) ! square matrix (size n x n)
REAL(kind=dp), allocatable, intent(out) :: &
    EIGV(:) ! eigenvalues

!! Internal variables

INTEGER :: &
    lwork, liwork, & ! The size of the work arrays
    info ! status argument
INTEGER, ALLOCATABLE :: iwork(:) ! The size of the work array (lwork>=n)
REAL(kind=dp), ALLOCATABLE :: &
    A(:,:), & ! used for H matrix and then overwritten by eigenvectors
    work(:) ! workspace array

!! Parameters

CHARACTER(len=1), PARAMETER :: &
    uplo='U', & ! stores the upper triangular part of A
    jobz = 'V' ! compute eigenvectors too

!! External subroutines

EXTERNAL :: &
    dsyevd ! Intel MKL

ALLOCATE( EIGV(n) ) ! allocate intent(out) array, but not deallocate here
ALLOCATE( A(n,n) )

!! Reassign H to another variable to avoid it being overwritten

A = H

!! Determine size of work arrays

lwork = -1
liwork = -1
allocate( work(1), iwork(1) )
call dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! double precision (kind=8)
!call ssyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! single precision (kind=4)
IF (info .ne. 0) THEN
    write (*,*) "ERROR in ssyevd (1st call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in "//&
        "rows and columns info/(n+1) through mod(info,n+1)."
    error stop
END IF
lwork = INT( work(1) )
liwork = iwork(1)
deallocate ( work, iwork )

!! Solve eigenvalues problem

allocate( work(lwork), iwork(liwork) )
call dsyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! double precision (kind=8)
!call ssyevd(jobz, uplo, n, A, n, EIGV, work, lwork, iwork, liwork, info) ! single precision (kind=4)
deallocate ( work, iwork, A )
IF (info .ne. 0) THEN
    write (*,*) "ERROR in ssyevd (2nd call): info = ", info
    write (*,*) "If info = -i, the i-th parameter had an illegal value."
    write (*,*) "If info = i, and jobz = 'V', then the algorithm failed to "//&
        "compute an eigenvalue while working on the submatrix lying in rows "//&
        "and columns info/(n+1) through mod(info,n+1)."
    error stop
END IF
! eigenvectors are returned through A-matrix are not used

RETURN
END SUBROUTINE eigvd

END MODULE eigvd_mod