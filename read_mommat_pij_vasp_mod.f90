MODULE read_mommat_pij_vasp_mod

USE, INTRINSIC :: iso_fortran_env, ONLY: int32, int64, output_unit
USE precision_mod, ONLY: sp, dp

IMPLICIT NONE

PRIVATE
PUBLIC :: read_mommat_pij_vasp

CONTAINS

SUBROUTINE read_mommat_pij_vasp( &
        fnameinp, nstot, nktot, nbcder, nb, ene, &
        pijks, dEijks)

! Name of the VASP WAVEDER file
CHARACTER(LEN=*), INTENT(IN) :: fnameinp

! WAVEDER metadata
INTEGER, INTENT(IN) :: &
    nstot, &  ! total number of spin channels
    nktot, &  ! total number of k-points
    nbcder, & ! number of bands included in derivative overlaps
    nb         ! total number of bands

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=sp), INTENT(IN) :: ene(:,:,:)

! Rectangular momentum-matrix elements [a.u.]
! pijks(:,i,j,k,s) stores <i|p|j>, where i <= nbcder and j <= nb.
COMPLEX(KIND=sp), ALLOCATABLE, INTENT(OUT) :: pijks(:,:,:,:,:)

! Rectangular energy differences [Ha]
! dEijks(i,j,k,s) stores E_j-E_i, where i <= nbcder and j <= nb.
REAL(KIND=sp), ALLOCATABLE, INTENT(OUT) :: dEijks(:,:,:,:)

! WAVEDER metadata read from the first record
INTEGER(KIND=int32) :: &
    nb_file, &
    nbcder_file, &
    nktot_file, &
    nstot_file

! WAVEDER records not used by BerryCPT
REAL(KIND=dp) :: &
    nodes_in_dielectric_function, &
    wplasmon(3,3)

! Input unit and I/O status
INTEGER :: unitinp, ierr
CHARACTER(LEN=256) :: iomsg

! Band, k-point, spin, and Cartesian counters
INTEGER :: &
    iband, iband_cder, ikpt, ispin, idim

! Energy difference for the current band pair [Ha]
REAL(KIND=sp) :: energy_difference

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64) :: &
    pijks_bytes, &
    dEijks_bytes

! Unit conversion from Angstrom to bohr
REAL(KIND=sp), PARAMETER :: angstrom_to_bohr = 1.8897261_sp

CALL validate_dimensions( &
    fnameinp, nstot, nktot, nbcder, nb, ene) ! in

CALL validate_waveder_kinds()

OPEN(NEWUNIT=unitinp,FILE=TRIM(fnameinp),FORM='UNFORMATTED', &
    ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ', &
    IOSTAT=ierr,IOMSG=iomsg)
IF (ierr /= 0) THEN
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to open the WAVEDER file: '//TRIM(iomsg)) ! in
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) &
    nb_file, nbcder_file, nktot_file, nstot_file
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the WAVEDER metadata: '//TRIM(iomsg)) ! in
END IF

IF (nb_file /= nb .OR. nbcder_file /= nbcder .OR. &
        nktot_file /= nktot .OR. nstot_file /= nstot) THEN
    CLOSE(unitinp)
    WRITE(*,'(/,A)') ' ERROR: Inconsistent WAVEDER metadata'
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of bands: expected ', nb, ', found ', nb_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of derivative bands: expected ', nbcder, &
        ', found ', nbcder_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of k-points: expected ', nktot, ', found ', nktot_file
    WRITE(*,'(A,I0,A,I0)') &
        ' Number of spin channels: expected ', nstot, &
        ', found ', nstot_file
    ERROR STOP 1
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) nodes_in_dielectric_function
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the second WAVEDER record: '//TRIM(iomsg)) ! in
END IF

READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) wplasmon
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the third WAVEDER record: '//TRIM(iomsg)) ! in
END IF

CALL estimate_output_memory( &
    nstot, nktot, nbcder, nb, & ! in
    pijks_bytes, dEijks_bytes)   ! out

CALL report_output_memory( &
    pijks_bytes, dEijks_bytes) ! in

ALLOCATE( &
    pijks(3,nbcder,nb,nktot,nstot), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL allocation_error( &
        fnameinp, 'pijks', pijks_bytes, iomsg) ! in
END IF

ALLOCATE( &
    dEijks(nbcder,nb,nktot,nstot), &
    STAT=ierr, ERRMSG=iomsg)
IF (ierr /= 0) THEN
    DEALLOCATE(pijks)
    CLOSE(unitinp)
    CALL allocation_error( &
        fnameinp, 'dEijks', dEijks_bytes, iomsg) ! in
END IF

! WAVEDER stores D^a_ji with the all-band index j first and the
! derivative-band index i second. Read the complete record directly into
! pijks(:,i,j,k,s), avoiding a second full-size derivative-overlap array.
READ(unitinp,IOSTAT=ierr,IOMSG=iomsg) &
    (((((pijks(idim,iband_cder,iband,ikpt,ispin), &
        iband = 1, nb), &
        iband_cder = 1, nbcder), &
        ikpt = 1, nktot), &
        ispin = 1, nstot), &
        idim = 1, 3)
IF (ierr /= 0) THEN
    CLOSE(unitinp)
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to read the WAVEDER derivative overlaps: '//TRIM(iomsg)) ! in
END IF

CLOSE(unitinp,IOSTAT=ierr,IOMSG=iomsg)
IF (ierr /= 0) THEN
    CALL waveder_file_error( &
        fnameinp, & ! in
        'Unable to close the WAVEDER file: '//TRIM(iomsg)) ! in
END IF

! Convert D^a_ji [Angstrom] to p^a_ij [a.u.] using Hermiticity:
!
!   p^a_ji = (E_i-E_j) D^a_ji
!   p^a_ij = CONJG(p^a_ji)
!
! Only the rectangular block i <= nbcder is retained. For pairs entirely
! within that block, the reverse element and energy difference are set
! explicitly rather than taken independently from WAVEDER.
dEijks = 0.0_sp

DO ispin = 1, nstot
    DO ikpt = 1, nktot
        DO iband_cder = 1, nbcder
            pijks(:,iband_cder,iband_cder,ikpt,ispin) = &
                CMPLX(0.0_sp,0.0_sp,KIND=sp)

            DO iband = iband_cder + 1, nb
                energy_difference = ene(iband,ikpt,ispin) - &
                    ene(iband_cder,ikpt,ispin)

                dEijks(iband_cder,iband,ikpt,ispin) = &
                    energy_difference

                DO idim = 1, 3
                    pijks(idim,iband_cder,iband,ikpt,ispin) = &
                        -angstrom_to_bohr * energy_difference * &
                        CONJG(pijks(idim,iband_cder,iband,ikpt,ispin))
                END DO

                IF (iband <= nbcder) THEN
                    dEijks(iband,iband_cder,ikpt,ispin) = &
                        -energy_difference
                    pijks(:,iband,iband_cder,ikpt,ispin) = &
                        CONJG(pijks(:,iband_cder,iband,ikpt,ispin))
                END IF
            END DO
        END DO
    END DO
END DO

WRITE(*,'(A)') ' VASP momentum matrix elements were read successfully.'
WRITE(*,'(A,5(I0,1X),A)') &
    ' Shape of pijks = (', SHAPE(pijks), &
    ') (dim_3D, nbands_smaller, nbands_larger, nkpts, nspins)'
WRITE(*,'(A,4(I0,1X),A)') &
    ' Shape of dEijks = (', SHAPE(dEijks), &
    ') (nbands_smaller, nbands_larger, nkpts, nspins)'

END SUBROUTINE read_mommat_pij_vasp


SUBROUTINE validate_waveder_kinds()

! WAVEDER uses 32-bit integers, 64-bit real header values, and
! 32-bit real components in the complex derivative-overlap record.
IF (STORAGE_SIZE(0_int32) /= 32) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported integer storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' INTEGER(KIND=int32) storage size [bits] = ', &
        STORAGE_SIZE(0_int32)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(0.0_dp) /= 64) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported double-precision storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' REAL(KIND=dp) storage size [bits] = ', &
        STORAGE_SIZE(0.0_dp)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(0.0_sp) /= 32) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported single-precision storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' REAL(KIND=sp) storage size [bits] = ', &
        STORAGE_SIZE(0.0_sp)
    ERROR STOP 1
END IF

IF (STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp)) /= 64) THEN
    WRITE(*,'(/,A)') ' ERROR: Unsupported complex storage size for WAVEDER'
    WRITE(*,'(A,I0)') ' COMPLEX(KIND=sp) storage size [bits] = ', &
        STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp))
    ERROR STOP 1
END IF

END SUBROUTINE validate_waveder_kinds


SUBROUTINE estimate_output_memory( &
        nstot, nktot, nbcder, nb, &
        pijks_bytes, dEijks_bytes)

! WAVEDER dimensions
INTEGER, INTENT(IN) :: nstot, nktot, nbcder, nb

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64), INTENT(OUT) :: &
    pijks_bytes, &
    dEijks_bytes

! Number of band-pair, k-point, and spin entries
INTEGER(KIND=int64) :: matrix_entries

! Storage required by one array element [bytes]
INTEGER(KIND=int64) :: &
    complex_sp_bytes, &
    real_sp_bytes

matrix_entries = &
    INT(nbcder,KIND=int64) * INT(nb,KIND=int64) * &
    INT(nktot,KIND=int64) * INT(nstot,KIND=int64)

complex_sp_bytes = INT( &
    STORAGE_SIZE(CMPLX(0.0_sp,0.0_sp,KIND=sp)) / 8, &
    KIND=int64)
real_sp_bytes = INT(STORAGE_SIZE(0.0_sp) / 8,KIND=int64)

pijks_bytes = 3_int64 * matrix_entries * complex_sp_bytes
dEijks_bytes = matrix_entries * real_sp_bytes

END SUBROUTINE estimate_output_memory


SUBROUTINE report_output_memory( &
        pijks_bytes, dEijks_bytes)

! Memory required by the returned arrays [bytes]
INTEGER(KIND=int64), INTENT(IN) :: &
    pijks_bytes, &
    dEijks_bytes

! Combined memory required by the returned arrays [bytes]
INTEGER(KIND=int64) :: total_bytes

! Conversion from bytes to gibibytes
REAL(KIND=dp), PARAMETER :: bytes_per_gib = 1024.0_dp**3

total_bytes = pijks_bytes + dEijks_bytes

WRITE(*,'(/,A)') ' VASP momentum-matrix memory estimate'
WRITE(*,'(A,F12.3,A)') &
    ' pijks   = ', REAL(pijks_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A,F12.3,A)') &
    ' dEijks  = ', REAL(dEijks_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A,F12.3,A)') &
    ' combined = ', REAL(total_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A)') &
    ' This estimate excludes eigenvalues, occupations, work arrays, and runtime overhead.'
FLUSH(OUTPUT_UNIT)

END SUBROUTINE report_output_memory


SUBROUTINE allocation_error( &
        fnameinp, array_name, requested_bytes, message)

! File name, array name, and allocation error description
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, array_name, message

! Requested allocation size [bytes]
INTEGER(KIND=int64), INTENT(IN) :: requested_bytes

! Conversion from bytes to gibibytes
REAL(KIND=dp), PARAMETER :: bytes_per_gib = 1024.0_dp**3

WRITE(*,'(/,A)') ' ERROR: Unable to allocate a VASP momentum-matrix array'
WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
WRITE(*,'(A)') ' Array = '//TRIM(array_name)
WRITE(*,'(A,F12.3,A)') &
    ' Requested allocation = ', &
    REAL(requested_bytes,KIND=dp) / bytes_per_gib, ' GiB'
WRITE(*,'(A)') ' '//TRIM(message)
ERROR STOP 1

END SUBROUTINE allocation_error


SUBROUTINE validate_dimensions( &
        fnameinp, nstot, nktot, nbcder, nb, ene)

! File name and WAVEDER metadata
CHARACTER(LEN=*), INTENT(IN) :: fnameinp
INTEGER, INTENT(IN) :: nstot, nktot, nbcder, nb

! Band-, k-point-, and spin-dependent eigenvalues [Ha]
REAL(KIND=sp), INTENT(IN) :: ene(:,:,:)

IF (nstot < 1 .OR. nstot > 2 .OR. nktot <= 0 .OR. nb <= 0 .OR. &
        nbcder <= 0 .OR. nbcder > nb) THEN
    WRITE(*,'(/,A)') ' ERROR: Invalid WAVEDER dimensions were supplied'
    WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
    WRITE(*,'(A,I0)') ' Number of spin channels = ', nstot
    WRITE(*,'(A,I0)') ' Number of k-points = ', nktot
    WRITE(*,'(A,I0)') ' Number of derivative bands = ', nbcder
    WRITE(*,'(A,I0)') ' Number of bands = ', nb
    ERROR STOP 1
END IF

IF (ANY(SHAPE(ene) /= [nb,nktot,nstot])) THEN
    WRITE(*,'(/,A)') ' ERROR: Inconsistent dimensions of ene'
    WRITE(*,'(A,3(I0,1X))') ' Actual shape = ', SHAPE(ene)
    WRITE(*,'(A,3(I0,1X))') &
        ' Expected shape = ', nb, nktot, nstot
    ERROR STOP 1
END IF

END SUBROUTINE validate_dimensions


SUBROUTINE waveder_file_error( &
        fnameinp, message)

! WAVEDER file name and error description
CHARACTER(LEN=*), INTENT(IN) :: fnameinp, message

WRITE(*,'(/,A)') ' ERROR: Unable to read VASP momentum matrix elements'
WRITE(*,'(A)') ' File = '//TRIM(fnameinp)
WRITE(*,'(A)') ' '//TRIM(message)
ERROR STOP 1

END SUBROUTINE waveder_file_error

END MODULE read_mommat_pij_vasp_mod
