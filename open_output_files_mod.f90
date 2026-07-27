MODULE open_output_files_mod

    USE version_mod, ONLY : version_string

    IMPLICIT NONE

    PRIVATE

    ! Output-file units used by the k-point writing routines
    TYPE, PUBLIC :: output_units_type
        INTEGER :: bcurv = -1
        INTEGER :: bcurv_up = -1
        INTEGER :: bcurv_dn = -1
        INTEGER :: bcurv_spin = -1
        INTEGER :: oam = -1
        INTEGER :: oam_up = -1
        INTEGER :: goam = -1
    END TYPE output_units_type

    ! Description retained for every output file until the final summary
    TYPE :: output_record_type
        CHARACTER(LEN=512) :: filename = ''
        CHARACTER(LEN=160) :: quantity = ''
        CHARACTER(LEN=240) :: values = ''
        CHARACTER(LEN=32) :: units = ''
    END TYPE output_record_type

    ! Output-file manifest accumulated over all spin channels
    TYPE(output_record_type), ALLOCATABLE, SAVE :: output_records(:)

    PUBLIC :: &
        open_output_files, &
        close_output_files, &
        write_output_summary

CONTAINS

    SUBROUTINE open_output_files( &
        wien2k, soc_sp_resolv_pij, charspin, & ! in
        output_units)                         ! out

        IMPLICIT NONE

        ! Calculation type and output-file suffix
        LOGICAL, INTENT(IN) :: &
            wien2k, &
            soc_sp_resolv_pij
        CHARACTER(LEN=*), INTENT(IN) :: charspin

        ! Units of the opened output files
        TYPE(output_units_type), INTENT(OUT) :: output_units

        ! Output file name
        CHARACTER(LEN=:), ALLOCATABLE :: filename

        output_units = output_units_type()

        !! Ordinary Berry curvature

        filename = 'bcurv_ij' // TRIM(charspin) // '.dat'
        CALL open_formatted_output_file( &
            filename, &                 ! in
            output_units%bcurv)         ! out
        CALL write_bcurv_header( &
            output_units%bcurv, filename, &             ! in
            'the Berry curvature tensor Omega_ij')      ! in

        !! Spin-resolved Berry-curvature contributions

        IF (wien2k .AND. soc_sp_resolv_pij) THEN
            filename = 'bcurv_ij-up.dat'
            CALL open_formatted_output_file( &
                filename, &                      ! in
                output_units%bcurv_up)           ! out
            CALL write_bcurv_header( &
                output_units%bcurv_up, filename, &                 ! in
                'the spin UP Berry curvature tensor Omega_ij')     ! in

            filename = 'bcurv_ij-dn.dat'
            CALL open_formatted_output_file( &
                filename, &                      ! in
                output_units%bcurv_dn)           ! out
            CALL write_bcurv_header( &
                output_units%bcurv_dn, filename, &                 ! in
                'the spin DN Berry curvature tensor Omega_ij')     ! in

            filename = 'bcurv_ij-up-dn.dat'
            CALL open_formatted_output_file( &
                filename, &                      ! in
                output_units%bcurv_spin)         ! out
            CALL write_bcurv_header( &
                output_units%bcurv_spin, filename, &                ! in
                'the spin UP-DN Berry curvature tensor Omega_ij')   ! in
        END IF

        !! Ordinary orbital angular momentum

        filename = 'oam_ij' // TRIM(charspin) // '.dat'
        CALL open_formatted_output_file( &
            filename, &              ! in
            output_units%oam)        ! out
        CALL write_oam_header( &
            output_units%oam, filename) ! in

        !! Spin-up orbital angular momentum

        IF (wien2k .AND. soc_sp_resolv_pij) THEN
            filename = 'oam_ij-sigma_z-up.dat'
            CALL open_formatted_output_file( &
                filename, &                 ! in
                output_units%oam_up)        ! out
            CALL write_oam_up_header( &
                output_units%oam_up, filename) ! in
        END IF

        !! Generalized orbital angular momentum

        filename = 'goam_ij_nm' // TRIM(charspin) // '.dat'
        CALL open_formatted_output_file( &
            filename, &               ! in
            output_units%goam)        ! out
        CALL write_goam_header( &
            output_units%goam, filename) ! in

    END SUBROUTINE open_output_files


    SUBROUTINE open_formatted_output_file( &
        filename, & ! in
        unit)       ! out

        IMPLICIT NONE

        ! Output file name
        CHARACTER(LEN=*), INTENT(IN) :: filename

        ! Connected Fortran unit
        INTEGER, INTENT(OUT) :: unit

        ! Input/output status information
        INTEGER :: ierr
        CHARACTER(LEN=512) :: iomsg

        OPEN( &
            NEWUNIT=unit, &
            FILE=TRIM(filename), &
            STATUS='REPLACE', &
            FORM='FORMATTED', &
            ACTION='WRITE', &
            IOSTAT=ierr, &
            IOMSG=iomsg)

        IF (ierr /= 0) THEN
            WRITE(*,'(/,A)') ' ERROR: Unable to open output file'
            WRITE(*,'(A)') ' File: ' // TRIM(filename)
            WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
            ERROR STOP 1
        END IF

    END SUBROUTINE open_formatted_output_file


    SUBROUTINE write_bcurv_header( &
        unit, filename, description) ! in

        IMPLICIT NONE

        ! Output-file unit, name, and quantity description
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: &
            filename, &
            description

        ! Header lines
        CHARACTER(LEN=160) :: lines(11)

        lines(1) = '# This file is generated by berrycpt. ' // TRIM(version_string)
        lines(2) = '# The output contains components of ' // TRIM(description) // '.'
        lines(3) = '# Data are grouped by k-point index and then by band index.'
        lines(4) = '# Each band record also gives its local degenerate-block ID.'
        lines(5) = '# Bands with the same positive block ID were treated together'
        lines(6) = '# by degenerate perturbation theory at that k-point and spin.'
        lines(7) = '# The final record for each k-point uses band=0 and block=0 and'
        lines(8) = '# gives the occupation-weighted total Berry curvature.'
        lines(9) = '# Columns: band  block  Omega_yz  Omega_zx  Omega_xy'
        lines(10) = '# Tensor labels: yz=Voigt 4, zx=Voigt 5, xy=Voigt 6'
        lines(11) = '# Units: bohr^2.'

        CALL write_output_lines( &
            unit, filename, lines) ! in

    END SUBROUTINE write_bcurv_header


    SUBROUTINE write_oam_header( &
        unit, filename) ! in

        IMPLICIT NONE

        ! Output-file unit and name
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: filename

        ! Header lines
        CHARACTER(LEN=160) :: lines(9)

        lines(1) = '# This file is generated by berrycpt. ' // TRIM(version_string)
        lines(2) = '# The output contains components of the orbital angular momentum tensor L_ij.'
        lines(3) = '# Data are grouped by k-point index and then by band index.'
        lines(4) = '# Each band record also gives its local degenerate-block ID.'
        lines(5) = '# Bands with the same positive block ID were treated together'
        lines(6) = '# by degenerate perturbation theory at that k-point and spin.'
        lines(7) = '# Columns: band  block  L_yz  L_zx  L_xy'
        lines(8) = '# Tensor labels: yz=Voigt 4, zx=Voigt 5, xy=Voigt 6'
        lines(9) = '# Units: hbar.'

        CALL write_output_lines( &
            unit, filename, lines) ! in

    END SUBROUTINE write_oam_header


    SUBROUTINE write_oam_up_header( &
        unit, filename) ! in

        IMPLICIT NONE

        ! Output-file unit and name
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: filename

        ! Header lines
        CHARACTER(LEN=160) :: lines(12)

        lines(1) = '# This file is generated by berrycpt. ' // TRIM(version_string)
        lines(2) = '# The output contains components of the spin-up-projected'
        lines(3) = '# orbital angular momentum tensor expressed as'
        lines(4) = '# L^{sigma_z,up}_{ij,n} = <u_n|L_ij P_up|u_n>,'
        lines(5) = '# where P_up = 1/2*(I+sigma_z).'
        lines(6) = '# Data are grouped by k-point index and then by band index.'
        lines(7) = '# Each band record also gives its local degenerate-block ID.'
        lines(8) = '# Bands with the same positive block ID were treated together'
        lines(9) = '# by degenerate perturbation theory at that k-point and spin.'
        lines(10) = '# Columns: band  block  L_yz  L_zx  L_xy'
        lines(11) = '# Tensor labels: yz=Voigt 4, zx=Voigt 5, xy=Voigt 6'
        lines(12) = '# Units: hbar.'

        CALL write_output_lines( &
            unit, filename, lines) ! in

    END SUBROUTINE write_oam_up_header


    SUBROUTINE write_goam_header( &
        unit, filename) ! in

        IMPLICIT NONE

        ! Output-file unit and name
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: filename

        ! Header lines
        CHARACTER(LEN=180) :: lines(18)

        lines(1) = '# This file is generated by berrycpt. ' // TRIM(version_string)
        lines(2) = '# The output contains generalized orbital-angular-momentum matrices.'
        lines(3) = '# Matrix convention: L_c(i,j) = <u_i|L_c|u_j>.'
        lines(4) = '# The first matrix index is the bra band and the second is the ket band.'
        lines(5) = '# Both matrix indices span bands 1:NBANDS_OUT.'
        lines(6) = '# Intermediate-state sums used bands 1:NBANDS_TRANS.'
        lines(7) = '# Components are yz=Voigt 4, zx=Voigt 5, and xy=Voigt 6.'
        lines(8) = '# Each component is a Hermitian matrix. The lower triangle was'
        lines(9) = '# reconstructed as the complex conjugate of the upper triangle.'
        lines(10) = '# Complex values are printed as (Re,Im).'
        lines(11) = '# Each matrix row begins with the bra-band index and its local'
        lines(12) = '# degenerate-block ID. Remaining entries correspond, in order,'
        lines(13) = '# to ket bands 1:NBANDS_OUT.'
        lines(14) = '# Degenerate-block IDs identify subspaces excluded from singular'
        lines(15) = '# intermediate-state couplings. The matrix remains in the original'
        lines(16) = '# DFT eigenstate basis and is not diagonalized within those blocks.'
        lines(17) = '# Row format: band  block  (Re L_i1,Im L_i1) ... (Re L_iN,Im L_iN)'
        lines(18) = '# Units: hbar.'

        CALL write_output_lines( &
            unit, filename, lines) ! in

    END SUBROUTINE write_goam_header


    SUBROUTINE write_output_lines( &
        unit, filename, lines) ! in

        IMPLICIT NONE

        ! Output-file unit, name, and header lines
        INTEGER, INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: &
            filename, &
            lines(:)

        ! Input/output status information
        INTEGER :: &
            ierr, &
            iline
        CHARACTER(LEN=512) :: iomsg

        DO iline = 1, SIZE(lines)
            WRITE(unit,'(A)',IOSTAT=ierr,IOMSG=iomsg) TRIM(lines(iline))

            IF (ierr /= 0) THEN
                WRITE(*,'(/,A)') ' ERROR: Unable to write an output-file header'
                WRITE(*,'(A)') ' File: ' // TRIM(filename)
                WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
                ERROR STOP 1
            END IF
        END DO

    END SUBROUTINE write_output_lines


    SUBROUTINE close_output_files( &
        output_units) ! in/out

        IMPLICIT NONE

        ! Units of output files opened for the current spin channel
        TYPE(output_units_type), INTENT(INOUT) :: output_units

        CALL close_and_record_output_file( &
            output_units%bcurv, & ! in/out
            'Ordinary Berry curvature', & ! in
            'Omega_yz, Omega_zx, Omega_xy by band, plus the occupation-weighted total for each k-point', & ! in
            'bohr^2') ! in

        CALL close_and_record_output_file( &
            output_units%bcurv_up, & ! in/out
            'Spin-up contribution to the Berry curvature', & ! in
            'Omega_yz, Omega_zx, Omega_xy by band, plus the occupation-weighted total for each k-point', & ! in
            'bohr^2') ! in

        CALL close_and_record_output_file( &
            output_units%bcurv_dn, & ! in/out
            'Spin-down contribution to the Berry curvature', & ! in
            'Omega_yz, Omega_zx, Omega_xy by band, plus the occupation-weighted total for each k-point', & ! in
            'bohr^2') ! in

        CALL close_and_record_output_file( &
            output_units%bcurv_spin, & ! in/out
            'Sigma_z-normalized spin Berry curvature', & ! in
            'Omega_yz, Omega_zx, Omega_xy by band, plus the occupation-weighted total for each k-point', & ! in
            'bohr^2') ! in

        CALL close_and_record_output_file( &
            output_units%oam, & ! in/out
            'Ordinary orbital angular momentum', & ! in
            'L_yz, L_zx, L_xy by band at each k-point', & ! in
            'hbar') ! in

        CALL close_and_record_output_file( &
            output_units%oam_up, & ! in/out
            'Spin-up-projected orbital angular momentum', & ! in
            'L_yz, L_zx, L_xy by band at each k-point', & ! in
            'hbar') ! in

        CALL close_and_record_output_file( &
            output_units%goam, & ! in/out
            'Generalized orbital-angular-momentum matrices', & ! in
            'complex <u_i|L_c|u_j> for c = yz, zx, xy in the original DFT eigenstate basis', & ! in
            'hbar') ! in

    END SUBROUTINE close_output_files


    SUBROUTINE close_and_record_output_file( &
        unit, quantity, values, units) ! in/out, in, in, in

        IMPLICIT NONE

        ! Connected output-file unit, reset after successful closing
        INTEGER, INTENT(INOUT) :: unit

        ! Description printed in the final output-file manifest
        CHARACTER(LEN=*), INTENT(IN) :: &
            quantity, &
            values, &
            units

        ! Actual file name associated with the Fortran unit
        CHARACTER(LEN=512) :: filename

        ! Input/output status information
        INTEGER :: ierr
        LOGICAL :: opened
        CHARACTER(LEN=512) :: iomsg

        IF (unit == -1) RETURN

        filename = ''
        INQUIRE( &
            UNIT=unit, &
            OPENED=opened, &
            NAME=filename, &
            IOSTAT=ierr, &
            IOMSG=iomsg)

        IF (ierr /= 0) THEN
            WRITE(*,'(/,A)') ' ERROR: Unable to inquire about an output file'
            WRITE(*,'(A,I0)') ' Unit: ', unit
            WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
            ERROR STOP 1
        END IF

        IF (.NOT. opened) THEN
            WRITE(*,'(/,A)') ' ERROR: Output-file unit is not connected'
            WRITE(*,'(A,I0)') ' Unit: ', unit
            ERROR STOP 1
        END IF

        CLOSE( &
            UNIT=unit, &
            STATUS='KEEP', &
            IOSTAT=ierr, &
            IOMSG=iomsg)

        IF (ierr /= 0) THEN
            WRITE(*,'(/,A)') ' ERROR: Unable to close an output file'
            WRITE(*,'(A)') ' File: ' // TRIM(filename)
            WRITE(*,'(A)') ' Reason: ' // TRIM(iomsg)
            ERROR STOP 1
        END IF

        unit = -1

        CALL append_output_record( &
            filename, quantity, values, units) ! in

    END SUBROUTINE close_and_record_output_file


    SUBROUTINE append_output_record( &
        filename, quantity, values, units) ! in

        IMPLICIT NONE

        ! Output-file information retained for the final manifest
        CHARACTER(LEN=*), INTENT(IN) :: &
            filename, &
            quantity, &
            values, &
            units

        ! Expanded copy of the output-file manifest
        TYPE(output_record_type), ALLOCATABLE :: expanded_records(:)

        ! Number of records already stored
        INTEGER :: number_of_records

        IF (.NOT. ALLOCATED(output_records)) THEN
            ALLOCATE(output_records(1))
        ELSE
            number_of_records = SIZE(output_records)
            ALLOCATE(expanded_records(number_of_records + 1))
            expanded_records(1:number_of_records) = output_records
            CALL MOVE_ALLOC(expanded_records, output_records)
        END IF

        output_records(SIZE(output_records))%filename = TRIM(filename)
        output_records(SIZE(output_records))%quantity = TRIM(quantity)
        output_records(SIZE(output_records))%values = TRIM(values)
        output_records(SIZE(output_records))%units = TRIM(units)

    END SUBROUTINE append_output_record


    SUBROUTINE write_output_summary()

        IMPLICIT NONE

        ! Output-file record index
        INTEGER :: irecord

        IF (.NOT. ALLOCATED(output_records)) THEN
            WRITE(*,'(/,A)') ' No output files were recorded.'
            RETURN
        END IF

        WRITE(*,'(/,A)') ' Output files created and closed successfully:'

        DO irecord = 1, SIZE(output_records)
            WRITE(*,'(/,2X,A)') TRIM(output_records(irecord)%filename)
            WRITE(*,'(4X,A)') &
                'Quantity: ' // TRIM(output_records(irecord)%quantity)
            WRITE(*,'(4X,A)') &
                'Values: ' // TRIM(output_records(irecord)%values)
            WRITE(*,'(4X,A)') &
                'Units: ' // TRIM(output_records(irecord)%units)
        END DO

        WRITE(*,'(/,A)') &
            ' Detailed k-point, band, degenerate-block, and matrix layouts'
        WRITE(*,'(A)') &
            ' are documented in the header of each output file.'

    END SUBROUTINE write_output_summary

END MODULE open_output_files_mod
