# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to select the Fortran compiler and compilation flags.

# Intel Fortran

FC = ifort

# Intel release build
FCFLAGS = -I${MKLROOT}/include -O2 -warn all
FLFLAGS = -L${MKLROOT}/lib/intel64 \
          -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
          -lpthread -ldl

# Intel debugging build
# FCFLAGS = -I${MKLROOT}/include -O0 -g -traceback -check all \
#           -warn all -debug all
# FLFLAGS = -L${MKLROOT}/lib/intel64 \
#           -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
#           -lpthread -ldl -g -traceback

# GNU Fortran

# FC = gfortran

# GNU Fortran release build
# FCFLAGS = -I${MKLROOT}/include -O2 -Wall -Wextra \
#           -Wsurprising -Wcharacter-truncation
# FLFLAGS = -L${MKLROOT}/lib/intel64 \
#           -lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
#           -lpthread -ldl

# GNU Fortran debugging build
# FCFLAGS = -I${MKLROOT}/include -O0 -g -fbacktrace -fcheck=all \
#           -Wall -Wextra -Wsurprising -Wcharacter-truncation \
#           -Wintrinsics-std -fno-omit-frame-pointer \
#           -ffpe-summary=none
# FLFLAGS = -L${MKLROOT}/lib/intel64 \
#           -lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
#           -lpthread -ldl -g -fbacktrace -ffpe-summary=none


# ~~~ Do not edit after that line ~~~

ifndef MKLROOT
    $(error MKLROOT environment variable is not set. Please check your MKL setup.)
endif

PROGRAM = berrycpt

MODULE_OBJS = \
	precision_mod.o \
	is_hermitian_mod.o \
	eigvz_mod.o \
	degenbc_mod.o \
	degenoam_mod.o \
	calculate_bcurv_kpoint_mod.o \
	calculate_goam_kpoint_mod.o \
	calculate_oam_kpoint_mod.o \
	command_line_args_mod.o \
	construct_occupations_mod.o \
	filename_contains_energyso_mod.o \
	find_degenerate_groups_mod.o \
	open_mommat_files_mod.o \
	open_output_files_mod.o \
	read_eigenvalues_vasp_mod.o \
	read_eigenvalues_wien2k_mod.o \
	read_mommat_pij_vasp_mod.o \
	read_mommat_pij_wien2k_kpoint_mod.o \
	read_occupations_wien2k_mod.o \
	read_waveder_header_mod.o \
	set_spin_suffix_mod.o \
	validate_input_files_mod.o \
	write_bcurv_kpoint_mod.o \
	write_goam_kpoint_mod.o \
	write_oam_kpoint_mod.o \
	write_progress_mod.o

OBJS = $(MODULE_OBJS) berrycpt.o

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(FLFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Module dependencies

eigvz_mod.o: precision_mod.o

degenbc_mod.o: precision_mod.o eigvz_mod.o is_hermitian_mod.o
degenoam_mod.o: precision_mod.o eigvz_mod.o is_hermitian_mod.o

calculate_bcurv_kpoint_mod.o: precision_mod.o degenbc_mod.o
calculate_goam_kpoint_mod.o: precision_mod.o
calculate_oam_kpoint_mod.o: precision_mod.o degenoam_mod.o

command_line_args_mod.o: precision_mod.o
construct_occupations_mod.o: precision_mod.o
find_degenerate_groups_mod.o: precision_mod.o

read_eigenvalues_vasp_mod.o: precision_mod.o
read_eigenvalues_wien2k_mod.o: precision_mod.o
read_mommat_pij_vasp_mod.o: precision_mod.o
read_mommat_pij_wien2k_kpoint_mod.o: precision_mod.o
read_occupations_wien2k_mod.o: precision_mod.o
read_waveder_header_mod.o: precision_mod.o

write_bcurv_kpoint_mod.o: precision_mod.o
write_goam_kpoint_mod.o: precision_mod.o
write_oam_kpoint_mod.o: precision_mod.o
write_progress_mod.o: precision_mod.o

berrycpt.o: $(MODULE_OBJS)

# Utility targets
.PHONY: all clean veryclean

clean:
	rm -f *.o *.mod *.MOD *__genmod.*

veryclean: clean
	rm -f *~ $(PROGRAM)
