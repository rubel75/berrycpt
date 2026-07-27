# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean
#
# Version information is generated automatically at build time.
#
# For a cloned Git repository:
#
#     Git revision: YYYY-MM-DD-<12-character commit hash>
#
# If tracked files differ from HEAD, "-dirty" is appended.
#
# For a source archive or other copy without Git metadata:
#
#     Git revision: unavailable; build date: YYYY-MM-DD


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


PROGRAM = berrycpt

VERSION_SRC = version_mod.f90


MODULE_OBJS = \
	precision_mod.o \
	version_mod.o \
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


$(PROGRAM): $(OBJS) | check-env
	$(FC) -o $@ $^ $(FLFLAGS)


%.o: %.f90 | check-env
	$(FC) $(FCFLAGS) -c $<


# Version information
#
# If this is a cloned Git repository, obtain the commit date and
# abbreviated commit hash directly from Git.  If tracked files have
# been modified relative to HEAD, append "-dirty".
#
# If Git is unavailable, or the source tree does not contain Git
# metadata (for example, a ZIP downloaded from GitHub), report that
# the Git revision is unavailable and record the build date instead.
#
# A temporary source file is generated and compared with the existing
# version_mod.f90.  The file is replaced only when its contents change,
# avoiding unnecessary recompilation when the version is unchanged.

$(VERSION_SRC): FORCE
	@version=""; \
	if command -v git >/dev/null 2>&1 && \
	   [ -e .git ] && \
	   git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		git_date=$$(git log -1 --format=%cs 2>/dev/null); \
		git_hash=$$(git rev-parse --short=12 HEAD 2>/dev/null); \
		if [ -n "$$git_date" ] && [ -n "$$git_hash" ]; then \
			dirty=""; \
			if ! git diff --quiet HEAD -- 2>/dev/null; then \
				dirty="-dirty"; \
			fi; \
			version="Git revision: $$git_date-$$git_hash$$dirty"; \
		fi; \
	fi; \
	if [ -z "$$version" ]; then \
		build_date=$$(date +%Y-%m-%d); \
		version="Git revision: unavailable; build date: $$build_date"; \
	fi; \
	{ \
		echo 'module version_mod'; \
		echo '    implicit none'; \
		echo '    private'; \
		echo '    public :: version_string'; \
		echo ''; \
		printf '    character(len=*), parameter :: version_string = "%s"\n' "$$version"; \
		echo ''; \
		echo 'end module version_mod'; \
	} > $(VERSION_SRC).tmp; \
	if [ ! -f $(VERSION_SRC) ] || \
	   ! cmp -s $(VERSION_SRC).tmp $(VERSION_SRC); then \
		mv $(VERSION_SRC).tmp $(VERSION_SRC); \
	else \
		rm -f $(VERSION_SRC).tmp; \
	fi


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

.PHONY: all clean veryclean check-env FORCE

FORCE:


clean:
	rm -f *.o *.mod *.MOD *__genmod.*


veryclean: clean
	rm -f *~ $(PROGRAM) $(VERSION_SRC) $(VERSION_SRC).tmp


check-env:
	@if [ -z "$(strip $(FC))" ]; then \
		echo "ERROR: FC is not set."; \
		exit 1; \
	fi
	@if ! command -v "$(firstword $(FC))" >/dev/null 2>&1; then \
		echo "ERROR: Fortran compiler '$(firstword $(FC))' was not found."; \
		echo "       Set FC to an available compiler."; \
		exit 1; \
	fi
	@if ! $(FC) --version >/dev/null 2>&1; then \
		echo "ERROR: Fortran compiler command '$(FC)' cannot be executed."; \
		exit 1; \
	fi
	@if [ -z "$(MKLROOT)" ]; then \
		echo "ERROR: MKLROOT environment variable is not set."; \
		exit 1; \
	fi
	@if [ ! -d "$(MKLROOT)/include" ]; then \
		echo "ERROR: MKL include directory not found:"; \
		echo "       $(MKLROOT)/include"; \
		exit 1; \
	fi
	@if [ ! -d "$(MKLROOT)/lib/intel64" ]; then \
		echo "ERROR: MKL library directory not found:"; \
		echo "       $(MKLROOT)/lib/intel64"; \
		exit 1; \
	fi