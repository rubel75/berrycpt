# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags.

# Intel Fortran

FC = ifort

# Intel debuging
FCFLAGS = -I${MKLROOT}/include -g -traceback -check all -warn all -debug all -qopenmp -O0 #-heap-arrays
FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl -g -traceback -check all -warn all -debug all -O0

# Intel performance
#FCFLAGS =  -I${MKLROOT}/include -qopenmp -heap-arrays
#FLFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl

# GNU Fortran

#FC = gfortran
#FCFLAGS = -I${MKLROOT}/include -fopenmp -g -fbacktrace -ffpe-summary=none -fno-automatic
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -fopenmp -ldl -g -fbacktrace -ffpe-summary=none

# ~~~ Do not edit after that line ~~~

# Check env variables
ifndef MKLROOT
    $(error MKLROOT environment variable is not set. Please check your MKL setup.)
endif

PROGRAM = berrycpt

# Get all .f90 source files
F90S := $(filter-out %__genmod.f90, $(wildcard *.f90))

# Convert to .o files
ALL_OBJS := $(patsubst %.f90, %.o, $(F90S))

# Remove berrycpt.o from the list
OBJS := $(filter-out berrycpt.o precision_mod.o eigvd_mod.o, $(ALL_OBJS))

# Add precision_mod.o first and berrycpt.o last
OBJS := precision_mod.o eigvd_mod.o $(OBJS) berrycpt.o


all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

#%.mod: %.h
#	$(FC) $(FCFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *__genmod.*

veryclean: clean
	rm -f *~ $(PROGRAM)
