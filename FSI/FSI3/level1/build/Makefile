TARGET := fsi3l1.x
BUILD_DIR := 
SRC_DIR := ../source
# SRCS := $(shell find $(SRC_DIR) -name '*.f90')
# SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRCS := \
$(SRC_DIR)/modInput.f90 \
$(SRC_DIR)/modFEM.f90 \
$(SRC_DIR)/modMath.f90 \
$(SRC_DIR)/modLBM.f90 \
$(SRC_DIR)/modSolver.f90 \
$(SRC_DIR)/turekFSI.f90 
OBJS := ${SRCS:$(SRC_DIR)/%.f90=%.o} #substitute .f90 with .o

# FC=ifort #fortran compiler
# FFLAGS=-mkl -xCORE-AVX2# -qopenmp #-CB -traceback
FC := gfortran
FFLAGS := -g -fopenmp #-CB -traceback
LIBS := -llapack -lblas
# LIBS= -L/usr/local/lib/ -llapack -lblas


$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: $(SRC_DIR)/%.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LIBS)

.PHONY: clean
clean:
	rm *mod *o *x
