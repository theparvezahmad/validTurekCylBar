TARGET := fsi3l1.x
BUILD_DIR :=
SRC_DIR := ../source/

SRCS := \
$(SRC_DIR)modInput.f90 \
$(SRC_DIR)modFEM.f90 \
$(SRC_DIR)modMath.f90 \
$(SRC_DIR)modLBM.f90 \
$(SRC_DIR)modModel.f90 \
$(SRC_DIR)modSolver.f90 \
$(SRC_DIR)turekFSI.f90

OBJS := ${SRCS:$(SRC_DIR)%.f90=$(BUILD_DIR)%.o}

FC := ifort #fortran compiler
FFLAGS :=-mkl# -qopenmp #-CB -traceback
# FC := gfortran
# FFLAGS := -g# -fopenmp #-CB -traceback
# LIBS := -llapack -lblas
# LIBS= -L/usr/local/lib/ -llapack -lblas


$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

$(BUILD_DIR)%.o: $(SRC_DIR)%.f90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LIBS)

.PHONY: clean
clean:
	@rm -f -v $(BUILD_DIR)*.mod *.o *.x
