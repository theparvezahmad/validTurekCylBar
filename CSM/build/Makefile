TARGET := csm.x
BUILD_DIR :=
SRC_DIR := ../src/

SRCS := \
$(SRC_DIR)modUserInput.f90 \
$(SRC_DIR)modShapeFun_v2.f90 \
$(SRC_DIR)modMathOp.f90 \
$(SRC_DIR)turekDynamic_v2.f90

OBJS := ${SRCS:$(SRC_DIR)%.f90=$(BUILD_DIR)%.o}

FC := ifort #fortran compiler
FFLAGS :=-qmkl# -CB -traceback# -xCORE-AVX2# -qopenmp #-CB -traceback
# FFLAGS :=-qmkl -g -check all -fpe0 -warn -traceback -debug extended
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
