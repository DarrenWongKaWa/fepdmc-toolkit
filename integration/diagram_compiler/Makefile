# Validation-only. Does not replace FEP-DMC update_swap or public P1.
# Compiles stage-5 generated kernel + native_kernel_v1 runtime.
FC ?= gfortran
FFLAGS ?= -O0 -std=f2008 -ffree-line-length-none -fcheck=bounds -J.

.PHONY: all clean

all: dc_validate.x

dc_validate.x: native_kernel_runtime_v1.f90 generated_da_kernel.f90 dc_validate.f90
	$(FC) $(FFLAGS) native_kernel_runtime_v1.f90 generated_da_kernel.f90 dc_validate.f90 -o dc_validate.x

clean:
	rm -f dc_validate.x *.mod *.o
