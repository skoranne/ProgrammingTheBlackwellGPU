### Just something so we dont have to retype library paths
FC=nvfortran
FOPTS=-O3
FOPTS+=-stdpar -gpu=nvlamath -cudalib=nvlamath
MATH_LIBS=/opt/nvidia/hpc_sdk/Linux_aarch64/25.9/math_libs/lib64
xt.exe:	xt.f90
	$(FC) $(FOPTS) xt.f90 -L $(MATH_LIBS) -lcublas -o xt.exe

