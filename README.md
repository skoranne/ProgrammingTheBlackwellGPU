This repository contains nvfortran based experiments on the GB10 Blackwell
also called the DGX Spark.

The problems are chosen to explain ISO language features which map well
onto the GPU.

-- Sandeep

Steps to setup the compiler toolchain

export MODULEPATH=/opt/nvidia/hpc_sdk/modulefiles/:$MODULEPATH

---------------------------------------------- /opt/nvidia/hpc_sdk/modulefiles -----------------------------------------------
nvhpc-byo-compiler/25.9  nvhpc-hpcx-cuda13/25.9  nvhpc-hpcx/25.9  nvhpc-nompi/25.9  nvhpc/25.9  

----------------------------------------------- /usr/share/modules/modulefiles -----------------------------------------------
dot  module-git  module-info  modules  null  use.own  
module load nvhpc-hpcx-cuda13/25.9

nvfortran 25.9-0 linuxarm64 target on aarch64 Linux -tp generic 
NVIDIA Compilers and Tools
Copyright (c) 2025, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
