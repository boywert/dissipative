
# AREPO Makefile
#   see documentation/getting_started.md
#
# If you add a new system below, also add that systype to Template-Makefile.systype

EXEC   = Arepo
LIBRARY = arepo
CONFIG   = Config.sh
BUILD_DIR = build
SRC_DIR = src

###################
#determine SYSTYPE#
###################
ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

MAKEFILES = Makefile config-makefile
ifeq ($(wildcard Makefile.systype), Makefile.systype)
MAKEFILES += Makefile.systype
endif

$(info Build configuration:)
$(info SYSTYPE: $(SYSTYPE))
$(info CONFIG: $(CONFIG))
$(info EXEC: $(EXEC))
$(info )

PYTHON = python
PERL	 = /usr/bin/perl
RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) BUILD_DIR=$(BUILD_DIR) make -f config-makefile)
CONFIGVARS := $(shell cat $(BUILD_DIR)/arepoconfig.h)
RESULT     := $(shell SRC_DIR=$(SRC_DIR) BUILD_DIR=$(BUILD_DIR) ./git_version.sh)

MPICHLIB  = #-lmpich
GMPLIB    = -lgmp
GSLLIB    = -lgsl -lgslcblas
MATHLIB   = -lm -lstdc++
HWLOC_LIB = -lhwloc

ifeq (MCMA,$(findstring MCMA,$(CONFIGVARS)))
LAPACK_LIB = -llapack
endif

#OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems

ifeq ($(SYSTYPE),"Curie")
CC         =  mpicc  

OPTIMIZE   =  -std=c11 -g -O2 #-xW -ipo -Wall
GMP_INCL   = -I/ccc/work/cont005/ra0844/springev/libs/gmp/include
GMP_LIBS   = -L/ccc/work/cont005/ra0844/springev/libs/gmp/lib  -Xlinker -R -Xlinker /ccc/work/cont005/ra0844/springev/libs/gmp/lib 
GSL_INCL   = -I$(GSL_INC_DIR)
GSL_LIBS   = -L$(GSL_LIB_DIR)
FFTW_INCL  = -I$(FFTW2_INC_DIR)
FFTW_LIBS  = -L$(FFTW2_LIB_DIR)
HDF5INCL   = -I$(HDF5_INC_DIR) -DH5_USE_16_API
HDF5LIB    = -L$(HDF5_LIB_DIR) -lhdf5 -lz
HWLOC_INCL = -I/ccc/work/cont005/ra0844/springev/libs/hwloc/include/
HWLOC_LIB  = -L/ccc/work/cont005/ra0844/springev/libs/hwloc/lib/ -lhwloc -Xlinker -R -Xlinker /ccc/work/cont005/ra0844/springev/libs/hwloc/lib/

MPICHLIB   =
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
endif


ifeq ($(SYSTYPE),"Gordon") #module load fftw/2.1.5 hdf5/1.8.11 gsl/1.15
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DNOCALLSOFSYSTEM  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =  -std=c99 -O3 -g -Wall #-mavx #-xHOST
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GSL_INCL = -I/opt/gsl/gnu/include -I/opt/gnu/gmp/4.3.2/include
GSL_LIBS = -L/opt/gsl/gnu/lib -L/opt/gnu/gmp/4.3.2/lib
FFTW_INCL= -I/opt/fftw/2.1.5/gnu/mvapich2/include
FFTW_LIBS= -L/opt/fftw/2.1.5/gnu/mvapich2/lib
MPICHLIB =
HDF5INCL = -I/opt/hdf5/intel/mvapich2/ib/include -DH5_USE_16_API
HDF5LIB  = -L/opt/hdf5/intel/mvapich2/ib/lib -lhdf5
endif



ifeq ($(SYSTYPE),"bwforcluster")
CC         =  mpicc     # sets the C-compiler
OPTIMIZE   =  -std=c99 -O3 -g -Wall
GSL_INCL   =  -I/home/hd/hd_hd/hd_mm002/libraries/include
GSL_LIBS   =  -L/home/hd/hd_hd/hd_mm002/libraries/lib  -Xlinker -R -Xlinker /home/hd/hd_hd/hd_mm002/libraries/lib
FFTW_INCL  =  -I/home/hd/hd_hd/hd_mm002/libraries/include
FFTW_LIBS  =  -L/home/hd/hd_hd/hd_mm002/libraries/lib  -Xlinker -R -Xlinker /home/hd/hd_hd/hd_mm002/libraries/lib
HDF5INCL   =  -I/home/hd/hd_hd/hd_mm002/libraries/include   -DH5_USE_16_API
HDF5LIB    =  -lhdf5 -L/home/hd/hd_hd/hd_mm002/libraries/lib  -Xlinker -R -Xlinker /home/hd/hd_hd/hd_mm002/libraries/lib 
HWLOC_INCL =  -I/home/hd/hd_hd/hd_mm002/libraries/include
HWLOC_LIBS =  -lhwloc -L/home/hd/hd_hd/hd_mm002/libraries/lib  -Xlinker -R -Xlinker /home/hd/hd_hd/hd_mm002/libraries/lib
endif




ifeq ($(SYSTYPE),"Odin")
CC       =  mpiicc  
OPTIMIZE =   -std=c11 -g -O2
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GMP_INCL = 
GMP_LIBS = 
GSL_INCL = -I$(GSL_INCDIR)
GSL_LIBS = -L$(GSL_LIBDIR)  -Xlinker -R -Xlinker  $(GSL_LIBDIR)
FFTW_INCL= -I$(FFTW_HOME)/include
FFTW_LIBS= -L$(FFTW_HOME)/lib  -Xlinker -R -Xlinker $(FFTW_HOME)/lib
HDF5INCL = -I$(HDF5_HOME)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_HOME)/lib  -lhdf5 -lz  -Xlinker -R -Xlinker $(HDF5_HOME)/lib
HWLOC_INCL= -I/u/vrs/Libs/include 
HWLOC_LIB = -L/u/vrs/Libs/lib -lhwloc    -Xlinker -R -Xlinker /u/vrs/Libs/lib
MPICHLIB =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"DracoAndHydra")
CC       =  mpiicc
OPTIMIZE =  -std=c11 -march=native -g -O3 -ipo
GSL_INCL =  -I$(GSL_HOME)/include/
GSL_LIBS =  -L$(GSL_HOME)/lib/
FFTW_INCL=  -I$(FFTW_HOME)/include/
FFTW_LIBS=  -L$(FFTW_HOME)/lib/
HWLOC_INCL =
HWLOC_LIB = -lhwloc -L/u/dnelson/.local/lib/
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API=1 -I$(HDF5_HOME)/include/
HDF5LIB  =  -lhdf5 -L$(HDF5_HOME)/lib/
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -qopenmp
endif
endif

ifeq ($(SYSTYPE),"Hornet")
CC       =  cc -std=c11 -g -O2  #-xW -ipo -Wall                                                                                                                    OPTIMIZE =
GMP_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
GMP_LIBS = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib  -Xlinker -R -Xlinker /zhome/academic/HLRS/lha/zahrpakm/libs/lib
GSL_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
GSL_LIBS = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib  -Xlinker -R -Xlinker /zhome/academic/HLRS/lha/zahrpakm/libs/lib
HDF5INCL = -I$(HDF5_DIR)/include -DH5_USE_16_API # need to do module load cray-hdf5                                                                                HDF5LIB  = -L$(HDF5_LIB)/lib -lhdf5 -lz
HWLOC_INCL = -I/zhome/academic/HLRS/lha/zahrpakm/libs/include
HWLOC_LIB  = -L/zhome/academic/HLRS/lha/zahrpakm/libs/lib -lhwloc
MPICHLIB =
OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"darwinUK")
CC       =  mpiicc -g -Wall -Wno-unknown-pragmas -std=c99 -DH5Dopen_vers=1 -DH5Dcreate_vers=1 -DH5Gopen_vers=1 -DH5Gcreate_vers=1 -DH5Acreate_vers=1
OPTIMIZE = -O3
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas
FFTW_INCL=
FFTW_LIBS= -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
MPICHLIB =
HDF5INCL =
HDF5LIB  = -lhdf5 -lm -lz
LIBS = $(GSL_LIBS) $(FFTW_LIBS)
endif

ifeq ($(SYSTYPE),"lonestar")
CC	  = mpicc
OPTIMIZE  = -std=c11 -O2 -m64 -Wno-uninitialized
GSL_INCL  = -I$(TACC_GSL_INC)
GSL_LIBS  = -L$(TACC_GSL_LIB)
FFTW_INCL = -I$(TACC_FFTW2_INC)
FFTW_LIBS = -L$(TACC_FFTW2_LIB)
MPICHLIB  =
HDF5INCL  = -I$(TACC_HDF5_INC)
HDF5LIB   = -L$(TACC_HDF5_LIB) -lhdf5 -lz
HWLOC_INCL= -I/home1/00025/tgreif/libs/hwloc/include
HWLOC_LIB = -L/home1/00025/tgreif/libs/hwloc/lib -lhwloc
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL = -I/home1/00025/tgreif/libs/cvode/include
CVODE_LIB  = -L/home1/00025/tgreif/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
HEALPIX_INCL = -I/home1/00025/tgreif/libs/healpix/include
HEALPIX_LIB  = -L/home1/00025/tgreif/libs/healpix/lib -lchealpix
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE += -openmp
else
OPTIMIZE += -Wno-unknown-pragmas
endif
OPT      += -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"aurora")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -g -O3 -Wall -ipo -vec_report0
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS= -L/home/spb/.local/fftw/lib
MPICHLIB = -lmpi
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 
endif

ifeq ($(SYSTYPE),"hecate")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -g -O3 -w1 -ipo -vec_report0 -xHOST
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = -lmpi
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 
endif

ifeq ($(SYSTYPE),"Darwin")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
CC       =  /opt/local/bin/gcc-mp-4.7  $(MPI_COMPILE_FLAGS)      # to replace clang with gcc (mpicc uses clang for some reason)
endif
GSL_INCL = -I/opt/local/include 
GSL_LIBS = -L/opt/local/lib
FFTW_INCL= -I/opt/local/include -I/usr/local/include
FFTW_LIBS= -L/opt/local/lib -I/usr/local/lib
HDF5INCL = -I/opt/local/include -DH5_USE_16_API 
HDF5LIB  = -L/opt/local/lib  -lhdf5 -lz
HWLOC_INCL= -I/opt/local/include 
HWLOC_LIB = -L/opt/local/lib -lhwloc
HYPRE_INCL = -I/Users/varan/codes/libs/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/Users/varan/codes/libs/hypre-2.10.0b/src/hypre/lib -lHYPRE
MPICHLIB = -lmpi
endif

ifeq ($(SYSTYPE),"Darwin-mpich")
CC       = mpicc   # sets the C-compiler
LINKER   = mpicc
OPTIMIZE = -std=c11 -m64 -ggdb -O3 -Wall -Wno-format-security -Wno-unknown-pragmas
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE+= -fopenmp
endif
GSL_INCL = -I/opt/local/include
GSL_LIBS = -L/opt/local/lib
FFTW_INCL= 
FFTW_LIBS= 
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lz
HWLOC_INCL =
HWLOC_LIB = -lhwloc
CUDA_INCL= -I/Developer/NVIDIA/CUDA-5.0/include
CUDA_LIBS= -Xlinker -rpath /Developer/NVIDIA/CUDA-5.0/lib -L/Developer/NVIDIA/CUDA-5.0/lib -lcudart -lnvToolsExt -framework CUDA
NVCC     = /Developer/NVIDIA/CUDA-5.0/bin/nvcc
CUDA_OPTIMIZE = -g -G -O3 -m64 --ptxas-options=-v -Xptxas="-v" --maxrregcount=32 -arch=sm_30 $(filter -I%,$(shell mpicc -show))
endif

# modules for Magny
# module add mvapich2/gcc/64/1.6-qlc

ifeq ($(SYSTYPE),"Magny") 
CC       =   mpicc
FC	   = mpif90
OPTIMIZE =   -std=c11 -O3 -msse3 -g -Wall -m64
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
else
OPTIMIZE +=  -Wno-unknown-pragmas -Wno-unused-function
endif
GMP_INCL =  -I/hits/tap/sw/libs/gmp-5.0.5/include
GMP_LIBS =  -L/hits/tap/sw/libs/gmp-5.0.5/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/gmp-5.0.5/lib
GSL_INCL =  -I/hits/tap/sw/libs/gsl-1.15/include
GSL_LIBS =  -L/hits/tap/sw/libs/gsl-1.15/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/gsl-1.15/lib
FFTW_INCL=  -I/hits/tap/sw/libs/fftw-3.3.4/include
FFTW_LIBS=  -L/hits/tap/sw/libs/fftw-3.3.4/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/fftw-3.3.4/lib
MPICHLIB =
HDF5INCL =  -I/hits/tap/sw/libs/hdf5-1.8.10/include -DH5_USE_16_API
HDF5LIB  =  -L/hits/tap/sw/libs/hdf5-1.8.10/lib -lhdf5  -Xlinker -R -Xlinker /hits/tap/sw/libs/hdf5-1.8.10/lib
HWLOC_INCL=  
HWLOC_LIB =  -lhwloc
NBC_INCL  =  -I/hits/tap/sw/libs/libNBC-1.1.1/include
NBC_LIB   =  -L/hits/tap/sw/libs/libNBC-1.1.1/lib -lnbc
LINKER   = $(CC) -lgfortran #-nofor-main
HYPRE_LIB = -lHYPRE
VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify
#LINKER   = $(FC) #-nofor-main
#OPT      +=  -DNOCALLSOFSYSTEM
#OPT      +=  -DUSE_SSE
endif

 
# modules for Magny-Intel
# module load intel/compiler
# module load mvapich2/intel/64/1.6-qlc

ifeq ($(SYSTYPE),"Magny-Intel") 
CC       =   mpiicc
OPTIMIZE =   -std=c11 -O3 -g -Wall -m64
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GSL_INCL =  -I/hits/tap/sw/libs/include
GSL_LIBS =  -L/hits/tap/sw/libs/lib  -Xlinker -R -Xlinker /hits/tap/sw/libs/lib
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
MATHLIB  = -limf -lm
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5
VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify

#OPT      +=  -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"Haswell")
CC       =   mpicc
FC         = mpif90
OPTIMIZE =   -std=c11 -O3 -msse3 -g -Wall -m64
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
else
OPTIMIZE +=  -Wno-unknown-pragmas -Wno-unused-function
endif
GMP_INCL =  -I/hits/basement/tap/sw/libs/gmp-6.1.1/include
GMP_LIBS =  -L/hits/basement/tap/sw/libs/gmp-6.1.1/lib  -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/gmp-6.1.1/lib
GSL_INCL =  -I/hits/basement/tap/sw/libs/gsl-2.1/include
GSL_LIBS =  -L/hits/basement/tap/sw/libs/gsl-2.1/lib  -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/gsl-2.1/lib
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB =
HDF5INCL =  -I/hits/basement/tap/sw/libs/hdf5-1.8.17/include -DH5_USE_16_API
HDF5LIB  =  -L/hits/basement/tap/sw/libs/hdf5-1.8.17/lib -lhdf5  -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/hdf5-1.8.17/lib
HWLOC_INCL= -I/hits/basement/tap/sw/libs/hwloc-1.11.3/include
HWLOC_LIB = -L/hits/basement/tap/sw/libs/hwloc-1.11.3/lib -lhwloc -Xlinker -R -Xlinker /hits/basement/tap/sw/libs/hwloc-1.11.3/lib
NBC_INCL  =  
NBC_LIB   =  
LINKER   = $(CC) -lgfortran #-nofor-main                                                                                                                          
HYPRE_LIB = -lHYPRE
VTUNE_INCL = -I/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/include
VTUNE_LIBS = -L/cm/shared/apps/intel/vtune_u13/vtune_amplifier_xe_2015/lib64 -littnotify
#LINKER   = $(FC) #-nofor-main                                                                                                                                    
#OPT      +=  -DNOCALLSOFSYSTEM                                                                                                                                   
#OPT      +=  -DUSE_SSE                                                                                                                                           
endif

ifeq ($(SYSTYPE),"OpenSuse")
CC       =  mpicc -g -O2 -Wall
OPTIMIZE =  -std=c11 
GSL_INCL =  
GSL_LIBS =  
FFTW_INCL=  
FFTW_LIBS=  
MPICHLIB = -L/usr/lib/mpi/gcc/openmpi/lib -Xlinker -R -Xlinker /usr/lib/mpi/gcc/openmpi/lib -lmpi 
HDF5INCL = -DH5_USE_16_API
HDF5LIB  =
endif

ifeq ($(SYSTYPE),"OpenSuse64")
CC       = mpicc
OPTIMIZE = -std=c11 -g -O3 -Wall 
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GSL_INCL = -I/sw/tap/include
GSL_LIBS = -L/sw/tap/lib -Xlinker -R -Xlinker /sw/tap/lib 
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB = 
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5 -lhdf5_hl
endif

ifeq ($(SYSTYPE),"OpenSuse64-cuda")
CC       = mpicc
LINKER   = mpicxx
OPTIMIZE = -std=c11 -g -O3 -Wall  
GSL_INCL = -I/sw/tap/include
GSL_LIBS = -L/sw/tap/lib -Xlinker -R -Xlinker /sw/tap/lib
FFTW_INCL= 
FFTW_LIBS=  
MPICHLIB =  
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5
CUDA_INCL= -I/usr/local/cuda/include
CUDA_LIBS= -L/usr/local/cuda/lib64 -lcuda -lcudart
NVCC     = /usr/local/cuda/bin/nvcc
CUDA_OPTIMIZE = -g -G -O0 --ptxas-options=-v -Xptxas="-v" -arch=sm_20 $(filter -I%,$(shell mpicc -show))
endif


ifeq ($(SYSTYPE),"Judge")
CC       =  mpicc
LINKER   =  mpicxx
OPTIMIZE =  -std=c11 -g -O3 -Wall
GSL_INCL =  -I/homeb/zam/baueras/libs/include 
GSL_LIBS =  -L/homeb/zam/baueras/libs/lib  -Xlinker -R -Xlinker /homeb/zam/baueras/libs/lib  
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5
CUDA_INCL= -I/opt/cuda/include/
CUDA_LIBS= -L/opt/cuda/lib64 -Xlinker -R -Xlinker /opt/cuda/lib64 -lcuda -lcudart
NVCC     = nvcc
CUDA_OPTIMIZE = -g -G -O3 --ptxas-options=-v -Xptxas="-v" -arch=sm_20 $(filter -I%,$(shell mpicc -show))
endif


ifeq ($(SYSTYPE),"VIP")
CC       =  mpcc_r -g # -qflttrap=enable:zerodivide:nanq # sets the C-compiler
OPTIMIZE = -qstrict -q64 -qcpluscmt  -O3 -qipa
GSL_INCL = -I/afs/rzg/u/vrs/Libs/regatta/include
GSL_LIBS = -L/afs/rzg/u/vrs/Libs/regatta/lib
FFTW_INCL= -I/afs/rzg/u/vrs/Libs/regatta/include
FFTW_LIBS= -L/afs/rzg/u/vrs/Libs/regatta/lib  -q64 # -qipa
MPICHLIB =
HDF5INCL = -I/afs/rzg/u/vrs/Libs/regatta/include
HDF5LIB  = -L/afs/rzg/u/vrs/Libs/regatta/lib  -lhdf5 -lz
OPT	 += -DNOINLINE
endif

ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unused-function
GSL_INCL =
GSL_LIBS =
FFTW_INCL=  -I/n/home12/mvogelsberger/opt/include/
FFTW_LIBS=  -L/n/home12/mvogelsberger/opt/lib/
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL =  -I/n/home01/tgreif/libs/cvode/include
CVODE_LIB  =  -L/n/home01/tgreif/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
else
OPTIMIZE += -Wno-unknown-pragmas
endif
endif

ifeq ($(SYSTYPE),"odyssey2")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -mprefer-avx128 -march=native
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (MONOTONE_CONDUCTION,$(findstring MONOTONE_CONDUCTION,$(CONFIGVARS)))
HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL =  -I/n/home01/tgreif/libs/cvode/include
CVODE_LIB  =  -L/n/home01/tgreif/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
CFITSIO_INCL = -I/n/home01/tgreif/libs/cfitsio/include
CFITSIO_LIB  = -L/n/home01/tgreif/libs/cfitsio/lib -lcfitsio
HEALPIX_INCL = -I/n/home01/tgreif/libs/healpix/include
HEALPIX_LIB  = -L/n/home01/tgreif/libs/healpix/lib -lchealpix
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"PhilipNUC")
CC = mpicc
CPP      =  mpicxx   # sets the C++-compiler
OPTIMIZE = -std=c11 -O3 -Wall -m64 -g -Wno-error=format-security -Wno-unknown-pragmas
GSL_INCL = -I/usr/include
GSL_LIBS = -L/usr/lib -lgsl -lm
FFTW_INCL= -I/usr/include
FFTW_LIBS= -L/usr/lib -lfftw3
MPICHLIB = -lmpich
HDF5INCL = -I/usr/local/include
HDF5LIB  = -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lz
HWLOC_INCL=  -I/usr/include/hdf5/serial -DH5_USE_16_API
HWLOC_LIB =  -lhwloc -Xlinker -R -Xlinker /usr/lib
CPV      =  $(CC)
LINKER = mpicxx
endif

ifeq ($(SYSTYPE),"odyssey-gcc")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c11 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -mprefer-avx128 -march=native
GSL_INCL =
GSL_LIBS =
FFTW_INCL= -I/n/home13/kannan/fftw-3.3.4/include
FFTW_LIBS= -L/n/home13/kannan/fftw-3.3.4/lib
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (MONOTONE_CONDUCTION,$(findstring MONOTONE_CONDUCTION,$(CONFIGVARS)))
HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif
ifeq (IMPLICIT_OHMIC_DIFFUSION,$(findstring IMPLICIT_OHMIC_DIFFUSION,$(CONFIGVARS)))
HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL =  -I/n/home02/fbecerra/libs/cvode/include
CVODE_LIB  =  -L/n/home02/fbecerra/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
CFITSIO_INCL = -I/n/home01/tgreif/libs/cfitsio/include
CFITSIO_LIB  = -L/n/home01/tgreif/libs/cfitsio/lib -lcfitsio
HEALPIX_INCL = -I/n/home01/tgreif/libs/healpix/include
HEALPIX_LIB  = -L/n/home01/tgreif/libs/healpix/lib -lchealpix
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"odyssey-opteron")
CC       =  mpiicc
OPT      +=
OPTIMIZE =  -std=c11 -O3 -m64 -Wno-uninitialized -Wno-unknown-pragmas
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
MATHLIB  =
HDF5INCL =  -I/n/home12/mvogelsberger/optmine/include
HDF5LIB  =  -L/n/home12/mvogelsberger/optmine/lib -Xlinker -R -Xlinker /n/home12/mvogelsberger/optmine/lib -lhdf5 -lz
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL =  -I/n/home02/fbecerra/libs/cvode/include
CVODE_LIB  =  -L/n/home02/fbecerra/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
HEALPIX_INCL =  -I/n/home01/tgreif/libs/healpix/include
HEALPIX_LIB  =  -L/n/home01/tgreif/libs/healpix/lib -lchealpix
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
endif

ifeq ($(SYSTYPE),"odyssey-intel")
CC       =  mpiicc
OPT      += -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =  -std=c99 -O2 -parallel -ipo -funroll-loops
GSL_INCL =
GSL_LIBS =
FFTW_INCL= -I/n/home03/cpopa/fftw3/include
FFTW_LIBS= -L/n/home03/cpopa/fftw3/lib
MPICHLIB =
HDF5INCL = -DH5_USE_16_API 
HDF5LIB  = -lhdf5 -lz
ifeq (MONOTONE_CONDUCTION,$(findstring MONOTONE_CONDUCTION,$(CONFIGVARS)))
HYPRE_INCL = -I/n/home13/kannan/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/n/home13/kannan/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL =  -I/n/home02/fbecerra/libs/cvode/include
CVODE_LIB  =  -L/n/home02/fbecerra/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE += -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
else
OPTIMIZE += -diag-disable 3180
endif
endif

ifeq ($(SYSTYPE),"Ubuntu")
CC       =  mpicc -O1   # sets the C-compiler
OPTIMIZE = -std=gnu11 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-format-security -Wno-unused-result
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL = -I/usr/include/hdf5/serial/ -DH5_USE_16_API 
HDF5LIB  = -lhdf5  
OPT     += -DX86FIX
HYPRE_INCL = -I/home/kannan/codes/hypre-2.10.0b/src/hypre/include
HYPRE_LIB = -L/home/kannan/codes/hypre-2.10.0b/src/hypre/lib -lHYPRE
endif

ifeq ($(SYSTYPE),"Centos5-Intel")
CC       =  icc   # sets the C-compiler
OPTIMIZE =  -std=c11 -g -O3 -ipo  -mtune=host -mcpu=host -mp
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB = -lmpi
HDF5INCL = 
HDF5LIB  = 
endif

ifeq ($(SYSTYPE),"Centos5-Gnu")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= -L/usr/local/lib
FFTW_LIBS= 
MPICHLIB = -lmpi
HDF5INCL = 
HDF5LIB  = -lhdf5
endif

ifeq ($(SYSTYPE),"stampede")
CC	  = mpicc
OPTIMIZE  = -std=c99 -O2 -xhost -Wno-uninitialized -Wno-unknown-pragmas
GMP_INCL  = 
GMP_LIBS  = 
GSL_INCL  = -I$(TACC_GSL_INC)
GSL_LIBS  = -L$(TACC_GSL_LIB)
FFTW_INCL = -I$(TACC_FFTW3_INC)
FFTW_LIBS = -L$(TACC_FFTW3_LIB)
MPICHLIB  =
HDF5INCL  = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB   = -L$(TACC_HDF5_LIB) -lhdf5 -lz
HWLOC_INCL= 
HWLOC_LIB = 
ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
CVODE_INCL = -I/home1/00025/tgreif/libs/cvode/include
CVODE_LIB  = -L/home1/00025/tgreif/libs/cvode/lib -lsundials_cvode -lsundials_nvecserial -Xlinker -R -Xlinker /home1/00025/tgreif/libs/cvode/lib
endif
ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
OPTIMIZE += -openmp
CFITSIO_INCL = -I/home1/00025/tgreif/libs/cfitsio/include
CFITSIO_LIB  = -L/home1/00025/tgreif/libs/cfitsio/lib -lcfitsio
HEALPIX_INCL = -I/home1/00025/tgreif/libs/healpix/include
HEALPIX_LIB  = -L/home1/00025/tgreif/libs/healpix/lib -lchealpix
ifeq (HEALRAY_MIC,$(findstring HEALRAY_MIC,$(CONFIGVARS)))
#OPTIMIZE += -offload-attribute-target=mic
endif
endif
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE += -openmp
endif
OPT      += -DNOCALLSOFSYSTEM
endif

ifeq ($(SYSTYPE),"stampede_illustris")
CC        = mpicc
OPTIMIZE  = -std=c11 -O2 -m64 -Wno-uninitialized
GMP_INCL  = -I/home1/01637/mvogelsb/libs/gmp/include
GMP_LIBS  = -L/home1/01637/mvogelsb/libs/gmp/lib -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/gmp/lib
GSL_INCL  = -I$(TACC_GSL_INC)
GSL_LIBS  = -L$(TACC_GSL_LIB)
FFTW_INCL = -I$(TACC_FFTW2_INC)
FFTW_LIBS = -L$(TACC_FFTW2_LIB)
MPICHLIB  =
HDF5INCL  = -I/home1/01637/mvogelsb/libs/hdf5/include/
HDF5LIB   = -L/home1/01637/mvogelsb/libs/hdf5/lib/ -lhdf5 -lz -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/hdf5/lib/
HWLOC_INCL= -I/home1/01637/mvogelsb/libs/hwloc/include
HWLOC_LIB = -L/home1/01637/mvogelsb/libs/hwloc/lib -lhwloc -Xlinker -R -Xlinker /home1/01637/mvogelsb/libs/hwloc/lib
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE += -openmp
else
OPTIMIZE +=  -Wno-unknown-pragmas
endif
endif

ifeq ($(SYSTYPE),"ScientificLinux")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O2 -Wno-uninitialized
GSL_INCL =  -L/usr/lib64/
GSL_LIBS = -L/usr/lib64/
FFTW_INCL= -I/home/djmunoz/opt/include
FFTW_LIBS= -L/home/djmunoz/opt/lib
MPICHLIB = -L/home/djmunoz/opt/lib -lmpi
HDF5INCL = -L/home/djmunoz/opt/lib/ -DH5_USE_16_API
HDF5LIB  = -L/home/djmunoz/opt/lib/ -lhdf5
endif


ifeq ($(SYSTYPE),"Ranger_intel" )
AR = xiar
else
ifeq ($(SYSTYPE),"aurora")
AR = xiar
else
ifeq ($(SYSTYPE),"odyssey")
AR = ar
else
AR = ar
endif
endif
endif

ifeq ($(SYSTYPE),"bocksbeutel")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -funroll-all-loops -fprefetch-loop-arrays --param prefetch-latency=300 -ftree-vectorize -march=bdver1 -mprefer-avx128 -mcx16 -msahf -mno-movbe -maes -mpclmul -mpopcnt -mabm -mlwp -mno-fma -mfma4 -mxop -mno-bmi -mno-bmi2 -mno-tbm -mavx -mno-avx2 -msse4.2 -msse4.1 -mlzcnt -mno-rdrnd -mno-f16c -mno-fsgsbase --param l1-cache-size=16 --param l1-cache-line-size=64 --param l2-cache-size=2048 -mtune=bdver1
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB = -lmpi -lhwloc
HDF5INCL = -DH5_USE_16_API=1
HDF5LIB  = -lhdf5
endif

ifeq ($(SYSTYPE),"wuerzburg")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -march=native -ftree-vectorize
GSL_INCL = 
GSL_LIBS = 
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB = -lmpi -lhwloc
HDF5INCL = -DH5_USE_16_API=1
HDF5LIB  = -lhdf5
endif

# module load fftw   needed on SuperMUC
ifeq ($(SYSTYPE),"SuperMuc")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -xAVX -g -O3 -Wall
GSL_INCL =  $(GSL_INC)
GSL_LIBS =  -L/lrz/sys/libraries/gsl/1.16/lib/
FFTW_INCL=  $(FFTW_INC)
FFTW_LIBS=  -L$(FFTW_LIBDIR)
HWLOC_INCL = $(HWLOC_INC)
HWLOC_LIB = -L/lrz/sys/tools/hwloc/1.8.1/lib -lhwloc
MPICHLIB =
HDF5INCL =  $(HDF5_INC) -DH5_USE_16_API=1
HDF5LIB  =  $(HDF5_LIB)  -lhdf5
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE += -diag-disable 3180
endif
endif

ifeq ($(SYSTYPE),"SuperMucIntel")
CC       =  mpicc   # sets the C-compiler
FC       =  mpif90
OPTIMIZE =  -std=c11 -xAVX -g -O3 -Wall
GSL_INCL =  $(GSL_INC)
GSL_LIBS =  -L/lrz/sys/libraries/gsl/1.16/lib/
FFTW_INCL=  $(FFTW_INC)
FFTW_LIBS=  -L$(FFTW_LIBDIR)
HWLOC_INCL = $(HWLOC_INC)
HWLOC_LIB = -L/lrz/sys/tools/hwloc/1.8.1/lib -lhwloc
MPICHLIB =
HDF5INCL =  $(HDF5_INC) -DH5_USE_16_API=1
HDF5LIB  =  -L/lrz/sys/libraries/hdf5/1.8.15/serial/lib -lhdf5_hl_cpp -lhdf5_cpp #-L/lrz/sys/libraries/hdf5/1.8.12/serial_gpfs/lib -lhdf5
LINKER   = $(FC) -nofor-main
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE += -diag-disable 3180
endif
endif


ifeq ($(SYSTYPE),"JuQueen")
CC       =  mpixlc_r -g
OPTIMIZE = -std=c11 -qstrict -q64 -qcpluscmt  -O3 -qipa
GMP_INCL = -I/bgsys/local/gmp/5.0.5/include
GMP_LIBS = -L//bgsys/local/gmp/5.0.5/lib
GSL_INCL = -I/bgsys/local/gsl/1.15_O3g/include
GSL_LIBS = -L/bgsys/local/gsl/1.15_O3g/lib
FFTW_INCL= -I/bgsys/local/fftw2/2.1.5/include
FFTW_LIBS= -L/bgsys/local/fftw2/2.1.5/lib  -q64 # -qipa                                                                                                                                                              
MPICHLIB =
HDF5INCL = -I/bgsys/local/hdf5/include
HDF5LIB  = -L/bgsys/local/hdf5/lib  -lhdf5 -lz
endif

ifeq ($(SYSTYPE),"monstrum")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c11 -O3 -g -Wall #-Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function
GSL_INCL =
GSL_LIBS =
FFTW_INCL= 
FFTW_LIBS=  
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"mira")
CC       =  mpixlc   # sets the C-compiler
OPTIMIZE =  -std=c11 -g -O2 -qarch=qp -qtune=qp
GMP_INCL =  -I/home/mvogelsb/gmp/include
GMP_LIBS =  -L/home/mvogelsb/gmp/lib
GSL_INCL =  -I/soft/libraries/unsupported/gsl/1.9/xl/include 
GSL_LIBS =  -L/soft/libraries/unsupported/gsl/1.9/xl/lib -lgsl
FFTW_INCL=  -I/home/mvogelsb/fftw/include 
FFTW_LIBS=  -L/home/mvogelsb/fftw/lib
HWLOC_INCL=
MPICHLIB =   -Wl,-dy -L/bgsys/drivers/ppcfloor/gnu-linux/powerpc64-bgq-linux/lib -lstdc++
HDF5INCL =  -I/soft/libraries/hdf5/1.8.10/cnk-xl/current/include -DH5_USE_16_API
HDF5LIB  =  -L/soft/libraries/hdf5/1.8.10/cnk-xl/current/lib -lhdf5 -L/soft/libraries/alcf/current/xl/ZLIB/lib -lz
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE +=
endif
endif

ifeq ($(SYSTYPE),"cosma4")
CC       = $(MPIROOT)/bin/mpicc
OPTIMIZE = -std=c11 -DUSE_IRECV -g -O3 -ip -DH5_USE_16_API -Wno-unknown-pragmas
HDF5LIB  = -lhdf5
MPICHLIB = -lmpi
EXTRA_LDFLAGS = -shared-intel $(OMPI_LDFLAGS) # only for intel compiler
endif

ifeq ($(SYSTYPE),"garpur")
CC       = mpicc
OPTIMIZE = -std=c11 -DUSE_IRECV -g -O3 -ip -DH5_USE_16_API -Wno-unknown-pragmas
HDF5LIB  = -lhdf5
MPICHLIB = -lmpi
#I try to use GCC
#EXTRA_LDFLAGS = -shared-intel $(OMPI_LDFLAGS) # only for intel compiler
endif

ifeq ($(SYSTYPE),"cosma5")
CC       = $(MPIROOT)/bin/mpicc
OPTIMIZE = -std=c11 -DUSE_IRECV -g -O3 -ip -DH5_USE_16_API -Wno-unknown-pragmas
HDF5LIB  = -lhdf5
MPICHLIB = -lmpi
EXTRA_LDFLAGS = -shared-intel $(OMPI_LDFLAGS) # only for intel compiler
endif

ifeq ($(SYSTYPE),"yeti")
CC       =   mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
GSL_INCL =
GSL_LIBS =
FFTW_INCL= -I/vega/opt/fftw2-2.1.5/include
FFTW_LIBS= -L/vega/opt/fftw2-2.1.5/lib
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"Flatiron")
CC       =   mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
OPTIMIZE =   -std=c99 -O3 -g -Wall -Wno-unused-but-set-variable -Wno-uninitialized -Wno-unknown-pragmas -Wno-unused-function -march=native
GSL_INCL =
GSL_LIBS =
FFTW_INCL= 
FFTW_LIBS= 
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -fopenmp
OPT      += -DIMPOSE_PINNING -DSOCKETS=4 -DMAX_CORES=16
endif
endif

ifeq ($(SYSTYPE),"OPA")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -ftree-vectorize -march=core2 -mcx16 -msahf -mno-movbe -mno-aes -mno-pclmul -mno-popcnt -mno-abm -mno-lwp -mno-fma -mno-fma4 -mno-xop -mno-bmi -mno-bmi2 -mno-tbm -mno-avx -mno-avx2 -mno-sse4.2 -mno-sse4.1 -mno-lzcnt -mno-rdrnd -mno-f16c -mno-fsgsbase --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=4096 -mtune=core2
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = -lmpi -lhwloc
HDF5INCL = -DH5_USE_16_API=1
HDF5LIB  = -lhdf5
endif

# module load compiler/gnu/5 use.own lib/gmp lib/hdf5/1.8.14-gnu-4.4
ifeq ($(SYSTYPE),"ForHLRI")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -ftree-vectorize -march=native
#OPTIMIZE =  -std=c11 -xHost -g -O3 -Wall
GMP_INCL =  -I/opt/gmp/6.0/include
#GMP_INCL =  -I/pfs/data1/home/kit/scc/fk0832/software/gmp/include
#GMP_LIBS   = -L/pfs/data1/home/kit/scc/fk0832/software/gmp/lib
GMP_LIBS   = -L/opt/gmp/6.0/lib
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = -lmpi
#HDF5INCL = -DH5_USE_16_API=1 -I/opt/bwhpc/common/lib/hdf5/1.8.15-intel-15.0/include
#HDF5LIB  = -L/opt/bwhpc/common/lib/hdf5/1.8.15-intel-15.0/lib -lhdf5
HDF5INCL = -DH5_USE_16_API=1 -I/opt/bwhpc/common/lib/hdf5/1.8.14-gnu-4.4/include
HDF5LIB  = -L/opt/bwhpc/common/lib/hdf5/1.8.14-gnu-4.4/lib -lhdf5
endif

ifeq ($(SYSTYPE),"ARCCA")
MPICHLIB =
CC       =  mpicc -std=c99
FC       =  mpif90 -nofor_main
OPTIMIZE =  -O3
GSL_INCL =  -I/software/libraries/gsl/1.15/gnu-4.1.2/include/
GSL_LIBS =  -L/software/libraries/gsl/1.15/gnu-4.1.2/lib/ -lgsl -lgslcblas -lm
GMP_INCL =  -I/software/libraries/gmp/5.0.2/include/
GMP_LIBS =  -L/software/libraries/gmp/5.0.2/lib/
FFTW_INCL=
FFTW_LIBS=
HDF5INCL =
HDF5LIB  =
HWLOC_INCL= -I/usr/include/
HWLOC_LIB = -L/usr/local/lib64 -lhwloc -Xlinker -R -Xlinker /usr/local/lib64
LINKER   = $(FC)
ifeq (PRIMCHEM,$(findstring PRIMCHEM,$(CONFIGVARS)))
ODE_INCL =  -I/usr/local/include
ODE_LIB  =  -L/usr/local/lib -lsundials_cvode -lsundials_nvecserial
endif
endif

ifeq ($(SYSTYPE),"Raijin")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -g -O3 -ipo -m64 -Wall -xCORE-AVX2
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = -lmpi
HDF5INCL = -DH5_USE_16_API
HDF5LIB  = -lhdf5
endif

ifeq ($(SYSTYPE),"bwFor")
CC        = mpicc
FC        = mpif90
OPTIMIZE  = -std=c11 -O3 -g -xhost
GSL_INCL  = -I$(GSL_INC_DIR)
GSL_LIBS   = -L$(GSL_LIB_DIR)
GMP_INCL  = 
GMP_LIBS  = 
FFTW_INCL = -I$(FFTW_INC_DIR)
FFTW_LIBS = -L$(FFTW_LIB_DIR)
MKL_INCL  = -I$(MKL_LIB_DIR)
MKL_LIBS  = -L$(MKL_INC_DIR)
HDF5INCL   = -I$(HDF5_INC_DIR) -DH5_USE_16_API
HDF5LIB    = -L$(HDF5_LIB_DIR) -lhdf5 -lz
MPICHLIB = -lmpi
LINKER   = $(FC) -nofor-main
endif

ifeq ($(SYSTYPE),"Milkyway")
CC        = mpicc
FC        = mpif90
OPTIMIZE  = -std=c11 -O3 -g
MPICHLIB = -lmpi
LINKER   = $(FC)
endif

ifeq ($(SYSTYPE),"Leibniz")
CC       =  mpiicc
OPTIMIZE =  -std=c11 -g -O3 -Wall
GSL_INCL =  $(GSL_INC)
GSL_LIBS =  $(GSL_LIB)
FFTW_INCL=  $(FFTW_INC)
FFTW_LIBS=  -L$(FFTW_LIBDIR)
HWLOC_INCL = $(HWLOC_INC)
HWLOC_LIB =
MPICHLIB =
HYPRE_INCL = -I/opt/hypre/2.11.2/include
HYPRE_LIB = -L/opt/hypre/2.11.2/lib -lHYPRE
HDF5INCL =  $(HDF5_INC) -DH5_USE_16_API=1
HDF5LIB  =  $(HDF5_LIB)  -lhdf5
ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
OPTIMIZE +=  -openmp
else
OPTIMIZE += -diag-disable 3180
endif
endif

ifeq ($(SYSTYPE),"elgato")
CC       =  mpicc   # sets the C-compiler
OPTIMIZE =  -std=c11 -ggdb -O3 -m64 -Wall -Wno-format-security -Wno-unknown-pragmas -Wno-unused-function
GMP_INCL   = -I$(GMP_INC_DIR)
GMP_LIBS   = -L$(GMP_LIB_DIR) # -Xlinker -R -Xlinker $(GMP_LIB_DIR) 
GSL_INCL   = -I$(GSL_INC_DIR)
GSL_LIBS   = -L$(GSL_LIB_DIR)
#FFTW_INCL  = -I$(FFTW2_INC_DIR)
#FFTW_LIBS  = -L$(FFTW2_LIB_DIR)
HDF5INCL   = -I$(HDF5_INC_DIR) -DH5_USE_16_API
HDF5LIB    = -L$(HDF5_LIB_DIR) -lhdf5 -lz #-Xlinker -R -Xlinker $(HDF5_LIB_DIR)
FFTW_INCL=
FFTW_LIBS=
MPICHLIB = -lmpi
endif

ifeq ($(SYSTYPE),"HI-Garpur")
CC       =  mpicc
OPTIMIZE =  -O3
GSL_INCL =  -I/users/home/boyd/local/include
GSL_LIBS =  -L/users/home/boyd/local/lib #-static
FFTW_INCL=  -I/users/home/boyd/local/include
FFTW_LIBS=  -L/users/home/boyd/local/lib
#MPICHLIB =  -L/opt/local/lib/mpich-mp -lmpich
HDF5INCL =  -I/users/home/boyd/local/include -DH5_USE_16_API
HDF5LIB  =  -L/users/home/boyd/local/lib  -lhdf5 #-lz
#OPT     += -DNOTYPEPREFIX_FFTW
endif

ifndef LINKER
LINKER = $(CC)
endif

ifeq ($(SYSTYPE),"LAPTOP")
CC       =  mpicc #-arch x86_64  
OPTIMIZE =  -O3    
GSL_INCL =  -I/usr/local/include
GSL_LIBS =  -L/usr/local/lib #-static
FFTW_INCL=  -I/usr/local/include
FFTW_LIBS=  -L/usr/local/lib 
#MPICHLIB =  -L/opt/local/lib/mpich-mp -lmpich 
HDF5INCL =  -I/home/jzavala/HDF5/hdf5-1.10.0-patch1/hdf5/include -DH5_USE_16_API
HDF5LIB  =  -L/home/jzavala/HDF5/hdf5-1.10.0-patch1/hdf5/lib  -lhdf5 #-lz 
#OPT     += -DNOTYPEPREFIX_FFTW
endif


##########################################
#determine the needed object/header files#
##########################################

SUBDIRS = . mpi_utils pm gitversion

OBJS = powerspec_vel.o  \
         forcetree.o forcetree_walk.o forcetree_ewald.o  \
         ngbtree.o ngbtree_walk.o forcetree_optimizebalance.o \
         domain.o domain_balance.o domain_box.o domain_counttogo.o \
         domain_DC_update.o domain_exchange.o domain_toplevel.o \
         domain_rearrange.o domain_sort_kernels.o domain_vars.o domain_checks.o	\
         voronoi_dynamic_update.o voronoi_ghost_search.o logs.o \
         timestep_treebased.o	\
         external_disk.o riemann.o riemann_rosunov.o riemann_hll.o	\
         riemann_hllc.o riemann_hlld.o riemann_gamma.o finite_volume_solver.o 	\
         set_vertex_velocities.o do_gravity_hydro.o voronoi.o voronoi_utils.o	\
         voronoi_3d.o voronoi_1d.o voronoi_1d_spherical.o \
         voronoi_exchange.o voronoi_check.o \
         voronoi_makeimage.o voronoi_makeimage_new.o coffee.o		\
         windtunnel.o growing_disk_potential.o	voronoi_3d_test.o \
         allvars.o gradients.o scalars.o hdf5_util.o  	\
         pinning.o mpz_extension.o run.o predict.o begrun.o	\
         global.o pm_periodic2d.o timestep.o init.o restart.o io.o io_fields.o \
         starformation.o mpi_utils/checksummed_sendrecv.o accel.o read_ic.o 	\
         parallel_sort.o second_order.o system.o allocate.o mpi_utils/myalltoall.o 	\
         density.o noh.o mpi_utils/sizelimited_sendrecv.o mpi_utils/hypercube_allgatherv.o	\
         gravtree.o gravdirect.o grav_softening.o grav_external.o gravtree_forcetest.o driftfac.o darkenergy.o 	\
         peano.o pm/pm_periodic.o pm/pm_mpi_fft.o pm/pm_nonperiodic.o longrange.o	\
         mymalloc.o helm_eos.o \
         second_derivatives.o\
         diffusion_fluxes.o extrapolate_quantities.o boundaries.o	\
         tracer_particle.o tracer_mc.o tracer_trajectory.o debug.o  \
         update_primitive_variables.o voronoi_derefinement.o		\
         voronoi_refinement.o voronoi_gradients.o			\
         refinement.o criterion_refinement.o parallel_logs.o	\
         criterion_derefinement.o mpi_utils/mpi_util.o                 \
         relaxobject.o voronoi_proj.o voronoi_derefinement_pairs.o      \
         sfr_eEOS.o sfr_quicklyalpha.o parameters.o voronoi_gradients_onedims.o \
         twopoint.o diffusion_general.o runge_kutta_full.o inspiral.o \
         ngbtree_search.o voronoi_gradients_lsf.o main.o

INCL += allvars.h proto.h forcetree.h domain.h dd.h dtypes.h	\
        voronoi.h mesh.h helm_eos.h voronoi_proj.h gitversion/version.h\
        timer.h timestep.h generic_comm_helpers2.h  generic_comm_helpers_async.h\
        runge_kutta_full.h parallel_logs.h


ifeq (VORONOI,$(findstring VORONOI,$(CONFIGVARS)))
ifeq (TWODIMS,$(findstring TWODIMS,$(CONFIGVARS)))
OBJS    += voronoi_2d.o
endif
endif

ifeq (NEW_FFT,$(findstring NEW_FFT,$(CONFIGVARS)))
OBJS    += my_fft/my_fft.o
SUBDIRS += my_fft
endif

ifeq (MYIBARRIER,$(findstring MYIBARRIER,$(CONFIGVARS)))
OBJS    += mpi_utils/myIBarrier.o
INCL    += mpi_utils/myIBarrier.h
endif

ifeq (MHD,$(findstring MHD,$(CONFIGVARS)))
OBJS    += mhd.o
endif

ifeq (ADDBACKGROUNDGRID,$(findstring ADDBACKGROUNDGRID,$(CONFIGVARS)))
OBJS    += add_backgroundgrid/add_bggrid.o add_backgroundgrid/calc_weights.o add_backgroundgrid/distribute.o
INCL    += add_backgroundgrid/add_bggrid.h
SUBDIRS += add_backgroundgrid
endif

ifeq (DVR_RENDER,$(findstring DVR_RENDER,$(CONFIGVARS)))
OBJS    += dvr_render/dvr_render.o
INCL    += dvr_render/dvr_render.h
SUBDIRS += dvr_render
endif

ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling/cooling.o cooling/simple_cooling.o
INCL    += cooling/cooling_vars.h cooling/cooling_proto.h
SUBDIRS += cooling 
endif

ifeq (RT_ENABLE,$(filter RT_ENABLE,$(CONFIGVARS)))            
OBJS    += rt/rt_CGmethod.o rt/rt_gradients.o rt/rt_voronoi.o  \
           rt/rt_stellar_sources.o rt/rt_advect.o rt/rt_cooling.o rt/rt_chem.o	\
           rt/rt_inject_photons_sfr.o rt/rt_inject_photons.o\
           rt/rt_optical_depth.o rt/pix2vec_ring.o
SUBDIRS += rt
endif

ifeq (BLACK_HOLES,$(findstring BLACK_HOLES,$(CONFIGVARS)))     
OBJS    += blackhole/blackhole_neighbors.o blackhole/blackhole_mergers.o blackhole/blackhole_bubbles.o blackhole/blackhole_adios_wind.o \
           blackhole/blackhole_swallowgas.o blackhole/blackhole.o blackhole/blackhole_density.o blackhole/blackhole_adios_wind_randomized.o \
           blackhole/blackhole_disk_vorticity.o blackhole/blackhole_bubbles_nf.o blackhole/blackhole_friction.o \
           blackhole/blackhole_refinement.o blackhole/blackhole_mdot.o  blackhole/blackhole_feedback.o  blackhole/blackhole_centering.o \
		   blackhole/blackhole_bipolar.o blackhole/blackhole_spin.o
INCL    += blackhole/blackhole_proto.h
SUBDIRS += blackhole
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += fof/fof.o fof/fof_vars.o fof/fof_distribute.o fof/fof_findgroups.o fof/fof_nearest.o fof/fof_io.o fof/fof_sort_kernels.o \
           fof/fof_fuzz.o fof/fof_gfm.o fof/fof_bh.o fof/fof_spinmeasurement.o fof/fof_massiveseeds.o                
INCL    += fof/fof.h
SUBDIRS += fof
endif

ifeq (GFM_AGN_RADIATION,$(findstring GFM_AGN_RADIATION,$(CONFIGVARS)))
OBJS    += GFM/agn_radiation.o
SUBDIRS += GFM
endif

ifeq (GFM_STELLAR_PHOTOMETRICS,$(findstring GFM_STELLAR_PHOTOMETRICS,$(CONFIGVARS)))
OBJS    += GFM/stellar_photometrics.o GFM/stellar_photometrics_vars.o
INCL    += GFM/stellar_photometrics_vars.h GFM/stellar_photometrics_proto.h
SUBDIRS += GFM
endif

ifeq (GFM_STELLAR_EVOLUTION,$(findstring GFM_STELLAR_EVOLUTION,$(CONFIGVARS)))
OBJS    += GFM/stellar_evolution_init.o GFM/stellar_evolution_main.o GFM/stellar_evolution_evolve.o \
           GFM/stellar_evolution_util.o GFM/stellar_evolution_vars.o GFM/stellar_density.o          \
           GFM/stellar_evolution_logs.o GFM/stellar_evolution_dust.o GFM/stellar_evolution_enrich.o
INCL    += GFM/stellar_evolution_proto.h GFM/stellar_evolution_vars.h
SUBDIRS += GFM
endif

ifeq (GFM_COOLING_METAL,$(findstring GFM_COOLING_METAL,$(CONFIGVARS)))
OBJS    += GFM/cooling_metal.o GFM/cooling_metal_vars.o
INCL    += GFM/cooling_metal_proto.h GFM/cooling_metal_vars.h
SUBDIRS += GFM
endif

ifeq (GFM_DUST_COOLING,$(findstring GFM_DUST_COOLING,$(CONFIGVARS)))
OBJS    += GFM/cooling_dust.o 
INCL    += GFM/cooling_dust_proto.h 
SUBDIRS += GFM
endif

ifeq (GFM_WINDS,$(findstring GFM_WINDS,$(CONFIGVARS)))
OBJS    += GFM/winds.o GFM/winds_variable.o GFM/winds_vars.o  GFM/winds_local.o GFM/winds_recouple.o GFM/winds_findcells.o
INCL    += GFM/winds_proto.h 
SUBDIRS += GFM
endif

ifeq (GFM,$(findstring GFM,$(CONFIGVARS)))
OBJS    += GFM/helper.o 
INCL    += GFM/helper_proto.h
SUBDIRS += GFM
endif

ifeq (DMPIC,$(findstring DMPIC,$(CONFIGVARS)))
OBJS    += dmpic/project.o
INCL    += dmpic/dmpic.h
SUBDIRS += dmpic
endif

ifeq (ISM,$(findstring ISM,$(CONFIGVARS)))
OBJS    += ism/ism_sfr.o
INCL    += ism/ism_proto.h
SUBDIRS += ism 
endif

ifeq (ISM_LOCAL_RADIATION_PRESSURE,$(findstring ISM_LOCAL_RADIATION_PRESSURE,$(CONFIGVARS)))
OBJS    += ism/ism_local_radiation_pressure.o ism/ism_neighbors.o
INCL    += ism/ism_local_radiation_pressure_proto.h
SUBDIRS += ism
endif

ifeq (ISM_HII_HEATING,$(findstring ISM_HII_HEATING,$(CONFIGVARS)))
OBJS   += ism/ism_HII_heating.o fm_star_formation/radiation_stellar_feedback_util.o
INCL   += ism/ism_HII_heating.h fm_star_formation/radiation_stellar_feedback_proto.h
SUBDIRS += ism fm_star_formation 
endif

inc_stellar_feedback = 
ifeq (GFM_STELLAR_FEEDBACK,$(findstring GFM_STELLAR_FEEDBACK,$(CONFIGVARS)))
inc_stellar_feedback = yes
endif

ifeq (GFM_WINDS_LOCAL,$(findstring GFM_WINDS_LOCAL,$(CONFIGVARS)))
inc_stellar_feedback = yes
endif

ifdef inc_stellar_feedback
OBJS    += GFM/stellar_feedback.o GFM/stellar_feedback_kernels.o
INCL    += GFM/stellar_feedback_kernels.h
SUBDIRS += GFM
endif


ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o subfind/subfind_vars.o  subfind/subfind_serial.o  subfind/subfind_coll_tree.o \
           subfind/subfind_properties.o subfind/subfind_so.o  subfind/subfind_distribute.o \
           subfind/subfind_collective.o subfind/subfind_findlinkngb.o subfind/subfind_nearesttwo.o \
           subfind/subfind_loctree.o subfind/subfind_coll_domain.o  subfind/subfind_coll_treewalk.o \
           subfind/subfind_density.o subfind/subfind_io.o subfind/subfind_sort_kernels.o \
           subfind/subfind_reprocess.o subfind/subfind_so_potegy.o
INCL	+= subfind/subfind.h
SUBDIRS += subfind
endif

ifeq (AMR,$(findstring AMR,$(CONFIGVARS)))
OBJS	+= amr/amr.o amr/amr_mesh.o amr/amr_refinement.o amr/amr_exchange.o amr/amr_ngb.o amr/amr_refinement_criterion.o amr/amr_validate.o amr/amr_remap.o amr/amr_gradients.o amr/amr_generate_gas.o amr/amr_update_nodes.o amr/amr_walk.o
#OBJS	+= amr/amr_walk.o  amr/amr_search.o
ifeq (TWODIMS,$(findstring TWODIMS,$(CONFIGVARS)))
OBJS    += amr/amr_2d.o
endif
INCL	+= amr/amr.h amr/amr_proto.h
SUBDIRS += amr
endif

ifeq (DG,$(filter DG,$(CONFIGVARS)))
OBJS	+= dg/dg_vars.o dg/dg_legendre.o dg/dg_core.o dg/dg_debug.o dg/dg_set_get.o dg/dg_limiter.o dg/dg_fluxes.o dg/dg_recompute.o dg/dg_projection.o dg/dg_time_integration.o dg/dg_limiter_special.o
OBJS    += dg/dg_io.o dg/dg_refinement.o
INCL	+= dg/dg_vars.h dg/dg_proto.h dg/dg_defines.h dg/dg_core_inline.h
SUBDIRS += dg
EXEC = Tenet
endif

ifeq (DG_TEST_PROBLEM,$(filter DG_TEST_PROBLEM,$(CONFIGVARS)))
OBJS	+= dg/dg_test_problems.o dg/cell_projection.o
INCL	+= dg/dg_test_problems.h 
SUBDIRS += dg
endif

ifeq (VS_TURB,$(findstring VS_TURB,$(CONFIGVARS)))
OBJS	+= turb/turb_driving.o turb/turb_powerspectra.o
SUBDIRS += turb
endif

ifeq (AB_TURB,$(findstring AB_TURB,$(CONFIGVARS)))
OBJS	+= turb/ab_turb.o  turb/turb_driving.o turb/turb_powerspectra.o
SUBDIRS += turb
endif

ifeq (MHD_SEEDPSPEC,$(findstring MHD_SEEDPSPEC,$(CONFIGVARS)))
OBJS	+= constrained_transport/mhd_seedpspec.o
SUBDIRS += constrained_transport
endif

ifeq (REGULARIZE_MESH_LLOYD,$(findstring REGULARIZE_MESH_LLOYD,$(CONFIGVARS)))
OBJS	+= constrained_transport/regularize_mesh_more_options.o
SUBDIRS += constrained_transport
endif

ifeq (REGULARIZE_MESH_SMOOTH,$(findstring REGULARIZE_MESH_SMOOTH,$(CONFIGVARS)))
OBJS	+= constrained_transport/regularize_mesh_more_options.o
SUBDIRS += constrained_transport
endif

ifeq (TGSET,$(findstring TGSET,$(CONFIGVARS)))
OBJS    += tgset/tgset.o
INCL	+= tgset/tgset.h tgset/tgset_proto.h
SUBDIRS += tgset
endif

ifeq (TGCHEM,$(findstring TGCHEM,$(CONFIGVARS)))
OBJS    += tgchem/tgchem.o tgchem/tgchem_chem.o tgchem/tgchem_cool.o tgchem/tgchem_init.o tgchem/tgchem_photo.o tgchem/tgchem_step.o tgchem/tgchem_utils.o
INCL	+= tgchem/tgchem.h tgchem/tgchem_proto.h
SUBDIRS	+= tgchem
endif

ifeq (HEALRAY,$(findstring HEALRAY,$(CONFIGVARS)))
OBJS	+= healray/healray.o healray/healray_comm.o healray/healray_experimental.o healray/healray_finish.o healray/healray_init.o healray/healray_rayout.o healray/healray_sources.o
INCL	+= healray/healray.h healray/healray_proto.h
SUBDIRS += healray
endif

ifeq (SGCHEM,$(findstring SGCHEM,$(CONFIGVARS)))
OBJS   += SGChem/calc_shield.o SGChem/cheminmo.o SGChem/cma-setup.o SGChem/cma-util.o \
          SGChem/const_rates.o SGChem/cool_func.o SGChem/cool_util.o \
          SGChem/coolinmo.o SGChem/dvode.o SGChem/evolve_abundances.o SGChem/init_chemistry_parameters.o \
          SGChem/jac.o SGChem/photoinit_ism.o SGChem/photoinit_lowZ.o SGChem/rate_eq_simple.o SGChem/rate_eq_primordial.o \
          SGChem/rate_eq_nl99.o SGChem/sgchem.o SGChem/spline.o SGChem/validate_rates.o SGChem/set_local_abundances.o SGChem/lwbg.o
INCL   += SGChem/cma.h SGChem/cool.h SGChem/f2c.h SGChem/fs_data.h SGChem/h2heat.h  SGChem/isrf.h \
          SGChem/mol_data.h SGChem/non_eq.h SGChem/sgchem_def.h SGChem/sgchem_proto.h SGChem/shield_data.h
SUBDIRS += SGChem
endif

ifeq (TREE_RAD,$(findstring TREE_RAD,$(CONFIGVARS)))
OBJS   += TreeCol/healpix_types.o TreeCol/bit_manipulation.o TreeCol/healpix.o TreeCol/project_column.o
INCL   += TreeCol/treecol_proto.h
SUBDIRS += TreeCol
endif

ifeq (SNE_FEEDBACK,$(findstring SNE_FEEDBACK,$(CONFIGVARS)))
OBJS   += sne/sne.o sne/sne_injection_criteria.o sne/sne_utility.o sne/sne_time_stepping.o
INCL   += sne/sne.h sne/sne_proto.h 
SUBDIRS += sne
endif

ifeq (SINK_PARTICLES,$(findstring SINK_PARTICLES,$(CONFIGVARS)))
OBJS    += sink_particles/accrete_onto_sink_particles.o sink_particles/dump_sink_particle_info.o \
           sink_particles/init_sink_particles.o  sink_particles/sink_particles.o \
           sink_particles/create_sink_particles.o sink_particles/get_all_sink_particle_info.o
INCL    += sink_particles/proto_sink_particles.h sink_particles/sink_particles.h
SUBDIRS += sink_particles
endif

ifeq (SINKS,$(findstring SINKS,$(CONFIGVARS)))
OBJS	+= sinks/sinks.o sinks/sinks_accrete.o sinks/sinks_aux.o sinks/sinks_create.o sinks/sinks_init.o sinks/sinks_dmass.o
INCL	+= sinks/sinks.h sinks/sinks_proto.h
SUBDIRS += sinks
endif

ifeq (SIDM,$(findstring SIDM,$(CONFIGVARS)))
OBJS    += sidm/sidm_hsml.o sidm/sidm_ngb.o sidm/sidm_scatter.o sidm/sidm_cross.o sidm/sidm.o sidm/sidm_vars.o sidm/sidm_scatter_process.o
INCL    += sidm/sidm_vars.h sidm/sidm_proto.h
SUBDIRS += sidm
endif

ifeq (DISSIPATIVE,$(findstring DISSIPATIVE,$(CONFIGVARS)))
OBJS    += dissipative/dissipative_eos.o dissipative/dissipative_cooling.o
INCL    += dissipative/dissipative_proto.h
SUBDIRS += dissipative
endif

ifeq (ADJ_BOX_POWERSPEC,$(findstring ADJ_BOX_POWERSPEC,$(CONFIGVARS)))
OBJS    += power_spec/adj_box_powerspec.o
INCL    += power_spec/adj_box_powerspec_proto.h
SUBDIRS += power_spec
endif

ifeq (FM_STAR_FEEDBACK,$(findstring FM_STAR_FEEDBACK,$(CONFIGVARS)))
OBJS    += GFM/stellar_feedback_kernels.o GFM/stellar_feedback.o GFM/stellar_feedback_find_cells.o fm_star_formation/stellar_feedback_util.o
INCL    += GFM/stellar_feedback_kernels.h fm_star_formation/stellar_feedback_proto.h 
SUBDIRS += fm_star_formation GFM
endif

ifeq (FM_SFR,$(findstring FM_SFR,$(CONFIGVARS)))
OBJS    += fm_star_formation/sfr.o 
SUBDIRS += fm_star_formation
endif

ifeq (LOCAL_FEEDBACK,$(findstring LOCAL_FEEDBACK,$(CONFIGVARS)))
OBJS    += local_feedback/compute_sfr.o local_feedback/inject_feedback.o local_feedback/create_star_particles.o
SUBDIRS += local_feedback
endif

ifeq (COMPUTE_SFR_FROM_H2,$(findstring COMPUTE_SFR_FROM_H2,$(CONFIGVARS)))
OBJS    += fm_star_formation/H2_frac.o
INCL    += fm_star_formation/H2_frac_proto.h
SUBDIRS += fm_star_formation
endif

ifeq (TEST_SFR,$(findstring TEST_SFR,$(CONFIGVARS)))
OBJS    += fm_star_formation/test_sfr.o
SUBDIRS += fm_star_formation
endif

ifeq (OTVET,$(findstring OTVET,$(CONFIGVARS)))
####OBJS   += OTVET/do_otvet.o OTVET/otvet_eddington.o OTVET/otvet_star_lum.o OTVET/otvet_CGmethod.o OTVET/otvet_chem.o OTVET/otvet_cooling.o OTVET/otvet_stardensity.o 
OBJS   += OTVET/do_otvet.o OTVET/otvet_eddington.o OTVET/otvet_star_lum.o OTVET/otvet_CGmethod.o OTVET/otvet_chem_ps2009.o OTVET/otvet_chem_ps2011.o OTVET/otvet_cooling.o 
INCL   += OTVET/otvet_proto.h
SUBDIRS += OTVET
endif

ifeq (RADPRESS_OPT_THIN,$(findstring RADPRESS_OPT_THIN,$(CONFIGVARS)))
OBJS   += fm_star_formation/radpressthin.o 
INCL   += fm_star_formation/radpressthin_proto.h
SUBDIRS += fm_star_formation 
endif

ifeq (RADPRESS_OPT_THICK,$(findstring RADPRESS_OPT_THICK,$(CONFIGVARS)))
OBJS   += fm_star_formation/radpressthick.o 
INCL   += fm_star_formation/radpressthick_proto.h
SUBDIRS += fm_star_formation 
endif

ifeq (FM_RADIATION_FEEDBACK,$(findstring FM_RADIATION_FEEDBACK,$(CONFIGVARS)))
OBJS   += GFM/radiation_stellar_feedback_find_cells.o GFM/radiation_stellar_feedback.o fm_star_formation/radiation_stellar_feedback_util.o
INCL   += fm_star_formation/radiation_stellar_feedback_proto.h
SUBDIRS += fm_star_formation GFM 
endif

ifeq (FM_EARLY_STAR_FEEDBACK,$(findstring FM_EARLY_STAR_FEEDBACK,$(CONFIGVARS)))
OBJS   += fm_star_formation/early_stellar_feedback.o fm_star_formation/early_stellar_feedback_util.o
INCL   += fm_star_formation/early_stellar_feedback_proto.h
SUBDIRS += fm_star_formation  
endif


OBJS    += debug_md5/Md5.o  debug_md5/calc_checksum.o
INCL    += debug_md5/Md5.h
SUBDIRS += debug_md5


ifeq (TEST_COOLING_METAL,$(findstring TEST_COOLING_METAL,$(CONFIGVARS)))
OBJS    += Test_Cooling_Metal/test_cooling_metal.o
INCL    += Test_Cooling_Metal/test_cooling_metal_proto.h
SUBDIRS += Test_Cooling_Metal
endif

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS    += network/network.o network/network_solver.o network/utilities.o network/integrate.o network/network_nse.o
INCL    += network/network.h network/network_solver.h network/utilities.h network/integrate.h network/network_nse.h
SUBDIRS += network
endif

ifeq (MRT,$(findstring MRT,$(CONFIGVARS)))
OBJS += MRT/RT.o MRT/RT_set_VET.o MRT/RT_cooling.o MRT/RT_chem_ps2011.o MRT/RT_chem_ps2009.o MRT/RT_riemann_rosunov.o \
        MRT/RT_IR.o MRT/RT_setup.o MRT/RT_riemann_HLLE.o MRT/RT_comoving.o MRT/RT_add_sources_stellar.o \
	MRT/RT_add_sources_blackholes.o MRT/RT_add_radiation_stellar.o MRT/RT_add_radiation_blackholes.o \
        MRT/RT_source_utils.o MRT/RT_finite_volume_solver.o MRT/RT_gradients.o MRT/RT_gradients_lsf.o MRT/RT_exchange.o \
	MRT/RT_run.o MRT/RT_update_primitive_variables.o
INCL += MRT/RT.h MRT/RT_proto.h
SUBDIRS += MRT
endif

ifeq (MONOTONE_CONDUCTION,$(findstring MONOTONE_CONDUCTION,$(CONFIGVARS)))                # Conduction module 
ifeq (IMPLICIT_TI,$(findstring IMPLICIT_TI,$(CONFIGVARS)))
ifeq (SEMI_IMPLICIT_TI,$(findstring SEMI_IMPLICIT_TI,$(CONFIGVARS)))
OBJS     += conduction/conduction_semi_HYPRE.o
else
OBJS     += conduction/conduction_HYPRE.o
endif
else
OBJS     += conduction/conduction.o
endif
INCL     += conduction/conduction.h
SUBDIRS  += conduction
endif

ifeq (CALCULATE_QUANTITIES_IN_POSTPROCESS,$(findstring CALCULATE_QUANTITIES_IN_POSTPROCESS,$(CONFIGVARS)))       
OBJS     += calculate_quantities_in_postprocess/calculate_quantities.o
INCL     += calculate_quantities_in_postprocess/calculate_quantities.h
SUBDIRS  += calculate_quantities_in_postprocess
endif


ifeq (NETWORK_PARDISO,$(findstring NETWORK_PARDISO,$(CONFIGVARS)))
MKL_INCL =  -I/cm/shared/apps/intel/composer_xe/2011_sp1.12.361/mkl/include
MKL_LIBS =  -L/cm/shared/apps/intel/composer_xe/2011_sp1.12.361/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
else
MKL_INCL =
MKL_LIBS =
endif

ifeq (CIRCUMSTELLAR,$(findstring CIRCUMSTELLAR,$(CONFIGVARS)))
OBJS    += circumstellar/circumstellar.o  circumstellar/circumstellar_sinks.o circumstellar/circumstellar_gravity.o circumstellar/circumstellar_kepler.o
INCL    += circumstellar/circumstellar_proto.h
SUBDIRS += circumstellar
endif

ifeq (ATOMIC_DM,$(findstring ATOMIC_DM,$(CONFIGVARS)))
OBJS    += atomic_dm/cooling_atomic_dm.o 
INCL    += atomic_dm/cooling_atomic_dm_proto.h atomic_dm/cooling_atomic_dm_vars.h 
SUBDIRS += atomic_dm
endif

ifeq (FLD,$(filter FLD,$(CONFIGVARS)))
OBJS    += fld/fld.o fld/fld_MGmethod.o fld/fld_HYPRE.o fld/fld_HYPRE_IJ.o
INCL    += fld/fld_proto.h fld/fld.h
SUBDIRS += fld
endif


ifeq (CUDA,$(findstring CUDA,$(CONFIGVARS)))
OBJS    += cuda_util.o
INCL    += cuda_util.h
endif

ifeq (SHOCK_FINDER,$(findstring SHOCK_FINDER,$(CONFIGVARS)))     
OBJS    += shock_finder/shock_finder.o shock_finder/shock_finder_ryu.o shock_finder/shock_finder_skillman.o shock_finder/shock_finder_arepo.o 
INCL    += shock_finder/shock_finder.h shock_finder/shock_finder_fields.h shock_finder/shock_finder_rays.h 
SUBDIRS += shock_finder 
endif

ifeq (MHD_CT,$(findstring MHD_CT,$(CONFIGVARS)))     
OBJS    += constrained_transport/constrained_transport.o 
INCL    += constrained_transport/constrained_transport.h
SUBDIRS += constrained_transport 
endif

ifeq (BECDM,$(findstring BECDM,$(CONFIGVARS)))     
OBJS    += becdm/becdm.o 
INCL    += becdm/becdm.h
SUBDIRS += becdm 
endif

ifeq (COSMIC_RAYS,$(findstring COSMIC_RAYS,$(CONFIGVARS)))
OBJS    += cosmic_rays/cosmic_rays.o
INCL    += cosmic_rays/cosmic_rays.h
SUBDIRS += cosmic_rays
endif

ifeq (COSMIC_RAYS_SN_INJECTION,$(findstring COSMIC_RAYS_SN_INJECTION,$(CONFIGVARS)))
OBJS    += cosmic_rays/cosmic_rays_injection.o cosmic_rays/cosmic_rays_find_ngbs.o
endif

ifeq (COSMIC_RAYS_DIFFUSION,$(findstring COSMIC_RAYS_DIFFUSION,$(CONFIGVARS)))
OBJS    += cosmic_rays/cosmic_rays_diffusion.o cosmic_rays/cosmic_rays_set_diffusion_coefficients.o
endif

ifeq (COSMIC_RAYS_STREAMING,$(findstring COSMIC_RAYS_STREAMING,$(CONFIGVARS)))
OBJS    += cosmic_rays/cosmic_rays_streaming.o
endif

ifeq (COSMIC_RAYS_SHOCK_ACCELERATION,$(findstring COSMIC_RAYS_SHOCK_ACCELERATION,$(CONFIGVARS)))
OBJS    += cosmic_rays/cosmic_rays_shock_acceleration.o
endif

ifeq (NON_IDEAL_MHD,$(findstring NON_IDEAL_MHD,$(CONFIGVARS)))
OBJS    += mhd_nonideal/ohmic_diffusion.o
OBJS    += mhd_nonideal/ambipolar_diffusion.o
INCL    += mhd_nonideal/proto_non_ideal_mhd.h
SUBDIRS += mhd_nonideal
endif

ifeq (IMPLICIT_OHMIC_DIFFUSION,$(findstring IMPLICIT_OHMIC_DIFFUSION,$(CONFIGVARS)))
ifeq (MHD_CT,$(findstring MHD_CT,$(CONFIGVARS)))
ifeq (OHM_CRANK_NICHOLSON,$(findstring OHM_CRANK_NICHOLSON,$(CONFIGVARS)))
OBJS    += mhd_nonideal/implicit_ohmic_diffusion_ct_cn.o
else
OBJS    += mhd_nonideal/implicit_ohmic_diffusion_ct.o
endif
else
ifeq (OHM_CRANK_NICHOLSON,$(findstring OHM_CRANK_NICHOLSON,$(CONFIGVARS)))
OBJS    += mhd_nonideal/implicit_ohmic_diffusion_cn.o
else
OBJS    += mhd_nonideal/implicit_ohmic_diffusion.o
endif
endif
INCL    += mhd_nonideal/proto_non_ideal_mhd.h
SUBDIRS += mhd_nonideal
endif

ifeq (EOS_OPAL,$(findstring EOS_OPAL,$(CONFIGVARS)))
OBJS    += opal_eos.o
INCL    += opal_eos.h
GSL_LIBS +=  -lz
endif

ifeq (SPECIAL_RELATIVITY,$(findstring SPECIAL_RELATIVITY,$(CONFIGVARS)))
OBJS    += special_relativity.o
ifeq (SPECIAL_RELATIVITY_HLLC,$(findstring SPECIAL_RELATIVITY_HLLC,$(CONFIGVARS)))
OBJS    += riemann_hllc_special_relativity.o
else
OBJS    += riemann_hlle_special_relativity.o
endif
endif

ifeq (GENERAL_RELATIVITY,$(findstring GENERAL_RELATIVITY,$(CONFIGVARS)))
OBJS    += general_relativity.o
ifeq (SPECIAL_RELATIVITY_HLLC,$(findstring SPECIAL_RELATIVITY_HLLC,$(CONFIGVARS)))
OBJS    += riemann_hllc_general_relativity.o
else
OBJS    += riemann_hlle_general_relativity.o
endif
endif

ifeq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
OBJS    += grackle/grackle.o
INCL    += grackle/grackle_def.h
SUBDIRS += grackle
endif

ifeq (GRAVITY_TABLE,$(findstring GRAVITY_TABLE,$(CONFIGVARS)))
OBJS    += gravity_table.o
endif

ifeq (MODGRAV,$(findstring MODGRAV,$(CONFIGVARS)))
OBJS    += modgrav/modgrav_forcetree.o modgrav/modgrav_pm_periodic.o  modgrav/modgrav_pm_nonperiodic.o 
INCL    += modgrav/modgrav_forcetree.h modgrav/modgrav_pm.h
SUBDIRS += modgrav
endif

ifeq (DUST_LIVE,$(findstring DUST_LIVE,$(CONFIGVARS)))
OBJS    += dust_live/drag_kernels.o dust_live/drag_kicks.o dust_live/dust_density.o dust_live/dust_init.o dust_live/dust_neighbors.o dust_live/dust_production.o dust_live/dust_timestep.o dust_live/dust_transfer.o dust_live/dust_transfer_kernels.o dust_live/dust_util.o dust_live/dust_vars.o dust_live/grain_sizes.o
INCL    += dust_live/dust_proto.h dust_live/dust_vars.h
SUBDIRS += dust_live
endif

ifeq (BAROTROPIC,$(findstring BAROTROPIC,$(CONFIGVARS)))
OBJS    += barotropic/barotropic.o
INCL    +=
SUBDIRS += barotropic
endif

ifeq (PERTURB_VELOCITIES,$(findstring PERTURB_VELOCITIES,$(CONFIGVARS)))
OBJS    += perturb_velocities.o
endif

ifeq (OPACITIES,$(findstring OPACITIES,$(CONFIGVARS)))
OBJS    += opacities/opacities_combined.o opacities/xztrin21.o
INCL    += opacities/opacities_combined.h
SUBDIRS += opacities
endif

################################
#determine the needed libraries#
################################

# we only need fftw if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIB = $(FFTW_LIBS) -lfftw3 
else
FFTW_LIB = $(FFTW_LIBS) -lfftw3f 
endif
else

# or if POWERSPEC_GRID is activated
ifeq (POWERSPEC_GRID, $(findstring POWERSPEC_GRID, $(CONFIGVARS)))
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
FFTW_LIB = $(FFTW_LIBS) -lfftw3 
else
FFTW_LIB = $(FFTW_LIBS) -lfftw3f 
endif
else
FFTW_LIB = 
endif

endif

ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5LIB  = 
endif

ifeq (NUM_THREADS,$(findstring NUM_THREADS,$(CONFIGVARS)))
THREAD_LIB = 
endif

ifneq (CUDA,$(findstring CUDA,$(CONFIGVARS)))
CUDA_INCL =
CUDA_LIBS = 
endif

ifneq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
HWLOC_INCL = 
HWLOC_LIB = 
endif


ifneq (USE_NBC_FOR_IBARRIER,$(findstring USE_NBC_FOR_IBARRIER,$(CONFIGVARS)))
NBC_INCL = 
NBC_LIB = 
endif

ifneq (FLD_HYPRE,$(findstring FLD_HYPRE,$(CONFIGVARS)))
ifneq (COSMIC_RAYS_DIFFUSION,$(findstring COSMIC_RAYS_DIFFUSION,$(CONFIGVARS)))
ifneq (COSMIC_RAYS_STREAMING,$(findstring COSMIC_RAYS_STREAMING,$(CONFIGVARS)))
ifneq (IMPLICIT_TI,$(findstring IMPLICIT_TI,$(CONFIGVARS)))
ifneq (IMPLICIT_OHMIC_DIFFUSION,$(findstring IMPLICIT_OHMIC_DIFFUSION,$(CONFIGVARS)))
HYPRE_INCL = 
HYPRE_LIB = 
endif
endif
endif
endif
endif

ifneq (VTUNE_INSTRUMENT,$(findstring VTUNE_INSTRUMENT,$(CONFIGVARS)))
VTUNE_INCL =
VTUNE_LIBS =
endif

ifneq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
GRACKLE_INCL =
GRACKLE_LIBS =
endif

ifeq (GALPOT,$(findstring GALPOT,$(CONFIGVARS)))
OBJS   += galpot/galpot.o galpot/potential.o
INCL   += galpot/galpot.h galpot/common.h galpot/potential.h
SUBDIRS += galpot
EIGEN_INCL = 
CPPC = mpic++ -std=c++11
LINKER   = $(CPPC) -lgfortran
endif


##########################
#combine compiler options#
##########################

CFLAGS = $(OPTIMIZE) $(OPT) $(HDF5INCL) $(GSL_INCL) $(EIGEN_INCL) $(FFTW_INCL) $(CVODE_INCL) $(CFITSIO_INCL) $(HEALPIX_INCL) $(GMP_INCL) $(MKL_INCL) $(CUDA_INCL) $(HWLOC_INCL) -I$(BUILD_DIR)  $(NBC_INCL) $(HYPRE_INCL) $(VTUNE_INCL) $(GRACKLE_INCL)

CFLAGS_CUDA = $(CUDA_OPTIMIZE) $(OPT) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(CVODE_INCL) $(CFITSIO_INCL) $(HEALPIX_INCL) $(GMP_INCL) $(MKL_INCL) $(CUDA_INCL)  $(VTUNE_INCL) $(GRACKLE_INCL) $(EIGEN_INCL) -I$(BUILD_DIR)

LIBS = $(MATHLIB) $(HDF5LIB) $(MPICHLIB) $(GSL_LIBS) $(GSLLIB) $(FFTW_LIB) $(GMP_LIBS) $(GMPLIB) $(CVODE_LIB) $(CFITSIO_LIB) $(HEALPIX_LIB) $(MKL_LIBS) $(THREAD_LIB) $(CUDA_LIBS) $(HWLOC_LIB) $(NBC_LIB) $(HYPRE_LIB) $(VTUNE_LIBS) $(GRACKLE_LIBS) $(LAPACK_LIB)

FOPTIONS = $(OPTIMIZE)
FFLAGS = $(FOPTIONS)



SUBDIRS := $(addprefix $(BUILD_DIR)/,$(SUBDIRS))
OBJS := $(addprefix $(BUILD_DIR)/,$(OBJS)) $(BUILD_DIR)/compile_time_info.o $(BUILD_DIR)/compile_time_info_hdf5.o $(BUILD_DIR)/version.o
INCL := $(addprefix $(SRC_DIR)/,$(INCL)) $(BUILD_DIR)/arepoconfig.h

TO_CHECK := $(addsuffix .check, $(OBJS) $(patsubst $(SRC_DIR)%, $(BUILD_DIR)%, $(INCL)) )
TO_CHECK +=  $(BUILD_DIR)/Makefile.check
CONFIG_CHECK = $(BUILD_DIR)/$(notdir $(CONFIG)).check

DOCS_CHECK = $(BUILD_DIR)/README.check

################
#create subdirs#
################
RESULT := $(shell mkdir -p $(SUBDIRS)  )



#############
#build rules#
#############

all: check build 

build: $(EXEC)

$(EXEC): $(OBJS)
	$(LINKER) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)

lib$(LIBRARY).a: $(filter-out $(BUILD_DIR)/main.o,$(OBJS))
	$(AR) -rcs lib$(LIBRARY).a $(OBJS)

clean:
	@echo Cleaning all build files...
	@rm -f $(OBJS) $(EXEC) lib$(LIBRARY).a
	@rm -f $(BUILD_DIR)/compile_time_info.c $(BUILD_DIR)/compile_time_info_hdf5.c $(BUILD_DIR)/arepoconfig.h
	@rm -f $(BUILD_DIR)/version.c
	@rm -f $(TO_CHECK) $(CONFIG_CHECK)
	@rm -rf doxygen/

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.F90 $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 $(INCL) $(MAKEFILES)
	$(FC) -c $(FFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCL) $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cc $(INCL) $(MAKEFILES)
	$(CPPC) $(CFLAGS)-c $< -o $@

$(BUILD_DIR)/compile_time_info.o: $(BUILD_DIR)/compile_time_info.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/compile_time_info_hdf5.o: $(BUILD_DIR)/compile_time_info_hdf5.c $(MAKEFILES)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu $(INCL) $(MAKEFILES)
	$(NVCC)  $(CFLAGS_CUDA)  -c $< -o $@

# sanity checks:

check: $(CONFIG_CHECK)

check_docs: $(DOCS_CHECK)

$(CONFIG_CHECK): $(TO_CHECK) $(CONFIG) check.py 
	@$(PYTHON) check.py 2 $(CONFIG) $(CONFIG_CHECK) defines_extra $(TO_CHECK)

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.f90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.F90
	touch $@

$(BUILD_DIR)/%.o.check: $(SRC_DIR)/%.cc
	touch $@

$(BUILD_DIR)/%.h.check: $(SRC_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.o.check: $(BUILD_DIR)/%.c Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/%.h.check: $(BUILD_DIR)/%.h Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 1 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Makefile.check: Makefile Template-Config.sh defines_extra check.py
	@$(PYTHON) check.py 3 $< $@ Template-Config.sh defines_extra

$(BUILD_DIR)/Config.check: Template-Config.sh check.py
	@$(PYTHON) check.py 4 Template-Config.sh $@

