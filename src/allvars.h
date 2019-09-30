/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/allvars.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <mpi.h>

#include <stddef.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

#include "arepoconfig.h"

#ifdef IMPOSE_PINNING
#include <hwloc.h>
#endif

#include "tags.h"
#include "dtypes.h"
#include "timestep.h"

#define  AREPO_VERSION   "Arepo 0.1"    /* code version string */
#define  TENET_VERSION   "0.1"    /* code version string */

/* default values for unspecified config options */

#ifndef LOAD_TYPES
#define LOAD_TYPES 0xff
#endif

#if !defined(AMR) && !defined(VORONOI)
#define GRADIENTS_GREEN_GAUSS
#endif

#ifndef AMR
#define VORONOI
#endif

#if defined (REFINEMENT_SPLIT_CELLS) || defined (REFINEMENT_MERGE_CELLS)
#define REFINEMENT
#else
#undef REFINEMENT
#endif

#ifdef GFM_STELLAR_EVOLUTION
#if !defined(GFM_CONST_IMF) && !defined(GFM_VARIABLE_IMF)
#define GFM_CONST_IMF 0
#endif
#endif

#ifndef NTYPES
#define NTYPES 6
#endif

#ifndef NSOFTTYPES
#define NSOFTTYPES NTYPES
#endif

#ifdef ADAPTIVE_HYDRO_SOFTENING
#ifndef NSOFTTYPES_HYDRO
#define NSOFTTYPES_HYDRO 64
#endif
#else
#undef NSOFTTYPES_HYDRO
#define NSOFTTYPES_HYDRO 0
#endif

#if defined(SAVE_HSML_IN_SNAPSHOT) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL) || defined(DMPIC)
#define SUBFIND_CALC_MORE
#endif

/* restrictions on config option combinations */

#if (NSOFTTYPES + NSOFTTYPES_HYDRO) >= 254
#error "(NSOFTTYPES + NSOFTTYPES_HYDRO) >= 254"
#endif

#if NSOFTTYPES < 2
#error "NSOFTTYPES < 2"
#endif

#if defined(HOST_MEMORY_REPORTING) && !defined(__linux__)
#error "HOST_MEMORY_REPORTING only works under Linux."
#endif

#if defined(USE_DIRECT_IO_FOR_RESTARTS) && !defined(__linux__)
#error "USE_DIRECT_IO_FOR_RESTARTS only works under Linux."
#endif


#ifdef INDIVIDUAL_GRAVITY_SOFTENING
#if !((INDIVIDUAL_GRAVITY_SOFTENING+0) >= 1)
#error "set INDIVIDUAL_GRAVITY_SOFTENING to a bitmask of particle types"
#endif
#endif

#if (defined(REPOSITION_ON_POTMIN) || defined(BH_FRICTION)) && !defined(MEASURE_POTMIN_AROUND_BH)
#define MEASURE_POTMIN_AROUND_BH
#endif

#ifdef OUTPUTPOTENTIAL
#ifndef EVALPOTENTIAL
#error "the option OUTPUTPOTENTIAL requires EVALPOTENTIAL"
#endif
#endif

#ifdef MEASURE_POTMIN_AROUND_BH
#ifndef EVALPOTENTIAL
#error "the option REPOSITION_ON_POTMIN requires EVALPOTENTIAL"
#endif
#endif

#ifdef COMPUTE_POTENTIAL_ENERGY
#ifndef EVALPOTENTIAL
#error "the option COMPUTE_POTENTIAL_ENERGY requires EVALPOTENTIAL"
#endif
#endif

#if defined(CELL_CENTER_GRAVITY) && defined(SELFGRAVITY)
#ifndef HIERARCHICAL_GRAVITY
#error "use of option CELL_CENTER_GRAVITY requires HIERARCHICAL_GRAVITY"
#endif
#endif

#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
#if !defined(EOS_DEGENERATE) || !defined(NUCLEAR_NETWORK)
#error "use of option TRACER_TRAJECTORY_EXTENDED_OUTPUT requires EOS_DEGENERATE and NUCLEAR_NETWORK"
#endif
#endif

#ifdef NUCLEAR_NETWORK
#ifndef OUTPUT_DIVVEL
#error "use of option NUCLEAR_NETWORK requires OUTPUT_DIVVEL"
#endif
#endif

#ifdef COSMIC_RAYS_DIFFUSION
#ifndef TETRA_INDEX_IN_FACE
#error "use of option COSMIC_RAYS_DIFFUSION requires TETRA_INDEX_IN_FACE"
#endif
#endif

#ifdef BECDM
#ifndef EVALPOTENTIAL
#error "use of option BECDM requires EVALPOTENTIAL"
#endif
#ifndef PMGRID
#error "use of option BECDM requires PMGRID"
#endif
#endif

#ifdef MONOTONE_CONDUCTION
#ifndef TETRA_INDEX_IN_FACE
#error "use of option MONOTONE_CONDUCTION requires TETRA_INDEX_IN_FACE"
#endif
#endif

#ifdef SIDM
#if defined(HIERARCHICAL_GRAVITY) || defined(ALLOW_DIRECT_SUMMATION)
#error "use of option SIDM does not work with HIERARCHICAL_GRAVITY or ALLOW_DIRECT_SUMMATION since we need a full gravity tree for neighbor searches at each time"
#endif
#endif

#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
#ifndef FORCE_EQUAL_TIMESTEPS
#error "NUCLEAR_NETWORK_TIMESTEP_LIMITER only works with FORCE_EQUAL_TIMESTEPS at the moment"
#endif
#endif

#ifdef NUCLEAR_NETWORK_USE_SHOCKFINDER
#ifndef SHOCK_FINDER_ON_THE_FLY
#error "NUCLEAR_NETWORK_USE_SHOCKFINDER requires SHOCK_FINDER_ON_THE_FLY"
#endif
#endif

#ifdef MRT

#ifdef MRT_NO_UV
#define UV_BINS 0
#else
#ifdef MRT_MULTI_FREQUENCY
#define UV_BINS 1
#else
#define UV_BINS 1
#endif
#endif

#ifdef MRT_IR
#define IR_BINS 1
#else
#define IR_BINS 0
#endif

#define MRT_BINS (UV_BINS+IR_BINS)

#ifdef MRT_IR_LTE
#define RADIATION_CONSTANT 7.5657e-15  /* erg cm^-3 K^-4 */
#endif

#ifdef MRT_BH
extern struct bh_particle
{
  int index;
  MyDouble NumNgb;
  MyDouble NormSph;
  MyDouble Dhsmlrho;
  MyFloat TotalPhotReleased[MRT_BINS];
} *BHParticle;
#endif

#endif

/* optional additional headers based on config options */

#include "timer.h"

#if defined(COOLING) && !defined(GRACKLE)
#include "cooling/cooling_vars.h"
#endif

#ifdef ATOMIC_DM
#include "atomic_dm/cooling_atomic_dm_vars.h"
#endif

#ifdef NUCLEAR_NETWORK
#include "network/network.h"
#include "network/integrate.h"
#endif

#ifdef NETWORK_NSE
#include "network/network_nse.h"
#endif

#ifdef GFM_STELLAR_EVOLUTION
#include "GFM/stellar_evolution_vars.h"
#endif

#ifdef GFM_COOLING_METAL
#include "GFM/cooling_metal_vars.h"
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
#include "GFM/stellar_photometrics_vars.h"
#endif

#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS_VARIABLE)
#include "GFM/winds_vars.h"
#endif

#ifdef DUST_LIVE
#include "dust_live/dust_vars.h"
#endif

#ifdef DVR_RENDER
#include "dvr_render/dvr_render.h"
#endif

#ifdef SIDM
#include "sidm/sidm_vars.h"
#endif

#ifdef SGCHEM
#include "SGChem/sgchem_def.h"
#endif

#ifdef SNE_FEEDBACK
#include "sne/sne.h"
#endif

#ifdef FLD
#include "fld/fld.h"
#endif

#ifdef COSMIC_RAYS
#include "cosmic_rays/cosmic_rays.h"
#endif

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
#include "shock_finder/shock_finder_fields.h"
#endif

#ifdef ADDBACKGROUNDGRID
#include "add_backgroundgrid/add_bggrid.h"
#endif

/* function mappings and macros */

#ifdef MPI_HYPERCUBE_ALLGATHERV
#define MPI_Allgatherv MPI_hypercube_Allgatherv
#endif

#ifdef MPISENDRECV_CHECKSUM
#define MPI_Sendrecv MPI_Check_Sendrecv
#endif

#ifdef MPISENDRECV_SIZELIMIT
#define MPI_Sendrecv MPI_Sizelimited_Sendrecv
#endif

#define  terminate(...) {if(FlagNyt==0){char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "TERMINATE: ******!!!!!******  Code termination on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s: %s\n", termbuf1, termbuf2); fflush(stdout); FlagNyt=1; MPI_Abort(MPI_COMM_WORLD, 1);} exit(0);}

#define  warn(...) {char termbuf1[1000], termbuf2[1000]; sprintf(termbuf1, "WARNING: Code warning on task=%d, function %s(), file %s, line %d", ThisTask, __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s: %s\n", termbuf1, termbuf2); myflush(stdout); FILE *fd=fopen("WARNINGS", "a"); fprintf(fd, "%s: %s\n", termbuf1, termbuf2); fclose(fd);}

/* define an "assert" macro which outputs MPI task (we do NOT want to
   call MPI_Abort, because then the assertion failure isn't caught in
   the debugger) */
#ifndef NDEBUG
#define myassert(cond)							\
  if(!(cond)) {								\
    char termbuf[1000];							\
    sprintf(termbuf, "Assertion failure!\n\ttask=%d, function %s(), file %s, line %d:\n\t%s\n", ThisTask, __FUNCTION__, __FILE__, __LINE__, #cond); \
    printf("%s", termbuf); myflush(stdout);				\
    assert(0);								\
  }
#else
#define myassert(cond)
#endif

#ifndef DISABLE_MEMORY_MANAGER
#define  mymalloc(x, y)              mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 0, NULL)
#define  mymalloc_g(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 0, callorigin)
#define  mymalloc_clear(x, y)        mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__, 1, NULL)
#define  mymalloc_movable(x, y, z)   mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, NULL)
#define  mymalloc_movable_g(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__, callorigin)

#define  myrealloc(x, y)             myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define  myrealloc_movable(x, y)     myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)

#define  myfree(x)                   myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define  myfree_movable(x)           myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)

#else
#define  mymalloc(x, y)              malloc(y)
#define  mymalloc_g(x, y)            malloc(y)
#define  mymalloc_clear(x, y)        calloc(1, y)
#define  mymalloc_movable(x, y, z)   malloc(z)
#define  mymalloc_movable_g(x, y, z) malloc(z)

#define  myrealloc(x, y)             realloc(x, y)
#define  myrealloc_movable(x, y)     realloc(x, y)

#define  myfree(x)                   free(x)
#define  myfree_movable(x)           free(x)

#endif

#define  MAX_FIRST_ELEMENTS_CONSIDERED     5    /* This sets the number of lowest loaded tasks to be considered for assignment of next domain patch */

#define  NUMBER_OF_MEASUREMENTS_TO_RECORD  6

#ifndef  GRAVCOSTLEVELS
#ifdef TGSET
#define  GRAVCOSTLEVELS      30
#else
#define  GRAVCOSTLEVELS      6
#endif
#endif

#define MODE_LOCAL_NO_EXPORT   -1
#define MODE_LOCAL_PARTICLES    0
#define MODE_IMPORTED_PARTICLES 1
#define MODE_FINISHED           2


#ifndef DIRECT_SUMMATION_THRESHOLD
#define DIRECT_SUMMATION_THRESHOLD 3000
#endif


#define MODE_FIRST_HALFSTEP   0
#define MODE_SECOND_HALFSTEP  1


#define FLAG_PARTIAL_TREE     0
#define FLAG_FULL_TREE        1


#ifndef MPI_MESSAGE_SIZELIMIT_IN_MB
#define MPI_MESSAGE_SIZELIMIT_IN_MB 200
#endif

#define MPI_MESSAGE_SIZELIMIT_IN_BYTES ((MPI_MESSAGE_SIZELIMIT_IN_MB) * 1024LL * 1024LL)

#define COMMBUFFERSIZE (32 * 1024LL * 1024LL)

#ifdef BLACK_HOLES
#define BH_QUASAR_MODE 0
#define BH_RADIO_MODE 1
#endif


#if !defined(NUM_THREADS)
#define NUM_THREADS 1
#endif

#if (NUM_THREADS > 1)
#include <omp.h>
#endif

#if (NUM_THREADS > 1)
extern omp_lock_t *ParticleLocks;
extern omp_lock_t *Ngb_NodeLocks;
#endif

extern int Nforces;
extern int *TargetList;

extern struct thread_data
{
  int Nexport              __attribute__ ((__aligned__(64)));  /* to align on different cache lines */
  int NexportNodes;
  int Interactions;
  int dummy;
  double Cost;

  double Costtotal;             /*!< The total cost of the particles/nodes processed by each thread */
  double Ewaldcount;            /*!< The total cost for the Ewald correction per thread */
  int FirstExec;                /*!< Keeps track, if a given thread executes the gravity_primary_loop() for the first time */

  size_t ExportSpace;
  size_t InitialSpace;
  size_t ItemSize;

  int *P_CostCount;
  int *TreePoints_CostCount;
  int *Node_CostCount;

  struct data_partlist *PartList;

  int *Ngblist;
  double *R2list;
  int *Exportflag;
  int *toGoDM;
  int *toGoSph;

#ifdef GENERIC_ASYNC
  struct data_partlist *PartListOld;
  int NexportOld;
#endif

} Thread[NUM_THREADS];


#ifdef SIDM
#include "sidm/sidm_vars.h"
#endif

/* If we use a static Voronoi mesh with local timestepping and no rebuild of
 * the static mesh, then we need to backup the face areas before calling
 * compute_interface_fluxes(), because this function calls face_get_normals()
 * which sets some face area to 0 under some circumstances */
#if defined(VORONOI_STATIC_MESH) && !defined(FORCE_EQUAL_TIMESTEPS) && !defined(VORONOI_STATIC_MESH_DO_DOMAIN_DECOMPOSITION)
#define VORONOI_BACKUP_RESTORE_FACE_AREAS
#else
#undef VORONOI_BACKUP_RESTORE_FACE_AREAS
#endif

#ifdef IMPOSE_PINNING
extern hwloc_cpuset_t cpuset_thread[NUM_THREADS];
#endif

#ifdef ONEDIMS
#define ALLOC_TOLERANCE      0.3
#else
#define ALLOC_TOLERANCE      0.1
#endif
#define ALLOC_STARBH_ROOM    0.02



#ifdef TOLERATE_WRITE_ERROR
#define IO_TRIALS 20
#define IO_SLEEP_TIME 10
#endif

/* calculate appropriate value of MAXSCALARS */

#if defined(EOS_DEGENERATE) || defined(REFINEMENT_HIGH_RES_GAS) || defined(METALS) || defined(GFM_STELLAR_EVOLUTION) || \
  defined(PASSIVE_SCALARS) || defined(DELAYED_COOLING_TURB) || defined(REFINEMENT_RPS) || defined(SGCHEM) || defined(EOS_OPAL) || defined (COSMIC_RAYS)

#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
#define COUNT_EOS EOS_NSPECIES
#else
#define COUNT_EOS 0
#endif

#ifdef MHD_THERMAL_ENERGY_SWITCH
#define COUNT_MHD_TES 1
#else
#define COUNT_MHD_TES 0
#endif

#ifdef  REFINEMENT_HIGH_RES_GAS
#define COUNT_REFINE 1
#else
#define COUNT_REFINE 0
#endif

#ifdef  METALS
#define COUNT_METALS 1
#else
#define COUNT_METALS 0
#endif

#ifdef  GFM_STELLAR_EVOLUTION
#define COUNT_STELLAR_EVOLUTION (GFM_N_CHEM_ELEMENTS+1+GFM_DUST_COUNT_SCALARS)
#else
#define COUNT_STELLAR_EVOLUTION 0
#endif

#ifdef SGCHEM
#ifdef SGCHEM_VARIABLE_Z
#define COUNT_SGCHEM SGCHEM_NUM_ADVECTED_SPECIES+SGCHEM_NUM_ELEMS
#else
#define COUNT_SGCHEM SGCHEM_NUM_ADVECTED_SPECIES
#endif
#else
#define COUNT_SGCHEM 0
#endif

#ifdef PASSIVE_SCALARS
#define COUNT_PASSIVE_SCALARS PASSIVE_SCALARS
#else
#define COUNT_PASSIVE_SCALARS 0
#endif

#ifdef DELAYED_COOLING_TURB
#define COUNT_DELAYED_COOLING 1
#else
#define COUNT_DELAYED_COOLING 0
#endif

#ifdef REFINEMENT_RPS
#define COUNT_RPS 1
#else
#define COUNT_RPS 0
#endif

#ifdef COSMIC_RAYS
#define COUNT_CR 1
#else
#define COUNT_CR 0
#endif

#define MAXSCALARS (COUNT_EOS + COUNT_MHD_TES + COUNT_REFINE + COUNT_METALS + COUNT_STELLAR_EVOLUTION + COUNT_PASSIVE_SCALARS + COUNT_DELAYED_COOLING + COUNT_RPS + COUNT_SGCHEM + COUNT_CR)

#endif

#ifdef RUNGE_KUTTA_FULL_UPDATE
#include "runge_kutta_full.h"
#endif

#if defined(COOLING) && defined(CIRCUMSTELLAR_IRRADIATION)
#define SIMPLE_COOLING
#endif

#ifdef RT_HEALPIX_NSIDE
#ifdef TWODIMS
/* here we keep RT_N_DIR in  case RT_HEALPIX_NSIDE is enabled */
#else
#undef RT_N_DIR
#define RT_N_DIR  (12*RT_HEALPIX_NSIDE*RT_HEALPIX_NSIDE)
#endif
#else

/* in this case, RT_N_DIR should be set in the Makefile */
#endif

/* We define RT_N_DIR=0 if no RT_ADVECT because it simplifies slice code */
#ifndef RT_ADVECT
#define RT_N_DIR 0
#endif

/* calculate appropriate value of MAXGRADIENTS */

#define COUNT_GRAD_DEFAULT 5

#ifdef RT_ADVECT
#define COUNT_GRAD_RT RT_N_DIR
#else
#define COUNT_GRAD_RT 0
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
#define COUNT_GRAD_ENTR 1
#else
#define COUNT_GRAD_ENTR 0
#endif

#ifdef TGCHEM
#define COUNT_GRAD_TGCHEM 1
#else
#define COUNT_GRAD_TGCHEM 0
#endif

#ifdef VARIABLE_GAMMA
#define COUNT_GRAD_GAMMA 2
#else
#define COUNT_GRAD_GAMMA 0
#endif

#ifdef MHD
#define COUNT_GRAD_MHD 3
#else
#define COUNT_GRAD_MHD 0
#endif

#ifdef MHD_CT
#define COUNT_GRAD_MHD_CT 3
#else
#define COUNT_GRAD_MHD_CT 0
#endif

#ifdef MAXSCALARS
#define COUNT_GRAD_SCALARS MAXSCALARS
#else
#define COUNT_GRAD_SCALARS 0
#endif

#ifdef TRACER_FIELD
#define COUNT_GRAD_TRACER 1
#else
#define COUNT_GRAD_TRACER 0
#endif

#if defined(GLOBAL_VISCOSITY) || defined(USE_KINEMATIC_VISCOSITY) || defined(ALPHA_VISCOSITY) || defined(LOCAL_VISCOSITY)
#define VISCOSITY
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_PLANET_GROWTH)
#define BLACK_HOLES
#define BH_BONDI_DISK_VORTICITY
#define DRAINGAS 3
#endif

#if defined(SPECIAL_RELATIVITY) || defined(DEREFINE_GENTLY) || defined(GENERAL_RELATIVITY) || defined(VORONOI_PROJ_TAU)
#define COUNT_GRAD_UTHERM 1
#else
#define COUNT_GRAD_UTHERM 0
#endif

#if defined(CONDUCTION_SATURATION) || defined(NON_LINEAR_SLOPE_LIMITERS) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
#define COUNT_GRAD_CONDUCTION 1
#else
#define COUNT_GRAD_CONDUCTION 0
#endif

#if defined(MRT) && defined(MRT_LSF_GRADIENTS)

#ifndef MRT_TIME_EXTRAPOLATION
#define COUNT_GRAD_MRT (5*MRT_BINS)
#else
#define COUNT_GRAD_MRT (14*MRT_BINS)
#endif

#else
#define COUNT_GRAD_MRT 0
#endif

#ifdef DVR_RENDER
#define COUNT_GRAD_DVR_RENDER DVR_NUM_FIELDS
#else
#define COUNT_GRAD_DVR_RENDER 0
#endif

#ifdef COSMIC_RAYS
#define COUNT_CR 1
#else
#define COUNT_CR 0
#endif

#ifdef FLD
#define COUNT_FLD 1
#else
#define COUNT_FLD 0
#endif

#define MAXGRADIENTS (COUNT_GRAD_DEFAULT + COUNT_GRAD_ENTR + COUNT_GRAD_UTHERM + COUNT_GRAD_TGCHEM + COUNT_GRAD_GAMMA + COUNT_GRAD_MHD + COUNT_GRAD_MHD_CT + COUNT_GRAD_SCALARS + COUNT_GRAD_TRACER + COUNT_GRAD_DVR_RENDER + COUNT_GRAD_CONDUCTION + COUNT_CR + COUNT_FLD)
#define RT_MAXGRADIENTS COUNT_GRAD_RT

#define MAXRTGRADIENTS COUNT_GRAD_MRT

#if defined(OTVET) 
#if !defined(OTVET_MULTI_FREQUENCY)
#define OT_N_BINS 1
#else
#define OT_N_BINS 4
#endif
#define ASSIGN_ADD(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#endif


/*************************************/


/** For Peano-Hilbert order.
 *  Note: Maximum is 10 to fit in 32-bit integer,
 *  maximum is 21 to fit into 64-bit integer,
 *  and 42 is the absolute maximum, for which 128-bit integers are needed
 */
#ifndef  BITS_PER_DIMENSION
#define  BITS_PER_DIMENSION 42
#endif
#if (BITS_PER_DIMENSION <= 21)
typedef unsigned long long peanokey;
#else
typedef __int128 peanokey;
#endif
#if (BITS_PER_DIMENSION <= 31)
typedef unsigned int peano1D;
#else
#if (BITS_PER_DIMENSION <= 42)
typedef unsigned long long peano1D;
#else
#error "BITS_PER_DIMENSION can be at most 42"
#endif
#endif

#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))


#define  MAX_FLOAT_NUMBER 1e37
#define  MIN_FLOAT_NUMBER 1e-37
#define  MAX_DOUBLE_NUMBER 1e306
#define  MIN_DOUBLE_NUMBER 1e-306

#define  BHPOTVALUEINIT 1.0e30

#ifdef DOUBLEPRECISION
#if (DOUBLEPRECISION==2)
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#else
#define  MAX_REAL_NUMBER  MAX_DOUBLE_NUMBER
#define  MIN_REAL_NUMBER  MIN_DOUBLE_NUMBER
#endif
#else
#define  MAX_REAL_NUMBER  MAX_FLOAT_NUMBER
#define  MIN_REAL_NUMBER  MIN_FLOAT_NUMBER
#endif


#ifndef  GAMMA
#define  GAMMA         (5. / 3.)  /**< adiabatic index of simulated gas */
#endif

#define  GAMMA_MINUS1  (GAMMA - 1.)
#define  GAMMA_PLUS1   (GAMMA + 1.)

#ifndef  HYDROGEN_ONLY
#define  HYDROGEN_MASSFRAC 0.76 /**< mass fraction of hydrogen, relevant only for radiative cooling */
#define  HE_ABUND ((1. / HYDROGEN_MASSFRAC - 1.) / 4.)
#else
#define  HYDROGEN_MASSFRAC 1.    /**< mass fraction of hydrogen, relevant only for radiative cooling */
#define  HE_ABUND 0.
#endif


#define  METAL_YIELD       0.02 /**< effective metal yield for star formation */

/* ... often used physical constants (cgs units; NIST 2010) */

#define  GRAVITY         6.6738e-8
#define  SOLAR_MASS      1.989e33
#define  SOLAR_LUM       3.826e33
#define  SOLAR_EFF_TEMP  5.780e3
#define  RAD_CONST       7.5657e-15
#define  AVOGADRO        6.02214e23
#define  BOLTZMANN       1.38065e-16
#define  GAS_CONST       8.31446e7
#if defined(RT_SLOWLIGHT) || defined(MRT_SLOWLIGHT)
#define  CLIGHT          2.99792458e9
#else
#define  CLIGHT          2.99792458e10
#endif

#ifdef MRT
#define MINDENSPHOT      1e-25
#endif

#define  PLANCK          6.6260695e-27
#define  PARSEC          3.085678e18
#define  KILOPARSEC      3.085678e21
#define  MEGAPARSEC      3.085678e24
#define  ASTRONOMICAL_UNIT   1.49598e13
#define  PROTONMASS      1.67262178e-24
#define  ELECTRONMASS    9.1093829e-28
#define  THOMPSON        6.65245873e-25
#define  ELECTRONCHARGE  4.8032042e-10
#define  HUBBLE          3.2407789e-18  /* in h/sec */
#define  LYMAN_ALPHA     1215.6e-8      /* 1215.6 Angstroem */
#define  LYMAN_ALPHA_HeII  303.8e-8     /* 303.8 Angstroem */
#define  OSCILLATOR_STRENGTH       0.41615
#define  OSCILLATOR_STRENGTH_HeII  0.41615
#define  ELECTRONVOLT_IN_ERGS      1.60217656e-12

#define  SEC_PER_GIGAYEAR   3.15576e16
#define  SEC_PER_MEGAYEAR   3.15576e13
#define  SEC_PER_YEAR       3.15576e7

#ifndef FOF_PRIMARY_LINK_TYPES
#define FOF_PRIMARY_LINK_TYPES 2
#endif

#ifndef FOF_SECONDARY_LINK_TYPES
#define FOF_SECONDARY_LINK_TYPES 0
#endif


#ifndef ASMTH
/** ASMTH gives the scale of the short-range/long-range force split in units of FFT-mesh cells */
#define ASMTH 1.25
#endif
#ifndef RCUT
/** RCUT gives the maximum distance (in units of the scale used for the force split) out to which short-range
 * forces are evaluated in the short-range tree walk.
 */
#define RCUT  4.5
#endif

#ifdef DECOUPLE_TIMESTEPS
#define MAX_TIMEBIN_DIFFERENCE 3
#else
#define MAX_TIMEBIN_DIFFERENCE 0
#endif

#ifdef GFM_AGN_RADIATION
#define GFM_MAX_TIMEBINS_WITHOUT_AGN_RAD 5
#endif

#ifdef TRACER_PARTICLE
#ifndef TRACER_PART_STORE_WHAT
#define TRACER_PART_STORE_WHAT 0
#endif
#define TRACER_PART_TMAX ((TRACER_PART_STORE_WHAT) & 1)
#define TRACER_PART_TMAX_TIME ((TRACER_PART_STORE_WHAT) & 2)
#define TRACER_PART_TMAX_RHO ((TRACER_PART_STORE_WHAT) & 4)
#define TRACER_PART_RHOMAX ((TRACER_PART_STORE_WHAT) & 8)
#define TRACER_PART_RHOMAX_TIME ((TRACER_PART_STORE_WHAT) & 16)
#define TRACER_PART_MACHMAX ((TRACER_PART_STORE_WHAT) & 32)
#define TRACER_PART_ENTMAX ((TRACER_PART_STORE_WHAT) & 64)
#define TRACER_PART_ENTMAX_TIME ((TRACER_PART_STORE_WHAT) & 128)
#endif

#ifdef TRACER_MC
#ifndef TRACER_MC_STORE_WHAT
#define TRACER_MC_STORE_WHAT 0
#endif
#define TRACER_MC_TMAX ((TRACER_MC_STORE_WHAT) & 1)
#define TRACER_MC_TMAX_TIME ((TRACER_MC_STORE_WHAT) & 2)
#define TRACER_MC_TMAX_RHO ((TRACER_MC_STORE_WHAT) & 4)
#define TRACER_MC_RHOMAX ((TRACER_MC_STORE_WHAT) & 8)
#define TRACER_MC_RHOMAX_TIME ((TRACER_MC_STORE_WHAT) & 16)
#define TRACER_MC_MACHMAX ((TRACER_MC_STORE_WHAT) & 32)
#define TRACER_MC_ENTMAX ((TRACER_MC_STORE_WHAT) & 64)
#define TRACER_MC_ENTMAX_TIME ((TRACER_MC_STORE_WHAT) & 128)
#define TRACER_MC_LAST_STAR_TIME ((TRACER_MC_STORE_WHAT) & 256)
#define TRACER_MC_WIND_COUNTER ((TRACER_MC_STORE_WHAT) & 512)
#define TRACER_MC_EXCHANGE_COUNTER ((TRACER_MC_STORE_WHAT) & 1024)
#define TRACER_MC_EXCHANGE_DISTANCE ((TRACER_MC_STORE_WHAT) & 2048)
#define TRACER_MC_EXCHANGE_DISTANCE_ERROR ((TRACER_MC_STORE_WHAT) & 4096)
#define TRACER_MC_SHOCKMACHNUM_MAX ((TRACER_MC_STORE_WHAT) & 8192)
#endif

#if TRACER_MC_SHOCKMACHNUM_MAX && !defined(SHOCK_FINDER_ON_THE_FLY)
#error "Enabling TRACER_MC_SHOCKMACHNUM_MAX (8192) requires SHOCK_FINDER_ON_THE_FLY."
#endif

#define MAXLEN_OUTPUTLIST 1100  /**< maxmimum number of entries in output list */

#define MAXLEN_PATH       256   /**< maximum length of various filenames (full path) */

#define MAXLEN_PARAM_TAG  50    /**< maximum length of the tag of a parameter in the parameter file */
#define MAXLEN_PARAM_VALUE 200  /**< maximum length of the value of a parameter in the parameter file */
#define MAX_PARAMETERS 300      /**< maximum number of parameters in the parameter file */

#define DRIFT_TABLE_LENGTH  1000        /**< length of the lookup table used to hold the drift and kick factors */

#define BASENUMBER          100

#define HIGHRESMASSFAC 0.5

#ifndef CSND_FRAC_BH_MERGE
#define CSND_FRAC_BH_MERGE 0.5
#endif

#ifndef BH_REPOSITION_POTMIN_TRUST_THRESHOLD
#define BH_REPOSITION_POTMIN_TRUST_THRESHOLD  0.9
#endif

#define MAXITER 300000

#ifndef FOF_LINKLENGTH
#define FOF_LINKLENGTH 0.2
#endif

#ifndef FOF_GROUP_MIN_LEN
#define FOF_GROUP_MIN_LEN 32
#endif

#ifdef RT_ADVECT

#ifndef RT_SFR_SOURCES
#define N_SOURCES 1
double Source_Pos[N_SOURCES][3];
double Source_Lum[N_SOURCES];
int Source_ID[N_SOURCES];
#endif

#endif /* end of RT_ADVECT */

typedef struct
{
  double r;
  double mass;
}
sort_r2list;

typedef struct
{
  MyFloat r2;
  int index;
}
r2type;

#ifdef TGSET
#include "tgset/tgset.h"
#endif

#ifdef TGCHEM
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include "tgchem/tgchem.h"
#endif

#ifdef HEALRAY
#include "healray/healray.h"
#endif

#ifdef SINKS
#include "sinks/sinks.h"
#endif

#ifdef SINK_PARTICLES
#include "sink_particles/sink_particles.h"
#endif

#include "mesh.h"

#ifdef VORONOI
#include "voronoi.h"
#endif

struct unbind_data
{
  int index;
};

#ifdef FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
extern MPI_Status mpistat;
#undef MPI_STATUS_IGNORE
#define MPI_STATUS_IGNORE &mpistat
#endif

#define FLT(x) (x)


#ifndef M_PI
#define M_PI      3.14159265358979323846
#endif

#define TO_MBYTE_FAC  (1.0 / (1024.0 * 1024.0))


#ifdef   ONEDIMS
#define  NUMDIMS 1
#define  KERNEL_COEFF_1  (4.0/3)
#define  KERNEL_COEFF_2  (8.0)
#define  KERNEL_COEFF_3  (24.0)
#define  KERNEL_COEFF_4  (16.0)
#define  KERNEL_COEFF_5  (8.0/3)
#define  KERNEL_COEFF_6  (-8.0)
#define  NORM_COEFF      2.0
#else
#ifndef  TWODIMS
#define  NUMDIMS 3              /**< For 3D-normalized kernel */
#define  KERNEL_COEFF_1  2.546479089470 /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  15.278874536822
#define  KERNEL_COEFF_3  45.836623610466
#define  KERNEL_COEFF_4  30.557749073644
#define  KERNEL_COEFF_5  5.092958178941
#define  KERNEL_COEFF_6  (-15.278874536822)
#define  NORM_COEFF      4.188790204786 /**< Coefficient for kernel normalization. Note:  4.0/3 * PI = 4.188790204786 */
#else
#define  NUMDIMS 2              /**< For 2D-normalized kernel */
#define  KERNEL_COEFF_1  (5.0/7*2.546479089470) /**< Coefficients for SPH spline kernel and its derivative */
#define  KERNEL_COEFF_2  (5.0/7*15.278874536822)
#define  KERNEL_COEFF_3  (5.0/7*45.836623610466)
#define  KERNEL_COEFF_4  (5.0/7*30.557749073644)
#define  KERNEL_COEFF_5  (5.0/7*5.092958178941)
#define  KERNEL_COEFF_6  (5.0/7*(-15.278874536822))
#define  NORM_COEFF      M_PI   /**< Coefficient for kernel normalization. */
#endif
#endif /* ONEDIMS */


#define  SOFTFAC1  10.666666666667 /**< Coefficients for gravitational softening */
#define  SOFTFAC2  32.0
#define  SOFTFAC3  (-38.4)
#define  SOFTFAC4  (-2.8)
#define  SOFTFAC5  5.333333333333
#define  SOFTFAC6  6.4
#define  SOFTFAC7  (-9.6)
#define  SOFTFAC8  21.333333333333
#define  SOFTFAC9  (-48.0)
#define  SOFTFAC10 38.4
#define  SOFTFAC11 (-10.666666666667)
#define  SOFTFAC12 (-0.066666666667)
#define  SOFTFAC13 (-3.2)
#define  SOFTFAC14  0.066666666667
#define  SOFTFAC15  (-16.0)
#define  SOFTFAC16  9.6
#define  SOFTFAC17  (-2.133333333333)
#define  SOFTFAC18  128.0
#define  SOFTFAC19  (-115.2)
#define  SOFTFAC20  21.333333333333
#define  SOFTFAC21  (-96.0)
#define  SOFTFAC22  115.2
#define  SOFTFAC23  (-42.666666666667)
#define  SOFTFAC24  0.1333333333333


#ifdef PERIODIC
extern MyDouble boxSize, boxHalf;
#ifdef LONG_X
extern MyDouble boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
extern MyDouble boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
extern MyDouble boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif

#ifdef AMR
#include "amr/amr.h"
#endif

#ifdef DG
#include "dg/dg_vars.h"
#endif

#ifdef DG_TEST_PROBLEM
#include "dg/dg_test_problems.h"
#endif

#ifdef GRACKLE
#define CONFIG_BFLOAT_8
#include <grackle.h>
#include "grackle/grackle_def.h"
#endif

#if !defined(GRAVITY_NOT_PERIODIC) || defined(GRAVITY_TALLBOX)
#define GRAVITY_NEAREST_X(x) (xtmp=(x),(xtmp>boxHalf_X)?(xtmp-boxSize_X):( (xtmp< -boxHalf_X)?(xtmp+boxSize_X):(xtmp) ) )
#define GRAVITY_NEAREST_Y(x) (ytmp=(x),(ytmp>boxHalf_Y)?(ytmp-boxSize_Y):( (ytmp< -boxHalf_Y)?(ytmp+boxSize_Y):(ytmp) ) )
#define GRAVITY_NEAREST_Z(x) (ztmp=(x),(ztmp>boxHalf_Z)?(ztmp-boxSize_Z):( (ztmp< -boxHalf_Z)?(ztmp+boxSize_Z):(ztmp) ) )
#else
#define GRAVITY_NEAREST_X(x) (x)
#define GRAVITY_NEAREST_Y(x) (x)
#define GRAVITY_NEAREST_Z(x) (x)
#endif

#ifdef GRAVITY_TALLBOX
#if GRAVITY_TALLBOX==0
#undef  GRAVITY_NEAREST_X
#define GRAVITY_NEAREST_X(x) (x)
#endif
#if GRAVITY_TALLBOX==1
#undef  GRAVITY_NEAREST_Y
#define GRAVITY_NEAREST_Y(x) (x)
#endif
#if GRAVITY_TALLBOX==2
#undef  GRAVITY_NEAREST_Z
#define GRAVITY_NEAREST_Z(x) (x)
#endif
#endif


#if !defined(GRAVITY_NOT_PERIODIC)
#define FOF_NEAREST_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define FOF_NEAREST_LONG_Y(x) (ytmp=fabs(x),(ytmp>boxHalf_Y)?(boxSize_Y-ytmp):ytmp)
#define FOF_NEAREST_LONG_Z(x) (ztmp=fabs(x),(ztmp>boxHalf_Z)?(boxSize_Z-ztmp):ztmp)
#else
#define FOF_NEAREST_LONG_X(x) fabs(x)
#define FOF_NEAREST_LONG_Y(x) fabs(x)
#define FOF_NEAREST_LONG_Z(x) fabs(x)
#endif


#ifdef PERIODIC
#define NGB_PERIODIC_LONG_X(x) (xtmp=fabs(x),(xtmp>boxHalf_X)?(boxSize_X-xtmp):xtmp)
#define NGB_PERIODIC_LONG_Y(x) (ytmp=fabs(x),(ytmp>boxHalf_Y)?(boxSize_Y-ytmp):ytmp)
#define NGB_PERIODIC_LONG_Z(x) (ztmp=fabs(x),(ztmp>boxHalf_Z)?(boxSize_Z-ztmp):ztmp)

#define NEAREST_X(x) (xtmp=(x),(xtmp>boxHalf_X)?(xtmp-boxSize_X):( (xtmp< -boxHalf_X)?(xtmp+boxSize_X):(xtmp) ) )
#define NEAREST_Y(x) (ytmp=(x),(ytmp>boxHalf_Y)?(ytmp-boxSize_Y):( (ytmp< -boxHalf_Y)?(ytmp+boxSize_Y):(ytmp) ) )
#define NEAREST_Z(x) (ztmp=(x),(ztmp>boxHalf_Z)?(ztmp-boxSize_Z):( (ztmp< -boxHalf_Z)?(ztmp+boxSize_Z):(ztmp) ) )

#define WRAP_X(x) (xtmp=(x),(xtmp>boxSize_X)?(xtmp-boxSize_X):( (xtmp< 0)?(xtmp+boxSize_X):(xtmp) ) )
#define WRAP_Y(x) (ytmp=(x),(ytmp>boxSize_Y)?(ytmp-boxSize_Y):( (ytmp< 0)?(ytmp+boxSize_Y):(ytmp) ) )
#define WRAP_Z(x) (ztmp=(x),(ztmp>boxSize_Z)?(ztmp-boxSize_Z):( (ztmp< 0)?(ztmp+boxSize_Z):(ztmp) ) )

#else

#define NGB_PERIODIC_LONG_X(x) fabs(x)
#define NGB_PERIODIC_LONG_Y(x) fabs(x)
#define NGB_PERIODIC_LONG_Z(x) fabs(x)

#define NEAREST_X(x) (x)
#define NEAREST_Y(x) (x)
#define NEAREST_Z(x) (x)

#define WRAP_X(x) (x)
#define WRAP_Y(x) (x)
#define WRAP_Z(x) (x)
#endif

#define FACT1 0.366025403785    /* FACT1 = 0.5 * (sqrt(3)-1) */
#define FAC_TWO_TO_TWO_THIRDS   1.5874011

/*********************************************************/
/*  Global variables                                     */
/*********************************************************/

extern int TimeBinSynchronized[TIMEBINS];
extern struct TimeBinData TimeBinsHydro, TimeBinsGravity;

#ifdef TRACER_PARTICLE
extern struct TimeBinData TimeBinsTracer;
#endif

#ifdef BLACK_HOLES
extern struct TimeBinData TimeBinsBHAccretion;
#endif

#ifdef SINKS
extern struct TimeBinData TimeBinsSinksAccretion;
#endif

#ifdef DUST_LIVE
extern struct TimeBinData TimeBinsDust;
#endif

#ifdef USE_SFR
extern double TimeBinSfr[TIMEBINS];
#endif

#ifdef BLACK_HOLES
extern double TimeBin_BH_mass[TIMEBINS];
extern double TimeBin_BH_dynamicalmass[TIMEBINS];
extern double TimeBin_BH_Mdot[TIMEBINS];
extern double TimeBin_BH_Medd[TIMEBINS];
#endif


extern int ThisTask;            /**< the number of the local processor  */
extern int NTask;               /**< number of processors */
extern int PTask;               /**< note: NTask = 2^PTask */

extern int ThisNode;            /**< the rank of the current compute node  */
extern int NumNodes;            /**< the number of compute nodes used  */
extern int MinTasksPerNode;     /**< the minimum number of MPI tasks that is found on any of the nodes  */
extern int MaxTasksPerNode;     /**< the maximum number of MPI tasks that is found on any of the nodes  */
extern int TasksInThisNode;     /**< number of MPI tasks on  current compute node */
extern int RankInThisNode;      /**< rank of the MPI task on the current compute node */
extern long long MemoryOnNode;

extern double CPUThisRun;       /**< Sums CPU time of current process */


extern int MaxTopNodes;         /**< Maximum number of nodes in the top-level tree used for domain decomposition */

extern int RestartFlag;         /**< taken from command line used to start code. 0 is normal start-up from
				   initial conditions, 1 is resuming a run from a set of restart files, while 2
				   marks a restart from a snapshot file. */
extern int RestartSnapNum;

extern int TakeLevel;

extern int TagOffset;

extern int Argc;
extern char **Argv;

extern double CPU_Step[CPU_LAST];
extern double CPU_Step_Stored[CPU_LAST];


extern double WallclockTime;    /**< This holds the last wallclock time measurement for timings measurements */
extern double StartOfRun;       /**< This stores the time of the start of the run for evaluating the elapsed time */

extern size_t AllocatedBytes;
extern size_t FreeBytes;

extern char DumpFlag;
extern char DumpFlagNextSnap;

extern int  FlagNyt;

extern int NumPart;             /**< number of particles on the LOCAL processor */
extern int NumGas;              /**< number of gas particles on the LOCAL processor  */
#ifdef TRACER_MC
extern int N_tracer;            /**< number of tracer particles on the LOCAL processor  */
#endif
#ifdef GFM
extern int N_star;              /**< number of star particles on the LOCAL processor  */
#endif
#ifdef BLACK_HOLES
extern int NumBHs;               /**< number of BH particles on the LOCAL processor  */
#endif
#ifdef DUST_LIVE
extern int N_dust;               /**< number of dust particles on the LOCAL processor */
#endif

#ifdef INJECT_TRACER_INTO_SN
extern MyIDType MaxTracerID;   /**< stores the maximum ID the tracer reached on any processor */
#endif

extern gsl_rng *random_generator;          /**< a random number generator  */
extern gsl_rng *random_generator_aux;      /**< an auxialiary random number generator for use if one doesn't want to influence the main code's random numbers  */

#ifdef USE_SFR
extern int Stars_converted;     /**< current number of star particles in gas particle block */
#endif

#ifdef TOLERATE_WRITE_ERROR
extern int WriteErrorFlag;
extern char AlternativeOutputDir[MAXLEN_PATH];
#endif

#ifdef GFM_AGN_RADIATION
extern int CellsWithAGNRadiation;
#endif

#if defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS)
extern double WindEnergy_Should, WindEnergy_Is;
#endif
#ifdef BLACK_HOLES
extern double AGNEnergyEM_Is, AGNEnergyEMobs_Is, AGNEnergyM_Should, AGNEnergyM_Is, AGNEnergyT_Should, AGNEnergyT_Is;
#endif

#ifdef SUBBOX_SNAPSHOTS
extern int Nsubboxes;                   /**< the number of subboxes for frequent snapshots */
extern double *SubboxXmin, *SubboxXmax, *SubboxYmin, *SubboxYmax, *SubboxZmin, *SubboxZmax;  /**< the coordinates of subboxes */
#endif

#ifdef EOS_OPAL
extern double opal_rhomax, opal_rhomin;  /**< boundaries of eos table */
#endif

#ifdef GRACKLE
extern code_units my_grackle_units;
#endif

extern double EgyInjection;

extern double TimeOfLastDomainConstruction;     /**< holds what it says */

/** Buffer to hold indices of neighbours retrieved by the neighbour
  search routines. Usually of size Numpart. */
//extern int *Ngblist;


extern double DomainCorner[3], DomainCenter[3], DomainLen, DomainFac;
extern double DomainInverseLen, DomainBigFac;
extern int *DomainStartList, *DomainEndList;
extern double *DomainCost, *TaskCost;
extern int *DomainCount, *TaskCount;
extern struct no_list_data
{
  int task;
  int no;
  int domainCount;
  double domainCost;
}
 *ListNoData;


extern int domain_bintolevel[TIMEBINS];
extern int domain_refbin[TIMEBINS];
extern int domain_grav_weight[TIMEBINS];
extern int domain_hydro_weight[TIMEBINS];
extern int domain_to_be_balanced[TIMEBINS];

/** Array of task numbers holding the respective top-level nodes. For
    the topnodes entries, it is indexed by the Leaf member, for
    pseudoparticles it is indexed by the node
    number-MaxPart-MaxNodes.  */
extern int *DomainTask;
extern int *DomainNewTask;

/** Array of indices of the main tree nodes that are identical to the
    top-level nodes. For the topnodes entries, it is indexed by the
    Leaf member, for pseudoparticles it is indexed by the node
    number-MaxPart-MaxNodes. */
extern int *DomainNodeIndex;

#ifdef RT_ADVECT
double rt_vec[RT_N_DIR][3];
double rt_vec_new[RT_N_DIR][3];
#endif

#ifdef OTVET
double otvet_sigma_HI[OT_N_BINS];
double otvet_sigma_HeI[OT_N_BINS];
double otvet_sigma_HeII[OT_N_BINS];
double G_HI[OT_N_BINS];
double G_HeI[OT_N_BINS];
double G_HeII[OT_N_BINS];
double lum[OT_N_BINS];
double nu[OT_N_BINS];
#endif


#ifdef MRT
#ifndef MRT_NO_UV
double mrt_sigma_HI[UV_BINS];
double mrt_sigma_HeI[UV_BINS];
double mrt_sigma_HeII[UV_BINS];
double G_HI[UV_BINS];
double G_HeI[UV_BINS];
double G_HeII[UV_BINS];
double P_HI[UV_BINS];
double P_HeI[UV_BINS];
double P_HeII[UV_BINS];
double lum[UV_BINS];
double nu[UV_BINS+1];
#endif
double c_internal_units ;
#ifdef MRT_IR_LTE
double radiation_constant ;
#endif
#endif



extern peanokey *Key, *KeySorted;

/** The top node structure is an octree used for encoding the domain
    decomposition. Its leaf nodes are the units into which the domain
    is decomposed. */
extern struct topnode_data
{
  peanokey Size;
  peanokey StartKey;
  long long Count;
  /** The index of the first daughter node. The remaining 7 follow
      sequentially, I think. */
  int Daughter;
  /** The index of this topnode in the DomainTask etc arrays. Is this
      only valid for topnodes that have daughter=-1, i.e. the actual
      leaves? */
  int Leaf;
  unsigned char MortonToPeanoSubnode[8];
} *TopNodes;

extern int NTopnodes, NTopleaves;


/** Variables for gravitational tree */
extern int Tree_MaxPart;
extern int Tree_NumNodes;
extern int Tree_MaxNodes;
extern int Tree_FirstNonTopLevelNode;
extern int Tree_NumPartImported;
extern int Tree_NumPartExported;
extern int Tree_ImportedNodeOffset;
extern int Tree_NextFreeNode;

extern int *Tree_ResultIndexList;
extern int *Tree_Task_list;
extern MyDouble *Tree_Pos_list;
extern unsigned long long *Tree_IntPos_list;

extern struct treepoint_data
{
  MyDouble Pos[3];
  unsigned long long IntPos[3];
  MyDouble Mass;
  float OldAcc;
  int index;
  int th;
  unsigned char level;
  unsigned char Type;
  unsigned char SofteningType : 7;
#ifndef HIERARCHICAL_GRAVITY
  unsigned char ActiveFlag : 1;
#endif

#if defined(SUBFIND) && defined(SUBFIND_EXTENDED_PROPERTIES)
  MyFloat GroupRad;
  int GrNr;
#endif
#ifdef BLACK_HOLES
  int AuxDataIndex;
#if defined(BLACK_HOLES) && (defined(MEASURE_POTMIN_AROUND_BH) || defined(BH_NEW_CENTERING))
  MySingle Potential;
#endif
#ifdef TRACER_MC
  int origTask; /** the task on which this particle resides in P */
#endif
#endif
#if defined(OTVET) || defined(TREE_RAD) || defined(TREE_RAD_H2) || defined (TREE_RAD_CO)
  MyFloat Hsml;
#endif
#if defined(TREE_RAD_H2) || defined (TREE_RAD_CO) || defined(TREE_RAD_VEL)
  double TracAbund[SGCHEM_NUM_ADVECTED_SPECIES];
#endif
#ifdef TREE_RAD_VEL
  MyFloat Vel[3];
#endif
#if defined(RADPRESS_OPT_THIN) && defined(GFM)
  MyFloat BirthTime;
#endif
#if defined(RADCOOL)
  MyFloat BirthTime;
  double young_stellar_mass;
  double old_stellar_mass;
  double young_stellar_s[3];
  double old_stellar_s[3];
#ifdef RADCOOL_HOTHALO
  MyFloat Utherm;
  MyFloat Density;
  MyFloat Temperature;
#ifdef COOLING
  MyFloat Ne;
#endif
#ifdef RADCOOL_HOTHALO_METAL_BOOST
  MyFloat Metallicity;
#endif
  double T6_gas_mass;
  double T7_gas_mass;
  double T8_gas_mass;
  double T6_gas_s[3];
  double T7_gas_s[3];
  double T8_gas_s[3];
#endif
#endif

#ifdef SIDM
  MyFloat Vel[3];
  MyDouble sidm_Hsml;
  unsigned char sidm_TimeBin;
  MyIDType sidm_ID;
  int sidm_State;
#endif

#ifdef ISM
  MyFloat Metallicity;
  MyFloat BirthTime;
  MyFloat Sfr;
#endif

#ifdef DUST_LIVE
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  MyFloat DustHsml;
#endif
#endif
}
 *Tree_Points;


#ifdef BLACK_HOLES
extern struct treepoint_aux_bh_data
{
  MyFloat BH_Mass;
  MyFloat BH_CumMass_QM;
  MyFloat BH_CumEgy_QM;
  MyFloat BH_CumMass_RM;
  MyFloat BH_CumEgy_RM;
  MyFloat BH_Hsml;
  MyFloat BH_U;
  MyFloat Vel[3];
  MyIDType ID;
  MyIDType SwallowID;
  unsigned char TimeBin;
  int BH_CountProgs;
#ifdef DRAINGAS
  MyDouble DrainBucketMass;
#endif
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#endif
#ifdef BH_NF_RADIO
  MyFloat BH_RadioEgyFeedback;
  MyFloat BH_HaloVvir;
  MyFloat BH_Mdot_quasar;
  MyFloat BH_Mdot_radio;
  MyFloat BH_XrayLum;
  MyFloat BH_RadioLum;
#endif
#if defined(MASSIVE_SEEDS_MERGER)
  MyFloat HostHaloMass;
#endif
#ifdef BH_SPIN_EVOLUTION
  MyFloat BH_SpinParameter;
  MyFloat BH_SpinOrientation[3];
#endif
}
 *Tree_AuxBH_Points;

extern int Tree_NumBHExported, Tree_NumBHImported;

#define TBPP(i) Tree_AuxBH_Points[Tree_Points[(i)].AuxDataIndex]

#endif



extern struct resultsactiveimported_data
{
  MyFloat GravAccel[3];
#ifdef EVALPOTENTIAL
  MyFloat Potential;
#endif
#ifdef RADCOOL
  double Phios, Phins;
#ifdef RADCOOL_HOTHALO
  double PhiT6, PhiT7, PhiT8;
#endif
#endif
#ifdef TREE_RAD
#define NPIX 12*NSIDE*NSIDE
  double Projection[NPIX];
#ifdef TREE_RAD_H2
  double ProjectionH2[NPIX];
#endif
#ifdef TREE_RAD_CO
  double ProjectionCO[NPIX];
#endif
#endif /* TREE_RAD */
#ifdef OUTPUTGRAVINTERACTIONS
  int GravInteractions;
#endif
#ifdef MODGRAV
  MyFloat ModgravAccel[3];
#endif
  int index;
}
 *Tree_ResultsActiveImported;







extern char ParameterFile[MAXLEN_PATH]; /**< file name of parameterfile used for starting the simulation */

extern FILE *FdInfo,            /**< file handle for info.txt log-file. */
 *FdEnergy,                     /**< file handle for energy.txt log-file. */
 *FdTimings,                    /**< file handle for timings.txt log-file. */
 *FdBalance,                    /**< file handle for balance.txt log-file. */
 *FdTimebin,                    /**< file handle for timebins.txt log-file. */
 *FdDomain,                     /**< file handle for domain.txt log-file. */
 *FdMemory,                     /**< file handle for memory.txt log-file. */
 *FdCPU;                        /**< file handle for cpu.txt log-file. */


#if defined(LOCAL_FEEDBACK)
extern FILE *FdLocalFeedback;
#endif

#ifdef DETAILEDTIMINGS
extern FILE *FdDetailed;
#endif


#ifdef OUTPUT_CPU_CSV
extern FILE *FdCPUCSV;          /**< file handle for cpu.csv log-file. Used if the cpu log is printed in csv format as well. */
#endif

#ifdef RESTART_DEBUG
extern FILE *FdRestartTest;
#endif

#ifdef USE_SFR
extern FILE *FdSfr;             /**< file handle for sfr.txt log-file. */
#endif

#ifdef BLACKHOLE_POTMIN_DIAGNOSTIC
extern FILE *FdBHDiag;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
extern FILE *FdTurb;
#endif

#ifdef RT_ADVECT
extern FILE *FdRad;             /**< file handle for radtransfer.txt log-file. */
#endif


#ifdef BLACK_HOLES
extern FILE *FdBlackHoles;      /**< file handle for blackholes.txt log-file. */
extern FILE *FdBlackHolesDetails;
extern FILE *FdBlackHolesMergers;
#ifdef BH_SPIN_EVOLUTION
extern FILE *FdBlackHolesSpin;
#endif
#ifdef  BH_NEW_CENTERING
extern FILE *FdBlackHolesRepos;
#endif
#ifdef BH_BIPOLAR_FEEDBACK
extern FILE *FdBlackHolesBipolar;
#endif
#endif

#ifdef GFM_STELLAR_EVOLUTION
extern FILE *FdMetalsGas;                      /**< file handle for metals_gas.txt log-file. */
extern FILE *FdMetalsStars;                    /**< file handle for metals_star.txt log-file. */
extern FILE *FdMetalsTot;                      /**< file handle for metals_tot.txt log-file. */
extern FILE *FdSN;                             /**< file handle for SN.txt log-file. */
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
extern FILE *FdFeedback;        /*!< file handle for stellar-feedback.txt log-file. */
extern FILE *FdGasFeed;         /*!< file handle for gas-feedback.txt log-file. */
#endif

#ifdef FORCETEST
extern FILE *FdForceTest;       /*!< file handle for forcetest.txt log-file. */
#endif


#ifdef DARKENERGY
extern FILE *FdDE;              /**< file handle for darkenergy.txt log-file. */
#endif

#ifdef OTVET
extern FILE *FdOTVET;           /*!< file handle for otvet.txt log-file. */
#ifdef OTVET_MULTI_FREQUENCY
extern FILE *FdOTVETStar;       /*!< file handle for lum_star.txt log-file. */
#endif
#endif

#ifdef MRT
extern FILE *FdMRT;           /*!< file handle for otvet.txt log-file. */
#ifdef MRT_MULTI_FREQUENCY
extern FILE *FdMRTStar;       /*!< file handle for lum_star.txt log-file. */
#endif
#endif


#ifdef BINARYLOG
extern FILE *FdBinary;		/**< file handle for binary.txt log-file. */
#endif

#ifdef NUCLEAR_NETWORK
extern FILE *FdNetwork;
#endif

#ifdef GENERAL_RELATIVITY
extern FILE *FdGR;
#endif

#ifdef SINK_PARTICLES
extern FILE *FdSinkPart;
#endif

#ifdef COSMIC_RAYS
extern FILE *FdCREnergy;
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
extern FILE *FdDust;
#endif

/** Determines whether various dump files are written. Normally true,
    set to false by Sunrise to avoid creating them. */
extern int WriteMiscFiles;







extern void *CommBuffer;        /**< points to communication buffer, which is used at a few places */

#ifdef TRACER_PARTICLE
extern int TracerPartTmaxIndex;         // store maximum temperature for velocity tracers
extern int TracerPartTmaxTimeIndex;     // store time of maximum temperature for velocity tracers
extern int TracerPartTmaxRhoIndex;      // store density at time of maximum temperature for velocity tracers
extern int TracerPartRhomaxIndex;       // store maximum density for velocity tracers
extern int TracerPartRhomaxTimeIndex;   // store time of maximum density for velocity tracers
extern int TracerPartMachmaxIndex;      // store maximum Mach from the riemann solution in SphP for the velocity tracers
extern int TracerPartEntmaxIndex;       // store maximum entropy (= P/rho^gamma) for velocity tracers
extern int TracerPartEntmaxTimeIndex;   // store time of maximum entropy (= P/rho^gamma) for velocity tracers
#endif

#ifdef TRACER_MC
extern int TracerMCTmaxIndex;           // store maximum temperature for MC tracers
extern int TracerMCTmaxTimeIndex;       // store time of maximum temperature for MC tracers
extern int TracerMCTmaxRhoIndex;        // store density at time of maximum temperature for MC tracers
extern int TracerMCRhomaxIndex;         // store maximum density for MC tracers
extern int TracerMCRhomaxTimeIndex;     // store time of maximum density for MC tracers
extern int TracerMCMachmaxIndex;        // store maximum Mach from the riemann solution in SphP for the MC tracers
extern int TracerMCEntmaxIndex;         // store maximum entropy (= P/rho^gamma) for MC tracers
extern int TracerMCEntmaxTimeIndex;     // store time of maximum entropy (= P/rho^gamma) for MC tracers
extern int TracerMCLastStarTimeIndex;   // store the last time a MC tracer belonged to a star/wind particle
extern int TracerMCWindCounterIndex;    // store the number of times a tracer has been kicked to the wind (GFM_WINDS)
extern int TracerMCExchangeCounterIndex;        // store the number of times a tracer has been moved between particles/cells
extern int TracerMCExchangeDistanceIndex;       // store the sum of the radii of the cells between which the tracer was exchanged
extern int TracerMCExchangeDistanceErrorIndex;  // store an estimate of the error of the trajectory of the tracer
extern int TracerMCShockMachMaxIndex;   // store maximum Mach number encountered from the on-the-fly shock finder
#endif

/** Data which is the SAME for all tasks (mostly code parameters read
 * from the parameter file).  Holding this data in a structure is
 * convenient for writing/reading the restart file, and it allows the
 * introduction of new global variables in a simple way. The only
 * thing to do is to introduce them into this structure.
 */
extern struct global_data_all_processes
{
  long long TotNumPart;         /**<  total particle numbers (global value) */
  long long TotNumGas;          /**<  total gas particle number (global value) */

#ifdef TRACER_MC
  long long N_alltracer_global;   /**< global number of tracer particles (belonging to all particle types)  */
#endif

#ifdef BLACK_HOLES
  int TotNumBHs;
#endif
#ifdef DUST_LIVE
  long long TotN_dust;               /**< global number of dust particles */
#endif

  int MaxPart;                  /**< This gives the maxmimum number of particles that can be stored on one
				   processor. */
  int MaxPartSph;               /**< This gives the maxmimum number of SPH particles that can be stored on one
				   processor. */

#ifdef MRT
  int RTNumSubCycles ;
#ifdef MRT_RIEMANN_HLLE
  char HLLEFile[MAXLEN_PATH];
#endif
#ifdef MRT_IR_GRAIN_KAPPA
  char GrainKappaPath[MAXLEN_PATH];
#endif
#ifdef MRT_BH
  double LogAGNLuminosity;
  double UVLuminosityFraction;
#endif
#endif

#if defined(COOLING) && !defined(GRACKLE)
  char TreecoolFile[MAXLEN_PATH];
#ifdef GFM_AGN_RADIATION
  char TreecoolFileAGN[MAXLEN_PATH];
#endif
#ifdef RADCOOL
  char TreecoolFileRAD[MAXLEN_PATH];
  int NewStarsOn;
  int OldStarsOn;
#ifdef RADCOOL_HOTHALO
  int HotHaloOn;
  double Redshift;
#endif
  int SelfShieldingOn;
  MyFloat SelfShieldingDensity;
#endif
#ifdef UVB_SELF_SHIELDING
  char SelfShieldingFile[MAXLEN_PATH];
#endif
#endif

#ifdef GRACKLE
  int GrackleOn;
  int GrackleRadiativeCooling;
  int GrackleMetalCooling;
  int GrackleUVB;
  char GrackleDataFile[MAXLEN_PATH];
#ifndef METALS
  MyFloat GrackleInitialMetallicity;
#endif
#endif

#ifndef COOLING
#ifdef ATOMIC_DM
  char TreecoolFile[MAXLEN_PATH];
#endif
#endif

#ifdef GFM_AGN_RADIATION
  double SelfShieldingDensity;    /**< hydrogen number density limit for shielding */
  double ObscurationFactor;       /**< obscuration prefactor */
  double ObscurationSlope;        /**< obscuration slope */
#endif

#ifdef TRACER_TRAJECTORY
#ifdef TRACER_TRAJECTORY_GENERATE
  int NumberOfTracersToGenerate;
#else
  char TracerInitFile[MAXLEN_PATH];
#endif
  char TracerOutputFile[MAXLEN_PATH];
  char TracerOutputConfFile[MAXLEN_PATH];

  int NTracerTot;               /* total number of tracers on all tasks */
  int NTracerOutputSteps;
  double *NTracerOutputTimes;
  double *NTracerOutputTimePeriod;
  double TracerNextOutput;
#endif

#ifdef TRACER_MC
  int MaxPartTracer;            /**< This gives the maxmimum number of tracer particles that can be stored on one processor */
  int TracerMCPerCell;          /**< The initial number of MC-tracers per gas cell, specified in the parameter file. */
  double ReferenceTracerMCMass; /**< The mass associated with each MC tracer */
#endif
#ifdef TRACER_PARTICLE
  double MinimumTracerHsml;     /**< The minimum distance used in search of the nearest cell - should be an estimate for the minimum cell size */
#endif
#ifdef GFM
  long long TotN_star;
  int MaxPartStar;              /**< This gives the maxmimum number of star particles that can be stored on one
                                   processor. */
#endif

#ifdef BLACK_HOLES
  int MaxPartBHs;               /**< This gives the maxmimum number of BH particles that can be stored on one
                                   processor. */
#endif

#ifdef DUST_LIVE
  int MaxPartDust;              /**< Maximum number of dust particles that can be stored on one processor. */
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  int TotPartSpecial, MaxPartSpecial;
#endif
#ifdef REFINEMENT_AROUND_DM
  int TotPartDM;
#endif
#ifdef CIRCUMSTELLAR
#ifdef CIRCUMSTELLAR_REFINEMENTS
  double CircumstellarDerefinementDistance;
#endif
#ifdef CIRCUMSTELLAR_IRRADIATION
  int TotPartSources, MaxPartSources;
  double IrradiationTempScaling;
#endif
#endif

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_SINKS)
  double CircumstellarSinkRadius;
#endif

#ifdef BH_THERMALFEEDBACK_ACC
  double BlackholeDeltaTemp;
  double BlackholeDeltaTime;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
  double BHBipolarTheta;
  double BHBipolarEfficiency;
  double BHBipolarColdTemp;
  double BHBipolarColdFraction;
#endif

#ifdef BH_NF_RADIO
  double RadioModeMachnumber;
  double RadioRelativeBubbleSize;
  double RadioRelativeBubbleEnergy;
  double RadioRelativeMaxDist;
  double RadioModeMetallicityInSolar;
#endif

#ifdef GFM_UVB_CORRECTIONS
  double UV_HeII_alpha;
  double UV_HeII_beta;
  double UV_HeII_threshold;
  double UV_HeatBoost;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  int DesNumNgbEnrichment;
  double MaxNumNgbDeviationEnrichment;
  double IMF_MinMass_Msun;      /**< IMF min mass */
  double IMF_MaxMass_Msun;      /**< IFM max mass */
  /* SNIa */
  double SNIa_Rate_Norm;
  double SNIa_Rate_TAU;
  int SNIa_MassTransferOn;
  /* SNII */
  int SNII_MassTransferOn;
  /* AGB flag */
  int AGB_MassTransferOn;
  /*NSNS*/
#ifdef GFM_RPROCESS
  double NSNS_MassPerEvent;
  double NSNS_per_SNIa;
  double NSNS_Rate_TAU;
  int NSNS_MassTransferOn;
#endif
  /* Yield tables path */
  char YieldTablePath[MAXLEN_PATH];
  double SNII_MinMass_Msun;
  double SNII_MaxMass_Msun;
#endif

#ifdef GFM_WIND_ENERGY_METAL_DEPENDENCE
  double WindEnergyReductionFactor;
  double WindEnergyReductionMetallicity;
  double WindEnergyReductionExponent;
#endif

#ifdef GFM_DUST
  /* Dust production efficiency factors. */
  double AGB_Dust_Delta_C;
  double AGB_Dust_Delta_Metal;
  double SN_Dust_Delta_C;
  double SN_Dust_Delta_Metal;

  double Dust_Growth_Tau;       /**< Dust accretion timescale, in Gyr. */

#ifndef GFM_DUST_MRN
  double DustSingleGrainSize; // grain size in mu
#endif

#ifdef GFM_DUST_SPUTTERING
#if (GFM_DUST_SPUTTERING==1)
  double Dust_Sputter_Tau_Fac;
#endif
#endif

#if (GFM_DUST_DESTMODE==1)
  double Dust_Destruction_Tau;  /**< Dust destruction timescale, in Gyr. */
#endif
#endif

#ifdef GFM_PREENRICH
  double PreEnrichTime;
  char PreEnrichAbundanceFile[MAXLEN_PATH];


#ifdef WINDTUNNEL_FIXVARIABLESININJECTIONREGION
  double mass_fractions[GFM_N_CHEM_ELEMENTS];
  double metallicity;
#endif

#endif

#ifdef GFM_COOLING_METAL
  char CoolingTablePath[MAXLEN_PATH];
  double MinMetalTemp;                       /**< no metal line cooling below this value */
#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
  char PhotometricsTablePath[MAXLEN_PATH];
#endif

#if defined(REFINEMENT) || defined(GFM_STELLAR_EVOLUTION) || defined(BLACK_HOLES)
  double ReferenceGasPartMass;
#endif

#ifdef REFINEMENT
  double TargetGasMass;
  double TargetGasMassFactor;
  int RefinementCriterion;
  int DerefinementCriterion;
#endif

#ifdef GMC_REFINEMENT
  int GMCRefCellsPerJeansLength;
  double GMCRefMinDensity;
  double GMCRefMaxDensity;
  double GMCDerefMinDensity;
#endif


#ifdef STICKYFLAGS
#if (REFLECTIVE_X == 2) || (REFLECTIVE_Y == 2) || (REFLECTIVE_Z == 2)
  double StickyLayerMaxDist;
#endif
#endif

  double TotGravCost;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  double AvgType1Mass;
#endif

  double MeanVolume;

  int MultipleDomains;
  double TopNodeFactor;

  int ICFormat;                 /**< selects different versions of IC file-format */

  int SnapFormat;               /**< selects different versions of snapshot file-formats */

  int NumFilesPerSnapshot;      /**< number of files in multi-file snapshot dumps */
  int NumFilesWrittenInParallel;/**< maximum number of files that may be written/read simultaneously when
				   writing/reading restart-files, or when writing snapshot files */

  double TreeAllocFactor;       /**< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double TopNodeAllocFactor;    /**< Each processor allocates a number of nodes which is TreeAllocFactor times
				   the maximum(!) number of particles.  Note: A typical local tree for N
				   particles needs usually about ~0.65*N nodes. */

  double NgbTreeAllocFactor;   /**< Each processor allocates a number of nodes for the neighbor search which is NgbTreeAllocFactor times
                                   the maximum(!) number of gas particles.  Note: A typical local tree for N
                                   particles needs usually about ~0.65*N nodes. */

  int MaxMemSize;               /**< size of maximum memory consumption in MB */

  /* some SPH parameters */

  int DesNumNgb;                /**< Desired number of SPH neighbours */



#if defined (VS_TURB) || defined (AB_TURB)
  double RefDensity;
  double RefEntropy;
  double TurbInjectedEnergy;
  double TurbDissipatedEnergy;
  double SetLastTime;
  double TimeBetTurbSpectrum;
  double TimeNextTurbSpectrum;
  int FileNumberTurbSpectrum;
#endif

#ifdef SUBFIND
  int DesLinkNgb;
  double ErrTolThetaSubfind;
#endif

  double TotCountReducedFluxes;
  double TotCountFluxes;
#ifdef COSMIC_RAYS
  double TotCountCRLimiter;
  double MinimumCREnergyDensity;
#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
  double AccelerationEfficiency;
  double CriticalMachnumber;
#endif

  double TotalCREnergyInjected;
  double TotalCREnergyCooledHadronic;
  double TotalCREnergyCooledCoulomb;
  double TotalCREnergyCooled;

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double TotalCREnergyErrorDiffusion;
  double TotalCREnergyChangeAdiabatic;
  double TotalCREnergyLossSfr;
  double TotalCREnergyUpdatePrims;
#endif
#endif

  double DtDisplacement;

  double MaxNumNgbDeviation;    /**< Maximum allowed deviation neighbour number */

  double InitGasTemp;           /**< may be used to set the temperature in the IC's */
  double InitGasU;              /**< the same, but converted to thermal energy per unit mass */
  double MinGasTemp;            /**< may be used to set a floor for the gas temperature */
  double MinEgySpec;            /**< the minimum allowed temperature expressed as energy per unit mass */

  double MinimumDensityOnStartUp;

#ifdef MIN_METALLICITY_ON_STARTUP
  double MinimumMetallicityOnStartUp;
#endif

  double GasSoftFactor;

  double LimitUBelowThisDensity;
  double LimitUBelowCertainDensityToThisValue;


  /* some force counters  */

  long long TotNumOfForces;     /**< counts total number of force computations  */

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)
  int PicXpixels, PicYpixels;
  int PicXaxis, PicYaxis, PicZaxis;
  double PicXmin, PicXmax, PicYmin, PicYmax, PicZmin, PicZmax;
#endif

#ifdef SUBBOX_SNAPSHOTS
  char SubboxCoordinatesPath[MAXLEN_PATH];
  double SubboxMinTime, SubboxMaxTime;
  int SubboxSyncCounter;
  int SubboxSyncModulo;
  int SubboxNumFilesPerSnapshot;
  int SubboxNumFilesWrittenInParallel;
  int SubboxNumber;
#endif

#ifdef MULTIPLE_RESTARTS
  int RestartFileCount;
#endif

  /* various cosmological factors that are only a function of the current scale factor, and in non-comoving runs are set to 1 */
  double cf_atime, cf_a2inv, cf_a3inv, cf_afac1, cf_afac2, cf_afac3, cf_hubble_a, cf_time_hubble_a, cf_redshift;
  /* Hubble rate at the current time, valid both for comoving and non-comoving integration */
  double cf_H;
  /* Hubble expansion rate, but in non-comoving integration set to zero */
  double cf_Hrate;

  /* system of units  */
  double UnitTime_in_s,         /**< factor to convert internal time unit to seconds/h */
    UnitMass_in_g,              /**< factor to convert internal mass unit to grams/h */
    UnitVelocity_in_cm_per_s,   /**< factor to convert internal velocity unit to cm/sec */
    UnitLength_in_cm,           /**< factor to convert internal length unit to cm/h */
    UnitPressure_in_cgs,        /**< factor to convert internal pressure unit to cgs units (little 'h' still
				   around!) */
    UnitDensity_in_cgs,         /**< factor to convert internal mass density unit to g/cm^3*h^2 */
    UnitCoolingRate_in_cgs,     /**< factor to convert internal cooling rate to cgs units */
    UnitEnergy_in_cgs,          /**< factor to convert internal energy to cgs units */
    UnitTime_in_Megayears,      /**< factor to convert internal time to megayears/h */
    GravityConstantInternal,    /**< If set to zero in the parameterfile, the internal value of the
				   gravitational constant is set to the Newtonian value based on the system of
				   units specified. Otherwise the value provided is taken as internal gravity
				   constant G. */
    G;                          /**< Gravity-constant in internal units */

#ifdef BECDM
  double hbar;
  double mAxion;
  double AxionMassEv;
#endif

  /* Cosmology */

  double Hubble;                /**< Hubble-constant in internal units */
  double Omega0,                /**< matter density in units of the critical density (at z=0) */
    OmegaLambda,                /**< vaccum energy density relative to crictical density (at z=0) */
    OmegaBaryon,                /**< baryon density in units of the critical density (at z=0) */
    HubbleParam;                /**< little `h', i.e. Hubble constant in units of 100 km/s/Mpc.  Only needed to get absolute
				 * physical values for cooling physics
				 */

  double BoxSize;               /**< Boxsize in case periodic boundary conditions are used */

  /* Code options */

  int ComovingIntegrationOn;    /**< flags that comoving integration is enabled */
  int PeriodicBoundariesOn;     /**< flags that periodic boundaries are enabled for gravity */
  int ResubmitOn;               /**< flags that automatic resubmission of job to queue system is enabled */
  int TypeOfOpeningCriterion;   /**< determines tree cell-opening criterion: 0 for Barnes-Hut, 1 for relative
				   criterion */
  int TypeOfTimestepCriterion;  /**< gives type of timestep criterion (only 0 supported right now - unlike
				   gadget-1.1) */
  int OutputListOn;             /**< flags that output times are listed in a specified file */
  int CoolingOn;                /**< flags that cooling is enabled */
  int StarformationOn;          /**< flags that star formation is enabled */

  int NParameters;

  int LowestActiveTimeBin;
  int HighestActiveTimeBin;
  int LowestOccupiedTimeBin;
  int HighestOccupiedTimeBin;
  int LowestOccupiedGravTimeBin;
  int HighestOccupiedGravTimeBin;
  int HighestSynchronizedTimeBin;
  int SmallestTimeBinWithDomainDecomposition;
  double ActivePartFracForNewDomainDecomp;

#ifdef GFM_AGN_RADIATION
  int SmallestTimeBinWithAGNRad;
#endif
  /* parameters determining output frequency */

  int SnapshotFileCount;        /**< number of snapshot that is written next */
#ifdef SUBBOX_SNAPSHOTS
  int SubboxSnapshotFileCount;  /**< number of subbox snapshot that is written next */
#endif
  double TimeBetSnapshot,       /**< simulation time interval between snapshot files */
    TimeOfFirstSnapshot,        /**< simulation time of first snapshot files */
    CpuTimeBetRestartFile,      /**< cpu-time between regularly generated restart files */
    TimeLastRestartFile,        /**< cpu-time when last restart-file was written */
    TimeBetStatistics,          /**< simulation time interval between computations of energy statistics */
    TimeLastStatistics;         /**< simulation time when the energy statistics was computed the last time */
  int NumCurrentTiStep;         /**< counts the number of system steps taken up to this point */

  /* Current time of the simulation, global step, and end of simulation */

  double Time,                  /**< current time of the simulation */
    TimeBegin,                  /**< time of initial conditions of the simulation */
    TimeStep,                   /**< difference between current times of previous and current timestep */
    TimeMax;                    /**< marks the point of time until the simulation is to be evolved */

  /* variables for organizing discrete timeline */

  double Timebase_interval;     /**< factor to convert from floating point time interval to integer timeline */
  integertime Ti_Current;       /**< current time on integer timeline */
  integertime Previous_Ti_Current;
  integertime Ti_nextoutput;    /**< next output time on integer timeline */
  integertime Ti_lastoutput;

  integertime Ti_begstep[TIMEBINS];    /**< marks start of current step of each timebin on integer timeline */

#ifdef RT_ADVECT
  integertime Ti_LastRadTransfer;
#endif

#ifdef RT_SHORT_CHARACTERISTICS
  integertime Radiation_Ti_begstep;
  integertime Radiation_Ti_endstep;
#endif

#ifdef PMGRID
  integertime PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], UpperCorner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];
#if defined(EVALPOTENTIAL) && defined(PMGRID) && defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double MassPMregions[2];
#endif
#endif

  long long GlobalNSynchronizedHydro;
  long long GlobalNSynchronizedGravity;

  int LevelToTimeBin[GRAVCOSTLEVELS];
  int LevelHasBeenMeasured[GRAVCOSTLEVELS];

  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_Sum[CPU_LAST];     /**< sums wallclock time/CPU consumption in whole run */

  /* tree code opening criterion */

  double ErrTolTheta;           /**< BH tree opening angle */
  double ErrTolForceAcc;        /**< parameter for relative opening criterion in tree walk */

  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;     /**< accuracy tolerance parameter \f$ \eta \f$ for timestep criterion. The
				   timesteps is \f$ \Delta t = \sqrt{\frac{2 \eta eps}{a}} \f$ */

  double MinSizeTimestep,       /**< minimum allowed timestep. Normally, the simulation terminates if the
				   timestep determined by the timestep criteria falls below this limit. */
    MaxSizeTimestep;            /**< maximum allowed timestep */

#ifdef TIMESTEP_OUTPUT_LIMIT
  double TimestepOutputLimit;
#endif

#ifdef FORCE_EQUAL_TIMESTEPS
  integertime GlobalTimeStep;
#endif

#ifdef LEGACY_DISPLACEMENT_CONSTRAINT
  double MaxRMSDisplacementFac; /**< this determines a global timestep criterion for cosmological simulations
				   in comoving coordinates.  To this end, the code computes the rms velocity
				   of all particles, and limits the timestep such that the rms displacement
				   is a fraction of the mean particle separation (determined from the
				   particle mass and the cosmological parameters). This parameter specifies
				   this fraction. */
#endif

  double IsoSoundSpeed;

  double CourantFac;            /**< SPH-Courant factor */


#ifdef REGULARIZE_MESH_FACE_ANGLE
  double CellMaxAngleFactor;
#else
  double CellShapingFactor;
#endif
  double CellShapingSpeed;


  int CPU_TimeBinCountMeasurements[TIMEBINS];
  double CPU_TimeBinMeasurements[TIMEBINS][NUMBER_OF_MEASUREMENTS_TO_RECORD];


  /* gravitational and hydrodynamical softening lengths (given in terms of an `equivalent' Plummer softening
   * length)
   *
   */

  int SofteningTypeOfPartType[NTYPES];

  double SofteningComoving[NSOFTTYPES];     /**< comoving gravitational softening lengths for each softeniung type */
  double SofteningMaxPhys[NSOFTTYPES];      /**< maximum physical gravitational softening lengths for each softening type */

  double SofteningTable[NSOFTTYPES + NSOFTTYPES_HYDRO];     /**< current (comoving) gravitational softening lengths for each softening type */
  double ForceSoftening[NSOFTTYPES + NSOFTTYPES_HYDRO + 1];     /**<  current (comoving) gravitational softening lengths, multiplied by a factor 2.8 - at that scale the force is Newtonian */

  /** If particle masses are all equal for one type, the corresponding entry in MassTable is set to this
   *  value, * allowing the size of the snapshot files to be reduced
   */
  double MassTable[NTYPES];

#ifdef ADAPTIVE_HYDRO_SOFTENING
  double MinimumComovingHydroSoftening;
  double AdaptiveHydroSofteningSpacing;
#endif


  /* some filenames */
  char InitCondFile[MAXLEN_PATH], OutputDir[MAXLEN_PATH], SnapshotFileBase[MAXLEN_PATH], ResubmitCommand[MAXLEN_PATH], OutputListFilename[MAXLEN_PATH];
  
#ifdef POWERSPECTRUM_IN_POSTPROCESSING
  char InputFileName[MAXLEN_PATH];
#endif

#ifdef VEL_POWERSPEC_BOX
  double PSCenterX;
  double PSCenterY;
  double PSCenterZ;
  double PSRadius;

  double PSBoxMinX;
  double PSBoxMinY;
  double PSBoxMinZ;
  double PSCellSize;
#endif

  /** table with desired output times */
  double OutputListTimes[MAXLEN_OUTPUTLIST];
  char OutputListFlag[MAXLEN_OUTPUTLIST];
  int OutputListLength;         /**< number of times stored in table of desired output times */

#ifdef VORONOI_FREQUENT_IMAGES
  int ImageCount;
  int ImageTimeBin;
  double TimeBetweenImages;
#endif

#ifdef BLACK_HOLES
  double CumAGNEnergyEM_Is, CumAGNEnergyEMobs_Is, CumAGNEnergyM_Should, CumAGNEnergyM_Is, CumAGNEnergyT_Should, CumAGNEnergyT_Is;
#endif

#ifdef USE_SFR			              /* star formation and feedback sector */
#ifdef QUICK_LYALPHA_LATETIMEONLY
  double TimeOfCoolingStart;
#endif
#if !defined(FM_SFR) && !defined(ISM) && !defined(LOCAL_FEEDBACK)        /* enable Springel & Hernquist model */
  double OverDensThresh;
  double CritOverDensity;
  double TemperatureThresh;
  double CritPhysDensity;
  double PhysDensThresh;
  double EgySpecSN;
  double EgySpecCold;
  double FactorEVP;
  double TempSupernova;
  double TempClouds;
#ifdef STEEPER_SFR_FOR_STARBURST
  double PhysDensThreshStarburst;
  double StarburstPowerLawIndex;
#endif
#ifdef SOFTEREQS
  double FactorForSofterEQS;
  double TempForSofterEQS;
  double UForSofterEQS;
#endif
#ifdef MODIFIED_EOS
  double UthermAtThresh;
  double JoinDens;
  double FactorDensThresh;
  double FactorUthermAtThresh;
  double FactorUthermJoin;
#endif
  double MaxSfrTimescale;
  double FactorSN;
#endif

#ifdef GFM_STELLAR_EVOLUTION
  double AGBMassReleased, SNIaMassReleased, SNIIMassReleased;
#endif

#if defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS)
  double CumWindEnergy_Should, CumWindEnergy_Is;
#endif

#ifdef GFM_WINDS_LOCAL
  double VariableWindVelFactor;
  double WindEnergyIn1e51erg;
  double WindEnergyFactor;
  double WindFreeTravelMaxTimeFactor;   /**< maximum free travel time in units of the Hubble time at the current simulation redshift */
  double WindFreeTravelDensFac;
  double MinWindVel;                    /**< the minimum wind velocity in physical code units */
#endif
#ifdef GFM_WINDS_THERMAL
  double ThermalWindFactor;
#endif
#ifdef GFM_WINDS_THERMAL_NEWDEF
  double ThermalWindFraction;
#endif
#ifdef GFM_WINDS
#ifndef GFM_WINDS_VARIABLE
  double WindEfficiency;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  double WindEnergyIn1e51erg;
#else
  double WindEnergyFraction;
#endif
  double WindEnergyFactor;
#ifdef GFM_CONST_IMF
  double WindEgySpecSN;
#endif
  double WindFreeTravelMaxTimeFactor;   /**< maximum free travel time in units of the Hubble time at the current simulation redshift */
  double WindFreeTravelDensFac;
#ifdef GFM_WINDS_VARIABLE
  double VariableWindVelFactor;         /**< wind velocity in units of the halo escape velcoity */
  double VariableWindSpecMomentum;      /**< momentum available for wind per unit mass of stars formed, in physical code units */
  double MinWindVel;                    /**< the minimum wind velocity in physical code units */
  double HaloConcentrationNorm;         /**< concentration c0 of a halo of unit mass */
  double HaloConcentrationSlope;        /**< slope n of mass concentration relation, namely c = c0 * M_200,crit^n */
#ifdef GFM_WINDS_MASSSCALING
  double VariableWindMassScale;
#endif
#ifdef GFM_WINDS_HUBBLESCALING
  double WindSuppressionRedshift;
#endif
#endif
#endif
#endif

#if defined(GFM_CONST_IMF) && (GFM_CONST_IMF == 1)
  double IMFslope;
#endif

#ifdef GFM_WINDS_STRIPPING
  double WindDumpFactor;
#endif

#ifdef GFM_STELLAR_FEEDBACK
  double EnergyPerSNIa;
  double AGBWindVelocity;
#endif
#if defined(FOF) && (defined(BLACK_HOLES) || defined(GFM_WINDS_VARIABLE) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_WINDS_LOCAL))
  double TimeNextOnTheFlyFoF;
  double TimeBetOnTheFlyFoF;
#endif

#ifdef BLACK_HOLES
#ifdef BH_BONDI_DENSITY
  double BlackHoleAccretionSlope;
#else
  double BlackHoleAccretionFactor;      /**< Fraction of BH bondi accretion rate */
#endif
  double BlackHoleFeedbackFactor;       /**< Fraction of the black luminosity feed into thermal feedback */
  double SeedBlackHoleMass;     /**< Seed black hole mass */
  double MinFoFMassForNewSeed;  /**< Halo mass required before new seed is put in */
  int DesNumNgbBlackHole;
  double BlackHoleMaxAccretionRadius;
  double BlackHoleEddingtonFactor;      /**< Factor above Eddington */
  double BlackHoleRadiativeEfficiency;  /**< Radiative efficiency determined by the spin value, default value is 0.1 */
#ifdef BH_NEW_CENTERING
  double BlackHoleCenteringMassMultiplier;
#endif
#ifdef BH_FRICTION
    double BHFrictionCoefficient;
    double BHFrictionAvgTime;
#endif
#ifdef BH_SPIN_EVOLUTION
    double BHInitialSpin;
    double ShakuraSunyaevParameter;
#endif

#ifdef BH_ADIOS_WIND
  double BlackHoleWindOverBubbleFraction;
  double BlackHoleWindAccretionFactor;
#ifdef BH_ADIOS_WIND_WITH_QUASARTHRESHOLD
  double QuasarThreshold;
#endif
#ifdef BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY
  double RadioFeedbackMinDensityFactor;
#endif
#endif

#if defined(BH_PRESSURE_CRITERION)
  double Ref_BH_Pressure;
  double Ref_BH_Mass;
#endif

#ifdef FOF
  double massDMpart;
#endif
#ifdef BH_BONDI_DISK_VORTICITY
  double DiskVorticityRadius;
#endif

#if defined(BH_ADIOS_WIND) || defined(BH_BUBBLES)
  double RadioFeedbackFactor;
#endif

#if defined(BH_ADIOS_WIND) && defined(BH_ADIOS_DENS_DEP_EFFICIANCY)
  double RadioFeedbackFactorPivotMass;
  double RadioFeedbackFactorSlope;
  double RadioFeedbackFactorRefDensityFactor;
  double RadioFeedbackFactorMaxEfficiency;
#endif

#ifdef BH_ADIOS_RANDOMIZED
  double RadioFeedbackReiorientationFactor;
#endif

#ifdef BH_BUBBLES
  double BubbleDistance;
  double BubbleRadius;
  double BubbleEnergy;
  double BlackHoleRadioTriggeringFactor;
  double DefaultICMDensity;
#ifdef UNIFIED_FEEDBACK
  double RadioThreshold;
#endif
#ifdef BH_MAGNETIC_BUBBLES
  double MagneticEnergyFraction;
#endif
#endif
#ifdef MASSIVE_SEEDS
  double DesNumNgbSeed;
  double SeedMaxAccretionRadius;
#endif
#endif

#ifdef DARKENERGY
  double DarkEnergyParam;       /**< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[MAXLEN_PATH];     /**< tabelized w for equation of state */
#ifdef TIMEDEPGRAV
  double Gini;
#endif
#endif
#endif

#ifdef RELAXOBJECT
  double RelaxBaseFac;
  double RelaxFac;
#endif

#ifdef RELAXOBJECT_COOLING
  double RelaxTemperature;
#endif

#ifdef RELAXOBJECT_COOLING2
  double TempCore;
  double TempShell;
  double ShellBaseRadius;
  double TempMin; // minimal temperature at outer bound of shell
  double RTempMin; // position of minimal temperature
  double BackgroundTemp;
#endif

#ifdef RELAXOBJECT_BINARY
  double Omega;
#endif

#ifdef INSPIRAL
  double InspiralVelocity;
  double WD1_Acc[3];
  double WD2_Acc[3];
#endif

#ifdef EOS_DEGENERATE
  char EosTable[MAXLEN_PATH];
  char EosSpecies[MAXLEN_PATH];
#endif

#ifdef EOS_OPAL
  char EosOpalTable[MAXLEN_PATH];
#endif

#ifdef SGCHEM
  int SGChemConstInitAbundances;
  double SGChemInitH2Abund;
  double SGChemInitHPAbund;
  double SGChemInitCPAbund;
  double SGChemInitCOAbund;
  double SGChemInitCHxAbund;
  double SGChemInitOHxAbund;
  double SGChemInitHCOPAbund;
  double SGChemInitHePAbund;
  double SGChemInitMPAbund;
  double SGChemInitDIIAbund;
  double SGChemInitHDAbund;
  double SGChemInitHeIIIAbund;

#ifndef SGCHEM_VARIABLE_Z
  double CarbAbund;
  double OxyAbund;
  double MAbund;
  double ZAtom;
  double DustToGasRatio;
#endif
  double DeutAbund;
  double InitDustTemp;
  double UVFieldStrength;
  long LWBGType;
  double LWBGStartRedsh;
  double CosmicRayIonRate;
  double InitRedshift;
  double ExternalDustExtinction;
  double H2FormEx;
  double H2FormKin;
  long PhotoApprox;
  long ISRFOption;
  long AtomicCoolOption;
  long H2OpacityOption;
#endif

#ifdef SNE_FEEDBACK
  long SNEInjectionCriterion;
  double SNETargetMass;
  double SNEPeriodInYears;
  long SNESeed;
  int SNEMinimalParticleNumber;

#ifdef INJECT_TRACER_INTO_SN
  long SNETracersForEachSn;
  int SNETracerBitmask;
#endif

#ifdef CLUSTERED_SNE 
  double SNEClusterTfin;
  double SNEClusterTinit;
  long SNENumber;
#endif

#endif /*SNE_FEEDBACK*/

#ifdef REFINE_ONLY_WITH_TRACER
  double MaxTracerVolume;
  double MinTracerVolume;
#endif

#ifdef TREE_RAD
  double TreeColMaxDistance;  /* Maximum distance to consider when computing column densities using TreeCol (in code units) */
#endif

#ifdef TREE_RAD_VEL
  double FracOverlap;
#endif

#ifdef NUCLEAR_NETWORK
  char NetworkRates[MAXLEN_PATH];
  char NetworkPartFunc[MAXLEN_PATH];
  char NetworkMasses[MAXLEN_PATH];
  char NetworkWeakrates[MAXLEN_PATH];
  struct network_data nd;
  struct network_workspace nw[NUM_THREADS];
  double NetworkTempThreshold;
#ifdef NUCLEAR_NETWORK_TIMESTEP_LIMITER
  double NuclearNetworkMaxEnergyDiff;
#endif
#ifdef NUCLEAR_NETWORK_LIMIT_COMPOSITION_CHANGE
  double NuclearNetworkMaxCompositionChange;
#endif
#endif

#ifdef NETWORK_NSE
  double NetworkNSEThreshold;
  struct network_data nd_nse;
  struct network_workspace nw_nse[NUM_THREADS];
#endif

#ifdef MHD_POWELL
  double Powell_Momentum[3];
  double Powell_Angular_Momentum[3];
  double Powell_Energy;
#endif

#ifdef MHD_CT
  double CT_mean_Bx;
  double CT_mean_By;
  double CT_mean_Bz;
#endif

#ifdef MHD_SEEDFIELD
  int B_dir;                    /* flags for direction: x = 1, y = 2, z = 4 */
  double B_value;               /* value for the chosen component(s) of the magnetic field */
#endif

#ifdef MHD_SEEDPSPEC
  double B_pspec_slope;
  double B_pspec_ampl;
  double B_pspec_kcut;           /* k=1 corresponds to box length */
  double B_pspec_helical;        /* 0 non-helical, 1 max helical */
  double B_spec_kmax; /*maxminum k */
#endif

#ifdef VISCOSITY
#ifdef GLOBAL_VISCOSITY
  double dyn_visc;
  double bulk_visc;
#else
#ifdef USE_KINEMATIC_VISCOSITY
  double KinematicViscosity;
#else
#ifdef ALPHA_VISCOSITY
  double AlphaCoefficient;
#endif
#endif
#endif
#endif

#ifdef THERMAL_CONDUCTION
  double ThermalConductivity;
#endif


#ifdef TRACER_DIFFUSION
  double TracerDiffusivity;
#endif

MyIDType MaxID;

#ifdef WINDTUNNEL
  double InjectionDensity;
  double InjectionVelocity;
  double InjectionUtherm;
  double InjectionRegion;
  double InjectionVolume;
#ifdef WINDTUNNEL_EXTERNAL_SOURCE
  char WindTunnelExternalSourceFile[MAXLEN_PATH];
#endif
#endif

#ifdef SPECIAL_BOUNDARY
  double BoundaryLayerScaleFactor;
  double SpecialBoundarySpeed;
  int SpecialBoundaryMotion;
  int SpecialBoundaryType;
  double OutflowPressure;
#ifdef THERMAL_CONDUCTION
  double BoundaryTemperature;
#endif
#endif


#if defined(COAXIAL_BOUNDARIES) && !defined(CIRCUMSTELLAR)
  double inner_radius;
  double outer_radius;
  double omega_in;
  double omega_out;
#endif

#if defined(SPECIAL_BOUNDARY) && defined(CIRCUMSTELLAR_WBOUNDARIES)
  double inner_radius;
  double outer_radius;
  double EvanescentBoundaryStrength;
  double CircumstellarBoundaryDensity;
#endif

#ifdef CENTRAL_MASS_POTENTIAL
  double CentralMass;
  double SofteningCentral;
#endif

#if defined(ACCRETE_ONTO_CENTRAL_POTENTIAL) && defined(CENTRAL_MASS_POTENTIAL)
  double CentralAccretionRadius;
#endif

#ifdef PERTURB_VELOCITIES
  double VelocityPerturbation;
#endif

#ifdef STAR_PLANET_POTENTIAL
  double MassRatio;
  double SofteningPlanet;
  double PlanetGrowthTime;
#endif

#ifdef BINARY_POTENTIAL
  double BinaryMassRatio;
  double BinarySoftening;
  double BinaryGrowthTime;
  double BinaryEccentricity;
  int BinaryBarycentricCoord;
#endif

#ifdef MHD_DIVBCLEANING
  double c_h;
#endif

#ifdef AMR
  int MinRefLevel;	//minimum allowed refinement level
  int MaxRefLevel;
  int AMRMeshSmoothing;
#endif

#ifdef AB_TURB
  double StDecay;
  double StEnergy;
  double StDtFreq;
  double StKmin;
  double StKmax;
  double StSolWeight;
  double StAmplFac;

  int StSpectForm;
  int StSeed;

#endif

#ifdef PREHEATING
  unsigned char FlagPreheating;
  double TimePreheating;
  double TempPreheating;
#endif

#ifdef ADJ_BOX_POWERSPEC
  double BoxWidth;              /*lenght of the transformation box side */
  double BoxCenter_x;           /*x coordinate of the box center */
  double BoxCenter_y;           /*y coordinate of the box center */
  double BoxCenter_z;           /*z coordinate of the box center */
  int FourierGrid;              /*dimension of the Fourier transform (actual size is FourierGrid^3) */
#endif

#ifdef REFINEMENT_VOLUME_LIMIT
  double MaxVolumeDiff;
  double MinVolume;
  double MaxVolume;
#endif

#ifdef REFINEMENT_BY_DENSITY
  double MinimumDensityForRefinement;
  double MinimumVolumeForDensityRefinement;
#endif

#ifdef REFINEMENT_AROUND_BH
#ifdef REFINEMENT_AROUND_BH_FIXED
  double RefBHRadius;           /* refinement region in code units */
  double RefBHMaxCellRadius;        /* in code units */
  double RefBHMinCellRadius;        /* in code units */
#else
  double RefBHRadiusHSML;           /* refinement region in units of hsml */
  double RefBHMaxCellRadiusHSML;    /* in units of hsml */
  double RefBHMinCellRadiusRBondi;  /* in units of rbondi */
#endif
  double RefBHMinCellMass;          /* do not refine below this mass */
  double RefBHLowerFactorC;         /* provides lower bound for refined cell sizes */
#endif

#ifdef REFINEMENT_AROUND_DM
  double RefinementCellsPerSoftening;
#endif

#ifdef FM_SFR
  double DensThreshold;          /**< density threshold for star formation in cm^-3 h^2*/
  double CritOverDensity;
  double OverDensThresh;
  double SfrEfficiency;
  double FactorSN;
  double MaxSfrTimescale;
#ifdef USE_POLYTROPIC_EQSTATE
  double UthermThreshold;
#endif
#endif
#ifdef FM_STAR_FEEDBACK
  double PhotoionizationGasTemp; /**< Minimum Temperature for gas that is photoionized    */
  double PhotoionizationEgySpec; /**< Equiv minimum energy for gas that is photoionized   */
  double FeedbackEfficiency;     /**< fraction of the total SN energy that goes in feedback   */
  double EtaKineticEnergy;       /**< fraction of feedback energy that goes in kinetic energy */
#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB) && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK)
  double FeedbackVelocity;       /**< the (constant) kick velocity indiced by stellar feedback */
#endif
#ifdef INSTANTANEOUS_DEPOSITION
  double MinFeedAge;
  double SNaePerUnitMass;
  int FeedbackInjectionEvents;
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
  double MomentumFactor;
#endif
#endif
#ifdef FM_EARLY_STAR_FEEDBACK
  double EarlyFeedbackEfficiency;     /**< fraction of the total early feedback energy that goes in feedback   */
  double EarlyEtaKineticEnergy;       /**< fraction of feedback energy that goes in kinetic energy */
                                          /*  double EarlyFeedbackVelocity;  *//**< the (constant) kick velocity indiced by early stellar feedback,deactivated now, using local escape speed */
#endif


#if defined(CONDUCTION) || defined(MONOTONE_CONDUCTION)                //VITALI
  double ConductionCoeff;       /*!< Thermal Conductivity */
  double ConductionEfficiency;  /* Fraction factor f of Spitzer Conductivity  */
  int MaxConductionSubCycles;
#ifdef CONDUCTION_SATURATION
  double ElectronFreePathFactor;        /*!< Factor to get electron mean free path */
#endif
  integertime conduction_Ti_endstep, conduction_Ti_begstep;
  double dt_conduction, dt_max_conduction;
#ifdef RESTRICT_KAPPA
  double Chi_max_int, MaxDiffusivity ;
#endif
#endif

#ifdef TURBULENT_METALDIFFUSION
#ifndef GFM_STELLAR_EVOLUTION
#define GFM_N_CHEM_ELEMENTS 9
#endif
  integertime metaldiff_Ti_endstep, metaldiff_Ti_begstep;
  double dt_metaldiff;
  double metaldiff_coeff;
  double Courantmetaldiff;
  double metaldiff_kappamax;
#endif

#ifdef DELAYED_COOLING_TURB
  double DissipationTime;
  double SigmaThreshold;
#endif

#ifdef FM_RADIATION_FEEDBACK
  double DustOppacityRadiationFeedback; /* Dust oppacity in gr cm^-2, typically ~5  */
  double InputTimeHeatRadiationFeedback;        /* Maximum time for input of heat due to young stars. Typically 10 Myr */
  double InputTimeMomRadiationFeedback; /* Momentum due to young stars is inputed only once at this specific time. Typically 3 Myr */
  double LumToMassRatioRadiationFeedback;       /* Luminosity to mass ratio for young stars. Typically 1000 */
  double RadiationFeedbackAvgPhotonEnergyineV;  /* Average photon energy for emission from young stars... ~17-20 eV */
#endif

#ifdef ISM
  double DensThreshold;          /**< density threshold for star formation in cm^-3 h^2*/
  double CritOverDensity;
  double OverDensThresh;
  double SfrEfficiency;

  double Kappa_IR;
  double WindMomentumBoost;
#endif

#ifdef ISM_HII_HEATING
  double HIIRegion_fLum_Coupled;
#endif

#ifdef REDUCE_FLUSH
  double FlushCpuTimeDiff;
  double FlushLast;
#endif

#ifdef TILE_ICS
  int TileICsFactor;
#endif

#ifdef DVR_RENDER
#if (DVR_RENDER==1)
  double DvrRenderTimeIntverallInGyr;
  int DvrPixelX, DvrPixelY;
#endif
  int DvrFrameFileNum;
  double DvrTauScaleFactor;
  double DvrTauFloor;
  char DvrOutputDir[MAXLEN_PATH];
#endif

#ifdef FLD
  integertime fld_Radiation_Ti_begstep;
  integertime fld_Radiation_Ti_endstep;

#ifdef FLD_TEST_BOUNDARY
  double fld_Flux;
  double fld_n_gamma;
  double fld_density;
  double fld_u;
#endif

#ifdef FLD_CONST_KAPPA
  double Kappa_R;
  double Kappa_P;
#endif

#ifdef FLD_MARSHAK
  double Epsilon;
#endif
#endif

#ifdef OTVET
  double star_Teff;
  integertime otvet_Radiation_Ti_begstep;
  integertime otvet_Radiation_Ti_endstep;
  double IonizingLumPerSolarMass;       /*For now works only on stars as sources */
  double OtvetNumNgbSource;
  double OtvetMaxNumNgbDeviationSource;
#endif
#ifdef RADPRESS_OPT_THIN
//  double IonizingLumPerSolarMass;          /*For now works only on stars as sources*/
  double RadPressure_MaxAge;    /* in Gyr */
#endif

#ifdef RADCOOL
#define TIMEON_NEWSTARS 0.01    /* Time (in Gyrs) till which the new stars provide ionizing radiation */
#define TIMEON_OLDSTARS 0.20    /*Time (in Gyrs) at which the ionizing radiation from old stars turn on */
#ifdef RADCOOL_HOTHALO
#define TLOGMIN6        5.5     /*Temperature ranges in log10(T [K]) */
#define TLOGMAX6        6.5
#define TLOGMIN7        6.5
#define TLOGMAX7        7.5
#define TLOGMIN8        7.5
#define TLOGMAX8        8.5
#ifdef RADCOOL_HOTHALO_METAL_BOOST
#define T6METALBOOSTFACTOR  28.37
#define T7METALBOOSTFACTOR  2.8435
#define T8METALBOOSTFACTOR  0.291
#define GFM_SOLAR_METALLICITY 0.0127
#endif
#endif
#endif


#ifdef SIDM
  double SIDMDesNumNgb, SIDMMaxNumNgbDeviation;
  double DtimeFac, DtimeFacLim;
  long long sidm_ShouldScatter[SIDM_REACTIONS], sidm_Scatters[SIDM_REACTIONS], sidm_Rejected[SIDM_REACTIONS], sidm_EnergyForbidden[SIDM_REACTIONS];
  double sidm_EnergyInjected[SIDM_REACTIONS], sidm_ScatteredMass[SIDM_REACTIONS];
  double sidm_EnergyInjectedCheck_sidm_parts, sidm_EnergyInjectedCheck_all_parts;
  double SIDM_clight, SIDM_GroundStateMass;
#if defined(SIDM_CONST_CROSS) || defined(SIDM_MAXWELLIAN)
  double CrossSectionPerMass_in_cgs;
#endif
#endif

#ifdef DUST_LIVE
  int DesNumNgbDust;
  double MaxNumNgbDeviationDust;
#ifndef DL_DRAG_SEMI_IMPLICIT
  double StoppingTimeFrac;
#endif
#ifdef DL_GRAIN_BINS
  double MinGrainSize;
  double MaxGrainSize;
  double GrainDensity;
  double MaxBinFracMassChg;
#if defined(DL_SHATTERING) || defined(DL_COAGULATION) || defined(DL_PRODUCTION)
  char GrainDataPath[MAXLEN_PATH];
#endif
#ifdef DL_PRODUCTION
  double DustTargetFrac;
  int NumDustPerSpawn;
#ifdef DL_REFINEMENT
  double DustMaxFrac;
#endif
#endif
#endif
#endif

#ifdef ATOMIC_DM
  double ADMProtonMassInkeV, ADMElectronMassInkeV, ADMFineStructureConstant;
#endif

#ifdef ADDBACKGROUNDGRID
  int GridSize;
#endif

#ifdef DMPIC
  int PicCount;
  int PicCurrentNr;
  char PicList[1024];
  char PicDataDir[1024];
  char PicDataName[1024];
  int PixelsX, PixelsY;
#endif


#ifdef BLACK_HOLES
#ifdef BH_RELATIVE_NGB_DEVIATION
  double DesNumNgbBlackHoleRelDeviationFactor;
#endif
  double MaxNumNgbDeviationBlackHole;
#endif

#if defined(TREE_RAD) || defined(TREE_RAD_H2) || defined(TREE_RAD_CO)
  integertime TreeRadUpdateTime;
  int TreeRadUpdate;
#endif

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  double RayMemFac;
  double MachMin;
  int RayStepsMax;

#ifdef SKIP_BORDER
  double DistToBorder;
#endif

#ifdef SHOCK_SUBVOLUME
  int ShockSubvolumeNum;
#endif

#ifdef SHOCK_FINDER_POST_PROCESSING
  int NumFilesPerOutput;
  char OutputDirShockFinder[1024];
#endif
#endif


#ifdef DG
  double DG_beta;
#ifdef DISCONTINUITY_DETECTION
  double DG_alpha;
#endif
#ifdef MINMOD_B
  double DG_M;
#endif
#ifdef MACHNUM_B
  double DG_lim_mach_min;
#endif
#ifdef ANGLE_BOUND
double DG_min_angle;
#endif
#ifdef RK2
  double DG_RK2_alpha;
#endif
#if defined (REFINEMENT_SPLIT_CELLS) || defined (REFINEMENT_MERGE_CELLS)
double DG_TargetSlope;
double DG_SlopeRangeFactor;
#endif
#endif

#ifdef COSMIC_RAYS
  double GammaCR; /**< fiducial value: 4./3. */
  double CREnergyInputPerSolarMassOfStarFormation; /**< fiducial value: 4e47 = 4e48 [erg] for supernova energy return * 0.1 for efficiency */
  double DesNumNgbCRInjection;
  double MaxNumNgbDeviationCRInjection;
#endif

#ifdef COSMIC_RAYS_DIFFUSION
  integertime cr_diffusion_Ti_endstep, cr_diffusion_Ti_begstep;
#ifdef COSMIC_RAYS_DIFFUSION_CONSTANT_TIMESTEP
  double CRDiffusionTimestep;
#else
  double CourantCRDiffusion;
#endif
  double CR_Diffusion_Coefficient;
  double dt_cr_diffusion;

#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_X
  double CR_Boundary_X_Lower, CR_Boundary_X_Upper;
#endif
#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Y
  double CR_Boundary_Y_Lower, CR_Boundary_Y_Upper;
#endif
#ifdef COSMIC_RAYS_DIFFUSION_BOUNDARY_Z
  double CR_Boundary_Z_Lower, CR_Boundary_Z_Upper;
#endif

#endif

#ifdef COSMIC_RAYS_STREAMING_GLOBAL_CHI
  double CR_Chi;
#endif

#ifdef EXTERNALSHEARBOX
  double ShearBoxSigma0;        // surface density of disk in code units
  double ShearBoxSigmaNow;        // surface density of disk in code units
  double ShearBoxFg;            // gas fraction
  double ShearBoxMu;            //this will be unnecessary in future work
  double ShearBoxB;
#endif

#ifdef FM_SN_COOLING_RADIUS_BOOST
  double one_SNe_energy;	/* energy of one SN in code units */
  double one_SNe_mass;		/* mass of one SN in code units   */
  double SNe_velocity;		/* resulting SN velocity          */
#endif

#ifdef LOCAL_FEEDBACK
  double LocalFeedbackSNEnergy;
  double LocalFeedbackSNMassReturn;
  double LocalFeedbackSNRate;

#if defined(COSMIC_RAYS) && !defined(COSMIC_RAYS_SHOCK_ACCELERATION)
  double LocalFeedbackCRInjectionFraction;
#endif

#ifdef LOCAL_FEEDBACK_PARTICLES
  double LocalFeedbackSNTimeDelay;
  double LocalFeedbackSNTimeSpread;
#endif

#ifdef LOCAL_KINETIC
  double LocalFeedbackKineticInjectionFraction;
#endif

#if !defined(EXTERNALSHEARBOX_KSRATE) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  double LocalFeedbackSFEff;
  double LocalFeedbackSFDenThresh;
#endif

#endif

#ifdef NON_IDEAL_MHD
#ifdef OHMIC_DIFFUSION
  double OhmicDiffusionCoefficient;
#endif
#ifdef AMBIPOLAR_DIFFUSION
  double AmbipolarDiffusionCoefficient;
#endif
#ifdef NON_IDEAL_MHD_EXPLICIT_LIMIT_TIMESTEP
 double NonidealMHDTimelimit;
#endif

#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
  integertime ohmdiffusion_Ti_endstep, ohmdiffusion_Ti_begstep;
  double dt_ohmdiffusion, dt_max_ohmdiffusion;
  double OhmicDiffusionCoefficient;
#endif

#ifdef ONEDIMS_SPHERICAL
  double CoreMass;
  double CoreRadius;
#endif

#ifdef GRAVITY_TABLE 
  char ExternalGravForcesFile[MAXLEN_PATH]; 
  double DeltaR, DeltaZ, MaxR, MaxZ, MinR, MinZ; 
  int NGravityTableBins;
#endif

#ifdef ROTATING_HIGHRES_REGION
  double Highres_x0, Highres_y0;
  double Highres_delta, Highres_deltaz;
  double Highres_vrot;
  double Highres_t0;
  double Highres_targetmass;
#endif 

#ifdef MODGRAV
  char ModifiedGravityFile[200];
  int MaxAMRLevel;
  int MinLevelTopLeaf;
  int MaxLevelTopLeaf;
  int MaxLevelFullTree;
#endif
}
All;




/** This structure holds all the information that is
 * stored for each particle of the simulation.
 */
extern struct particle_data
{
  MyDouble Pos[3];              /**< particle position at its current time */
  MyDouble Mass;                /**< particle mass */
  MyFloat Vel[3];               /**< particle velocity at its current time */
  MySingle GravAccel[3];         /**< particle acceleration due to gravity */
#ifdef PMGRID
  MySingle GravPM[3];             /**< particle acceleration due to long-range PM gravity force */
#endif
#ifdef BECDM
  MyDouble PsiRe;
  MyDouble PsiIm;
#endif
#ifdef FORCETEST
  MyFloat GravAccelDirect[3];   /*!< particle acceleration calculated by direct summation */
  MyFloat PotentialDirect;
  MyFloat DistToID1;
#ifdef PMGRID
  MyFloat GravAccelShortRange[3];
  MyFloat GravAccelLongRange[3];
  MyFloat PotentialShortRange;
  MyFloat PotentialLongRange;
#endif
#endif

#if defined(EVALPOTENTIAL) || defined (OUTPUTPOTENTIAL)
  MySingle Potential;            /**< gravitational potential */
#if defined(PMGRID)
  MySingle PM_Potential;
#endif
#endif

#ifdef OUTPUTGRAVINTERACTIONS
  int GravInteractions;
#endif

#ifdef STELLARAGE
  MyFloat StellarAge;           /**< formation time of star particle */
#ifdef LOCAL_FEEDBACK_PARTICLES
  short int FeedbackDone;
#endif
#endif

#ifdef METALS
  MyFloat Metallicity;          /**< metallicity of gas or star particle */
#endif

#ifdef EXTERNALGRAVITY
  MyFloat ExtPotential;
#endif

#if defined(GFM) || defined(BLACK_HOLES) || defined(DUST_LIVE)
  MyIDType AuxDataID;
#endif

#ifdef TRACER_PARTICLE
  MyFloat TracerHsml;
#endif
#if defined(TRACER_PARTICLE) && defined(TRACER_PART_NUM_FLUID_QUANTITIES)
  MyFloat fluid_quantities[TRACER_PART_NUM_FLUID_QUANTITIES];
#endif

#if defined(TRACER_TRAJECTORY)
  MyFloat tRho, tTemp, tUtherm;
#ifdef TRACER_TRAJECTORY_EXTENDED_OUTPUT
  MyFloat tComposition[EOS_NSPECIES], tDedt;
#endif
#endif

#ifdef TRACER_MC
  int TracerHead;
  int NumberOfTracers;
  int OriginTask;
#endif

#if defined(OTVET_SCATTER_SOURCE) && !defined(GFM)
  MyFloat OtvetHsml, OtvetGasDensity;
#endif

#ifdef SIDM
  MyDouble sidm_Density[SIDM_STATES];
  MyDouble sidm_VelDisp[SIDM_STATES];
  MyDouble sidm_PSum[SIDM_REACTIONS];
  MyDouble sidm_NumNgb;
  MyDouble sidm_Hsml;
  int sidm_NumTotalScatter[SIDM_REACTIONS];
  int sidm_State;
#endif

#ifdef DMPIC
  MyFloat DM_VelDisp;
  MyFloat DM_Rho;
  MyFloat DM_Hsml;
#endif

#ifdef COSMIC_RAYS
  short int CRInjection;
  MyFloat Hsml;
#endif

#ifdef MODGRAV
  MySingle ModgravAccel[3];
#endif

  MyIDType ID;

#if defined(ADD_GROUP_PROPERTIES) || defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
  MyIDType FileOrder;
#endif
#ifdef ADD_GROUP_PROPERTIES
  MyIDType MinID;
  int MinIDTask;
  int OriginalIndex;
  int OriginalTask;
  int OriginalGrNr;
  int OriginalSubNr;
#endif

  integertime Ti_Current;       /**< current time on integer timeline */

  float OldAcc;         /**< magnitude of old gravitational force. Used in relative opening criterion */

  float GravCost[GRAVCOSTLEVELS];       /**< weight factors used for balancing the work-load */

  unsigned char Type;               /**< flags particle type.  0=gas, 1=halo, 2=disk, 3=bulge, 4=stars, 5=bndry */
  unsigned char SofteningType;
  signed char TimeBinGrav;
  signed char TimeBinHydro;
#ifdef SINKS
  signed char TimeBinSink;
#endif
}
 *P,                            /**< holds particle data on local processor */
 *DomainPartBuf;                /**< buffer for particle data used in domain decomposition */


#ifdef GRAVITY_TABLE 
/** holds the R and z components of the gravitational accelerations defined on a lookup table 
 */
extern struct grav_table_data
{
  double R, z;
  double  acc_R;
  double  acc_z;
}
 *GravT;
#endif




extern struct subfind_data
{
  int OriginIndex, OriginTask;
  int TargetIndex, TargetTask;
  int GrNr;

#ifdef SUBFIND
  int SubNr;
  int OldIndex;
  int submark;
  int originindex, origintask;
  MyFloat Utherm;
  MyFloat Density;
  MyFloat Potential;
  MyFloat Hsml;
  MyFloat BindingEnergy;

#ifdef CELL_CENTER_GRAVITY
  MyDouble Center[3];
#endif

#ifdef SUBFIND_CALC_MORE
  MyFloat SubfindHsml;
  MyFloat SubfindDensity;       /* total matter density */
  MyFloat SubfindDMDensity;       /* dark matter density */
  MyFloat SubfindVelDisp;       /* 3D DM velocity dispersion */
#endif
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
  int GroupNr;
#endif
#endif
}
 *PS;




/** Holds data that is stored for each hydro mesh cell in addition to
    the collisionless variables.
 */
extern struct sph_particle_data
{
  /* conserved variables */
  MyFloat Energy;
  MyFloat Momentum[3];
  MyFloat Volume;
  MyFloat OldMass;

#ifdef MRT_LSF_GRADIENTS
  MyFloat OldCons_DensPhot[MRT_BINS] ;
#ifdef MRT_COMOVING
  MyFloat Old_Vel[3] ;
#endif
#endif

#ifdef GENERAL_RELATIVITY
  /* store source terms */
  MyFloat srcEnergy;
  MyFloat srcMomentum[3];
#endif

  /* primitive variables */
  MyFloat Density;
  MyFloat Pressure;             /**< current pressure */
  MySingle Utherm;
#if defined(USE_ENTROPY_FOR_COLD_FLOWS) || defined(OUTPUT_ENTROPY)
  MyFloat Entropy;
  MyFloat A;
#ifdef ENTROPY_MACH_THRESHOLD
  MyFloat MaxMach;
#endif
#endif

#ifdef HIERARCHICAL_GRAVITY
  MySingle FullGravAccel[3];
#endif

#ifdef STICKYFLAGS
  short int StickyFlag;
#endif

#if !(defined(USE_ENTROPY_FOR_COLD_FLOWS) && defined(ENTROPY_MACH_THRESHOLD))
#if defined(TRACER_MC_MACHMAX) || (TRACER_PART_MACHMAX)
  MyFloat MaxMach;
#endif
#endif

#ifdef OUTPUT_MACHNUM
  MyFloat MaxMachNumber;
#endif
#ifdef MEASURE_DISSIPATION_RATE
  MyFloat DuDt;
#endif

  /* variables for mesh  */
  MyDouble Center[3];            /* center of mass of cell */
  MySingle VelVertex[3];         /* current vertex velocity (primitive variable) */
#ifdef REGULARIZE_MESH_LLOYD
  MyFloat Center_predict[3];
  MyFloat Center_anchor[3];
  MyFloat Volume_predict;
#endif
#ifdef REGULARIZE_MESH_SMOOTH
  MyFloat VelVertexAvg[3];
  MyFloat VelVertexAvgNorm;
#endif
  MySingle  MaxDelaunayRadius;
  MySingle Hsml;                 /* auxiliary search radius for points around a delaunay triangle */
  MySingle SurfaceArea;
#if defined(REGULARIZE_MESH_FACE_ANGLE) || defined(OUTPUT_MESH_FACE_ANGLE)
  MySingle MaxFaceAngle;
#endif
  MySingle ActiveArea;

#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
  /* fields that are potentially written into the AREPO snapshots */
  MySingle Machnumber;
  MySingle EnergyDissipation;

#ifdef SHOCK_FINDER_ON_THE_FLY
  SHOCK_FINDER_FIELDS
#ifdef COSMIC_RAYS
  SHOCK_FINDER_CR_FIELDS
#else
  SHOCK_FINDER_IDEAL_HYDRO_FIELDS
#endif
#endif

#ifdef SHOCK_FINDER_BEFORE_OUTPUT_MORE
  SHOCK_FINDER_FIELDS_MORE
#ifndef COSMIC_RAYS
  SHOCK_FINDER_IDEAL_HYDRO_FIELDS
#endif
#endif
#endif

#if defined(OUTPUT_DIVVEL) || defined(TGCHEM) || defined(SGCHEM)
  MyFloat DivVel;               /* divergence of the velocity field */
#endif
#if defined(REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED) || defined(OUTPUT_CURLVEL)
  MySingle CurlVel;              /* magnitude of the curl of the velocity field */
#endif
#if defined(VS_TURB) || defined(AB_TURB) || defined(ADJ_BOX_POWERSPEC)
  MyFloat Vorticity[3];         /* keep vorticity */
#endif
#ifdef OUTPUT_VERTEX_VELOCITY_DIVERGENCE
  MyFloat DivVelVertex;
#endif

#ifdef OUTPUT_CELL_SPIN
  MyFloat Spin[3];
  MyFloat CenterOld[3];
  MyFloat CenterOffsetMass[3];
  MyFloat CenterOffset[3];
#endif
#ifdef ACTIVE_CELL_SPIN
  MyFloat Omega[3];
  MyFloat MomentOfInertia[3][3];
#endif

#ifdef TREE_BASED_TIMESTEPS
  MySingle CurrentMaxTiStep;
  MySingle Csnd;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  MyDouble DuDt_diss;
  MyDouble DuDt_drive;
  MyDouble EgyDiss;
  MyDouble EgyDrive;
  MyDouble TurbAccel[3];
#endif

#ifdef VARIABLE_GAMMA
  MyFloat GammaE;
  MyFloat GammaC;
#endif

#if defined(REFINEMENT_HIGH_RES_GAS) && !defined(TGSET)
  MyFloat HighResMass;
  MyFloat HighResDensity;
#endif

#ifdef REFINEMENT_RPS
  MyFloat RPSGalaxyMass;
  MyFloat RPSGalaxyDensity;
#endif

#ifdef REFINEMENT_AROUND_BH
  char RefBHFlag;
  MyFloat RefBHMaxRad;
#endif


#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
  double EOSTemperature;
  double Composition[EOS_NSPECIES];     /* mass fractions of nuclei */
  double MassComposition[EOS_NSPECIES]; /* mass of a nuclei in the cell */
#endif

#ifdef NUCLEAR_NETWORK
  double dedt;
#endif

#ifdef SGCHEM
  double TracAbund[SGCHEM_NUM_ADVECTED_SPECIES];
  double MassTracAbund[SGCHEM_NUM_ADVECTED_SPECIES];
  double DustTemp;

#ifdef SGCHEM_VARIABLE_Z
  double CarbAbund;
  double OxyAbund;
  double MAbund;
  double ZAtom;
  /* Element abundance * cell mass - needed for advection */
  double CarbMass;
  double OxyMass;
  double MMass;
  double ZMass;

  double DustToGasRatio;
  /* Dust-to-gas ratio (in units of solar) times gas mass; to get true dust mass,
   * multiply by actual value of solar dust to gas ratio
   */
  double ScaledDustMass;
#endif

#ifdef SGCHEM_OUTPUT_COOLTIME
  double CoolTime;
#endif
#ifdef SGCHEM_DUMP_THERMAL_RATES
  double HeatCoolRates[SGCHEM_NUM_THERMAL_RATES];
#endif
#ifdef TREE_RAD
#define NPIX 12*NSIDE*NSIDE
  double Projection[NPIX];
#ifdef TREE_RAD_H2
  double ProjectionH2[NPIX];
#endif
#ifdef TREE_RAD_CO
  double ProjectionCO[NPIX];
#endif
#else                           /* TREE_RAD */
#define NPIX 1
#endif

#ifdef TREE_RAD_VEL
  MyFloat Vth2;   /* Square of the thermal velocity */
#endif

#endif                          /* SGCHEM */

#ifdef MHD
  MyFloat B[3];
  MyFloat BConserved[3];
  MyFloat DivB;
#ifdef MHD_THERMAL_ENERGY_SWITCH
  MyFloat Etherm;
#endif
#ifdef MHD_CT
  MyFloat AConserved[3];
  MyFloat A[3];
  double TimeLastBUpdate;
#endif
  MyFloat CurlB[3];
#endif

#ifdef METALS
  MyFloat Metallicity;
  MyFloat MassMetallicity;
#endif

#ifdef CONDUCTION_SATURATION    //VITALI
  MyFloat GradEntr[3];
#endif

#ifdef CONDUCTION_ANISOTROPIC
  double gradU[3];
#endif

#ifdef TURBULENT_METALDIFFUSION
#ifndef GFM_STELLAR_EVOLUTION
  MyFloat Metallicity;
  MyFloat MassMetallicity;
  MyFloat MetalsFraction[GFM_N_CHEM_ELEMENTS];
  MyFloat MassMetals[GFM_N_CHEM_ELEMENTS];
#endif
  double gradm[3];
  double st[3][3];
#endif

#ifdef GFM_LAMBDA
  MyFloat Lambda[6];
#endif

#ifdef GFM_STELLAR_EVOLUTION
  MyFloat Metallicity;
  MyFloat MassMetallicity;
  MyFloat MetalsFraction[GFM_N_CHEM_ELEMENTS];
  MyFloat MassMetals[GFM_N_CHEM_ELEMENTS];
#ifdef GFM_DUST
#ifdef GFM_DUST_CAP
  MyFloat DustMassCap;
#endif
  MyFloat MetalsDustFraction[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
  MyFloat MassMetalsDust[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
  MyFloat DustTauGrowth;
#ifdef GFM_DUST_SPUTTERING
#if (GFM_DUST_SPUTTERING==1)
  MyFloat DustTauSputter;
#endif
#endif
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
  MyFloat NumSNII;
#endif
#if defined(DUST_LIVE) && defined(DL_SNE_DESTRUCTION)
  MyFloat SNRate; /* local supernova rate in Gyr^-1 */
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MassMetalsChemTags[GFM_N_CHEM_TAGS];
  MyFloat MassMetalsChemTagsFraction[GFM_N_CHEM_TAGS];
#endif
#endif

#ifdef RT_SHORT_CHARACTERISTICS
  MyFloat ColumnDensity0;
  MyFloat ColumnDensity1;
  MyFloat v_shell;
  MyFloat SourceDistance;
  MyFloat nH;
  MyFloat nHI;
  MyFloat nHII;
  MyFloat n_elec;
#endif

#if defined(RT_ADVECT) || defined(MRT)
#ifndef MRT
  struct rt_grad_data rt_Grad;
#else
  struct rt_grad_data RTGrad ;
#endif
  MyFloat nHI;
  MyFloat nHII;
  MyFloat ne ;
#if defined(RT_INCLUDE_HE) || defined(MRT_INCLUDE_HE)
  MyFloat nHeI;
  MyFloat nHeII;
  MyFloat nHeIII;
#endif
  MyFloat n_elec;
  MyFloat Photons[RT_N_DIR];

#ifdef RT_ADVECT
  MyFloat DensPhot[RT_N_DIR];
#else
  MyFloat DensPhot[MRT_BINS];
  MyFloat Cons_DensPhot[MRT_BINS];
  MyFloat Cons_DensPhot_absorbed[MRT_BINS];
#endif

  MyFloat RT_F[MRT_BINS][3] ;
  MyFloat Cons_RT_F[MRT_BINS][3] ;
  MyFloat FN[MRT_BINS][3] ;
  MyFloat modFN[MRT_BINS] ;

#ifdef MRT_IR
  MyFloat KappaIR_R[IR_BINS] ; /*Rosseland mean opacity*/
  MyFloat KappaIR_P[IR_BINS] ; /*Planck mean opacity*/
#ifdef MRT_IR_PHOTON_TRAPPING
  MyFloat Trapped_Cons_DensPhot[IR_BINS] ;
  MyFloat Trapped_DensPhot[IR_BINS] ;
#endif
#endif

#ifndef RT_HEALPIX_NSIDE
  int SourceID[RT_N_DIR];
  MyFloat SourcePos[RT_N_DIR][3];
#endif
#endif


#if defined(GRACKLE) && !defined(GRACKLE_TAB)
  MyFloat e_frac; /*fraction of mass in electrons*/
  MyFloat e_mass;    /*mass of electrons*/
  MyFloat GrackleSpeciesFraction[GRACKLE_SPECIES_NUMBER];
  MyFloat GrackleSpeciesMass[GRACKLE_SPECIES_NUMBER];
#endif

#ifdef COFFEE_PROBLEM
  int Flag;
#endif

#ifdef TRACER_FIELD
  MyFloat Tracer;
  MyFloat ConservedTracer;
#endif

#ifdef PASSIVE_SCALARS
  MyFloat PScalars[PASSIVE_SCALARS];
  MyFloat PConservedScalars[PASSIVE_SCALARS];
#endif

#ifdef OUTPUT_SURFACE_AREA
  int CountFaces;
#endif
#if defined(REFINEMENT_SPLIT_CELLS)
  MySingle MinimumEdgeDistance;
#endif

#if defined(COOLING) && !defined(GRACKLE)
  MyFloat Ne;                   /* electron fraction, expressed as local electron number
                                   density normalized to the hydrogen number density. Gives
                                   indirectly ionization state and mean molecular weight. */
#endif

#ifndef COOLING
#ifdef ATOMIC_DM
  MyFloat Ne;
#endif
#endif

#ifdef USE_SFR
  MySingle Sfr;
#ifdef DELAYED_COOLING
  MyFloat CoolShutoffTime;
#endif
#ifdef DELAYED_COOLING_TURB
  MyFloat Uturb;
  MyFloat MassUturb;
#endif
#endif

#if defined(FM_RADIATION_FEEDBACK) || defined(ISM_HII_HEATING)
  MyFloat GasRadCoolShutoffTime;
#endif

#if defined(GFM_WINDS_VARIABLE) || defined(GFM_WINDS_LOCAL) || (defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF==0))
  union
  {
    MySingle HostHaloMass;
    MySingle DMVelDisp;
  } w;
#endif

#ifdef AMR
  int Level;

#ifdef OUTPUT_AMR_REFFLAG
  int refflag;
#endif
#endif

#ifdef GFM_AGN_RADIATION
  MySingle AGNBolIntensity;
#endif

#ifdef GFM_WINDS_LOCAL
  MyFloat WindEnergyReceived;
#endif

#ifdef GFM_BIPOLAR_WINDS
#if (GFM_BIPOLAR_WINDS == 3)
  MyFloat DensGasAngMomentum[3];
#else
  MyFloat GroupVel[3];
  MyFloat GroupGravAcc[3];
#endif
#endif

#ifdef OUTPUT_COOLHEAT
  MyFloat CoolHeat;
#endif

#ifdef BH_THERMALFEEDBACK
  MySingle Injected_BH_Energy;
#endif
#ifdef BH_ADIOS_WIND
  MyFloat Injected_BH_Wind_Momentum[3];
#endif

#ifdef TGCHEM
  double Abund[TGCHEM_NUM_ABUNDANCES];
  double Gamma;
  double HydroHeatRate;
#endif

#ifdef HEALRAY
  double HeatRate;
#endif

#if defined(TGCHEM) || defined(HEALRAY)
  double EscFrac;
#endif

  struct grad_data Grad;

#ifdef RUNGE_KUTTA_FULL_UPDATE
  struct conservative_variables rk;
#endif

#ifdef TVD_SLOPE_LIMITER
  struct grad_data GradUl;
#endif

#ifdef SECOND_DERIVATIVES
  struct hessian_data Hessian;
#endif

#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
  int first_connection;
  int last_connection;
#endif

#ifdef DEREFINE_GENTLY
  int DoNotDerefFlag;
#endif

#ifdef REFINEMENT_HIGH_RES_GAS
  int AllowRefinement;
#endif

#ifdef REFINEMENT_SPLIT_CELLS
  MySingle SepVector[3];
#endif

#ifdef REFINEMENT_MERGE_PAIRS
  MyIDType DerefPartnerId;      /* Id of partner cell we want to merge with */
  int DerefPartnerIndex;        /* local index of partner cell. This value is not always valid! */
#endif

#ifdef SPECIAL_BOUNDARY
  MyFloat MinDistBoundaryCell;
#endif

#ifdef REFINEMENT_VOLUME_LIMIT
  MyFloat MinNgbVolume;
#endif

#ifdef SINKS
  short int InAccrRadius;
#endif

#ifdef OUTPUT_MOLECULAR_FRACTION
  MyFloat MolecularFrac;
#endif

#ifdef OUTPUT_OPTICAL_DEPTH
  MyFloat OpticalDepth;
#endif

#ifdef FLD
  MyFloat Temperature;

  MyFloat n_gamma;
  MyFloat Lambda;
  MyFloat Kappa_P;
  MyFloat Kappa_R;
  MyFloat Kappa_diff;
  MyFloat R2;

  double b;

#ifdef FLD_CONES
  MyFloat gammas[FLD_NCONES];
#endif

#ifndef FLD_HYPRE_IJ1

#ifndef FLD_CONES
  double w[2*NUMDIMS+1];
#else
  double w[9];
#endif

#else
  MyFloat w;
#endif

#ifdef FLD_MG
  MyFloat residuum;
#endif

#endif

#ifdef OTVET
  MyFloat ET[6];                /* eddington tensor - symmetric -> only 6 elements needed */
  MyFloat Je[OT_N_BINS];        /* emmisivity */
  MyFloat n_gamma[OT_N_BINS];
  MyFloat HI;
  MyFloat HII;
#if !defined(COOLING)
  MyFloat Ne;
#endif
#ifdef RADPRESS_OPT_THICK
  MyFloat vector_n[3];          /* computes the direction normal to stellar sources, uses loop of eddington tensor calculation */
  MyFloat RadPress[3];          /* notice same name than optically thin case, but both can't be activated at once */
  MyFloat n_gamma_abs[OT_N_BINS];       /* number density of photons absorved in a given timestep */
#endif

#ifdef OTVET_INCLUDE_HE
  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;
#endif

#ifdef OTVET_FLUXLIMITER
  MyFloat Grad_ngamma[3][OT_N_BINS];
#endif
#endif                          /* OTVET */

#ifdef RADPRESS_OPT_THIN
  MyFloat RadPress[3];
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  MyFloat TotEgyFeed;                           /**< cumulative energy received by stellar feedback */
  MyFloat IntEgyFeed;                           /**< cumulative internal energy received by stellar feedback */
  MyFloat KinEgyFeed;                           /**< cumulative kinetic energy received by stellar feedback */
#endif

#ifdef DVR_RENDER
  char DoMesh;
  MyFloat DvrFields[DVR_NUM_FIELDS];
#endif

#ifdef COSMIC_RAYS
  MyFloat CR_Energy;              /**< Cosmic ray energy in a cell */
  MyFloat CR_SpecificEnergy;      /**< Cosmic ray energy per unit gas mass in a cell */
  MyFloat CR_Pressure;            /**< Cosmic ray pressure in a cell */
#ifdef COSMIC_RAYS_STREAMING
  MyFloat CR_Chi;
#endif
#endif

/*
#ifdef FM_STAR_FEEDBACK
#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB)
  short IsNotKicked;
#endif
#endif
*/


#ifdef DG
  MyDouble Weights[NOF_BASE_FUNCTIONS][5]; /** < Stores the degrees of freedom associated with the cell */
#ifdef DISCONTINUITY_DETECTION
  int Inflow_boundaries;
#ifdef OUTPUT_DG_DISCONTINUITIES
  double Discontinuity;
#endif
#endif
#ifdef MACHNUM_JUMP_DETECTION
#ifdef OUTPUT_DG_MACHNUMS
  double Machnum;
#endif
#endif
#ifdef OUTPUT_MIN_ANGLES
  double min_angles_x[5];
  double min_angles_y[5];
  double min_angles_z[5];
#endif
#ifdef DG_TURBULENCE
  MyDouble WeightsDiss[NOF_BASE_FUNCTIONS];
  MyDouble WeightsDeDt[NOF_BASE_FUNCTIONS];
#endif
#endif




#ifdef ISM
#ifdef ISM_LOCAL_RADIATION_PRESSURE
  MyDouble ClumpPos[3];
  MyFloat ClumpDensity;
  char ClumpCenterFound;
  MyFloat ClumpRadius;
  MyIDType ClumpID;
  MyFloat ClumpMass;
  MyFloat ClumpLuminosity;
#endif
#endif

  double TimeLastPrimUpdate;

#ifdef RADCOOL
  double Phios;                 /*Normalization or amount of flux from new stars weighted by the inverse distance squared */
  double Phins;                 /*Normalization of flux from old stars weighted by inverse disctance squared */
#ifdef RADCOOL_HOTHALO
  double PhiT6;                 /*Normalization or amount of flux from gas with 5.5<logT<6.5  weighted by the inverse distance squared */
  double PhiT7;                 /*Normalization or amount of flux from gas with 6.5<logT<7.5  weighted by the inverse distance squared */
  double PhiT8;                 /*Normalization or amount of flux from gas with 7.5<logT<8.5  weighted by the inverse distance squared */
  MyFloat Temperature;
#endif
#endif

#ifdef ADDBACKGROUNDGRID
  MyFloat Weight;
#endif


#ifdef MRT
#ifdef MRT_LSF_GRADIENTS
  MyFloat PT[MRT_BINS][DIMS][DIMS] ;
#endif
  MyFloat HI;
  MyFloat HII;
#if !defined(COOLING)
  MyFloat Ne;
#endif
#ifdef MRT_INCLUDE_HE
  MyFloat HeI;
  MyFloat HeII;
  MyFloat HeIII;
#endif
#endif                          /* MONOTONE_RT */

#if defined(DUST_LIVE) && defined(DL_GRAIN_BINS) && (defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION))
  MyFloat CloudFrac;
#endif

#ifdef CALCULATE_QUANTITIES_IN_POSTPROCESS
  MyFloat Vorticity[DIMS] ;
#endif
}
 *SphP,                         /**< holds SPH particle data on local processor */
 *DomainSphBuf;                 /**< buffer for SPH particle data in domain decomposition */

#ifdef TRACER_MC
extern struct tracer_linked_list
{
  MyIDType ID;                  /**< unique ID of the tracer particles */
#ifdef TRACER_MC_CHECKS
  MyIDType ParentID;
#endif
#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  MyFloat fluid_quantities[TRACER_MC_NUM_FLUID_QUANTITIES];  /**< recorded properties of parent gas cell */
#endif
  int Next;                     /**< list index of next tracer, or -1 if end */
  int Prev;                     /**< list index of prev tracer, or (-index-1) of SphP if first child */
}
 *TracerLinkedList;

extern int *TracerLinkedListHeap;

extern MyIDType *tracer_cellids;
#endif /* TRACER_MC */

#ifdef GFM
extern struct star_particle_data
{
  unsigned int PID;
  MyFloat BirthTime;
  MyFloat BirthPos[3];
  MyFloat BirthVel[3];
  MyFloat BirthDensity;
#ifdef GFM_STELLAR_EVOLUTION
  MyDouble InitialMass;
  MyFloat MassMetals[GFM_N_CHEM_ELEMENTS];
  MyFloat Metallicity;
  MyFloat SNIaRate;
  MyFloat SNIIRate;
#endif

#ifdef GFM_DISCRETE_ENRICHMENT
  // To accumulate mass/metals for discrete enrichment events, need to keep track of last enrichment time
  MyDouble lastEnrichTime;
#endif

#ifdef GFM_DUST
  /* To ensure wind particles return the gas-phase metal and dust content of */
  /* the ISM at their creation when recoupled (or stripped), store the */
  /* relative ratios of gas-phase and dust channel metals for each species. */
  /* For each species, the sum of the gas metal fraction plus all dust */
  /* channel fractions should sum to 1. */
  MyFloat InitialMetalFractions[GFM_N_CHEM_ELEMENTS];
  MyFloat InitialDustFractions[GFM_DUST_N_CHANNELS][GFM_N_CHEM_ELEMENTS];
#endif
#ifdef GFM_CHEMTAGS
  MyFloat MassMetalsChemTags[GFM_N_CHEM_TAGS];
#endif
#if defined(GFM_STELLAR_EVOLUTION) || defined(GFM_WINDS) || defined(MRT_SOURCES)
  MyFloat Hsml;
#endif
#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)
  MyFloat Utherm;
#endif
#if defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF==0)
  MySingle DMVelDisp;
#endif
#ifdef FM_STAR_FEEDBACK
#ifdef OUTPUT_STELLAR_FEEDBACK
  MyFloat SNII_Num;
  MyFloat FeedbackEnergy;
  MyFloat FeedbackKinEnergy;
  MyFloat FeedbackThEnergy;
  MyFloat FeedbackMomentum;
  MyFloat TotalMassReleased;
  MyFloat TotalMassToEnrich;
  MyDouble Cum_SNII_Num;
  MyDouble Cum_FeedbackMomentum;
#endif
#ifdef DELAYED_COOLING
  MyFloat BlastRadius;
#endif
#ifdef INSTANTANEOUS_DEPOSITION
  short int FeedbackFlag;
#endif
#endif
#if defined(FM_RADIATION_FEEDBACK) || defined(ISM_HII_HEATING)
  MyFloat StromgrenRadius;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat StromgrenMass;
#endif
  MyFloat RadFeedTau;
  int RadFeed_NumNgb;
  int RadFeed_Flag;
  MyFloat RadFeed_NormSphRad;

#ifdef FM_RADIATION_FEEDBACK_DEBUG      /* #LVS: TODO --> once tested, keep variables in "debug" for only StarParticle */
  MyDouble RadiationMomentumReleased;
  MyFloat RadVelocityKick;
  MyFloat NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat NormSphRadFeedback_cold;
#endif
  MyFloat RadCoolShutoffTime;
#endif
#endif
#ifdef FM_EARLY_STAR_FEEDBACK
  MyDouble EarlyCumulativeFeedback;     /* stores the total(cumulative) amount of energy input as early feedback throughout lifetime of star */
#endif
#ifdef OUTPUT_EARLY_STELLAR_FEEDBACK
  MyFloat EarlyFeedbackMomentum;
  MyDouble Cum_EarlyFeedbackMomentum;
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
  MyFloat HighResMass;
#endif
#ifdef ISM_LOCAL_RADIATION_PRESSURE
  MyFloat Luminosity;
#endif

#if defined(FM_SN_COOLING_RADIUS_BOOST) || defined(FM_RADIATION_FEEDBACK)
MyDouble LocISMdens;			/* stores the local ISM density (and metal density) for the cooling radius calculation */
MyDouble LocISMZdens;
#endif

#ifdef FM_VAR_SN_EFF
MyDouble AvgMetalNgb;
#endif
#ifdef FM_MASS_WEIGHT_SN
MyDouble TotNgbMass;
#endif

#if defined(DUST_LIVE) && defined(DL_PRODUCTION)
  MyDouble DeltaDustMassTot;
  MyDouble DeltaDustMass[GFM_N_CHEM_ELEMENTS];
  enum gsd_dnda_type DndaType;
#endif
}
 *StarP, *DomainStarBuf;
#endif

#ifdef BLACK_HOLES
extern struct bh_particle_data
{
  unsigned int PID;
  int BH_CountProgs;
  MyFloat BH_NumNgb;
  MyFloat BH_Hsml;
  MyFloat BH_Mass;
  MyFloat BH_Mdot;
  MyFloat BH_MdotBondi;
  MyFloat BH_MdotEddington;
  MyFloat BH_CumMass_QM;
  MyFloat BH_CumEgy_QM;
  MyFloat BH_CumMass_RM;
  MyFloat BH_CumEgy_RM;
  MyFloat BH_DtGasNeighbor;
  MyFloat BH_VolSum;
  MyFloat BH_Density;
  MyFloat BH_U;
  MyFloat BH_Pressure;
  MyFloat BH_SurroundingGasVel[3];
  MyIDType SwallowID;
#ifdef BH_NEW_CENTERING
  MyFloat HsmlCentering;
#endif
#ifdef MEASURE_POTMIN_AROUND_BH
  MyDouble BH_MinPotPos[3];
  MyFloat BH_MinPot;
  MyFloat BH_MinPot_ActiveM;
  MyFloat BH_MinPot_TotalM;
#endif
#ifdef BH_FRICTION
  MyDouble BH_MinPotPos_Previous[3];
  integertime BH_MinPotTime;
  integertime BH_MinPotTime_Previous;
  MyFloat BH_MinPotCumAvgTime;
  MyFloat BH_RhoTot;
  MyFloat BH_MinPotVel[3];
  MyDouble BH_MinPotPos_Extended[3];
  MyFloat BH_MinPot_Extended;
#endif
#ifdef DRAINGAS
  MyFloat NearestDist;
  MyIDType DrainID;
  MyFloat CellDensity;
  MyFloat CellUtherm;
  MyDouble DrainBucketMass;
#endif
#ifdef BH_DYN_FRICTION
  MyFloat BH_DynFric_Hsml;
  MyFloat BH_DynFric_NumNgb;
#endif
#ifdef BH_BONDI_DISK_VORTICITY
  MyFloat BH_GasVort[3];
  MyFloat Gal_Mass;
#endif
#ifdef BH_BONDI_CAPTURE
  MyFloat BH_CaptureMass;
#endif
#ifdef BH_BUBBLES
  MyFloat BH_Mass_bubbles;
  MyFloat BH_Mass_ini;
#endif
#ifdef BH_ADIOS_WIND
  MyFloat BH_WindEnergy;
#ifndef BH_ADIOS_RANDOMIZED
  MyFloat Asum;
  MyFloat Bsum;
  MyFloat Msum;
  MyFloat Qsum[3];
#endif
#if defined(BH_ADIOS_WIND_DIRECTIONAL) || defined(BH_ADIOS_RANDOMIZED)
  MyFloat WindDir[3];
#endif

#endif
#if defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER)
  MyFloat HostHaloMass;
#endif
#if (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  MyFloat BH_DMVelDisp;
#endif
#ifdef BH_THERMALFEEDBACK
  MyFloat BH_ThermEnergy;
#endif
#ifdef BH_THERMALFEEDBACK_ACC
  MyFloat BH_AccEnergy;
  MyFloat BH_AccTime;
#endif

#ifdef BH_BIPOLAR_FEEDBACK
  MyFloat BH_BipolarJ[3];
  MyFloat BH_BipolarSum;
  MyFloat BH_BipolarColdFraction;
  MyFloat BH_BipolarColdMass;
  int BH_BipolarColdDisk;
#endif

#ifdef BH_NF_RADIO
  MyFloat BH_Mdot_quasar;
  MyFloat BH_Mdot_radio;
  MyFloat BH_RadioEgyFeedback;
  MyFloat BH_HaloVvir;
  MyFloat BH_XrayLum;
  MyFloat BH_RadioLum;
  MyIDType ID_Min_BH_Potential;
#endif
#ifdef BH_USE_ALFVEN_SPEED_IN_BONDI
  MyFloat BH_Bpress;
#endif

#ifdef BH_SPIN_EVOLUTION
  MyFloat BH_SpinParameter;
  MyFloat BH_SpinOrientation[3];

  int BH_SpinModel;
  int BH_FlagOngAccEpis;
  MyFloat BH_AngMomGasCells[3];
  MyFloat BH_TimeAccretion_Previous;
  MyFloat BH_DTimeAccretion_Current;
  MyFloat BH_Mass_Previous;
  MyFloat BH_DMass_Current;
  MyFloat BlackHoleRadiativeEfficiency;


#endif
}
 *BHP, *DomainBHBuf;
#endif

#ifdef BLACK_HOLES
#define BPP(i) BHP[P[(i)].AuxDataID]
#endif

#ifdef GFM
#define STP(i) StarP[P[(i)].AuxDataID]
#endif

#ifdef DUST_LIVE
extern struct dust_particle_data
{
  unsigned int PID;
  MyFloat Hsml;
  MyFloat LocalGasDensity;
  MyFloat LocalGasVelocity[3];
  MyFloat LocalSoundSpeed;
  MyFloat DragAccel[3];
  MyFloat StoppingTime;
  int MinGasTimeBin;
#ifdef DL_GRAIN_BINS
  MyFloat NumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  MyFloat BinSlopes[DL_GRAIN_BINS];
#endif
  MyFloat MetalFractions[GFM_N_CHEM_ELEMENTS];
  MyFloat BinMassChgTau;
  MyFloat OrigMass;
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  MyFloat DustHsml;
  MyFloat DustNumNgb;
  MyFloat DustDensity;
  MyFloat LocalCloudFrac;
#endif
#endif
} *DustP, *DomainDustP;

#define DTP(i) DustP[P[(i)].AuxDataID]
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
extern struct special_particle_data
{
  MyIDType ID;
  double pos[3];
  double mass;
}
 *PartSpecialListGlobal;
#endif

#ifdef REFINEMENT_AROUND_DM
extern struct refine_dm_data
{
  MyIDType ID;
  double pos[3];
  double softening;
}
 *DMPartListGlobal;
#endif

#ifdef TRACER_MC
/* tracer communication structure for export/import */
extern struct tracer_flux_list_data
{
  int task, index;
  MyIDType ID;
#ifdef TRACER_MC_NUM_FLUID_QUANTITIES
  MyFloat fluid_quantities[TRACER_MC_NUM_FLUID_QUANTITIES];
#endif
#ifdef TRACER_MC_CHECKS
  MyIDType ParentID;
#endif
}
 *TracerFluxListIn, *TracerFluxListGet;

extern int Ntracerflux;                    /**< number of tracers in local export buffer */
extern int TracerFluxListGet_index; /**< walk through the TracerFluxListGet after exchange */
#endif

extern peanokey *DomainKeyBuf;

/* global state of system
*/
extern struct state_of_system
{
  double Mass,
    EnergyKin,
    EnergyPot,
    EnergyInt,
    EnergyTot,
    Momentum[4],
    AngMomentum[4],
    CenterOfMass[4],
    MassComp[NTYPES], EnergyKinComp[NTYPES], EnergyPotComp[NTYPES], EnergyIntComp[NTYPES], EnergyTotComp[NTYPES], MomentumComp[NTYPES][4], AngMomentumComp[NTYPES][4], CenterOfMassComp[NTYPES][4];
}
SysState, SysStateAtStart, SysStateAtEnd;



/** Struct used for passing the parameters during the mesh cell search. */
typedef struct
{
  MyDouble Pos[3];
  int Task;
  union
  {
    int Index;
    float hsmlguess;
  } u;

#ifdef TRACER_TRAJECTORY
  MyDouble Density;
#endif
#if defined(VORONOI_FIELD_DUMP_PIXELS_X) && defined(VORONOI_FIELD_DUMP_PIXELS_Y)
  int pixel_index;
#endif
} mesh_search_data;

/** Struct used for sending positions to other tasks during the
    mesh cell search. */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Distance;
} mesh_search_request;

/** Struct used for receiving the results from other tasks during the
    mesh cell search. */
typedef struct
{
  MyDouble Distance;
  int Task;
  int Index;
#ifdef TRACER_TRAJECTORY
  MyDouble Density;
#endif
} mesh_search_response;



extern struct data_partlist
{
  int Task;          /** The task the item was exported to. */
  int Index;         /** The particle index of the item on the sending task. */
}
 *PartList;


extern struct datanodelist
{
  int Task;   /** target process */
  int Index;  /** local index that wants to open this node */
  int Node;   /** node to be opened on foreign process */
}
 *NodeList;


#define FAC_AVG_NODES_PER_EXPORT 4.0    /* default choice for estimated average number of exported nodes per exported particle */

extern struct directdata
{
  MyDouble Pos[3];
  MyDouble Mass;
  unsigned char Type;
  unsigned char SofteningType;
}
 *DirectDataIn, *DirectDataAll;

extern struct accdata
{
  MyFloat Acc[3];
#ifdef EVALPOTENTIAL
  MyFloat Potential;
#endif
}
 *DirectAccOut, *DirectAccIn;





#if defined (EVALPOTENTIAL) || defined (OUTPUTPOTENTIAL) || defined(SUBFIND)
extern struct potdata_out
{
  MyFloat Potential;
}
/** Holds the partial results computed for imported particles. Note:
    We use GravDataResult = GravDataGet, such that the result replaces
    the imported data */
 *PotDataResult,
/** Holds partial results received from other processors. This will
    overwrite the GravDataIn array */
 *PotDataOut;
#endif

/** Buffer of size NTask used for flagging whether a particle needs to
  be exported to the other tasks. */
extern int *Exportflag;
/** Buffer of size NTask used for counting how many nodes are to be
    exported to the other tasks? */
extern int *Exportnodecount;
/** Buffer of size NTask used for holding the index into the
    DataIndexTable. */
extern int *Exportindex;

/** Array of NTask size of the offset into the send array where the
    objects to be sent to the specified task starts. */
extern int *Send_offset,
/** Array of NTask size of the number of objects to send to the
    tasks. */
 *Send_count,
/** Array of NTask size of the number of objects to receive from the
    tasks. */
 *Recv_count,
/** Array of NTask size of the offset into the receive array where the
    objects from the specified task starts. */
 *Recv_offset;


extern int *TasksThatSend, *TasksThatRecv, NSendTasks, NRecvTasks;

extern struct send_recv_counts
{
  int Count;
  int CountNodes;
}
*Send, *Recv;

extern int *Send_offset_nodes, *Send_count_nodes, *Recv_count_nodes, *Recv_offset_nodes;


extern int Mesh_nimport, Mesh_nexport, *Mesh_Send_offset, *Mesh_Send_count, *Mesh_Recv_count, *Mesh_Recv_offset;

extern int Force_nimport, Force_nexport, *Force_Send_offset, *Force_Send_count, *Force_Recv_count, *Force_Recv_offset;

/** Header for the standard file format.
 */
#if (NTYPES==7 || NTYPES==8)
#define NTYPES_INT_HEADER 8
#else
#define NTYPES_INT_HEADER NTYPES
#endif
extern struct io_header
{
  int npart[NTYPES_INT_HEADER];                 /**< number of particles of each type in this file */
  double mass[NTYPES];          /**< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                  /**< time of snapshot file */
  double redshift;              /**< redshift of snapshot file */
  int flag_sfr;                 /**< flags whether the simulation was including star formation */
  int flag_feedback;            /**< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[NTYPES_INT_HEADER];   /**< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;             /**< flags whether cooling was included  */
  int num_files;                /**< number of files in multi-file snapshot */
  double BoxSize;               /**< box-size of simulation in case periodic boundaries were used */
  double Omega0;                /**< matter density in units of critical density */
  double OmegaLambda;           /**< cosmological constant parameter */
  double HubbleParam;           /**< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;          /**< flags whether the file contains formation times of star particles */
  int flag_metals;              /**< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[NTYPES_INT_HEADER];   /**< High word of the total number of particles of each type */
  int flag_entropy_instead_u;   /**< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;     /**< flags that snapshot contains double-precision instead of single precision */

  int flag_lpt_ics;             /**< flag to signal that IC file contains 2lpt initial conditions */
  float lpt_scalingfactor;      /**< scaling factor for 2lpt initial conditions */

  int flag_tracer_field;        /**< flags presence of a tracer field */

  int composition_vector_length;        /**< specifies the length of the composition vector (0 if not present)  */


#if (NTYPES==6)
  char fill[40];                /**< fills to 256 Bytes */
#elif (NTYPES==7)
  char fill[8];                 /**< fills to 256 Bytes */
#endif
}
header;                         /**< holds header for snapshot files */

/** Header for the ICs file format, if NTYPES does not match
 */
#ifdef NTYPES_ICS
extern struct io_header_ICs
{
  int npart[NTYPES_ICS];        /**< number of particles of each type in this file */
  double mass[NTYPES_ICS];      /**< mass of particles of each type. If 0, then the masses are explicitly
				   stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                  /**< time of snapshot file */
  double redshift;              /**< redshift of snapshot file */
  int flag_sfr;                 /**< flags whether the simulation was including star formation */
  int flag_feedback;            /**< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[NTYPES_ICS];  /**< total number of particles of each type in this snapshot. This can be
				   different from npart if one is dealing with a multi-file snapshot. */
  int flag_cooling;             /**< flags whether cooling was included  */
  int num_files;                /**< number of files in multi-file snapshot */
  double BoxSize;               /**< box-size of simulation in case periodic boundaries were used */
  double Omega0;                /**< matter density in units of critical density */
  double OmegaLambda;           /**< cosmological constant parameter */
  double HubbleParam;           /**< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;          /**< flags whether the file contains formation times of star particles */
  int flag_metals;              /**< flags whether the file contains metallicity values for gas and star
				   particles */
  unsigned int npartTotalHighWord[NTYPES_ICS];  /**< High word of the total number of particles of each type */
  int flag_entropy_instead_u;   /**< flags that IC-file contains entropy instead of u */
  int flag_doubleprecision;     /**< flags that snapshot contains double-precision instead of single precision */

  int flag_lpt_ics;             /**< flag to signal that IC file contains 2lpt initial conditions */
  float lpt_scalingfactor;      /**< scaling factor for 2lpt initial conditions */

  int flag_tracer_field;        /**< flags presence of a tracer field */

  int composition_vector_length;        /**< specifies the length of the composition vector (0 if not present)  */


#if (NTYPES_ICS==6)
  char fill[40];                /**< fills to 256 Bytes */
#else
    terminate("NTYPES_ICS != 6")
#endif
}
header_ICs;                             /**< holds header for IC files */
#endif


enum iofields
{ IO_POS,
  IO_VEL,
  IO_ID,
  IO_MASS,
  IO_TRACER,
  IO_SECONDORDERMASS,
  IO_U,
  IO_RHO,
  IO_VORT,
  IO_VOL,
  IO_CM,
  IO_VERTEXVEL,
  IO_FACEANGLE,
  IO_DIVVERTEXVEL,
  IO_SAREA,
  IO_NFACES,
  IO_SPIN,

  IO_TRACER_PARTICLE,
  IO_TRACER_MC_NumTracers,
  IO_TRACER_MC_ID,
  IO_TRACER_MC_ParentID,
  IO_TRACER_MC_FluidQuantities,

  IO_HIGHRESMASS,
  IO_COMPOSITION,
  IO_PRESSURE,
  IO_ENTROPY,
  IO_CSND,
  IO_NE,
  IO_NH,
  IO_SFR,
  IO_AGE,
  IO_LOCFBEVENT,
  IO_STKY,
  IO_Z,
  IO_BHMASS,
  IO_BHMDOT,
  IO_BHMDOTBONDI,
  IO_BHMDOTEDDIN,
  IO_BHHSML,
  IO_BHU,
  IO_BHRHO,
  IO_BHPRESS,
  IO_BHMBUB,
  IO_BHMINI,
  IO_BHPROGS,
  IO_BHCMQM,
  IO_BHCEQM,
  IO_BHRADIOLUM,
  IO_BHXRAYLUM,
  IO_BHCMRM,
  IO_BHCERM,

  IO_BHMDOTQUASAR,
  IO_BHMDOTRADIO,
  IO_BHVVIR,
  IO_BHTIMESTEP,
  IO_BHBPRESS,
  IO_BHHOSTHALOMASS,

  IO_BH_FRC_MINPOT,
  IO_BH_FRC_POTPOS,
  IO_BH_FRC_MINPOT_EXT,
  IO_BH_FRC_POTPOS_EXT,
  IO_BH_FRC_POTVEL,
  IO_BH_FRC_RHOTOT,

  IO_BHSPINPARAMETER,
  IO_BHSPINORIENTATION,
  IO_BHSPINMODEL,
  IO_BHANGMOMGASCELLS,

  IO_POT,
  IO_ACCEL,
  IO_GRADP,
  IO_GRADR,
  IO_GRADPCR,
  IO_GRADV,
  IO_GRADB,

  IO_POT_MINI,
  IO_POS_MINI,

  IO_CR_C0,
  IO_CR_Q0,
  IO_CR_P0,
  IO_CR_E0,
  IO_CR_n0,
  IO_CR_ThermalizationTime,
  IO_CR_DissipationTime,

  IO_ELECT,
  IO_HI,
  IO_HII,
  IO_HeI,
  IO_HeII,
  IO_HeIII,
  IO_H2I,
  IO_H2II,
  IO_HM,

  IO_TSTP,

  IO_BFLD,
  IO_AFLD,
  IO_CURLB,
  IO_BSMTH,
  IO_DIVB,
  IO_COOLRATE,
  IO_MACH,
  IO_PRESHOCK_DENSITY,
  IO_PRESHOCK_ENERGY,
  IO_PRESHOCK_XCR,
  IO_DENSITY_JUMP,
  IO_ENERGY_JUMP,

  IO_PSIR,
  IO_PSII,

  IO_nHI,
  IO_nHeI,
  IO_nHeII,
  IO_ColDens,
  IO_PHOTONS,
  IO_PHOTDENSITY,
  IO_SGCHEM_METALS,  /* Must always come before IO_CHEM */
  IO_METL,
  IO_CHEM,
  IO_COOLTIME,
  IO_TREECOLUMN,
  IO_GAMMA,
  IO_HEATRATE,
  IO_ESCFRAC,
  IO_EOSXNUC,

  IO_ALLOWREFINEMENT,

  IO_GFM_AGE,
  IO_GFM_INITIAL_MASS,
  IO_GFM_METALLICITY,
  IO_GFM_METALLICITY2,
  IO_GFM_METALS,
  IO_GFM_METALS2,
  IO_GFM_DUST_AGB,
  IO_GFM_DUST_AGB2,
  IO_GFM_DUST_SNII,
  IO_GFM_DUST_SNII2,
  IO_GFM_DUST_SNIa,
  IO_GFM_DUST_SNIa2,
  IO_GFM_DUST_METALLICITY,
  IO_GFM_DUST_METALLICITY2,
  IO_GFM_DUST_CAPMASS,
  IO_GFM_DUST_TAUGROWTH,
  IO_GFM_DUST_TAUSPUTTER,
  IO_GFM_LAMBDA,
  IO_GFM_CHEM_TAGS,
  IO_GFM_CHEM_TAGS2,
  IO_GFM_STELLAR_PHOTOMETRICS,
  IO_GFM_WINDHOSTMASS,
  IO_GFM_WINDHOSTDISP,
  IO_GFM_COOLRATE,
  IO_GFM_AGN_RADIATION,
  IO_ISM_CLUMP_POS,
  IO_ISM_CLUMP_DENSITY,
  IO_ISM_CLUMP_RADIUS,

  IO_DIVVEL,
  IO_CURLVEL,
  IO_COOLHEAT,
  IO_DUDT,
  IO_TEMP,
  IO_DUSTTEMP,
  IO_PASS,

  IO_VSTURB_DISS,
  IO_VSTURB_DRIVE,

  IO_SUBFINDHSML,
  IO_SUBFINDDENSITY,
  IO_SUBFINDDMDENSITY,
  IO_SUBFINDVELDISP,

  IO_COOL_DELAY,
  IO_BLAST_RADIUS,
  IO_HSML_STARS,
  IO_GFM_SNII_NUM,
  IO_GFM_CUM_SNII_NUM,
  IO_GFM_FEED_ENERGY,
  IO_GFM_FEED_MOMENTUM,
  IO_GFM_CUM_FEED_MOMENTUM,
  IO_GFM_MASS_RELEASED,

  IO_GFM_BIRTH_POS,
  IO_GFM_BIRTH_VEL,
  IO_GFM_BIRTH_RHO,

  IO_TURBULENT_ENERGY,
  IO_MOLECULAR_FRAC,
  IO_RADIATION_MOMENTUM,
  IO_RADIATION_KICK,
  IO_STROMGREN_RADIUS,
  IO_EARLY_STAR_FEEDBACK,
  IO_EARLY_STAR_FEEDBACK_MOMENTUM,
  IO_EARLY_STAR_FEEDBACK_CUM_MOMENTUM,
  IO_VAREFF_SN_METAL,
  IO_MASS_WEIGHT_SN,

  IO_NORMSPH_RADFEED,
  IO_COOL_DELAY_RADFEED,
  IO_TAU_RADFEED,
  IO_NUMNGB_RADFEED,
  IO_COOL_DELAY_GAS_RADFEED,
  IO_FEEDBACK_EVENTS,

  IO_OTVET_GAMMA,
  IO_OTVET_HI,
  IO_OTVET_HII,
  IO_OTVET_HeI,
  IO_OTVET_HeII,
  IO_OTVET_HeIII,
  IO_OTVET_ET,
  IO_OTVET_SRC_HSML,
  IO_OTVET_SRC_RHO,

  IO_FLD,
  IO_FLD_GRAD,
  IO_FLD_LAMBDA,
  IO_FLD_KAPPA_P,
  IO_FLD_KAPPA_R,
  IO_FLD_GAMMAS,

  IO_RP_MOM,

  IO_SIDM_PSUM,
  IO_SIDM_NUMNGB,
  IO_SIDM_NUMTOTALSCATTER,
  IO_SIDM_HSML,
  IO_SIDM_DENSITY,
  IO_SIDM_VELDISP,
  IO_SIDM_STATE,

  IO_GROUPNR,

  IO_STAR_INDEX,

  IO_RPS,
  IO_CRCHI,
  IO_CRENERGY,

  IO_SOFTENING,
  IO_GRAVINTERACTIONS,

  IO_SF_PROBABILITY,

  IO_LPOS,

  IO_PHIOS,
  IO_PHINS,
  IO_PHIT6,
  IO_PHIT7,
  IO_PHIT8,

  IO_TASK,
  IO_TIMEBIN_HYDRO,

  IO_DG_W0,
  IO_DG_W1,
  IO_DG_W2,
  IO_DG_W3,
  IO_DG_W4,
  IO_DG_T,
  IO_DG_U,
  IO_DG_ACCEL,
  IO_DG_DT,
  IO_DG_IB,
  IO_DG_DC,
  IO_DG_AM,
  IO_DG_SPIN,
  IO_DG_MN,
  IO_DG_MA_X,
  IO_DG_MA_Y,
  IO_DG_MA_Z,
  IO_DG_L1,

  IO_AMR_LEVEL,
  IO_AMR_REFFLAG,

  IO_SHOCK_MACHNUM,
  IO_SHOCK_EDISS,
  IO_SHOCK_GRAVFRAC,
  IO_SHOCK_AREA,
  IO_SHOCK_C_PRE,
  IO_SHOCK_RHO_PRE,
  IO_SHOCK_P_PRE,
  IO_SHOCK_RHO_POST,
  IO_SHOCK_P_POST,
  IO_SHOCK_T_PRE,
  IO_SHOCK_ZONE,
  IO_SHOCK_SURFACE,
  IO_SHOCK_ZONE_FLAG,
  IO_SHOCK_DIVVEL,
  IO_SHOCK_DIR,
  IO_SHOCK_V_POST,
  IO_SHOCK_V_PRE,

  IO_NUC_DEDT,

  IO_GAMMAC,
  IO_GAMMAE,

  IO_GRACKLE_TEMP,
  IO_GRACKLE_COOL_TIME,
  IO_GRACKLE_E,
  IO_GRACKLE_HI,
  IO_GRACKLE_HII,
  IO_GRACKLE_HeI,
  IO_GRACKLE_HeII,
  IO_GRACKLE_HeIII,
  IO_GRACKLE_HM,
  IO_GRACKLE_H2I,
  IO_GRACKLE_H2II,
  IO_GRACKLE_DI,
  IO_GRACKLE_DII,
  IO_GRACKLE_HDI,

  IO_SGCHEM_THERMAL_RATES,
  IO_SGCHEM_DUSTTOGAS,

  IO_DUST_GASHSML,
  IO_DUST_GASDENSITY,
  IO_DUST_GSD,
  IO_DUST_GSD_SLOPES,
  IO_DUST_METALFRACS,
  IO_DUST_DENSITY,
  IO_DUST_DUSTHSML,

  IO_LASTENTRY                  /* This should be kept - it signals the end of the list */
};

enum arrays
{
  A_NONE,
  A_SPHP,
#ifdef TRACER_MC
  A_TLL,
#endif
#ifdef GFM
  A_STARP,
#endif
#ifdef BLACK_HOLES
  A_BHP,
#endif
#ifdef DUST_LIVE
  A_DUSTP,
#endif
  A_P,
  A_PS
};

enum types_in_file
{
  FILE_NONE = -1,
  FILE_INT = 0,
  FILE_MY_ID_TYPE = 2,
  FILE_MY_IO_FLOAT = 1,
  FILE_DOUBLE = 3,
  FILE_FLOAT = 4
};

enum types_in_memory
{
  MEM_INT,
  MEM_MY_ID_TYPE,
  MEM_FLOAT,
  MEM_DOUBLE,
  MEM_MY_SINGLE,
  MEM_MY_FLOAT,
  MEM_MY_DOUBLE,
  MEM_NONE
};

enum e_typelist
{
  GAS_ONLY = 1,
  STARS_ONLY = 16,
  GAS_AND_STARS = 17,
  BHS_ONLY = 32,
#ifdef DUST_LIVE
  DUST_ONLY = (1 << DUST_LIVE),
#endif
#ifdef TRACER_MC
  TRACER_MC_PARENTS = 49, // gas, stars, BHs
  TRACER_MC_ONLY = 1<<TRACER_MC,
  ALL_TYPES = ((1<<NTYPES)-1) & ~(1<<TRACER_MC), // exclude MC tracers from e.g. IO_POS
#else
  ALL_TYPES = ((1<<NTYPES)-1),
#endif
  SET_IN_GET_PARTICLES_IN_BLOCK = 0
};

enum sn_type
{
  SN_FULL = 0,
  SN_MINI = 1,
  SN_MINI_ONLY = 2,
  SN_NO_SUBBOX = 3
};


typedef struct
{

  enum iofields field;
  enum types_in_memory type_in_memory;
  enum types_in_file type_in_file_input;
  enum types_in_file type_in_file_output;
  int values_per_block;
  char label[4];
  char datasetname[256];
  void (*io_func) (int, int, void*, int);
  int typelist;
  enum arrays array;
  size_t offset;
  enum sn_type snap_type;

  char hasunit;
  double a;
  double h;
  double L;
  double M;
  double V;
  double c;
}
IO_Field;

extern IO_Field *IO_Fields;
extern int N_IO_Fields;
extern int Max_IO_Fields;


extern char (*Parameters)[MAXLEN_PARAM_TAG];
extern char (*ParametersValue)[MAXLEN_PARAM_VALUE];
extern char *ParametersType;

#ifdef MODGRAV
typedef struct
 {
   unsigned int level:8; /*!< the level of the node, starting with 1 for the largest node, works for MaxAMRLevel < 128 */
   unsigned int num_daughters:4; /*!< number of regular daughter nodes  */
   unsigned int amr_node:2; /*!< 1 for lowest level nodes that have 7 siblings, 0 for nodes on coarser levels, 2 for nodes on finer levels */
   unsigned int red:1; /*!< flag for red-black sweep */
   unsigned int always_open:1;
 } MG_Bitflag_Data;
#endif


/** The tree data structure. Nodes points to the actual memory
    allocated for the internal nodes, but is shifted such that
    Nodes[All.MaxPart] gives the first allocated node. Note that node
    numbers less than All.MaxPart are the leaf nodes that contain a
    single particle, and node numbers >= MaxPart+MaxNodes are "pseudo
    particles" that hang off the toplevel leaf nodes belonging to
    other tasks. These are not represented by this structure. Instead,
    the tree traversal for these are saved in the Nextnode, Prevnode
    and Father arrays, indexed with the node number in the case of
    real particles and by nodenumber-MaxNodes for pseudo
    particles.  */
extern struct NODE
{
  union
  {
    int suns[8];                /**< temporary pointers to daughter nodes */
    struct
    {
      MyDouble s[3];            /**< center of mass of node */
      MyDouble mass;            /**< mass of node */
      /** The next node in the tree walk in case the current node does
          not need to be opened. This means that it traverses the 8
          subnodes of a node in a breadth-first fashion, and then goes
          to father->sibling. */
      int sibling;
      /** The next node in case the current node needs to be
          opened. Applying nextnode repeatedly results in a pure
          depth-first traversal of the tree. */
      int nextnode;
      /** The parent node of the node. (Is -1 for the root node.) */
      int father;
#if (NSOFTTYPES > 1)
      unsigned char maxsofttype;    /**< hold the maximum gravitational softening of particles */
#if defined(MULTIPLE_NODE_SOFTENING) && defined(ADAPTIVE_HYDRO_SOFTENING)
      unsigned char maxhydrosofttype;
      unsigned char minhydrosofttype;
#endif
#endif

#ifdef TREE_RAD_H2
      MyFloat h2mass;           /*!< mass of h2 in node */
      MyFloat comass;           /*!< mass of co in node */
#endif
    }
    d;
  }
  u;

#ifdef TREE_RAD_VEL
  MyFloat Vel[3];              /*!< centre of mass velocity of node */
#endif

#ifndef SIDM
  MyDouble center[3];           /**< geometrical center of node */
  MyFloat  len;                 /**< sidelength of treenode */
#else
  double center[3];           /**< geometrical center of node */
  double len;                 /**< sidelength of treenode */
#endif

#ifdef OTVET
  MyFloat stellar_s[3];         /*!< center of mass for the stars in the node */
  MyFloat stellar_mass;         /*!< mass in stars in the node */
#endif
#ifdef RADPRESS_OPT_THIN
  MyFloat young_stellar_s[3];   /*!< center of mass for the stars in the node */
  MyFloat young_stellar_mass;   /*!< mass in stars in the node */
#endif
#ifdef RADCOOL
  MyFloat young_stellar_s[3];   /*COM of young stars in the node */
  MyFloat young_stellar_mass;   /*mass of young stars in the node */
  MyFloat old_stellar_s[3];     /*COM of old stars in the node */
  MyFloat old_stellar_mass;     /*mass of old stars in the node */
#ifdef RADCOOL_HOTHALO
  MyFloat T6_gas_mass;          /*mass of gas with 5.5<logT<6.5 in the node */
  MyFloat T6_gas_s[3];          /*COM of gas */
  MyFloat T7_gas_mass;
  MyFloat T7_gas_s[3];
  MyFloat T8_gas_mass;
  MyFloat T8_gas_s[3];
#endif
#endif

#ifdef MODGRAV
  float mass_cic;

  double phi;
  double exp_phi;
#ifdef MODGRAV_INTERPOLATE_PHI
  MyFloat grad_phi[3];
#endif

  /* union which holds either effective mass / COM or variables for the field solver (only used in mg_fieldsolve())*/
  union
  {
    struct
    {
      double s_eff[3];
      float eff_mass;
    }
    eff; /* effective mass quantities */

    struct
    {
      double phi_mapped_up;   /*!< phi value that was mapped up from the finer level */
      double phi_orig;   /*!< original phi value at this level, not the value that comes out of the multigrid method at the coarse levels */
      double fsource;   /*!< source term for linear equation that is iteratively solved */
      float iterate;   /*!< whether this cell should be updated in iteration */
    }
    fs; /* variables for mg_fieldsolve() */
  }
  mg;

  double Lphi;   /*!< value of L phi */
  float len_orig;
  int nxm1, nxp1, nym1, nyp1, nzm1, nzp1;   /*!< neighbor nodes on same or coarser grid level */
  MG_Bitflag_Data mg_bitflag;
#endif
}
 *Nodes;


#ifdef MULTIPLE_NODE_SOFTENING
extern struct ExtNODE
{
  MyDouble mass_per_type[NSOFTTYPES];
}
 *ExtNodes;
#endif



/** Gives next node in tree walk for the "particle" nodes. Entries 0
    -- MaxPart-1 are the real particles, and the "pseudoparticles" are
    indexed by the node number-MaxNodes. */
extern int *Nextnode;
/** Gives previous node in tree walk for the leaf (particle)
    nodes. Entries 0 -- MaxPart-1 are the real particles, and the
    "pseudoparticles" are indexed by the node number-MaxNodes. */
extern int *Father;


/** Variables for neighbor tree */
extern int Ngb_MaxPart;
extern int Ngb_NumNodes;
extern int Ngb_MaxNodes;
extern int Ngb_FirstNonTopLevelNode;
extern int Ngb_NextFreeNode;
extern int *Ngb_Father;
extern int *Ngb_Marker;
extern int Ngb_MarkerValue;

extern int *Ngb_DomainNodeIndex;
extern int *DomainListOfLocalTopleaves;
extern int *DomainNLocalTopleave;
extern int *DomainFirstLocTopleave;
#ifdef AMR
extern int *Ngb_DomainTask;
#endif
extern int *Ngb_Nextnode;



/** The ngb-tree data structure
 */
extern struct NgbNODE
{
#ifndef AMR
  union
#else
  struct
#endif
  {
    int suns[8];                /**< temporary pointers to daughter nodes */
    struct
    {
      int sibling;
      int nextnode;
      MyNgbTreeFloat range_min[3];
      MyNgbTreeFloat range_max[3];
    }
    d;
  }
  u;

  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];

  int father;

  integertime Ti_Current;

#ifdef AMR
  int neighbors[2*NUMDIMS];
  int level;
  int nextinlevel;
  int previnlevel;
  MyFloat Center[3];
  amr_node_data hydro;
#endif
}
 *Ngb_Nodes;


extern struct ExtNgbNODE
{
  float vmin[3];
  float vmax[3];
  float MaxCsnd;
}
 *ExtNgb_Nodes;




#ifdef STATICNFW
extern double Rs, R200;
extern double Dc;
extern double RhoCrit, V200;
extern double fac;
#endif

#ifdef CIRCUMSTELLAR
#if (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined (EXTERNALGRAVITY)
extern struct source_particle_data
{
  unsigned char SoftningType;
  double Pos[3];
  double Mass;
  double SurfaceTemp;
  double SourceID;
}
 *SourcePartListGlobal;
#endif
#endif

#ifdef AB_TURB
//parameters
extern double StDecay;
extern double StEnergy;
extern double StDtFreq;
extern double StKmin;
extern double StKmax;
extern double StSolWeight;
extern double StAmplFac;

//Ornstein-Uhlenbeck variables
extern int StSeed;
extern double StOUVar;
extern double *StOUPhases;
extern gsl_rng *StRng;


//forcing field in fourie space
extern double *StAmpl;
extern double *StAka;           //phases (real part)
extern double *StAkb;           //phases (imag part)
extern double *StMode;
extern int StNModes;

extern integertime StTPrev;
extern double StSolWeightNorm;

extern int StSpectForm;

//extern FILE *FdTurb;
//extern double StEnergyAcc;
//extern double StEnergyDeacc;
//extern double StLastStatTime;
#endif


/** Struct holding the data for the rays used for making projections. */
typedef struct
{
  /** Index of the primary Voronoi cell in which the ray is located. */
  int index;
  /** Index of the previous cell the particle traversed. */
  int prev;
  /** Task of the primary Voronoi cell in which the ray is located. */
  int task;
  /** The index in the pixel array this ray is associated with. */
  int pixel;
  /** The length the ray has traveled. */
  double len;
  /** The total length the ray will travel when done. */
  double target_len;
  /** The current location of the ray. */
  MyDouble pos[3];
  /** The direction the ray is travelling in. */
  double dir[3];
} ray_data;

extern ray_data *Ray;

/** Local number of rays */
extern int Nray;
/** Maximum number of rays (equal to the number of pixels). */
extern int MaxNray;



extern int MaxThreads;


#endif
