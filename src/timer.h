/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/timer.h
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

#if !defined(TIMER_H) || defined(TIMER_STRUCT)
#define TIMER_H

#define DETAILED_TIMING_GRAVWALK         0
#define DETAILED_TIMING_STELLARDENSITY   1

#if defined(CUDA_INSTRUMENT)
#include <nvToolsExtCuda.h>
#define TIMER_INSTRUMENT_START(counter) nvtxRangePush(#counter);
#define TIMER_INSTRUMENT_STOP(counter)  nvtxRangePop();
#define TIMER_INSTRUMENT_CREATE(name, desc) ;
#elif defined(VTUNE_INSTRUMENT)
#include <ittnotify.h>
extern __itt_domain* vtune_domain;
#define TIMER_INSTRUMENT_START(counter) __itt_task_begin(vtune_domain, __itt_null, __itt_null, Timer_data[counter].vtune_string);
#define TIMER_INSTRUMENT_STOP(counter)  __itt_task_end(vtune_domain);
#define TIMER_INSTRUMENT_CREATE(name, descr) Timer_data[name].vtune_string = __itt_string_handle_create(descr);
#else
#define TIMER_INSTRUMENT_START(counter)
#define TIMER_INSTRUMENT_STOP(counter)
#define TIMER_INSTRUMENT_CREATE(name,descr);
#endif




#ifdef TIMER_STRUCT
#undef TIMER_CREATE
/*! \def TIMER_CREATE(name,desc, par, symba, symbb )
 * \brief creates a new CPU timer
 *
 * \param name name used in the code to reference this timer
 * \param desc description string used in output files
 * \param parent parent of this timer to build a tree-like hierarchy of timers
 * \param symba character used for active time in balance.txt
 * \param symbb character used for imbalance in balance.txt
 *
 */
#define TIMER_CREATE(name,desc, par, symba, symbb )  Timer_data[name].parent = par; \
         strncpy(Timer_data[name].shortname, #name, 40); \
         strncpy(Timer_data[name].longname,  (desc), 40); \
         Timer_data[name].symb = (symba); \
         Timer_data[name].symbImbal = (symbb); \
         TIMER_INSTRUMENT_CREATE(name, desc)

#else

#define TIMER_STACK_DEPTH 30

#define TIMER_CREATE(name,desc, parent,symba,symbb )  name ,


/*! \def  TIMER_START(counter)
 * \brief Starts the timer counter
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual  instrumentation APIs can be attached.
 *
 * \param counter Name of the timer to start
 */
#define TIMER_START_INTERNAL(counter) {\
    TIMER_INSTRUMENT_START(counter);\
    CPU_Step[TimerStack[TimerStackPos]] += measure_time();\
    int itimer; for(itimer=0;itimer<=TimerStackPos;itimer++) \
    if(counter==TimerStack[itimer]) \
      {printf("Try to start timer %d, but it is already running.\n",counter); terminate("fail")}; \
    if(++TimerStackPos >= TIMER_STACK_DEPTH)\
      {\
        terminate("Run out of timer stack space, increase TIMER_STACK_DEPTH");\
      }\
    else\
      {\
        TimerStack[TimerStackPos] = (counter);\
      }}

#define TIMER_START(counter) TIMER_START_INTERNAL(counter)

/*! \def TIMER_STOP(counter)
 * \brief Stops the timer counter
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual instrumentation APIs can be attached.
 *
 * \param counter Name of the timer to stop
 */
#define TIMER_STOP_INTERNAL(counter) {\
  if(TimerStack[TimerStackPos] != (counter))\
    {\
      terminate("Wrong use of TIMER_STOP, you must stop the timer started last");\
    }\
  CPU_Step[TimerStack[TimerStackPos--]] += measure_time();\
  if(TimerStackPos < 0)\
    {\
      terminate("Do not stop the out CPU_MISC timer");\
    }\
  TIMER_INSTRUMENT_STOP(counter); }

#define TIMER_STOP(counter) TIMER_STOP_INTERNAL(counter)


/*! \def TIMER_STOPSTART(stop, start)
 * \brief Stops the timer 'stop' and starts the timer 'start'
 *
 * Use this macro instead of directly accessing the CPU_Step array,
 * so manual instrumentation APIs can be attached.
 *
 * \param stop Name of the timer to stop
 * \param start Name of the timer to start
 */
#define TIMER_STOPSTART(stop,start) { TIMER_STOP_INTERNAL(stop); TIMER_START_INTERNAL(start); }


/*! \def TIMER_ADD(counter, amount)
 * \brief Adds amount to the timer counter

 * \param counter Name of the timer to add to
 * \param amount amount to add to timer counter
 */
#define TIMER_ADD(counter, amount)  CPU_Step[counter] += (amount);


/*! \def TIMER_DIFF(counter)
 * \brief Returns amount elapsed for the timer since last save with TIMER_STORE

 * \param counter Name of the timer to add to
 */
#define TIMER_DIFF(counter)  (CPU_Step[counter] - CPU_Step_Stored[counter])


/*! \def TIMER_STORE
 * \brief Copies the current value of CPU times to a stored variable, such that differences with respect to this reference can be calculated
 *
 */
#define TIMER_STORE  memcpy(CPU_Step_Stored, CPU_Step, sizeof(CPU_Step));


enum timers
{
  CPU_NONE = -2,                /*!< used for counters without a parent */
  CPU_ROOT = -1,                /*!< root node of the tree */
#endif

/* possible characters to use for marking the parts:
 *
 *   abdefghijklmnopqrstuvABCDEFGHHIJKLMNOPQRSTUV
 *   0123456789
 *   -:.*=[]^&;~/_$()?+"<>@#!|\
 */


/*add your counter here, they must appear in the right order*/

  TIMER_CREATE(CPU_ALL, "total", CPU_ROOT, '-', '-')    /*!< root timer, everything should be below this timer */
    TIMER_CREATE(CPU_TREE, "treegrav", CPU_ALL, 'a', ')')
    TIMER_CREATE(CPU_TREEBUILD, "treebuild", CPU_TREE, 'b', '(')
    TIMER_CREATE(CPU_TREEBUILD_INSERT, "insert", CPU_TREEBUILD, 'c', '*')
    TIMER_CREATE(CPU_TREEBUILD_BRANCHES, "branches", CPU_TREEBUILD, 'd', '&')
    TIMER_CREATE(CPU_TREEBUILD_TOPLEVEL, "toplevel", CPU_TREEBUILD, 'e', '^')
    TIMER_CREATE(CPU_TREECOSTMEASURE, "treecostm", CPU_TREE, 'f', '%')
    TIMER_CREATE(CPU_TREEWALK, "treewalk", CPU_TREE, 'g', '$')
    TIMER_CREATE(CPU_TREEWALK1, "treewalk1", CPU_TREEWALK, 'h', '#')
    TIMER_CREATE(CPU_TREEWALK2, "treewalk2", CPU_TREEWALK, 'i', '@')
    TIMER_CREATE(CPU_TREEBALSNDRCV, "treebalsndrcv", CPU_TREE, 'j', '!')
    TIMER_CREATE(CPU_TREESENDBACK, "treeback", CPU_TREE, 'm', '7')
    TIMER_CREATE(CPU_TREEDIRECT, "treedirect", CPU_TREE, 'r', '2')
#ifdef PMGRID
    TIMER_CREATE(CPU_PM_GRAVITY, "pm_grav", CPU_ALL, 's', '1')
#endif
    TIMER_CREATE(CPU_NGBTREEBUILD, "ngbtreebuild", CPU_ALL, 't', 'Z')
    TIMER_CREATE(CPU_NGBTREEUPDATEVEL, "ngbtreevelupdate", CPU_ALL, 'u', 'Y')
#ifdef VORONOI
    TIMER_CREATE(CPU_MESH, "voronoi", CPU_ALL, 'v', 'X')
    TIMER_CREATE(CPU_MESH_INSERT, "insert", CPU_MESH, 'w', 'W')
    TIMER_CREATE(CPU_MESH_FIND_DP, "findpoints", CPU_MESH, 'x', 'V')
    TIMER_CREATE(CPU_MESH_CELLCHECK, "cellcheck", CPU_MESH, 'y', 'U')
    TIMER_CREATE(CPU_MESH_GEOMETRY, "geometry", CPU_MESH, 'z', 'T')
    TIMER_CREATE(CPU_MESH_EXCHANGE, "exchange", CPU_MESH, 'A', 'S')
    TIMER_CREATE(CPU_MESH_DYNAMIC, "dynamic", CPU_MESH, 'B', 'R')
#endif
#ifdef AMR
    TIMER_CREATE(CPU_AMR, "amr", CPU_ALL, 'a', 'A')
    TIMER_CREATE(CPU_AMR_MESH, "amr mesh", CPU_AMR, 'm', 'M')
    TIMER_CREATE(CPU_AMR_EXCHANGE_CELLS, "amr exch cells", CPU_AMR, 'c', 'C')
    TIMER_CREATE(CPU_AMR_EXCHANGE_NODES, "amr exch nodes", CPU_AMR, 'n', 'N')
    TIMER_CREATE(CPU_AMR_EXCHANGE, "amr exchange", CPU_AMR, 'e', 'E')
    TIMER_CREATE(CPU_AMR_LINK_NGB, "amr link ngb", CPU_AMR, 'l', 'L')
    TIMER_CREATE(CPU_AMR_REFINEMENT, "amr refinement", CPU_AMR, 'l', 'L')
    TIMER_CREATE(CPU_AMR_UPDATE_NODES, "amr update nodes", CPU_AMR, 'u', 'U')
#ifdef AMR_CONNECTIONS
    TIMER_CREATE(CPU_MESH_DYNAMIC, "dynamic", CPU_AMR, 'B', 'R')
#endif
#endif
    TIMER_CREATE(CPU_HYDRO, "hydro", CPU_ALL, 'C', 'Q')
    TIMER_CREATE(CPU_GRADIENTS, "gradients", CPU_HYDRO, 'D', 'P')
#ifdef SECOND_DERIVATIVES
    TIMER_CREATE(CPU_HESSIAN, "hessians", CPU_HYDRO, 'E', 'O')
#endif
    TIMER_CREATE(CPU_FLUXES, "fluxes", CPU_HYDRO, 'F', 'N')
    TIMER_CREATE(CPU_FLUXES_COMM, "fluxcomm", CPU_HYDRO, 'H', 'L')
#ifdef VISCOSITY
    TIMER_CREATE(CPU_VISCOUS_FLUXES, "viscfluxes", CPU_HYDRO, 'I', 'K')
#endif
    TIMER_CREATE(CPU_CELL_UPDATES, "updates", CPU_HYDRO, 'J', 'j')
    TIMER_CREATE(CPU_SET_VERTEXVELS, "vertex vel", CPU_HYDRO, 'K', 'I')
    TIMER_CREATE(CPU_MHD, "mhd", CPU_HYDRO, '4', 'p')
#ifdef NUCLEAR_NETWORK
    TIMER_CREATE(CPU_NETWORK, "network", CPU_ALL, 'L', 'H')
    TIMER_CREATE(CPU_NETWORK_INTEGRATION, "integration", CPU_NETWORK, 'M', 'G')
    TIMER_CREATE(CPU_NETWORK_IMBALANCE, "imbalance", CPU_NETWORK, 'N', 'F')
#endif
#ifdef GFM
    TIMER_CREATE(CPU_GFM, "GFM", CPU_ALL, 'O', 'E')
    TIMER_CREATE(CPU_GFM_WINDS, "winds", CPU_GFM, 'P', 'D')
    TIMER_CREATE(CPU_GFM_ENRICH, "enrich", CPU_GFM, 'Q', 'C')
    TIMER_CREATE(CPU_GFM_AGNRAD, "agnrad", CPU_GFM, 'R', 'B')
    TIMER_CREATE(CPU_GFM_FEEDBACK, "stellarfeed", CPU_GFM, 'S', 'A')
#endif
#ifdef TRACER_PARTICLE
    TIMER_CREATE(CPU_TRACERS, "vel-tracers", CPU_ALL, 'T', 'z')
#endif
    TIMER_CREATE(CPU_DOMAIN, "domain", CPU_ALL, 'U', 'y')
    TIMER_CREATE(CPU_PEANO, "peano", CPU_ALL, 'V', 'x')
    TIMER_CREATE(CPU_DRIFTS, "drift/kicks", CPU_ALL, 'W', 'w')
    TIMER_CREATE(CPU_TIMELINE, "timeline", CPU_ALL, 'X', 'v')
#ifdef TREE_BASED_TIMESTEPS
    TIMER_CREATE(CPU_TREE_TIMESTEPS, "treetimesteps", CPU_ALL, 'Y', 'u')
#endif
    TIMER_CREATE(CPU_SNAPSHOT, "i/o", CPU_ALL, 'Z', 't')
    TIMER_CREATE(CPU_LOGS, "logs", CPU_ALL, '1', 's')
    TIMER_CREATE(CPU_COOLINGSFR, "sfrcool", CPU_ALL, '2', 'r')
#ifdef BLACK_HOLES
    TIMER_CREATE(CPU_BH, "blackholes", CPU_ALL, '3', 'q')
    TIMER_CREATE(CPU_BH_NGB, "ngb search", CPU_BH, '4', 'p')
    TIMER_CREATE(CPU_BH_DENSITY, "density", CPU_BH, '5', 'o')
    TIMER_CREATE(CPU_BH_ACCRETION, "accretion", CPU_BH, '6', 'n')
    TIMER_CREATE(CPU_BH_MERGERS, "mergers", CPU_BH, '4', 'p')
#endif
#ifdef DUST_LIVE
    TIMER_CREATE(CPU_DUST, "dust", CPU_ALL, 'd', 'D')
    TIMER_CREATE(CPU_DUST_DRAG, "drag", CPU_DUST, 'e', 'E')
    TIMER_CREATE(CPU_DUST_TRANSFER, "transfer", CPU_DUST, 'f', 'F')
    TIMER_CREATE(CPU_DUST_SHATTER, "shatter", CPU_DUST, 'g', 'G')
    TIMER_CREATE(CPU_DUST_PRODUCTION, "production", CPU_DUST, 'p', 'P')
#endif
#ifdef TGCHEM
    TIMER_CREATE(CPU_TGCHEM, "tgchem", CPU_ALL, 's', 'O')
#endif
#ifdef HEALRAY
    TIMER_CREATE(CPU_HEALRAY, "healray", CPU_ALL, 'R', 'S')
#endif
#if defined(SINKS) || defined(SINK_PARTICLES)
    TIMER_CREATE(CPU_SINKS, "sinks", CPU_ALL, 'M', 'm')
#endif
#ifdef MONOTONE_CONDUCTION
    TIMER_CREATE(CPU_CONDUCTION, "conduction", CPU_ALL,  '{', '}')
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
    TIMER_CREATE(CPU_OHM, "ohm diffusion", CPU_ALL,  '{', '}')
#endif
#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT) || defined(SHOCK_FINDER_ON_THE_FLY)
    TIMER_CREATE(CPU_SHOCK_FINDER, "shock-finder", CPU_ALL, '7', 'm')
    TIMER_CREATE(CPU_SF_OUTPUT, "output", CPU_SHOCK_FINDER, '8', 'l')
    TIMER_CREATE(CPU_SF_ZONE, "shock zone", CPU_SHOCK_FINDER, '9', 'k')
    TIMER_CREATE(CPU_SF_SURFACE, "shock surface", CPU_SHOCK_FINDER, '!', 'j')
    TIMER_CREATE(CPU_SF_MISC, "misc", CPU_SHOCK_FINDER, '@', 'i')
#endif
#ifdef DG
    TIMER_CREATE(CPU_DISCONTINUOUS_GALERKIN,"discont. galerkin",CPU_ALL,'S','s')
    TIMER_CREATE(CPU_DG_INNER,"R inner",CPU_DISCONTINUOUS_GALERKIN,'T','t')
    TIMER_CREATE(CPU_DG_OUTER,"R_outer",CPU_DISCONTINUOUS_GALERKIN,'U','u')
    TIMER_CREATE(CPU_DG_FLEXCH,"flux exch.",CPU_DG_OUTER,'F','f')
    TIMER_CREATE(CPU_DG_FLIMBAL,"flux imbal.",CPU_DG_OUTER,'F','f')
    TIMER_CREATE(CPU_DG_LIMITER,"limiter",CPU_DISCONTINUOUS_GALERKIN,'V','v')
    TIMER_CREATE(CPU_DG_RECOMP,"recomp node",CPU_DISCONTINUOUS_GALERKIN,'R','r')
    TIMER_CREATE(CPU_DG_EXCHANGE,"exchange",CPU_DISCONTINUOUS_GALERKIN,'E','e')
#endif
#ifdef AB_TURB
    TIMER_CREATE(CPU_TURB,"turb",CPU_ALL,'T','t')
    TIMER_CREATE(CPU_TURB_RESET,"turb_reset",CPU_TURB,'T','t')
    TIMER_CREATE(CPU_TURB_UPDATE,"turb_update",CPU_TURB,'T','t')
    TIMER_CREATE(CPU_TURB_FORCE,"turb_force",CPU_TURB,'T','t')
#endif
#ifdef FLD
    TIMER_CREATE(CPU_FLD,"fld",CPU_ALL,'F','f')
    TIMER_CREATE(CPU_FLD_MATRIX,"fld_matrix",CPU_FLD,'M','m')
#endif
#ifdef COSMIC_RAYS_DIFFUSION
    TIMER_CREATE(CPU_CR_DIFFUSION, "CR diffusion", CPU_ALL, '{', '}')
    TIMER_CREATE(CPU_CR_DIFFUSION_PREPARE, "preparation", CPU_CR_DIFFUSION, '{', '}')
    TIMER_CREATE(CPU_CR_DIFFUSION_SOLVE, "solver", CPU_CR_DIFFUSION, '{', '}')
    TIMER_CREATE(CPU_CR_DIFFUSION_SOLVE_EXPLICIT, "explicit flux", CPU_CR_DIFFUSION_SOLVE, '{', '}')
    TIMER_CREATE(CPU_CR_DIFFUSION_SOLVE_MATRIX, "matrix setup", CPU_CR_DIFFUSION_SOLVE, '{', '}')
    TIMER_CREATE(CPU_CR_DIFFUSION_SET_COEFF, "set coeff", CPU_CR_DIFFUSION_SOLVE, '{', '}' )
    TIMER_CREATE(CPU_CR_DIFFUSION_SOLVE_IMPLICIT, "implicit solve", CPU_CR_DIFFUSION_SOLVE, '{', '}')
#endif
#ifdef COSMIC_RAYS_STREAMING
    TIMER_CREATE(CPU_CR_STREAMING, "CR streaming", CPU_ALL, '{', '}')
#endif
#ifdef FOF
    TIMER_CREATE(CPU_FOF, "fof", CPU_ALL, '#', 'h')
#endif
#ifdef SUBFIND
    TIMER_CREATE(CPU_SUBFIND, "subfind", CPU_ALL, '$', 'g')
#endif
#ifdef VORONOI
    TIMER_CREATE(CPU_REFINE, "refine", CPU_ALL, '%', 'f')
    TIMER_CREATE(CPU_DEREFINE, "mesh_derefine", CPU_ALL, '^', 'e')
#endif
#ifdef SIDM
    TIMER_CREATE(CPU_SIDM, "sidm", CPU_ALL, '+', '+')
    TIMER_CREATE(CPU_SIDM_ALLOCFREE, "allocfree", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_HSML, "hsml", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_CHECK, "checkscatter", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_NGB, "ngb", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_ASSIGN, "assignscatter", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_SCATTER, "scatter", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_EXCHANGE, "exchange", CPU_SIDM, '+', '+')
    TIMER_CREATE(CPU_SIDM_STATS, "stats", CPU_SIDM, '+', '+')
#endif
    TIMER_CREATE(CPU_MAKEIMAGES, "images", CPU_ALL, '&', 'd')
    TIMER_CREATE(CPU_INIT, "initializ.", CPU_ALL, '*', 'c')
    TIMER_CREATE(CPU_RESTART, "restart", CPU_ALL, '(', 'b')
    TIMER_CREATE(CPU_MISC, "misc", CPU_ALL, ')', 'a')
    TIMER_CREATE(CPU_LAST, "LAST", CPU_NONE, ' ', ' ')   /*!<last item, do not use! */
#ifndef TIMER_STRUCT
};

extern enum timers TimerStack[TIMER_STACK_DEPTH];
extern int TimerStackPos;

/*! \brief struct containing the information of a CPU timer
 *
 */
struct timer_d
{

  int parent;                   /*!< id of the parent timer */
  char shortname[40];           /*!< string containing the internal name of the timer */
  char longname[40];            /*!< name of the timer */
  char symb;                    /*!< symbol used in balance.txt for the active part */
  char symbImbal;               /*!< symbol used in balance.txt for imbalance */
  char depth;                   /*!< depth in the tree-like structure of this timer */
#ifdef VTUNE_INSTRUMENT
  __itt_string_handle* vtune_string;
#endif
};
extern struct timer_d Timer_data[CPU_LAST + 1];
#else
#undef TIMER_STRUCT
#endif
#endif
