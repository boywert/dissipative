/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/begrun.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"
#include "domain.h"
#include "MRT/RT.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
herr_t my_hdf5_error_handler(void *unused);
#endif

/*! \file begrun.c
 *  \brief Initial set-up of a simulation run
 *
 *  This file contains various functions to initialize a simulation run. In
 *  particular, the parameter file is read in and parsed and global variables
 *  are initialized to their proper values.
 */

void hello(void)
{
#ifndef DG
  mpi_printf("\n   __    ____  ____  ____  _____\n  /__\\  (  _ \\( ___)(  _ \\(  _  )\n /(__)\\  )   / )__)  )___/ )(_)(\n(__)(__)(_)\\_)(____)(__)  (_____)\n\n");
#else
  mpi_printf("  _____   U _____ u   _   _     U _____ u   _____      \n |_ \" _|  \\| ___\"|/  | \\ |\"|    \\| ___\"|/  |_ \" _|     \n"
             "   | |     |  _|\"   <|  \\| |>    |  _|\"      | |       \n  /| |\\    | |___   U| |\\  |u    | |___     /| |\\      \n"
             " u |_|U    |_____|   |_| \\_|     |_____|   u |_|U      \n _//\\\\_    <<   >>   ||   \\\\,-.  <<   >>   _// \\\\_     \n" "(__) (__) (__) (__)  (_\")  (_/  (__) (__) (__) (__)    \n");
#endif
}



/*! \brief Print a welcome message and used compile options
 *
 */
void begrun0(void)
{
#ifndef DG
  mpi_printf
    ("\nThis is Arepo, version %s.\n\nRunning with %d MPI tasks.\n\nApparently we're using %d compute nodes (we have a minimum of %d MPI tasks per node, and a maximum of %d)\n\nCode was compiled with settings:\n\n",
     AREPO_VERSION, NTask, NumNodes, MinTasksPerNode, MaxTasksPerNode);
#else
  mpi_printf
    ("\nThis is Tenet, version %s.\n\nRunning with %d MPI tasks.\n\nApparently we're using %d compute nodes (we have a minimum of %d MPI tasks per node, and a maximum of %d)\nCode was compiled with settings:\n\n",
     TENET_VERSION, NTask, NumNodes, MinTasksPerNode, MaxTasksPerNode);
#endif

  if(ThisTask == 0)
    output_compile_time_options();
}

/*! \brief This function performs the initial set-up of the simulation.
 *
 *  First, the parameter file is read by read_parameter_file(),
 *  then routines for setting units, etc are called. This function only does
 *  the setup necessary to load the IC file. After the IC file has been loaded
 *  and prepared by init(), setup continues with begrun2(). This splitting is
 *  done so that we can return cleanly from operations that don't actually
 *  start the simulation (converting snapshots, making projected images, etc.)
 */
void begrun1(void)
{
#if defined(X86FIX) && defined(SOFTDOUBLEDOUBLE)
  x86_fix();                    /* disable 80bit treatment of internal FPU registers in favour of proper IEEE 64bit double precision arithmetic */
#endif

  read_parameter_file(ParameterFile);   /* ... read in parameters for this run */

  check_parameters(); /* consistency check of parameters */

#ifdef HAVE_HDF5
  H5Eset_auto(my_hdf5_error_handler, NULL);
#endif

  gsl_set_error_handler(my_gsl_error_handler);

#ifdef CUDA
  cuda_init();
#endif

#if defined(WINDTUNNEL) && defined(WINDTUNNEL_EXTERNAL_SOURCE)
  read_windtunnel_file();
#endif

#ifdef DMPIC
  dmpic_init();
#endif

#ifdef DEBUG
  enable_core_dumps_and_fpu_exceptions();
#endif

  mpi_printf("BEGRUN: Size of particle structure       %3d  [bytes]\n", (int) sizeof(struct particle_data));
  mpi_printf("BEGRUN: Size of sph particle structure   %3d  [bytes]\n", (int) sizeof(struct sph_particle_data));
  mpi_printf("BEGRUN: Size of gravity tree node        %3d  [bytes]\n", (int) sizeof(struct NODE));
#ifdef MULTIPLE_NODE_SOFTENING
  mpi_printf("BEGRUN: Size of auxiliary gravity node   %3d  [bytes]\n", (int) sizeof(struct ExtNODE));
#endif
#ifdef GFM
  mpi_printf("BEGRUN: Size of star particle structure  %3d  [bytes]\n", (int) sizeof(struct star_particle_data));
#endif
#ifdef BLACK_HOLES
  mpi_printf("BEGRUN: Size of BH particle structure    %3d  [bytes]\n\n", (int) sizeof(struct bh_particle_data));
#endif
#ifdef DUST_LIVE
  mpi_printf("BEGRUN: Size of dust particle structure  %3d  [bytes]\n\n", (int) sizeof(struct dust_particle_data));
#endif


#ifdef DARKENERGY
#ifdef TIMEDEPDE
/* set up table needed for hubble functions with time dependent w */
  fwa_init();
#endif
#endif

  set_units();

  if(RestartFlag == 1) /* this is needed here to allow domain decomposition right after restart */
    if(All.ComovingIntegrationOn)
      init_drift_table();

  init_io_fields();

  force_short_range_init();

#if defined (FORCETEST) && !defined(FORCETEST_TESTFORCELAW) && !defined(GRAVITY_TALLBOX)
  forcetest_ewald_init();
#endif

  /* set up random number generators */
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  random_generator_aux = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* individual start-up seed */
  gsl_rng_set(random_generator, 42 + ThisTask);
  gsl_rng_set(random_generator_aux, 31452 + ThisTask);


#ifdef PERFORMANCE_TEST_SPARSE_MPI_ALLTOALL
  test_mpi_alltoall_performance();
#endif

  timebins_init(&TimeBinsHydro, "Hydro", &All.MaxPartSph);
  timebins_init(&TimeBinsGravity, "Gravity", &All.MaxPart);

#ifdef TRACER_PARTICLE
  timebins_init(&TimeBinsTracer, "Tracer", &All.MaxPart);       /* assume that there are less tracers than cells */
#endif

#ifdef BLACK_HOLES
  timebins_init(&TimeBinsBHAccretion, "BHAccretion", &All.MaxPart);
#endif

#ifdef SNE_FEEDBACK
  sne_init();
#endif

#ifdef SINKS
  timebins_init(&TimeBinsSinksAccretion, "SinksAccretion", &All.MaxPart);
#endif

#ifdef DUST_LIVE
  timebins_init(&TimeBinsDust, "Dust", &All.MaxPart);
#endif

#ifdef SUBBOX_SNAPSHOTS
  read_subbox_coordinates(All.SubboxCoordinatesPath);
#endif


#ifdef GFM_STELLAR_EVOLUTION
  /* copy element names to ElementNames */
  strcpy(ElementNames[0], "Hydrogen");
  strcpy(ElementNames[1], "Helium");
  strcpy(ElementNames[2], "Carbon");
  strcpy(ElementNames[3], "Nitrogen");
  strcpy(ElementNames[4], "Oxygen");
  strcpy(ElementNames[5], "Neon");
  strcpy(ElementNames[6], "Magnesium");
  strcpy(ElementNames[7], "Silicon");
  strcpy(ElementNames[8], "Iron");
#ifdef GFM_NORMALIZED_METAL_ADVECTION
  strcpy(ElementNames[9], "OtherMetals");
#endif




  init_imf();
  mpi_printf("GFM_STELLAR_EVOLUTION: IMF initialized.\n");

  init_yields();
  mpi_printf("GFM_STELLAR_EVOLUTION: Yields initialized.\n");

#ifdef RADCOOL
  init_radcool_units();
  mpi_printf("RADCOOL, RADCOOL_HOTHALO : Mass and length units converted to Msun and kpc respectively.\n");
#endif

#ifdef GFM_COOLING_METAL
  init_cooling_metal();
  mpi_printf("GFM_COOLING_METAL: Metal line cooling rates initialized.\n");
#endif

#ifdef GFM_WINDS
  init_winds();
#endif

#ifdef GFM_WINDS_VARIABLE
  if(ThisTask == 0)
    init_variable_winds();
#endif

#ifdef GFM_PREENRICH
  gfm_read_preenrich_table(All.PreEnrichAbundanceFile);
#endif

#endif

#ifdef GFM_STELLAR_PHOTOMETRICS
  init_stellar_photometrics();
  mpi_printf("GFM_STELLAR_PHOTOMETRICS: Stellar photometrics initialized.\n");
#endif

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)
  init_agn_radiation();
  mpi_printf("GFM_AGN_RADIATION: initialized.\n");
#endif

#if defined(TRACER_PART_NUM_FLUID_QUANTITIES) || defined(TRACER_MC_NUM_FLUID_QUANTITIES)
  set_tracer_part_indices();
#endif

#ifdef GROWING_DISK_POTENTIAL
  growing_disk_init();
#endif

#if defined(COOLING) & !defined(SIMPLE_COOLING) & !defined(GRACKLE)
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  InitCool();
#endif

#if defined(COOLING) && defined(GRACKLE)
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  initialise_grackle();
#endif

#ifdef ATOMIC_DM
  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
  ADM_InitCool();
#endif

#ifdef TGSET
  tgset_begrun();
#endif

#ifdef HEALRAY
  healray_begrun();             // Needs to be before TGCHEM!
#endif

#ifdef TGCHEM
  tgchem_begrun();
#endif

#ifdef SGCHEM
  init_chemistry();
#endif

#ifdef TREE_RAD

  mpi_printf("\n Calculating pixel centres \n");
  CALCULATE_PIXEL_CENTRES();
  CREATETRIGLOOKUP();

#if defined(TREE_RAD_H2) || defined(TREE_RAD_CO)
#ifndef SGCHEM
  fprintf(stderr, "Error: must define SGCHEM\n");
  endrun(1213);
#endif
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 15
#define CO_OK
#endif
  /* Check that Treecol setup makes sense */
#ifdef TREE_RAD_CO
#ifndef CO_OK
  fprintf(stderr, "Error: chemistry network does not contain CO, so TreeCol cannot compute CO columns\n");
  endrun(1213);
#endif
#endif /* TREE_RAD_CO */
#undef CO_OK
#endif /* TREE_RAD_H2 || TREE_RAD_CO */
#endif /* TREE_RAD */

#ifdef SINKS
  sinks_begrun();
#endif

#ifdef SINK_PARTICLES
  init_sink_particles(0);
#endif

#if !defined(PMGRID) && defined(SELFGRAVITY) && !defined(GRAVITY_NOT_PERIODIC) && defined(PERIODIC) && !defined(ONEDIMS_SPHERICAL)
  ewald_init();
#endif

#ifdef TILE_ICS
  All.BoxSize *= All.TileICsFactor;
#endif

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

  EgyInjection = 0;

#ifdef PMGRID
  if((RestartFlag != 3) && (RestartFlag != 4) && (RestartFlag != 6) && (RestartFlag != 7) && (RestartFlag != 11))
    long_range_init();
#endif

  if(RestartFlag <= 2)
    open_logfiles();

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif

#ifdef EOS_DEGENERATE
#ifndef VARIABLE_GAMMA
#error "EOS_DEGENERATE requires VARIABLE_GAMMA"
#endif
  eos_init(All.EosTable, All.EosSpecies);
#endif

#ifdef NUCLEAR_NETWORK
  network_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses, All.NetworkWeakrates, &All.nd);

  {
    int k;
    for(k = 0; k < NUM_THREADS; k++)
      network_workspace_init(&All.nd, &All.nw[k]);
  }

#ifdef NETWORK_NSE
  network_nse_init(All.EosSpecies, All.NetworkRates, All.NetworkPartFunc, All.NetworkMasses, All.NetworkWeakrates, &All.nd_nse, All.nw_nse);
#endif

#endif

#ifdef EOS_OPAL
#ifndef VARIABLE_GAMMA
#error "EOS_OPAL requires VARIABLE_GAMMA"
#endif
  if(opaleos_init(All.EosOpalTable) < 0)
    terminate("Error in OPAL EOS initialization!\n");
  /* get rho boundaries */
  opaleos_get_rho_limits(&opal_rhomin, &opal_rhomax);
#endif

#ifdef GENERAL_RELATIVITY
#if METRIC_TYPE==3 || METRIC_TYPE==4
 if ( (read_fixed_numerical_1d_metric()) != 0 )
 {
      printf("problems with reading fixed metric");
      terminate("bled");
 }

#endif
#endif

#ifdef COSMIC_RAYS
  init_cosmic_rays();
#endif

#if defined(DUST_LIVE) && defined(DL_GRAIN_BINS)
  init_dust();
#endif

  init_scalars();

  init_gradients();
#ifdef MRT
  init_gradients_RT() ;
#endif

#ifdef RT_ADVECT
  rt_init_gradients();
#endif

#ifdef SECOND_DERIVATIVES
  init_hessians();
#endif

#ifdef AMR
  amr_init();
#endif

#if defined (VS_TURB) || defined (AB_TURB)
  init_turb();
#ifdef POWERSPEC_GRID
  powersepc_turb_init();
#endif
#endif

#ifdef SIDM
  sidm_Init_CrossSection();
#endif

#ifdef DG
  dg_initialize();
#endif

#ifdef GRAVITY_TABLE  
  grav_table_init();
#endif

#ifdef RELAXOBJECT_COOLING2
  load_temperature_profil();
#endif

}


/*! \brief This function does late setup, after the IC file has been loaded
 *  but before run() is called.
 *
 *  The output files are opened and various modules are initialized. The next output
 *  time is determined by find_next_outputtime() and various timers are set.
 *
 */
void begrun2(void)
{
  char contfname[1000];
  sprintf(contfname, "%scont", All.OutputDir);
  unlink(contfname);

  if(RestartFlag > 2)
    open_logfiles();

#if defined(DG_SET_IC_FROM_AVERAGES) && defined(DG)
  load_weights_from_averages();
#endif

#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK)
   sfr_init();
#endif

#ifdef GFM_STELLAR_EVOLUTION
  init_SNIa_rates();
  mpi_printf("GFM_STELLAR_EVOLUTION: Type Ia rates initialized.\n");
#endif

#ifdef PMGRID
  long_range_init_regionsize();
#endif

#ifdef RT_ADVECT
  if(RestartFlag == 0)
    {
      rt_set_simple_inits();
      rt_init_sourceid();
#ifdef RT_STELLAR_SOURCES
      rt_create_source_list();
#endif
    }
#ifdef RT_HEALPIX_NSIDE
  rt_get_vectors();
#endif
#endif

#ifdef EXACT_GRAVITY_FOR_PARTICLE_TYPE
  special_particle_create_list();
#endif

#ifdef REFINEMENT_AROUND_DM
  dm_particle_create_list();
#endif

#if (defined(CIRCUMSTELLAR_IRRADIATION) || defined(ALPHA_VISCOSITY) || defined(CIRCUMSTELLAR_REFINEMENTS)) && !defined (EXTERNALGRAVITY)
  source_particle_create_list();
#endif

  if(RestartFlag != 1) /* this needs to be done here because here All.TimeBegin has the correct value */
    if(All.ComovingIntegrationOn)
      init_drift_table();

#ifdef TGSET
  if(All.TimeBetSnapshot)
#endif
    {
      if(RestartFlag == 2)
        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 100);
      else
        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);
    }

#ifdef OTVET
  if(RestartFlag == 0)
    otvet_set_simple_inits();

  ot_get_sigma();

#ifdef OTVET_MULTI_FREQUENCY
#if defined(EDDINGTON_TENSOR_STARS)
  ot_get_lum_stars();
#endif
#endif
#endif


#ifdef MRT
  init_RT() ;
#endif


#if defined(MRT_CHEMISTRY_PS2011) || defined(MRT_CHEMISTRY_PS2009)
  mrt_get_sigma();
#endif

#ifdef TRACER_TRAJECTORY
  tracer_init_output_configuration();
#endif

  All.TimeLastRestartFile = CPUThisRun;

#ifdef REDUCE_FLUSH
  All.FlushLast = CPUThisRun;
#endif

#if defined(FORCETEST) && defined(FORCETEST_TESTFORCELAW)
  gravity_forcetest_testforcelaw();
#endif

#ifdef SHOCK_FINDER_ON_THE_FLY  //create full mesh and run shock finder
  if(RestartFlag == 1 || RestartFlag == 2)
    {
      int k;
      short int *buTimeBin = mymalloc_movable(&buTimeBin, "buTimeBin", NumPart * sizeof(short int));
      static int buTimeBinActive[TIMEBINS];

      for(k = 0; k < NumPart; k++)
        {
          buTimeBin[k] = P[k].TimeBinHydro;
          P[k].TimeBinHydro = 0;
        }

      for(k = 0; k < TIMEBINS; k++)
        {
          buTimeBinActive[k] = TimeBinSynchronized[k];

          TimeBinSynchronized[k] = 1;
        }

      reconstruct_timebins();

      create_mesh();
      mesh_setup_exchange();
      shock_finder_on_the_fly();
      free_mesh();

      for(k = 0; k < TIMEBINS; k++)
        TimeBinSynchronized[k] = buTimeBinActive[k];

      for(k = 0; k < NumPart; k++)
        P[k].TimeBinHydro = buTimeBin[k];

      reconstruct_timebins();

      myfree_movable(buTimeBin);
    }
#endif
}




/*! \brief Computes conversion factors between internal code units and the
 *  cgs-system.
 *
 *  In addition constants like the gravitation constant are set.
 */
void set_units(void)
{
  double meanweight;

#ifdef STATICNFW
  double Mtot;
#endif

  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitTime_in_Megayears = All.UnitTime_in_s / SEC_PER_MEGAYEAR;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;

#ifdef BECDM_H
    All.hbar = 1.0545718e-27 / pow(All.UnitLength_in_cm, 2) / All.UnitMass_in_g * All.UnitTime_in_s;
    All.mAxion = All.AxionMassEv * 1.78266191e-33/ All.UnitMass_in_g;
#endif

  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitPressure_in_cgs = All.UnitMass_in_g / All.UnitLength_in_cm / pow(All.UnitTime_in_s, 2);
  All.UnitCoolingRate_in_cgs = All.UnitPressure_in_cgs / All.UnitTime_in_s;
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  mpi_printf("BEGRUN: Hubble (internal units)   = %g\n", All.Hubble);
  mpi_printf("BEGRUN: G (internal units)        = %g\n", All.G);
  mpi_printf("BEGRUN: UnitMass_in_g             = %g\n", All.UnitMass_in_g);
  mpi_printf("BEGRUN: UnitTime_in_s             = %g\n", All.UnitTime_in_s);
  mpi_printf("BEGRUN: UnitVelocity_in_cm_per_s  = %g\n", All.UnitVelocity_in_cm_per_s);
  mpi_printf("BEGRUN: UnitDensity_in_cgs        = %g\n", All.UnitDensity_in_cgs);
  mpi_printf("BEGRUN: UnitEnergy_in_cgs         = %g\n", All.UnitEnergy_in_cgs);
  mpi_printf("\n");

  meanweight = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC);       /* note: assuming NEUTRAL GAS */

  if(All.MinEgySpec == 0)
    {
      All.MinEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.MinGasTemp;
      All.MinEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;

      mpi_printf("BEGRUN: MinEgySpec set to %g based on MinGasTemp=%g\n", All.MinEgySpec, All.MinGasTemp);
    }

#ifdef FM_STAR_FEEDBACK
    All.PhotoionizationEgySpec = 1 / meanweight * (1.0 / GAMMA_MINUS1) * (BOLTZMANN / PROTONMASS) * All.PhotoionizationGasTemp;
    All.PhotoionizationEgySpec *= All.UnitMass_in_g / All.UnitEnergy_in_cgs;
#endif

#if defined(USE_SFR) && !defined(FM_SFR) && !defined(ISM) && !defined(LOCAL_FEEDBACK)
  set_units_sfr();
#endif

#ifdef CONDUCTION               //VITALI
  init_spitzer_conductivity();
#endif

#ifdef MRT
  init_RT() ;
#endif


#ifdef MONOTONE_CONDUCTION
  init_conductivity();
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
  init_ohm_conductivity();
#endif


#ifdef STATICNFW
  R200 = pow(NFW_M200 * All.G / (100 * All.Hubble * All.Hubble), 1.0 / 3);
  Rs = R200 / NFW_C;
  Dc = 200.0 / 3 * NFW_C * NFW_C * NFW_C / (log(1 + NFW_C) - NFW_C / (1 + NFW_C));
  RhoCrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
  V200 = 10 * All.Hubble * R200;
  mpi_printf("V200= %g\n", V200);

  fac = 1.0;
  Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
  fac = V200 * V200 * V200 / (10 * All.G * All.Hubble) / Mtot;
  Mtot = enclosed_mass(R200);
  mpi_printf("M200= %g\n", Mtot);
#endif
}
