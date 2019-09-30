/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/init.c
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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"
#include "domain.h"

#ifdef CALCULATE_QUANTITIES_IN_POSTPROCESS
#include "calculate_quantities_in_postprocess/calculate_quantities.h"
#endif

/*! \file init.c
 *  \brief Code for initialization of a simulation from initial conditions
 */

/*! \brief Prepares the loaded initial conditions for the run
 *
 *  It is only called if RestartFlag !=1. Various counters and variables are initialized.
 *  Entries of the particle data structures not read from initial conditions are
 *  initialized or converted and a initial domain decomposition is performed.
 *  If SPH particles are present, the initial SPH smoothing lengths are determined.
 *
 *  \return status code: <0 if finished without errors and run can start, 0 code ends after calling init() (used for
 *  special run modes like image generation), > 0 an error occurred, terminate
 */
int init(void)
{
  int i, j;

#ifdef BLACK_HOLES
#ifdef BH_RELATIVE_NGB_DEVIATION
  All.MaxNumNgbDeviationBlackHole = All.DesNumNgbBlackHole * All.DesNumNgbBlackHoleRelDeviationFactor;
#else
  All.MaxNumNgbDeviationBlackHole = All.MaxNumNgbDeviation;
#endif
#endif

#ifdef BLACK_HOLES
  int count_holes = 0;
#endif
#ifdef GFM_STELLAR_EVOLUTION
  MyDouble mass_fractions[GFM_N_CHEM_ELEMENTS];
#endif

  double mass;

  assert(RestartFlag != 1);


#ifndef TGSET
  if(All.ComovingIntegrationOn)
    if(All.PeriodicBoundariesOn == 1)
      {
        if(RestartFlag < 3)
          /* can't do this check when not all particles are loaded */
          check_omega();
        else
          mpi_printf("INIT: Skipping Omega check since we are not doing a dynamical evolution (not all particles may be loaded)\n");
      }
#endif

#ifdef GALPOT
  mpi_printf("Initializing potential\n");
  galpot_init();
#endif

#if defined(COOLING) && !defined(GRACKLE)
  IonizeParams();
#ifdef GFM_COOLING_METAL
  read_cooling_tables_current_time();
#endif
#endif

#ifdef GFM_STELLAR_EVOLUTION 
  get_initial_mass_fractions(&mass_fractions[0]);
#endif

  if(All.ComovingIntegrationOn)
    {
      All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
      All.Ti_Current = 0;
    }
  else
    {
      All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
      All.Ti_Current = 0;
    }

  set_cosmo_factors_for_current_time();

  All.NumCurrentTiStep = 0;     /* setup some counters */
  All.SnapshotFileCount = 0;
#ifdef SUBBOX_SNAPSHOTS
  All.SubboxSnapshotFileCount = 0;
#endif
#ifdef MULTIPLE_RESTARTS
  All.RestartFileCount = 0;
#ifdef SUBBOX_SNAPSHOTS
  All.SubboxSyncCounter = 0;
#endif
#endif

  if(RestartFlag == 2)
    {
      if(RestartSnapNum < 0)
        All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;
      else
        All.SnapshotFileCount = RestartSnapNum + 1;
    }

  if(RestartFlag == 10)
    All.SnapshotFileCount = RestartSnapNum;

  All.TotNumOfForces = 0;
  All.TopNodeAllocFactor = 0.08;
  All.TreeAllocFactor = 0.7;
  All.NgbTreeAllocFactor = 0.7;

  if(NumPart < 1000)
    All.TreeAllocFactor = 10.0;

#ifdef VORONOI                  /* TODO: move to voronoi_init() */
  DeRefMesh.Indi.AllocFacNdp = MIN_ALLOC_NUMBER;
  DeRefMesh.Indi.AllocFacNdt = MIN_ALLOC_NUMBER;
#endif


#ifdef VORONOI
  Mesh.Indi.AllocFacNdp = 1.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNdt = 8.0 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNvf = 8.0 * NumGas + MIN_ALLOC_NUMBER;
#if defined(VORONOI_DYNAMIC_UPDATE)
  Mesh.Indi.AllocFacNvc = 16.0 * NumGas + MIN_ALLOC_NUMBER;
  Nvc = 0;
#endif
  Mesh.Indi.AllocFacNinlist = 1.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacN_DP_Buffer = 0.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNflux = 0.01 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNradinflux = 0.01 * NumGas + MIN_ALLOC_NUMBER;
#endif

#ifdef AMR
  Mesh.Indi.AllocFacNdp = 1.2;
  Mesh.Indi.AllocFacNdt = 8.0 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNvf = 12;   //FIXME too much
  Mesh.Indi.AllocFacNnodes = 1.2;
  Mesh.Indi.AllocFacNgnodes = 0.2;

  //Mesh.Indi.AllocFacNinlist = 1.2 * NumGas + MIN_ALLOC_NUMBER;
  //Mesh.Indi.AllocFacN_DP_Buffer = 0.2 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNflux = 0.01 * NumGas + MIN_ALLOC_NUMBER;
  Mesh.Indi.AllocFacNradinflux = 0.01 * NumGas + MIN_ALLOC_NUMBER;

#if defined(AMR_CONNECTIONS)
  Mesh.Indi.AllocFacNvc = 16.0 * NumGas + MIN_ALLOC_NUMBER;
  Nvc = 0;
#endif
#endif

#ifdef MHD_POWELL
  for(j = 0; j < 3; j++)
    {
      All.Powell_Momentum[j] = 0;
      All.Powell_Angular_Momentum[j] = 0;
    }
  All.Powell_Energy = 0;
#endif

#ifdef NON_IDEAL_MHD
#if defined(OHMIC_DIFFUSION) || defined(IMPLICIT_OHMIC_DIFFUSION)
  All.OhmicDiffusionCoefficient /= All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam * All.UnitTime_in_s);
#endif
#ifdef AMBIPOLAR_DIFFUSION
  All.AmbipolarDiffusionCoefficient /= All.UnitLength_in_cm * All.UnitLength_in_cm / (All.HubbleParam * All.UnitTime_in_s);
#endif
#endif

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

#if defined(FOF) && (defined(BLACK_HOLES) || defined(GFM_WINDS_VARIABLE) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_WINDS_LOCAL))
  All.TimeNextOnTheFlyFoF = All.TimeBegin;
#endif

#if defined(GFM_WINDS_LOCAL) || defined(GFM_WINDS)
  All.CumWindEnergy_Should = 0.0;
  All.CumWindEnergy_Is = 0.0;
#endif

#if defined(BLACK_HOLES) && defined(USE_SFR)
  All.CumAGNEnergyEM_Is = 0.0;
  All.CumAGNEnergyEMobs_Is = 0.0;
  All.CumAGNEnergyM_Should = 0.0;
  All.CumAGNEnergyM_Is = 0.0;
  All.CumAGNEnergyT_Should = 0.0;
  All.CumAGNEnergyT_Is = 0.0;
#endif

#ifdef OTVET
  All.otvet_Radiation_Ti_begstep = 0;
#endif

#ifdef HEALRAY
  healray_init();
#endif


  set_softenings();

#ifdef ADAPTIVE_HYDRO_SOFTENING
  mpi_printf("INIT: Adaptive hydro softening, minimum gravitational softening for cells: %g\n", All.MinimumComovingHydroSoftening);
  mpi_printf("INIT: Adaptive hydro softening, maximum gravitational softening for cells: %g\n", All.MinimumComovingHydroSoftening * pow(All.AdaptiveHydroSofteningSpacing, NSOFTTYPES_HYDRO - 1));
  mpi_printf("INIT: Adaptive hydro softening, number of softening values: %d\n", NSOFTTYPES_HYDRO);
#endif
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  init_individual_softenings();
#endif

#ifdef SPECIAL_SOFTENINGS
  /* set special softenings for RG core and companion; companion ID must
   * be one larger than RG core ID; Softening types are 1 and 2 */
#ifndef ID_RGCORE
#error ID_RGCORE not defined!
#endif
  int firstID = ID_RGCORE;
  int secondID = ID_RGCORE + 1;
  int firstType = 1;
  int secondType = 2;
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].ID == firstID)
        {
          P[i].SofteningType = firstType;
          printf("INIT: Setting softening type of particle with ID %d and mass %g to %d\n", P[i].ID, P[i].Mass, P[i].SofteningType);
        }
      if(P[i].ID == secondID)
        {
          P[i].SofteningType = secondType;
          printf("INIT: Setting softening type of particle with ID %d and mass %g to %d\n", P[i].ID, P[i].Mass, P[i].SofteningType);
        }
    }
#endif

#ifdef SHIFT_BY_HALF_BOX
  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Pos[j] += 0.5 * All.BoxSize;
#endif

  for(i = 0; i < GRAVCOSTLEVELS; i++)
    All.LevelToTimeBin[i] = -1;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < GRAVCOSTLEVELS; j++)
      P[i].GravCost[j] = 0;

#ifdef TWODIMS
  for(i = 0; i < NumPart; i++)
    {
      P[i].Pos[2] = 0;
      P[i].Vel[2] = 0;
    }
#endif

  if(All.ComovingIntegrationOn) /*  change to new velocity variable */
    {
      for(i = 0; i < NumPart; i++)
        {
          for(j = 0; j < 3; j++)
            P[i].Vel[j] *= sqrt(All.Time) * All.Time;   /* for dm/gas particles, p = a^2 xdot */
        }
    }

  /* measure mean cell mass */
  int num = 0;
  long long glob_num;
  double glob_mass;
  mass = 0;

  for(i = 0; i < NumGas; i++)
#ifdef REFINEMENT_HIGH_RES_GAS
    if(SphP[i].AllowRefinement != 0)
#endif
      {
        num += 1;
        mass += P[i].Mass;
      }

  sumup_large_ints(1, &num, &glob_num);
  MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifndef REFINEMENT_HIGH_RES_GAS
  if(glob_num != All.TotNumGas)
    terminate("glob_num(=%lld) != All.TotNumGas(=%lld)", glob_num, All.TotNumGas);
#endif

  if(All.TotNumGas > 0 && (glob_num == 0 || glob_mass == 0))
    terminate("All.TotNumGas(=%lld) > 0 && (glob_num(=%lld) == 0 || glob_mass(=%g) == 0)", All.TotNumGas, glob_num, glob_mass);

  /* assign global variables that depend on the mean cell mass */
#if defined(REFINEMENT) || defined(GFM_STELLAR_EVOLUTION) || defined(BLACK_HOLES)
  if(All.ReferenceGasPartMass == 0)
    {
      if(!All.ComovingIntegrationOn)
        terminate("In non-comoving runs, ReferenceGasPartMass must be set to a non-zero value");

      All.ReferenceGasPartMass = glob_mass / glob_num;

      mpi_printf("REFINEMENT: The mean cell mass, which is used as a reference, is %g\n", All.ReferenceGasPartMass);
    }
  else
    mpi_printf("REFINEMENT: The given reference cell mass is %g\n", All.ReferenceGasPartMass);
#ifdef REFINEMENT
  All.TargetGasMass = All.TargetGasMassFactor * All.ReferenceGasPartMass;
  mpi_printf("REFINEMENT: setting All.TargetGasMass=%g\n", All.TargetGasMass);
#endif
#endif

#ifdef TRACER_MC
  All.MassTable[TRACER_MC] = (glob_mass / glob_num) / All.TracerMCPerCell;

  if(RestartFlag >= 2)
    {
      /* we have loaded tracers from a snapshot file, so verify counts */
      if(All.N_alltracer_global != get_total_number_of_tracers(-1))
        terminate("TRACER_MC: N_alltracer_global in init() after load differs! N_alltracer_global=%lld but current number of tracers is %lld", All.N_alltracer_global, get_total_number_of_tracers(-1));
    }
  else
    {
      /* we have generated tracers after IC start, so set global count */
      All.N_alltracer_global = get_total_number_of_tracers(-1);
    }

  mpi_printf("TRACER_MC: N_alltracer_global=%lld\n", All.N_alltracer_global);
#endif

#ifdef GENERATE_TRACER_MC_IN_ICS
  All.ReferenceTracerMCMass = (glob_mass / glob_num) / All.TracerMCPerCell;
#endif

  for(i = 0; i < TIMEBINS; i++)
    All.Ti_begstep[i] = 0;

  for(i = 0; i < NumPart; i++)  /*  start-up initialization */
    {
      for(j = 0; j < 3; j++)
        P[i].GravAccel[j] = 0;

#ifdef PMGRID
      for(j = 0; j < 3; j++)
        P[i].GravPM[j] = 0;
#endif
      P[i].TimeBinHydro = 0;
      P[i].TimeBinGrav = 0;
#ifdef SINKS
      P[i].TimeBinSink = 0;
#endif
#ifndef SECOND_ORDER_ICS
      P[i].OldAcc = 0;          /* Do not zero as masses are stored here */
#endif

#ifdef SELFGRAVITY
#ifdef EVALPOTENTIAL
      if(RestartFlag == 0)
        P[i].Potential = 0;
#endif
#endif

#ifdef STELLARAGE
      // For hdf5 ICs, these are read or set to zero in read_ic
      if(RestartFlag == 0 && All.ICFormat != 3)
        P[i].StellarAge = 0;
#endif

#if defined(METALS) && !defined(ADDBACKGROUNDGRID)
      if(RestartFlag == 0 && All.ICFormat != 3)
        P[i].Metallicity = 0;
#endif

#ifdef USE_SFR
      if(RestartFlag == 0 && P[i].Type == 0)
        SphP[i].Sfr = 0;
#endif

#ifdef REFINEMENT_AROUND_BH
      if(P[i].Type == 0)
        {
          SphP[i].RefBHFlag = 0;
#ifdef OUTPUT_REFBHCOUNTER
          SphP[i].RefBHCounter = 0;
#endif
        }
#endif

#ifdef SINKS
      if(P[i].Type == 0)
        SphP[i].InAccrRadius = 0;
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
        {
          count_holes++;

          if(RestartFlag == 0)
            {
              BPP(i).BH_Mass = All.SeedBlackHoleMass;

              if(P[i].Mass == 0)
                P[i].Mass = All.ReferenceGasPartMass;

#ifdef BH_FRICTION
              BPP(i).BH_MinPotTime = -1;
              BPP(i).BH_MinPotTime_Previous = -1;
              BPP(i).BH_MinPotCumAvgTime = 0;
#endif
              if(BPP(i).BH_DtGasNeighbor == 0)
                BPP(i).BH_DtGasNeighbor = 1.0;
#ifdef DRAINGAS
#if (DRAINGAS > 1)
              BPP(i).DrainBucketMass = 0.0;
#endif
#endif
#ifdef BH_BUBBLES
              BPP(i).BH_Mass_bubbles = All.SeedBlackHoleMass;
              BPP(i).BH_Mass_ini = All.SeedBlackHoleMass;
#endif
#if defined(GFM_AGN_RADIATION) || defined(MASSIVE_SEEDS_MERGER)
              BPP(i).HostHaloMass = 0.0;
#endif
#ifdef BH_THERMALFEEDBACK_ACC
              BPP(i).BH_AccEnergy = 0.0;
              BPP(i).BH_AccTime = 0.0;
#endif

#ifdef BH_SPIN_EVOLUTION
              double Phi, CosTheta, SinTheta;

              BPP(i).BH_SpinParameter = All.BHInitialSpin;
              BPP(i).BlackHoleRadiativeEfficiency = All.BlackHoleRadiativeEfficiency;
              BPP(i).BH_FlagOngAccEpis = 0;

              Phi = 2*M_PI*get_random_number();
              CosTheta = 1-2*get_random_number();
              SinTheta = sqrt(1-CosTheta*CosTheta);
              BPP(i).BH_SpinOrientation[0] = SinTheta*cos(Phi);
              BPP(i).BH_SpinOrientation[1] = SinTheta*sin(Phi);
              BPP(i).BH_SpinOrientation[2] = CosTheta;
#endif
            }
        }
#endif

#if defined(GFM_WINDS_VARIABLE) || defined(GFM_WINDS_LOCAL)
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].w.HostHaloMass = 0;     /* union! */
#endif

#ifdef GFM_STELLAR_EVOLUTION
      if(P[i].Type == 4 && RestartFlag == 0)
        {
          StarP[P[i].AuxDataID].SNIaRate = 0;
          StarP[P[i].AuxDataID].SNIIRate = 0;
        }
#ifdef GFM_CHEMTAGS
      for(int j = 0; j < GFM_N_CHEM_TAGS; j++)
        StarP[P[i].AuxDataID].MassMetalsChemTags[j] = 0.0;
#endif

      if(P[i].Type == 0 && RestartFlag == 0)
        {
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            SphP[i].MetalsFraction[j] = mass_fractions[j];
          SphP[i].Metallicity = GFM_INITIAL_METALLICITY;
#ifdef GFM_DUST
          int k;
          for(j = 0; j < GFM_DUST_N_CHANNELS; j++)
            {
              for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                {
                  SphP[i].MetalsDustFraction[j][k] = 0.0;
                }
            }
          SphP[i].NumSNII = 0.0;
#ifdef GFM_DUST_CAP
          SphP[i].DustMassCap = 0.0;
#endif
#endif
        }
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
      if(P[i].Type == 0 && RestartFlag == 0)
        {
          SphP[i].TotEgyFeed = 0.0;
          SphP[i].IntEgyFeed = 0.0;
          SphP[i].KinEgyFeed = 0.0;
        }
#endif

#ifdef DELAYED_COOLING
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].CoolShutoffTime = 0;

      if(P[i].Type == 4 && RestartFlag == 0)
        StarP[P[i].AuxDataID].BlastRadius = 0;
#endif

#ifdef DELAYED_COOLING_TURB
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].Uturb = 0;
#endif

#ifdef OUTPUT_STELLAR_FEEDBACK
      if(P[i].Type == 4 && RestartFlag == 0)
        {
          StarP[P[i].AuxDataID].SNII_Num = 0;
          StarP[P[i].AuxDataID].FeedbackEnergy = 0;
          StarP[P[i].AuxDataID].FeedbackMomentum = 0;
          StarP[P[i].AuxDataID].TotalMassReleased = 0;
          StarP[P[i].AuxDataID].TotalMassToEnrich = 0;
          StarP[P[i].AuxDataID].Cum_SNII_Num = 0;
          StarP[P[i].AuxDataID].Cum_FeedbackMomentum = 0;
        }
#endif

#ifdef OUTPUT_MOLECULAR_FRACTION
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].MolecularFrac = 0;
#endif

#ifdef OUTPUT_OPTICAL_DEPTH
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].OpticalDepth = 0;
#endif

#ifdef FM_RADIATION_FEEDBACK
      if(P[i].Type == 0 && RestartFlag == 0)
        SphP[i].GasRadCoolShutoffTime = 0;

      if(P[i].Type == 4 && RestartFlag == 0)
        {
          StarP[P[i].AuxDataID].StromgrenRadius = 0;
          StarP[P[i].AuxDataID].RadFeedTau = 0;
          StarP[P[i].AuxDataID].RadFeed_NumNgb = 0;
          StarP[P[i].AuxDataID].RadFeed_Flag = 0;
#ifdef FM_RADIATION_FEEDBACK_DEBUG
          StarP[P[i].AuxDataID].RadiationMomentumReleased = 0;
          StarP[P[i].AuxDataID].RadVelocityKick = 0;
          StarP[P[i].AuxDataID].NormSphRadFeedback = 0;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
          StarP[P[i].AuxDataID].NormSphRadFeedback_cold = 0;
#endif
          StarP[P[i].AuxDataID].RadCoolShutoffTime = 0;
#endif
        }

#endif

#ifdef FM_EARLY_STAR_FEEDBACK
      if(P[i].Type == 4 && RestartFlag == 0)
        StarP[P[i].AuxDataID].EarlyCumulativeFeedback = 0;
#endif
#ifdef OUTPUT_EARLY_STELLAR_FEEDBACK
      if(P[i].Type == 4 && RestartFlag == 0)
        {
          StarP[P[i].AuxDataID].EarlyFeedbackMomentum = 0;
          StarP[P[i].AuxDataID].Cum_EarlyFeedbackMomentum = 0;
        }
#endif
    }

#ifdef GFM_PREENRICH
  if(RestartFlag == 2)
    if(All.Time >= All.PreEnrichTime)
      {
        /* flag as done already */
        All.PreEnrichTime = -1;
      }
#endif

#ifdef SIDM
  sidm_SetGroundStateMass();
#endif

#if defined(BLACK_HOLES) && defined(FOF)
  set_DMmass();
#endif

#ifdef BLACK_HOLES
  MPI_Allreduce(&count_holes, &All.TotNumBHs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

#ifdef SINK_PARTICLES
  init_sink_particles(1);
#endif

  for(i = 0; i < TIMEBINS; i++)
    TimeBinSynchronized[i] = 1;

  reconstruct_timebins();

#ifdef PMGRID
  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif

#ifdef MONOTONE_CONDUCTION
  All.conduction_Ti_endstep = All.conduction_Ti_begstep = 0 ;
#endif

#ifdef TURBULENT_METALDIFFUSION
  All.metaldiff_Ti_endstep = All.metaldiff_Ti_begstep = 0;
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION
  All.ohmdiffusion_Ti_endstep = All.ohmdiffusion_Ti_begstep = 0 ;
#endif

  for(i = 0; i < NumGas; i++)   /* initialize sph_properties */
    {
      if(RestartFlag == 2 || RestartFlag == 3 || RestartFlag == 15)
        for(j = 0; j < 3; j++)
          SphP[i].Center[j] = P[i].Pos[j];

#if defined(CELL_CENTER_GRAVITY) && !defined(OUTPUT_CENTER_OF_MASS)
      if(RestartFlag == 18)
        for(j = 0; j < 3; j++)
          SphP[i].Center[j] = P[i].Pos[j];
#endif

      if(RestartFlag == 0)
        {
          for(j = 0; j < 3; j++)
            SphP[i].Center[j] = P[i].Pos[j];

          SphP[i].Hsml = 0;
#if defined(COOLING) && !defined(GRACKLE)
          SphP[i].Ne = 1.0;
#endif

#ifdef ATOMIC_DM
          SphP[i].Ne = 1.0;
#endif

#ifdef RADCOOL
          SphP[i].Phios = 0.0;
          SphP[i].Phins = 0.0;
#ifdef RADCOOL_HOTHALO
          SphP[i].PhiT6 = 0.0;
          SphP[i].PhiT7 = 0.0;
          SphP[i].PhiT8 = 0.0;
#endif
#endif

#ifdef TGCHEM
          if(TGCD.ChemMode >= 0)
            {
              SphP[i].Abund[0] = 0.;
              SphP[i].Abund[1] = TGCD.ChemInitAbH2;
              SphP[i].Abund[2] = TGCD.ChemInitAbHII;
            }
          else
            {
              for(j = 0; j < TGCHEM_NUM_ABUNDANCES; j++)
                SphP[i].Abund[j] = 0.;
            }

          SphP[i].Gamma = GAMMA;
#endif

#ifdef SGCHEM
#ifdef SGCHEM_OUTPUT_COOLTIME
          SphP[i].CoolTime = 0;
#endif
          if (All.SGChemConstInitAbundances) 
            {
#if CHEMISTRYNETWORK == 1
              SphP[i].TracAbund[IH2] = All.SGChemInitH2Abund;
              SphP[i].TracAbund[IHP] = All.SGChemInitHPAbund;
              SphP[i].TracAbund[IDP]   = All.SGChemInitDIIAbund;
              SphP[i].TracAbund[IHD]   = All.SGChemInitHDAbund;
              SphP[i].TracAbund[IHEP]  = All.SGChemInitHePAbund;
              SphP[i].TracAbund[IHEPP] = All.SGChemInitHeIIIAbund;
              SphP[i].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[i].TracAbund[IH2] - SphP[i].TracAbund[IHP] - SphP[i].TracAbund[IHD];
              SphP[i].TracAbund[IHEATOM] = ABHE - SphP[i].TracAbund[IHEP] - SphP[i].TracAbund[IHEPP];
              SphP[i].TracAbund[IDATOM] = All.DeutAbund - SphP[i].TracAbund[IDP] - SphP[i].TracAbund[IHD];
#endif
#if CHEMISTRYNETWORK == 4
              SphP[i].TracAbund[IH2] = All.SGChemInitH2Abund;
              SphP[i].TracAbund[IHP] = All.SGChemInitHPAbund;
              SphP[i].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[i].TracAbund[IH2] - SphP[i].TracAbund[IHP];
#endif
#if CHEMISTRYNETWORK == 5
              SphP[i].TracAbund[IH2] = All.SGChemInitH2Abund;
              SphP[i].TracAbund[IHP] = All.SGChemInitHPAbund;
              SphP[i].TracAbund[ICP] = All.SGChemInitCPAbund;
              SphP[i].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[i].TracAbund[IH2] - SphP[i].TracAbund[IHP];
#ifdef SGCHEM_VARIABLE_Z
              SphP[i].TracAbund[ICO] = SphP[i].CarbAbund - SphP[i].TracAbund[ICP];
#else
              SphP[i].TracAbund[ICO] = All.CarbAbund - SphP[i].TracAbund[ICP];
#endif
#endif /* CHEMISTRYNETWORK == 5 */
#if CHEMISTRYNETWORK == 15
              SphP[i].TracAbund[IH2] = All.SGChemInitH2Abund;
              SphP[i].TracAbund[IHP] = All.SGChemInitHPAbund;
              SphP[i].TracAbund[ICP] = All.SGChemInitCPAbund;
              SphP[i].TracAbund[ICHX]  = All.SGChemInitCHxAbund;
              SphP[i].TracAbund[IOHX]  = All.SGChemInitOHxAbund;
              SphP[i].TracAbund[ICO]   = All.SGChemInitCOAbund;
              SphP[i].TracAbund[IHCOP] = All.SGChemInitHCOPAbund;
              SphP[i].TracAbund[IHEP]  = All.SGChemInitHePAbund;
              SphP[i].TracAbund[IMP]   = All.SGChemInitMPAbund;
              SphP[i].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[i].TracAbund[IH2] - SphP[i].TracAbund[IHP] - SphP[i].TracAbund[ICHX] - SphP[i].TracAbund[IOHX] - SphP[i].TracAbund[IHCOP];
              SphP[i].TracAbund[IHEATOM] = HE_ABUND - SphP[i].TracAbund[IHEP];
#ifdef SGCHEM_VARIABLE_Z
              SphP[i].TracAbund[ICATOM] = SphP[i].CarbAbund - SphP[i].TracAbund[ICP] - SphP[i].TracAbund[ICHX] 
                                        - SphP[i].TracAbund[ICO] - SphP[i].TracAbund[IHCOP];
              SphP[i].TracAbund[IOATOM] = SphP[i].OxyAbund - SphP[i].TracAbund[IOHX] - SphP[i].TracAbund[ICO] - SphP[i].TracAbund[IHCOP];
              SphP[i].TracAbund[IMATOM] = SphP[i].MAbund - SphP[i].TracAbund[IMP];
#else
              SphP[i].TracAbund[ICATOM] = All.CarbAbund - SphP[i].TracAbund[ICP] - SphP[i].TracAbund[ICHX] - SphP[i].TracAbund[ICO] - SphP[i].TracAbund[IHCOP];
              SphP[i].TracAbund[IOATOM] = All.OxyAbund - SphP[i].TracAbund[IOHX] - SphP[i].TracAbund[ICO] - SphP[i].TracAbund[IHCOP];
              SphP[i].TracAbund[IMATOM] = All.MAbund - SphP[i].TracAbund[IMP];
#endif
#endif /* CHEMISTRYNETWORK == 15 */
	    }
          SphP[i].DustTemp = All.InitDustTemp;
#endif /* SGCHEM */
        }

#ifdef SGCHEM
      for(j = 0; j < SGCHEM_NUM_ADVECTED_SPECIES; j++)
        {
          SphP[i].MassTracAbund[j] = SphP[i].TracAbund[j] * P[i].Mass;
        }
#ifdef SGCHEM_VARIABLE_Z
      SphP[i].CarbMass       = SphP[i].CarbAbund * P[i].Mass;
      SphP[i].OxyMass        = SphP[i].OxyAbund  * P[i].Mass;
      SphP[i].MMass          = SphP[i].MAbund    * P[i].Mass;
      SphP[i].ZMass          = SphP[i].ZAtom     * P[i].Mass;
      SphP[i].ScaledDustMass = SphP[i].DustToGasRatio * P[i].Mass;
#endif

#endif

#ifdef FLD
      SphP[i].Lambda = 1. / 3.;
#endif

#ifdef MRT
      //      double nH_times_volume = HYDROGEN_MASSFRAC * P[i].Mass / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
      double nH_times_volume = P[i].Mass*SphP[i].Volume ;
      double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
      double nh_fac = 1.0 / (1.0 + y_fac) ;
      SphP[i].HI = 0.999999 ;
      SphP[i].HII = 0.000001 ;
      SphP[i].nHI = SphP[i].HI * nH_times_volume ;
      SphP[i].nHII = SphP[i].HII * nH_times_volume ;
#ifdef MRT_INCLUDE_HE
      SphP[i].HeI = 0.999998*y_fac ;
      SphP[i].HeII = 0.000001*y_fac ; 
      SphP[i].HeIII = 0.000001*y_fac ;
      SphP[i].Ne = SphP[i].HII + SphP[i].HeII + 2.0 * SphP[i].HeIII ;
      SphP[i].nHeI = SphP[i].HeI * nH_times_volume ;
      SphP[i].nHeII = SphP[i].HeII * nH_times_volume ;
      SphP[i].nHeIII = SphP[i].HeIII * nH_times_volume ;
#else
      SphP[i].Ne = SphP[i].HII ;
#endif
      SphP[i].ne = SphP[i].Ne * nH_times_volume ;
#endif

#ifdef MEASURE_DISSIPATION_RATE
      SphP[i].DuDt = 0;
#endif

#ifdef BH_THERMALFEEDBACK
      SphP[i].Injected_BH_Energy = 0;
#endif

#ifdef TGCHEM
      SphP[i].HydroHeatRate = 0.;
#endif

#if defined(TGCHEM) || defined(HEALRAY)
      SphP[i].EscFrac = -1.;
#endif

#ifdef HEALRAY
      SphP[i].HeatRate = 0.;
#endif

#ifdef AMR
      SphP[i].Level = All.MinRefLevel;  //FIXME a better estimate might be needed
#endif
    }

#ifndef NODEREFINE_BACKGROUND_GRID
  double mvol = 0;
  if(All.TotNumGas)
    {
#ifdef TWODIMS
      mvol = boxSize_X * boxSize_Y / All.TotNumGas;
#else
#ifdef ONEDIMS
      mvol = boxSize_X / All.TotNumGas;
#else
      mvol = boxSize_X * boxSize_Y * boxSize_Z / All.TotNumGas;
#endif
#endif
    }

  All.MeanVolume = mvol;
#endif

  mpi_printf("INIT: MeanVolume=%g\n", All.MeanVolume);

#if defined(BOUNDARY_STICKY_MINID) && defined(BOUNDARY_STICKY_MAXID) && defined(EXTERNALSHEARBOX)

  terminate("INIT: You are using an old EXTERNALSHEARBOX Config file.  Replace BOUNDARY_STICKY_MINID and BOUNDARY_STICKY_MAXID with STICKYFLAGS.");

#endif


#if defined(STICKYFLAGS)
  if(RestartFlag == 0){
    for(i = 0; i < NumPart; i++)
      {
        /* check whether particle is within stickly layer MaxDist
           if yes, check whether id is within sticky id range
           if no, change id
        */
        double boundary_dist;
        boundary_dist = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));

#if (REFLECTIVE_X == 2)
        boundary_dist = fmin(boundary_dist, fmin(boxSize_X - P[i].Pos[0], P[i].Pos[0]));
#endif
#if (REFLECTIVE_Y == 2)
        boundary_dist = fmin(boundary_dist, fmin(boxSize_Y - P[i].Pos[1], P[i].Pos[1]));
#endif
#if (REFLECTIVE_Z == 2)
        boundary_dist = fmin(boundary_dist, fmin(boxSize_Z - P[i].Pos[2], P[i].Pos[2]));
#endif

      if(boundary_dist < All.StickyLayerMaxDist)
        {
	  SphP[i].StickyFlag = 1;
	  
        }
      }
  }
#endif

#ifndef NO_ID_UNIQUE_CHECK
  test_id_uniqueness();
#ifdef TRACER_MC
  test_tracer_mc_id_uniqueness();
#endif
#endif

#ifdef INJECT_TRACER_INTO_SN
  calculate_max_tracer_mc_id();
#endif

#ifdef REFINEMENT_MERGE_CELLS
  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 0 && P[i].ID == 0)
      terminate("INIT: Cannot use ID==0 for gas in ICs with derefinement enabled.");
#endif

#if defined(VORONOI_DYNAMIC_UPDATE) || defined(AMR_CONNECTIONS)
  voronoi_init_connectivity(&Mesh);
#endif

#ifdef SUBBOX_SNAPSHOTS
  if(RestartFlag == 5 || RestartFlag == 10 || RestartFlag == 15)
    {
      for(i = 0; i < NumPart; i++)
        {
          P[i].Pos[0] -= SubboxXmin[All.SubboxNumber];
          P[i].Pos[1] -= SubboxYmin[All.SubboxNumber];
          P[i].Pos[2] -= SubboxZmin[All.SubboxNumber];
        }
      mpi_printf("SUBBOX_SNAPSHOTS: shifted for image generation -> lower left corner = (%g|%g|%g)\n", SubboxXmin[All.SubboxNumber], SubboxYmin[All.SubboxNumber], SubboxZmin[All.SubboxNumber]);
    }
#endif

#ifdef ADDBACKGROUNDGRID
  prepare_domain_backgroundgrid();
#endif


  domain_Decomposition();       /* do initial domain decomposition (gives equal numbers of particles) */

  if(RestartFlag == 18)  /* recalculation of potential */
    {
      mark_active_timebins();
      open_logfiles();
#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK) 
      sfr_init();
#endif
      set_non_standard_physics_for_current_time();

#ifdef PMGRID
      long_range_init_regionsize();
#endif

      compute_grav_accelerations(All.HighestActiveTimeBin, FLAG_FULL_TREE);

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) && defined(FOF)
      PS = (struct subfind_data *) mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));
      fof_prepare_output_order();   /* sort by type and Fileorder */
      fof_subfind_exchange(MPI_COMM_WORLD);
#endif

      sprintf(All.SnapshotFileBase, "%s_potupdated", All.SnapshotFileBase);
      mpi_printf("Start writing file %s\nRestartSnapNum %d\n", All.SnapshotFileBase, RestartSnapNum);
      savepositions(RestartSnapNum, 0);

      endrun();
    }

  if(RestartFlag == 12)
    {
      force_treeallocate(NumPart, All.MaxPart);
      force_treebuild(NumPart, 0, 0, 0);
      /*twopoint(); */
      myfree(Father);
      myfree(Nextnode);
#ifdef BLACK_HOLES
      myfree(Tree_AuxBH_Points);
#endif
      myfree(Tree_Points);
      force_treefree();
      endrun();
    }

  if(RestartFlag == 13)
    {
#ifdef PMGRID
      long_range_init_regionsize();
#ifndef GRAVITY_NOT_PERIODIC
      int n, n_type[NTYPES];
      long long ntot_type_all[NTYPES];
      /* determine global and local particle numbers */
      for(n = 0; n < NTYPES; n++)
        n_type[n] = 0;
      for(n = 0; n < NumPart; n++)
        n_type[P[n].Type]++;
      sumup_large_ints(NTYPES, n_type, ntot_type_all);

      calculate_power_spectra(RestartSnapNum, ntot_type_all);
#endif
#endif
      endrun();
    }

#ifdef AMR
#ifdef AMR_REMAP
  amr_amrtree = 0;
#endif
#endif
  /* will build tree */
  ngb_treeallocate();
  ngb_treebuild(NumGas);

  if(RestartFlag == 3)
    {
#ifdef FOF
      fof_fof(RestartSnapNum);
      DumpFlag = 1;
      savepositions(RestartSnapNum, 0);
#endif
      return (0);
    }

  if(RestartFlag == 7)
    {
#if defined(PMGRID) && defined(PERIODIC) && defined(VEL_POWERSPEC)
      powerspec_vel(RestartSnapNum);
#endif

#if defined(PERIODIC) && defined(ADJ_BOX_POWERSPEC)
      adj_box_powerspec();
#endif
      return (0);
    }

  All.Ti_Current = 0;

  if(RestartFlag == 0 || RestartFlag == 2 || RestartFlag == 4 || RestartFlag == 5 || RestartFlag == 8
     || RestartFlag == 9 || RestartFlag == 10 || RestartFlag == 11 || RestartFlag == 14 || RestartFlag == 15 || RestartFlag == 16 || RestartFlag == 17 || RestartFlag == 19)
    setup_smoothinglengths();

#ifdef RT_ADVECT
  All.Ti_LastRadTransfer = 0;
#endif


#ifdef BLACK_HOLES
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == 5)
          {
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
            BPP(i).BH_Hsml = All.SofteningTable[get_softening_type_from_mass(P[i].Mass)];
#else
            BPP(i).BH_Hsml = get_default_softening_of_particletype(P[i].Type);
#endif
            if(BPP(i).BH_DtGasNeighbor == 0)
              BPP(i).BH_DtGasNeighbor = 1.0;

          }
    }
#endif

#if defined(GFM_STELLAR_EVOLUTION) || defined(MRT_STARS)
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == 4)
          {
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
            StarP[P[i].AuxDataID].Hsml = All.SofteningTable[get_softening_type_from_mass(P[i].Mass)];
#else
            StarP[P[i].AuxDataID].Hsml = get_default_softening_of_particletype(4);
#endif
            if(All.StarformationOn == 0)
              {
                StarP[P[i].AuxDataID].BirthTime = All.Time;
 
                if(All.ComovingIntegrationOn == 0)
                  StarP[P[i].AuxDataID].BirthTime += 1e-15;
              }
          }
    }
#endif

#ifdef DUST_LIVE
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == DUST_LIVE)
          {
#ifdef INDIVIDUAL_GRAVITY_SOFTENING
            DustP[P[i].AuxDataID].Hsml = All.SofteningTable[get_softening_type_from_mass(P[i].Mass)];
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
            DustP[P[i].AuxDataID].DustHsml = All.SofteningTable[get_softening_type_from_mass(P[i].Mass)];
#endif
#else
            DustP[P[i].AuxDataID].Hsml = get_default_softening_of_particletype(DUST_LIVE);
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
            DustP[P[i].AuxDataID].DustHsml = get_default_softening_of_particletype(DUST_LIVE);
#endif
#endif
          }
    }
#endif

#if defined(OTVET_SCATTER_SOURCE) && !defined(GFM)
  if(RestartFlag == 0 || RestartFlag == 2)
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == 4)
          {
            P[i].OtvetHsml = get_default_softening_of_particletype(4);
            P[i].OtvetGasDensity = 0;
          }
    }
#endif

#ifdef SIDM
  sidm_Init_Particles();
#endif


#ifdef ADDBACKGROUNDGRID
  init_aux_fields();
  // This return more clearly shows that this function terminates the run
  return add_backgroundgrid();
#endif

#ifdef AMR
#ifdef AMR_REMAP
  init_aux_fields();
  return amr_remap();
#endif
#endif

  /* at this point, the entropy variable actually contains the
   * internal energy, read in from the initial conditions file.
   * Once the density has been computed, we can convert to entropy.
   */

#if defined (VS_TURB) || defined (AB_TURB)
  {
    double mass = 0, glob_mass;
    int i;
    for(i = 0; i < NumGas; i++)
      mass += P[i].Mass;
    MPI_Allreduce(&mass, &glob_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    All.RefDensity = glob_mass / pow(All.BoxSize, 3);
    All.RefEntropy = All.IsoSoundSpeed * All.IsoSoundSpeed / pow(All.RefDensity, GAMMA_MINUS1);
  }
#endif

  create_mesh();
  mesh_setup_exchange();

  if(RestartFlag == 14)
    {
#ifdef VORONOI
      char tess_name[1024];
      sprintf(tess_name, "%s/tess_%03d", All.OutputDir, RestartSnapNum);
      write_voronoi_mesh(&Mesh, tess_name, 0, NTask - 1);
#else
      terminate("termination in init(): set flag VORONOI in Config.sh.");
#endif
      return 0;
    }




  if(RestartFlag == 16)
    {
      if(Argc != 10)
        {
        terminate("termination in init(): 10 arguments needed")}


#if defined (TETRA_INDEX_IN_FACE) && !defined(TWODIMS) && !defined(ONEDIMS)

      double center[3];
      center[0] = atof(Argv[4]);
      center[1] = atof(Argv[5]);
      center[2] = atof(Argv[6]);

      double normal[3];

      normal[0] = atof(Argv[7]);
      normal[1] = atof(Argv[8]);
      normal[2] = atof(Argv[9]);

      char tess_name[1024];
      sprintf(tess_name, "%s/tess_section%03d", All.OutputDir, RestartSnapNum);

      intersection_plane_grid(center, normal, tess_name);
#else
      terminate("termination in init(): set flag TETRA_INDEX_IN_FACE and unset TWODIMS/ONEDIMS in Config.sh.");
#endif
      return 0;
    }

  if(RestartFlag == 11)
    {
#ifdef TRACER_PARTICLE
      update_primitive_variables();
      exchange_primitive_variables();
      calculate_gradients();
      exchange_primitive_variables_and_gradients();

#if (defined (VS_TURB) || defined (AB_TURB)) && defined(POWERSPEC_GRID)
      powersepc_turb_init();
      powerspec_turb(RestartSnapNum, TRACER_PARTICLE);
#endif

#if defined(VORONOI_IMAGES_FOREACHSNAPSHOT) || defined(VORONOI_FREQUENT_IMAGES)
      make_voronoi_image(RestartSnapNum);
      make_voronoi_image_slice(RestartSnapNum);
#endif

#if defined(TWODIMS) && defined(VORONOI_FIELD_DUMP_PIXELS_X)
      do_special_dump(RestartSnapNum, 1);
#ifdef VORONOI_NOGRADS
      do_special_dump(RestartSnapNum, 0);
#endif
#endif

      DumpFlag = 1;
      savepositions(RestartSnapNum, 0);
#endif
      return (0);
    }


#ifdef MHD_CT
#ifdef MHD_CT_IC
  set_A_ICs(MHD_CT_IC);
#else
  for(i = 0, mass = 0; i < NumGas; i++)
    {
      SphP[i].AConserved[0] = SphP[i].A[0] * SphP[i].Volume;
      SphP[i].AConserved[1] = SphP[i].A[1] * SphP[i].Volume;
      SphP[i].AConserved[2] = SphP[i].A[2] * SphP[i].Volume;
    }
  correct_ctr_b();
#endif
#endif



#ifdef MHD_SEEDPSPEC
  set_pscpec_ICs();
#endif

#if defined(DG) && defined(DG_TEST_PROBLEM)
  dg_set_initial_conditions();
#elif defined(DG_TEST_PROBLEM)
  arepo_set_initial_conditions();
#endif

#ifdef COSMIC_RAYS
  All.TotCountCRLimiter = 0.;
  All.TotalCREnergyInjected = 0.;
  All.TotalCREnergyCooledHadronic = 0.;
  All.TotalCREnergyCooledCoulomb = 0.;
  All.TotalCREnergyCooled = 0;
  
#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  All.TotalCREnergyErrorDiffusion = 0;
  All.TotalCREnergyChangeAdiabatic = 0;
  All.TotalCREnergyLossSfr = 0;
  All.TotalCREnergyUpdatePrims = 0;
#endif
#endif

  for(i = 0, mass = 0; i < NumGas; i++)
    {
      if(RestartFlag == 0)
        {
#ifdef MESHRELAX_DENSITY_IN_INPUT
          P[i].Mass *= SphP[i].Volume;
#endif
#ifdef TRACER_FIELD
          SphP[i].ConservedTracer *= P[i].Mass;
#endif
        }

#ifdef REFINEMENT_MERGE_PAIRS
      SphP[i].DerefPartnerId = 0;
      SphP[i].DerefPartnerIndex = -1;
#endif

#ifdef TGSET
      tgset_dens_init(i);
#endif

      SphP[i].Density = P[i].Mass / SphP[i].Volume;

      if(SphP[i].Density < All.MinimumDensityOnStartUp)
        {
          SphP[i].Density = All.MinimumDensityOnStartUp;

          P[i].Mass = SphP[i].Volume * SphP[i].Density;
        }

      SphP[i].Momentum[0] = P[i].Mass * P[i].Vel[0];
      SphP[i].Momentum[1] = P[i].Mass * P[i].Vel[1];
      SphP[i].Momentum[2] = P[i].Mass * P[i].Vel[2];


#ifdef MRT
      double X = P[i].Pos[0] - All.BoxSize/2.0 ;
      double Y = P[i].Pos[1] - All.BoxSize/2.0 ;
      double Z = P[i].Pos[2] ;

      double R = sqrt(X*X + Y*Y + Z*Z) ;

      if(R==0.0)
	R = 1e-8 ;

      //      printf("R = %g \n", R) ;
      
      for(int num1=0;num1<MRT_BINS;num1++)
	{
	  //if(R<0.2)
	  //{
	  SphP[i].DensPhot[num1] = 40.0*exp(-(R-1.2)*(R-1.2)/2.0/0.01) ;
	      //      printf("Entered Here\n") ;

	      //}
	      //else
	      //{
	  //SphP[i].DensPhot[num1] = MINDENSPHOT ; 
	      //printf("Only entered here\n") ;
	      //}
	  SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume ;
	  if(isnan(SphP[i].Cons_DensPhot[num1]))
	    terminate("Something went wrong\n") ;
	  
	  for(j=0;j<3;j++)
	    {
	      SphP[i].RT_F[num1][j] = MINDENSPHOT;
	      
	      SphP[i].Cons_RT_F[num1][j] = SphP[i].RT_F[num1][j] * SphP[i].Volume ;
	    }
	  SphP[i].RT_F[num1][0] = 0.99 * SphP[i].DensPhot[num1] * X/R;
	  SphP[i].RT_F[num1][1] = 0.99 * SphP[i].DensPhot[num1] *Y/R;
	  SphP[i].RT_F[num1][2] = 0.99 * SphP[i].DensPhot[num1] *Z/R;

	  SphP[i].Cons_RT_F[num1][0] = SphP[i].RT_F[num1][0] * SphP[i].Volume ;
	  SphP[i].Cons_RT_F[num1][1] = SphP[i].RT_F[num1][1] * SphP[i].Volume ;
	  SphP[i].Cons_RT_F[num1][2] = SphP[i].RT_F[num1][2] * SphP[i].Volume ;
	}

#ifdef MRT_COMOVING
      SphP[i].Old_Vel[0] = 0.6 ;
      SphP[i].Old_Vel[0] = 0.0 ;
      SphP[i].Old_Vel[0] = 0.0 ;
#endif

      if(i==0)
	mpi_printf("RT: INIT: Initialized conservative quantities\n") ;
#endif




#ifdef MHD
#ifndef MHD_CT                  // ICs for B, BConserved are already set for CT by this time in set_A_ICs()
      /* The initial magnetic field is assumed to be the comoving magnetic field ( B_c = B_phys * a^2 ) for cosmological simulations */
#ifdef MHD_SEEDFIELD
if (RestartFlag == 0)
{
      if(i == 0)
        {
          mpi_printf("MHD Seed field=%g, direction=%d\n", All.B_value, All.B_dir);
        }

      int k;
      double bfac = 1. / (sqrt(All.UnitMass_in_g / All.UnitLength_in_cm) / (All.UnitTime_in_s / All.HubbleParam));

#if defined(MONOTONE_CONDUCTION) && defined(THERMAL_INSTABILITY_TEST)   /*&& defined(DO_NOT_MOVE_GAS) */
      /*if(P[i].Pos[0] < All.BoxSize/2.0)
         {
         SphP[i].BConserved[0] = -1.0 ;
         SphP[i].BConserved[1] = 1.0 ;
         }
         else
         {
         SphP[i].BConserved[0] = 1.0 ;
         SphP[i].BConserved[1] = 1.0 ;
         } */
      /*

         double  mcx = P[i].Pos[0] - All.BoxSize/2.0 ;
         double  mcy = P[i].Pos[1] - All.BoxSize/2.0 ;

         double   mcr = sqrt(mcx*mcx + mcy*mcy) ;



         double b = 1e-6;

       */
/*
e    if( mcr < 5 )
      b = 0.0 ;
    else if (mcr < 15)
      b = 10. * (mcr - 5);
    else if (mcr < 35)
      b = 1.;
    else if (mcr < 45)
      b = 10. * (45 - mcr);
    else
    b = 0;*/


      //double b = 1.0 ;
      /*      if(mcr==0.0)
         mcr = 1.0 ;

         double sintheta = b*mod(mcy/mcr) ;
         double costheta = b*mod(mcx/mcr) ;


         //    double modtheta = sqrt(sintheta*sintheta + costheta*costheta) ;

         //if(modtheta == 0.0)
         //modtheta = 1.0 ;


         //sintheta /= modtheta ;
         //costheta /= modtheta ;

         if((sintheta==0.0) && (costheta ==0.0))
         {
         sintheta = 1.0 ;
         costheta = 0.0 ;
         }

         if(mcx>=0.0 && mcy>=0.0)
         {
         SphP[i].BConserved[0] = sintheta ;
         SphP[i].BConserved[1] = -costheta ;
         }
         else if(mcx>0.0 && mcy<0.0)
         {
         SphP[i].BConserved[0] = -sintheta ;
         SphP[i].BConserved[1] = -costheta ;
         }
         else if(mcx<0.0 && mcy<0.0)
         {
         SphP[i].BConserved[0] = -sintheta ;
         SphP[i].BConserved[1] = costheta ;
         }
         else
         {
         SphP[i].BConserved[0] = sintheta ;
         SphP[i].BConserved[1] = costheta ;
         }
       */
      /*if(P[i].Pos[0]<All.BoxSize/2.0)
         {
         SphP[i].BConserved[0] = -1.0 ;
         SphP[i].BConserved[1] = 1.0 ;
         }
         else
         {
         SphP[i].BConserved[0] = 1.0 ;
         SphP[i].BConserved[1] = 1.0 ;
         } */


      //SphP[i].BConserved[2] = SphP[i].BConserved[0]  ;
      //SphP[i].BConserved[0] = 0.0 ;

      double b = 1e-9;

      SphP[i].BConserved[2] = b * 0.0;
      SphP[i].BConserved[1] = b * 1.0;
      SphP[i].BConserved[0] = b * 1.0;

      for(k = 0; k < 3; k++)
        SphP[i].B[k] = SphP[i].BConserved[k] / SphP[i].Volume;


#else
#ifdef EXTERNALSHEARBOX
#ifdef EQUIPARTITION
      double eunit=All.UnitLength_in_cm*All.UnitLength_in_cm/All.UnitTime_in_s/All.UnitTime_in_s;
      double dunit=All.UnitMass_in_g/All.UnitLength_in_cm/All.UnitLength_in_cm/All.UnitLength_in_cm;
      double B_value = sqrt(EQUIPARTITION * SphP[i].Utherm*eunit*SphP[i].Density*dunit * 8.0 * M_PI * GAMMA_MINUS1);
      SphP[i].Utherm /= (1.+EQUIPARTITION);
#else
      double bscale = 61.0 * (All.ShearBoxFg / 0.1 / All.ShearBoxMu) / (All.ShearBoxSigma0 / 10.0);
      double midplane_density = All.ShearBoxSigma0 / 2.0 / bscale;
      double B_value = All.B_value * pow(SphP[i].Density / midplane_density, 2.0 / 3.0);
#endif

#else
      double B_value = All.B_value;
#endif
      for(k = 0; k < 3; k++)
        if(All.B_dir & (1 << k))
          {
            SphP[i].BConserved[k] = B_value * SphP[i].Volume * bfac;
            SphP[i].B[k] = SphP[i].BConserved[k] / SphP[i].Volume;
          }
        else
          {
            SphP[i].BConserved[k] = 0;
            SphP[i].B[k] = SphP[i].BConserved[k] / SphP[i].Volume;
          }

#endif
      if(i == 0)
        {
          mpi_printf("BConserved[0] = %g|%g|%g\n", SphP[i].BConserved[0], SphP[i].BConserved[1], SphP[i].BConserved[2]);
          mpi_printf("Volume[0] %g bfac %g\n", SphP[i].Volume, bfac);
        }
      /* convert Gauss-cgs to heavyside - lorentz */
      {
        int kk;
        for(kk = 0; kk < 3; kk++)
          {
            SphP[i].BConserved[kk] /= sqrt(4. * M_PI);
            SphP[i].B[kk] /= sqrt(4. * M_PI);
          }
      }
}
else
{
      SphP[i].BConserved[0] = SphP[i].B[0] * SphP[i].Volume;
      SphP[i].BConserved[1] = SphP[i].B[1] * SphP[i].Volume;
      SphP[i].BConserved[2] = SphP[i].B[2] * SphP[i].Volume;

}
#else
      SphP[i].BConserved[0] = SphP[i].B[0] * SphP[i].Volume;
      SphP[i].BConserved[1] = SphP[i].B[1] * SphP[i].Volume;
      SphP[i].BConserved[2] = SphP[i].B[2] * SphP[i].Volume;

#endif

#endif
#endif

#ifdef EOS_DEGENERATE
#ifdef EOS_NSPECIES
#endif
#ifdef EOS_INITTEMP
      {
        if(SphP[i].Density > 1e-4)
          {
            struct eos_result res;
            SphP[i].EOSTemperature = EOS_INITTEMP;
            eos_calc_tgiven(SphP[i].Density, SphP[i].Composition, SphP[i].EOSTemperature, &res);
            SphP[i].Utherm = res.e.v;
          }
      }
#else
      SphP[i].EOSTemperature = -1;
#endif
#endif

#ifdef EOS_OPAL
      SphP[i].EOSTemperature = -1.0;
#endif

#ifdef COSMIC_RAYS
#ifdef COSMIC_RAYS_IN_ICS
      SphP[i].CR_Energy = P[i].Mass * SphP[i].CR_SpecificEnergy;
#else
#ifdef EQUIPARTITION
      //      double Utherm_0 = SphP[i].Utherm/(1.0-EQUIPARTITION);
      // SphP[i].CR_Energy = EQUIPARTITION * Utherm_0 * P[i].Mass * (GAMMA_MINUS1/(All.GammaCR-1.0));
      //SphP[i].Utherm = Utherm_0 * (1.0 - 2.0*EQUIPARTITION);

      //SphP[i].CR_Energy = (GAMMA_MINUS1*EQUIPARTITION*SphP[i].Utherm*midplane_density*P[i].Mass)/((All.GammaCR - 1.0) * SphP[i].Density);
      SphP[i].CR_Energy = 0.;
#else
      SphP[i].CR_Energy = 0.;
#endif
#endif
      if(SphP[i].CR_Energy < All.MinimumCREnergyDensity * SphP[i].Volume)
        SphP[i].CR_Energy = All.MinimumCREnergyDensity * SphP[i].Volume;
#endif /* COSMIC_RAYS */

      /* utherm has been loaded from IC file */
#ifdef MESHRELAX
      SphP[i].Energy = P[i].Mass * SphP[i].Utherm;
#else

#if defined (VS_TURB) || defined (AB_TURB)
      SphP[i].Utherm = All.RefEntropy * pow(SphP[i].Density, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif

      SphP[i].Energy = P[i].Mass * All.cf_atime * All.cf_atime * SphP[i].Utherm + 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
#endif

#if defined (VS_TURB) || defined (AB_TURB)
      SphP[i].Energy += P[i].Mass * get_turb_pot(P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
#endif

#ifdef MHD
      SphP[i].Energy += 0.5 * (SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) * SphP[i].Volume * All.cf_atime;
#ifdef MHD_THERMAL_ENERGY_SWITCH
      SphP[i].Etherm = SphP[i].Utherm * P[i].Mass;
#endif
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
#ifdef TGCHEM
      SphP[i].A = (SphP[i].Gamma - 1) * SphP[i].Utherm / pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
#else
      SphP[i].A = GAMMA_MINUS1 * SphP[i].Utherm / pow(SphP[i].Density * All.cf_a3inv, GAMMA_MINUS1);
#endif
      SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;

#endif

      for(j = 0; j < 3; j++)
        SphP[i].VelVertex[j] = P[i].Vel[j];

#ifdef OUTPUT_CELL_SPIN
      for(j = 0; j < 3; j++)
        {
          SphP[i].Spin[j] = 0;
          SphP[i].CenterOld[j] = SphP[i].Center[j];
          SphP[i].CenterOffset[j] = 0;
          SphP[i].CenterOffsetMass[j] = 0;
        }
#endif

#ifdef TRACER_FIELD
      if(P[i].Mass > 0)
        SphP[i].Tracer = SphP[i].ConservedTracer / P[i].Mass;
      else
        SphP[i].Tracer = 0;
#endif
      mass += P[i].Mass;
    }


#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
  for(i = 0; i < NumGas; i++)
    {
      for(j = 0; j < EOS_NSPECIES; j++)
        SphP[i].MassComposition[j] = SphP[i].Composition[j] * P[i].Mass;
      SphP[i].EOSTemperature = -1.0;       /* initial temperature for EOS */
    }
#endif

#ifdef PASSIVE_SCALARS
  for(i = 0; i < NumGas; i++)
    {
      for(j = 0; j < PASSIVE_SCALARS; j++)
        SphP[i].PConservedScalars[j] = SphP[i].PScalars[j] * P[i].Mass;
    }

#endif

#ifdef NUCLEAR_NETWORK
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].dedt = 0;
    }
#endif

#ifdef METALS
  /* initialize absolute masses in materials */
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Metallicity = P[i].Metallicity;   // set above

#ifdef MIN_METALLICITY_ON_STARTUP
      if(SphP[i].Metallicity < All.MinimumMetallicityOnStartUp)
        {
          SphP[i].Metallicity = All.MinimumMetallicityOnStartUp;
          P[i].Metallicity = All.MinimumMetallicityOnStartUp;
        }
#endif

      SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
    }
#endif


#ifdef CLOUD_ENRICH_METALS
  mpi_printf("INIT:CLOUD_ENRICH_METALS entered\n");
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].Metallicity = 0.3 * GFM_SOLAR_METALLICITY;
      //printf("Metallicity = %g\t solar met = %g\t norm met = %g\n", SphP[i].Metallicity, GFM_SOLAR_METALLICITY, SphP[i].Metallicity/GFM_SOLAR_METALLICITY) ;
    }
#endif



#if defined(GFM_STELLAR_EVOLUTION) || defined(TURBULENT_METALDIFFUSION)
  /* initialize absolute masses in metals */
  /* note: for GFM_STELLAR_EVOLUTION==1 MassMetallicity and MassMetals are inconsistent with */
  /* Metallicity and MetalsFraction, so their values are not written to the snapshots, */
  /* therefore if RestartFlag==2 here, MassMetallicity and MassMetals will not have their old values */
  for(i = 0; i < NumGas; i++)
    {
      SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        SphP[i].MassMetals[j] = SphP[i].MetalsFraction[j] * P[i].Mass;
#ifdef GFM_DUST
      int k;
      for(j = 0; j < GFM_DUST_N_CHANNELS; j++)
        {
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            {
              SphP[i].MassMetalsDust[j][k] = SphP[i].MetalsDustFraction[j][k] * P[i].Mass;
            }
        }
#endif
    }
#endif

#ifdef SPECIAL_RELATIVITY
  for(i = 0; i < NumGas; i++)
    {
      compute_conserved_quantities_from_ICs_special_relativity(i);
    }
#endif

#if defined(GRACKLE) && !defined(GRACKLE_TAB)
  if((RestartFlag == 0) || (RestartFlag == 2) || (RestartFlag == 3))
  {
#ifndef GRACKLE_ABUNDANCE_IN_ICS
    grackle_initialise_abundances();
#endif
    grackle_converge_abundances();
  }
#endif

#ifdef GENERAL_RELATIVITY
#if ATMOSPHERE_GENERAL_RELATIVITY == 1

 general_relativity_determine_maximum_density();
#endif

 // we should include a first atmosphere treatment here !
#ifdef MESHRELAX_DENSITY_IN_INPUT
    for(i = 0; i < NumGas; i++)
    {
#if ATMOSPHERE_GENERAL_RELATIVITY == 1
      compute_conserved_quantities_from_ICs_general_relativity_apply_atmosphere(i);
#else
      compute_conserved_quantities_from_ICs_general_relativity(i);
#endif
    }
#else
  for(i = 0; i < NumGas; i++)
    {
      compute_conserved_quantities_from_ICs_general_relativity(i);
    }
  //terminate("with this prepro flag no atmosphere treatment possible\n");
#endif
#endif

  if(RestartFlag == 17)
    {
      update_primitive_variables();
      exchange_primitive_variables();
      calculate_gradients();
      exchange_primitive_variables_and_gradients();
      DumpFlag = 1;
      savepositions(RestartSnapNum + 1, 0);
      return (0);
    }


 if(RestartFlag == 19)
    {

#ifdef CALCULATE_QUANTITIES_IN_POSTPROCESS

      update_primitive_variables();
      exchange_primitive_variables();
      calculate_gradients();
      exchange_primitive_variables_and_gradients();


      calculate_quantities_in_postprocess();
#endif

#if defined(FOF)
      PS = (struct subfind_data *) mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));
      fof_prepare_output_order();   /* sort by type and Fileorder */
      fof_subfind_exchange(MPI_COMM_WORLD);
#endif

#ifdef CALCULATE_QUANTITIES_IN_POSTPROCESS
      save_additional_quantities_data(RestartSnapNum);
#endif
      return (0) ;
    }

#ifdef NUCLEAR_NETWORK_DETONATE
  detonate();
#endif

  update_primitive_variables();

#ifdef TREE_BASED_TIMESTEPS
  tree_based_timesteps_setsoundspeeds();
#endif


  /* initialize star formation rate */
#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK) //&& !defined(QUICK_LYALPHA)
  sfr_init();
#endif
  if(RestartFlag != 15) // for the shock finder the sfr is read from the snapshot file
    {
#if defined(USE_SFR) && !defined(LOCAL_FEEDBACK)
      for(i = 0; i < NumGas; i++)
        SphP[i].Sfr = get_starformation_rate(i);
#endif
    }


#ifdef BECDM
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != 1)
        continue;
      P[i].PsiRe = sqrt( P[i].Mass / pow(All.BoxSize/((double) PMGRID), 3) );
      P[i].PsiIm = 0.0;
#ifdef BECDM_INPUT_PHASE_AS_VX
      double fac = 1;
      if(All.ComovingIntegrationOn) 
        fac = sqrt(All.Time) * All.Time;
      if (RestartFlag == 0) {
        P[i].PsiRe = sqrt( P[i].Mass / pow(All.BoxSize/((double) PMGRID), 3) ) * cos(P[i].Vel[0]/fac);
        P[i].PsiIm = sqrt( P[i].Mass / pow(All.BoxSize/((double) PMGRID), 3) ) * sin(P[i].Vel[0]/fac);
      }
#endif   
    P[i].Vel[0] = 0.0;   
    P[i].Vel[1] = 0.0;   
    P[i].Vel[2] = 0.0;   
    }
#endif


#ifdef TRACER_PARTICLE
  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type != TRACER_PARTICLE)
        continue;
      P[i].TracerHsml = pow(mvol, 1.0 / NUMDIMS);
    }

  All.MassTable[TRACER_PARTICLE] = 1.0; /* do not output tracer masses in snapshots, nor expect them in ICs */
#endif

  update_primitive_variables();

  exchange_primitive_variables();

  calculate_gradients();

  exchange_primitive_variables_and_gradients();

#ifdef RT_ADVECT
  mpi_printf("rt_advect grad init\n");
  rt_voronoi_exchange_primitive_variables();
  rt_calculate_green_gauss_gradients();
  rt_voronoi_exchange_primitive_variables_and_gradients();
#endif

#ifdef SECOND_DERIVATIVES
  calculate_green_gauss_hessian(&Mesh);
  voronoi_exchange_hessians();
#endif

#ifdef TRACER_TRAJECTORY
  tracer_init();
  reconstruct_timebins();
#endif

#if !defined(ONEDIMS) && !defined(TWODIMS)
  int xaxis, yaxis, zaxis, weight_flag = 0;
  double xmin, xmax, ymin, ymax, zmin, zmax;
#endif

  if(RestartFlag == 4)
    {
#ifdef VORONOI
      int pixels_x, pixels_y;
#ifdef TWODIMS
      if(Argc == 6)
        {
          pixels_x = atoi(Argv[4]);
          pixels_y = atoi(Argv[5]);
          make_2d_voronoi_image(RestartSnapNum, pixels_x, pixels_y);
        }
      if(Argc > 6)
        {
          double xmin, xmax, ymin, ymax;
          pixels_x = atoi(Argv[4]);
          pixels_y = atoi(Argv[5]);
          xmin = atof(Argv[9]);
          xmax = atof(Argv[10]);
          ymin = atof(Argv[11]);
          ymax = atof(Argv[12]);

          make_2d_voronoi_image_zoomed(&Mesh, RestartSnapNum, pixels_x, pixels_y, xmin, xmax, ymin, ymax);
        }
#elif !defined(ONEDIMS) && !defined(TWODIMS)
      /*   3D */

      if(!get_image_limits(Argc, Argv, RestartFlag, &pixels_x, &pixels_y, &xaxis, &yaxis, &zaxis, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &weight_flag))
        return 0;

      // zmin is zval in this case
      make_3d_voronoi_slice_image(RestartSnapNum, 1, pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin);
#ifdef VORONOI_NOGRADS
      make_3d_voronoi_slice_image(RestartSnapNum, 0, pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin);
#endif
      make_3d_voronoi_listfaces(&Mesh, RestartSnapNum, xaxis, yaxis, zaxis, zmin);
#endif
      mpi_printf("\ndensity slice generated. Stopping.\n");
      return (0);
#endif /* VORONOI */
    }

  if(RestartFlag == 5)
    {
#ifdef VORONOI
#if !defined(ONEDIMS) && !defined(TWODIMS)
      int pixels_x, pixels_y;

      if(!get_image_limits(Argc, Argv, RestartFlag, &pixels_x, &pixels_y, &xaxis, &yaxis, &zaxis, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &weight_flag))
        return 0;

      make_3d_voronoi_projected_image(RestartSnapNum, 1, pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax, weight_flag);
#ifdef VORONOI_NOGRADS
      make_3d_voronoi_projected_image(RestartSnapNum, 0, pixels_x, pixels_y, xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax, weight_flag);
#endif
#ifdef SINKS
      write_sink_data(xaxis, yaxis, zaxis, xmin, xmax, ymin, ymax, zmin, zmax);
#endif
#endif
      mpi_printf("\ndensity projection generated. Stopping.\n");

      return (0);
#endif /* VORONOI */
    }

  if(RestartFlag == 8)
    {
#ifdef VORONOI
#if !defined(ONEDIMS) && !defined(TWODIMS)
      /*   3D */
      int pixels_x, pixels_y;
      if(!get_image_limits(Argc, Argv, RestartFlag, &pixels_x, &pixels_y, &xaxis, &yaxis, &zaxis, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, &weight_flag))
        return 0;

      make_3d_voronoi_grid(RestartSnapNum, pixels_x, pixels_y, xaxis, xmin, xmax, ymin, ymax, zmin, zmax);

#endif
      mpi_printf("\ndensity grid generated. Stopping.\n");

      return (0);
#endif /* VORONOI */
    }

  if(RestartFlag == 9)
    {
#if defined(VORONOI) && defined(VORONOI_PROJ)
#if !defined(ONEDIMS) && !defined(TWODIMS)
      struct projection_data pdata;
      if(!voronoi_proj_setup(Argc, Argv, &pdata))
        return 0;

      voronoi_proj_init(&pdata, 0);
      voronoi_proj_run(&pdata);
      voronoi_proj_deinit(&pdata);

      voronoi_proj_init(&pdata, 1);
      voronoi_proj_run(&pdata);
      voronoi_proj_deinit(&pdata);
#endif
      mpi_printf("\nProjections done. Stopping.\n");

      return (0);
#endif /* VORONOI */
    }

#if defined(VORONOI) && defined(DVR_RENDER) && !defined(ONEDIMS) && !defined(TWODIMS)
#if (DVR_RENDER==1)
  All.DvrFrameFileNum = 0;
#endif

  if(RestartFlag == 10)
    {
      if(Argc < 7)
        {
          printf("DVR_RENDER: Not enough parameters, aborting.\n");
          return 0;
        }

      dvr_render_loop(Argv[6], atoi(Argv[4]), atoi(Argv[5]), 0);
      mpi_printf("\nDVR_RENDER: Rendering done.\n");
      return (0);
    }
#endif

  if(RestartFlag == 15)
    {
      open_logfiles();

#ifdef SHOCK_FINDER_POST_PROCESSING
      All.SnapshotFileCount = RestartSnapNum;
      shock_finder();
      write_cpu_log();
      close_logfiles();
#elif defined(SHOCK_FINDER_BEFORE_OUTPUT)
      All.SnapshotFileCount = RestartSnapNum;
      char mybase[]="shocks";
      strncpy(All.SnapshotFileBase, mybase, MAXLEN_PARAM_VALUE);
      produce_dump();
#else
      terminate("Compile with flag SHOCK_FINDER_POST_PROCESSING or SHOCK_FINDER_BEFORE_OUTPUT!\n");
#endif
#ifdef DVR_RENDER
      dvr_render_loop("camerafile", 1024, 1024, 0);
#endif
      return 0;

    }
#ifdef VORONOI
  free_mesh();
#endif

  return -1;                    // return -1 means we ran to completion, i.e. not an endrun code
}


/*! \brief This routine computes the mass content of the box and compares it to the
 * specified value of Omega-matter.
 *
 * If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  double mass_b = 0, masstot_b, omega_b;
  int i, n_b = 0;

  for(i = 0; i < NumPart; i++)
    {
      mass += P[i].Mass;
      if(P[i].Type == 0)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#ifdef USE_SFR
      if(P[i].Type == 4)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#endif
#ifdef BLACK_HOLES
      if(P[i].Type == 5)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#endif
#ifdef DUST_LIVE
      if(P[i].Type == DUST_LIVE)
        {
          mass_b += P[i].Mass;
          n_b += 1;
        }
#endif
    }
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mass_b, &masstot_b, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega = masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));
  omega_b = masstot_b / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(n_b > 0)
    {
      if(fabs((omega - All.Omega0) / omega) > 1.0e-1 || fabs((omega_b - All.OmegaBaryon) / omega_b) > 1.0e-1)
        {
          mpi_printf("\n\nI've found something odd!\nThe mass content accounts for Omega=%g and OmegaBaryon=%g,\nbut you specified Omega=%g and OmegaBaryon=%g in the parameterfile.\n\nI better stop.\n", omega, omega_b, All.Omega0, All.OmegaBaryon);
#ifndef TWODIMS
          endrun();
#endif
        }

      if(fabs((omega - All.Omega0) / omega) > 1.0e-3 || fabs((omega_b - All.OmegaBaryon) / omega_b) > 1.0e-3)
	if(ThisTask == 0)
	  warn("I've found something odd! The mass content accounts for Omega=%g and OmegaBaryon=%g, but you specified Omega=%g and OmegaBaryon=%g in the parameterfile.", omega, omega_b, All.Omega0, All.OmegaBaryon);

    }
  else
    {
      if(All.OmegaBaryon != 0)
	if(ThisTask == 0)
	  warn("We are running with no baryons, even though you have specified OmegaBaryon=%g in the parameterfile. Please make sure you really want this.\n\n",All.OmegaBaryon);

      if(fabs((omega - All.Omega0) / omega) > 1.0e-1)
        {
          mpi_printf("\n\nI've found something odd!\nThe mass content accounts for Omega=%g and OmegaBaryon=%g,\nbut you specified Omega=%g and OmegaBaryon=%g in the parameterfile.\n\nI better stop.\n", omega, omega_b, All.Omega0, All.OmegaBaryon);
#ifndef TWODIMS
          endrun();
#endif
        }

      if(fabs((omega - All.Omega0) / omega) > 1.0e-3)
	if(ThisTask == 0)
	  warn("I've found something odd! The mass content accounts for Omega=%g and OmegaBaryon=%g, but you specified Omega=%g and OmegaBaryon=%g in the parameterfile.", omega, omega_b, All.Omega0, All.OmegaBaryon);
    }
}


/*! \brief This function is used to find an initial smoothing length for each SPH
 *  particle.
 *
 *  It guarantees that the number of neighbours will be between
 *  desired_ngb-MAXDEV and desired_ngb+MAXDEV. For simplicity, a first guess
 *  of the smoothing length is provided to the function density(), which will
 *  then iterate if needed to find the right smoothing length.
 */
void setup_smoothinglengths(void)
{
  int i, no, p;
  double *save_masses = mymalloc("save_masses", NumGas * sizeof(double));

  for(i = 0; i < NumGas; i++)
    {
#ifdef NO_GAS_SELFGRAVITY
/* This is needed otherwise the force tree will not be constructed for gas particles */
      P[i].Type = -1;
#endif
      save_masses[i] = P[i].Mass;
      P[i].Mass = 1.0;
    }

#ifdef HIERARCHICAL_GRAVITY
  TimeBinsGravity.NActiveParticles = 0;
  for(i = 0; i < NumGas; i++)
    {
      TimeBinsGravity.ActiveParticleList[TimeBinsGravity.NActiveParticles] = i;
      TimeBinsGravity.NActiveParticles++;
    }
#endif

  force_treeallocate(NumGas, All.MaxPart);
  force_treebuild(NumGas, 1, 0, 0);

  for(i = 0; i < NumGas; i++)
    {
      no = Father[i];

      if(no < 0)
        terminate("i=%d no=%d\n", i, no);

      while(10 * All.DesNumNgb * P[i].Mass > Nodes[no].u.d.mass)
        {
          p = Nodes[no].u.d.father;

          if(p < 0)
            break;

          no = p;
        }
#ifndef TWODIMS
      SphP[i].Hsml = pow(3.0 / (4 * M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 3) * Nodes[no].len;
#else
      SphP[i].Hsml = pow(1.0 / (M_PI) * All.DesNumNgb * P[i].Mass / Nodes[no].u.d.mass, 1.0 / 2) * Nodes[no].len;
#endif
#ifdef NO_GAS_SELFGRAVITY
/* Reset the original particle type */
      P[i].Type = 0;
#endif
    }

  myfree(Father);
  myfree(Nextnode);

#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
  myfree(Tree_Points);
  force_treefree();

  density();

  for(i = 0; i < NumGas; i++)
    P[i].Mass = save_masses[i];

  myfree(save_masses);

  for(i = 0; i < NumGas; i++)
    SphP[i].MaxDelaunayRadius = SphP[i].Hsml;


#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  ngb_treefree();
  domain_free();
  domain_Decomposition();
  ngb_treeallocate();
  ngb_treebuild(NumGas);
#endif
}

/*! \brief This function checks for unique particle  IDs
 *
 *  The particle IDs are copied to an array and then sorted among all tasks.
 *  This array is then checked for duplicates. In that case the code terminates.
 */
void test_id_uniqueness(void)
{
  int i;
  double t0, t1;
  MyIDType *ids, *ids_first;

  mpi_printf("INIT: Testing ID uniqueness...\n");

  if(NumPart == 0)
    terminate("need at least one particle per cpu\n");

  t0 = second();

#ifndef SPH_BND_PARTICLES
  ids = (MyIDType *) mymalloc("ids", NumPart * sizeof(MyIDType));
  ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

  for(i = 0; i < NumPart; i++)
    ids[i] = P[i].ID;

  parallel_sort(ids, NumPart, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < NumPart; i++)
    {
#ifdef TGSET
      if(ids[i] == 0)
        continue;
#endif
#ifdef SPECIAL_BOUNDARY
      if(ids[i] == -2 || ids[i] == -1)
        continue;
#endif
      if(ids[i] == ids[i - 1])
        terminate("non-unique ID=%lld found on task=%d (i=%d NumPart=%d)\n", (long long) ids[i], ThisTask, i, NumPart);
    }
  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
#ifdef SPECIAL_BOUNDARY
    if(ids_first[ThisTask + 1] != -2 && ids_first[ThisTask + 1] != -1)
#endif
      {
        if(ids[NumPart - 1] == ids_first[ThisTask + 1])
          terminate("non-unique ID=%lld found on task=%d\n", (long long) ids[NumPart - 1], ThisTask);
      }
  myfree(ids_first);
  myfree(ids);
#endif

  t1 = second();

  mpi_printf("INIT: success.  took=%g sec\n", timediff(t0, t1));
}


#ifdef TRACER_MC
/*! \brief This function checks for unique monte carlo tracer IDs
 *
 *  The tracer IDs are copied to an array and then sorted among all tasks.
 *  This array is then checked for duplicates. In that case the code terminates.
 */
void test_tracer_mc_id_uniqueness(void)
{
  int i, j;
  double t0, t1;
  MyIDType *ids, *ids_first;
  int new_n_tracer;  

  mpi_printf("TRACERS: Testing tracer ID uniqueness...\n");

  t0 = second();

  if(N_tracer == 0)
    {
      ids = (MyIDType *) mymalloc("ids", sizeof(MyIDType));
      ids[0] = -1;
      new_n_tracer = 1;
    }
  else
    {  

      ids = (MyIDType *) mymalloc("ids", N_tracer * sizeof(MyIDType));

      for(i = 0, j = 0; i < NumPart; i++)
        {
          int next = P[i].TracerHead;
          while(next >= 0)
            {
              ids[j++] = TracerLinkedList[next].ID;
              next = TracerLinkedList[next].Next;
            }
        }
      new_n_tracer = N_tracer;
    }

  ids_first = (MyIDType *) mymalloc("ids_first", NTask * sizeof(MyIDType));

  parallel_sort(ids, new_n_tracer, sizeof(MyIDType), compare_IDs);

  for(i = 1; i < new_n_tracer; i++)
    {
      if(ids[i] != -1 && ids[i] == ids[i - 1])
        terminate("non-unique ID=%lld found on task=%d (i=%d NumTracer=%d)\n", (long long) (ids[i]), ThisTask, i, N_tracer);
    }
  MPI_Allgather(&ids[0], sizeof(MyIDType), MPI_BYTE, ids_first, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  if(ThisTask < NTask - 1)
    {
      if(ids_first[ThisTask + 1] != -1 && ids[new_n_tracer - 1] == ids_first[ThisTask + 1])
        terminate("non-unique ID=%lld found on task=%d\n", (long long) ids[new_n_tracer - 1], ThisTask);
    }

  myfree(ids_first);
  myfree(ids);

  t1 = second();

  mpi_printf("TRACERS: success.  took=%g sec\n", timediff(t0, t1));
}
#endif

#ifdef INJECT_TRACER_INTO_SN
void calculate_max_tracer_mc_id(void)
{
  int i,j;
  MyIDType local_max = 0;

  if(N_tracer != 0)
    {  
      for(i = 0, j = 0; i < NumPart; i++)
      {
        int next = P[i].TracerHead;
        while(next >= 0)
          {
            if(TracerLinkedList[next].ID > local_max)
              local_max = TracerLinkedList[next].ID;
            next = TracerLinkedList[next].Next;
          }
      }
    }

  MPI_Allreduce(&local_max, &MaxTracerID, 1, MPI_MYIDTYPE, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("MaxTracerID = %lld \n",(long long) MaxTracerID);

}
#endif

void calculate_maxid(void)
{
  /* determine maximum ID */
  MyIDType maxid, *tmp;
  int i;

  for(i = 0, maxid = 0; i < NumPart; i++)
    if(P[i].ID > maxid)
      {
#ifdef SPECIAL_BOUNDARY
        if(P[i].ID == -1 || P[i].ID == -2)
          continue;
#endif
#ifdef BOUNDARY_INFLOWOUTFLOW_MINID
        if(P[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID)
          continue;
#endif
#ifdef BOUNDARY_REFL_FLUIDSIDE_MINID
        if(P[i].ID >= BOUNDARY_REFL_FLUIDSIDE_MINID)
          continue;
#endif
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
        if(P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID)
          continue;
#endif
#ifdef BOUNDARY_STICKY_MINID
        if(P[i].ID >= BOUNDARY_STICKY_MINID)
          continue;
#endif
        maxid = P[i].ID;
      }

  tmp = mymalloc("tmp", NTask * sizeof(MyIDType));

  MPI_Allgather(&maxid, sizeof(MyIDType), MPI_BYTE, tmp, sizeof(MyIDType), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    if(tmp[i] > maxid)
      maxid = tmp[i];

#if defined(REFINEMENT_SPLIT_CELLS) || defined(USE_SFR)
  All.MaxID = maxid;
#endif

  myfree(tmp);
}


#if defined(BLACK_HOLES) && defined(FOF)
void set_DMmass(void)
{
  int i;
  double MaxmassDMpart = 0.0;
  double massDMpart = 0.0;

  if(All.MassTable[1] > 0)
    {
      massDMpart = All.MassTable[1];
    }
  else
    {
      for(i = 0; i < NumPart; i++)
        if(P[i].Type == 1)
          {
            if(massDMpart < P[i].Mass)
              massDMpart = P[i].Mass;
          }
    }

  MPI_Allreduce(&massDMpart, &MaxmassDMpart, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  All.massDMpart = MaxmassDMpart;
}
#endif

/*! \brief The sorting kernel used during the ID uniqueness check
 *
 */
int compare_IDs(const void *a, const void *b)
{
  if(*((MyIDType *) a) < *((MyIDType *) b))
    return -1;

  if(*((MyIDType *) a) > *((MyIDType *) b))
    return +1;

  return 0;
}
