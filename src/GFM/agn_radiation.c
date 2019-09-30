/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/agn_radiation.c
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
#include <gsl/gsl_math.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(BLACK_HOLES) && defined(GFM_AGN_RADIATION)

#define GFM_AGN_RAD_RVIR_FAC   3.0



static MyDouble ADAFThresholdMdotInEddington, ADAFparam1, ADAFparam2;
static MyDouble *Radius;

static int LocalRadiativeInefficientCount, LocalRadiativeEfficientCount;
static double LocalMinMdotInEddington, LocalMaxMdotInEddington;
static int LocalFlagActive;

static int assign_agn_radiation_evaluate(int target, int mode, int threadid);
static double get_agn_luminosity(MyFloat bh_mdot, MyFloat bh_mass);


typedef struct
{
  MyDouble Pos[3];
  MyDouble Radius;
  double AGNLuminosity;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;



/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];

  in->AGNLuminosity = get_agn_luminosity(BPP(i).BH_Mdot, BPP(i).BH_Mass);
  in->Radius = Radius[i];

  in->Firstnode = firstnode;
}


/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult, *DataOut;



 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
    }
  else                          /* merge */
    {
    }
}



#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, idx)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int i = TimeBinsBHAccretion.ActiveParticleList[idx];
        if(i < 0)
          continue;

        assign_agn_radiation_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }


}

static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        assign_agn_radiation_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}


static double get_agn_luminosity(MyFloat bh_mdot, MyFloat bh_mass)
{
  double Lbol = 0;

  if((bh_mass == 0) || (bh_mdot == 0))
    return Lbol;

  double meddington = (4 * M_PI * GRAVITY * CLIGHT * PROTONMASS / (All.BlackHoleRadiativeEfficiency * CLIGHT * CLIGHT * THOMPSON)) * bh_mass * All.UnitTime_in_s / All.HubbleParam;

  /* accretion rate in units of Eddingtion */
  double mdot_in_edd = bh_mdot / meddington;

  /* accretion rate in cgs units */
  double mdot_in_cgs = bh_mdot * All.UnitMass_in_g / All.UnitTime_in_s;

  /* statistics */
  if(mdot_in_edd > 0)
    {
      LocalFlagActive = 1;
      if(mdot_in_edd < LocalMinMdotInEddington)
        LocalMinMdotInEddington = mdot_in_edd;
      if(mdot_in_edd > LocalMaxMdotInEddington)
        LocalMaxMdotInEddington = mdot_in_edd;
    }

  /* radiatively efficient */
  if(mdot_in_edd > ADAFThresholdMdotInEddington)
    {
      Lbol = (1. - All.BlackHoleFeedbackFactor) * All.BlackHoleRadiativeEfficiency * mdot_in_cgs * CLIGHT * CLIGHT;
      LocalRadiativeEfficientCount++;
    }
  /* ADAF efficiency (Narayan & Yi 1994) with continuous transition to All.BlackHoleRadiativeEfficiency */
  else
    {
      Lbol = ADAFparam1 * ADAFparam2 * mdot_in_edd / (1.0 + ADAFparam2 * mdot_in_edd) * All.BlackHoleRadiativeEfficiency * mdot_in_cgs * CLIGHT * CLIGHT;
      LocalRadiativeInefficientCount++;
    }

  AGNEnergyEM_Is += Lbol;

  /* Hopkins et al. 2007 obscuration scaling */
  Lbol *= All.ObscurationFactor * pow((Lbol / 1.0e46), All.ObscurationSlope);

  AGNEnergyEMobs_Is += Lbol;

  if(!gsl_finite(Lbol))
    terminate("GFM_AGN_RADIATION: Lbol=%g mdot_in_edd=%g mdot_in_cgs=%g bh_mdot=%g meddington=%g bh_mass=%g", Lbol, mdot_in_edd, mdot_in_cgs, bh_mdot, meddington, bh_mass);

  return Lbol;
}


void assign_agn_radiation(void)
{
  int i, idx;
  double nH_mean, agnlum, radius_by_ionization_param;
  double ionization_param_min = GFM_MIN_IONIZATION_PARAMETER;
  double r200c, rhocrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  /* physical critical density at current redshift */
  rhocrit *= (All.Omega0 * pow(1 + All.cf_redshift, 3) + (1 - All.Omega0 - All.OmegaLambda) * pow(1 + All.cf_redshift, 2) + All.OmegaLambda);

  for(i = 0; i < NumGas; i++)
    SphP[i].AGNBolIntensity = 0;

  /* mean hydrogen number density */
  nH_mean = HYDROGEN_MASSFRAC * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G) / PROTONMASS * All.UnitDensity_in_cgs;
  nH_mean /= All.cf_atime * All.cf_atime * All.cf_atime;

  /* set up sphere of influence for each BH based on ionization parameter */

  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      i = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(i < 0)
        continue;

      agnlum = get_agn_luminosity(BPP(i).BH_Mdot, BPP(i).BH_Mass);

      r200c = pow(BPP(i).HostHaloMass / (4 * M_PI / 3.0 * 200 * rhocrit), 1.0 / 3.0);   /* physical r_200,crit value, assuming FoF mass = M_200,crit */
      radius_by_ionization_param = sqrt(agnlum / (ionization_param_min * nH_mean)) / All.UnitLength_in_cm;
      Radius[i] = dmin(radius_by_ionization_param, BPP(i).HostHaloMass == 0 ? radius_by_ionization_param : GFM_AGN_RAD_RVIR_FAC * r200c);
    }

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);

}


int assign_agn_radiation_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  double dx, dy, dz, r2, h, h2;
  MyDouble *pos;
  double agnlum, local_intensity = 0;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  h = target_data->Radius;
  agnlum = target_data->AGNLuminosity;

  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2 && r2 > 0)
        {
          local_intensity = agnlum / r2;

          /* assign bolometric intensity erg/s/cm^2 in physical units */
          local_intensity /= pow(All.cf_atime * All.UnitLength_in_cm / All.HubbleParam, 2);

          SphP[j].AGNBolIntensity += local_intensity;
        }
    }

  return 0;
}

/*! \brief setup ADAF regime
 */
void init_agn_radiation(void)
{
#if defined(UNIFIED_FEEDBACK) && defined(BH_BUBBLES)
  ADAFThresholdMdotInEddington = All.RadioThreshold;
  if(All.RadioFeedbackFactor < 1.0)
    ADAFparam1 = 2.0 * (1.0 - All.RadioFeedbackFactor);
  else
    {
      ADAFparam1 = 0.0;
    }
  ADAFparam2 = 1.0 / All.RadioThreshold;
#else

#if defined(BH_ADIOS_WIND)

#if defined(BH_ADIOS_WIND_WITH_QUASARTHRESHOLD)
  ADAFThresholdMdotInEddington = All.QuasarThreshold;
  ADAFparam2 = 1.0 / All.QuasarThreshold;
#else
  ADAFThresholdMdotInEddington = 0.01;
  ADAFparam2 = 100.0;
#endif
  if(All.RadioFeedbackFactor < 1.0)
    ADAFparam1 = 2.0 * (1.0 - All.RadioFeedbackFactor);
  else
    {
      ADAFparam1 = 0.0;
    }

#else

  ADAFThresholdMdotInEddington = 0.01;
  ADAFparam1 = 2.0;
  ADAFparam2 = 100.0;

#endif


#endif
}

void agn_radiation_info(void)
{
  TIMER_START(CPU_GFM_AGNRAD);

  long long totCellsWithAGNRadiation;
  sumup_large_ints(1, &CellsWithAGNRadiation, &totCellsWithAGNRadiation);
  mpi_printf("GFM_AGN_RADIATION: total number of cells with AGN background = %lld\n", totCellsWithAGNRadiation);

  TIMER_STOP(CPU_GFM_AGNRAD);
}

void do_agn_radiation(void)
{
  TIMER_START(CPU_GFM_AGNRAD);

  int totRadiativeEfficientCount, totRadiativeInefficientCount;
  int totFlagActive;
  double totMinMdotInEddington, totMaxMdotInEddington;
  int total_count;

  if(All.HighestSynchronizedTimeBin >= All.HighestOccupiedTimeBin)
    {
      /* we only do this on steps where all BHs are synchronized, and to avoid a global operation for the whole system at all times */

      /* for feedback statistics */
      AGNEnergyEM_Is = 0.0;
      AGNEnergyEMobs_Is = 0.0;

      /* allocate temp arrays */
      Radius = (MyDouble *) mymalloc("Radius", NumPart * sizeof(MyDouble));

      /* reset */
      LocalRadiativeEfficientCount = 0;
      LocalRadiativeInefficientCount = 0;
      LocalMinMdotInEddington = +MAX_REAL_NUMBER;
      LocalMaxMdotInEddington = -MAX_REAL_NUMBER;
      LocalFlagActive = 0;

      /* assign AGN luminosty to cells */
      assign_agn_radiation();

      /* collect accretion statistics */
      MPI_Reduce(&LocalRadiativeEfficientCount, &totRadiativeEfficientCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocalRadiativeInefficientCount, &totRadiativeInefficientCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocalFlagActive, &totFlagActive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocalMinMdotInEddington, &totMinMdotInEddington, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&LocalMaxMdotInEddington, &totMaxMdotInEddington, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      total_count = totRadiativeEfficientCount + totRadiativeInefficientCount;

      if(totFlagActive > 0)
        mpi_printf("GFM_AGN_RADIATION: FlagActive=%d  RadiativeEfficient=%d (%g%)  RadiativeInefficient=%d (%g%)  MinMdotInEddington=%g  MaxMdotInEddington=%g\n",
                   totFlagActive, totRadiativeEfficientCount, 100.0 * totRadiativeEfficientCount / total_count, totRadiativeInefficientCount,
                   100.0 * totRadiativeInefficientCount / total_count, totMinMdotInEddington, totMaxMdotInEddington);
      else
        mpi_printf("GFM_AGN_RADIATION: non active\n");

      /* free temp arrays */
      myfree(Radius);
    }

  TIMER_STOP(CPU_GFM_AGNRAD);
}

#endif
