/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/radiation_stellar_feedback_find_cells.c
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


#if defined(FM_RADIATION_FEEDBACK)

/* communication structures */
typedef struct
{
  MyDouble Pos[3];
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  MyFloat Hsml;
#else
  MyFloat StromgrenRadius;
#endif
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  MyDouble NormSphRadFeedback;
  int RadFeed_NumNgb;
  MyFloat RadFeed_MinGasDist;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  MyDouble NormSphRadFeedback_cold;
  MyFloat LowestDensity;
  MyFloat LowestDensityDirection_x;
  MyFloat LowestDensityDirection_y;
  MyFloat LowestDensityDirection_z;
#endif
} data_out;

static data_out *DataResult, *DataOut;

static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[StarParticle[i].index].Pos[0];
  in->Pos[1] = P[StarParticle[i].index].Pos[1];
  in->Pos[2] = P[StarParticle[i].index].Pos[2];
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  if( STP(StarParticle[i].index).Hsml > STP(StarParticle[i].index).StromgrenRadius)
      in->Hsml = STP(StarParticle[i].index).Hsml;
  else
      in->Hsml = STP(StarParticle[i].index).StromgrenRadius;
#else
  in->StromgrenRadius = STP(StarParticle[i].index).StromgrenRadius;
#endif
  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)
    {
      StarParticle[i].NormSphRadFeedback = out->NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      StarParticle[i].NormSphRadFeedback_cold = out->NormSphRadFeedback_cold;
#endif
      StarParticle[i].RadFeed_NumNgb = out->RadFeed_NumNgb;

      STP(StarParticle[i].index).RadFeed_NumNgb = StarParticle[i].RadFeed_NumNgb;
      StarParticle[i].RadFeed_MinGasDist = dmin(out->RadFeed_MinGasDist, StarParticle[i].RadFeed_MinGasDist);
#ifdef FM_RADIATION_FEEDBACK_DEBUG
      STP(StarParticle[i].index).NormSphRadFeedback = StarParticle[i].NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      STP(StarParticle[i].index).NormSphRadFeedback_cold = StarParticle[i].NormSphRadFeedback_cold;
#endif
#endif
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      StarParticle[i].LowestDensity = out->LowestDensity;
      StarParticle[i].LowestDensityDirection_x = out->LowestDensityDirection_x;
      StarParticle[i].LowestDensityDirection_y = out->LowestDensityDirection_y;
      StarParticle[i].LowestDensityDirection_z = out->LowestDensityDirection_z;
#endif
    }
  else
    {
      StarParticle[i].NormSphRadFeedback += out->NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      StarParticle[i].NormSphRadFeedback_cold += out->NormSphRadFeedback_cold;
#endif
      StarParticle[i].RadFeed_NumNgb += out->RadFeed_NumNgb;

      STP(StarParticle[i].index).RadFeed_NumNgb = StarParticle[i].RadFeed_NumNgb;
      StarParticle[i].RadFeed_MinGasDist = dmin(out->RadFeed_MinGasDist, StarParticle[i].RadFeed_MinGasDist);
#ifdef FM_RADIATION_FEEDBACK_DEBUG
      STP(StarParticle[i].index).NormSphRadFeedback += StarParticle[i].NormSphRadFeedback;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      STP(StarParticle[i].index).NormSphRadFeedback_cold += StarParticle[i].NormSphRadFeedback_cold;
#endif
#endif
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
      if(out->LowestDensity < StarParticle[i].LowestDensity)
      {
          StarParticle[i].LowestDensity = out->LowestDensity;
          StarParticle[i].LowestDensityDirection_x = out->LowestDensityDirection_x;
          StarParticle[i].LowestDensityDirection_y = out->LowestDensityDirection_y;
          StarParticle[i].LowestDensityDirection_z = out->LowestDensityDirection_z;
      }
#endif
    }
}


#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, Nstar))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nstar)
          break;

#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
        if(StarParticle[i].StromgrenMass >= 0.0)
          find_radiation_feedback_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
#else
        if(StarParticle[i].StromgrenRadius >= 0.0)
          find_radiation_feedback_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
#endif
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

        find_radiation_feedback_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void find_radiation_feedback_cells(void)
{
  long long ntot;

  sumup_large_ints(1, &Nstar, &ntot);
  if(ntot == 0)
    return;


  mpi_printf("LVS: start find_radiation_feedback_cells\n");

  generic_set_MaxNexport();

  double t0 = second();

  generic_comm_pattern(Nstar, kernel_local, kernel_imported);

#if !defined(FM_STOCHASTIC_HII_PHOTOIONIZATION)
  for(int j = 0; j < Nstar; j++)
    {
      if(StarParticle[j].StromgrenRadius > 0)
        {
          if(StarParticle[j].NormSphRadFeedback == 0)
            {
              //printf("Problem LVS. Task=%d ID=%d h=%g NormSph=%g NumNgb=%d \n",ThisTask, j, StarParticle[j].StromgrenRadius, StarParticle[j].NormSphRadFeedback, StarParticle[j].RadFeed_NumNgb);    
              //terminate("sucks");
              printf("Problem LVS0 Stromgren rad. too small. rs=%g, min gas dist=%g \n", StarParticle[j].StromgrenRadius, StarParticle[j].RadFeed_MinGasDist);
              STP(StarParticle[j].index).RadFeed_NumNgb = 1;
              StarParticle[j].NormSphRadFeedback = All.TargetGasMass;
              StarParticle[j].StromgrenRadius = 1.1 * StarParticle[j].RadFeed_MinGasDist;
            }
        }
    }
#endif

  double t1 = second();

  mpi_printf("LVS: Done! Calculation took %g sec\n", timediff(t0, t1));
}


/* Modified to compute the total mass of cells within the stromgren radius  
 *
 */
int find_radiation_feedback_cells_evaluate(int target, int mode, int thread_id)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;

  double h, h2, hinv, hinv3;
  double dx, dy, dz, r2;
  double minDistGas = MAX_REAL_NUMBER;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  MyFloat normsph=0;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  double lowest_density = MAX_REAL_NUMBER;
  double lowest_density_direction_x,  lowest_density_direction_y,  lowest_density_direction_z;
  MyFloat normsph_cold=0;
#endif
  MyDouble *pos;

  int countcells = 0;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;
      numnodes = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }
  pos = in->Pos;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  h = in->Hsml;
  if(h > 0.1) h = 0.1; // ToDo:  At large distances the ionizing photon input rate needs to account for the geometric incidence of radiation
#else
  h = in->StromgrenRadius;
#endif

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, thread_id, numnodes, firstnode);

  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[thread_id].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0) /* skip cells that have been swallowed or dissolved */
        {
          dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
          dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
          dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

          r2 = dx * dx + dy * dy + dz * dz;

          if(r2 < minDistGas)
            minDistGas = r2;    /*is the square of r */

          if(r2 < h2)
            {
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
              if(SphP[j].Density < lowest_density)
              {
                lowest_density = SphP[j].Density;
                lowest_density_direction_x = dx;
                lowest_density_direction_y = dy;
                lowest_density_direction_z = dz;
              }
              /* this is a more liberal cut to identify gas that could be photoionized to avoid.  Avoids photoionization flickering. */
              if(SphP[j].Utherm < 2.0*All.PhotoionizationEgySpec && SphP[j].GasRadCoolShutoffTime==0.0)  normsph_cold += P[j].Mass;
#endif
              normsph += P[j].Mass;     /*notice that we store mass and not volume */
              countcells++;
            }

        }
    }

  out.NormSphRadFeedback = normsph;
  out.RadFeed_NumNgb = countcells;
  out.RadFeed_MinGasDist = sqrt(minDistGas);

#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  out.NormSphRadFeedback_cold = normsph_cold;
  out.LowestDensity = lowest_density;
  out.LowestDensityDirection_x = lowest_density_direction_x;
  out.LowestDensityDirection_y = lowest_density_direction_y;
  out.LowestDensityDirection_z = lowest_density_direction_z;
#endif
    
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
