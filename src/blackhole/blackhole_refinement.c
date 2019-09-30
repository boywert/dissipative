/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/blackhole/blackhole_refinement.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *      01/12/2015 Mike Curtis - updated with new changes to black hole refinement scheme.
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

#if defined(BLACK_HOLES) && defined(REFINEMENT_AROUND_BH)


static int blackhole_mark_cells_for_refinement_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyFloat Rbondi;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->BH_Hsml = BPP(i).BH_Hsml;

  double soundspeed = sqrt(GAMMA * GAMMA_MINUS1 * BPP(i).BH_U);
  in->Rbondi = 50.0 * (PARSEC / All.UnitLength_in_cm) * ((BPP(i).BH_Mass * All.UnitMass_in_g / All.HubbleParam) / (1.e7 * SOLAR_MASS)) *
               pow(((soundspeed * All.UnitVelocity_in_cm_per_s) / (30.0 * 1.e5)), -2);

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
} data_out;

static data_out *DataResult, *DataOut;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{

}


#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int idx;
  #pragma omp parallel private(i, idx)
  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        #pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int i = TimeBinsBHAccretion.ActiveParticleList[idx];
        if(i < 0)
          continue;

        if(BPP(i).SwallowID == 0)
          blackhole_mark_cells_for_refinement_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, count = 0;
  #pragma omp parallel private(i)
  {
    int threadid = get_thread_num();

    while(1)
      {
        #pragma omp atomic capture
        i = count++;

        if(i >= Nimport)
          break;

        blackhole_mark_cells_for_refinement_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void blackhole_mark_cells_for_refinement(void)
{

  int idx, n;
  mpi_printf("REFINEMENT_AROUND_BH: Begin marking cells for refinement.\n");

  /* reset cell refinement flag */
  for(int i = 0; i < NumGas; i++)
    {
      SphP[i].RefBHFlag = 0;
      SphP[i].RefBHMaxRad = MAX_REAL_NUMBER;
    }

  generic_set_MaxNexport();
  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);

  int refcount_can, refcount;
  long long totrefcount, totrefcount_can;

  for(int i = 0, refcount = refcount_can = 0; i < NumGas; i++)
    if(SphP[i].RefBHFlag)
      {
        refcount++;
#if (REFINEMENT_AROUND_BH==0)
        if(can_this_cell_be_split(i))
          refcount_can++;
#endif
#if (REFINEMENT_AROUND_BH==1)
        refcount_can++;
#endif
      }

  sumup_large_ints(1, &refcount, &totrefcount);
  sumup_large_ints(1, &refcount_can, &totrefcount_can);
  mpi_printf("REFINEMENT_AROUND_BH: all cells -> %lld/%lld want/can refine near BHs\n", totrefcount, totrefcount_can);

  refcount = refcount_can = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].RefBHFlag)
        {
          refcount++;
#if (REFINEMENT_AROUND_BH==0)
          if(can_this_cell_be_split(i))
            refcount_can++;
#endif
#if (REFINEMENT_AROUND_BH==1)
          refcount_can++;
#endif
        }
    }

  sumup_large_ints(1, &refcount, &totrefcount);
  sumup_large_ints(1, &refcount_can, &totrefcount_can);
  mpi_printf("REFINEMENT_AROUND_BH: active cells -> %lld/%lld want/can refine near BHs\n", totrefcount, totrefcount_can);
  mpi_printf("REFINEMENT_AROUND_BH: Cells now marked for refinement.");
}


int blackhole_mark_cells_for_refinement_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode, j, k, n;

  MyDouble dx, dy, dz, r2;

  double xtmp, ytmp, ztmp;   // required for NGB_PERIODIC_LONG_X

  // for data in from BH
  MyDouble *pos;
  MyFloat hsml, rbondi;

  data_in local, *in;
  data_out out;

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
  hsml = in->BH_Hsml;
  rbondi = in->Rbondi;

#ifdef REFINEMENT_AROUND_BH_FIXED
  double refinement_radius =  All.RefBHRadius;
#else
  double refinement_radius =  All.RefBHRadiusHSML * hsml;
#endif
  int nfound = ngb_treefind_variable_threads(pos, refinement_radius, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < refinement_radius * refinement_radius && P[j].Mass != 0 && P[j].ID != 0)
        {
#ifdef REFINEMENT_AROUND_BH_FIXED
          MyFloat max_rad = All.RefBHMaxCellRadius;
          MyFloat min_rad = All.RefBHMinCellRadius;
#else
          MyFloat max_rad = All.RefBHMaxCellRadiusHSML * hsml;
          MyFloat min_rad = All.RefBHMinCellRadiusRBondi * rbondi;
#endif
          if(min_rad > max_rad || !isfinite(rbondi))
            {
              warn("The bondi radius for a black hole is non finite. Setting min_rad = max_rad/10000.");
              min_rad = max_rad/10000;
            }

          double r_norm = sqrt(r2) / refinement_radius;
          MyFloat cell_max_rad = min_rad + r_norm * (max_rad-min_rad);

          if(SphP[j].RefBHFlag)
            {
              if(cell_max_rad < SphP[j].RefBHMaxRad)
                {
                  SphP[j].RefBHMaxRad = cell_max_rad;
                }
            }
          else
            {
              SphP[j].RefBHFlag = 1;
              SphP[j].RefBHMaxRad = cell_max_rad;
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
  return 0;
}

#endif
