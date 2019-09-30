/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/winds.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)

static int find_wind_cells_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];              /* wind particle position */
  MyFloat hsml;                 /* current search radius */

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[WindParticle[i].index].Pos[0];
  in->Pos[1] = P[WindParticle[i].index].Pos[1];
  in->Pos[2] = P[WindParticle[i].index].Pos[2];
  in->hsml = STP(WindParticle[i].index).Hsml;

  in->Firstnode = firstnode;
}


typedef struct
{
  MyFloat NearestDist;          /* distance to closest cell on task */
  MyFloat Density;              /* density of nearest cell */
  MyIDType CellID;              /* ID of cell */
} data_out;

static data_out *DataResult, *DataOut;



/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      WindParticle[i].NearestDist = out->NearestDist;
      WindParticle[i].CellID = out->CellID;
      WindParticle[i].Density = out->Density;
    }
  else                          /* merge */
    {
      if(out->NearestDist < WindParticle[i].NearestDist)
        {
          WindParticle[i].NearestDist = out->NearestDist;
          WindParticle[i].CellID = out->CellID;
          WindParticle[i].Density = out->Density;
        }
    }
}



#include "../generic_comm_helpers2.h"

static int Nwind;

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i)
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
              if(generic_polling_primary(count, Nwind))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Nwind)
          break;

        if(WindParticle[i].CellID == 0) /* do we have the cell already for that wind particle? */
          find_wind_cells_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        find_wind_cells_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}



void find_wind_cells(int nwind)
{
  Nwind = nwind;

  long long ntotWind;
  sumup_large_ints(1, &Nwind, &ntotWind);
  if(ntotWind == 0)
    return;

  int i;
  for(i = 0; i < Nwind; i++)
    {
      WindParticle[i].NearestDist = MAX_REAL_NUMBER;
      WindParticle[i].CellID = 0;
    }

  generic_set_MaxNexport();

  int ntot, iter = 0;

  /* we will repeat the whole thing for those points where we did not find a nearest neighbor */
  do
    {
      generic_comm_pattern(Nwind, kernel_local, kernel_imported);


      int npleft;

      /* do final operations on results */
      for(i = 0, npleft = 0; i < Nwind; i++)
        {
          if(WindParticle[i].CellID == 0)
            {
              npleft++;
              STP(WindParticle[i].index).Hsml *= 2.0;
              if(iter >= MAXITER - 10)
                {
                  printf("i=%d task=%d hsml=%g nearest dist=%g pos=(%g|%g|%g)\n", i, ThisTask,
                         STP(WindParticle[i].index).Hsml, WindParticle[i].NearestDist, P[WindParticle[i].index].Pos[0], P[WindParticle[i].index].Pos[1], P[WindParticle[i].index].Pos[2]);
                  myflush(stdout);
                }
              if(iter > MAXITER)
                terminate("iter > MAXITER");
            }
          else
            /* note: NearestDist can be zero if we start from a snapshot, since a snapshot can contain wind and SPH particles at identical coordinates, since
               after wind particle creation only a gravity kick is applied, but no drift before the snapshot is written
             */
            STP(WindParticle[i].index).Hsml = (WindParticle[i].NearestDist > 0) ? (1.5 * WindParticle[i].NearestDist) : get_default_softening_of_particletype(4);
        }

      /* sum up the left overs */
      MPI_Allreduce(&npleft, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(ntot > 0)              /* ok, we need to repeat for a few wind particles */
        {
          iter++;
          if(iter > 0)
            mpi_printf("GFM_WINDS: wind-nearest iteration %3d: need to repeat for %010d particles.\n", iter, ntot);

          if(iter > MAXITER)
            terminate("GFM_WINDS: failed to converge in wind-nearest\n");
        }
    }
  while(ntot > 0);
}


static int find_wind_cells_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode;
  double h;
  MyDouble dx, dy, dz, r;
  MyDouble *pos;
  MyIDType cellID = 0;
  MyFloat density = MAX_FLOAT_NUMBER;
  double nearestDist = MAX_REAL_NUMBER;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;
  data_out out;

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
  h = target_data->hsml;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
      dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

      r = sqrt(dx * dx + dy * dy + dz * dz);

      if(r < nearestDist && r < h && P[j].ID != 0 && P[j].Mass > 0)     /* check for distance, and ignore merged cells */
        {
          nearestDist = r;
          cellID = P[j].ID;
          density = SphP[j].Density;
        }
    }

  out.NearestDist = nearestDist;
  out.CellID = cellID;
  out.Density = density;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#endif
