/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_3d.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gmp.h>

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#if defined (VORONOI) && !defined(TWODIMS) && !defined(ONEDIMS)
#ifdef VORONOI_TEST

void check_circumcircle(tessellation * T, int tt);
void re_compute_circumcircles(tessellation * T);
void get_circumcircle_exact(tessellation * T, int tt, double *x, double *y, double *z);
int voronoi_search_check(tessellation * TT);
static int voronoi_test_evaluate(int target, int mode, int threadid);

static int *NgbList;

void voronoi_test(void)
{
  mpi_printf("VORONOI_TEST:\n");

  re_compute_circumcircles(&Mesh);

  voronoi_search_check(&Mesh);

  mpi_printf("VORONOI_TEST completed.\n");
}


void re_compute_circumcircles(tessellation * T)
{
  tetra *DT = T->DT;
  int i;

  for(i = 0; i < T->Ndt; i++)
    {
      if(DT[i].t[0] < 0)        /* deleted ? */
        continue;

      if(DT[i].p[0] == DPinfinity)
        continue;
      if(DT[i].p[1] == DPinfinity)
        continue;
      if(DT[i].p[2] == DPinfinity)
        continue;
      if(DT[i].p[3] == DPinfinity)
        continue;

      check_circumcircle(T, i);
    }
}

void check_circumcircle(tessellation * T, int tt)
{
  tetra *DT = T->DT;
  tetra_center *DTC = T->DTC;
  point *DP = T->DP;
  tetra *t = &DT[tt];
  tetra_center *tc = &DTC[tt];

  if(t->t[0] < 0)               /* deleted ? */
    return;

  point *p0 = &DP[t->p[0]];
  point *p1 = &DP[t->p[1]];
  point *p2 = &DP[t->p[2]];
  point *p3 = &DP[t->p[3]];

  if(isInfinity(p0) || isInfinity(p1) || isInfinity(p2) || isInfinity(p3))
    return;

  int orient = Orient3d_Exact(p0, p1, p2, p3);
  double vol = calculate_tetra_volume(p0, p1, p2, p3);

  if(orient != 1)
    terminate("orientation is wrong\n");


  double xc, yc, zc;
  get_circumcircle_exact(T, tt, &xc, &yc, &zc);

  double dx = tc->cx - xc;
  double dy = tc->cy - yc;
  double dz = tc->cz - zc;

  double r2 = dx * dx + dy * dy + dz * dz;

  double r = sqrt(r2);
  double len = pow(6 * vol, 1.0 / 3);

  if(r > 1.0e-6 * len)
    printf("ThisTask=%d: tt=%d  r/len=%g\n", ThisTask, tt, r / len);
}








static tessellation *T;

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Rad;

  int Firstnode;

} data_in;

static data_in *DataGet, *DataIn;



/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  point *DP = T->DP;
  tetra *DT = T->DT;
  tetra_center *DTC = T->DTC;

  in->Pos[0] = DTC[i].cx;
  in->Pos[1] = DTC[i].cy;
  in->Pos[2] = DTC[i].cz;

  double dx = DTC[i].cx - DP[DT[i].p[0]].x;
  double dy = DTC[i].cy - DP[DT[i].p[0]].y;
  double dz = DTC[i].cz - DP[DT[i].p[0]].z;

  in->Rad = sqrt(dx * dx + dy * dy + dz * dz);


  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  int NumNgb;
  /* counts how many have been found */
} data_out;

static data_out *DataResult, *DataOut;


/* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      NgbList[i] = out->NumNgb;
    }
  else                          /* merge */
    {
      NgbList[i] += out->NumNgb;
    }
}

#include "generic_comm_helpers2.h"


static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
  /* do local particles and prepare export list */
#pragma omp parallel private(i)
  {
    int thread_id = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    for(int j = 0; j < NTask; j++)
      Thread[thread_id].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[thread_id].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, T->Ndt))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= T->Ndt)
          break;

        tetra *DT = T->DT;
        point *DP = T->DP;

        if(DT[i].t[0] < 0)      /* deleted ? */
          continue;

        if(DT[i].p[0] == DPinfinity || DT[i].p[1] == DPinfinity || DT[i].p[2] == DPinfinity)
          continue;

        if(DT[i].p[3] == DPinfinity)
          continue;

        int j;
        for(j = 0; j < (NUMDIMS + 1); j++)
          {
            if(DP[DT[i].p[j]].task == ThisTask)
              if(DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < NumGas)
                if(TimeBinSynchronized[P[DP[DT[i].p[j]].index].TimeBinHydro])
                  break;
          }

        if(j == (NUMDIMS + 1))  /* this triangle does not a local point */
          continue;

        voronoi_test_evaluate(i, MODE_LOCAL_PARTICLES, thread_id);
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
        voronoi_test_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}



int voronoi_search_check(tessellation * TT)
{
  T = TT;
  int i, j;

  /* allocate buffers to arrange communication */
  NgbList = (int *) mymalloc_movable(&NgbList, "NgbList", T->Ndt * sizeof(int));


  generic_set_MaxNexport();

  generic_comm_pattern(T->Ndt, kernel_local, kernel_imported);


  for(i = 0; i < T->Ndt; i++)
    {
      tetra *DT = T->DT;
      point *DP = T->DP;

      if(DT[i].t[0] < 0)        /* deleted ? */
        continue;

      if(DT[i].p[0] == DPinfinity || DT[i].p[1] == DPinfinity || DT[i].p[2] == DPinfinity)
        continue;

      if(DT[i].p[3] == DPinfinity)
        continue;

      for(j = 0; j < (NUMDIMS + 1); j++)
        {
          if(DP[DT[i].p[j]].task == ThisTask)
            if(DP[DT[i].p[j]].index >= 0 && DP[DT[i].p[j]].index < NumGas)
              if(TimeBinSynchronized[P[DP[DT[i].p[j]].index].TimeBinHydro])
                break;
        }

      if(j == (NUMDIMS + 1))    /* this triangle does not a local point */
        continue;

      if(NgbList[i] != 4)
        {
          printf("TASK=%d: Triangle=%d NumNgb=%d\n", ThisTask, i, NgbList[i]);
        }
    }

  myfree(NgbList);

  return 0;
}

static int voronoi_test_evaluate(int target, int mode, int threadid)
{
  MyDouble *pos, h;
  int numnodes, *firstnode;
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
  h = 1.000001 * target_data->Rad;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  out.NumNgb = nfound;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif
