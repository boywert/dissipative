/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_neighbors.c
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

#ifdef DUST_LIVE
#ifdef DL_GRAIN_BINS
/* DL_SNE_DESTRUCTION uses dust density calculations to possibly enhance the
 * supernova rate at high density. */
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)

#if defined(HIERARCHICAL_GRAVITY) || defined(ALLOW_DIRECT_SUMMATION)
#error "Use of options DL_SNE_DESTRUCTION and DL_SHATTERING and DL_COAGULATION does not work with HIERARCHICAL_GRAVITY or ALLOW_DIRECT_SUMMATION since we need a full gravity tree for neighbor searches at each time!"
#endif

typedef struct
{
  /* your fields go here */
  MyDouble Pos[3];
  MyFloat Hsml;
  //double Mass;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

typedef struct
{
  /* your fields go here */
  MyFloat DustDensity;
  MyFloat DustNumNgb;
} data_out;

static data_out *DataResult, *DataOut;

static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

  int idx = TargetList[i];

  if(idx < NumPart)
    {
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = P[idx].Pos[k];
        }
    }
  else
    {
      idx -= Tree_ImportedNodeOffset;
      for(k = 0; k < 3; k++)
        {
          in->Pos[k] = Tree_Points[idx].Pos[k];
        }
    }
  in->Hsml = PShatter[i].DustHsml;
  in->Firstnode = firstnode;
}

static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      PShatter[i].DustDensity = out->DustDensity;
      PShatter[i].DustNumNgb = out->DustNumNgb;
    }
  else                          /* combine */
    {
      PShatter[i].DustDensity += out->DustDensity;
      PShatter[i].DustNumNgb += out->DustNumNgb;
    }
}

#include "../generic_comm_helpers2.h"

static unsigned char *Todo;

static void kernel_local(void)
{
  int idx;
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
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= Nforces)
          break;

        if(Todo[idx] > 0)
          dust_findHsml_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
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

        dust_findHsml_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void dust_findHsml(void)
{
  int idx, i;
  int iter = 0, npleft;
  MyDouble *Left, *Right;
  double t0, t1;
  long long ntot, npartall;

  sumup_large_ints(1, &Nforces, &npartall);
  if(npartall == 0)
    return;

  mpi_printf("DUST_LIVE: Finding dust-dust hsml values\n");

  /* Nforces is set in begin_shattering() */
  Left = mymalloc("Left", Nforces * sizeof(MyDouble));
  Right = mymalloc("Right", Nforces * sizeof(MyDouble));
  Todo = mymalloc("Todo", Nforces * sizeof(unsigned char));

  int nforces = 0;
  for(idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == DUST_LIVE) && (Tree_Task_list[i] == ThisTask) && (P[i].Mass > 0.0))
        {
          Left[nforces] = 0;
          Right[nforces] = 0;
          Todo[nforces] = 1;
          PShatter[nforces].DustHsml = DTP(i).DustHsml;
          /* Reset dust smoothing length if previously there were not enough
           * dust particles. */
          if(PShatter[nforces].DustHsml >= 0.1*All.BoxSize)
            PShatter[nforces].DustHsml = get_default_softening_of_particletype(DUST_LIVE);
          nforces++;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if((Tree_Points[i].Type == DUST_LIVE) && (Tree_Points[i].Mass > 0.0))
      {
        Left[nforces] = 0;
        Right[nforces] = 0;
        Todo[nforces] = 1;
        PShatter[nforces].DustHsml = Tree_Points[i].DustHsml;
        /* Reset dust smoothing length if previously there were not enough
         * dust particles. */
        if(PShatter[nforces].DustHsml >= 0.1*All.BoxSize)
          PShatter[nforces].DustHsml = get_default_softening_of_particletype(DUST_LIVE);
        nforces++;
      }

  generic_set_MaxNexport();

  do
    {
      t0 = second();

      generic_comm_pattern(Nforces, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < Nforces; i++)
        {
          if(Todo[i])
            {
              if(((PShatter[i].DustNumNgb < (All.DesNumNgbDust - All.MaxNumNgbDeviationDust)) || (PShatter[i].DustNumNgb > (All.DesNumNgbDust + All.MaxNumNgbDeviationDust))) && (PShatter[i].DustHsml < All.BoxSize/10.0))
                {
                  /* need to redo this particle */
                  npleft++;

                  if(PShatter[i].DustNumNgb < All.DesNumNgbDust - All.MaxNumNgbDeviationDust)
                    Left[i] = dmax(PShatter[i].DustHsml, Left[i]);
                  else
                    {
                      if(Right[i] != 0)
                        {
                          if(PShatter[i].DustHsml < Right[i])
                            Right[i] = PShatter[i].DustHsml;
                        }
                      else
                        Right[i] = PShatter[i].DustHsml;
                    }

                  if(iter >= MAXITER - 10)
                    {
                      printf("DUST_LIVE: i=%d task=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n", i, ThisTask, PShatter[i].DustHsml, Left[i], Right[i], (double) PShatter[i].DustNumNgb, Right[i] - Left[i]);
                      fflush(stdout);
                    }

                  if(Right[i] > 0 && Left[i] > 0)
                    PShatter[i].DustHsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                    {
                      if(Right[i] == 0 && Left[i] == 0)
                        terminate("BAD");       /* can't occur */

                      if(Right[i] == 0 && Left[i] > 0)
                        PShatter[i].DustHsml *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                        PShatter[i].DustHsml /= 1.26;
                    }
                }
              else
                {
                  Todo[i] = 0;
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0 && ThisTask == 0)
            {
              printf("DUST_LIVE: ngb iteration %d: need to repeat for %llu particles. (took %g sec)\n", iter, (unsigned long long) ntot, timediff(t0, t1));
              fflush(stdout);
            }

          if(iter > MAXITER)
            {
              printf("DUST_LIVE: failed to converge in dust-dust neighbour iteration\n");
              fflush(stdout);
              terminate("BAD");
            }
        }
    }
  while(ntot > 0);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

  mpi_printf("DUST_LIVE: Done with dust-dust hsml values\n");
}

void dust_findHsml_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;
  data_in local, *in;
  data_out out;
  int k;
  int no;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double wk, h, h2, h3, hinv, hinv3;
  double r, r2, u;
  MyDouble *pos;
  MyDouble dx, dy, dz;

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

  memset(&out, 0, sizeof(data_out));

  pos = in->Pos;
  h = in->Hsml;
  h2 = h * h;
  h3 = h2 * h;
  hinv = 1. / h;
  hinv3 = hinv * hinv * hinv;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart;    /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode;  /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if((r2 < h2) && (P[no].Type == DUST_LIVE) && (P[no].Mass > 0.0))
                {
                  //hinv = 1. / h;
                  //hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);
                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }

                  out.DustDensity += P[no].Mass * wk;
                  out.DustNumNgb += NORM_COEFF * wk * h3;

                }
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];
              no = current->u.d.sibling;        /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }

          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if((r2 < h2) && (Tree_Points[n].Type == DUST_LIVE) && (Tree_Points[n].Mass > 0.0))
                {
                  //hinv = 1. / h;
                  //hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);
                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }

                  out.DustDensity += Tree_Points[n].Mass * wk;
                  out.DustNumNgb += NORM_COEFF * wk * h3;
                }

              no = Nextnode[no - Tree_MaxNodes];

            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }
    }

  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
}

#endif
#endif
#endif
