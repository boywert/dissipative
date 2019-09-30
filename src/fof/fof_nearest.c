/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fof/fof_nearest.c
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
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_math.h>
#include <inttypes.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "fof.h"
#include "../subfind/subfind.h"

/*! \file fof.c
 *  \brief parallel FoF group finder
 */

#ifdef FOF


static MyFloat *fof_nearest_distance;
static MyFloat *fof_nearest_hsml;

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

static int fof_find_nearest_dmparticle_evaluate(int target, int mode, int threadid);


/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->Hsml = fof_nearest_hsml[i];

  in->Firstnode = firstnode;
}


 /* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat Distance;
  MyIDType MinID;
  int MinIDTask;
#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
  MyFloat DM_Hsml;
#endif
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(out->Distance < fof_nearest_distance[i])
    {
      fof_nearest_distance[i] = out->Distance;
      MinID[i] = out->MinID;
      MinIDTask[i] = out->MinIDTask;
#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      PS[i].Hsml = out->DM_Hsml;
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

  /* do local particles */
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
              if(generic_polling_primary(count, NumPart))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPart)
          break;

#ifdef TRACER_PARTICLE
        if(((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES)) && (P[i].Mass > 0 || P[i].Type == TRACER_PARTICLE) && P[i].ID > 0)
#else
        if((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES))
#endif
          {
            if(fof_nearest_distance[i] > 1.0e29)        /* we haven't found any neighbor yet */
              {
                fof_find_nearest_dmparticle_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
              }
          }
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

        fof_find_nearest_dmparticle_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}

double fof_find_nearest_dmparticle(MyIDType * vMinID, int *vHead, int *vLen, int *vNext, int *vTail, int *vMinIDTask)
{
  MinID = vMinID;
  Head = vHead;
  Len = vLen;
  Next = vNext;
  Tail = vTail;
  MinIDTask = vMinIDTask;

  int i, n, npleft, iter;
  long long ntot;
  double tstart = second();

  mpi_printf("FOF: Start finding nearest dm-particle (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));


  fof_nearest_distance = (MyFloat *) mymalloc("fof_nearest_distance", sizeof(MyFloat) * NumPart);
  fof_nearest_hsml = (MyFloat *) mymalloc("fof_nearest_hsml", sizeof(MyFloat) * NumPart);

  for(n = 0; n < NumPart; n++)
    {
#ifdef TRACER_PARTICLE
      if(((1 << P[n].Type) & (FOF_SECONDARY_LINK_TYPES)) && (P[n].Mass > 0 || P[n].Type == TRACER_PARTICLE) && P[n].ID > 0)
#else
      if((1 << P[n].Type) & (FOF_SECONDARY_LINK_TYPES))
#endif
        {
          fof_nearest_distance[n] = 1.0e30;
          if(P[n].Type == 0)
#ifdef USE_AREPO_FOF_WITH_GADGET_FIX
            fof_nearest_hsml[n] = SphP[n].Hsml;
#else
            fof_nearest_hsml[n] = get_cell_radius(n);
#endif
          else
            fof_nearest_hsml[n] = 0.1 * LinkL;
        }
    }

  generic_set_MaxNexport();

  iter = 0;
  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      double t0 = second();

      generic_comm_pattern(NumPart, kernel_local, kernel_imported);


      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
        {
#ifdef TRACER_PARTICLE
          if(((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES)) && (P[i].Mass > 0 || P[i].Type == TRACER_PARTICLE) && P[i].ID > 0)
#else
          if((1 << P[i].Type) & (FOF_SECONDARY_LINK_TYPES))
#endif
            {
              if(fof_nearest_distance[i] > 1.0e29)
                {
                  if(fof_nearest_hsml[i] < 4 * LinkL)   /* we only search out to a maximum distance */
                    {

                      /* need to redo this particle */
                      npleft++;
                      fof_nearest_hsml[i] *= 2.0;
                      if(iter >= MAXITER - 10)
                        {
                          printf("FOF: i=%d task=%d ID=%d P[i].Type=%d Hsml=%g LinkL=%g nearest=%g pos=(%g|%g|%g)\n",
                                 i, ThisTask, (int) P[i].ID, P[i].Type, fof_nearest_hsml[i], LinkL, fof_nearest_distance[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                          myflush(stdout);
                        }
                    }
                  else
                    {
                      fof_nearest_distance[i] = 0;      /* we do not continue to search for this particle */
                    }
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      double t1 = second();
      if(ntot > 0)
        {
          iter++;
          if(iter > 0)
            mpi_printf("FOF: fof-nearest iteration %d: need to repeat for %lld particles. (took = %g sec)\n", iter, ntot, timediff(t0, t1));

          if(iter > MAXITER)
            terminate("FOF: failed to converge in fof-nearest\n");
        }
    }
  while(ntot > 0);


  myfree(fof_nearest_hsml);
  myfree(fof_nearest_distance);

  mpi_printf("FOF: done finding nearest dm-particle\n");

  double tend = second();
  return timediff(tstart, tend);
}


static int fof_find_nearest_dmparticle_evaluate(int target, int mode, int threadid)
{
  int k, no, index, numnodes, *firstnode;
  double h, r2max, dist;
  double dx, dy, dz, r2;
  MyDouble *pos;
  data_in local, *target_data;
  data_out out;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

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
  h = target_data->Hsml;

  index = -1;
  r2max = 1.0e30;

  /* Now start the actual tree-walk computation for this particle */

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
              int p = no;
              no = Nextnode[no];

              if(!((1 << P[p].Type) & (FOF_SECONDARY_LINK_TARGET_TYPES)))
                continue;

              dist = h;
              dx = FOF_NEAREST_LONG_X(Tree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(Tree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(Tree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > dist)
                continue;

              r2 = dx * dx + dy * dy + dz * dz;
              if(r2 < r2max && r2 < h * h)
                {
                  index = p;
                  r2max = r2;
                }
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

              dist = h + 0.5 * current->len;
              dx = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
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
              terminate("do not expect imported points here");
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  if(index >= 0)
    {
      out.Distance = sqrt(r2max);
      out.MinID = MinID[Head[index]];
      out.MinIDTask = MinIDTask[Head[index]];
#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      out.DM_Hsml = PS[index].Hsml;
#endif
    }
  else
    {
      out.Distance = 2.0e30;
      out.MinID = 0;
      out.MinIDTask = -1;
#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==1)) || defined(GFM_WINDS_LOCAL)
      out.DM_Hsml = 0;
#endif
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif /* of FOF */
