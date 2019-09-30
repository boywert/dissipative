/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/subfind/subfind_nearesttwo.c
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

#ifdef SUBFIND
#include "subfind.h"


static int subfind_nearesttwo_evaluate(int target, int mode, int threadid);

/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyIDType ID;
  MyFloat Hsml;
  MyFloat Density;
  MyFloat Dist[2];
  int Count;
  long long Index[2];

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;


 /* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in * in, int i, int firstnode)
{
  int k;

#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = PS[i].Center[0];
      in->Pos[1] = PS[i].Center[1];
      in->Pos[2] = PS[i].Center[2];
    }
  else
#endif
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }

  in->Hsml = PS[i].Hsml;
  in->ID = P[i].ID;
  in->Density = PS[i].Density;
  in->Count = NgbLoc[i].count;
  for(k = 0; k < NgbLoc[i].count; k++)
    {
      in->Dist[k] = R2Loc[i].dist[k];
      in->Index[k] = NgbLoc[i].index[k];
    }
  in->Firstnode = firstnode;
}


typedef struct
{
  MyFloat Dist[2];
  long long Index[2];
  int Count;
} data_out;

static data_out *DataResult, *DataOut;


 /* routine to store or combine result data */
static void out2particle(data_out * out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
    {
      int k;

      NgbLoc[i].count = out->Count;

      for(k = 0; k < out->Count; k++)
        {
          R2Loc[i].dist[k] = out->Dist[k];
          NgbLoc[i].index[k] = out->Index[k];
        }
    }
  else                          /* combine */
    {
      int k, l;

      for(k = 0; k < out->Count; k++)
        {
          if(NgbLoc[i].count >= 1)
            if(NgbLoc[i].index[0] == out->Index[k])
              continue;

          if(NgbLoc[i].count == 2)
            if(NgbLoc[i].index[1] == out->Index[k])
              continue;

          if(NgbLoc[i].count < 2)
            {
              l = NgbLoc[i].count;
              NgbLoc[i].count++;
            }
          else
            {
              if(R2Loc[i].dist[0] > R2Loc[i].dist[1])
                l = 0;
              else
                l = 1;

              if(out->Dist[k] >= R2Loc[i].dist[l])
                continue;
            }

          R2Loc[i].dist[l] = out->Dist[k];
          NgbLoc[i].index[l] = out->Index[k];

          if(NgbLoc[i].count == 2)
            if(NgbLoc[i].index[0] == NgbLoc[i].index[1])
              terminate("this is not supposed to happen");
        }
    }
}

#define USE_SUBCOMM_COMMUNICATOR
#include "../generic_comm_helpers2.h"

static double *Dist2list;
static int *Ngblist;


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

    for(j = 0; j < SubNTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, NumPartGroup))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= NumPartGroup)
          break;

        subfind_nearesttwo_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        subfind_nearesttwo_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }

}

void subfind_find_nearesttwo(void)
{
  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Start finding nearest two.\n", ThisTask);

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPartGroup * sizeof(int));
  Dist2list = (double *) mymalloc("Dist2list", NumPartGroup * sizeof(double));

  generic_set_MaxNexport();


  for(int i = 0; i < NumPartGroup; i++)
    NgbLoc[i].count = 0;


  generic_comm_pattern(NumPartGroup, kernel_local, kernel_imported);


  myfree(Dist2list);
  myfree(Ngblist);

  if(SubThisTask == 0)
    printf("SUBFIND-COLLECTIVE, root-task=%d: Done with nearest two.\n", ThisTask);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int subfind_nearesttwo_evaluate(int target, int mode, int threadid)
{
  int j, k, n, no, count;
  MyIDType ID;
  long long index[2];
  double dist[2];
  int numngb, numnodes, *firstnode;
  double hsml;
  double density;
  MyDouble *pos;
  struct NODE *current;
  double dx, dy, dz, disthsml, r2;
#ifdef PERIODIC
  MyDouble xtmp, ytmp, ztmp;
#endif

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

  ID = in->ID;
  density = in->Density;
  pos = in->Pos;
  hsml = in->Hsml;
  count = in->Count;
  for(k = 0; k < count; k++)
    {
      dist[k] = in->Dist[k];
      index[k] = in->Index[k];
    }

  if(count == 2)
    if(index[0] == index[1])
      {
        terminate("task=%d/%d target=%d mode=%d  index_0=%lld  index_1=%lld\n", SubThisTask, ThisTask, target, mode, index[0], index[1]);
      }

  numngb = 0;
  count = 0;

  hsml *= 1.00001;              /* prevents that the most distant neighbour on the edge of the search region may not be found.
                                 * (needed for consistency with serial algorithm)
                                 */

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = SubTree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = SubNodes[no].u.d.nextnode;       /* open it */
        }
      while(no >= 0)
        {
          if(no < SubTree_MaxPart)      /* single particle */
            {
              int p = no;
              no = SubNextnode[no];

              disthsml = hsml;
              dx = FOF_NEAREST_LONG_X(SubTree_Pos_list[3 * p + 0] - pos[0]);
              if(dx > disthsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(SubTree_Pos_list[3 * p + 1] - pos[1]);
              if(dy > disthsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(SubTree_Pos_list[3 * p + 2] - pos[2]);
              if(dz > disthsml)
                continue;
              if((r2 = (dx * dx + dy * dy + dz * dz)) > disthsml * disthsml)
                continue;

              Dist2list[numngb] = r2;
              Ngblist[numngb++] = p;
            }
          else if(no < SubTree_MaxPart + SubTree_MaxNodes)      /* internal node */
            {
              if(mode == 1)
                {
                  if(no < SubTree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      break;
                    }
                }

              current = &SubNodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              disthsml = hsml + 0.5 * current->len;

              dx = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > disthsml)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > disthsml)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > disthsml)
                continue;
              /* now test against the minimal sphere enclosing everything */
              disthsml += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > disthsml * disthsml)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }
          else if(no >= SubTree_ImportedNodeOffset)     /* point from imported nodelist */
            {
              terminate("do not expect imported points here");
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)   /* note: if no target is given, export will not occur */
                subfind_treefind_collective_export_node_threads(no, target, threadid);

              no = SubNextnode[no - SubTree_MaxNodes];
            }
        }
    }

  for(n = 0; n < numngb; n++)
    {
      j = Ngblist[n];
      r2 = Dist2list[n];

      if(P[j].ID != ID)         /* exclude the self-particle */
        {
          if(PS[j].Density > density)   /* we only look at neighbours that are denser */
            {
              if(count < 2)
                {
                  dist[count] = r2;
                  index[count] = (((long long) SubThisTask) << 32) + j;
                  count++;
                }
              else
                {
                  if(dist[0] > dist[1])
                    k = 0;
                  else
                    k = 1;

                  if(r2 < dist[k])
                    {
                      dist[k] = r2;
                      index[k] = (((long long) SubThisTask) << 32) + j;
                    }
                }
            }
        }
    }

  out.Count = count;
  for(k = 0; k < count; k++)
    {
      out.Dist[k] = dist[k];
      out.Index[k] = index[k];
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}



#endif
