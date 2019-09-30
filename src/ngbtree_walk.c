/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/ngbtree_walk.c
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
#include <time.h>

#include "allvars.h"
#include "proto.h"



/*! This function returns the number of neighbours with distance <=
    hsml, and returns the particle indices in the global buffer
    Ngblist. (Actually, particles in a box of half side length hsml
    are returned, i.e. the reduction to a sphere still needs to be
    done in the calling routine.) The return value is the number of
    neighbours found. The tree traversal starts at startnode.
    Target seems to be the index of the particle at the center.
    See ngb_treefind_export_node for definitions of nexport and nsend_local.
    Don't know what mode is...
 */
int ngb_treefind_variable_threads(MyDouble searchcenter[3], MyFloat hsml, int target, int mode, int thread_id, int numnodes, int *firstnode)
{
  MyDouble search_min[3], search_max[3], search_max_Lsub[3], search_min_Ladd[3];

  for(int i = 0; i < 3; i++)
    {
      search_min[i] = searchcenter[i] - 1.001 * hsml;
      search_max[i] = searchcenter[i] + 1.001 * hsml;
    }

  search_max_Lsub[0] = search_max[0] - boxSize_X;
  search_max_Lsub[1] = search_max[1] - boxSize_Y;
  search_max_Lsub[2] = search_max[2] - boxSize_Z;

  search_min_Ladd[0] = search_min[0] + boxSize_X;
  search_min_Ladd[1] = search_min[1] + boxSize_Y;
  search_min_Ladd[2] = search_min[2] + boxSize_Z;

  int numngb = 0;
  double xtmp, ytmp, ztmp;
  double hsml2 = hsml * hsml;

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Ngb_MaxPart;     /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Ngb_Nodes[no].u.d.nextnode;      /* open it */
        }

      while(no >= 0)
        {
          if(no < Ngb_MaxPart)  /* single particle */
            {
              int p = no;
              no = Ngb_Nextnode[no];

              if(P[p].Type > 0)
                continue;

              if(P[p].Ti_Current != All.Ti_Current)
                {
#if (NUM_THREADS > 1)
                  omp_set_lock(&ParticleLocks[p]);

                  if(P[p].Ti_Current != All.Ti_Current)
                    {
#endif
                      drift_particle(p, All.Ti_Current);
#if (NUM_THREADS > 1)
                    }
                  omp_unset_lock(&ParticleLocks[p]);
#endif
                }

              double dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
              if(dx > hsml)
                continue;
              double dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
              if(dy > hsml)
                continue;
              double dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
              if(dz > hsml)
                continue;

              double r2 = dx * dx + dy * dy + dz * dz;
              if(r2 > hsml2)
                continue;

              Thread[thread_id].R2list[numngb] = r2;
              Thread[thread_id].Ngblist[numngb++] = p;
            }
          else if(no < Ngb_MaxPart + Ngb_MaxNodes)      /* internal node */
            {
              struct NgbNODE *current = &Ngb_Nodes[no];

              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no < Ngb_FirstNonTopLevelNode)     /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

#if (NUM_THREADS > 1)
              int no_old = no;
#endif

              no = current->u.d.sibling;        /* in case the node can be discarded */

              if(current->Ti_Current != All.Ti_Current)
                {
#if (NUM_THREADS > 1)
                  omp_set_lock(&Ngb_NodeLocks[no_old]);

                  if(current->Ti_Current != All.Ti_Current)
                    {
#endif
                      drift_node(current, All.Ti_Current);

#if (NUM_THREADS > 1)
                    }
                  omp_unset_lock(&Ngb_NodeLocks[no_old]);
#endif
                }

              if(search_min[0] > current->u.d.range_max[0] && search_max_Lsub[0] < current->u.d.range_min[0])
                continue;
              if(search_min_Ladd[0] > current->u.d.range_max[0] && search_max[0] < current->u.d.range_min[0])
                continue;

              if(search_min[1] > current->u.d.range_max[1] && search_max_Lsub[1] < current->u.d.range_min[1])
                continue;
              if(search_min_Ladd[1] > current->u.d.range_max[1] && search_max[1] < current->u.d.range_min[1])
                continue;

              if(search_min[2] > current->u.d.range_max[2] && search_max_Lsub[2] < current->u.d.range_min[2])
                continue;
              if(search_min_Ladd[2] > current->u.d.range_max[2] && search_max[2] < current->u.d.range_min[2])
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }
          else                  /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES should not occur here");

              if(target >= 0)   /* if no target is given, export will not occur */
                if(ngb_treefind_export_node_threads(no, target, thread_id, 0))
                  return -1;

              no = Ngb_Nextnode[no - Ngb_MaxNodes];
              continue;
            }

        }
    }
  return numngb;
}



/** Exports node when a pseudoparticle node no is found during the
    tree traversal in ngb_treefind_variable. target is the index
    written into the export arrays to that it can be deduced what it
    refers to. Nsend_local is the send count array that indicates how
    many items to send to the respective task. It is increased by one
    for the task the pseudoparticle is on as the item is added to the
    DataIndexTable.

    nexport is the total number of items to send, and nexport_save is the
    value of nexport to revert to if we run out of space in the send buffer.

    A return code -1 indicates that the exchange buffers are full and
    we need to do an exchange before continuing.

    The function writes to the global buffers Exportflag,
    Exportnodecount, Exportindex, DataIndexTable, DataNodeList.
  */
int ngb_treefind_export_node_threads(int no, int target, int thread_id, int image_flag)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != target)
    {
      Thread[thread_id].Exportflag[task] = target;
      int nexp = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task = task;
      Thread[thread_id].PartList[nexp].Index = target;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp = Thread[thread_id].NexportNodes++;
  nexp = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *) (((char *) Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task = task;
  nodelist[nexp].Index = target;
  nodelist[nexp].Node = Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];
#ifdef EXTENDED_GHOST_SEARCH
  nodelist[nexp].BitFlags = image_flag;
#endif
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}
