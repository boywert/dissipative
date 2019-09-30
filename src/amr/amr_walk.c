/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_walk.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"

#ifdef AMR
static inline unsigned long long ngb_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}


int amr_treefind_single_threads(MyDouble searchcenter[3], int target, int mode, int thread_id, int numnodes, int *firstnode)
{
  int numngb = 0;

  if(mode == MODE_IMPORTED_PARTICLES)
    {
      assert(numnodes == 1);
    }

  unsigned long long xxb = ngb_double_to_int(((searchcenter[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
  unsigned long long yyb = ngb_double_to_int(((searchcenter[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
  unsigned long long zzb = ngb_double_to_int(((searchcenter[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);
  unsigned long long mask = ((unsigned long long) 1) << (52 - 1);
  unsigned char shiftx = (52 - 1);
  unsigned char shifty = (52 - 2);
  unsigned char shiftz = (52 - 3);
  unsigned char levels = 0;

  int no = 0;
  while(TopNodes[no].Daughter >= 0)     /* walk down top tree to find correct leaf */
    {
      unsigned char subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

      mask >>= 1;
      levels++;

      no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
    }

  no = TopNodes[no].Leaf;

  int th = Ngb_DomainNodeIndex[no];

  if(mode == MODE_IMPORTED_PARTICLES)
    {
      assert(numnodes == 1);
      assert(th == firstnode[0]);
    }

  if(Ngb_DomainTask[th] != ThisTask)
    {
      int task = Ngb_DomainTask[th];

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
      nodelist[nexp].Node = th;

      Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
      return 0;
   }

  signed long long centermask = (0xFFF0000000000000llu) >> levels;

  unsigned char subnode = 0;

  while(1)
    {
      if(th < Ngb_MaxPart)      /* single particle */
        {
          Thread[thread_id].Ngblist[numngb++] = th;
          break;
        }
      else if(no < Ngb_MaxPart + Ngb_MaxNodes)  /* internal node */
        {
          subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

          centermask >>= 1;
          mask >>= 1;
          levels++;

          if(levels > MAX_TREE_LEVEL)
            {
              terminate("bad");
            }
          th = Ngb_Nodes[th].u.suns[subnode];
        }
      else
        {
          if(mode == MODE_IMPORTED_PARTICLES)
            terminate("mode == MODE_IMPORTED_PARTICLES should not occur here");
          terminate("should not happen");
        }

    }

  return numngb;
}
#endif
