/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/domain_counttogo.c
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
#include <strings.h>
#include <math.h>


#include "allvars.h"
#include "proto.h"
#include "domain.h"
#include "voronoi.h"



/*! This function determines how many particles that are currently stored
 *  on the local CPU have to be moved off according to the domain
 *  decomposition.
 */
int domain_countToGo(void)
{
  for(int n = 0; n < NTask; n++)
    {
      toGo[n] = 0;
      toGoSph[n] = 0;
#ifdef TRACER_MC
      toGoTracer[n] = 0;
#endif
#ifdef GFM
      toGoStar[n] = 0;
#endif
#ifdef BLACK_HOLES
      toGoBHs[n] = 0;
#endif
#ifdef DUST_LIVE
      toGoDust[n] = 0;
#endif
    }

  for(int n = 0; n < NumPart; n++)
    {
      int no = 0;

      while(topNodes[no].Daughter >= 0)
        no = topNodes[no].Daughter + (Key[n] - topNodes[no].StartKey) / (topNodes[no].Size >> 3);

      no = topNodes[no].Leaf;

      if(DomainTask[no] != ThisTask)
        {
          toGo[DomainTask[no]] += 1;

          if(P[n].Type == 0)
            toGoSph[DomainTask[no]] += 1;
#ifdef GFM
          if(P[n].Type == 4)
            toGoStar[DomainTask[no]] += 1;
#endif
#ifdef BLACK_HOLES
          if(P[n].Type == 5)
            toGoBHs[DomainTask[no]] += 1;
#endif
#ifdef DUST_LIVE
          if(P[n].Type == DUST_LIVE)
            toGoDust[DomainTask[no]] += 1;
#endif

#ifdef TRACER_MC
          toGoTracer[DomainTask[no]] += get_number_of_tracers(n);
#endif
        }
    }

  MPI_Alltoall(toGo, 1, MPI_INT, toGet, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(toGoSph, 1, MPI_INT, toGetSph, 1, MPI_INT, MPI_COMM_WORLD);
#ifdef TRACER_MC
  MPI_Alltoall(toGoTracer, 1, MPI_INT, toGetTracer, 1, MPI_INT, MPI_COMM_WORLD);
#endif
#ifdef GFM
  MPI_Alltoall(toGoStar, 1, MPI_INT, toGetStar, 1, MPI_INT, MPI_COMM_WORLD);
#endif
#ifdef BLACK_HOLES
  MPI_Alltoall(toGoBHs, 1, MPI_INT, toGetBHs, 1, MPI_INT, MPI_COMM_WORLD);
#endif
#ifdef DUST_LIVE
  MPI_Alltoall(toGoDust, 1, MPI_INT, toGetDust, 1, MPI_INT, MPI_COMM_WORLD);
#endif

  return 0;
}
