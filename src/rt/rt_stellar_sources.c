/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_stellar_sources.c
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

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef RT_ADVECT
#ifdef RT_STELLAR_SOURCES

static struct source_list
{
  double pos[3], lum;
}
 *SourceList, *SourceListGet;

void rt_create_source_list()
{
  int i, j, nsrc, nimport, ngrp;

  SourceList = mymalloc("SourceList", N_SOURCES * sizeof(struct source_list));

  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if(P[i].Type == 4)
        {
          SourceList[nsrc].pos[0] = P[i].Pos[0];
          SourceList[nsrc].pos[1] = P[i].Pos[1];
          SourceList[nsrc].pos[2] = P[i].Pos[2];

          SourceList[nsrc++].lum = P[i].Mass * 1e52;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  SourceListGet = mymalloc("SourceListGet", nimport * sizeof(struct source_list));

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SourceList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct source_list), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &SourceListGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct source_list), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  mysort(SourceListGet, nimport, sizeof(struct source_list), source_list_lum_compare);

  for(i = 0; i < nimport; i++)
    {
      Source_Pos[i][0] = SourceListGet[i].pos[0];
      Source_Pos[i][1] = SourceListGet[i].pos[1];
      Source_Pos[i][2] = SourceListGet[i].pos[2];

      Source_Lum[i] = SourceListGet[i].lum;
      Source_ID[i] = i;

      if(ThisTask == 0)
        printf("%d %g \n", Source_ID[i], Source_Lum[i]);
    }

  myfree(SourceListGet);
  myfree(SourceList);

}

int source_list_lum_compare(const void *a, const void *b)
{
  if(((struct source_list *) a)->lum < (((struct source_list *) b)->lum))
    return -1;

  if(((struct source_list *) a)->lum > (((struct source_list *) b)->lum))
    return +1;

  return 0;
}

#endif
#endif
