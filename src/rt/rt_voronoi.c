/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_voronoi.c
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
#include "../voronoi.h"


#ifdef RT_ADVECT
void rt_voronoi_exchange_primitive_variables()
{
  CPU_Step[CPU_MISC] += measure_time();

  int listp;
  struct rt_primexch *rt_tmpPrimExch;
  int j, p, task, off;
  int ngrp, sendTask, recvTask, place;

  rt_tmpPrimExch = (struct rt_primexch *) mymalloc("rt_tmpPrimExch", Mesh_nexport * sizeof(struct rt_primexch));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  double sum;

  for(p = 0, sum = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;
          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  for(j = 0; j < RT_N_DIR; j++)
                    {
                      rt_tmpPrimExch[off].DensPhot[j] = SphP[place].DensPhot[j];
                      sum += rt_tmpPrimExch[off].DensPhot[j];
#ifndef RT_HEALPIX_NSIDE
                      rt_tmpPrimExch[off].SourceID[j] = SphP[place].SourceID[j];
                      rt_tmpPrimExch[off].SourcePos[j][0] = SphP[place].SourcePos[j][0];
                      rt_tmpPrimExch[off].SourcePos[j][1] = SphP[place].SourcePos[j][1];
                      rt_tmpPrimExch[off].SourcePos[j][2] = SphP[place].SourcePos[j][2];
#endif
                    }
                }
              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&rt_tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &RTPrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(rt_tmpPrimExch);

  CPU_Step[CPU_MESH_EXCHANGE] += measure_time();
}

void rt_voronoi_exchange_primitive_variables_and_gradients()
{
  CPU_Step[CPU_MISC] += measure_time();

  int listp;
  struct rt_grad_data *tmpGradExch;
  struct rt_primexch *tmpPrimExch;
  int j, p, task, off;
  int ngrp, sendTask, recvTask, place;

  tmpPrimExch = (struct rt_primexch *) mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct rt_primexch));
  tmpGradExch = (struct rt_grad_data *) mymalloc("tmpGradExch", Mesh_nexport * sizeof(struct rt_grad_data));

  /* prepare data for export */
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(p = 0; p < NumGas; p++)
    {
      if(P[p].Type == 0)
        {
          listp = List_P[p].firstexport;
          while(listp >= 0)
            {
              if((task = ListExports[listp].origin) != ThisTask)
                {
                  place = ListExports[listp].index;
                  off = Mesh_Send_offset[task] + Mesh_Send_count[task]++;

                  for(j = 0; j < RT_N_DIR; j++)
                    {
#ifndef RT_HEALPIX_NSIDE
                      tmpPrimExch[off].SourceID[j] = SphP[place].SourceID[j];
                      tmpPrimExch[off].SourcePos[j][0] = SphP[place].SourcePos[j][0];
                      tmpPrimExch[off].SourcePos[j][1] = SphP[place].SourcePos[j][1];
                      tmpPrimExch[off].SourcePos[j][2] = SphP[place].SourcePos[j][2];
#endif
                      tmpPrimExch[off].DensPhot[j] = SphP[place].DensPhot[j];

                      tmpGradExch[off].ddensphot[j][0] = SphP[place].rt_Grad.ddensphot[j][0];
                      tmpGradExch[off].ddensphot[j][1] = SphP[place].rt_Grad.ddensphot[j][1];
                      tmpGradExch[off].ddensphot[j][2] = SphP[place].rt_Grad.ddensphot[j][2];
                    }

#ifdef RT_HEALPIX_NSIDE
                  tmpGradExch[off].ddensphot_unlimited[0] = SphP[place].rt_Grad.ddensphot_unlimited[0];
                  tmpGradExch[off].ddensphot_unlimited[1] = SphP[place].rt_Grad.ddensphot_unlimited[1];
                  tmpGradExch[off].ddensphot_unlimited[2] = SphP[place].rt_Grad.ddensphot_unlimited[2];
#endif
                }
              listp = ListExports[listp].nextexport;
            }
        }
    }

  /* exchange data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &RTPrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct rt_grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &RTGradExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct rt_grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpGradExch);
  myfree(tmpPrimExch);

  CPU_Step[CPU_MESH_EXCHANGE] += measure_time();
}

#endif
