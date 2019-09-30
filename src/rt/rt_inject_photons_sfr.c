/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_inject_photons_sfr.c
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

/* spread sfr sources */
#if defined(RT_SPREAD_SOURCE) && defined(USE_SFR) && defined(RT_HEALPIX_NSIDE)
static struct sourcedata_in
{
  MyDouble Pos[3];
  double Volume;
  int NodeList[NODELISTLENGTH];
}
 *SourceDataIn, *SourceDataGet;


void rt_inject_photons_spread(tessellation * T, double dt)
{
  int i, j, target, dummy;
  int ngrp, sendTask, recvTask, place, nexport, nimport, ndone, ndone_flag;
  double dphotons, hubble_a;
  MPI_Status status;

  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      hubble_a = 1.0;
    }


  CPU_Step[CPU_MISC] += measure_time();

  /* allocate buffers to arrange communication */


  Ngblist = (int *) mymalloc("NgbList", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct sourcedata_in) + sizeof(struct sourcedata_in) + sizemax(sizeof(struct sourcedata_in), sizeof(struct sourcedata_in))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  i = 0;

  do                            /* communication loop */
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      for(nexport = 0; i < NumGas; i++)
        {
          dphotons = SphP[i].Sfr * 1e53 * dt / hubble_a * All.UnitTime_in_s;

          if(SphP[i].Sfr > 0)
            if(rt_source_evaluate_sfr(i, 0, &nexport, Send_count, dphotons) < 0)
              break;
        }

      mysort_dataindex(DataIndexTable, nexport, sizeof(data_index), data_index_compare);

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
          nimport += Recv_count[j];

          if(j > 0)
            {
              Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
              Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }

      SourceDataGet = (struct sourcedata_in *) mymalloc("SourceDataGet", nimport * sizeof(struct sourcedata_in));
      SourceDataIn = (struct sourcedata_in *) mymalloc("SourceDataIn", nexport * sizeof(struct sourcedata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          SourceDataIn[j].Pos[0] = P[place].Pos[0];
          SourceDataIn[j].Pos[1] = P[place].Pos[1];
          SourceDataIn[j].Pos[2] = P[place].Pos[2];
          SourceDataIn[j].Volume = SphP[place].Volume;

          memcpy(SourceDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          sendTask = ThisTask;
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* get the particles */
                  MPI_Sendrecv(&SourceDataIn[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct sourcedata_in), MPI_BYTE,
                               recvTask, TAG_INJECT_A,
                               &SourceDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct sourcedata_in), MPI_BYTE, recvTask, TAG_INJECT_A, MPI_COMM_WORLD, &status);
                }
            }

        }

      myfree(SourceDataIn);

      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
        rt_source_evaluate_sfr(j, 1, &dummy, &dummy, dphotons);


      if(i < NumGas)
        ndone_flag = 0;
      else
        ndone_flag = 1;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


      myfree(SourceDataGet);

    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int rt_source_evaluate_sfr(int target, int mode, int *nexport, int *nsend_local, double dphotons)
{
  int i, k, j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3;
  double wk, dphot;
  double dx, dy, dz, r2, *pos, u, r;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = 4.0 * pow(3.0 / 4.0 / M_PI * SphP[target].Volume, 0.33);
    }
  else
    {
      pos = SourceDataGet[target].Pos;
      h = 4.0 * pow(3.0 / 4.0 / M_PI * SourceDataGet[target].Volume, 0.33);
    }


  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif


  if(mode == 0)
    {
      startnode = All.MaxPart;  /* root node */
    }
  else
    {
      startnode = SourceDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;        /* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb_inbox = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

          if(numngb_inbox < 0)
            return -1;

          for(n = 0; n < numngb_inbox; n++)
            {
              j = Ngblist[n];

              dx = pos[0] - P[j].Pos[0];
              dy = pos[1] - P[j].Pos[1];
              dz = pos[2] - P[j].Pos[2];

#ifdef SOURCE_PERIODIC          /*  now find the closest image in the given box size  */
              if(dx > boxHalf_X)
                dx -= boxSize_X;
              if(dx < -boxHalf_X)
                dx += boxSize_X;
              if(dy > boxHalf_Y)
                dy -= boxSize_Y;
              if(dy < -boxHalf_Y)
                dy += boxSize_Y;
              if(dz > boxHalf_Z)
                dz -= boxSize_Z;
              if(dz < -boxHalf_Z)
                dz += boxSize_Z;
#endif
              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  numngb++;

                  r = sqrt(r2);

                  u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                  dphot = P[j].Mass / SphP[j].Density * wk * dphotons;

                  for(i = 0; i < RT_N_DIR; i++)
                    {
                      SphP[j].Photons[i] += dphot / RT_N_DIR;
                      SphP[j].DensPhot[i] = SphP[j].Photons[i] / SphP[j].Volume;
                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = SourceDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Nodes[startnode].u.d.nextnode;      /* open it */
            }
        }

    }

  return 0;
}

#endif
#endif
