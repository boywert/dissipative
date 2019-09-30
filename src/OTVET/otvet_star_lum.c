/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/OTVET/otvet_star_lum.c
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
#include "otvet_proto.h"


#ifdef OTVET

static struct stardata_in
{
  MyDouble Pos[3], Mass;
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
 *StarDataIn, *StarDataGet;

void otvet_star_lum(void)
{
  int j;
  int idx, i, dummy;
  int ngrp, recvTask, place, nexport, nimport, ndone, ndone_flag;

  /* clear Je in all gas particles */

  for(j = 0; j < NumGas; j++)
    if(P[j].Type == 0)
      for(i = 0; i < OT_N_BINS; i++)
        SphP[j].Je[i] = 0;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) + 2 * sizeof(struct stardata_in)));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  idx = 0;                      /* beginn with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      for(nexport = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(P[i].Type == 4)
            {
              if(star_lum_evaluate(i, 0, &nexport, Send_count) < 0)
                break;
            }
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

      StarDataGet = (struct stardata_in *) mymalloc("StarDataGet", nimport * sizeof(struct stardata_in));
      StarDataIn = (struct stardata_in *) mymalloc("StarDataIn", nexport * sizeof(struct stardata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          StarDataIn[j].Pos[0] = P[place].Pos[0];
          StarDataIn[j].Pos[1] = P[place].Pos[1];
          StarDataIn[j].Pos[2] = P[place].Pos[2];
          StarDataIn[j].Mass = P[place].Mass;
#if defined(EDDINGTON_TENSOR_STARS) && defined(OTVET_SCATTER_SOURCE)
          StarDataIn[j].Hsml = P[place].OtvetHsml;
#else
          StarDataIn[j].Hsml = All.FoeceSoftening[P[place].SofteningType];
#endif

          memcpy(StarDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

        }

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  /* get the particles */
                  MPI_Sendrecv(&StarDataIn[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct stardata_in), MPI_BYTE,
                               recvTask, TAG_DENS_A,
                               &StarDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct stardata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(StarDataIn);


      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
        star_lum_evaluate(j, 1, &dummy, &dummy);

      /* check whether this is the last iteration */
      if(idx == TimeBinsGravity.NActiveParticles)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      myfree(StarDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}


int star_lum_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int i, j, n, numngb;
  int startnode, listindex = 0;
  double h, hinv, h2, hinv3;
  double wk, mass, fac;
  double dx, dy, dz, r, r2, u;
  MyDouble *pos;
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif

  if(mode == 0)
    {
      pos = P[target].Pos;
#if defined(EDDINGTON_TENSOR_STARS) && defined(OTVET_SCATTER_SOURCE)
      h = P[target].OtvetHsml;
#else
      h = All.ForceSoftening[P[target].SofteningType];
#endif
      mass = P[target].Mass;
    }
  else
    {
      pos = StarDataGet[target].Pos;
      h = StarDataGet[target].Hsml;
      mass = StarDataGet[target].Mass;
    }

  h2 = h * h;
  hinv = 1.0 / h;
  hinv3 = hinv * hinv * hinv;

  fac = mass * All.UnitMass_in_g / SOLAR_MASS;  /* LVS: for RADPRESS_THIN we have All.HubbleParam multiplying, check! */
#ifndef OTVET_MULTI_FREQUENCY
  lum[0] = All.IonizingLumPerSolarMass * All.UnitTime_in_s / All.HubbleParam;
//  lum[0] = fac * All.IonizingLumPerSolarMass * All.UnitTime_in_s / All.HubbleParam;   //LVS: original in Gadget, but "fac" is used again in Je??
#endif


  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = StarDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

          if(numngb < 0)
            return -1;

          for(n = 0; n < numngb; n++)
            {
              j = Ngblist[n];

              dx = GRAVITY_NEAREST_X(pos[0] - P[j].Pos[0]);
              dy = GRAVITY_NEAREST_Y(pos[1] - P[j].Pos[1]);
              dz = GRAVITY_NEAREST_Z(pos[2] - P[j].Pos[2]);
              r2 = dx * dx + dy * dy + dz * dz;
              r = sqrt(r2);

              if(r2 < h2)
                {
                  u = r * hinv;

                  if(u < 0.5)
                    wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                  else
                    wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

                }
              else
                wk = 0;

              for(i = 0; i < OT_N_BINS; i++)
                SphP[j].Je[i] += lum[i] * wk * fac;     /* Je is in "photon s^-1 */
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = StarDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  return 0;
}

#endif
