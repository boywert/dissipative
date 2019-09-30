/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_disk_vorticity.c
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


#if defined(BLACK_HOLES) && defined(BH_BONDI_DISK_VORTICITY)

static struct blackholedata_in
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  int NodeList[NODELISTLENGTH];
}
 *BlackholeDataIn, *BlackholeDataGet;

static struct blackholedata_out
{
  MyFloat Gal_Mass;
  MyFloat GasVort[3];
}
 *BlackholeDataResult, *BlackholeDataOut;

static int blackhole_evaluate_vorticity(int target, int mode, int *nexport, int *nSend_local);

void blackhole_disk_vorticity(void)
{
  int idx, i, j, k, n, nexport, nimport, ndone, ndone_flag;
  int ngrp, recvTask, dummy;

  mpi_printf("BLACK_HOLES: Begin disk vorticity evaluation for each BH.\n");

  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct blackholedata_in) + sizeof(struct blackholedata_out) + sizemax(sizeof(struct blackholedata_in), sizeof(struct blackholedata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));



  idx = 0;                      /* firstparticle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */

      for(nexport = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
        {
          i = TimeBinsBHAccretion.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(blackhole_evaluate_vorticity(i, 0, &nexport, Send_count) < 0)
            break;
        }


      mysort(DataIndexTable, nexport, sizeof(data_index), data_index_compare);

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

      BlackholeDataGet = (struct blackholedata_in *) mymalloc("BlackholeDataGet", nimport * sizeof(struct blackholedata_in));
      BlackholeDataIn = (struct blackholedata_in *) mymalloc("BlackholeDataIn", nexport * sizeof(struct blackholedata_in));

      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          for(k = 0; k < 3; k++)
            BlackholeDataIn[j].Pos[k] = P[place].Pos[k];

          BlackholeDataIn[j].BH_Hsml = BPP(place).BH_Hsml;

          memcpy(BlackholeDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              MPI_Sendrecv(&BlackholeDataIn[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &BlackholeDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct blackholedata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      myfree(BlackholeDataIn);
      BlackholeDataResult = (struct blackholedata_out *) mymalloc("BlackholeDataResult", nimport * sizeof(struct blackholedata_out));
      BlackholeDataOut = (struct blackholedata_out *) mymalloc("BlackholeDataOut", nexport * sizeof(struct blackholedata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        blackhole_evaluate_vorticity(j, 1, &dummy, &dummy);

      if(idx == TimeBinsBHAccretion.NActiveParticles)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              MPI_Sendrecv(&BlackholeDataResult[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct blackholedata_out),
                           MPI_BYTE, recvTask, TAG_DENS_B,
                           &BlackholeDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct blackholedata_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;
          BPP(place).Gal_Mass += BlackholeDataOut[j].Gal_Mass;
          BPP(place).BH_GasVort[0] += BlackholeDataOut[j].GasVort[0];
          BPP(place).BH_GasVort[1] += BlackholeDataOut[j].GasVort[1];
          BPP(place).BH_GasVort[2] += BlackholeDataOut[j].GasVort[2];
        }

      myfree(BlackholeDataOut);
      myfree(BlackholeDataResult);
      myfree(BlackholeDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


  for(idx = 0; idx < TimeBinsBHAccretion.NActiveParticles; idx++)
    {
      n = TimeBinsBHAccretion.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(BPP(n).Gal_Mass > 0)
        {
          BPP(n).BH_GasVort[0] /= BPP(n).Gal_Mass;
          BPP(n).BH_GasVort[1] /= BPP(n).Gal_Mass;
          BPP(n).BH_GasVort[2] /= BPP(n).Gal_Mass;
        }
    }

  mpi_printf("\nBLACK_HOLES: Disk vorticity done.\n");
}




static int blackhole_evaluate_vorticity(int target, int mode, int *nexport, int *nSend_local)
{
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  int startnode, numngb, j, k, n, listindex = 0;
  double dx, dy, dz, r2;
  MyDouble *pos;
  double h;
  MyFloat galmass = 0;
  MyFloat gasvort[3];
  gasvort[0] = gasvort[1] = gasvort[2] = 0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = BPP(target).BH_Hsml;
    }
  else
    {
      pos = BlackholeDataGet[target].Pos;
      h = BlackholeDataGet[target].BH_Hsml;
    }

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = BlackholeDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          if(All.DiskVorticityRadius == 0)
            numngb = ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nSend_local);
          else
            numngb = ngb_treefind_variable(pos, All.DiskVorticityRadius, target, &startnode, mode, nexport, nSend_local);

          if(numngb < 0)
            return -1;

          for(n = 0; n < numngb; n++)
            {
              j = Ngblist[n];

              if(P[j].Type == 0)
                {
                  dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
                  dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
                  dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

                  r2 = dx * dx + dy * dy + dz * dz;

                  if(All.DiskVorticityRadius == 0)
                    if(r2 < h * h)
                      {
                        galmass += P[j].Mass;
                        gasvort[0] += (P[j].Mass * (SphP[j].Grad.dvel[2][1] - SphP[j].Grad.dvel[1][2]));
                        gasvort[1] += (P[j].Mass * (SphP[j].Grad.dvel[0][2] - SphP[j].Grad.dvel[2][0]));
                        gasvort[2] += (P[j].Mass * (SphP[j].Grad.dvel[1][0] - SphP[j].Grad.dvel[0][1]));
                      }

                  if(All.DiskVorticityRadius != 0)
                    if(r2 < All.DiskVorticityRadius * All.DiskVorticityRadius)
                      {
                        galmass += P[j].Mass;
                        gasvort[0] += (P[j].Mass * (SphP[j].Grad.dvel[2][1] - SphP[j].Grad.dvel[1][2]));
                        gasvort[1] += (P[j].Mass * (SphP[j].Grad.dvel[0][2] - SphP[j].Grad.dvel[2][0]));
                        gasvort[2] += (P[j].Mass * (SphP[j].Grad.dvel[1][0] - SphP[j].Grad.dvel[0][1]));
                      }
                }
            }
        }
      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = BlackholeDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      BPP(target).Gal_Mass = galmass;
      for(k = 0; k < 3; k++)
        BPP(target).BH_GasVort[k] = gasvort[k];
    }
  else
    {
      BlackholeDataResult[target].Gal_Mass = galmass;
      for(k = 0; k < 3; k++)
        BlackholeDataResult[target].GasVort[k] = gasvort[k];
    }

  return 0;
}


#endif
