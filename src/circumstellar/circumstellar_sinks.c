/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/circumstellar/circumstellar_sinks.c
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

#if defined(CIRCUMSTELLAR) && defined(CIRCUMSTELLAR_SINKS)


static struct circumstellardata_in
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble Mass;
  MyFloat Radius;

  /** following needed for bookkeeping */
  int Index;
  int NodeList[NODELISTLENGTH];
}
 *CircumstellarDataIn, *CircumstellarDataGet;

static struct circumstellardata_out
{
  MyFloat AccretedMass;
  MyFloat AccretedMomentum[3];
}
 *CircumstellarDataResult, *CircumstellarDataOut;

static int NumGas_swallowed, Ntot_gas_swallowed;

static int circumstellar_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local);

static double hubble_a, ascale;



void circumstellar_swallow_gas(void)
{
  int idx, i, j, k, nexport, nimport, ndone, ndone_flag;
  int ngrp, recvTask, dummy;

  if(All.Time == All.TimeBegin)
    return;

  mpi_printf("CIRCUMSTELLAR: Begin sink swallowing.\n");


  if(All.ComovingIntegrationOn)
    {
      ascale = All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    hubble_a = ascale = 1;

  /* Now do the swallowing of particles */
  NumGas_swallowed = 0;


  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct circumstellardata_in) + sizeof(struct circumstellardata_out) +
                                             sizemax(sizeof(struct circumstellardata_in), sizeof(struct circumstellardata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));



  idx = 0;                      /* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */

      for(nexport = 0; i <= TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(P[i].Type == 5 || P[i].Type == 4)
            if(circumstellar_evaluate_swallow(i, 0, &nexport, Send_count) < 0)
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

      CircumstellarDataGet = (struct circumstellardata_in *) mymalloc("CircumstellarDataGet", nimport * sizeof(struct circumstellardata_in));
      CircumstellarDataIn = (struct circumstellardata_in *) mymalloc("CircumstellarDataIn", nexport * sizeof(struct circumstellardata_in));

      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          for(k = 0; k < 3; k++)
            CircumstellarDataIn[j].Pos[k] = P[place].Pos[k];

          CircumstellarDataIn[j].Radius = All.CircumstellarSinkRadius;

          memcpy(CircumstellarDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }


      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              MPI_Sendrecv(&CircumstellarDataIn[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct circumstellardata_in), MPI_BYTE,
                           recvTask, TAG_DENS_A,
                           &CircumstellarDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct circumstellardata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      myfree(CircumstellarDataIn);
      CircumstellarDataResult = (struct circumstellardata_out *) mymalloc("CircumstellarDataResult", nimport * sizeof(struct circumstellardata_out));
      CircumstellarDataOut = (struct circumstellardata_out *) mymalloc("CircumstellarDataOut", nexport * sizeof(struct circumstellardata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        circumstellar_evaluate_swallow(j, 1, &dummy, &dummy);

      if(idx == TimeBinsGravity.NActiveParticles)
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
              MPI_Sendrecv(&CircumstellarDataResult[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct circumstellardata_out),
                           MPI_BYTE, recvTask, TAG_DENS_B,
                           &CircumstellarDataOut[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct circumstellardata_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      /* add the result to the particles */
      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;

          for(k = 0; k < 3; k++)
            P[place].Vel[k] = (P[place].Vel[k] * P[place].Mass + CircumstellarDataOut[j].AccretedMomentum[k]) / (P[place].Mass + CircumstellarDataOut[j].AccretedMass);

          P[place].Mass += CircumstellarDataOut[j].AccretedMass;


        }

      myfree(CircumstellarDataOut);
      myfree(CircumstellarDataResult);
      myfree(CircumstellarDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  MPI_Reduce(&NumGas_swallowed, &Ntot_gas_swallowed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  mpi_printf("CIRCUMSTELLAR: Accretion done: %d gas particles swallowed.\n", Ntot_gas_swallowed);
}




static int circumstellar_evaluate_swallow(int target, int mode, int *nexport, int *nSend_local)
{
  int startnode, numngb, j, k, n, listindex = 0;
  MyFloat accreted_mass, dmass, fac;
  MyFloat accreted_momentum[3];
  MyDouble *pos, *vel, mass;
  MyFloat h_i;

  if(mode == 0)
    {
      /* local sink */
      pos = P[target].Pos;
      vel = P[target].Vel;
      mass = P[target].Mass;
      h_i = All.CircumstellarSinkRadius;
    }
  else
    {
      /* imported sink */
      pos = CircumstellarDataGet[target].Pos;
      vel = CircumstellarDataGet[target].Vel;
      mass = CircumstellarDataGet[target].Mass;
      h_i = CircumstellarDataGet[target].Radius;
    }


  accreted_mass = 0;
  accreted_momentum[0] = accreted_momentum[1] = accreted_momentum[2] = 0;

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = CircumstellarDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb = ngb_treefind_variable(pos, h_i, target, &startnode, mode, nexport, nSend_local);

          if(numngb < 0)
            return -1;

          for(n = 0; n < numngb; n++)
            {
              j = Ngblist[n];
              /* gas cell -> drain or swallow */
              if(P[j].Type == 0)
                {
#ifndef CIRCUMSTELLAR_SINKS_ALTERNATIVE
                  double sep = sqrt((pos[0] - P[j].Pos[0]) * (pos[0] - P[j].Pos[0]) + (pos[1] - P[j].Pos[1]) * (pos[1] - P[j].Pos[1]) + (pos[2] - P[j].Pos[2]) * (pos[2] - P[j].Pos[2]));
                  double E_bind = 0.5 * ((vel[0] - P[j].Vel[0]) * (vel[0] - P[j].Vel[0]) +
                                         (vel[1] - P[j].Vel[1]) * (vel[1] - P[j].Vel[1]) + (vel[2] - P[j].Vel[2]) * (vel[2] - P[j].Vel[2])) + -All.G * (mass + P[j].Mass) / sep / sep / sep;

                  if(E_bind > 0)
                    continue;

                  double cellrad = get_cell_radius(j);
                  if(cellrad > 0.5 * All.CircumstellarSinkRadius)
                    continue;
                  if(SphP[j].MaxFaceAngle > 1.5 * All.CellMaxAngleFactor)
                    continue;
                  if(P[j].Mass <= All.TargetGasMass)
                    {
                      dmass = P[j].Mass;
                      accreted_mass += dmass;
                      for(k = 0; k < 3; k++)
                        accreted_momentum[k] += dmass * P[j].Vel[k];

                      P[j].Mass = 0;
                      P[j].ID = 0;
                    }
                  else
                    {
                      dmass = 0.9 * All.TargetGasMass;
                      fac = (P[j].Mass - dmass) / P[j].Mass;

                      accreted_mass += dmass;
                      for(k = 0; k < 3; k++)
                        accreted_momentum[k] += dmass * P[j].Vel[k];

                      P[j].Mass *= fac;
                      SphP[j].Energy *= fac;
                      SphP[j].Momentum[0] *= fac;
                      SphP[j].Momentum[1] *= fac;
                      SphP[j].Momentum[2] *= fac;
                    }
#else
                  accreted_mass += P[j].Mass;
                  accreted_momentum[0] += SphP[j].Momentum[0];
                  accreted_momentum[1] += SphP[j].Momentum[1];
                  accreted_momentum[2] += SphP[j].Momentum[2];

                  P[j].Mass = 0;
                  P[j].ID = 0;

#ifdef VORONOI_DYNAMIC_UPDATE
                  voronoi_remove_connection(j);
#endif

#endif
                  /* count drained as swallowed for statistics */

                  NumGas_swallowed++;
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = CircumstellarDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    {
      for(k = 0; k < 3; k++)
        P[target].Vel[k] = (P[target].Vel[k] * P[target].Mass + accreted_momentum[k]) / (P[target].Mass + accreted_mass);
      P[target].Mass += accreted_mass;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[target].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[target].SofteningType = get_softening_type_from_mass(P[target].Mass);
#endif

    }
  else
    {
      CircumstellarDataResult[target].AccretedMass = accreted_mass;

      for(k = 0; k < 3; k++)
        CircumstellarDataResult[target].AccretedMomentum[k] = accreted_momentum[k];
    }

  return 0;
}


#endif
