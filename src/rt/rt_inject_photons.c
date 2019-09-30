/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/rt/rt_inject_photons.c
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


void rt_advect_main(void)
{
  /* get the time step in internal units */
  double dt = rt_get_advect_tistep();
  double delta = (All.Ti_Current - All.Ti_LastRadTransfer) * All.Timebase_interval;
  double steps = delta / dt;

#ifdef RT_CGMETHOD
  dt = delta;
  steps = 1;
#endif

  if(steps > 0)
    {
      int i, isteps = (int) (steps + 1);

      dt = delta / isteps;

      mpi_printf("----> Subcycling radiative transfer in %d cycles, dt=%g\n", isteps, dt);

      for(i = 0; i < isteps; i++)
        {
          mpi_printf("cycle =%d|%d\n", i + 1, isteps);
#ifndef RT_SPREAD_SOURCE
          rt_inject_photons_single(&Mesh, dt);
#else
          rt_inject_photons_spread(&Mesh, dt);
#endif
          /* update mesh variables and compute gradients */
          rt_voronoi_exchange_primitive_variables();
          rt_calculate_green_gauss_gradients();
          rt_voronoi_exchange_primitive_variables_and_gradients();

#ifndef RT_CGMETHOD
          rt_advect_radiation(&Mesh, dt);
#else
          rt_cgmethod(&Mesh, dt);
#endif

          rt_update_chemistry(dt);
        }
    }
  All.Ti_LastRadTransfer = All.Ti_Current;
}


struct diff_list
{
  int sourceid;
  double dPhotons;
  MyFloat sourcepos[3];
}
 *diff;

#if !defined(RT_SPREAD_SOURCE) && !defined(USE_SFR)
void rt_inject_photons_single(tessellation * T, double dt)
{
  int i, j, p;
  double dphot;

  mesh_search_data search[N_SOURCES];

  for(i = 0; i < N_SOURCES; i++)
    {
      search[i].Pos[0] = Source_Pos[i][0];
      search[i].Pos[1] = Source_Pos[i][1];
      search[i].Pos[2] = Source_Pos[i][2];
    }

  find_nearest_meshpoint_global(search, N_SOURCES, 0, 1);


  for(i = 0; i < N_SOURCES; i++)
    {

      dphot = Source_Lum[i] * dt * All.UnitTime_in_s;

      if(search[i].Task == ThisTask)    /* point lies in local cell p on this CPU */
        {
          p = search[i].u.Index;

#ifdef RT_HEALPIX_NSIDE
          for(j = 0; j < RT_N_DIR; j++)
            SphP[p].Photons[j] += dphot / RT_N_DIR;
#else

          for(j = 1; j < RT_N_DIR; j++)
            {
              if(SphP[p].SourceID[j] == Source_ID[i])
                {
                  printf("sourceid=%d p=%d j=%d dphot=%g  SphP[p].Photons[j]=%g\n", Source_ID[i], p, j, dphot, SphP[p].Photons[j]);

                  SphP[p].Photons[j] += dphot;
                  break;
                }
            }
          if(j >= RT_N_DIR)     /* we could not find this source ID */
            {
              int k;
              /* look for the weakest source */
              for(j = 1, k = 1; j < RT_N_DIR; j++)
                {
                  if(SphP[p].Photons[j] < SphP[p].Photons[k])
                    k = j;
                }

              if(SphP[p].Photons[k] < dphot)
                {
                  /* ok, replace this source */
                  SphP[p].Photons[0] += SphP[p].Photons[k];
                  SphP[p].Photons[k] = dphot;
                  SphP[p].SourceID[k] = Source_ID[i];

                  /* let's assign the cell centre as source position to avoid discreteness effects in the first
                     advection step out of the source cell */
                  SphP[p].SourcePos[k][0] = SphP[p].Center[0];
                  SphP[p].SourcePos[k][1] = SphP[p].Center[1];
                  SphP[p].SourcePos[k][2] = SphP[p].Center[2];


                  /* note: this may have destroyed the order according to source ID (needed for gradient estimate),
                     so let's recreate this */

                  struct diff_list *data_for_sort = mymalloc("to_sort", RT_N_DIR * sizeof(struct diff_list));

                  for(j = 1; j < RT_N_DIR; j++)
                    {
                      data_for_sort[j].sourceid = SphP[p].SourceID[j];
                      data_for_sort[j].dPhotons = SphP[p].Photons[j];
                      memcpy(data_for_sort[j].sourcepos, &SphP[p].SourcePos[j][0], 3 * sizeof(MyFloat));
                    }
                  /* let's sort according to sourceid, which is adopted in the gradient estimation */
                  mysort(data_for_sort + 1, RT_N_DIR - 1, sizeof(struct diff_list), rt_compare_difflist_sourceid);

                  for(j = 1; j < RT_N_DIR; j++)
                    {
                      SphP[p].SourceID[j] = data_for_sort[j].sourceid;
                      SphP[p].Photons[j] = data_for_sort[j].dPhotons;
                      memcpy(&SphP[p].SourcePos[j][0], data_for_sort[j].sourcepos, 3 * sizeof(MyFloat));
                    }

                  myfree(data_for_sort);
                }
              else
                {
                  /* source seems to be so weak that we add it to the diffuse component right away */
                  SphP[p].Photons[0] += dphot;
                }
            }
#endif
          for(j = 0; j < RT_N_DIR; j++)
            SphP[p].DensPhot[j] = SphP[p].Photons[j] / SphP[p].Volume;
        }
    }
}
#endif


/* spread a predefined list of sources */
#if defined(RT_SPREAD_SOURCE) && !defined(USE_SFR)
static struct sourcedata_in
{
  MyDouble Pos[3];
  double Center[3];
  double Volume;
  int NodeList[NODELISTLENGTH];
}
 *SourceDataIn, *SourceDataGet;


void rt_inject_photons_spread(tessellation * T, double dt)
{
  int i, j, target, dummy;
  int ngrp, sendTask, recvTask, place, nexport, nimport;
  double dphotons, hubble_a;

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

  for(i = 0; i < N_SOURCES; i++)
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      nexport = 0;

      /* first check whether source position is in local cell */
      target = rt_get_cell(T, Source_Pos[i][0], Source_Pos[i][1], Source_Pos[i][2]);

      dphotons = Source_Lum[i] * dt / hubble_a * All.UnitTime_in_s;

      if(target >= 0 && target < NumGas)        /* point lies in local cell p on this CPU */
        if(rt_source_evaluate(target, 0, &nexport, Send_count, Source_ID[i], dphotons) < 0)
          terminate("An error occurred int the source evaluation");

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
          SourceDataIn[j].Center[0] = SphP[place].Center[0];
          SourceDataIn[j].Center[1] = SphP[place].Center[1];
          SourceDataIn[j].Center[2] = SphP[place].Center[2];
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
                               recvTask, TAG_DENS_A,
                               &SourceDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct sourcedata_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(SourceDataIn);

      /* now do the particles that were sent to us */

      for(j = 0; j < nimport; j++)
        rt_source_evaluate(j, 1, &dummy, &dummy, Source_ID[i], dphotons);

      myfree(SourceDataGet);
    }

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

}


int rt_source_evaluate(int target, int mode, int *nexport, int *nsend_local, int source_id, double dphotons)
{
  int i, k, j, n;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3;
  double wk, dphot, cx, cy, cz;
  double dx, dy, dz, r2, *pos, u, r;

  if(mode == 0)
    {
      pos = P[target].Pos;
      h = 4.0 * pow(3.0 / 4.0 / M_PI * SphP[target].Volume, 0.33);
      cx = SphP[target].Center[0];
      cy = SphP[target].Center[1];
      cz = SphP[target].Center[2];
    }
  else
    {
      pos = SourceDataGet[target].Pos;
      h = 4.0 * pow(3.0 / 4.0 / M_PI * SourceDataGet[target].Volume, 0.33);
      cx = SourceDataGet[target].Center[0];
      cy = SourceDataGet[target].Center[1];
      cz = SourceDataGet[target].Center[2];
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

#ifdef RT_HEALPIX_NSIDE
                  for(i = 0; i < RT_N_DIR; i++)
                    SphP[j].Photons[i] += dphot / RT_N_DIR;
#else
                  for(i = 1; i < RT_N_DIR; i++)
                    {
                      if(SphP[j].SourceID[i] == source_id)
                        {
                          SphP[j].Photons[i] += dphot;
                          break;
                        }
                    }
                  if(i >= RT_N_DIR)     /* we could not find this source ID */
                    {
                      /* look for the weakest source */
                      for(i = 1, k = 1; i < RT_N_DIR; i++)
                        {
                          if(SphP[j].Photons[i] < SphP[j].Photons[k])
                            k = i;
                        }

                      if(SphP[j].Photons[k] < dphot)
                        {
                          /* ok, replace this source */
                          SphP[j].Photons[0] += SphP[j].Photons[k];
                          SphP[j].Photons[k] = dphot;
                          SphP[j].SourceID[k] = source_id;


                          /* let's assign the cell centre as source position to avoid discreteness effects in the first
                             advection step out of the source cell */
                          SphP[j].SourcePos[k][0] = cx;
                          SphP[j].SourcePos[k][1] = cy;
                          SphP[j].SourcePos[k][2] = cz;

                          /* note: this may have destroyed the order according to source ID (needed for gradient estimate),
                             so let's recreate this */

                          struct diff_list *data_for_sort = mymalloc("data_for_sort", RT_N_DIR * sizeof(struct diff_list));

                          for(i = 1; i < RT_N_DIR; i++)
                            {
                              data_for_sort[i].sourceid = SphP[j].SourceID[i];
                              data_for_sort[i].dPhotons = SphP[j].Photons[i];
                              memcpy(data_for_sort[i].sourcepos, &SphP[j].SourcePos[i][0], 3 * sizeof(MyFloat));
                            }
                          /* let's sort according to sourceid, which is adopted in the gradient estimation */
                          mysort(data_for_sort + 1, RT_N_DIR - 1, sizeof(struct diff_list), rt_compare_difflist_sourceid);

                          for(i = 1; i < RT_N_DIR; i++)
                            {
                              SphP[j].SourceID[i] = data_for_sort[i].sourceid;
                              SphP[j].Photons[i] = data_for_sort[i].dPhotons;
                              memcpy(&SphP[j].SourcePos[i][0], data_for_sort[i].sourcepos, 3 * sizeof(MyFloat));
                            }

                          myfree(data_for_sort);
                        }
                      else
                        {
                          /* source seems to be so weak that we add it to the diffuse component right away */
                          SphP[j].Photons[0] += dphot;
                        }
                    }
#endif
                  for(i = 0; i < RT_N_DIR; i++)
                    SphP[j].DensPhot[i] = SphP[j].Photons[i] / SphP[j].Volume;
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
