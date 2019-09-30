/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/ism/ism_local_radiation_pressure.c
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


#ifdef ISM_LOCAL_RADIATION_PRESSURE


/* communication structures */
static struct clumplocation_data_in
{
  MyDouble Pos[3];              /* need to pass: gas cells position; gas smoothing length; current clump position; current max density */
  MyDouble ClumpPos[3];
  MyFloat h;
  MyFloat ClumpDensity;
  MyIDType ClumpID;

  int NodeList[NODELISTLENGTH];
}
 *ClumpLocationDataIn, *ClumpLocationDataGet;

static struct clumplocation_data_out
{
  MyDouble ClumpPos[3];         /* need to get clump position, clump density, flag for repeat needed if new center found */
  MyFloat ClumpDensity;
  MyFloat ClumpRadius;
  MyIDType ClumpID;

  int found_new_center;
}
 *ClumpLocationDataResult, *ClumpLocationDataOut;



void ism_find_clump_centers(void)       /* loop over all (star forming) gas particles, set field SphP[i].ClumpPos to clump center location                   */
{
  int idx, i, j, ndone, ndone_flag, dummy;
  int ngrp, sendTask, recvTask, place, nexport, nimport;
  double xtmp, ytmp, ztmp;
  int npleft = 0, iter = 0;
  long long int ntot;
  float distance = 0;

  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct clumplocation_data_in) + sizeof(struct clumplocation_data_out) +
                                             sizemax(sizeof(struct clumplocation_data_in), sizeof(struct clumplocation_data_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].ClumpPos[0] = P[i].Pos[0];
      SphP[i].ClumpPos[1] = P[i].Pos[1];
      SphP[i].ClumpPos[2] = P[i].Pos[2];
      SphP[i].ClumpDensity = 0.0;
      SphP[i].ClumpRadius = 0.0;
      SphP[i].ClumpMass = 0.0;
      SphP[i].ClumpLuminosity = 0.0;

#ifdef USE_SFR
      if(SphP[i].Sfr > 0)
#else
      if(SphP[i].Density * All.cf_a3inv >= All.DensThreshold)
#endif
        SphP[i].ClumpCenterFound = 0;
      else
        SphP[i].ClumpCenterFound = 1;
    }



  do
    {
      idx = 0;                  /* first particle for this task */
      ndone = 0;

      do
        {
          for(j = 0; j < NTask; j++)
            {
              Send_count[j] = 0;
              Exportflag[j] = -1;
            }

          /* do local particles and prepare export list */

          for(nexport = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
                continue;

              if(SphP[i].ClumpCenterFound == 0)
                {
                  SphP[i].ClumpCenterFound = 1; /* assume SphP[i].ClumpPos is already properly set because this cell is the center or we found the center during the last iteration */
                  xtmp = NGB_PERIODIC_LONG_X(P[i].Pos[0] - SphP[i].ClumpPos[0]);
                  ytmp = NGB_PERIODIC_LONG_X(P[i].Pos[1] - SphP[i].ClumpPos[1]);
                  ztmp = NGB_PERIODIC_LONG_X(P[i].Pos[2] - SphP[i].ClumpPos[2]);

                  distance = sqrt(xtmp * xtmp + ytmp * ytmp + ztmp * ztmp);

                  if(distance < 100.0 * get_cell_radius(i))     //FIXME: make this check more elegant
                    if(find_clump_centers_evaluate(i, 0, &nexport, Send_count) < 0)
                      break;
                }
            }

#ifdef MYSORT
          mysort_dataindex(DataIndexTable, nexport, sizeof(data_index), data_index_compare);
#else
          qsort(DataIndexTable, nexport, sizeof(data_index), data_index_compare);
#endif
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

          ClumpLocationDataGet = (struct clumplocation_data_in *) mymalloc("ClumpLocationDataGet", nimport * sizeof(struct clumplocation_data_in));
          ClumpLocationDataIn = (struct clumplocation_data_in *) mymalloc("ClumpLocationDataIn", nexport * sizeof(struct clumplocation_data_in));

          // prepare particle data for export 
          for(j = 0; j < nexport; j++)
            {
              place = DataIndexTable[j].Index;

              ClumpLocationDataIn[j].Pos[0] = P[place].Pos[0];
              ClumpLocationDataIn[j].Pos[1] = P[place].Pos[1];
              ClumpLocationDataIn[j].Pos[2] = P[place].Pos[2];

              ClumpLocationDataIn[j].ClumpPos[0] = SphP[place].ClumpPos[0];
              ClumpLocationDataIn[j].ClumpPos[1] = SphP[place].ClumpPos[1];
              ClumpLocationDataIn[j].ClumpPos[2] = SphP[place].ClumpPos[2];

              ClumpLocationDataIn[j].ClumpDensity = SphP[place].ClumpDensity;
              ClumpLocationDataIn[j].ClumpID = SphP[place].ClumpID;

              ClumpLocationDataIn[j].h = get_cell_radius(place);

              memcpy(ClumpLocationDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
            }

          // exchange particle data 
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              sendTask = ThisTask;
              recvTask = ThisTask ^ ngrp;

              if(recvTask < NTask)
                {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                      // get the particles 
                      MPI_Sendrecv(&ClumpLocationDataIn[Send_offset[recvTask]],
                                   Send_count[recvTask] * sizeof(struct clumplocation_data_in), MPI_BYTE,
                                   recvTask, TAG_DENS_A,
                                   &ClumpLocationDataGet[Recv_offset[recvTask]],
                                   Recv_count[recvTask] * sizeof(struct clumplocation_data_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

          myfree(ClumpLocationDataIn);

          ClumpLocationDataResult = (struct clumplocation_data_out *) mymalloc("ClumpLocationDataResult", nimport * sizeof(struct clumplocation_data_out));
          ClumpLocationDataOut = (struct clumplocation_data_out *) mymalloc("ClumpLocationDataOut", nexport * sizeof(struct clumplocation_data_out));

          // now do the particles that were sent to us
          for(j = 0; j < nimport; j++)
            find_clump_centers_evaluate(j, 1, &dummy, &dummy);

          // get the result 
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              sendTask = ThisTask;
              recvTask = ThisTask ^ ngrp;
              if(recvTask < NTask)
                {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                      // send the results 
                      MPI_Sendrecv(&ClumpLocationDataResult[Recv_offset[recvTask]],
                                   Recv_count[recvTask] * sizeof(struct clumplocation_data_out),
                                   MPI_BYTE, recvTask, TAG_DENS_B,
                                   &ClumpLocationDataOut[Send_offset[recvTask]],
                                   Send_count[recvTask] * sizeof(struct clumplocation_data_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }

          // update the local particles 

          for(j = 0; j < nexport; j++)
            {
              place = DataIndexTable[j].Index;
              if(ClumpLocationDataOut[j].ClumpDensity > SphP[place].ClumpDensity)
                {
                  SphP[place].ClumpPos[0] = ClumpLocationDataOut[j].ClumpPos[0];
                  SphP[place].ClumpPos[1] = ClumpLocationDataOut[j].ClumpPos[1];
                  SphP[place].ClumpPos[2] = ClumpLocationDataOut[j].ClumpPos[2];
                  SphP[place].ClumpDensity = ClumpLocationDataOut[j].ClumpDensity;
                  SphP[place].ClumpID = ClumpLocationDataOut[j].ClumpID;
                  SphP[place].ClumpCenterFound = 0;
                }
            }

          myfree(ClumpLocationDataOut);
          myfree(ClumpLocationDataResult);
          myfree(ClumpLocationDataGet);

          // only done if all local star particles got assigned to a hsml value 
          if(i < 0)
            ndone_flag = 1;
          else
            ndone_flag = 0;

          MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      while(ndone < NTask);

      npleft = 0;
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(SphP[i].ClumpCenterFound == 0)
            npleft++;
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("ISM: clump finding iteration %d: need to repeat for %lld particles.\n", iter, ntot);

          if(iter > 10)
            {
              mpi_printf("ISM: 10 iterations reached.\n");
              npleft = 0;
              ntot = 0;
            }
        }
    }
  while(ntot > 0);



  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

}



int find_clump_centers_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb_inbox, listindex = 0;
  double h, h2, hinv, hinv3, r2;

  MyIDType clumpID;

#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  char found_new_center = 0;

  MyFloat clumpDensity, newclumpDensity;
  MyDouble pos[3], clumpPos[3], newclumpPos[3], rclump;

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
      clumpPos[0] = SphP[target].ClumpPos[0];
      clumpPos[1] = SphP[target].ClumpPos[1];
      clumpPos[2] = SphP[target].ClumpPos[2];
      clumpDensity = SphP[target].ClumpDensity;
      clumpID = SphP[target].ClumpID;

      h = get_cell_radius(target) * 5.0;
    }
  else
    {
      pos[0] = ClumpLocationDataGet[target].Pos[0];
      pos[1] = ClumpLocationDataGet[target].Pos[1];
      pos[2] = ClumpLocationDataGet[target].Pos[2];
      clumpPos[0] = ClumpLocationDataGet[target].ClumpPos[0];
      clumpPos[1] = ClumpLocationDataGet[target].ClumpPos[1];
      clumpPos[2] = ClumpLocationDataGet[target].ClumpPos[2];
      clumpDensity = ClumpLocationDataGet[target].ClumpDensity;
      clumpID = ClumpLocationDataGet[target].ClumpID;
      h = ClumpLocationDataGet[target].h * 5.0;
    }

  newclumpDensity = clumpDensity;
  newclumpPos[0] = clumpPos[0];
  newclumpPos[1] = clumpPos[1];
  newclumpPos[2] = clumpPos[2];

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef  TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif

  if(mode == 0)
    {
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      startnode = ClumpLocationDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb_inbox = ngb_treefind_variable(clumpPos, h, target, &startnode, mode, nexport, nsend_local);

          if(numngb_inbox < 0)  // what is this doing here? will it cause an error?
            return -1;

          for(n = 0; n < numngb_inbox; n++)
            {
              j = Ngblist[n];

              xtmp = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
              ytmp = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
              ztmp = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

              r2 = xtmp * xtmp + ytmp * ytmp + ztmp * ztmp;

              if(r2 < h2 && P[j].Mass > 0 && P[j].ID != 0)
                {
                  if(SphP[j].Density > clumpDensity)
                    {
                      found_new_center = 1;
                      clumpDensity = SphP[j].Density;
                      newclumpPos[0] = P[j].Pos[0];
                      newclumpPos[1] = P[j].Pos[1];
                      newclumpPos[2] = P[j].Pos[2];
                      clumpID = P[j].ID;

                      rclump = sqrt(r2);
                      if(rclump > SphP[j].ClumpRadius)
                        SphP[j].ClumpRadius = rclump;   /*   SphP[j].ClumpRadius should now contain the 
                                                           distance from the clump center to the most distant particle, regardless 
                                                           of which Task the most distant particle resides on.  However, the 
                                                           clump constituent particles do not yet know 
                                                           the size of the clump (i.e. rclump) */

                    }
                }
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = ClumpLocationDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  if(mode == 0)
    {
      if(found_new_center)
        {
          SphP[target].ClumpPos[0] = newclumpPos[0];
          SphP[target].ClumpPos[1] = newclumpPos[1];
          SphP[target].ClumpPos[2] = newclumpPos[2];
          SphP[target].ClumpDensity = clumpDensity;
          SphP[target].ClumpID = clumpID;
          SphP[target].ClumpCenterFound = 0;    // set this to zero because we need to repeat search for new clump location
        }
    }
  else
    {
      ClumpLocationDataResult[target].ClumpPos[0] = newclumpPos[0];
      ClumpLocationDataResult[target].ClumpPos[1] = newclumpPos[1];
      ClumpLocationDataResult[target].ClumpPos[2] = newclumpPos[2];
      ClumpLocationDataResult[target].ClumpDensity = clumpDensity;
      ClumpLocationDataResult[target].ClumpID = clumpID;
      ClumpLocationDataResult[target].found_new_center = found_new_center;
    }

  return 0;
}














void ism_find_clump_radii(void) /* Extract field SphP[j].ClumpRadius for all clump particles                   */
{
  int idx, i, j, ndone, ndone_flag, dummy;
  int ngrp, sendTask, recvTask, place, nexport, nimport;

  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct clumplocation_data_in) + sizeof(struct clumplocation_data_out) +
                                             sizemax(sizeof(struct clumplocation_data_in), sizeof(struct clumplocation_data_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  idx = 0;                      /* first particle for this task */
  ndone = 0;

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      for(nexport = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(find_clump_radii_evaluate(i, 0, &nexport, Send_count) < 0)
            break;
        }

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(data_index), data_index_compare);
#endif
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

      ClumpLocationDataGet = (struct clumplocation_data_in *) mymalloc("ClumpLocationDataGet", nimport * sizeof(struct clumplocation_data_in));
      ClumpLocationDataIn = (struct clumplocation_data_in *) mymalloc("ClumpLocationDataIn", nexport * sizeof(struct clumplocation_data_in));

      // prepare particle data for export 
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          ClumpLocationDataIn[j].ClumpPos[0] = SphP[place].ClumpPos[0];
          ClumpLocationDataIn[j].ClumpPos[1] = SphP[place].ClumpPos[1];
          ClumpLocationDataIn[j].ClumpPos[2] = SphP[place].ClumpPos[2];

          ClumpLocationDataIn[j].ClumpID = SphP[place].ClumpID;

          ClumpLocationDataIn[j].h = get_cell_radius(place);

          memcpy(ClumpLocationDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }

      // exchange particle data 
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          sendTask = ThisTask;
          recvTask = ThisTask ^ ngrp;

          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  // get the particles 
                  MPI_Sendrecv(&ClumpLocationDataIn[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct clumplocation_data_in), MPI_BYTE,
                               recvTask, TAG_DENS_A,
                               &ClumpLocationDataGet[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct clumplocation_data_in), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      myfree(ClumpLocationDataIn);

      ClumpLocationDataResult = (struct clumplocation_data_out *) mymalloc("ClumpLocationDataResult", nimport * sizeof(struct clumplocation_data_out));
      ClumpLocationDataOut = (struct clumplocation_data_out *) mymalloc("ClumpLocationDataOut", nexport * sizeof(struct clumplocation_data_out));

      // now do the particles that were sent to us
      for(j = 0; j < nimport; j++)
        find_clump_radii_evaluate(j, 1, &dummy, &dummy);

      // get the result 
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          sendTask = ThisTask;
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            {
              if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                {
                  // send the results 
                  MPI_Sendrecv(&ClumpLocationDataResult[Recv_offset[recvTask]],
                               Recv_count[recvTask] * sizeof(struct clumplocation_data_out),
                               MPI_BYTE, recvTask, TAG_DENS_B,
                               &ClumpLocationDataOut[Send_offset[recvTask]],
                               Send_count[recvTask] * sizeof(struct clumplocation_data_out), MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      // update the local particles 
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          if(ClumpLocationDataOut[j].ClumpRadius > 0)
            SphP[place].ClumpRadius = ClumpLocationDataOut[j].ClumpRadius;
        }

      myfree(ClumpLocationDataOut);
      myfree(ClumpLocationDataResult);
      myfree(ClumpLocationDataGet);

      // only done if all local star particles got assigned to a hsml value 
      if(idx == TimeBinsHydro.NActiveParticles)
        ndone_flag = 1;
      else
        ndone_flag = 0;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

}



int find_clump_radii_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb_inbox, listindex = 0;
  MyFloat h;
  MyDouble clumpPos[3], rclump = 0;
  MyIDType clumpID;

  if(mode == 0)
    {
      clumpPos[0] = SphP[target].ClumpPos[0];
      clumpPos[1] = SphP[target].ClumpPos[1];
      clumpPos[2] = SphP[target].ClumpPos[2];
      clumpID = SphP[target].ClumpID;
      h = get_cell_radius(target);
      startnode = Ngb_MaxPart;  /* root node */
    }
  else
    {
      clumpPos[0] = ClumpLocationDataGet[target].ClumpPos[0];
      clumpPos[1] = ClumpLocationDataGet[target].ClumpPos[1];
      clumpPos[2] = ClumpLocationDataGet[target].ClumpPos[2];
      clumpID = ClumpLocationDataGet[target].ClumpID;
      h = ClumpLocationDataGet[target].h;
      startnode = ClumpLocationDataGet[target].NodeList[0];
      startnode = Ngb_Nodes[startnode].u.d.nextnode;    /* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
        {
          numngb_inbox = ngb_treefind_variable(clumpPos, h, target, &startnode, mode, nexport, nsend_local);

          if(numngb_inbox < 0)  // what is this doing here? will it cause an error?
            return -1;

          for(n = 0; n < numngb_inbox; n++)
            {
              j = Ngblist[n];

              if(P[j].ID == clumpID)
                rclump = SphP[j].ClumpRadius;
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              startnode = ClumpLocationDataGet[target].NodeList[listindex];
              if(startnode >= 0)
                startnode = Ngb_Nodes[startnode].u.d.nextnode;  /* open it */
            }
        }
    }

  if(mode == 0)
    {
      if(rclump > 0)
        SphP[target].ClumpRadius = rclump;
    }
  else
    ClumpLocationDataResult[target].ClumpRadius = rclump;

  return 0;
}




void ism_kick_particles(void)   /* take some action */
{
  int idx, i, k;
  MyFloat dr, dt, dtime;        //, dtime_in_Gyr;
  MyFloat EKin;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].ClumpDensity > 0 && SphP[i].ClumpLuminosity > 0)
        {
          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;

          if(All.ComovingIntegrationOn)
            dtime = All.Time * dt / All.cf_time_hubble_a;
          else
            dtime = dt;


          dtime /= All.HubbleParam;

          dr = sqrt((P[i].Pos[0] - SphP[i].ClumpPos[0]) * (P[i].Pos[0] - SphP[i].ClumpPos[0]) +
                    (P[i].Pos[1] - SphP[i].ClumpPos[1]) * (P[i].Pos[1] - SphP[i].ClumpPos[1]) + (P[i].Pos[2] - SphP[i].ClumpPos[2]) * (P[i].Pos[2] - SphP[i].ClumpPos[2]));

          /* Compute internal energy */
          EKin = 0.5 * (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] + SphP[i].Momentum[2] * SphP[i].Momentum[2]) / P[i].Mass;
          SphP[i].Energy -= EKin;

//        if(dr > 0 && SphP[i].ClumpMass > 0)
          //          for(k=0; k<3; k++)
          //          SphP[i].Momentum[k] += (P[i].Pos[k] - SphP[i].ClumpPos[k]) / dr * 
          //                               (SphP[i].ClumpLuminosity * SOLAR_LUM / (All.UnitEnergy_in_cgs / All.UnitTime_in_s)) / (CLIGHT / All.UnitVelocity_in_cm_per_s) *
          //                           (P[i].Mass / SphP[i].ClumpMass) * dtime;

          if(dr > 0 && SphP[i].ClumpMass > 0)
            for(k = 0; k < 3; k++)
              SphP[i].Momentum[k] +=
                P[i].Mass * (1.0 / SphP[i].ClumpMass +
                             All.Kappa_IR * (0.1 +
                                             SphP[i].Metallicity / 0.02) / (4 * 3.14159 * dr * dr)) * (SphP[i].ClumpLuminosity * SOLAR_LUM /
                                                                                                       (All.UnitEnergy_in_cgs / All.UnitTime_in_s)) /
                (CLIGHT / All.UnitVelocity_in_cm_per_s) * (P[i].Pos[k] - SphP[i].ClumpPos[k]) / dr * dtime;


          /* Compute new kinetic energy and add to the internal energy */
          EKin = 0.5 * (SphP[i].Momentum[0] * SphP[i].Momentum[0] + SphP[i].Momentum[1] * SphP[i].Momentum[1] + SphP[i].Momentum[2] * SphP[i].Momentum[2]) / P[i].Mass;
          SphP[i].Energy += EKin;

#ifdef ISM_VERBOSE
          printf("stats for particle i=%d:\n", i);
          printf("  dx/dr                   = %f\n", (P[i].Pos[0] - SphP[i].ClumpPos[0]) / dr);
          printf("  dy/dr                   = %f\n", (P[i].Pos[1] - SphP[i].ClumpPos[1]) / dr);
          printf("  dz/dr                   = %f\n", (P[i].Pos[2] - SphP[i].ClumpPos[2]) / dr);
          printf("  fair share luminosity   = %f\n\n", (10.0 * SphP[i].ClumpLuminosity / SphP[i].ClumpMass * P[i].Mass));

          printf("  dr                      = %f\n", dr);
          printf("  SphP[i].ClumpLuminosity = %f\n", SphP[i].ClumpLuminosity);
          printf("  SphP[i].ClumpMass       = %f\n", SphP[i].ClumpMass);
          printf("  P[i].Mass               = %f\n", P[i].Mass);
          printf("  est number of particles = %f\n", SphP[i].ClumpMass / P[i].Mass);
          printf("  dtime_in_Gyr            = %f\n\n", dtime_in_Gyr);

          printf("  P[i].Pos[0]             = %f\n", P[i].Pos[0]);
          printf("  P[i].Pos[1]             = %f\n", P[i].Pos[1]);
          printf("  P[i].Pos[2]             = %f\n\n", P[i].Pos[2]);

          printf("  SphP[i].ClumpPos[0]             = %f\n", SphP[i].ClumpPos[0]);
          printf("  SphP[i].ClumpPos[1]             = %f\n", SphP[i].ClumpPos[1]);
          printf("  SphP[i].ClumpPos[2]             = %f\n\n", SphP[i].ClumpPos[2]);

          printf("  dvx                     = %f\n", 1e9 * (P[i].Pos[0] - SphP[i].ClumpPos[0]) / dr * (10.0 * SphP[i].ClumpLuminosity / SphP[i].ClumpMass * P[i].Mass) * dtime_in_Gyr);
          printf("  dvy                     = %f\n", 1e9 * (P[i].Pos[1] - SphP[i].ClumpPos[1]) / dr * (10.0 * SphP[i].ClumpLuminosity / SphP[i].ClumpMass * P[i].Mass) * dtime_in_Gyr);
          printf("  dvz                     = %f\n\n", 1e9 * (P[i].Pos[2] - SphP[i].ClumpPos[2]) / dr * (10.0 * SphP[i].ClumpLuminosity / SphP[i].ClumpMass * P[i].Mass) * dtime_in_Gyr);

#endif

        }
    }
}


#endif
