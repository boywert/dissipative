/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/ism/ism_neighbors.c
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

#ifdef ISM

static int ism_tree_ngb_evaluate(int i, int mode, int *nexport, int *nsend_local);

static struct ism_ngbdata_in
{
  MyDouble SearchPos[3];        /* Clump center         */
  MyFloat SearchRadius;         /* Clump size           */
  int NodeList[NODELISTLENGTH];
}
 *ISMngbDataIn, *ISMngbDataGet;

static struct ism_ngbdata_out
{
  MyFloat ClumpMass;            /* Enclosed gas mass            */
  MyFloat ClumpLuminosity;      /* Enclosed stellar luminosity  */
}
 *ISMngbDataResult, *ISMngbDataOut;

extern int *TargetList;
extern int Nforces;

/* ================================================ */
/* ================================================ */
/* ================================================ */
/* ========  STILL UNDER CONSTRUCTION!!! ========== */
/* ================================================ */
/* ================================================ */
/* ================================================ */

/* This routines uses the gravitational tree to search in the R_clump neighborhood around each
 * clump center for gas and star particles to determine the enclosed gas mass and stellar luminosity.
 */
void ism_find_clump_interior_properties(void)
{
  int idx, i, j, nexport, nimport;      // ncount;
  int ngrp, recvTask, dummy, ndone, ndone_flag;

  mpi_printf("ISM: Begin clump property evaluation\n");

  /* allocate buffers to arrange communication */
  All.BunchSize = (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) + sizeof(struct ism_ngbdata_in) +
                                                           sizeof(struct ism_ngbdata_out) + sizemax(sizeof(struct ism_ngbdata_in), sizeof(struct ism_ngbdata_out))));

  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  /* Create list of targets. We do this here to simplify the treatment of the two possible locations of source points */
  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

/* ============================================================================================= */
/* ============================================================================================= */
/* ============================================================================================= */
  /* Populate the target list with all active, star forming gas, on this processor */

#ifdef USE_SFR
  Nforces = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].Sfr > 0)
        TargetList[Nforces++] = i;
    }
#else

  printf("All.DensThreshold = %f  All.cf_a3inv = %f  SphP[i].Density = %f\n", All.DensThreshold, All.cf_a3inv, SphP[i].Density);
  Nforces = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].Density * All.cf_a3inv >= All.DensThreshold)
        TargetList[Nforces++] = i;
    }
#endif


/* ============================================================================================= */
/* ============================================================================================= */
/* ============================================================================================= */

  i = 0;                        /* begin with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }

      /* do local particles and prepare export list */
      for(nexport = 0; i < Nforces; i++)
        if(ism_tree_ngb_evaluate(i, 0, &nexport, Send_count) < 0)
          break;

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

      ISMngbDataGet = (struct ism_ngbdata_in *) mymalloc("ISMngbDataGet", nimport * sizeof(struct ism_ngbdata_in));
      ISMngbDataIn = (struct ism_ngbdata_in *) mymalloc("ISMngbDataIn", nexport * sizeof(struct ism_ngbdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;
          int target = TargetList[place];

          ISMngbDataIn[j].SearchPos[0] = SphP[target].ClumpPos[0];
          ISMngbDataIn[j].SearchPos[1] = SphP[target].ClumpPos[1];
          ISMngbDataIn[j].SearchPos[2] = SphP[target].ClumpPos[2];
          ISMngbDataIn[j].SearchRadius = sqrt((P[target].Pos[0] - SphP[target].ClumpPos[0]) * (P[target].Pos[0] - SphP[target].ClumpPos[0]) +
                                              (P[target].Pos[1] - SphP[target].ClumpPos[1]) * (P[target].Pos[1] - SphP[target].ClumpPos[1]) +
                                              (P[target].Pos[2] - SphP[target].ClumpPos[2]) * (P[target].Pos[2] - SphP[target].ClumpPos[2]));
          memcpy(ISMngbDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));

//        if(target < NumPart)
//          {
//            ISMngbDataIn[j].ClumpRadius = SphP[target].ClumpRadius;
//            for(k = 0; k < 3; k++)
//              ISMngbDataIn[j].ClumpPos[k] = SphP[target].ClumpPos[k];
//          }
//        else
//          {
//            target -= Tree_ImportedNodeOffset;
//            ISMngbDataIn[j].Hsml = TSphP(target).ClumpRadius;
//            for(k = 0; k < 3; k++)
//              ISMngbDataIn[j].ClumpPos[k] = TSphP(target).ClumpPos[k];
//          }

        }

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              MPI_Sendrecv(&ISMngbDataIn[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct ism_ngbdata_in), MPI_BYTE,
                           recvTask, TAG_GRAV_A,
                           &ISMngbDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct ism_ngbdata_in), MPI_BYTE, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      myfree(ISMngbDataIn);

      ISMngbDataResult = (struct ism_ngbdata_out *) mymalloc("ISMngbDataIn", nimport * sizeof(struct ism_ngbdata_out));
      ISMngbDataOut = (struct ism_ngbdata_out *) mymalloc("ISMngbDataOut", nexport * sizeof(struct ism_ngbdata_out));

      memset(ISMngbDataResult, 0, nimport * sizeof(struct ism_ngbdata_out));
      memset(ISMngbDataOut, 0, nexport * sizeof(struct ism_ngbdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        ism_tree_ngb_evaluate(j, 1, &dummy, &dummy);

      if(i >= Nforces)
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
              MPI_Sendrecv(&ISMngbDataResult[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct ism_ngbdata_out),
                           MPI_BYTE, recvTask, TAG_GRAV_B,
                           &ISMngbDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct ism_ngbdata_out), MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      for(j = 0; j < nexport; j++)
        {
          int place = DataIndexTable[j].Index;
          int target = TargetList[place];

//        if(target < NumPart)
//          {
          SphP[target].ClumpMass += ISMngbDataOut[j].ClumpMass;
          SphP[target].ClumpLuminosity += ISMngbDataOut[j].ClumpLuminosity;
//          }
//        else
//          {
//            int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
//
//            SphP[idx].ClumpMass += ISMngbDataOut[j].ClumpMass;
//              SphP[idx].ClumpLuminosity += ISMngbDataOut[j].ClumpLuminosity;
//          }
        }
      myfree(ISMngbDataOut);
      myfree(ISMngbDataResult);
      myfree(ISMngbDataGet);
    }
  while(ndone < NTask);

  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  mpi_printf("ISM: done with ISM neighbours\n");
}


static int ism_tree_ngb_evaluate(int i, int mode, int *nexport, int *nsend_local)
{
  int target, task;
  int no, listindex = 0;
  double h, h2;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif
  double dx, dy, dz, r2;
//  MyIDType id;
  MyDouble *pos;
//  MyFloat *vel;
  MyFloat ageInGyr, lm_ssp;
  int nexport_save = *nexport;
  double enclosed_mass = 0, enclosed_luminosity = 0;

  if(mode == 0)
    {
      target = TargetList[i];
      pos = P[target].Pos;
      h = SphP[target].ClumpRadius;

//      if(target < NumPart)                    /* This is a real particle */
//      {
//        pos = P[target].Pos;
//        h = SphP[target].ClumpRadius;
//      }
//      else                                    /* This is a ghost particle -- SHOULDN'T HAPPEN!!! */
//      {
//        int n = target - Tree_ImportedNodeOffset;
//
//        pos = Tree_Points[n].Pos;
//        vel = TBPP(n).Vel;
//        h = TBPP(n).Hsml;
//        id = TBPP(n).ID;
//      }


    }
  else                          /* This is a particle on another processor */
    {
      target = i;
      pos = ISMngbDataGet[target].SearchPos;
      h = ISMngbDataGet[target].SearchRadius;
    }

  h2 = h * h;

  if(mode == 0)
    {
      no = Tree_MaxPart;        /* gravity tree root node */
    }
  else
    {
      no = ISMngbDataGet[target].NodeList[0];
      no = Nodes[no].u.d.nextnode;      /* open it */
    }

  while(no >= 0)
    {
      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle -- local processor -- easy. */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
#ifdef GFM
                  if(P[no].Type == 4)
                    {
                      ageInGyr = get_time_difference_in_Gyr(STP(no).BirthTime, All.Time);
                      lm_ssp = evaluate_l_over_m_ssp(ageInGyr);
                      enclosed_luminosity += lm_ssp * P[no].Mass * (All.UnitMass_in_g / SOLAR_MASS);
                    }

                  if(P[no].Type == 0)
                    if(SphP[no].Sfr > 0.0)
                      enclosed_mass += P[no].Mass;
#else
                  if(P[no].Type == 4)
                    {
                      ageInGyr = get_time_difference_in_Gyr(0.0, All.Time);
                      lm_ssp = evaluate_l_over_m_ssp(ageInGyr);
                      enclosed_luminosity += lm_ssp * P[no].Mass * (All.UnitMass_in_g / SOLAR_MASS);
                    }

                  if(P[no].Type == 0)
                    if(SphP[i].Density * All.cf_a3inv >= All.DensThreshold)
                      enclosed_mass += P[no].Mass;
#endif
                }

              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node -- decide if it should be opened. */
            {
              if(mode == 1)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              struct NODE *current = &Nodes[no];

              no = current->u.d.sibling;        /* in case the node can be discarded */

              double dist = h + 0.5 * current->len;
              dx = NGB_PERIODIC_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = NGB_PERIODIC_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = NGB_PERIODIC_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;
              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode;       /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist -- this is tricky.  We've found a star particle on this task, but its info is on another task */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;

              if(r2 < h2)
                {
                  if(Tree_Points[n].Type == 4)
                    {
                      ageInGyr = get_time_difference_in_Gyr(Tree_Points[n].BirthTime, All.Time);
                      lm_ssp = evaluate_l_over_m_ssp(ageInGyr);
                      enclosed_luminosity += lm_ssp * Tree_Points[n].Mass * (All.UnitMass_in_g / SOLAR_MASS);
                    }

                  if(Tree_Points[n].Type == 0)
                    if(Tree_Points[n].Sfr > 0)
                      enclosed_mass += Tree_Points[n].Mass;

                }

              no = Nextnode[no - Tree_MaxNodes];
            }
          else                  /* pseudo particle */
            {
              if(mode == 0)
                {
                  if(Exportflag[task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)]] != i)
                    {
                      Exportflag[task] = i;
                      Exportnodecount[task] = NODELISTLENGTH;
                    }

                  if(Exportnodecount[task] == NODELISTLENGTH)
                    {
                      if(*nexport >= All.BunchSize)
                        {
                          *nexport = nexport_save;
                          for(task = 0; task < NTask; task++)
                            nsend_local[task] = 0;
                          for(no = 0; no < nexport_save; no++)
                            nsend_local[DataIndexTable[no].Task]++;
                          return -1;
                        }
                      Exportnodecount[task] = 0;
                      Exportindex[task] = *nexport;
                      DataIndexTable[*nexport].Task = task;
                      DataIndexTable[*nexport].Index = i;
                      DataIndexTable[*nexport].IndexGet = *nexport;
                      *nexport = *nexport + 1;
                      nsend_local[task]++;
                    }

                  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];

                  if(Exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }
              else
                terminate("stop");

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
        }

      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              no = ISMngbDataGet[target].NodeList[listindex];
              if(no >= 0)
                {
                  no = Nodes[no].u.d.nextnode;  /* open it */
                }
            }
        }
    }

  /* store result at the proper place */
  if(mode == 0)
    {
      if(target < NumPart)      /* real particle -- local processor */
        {
          SphP[target].ClumpMass += enclosed_mass;
          SphP[target].ClumpLuminosity += enclosed_luminosity;
        }
      else
        {
//          int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
//        Tree_ISMngbResultsImported[idx].ClumpMass += enclosed_mass;
//        Tree_ISMngbResultsImported[idx].ClumpLuminosity += enclosed_luminosity;
        }
    }
  else                          /* imported data */
    {
      ISMngbDataResult[target].ClumpMass += enclosed_mass;
      ISMngbDataResult[target].ClumpLuminosity += enclosed_luminosity;
    }

  return 0;
}


/* return the (solar-scaled) light-to-mass ratio of an SSP with a given age; used throughout */
double evaluate_l_over_m_ssp(double stellar_age_in_gyr)
{
  if(stellar_age_in_gyr < 0.0029)
    {
      return 1136.59;
    }
  else
    {
      double log_age = log10(stellar_age_in_gyr) - (-2.2681);
      return 478.63 * pow(10., -1.3625 * log_age + 0.115765 * log_age * log_age);
      /* could replace with piecewise linear functions; if this call in forcetree gets expensive */
    }
}


#endif // ISM
