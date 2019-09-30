/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/OTVET/otvet_eddington.c
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


#ifdef OTVET


int eddington_treeevaluate(int i, int mode, int *nexport, int *nsend_local);

extern int *TargetList;
extern int Nforces;

/*structures for eddington tensor*/
struct eddingtondata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  MyFloat Hsml;
}
 *EddingtonDataIn, *EddingtonDataGet;


struct eddingtondata_out
{
  MyFloat ET[6];
#ifdef RADPRESS_OPT_THICK
  MyFloat vector_n[3];
#endif
}
 *EddingtonDataResult, *EddingtonDataOut;

struct eddingtonresultsimported_data
{
  int index;
  MyFloat ET[6];
#ifdef RADPRESS_OPT_THICK
  MyFloat vector_n[3];
#endif
}
 *EddingtonResultsImported;


/* eddington tensor computation and arrangement of particle communication*/
void otvet_eddington(void)
{
  int idx, i, j, k, ngrp, dummy;
  int recvTask, nexport, nimport, place, ndone, ndone_flag;
  MyDouble trace;
  int target, ncount, n;
#ifdef RADPRESS_OPT_THICK
  MyDouble trace1;
#endif

  /* allocate buffers to arrange communication */

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct eddingtondata_in) + sizeof(struct eddingtondata_out) + sizemax(sizeof(struct eddingtondata_in), sizeof(struct eddingtondata_out))));
  DataIndexTable = (data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(data_index));
  DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));



  /* Create list of targets. We do this here to simplify the treatment of the two possible sources of points */

  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;


  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(Tree_Task_list[i] == ThisTask)
        TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(Tree_Points[i].Type == 0)
      {
        Tree_ResultIndexList[i] = ncount++;
        TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
      }

  EddingtonResultsImported = mymalloc("EddingtonResultsImported", ncount * sizeof(struct eddingtonresultsimported_data));

  i = 0;                        /* beginn with this index */

  do
    {
      for(j = 0; j < NTask; j++)
        {
          Send_count[j] = 0;
          Exportflag[j] = -1;
        }


      /* do local particles and prepare export list */
      for(nexport = 0; i < Nforces; i++)
        {
          if(eddington_treeevaluate(i, 0, &nexport, Send_count) < 0)
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

      EddingtonDataGet = (struct eddingtondata_in *) mymalloc("EddingtonDataGet", nimport * sizeof(struct eddingtondata_in));
      EddingtonDataIn = (struct eddingtondata_in *) mymalloc("EddingtonDataIn", nexport * sizeof(struct eddingtondata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          target = TargetList[place];

          if(target < NumPart)
            {
              for(k = 0; k < 3; k++)
                EddingtonDataIn[j].Pos[k] = P[target].Pos[k];
              EddingtonDataIn[j].Hsml = SphP[target].Hsml;      /* Sph-based Hsml obtained from density() */
            }
          else
            {
              target -= Tree_ImportedNodeOffset;

              for(k = 0; k < 3; k++)
                EddingtonDataIn[j].Pos[k] = Tree_Points[target].Pos[k];
              EddingtonDataIn[j].Hsml = Tree_Points[target].Hsml;
            }

          memcpy(EddingtonDataIn[j].NodeList, DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
        }



      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
        {
          recvTask = ThisTask ^ ngrp;
          if(recvTask < NTask)
            if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
              /* get the particles */
              MPI_Sendrecv(&EddingtonDataIn[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE,
                           recvTask, TAG_HYDRO_A,
                           &EddingtonDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct eddingtondata_in), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

      myfree(EddingtonDataIn);
      EddingtonDataResult = (struct eddingtondata_out *) mymalloc("EddingtonDataResult", nimport * sizeof(struct eddingtondata_out));
      EddingtonDataOut = (struct eddingtondata_out *) mymalloc("EddingtonDataOut", nexport * sizeof(struct eddingtondata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
        eddington_treeevaluate(j, 1, &dummy, &dummy);

      /* check whether this is the last iteration */
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
              MPI_Sendrecv(&EddingtonDataResult[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct eddingtondata_out),
                           MPI_BYTE, recvTask, TAG_HYDRO_B,
                           &EddingtonDataOut[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct eddingtondata_out), MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;
          target = TargetList[place];

          if(target < NumPart)
            for(k = 0; k < 6; k++)
              SphP[target].ET[k] += EddingtonDataOut[j].ET[k];

          else
            {
              int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];

              for(k = 0; k < 6; k++)
                EddingtonResultsImported[idx].ET[k] = EddingtonDataOut[j].ET[k];
            }

#ifdef RADPRESS_OPT_THICK
          if(target < NumPart)
            for(k = 0; k < 3; k++)
              SphP[target].vector_n[k] += EddingtonDataOut[j].vector_n[k];

          else
            {
              int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];

              for(k = 0; k < 3; k++)
                EddingtonResultsImported[idx].vector_n[k] = EddingtonDataOut[j].vector_n[k];
            }
#endif
        }

      myfree(EddingtonDataOut);
      myfree(EddingtonDataResult);
      myfree(EddingtonDataGet);
    }
  while(ndone < NTask);


  /* now communicate the results in EddingtonResultsImported */

  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Mesh_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
        if(Tree_Points[n].Type == 0)
          {
            EddingtonResultsImported[k].index = Tree_Points[n].index;
            Recv_count[i]++;
            k++;
          }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nexport = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  struct eddingtonresultsimported_data *tmp_results = mymalloc("tmp_results", nexport * sizeof(struct eddingtonresultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct eddingtonresultsimported_data));

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&EddingtonResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct eddingtonresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct eddingtonresultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      for(k = 0; k < 6; k++)
        SphP[target].ET[k] = tmp_results[i].ET[k];
#ifdef RADPRESS_OPT_THICK
      for(k = 0; k < 3; k++)
        SphP[target].vector_n[k] = tmp_results[i].vector_n[k];
#endif
    }

  myfree(tmp_results);

  myfree(EddingtonResultsImported);
  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  /* do final operations divide by the trace */

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      trace = SphP[i].ET[0] + SphP[i].ET[1] + SphP[i].ET[2];

      if(trace)
        for(k = 0; k < 6; k++)
          SphP[i].ET[k] /= trace;
      else
        {
          for(k = 0; k < 6; k++)
            SphP[i].ET[k] = 0.0;
          SphP[i].ET[0] = SphP[i].ET[1] = SphP[i].ET[2] = 0.333;
        }
#ifdef RADPRESS_OPT_THICK
      trace1 = SphP[i].vector_n[0] * SphP[i].vector_n[0] + SphP[i].vector_n[1] * SphP[i].vector_n[1] + SphP[i].vector_n[2] * SphP[i].vector_n[2];
      trace1 = sqrt(trace1);

      if(trace1)
        for(k = 0; k < 3; k++)
          SphP[i].vector_n[k] /= trace1;
      else
        {
          for(k = 0; k < 3; k++)
            SphP[i].vector_n[k] = 0.0;
          SphP[i].vector_n[0] = SphP[i].vector_n[1] = SphP[i].vector_n[2] = 0.333;
        }
#endif
    }

}



/*! This routine computes the eddington tensor ET for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int eddington_treeevaluate(int i, int mode, int *nexport, int *nsend_local)
{
  struct NODE *nop = 0;
  int k, no, nodesinlist, ninteractions, task, ptype;
  int target;
  int startnode, numngb, numngb_inbox, listindex = 0;
  double r4, r2, r4_invfac, dx = 0, dy = 0, dz = 0, stellar_mass = 0;
  double pos_x, pos_y, pos_z, h_i, h_j, hmax;
  MyDouble ET[6] = { 0, 0, 0, 0, 0, 0 };
#if defined(PERIODIC) && !defined(GRAVITY_NOT_PERIODIC)
  double xtmp, ytmp, ztmp;
#endif
  int nexport_save;
  nexport_save = *nexport;
#ifdef RADPRESS_OPT_THICK
  MyDouble vector_n[3] = { 0, 0, 0 };
  double r1_invfac;
  double stellar_mass_acum = 0;
#endif

  ninteractions = 0;
  nodesinlist = 0;
  ptype = 0;                    /* we only deal with gas particles */


  if(mode == 0)
    {
      target = TargetList[i];

      if(target < NumPart)
        {
          pos_x = Tree_Pos_list[3 * target + 0];
          pos_y = Tree_Pos_list[3 * target + 1];
          pos_z = Tree_Pos_list[3 * target + 2];
          h_i = SphP[target].Hsml;
        }
      else
        {
          int n = target - Tree_ImportedNodeOffset;

          pos_x = Tree_Points[n].Pos[0];
          pos_y = Tree_Points[n].Pos[1];
          pos_z = Tree_Points[n].Pos[2];
          h_i = Tree_Points[n].Hsml;
        }
    }
  else
    {
      target = i;
      pos_x = EddingtonDataGet[target].Pos[0];
      pos_y = EddingtonDataGet[target].Pos[1];
      pos_z = EddingtonDataGet[target].Pos[2];
      h_i = EddingtonDataGet[target].Hsml;
    }

  if(mode == 0)
    no = Tree_MaxPart;          /* root node */
  else
    {
      nodesinlist++;
      no = EddingtonDataGet[target].NodeList[0];
      no = Nodes[no].u.d.nextnode;      /* open it */
    }

  while(no >= 0)
    {
      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              /* the index of the node is the index of the particle */
              /* observe the sign */

#ifdef EDDINGTON_TENSOR_STARS
              if(P[no].Type == 4)
                {
                  dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);
                  dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
                  dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);
                  stellar_mass = P[no].Mass;
                }
              else
                {
                  dx = 0;
                  dy = 0;
                  dz = 0;
                  stellar_mass = 0;
                }
#endif

              r2 = dx * dx + dy * dy + dz * dz;
              if(P[no].Type == 0)
                h_j = SphP[no].Hsml;
              else
                h_j = All.ForceSoftening[P[no].SofteningType];
              if(h_j > h_i)
                hmax = h_j;
              else
                hmax = h_i;
              no = Nextnode[no];

            }
          else if(no < Tree_MaxPart + Tree_MaxNodes)    /* internal node */
            {
              if(mode == 1)
                {
                  if(no < Tree_FirstNonTopLevelNode)    /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];
              stellar_mass = nop->stellar_mass;

              dx = GRAVITY_NEAREST_X(nop->stellar_s[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(nop->stellar_s[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(nop->stellar_s[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              /* we have an  internal node. Need to check opening criterion */

              if(All.ErrTolTheta)       /* check Barnes-Hut opening criterion */
                {
                  if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                }
              else
                {
                  /* check in addition whether we lie inside the cell */
                  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
                    {
                      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
                        {
                          if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
                            {
                              no = nop->u.d.nextnode;
                              continue;
                            }
                        }
                    }
                }

              h_j = nop->u.d.maxsoft;

              if(h_j > h_i)
                {
                  if(r2 < h_j * h_j)
                    {
                      /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
                    }
                  hmax = h_j;
                }
              else
                hmax = h_i;

              no = nop->u.d.sibling;
            }

          else if(no >= Tree_ImportedNodeOffset)        /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;


/* #ifdef EDDINGTON_TENSOR_STARS */
              if(Tree_Points[n].Type == 4)
                {
                  dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
                  dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
                  dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);
                  stellar_mass = Tree_Points[n].Mass;
                }
              else
                {
                  dx = 0;
                  dy = 0;
                  dz = 0;
                  stellar_mass = 0;
                }
/* #endif */


              r2 = dx * dx + dy * dy + dz * dz;
              if(Tree_Points[n].Type == 0)
                h_j = Tree_Points[n].Hsml;
              else
                h_j = All.ForceSoftening[Tree_Points[n].SofteningType];

              if(h_j > h_i)
                hmax = h_j;
              else
                hmax = h_i;

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
                          /* out of buffer space. Need to discard work for this particle and interrupt */
                          *nexport = nexport_save;
                          if(nexport_save == 0)
                            terminate("error eddington, buffer too small\n");   /* in this case, the buffer is too small to process even a single particle */

                          nsend_local[task] = 0;
                          for(no = 0; no < nexport_save; no++)
                            nsend_local[DataIndexTable[no].Task]++;
                          return -1;
                        }
                      Exportnodecount[task] = 0;
                      Exportindex[task] = *nexport;
                      DataIndexTable[*nexport].Task = task;
                      DataIndexTable[*nexport].Index = target;
                      DataIndexTable[*nexport].IndexGet = *nexport;
                      *nexport = *nexport + 1;
                      nsend_local[task]++;
                    }

                  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];

                  if(Exportnodecount[task] < NODELISTLENGTH)
                    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }


          r4 = r2 * r2;
          r4_invfac = 1.0 / (r4 + hmax * hmax * hmax * hmax);

          ET[0] += stellar_mass * dx * dx * r4_invfac;
          ET[1] += stellar_mass * dy * dy * r4_invfac;
          ET[2] += stellar_mass * dz * dz * r4_invfac;
          ET[3] += stellar_mass * dx * dy * r4_invfac;
          ET[4] += stellar_mass * dy * dz * r4_invfac;
          ET[5] += stellar_mass * dz * dx * r4_invfac;


#ifdef RADPRESS_OPT_THICK

/*
          vector_n[0] += stellar_mass * -dx;  // -dx is radially outwards from star 
          vector_n[1] += stellar_mass * -dy;  
          vector_n[2] += stellar_mass * -dz;  
*/
          r1_invfac = 1.0 / (sqrt(r2) + hmax);
          vector_n[0] += stellar_mass * -dx * r1_invfac;        // -dx is radially outwards from star 
          vector_n[1] += stellar_mass * -dy * r1_invfac;
          vector_n[2] += stellar_mass * -dz * r1_invfac;
          stellar_mass_acum += stellar_mass;
#endif

          if(stellar_mass > 0)
            ninteractions++;
        }
      if(mode == 1)
        {
          listindex++;
          if(listindex < NODELISTLENGTH)
            {
              no = EddingtonDataGet[target].NodeList[listindex];
              if(no >= 0)
                {
                  nodesinlist++;
                  no = Nodes[no].u.d.nextnode;  /* open it */
                }
            }
        }
    }

  /* store result at the proper place */
  if(mode == 0)
    {
      if(target < NumPart)
        for(k = 0; k < 6; k++)
          SphP[target].ET[k] = ET[k];
      else
        {
          int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
          for(k = 0; k < 6; k++)
            EddingtonResultsImported[idx].ET[k] = ET[k];
        }
    }
  else
    {
      for(k = 0; k < 6; k++)
        EddingtonDataResult[target].ET[k] = ET[k];
      *nexport = nodesinlist;
    }

#ifdef RADPRESS_OPT_THICK
  if(stellar_mass_acum > 0)
    for(k = 0; k < 3; k++)
      vector_n[k] /= stellar_mass_acum;

  /* store result at the proper place */
  if(mode == 0)
    {
      if(target < NumPart)
        for(k = 0; k < 3; k++)
          SphP[target].vector_n[k] = vector_n[k];
      else
        {
          int idx = Tree_ResultIndexList[target - Tree_ImportedNodeOffset];
          for(k = 0; k < 3; k++)
            EddingtonResultsImported[idx].vector_n[k] = vector_n[k];
        }
    }
  else
    {
      for(k = 0; k < 3; k++)
        EddingtonDataResult[target].vector_n[k] = vector_n[k];
    }
#endif


  return ninteractions;
}

#endif
