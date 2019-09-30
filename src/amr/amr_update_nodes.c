/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/amr/amr_mesh.c
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
#include <time.h>

#include "../allvars.h"
#include "../proto.h"

#ifdef AMR


void amr_update_nodes_recursive(int no, int mode)
{
  int j, p;

  if(mode == 1)
    {
      if(no >= Ngb_FirstNonTopLevelNode)
        {
          return;
        }
      else if(Ngb_DomainTask[no] >= 0)
        {
          return;
        }
    }


  if(no >= Ngb_MaxPart && no < Ngb_MaxPart + Ngb_MaxNodes)      /* internal node */
    {
      amr_reset_node_data(no);

      for(j = 0; j < (1 << NUMDIMS); j++)
        {
          p = Ngb_Nodes[no].u.suns[j];
          if(p >= 0)
            {
#ifndef AMR_STATIC_MESH
              if(p >= Ngb_MaxPart)
                {
                  if(p < Ngb_MaxPart + Ngb_MaxNodes)    // not a pseudo particle
                    {
                      int d;
                      if(p >= Ngb_FirstNonTopLevelNode)
                        {
                          amr_link_ngb_node(p, no, j);
                          for(d = 0; d < 2 * NUMDIMS; d++)
                            {
                              int ngb = Ngb_Nodes[p].neighbors[d];

                              if(ngb >= Ngb_MaxPart && ngb < Ngb_FirstNonTopLevelNode)  //a top level node, must be on an other task

                                {
                                  amr_export_node(p, d, ngb);
                                }
                            }
                        }
                      else
                        {
                          if(j < (1 << NUMDIMS))
                            {
                              amr_link_ngb_node(p, no, j);
                            }
                        }
                    }
                }
              else
                {
                  amr_link_ngb_particle(p, no, j);
                  int d;
                  for(d = 0; d < 2 * NUMDIMS; d++)
                    {

                      int ngb = Mesh.DP[p].neighbors[d];

                      if(ngb >= Ngb_MaxPart && ngb < Ngb_FirstNonTopLevelNode && Ngb_DomainTask[ngb] != ThisTask)       //a top level node on an other task
                        {
                          if(Mesh.NexpList == Mesh.MaxNexpList)
                            {
                              Mesh.MaxNexpList *= 1.2 + 1;      //FIXME better estiamte
                              Mesh.expList = myrealloc_movable(Mesh.expList, Mesh.MaxNexpList * sizeof(list_export_data));
                            }
                          // write this particle to an export list, so we can send it as a ghost to the other side.
                          Mesh.expList[Mesh.NexpList].task = Ngb_DomainTask[ngb];
                          Mesh.expList[Mesh.NexpList].index = p;
                          Mesh.NexpList++;
                        }
                    }
                }
#endif
              if(p >= Ngb_MaxPart && no < Ngb_MaxPart + Ngb_MaxNodes)
                {
                  amr_update_nodes_recursive(p, mode);
                }


              if(p >= 0 && p < Ngb_MaxPart + Ngb_MaxNodes)
                {
                  amr_accumulate_node_data(no, p);
                }
            }
          else if(j < (1 << NUMDIMS))
            {
              if(no > Ngb_FirstNonTopLevelNode)
                {
                  terminate("invalid amr mesh: not all subnodes are occupied in node %d\n", no);
                }
            }
        }
    }

}


void amr_exchange_topleafdata(void)
{
  int n, no, idx, task;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;
  struct DomainNODE
  {
    amr_node_data hydro;
  }
   *DomainMoment, *loc_DomainMoment;

  DomainMoment = (struct DomainNODE *) mymalloc("DomainMoment", NTopleaves * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  recvcounts = (int *) mymalloc("recvcounts", sizeof(int) * NTask);
  recvoffset = (int *) mymalloc("recvoffset", sizeof(int) * NTask);
  bytecounts = (int *) mymalloc("bytecounts", sizeof(int) * NTask);
  byteoffset = (int *) mymalloc("byteoffset", sizeof(int) * NTask);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    recvcounts[DomainTask[n]]++;

  for(task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_DomainMoment = (struct DomainNODE *) mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));


  for(n = 0, idx = 0; n < NTopleaves; n++)
    {
      if(DomainTask[n] == ThisTask)
        {
          no = Ngb_DomainNodeIndex[n];

          loc_DomainMoment[idx].hydro = Ngb_Nodes[no].hydro;

          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(n = 0; n < NTopleaves; n++)
    {
      task = DomainTask[n];
      if(task != ThisTask)
        {
          no = Ngb_DomainNodeIndex[n];
          idx = recvoffset[task] + recvcounts[task]++;

          Ngb_Nodes[no].hydro = DomainMoment[idx].hydro;

          Ngb_Nodes[no].Ti_Current = All.Ti_Current;
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}


void amr_update_nodes()
{
  int i;

  TIMER_START(CPU_AMR_UPDATE_NODES) mpi_printf("AMR: updating nodes and linkage\n");

#ifndef AMR_STATIC_MESH
  Mesh.NexpList = 0;
  Mesh.NexpListNode = 0;
  Mesh.Nvf = 0;

  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] != ThisTask)
        {
          int index = Ngb_DomainNodeIndex[i];
          int j;
          for(j = 0; j < 8; j++)
            {
              Ngb_Nodes[index].u.suns[j] = -1;
            }
        }
    }

  amr_link_toptree(Ngb_MaxPart);
#endif

  for(i = 0; i < NTopleaves; i++)
    {
      int no = Ngb_DomainNodeIndex[i];

      if(no < Ngb_MaxPart || no >= Ngb_MaxPart + Ngb_MaxNodes)
        terminate("i=%d no=%d  task=%d \n", i, no, DomainTask[i]);

      amr_update_nodes_recursive(no, 0);
    }

  amr_reset_node_data(Ngb_MaxPart);

  amr_exchange_topleafdata();

  amr_update_nodes_recursive(Ngb_MaxPart, 1);

TIMER_STOP(CPU_AMR_UPDATE_NODES)}

#endif
