/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_recompute.c
 * \date        12/2014
 * \author      Kevin Schaal
 * \brief       Recompute hydro data on the nodes
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

//3d done

#include "../allvars.h"
#include "../proto.h"


#ifdef DG


void dg_recompute_nodes(void)
{
  TIMER_START(CPU_DG_RECOMP);
  if(All.TotNumGas == 0)
    return;

  dg_recompute_nodes_normal();
  TIMER_STOP(CPU_DG_RECOMP);
}

void dg_recompute_nodes_normal(void)
{
  int i, no;

  for(i = 0; i < NTopleaves; i++)
    {
      if(DomainTask[i] == ThisTask)
        {
          no = Ngb_DomainNodeIndex[i];

          dg_recompute_node_recursive(no, 0);
        }
    }

  dg_exchange_topleafdata();

  dg_recompute_node_recursive(Ngb_MaxPart, 1);
}




void dg_recompute_node_recursive(int no, int mode)
{
  int p;


  if(no >= Ngb_MaxPart && no < Ngb_MaxPart + Ngb_MaxNodes)      /* internal node */
    {
      p = Ngb_Nodes[no].u.d.nextnode;

#ifdef AMR
      amr_reset_node_data(no);
#endif

      while(p != Ngb_Nodes[no].u.d.sibling)
        {
          if(p < Ngb_MaxPart)   /* a particle */
            {




#ifdef AMR
              amr_accumulate_node_data(no, p);
#endif

              p = Ngb_Nextnode[p];
            }
          else if(p < Ngb_MaxPart + Ngb_MaxNodes)       /* an internal node  */
            {
              if(mode == 1)
                {
                  if(p < Ngb_FirstNonTopLevelNode)      /* we still have a top-level node */
                    dg_recompute_node_recursive(p, mode);
                }
              else
                {

                  dg_recompute_node_recursive(p, mode);

#ifdef AMR
                  amr_accumulate_node_data(no, p);
#endif
                }


              p = Ngb_Nodes[p].u.d.sibling;
            }
          else                  /* a pseudo particle */
            {
              if(mode == 0)
                terminate("should not happen");

              p = Ngb_Nextnode[p - Ngb_MaxNodes];
            }
        }


    }
  else
    {
      terminate("should not happen");
    }
}

void dg_exchange_topleafdata(void)
{
  int n, no, idx, task;
  int *recvcounts, *recvoffset, *bytecounts, *byteoffset;

  amr_node_data *DomainMoment, *loc_DomainMoment;

  DomainMoment = (amr_node_data *) mymalloc("DomainMoment", NTopleaves * sizeof(amr_node_data));

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
    bytecounts[task] = recvcounts[task] * sizeof(amr_node_data);

  for(task = 1, recvoffset[0] = 0, byteoffset[0] = 0; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  loc_DomainMoment = (amr_node_data *) mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(amr_node_data));


  for(n = 0, idx = 0; n < NTopleaves; n++)
    {
      if(DomainTask[n] == ThisTask)
        {
          no = Ngb_DomainNodeIndex[n];

          loc_DomainMoment[idx] = Ngb_Nodes[no].hydro;

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

          Ngb_Nodes[no].hydro = DomainMoment[idx];
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}


#endif
