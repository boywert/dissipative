/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_refinement.c
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

#include <string.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef AMR
#ifdef REFINEMENT

void amr_refmap();
static void amr_refmap_loop(int mode);

static int amr_compare_partlist_task_index(const void *a, const void *b);
static int amr_compare_abort_task_index(const void *a, const void *b);

static int amr_check_refinement(int node);
static int amr_check_refinement_smooth(int node);
static void amr_smooth_node(int node, int status);
static void amr_smooth_0(int node, int reason);
static void amr_smooth_1(int node, int side);
#if NUMDIMS >=2
static void amr_smooth_2(int node, int side);
#endif
#if NUMDIMS >=3
static void amr_smooth_3(int node, int side);
#endif


#ifdef REFINEMENT_SPLIT_CELLS
static void amr_refine_cells();
static int amr_refine_cell(int i, MyIDType * newid, int *count);
static void amr_set_cell_pos(int cell, int node, int subnode);
static void amr_split_hydro(int i, double fac);
#endif

#ifdef REFINEMENT_MERGE_CELLS
static void amr_derefine_cells();
static int amr_derefine_cell(int no);
static void amr_sum_hydro(int new_cell, int i);
static void amr_remove_nextinlevel(int no);
#endif

#ifndef FORCE_EQUAL_TIMESTEPS
static int amr_smooth_abort(int node, int mode);

static void amr_smooth_abort_1(int node, int ngb, int reason);
#if NUMDIMS >=2
static void amr_smooth_abort_2(int node, int ngb, int reason);
#endif
#if NUMDIMS >=3
static void amr_smooth_abort_3(int node, int ngb, int reason);
#endif

static int amr_mark_abort(int node, int ngb, int reason);
#endif

static int status_smooth_a;
static int status_smooth_b;

static int Nexport = 0;
static int Nimport = 0;

static struct abort_refinement
{
  int Task;
  int Index;
  int Reason;
}
 *AbortRefinement;


void do_derefinements_and_refinements()
{
  TIMER_START(CPU_AMR_REFINEMENT) if(All.TotNumGas == 0)
    return;

  mpi_printf("AMR: generate new refinement map\n");

  //amr_validate_mesh(Ngb_MaxPart, -1, -1);
  //amr_validate_lists();

  Mesh.Nrefine = 0;
  Mesh.refcells = (int *) mymalloc_movable(&Mesh.refcells, "amr_refcells", Mesh.Ndp * sizeof(int));

  Mesh.Nderefine = 0;
  Mesh.derefcand = (int *) mymalloc_movable(&Mesh.derefcand, "amr_derefcand", Ngb_NumNodes * sizeof(int));

  Mesh.Refflag = (int *) mymalloc_movable(&Mesh.Refflag, "Refflag", (Ngb_MaxPart + Mesh.nodes_total + Mesh.Nghost_dp) * sizeof(int));
  memset(Mesh.Refflag, 0, (Ngb_MaxPart + Mesh.nodes_total + Mesh.Nghost_dp) * sizeof(int));
  Mesh.allocated = 3;

  amr_refmap();

  //remove ghost part/old interface list
  Mesh.Ndp = Mesh.LastDP;

#ifdef REFINEMENT_SPLIT_CELLS
  amr_refine_cells();
#endif
#ifdef REFINEMENT_MERGE_CELLS
  amr_derefine_cells();
#endif

  MPI_Allreduce(MPI_IN_PLACE, &Mesh.maxlevel, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &Mesh.minlevel, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  Mesh.LastDP = Mesh.Ndp;

  myfree(Mesh.Refflag);
  myfree(Mesh.derefcand);
  myfree(Mesh.refcells);
  Mesh.allocated = 2;

  Mesh.Nrefine = 0;
  Mesh.Nderefine = 0;
TIMER_STOP(CPU_AMR_REFINEMENT)}


struct refdata_in
{
  int index;
  int refflag;
};

void amr_refmap()
{
  PartList = (struct data_partlist *) mymalloc("PartList", Mesh.Nghost_nodes + Mesh.Nghost_dp * sizeof(struct data_partlist));
#ifndef FORCE_EQUAL_TIMESTEPS
  AbortRefinement = (struct abort_refinement *) mymalloc("AbortRefinement", Mesh.Nghost_nodes + Mesh.Nghost_dp * sizeof(struct abort_refinement));
#endif
  Nexport = 0;

  status_smooth_a = AMR_REFINE_SMOOTH_A;
  status_smooth_b = AMR_REFINE_SMOOTH_A;

  amr_refmap_loop(0);

  status_smooth_a = AMR_REFINE_SMOOTH_A;
  status_smooth_b = AMR_REFINE_SMOOTH_B;

  int i;
  for(i = All.AMRMeshSmoothing; i > 0; i--)
    {
      amr_refmap_loop(i);

      int tmp = status_smooth_a;
      status_smooth_a = status_smooth_b;
      status_smooth_b = tmp;
    }

#ifndef FORCE_EQUAL_TIMESTEPS
  myfree(AbortRefinement);
#endif
  myfree(PartList);
}

static void amr_refmap_loop(int mode)
{
  int level;
  int j;

  for(level = Mesh.maxlevel; level >= 0; level--)
    {
      int node;

      node = Mesh.lastinlevel[level];
      while(node >= 0)
        {
          int status;
          if(mode == 0)
            {
              status = amr_check_refinement(node);
            }
          else
            {
              status = amr_check_refinement_smooth(node);
            }

          if(status & AMR_REFINE_EXPAND)
            {
              amr_smooth_node(node, status);
            }

          if(node >= Ngb_MaxPart)
            node = Ngb_Nodes[node].nextinlevel;
          else
            node = Mesh.DP[node].nextinlevel;
        }

      int count;

      for(count = 0; count < NUMDIMS; count++)
        {

          qsort(PartList, Nexport, sizeof(struct data_partlist), amr_compare_partlist_task_index);

          for(j = 0; j < NTask; j++)
            {
              Send_count[j] = 0;
              Send_count_nodes[j] = 0;
            }

          for(j = 0; j < Nexport; j++)
            Send_count[PartList[j].Task]++;

          MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

          for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
            {
              Nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          struct refdata_in *RefDataGet = (struct refdata_in *) mymalloc("RefDataGet", Nimport * sizeof(struct refdata_in));
          struct refdata_in *RefDataIn = (struct refdata_in *) mymalloc("RefDataIn", Nexport * sizeof(struct refdata_in));

          /* prepare particle data for export */
          for(j = 0; j < Nexport; j++)
            {
              int index = PartList[j].Index;
              if(index < Ngb_MaxPart + Mesh.nodes_total)
                {
                  RefDataIn[j].index = Mesh.original_index[index];
                  RefDataIn[j].refflag = Mesh.Refflag[index];
                }
              else
                {
                  RefDataIn[j].refflag = Mesh.Refflag[index];
                  index -= Mesh.nodes_total;
                  RefDataIn[j].index = Mesh.DP[index].originalindex;
                }
            }

          int ngrp, recvTask;
          MPI_Status status;

          /* exchange particle data */
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              recvTask = ThisTask ^ ngrp;

              if(recvTask < NTask)
                {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                      /* get the particles */
                      MPI_Sendrecv(&RefDataIn[Send_offset[recvTask]],
                                   Send_count[recvTask] * sizeof(struct refdata_in), MPI_BYTE,
                                   recvTask, TAG_GRAV_A, &RefDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct refdata_in), MPI_BYTE, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
                    }
                }
            }

          myfree(RefDataIn);

          Nexport = 0;

          for(j = 0; j < Nimport; j++)
            {
              int refflag = RefDataGet[j].refflag;
              int index = RefDataGet[j].index;

              assert(index < Ngb_MaxPart + Ngb_MaxNodes);

              refflag = ((Mesh.Refflag[index] & AMR_REF_ANY) ^ refflag) & refflag;

              if(refflag & AMR_REF1)
                {
                  amr_smooth_0(index, refflag & AMR_REF1);
                }

#if NUMDIMS >= 2
              if(refflag & AMR_REF2)
                {
                  amr_smooth_0(index, refflag & AMR_REF2);

                  if(refflag & AMR_REF2_LEFT)
                    amr_smooth_1(index, AMR_REF1_LEFT);
                  if(refflag & AMR_REF2_RIGHT)
                    amr_smooth_1(index, AMR_REF1_RIGHT);
                }

#if NUMDIMS >=3
              if(refflag & AMR_REF3)
                {
                  amr_smooth_0(index, refflag & AMR_REF3);

                  if(refflag & AMR_REF3_LEFT)
                    {
                      amr_smooth_1(index, AMR_REF1_LEFT);
                    }
                  if(refflag & AMR_REF3_RIGHT)
                    {
                      amr_smooth_1(index, AMR_REF1_RIGHT);
                    }

                  if(refflag & AMR_REF3_FRONT)
                    {
                      int tmp = 0;
                      if(refflag & AMR_REF3_FRONT_LEFT)
                        tmp |= AMR_REF2_FRONT_LEFT;
                      if(refflag & AMR_REF3_FRONT_RIGHT)
                        tmp |= AMR_REF2_FRONT_RIGHT;

                      amr_smooth_2(index, tmp);
                    }
                  if(refflag & AMR_REF3_BACK)
                    {
                      int tmp = 0;
                      if(refflag & AMR_REF3_BACK_LEFT)
                        tmp |= AMR_REF2_BACK_LEFT;
                      if(refflag & AMR_REF3_BACK_RIGHT)
                        tmp |= AMR_REF2_BACK_RIGHT;

                      amr_smooth_2(index, tmp);
                    }
                }
#endif
#endif
            }
          myfree(RefDataGet);
        }
    }

#ifndef FORCE_EQUAL_TIMESTEPS
  for(level = 0; level <= Mesh.maxlevel; level++)
    {
      int node = Mesh.lastinlevel[level];
      while(node >= 0)
        {
          //abort surrounding cells
          if(Mesh.Refflag[node] & AMR_REFINE_ABORT)
            {
              amr_smooth_abort(node, mode);
            }

          if(node >= Ngb_MaxPart)
            node = Ngb_Nodes[node].nextinlevel;
          else
            node = Mesh.DP[node].nextinlevel;
        }

      int count;
      for(count = 0; count < NUMDIMS; count++)
        {

          qsort(AbortRefinement, Nexport, sizeof(struct abort_refinement), amr_compare_abort_task_index);

          for(j = 0; j < NTask; j++)
            {
              Send_count[j] = 0;
              Send_count_nodes[j] = 0;
            }

          for(j = 0; j < Nexport; j++)
            Send_count[AbortRefinement[j].Task]++;

          MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

          for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
            {
              Nimport += Recv_count[j];

              if(j > 0)
                {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
                }
            }

          struct refdata_in *RefDataGet = (struct refdata_in *) mymalloc("RefDataGet", Nimport * sizeof(struct refdata_in));
          struct refdata_in *RefDataIn = (struct refdata_in *) mymalloc("RefDataIn", Nexport * sizeof(struct refdata_in));

          /* prepare particle data for export */
          for(j = 0; j < Nexport; j++)
            {
              int index = AbortRefinement[j].Index;
              if(index < Ngb_MaxPart + Mesh.nodes_total)
                {
                  RefDataIn[j].index = Mesh.original_index[index];
                  RefDataIn[j].refflag = AbortRefinement[j].Reason;
                }
              else
                {
                  RefDataIn[j].refflag = AbortRefinement[j].Reason;
                  index -= Mesh.nodes_total;
                  RefDataIn[j].index = Mesh.DP[index].originalindex;
                }
            }

          int ngrp, recvTask;
          MPI_Status status;

          /* exchange particle data */
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
            {
              recvTask = ThisTask ^ ngrp;

              if(recvTask < NTask)
                {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                    {
                      /* get the particles */
                      MPI_Sendrecv(&RefDataIn[Send_offset[recvTask]],
                                   Send_count[recvTask] * sizeof(struct refdata_in), MPI_BYTE,
                                   recvTask, TAG_GRAV_A, &RefDataGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct refdata_in), MPI_BYTE, recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
                    }
                }
            }

          myfree(RefDataIn);

          Nexport = 0;

          for(j = 0; j < Nimport; j++)
            {
              int refflag = RefDataGet[j].refflag;
              int index = RefDataGet[j].index;


              assert(index < Ngb_MaxPart + Ngb_MaxNodes);


              if(refflag & AMR_REF1)
                {
                  amr_smooth_abort_1(-1, index, refflag);
                }

#if NUMDIMS >= 2
              if(refflag & AMR_REF2)
                {
                  amr_smooth_abort_2(-1, index, refflag);
                }

#if NUMDIMS >=3
              if(refflag & AMR_REF3)
                {
                  amr_smooth_abort_3(-1, index, refflag);
                }
#endif
#endif

            }
          myfree(RefDataGet);

        }

    }
#endif

}


static int amr_compare_partlist_task_index(const void *a, const void *b)
{
  if(((struct data_partlist *) a)->Task < (((struct data_partlist *) b)->Task))
    return -1;

  if(((struct data_partlist *) a)->Task > (((struct data_partlist *) b)->Task))
    return +1;

  if(((struct data_partlist *) a)->Index < (((struct data_partlist *) b)->Index))
    return -1;

  if(((struct data_partlist *) a)->Index > (((struct data_partlist *) b)->Index))
    return +1;

  return 0;
}

static int amr_compare_abort_task_index(const void *a, const void *b)
{
  if(((struct abort_refinement *) a)->Task < (((struct abort_refinement *) b)->Task))
    return -1;

  if(((struct abort_refinement *) a)->Task > (((struct abort_refinement *) b)->Task))
    return +1;

  if(((struct abort_refinement *) a)->Index < (((struct abort_refinement *) b)->Index))
    return -1;

  if(((struct abort_refinement *) a)->Index > (((struct abort_refinement *) b)->Index))
    return +1;

  return 0;
}

int amr_check_refinement(int node)
{
  if(node < Ngb_MaxPart)        // a cell
    {
      if(TimeBinSynchronized[P[node].TimeBinHydro] && Mesh.DP[node].level < All.MaxRefLevel)
        {
          int status = amr_should_this_node_be_split(Mesh.DP[node].index);

          if(status > 0)
            {
              Mesh.Refflag[node] |= AMR_REFINE;
              Mesh.Refflag[node] |= status_smooth_a;
              if((Mesh.Refflag[node] & AMR_REFM) == 0)
                {
                  Mesh.refcells[Mesh.Nrefine++] = Mesh.DP[node].index;
                  Mesh.Refflag[node] |= AMR_REFM;
                }
              return AMR_REFINE;
            }
        }
    }
  else                          //a local node
    {
      int status = 0;

      if(node < Ngb_FirstNonTopLevelNode ||  Ngb_Nodes[node].level < All.MinRefLevel)
        {
          status |= AMR_REFINE;
        }
      if(amr_should_this_node_be_split(node) > 0)
        {
          status |= AMR_REFINE;
        }

      if(node >= Ngb_FirstNonTopLevelNode)
        {
          int i;
          for(i = 0; i < (1 << NUMDIMS); i++)
            {
              int sun = Ngb_Nodes[node].u.suns[i];
              int s = 0;

              if(sun < Ngb_MaxPart)
                {
#ifdef OUTPUT_AMR_REFFLAG
                  if(amr_should_this_node_be_split(node) > 0)
                    {
                      SphP[sun].refflag |= AMR_REFINE_CRITERION;
                    }
#endif

                  s |= Mesh.Refflag[sun];

                  //avoid derefinement of inactive cells
                  if(!TimeBinSynchronized[P[Mesh.DP[sun].index].TimeBinHydro])
                    {
                      status |= AMR_REFINE;
                    }
                }
              else
                {
                  s |= Mesh.Refflag[sun];
                  status |= AMR_REFINE; //enforce slow derefinement
                }

              if(s > 0)
                {
                  if((i & 1) == 1)
                    status |= AMR_REFINE_EXPAND_RIGHT;
                  else
                    status |= AMR_REFINE_EXPAND_LEFT;
#if NUMDIMS >= 2
                  if((i & 2) == 2)
                    status |= AMR_REFINE_EXPAND_BACK;
                  else
                    status |= AMR_REFINE_EXPAND_FRONT;
#if NUMDIMS >=3
                  if((i & 4) == 4)
                    status |= AMR_REFINE_EXPAND_TOP;
                  else
                    status |= AMR_REFINE_EXPAND_BOTTOM;
#endif
#endif
                }
            }
        }

      if(status > 0)
        {
          Mesh.Refflag[node] |= AMR_REFINE;
          Mesh.Refflag[node] |= status_smooth_a;

          return status;
        }

      if(Mesh.Refflag[node] == 0)
        {
          Mesh.derefcand[Mesh.Nderefine++] = node;
          return 0;
        }
    }
  return 0;
}

static int amr_check_refinement_smooth(int node)
{
  if(node < Ngb_MaxPart)        // a cell
    {
      if((Mesh.Refflag[node] & status_smooth_a) && !(Mesh.Refflag[node] & AMR_SMOOTH))
        {
          Mesh.Refflag[node] |= AMR_SMOOTH;

          int status = 0;

          status |= AMR_REFINE_EXPAND_RIGHT;
          status |= AMR_REFINE_EXPAND_LEFT;
#if NUMDIMS >= 2
          status |= AMR_REFINE_EXPAND_BACK;
          status |= AMR_REFINE_EXPAND_FRONT;
#if NUMDIMS >=3
          status |= AMR_REFINE_EXPAND_TOP;
          status |= AMR_REFINE_EXPAND_BOTTOM;
#endif
#endif
          return status;
        }
    }
  else if(node < Ngb_MaxPart + Ngb_MaxNodes)    //a local node
    {
      int status = 0;
      if((Mesh.Refflag[node] & status_smooth_a) && !(Mesh.Refflag[node] & AMR_SMOOTH))
        {
          Mesh.Refflag[node] |= AMR_SMOOTH;

          status |= AMR_REFINE_EXPAND_RIGHT;
          status |= AMR_REFINE_EXPAND_LEFT;
#if NUMDIMS >= 2
          status |= AMR_REFINE_EXPAND_BACK;
          status |= AMR_REFINE_EXPAND_FRONT;
#if NUMDIMS >=3
          status |= AMR_REFINE_EXPAND_TOP;
          status |= AMR_REFINE_EXPAND_BOTTOM;
#endif
#endif
        }

      if(node >= Ngb_FirstNonTopLevelNode)
        {
          int i;
          for(i = 0; i < 1 << NUMDIMS; i++)
            {
              int sun = Ngb_Nodes[node].u.suns[i];

              int s = 0;

              if(sun < Ngb_MaxPart)
                {
                  s |= Mesh.Refflag[sun];
                  if(!TimeBinSynchronized[P[Mesh.DP[sun].index].TimeBinHydro])
                    {
                      status |= AMR_REFINE;
                    }
                }
              else
                {
                  s |= Mesh.Refflag[sun];
                  status |= AMR_REFINE; //enforce slow derefinement
                }

              if(s > 0)
                {
                  if((i & 1) == 1)
                    status |= AMR_REFINE_EXPAND_RIGHT;
                  else
                    status |= AMR_REFINE_EXPAND_LEFT;
#if NUMDIMS >= 2
                  if((i & 2) == 2)
                    status |= AMR_REFINE_EXPAND_BACK;
                  else
                    status |= AMR_REFINE_EXPAND_FRONT;
#if NUMDIMS >=3
                  if((i & 4) == 4)
                    status |= AMR_REFINE_EXPAND_TOP;
                  else
                    status |= AMR_REFINE_EXPAND_BOTTOM;
#endif
#endif
                }
            }
        }

      if(status > 0)
        {
          Mesh.Refflag[node] |= AMR_REFINE;
          return status;
        }
    }
  return 0;
}

void amr_smooth_node(int node, int status)
{
  if(status & AMR_REFINE_EXPAND_LEFT)
    {
      amr_smooth_1(node, AMR_REF1_LEFT);

#if NUMDIMS >=2
      if(status & AMR_REFINE_EXPAND_FRONT)
        {
          amr_smooth_2(node, AMR_REF2_FRONT_LEFT);

#if NUMDIMS >=3
          if(status & AMR_REFINE_EXPAND_BOTTOM)
            amr_smooth_3(node, AMR_REF3_BOTTOM_FRONT_LEFT);
          if(status & AMR_REFINE_EXPAND_TOP)
            amr_smooth_3(node, AMR_REF3_TOP_FRONT_LEFT);
#endif
        }
      if(status & AMR_REFINE_EXPAND_BACK)
        {
          amr_smooth_2(node, AMR_REF2_BACK_LEFT);

#if NUMDIMS >=3
          if(status & AMR_REFINE_EXPAND_BOTTOM)
            amr_smooth_3(node, AMR_REF3_BOTTOM_BACK_LEFT);
          if(status & AMR_REFINE_EXPAND_TOP)
            amr_smooth_3(node, AMR_REF3_TOP_BACK_LEFT);
#endif
        }
#endif
    }
  if(status & AMR_REFINE_EXPAND_RIGHT)
    {
      amr_smooth_1(node, AMR_REF1_RIGHT);

#if NUMDIMS >=2
      if(status & AMR_REFINE_EXPAND_FRONT)
        {
          amr_smooth_2(node, AMR_REF2_FRONT_RIGHT);

#if NUMDIMS >=3
          if(status & AMR_REFINE_EXPAND_BOTTOM)
            amr_smooth_3(node, AMR_REF3_BOTTOM_FRONT_RIGHT);
          if(status & AMR_REFINE_EXPAND_TOP)
            amr_smooth_3(node, AMR_REF3_TOP_FRONT_RIGHT);
#endif
        }
      if(status & AMR_REFINE_EXPAND_BACK)
        {
          amr_smooth_2(node, AMR_REF2_BACK_RIGHT);

#if NUMDIMS >=3
          if(status & AMR_REFINE_EXPAND_BOTTOM)
            amr_smooth_3(node, AMR_REF3_BOTTOM_BACK_RIGHT);
          if(status & AMR_REFINE_EXPAND_TOP)
            amr_smooth_3(node, AMR_REF3_TOP_BACK_RIGHT);
#endif
        }
#endif
    }
}

void amr_smooth_0(int node, int reason)
{

  if(node < Ngb_MaxPart)
    {
      if((Mesh.Refflag[node] & reason) == 0)
        {
          Mesh.Refflag[node] |= reason;

          Mesh.Refflag[node] |= status_smooth_b;


          if((Mesh.Refflag[node] & AMR_REFM) == 0)
            {
              if(TimeBinSynchronized[P[node].TimeBinHydro])
                {
                  Mesh.Refflag[node] |= AMR_REFINE;
                  Mesh.refcells[Mesh.Nrefine++] = Mesh.DP[node].index;
                  Mesh.Refflag[node] |= AMR_REFM;
                }
#ifndef FORCE_EQUAL_TIMESTEPS
              else
                {
                  Mesh.Refflag[node] |= AMR_REFINE_ABORT;
                  Mesh.Refflag[node] |= reason;
                }
#endif
            }
        }
    }
  else if(node < Ngb_MaxPart + Ngb_MaxNodes)
    {
      if((Mesh.Refflag[node] & reason) == 0)
        {
          Mesh.Refflag[node] |= AMR_REFINE;
          Mesh.Refflag[node] |= reason;
        }

      Mesh.Refflag[node] |= status_smooth_b;

    }
  else
    {
      terminate("error\n");
    }
}

void amr_smooth_1(int node, int reason)
{
  int ngb;
  int side;
  if(reason & AMR_REF1_LEFT)
    side = 0;
  else
    side = 1;

  if(node < Ngb_MaxPart)
    {
      ngb = Mesh.DP[node].neighbors[side];
    }
  else
    {
      ngb = Ngb_Nodes[node].neighbors[side];
    }

  if(ngb < Ngb_MaxPart)
    {
      if((Mesh.Refflag[ngb] & reason) == 0)
        {
          if((Mesh.Refflag[ngb] & AMR_REFM) == 0)
            {
              if(TimeBinSynchronized[P[ngb].TimeBinHydro])
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE;
                  Mesh.refcells[Mesh.Nrefine++] = Mesh.DP[ngb].index;
                  Mesh.Refflag[ngb] |= AMR_REFM;
                }
#ifndef FORCE_EQUAL_TIMESTEPS
              else
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE_ABORT;
                  Mesh.Refflag[ngb] |= reason;
                  return;
                }
#endif
            }

          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;


        }
    }
  else if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
    {
      if((Mesh.Refflag[ngb] & reason) == 0)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;
        }
    }
  //ghost nodes/cells
  else if(ngb < Ngb_MaxPart + Mesh.nodes_total)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          PartList[Nexport].Task = Mesh.original_task[ngb];
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.original_task[ngb] < NTask);
          assert(Mesh.original_task[ngb] > -1);
          assert(Mesh.original_task[ngb] != ThisTask);
        }
    }
  else
    {
      int index = ngb - Mesh.nodes_total;
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          PartList[Nexport].Task = Mesh.DP[index].task;
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.DP[index].task < NTask);
          assert(Mesh.DP[index].task > -1);
          assert(Mesh.DP[index].task != ThisTask);
        }
    }
}

#if NUMDIMS >=2
void amr_smooth_2(int node, int reason)
{
  int ngb;
  int side;

  if(reason & AMR_REF2_FRONT)
    side = 2;
  else
    side = 3;

  if(node < Ngb_MaxPart)
    {
      ngb = Mesh.DP[node].neighbors[side];
    }
  else
    {
      ngb = Ngb_Nodes[node].neighbors[side];
    }


  if(ngb < Ngb_MaxPart)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          if((Mesh.Refflag[ngb] & AMR_REFM) == 0)
            {
              if(TimeBinSynchronized[P[ngb].TimeBinHydro])
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE;
                  Mesh.refcells[Mesh.Nrefine++] = Mesh.DP[ngb].index;
                  Mesh.Refflag[ngb] |= AMR_REFM;
                }
#ifndef FORCE_EQUAL_TIMESTEPS
              else
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE_ABORT;
                  Mesh.Refflag[ngb] |= reason;
                  return;
                }
#endif
            }

          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          if(reason & AMR_REF2_LEFT)
            amr_smooth_1(ngb, AMR_REF1_LEFT);
          if(reason & AMR_REF2_RIGHT)
            amr_smooth_1(ngb, AMR_REF1_RIGHT);
        }
    }
  else if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          if(reason & AMR_REF2_LEFT)
            amr_smooth_1(ngb, AMR_REF1_LEFT);
          if(reason & AMR_REF2_RIGHT)
            amr_smooth_1(ngb, AMR_REF1_RIGHT);
        }
    }
  //ghost nodes/cells
  else if(ngb < Ngb_MaxPart + Mesh.nodes_total)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          PartList[Nexport].Task = Mesh.original_task[ngb];
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.original_task[ngb] < NTask);
          assert(Mesh.original_task[ngb] > -1);
          assert(Mesh.original_task[ngb] != ThisTask);
        }
    }
  else
    {
      int index = ngb - Mesh.nodes_total;
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          PartList[Nexport].Task = Mesh.DP[index].task;
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.DP[index].task < NTask);
          assert(Mesh.DP[index].task > -1);
          assert(Mesh.DP[index].task != ThisTask);
        }
    }
}
#endif

#if NUMDIMS >=3
void amr_smooth_3(int node, int reason)
{
  int ngb;
  int side;

  if(reason & AMR_REF3_BOTTOM)
    side = 4;
  else
    side = 5;

  if(node < Ngb_MaxPart)
    {
      ngb = Mesh.DP[node].neighbors[side];
    }
  else
    {
      ngb = Ngb_Nodes[node].neighbors[side];
    }


  if(ngb < Ngb_MaxPart)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          if((Mesh.Refflag[ngb] & AMR_REFM) == 0)
            {
              if(TimeBinSynchronized[P[ngb].TimeBinHydro])
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE;
                  Mesh.refcells[Mesh.Nrefine++] = Mesh.DP[ngb].index;
                  Mesh.Refflag[ngb] |= AMR_REFM;
                }
#ifndef FORCE_EQUAL_TIMESTEPS
              else
                {
                  Mesh.Refflag[ngb] |= AMR_REFINE_ABORT;
                  Mesh.Refflag[ngb] |= reason;
                  return;
                }
#endif
            }

          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          if(reason & AMR_REF3_LEFT)
            amr_smooth_1(ngb, AMR_REF1_LEFT);
          if(reason & AMR_REF3_RIGHT)
            amr_smooth_1(ngb, AMR_REF1_RIGHT);
          if(reason & AMR_REF3_FRONT)
            {
              int tmp = 0;
              if(reason & AMR_REF3_LEFT)
                tmp |= AMR_REF2_FRONT_LEFT;
              if(reason & AMR_REF3_RIGHT)
                tmp |= AMR_REF2_FRONT_RIGHT;

              amr_smooth_2(ngb, tmp);
            }
          if(reason & AMR_REF3_BACK)
            {
              int tmp = 0;
              if(reason & AMR_REF3_LEFT)
                tmp |= AMR_REF2_BACK_LEFT;
              if(reason & AMR_REF3_RIGHT)
                tmp |= AMR_REF2_BACK_RIGHT;

              amr_smooth_2(ngb, tmp);
            }
        }
    }
  else if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          if(reason & AMR_REF3_LEFT)
            amr_smooth_1(ngb, AMR_REF1_LEFT);
          if(reason & AMR_REF3_RIGHT)
            amr_smooth_1(ngb, AMR_REF1_RIGHT);
          if(reason & AMR_REF3_FRONT)
            {
              int tmp = 0;
              if(reason & AMR_REF3_LEFT)
                tmp |= AMR_REF2_FRONT_LEFT;
              if(reason & AMR_REF3_RIGHT)
                tmp |= AMR_REF2_FRONT_RIGHT;

              amr_smooth_2(ngb, tmp);
            }
          if(reason & AMR_REF3_BACK)
            {
              int tmp = 0;
              if(reason & AMR_REF3_LEFT)
                tmp |= AMR_REF2_BACK_LEFT;
              if(reason & AMR_REF3_RIGHT)
                tmp |= AMR_REF2_BACK_RIGHT;

              amr_smooth_2(ngb, tmp);
            }

        }
    }
  //ghost nodes/cells
  else if(ngb < Ngb_MaxPart + Mesh.nodes_total)
    {
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;


          PartList[Nexport].Task = Mesh.original_task[ngb];
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.original_task[ngb] < NTask);
          assert(Mesh.original_task[ngb] > -1);
          assert(Mesh.original_task[ngb] != ThisTask);
        }
    }
  else
    {
      int index = ngb - Mesh.nodes_total;
      if((Mesh.Refflag[ngb] & reason) != reason)
        {
          Mesh.Refflag[ngb] |= AMR_REFINE;
          Mesh.Refflag[ngb] |= reason;

          Mesh.Refflag[ngb] |= status_smooth_b;

          PartList[Nexport].Task = Mesh.DP[index].task;
          PartList[Nexport].Index = ngb;
          Nexport++;

          assert(Mesh.DP[index].task < NTask);
          assert(Mesh.DP[index].task > -1);
          assert(Mesh.DP[index].task != ThisTask);
        }
    }
}
#endif

#ifndef FORCE_EQUAL_TIMESTEPS

int amr_smooth_abort(int node, int modus)
{

  assert(Mesh.Refflag[node] & AMR_REFINE_ABORT);
  assert(node < Ngb_MaxPart + Ngb_MaxNodes);

  int* ngb;

  if(node < Ngb_MaxPart)
    {
      ngb = Mesh.DP[node].neighbors;
    }
  else
    {
      ngb = Ngb_Nodes[node].neighbors;
    }

  if(Mesh.Refflag[node] & AMR_REF1)
    {
      if(Mesh.Refflag[node] & AMR_REF1_LEFT)
        {

          amr_smooth_abort_1(node, ngb[1], AMR_REF1_LEFT);
        }

      if(Mesh.Refflag[node] & AMR_REF1_RIGHT)
        {
          amr_smooth_abort_1(node, ngb[0], AMR_REF1_RIGHT);
        }
    }

#if NUMDIMS >=2
  if(Mesh.Refflag[node] & AMR_REF2)
    {
      if(Mesh.Refflag[node] & AMR_REF2_FRONT_LEFT)
        {
          amr_smooth_abort_2(node, ngb[3], AMR_REF2_FRONT_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF2_FRONT_RIGHT)
        {
          amr_smooth_abort_2(node, ngb[3], AMR_REF2_FRONT_RIGHT);
        }

      if(Mesh.Refflag[node] & AMR_REF2_BACK_LEFT)
        {
          amr_smooth_abort_2(node, ngb[2], AMR_REF2_BACK_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF2_BACK_RIGHT)
        {
          amr_smooth_abort_2(node, ngb[2], AMR_REF2_BACK_RIGHT);
        }
    }
#endif

#if NUMDIMS >=3
  if(Mesh.Refflag[node] & AMR_REF3)
    {
      if(Mesh.Refflag[node] & AMR_REF3_BOTTOM_FRONT_LEFT)
        {
          amr_smooth_abort_3(node, ngb[5], AMR_REF3_BOTTOM_FRONT_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_BOTTOM_FRONT_RIGHT)
        {
          amr_smooth_abort_3(node, ngb[5], AMR_REF3_BOTTOM_FRONT_RIGHT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_BOTTOM_BACK_LEFT)
        {
          amr_smooth_abort_3(node, ngb[5], AMR_REF3_BOTTOM_BACK_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_BOTTOM_BACK_RIGHT)
        {
          amr_smooth_abort_3(node, ngb[5], AMR_REF3_BOTTOM_BACK_RIGHT);
        }

      if(Mesh.Refflag[node] & AMR_REF3_TOP_FRONT_LEFT)
        {
          amr_smooth_abort_3(node, ngb[4], AMR_REF3_TOP_FRONT_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_TOP_FRONT_RIGHT)
        {
          amr_smooth_abort_3(node, ngb[4], AMR_REF3_TOP_FRONT_RIGHT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_TOP_BACK_LEFT)
        {
          amr_smooth_abort_3(node, ngb[4], AMR_REF3_TOP_BACK_LEFT);
        }
      if(Mesh.Refflag[node] & AMR_REF3_TOP_BACK_RIGHT)
        {
          amr_smooth_abort_3(node, ngb[4], AMR_REF3_TOP_BACK_RIGHT);
        }
    }
#endif

  return 0;
}

void amr_smooth_abort_1(int node, int ngb, int reason)
{
  if(amr_mark_abort(node, ngb, reason))
    {
      return;
    }

#if NUMDIMS >=2
  int* ngb2;
  if(ngb < Ngb_MaxPart)
    {
      ngb2 = Mesh.DP[ngb].neighbors;
    }
  else
    {
      ngb2 = Ngb_Nodes[ngb].neighbors;
    }

  if(reason & AMR_REF1_LEFT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF2_FRONT_LEFT)
        {
          amr_smooth_abort_2(ngb, ngb2[2], AMR_REF2_FRONT_LEFT);
        }

      if(Mesh.Refflag[ngb] & AMR_REF2_BACK_LEFT)
        {
          amr_smooth_abort_2(ngb, ngb2[3], AMR_REF2_BACK_LEFT);
        }
    }

  if(reason & AMR_REF1_RIGHT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF2_FRONT_RIGHT)
        {
          amr_smooth_abort_2(ngb, ngb2[2], AMR_REF2_FRONT_RIGHT);
        }

      if(Mesh.Refflag[ngb] & AMR_REF2_BACK_RIGHT)
        {
          amr_smooth_abort_2(ngb, ngb2[3], AMR_REF2_BACK_RIGHT);
        }
    }
#endif
}

#if NUMDIMS >=2
void amr_smooth_abort_2(int node, int ngb, int reason)
{
  if(amr_mark_abort(node, ngb, reason))
    {
      return;
    }

#if NUMDIMS >=3
  int* ngb2;
  if(ngb < Ngb_MaxPart)
    {
      ngb2 = Mesh.DP[ngb].neighbors;
    }
  else
    {
      ngb2 = Ngb_Nodes[ngb].neighbors;
    }

  if(reason & AMR_REF2_FRONT_LEFT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF3_BOTTOM_FRONT_LEFT)
        {
          amr_smooth_abort_3(ngb, ngb2[4], AMR_REF3_BOTTOM_FRONT_LEFT);
        }
      if(Mesh.Refflag[ngb] & AMR_REF3_TOP_FRONT_LEFT)
        {
          amr_smooth_abort_3(ngb, ngb2[5], AMR_REF3_TOP_FRONT_LEFT);
        }
    }

  if(reason & AMR_REF2_FRONT_RIGHT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF3_BOTTOM_FRONT_RIGHT)
        {
          amr_smooth_abort_3(ngb, ngb2[4], AMR_REF3_BOTTOM_FRONT_RIGHT);
        }
      if(Mesh.Refflag[ngb] & AMR_REF3_TOP_FRONT_RIGHT)
        {
          amr_smooth_abort_3(ngb, ngb2[5], AMR_REF3_TOP_FRONT_RIGHT);
        }
    }

  if(reason & AMR_REF2_BACK_LEFT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF3_BOTTOM_BACK_LEFT)
        {
          amr_smooth_abort_3(ngb, ngb2[4], AMR_REF3_BOTTOM_BACK_LEFT);
        }
      if(Mesh.Refflag[ngb] & AMR_REF3_TOP_BACK_LEFT)
        {
          amr_smooth_abort_3(ngb, ngb2[5], AMR_REF3_TOP_BACK_LEFT);
        }
    }

  if(reason & AMR_REF2_BACK_RIGHT)
    {
      if(Mesh.Refflag[ngb] & AMR_REF3_BOTTOM_BACK_RIGHT)
        {
          amr_smooth_abort_3(ngb, ngb2[4], AMR_REF3_BOTTOM_BACK_RIGHT);
        }
      if(Mesh.Refflag[ngb] & AMR_REF3_TOP_BACK_RIGHT)
        {
          amr_smooth_abort_3(ngb, ngb2[5], AMR_REF3_TOP_BACK_RIGHT);
        }
    }



#endif
}
#endif

#if NUMDIMS >=3
void amr_smooth_abort_3(int node, int ngb, int reason)
{
  if(amr_mark_abort(node, ngb, reason))
    {
      return;
    }
}
#endif


int amr_mark_abort(int node, int ngb, int reason)
{


  if(node >=0)
    Mesh.Refflag[node] &= ~reason;

  if(ngb < Ngb_MaxPart)
    {
      return 0; //should only happen during smoothing
    }
  else if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart+Ngb_MaxNodes)
    {
      if(reason & (AMR_REF1_LEFT | AMR_REF2_BACK_LEFT | AMR_REF3_BOTTOM_BACK_LEFT))
        {
          if(Ngb_Nodes[ngb].u.suns[2] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[2]] |= AMR_REFINE_ABORT;
        }
      if(reason & (AMR_REF1_LEFT | AMR_REF2_FRONT_LEFT | AMR_REF3_BOTTOM_FRONT_LEFT))
        {
          if(Ngb_Nodes[ngb].u.suns[0] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[0]] |= AMR_REFINE_ABORT;
        }

#if NUMDIMS >=2
      if(reason & (AMR_REF1_RIGHT | AMR_REF2_BACK_RIGHT | AMR_REF3_BOTTOM_BACK_RIGHT))
        {
          if(Ngb_Nodes[ngb].u.suns[3] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[3]] |= AMR_REFINE_ABORT;
        }

      if(reason & (AMR_REF1_RIGHT | AMR_REF2_FRONT_RIGHT | AMR_REF3_BOTTOM_FRONT_RIGHT))
        {
          if(Ngb_Nodes[ngb].u.suns[1] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[1]] |= AMR_REFINE_ABORT;
        }
#endif

#if NUMDIMS >=3
      if(reason & (AMR_REF1_LEFT | AMR_REF2_BACK_LEFT | AMR_REF3_TOP_BACK_LEFT))
        {
          if(Ngb_Nodes[ngb].u.suns[6] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[6]] |= AMR_REFINE_ABORT;
        }
      if(reason & (AMR_REF1_LEFT | AMR_REF2_FRONT_LEFT | AMR_REF3_TOP_FRONT_LEFT))
        {
          if(Ngb_Nodes[ngb].u.suns[4] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[4]] |= AMR_REFINE_ABORT;
        }

      if(reason & (AMR_REF1_RIGHT | AMR_REF2_BACK_RIGHT | AMR_REF3_TOP_BACK_RIGHT))
        {
          if(Ngb_Nodes[ngb].u.suns[7] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[7]] |= AMR_REFINE_ABORT;
        }

      if(reason & (AMR_REF1_RIGHT | AMR_REF2_FRONT_RIGHT | AMR_REF3_TOP_FRONT_RIGHT))
        {
          if(Ngb_Nodes[ngb].u.suns[5] < Ngb_MaxPart)
          Mesh.Refflag[Ngb_Nodes[ngb].u.suns[5]] |= AMR_REFINE_ABORT;
        }
#endif

      return 0;
    }

  //ghost nodes/cells
  else if(ngb < Ngb_MaxPart + Mesh.nodes_total)
    {
      AbortRefinement[Nexport].Task = Mesh.original_task[ngb];
      AbortRefinement[Nexport].Index = ngb;
      AbortRefinement[Nexport].Reason = reason;
      Nexport++;

      assert(Mesh.original_task[ngb] < NTask);
      assert(Mesh.original_task[ngb] > -1);
      assert(Mesh.original_task[ngb] != ThisTask);
    }
  else
    {
      int index = ngb - Mesh.nodes_total;

      AbortRefinement[Nexport].Task = Mesh.DP[index].task;
      AbortRefinement[Nexport].Index = ngb;
      AbortRefinement[Nexport].Reason = reason;
      Nexport++;

      assert(Mesh.DP[index].task < NTask);
      assert(Mesh.DP[index].task > -1);
      assert(Mesh.DP[index].task != ThisTask);
    }

  return 1;
}
#endif

#ifdef REFINEMENT_SPLIT_CELLS
void amr_refine_cells()
{
  int i;
  int count, countall;
  int countabort, countabortall;
  MyIDType newid = 0;

  count = Mesh.Nrefine * ((1 << NUMDIMS) - 1);

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("AMR: want to refine %d cells\n", countall / ((1 << NUMDIMS) - 1));

  if(countall)
    {
      domain_resize_storage(count, count, 2);

      if(Mesh.Ndp + count >= Mesh.MaxNdp)
        {
          Mesh.DP = myrealloc_movable(Mesh.DP, (Mesh.MaxNdp + count) * sizeof(point));
          Mesh.MaxNdp += count;
        }

      if(NumPart + count >= All.MaxPart)
        {
          terminate("On Task=%d with NumPart=%d we try to produce %d cells. Sorry, no space left...(All.MaxPart=%d)\n", ThisTask, NumPart, count, All.MaxPart);
        }

      if(NumGas + count >= All.MaxPartSph)
        {
          terminate("On Task=%d with NumGas=%d we try to produce %d cells. Sorry, no space left...(All.MaxPartSph=%d)\n", ThisTask, NumGas, count, All.MaxPartSph);
        }


      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      int *list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&count, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      newid = All.MaxID + 1;

      for(i = 0; i < ThisTask; i++)
        newid += list[i];

      All.MaxID += countall;

      myfree(list);
    }

  if(Ngb_NextFreeNode + count > Ngb_MaxPart + Ngb_MaxNodes)
    {
      ngb_treerealloc(Ngb_NextFreeNode + count - (Ngb_MaxPart + Ngb_MaxNodes));
    }
  count = 0;
  countabort = 0;

  int refined = 0;

  for(i = 0; i < Mesh.Nrefine; i++)
    {
      int cell = Mesh.refcells[i];
      if((Mesh.Refflag[cell] & AMR_REFINE) && ((Mesh.Refflag[cell] & AMR_REFINE_ABORT)==0))
        {
          assert(TimeBinSynchronized[P[cell].TimeBinHydro]);
          int node = amr_refine_cell(cell, &newid, &count);
          Mesh.refcells[i] = node;
          refined++;
        }
      else
        {
          countabort++;
        }
    }

  MPI_Allreduce(&countabort, &countabortall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(countabortall > 0)
    {
      mpi_printf("AMR: aborted refinement of %d cells, refined %d cells\n", countabortall, countall - countabortall);
    }

  MPI_Allreduce(&refined, &Mesh.Refined, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  NumPart += count;
  NumGas += count;
  All.TotNumPart += countall;
  All.TotNumGas += countall;
}

/** This function does the refinement of the cell at index cell
 *
 */
int amr_refine_cell(int i, MyIDType * newid, int *count)
{
  int father;
  int no;
  int j, k;

  father = Mesh.DP[i].father;

  Ngb_NumNodes++;
  no = Ngb_NextFreeNode++;
  assert(no < (Ngb_MaxNodes + Ngb_MaxPart));
  assert(no >= Ngb_MaxPart);

  int subnode = 0;
  if(P[i].Pos[0] > Ngb_Nodes[father].Center[0])
    subnode += 1;
#if NUMDIMS >=2
  if(P[i].Pos[1] > Ngb_Nodes[father].Center[1])
    subnode += 2;
#if NUMDIMS >=3
  if(P[i].Pos[2] > Ngb_Nodes[father].Center[2])
    subnode += 4;
#endif
#endif



  //create  a new node
  Ngb_Nodes[father].u.suns[subnode] = no;

  Ngb_Nodes[no].level = Mesh.DP[i].level;

  amr_set_node_data(no, father, subnode, 1, &Mesh);

  amr_remove_nextinlevel(i);

  //insert that node into nextnode list
  int nnext = father;
  int next;
  do
    {
      next = nnext;
      if(next >= Ngb_MaxPart)
        {
          nnext = Ngb_Nodes[next].u.d.nextnode;
        }
      else
        {
          nnext = Ngb_Nextnode[next];
        }
    }
  while(nnext != i);

  if(next >= Ngb_MaxPart)
    {
      Ngb_Nodes[next].u.d.nextnode = no;
    }
  else
    {
      Ngb_Nextnode[next] = no;
    }

  next = father;
  do
    {
      nnext = Ngb_Nodes[next].u.d.nextnode;
      while(nnext != i && nnext != no && nnext >= 0)
        {
          next = nnext;
          if(next >= Ngb_MaxPart)
            {
              nnext = Ngb_Nodes[next].u.d.sibling;
            }
          else
            {
              nnext = Ngb_Nextnode[next];
            }
        }

      if(next >= Ngb_MaxPart && nnext == i)
        {
          Ngb_Nodes[next].u.d.sibling = no;
        }


    }
  while(next >= Ngb_MaxPart && nnext == i);


  Ngb_Nodes[no].u.d.nextnode = i;
  Ngb_Nodes[no].u.d.sibling = Ngb_Nextnode[i];

  //copy hydro quantities to new cells
  int num;
  int last = i;
  for(num = 1; num < (1 << NUMDIMS); num++)
    {
      int addToGravList = TimeBinSynchronized[P[i].TimeBinGrav];
      if(NumPart > NumGas)
        {
          move_collisionless_particle(NumPart + *count, NumGas + *count);
          if(TimeBinSynchronized[P[NumPart + *count].TimeBinGrav])
            addToGravList = 0;
#ifdef TRACER_PARTICLE
          if(P[NumPart + count].Type == TRACER_PARTICLE)
            addToGravList = 1;
#endif
          /* there is already an entry in the list of active particles for 
             gravity that points to the index that we will use for our new cell */
        }

      j = NumGas + *count;

      assert(j < Ngb_MaxPart);

      P[j] = P[i];
      SphP[j] = SphP[i];

      P[j].ID = *newid;
      *newid = *newid + 1;

#ifdef AMR_CONNECTIONS
      SphP[j].first_connection = SphP[j].last_connection = -1;
#endif

      Ngb_Nodes[no].u.suns[num] = j;
      amr_set_cell_pos(j, no, num);
      amr_set_cell_data(j, no, 1, &Mesh);
      Ngb_Father[j] = no;
#ifdef DG
      dg_split_hydro(num, j);
#else

      amr_split_hydro(j, 1. / (1 << NUMDIMS));
#endif
#ifdef TRACER_MC
      P[j].TracerHead = -1;     /* P[j] was copied from P[i] so reset its TLL to empty */
      P[j].NumberOfTracers = 0;
      consider_moving_tracers_local(i, j, 1. / ((1 << NUMDIMS) - (num - 1)));
#endif
      Mesh.Ndp++;

      /* add the new particle into the neighbour tree */
      int tmp = Ngb_Nextnode[last];
      Ngb_Nextnode[last] = j;
      Ngb_Nextnode[j] = tmp;

      /* now add the new particle into the link-lists for the time integration */

      timebin_add_particle(&TimeBinsHydro, j, i, P[i].TimeBinHydro, 1);
      timebin_add_particle(&TimeBinsGravity, j, i, P[i].TimeBinGrav, addToGravList);

      *count = *count + 1;
      last = j;
    }

  //recycle old cell as first new cell
  Ngb_Nodes[no].u.suns[0] = i;
  amr_set_cell_pos(i, no, 0);
  amr_set_cell_data(i, no, 1, &Mesh);
  Ngb_Father[i] = no;
#ifdef DG
  dg_split_hydro(0, i);
#else
  amr_split_hydro(i, 1. / (1 << NUMDIMS));
#endif


  for(k = 0; k < 2 * NUMDIMS; k++)
    {
      Ngb_Nodes[no].neighbors[k] = -1;  //Mesh.DP[i].neighbors[k];
    }

  return no;
}

void amr_set_cell_pos(int cell, int node, int subnode)
{
  if(subnode & 1)
    P[cell].Pos[0] = Ngb_Nodes[node].Center[0] + amr_length[Ngb_Nodes[node].level + 1] / 2.;
  else
    P[cell].Pos[0] = Ngb_Nodes[node].Center[0] - amr_length[Ngb_Nodes[node].level + 1] / 2.;
#if NUMDIMS >=2
  if(subnode & 2)
    P[cell].Pos[1] = Ngb_Nodes[node].Center[1] + amr_length[Ngb_Nodes[node].level + 1] / 2.;
  else
    P[cell].Pos[1] = Ngb_Nodes[node].Center[1] - amr_length[Ngb_Nodes[node].level + 1] / 2.;

#if NUMDIMS >=3
  if(subnode & 4)
    P[cell].Pos[2] = Ngb_Nodes[node].Center[2] + amr_length[Ngb_Nodes[node].level + 1] / 2.;
  else
    P[cell].Pos[2] = Ngb_Nodes[node].Center[2] - amr_length[Ngb_Nodes[node].level + 1] / 2.;
#endif
#endif
}



void amr_split_hydro(int i, double fac)
{
  int k;

  P[i].Mass *= fac;
  SphP[i].OldMass *= fac;
  SphP[i].Energy *= fac;

#ifdef MHD
  for(k = 0; k < 3; k++)
    {
      SphP[i].B[k] = SphP[i].BConserved[k] / (voli + volj);
      SphP[j].B[k] = SphP[i].B[k] + SphP[i].Grad.dB[k][0] * (P[j].Pos[0] - P[i].Pos[0]) + SphP[i].Grad.dB[k][1] * (P[j].Pos[1] - P[i].Pos[1]) + SphP[i].Grad.dB[k][2] * (P[j].Pos[2] - P[i].Pos[2]);    /* extrapolate B to the position of the new cell */

      /* update conserved variables */
      SphP[i].BConserved[k] = SphP[i].B[k] * voli;
      SphP[j].BConserved[k] = SphP[j].B[k] * volj;
    }
#endif
#ifdef MHD_DIVBCLEANING
  SphP[i].Psi = SphP[i].PsiConserved / (voli + volj);
  SphP[j].Psi = SphP[i].Psi + SphP[i].Grad.dPsi[0] * (P[j].Pos[0] - P[i].Pos[0]) + SphP[i].Grad.dPsi[1] * (P[j].Pos[1] - P[i].Pos[1]) + SphP[i].Grad.dPsi[2] * (P[j].Pos[2] - P[i].Pos[2]);

  SphP[i].PsiConserved = SphP[i].Psi * voli;
  SphP[j].PsiConserved = SphP[j].Psi * volj;
#endif

  for(k = 0; k < 3; k++)
    {
      SphP[i].Momentum[k] *= fac;
    }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[i].Entropy *= fac;
#endif
#ifdef USE_SFR
  SphP[i].Sfr *= fac;
#endif
#ifdef BH_THERMALFEEDBACK
  SphP[i].Injected_BH_Energy *= fac;
#endif
#ifdef GFM_WINDS_LOCAL
  SphP[i].WindEnergyReceived *= fac;
#endif

#ifdef MAXSCALARS
  int s;
  for(s = 0; s < N_Scalar; s++) /* Note, the changes in MATERIALS, HIGHRESGASMASS, etc., are treated as part of the Scalars */
    {
      *(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[s].offset_mass) *= fac;
    }
#endif
#ifdef TRACER_FIELD
  SphP[i].ConservedTracer *= fac;
#endif

#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  SphP[i].TotEgyFeed *= fac;
  SphP[i].IntEgyFeed *= fac;
  SphP[i].KinEgyFeed *= fac;
#endif
}
#endif

#ifdef REFINEMENT_MERGE_CELLS
void amr_derefine_cells()
{
  int i;

  int filtered = 0;
  int global_filtered;

  for(i = 0; i < Mesh.Nderefine; i++)
    {
      int cand = Mesh.derefcand[i];
      if(Mesh.Refflag[cand] == 0)
        {
          int cell = amr_derefine_cell(cand);
          Mesh.derefcand[i] = cell;
        }
      else
        {
          filtered++;
        }
    }

  MPI_Allreduce(&Mesh.Nderefine, &Mesh.Derefined, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&filtered, &global_filtered, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("AMR: derefined %d nodes\n", Mesh.Derefined - global_filtered);


  timebin_cleanup_list_of_active_particles(&TimeBinsHydro);
  timebin_cleanup_list_of_active_particles(&TimeBinsGravity);
}

/** This function does the derefinement of a node
 *
 */
int amr_derefine_cell(int no)
{
  int father;
  int j;
  int new_cell;

  father = Ngb_Nodes[no].father;
  new_cell = Ngb_Nodes[no].u.suns[0];


  amr_remove_nextinlevel(no);
  for(j = 0; j < (1 << NUMDIMS); j++)
    {
      amr_remove_nextinlevel(Ngb_Nodes[no].u.suns[j]);
    }

  int subnode = 0;

  if(Ngb_Nodes[no].Center[0] > Ngb_Nodes[father].Center[0])
    subnode += 1;
#if NUMDIMS >=2
  if(Ngb_Nodes[no].Center[1] > Ngb_Nodes[father].Center[1])
    subnode += 2;
#if NUMDIMS >=3
  if(Ngb_Nodes[no].Center[2] > Ngb_Nodes[father].Center[2])
    subnode += 4;
#endif
#endif

  //insert the new cell
  Ngb_Nodes[father].u.suns[subnode] = new_cell;

  P[new_cell].Pos[0] = Ngb_Nodes[no].Center[0];
  P[new_cell].Pos[1] = Ngb_Nodes[no].Center[1];
  P[new_cell].Pos[2] = Ngb_Nodes[no].Center[2];

  amr_set_cell_data(new_cell, father, 1, &Mesh);
  Ngb_Father[new_cell] = father;

  //fix nextnode pointers towards that cell
  int nnext = father;
  int next;
  do
    {
      next = nnext;
      if(next >= Ngb_MaxPart)
        {
          nnext = Ngb_Nodes[next].u.d.nextnode;
        }
      else
        {
          nnext = Ngb_Nextnode[next];
        }
    }
  while(nnext != no);

  if(next >= Ngb_MaxPart)
    {
      Ngb_Nodes[next].u.d.nextnode = new_cell;
    }
  else
    {
      Ngb_Nextnode[next] = new_cell;
    }

  next = father;
  do
    {
      nnext = Ngb_Nodes[next].u.d.nextnode;
      while(nnext != new_cell && nnext != no && nnext >= 0)
        {
          next = nnext;
          if(next >= Ngb_MaxPart)
            {
              nnext = Ngb_Nodes[next].u.d.sibling;
            }
          else
            {
              nnext = Ngb_Nextnode[next];
            }
        }

      if(next >= Ngb_MaxPart && nnext == no)
        {
          Ngb_Nodes[next].u.d.sibling = new_cell;
        }


    }
  while(next >= Ngb_MaxPart && nnext == no);

  int num;
  int cell = new_cell;
  for(num = 0; num < (1 << NUMDIMS); num++)
    {
      assert(cell < Ngb_MaxPart);
      cell = Ngb_Nextnode[cell];
    }
  Ngb_Nextnode[new_cell] = cell;


  //sum up hydro quantities to new cells
#ifdef DG
  dg_sum_hydro(new_cell, 0, new_cell);
#endif

  for(num = 1; num < (1 << NUMDIMS); num++)
    {
      j = Ngb_Nodes[no].u.suns[num];

#ifdef DG
      dg_sum_hydro(new_cell, num, j);
#else
      amr_sum_hydro(new_cell, j);
#endif

#ifdef TRACER_MC
      consider_moving_tracers_local(j, new_cell, 1.0);
#endif

      //mark cell as inactive
      P[j].ID = 0;
      P[j].Mass = 0;
      P[j].Vel[0] = 0;
      P[j].Vel[1] = 0;
      P[j].Vel[2] = 0;

#ifdef AMR_CONNECTIONS
      voronoi_remove_connection(j);
#endif
    }

  return new_cell;
}



void amr_sum_hydro(int new_cell, int i)
{
  P[new_cell].Mass += P[i].Mass;
  SphP[new_cell].Momentum[0] += SphP[i].Momentum[0];
  SphP[new_cell].Momentum[1] += SphP[i].Momentum[1];
  SphP[new_cell].Momentum[2] += SphP[i].Momentum[2];

#ifdef MHD
  SphP[new_cell].BConserved[0] += SphP[i].BConserved[0];
  SphP[new_cell].BConserved[1] += SphP[i].BConserved[1];
  SphP[new_cell].BConserved[2] += SphP[i].BConserved[2];
#ifdef MHD_DIVBCLEANING
  SphP[new_cell].PsiConserved += SphP[i].PsiConserved;
#endif
#endif

#ifndef ISOTHERM_EQS
  SphP[new_cell].Energy += SphP[i].Energy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[new_cell].Entropy += SphP[i].Entropy;
#endif
#ifdef MAXSCALARS
  int s;
  for(s = 0; s < N_Scalar; s++)
    *(MyFloat *) (((char *) (&SphP[new_cell])) + scalar_elements[s].offset_mass) += (*(MyFloat *) (((char *) (&SphP[i])) + scalar_elements[s].offset_mass));
#endif
#ifdef TRACER_FIELD
  SphP[new_cell].ConservedTracer += SphP[i].ConservedTracer;
#endif
#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  SphP[new_cell].TotEgyFeed += SphP[i].TotEgyFeed;
  SphP[new_cell].IntEgyFeed += SphP[i].IntEgyFeed;
  SphP[new_cell].KinEgyFeed += SphP[i].KinEgyFeed;
#endif
#ifdef BH_THERMALFEEDBACK
  SphP[new_cell].Injected_BH_Energy += SphP[i].Injected_BH_Energy;
#endif
#ifdef GFM_WINDS_LOCAL
  SphP[new_cell].WindEnergyReceived += SphP[i].WindEnergyReceived;
#endif
#ifdef TGCHEM
  int m;
  for(m = 0; m < TGCHEM_NUM_ABUNDANCES; m++)
    SphP[p].Abund[m] = (SphP[p].Abund[m] * SphP[p].Volume + SphP[i].Abund[m] * volume[q]) / (SphP[p].Volume + volume[q]);
#endif
}


void amr_remove_nextinlevel(int no)
{
  int next, prev, level;
  //remove node from nextinlevel linked list
  if(no < Ngb_MaxPart)
    {
      prev = Mesh.DP[no].previnlevel;
      next = Mesh.DP[no].nextinlevel;
      level = Mesh.DP[no].level;
    }
  else
    {
      prev = Ngb_Nodes[no].previnlevel;
      next = Ngb_Nodes[no].nextinlevel;
      level = Ngb_Nodes[no].level;
    }



  if(prev >= 0)
    {
      if(prev < Ngb_MaxPart && prev >= 0)
        {
          Mesh.DP[prev].nextinlevel = next;
        }
      else if(prev >= Ngb_MaxPart)
        {
          Ngb_Nodes[prev].nextinlevel = next;
        }

      if(next < Ngb_MaxPart && next >= 0)
        {
          Mesh.DP[next].previnlevel = prev;
        }
      else if(next >= Ngb_MaxPart)
        {
          Ngb_Nodes[next].previnlevel = prev;
        }
    }
  else
    {
      Mesh.lastinlevel[level] = next;
      if(next < Ngb_MaxPart && next >= 0)
        {
          Mesh.DP[next].previnlevel = prev;
        }
      else if(next >= Ngb_MaxPart)
        {
          Ngb_Nodes[next].previnlevel = prev;
        }
    }
}
#endif

#endif
#endif
