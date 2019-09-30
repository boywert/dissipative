/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr_exchange.c
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
#include "../domain.h"

#ifdef AMR

void mesh_setup_exchange()
{
}

void amr_export_node(int p, int dir, int ngb)
{
  assert(ngb >= Ngb_MaxPart);
  assert(ngb < Ngb_FirstNonTopLevelNode);

  if(Ngb_DomainTask[ngb] >= 0 && Ngb_DomainTask[ngb] != ThisTask)
    {
      if(Mesh.NexpListNode == Mesh.MaxNexpListNode)
        {
          Mesh.MaxNexpListNode *= 1.2 + 1;      //FIXME better estiamte
          Mesh.expListNode = myrealloc_movable(Mesh.expListNode, Mesh.MaxNexpListNode * sizeof(list_export_data));
        }
      // write this node to an export list, so we can send it as a ghost to the other side.
      Mesh.expListNode[Mesh.NexpListNode].task = Ngb_DomainTask[ngb];
      Mesh.expListNode[Mesh.NexpListNode].index = p;
      Mesh.NexpListNode++;
    }
  else if(Ngb_DomainTask[ngb] == -1)    //not a top level leaf node
    {
      int d = dir / 2;
      int m = !(dir % 2);

      int j = 0;
      int k = 0;
#if NUMDIMS >= 2
      for(j = 0; j < 2; j++)
#endif
        {
#if NUMDIMS == 3
          for(k = 0; k < 2; k++)
#endif
            {
              int subnode = (j << ((d + 1) % NUMDIMS)) + (k << ((d + 2) % NUMDIMS)) + (m << d);
              int n = Ngb_Nodes[ngb].u.suns[subnode];
              amr_export_node(p, dir, n);
            }
        }

    }
}

static inline unsigned long long ngb_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

void amr_exchange_ghost_nodes(tessellation * T)
{
  int j;
  int index;
  node_exchange *nodeExchOut;
  node_exchange *nodeExchIn;

  int ngrp, recvTask;

  TIMER_START(CPU_AMR_EXCHANGE_NODES)
    //sort our export list
    mysort(T->expListNode, T->NexpListNode, sizeof(list_export_data), compare_list_export_data);


  if(T->NexpListNode > 1)
    {
      int new_nexp = 1;
      for(j = 1; j < T->NexpListNode; j++)
        {
          if(T->expListNode[j].index != T->expListNode[new_nexp - 1].index || T->expListNode[j].task != T->expListNode[new_nexp - 1].task)
            {
              T->expListNode[new_nexp++] = T->expListNode[j];
            }
        }
      T->NexpListNode = new_nexp;
    }

  //now find send/recv offset and count
  for(j = 0; j < NTask; j++)
    T->Node_Send_count[j] = 0;

  for(j = 0; j < T->NexpListNode; j++)
    T->Node_Send_count[T->expListNode[j].task]++;

  //TODO safe some space and just store the index (of expList)

  MPI_Alltoall(T->Node_Send_count, 1, MPI_INT, T->Node_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, T->Node_nimport = 0, T->Node_nexport = 0, T->Node_Recv_offset[0] = 0, T->Node_Send_offset[0] = 0; j < NTask; j++)
    {
      T->Node_nimport += T->Node_Recv_count[j];
      T->Node_nexport += T->Node_Send_count[j];

      if(j > 0)
        {
          T->Node_Send_offset[j] = T->Node_Send_offset[j - 1] + T->Node_Send_count[j - 1];
          T->Node_Recv_offset[j] = T->Node_Recv_offset[j - 1] + T->Node_Recv_count[j - 1];
        }
    }

  assert(T->Node_nexport == T->NexpListNode);

  T->nodes_total = Ngb_MaxNodes + 1 + NTopleaves + T->Node_nimport;
  T->Nghost_nodes = T->Node_nimport;

  Ngb_Nodes += Ngb_MaxPart;
  Ngb_Nodes = (struct NgbNODE *) myrealloc_movable(Ngb_Nodes, sizeof(struct NgbNODE) * T->nodes_total);
  Ngb_Nodes -= Ngb_MaxPart;

  T->original_index = (int *) mymalloc_movable(&T->original_index, "amr_original_index", T->Node_nimport * sizeof(int));
  T->original_task = (int *) mymalloc_movable(&T->original_task, "amr_original_task", T->Node_nimport * sizeof(int));

  T->original_index -= Ngb_MaxPart + Ngb_MaxNodes + NTopleaves;
  T->original_task -= Ngb_MaxPart + Ngb_MaxNodes + NTopleaves;

  //prepare nodes for exchange
  nodeExchOut = (node_exchange *) mymalloc("nodeExchOut", T->Node_nexport * sizeof(node_exchange));
  nodeExchIn = (node_exchange *) mymalloc("nodeExchIn", T->Node_nimport * sizeof(node_exchange));

  for(j = 0; j < T->Node_nexport; j++)
    {
      index = T->expListNode[j].index;

      nodeExchOut[j].originalindex = index;
      nodeExchOut[j].x = Ngb_Nodes[index].Center[0];
      nodeExchOut[j].y = Ngb_Nodes[index].Center[1];
      nodeExchOut[j].z = Ngb_Nodes[index].Center[2];
      nodeExchOut[j].level = Ngb_Nodes[index].level;
    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(T->Node_Send_count[recvTask] > 0 || T->Node_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&nodeExchOut[T->Node_Send_offset[recvTask]], T->Node_Send_count[recvTask]
                           * sizeof(node_exchange), MPI_BYTE, recvTask, TAG_DENS_A,
                           &nodeExchIn[T->Node_Recv_offset[recvTask]], T->Node_Recv_count[recvTask] * sizeof(node_exchange), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  //FIXME MPI Tag



  //remove pseudo particles from Ngb Topleaves
  int no;
  for(j = 0; j < NTopleaves; j++)
    {
      if(DomainTask[j] != ThisTask)
        {
          no = Ngb_DomainNodeIndex[j];

          int i;
          for(i = 0; i < 8; i++)
            {
              Ngb_Nodes[no].u.suns[i] = -1;
            }
        }
    }

  //now insert ghost nodes
  int currTask = 0;
  for(j = 0; j < T->Node_nimport; j++)
    {
      while((currTask < NTask - 1) && (T->Node_Recv_offset[currTask + 1] == j))
        {
          currTask++;
        }
      if(currTask == ThisTask)
        terminate("strange");


      unsigned long long xxb = ngb_double_to_int(((nodeExchIn[j].x - DomainCorner[0]) * DomainInverseLen) + 1.0);
      unsigned long long yyb = ngb_double_to_int(((nodeExchIn[j].y - DomainCorner[1]) * DomainInverseLen) + 1.0);
      unsigned long long zzb = ngb_double_to_int(((nodeExchIn[j].z - DomainCorner[2]) * DomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long) 1) << (52 - 1);
      unsigned char shiftx = (52 - 1);
      unsigned char shifty = (52 - 2);
      unsigned char shiftz = (52 - 3);
      unsigned char levels = 0;

      int no = 0;
      while(TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
        {
          unsigned char subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

          mask >>= 1;
          levels++;

          no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
        }

      no = TopNodes[no].Leaf;

      int th = Ngb_DomainNodeIndex[no];

      signed long long centermask = (0xFFF0000000000000llu) >> levels;

      unsigned char subnode = 0;

      while(1)
        {
          if(th >= Ngb_MaxPart) /* we are dealing with an internal node */
            {
              subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

              centermask >>= 1;
              mask >>= 1;
              levels++;

              if(levels > MAX_TREE_LEVEL)
                {
                  terminate("bad");
                }

              int nn = Ngb_Nodes[th].u.suns[subnode];

              if(nn >= 0)       /* ok, something is in the daughter slot already, need to continue */
                {
                  th = nn;
                }
              else
                {
                  /* here we have found an empty slot where we can attach
                   * the new ghost node.
                   */
                  int new_node = j + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves;

                  Ngb_Nodes[th].u.suns[subnode] = new_node;
                  int k;
                  for(k = 0; k < 8; k++)
                    Ngb_Nodes[new_node].u.suns[k] = -1;

                  for(k = 0; k < 2 * NUMDIMS; k++)
                    Ngb_Nodes[new_node].neighbors[k] = -1;

                  Ngb_Nodes[new_node].father = th;
                  Ngb_Nodes[new_node].Center[0] = nodeExchIn[j].x;
                  Ngb_Nodes[new_node].Center[1] = nodeExchIn[j].y;
                  Ngb_Nodes[new_node].Center[2] = nodeExchIn[j].z;
                  Ngb_Nodes[new_node].level = Ngb_Nodes[th].level + 1;
                  assert(Ngb_Nodes[new_node].level == nodeExchIn[j].level);

                  T->original_index[new_node] = nodeExchIn[j].originalindex;
                  T->original_task[new_node] = currTask;
                  break;        /* done for this node */
                }
            }
          else
            {
              terminate("error");
            }
        }
    }

  //free
  myfree(nodeExchIn);
  myfree(nodeExchOut);

TIMER_STOP(CPU_AMR_EXCHANGE_NODES)}

void amr_exchange_ghost_cells(tessellation * T)
{
  int j;
  int index;
  point_exchange *dpExchOut;
  point_exchange *dpExchIn;

  int ngrp, recvTask;

  TIMER_START(CPU_AMR_EXCHANGE_CELLS)
    //sort our export list
    mysort(T->expList, T->NexpList, sizeof(list_export_data), compare_list_export_data);

  if(T->NexpList > 1)
    {
      int new_nexp = 1;
      for(j = 1; j < T->NexpList; j++)
        {
          if(T->expList[j].index != T->expList[new_nexp - 1].index || T->expList[j].task != T->expList[new_nexp - 1].task)
            {
              T->expList[new_nexp++] = T->expList[j];
            }
        }
      T->NexpList = new_nexp;
    }


  //now find send/recv offset and count
  for(j = 0; j < NTask; j++)
    Mesh_Send_count[j] = 0;

  for(j = 0; j < T->NexpList; j++)
    Mesh_Send_count[T->expList[j].task]++;

  MPI_Alltoall(Mesh_Send_count, 1, MPI_INT, Mesh_Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Mesh_nimport = 0, Mesh_nexport = 0, Mesh_Recv_offset[0] = 0, Mesh_Send_offset[0] = 0; j < NTask; j++)
    {
      Mesh_nimport += Mesh_Recv_count[j];
      Mesh_nexport += Mesh_Send_count[j];

      if(j > 0)
        {
          Mesh_Send_offset[j] = Mesh_Send_offset[j - 1] + Mesh_Send_count[j - 1];
          Mesh_Recv_offset[j] = Mesh_Recv_offset[j - 1] + Mesh_Recv_count[j - 1];
        }
    }

  assert(Mesh_nexport == T->NexpList);

  //prepare DP for exchange
  dpExchIn = (point_exchange *) mymalloc_movable(&dpExchIn, "dpExchIn", Mesh_nimport * sizeof(point_exchange));
  dpExchOut = (point_exchange *) mymalloc_movable(&dpExchIn, "dpExchOut", Mesh_nexport * sizeof(point_exchange));

  for(j = 0; j < Mesh_nexport; j++)
    {
      index = T->expList[j].index;
      dpExchOut[j].x = T->DP[index].x;
      dpExchOut[j].y = T->DP[index].y;
      dpExchOut[j].z = T->DP[index].z;

      dpExchOut[j].originalindex = T->DP[index].index;
      dpExchOut[j].level = T->DP[index].level;
      dpExchOut[j].ID = T->DP[index].ID;
      dpExchOut[j].timebin = T->DP[index].timebin;
    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&dpExchOut[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(point_exchange), MPI_BYTE, recvTask, TAG_DENS_A,
                           &dpExchIn[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(point_exchange), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
  myfree(dpExchOut);

  T->Nghost_dp = Mesh_nimport;

  assert(T->Ndp <= Ngb_MaxPart);

  T->Ndp = Ngb_MaxPart;         //skip some dp

  if(T->Ndp + Mesh_nimport > T->MaxNdp)
    {
      T->MaxNdp = T->Ndp + Mesh_nimport;
      T->DP = myrealloc_movable(T->DP, T->MaxNdp * sizeof(point));
    }

  //now create ghost DP
  int currTask = 0;
  for(j = 0; j < Mesh_nimport; j++)
    {
      while((currTask < NTask - 1) && (Mesh_Recv_offset[currTask + 1] == j))
        {
          currTask++;
        }

      assert(T->Ndp < T->MaxNdp);

      T->DP[T->Ndp].x = dpExchIn[j].x;
      T->DP[T->Ndp].y = dpExchIn[j].y;
      T->DP[T->Ndp].z = dpExchIn[j].z;
      T->DP[T->Ndp].timebin = dpExchIn[j].timebin;
      T->DP[T->Ndp].ID = dpExchIn[j].ID;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
      T->DP[T->Ndp].image_flags = 1;
#endif

      T->DP[T->Ndp].originalindex = dpExchIn[j].originalindex;
      if(currTask == ThisTask)
        {
          assert(0);
          T->DP[T->Ndp].index = T->DP[T->Ndp].originalindex;    //this is a local ghost, we can link directly to P/SphP
        }
      else
        {
          T->DP[T->Ndp].index = j;      //index in primexch list
        }
      T->DP[T->Ndp].task = currTask;

      //now insert dp into tree
      unsigned long long xxb = ngb_double_to_int(((dpExchIn[j].x - DomainCorner[0]) * DomainInverseLen) + 1.0);
      unsigned long long yyb = ngb_double_to_int(((dpExchIn[j].y - DomainCorner[1]) * DomainInverseLen) + 1.0);
      unsigned long long zzb = ngb_double_to_int(((dpExchIn[j].z - DomainCorner[2]) * DomainInverseLen) + 1.0);
      unsigned long long mask = ((unsigned long long) 1) << (52 - 1);
      unsigned char shiftx = (52 - 1);
      unsigned char shifty = (52 - 2);
      unsigned char shiftz = (52 - 3);
      unsigned char levels = 0;

      int no = 0;
      while(TopNodes[no].Daughter >= 0) /* walk down top tree to find correct leaf */
        {
          unsigned char subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

          mask >>= 1;
          levels++;

          no = TopNodes[no].Daughter + TopNodes[no].MortonToPeanoSubnode[subnode];
        }

      no = TopNodes[no].Leaf;

      int th = Ngb_DomainNodeIndex[no];

      signed long long centermask = (0xFFF0000000000000llu) >> levels;

      unsigned char subnode = 0;

      while(1)
        {
          if(th >= Ngb_MaxPart + Ngb_MaxNodes || (th >= Ngb_MaxPart && th < Ngb_FirstNonTopLevelNode))  /* we are dealing with a ghost node or top node */
            {
              subnode = (((unsigned char) ((xxb & mask) >> (shiftx--))) | ((unsigned char) ((yyb & mask) >> (shifty--))) | ((unsigned char) ((zzb & mask) >> (shiftz--))));

              centermask >>= 1;
              mask >>= 1;
              levels++;

              if(levels > MAX_TREE_LEVEL)
                {
                  terminate("bad");
                }

              int nn = Ngb_Nodes[th].u.suns[subnode];

              if(nn >= 0)       /* ok, something is in the daughter slot already, need to continue */
                {
                  if(nn >= T->nodes_total + Ngb_MaxPart && Mesh.DP[nn - T->nodes_total].originalindex == dpExchIn[j].originalindex)
                    {
                      terminate("Error");
                    }

                  th = nn;
                }
              else
                {
                  /* here we have found an empty slot where we can attach
                   * the new particle as a leaf.
                   */
                  Ngb_Nodes[th].u.suns[subnode] = T->Ndp + T->nodes_total;

                  int k;
                  for(k = 0; k < 2 * NUMDIMS; k++)
                    T->DP[T->Ndp].neighbors[k] = -1;

                  T->DP[T->Ndp].father = th;
                  T->DP[T->Ndp].level = Ngb_Nodes[th].level + 1;
                  assert(T->DP[T->Ndp].level == dpExchIn[j].level);
                  break;        /* done for this particle */
                }
            }
          else
            {
              terminate("error %d", th);
            }
        }

      T->Ndp++;
    }


  //free
  myfree(dpExchIn);


  //allocate space for primexch
  PrimExch = (struct primexch *) mymalloc_movable(&(PrimExch), "PrimExch", Mesh_nimport * sizeof(struct primexch));
  GradExch = (struct grad_data *) mymalloc_movable(&(GradExch), "GradExch", Mesh_nimport * sizeof(struct grad_data));
#ifdef TVD_SLOPE_LIMITER
  GradExchUl = (struct grad_data *) mymalloc_movable(&(GradExchUl), "GradExchUl", Mesh_nimport * sizeof(struct grad_data));
#endif

  Mesh.allocated = 2;

TIMER_STOP(CPU_AMR_EXCHANGE_CELLS)}

void exchange_primitive_variables()
{
  TIMER_START(CPU_AMR_EXCHANGE) if(All.TotNumGas == 0)
    return;

  struct primexch *tmpPrimExch;
  int i, j;
  int place;
  int ngrp, recvTask;
  tessellation *T = &Mesh;

  tmpPrimExch = (struct primexch *) mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));

  for(i = 0; i < Mesh_nexport; i++)
    {
      place = T->expList[i].index;

      tmpPrimExch[i].Volume = SphP[place].Volume;
      tmpPrimExch[i].Density = SphP[place].Density;
      tmpPrimExch[i].Pressure = SphP[place].Pressure;

#ifdef DG
      int k, l;

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              tmpPrimExch[i].Weights[l][k] = SphP[place].Weights[l][k];
            }
        }
#endif
#ifdef TGCHEM
      tmpPrimExch[i].Gamma = SphP[place].Gamma;
#endif
#ifdef MHD
      tmpPrimExch[i].B[0] = SphP[place].B[0];
      tmpPrimExch[i].B[1] = SphP[place].B[1];
      tmpPrimExch[i].B[2] = SphP[place].B[2];
#ifdef MHD_DIVBCLEANING
      tmpPrimExch[i].Psi = SphP[place].Psi;
#endif
#endif
#ifdef VARIABLE_GAMMA
      tmpPrimExch[i].GammaC = SphP[place].GammaC;
      tmpPrimExch[i].GammaE = SphP[place].GammaE;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      tmpPrimExch[i].A = SphP[place].A;
      tmpPrimExch[i].Mass = P[place].Mass;
#endif
      tmpPrimExch[i].OldMass = P[place].Mass;
      tmpPrimExch[i].SurfaceArea = SphP[place].SurfaceArea;
      tmpPrimExch[i].ActiveArea = SphP[place].ActiveArea;

      tmpPrimExch[i].TimeBinHydro = P[place].TimeBinHydro;

      tmpPrimExch[i].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        tmpPrimExch[i].Scalars[j] = *(MyFloat *) (((char *) (&SphP[place])) + scalar_elements[j].offset);
#endif
#ifdef TRACER_FIELD
      tmpPrimExch[i].Tracer = SphP[place].Tracer;
#endif

#ifdef FLD
      tmpPrimExch[i].n_gamma = SphP[place].n_gamma;
      tmpPrimExch[i].Kappa_diff = SphP[place].Kappa_diff;
      tmpPrimExch[i].R2 = SphP[place].R2;

#ifdef FLD_CONES
      int cone;
      for(cone = 0; cone < FLD_NCONES; cone++)
        {
          tmpPrimExch[i].gammas[cone] = SphP[place].gammas[cone];
        }
#endif
#endif

      for(j = 0; j < 3; j++)
        {
          tmpPrimExch[i].VelGas[j] = P[place].Vel[j];
          tmpPrimExch[i].Center[j] = SphP[place].Center[j];
          tmpPrimExch[i].VelVertex[j] = SphP[place].VelVertex[j];
        }
    }




  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &PrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpPrimExch);

  TIMER_STOP(CPU_AMR_EXCHANGE) exchange_node_data();
  return;
}

void exchange_primitive_variables_and_gradients()
{
  TIMER_START(CPU_AMR_EXCHANGE) struct grad_data *tmpGradExch;
  struct primexch *tmpPrimExch;
#ifdef TVD_SLOPE_LIMITER
  struct grad_data *tmpGradExchUl;
#endif
  int i, j;
  int place;
  int ngrp, recvTask;
  tessellation *T = &Mesh;

  tmpPrimExch = (struct primexch *) mymalloc("tmpPrimExch", Mesh_nexport * sizeof(struct primexch));
  tmpGradExch = (struct grad_data *) mymalloc("tmpGradExch", Mesh_nexport * sizeof(struct grad_data));
#ifdef TVD_SLOPE_LIMITER
  tmpGradExchUl = (struct grad_data *) mymalloc("tmpGradExchUl", Mesh_nexport * sizeof(struct grad_data));
#endif

  for(i = 0; i < Mesh_nexport; i++)
    {
      place = T->expList[i].index;

      tmpPrimExch[i].Volume = SphP[place].Volume;
      tmpPrimExch[i].Density = SphP[place].Density;
      tmpPrimExch[i].Pressure = SphP[place].Pressure;

#ifdef DG
      int k, l;

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              tmpPrimExch[i].Weights[l][k] = SphP[place].Weights[l][k];
            }
        }
#endif
#ifdef TGCHEM
      tmpPrimExch[i].Gamma = SphP[place].Gamma;
#endif
#ifdef MHD
      tmpPrimExch[i].B[0] = SphP[place].B[0];
      tmpPrimExch[i].B[1] = SphP[place].B[1];
      tmpPrimExch[i].B[2] = SphP[place].B[2];
#ifdef MHD_DIVBCLEANING
      tmpPrimExch[i].Psi = SphP[place].Psi;
#endif
#endif
#ifdef VARIABLE_GAMMA
      tmpPrimExch[i].GammaC = SphP[place].GammaC;
      tmpPrimExch[i].GammaE = SphP[place].GammaE;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      tmpPrimExch[i].A = SphP[place].A;
      tmpPrimExch[i].Mass = P[place].Mass;
#endif
      tmpPrimExch[i].OldMass = P[place].Mass;
      tmpPrimExch[i].SurfaceArea = SphP[place].SurfaceArea;
      tmpPrimExch[i].ActiveArea = SphP[place].ActiveArea;

      tmpPrimExch[i].TimeBinHydro = P[place].TimeBinHydro;

      tmpPrimExch[i].TimeLastPrimUpdate = SphP[place].TimeLastPrimUpdate;

#ifdef FLD
      tmpPrimExch[i].n_gamma = SphP[place].n_gamma;
      tmpPrimExch[i].Kappa_diff = SphP[place].Kappa_diff;
      tmpPrimExch[i].R2 = SphP[place].R2;

#ifdef FLD_CONES
      int cone;
      for(cone = 0; cone < FLD_NCONES; cone++)
        {
          tmpPrimExch[i].gammas[cone] = SphP[place].gammas[cone];
        }
#endif
#endif

#ifdef MAXSCALARS
      for(j = 0; j < N_Scalar; j++)
        tmpPrimExch[i].Scalars[j] = *(MyFloat *) (((char *) (&SphP[place])) + scalar_elements[j].offset);
#endif
#ifdef TRACER_FIELD
      tmpPrimExch[i].Tracer = SphP[place].Tracer;
#endif
      for(j = 0; j < 3; j++)
        {
          tmpPrimExch[i].VelGas[j] = P[place].Vel[j];
          tmpPrimExch[i].Center[j] = SphP[place].Center[j];
          tmpPrimExch[i].VelVertex[j] = SphP[place].VelVertex[j];
        }

      tmpGradExch[i] = SphP[place].Grad;


#ifdef TVD_SLOPE_LIMITER
      tmpGradExchUl[i] = SphP[place].GradUl;
#endif

      tmpPrimExch[i].Csnd = get_sound_speed(place);


    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpPrimExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A,
                           &PrimExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct primexch), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef TVD_SLOPE_LIMITER
              MPI_Sendrecv(&tmpGradExchUl[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExchUl[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
            }
        }
    }

#ifdef TVD_SLOPE_LIMITER
  myfree(tmpGradExchUl);
#endif
  myfree(tmpGradExch);
  myfree(tmpPrimExch);

TIMER_STOP(CPU_AMR_EXCHANGE)}


void exchange_gradients()
{
  TIMER_START(CPU_AMR_EXCHANGE) struct grad_data *tmpGradExch;
#ifdef TVD_SLOPE_LIMITER
  struct grad_data *tmpGradExchUl;
#endif
  int i;
  int place;
  int ngrp, recvTask;
  tessellation *T = &Mesh;

  tmpGradExch = (struct grad_data *) mymalloc("tmpGradExch", Mesh_nexport * sizeof(struct grad_data));
#ifdef TVD_SLOPE_LIMITER
  tmpGradExchUl = (struct grad_data *) mymalloc("tmpGradExchUl", Mesh_nexport * sizeof(struct grad_data));
#endif

  for(i = 0; i < Mesh_nexport; i++)
    {
      place = T->expList[i].index;

      tmpGradExch[i] = SphP[place].Grad;

#ifdef TVD_SLOPE_LIMITER
      tmpGradExchUl[i] = SphP[place].GradUl;
#endif
    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&tmpGradExch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

#ifdef TVD_SLOPE_LIMITER
              MPI_Sendrecv(&tmpGradExchUl[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A,
                           &GradExchUl[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * sizeof(struct grad_data), MPI_BYTE, recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
            }
        }
    }

#ifdef TVD_SLOPE_LIMITER
  myfree(tmpGradExchUl);
#endif
  myfree(tmpGradExch);

TIMER_STOP(CPU_AMR_EXCHANGE)}


void exchange_node_data()
{
  TIMER_START(CPU_AMR_EXCHANGE) amr_node_data *tmpExchExport, *tmpExchImport;
  int i;
  int place;
  int ngrp, recvTask;
  tessellation *T = &Mesh;

  tmpExchExport = (amr_node_data *) mymalloc("tmpExchExport", Mesh.Node_nexport * sizeof(amr_node_data));
  tmpExchImport = (amr_node_data *) mymalloc("tmpExchImport", Mesh.Node_nimport * sizeof(amr_node_data));

  for(i = 0; i < Mesh.Node_nexport; i++)
    {
      place = T->expListNode[i].index;

      tmpExchExport[i] = Ngb_Nodes[place].hydro;
    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh.Node_Send_count[recvTask] > 0 || Mesh.Node_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpExchExport[Mesh.Node_Send_offset[recvTask]], Mesh.Node_Send_count[recvTask]
                           * sizeof(amr_node_data), MPI_BYTE, recvTask, TAG_DENS_A,
                           &tmpExchImport[Mesh.Node_Recv_offset[recvTask]], Mesh.Node_Recv_count[recvTask] * sizeof(amr_node_data), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < Mesh.Node_nimport; i++)
    {
      int node = i + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves;

      Ngb_Nodes[node].hydro = tmpExchImport[i];
    }

  myfree(tmpExchImport);
  myfree(tmpExchExport);

  TIMER_STOP(CPU_AMR_EXCHANGE) return;
}

int compare_list_export_data(const void *a, const void *b)
{
  if(((list_export_data *) a)->task < ((list_export_data *) b)->task)
    return -1;

  if(((list_export_data *) a)->task > ((list_export_data *) b)->task)
    return +1;

  if(((list_export_data *) a)->index < ((list_export_data *) b)->index)
    return -1;

  if(((list_export_data *) a)->index > ((list_export_data *) b)->index)
    return +1;

  return 0;
}


#endif //AMR
