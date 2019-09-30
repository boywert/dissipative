/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/amr/amr.c
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

tessellation Mesh;

double amr_length[AMR_MAX_REFLEVEL + 1];
double amr_area[AMR_MAX_REFLEVEL + 1];
double amr_volume[AMR_MAX_REFLEVEL + 1];

int amr_amrtree = 1;

#if defined(LONG_X) || defined(LONG_Y) ||defined(LONG_Z)
int amr_last_long_level[3];
#endif


void amr_init()
{
  int i = 0;

  Mesh.allocated = 0;

#ifdef REFINEMENT
  Mesh.Nrefine = 0;
  Mesh.Nderefine = 0;
#endif

  double len = fmax(boxSize_X, fmax(boxSize_Y, boxSize_Z));

#if defined(LONG_X) || defined(LONG_Y) ||defined(LONG_Z)
  amr_last_long_level[0] = 0;
  amr_last_long_level[1] = 0;
  amr_last_long_level[2] = 0;

  double long_x = boxSize_X / len;

  while(long_x < 0.99)
    {
      amr_last_long_level[0]++;

      long_x *= 2.;
    }

  double long_y = boxSize_Y / len;

  while(long_y < 0.99)
    {
      amr_last_long_level[1]++;

      long_y *= 2.;
    }

  double long_z = boxSize_Z / len;

  while(long_z < 0.99)
    {
      amr_last_long_level[2]++;

      long_z *= 2.;
    }
#endif

  amr_length[0] = len;

#ifdef ONEDIMS
  amr_area[0] = 1.;
  amr_volume[0] = amr_length[0];
#else
#ifdef TWODIMS
  amr_area[0] = amr_length[0];
  amr_volume[0] = amr_length[0] * amr_length[0];
#else
  amr_area[0] = amr_length[0] * amr_length[0];
  amr_volume[0] = amr_area[0] * amr_length[0];
#endif
#endif

  for(i = 1; i < AMR_MAX_REFLEVEL + 1; i++)
    {
      amr_length[i] = amr_length[i - 1] / 2;
#ifdef ONEDIMS
      amr_area[i] = 1.;
      amr_volume[i] = amr_length[i];
#else
#ifdef TWODIMS
      amr_area[i] = amr_length[i];
      amr_volume[i] = amr_length[i] * amr_length[i];
#else
      amr_area[i] = amr_length[i] * amr_length[i];
      amr_volume[i] = amr_area[i] * amr_length[i];
#endif
#endif
    }

  if(All.MinRefLevel < 2)
    terminate("Please set the parameter MinRefLevel to at least 2");
}

int amr_check_domain_decomposition()
{
#ifdef AMR_REDUCE_DOMAIN_DECOMPOISTION

#ifdef OUTPUT_EVERY_STEP
  return 1;
#endif

  if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      return 1;
    }

  return 0;

#else

  return 1;

#endif
}

void amr_set_cell_data(int cell, int father, int lists, tessellation * T)
{
  T->DP[cell].father = father;
  T->DP[cell].index = cell;
  T->DP[cell].task = ThisTask;
  T->DP[cell].originalindex = cell;
  T->DP[cell].x = P[cell].Pos[0];
  T->DP[cell].y = P[cell].Pos[1];
  T->DP[cell].z = P[cell].Pos[2];
  T->DP[cell].ID = P[cell].ID;
  T->DP[cell].timebin = P[cell].TimeBinHydro;
  T->DP[cell].level = Ngb_Nodes[father].level + 1;

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  T->DP[cell].image_flags = 1;
#endif

  assert(T->DP[cell].level >= All.MinRefLevel && T->DP[cell].level <= All.MaxRefLevel);

  if(T->minlevel > T->DP[cell].level)
    T->minlevel = T->DP[cell].level;

  if(T->maxlevel < T->DP[cell].level)
    T->maxlevel = T->DP[cell].level;

  SphP[cell].Volume = amr_volume[T->DP[cell].level];
  SphP[cell].SurfaceArea = 2 * NUMDIMS * amr_area[T->DP[cell].level];
  SphP[cell].ActiveArea = 2 * NUMDIMS * amr_area[T->DP[cell].level];    //TODO this is just a rough estimate

  SphP[cell].Level = T->DP[cell].level;

  SphP[cell].Center[0] = P[cell].Pos[0];
  SphP[cell].Center[1] = P[cell].Pos[1];
  SphP[cell].Center[2] = P[cell].Pos[2];

  if(lists)
    {
      int last = T->lastinlevel[T->DP[cell].level];
      T->DP[cell].nextinlevel = last;
      T->DP[cell].previnlevel = -1;

      if(last >= Ngb_MaxPart)
        {
          Ngb_Nodes[last].previnlevel = cell;
        }
      else if(last >= 0)
        {
          T->DP[last].previnlevel = cell;
        }

      T->lastinlevel[T->DP[cell].level] = cell;
    }
}


void amr_set_node_data(int target, int parent, int subnode, int lists, tessellation * T)
{
  int level = Ngb_Nodes[target].level;
  assert(level < All.MaxRefLevel);
  if(parent != -1)
    {
      assert(level == Ngb_Nodes[parent].level + 1);
      assert(level > 0);
    }
  else
    {
      assert(level == 0);
    }

  Ngb_Nodes[target].father = parent;
  //set center of target
  if(level != 0)
    {
      if(subnode & 1)
        {
          Ngb_Nodes[target].Center[0] = Ngb_Nodes[parent].Center[0] + amr_length[level + 1];
        }
      else
        {
          Ngb_Nodes[target].Center[0] = Ngb_Nodes[parent].Center[0] - amr_length[level + 1];
        }
#if NUMDIMS >= 2
      if(subnode & 2)
        {
          Ngb_Nodes[target].Center[1] = Ngb_Nodes[parent].Center[1] + amr_length[level + 1];
        }
      else
        {
          Ngb_Nodes[target].Center[1] = Ngb_Nodes[parent].Center[1] - amr_length[level + 1];
        }
#if NUMDIMS == 3
      if(subnode & 4)
        {
          Ngb_Nodes[target].Center[2] = Ngb_Nodes[parent].Center[2] + amr_length[level + 1];
        }
      else
        {
          Ngb_Nodes[target].Center[2] = Ngb_Nodes[parent].Center[2] - amr_length[level + 1];
        }
#else
      Ngb_Nodes[target].Center[2] = 0.;
#endif
#else
      Ngb_Nodes[target].Center[1] = 0.;
      Ngb_Nodes[target].Center[2] = 0.;
#endif
    }

  if(lists)
    {
      if(subnode < (1 << NUMDIMS))
        {
          int last = Mesh.lastinlevel[Ngb_Nodes[target].level];
          Ngb_Nodes[target].nextinlevel = last;
          Ngb_Nodes[target].previnlevel = -1;

          if(last >= Ngb_MaxPart)
            {
              Ngb_Nodes[last].previnlevel = target;
            }
          else if(last >= 0)
            {
              Mesh.DP[last].previnlevel = target;
            }

          Mesh.lastinlevel[Ngb_Nodes[target].level] = target;
        }
    }
}


void amr_free_mesh(tessellation * T)
{
  if(T->allocated)
    {
      if(T->allocated > 1)
        {
#ifdef TVD_SLOPE_LIMITER
          myfree(GradExchUl);
#endif
          myfree(GradExch);
          myfree(PrimExch);

          myfree(T->original_task + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves);
          myfree(T->original_index + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves);
        }

      myfree(T->Node_Recv_offset);
      myfree(T->Node_Recv_count);
      myfree(T->Node_Send_offset);
      myfree(T->Node_Send_count);


      myfree(T->expListNode);
      myfree(T->expList);
      myfree(T->VF);
      myfree(T->DP);

    }
  T->allocated = 0;

}

void amr_allocate_mesh(tessellation * T)
{
  if(T->allocated)
    {
      terminate("trying to allocate an already allocated mesh.");
    }

  int expectedGhosts = All.MaxPartSph * 2.0;    // TODO: use better value

  T->MaxNdp = T->Indi.AllocFacNdp * NumGas + expectedGhosts;    //add space for ghosts
  T->DP = (point *) mymalloc_movable(&(T->DP), "DP", (T->MaxNdp) * sizeof(point));
  T->Ndp = 0;

  T->MaxNvf = T->Indi.AllocFacNvf * NumGas;
  T->VF = (face *) mymalloc_movable(&(T->VF), "VF", T->MaxNvf * sizeof(face));
  T->Nvf = 0;

  T->MaxNexpList = expectedGhosts;
  T->expList = (list_export_data *) mymalloc_movable(&T->expList, "expList", T->MaxNexpList * sizeof(list_export_data));
  T->NexpList = 0;

  T->MaxNexpListNode = expectedGhosts;
  T->expListNode = (list_export_data *) mymalloc_movable(&T->expListNode, "expListNode", T->MaxNexpListNode * sizeof(list_export_data));
  T->NexpListNode = 0;

  T->Node_Send_count = (int *) mymalloc_movable(&T->Node_Send_count, "Node_Send_count", sizeof(int) * NTask);
  T->Node_Send_offset = (int *) mymalloc_movable(&T->Node_Send_offset, "Node_Send_offset", sizeof(int) * NTask);
  T->Node_Recv_count = (int *) mymalloc_movable(&T->Node_Recv_count, "Node_Recv_count", sizeof(int) * NTask);
  T->Node_Recv_offset = (int *) mymalloc_movable(&T->Node_Recv_offset, "Node_Recv_offset", sizeof(int) * NTask);

#ifdef REFINEMENT
  Mesh.Nrefine = 0;
  Mesh.Nderefine = 0;
#endif

  T->allocated = 1;
}

void amr_initmesh()
{
  int i = 0;

  Mesh.maxlevel = 0;
  Mesh.minlevel = INT_MAX;

  for(i = 0; i < AMR_MAX_REFLEVEL + 1; i++)
    {
      Mesh.lastinlevel[i] = -1;
    }
}

void amr_reset_node_data(int no)
{
  Ngb_Nodes[no].hydro.mass = 0.;
  Ngb_Nodes[no].hydro.momentum[0] = 0.;
  Ngb_Nodes[no].hydro.momentum[1] = 0.;
  Ngb_Nodes[no].hydro.momentum[2] = 0.;
  Ngb_Nodes[no].hydro.etherm = 0.;

#ifdef DG
  int l, k;

  for(l = 0; l < Nof_base_functions; l++)
    {
      for(k = 0; k < 5; k++)
        {
          Ngb_Nodes[no].hydro.Weights[l][k] = 0;
        }
    }
#endif
}


int amr_get_subnode(int father, int node)
{
  double node_x = Ngb_Nodes[father].Center[0];
  double node_y = Ngb_Nodes[father].Center[1];
  double node_z = Ngb_Nodes[father].Center[2];

  double sun_x;
  double sun_y;
  double sun_z;

  if(node < Ngb_MaxPart)
    {
      sun_x = P[node].Pos[0];
      sun_y = P[node].Pos[1];
      sun_z = P[node].Pos[2];
    }
  else if(node < Ngb_MaxPart + Ngb_MaxNodes)
    {
      sun_x = Ngb_Nodes[node].Center[0];
      sun_y = Ngb_Nodes[node].Center[1];
      sun_z = Ngb_Nodes[node].Center[2];
    }
  else
    {
      terminate("ERROR in amr_get_subnode: invalid node index!\n");
    }

  int result = 0;

  if(sun_x > node_x)
    result += 1;

#if (NUMDIMS > 1)
  if(sun_y > node_y)
    result += 2;

#if (NUMDIMS > 2)
  if(sun_z > node_z)
    result += 4;
#endif
#endif

  return result;

}

void amr_accumulate_node_data(int no, int p)
{
  if(p < Ngb_MaxPart)
    {
      Ngb_Nodes[no].hydro.mass += P[p].Mass;
      Ngb_Nodes[no].hydro.momentum[0] += SphP[p].Momentum[0];
      Ngb_Nodes[no].hydro.momentum[1] += SphP[p].Momentum[1];
      Ngb_Nodes[no].hydro.momentum[2] += SphP[p].Momentum[2];
      Ngb_Nodes[no].hydro.etherm += SphP[p].Utherm * P[p].Mass;
    }
  else
    {
      Ngb_Nodes[no].hydro.mass += Ngb_Nodes[p].hydro.mass;
      Ngb_Nodes[no].hydro.momentum[0] += Ngb_Nodes[p].hydro.momentum[0];
      Ngb_Nodes[no].hydro.momentum[1] += Ngb_Nodes[p].hydro.momentum[1];
      Ngb_Nodes[no].hydro.momentum[2] += Ngb_Nodes[p].hydro.momentum[2];
      Ngb_Nodes[no].hydro.etherm += Ngb_Nodes[p].hydro.etherm;
    }

#ifdef DG

  //Do a L2 projection
  //

  double (*P_X)[NOF_BASE_FUNCTIONS];

  switch (amr_get_subnode(no, p))
    {
    case 0:
      P_X = P_A;
      break;
    case 1:
      P_X = P_B;
      break;
    case 2:
      P_X = P_C;
      break;
    case 3:
      P_X = P_D;
      break;
#ifndef TWODIMS
    case 4:
      P_X = P_E;
      break;
    case 5:
      P_X = P_F;
      break;
    case 6:
      P_X = P_G;
      break;
    case 7:
      P_X = P_H;
      break;
#endif
    default:
      terminate("ERROR in amr_accumulate_node_data: invalid sun position!\n");
    }

  int l, j;

  if(p < Ngb_MaxPart)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(j = 0; j < Nof_base_functions; j++)
            {
              Ngb_Nodes[no].hydro.Weights[l][0] += 1. / DG_PROJ_NORM * SphP[p].Weights[j][0] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][1] += 1. / DG_PROJ_NORM * SphP[p].Weights[j][1] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][2] += 1. / DG_PROJ_NORM * SphP[p].Weights[j][2] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][3] += 1. / DG_PROJ_NORM * SphP[p].Weights[j][3] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][4] += 1. / DG_PROJ_NORM * SphP[p].Weights[j][4] * P_X[l][j];
            }
        }
    }
  else
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(j = 0; j < Nof_base_functions; j++)
            {
              Ngb_Nodes[no].hydro.Weights[l][0] += 1. / DG_PROJ_NORM * Ngb_Nodes[p].hydro.Weights[j][0] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][1] += 1. / DG_PROJ_NORM * Ngb_Nodes[p].hydro.Weights[j][1] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][2] += 1. / DG_PROJ_NORM * Ngb_Nodes[p].hydro.Weights[j][2] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][3] += 1. / DG_PROJ_NORM * Ngb_Nodes[p].hydro.Weights[j][3] * P_X[l][j];
              Ngb_Nodes[no].hydro.Weights[l][4] += 1. / DG_PROJ_NORM * Ngb_Nodes[p].hydro.Weights[j][4] * P_X[l][j];
            }
        }
    }
#endif
}

void amr_treemodifylength(int delta_NgbMaxPart)
{
  mpi_printf("ALLOCATE: Need to adjust AMR data because Ngb_MaxPart needs to grow by %d\n", delta_NgbMaxPart);

  int i;
  for(i = 0; i < Mesh.Ndp; i++)
    {
      if(Mesh.DP[i].father >= Ngb_MaxPart)
        Mesh.DP[i].father += delta_NgbMaxPart;

      if(Mesh.DP[i].nextinlevel >= Ngb_MaxPart)
        Mesh.DP[i].nextinlevel += delta_NgbMaxPart;

      if(Mesh.DP[i].previnlevel >= Ngb_MaxPart)
        Mesh.DP[i].previnlevel += delta_NgbMaxPart;

      int k;
      for(k = 0; k < 2 * NUMDIMS; k++)
        {
          if(Mesh.DP[i].neighbors[k] >= Ngb_MaxPart)
            Mesh.DP[i].neighbors[k] += delta_NgbMaxPart;
        }
    }

  for(i = 0; i < AMR_MAX_REFLEVEL + 1; i++)
    {
      if(Mesh.lastinlevel[i] >= Ngb_MaxPart)
        Mesh.lastinlevel[i] += delta_NgbMaxPart;
    }

  Mesh.original_index -= delta_NgbMaxPart;
  Mesh.original_task -= delta_NgbMaxPart;

#ifdef REFINEMENT
  for(i = 0; i < Mesh.Nderefine; i++)
    {
      Mesh.derefcand[i] += delta_NgbMaxPart;
    }

  if(Mesh.allocated > 2)
    {
      Mesh.Refflag = myrealloc_movable(Mesh.Refflag, (Ngb_MaxPart + delta_NgbMaxPart + Mesh.nodes_total + Mesh.Nghost_dp) * sizeof(int));
      memmove(&Mesh.Refflag[Ngb_MaxPart + delta_NgbMaxPart], &Mesh.Refflag[Ngb_MaxPart], sizeof(int) * Mesh.nodes_total + Mesh.Nghost_dp);
    }
#endif
}

void amr_treerealloc(int delta_Nodes)
{
  Mesh.nodes_total += delta_Nodes;

  Mesh.original_index -= delta_Nodes;
  Mesh.original_task -= delta_Nodes;
}

double get_cell_radius(int i)
{
  return 0.5 * amr_length[Mesh.DP[i].level];
}

double nearest_x(double d)
{
#if !defined(REFLECTIVE_X)
  if(d < -boxHalf_X)
    d += boxSize_X;
  if(d > boxHalf_X)
    d -= boxSize_X;
#endif
  return d;
}

double nearest_y(double d)
{
#if !defined(REFLECTIVE_Y)
  if(d < -boxHalf_Y)
    d += boxSize_Y;
  if(d > boxHalf_Y)
    d -= boxSize_Y;
#endif
  return d;
}

double nearest_z(double d)
{
#if !defined(REFLECTIVE_Z)
  if(d < -boxHalf_Z)
    d += boxSize_Z;
  if(d > boxHalf_Z)
    d -= boxSize_Z;
#endif
  return d;
}

int face_get_normals(tessellation * T, int i, struct geometry *geom)
{
  int li, ri;
  double mm;

  face *VF = T->VF;
  point *DP = T->DP;

  li = DP[VF[i].p1].index;
  ri = DP[VF[i].p2].index;

  if(li < 0 || ri < 0)
    return -1;


  /* center of face */
  geom->cx = VF[i].cx;
  geom->cy = VF[i].cy;
  geom->cz = VF[i].cz;


  double dx = fabs(DP[VF[i].p2].x - DP[VF[i].p1].x);
  double dy = fabs(DP[VF[i].p2].y - DP[VF[i].p1].y);
  double dz = fabs(DP[VF[i].p2].z - DP[VF[i].p1].z);
  /* normal vector pointing to "right" state */
  geom->nx = 0.;
  geom->ny = 0.;
  geom->nz = 0.;

  double max = fmax(dx, fmax(dy, dz));

  if(dx == max)
    geom->nx = DP[VF[i].p2].x - DP[VF[i].p1].x;
  else if(dy == max)
    geom->ny = DP[VF[i].p2].y - DP[VF[i].p1].y;
  else if(dz == max)
    geom->nz = DP[VF[i].p2].z - DP[VF[i].p1].z;

  geom->nn = sqrt(geom->nx * geom->nx + geom->ny * geom->ny + geom->nz * geom->nz);
  geom->nx /= geom->nn;
  geom->ny /= geom->nn;
  geom->nz /= geom->nn;

  /* need an ortonormal basis */
  if(geom->nx != 0 || geom->ny != 0)
    {
      geom->mx = -geom->ny;
      geom->my = geom->nx;
      geom->mz = 0;
    }
  else
    {
      geom->mx = 1;
      geom->my = 0;
      geom->mz = 0;
    }

  mm = sqrt(geom->mx * geom->mx + geom->my * geom->my + geom->mz * geom->mz);
  geom->mx /= mm;
  geom->my /= mm;
  geom->mz /= mm;

  geom->px = geom->ny * geom->mz - geom->nz * geom->my;
  geom->py = geom->nz * geom->mx - geom->nx * geom->mz;
  geom->pz = geom->nx * geom->my - geom->ny * geom->mx;

  return 0;
}

#endif /* AMR */
