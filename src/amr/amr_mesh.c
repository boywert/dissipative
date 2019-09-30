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



/*! This function is the main driver of the amr mesh
 * construction. It also does the ghost cell search/exchange
 *
 */
void create_mesh()
{
  TIMER_START(CPU_AMR) if(All.TotNumGas == 0)
    {
      return;
    }

  mpi_printf("AMR: building amr mesh\n");

  if(Mesh.allocated > 1)
    {
#ifdef TVD_SLOPE_LIMITER
      myfree(GradExchUl);
#endif
      myfree(GradExch);
      myfree(PrimExch);

      myfree(Mesh.original_task + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves);
      myfree(Mesh.original_index + Ngb_MaxPart + Ngb_MaxNodes + NTopleaves);
      Mesh.allocated = 1;
    }

  //remove old ghost points
  Mesh.Ndp = Mesh.LastDP;

  Mesh.Nghost_dp = 0;
  Mesh.Nghost_nodes = 0;
  Mesh.nodes_total = Ngb_MaxNodes + 1 + NTopleaves;

  amr_exchange_ghost_nodes(&Mesh);
  amr_exchange_ghost_cells(&Mesh);

  TIMER_START(CPU_AMR_LINK_NGB);
  amr_link_ngb(Ngb_MaxPart);
  TIMER_STOP(CPU_AMR_LINK_NGB);

  amr_create_faces(&Mesh);

#ifdef AMR_CONNECTIONS
  voronoi_update_connectivity(&Mesh);
#endif

  MPI_Allreduce(MPI_IN_PLACE, &Mesh.maxlevel, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &Mesh.minlevel, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

TIMER_STOP(CPU_AMR)}

/* This function creates all interfaces
 */
int amr_create_faces(tessellation * T)
{
  int i, j, k;
  int d;

  TIMER_START(CPU_AMR_MESH) Mesh.Nvf = 0;

  int allFaces = 0;

  int idx;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;


#ifdef AMR_REMAP
      if(P[i].ID < amr_IDNew)
        {
          continue;
        }
#endif

      for(d = 0; d < NUMDIMS; d++)
        {
          int ngb = T->DP[i].neighbors[2 * d];

          if(ngb < 0)
            terminate("broken amr mesh");

          if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart + T->nodes_total)  // a node
            {
              j = 0;
              k = 0;
#if NUMDIMS >=2
              for(j = 0; j < 2; j++)
#endif
                {
#if NUMDIMS == 3
                  for(k = 0; k < 2; k++)
#endif
                    {
                      int subnode = (j << ((d + 1) % NUMDIMS)) + (k << ((d + 2) % NUMDIMS)) + (1 << d);
                      int p2 = Ngb_Nodes[ngb].u.suns[subnode];

                      assert(p2 < Ngb_MaxPart || p2 >= Ngb_MaxPart + T->nodes_total);

                      if(p2 >= Ngb_MaxPart + T->nodes_total)
                        {
                          p2 -= T->nodes_total;
                        }
                      amr_create_face(i, p2, 2 * d, T);
                    }
                }
            }
          else
            {
              int p2 = ngb;

              assert(p2 < Ngb_MaxPart || p2 >= Ngb_MaxPart + T->nodes_total);

              if(p2 >= Ngb_MaxPart + T->nodes_total)
                {
                  p2 -= T->nodes_total;
                }
              amr_create_face(i, p2, 2 * d, T);
            }




          ngb = T->DP[i].neighbors[2 * d + 1];
          if(ngb < 0)
            terminate("broken amr mesh");


          if(ngb > Ngb_MaxPart + Ngb_MaxNodes && ngb < Ngb_MaxPart + T->nodes_total)    //a ghost node
            {
              j = 0;
              k = 0;
#if NUMDIMS >= 2
              for(j = 0; j < 2; j++)
#endif
                {
#if NUMDIMS == 3
                  for(k = 0; k < 2; k++)
#endif
                    {
                      int subnode = (j << ((d + 1) % NUMDIMS)) + (k << ((d + 2) % NUMDIMS));
                      int p2 = Ngb_Nodes[ngb].u.suns[subnode];

                      assert(p2 >= Ngb_MaxPart + T->nodes_total);
                      p2 -= T->nodes_total;
                      amr_create_face(i, p2, 2 * d + 1, T);
                    }
                }
            }
          else if(ngb >= Ngb_MaxPart + T->nodes_total)              // a ghost particle
            {
              int p2 = ngb;

              assert(p2 >= Ngb_MaxPart + T->nodes_total);
              p2 -= T->nodes_total;
              amr_create_face(i, p2, 2 * d + 1, T);
            }

#ifndef FORCE_EQUAL_TIMESTEPS
          else if(ngb < Ngb_MaxPart && T->DP[ngb].timebin > All.HighestActiveTimeBin)
            {
              int p2 = ngb;
              amr_create_face(i, p2, 2 * d + 1, T);
            }
          else if(ngb >= Ngb_MaxPart && ngb < Ngb_MaxPart + T->nodes_total)  // a node
            {
              j = 0;
              k = 0;
#if NUMDIMS >=2
              for(j = 0; j < 2; j++)
#endif
                {
#if NUMDIMS == 3
                  for(k = 0; k < 2; k++)
#endif
                    {
                      int subnode = (j << ((d + 1) % NUMDIMS)) + (k << ((d + 2) % NUMDIMS));
                      int p2 = Ngb_Nodes[ngb].u.suns[subnode];

                      assert(p2 < Ngb_MaxPart || p2 >= Ngb_MaxPart + T->nodes_total);

                      if(p2 >= Ngb_MaxPart + T->nodes_total)
                        {
                          p2 -= T->nodes_total;
                        }
                      if(T->DP[p2].timebin > All.HighestActiveTimeBin)
                        {
                          amr_create_face(i, p2, 2 * d, T);
                        }
                    }
                }
            }
#endif


#if  defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
          if(ngb == i)
            {
              amr_create_face(i, i, 2 * d + 1, T);
            }
#endif
        }
    }

  MPI_Reduce(&(T->Nvf), &allFaces, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  mpi_printf("AMR: %d faces created\n", allFaces);
  TIMER_STOP(CPU_AMR_MESH) return 0;
}

void amr_create_face(int p1, int ngb, int direction, tessellation * T)
{
  int p2 = ngb;

  assert(abs(T->DP[p1].level - T->DP[ngb].level) <= 1);

  if(T->DP[p1].task == ThisTask)
    {
      int particle = T->DP[p1].index;
      if(P[particle].Ti_Current != All.Ti_Current)
        {
          drift_particle(particle, All.Ti_Current);
        }
    }

  if(T->DP[ngb].task == ThisTask)
    {
      int particle = T->DP[ngb].index;
      if(P[particle].Ti_Current != All.Ti_Current)
        {
          drift_particle(particle, All.Ti_Current);
        }
    }

  switch (direction)
    {
    case 0:
      if(T->DP[p1].x < T->DP[ngb].x)
        {
          if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
            {
              amr_create_face(ngb, p1, direction + 1, T);
            }
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].x -= boxSize_X;
        }
#ifdef REFLECTIVE_X
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].x = -T->DP[p2].x;
          T->DP[p2].image_flags = 1 << 1;
#if (REFLECTIVE_X == 2)
          T->DP[p2].image_flags |= OUTFLOW_X;
#endif
        }
#endif
      break;

    case 1:
      if(T->DP[p1].x > T->DP[ngb].x)
        {
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].x += boxSize_X;
        }
#ifdef REFLECTIVE_X
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].x = -T->DP[p1].x + 2 * boxSize_X;
          T->DP[p2].image_flags = 1 << 2;
#if (REFLECTIVE_X == 2)
          T->DP[p2].image_flags |= OUTFLOW_X;
#endif
        }
#endif
      break;

    case 2:
      if(T->DP[p1].y < T->DP[ngb].y)
        {
          if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
            {
              amr_create_face(ngb, p1, direction + 1, T);
            }
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].y -= boxSize_Y;
        }
#ifdef REFLECTIVE_Y
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].y = -T->DP[p1].y;
          T->DP[p2].image_flags = 1 << 1 * 3;
#if (REFLECTIVE_Y == 2)
          T->DP[p2].image_flags |= OUTFLOW_Y;
#endif
        }
#endif
      break;

    case 3:
      if(T->DP[p1].y > T->DP[ngb].y)
        {
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].y += boxSize_Y;
        }
#ifdef REFLECTIVE_Y
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].y = -T->DP[p1].y + 2 * boxSize_Y;
          T->DP[p2].image_flags = 1 << 2 * 3;
#if (REFLECTIVE_Y == 2)
          T->DP[p2].image_flags |= OUTFLOW_Y;
#endif
        }
#endif
      break;

    case 4:
      if(T->DP[p1].z < T->DP[ngb].z)
        {
          if(ngb < Ngb_MaxPart + Ngb_MaxNodes)
            {
              amr_create_face(ngb, p1, direction + 1, T);
            }
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].z -= boxSize_Z;
        }
#ifdef REFLECTIVE_Z
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].z = -T->DP[p1].z;
          T->DP[p2].image_flags = 1 << 1 * 9;
#if (REFLECTIVE_Z == 2)
          T->DP[p2].image_flags |= OUTFLOW_Z;
#endif
        }
#endif
      break;

    case 5:
      if(T->DP[p1].z > T->DP[ngb].z)
        {
          p2 = amr_copy_dp(ngb, T);
          T->DP[p2].z += boxSize_Z;
        }
#ifdef REFLECTIVE_Z
      else if(p1 == ngb)
        {
          p2 = amr_copy_dp(p1, T);
          T->DP[p2].z = -T->DP[p1].z + 2 * boxSize_Z;
          T->DP[p2].image_flags = 1 << 2 * 9;
#if (REFLECTIVE_Z == 2)
          T->DP[p2].image_flags |= OUTFLOW_Z;
#endif
        }
#endif
      break;
    }

  if(T->Nvf >= T->MaxNvf)
    {
      T->MaxNvf = T->MaxNvf * 1.5 + 1;  //FIXME better estimate
      T->VF = myrealloc_movable(T->VF, T->MaxNvf * sizeof(face));
    }

  T->VF[T->Nvf].p1 = p1;
  T->VF[T->Nvf].p2 = p2;

  int maxlevel = imax(T->DP[p1].level, T->DP[p2].level);

  T->VF[T->Nvf].area = amr_area[maxlevel];

  int pmax;
  int d;
  if(direction % 2 == 0)
    d = -1;
  else
    d = 1;

  if(T->DP[p2].level > T->DP[p1].level)
    {
      pmax = p2;
      d = (-1) * d;
    }
  else
    {
      pmax = p1;
    }



  T->VF[T->Nvf].cx = T->DP[pmax].x;
  T->VF[T->Nvf].cy = T->DP[pmax].y;
  T->VF[T->Nvf].cz = T->DP[pmax].z;

  if(direction / 2 == 0)
    T->VF[T->Nvf].cx += d * amr_length[maxlevel + 1];
  else if(direction / 2 == 1)
    T->VF[T->Nvf].cy += d * amr_length[maxlevel + 1];
  else if(direction / 2 == 2)
    T->VF[T->Nvf].cz += d * amr_length[maxlevel + 1];

  T->Nvf++;
}

int amr_copy_dp(int i, tessellation * T)
{
  if(T->Ndp >= T->MaxNdp)
    {
      T->MaxNdp = T->MaxNdp * T->Indi.AllocFacNdp + 1;  //FIXME better estimate
      T->DP = myrealloc_movable(T->DP, T->MaxNdp * sizeof(point));
    }

  int new = T->Ndp;

  memcpy(T->DP + new, T->DP + i, sizeof(point));

  int j;
  for(j = 0; j < 2 * NUMDIMS; j++)
    {
      T->DP[new].neighbors[j] = -1;
    }

  if(T->DP[new].task == ThisTask)
    {
      T->DP[new].index += NumGas;
    }

  T->Ndp++;

  return new;
}

#endif
