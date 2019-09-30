/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/voronoi_derefinement_pairs.c
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

#include "allvars.h"
#include "proto.h"
#include "voronoi.h"

#ifdef VORONOI

#if defined(REFINEMENT_MERGE_CELLS) && !defined(ONEDIMS) && defined(REFINEMENT_MERGE_PAIRS)

#ifdef GFM_CHEMTAGS
#error "does presently not work with GFM_CHEMTAGS"
#endif

int do_derefinements(void);
void voronoi_derefinement_find_partner(int i);
void voronoi_derefinement_merge_pair(int i, int j);
void voronoi_derefinement_update_partner_index(int i);

#ifdef REFINEMENT_SPLIT_CELLS
extern char *FlagDoNotRefine;
#endif

int do_derefinements(void)
{
  int idx, i, count = 0, countall = 0, pcount = 0, pcountall = 0, fcount = 0, fcountall = 0;
  int lcount = 0, lcountall = 0;

  CPU_Step[CPU_MISC] += measure_time();

  /* first, check whether we have cells to derefine */
  int NActiveParticles = TimeBinsHydro.NActiveParticles;
  for(idx = 0; idx < NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;


#ifdef REFINEMENT_SPLIT_CELLS
      FlagDoNotRefine[i] = 0;
#endif
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip cells that have been swallowed or dissolved */

      /* check if this cell already has a partner to be merged with */
      if(SphP[i].DerefPartnerId > 0)
        {
          /* check partner index */
          if(SphP[i].DerefPartnerIndex >= NumGas || SphP[i].DerefPartnerIndex < 0 || P[SphP[i].DerefPartnerIndex].ID != SphP[i].DerefPartnerId)
            /* we have to find our partner again */
            voronoi_derefinement_update_partner_index(i);

          if(SphP[i].DerefPartnerIndex < 0)
            {
              /* our partner got lost, i.e. by a domain decomposition, so discard him */
              lcount++;
              SphP[i].DerefPartnerId = 0;
            }
        }

      if(SphP[i].DerefPartnerId == 0 && derefine_should_this_cell_be_merged(i, 0))
        {
          /* find a partner to merge with, this also sets the index already */
          voronoi_derefinement_find_partner(i);
          if(SphP[i].DerefPartnerId == 0)
            fcount++;
        }

      /* check if the pair is close enough to dissolve one of the cells */
      if(SphP[i].DerefPartnerId > 0 && SphP[i].DerefPartnerIndex >= 0)
        {
          FlagDoNotRefine[i] = 1;

          int j = SphP[i].DerefPartnerIndex;
          double dx = P[i].Pos[0] - P[j].Pos[0];
          double dy = P[i].Pos[1] - P[j].Pos[1];
          double dz = P[i].Pos[2] - P[j].Pos[2];

          double dist = sqrt(dx * dx + dy * dy + dz * dz);
          printf("Task: %d, ID=%d, PID=%d, dist=%g, thres=%g\n", ThisTask, P[i].ID, SphP[i].DerefPartnerId, dist, 0.01 * (get_cell_radius(i) + get_cell_radius(j)));
          if(dist < 0.01 * (get_cell_radius(i) + get_cell_radius(j)))
            {
              voronoi_derefinement_merge_pair(i, j);
              count++;
            }
          else
            pcount++;
        }
    }

  MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&pcount, &pcountall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&fcount, &fcountall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&lcount, &lcountall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  mpi_printf("DEREFINE: Number of cells we de-refined: %d, Number of pairs left: %d (lost: %d), failed to find partner for %d particles.\n", countall, (pcountall + 1) / 2, lcountall, fcountall);

  if(count > 0)
    {
      for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
        {
          i = TimeBinsGravity.ActiveParticleList[idx];
          if(i < 0)
            continue;

          if(P[i].ID == 0 && P[i].Mass == 0)
            timebin_remove_particle(&TimeBinsGravity, idx, P[i].TimeBinGrav);
        }
    }

  return 0;
}

void voronoi_derefinement_update_partner_index(int i)
{
  SphP[i].DerefPartnerIndex = -1;
  int q = SphP[i].first_connection;
  while(q >= 0)
    {
      int particle = DC[q].index;

      /* discard cells on other tasks and mirrored cells */
      if(DC[q].task != ThisTask || DC[q].image_flags != 1 || particle >= NumGas || particle < 0)
        {
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
          continue;
        }

      if(P[particle].ID == SphP[i].DerefPartnerId)
        {
          /* we found our partner */
          SphP[i].DerefPartnerIndex = particle;
          break;
        }

      if(q == SphP[i].last_connection)
        break;

      q = DC[q].next;
    }
}

void voronoi_derefinement_find_partner(int i)
{
  /* we look for the neighbour with the smallest mass,
     but only take into account cells which are on our
     task and do not already have a partner */
  double minmass = MAX_REAL_NUMBER;
  int partnerindex = -1;
  SphP[i].DerefPartnerIndex = -1;

  int q = SphP[i].first_connection;
  while(q >= 0)
    {
      int particle = DC[q].index;

      /* discard cells on other tasks and mirrored cells */
      if(DC[q].task != ThisTask || DC[q].image_flags != 1 || particle >= NumGas || particle < 0 || SphP[particle].DerefPartnerId != 0 || !TimeBinSynchronized[P[particle].TimeBinHydro])
        {
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
          continue;
        }

      if(P[particle].Mass < minmass)
        {
          partnerindex = particle;
          minmass = P[particle].Mass;
        }

      if(q == SphP[i].last_connection)
        break;

      q = DC[q].next;
    }

  if(partnerindex >= 0)
    {
      int j = partnerindex;
      /* we found a partner, lets register with him */
      SphP[i].DerefPartnerId = P[j].ID;
      SphP[j].DerefPartnerId = P[i].ID;

      SphP[i].DerefPartnerIndex = j;
      SphP[j].DerefPartnerIndex = i;
    }
}

void voronoi_derefinement_merge_pair(int i, int j)
{
  int keep, del;
  if(P[i].Mass > P[j].Mass)
    {
      keep = i;
      del = j;
    }
  else
    {
      keep = j;
      del = i;
    }

  P[keep].Mass += P[del].Mass;
  SphP[keep].Momentum[0] += SphP[del].Momentum[0];
  SphP[keep].Momentum[1] += SphP[del].Momentum[1];
  SphP[keep].Momentum[2] += SphP[del].Momentum[2];

#ifdef MHD
  /* it does not make any sense to "conserve" the magnetic field, so lets take the average */
  SphP[keep].BConserved[0] = 0.5 * (SphP[keep].BConserved[0] + SphP[del].BConserved[0]);
  SphP[keep].BConserved[1] = 0.5 * (SphP[keep].BConserved[1] + SphP[del].BConserved[1]);
  SphP[keep].BConserved[2] = 0.5 * (SphP[keep].BConserved[2] + SphP[del].BConserved[2]);
#endif

#ifdef MRT
  for(int num1=0; num1<MRT_BINS; num1++)
    {
      SphP[keep].Cons_DensPhot[num1] += SphP[del].Cons_DensPhot[num1] ;
      SphP[keep].Cons_RT_F[num1][0] += SphP[del].Cons_RT_F[num1][0] ;
      SphP[keep].Cons_RT_F[num1][1] += SphP[del].Cons_RT_F[num1][1] ;
      SphP[keep].Cons_RT_F[num1][2] += SphP[del].Cons_RT_F[num1][2] ;
    }
#endif

#ifndef ISOTHERM_EQS
  SphP[keep].Energy += SphP[del].Energy;
#endif
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[keep].Entropy += SphP[del].Entropy;
#endif
#ifdef MAXSCALARS
  int s;
  for(s = 0; s < N_Scalar; s++)
    *(MyFloat *) (((char *) (&SphP[keep])) + scalar_elements[s].offset_mass) += *(MyFloat *) (((char *) (&SphP[del])) + scalar_elements[s].offset_mass);
#endif
#ifdef TRACER_FIELD
  SphP[keep].ConservedTracer += SphP[del].ConservedTracer;
#endif
#if defined(FM_STAR_FEEDBACK) && defined(OUTPUT_STELLAR_FEEDBACK)
  SphP[keep].TotEgyFeed += SphP[del].TotEgyFeed;
  SphP[keep].IntEgyFeed += SphP[del].IntEgyFeed;
  SphP[keep].KinEgyFeed += SphP[del].KinEgyFeed;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
  voronoi_remove_connection(del);
#endif

  /* mark cell as deleted */
  P[del].ID = 0;
  P[del].Mass = 0;
  P[del].Vel[0] = 0;
  P[del].Vel[1] = 0;
  P[del].Vel[2] = 0;

  SphP[del].DerefPartnerId = 0;
  SphP[del].DerefPartnerIndex = -1;

  SphP[keep].DerefPartnerId = 0;
  SphP[keep].DerefPartnerIndex = -1;
}

void voronoi_derefinement_pairs_velvertex_corrections()
{
  int idx, i;
  double atime, hubble_a;
  if(All.ComovingIntegrationOn)
    {
      hubble_a = hubble_function(All.Time);
      atime = All.Time;
    }
  else
    atime = hubble_a = 1.0;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(SphP[i].DerefPartnerId > 0)
        {
          voronoi_derefinement_update_partner_index(i);
          if(SphP[i].DerefPartnerIndex >= 0)
            {
              int j = SphP[i].DerefPartnerIndex;

              /* if both partners are active (that should be the case), the one with the smaller ID
                 calculates both velocities */
              if(TimeBinSynchronized[P[j].TimeBinHydro] && P[j].ID < P[i].ID)
                continue;

              double dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / hubble_a;

              double vx = 0.5 * (SphP[i].VelVertex[0] + SphP[j].VelVertex[0]);
              double vy = 0.5 * (SphP[i].VelVertex[1] + SphP[j].VelVertex[1]);
              double vz = 0.5 * (SphP[i].VelVertex[2] + SphP[j].VelVertex[2]);

              double dvx = 0.3 * (P[j].Pos[0] - P[i].Pos[0]) / dt * atime * atime;
              double dvy = 0.3 * (P[j].Pos[1] - P[i].Pos[1]) / dt * atime * atime;
              double dvz = 0.3 * (P[j].Pos[2] - P[i].Pos[2]) / dt * atime * atime;

              SphP[i].VelVertex[0] = vx + dvx;
              SphP[i].VelVertex[1] = vy + dvy;
              SphP[i].VelVertex[2] = vz + dvz;

              if(TimeBinSynchronized[P[j].TimeBinHydro])
                {
                  SphP[j].VelVertex[0] = vx - dvx;
                  SphP[j].VelVertex[1] = vy - dvy;
                  SphP[j].VelVertex[2] = vz - dvz;
                }
            }
        }
    }
}
#endif

#endif /* VORONOI */
