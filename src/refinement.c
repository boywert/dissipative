/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/refinement.c
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

#include "allvars.h"

#ifdef REFINEMENT

#include "proto.h"

#if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS)
char *FlagDoNotRefine;
#endif

#ifdef VORONOI
static void refinement_prepare();
static void refinement_cleanup();

void do_derefinements_and_refinements()
{
  refinement_prepare();

#ifdef REFINEMENT_MERGE_CELLS
  do_derefinements();
#endif

#ifdef REFINEMENT_SPLIT_CELLS
  do_refinements();
#endif

  refinement_cleanup();
}

void refinement_prepare()
{
  TIMER_START(CPU_REFINE);

#ifdef REFINEMENT_VOLUME_LIMIT
  int idx, i;
#endif

#if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS)
  FlagDoNotRefine = mymalloc_movable(&FlagDoNotRefine, "FlagDoNotRefine", NumGas * sizeof(char));
#endif

#ifdef REFINEMENT_VOLUME_LIMIT
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      SphP[i].MinNgbVolume = MAX_REAL_NUMBER;

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            {
              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
              continue;
            }

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

          double Volume;
          if(DC[q].task == ThisTask)
            Volume = SphP[particle].Volume;
          else
            {
#ifndef OPTIMIZE_MESH_MEMORY_FOR_REFINEMENT
              Volume = PrimExch[particle].Volume;
#else
              Volume = RefExch[particle].Volume;
#endif
            }

          if(Volume < SphP[i].MinNgbVolume)
            SphP[i].MinNgbVolume = Volume;

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

    }
#endif

#ifdef REFINEMENT_AROUND_BH
  blackhole_mark_cells_for_refinement();
#endif

#ifdef REFINEMENT_AROUND_DM
  dm_particle_update_list();
#endif

  TIMER_STOP(CPU_REFINE);
}

void refinement_cleanup()
{
#if defined (REFINEMENT_MERGE_CELLS) && defined (REFINEMENT_SPLIT_CELLS)
  myfree(FlagDoNotRefine);
#endif
}

#endif

void move_collisionless_particle(int new_i, int old_i)
{
  int prev, next, bin;
  struct TimeBinData *tbData;

  P[new_i] = P[old_i];

#ifdef GFM
  if(P[new_i].Type == 4)
    StarP[P[new_i].AuxDataID].PID = new_i;
#endif
#ifdef BLACK_HOLES
  if(P[new_i].Type == 5)
    BHP[P[new_i].AuxDataID].PID = new_i;
#endif
#ifdef DUST_LIVE
  if(P[new_i].Type == DUST_LIVE)
    DustP[P[new_i].AuxDataID].PID = new_i;
#endif

  if(P[old_i].Mass == 0 && P[old_i].ID == 0)
    return;

  if(P[old_i].Mass == 0 && P[old_i].Type == 4)
    return;

#ifdef DUST_LIVE
  if(P[old_i].Mass == 0 && P[old_i].Type == DUST_LIVE)
    return;
#endif

  tbData = &TimeBinsGravity;
  bin = P[old_i].TimeBinGrav;

#ifdef TRACER_PARTICLE
  if(P[old_i].Type == TRACER_PARTICLE)
    {
      tbData = &TimeBinsTracer;
      bin = P[old_i].TimeBinHydro;
    }
#endif

  if(TimeBinSynchronized[bin])
    {
      /* particle is active, need to add it to the list of active particles again
         we assume here, that the new particle at the old index in this list is also active! */
      tbData->ActiveParticleList[tbData->NActiveParticles] = new_i;
      tbData->NActiveParticles++;
    }

  /* now move it in the link list of its timebin 
     we only need to change the gravity timebin here */

  tbData->NextInTimeBin[new_i] = tbData->NextInTimeBin[old_i];
  tbData->PrevInTimeBin[new_i] = tbData->PrevInTimeBin[old_i];

  prev = tbData->PrevInTimeBin[old_i];
  next = tbData->NextInTimeBin[old_i];

  if(prev >= 0)
    tbData->NextInTimeBin[prev] = new_i;
  else
    {
      if(tbData->FirstInTimeBin[bin] != old_i)
        terminate("strange");
      tbData->FirstInTimeBin[bin] = new_i;
    }

  if(next >= 0)
    tbData->PrevInTimeBin[next] = new_i;
  else
    {
      if(tbData->LastInTimeBin[bin] != old_i)
        terminate("strange");
      tbData->LastInTimeBin[bin] = new_i;
    }
    
#ifdef BLACK_HOLES
  /* move in BH accretion timebin */
  if(P[old_i].Type == 5)
  {
	  if(TimeBinSynchronized[bin])
	    {
	      /* particle is active, need to add it to the list of active particles again
	         we assume here, that the new particle at the old index in this list is also active! */
		  TimeBinsBHAccretion.ActiveParticleList[TimeBinsBHAccretion.NActiveParticles] = new_i;
		  TimeBinsBHAccretion.NActiveParticles++;
	    }

	  /* now move it in the link list of its timebin
	     we only need to change the gravity timebin here */

	  TimeBinsBHAccretion.NextInTimeBin[new_i] = TimeBinsBHAccretion.NextInTimeBin[old_i];
	  TimeBinsBHAccretion.PrevInTimeBin[new_i] = TimeBinsBHAccretion.PrevInTimeBin[old_i];

	  prev = TimeBinsBHAccretion.PrevInTimeBin[old_i];
	  next = TimeBinsBHAccretion.NextInTimeBin[old_i];

	  if(prev >= 0)
		  TimeBinsBHAccretion.NextInTimeBin[prev] = new_i;
	  else
	    {
	      if(TimeBinsBHAccretion.FirstInTimeBin[bin] != old_i)
	        terminate("strange");
	      TimeBinsBHAccretion.FirstInTimeBin[bin] = new_i;
	    }

	  if(next >= 0)
		  TimeBinsBHAccretion.PrevInTimeBin[next] = new_i;
	  else
	    {
	      if(TimeBinsBHAccretion.LastInTimeBin[bin] != old_i)
	        terminate("strange");
	      TimeBinsBHAccretion.LastInTimeBin[bin] = new_i;
	    }
  }
#endif

#ifdef SINKS
  /* move in Sinks accretion timebin */
  if(P[old_i].Type == 5)
  {
	  if(TimeBinSynchronized[bin])
	    {
	      /* particle is active, need to add it to the list of active particles again
	         we assume here, that the new particle at the old index in this list is also active! */
		  TimeBinsSinksAccretion.ActiveParticleList[TimeBinsSinksAccretion.NActiveParticles] = new_i;
		  TimeBinsSinksAccretion.NActiveParticles++;
	    }

	  /* now move it in the link list of its timebin
	     we only need to change the gravity timebin here */

	  TimeBinsSinksAccretion.NextInTimeBin[new_i] = TimeBinsSinksAccretion.NextInTimeBin[old_i];
	  TimeBinsSinksAccretion.PrevInTimeBin[new_i] = TimeBinsSinksAccretion.PrevInTimeBin[old_i];

	  prev = TimeBinsSinksAccretion.PrevInTimeBin[old_i];
	  next = TimeBinsSinksAccretion.NextInTimeBin[old_i];

	  if(prev >= 0)
		  TimeBinsSinksAccretion.NextInTimeBin[prev] = new_i;
	  else
	    {
	      if(TimeBinsSinksAccretion.FirstInTimeBin[bin] != old_i)
	        terminate("strange");
	      TimeBinsSinksAccretion.FirstInTimeBin[bin] = new_i;
	    }

	  if(next >= 0)
		  TimeBinsSinksAccretion.PrevInTimeBin[next] = new_i;
	  else
	    {
	      if(TimeBinsSinksAccretion.LastInTimeBin[bin] != old_i)
	        terminate("strange");
	      TimeBinsSinksAccretion.LastInTimeBin[bin] = new_i;
	    }
  }
#endif

#ifdef DUST_LIVE
  /* move in dust timebin */
  if(P[old_i].Type == DUST_LIVE)
    {
      if(TimeBinSynchronized[bin])
        {
          /* particle is active, need to add it to the list of active particles again
             we assume here, that the new particle at the old index in this list is also active! */
          TimeBinsDust.ActiveParticleList[TimeBinsDust.NActiveParticles] = new_i;
          TimeBinsDust.NActiveParticles++;
        }

      /* now move it in the link list of its timebin */
      TimeBinsDust.NextInTimeBin[new_i] = TimeBinsDust.NextInTimeBin[old_i];
      TimeBinsDust.PrevInTimeBin[new_i] = TimeBinsDust.PrevInTimeBin[old_i];

      prev = TimeBinsDust.PrevInTimeBin[old_i];
      next = TimeBinsDust.NextInTimeBin[old_i];

      if(prev >= 0)
        TimeBinsDust.NextInTimeBin[prev] = new_i;
      else
        {
          if(TimeBinsDust.FirstInTimeBin[bin] != old_i)
            terminate("strange");
          TimeBinsDust.FirstInTimeBin[bin] = new_i;
        }

      if(next >= 0)
        TimeBinsDust.PrevInTimeBin[next] = new_i;
      else
        {
          if(TimeBinsDust.LastInTimeBin[bin] != old_i)
            terminate("strange");
          TimeBinsDust.LastInTimeBin[bin] = new_i;
        }
    }
#endif
}

#endif /* REFINEMENT */
