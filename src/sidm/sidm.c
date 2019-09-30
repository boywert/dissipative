/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm.c
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
#include "sidm_vars.h"


/* first init from init.c only called on startup not for restarts */
void sidm_Init_Particles(void)
{
  int i;
  unsigned char state;
  unsigned char reaction;
  double psum, xran;

  set_cosmo_factors_for_current_time();

  SIDM_arho = 1. / (All.cf_atime * All.cf_atime * All.cf_atime);
  SIDM_arho *= All.HubbleParam * All.HubbleParam;
  SIDM_avel = 1.0 / All.cf_atime;

  mpi_printf("SIDM: Init Particles...\n");
  if(All.SIDMDesNumNgb + All.SIDMMaxNumNgbDeviation > SIDM_MAX_NGBS)
    {
      terminate("All.SIDMDesNumNgb+All.SIDMMaxNumNgbDeviation>SIDM_MAX_NGBS\n");
    }
  if(All.SIDMDesNumNgb - All.SIDMMaxNumNgbDeviation < SIDM_MIN_NGBS)
    {
      terminate("All.SIDMDesNumNgb-All.SIDMMaxNumNgbDeviation<SIDM_MIN_NGBS\n");
    }

  /* set internal SIDM variables for each particle, we set all auxiliary SIDM variables here for all particles */
  for(i = 0; i < NumPart; i++)
    {
      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          P[i].sidm_PSum[reaction] = 0;
          P[i].sidm_NumTotalScatter[reaction] = 0;
        }
      P[i].sidm_Hsml = All.SofteningTable[P[i].SofteningType];
      for(state = 0; state < SIDM_STATES; state++)
        {
          P[i].sidm_Density[state] = 0.0;
          P[i].sidm_VelDisp[state] = 0.0;
        }
      P[i].sidm_NumNgb = 0;
    }

  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      All.sidm_ShouldScatter[reaction] = All.sidm_Scatters[reaction] = All.sidm_Rejected[reaction] = All.sidm_EnergyForbidden[reaction] = 0;
      All.sidm_EnergyInjected[reaction] = All.sidm_ScatteredMass[reaction] = 0.0;
    }

  All.sidm_EnergyInjectedCheck_sidm_parts = All.sidm_EnergyInjectedCheck_all_parts = 0.0;

  /* set up initial state population (i.e. set the state variable and change the mass accordingly) */
  for(i = 0; i < NumPart; i++)
    {
      /* here we have to filter by type to get the initial fractions and masses right */
      if ((1 << P[i].Type) & (SIDM))
      {
        xran = get_random_number();

        for (psum = 0.0, state = 0; state < SIDM_STATES; state++)
          {
            psum += STSIDM[state].InitialFraction;
            if (xran < psum)
              break;
          }

        if ((state >= SIDM_STATES) || (state < 0))
          terminate("SIDM: wrong state initialization state=%d  SIDM_STATES=%d\n", state, SIDM_STATES);

        P[i].sidm_State = state;

        /* this mass change is tiny; note that the check_omega routine is called before this mass change in init.c */ 
        P[i].Mass = (1.0 + STSIDM[state].DeltaMass) * All.SIDM_GroundStateMass;
      }
    }

  /* check the state stats */
  state_stats();

  mpi_printf("SIDM: done.\n");

}

void sidm_SetGroundStateMass()
{
  unsigned char reaction;
  double delta_E_1, delta_E_2, delta_E;
  double vmin;

  All.SIDM_clight = (2.99792458e10 / All.UnitVelocity_in_cm_per_s);
  All.SIDM_GroundStateMass = All.MassTable[1];                       //FIXME: particle mass
  
  mpi_printf("SIDM: scatter summary:\n");

  mpi_printf("SIDM: All.SIDM_GroundStateMass = %g\n", All.SIDM_GroundStateMass);
  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      sidm_get_delta_energy(reaction, &delta_E_1, &delta_E_2);
      delta_E = delta_E_1 + delta_E_2;
      vmin = sqrt(fabs(delta_E / (0.5 * 0.5 * All.SIDM_GroundStateMass))); //assuming mu=0.5 (equal mass particles) 
      vmin = (((delta_E)>(0)) ? (vmin) : (-vmin));
      mpi_printf("SIDM: Scatter Matrix: reaction=%d:  %d %d --> %d %d   delta_mass_1=%g   delta_mass_2=%g   delta_E_1=%g  delta_E_2=%g   delta_E=%g (required min. relative velocity=%g)\n", reaction, 
                  SMSIDM[reaction].In1, SMSIDM[reaction].In2, 
                  SMSIDM[reaction].Out1, SMSIDM[reaction].Out2, 
                  (STSIDM[SMSIDM[reaction].Out1].DeltaMass - STSIDM[SMSIDM[reaction].In1].DeltaMass), (STSIDM[SMSIDM[reaction].Out2].DeltaMass - STSIDM[SMSIDM[reaction].In2].DeltaMass), 
                  delta_E_1, delta_E_2, 
                  delta_E,
                  vmin);
    }
    mpi_printf("SIDM: init done.\n");
}

void sidm_AllocTargets(void)
{
  int idx, i;
  unsigned char reaction;

  Nforces = 0;
  TargetList = mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(((1 << P[i].Type) & (SIDM)) && Tree_Task_list[i] == ThisTask)
        {
          TargetList[Nforces++] = i;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if(((1 << (Tree_Points[i].Type)) & (SIDM)))
        {
          TargetList[Nforces++] = i + Tree_ImportedNodeOffset;
        }

  PSIDM = (struct sidm_data_p *) mymalloc("PSIDM", Nforces * sizeof(struct sidm_data_p));
  TSIDM = (struct sidm_data_t *) mymalloc("TSIDM", (NumPart + Tree_NumPartImported) * sizeof(struct sidm_data_t));
  LISIDM = (struct sidm_list *) mymalloc("LISIDM", Nforces * sizeof(struct sidm_list));


  memset(&PSIDM[0], 0, Nforces * sizeof(struct sidm_data_p));
  memset(&TSIDM[0], 0, (NumPart + Tree_NumPartImported) * sizeof(struct sidm_data_t));
  memset(&LISIDM[0], 0, Nforces * sizeof(struct sidm_list));

  for (reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    SIDM_Ekin_before_total[reaction] = SIDM_Ekin_after_total[reaction] = SIDM_Scattered_Mass[reaction] = 0.0;

  for(i = 0; i < NumPart; i++)
    {
      TSIDM[i].NewState  = P[i].sidm_State;
      TSIDM[i].NewMass   = P[i].Mass;
      TSIDM[i].NewVel[0] = P[i].Vel[0];
      TSIDM[i].NewVel[1] = P[i].Vel[1];
      TSIDM[i].NewVel[2] = P[i].Vel[2];
    }

  for(i = 0; i < Tree_NumPartImported; i++)
    {
      TSIDM[NumPart + i].NewState  = Tree_Points[i].sidm_State;
      TSIDM[NumPart + i].NewMass   = Tree_Points[i].Mass;
      TSIDM[NumPart + i].NewVel[0] = Tree_Points[i].Vel[0];
      TSIDM[NumPart + i].NewVel[1] = Tree_Points[i].Vel[1];
      TSIDM[NumPart + i].NewVel[2] = Tree_Points[i].Vel[2];
    }
}

void sidm_FreeTargets(void)
{
  myfree(LISIDM);
  myfree(TSIDM);
  myfree(PSIDM);
  myfree(TargetList);
}

/* main loop for scattering, called in run.c */
void sidm_DoScatter()
{
  mpi_printf("SIDM: starting...\n");

  TIMER_START(CPU_SIDM);

  //set_cosmo_factors_for_current_time();

  SIDM_arho = 1. / (All.cf_atime * All.cf_atime * All.cf_atime);
  SIDM_arho *= All.HubbleParam * All.HubbleParam;
  SIDM_avel = 1.0 / All.cf_atime;

  TIMER_STOPSTART(CPU_SIDM, CPU_SIDM_ALLOCFREE);

  sidm_AllocTargets();

  TIMER_STOPSTART(CPU_SIDM_ALLOCFREE, CPU_SIDM_HSML);

  sidm_findHsml();

  TIMER_STOPSTART(CPU_SIDM_HSML, CPU_SIDM_CHECK); 

  sidm_check_particle_scatter();

  TIMER_STOPSTART(CPU_SIDM_CHECK, CPU_SIDM_NGB);

  sidm_NgbList();

  TIMER_STOPSTART(CPU_SIDM_NGB, CPU_SIDM_ASSIGN);

  sidm_AssignScatterPartner();

  TIMER_STOPSTART(CPU_SIDM_ASSIGN, CPU_SIDM_SCATTER);

#ifndef SIDM_NO_SCATTER
  sidm_Scatter();
#endif

  TIMER_STOPSTART(CPU_SIDM_SCATTER, CPU_SIDM_STATS);

  scatter_stats();

  TIMER_STOPSTART(CPU_SIDM_STATS, CPU_SIDM_EXCHANGE);

  sidm_SetAndExchangeData();

  TIMER_STOPSTART(CPU_SIDM_EXCHANGE, CPU_SIDM_STATS);

  state_stats();

  TIMER_STOPSTART(CPU_SIDM_STATS, CPU_SIDM_ALLOCFREE);

  sidm_FreeTargets();

  TIMER_STOPSTART(CPU_SIDM_ALLOCFREE, CPU_SIDM);

  TIMER_STOP(CPU_SIDM);

  mpi_printf("SIDM: done.\n");
}

/* check which particle to scatter */
void sidm_check_particle_scatter(void)
{
  int idx, i;
  MyDouble dtime, xran, dt, time_hubble_a, hubble_a;
  MyDouble PSum[SIDM_REACTIONS], PSum_allreaction, PSum_partial;
  int done_flag;
  unsigned char reaction;

  SIDM_NumScatterParticles = 0;

  int nforces = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(((1 << P[i].Type) & (SIDM)) && Tree_Task_list[i] == ThisTask)
        {
          dt = (P[i].TimeBinGrav ? (((integertime) 1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;

          if(All.ComovingIntegrationOn)
            {
              hubble_a = hubble_function(All.Time);
              time_hubble_a = All.Time * hubble_a;
              dtime = All.Time * dt / time_hubble_a;
            }
          else
            {
              dtime = dt;
            }

          dtime *= 1.0 / All.HubbleParam;

          xran = get_random_number();
          for(reaction = 0, PSum_allreaction = 0.0; reaction < SIDM_REACTIONS; reaction++)
            {
              PSum[reaction] = PSIDM[nforces].PSum[reaction] * dtime;
              PSum_allreaction += PSum[reaction];
            }

          /* scatter in any reaction? */
          if(xran < PSum_allreaction)
            {
              PSIDM[nforces].RandX = xran / dtime;
              LISIDM[SIDM_NumScatterParticles].List1 = nforces;
              LISIDM[nforces].List2 = SIDM_NumScatterParticles;
              SIDM_NumScatterParticles++;
              /* scatter in which reaction? */
              for(reaction = 0, done_flag = 0, PSum_partial = 0.0; reaction < SIDM_REACTIONS; reaction++)
                {
                  PSum_partial += PSum[reaction];
                  if(PSum_partial > xran)
                    {
                      PSIDM[nforces].ScatterReaction = reaction;
                      PSIDM[nforces].ShouldScatterInStep[reaction] = 1;
                      done_flag = 1;
                      break;
                    }
                }
              if(done_flag == 0)
                terminate("no scatter reaction found PSum_allreaction=%g xran=%g PSum_partial=%g\n", PSum_allreaction, xran, PSum_partial);

            }

          nforces++;
        }
    }

  for(i = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(((1 << (Tree_Points[i].Type)) & (SIDM)))
      {
        dt = (Tree_Points[i].sidm_TimeBin ? (((integertime) 1) << Tree_Points[i].sidm_TimeBin) : 0) * All.Timebase_interval;

        if(All.ComovingIntegrationOn)
          {
            hubble_a = hubble_function(All.Time);
            time_hubble_a = All.Time * hubble_a;
            dtime = All.Time * dt / time_hubble_a;
          }
        else
          {
            dtime = dt;
          }

        dtime *= 1.0 / All.HubbleParam;

        xran = get_random_number();
        for(reaction = 0, PSum_allreaction = 0.0; reaction < SIDM_REACTIONS; reaction++)
          {
            PSum[reaction] = PSIDM[nforces].PSum[reaction] * dtime;
            PSum_allreaction += PSum[reaction];
          }

        /* scatter in any reaction? */
        if(xran < PSum_allreaction)
          {
            PSIDM[nforces].RandX = xran / dtime;
            LISIDM[SIDM_NumScatterParticles].List1 = nforces;
            LISIDM[nforces].List2 = SIDM_NumScatterParticles;
            SIDM_NumScatterParticles++;
            /* scatter in which reaction? */
            for(reaction = 0, done_flag = 0, PSum_partial = 0.0; reaction < SIDM_REACTIONS; reaction++)
              {
                PSum_partial += PSum[reaction];
                if(PSum_partial > xran)
                  {
                    PSIDM[nforces].ScatterReaction = reaction;
                    PSIDM[nforces].ShouldScatterInStep[reaction] = 1;
                    done_flag = 1;
                    break;
                  }
              }
            if(done_flag == 0)
              terminate("no scatter reaction found PSum_allreaction=%g xran=%g PSum_partial=%g\n", PSum_allreaction, xran, PSum_partial);
          }

        nforces++;
      }

}

/* ID compare function */
int ID_cmp(const void *a, const void *b)
{
  const MyIDType *ia = (const MyIDType *) a;
  const MyIDType *ib = (const MyIDType *) b;
  if(*ia < *ib)
    return -1;
  if(*ia > *ib)
    return +1;
  return 0;
}


/* ngb struct compare, increasing distance */
int ngb_entry_distance_cmp(const void *a, const void *b)
{
  const ngb_entry *ia = (const ngb_entry *) a;
  const ngb_entry *ib = (const ngb_entry *) b;
  if(ia->Distance < ib->Distance)
    return -1;
  if(ia->Distance > ib->Distance)
    return +1;
  return 0;
}


/* assign scatter partner to particles supposed to scatter */
void sidm_AssignScatterPartner(void)
{
  int n, p, done_flag, i, numdup;
  MyDouble PSum_partial;
  MyIDType *scatterlist_local, *scatterlist_global, *duplicate;
  int *count, *offset, tot_count;
  int TotNumScatterParticles = 2 * SIDM_NumScatterParticles;
  unsigned char scatter_reaction, reaction;

  scatterlist_local = mymalloc("scatterlist_local", TotNumScatterParticles * sizeof(MyIDType));


  /* select scatter partner */
  for(n = 0; n < SIDM_NumScatterParticles; n++)
    {
      int sidx = LISIDM[n].List1;

      if((PSIDM[n].ngb_Offset != PSIDM[sidx].NumNgb) || (PSIDM[sidx].NumNgb >= SIDM_MAX_NGBS) || (PSIDM[sidx].NumNgb <= SIDM_MIN_NGBS))
        {
          terminate("SIDM: offset=%d NumNgb=%d S1=%d S2=%d n=%d sidx=%d\n", PSIDM[n].ngb_Offset, PSIDM[sidx].NumNgb, LISIDM[n].List1, LISIDM[sidx].List2, n, sidx);
        }

      //sort by increasing distance
      qsort(&PSIDM[n].ngb_Entry[0], PSIDM[n].ngb_Offset, sizeof(ngb_entry), ngb_entry_distance_cmp);

      done_flag = 0;

      PSum_partial = 0.0;

      /* this is the reaction that is going to scatter */
      scatter_reaction = PSIDM[sidx].ScatterReaction;


      /* sum up everything before actuall scatter reaction */
      for(p = 0; p < PSIDM[n].ngb_Offset; p++)
        for(reaction = 0; reaction < scatter_reaction; reaction++)
          PSum_partial += PSIDM[n].ngb_Entry[p].P0j_half[reaction];


      /* now add the rest from scatter reaction until we cross the random value to select scatter partner */
      for(p = 0; p < PSIDM[n].ngb_Offset; p++)
        {
          PSum_partial += PSIDM[n].ngb_Entry[p].P0j_half[scatter_reaction];
          if(PSum_partial > PSIDM[sidx].RandX)
            {
              PSIDM[sidx].ScatterID = PSIDM[n].ngb_Entry[p].NgbIDs;
              scatterlist_local[n] = PSIDM[n].ngb_Entry[p].NgbIDs;
              scatterlist_local[n + SIDM_NumScatterParticles] = PSIDM[sidx].ID;
              done_flag = 1;
              break;
            }
        }

#ifdef SIDM_NO_NGB_SEL
      /* take closest particle to scatter with; will select particles in wrong state for multiple states; scatter state check is then turned off */
      PSIDM[sidx].ScatterID = PSIDM[n].ngb_Entry[1].NgbIDs;
      scatterlist_local[n] = PSIDM[n].ngb_Entry[1].NgbIDs;
      scatterlist_local[n + SIDM_NumScatterParticles] = PSIDM[sidx].ID;
#endif
      if(done_flag == 0)
        terminate("no scatter partner found offset=%d numbgb=%d scatter_reaction=%d PSum_partial=%g p=%d  PSIDM[sidx].RandX=%g\n", PSIDM[n].ngb_Offset, PSIDM[sidx].NumNgb, scatter_reaction, PSum_partial, p, PSIDM[sidx].RandX);
    }


  count = mymalloc("count", sizeof(int) * NTask);
  offset = mymalloc("offset", sizeof(int) * NTask);

  MPI_Allgather(&TotNumScatterParticles, 1, MPI_INT, count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, tot_count = 0, offset[0] = 0; i < NTask; i++)
    {
      tot_count += count[i];
      if(i > 0)
        offset[i] = offset[i - 1] + count[i - 1];
    }

  for(i = 0; i < NTask; i++)
    {
      count[i] *= sizeof(MyIDType);
      offset[i] *= sizeof(MyIDType);
    }

  scatterlist_global = mymalloc("scatterlist_global", tot_count * sizeof(MyIDType));
  MPI_Allgatherv(scatterlist_local, TotNumScatterParticles * sizeof(MyIDType), MPI_BYTE, scatterlist_global, count, offset, MPI_BYTE, MPI_COMM_WORLD);


  duplicate = mymalloc("duplicate", tot_count * sizeof(MyIDType));
  qsort(scatterlist_global, tot_count, sizeof(MyIDType), ID_cmp);

  numdup = 0;
  for(i = 1; i < tot_count; i++)
    {
      if(scatterlist_global[i - 1] == scatterlist_global[i])
        {
          duplicate[numdup++] = scatterlist_global[i];
        }
    }

  for(n = 0; n < SIDM_NumScatterParticles; n++)
    {
      int sidx = LISIDM[n].List1;
      unsigned char reaction = PSIDM[sidx].ScatterReaction;

      for(i = 0; i < numdup; i++)
        {
          if((PSIDM[sidx].ScatterID == duplicate[i]) || (PSIDM[sidx].ID == duplicate[i]))
            PSIDM[sidx].ShouldScatterInStep[reaction] = -1;
        }
    }

  myfree(duplicate);
  myfree(scatterlist_global);
  myfree(offset);
  myfree(count);
  myfree(scatterlist_local);
}

/* some state stats */
void state_stats(void)
{
  int i;
  unsigned char state;
  int input_ints[SIDM_STATES];
  long long output_longs[SIDM_STATES], tot_particles;
  double local_mass_states[SIDM_STATES], tot_mass_states[SIDM_STATES];
  

  for (state = 0; state < SIDM_STATES; state++)
    local_mass_states[state] = 0;

  for (i = 0; i < NumPart; i++)
    {
      if((1 << P[i].Type) & (SIDM)) 
        {
          local_mass_states[P[i].sidm_State] += P[i].Mass;
        }
    }

  MPI_Allreduce(&local_mass_states[0], &tot_mass_states[0], SIDM_STATES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(state = 0; state < SIDM_STATES; state++)
    input_ints[state] = 0;

  for(i = 0; i < NumPart; i++)
    if((1 << P[i].Type) & (SIDM))
      input_ints[P[i].sidm_State]++;

  sumup_large_ints(SIDM_STATES, input_ints, output_longs);

  for(tot_particles = 0, state = 0; state < SIDM_STATES; state++)
    tot_particles += output_longs[state];

  if(ThisTask == 0)
    for(state = 0; state < SIDM_STATES; state++)
      {
        printf("SIDM: %06.6f percent of particles (%012llu) in state %02d\n", (100.0 * output_longs[state] / tot_particles), output_longs[state], state);
        printf("SIDM: average particle mass in that state: (mass - groundstatemass)/groundstatemass: %g   (groundstatemass=%g)\n", (tot_mass_states[state]>0)?((tot_mass_states[state]/output_longs[state]-All.SIDM_GroundStateMass)/All.SIDM_GroundStateMass):0, All.SIDM_GroundStateMass);
      }


}



/* some scatter statistics */
void scatter_stats(void)
{
  int global_ShouldScatterInStep[SIDM_REACTIONS], local_ShouldScatterInStep[SIDM_REACTIONS];
  int global_ScattersInStep[SIDM_REACTIONS], local_ScattersInStep[SIDM_REACTIONS];
  int global_rejected[SIDM_REACTIONS], local_rejected[SIDM_REACTIONS];
  int global_EnergyForbidden[SIDM_REACTIONS], local_EnergyForbidden[SIDM_REACTIONS];
  double global_Ekin_before_total[SIDM_REACTIONS], global_Ekin_after_total[SIDM_REACTIONS], global_Scattered_Mass[SIDM_REACTIONS];
  int total_ShouldScatter = 0, total_Scatters = 0, total_Rejected = 0, total_EnergyForbidden = 0;
  double total_EnergyInjected = 0.0, total_ScatteredMass= 0.0;
  int i;
  unsigned char reaction;
  

  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      global_ShouldScatterInStep[reaction] = 0;
      local_ShouldScatterInStep[reaction] = 0;

      global_ScattersInStep[reaction] = 0;
      local_ScattersInStep[reaction] = 0;

      global_rejected[reaction] = 0;
      local_rejected[reaction] = 0;

      global_EnergyForbidden[reaction] = 0;
      local_EnergyForbidden[reaction] = 0;
    }

  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      for(i = 0; i < Nforces; i++)
        {
          local_EnergyForbidden[reaction] += PSIDM[i].EnergyForbidden[reaction] ;
        }
    }

  for(i = 0; i < Nforces; i++)
    {
      int target = TargetList[i];
      reaction = PSIDM[i].ScatterReaction;

      if(target > NumPart)
        target = NumPart + target - Tree_ImportedNodeOffset;

      if(PSIDM[i].ShouldScatterInStep[reaction] != 0)
        local_ShouldScatterInStep[reaction] += 1;

      if(PSIDM[i].ShouldScatterInStep[reaction] < 0)
        local_rejected[reaction] += 1;

#ifndef SIDM_NO_SCATTER
      /* note that the reverse check is not necessary since a scatter partner is usually not marked as a particle that is supposed to scatter in the first place */
      if((PSIDM[i].ShouldScatterInStep[reaction] > 0) && (TSIDM[target].ScattersInStep[reaction] == 0))
        {
          terminate("failed scatter: %d %d | local_index=%d global_index=%d Nforces=%d\n", TSIDM[target].ScattersInStep[reaction], PSIDM[i].ShouldScatterInStep[reaction], i, target, Nforces);
        }
#endif
    }

  for(i = 0; i < NumPart + Tree_NumPartImported; i++)
    {
      reaction = TSIDM[i].ScatterReaction;
      local_ScattersInStep[reaction] += (int) TSIDM[i].ScattersInStep[reaction];

      if(TSIDM[i].ScattersInStep[reaction] > 1)
        terminate("multiple scatter: %d %d/%d\n", TSIDM[i].ScattersInStep[reaction], i, NumPart);
    }

  for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
    {
      MPI_Allreduce(&local_ShouldScatterInStep[reaction], &global_ShouldScatterInStep[reaction], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_ScattersInStep[reaction], &global_ScattersInStep[reaction], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_rejected[reaction], &global_rejected[reaction], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_EnergyForbidden[reaction], &global_EnergyForbidden[reaction], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&SIDM_Ekin_before_total[reaction], &global_Ekin_before_total[reaction], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&SIDM_Ekin_after_total[reaction], &global_Ekin_after_total[reaction], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&SIDM_Scattered_Mass[reaction], &global_Scattered_Mass[reaction], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

 
  if(ThisTask == 0)
    {
      unsigned char reaction;
      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          printf("SIDM: STEP: SCATTERS:    reaction=%02d ShouldScatterInStep=%010d  ScattersInStep=%010d  rejected=%010d  (%06.2f percent)   energ_forbidden=%012d\n", reaction, 
                 global_ShouldScatterInStep[reaction], global_ScattersInStep[reaction], global_rejected[reaction],
                 100.0 * global_rejected[reaction] / (2.0 * global_ShouldScatterInStep[reaction] + 1e-5), global_EnergyForbidden[reaction]);
          printf("SIDM: STEP: KIN. ENERGY: reaction=%02d before scattering=%g  after scatter=%g absdel=%g reldelta=%g \n", reaction, global_Ekin_before_total[reaction], global_Ekin_after_total[reaction], global_Ekin_after_total[reaction]-global_Ekin_before_total[reaction], (global_Ekin_after_total[reaction]-global_Ekin_before_total[reaction])/(1e-20 + global_Ekin_before_total[reaction]));          
         }

      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        {
          All.sidm_ShouldScatter[reaction] += global_ShouldScatterInStep[reaction];
          All.sidm_Scatters[reaction] += global_ScattersInStep[reaction];
          All.sidm_Rejected[reaction] += global_rejected[reaction];
          All.sidm_EnergyForbidden[reaction] += global_EnergyForbidden[reaction];
          All.sidm_EnergyInjected[reaction] += (global_Ekin_after_total[reaction] - global_Ekin_before_total[reaction]);
          All.sidm_ScatteredMass[reaction] += global_Scattered_Mass[reaction]; 
          total_ShouldScatter += All.sidm_ShouldScatter[reaction];
          total_Scatters += All.sidm_Scatters[reaction];
          total_Rejected += All.sidm_Rejected[reaction];
          total_EnergyForbidden += All.sidm_EnergyForbidden[reaction];
	  total_EnergyInjected += All.sidm_EnergyInjected[reaction];
          total_ScatteredMass += All.sidm_ScatteredMass[reaction];
          printf("SIDM: TOTAL: SCATTERS:    reaction=%02d ShouldScatter=%010llu  Scatters=%010llu  Rejected=%010llu (%06.2f percent)  energ_forbidden=%012llu\n", reaction, 
                 All.sidm_ShouldScatter[reaction], All.sidm_Scatters[reaction], All.sidm_Rejected[reaction], 100.0 * All.sidm_Rejected[reaction] / (All.sidm_ShouldScatter[reaction] + 1e-20), All.sidm_EnergyForbidden[reaction]);

          printf("SIDM: TOTAL: KIN. ENERGY: reaction=%02d energy injected=%g  scattered mass=%g  velkick=%g\n", reaction, All.sidm_EnergyInjected[reaction], All.sidm_ScatteredMass[reaction], sqrt(fabs(All.sidm_EnergyInjected[reaction])/(1e-20 + All.sidm_ScatteredMass[reaction])));

        }
      printf("SIDM: total over all reactions\n");
      printf("SIDM: total: scatters sum over all reaction channels: ShouldScatter = %010d  Scatters = %010d  Rejected=%010d  (%06.2f percent)\n",
             total_ShouldScatter, total_Scatters, total_Rejected, 100.0 * total_Rejected / (total_ShouldScatter + 1e-20));
      printf("SIDM: total: injected energy over all reaction channels: energy injected=%g scattered mass=%g\n", total_EnergyInjected, total_ScatteredMass);
      fflush(stdout);
    }
}


static struct resultsimported_data
{
  int index;
  MyDouble sidm_Density[SIDM_STATES];
  MyDouble sidm_VelDisp[SIDM_STATES];
  MyDouble sidm_PSum[SIDM_REACTIONS];
  MyDouble sidm_NumNgb;
  MyDouble sidm_Hsml;
  double sidm_Mass;
  unsigned char sidm_State;
  int ScattersInStep[SIDM_REACTIONS];
  MyFloat Vel[3];
}
 *Tree_ResultsImported;


void sidm_SetAndExchangeData(void)
{
  int idx, i, j, k, nexport, nimport;
  int ngrp, recvTask, n;
  int ncount, ncount1, ncount2;
  struct resultsimported_data *tmp_results;
  unsigned char state;
  unsigned char reaction;

  /* first all active particles -> update all properties that are always calculated; i.e. those that also need updates even without a scattering happening */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      if(((1 << (Tree_Points[i].Type)) & (SIDM)))
        ncount++;

  Tree_ResultsImported = mymalloc("Tree_HsmlResultsImported", ncount * sizeof(struct resultsimported_data));

  memset(Tree_ResultsImported, 0, ncount * sizeof(struct resultsimported_data));

  for(idx = 0, ncount1 = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(((1 << P[i].Type) & (SIDM)) && Tree_Task_list[i] == ThisTask)
        {
          P[i].sidm_NumNgb = PSIDM[ncount1].NumNgb;
          P[i].sidm_Hsml = PSIDM[ncount1].Hsml;
          for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
            P[i].sidm_PSum[reaction] = PSIDM[ncount1].PSum[reaction];
          for(state = 0; state < SIDM_STATES; state++)
            {
              P[i].sidm_Density[state] = PSIDM[ncount1].Density[state];
              if(PSIDM[ncount1].NumNgbState[state] > 0)
                {
                  PSIDM[ncount1].Vx[state] /= PSIDM[ncount1].NumNgbState[state];
                  PSIDM[ncount1].Vy[state] /= PSIDM[ncount1].NumNgbState[state];
                  PSIDM[ncount1].Vz[state] /= PSIDM[ncount1].NumNgbState[state];
                  PSIDM[ncount1].VelDisp[state] /= PSIDM[ncount1].NumNgbState[state];
                  P[i].sidm_VelDisp[state] =
                    sqrt(PSIDM[ncount1].VelDisp[state] - PSIDM[ncount1].Vx[state] * PSIDM[ncount1].Vx[state] - PSIDM[ncount1].Vy[state] * PSIDM[ncount1].Vy[state] -
                         PSIDM[ncount1].Vz[state] * PSIDM[ncount1].Vz[state]);
                }
            }
          ncount1++;
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
    if(((1 << (Tree_Points[i].Type)) & (SIDM)))
      {
        Tree_ResultsImported[ncount2].sidm_NumNgb = PSIDM[ncount1].NumNgb;
        Tree_ResultsImported[ncount2].sidm_Hsml = PSIDM[ncount1].Hsml;
        for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
          Tree_ResultsImported[ncount2].sidm_PSum[reaction] = PSIDM[ncount1].PSum[reaction];
        for(state = 0; state < SIDM_STATES; state++)
          {
            Tree_ResultsImported[ncount2].sidm_Density[state] = PSIDM[ncount1].Density[state];
            if(PSIDM[ncount1].NumNgbState[state] > 0)
              {
                PSIDM[ncount1].Vx[state] /= PSIDM[ncount1].NumNgbState[state];
                PSIDM[ncount1].Vy[state] /= PSIDM[ncount1].NumNgbState[state];
                PSIDM[ncount1].Vz[state] /= PSIDM[ncount1].NumNgbState[state];
                PSIDM[ncount1].VelDisp[state] /= PSIDM[ncount1].NumNgbState[state];
                Tree_ResultsImported[ncount2].sidm_VelDisp[state] =
                  sqrt(PSIDM[ncount1].VelDisp[state] - PSIDM[ncount1].Vx[state] * PSIDM[ncount1].Vx[state] - PSIDM[ncount1].Vy[state] * PSIDM[ncount1].Vy[state] -
                       PSIDM[ncount1].Vz[state] * PSIDM[ncount1].Vz[state]);
              }
          }
        ncount1++;
        ncount2++;
      }


  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          if(((1 << (Tree_Points[n].Type)) & (SIDM)))
            {
              Tree_ResultsImported[k].index = Tree_Points[n].index;
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

  tmp_results = mymalloc("tmp_results", nexport * sizeof(struct resultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct resultsimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_ResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct resultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct resultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;
      for(state = 0; state < SIDM_STATES; state++)
        {
          P[target].sidm_Density[state] = tmp_results[i].sidm_Density[state];
          P[target].sidm_VelDisp[state] = tmp_results[i].sidm_VelDisp[state];
        }
      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        P[target].sidm_PSum[reaction] = tmp_results[i].sidm_PSum[reaction];
      P[target].sidm_NumNgb = tmp_results[i].sidm_NumNgb;
      P[target].sidm_Hsml = tmp_results[i].sidm_Hsml;
    }

  myfree(tmp_results);
  myfree(Tree_ResultsImported);


  /* now all scatter particles -> update scatter quantities */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
    {
      reaction = TSIDM[NumPart + i].ScatterReaction;
      if(TSIDM[NumPart + i].ScattersInStep[reaction])
        ncount++;
    }

  Tree_ResultsImported = mymalloc("Tree_HsmlResultsImported", ncount * sizeof(struct resultsimported_data));
  memset(Tree_ResultsImported, 0, ncount * sizeof(struct resultsimported_data));

  for(i = 0; i < NumPart; i++)
    {
      reaction = TSIDM[i].ScatterReaction;
      if(TSIDM[i].ScattersInStep[reaction])
        {
          P[i].sidm_NumTotalScatter[reaction] += TSIDM[i].ScattersInStep[reaction];
          P[i].sidm_State = TSIDM[i].NewState;
          P[i].Mass = TSIDM[i].NewMass;
#ifndef SIDM_NO_KINEMATICS
          P[i].Vel[0] = TSIDM[i].NewVel[0];
          P[i].Vel[1] = TSIDM[i].NewVel[1];
          P[i].Vel[2] = TSIDM[i].NewVel[2];
          if ((isnan(P[i].Vel[0])) || (isnan(P[i].Vel[1])) || (isnan(P[i].Vel[2])) || (isinf(P[i].Vel[0])) || (isinf(P[i].Vel[1])) || (isinf(P[i].Vel[2])))
            terminate("SIDM: BAD VELOCITIES AFTER SCATTER: (%g|%g|%g)  %d\n", P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], i);
#endif
        }
    }

  for(i = 0, ncount2 = 0; i < Tree_NumPartImported; i++)
    {
      reaction = TSIDM[NumPart + i].ScatterReaction;
      if(TSIDM[NumPart + i].ScattersInStep[reaction])
        {
          Tree_ResultsImported[ncount2].ScattersInStep[reaction] = TSIDM[NumPart + i].ScattersInStep[reaction];
          Tree_ResultsImported[ncount2].sidm_State = TSIDM[NumPart + i].NewState;
          Tree_ResultsImported[ncount2].sidm_Mass = TSIDM[NumPart + i].NewMass;
#ifndef SIDM_NO_KINEMATICS
          Tree_ResultsImported[ncount2].Vel[0] = TSIDM[NumPart + i].NewVel[0];
          Tree_ResultsImported[ncount2].Vel[1] = TSIDM[NumPart + i].NewVel[1];
          Tree_ResultsImported[ncount2].Vel[2] = TSIDM[NumPart + i].NewVel[2];
#endif
          ncount2++;
        }
    }
  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      {
        reaction = TSIDM[NumPart + n].ScatterReaction;
        if(TSIDM[NumPart + n].ScattersInStep[reaction])
          {
            Tree_ResultsImported[k].index = Tree_Points[n].index;
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

  tmp_results = mymalloc("tmp_results", nexport * sizeof(struct resultsimported_data));
  memset(tmp_results, -1, nexport * sizeof(struct resultsimported_data));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Tree_ResultsImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct resultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                       &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct resultsimported_data), MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  for(i = 0; i < nexport; i++)
    {
      int target = tmp_results[i].index;

      for(reaction = 0; reaction < SIDM_REACTIONS; reaction++)
        P[target].sidm_NumTotalScatter[reaction] += tmp_results[i].ScattersInStep[reaction];
      P[target].sidm_State = tmp_results[i].sidm_State;
      P[target].Mass = tmp_results[i].sidm_Mass;
#ifndef SIDM_NO_KINEMATICS
      P[target].Vel[0] = tmp_results[i].Vel[0];
      P[target].Vel[1] = tmp_results[i].Vel[1];
      P[target].Vel[2] = tmp_results[i].Vel[2];
          if ((isnan(P[target].Vel[0])) || (isnan(P[target].Vel[1])) || (isnan(P[target].Vel[2])) || (isinf(P[target].Vel[0])) || (isinf(P[target].Vel[1])) || (isinf(P[target].Vel[2])))
            terminate("SIDM: BAD VELOCITIES AFTER SCATTER: (%g|%g|%g)  %d\n", P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], target);

#endif
    }

  myfree(tmp_results);
  myfree(Tree_ResultsImported);

}


double sidm_calc_kinetic_energy_sidm_parts(void)
{
 int i;
 double local_kin_energy = 0.0, global_kin_energy = 0.0;

 SIDM_avel = 1.0 / All.cf_atime;

 for (i = 0; i < NumPart; i++)
   if ((1 << P[i].Type) & (SIDM))
     local_kin_energy +=  SIDM_avel * SIDM_avel * (0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]));

  MPI_Allreduce(&local_kin_energy, &global_kin_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 return global_kin_energy;
}

double sidm_calc_kinetic_energy_all_parts(void)
{
 int i;
 double local_kin_energy = 0.0, global_kin_energy = 0.0;

 SIDM_avel = 1.0 / All.cf_atime;

 for (i = 0; i < NumPart; i++)
   local_kin_energy +=  SIDM_avel * SIDM_avel * (0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]));

  MPI_Allreduce(&local_kin_energy, &global_kin_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 return global_kin_energy;
}

