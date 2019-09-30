/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_production.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../proto.h"
#include "../allvars.h"

#ifdef DUST_LIVE
#ifdef DL_PRODUCTION

#ifndef DL_GRAIN_BINS
#error "Creation of dust particles also requires grain size evolution from DL_GRAIN_BINS!"
#endif

#if defined(DL_REFINEMENT) && (defined(ONEDIMS) || defined(TWODIMS))
#error "Dust refinement currently only designed for three dimensions!"
#endif

static double cum_mass_dust = 0.0;      /**< cumulative mass of dust created */
static double cum_mass_dust_exp = 0.0;  /**< cumulative mass of expected to be created */

void create_dust_particles(void)
{
  TIMER_START(CPU_DUST_PRODUCTION);

  int dust_spawned = 0;
  double sum_mass_dust = 0.0, sum_mass_dust_exp = 0.0;

  for(int idx = 0; idx < Nstar; idx++)
    {
      int i = StarParticle[idx].index;

      /* m_dust is the total mass of dust desired to be spawned, which may
       * be subdivided into multiple dust particles */
      double m_dust_each = All.DustTargetFrac * STP(i).InitialMass;
      double m_dust = All.NumDustPerSpawn * m_dust_each;
      double m_frac = STP(i).DeltaDustMassTot / P[i].Mass;
      double prob = P[i].Mass / m_dust * (1.0 - exp(-m_frac));
      sum_mass_dust_exp += STP(i).DeltaDustMassTot;

      if(prob == 0)
        continue;

      if(prob < 0)
        terminate("prob < 0, P[i].Mass=%g m_dust=%g m_frac=%g STP(i).DeltaDustMassTot=%g", P[i].Mass, m_dust, m_frac, STP(i).DeltaDustMassTot);

      if(prob > 1)
        {
          printf
            ("DUST_LIVE: Warning, need to make a dust particle heavier than desired. Task=%d prob=%g P[i].Mass=%g m_dust=%g m_dust_new=%g\n",
             ThisTask, prob, P[i].Mass, m_dust, P[i].Mass * (1.0 - exp(-m_frac)));
          m_dust = P[i].Mass * (1.0 - exp(-m_frac));
          prob = 1.0;
        }

      double p_decide = get_random_number();

      if(p_decide < prob)
        {
          /* Allow for spawning multiple dust particles. */
          for(int j = 0; j < All.NumDustPerSpawn; j++)
            {
              if(NumPart + dust_spawned >= All.MaxPart)
                terminate("For NumPart=%d and All.MaxPart=%d, no space left after spawning %d dust particles", NumPart, All.MaxPart, dust_spawned);

              int idust = NumPart + dust_spawned;

              spawn_dust_from_star(i, idust, m_dust, m_dust_each, j);

              dust_spawned++;
              sum_mass_dust += m_dust_each;
            }
        }

    } /* End of main loop over active star particles */

  int tot_dust_spawned;
  MPI_Allreduce(&dust_spawned, &tot_dust_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  All.TotN_dust += tot_dust_spawned;

  if(tot_dust_spawned > 0)
    mpi_printf("DUST_LIVE: spawned %d dust particles\n", tot_dust_spawned);

  if(tot_dust_spawned)
    {
      /* need to assign new unique IDs to the spawned dust */
      int *list;

      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&dust_spawned, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(int i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(int i = 0; i < dust_spawned; i++)
        {
          P[NumPart + i].ID = newid;
          newid++;
        }

      All.MaxID += tot_dust_spawned;
    }

  if(tot_dust_spawned > 0)
    {
      All.TotNumPart += tot_dust_spawned;
      NumPart += dust_spawned;
    }

  double din[2] = {sum_mass_dust, sum_mass_dust_exp}, dout[2];
  MPI_Reduce(din, dout, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      double tot_sum_mass_dust = dout[0];
      double tot_sum_mass_dust_exp = dout[1];

      cum_mass_dust += tot_sum_mass_dust;
      cum_mass_dust_exp += tot_sum_mass_dust_exp;

      fprintf(FdDust, "%14e %14e %14e\n", All.Time, cum_mass_dust, cum_mass_dust_exp);
      myflush(FdDust);
    }

  TIMER_STOP(CPU_DUST_PRODUCTION);
}

/* Take a grain size distribution of total mass m_dust = N*m_dust_each, and
 * modify it so the grain size distribution ignores the first j*m_dust_each of
 * mass, then contains m_dust_each of mass, then ignores the last
 * (N-j-1)*m_dust_each of mass.  This selects a contiguous section of the grain
 * size distribution.
 *
 * Summing the grain size distributions obtained in this way for all values of
 * j gives the original grain size distribution. */
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
void restrict_gsd(double m_dust_each, int j, double *num_grains, double *bin_slopes)
#else
void restrict_gsd(double m_dust_each, int j, double *num_grains)
#endif
{
  double bin_mass[DL_GRAIN_BINS];
  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_mass[i] = num_grains[i] * bin_avg_mass(i, num_grains[i], bin_slopes[i]);
#else
      bin_mass[i] = num_grains[i] * GSD.AvgMasses[i];
#endif
    }

  /* Entry i is the dust mass in bins 0, 1, ..., i-1, that is, up to the
   * left edge of bin i. */
  double cumulative_mass[DL_GRAIN_BINS+1];
  cumulative_mass[0] = 0.0;
  for(int i = 1; i <= DL_GRAIN_BINS; i++)
    {
      cumulative_mass[i] = cumulative_mass[i-1] + bin_mass[i-1];
    }

  for(int i = 0; i < DL_GRAIN_BINS; i++)
    {
      /* We're only interested in the bins containing cumulative mass between
       * j*m_dust_each and (j+1)*m_dust_each.  If a bin overlaps partially, use
       * the fractional overlap. */
      double overlap = interval_overlap(j*m_dust_each, (j+1)*m_dust_each, cumulative_mass[i], cumulative_mass[i+1]);
      double frac_overlap = overlap / bin_mass[i];

      /* Mass in a bin is a linear function of number of grains and slope, so
       * we can just scale down. */
      num_grains[i] *= frac_overlap;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_slopes[i] *= frac_overlap;
#endif
    }
}

/* Spawn a dust particle with index idust and mass m_dust_each from a star
 * particle with index istar.  The star particle may create multiple dust
 * particles during this timestep and subdivide the grain size distribution: in
 * that case, m_dust is the total mass of dust released (a multiple of
 * m_dust_each) and ispawn >= 0 denotes the ordering of the simultaneously
 * spawned dust particles. */
void spawn_dust_from_star(int istar, int idust, double m_dust, double m_dust_each, int ispawn)
{
  P[idust] = P[istar];
  P[idust].Type = DUST_LIVE;
  P[idust].SofteningType = All.SofteningTypeOfPartType[P[idust].Type];
  P[idust].Mass = m_dust_each;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[idust].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    P[idust].SofteningType = get_softening_type_from_mass(P[idust].Mass);
#endif

  // TODO: make a note of this factor in the startup
  /* Drop the dust particle a few timebins in case the grain size evolution is
   * more restrictive. */
  P[idust].TimeBinGrav -= 2;
  /* Initialize the location in TimeBinsDust using the gravity time bin.  The
   * dust time bin will be properly calculated at the start of the next
   * time-step. */
  P[idust].TimeBinHydro = P[idust].TimeBinGrav;
  timebin_add_particle(&TimeBinsGravity, idust, -1, P[idust].TimeBinGrav, TimeBinSynchronized[P[idust].TimeBinGrav]);
  timebin_add_particle(&TimeBinsDust, idust, -1, P[idust].TimeBinHydro, TimeBinSynchronized[P[idust].TimeBinHydro]);

  if(N_dust >= All.MaxPartDust)
    terminate("There is no space left to create new dust! N_dust = %d, MaxPartDust = %d", N_dust, All.MaxPartDust);

  /* Create new entry in DustP array, update indices as necessary */
  memset(&DustP[N_dust], 0, sizeof(struct dust_particle_data));
  P[idust].AuxDataID = N_dust;
  DustP[N_dust].PID = idust;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  DustP[N_dust].Hsml = All.SofteningTable[get_softening_type_from_mass(P[idust].Mass)];
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  DustP[N_dust].DustHsml = All.SofteningTable[get_softening_type_from_mass(P[idust].Mass)];
#endif
#else
  DustP[N_dust].Hsml = get_default_softening_of_particletype(DUST_LIVE);
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  DustP[N_dust].DustHsml = get_default_softening_of_particletype(DUST_LIVE);
#endif
#endif

  DustP[N_dust].BinMassChgTau = MAX_REAL_NUMBER;
  /* Needed for timestep criterion. */
  DustP[N_dust].OrigMass = P[idust].Mass;

  /* Initialize grain size distribution based on stellar type */
  double *num_grains;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  double *bin_slopes;
#endif

  if(STP(istar).DndaType == GSD_DNDA_AGB)
    {
      num_grains = GSD.AGB_NumGrains;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_slopes = GSD.AGB_BinSlopes;
#endif
    }
  else
    {
      num_grains = GSD.SNII_NumGrains;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      bin_slopes = GSD.SNII_BinSlopes;
#endif
    }

  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      /* Precomputed initial grain size distributions are normalized to unit
       * mass in internal units.  We scale them up by the total mass of dust
       * being spawned by the star, across possibly many dust particles. */
      DustP[N_dust].NumGrains[j] = m_dust * num_grains[j];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      DustP[N_dust].BinSlopes[j] = m_dust * bin_slopes[j];
#endif
    }

  /* Restrict grain size distribution so that it only contains the ispawn-th
   * contiguous chunk of the grain size distribution of mass m_dust_each. */
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  restrict_gsd(m_dust_each, ispawn, DustP[N_dust].NumGrains, DustP[N_dust].BinSlopes);
#else
  restrict_gsd(m_dust_each, ispawn, DustP[N_dust].NumGrains);
#endif

  /* Initialize element fractions based on yields from stars */
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      DustP[N_dust].MetalFractions[k] = STP(istar).DeltaDustMass[k] / STP(istar).DeltaDustMassTot;
    }

  /* During stellar evolution, we have not been reducing the star particle's
   * mass to account for dust production.  We need to do it here and keep
   * the star's metallicity the same. */
  double frac = 1.0 - m_dust_each / P[istar].Mass;
  if(frac <= 0.0)
    terminate("Spawning a dust particle would take away too much mass from a star particle! m_dust_each=%g, P[i].Mass=%g\n", m_dust_each, P[istar].Mass);
  P[istar].Mass *= frac;
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      STP(istar).MassMetals[k] *= frac;
    }

  N_dust++;
}

#ifdef DL_REFINEMENT
void refine_dust_particles(void)
{
  TIMER_START(CPU_DUST_PRODUCTION);

  int dust_created = 0;

  for(int idx = 0; idx < TimeBinsDust.NActiveParticles; idx++)
    {
      int i = TimeBinsDust.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(!((P[i].Type == DUST_LIVE) && (P[i].Mass > 0)))
        continue;

      /* Only refine sufficiently large dust particles. */
      if(P[i].Mass < All.DustMaxFrac * All.TargetGasMass)
        continue;

      if(NumPart + dust_created >= All.MaxPart)
        terminate("For NumPart=%d and All.MaxPart=%d, no space left after creating %d dust particles", NumPart, All.MaxPart, dust_created);

      int idust = NumPart + dust_created;

      //spawn_dust_from_star(i, idust, m_dust, m_dust_each, j);
      subdivide_dust_particle(i, idust);

      dust_created++;

    } /* End of main loop over active dust particles */

  int tot_dust_created;
  MPI_Allreduce(&dust_created, &tot_dust_created, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  All.TotN_dust += tot_dust_created;

  if(tot_dust_created > 0)
    mpi_printf("DUST_LIVE: refinement created %d dust particles\n", tot_dust_created);

  if(tot_dust_created)
    {
      /* need to assign new unique IDs to the created dust */
      int *list;

      if(All.MaxID == 0)        /* MaxID not calculated yet */
        calculate_maxid();

      list = mymalloc("list", NTask * sizeof(int));

      MPI_Allgather(&dust_created, 1, MPI_INT, list, 1, MPI_INT, MPI_COMM_WORLD);

      MyIDType newid = All.MaxID + 1;

      for(int i = 0; i < ThisTask; i++)
        newid += list[i];

      myfree(list);

      for(int i = 0; i < dust_created; i++)
        {
          P[NumPart + i].ID = newid;
          newid++;
        }

      All.MaxID += tot_dust_created;
    }

  if(tot_dust_created > 0)
    {
      All.TotNumPart += tot_dust_created;
      NumPart += dust_created;
    }

  TIMER_STOP(CPU_DUST_PRODUCTION);
}

/* Split off part of a dust particle with index iorig into a new dust particle
 * with index idust. */
void subdivide_dust_particle(int iorig, int idust)
{
  P[idust] = P[iorig];

  double oldmass = P[iorig].Mass;
  P[iorig].Mass = P[idust].Mass = oldmass/2.0;

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
  if(((1 << P[idust].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
    {
      P[iorig].SofteningType = get_softening_type_from_mass(P[iorig].Mass);
      P[idust].SofteningType = get_softening_type_from_mass(P[idust].Mass);
    }
#endif

  /* The new dust particle assumes the particle type and timebins of the
   * original dust particle. */
  timebin_add_particle(&TimeBinsGravity, idust, iorig, P[idust].TimeBinGrav, TimeBinSynchronized[P[idust].TimeBinGrav]);
  timebin_add_particle(&TimeBinsDust, idust, iorig, P[idust].TimeBinHydro, TimeBinSynchronized[P[idust].TimeBinHydro]);

  if(N_dust >= All.MaxPartDust)
    terminate("There is no space left to create new dust! N_dust = %d, MaxPartDust = %d", N_dust, All.MaxPartDust);

  /* Create new entry in DustP array, update indices as necessary */
  memset(&DustP[N_dust], 0, sizeof(struct dust_particle_data));
  P[idust].AuxDataID = N_dust;
  DustP[N_dust].PID = idust;

  /* The new dust particle also assumes the same smoothing lengths. */
  DustP[N_dust].Hsml = DTP(iorig).Hsml;
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
  DustP[N_dust].DustHsml = DTP(iorig).DustHsml;
#endif

  DustP[N_dust].BinMassChgTau = DTP(iorig).BinMassChgTau;
  /* Needed for timestep criterion. */
  DTP(iorig).OrigMass = P[iorig].Mass;
  DustP[N_dust].OrigMass = P[idust].Mass;

  /* Split the grain size distributions evenly. */
  for(int j = 0; j < DL_GRAIN_BINS; j++)
    {
      double oldnum = DTP(iorig).NumGrains[j];
      DTP(iorig).NumGrains[j] = DustP[N_dust].NumGrains[j] = oldnum/2.0;
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
      double oldslope = DTP(iorig).BinSlopes[j];
      DTP(iorig).BinSlopes[j] = DustP[N_dust].BinSlopes[j] = oldslope/2.0;
#endif
    }

  /* Offset the particles to preserve the center of mass. */
  double dir[3];
  double theta = acos(2 * get_random_number() - 1);
  double phi = 2 * M_PI * get_random_number();

  dir[0] = sin(theta) * cos(phi);
  dir[1] = sin(theta) * sin(phi);
  dir[2] = cos(theta);

  double fac = 0.025 * DTP(iorig).Hsml;

  P[iorig].Pos[0] += fac * dir[0];
  P[iorig].Pos[1] += fac * dir[1];
  P[iorig].Pos[2] += fac * dir[2];

  P[idust].Pos[0] -= fac * dir[0];
  P[idust].Pos[1] -= fac * dir[1];
  P[idust].Pos[2] -= fac * dir[2];

  /* Assume that the new dust particle and original dust particle share the
   * same metal fractions. */
  for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      DustP[N_dust].MetalFractions[k] = DTP(iorig).MetalFractions[k];
    }

  N_dust++;
}
#endif

#endif
#endif
