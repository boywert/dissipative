/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/helper.c
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

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


void check_AuxDataID_references(void)
{
  mpi_printf("GFM: GFM Checks...\n");

  int i;
  for(i = 0; i < NumPart; i++)
    {
#ifdef GFM
      if(P[i].Type == 4)
        if(StarP[P[i].AuxDataID].PID != i)
          {
            printf("StarP broken: %llu %d %d %d %d %d\n", (long long) P[i].AuxDataID, i, StarP[P[i].AuxDataID].PID, NumGas, N_star, NumPart);
            terminate("StarP[P[i].AuxDataID].PID!=i\n");
          }
#endif

#ifdef BLACK_HOLES
      if(P[i].Type == 5)
        if(BHP[P[i].AuxDataID].PID != i)
          {
            printf("BHP broken: %llu %d %d %d %d %d\n", (long long) P[i].AuxDataID, i, BHP[P[i].AuxDataID].PID, NumGas, NumBHs, NumPart);
            terminate("BHP[P[i].AuxDataID].PID!=i\n");
          }
#endif
    }

  mpi_printf("GFM: done.\n");
}

#ifdef GFM
void gfm_inject_into_cell(int j, double inj_mass, double inj_thermalenergy, double inj_mom[3])
{

  if(inj_mass <= 0)
    terminate("GFM: inj_mass<=0");

  /* add injected momentum to cell */
  SphP[j].Momentum[0] += inj_mom[0];
  SphP[j].Momentum[1] += inj_mom[1];
  SphP[j].Momentum[2] += inj_mom[2];

  /* add injected mass to cell */
  P[j].Mass += inj_mass;

  /* add injected thermal energy to cell */
  SphP[j].Energy += inj_thermalenergy;

  /* add injected kinetic energy to cell */
  SphP[j].Energy += 0.5 * (inj_mom[0] * inj_mom[0] + inj_mom[1] * inj_mom[1] + inj_mom[2] * inj_mom[2]) / inj_mass;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  SphP[j].Utherm = (SphP[j].Energy -
                    0.5 * (SphP[j].Momentum[0] * SphP[j].Momentum[0] + SphP[j].Momentum[1] * SphP[j].Momentum[1] +
                           SphP[j].Momentum[2] * SphP[j].Momentum[2]) / P[j].Mass) / P[j].Mass / (All.cf_atime * All.cf_atime);
  SphP[j].A = (GAMMA - 1.0) * SphP[j].Utherm / pow(SphP[j].Density * All.cf_a3inv, GAMMA - 1);
  SphP[j].Entropy = log(SphP[j].A) * P[j].Mass;
#endif
}
#endif


#ifdef GFM
/* add a star particle to StarP */
void gfm_add_star(int i, int j, MyDouble mass_of_star, MyFloat birthtime, MyFloat hsml_guess)
{
  if(N_star >= All.MaxPartStar)
    terminate("There is no space left to creat new stars. N_star = %d, MaxPartStar = %d", N_star, All.MaxPartStar);

#ifdef GFM_STELLAR_EVOLUTION
  int k_elem;
#endif
  P[i].AuxDataID = N_star;

  /* zero StarP[] entries */
  memset(&StarP[N_star], 0, sizeof(struct star_particle_data));

  /* set values */
  StarP[N_star].PID = i;
  StarP[N_star].BirthTime = birthtime;

  StarP[N_star].BirthPos[0] = P[i].Pos[0];
  StarP[N_star].BirthPos[1] = P[i].Pos[1];
  StarP[N_star].BirthPos[2] = P[i].Pos[2];
  StarP[N_star].BirthVel[0] = P[i].Vel[0];
  StarP[N_star].BirthVel[1] = P[i].Vel[1];
  StarP[N_star].BirthVel[2] = P[i].Vel[2];
  StarP[N_star].BirthDensity = SphP[j].Density;

#ifdef GFM_STELLAR_EVOLUTION
  StarP[N_star].InitialMass = mass_of_star;
  StarP[N_star].Hsml = hsml_guess;

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    StarP[N_star].MassMetals[k_elem] = (mass_of_star / P[j].Mass) * SphP[j].MassMetals[k_elem];

  double metmass = 0;
  for(k_elem = 2; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    metmass += StarP[N_star].MassMetals[k_elem];

  StarP[N_star].Metallicity = metmass / mass_of_star;

#ifdef GFM_CHEMTAGS
  for(k_elem = 0; k_elem < GFM_N_CHEM_TAGS; k_elem++)
    StarP[N_star].MassMetalsChemTags[k_elem] = (mass_of_star / P[j].Mass) * SphP[j].MassMetalsChemTags[k_elem];
#endif

#ifdef GFM_DUST
  /* Assume that a star particle also starts with initial metals */
  /* due to dust in the ISM. */
  int chan;
  for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
    {
      for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
        {
          StarP[N_star].MassMetals[k_elem] += mass_of_star * SphP[j].MetalsDustFraction[chan][k_elem];
          StarP[N_star].Metallicity += SphP[j].MetalsDustFraction[chan][k_elem];
        }
    }

  for(k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      double tot_mass = 0.0;
      tot_mass += SphP[j].MassMetals[k_elem];
      for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          tot_mass += SphP[j].MassMetalsDust[chan][k_elem];
        }

      /* Normalize the star's initial metal and dust fractions.  If not */
      /* normalizable, just set to zero: the  ISM had no gas-phase or dust */
      /* metals. */
      if(tot_mass > 0.0)
        {
          StarP[N_star].InitialMetalFractions[k_elem] = SphP[j].MassMetals[k_elem] / tot_mass;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              StarP[N_star].InitialDustFractions[chan][k_elem] = SphP[j].MassMetalsDust[chan][k_elem] / tot_mass;
            }
        }
      else
        {
          StarP[N_star].InitialMetalFractions[k_elem] = 0.0;
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              StarP[N_star].InitialDustFractions[chan][k_elem] = 0.0;
            }
        }
    }
#endif

#ifdef GFM_DISCRETE_ENRICHMENT
  StarP[N_star].lastEnrichTime = birthtime;
#endif

#if defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF==0)
  StarP[N_star].DMVelDisp = SphP[j].w.DMVelDisp;
#endif
#endif

#ifdef INSTANTANEOUS_DEPOSITION
  StarP[N_star].FeedbackFlag = All.FeedbackInjectionEvents;     /* star is active for FeedbackInjectionEvents when it is born */
#endif

#ifdef FM_RADIATION_FEEDBACK
  StarP[N_star].RadFeed_Flag = 0;       /* it changes to -1 once the momentum due to radiation is input */
  StarP[N_star].StromgrenRadius = 0;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  StarP[N_star].StromgrenMass = 0;
#endif
  StarP[N_star].RadFeedTau = 0;
  StarP[N_star].RadFeed_NumNgb = 0;
#ifdef FM_RADIATION_FEEDBACK_DEBUG      /* #LVS: TODO --> once tested, keep variables in "debug" for only StarParticle */
  StarP[N_star].RadiationMomentumReleased = 0;
  StarP[N_star].NormSphRadFeedback = 0;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
  StarP[N_star].NormSphRadFeedback_cold = 0;
#endif
  StarP[N_star].RadCoolShutoffTime = 0;
#endif
#endif

#ifdef FM_EARLY_STAR_FEEDBACK
  StarP[N_star].EarlyCumulativeFeedback = 0.;
#endif

  N_star++;
}


#ifdef GFM_DUST_CAP
void gfm_cap_dust(void)
{
  for (int i = 0; i < NumGas; i++)
    {
      double sum = 0.0;
      for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              if (SphP[i].MetalsDustFraction[chan][k_elem] < 0)
                {
                  SphP[i].MetalsDustFraction[chan][k_elem] = 0.0;
                }
              if (SphP[i].MassMetalsDust[chan][k_elem] < 0)
                {
                  sum += -SphP[i].MassMetalsDust[chan][k_elem];
                  SphP[i].MassMetalsDust[chan][k_elem] = 0.0; 
                }
            }
        }
       SphP[i].DustMassCap += sum;
    }
}

void gfm_check_dust(int checkid)
{
  for (int i = 0; i < NumGas; i++)
    {
      for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
        {
          for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
            {
              if (SphP[i].MetalsDustFraction[chan][k_elem] < 0)
                terminate("GFM_DUST: %d negative dust %g\n", checkid, SphP[i].MetalsDustFraction[chan][k_elem]);
              if (SphP[i].MassMetalsDust[chan][k_elem] < 0)
                terminate("GFM_DUST: %d negative dust %g\n", checkid, SphP[i].MassMetalsDust[chan][k_elem]);
            }
        }
    }
}

#endif

#endif // GFM
