/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/stellar_evolution_util.c
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

#ifdef GFM_STELLAR_EVOLUTION


void start_enrichment()
{
  int idx, i, k;
#ifdef GFM_DUST
  int chan;
#endif

  TIMER_START(CPU_GFM_ENRICH);

  //  mpi_printf("GFM_STELLAR_EVOLUTION: GFM_STELLAR_EVOLUTION...\n");

  StarParticle = mymalloc("StarParticle", N_star * sizeof(struct star_particle));

  Nstar = 0;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Ti_Current != All.Ti_Current)
        {
          terminate("how can this be?");
          drift_particle(i, All.Ti_Current);
        }

      if(P[i].Type == 4 && P[i].Mass > 0 && STP(i).BirthTime > 0)
        {
          StarParticle[Nstar].index = i;
          StarParticle[Nstar].NumNgb = 0;
          StarParticle[Nstar].NormSph = 0;
          StarParticle[Nstar].TotalMassReleased = -1;   /* indicates that we have not found its hsml value */
          StarParticle[Nstar].TotalMetalMassReleased = 0;
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            StarParticle[Nstar].MetalMassReleased[k] = 0;

#ifdef GFM_DUST
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                {
                  StarParticle[Nstar].DustMassReleased[chan][k] = 0.0;
                }
            }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
          StarParticle[Nstar].NumSNII = 0.0;
#endif
#ifdef GFM_CHEMTAGS
          for(k = 0; k < GFM_N_CHEM_TAGS; k++)
            StarParticle[Nstar].MetalMassReleasedChemTags[k] = 0;
#endif

#ifdef FM_STAR_FEEDBACK
          StarParticle[Nstar].TotalEnergyReleased = 0;
          StarParticle[Nstar].TotalMomentumInjected = 0;
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
          StarParticle[Nstar].TotalMomentumReleased = 0;
#endif
#ifdef DELAYED_COOLING
          StarParticle[Nstar].AvgPress = 0;
          StarParticle[Nstar].AvgHDens = 0;
          StarParticle[Nstar].BlastRadius = 0;
          StarParticle[Nstar].CoolShutoffTime = 0;
#else
#ifndef DELAYED_COOLING_TURB
          StarParticle[Nstar].deltaEKin = 0;
          StarParticle[Nstar].deltaMomentum[0] = 0;
          StarParticle[Nstar].deltaMomentum[1] = 0;
          StarParticle[Nstar].deltaMomentum[2] = 0;
#endif
#endif
#ifdef FM_RADIATION_FEEDBACK
          StarParticle[Nstar].RadiationMomentumReleased = 0;
          StarParticle[Nstar].StromgrenRadius = 0;
#ifdef FM_STOCHASTIC_HII_PHOTOIONIZATION
          StarParticle[Nstar].StromgrenMass = 0;
          StarParticle[Nstar].NormSphRadFeedback_cold = 0;
#endif
          StarParticle[Nstar].NormSphRadFeedback = 0;
          StarParticle[Nstar].RadCoolShutoffTime = 0;
          StarParticle[Nstar].RadFeedTau = 0;
          StarParticle[Nstar].RadFeed_MinGasDist = 0;
          StarParticle[Nstar].RadFeed_NumNgb = 0;
#endif
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
          StarParticle[Nstar].ClosestNeighbourID = 0;
#endif
#endif
#ifdef FM_EARLY_STAR_FEEDBACK
          StarParticle[Nstar].EarlyTotalEnergyReleased = 0;
          StarParticle[Nstar].deltaEarlyEKin = 0;
          StarParticle[Nstar].deltaEarlyMomentum[0] = 0;
          StarParticle[Nstar].deltaEarlyMomentum[1] = 0;
          StarParticle[Nstar].deltaEarlyMomentum[2] = 0;
          StarParticle[Nstar].EarlyMomentumInjected = 0;
#endif
#ifdef FM_VAR_SN_EFF
          StarParticle[Nstar].AvgMetalNgb = 0;   
#endif
#ifdef FM_MASS_WEIGHT_SN
          StarParticle[Nstar].TotNgbMass = 0;   
#endif
          Nstar++;
        }
    }

}

void end_enrichment()
{
  myfree(StarParticle);

  mpi_printf("GFM_STELLAR_EVOLUTION: done.\n");

  TIMER_STOP(CPU_GFM_ENRICH);
}

#ifdef GFM_WINDS_STRIPPING
void start_stripping()
{
  TIMER_START(CPU_GFM_ENRICH);

  mpi_printf("GFM_WINDS_STRIPPING: starting...\n");

  int idx, i, k;
#ifdef GFM_DUST
  int chan;
#endif


  float targetAge = (float) (-1.0 * All.WindFreeTravelMaxTimeFactor / All.cf_hubble_a);
  float currentAge = 0.0;

  Nwinds_to_strip = 0;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      currentAge = (float) STP(i).BirthTime;
      if((P[i].Type == 4 && P[i].Mass > 0) && currentAge == targetAge)
        Nwinds_to_strip++;
    }

  mpi_printf("GFM_WINDS_STRIPPING: found %d wind particles to strip.\n", Nwinds_to_strip);


  StarParticle = mymalloc("StarParticle", Nwinds_to_strip * sizeof(struct star_particle));

  Nwinds_to_strip = 0;
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      currentAge = (float) STP(i).BirthTime;
      if((P[i].Type == 4 && P[i].Mass > 0) && currentAge == targetAge)
        {
          StarParticle[Nwinds_to_strip].index = i;
          StarParticle[Nwinds_to_strip].NumNgb = 0;
          StarParticle[Nwinds_to_strip].NormSph = 0;
          StarParticle[Nwinds_to_strip].TotalMassReleased = -1; /* indicates that we have not found its hsml value */
          StarParticle[Nwinds_to_strip].TotalMetalMassReleased = 0;
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            StarParticle[Nwinds_to_strip].MetalMassReleased[k] = 0;
#ifdef GFM_DUST
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            {
              for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
                {
                  StarParticle[Nwinds_to_strip].DustMassReleased[chan][k] = 0.0;
                }

            }
#endif
#if defined(GFM_DUST) || defined(DUST_LIVE)
          StarParticle[Nwinds_to_strip].NumSNII = 0.0;
#endif
#ifdef GFM_CHEMTAGS
          for(k = 0; k < GFM_N_CHEM_TAGS; k++)
            StarParticle[Nwinds_to_strip].MetalMassReleasedChemTags[k] = 0;
#endif

          Nwinds_to_strip++;
        }
    }
}

void end_stripping()
{
  myfree(StarParticle);

  mpi_printf("GFM_WINDS_STRIPPING: done.\n");

  TIMER_STOP(CPU_GFM_ENRICH);
}
#endif


#ifdef GFM_PREENRICH

static double mass_fractions[GFM_N_CHEM_ELEMENTS];
static double metallicity;

void gfm_read_preenrich_table(char *fname)
{
  int n_preenrich_tab, elem;
  FILE *fdpreenrich;
  char elementname[100];
  double massfrac;

  mpi_printf("GFM_PREENRICH: reading preenrichment file '%s'...\n", fname);
  if(!(fdpreenrich = fopen(fname, "r")))
    terminate("GFM_PREENRICH: Cannot read preenrichment file '%s'\n", fname);

  n_preenrich_tab = 0;
  while(fscanf(fdpreenrich, "%s %lg", elementname, &massfrac) != EOF)
    {
      if(strcmp(elementname, "Metallicity") == 0)
        {
          mpi_printf("GFM_PREENRICH: preenrich Metallicity=%g\n", massfrac);
          metallicity = massfrac;
        }
      else
        {
          if((elem = element_index(elementname)) < 0)
            terminate("GFM_PREENRICH: element '%s' not found\n", elementname);

          mass_fractions[elem] = massfrac;
          mpi_printf("GFM_PREENRICH: preenrich element=%-10s (index=%02d) with massfrac=%g\n", elementname, elem, massfrac);
          if(n_preenrich_tab == GFM_N_CHEM_ELEMENTS)
            terminate("GFM_PREENRICH: too many elements\n");
          n_preenrich_tab++;
        }
    }
  fclose(fdpreenrich);

#ifdef GFM_NORMALIZED_METAL_ADVECTION
  double sum = 0;
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    sum += mass_fractions[i];
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    mass_fractions[i] /= sum;

  metallicity = 0;

  for(int i = 2; i < GFM_N_CHEM_ELEMENTS; i++)
    metallicity += mass_fractions[i];


  #ifdef WINDTUNNEL_FIXVARIABLESININJECTIONREGION
  All.metallicity = metallicity;
  for(int i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    All.mass_fractions[i]= mass_fractions[i];
  
  if(All.PreEnrichTime > All.TimeBegin)
    terminate("PreEnrichTime>All.TimeBegin, this is not allowed when using WINDTUNNEL_FIXVARIABLESININJECTIONREGION\n");
  #endif

  mpi_printf("GFM_PREENRICH: Metallicity renormalized to %g, hydrogen=%18.12g helium=%18.12g\n", metallicity, mass_fractions[0], mass_fractions[1]);

#else
  if(metallicity == 0)
    terminate("GFM_PREENRICH: Metallicity not set or set to zero in preenrich.txt!\n");
#endif

  mpi_printf("GFM_PREENRICH: read preenrichment file.\n");
}


void gfm_preenrich_gas(void)
{
  int i, j;

  mpi_printf("GFM_PREENRICH: now at time %g we carry out the preenrichment by imposing the specified abundance pattern.\n", All.Time);

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type == 0)
        {
          SphP[i].Metallicity = metallicity;
          SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            {
              SphP[i].MetalsFraction[j] = mass_fractions[j];
              SphP[i].MassMetals[j] = SphP[i].MetalsFraction[j] * P[i].Mass;
            }
#ifdef GFM_CHEMTAGS
	  double summetals = 0.0;
	  for(int j = 2; j < GFM_N_CHEM_ELEMENTS; j++)
	    summetals += SphP[i].MassMetals[j];
	  
          for(int j = 0; j < GFM_N_CHEM_TAGS; j++)
            {
              SphP[i].MassMetalsChemTags[j] = summetals / 3.0; // 3: SNIA, SNII, AGB
              SphP[i].MassMetalsChemTagsFraction[j] = SphP[i].MassMetalsChemTags[j] / P[i].Mass;
            }

#ifdef GFM_SPLITFE
          double massiron = SphP[i].MassMetals[element_index("Iron")];

          SphP[i].MassMetalsChemTags[GFM_FESNIA_CHEMTAG] = 0.5 * massiron;
          SphP[i].MassMetalsChemTags[GFM_FESNII_CHEMTAG] = 0.5 * massiron;
          SphP[i].MassMetalsChemTagsFraction[GFM_FESNIA_CHEMTAG] = SphP[i].MassMetalsChemTags[GFM_FESNIA_CHEMTAG] / P[i].Mass;
          SphP[i].MassMetalsChemTagsFraction[GFM_FESNII_CHEMTAG] = SphP[i].MassMetalsChemTags[GFM_FESNIA_CHEMTAG] / P[i].Mass;
#endif

#ifdef GFM_RPROCESS
          SphP[i].MassMetalsChemTags[GFM_NSNS_CHEMTAG] = SphP[i].MassMetals[element_index("Iron")];  /* set this to iron abundance in lack of other info */
          SphP[i].MassMetalsChemTagsFraction[GFM_NSNS_CHEMTAG] = SphP[i].MassMetalsChemTags[GFM_NSNS_CHEMTAG] / P[i].Mass;
#endif

#endif // end of GFM_CHEMTAGS
        }
    }
}
#endif



int element_index(char *element_name)
{
  int i;

  for(i = 0; i < GFM_N_CHEM_ELEMENTS; i++)
    if(strcmp(ElementNames[i], element_name) == 0)
      return i;

  /* element not found */
  return -1;
}



int get_element_index(char *table[20], int size, char *element_name)
{
  int i;

  for(i = 0; i < size; i++)
    if(strcmp(table[i], element_name) == 0)
      return i;

  /* element not found */
  return -1;
}



void get_initial_mass_fractions(MyDouble * mass_fractions)
{
  int index = element_index("Hydrogen");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_HYDROGEN;

  index = element_index("Helium");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_HELIUM;

  index = element_index("Carbon");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_CARBON;

  index = element_index("Nitrogen");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_NITROGEN;

  index = element_index("Oxygen");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_OXYGEN;

  index = element_index("Neon");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_NEON;

  index = element_index("Magnesium");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_MAGNESIUM;

  index = element_index("Silicon");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_SILICON;

  index = element_index("Iron");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_IRON;

#ifdef GFM_NORMALIZED_METAL_ADVECTION
  index = element_index("OtherMetals");
  if(index > -1)
    mass_fractions[index] = GFM_INITIAL_ABUNDANCE_OTHER;
#endif
}

void convert_metal_abundances_to_masses()
{
  int i, k;
  mpi_printf("GFM: Converting metal abundances to masses ... \n");

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)        /* gas particle */
        {
          SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            SphP[i].MassMetals[k] = SphP[i].MetalsFraction[k] * P[i].Mass;

#ifdef GFM_CHEMTAGS
          for(k = 0; k < GFM_N_CHEM_TAGS; k++)
            {
              SphP[i].MassMetalsChemTagsFraction[k] = SphP[i].MassMetalsChemTags[k];
              SphP[i].MassMetalsChemTags[k] *= P[i].Mass;
            }
#endif
        }

      if(P[i].Type == 4)        /* star particle */
        {
          for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
            STP(i).MassMetals[k] *= P[i].Mass;

#ifdef GFM_CHEMTAGS
          for(k = 0; k < GFM_N_CHEM_TAGS; k++)
            STP(i).MassMetalsChemTags[k] *= P[i].Mass;
#endif
        }
    }
}

MyFloat ejecta_mass_for_bin(MyFloat ** ejecta_spline, MyFloat *** yield_spline, int elem_index, MyFloat * initial_metals, int iz_low, int iz_high, int imass, MyFloat dz)
{
  return initial_metals[elem_index] * ((1 - dz) * ejecta_spline[iz_low][imass]
                                       + dz * ejecta_spline[iz_high][imass]) + (1 - dz) * yield_spline[iz_low][elem_index][imass] + dz * yield_spline[iz_high][elem_index][imass];
}

#endif
