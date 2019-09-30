#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "../allvars.h"
#include "../proto.h"
#include "sgchem_proto.h"
#include "f2c.h"

void init_chemistry(void)
{
  /* Initialize parameters stored in coolr and cooli common blocks */
  INIT_CHEMISTRY_PARAMETERS(&All.CarbAbund, &All.OxyAbund, &All.MAbund, &All.InitDustTemp,
                            &All.UVFieldStrength, &All.DustToGasRatio, &All.CosmicRayIonRate,
                            &All.InitRedshift, &All.ExternalDustExtinction, &All.H2FormEx, &All.H2FormKin, &All.PhotoApprox, &All.ISRFOption);

  /* Initialize chemical rates */
  COOLINMO();
  CHEMINMO();

  /* ODE integrator tolerances */
  INIT_TOLERANCES();
  return;
}

void evolve_chemistry(void)
{
  double dt;
  int idx, i;


  /* Loop over active particles, and do chemical evolution for each one */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      evolve_chemistry_for_single_cell(dt, i);

    }

  return;
}

/* Chemical evolution, radiative heating & cooling */
void evolve_chemistry_for_single_cell(double dt, int index)
{
  double time, rho, energy, divv, dl, yn;
  double rho_cgs, redshift, dust_temp;
  double non_eq_abundances[SGCHEM_NUM_SPECIES];
  double column_density_projection[NPIX];
  double column_density_projection_H2[NPIX];
  double column_density_projection_CO[NPIX];
  double a, a3inv, hubble_a, hubble_param, hubble_param2;
  double column_correction_cosmic, column_correction_cgs;
  int i, id;

#if NO_CHEMISTRY
  return;
#else
  /* Check size of timestep */
  if(dt < 0)
    {
      printf("Error: negative timestep in evolve_chemistry\n");
      endrun(801);
    }
  else if(dt == 0)
    {
      return;
    }

  /* Set current index and ID in cool.h */
  id = P[index].ID;
  SET_INDEX_ID_FOR_CHEMISTRY(&index, &id);

  /* Starting values */
  time = dt;
  rho = SphP[index].Density;
  divv = SphP[index].DivVel;

  /* Estimate of typical distance from centre to edge of local grid-cell, used for computing
   * local contribution to shielding.
   *
   * XXX: should we make it possible to switch this on or off?
   */
  dl = 0.5 * pow(SphP[index].Volume, 1. / 3.);

  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      non_eq_abundances[i] = SphP[index].TracAbund[i];
    }

  /* If TREE_RAD and friends are not in use, then these are all single-element arrays */
  column_density_projection[0] = column_density_projection_H2[0] = column_density_projection_CO[0] = 0.0;

#ifdef TREE_RAD
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] = SphP[index].Projection[i];
#ifdef TREE_RAD_H2
      column_density_projection_H2[i] = SphP[index].ProjectionH2[i];
#endif
#ifdef TREE_RAD_CO
      column_density_projection_CO[i] = SphP[index].ProjectionCO[i];
#endif
    }
#endif

  /* Convert from comoving to physical units */
  if(All.ComovingIntegrationOn)
    {
      a = All.Time;
      a3inv = 1 / (a * a * a);
      hubble_a = hubble_function(a);
      hubble_param = All.HubbleParam;
      hubble_param2 = hubble_param * hubble_param;
      redshift = 1 / (a - 1);
    }
  else
    {
      a = a3inv = hubble_a = hubble_param = hubble_param2 = 1;
      redshift = All.InitRedshift;
    }

  dt *= 1 / hubble_a / hubble_param;
  rho *= a3inv * hubble_param2;
  energy = SphP[index].Utherm * rho;    /* Energy here is energy density, AREPO evolves specific energy */
  divv *= 1.0 / a;
  if(All.ComovingIntegrationOn)
    divv += 3 * hubble_a;
  dl *= a;

  column_correction_cosmic = hubble_param2 / (a * a);
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] *= column_correction_cosmic;
      column_density_projection_H2[i] *= column_correction_cosmic;
      column_density_projection_CO[i] *= column_correction_cosmic;
    }

  /* Convert from AREPO code units to cgs */
  dt *= All.UnitTime_in_s;
  rho_cgs = rho * All.UnitDensity_in_cgs;       /* We keep rho in code units as we need it later */
  energy *= All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  divv *= 1.0 / All.UnitTime_in_s;
  dl *= All.UnitLength_in_cm;

  column_correction_cgs = All.UnitDensity_in_cgs * All.UnitLength_in_cm;
  for(i = 0; i < NPIX; i++)
    {
      column_density_projection[i] *= column_correction_cgs;
      column_density_projection_H2[i] *= column_correction_cgs;
      column_density_projection_CO[i] *= column_correction_cgs;
    }

  /* Compute derived units */
  yn = rho_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);

  EVOLVE_ABUNDANCES(&time, &dl, &yn, &divv, &energy, &redshift, non_eq_abundances, column_density_projection, column_density_projection_H2, column_density_projection_CO, &dust_temp);

  /* Convert from cgs to AREPO code units */
  energy /= All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);       /* Energy density in Arepo code units */
  SphP[index].Utherm = energy / rho;

  if(id == 1000)
    {
      printf("Energy: u %g energy % g \n", SphP[index].Utherm, energy);

    }


  /* Set new dust temperature */
  SphP[index].DustTemp = dust_temp;

  /* Update evolved values */
  SphP[index].Utherm = energy;

  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      SphP[index].TracAbund[i] = non_eq_abundances[i];
      SphP[index].MassTracAbund[i] = P[index].Mass * SphP[index].TracAbund[i];
    }
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  SphP[index].TracAbund[IHATOM] = 1.0 - 2.0 * SphP[index].TracAbund[IH2] - SphP[index].TracAbund[IHP];
  if(SphP[index].TracAbund[IHATOM] < 0.0)
    {
      SphP[index].TracAbund[IHATOM] = 0.0;
    }
  SphP[index].MassTracAbund[IHATOM] = SphP[index].TracAbund[IHATOM] * P[index].Mass;
#endif
#if CHEMISTRYNETWORK == 5
  SphP[index].TracAbund[ICP] = All.CarbAbund - SphP[index].TracAbund[ICO];
  if(SphP[index].TracAbund[ICP] < 0.0)
    {
      SphP[index].TracAbund[ICP] = 0.0;
    }
  SphP[index].MassTracAbund[ICP] = SphP[index].TracAbund[ICP] * P[index].Mass;
#endif

  /* Update pressure */
  set_pressure_of_cell(index);

  return;
#endif /* No CHEMISTRY */
}
