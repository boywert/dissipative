/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/update_primitive_variables.c
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
#include "proto.h"

#include "helm_eos.h"
#include "opal_eos.h"
#include <gsl/gsl_linalg.h>

void update_primitive_variables(void)
{
  TIMER_START(CPU_CELL_UPDATES);
  
#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double InitialCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      InitialCREnergy += SphP[i].CR_Energy;
#endif

  struct pv_update_data pvd;
  int idx, i;

  if(All.ComovingIntegrationOn)
    {
      pvd.atime = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  pvd.count_keep_entropy = pvd.count_update_entropy = 0;
#endif

#pragma omp parallel for private(idx,i)
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      do_validity_checks(P, SphP, i, &pvd);

      do_special_boundaries(P, SphP, i, &pvd);


#ifdef SPECIAL_RELATIVITY
      update_primitive_variables_special_relativity(P, SphP, i, &pvd);
#else
#ifdef GENERAL_RELATIVITY
      update_primitive_variables_general_relativity(P, SphP, i, &pvd);
#else
      update_primitive_variables_single(P, SphP, i, &pvd);

      update_internal_energy(P, SphP, i, &pvd);

      set_pressure_of_cell_internal(P, SphP, i);        /* calculate the pressure from Density and Utherm (and composition) */
#endif
#endif

#ifdef SGCHEM
      do_chemical_abund_checks(P, SphP, i, &pvd);
#endif

      SphP[i].OldMass = P[i].Mass;

#ifdef MRT_LSF_GRADIENTS
      //      for(int num1=0;num1<MRT_BINS;num1++)
      //SphP[i].OldCons_DensPhot[num1] = SphP[i].Cons_DensPhot[num1] ;
#ifdef MRT_COMOVING
      SphP[i].Old_Vel[0] = P[i].Vel[0] ;
      SphP[i].Old_Vel[1] = P[i].Vel[1] ;
      SphP[i].Old_Vel[2] = P[i].Vel[2] ;
#endif
#endif

      SphP[i].TimeLastPrimUpdate = All.Time;
    }

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
  long long tot_count_keep_entropy, tot_count_update_entropy;

  sumup_large_ints(1, &pvd.count_keep_entropy, &tot_count_keep_entropy);
  sumup_large_ints(1, &pvd.count_update_entropy, &tot_count_update_entropy);

  mpi_printf("Keep-Entropy=%llu  Update-Entropy=%llu\n", tot_count_keep_entropy, tot_count_update_entropy);
#endif

#ifdef COSMIC_RAYS_EXTRA_DIAGNOSTICS
  double FinalCREnergy = 0;
  for(int i=0; i<NumGas; i++)
    if(P[i].Mass != 0 && P[i].ID != 0 && P[i].Type == 0)
      FinalCREnergy += SphP[i].CR_Energy;
  
  double dCREnergy = FinalCREnergy - InitialCREnergy;
  All.TotalCREnergyUpdatePrims += dCREnergy;
#endif

  TIMER_STOP(CPU_CELL_UPDATES);
}

void set_pressure_of_cell(int i)
{
  set_pressure_of_cell_internal(P, SphP, i);
}

void set_pressure_of_cell_internal(struct particle_data *localP, struct sph_particle_data *localSphP, int i)
{
#ifdef DG
  localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#ifndef FIX_MEAN_VALUES
  if(localSphP[i].Pressure < Epsilon_p)
    {
      printf("cell: %d, pressure: %g\n", i, localSphP[i].Pressure);
      terminate("Negative pressure detected!\n");
    }
#endif

  return;
#endif

#ifdef DISSIPATIVE
  double alpha_D = 0.01;
  double dark_photon_mass = 1e-10;
  printf("update pressure %d U = %f\n",i, localSphP[i].Utherm);
  
  localSphP[i].Pressure = localSphP[i].Density * localSphP[i].EOSTemperature
    + 2.0 * M_PI * alpha_D * localSphP[i].Density
    * localSphP[i].Density / (dark_photon_mass*dark_photon_mass);
  return;
#endif  // DISSIPATIVE
#ifdef ISOTHERM_EQS
  localSphP[i].Pressure = localSphP[i].Density * All.IsoSoundSpeed * All.IsoSoundSpeed;
#else

#ifdef LOCALLY_ISOTHERM_DISK
  double csnd;
  csnd = get_isotherm_disk_sound_speed(i);
  localSphP[i].Pressure = localSphP[i].Density * csnd * csnd;
#else

  if(localSphP[i].Utherm >= 0)
    {
#if !defined(USE_SFR) || defined(LOCAL_FEEDBACK)

#ifdef VARIABLE_GAMMA
      localSphP[i].GammaE = GAMMA;
      localSphP[i].GammaC = GAMMA;
#endif

#ifdef TGCHEM
      localSphP[i].Pressure = (localSphP[i].Gamma - 1) * localSphP[i].Density * localSphP[i].Utherm;
#else
#ifdef EOS_DEGENERATE
      struct eos_result res;

      if(eos_calc_egiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].Utherm, &localSphP[i].EOSTemperature, &res) == -2)
        {
          localSphP[i].EOSTemperature = -1.;
          if(eos_calc_egiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].Utherm, &localSphP[i].EOSTemperature, &res) != 0)
            {
              if(res.temp < 1000.)
                localSphP[i].EOSTemperature = 1000;
              else if(res.temp > 1e10)
                localSphP[i].EOSTemperature = 1e10;

              eos_calc_tgiven(localSphP[i].Density, localSphP[i].Composition, localSphP[i].EOSTemperature, &res);
            }
        }

      if(localSphP[i].EOSTemperature <= 1.001 * 1000. || localSphP[i].EOSTemperature >= 1e10 * 0.999)
        {
          /* looks like the equation of state corrected the internal energy, as we run out of the table */
          EgyInjection += res.e.v - localSphP[i].Utherm;
          localSphP[i].Energy += (res.e.v - localSphP[i].Utherm) * localP[i].Mass;
          localSphP[i].Utherm = res.e.v;
        }

#ifdef OUTPUT_ENTROPY
      localSphP[i].Entropy = res.s.v;
#endif

      localSphP[i].Pressure = res.p.v;
      localSphP[i].GammaE = localSphP[i].Pressure / localSphP[i].Utherm / localSphP[i].Density + 1.0;
      localSphP[i].GammaC = (res.p.drho + res.temp * gsl_pow_2(res.p.dtemp / localSphP[i].Density) / res.e.dtemp) * localSphP[i].Density / localSphP[i].Pressure;

      if(localSphP[i].GammaC <= 0)
        {
          printf("Pressure=%g, Temperature=%g, Density=%g, dpdr=%g, dedt=%g.\n", res.p.v, res.temp, localSphP[i].Density, res.p.drho, res.e.dtemp);
          print_particle_info(i);
          terminate("bla");
        }
#else
#ifdef EOS_OPAL
      struct opal_eos_result res;

      /* limit and renormalize composition */
      double xsum = 0.0;
      int j;
      /* X can be 0.8 maximum */
      localSphP[i].Composition[0] = dmin(0.8, localSphP[i].Composition[0]);
      for(j = 0; j < EOS_NSPECIES; j++)
        {
          localSphP[i].Composition[j] = dmin(1.0, dmax(0.0, localSphP[i].Composition[j]));
          xsum += localSphP[i].Composition[j];
        }
      for(j = 0; j < EOS_NSPECIES; j++)
        localSphP[i].Composition[j] /= xsum;
      /* check rho boundaries */
      /*
      if(localSphP[i].Density < opal_rhomin)
        {
          //printf("Error in OPAL EOS: density too small at %e.\n", localSphP[i].Density);
          localSphP[i].Density = opal_rhomin;
        }
      if(localSphP[i].Density > opal_rhomax)
        {
          //printf("Error in OPAL EOS: density too large at %e.\n", localSphP[i].Density);
          localSphP[i].Density = opal_rhomax;
        }
        */

      /* call EOS */
      if(opaleos_egiven(localSphP[i].Composition[0], localSphP[i].Density,
                        localSphP[i].Utherm, &localSphP[i].EOSTemperature, EOS_PRESSURE | EOS_ENERGY | EOS_CHIT | EOS_CHIR | EOS_CV | EOS_GAMMA1 | EOS_RADIATION, &res) < 0)
        {
          /* printf("Correcting internal energy from %e to %e (total energy: %e;\
            rel. change: %e).\n", localSphP[i].Utherm, res.e, localSphP[i].Energy, res.e / localSphP[i].Energy); */
          localSphP[i].Utherm = res.e;
          /* print_particle_info(i); */
        }

      localSphP[i].Pressure = res.p;

      localSphP[i].GammaE = localSphP[i].Pressure / localSphP[i].Utherm / localSphP[i].Density + 1.0;
      localSphP[i].GammaC = res.gamma1;

      /* check results */
      if(!gsl_finite(localSphP[i].GammaC) || !gsl_finite(localSphP[i].GammaE))
        {
          printf("Error in EOS for particle %d.\n", i);
          print_particle_info(i);
        terminate("Infinity encountered in EOS.")}

#else
      /* this block applies to ordinary gas dynamics, without SFR */
      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#if defined(MRT_IR_PHOTON_TRAPPING) && defined(MRT_RADIATION_PRESSURE)
      for(int num1=UV_BINS;num1<(UV_BINS+IR_BINS);num1++)
	localSphP[i].Pressure += localSphP[i].Trapped_DensPhot[num1]*c_internal_units/(2.99792458e10 / All.UnitVelocity_in_cm_per_s)/3.0 ;
#endif

#endif /* EOS_OPAL */
#endif
#endif // TGCHEM

#ifdef JEANS_PRESSURE_LIMIT
      double jeans_pressure, celldim, ncells;
      ncells = JEANS_PRESSURE_LIMIT;
      //celldim = 2.0 * get_cell_radius(i);
      celldim = 2.0 * pow(All.MinVolume * 3.0 / 4.0 / M_PI / 2.0, 1.0 / 3.0);
      jeans_pressure = ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / GAMMA;
      localSphP[i].Pressure = fmax(localSphP[i].Pressure, jeans_pressure);

#ifdef MAKE_PRES_UTHERM_CONSISTENT
      double utherm_from_pres = localSphP[i].Pressure / GAMMA_MINUS1 / localSphP[i].Density;
      double utherm_diff = utherm_from_pres - localSphP[i].Utherm;
      localSphP[i].Utherm += utherm_diff;
      localSphP[i].Energy += utherm_diff * localP[i].Mass;
#endif

#endif

#ifdef JEANS_TOTPRESSURE_LIMIT

      double jeans_pressure, celldim, ncells;
      ncells = JEANS_PRESSURE_LIMIT;
      //celldim = 2.0 * get_cell_radius(i);
      celldim = 2.0 * pow(All.MinVolume * 3.0 / 4.0 / M_PI / 2.0, 1.0 / 3.0);
      jeans_pressure = ncells * ncells * All.G * celldim * celldim * localSphP[i].Density * localSphP[i].Density / M_PI / GAMMA;


      double tot_pressure = GAMMA_MINUS1 * localSphP[i].Energy / localSphP[i].Volume;
      tot_pressure = fmax(jeans_pressure, tot_pressure);




#endif


#else /* now comes the SFR treatment */

      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#endif // USE_SFR
    }                           //end utherm >= 0
  else
    localSphP[i].Pressure = 0;
#endif // LOCAL_ISOTHERM_DISK
#endif // ISOTHERM_EQS

#ifdef COSMIC_RAYS
  localSphP[i].CR_Pressure = (All.GammaCR - 1.0) * localSphP[i].CR_SpecificEnergy * localSphP[i].Density / All.cf_atime;
#ifdef COSMIC_RAYS_STREAMING
  localSphP[i].CR_Chi = compute_chi(i);
#endif
#endif

#ifdef ENFORCE_JEANS_STABILITY_OF_CELLS
#ifndef ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
#if defined(USE_SFR) && !defined(FM_SFR) && !defined(LOCAL_FEEDBACK)
  if(get_starformation_rate(i) == 0)
#endif
    {
#endif
      localSphP[i].Pressure = dmax(localSphP[i].Pressure, GAMMA_MINUS1 * localSphP[i].Density * 2 * All.G * localP[i].Mass / (All.cf_atime * All.ForceSoftening[localP[i].SofteningType]));
#ifndef ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS
    }
#endif
#endif
}

void do_validity_checks(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
  if(localP[i].Mass < 0)
    {
      printf("very bad...i=%d ID=%d mass=%g oldMass=%g utherm=%g pos=%g|%g|%g\n",
             i, (int) localP[i].ID, localP[i].Mass, localSphP[i].OldMass, localSphP[i].Utherm, localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2]);


      terminate("stop");
    }
}

#ifdef SGCHEM
void do_chemical_abund_checks(struct particle_data *localP, struct sph_particle_data *localSphP, int index, struct pv_update_data *pvd)
{
  double non_eq_abundances_i, carb_abund, oxy_abund;
  int i;

#ifdef SGCHEM_VARIABLE_Z
  carb_abund = localSphP[index].CarbAbund;
  oxy_abund  = localSphP[index].OxyAbund;
#else
  carb_abund = All.CarbAbund;
  oxy_abund  = All.OxyAbund;
#endif

  // Check that the chemistry has sane values after the advection.
  for(i = 0; i < SGCHEM_NUM_SPECIES; i++)
    {
      non_eq_abundances_i = SphP[index].TracAbund[i];

      if(non_eq_abundances_i < 0.0)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c negative abundance from advection, species = %d abundance = %g\n", i, non_eq_abundances_i);
          printf("update_primitive_variables.c: Setting abundance to +1e-20 (this might not help!)\n");
#endif
          non_eq_abundances_i = 1e-20;
        }

      if(i == IH2 && non_eq_abundances_i > 0.5)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c H2 abundance greater than 0.5; abundance = %g\n", non_eq_abundances_i);
          non_eq_abundances_i = 0.5;
#endif
        }

      if(i == IHP && non_eq_abundances_i > 1.0)
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c HP abundance greater than 1.0; abundance = %g\n", non_eq_abundances_i);
#endif
          non_eq_abundances_i = 1.0;
        }

      if(i == ICO && non_eq_abundances_i > fmin(carb_abund, oxy_abund))
        {
#ifdef DEBUG_SGCHEM
          printf("update_primitive_variables.c CO abundance greater than Carbon & Oxygen abundances; abundance = %g\n", non_eq_abundances_i);
#endif
          non_eq_abundances_i = fmin(carb_abund, oxy_abund);
        }

      SphP[index].TracAbund[i] = non_eq_abundances_i;
      SphP[index].MassTracAbund[i] = P[index].Mass * SphP[index].TracAbund[i];

    }

}
#endif

void do_special_boundaries(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifdef SPECIAL_BOUNDARY
  double dt;

  // ID <-3 cells are buffer cells and are never updated
  if((localP[i].ID <= -3) && (All.Time > All.TimeBegin))
    return;

  if((localP[i].ID == -2) && (All.SpecialBoundaryType == 3 || All.SpecialBoundaryType == 4))
    {
      dt = (localP[i].TimeBinHydro ? (((integertime) 1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval;
      boundary_get_velocity(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2], &localP[i].Vel[0], &localP[i].Vel[1], &localP[i].Vel[2], dt);
      localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
      localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
      localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
      /*update vectorial quantities only, not the scalars */
    }
#endif

 

#ifdef WINDTUNNEL_FIXVARIABLESININJECTIONREGION

#ifdef TWODIMS
  if(   (localP[i].Pos[0]<All.InjectionRegion) || (localP[i].Pos[1]<All.InjectionRegion)  ) //Create an injection coordinate at the upstream boundary and at the edge of the order coordinates (the latter is done to avoid shock reflection or crossing at boundary )
#else
  if(   (localP[i].Pos[0]<All.InjectionRegion) || (localP[i].Pos[1]<All.InjectionRegion) || (localP[i].Pos[2]<All.InjectionRegion) )
#endif

  {
      localSphP[i].Density = All.InjectionDensity;
      localP[i].Vel[0] = 0;
      localP[i].Vel[1] = 0;
      localP[i].Vel[2] = 0;
      localP[i].Vel[WINDTUNNEL_COORD] = All.InjectionVelocity;
      localSphP[i].Utherm = All.InjectionUtherm;

      localP[i].Mass = localSphP[i].Density * localSphP[i].Volume;
      localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
      localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
      localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
      localSphP[i].Energy =
          localP[i].Mass * pvd->atime * pvd->atime * localSphP[i].Utherm +
          0.5 * localP[i].Mass * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);

      localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;

#ifdef GFM_NORMALIZED_METAL_ADVECTION
      localSphP[i].Metallicity = All.metallicity;
      localSphP[i].MassMetallicity = SphP[i].Metallicity * localP[i].Mass;
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
      {
          localSphP[i].MetalsFraction[j] = All.mass_fractions[j];
          localSphP[i].MassMetals[j] = localSphP[i].MetalsFraction[j] * localP[i].Mass;
      }      
#endif

  }
#endif




  /*Check if we have buffer (inflow-outflow)cells with prescribed values */
#if defined(BOUNDARY_INFLOWOUTFLOW_MINID) && defined(BOUNDARY_INFLOWOUTFLOW_MAXID)
  if(localP[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && localP[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
    {
#ifdef WINDTUNNEL
      localSphP[i].Density = All.InjectionDensity;
      localP[i].Vel[0] = 0;
      localP[i].Vel[1] = 0;
      localP[i].Vel[2] = 0;
      localP[i].Vel[WINDTUNNEL_COORD] = All.InjectionVelocity;
      localSphP[i].Utherm = All.InjectionUtherm;

 #ifdef GFM_NORMALIZED_METAL_ADVECTION
      localSphP[i].Metallicity = All.metallicity;
      localSphP[i].MassMetallicity = SphP[i].Metallicity * localP[i].Mass;
      for(int j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
      {
            localSphP[i].MetalsFraction[j] = All.mass_fractions[j];
                localSphP[i].MassMetals[j] = localSphP[i].MetalsFraction[j] * localP[i].Mass;
      }      
 #endif

#endif

#ifdef SPECIAL_BOUNDARY
      /*If the special boundaries are at the same time buffer cells */
      if((localP[i].ID == -1) || (localP[i].ID == -2))
        {
          dt = (localP[i].TimeBinHydro ? (((integertime) 1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval;
          boundary_get_velocity(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2], &localP[i].Vel[0], &localP[i].Vel[1], &localP[i].Vel[2]);
        }
#endif

#if defined(WINDTUNNEL) || defined(SPECIAL_BOUNDARY)
      if(localP[i].ID >= BOUNDARY_INFLOWOUTFLOW_MINID && localP[i].ID < BOUNDARY_INFLOWOUTFLOW_MAXID)
        {
          localP[i].Mass = localSphP[i].Density * localSphP[i].Volume;
          localSphP[i].Momentum[0] = localP[i].Vel[0] * localP[i].Mass;
          localSphP[i].Momentum[1] = localP[i].Vel[1] * localP[i].Mass;
          localSphP[i].Momentum[2] = localP[i].Vel[2] * localP[i].Mass;
#ifndef ISOTHERM_EQS
          localSphP[i].Energy =
            pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
            0.5 * localP[i].Mass * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
          localSphP[i].Pressure = GAMMA_MINUS1 * localSphP[i].Density * localSphP[i].Utherm;
#endif
        }
#endif
    }
#endif
}

void update_primitive_variables_single(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifdef DG

  localSphP[i].Density = localP[i].Mass / localSphP[i].Volume;

#ifndef FIX_MEAN_VALUES
  if(localSphP[i].Density < Epsilon_rho)
    {
      printf("cell: %d, density: %g\n", i, localSphP[i].Density);
      terminate("Negative density detected!\n");
    }
#endif

  localP[i].Vel[0] = localSphP[i].Momentum[0] / localP[i].Mass;
  localP[i].Vel[1] = localSphP[i].Momentum[1] / localP[i].Mass;
  localP[i].Vel[2] = localSphP[i].Momentum[2] / localP[i].Mass;

  return;
#endif
  localSphP[i].Density = localP[i].Mass / localSphP[i].Volume;

#ifdef RT_ADVECT
  int j;
  for(j = 0; j < RT_N_DIR; j++)
    localSphP[i].DensPhot[j] = localSphP[i].Photons[j] / localSphP[i].Volume;
#endif

#ifdef MRT
  double nH_times_volume ;

  nH_times_volume  = localP[i].Mass ;

  if(localSphP[i].nHI == 0.0 && localSphP[i].nHII == 0.0)
    {
      double y_fac = (1.0 - HYDROGEN_MASSFRAC) / 4.0 / HYDROGEN_MASSFRAC;
      double nh_fac = 1.0 / (1.0 + y_fac) ;
      localSphP[i].HI = 0.999999 ;
      localSphP[i].HII = 0.000001 ;
      localSphP[i].nHI = localSphP[i].HI * nH_times_volume ;
      localSphP[i].nHII = localSphP[i].HII * nH_times_volume ;
#ifdef MRT_INCLUDE_HE
      localSphP[i].HeI = 0.999998*y_fac ;
      localSphP[i].HeII = 0.000001*y_fac ; 
      localSphP[i].HeIII = 0.000001*y_fac ;
      localSphP[i].Ne = localSphP[i].HII + localSphP[i].HeII + 2.0 * localSphP[i].HeIII ;
      localSphP[i].nHeI = localSphP[i].HeI * nH_times_volume ;
      localSphP[i].nHeII = localSphP[i].HeII * nH_times_volume ;
      localSphP[i].nHeIII = localSphP[i].HeIII * nH_times_volume ;
#else
      localSphP[i].Ne = localSphP[i].HII ;
#endif
      localSphP[i].ne = localSphP[i].Ne * nH_times_volume ;
    }
  else
    {
      localSphP[i].HI = localSphP[i].nHI / nH_times_volume ;
      localSphP[i].HII = localSphP[i].nHII / nH_times_volume ;
      localSphP[i].Ne = localSphP[i].ne / nH_times_volume ;
#ifdef MRT_INCLUDE_HE
      localSphP[i].HeI = localSphP[i].nHeI / nH_times_volume ;
      localSphP[i].HeII = localSphP[i].nHeII / nH_times_volume ;
      localSphP[i].HeIII = localSphP[i].nHeIII / nH_times_volume ;
#endif
    }

  if(isnan(localSphP[i].nHI))
    terminate("Nan nHI\n") ;



  if((localSphP[i].HI + localSphP[i].HII > 1.1) || (localSphP[i].HI + localSphP[i].HII < 0.95))
    terminate("i = %d task = %d \n\n Passive advection gone wrong HI = %g \n HII = %g \t total = %g\n nhI, nhII, Mass = %g %g %g \n", i, ThisTask, localSphP[i].HI, localSphP[i].HII, localSphP[i].HI + localSphP[i].HII, localSphP[i].nHI, localSphP[i].nHII, localP[i].Mass) ;


#ifdef MRT_IR
  set_kappa_times_rho_IR(i, localSphP) ;
#endif


  
#ifdef MRT_COMOVING
  cell_do_lorentz_boost(i, localSphP, P[i].Vel[0]-localSphP[i].Old_Vel[0], P[i].Vel[1]-localSphP[i].Old_Vel[1], P[i].Vel[2]-localSphP[i].Old_Vel[2]) ;
#endif

#endif

  if(localP[i].Mass > 0)
    {
      localP[i].Vel[0] = localSphP[i].Momentum[0] / localP[i].Mass;
      localP[i].Vel[1] = localSphP[i].Momentum[1] / localP[i].Mass;
      localP[i].Vel[2] = localSphP[i].Momentum[2] / localP[i].Mass;

#ifdef COSMIC_RAYS
      if(SphP[i].CR_Energy / SphP[i].Volume < All.MinimumCREnergyDensity)
        {
          SphP[i].CR_Energy = All.MinimumCREnergyDensity * SphP[i].Volume;
        }
#endif

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        {
#ifdef MHD_THERMAL_ENERGY_SWITCH
          if(k != ScalarIndex.Etherm)   /* don't always update thermal energy */
#endif
            *(MyFloat *) (((char *) (&localSphP[i])) + scalar_elements[k].offset) = *(MyFloat *) (((char *) (&localSphP[i])) + scalar_elements[k].offset_mass) / localP[i].Mass;
        }
#endif

#ifdef GFM_NORMALIZED_METAL_ADVECTION
      double mass_metals = 0;

      for(int j = 2; j < GFM_N_CHEM_ELEMENTS; j++)
        mass_metals += localSphP[i].MassMetals[j];

      localSphP[i].MassMetallicity = mass_metals;
      localSphP[i].Metallicity = mass_metals / localP[i].Mass;
#endif


#ifdef GFM_CHEMTAGS
      for(int k = 0; k < GFM_N_CHEM_TAGS; k++)
        localSphP[i].MassMetalsChemTagsFraction[k] = localSphP[i].MassMetalsChemTags[k] / localP[i].Mass;
#endif  // en


#ifdef ACTIVE_CELL_SPIN
      {
        /* Spin = I * Omega -> solve for Omega */
        int j, k;
        double momInertia[9];

        /* we have to move the tensor from the center of mass to the center of rotation of the cell */
        {
          double xtmp, ytmp, ztmp;

          double dx = NEAREST_X(localSphP[i].CenterOffset[0]);
          double dy = NEAREST_Y(localSphP[i].CenterOffset[1]);
          double dz = NEAREST_Z(localSphP[i].CenterOffset[2]);

          momInertia[0] = localSphP[i].MomentOfInertia[0][0] + localSphP[i].Volume * (dy * dy + dz * dz);
          momInertia[4] = localSphP[i].MomentOfInertia[1][1] + localSphP[i].Volume * (dx * dx + dz * dz);
          momInertia[8] = localSphP[i].MomentOfInertia[2][2] + localSphP[i].Volume * (dx * dx + dy * dy);
          momInertia[1] = localSphP[i].MomentOfInertia[0][1] - localSphP[i].Volume * dx * dy;
          momInertia[3] = localSphP[i].MomentOfInertia[1][0] - localSphP[i].Volume * dx * dy;
          momInertia[2] = localSphP[i].MomentOfInertia[0][2] - localSphP[i].Volume * dx * dz;
          momInertia[6] = localSphP[i].MomentOfInertia[2][0] - localSphP[i].Volume * dx * dz;
          momInertia[5] = localSphP[i].MomentOfInertia[1][2] - localSphP[i].Volume * dy * dz;
          momInertia[7] = localSphP[i].MomentOfInertia[2][1] - localSphP[i].Volume * dy * dz;
        }

        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            momInertia[j * 3 + k] *= localSphP[i].Density;

        gsl_matrix_view A = gsl_matrix_view_array(momInertia, 3, 3);
        gsl_vector_view b = gsl_vector_view_array(SphP[i].Spin, 3);
        gsl_vector *x = gsl_vector_alloc(3);

        int s;
        gsl_permutation *p = gsl_permutation_alloc(3);
        gsl_linalg_LU_decomp(&A.matrix, p, &s);
        gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);

        for(j = 0; j < 3; j++)
          localSphP[i].Omega[j] = gsl_vector_get(x, j);

        gsl_permutation_free(p);
        gsl_vector_free(x);
      }
#endif

#if defined(GFM_STELLAR_EVOLUTION) && (GFM_STELLAR_EVOLUTION==1)
      /* re-normalizing the primitive metallicity variables, because in this mode
         the gas MassMetallicity & MassMetals are updated (from stellar evolution) but not gas Mass */
      double mass_by_metals = 0;
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          mass_by_metals += localSphP[i].MassMetals[k];
        }
      localSphP[i].Metallicity *= (localP[i].Mass / mass_by_metals);
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          localSphP[i].MetalsFraction[k] *= (localP[i].Mass / mass_by_metals);
        }
#ifdef GFM_DUST
      /* TODO? */
      for(int k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
        {
          localSphP[i].MetalsDustFraction[k] *= (localP[i].Mass / mass_by_metals);
        }
#endif

#endif

#ifdef MHD
      localSphP[i].B[0] = localSphP[i].BConserved[0] / localSphP[i].Volume;
      localSphP[i].B[1] = localSphP[i].BConserved[1] / localSphP[i].Volume;
      localSphP[i].B[2] = localSphP[i].BConserved[2] / localSphP[i].Volume;
#endif

#ifdef MHD_CT
      localSphP[i].A[0] = localSphP[i].AConserved[0] / localSphP[i].Volume;
      localSphP[i].A[1] = localSphP[i].AConserved[1] / localSphP[i].Volume;
      localSphP[i].A[2] = localSphP[i].AConserved[2] / localSphP[i].Volume;
#endif

#ifdef TRACER_FIELD
      localSphP[i].Tracer = localSphP[i].ConservedTracer / localP[i].Mass;
#endif
    }
  else                          /* P[i].Mass <= 0 */
    {
      localP[i].Vel[0] = 0;
      localP[i].Vel[1] = 0;
      localP[i].Vel[2] = 0;

#ifdef MAXSCALARS
      for(int k = 0; k < N_Scalar; k++)
        *(MyFloat *) (((char *) (&localSphP[i])) + scalar_elements[k].offset) = 0;
#endif
#ifdef TRACER_FIELD
      localSphP[i].Tracer = 0;
#endif
    }
}


void update_internal_energy(struct particle_data *localP, struct sph_particle_data *localSphP, int i, struct pv_update_data *pvd)
{
#ifdef DG
  //rho>=Epsilon_rho
  localSphP[i].Utherm = (localSphP[i].Energy / localP[i].Mass - 0.5 * (localP[i].Vel[0] * localP[i].Vel[0] +
                                                                       localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2])) / (pvd->atime * pvd->atime);
  return;
#endif


#ifndef ISOTHERM_EQS
  double ulimit;
#if (defined(MEASURE_DISSIPATION_RATE) && defined(USE_ENTROPY_FOR_COLD_FLOWS)) || defined(TGCHEM)
  double dt = (localP[i].TimeBinHydro ? (((integertime) 1) << localP[i].TimeBinHydro) : 0) * All.Timebase_interval / pvd->hubble_a;
#endif

#ifdef TGCHEM
  double uold = localSphP[i].Utherm;
#endif

  if(localP[i].Mass > 0)
    {
#ifdef MESHRELAX
      localSphP[i].Utherm = localSphP[i].Energy / localP[i].Mass;
#else
      localSphP[i].Utherm = (localSphP[i].Energy / localP[i].Mass - 0.5 * (localP[i].Vel[0] * localP[i].Vel[0] +
                                                                           localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2])) / (pvd->atime * pvd->atime);
#endif

#ifdef FLD
      MyFloat mu = 2.33;
      MyFloat u = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
      MyFloat temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;
      localSphP[i].Temperature = temp;
#endif

#if defined (VS_TURB) || defined (AB_TURB)
      localSphP[i].Utherm -= get_turb_pot(localP[i].Pos[0], localP[i].Pos[1], localP[i].Pos[2]);
#endif

#ifdef MHD
      localSphP[i].Utherm -= 0.5 * (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) / localSphP[i].Density / pvd->atime;
#endif

#ifdef MHD_THERMAL_ENERGY_SWITCH
      if(localSphP[i].Utherm < 0.)
        {
          EgyInjection -= localSphP[i].Energy;
          localSphP[i].Energy -= localSphP[i].Utherm * localP[i].Mass * pvd->atime * pvd->atime;
          localSphP[i].Utherm = localSphP[i].Etherm / localP[i].Mass / (pvd->atime * pvd->atime);
          localSphP[i].Energy += localSphP[i].Utherm * localP[i].Mass * pvd->atime * pvd->atime;
          EgyInjection += localSphP[i].Energy;
        }
#endif

#ifdef ACTIVE_CELL_SPIN
      {
        int j, k;
        double Erot = 0;
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            Erot += localSphP[i].Omega[j] * localSphP[i].Omega[k] * localSphP[i].MomentOfInertia[j][k] / localSphP[i].Volume;
        Erot *= 0.5;

        /* still undecided on including this term -> more testing required */
        //localSphP[i].Utherm -= Erot;
      }
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
      double Anew;
      int flag = 1;

#ifdef ENTROPY_MACH_THRESHOLD
      if(localSphP[i].MaxMach > ENTROPY_MACH_THRESHOLD)
        flag = 0;
#else
      terminate("not implemented any more");
#endif

      localSphP[i].A = exp(localSphP[i].Entropy / localP[i].Mass);

#ifdef TGCHEM
      Anew = (localSphP[i].Gamma - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
      Anew = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif

#ifndef MEASURE_DISSIPATION_RATE
      if(flag != 0 || Anew < 0.25 * localSphP[i].A)     /* here we keep the entropy, and reinitialize the thermal energy */
        {
#ifdef TGCHEM
          localSphP[i].Utherm = localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / (localSphP[i].Gamma - 1);
#else
          localSphP[i].Utherm = localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
          localSphP[i].Energy =
            localP[i].Mass * pvd->atime * pvd->atime * localSphP[i].Utherm +
            0.5 * localP[i].Mass * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);

          pvd->count_keep_entropy++;
        }
      else                      /* here the entropy is reinitialized, the thermal energy is kept */
#endif /* end of not MEASURE_DISSIPATION_RATE */
        {
#ifdef MEASURE_DISSIPATION_RATE
          if(dt)
            localSphP[i].DuDt = (localSphP[i].Utherm - localSphP[i].A * pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1) / GAMMA_MINUS1) / dt;
          else
            localSphP[i].DuDt = 0;
#endif
          localSphP[i].A = Anew;

          localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;

          pvd->count_update_entropy++;

          if(!gsl_finite(localSphP[i].A))
            {
              printf("i=%d BUMMER: %g %g egy=%g u=%g dens=%g m=%g\n", i, localSphP[i].A, localSphP[i].Entropy, localSphP[i].Energy, localSphP[i].Utherm, localSphP[i].Density, localP[i].Mass);
              terminate("stop");
            }
        }
#endif /* end of USE_ENTROPY_FOR_COLD_FLOWS */

      ulimit = All.MinEgySpec;
#ifdef JEANS_UTHERM_LIMIT
      //double celldim = 2.0 * get_cell_radius(i);
      //double jeans_length = JEANS_UTHERM_LIMIT*celldim;
      //double ulimit_jeans = All.G * localSphP[i].Density*jeans_length*jeans_length/M_PI/GAMMA/GAMMA_MINUS1;
      //ulimit = fmax(ulimit,ulimit_jeans);
      double max_rho = 4.0 * All.ReferenceGasPartMass / All.MinVolume;
      ulimit *= pow(localSphP[i].Density / max_rho, 2.0 / 3.0);
#endif


      if(localSphP[i].Utherm < ulimit)
        {
          EgyInjection -= localSphP[i].Energy;

          localSphP[i].Utherm = ulimit;

#ifdef MESHRELAX
          localSphP[i].Energy = localP[i].Mass * localSphP[i].Utherm;
#else
          localSphP[i].Energy =
            pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm +
            0.5 * localP[i].Mass * (localP[i].Vel[0] * localP[i].Vel[0] + localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#endif

#ifdef MHD
          localSphP[i].Energy += 0.5 * (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) * localSphP[i].Volume * pvd->atime;
#endif

#ifdef ACTIVE_CELL_SPIN
          {
            int j, k;
            double Erot = 0;
            for(j = 0; j < 3; j++)
              for(k = 0; k < 3; k++)
                Erot += localSphP[i].Omega[j] * localSphP[i].Omega[k] * localSphP[i].MomentOfInertia[j][k] * localSphP[i].Density;
            Erot *= 0.5;

            /* still undecided on including this term -> more testing required */
            //localSphP[i].Energy += Erot;
          }
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS

#ifdef TGCHEM
          localSphP[i].A = (localSphP[i].Gamma - 1) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
          localSphP[i].A = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif
          localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;

          if(!gsl_finite(localSphP[i].A))
            {
              printf("i=%d ouch: %g %g %g\n", i, localSphP[i].A, localSphP[i].Entropy, localP[i].Mass);
              terminate("stop");
            }
#endif
          EgyInjection += localSphP[i].Energy;
        }
    }
  else
    localSphP[i].Utherm = 0;

  if(localSphP[i].Density < All.LimitUBelowThisDensity && localSphP[i].Utherm > All.LimitUBelowCertainDensityToThisValue)
    {
      localSphP[i].Utherm = All.LimitUBelowCertainDensityToThisValue;
      localSphP[i].Energy =
        pvd->atime * pvd->atime * localP[i].Mass * localSphP[i].Utherm + 0.5 * localP[i].Mass * (localP[i].Vel[0] * localP[i].Vel[0] +
                                                                                                 localP[i].Vel[1] * localP[i].Vel[1] + localP[i].Vel[2] * localP[i].Vel[2]);
#ifdef MHD
      localSphP[i].Energy += 0.5 * (localSphP[i].B[0] * localSphP[i].B[0] + localSphP[i].B[1] * localSphP[i].B[1] + localSphP[i].B[2] * localSphP[i].B[2]) * localSphP[i].Volume * pvd->atime;
#endif

#ifdef USE_ENTROPY_FOR_COLD_FLOWS
#ifdef TGCHEM
      localSphP[i].A = (localSphP[i].Gamma - 1.0) * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#else
      localSphP[i].A = GAMMA_MINUS1 * localSphP[i].Utherm / pow(localSphP[i].Density * pvd->a3inv, GAMMA_MINUS1);
#endif
      localSphP[i].Entropy = log(localSphP[i].A) * localP[i].Mass;
#endif
    }

  if(localSphP[i].Utherm < 0)
    {
      printf("negative utherm %g\n", localSphP[i].Utherm);
      terminate("stop");
    }

#ifdef ENTROPY_MACH_THRESHOLD
#ifdef OUTPUT_MACHNUM
  localSphP[i].MaxMachNumber = localSphP[i].MaxMach;
#endif
  localSphP[i].MaxMach = 0;     /* reset */
#endif

#ifdef TGCHEM
  if(dt)
    localSphP[i].HydroHeatRate = (localSphP[i].Utherm - uold) / dt;
  else
    localSphP[i].HydroHeatRate = 0.;


  //if(localP[i].ID == 1068839580)
  //printf("unew = %d %g %g %g %g\n", localP[i].ID, uold, localSphP[i].Utherm, dt, localSphP[i].HydroHeatRate);
#endif

#endif /* end of not ISOTHERM_EQS */

  //temperature calculation

#ifdef RADCOOL_HOTHALO
#ifdef COOLING
  double meanweight = 4. / (3 * HYDROGEN_MASSFRAC + 1 + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
#else
  double meanweight = 0.5882352941176471;       //fully ionized
#endif
  localSphP[i].Temperature = localSphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g * GAMMA_MINUS1 * meanweight * PROTONMASS / BOLTZMANN;
#endif
}


double get_sound_speed(int p)
{
  double csnd;

#ifdef ISOTHERM_EQS
  csnd = All.IsoSoundSpeed;
#else
#ifdef LOCALLY_ISOTHERM_DISK
  csnd = get_isotherm_disk_sound_speed(p);
#else

  double gamma;
#ifdef TGCHEM
  gamma = SphP[p].Gamma;
#else
#ifdef VARIABLE_GAMMA
  gamma = SphP[p].GammaC;
#else /* This corresponds to ordinary gas dynamics, not ISOTHERM   */
  gamma = GAMMA;
#endif /* VARIABLE_GAMMA */
#endif /* TGCHEM */

  if(SphP[p].Density > 0)
    csnd = sqrt(gamma * SphP[p].Pressure / SphP[p].Density);
  else
    csnd = 0;

#endif /* LOCALLY_ISOTHERM_DISK */
#endif /* ISOTHERM_EQS */

#ifdef MHD
  /* for MHD, this is an upper bound to the signal velocity
     to do it more precisely, the magnet field in normal direction to the interfaces
     has to be taken into account */
  double Bsqr = SphP[p].B[0] * SphP[p].B[0] + SphP[p].B[1] * SphP[p].B[1] + SphP[p].B[2] * SphP[p].B[2];
  if(All.ComovingIntegrationOn)
    Bsqr /= All.Time;
  csnd = sqrt(csnd * csnd + Bsqr / SphP[p].Density);
#endif

#ifdef COSMIC_RAYS
  csnd = sqrt(csnd * csnd + All.GammaCR * SphP[p].CR_Pressure / SphP[p].Density);
#endif

#ifdef FLD
  csnd = sqrt(csnd * csnd + (4. / 9.) * SphP[p].n_gamma * (1. - exp(-SphP[p].Kappa_R * amr_length[Mesh.DP[p].level])) / SphP[p].Density);
#endif

  return csnd;
}


#ifdef LOCALLY_ISOTHERM_DISK
double get_isotherm_disk_sound_speed(int p)
{
  double dx, dy, r;
  double csnd = 0;
  dx = P[p].Pos[0] - boxHalf_X;
  dy = P[p].Pos[1] - boxHalf_Y;

  r = sqrt(dx * dx + dy * dy);

  if(r <= All.inner_radius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / All.inner_radius);
  if(r > All.inner_radius && r < All.outer_radius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / r);
  if(r >= All.inner_radius)
    csnd = All.AspectRatio * sqrt(All.G * All.CentralMass / All.outer_radius);

  return csnd;
}

int get_isotherm_disk_flag(int p)
{
  double dx, dy, r;
  dx = P[p].Pos[0] - boxHalf_X;
  dy = P[p].Pos[1] - boxHalf_Y;

  r = sqrt(dx * dx + dy * dy);

  if(r > All.inner_radius && r < All.outer_radius)
    return 1;
  else
    return 0;
}
#endif
