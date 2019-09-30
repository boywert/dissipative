/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays.c
 * \date        11/2013
 * \author      R. Pakmor & C. Pfrommer
 * \brief       A simple cosmic ray implementation
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "string.h"

#ifdef COSMIC_RAYS

void init_cosmic_rays(void)
{
  if(All.ComovingIntegrationOn && fabs(All.GammaCR - 4. / 3.) > 1e-3)
    terminate("Comoving integration with cosmic rays requires GammaCR = 4./3.");

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION

  if(All.AccelerationEfficiency == 0)
    All.AccelerationEfficiency = 0.1;   /* fiducial value of Caprioli D., Spitkovsky A., 2014, ApJ, 783, 91 */

  if(All.CriticalMachnumber == 0)
    All.CriticalMachnumber = 3.;
#endif
}

#ifdef COSMIC_RAYS_COOLING
/* adopting CR cooling of an equilibrium distribution; assumptions:
   alpha=2.2, q_min=0.1 (only relevant for Coulomb heating/cooling), ne = 1/cm^3 (and then scaled to different density):
   Gamma = Gamma_hadronic + Gamma_Coulomb, 
   Gamma_Coulomb  = 2.78e-16 / s, 
   Gamma_hadronic = 7.44e-16 / s ~ 0.5 * SIGMA_PP * CLIGHT / (1. - 0.5 * (1. - HYDROGEN_MASSFRAC))
   with SIGMA_PP  = 4.4e-26 as effective pp cross section in cm^2 for a CR population f(p) = C p^-alpha with alpha=2.2
   collisional heating rate due to Coulomb and hadronic interactions (5/6 of the latter is emitted in neutrinos and gama-rays):
   Gamma_heating = Gamma_Coulomb + Gamma_hadronic / 6 = 4.02e-16 / s */

void do_cosmic_ray_cooling(void)
{
  int idx, i;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;    /* timestep of the cell in Arepo units */
      double rho = SphP[i].Density * All.UnitDensity_in_cgs * All.cf_a3inv * All.HubbleParam * All.HubbleParam;         /* physical mass density in cgs */
      double Gamma_heating = 4.02e-16 * rho / PROTONMASS * (1. - 0.5 * (1. - HYDROGEN_MASSFRAC)) * All.UnitTime_in_s;   /* collisional heating rate in Arepo units */
      double Gamma = 1.022e-15 * rho / PROTONMASS * (1. - 0.5 * (1. - HYDROGEN_MASSFRAC)) * All.UnitTime_in_s;          /* total CR cooling rate in Arepo units */

      double Gamma_Coulomb  = 2.78e-16 * rho / PROTONMASS * (1. - 0.5 * (1. - HYDROGEN_MASSFRAC)) * All.UnitTime_in_s;
      double Gamma_Hadronic = 7.44e-16 * rho / PROTONMASS * (1. - 0.5 * (1. - HYDROGEN_MASSFRAC)) * All.UnitTime_in_s;
      All.TotalCREnergyCooledCoulomb  += SphP[i].CR_Energy * ( exp(-Gamma_Coulomb  * dt_cell) - 1. ) / All.cf_atime;
      All.TotalCREnergyCooledHadronic += SphP[i].CR_Energy * ( exp(-Gamma_Hadronic * dt_cell) - 1. ) / All.cf_atime;
      
      double CR_EnergyOld = SphP[i].CR_Energy;
      SphP[i].Energy += SphP[i].CR_Energy * (1. - exp(-Gamma_heating * dt_cell)) * All.cf_atime;
      SphP[i].CR_Energy *= exp(-Gamma * dt_cell);
      SphP[i].CR_Energy = dmax(SphP[i].CR_Energy, All.MinimumCREnergyDensity * SphP[i].Volume);
      
      double dCREnergyEff = CR_EnergyOld - SphP[i].CR_Energy;
      All.TotalCREnergyCooled -= dCREnergyEff;

      SphP[i].CR_SpecificEnergy = SphP[i].CR_Energy / P[i].Mass;
    }
  update_primitive_variables();
}
#endif

#ifdef COSMIC_RAYS_ALFVEN_COOLING
/* isotropic CR transport: Alfven cooling term is an upper limit to the (projected) anisotropic Alfven cooling term */
void do_cosmic_ray_Alfven_cooling(void)
{
  int idx, i;
  /* update gradients */
  exchange_primitive_variables();
  calculate_gradients();
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt_cell = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval / All.cf_hubble_a;    /* timestep of the cell in Arepo units */

      double vel_Alfven = sqrt((SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) / SphP[i].Density);
      double b = sqrt(SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]);
      double prod = (SphP[i].B[0] * SphP[i].Grad.dcrPressure[0] + SphP[i].B[1] * SphP[i].Grad.dcrPressure[1] + SphP[i].B[2] * SphP[i].Grad.dcrPressure[2]) / b;

      double dCR_Energy = dmax(0., dmin(dt_cell * vel_Alfven * fabs(prod) * SphP[i].Volume / sqrt(All.cf_atime), SphP[i].CR_Energy - All.MinimumCREnergyDensity * SphP[i].Volume));

      SphP[i].CR_Energy -= dCR_Energy;
      SphP[i].CR_SpecificEnergy = SphP[i].CR_Energy / P[i].Mass;
      SphP[i].Energy += dCR_Energy * All.cf_atime;
    }
  update_primitive_variables();
}
#endif

#ifdef COSMIC_RAYS_DIFFUSION_OLD
void do_cr_diffusion(void)
{
  double *diffusion_coeff, *cr_energy_density;
  double energy_sum, energy_sum_all;
  int i;
#ifdef COSMIC_RAYS_DIFFUSION_LIMITER
  double diffusion_coeff_max, diffusion_coeff_0, vel_Alfven, grad_pcr;
  double *diffusion_coeff_ratio;
  double diffusion_coeff_max_new, vel_Alfven_new, grad_pcr_new;
  double diffusion_coeff_ratio_max, diffusion_coeff_ratio_max_all;
#endif

  if(All.CR_Diffusion_Coefficient == 0)
    All.CR_Diffusion_Coefficient = 1e28;        /* fiducial value from Salem & Bryan */
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
    if(All.Time >= All.CR_Diffusion_Time + All.CR_Diffusion_TimeStep || 1)
      {
        mpi_printf("COSMIC_RAYS: Doing CR diffusion.\n");
        cr_energy_density = (double *) mymalloc("cr_energy_density", NumGas * sizeof(double));
        diffusion_coeff = (double *) mymalloc("cr_diffusion_coeff", NumGas * sizeof(double));

#ifdef COSMIC_RAYS_DIFFUSION_LIMITER
        /* diagnostics for change of kappa during diffusion time step */
        diffusion_coeff_ratio = (double *) mymalloc("cr_diffusion_coeff_ratio", NumGas * sizeof(double));
        /* update gradients */
        exchange_primitive_variables();
        calculate_gradients();
#endif

        energy_sum = 0;
        for(i = 0; i < NumGas; i++)
          if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
            {
              cr_energy_density[i] = SphP[i].CR_Energy / SphP[i].Volume;
#if !defined(COSMIC_RAYS_DIFFUSION_LIMITER)
              diffusion_coeff[i] = All.CR_Diffusion_Coefficient / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;
#else
              /* limiting the diffusion velocity to the Alfven speed using a smooth interpolating function */
              /* kappa = [1 - exp(-kappa_0/kappa_max)] * kappa_max, where kappa_0 = global constant variable and kappa_max = v_A P_cr / \nabla P_cr */
              vel_Alfven = sqrt((SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) / SphP[i].Density);
              grad_pcr = sqrt(SphP[i].Grad.dcrPressure[0] * SphP[i].Grad.dcrPressure[0] +
                              SphP[i].Grad.dcrPressure[1] * SphP[i].Grad.dcrPressure[1] + SphP[i].Grad.dcrPressure[2] * SphP[i].Grad.dcrPressure[2]);
              diffusion_coeff_0 = All.CR_Diffusion_Coefficient / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;

              if(grad_pcr > 0 && vel_Alfven > 0)
                {
                  diffusion_coeff_max = vel_Alfven * SphP[i].CR_Pressure / grad_pcr;
                  diffusion_coeff[i] = (1. - exp(-diffusion_coeff_0 / diffusion_coeff_max)) * diffusion_coeff_max;
                }
              else
                {
                  diffusion_coeff[i] = diffusion_coeff_0;
                }
#endif
              energy_sum += SphP[i].CR_Energy;
            }
        MPI_Reduce(&energy_sum, &energy_sum_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        mpi_printf("COSMIC_RAYS: Total CR Energy before diffusion %g erg.\n", energy_sum_all * All.UnitEnergy_in_cgs);

        double dt = (All.HighestActiveTimeBin ? (((integertime) 1) << All.HighestActiveTimeBin) : 0) * All.Timebase_interval;
        diffuse(cr_energy_density, diffusion_coeff, dt);

        energy_sum = 0;
        for(i = 0; i < NumGas; i++)
          if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
            {
              SphP[i].CR_SpecificEnergy = cr_energy_density[i] * SphP[i].Volume / P[i].Mass;
              SphP[i].CR_Energy = SphP[i].CR_SpecificEnergy * P[i].Mass;
              energy_sum += SphP[i].CR_Energy;
            }

#ifdef COSMIC_RAYS_DIFFUSION_LIMITER
        /* diagnostics for change of kappa during diffusion time step */
        update_primitive_variables();
        exchange_primitive_variables();
        calculate_gradients();
        diffusion_coeff_ratio_max = 0;

        for(i = 0; i < NumGas; i++)
          if(P[i].Type == 0 && P[i].ID != 0 && P[i].Mass > 0)
            {
              vel_Alfven_new = sqrt((SphP[i].B[0] * SphP[i].B[0] + SphP[i].B[1] * SphP[i].B[1] + SphP[i].B[2] * SphP[i].B[2]) / SphP[i].Density);
              grad_pcr_new = sqrt(SphP[i].Grad.dcrPressure[0] * SphP[i].Grad.dcrPressure[0] +
                                  SphP[i].Grad.dcrPressure[1] * SphP[i].Grad.dcrPressure[1] + SphP[i].Grad.dcrPressure[2] * SphP[i].Grad.dcrPressure[2]);
              diffusion_coeff_0 = All.CR_Diffusion_Coefficient / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;

              if(grad_pcr_new > 0 && vel_Alfven_new > 0)
                {
                  diffusion_coeff_max_new = vel_Alfven_new * SphP[i].CR_Pressure / grad_pcr_new;
                  diffusion_coeff_ratio[i] = (1. - exp(-diffusion_coeff_0 / diffusion_coeff_max_new)) * diffusion_coeff_max_new;
                }
              else
                {
                  diffusion_coeff_ratio[i] = diffusion_coeff_0;
                }
              diffusion_coeff_ratio[i] = fabs(diffusion_coeff_ratio[i] / diffusion_coeff[i] - 1.);
              if(diffusion_coeff_ratio[i] > diffusion_coeff_ratio_max)
                {
                  diffusion_coeff_ratio_max = diffusion_coeff_ratio[i];
                }
            }
        MPI_Reduce(&diffusion_coeff_ratio_max, &diffusion_coeff_ratio_max_all, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        mpi_printf("COSMIC_RAYS: max(kappa_new / kappa) = %g.\n", diffusion_coeff_ratio_max_all);
        //mpi_printf( "COSMIC_RAYS: max(kappa_new / kappa) = %g.\n", max_array(diffusion_coeff_ratio, NumGas));
        myfree(diffusion_coeff_ratio);
#endif
        MPI_Reduce(&energy_sum, &energy_sum_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        mpi_printf("COSMIC_RAYS: Total CR Energy after diffusion %g erg.\n", energy_sum_all * All.UnitEnergy_in_cgs);

        myfree(diffusion_coeff);
        myfree(cr_energy_density);

        All.CR_Diffusion_Time = All.Time;
      }
}
#endif

#endif
