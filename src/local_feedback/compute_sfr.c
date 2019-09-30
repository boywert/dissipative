/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        
 * \date        10/2014
 * \author      Christine Simpson
 * \brief
 * \details     Compute sfr either from KS relation or from the local dynamical time
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include "../allvars.h"
#include "../proto.h"

#if defined(USE_SFR) && defined(LOCAL_FEEDBACK) && !defined(EXTERNALSHEARBOX_KSRATE_RANDOM)

void compute_sfr(void)
{
  CPU_Step[CPU_MISC] += measure_time();
  
  int i, idx;
  double sfr, dt, deltaM;
  double sum_sfr = 0.0;
  double sum_sm = 0.0;
  double sum_mass_stars = 0.0;

#ifdef EXTERNALSHEARBOX
  double dz, sech, rho0, rho, bscale, bscale0, Sigma0, SigmaNow;
  double fg = All.ShearBoxFg;
  double mu = All.ShearBoxMu;
  Sigma0 = All.ShearBoxSigma0;
  bscale0 = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);
  Sigma0 *= (1.9891e33 / All.UnitMass_in_g) * (All.UnitLength_in_cm / 3.0857e18) * (All.UnitLength_in_cm / 3.0857e18);
  bscale0 *= (3.0857e18 / All.UnitLength_in_cm);
#ifdef EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM
  bscale = All.ShearBoxB;
  SigmaNow = All.ShearBoxSigmaNow;
#else
  SigmaNow = All.ShearBoxSigma0;
  bscale = 61.0 * (fg / 0.1 / mu) / (Sigma0 / 10.0);
  // convert Sigma0 and bscale to code units from Msun/pc^2 and pc
  SigmaNow *= (1.9891e33 / All.UnitMass_in_g) * (All.UnitLength_in_cm / 3.0857e18) * (All.UnitLength_in_cm / 3.0857e18);
  bscale *= (3.0857e18 / All.UnitLength_in_cm);
#endif
#endif

#if !defined(EXTERNALSHEARBOX_KSRATE) || defined(EXTERNALSHEARBOX_MIXED_INJECTION)
  double SFEff = All.LocalFeedbackSFEff;
  double SFDenThresh = All.LocalFeedbackSFDenThresh * PROTONMASS * All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
  double density_threshold = SFDenThresh;
#endif

#ifdef EXTERNALSHEARBOX_KSRATE /* compute SFR from surface density */
  double sigma_star = 2.5e-4 * pow(SigmaNow, 1.4) * (All.UnitTime_in_s / 3.15569e7) * (1.9891e33 / All.UnitMass_in_g) * (All.UnitLength_in_cm / 3.0857e21) * (All.UnitLength_in_cm / 3.0857e21);        /*KS relation in code units */
  mpi_printf("CMS_FEEDBACK: compute_sfr bscale0=%g Sigma0=%g bscale=%g SigmaNow=%g sigma_star=%g\n", bscale0, Sigma0, bscale, SigmaNow, sigma_star);
#endif

  double prob_factor = 1.0;
#ifdef EXTERNALSHEARBOX_MIXED_INJECTION
  prob_factor = EXTERNALSHEARBOX_MIXED_INJECTION;
#endif

  /* loop over all active cells */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* just consider gas cells */
      if(P[i].Type == 0)
        {
          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;

          SphP[i].Sfr = 0.0;

#ifdef EXTERNALSHEARBOX
          dz = P[i].Pos[2] - boxHalf_Z;
          sech = 1.0 / cosh(dz / bscale);
          rho0 = Sigma0 * sech * sech / 2.0 / bscale0;
          rho = SigmaNow * sech * sech / 2.0 / bscale;
#endif
          /* compute sfr */

#if !defined(EXTERNALSHEARBOX_KSRATE) || defined(EXTERNALSHEARBOX_MIXED_INJECTION) /* compute sfr from dynamical time */

          //if(SphP[i].DivVel >= 0.0)
          //  continue;

          if(SphP[i].Density < density_threshold)
            continue;

          double tdyn, tot_rho;

          tot_rho = SphP[i].Density;
#ifdef EXTERNALSHEARBOX
          tot_rho += (1.0 / fg - 1.0) * rho0;
#endif
          tdyn = sqrt(3.0 * M_PI / 32.0 / All.G / tot_rho);
          sfr = SFEff * P[i].Mass / tdyn;
	  sfr *= prob_factor;

#else /* compute sfr from ks law for shearing box ic's */

#ifdef EXTERNALSHEARBOX_KSRATE_IC_MASS  /* Use cell mass assuming isothermal profile */

#ifdef EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM     /* Use isothermal profile of current gas distribution */
          sfr = sigma_star * rho * SphP[i].Volume / SigmaNow;
#else /* Use isothermal profile from initial conditions */
          sfr = sigma_star * rho0 * SphP[i].Volume / Sigma0;
#endif

#else /* Use current cell mass */

          sfr = sigma_star * P[i].Mass / SigmaNow;

#endif /* end EXTERNALSHEARBOX_KSRATE_IC_MASS */

#endif
          SphP[i].Sfr = sfr;
	  sum_sfr += sfr;
        }
    } /* end loop over active cells */

#if defined(EXTERNALSHEARBOX_KSRATE) & defined(EXTERNALSHEARBOX_MIXED_INJECTION)

  // sum sfr
  double totpeaksfrrate;
  MPI_Allreduce(&sum_sfr, &totpeaksfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  double scale_factor;
  scale_factor = totpeaksfrrate/(sigma_star*boxSize_X*boxSize_Y*prob_factor);
  mpi_printf("CMS_FEEDBACK: scale_factor: %e sigma_star: %e totpeaksfrrate: %e\n",scale_factor,sigma_star,totpeaksfrrate);
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
	continue;
      
      /* just consider gas cells */
      if(P[i].Type == 0)
	{
	  /* skip cells that have been swallowed or eliminated */
	  if(P[i].Mass == 0 && P[i].ID == 0)
	    continue;
	  
	  SphP[i].Sfr /= scale_factor;
	}
      }	  
#endif
}
#endif


#if defined(EXTERNALSHEARBOX_KSRATE_UPDATE_PARAM)

void update_shearbox_param(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  mpi_printf("CMS_FEEDBACK: Starting Update Shearbox Param\n");

  int i, idx;
  double arctanh_half = 0.54930614433405;       // arctanh(0.5)
  double sum_mass;
  double sum_mass_local = 0.0;

  int count = 0;
  /* Compute total gas mass in the box */
  /* loop over all active cells */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* just consider gas cells */
      if(P[i].Type == 0)
        {
          /* skip cells that have been swallowed or eliminated */
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;
          count++;
          sum_mass_local += P[i].Mass;
        }
    }

  MPI_Allreduce(&sum_mass_local, &sum_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  int sum_gas_particles;
  MPI_Allreduce(&count, &sum_gas_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  double box_area = boxSize * boxSize;
  All.ShearBoxSigmaNow = sum_mass / box_area;

  double h;
  double total_mass = 1.0 * sum_mass;
  double half_total_mass = 0.5 * sum_mass;
  double dh = 0.5 * boxSize_Z;
  double hsum = 0.5 * boxSize_Z;
  int niter = 0;
  int new_sum_gas_particles = 10 * sum_gas_particles;
  double all_gas_particles = sum_gas_particles * 1.0;

  /* Compute height that contains half the mass; this is related to the scale height */
  while(fabs(sum_mass - half_total_mass) / half_total_mass > 0.01)
    {
      if(fabs(sum_gas_particles - new_sum_gas_particles) < 5)
        break;                  //if subsequent bins are only a few particles different in size, break

      sum_gas_particles = new_sum_gas_particles * 1;

      if(sum_mass > half_total_mass)
        {
          dh *= 0.5;
          hsum -= dh;
        }
      else
        {
          hsum += dh * 0.5;
        }

      if(niter > 100)
        {
          terminate("CMS_FEEDBACK: Stratified box scale height did not converge: hsum %g; sum_mass %g; half_total_mass %g\n", hsum, sum_mass, half_total_mass);
        }

      if(hsum > 0.5 * boxSize_Z)
        {
          terminate("CMS_FEEDBACK: half_mass height cannot exceed box height: niter %d; hsum %g; sum_mass %g; half_total_mass %g\n", niter, hsum, sum_mass, half_total_mass);
        }

      double sum_mass_local = 0.0;
      int count = 0;

      /* loop over all active cells */
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;

          /* just consider gas cells */
          if(P[i].Type == 0)
            {
              /* skip cells that have been swallowed or eliminated */
              if(P[i].Mass == 0 && P[i].ID == 0)
                continue;

              h = P[i].Pos[2] - boxHalf_Z;
              if(fabs(h) < hsum)
                {
                  sum_mass_local += P[i].Mass;
                  count++;
                }
            }
        }
      
      MPI_Allreduce(&sum_mass_local, &sum_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&count, &new_sum_gas_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      niter++;
    }

  All.ShearBoxB = hsum / arctanh_half;

  mpi_printf("CMS_FEEDBACK: Sigma0 = %g; bscale = %g; niter = %d\n", All.ShearBoxSigma0, All.ShearBoxB, niter);

}
#endif
