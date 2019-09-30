/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/atomic_dm/cooling_atomic_dm.c
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

#include "../allvars.h"
#include "../proto.h"

/** \file cooling_atomic_dm.c
 */

#ifdef ATOMIC_DM


static ADM_GasState ADM_gs;           /**< gas state */
static ADM_DoCoolData ADM_DoCool;     /**< cooling data */


static double ADMProtonMassInCgs, ADMElectronMassInCgs;



/** \brief Compute the new internal energy per unit mass.
 * 
 *   The function solves for the new internal energy per unit mass of the gas by integrating the equation
 *   for the internal energy with an implicit Euler scheme. The root of resulting non linear equation,
 *   which gives tnew internal energy, is found with the bisection method.
 *   Arguments are passed in code units.
 *   
 *   \param u_old the initial (before cooling is applied) internal energy per unit mass of the gas cell
 *   \param rho   the proper density of the gas cell
 *   \param dt    the duration of the time step
 *   \param ne_guess electron number density relative to hydrogen number density (for molecular weight computation)
 *   \return the new internal energy per unit mass of the gas cell
 */
double ADM_DoCooling(double u_old, double rho, double dt, double *ne_guess)
{
  double u, du;
  double u_lower, u_upper;
  double ratefact;
  double LambdaNet;

  int iter = 0;

  ADM_DoCool.u_old_input = u_old;
  ADM_DoCool.rho_input = rho;
  ADM_DoCool.dt_input = dt;
  ADM_DoCool.ne_guess_input = *ne_guess;

  if(!gsl_finite(u_old))
    terminate("invalid input: u_old=%g\n", u_old);

  if(u_old < 0 || rho < 0)
    terminate("invalid input: task=%d u_old=%g  rho=%g  dt=%g  All.MinEgySpec=%g\n", ThisTask, u_old, rho, dt, All.MinEgySpec);

  rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;    /* convert to physical cgs units */
  u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
  dt *= All.UnitTime_in_s / All.HubbleParam;

  ADM_gs.nHcgs = rho / ADMProtonMassInCgs;      /* hydrogen number dens in cgs units */
  ratefact = ADM_gs.nHcgs * ADM_gs.nHcgs / rho;

  u = u_old;
  u_lower = u;
  u_upper = u;

  LambdaNet = ADM_CoolingRateFromU(u, rho, ne_guess);

  /* solve for new u; implicit Euler for du/dt = Lambda(u) */
  /* bracketing */
  if(u - u_old - ratefact * LambdaNet * dt < 0) /* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while(u_upper - u_old - ratefact * ADM_CoolingRateFromU(u_upper, rho, ne_guess) * dt < 0)
        {
          u_upper *= 1.1;
          u_lower *= 1.1;
        }
    }

  if(u - u_old - ratefact * LambdaNet * dt > 0) /* cooling */
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while(u_lower - u_old - ratefact * ADM_CoolingRateFromU(u_lower, rho, ne_guess) * dt > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;
        }
    }

  /* bisection algorithm */
  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = ADM_CoolingRateFromU(u, rho, ne_guess);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
        {
          u_upper = u;
        }
      else
        {
          u_lower = u;
        }

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
        printf("u= %g\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  if(iter >= MAXITER)
    terminate("Failed to converge in ADM_DoCooling(): u_old_input=%g  rho_input= %g  cool.dt_input= %g   ne_guess_input= %g\n",
              ADM_DoCool.u_old_input, ADM_DoCool.rho_input, ADM_DoCool.dt_input, ADM_DoCool.ne_guess_input);

  u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs;        /* to internal units */

  return u;
}



/** \brief Get cooling rate from gas internal energy.
 * 
 *  This function first computes the self-consistent temperature
 *  and abundance ratios, and then it calculates
 *  (heating rate-cooling rate)/n_h^2 in cgs units
 * 
 *  \param u   gas internal energy per unit mass
 *  \param rho gas density
 *  \param ne_guess electron number density relative to hydrogen number density
 */
double ADM_CoolingRateFromU(double u, double rho, double *ne_guess)
{
  double Tplasma, mu;
  double Lambda, LambdaFF, LambdaCmptn;
  double Heat;
  double gff = 1;
  double TdarkCMB, redshift;

  /* set ionization (all in units of nH) */
  ADM_gs.nH0 = 0;
  ADM_gs.nHp = 1.0;
  ADM_gs.ne = ADM_gs.nHp;
  *ne_guess = ADM_gs.ne;

  /* nH in cgs */
  ADM_gs.nHcgs = rho / ADMProtonMassInCgs;

  /* calculate mu */
  mu = 1. / (1 + ADM_gs.ne);

  /* calculate Tplasma */
  Tplasma = GAMMA_MINUS1 / BOLTZMANN * u * ADMProtonMassInCgs * mu;

  /* dark photon temperature */
  TdarkCMB = 1.35;              //keep this constant for non-comoving runs FIXME

  if(All.ComovingIntegrationOn)
    {
      redshift = 1 / All.Time - 1;
      TdarkCMB = (1 + redshift) * 1.35;
    }
  else
    TdarkCMB = 1.35;

  /* cooling rates */

  /* Bremsstrahlung */
  LambdaFF = gff * ADM_gs.nHp * ADM_gs.ne * 3.66e-27 * pow(All.ADMFineStructureConstant / 1e-2, 3) * pow(511.0 / All.ADMElectronMassInkeV, 1.5) * pow(Tplasma, 0.5);

  /* Compton cooling off the dark CMB */
  LambdaCmptn = ADM_gs.ne * 1.915e-37 * (Tplasma - TdarkCMB) * pow(511 / All.ADMElectronMassInkeV, 3.0) * pow(All.ADMFineStructureConstant / 1e-2, 2) * pow(TdarkCMB, 4);

  Lambda = LambdaFF + LambdaCmptn;

  Heat = 0;

  return (Heat - Lambda);
}



/** \brief Initialize the cooling module.
 * 
 *   This function initializes the cooling module. In particular,
 *   it allocates the memory for the cooling rate and ionization tables
 *   and initializes them.
 */
void ADM_InitCool(void)
{
  ADMProtonMassInCgs = All.ADMProtonMassInkeV * keV_to_gramm;
  ADMElectronMassInCgs = All.ADMElectronMassInkeV * keV_to_gramm;

  All.Time = All.TimeBegin;
  set_cosmo_factors_for_current_time();
}

/** \brief Apply the isochoric cooling to all the active gas cells.
 * 
 */
void ADM_cooling(void)          /* normal cooling routine when star formation is disabled */
{
  int idx, i;
  mpi_printf("ATOMIC_DM: ADM Cooling (sync-point %d).\n", All.NumCurrentTiStep);
  CPU_Step[CPU_MISC] += measure_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      ADM_cool_cell(i);
    }

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/** \brief Apply the isochoric cooling to a given gas cell.
 * 
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models implemented, 
 *  the internal energy per unit mass, the total energy and the pressure of the cell are updated. 
 * 
 *  \param i index of the gas cell to which cooling is applied
 */
void ADM_cool_cell(int i)
{
  double dt, dtime, ne = 1;
  double unew, dens, dtcool;

  set_cosmo_factors_for_current_time();

  dens = SphP[i].Density;

  dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;
  dtcool = dtime;

  ne = SphP[i].Ne;              /* electron abundance (gives ionization state and mean molecular weight) */
  unew = ADM_DoCooling(dmax(All.MinEgySpec, SphP[i].Utherm), dens * All.cf_a3inv, dtcool, &ne);
  SphP[i].Ne = ne;

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

  set_pressure_of_cell(i);
}

#endif
