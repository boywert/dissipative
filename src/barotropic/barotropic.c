/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/barotropic/barotropic.c
 * \date        11/2014
 * \author      Paul C. Clark
 * \brief       Provides the kind of barotropic cooling used in SF collapse simulations
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



/** \file barotropic.c
 *  \brief The barotropic EOS used in SF collapse simulations
 */

#ifdef BAROTROPIC
/** \brief Apply the new barotropic EOS to all gas cells.
 * 
 */
void apply_barotropic_eos(void)		/* normal cooling routine when star formation is disabled */
{
  int i, idx;
  double minus_two_thirds = -2.0/3.0;
  double four_thirds = 1.4;
  double Bhattal_rho_0 = 1e-21 / All.UnitDensity_in_cgs;
  double Bhattal_rho_1 = 1e-14 / All.UnitDensity_in_cgs;
  double iso_rho_2 = 1e-9 / All.UnitDensity_in_cgs;
  
  double c_0 = 0.44721 * 1e5 / All.UnitVelocity_in_cm_per_s; /* 0.45 km/s eqiv to 100 K  */
  double c_s = 0.2 * 1e5 / All.UnitVelocity_in_cm_per_s;     /* 0.2 km/s eqiv to 10 K  */
  double c_0_2, c_s_2;

  double density;
  double t1, t2;

  c_0_2 = c_0 * c_0;
  c_s_2 = c_s * c_s;

  int icount_baro = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0) 
	{
	  if((P[i].Mass == 0 && P[i].ID == 0))
	    continue;		/* skip cells that have been swallowed or eliminated */
          
#ifdef BARO_CONSTANT_GAMMA_EOS          
          SphP[i].Utherm = SphP[i].InitialAdiabat * pow(SphP[i].Density, BARO_GAMMA-1) / (GAMMA-1);
#else
#ifdef BHATTAL98
          if (icount_baro == 0)
             mpi_printf("BAROTROPIC: Utherm before EOS %g Density %g A %g Entropy %g\n", SphP[i].Utherm, SphP[i].Density, SphP[i].A, SphP[i].Entropy);
          /* Using the Bhattal '98 parameterisation */
          if(SphP[i].Density <= Bhattal_rho_0)
            SphP[i].Utherm = c_0_2;
          else
            {
              /* Allow an isothermal regime to help sinks form */
              if (SphP[i].Density < iso_rho_2)
                density = SphP[i].Density;
              else
                density = iso_rho_2;
              t1 = ((c_0_2 - c_s_2) * pow(density/Bhattal_rho_0, minus_two_thirds) + c_s_2);
              t2 = ((c_0_2 - c_s_2) * pow(Bhattal_rho_1/Bhattal_rho_0, minus_two_thirds) + c_s_2);
              t2 = t2 * pow(density/Bhattal_rho_1, four_thirds - 1);
              SphP[i].Utherm = sqrt(t1*t1 + t2*t2);
            }
#else
           /* using a standard smooth isothermal / adiabatic switch 
           */
           t1 = 10.0 * (1.0 + pow(SphP[i].Density/Bhattal_rho_1, 0.4));
           SphP[i].Utherm = t1 * BOLTZMANN / PROTONMASS / 2.33 / (GAMMA - 1) / (All.UnitEnergy_in_cgs / All.UnitMass_in_g); 
#endif
#endif
          SphP[i].Energy = All.cf_atime * All.cf_atime * SphP[i].Utherm * P[i].Mass;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          SphP[i].A = (GAMMA - 1.0) * SphP[i].Utherm / pow(SphP[i].Density * All.cf_a3inv, GAMMA - 1);
          SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif
          if (icount_baro == 0)
            mpi_printf("BAROTROPIC: Utherm after EOS %g A %g Entropy %g \n", SphP[i].Utherm, SphP[i].A, SphP[i].Entropy);         
 
          set_pressure_of_cell(i);
          icount_baro++;
	}
    }

}

#endif
