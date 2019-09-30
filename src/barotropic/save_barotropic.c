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

#ifdef BARO_CONSTANT_GAMMA_EOS
void calc_initial_adiabat( void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {
      if(P[i].Type == 0)
        SphP[i].InitialAdiabat = (GAMMA-1) * pow(SphP[i].Density, 1 - BARO_GAMMA) * SphP[i].Utherm;
    }

}
#endif


/** \brief Apply the new barotropic EOS to all gas cells.
 * 
 */
void apply_barotropic_eos(void)		/* normal cooling routine when star formation is disabled */
{
  int i;
  double minus_two_thirds = -2.0/3.0;
  double four_thirds = 4.0/3.0;
  double Bhattal_rho_0 = 1e-21 / All.UnitDensity_in_cgs;
  double Bhattal_rho_1 = 1e-13 / All.UnitDensity_in_cgs;
  double iso_rho_2 = 1e-12 / All.UnitDensity_in_cgs;
  
  double c_0 = 0.44721 * 1e5 / All.UnitVelocity_in_cm_per_s; /* 0.45 km/s eqiv to 100 K  */
  double c_s = 0.2 * 1e5 / All.UnitVelocity_in_cm_per_s;     /* 0.2 km/s eqiv to 10 K  */
  double c_0_2, c_s_2;

  c_0_2 = c_0 * c_0;
  c_s_2 = c_s * c_s;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  if(P[i].Mass == 0 && P[i].ID == 0)
	    continue;		/* skip cells that have been swallowed or eliminated */

#ifdef BARO_CONSTANT_GAMMA_EOS          
          SphP[i].Utherm = SphP[i].InitialAdiabat * pow(SphP[i].Density, BARO_GAMMA-1) / (GAMMA-1);
#else
          /* Using the Bhattal '98 parameterisation */
          if(SphP[i].Density <= Bhattal_rho_0)
            SphP[i].Utherm = c_0_2 * c_0_2 / (GAMMA-1);
          else
            {
              if (SphP[i].Density <= iso_rho_2 )
                {
                  SphP[i].Utherm = ((c_0_2 - c_s_2) * pow(SphP[i].Density/Bhattal_rho_0, minus_two_thirds) + c_s_2);
                  SphP[i].Utherm *= sqrt(1 + pow(SphP[i].Density/Bhattal_rho_1, four_thirds)) / (GAMMA-1);
                }
              else
                {  /* We hack an isothermal regime so that sinks can form */
                  SphP[i].Utherm = ((c_0_2 - c_s_2) * pow(iso_rho_2/Bhattal_rho_0, minus_two_thirds) + c_s_2);
                  SphP[i].Utherm *= sqrt(1 + pow(iso_rho_2/Bhattal_rho_1, four_thirds)) / (GAMMA-1);
                }
            }
#endif

          SphP[i].Energy = All.cf_atime * All.cf_atime * SphP[i].Utherm * P[i].Mass;
#ifdef USE_ENTROPY_FOR_COLD_FLOWS
          SphP[i].A = (GAMMA - 1.0) * SphP[i].Utherm / pow(SphP[i].Density * All.cf_a3inv, GAMMA - 1);
          SphP[i].Entropy = log(SphP[i].A) * P[i].Mass;
#endif
          set_pressure_of_cell(i);
	}
    }

}

#endif
