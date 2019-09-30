/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/cooling_dust.c
 * \date        MM/YYYY
 * \author      Mark Vogelsberger & Ryan McKinnon
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
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "../allvars.h"
#include "../proto.h"

#ifdef GFM_DUST
#ifdef GFM_DUST_COOLING

double get_CoolingDustRate(double logT, double dust_to_gas_ratio, double mu, char dust_cool)
{
/* MRN grain distribution */
#ifdef GFM_DUST_MRN
  if (dust_cool == 1)
    {
      if ((logT > 5) && (logT < 9.0))
        return (dust_to_gas_ratio/0.0075) * pow(10.0, 0.01575072*logT*logT*logT*logT - 0.41930701*logT*logT*logT + 3.87861885*logT*logT - 13.66242274*logT - 9.09918663);
      else
        return 0.0;
    }
  else 
    return 0.0;
#endif

/* single grainsize */
#ifndef GFM_DUST_MRN
    if (dust_cool == 1)
    {
    double T = pow(10.0, logT);
    double a = All.DustSingleGrainSize; // grain size in mu
    double x = 2.71e8 * pow(a, (2.0/3.0)) / T;
    double H_ne;

    if (x > 4.5)
        H_ne = 5.38e-18 * a*a * pow(T,1.5);
    else if  ((x < 4.5) && (x > 1.5))
        H_ne = 3.37e-13 * pow(a,2.41) * pow(T,0.88);
    else
        H_ne = 6.48e-6 * a*a*a;

    double m_p = 1.67e-24; // g
    double a_cm = a/1.0e6*100;// cm
    double rho_grain = 3.0; // g / cm^3
    double fac = mu * m_p / (4.0/3.0 * 3.141592 * a_cm*a_cm*a_cm * rho_grain) * dust_to_gas_ratio;  
    return H_ne * fac;
    }
    else
      return 0.0;
#endif

}

void get_CoolingDustState(int i, double *dust_to_gas_ratio, char *dust_cool)
{
  *dust_to_gas_ratio = 0.0;
  for(int chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
  {
    for(int k_elem = 0; k_elem < GFM_N_CHEM_ELEMENTS; k_elem++)
    {
      *dust_to_gas_ratio += SphP[i].MetalsDustFraction[chan][k_elem];
    }
  }
  if (SphP[i].Sfr == 0.0)
    *dust_cool = 1;
  else
    *dust_cool = 0;
}
#endif
#endif
