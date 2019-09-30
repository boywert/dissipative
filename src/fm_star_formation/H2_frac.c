/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/H2_frac.c
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

#include <math.h>
#include "../allvars.h"

#ifdef COMPUTE_SFR_FROM_H2

#define SOLAR_METALLICITY 0.0127        /* the same as GFM */

/* computes the molecular Hydrogen mass fraction see Kuhlen 2011 */
double MolecularHFrac(int i)
{
  double fH2 = 1.0;
  double chi, s = 10000.;       //initialize big to ensure fH2 if necessary 
  double density, a3inv = 1.0;
//double radius; 
  double gradD;
  double metallicity;
  double sigma = 0.;
  double tau = 0.;


  if(All.ComovingIntegrationOn)
    a3inv = 1 / (All.Time * All.Time * All.Time);

  density = SphP[i].Density * a3inv;

#ifdef GFM_STELLAR_EVOLUTION
  if(SphP[i].Metallicity > 1.0e-2 * SOLAR_METALLICITY)
    metallicity = SphP[i].Metallicity;
  else
    metallicity = 1.e-2 * SOLAR_METALLICITY;
#else
  metallicity = 1.e-2 * SOLAR_METALLICITY;
#endif

/* There is also the possibility of obtaining the cell radius with the
   function get_cell_radius(i) defined in voronoi.c */

/*  Alternative for the column mass density \Sigma in Kuhlen+2011 */
/*
radius = pow(0.75 * SphP[i].Volume / (a3inv * M_PI), 1./3);
sigma = density * All.HubbleParam * radius;
tau = 0.067 * ((All.UnitMass_in_g / SOLAR_MASS) * (pow(PARSEC / All.UnitLength_in_cm, 2.0))) * metallicity / SOLAR_METALLICITY * sigma;
*/

/*  Alternative for the column density \Sigma in Gnedin, Tassi & Kravtsov 2009 */

  gradD = sqrt(SphP[i].Grad.drho[0] * SphP[i].Grad.drho[0] + SphP[i].Grad.drho[1] * SphP[i].Grad.drho[1] + SphP[i].Grad.drho[2] * SphP[i].Grad.drho[2]) * pow(a3inv, 0.66);     // LVS this shouldn't be ( )*a3inv * a3inv, but ainv^2;
  if(gradD > 0)
    {
      sigma = density * All.HubbleParam * (density / gradD);
      tau = 0.067 * ((All.UnitMass_in_g / SOLAR_MASS) * (pow(PARSEC / All.UnitLength_in_cm, 2.0))) * metallicity / SOLAR_METALLICITY * sigma;
    }

//printf("gradD=%f sigma=%f tau=%f\n",gradD,sigma,tau);

  chi = 2.3 * (1.0 + 3.1 * pow(metallicity / SOLAR_METALLICITY, 0.365)) / 3.0;
  if(tau > 0)
    {
      s = log(1.0 + chi * (0.6 + 0.01 * chi)) / (0.6 * tau);
    }
  fH2 = 1.0 - 0.75 * (s / (1.0 + 0.25 * s));
  if(fH2 < 0.0)
    fH2 = 0.0;

/*printf("radius %e, tau %e, Fh2 %e, s %e dens %e\n",radius, tau, fH2, s, density);*/
#ifdef OUTPUT_OPTICAL_DEPTH
  SphP[i].OpticalDepth = tau;
#endif
//printf("tau %e\n", SphP[i].MolecularFrac);
  return fH2;
}

#endif /* closes COMPUTE_SFR_FROM_H2 */
