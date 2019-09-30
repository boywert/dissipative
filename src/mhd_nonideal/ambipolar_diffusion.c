/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MHD_non_ideal/ambipolar_diffusion.c
 * \date        04/2016
 * \author      Mark Vogelsberger and Federico Marinacci 
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"

#ifdef AMBIPOLAR_DIFFUSION

#ifdef MHD_POWELL
void calc_ambipolar_fluxes(double B[3], double J[3], double K_dot_J[3])
{
  double K[3][3];
  double B_normed[3] = {0.0, 0.0, 0.0};
  double B_norm;

  B_norm = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

  if(B_norm > 0.0)
    {
      B_normed[0] = B[0] / B_norm;
      B_normed[1] = B[1] / B_norm;
      B_normed[2] = B[2] / B_norm;
    }
 
  K[0][0] = -1.0 + B_normed[0] * B_normed[0];
  K[1][1] = -1.0 + B_normed[1] * B_normed[1];
  K[2][2] = -1.0 + B_normed[2] * B_normed[2];
  K[0][1] = B_normed[0] * B_normed[1];
  K[1][2] = B_normed[1] * B_normed[2];
  K[1][0] = B_normed[1] * B_normed[0];
  K[1][2] = B_normed[1] * B_normed[2];
  K[2][0] = B_normed[2] * B_normed[0];
  K[2][1] = B_normed[2] * B_normed[1];

  K_dot_J[0] = K[0][0] * J[0] + K[0][1] * J[1] + K[0][2] * J[2];
  K_dot_J[1] = K[1][0] * J[0] + K[1][1] * J[1] + K[1][2] * J[2];
  K_dot_J[2] = K[2][0] * J[0] + K[2][1] * J[1] + K[2][2] * J[2];
}



void face_add_ambipolar_fluxes(const struct state *state_L, const struct state *state_R, const struct state *delta_time_L,
                           const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime,
                           double sqrtatime)
{
  double Bleft[3], Bright[3], Cleft[3], Cright[3];
  double K_dot_J[3], erg_flux[3], Bavg[3], Cavg[3], Bnormsq;
 
#ifdef ONLY_AMBIPOLAR_DIFFUSION
  face_clear_fluxes(flux);
#endif

  Bleft[0] = state_L->Bx + delta_time_L->Bx;
  Bleft[1] = state_L->By + delta_time_L->By;
  Bleft[2] = state_L->Bz + delta_time_L->Bz;

  Bright[0] = state_R->Bx + delta_time_R->Bx;
  Bright[1] = state_R->By + delta_time_R->By;
  Bright[2] = state_R->Bz + delta_time_R->Bz;

  Bavg[0] = 0.5 * (Bleft[0] + Bright[0]);
  Bavg[1] = 0.5 * (Bleft[1] + Bright[1]);
  Bavg[2] = 0.5 * (Bleft[2] + Bright[2]);

  Bnormsq = Bavg[0] * Bavg[0] + Bavg[1] * Bavg[1] + Bavg[2] * Bavg[2];

  Cleft[0] = state_L->CurlB[0];
  Cleft[1] = state_L->CurlB[1];
  Cleft[2] = state_L->CurlB[2];

  Cright[0] = state_R->CurlB[0];
  Cright[1] = state_R->CurlB[1];
  Cright[2] = state_R->CurlB[2];

  Cavg[0] = 0.5 * (Cleft[0] + Cright[0]);
  Cavg[1] = 0.5 * (Cleft[1] + Cright[1]);
  Cavg[2] = 0.5 * (Cleft[2] + Cright[2]);

  calc_ambipolar_fluxes(Bavg, Cavg, K_dot_J); 

  flux->B[0] -= All.AmbipolarDiffusionCoefficient * Bnormsq * (K_dot_J[1] * geom->nz - K_dot_J[2] * geom->ny);
  flux->B[1] -= All.AmbipolarDiffusionCoefficient * Bnormsq * (K_dot_J[2] * geom->nx - K_dot_J[0] * geom->nz);
  flux->B[2] -= All.AmbipolarDiffusionCoefficient * Bnormsq * (K_dot_J[0] * geom->ny - K_dot_J[1] * geom->nx);

  erg_flux[0] = K_dot_J[1] * Bavg[2] - K_dot_J[2] * Bavg[1];
  erg_flux[1] = K_dot_J[2] * Bavg[0] - K_dot_J[0] * Bavg[2]; 
  erg_flux[2] = K_dot_J[0] * Bavg[1] - K_dot_J[1] * Bavg[0];

  flux->energy += All.AmbipolarDiffusionCoefficient * Bnormsq * (erg_flux[0] * geom->nx + erg_flux[1] * geom->ny + erg_flux[2] * geom->nz); 
}
#endif


#endif
