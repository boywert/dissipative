/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MHD_non_ideal/ohmic_diffusion.c
 * \date        04/2015
 * \author      Federico Marinacci
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

#ifdef OHMIC_DIFFUSION

#ifdef MHD_POWELL
void face_add_ohmic_fluxes(const struct state *state_L, const struct state *state_R, const struct state *delta_time_L,
                           const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime,
                           double sqrtatime)
{
  double gradB, heat_L[3], heat_R[3], heat[3], Bleft[3], Bright[3], deltaB[3];

#ifdef ONLY_OHMIC_DIFFUSION
  face_clear_fluxes(flux);
#endif

  //deltaB[0] = state_L->dcp[0] * state_L->grad->dB[0][0] + state_L->dcp[1] * state_L->grad->dB[0][1] + state_L->dcp[2] * state_L->grad->dB[0][2];
  //deltaB[1] = state_L->dcp[0] * state_L->grad->dB[1][0] + state_L->dcp[1] * state_L->grad->dB[1][1] + state_L->dcp[2] * state_L->grad->dB[1][2];
  //deltaB[2] = state_L->dcp[0] * state_L->grad->dB[2][0] + state_L->dcp[1] * state_L->grad->dB[2][1] + state_L->dcp[2] * state_L->grad->dB[2][2];
  deltaB[0] = 0;
  deltaB[1] = 0;
  deltaB[2] = 0;

  /* extrapolate B in time */
  Bleft[0] = state_L->Bx + delta_time_L->Bx + deltaB[0];
  Bleft[1] = state_L->By + delta_time_L->By + deltaB[1];
  Bleft[2] = state_L->Bz + delta_time_L->Bz + deltaB[2];

  //deltaB[0] = state_R->dcp[0] * state_R->grad->dB[0][0] + state_R->dcp[1] * state_R->grad->dB[0][1] + state_R->dcp[2] * state_R->grad->dB[0][2];
  //deltaB[1] = state_R->dcp[0] * state_R->grad->dB[1][0] + state_R->dcp[1] * state_R->grad->dB[1][1] + state_R->dcp[2] * state_R->grad->dB[1][2];
  //deltaB[2] = state_R->dcp[0] * state_R->grad->dB[2][0] + state_R->dcp[1] * state_R->grad->dB[2][1] + state_R->dcp[2] * state_R->grad->dB[2][2];
  deltaB[0] = 0;
  deltaB[1] = 0;
  deltaB[2] = 0;

  Bright[0] = state_R->Bx + delta_time_R->Bx + deltaB[0];
  Bright[1] = state_R->By + delta_time_R->By + deltaB[1];
  Bright[2] = state_R->Bz + delta_time_R->Bz + deltaB[2];

  /* NOTE: both state_L and state_R are not rescaled by sqrtatime so it is done here */
  /* diffusion of Bx component */
  gradB = (Bright[0] - Bleft[0]) / geom->nn;
  flux->B[0] -= All.OhmicDiffusionCoefficient * gradB / (atime * sqrtatime);

  /* diffusion of By component */
  gradB = (Bright[1] - Bleft[1]) / geom->nn;
  flux->B[1] -= All.OhmicDiffusionCoefficient * gradB / (atime * sqrtatime);

  /* diffusion of Bz component */
  gradB = (Bright[2] - Bleft[2]) / geom->nn;
  flux->B[2] -= All.OhmicDiffusionCoefficient * gradB / (atime * sqrtatime);

  /* now ohmic heating */
  heat_L[0] = state_L->CurlB[1] * Bleft[2] - state_L->CurlB[2] * Bleft[1];
  heat_L[1] = state_L->CurlB[2] * Bleft[0] - state_L->CurlB[0] * Bleft[2];
  heat_L[2] = state_L->CurlB[0] * Bleft[1] - state_L->CurlB[1] * Bleft[0];

  heat_R[0] = state_R->CurlB[1] * Bright[2] - state_R->CurlB[2] * Bright[1];
  heat_R[1] = state_R->CurlB[2] * Bright[0] - state_R->CurlB[0] * Bright[2];
  heat_R[2] = state_R->CurlB[0] * Bright[1] - state_R->CurlB[1] * Bright[0];

  /* take the average (should be second order accurate) */
  heat[0] = 0.5 * (heat_L[0] + heat_R[0]);
  heat[1] = 0.5 * (heat_L[1] + heat_R[1]);
  heat[2] = 0.5 * (heat_L[2] + heat_R[2]);

  /* NOTE: both curls and B field are not rescaled by sqrtatime so the cosmological factor is a^2 instead of a */
  flux->energy -= All.OhmicDiffusionCoefficient * (heat[0] * geom->nx + heat[1] * geom->ny + heat[2] * geom->nz) / (atime * atime);
}
#endif

#ifdef MHD_CT
void face_add_ohmic_A_fluxes(int q, int qother, const struct state *st_L, const struct state *st_R, struct fluxes *flux,
                             double atime, double sqrtatime)
{
  double pUpdate[3], pOther[3];
  double Aleft[3], Aright[3];
  double dx, dy, dz, dr;

#ifdef ONLY_OHMIC_DIFFUSION
  face_clear_fluxes(flux);
#endif

  /* computing distance between the two DP points */
  pUpdate[0] = Mesh.DP[q].x;
  pUpdate[1] = Mesh.DP[q].y;
  pUpdate[2] = Mesh.DP[q].z;

  pOther[0] = Mesh.DP[qother].x;
  pOther[1] = Mesh.DP[qother].y;
  pOther[2] = Mesh.DP[qother].z;

  dx = pUpdate[0] - pOther[0];
  dy = pUpdate[1] - pOther[1];
  dz = pUpdate[2] - pOther[2];
  dr = sqrt(dx * dx + dy * dy + dz * dz);


  /* get states */
  Aleft[0] = st_L->Ax;
  Aleft[1] = st_L->Ay;
  Aleft[2] = st_L->Az;

  Aright[0] = st_R->Ax;
  Aright[1] = st_R->Ay;
  Aright[2] = st_R->Az;

  for(int s = 0; s < 3; s++)
    flux->A[s] -= All.OhmicDiffusionCoefficient * (Aright[s] - Aleft[s]) / dr;
}

void face_add_ohmic_heating(const struct state *st_L, const struct state *st_R, const struct state *delta_time_L,
                            const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime,
                            double sqrtatime)
{
  double heat_L[3], heat_R[3], heat[3], Bleft[3], Bright[3];

  /* extrapolate B in time */
  Bleft[0] = st_L->Bx + delta_time_L->Bx;
  Bleft[1] = st_L->By + delta_time_L->By;
  Bleft[2] = st_L->Bz + delta_time_L->Bz;

  Bright[0] = st_R->Bx + delta_time_R->Bx;
  Bright[1] = st_R->By + delta_time_R->By;
  Bright[2] = st_R->Bz + delta_time_R->Bz;

  /* get heating for left and right states */
  heat_L[0] = st_L->CurlB[1] * Bleft[2] - st_L->CurlB[2] * Bleft[1];
  heat_L[1] = st_L->CurlB[2] * Bleft[0] - st_L->CurlB[0] * Bleft[2];
  heat_L[2] = st_L->CurlB[0] * Bleft[1] - st_L->CurlB[1] * Bleft[0];

  heat_R[0] = st_R->CurlB[1] * Bright[2] - st_R->CurlB[2] * Bright[1];
  heat_R[1] = st_R->CurlB[2] * Bright[0] - st_R->CurlB[0] * Bright[2];
  heat_R[2] = st_R->CurlB[0] * Bright[1] - st_R->CurlB[1] * Bright[0];

  /* take the average (should be second order accurate) */
  heat[0] = 0.5 * (heat_L[0] + heat_R[0]);
  heat[1] = 0.5 * (heat_L[1] + heat_R[1]);
  heat[2] = 0.5 * (heat_L[2] + heat_R[2]);

  /* NOTE: both curls and B field are not rescaled by sqrtatime so the cosmological factor is a^2 instead of a */
  flux->energy -= All.OhmicDiffusionCoefficient * (heat[0] * geom->nx + heat[1] * geom->ny + heat[2] * geom->nz) / (atime * atime);
}
#endif

#endif
