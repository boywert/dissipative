/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays/cosmic_rays.h
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

#ifndef COSMICRAYS_H
#define COSMICRAYS_H

#include "../allvars.h"
#include "../shock_finder/shock_finder_rays.h"

void init_cosmic_rays(void);

#ifdef COSMIC_RAYS_DIFFUSION
struct diff_face_data {
  int cornerFirst;
  int cornerCount;
#ifdef COSMIC_RAYS_DIFFUSION_ANISOTROPIC
  double bfld[3];
#endif
  int sng;

  int active;
  double dt;
  double area;
  double kappa_iso;
#ifdef COSMIC_RAYS_DIFFUSION_ANISOTROPIC
  double kappa_aniso;
#endif
  
  double nx, ny, nz;
  double mx, my, mz;
  double px, py, pz;

  double failWeight;
};
#endif

/* Diffusion */
#ifdef COSMIC_RAYS_DIFFUSION
void do_cr_diffusion(void);
void set_diffusion_coefficients( struct diff_face_data *diff_face_data );
#endif

/* Streaming */
#ifdef COSMIC_RAYS_STREAMING
void cosmic_rays_do_streaming(void);
double get_timestep_streaming(int i);
double compute_chi( int i );
#endif

/* Various source terms */
#ifdef COSMIC_RAYS_SN_INJECTION
void deposit_cosmic_rays_from_supernovae(void);
void cosmic_rays_find_ngbs( int *ngbCount, MyFloat *numNgb, MyFloat *normSph );
#endif

#ifdef COSMIC_RAYS_COOLING
void do_cosmic_ray_cooling(void);
#endif
#ifdef COSMIC_RAYS_ALFVEN_COOLING
void do_cosmic_ray_Alfven_cooling(void);
#endif

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
void accelerate_cosmic_rays_at_shock(shock_ray* first_ray, shock_ray* last_ray);
#endif

#endif /* COSMICRAYS_H */
