/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/winds_vars.h
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

#ifndef WINDS_VARS_H
#define WINDS_VARS_H

#define GFM_FAC_SIGMA 0.57735   /* 1/sqrt(3) */

typedef struct
{
  double wind_mass;
  double wind_velocity;
  double wind_utherm;
  double v_esc_halo;
  double v_vir_halo;
#if defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE==0)
  double c_halo;
#endif
} wind_parameter;

/* data for nearest ngb search */
extern struct wind_particle
{
  MyIDType CellID;              /* ID of cell to recouple with */
  int index;                    /* index (local) of wind particle */
  MyFloat NearestDist;          /* nearest distance to ngb */
  MyFloat Density;              /* density of nearest cell */
} *WindParticle;

#endif
