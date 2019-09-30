/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cosmic_rays_diffusion.c
 * \date        06/2015
 * \author      R. Pakmor
 * \brief       (anisotropic) Cosmic Ray diffusion
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - 02/09/2015 Added limiting for explicit and implicit time integration
 */

#include "../allvars.h"
#include "../proto.h"

#if defined(COSMIC_RAYS_DIFFUSION) && !defined(COSMIC_RAYS_DIFFUSION_OLD)

void set_diffusion_coefficients( struct diff_face_data *diff_face_data )
{
  double coeff = 1e28;
  if(All.CR_Diffusion_Coefficient > 0)
    coeff = All.CR_Diffusion_Coefficient;
  coeff = coeff / (All.UnitLength_in_cm * All.UnitLength_in_cm) * All.UnitTime_in_s;

  if(All.ComovingIntegrationOn)
    coeff *= All.cf_atime * All.cf_atime / All.HubbleParam;
  
  int iface;
  for(iface = 0; iface < Mesh.Nvf; iface++)
    {
      struct diff_face_data *fd = &diff_face_data[iface];

      if(fd->active)
        {
#ifndef COSMIC_RAYS_DIFFUSION_ANISOTROPIC
          fd->kappa_iso = coeff ;
#else
          fd->kappa_iso = 0.;
          fd->kappa_aniso = coeff;
#endif
        }
    }
}

#endif
