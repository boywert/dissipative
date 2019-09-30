/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/riemann_rosunov.c
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

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"
//#include "../mesh.h"

#ifdef MRT_RIEMANN_ROSUNOV


static void rosunov_get_fluxes_from_state_RT(struct state *st, struct fluxes *flux)
{
  for(int num1=0;num1<MRT_BINS;num1++)
    {
#ifndef MRT_COMOVING
      flux->DensPhot[num1] = st->RT_F[num1][0] ;
      
      flux->RT_F[num1][0] = c_internal_units*c_internal_units*st->PT[num1][0][0] ;
      flux->RT_F[num1][1] = c_internal_units*c_internal_units*st->PT[num1][1][0] ;
      flux->RT_F[num1][2] = c_internal_units*c_internal_units*st->PT[num1][2][0] ;
#else
      flux->DensPhot[num1] = st->RT_F[num1][0] + st->velx*st->DensPhot[num1];
      
      flux->RT_F[num1][0] = c_internal_units*c_internal_units*st->PT[num1][0][0] + st->velx*st->RT_F[num1][0] ;
      flux->RT_F[num1][1] = c_internal_units*c_internal_units*st->PT[num1][1][0] + st->velx*st->RT_F[num1][1];
      flux->RT_F[num1][2] = c_internal_units*c_internal_units*st->PT[num1][2][0] + st->velx*st->RT_F[num1][2];

#endif
    }

}

static void get_rosunov_fluxes_RT(const struct state *st_L, const struct state *st_R, const struct fluxes *flux_L, const struct fluxes *flux_R, struct fluxes *rosunov_flux, double S)
{
  for(int num1=0;num1<MRT_BINS;num1++)
    {
      rosunov_flux->DensPhot[num1] = 0.5 * (flux_L->DensPhot[num1] + flux_R->DensPhot[num1] - S * (st_R->DensPhot[num1] - st_L->DensPhot[num1]));
      
      rosunov_flux->RT_F[num1][0] = 0.5 * (flux_L->RT_F[num1][0] + flux_R->RT_F[num1][0] - S * (st_R->RT_F[num1][0] - st_L->RT_F[num1][0]));
      rosunov_flux->RT_F[num1][1] = 0.5 * (flux_L->RT_F[num1][1] + flux_R->RT_F[num1][1] - S * (st_R->RT_F[num1][1] - st_L->RT_F[num1][1]));
      rosunov_flux->RT_F[num1][2] = 0.5 * (flux_L->RT_F[num1][2] + flux_R->RT_F[num1][2] - S * (st_R->RT_F[num1][2] - st_L->RT_F[num1][2]));
    }

}

double godunov_flux_3d_rosunov_RT(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{


  int calc_flag = 1 ;
  for(int num1=0;num1<MRT_BINS;num1++)
    {
      if(st_L->DensPhot[num1]>=0.0 && st_R->DensPhot[num1]>=0.0)
	calc_flag *= 1 ;
      else
	calc_flag *= 0 ;
    }

  if(calc_flag)
    {
      double S_plus;
      struct fluxes flux_L, flux_R;

      //#if defined(MRT_COMOVING) && defined(VORONOI_STATIC_MESH)
      S_plus = dmax(dmax(fabs(st_L->velx - c_internal_units), fabs(st_R->velx - c_internal_units)), dmax(fabs(st_L->velx + c_internal_units), fabs(st_R->velx + c_internal_units)));
      //#else
	// S_plus = c_internal_units ;
      //#endif

      /* compute fluxes for the left and right states */
      rosunov_get_fluxes_from_state_RT(st_L, &flux_L);
      rosunov_get_fluxes_from_state_RT(st_R, &flux_R);

      /* compute Rosunov fluxes */
      get_rosunov_fluxes_RT(st_L, st_R, &flux_L, &flux_R, flux, S_plus);
      
      }
  else
    {
      terminate("Ngamma is negative\n");
    }

  return 0 ;
}

#endif
