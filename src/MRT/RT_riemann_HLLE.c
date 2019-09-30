/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_riemann_HLLE.c
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
#include "RT_proto.h"

#ifdef MRT_RIEMANN_HLLE

void HLLE_interpolate(double, double, double *, double *) ;
void HLLE_get_minmax_eigenvalues(struct state *state, double *S_plus, double *S_minus) ;

void readin_hlle_eingenvalues() 
{
  FILE *fp ;
  fp = fopen(All.HLLEFile, "r") ;
  int f, theta ;
  double l1, l2, l3, l4 ;
  for(int num1=0; num1<101; num1++)
    {
      for(int num2=0; num2<101; num2++)
	{
	  fscanf(fp, "%d %d %lf %lf %lf %lf", &f, &theta, &l1, &l2, &l3, &l4) ;
	  lambda1[num1][num2] = l1 ;
	  lambda2[num1][num2] = l2 ;
	  lambda3[num1][num2] = l3 ;
	  lambda4[num1][num2] = l4 ;

	  //	  printf("%d\t%d\t%lf\t%lf\t%lf\t%lf\n", num1, num2, lambda1[num1][num2], lambda2[num1][num2], lambda3[num1][num2], lambda4[num1][num2]) ;
	}
    }
  fclose(fp) ;
  mpi_printf("RT: Read in the HLL eingenvalues\n") ;

  /*   mpi_printf("RT:testing HLLE interpolation scheme\n") ;

  double fd, td ;
  fd = 0.0 ;
  td = M_PI ;

  printf("\n\n\n\n") ;
  
  double min, max ;
  for(int k=0;k<1000;k++)
    {
      td =  ((double) (k))/1000.0  ;
      double m = 1.0*100.0 ;
      double n = td*100.0 ;
      HLLE_interpolate(m,n, &min, &max) ;      
      printf("%g %g %g \t %g \n", fd, td, min, max) ;
    }

  printf("\n\n\n\n") ;
  
  
  terminate("BBBB\n") ;*/
}



static void HLLE_get_fluxes_from_state_RT(struct state *st, struct fluxes *flux)
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

static void get_HLLE_fluxes_RT(const struct state *st_L, const struct state *st_R, const struct fluxes *flux_L, const struct fluxes *flux_R, struct fluxes *HLLE_flux, double *S_plus_L, double *S_minus_L, double *S_plus_R, double *S_minus_R)
{
  for(int num1=0;num1<MRT_BINS;num1++)
    {

      double lpl = S_plus_L[num1] ; 
      double lml = S_minus_L[num1] ;

      double lpr = S_plus_R[num1] ;
      double lmr = S_minus_R[num1] ;

      double lp = 0.0 ;

      if(lpl>lp)
	lp=lpl;
      if(lpr>lp)
	lp=lpr ;

      double lm = 0.0 ;

      if(lml<lm)
	lm=lml ;
      if(lmr<lm)
	lm=lmr;

#if defined(BOUNDARY_REFL_FLUIDSIDE_MINID) && defined(BOUNDARY_REFL_FLUIDSIDE_MAXID) && defined(BOUNDARY_REFL_SOLIDSIDE_MINID) && defined(BOUNDARY_REFL_SOLIDSIDE_MAXID)
      MyIDType idL = st_L->ID;
      MyIDType idR = st_R->ID;

      int nn1 ;
      if(idR >= BOUNDARY_REFL_FLUIDSIDE_MINID && idR < BOUNDARY_REFL_FLUIDSIDE_MAXID && idL >= BOUNDARY_REFL_SOLIDSIDE_MINID && idL < BOUNDARY_REFL_SOLIDSIDE_MAXID)
        {
	  lp = 1.0 ;
	  lm = -1.0 ;
	}
      else if(idL >= BOUNDARY_REFL_FLUIDSIDE_MINID && idL < BOUNDARY_REFL_FLUIDSIDE_MAXID && idR >= BOUNDARY_REFL_SOLIDSIDE_MINID && idR < BOUNDARY_REFL_SOLIDSIDE_MAXID)
	{
	  lp = 1.0 ;
	  lm = -1.0 ;
	}
#endif


      //      if((st_L->ID == 561 && st_R->ID == 611) || (st_L->ID == 611 && st_R->ID == 561))
      //printf("lid, rid, ld, rd, lp, lm = %d %d %g %g %g %g %g %g %g %g\n", st_L->ID, st_R->ID, st_L->DensPhot[num1], st_R->DensPhot[num1], lp, lm, lpl, lml, lpr, lmr) ;


      if(lp-lm == 0.0)
	{
	  HLLE_flux->DensPhot[num1] = 0.0;
	  HLLE_flux->RT_F[num1][0] = 0.0 ;
	  HLLE_flux->RT_F[num1][1] = 0.0 ;
	  HLLE_flux->RT_F[num1][2] = 0.0 ;
	}
      else
	{
	  HLLE_flux->DensPhot[num1] = (lp*flux_L->DensPhot[num1] - lm*flux_R->DensPhot[num1] + c_internal_units*lm*lp*(st_R->DensPhot[num1] - st_L->DensPhot[num1]))/(lp-lm);
      
	  HLLE_flux->RT_F[num1][0] = (lp*flux_L->RT_F[num1][0] - lm*flux_R->RT_F[num1][0] + c_internal_units*lp*lm * (st_R->RT_F[num1][0] - st_L->RT_F[num1][0]))/(lp-lm);
	  HLLE_flux->RT_F[num1][1] = (lp*flux_L->RT_F[num1][1] - lm*flux_R->RT_F[num1][1] + c_internal_units*lp*lm * (st_R->RT_F[num1][1] - st_L->RT_F[num1][1]))/(lp-lm);
	  HLLE_flux->RT_F[num1][2] = (lp*flux_L->RT_F[num1][2] - lm*flux_R->RT_F[num1][2] + c_internal_units*lp*lm * (st_R->RT_F[num1][2] - st_L->RT_F[num1][2]))/(lp-lm);
	}
    }
}



void HLLE_get_minmax_eigenvalues(struct state *st, double S_plus[MRT_BINS], double S_minus[MRT_BINS])
{
  //  printf("Here func 1\n") ;

  for(int num1=0;num1<MRT_BINS;num1++)
    {
      double modF = sqrt(st->RT_F[num1][0]*st->RT_F[num1][0] + st->RT_F[num1][1]*st->RT_F[num1][1] + st->RT_F[num1][2]*st->RT_F[num1][2]) ;
      //      if(modF == 0.0)
      //modF = 1.0 ;
      
      double red_flux = modF / c_internal_units / st->DensPhot[num1] ;

      
      double theta ;
      if(modF > MINDENSPHOT)
	theta = acos(st->RT_F[num1][0] / modF) ;
      else
	theta = 0.0 ;

      //  printf("%g\t%g\t%g\n", modF, c_internal_units, st->DensPhot[num1]) ;

      double j = red_flux*100.0 ;
      double k = theta*100.0/M_PI ;

      if(isnan(k))
	terminate("%g\t%g\t%g\t%g\t%g\t%g\n", j, k, red_flux, theta, st->RT_F[num1][0], modF) ;

      double min, max ;

      // printf("Here func 2\n") ;
      HLLE_interpolate(j,k, &min, &max) ;
      //      if(st->DensPhot[num1] < 10.0*MINDENSPHOT)
      //min = max = 0.0 ;
      
      S_plus[num1] = max ;
      S_minus[num1] = min ;

      //      if(st->ID == 561 || st->ID == 611)
      //printf("ID, red_flux, theta j, k, max, min = %d %g %g %g %g %g %g\n", st->ID, red_flux, theta, j, k, max, min) ;

	    
      //printf("Here func 3\n") ;
    }
  return ;
}


void HLLE_interpolate(double j, double k, double *min, double *max)
{

  //printf("%g\t%g\n", j, k ) ;

  int num1 = ((int) (j)) ;
  int num2 = ((int) (k)) ;

  //  printf("num1 = %d, num2 = %d \n", num1, num2) ;
  
  double rem1 = j - ((double) (num1)) ;
  double rem2 = k - ((double) (num2)) ;

  if(num1<0 || num1>100 || num2<0 || num2>100)
    terminate("What!! %g \t %g \t \n\n", j, k) ;


  double l1 = interpolate_2d(lambda1, num1, num2, rem1, rem2) ; //(lambda1[num1][num2] * (1.0-rem1) + lambda1[num1+1][num2] * rem1  + lambda1[num1][num2] * (1.0-rem2) + lambda1[num1][num2+1] * rem2)/2.0 ;
  double l2 = interpolate_2d(lambda2, num1, num2, rem1, rem2) ;
  double l3 = interpolate_2d(lambda3, num1, num2, rem1, rem2) ;
  double l4 = interpolate_2d(lambda4, num1, num2, rem1, rem2) ;
  /*
  double l2 = (lambda2[num1][num2] * (1.0-rem1) + lambda2[num1+1][num2] * rem1  + lambda2[num1][num2] * (1.0-rem2) + lambda2[num1][num2+1] * rem2)/2.0 ;
  double l3 = (lambda3[num1][num2] * (1.0-rem1) + lambda3[num1+1][num2] * rem1  + lambda3[num1][num2] * (1.0-rem2) + lambda3[num1][num2+1] * rem2)/2.0 ;
  double l4 = (lambda4[num1][num2] * (1.0-rem1) + lambda4[num1+1][num2] * rem1  + lambda4[num1][num2] * (1.0-rem2) + lambda4[num1][num2+1] * rem2)/2.0 ;
  */

  *min = 1000.0 ;
  *max = -1000.0 ;

  if(l1<*min)
    *min=l1 ;
  if(l2<*min)
    *min=l2 ;
  if(l3<*min)
    *min=l3 ;
  if(l4<*min)
    *min=l4 ;

  if(l1>*max)
    *max=l1 ;
  if(l2>*max)
    *max=l2 ;
  if(l3>*max)
    *max=l3;
  if(l4>*max)
    *max=l4;
}

double godunov_flux_3d_HLLE_RT(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
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
      //      double S_plus;
      struct fluxes flux_L, flux_R;

      //      printf("Here\n") ;
      double S_plus_R[MRT_BINS], S_minus_R[MRT_BINS] ;
      double S_plus_L[MRT_BINS], S_minus_L[MRT_BINS] ;

      HLLE_get_minmax_eigenvalues(st_L, S_plus_L, S_minus_L) ;
      HLLE_get_minmax_eigenvalues(st_R, S_plus_R, S_minus_R) ;


      /* compute fluxes for the left and right states */
      HLLE_get_fluxes_from_state_RT(st_L, &flux_L);
      HLLE_get_fluxes_from_state_RT(st_R, &flux_R);

      /* compute Rosunov fluxes */
      get_HLLE_fluxes_RT(st_L, st_R, &flux_L, &flux_R, flux, S_plus_L, S_minus_L, S_plus_R, S_minus_R);
      
      }
  else
    {
      terminate("Ngamma is negative\n");
    }

  return 0 ;
}

#endif
