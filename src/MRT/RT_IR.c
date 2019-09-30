/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_cooling.c
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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"

#ifdef MRT_IR_ONLY_CHEMISTRY

void mrt_IR_chemistry(void) 
{
  int idx, i ;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      double dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ;
      
      for(int num1=UV_BINS;num1<(UV_BINS+IR_BINS);num1++)
	{
	  double A = dt * SphP[i].KappaIR_R[num1-UV_BINS] * c_internal_units ;
	  double ratio = 1.0 / (1 + A) ;
	  //double ratio = exp(-A) ;

	  
	  //SphP[i].Cons_DensPhot[num1] *=ratio ;
	  for(int j=0;j<3;j++)
	    {
	      SphP[i].Cons_RT_F[num1][j] *= ratio ;
	      if(isnan(SphP[i].Cons_RT_F[num1][j]))
		terminate("IR Cooling - What!!\n") ;
	    }
	}
    }
  return ;
}
#endif

#ifdef MRT_IR_LTE
#ifdef MRT_IR_LTE_SEMI_IMPLICIT
double mrt_update_IR_cooling(int i, double dt) 
{
  if(dt<=0.0 || (P[i].ID>=40000000 && P[i].ID <= 49999999))
    return 0.0 ;
  else
    {
      double cspeed =  2.99792458e10 / All.UnitVelocity_in_cm_per_s ;
      double molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);
      
      double unew = SphP[i].Utherm*SphP[i].Density ;
      double uold = unew ;

      molecular_weight = 1.0 ;
      
      double inv_cv =  GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam) / SphP[i].Density ;

      double temp = uold * inv_cv ;
      
      for(int num1=UV_BINS;num1<(UV_BINS+IR_BINS);num1++)
	{
	  /*First take care of flux absorption*/
	  double A = dt * SphP[i].KappaIR_R[num1-UV_BINS] * c_internal_units ;
	  double ratio = 1.0 / (1.0 + A) ;
	  
	  for(int j=0;j<3;j++)
	    {
	      SphP[i].Cons_RT_F[num1][j] *= ratio ;
	      if(isnan(SphP[i].Cons_RT_F[num1][j]))
		terminate("IR Cooling - What!!\n") ;
	    }

	  double Enew = SphP[i].Cons_DensPhot[num1]/SphP[i].Volume ;
	  double Eold = Enew ;


	  int times = 1 ;
	  int flag = 1 ;
	  do 
	    {
	      for(int kk=0; kk<times; kk++)
		{
		  double de = (cspeed*radiation_constant*pow(temp, 4) - c_internal_units*Enew)/( 1.0/(dt * SphP[i].KappaIR_P[num1-UV_BINS]) + c_internal_units + 4.0*cspeed*radiation_constant*pow(temp,3)*inv_cv) ;
		  if(fabs(de/Enew) >0.1 || fabs(de/uold) > 0.1)
		    {
		      times = times * 2 ;
		      dt = dt/2.0 ;
		      flag = 1 ;
		      unew = SphP[i].Utherm*SphP[i].Density ;
		      Enew = SphP[i].Cons_DensPhot[num1]/SphP[i].Volume ;
		      temp = unew * inv_cv ;
		      break ;
		    }
		  else
		    {
		      unew = unew - de ;
		      Enew = Enew + de ;
		      temp = unew * inv_cv ;
		      flag = 0 ;
		    }
		}
	    }
	  while(flag) ;

	  SphP[i].Cons_DensPhot[num1] = Enew * SphP[i].Volume ;
	}
      return (unew - uold)/SphP[i].Density ;
    }
}
#endif




#ifdef MRT_IR_LTE_GSL
double mrt_update_IR_cooling(int i, double dt)
{
  double du = 0.0;
  if(dt<=0.0 
#ifdef BOUNDARY_REFL_SOLIDSIDE_MINID
     || (P[i].ID >= BOUNDARY_REFL_SOLIDSIDE_MINID && P[i].ID <= BOUNDARY_REFL_SOLIDSIDE_MAXID)
#endif
     )
    return du ;
  else
    {
      int num1 = UV_BINS ;
      /*First take care of flux absorption*/
      double A = dt * SphP[i].KappaIR_R[num1-UV_BINS] * c_internal_units ;
      //double ratio = 1.0 / (1.0 + A) ;
      double ratio = exp(-A) ;
      for(int j=0;j<3;j++)
	{
	  SphP[i].Cons_RT_F[num1][j] *= ratio ;
	  if(isnan(SphP[i].Cons_RT_F[num1][j]))
	    terminate("IR Cooling - What!!\n") ;
	}


      /*Solve the stiff equation*/

      double cspeed =  2.99792458e10 / All.UnitVelocity_in_cm_per_s ;
      double molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].Ne);

      //double molecular_weight = 2.33 ;

      double unew = SphP[i].Utherm*SphP[i].Density ;
      double uold = unew ;
    
#ifdef MRT_IR_PHOTON_TRAPPING
      double Enew = (SphP[i].Cons_DensPhot[num1] + SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS])/SphP[i].Volume ;
      double trap_frac = SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] / ( SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] + SphP[i].Cons_DensPhot[num1]) ;
#else
      double Enew = SphP[i].Cons_DensPhot[num1]/SphP[i].Volume ;
#endif
      double Eold = Enew ;


      if(uold <0.0 || Eold < 0.0)
	terminate("Something Wrong!!\n") ;


      double inv_cv =  GAMMA_MINUS1 * (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) / (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam) / SphP[i].Density ;

      double temp = uold * inv_cv ;

      double de = (cspeed*radiation_constant*pow(temp, 4) - c_internal_units*Enew)/( 1.0/(dt * SphP[i].KappaIR_P[num1-UV_BINS]) + c_internal_units + 4.0*cspeed*radiation_constant*pow(temp,3)*inv_cv) ; /*Semi-implicit estimate*/

      if(fabs(de/unew) < 0.1 && fabs(de/Enew)<0.1)
	{
	  /*Keep the semi-implicit solution*/
	  unew = unew - de ;
	  Enew = Enew + de ;

	  if(unew<0.0 || Enew <0.0)
	    terminate("Something wrong!!\n") ;


#ifdef MRT_IR_PHOTON_TRAPPING
	  SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] = trap_frac * Enew * SphP[i].Volume ;
	  SphP[i].Cons_DensPhot[num1] = (1.0 - trap_frac) * Enew * SphP[i].Volume ;
#else	  
	  SphP[i].Cons_DensPhot[num1]  = Enew * SphP[i].Volume ;
#endif
	}
      else
	{
	  /*Solve using the implicit Bulirsch-Stoer method of Bader and Deuflhard*/
	  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_bsimp;
	  gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(T, 2);
	  gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(1e-6, 0.0);
	  gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(2);
	  double coeff[2] = { SphP[i].KappaIR_P[num1-UV_BINS]*c_internal_units, SphP[i].KappaIR_P[num1-UV_BINS]*cspeed*radiation_constant*pow(inv_cv, 4.0)};
	  gsl_odeiv2_system sys = { mrt_IR_rate_ODEs, jacobian, 2, coeff };
	  double t = 0.0, t1 = dt ;
	  double h = 0.00025 * dt ;
	  double y[2] = {Enew, unew} ;

	  while(t<t1)
	    {
	      int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t1, &h, y);
	      if(status != GSL_SUCCESS)
		terminate("we seem to have failed in the integration ----  Tgas = %g \t Trad = %g \t\t Position of particle = %g %g %g\t\t Coeffss = %g %g \t\t ratio = %g\n", uold*inv_cv, pow(Eold/radiation_constant, (1./4.)), P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], coeff[0], coeff[1], coeff[0]*y[0]/(coeff[1]*y[1]*y[1]*y[1]*y[1]));
	    }



	  gsl_odeiv2_evolve_free(e);
	  gsl_odeiv2_control_free(c);
	  gsl_odeiv2_step_free(s);
	  
	  Enew = y[0] ;
	  unew = y[1] ;


	  /*	  if(uold>Eold && unew>uold)
	    {
	      unew = uold ;
	      Enew = Eold ;
	    }

	  if(uold<Eold && unew<uold)
	    {
	      unew = uold ;
	      Enew = Eold ;
	    }
	  */
	  if(isnan(Enew) || isnan(unew) || unew/SphP[i].Density<0.0 || Enew<0.0)
	    terminate("Enew/unew is a nan/negative uold=%g \t Eold=%g unew = %g \t Enew = %g \n", uold, Eold, unew, Enew) ;
	  
	  
#ifdef MRT_IR_PHOTON_TRAPPING
	  SphP[i].Trapped_Cons_DensPhot[num1-UV_BINS] = trap_frac * Enew * SphP[i].Volume ;
	  SphP[i].Cons_DensPhot[num1] = (1.0 - trap_frac) * Enew * SphP[i].Volume ;
#else	  
	  SphP[i].Cons_DensPhot[num1]  = Enew * SphP[i].Volume ;
#endif
	  
	}
      return (unew-uold)/SphP[i].Density ;
    }
}

int mrt_IR_rate_ODEs(double t, const double y[], double f[], void *params)
{
  double *coeff = (double *) params;
  f[0] = -coeff[0]*y[0] + coeff[1]*y[1]*y[1]*y[1]*y[1] ;
  f[1] = coeff[0]*y[0] - coeff[1]*y[1]*y[1]*y[1]*y[1] ;
  return GSL_SUCCESS;
}

int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  double *coeff = (double *) params;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);

  gsl_matrix *m_ptr = &dfdy_mat.matrix;/* m_ptr points to the matrix */

  /* fill the Jacobian matrix as shown */
  gsl_matrix_set (m_ptr, 0, 0, -coeff[0]);/* df[0]/dy[0] = 0 */
  gsl_matrix_set (m_ptr, 0, 1, coeff[1]*4*y[1]*y[1]*y[1]);/* df[0]/dy[1] = 1 */
  gsl_matrix_set (m_ptr, 1, 0, coeff[0]); /* df[1]/dy[0] */
  gsl_matrix_set (m_ptr, 1, 1, -coeff[1]*4*y[1]*y[1]*y[1]);     /* df[1]/dy[1] */

  /* set explicit t dependence of f[i] */
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;

  return GSL_SUCCESS ;
}
#endif
#endif


