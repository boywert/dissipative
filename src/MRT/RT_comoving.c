/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT_comoving.c
 * \date        MM/YYYY
 * \author      Rahul Kannan
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */


/*! \file RT_comoving.c
 *  \brief functions to solve the RT equations in the comoving frame 
 *
 *  
 */





#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "../allvars.h"
#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"


#ifdef MRT_COMOVING

void cell_do_lorentz_boost(int i, struct sph_particle_data *sphp, double vx, double vy, double vz)
{
  double cspeed ;

#ifdef MRT_CONSTANT_KAPPA
  cspeed = 10.0 ;
#else
  cspeed = 2.99792458e10 / All.UnitVelocity_in_cm_per_s ;
#endif

  for(int num1=0;num1<MRT_BINS;num1++)
    {

      double pt[3][3] ;
      double n[3] ;
      double f ;

      double modF = sqrt(sphp[i].Cons_RT_F[num1][0]*sphp[i].Cons_RT_F[num1][0] + sphp[i].Cons_RT_F[num1][1]*sphp[i].Cons_RT_F[num1][1] + sphp[i].Cons_RT_F[num1][2]*sphp[i].Cons_RT_F[num1][2]) ;

      if(sphp[i].Cons_DensPhot[num1]==0.0)
	f = modF/c_internal_units ;
      else
	f = modF/c_internal_units/sphp[i].Cons_DensPhot[num1] ;

      //      printf("Cons_DensPhot = %g Fx = %g Fy = %g Fz = %g\n", sphp[i].Cons_DensPhot[num1], sphp[i].Cons_RT_F[num1][0], sphp[i].Cons_RT_F[num1][1], sphp[i].Cons_RT_F[num1][2]) ;

      if(f>1.0)
        {
          double ratio = 0.99/f ;
          f = 0.99 ;
          sphp[i].Cons_RT_F[num1][0] *= ratio ;
          sphp[i].Cons_RT_F[num1][1] *= ratio ;
          sphp[i].Cons_RT_F[num1][2] *= ratio ;
	}

      double chi = (3.0 + 4.0*f*f)/(5.0 + 2.0*sqrt(4.0 - 3.0*f*f)) ;

      if(modF == 0.0)
        modF = 1.0 ;

      for(int l=0;l<3;l++)
        n[l] = sphp[i].Cons_RT_F[num1][l]/modF ;


      for(int j=0;j<3;j++)
	{
          for(int k=0;k<3;k++)
            {
              if(j==k)
                pt[j][k] =  1.0 + ((1.0 - chi)/2.0) + ((3.0*chi - 1.0)/2.0)*n[j]*n[k] ; /*Addition 1.0 + comes due to the fact that you have I+D in transformation*/
              else
                pt[j][k] = ((3.0*chi - 1.0)/2.0)*n[j]*n[k] ;

              pt[j][k] *= sphp[i].Cons_DensPhot[num1] ;
            }
	}

      double Ec = sphp[i].Cons_DensPhot[num1] - 2.0 * (vx*sphp[i].Cons_RT_F[num1][0] + vy*sphp[i].Cons_RT_F[num1][1] + vz*sphp[i].Cons_RT_F[num1][2])  / c_internal_units / cspeed ;
      double Fcx = sphp[i].Cons_RT_F[num1][0] - (c_internal_units/cspeed) * (vx*pt[0][0] + vy*pt[1][0] + vz*pt[2][0]) ;
      double Fcy = sphp[i].Cons_RT_F[num1][1] - (c_internal_units/cspeed) * (vx*pt[0][1] + vy*pt[1][1] + vz*pt[2][1]) ;
      double Fcz = sphp[i].Cons_RT_F[num1][2] - (c_internal_units/cspeed) * (vx*pt[0][2] + vy*pt[1][2] + vz*pt[2][2]) ;
      sphp[i].Cons_DensPhot[num1] = Ec ;
      sphp[i].Cons_RT_F[num1][0] = Fcx ;
      sphp[i].Cons_RT_F[num1][1] = Fcy ;
      sphp[i].Cons_RT_F[num1][2] = Fcz ;

      if(isnan(sphp[i].Cons_RT_F[num1][0]) || isnan(sphp[i].Cons_DensPhot[num1]))
	terminate("NAN ||| vx = %g \t vy = %g \t vz = %g \t pt = %g %g %g %g %g %g %g %g %g \n modF = %g, f = %g, chi = %g, cons_densphot = %g \n c_int = %g, c_tot = %g\n", vx, vy, vz, pt[0][0], pt[0][1], pt[0][2], pt[1][0], pt[1][1], pt[1][2], pt[2][0], pt[2][1], pt[2][2], modF, f, chi, sphp[i].Cons_DensPhot[num1], c_internal_units, cspeed) ;
    }
}


void do_comoving_frame_source_terms()
{
  mpi_printf("RT: Adding the comoving frame source terms...\n") ;
  int idx, i ;
 
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
	continue;

      double dt = 0.5 * (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval ; /*factor of 0.5 included because of spilt for RK time integration*/

     
      for(int num1=0;num1<MRT_BINS;num1++)
	{
	  double sumE, sumFx, sumFy, sumFz ;
	  sumE = sumFx = sumFy = sumFz = 0.0 ;
	  for(int k=0;k<3;k++)
	    {
	      sumFx += SphP[i].RT_F[num1][0]*SphP[i].Grad.dvel[k][0] ;
	      sumFy += SphP[i].RT_F[num1][1]*SphP[i].Grad.dvel[k][1] ;
	      sumFz += SphP[i].RT_F[num1][2]*SphP[i].Grad.dvel[k][2] ;
	      for(int l=0;l<3;l++)
		  sumE += SphP[i].PT[num1][k][l]*SphP[i].Grad.dvel[l][k] ;

	    }
	  SphP[i].Cons_DensPhot[num1] += SphP[i].Volume * (-dt) * sumE ;
	  SphP[i].Cons_RT_F[num1][0] += SphP[i].Volume * (-dt) * sumFx ;
	  SphP[i].Cons_RT_F[num1][1] += SphP[i].Volume * (-dt) * sumFy ;
	  SphP[i].Cons_RT_F[num1][2] += SphP[i].Volume * (-dt) * sumFz ;
	}
    }

}

#endif
