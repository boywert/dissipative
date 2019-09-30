/*! * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MRT/RT.c
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


/*! \file RT_semi_HYPRE.c
 *  \brief main driver for an moment based RT with the VET formalism
 *
 *  
 */



#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "../allvars.h"
//#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"


#ifdef MRT

void update_primitive_variables_RT()
{

  mpi_printf("RT: Updating primitive variables...\n") ;

  int idx, i;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;


      double nH_times_volume  = P[i].Mass ;
      
      SphP[i].HI = SphP[i].nHI / nH_times_volume ;
      SphP[i].HII = SphP[i].nHII / nH_times_volume ;
      SphP[i].Ne = SphP[i].ne / nH_times_volume ;
#ifdef MRT_INCLUDE_HE
      SphP[i].HeI = SphP[i].nHeI / nH_times_volume ;
      SphP[i].HeII = SphP[i].nHeII / nH_times_volume ;
      SphP[i].HeIII = SphP[i].nHeIII / nH_times_volume ;
#endif
      
  
      for(int num1=0;num1<MRT_BINS;num1++)
	{
	  SphP[i].OldCons_DensPhot[num1] = SphP[i].Cons_DensPhot[num1] ;

	  SphP[i].DensPhot[num1] = SphP[i].Cons_DensPhot[num1] / SphP[i].Volume ;
	  
	  if(SphP[i].DensPhot[num1] < MINDENSPHOT)
	    {
	      SphP[i].DensPhot[num1] = MINDENSPHOT ;
	      SphP[i].Cons_DensPhot[num1] = SphP[i].DensPhot[num1] * SphP[i].Volume ;
	    }
	  
	  for(int j=0;j<3;j++)
	    {
	      SphP[i].RT_F[num1][j] = SphP[i].Cons_RT_F[num1][j] / SphP[i].Volume ;
	      SphP[i].FN[num1][j] = SphP[i].RT_F[num1][j] / SphP[i].DensPhot[num1] ;
	    }
	}
      
      
      set_VET_single(i, SphP) ;
      
      for(int num1=0;num1<MRT_BINS;num1++)
	{
	  double modF = sqrt(SphP[i].RT_F[num1][0]*SphP[i].RT_F[num1][0] + SphP[i].RT_F[num1][1]*SphP[i].RT_F[num1][1] + SphP[i].RT_F[num1][2]*SphP[i].RT_F[num1][2]) ;
	  SphP[i].modFN[num1] = modF/SphP[i].DensPhot[num1] ;
	  if(!(SphP[i].modFN[num1] <= c_internal_units*fTOLERENCE))
	    terminate("%g\t%g\t%g\t%g\t%g\t%g\n", SphP[i].modFN[num1], SphP[i].RT_F[num1][0], SphP[i].RT_F[num1][1], SphP[i].RT_F[num1][2], SphP[i].DensPhot[num1], c_internal_units) ; 
	  
	  if(SphP[i].DensPhot[num1]< 0.0)
	    terminate("UPDATE PRIMITIVE VARIABLES Ngamma < 0.0\n") ;
	  
	  if(isnan(SphP[i].RT_F[0][0]))
	    terminate("Nan F\n") ;
	}
    }
}

#endif
