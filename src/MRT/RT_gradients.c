/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/gradients.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "../voronoi.h"


int N_Grad_RT = 0;

struct grad_elements grad_elements_RT[MAXRTGRADIENTS] ;

void init_gradients_RT()
{
  int j, num1 ;

  for(num1=0;num1<MRT_BINS;num1++)
    {
      for(j=0;j<DIMS;j++)
	{
	  gradient_init_RT(&SphP[0].FN[num1][j], &RTPrimExch[0].FN[num1][j], SphP[0].RTGrad.dFN[num1][j], GRADIENT_TYPE_NORMAL) ;
#ifdef MRT_TIME_EXTRAPOLATION
	  for(int k=0;k<DIMS;k++)
	    gradient_init_RT(&SphP[0].PT[num1][j][k], &RTPrimExch[0].PT[num1][j][k], SphP[0].RTGrad.dPT[num1][j][k], GRADIENT_TYPE_NORMAL) ;
#endif
	}
      
      gradient_init_RT(&SphP[0].modFN[num1], &RTPrimExch[0].modFN[num1], SphP[0].RTGrad.dmodFN[num1], GRADIENT_TYPE_NORMAL) ;
      
      gradient_init_RT(&SphP[0].DensPhot[num1], &RTPrimExch[0].DensPhot[num1], SphP[0].RTGrad.dDensPhot[num1], GRADIENT_TYPE_NORMAL) ;
    }
  mpi_printf("INIT RT: %d/%d Gradients used.\n", N_Grad_RT, MAXRTGRADIENTS);
}

void gradient_init_RT(MyFloat * addr, MyFloat * addr_exch, MySingle * addr_grad, int type)
{
  if(N_Grad_RT == MAXRTGRADIENTS)
    {
      mpi_printf("Failed to register gradient, maximum of %d already reached\n", MAXGRADIENTS);
      terminate("MAXRTGRADIENTS reached");
    }

  grad_elements_RT[N_Grad_RT].type = type;
  
  grad_elements_RT[N_Grad_RT].offset = ((char *) addr) - ((char *) &SphP[0]);

  grad_elements_RT[N_Grad_RT].offset_exch = ((char *) addr_exch) - ((char *) &RTPrimExch[0]);
  grad_elements_RT[N_Grad_RT].offset_grad = ((char *) addr_grad) - ((char *) &(SphP[0].RTGrad));


  N_Grad_RT++;
}
