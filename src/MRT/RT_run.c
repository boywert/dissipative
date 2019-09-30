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
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "../allvars.h"
//#include "../proto.h"
#include "../domain.h"
#include "../voronoi.h"
#include "RT.h"
#include "RT_proto.h"


#ifdef MRT

void mrt_run()
{
  for(int i;i<All.RTNumSubCycles;i++)
    {
      mpi_printf("RT: Sub Cycle %d\n", i) ;
      update_primitive_variables_RT() ;
      exchange_primitive_variables_RT() ;
      calculate_gradients_RT() ;
      exchange_primitive_variables_and_gradients_RT() ;
      compute_interface_fluxes_RT(&Mesh) ;
    }
}

#endif
