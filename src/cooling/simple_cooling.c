/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cooling/simple_cooling.c
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

#include "../allvars.h"
#include "../proto.h"

/* User-defined "simple" cooling routines of arbitrary parametric form that replace
   standard cooling functions" */

#ifdef COOLING
#ifdef SIMPLE_COOLING
double DoSimpleCoolingHeating(int i, double dtcool)
{

#ifdef CIRCUMSTELLAR_IRRADIATION
  return circumstellar_DoCoolingHeating(i, dtcool);
  //return 0;
#endif
  return 0;

}
#endif
#endif
