/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/runge_kutta_full.h
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

#ifndef RUNGE_KUTTA_FULL_UPDATE_H
#define RUNGE_KUTTA_FULL_UPDATE_H

#include "arepoconfig.h"

#ifdef RUNGE_KUTTA_FULL_UPDATE

#include "allvars.h"

struct conservative_variables {
  MyFloat Mass;
  MyFloat Momentum[3];
  MyFloat Energy;
#ifdef MAXSCALARS
  MyFloat Scalars[MAXSCALARS];
#endif
#ifdef MHD
  MyFloat BConserved[3];
#endif
};

void rk_save_conservative_variables();
void rk_finish_step();

void rk_derefinement_add( struct conservative_variables *target, struct conservative_variables *source, double fac );
void rk_derefinement_set( struct conservative_variables *target, struct conservative_variables *source, double fac );
void rk_multiply( struct conservative_variables *target, double fac );

#endif

#endif