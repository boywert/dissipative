/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/atomic_dm/cooling_atomic_dm_proto.h
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

double ADM_CoolingRateFromU(double u, double rho, double *ne_guess);
double ADM_DoCooling(double u_old, double rho, double dt, double *ne_guess);
void ADM_InitCool(void);
void ADM_cool_cell(int i);
