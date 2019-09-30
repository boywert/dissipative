/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/network/integrate.h
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

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "../allvars.h"
#include "network.h"

void network_main(int *timeBin);
void network_normalize(double *x, double *e, const struct network_data *nd, struct network_workspace *nw);
int network_integrate(double temp, double rho, double *x, double dt, double *dedt, const struct network_data *nd, struct network_workspace *nw);
void network_composition_statistics();

#endif
