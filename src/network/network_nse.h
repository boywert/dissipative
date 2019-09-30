/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/network/network_nse.h
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

#ifndef NETWORK_NSE_H
#define NETWORK_NSE_H

#include "arepoconfig.h"

#if defined (NUCLEAR_NETWORK) && defined (NETWORK_NSE)

#define NSE_MAXITER 500
#define NSE_DIFFVAR 1e-12       /* variation for numerical derivatives */

#define ELECTRON_CHARGE_ESU (4.80320427e-10)

#include "network.h"
#include "network_solver.h"

struct nse_params_ene
{
  struct network_data *nd;
  struct network_workspace *nw;
  struct network_solver_trajectory *traj;
  double temp, rho;
};

struct nse_params
{
  struct network_data *nd;
  struct network_workspace *nw;
  double temp;
  double rho;
};

struct nse_params_iter
{
  struct network_data *nd;
  struct network_workspace *nw;
  double rho;
  double ye;
  double *x;
  double energy;
};

struct network_nse_eparams
{
  struct network_data *nd;
  struct network_workspace *nw;
  double rho;
  double energy;
  double ye;
};

int network_nse_init(char *speciesfile, char *ratesfile, char *partfile, char *massesfile, char *weakratesfile, struct network_data *nd, struct network_workspace *nw);
void network_nse_deinit(struct network_data *nd, struct network_workspace *nw);
int network_nse(double temp, double rho, double ye, double *x, struct network_data *nd, struct network_workspace *nw);
double network_nse_iter(double energy, double rho, double x[], struct network_data *nd, struct network_workspace *nw);
int network_nse_egiven(double energy, double rho, double ye, double *xnuc, double *tempguess, struct network_data *nd, struct network_workspace *nw);
void network_nse_integrate_ye(double energy, double rho, double x[], double dt, double *dedt, double *temp, struct network_data *nd, struct network_workspace *nw);
void network_nse_integrate(double temp, double rho, double x[], double dt, double *dedt, struct network_data *nd, struct network_workspace *nw, int verbose);
void network_nse_integrate_traj(struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj);
void network_nse_integrate_traj_temp(struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj);

#endif /* NUCLEAR_NETWORK */

#endif /* NETWORK_NSE_H */
