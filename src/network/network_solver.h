/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/network/network_solver.h
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

#ifndef NETWORK_SOLVER_H
#define NETWORK_SOLVER_H

#include "arepoconfig.h"

#if defined(NETWORK_SUPERLU) || defined(NETWORK_PARDISO)
#define NETWORK_SPARSE 1
#else
#define NETWORK_SPARSE 0
#endif

#include "network.h"
#define NETWORK_T_MAX 1.0e10

struct network_solver_data
{
  int nsteps;
  int maxstep;
  int *steps;
  int maxiter;
  int matrixsize;
  int matrixsize2;
  int nelements;
  double *aion;
  double tolerance;
  double *dy;
  double *ynew;
  double *yscale;
  network_var *rhs_deriv;
  double *first_rhs, *rhs;
  jacob_t jacob;
  jacob_t mod_jacob;
  double *x;
  double *err;
  double *qcol;
  double *a;
  double *alf;
#if NETWORK_VAR_RHO_T
  int iTemp;
#endif
#if NETWORK_VAR_RHO_T == NETWORK_VAR_RHO
  int iRho;
#endif
#ifdef NETWORK_OUTPUT
  FILE *fp;
#endif
};

struct network_solver_trajectory
{
  int ntimesteps, timestep;
  double *timesteps;
  double *rho;
  double *energy;
  double *x;
  double time, maxtime;
  double mintemp, maxtemp;
};

struct network_solver_data *network_solver_init(double tolerance, int matrixsize, int nelements, const struct network_nucdata *nucdata);
void network_solver_deinit(struct network_solver_data *nsd);
void network_solver_interpolate_trajectory(struct network_solver_trajectory *traj, double time, double *rho, double *energy);
void network_solver_integrate_traj(const struct network_data *nd, struct network_workspace *nw, struct network_solver_trajectory *traj);
void network_solver_integrate(double temp, double rho, double *y, double dt, const struct network_data *nd, struct network_workspace *nw);
#endif /* NETWORK_SOLVER_H */
