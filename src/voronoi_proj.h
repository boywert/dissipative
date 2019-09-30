/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/voronoi_proj.h
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

#ifndef VORONOI_PROJ_H
#define VORONOI_PROJ_H

#include "arepoconfig.h"
#include "allvars.h"

#ifdef VORONOI_PROJ
struct projection_data
{
  double center[3];
  double dir[3];
  double normal1[3];
  double normal2[3];

  int nx, ny, npixels, pflag;
  double boxx, boxy, boxz;
  double dx, dy;

  double *rho;
#ifdef MHD
  double *bsqr;
#endif
  double *ptot;
#ifdef GFM_DUST
  double *metal_rho, *dust_rho;
#endif
#ifdef VORONOI_PROJ_TAU
  double *tau, *oldk, *tps, *rhops, *len;
  int *finished;
#endif

  /* Rays */
  int Nray, MaxNray;
  ray_data *Ray;
};

int voronoi_proj_setup(int argc, char **argv, struct projection_data *pdata);
int voronoi_proj_init(struct projection_data *pdata, int flag);
int voronoi_proj_run(struct projection_data *pdata);
int voronoi_proj_deinit(struct projection_data *pdata);
#endif

#endif /* VORONOI_PROJ_H */
