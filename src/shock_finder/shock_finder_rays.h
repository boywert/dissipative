/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/shock_finder_rays.h
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

#ifndef SHOCK_FINDER_RAYS_H
#define SHOCK_FINDER_RAYS_H

//Representation of a ray
typedef struct
{
  double startpos[3];           //!< the starting position of the ray
  double pos[3];                //!< the current position of the ray
  double dir[3];                //!< the direction of the ray
  double machnumber;            //!< the machnumber of the original cell
  double surface_value;         //!< the divergence of the surface_cell
  int current_cell;             //!< the cell in which the ray currently is
  int previous_cell;            //!< the cell from which the ray came from (needed for stability in find_next_voronoi_cell)
  int previous_cell_task;
  int original_cell;            //!< the cell from which the ray started
  int original_cell_task;
  int steps;                    //!< the number of steps the ray is on its way
  char preshock_dir;            //!< flags whether the ray goes currently in the pre- or postshock direction
  char quantities_set;          //!< 0: no quantities set, 1: post-shock quantity set, 2: pre-shock quantity set
  double rho_post_shock;        //!< post shock density
  double rho_pre_shock;         //!< the preshock density (needed for thermal energy flux calculation)
  double p_post_shock;          //!< post shock pressure
  double v_pre_shock[3];        //!< pre-shock velocity
  double v_post_shock[3];       //!< post shock velocity

#ifdef COSMIC_RAYS
  double pth_post_shock;       //!< post shock thermal pressure
  double pcr_post_shock;       //!< post shock cosmic ray pressure
  double pth_pre_shock;        //!< pre shock thermal pressure
  double pcr_pre_shock;        //!< pre shock cosmic ray pressure
  double eth_post_shock;       //!< post shock thermal energy density
  double ecr_post_shock;       //!< post shock cosmic ray energy density
  double gte;                  //!< generated thermal energy flux

#ifdef COSMIC_RAYS_SHOCK_ACCELERATION
#ifdef COSMIC_RAYS_MAGNETIC_OBLIQUITY
  double magnetic_obliquity_pre_shock; //!< the pre shock magnetic obliquity angle
#endif
  double eth_tot_pre_shock;      //!< pre shock thermal energy (total, no density)
  double post_shock_eth_tot_sum; //!< sum of the post shock thermal energies (including the shocked cell)
  int cell_count;                //!< number of elements written to post_shock_cells/post_shock_cells_tasks
  int post_shock_cells[3];       //!< indices of the post shock cells (shocked cell not included)
  int post_shock_cells_tasks[3]; //!< tasks of the post shock cells (shocked cell not included)
#endif
#else
  double post_shock_quant;      //!< the post shock quantity
  double c_pre_shock;           //!< the preshock sound speed (needed for thermal energy flux calculation)
  double p_pre_shock;           //!< pre shock pressure
  double t_pre_shock;           //!< pre-shock temperature
#endif
} shock_ray;

#endif
