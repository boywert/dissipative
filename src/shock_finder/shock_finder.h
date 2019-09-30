/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/shock_finder/shock_finder.h
 * \date        2013
 * \author      Kevin Schaal
 * \brief       A (post-processing) shock finder for Arepo
 * \details
 *
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef SHOCK_FINDER_H
#define SHOCK_FINDER_H

#include <stdio.h>
#include "shock_finder_fields.h"
#include "shock_finder_rays.h"

//filtering of the starforming region (effective equation of state)
#ifdef USE_SFR
#if !defined(FM_SFR) && !defined(ISM) && !defined(LOCAL_FEEDBACK) && !defined(QUICK_LYALPHA)
#define SHOCK_FINDER_FILTER_SFR
#endif
#endif

//filtering of inconsistent jumps
#if (defined(COOLING) || defined(TGCHEM))
#define SHOCK_FINDER_FILTER_INCONSISTENT
#endif

//post-processing shock finder
#ifdef SHOCK_FINDER_POST_PROCESSING
#define SDATA(i)  SData[i]
#define SDATAEXCH(i) SDataExch[i]
#endif

//shock finding before output
#ifdef SHOCK_FINDER_BEFORE_OUTPUT
#define SHOCK_FINDER_AREPO
#ifdef COSMIC_RAYS
#define ZONE_JUMP_CR
#define SHOCK_JUMP_CR
#define SHOCK_DIR_GRAD_T_CR
#else
#define ZONE_JUMP_P
#define ZONE_JUMP_T
#define SHOCK_DIR_GRAD_T
#define SHOCK_JUMP_T
#endif
#ifdef TWODIMS
#define SURFACE_SPHERE_APPROX
#else
#define SURFACE_ANGLE_APPROX
#endif
#define RESET_WRONG_JUMPS
#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_ON_THE_FLY)
#error Compile only with one of the flags SHOCK_FINDER_POST_PROCESSING / SHOCK_FINDER_BEFORE_OUTPUT / SHOCK_FINDER_ON_THE_FLY
#endif
#define SDATA(i)  SData[i]
#define SDATAEXCH(i) SDataExch[i]
#endif

//shock finder on-the-fly
#ifdef SHOCK_FINDER_ON_THE_FLY
#define SHOCK_FINDER_AREPO
#ifdef COSMIC_RAYS
#define ZONE_JUMP_CR
#define SHOCK_JUMP_CR
#define SHOCK_DIR_GRAD_T_CR
#else
#define ZONE_JUMP_P
#define ZONE_JUMP_T
#define SHOCK_DIR_GRAD_T
#define SHOCK_JUMP_T
#endif
#ifdef VORONOI_STATIC_MESH
#define SURFACE_CUBE_APPROX
#else
#ifdef TWODIMS
#define SURFACE_SPHERE_APPROX
#else
#define SURFACE_ANGLE_APPROX
#endif
#endif
#define RESET_WRONG_JUMPS
#if defined(SHOCK_FINDER_POST_PROCESSING) || defined(SHOCK_FINDER_BEFORE_OUTPUT)
#error Compile only with one of the flags SHOCK_FINDER_POST_PROCESSING / SHOCK_FINDER_BEFORE_OUTPUT / SHOCK_FINDER_ON_THE_FLY
#endif
#define SDATA(i)  SphP[i]
#define SDATAEXCH(i) PrimExch[i]
#endif

//debug
#ifdef SHOCK_FINDER_VERBOSE2
#define VPRINTF(...) mpi_printf(__VA_ARGS__)
#else
#define VPRINTF(...)
#endif


//write shock direction in hdf5 output
#define SHOCK_DIR_IN_OUTPUT

#ifdef SHOCK_DIR_IN_OUTPUT
#define SDIR(...) __VA_ARGS__
#else
#define SDIR(...)
#endif

//write maximum face angle in hdf5 output
#ifdef MAX_FACE_ANGLE_IN_OUTPUT
#define FANGLE(...) __VA_ARGS__
#else
#define FANGLE(...)
#endif

//shock finder global variables
struct ShockFinderVariables
{
  FILE *pFile;                  //!< file available for additional output (for each task)
  long long int nof_gas_part_total;
  long long int nof_shockzone_part_total;
  long long int nof_shocksurface_part_total;

  long long int nof_shockzone_part_local;
  long long int nof_shocksurface_part_local;

  //flags for shock finder before output
  int limit_gradients;

  //memory reservation for shock finder output
  float *coord_buffer;
  float *velocity_buffer;
  SDIR(float *shockdir_buffer;)
  float *vpre_buffer;
  float *vpost_buffer;
  FANGLE(float *fangle_buffer;)
  float *machnum_buffer;
  float *surface_buffer;
  float *volume_buffer;
  float *density_buffer;
  float *genUflux_buffer;
  float *cpre_buffer;
  float *rhopre_buffer;
  float *ppre_buffer;
  float *rhopost_buffer;
  float *ppost_buffer;
  float *tpre_buffer;
  float *intenergy_buffer;
  float *temperature_buffer;
  int *zoneflag_buffer;
  unsigned long long int *id_buffer;
};

extern struct ShockFinderVariables SfVars;

#ifdef SHOCK_FINDER_POST_PROCESSING
void shock_finder();
void shock_finder_output();
#endif

void shock_finder_on_the_fly();
void exchange_shock_data();
void calculate_additional_shock_data();
void print_ini_info();
void reset_shock_data();
void shock_finder_alloc_memory();
void shock_finder_free_memory();
void shock_finder_free_storage();
void shock_finder_arepo();
void shock_finder_ryu();
void shock_finder_skillman();
double distance_to_border_x(int cell);
double distance_to_border_y(int cell);
double distance_to_border_z(int cell);

/*!
 *  Shock finder data structure for shock finder before output and post-processing shock finder.
 *  In this case the fields are not stored on SphP, in order to save memory.
 */
typedef struct
{
  SHOCK_FINDER_FIELDS
#ifdef COSMIC_RAYS
  SHOCK_FINDER_CR_FIELDS
#else
  SHOCK_FINDER_IDEAL_HYDRO_FIELDS
#endif
} ShockData;

#ifndef SHOCK_FINDER_ON_THE_FLY
extern ShockData *SData;
extern ShockData *SDataExch;
#endif

//helper functions
//

void calculate_temperatures();
void calculate_gradT();
void calculate_divvel();

double shock_surface_area(int k);

double delta_M(double M);
double generated_thermal_energy_flux(int k);
double dissipated_energy_in_erg_per_s(int k);

double calculate_machnumber_T_jump(double T_preshock, double T_postshock);
double calculate_machnumber_rho_jump(double rho_preshock, double rho_postshock);
double calculate_machnumber_p_jump(double p_preshock, double p_postshock);
double calculate_machnumber_vdiff(double dv, double c_preshock);
double calculate_machnumber_S_jump(double S_preshock, double S_postshock, double guess);

//Mach number and energy dissipation for cosmic rays
#ifdef COSMIC_RAYS
double calculate_machnumber_cosmic_rays(double pth1, double pth2, double pcr1, double pcr2, double rho1, double rho2);
double calculate_machnumber_cosmic_rays_zone(double pth1,double pth2,double pcr1,double pcr2,double gamma1,double gamma2);
double calculate_machnumber_cosmic_rays_from_energy_equation(double pth1, double pth2, double pcr1, double pcr2, double rho1, double rho2);
double generated_thermal_energy_flux_cosmic_rays(double eth1, double eth2, double ecr1, double ecr2, double rho1, double rho2, double M, double c1);
double calculate_gamma_effective(double pcr, double pth);
double sound_speed_cosmic_rays(double pcr1, double pth1, double rho1);
double magnetic_shock_obliquity(double* sdir, double* bfield);
double injection_weight(int cell);
double injection_weight_exch(int cell);
#endif

//helper functions
int find_neighbouring_cell(int k, double dir[3]);
int find_neighbouring_cell_pos(int k, double start_point[3], double dir[3], int previous, double *end_point);
int find_neighbouring_cell_para(int k, double dir[3], int *task, int *original_index, int *dp_index, int *boundary_reached);
int find_neighbouring_cell_pos_para(int k, double start_point[3], double dir[3], int previous, int task_of_previous, double *end_point, int *task, int *original_index, int *edge, int *boundary_reached);

void solve_quadratic_equation(double a, double b, double c, double *solution_plus, double *solution_minus);
void solve_quadratic_equation_pq(double p, double q, double *solution_plus, double *solution_minus);

//calculate the divergence of a cell
static double inline divergence(int k)
{
#ifdef LOCAL_DIVERGENCE
  return SphP[k].Grad.dvel[0][0] + SphP[k].Grad.dvel[1][1] + SphP[k].Grad.dvel[2][2];
#else
  return SDATA(k).Divvel;
#endif
}

//calculate dhe divergence of a cell if it is imported
static double inline divergence_exchange(int k)
{
#ifdef LOCAL_DIVERGENCE
  return GradExch[k].dvel[0][0] + GradExch[k].dvel[1][1] + GradExch[k].dvel[2][2];
#else
  return SDATAEXCH(k).Divvel;
#endif
}


#endif
