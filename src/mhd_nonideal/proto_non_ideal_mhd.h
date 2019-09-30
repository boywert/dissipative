/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/MHD_non_ideal/proto_non_ideal_MHD.h
 * \date        04/2015
 * \author      Federico Marinacci
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef NON_IDEAL_MHD_H
#define NON_IDEAL_MHD_H

#ifdef NON_IDEAL_MHD
#ifdef MHD_POWELL
#ifdef AMBIPOLAR_DIFFUSION
void face_add_ambipolar_fluxes(const struct state *state_L, const struct state *state_R, const struct state *delta_time_L,
                           const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime,
                           double sqrtatime);
#endif
#ifdef OHMIC_DIFFUSION
void face_add_ohmic_fluxes(const struct state *state_L, const struct state *state_R, const struct state *delta_time_L,
                           const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime,
                           double sqrtatime);
#endif
#endif

#ifdef MHD_CT
#ifdef OHMIC_DIFFUSION
void face_add_ohmic_A_fluxes(int q, int qother, const struct state *st_L, const struct state *st_R, struct fluxes *flux,
                             double atime, double sqrtatime);
void face_add_ohmic_heating(const struct state *state_L, const struct state *state_R, const struct state *delta_time_L,
                            const struct state *delta_time_R, const struct geometry *geom, struct fluxes *flux, double atime, 
                            double sqrtatime);
#endif
#endif

#ifdef IMPLICIT_OHMIC_DIFFUSION

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

void diffusion_exchange_vector(int);
void update_B_and_Energy(void);
void init_ohm_conductivity(void);
void flux_exchange_vector(void);
void ohmic_diffusion(void);
double get_time_step(int *);
double get_subcycle_timestep(void);
void initialize_A_init(int);
void initialize_eta(int);
void calculate_unidirectional_fluxes_3D(int);
void set_initial_coeff(double, int);
void set_coeff(HYPRE_IJMatrix *, HYPRE_IJVector *, HYPRE_IJVector *, int *, double, int);
void prepare_linear_iterations(double, int);
void linear_iterations(double, int);
void get_xval(HYPRE_IJVector *, int *, int);
double square(double);
double mod(double);
double sign(double);
double MC(double, double);
double minmod(double, double);
double min(double, double);
double max(double, double);

#endif
#endif

#endif
