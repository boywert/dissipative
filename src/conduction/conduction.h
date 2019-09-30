/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/conduction/conduction.h
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

#ifdef IMPLICIT_TI
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#endif



#ifdef MONOTONE_CONDUCTION
#ifndef IMPLICIT_TI
void conduction_exchange_vector(void) ;
void init_conductivity(void) ;
void flux_exchange_vector(void) ;
void monotone_conduction(void) ;
double get_time_step(void) ;
void initialize_kappa(void) ;
void calculate_unidirectional_fluxes_2D(void) ;
void calculate_unidirectional_fluxes_3D(void) ;
void calculate_incremental_u(double) ;
double square(double) ;
double mod(double) ;
double sign(double) ;
double MC(double, double) ;
double minmod(double, double) ;
double min(double, double) ;
double max(double, double) ;
#else
#ifndef SEMI_IMPLICIT_TI
void conduction_exchange_vector(void) ;
void U_exchange_vector(void) ;
void F2_exchange_vector(void) ;
void update_U_and_Energy(void) ;
void init_conductivity(void) ;
void flux_exchange_vector(void) ;
void monotone_conduction(void) ;
double get_time_step(void) ;
void initialize_kappa(void) ;
void calculate_unidirectional_fluxes_2D(void) ;
void calculate_unidirectional_fluxes_3D(void) ;
void nonlinear_iterations(double) ;
double calculate_nonlinear_error(double) ;
void calculate_F2(void) ;
double get_Utherm(int) ;
void set_initial_coeff(void) ;
double set_coeff(HYPRE_IJMatrix*, HYPRE_IJVector*, HYPRE_IJVector*, int*, double) ;
int calculate_cols(int) ;
double get_LF2(int, int, int, int, int, double*) ;
double linear_iterations(double) ;
void get_xval(HYPRE_IJVector*, int*) ;
double square(double) ;
double mod(double) ;
double sign(double) ;
double MC(double, double) ;
double minmod(double, double) ;
double min(double, double) ;
double max(double, double) ;
double compute_timestep(double, int, int) ;
#else
void conduction_exchange_vector(void) ;
void update_U_and_Energy(void) ;
void init_conductivity(void) ;
void flux_exchange_vector(void) ;
void monotone_conduction(void) ;
double get_time_step(int*) ;
double get_subcycle_timestep(void) ;
void initialize_u_init(void) ;
void initialize_kappa(void) ;
#ifdef TWODIMS
void calculate_unidirectional_fluxes_2D(void) ;
#else
void calculate_unidirectional_fluxes_3D(void) ;
#endif
void set_initial_coeff(double) ;
void set_coeff(HYPRE_IJMatrix*, HYPRE_IJVector*, HYPRE_IJVector*, int*, double) ;
void prepare_linear_iterations(double) ;
void linear_iterations(double) ;
void get_xval(HYPRE_IJVector*, int*) ;
double square(double) ;
double mod(double) ;
double sign(double) ;
double MC(double, double) ;
double minmod(double, double) ;
double min(double, double) ;
double max(double, double) ;
#endif
#endif
#endif
