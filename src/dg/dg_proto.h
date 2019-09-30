/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_proto.h
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Prototypes of DG functions
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef DG_PROTO_H
#define DG_PROTO_H

#include "dg_vars.h"

/*
 * dg_core.c
 */

void dg_initialize();
void dg_finalize();
void dg_set_initial_conditions();
void dg_compute_step();
void dg_setup_step();
void dg_end_step();

void delete_duplicated_union_quad_points();
void calc_R_inner(CBV(a), int cell, double time, double dt, double R_inner[NOF_BASE_FUNCTIONS][5], double rk_weight);
enum E_ngb get_neighbour_type(int index_ngb);
void get_neighbour_weights(int dp_center, enum E_interface neighbour, int base_fct, CBV(a), double w[5]);
void calc_state_at_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int quad_point, double w[5]);
void calc_state_at_union_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int quad_point, double w[5]);
void calc_state_at_outer_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int interface, int quad_point, double w[5], double (*P_X)[NOF_BASE_FUNCTIONS]);
void calc_state_at_outer_fine_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int interface, int fine_quad_point, double w[5]);
void calc_state_at_pos(double a_cell[NOF_BASE_FUNCTIONS][5], double x, double y, double z, double w[5]);
void flux1(double w[5], double f[5]);
void flux2(double w[5], double f[5]);
void flux3(double w[5], double f[5]);

void copy_cell_weights_to_array(CBV(w));
void copy_array_to_cell_weights(CBV(w));
void copy_array_to_array(CBV(a),CBV(a_copy));
void exchange_weights(CBV(a));
void assign_projection_matrix(int dp_high_lvl, int dp_low_lvl, double (**P_X)[NOF_BASE_FUNCTIONS]);

void dg_update_conserved_variables();
double dg_time_step_cell(int cell);
double dg_time_step_cell_gravity(int cell);
double angular_momentum(double a_cell[NOF_BASE_FUNCTIONS][5], int cell);
double spin(double a_cell[NOF_BASE_FUNCTIONS][5], int cell);

int opposite_interface(int e);

void multiply_matrix_matrix(double m1[5][5],double m2[5][5],double m_out[5][5]);
void multiply_matrix_vector(double m[5][5],double v[5],double v_out[5]);
void multiply_matrix_vector_dims(double* m,double* v,double* v_out, int dims);

void solve_quadratic_equation(double a, double b, double c, double *solution_plus, double *solution_minus);

/*
 * dg_fluxes.c
 */

void subtract_R_outer(tessellation *T, double c, CBV(a1), CBV(a2));
int dg_face_get_state(tessellation * T, CBV(a), int interface, int q, int p, int i, struct state *st, double (*P_X)[NOF_BASE_FUNCTIONS]);
void dg_apply_flux_list(CBV(a));

/*
 * dg_limiter.c
 */

double minmod(double a, double b, double c);
double minmodB(double a, double b, double c, double h);
void minmod_limiter(CBV(a));
void positivity_limiter(double a_cell[NOF_BASE_FUNCTIONS][5], int cell);

/*
 *   dg_legendre.c
 */

double P_0(double x);
double P_1(double x);
double P_2(double x);
double P_3(double x);
double P_4(double x);
double P_5(double x);
double dP_0(double x);
double dP_1(double x);
double dP_2(double x);
double dP_3(double x);
double dP_4(double x);
double dP_5(double x);

double Legendre(int n, double x);
double scaled_Legendre(int n, double x);
double deriv_Legendre(int n, double x);
double deriv_scaled_Legendre(int n, double x);
double legendre_root(int n, int k);
double gauss_weight(int n, double root);
void index_to_base_function(int k, int *Px, int *Py, int *Pz);
double base_function(int k, double x, double y, double z);
double deriv_x_base_function(int k, double x, double y, double z);
double deriv_y_base_function(int k, double x, double y, double z);
double deriv_z_base_function(int k, double x, double y, double z);

void abs_pos_inner_quad_points(int q, double cell_dl, double cell_center[3], double pos_abs[3]);
void abs_pos_outer_quad_points(int e, int q, double cell_dl, double cell_center[3], double pos_abs[3]);

void ini_1d_quadrature_points();
void ini_multi_d_quadrature_points();
void ini_base_function_values();
void ini_refinement_matrices();

void print_quad_info();
void print_base_info();
void print_ref_matrix_info();

/*
 * dg_time_integration.c
 */

//void dg_compute_step();
void print_time_integration_info();
void ini_rk_coefficients();

/*
 * dg_recompute.c
 */

void dg_recompute_nodes(void);
void dg_recompute_nodes_normal(void);
void dg_recompute_node_recursive(int no, int mode);
void dg_exchange_topleafdata(void);

/*
 * dg_projection.c
 */
void ini_higher_order_quadrature();
double calc_L1_norm(int cell);
double integrate_density(int cell);
void calc_temperature_weights(double a_cell[NOF_BASE_FUNCTIONS][5], double weights_out[NOF_BASE_FUNCTIONS]);
void calc_u_weights(double a_cell[NOF_BASE_FUNCTIONS][5], double weights_out[NOF_BASE_FUNCTIONS]);

/*
 * dg_set_get.c
 */
void set_lobatto_points_1d(int q,double value);
void set_lobatto_points_weights_1d(int q,double value);
void set_quad_points_1d(int q,double value);
void set_quad_points_weights_1d(int q,double value);
void set_union_quad_points_x(int q,double value);
void set_union_quad_points_y(int q,double value);
void set_union_quad_points_z(int q,double value);
void set_inner_quad_points_x(int q,double value);
void set_inner_quad_points_y(int q,double value);
void set_inner_quad_points_z(int q,double value);
void set_inner_quad_points_weights(int q,double value);
void set_outer_quad_points_x(int e,int q,double value);
void set_outer_quad_points_y(int e,int q,double value);
void set_outer_quad_points_z(int e,int q,double value);
void set_outer_quad_points_weights(int e,int q,double value);
void set_union_base_values(int q,int l,double value);
void set_inner_base_values(int q,int l,double value);
void set_inner_base_dx_values(int q,int l,double value);
void set_inner_base_dy_values(int q,int l,double value);
void set_inner_base_dz_values(int q,int l,double value);
void set_outer_base_values(int e,int q, int l,double value);
void set_outer_fine_base_values(int i,int q, int l,double value);
void set_p_a(int j,int l,double value);
void set_p_b(int j,int l,double value);
void set_p_c(int j,int l,double value);
void set_p_d(int j,int l,double value);

double get_lobatto_points_1d(int q);
double get_lobatto_points_weights_1d(int q);
double get_quad_points_1d(int q);
double get_quad_points_weights_1d(int q);
double get_union_quad_points_x(int q);
double get_union_quad_points_y(int q);
double get_union_quad_points_z(int q);
double get_inner_quad_points_x(int q);
double get_inner_quad_points_y(int q);
double get_inner_quad_points_z(int q);
double get_inner_quad_points_weights(int q);
double get_outer_quad_points_x(int e,int q);
double get_outer_quad_points_y(int e,int q);
double get_outer_quad_points_z(int e,int q);
double get_outer_quad_points_weights(int e,int q);
double get_union_base_values(int q, int l);
double get_inner_base_values(int q, int l);
double get_inner_base_dx_values(int q, int l);
double get_inner_base_dy_values(int q, int l);
double get_inner_base_dz_values(int q, int l);
double get_outer_base_values(int e, int q, int l);
double get_outer_fine_base_values(int i, int q, int l);
double get_p_a(int j,int l);
double get_p_b(int j,int l);
double get_p_c(int j,int l);
double get_p_d(int j,int l);

/*
 * dg_debug.c
 */


void cell_info_a(const char* location, CBV(a));
void cell_info(const char* location);
void print_outer_quad_points_weights();
void print_weights(int cell);
void print_imported_weights();
void assert_imported_weights_positivity(const char* error);
void print_timebins();
void print_weights_array(CBV(a));
void print_weights_a(int cell, CBV(a));
void assert_weights_positivity(CBV(a), const char* error);
void assert_sphp_weights_positivity(const char* error);
void assert_state_positivity(struct state *s, const char* error);
int is_id(double m[5][5]);
void check_time_steps();

/*
 * External Acceleration 
 */

#ifdef DG_EXTERNAL_ACCELERATION
void dg_acceleration(double x, double y, double z, double* acc);
#endif


/*
 * dg_io.c
 */
#ifdef DG
void io_func_dgw0(int particle, int components, void* buffer, int mode);
void io_func_dgw1(int particle, int components, void* buffer, int mode);
void io_func_dgw2(int particle, int components, void* buffer, int mode);
void io_func_dgw3(int particle, int components, void* buffer, int mode);
void io_func_dgw4(int particle, int components, void* buffer, int mode);
#ifdef DG_SET_IC_FROM_AVERAGES
void load_weights_from_averages();
#endif
#endif
#ifdef OUTPUT_DG_ACCELERATION
void io_func_dg_accel(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_ANGULAR_MOMENTUM
void io_func_dg_angular_momentum(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_SPIN
void io_func_dg_spin(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_TIMESTEP
void io_func_dg_timestep(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_TEMPERATURE
void io_func_dgw_temperature(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_U
void io_func_dgw_u(int particle, int components, void* buffer, int mode);
#endif
#ifdef OUTPUT_DG_L1_NORM
void io_func_dg_norm(int particle, int components, void* buffer, int mode);
#endif

/*
 * dg_refinement.c
 */

#ifdef DG
void dg_split_hydro(int sun_num, int sun_index);
void dg_sum_hydro(int new_cell, int sun_num, int sun_index);
#endif
#endif
