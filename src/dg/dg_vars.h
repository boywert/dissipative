/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_vars.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Variables for the DG module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef DG_VARS_H
#define DG_VARS_H

#include "dg_defines.h"

// three dimensional weights array: cells, base functions, variables
#define CBV(name) double (*name)[NOF_BASE_FUNCTIONS][5]

//#define INI_MIDPOINT_RULE
#define MIDPOINT_RULE_STEPS	20

/*
 * Time integration scheme
 */

#if (DG_ORDER == 1)
#define RK1
#define RK_STAGES 1

#elif (DG_ORDER == 2)
#define RK2
#define RK_STAGES 2

#elif (DG_ORDER == 3)
#define RK3
#define RK_STAGES 3

#else
#define RK4
#define RK_STAGES 5
#endif



/*
 * Verbose output
 */

#ifdef DG_VERBOSE
#define	DG_PRINTF(...) mpi_printf(__VA_ARGS__)
#else
#define	DG_PRINTF(...)
#endif

#define LEGENDRE_PRINTF(...) mpi_printf(__VA_ARGS__)
//#define LEGENDRE_PRINTF(...)


/*
 * Data access of the global arrays
 */

#if defined(SAVE_ACCESS) || defined(DG_DEBUG)

#define SET_Lobatto_points_1d(q,value)               set_lobatto_points_1d(q,value)
#define SET_Lobatto_points_weights_1d(q,value)       set_lobatto_points_weights_1d(q,value)
#define SET_Quad_points_1d(q,value)                  set_quad_points_1d(q,value)
#define SET_Quad_points_weights_1d(q,value)          set_quad_points_weights_1d(q,value)
#define SET_Union_quad_points_x(q,value)             set_union_quad_points_x(q,value)
#define SET_Union_quad_points_y(q,value)             set_union_quad_points_y(q,value)
#define SET_Union_quad_points_z(q,value)             set_union_quad_points_z(q,value)
#define SET_Inner_quad_points_x(q,value)             set_inner_quad_points_x(q,value)
#define SET_Inner_quad_points_y(q,value)             set_inner_quad_points_y(q,value)
#define SET_Inner_quad_points_z(q,value)             set_inner_quad_points_z(q,value)
#define SET_Inner_quad_points_weights(q,value)       set_inner_quad_points_weights(q,value)
#define SET_Outer_quad_points_x(e,q,value)	         set_outer_quad_points_x(e,q,value)
#define SET_Outer_quad_points_y(e,q,value)	         set_outer_quad_points_y(e,q,value)
#define SET_Outer_quad_points_z(e,q,value)           set_outer_quad_points_z(e,q,value)
#define SET_Outer_quad_points_weights(e,q,value)     set_outer_quad_points_weights(e,q,value)
#define SET_Union_base_values(q,l,value)             set_union_base_values(q,l,value)
#define SET_Inner_base_values(q,l,value)             set_inner_base_values(q,l,value)
#define SET_Inner_base_dx_values(q,l,value)          set_inner_base_dx_values(q,l,value)
#define SET_Inner_base_dy_values(q,l,value)          set_inner_base_dy_values(q,l,value)
#define SET_Inner_base_dz_values(q,l,value)          set_inner_base_dz_values(q,l,value)
#define SET_Outer_base_values(e,q,l,value)           set_outer_base_values(e,q,l,value)
#define SET_Outer_fine_base_values(e,q,l,value)      set_outer_fine_base_values(e,q,l,value)
#define SET_P_A(j,l,value)                           set_p_a(j,l,value)
#define SET_P_B(j,l,value)                           set_p_b(j,l,value)
#define SET_P_C(j,l,value)                           set_p_c(j,l,value)
#define SET_P_D(j,l,value)                           set_p_d(j,l,value)
#define SET_P_E(j,l,value)                           set_p_e(j,l,value)
#define SET_P_F(j,l,value)                           set_p_f(j,l,value)
#define SET_P_G(j,l,value)                           set_p_g(j,l,value)
#define SET_P_H(j,l,value)                           set_p_h(j,l,value)

#define GET_Lobatto_points_1d(q)                     get_lobatto_points_1d(q)
#define GET_Lobatto_points_weights_1d(q)             get_lobatto_points_weights_1d(q)
#define GET_Quad_points_1d(q)                        get_quad_points_1d(q)
#define GET_Quad_points_weights_1d(q)                get_quad_points_weights_1d(q)
#define GET_Union_quad_points_x(q)                   get_union_quad_points_x(q)
#define GET_Union_quad_points_y(q)                   get_union_quad_points_y(q)
#define GET_Union_quad_points_z(q)                   get_union_quad_points_z(q)
#define GET_Inner_quad_points_x(q)                   get_inner_quad_points_x(q)
#define GET_Inner_quad_points_y(q)                   get_inner_quad_points_y(q)
#define GET_Inner_quad_points_z(q)                   get_inner_quad_points_z(q)
#define GET_Inner_quad_points_weights(q)             get_inner_quad_points_weights(q)
#define GET_Outer_quad_points_x(e,q)                 get_outer_quad_points_x(e,q)
#define GET_Outer_quad_points_y(e,q)	               get_outer_quad_points_y(e,q)
#define GET_Outer_quad_points_z(e,q)                 get_outer_quad_points_z(e,q)
#define GET_Outer_quad_points_weights(e,q)           get_outer_quad_points_weights(e,q)
#define GET_Union_base_values(q,l)                   get_union_base_values(q,l)
#define GET_Inner_base_values(q,l)                   get_inner_base_values(q,l)
#define GET_Inner_base_dx_values(q,l)                get_inner_base_dx_values(q,l)
#define GET_Inner_base_dy_values(q,l)                get_inner_base_dy_values(q,l)
#define GET_Inner_base_dz_values(q,l)                get_inner_base_dz_values(q,l)
#define GET_Outer_base_values(e,q,l)                 get_outer_base_values(e,q,l)
#define GET_Outer_fine_base_values(e,q,l)            get_outer_fine_base_values(e,q,l)
#define GET_P_A(j,l)                                 get_p_a(j,l)
#define GET_P_B(j,l)                                 get_p_b(j,l)
#define GET_P_C(j,l)                                 get_p_c(j,l)
#define GET_P_D(j,l)                                 get_p_d(j,l)
#define GET_P_E(j,l)                                 get_p_e(j,l)
#define GET_P_F(j,l)                                 get_p_f(j,l)
#define GET_P_G(j,l)                                 get_p_g(j,l)
#define GET_P_H(j,l)                                 get_p_h(j,l)

#else

#define SET_Lobatto_points_1d(q,value)               Lobatto_points_1d[q]=value
#define SET_Lobatto_points_weights_1d(q,value)       Lobatto_points_weights_1d[q]=value
#define SET_Quad_points_1d(q,value)                  Quad_points_1d[q]=value
#define SET_Quad_points_weights_1d(q,value)          Quad_points_weights_1d[q]=value
#define SET_Union_quad_points_x(q,value)             Union_quad_points_x[q]=value
#define SET_Union_quad_points_y(q,value)             Union_quad_points_y[q]=value
#define SET_Union_quad_points_z(q,value)             Union_quad_points_z[q]=value
#define SET_Inner_quad_points_x(q,value)             Inner_quad_points_x[q]=value
#define SET_Inner_quad_points_y(q,value)             Inner_quad_points_y[q]=value
#define SET_Inner_quad_points_z(q,value)             Inner_quad_points_z[q]=value
#define SET_Inner_quad_points_weights(q,value)       Inner_quad_points_weights[q]=value
#define SET_Outer_quad_points_x(e,q,value)	         Outer_quad_points_x[e][q]=value
#define SET_Outer_quad_points_y(e,q,value)	         Outer_quad_points_y[e][q]=value
#define SET_Outer_quad_points_z(e,q,value)           Outer_quad_points_z[e][q]=value
#define SET_Outer_quad_points_weights(e,q,value)     Outer_quad_points_weights[e][q]=value
#define SET_Union_base_values(q,l,value)             Union_base_values[q][l]=value
#define SET_Inner_base_values(q,l,value)             Inner_base_values[q][l]=value
#define SET_Inner_base_dx_values(q,l,value)          Inner_base_dx_values[q][l]=value
#define SET_Inner_base_dy_values(q,l,value)          Inner_base_dy_values[q][l]=value
#define SET_Inner_base_dz_values(q,l,value)          Inner_base_dz_values[q][l]=value
#define SET_Outer_base_values(e,q,l,value)           Outer_base_values[e][q][l]=value
#define SET_Outer_fine_base_values(e,q,l,value)      Outer_fine_base_values[e][q][l]=value
#define SET_P_A(j,l,value)                           P_A[j][l]=value
#define SET_P_B(j,l,value)                           P_B[j][l]=value
#define SET_P_C(j,l,value)                           P_C[j][l]=value
#define SET_P_D(j,l,value)                           P_D[j][l]=value
#define SET_P_E(j,l,value)                           P_E[j][l]=value
#define SET_P_F(j,l,value)                           P_F[j][l]=value
#define SET_P_G(j,l,value)                           P_G[j][l]=value
#define SET_P_H(j,l,value)                           P_H[j][l]=value

#define GET_Lobatto_points_1d(q)                     Lobatto_points_1d[q]
#define GET_Lobatto_points_weights_1d(q)             Lobatto_points_weights_1d[q]
#define GET_Quad_points_1d(q)                        Quad_points_1d[q]
#define GET_Quad_points_weights_1d(q)                Quad_points_weights_1d[q]
#define GET_Union_quad_points_x(q)                   Union_quad_points_x[q]
#define GET_Union_quad_points_y(q)                   Union_quad_points_y[q]
#define GET_Union_quad_points_z(q)                   Union_quad_points_z[q]
#define GET_Inner_quad_points_x(q)                   Inner_quad_points_x[q]
#define GET_Inner_quad_points_y(q)                   Inner_quad_points_y[q]
#define GET_Inner_quad_points_z(q)                   Inner_quad_points_z[q]
#define GET_Inner_quad_points_weights(q)             Inner_quad_points_weights[q]
#define GET_Outer_quad_points_x(e,q)	               Outer_quad_points_x[e][q]
#define GET_Outer_quad_points_y(e,q)	               Outer_quad_points_y[e][q]
#define GET_Outer_quad_points_z(e,q)                 Outer_quad_points_z[e][q]
#define GET_Outer_quad_points_weights(e,q)           Outer_quad_points_weights[e][q]
#define GET_Union_base_values(q,l)                   Union_base_values[q][l]
#define GET_Inner_base_values(q,l)                   Inner_base_values[q][l]
#define GET_Inner_base_dx_values(q,l)                Inner_base_dx_values[q][l]
#define GET_Inner_base_dy_values(q,l)                Inner_base_dy_values[q][l]
#define GET_Inner_base_dz_values(q,l)                Inner_base_dz_values[q][l]
#define GET_Outer_base_values(e,q,l)                 Outer_base_values[e][q][l]
#define GET_Outer_fine_base_values(e,q,l)            Outer_fine_base_values[e][q][l]
#define GET_P_A(j,l)                                 P_A[j][l]
#define GET_P_B(j,l)                                 P_B[j][l]
#define GET_P_C(j,l)                                 P_C[j][l]
#define GET_P_D(j,l)                                 P_D[j][l]
#define GET_P_E(j,l)                                 P_E[j][l]
#define GET_P_F(j,l)                                 P_F[j][l]
#define GET_P_G(j,l)                                 P_G[j][l]
#define GET_P_H(j,l)                                 P_H[j][l]

#endif



enum E_interface {LEFT=0,RIGHT=1,BACK=2,FRONT=3,BOTTOM=4,TOP=5};
enum E_prim_var {RHO=0, PRES=1, VX=2, VY=3, VZ=4};
enum E_con_var {W_RHO=0, W_PX=1, W_PY=2, W_PZ=3, W_E=4};
enum E_ngb {LOCAL_CELL=0, LOCAL_NODE=1, GHOST_NODE=2, GHOST_CELL=3};

/*
 * global variables
 */


extern const int Nof_interfaces;
extern const int Nof_base_functions;
extern const int Nof_inner_quad_points;
extern const int Nof_outer_quad_points;
extern const int Nof_outer_fine_quad_points;
extern const int Nof_quad_points_1d;
extern int Nof_lobatto_points_1d;
extern int Nof_union_quad_points;

extern double* Union_quad_points_x;
extern double* Union_quad_points_y;
extern double* Union_quad_points_z;
extern double (*Union_base_values)[NOF_BASE_FUNCTIONS];

extern double* Lobatto_points_1d;
extern double* Lobatto_points_weights_1d;

extern double* Quad_points_1d;
extern double* Quad_points_weights_1d;

extern double* Inner_quad_points_x;
extern double* Inner_quad_points_y;
extern double* Inner_quad_points_z;
extern double* Inner_quad_points_weights;

extern double (*Outer_quad_points_x)[NOF_OUTER_QUAD_POINTS];
extern double (*Outer_quad_points_y)[NOF_OUTER_QUAD_POINTS];
extern double (*Outer_quad_points_z)[NOF_OUTER_QUAD_POINTS];
extern double (*Outer_quad_points_weights)[NOF_OUTER_QUAD_POINTS];

extern double (*Inner_base_values)[NOF_BASE_FUNCTIONS];
extern double (*Inner_base_dx_values)[NOF_BASE_FUNCTIONS];
extern double (*Inner_base_dy_values)[NOF_BASE_FUNCTIONS];
extern double (*Inner_base_dz_values)[NOF_BASE_FUNCTIONS];

extern double (*Outer_base_values)[NOF_OUTER_QUAD_POINTS][NOF_BASE_FUNCTIONS];
extern double (*Outer_fine_base_values)[NOF_OUTER_FINE_QUAD_POINTS][NOF_BASE_FUNCTIONS];

/*
 * Refinement / Derefinement Matrices
 */

extern double (*P_A)[NOF_BASE_FUNCTIONS];
extern double (*P_B)[NOF_BASE_FUNCTIONS];
extern double (*P_C)[NOF_BASE_FUNCTIONS];
extern double (*P_D)[NOF_BASE_FUNCTIONS];

#ifndef TWODIMS
extern double (*P_E)[NOF_BASE_FUNCTIONS];
extern double (*P_F)[NOF_BASE_FUNCTIONS];
extern double (*P_G)[NOF_BASE_FUNCTIONS];
extern double (*P_H)[NOF_BASE_FUNCTIONS];
#endif

#ifdef CHARACTERISTIC_LIMITER

/*
 * Flux Jacobian Eigenvector Matrices
 */

extern double (*Lx)[5];
extern double (*Ly)[5];
DG3D(extern double (*Lz)[5];)
extern double (*Rx)[5];
extern double (*Ry)[5];
DG3D(extern double (*Rz)[5];)

#endif

/*
 * Temporary storage of SphP weights
 */
extern double (*Temp_weights)[5];

/*
 * array for weights communication
 */

extern double (*Aexch)[NOF_BASE_FUNCTIONS][5];

/*
 * minimum values for pressure and density
 */
const extern double Epsilon_rho;
const extern double Epsilon_p;

/*
 * Additional text output file handels
 */
#ifdef DG_CON_VARS_SUM_TO_FILE
extern FILE* FdAngularMomentumDG;
extern FILE* FdMassDG;
extern FILE* FdEnergyDG;
#endif

extern FILE* FdInfoDG;

#endif
