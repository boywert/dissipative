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

//3d done

#include "../allvars.h"
#include "../proto.h"

#ifdef DG

const int Nof_interfaces = NOF_INTERFACES;      /*!< number of interfaces of a cell */
const int Nof_base_functions = NOF_BASE_FUNCTIONS;      /*!< number of base functions */
const int Nof_inner_quad_points = NOF_INNER_QUAD_POINTS;        /*!< number of quadrature points in the interior of a cell */
const int Nof_outer_quad_points = NOF_OUTER_QUAD_POINTS;        /*!< number of quadrature points per interface */
const int Nof_outer_fine_quad_points = NOF_OUTER_FINE_QUAD_POINTS;      /*!< number of fine quadrature points per interface for AMR boundary cells */
const int Nof_quad_points_1d = NOF_QUAD_POINTS_1D;      /*!< number of quadrature points in one dimension */

int Nof_lobatto_points_1d;      /*! < number of lobatto quadrature points in one dimension */
int Nof_union_quad_points;      /*! < number of tensored gauss and lobatto points */

double *Union_quad_points_x;    /*! < array of the x-coordinates of the union points */
double *Union_quad_points_y;    /*! < array of the y-coordinates of the union points */
double *Union_quad_points_z;    /*! < array of the z-coordinates of the union points */
double (*Union_base_values)[NOF_BASE_FUNCTIONS];        /*! < 2d array [union point][base function]: values of the base functions at the union points */

double *Lobatto_points_1d;      /*!< array of the position of the Lobatto points in one dimension */
double *Lobatto_points_weights_1d;      /*! < array of the weights of the 1d Lobatto points */

double *Quad_points_1d;         /*!< array of the position of the quadrature points in one dimension */
double *Quad_points_weights_1d; /*! < array of the weights of the 1d quadrature points */

double *Inner_quad_points_x;    /*!< array of the x-coordinates of the inner quadrature points */
double *Inner_quad_points_y;    /*!< array of the y-coordinates of the inner quadrature points */
double *Inner_quad_points_z;    /*!< array of the z-coordinates of the inner quadrature points */
double *Inner_quad_points_weights;      /*!< array of the quadrature weights of the inner quadrature points */

double (*Outer_quad_points_x)[NOF_OUTER_QUAD_POINTS];   /*!< 2d array [interface][quad point]: the x-coordinates of the outer quadrature points */
double (*Outer_quad_points_y)[NOF_OUTER_QUAD_POINTS];   /*!< 2d array [interface][quad point]: the y-coordinates of the outer quadrature points */
double (*Outer_quad_points_z)[NOF_OUTER_QUAD_POINTS];   /*!< 2d array [interface][quad point]: the z-coordinates of the outer quadrature points */
double (*Outer_quad_points_weights)[NOF_OUTER_QUAD_POINTS];     /*!< 2d array [interface][quad point]: the quadrature weights of the outer quadrature points */

double (*Inner_base_values)[NOF_BASE_FUNCTIONS];        /*!< 2d array [quad point][base function]: values of the base functions at the quadrature points */
double (*Inner_base_dx_values)[NOF_BASE_FUNCTIONS];     /*!< 2d array [quad point][base function]: values of the x-derivative of the base functions at the quadrature points */
double (*Inner_base_dy_values)[NOF_BASE_FUNCTIONS];     /*!< 2d array [quad point][base function]: values of the y-derivative of the base functions at the quadrature points */
double (*Inner_base_dz_values)[NOF_BASE_FUNCTIONS];     /*!< 2d array [quad point][base function]: values of the z-derivative of the base functions at the quadrature points */

double (*Outer_base_values)[NOF_OUTER_QUAD_POINTS][NOF_BASE_FUNCTIONS]; /*!< 3d array [interface][quad point][base function]: values of the base functions at the interface quadrature points */
double (*Outer_fine_base_values)[NOF_OUTER_FINE_QUAD_POINTS][NOF_BASE_FUNCTIONS];       /*!< 3d array [interface][fine quad point][base function]: values of the base functions at the interface part quadrature points */

double (*P_A)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_B)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_C)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_D)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */

#ifndef TWODIMS
double (*P_E)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_F)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_G)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
double (*P_H)[NOF_BASE_FUNCTIONS];      /*!< 2d array [base function][base function]: values of the refinement/derefinement matrix */
#endif

#ifdef CHARACTERISTIC_LIMITER
double (*Lx)[5];                /*!< 2d array: values of first left Eigenvector Matrix of the flux jacobian */
double (*Ly)[5];                /*!< 2d array: values of second left Eigenvector Matrix of the flux jacobian */
DG3D(double (*Lz)[5];)          /*!< 2d array: values of third left Eigenvector Matrix of the flux jacobian */
     double (*Rx)[5];           /*!< 2d array: values of first right Eigenvector Matrix of the flux jacobian */
     double (*Ry)[5];           /*!< 2d array: values of second reft Eigenvector Matrix of the flux jacobian */
DG3D(double (*Rz)[5];)          /*!< 2d array: values of third reft Eigenvector Matrix of the flux jacobian */
#endif
     double (*Temp_weights)[5]; /*!< 2d array [base function][variable]: temporary storage of SphP weights. */
     double (*Aexch)[NOF_BASE_FUNCTIONS][5];    /*!< Array [import cell][base function][variable] for imported weights */

/*
 * minimum values for pressure and density
 */
     const double Epsilon_rho = 1e-10;  /*!< Minimum mean density which is allowed. */
     const double Epsilon_p = 1e-10;    /*!< Minimum mean pressure which is allowed. */

#ifdef DG_CON_VARS_SUM_TO_FILE
     FILE *FdAngularMomentumDG; /*!< File handle for the angular momentum text file */
     FILE *FdMassDG;            /*!< File handle for the mass text file */
     FILE *FdEnergyDG;          /*!< File handle for the energy text file */
#endif

     FILE *FdInfoDG;            /*!< File handle for the info file */

#endif /* DG */
