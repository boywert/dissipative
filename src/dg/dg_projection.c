/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_projection.c
 * \date        02/2015
 * \author      Kevin Schaal
 * \brief       Projection functions
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */


#include "../allvars.h"
#include "../proto.h"
#include "dg_core_inline.h"

#ifdef DG

//higher order quadrature
#define NOF_QP_1D_MORE (DEGREE_K+3)

#ifdef TWODIMS
#define NOF_QP_MORE (NOF_QP_1D_MORE * NOF_QP_1D_MORE)
#else
#define NOF_QP_MORE (NOF_QP_1D_MORE * NOF_QP_1D_MORE * NOF_QP_1D_MORE)
#endif

static const int s_nof_qp_1d_more = NOF_QP_1D_MORE;
static double s_qp_1d_more[NOF_QP_1D_MORE];
static double s_qp_weights_1d_more[NOF_QP_1D_MORE];

static const int s_nof_qp_more = NOF_QP_MORE;
static double s_qp_x[NOF_QP_MORE];
static double s_qp_y[NOF_QP_MORE];
static double s_qp_z[NOF_QP_MORE];
static double s_qp_weight[NOF_QP_MORE];

/*!
 * Initialize the higher order quadrature
 */
void ini_higher_order_quadrature()
{
  int i, j, k;

  DG3D(int l;) for(i = 0; i < s_nof_qp_1d_more; i++)
    {
      s_qp_1d_more[s_nof_qp_1d_more - 1 - i] = legendre_root(s_nof_qp_1d_more, i + 1);
      s_qp_weights_1d_more[i] = gauss_weight(s_nof_qp_1d_more, s_qp_1d_more[s_nof_qp_1d_more - 1 - i]);
    }

  k = 0;

  for(i = 0; i < s_nof_qp_1d_more; i++)
    {
      for(j = 0; j < s_nof_qp_1d_more; j++)
        {
          DG3D(for(l = 0; l < s_nof_qp_1d_more; l++))
            {
              s_qp_x[k] = s_qp_1d_more[i];
              s_qp_y[k] = s_qp_1d_more[j];

#ifdef TWODIMS
              s_qp_z[k] = 0;
#else
              s_qp_z[k] = s_qp_1d_more[l];
#endif
              s_qp_weight[k] = s_qp_weights_1d_more[i] * s_qp_weights_1d_more[j] DG3D(*s_qp_weights_1d_more[l]);

              k++;
            }
        }
    }

  assert(k == s_nof_qp_more);
}


double integrate_density(int cell)
{
  int i;
  double x_cell, y_cell, z_cell;

  double w[5];
  double pv[5];
  double result = 0;

  for(i = 0; i < s_nof_qp_more; i++)
    {
      x_cell = s_qp_x[i];
      y_cell = s_qp_y[i];
      z_cell = s_qp_z[i];

      calc_state_at_pos(SphP[cell].Weights, x_cell, y_cell, z_cell, w);
      w_to_primvars(w, pv);

      result += pv[RHO] * s_qp_weight[i];
    }

  result *= SphP[cell].Volume / DG_PROJ_NORM;

  return result;

}

/*!
 * Calculate the L1 Norm Integral(|rho-rho0|)dV for the cell
 */
#ifdef DG_TEST_PROBLEM
double calc_L1_norm(int cell)
{
  int i;
  double x, y, z;

  double x_cell, y_cell, z_cell;

  double w[5];
  double norm = 0;

  for(i = 0; i < s_nof_qp_more; i++)
    {

      x_cell = s_qp_x[i];
      y_cell = s_qp_y[i];
      z_cell = s_qp_z[i];

      cell_coords_to_lab_coords(cell, x_cell, y_cell, z_cell, &x, &y, &z);

      calc_state_at_pos(SphP[cell].Weights, x_cell, y_cell, z_cell, w);

      norm += fabs(w[W_RHO] - ic_density(x, y, z)) * s_qp_weight[i];
    }

  norm *= SphP[cell].Volume / DG_PROJ_NORM;

  return norm;
}
#endif

/*!
 * \brief calculate the weights for the temperature polynomials
 * \param a_cell[base function][conserved variable]: the state of the cell
 * \param weights_out[base function][conserved varialbe]: output weights
 */
void calc_temperature_weights(double a_cell[NOF_BASE_FUNCTIONS][5], double weights_out[NOF_BASE_FUNCTIONS])
{
  double w_cell[5];
  double pv_cell[5];
  double temperature;

  int l, j;

  for(l = 0; l < Nof_base_functions; l++)       //loop over base functions
    {
      weights_out[l] = 0;

      for(j = 0; j < Nof_inner_quad_points; j++)        //loop over quadrature points
        {
          calc_state_at_quad_point(a_cell, j, w_cell);
          w_to_primvars(w_cell, pv_cell);

          temperature = pv_cell[PRES] / pv_cell[RHO];

          weights_out[l] += temperature * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
        }

      weights_out[l] /= DG_PROJ_NORM;
    }
}

/*!
 * \brief calculate the weights for the internal energy
 * \param a_cell[base function][conserved variable]: the state of the cell
 * \param weights_out[base function][conserved varialbe]: output weights
 */
void calc_u_weights(double a_cell[NOF_BASE_FUNCTIONS][5], double weights_out[NOF_BASE_FUNCTIONS])
{
  double w_cell[5];
  double pv_cell[5];
  double u;

  int l, j;

  for(l = 0; l < Nof_base_functions; l++)       //loop over base functions
    {
      weights_out[l] = 0;

      for(j = 0; j < Nof_inner_quad_points; j++)        //loop over quadrature points
        {
          calc_state_at_quad_point(a_cell, j, w_cell);
          w_to_primvars(w_cell, pv_cell);

          u = pv_cell[PRES] / (pv_cell[RHO] * (GAMMA - 1));

          weights_out[l] += u * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
        }

      weights_out[l] /= DG_PROJ_NORM;
    }
}





















#endif
