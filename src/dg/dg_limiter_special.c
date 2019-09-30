/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_limiter_special.c
 * \date        10/2014
 * \author      Kevin Schaal
 * \brief       Functions for limiting the higher order terms
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"


//Before using functions from this file adapt them to the 3d code version!
#ifdef DG_SPECIAL_LIMITERS

#ifdef DISCONTINUITY_DETECTION
inline static int detect_discontinuity(CBV(a), int cell);
#endif

#ifdef MACHNUM_JUMP_DETECTION
inline static double calculate_mach_number(CBV(a), int cell);
static double calculate_machnumber_T_jump(double T_jump);
static double calculate_machnumber_rho_jump(double rho_jump);
static double calculate_T_jump(double M);
#endif


#ifdef ANGULAR_MOMENTUM_LIMITER
static void angular_momentum_limiter(CBV(a), int cell, double target_spin);
#endif

#ifdef NEIGHBOUR_BOUND
static void min_max_conserved_values_on_interface(double a_cell[NOF_BASE_FUNCTIONS][4], int cell, double w_max[4], double w_min[4])
#endif
#ifdef NEIGHBOUR_BOUND
     static void min_max_conserved_values_on_interface(double a_cell[NOF_BASE_FUNCTIONS][4], int cell, double w_max[4], double w_min[4])
{
  int e, q, k;
  double w[4];

  w_max[0] = -DBL_MAX;
  w_max[1] = -DBL_MAX;
  w_max[2] = -DBL_MAX;
  w_max[3] = -DBL_MAX;

  w_min[0] = DBL_MAX;
  w_min[1] = DBL_MAX;
  w_min[2] = DBL_MAX;
  w_min[3] = DBL_MAX;

  for(e = 0; e < 4; e++)
    {
      for(q = 0; q < Nof_outer_quad_points; q++)
        {
          calc_state_at_outer_quad_point(a_cell, e, q, w, 0);

          for(k = 0; k < 4; k++)
            {
              w_max[k] = fmax(w[k], w_max[k]);
              w_min[k] = fmin(w[k], w_min[k]);
            }
        }
    }
}

static void min_max_mean_conserved_values_of_neighbours(double w0_left[4], double w0_right[4], double w0_back[4], double w0_front[4], double w_max[4], double w_min[4])
{
  int k;

  for(k = 0; k < 4; k++)
    {
      w_max[k] = fmax(fmax(fmax(w0_left[k], w0_right[k]), w0_back[k]), w0_front[k]);
      w_min[k] = fmin(fmin(fmin(w0_left[k], w0_right[k]), w0_back[k]), w0_front[k]);
    }
}
#endif


#ifdef ANGULAR_MOMENTUM_LIMITER
static void angular_momentum_limiter(CBV(a), int cell, double target_spin)
{
  DG_PRINTF("Angular Momentum Limiter\n");

  //constant factors
  double dl = amr_length[Mesh.DP[cell].level];
  double c = 0.5 / sqrt(3.) * dl * dl * dl;

  //limited slopes
  double a1_tilde = a[cell][1][2];
  double a2_tilde = a[cell][2][1];

#ifdef PRAVEENS_VERSION
  //recover the initial angular momentum by adjusting the velocity slopes

  a[cell][1][2] = 0.5 * (a1_tilde + a2_tilde + target_spin / c);
  a[cell][2][1] = a[cell][1][2] - target_spin / c;

  assert(fabs(target_spin - spin(a[cell], cell) < 1e-10));

#else
  //improve angular momentum conservation by limiting the velocity slopes further

  //spin after limiting
  double limited_spin = spin(a[cell], cell);

  //best slopes
  double a1_best = a1_tilde;
  double a2_best = a2_tilde;

  double delta_spin_best = fabs(target_spin - limited_spin);

  //temporary slopes
  double a1_temp;
  double a2_temp;
  double spin_temp;
  double delta_spin_temp;

  //slope arrays
  double a1_array[4];
  double a2_array[4];

  //solutions along the boundary lines of the allowed R^2 region
  a1_array[0] = 0;
  a2_array[0] = -target_spin / c;

  a1_array[1] = target_spin / c;
  a2_array[1] = 0;

  a1_array[2] = a1_tilde;
  a2_array[2] = a1_tilde - target_spin / c;

  a1_array[3] = a2_tilde + target_spin / c;
  a2_array[3] = a2_tilde;

  //check for minimum
  int i;

  for(i = 0; i < 4; i++)
    {
      a1_temp = a1_array[i];
      a2_temp = a2_array[i];

      if(!(a1_temp * a1_tilde >= 0 && fabs(a1_temp) <= fabs(a1_tilde)))
        {
          a1_temp = 0;

          spin_temp = c * (a1_temp - a2_temp);
          //printf("\t\ta1_temp:%e, a2_temp:%e, spin_temp:%e\n", a1_temp, a2_temp, spin_temp);
          delta_spin_temp = fabs(spin_temp - target_spin);
          if(delta_spin_temp < delta_spin_best)
            {
              delta_spin_best = delta_spin_temp;
              a1_best = a1_temp;
              a2_best = a2_temp;
            }

          a1_temp = a1_tilde;

          spin_temp = c * (a1_temp - a2_temp);
          //printf("\t\ta1_temp:%e, a2_temp:%e, spin_temp:%e\n", a1_temp, a2_temp, spin_temp);
          delta_spin_temp = fabs(spin_temp - target_spin);
          if(delta_spin_temp < delta_spin_best)
            {
              delta_spin_best = delta_spin_temp;
              a1_best = a1_temp;
              a2_best = a2_temp;
            }
        }
      else if(!(a2_temp * a2_tilde >= 0 && fabs(a2_temp) <= fabs(a2_tilde)))
        {
          a2_temp = 0;

          spin_temp = c * (a1_temp - a2_temp);
          //printf("\t\ta1_temp:%e, a2_temp:%e, spin_temp:%e\n", a1_temp, a2_temp, spin_temp);
          delta_spin_temp = fabs(spin_temp - target_spin);
          if(delta_spin_temp < delta_spin_best)
            {
              delta_spin_best = delta_spin_temp;
              a1_best = a1_temp;
              a2_best = a2_temp;
            }

          a2_temp = a2_tilde;

          spin_temp = c * (a1_temp - a2_temp);
          //printf("\t\ta1_temp:%e, a2_temp:%e, spin_temp:%e\n", a1_temp, a2_temp, spin_temp);
          delta_spin_temp = fabs(spin_temp - target_spin);
          if(delta_spin_temp < delta_spin_best)
            {
              delta_spin_best = delta_spin_temp;
              a1_best = a1_temp;
              a2_best = a2_temp;
            }
        }
      else                      //calculated solution is in the allowed range
        {
          //printf("solution is in allowed range!\n");
          spin_temp = c * (a1_temp - a2_temp);
          //printf("\t\ta1_temp:%e, a2_temp:%e, spin_temp:%e\n", a1_temp, a2_temp, spin_temp);
          delta_spin_temp = fabs(spin_temp - target_spin);
          if(delta_spin_temp < delta_spin_best)
            {
              delta_spin_best = delta_spin_temp;
              a1_best = a1_temp;
              a2_best = a2_temp;
            }
        }
    }                           //end loop

  a[cell][1][2] = a1_best;
  a[cell][2][1] = a2_best;

  assert(a1_best * a1_tilde >= 0 && fabs(a1_best) <= fabs(a1_tilde));
  assert(a2_best * a2_tilde >= 0 && fabs(a2_best) <= fabs(a2_tilde));

  //printf("cell: %d, old slopes:%e, %e, new slopes: %e, %e\n",cell, a1_tilde, a2_tilde, a1_best, a2_best);
  //printf("\tlimited_spin:%e\n\t new_spin: %e\n\t target spin:  %e\n\n\n", limited_spin, spin(a[cell],cell), target_spin);

#endif
}
#endif

#ifdef DISCONTINUITY_DETECTION
/*!
 * Check whether a discontinuity exists between this cell and the neighbouring cells
 */
inline static int detect_discontinuity(CBV(a), int cell)
{
  DG_PRINTF("Discontinuity Detection\n");

  int e, eo, q;
  int cell_ngb;
  double rho_max;
  double sum;
  double T;
  CBV(a_general);
  int cell_ngb_general;

  //states
  double w1[4];
  double w2[4];

  //cell length
  double dl = amr_length[Mesh.DP[cell].level];

  rho_max = 0;
  T = 0;
  sum = 0;

  for(e = 0; e < 4; e++)        //loop over edges
    {

      if(!(SphP[cell].Inflow_boundaries & (1 << e)))
        continue;               //edge is no inflow boundary

      eo = opposite_interface(e);
      cell_ngb = Mesh.DP[cell].neighbors[e];

      if(cell_ngb < Ngb_MaxPart)        //local cell
        {
          cell_ngb_general = cell_ngb;
          a_general = a;
        }
      else                      //ghost cell
        {
          cell_ngb_general = Mesh.DP[cell_ngb - Mesh.nodes_total].index;
          a_general = Aexch;
        }

      //integrate density difference along the edge
      for(q = 0; q < Nof_outer_quad_points; q++)
        {
          calc_state_at_outer_quad_point(a[cell], e, q, w1, 0);
          calc_state_at_outer_quad_point(a_general[cell_ngb_general], eo, q, w2, 0);

          if(w1[0] > rho_max)
            {
              rho_max = w1[0];
            }

          sum += 0.5 * (w1[0] - w2[0]) * GET_Outer_quad_points_weights(e, q);
        }
    }

  if(rho_max == 0)              //no inflow boundaries
    {
      T = 0;
    }
  else
    {
      T = fabs(sum) / rho_max / pow(dl, 1 + All.DG_alpha * (DEGREE_K + 1));
    }

#ifdef OUTPUT_DG_DISCONTINUITIES
  SphP[cell].Discontinuity = T;
#endif

  return (T > 1);               // T>1 => discontinuity found
}
#endif

#ifdef MACHNUM_JUMP_DETECTION
/*!
 * Calculate the largest discontinuity jump relative to neighbouring cells
 */
inline static double calculate_mach_number(CBV(a), int cell)
{
  int e, eo, q;
  int cell_ngb;

  CBV(a_general);
  int cell_ngb_general;

  //states
  double w1[4];
  double w2[4];
  double pv1[4];
  double pv2[4];

  double jump;
  double jump_max = 0;

  for(e = 0; e < 4; e++)        //loop over edges
    {
      eo = opposite_interface(e);
      cell_ngb = Mesh.DP[cell].neighbors[e];

      if(cell_ngb < Ngb_MaxPart)        //local cell
        {
          cell_ngb_general = cell_ngb;
          a_general = a;
        }
      else                      //ghost cell
        {
          cell_ngb_general = Mesh.DP[cell_ngb - Mesh.nodes_total].index;
          a_general = Aexch;
        }

      //find largest temperautre jump among quadrature points
      for(q = 0; q < Nof_outer_quad_points; q++)
        {
          calc_state_at_outer_quad_point(a[cell], e, q, w1, 0);
          calc_state_at_outer_quad_point(a_general[cell_ngb_general], eo, q, w2, 0);

          w_to_primvars(w1, pv1);
          w_to_primvars(w2, pv2);

          if(w1[0] < 0 || w2[0] < 0)    //can happen since the solution is not yet positivity limited
            {
              return 3;
            }
          else if(w1[0] > w2[0])
            {
              jump = w1[0] / w2[0];
            }
          else
            {
              jump = w2[0] / w1[0];
            }

          if(jump > jump_max)
            {
              jump_max = jump;
            }
        }
    }
#ifdef OUTPUT_DG_MACHNUMS
  SphP[cell].Machnum = calculate_machnumber_rho_jump(jump_max);
#endif

  return calculate_machnumber_rho_jump(jump_max);
}

/*!
 * Given a temperature jump, calculate the Mach number
 */
static double calculate_machnumber_T_jump(double T_jump)
{
  double solution1;
  double solution2;

  //solve_quadratic_equation(5, 14 - 16 * T_postshock / T_preshock, -3, &solution1, &solution2);  //gamma = 5/3
  solve_quadratic_equation(2 * GAMMA * (GAMMA - 1), 4 * GAMMA - (GAMMA - 1) * (GAMMA - 1) - (GAMMA + 1) * (GAMMA + 1) * T_jump, -2 * (GAMMA - 1), &solution1, &solution2);      //general gamma

  return sqrt(solution1);
}

static double calculate_machnumber_rho_jump(double rho_jump)
{
  if(rho_jump > 3)
    {
      return 3;
    }
  else
    {
      return sqrt(2. / (1. / rho_jump * (GAMMA + 1) - (GAMMA - 1)));
    }
}

/*!
 * Given a Mach number, calculate the density jump
 */
static double calculate_T_jump(double M)
{
  return (2 * GAMMA * M * M - (GAMMA - 1)) * ((GAMMA - 1) * M * M + 2) / ((GAMMA + 1) * (GAMMA + 1) * M * M);
}

#endif
#endif
