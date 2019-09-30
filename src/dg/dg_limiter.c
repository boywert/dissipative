/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_limiter.c
 * \date        10/2014
 * \author		  Kevin Schaal
 * \brief		    Functions for limiting the higher order terms
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

//3d: positivity limiter done

#include "../allvars.h"
#include "../proto.h"
#include "dg_core_inline.h"

#ifdef DG

//general functions
static double tripple_min(double a, double b, double c);
static double tripple_max(double a, double b, double c);


#ifdef POSITIVITY_LIMITER
static void info_positivity_limiter(double *w_mean, double *w, double sol1, double sol2, double a_, double b, double c, int t);
static double calc_tau(double w_mean[5], double pv_mean[5], double w[5]);
inline static int fix_mean_value(double a_cell[NOF_BASE_FUNCTIONS][5], int cell, double eps_rho, double eps_p);
#endif

#ifdef CHARACTERISTIC_LIMITER
static void calc_flux_jacobian_matrices(double w[5]);
#endif

#ifdef ANGLE_BOUND
static void min_angles(CBV(a), int cell, int dim, double min_angle[5]);
#endif

#if defined(ANGLE_BOUND) && defined(CHARACTERISTIC_LIMITER)
static void calculate_characteristic_slopes(CBV(a));
static double (*s_char_slopes)[NUMDIMS][5];     // char slopes[cell][dim][char_var]
#endif


/*!
 * minimum of three numbers
 */
static double tripple_min(double a, double b, double c)
{
  return fmin(fmin(a, b), c);
}

/*!
 * maximum of three numbers
 */
static double tripple_max(double a, double b, double c)
{
  return fmax(fmax(a, b), c);
}


/*!
 * multivariable minmod function
 */
double minmod(double a, double b, double c)
{
  if(a > 0 && b > 0 && c > 0)
    {
      return tripple_min(a, b, c);
    }
  else if(a < 0 && b < 0 && c < 0)
    {
      return tripple_max(a, b, c);
    }
  else
    {
      return 0;
    }
}

/*!
 * total variation bounded version
 */
#ifdef MINMOD_B
double minmodB(double a, double b, double c, double h)
{
  //All.DG_M corresponds to M_tilde, not the M in Praveen's report!
  if(fabs(a) < All.DG_M * h / sqrt(3))
    {
      return a;
    }
  else
    {
      return minmod(a, b, c);
    }
}
#endif

#ifdef ANGLE_BOUND
/*!
 * Returns the angle difference to neighbouring cells
 * \param CBV(a) the state of all cells
 * \param cell the index of the central cell
 * \param dim indicates the dimension, x=0, y=1, z=2
 * \param min_angle output parameter, the minimum angle difference for each conserved variable is stored
 */
static void min_angles(CBV(a), int cell, int dim, double min_angle[5])
{
#ifdef CHARACTERISTIC_LIMITER
  //neighbour indices
  int index_ngb_neg;            //negative neighbour (left, back, bottom)
  int index_ngb_pos;            //positive neighbour (right, front, top)
#endif

  enum E_interface interface_pos;
  enum E_interface interface_neg;

  //slopes
  double a_wx_ngb_pos[5];
  double a_wx_ngb_neg[5];

  double *wx_ngb_pos = 0;
  double *wx_ngb_neg = 0;
  double *wx_center = 0;

#ifdef CHARACTERISTIC_LIMITER
  double temp[5];

  //mean values
  double mean_ngb_pos[5];
  double mean_ngb_neg[5];

  //transformation matrix
  double (*LL)[5];
#endif

  if(dim == 0)
    {
      interface_neg = LEFT;
      interface_pos = RIGHT;

#ifdef CHARACTERISTIC_LIMITER
      index_ngb_neg = Mesh.DP[cell].neighbors[LEFT];
      index_ngb_pos = Mesh.DP[cell].neighbors[RIGHT];
      LL = Lx;
#endif
    }
  else if(dim == 1)
    {
      interface_neg = BACK;
      interface_pos = FRONT;

#ifdef CHARACTERISTIC_LIMITER
      index_ngb_neg = Mesh.DP[cell].neighbors[BACK];
      index_ngb_pos = Mesh.DP[cell].neighbors[FRONT];
      LL = Ly;
#endif
    }
#ifndef TWODIMS
  else if(dim == 2)
    {
      interface_neg = BOTTOM;
      interface_pos = TOP;

#ifdef CHARACTERISTIC_LIMITER
      index_ngb_neg = Mesh.DP[cell].neighbors[BOTTOM];
      index_ngb_pos = Mesh.DP[cell].neighbors[TOP];
      LL = Lz;
#endif
    }
#endif
  else
    {
      interface_neg = 0;
      interface_pos = 0;
#ifdef CHARACTERISTIC_LIMITER
      index_ngb_neg = 0;
      index_ngb_pos = 0;
      LL = 0;
#endif

      terminate("invalid dim!\n");
    }

#ifdef CONSERVED_LIMITER
  //get slopes in dim-direction
  get_neighbour_weights(cell, interface_pos, dim + 1, a, a_wx_ngb_pos);
  get_neighbour_weights(cell, interface_neg, dim + 1, a, a_wx_ngb_neg);

  wx_center = a[cell][dim + 1];
  wx_ngb_pos = &a_wx_ngb_pos[0];
  wx_ngb_neg = &a_wx_ngb_neg[0];

#elif defined(CHARACTERISTIC_LIMITER)
  //transform to slopes of characteristic variables

  wx_center = s_char_slopes[cell][dim];

  if(index_ngb_pos < Ngb_MaxPart && Mesh.DP[cell].level == Mesh.DP[index_ngb_pos].level)        //neighbour is local cell
    {
      wx_ngb_pos = s_char_slopes[index_ngb_pos][dim];
    }
  else
    {
      get_neighbour_weights(cell, interface_pos, 0, a, mean_ngb_pos);
      get_neighbour_weights(cell, interface_pos, dim + 1, a, a_wx_ngb_pos);

      calc_flux_jacobian_matrices(mean_ngb_pos);
      copy_vector_vector(a_wx_ngb_pos, temp, 5);
      multiply_matrix_vector(LL, temp, a_wx_ngb_pos);

      wx_ngb_pos = &a_wx_ngb_pos[0];
    }

  if(index_ngb_neg < Ngb_MaxPart && Mesh.DP[cell].level == Mesh.DP[index_ngb_neg].level)        //neighbour is local cell
    {
      wx_ngb_neg = s_char_slopes[index_ngb_neg][dim];
    }
  else
    {
      get_neighbour_weights(cell, interface_neg, 0, a, mean_ngb_neg);
      get_neighbour_weights(cell, interface_neg, dim + 1, a, a_wx_ngb_neg);

      calc_flux_jacobian_matrices(mean_ngb_neg);
      copy_vector_vector(a_wx_ngb_neg, temp, 5);
      multiply_matrix_vector(LL, temp, a_wx_ngb_neg);

      wx_ngb_neg = &a_wx_ngb_neg[0];
    }
#endif


  //cell size
  double dx_center = amr_length[Mesh.DP[cell].level];

  //the neighbours have the same cell size, since
  //get neighbour weights projects in case of different cell sizes
  double dx_ngb_pos = amr_length[Mesh.DP[cell].level];
  double dx_ngb_neg = amr_length[Mesh.DP[cell].level];

  int k;
  double cos_phi;               //has to be in [-1,1]
  double phi;                   //angle, has to be in [0,pi/2.]

  for(k = 0; k < 5; k++)
    {
      //slope difference to positive neighbour (center cell m gets a minus sign)
      cos_phi = 3 * wx_ngb_pos[k] * (-wx_center[k]) + (-dx_center) * dx_ngb_pos / 4.;
      cos_phi /= sqrt(3 * wx_ngb_pos[k] * wx_ngb_pos[k] + dx_ngb_pos * dx_ngb_pos / 4.) * sqrt(3 * wx_center[k] * wx_center[k] + dx_center * dx_center / 4.);   //minus sign irrelevant

      //correct rounding errors, cos_phi has to be in [-1,1]
      cos_phi = fmin(fmax(-1, cos_phi), 1);

      phi = acos(cos_phi);
      min_angle[k] = phi;

      //slope difference to negative neighbour (negative neighbour gets minus sign)
      cos_phi = 3 * (-wx_ngb_neg[k]) * wx_center[k] + dx_center * (-dx_ngb_neg) / 4.;
      cos_phi /= sqrt(3 * wx_ngb_neg[k] * wx_ngb_neg[k] + dx_ngb_neg * dx_ngb_neg / 4.) * sqrt(3 * wx_center[k] * wx_center[k] + dx_center * dx_center / 4.);   //minus sign irrelevant

      //correct rounding errors, cos_phi has to be in [-1,1]
      cos_phi = fmin(fmax(-1, cos_phi), 1);

      phi = acos(cos_phi);
      min_angle[k] = fmin(min_angle[k], phi);
    }
}
#endif


/*!
 * apply the componentwise minmod limiter to a weights array
 */
void minmod_limiter(CBV(a))
{
  mpi_printf("DG: Minmod limiter\n");

#if (DEGREE_K==0)
  return;
#endif

  TIMER_START(CPU_DG_LIMITER);

  double (*a_temp)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("a_temp", NumGas * Nof_base_functions * 5 * sizeof(double));
  copy_array_to_array(a, a_temp);

  int i;

#ifdef MINMOD_B
  double dl;
#endif

#ifdef MACHNUM_B
  double jump = calculate_rho_jump(All.DG_lim_mach_min);
#endif

#ifdef NEIGHBOUR_BOUND
  double w_max_interface[4];
  double w_min_interface[4];
  double w_max_neighbours[4];
  double w_min_neighbours[4];
  int nb_limiting;
#endif

#ifdef ANGLE_BOUND

#if !(defined(CONSERVED_LIMITER) || defined(CHARACTERISTIC_LIMITER))
#error Define a limiter!
#endif

  double min_angles_x[5];
  double min_angles_y[5];
  DG3D(double min_angles_z[5];)
#ifdef CHARACTERISTIC_LIMITER
    s_char_slopes = (double (*)[NUMDIMS][5]) mymalloc("s_char_slopes", NumGas * NUMDIMS * 5 * sizeof(double));  // char slopes[cell][dim][char_var]
  calculate_characteristic_slopes(a);
#endif

#endif

#if defined(CHARACTERISTIC_LIMITER) || defined(CONSERVED_LIMITER)
  double eps = 1.0e-8;
  int k, t;

  double w0_left[5];
  double w0_right[5];
  double w0_back[5];
  double w0_front[5];
#ifndef TWODIMS
  double w0_bottom[5];
  double w0_top[5];
#endif
#endif

#ifdef CHARACTERISTIC_LIMITER
  double dw_left[5], dw_right[5], dw_front[5], dw_back[5];
  double *c2, *c3;
  double c2_new[5], c3_new[5], c2_left[5], c2_right[5], c3_back[5], c3_front[5];
  int triggered;

#ifndef TWODIMS
  double dw_top[5], dw_bottom[5];
  double *c4;
  double c4_new[5], c4_bottom[5], c4_top[5];
#endif

#ifndef ANGLE_BOUND
  double c2_[5], c3_[5];
  DG3D(double c4_[5];)
#endif
#endif
#ifdef CONSERVED_LIMITER
  double w1new, w2new;
  DG3D(double w3new;)
#endif
    //loop over cells
    for(i = 0; i < NumGas; i++)
    {

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

#ifdef ANGLE_BOUND
      min_angles(a, i, 0, min_angles_x);
      min_angles(a, i, 1, min_angles_y);
      DG3D(min_angles(a, i, 2, min_angles_z);)
#ifdef OUTPUT_MIN_ANGLES
        SphP[i].min_angles_x[0] = min_angles_x[0];
      SphP[i].min_angles_x[1] = min_angles_x[1];
      SphP[i].min_angles_x[2] = min_angles_x[2];
      SphP[i].min_angles_x[3] = min_angles_x[3];
      SphP[i].min_angles_x[4] = min_angles_x[4];

      SphP[i].min_angles_y[0] = min_angles_y[0];
      SphP[i].min_angles_y[1] = min_angles_y[1];
      SphP[i].min_angles_y[2] = min_angles_y[2];
      SphP[i].min_angles_y[3] = min_angles_y[3];
      SphP[i].min_angles_y[4] = min_angles_y[4];

#ifndef TWODIMS
      SphP[i].min_angles_z[0] = min_angles_z[0];
      SphP[i].min_angles_z[1] = min_angles_z[1];
      SphP[i].min_angles_z[2] = min_angles_z[2];
      SphP[i].min_angles_z[3] = min_angles_z[3];
      SphP[i].min_angles_z[4] = min_angles_z[4];
#endif
#endif
#endif

#ifdef MACHNUM_JUMP_DETECTION
      calculate_mach_number(a, i);
#endif

#ifdef DISCONTINUITY_DETECTION
      if(!detect_discontinuity(a, i))   //no discontinuity detected, apply only the positivity limiter
        {
#ifdef POSITIVITY_LIMITER
          positivity_limiter(a_temp[i], i);
#endif
          continue;
        }
#endif //end DISCONTINUITY_DETECTION

#ifdef MINMOD_B
      dl = amr_length[Mesh.DP[i].level];
#endif
//done

#if defined(CHARACTERISTIC_LIMITER) || defined(CONSERVED_LIMITER)
      //get the mean states of the neighbours
      get_neighbour_weights(i, LEFT, 0, a, w0_left);
      get_neighbour_weights(i, RIGHT, 0, a, w0_right);
      get_neighbour_weights(i, BACK, 0, a, w0_back);
      get_neighbour_weights(i, FRONT, 0, a, w0_front);
#ifndef TWODIMS
      get_neighbour_weights(i, BOTTOM, 0, a, w0_bottom);
      get_neighbour_weights(i, TOP, 0, a, w0_top);
#endif
#endif

#ifdef NEIGHBOUR_BOUND
      min_max_conserved_values_on_interface(a[i], i, w_max_interface, w_min_interface);
      min_max_mean_conserved_values_of_neighbours(w0_left, w0_right, w0_back, w0_front, w_max_neighbours, w_min_neighbours);
#endif

#ifdef CHARACTERISTIC_LIMITER

      DG_PRINTF("Characteristic Limiter\n");

#ifdef CONSERVED_LIMITER
#error Limit either conserved OR characteristic variables!
#endif

      /*
       * w_mean: a[i][0]
       * w_slope_x: a[i][1]
       * w_slope_y: a[i][2]
       * w_slope_z: a[i][3]
       *
       */

      sub_vector_vector(w0_right, a[i][0], dw_right);
      sub_vector_vector(a[i][0], w0_left, dw_left);

      sub_vector_vector(w0_front, a[i][0], dw_front);
      sub_vector_vector(a[i][0], w0_back, dw_back);
#ifndef TWODIMS
      sub_vector_vector(w0_top, a[i][0], dw_top);
      sub_vector_vector(a[i][0], w0_bottom, dw_bottom);
#endif


      calc_flux_jacobian_matrices(a[i][0]);

#ifdef ANGLE_BOUND

      c2 = s_char_slopes[i][0];
      c3 = s_char_slopes[i][1];
      DG3D(c4 = s_char_slopes[i][2];)
#else
      multiply_matrix_vector(Lx, a[i][1], c2_);
      multiply_matrix_vector(Ly, a[i][2], c3_);
      DG3D(multiply_matrix_vector(Lz, a[i][3], c4_);) c2 = &c2_[0];
      c3 = &c3_[0];
      DG3D(c4 = &c4_[0];)
#endif
        multiply_matrix_vector(Lx, dw_left, c2_left);
      multiply_matrix_vector(Lx, dw_right, c2_right);

      multiply_matrix_vector(Ly, dw_back, c3_back);
      multiply_matrix_vector(Ly, dw_front, c3_front);

#ifndef TWODIMS
      multiply_matrix_vector(Lz, dw_bottom, c4_bottom);
      multiply_matrix_vector(Lz, dw_top, c4_top);
#endif

      //loop over components
      triggered = 0;

#ifdef NEIGHBOUR_BOUND
      nb_limiting = 1;

      if(w_max_interface[0] <= w_max_neighbours[0] && w_min_interface[0] >= w_min_neighbours[0])        //no limiting!
        {
          nb_limiting = 0;
        }
#endif

      for(k = 0; k < 5; k++)
        {
#ifdef MINMOD_B
          c2_new[k] = minmodB(c2[k], All.DG_beta / sqrt(3.) * c2_left[k], All.DG_beta / sqrt(3.) * c2_right[k], dl);
          c3_new[k] = minmodB(c3[k], All.DG_beta / sqrt(3.) * c3_back[k], All.DG_beta / sqrt(3.) * c3_front[k], dl);
          DG3D(c4_new[k] = minmodB(c4[k], All.DG_beta / sqrt(3.) * c4_bottom[k], All.DG_beta / sqrt(3.) * c4_top[k], dl);)
#elif defined(ANGLE_BOUND)
          if(min_angles_x[k] > All.DG_min_angle)        //no limiting
            {
              c2_new[k] = c2[k];
            }
          else                  //limiting
            {
              c2_new[k] = minmod(c2[k], All.DG_beta / sqrt(3.) * c2_left[k], All.DG_beta / sqrt(3.) * c2_right[k]);
            }

          if(min_angles_y[k] > All.DG_min_angle)        //no limiting
            {
              c3_new[k] = c3[k];
            }
          else                  //limiting
            {
              c3_new[k] = minmod(c3[k], All.DG_beta / sqrt(3.) * c3_back[k], All.DG_beta / sqrt(3.) * c3_front[k]);
            }
#ifndef TWODIMS
          if(min_angles_z[k] > All.DG_min_angle)        //no limiting
            {
              c4_new[k] = c4[k];
            }
          else                  //limiting
            {
              c4_new[k] = minmod(c4[k], All.DG_beta / sqrt(3.) * c4_bottom[k], All.DG_beta / sqrt(3.) * c4_top[k]);
            }
#endif
#elif defined(NEIGHBOUR_BOUND)
          if(!nb_limiting)
            {
              c2_new[k] = c2[k];
              c3_new[k] = c3[k];
            }
          else
            {
              c2_new[k] = minmod(c2[k], All.DG_beta / sqrt(3.) * c2_left[k], All.DG_beta / sqrt(3.) * c2_right[k]);
              c3_new[k] = minmod(c3[k], All.DG_beta / sqrt(3.) * c3_back[k], All.DG_beta / sqrt(3.) * c3_front[k]);
            }
#elif defined(MACHNUM_B)
          if((a[i][0][0] + fabs(a[i][1][0]) * sqrt(3.)) / (a[i][0][0] - fabs(a[i][1][0]) * sqrt(3.)) > jump)
            {
              c2_new[k] = minmod(c2[k], All.DG_beta / sqrt(3.) * c2_left[k], All.DG_beta / sqrt(3.) * c2_right[k]);
            }
          else
            {
              c2_new[k] = c2[k];
            }

          if((a[i][0][0] + fabs(a[i][2][0]) * sqrt(3.)) / (a[i][0][0] - fabs(a[i][2][0]) * sqrt(3.)) > jump)
            {
              c3_new[k] = minmod(c3[k], All.DG_beta / sqrt(3.) * c3_back[k], All.DG_beta / sqrt(3.) * c3_front[k]);
            }
          else
            {
              c3_new[k] = c3[k];
            }
#else
          c2_new[k] = minmod(c2[k], All.DG_beta / sqrt(3.) * c2_left[k], All.DG_beta / sqrt(3.) * c2_right[k]);
          c3_new[k] = minmod(c3[k], All.DG_beta / sqrt(3.) * c3_back[k], All.DG_beta / sqrt(3.) * c3_front[k]);
          DG3D(c4_new[k] = minmod(c4[k], All.DG_beta / sqrt(3.) * c4_bottom[k], All.DG_beta / sqrt(3.) * c4_top[k]);)
#endif
            if(fabs(c2_new[k] - c2[k]) > eps || fabs(c3_new[k] - c3[k]) > eps DG3D(||fabs(c4_new[k] - c4[k]) > eps))
            {
              triggered = 1;
            }
        }

      if(triggered)
        {
          //discard higher order (non-linear) terms
          for(t = NUMDIMS + 1; t < Nof_base_functions; t++)
            {
              a_temp[i][t][0] = 0;
              a_temp[i][t][1] = 0;
              a_temp[i][t][2] = 0;
              a_temp[i][t][3] = 0;
              a_temp[i][t][4] = 0;
            }

          //calculate new slopes
          multiply_matrix_vector(Rx, c2_new, a_temp[i][1]);
          multiply_matrix_vector(Ry, c3_new, a_temp[i][2]);
        DG3D(multiply_matrix_vector(Rz, c4_new, a_temp[i][3]);)}
#endif


#ifdef CONSERVED_LIMITER

      DG_PRINTF("Conserved Limiter\n");

      //loop over conserved variables
      for(k = 0; k < 5; k++)
        {
#ifdef NEIGHBOUR_BOUND
          if(w_max_interface[k] <= w_max_neighbours[k] && w_min_interface[k] >= w_min_neighbours[k])    //no limiting!
            {
              continue;
            }
#elif defined(ANGLE_BOUND)
          if(min_angles_x[k] > All.DG_min_angle)        //no limiting
            {
              w1new = a[i][1][k];
            }
          else                  //limiting
            {
              w1new = minmod(a[i][1][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_left[k]), All.DG_beta / sqrt(3.) * (w0_right[k] - a[i][0][k]));
            }

          if(min_angles_y[k] > All.DG_min_angle)        //no limiting
            {
              w2new = a[i][2][k];
            }
          else                  //limiting
            {
              w2new = minmod(a[i][2][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_back[k]), All.DG_beta / sqrt(3.) * (w0_front[k] - a[i][0][k]));
            }
#ifndef TWODIMS
          if(min_angles_z[k] > All.DG_min_angle)        //no limiting
            {
              w3new = a[i][3][k];
            }
          else                  //limiting
            {
              w3new = minmod(a[i][3][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_bottom[k]), All.DG_beta / sqrt(3.) * (w0_top[k] - a[i][0][k]));
            }
#endif
#elif defined(MINMOD_B)
          w1new = minmodB(a[i][1][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_left[k]), All.DG_beta / sqrt(3.) * (w0_right[k] - a[i][0][k]), dl);
          w2new = minmodB(a[i][2][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_back[k]), All.DG_beta / sqrt(3.) * (w0_front[k] - a[i][0][k]), dl);
          DG3D(w3new = minmodB(a[i][3][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_bottom[k]), All.DG_beta / sqrt(3.) * (w0_top[k] - a[i][0][k]), dl);)
#else
          w1new = minmod(a[i][1][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_left[k]), All.DG_beta / sqrt(3.) * (w0_right[k] - a[i][0][k]));
          w2new = minmod(a[i][2][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_back[k]), All.DG_beta / sqrt(3.) * (w0_front[k] - a[i][0][k]));
          DG3D(w3new = minmod(a[i][3][k], All.DG_beta / sqrt(3.) * (a[i][0][k] - w0_bottom[k]), All.DG_beta / sqrt(3.) * (w0_top[k] - a[i][0][k]));)
#endif
            if(fabs(w1new - a[i][1][k]) > eps || fabs(w2new - a[i][2][k]) > eps DG3D(||fabs(w3new - a[i][3][k]) > eps)) //limiter triggered
            {
              a_temp[i][1][k] = w1new;
              a_temp[i][2][k] = w2new;
              DG3D(a_temp[i][3][k] = w3new;)
                //discard higher order terms
                for(t = NUMDIMS + 1; t < Nof_base_functions; t++)
                {
                  a_temp[i][t][k] = 0;
                }
            }
        }
#endif //end normal limiter

#ifdef ANGULAR_MOMENTUM_LIMITER
      angular_momentum_limiter(a_temp, i, spin(a[i], i));
#endif
    }                           //end loop over cells

  copy_array_to_array(a_temp, a);

#ifdef POSITIVITY_LIMITER
  for(i = 0; i < NumGas; i++)
    {
      positivity_limiter(a[i], i);
    }
#endif

  //free memory

#if defined(ANGLE_BOUND) && defined(CHARACTERISTIC_LIMITER)
  myfree(s_char_slopes);
#endif

  myfree(a_temp);

  TIMER_STOP(CPU_DG_LIMITER);
}

#if defined(ANGLE_BOUND) && defined(CHARACTERISTIC_LIMITER)
static void calculate_characteristic_slopes(CBV(a))
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      calc_flux_jacobian_matrices(a[i][0]);
      multiply_matrix_vector(Lx, a[i][1], s_char_slopes[i][0]);
      multiply_matrix_vector(Ly, a[i][2], s_char_slopes[i][1]);
    DG3D(multiply_matrix_vector(Lz, a[i][3], s_char_slopes[i][2]);)}
}
#endif

#ifdef POSITIVITY_LIMITER

static double generated_mass = 0;
static double generated_internal_energy = 0;

inline static int fix_mean_value(double a_cell[NOF_BASE_FUNCTIONS][5], int cell, double eps_rho, double eps_p)
{
  double w_mean[5], pv_mean[5];
  double dp = 0;
  double drho = 0;

  w_mean[0] = a_cell[0][0];
  w_mean[1] = a_cell[0][1];
  w_mean[2] = a_cell[0][2];
  w_mean[3] = a_cell[0][3];
  w_mean[4] = a_cell[0][4];

  w_to_primvars(w_mean, pv_mean);

  int fixed = 0;

  if(pv_mean[PRES] < eps_p)
    {
      dp = eps_p - pv_mean[PRES];
      pv_mean[PRES] = eps_p;
      fixed = 1;
    }

  if(pv_mean[RHO] < eps_rho)
    {
      drho = eps_rho - pv_mean[RHO];
      pv_mean[RHO] = eps_rho;
      fixed = 1;
    }

  if(fixed)
    {
      //set new mean values
      primvars_to_w(pv_mean, a_cell[0]);

      int l;

      //discard non-constant terms
      for(l = 1; l < Nof_base_functions; l++)
        {
          a_cell[l][0] = 0;
          a_cell[l][1] = 0;
          a_cell[l][2] = 0;
          a_cell[l][3] = 0;
          a_cell[l][4] = 0;
        }

      generated_mass += drho * SphP[cell].Volume;
      generated_internal_energy += dp / (GAMMA - 1) * SphP[cell].Volume;
    }

  return fixed;
}

static double calc_tau(double w_mean[5], double pv_mean[5], double w[5])
{
  double pv[5];

#ifdef DG_DEBUG
  double pv_tau[5];
  double w_tau[5];
#endif

  w_to_primvars(w, pv);

  double drho, dm1, dm2, dm3, dE, a_, b, c, tau, sol1, sol2;
  int t;


  if(pv[RHO] < Epsilon_rho)
    {
      printf("rho_eps=%.16e, rho=%.16e, rho_mean=%.16e\n", Epsilon_rho, pv[RHO], pv_mean[RHO]);
      terminate("Error in positivity limiter: density not limited correctly!");
    }


  if(pv[PRES] >= Epsilon_p || same_w(w, w_mean))        //no limiting necessary
    {
      return 1;
    }
  else                          //solve quadratic equation
    {
      drho = w[0] - w_mean[0];
      dm1 = w[1] - w_mean[1];
      dm2 = w[2] - w_mean[2];
      dm3 = w[3] - w_mean[3];
      dE = w[4] - w_mean[4];

      a_ = dE * drho - 0.5 * (dm1 * dm1 + dm2 * dm2 + dm3 * dm3);
      b = (w_mean[4] - Epsilon_p / (GAMMA - 1.)) * drho + w_mean[0] * dE - (w_mean[1] * dm1 + w_mean[2] * dm2 + w_mean[3] * dm3);
      c = (w_mean[4] - Epsilon_p / (GAMMA - 1.)) * w_mean[0] - 0.5 * (w_mean[1] * w_mean[1] + w_mean[2] * w_mean[2] + w_mean[3] * w_mean[3]);

      if(a_ == 0)               //can happen due to interaction with TVD limiter
        {
          tau = -c / b;
          sol1 = 0;
          sol2 = 0;
          t = 0;

          assert(tau >= 0 && tau <= 1);

          return tau;
        }
      else
        {
          solve_quadratic_equation(a_, b, c, &sol1, &sol2);

          t = 0;

          //improve the solution with Newton iterations
          {
            int converged = 0;
            double change;
            double sol_new;

            for(t = 0; t < 100; t++)
              {
                sol_new = sol2 - (a_ * sol2 * sol2 + b * sol2 + c) / (2 * sol2 * a_ + b);

                change = 2 * fabs((sol_new - sol2) / (sol_new + sol2));

                sol2 = sol_new;

                if(change < 1e-12 && t >= 2)
                  {
                    converged = 1;
                    t++;
                    break;
                  }
              }

            if(!converged)
              {
                info_positivity_limiter(w_mean, w, sol1, sol2, a_, b, c, t);
                terminate("Error in positivity limiter: Newton iterations did not converge!\n");
              }
          }

          tau = sol2;

          //consistency check
          if(!(tau > -1e-10 && tau < 1 + 1e-10))
            {
              info_positivity_limiter(w_mean, w, sol1, sol2, a_, b, c, t);
              terminate("Error in positivity limiter: Tau is not in the expected range!\n");
            }

          tau = fmax(0, fmin(1.0, tau));        //enforce tau to be in [0,1]

          //safety factor
#ifdef OLD_VERSION
          tau = tau * (1. - 1e-6);
#else
          tau = tau * (1. - 1e-10);
#endif

#ifdef DG_DEBUG
          //consistency check
          w_tau[0] = w_mean[0] + tau * (w[0] - w_mean[0]);
          w_tau[1] = w_mean[1] + tau * (w[1] - w_mean[1]);
          w_tau[2] = w_mean[2] + tau * (w[2] - w_mean[2]);
          w_tau[3] = w_mean[3] + tau * (w[3] - w_mean[3]);
          w_tau[4] = w_mean[4] + tau * (w[4] - w_mean[4]);

          w_to_primvars(w_tau, pv_tau);

          if(pv_tau[PRES] < Epsilon_p)
            {
              info_positivity_limiter(w_mean, w, sol1, sol2, a_, b, c, t);
              printf("p_eps=%.16e, p=%.16e, p_mean=%.16e\n", Epsilon_p, pv_tau[PRES], pv_mean[PRES]);
              terminate("Error in positivity limiter, pressure not limited correctly!\n");
            }
#endif

          return tau;
        }
    }
}

void positivity_limiter(double a_cell[NOF_BASE_FUNCTIONS][5], int cell)
{
  if(P[cell].Mass == 0 && P[cell].ID == 0)
    return;                     /* skip dissolved cells */

#ifdef FIX_MEAN_VALUES
  if(fix_mean_value(a_cell, cell, Epsilon_rho, Epsilon_p))      //mean values are below minimum values
    {
      return;
    }
#endif

#ifdef DG_REFINEMENT
  int i;
#endif

  double rho, rho_mean, rho_min, lambda_rho;
  double tau, tau_min;

  double w[5];
  double w_mean[5], pv_mean[5];

  int q, l;

  rho_min = DBL_MAX;
  tau_min = DBL_MAX;

  //limit the density
  //

  rho_mean = a_cell[0][0];

  if(rho_mean < Epsilon_rho)
    {
      lambda_rho = 0;
    }
  else                          //compute minimum density of union quadrature points
    {
      for(q = 0; q < Nof_union_quad_points; q++)
        {
          rho = 0;

          for(l = 0; l < Nof_base_functions; l++)
            {
              rho += a_cell[l][0] * GET_Union_base_values(q, l);
            }

          if(rho < rho_min)
            {
              rho_min = rho;
            }
        }

#ifdef DG_REFINEMENT

      for(i = 0; i < Nof_interfaces; i++)       //speedup: only ngb interfaces with smaller cells
        {
          for(q = 0; q < Nof_outer_fine_quad_points; q++)
            {
              rho = 0;

              for(l = 0; l < Nof_base_functions; l++)
                {
                  rho += a_cell[l][0] * GET_Outer_fine_base_values(i, q, l);
                }

              if(rho < rho_min)
                {
                  rho_min = rho;
                }
            }
        }
#endif

      lambda_rho = fmin(fabs((rho_mean - Epsilon_rho) / (rho_mean - rho_min)), 1);
    }

  //limit
  if(lambda_rho < 1)
    {
      //safety factor
#ifdef OLD_VERSION
      lambda_rho = lambda_rho * (1 - 1e-6);
#else
      lambda_rho = lambda_rho * (1 - 1e-10);
#endif

      for(l = 1; l < Nof_base_functions; l++)
        {
          a_cell[l][0] = lambda_rho * a_cell[l][0];
        }
    }

  //limit the pressure
  //

  w_mean[0] = a_cell[0][0];
  w_mean[1] = a_cell[0][1];
  w_mean[2] = a_cell[0][2];
  w_mean[3] = a_cell[0][3];
  w_mean[4] = a_cell[0][4];

  w_to_primvars(w_mean, pv_mean);

  if(pv_mean[PRES] < Epsilon_p || pv_mean[RHO] < Epsilon_rho)
    {
      printf("cell (%f, %f, %f), mean rho:%g, mean p:%g\n", P[cell].Pos[0], P[cell].Pos[1], P[cell].Pos[2], pv_mean[RHO], pv_mean[PRES]);
      terminate("Error in positivity limiter: negative mean values!\n");
    }

  for(q = 0; q < Nof_union_quad_points; q++)
    {
      calc_state_at_union_quad_point(a_cell, q, w);

      tau = calc_tau(w_mean, pv_mean, w);

      if(tau < tau_min)
        {
          tau_min = tau;
        }

      if(tau_min == 0)          //minimum found
        {
          break;
        }
    }                           //end loop over quadrature points, tau_min set

#ifdef DG_REFINEMENT
  int stop = (tau_min == 0);

  for(i = 0; (i < Nof_interfaces) && !stop; i++)        //speedup: only ngb interfaces with smaller cells
    {
      for(q = 0; q < Nof_outer_fine_quad_points; q++)
        {
          calc_state_at_outer_fine_quad_point(a_cell, i, q, w);

          tau = calc_tau(w_mean, pv_mean, w);

          if(tau < tau_min)
            {
              tau_min = tau;
            }

          if(tau_min == 0)      //minimum found
            {
              stop = 1;
              break;
            }
        }
    }
#endif

  if(tau_min < 1)
    {
      for(l = 1; l < Nof_base_functions; l++)
        {
          a_cell[l][0] = tau_min * a_cell[l][0];
          a_cell[l][1] = tau_min * a_cell[l][1];
          a_cell[l][2] = tau_min * a_cell[l][2];
          a_cell[l][3] = tau_min * a_cell[l][3];
          a_cell[l][4] = tau_min * a_cell[l][4];
        }
    }
}

static void info_positivity_limiter(double *w_mean, double *w, double sol1, double sol2, double a_, double b, double c, int t)
{
  printf("Task: %d\n", ThisTask);

  double w_tau_plus[5];
  double w_tau_minus[5];
  double dw[5];
  double pv_plus[5];
  double pv_minus[5];
  double pv[5];
  double pv_mean[5];

  dw[0] = w[0] - w_mean[0];
  dw[1] = w[1] - w_mean[1];
  dw[2] = w[2] - w_mean[2];
  dw[3] = w[3] - w_mean[3];
  dw[4] = w[4] - w_mean[4];

  w_tau_plus[0] = w_mean[0] + sol1 * (w[0] - w_mean[0]);
  w_tau_plus[1] = w_mean[1] + sol1 * (w[1] - w_mean[1]);
  w_tau_plus[2] = w_mean[2] + sol1 * (w[2] - w_mean[2]);
  w_tau_plus[3] = w_mean[3] + sol1 * (w[3] - w_mean[3]);
  w_tau_plus[4] = w_mean[4] + sol1 * (w[4] - w_mean[4]);

  w_tau_minus[0] = w_mean[0] + sol2 * (w[0] - w_mean[0]);
  w_tau_minus[1] = w_mean[1] + sol2 * (w[1] - w_mean[1]);
  w_tau_minus[2] = w_mean[2] + sol2 * (w[2] - w_mean[2]);
  w_tau_minus[3] = w_mean[3] + sol2 * (w[3] - w_mean[3]);
  w_tau_minus[4] = w_mean[4] + sol2 * (w[4] - w_mean[4]);

  w_to_primvars(w_tau_plus, pv_plus);
  w_to_primvars(w_tau_minus, pv_minus);
  w_to_primvars(w, pv);
  w_to_primvars(w_mean, pv_mean);

  printf("w: %e, %e, %e, %e, %e\n", w[0], w[1], w[2], w[3], w[4]);
  printf("p:%.16e, rho:%.16e\n", pv[PRES], pv[RHO]);
  printf("p_mean:%.16e, rho_mean:%.16e\n", pv_mean[PRES], pv_mean[RHO]);
  printf("w_mean: %e, %e, %e, %e, %e\n", w_mean[0], w_mean[1], w_mean[2], w_mean[3], w_mean[4]);
  printf("dw: %e, %e, %e, %e, %e\n", dw[0], dw[1], dw[2], dw[3], dw[4]);

  printf("sol1: ax^2+bx+c=%.16e\n", a_ * sol1 * sol1 + b * sol1 + c);
  printf("sol2: ax^2+bx+c=%.16e\n", a_ * sol2 * sol2 + b * sol2 + c);
  printf("newton iterations: %d\n", t);
  printf("pressure sol1: %e, pressure sol2: %.16e\n", pv_plus[PRES], pv_minus[PRES]);
  printf("density sol1: %e, density sol2: %.16e\n", pv_plus[RHO], pv_minus[RHO]);
  printf("a:%e,b:%e,c:%e,sol1:%.16e,sol2:%.16e\n", a_, b, c, sol1, sol2);
}

#endif

#ifdef CHARACTERISTIC_LIMITER
static void calc_flux_jacobian_matrices(double state_w[5])
{
  double pv[5];
  w_to_primvars(state_w, pv);

  double g1 = GAMMA_MINUS1;
  double c = sqrt(GAMMA * pv[PRES] / pv[RHO]);
  double c2 = GAMMA * pv[PRES] / pv[RHO];
  double beta = 0.5 / c2;
  double k = 0.5 * (pv[VX] * pv[VX] + pv[VY] * pv[VY] + pv[VZ] * pv[VZ]);
  double phi = g1 * k;
  double h = c2 / g1 + k;

  double u = pv[VX];
  double v = pv[VY];
  double w = pv[VZ];

  //Lx
  //

  Lx[0][0] = 1 - phi / (c2);
  Lx[0][1] = g1 * u / (c2);
  Lx[0][2] = g1 * v / (c2);
  Lx[0][3] = 2 * beta * g1 * w;
  Lx[0][4] = -g1 / (c2);
  Lx[1][0] = v;
  Lx[1][1] = 0;
  Lx[1][2] = -1;
  Lx[1][3] = 0;
  Lx[1][4] = 0;
  Lx[2][0] = beta * (phi - c * u);
  Lx[2][1] = beta * (c - g1 * u);
  Lx[2][2] = -beta * g1 * v;
  Lx[2][3] = -beta * g1 * w;
  Lx[2][4] = beta * g1;
  Lx[3][0] = -w;
  Lx[3][1] = 0;
  Lx[3][2] = 0;
  Lx[3][3] = 1;
  Lx[3][4] = 0;
  Lx[4][0] = beta * (phi + c * u);
  Lx[4][1] = -beta * (c + g1 * u);
  Lx[4][2] = -beta * g1 * v;
  Lx[4][3] = -beta * g1 * w;
  Lx[4][4] = beta * g1;

  //Ly
  //

  Ly[0][0] = 1 - phi / (c2);
  Ly[0][1] = g1 * u / (c2);
  Ly[0][2] = g1 * v / (c2);
  Ly[0][3] = 2 * beta * g1 * w;
  Ly[0][4] = -g1 / (c2);
  Ly[1][0] = -u;
  Ly[1][1] = 1;
  Ly[1][2] = 0;
  Ly[1][3] = 0;
  Ly[1][4] = 0;
  Ly[2][0] = beta * (phi - c * v);
  Ly[2][1] = -beta * g1 * u;
  Ly[2][2] = beta * (c - g1 * v);
  Ly[2][3] = -beta * g1 * w;
  Ly[2][4] = beta * g1;
  Ly[3][0] = w;
  Ly[3][1] = 0;
  Ly[3][2] = 0;
  Ly[3][3] = -1;
  Ly[3][4] = 0;
  Ly[4][0] = beta * (phi + c * v);
  Ly[4][1] = -beta * g1 * u;
  Ly[4][2] = -beta * (c + g1 * v);
  Ly[4][3] = -beta * g1 * w;
  Ly[4][4] = beta * g1;

#ifndef TWODIMS
  //Lz
  //

  Lz[0][0] = 1 - 2 * beta * phi;
  Lz[0][1] = 2 * beta * g1 * u;
  Lz[0][2] = 2 * beta * g1 * v;
  Lz[0][3] = 2 * beta * g1 * w;
  Lz[0][4] = -2 * beta * g1;
  Lz[1][0] = u;
  Lz[1][1] = -1;
  Lz[1][2] = 0;
  Lz[1][3] = 0;
  Lz[1][4] = 0;
  Lz[2][0] = beta * (phi - c * w);
  Lz[2][1] = -beta * g1 * u;
  Lz[2][2] = -beta * g1 * v;
  Lz[2][3] = -beta * (g1 * w - c);
  Lz[2][4] = beta * g1;
  Lz[3][0] = -v;
  Lz[3][1] = 0;
  Lz[3][2] = 1;
  Lz[3][3] = 0;
  Lz[3][4] = 0;
  Lz[4][0] = beta * (phi + c * w);
  Lz[4][1] = -beta * g1 * u;
  Lz[4][2] = -beta * (g1 * v);
  Lz[4][3] = -beta * (g1 * w + c);
  Lz[4][4] = beta * g1;
#endif

  //Rx
  //

  Rx[0][0] = 1;
  Rx[0][1] = 0;
  Rx[0][2] = 1;
  Rx[0][3] = 0;
  Rx[0][4] = 1;
  Rx[1][0] = u;
  Rx[1][1] = 0;
  Rx[1][2] = u + c;
  Rx[1][3] = 0;
  Rx[1][4] = u - c;
  Rx[2][0] = v;
  Rx[2][1] = -1;
  Rx[2][2] = v;
  Rx[2][3] = 0;
  Rx[2][4] = v;
  Rx[3][0] = w;
  Rx[3][1] = 0;
  Rx[3][2] = w;
  Rx[3][3] = 1;
  Rx[3][4] = w;
  Rx[4][0] = k;
  Rx[4][1] = -v;
  Rx[4][2] = h + c * u;
  Rx[4][3] = w;
  Rx[4][4] = h - c * u;

  //Ry
  //

  Ry[0][0] = 1;
  Ry[0][1] = 0;
  Ry[0][2] = 1;
  Ry[0][3] = 0;
  Ry[0][4] = 1;
  Ry[1][0] = u;
  Ry[1][1] = 1;
  Ry[1][2] = u;
  Ry[1][3] = 0;
  Ry[1][4] = u;
  Ry[2][0] = v;
  Ry[2][1] = 0;
  Ry[2][2] = v + c;
  Ry[2][3] = 0;
  Ry[2][4] = v - c;
  Ry[3][0] = w;
  Ry[3][1] = 0;
  Ry[3][2] = w;
  Ry[3][3] = -1;
  Ry[3][4] = w;
  Ry[4][0] = k;
  Ry[4][1] = u;
  Ry[4][2] = h + c * v;
  Ry[4][3] = -w;
  Ry[4][4] = h - c * v;

#ifndef TWODIMS
  //Rz
  //

  Rz[0][0] = 1;
  Rz[0][1] = 0;
  Rz[0][2] = 1;
  Rz[0][3] = 0;
  Rz[0][4] = 1;
  Rz[1][0] = u;
  Rz[1][1] = -1;
  Rz[1][2] = u;
  Rz[1][3] = 0;
  Rz[1][4] = u;
  Rz[2][0] = v;
  Rz[2][1] = 0;
  Rz[2][2] = v;
  Rz[2][3] = 1;
  Rz[2][4] = v;
  Rz[3][0] = w;
  Rz[3][1] = 0;
  Rz[3][2] = w + c;
  Rz[3][3] = 0;
  Rz[3][4] = w - c;
  Rz[4][0] = k;
  Rz[4][1] = -u;
  Rz[4][2] = h + c * w;
  Rz[4][3] = v;
  Rz[4][4] = h - c * w;
#endif

#ifdef DG_DEBUG
  double idx[5][5];
  double idy[5][5];
  DG3D(double idz[5][5];) multiply_matrix_matrix(Lx, Rx, idx);
  multiply_matrix_matrix(Ly, Ry, idy);
  DG3D(multiply_matrix_matrix(Lz, Rz, idz);) if(!is_id(idx))
    {
      printf("u:%f, v:%f, c:%f, beta:%f, phi:%f,h:%f,k:%f,g1:%f\n", u, v, c, beta, phi, h, k, g1);
      printf("Lx*Rx!=Id\n");

      printf("Lx:\n");
      printf("%f\t\%f\t%f\t%f\t%f\n", Lx[0][0], Lx[0][1], Lx[0][2], Lx[0][3], Lx[0][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Lx[1][0], Lx[1][1], Lx[1][2], Lx[1][3], Lx[1][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Lx[2][0], Lx[2][1], Lx[2][2], Lx[2][3], Lx[2][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Lx[3][0], Lx[3][1], Lx[3][2], Lx[3][3], Lx[3][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Lx[4][0], Lx[4][1], Lx[4][2], Lx[4][3], Lx[4][4]);

      printf("Rx:\n");
      printf("%f\t\%f\t%f\t%f\t%f\n", Rx[0][0], Rx[0][1], Rx[0][2], Rx[0][3], Rx[0][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Rx[1][0], Rx[1][1], Rx[1][2], Rx[1][3], Rx[1][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Rx[2][0], Rx[2][1], Rx[2][2], Rx[2][3], Rx[2][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Rx[3][0], Rx[3][1], Rx[3][2], Rx[3][3], Rx[3][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", Rx[4][0], Rx[4][1], Rx[4][2], Rx[4][3], Rx[4][4]);


      printf("id:\n");
      printf("%f\t\%f\t%f\t%f\t%f\n", idx[0][0], idx[0][1], idx[0][2], idx[0][3], idx[0][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idx[1][0], idx[1][1], idx[1][2], idx[1][3], idx[1][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idx[2][0], idx[2][1], idx[2][2], idx[2][3], idx[2][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idx[3][0], idx[3][1], idx[3][2], idx[3][3], idx[3][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idx[4][0], idx[4][1], idx[4][2], idx[4][3], idx[4][4]);
      terminate("Error in characteristic limiter!\n");
    }

  if(!is_id(idy))
    {
      printf("Ly*Ry!=Id\n");

      printf("id:\n");
      printf("%f\t\%f\t%f\t%f\t%f\n", idy[0][0], idy[0][1], idy[0][2], idy[0][3], idy[0][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idy[1][0], idy[1][1], idy[1][2], idy[1][3], idy[1][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idy[2][0], idy[2][1], idy[2][2], idy[2][3], idy[2][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idy[3][0], idy[3][1], idy[3][2], idy[3][3], idy[3][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idy[4][0], idy[4][1], idy[4][2], idy[4][3], idy[4][4]);
      terminate("Error in characteristic limiter!\n");
    }

#ifndef TWODIMS
  if(!is_id(idz))
    {
      printf("Lz*Rz!=Id\n");

      printf("id:\n");
      printf("%f\t\%f\t%f\t%f\t%f\n", idz[0][0], idz[0][1], idz[0][2], idz[0][3], idz[0][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idz[1][0], idz[1][1], idz[1][2], idz[1][3], idz[1][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idz[2][0], idz[2][1], idz[2][2], idz[2][3], idz[2][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idz[3][0], idz[3][1], idz[3][2], idz[3][3], idz[3][4]);
      printf("%f\t\%f\t%f\t%f\t%f\n", idz[4][0], idz[4][1], idz[4][2], idz[4][3], idz[4][4]);
      terminate("Error in characteristic limiter!\n");
    }
#endif

#endif

  return;
}
#endif

#endif //endif DG
