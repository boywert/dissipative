/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_core.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		The main functions of the Discontinuous Galerkin (DG) module
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"
#include <string.h>
#include "dg_core_inline.h"

//3d basic features done


#ifdef DG


/*!
 * Reserve memory, calculate quadrature data and base functions
 */
void dg_initialize()
{                               //todo: movable

  mpi_printf("\n\n===== THIS IS TENET =====\n\n");

  //checks
  //

#ifndef FORCE_EQUAL_TIMESTEPS
  terminate("ERROR in dg_initialize: Compile with FORCE_EQUAL_TIMESTEPS!\n");
#endif

#if (defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)) && (DEGREE_K==0)
  terminate("ERROR in dg_initialize: For DG AMR use a higher degree than DEGREE_K=0\n");
#endif


  //set global non-const variables and reserve memory
  //

  Nof_lobatto_points_1d = NOF_LOBATTO_POINTS_1D;

  Lobatto_points_1d = (double *) mymalloc("Lobatto_points_1d", Nof_lobatto_points_1d * sizeof(double));
  Lobatto_points_weights_1d = (double *) mymalloc("Lobatto_points_weights_1d", Nof_lobatto_points_1d * sizeof(double));

  Quad_points_1d = (double *) mymalloc("Quad_points_1d", Nof_quad_points_1d * sizeof(double));
  Quad_points_weights_1d = (double *) mymalloc("Quad_points_weights_1d", Nof_quad_points_1d * sizeof(double));

  Inner_quad_points_x = (double *) mymalloc("Inner_quad_points_x", Nof_inner_quad_points * sizeof(double));
  Inner_quad_points_y = (double *) mymalloc("Inner_quad_points_y", Nof_inner_quad_points * sizeof(double));
  Inner_quad_points_z = (double *) mymalloc("Inner_quad_points_y", Nof_inner_quad_points * sizeof(double));
  Inner_quad_points_weights = (double *) mymalloc("Inner_quad_points_weights", Nof_inner_quad_points * sizeof(double));

  Outer_quad_points_x = (double (*)[NOF_OUTER_QUAD_POINTS]) mymalloc("Outer_quad_points_x", 2 * NUMDIMS * Nof_outer_quad_points * sizeof(double));
  Outer_quad_points_y = (double (*)[NOF_OUTER_QUAD_POINTS]) mymalloc("Outer_quad_points_y", 2 * NUMDIMS * Nof_outer_quad_points * sizeof(double));
  Outer_quad_points_z = (double (*)[NOF_OUTER_QUAD_POINTS]) mymalloc("Outer_quad_points_z", 2 * NUMDIMS * Nof_outer_quad_points * sizeof(double));
  Outer_quad_points_weights = (double (*)[NOF_OUTER_QUAD_POINTS]) mymalloc("Outer_quad_points_weights", 2 * NUMDIMS * Nof_outer_quad_points * sizeof(double));

  Inner_base_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Inner_base_values", Nof_base_functions * Nof_inner_quad_points * sizeof(double));
  Inner_base_dx_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Inner_base_dx_values", Nof_base_functions * Nof_inner_quad_points * sizeof(double));
  Inner_base_dy_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Inner_base_dy_values", Nof_base_functions * Nof_inner_quad_points * sizeof(double));
  Inner_base_dz_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Inner_base_dz_values", Nof_base_functions * Nof_inner_quad_points * sizeof(double));

  Outer_base_values = (double (*)[NOF_OUTER_QUAD_POINTS][NOF_BASE_FUNCTIONS]) mymalloc("Outer_base_values", Nof_base_functions * 2 * NUMDIMS * Nof_outer_quad_points * sizeof(double));
  Outer_fine_base_values =
    (double (*)[NOF_OUTER_FINE_QUAD_POINTS][NOF_BASE_FUNCTIONS]) mymalloc("Outer_fine_base_values", NOF_INTERFACES * NOF_OUTER_FINE_QUAD_POINTS * Nof_base_functions * sizeof(double));

  P_A = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_A", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_B = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_B", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_C = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_C", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_D = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_D", Nof_base_functions * Nof_base_functions * sizeof(double));

#ifndef TWODIMS
  P_E = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_E", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_F = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_F", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_G = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_G", Nof_base_functions * Nof_base_functions * sizeof(double));
  P_H = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("P_H", Nof_base_functions * Nof_base_functions * sizeof(double));
#endif

#ifdef CHARACTERISTIC_LIMITER
  Lx = (double (*)[5]) mymalloc("Lx", 5 * 5 * sizeof(double));
  Ly = (double (*)[5]) mymalloc("Ly", 5 * 5 * sizeof(double));
  DG3D(Lz = (double (*)[5]) mymalloc("Lz", 5 * 5 * sizeof(double));) Rx = (double (*)[5]) mymalloc("Rx", 5 * 5 * sizeof(double));
  Ry = (double (*)[5]) mymalloc("Ry", 5 * 5 * sizeof(double));
  DG3D(Rz = (double (*)[5]) mymalloc("Rz", 5 * 5 * sizeof(double));)
#endif
    Temp_weights = (double (*)[5]) mymalloc("Temp_weights", Nof_base_functions * 5 * sizeof(double));

#if (NUMDIMS==2)
  Nof_union_quad_points = (NOF_LOBATTO_POINTS_1D * NOF_QUAD_POINTS_1D * 2);
#else
  Nof_union_quad_points = (NOF_LOBATTO_POINTS_1D * NOF_QUAD_POINTS_1D * NOF_QUAD_POINTS_1D * 3);
#endif

  Union_quad_points_x = (double *) mymalloc("Union_quad_points_x", Nof_union_quad_points * sizeof(double));
  Union_quad_points_y = (double *) mymalloc("Union_quad_points_y", Nof_union_quad_points * sizeof(double));
  Union_quad_points_z = (double *) mymalloc("Union_quad_points_z", Nof_union_quad_points * sizeof(double));
  Union_base_values = (double (*)[NOF_BASE_FUNCTIONS]) mymalloc("Union_base_values", Nof_base_functions * Nof_union_quad_points * sizeof(double));

  //initialize base functions and quadrature data
  //

  ini_1d_quadrature_points();
  ini_multi_d_quadrature_points();

  ini_base_function_values();
  ini_refinement_matrices();

  print_quad_info();
  print_base_info();
  print_ref_matrix_info();
  print_time_integration_info();


  //initialize time integration
  ini_rk_coefficients();

  //initialize higher order quadrature
  ini_higher_order_quadrature();
}

/*!
 * Free memory and clean up
 * Note: this function is never called since AREPO is finished with exit(0).
 */
void dg_finalize()
{
  myfree(Union_base_values);
  myfree(Union_quad_points_z);
  myfree(Union_quad_points_y);
  myfree(Union_quad_points_x);

  myfree(Temp_weights);
#ifdef CHARACTERISTIC_LIMITER
  DG3D(myfree(Rz);) myfree(Ry);
  myfree(Rx);
  DG3D(myfree(Lz);) myfree(Ly);
  myfree(Lx);
#endif
#ifndef TWODIMS
  myfree(P_H);
  myfree(P_G);
  myfree(P_F);
  myfree(P_E);
#endif
  myfree(P_D);
  myfree(P_C);
  myfree(P_B);
  myfree(P_A);
  myfree(Outer_fine_base_values);
  myfree(Outer_base_values);
  myfree(Inner_base_dy_values);
  myfree(Inner_base_dx_values);
  myfree(Inner_base_values);
  myfree(Union_base_values);
  myfree(Outer_quad_points_weights);
  myfree(Outer_quad_points_y);
  myfree(Outer_quad_points_x);
  myfree(Inner_quad_points_weights);
  myfree(Inner_quad_points_y);
  myfree(Inner_quad_points_x);
  myfree(Quad_points_weights_1d);
  myfree(Quad_points_1d);
  myfree(Lobatto_points_weights_1d);
  myfree(Lobatto_points_1d);
}

#ifdef DG_TEST_PROBLEM
/*!
 * Set the initial conditions from density,pressure and velocity specified in dg_test_problems.c
 * The projection integral is solved with the midpoint rule.
 */
void dg_set_initial_conditions()
{
  DG_PRINTF("Setting initial conditions...\n");

#ifndef MESHRELAX_DENSITY_IN_INPUT
  terminate("Error in dg_set_initial_conditions: compile with MESHRELAX_DENSITY_IN_INPUT!\n");
#endif

  double w[5];
  double pv[5];

  MyDouble *w0l, *w1l, *w2l, *w3l, *w4l;

  //absolute location in the box
  double abs_x, abs_y, abs_z;

  //cell length
  double dl;

#ifdef INI_MIDPOINT_RULE
#ifndef TWODIMS
#error Not yet implemented in 3D!
#endif

  double hx = 2. / MIDPOINT_RULE_STEPS;
  double hy = 2. / MIDPOINT_RULE_STEPS;
  double hh = hx * hy;

  double base_function_value;

  //location in [-1,1]^2
  double rel_x, rel_y;

  int i;

#endif

  int j, k, l;

  for(k = 0; k < NumGas; k++)
    {

      dl = amr_length[Mesh.DP[k].level];


      for(l = 0; l < Nof_base_functions; l++)
        {

          w0l = &SphP[k].Weights[l][0];
          w1l = &SphP[k].Weights[l][1];
          w2l = &SphP[k].Weights[l][2];
          w3l = &SphP[k].Weights[l][3];
          w4l = &SphP[k].Weights[l][4];

          *w0l = 0;
          *w1l = 0;
          *w2l = 0;
          *w3l = 0;
          *w4l = 0;

#ifdef INI_MIDPOINT_RULE
#ifndef TWODIMS
#error Not yet implemented in 3D!
#endif

          //projection (integral)
          for(i = 0; i < MIDPOINT_RULE_STEPS; i++)      //loop in x-direction
            {
              for(j = 0; j < MIDPOINT_RULE_STEPS; j++)  //loop in y-direction
                {
                  rel_x = -1 + hx * (i + 0.5);  //in [-1,1]
                  rel_y = -1 + hy * (j + 0.5);  //in [-1,1]

                  abs_x = 0.5 * rel_x * dl + P[k].Pos[0];
                  abs_y = 0.5 * rel_y * dl + P[k].Pos[1];

                  pv[RHO] = ic_density(abs_x, abs_y);
                  pv[PRES] = ic_pressure(abs_x, abs_y);
                  pv[VX] = ic_velocity_x(abs_x, abs_y);
                  pv[VY] = ic_velocity_y(abs_x, abs_y);

                  primvars_to_w(pv, w);

                  base_function_value = base_function(l, rel_x, rel_y);

                  *w0l += w[0] * base_function_value * hh;
                  *w1l += w[1] * base_function_value * hh;
                  *w2l += w[2] * base_function_value * hh;
                  *w3l += w[3] * base_function_value * hh;
                }
            }
#else //gaussian quadrature

          double pos_abs[3];

          for(j = 0; j < Nof_inner_quad_points; j++)    //loop over quadrature points
            {
              abs_pos_inner_quad_points(j, dl, P[k].Pos, pos_abs);

              abs_x = pos_abs[0];
              abs_y = pos_abs[1];
              abs_z = pos_abs[2];

              pv[RHO] = ic_density(abs_x, abs_y, abs_z);
              pv[PRES] = ic_pressure(abs_x, abs_y, abs_z);
              pv[VX] = ic_velocity_x(abs_x, abs_y, abs_z);
              pv[VY] = ic_velocity_y(abs_x, abs_y, abs_z);
              pv[VZ] = ic_velocity_z(abs_x, abs_y, abs_z);

              primvars_to_w(pv, w);

              *w0l += w[0] * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
              *w1l += w[1] * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
              *w2l += w[2] * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
              *w3l += w[3] * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
              *w4l += w[4] * GET_Inner_base_values(j, l) * GET_Inner_quad_points_weights(j);
            }
#endif

#ifdef TWODIMS
          const double norm_fac = 4;
#else
          const double norm_fac = 8;
#endif
          *w0l /= norm_fac;
          *w1l /= norm_fac;
          *w2l /= norm_fac;
          *w3l /= norm_fac;
          *w4l /= norm_fac;
        }

      //set values of SphP to cell averages (w[...][0])

      w[0] = SphP[k].Weights[0][0];
      w[1] = SphP[k].Weights[0][1];
      w[2] = SphP[k].Weights[0][2];
      w[3] = SphP[k].Weights[0][3];
      w[4] = SphP[k].Weights[0][4];


      w_to_primvars(w, pv);

      SphP[k].Utherm = ideal_gas_utherm(pv[RHO], pv[PRES]);
      P[k].Mass = pv[RHO];
      P[k].Vel[0] = pv[VX];
      P[k].Vel[1] = pv[VY];
      P[k].Vel[2] = pv[VZ];

#ifdef DISCONTINUITY_DETECTION
#ifndef TWODIMS
#error Not yet implemented in 3D!
#endif

      //set inflow boundaries
      int e, q;
      double sum;
      double w[4];
      int flux_component;

      SphP[k].Inflow_boundaries = 0;

      for(e = 0; e < 4; e++)    //loop over interfaces
        {
          sum = 0;

          if(e == LEFT || e == RIGHT)
            {
              flux_component = 1;
            }
          else
            {
              flux_component = 2;
            }

          for(q = 0; q < Nof_outer_quad_points; q++)    //integrate the mass flux along the interface
            {
              calc_state_at_outer_quad_point(SphP[k].Weights, e, q, w, 0);

              sum += 0.5 * w[flux_component] * GET_Outer_quad_points_weights(e, q);;
            }

          switch (e)
            {
            case LEFT:
              if(sum > 0)
                SphP[k].Inflow_boundaries += (1 << e);
              break;

            case RIGHT:
              if(sum < 0)
                SphP[k].Inflow_boundaries += (1 << e);
              break;

            case BACK:
              if(sum > 0)
                SphP[k].Inflow_boundaries += (1 << e);
              break;

            case FRONT:
              if(sum < 0)
                SphP[k].Inflow_boundaries += (1 << e);
              break;

            default:
              terminate("interface not valid!\n");
              break;
            }

        }                       //end loop over interface
#endif //end DISCONTINUITY_DETECTION

#ifdef POSITIVITY_LIMITER
      positivity_limiter(SphP[k].Weights, k);
#endif
    }                           //end loop over gas cells

  //slope limit the initial conditions
  Aexch = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc_movable(&(Aexch), "Aexch", Mesh_nimport * Nof_base_functions * 5 * sizeof(double));
  double (*a1)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("a1", NumGas * Nof_base_functions * 5 * sizeof(double));

  copy_cell_weights_to_array(a1);

  exchange_weights(a1);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a1);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a1);

  copy_array_to_cell_weights(a1);

  myfree(a1);

  myfree(Aexch);

  DG_PRINTF("Initial conditions set!\n");
}
#endif

#ifdef OLD_TIME_INTEGRATION

/*!
 *  Compute a full Discontinuous Galerkin step with Runge Kutta 3 and update the conserved variables
 */
void dg_compute_step()
{

  TIMER_START(DISCONTINUOUS_GALERKIN);

  mpi_printf("DG: Starting a step\n");

#ifdef DG_DEBUG
  check_time_steps();
#endif

  dg_setup_step();

  //memory reservation
  Aexch = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc_movable(&(Aexch), "Aexch", Mesh_nimport * Nof_base_functions * 5 * sizeof(double));
  double (*w0)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("w0", NumGas * Nof_base_functions * 5 * sizeof(double));
  double (*a1)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("a1", NumGas * Nof_base_functions * 5 * sizeof(double));
  double (*a2)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("a2", NumGas * Nof_base_functions * 5 * sizeof(double));

  double (*R_inner)[5] = (double (*)[5]) mymalloc("R_inner", Nof_base_functions * 5 * sizeof(double));

  copy_cell_weights_to_array(w0);       //store the original weights

  //slope limit the original weights
  exchange_weights(w0);
#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  //copy_array_to_cell_weights(w0); not needed here
  dg_recompute_nodes();
  exchange_node_data();
#endif
  minmod_limiter(w0);


#ifdef RK1

  mpi_printf("DG: Propagation with an Euler step\n");

  copy_array_to_array(w0, a1);  // a1 = w0

  int i, k, l;
  double dt;

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a1, i, All.Time, dt, R_inner, 1);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a2[i][l][k] = a1[i][l][k] - dt * R_inner[l][k];   // a_2 = a_1 - dt * R_inner(a1)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a1);

  subtract_R_outer(&Mesh, 1, a1, a2);   //a2 = a2 - 1 * dt * R_outer(a1)

  exchange_weights(a2);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a2);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a2);

  copy_array_to_cell_weights(a2);

#elif defined(RK2)              //RK2 (Wikipedia: Runge-Kutta methods)

  mpi_printf("DG: Propagation with a RK2 step\n");

  double factor = (1 - 0.5 / All.DG_RK2_alpha) / All.DG_RK2_alpha;      // (1-(1/(2alpha))/alpha

  copy_array_to_array(w0, a1);  // a1 = w0

  int i, k, l;
  double dt;

  // a_2 = a1 - alpha * dt * R(a1)
  //

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a1, i, All.Time, dt, R_inner, 0.5);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a2[i][l][k] = a1[i][l][k] - All.DG_RK2_alpha * dt * R_inner[l][k];        // a_2 = a_1 - alpha * dt * R_inner(a1)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a1);

  subtract_R_outer(&Mesh, All.DG_RK2_alpha, a1, a2);    //a2 = a2 - alpha * dt * R_outer(a1)

  exchange_weights(a2);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a2);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a2);

  //
  ////


  // a_1 = ((1/(2alpha)-1)/alpha+1) * w_0 - (1/(2alpha)-1)/alpha * a2 - 1/(2alpha) * dt * R(a2)
  //

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a2, i, All.Time + dt, dt, R_inner, 0.5);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a1[i][l][k] = (1 - factor) * w0[i][l][k] + factor * a2[i][l][k] - 0.5 / All.DG_RK2_alpha * dt * R_inner[l][k];    // a_1 = (1-(1-1/(2alpha))/alpha) * w_0 + (1-1/(2alpha))/alpha * a2 - 1/(2alpha) * dt * R_inner(a2)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a2);

  subtract_R_outer(&Mesh, 0.5 / All.DG_RK2_alpha, a2, a1);      //a1 = a1 - 1/(2alpha) * dt * R_outer(a2)

  exchange_weights(a1);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a1);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a1);

  //
  ////

  copy_array_to_cell_weights(a1);


#elif defined(RK3)              //RK3

  mpi_printf("DG: Propagation with a RK3 step\n");

  copy_array_to_array(w0, a1);  // a1 = w0

  int i, k, l;
  double dt;


  // a_2 = a1 - 1 * dt * R(a1)
  //

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a1, i, All.Time, dt, R_inner, 1. / 6.);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a2[i][l][k] = a1[i][l][k] - dt * R_inner[l][k];   // a_2 = a_1 - dt * R_inner(a1)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a1);

  subtract_R_outer(&Mesh, 1, a1, a2);   //a2 = a2 - 1 * dt * R_outer(a1)

  exchange_weights(a2);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a2);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a2);

  //
  ////

  // a_1 = 3./4. * w_0 + 1./4. * a2  - 1./4. * dt * R(a2)
  //

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a2, i, All.Time + dt, dt, R_inner, 1. / 6.);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {

              a1[i][l][k] = 3. / 4. * w0[i][l][k] + 1. / 4. * a2[i][l][k] - 1. / 4. * dt * R_inner[l][k];       //  a_1 = 3./4. * w_0 + 1./4. * a2  - 1./4. * dt * R_inner(a2)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a2);

  subtract_R_outer(&Mesh, 1. / 4., a2, a1);     //a1 = a1 - 1./4. * dt * R_outer(a2)

  exchange_weights(a1);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a1);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a1);

  //
  ////



  // a_2 = 1./3. * w_0 + 2./3. * a_1 - 2./3. * dt * R(a1)
  //

  TIMER_START(DG_INNER);
  mpi_printf("DG: Calculating inner integral\n");

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;  /* time step of the cell */

      calc_R_inner(a1, i, All.Time + 0.5 * dt, dt, R_inner, 2. / 3.);

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a2[i][l][k] = 1. / 3. * w0[i][l][k] + 2. / 3. * a1[i][l][k] - 2. / 3. * dt * R_inner[l][k];       //  a_2 = 1./3. * w_0 + 2./3. * a_1 - 2./3. * dt * R_inner(a1)
            }
        }
    }

  TIMER_STOP(DG_INNER);

  exchange_weights(a1);

  subtract_R_outer(&Mesh, 2. / 3., a1, a2);     //a2 = a2 - 2./3. * dt * R_outer(a1)

  exchange_weights(a2);

#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(a2);
  dg_recompute_nodes();
  exchange_node_data();
#endif

  minmod_limiter(a2);

  //
  ////

  copy_array_to_cell_weights(a2);
#else

#error Specify a time integrator!

#endif

  dg_update_conserved_variables();

  dg_end_step();

  //free memory
  myfree(R_inner);
  myfree(a2);
  myfree(a1);
  myfree(w0);
  myfree(Aexch);

  mpi_printf("DG: step done!\n");

  TIMER_STOP(DISCONTINUOUS_GALERKIN);
}

#endif

/*!
 * set values before a dg step
 */
void dg_setup_step()
{

}

/*!
 * end a dg step
 */

#ifdef DG_CON_VARS_SUM_TO_FILE
static double s_time_next_to_file = 0;
#endif

void dg_end_step()
{
#ifdef DG_CON_VARS_SUM_TO_FILE

  if(All.Time >= s_time_next_to_file)
    {

      int i;
      double local_angular_momentum = 0;
      double tot_angular_momentum = 0;
      double local_mass = 0;
      double tot_mass = 0;
      double local_energy = 0;
      double tot_energy = 0;

      for(i = 0; i < NumGas; i++)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;           /* skip dissolved cells */

          local_angular_momentum += angular_momentum(SphP[i].Weights, i);
          local_mass += P[i].Mass;
          local_energy += SphP[i].Energy;
        }

      MPI_Allreduce(&local_angular_momentum, &tot_angular_momentum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_mass, &tot_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&local_energy, &tot_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      if(ThisTask == 0)
        {
          fprintf(FdAngularMomentumDG, "%e\t%e\n", All.Time, tot_angular_momentum);
          fprintf(FdMassDG, "%e\t%e\n", All.Time, tot_mass);
          fprintf(FdEnergyDG, "%e\t%e\n", All.Time, tot_energy);

          myflush(FdAngularMomentumDG);
          myflush(FdMassDG);
          myflush(FdEnergyDG);
        }

      s_time_next_to_file += DG_CON_VARS_SUM_TO_FILE;
    }

#endif
}


/*!
 * Calculate the contribution of the inner Integral to the change of the weights of a cell
 */
void calc_R_inner(CBV(a), int cell, double time, double dt, double R_inner[NOF_BASE_FUNCTIONS][5], double rk_weight)
{
  int k, l, q;

  double w[5];                  //state in specific cell at specific quad point (rho, rho*vx, rho*vy, E_density)
  double f1[5];                 //flux in x-direction
  double f2[5];                 //flux in y-direction
  DG3D(double f3[5];)           //flux in z-direction
#ifdef DG_EXTERNAL_ACCELERATION
  double src_term[5];
#endif

  double dl = amr_length[Mesh.DP[cell].level];

#ifdef TWODIMS
  double volume = dl * dl;
  double norm_fac = dl / 2. / volume;   //todo: change to SphP[cell].Volume
#else
  double volume = dl * dl * dl;
  double norm_fac = dl * dl / 4. / volume;
#endif

  memset(R_inner, 0, sizeof(double) * NOF_BASE_FUNCTIONS * 5);

#ifdef DG_EXTERNAL_ACCELERATION
  double pv[5];
  double xlab, ylab, zlab;
#ifdef TWODIMS
  double norm_fac_grav = 1. / 4.;
#else
  double norm_fac_grav = 1. / 8;
#endif
#endif


  //contributions from interior quad points (first & second term)
  for(q = 0; q < Nof_inner_quad_points; q++)
    {
      //calculate conserved quantities at quadrature point q
      calc_state_at_quad_point(a[cell], q, w);

      //calculate the flux
      flux1(w, f1);
      flux2(w, f2);
      DG3D(flux3(w, f3);)
#ifdef DG_EXTERNAL_ACCELERATION
        cell_coords_to_lab_coords(cell, GET_Inner_quad_points_x(q), GET_Inner_quad_points_y(q), GET_Inner_quad_points_z(q), &xlab, &ylab, &zlab);

      double acc[3];
      dg_acceleration(xlab, ylab, zlab, acc);

      w_to_primvars(w, pv);

#ifdef DG_TURBULENCE
      SphP[cell].EgyDrive += pv[RHO] * (pv[VX] * acc[0] + pv[VY] * acc[1] + pv[VZ] * acc[2]) * rk_weight * dt * SphP[cell].Volume * norm_fac_grav * GET_Inner_quad_points_weights(q);
#endif
#endif

      //loop over base functions
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              R_inner[l][k] += -norm_fac * f1[k] * GET_Inner_base_dx_values(q, l) * GET_Inner_quad_points_weights(q);   //term1
              R_inner[l][k] += -norm_fac * f2[k] * GET_Inner_base_dy_values(q, l) * GET_Inner_quad_points_weights(q);   //term2
              DG3D(R_inner[l][k] += -norm_fac * f3[k] * GET_Inner_base_dz_values(q, l) * GET_Inner_quad_points_weights(q);)     //term3
            }

#ifdef DG_EXTERNAL_ACCELERATION
          src_term[W_RHO] = 0;
          src_term[W_PX] = pv[RHO] * acc[0] * GET_Inner_base_values(q, l) * GET_Inner_quad_points_weights(q);
          src_term[W_PY] = pv[RHO] * acc[1] * GET_Inner_base_values(q, l) * GET_Inner_quad_points_weights(q);
          src_term[W_PZ] = pv[RHO] * acc[2] * GET_Inner_base_values(q, l) * GET_Inner_quad_points_weights(q);
          src_term[W_E] = pv[RHO] * (pv[VX] * acc[0] + pv[VY] * acc[1] + pv[VZ] * acc[2]) * GET_Inner_base_values(q, l) * GET_Inner_quad_points_weights(q);

          for(k = 0; k < 5; k++)
            {
              R_inner[l][k] -= norm_fac_grav * src_term[k];
            }
#endif
        }
    }
}

/*!
 * \brief returns the type of the neighbour
 */
enum E_ngb get_neighbour_type(int index_ngb)
{
  if(index_ngb < Ngb_MaxPart)   //local cell
    {
      return LOCAL_CELL;
    }
  else if(index_ngb < Ngb_MaxPart + Ngb_MaxNodes)       //local node
    {
      return LOCAL_NODE;
    }
  else if(index_ngb < Ngb_MaxPart + Mesh.nodes_total)   //ghost node
    {
      return GHOST_NODE;
    }
  else                          //ghost cell
    {
      return GHOST_CELL;
    }
}

/*!
 * Returns the weights of a neighbouring cell, corresponding to base function l
 * \param dp_center index of the central cell
 * \param neighbour indicates the neighbour: LEFT, RIGHT, ...
 * \param base_fct index of the base function (0: weights for cell average, 1,2,(3): weights for the slopes)
 * \param CBV(a) degrees of freedom of all cells
 * \param w output parameter, weights of neighbouring cell is stored
 */
void get_neighbour_weights(int dp_center, enum E_interface neighbour, int base_fct, CBV(a), double w[5])
{
  assert(neighbour <= 5 && neighbour >= 0);

  int index, index_new, j;
  index = Mesh.DP[dp_center].neighbors[neighbour];

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  if(index < Ngb_MaxPart)       //local cell
    {
      if(Mesh.DP[dp_center].level == Mesh.DP[index].level)
        {
          w[0] = a[index][base_fct][0];
          w[1] = a[index][base_fct][1];
          w[2] = a[index][base_fct][2];
          w[3] = a[index][base_fct][3];
          w[4] = a[index][base_fct][4];
        }
      else                      //calculate a L2 projection
        {
          assert(Mesh.DP[dp_center].level > Mesh.DP[index].level);

          double (*P_X)[NOF_BASE_FUNCTIONS];

          assign_projection_matrix(dp_center, index, &P_X);

          for(j = 0; j < Nof_base_functions; j++)
            {
              w[0] += a[index][j][0] * P_X[j][base_fct];
              w[1] += a[index][j][1] * P_X[j][base_fct];
              w[2] += a[index][j][2] * P_X[j][base_fct];
              w[3] += a[index][j][3] * P_X[j][base_fct];
              w[4] += a[index][j][4] * P_X[j][base_fct];
            }
        }
    }
  else if(index < Ngb_MaxPart + Ngb_MaxNodes)   //local node
    {
      w[0] = Ngb_Nodes[index].hydro.Weights[base_fct][0];
      w[1] = Ngb_Nodes[index].hydro.Weights[base_fct][1];
      w[2] = Ngb_Nodes[index].hydro.Weights[base_fct][2];
      w[3] = Ngb_Nodes[index].hydro.Weights[base_fct][3];
      w[4] = Ngb_Nodes[index].hydro.Weights[base_fct][4];
    }
  else if(index < Ngb_MaxPart + Mesh.nodes_total)       //ghost node
    {
      w[0] = Ngb_Nodes[index].hydro.Weights[base_fct][0];
      w[1] = Ngb_Nodes[index].hydro.Weights[base_fct][1];
      w[2] = Ngb_Nodes[index].hydro.Weights[base_fct][2];
      w[3] = Ngb_Nodes[index].hydro.Weights[base_fct][3];
      w[4] = Ngb_Nodes[index].hydro.Weights[base_fct][4];
    }
  else                          //ghost cell
    {
      index_new = Mesh.DP[index - Mesh.nodes_total].index;

      assert(index_new >= 0 && index_new < Mesh_nimport);

      if(Mesh.DP[dp_center].level == Mesh.DP[index - Mesh.nodes_total].level)
        {
          w[0] = Aexch[index_new][base_fct][0];
          w[1] = Aexch[index_new][base_fct][1];
          w[2] = Aexch[index_new][base_fct][2];
          w[3] = Aexch[index_new][base_fct][3];
          w[4] = Aexch[index_new][base_fct][4];
        }
      else                      //calculate a L1 projection
        {
          assert(Mesh.DP[dp_center].level > Mesh.DP[index - Mesh.nodes_total].level);

          double (*P_X)[NOF_BASE_FUNCTIONS];

          assign_projection_matrix(dp_center, index - Mesh.nodes_total, &P_X);

          for(j = 0; j < Nof_base_functions; j++)
            {
              w[0] += Aexch[index_new][j][0] * P_X[j][base_fct];
              w[1] += Aexch[index_new][j][1] * P_X[j][base_fct];
              w[2] += Aexch[index_new][j][2] * P_X[j][base_fct];
              w[3] += Aexch[index_new][j][3] * P_X[j][base_fct];
              w[4] += Aexch[index_new][j][4] * P_X[j][base_fct];
            }
        }
    }

//handle non-periodic boundaries, set neighbour weights to zero (except the mean values)

#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)

  if(index != dp_center)        //no boundary cell
    {
      return;
    }

  assert(index < Ngb_MaxPart);

  if(base_fct > 0)
    {
      w[0] = 0;
      w[1] = 0;
      w[2] = 0;
      w[3] = 0;
      w[4] = 0;
    }
#endif

}

/*!
 * Given a weights array of a cell and a quadrature point, calculate the state w
 */
void calc_state_at_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int quad_point, double w[5])
{
  int l;

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  //loop over base functions
  for(l = 0; l < Nof_base_functions; l++)
    {
      w[0] += a_cell[l][0] * GET_Inner_base_values(quad_point, l);
      w[1] += a_cell[l][1] * GET_Inner_base_values(quad_point, l);
      w[2] += a_cell[l][2] * GET_Inner_base_values(quad_point, l);
      w[3] += a_cell[l][3] * GET_Inner_base_values(quad_point, l);
      w[4] += a_cell[l][4] * GET_Inner_base_values(quad_point, l);
    }
}

/*!
 * Given a weights array of a cell and a union quadrature point, calculate the state w
 */
void calc_state_at_union_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int quad_point, double w[5])
{
  int l;

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  //loop over base functions
  for(l = 0; l < Nof_base_functions; l++)
    {
      w[0] += a_cell[l][0] * GET_Union_base_values(quad_point, l);
      w[1] += a_cell[l][1] * GET_Union_base_values(quad_point, l);
      w[2] += a_cell[l][2] * GET_Union_base_values(quad_point, l);
      w[3] += a_cell[l][3] * GET_Union_base_values(quad_point, l);
      w[4] += a_cell[l][4] * GET_Union_base_values(quad_point, l);
    }
}

/*!
 * Given a weights array of a cell and an interface and a quadrature point, calculate the state w
 */
void calc_state_at_outer_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int interface, int quad_point, double w[5], double (*P_X)[NOF_BASE_FUNCTIONS])
{
  int l, j;

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  if(P_X == 0)                  //no projection
    {
      //loop over base functions
      for(l = 0; l < Nof_base_functions; l++)
        {
          w[0] += a_cell[l][0] * GET_Outer_base_values(interface, quad_point, l);
          w[1] += a_cell[l][1] * GET_Outer_base_values(interface, quad_point, l);
          w[2] += a_cell[l][2] * GET_Outer_base_values(interface, quad_point, l);
          w[3] += a_cell[l][3] * GET_Outer_base_values(interface, quad_point, l);
          w[4] += a_cell[l][4] * GET_Outer_base_values(interface, quad_point, l);
        }
    }
  else                          //project bigger cell down, amr boundary
    {
      //calculate projected weights
      for(l = 0; l < Nof_base_functions; l++)
        {
          Temp_weights[l][0] = 0;
          Temp_weights[l][1] = 0;
          Temp_weights[l][2] = 0;
          Temp_weights[l][3] = 0;
          Temp_weights[l][4] = 0;

          for(j = 0; j < Nof_base_functions; j++)
            {
              Temp_weights[l][0] += a_cell[j][0] * P_X[j][l];
              Temp_weights[l][1] += a_cell[j][1] * P_X[j][l];
              Temp_weights[l][2] += a_cell[j][2] * P_X[j][l];
              Temp_weights[l][3] += a_cell[j][3] * P_X[j][l];
              Temp_weights[l][4] += a_cell[j][4] * P_X[j][l];
            }

          w[0] += Temp_weights[l][0] * GET_Outer_base_values(interface, quad_point, l);
          w[1] += Temp_weights[l][1] * GET_Outer_base_values(interface, quad_point, l);
          w[2] += Temp_weights[l][2] * GET_Outer_base_values(interface, quad_point, l);
          w[3] += Temp_weights[l][3] * GET_Outer_base_values(interface, quad_point, l);
          w[4] += Temp_weights[l][4] * GET_Outer_base_values(interface, quad_point, l);
        }
    }
}


/*!
 * Given a weights array of a cell, an interface and a fine quadrature point, calculate the state w
 */
void calc_state_at_outer_fine_quad_point(double a_cell[NOF_BASE_FUNCTIONS][5], int interface, int fine_quad_point, double w[5])
{
  int l;

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  //loop over base functions
  for(l = 0; l < Nof_base_functions; l++)
    {
      w[0] += a_cell[l][0] * GET_Outer_fine_base_values(interface, fine_quad_point, l);
      w[1] += a_cell[l][1] * GET_Outer_fine_base_values(interface, fine_quad_point, l);
      w[2] += a_cell[l][2] * GET_Outer_fine_base_values(interface, fine_quad_point, l);
      w[3] += a_cell[l][3] * GET_Outer_fine_base_values(interface, fine_quad_point, l);
      w[4] += a_cell[l][4] * GET_Outer_fine_base_values(interface, fine_quad_point, l);
    }
}

/*!
 * Calculate the state at an arbitrary position within a cell
 */

void calc_state_at_pos(double a_cell[NOF_BASE_FUNCTIONS][5], double x, double y, double z, double w[5])
{
  int l;

  w[0] = 0;
  w[1] = 0;
  w[2] = 0;
  w[3] = 0;
  w[4] = 0;

  //loop over base functions
  for(l = 0; l < Nof_base_functions; l++)
    {
      w[0] += a_cell[l][0] * base_function(l, x, y, z);
      w[1] += a_cell[l][1] * base_function(l, x, y, z);
      w[2] += a_cell[l][2] * base_function(l, x, y, z);
      w[3] += a_cell[l][3] * base_function(l, x, y, z);
      w[4] += a_cell[l][4] * base_function(l, x, y, z);
    }
}

#ifdef OLD_VERSION


void flux1(double w[5], double f[5])
{
  double u = w[1] / w[0];
  double v = w[2] / w[0];

  f[0] = w[1];
  f[1] = w[4] * (GAMMA - 1) - 0.5 * w[0] * (u * u + v * v) * (GAMMA - 1) + w[0] * u * u;
  f[2] = (w[1] * w[2]) / w[0];
  f[3] = 0;
  f[4] = w[4] * GAMMA * u - 0.5 * w[0] * (u * u + v * v) * u * (GAMMA - 1);
}

/*!
 * calculate the flux1 from a state w
 */
void flux2(double w[5], double f[5])
{
  double u = w[1] / w[0];
  double v = w[2] / w[0];

  f[0] = w[2];
  f[1] = (w[1] * w[2]) / w[0];
  f[2] = w[4] * (GAMMA - 1) - 0.5 * w[0] * (u * u + v * v) * (GAMMA - 1) + w[0] * v * v;
  f[3] = 0;
  f[4] = w[4] * GAMMA * v - 0.5 * w[0] * (u * u + v * v) * v * (GAMMA - 1);
}

void flux3(double w[5], double f[5])
{
  return;
}




#else

/*!
 * calculate the flux1 from a state w
 */
void flux1(double w[5], double f[5])    //speedup: if possible calc all fluxes (1-3) in the same function
{
  double vx = w[W_PX] / w[W_RHO];
  double p = (w[W_E] - 0.5 * (w[W_PX] * w[W_PX] + w[W_PY] * w[W_PY] + w[W_PZ] * w[W_PZ]) / w[W_RHO]) * (GAMMA - 1);

  f[0] = w[W_PX];
  f[1] = p + w[W_PX] * vx;
  f[2] = (w[W_PX] * w[W_PY]) / w[W_RHO];
  f[3] = (w[W_PX] * w[W_PZ]) / w[W_RHO];
  f[4] = (w[W_E] + p) * vx;
}

/*!
 * calculate the flux2 from a state w
 */
void flux2(double w[5], double f[5])
{
  double vy = w[W_PY] / w[W_RHO];
  double p = (w[W_E] - 0.5 * (w[W_PX] * w[W_PX] + w[W_PY] * w[W_PY] + w[W_PZ] * w[W_PZ]) / w[W_RHO]) * (GAMMA - 1);

  f[0] = w[W_PY];
  f[1] = (w[W_PX] * w[W_PY]) / w[W_RHO];
  f[2] = p + w[W_PY] * vy;
  f[3] = (w[W_PY] * w[W_PZ]) / w[W_RHO];
  f[4] = (w[W_E] + p) * vy;
}

/*!
 * calculate the flux1 from a state w
 */
void flux3(double w[5], double f[5])
{
  double vz = w[W_PZ] / w[W_RHO];
  double p = (w[W_E] - 0.5 * (w[W_PX] * w[W_PX] + w[W_PY] * w[W_PY] + w[W_PZ] * w[W_PZ]) / w[W_RHO]) * (GAMMA - 1);

  f[0] = w[W_PZ];
  f[1] = (w[W_PX] * w[W_PZ]) / w[W_RHO];
  f[2] = (w[W_PY] * w[W_PZ]) / w[W_RHO];
  f[3] = p + w[W_PZ] * vz;
  f[4] = (w[W_E] + p) * vz;
}

#endif

/*!
 * Copy the weights of every cell to an array
 */
void copy_cell_weights_to_array(double w[][NOF_BASE_FUNCTIONS][5])
{
  int i, k, l;

  for(i = 0; i < NumGas; i++)   //loop over cells
    {
      if(P[i].Mass == 0 && P[i].ID == 0)        /* dissolved cells */
        {
          for(l = 0; l < Nof_base_functions; l++)       //loop over base functions
            {
              for(k = 0; k < 5; k++)    //loop over variables
                {
                  w[i][l][k] = 0;
                }
            }
        }

      for(l = 0; l < Nof_base_functions; l++)   //loop over base functions
        {
          for(k = 0; k < 5; k++)        //loop over variables
            {
              w[i][l][k] = SphP[i].Weights[l][k];
            }
        }
    }
}

/*!
 * Copy a weights array to the cells
 */
void copy_array_to_cell_weights(double w[][NOF_BASE_FUNCTIONS][5])
{
  int i, k, l;

  for(i = 0; i < NumGas; i++)   //loop over cells
    {
      if(P[i].Mass == 0 && P[i].ID == 0)        /* skip dissolved cells */
        {
          continue;
        }


      for(k = 0; k < 5; k++)    //loop over variables
        {
          for(l = 0; l < Nof_base_functions; l++)       //loop over base functions
            {
              SphP[i].Weights[l][k] = w[i][l][k];
            }
        }
    }
}

/*!
 * Copy a weights array into a second one
 */

void copy_array_to_array(CBV(a), CBV(a_copy))
{
  int i, l, k;

  for(i = 0; i < NumGas; i++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              a_copy[i][l][k] = a[i][l][k];
            }
        }
    }
}


void exchange_weights(CBV(a))
{
  double (*tmpAexch)[NOF_BASE_FUNCTIONS][5];
  int i;
  int place;
  int ngrp, recvTask;
  tessellation *T = &Mesh;

  int k, l;

  TIMER_START(CPU_DG_EXCHANGE) tmpAexch = (double (*)[(NOF_BASE_FUNCTIONS)][5]) mymalloc("tmpAexch", Mesh_nexport * Nof_base_functions * 5 * sizeof(double));

  for(i = 0; i < Mesh_nexport; i++)
    {
      place = T->expList[i].index;

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              tmpAexch[i][l][k] = a[place][l][k];
            }
        }
    }

  /* exchange data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Mesh_Send_count[recvTask] > 0 || Mesh_Recv_count[recvTask] > 0)
            {
              /* exchange the data */
              MPI_Sendrecv(&tmpAexch[Mesh_Send_offset[recvTask]], Mesh_Send_count[recvTask]
                           * Nof_base_functions * 5 * sizeof(double), MPI_BYTE, recvTask, 0,
                           &Aexch[Mesh_Recv_offset[recvTask]], Mesh_Recv_count[recvTask] * Nof_base_functions * 5 * sizeof(double), MPI_BYTE, recvTask, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(tmpAexch);

  TIMER_STOP(CPU_DG_EXCHANGE) return;
}

//helper function for level difference
void assign_projection_matrix(int dp_high_lvl, int dp_low_lvl, double (**P_X)[NOF_BASE_FUNCTIONS])
{
  double highx = Mesh.DP[dp_high_lvl].x;
  double highy = Mesh.DP[dp_high_lvl].y;
  DG3D(double highz = Mesh.DP[dp_high_lvl].z;) double lowx = Mesh.DP[dp_low_lvl].x;
  double lowy = Mesh.DP[dp_low_lvl].y;
  DG3D(double lowz = Mesh.DP[dp_low_lvl].z;) DG3D(if(highz < lowz))
    {
      if(highx < lowx && highy < lowy)
        {
          *P_X = P_A;
          return;
        }
      else if(highx > lowx && highy < lowy)
        {
          *P_X = P_B;
          return;
        }
      else if(highx < lowx && highy > lowy)
        {
          *P_X = P_C;
          return;
        }
      else if(highx > lowx && highy > lowy)
        {
          *P_X = P_D;
          return;
        }
      else
        {
          terminate("ERROR in assign_projection matrix: Could not locate Delaunay points");
        }
    }
#ifndef TWODIMS
  else
    {
      if(highx < lowx && highy < lowy)
        {
          *P_X = P_E;
          return;
        }
      else if(highx > lowx && highy < lowy)
        {
          *P_X = P_F;
          return;
        }
      else if(highx < lowx && highy > lowy)
        {
          *P_X = P_G;
          return;
        }
      else if(highx > lowx && highy > lowy)
        {
          *P_X = P_H;
          return;
        }
      else
        {
          terminate("ERROR in assign_projection matrix: Could not locate Delaunay points");
        }
    }
#endif

  assert(0);
}


/*!
 * Update the conserved variables from the weights
 */
void dg_update_conserved_variables()
{
  int i;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      P[i].Mass = SphP[i].Weights[0][W_RHO] * SphP[i].Volume;
      SphP[i].Momentum[0] = SphP[i].Weights[0][W_PX] * SphP[i].Volume;
      SphP[i].Momentum[1] = SphP[i].Weights[0][W_PY] * SphP[i].Volume;
      SphP[i].Momentum[2] = SphP[i].Weights[0][W_PZ] * SphP[i].Volume;
      SphP[i].Energy = SphP[i].Weights[0][W_E] * SphP[i].Volume;
    }
}

/*!
 * Calculate the time step of a cell
 * \param cell the SphP index
 */
double dg_time_step_cell(int cell)
{
  double c = dg_sound_speed(cell);      //sound speed in the cell
  double dl = amr_length[Mesh.DP[cell].level];  //cell length


#if defined(RK4) && (DG_ORDER > 4)      //reduce timestep
  if(dl < 1)
    {
      dl = pow(dl, DG_ORDER / 4.);
    }
#endif

  double fac = 1. / (2 * DEGREE_K + 1);

#ifdef POSITIVITY_LIMITER
  fac = fmin(fac, 0.5 * GET_Lobatto_points_weights_1d(0));
#endif

  return fac * All.CourantFac * dl / (NUMDIMS * c + fabs(P[cell].Vel[0]) + fabs(P[cell].Vel[1]) + fabs(P[cell].Vel[2]));
}

/*!
 * Calculate the gravity time step of a cell
 * \param cell the SphP index
 */
double dg_time_step_cell_gravity(int cell)
{
#ifdef DG_EXTERNAL_ACCELERATION
  double acc[3];
  dg_acceleration(SphP[cell].Center[0], SphP[cell].Center[1], SphP[cell].Center[2], acc);

  double acc_abs = sqrt(acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]);

  double c = dg_sound_speed(cell);      //sound speed in the cell

  return 1. / sqrt(2 * GAMMA * (GAMMA - 1)) * c / acc_abs;      // src: Zhang's thesis

#else
  terminate("Called dg_time_step_cell_gravity, but DG_EXTERNAL_ACCELERATION is not active!\n");
  return 0;
#endif
}

/*!
 * Calculate the angular momentum of a cell
 */

double angular_momentum(double a_cell[NOF_BASE_FUNCTIONS][5], int cell)
{
#ifndef TWODIMS
  terminate("implement me");
#endif
  double dl = amr_length[Mesh.DP[cell].level];

  //transform cell center to coord system attached to the middle of the box
  double cell_x = SphP[cell].Center[0] - boxHalf_X;
  double cell_y = SphP[cell].Center[1] - boxHalf_Y;

  return spin(a_cell, cell) + dl * dl * (cell_x * a_cell[0][2] - cell_y * a_cell[0][1]);
}

/*!
 * Calculate the spin of a cell
 */

double spin(double a_cell[NOF_BASE_FUNCTIONS][5], int cell)
{
#ifndef TWODIMS
  terminate("implement me");
#endif
  double dl = amr_length[Mesh.DP[cell].level];

  return 1. / sqrt(3.) * dl * dl * dl / 2. * (a_cell[1][2] - a_cell[2][1]);
}


/*!
 * return the opposite interface
 */

int opposite_interface(int e)
{
  switch (e)
    {
    case LEFT:
      return RIGHT;
      break;

    case RIGHT:
      return LEFT;
      break;

    case FRONT:
      return BACK;
      break;

    case BACK:
      return FRONT;
      break;

    case BOTTOM:
      return TOP;
      break;

    case TOP:
      return BOTTOM;
      break;

    default:
      terminate("opposite interface not found!\n");
      return 0;
    }
}

/*!
 * Matrix-Matrix multiplication
 */

void multiply_matrix_matrix(double m1[5][5], double m2[5][5], double m_out[5][5])
{
  int i, j, k;

  for(i = 0; i < 5; i++)
    {
      for(j = 0; j < 5; j++)
        {
          m_out[i][j] = 0;

          for(k = 0; k < 5; k++)
            {
              m_out[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

/*!
 * Matrix-Vector multiplication
 */

void multiply_matrix_vector(double m[5][5], double v[5], double v_out[5])
{
  int i, k;

  for(i = 0; i < 5; i++)
    {
      v_out[i] = 0;

      for(k = 0; k < 5; k++)
        {
          v_out[i] += m[i][k] * v[k];
        }
    }
}

void multiply_matrix_vector_dims(double *m, double *v, double *v_out, int dims)
{
  double (*m_)[dims] = (double (*)[dims]) m;

  int i, k;

  for(i = 0; i < dims; i++)
    {
      v_out[i] = 0;

      for(k = 0; k < dims; k++)
        {
          v_out[i] += m_[i][k] * v[k];
        }
    }
}

/*!
 * stable solving of a quadratic equation (ax^2+bx+c=0)
 */

void solve_quadratic_equation(double a, double b, double c, double *solution_plus, double *solution_minus)
{
  assert(a != 0);

  double p = -b / a;
  double q = c / a;

  if(p > 0)
    {
      if(a > 0)
        {
          *solution_plus = p * 0.5 + sqrt(p * p * 0.25 - q);
          *solution_minus = q / (*solution_plus);
        }
      else
        {
          *solution_minus = p * 0.5 + sqrt(p * p * 0.25 - q);
          *solution_plus = q / (*solution_minus);
        }
    }

  else                          //p<=0
    {
      if(a > 0)
        {
          *solution_minus = p * 0.5 - sqrt(p * p * 0.25 - q);
          *solution_plus = q / (*solution_minus);
        }
      else
        {
          *solution_plus = p * 0.5 - sqrt(p * p * 0.25 - q);
          *solution_minus = q / (*solution_plus);
        }
    }
}

#endif /* DG */
