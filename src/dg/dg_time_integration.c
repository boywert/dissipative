/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_time_integration.c
 * \date        10/2014
 * \author    Kevin Schaal
 * \brief     Excplicit Runge Kutta time integration
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"


#ifdef DG


//Butcher tableau coefficients

static const int _rk_stages = RK_STAGES;
static double _c[RK_STAGES];
static double _b[RK_STAGES];
static double _a[RK_STAGES][RK_STAGES];


void ini_rk_coefficients()
{

#ifdef RK1                      // explicit Euler
  _c[0] = 0;

  _b[0] = 1;
#endif

#ifdef RK2                      // Heun's method
  _c[0] = 0;

  _c[1] = 1;
  _a[1][0] = 1;

  _b[0] = 0.5;
  _b[1] = 0.5;
#endif

#ifdef RK3                      // SSP RK3
  _c[0] = 0;

  _c[1] = 1;
  _a[1][0] = 1.;

  _c[2] = 0.5;
  _a[2][0] = 1. / 4.;
  _a[2][1] = 1. / 4.;

  _b[0] = 1. / 6.;
  _b[1] = 1. / 6.;
  _b[2] = 2. / 3.;
#endif

#ifdef RK4                      // SSP RK4

  _c[0] = 0;

  _c[1] = 0.39175222700392;
  _a[1][0] = 0.39175222700392;

  _c[2] = 0.58607968896779;
  _a[2][0] = 0.21766909633821;
  _a[2][1] = 0.36841059262959;

  _c[3] = 0.47454236302687;
  _a[3][0] = 0.08269208670950;
  _a[3][1] = 0.13995850206999;
  _a[3][2] = 0.25189177424738;

  _c[4] = 0.93501063100924;
  _a[4][0] = 0.06796628370320;
  _a[4][1] = 0.11503469844438;
  _a[4][2] = 0.20703489864929;
  _a[4][3] = 0.54497475021237;

  _b[0] = 0.14681187618661;
  _b[1] = 0.24848290924556;
  _b[2] = 0.10425883036650;
  _b[3] = 0.27443890091960;
  _b[4] = 0.22600748319395;
#endif


  //check Runge Kutta table for validity.
  int error = 0;

  int i, s;
  double sum = 0;
  double eps = 1e-10;

  if(_c[0] != 0)
    {
      error = 1;
    }

  //the weights should sum up to 1
  for(i = 0; i < _rk_stages; i++)
    {
      sum += _b[i];
    }

  if(fabs(sum - 1) > eps)
    {
      error = 1;
    }

  //The RK-Matrix rows should sum up to the nodes
  for(s = 1; s < _rk_stages; s++)
    {
      sum = 0;

      for(i = 0; i < s; i++)
        {
          sum += _a[s][i];
        }

      if(fabs(sum - _c[s]) > eps)
        {
          error = 1;
        }
    }

  if(error)
    {
      terminate("Runge Kutta table not valid!\n");
    }
}

#ifndef OLD_TIME_INTEGRATION

void dg_compute_step()
{
  TIMER_START(CPU_DISCONTINUOUS_GALERKIN);

  mpi_printf("DG: Starting a step\n");

#ifdef DG_DEBUG
  check_time_steps();
#endif

  dg_setup_step();

  Aexch = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc_movable(&(Aexch), "Aexch", Mesh_nimport * Nof_base_functions * 5 * sizeof(double));
  double (*_mRdt)[NumGas][NOF_BASE_FUNCTIONS][5] = (double (*)[NumGas][NOF_BASE_FUNCTIONS][5]) mymalloc("a2", NumGas * _rk_stages * Nof_base_functions * 5 * sizeof(double));
  double (*w0)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("w0", NumGas * Nof_base_functions * 5 * sizeof(double));
  double (*w_temp)[NOF_BASE_FUNCTIONS][5] = (double (*)[NOF_BASE_FUNCTIONS][5]) mymalloc("w0", NumGas * Nof_base_functions * 5 * sizeof(double));
  double (*R_inner)[5] = (double (*)[5]) mymalloc("R_inner", Nof_base_functions * 5 * sizeof(double));


  mpi_printf("DG: Propagation with a %d-stage RK\n", _rk_stages);

  copy_cell_weights_to_array(w0);       //store the original weights


  //slope limit the original weights
  exchange_weights(w0);
#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  //copy_array_to_cell_weights(w0); not needed here
  dg_recompute_nodes();
  exchange_node_data();
#endif
  minmod_limiter(w0);


  int s;
  int i, l, k;
  int m;
  double dt;

  //calculate the r.h.s of the differential equation
  for(s = 0; s < _rk_stages; s++)
    {

      copy_array_to_array(w0, w_temp);

      //calculate argument w for R(t,w) (current solution)
      if(s != 0)                //not for the first stage
        {
          for(i = 0; i < NumGas; i++)
            {
              if(P[i].Mass == 0 && P[i].ID == 0)
                continue;       /* skip dissolved cells */

              for(m = 0; m < s; m++)
                {
                  for(l = 0; l < Nof_base_functions; l++)
                    {
                      for(k = 0; k < 5; k++)
                        {
                          w_temp[i][l][k] += _a[s][m] * _mRdt[m][i][l][k];
                        }
                    }
                }
            }
        }

      //slope limit current solution
      if(s != 0)                //w0 is already slope limited
        {
          exchange_weights(w_temp);
#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
          copy_array_to_cell_weights(w_temp);
          dg_recompute_nodes();
          exchange_node_data();
#endif
          minmod_limiter(w_temp);
        }

      TIMER_START(CPU_DG_INNER);

      mpi_printf("DG: Calculating inner integral\n");

      for(i = 0; i < NumGas; i++)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue;           /* skip dissolved cells */

          dt = (P[i].TimeBinHydro ? (((integertime) 1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;      /* time step of the cell */

          calc_R_inner(w_temp, i, All.Time + _c[s] * dt, dt, R_inner, _b[s]);

          for(l = 0; l < Nof_base_functions; l++)
            {
              for(k = 0; k < 5; k++)
                {
                  _mRdt[s][i][l][k] = -dt * R_inner[l][k];
                }
            }
        }

      TIMER_STOP(CPU_DG_INNER);

      //outer integral
      exchange_weights(w_temp);
      subtract_R_outer(&Mesh, 1, w_temp, _mRdt[s]);     //_mRdt[s] = _mRdt[s] - 1 * dt * R_outer(w_temp)
    }

  //update SphP Weights
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;               /* skip dissolved cells */

      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 5; k++)
            {
              for(s = 0; s < _rk_stages; s++)
                {
                  w0[i][l][k] += _b[s] * _mRdt[s][i][l][k];
                }
            }
        }
    }

  exchange_weights(w0);
#if defined(REFINEMENT_SPLIT_CELLS) || defined(REFINEMENT_MERGE_CELLS)
  copy_array_to_cell_weights(w0);
  dg_recompute_nodes();
  exchange_node_data();
#endif
  minmod_limiter(w0);

  copy_array_to_cell_weights(w0);

  dg_update_conserved_variables();

  dg_end_step();

  //free memory
  myfree(R_inner);
  myfree(w_temp);
  myfree(w0);
  myfree(_mRdt);
  myfree(Aexch);

  mpi_printf("DG: step done!\n");

  TIMER_STOP(CPU_DISCONTINUOUS_GALERKIN);
}

#endif

/*!
 * Print which integration scheme is used
 */

void print_time_integration_info()
{
  int counter = 0;

#ifdef RK1
  mpi_printf("Using RK1 (explicit Euler) time integration!\n");
  counter++;
#endif

#ifdef RK2
  mpi_printf("Using explicit Runge-Kutta 2 time integration with alpha=%f!\n", All.DG_RK2_alpha);
  counter++;
#endif

#ifdef RK3
  mpi_printf("Using explicit Runge-Kutta 3 time integration!\n");
  counter++;
#endif

#ifdef RK4
  mpi_printf("Using explicit Runge-Kutta 4 time integration!\n");
  counter++;
#endif

  if(counter != 1)
    {
      terminate("Time integrator not specified correctly!\n");
    }

}


#endif
