/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_debug.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Functions for testing and debugging the dg code
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef DG_DEBUG

void cell_info(const char *location)
{

  int i;

  for(i = 0; i < NumGas; i++)
    {
      assert(P[i].ID == Mesh.DP[i].ID);
      if(P[i].ID == 258112)
        {
          printf("w[2]=%g, %s\n", SphP[i].Weights[0][2], location);
          printf("\t index: %d, task: %d\n", i, ThisTask);
          printf("\t position: (%g, %g, %g)\n", SphP[i].Center[0], SphP[i].Center[1], SphP[i].Center[2]);
        }
    }
}

void cell_info_a(const char *location, CBV(a))
{

  int i;

  for(i = 0; i < NumGas; i++)
    {
      if(P[i].ID == 258112)
        {
          printf("w[2]=%g, %s\n", a[i][0][2], location);
          printf("\t index: %d, task: %d\n", i, ThisTask);
          printf("\t position: (%g, %g, %g)\n", SphP[i].Center[0], SphP[i].Center[1], SphP[i].Center[2]);
        }
    }
}

/*!
 * Check whether the matrix is the identity
 */
int is_id(double m[5][5])
{
  int i, j;

  double epsilon = 1e-10;

  for(i = 0; i < 5; i++)
    {
      for(j = 0; j < 5; j++)
        {
          if(i == j && !(m[i][j] < 1 + epsilon && m[i][j] > 1 - epsilon))
            {
              printf("i:%d,j:%d,val:%e\n", i, j, m[i][j]);
              return 0;
            }

          if(i != j && !(m[i][j] < epsilon && m[i][j] > -epsilon))
            {
              printf("i:%d,j:%d,val:%e\n", i, j, m[i][j]);
              return 0;
            }
        }
    }

  return 1;
}

/*!
 * Check the timestep of the cells
 */

void check_time_steps()
{
  int i;
  double dt_hydro = DBL_MAX;
  double dt_gravity = DBL_MAX;

  double dt_hydro_temp;
  double dt_gravity_temp;

  double dt_hydro_global;
  double dt_gravity_global;

  for(i = 0; i < NumGas; i++)
    {
      dt_hydro_temp = dg_time_step_cell(i);

#ifdef DG_EXTERNAL_ACCELERATION
      dt_gravity_temp = dg_time_step_cell_gravity(i);
#else
      dt_gravity_temp = DBL_MAX;
#endif
      if(dt_hydro_temp < dt_hydro)
        {
          dt_hydro = dt_hydro_temp;
        }

      if(dt_gravity_temp < dt_gravity)
        {
          dt_gravity = dt_gravity_temp;
        }
    }

  if(NTask > 1)
    {
      MPI_Allreduce(&dt_hydro, &dt_hydro_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&dt_gravity, &dt_gravity_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
  else
    {
      dt_hydro_global = dt_hydro;
      dt_gravity_global = dt_gravity;
    }

  double time_step = fmin(dt_hydro_global, dt_gravity_global);

  mpi_printf("hydro time step: %.16e, gravity time step: %.16e, time step: %.16e\n", dt_hydro_global, dt_gravity_global, time_step);
}

#endif

#if defined(DG_DEBUG) && defined(DG_DEBUG_OLD)


void print_outer_quad_points_weights()
{
  mpi_printf("outer quadrature points weights:\n");

  int i, j;

  for(i = 0; i < 4; i++)
    {
      for(j = 0; j < Nof_outer_quad_points; j++)
        {
          mpi_printf("%f\n", GET_Outer_quad_points_weights(i, j));
        }
    }
}

void print_weights(int cell)
{
  int i, l;

  int Px, Py;

  for(l = 0; l < Nof_base_functions; l++)
    {
      for(i = 0; i < 4; i++)
        {
          index_to_base_function(l, &Px, &Py);
          mpi_printf("base function: %d (Px:%d, Py:%d), w0:%f, w1:%f, w2:%f, w3:%f\n", l, Px, Py, 0, 0, 0, 0, SphP[cell].Weights[l][0], SphP[cell].Weights[l][1], SphP[cell].Weights[l][2],
                     SphP[cell].Weights[l][3]);
        }
    }
}

void print_imported_weights()
{
  int i, l;

  for(i = 0; i < Mesh_nimport; i++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          printf("weight: %f\n", Aexch[i][l][0]);
          printf("weight: %f\n", Aexch[i][l][1]);
          printf("weight: %f\n", Aexch[i][l][2]);
          printf("weight: %f\n", Aexch[i][l][3]);
          printf("\n");
        }
    }
}

void assert_imported_weights_positivity(const char *error)
{
  int i, l;

  for(i = 0; i < Mesh_nimport; i++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          if(Aexch[i][l][0] < 0 || Aexch[i][l][1] < 0 || Aexch[i][l][2] < 0 || Aexch[i][l][3] < 0)
            {
              printf("weight: %f\n", Aexch[i][l][0]);
              printf("weight: %f\n", Aexch[i][l][1]);
              printf("weight: %f\n", Aexch[i][l][2]);
              printf("weight: %f\n", Aexch[i][l][3]);
              printf("\n");

              terminate(error);
            }
        }
    }
}



void print_timebins()
{
  int n;

  for(n = 0; n < NumGas; n++)
    {
      printf("bin:%d\n", P[n].TimeBinHydro);
    }
}

void print_weights_array(CBV(a))
{
  int i, l, k;

  for(i = 0; i < NumGas; i++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 4; k++)
            {
              printf("task:%d, a[%d][%d][%d]: %f\n", ThisTask, i, l, k, a[i][l][k]);
            }
        }
    }
}

void print_weights_a(int cell, CBV(a))
{
  int l, k;

  for(l = 0; l < Nof_base_functions; l++)
    {
      for(k = 0; k < 4; k++)
        {
          printf("task:%d, a[%d][%d][%d]: %f\n", ThisTask, cell, l, k, a[cell][l][k]);
        }
    }
}

void assert_weights_positivity(CBV(a), const char *error)
{
  int cell, l, k;

  for(cell = 0; cell < NumGas; cell++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 4; k++)
            {
              if(a[cell][l][k] < 0)
                {
                  terminate(error);
                }
            }
        }
    }
}

void assert_sphp_weights_positivity(const char *error)
{
  int cell, l, k;

  for(cell = 0; cell < NumGas; cell++)
    {
      for(l = 0; l < Nof_base_functions; l++)
        {
          for(k = 0; k < 4; k++)
            {
              if(SphP[cell].Weights[l][k] < 0)
                {
                  terminate(error);
                }
            }
        }
    }
}

void assert_state_positivity(struct state *s, const char *error)
{
  if(s->Energy < 0 || s->press < 0 || s->velx < 0 || s->vely < 0 || s->velz < 0)
    {
      terminate(error);
    }

}
#endif /* DG_DEBUG */
