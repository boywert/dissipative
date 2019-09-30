/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/cell_projection.c
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

#ifndef DG

/*!
 *  Arepo: Calculate the Integral(func(x,y,z))dV for the cell.
 */
double integrate(int cell, double (*func) (int cell, double x, double y, double z))
{
  int i, j;
  int k = 0;
  double x, y, z;

  int nof_points = 10;

  double x_cell = SphP[cell].Center[0];
  double y_cell = SphP[cell].Center[1];
  double z_cell = SphP[cell].Center[2];

  double l_cell = amr_length[Mesh.DP[cell].level];

  double dd = l_cell / nof_points;

  double px, py, pz;

  px = x_cell + 0.5 * (dd - l_cell);
  py = y_cell + 0.5 * (dd - l_cell);
  pz = z_cell + 0.5 * (dd - l_cell);


  double result = 0;
  int counter = 0;


  for(i = 0; i < nof_points; i++)
    {
      for(j = 0; j < nof_points; j++)
        {
#ifndef TWODIMS
          for(k = 0; k < nof_points; k++)
            {
#endif
              x = px + i * dd;
              y = py + j * dd;
              z = pz + k * dd;

              result += func(cell, x, y, z);

              counter++;
#ifndef TWODIMS
            }
#endif
        }
    }

  result *= SphP[cell].Volume / counter;

  return result;
}

/*!
 *  Arepo: Calculate the L1 Norm Integral(|rho0-rho|)dV for the cell.
 */

static double L1_density_error(int cell, double x, double y, double z)
{
  return fabs(SphP[cell].Density - ic_density(x, y, z));
}

double calc_L1_norm(int cell)
{
  double result = integrate(cell, &L1_density_error);

  return result;
}

#endif
