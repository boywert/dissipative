/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_set_get.c
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief		Setter and getter methods for the global arrays
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

#if ((defined DG) && ((defined SAVE_ACCESS) || defined(DG_DEBUG)))

/*
 * setter methods which check whether the indices are valid
 */

void set_lobatto_points_1d(int q, double value)
{
  assert(q < Nof_lobatto_points_1d);

  Lobatto_points_1d[q] = value;

}

void set_lobatto_points_weights_1d(int q, double value)
{
  assert(q < Nof_lobatto_points_1d);

  Lobatto_points_weights_1d[q] = value;
}

void set_quad_points_1d(int q, double value)
{
  assert(q < Nof_quad_points_1d);

  Quad_points_1d[q] = value;

}

void set_quad_points_weights_1d(int q, double value)
{
  assert(q < Nof_quad_points_1d);

  Quad_points_weights_1d[q] = value;
}

void set_union_quad_points_x(int q, double value)
{
  assert(q < Nof_union_quad_points);

  Union_quad_points_x[q] = value;
}

void set_union_quad_points_y(int q, double value)
{
  assert(q < Nof_union_quad_points);

  Union_quad_points_y[q] = value;
}

void set_union_quad_points_z(int q, double value)
{
  assert(q < Nof_union_quad_points);

  Union_quad_points_z[q] = value;
}

void set_inner_quad_points_x(int q, double value)
{
  assert(q < Nof_inner_quad_points);

  Inner_quad_points_x[q] = value;
}

void set_inner_quad_points_y(int q, double value)
{
  assert(q < Nof_inner_quad_points);

  Inner_quad_points_y[q] = value;
}

void set_inner_quad_points_z(int q, double value)
{
  assert(q < Nof_inner_quad_points);

  Inner_quad_points_z[q] = value;
}

void set_inner_quad_points_weights(int q, double value)
{
  assert(q < Nof_inner_quad_points);

  Inner_quad_points_weights[q] = value;
}

void set_outer_quad_points_x(int e, int q, double value)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  Outer_quad_points_x[e][q] = value;
}

void set_outer_quad_points_y(int e, int q, double value)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  Outer_quad_points_y[e][q] = value;
}

void set_outer_quad_points_z(int e, int q, double value)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  Outer_quad_points_z[e][q] = value;
}

void set_outer_quad_points_weights(int e, int q, double value)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  Outer_quad_points_weights[e][q] = value;
}

void set_union_base_values(int q, int l, double value)
{
  assert(q < Nof_union_quad_points);
  assert(l < Nof_base_functions);

  Union_base_values[q][l] = value;
}

void set_inner_base_values(int q, int l, double value)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  Inner_base_values[q][l] = value;
}

void set_inner_base_dx_values(int q, int l, double value)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  Inner_base_dx_values[q][l] = value;
}

void set_inner_base_dy_values(int q, int l, double value)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  Inner_base_dy_values[q][l] = value;
}

void set_inner_base_dz_values(int q, int l, double value)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  Inner_base_dz_values[q][l] = value;
}

void set_outer_base_values(int e, int q, int l, double value)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);
  assert(l < Nof_base_functions);

  Outer_base_values[e][q][l] = value;
}

void set_outer_fine_base_values(int e, int q, int l, double value)
{
  assert(e < NOF_INTERFACES);   //interface parts
  assert(q < NOF_OUTER_FINE_QUAD_POINTS);
  assert(l < Nof_base_functions);

  Outer_fine_base_values[e][q][l] = value;
}

void set_p_a(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_A[j][l] = value;
}

void set_p_b(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_B[j][l] = value;
}

void set_p_c(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_C[j][l] = value;
}

void set_p_d(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_D[j][l] = value;
}

#ifndef TWODIMS
void set_p_e(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_E[j][l] = value;
}

void set_p_f(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_F[j][l] = value;
}

void set_p_g(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_G[j][l] = value;
}

void set_p_h(int j, int l, double value)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  P_H[j][l] = value;
}
#endif

/*
 * getter methods which check whether the indices are valid
 */

double get_lobatto_points_1d(int q)
{
  assert(q < Nof_lobatto_points_1d);

  return Lobatto_points_1d[q];

}

double get_lobatto_points_weights_1d(int q)
{
  assert(q < Nof_lobatto_points_1d);

  return Lobatto_points_weights_1d[q];
}

double get_quad_points_1d(int q)
{
  assert(q < Nof_quad_points_1d);

  return Quad_points_1d[q];

}

double get_quad_points_weights_1d(int q)
{
  assert(q < Nof_quad_points_1d);

  return Quad_points_weights_1d[q];
}

double get_union_quad_points_x(int q)
{
  assert(q < Nof_union_quad_points);

  return Union_quad_points_x[q];
}

double get_union_quad_points_y(int q)
{
  assert(q < Nof_union_quad_points);

  return Union_quad_points_y[q];
}

double get_union_quad_points_z(int q)
{
  assert(q < Nof_union_quad_points);

  return Union_quad_points_z[q];
}

double get_inner_quad_points_x(int q)
{
  assert(q < Nof_inner_quad_points);

  return Inner_quad_points_x[q];
}

double get_inner_quad_points_y(int q)
{
  assert(q < Nof_inner_quad_points);

  return Inner_quad_points_y[q];
}

double get_inner_quad_points_z(int q)
{
  assert(q < Nof_inner_quad_points);

  return Inner_quad_points_z[q];
}

double get_inner_quad_points_weights(int q)
{
  assert(q < Nof_inner_quad_points);

  return Inner_quad_points_weights[q];
}

double get_outer_quad_points_x(int e, int q)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  return Outer_quad_points_x[e][q];
}

double get_outer_quad_points_y(int e, int q)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  return Outer_quad_points_y[e][q];
}

double get_outer_quad_points_z(int e, int q)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  return Outer_quad_points_z[e][q];
}

double get_outer_quad_points_weights(int e, int q)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);

  return Outer_quad_points_weights[e][q];
}

double get_union_base_values(int q, int l)
{
  assert(q < Nof_union_quad_points);
  assert(l < Nof_base_functions);

  return Union_base_values[q][l];
}

double get_inner_base_values(int q, int l)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  return Inner_base_values[q][l];
}

double get_inner_base_dx_values(int q, int l)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  return Inner_base_dx_values[q][l];
}

double get_inner_base_dy_values(int q, int l)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  return Inner_base_dy_values[q][l];
}

double get_inner_base_dz_values(int q, int l)
{
  assert(q < Nof_inner_quad_points);
  assert(l < Nof_base_functions);

  return Inner_base_dz_values[q][l];
}

double get_outer_base_values(int e, int q, int l)
{
  assert(e < Nof_interfaces);   //interfaces
  assert(q < Nof_outer_quad_points);
  assert(l < Nof_base_functions);

  return Outer_base_values[e][q][l];
}

double get_outer_fine_base_values(int e, int q, int l)
{
  assert(e < NOF_INTERFACES);   //interfaces
  assert(q < NOF_OUTER_FINE_QUAD_POINTS);
  assert(l < Nof_base_functions);

  return Outer_fine_base_values[e][q][l];
}

double get_p_a(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_A[j][l];
}

double get_p_b(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_B[j][l];
}

double get_p_c(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_C[j][l];
}

double get_p_d(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_D[j][l];
}

#ifndef TWODIMS
double get_p_e(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_E[j][l];
}

double get_p_f(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_F[j][l];
}

double get_p_g(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_G[j][l];
}

double get_p_h(int j, int l)
{
  assert(j < Nof_base_functions);
  assert(l < Nof_base_functions);

  return P_H[j][l];
}
#endif

#endif /* DG && SAVE_ACCESS */
