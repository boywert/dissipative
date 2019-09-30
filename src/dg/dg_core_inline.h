/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_core_inline.h
 * \date        06/2015
 * \author      Kevin Schaal
 * \brief       Inline functions
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

/*!
 *  Inline functions
 */

inline static double dg_sound_speed(int cell);
inline static void primvars_to_w(double *pv, double *w);
inline static int same_w(double *w1, double *w2);
inline static double ideal_gas_pressure(double density, double utherm);
inline static double ideal_gas_density(double pressure, double utherm);
inline static double ideal_gas_utherm(double density, double pressure);
inline static void w_to_primvars(double *w, double *pv);
inline static void lab_coords_to_cell_coords(int cell, double x_lab, double y_lab, double z_lab, double *x_cell, double *y_cell, double *z_cell);
inline static void cell_coords_to_lab_coords(int cell, double x_cell, double y_cell, double z_cell, double *x_lab, double *y_lab, double *z_lab);
inline static void sub_vector_vector(double v1[5], double v2[5], double v_out[5]);
inline static void copy_vector_vector(double *v_origin, double *v_destin, int dims);

/*!
 * Calculate the sound speed of a cell
 */
inline static double dg_sound_speed(int cell)
{
  return sqrt(GAMMA * SphP[cell].Pressure / SphP[cell].Density);
}

/*!
 * Calculate the w vector from the primitive variables
 */
inline static void primvars_to_w(double *pv, double *w)
{
  w[0] = pv[RHO];
  w[1] = pv[RHO] * pv[VX];
  w[2] = pv[RHO] * pv[VY];
  w[3] = pv[RHO] * pv[VZ];
  w[4] = pv[PRES] / (GAMMA - 1) + 0.5 * pv[RHO] * (pv[VX] * pv[VX] + pv[VY] * pv[VY]+ pv[VZ] * pv[VZ]);
}

/*!
 * Check whether two state are the same
 */
inline static int same_w(double *w1, double *w2)
{
  return (w1[0] == w2[0] && w1[1] == w2[1] && w1[2] == w2[2] && w1[3] == w2[3] && w1[4] == w2[4]);
}

/*!
 * ideal gas equation
 */

inline static double ideal_gas_pressure(double density, double utherm)
{
  return density * utherm * (GAMMA - 1);
}

inline static double ideal_gas_density(double pressure, double utherm)
{
  return pressure / (utherm * (GAMMA - 1));
}

inline static double ideal_gas_utherm(double density, double pressure)
{
  return pressure / (density * (GAMMA - 1));
}

/*!
 * Calculate the primitive variables from the w vector
 */
inline static void w_to_primvars(double *w, double *pv)
{
  pv[RHO] = w[W_RHO];
  pv[VX] = w[W_PX] / w[W_RHO];
  pv[VY] = w[W_PY] / w[W_RHO];
  pv[VZ] = w[W_PZ] / w[W_RHO];

#ifdef OLD_VERSION
  pv[PRES] = (w[W_E] - 0.5 * (w[W_PX]*w[W_PX]/w[W_RHO]+w[W_PY]*w[W_PY]/w[W_RHO]+w[W_PZ]*w[W_PZ]/w[W_RHO])) * (GAMMA - 1);
#else
  pv[PRES] = (w[W_E] - 0.5 * (w[W_PX]*w[W_PX]+w[W_PY]*w[W_PY]+w[W_PZ]*w[W_PZ])/w[W_RHO]) * (GAMMA - 1);
#endif
}


/*!
 * Transformation between coordinates in the lab system and coordinates in the cell system ([-1,1])
 */


inline static void lab_coords_to_cell_coords(int cell, double x_lab, double y_lab, double z_lab, double *x_cell, double *y_cell, double *z_cell)
{
  double dl = amr_length[Mesh.DP[cell].level];

  *x_cell = (x_lab - SphP[cell].Center[0]) * 2. / dl;
  *y_cell = (y_lab - SphP[cell].Center[1]) * 2. / dl;
  *z_cell = (z_lab - SphP[cell].Center[2]) * 2. / dl;
}


inline static void cell_coords_to_lab_coords(int cell, double x_cell, double y_cell, double z_cell, double *x_lab, double *y_lab, double *z_lab)
{
  double dl = amr_length[Mesh.DP[cell].level];

  *x_lab = x_cell * dl / 2. + SphP[cell].Center[0];
  *y_lab = y_cell * dl / 2. + SphP[cell].Center[1];
  *z_lab = z_cell * dl / 2. + SphP[cell].Center[2];
}

/*!
 * Subtract two vectors
 */

inline static void sub_vector_vector(double v1[5], double v2[5], double v_out[5])
{
  v_out[0] = v1[0] - v2[0];
  v_out[1] = v1[1] - v2[1];
  v_out[2] = v1[2] - v2[2];
  v_out[3] = v1[3] - v2[3];
  v_out[4] = v1[4] - v2[4];
}

/*!
 * Copy a vector to a vector
 */

inline static void copy_vector_vector(double *v_origin, double *v_destin, int dims)
{
  int i;

  for(i = 0; i < dims; i++)
    {
      v_destin[i] = v_origin[i];
    }
}
