/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/modgrav/modgrav_pm.h
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

struct pm_amr_data
{
  int task;
  MyFloat center[3];
  float len;
  float eff_mass;
};

void modgrav_add_amr_nodes_to_periodic_grid(fft_real *rhogrid, MyFloat to_slab_fac, int fftsize, fft_plan myplan);
void modgrav_add_amr_nodes_to_nonperiodic_grid(int grnr, fft_real *rhogrid, double to_slab_fac, fft_plan myplan,  int fftsize);

