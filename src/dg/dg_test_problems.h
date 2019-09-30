/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dg/dg_test_problems.h
 * \date        10/2014
 * \author		Kevin Schaal
 * \brief	    functions for dg test problems	
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef DG_TEST_PROBLEMS_H
#define DG_TEST_PROBLEMS_H

#ifndef DG
void arepo_set_initial_conditions();
#endif

/*
 * dg_test_problems.c
 */

double ic_pressure(double x, double y, double z);
double ic_density(double x, double y, double z);
double ic_velocity_x(double x, double y, double z);
double ic_velocity_y(double x, double y, double z);
double ic_velocity_z(double x, double y, double z);

/*
 * cell_projection.c
 */

#ifndef DG
double integrate(int cell, double (*func)(int cell, double x, double y, double z));
double calc_L1_norm(int cell);
#endif

#endif
