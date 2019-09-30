/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_test_utils.c
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include "tgchem_test.h"


double dabs(double a)
{
  if(a < 0)
    return -a;
  else
    return a;
}


double dmax(double a, double b)
{
  if(a > b)
    return a;
  else
    return b;
}


double dmin(double a, double b)
{
  if(a < b)
    return a;
  else
    return b;
}


int imax(int a, int b)
{
  if(a > b)
    return a;
  else
    return b;
}


int imin(int a, int b)
{
  if(a < b)
    return a;
  else
    return b;
}


double second(void)
{
  return MPI_Wtime();
}


double measure_time(void)
{
  double t, dt;

  t = second();

  dt = t - WallClockTime;

  WallClockTime = t;

  return dt;
}


double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)
    dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;

  return dt;
}
