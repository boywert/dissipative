/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_test_proto.h
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

void tgchem_test_pars();
void tgchem_test_init();
double dabs(double a);
double dmax(double a, double b);
double dmin(double a, double b);
int imax(int a, int b);
int imin(int a, int b);
double second(void);
double measure_time(void);
double timediff(double t0, double t1);
