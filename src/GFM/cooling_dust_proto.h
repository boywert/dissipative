/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/cooling_dust_proto.h
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

double get_CoolingDustRate(double logT, double dust_to_gas_ratio, double mu, char dust_cool);
void get_CoolingDustState(int i, double *dust_to_gas_ratio, char *dust_cool);
