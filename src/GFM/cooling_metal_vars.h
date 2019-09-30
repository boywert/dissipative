/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/cooling_metal_vars.h
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

#ifndef COOLING_METAL_VARS_H
#define COOLING_METAL_VARS_H


/* solar metal mass fraction */
#ifndef RADCOOL_HOTHALO_METAL_BOOST
#define GFM_SOLAR_METALLICITY 0.0127
#endif

#define GFM_SOLAR_HE_ABUNDANCE 0.1
#ifdef RADCOOL_HOTHALO
#define HOTHALO_ON_zFACTOR 3.0
#endif
extern struct CoolingMetalTable
{
  unsigned int N_Redshift;
#ifdef GFM_AGN_RADIATION
  unsigned int N_BolFlux;
#endif
#ifdef RADCOOL
  unsigned int N_Phios;
  unsigned int N_Phins;
#ifdef RADCOOL_HOTHALO
  unsigned int N_Redshift_lz;
  unsigned int N_Phios_lz;
  unsigned int N_Phins_lz;
  unsigned int N_PhiT6;
  unsigned int N_PhiT7;
  unsigned int N_PhiT8;
#endif
#endif
#ifndef RADCOOL
  unsigned int N_MetallicityInSolar;
#endif
  unsigned int N_HydrogenNumberDensity;
  unsigned int N_Temperature;

  MyFloat *Redshift_bins;
#ifdef GFM_AGN_RADIATION
  MyFloat *BolFlux_bins;
#endif
#ifdef RADCOOL
  MyFloat *Phios_bins;
  MyFloat *Phins_bins;
#ifdef RADCOOL_HOTHALO
  MyFloat *Redshift_bins_lz;
  MyFloat *Phios_bins_lz;
  MyFloat *Phins_bins_lz;
  MyFloat *PhiT6_bins;
  MyFloat *PhiT7_bins;
  MyFloat *PhiT8_bins;
#endif
#endif
#ifndef RADCOOL
  MyFloat *MetallicityInSolar_bins;
#endif
  MyFloat *HydrogenNumberDensity_bins;
  MyFloat *Temperature_bins;


#if defined(GFM_AGN_RADIATION)
  MyFloat *****NetCoolingRate;
#else
#ifdef RADCOOL
  float *****NetCoolingRate;
  float *****NetHeatingRate;
#ifdef RADCOOL_HOTHALO
  float ********NetCoolingRateLowZ;
  float ********NetHeatingRateLowZ;
#endif
#else
  MyFloat ****NetCoolingRate;
#endif
#endif
}
coolingMetalTable;

extern int dummyvar;

#endif
