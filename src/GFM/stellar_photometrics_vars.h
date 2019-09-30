/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_photometrics_vars.h
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

#ifndef STELLAR_PHOTOMETRICS_VARS_H
#define STELLAR_PHOTOMETRICS_VARS_H

#define GFM_STELLAR_PHOTOMETRICS_BANDS 8        /* number of bands */
#define GFM_STELLAR_PHOTOMETRICS_DIRECTIONS 100 /* number of directions for surface brightness calculation */
#define GFM_STELLAR_PHOTOMETRICS_RADII 100      /* number of radii (rings) for surface brightness calculation (must be >=3) */
#define GFM_STELLAR_PHOTOMETRICS_K_LIMIT 20.7   /* limiting surface brightness determining 'detectable radius' */

typedef struct
{
  MyFloat Magnitude_U, Magnitude_B, Magnitude_V, Magnitude_K;
  MyFloat Magnitude_g, Magnitude_r, Magnitude_i, Magnitude_z;
} stellar_photometrics;

extern struct StellarLuminosityTable
{
  unsigned int N_LogMetallicity;        /* metallicity NOT in solar */
  unsigned int N_LogAgeInGyr;

  MyFloat *LogMetallicity_bins;
  MyFloat *LogAgeInGyr_bins;
  MyFloat **Magnitude_U, **Magnitude_B, **Magnitude_V, **Magnitude_K;
  MyFloat **Magnitude_g, **Magnitude_r, **Magnitude_i, **Magnitude_z;
} stellarLuminosityTable;

extern float StellarPhotometricsRandomAngles[GFM_STELLAR_PHOTOMETRICS_DIRECTIONS][2];

extern int dummyvarstellarphotometrics;

#endif
