/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_vars.h
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

#ifndef DUST_VARS_H
#define DUST_VARS_H

extern int Ndust;

extern struct dust_particle
{
  int index;
  int active_idx;
  MyFloat NumNgb;
  MyFloat NormSph;
  MyFloat TotNgbMass;
  MyFloat Dhsmlrho;
  MyFloat LocalGasAccel[3];
  MyFloat LocalGradP[3];
#ifdef DL_GRAIN_BINS
  MyFloat LocalGasDensityH;
  MyFloat LocalGasZ;
  MyFloat LocalGasTemp;
  MyFloat NewNumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  MyFloat NewBinSlopes[DL_GRAIN_BINS];
#endif
  MyFloat DeltaMassExpected;
  MyFloat DeltaMassActual;
  MyFloat DeltaMetalMasses[GFM_N_CHEM_ELEMENTS];
  MyFloat DeltaMomentum[3];
#ifdef DL_SNE_DESTRUCTION
  MyFloat LocalSNPrefactor;
#endif
#endif
} *DustParticle;

#ifdef DL_GRAIN_BINS
extern struct grain_size_distribution
{
  double InternalDensity; /* Grain density in internal mass / um^3 */
  double *Edges;
  double *Midpoints;
  double *Widths;
  double *AvgMasses;
  double *EdgeMasses;
#if defined(DL_SHATTERING) || defined(DL_COAGULATION)
  double *VelocitiesCNM; /* CNM grain velocities in um / s */
  double *VelocitiesWIM; /* WIM grain velocities in um / s */
  double I_2_kj_Prefac[DL_GRAIN_BINS][DL_GRAIN_BINS]; /* prefactors for shattering and coagulation integrals */
#ifdef DL_SHATTERING
  double VelShat; /* Shattering threshold velocity in um / s */
#endif
#ifdef DL_COAGULATION
  double VelCoag[DL_GRAIN_BINS][DL_GRAIN_BINS]; /* Coagulation threshold velocity in um / s */
#endif
#endif
#ifdef DL_SNE_DESTRUCTION
  /* Consult equation 11 of Asano+ (2013) for the full definition of this xi
   * fraction.  In short, XiFrac[i][j] describes the differential number
   * fraction of grains shifted from bin j to bin i as a result of one
   * supernova shock. */
  double XiFrac[DL_GRAIN_BINS][DL_GRAIN_BINS];
#endif
#ifdef DL_PRODUCTION
  /* Initial grain size distributions on a relative scale. */
  double AGB_NumGrains[DL_GRAIN_BINS];
  double SNII_NumGrains[DL_GRAIN_BINS];
#ifdef DL_GRAIN_BINS_PIECEWISE_LINEAR
  double AGB_BinSlopes[DL_GRAIN_BINS];
  double SNII_BinSlopes[DL_GRAIN_BINS];
#endif
  /* Condensation efficiencies for SNII */
  double SNII_CondEff[GFM_N_CHEM_ELEMENTS];
  double ***AGB_CondEff;
  double ***AGB_CondEff_spline;
  int AGB_CondEff_NZ;
  int AGB_CondEff_NM;
  double *AGB_CondEff_Mass;
  double *AGB_CondEff_Metallicity;
#endif
} GSD;
#endif

enum gsd_collision_type
{
  GSD_SHATTERING,
  GSD_COAGULATION
};

#ifdef DL_GRAIN_BINS
#if defined(DL_SNE_DESTRUCTION) || defined(DL_SHATTERING) || defined(DL_COAGULATION)
extern struct shatter_data_p
{
  MyFloat DustHsml;
  MyFloat DustNumNgb;
  MyFloat DustDensity;
} *PShatter;
#endif
#endif

#ifdef DL_PRODUCTION
enum gsd_dnda_type
{
  GSD_DNDA_AGB,
  GSD_DNDA_SNII,
  GSD_DNDA_SNIA
};
#endif

#endif
