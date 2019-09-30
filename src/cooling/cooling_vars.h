/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/cooling/cooling_vars.h
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

#define NCOOLTAB       2000
#define SMALLNUM       1.0e-60
#define COOLLIM        0.1
#define HEATLIM        20.0
#define eV_to_K        11606.0
#define eV_to_erg      1.60184e-12
#define MAX_TABLESIZE  250      /* Max # of lines in TREECOOL */

#ifdef EXPLICIT_COOLING
#define COOLING_TOLERANCE  0.05 /* tolerance on the new value of u for which to use explicit cooling */
#endif

#ifdef GFM_AGN_RADIATION
#define GFM_MIN_AGN_BOL_INTENSITY      1e-35
#define GFM_MIN_IONIZATION_PARAMETER   1e-3
#define GFM_MAX_IONIZATION_PARAMETER   1e8
#endif


#ifdef RADCOOL
#define PHIOS_MIN         1e4
#define PHIOS_MAX         1e12
#define PHINS_MIN         1e-7
#define PHINS_MAX         1e3
#define PHIT_MIN          1e14
#define PHIT_MAX          3.162277e23
#define PHIOS_MINVAL      1e5
#define PHINS_MINVAL      1e-4
#define PHIT_MINVAL       3.162277e15
#define PHIOS_MINVAL_LZ   1e6
#define PHINS_MINVAL_LZ   1e-5
#define EPS               1e-2  /*hack to get the values within range */
#define SOLARMASS_in_g    1.989e33
#define KPC_in_cm         3.0856776e21
#define GYR_to_YR         1e9

#ifdef RADCOOL_HOTHALO
#define HOTHALO_ON_zFACTOR 3.0  /*Redshift at which the hothalo emission turns on || z=3 */
#define FACTOR_Mp2_cm5     5055510.5    /*conversion from solarmass^2/kpc^5 to protonmass^2/cm^5 */
#define INVERSE_mu2        2.89 /*1/meanweight^2 meanweight ~ 0.58 fully ionized */
#endif


#endif

/* data for gas state */
typedef struct
{
  double ne, necgs, nHcgs;
  double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
  double gJH0ne, gJHe0ne, gJHepne;
  double nH0, nHp, nHep, nHe0, nHepp;
  double XH, yhelium;
  double mhboltz;
  double ethmin;                /* minimum internal energy for neutral gas */
  double mu;
#ifdef GFM_COOLING_METAL
  double log_MetallicityInSolar;
  double log_HydrogenNumberDensity;
#endif
#ifdef GFM_AGN_RADIATION
  int FlagAGNBackground;
  double LogAGNBolIntensity, AGNBolIntensity;
#endif
#ifdef RADCOOL
  int FlagRadcool;
  double LogPhios, LogPhins;
  double Phios, Phins;
#ifdef RADCOOL_HOTHALO
  double LogPhiT6, LogPhiT7, LogPhiT8, PhiT6, PhiT7, PhiT8;
#endif
#endif
#ifdef GFM_UVB_CORRECTIONS
  double UV_HeII_Factor;
  double RhoGasMean;
#endif
#ifdef GFM_DUST_COOLING
  double dust_to_gas_ratio;
  char dust_cool;
#endif
#ifdef GFM_LAMBDA
  int partindex;
#endif
} GasState;


/* tabulated rates */
typedef struct
{
  double BetaH0, BetaHep, Betaff;
  double AlphaHp, AlphaHep, Alphad, AlphaHepp;
  double GammaeH0, GammaeHe0, GammaeHep;
} RateTable;

/* photo-ionization/heating rate table */
typedef struct
{
  float variable;               /* logz for UVB */
  float gH0, gHe, gHep;         /* photo-ionization rates */
  float eH0, eHe, eHep;         /* photo-heating rates */
} PhotoTable;

/* current interpolated photo-ionization/heating rates */
typedef struct
{
  char J_UV;
  double gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
} PhotoCurrent;

/* cooling data */
typedef struct
{
  double u_old_input, rho_input, dt_input, ne_guess_input;
} DoCoolData;

#ifdef UVB_SELF_SHIELDING
typedef struct
{
  float redshift;
  float alpha1, alpha2, beta, xi, n0, f;
} SelfShieldingTable;
#endif
