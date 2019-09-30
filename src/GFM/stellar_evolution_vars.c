/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_evolution_vars.c
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

#include "../allvars.h"

#ifdef GFM_STELLAR_EVOLUTION

char ElementNames[GFM_N_CHEM_ELEMENTS][GFM_EL_NAME_LENGTH];

int element_index_Hydrogen;
int element_index_Helium;
int element_index_Carbon;
int element_index_Nitrogen;
int element_index_Oxygen;
int element_index_Neon;
int element_index_Magnesium;
int element_index_Silicon;
int element_index_Iron;


double imf_by_number[GFM_N_MASS_BINS], imf_mass_bin[GFM_N_MASS_BINS], imf_mass_bin_log10[GFM_N_MASS_BINS];
double yield_mass_bin[GFM_N_MASS_BINS];
double imf_dlog10_Msun;

struct YieldTypeIa yieldsSNIa;
struct YieldTypeII_and_AGB yieldsSNII, yieldsAGB;
struct Lifetime_Table Lifetimes;

struct star_particle *StarParticle;

/* this variable initialization is just added to fix a strange linking problem on Mac OSX - without it, the file cooling_metal_vars.c is somehow not linked in */
int Nstar = 0;
#endif

#ifdef GFM_WINDS_STRIPPING
int Nwinds_to_strip = 0;
#endif
