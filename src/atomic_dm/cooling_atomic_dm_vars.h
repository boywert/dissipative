/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/atomic_dm/cooling_atomic_dm_vars.h
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

#define SMALLNUM       1.0e-60
#define eV_to_K        11606.0
#define eV_to_erg      1.60184e-12
#define MAX_TABLESIZE  250      /* Max # of lines in TREECOOL */
#define keV_to_gramm   1.782662e-30

/* data for gas state */
typedef struct
{
  double ne, nHcgs;
  double nH0, nHp;
  double ethmin;                /* minimum internal energy for neutral gas */
} ADM_GasState;

/* cooling data */
typedef struct
{
  double u_old_input, rho_input, dt_input, ne_guess_input;
} ADM_DoCoolData;
