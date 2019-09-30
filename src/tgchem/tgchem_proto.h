/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_proto.h
 * \date        01/2013
 * \author      Primordial chemistry and cooling network
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

void tgchem();
double tgchem_get_timestep(int i, double dt);

void tgchem_chem(N_Vector dspecies);

void tgchem_cool(N_Vector dspecies);

double tgchem_h2_line_fit_fesc();
double tgchem_h2_line_sob_fesc();

void tgchem_begrun();
void tgchem_init_cvode();
void tgchem_free_cvode();
void tgchem_init_rates();

void tgchem_photo();

void tgchem_step(double dt);

int tgchem_rates(double time, N_Vector species, N_Vector dspecies, void *user_data);
void tgchem_check_for_nan(int mode, N_Vector species);
void tgchem_compute_vars(int mode, N_Vector species);
void tgchem_compute_rates();
void tgchem_compute_aux();
int tgchem_check_species();
void tgchem_debug_rates(int mode, N_Vector species, N_Vector dspecies);
