/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/OTVET/otvet_proto.h
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

int do_otvet(void);
void otvet_eddington(void);
int eddington_treeevaluate(int target, int mode, int *nexport, int *nsend_local);
void otvet_star_lum(void);
int star_lum_evaluate(int target, int mode, int *nexport, int *nsend_local);
void otvet_compute_sphdensity(void);
void otvet_radtransfer(void);
void ot_get_sigma(void);
void ot_get_lum_stars(void);
void otvet_set_simple_inits(void);
double otvet_vector_multiply(double *a, double *b);
double otvet_vector_sum(double *a);
void otvet_matrix_multiply(double *in, double *out, double *sum);
void otvet_update_chemistry(void);
double otvet_DoHeating(int i, double dt_internal);
double otvet_DoCooling(int i, double dt_internal);
double otvet_get_heating_rate(int i);
double otvet_get_cooling_rate(int i, double utherm);
double otvet_GetCoolingTime(int i, double u_old, double rho, double *ne_guess);
void otvet_write_stats(void);
void otvet_update_chemistry(void);
void otvet_stardensity(void);
