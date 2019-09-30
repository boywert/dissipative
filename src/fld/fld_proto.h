/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fld/fld_proto.h
 * \date        09/2014
 * \author      Andreas Bauer (andreas.bauer@h-its.org)
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

int fld(void);

void fld_radtransfer(double *XVec, double dt);

double fld_vector_multiply(double *a, double *b);
double fld_vector_sum(double *a);
void fld_matrix_multiply(double *in, double *out, double *sum);


void fld_update_kappa();
void fld_setup_matrix(int cone, double dt);
void fld_compute_flux_limiter();
void fld_source(double dt);
void fld_add_emission(double dt);


void fld_update_gas(double dt, double gamma_new[FLD_NCONES][NumGas], MyFloat *u_new);
void fld_explicit_term(double dt);
double fld_update_state(double dt, double *gamma_new, double *u_new, double *u_old, double tau, double *change);
void fld_compute_coeff(int cone, double dt, double* gamma_new, double* u_new, double* u_old, double tau);
void fld_update_kappa_P(MyFloat *u_new);
void fld_diffusion_operator(double *gamma_new, double *out);
double fld_sigma(double *gamma, double *u_new, double dt);

#ifdef FLD_CONES
void fld_update_n_gamma();
void fld_limit_cone(double out[], double grad[], int cone);
void fld_get_vectors();
#endif
