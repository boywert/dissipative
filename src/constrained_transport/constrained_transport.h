/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/constrained_transport/constrained_transport.h
 * \date        06/2015
 * \author      Philip Mocz
 * \brief       includes functions for constrained transport algorithm
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef CONSTRAINED_TRANSPORT_H
#define CONSTRAINED_TRANSPORT_H
 
#ifdef MHD
#ifdef MHD_CT

/* Constrained Transport Main Methods */
void set_A_ICs(int ic_flag);
void correct_ctr_b(void);
void update_A(int i, double dt);
void do_mhd_ct_source_terms();
void get_A_fluxes(int k, int q, int qother, double face_normal_vel, double face_velx, struct state *st_L, struct state *st_R, struct fluxes *flux);

/* Constrained Transport Helpers */
void get_init_A(double x, double y, double z, double *Ax, double *Ay, double *Az, int ic_flag);
void get_face_area(point * dp0, point * dp1, point * dp2, double *Ax, double *Ay, double *Az);
double get_phi(point * dp0, point * dp1, point * dp2, point * dp_ref);
void solve3x3(double A[3][3], double b[3], double *x, double *y, double *z);
void A_periodic_shift_correction(point *dp0, point *dp1, double shift[3]);
void switch_to_CM(point** p0, point* p1);
void time_extrapolate_A_if_needed(double A[3], point *dp, point *dp_ref);
int not_both_active(point *dp_ref, point *dp0);

#endif
#endif

#endif /* CONSTRAINED_TRANSPORT_H */

