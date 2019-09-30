/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/anisotropic_RT/RT.h
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
//#include "../mesh.h"
#include "../proto.h"

#ifdef MRT

void mrt_run(void) ;
void exchange_primitive_variables_RT(void) ;
void exchange_primitive_variables_and_gradients_RT(void) ;
void update_primitive_variables_RT(void) ;
void calculate_gradients_RT(void) ;
void compute_interface_fluxes_RT(tessellation * T);
int face_get_state_RT(tessellation * T, int p, int i, struct state *st);
void state_convert_to_local_frame_RT(struct state *st, double *vel_face, double hubble_a, double atime);
void face_do_time_extrapolation_RT(struct state *delta, struct state *st, double atime);
void face_do_spatial_extrapolation_RT(struct state *delta, struct state *st, struct state *st_other);
void face_add_extrapolations_RT(struct state *st_face, struct state *delta_time, struct state *delta_space);
void face_add_extrapolation_RT(struct state *st_face, struct state *delta);
void face_turn_velocities_RT(struct state *st, struct geometry *geom);
void face_turnback_velocities_RT(struct state_face *st_face, struct geometry *geom);
void face_limit_fluxes_RT(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R, struct fluxes *flux, double dt, double *count, double *count_reduced);
void apply_flux_list_RT(void); 
double face_timestep_RT(struct state *state_L, struct state *state_R, double *hubble_a, double *atime);


static inline double interpolate_2d(double (*table)[101], int i, int j, double dx, double dy)
{
  return (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] + dx * (1 - dy) * table[i + 1][j] + dx * dy * table[i + 1][j + 1];
}


/*Required structures*/


extern struct rt_face_data {
struct geometry geom ;
int face_active ;
int face_responsibility ;
double face_timestep ;
double vel_face[3] ;
double vel_face_turned[3] ;
//double Ldx , Ldy, Ldz, Rdx, Rdy, Rdz ;
} *rtfacedata ;


#endif
