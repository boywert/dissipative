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


#ifdef MRT


#define fTOLERENCE 1.000000001

#ifdef MRT_RIEMANN_HLLE
double lambda1[101][101] ;
double lambda2[101][101] ;
double lambda3[101][101] ;
double lambda4[101][101] ;

void readin_hlle_eingenvalues(void) ;

#endif

void init_RT(void);
void add_source_fluxes(void);
void set_VET_single(int, struct sph_particle_data *) ;
void mrt_update_chemistry(void);

#ifdef MRT_COMOVING
void cell_do_lorentz_boost(int, struct sph_particle_data *, double, double, double) ;
void do_comoving_frame_source_terms(void);
#endif


#ifdef MRT_IR_LTE
#ifdef MRT_IR_GRAIN_KAPPA
void read_grain_kappa_data(void);
extern int IR_N_pts;
extern double *IR_logT, *IR_logkappaP, *IR_logkappaR;
extern gsl_interp_accel *accIR_kappaP;
extern gsl_interp_accel *accIR_kappaR;
extern gsl_spline *splineIR_kappaP;
extern gsl_spline *splineIR_kappaR;
#endif
void set_kappa_times_rho_IR(int, struct sph_particle_data *) ;
double mrt_update_IR_cooling(int, double) ;
#ifdef MRT_IR_LTE_GSL
int mrt_IR_rate_ODEs(double, const double *, double *, void *) ;
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) ;
#endif
#endif



#ifdef MRT_RADIATION_PRESSURE
void do_radiation_pressure_source_terms(void);
#endif

#if defined(MRT_COOLING_HEATING)
double mrt_DoHeating(int, double);
double mrt_DoCooling(int, double);
double mrt_get_cooling_rate(int, double);
double mrt_get_heating_rate(int);
double mrt_GetCoolingTime(int, double, double, double*);
#endif

#ifdef MRT_CHEMISTRY_PS2009
void mrt_update_chemistry_ps2009(void);
void mrt_write_stats(void);
#endif

#ifdef MRT_CHEMISTRY_PS2009
void mrt_IR_chemistry(void);
#endif

#ifdef MRT_IR_ONLY_CHEMISTRY
void mrt_IR_chemistry(void) ;
#endif

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
void mrt_setup(void) ;
void exchange_vector(void) ;
#endif

#ifdef MRT_CHEMISTRY_PS2011
int mrt_rate_ODEs(double, const double *, double *, void *) ;
void mrt_update_chemistry_ps2011(void) ;
void mrt_get_sigma(void) ;
#endif

#ifdef MRT_TIME_EXTRAPOLATION
void calculate_div_F(struct state *dl, struct state *st) ;
#endif

#ifdef MRT_SOURCES
extern int Nsource;

#ifdef MRT_STARS
void start_stellar_sources(void);
void end_stellar_sources(void);
void do_ionizing_stellar_sources(void);
void add_ionizing_stellar_radiation(void);
#endif

#ifdef MRT_BH
void start_blackhole_sources(void);
void end_blackhole_sources(void);
void do_ionizing_blackhole_sources(void);
void add_ionizing_blackhole_radiation(void);
#endif

#endif /* MRT_SOURCES */

#endif /* MRT */
