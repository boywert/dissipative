/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/winds_proto.h
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

#if defined(GFM_WINDS) || defined(GFM_WINDS_LOCAL)
void do_winds(void);
int data_index_compare_wind(const void *a, const void *b);
void setup_windparticles(void);
void find_wind_cells(int Nwind);
void recouple_wind_particles(int Nwind, int *ret_recoupled, double *ret_mass_recoupled);
void cool_and_check_recouple(void);
#endif
#ifdef GFM_WINDS
void init_winds(void);
#endif
#ifdef GFM_WINDS_LOCAL
void create_winds_local(void);
#endif
#ifdef GFM_WINDS_VARIABLE
void gfm_calc_variable_wind_parameters(double stellar_mass, double halo_mass, double vel_disp, double metallicity, wind_parameter * wp);
void init_variable_winds(void);
#endif
