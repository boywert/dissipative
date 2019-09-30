/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm_proto.h
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

void sidm_SetAndExchangeData(void);
void sidm_Scatter(void);
void scatter_stats(void);
void state_stats(void);
void sidm_AssignScatterPartner(void);
void sidm_NgbList(void);
void sidm_DoScatter(void);
void sidm_findHsml(void);
void sidm_findHsml_evaluate(int target, int mode, int threadid);
#if defined(SIDM_CONST_CROSS) || defined(SIDM_MAXWELLIAN)
MyDouble sidm_cross_sigma(MyDouble rel_vel, unsigned char reaction);
#else
MyDouble sidm_cross_sigma(MyDouble rel_vel, unsigned char reaction);
#endif
double sidm_scatter_P(double phys_rho, double phys_rel_vel, double Ekin, unsigned char reaction, int *retval);
void init_cross_table(void);
void sidm_Init_Particles(void);
void sidm_Init_CrossSection(void);
void sidm_ReInit(void);
void sidm_check_particle_scatter(void);
void sidm_NgbList(void);
void sidm_Scatter_evaluate(int target, int mode, int threadid);
void sidm_NgbList_evaluate(int target, int mode, int threadid);
void load_scatter_matrix_and_state(void);
void sidm_Init_States(void);
void sidm_get_random_unit_sphere_vector(double randnum1, double randnum2, double *xunit);
void sidm_set_velocities(scatter_process_data_in *sprdata_in, scatter_process_data_out *sprdata_out, double dvabs_in, double restitution, double *vcm_in, double *xunit);
void sidm_get_velocities(scatter_process_data_in *sprdata_in, scatter_process_data_out *sprdata_out, double *vcm_in, double *dv_in, double *dvabs_in);
void sidm_evaluate_scatter_process(scatter_process_data_in *sprdata_in, scatter_process_data_out *sprdata_out, MyIDType pID);
void sidm_scatter_in_to_out(scatter_process_data_in *sprdata_in, scatter_process_data_out *sprdata_out, MyIDType pID);
void sidm_get_delta_energy(unsigned char reaction, double *delta_E_1, double *delta_E_2);
void sidm_get_delta_mass(unsigned char reaction, double in_mass1, double in_mass2, double *out_mass1, double *out_mass2);
void sidm_SetGroundStateMass();
double sidm_calc_kinetic_energy_sidm_parts(void);
double sidm_calc_kinetic_energy_all_parts(void);
