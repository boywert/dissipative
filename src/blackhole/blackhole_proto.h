/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/blackhole/blackhole_proto.h
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

double blackhole_mdot_eddington(double bh_mass);
double blackhole_luminosity_eddington(double bh_mass);
double blackhole_bondi_rate(int n, int type);
void blackhole_drag_force(void);
void blackhole_mark_cells_for_refinement(void);
void blackhole_calculate_mdot(void);
void blackhole_reset_energy(void);
void blackhole_accumulate_energy(void);
void blackhole_do_bubbles(int num_activebh);
void blackhole_massiveseeds(int nseed, int *seed_indices);
void remove_swallowed_gas(int nseed, int *seed_indices);
int blackhole_gasngb_evaluate(int target, int mode, int *nexport, int *nsend_local);
void blackhole_do_mergers(void);
void blackhole_swallow_gas(void);
void blackhole_disk_vorticity(void);
void blackhole_assign_feedback(void);
void blackhole_find_neighboring_holes_and_potmin(void);
void blackhole_calculate_mdot_radiomode(void);
void blackhole_do_bubbles_newradio(int num_activebh);
double blackhole_get_bubble_energy_thresh(int n);
void bh_newradio_bubble(MyDouble bh_egy, MyDouble center[3], MyDouble rvir_comov, MyFloat bubble_szie, MyIDType BH_id);
void blackhole_place_newradio_bubbles(int num_activebh);
int blackhole_place_newradio_bubbles_evaluate(int target, int mode, int rep, int *nexport, int *nsend_local);
void bh_bubble(MyDouble bh_dmass, MyDouble center[3], MyIDType BH_id);
void blackhole_energy_log_info(void);
void blackhole_reposition(void);
int blackhole_reposition_evaluate(int target, int mode, int *nexport, int *nsend_local);
int ngb_treefind_dyn_friction(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);
void blackhole_dyn_friction(void);
int blackhole_isactive(int n);
void blackhole_density(void);
void blackhole_accretion(void);
int blackhole_check_duty_cycle(double delta_temp, double delta_time);
void blackhole_potmin_diagnostic(void);
void blackhole_friction_update_vel_pot_minimum(void);
void blackhole_friction_store_previous_minimum(void);
// Blackhole Bipolar Feedback
void blackhole_bipolar(void);
int is_cell_in_bipolar_cone(MyIDType i, MyDouble *pos, MyDouble *bipolar_j);
int is_cell_in_bipolar_cold_disk(MyIDType i, int ColdDisk);
// Blackhole Spin Evolution
void cross_product(MyFloat *A, MyFloat *B, MyFloat *C);
void blackhole_spin_evolution(void);
void blackhole_new_accretion_episode(int i, MyFloat CurrentTime, MyFloat BH_Mass);
