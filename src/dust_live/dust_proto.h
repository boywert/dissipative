/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/dust_live/dust_proto.h
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

#ifndef DUST_PROTO_H
#define DUST_PROTO_H

void init_dust(void);
void init_grain_size_distribution(void);

void start_dust(void);
void end_dust(void);

void update_timebins_dust(void);
integertime get_timestep_dust(int p);

void find_drag_cells(int npart);
int find_drag_cells_evaluate(int target, int mode, int thread_id);

void drag_kernel(void);
int drag_kernel_evaluate(int target, int mode, int thread_id);

void compute_drag_acceleration(void);
void do_drag_step(void);
void do_drag_step_first_half(void);
void do_drag_step_second_half(void);

void dust_transfer_kernel(void);
int dust_transfer_kernel_evaluate(int target, int mode, int thread_id);

void dust_findHsml(void);
void dust_findHsml_evaluate(int target, int mode, int threadid);

void do_dust_transfer(void);
void do_dust_transfer_sne(void);
void update_dust_element_fractions(int i);
void update_grain_sizes(int k, double adot, double dt, double dt_code);
void update_sn_rates(void);
void update_grain_sizes_sne(int k, double dt, double dt_code);
void correct_grain_conserved_quantities(int k);
void check_dust_for_removal(int i);

void begin_shattering(void);
void end_shattering(void);
void exchange_shattering_results(void);
void do_shattering_coagulation(void);
void update_grain_sizes_shattering(int l, double dt, double dt_code, enum gsd_collision_type type);

double interval_overlap(double a_1, double b_1, double a_2, double b_2);
double bin_avg_mass(int j, double num_grains, double slope);
int elem_can_be_dust(int k);

void create_dust_particles(void);
void spawn_dust_from_star(int istar, int idust, double m_dust, double m_dust_each, int ispawn);

void refine_dust_particles(void);
void subdivide_dust_particle(int iorig, int idust);

#endif
