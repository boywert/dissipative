/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/stellar_feedback_proto.h
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

#ifndef FM_STELLAR_FEED_PROTO_H
#define FM_STELLAR_FEED_PROTO_H

MyDouble compute_SN_energy(MyDouble Number_of_SN);
MyDouble quadratic_equation_solver(MyDouble b, MyDouble c);

#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
MyDouble compute_SN_momentum(MyDouble Number_of_SN);
#endif

#ifdef DELAYED_COOLING
void init_delayed_cooling();
void set_cooling_shutoff_time(int i, MyFloat CoolShutoffTime);
MyFloat compute_blast_radius(MyFloat energy, MyFloat avg_H_density, MyFloat avg_pressure);
MyFloat compute_cooling_shutoff_time(MyFloat energy, MyFloat avg_H_density, MyFloat avg_pressure);
MyFloat update_cooling_shutoff_time(int i, MyFloat dt);
#endif

#ifdef DELAYED_COOLING_TURB
void init_turbulent_energy();
void convert_specific_turbulent_energy_to_turbulent_energy();
int check_if_turbulent_energy_above_threshold(int i);
MyFloat dissipate_turbulent_energy(MyFloat Uturb_old, MyFloat dt);
#endif

#ifdef INSTANTANEOUS_DEPOSITION
void init_instantaneous_deposition();
MyDouble compute_total_SN_energy(int index, MyDouble Age, MyFloat * SN_Number);
#endif

#endif
