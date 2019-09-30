/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_feedback_kernels.h
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

#ifndef STELLAR_FEEDBACK_KERNELS_H
#define STELLAR_FEEDBACK_KERNELS_H

typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat NormSph;
  MyFloat TotalMassReleased;
#if defined(FM_MASS_WEIGHT_SN) || defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK)
  MyFloat TotNgbMass; 
#endif
#ifdef GFM_WINDS_LOCAL
  MyFloat WindEnergyReleased;
#endif
#ifdef GFM_STELLAR_FEEDBACK
  MyDouble SNIaEnergyReleased;
  MyDouble AGBMomentumReleased;
#endif
#ifdef FM_STAR_FEEDBACK
  MyDouble TotalEnergyReleased;
#ifdef FM_SN_COOLING_RADIUS_BOOST
  MyFloat n_SNII;		/* these are floats; sampled to ints during feedback routine */
  MyFloat n_SNIa;
  MyFloat LocISMdens;
  MyFloat LocISMZdens;
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
  MyDouble TotalMomentumReleased;
#endif
#if defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && defined(INJECT_INTO_SINGLE_CELL)
  MyIDType ClosestNeighbourID;
#endif
#ifdef DELAYED_COOLING
  MyFloat NormSphFeedback;
  MyFloat BlastRadius;
  MyFloat CoolShutoffTime;
#else
#ifndef DELAYED_COOLING_TURB
  MyFloat NumNgb;
#endif
#endif
#endif
  int Firstnode;
} data_in;

typedef struct
{
#ifndef DELAYED_COOLING_TURB
#if !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
  int kicked_cells;
#ifdef DELAYED_COOLING
  MyFloat BlastRadius;
  MyFloat ShutoffTime;
#else
  MyDouble mass_to_kick;
  MyDouble mass_kicked;
  MyFloat kick_vel;
  MyDouble deltaEKin;
  MyDouble deltaMomentum[3];
#endif
#endif
#endif
#if !defined(DELAYED_COOLING_TURB) && !defined (DELAYED_COOLING)
  MyDouble TotalMomentumInjected;
#endif
} data_out;

int is_doing_stellar_feedback(int i);

#ifdef GFM_STELLAR_FEEDBACK
void GFM_stellar_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in);
#endif

#ifdef GFM_WINDS_LOCAL
void GFM_winds_local(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in);
#endif

#ifdef FM_STAR_FEEDBACK

#ifdef DELAYED_COOLING
void delayed_cooling(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif
#ifdef DELAYED_COOLING_TURB
void delayed_cooling_turbulence(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in);
#endif
#ifdef NON_STOCHASTIC_MOMENTUM_FEEDBACK
void momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif
#ifdef DIRECT_MOMENTUM_INJECTION_FEEDBACK
void direct_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif
#ifdef FM_SN_COOLING_RADIUS_BOOST
void cooling_radius_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif
#if !defined(DELAYED_COOLING) && !defined(DELAYED_COOLING_TURB) && !defined(NON_STOCHASTIC_MOMENTUM_FEEDBACK) && !defined(DIRECT_MOMENTUM_INJECTION_FEEDBACK)
void stochastic_momentum_feedback(int target, int mode, int thread_id, int numnodes, int *firstnode, data_in *in, data_out *out);
#endif

#endif

#endif
