/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm_vars.h
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

#ifndef SIDM_VARS_H
#define SIDM_VARS_H

#define SIDM_MAX_NGBS 38 
#define SIDM_MIN_NGBS 10

#ifndef SIDM_STATES
#define SIDM_STATES 1
#endif

#ifndef SIDM_REACTIONS
#define SIDM_REACTIONS 1
#endif

extern int *TargetList;

extern int Nforces;

extern int SIDM_NumScatterParticles;

extern double SIDM_Ekin_before_total[SIDM_REACTIONS], SIDM_Ekin_after_total[SIDM_REACTIONS], SIDM_Scattered_Mass[SIDM_REACTIONS];

extern MyFloat SIDM_arho, SIDM_avel;

extern double SIDM_clight;

extern double SIDM_GroundStateMass;

typedef struct
{
  MyIDType NgbIDs;
  double P0j_half[SIDM_REACTIONS];
  float Distance;
  unsigned char State;
} ngb_entry;

typedef struct
{
  MyDouble Pos1[3], Pos2[3], Vel1[3], Vel2[3];
  MyDouble Mass1, Mass2;
  unsigned char State1, State2;
  unsigned char Reaction;
} scatter_process_data_in;

typedef struct
{
  MyDouble Vel1[3], Vel2[3];
  MyDouble Mass1, Mass2;
  unsigned char OutState1, OutState2;
} scatter_process_data_out; 


extern struct sidm_scatter_matrix
{
  int Out1, Out2, In1, In2;
}
 SMSIDM[SIDM_REACTIONS];

extern struct sidm_scatter_states
{
  double DeltaMass;
  double InitialFraction;
}
 STSIDM[SIDM_STATES];

extern struct sidm_list
{
  int List1, List2;
}
 *LISIDM;

extern struct sidm_data_p
{
  char ShouldScatterInStep[SIDM_REACTIONS];
  MyIDType ScatterID, ID;
  float Hsml, VelDisp[SIDM_STATES];
  double RandX;
  int NumNgb, NumNgbState[SIDM_STATES];
  float Vx[SIDM_STATES], Vy[SIDM_STATES], Vz[SIDM_STATES];
  float Density[SIDM_STATES];
  double PSum[SIDM_REACTIONS];
  unsigned char ScatterReaction;
  int EnergyForbidden[SIDM_REACTIONS];
  int ngb_Offset;
  ngb_entry ngb_Entry[SIDM_MAX_NGBS];
}
 *PSIDM;


extern struct sidm_data_t
{
  MyFloat NewVel[3];
  double NewMass;
  unsigned char NewState;
  unsigned char ScatterReaction;
  int ScattersInStep[SIDM_REACTIONS];
}
 *TSIDM;


#endif
