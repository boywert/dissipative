/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sidm/sidm_vars.c
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

#include "sidm_vars.h"

int *TargetList;

int Nforces;

int SIDM_NumScatterParticles;

double SIDM_Ekin_before_total[SIDM_REACTIONS], SIDM_Ekin_after_total[SIDM_REACTIONS], SIDM_Scattered_Mass[SIDM_REACTIONS];

MyFloat SIDM_arho, SIDM_avel;

double SIDM_clight;

double SIDM_GroundStateMass;

struct sidm_scatter_matrix SMSIDM[SIDM_REACTIONS];

struct sidm_scatter_states STSIDM[SIDM_STATES];

struct sidm_list *LISIDM;

struct sidm_data_p *PSIDM;

struct sidm_data_t *TSIDM;
