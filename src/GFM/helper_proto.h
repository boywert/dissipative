/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/helper_proto.h
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

void gfm_add_star(int i, int j, MyDouble mass_of_star, MyFloat birthtime, MyFloat hsml_guess);
void gfm_calc_wind_parameters(int i, double p, double *p_wind, double *v_wind, double *u_wind);
void gfm_add_wind(int i, double v, double utherm);
void gfm_spawn_wind_from_cell(int igas, double v, double u, int istar, MyFloat mass_of_wind);
void check_AuxDataID_references(void);
void gfm_inject_into_cell(int j, double mass, double inj_thermalenergy, double *inj_mom);
void gfm_cap_dust(void);
void gfm_check_dust(int checkid);
