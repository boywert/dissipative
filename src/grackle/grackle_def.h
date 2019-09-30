/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/grackle/grackle_proto.h
 * \date        MM/YYYY
 * \author     	Matthew C Smith
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifdef COOLING
#ifdef GRACKLE

#ifdef GRACKLE_D
#define GRACKLE_SPECIES_NUMBER 11
#define GRACKLE_H2
#elif defined(GRACKLE_H2)
#define GRACKLE_SPECIES_NUMBER 8
#elif !defined(GRACKLE_TAB)
#define GRACKLE_SPECIES_NUMBER 5
#endif
 
void initialise_grackle(void);
void cooling_only(void);
void cool_active_cells(void);
double get_temp_individual_cell_grackle(int i);
double get_cooling_time_individual_cell_grackle(int i);
#ifndef GRACKLE_TAB
void grackle_initialise_abundances(void);
void grackle_converge_abundances(void);
#endif
#endif
#endif
