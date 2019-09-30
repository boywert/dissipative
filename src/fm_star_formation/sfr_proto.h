/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/sfr_proto.h
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

#ifndef FM_STAR_FORM_PROTO_H
#define FM_STAR_FORM_PROTO_H

void init_star_formation();

#ifdef OUTPUT_SF_PROBABILITY
double compute_sf_probability(int index);
#endif

#endif
