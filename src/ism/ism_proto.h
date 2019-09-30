/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/ism/ism_proto.h
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

#ifndef ISM_PROTO_H
#define ISM_PROTO_H

void init_star_formation();
double get_ism_h2_frac(int i);

double evaluate_l_over_m_ssp(double stellar_age_in_gyr);




void ism_find_enclosed_gas_mass(void);  /* find total star forming gas mass enclosed between P[i].pos and SphP[i].ClumpCenter. Store in SphP[i].ClumpMass       */
void ism_find_enclosed_stellar_luminosity(void);        /* find total stellar luminosity enclosed in same region.  Store in SphP[i].ClumpLuminosity                             */
void ism_kick_particles(void);


void ism_find_clump_centers(void);      /* loop over all (star forming) gas particles, set field SphP[i].ClumpCenter to clump center location                   */
int find_clump_centers_evaluate(int target, int mode, int *nexport, int *nsend_local);

void ism_find_clump_radii(void);        /* loop over all (star forming) gas particles, set field SphP[i].ClumpCenter to clump center location                   */
int find_clump_radii_evaluate(int target, int mode, int *nexport, int *nsend_local);

#endif
