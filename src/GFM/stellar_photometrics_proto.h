/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/stellar_photometrics_proto.h
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

void get_stellar_photometrics(MyFloat ageInGyr, MyFloat metallicity, MyFloat mass, stellar_photometrics * st_photo);
void init_stellar_photometrics(void);
void assign_stellar_photometrics(int i, stellar_photometrics * st_photo);
void test_stellar_photometrics(double mass, double metallicity, double log_AgeInGyr_min, double log_AgeInGyr_max, int log_AgeInGyr_bins);
void generate_random_directions(void);
