/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/GFM/cooling_metal_proto.h
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

void read_cooling_tables_init(void);
void read_cooling_tables_current_time(void);
void init_cooling_metal(void);
void dump_cooling_table();
#ifdef GFM_AGN_RADIATION
MyFloat get_CoolingMetalRate(MyFloat log_BolFlux, MyFloat log_MetallicityInSolar, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature);
#else
#ifdef RADCOOL
MyFloat get_CoolingMetalRate(MyFloat log_Phios, MyFloat log_Phins
#ifdef RADCOOL_HOTHALO
                             , MyFloat log_PhiT6, MyFloat log_PhiT7, MyFloat log_PhiT8
#endif
                             , MyFloat log_MetallicityInSolar, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature);
#else
MyFloat get_CoolingMetalRate(MyFloat log_Metallicity, MyFloat log_HydrogenNumberDensity, MyFloat log_Temperature);
#endif
#endif
