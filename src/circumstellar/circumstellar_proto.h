/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/circumstellar/circumstellar_proto.h
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

#ifdef CIRCUMSTELLAR


double get_circumstellar_distance(int i);
double circumstellar_DoCoolingHeating(int i, double dtime);
double circumstellar_irradiation_utherm(int i);
double circumstellar_convert_temp_to_u(double T);
void source_particle_create_list();
void source_particle_update_list();

void circumstellar_calc_gravity_from_stars_planets_only(void);

double get_circumstellar_alpha_viscosity(double x, double y, double z, double rho, double press);

#ifdef CIRCUMSTELLAR_WBOUNDARIES
void do_circumstellar_disk_update(struct particle_data *localP, struct sph_particle_data *localSphP, int p);
int set_vertex_vel_circumstellar(struct particle_data *localP, int p);
#endif


#endif
