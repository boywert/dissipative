/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/sinks/sinks_proto.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Sink particles
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

void sinks(void);

void sinks_dmass(void);
void sinks_accrete(void);

void sinks_set_constants(void);
void write_sink_data(int xaxis, int yaxis, int zaxis, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

void sinks_create(void);

void sinks_begrun(void);

void sinks_get_num_sinks(void);
void sinks_get_active_sinks(void);
void sinks_begin(void);
void sinks_end(void);

