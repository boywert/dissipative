/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/healray/healray_proto.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Adaptive ray-tracing
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

void healray();

void healray_comm(int mode);
void healray_split_comm();

void healray_merge();
void healray_mic_test();

void healray_finish_step();

void healray_begrun();
void healray_init();
void healray_init_step();
void healray_init_rays();
void healray_init_rays_alt();

void healray_init_rayout();
void healray_finish_rayout();

void healray_read_sources();
void healray_update_sources();
void healray_init_sources();
