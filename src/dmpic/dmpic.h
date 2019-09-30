/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/dmpic/dmpic.h
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

#ifndef DMPIC_H
#define DMPIC_H

void dmpic_init(void);
void dmpic_geometry(int nr);
void dmpic_save_picture(int nr);
void dmpic_project_and_smooth(int nr);
integertime dmpic_get_next_pictime(void);
void dmpic_make_image(void);

#endif
