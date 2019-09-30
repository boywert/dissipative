/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/becdm/becdm.h
 * \date        02/2017
 * \author      Philip Mocz
 * \brief       includes functions for BECDM pseudo-spectral algorithm
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef BECDM_H
#define BECDM_H
 
#ifdef BECDM

/* Becdm Main Methods */
void do_becdm_potential_kick(int i, double dt_gravkick);
void do_becdm_kinetic_drift(void);

/* helpers */
void particle2fftgrid(fft_plan *myplan);
void fftgrid2particle(fft_plan *myplan);

#endif

#endif /* BECDM_H */

