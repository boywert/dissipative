/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/tgchem/tgchem_test.h
 * \date        01/2013
 * \author      Thomas Greif
 * \brief       Primordial chemistry and cooling network
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>



#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include "tgchem.h"
#include "tgchem_proto.h"
#include "tgchem_test_proto.h"

#define MAX_STRING_LEN 1000

#define GAMMA (5. / 3.)
#define GAMMA_MINUS1 (GAMMA - 1.)
#define HYDROGEN_MASSFRAC 0.76
#define HE_ABUND ((1. / HYDROGEN_MASSFRAC - 1.) / 4.)
#define GRAVITY 6.6738e-8
#define BOLTZMANN 1.38065e-16
#define CLIGHT 2.99792458e10
#define PLANCK 6.6260695e-27
#define PROTONMASS 1.67262178e-24
#define ELECTRONMASS 9.1093829e-28

#define mymalloc(x, y) malloc(y)
#define mymalloc_movable(x, y, z) malloc(z)

#define myrealloc(x, y) realloc(x, y)
#define myrealloc_movable(x, y) realloc(x, y)

#define myfree(x) free(x)
#define myfree_movable(x) free(x)

#define terminate(...) {char termbuf1[MAX_STRING_LEN], termbuf2[MAX_STRING_LEN]; sprintf(termbuf1, "Code termination in function %s(), file %s, line %d", __FUNCTION__, __FILE__, __LINE__); sprintf(termbuf2, __VA_ARGS__); printf("%s\n%s\n", termbuf1, termbuf2); fflush(stdout); exit(0);}

extern int ThisTask;
extern int NTask;

extern double WallClockTime;
