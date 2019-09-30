#ifndef ADD_BGGRID_H
#define ADD_BGGRID_H

#include "../allvars.h"

#ifdef ADDBACKGROUNDGRID

#define ADDBACKGROUNDGRIDMAX  256
#define FACTOR_MAX_BOX_SIZE   15.0
#define FACTOR_MIN_BOX_SIZE   2.0

extern MyIDType IDNew;

int add_backgroundgrid(void);
void prepare_domain_backgroundgrid(void);
void calculate_weights();
void distribute_particles();

#endif

#endif /* ADD_BGGRID_H */
