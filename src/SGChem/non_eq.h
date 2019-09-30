c
c Written by S. Glover, AMNH, 2004-2005, AIP, 2006-2007
c
#ifdef SGCHEM
#include "sgchem_def.h"

#ifdef CHEMISTRYNETWORK
      integer nchem_network
      parameter (nchem_network = CHEMISTRYNETWORK)
#endif
c
c Set up quantities (such as the absolute tolerances) that are used in 
c multiple places in the non-equilibrium chemistry code. Note that most 
c DVODE-specific setup should go in evolve_abundances.F -- nrpar & nipar
c are exceptions, as they are used elsewhere, so it is useful to define
c them here
c 
      integer nrpar, nipar
      parameter (nrpar = 9)
      parameter (nipar = 2)
c
      integer num_non_eq_species
      parameter (num_non_eq_species = SGCHEM_NUM_SPECIES)
c
      integer nspec 
      parameter (nspec = num_non_eq_species+1)

#if CHEMISTRYNETWORK == 1
      integer num_eqb_species
      parameter (num_eqb_species = 3)
#endif

c
c Amount by which abundances are allowed to stray over their theoretical
c maximum before triggering an error in rate_eq -- set to a blanket value
c of 1d-4 for the time being...
c 
      REAL eps_max
      parameter (eps_max = 1d-4)
c
      REAL atol(nspec)
c
#define  RTOL      1d-4
#define  ATOL_H2   1d-9
#define  ATOL_HP   1d-9
#define  ATOL_DP   1d-10
#define  ATOL_HD   1d-10
#define  ATOL_HEPP 1d-14
#define  ATOL_CP   1d-16
#define  ATOL_HEP  1d-14
#define  ATOL_CO   1d-14
#define  ATOL_HCOP 1d-18
#define  ATOL_CHX  1d-14
#define  ATOL_OHX  1d-14
#define  ATOL_MP   1d-14
#define  ATOL_TMP  0d0
c
      common /tolerance/ atol
c
#endif /* SGCHEM */
