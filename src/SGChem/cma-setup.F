#include "../../build/arepoconfig.h"
#ifdef MCMA      
      blockdata cma_setup
      implicit NONE
      integer nelem_cma, nspec_cma

#if CHEMISTRYNETWORK == 15
      parameter (nelem_cma = 3)
      parameter (nspec_cma = 10)

      integer cma_atoms_y_in_x(nspec_cma, nelem_cma)
c
c H2, H+, C+, CHx, OHx, CO, HCO+, H, C, O
c
      DATA cma_atoms_y_in_x /2,1,0,1,1,0,1,1,0,0,
     $                       0,0,1,1,0,1,1,0,1,0,
     $                       0,0,0,0,1,1,1,0,0,1/

      integer cma_total_atoms(nspec_cma)
      DATA cma_total_atoms  /2,1,1,2,2,2,3,1,1,1/
#endif

      real*8 cma_weight(nspec_cma, nelem_cma, nelem_cma)
      real*8 cma_nonzero(nelem_cma, nelem_cma)

      common /cma_data/ cma_atoms_y_in_x, cma_total_atoms, 
     $                  cma_weight, cma_nonzero

      end
#endif /* MCMA */
