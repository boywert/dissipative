#include "../../build/arepoconfig.h"
#ifdef MCMA
      integer nelem_cma, nspec_cma

      parameter (nelem_cma = NELEM_CMA)
      parameter (nspec_cma = NSPEC_CMA)

c NB cma_atoms_y_in_x and cma_total_atoms are initialized
c in the cma-setup BLOCKDATA subroutine
      integer cma_atoms_y_in_x(nspec_cma, nelem_cma)
      integer cma_total_atoms(nspec_cma)
      real*8 cma_weight(nspec_cma, nelem_cma, nelem_cma)
      real*8 cma_nonzero(nelem_cma, nelem_cma)

      common /cma_data/ cma_atoms_y_in_x, cma_total_atoms, 
     $                  cma_weight, cma_nonzero

#endif /* MCMA */
