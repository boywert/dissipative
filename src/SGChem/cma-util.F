#include "../../build/arepoconfig.h"
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            C M A _ I N I T                \\\\\\\\\\
c
c=======================================================================
      subroutine cma_init
#ifdef MCMA
      implicit none
#include "sgchem_def.h"
#include "cma.h"
#include "cool.h"
      integer I, J, K

      do I = 1, nelem_cma
        do J = 1, nelem_cma
          cma_nonzero(I,J) = 0
          do K = 1, nspec_cma
            if (cma_atoms_y_in_x(K,J) .gt. 0d0) then
              cma_weight(K,I,J) = (cma_atoms_y_in_x(K,I)
     $                          / real(cma_total_atoms(K)))
     $                          * cma_atoms_y_in_x(K,J)
              cma_nonzero(I,J) = 1
            else
              cma_weight(K,I,J) = 0d0
            endif
          enddo
        enddo
      enddo
#endif
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////             C M A _ I N I T             \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            C M A _ M A T R I X            \\\\\\\\\\
c
c=======================================================================
      subroutine cma_matrix(flux, M)
#ifdef MCMA
      implicit none
#include "sgchem_def.h"
#include "cma.h"
      real*8 flux(nspec_cma)
      real*8 M(nelem_cma, nelem_cma)
      integer I, J, K
c
c  Mij = Sum_k [N(k,i) / N(k)] * flux(k,j)
c  flux(k,j) = flux(k)   if N(k,j) > 0
c            = 0         if N(k,j) = 0 
c
c  We define a function cma_weight(k,i,j) such that:
c
c  cma_weight(k,i,j) = [N(k,i) / N(k)] * truth(N(k,j) > 0) 
c
c where truth(X) = 1 if X is true and 0 if X is false
c
      do I = 1, nelem_cma
        do J = 1, nelem_cma
          M(I,J) = 0d0
          if (cma_nonzero(I,J) .gt. 0) then
            do K = 1, nspec_cma
              M(I,J) = M(I,J) + cma_weight(K,I,J) * flux(K)
            enddo
          endif
        enddo
      enddo
#endif
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////            C M A _ M A T R I X            \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////           C M A _ C O R R E C T           \\\\\\\\\\
c
c=======================================================================
      subroutine cma_correct(species_flux, element_flux, nspec_input)
#ifdef MCMA
      implicit none
#include "sgchem_def.h"
#include "cma.h"
      REAL species_flux(nspec_cma)
      REAL element_flux(nelem_cma)
      REAL m_cma(nelem_cma, nelem_cma)
      REAL flux_cma(nelem_cma)
      integer I, info, ipiv(nelem_cma), nspec_input

      if (nspec_input .ne. nspec_cma) then
        print*, "Error: wrong number of species in CMA"
        print*, "Got ", nspec_input, " Expected ", nspec_cma
        stop
      endif

      call cma_matrix(species_flux, m_cma)
      do I = 1, nelem_cma
        flux_cma(I) = element_flux(I)
      enddo

      call DGESV(nelem_cma, 1, m_cma, nelem_cma, ipiv, flux_cma, 
     $           nelem_cma, info)
      if (info .ne. 0) then
        print*, "Error: MCMA failed, error code ", info
        print*, "Species flux: ", species_flux
        print*, "CMA matrix: ", m_cma
        stop
      endif             

      call cma_flux_correct(flux_cma, species_flux)
#endif
      return
      end
c=======================================================================
c
c    \\\\\\\\\\         E N D   S U B R O U T I N E       //////////
c    //////////           C M A _ C O R R E C T           \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      C M A _ F L U X _ C O R R E C T      \\\\\\\\\\
c
c=======================================================================
      subroutine cma_flux_correct(flux_cma, species_flux)
#ifdef MCMA
      implicit none
#include "sgchem_def.h"
#include "cma.h"
      real*8 flux_cma(nelem_cma), species_flux(nspec_cma)
      real*8 eta_H, eta_C, eta_O, correction
      integer I

#if CHEMISTRYNETWORK == 15
      eta_H  = flux_cma(1)
      eta_C  = flux_cma(2)
      eta_O  = flux_cma(3)

      do I = 1, nspec_cma
        correction = (eta_H  * cma_atoms_y_in_x(I,1)
     $             +  eta_C  * cma_atoms_y_in_x(I,2)
     $             +  eta_O  * cma_atoms_y_in_x(I,3))
     $             /  cma_total_atoms(I)
        species_flux(I) = species_flux(I) * correction
      enddo
#endif
c
#endif
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////      C M A _ F L U X _ C O R R E C T      \\\\\\\\\\
c
c=======================================================================
