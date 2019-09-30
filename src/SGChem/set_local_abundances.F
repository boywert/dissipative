#include "../../build/arepoconfig.h"

      subroutine set_local_dust_abundance(DustToGasRatio)
#ifdef SGCHEM

#include "cool.h"
      REAL DustToGasRatio

#ifdef SGCHEM_VARIABLE_Z
c Set value in coolr common block using value passed in from sgchem.c
      dust_to_gas_ratio = DustToGasRatio

#endif

#endif 
      return
      end


      subroutine set_local_element_abundances(CarbAbund, OxyAbund, 
     $                                        MAbund, ZAtom)
#ifdef SGCHEM

#include "cool.h"
      REAL CarbAbund, OxyAbund, MAbund, ZAtom

#ifdef SGCHEM_VARIABLE_Z
      abundc = CarbAbund
      abundo = OxyAbund
      abundM = MAbund
      Z_atom = ZAtom

c For Si, take local ISM value and then scale by gas-phase metallicity
      abundsi = 1.5d-5 * Z_atom
#endif

#endif
      return
      end


      subroutine set_local_isrf(local_ISRF)
#ifdef SGCHEM

#include "cool.h"
      REAL local_ISRF

#ifdef SGCHEM_VARIABLE_ISRF
      G0 = local_ISRF
#endif

#endif
      return
      end


      subroutine set_local_cr_ion_rate(CR_energy_density,
     $                                 local_CR_ion_rate)
#ifdef SGCHEM

#include "cool.h"
      REAL local_CRIR, CR_energy_density

#ifdef SGCHEM_VARIABLE_CRION
#ifdef COSMIC_RAYS
c When running with cosmic rays enabled, we use the local energy density
c of cosmic rays to estimate the local cosmic ray ionization rate
      REAL energy_ratio
      energy_ratio = CR_energy_density / eV / CRenden_solcirc
      cosmic_ray_ion_rate = CRionrate_solcirc * energy_ratio
#else
      cosmic_ray_ion_rate = local_CR_ion_rate
#endif
#endif
#endif
      return
      end
