#include "../../build/arepoconfig.h"
      subroutine init_chemistry_parameters(DeutAbund,
#ifndef SGCHEM_VARIABLE_Z
     $                                     CarbAbund,
     $                                     OxyAbund,
     $                                     MAbund,
     $                                     ZAtom,
#endif
     $                                     InitDustTemp, 
#ifndef SGCHEM_VARIABLE_ISRF
     $                                     UVFieldStrength,
#endif
     $                                     LWBGType,
     $                                     LWBGStartRedsh,
#ifndef SGCHEM_VARIABLE_Z
     $                                     DustToGasRatio, 
#endif 
#ifndef SGCHEM_VARIABLE_CRION
     $                                     CosmicRayIonRate, 
#endif
     $                                     InitRedshift,
     $                                     ExternalDustExtinction, 
     $                                     H2FormEx, H2FormKin,
     $                                     PhotoApprox, ISRFOption,
     $                                     AtomicCoolOption,
     $                                     H2OpacityOption)
#ifdef SGCHEM
#include "cool.h"
      REAL InitDustTemp, DeutAbund
      REAL LWBGStartRedsh
#ifndef SGCHEM_VARIABLE_Z
      REAL CarbAbund, OxyAbund, ZAtom, MAbund
      REAL DustToGasRatio
#endif
#ifndef SGCHEM_VARIABLE_ISRF
      REAL UVFieldStrength
#endif
#ifndef SGCHEM_VARIABLE_CRION
      REAL CosmicRayIonRate
#endif
      REAL InitRedshift
      REAL ExternalDustExtinction, H2FormEx, H2FormKin
      integer PhotoApprox, ISRFOption, AtomicCoolOption, LWBGType
      integer H2OpacityOption

c Coolr parameters initialized via argument list:
      tdust  = InitDustTemp
#ifndef SGCHEM_VARIABLE_ISRF
      G0     = UVFieldStrength
#endif
      lwstartz = LWBGStartRedsh
#ifndef SGCHEM_VARIABLE_Z
      abundc = CarbAbund
      abundo = OxyAbund
      abundM = MAbund
      Z_atom = ZAtom
      dust_to_gas_ratio   = DustToGasRatio
#endif
#ifndef SGCHEM_VARIABLE_CRION
      cosmic_ray_ion_rate = CosmicRayIonRate
#endif
      redshift            = InitRedshift
      AV_ext              = ExternalDustExtinction
      h2_form_ex          = H2FormEx
      h2_form_kin         = H2FormKin

c Coolr parameters initialized here:

      abundD  = DeutAbund
      abundN  = 0d0
#ifndef SGCHEM_VARIABLE_Z
      abundsi = 1.5d-5 * Z_atom   ! Scaled from local ISM value
#endif
      f_rsc   = 1d0
      deff    = 1d0
      phi_pah = 0.5d0
      AV_conversion_factor = 5.348d-22

c Cooli parameters initialized via argument list
      iphoto = PhotoApprox
      isrf_option = ISRFOption
      iflag_atom  = AtomicCoolOption
      lwtype = LWBGType
      iflag_H2_opacity = H2OpacityOption

c Cooli parameters initialized here
c ID and index parameters will be overridden later
      id_current    = -1
      index_current = -1

c Select the correct set of rate coefficients for primordial
c runs. Note that these reactions aren't used in the ISM networks,
c so the values we choose here generally don't matter in that case.
c The exception is iflag_3bh2a - we set this to 5 to ensure that the
c MSM96 rate is selected
      iflag_mn = 4
      iflag_ad = 6
#if CHEMISTRYNETWORK != 1
      iflag_3bh2a = 5
#else
      iflag_3bh2a = 6 
#endif
      iflag_3bh2b = 4 
      iflag_h3pra = 1

#endif
      return
      end

      subroutine set_index_id_for_chemistry(index, id)
#ifdef SGCHEM
#include "cool.h"
      integer index, id

      id_current    = id
      index_current = index
#endif
      return
      end
