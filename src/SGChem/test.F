#include "../../build/arepoconfig.h"
      program test
#ifdef SGCHEM
#define NSIDE 2
#define NPIX  12*NSIDE*NSIDE
#include "sgchem_def.h"
      implicit NONE
      real*8 abundD, abundc, abundo, abundM, Zatom, tdust, G0,
     $       LWBGz, D2G, crion, initz, h2ex, h2kin, AVeff, time,
     $       dl, yn, divv, energy, current_redshift
      real*8 non_eq_abundances(SGCHEM_NUM_SPECIES)
      real*8 column_density_projection(NPIX)
      real*8 column_density_projection_H2(NPIX)
      real*8 column_density_projection_CO(NPIX)
      real*8 thermal_rates_out(SGCHEM_NUM_THERMAL_RATES)

      integer LWBG, iphoto, isrf, iatom, id, I

      abundD = 2.6d-5
      abundc = 5.66596682d-05
      abundo = 1.29507811d-04
      abundM = 4.04711891d-08
      Zatom  = 4.04711932d-01
      tdust  = 15.21233845
      G0     = 1.7d0
      LWBG   = 1
      LWBGz  = 0d0
      D2G    = 0.40471193
      crion  = 3d-17
      initz  = 0d0
      AVeff  = 0d0
      h2ex   = 0.84d0
      h2kin  = 0.12d0
      iphoto = 6
      isrf   = 1
      iatom  = 1

      call init_chemistry_parameters(abundD,
#ifndef SGCHEM_VARIABLE_Z
     $ abundc, abundo, abundM, Zatom,
#endif
     $ tdust, G0, LWBG, LWBGz,
#ifndef SGCHEM_VARIABLE_Z
     $ D2G,
#endif
#ifndef SGCHEM_VARIABLE_CRION
     $ crion,
#endif
     $ initz, AVeff, h2ex, h2kin, iphoto, isrf, iatom)

      call coolinmo
      call cheminmo
      call init_tolerances

#ifdef SGCHEM_VARIABLE_Z
      call set_local_dust_abundance(D2G)
      call set_local_element_abundances(abundc, abundo, abundM,
     $                                  Zatom)
#endif

      time = 461205148662.6d0
      dl   = 1.035280993819809d20
      yn   = 4.834646667854838d-03
      divv = 8.049198965407391d-15
      energy   = 4.209412687868616d-13
      current_redshift = 0d0
c Current abundances
      non_eq_abundances(IH2+1) = 9.999999999999999D-21
      non_eq_abundances(IHP+1) = 1d0
      non_eq_abundances(ICO+1) = -4.620165234693153D-06
c Columns from Treecol
      do I = 1, NPIX
        column_density_projection(I) = 0d0
        column_density_projection_H2(I) = 0d0
        column_density_projection_CO(I) = 0d0
      enddo
      id   = 1

      call evolve_abundances(time, dl, yn, divv, energy,
     $                       current_redshift, non_eq_abundances,
     $                       thermal_rates_out,
     $                       column_density_projection,
     $                       column_density_projection_H2,
     $                       column_density_projection_co,
     $                       tdust,id)

      print*, non_eq_abundances, energy

#endif
      end
