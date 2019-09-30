#define REAL real*8
#define UNUSED_PARAM(x)  x = x
#define ABORT(x) stop

#include "sgchem_def.h"
#ifdef SGCHEM
c
c He:H ratio by number (=> ratio by mass is 4*abhe)
c
      REAL abhe
      parameter(abhe = ABHE)
c
c Symbolic constants representing position of each species in
c the abundance array. Note that these values are for use in
c Fortran code and hence start from 1; the corresponding 
c #defines are for use in C code and start from 0.

      integer ih2
      parameter(ih2 = IH2+1)

      integer ihp
      parameter(ihp = IHP+1)

#if CHEMISTRYNETWORK == 1
      integer idp
      parameter(idp = IDP+1)

      integer ihd
      parameter(ihd = IHD+1)

      integer ihepp
      parameter(ihepp = IHEPP+1)
#endif

      integer icp
      parameter(icp  = ICP+1)

      integer ichx
      parameter(ichx = ICHX+1)

      integer iohx
      parameter(iohx = IOHX+1)

      integer ico
      parameter(ico  = ICO+1)

      integer ihcop
      parameter(ihcop = IHCOP+1)

      integer ihep
      parameter(ihep = IHEP+1)

      integer imp
      parameter(imp = IMP+1)

      integer itmp
      parameter(itmp = ITMP+1)

c Number of entries in cooling table
      integer nmd
      parameter(nmd = 10000)

c Number of cooling / heating rates computed in cooling fn.
      integer nrates
      parameter(nrates = 30)

c Number of cooling / heating rates computed in chemistry routines
      integer nrates_chem
      parameter(nrates_chem = 21)

c Total number of cooling / heating rates
      integer nrates_tot
      parameter(nrates_tot = nrates + nrates_chem)

c Number of abundances passed to cooling function
c (Note that this is not necessarily the same as the number
c that we actually track as field variables)
      integer nabn
      parameter(nabn = 17)

c Boltzmann constant
      REAL kboltz
      parameter (kboltz = 1.38066e-16)

c Hydrogen mass
      REAL mh
      parameter (mh = 1.6726e-24)

c One electron volt, in ergs
      REAL eV
      parameter (eV = 1.60219e-12)

c Gravitational constant
      REAL Gn
      parameter (Gn = 6.672e-8)

c Number of different quantities stored in cooling look-up table
      integer ncltab
      parameter (ncltab = 86) 

c Number of different quantities stored in chemistry look-up table
      integer nchtab
      parameter (nchtab = 173) 

c Number of cosmic ray ionization rates tabulated
      integer ncrtab
      parameter(ncrtab = 13)

c Number of cosmic ray induced photoionizations/photodissociations tabulated
      integer ncrphot
      parameter(ncrphot = 12)

c Number of photochemical rates tabulated
      integer nphtab
      parameter(nphtab = 54)

c Size of rate array returned by calc_photo - NB not all of the entries are filled
      integer npr
      parameter (npr = 54)

c Number of constant rate coefficient initialized in const_rates
      integer nconst
      parameter(nconst = 91)

c These variables are initialized in cheminmo
      REAL chtab(nchtab,nmd), dtchtab(nchtab,nmd),
     $     crtab(ncrtab), crphot(ncrphot)

c These variables are initialized in photoinit
      REAL phtab(nphtab), f_rsc

c This is initialized in const_rates
      REAL cst(nconst)

c These variables are initialized in coolinmo
      REAL temptab(nmd)
      REAL cltab(ncltab,nmd), dtcltab(ncltab,nmd)
      REAL dtlog,tmax,tmin
c
c CO rotational cooling
       integer nTco
       parameter (nTco = 1996)

       integer ncdco
       parameter (ncdco = 46)

       REAL co_temptab(nTco), co_colntab(ncdco)

       REAL co_L0(nTco), dTco_L0(nTco)
       REAL co_lte(ncdco,nTco), co_n05(ncdco,nTco), co_alp(ncdco,nTco)
       REAL dTco_lte(ncdco,nTco), dTco_n05(ncdco,nTco), 
     $      dTco_alp(ncdco,nTco)

c CO vibrational cooling
       integer nTco_vib
       parameter (nTco_vib = 3901)

       integer ncdco_vib
       parameter (ncdco_vib = 61)

       REAL co_vib_temptab(nTco_vib), co_vib_colntab(ncdco_vib)
       REAL co_vib_LTE_final(ncdco_vib, nTco_vib)
       REAL dTco_vib_LTE(ncdco_vib, nTco_vib)

       common /co_data/ co_temptab, co_colntab, co_L0, dTco_L0,
     $                  co_lte, co_n05, co_alp, dTco_lte, dTco_n05, 
     $                  dTco_alp, co_vib_temptab, co_vib_colntab,
     $                  co_vib_LTE_final, dTco_vib_LTE

c H2O rotational cooling
       integer nTh2o
       parameter (nTh2o = 3991)

       integer ncdh2o
       parameter (ncdh2o = 91)

       REAL h2o_temptab(nTh2o), h2o_colntab(ncdh2o)

       REAL h2o_L0_ortho(nTh2o), dTh2o_L0_ortho(nTh2o),
     $      h2o_L0_para(nTh2o),  dTh2o_L0_para(nTh2o)

       REAL h2o_LTE_ortho(ncdh2o,nTh2o), 
     $      h2o_n05_ortho(ncdh2o,nTh2o), 
     $      h2o_alp_ortho(ncdh2o,nTh2o),
     $      h2o_LTE_para(ncdh2o,nTh2o), 
     $      h2o_n05_para(ncdh2o,nTh2o), 
     $      h2o_alp_para(ncdh2o,nTh2o)

       REAL dTh2o_LTE_ortho(ncdh2o,nTh2o), 
     $      dTh2o_n05_ortho(ncdh2o,nTh2o), 
     $      dTh2o_alp_ortho(ncdh2o,nTh2o),
     $      dTh2o_LTE_para(ncdh2o,nTh2o), 
     $      dTh2o_n05_para(ncdh2o,nTh2o), 
     $      dTh2o_alp_para(ncdh2o,nTh2o)

c H2O vibrational cooling
       integer nTh2o_vib
       parameter (nTh2o_vib = 3901)

       integer ncdh2o_vib
       parameter (ncdh2o_vib = 61)

       REAL h2o_vib_temptab(nTh2o_vib), h2o_vib_colntab(ncdh2o_vib)
       REAL h2o_vib_LTE_final(ncdh2o_vib, nTh2o_vib)
       REAL dTh2o_vib_LTE(ncdh2o_vib, nTh2o_vib)

       common /h2o_data/ h2o_temptab, h2o_colntab, h2o_L0_ortho,
     $                   dTh2o_L0_ortho, h2o_L0_para, dTh2o_L0_para,
     $                   h2o_LTE_ortho, h2o_n05_ortho, h2o_alp_ortho,
     $                   h2o_LTE_para, h2o_n05_para, h2o_alp_para,
     $                   dTh2o_LTE_ortho, dTh2o_n05_ortho,
     $                   dTh2o_alp_ortho, dTh2o_LTE_para,
     $                   dTh2o_n05_para, dTh2o_alp_para,
     $                   h2o_vib_temptab, h2o_vib_colntab,
     $                   h2o_vib_LTE_final, dTh2o_vib_LTE
c
c These variables are initialized during problem setup
c 
      REAL deff, abundc, abundo, abundsi, abundD, abundM, 
     $     abundN, G0, phi_pah, tdust, dust_to_gas_ratio, lwstartz,
     $     AV_conversion_factor, cosmic_ray_ion_rate, redshift, 
     $     AV_ext, h2_form_ex, h2_form_kin, Z_atom, G0_curr
      integer iphoto, iflag_mn, iflag_ad, iflag_atom, 
     $        iflag_3bh2a, iflag_3bh2b, iflag_h3pra,
     $        id_current, index_current, isrf_option, lwtype,
     $        iflag_h2_opacity

      common /coolr/ temptab, cltab, chtab, dtcltab, dtchtab, 
     $               crtab, crphot, lwstartz, G0_curr, 
     $               phtab, cst, dtlog, tdust, tmax, tmin, 
     $               deff, abundc, abundo, abundsi, abundD, 
     $               abundM, abundN, G0, f_rsc, phi_pah, 
     $               dust_to_gas_ratio, AV_conversion_factor,
     $               cosmic_ray_ion_rate, redshift, AV_ext,
     $               h2_form_ex, h2_form_kin, Z_atom

      common /cooli/ iphoto, iflag_mn, iflag_ad, iflag_atom
     $,              iflag_3bh2a, iflag_3bh2b, iflag_h3pra
     $,              id_current, index_current, isrf_option
     $,              lwtype, iflag_h2_opacity


      REAL thermal_rates(SGCHEM_NUM_THERMAL_RATES)
      common /thermal_info/ thermal_rates

#ifdef DEBUG_COOLING_RATES
      REAL radiative_rates(nrates), chemical_rates(nrates_chem)
      REAL newdt
      common /rate_block/ radiative_rates, chemical_rates, newdt
#endif

#if CHEMISTRYNETWORK == 1
      integer no_dchem
      common /dchem_block/ no_dchem
#endif

#ifdef TREE_RAD
#define NPIX  12*NSIDE*NSIDE
#else
#define NPIX  1
#endif

#define DTCOOL_SCALE_FACTOR 0.1

#if defined(SGCHEM_VARIABLE_CRION) && defined(COSMIC_RAYS)

c CR ionization rate at solar circle
      REAL CRionrate_solcirc
      parameter(CRionrate_solcirc = 3e-17)

c CR energy density at solar circle in eV cm^-3 (Webber 1998)
      REAL CRenden_solcirc
      parameter(CRenden_solcirc = 1.8)

#endif

#endif /* SGCHEM */
