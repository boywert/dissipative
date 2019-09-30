#include "../../build/arepoconfig.h"
#if CHEMISTRYNETWORK == 15
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "non_eq.h"
      integer neq
c
      REAL time, y(nspec), ydot(nspec), rpar(nrpar)
      integer ipar(nipar)
c
      integer itemp, i
      REAL abundances(nabn), rates(nrates), 
     $     rates_chem(nrates_chem)
c
      REAL abe     , abhp    , abHI    , abh2    ,
     $     abcp    , abco    , abheI   , abh3p   ,
     $     abcI    , aboI    , abchx   , abohx   ,
     $     abmp    , abMI    , abhep   , abhcop  ,
     $     ab13co  , ab13cp  , ab13chx , abh13cop,
     $     ab13cI
c
      REAL ylam    , temp    , NH      ,
     $     yn      , dl      , divv    , yntot   , gamma , 
     $     energy  , dtemp   , ncrinv  , h2var0  , h2var1, 
     $     h2var2  , ykdh    , ykdh2   , ykdis   ,
     $     ylam_chem, ekn    , f_ex    , f_tot   , ynsq
c
      REAL rho, cs2, L_jeans
c
      REAL fshield_H2, fshield_CO, fshield_dust, fshield_dust_mol,
     $     chi_mean, AV_mean
c
      REAL G_dust
c
      REAL ph2, ph10, ph48, ph49, ph50, ph51
#ifdef _13CO_
      REAL ph55  ! 13CO
#endif
c
      REAL fh3p
c
      REAL num, denom
c
      REAL Rdust, RH2
c
      REAL phi     , h_gr    , hgrvar1 , hgrvar2
c
      REAL cr1     , cr2     , cr3     , cr6     , cp1  ,
     $     cp11    , cr13
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch9     , ch11    , ch13 , 
     $     ch17    , ch25    , ch34    , ch35    , ch40 ,
     $     ch87    , ch92    , ch105   , ch106   , ch146,
     $     ch147   , ch153   , ch154
c
      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch9   , dtch11  , dtch13 , 
     $     dtch17  , dtch25  , dtch34  , dtch35  , dtch40 ,
     $     dtch87  , dtch92  , dtch105 , dtch106 , dtch146,
     $     dtch147 , dtch153 , dtch154
c
      REAL ch87a, ch87b, ch106a
c
      REAL cdex, rdex, fpheat
c
      REAL energy_dummy
c
      REAL pi
      parameter (pi = 3.1415927d0)
c
c Unpack parameters
c
      yn   = rpar(1)
      dl   = rpar(2)
      divv = rpar(3)
      fshield_H2   = rpar(4)
      fshield_CO   = rpar(5)
      fshield_dust = rpar(6)
      chi_mean     = rpar(7)
      fshield_dust_mol = rpar(8)
      AV_mean      = rpar(9)
c
      NH   = 0.5d0 * yn * dl
      ynsq = yn * yn
c
c Chemical species:
c
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
c Return with zeroed ydot array if input abundances invalid
        do I = 1, nspec
          ydot(I) = 0d0
        enddo
        return
      endif
c
c Tracked, non-equilibrium species
c
      abh2   = max(0d0, y(ih2))
      abhp   = max(0d0, y(ihp))
      abcp   = max(0d0, y(icp))
      abchx  = max(0d0, y(ichx))
      abohx  = max(0d0, y(iohx))
      abco   = max(0d0, y(ico))
      abhcop = max(0d0, y(ihcop))
      abhep  = max(0d0, y(ihep))
      abmp   = max(0d0, y(imp))
#ifdef _13CO_
      ab13cp   = max(0d0, y(i13cp))
      ab13co   = max(0d0, y(i13co))
      ab13chx  = max(0d0, y(i13chx))
      abh13cop = max(0d0, y(ih13cop))
#else
      ab13cp   = 0d0
      ab13co   = 0d0
      ab13chx  = 0d0
      abh13cop = 0d0
#endif
c
c Non-equilibrium species computed using conservation laws
c
      abe   = max(abhp + abcp + abhcop + abhep + abmp + ab13cp 
     $      + abh13cop, 0d0)
      abHI  = max(1d0 - 2d0 * abh2 - abhp - abchx - abohx - ab13chx, 
     $        0d0)
      abheI = max(abhe - abhep, 0d0)
      abcI  = max(abundc - abcp - abchx - abco - abhcop, 0d0)
      aboI  = max(abundo - abohx - abco, 0d0)
      abmI  = max(abundM - abmp, 0d0)
#ifdef _13CO_
      ab13cI = max(abund13c - ab13cp - ab13chx - ab13co - abh13cop, 0d0)
#else
      ab13cI = 0d0
#endif
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance 
c     molecular species (e.g. H-, H2+) and metals]
c
      yntot = (1d0 + abhe - abh2 + abe) * yn
      ekn   = energy / (yntot * kboltz)
c For time being, assume gamma = 5/3
      gamma = 5d0 / 3d0
      temp = (2d0 / 3d0) * ekn
c
c Bounds checking
c
c Temp:
c
      if (temp .le. 0d0) then
        itemp   = 1
        dtemp   = 0d0
        temp    = tmin
        ipar(1) = 1
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = 0d0
      else
        itemp = int(dlog10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in rate_eq.F', itemp, temp, energy
          ABORT(1)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
c Read in coefficients
c
      dtch1  = dtchtab(1, itemp)
      dtch2  = dtchtab(2, itemp)
      dtch3  = dtchtab(3, itemp)
      dtch4  = dtchtab(4, itemp)
      dtch5  = dtchtab(5, itemp)
      dtch6  = dtchtab(6, itemp)
      dtch9  = dtchtab(9, itemp)
      dtch11 = dtchtab(11, itemp)
      dtch13 = dtchtab(13, itemp)
      dtch17 = dtchtab(17, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
      dtch40 = dtchtab(40, itemp)
      dtch87 = dtchtab(87, itemp)
      dtch92 = dtchtab(92, itemp)
      dtch105 = dtchtab(105, itemp)
      dtch106 = dtchtab(106, itemp)
      dtch146 = dtchtab(146, itemp)
      dtch147 = dtchtab(147, itemp)
      dtch153 = dtchtab(153, itemp)
      dtch154 = dtchtab(154, itemp)
c      
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch9  = chtab(9, itemp) + dtemp * dtch9
      ch11 = chtab(11, itemp) + dtemp * dtch11
      ch13 = chtab(13, itemp) + dtemp * dtch13
      ch17 = chtab(17, itemp) + dtemp * dtch17
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
      ch40 = chtab(40, itemp) + dtemp * dtch40
      ch87 = chtab(87, itemp) + dtemp * dtch87
      ch92 = chtab(92, itemp) + dtemp * dtch92
      ch105 = chtab(105, itemp) + dtemp * dtch105
      ch106 = chtab(106, itemp) + dtemp * dtch106
      ch146 = chtab(146, itemp) + dtemp * dtch146
      ch147 = chtab(147, itemp) + dtemp * dtch147
      ch153 = chtab(153, itemp) + dtemp * dtch153
      ch154 = chtab(154, itemp) + dtemp * dtch154
c
c H3+ dissociative recombination: branching fractions from UMIST06,
c based on MAC04
c
      ch106a = (2.34d0 / 6.7d0) * ch106   ! H2 + H
c
c HCO+ dissociative recombination: rate and branching ratios from
c Geppert et al (2005)
c
      ch87a = 0.92d0 * ch87       ! CO + H
      ch87b = 0.08d0 * ch87       ! OH + C
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      if (temp .lt. 5d2) then
c At this temperature, H2 collisional dissociation should be negligible                                                                                       
        ykdh  = 0d0
        ykdh2 = 0d0
      else
        dtch1  = dtchtab(1, itemp)
        dtch2  = dtchtab(2, itemp)
        dtch3  = dtchtab(3, itemp)
        dtch4  = dtchtab(4, itemp)
        dtch5  = dtchtab(5, itemp)
        dtch6  = dtchtab(6, itemp)
        ch1  = chtab(1, itemp) + dtemp * dtch1
        ch2  = chtab(2, itemp) + dtemp * dtch2
        ch3  = chtab(3, itemp) + dtemp * dtch3
        ch4  = chtab(4, itemp) + dtemp * dtch4
        ch5  = chtab(5, itemp) + dtemp * dtch5
        ch6  = chtab(6, itemp) + dtemp * dtch6
c                                                                                                                                                             
        ncrinv   = (2d0 * abh2 * (ch6 - ch5) + ch5)
        h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
        h2var1   = ch3**h2var0
        h2var2   = ch4**h2var0
        ykdh     = ch1 * h2var1
        ykdh2    = ch2 * h2var2
      endif
c
c Photodissociation:
c
c Photoelectric effect:
c
      G_dust  = G0 * fshield_dust
c
c Photodissociation:
c
      ykdis = phtab(1)  * G0 * fshield_H2 * fshield_dust_mol
      ph2   = phtab(2)  * G0 * dexp(-3.0d0 * AV_mean)
      ph10  = phtab(10) * G0 * fshield_CO * fshield_dust_mol
#ifdef _13CO_
      ph55  = phtab(55) * G0 * fshield_13CO * fshield_dust_mol
#endif
      ph48  = phtab(48) * G0 * dexp(-1.5d0 * AV_mean)
      ph49  = phtab(49) * G0 * dexp(-1.7d0 * AV_mean)
      ph50  = phtab(50) * G0 * fshield_dust
      ph51  = phtab(51) * G0 * dexp(-1.9d0 * AV_mean)
c
c Rates for recombination on grain surfaces. These rates must be 
c computed here because of their dependence on the electron number
c density. 
c
      if (abe .eq. 0d0) then
c If the fractional ionization is zero, then there won't be any recombination,
c so the value we use for phi doesn't matter too much -- 1d20 is simply an 
c arbitrary large number
c 
        phi = 1d20
      else
        phi = G_dust * sqrt(temp) / (yn * abe)
      endif
c
c H+:
c
      if (phi .lt. 1d-6) then
        h_gr = 1.225d-13 * dust_to_gas_ratio
      else
        hgrvar1  = 8.074d-6 * phi**1.378d0
        hgrvar2  = (1d0 + ch34 * phi**ch35)
        h_gr     = 1.225d-13 * dust_to_gas_ratio /
     $             (1d0 + hgrvar1 * hgrvar2)
      endif
c Exponentially reduce rate at high temperature to account for fact that we
c expect hot gas to be largely dust free. This also prevents us from getting
c an anomalously large grain recombination cooling rate when T is large
c and G_dust is very small or zero
      h_gr = h_gr * dexp(-temp / 3d4)
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c                  
      if (yn .gt. 1d8 / dust_to_gas_ratio**1.5d0) then
        tdust = temp
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, 1)
      else
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, 0)
      endif
c
c Cosmic ray ionization rates: precomputed in cheminmo.
c
c HI:
      cr1 = crtab(1)
      cr2 = crtab(2)
c HeI:
      cr6 = crtab(6)
c
      cr3  = crtab(3)   ! direct CR photoionization of carbon
c CO: CO + CR -> C + O
      cr13 = crtab(13)
      cp1  = crphot(1)  ! CR-induced photoionization of carbon
      cp11 = crphot(11) * dsqrt(temp) * abH2 * dsqrt(abco) ! CO
c
c H3+ equilibrium abundance
c
      if (abh2 .gt. 0d0) then
        fh3p = ch105 * abh2 / (ch105 * abh2 + cst(1) * abHI)
      else
        fh3p = 0d0
      endif
c
      num   = fh3p * ((cr2 / yn) * abh2 + cst(75) * abhep * abh2)
      denom = cst(84) * abcI + cst(86) * aboI + cst(59) * abco 
     $      + ch106 * abe + cst(81) * abmI
     $      + cst(59) * ab13co + cst(84) * ab13cI
c
      abh3p = num / denom
      abe   = abe + abh3p
c
c Protons
c 
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb. (case B)
     $          -  h_gr * abhp  * 1d0           ! Grain H recomb.
     $          + (1d0 - fh3p)  * cst(75) * abhep * abh2
     $          +  ch9  * abhep * abh2) * yn    ! DCT: He+ + H2 
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + (1d0 - fh3p) * cr2 * abh2     ! H2+ + H -> H2 + H+
     $          +  ph50 * abhcop                ! Photodissoc. of HCO+
     $          +  ph50 * abh13cop              ! Photodissoc. of H13CO+
c
c Molecular hydrogen
c
      ydot(ih2) = (-ch25  * abe  * abh2          ! Collisional dissoc.: H2 + e- 
     $          -   ykdh  * abHI * abh2          ! Collisional dissoc.: H2 + H 
     $          -   ykdh2 * abh2 * abh2          ! Collisional dissoc.: H2 + H2
     $          +   Rdust * abHI * 1d0 
     $          - cst(82) * abh2 * abcp
     $          - cst(82) * abh2 * ab13cp
     $          - cst(54) * abcI * abh2
     $          - cst(54) * ab13cI * abh2
     $          + cst(84) * abh3p * abcI
     $          + cst(84) * abh3p * ab13cI
     $          + cst(86) * abh3p * aboI
     $          + cst(59) * abh3p * abco
     $          + cst(59) * abh3p * ab13co
     $          -     ch9 * abhep * abh2
     $          - cst(75) * abhep * abh2 * (1d0 - fh3p)
     $          +  ch106a * abh3p * abe
     $          + cst(81) * abh3p * abmI) * yn
     $          - 2d0 * fh3p * cr2 * abh2        ! Cosmic rays
     $          -   ykdis * abh2                 ! Photodissoc:         H2 + PH
c
c C+
c
      ydot(icp) = ( ch92 * abhep  * abco
     $          -   ch40 * abcp   * abe
     $          -cst(82) * abcp   * abh2
     $          -cst(83) * abcp   * abohx
     $          +  ch153 * ab13cp * abco
     $          -  ch154 * abcp   * ab13co) * yn
     $          +    ph2 * abcI
     $          +    cr3 * abcI
     $          +    cp1 * abcI
c
c CHx
c
      ydot(ichx) = (cst(84) * abcI   * abh3p
     $           +  cst(82) * abcp   * abh2
     $           +  cst(54) * abcI   * abh2
     $           +    ch87b * abhcop * abe
     $           -  cst(85) * aboI   * abchx) * yn
     $           -     ph48 * abchx
c
c OHx
c
      ydot(iohx) = (cst(86) * aboI * abh3p
     $           -    ch147 * abcI * abohx
     $           -    ch147 * ab13cI * abohx
     $           -  cst(83) * abcp * abohx
     $           -  cst(83) * ab13cp * abohx) * yn
     $           -     ph49 * abohx
c
c HCO+
c
      ydot(ihcop) = (cst(59) * abh3p  * abco
     $            +  cst(83) * abcp   * abohx
     $            -     ch87 * abhcop * abe) * yn
     $            -     ph50 * abhcop
c
c CO
c
      ydot(ico) = (cst(85) * abchx  * aboI
     $          +    ch147 * abohx  * abcI
     $          +    ch87a * abhcop * abe
     $          -  cst(59) * abh3p  * abco
     $          -     ch92 * abhep  * abco
     $          -    ch153 * ab13cp * abco
     $          +    ch154 * abcp   * ab13co) * yn
     $          -     ph10 * abco
     $          +     ph50 * abhcop      ! produces H+       
     $          -     cr13 * abco
     $          -     cp11
c
#ifdef _13CO_
c 13C+
c
      ydot(i13cp)  = (ch92 * abhep  * ab13co
     $          -     ch40 * ab13cp * abe
     $          -  cst(82) * ab13cp * abh2
     $          -  cst(83) * ab13cp * abohx
     $          -    ch153 * ab13cp * abco
     $          +    ch154 * abcp   * ab13co) * yn
     $          +      ph2 * ab13cI
     $          +      cr3 * ab13cI
     $          +      cp1 * ab13cI

c
c 13CHx
c
      ydot(i13chx) = (cst(84) * ab13cI   * abh3p
     $           +  cst(82)   * ab13cp   * abh2
     $           +  cst(54)   * ab13cI   * abh2
     $           +    ch87b   * abh13cop * abe
     $           -  cst(85)   * aboI     * ab13chx) * yn
     $           -     ph48   * ab13chx

c
c H13CO+
c
      ydot(ih13cop) = (cst(59) * abh3p    * ab13co
     $              +  cst(83) * ab13cp   * abohx
     $              -     ch87 * abh13cop * abe) * yn
     $              -     ph50 * abh13cop
c
c 13CO
c
      ydot(i13co) = (cst(85) * ab13chx  * aboI
     $            +    ch147 * abohx    * ab13cI
     $            +    ch87a * abh13cop * abe
     $            -  cst(59) * abh3p    * ab13co
     $            -     ch92 * abhep    * ab13co
     $            +    ch153 * ab13cp   * abco
     $            -    ch154 * abcp     * ab13co) * yn
     $            -     ph55 * ab13co
     $            +     ph50 * abh13cop      ! produces H+
     
      if (abco .gt. 1d-12) then
        ydot(i13co) = ydot(i13co) - cp11 * (ab13co / abco)
      endif
#endif
c
c M+
c
      ydot(imp) = (cst(81) * abh3p * abmI
     $          -    ch146 * abmp  * abe) * yn
     $          +     ph51 * abmI
c
c He+
c
      ydot(ihep) =  cr6  * abheI
     $           - ( ch9 * abhep * abh2
     $           +  ch92 * abhep * abco
     $           +  ch92 * abhep * ab13co
     $           +  ch17 * abhep * abe
     $         + cst(75) * abhep * abh2) * yn
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI * abh2 + 
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) * ynsq
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(2) = - 4d-1 * eV * abH2 * yn * ykdis
c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990),
c and adopt a pumping rate 6.94 times larger than the dissociation rate
c (Draine & Bertoldi 1996, table 2). 
c
      cdex = (1.4 * dexp(-1.81d4 / (temp + 1200)) * abH2
     $     +  1.0 * dexp(-1d3 / temp) * abHI) * 1d-12 * dsqrt(temp) * yn
      rdex = 2d-7
      fpheat = cdex / (cdex + rdex)

      rates_chem(3) = - 2d0 * eV * abH2 * yn * ykdis * 6.94d0 * fpheat
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c Note that energy input depends on formation mechanism:
c
c Dust -- 4.48eV 
c
      if (RH2 .gt. 0d0) then
        call compute_h2_heating_fraction(temp, yn, abh2, f_ex)
        f_ex = f_ex * h2_form_ex
      else
        f_ex = 0d0
      endif
      f_tot = h2_form_kin + f_ex
c
      rates_chem(4) = - 4.48d0 * eV * f_tot * RH2
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(5) = 13.6d0 * eV * ch11 * abe * abHI * ynsq
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(6) = kboltz * temp * abhp * (ch13 * abe + h_gr) * ynsq
c
      ylam_chem = rates_chem(1)  + rates_chem(2)  + rates_chem(3)  +
     $            rates_chem(4)  + rates_chem(5)  + rates_chem(6)
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      abundances(1)  = abh2
      abundances(2)  = aboI
      abundances(5)  = abco
      abundances(6)  = abcI + ab13cI
      abundances(7)  = abcp + ab13cp
      abundances(9)  = abundsi
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
      abundances(14) = abheI
      abundances(15) = abhep
      abundances(17) = abhp * aboI
c
      energy_dummy = energy
c
c Compute local Jeans length (for use in cool_func)
c
      rho = (1d0 + 4d0 * abhe) * mh * yn
      cs2 = (gamma * (gamma - 1d0) * energy / rho)
      L_jeans  = dsqrt(pi * cs2 /
     $           (Gn * rho))
c
      call cool_func(temp, yn, dl, divv, G_dust, abundances, ylam,
     $               rates, 0, energy_dummy)
c
      ydot(itmp) = -(ylam + ylam_chem)
c
#ifdef THERMAL_INFO_DUMP
c
      heat_accr   = 0.
      heat_cray    = rates(11)
      heat_photo_elec = rates(12)
      heat_h2_form = rates_chem(4)
      heat_h2_pdiss = rates_chem(2)
      heat_h2_pump  = rates_chem(3)
      cool_h2_line = rates(2)
      cool_cplus   = rates(16)
      cool_oxy     = rates(13)
      cool_co      = rates(7) + rates(21) + rates(22) +
     $               rates(8) + rates(9) + rates(24)
      cool_13co     = rates(9)  + rates(21)
      cool_c18o     = rates(22) + rates(24)
      cool_gas_dust = rates(1)
      cool_c        = rates(14)
#endif
c
#endif /* SGCHEM */
      return      
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              R A T E _ E Q                \\\\\\\\\\
c
c=======================================================================
#endif /* CHEMISTRYNETWORK == 15 */
