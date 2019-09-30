#include "../../build/arepoconfig.h"
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
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
     $     rates_chem(nrates_chem),
     $     y_in(nspec)
c
      REAL abe     , abhp    , abHI    , abh2    ,
     $     abcp    , abco    , abo     , abcos
c
      REAL ylam    , temp    ,
     $     yn      , dl      , divv    , yntot   , gamma , 
     $     energy  , dtemp   , ncrinv  , h2var0  , h2var1, 
     $     h2var2  , ykdh    , ykdh2   , ykdis   , NH    ,
     $     ylam_chem, f_ex   , rho     , cs2     , L_jeans,
     $     f_tot   , RH2     , Rdust   , ekn, energy_dummy
c
      REAL fshield_H2, fshield_CO, fshield_dust, chi_mean, 
     $     fshield_dust_H2
c
      REAL G_dust
c
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      REAL k0, k1, gamma_chx, gamma_co, beta
      REAL kd1, kd2, kd3
#endif
c
      REAL phi     , h_gr    , hgrvar1 , hgrvar2
c
      REAL cr1, cr2, cr10, cr11
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch11    , ch13    , ch25    ,
     $     ch34    , ch35
c
      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch11  , dtch13  , dtch25  ,
     $     dtch34  , dtch35
c
      REAL cdex, rdex, fpheat 
c
c Unpack parameters
c

      yn   = rpar(1)
      dl   = rpar(2)
      divv = rpar(3)
      fshield_H2   = rpar(4)
      if (isnan(fshield_H2)) print*,"rate_eq_simple fshield_H2 is NaN"
      fshield_CO   = rpar(5)
      fshield_dust = rpar(6)
      chi_mean     = rpar(7)
      fshield_dust_H2 = rpar(8)
c
c Chemical species:
c
      do I = 1, nspec
        y_in(I) = y(I)
      enddo
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
      abhp = max(0d0, y(ihp))
      abh2 = max(0d0, y(ih2))
      abco  = 0d0
#if CHEMISTRYNETWORK == 5
      abco = max(0d0, y(ico))
      abcos=0d0
#endif
c
c Non-equilibrium species computed using conservation laws
c
      abHI = max(1d0 - 2d0 * abh2 - abhp, 0d0)
#if CHEMISTRYNETWORK == 5
      abcp = max(0d0, abundc - abco)
      abe  = abhp + abcp + abundsi
c     if (id_current .eq. 1) print*, "Fortran co", abco,abcp,abe
#else
      abcp = abundc
      abe  = abhp + abundc + abundsi
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
c     call lookup_temp_gamma(abh2, ekn, temp, gamma)
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
          print*, 'Fatal error in rate_eq.F', itemp, temp
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
      dtch11 = dtchtab(11, itemp)
      dtch13 = dtchtab(13, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
c
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch11 = chtab(11, itemp) + dtemp * dtch11
      ch13 = chtab(13, itemp) + dtemp * dtch13
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      ncrinv   = (2d0 * abh2 * (ch6 - ch5) + ch5)
      h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
      h2var1   = ch3**h2var0
      h2var2   = ch4**h2var0
      ykdh     = ch1 * h2var1
      ykdh2    = ch2 * h2var2
c
c Photoelectric effect:
c
      G_dust  = G0 * fshield_dust
c
c Photodissociation:
c
      ykdis = phtab(1) * G0 * fshield_H2 * fshield_dust_H2
#if CHEMISTRYNETWORK == 5
      gamma_co = phtab(10) * G0 * fshield_dust * fshield_CO
      gamma_chx = 2.94d-10 * G_dust
c      if (id_current .eq. 1) print*, "Fortran gamma co", gamma_co
#endif
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
c HII:
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
      call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, 0)


c
c Cosmic ray ionization rates: precomputed in cheminmo.
c
c HI:
      cr1  = crtab(1)  
      cr2  = crtab(2)   ! H2 -> H2+ + e-
      cr10 = crtab(10)  ! H2 -> H + H+ + e
      cr11 = crtab(11)  ! H2 -> H + H
c
c Protons
c 
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $          -  h_gr * abhp  * 1d0) * yn     ! Grain H recomb.
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + cr10  * abh2
c
c Molecular hydrogen
c
      ydot(ih2) = (-ch25  * abe  * abh2          ! Collisional dissoc.: H2 + e- 
     $          -   ykdh  * abHI * abh2          ! Collisional dissoc.: H2 + H 
     $          -   ykdh2 * abh2 * abh2          ! Collisional dissoc.: H2 + H2
     $          +   Rdust * abHI * 1d0 ) * yn    ! H2 formation on dust
     $          -   ykdis * abh2                 ! Photodissoc:         H2 + PH
     $      - 0.5d0 * cr2 * abh2   ! Fudge factor of two gives correct HI production
     $          -   cr10  * abh2
     $          -   cr11  * abh2
c
#if CHEMISTRYNETWORK == 5
c
c CO -- from NL97
c
      k0 = 5d-16
      k1 = 5d-10
#if CHEMISTRYNETWORK == 5
      abo = max(0.0, abundo - abco)
#endif
#if CHEMISTRYNETWORK == 6
      abo = max(0.0, abundo - abco - abcos)
#endif
      if (abo .eq. 0.0) then
        beta = 0.0
      elseif (gamma_chx .lt. 1d-30) then
        beta = 1.0
      else
        beta = k1 * abo / (k1 * abo + 
     $       gamma_chx / yn)
      endif
c
      ydot(ico) = k0 * abcp * abh2 * beta * yn 
     $          - gamma_co * abco
c      if (id_current .eq. 1) print*, "Fortran ydot", ydot(ico),ico
c
#if CHEMISTRYNETWORK == 6
c
      kd1 = 1.1d-17 * dsqrt(temp)  ! Accretion
      kd2 = 1.04d12 * dexp(-1.21d3 / tdust) ! Thermal evaporation
      kd3 = 9.8d-15 ! Cosmic ray desorption
c
      ydot(icos) = kd1 * abco * yn - (kd2 + kd3) * abcos
      ydot(ico) = ydot(ico) - ydot(icos)
c
#endif
#endif
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI * abh2 + 
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) * yn**2
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(2) = - 4d-1 * eV * abH2 * yn * ykdis
c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990),
c and that 85% of excitations end in pumping, so that the pumping
c rate is 0.85/0.15 ~ 5.67 times the dissociation rate
c (see Draine & Bertoldi 1996, table 2).
c
      cdex = (1.4 * dexp(-1.81d4 / (temp + 1200)) * abH2
     $     +  1.0 * dexp(-1d3 / temp) * abHI) * 1d-12 * dsqrt(temp) * yn
      rdex = 2d-7
      fpheat = cdex / (cdex + rdex)

      rates_chem(3) = - 2d0 * eV * abH2 * yn * ykdis * 5.67d0 * fpheat
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c We assume that a fraction f_ex of the energy goes into rotational &
c vibrational excitation of the H2, and a further fraction f_kin goes
c directly into kinetic energy (i.e. heat). The remaining energy goes
c directly into the dust 
c
c Dust -- 4.48eV 
c
      if (RH2 .gt. 0d0) then
c        call compute_h2_heating_fraction(temp, yn, abh2, f_ex)
c        f_ex = f_ex * h2_form_ex
        f_ex = 0d0
      else
        f_ex = 0d0
      endif
      f_tot = h2_form_kin + f_ex
c
      rates_chem(4) = - 4.48d0 * eV * f_tot * RH2
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(5) = 13.6d0 * eV * ch11 * abe * abHI * yn**2
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(6) = kboltz * temp * abhp * (ch13 * abe + h_gr) * yn**2
c
      ylam_chem = rates_chem(1)  + rates_chem(2)  + rates_chem(3)  +
     $            rates_chem(4)  + rates_chem(5)  + rates_chem(6)

c      if (id_current .eq. 1) print*, "Fortran calculated rates"
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      abundances(1)  = abh2

#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      abundances(2)  = max(0d0, abundo - abco - abcos)
      abundances(5)  = abco
      abundances(7)  = abcp
#else
      abundances(2)  = abundo
      abundances(7)  = abundc
#endif
      abundances(9)  = abundsi
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
      abundances(14) = abhe
      abundances(17) = abundo * abhp
c
c Compute local Jeans length (for use in cool_func)
c
      call compute_Jeans_length(energy, yn, gamma, L_jeans)
c
      energy_dummy = energy

      call cool_func(temp, yn, L_jeans, divv, G_dust, abundances, ylam,
     $               rates, 0, energy_dummy)

c            
c Store cooling rates for output
c
#ifdef SGCHEM_DUMP_THERMAL_RATES
c H2 collisional dissociation cooling: 
      thermal_rates(13) = rates_chem(1)
c H2 photodissociation heating: 
      thermal_rates(14) = rates_chem(2)
c UV pumping of H2: 
      thermal_rates(15) = rates_chem(3)
c H2 formation heating: 
      thermal_rates(16) = rates_chem(4)
c HI collisional ionization cooling:
      thermal_rates(17) = rates_chem(5)
c HII recombination cooling:
      thermal_rates(18) = rates_chem(6) 
#endif

      ydot(itmp) = -(ylam + ylam_chem)

c      if (id_current .eq. 1) print*, "Fortran: rate_eq  ydot",
c     &     ydot(itmp),itmp
c
#if 0
c Useful for debugging if we can isolate problems to a single
c particle. 
      if (id_current .eq. 6658924) then
        print*, "Particle number", id_current
        print*, '1 yn', yn, temp
        print*, '2 y', y
        print*, '3 ydot', ydot
        print*, '4 ylam', ylam, ylam_chem
        print*, '5G dust', G_dust
        print*, '6 rates_chem', rates_chem
        print*, '7 rates', rates
        print*, '8 tdust', tdust
      endif
#endif
c
      call validate_output(y,ydot,rpar)
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
#endif /* CHEMISTRYNETWORK == 4 || 5 || 6 */
