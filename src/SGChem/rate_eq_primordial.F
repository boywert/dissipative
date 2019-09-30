#include "../../build/arepoconfig.h"
#if CHEMISTRYNETWORK == 1
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
c    Written by S. Glover (2004-2006)
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
      integer niter
      parameter (niter = 25)
      REAL eqtol
      parameter (eqtol = 1d-3)
c
      REAL    dtemp
      integer itemp, I, J, K
c
      REAL abundances(nabn), rates(nrates),
     $     rates_chem(nrates_chem), ylam, ylam_chem,
     $     CMB_temp, ylam_cmb, rates_cmb(nrates)
      REAL temp, yn, dl, divv, energy, gamma, yntot,
     $     ekn

      REAL abhp, abh2, abhd, abdp, abhep, abhepp, abco 
      REAL abhm, abhm_old, abh2p, abh2p_old, abhehp, abhehp_old
c
      REAL num, denom
      REAL comp_equilib
      external comp_equilib

      REAL total_rad_heating, heating_rate, Lstar, distance

      REAL time_to_equilibrium(num_eqb_species)

      REAL abe, abHI, abDI, abheI

      REAL abM

      REAL ph4,  ph5, ph8, ph9, ph35, ph38, ph39,
     $     ph40, ph41, ph42, ph43, ph44, ph45
 
      REAL h2p_pd_low, h2p_pd_lte

      REAL cr1, cr1a, cr2, cr6

      REAL ykdis, NH, ncrinv, h2var0,
     $     h2var2, ynsq, rho, cs2, L_jeans

      REAL fshield_H2

      REAL cdex, rdex, fpheat

      REAL cdh2_h, cdh2_h_dt, cdh2_h2, cdh2_he

      REAL h2_ct_low

      REAL h2p_dr, h2p_cdh, ncrit_h2p_h, ncrit_h2p_e, ncrit_h2p_inv, 
     $     h2p_alpha

      REAL  ch1,  ch2,  ch3,  ch4,  ch5,  ch6,  ch8,  ch9,
     $     ch10, ch11, ch13, ch14, ch15, ch17, ch18, ch20,
     $     ch21, ch22, ch23, ch24, ch25, ch26, ch27, ch28, 
     $     ch29, ch30, ch31, ch32, ch33, ch98, ch99, 
     $     ch100, ch101, ch102, ch103, ch108, ch109, ch110,
     $     ch126, ch127, ch137, ch138, ch139, ch144,
     $     ch153, ch154, ch155, ch156, ch157, ch158, 
     $     ch159, ch160, ch161, ch162, ch163, ch164, ch165, 
     $     ch166, ch167, ch168, ch169, ch170, ch171, ch172, 
     $     ch173

      REAL  dtch1,  dtch2, dtch3,   dtch4,  dtch5,  dtch6,
     $      dtch8,  dtch9, dtch10, dtch11, dtch13, dtch14, 
     $     dtch15, dtch17, dtch18, dtch20, dtch21, dtch22, 
     $     dtch23, dtch24, dtch25, dtch26, dtch27, dtch28, 
     $     dtch29, dtch30, dtch31, dtch32, dtch33,
     $     dtch98, dtch99,  dtch100, dtch101, dtch102, 
     $     dtch103, dtch108, dtch109, dtch110,
     $     dtch126, dtch127, dtch137, dtch138, dtch139,
     $     dtch144, dtch153, dtch154, dtch155,
     $     dtch156, dtch157, dtch158, dtch159, 
     $     dtch160, dtch161, dtch162, dtch163, dtch164, 
     $     dtch165, dtch166, dtch167, dtch168, dtch169, 
     $     dtch170, dtch171, dtch172, dtch173
c
      no_dchem = 0
c
c Unpack parameters
c
      yn         = rpar(1)
      dl         = rpar(2)
      divv       = rpar(3)
      fshield_H2 = rpar(4)
c
      ynsq  = yn**2
c
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
        do I = 1, nspec
          ydot(I) = 0d0
        enddo
c We simply return here because we'll be redoing the step anyway
c (making any further work pointless) and because we may be doing 
c this on account of having negative chemical abundances, which 
c would otherwise cause serious problems later on.
        return
      endif
c
c Here, we call lwbg to calculate the current LWBG rate G0_curr:
c    
      call lwbg
c      print*,'G0_curr = ',G0_curr
c
c From this point on, we know that our tracked abundances cannot be negative
c
      abhp   = y(ihp)
      abh2   = y(ih2)
      abhd   = y(ihd)
      abdp   = y(idp)
      abhep  = y(ihep)
      abhepp = y(ihepp)
c
      abe  = abhp + abdp + abhep + 2d0 * abhepp
c
      abHI = 1d0 - 2d0 * abH2 - abhp - abhd
      abHI = max(0d0, abHI)
c
      abDI  = max(0d0, abundD - abdp - abhd)
      if (yn .gt. 9d7 .and. no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        abDI = max(abDI, atol(ihd))
      endif
c
      abheI = max(0d0, abhe - abhep - abhepp)
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute current adiabatic index and total number density
c
      yntot = (1d0 + abhe + abe - abh2) * yn
      ekn = energy / (kboltz * yntot)
c      call lookup_temp_gamma(abh2, ekn, temp, gamma) !XXX
c
c XXX: for high density runs, need to re-enable this. However,
c at low densities, H2 fraction is small and gamma = 5/3 to
c within a small fraction of a percent
c
      gamma = 5d0/3d0
      temp = (gamma - 1d0) * ekn
c
c Compute local Jeans length (for use in cool_func)
c
      call compute_Jeans_length(energy, yn, gamma, L_jeans)
c
c Compute CMB temperature
c
      CMB_temp = 2.726d0 * (1d0 + redshift)
c
c Bounds checking
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
c Rates:
c
       dtch1  = dtchtab(1, itemp)
       dtch2  = dtchtab(2, itemp)
       dtch3  = dtchtab(3, itemp)
       dtch4  = dtchtab(4, itemp)
       dtch5  = dtchtab(5, itemp)
       dtch6  = dtchtab(6, itemp)
c
       dtch8  = dtchtab(8, itemp)
       dtch9  = dtchtab(9, itemp)
c 
       dtch10 = dtchtab(10,itemp)
       dtch11 = dtchtab(11,itemp) 
c
       dtch13 = dtchtab(13,itemp)
       dtch14 = dtchtab(14,itemp)
       dtch15 = dtchtab(15,itemp)
c
       dtch17 = dtchtab(17,itemp)
       dtch18 = dtchtab(18,itemp)
c
       dtch20 = dtchtab(20,itemp)
       dtch21 = dtchtab(21,itemp)
       dtch22 = dtchtab(22,itemp)
       dtch23 = dtchtab(23,itemp)
       dtch24 = dtchtab(24,itemp)
       dtch25 = dtchtab(25,itemp)
       dtch26 = dtchtab(26,itemp)
       dtch27 = dtchtab(27,itemp)
       dtch28 = dtchtab(28,itemp)
       dtch29 = dtchtab(29,itemp)
       dtch30 = dtchtab(30,itemp)
       dtch31 = dtchtab(31,itemp)
       dtch32 = dtchtab(32,itemp)
       dtch33 = dtchtab(33,itemp)
c
       dtch98 = dtchtab(98,itemp)
       dtch99 = dtchtab(99,itemp)
       dtch100 = dtchtab(100,itemp)
       dtch101 = dtchtab(101,itemp)
       dtch102 = dtchtab(102,itemp)
       dtch103 = dtchtab(103,itemp)
c
       dtch108 = dtchtab(108,itemp)
       dtch109 = dtchtab(109,itemp)
       dtch110 = dtchtab(110,itemp)
c
       dtch126 = dtchtab(126,itemp)
       dtch127 = dtchtab(127,itemp)
c
       dtch137 = dtchtab(137,itemp)
       dtch138 = dtchtab(138,itemp)
       dtch139 = dtchtab(139,itemp)
c
       dtch144 = dtchtab(144,itemp)
c
       dtch153 = dtchtab(153,itemp)
       dtch154 = dtchtab(154,itemp)
       dtch155 = dtchtab(155,itemp)
       dtch156 = dtchtab(156,itemp)
       dtch157 = dtchtab(157,itemp)
       dtch158 = dtchtab(158,itemp)
       dtch159 = dtchtab(159,itemp)
       dtch160 = dtchtab(160,itemp)
       dtch161 = dtchtab(161,itemp)
       dtch162 = dtchtab(162,itemp)
       dtch163 = dtchtab(163,itemp)
       dtch164 = dtchtab(164,itemp)
       dtch165 = dtchtab(165,itemp)
       dtch166 = dtchtab(166,itemp)
       dtch167 = dtchtab(167,itemp)
       dtch168 = dtchtab(168,itemp)
       dtch169 = dtchtab(169,itemp)
       dtch170 = dtchtab(170,itemp)
       dtch171 = dtchtab(171,itemp)
       dtch172 = dtchtab(172,itemp)
       dtch173 = dtchtab(173,itemp)
c
       ch1  = chtab(1, itemp) + dtemp * dtch1
       ch2  = chtab(2, itemp) + dtemp * dtch2
       ch3  = chtab(3, itemp) + dtemp * dtch3
       ch4  = chtab(4, itemp) + dtemp * dtch4
       ch5  = chtab(5, itemp) + dtemp * dtch5
       ch6  = chtab(6, itemp) + dtemp * dtch6
c       
       ch8  = chtab(8, itemp) + dtemp * dtch8
       ch9  = chtab(9, itemp) + dtemp * dtch9
c
       ch10 = chtab(10,itemp) + dtemp * dtch10 
       ch11 = chtab(11,itemp) + dtemp * dtch11
c
       ch13 = chtab(13,itemp) + dtemp * dtch13
c
       ch14 = chtab(14,itemp) + dtemp * dtch14
       ch15 = chtab(15,itemp) + dtemp * dtch15
c
       ch17 = chtab(17,itemp) + dtemp * dtch17
       ch18 = chtab(18,itemp) + dtemp * dtch18
c
       ch20 = chtab(20,itemp) + dtemp * dtch20
       ch21 = chtab(21,itemp) + dtemp * dtch21
       ch22 = chtab(22,itemp) + dtemp * dtch22
       ch23 = chtab(23,itemp) + dtemp * dtch23
       ch24 = chtab(24,itemp) + dtemp * dtch24
       ch25 = chtab(25,itemp) + dtemp * dtch25
       ch26 = chtab(26,itemp) + dtemp * dtch26
       ch27 = chtab(27,itemp) + dtemp * dtch27
       ch28 = chtab(28,itemp) + dtemp * dtch28
       ch29 = chtab(29,itemp) + dtemp * dtch29
       ch30 = chtab(30,itemp) + dtemp * dtch30
       ch31 = chtab(31,itemp) + dtemp * dtch31
       ch32 = chtab(32,itemp) + dtemp * dtch32
       ch33 = chtab(33,itemp) + dtemp * dtch33
c
       ch98 = chtab(98,itemp) + dtemp * dtch98
       ch99 = chtab(99,itemp) + dtemp * dtch99
       ch100 = chtab(100,itemp) + dtemp * dtch100
       ch101 = chtab(101,itemp) + dtemp * dtch101
       ch102 = chtab(102,itemp) + dtemp * dtch102
       ch103 = chtab(103,itemp) + dtemp * dtch103
c
       ch108 = chtab(108,itemp) + dtemp * dtch108
       ch109 = chtab(109,itemp) + dtemp * dtch109
       ch110 = chtab(110,itemp) + dtemp * dtch110
c
       ch126 = chtab(126,itemp) + dtemp * dtch126
       ch127 = chtab(127,itemp) + dtemp * dtch127
c
       ch137 = chtab(137,itemp) + dtemp * dtch137
       ch138 = chtab(138,itemp) + dtemp * dtch138
       ch139 = chtab(139,itemp) + dtemp * dtch139
c
       ch144 = chtab(144,itemp) + dtemp * dtch144
c
       ch153 = chtab(153,itemp) + dtemp * dtch153
       ch154 = chtab(154,itemp) + dtemp * dtch154
       ch155 = chtab(155,itemp) + dtemp * dtch155
       ch156 = chtab(156,itemp) + dtemp * dtch156
       ch157 = chtab(157,itemp) + dtemp * dtch157
       ch158 = chtab(158,itemp) + dtemp * dtch158
       ch159 = chtab(159,itemp) + dtemp * dtch159
       ch160 = chtab(160,itemp) + dtemp * dtch160
       ch161 = chtab(161,itemp) + dtemp * dtch161
       ch162 = chtab(162,itemp) + dtemp * dtch162
       ch163 = chtab(163,itemp) + dtemp * dtch163
       ch164 = chtab(164,itemp) + dtemp * dtch164
       ch165 = chtab(165,itemp) + dtemp * dtch165
       ch166 = chtab(166,itemp) + dtemp * dtch166
       ch167 = chtab(167,itemp) + dtemp * dtch167
       ch168 = chtab(168,itemp) + dtemp * dtch168
       ch169 = chtab(169,itemp) + dtemp * dtch169
       ch170 = chtab(170,itemp) + dtemp * dtch170
       ch171 = chtab(171,itemp) + dtemp * dtch171
       ch172 = chtab(172,itemp) + dtemp * dtch172
       ch173 = chtab(173,itemp) + dtemp * dtch173
c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
c CID rate from MSM96
c
      cdh2_h = (ch161 - (ch161 - ch163)
     $       / (1d0 + (yn * abhI / ch165)**ch166))
     $       + (ch162 - (ch162 - ch164)
     $       / (1d0 + (yn * abhI / (1d1**(-1.164408e0)
     $       * ch165))**ch166))
      cdh2_h = 1d1**cdh2_h
c
c DT rate from MSM96
c
      cdh2_h_dt = (ch167 - (ch167 - ch169)
     $          / (1d0 + (yn * abhI / ch171)**ch172))
     $          + (ch168 - (ch168 - ch170)
     $          / (1d0 + (yn * abhI / (1d1**(-1.174242e0)
     $          * ch171))**ch172))
      cdh2_h_dt = 1d1**cdh2_h_dt
c
c Total dissociation rate
c
      cdh2_h = cdh2_h + cdh2_h_dt
c
      ncrinv   = ch5 + abheI * ch139
      h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
      h2var2   = ch4**h2var0
      cdh2_h2  = ch2 * h2var2
      if (ch138 .lt. 1d-100) then
        cdh2_he = 1d-100
      else
        cdh2_he = ch138 * (ch137 / ch138)**h2var0
      endif
c
c H2 destruction by charge transfer with H+
c
      if (dabs(ch24) .lt. 1d-100) then
        ch24 = 1d-100 !sign?!
      endif
      if(ch173 .lt. 1d-100) then
        ch173 = 1d-100
      endif
        h2_ct_low = ch24
        ch24 = ch173 * (h2_ct_low / ch173)**h2var0
c      if (ch24 .ne. ch24) then
c        print*,'RATES_EQ ch24 NaN:, ch173:',ch173,'h2_ct_low:',h2_ct_low
c        print*,'RATES_EQ ch24 NaN: h2var0',h2var0
c      endif
c
c Compute H2+ critical density
c
      ncrit_h2p_h = 4d2 * (1d4 / temp)
      ncrit_h2p_e = 5d1
      if (abHI .gt. 0d0 .or. abe .gt. 0d0) then
        ncrit_h2p_inv = (abHI / ncrit_h2p_h) + (abe / ncrit_h2p_e)
      else
        ncrit_h2p_inv = 1d0 / 4d2
      endif
      h2p_alpha = 1d0 / (1d0 + yn * ncrit_h2p_inv)
c
c H2+ dissociative recombination
c
      if (ch157 .lt. 1d-100) then
        ch157 = 1d-100
      endif
      h2p_dr   = ch157 * (ch30 / ch157)**h2p_alpha
c
c H2+ collisional dissociation by H 
c
      if (ch160 .lt. 1d-100) then
        ch160 = 1d-100
      endif
      h2p_cdh  = ch160 * (ch159 / ch160)**h2p_alpha
c
c XXX: shielding by HI currently neglected; only significant when NHI very large
c
      ykdis = phtab(1) * G0_curr * fshield_H2
      ph4   = phtab(4) * G0_curr
c At high redshift, we also need to account for dissociation of H- by the CMB.
c There are two contributions: one from the thermal spectrum and a second from
c the non-thermal photons produced during recombination (see Hirata & Padmanabhan, 2006).
c
c Thermal (fit from Galli & Palla 1998):
      ph4 = ph4 + 0.11d0 * CMB_temp**2.13d0 * exp(-8.823d3 / CMB_temp)
c Non-thermal (fit from Coppola et al 2011):
      ph4 = ph4 + 8.0d-8 * CMB_temp**1.3d0 * exp(-2.3d3 / CMB_temp)
c
c H2+ photodissociation rate depends on density, temperature (at high n)
c This prescription assumes a T5 spectrum
      if (G0_curr .eq. 0) then
        ph5 = 0
      else
        h2p_pd_low = phtab(5) * G0_curr
        h2p_pd_lte = 2.2d-11  * G0_curr * temp**(-0.22d0)
c
        if (temp .lt. 3.15d3 .and. h2p_pd_lte .gt. h2p_pd_low) then
          h2p_pd_lte = h2p_pd_low
        endif
c
        if (h2p_pd_lte .lt. 1d-100) then
          h2p_pd_lte = 1d-100
        endif
c
        ph5 = h2p_pd_lte * (h2p_pd_low / h2p_pd_lte)**h2p_alpha
      endif
c
      ph8   = phtab(8) * G0_curr
      ph9   = phtab(9) * G0_curr
c
c XXX: HD self-shielding and shielding by H2 currently neglected
c
      ph35  = phtab(35) * G0_curr
c
c Photoionization: currently not treated
c
      ph38  = phtab(38) * G0_curr  
      ph39  = phtab(39) * G0_curr
      ph40  = phtab(40) * G0_curr
      ph41  = phtab(41) * G0_curr
      ph42  = phtab(42) * G0_curr
      ph43  = phtab(43) * G0_curr
      ph44  = phtab(44) * G0_curr
      ph45  = phtab(45) * G0_curr

      cr1  = crtab(1)  ! HI
      cr1a = crtab(1)  ! DI, assumed same as HI
      cr2  = crtab(2)  ! H2
      cr6  = crtab(6)  ! HeI
c
c Estimate whether we're optically thick to ionizing radiation;
c if we're not, then we need to set ch18 to zero.
c
      NH = 0.5d0 * dl * yn
c
      if (2d0 * NH * abhI .lt. 1.6d17) then
        ch18 = 0d0
      endif
c
c Equilibrium species
c
c Since the equilibrium abundances are coupled, we have to iterate to find 
c accpetable solutions. We iterate niter times (or until all abundances have
c converged).
c
      abhm_old   = 0d0
      abh2p_old  = 0d0
      abhehp_old = 0d0
c
      abhm   = 0d0
      abh2p  = 0d0
      abhehp = 0d0
#ifdef COLLECT_EQB_DATA
      do I = 1, num_eqb_species
        time_to_equilibrium(I) = 0d0
      enddo
#endif
c
      do I = 1, niter
c
c H- :
c
        num   = ch21 * abHI * abe + ch126 * abh2 * abe
        denom = ch22 * abHI + ch26 * abe + ch27 * abHI + ch28 * abhp
     $        + ch29 * abhp + (ph4 / yn) + ch8 * abhep
     $        + ch31 * abheI + 2d0 * ch127 * abh2p
c
        if(num .ne. num .or. denom .ne. denom) then
          print*,'RATE_EQ1 num=',num,'  denom=',denom
        endif

        abhm = comp_equilib(num, denom, 1d0, 1, yn, time_to_equilibrium)
        if (abhm .lt. 0d0) then
           PRINT*,"H- abundance error"
           PRINT*,"Particle id ",id_current
           PRINT*,"T=",temp," no. density H", yn
           PRINT*,"abHI=",abHI," abe=",abe,"abh2= ",abh2
           PRINT*,"abhp=",abhp," abhep=",abhep," abheI=",abheI,
     &          "abh2p= ",abh2p
          ABORT(1000)
        endif
c
        abe  = max(0d0, abe  - (abhm - abhm_old))
        abHI = max(0d0, abHI - (abhm - abhm_old))
c
c H2+ :
c
        abM = abHI + abh2 / dsqrt(2d0) + abheI / 2d0
c
        num   = ch23 * abHI * abhp + ch24 * abh2 * abhp
     $        + ch29 * abhp * abhm + ch156 * abhehp * abHI
     $        + (cr2 * abh2 / yn)  + cst(75) * abhep * abh2 
     $        + (ph41 / yn) * abh2
c
        denom = cst(1) * abHI + h2p_dr * abe + (ph5 / yn)
     $        + 2d0 * ch127 * abhm + ch154 * abheI
     $        + h2p_cdh * abM
c
        if(num .ne. num .or. denom .ne. denom) then
          print*,'RATE_EQ2 num=',num,'  denom=',denom
          print*,'RATE_EQ denom1:',cst(1), abHI, h2p_dr, abe, ph5, yn
          print*,'RATE_EQ denom2:',ch127, abhm, ch154, abheI
          print*,'RATE_EQ denom3:',h2p_cdh, abM
          print*,'RATE_EQ denom4:',ph41
        endif
        abh2p = comp_equilib(num, denom, 0.5d0, 2, yn, 
     $                       time_to_equilibrium)
        if (abh2p .lt. 0d0) then
           PRINT*,"H2+ abundance error"
           PRINT*,"Particle id ",id_current
           PRINT*,"T=",temp," no. density H", yn
           PRINT*,"abHI=",abHI," abe=",abe,"abh2= ",abh2
           PRINT*,"abhp=",abhp," abhep=",abhep," abheI=",abheI,
     &          "abh2p= ",abh2p
           abh2p = 0.0 
        endif
c
        abe  = max(0d0, abe  + (abh2p - abh2p_old))
        abHI = max(0d0, abHI - 2d0*(abh2p - abh2p_old))
c
c HeH+
c
        num = ch155 * abheI * abhp
     $      + ch154 * abheI * abh2p

        denom = ch156 * abHI + ch153 * abe
        if(num .ne. num .or. denom .ne. denom) then
          print*,'RATE_EQ3 num=',num,'  denom=',denom
        endif
        abhehp = comp_equilib(num, denom, abhe, 3, yn,
     $                        time_to_equilibrium)
c
        abe   = max(0d0, abe + (abhehp - abhehp_old))
        abHI  = max(0d0, abhI - (abhehp - abhehp_old))
        abheI = max(0d0, abheI - (abhehp - abhehp_old))
c
c Check for convergence:
c
        if (I .ge. 2) then
          if ((abs(abhm_old - abhm)   .lt. eqtol * abhm)  .and.
     $        (abs(abh2p_old - abh2p) .lt. eqtol * abh2p) .and.
     $        (abs(abhehp_old - abhehp) .lt. eqtol * abhehp)) then
            goto 100
          endif
        endif
        abhm_old  = abhm
        abh2p_old = abh2p
      enddo
c
100   continue
c
c Protons
c
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          + ch144 * abHI  * abHI          ! Collisional ion: H + H
     $          + ch158 * abheI * abHI          ! Collisional ion: H + He
     $         + cst(1) * abh2p * abHI          ! Charge transfer: H  + H2+
     $           - ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $           - ch23 * abHI  * abhp          ! H2+ formation:   H  + H+
     $           - ch24 * abh2  * abhp          ! Charge transfer: H2 + H+
     $           - ch28 * abhp  * abhm          ! Neutralization:  H- + H+
     $           - ch29 * abhp  * abhm          ! H2+ formation:   H- + H+
     $           + ch18 * abhep * abe           ! Ion. from He+ rec: He+ + e
     $           + ch32 * abhep * abHI          ! Charge transfer: He+ + H
     $           - ch33 * abheI * abhp          ! Charge transfer: He  + H+
     $          + h2p_cdh * abh2p * abM         ! Collisional dissoc. of H2+
     $          + ch9 * abhep * abh2) * yn    ! DiCT: He+ + H2 -> He + H + H+
     $          +  ph5  * abh2p                 ! Photodissoc:    H2+ + PH
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + ph38  * abhI                  ! Photoion: H + PH
c
      if (no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        ydot(ihp) = ydot(ihp) 
     $            + (ch99 * abdp  * abHI          ! Charge transfer: H + D+
     $            -  ch98 * abhp  * abDI          ! Charge transfer: D + H+
     $            + ch100 * abdp  * abh2          ! HD form:         D+ + H2
     $            - ch101 * abhd  * abhp) * yn    ! HD dest:         HD + H+
      endif
c
c Molecular hydrogen
c
      ydot(ih2) = (ch108  * abHI**3           ! 3-body H2 form: 3H -> H2 + H 
     $           + ch109  * abHI**2 * abh2    ! "  : 2H + H2 -> H2 + H2
     $           + ch110  * abHI**2 * abheI) * yn**2 ! "  : 2H + He -> H2 + He 
     $          + (ch22   * abHI    * abhm    ! H2 formation:        H- + H
     $          + cst(1)  * abHI    * abh2p   ! Charge transfer:     H  + H2+
     $           -  ch24  * abhp    * abh2    ! Charge transfer:     H2 + H+
     $           -  ch25  * abe     * abh2    ! Collisional dissoc.: H2 + e- 
     $           - cdh2_h * abHI    * abh2    ! Collisional dissoc.: H2 + H 
     $           - cdh2_h2* abh2    * abh2    ! Collisional dissoc.: H2 + H2
     $           - cdh2_he* abheI   * abh2    ! Collisional dissoc.: H2 + He
     $           - ch9    * abhep   * abh2    ! DiCT: He+ + H2 -> He + H + H+
     $           - cst(75) * abhep  * abh2    ! CT: He+ + H2 -> He + H2+
     $           - ch126   * abh2   * abe     ! H2 + e -> H- + H
     $           + ch127   * abhm   * abh2p) * yn ! H2+ + H- -> H2 + H
     $           -  ykdis  * abh2             ! Photodissoc:         H2 + PH
     $           -  cr2    * abh2             ! Cosmic ray ion.:     H2 + CR
     $           - ph41    * abh2             ! Photoion:           H2  + PH

      if (no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        ydot(ih2) = ydot(ih2)
     $            + (ch101  * abhd    * abhp       ! HD destruction:      HD + H+
     $            -  ch100  * abh2    * abdp       ! HD formation:        D+ + H2
     $            -  ch102  * abh2    * abDI       ! HD formation:         D + H2
     $            +  ch103  * abhd    * abHI) * yn ! HD destruction:      HD + H
      endif
c
c  Ionized deuterium (D+)
c
      if (abundD .le. 0d0 .or. no_dchem .eq. 1) then
        ydot(idp) = 0d0
      else
        ydot(idp) = (ch98   * abhp  * abDI          ! Charge transfer: D + H+
     $             - ch99   * abdp  * abHI          ! Charge transfer: H + D+
     $             - ch100  * abdp  * abh2          ! HD form: D+ + H2
     $             + ch101  * abhd  * abhp          ! HD dest: HD + H+
     $             - ch13   * abdp  * abe           ! Case B recomb.
     $             + ch11   * abDI  * abe           ! Collisional ionization
     $             + ch32   * abhep * abDI          ! Charge transfer: He+ + D
     $             - ch33   * abheI * abdp) * yn    ! Charge transfer: He  + D+
     $             + cr1a   * abDI                  ! Cosmic ray ion: D + CR
     $             + ph38   * abDI                  ! Photoionization
      endif
c
c  Hydrogen deuteride (HD)
c
c  NB We follow Flower & Harris (2007) and assume that the three-body rates for
c  HD formation are the same as those for H2 formation
c
c  NB For simplicity, we assume that each photoionization leads to the destruction
c  of the molecule; this is incorrect when the fractional ionization is small, but
c  we do not expect HD photoionization to be important in those conditions
c
      if (abundD .le. 0d0 .or. no_dchem .eq. 1) then
        ydot(ihd) = 0d0
      else
        ydot(ihd) = (ch108  * abHI
     $            +  ch109  * abh2
     $            +  ch110  * abheI) * abHI * abDI * ynsq
     $            + (ch100  * abh2  * abdp          ! HD form: D+ + H2
     $             - ch101  * abhd  * abhp          ! HD dest: HD + H+
     $             + ch102  * abh2  * abDI          ! HD form: D + H2
     $             - ch103  * abhd  * abHI) * yn    ! HD dest: HD + H
     $             - ph35   * abhd                  ! Photodissociation
     $             - ph41   * abhd                  ! Photoionization
      endif
c
c Singly-ionized helium (He+)
c
      if (abhe .eq. 0d0) then
        ydot(ihep) = 0d0
      else
        ydot(ihep) = (ch14    * abheI * abe
     $              - ch15    * abhep * abe
     $              - ch17    * abhep * abe
     $              + ch20    * abhepp * abe
     $              - ch32    * abhep * abHI
     $              + ch33    * abheI * abhp
     $              -  ch8    * abhep * abhm
     $              -  ch9    * abhep * abh2
     $              - cst(75) * abhep * abh2) * yn
     $              + cr6     * abheI
     $              + ph39    * abheI
     $              - ph40    * abhep
      endif
c
      if (no_dchem .eq. 0 .and. abundD .gt. 0d0) then
        ydot(ihep) = ydot(ihep)
     $             + (ch33    * abheI * abdp
     $             -  ch32    * abhep * abDI) * yn
      endif
c
c Doubly-ionized helium (He++)
c
      if (abhe .eq. 0d0) then
        ydot(ihepp)  = 0d0
      else
        ydot(ihepp) = (ch15 * abhep  * abe
     $              -  ch20 * abhepp * abe) * yn
     $              +  ph40 * abhep
      endif
c
c Chemical heating & cooling
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (cdh2_h * abHI * abh2 + 
     $                               cdh2_h2* abh2 * abh2 +
     $                               cdh2_he* abheI * abh2 +
     $                               ch25   * abe  * abh2 ) * ynsq
c
c H2 destruction by charge transfer with HII -- endothermic by 1.83eV 
c
      rates_chem(2) = 1.83d0 * eV * ch24 * abh2 * abhp * ynsq
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(3) = - 4d-1 * eV * abH2 * yn * ykdis
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

      rates_chem(4) = - 2d0 * eV * abH2 * yn * ykdis * 5.67d0 * fpheat
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to 
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c Note that energy input depends on formation mechanism:
c
c H-   -- 3.73eV
c H2+  -- 1.83eV
c 3-body -- 4.48eV
c
      rates_chem(5) = - 3.73d0 * eV * ch22 * abHI * abhm * ynsq
     $                * (1d0 - h2var0)
c
      rates_chem(6) = - 1.83d0 * eV * cst(1) * abHI * abh2p * 
     $                ynsq * (1d0 - h2var0)
c
      rates_chem(7) = - 4.48d0 * eV * (ch108 * abHI**3 + 
     $                  ch109 * abHI**2 * abh2 + 
     $                  ch110 * abHI**2 * abheI) * ynsq * yn
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(8) = 13.6d0 * eV * (ch11 * abe + ch144 * abHI
     $              + ch158 * abheI) * abHI * ynsq
c
c HeI collisional ionization cooling -- 24.6 eV per ionization
c
      rates_chem(9) = 24.6d0 * eV * ch14 * abe * abheI * ynsq
c
c HeII collisional ionization cooling -- 54.4 eV per ionization
c
      rates_chem(10) = 54.4d0 * eV * ch15 * abe * abhep * ynsq
c
c HII recombination cooling -- we assume that kT ergs is radiated 
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
      rates_chem(11) = kboltz * temp * abhp * ch13 * abe * 
     $                 ynsq
c
c HeII recombination cooling. Dielectronic recombination results
c in the loss of ~ 40.8eV per recombination, while for standard 
c radiative recombination, we make the same assumption as for HII
c
c Dielectronic:
      rates_chem(12) = 40.8d0 * eV * ch10 * abhep * abe * ynsq
c
c Radiative:
      rates_chem(13) = kboltz * temp * (ch17 - ch10) * abhep * abe * 
     $                 ynsq
c
c HeIII recombination cooling -- see note for HII
c
      rates_chem(14) = kboltz * temp * abhepp * ch20 * abe * 
     $                 ynsq
c
c H- photodissociation heating
c
      rates_chem(15) = - ph8 * abhm * yn
c
c H2+ photodissociation heating
c
      rates_chem(16) = - ph9 * abh2p * yn
c
c HI photoionization heating
c 
      rates_chem(17) = - ph42 * abHI * yn
c
c HeI photoionization heating
c 
      rates_chem(18) = - ph43 * abheI * yn
c
c HeII photoionization heating
c 
      rates_chem(19) = - ph44 * abhepp * yn
c
c H2 photoionization heating
c
      rates_chem(20) = - ph45 * abh2 * yn
c
c H- formation cooling (NB kinetic contribution is approximate)
c
      rates_chem(21) = (0.75d0 * eV + kboltz * temp) 
     $               * ch21 * abHI * abe * ynsq 
c
      ylam_chem = 0d0
      do I = 1, nrates_chem
        ylam_chem = ylam_chem + rates_chem(I)
      enddo
c
c Initialize all abundances to zero, then just set the non-zero ones
c
      do i = 1, nabn
        abundances(i) = 0d0
      enddo
c
      if (no_dchem .eq. 1) then
        abhd = min(abundD, exp(4.88d2 / temp) * abh2)
      endif
c
      abundances(1)  = abh2
      abundances(10) = abe
      abundances(11) = abhp
      abundances(12) = abHI
      abundances(13) = abhd
      abundances(14) = abheI
      abundances(15) = abhep
      abundances(16) = abhepp
c
#ifdef ADIABATIC_DENSITY_THRESHOLD
      if (yn .ge. yn_adiabatic) then
        ydot(itmp) = 0d0
#ifdef DEBUG_RATE_EQ
        call validate_output(y,ydot,rpar)
#endif
        return
      endif
#endif
c
      call cool_func(temp, yn, L_jeans, divv, G0_curr, abundances, ylam,
     $               rates, 0, energy)
c
      if(ylam .ne. ylam) then
       print*,'RATE_EQ2: ylam NaN:', temp, yn, L_jeans, divv,
     $        'abund:', abundances, ylam,
     $        'rates:', rates,'energy:', energy
      endif
c
c Now correct for effects of CMB:
c
      if (redshift .gt. 0d0) then
        call cool_func(CMB_temp, yn, L_jeans, divv, G0_curr, abundances,
     $                 ylam_cmb, rates_cmb, 1, energy)
        if (ylam_cmb .ne. ylam_cmb) then
          print*,'RATE_EQ3: ylam_cmb NaN:', CMB_temp, yn, L_jeans, divv,
     $           'abund:', abundances,'ylam_cmb:', ylam_cmb,
     $           'rates:', rates_cmb,'energy:', energy
        endif
      else
        ylam_cmb = 0d0
      endif

      ydot(itmp) = -(ylam - ylam_cmb + ylam_chem) 

#ifdef THERMAL_INFO_DUMP
      dtcool_nopdv = rates(29) !hitchhiked H2_opacity_correction
      cool_h2_line = rates(2) - rates_cmb(2)
      cool_h2_cie  = rates(20) - rates_cmb(20)
      cool_h2_diss = rates_chem(1)
      heat_3b      = rates_chem(5) + rates_chem(6) + rates_chem(7)
#endif
c
#ifdef DEBUG_RATE_EQ
      if (ydot(itmp) .ne. ydot(itmp)) then
        print*, "Cooling rate is NaN: y=",y,' rpar=',rpar,' ylam=',ylam,
     $         ' ylam_cmb=',ylam_cmb,' ylam_chem=',ylam_chem,' rates='
     $        ,rates,' rates_chem=',rates_chem
      endif
      call validate_output(y,ydot,rpar)
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
c=======================================================================
c
c    \\\\\\\\\\        B E G I N   F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
c
#ifdef SGCHEM
      REAL function comp_equilib(num, denom, max_abundance, id, yn, 
     $                           time_to_equilibrium)
      implicit NONE
#include "cool.h"
#include "non_eq.h"
       REAL num, denom, max_abundance, yn
       REAL time_to_equilibrium(num_eqb_species)
       integer id
c
       if (num .eq. 0d0) then
         comp_equilib = 0d0
       elseif (denom .eq. 0d0) then
c If this happens, then our assumption that this species is in
c equilibrium has fatally broken down...
         print*, "Creation without destruction in comp_equilib"
         print*, "Species = ", id
c         ABORT(9751)
         comp_equilib = -1d0
      else
#ifdef COLLECT_EQB_DATA
         time_to_equilibrium(id) = 1d0 / (denom * yn)
#endif
        comp_equilib = num / denom
        if (comp_equilib .gt. max_abundance) then
          print*, "COMP_EQUI Error: impossible abundance!"
          print*, "Species ", id
          print*, num, denom, num/denom
c          ABORT(9752)
          comp_equilib = -1d0
        endif
      endif
c
      return
      end
#endif /* SGCHEM */
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////          C O M P _ E Q U I L I B          \\\\\\\\\\
c
c=======================================================================
#endif /* CHEMISTRYNETWORK == 1 */

