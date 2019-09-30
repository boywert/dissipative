#include "../../build/arepoconfig.h"
#ifdef SGCHEM
      subroutine calc_shielding(yn, dl, temp, abh2, abco,
     $                          column_density_projection,
     $                          column_density_projection_H2,
     $                          column_density_projection_co,
     $                          fshield_H2, fshield_CO, AV_mean, 
     $                          chi_mean)
      implicit NONE
#include "cool.h"
      REAL yn, dl, abh2, abco, temp
      DOUBLE PRECISION column_density_projection(NPIX)
      DOUBLE PRECISION column_density_projection_H2(NPIX)
      DOUBLE PRECISION column_density_projection_co(NPIX)
      REAL fshield_H2, fshield_CO, AV_mean, chi_mean, AV, chi
      REAL local_NH, AV_local, local_NH2, local_NCO, f_dust
      external chi
      integer I, nray, nray_start, nray_end
c
c Optically thin approx:
c
      if (iphoto .eq. 0) then
        fshield_H2 = 1.0
        fshield_CO = 1.0
        AV_mean    = 0.0
        chi_mean   = 1.0
        return
      endif
c
      if (iphoto .ne. 6) then
        print*, "Unknown photochemistry option: ", iphoto
        ABORT(3)
      endif

      if (isnan(AV)) print*,"calc_shield.F AV NaN",AV,yn
c
c Total column density of HI nuclei within characteristic scale-length dl
c
      local_NH = dl * yn
      AV_local = AV_conversion_factor * dust_to_gas_ratio 
     $         * local_NH
c
c H2 self-shielding
c
      local_NH2 = local_NH * abh2
      if (local_NH2.lt.0.0) then
         print*,"local_NH2 very wrong ",local_NH2,local_NH,dl,yn,
     &        abh2
         local_NH2=0.0
      endif
      call calc_H2_shielding_factor(temp, local_NH2, 
     $                              column_density_projection_h2, 
     $                              fshield_H2)

      if (isnan(fshield_H2)) print*,"calc_shield.F fshield_H2 NaN",
     $         fshield_H2,yn

c
c CO self-shielding, shielding of CO by H2
c
#if CHEMISTRYNETWORK != 1
      local_NCO = local_NH * abco
      call calc_CO_shielding_factor(column_density_projection_h2, 
     $                              column_density_projection_co, 
     $                              local_NCO, local_NH2, fshield_CO)
#else
      fshield_CO = 1d0
#endif
c
c Dust shielding -- we compute a mean AV, assuming a weighting factor exp(-2.5AV); this
c is appropriate for photoelectric heating and CO photodissociation, shouldn't be too
c much in error for other things
c
#if CHEMISTRYNETWORK != 1
      nray     = NPIX
      f_dust   = 0d0
      chi_mean = 0d0
      do I = 1, nray
        AV = AV_conversion_factor * column_density_projection(I) 
     $     * dust_to_gas_ratio + AV_local + AV_ext
        if (isnan(AV)) print*,"calc_shield.F AV loop NaN",AV,I,'/',nray,
     $       column_density_projection(I)
        f_dust   = f_dust + (1d0 / nray) * dexp(-2.5d0 * AV)
        chi_mean = chi_mean + (1d0 / nray) * chi(AV)
      enddo
      AV_mean  = -log(f_dust) / 2.5d0
#else
      chi_mean = 1d0
      AV_mean  = 0d0
#endif
#endif /* SGCHEM */
      return
      end
c=======================================================================
c
c    \\\\\\\\\\          E N D     F U N C T I O N        //////////
c    //////////             C A L C _ P H O T O           \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   F U N C T I O N          //////////
c    //////////                C H I                      \\\\\\\\\\
c
c=======================================================================
c
      REAL function chi(AV)
#ifdef SGCHEM
      implicit NONE
#include "cool.h"
#include "isrf.h"
      REAL AV, dchi, dAV
      integer index

      if (AV .ne. AV) then
        print*, "Error: AV is NaN"
        stop
        
      endif

      if (AV .eq. 0) then
        chi = 1d0
        return
      elseif (AV .gt. 2d2) then
c Assume chi remains constant for AV > 200
        if (isrf_option .eq. 1) then
          chi = chi_tab_mathis(200)
        else
          chi = chi_tab_black(200)
        endif
        return
      else
        index = int(AV)
        if ( index.eq. 200 ) then
           print *, "Index = 200... Whoops."
           stop
        end if 
        if (isrf_option .eq. 1) then
          if (index .eq. 0) then
            dchi = chi_tab_mathis(1) - 1d0
            dAV  = AV - index
            chi  = 1d0 + dAV * dchi
          else
            dchi = chi_tab_mathis(index+1) - chi_tab_mathis(index)
            dAV  = AV - index
            chi  = chi_tab_mathis(index) + dAV * dchi
          endif
        else
          if (index .eq. 0) then
            dchi = chi_tab_black(1) - 1d0
            dAV  = AV - index
            chi  = 1d0 + dAV * dchi
          else
            dchi = chi_tab_black(index+1) - chi_tab_black(index)
            dAV  = AV - index
            chi  = chi_tab_black(index) + dAV * dchi
          endif
        endif
        return
      endif
#else
      chi = 1.0
#endif /* SGCHEM */
      return
      end
c
c=======================================================================
c
c
c    \\\\\\\\\\          E N D   F U N C T I O N          //////////
c    //////////                  C H I                    \\\\\\\\\\
c
c=======================================================================
c
      subroutine calc_H2_shielding_factor(temp, local_NH2,
     $                                    column_density_projection_h2, 
     $                                    fshield_H2)

#ifdef SGCHEM
      implicit none
#include "cool.h"
      REAL fshield_H2, temp, local_NH2
      REAL column_density_projection_H2(NPIX)
      REAL b, b5, h2dis0, h2dis1, alpha
      integer I, nray
c
c alpha depends on choice of H2 self-shielding function
c In the original DB96 function, for rotationally cold H2, alpha = 2.0
c In the revised WGH11 function, for rotationally hot H2,  alpha = 1.1
c
      alpha = 2.0
c
      nray = NPIX
c
c To compute H2 self-shielding, we assume that the local b param. is 
c representative of the whole flow. This isn't true in general, but 
c is still a minor approximation compared to the others we are making
c
      b  = dsqrt(kboltz * temp / mh)
      b5 = b / 1d5
      fshield_H2 = 0d0
      do I = 1, NPIX
        h2dis0  = (local_NH2 + column_density_projection_h2(I)) / 5d14
        h2dis1  = dsqrt(1d0 + h2dis0)
        fshield_H2 = fshield_H2 + (1d0 / nray) * (9.65d-1 /
     $             ((1d0 + h2dis0 / b5)**alpha)
     $             + 3.5d-2 * dexp(-8.5d-4 * h2dis1) / h2dis1)
        if (isnan(fshield_H2)) print*, "Calc H2 shielding", I, local_NH2
     &       ,column_density_projection_h2(I),b5,h2dis0
      enddo
#endif /* SGCHEM */
      return
      end
c
c
c
      subroutine calc_CO_shielding_factor(column_density_projection_h2, 
     $                                    column_density_projection_co, 
     $                                    local_NCO, local_NH2, 
     $                                    fshield_CO)
#ifdef SGCHEM
      implicit none
#include "cool.h"
#include "shield_data.h"
      REAL column_density_projection_h2(NPIX)
      REAL column_density_projection_co(NPIX)
      REAL fshield_CO, theta_CO, theta_H2, NCO, NH2
      REAL dth, dNact, local_NCO, local_NH2
      integer nray, I, ipos
c
c CO self-shielding, shielding by H2.
c
c Our treatment here is based on Lee et al (1996)
c
      fshield_CO = 1d0

      nray = NPIX
      theta_CO   = 0d0
      theta_H2   = 0d0
      do I = 1, nray
        NCO = column_density_projection_CO(I) + local_NCO
        if (NCO .le. NCO_shield_min) then
          theta_CO = theta_CO + 1d0 / nray
        elseif (NCO .ge. NCO_shield_max) then
          theta_CO = theta_CO + (1d0 / nray)
     $             * CO_self_shielding(imax_COss)
        else
          ipos  = int(dlog10(NCO) / dNshield) - ioff_co
          if (ipos .eq. 0) then
            theta_CO = theta_CO + 1d0 / nray
          else
            dth = CO_self_shielding(ipos+1) - 
     $            CO_self_shielding(ipos)
            dNact = dlog10(NCO) - (ipos + ioff_co) * dNshield
            theta_CO = theta_CO + (1d0/nray)
     $               * (CO_self_shielding(ipos) 
     $               + dth * (dNact / dNshield))
          endif
        endif
c
        NH2 = column_density_projection_h2(I) + local_NH2
        if (NH2 .le. NH2_shield_min) then
          theta_H2 = theta_H2 + 1d0 / nray
        elseif (NH2 .ge. NH2_shield_max) then
          theta_H2 = theta_H2 + (1d0 / nray) 
     $             * CO_H2_shielding(imax_COH2)
        else
          ipos  = int(dlog10(NH2) / dNshield) - ioff_h2
          if (ipos .eq. 0) then
            theta_H2 = theta_H2 + 1d0 / nray
          else
            dth   = CO_H2_shielding(ipos+1) - 
     $              CO_H2_shielding(ipos)
            dNact = dlog10(NH2) - (ipos + ioff_h2) * dNshield
            theta_H2 = theta_H2 + (1d0 / nray)
     $               * (CO_H2_shielding(ipos) 
     $               + dth * (dNact / dNshield))
          endif
        endif
      enddo

      fshield_CO = theta_H2 * theta_CO

#endif /* SGCHEM */
      return
      end
