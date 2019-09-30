#include "../../build/arepoconfig.h"
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////                   L W B G                 \\\\\\\\\\
c
c=======================================================================
c
      subroutine lwbg
c
c    written by: A. Schauer, ITA, Dec 2016
c
c    PURPOSE: Compute a time-dependen LW BG
c             Taken from photoinit_lowZ.F : 
c             This version assumes an optical/UV spectrum corresponding
c             to a 10^5 K black body, cut off at 13.6eV. Scaling is 
c             controlled by G0: the rates are normalized so that G0 = 1
c             corresponds to a flux density of 10^-21 erg/s/cm^2/Hz/sr at
c             the Lyman limit.
c             lwtype = 0: uniform LWBG throughout the simulation with 
c                         value G0
c             lwtype = 1: uniform LWBG, turned on at lwstart with the 
c                         value G0
c
c    REFERENCES:
c
c
      implicit NONE
#include "cool.h"
c
c     this is the first option: lwtype = 0 
c     corresponds to a uniform LWBG from start 
c     to end of the simulation
c
      if (lwtype .eq. 0) then 
        G0_curr = G0
        return
      endif
c
c     this is the second option: lwtype = 1
c     corresponds to a LWBG that is turned on at 
c     redshift lwstart
c
      if (lwtype .eq. 1) then
        if (redshift .le. lwstartz) then
c        if (redshiftanna .le. lwstart) then
          G0_curr = G0
        else
          G0_curr = 0
        endif
      endif
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////                   L W B G                 \\\\\\\\\\
c
c=======================================================================
