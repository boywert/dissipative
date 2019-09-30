#include "../../build/arepoconfig.h"

#if defined TREE_RAD || defined TREE_RAD_H2 || defined TREE_RAD_CO
#define NPIX  12*NSIDE*NSIDE

 subroutine project_column(column, columnH2, columnCO, dx, dy, dz, angular_radius, projection, projectionH2, projectionCO, id)
 use healpix

 implicit none

 real(KIND=DP), intent(inout) :: column, columnH2, columnCO, dx, dy, dz, angular_radius
 real(KIND=DP), dimension(0:NPIX-1), intent(inout) :: projection, projectionH2, projectionCO
 real(KIND=DP) :: rp, tworp, one_over_rp2
 real(KIND=DP), dimension(0:2, 0:NPIX-1) :: pvector, pnode
 real(KIND=DP) :: atheta, aphi, Rphi, Rtheta, rn
 integer(kind=I4B) :: i
 real*8  :: px, py, pz
 real*8  :: PYx, PYy, PYz, PZx, PZy, PZz, PYn, PZn
 real*8 :: newx, newy, newz
 real*8 :: ax, ay, az, anorm
 real(KIND=DP) :: nodesum, tworn, areanorm, weighting
 real(KIND=DP) :: acos_pcc, asin_pcc
 integer :: id
 
 common / pixel_share / pvector
 common / healpix_res / rp, tworp, one_over_rp2

 rn = angular_radius
 tworn = 2.*rn
 nodesum = rp + rn
 if (rn > rp ) then
    areanorm = tworp
 else
    areanorm = tworn
 end if
	
! if (id.eq.7746) then	
!    print*,'NPIX fortran',NPIX
!    print*,'Size', size(projection),size(projectionH2),size(projectionCO)
!    print*,'Lower and upper bound ', lbound(projection),ubound(projection)	
! endif
!
! if (size(projection).ne.12*NSIDE*NSIDE) then
!	print*,"ERROR in size of projection ",size(projection),NPIX,NSIDE
! endif
!
! if (size(projectionH2).ne.12*NSIDE*NSIDE) then
!	print*,"ERROR in size of projectionH2 ",size(projectionH2),NPIX,NSIDE
! endif
!
! if (size(projectionCO).ne.12*NSIDE*NSIDE) then
!	print*,"ERROR in size of projectionCO ",size(projectionCO)
! endif

 !
 ! Find a vector that is close to dx, dy, dz, and use as the control
 ! vector for defining the new co-ordinate axis.
 !

 ax = dx + 2.*rp
 ay = dy
 az = dz + 2.*rp
 anorm = 1./sqrt(ax*ax + ay*ay + az*az)
 ax = ax*anorm
 ay = ay*anorm
 az = az*anorm

 !
 ! dx, dy, dz is the new x-axis, so create new y and z axis, using
 ! the control vector (ax, ay, az) as a starting point.
 !

 PYx = dy*az - dz*ay
 PYy = -(dx*az - dz*ax)
 PYz = dx*ay - dy*ax
 PYn = 1./sqrt(PYx*PYx + PYy*PYy + PYz*PYz)
 PYx = PYx*PYn
 PYy = PYy*PYn
 PYz = PYz*PYn
             
 PZx = dy*PYz - dz*PYy
 PZy = -(dx*PYz - dz*PYx)
 PZz = dx*PYy - dy*PYx

 !
 ! To save time, we don't use HEALPIX'S  query_disc subroutine anymore.
 ! We simply loop around all pixels and calculate the angle to them
 ! using the dot product. The angles here are computed using a lookup table. 
 ! 

 do i = 0, NPIX-1
    !
    ! Transform all the pixels into the node's coordinate system
    ! can skip if new x is < 0 since then phi > 90deg.
    !
    px = pvector(0, i)
    py = pvector(1, i)
    pz = pvector(2, i)
    newx = px*dx + py*dy + pz*dz

    if ( newx.lt.0 ) cycle
    newy = px*PYx + py*PYy + pz*PYz
    newz = px*PZx + py*PZy + pz*PZz

    !
    ! Calculate the orthogonal angular distances between pixel
    ! and node centres.
    !
    Rtheta = abs(asin_pcc(newz))
    if ( Rtheta .gt. nodesum ) cycle
    Rphi = acos_pcc(newx/sqrt(newx*newx + newy*newy))
    if ( Rphi .gt. nodesum ) cycle
    !
    ! There is some overlap between the node and the pixel
    !
    atheta = min(nodesum - Rtheta, areanorm)
    aphi = min(nodesum - Rphi, areanorm)
    weighting = 0.25*atheta*aphi*one_over_rp2
#ifdef TREE_RAD
    projection(i) = projection(i) + column*weighting
!crjs    
!    if (id.eq.7746) then	
!       print*,'F90 Projection ', i, projection(i), column, weighting
!    endif
#endif
#ifdef TREE_RAD_H2
    projectionH2(i) = projectionH2(i) + columnH2*weighting
#endif
#ifdef TREE_RAD_CO
    projectionCO(i) = projectionCO(i) + columnCO*weighting
#endif
 enddo

 !
 ! Done!
 !

 return
 end subroutine project_column 
!____________________________________________________________________________________
 
 subroutine get_angular_coords(ipix, theta, phi)
 use healpix
 integer, intent(inout) :: ipix
 real*8, intent(inout) :: theta, phi

 call pix2ang_ring(NSIDE, ipix, theta, phi)

 return
 end subroutine get_angular_coords

 subroutine get_pixels_for_xyz_axes(pixels)
 use healpix
 integer, dimension(6), intent(inout) :: pixels

 real*8, dimension(3) :: vector
 integer :: I, ipring

 do I = 1, 6
   vector(1) = 0.0
   vector(2) = 0.0
   vector(3) = 0.0
   select case (I)
   case (1)
     vector(1) = 1.0
   case (2)
     vector(1) = -1.0
   case (3)
     vector(2) = 1.0
   case (4)
     vector(2) = -1.0
   case (5)
     vector(3) = 1.0
   case (6)
     vector(3) = -1.0
   end select
   call vec2pix_ring(NSIDE, vector, ipring)
   pixels(I) = ipring
 enddo

 return
 end subroutine get_pixels_for_xyz_axes
!__________________________________________________________________________________
subroutine calculate_pixel_centres

 use healpix

 implicit none
 
 real(KIND=DP), dimension(1:3) :: pixel_centre
 real(KIND=DP), dimension(0:2, 0:NPIX-1) :: pvector
 real(KIND=DP) :: rp, tworp, one_over_rp2
 integer i

 common / pixel_share / pvector
 common / healpix_res / rp, tworp, one_over_rp2
 
 do i = 1, NPIX
    call pix2vec_ring(NSIDE, i-1, pixel_centre)
    pvector(0, i-1) = pixel_centre(1)
    pvector(1, i-1) = pixel_centre(2)
    pvector(2, i-1) = pixel_centre(3)
 end do

 !
 ! Set the angular distance between healpix pixels.
 !

 rp = sqrt(12.5663706143592/real(NPIX))/2.
 tworp = 2.*rp
 one_over_rp2 = 1./rp/rp

end subroutine calculate_pixel_centres
!_____________________________________________________________________________
subroutine createtriglookup

!
! create loop-up tables on each thread for various trig functions
!

 use healpix

 integer :: i
 real(KIND=DP) :: dx
 real(KIND=DP) :: acos_table(8001), asin_table(8001)
 common / trig_tables / acos_table, asin_table

 dx = 2.0/8000.

 do i = 1, 8000
    acos_table(i) = acos( -1 + real(i-1)*dx )
    asin_table(i) = asin( -1 + real(i-1)*dx )
 end do
 
 acos_table(8001) = 0.0
 asin_table(8001) = 1.5707963267949


end subroutine createtriglookup
!_____________________________________________________________________________
function acos_pcc(x)
 
!
! The function that does the look-up for the inverse cos.
! Should be fed with a number between -1 and 1
!

 use healpix

 integer index
 real(KIND=DP) :: acos_pcc
 real(KIND=DP) :: x
 real(KIND=DP) :: acos_table(8001), asin_table(8001)
 common / trig_tables / acos_table, asin_table
 
 index = int((x+1)*4000) + 1
 !index = (int(x)+1)*4000 + 1
 index = min(index, 8001)
 index = max(index, 1) 
 acos_pcc = acos_table(index)

end function acos_pcc

!____________________________________________________________________________
       
 function asin_pcc(x)
        
 !
 ! The function that does the look-up for the inverse cos.
 ! Should be fed with a number between -1 and 1
 !

 implicit none
 
 integer index
 real*8 :: asin_pcc
 real*8  :: x
 real*8 :: acos_table(8001), asin_table(8001)
 common / trig_tables / acos_table, asin_table
        
 index = int((x+1)*4000) + 1
 !index = (int(x)+1)*4000 + 1
 index = min(index, 8001)
 index = max(index, 1) 
 asin_pcc = asin_table(index)

 end function asin_pcc


#endif /* TREE_RAD || TREE_RAD_H2 || TREE_RAD_CO */
