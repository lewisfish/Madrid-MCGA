MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xcell,ycell,zcell,iseed, delta)

   use constants, only : nxg,nyg,nzg,pi,twopi,xmax,ymax,zmax
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real, intent(IN)    :: delta
   real                   :: ran2


      zp = zmax*(2.*ran2(iseed)-1.)
      xp = xmax*(2.*ran2(iseed)-1.)
      yp = ymax!*(2.*ran2(iseed)-1.)
      
      phi  = 3.*pi/2.!twopi*ran2(iseed)  
      cosp = cos(phi)
      sinp = sin(phi)         
      sint = 1.
      cost = 0.
 
   
   nxp = sint * cosp  
   nyp = sint * sinp
   nzp = cost
   
   ! print*,nxp,nyp,nzp
   ! call exit(0)

   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************
   end subroutine sourceph
end MODULE sourceph_mod
