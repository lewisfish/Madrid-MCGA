MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed, j)

   use constants, only : nxg,nyg,nzg,pi,twopi
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   integer, intent(IN)    :: j
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2


      zp = zmax-(1.e-5*(2.*zmax/nzg))
      xp = xmax*(2.*ran2(iseed)-1.)
      yp = ymax*(2.*ran2(iseed)-1.)
      
      phi  = twopi*ran2(iseed)
      cosp = cos(phi)
      sinp = sin(phi)         
      sint = 0.
      cost = -1.
 
   
   nxp = sint * cosp  
   nyp = sint * sinp
   nzp = cost
   
   
   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************
   end subroutine sourceph
end MODULE sourceph_mod
