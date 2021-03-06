module taufind2
   implicit none
   contains
subroutine tau2(xcell, ycell, zcell, delta, taurun)

      use iarray,      only : rhokap, xface, yface, zface, refrac
      use photon_vars, only : xp, yp, zp, nxp, nyp, nzp
      use constants,   only : xmax, ymax, zmax, nxg, nyg, nzg
      use opt_prop,    only : wavelength, material

      implicit none

      integer, intent(IN)  :: xcell, ycell, zcell
      real,    intent(IN)  :: delta
      real,    intent(OUT) :: taurun
 
      integer :: celli,cellj,cellk
      real :: xcur,ycur,zcur,dx,dy,dz,d,dcell,n1
      real :: smax,dsx,dsy,dsz,taucell,tau1


!***** set the cumulative distance and optical depth (d and taurun) 
!***** along the photon path to zero.  set the current photon coordinates.
!***** note that the origin of the (xcur,ycur,zcur) system is at the 
!***** bottom corner of the grid.
      taurun = 0.
      taucell = 0.
      d = 0.
      dcell = 0.
      xcur = xp + xmax
      ycur = yp + ymax
      zcur = zp + zmax

      celli = xcell
      cellj = ycell
      cellk = zcell

!***** calculate smax -- maximum distance photon can travel *******
      if(nxp .gt. 0.) then
         dsx = (2.*xmax-xcur)/nxp
      elseif(nxp.lt.0.) then
         dsx = -xcur/nxp
      elseif(nxp.eq.0.) then
         dsx = 1.e2*xmax
      endif

      if(nyp .gt. 0.) then
         dsy = (2.*ymax-ycur)/nyp
      elseif(nyp.lt.0.) then
         dsy = -ycur/nyp
      elseif(nyp.eq.0.) then
         dsy = 1.e2*ymax
      endif

      if(nzp .gt. 0.) then
         dsz = (2.*zmax-zcur)/nzp
      elseif(nzp.lt.0.) then
         dsz = -zcur/nzp
      elseif(nzp.eq.0.) then
         dsz = 1.e2*zmax
      endif

      smax = min(dsx,dsy,dsz)

      if(smax.lt.delta) then
         tau1 = smax*rhokap(celli,cellj,cellk,wavelength + material)
         return
      endif
!***** integrate through grid
      do while(d.lt.(.999*smax))

!***** optical depth to next cell wall is 
!***** taucell =  (distance to cell)*(opacity of current cell)

      n1 = refrac(celli,cellj,cellk)

!***** find distance to next x, y, and z cell walls.  
!***** note that dx is not the x-distance, but the actual distance along 
!*****the direction of travel to the next x-face, and likewise for dy and dz.
         if(nxp .gt. 0.) then
            dx = (xface(celli + 1)-xcur)/nxp
            if(dx.lt.delta) then
               xcur = xface(celli + 1)
               celli = celli + 1
               dx = (xface(celli + 1)-xcur)/nxp
            endif
         elseif(nxp.lt.0.) then
            dx = (xface(celli)-xcur)/nxp
            if(dx.lt.delta) then
               xcur = xface(celli)
               dx = (xface(celli-1)-xcur)/nxp
               celli = celli-1
            endif
         elseif(nxp.eq.0.) then
            dx = 1.e2*xmax
         endif

         if(nyp .gt. 0.) then
            dy = (yface(cellj + 1)-ycur)/nyp
            if(dy.lt.delta) then
               ycur = yface(cellj + 1)
               cellj = cellj + 1
               dy = (yface(cellj + 1)-ycur)/nyp
            endif
         elseif(nyp.lt.0.) then
            dy = (yface(cellj)-ycur)/nyp
            if(dy.lt.delta) then
               ycur = yface(cellj)
               dy = (yface(cellj-1)-ycur)/nyp
               cellj = cellj-1
            endif
         elseif(nyp.eq.0.) then
            dy = 1.e2*ymax
         endif

         if(nzp .gt. 0.) then
            dz = (zface(cellk + 1)-zcur)/nzp
            if(dz.lt.delta) then
               zcur = zface(cellk + 1)
               cellk = cellk + 1
               dz = (zface(cellk + 1)-zcur)/nzp
            endif
         elseif(nzp.lt.0.) then
            dz = (zface(cellk)-zcur)/nzp
            if(dz.lt.delta) then
               zcur = zface(cellk)
               dz = (zface(cellk-1)-zcur)/nzp
               cellk = cellk-1
            endif
         elseif(nzp.eq.0.) then
            dz = 1.e2*zmax
         endif

         if(n1 == 1.38)then
               material = 1
         else
               material = 3
         end if

!***** distances are only zero if photon is on cell wall.  if it is 
!***** on cell wall then set to arbitrary large distance, since we will
!***** in fact hit another wall
         if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx = 1.e2*xmax
         if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy = 1.e2*ymax
         if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz = 1.e2*zmax

!***** find distance to next cell wall -- minimum of dx, dy, and dz
         dcell = min(dx,dy,dz)
         if(dcell.le.0.) then
            print *,'taufind1: dcell < 0'
!            stop
         endif
         if(dx.lt.0.) dcell = min(dy,dz)
         if(dy.lt.0.) dcell = min(dx,dz)
         if(dz.lt.0.) dcell = min(dx,dy)

         d = d + dcell
         taucell = dcell*rhokap(celli,cellj,cellk,wavelength + material)
            
!***** if taurun + taucell<tau then photon moves distance dcell 
!***** (i.e. ends up on next cell wall) and update photon position
!***** and cell.
         taurun = taurun + taucell
         xcur = xcur + dcell*nxp
         ycur = ycur + dcell*nyp
         zcur = zcur + dcell*nzp

!*************** Linear Grid ************************
         celli = int(nxg*xcur/(2.*xmax)) + 1
         cellj = int(nyg*ycur/(2.*ymax)) + 1
         cellk = int(nzg*zcur/(2.*zmax)) + 1
!****************************************************
      end do

      tau1 = taurun
      end subroutine tau2
end module taufind2