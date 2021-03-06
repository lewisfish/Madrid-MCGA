module inttau2

   implicit none
   
   private
   public :: tauint1, peeling

CONTAINS

    subroutine tauint1(xcell,ycell,zcell,tflag,iseed,delta, dflag)
    !optical depth integration subroutine
    !
    !
        use constants,   only : xmax, ymax, zmax,nxg, nyg
        use photon_vars, only : xp, yp, zp
        use iarray,      only : jmean, rhokap, absorb

        use opt_prop,    only : material, wavelength, mua
        use vector_class
        use ch_opt
   
        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell, iseed
        logical, intent(INOUT) :: tflag
        logical, intent(IN)    :: dflag

        real                   :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell, ran2
        integer                :: celli, cellj, cellk
        logical                :: dir(3)

        xcur = xp + xmax
        ycur = yp + ymax
        zcur = zp + zmax

        celli = xcell
        cellj = ycell
        cellk = zcell

        taurun = 0.
        d = 0.
        dir = (/.FALSE., .FALSE., .FALSE./)

        tau = -log(ran2(iseed))

        do
            dir = (/.FALSE., .FALSE., .FALSE./)
            dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
            taucell = dcell * rhokap(celli,cellj,cellk,wavelength+material)

            if(taurun + taucell < tau)then
                taurun = taurun + taucell
                d = d + dcell

    jmean(celli, cellj, cellk, wavelength+material) = jmean(celli, cellj, cellk, wavelength+material) + dcell
                if(material+wavelength == 3 .and. .not. dflag)then
                    call init_opt3()
            absorb(celli, cellj, cellk) = absorb(celli, cellj, cellk) + dcell*mua
                end if

                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta, tflag, iseed, &
                                tau, taurun)
            else

                dcell = (tau - taurun) / rhokap(celli,cellj,cellk,wavelength+material)
                d = d + dcell

    jmean(celli, cellj, cellk, wavelength+material) = jmean(celli, cellj, cellk, wavelength+material) + dcell
                if(material+wavelength == 3 .and. .not. dflag)then
                    call init_opt3()
            absorb(celli, cellj, cellk) = absorb(celli, cellj, cellk) + dcell*mua
                end if

                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, delta, tflag, iseed, &
                                tau, taurun)
                exit
            end if

            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                if(celli == -1 .or. cellj == -1)then
                    call repeat_bounds(celli, cellj, xcur, ycur, xmax, ymax, nxg, nyg, delta)
                    tflag = .false.
                    if(celli == -1 .or. cellj == -1 .or. tflag)then
                       print*,'error',celli,cellj,tflag
                    end if
                else
                    tflag = .true.
                    exit
                end if
            end if

        end do
   
        xp = xcur - xmax
        yp = ycur - ymax
        zp = zcur - zmax
        xcell = celli
        ycell = cellj
        zcell = cellk

    end subroutine tauint1
   

    real function wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
    !funtion that returns distant to nearest wall and which wall that is (x,y or z)
    !
    !
        use iarray,      only : xface, yface, zface
        use photon_vars, only : nxp, nyp, nzp

        implicit none

        real,    intent(INOUT) :: xcur, ycur, zcur
        logical, intent(INOUT) :: dir(:)
        integer, intent(INOUT) :: celli, cellj, cellk
        real                   :: dx, dy, dz


        if(nxp > 0.)then
            dx = (xface(celli+1) - xcur)/nxp
        elseif(nxp < 0.)then
            dx = (xface(celli) - xcur)/nxp
        elseif(nxp == 0.)then
            dx = 100000.
        end if

        if(nyp > 0.)then
            dy = (yface(cellj+1) - ycur)/nyp
        elseif(nyp < 0.)then
            dy = (yface(cellj) - ycur)/nyp
        elseif(nyp == 0.)then
            dy = 100000.
        end if

        if(nzp > 0.)then
            dz = (zface(cellk+1) - zcur)/nzp
        elseif(nzp < 0.)then
            dz = (zface(cellk) - zcur)/nzp
        elseif(nzp == 0.)then
            dz = 100000.
        end if

        wall_dist = min(dx, dy, dz)
        if(wall_dist < 0.)print'(A,7F9.5)','dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
        if(wall_dist == dx)dir=(/.TRUE., .FALSE., .FALSE./)
        if(wall_dist == dy)dir=(/.FALSE., .TRUE., .FALSE./)
        if(wall_dist == dz)dir=(/.FALSE., .FALSE., .TRUE./)
        if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, wall_flag, dir, delta, tflag, iseed, &
                          tau, taurun, flag)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        use photon_vars, only : nxp, nyp, nzp, phi, cost, sint, cosp, sinp
        use opt_prop,    only : material
        use iarray,      only : xface, yface, zface, refrac
        use constants,   only : nzg
        use vector_class

        implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur, tau, taurun
      real,    intent(IN)    :: dcell, delta
      integer, intent(INOUT) :: celli, cellj, cellk, iseed
      logical, intent(IN)    :: wall_flag, dir(:)
      logical, intent(INOUT) :: tflag
      logical, optional      :: flag
      
      type(vector) :: norm, incd
      real         :: n1, n2, ran2
      integer      :: iold, jold, kold
      logical      :: rflag

      iold = celli
      jold = cellj
      kold = cellk

      rflag = .false.

      n1 = refrac(celli,cellj,cellk)
      
      if(wall_flag)then
      
         if(dir(1))then
            if(nxp > 0.)then
               xcur = xface(celli+1) + delta
            elseif(nxp < 0.)then
               xcur = xface(celli) - delta
            else
               print*,'Error in x dir in update_pos', dir, nxp, nyp, nzp
               error stop 
            end if
            ycur = ycur + nyp*dcell 
            zcur = zcur + nzp*dcell
         elseif(dir(2))then
            xcur = xcur + nxp*dcell
            if(nyp > 0.)then
                ycur = yface(cellj+1) + delta
            elseif(nyp < 0.)then
                ycur = yface(cellj) - delta
            else
                print*,'Error in y dir in update_pos', dir, nxp, nyp, nzp
                error stop 
            end if
            zcur = zcur + nzp*dcell
         elseif(dir(3))then
            xcur = xcur + nxp*dcell
            ycur = ycur + nyp*dcell 
            if(nzp > 0.)then
               zcur = zface(cellk+1) + delta
            elseif(nzp < 0.)then
               zcur = zface(cellk) - delta
            else
               print*,'Error in z dir in update_pos', dir, nxp, nyp, nzp
               error stop 
            end if
         else
            print*,'Error in update_pos...',dir
            error stop 
         end if
      else
      
        xcur = xcur + nxp*dcell
        ycur = ycur + nyp*dcell 
        zcur = zcur + nzp*dcell
      
      end if


      if(wall_flag)then
         call update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
         if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            tflag = .true.
            return
        end if
        n2 = refrac(celli,cellj,cellk)
      end if

      if(n1 == 1.38)then
            material = 1
      else
            material = 3
      end if

      if(.not. present(flag))then
      if(nzp .gt. 0. .and. cellk == -1)then
        cellk = nzg+1
        n2 = refrac(celli,cellj,cellk)
      end if

      if(wall_flag)then
        if(n1 /= n2)then
            incd = vector(nxp, nyp, nzp)
            incd = incd%magnitude()
            if(iold /= celli)then                           ! x-dir

                norm = vector(1., 0., 0.)
                call reflect_refract(incd, norm, n1, n2, iseed, rflag)
                if(rflag)then
                    celli = iold
                    if(nxp > 0.)then
                        xcur = xface(celli+1) - delta
                    elseif(nxp < 0.)then
                        xcur = xface(celli) + delta
                    end if
                end if
            elseif(jold /= cellj)then                       ! y-dir

                norm = vector(0., 1., 0.)
                call reflect_refract(incd, norm, n1, n2, iseed, rflag)
                if(rflag)then
                    cellj = jold
                    if(nyp > 0.)then
                        ycur = yface(cellj+1) - delta
                    elseif(nyp < 0.)then
                        ycur = yface(cellj) + delta
                    end if
                end if
            elseif(kold /= cellk)then                       ! z-dir
                norm = vector(0., 0., 1.)
                call reflect_refract(incd, norm, n1, n2, iseed, rflag)
                if(rflag)then
                    cellk = kold
                    if(nzp > 0.)then
                        zcur = zface(cellk+1) - delta
                    elseif(nzp < 0.)then
                        zcur = zface(cellk) + delta
                    end if
                elseif(cellk == nzg + 1)then
                    cellk = -1
                    tflag = .true.
                end if
            else
                print*,'Error in reflect/refract in update_pos!'
                error stop 
            end if
                nxp = incd%x
                nyp = incd%y
                nzp = incd%z

                phi = atan2(nyp, nxp)
                sinp = sin(phi)
                cosp = cos(phi)

                cost = nzp
                sint = sqrt(1.-cost*cost)
        end if
    end if
end if
    end subroutine update_pos


    subroutine update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        use iarray,    only : xface, yface, zface

        implicit none

        real,    intent(IN)    :: xcur, ycur, zcur
        integer, intent(INOUT) :: celli, cellj, cellk

        celli = find(xcur, xface) 
        cellj = find(ycur, yface)
        cellk = find(zcur, zface) 

    end subroutine update_voxels


    integer function find(val, a)
    !searchs for bracketing indicies for a value val in an array a
    !
    !
        implicit none

        real, intent(IN) :: val, a(:)
        integer          :: n, lo, mid, hi

        n = size(a)
        lo = 0
        hi = n + 1

        if (val == a(1)) then
            find = 1
        else if (val == a(n)) then
            find = n-1
        else if((val > a(n)) .or. (val < a(1))) then
            find = -1
        else
            do
                if (hi-lo <= 1) exit
                mid = (hi+lo)/2
                if (val >= a(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            find = lo
        end if
    end function find


    subroutine repeat_bounds(cella, cellb, acur, bcur, amax, bmax, nag, nbg, delta) 
    !   if photon leaves grid in a direction a or b, then photon is transported to otherside and continues being simulated
    !
    !
        implicit none

        real,    intent(INOUT) :: acur, bcur
        real,    intent(IN)    :: delta, amax, bmax
        integer, intent(IN)    :: nag, nbg
        integer, intent(INOUT) :: cella, cellb

        if(cella == -1)then
            if(acur < delta)then
                acur = 2.*amax  - delta
                cella = nag
            elseif(acur > 2. * amax - delta)then
                acur = delta
                cella = 1
            else
                print*,'Error in Repeat_bounds...'
                error stop 
            end if
        end if
        if(cellb == -1)then
            if(bcur < delta)then
                bcur = 2.*bmax-delta
                cellb = nbg
            elseif(bcur > 2.*bmax-delta)then
                bcur = delta
                cellb = 1
            else
                print*,'Error in Repeat_bounds...'
                error stop 
            end if
        end if
    end subroutine repeat_bounds
   

    subroutine reflect_refract(I, N, n1, n2, iseed, rflag)

        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(INOUT) :: N
        real,         intent(IN)    :: n1, n2
        integer,      intent(INOUT) :: iseed
        logical,      intent(OUT)   :: rflag

        real :: ran2

        rflag = .FALSE.

        if(ran2(iseed) <= fresnel(I, N, n1, n2))then
            call reflect(I, N)
            rflag = .true.
        else
            call refract(I, N, n1/n2)
        end if

    end subroutine reflect_refract


    subroutine reflect(I, N)
    !   get vector of reflected photon
    !
    !
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2. * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
    !   get vector of refracted photon
    !
    !
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N
        real,         intent(IN)    :: eta

        type(vector) :: T, Ntmp

        real :: c1, c2

        Ntmp = N

        c1 = (Ntmp .dot. I)
        if(c1 < 0.)then
            c1 = -c1
        else
            Ntmp = (-1.) * N
        end if
        c2 = sqrt(1. - (eta)**2 * (1.-c1**2))

        T = eta*I + (eta * c1 - c2) * Ntmp 

        I = T

    end subroutine refract


    function fresnel(I, N, n1, n2) result (tir)
    !   calculates the fresnel coefficents
    !
    !
        use vector_class

        implicit none

        real, intent(IN)         :: n1, n2
        type(vector), intent(IN) :: I, N

        real             ::  costt, sintt, sint2, cost2, tir, f1, f2

        costt = abs(I .dot. N)

        sintt = sqrt(1. - costt * costt)
        sint2 = n1/n2 * sintt
        if(sint2 > 1.)then
            tir = 1.0
            return
        elseif(costt == 1.)then
            tir = 0.
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1. - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5 * (f1 + f2)
        if(isnan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir, f1, f2, costt,sintt,cost2,sint2
            return
        end if
    end function fresnel


    subroutine taufind1(xcell,ycell,zcell,delta,taurun)
    !   routine to find tau from current position to edge of grid in a direction (nxp,nyp,nzp)
    !
    !
        use photon_vars, only : xp, yp, zp
        use iarray,      only : rhokap
        use opt_prop,    only : wavelength, material
        use constants,   only : xmax, ymax, zmax 
     
        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell

        real                   :: taurun, taucell, xcur, ycur, zcur, d, dcell, tau
        integer                :: celli, cellj, cellk, iseed
        logical                :: dir(3),tmp

        xcur = xp + xmax
        ycur = yp + ymax
        zcur = zp + zmax

        celli = xcell
        cellj = ycell
        cellk = zcell

        taurun = 0.
        taucell = 0.

        d = 0.
        dcell = 0.

        dir = (/.FALSE., .FALSE., .FALSE./)
        do
            dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
            taucell = dcell * rhokap(celli,cellj,cellk, wavelength+material)

            taurun = taurun + taucell
            d = d + dcell
            call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta, tmp, iseed, &
                            tau, taurun,.false.)

            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                exit
            end if
        end do
    end subroutine taufind1


    subroutine tauquick(xcell, ycell, zcell, zp, delta, taurun)

        use constants,   only : nzg, zmax
        use iarray,      only : refrac, rhokap, zface
        use opt_prop,    only : material, wavelength
        use photon_vars, only : nzp

        implicit none


        integer, intent(IN)  :: xcell, ycell, zcell
        real,    intent(IN)  :: zp, delta
        real,    intent(OUT) :: taurun

        integer :: cellk
        real    :: n1, tau, dcell, taucell, zcur

        taurun = 0.

        zcur = zp + zmax

        do cellk = zcell, nzg

            dcell = (zface(cellk+1) - zcur)/nzp
            taucell = dcell * rhokap(xcell, ycell, cellk, wavelength + material)
            taurun = taurun + taucell
            zcur = zface(cellk+1) + delta
            n1 = refrac(xcell, ycell, cellk)
            if(n1 == 1.38)then
                material = 1
            else
                material = 3
            end if
        end do
! error stop 
    end subroutine tauquick
   

    subroutine peeling(xcell,ycell,zcell,delta,flag)
   
        use iarray,      only : image
        use constants,   only : PI, nbins, xmax, ymax, zmax, v, costim, sintim, cospim, sinpim
        use photon_vars, only : xp, yp, zp, nxp, nyp, nzp
        use opt_prop,    only : hgg, g2, material, wavelength
        use taufind2

        implicit none


        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell
        logical, intent(IN)    :: flag
        real                   :: cosa, tau1, prob, xim, yim, bin_wid,xpold,ypold
        real                   :: nxpold, nypold, nzpold, hgfact,zpold, tau3
        integer                :: binx, biny, xcellold, ycellold, zcellold, matold, wavold


        ! print'(6(F9.5,1X),5(I3.3,1X))',nxp,nyp,nzp,xp,yp,zp,xcell,ycell,zcell,material,wavelength

        nxpold = nxp
        nypold = nyp
        nzpold = nzp

        xcellold = xcell
        ycellold = ycell
        zcellold = zcell

        matold = material
        wavold = wavelength

        xpold = xp
        ypold = yp
        zpold = zp
        cosa = nxp*v(1) + nyp*v(2) + nzp*v(3)!angle of peeled off photon

        nxp = v(1)
        nyp = v(2)
        nzp = v(3)
        xim = yp*cospim - xp*sinpim
        yim = zp*sintim - yp*costim*sinpim - xp*costim*cospim

        ! call taufind1(xcell, ycell, zcell, delta, tau3)
        call tauquick(xcell, ycell, zcell, zp, delta, tau1)
        ! call tau2(xcell,ycell,zcell,delta,tau3)
        ! if(tau3 /= tau1)print*,tau3, tau1
        ! prob = exp(-tau3)

        bin_wid = 4.*xmax/Nbins

        binx = floor(xim/bin_wid)
        biny = floor(yim/bin_wid)

        if(flag)then
            hgfact = 1./(4.*pi)
        else
            hgfact = (1.-g2) / ((4.*pi)*(1.+g2-2.*hgg*cosa)**(1.5))
        end if

        ! print*,prob*hgfact,prob,hgfact
        prob = hgfact * exp(-tau1)
        ! if(material==3)print*,material+wavelength
                material = matold
        wavelength = wavold
        image(binx, biny,material + wavelength) = image(binx, biny, material + wavelength) + prob

        nxp = nxpold
        nyp = nypold
        nzp = nzpold

        xcell = xcellold 
        ycell = ycellold 
        zcell = zcellold 

        material = matold
        wavelength = wavold
        ! call update_voxels(xp+xmax, yp+ymax, zp+zmax, xcell, ycell, zcell)
! print'(6(F9.5,1X),5(I3.3,1X))',nxp,nyp,nzp,xp,yp,zp,xcell,ycell,zcell,material,wavelength
!         print*,
    end subroutine peeling
end module inttau2