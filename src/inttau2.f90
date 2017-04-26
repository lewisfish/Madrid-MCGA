module inttau2

   implicit none
   
   private
   public :: tauint1, peeling

CONTAINS

    subroutine tauint1(xcell,ycell,zcell,tflag,iseed,delta)
    !optical depth integration subroutine
    !
    !
        use constants,   only : nxg, nyg, nzg, xmax, ymax, zmax
        use photon_vars, only : xp, yp, zp, nxp, nyp, nzp, cost, sint, cosp, sinp, phi
        use iarray,      only : rhokap, jmean
        use opt_prop,    only : wavelength, material, kappa
        use vector_class
   
        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell, iseed
        logical, intent(INOUT) :: tflag

        type(vector)           :: incd, norm
        real                   :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell, ran2
        integer                :: celli, cellj, cellk
        logical                :: dir(3), rflag

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
            taucell = dcell * kappa!rhokap(celli,cellj,cellk,1)

            if(taurun + taucell < tau)then
                taurun = taurun + taucell
                d = d + dcell
                jmean(celli,cellj,cellk,1)=jmean(celli,cellj,cellk,1) + dcell
                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)

            else

                dcell = (tau - taurun) / kappa!rhokap(celli,cellj,cellk,1)
                d = d + dcell
                jmean(celli,cellj,cellk,1)=jmean(celli,cellj,cellk,1) + dcell

                call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, delta)
                exit
            end if
         
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                if(celli == -1 .or. cellk == -1)then
                    call repeat_bounds(celli, cellk, xcur, zcur, xmax, zmax, nxg, nzg, delta)
                    if(celli == -1 .or. cellk == -1 .or. tflag)then
                       print*,'error',celli,cellj,tflag
                    end if
                elseif(cellj == -1)then
                    if(ycur >= 2.*ymax-delta)then
                        ycur = 2.*ymax
                        cellj = nyg
                        rflag = .false.
                        incd = vector(nxp,nyp,nzp)
                        norm = vector(0.,-1,0.)
                        ! print*,incd
                        incd = incd%magnitude()
                        ! print*,incd
                        ! call exit(0)
                        if(ran2(iseed) <= fresnel(incd, norm ,1.38, 1.0))then
                            call incd%print()
                            call reflect(incd, norm)
                            call incd%print()

                            print*,
                            nxp = incd%x
                            nyp = incd%y
                            nzp = incd%z
                            phi = atan2(nyp, nzp)
                            sinp = sin(phi)
                            cosp = cos(phi)
                            cost = nzp
                            sint = 1.-cost*cost
                        ! call reflect(zcur, zmax, nzp, cost, sint, 1.38, 1.0, iseed, delta, rflag)
                        ! if(rflag)then
                        !     phi = atan2(sinp, cosp)
                        !     nxp = sint * cosp
                        !     nyp = sint * sinp
                        !     nzp = cost
                            taurun = 0.
                            tau = -log(ran2(iseed))
                        else
                            tflag = .true.
                            exit
                        end if
                    else
                        tflag = .true.
                        exit
                    end if
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
        ! if(wall_dist < 0.)print'(A,7F9.5)','dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
        if(wall_dist == dx)dir=(/.TRUE., .FALSE., .FALSE./)
        if(wall_dist == dy)dir=(/.FALSE., .TRUE., .FALSE./)
        if(wall_dist == dz)dir=(/.FALSE., .FALSE., .TRUE./)
        if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, wall_flag, dir, delta, peelflag)
    !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
    !
    !
        use photon_vars, only : nxp, nyp, nzp
        use opt_prop,    only : material
        use iarray,      only : xface, yface, zface, refrac

        implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur
      real,    intent(IN)    :: dcell, delta
      integer, intent(INOUT) :: celli, cellj, cellk
      logical, intent(IN)    :: wall_flag, dir(:)   
      logical, optional, intent(IN) :: peelflag
      
      real :: n1

      n1 = refrac(celli,cellj,cellk)
      
      if(wall_flag)then
      
         if(dir(1))then
            if(nxp > 0.)then
               xcur = xface(celli+1) + delta
            elseif(nxp < 0.)then
               xcur = xface(celli) - delta
            else
               print*,'Error in x dir in update_pos', dir, nxp, nyp, nzp
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
            end if
         else
            print*,'Error in update_pos...',dir
            call exit(0)
         end if
      else
      
         xcur = xcur + nxp*dcell
         ycur = ycur + nyp*dcell 
         zcur = zcur + nzp*dcell
      
      end if
     if(present(peelflag))then
        if(.not.peelflag)then
            if(n1 == 1.38)then
            material = 3
            else
            material = 1
            end if
        end if
    else
        if(n1 == 1.38)then
            material = 3
            else
            material = 1
        end if
    end if

      if(wall_flag)then
         call update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
      end if
      
    end subroutine update_pos


    subroutine update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
    !updates the current voxel based upon position
    !
    !
        use iarray,    only : xface, yface, zface
        use constants, only : nxg, nyg, nzg, xmax, ymax, zmax

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
    !if photon leaves grid in a direction a or b, then photon is transported to otherside and continues being simulated
    !
    !
        implicit none

        real,    intent(INOUT) :: acur, bcur
        real,    intent(IN)    :: delta, amax, bmax
        integer, intent(IN)    :: nag, nbg
        integer, intent(INOUT) :: cella, cellb

        if(cella == -1)then
            if(acur < delta)then
                acur = 2.*amax  -delta
                cella = nag
            elseif(acur > 2.*amax-delta)then
                acur = delta
                cella = 1
            else
                print*,'Error in Repeat_bounds...'
                call exit(0)
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
                call exit(0)
            end if
        ! else
        ! tflag=.true.
        end if
    end subroutine repeat_bounds
   

    ! subroutine reflect(pcur, pmax, pdir, pdir_c, pdir_s, n1, n2, iseed, delta, rflag)
    ! !carries out fresnel reflection
    ! !
    ! !
    !     implicit none

    !     real,    intent(IN)    :: n1, n2, pcur, delta, pdir, pmax
    !     real,    intent(INOUT) :: pdir_c, pdir_s
    !     integer, intent(INOUT) :: iseed
    !     logical, intent(INOUT) :: rflag
    !     real                   :: ran2, tmp

    !     rflag = .false.
    !     tmp = pdir_s
    !     if(ran2(iseed) <= fresnel(pdir, n1, n2))then
    !         rflag = .true.
    !         pdir_c = -pdir_c
    !         pdir_s = (1. - pdir_c*pdir_c)
    !         if(pdir_s < 0.)then
    !             pdir_s = 0.
    !         else
    !             pdir_s = (sqrt(pdir_s))
    !         end if 
    !     end if
    ! end subroutine reflect


    subroutine reflect(I, N)

        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N

        type(vector) :: R

        R = I - 2. * (N .dot. I) * N
        I = R

    end subroutine reflect


    subroutine refract(I, N, eta)
        
        use vector_class

        implicit none

        type(vector), intent(INOUT) :: I
        type(vector), intent(IN)    :: N
        real,         intent(IN)    :: eta

        type(vector) :: T

        real :: c1, c2

        c1 = N .dot. I ! or cos(theta_1)
        c2 = sqrt(1. - (eta)**2 * (1.-c1**2))

        T = eta + (eta * c1 - c2) * N 

    end subroutine refract


    function fresnel(I, N, n1, n2) result (tir)
    !calculates the fresnel coefficents
    !
    !
        use vector_class

        implicit none

        real, intent(IN)         :: n1, n2
        type(vector), intent(IN) :: I, N
        real             :: crit, costt, sintt, sint2, cost2, tir, f1, f2

        crit = n2/n1

        costt = N .dot. I
        if(costt < 0.)then
            print*,costt
            call exit(0)
        end if
        sintt = sqrt(1. - costt * costt)

        if(sintt > crit)then
            tir = 1.0
            return
        else
            sint2 = (n1/n2)*sintt
            cost2 = sqrt(1. - sint2 * sint2)
            f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
            f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

            tir = 0.5 * (f1 + f2)
        if(isnan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir!, f1, f2, cost,sint,cost,sint2
            return
        end if
   
    end function fresnel

   
    ! function fresnel(pdir, n1, n2) result (tir)
    ! !calculates the fresnel coefficents
    ! !
    ! !
    !     implicit none

    !     real, intent(IN) :: n1, n2, pdir
    !     real             :: crit, costt, sintt, sint2, cost2, tir, f1, f2

    !     crit = n2/n1

    !     costt = abs(pdir)
    !     sintt = sqrt(1. - costt * costt)

    !     if(sintt > crit)then
    !         tir = 1.0
    !         return
    !     else
    !         sint2 = (n1/n2)*sintt
    !         cost2 = sqrt(1. - sint2 * sint2)
    !         f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
    !         f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

    !         tir = 0.5 * (f1 + f2)
    !     if(isnan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir!, f1, f2, cost,sint,cost,sint2
    !         return
    !     end if
   
    ! end function fresnel


    ! function fresnel(pdir, n1, n2) result (tir)
    ! !calculates the fresnel coefficents
    ! !
    ! !
    !     implicit none

    !     real, intent(IN) :: n1, n2, pdir
    !     real             :: crit, costt, sintt, sint2, cost2, tir, f1, f2

    !     crit = n2/n1

    !     costt = abs(pdir)
    !     sintt = sqrt(1. - costt * costt)

    !     if(sintt > crit)then
    !         tir = 1.0
    !         return
    !     else
    !         sint2 = (n1/n2)*sintt
    !         cost2 = sqrt(1. - sint2 * sint2)
    !         f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
    !         f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

    !         tir = 0.5 * (f1 + f2)
    !     if(isnan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir!, f1, f2, cost,sint,cost,sint2
    !         return
    !     end if
   
    ! end function fresnel


    subroutine taufind1(xcell,ycell,zcell,delta,taurun)
    !routine to find tau from current position to edge of grid in a direction (nxp,nyp,nzp)
    !
    !
        use photon_vars, only : xp, yp, zp
        use iarray,      only : rhokap
        use opt_prop,    only : wavelength, material
        use constants,   only : xmax, ymax, zmax 
     
        implicit none

        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell

        real                   :: taurun, taucell, xcur, ycur, zcur, d, dcell
        integer                :: celli, cellj, cellk
        logical                :: dir(3)

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
            call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta, .TRUE.)

            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                exit
            end if
        end do
    end subroutine taufind1
   

    subroutine peeling(xcell,ycell,zcell,delta)
   
        use iarray,      only : image
        use constants,   only : PI, nbins, xmax, ymax, zmax, v, costim, sintim, cospim, sinpim
        use photon_vars, only : xp, yp, zp, nxp, nyp, nzp
        use opt_prop,    only : hgg, g2, material, wavelength
        use taufind2

        implicit none


        real,    intent(IN)    :: delta
        integer, intent(INOUT) :: xcell, ycell, zcell
        real                   :: cosa, tau1, prob, xim, yim, bin_wid,xpold,ypold
        real                   :: nxpold, nypold, nzpold, hgfact,zpold, tau3
        integer                :: binx, biny, xcellold, ycellold, zcellold

        nxpold = nxp
        nypold = nyp
        nzpold = nzp

        xcellold = xcell
        ycellold = ycell
        zcellold = zcell

        xpold = xp
        ypold = yp
        zpold = zp
        cosa = nxp*v(1) + nyp*v(2) + nzp*v(3)!angle of peeled off photon

        nxp = v(1)
        nyp = v(2)
        nzp = v(3)
        xim = yp*cospim - xp*sinpim
        yim = zp*sintim - yp*costim*sinpim - xp*costim*cospim

        ! call taufind1(xcell,ycell,zcell,delta,tau3)
        call tau2(xcell,ycell,zcell,delta,tau3)


        !over 4 pi for 1st fluro photon
        prob = exp(-tau1)

        bin_wid = 4.*xmax/Nbins

        binx = floor(xim/bin_wid)
        biny = floor(yim/bin_wid)

        hgfact=(1.-g2)/(4.*pi*(1.+g2-2.*hgg*cosa)**(1.5))

        prob = exp(-tau3) * hgfact
        ! print*,hgfact,prob

        ! if(material==3)print*,material+wavelength
        image(binx, biny,material + wavelength) = image(binx, biny, material + wavelength) + prob

        nxp = nxpold
        nyp = nypold
        nzp = nzpold

        xcell = xcellold 
        ycell = ycellold 
        zcell = zcellold 

        ! call update_voxels(xp+xmax, yp+ymax, zp+zmax, xcell, ycell, zcell)
    
    end subroutine peeling

end module inttau2