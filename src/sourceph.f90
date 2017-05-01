MODULE sourceph_mod

    implicit none
    save

    contains
        subroutine sourceph(xcell,ycell,zcell,iseed, delta)
        !   generates stataring pos and dir of photon
        !
        !
            use constants, only : nxg, nyg, nzg, twopi, xmax, ymax, zmax
            use photon_vars

            implicit none


            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed
            real,    intent(IN)    :: delta
            
            real :: ran2

            zp = zmax - delta
            xp = xmax * (2. * ran2(iseed) - 1.)
            yp = ymax * (2. * ran2(iseed) - 1.)

            phi  = twopi * ran2(iseed)
            cosp = cos(phi)
            sinp = sin(phi)         
            cost = -1.
            sint = 0.


            nxp = sint * cosp  
            nyp = sint * sinp
            nzp = cost

            xcell = int(nxg * (xp + xmax) / (2. * xmax)) + 1
            ycell = int(nyg * (yp + ymax) / (2. * ymax)) + 1
            zcell = int(nzg * (zp + zmax) / (2. * zmax)) + 1
        end subroutine sourceph


        subroutine source_diffuse(xcell, ycell, zcell, i, j, k, iseed)

            use iarray,      only : xface, yface, zface
            use constants,   only : xmax, ymax, zmax, twopi, nxg, nyg, nzg
            use photon_vars, only : xp, yp, zp, cost, sint, cosp, sinp, phi, nxp, nyp, nzp

            implicit none

            integer, intent(IN)    :: i, j, k
            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed

            real :: ran2

            xp = xface(i) + ran2(iseed) * (xface(i+1) - xface(i)) - xmax
            yp = yface(j) + ran2(iseed) * (yface(j+1) - yface(j)) - ymax
            zp = zface(k) + ran2(iseed) * (zface(k+1) - zface(k)) - zmax

            cost = 2. * ran2(iseed) - 1.
            sint = sqrt(1. - cost**2)

            phi = twopi * ran2(iseed)
            sinp = sin(phi)
            cosp = cos(phi)

            nxp = sint * cosp
            nyp = sint * sinp
            nzp = cost

            xcell = int(nxg*(xp+xmax)/(2.*xmax))+1
            ycell = int(nyg*(yp+ymax)/(2.*ymax))+1
            zcell = int(nzg*(zp+zmax)/(2.*zmax))+1

        end subroutine source_diffuse
end MODULE sourceph_mod
