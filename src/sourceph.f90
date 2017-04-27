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
end MODULE sourceph_mod
