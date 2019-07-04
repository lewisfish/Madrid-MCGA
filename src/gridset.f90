MODULE gridset_mod

    implicit none
    save

    private
    public :: gridset

    contains
        subroutine gridset(id, height)

            use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
            use iarray, only    : rhokap,xface,yface,zface, rhokap,albedo_a, refrac
            use opt_prop, only  : kappa, albedo
            use ch_opt

            implicit none

            integer, intent(IN) :: id
            real,    intent(IN) :: height

            integer             :: i, j, k
            real                :: x, y, z, taueq1, taupole1, taueq2, taupole2

            if(id == 0)then
                print*, ' '
                print *, 'Setting up density grid....'
            end if

            ! setup grid faces
            do i = 1, nxg + 1
                xface(i) = (i - 1) * 2. * xmax/nxg
            end do

            do i = 1, nyg + 1
                yface(i) = (i - 1) * 2. * ymax/nyg
            end do

            do i = 1, nzg + 1
                zface(i) = (i - 1) * 2. * zmax/nzg
            end do

            call init_opt1
            refrac = 1.38
            refrac(:,:,nzg+1) = 1.0


            !set up optical properties grid 
            do i = 1, nxg
                x = xface(i) - xmax + xmax/nxg
                do j = 1, nyg
                    y = yface(j) - ymax + ymax/nyg
                    do k = 1, nzg
                        z = zface(k) - zmax + zmax/nzg

                        if(x >= -.33 .and. x <= 0.33 .and. y >= -0.51 .and. y <= .51 .and. z >= (-zmax+.013d0) .and.&
                         z <= ((-zmax + .013d0)+.29))then
                            refrac(i,j,k) = 1.8                        !crystal

                            call init_opt3                               !809 CRYSTAL         
                            albedo_a(i,j,k,3) = albedo
                            rhokap(i,j,k,3)   = kappa

                            call init_opt4                               !1064 CRYSTAL         
                            albedo_a(i,j,k,4) = albedo
                            rhokap(i,j,k,4)   = kappa
                        elseif(x <= (-xmax+.013d0) .or. x >= (xmax - .013d0)  .or. y <= (-ymax+.013d0) .or. y >= (ymax - .013d0)&
                         .or. z <= (-zmax+.013d0))then
                            refrac(i,j,k) = 1.5d0
                            rhokap(i,j,k,1:2)   = 0.000000001d0        
                            albedo_a(i,j,k,1:2) = 0.000000001d0
                        elseif(z <= ((-zmax + .013d0)+.29+height))then
                            refrac(i,j,k) = 1.38                         !tissue

                            call init_intra(809.d0)                      !809 TISSUE              
                            rhokap(i,j,k,1)   = kappa        
                            albedo_a(i,j,k,1) = albedo

                            call init_intra(1064.d0)                     !1064 TISSUE         
                            rhokap(i,j,k,2)   = kappa        
                            albedo_a(i,j,k,2) = albedo
                        else!air
                            refrac(i,j,k)=1.0d0
                            rhokap(i,j,k,:) = 0.000000001d0
                            albedo_a(i,j,k,:) = 0.000000001d0
                        end if
                    end do
                end do
            end do

            !****************** Calculate equatorial and polar optical depths ****
            taueq1   = 0.
            taupole1 = 0.
            taueq2   = 0.
            taupole2 = 0.

            do i = 1, nxg
                taueq1 = taueq1 + rhokap(i,nyg/2,nzg/2,1)
            end do

            do i = 1, nzg
                taupole1 = taupole1 + rhokap(nxg/2,nyg/2,i,1)
            end do

            taueq1 = taueq1 * 2. * xmax/nxg
            taupole1 = taupole1 * 2. * zmax/nzg
            if(id == 0)then
                print'(A,F9.5,A,F9.5)',' taueq1 = ',taueq1,'  taupole1 = ',taupole1
            end if

            if(id == 0)then
                print*,(2.d0*xmax) / nxg
                print*,(2.d0*ymax) / nyg
                print*,(2.d0*zmax) / nzg

                open(newunit=j,file='refrac.dat',access='stream',form='unformatted',status='replace')
                write(j)refrac
                close(j)
                open(newunit=j,file='rhokap.dat',access='stream',form='unformatted',status='replace')
                write(j)rhokap
                close(j)

                inquire(iolength=i)albedo_a
                open(newunit=j,file='albedo.dat',access='stream',form='unformatted',status='replace')
                write(j)albedo_a
                close(j)
            end if
            stop
        end subroutine gridset
end MODULE gridset_mod
