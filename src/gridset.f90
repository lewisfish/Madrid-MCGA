MODULE gridset_mod

implicit none
save

private
public :: gridset

CONTAINS
    subroutine gridset(id, depth)

    use constants, only : nxg, nyg, nzg, xmax, ymax, zmax
    use iarray, only    : rhokap,xface,yface,zface, rhokap,albedo_a, refrac
    use opt_prop, only  : kappa, albedo
    use ch_opt

    implicit none

    integer, intent(IN) :: id
    real,    intent(IN) :: depth

    integer             :: i, j, k, u
    real                :: x, y, z, taueq1, taupole1, taueq2, taupole2

    if(id == 0)then
        print*, ' '
        print *, 'Setting up density grid....'
    end if
   !**********  Linear Cartesian grid. Set up grid faces ****************
    do i=1,nxg+1
        xface(i) = (i - 1) * 2. * xmax/nxg
    end do

    do i=1,nyg+1
        yface(i) = (i - 1) * 2. * ymax/nyg
    end do

    do i=1,nzg+1
        zface(i) = (i - 1) * 2. * zmax/nzg
    end do

   call init_opt1
   !**************  Loop through x, y, and z to set up grid density and refractive index grid.  ****
    do i = 1, nxg
        x = xface(i) - xmax + xmax/nxg
        do j = 1, nyg
            y = yface(j) - ymax + ymax/nyg
            do k = 1, nzg
                z = zface(k) - zmax + zmax/nzg
!***********Call density setup subroutine 
                        rhokap(i,j,k,1)   = kappa
                        albedo_a(i,j,k,1) = albedo
                        ! if(x >= -.2 .and. x <= 0.2 .and. y >= -0.2 .and. y <= .2 .and. z >= -0.2 .and. z <= 0.2)then
                        !     refrac(i,j,k) = 1.5
                    ! if(x >= -.3 .and. x <= 0.3 .and. y >= -0.5 .and. y <= .5 .and. z >= -zmax .and. z <= (-zmax+.3))then
                    !     refrac(i,j,k)=1.38
                    !     albedo_a(i,j,k,3) = 0.001                   !809nm crystal
                    !     rhokap(i,j,k,3)   = 6.
                    
                    !     albedo_a(i,j,k,4) = 0.001                   !1064nm crystal  mus = 0.001 cm-1, mua = 0.004 cm-1
                    !     rhokap(i,j,k,4)   = 0.02       !900nm crystal   mus = 0.001 cm-1, mua = 0.02 cm-1
                    ! else
                    !     call init_opt1
                    !     rhokap(i,j,k,1)   = kappa        
                    !     albedo_a(i,j,k,1) = albedo
                    ! else
                    !     refrac(i,j,k) = 1.0
                    ! end if
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

    ! inquire(iolength=i)rhokap(:,:,:,3)
    ! open(newunit=u, file='rhokap.dat', access='direct', status='REPLACE',form='unformatted', &
    !  recl=i)
    !  write(u,rec=1) rhokap(:,:,:,3)
    !  close(u)

    end subroutine gridset
end MODULE gridset_mod
