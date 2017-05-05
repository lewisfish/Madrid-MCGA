module MCf

implicit none

private
public :: mcpolar

contains


subroutine mcpolar(depth, id, numproc, dflag)

use mpi

!shared data
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use gridset_mod
use ch_opt
use writer_mod

implicit none

integer, intent(IN) :: id, numproc
real,    intent(IN) :: depth
logical, intent(IN) :: dflag

integer           :: nphotons, iseed
double precision  :: nscatt, nscattGLOBAL
real              :: delta, phiim, thetaim

integer           :: error

call zarray

phiim   = 0. * pi/180.
thetaim = 0. * pi/180.

!image postion vector
!angle for vector
costim = cos(thetaim)
sintim = sin(thetaim)
sinpim = sin(phiim)
cospim = cos(phiim)

!vector
v(1) = sintim * cospim     
v(2) = sintim * sinpim
v(3) = costim  


!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
   read(10,*) nphotons
   read(10,*) xmax
   read(10,*) ymax
   read(10,*) zmax
   close(10)

   zmax = depth

! set seed for rnd generator. id to change seed for each process
iseed=-456123456+id

!****** setup up arrays and bin numbers/dimensions

!***** Set up constants, pi and 2*pi  ********************************

iseed=-abs(iseed)  ! Random number seed must be negative for ran2

call init_opt1

if(id == 0)then
   print*, ''      
   print*,'# of photons to run',dble(nphotons)*dble(numproc)
end if

!***** Set up density grid *******************************************
call gridset(id)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
   delta = 1.e-9*(2.*zmax/nzg)
!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',id


if(dflag)then
    call diffuse(nphotons, iseed, delta, nscatt, id, dflag)
else
    call emit(nphotons, iseed, delta, nscatt, id, dflag)
end if


call MPI_REDUCE(image,imageGLOBAL,size(image,1)*size(image,2)*size(image,3),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_REDUCE(absorb, absorbGLOBAL, nxg*nyg*nzg, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_REDUCE(jmean, jmeanGLOBAL, nxg*nyg*nzg*4, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_REDUCE(nscatt, nscattGLOBAL, 1, MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

if(id == 0)then
   call writer(nphotons,numproc,depth,dflag)
   print*,'write done',nscattGLOBAL/(numproc*nphotons)
end if
end subroutine mcpolar


subroutine emit(nphotons, iseed, delta, nscatt, id,dflag)

    use opt_prop,    only : material, wavelength
    use iarray,      only : albedo_a
    use sourceph_mod
    use inttau2
    use stokes_mod

    implicit none

    integer, intent(IN)    :: nphotons, id
    integer, intent(INOUT) :: iseed
    real,    intent(IN)    :: delta
    real,    intent(INOUT) :: nscatt
    logical, intent(IN)    :: dflag

    integer :: j, xcell, ycell, zcell
    real    :: ran2, ran
    logical :: tflag

    do j = 1, nphotons

        wavelength = 0
        material = 1

        tflag = .FALSE.

        if(mod(j,100000) == 0)then
            print *, j,' scattered photons completed on core: ',id
        end if

        call sourceph(xcell,ycell,zcell,iseed,delta)


        call tauint1(xcell,ycell,zcell,tflag,iseed,delta,dflag)

        do while(tflag.eqv..FALSE.)
            ran = ran2(iseed)
            if(ran < albedo_a(xcell,ycell,zcell,wavelength+material))then
                call stokes(iseed)
                nscatt = nscatt + 1
            else
                tflag=.true.
            end if

            call tauint1(xcell,ycell,zcell,tflag,iseed,delta,dflag)

        end do
    end do      ! end loop over nph photons
end subroutine emit


subroutine diffuse(nphotons, iseed, delta, nscatt, id, dflag)

    use constants,   only : nxg, nyg, nzg, fileplace
    use opt_prop,    only : wavelength, material
    use iarray,      only : absorb, albedo_a
    use inttau2
    use stokes_mod
    use sourceph_mod

    implicit none

    integer, intent(IN)    :: nphotons, id
    integer, intent(INOUT) :: iseed
    real,    intent(IN)    :: delta
    real,    intent(INOUT) :: nscatt
    logical, intent(IN)    :: dflag

    integer :: i, j, k, p, nph, jcount, xcell, ycell, zcell, u, ncount
    real :: tmp, diff, tot, ran2, ran
    logical :: tflag

    jcount = 0
    ncount = 0
    inquire(iolength=i)absorb
    open(newunit=u, file=trim(fileplace)//'jmean/absorb.dat',status='old',form='unformatted',access='direct',recl=i)
    read(u,rec=1)absorb
    close(u)

    tot = sum(absorb)


    do i = 1, nxg
        do j = 1, nyg
            do k = 1, nzg
                tmp = dble(nphotons) * absorb(i,j,k)/dble(tot)
                diff = tmp - int(tmp)
                if(ran2(iseed) < diff .and. diff > 0)then
                    nph = int(tmp) + 1
                else
                    nph = int(tmp)
                end if
                ncount = ncount + nph
                do p = 1, nph
                    jcount = jcount + 1
                    wavelength = 1
                    material = 3

                    tflag = .FALSE.

                    if(mod(jcount,100000) == 0)then
                        print *, jcount,' scattered photons completed on core: ',id
                    end if

                    call source_diffuse(xcell, ycell, zcell, i, j, k, iseed)
                    call peeling(xcell,ycell,zcell,delta,.true.)

                    call tauint1(xcell,ycell,zcell,tflag,iseed,delta,dflag)
                    if(.not. tflag)then
                        call peeling(xcell,ycell,zcell,delta,.false.)
                    ! else
                    !     exit
                    end if
                    do while(tflag .eqv. .FALSE.)
                        ran = ran2(iseed)
                        if(ran < albedo_a(xcell,ycell,zcell,wavelength+material))then!interacts with tissue
                            call stokes(iseed)
                            nscatt = nscatt + 1
                        else
                            tflag=.true.
                            exit
                        end if
                        call tauint1(xcell,ycell,zcell,tflag,iseed,delta,dflag)
                        if(.not. tflag)then
                            call peeling(xcell,ycell,zcell,delta,.false.)
                        end if
                    end do
                end do
            end do
        end do
    end do
    print*,jcount,ncount
end subroutine diffuse
end module MCf