module MCf

implicit none

contains


subroutine mcpolar(depth, id, numproc)

use mpi

!shared data
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use gridset_mod
use sourceph_mod
use inttau2
use ch_opt
use stokes_mod
use writer_mod

implicit none

integer, intent(IN) :: id, numproc
real,    intent(IN) :: depth

integer           :: nphotons, iseed, j, xcell, ycell, zcell
logical           :: tflag
double precision  :: nscatt
real              :: ran, delta, ran2, phiim, thetaim

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
   read(10,*) n1
   read(10,*) n2
   close(10)

   zmax = depth

! set seed for rnd generator. id to change seed for each process
iseed=-456123456+id

!****** setup up arrays and bin numbers/dimensions

!***** Set up constants, pi and 2*pi  ********************************

iseed=-abs(iseed)  ! Random number seed must be negative for ran2

call init_opt4

if(id == 0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(id, depth)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
   delta = 1.e-8*(2.*zmax/nzg)

!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',id
do j=1,nphotons
  
   ! call init_opt4
   wavelength = 0
   material = 1

   tflag = .FALSE.

   if(mod(j,10000) == 0)then
      print *, j,' scattered photons completed on core: ',id
   end if
    
!***** Release photon from point source *******************************
   call sourceph(xcell,ycell,zcell,iseed,j)

!****** Find scattering location

    call tauint1(xcell,ycell,zcell,tflag,iseed,delta)
    if(wavelength /= 0)then
        call peeling(xcell,ycell,zcell,delta)
    end if
!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.)
      ran = ran2(iseed)
      if(ran < albedo_a(xcell,ycell,zcell,wavelength +material))then!interacts with tissue
            call stokes(iseed)
            nscatt = nscatt + 1
         else
            if(material==3)then
               wavelength = 1
               hgg = 0.
               call stokes(iseed)
               hgg = .9
            else
               tflag=.true.
               exit
            end if
      end if



!************ Find next scattering location
        call tauint1(xcell,ycell,zcell,tflag,iseed,delta)

        if(wavelength /=0 .and. .not. tflag)then
            ! if(material==3)print*,material,tflag
            call peeling(xcell,ycell,zcell,delta)
        end if
   end do
end do      ! end loop over nph photons

call MPI_REDUCE(image,imageGLOBAL,size(image,1)*size(image,2)*size(image,3),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

if(id == 0)then
   call writer(depth)
   print*,'write done'
end if
end subroutine mcpolar
end module MCf