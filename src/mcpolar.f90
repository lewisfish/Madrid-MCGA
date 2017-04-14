program mcpolar

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

integer           :: nphotons, iseed, j, xcell, ycell, zcell
logical           :: tflag
double precision  :: nscatt
real              :: xmax, ymax, zmax, ran, delta
real              :: start, finish, ran2, depth


integer           :: id, error, numproc, i
real              :: nscattGLOBAL

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray


call MPI_init(error)

call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
   read(10,*) nphotons
   read(10,*) xmax
   read(10,*) ymax
   read(10,*) zmax
   read(10,*) n1
   read(10,*) n2
   close(10)

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


! do i = 1, 10
! call zarray()
depth = zmax - .5!real(i)/10.
if(id == 1)print*,depth
!***** Set up density grid *******************************************
call gridset(xmax,ymax,zmax,id,depth)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
   delta = 1.e-8*(2.*zmax/nzg)
!   deltay = 1.e-5*(2.*ymax/nyg)          !!!!!!!!1! impliment!!!!!!!!!!!!!!!!!!!!!!!!!
!   deltaz = 1.e-5*(2.*zmax/nzg)
nscatt=0

call cpu_time(start)

!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',id
do j=1,nphotons
  
   ! call init_opt4
   wavelength = 0
   material = 1

   tflag=.FALSE.

   if(mod(j,10000) == 0)then
      print *, j,' scattered photons completed on core: ',id
   end if
    
!***** Release photon from point source *******************************
   call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed,j)

!****** Find scattering location

      call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
    if(wavelength /= 0)then
        print*,wavelength,material
        call peeling(xmax,ymax,zmax,xcell,ycell,zcell,delta)
    end if
!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.)
      ran = ran2(iseed)
      ! if(material/=)print*,material,wavelength
      if(ran < albedo_a(xcell,ycell,zcell,wavelength +material))then!interacts with tissue
            call stokes(iseed)
            nscatt = nscatt + 1
         else
            ! if(material==3)print*,'here',material
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
        call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)

        if(wavelength /= 0 .and. .not. tflag)then
            call peeling(xmax,ymax,zmax,xcell,ycell,zcell,delta)
        end if
   end do
end do      ! end loop over nph photons


call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg)*4,MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_REDUCE(image,imageGLOBAL,size(image,1)*size(image,2)*size(image,3),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

if(id == 0)then
   print*,'Average # of scatters per photon:',nscattGLOBAL/(nphotons*numproc)
   !write out files
   call writer(xmax,ymax,zmax,nphotons, numproc,depth)
   print*,'write done'
end if
! end do
call MPI_Finalize(error)
end program mcpolar
