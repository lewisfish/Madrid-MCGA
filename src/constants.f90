MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.
!

implicit none
save

integer, parameter :: nxg =200, nyg=200, nzg=200, nbins=200
real,    parameter :: PI=4.*atan(1.), TWOPI=8.*atan(1.), OFFSET=1.e-2*(2.*.5/nxg)
real               :: xmax, ymax, zmax, v(3), costim, sintim, cospim, sinpim
character(len=255) :: cwd, homedir, fileplace, resdir

end MODULE constants
