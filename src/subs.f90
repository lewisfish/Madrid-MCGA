MODULE subs

    implicit none
    save

    contains
        
        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir = trim(homedir)//'res/'

        end subroutine directory


        subroutine zarray
        !   sets all arrays to zero
        !
        !
            use iarray

            implicit none

            xface = 0.
            yface = 0.
            zface = 0.
            rhokap = 0.
            albedo_a = 0.
            refrac = 1.
            image = 0.
            imageGLOBAL = 0.
            jmean = 0.
            jmeanGLOBAL = 0.
            absorb = 0.
            absorbGLOBAL = 0.
        end subroutine zarray


        subroutine alloc_array
        !   subroutine allocates allocatable arrays
        !   
        !   
            use iarray
            use constants, only : nxg, nyg, nzg, nbins

            implicit none

            allocate(xface(nxg+1), yface(nyg+1), zface(nzg+1))
            allocate(rhokap(nxg,nyg,nzg,4), albedo_a(nxg,nyg,nzg,4))
            allocate(refrac(nxg,nyg,nzg+1))
            allocate(image(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2), 4), &
            imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2), -((Nbins-1)/2):((Nbins-1)/2),4))
            allocate(jmean(nxg,nyg,nzg,4), jmeanGLOBAL(nxg,nyg,nzg,4))
            allocate(absorb(nxg, nyg, nzg), absorbGLOBAL(nxg, nyg, nzg))
        end subroutine alloc_array
end MODULE subs
