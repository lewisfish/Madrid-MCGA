MODULE writer_mod

implicit none
save

CONTAINS
    subroutine writer(xmax,ymax,zmax,nphotons, numproc, depth)

        use constants, only : nxg, nyg, nzg, fileplace, nbins
        use iarray,    only : jmeanGLOBAL, imageGLOBAL

        implicit none

        integer           :: nphotons, i, u, numproc, j
        real              :: xmax, ymax, zmax, depth
        character(len=10) :: fn


        !maybe not right
        jmeanGLOBAL =jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))


        inquire(iolength=i)jmeanGLOBAL
        write(fn,'(F6.3)') depth
        open(newunit=u,file=trim(fileplace)//'jmean/jmean-small.dat',access='direct',status='REPLACE',form='unformatted',&
        recl=i)
        write(u,rec=1) jmeanGLOBAL
        close(u)


        open(newunit=u,file=trim(fileplace)//'im/image_small.dat',status='REPLACE')

        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) (imageGLOBAL(j,i,1),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)
    end subroutine writer
end MODULE writer_mod
