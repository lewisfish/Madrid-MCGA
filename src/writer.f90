MODULE writer_mod

implicit none
save

CONTAINS
    subroutine writer(nphotons,numproc,depth,dflag)


        use constants, only : fileplace, xmax, ymax, zmax, nxg, nyg, nzg, Nbins
        use iarray,    only : jmeanGLOBAL, absorbGLOBAL, imageGLOBAL

        implicit none

        integer, intent(IN) :: numproc, nphotons
        real,    intent(IN) :: depth
        logical, intent(IN) :: dflag

        character(len=10) :: fn
        integer :: u, j, i

        write(fn,'(F4.2)') (2.*depth-.3)

        if(dflag)then
        ! jmeanGLOBAL = jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        ! open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat', access='stream',form='unformatted',status='replace')
        ! write(u)jmeanGLOBAL
        ! close(u)
        open(newunit=u,file=trim(fileplace)//'jmean/1064.dat', access='stream',form='unformatted',status='replace')
        write(u)jmeanGLOBAL(:,:,:,2)+jmeanGLOBAL(:,:,:,4)
        close(u)

        open(newunit=u,file=trim(fileplace)//'im/image-1064-'//trim(fn)//'.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) (imageGLOBAL(j,i,2)+imageGLOBAL(j,i,4),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)

        open(newunit=u,file=trim(fileplace)//'im/prob-scatt-slice-1064.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*)imageGLOBAL(i,0,2)+imageGLOBAL(i,0,4)
        end do
        
        ! open(newunit=u,file=trim(fileplace)//'im/plot.gp',status='replace')
        ! write(u,*) 'set terminal pngcairo size 1920,1080'
        ! write(u,*) 'set output "prob-scatt-1064-'//trim(fn)//'.png"'
        ! write(u,*) 'p"image-1064-'//trim(fn)//'.dat" matrix w image'
        ! close(u)
! #ifdef intel
!         i = chdir(trim(fileplace)//'im/')
! #else
!         call chdir(trim(fileplace)//'im/',i)
! #endif  
        ! call system('gnuplot plot.gp')

        else
        absorbGLOBAL = absorbGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        open(newunit=u,file=trim(fileplace)//'jmean/absorb.dat', access='stream',form='unformatted',status='replace')
        write(u)absorbGLOBAL
        close(u)

        ! jmeanGLOBAL = jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        ! open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat', access='stream',form='unformatted',status='replace')
        ! write(u)jmeanGLOBAL
        ! close(u)
        end if


    end subroutine writer
end MODULE writer_mod
