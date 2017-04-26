MODULE writer_mod

implicit none
save

CONTAINS
    subroutine writer(nphotons,numproc,depth)


        use constants, only : fileplace, xmax, ymax, zmax, nxg, nyg, nzg!, nbins
        use iarray,    only : jmeanGLOBAL,yface

        implicit none

        integer           :: i, u, j, numproc, nphotons,k
        real              :: depth, tmp(nzg)
        ! character(len=10) :: fn


        ! write(fn,'(F4.2)') (2.*depth-.3)

        inquire(iolength=i)jmeanGLOBAL(:,:,:,1)
        jmeanGLOBAL = jmeanGLOBAL * ((2.*xmax)**2./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

        open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat', access='direct',form='unformatted',status='replace',recl=i)
        write(u,rec=1)jmeanGLOBAL(:,:,:,1)
        close(u)

        ! tmp = 0.
        ! do i = 1, nzg
        !     do j = 1, nxg
        !         do k = 1, nyg
        !             tmp(i) = tmp(i) + jmeanGLOBAL(i,j,k,1)
        !         end do
        !     end do
        ! end do

        ! tmp = tmp / (nxg*nzg)

        ! open(newunit=u,file='x-plot.dat')
        ! do i = nyg, 1, -1
        !     write(u,*)yface(nzg-i+1),tmp(i)
        ! end do

        ! open(newunit=u,file=trim(fileplace)//'im/image-809-'//trim(fn)//'.dat',status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imageGLOBAL(j,i,1)+imageGLOBAL(j,i,3),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)

        ! open(newunit=u,file=trim(fileplace)//'im/image-1064-'//trim(fn)//'.dat',status='replace')
        ! do i = -((Nbins-1)/2), ((Nbins-1)/2)
        !     write(u,*) (imageGLOBAL(j,i,2)+imageGLOBAL(j,i,4),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        ! end do
        ! close(u)

    end subroutine writer
end MODULE writer_mod
