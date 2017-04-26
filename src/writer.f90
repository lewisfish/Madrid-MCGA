MODULE writer_mod

implicit none
save

CONTAINS
    subroutine writer(depth)

        use constants, only : fileplace, nbins, nzg, nyg, nxg
        use iarray,    only : imageGLOBAL, jmeanGLOBAL,zface

        implicit none

        integer           :: i, u, j, k
        real              :: depth, jmeanf(nzg)
        character(len=10) :: fn


        inquire(iolength=i)jmeanGLOBAL

        open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat',access='direct',form='unformatted',status='replace',recl=i)
        write(u,rec=1)jmeanGLOBAL
        close(u)

        jmeanf = 0.

        do i = 1, nzg
            do j = 1, nxg
                do k = 1 , nyg
                    jmeanf(i) = jmeanf(i) + jmeanGLOBAL(j,k,i,1)
                end do
            end do
        end do

        jmeanf = jmeanf / real(nyg*nxg)

        open(newunit=u,file='z-dep.dat',status='replace')

        do i = 1, nzg
            write(u,*) zface(nzg-i), jmeanf(i)
        end do

        close(u)

        write(fn,'(F4.2)') (2.*depth-.3)

        open(newunit=u,file=trim(fileplace)//'im/image-809-'//trim(fn)//'.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) (imageGLOBAL(j,i,1)+imageGLOBAL(j,i,3),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)

        open(newunit=u,file=trim(fileplace)//'im/image-1064-'//trim(fn)//'.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) (imageGLOBAL(j,i,2)+imageGLOBAL(j,i,4),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)

    end subroutine writer
end MODULE writer_mod
