MODULE writer_mod

implicit none
save

CONTAINS
    subroutine writer(depth)

        use constants, only : fileplace, nbins
        use iarray,    only : imageGLOBAL

        implicit none

        integer           :: i, u, j
        real              :: depth
        character(len=10) :: fn


        write(fn,'(F4.2)') depth


        open(newunit=u,file=trim(fileplace)//'im/im-1064-'//trim(fn)//'.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) imageGLOBAL(0,i,2)+imageGLOBAL(0,i,4)!,j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)

        open(newunit=u,file=trim(fileplace)//'im/im-809-'//trim(fn)//'.dat',status='replace')
        do i = -((Nbins-1)/2), ((Nbins-1)/2)
            write(u,*) (imageGLOBAL(j,i,1)+imageGLOBAL(j,i,3),j = ((Nbins-1)/2), -((Nbins-1)/2),-1)
        end do
        close(u)

    end subroutine writer
end MODULE writer_mod
