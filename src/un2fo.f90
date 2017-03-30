program converter

    implicit none
    
    integer :: i, u,j 
    double precision :: array(200,200,200,4), tmp(200,200,200)
    character(len=60) :: fileplace,fn


    fileplace = "/home/lewis/phdshizz/Madrid paper/data/jmean/"


    ! do j = 1, 5
        write(fn,'(i1)')j
        inquire(iolength=i)array
        open(newunit=u, file=trim(fileplace)//'jmean- 0.000.dat', access='direct', form='unformatted', recl=i)
        print*,trim(fileplace)//'jmean--0.'//trim(fn)//'00.dat'
        read(u,rec=1)array

        close(u)

        tmp = 0.
        tmp(:,:,:) = array(:,:,:,1) + array(:,:,:,3)

        inquire(iolength=i)tmp(:,:,:)
        open(newunit=u, file=trim(fileplace)//'809-large.dat', access='direct', form='unformatted', recl=i)
        write(u,rec=1)tmp(:,:,:)
        close(u)
        tmp = 0.
        tmp(:,:,:) = array(:,:,:,2) + array(:,:,:,4)

        inquire(iolength=i)tmp(:,:,:)
        open(newunit=u, file=trim(fileplace)//'1064-large.dat', access='direct', form='unformatted', recl=i)
        write(u,rec=1)tmp(:,:,:)
        close(u)
    ! end do
end program converter