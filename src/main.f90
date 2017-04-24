program main

    use mpi
    use MCf
    use subs
    use iarray, only : image

    implicit none
    
    integer :: error, numproc, id, i
    real    :: depth
        

    call MPI_init(error)

    call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)


    call directory
    call alloc_array

    depth = 0.5 !depth = zmax
    do i = 1, 1
        call mcpolar(depth, id, numproc)
        depth = depth - 0.1
    end do

    call mpi_finalize(error)
end program main