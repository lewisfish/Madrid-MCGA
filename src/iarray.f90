MODULE iarray
!
!  Contains all array var names.
!
    implicit none
    save

    real, allocatable :: xface(:), yface(:), zface(:)
    real, allocatable :: refrac(:,:,:)
    real, allocatable :: rhokap(:,:,:,:), albedo_a(:,:,:,:)
    real, allocatable :: image(:,:,:), imageGLOBAL(:,:,:)
    real, allocatable :: jmean(:,:,:,:), jmeanGLOBAL(:,:,:,:)
    real, allocatable :: absorb(:,:,:), absorbGLOBAL(:,:,:)

end MODULE iarray
