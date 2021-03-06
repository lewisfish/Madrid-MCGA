MODULE ch_opt

implicit none

private
public :: init_opt1, init_opt2, init_opt3, init_opt4, init_intra

CONTAINS
   
   subroutine init_opt1()
   !
   !  subroutine to set 809nm TISSUE optical properties
   !
      use opt_prop, only : hgg, g2, mua, mus, kappa, albedo
      
      implicit none

      hgg = 0.9
      g2  = hgg**2
      mua = 0.34 !from mua.dat
      mus = 13.453210015633305 / (1. - hgg) !from 96 jacques data + jacques formula

      kappa  = mus + mua 
      albedo = mus / kappa

   end subroutine init_opt1
   

   subroutine init_opt2()
   !
   !  subroutine to set 1064nm TISSUE optical properties
   !
      use opt_prop, only : hgg, g2, mua, mus, kappa, albedo

      implicit none

      hgg = 0.9
      g2  = hgg**2
      mua = 5.7 !steve jacques omlc webpage  !https://en.wikipedia.org/wiki/Near-infrared_window_in_biological_tissue
      mus = 9.488560164729575 / (1. - hgg) !from 96 jacques data + jacques formula

      kappa  = mus + mua
      albedo = mus / kappa

   end subroutine init_opt2
   

   subroutine init_opt3()
   !
   !  subroutine to set 809nm CRYSTAL optical properties
   !
      use opt_prop, only : hgg, g2, mua, mus, kappa, albedo

      implicit none

      hgg = 0.9
      g2  = hgg**2
      mua = 10.0 !from northrop grumman
      mus =  0.001 !from Scattering in Polycrystalline Nd:YAG Laser Ikesue et al. 1998

      kappa  = mus + mua
      albedo = mus / kappa

   end subroutine init_opt3
   

   subroutine init_opt4()
   !
   !  subroutine to set 1064nm CRYSTAL optical properties
   !
      use opt_prop, only : hgg, g2, mua, mus, kappa, albedo
      
      implicit none

      hgg = 0.7
      g2  = hgg**2
      mua = 0.001 !from http://users.unimi.it/aqm/wp-content/uploads/YAGBrochure.pdf 
      ! Scattering effect and laser performance for the Nd:YAG transparent ceramics  Li et al 2011
      mus = 0.004 !Scattering effect and laser performance for the Nd:YAG transparent ceramics  Li et al 2011

      kappa  = mus + mua 
      albedo = mus / kappa

   end subroutine init_opt4


   subroutine init_intra(wavecur)
   !
   !  subroutine to set intralipid optical properties
   !  from supercontinuum laser based optical characterization of intralipid phantoms in the 500-2250nm range

   use opt_prop
   
   implicit none

   real, intent(IN) :: wavecur

   real :: fact, phipr

   phipr = (2.d0 / 100.d0) * 22.7d0
   fact = phipr / .227d0
   mus = 1.868d10 * wavecur**(-2.59d0)
   mus = mus * fact

   hgg = 0.7d0
   g2  = hgg**2.
   mua = 0.d0

   kappa  = mus + mua 
   albedo = mus / kappa

   end subroutine init_intra


   subroutine sample(array, size_of, cdf, wave, iseed)
   !      
   !  samples a random value from an array based upon its cdf     
   !      
      implicit none
      
      integer, intent(IN)    :: iseed, size_of
      real,    intent(IN)    :: array(size_of, 2), cdf(size_of)
      real,    intent(OUT)   :: wave

      real :: ran2, value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(size(cdf), cdf, nlow, value)
      call lin_inter_1D(array, cdf, value, size(cdf), nlow, wave)
   
   end subroutine sample
   

   subroutine lin_inter_1D(array, cdf, value, length, nlow, y)
   !
   !  linear interpolates between values for an array and its cdf
   !   
      implicit none
   
      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2),cdf(length-1)
      integer, intent(IN)   :: nlow
   
      y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   

   subroutine lin_inter_2D(array,value,length,nlow,y)
   !
   !  linear interpolation for an array
   !
      implicit none

      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2)
      integer, intent(IN)   :: nlow
   
      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))
   
   end subroutine lin_inter_2D
   

   subroutine search_1D(length,array,nlow,value)
   !
   !  search by bisection for 1D array
   !
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_1D
   
   subroutine search_2D(length,array,nlow,value)
   !
   !  search by bisection for 2D array
   !
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length,2),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_2D
   
   subroutine mk_cdf(array,cdf,length)
   !
   !  subroutine that creates cdf for an array of values.
   !
      implicit none

      integer, intent(IN)    :: length
      real,    intent(IN)    :: array(length,2)
      real,    intent(INOUT) :: cdf(length)
      real                   :: summ
      integer                :: i,j
   
      do j = 1, length-1
         summ = 0.
         do i = 1, j   
            summ = summ + 0.5*(array(i+1, 2) + array(i, 2)) * (array(i+1, 1) - array(i, 1))
         end do
         cdf(j) = summ      
      end do
      cdf = cdf/cdf(length-1)
   
   end subroutine mk_cdf
end module ch_opt
