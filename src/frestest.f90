module types

    implicit none

    type :: point
        real :: x, y, z
    end type

    interface operator(.dot.)
        module procedure dot_prod
    end interface

    interface operator(-)
        module procedure vec_minus
    end interface

    interface operator(/)
        module procedure vec_scal_divide
    end interface


    contains

        type(point) function vec_scal_divide(a, b)

            implicit none

            type(point), intent(IN) :: a
            real ,       intent(IN) :: b

            vec_scal_divide = point(a%x/b, a%y/b, a%z/b)

        end function vec_scal_divide


        real function dot_prod(a, b)

            implicit none

        type(point), intent(IN) :: a, b

        dot_prod = a%x * b%x + a%y * b%y + a%z * b%z

        end function dot_prod


        type(point) function vec_minus(a, b)

            implicit none

            type(point), intent(IN) :: a, b

            vec_minus = point(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus


        real function mag(a)

            implicit none

            type(point), intent(IN) :: a

            mag = sqrt(a%x**2 + a%y**2 + a%z**2)

        end function mag


        function fresnel(pdir, n1, n2) result (tir)
        !calculates the fresnel coefficents
        !
        !
            implicit none

            real, intent(IN) :: n1, n2, pdir
            real             :: crit, costt, sintt, sint2, cost2, tir(2), f1, f2

            crit = n2/n1

            costt = abs(pdir)
            sintt = sqrt(1. - costt * costt)

            if(sintt > crit)then
                tir = 1.0
                return
            else
                sint2 = (n1/n2)*sintt
                cost2 = sqrt(1. - sint2 * sint2)
                f1 = abs((n1*costt - n2*cost2) / (n1*costt + n2*cost2))**2
                f2 = abs((n1*cost2 - n2*costt) / (n1*cost2 + n2*costt))**2

                tir = [f1,f2]!0.5 * (f1 + f2)
            ! if(isnan(tir) .or. tir > 1. .or. tir < 0.)print*,'TIR: ', tir!, f1, f2, cost,sint,cost,sint2
                return
            end if
       
        end function fresnel
end module types


program test

    use types   

    implicit none


    type(point) :: o, p, vec1, vec2
    real    :: dir, n1, n2, cost,x, y,pi, tir(2)
    integer :: i, u

    pi = 4.*atan(1.)

    o = point(0.,0.,0.)
    n1 = 1.
    n2 = 1.5

    vec1 = point(1.,0.,0.) - o
    vec1 = vec1/mag(vec1)
    dir = 0.
    x = 1.
    y = 0.

    open(newunit=u, file='fres_test.dat')

    do i = 0, 100

        p = point(x,y,0.)
        vec2 = p - o
        vec2 = vec2/mag(vec2)

        cost = vec1 .dot. vec2 
        
        tir = fresnel(cost, n2, n1)
        write(u,*) x,y,acos(cost)*180/pi,tir

        y = y - .01
        x = sqrt(.5 - y*y)
    end do
    close(u)
end program test