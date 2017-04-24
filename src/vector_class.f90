Module vector_class

    type :: vector
        real :: x, y, z
        contains

        procedure :: magnitude => magnitude_fn
        procedure :: print => print_sub
    end type vector

    interface operator (.dot.)
        module procedure vec_dot
    end interface

    interface operator (/)
        module procedure vec_div_scal
    end interface

    interface operator (*)
        module procedure vec_mult_vec
        module procedure vec_mult_scal
        module procedure scal_mult_vec
    end interface

    interface operator(+)
        module procedure vec_add_vec
        module procedure vec_add_scal
        module procedure scal_add_vec
    end interface

    interface operator(-)
        module procedure vec_minus_vec
    end interface

    private
    public :: magnitude, vector, print, operator(.dot.), operator(/), operator(*), operator(+), operator(-)

    contains

        type(vector) function vec_minus_vec(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vec_minus_vec = vector(a%x - b%x, a%y - b%y, a%z - b%z)

        end function vec_minus_vec


        type(vector) function vec_add_scal(a, b)

            implicit none

            type(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vec_add_scal = vector(a%x + b, a%y + b, a%z + b)

        end function vec_add_scal

        type(vector) function scal_add_vec(a, b)

            implicit none

            type(vector), intent(IN) :: b
            real,         intent(IN) :: a

            scal_add_vec = vector(b%x + a, b%y + a, b%z + a)

        end function scal_add_vec


        type(vector) function vec_add_vec(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vec_add_vec = vector(a%x + b%x, a%y + b%y, a%z + b%z)

        end function vec_add_vec


        real function vec_dot(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vec_dot = a%x * b%x + a%y * b%y + a%z * b%z

        end function vec_dot


        type(vector) function vec_mult_vec(a, b)

            implicit none

            type(vector), intent(IN) :: a, b

            vec_mult_vec = vector(a%x * b%x, a%y * b%y, a%z * b%z)

        end function vec_mult_vec


        type(vector) function vec_mult_scal(a, b)

            implicit none

            type(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vec_mult_scal = vector(a%x * b, a%y * b, a%z * b)

        end function vec_mult_scal


        type(vector) function scal_mult_vec(a, b)

            implicit none

            type(vector), intent(IN) :: b
            real,         intent(IN) :: a

            scal_mult_vec = vector(a * b%x, a * b%y, a * b%z)

        end function scal_mult_vec


        type(vector) function vec_div_scal(a, b)

            implicit none

            type(vector), intent(IN) :: a
            real,         intent(IN) :: b

            vec_div_scal = vector(a%x / b, a%y / b, a%z / b)

        end function vec_div_scal


        type(vector) function magnitude_fn(this)

            implicit none

            class(vector) :: this

            real :: tmp

            tmp = sqrt(this%x + this%y + this%z)
            magnitude_fn = this / tmp

        end function magnitude_fn

        subroutine print_sub(this)

            implicit none

            class(vector) :: this

                print*,this%x, this%y, this%z

        end subroutine
end Module vector_class