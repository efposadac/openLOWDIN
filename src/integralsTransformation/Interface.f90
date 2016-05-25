module Interface_
    use, intrinsic ::  iso_c_binding
    implicit none

    interface
        subroutine c_test(coeff, ints, nao) bind(C, name="c_test")
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: coeff
            type(c_ptr), value :: ints
            integer(c_int), value :: nao
        end subroutine c_test
    end interface


end module Interface_