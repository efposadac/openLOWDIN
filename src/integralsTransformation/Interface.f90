module Interface_
    use, intrinsic ::  iso_c_binding
    implicit none

    interface
        subroutine Interface_integralsTransform(coeff, ints, nao, lp, up, lq, uq, lr, ur, ls, us) bind(C, name="c_integrals_transform")
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), value :: coeff
            type(c_ptr), value :: ints
            integer(c_int), value :: nao
            integer(c_int), value :: lp, up, lq, uq, lr, ur, ls, us
        end subroutine Interface_integralsTransform
    end interface


end module Interface_