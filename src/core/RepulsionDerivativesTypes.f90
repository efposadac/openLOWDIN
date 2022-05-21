!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package
!!
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Libint 1.04 interface for Repulsion Derivative Module.
!!        This module contains all basic types of the repulsion derivatives calculations
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2015-03-17
!!
!! <b> History: </b>
!!
!!   - <tt> 2015-03-17 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs,
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module RepulsionDerivativesTypes_
  use Exception_
  use, intrinsic :: iso_c_binding
  implicit none

  INTEGER, PARAMETER :: libint_max_am = 5
  INTEGER, PARAMETER :: libderiv_max_am1 = 4
  INTEGER, PARAMETER :: prim_data_f_size = 4 * (libint_max_am - 1) + 1 
  INTEGER, PARAMETER :: libint_vrr_classes_size = 2 * (libint_max_am - 1) + 1  
  INTEGER, PARAMETER :: libint_dvrr_classes_size = 2 * (libderiv_max_am1 - 1) + 1 
  INTEGER, PARAMETER :: build_eri_size = libint_max_am - 1
  INTEGER, PARAMETER :: build_deriv1_eri_size = libderiv_max_am1 - 1 

  !> @brief libint.h structure
  type, bind(c) :: prim_data
     real(c_double) :: F(17)
     real(c_double) :: U(3,6)
     real(c_double) :: twozeta_a
     real(c_double) :: twozeta_b
     real(c_double) :: twozeta_c
     real(c_double) :: twozeta_d
     real(c_double) :: oo2z
     real(c_double) :: oo2n
     real(c_double) :: oo2zn
     real(c_double) :: poz
     real(c_double) :: pon
     real(c_double) :: oo2p
     real(c_double) :: ss_r12_ss
  end type prim_data

  !>
  !!@brief libint.h structure
  type, bind(c) :: lib_int
     type(c_ptr)  :: int_stack
     type(c_ptr)     :: PrimQuartet
     real(c_double)  :: AB(3)
     real(c_double)  :: CD(3)
     type(c_ptr)  :: vrr_classes(9,9)
     type(c_ptr)  :: vrr_stack
  end type lib_int

  !>
  !!@brief libderiv.h structure
  type, bind(c) :: lib_deriv
     type(c_ptr) :: int_stack
     type(c_ptr)    :: PrimQuartet
     type(c_ptr) :: zero_stack
     type(c_ptr) :: ABCD(12+144)
     real(c_double) :: AB(3)
     real(c_double) :: CD(3)
     type(c_ptr) :: deriv_classes(12,7,7)
     type(c_ptr) :: deriv2_classes(144,7,7)
     type(c_ptr) :: dvrr_classes(7,7)
     type(c_ptr) :: dvrr_stack
  end type lib_deriv
  
end module RepulsionDerivativesTypes_
