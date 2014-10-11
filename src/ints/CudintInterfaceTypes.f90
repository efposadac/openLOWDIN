!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Modulo que contiene los diferentes tipos a usarse en la interface Lowdin - LIBINT
!!
!! @author Edwin Posada
!!
!! <b> Creation data : </b> 10-07-10
!!
!! <b> History change: </b>
!!
!!   - <tt> 10-07-10 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# Modulo que contiene los diferentes tipos a usarse en la interface Lowdin - LIBINT
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
!! @todo Para tener total claridad sobre el tamaño de los arreglos mire los archivos libint.h, libderiv.h y libr12.h
!! @par info: Configurado con opciones de compilacion por omision de libint V.1.04
!! @version 2.0
module CudintInterfaceTypes_
  use Exception_
  use, intrinsic :: iso_c_binding
  implicit none
  
  !>
  !!@brief libint.h structure
  type, bind(c) :: lib_int
     type(c_ptr)     :: int_stack
     type(c_ptr)     :: PrimQuartet
     real(c_double)  :: AB(3)
     real(c_double)  :: CD(3)
     type(c_ptr)     :: vrr_classes(9,9)
     type(c_ptr)     :: vrr_stack
  end type lib_int
  
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
  
end module CudintInterfaceTypes_
