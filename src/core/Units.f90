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
!! @brief Modulo para definicion y ajuste de constantes globales
!!
!! Este modulo contiene la definicion de constantes globales utilizadas dentro de
!! modulos y sus procedimiento asociados, al igual que algunos metodos que ajustan
!! sus valores de forma adecuada en tiempo de ejecucion.
!!
!! @author Sergio Gonz\'alez
!!
!! <b> Creation data : </b> 02-12-11
!!
!! <b> History change: </b>
!!
!!   - <tt> 2007-05-15 </tt>: Sergio A. Gonz\'alez M. ( sagonzalezm@unal.edu.co )
!!        -# Creaci\'on de m\'odulo y procedimientos b\'asicos
module Units_
  implicit none
    
  !>
  !! @brief Definici\'on de factores de conversi\'on
  !!
  real(8) , parameter :: HARTREE = 1.0_8
  real(8) , parameter :: AMSTRONG = 0.52917724924_8
  real(8) , parameter :: DEBYES = 2.541764_8
  real(8) , parameter :: ELECTRON_REST = 1.0_8
  real(8) , parameter :: AMU = 1.0_8/1822.88850065855_8
  real(8) , parameter :: DALTON = 1822.88850065855_8
  real(8) , parameter :: kg = 9.109382616D-31
  real(8) , parameter :: DEGREES = 57.29577951_8
  real(8) , parameter :: CM_NEG1 = 219476.0_8
    
end module Units_
