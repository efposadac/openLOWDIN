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
  use Math_
  implicit none
    
  !>
  !! @brief Conversion factors from atomic units. Updated March-2025, according to CODATA values  F.M.
  !! CODATA Internationally recommended 2022 values of the Fundamental Physical Constants
  !!
  !a0=5.29177210544x10-11 m 
  real(8) , parameter :: METER = 5.29177210544E-11_8
  real(8) , parameter :: ANGSTROM = METER*1.0E10_8
  !e=1.602176634E-19 C
  real(8) , parameter :: COULOMB = 1.602176634E-19_8
  !Debye = (1/c)*1E-21 C m^2/s
  !c= 299792458 m s-1    
  real(8) , parameter :: DEBYE = METER*COULOMB*299792458.0_8/1.0E-21_8
  !me=9.1093837139x10-31 kg
  real(8) , parameter :: KG = 9.1093837139E-31_8
  !NA=6.02214076x10+23
  real(8) , parameter :: AMU = 1000_8*6.02214076E+23_8*KG
  real(8) , parameter :: DALTON = 1/AMU

  real(8) , parameter :: DEGREES = 180.0_8/Math_PI
  !Eh=2.1947463136314x10+7 m^-1
  real(8) , parameter :: CM_NEG1 = 219474.63136314_8
  !Eh=27.211386245981 eV
  real(8) , parameter :: EV=27.211386245981_8
  !Eh=4.3597447222060x10-18 J 
end module Units_
