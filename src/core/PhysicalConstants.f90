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
!! @brief  Modulo para definicion y ajuste de constantes globales
!!
!! Este modulo contiene la definicion de constantes globales utilizadas dentro de
!! modulos y sus procedimiento asociados, al igual que algunos metodos que ajustan
!! sus valores de forma adecuada en tiempo de ejecucion.
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2006-03-10
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-05-15 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapta al estandar de codificacion propuesto.
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapt the module to LOWDIN package.
module PhysicalConstants_
  implicit none
    
  !<
  !! Define algunas constantes numericas
  !! definicion de algunas constantes atomicas y nucleares en unidades atomicas.
  real(8) , parameter :: PhysicalConstants_ELECTRON_CHARGE = -1.0_8
  real(8) , parameter :: PhysicalConstants_ELECTRON_MASS = 1.0_8
  real(8) , parameter :: PhysicalConstants_PROTON_CHARGE = 1.0_8
  real(8) , parameter :: PhysicalConstants_PROTON_MASS = 1*1836.15267247_8
  real(8) , parameter :: PhysicalConstants_NEUTRON_CHARGE = 0.0_8
  real(8) , parameter :: PhysicalConstants_NEUTRON_MASS = 1*1838.6836605_8
  
  real(8) , parameter :: PhysicalConstants_SPIN_UP_ELECTRON = 0.5_8
  real(8) , parameter :: PhysicalConstants_SPIN_DOWN_ELECTRON = -0.5_8
  real(8) , parameter :: PhysicalConstants_SPIN_ELECTRON = 0.5_8
  real(8) , parameter :: PhysicalConstants_SPIN_FERMION = 0.5_8 !< solo se especifica este valor pero puede ser n*0.5 n=1,2,3,..., en U.A.
  real(8) , parameter :: PhysicalConstants_SPIN_BOSON = 1.0_8 !< solo se especifica este valor pero puede ser n=0,1,2,3,..., en U.A.
  !!***************************************************************************
  
end module PhysicalConstants_
