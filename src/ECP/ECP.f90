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
!! @brief ECP and RHF-ECP program.
!!        This module allows to make calculations using effective core potentials
!! @author  I. Ortiz-Verano
!!
!! <b> Creation date : </b> 2017-02-16
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-09-17 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos (RHF.f90)
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adpata el modulo para incluirlo en LOWDIN (RHF.90)
!!   - <tt> 2013-05-09 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Rewrite the module as a program, handles closed and open shell systems
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
Program ECP
  use CONTROL_
  use EffectiveCorePotentials_
  use Exception_

end program ECP
