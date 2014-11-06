
!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************



	!>
	!!
	!!  Este modulo define una seudoclase para calculo de propiedades derivadas de
	!! la funcion de onda como cargas, dipolos, polarizabilidades, etc.
	!!
	!! @author Sergio A. Gonzalez Monico
	!!
	!! <b> Fecha de creacion : </b> 2007-09-18
	!!
	!! <b> Historial de modificaciones: </b>
	!!
	!!   - <tt> 2007-09-18 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
	!!        -# Creacion de modulo y metodos basicos.
	!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
	!!   - <tt> 2011-11-23 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
	!!        -# Adds numerical integration properties, ADPT calculations and brings population analyses 
	!!   - <tt> 2014-01-23 </tt>: Matheus Rodriguez ( matrodriguezalv@unal.edu.co )
        !!        -# Reescribe y adapta el modulo de Calculate properties en Lowdin2
	!<

program CalcProp_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use Units_
  use WaveFunction_
  use ContractedGaussian_
  use CalculateProperties_ ! module name
  implicit none

  character(50) :: job

  job = ""
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))


  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

!  call CalculatePropertiesmod_constructor () !modificar
 ! call CalculatePropertiesmod_destructor ()

 !! Calculate properties subroutines
  call CalculateProperties_showPopulationAnalyses()
    
end program CalcProp_





