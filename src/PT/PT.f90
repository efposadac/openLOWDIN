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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
! this is a change
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program PT
  use CONTROL_
  use MolecularSystem_
  use InputCI_
  use Exception_
  use Vector_
  use PropagatorTheory_
  use String_
  implicit none

  character(50) :: job

  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  ! if ( .not. CONTROL_instance%LOCALIZE_ORBITALS) then
     !!Load the system in lowdin.sys format
     call MolecularSystem_loadFromFile( "LOWDIN.SYS" )
  ! else
  !    !!Load the system in lowdin.sys format
  !    call MolecularSystem_loadFromFile( "LOWDIN.SYS", "lowdin-subsystemA" )
  ! end if

  call InputCI_constructor( )
  call InputCI_load( MolecularSystem_getNumberOfQuantumSpecies() )
  
  call PropagatorTheory_constructor( CONTROL_instance%PT_ORDER )
  call PropagatorTheory_run()
  call PropagatorTheory_show()
  call PropagatorTheory_destructor()

  !!stop time
  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL CPU Time PT : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time PT : ", lowdin_stopwatch%elapsetWTime ," (s)"
  write(*, *) ""
  close(30)


end program PT
