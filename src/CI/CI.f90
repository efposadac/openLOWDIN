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
!! @brief Configuration Interaction (CI) program
!!        This module allows to make CI calculations in the APMO framework
!! @author  F. Moncada (fsmoncadaa@unal.edu.co)
!!
!! <b> Creation date : </b> 2012-07-24
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: J.Charry  ( jacharrym@unal.edu.co )
!!        -# Adapts this module to its inclusion in Lowdin2
!!
!! @warning This programs only works linked to lowdincore library, provided by LOWDIN quantum chemistry package
!!
program CI
  use CONTROL_
  use MolecularSystem_
  use Exception_
  use ConfigurationInteraction_
  use String_
  use InputCI_
  implicit none

  character(50) :: job
  integer :: numberOfSpeciesInCI

  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))
  read(job,"(I10)"), numberOfSpeciesInCI

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  call InputCI_constructor( )
  call InputCI_load( numberOfSpeciesInCI )

  call ConfigurationInteraction_constructor(CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL )
  call ConfigurationInteraction_run()
  call ConfigurationInteraction_show()
  call ConfigurationInteraction_destructor()

  !!stop time
  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time CI : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*, *) ""
  close(30)


end program CI
