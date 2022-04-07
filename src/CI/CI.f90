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
  use NonOrthogonalCI_
  use ConfigurationInteraction_
  use String_
  use InputCI_
  implicit none

  character(50) :: job
  integer :: numberOfSpeciesInCI

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

  if( trim(job) .eq. "NOCI") then

     call NonOrthogonalCI_constructor(NonOrthogonalCI_instance)
     call NonOrthogonalCI_displaceGeometries(NonOrthogonalCI_instance)
     call NonOrthogonalCI_buildOverlapAndHamiltonianMatrix(NonOrthogonalCI_instance)
     call NonOrthogonalCI_diagonalizeCImatrix(NonOrthogonalCI_instance)
     call NonOrthogonalCI_plotDensities(NonOrthogonalCI_instance)
  else
     
     read(job,"(I10)") numberOfSpeciesInCI

     call InputCI_constructor( )
     if(numberOfSpeciesInCI .ne. 0) then
        call InputCI_load( numberOfSpeciesInCI )
     else
        call InputCI_load( MolecularSystem_getNumberOfQuantumSpecies() )
     end if
     call ConfigurationInteraction_constructor(CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL )
     call ConfigurationInteraction_run()
     call ConfigurationInteraction_show()
     call ConfigurationInteraction_showEigenVectors()
     call ConfigurationInteraction_densityMatrices()
     call ConfigurationInteraction_destructor()

  end if
  
  !!stop time
  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL CPU Time CI : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time CI : ", lowdin_stopwatch%elapsetWTime ," (s)"
  write(*, *) ""


end program CI
