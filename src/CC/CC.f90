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
!! @brief Coupled Cluster (CC) program
!!        This module allows to make CC calculations in the APMO framework
!! @author  A. Peña (apenat@unal.edu.co)
!!
!! <b> Creation date : </b> 2015
!!
!! <b> History: </b>
!!
!!
program CC
  use CONTROL_
  use MolecularSystem_
  use Exception_
  use CoupledCluster_
  use CCD_
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

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

 
 print *," Control: ", CONTROL_instance%COUPLED_CLUSTER_LEVEL


  select case ( trim(CONTROL_instance%COUPLED_CLUSTER_LEVEL) )

    case ( "CCSD" )


       call CoupledCluster_constructor(CONTROL_instance%COUPLED_CLUSTER_LEVEL )
       call CoupledCluster_run()
       call CoupledCluster_show()
       call CoupledCluster_destructor()

       !!stop time
       call Stopwatch_stop(lowdin_stopwatch)
  
       write(*, *) ""
       write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time CC : ", lowdin_stopwatch%enlapsetTime ," (s)"
       write(*, *) ""

    case ( "CCD" )


 print *," Control: ", CONTROL_instance%COUPLED_CLUSTER_LEVEL
       call CCD_constructor(CONTROL_instance%COUPLED_CLUSTER_LEVEL )
       call CCD_run()
       call CCD_show()
       call CCD_destructor()

       !!stop time
       call Stopwatch_stop(lowdin_stopwatch)
  
       write(*, *) ""
       write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time CC : ", lowdin_stopwatch%enlapsetTime ," (s)"
       write(*, *) ""

 print *," Control: ", CONTROL_instance%COUPLED_CLUSTER_LEVEL

    case default

      call CoupledCluster_exception( ERROR, "Coupled interactor constructor", "Correction level not implemented")

  end select

  close(30)


end program CC
