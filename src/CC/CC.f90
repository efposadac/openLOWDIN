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
  use CCSD_
  use String_
  implicit none



  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)

  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  !! Loads General information
  call CoupledCluster_constructor()
  call CoupledCluster_init()
  
  select case(CONTROL_instance%COUPLED_CLUSTER_LEVEL)

    case("CCSD")

      print*, "CoupledCluster_"      
      call CCSD_constructor()


    case default

      ! call Exception_.....

  end select

  
  !!stop time

  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time Coupled Cluster : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*, *) ""
  close(30)


end program CC
