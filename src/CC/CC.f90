!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!    http://www.qcc.unal.edu.co/
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Coupled Cluster program
!!        This program allows to make calculations in Coupled Cluster over APMO approach (CC-APMO)
!!                from energy and amplitude equations obtained from SympyMaquinola (SymPyMaq see more: caraortizmah@una.edu.co) (based on SymPy). 
!! @author  Carlos Andres Ortiz Mahecha (CAOM) (caraortizmah@unal.edu.co)
!!
!! <b> Creation date : </b> 2016-10-26
!!
!! <b> History: </b>
!!
!!   - <tt> 2016-10-26 </tt>: (CAOM) ( caraortizmah@unal.edu.co )
!!        -# Development of Coupled Cluster (CC) program:
!!                This Program calls a number of modules that solve equations of CC-APMO.
!!   - <tt> data </tt>:  
!!
!!
!! @warning <em>  All characters and events in this Quantum chemistry package -- even those based on real source code -- are entirely fictional. </br>
!!                All celebrity lines are impersonated.....poorly. </br> 
!!                The following program contains corase language and due to it's cintent should not be viewed by anyone. </em>
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
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

      call CCSD_constructor()
      call CCSD_run()


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
