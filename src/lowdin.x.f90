!!******************************************************************************
!!      This code is part of LOWDIN Quantum chemistry package
!!
!!      this program has been developed under direction of:
!!
!!      Prof. A REYES' Lab. Universidad Nacional de Colombia
!!              http://www.qcc.unal.edu.co
!!      Prof. R. FLORES' Lab. Universidad de Guadalajara
!!              http://www.cucei.udg.mx/~robertof
!!
!!	(c) All rights reserved, 2013
!!
!!******************************************************************************

!>
!! @brief Source of the program, based on the apmo.x.f90 file.
!! @author E. F. Posada (efposadac@unal.edu.co)
!! <b> Creation data : </b> 02-11-2011
!! <b> History change: </b>
!!   - <tt> 02-11-2011 </tt>:  E. F. Posada ( efposadac@unal.edu.co )
!!        -# Creation and design of program based on APMO code
!!   - <tt> 03-24-2013 </tt>:  E. F. Posada ( efposadac@unal.edu.co )
!!        -# Rewritten to match with new lowdin core standard
!! @info: this program compiles by linking with lowdincore.a library
program lowdin_
 use CONTROL_
 use Stopwatch_
 use InputManager_
 use MolecularSystem_
 use MecanicProperties_
 use GeometryOptimizer_
 use Solver_
 implicit none

 character(50) :: strAuxNumber
 integer :: statusSystem

 !! Time Control
 call Stopwatch_constructor( lowdin_stopwatch )
 call Stopwatch_start( lowdin_stopwatch )
 
 !! Show credits
 write (6,*) "LOWDIN execution started at : ", trim( Stopwatch_getCurretData( lowdin_stopwatch ) )
 write (6,*) "---------------------------------------------------------------"
 write (6,*) ""
 write (6,*) "***************************************************************"
 write (6,*) "*                   LOWDIN 2.0  (May/2013)                    *"
 write (6,*) "*                                                             *"
 write (6,*) "*  R. FLORES-MORENO, E. F. POSADA, F. S. MONCADA, J.ROMERO,   *"
 write (6,*) "*  J. CHARRY, M. DIAZ-TINOCO, S. A. GONZALEZ, N. F. AGUIRRE,  *"
 write (6,*) "*  A. REYES                                                   *"
 write (6,*) "*                                                             *"
 write (6,*) "*  https://sites.google.com/site/lowdinproject/               *" 
 write (6,*) "*                                                             *"
 write (6,*) "***************************************************************"
 write (6,*) ""
 
 !!***************************************************************************
 !! Load input and build the molecular system
 !!
 write(6, "(1A)", advance="no") " PARSING INPUT..."


 !! Load info for system being calculated
 call InputManager_loadSystem() 

 !! Load CONTROL block 
 call InputManager_loadControl() 
 
 !! Load TASKS block
 call InputManager_loadTask()   

 !! Load GEOMETRY block
 call InputManager_loadGeometry() 

 !! Load potentials if any
 call InputManager_loadPotentials()

 !! Load OUTPUTY block
! call InputManager_loadOUTPUT() 
 
 write(6, "(1A)") " DONE!"

 !! Shows running parameters
 call CONTROL_show()   
 
 !! Builds the molecular system 
 call MolecularSystem_build()
 
 !! Shows some information related to molecular system
 call MolecularSystem_showInformation()  

 if (CONTROL_instance%METHOD/="MM") then 
 
    call MolecularSystem_showParticlesInformation()

    !!
    !!****************************************************************************

    !!***************************************************************************
    !!        Shows system's geometry
    !!
    write (6,"(T20,A30)") " INITIAL GEOMETRY: AMSTRONG"
    write (6,"(T18,A35)") "------------------------------------------"
 
    call MolecularSystem_showCartesianMatrix()
 
 end if 

 !! Transform to center of mass
 call MecanicProperties_constructor(MolecularSystem_instance%mechanicalProp)
  
 if (CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS .and. (.not. CONTROL_instance%ARE_THERE_DUMMY_ATOMS)) then
    
    call MolecularSystem_moveToCenterOfMass()
    call MolecularSystem_rotateOnPrincipalAxes()
    write (6,"(T20,A30)") " GEOMETRY IN C.M. : AMSTRONG"
    write (6,"(T18,A35)") "------------------------------------------"
    call MolecularSystem_showCartesianMatrix()
    
 end if
 
 if (CONTROL_instance%METHOD/="MM") then 
    call MolecularSystem_showDistanceMatrix()
 end if
  
 !! At this moment it is not relevant
 !! call MolecularSystem_showZMatrix( MolecularSystem_instance )  
 !!
 !!****************************************************************************
  
 !!***************************************************************************
 !!        Save checkpoint (lowdin.sys, lowdin.dat and lowdin.bas)
 !!
 call MolecularSystem_saveToFile()
 !! DEBUG
 !call MolecularSystem_loadFromFile("LOWDIN.DAT")
 !call CONTROL_show()
 !call MolecularSystem_loadFromFile("LOWDIN.SYS")
 !call MolecularSystem_showInformation()  
 !call MolecularSystem_showParticlesInformation()
 !call MolecularSystem_loadFromFile("LOWDIN.BAS")
 !call MolecularSystem_showInformation()  
 !call MolecularSystem_showParticlesInformation()
 !!
 !!****************************************************************************
   
  
 !!***************************************************************************
 !!        Running the properly solver for a selected method
 !!
 ! call Solver_run()
  
 ! if  (Parameters%OPTIMIZE .and. .not. Parameters%ELECTRONIC_WAVEFUNCTION_ANALYSIS &
 !      .and. .not.Parameters%ARE_THERE_DUMMY_ATOMS .and. .not. Parameters%CP_CORRECTION ) then
 

 !    call GeometryOptimizer_constructor( lowdin_geometryOptimizer, ssolver = lowdin_solver)
 !    call GeometryOptimizer_run( lowdin_geometryOptimizer )
 !    call GeometryOptimizer_destructor( lowdin_geometryOptimizer )
 
 ! else if ( Parameters%CP_CORRECTION ) then
 
 !    !          call ExternalSoftware_constructor( external_instance )
 !    !          call ExternalSoftware_makeBSSEwMolecularSystem_moveToCenterOfMassithCP( external_instance )
 
 ! else if  (Parameters%TDHF .and. .not. Parameters%ELECTRONIC_WAVEFUNCTION_ANALYSIS &
 !      .and. .not.Parameters%ARE_THERE_DUMMY_ATOMS .and. .not. Parameters%CP_CORRECTION ) then
 
 !    call TimeEvolution_constructor( lowdin_TimeEvolution, ssolver = lowdin_solver)
 !    call TimeEvolution_run( lowdin_TimeEvolution )
 !    call TimeEvolution_destructor( lowdin_TimeEvolution )

 if ( CONTROL_instance%OPTIMIZE ) then 
    call GeometryOptimizer_constructor( GeometryOptimizer_instance )
    call GeometryOptimizer_run( GeometryOptimizer_instance )
    call GeometryOptimizer_destructor( GeometryOptimizer_instance )
 else
    call Solver_run()
 end if

 !!
 !!******************************************************************************

  statusSystem = system ("lowdin-CalcProp.x")

  if ( CONTROL_instance%IS_THERE_OUTPUT ) then
    write(strAuxNumber,"(I10)") Input_instance%numberOfOutputs
!!  call system("lowdin-output.x" //trim(strAuxNumber))
    statusSystem = system ("lowdin-output.x" //trim(strAuxNumber))
  end if


 !!Cleaning
 call MolecularSystem_destroy()

 !!Shows time information 
 call Stopwatch_stop(lowdin_stopwatch)
 write(*, *)
 write(*,"(A,F10.3,A4)") "  TOTAL CPU Time: ", lowdin_stopwatch%enlapsetTime ," (s)"
 write(*,"(A,F10.3,A4)") " TOTAL Wall Time: ", lowdin_stopwatch%elapsetWTime ," (s)"
 write(6,"(A16,i3,A1,i3,A1,i3,A1,i4,A2)") &
      "Elapsed Time: ", &
      lowdin_stopwatch%endTime(5),"h", &
      lowdin_stopwatch%endTime(6),"m", &
      lowdin_stopwatch%endTime(7),"s", &
      lowdin_stopwatch%endTime(8),"ms"
 
 write (6,"(A, A)") "LOWDIN execution terminated normally at : ", trim( Stopwatch_getCurretData( lowdin_stopwatch ) )
 call Stopwatch_destructor( lowdin_stopwatch )

end program lowdin_
