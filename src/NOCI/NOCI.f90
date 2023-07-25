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
!! @brief Non-Orthogonal Configuration Interaction (CI) program
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
program NOCI
  use CONTROL_
  use InputManager_
  use MolecularSystem_
  use Exception_
  use NonOrthogonalCI_
  use String_
  use InputCI_
  use Stopwatch_
  use MecanicProperties_
  implicit none

  character(50) :: job
  character(50) :: strAuxNumber

  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  if( trim(job) .eq. "POSTSCF") then !POST SCF after a main lowdin2 calculation

     !!Start time
     call Stopwatch_constructor(lowdin_stopwatch)
     call Stopwatch_start(lowdin_stopwatch)

     write (6,*) "***************************************************************"
     write (6,*) "*                   LOWDIN-NOCI  (Dec/2022)                   *"
     write (6,*) "*                                                             *"
     write (6,*) "*               F. MONCADA, L.G.M. PETTERSON                  *"
     write (6,*) "*                                                             *"
     write (6,*) "***************************************************************"

     !!Load CONTROL Parameters
     call MolecularSystem_loadFromFile( "LOWDIN.DAT" )
     !!Load the system in lowdin.sys format
     call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

     call NonOrthogonalCI_constructor(NonOrthogonalCI_instance)
     if(CONTROL_instance%READ_NOCI_GEOMETRIES) then
        call NonOrthogonalCI_readGeometries(NonOrthogonalCI_instance)
     else
        call NonOrthogonalCI_displaceGeometries(NonOrthogonalCI_instance)
     end if
     call NonOrthogonalCI_runHFs(NonOrthogonalCI_instance)
     call NonOrthogonalCI_buildOverlapAndHamiltonianMatrix(NonOrthogonalCI_instance)
     call NonOrthogonalCI_diagonalizeCImatrix(NonOrthogonalCI_instance)
     call NonOrthogonalCI_generateDensities(NonOrthogonalCI_instance)

     !!stop time
     call Stopwatch_stop(lowdin_stopwatch)
    
     write(*, *) ""
     write(*,"(A,F10.3,A4)") "** TOTAL CPU Time NOCI : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time NOCI : ", lowdin_stopwatch%elapsetWTime ," (s)"
     write(*, *) ""

  else !STAND ALONE calculation

     !!Start time
     call Stopwatch_constructor(lowdin_stopwatch)
     call Stopwatch_start(lowdin_stopwatch)

     !! Show credits
     write (6,*) "LOWDIN-NOCI execution started at : ", trim( Stopwatch_getCurretData( lowdin_stopwatch ) )
     write (6,*) "---------------------------------------------------------------"
     write (6,*) ""
     write (6,*) "***************************************************************"
     write (6,*) "*                   LOWDIN-NOCI  (Dec/2022)                   *"
     write (6,*) "*                                                             *"
     write (6,*) "*               F. MONCADA, L.G.M. PETTERSON                  *"
     write (6,*) "*                                                             *"
     write (6,*) "*                       based on                              *"
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
     call MolecularSystem_showParticlesInformation()

     !!***************************************************************************
     !!        Shows system's geometry
     !!
     write (6,"(T20,A30)") " INITIAL GEOMETRY: AMSTRONG"
     write (6,"(T18,A35)") "------------------------------------------"
     call MolecularSystem_showCartesianMatrix(molecularSystem_instance)

     !! Transform to center of mass
     call MecanicProperties_constructor(MolecularSystem_instance%mechanicalProp)

     if (CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS .and. (.not. CONTROL_instance%ARE_THERE_DUMMY_ATOMS)) then

        call MolecularSystem_moveToCenterOfMass()
        call MolecularSystem_rotateOnPrincipalAxes()
        write (6,"(T20,A30)") " GEOMETRY IN C.M. : AMSTRONG"
        write (6,"(T18,A35)") "------------------------------------------"
        call MolecularSystem_showCartesianMatrix(molecularSystem_instance)

     end if

     if (CONTROL_instance%METHOD/="MM") then 
        call MolecularSystem_showDistanceMatrix()
     end if

     !!***************************************************************************
     !!        Save checkpoint (lowdin.sys, lowdin.dat and lowdin.bas)
     !!
     call MolecularSystem_saveToFile()
     !!
     !!
     !Do SCF
     select case ( trim(CONTROL_instance%METHOD) )

     case('RHF')
        call system("lowdin-SCF.x RHF")
     case('UHF')
        call system("lowdin-SCF.x UHF")
     case('RKS')
        call system("lowdin-SCF.x RKS")
     case('UKS')
        call system("lowdin-SCF.x UKS")
     case default
        STOP "The method: "//trim(CONTROL_instance%METHOD)//" is not implemented"
     end select

     ! if ( trim(CONTROL_instance%INTEGRAL_STORAGE) .ne. "DIRECT" ) then
     !    CONTROL_instance%INTEGRAL_STORAGE="DIRECT"
     !    write(*,*)"This program only works for DIRECT integral storage!"
     ! end if

     if ( CONTROL_instance%NONORTHOGONAL_CONFIGURATION_INTERACTION ) then
        call NonOrthogonalCI_constructor(NonOrthogonalCI_instance)
        call NonOrthogonalCI_displaceGeometries(NonOrthogonalCI_instance)
        call NonOrthogonalCI_runHFs(NonOrthogonalCI_instance)
        call NonOrthogonalCI_buildOverlapAndHamiltonianMatrix(NonOrthogonalCI_instance)
        call NonOrthogonalCI_diagonalizeCImatrix(NonOrthogonalCI_instance)
        call NonOrthogonalCI_generateDensities(NonOrthogonalCI_instance)
     end if

     call MolecularSystem_saveToFile()
     
     !!calculate CI density properties
     call system ("lowdin-CalcProp.x")

     if ( CONTROL_instance%IS_THERE_OUTPUT ) then
        write(strAuxNumber,"(I10)") Input_instance%numberOfOutputs
        call system("lowdin-output.x" //trim(strAuxNumber))
        ! statusSystem = system ("lowdin-output.x" //trim(strAuxNumber))
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

     write (6,"(A, A)") "LOWDIN-NOCI execution terminated normally at : ", trim( Stopwatch_getCurretData( lowdin_stopwatch ) )
     call Stopwatch_destructor( lowdin_stopwatch )
  end if
end program NOCI
