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
!!        (c) All rights reserved, 2013
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
  use Solver_,          only : Solver_main
  use Output_,          only : Output_main
  implicit none

  character(50) :: strAuxNumber
  type(Stopwatch) :: global_stopwatch

  !! Time Control
  call Stopwatch_constructor(global_stopwatch)
  call Stopwatch_start(global_stopwatch)

  write (6, *) "LOWDIN execution started at : ", trim(Stopwatch_getCurretData(global_stopwatch))
  write (6, *) "---------------------------------------------------------------------------------"

  !! Show logo
  write (6, *) "                                                                                 "
  write (6, *) "                  0000        00000      00000000000           000               " 
  write (6, *) "                 0000000         000000  000                   000               " 
  write (6, *) "                      000           000                 000    000               " 
  write (6, *) "           00          000          000  000000000       0000  000               " 
  write (6, *) "           000         000  000 00000    000          000  0000000               " 
  write (6, *) "            000       000   000                       000    00000               " 
  write (6, *) "             00000000000    000          00000000000  000      000               " 
  write (6, *) "                                                                                 " 
  write (6, *) "                                                                                 " 
  write (6, *) "  00                0000     00             00  0000000       000           000  " 
  write (6, *) "  00               0000000   000            00  00    00000   000           000  " 
  write (6, *) "  00                    000   000   00   0000   00       000         000    000  " 
  write (6, *) "  00         00          000   000  000  000              00  000     0000  000  " 
  write (6, *) "  00         000         000    000  000000              000  000  000   000000  " 
  write (6, *) "  00          0000     0000      00   0000              000   000  000     0000  " 
  write (6, *) "  000000000     000000000         0    00       0000000000    000  000       00  " 
  write (6, *) "                                                                                 "
  write (6, *) "---------------------------------------------------------------------------------"

  !! Show credits
  write (6, *) ""
  write (6, *) "---------------------------------------------------------------------------------"
  write (6, *) "                          OPEN LOWDIN   (May/2026)                               "
  write (6, *) "                                                                                 "
  write (6, *) "            R. FLORES-MORENO, E. F. POSADA, F. S. MONCADA, J.ROMERO,             "
  write (6, *) "            J. CHARRY, M. DIAZ-TINOCO, S. A. GONZALEZ, N. F. AGUIRRE,            "
  write (6, *) "            A. REYES                                                             "
  write (6, *) "                                                                                 "
  write (6, *) "            https://github.com/efposadac/openLOWDIN/                             "
  write (6, *) "                                                                                 "
  write (6, *) "---------------------------------------------------------------------------------"
  write (6, *) ""

  !!***************************************************************************
  !! Load input build the molecular system
  !!
  write (6, "(1A)", advance="no") " PARSING INPUT..."

  !! Load info for system being calculated
  call InputManager_loadSystem()

  !! Load CONTROL block
  call InputManager_loadControl()

  !! Load TASKS block
  call InputManager_loadTask()

  !! Load GEOMETRY block and initialize Molecular system
  call InputManager_loadGeometry()

  !! Load potentials if any
  call InputManager_loadPotentials()

  write (6, "(1A)") " DONE!"

  !! Shows running parameters
  call CONTROL_show()

  !!
  !!****************************************************************************

  !!***************************************************************************
  !! Builds the molecular system
  !!
  call MolecularSystem_build()

  !! Shows some information related to molecular system
  call MolecularSystem_showInformation()

  if (CONTROL_instance%METHOD /= "MM") then

    call MolecularSystem_showParticlesInformation()

    !! Shows system's geometry
    write (6, "(T20,A30)") " INITIAL GEOMETRY: ANGSTROM"
    write (6, "(T18,A35)") "------------------------------------------"

    call MolecularSystem_showCartesianMatrix(molecularSystem_instance)

  end if


  !! Transform to center of mass
  if (CONTROL_instance%TRANSFORM_TO_CENTER_OF_MASS .and. (.not. CONTROL_instance%ARE_THERE_DUMMY_ATOMS)) then

    call MolecularSystem_moveToCenterOfMass()
    call MolecularSystem_rotateOnPrincipalAxes()
    write (6, "(T20,A30)") " GEOMETRY IN C.M. : ANGSTROM"
    write (6, "(T18,A35)") "------------------------------------------"
    call MolecularSystem_showCartesianMatrix(molecularSystem_instance)

  end if

  if (CONTROL_instance%METHOD /= "MM") then
    call MolecularSystem_showDistanceMatrix()
  end if

  !! call MolecularSystem_showZMatrix( MolecularSystem_instance )

  !!
  !!****************************************************************************

  !!***************************************************************************
  !!        Save checkpoint (lowdin.sys, lowdin.dat and lowdin.bas)
  !!
  call MolecularSystem_saveToFile()
  !!
  !!****************************************************************************

  !!***************************************************************************
  !!        Running the properly solver for a selected method
  !!
  call Solver_main()
  !!
  !!******************************************************************************

  if (CONTROL_instance%IS_THERE_OUTPUT) then
    write (strAuxNumber, "(I10)") Input_instance%numberOfOutputs
    call Output_main(trim(strAuxNumber))
  end if

  !!Cleaning
  call MolecularSystem_destroy()

  !!Shows time information
  call Stopwatch_stop(global_stopwatch)
  write (*, *)
  write (*, "(A,F10.3,A4)") "  TOTAL CPU Time: ", global_stopwatch%enlapsetTime, " (s)"
  write (*, "(A,F10.3,A4)") " TOTAL Wall Time: ", global_stopwatch%elapsetWTime, " (s)"
  write (6, "(A16,i3,A1,i3,A1,i3,A1,i4,A2)") &
    "Elapsed Time: ", &
    global_stopwatch%endTime(5), "h", &
    global_stopwatch%endTime(6), "m", &
    global_stopwatch%endTime(7), "s", &
    global_stopwatch%endTime(8), "ms"

  write (6, "(A, A)") "LOWDIN execution terminated normally at : ", trim(Stopwatch_getCurretData(global_stopwatch))
  call Stopwatch_destructor(global_stopwatch)

end program lowdin_
