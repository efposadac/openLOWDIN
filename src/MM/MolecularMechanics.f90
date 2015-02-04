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
!! @brief Molecular Mechanics program.
!!        This module allows to make calculations using Force Fields (UFF)
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions using Universal Force Field has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program MolecularMechanics
  use CONTROL_
  use MolecularSystem_
  use String_
  use MMFunctions_
  use Exception_
  implicit none

       character(50) :: job
       character(50) :: ffmethod
       logical :: electrostaticEnergy
       logical :: printAllMM


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


  ffmethod = CONTROL_instance%FORCE_FIELD
  electrostaticEnergy = CONTROL_instance%ELECTROSTATIC_MM
  printAllMM = CONTROL_instance%PRINT_MM
  !! Initializes the object MolecularMechanics
  call MolecularMechanics_constructor()
  !! Initializes calculation using the force field selected
  call MolecularMechanics_run(ffmethod,electrostaticEnergy,printAllMM)
  ! call MolecularMechanics_show()
  call MolecularMechanics_destructor()

  !!stop time
  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time Molecular Mechanics : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*, *) ""
  close(30)

end program MolecularMechanics
