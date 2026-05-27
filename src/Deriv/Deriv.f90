!!******************************************************************************
!!        This code is part of LOWDIN Quantum chemistry package
!!
!!        this program has been developed under direction of:
!!
!!        Prof. A REYES' Lab. Universidad Nacional de Colombia
!!                http://www.qcc.unal.edu.co
!!        Prof. R. FLORES' Lab. Universidad de Guadalajara
!!                http://www.cucei.udg.mx/~robertof
!!
!!                Todos los derechos reservados, 2013
!!
!!******************************************************************************
!> @brief This program calculates nuclear derivatives of one-particle integrals for contracted gaussian functions representations

module Deriv_
  use CONTROL_
  use MolecularSystem_
  use EnergyGradients_
  implicit none

  public Deriv_main

contains

  subroutine Deriv_main(auxjob)
    implicit none
    character(len=*) :: auxjob
    character(50) :: job

    job = trim(String_getUppercase(auxjob))
  
    !!Start time
    call Stopwatch_constructor(lowdin_stopwatch)
    call Stopwatch_start(lowdin_stopwatch)
  
    !!Load CONTROL Parameters
    call MolecularSystem_loadFromFile("LOWDIN.DAT")
  
    !!Load the system in lowdin.sys format
    call MolecularSystem_loadFromFile("LOWDIN.SYS")
  
    select case (trim(job))

    case ("GET_GRADIENTS")
      call EnergyGradients_constructor()
      call EnergyGradients_getAnalyticDerivative()

    case default
  
      write (*, *) "USAGE: lowdin-ints.x job "
      write (*, *) "Where job can be: "
      write (*, *) "  GET_GRADIENTS"
      stop "ERROR"
  
    end select

  end subroutine Deriv_main

end module Deriv_
