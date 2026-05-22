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
module Solver_
  use CONTROL_
  use SinglePoint_
  use GeometryOptimizer_
  implicit none

  public Solver_main

contains

  subroutine Solver_main()
    implicit none

    if (CONTROL_instance%OPTIMIZE) then
      call GeometryOptimizer_constructor(GeometryOptimizer_instance)
      call GeometryOptimizer_run(GeometryOptimizer_instance)
      call GeometryOptimizer_destructor(GeometryOptimizer_instance)
    else
      call SinglePoint_run()
    end if

  end subroutine Solver_main


end module Solver_
