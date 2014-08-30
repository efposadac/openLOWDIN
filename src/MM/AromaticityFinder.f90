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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module AromaticityFinder_
  use MolecularSystem_
  use MatrixInteger_
  use Rings_
  use Vector_
  use MMCommons_
  use RingFinder_
  use Exception_
  implicit none


  public :: &
       AromaticityFinder_isAromatic
  
contains

  function AromaticityFinder_isAromatic( this, atomIdx ) result(output)
    implicit none
    type(Rings) :: this
    integer, intent(in) :: atomIdx
    logical :: output
    integer :: numberOfRings 
    integer :: numberOfColumns
    integer :: i, j

    output = .false.
    
    do i=1, this%numberOfRings
       do j=1, this%ringSize(i)
          if(this%connectionMatrix(i)%values(1,j)==atomIdx) then
             if(this%aromaticity(i) == 1) then
                output = .true.
             end if
          end if
       end do
    end do

  end function AromaticityFinder_isAromatic

end module AromaticityFinder_

