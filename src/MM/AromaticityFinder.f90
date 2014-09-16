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
!!        This module evaluates the aromaticity in a ring
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
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

  !>
  !! @brief This function evaluates in an atom belongs to an aromatic ring
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class with informations of the rings
  !! @param [in] atomIdx INTEGER atom to evaluate
  !! @return [out] output LOGICAL if the ring is aromatic returns .true.
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

