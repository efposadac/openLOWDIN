!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!		R. Flores (roberto.floresmoreno.qt@gmail.com)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module InputOutput_
  use Exception_
  use OutputBuilder_
  use Vector_
  implicit none

  !>
  !! @brief Description: Reads information from the input, as a namelist, for additional output files 
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 08-04-11
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 08-04-11 </tt>:  felix ( email@server )
  !!        -# description.
  !!   - <tt> 10-31-2014 </tt>:  Jorge Charry ( jacharry@unal.edu.co )
  !!        -# Adapts this module to Lowdin2
  !!   - <tt> 21-11-2024 </tt>:  Felix ( email@server )
  !!        -# Simplifies routines, adds more options
  !!
  !<
  character(50) :: Output_type
  character(50) :: Output_species
  character(2) :: Output_plane
  character(1) :: Output_axis
  integer :: Output_state
  integer :: Output_orbital
  integer :: Output_dimensions
  integer :: Output_pointsPerDim
  real(8) :: Output_scanStep
  real(8) :: Output_cubeSize
  real(8) :: Output_minValue
  real(8) :: Output_maxValue
  real(8) :: Output_offsetX
  real(8) :: Output_offsetY
  real(8) :: Output_offsetZ
  real(8) :: Output_limitX(2)
  real(8) :: Output_limitY(2)
  real(8) :: Output_limitZ(2)
  real(8) :: Output_center(3)
  real(8) :: Output_point1(3)
  real(8) :: Output_point2(3)
  real(8) :: Output_point3(3)
  
  NAMELIST /Output/ &
       Output_type, &
       Output_species, &
       Output_plane, &
       Output_axis, &
       Output_state, &
       Output_orbital, &
       Output_dimensions, &
       Output_pointsPerDim, &
       Output_scanStep, &
       Output_cubeSize, &
       Output_minValue, &
       Output_maxValue, &
       Output_offsetX, &
       Output_offsetY, &
       Output_offsetZ, &
       Output_limitX, &
       Output_limitY, &
       Output_limitZ, &
       Output_center, &
       Output_point1, &
       Output_point2, &
       Output_point3

  public :: &
       InputOutput_load

  private	

contains

  !>
  !! @brief Carga la informacion de potenciales externos desde el input
  !!
  !! @param this
  !<
  subroutine  InputOutput_load( outputObjects )
    implicit none
    type(OutputBuilder) :: outputObjects(:)
    integer :: i
    integer :: stat
    character(1000) :: line

    open (unit=4, file=trim(CONTROL_instance%INPUT_FILE)//"aux")
    rewind(4)
    do i=1, size(outputObjects)
       Output_type=""
       Output_species="ALL"
       Output_plane=""
       Output_axis=""
       Output_state=1
       Output_orbital=0
       Output_dimensions=0
       Output_pointsPerDim=0
       Output_scanStep=0.0_8
       Output_cubeSize=0.0_8
       Output_minValue=0.0_8
       Output_maxValue=0.0_8
       Output_offsetX=0.0_8
       Output_offsetY=0.0_8
       Output_offsetZ=0.0_8
       Output_limitX(:)=0.0_8
       Output_limitY(:)=0.0_8
       Output_limitZ(:)=0.0_8
       Output_center(:)=0.0_8
       Output_point1(:)=0.0_8
       Output_point2(:)=0.0_8
       Output_point3(:)=0.0_8
       read(4,NML=Output, iostat=stat)

       if( stat > 0 ) then       
          write (*,'(A)') 'Error reading Output block'
          backspace(4)
          read(4,fmt='(A)') line
          write(*,'(A)') 'Invalid line : '//trim(line)
          call Exception_stopError("Class object InputOutput in the load function", &
               "check the OUTPUTS block in your input file")
       end if

       call OutputBuilder_constructor( outputs_instance(i), i, &
            Output_type, &
            Output_species, &
            Output_plane, &
            Output_axis, &
            Output_state, &
            Output_orbital, &
            Output_dimensions, &
            Output_pointsPerDim, &
            Output_scanStep, &
            Output_cubeSize, &
            Output_minValue, &
            Output_maxValue, &
            Output_offsetX, &
            Output_offsetY, &
            Output_offsetZ, &
            Output_limitX, &
            Output_limitY, &
            Output_limitZ, &
            Output_center, &
            Output_point1, &
            Output_point2, &
            Output_point3)

    end do
    
    close(4)
    
  end subroutine InputOutput_load

end module InputOutput_
