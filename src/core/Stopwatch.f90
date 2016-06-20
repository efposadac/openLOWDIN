!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!> @author Sergio Gonzalez
!!
!! <b> Creation data : </b> 02-11-11
!!
!! <b> History change: </b>
!!
!!   - <tt> 02-11-11 </tt>:  edwin posada (efposadac@unal.edu.co)
!!        -# Implementation of split time function
module Stopwatch_
  implicit none

  type, public :: Stopwatch
     character(20) :: name
     logical :: isInstanced
     real(8) :: startTime
     real(8) :: enlapsetTime
     real(8) :: currentTime
     character(255) :: initialDate
     character(255) :: currentDate
     integer :: initTime(8)
     integer :: endTime(8)
  end type Stopwatch

  type(Stopwatch), public :: lowdin_stopwatch

  public :: &
       Stopwatch_constructor, &
       Stopwatch_destructor, &
       Stopwatch_show, &
       Stopwatch_start, &
       Stopwatch_stop, &
       Stopwatch_splitTime, &
       Stopwatch_getCurretData

  private

contains

  
  !>
  !! @brief Default constructor 
  !! @param this
  subroutine Stopwatch_constructor( this )
    implicit none
    type(Stopwatch) :: this

    this%name="UNTITLED"
    this%startTime = 0.0_8
    call FDATE(this%initialDate)
    this%currentDate = trim(this%initialDate)
    this%enlapsetTime = 0.0_8

  end subroutine Stopwatch_constructor


  !>
  !! @brief Default destructor 
  !! @param this
  subroutine Stopwatch_destructor( this )
    implicit none
    type(Stopwatch) :: this

    this%name=""
    this%startTime = 0.0_8
    this%initialDate = ""
    this%currentDate = ""
    this%enlapsetTime = 0.0_8

  end subroutine Stopwatch_destructor



  subroutine Stopwatch_start( this )
    implicit none
    type(Stopwatch) :: this

    call cpu_time(this%startTime)
    call DATE_AND_TIME(values=this%initTime)

    this%currentTime = this%startTime

  end subroutine Stopwatch_start

  !> @brief stops the global clock
  !! @author E. F. Posada, S. A. Gonzalez
  !! @version 2.0
  subroutine Stopwatch_stop( this )
    implicit none
    type(Stopwatch) :: this
    real(8) :: currentTime
    real(8) :: iinit, eend
    
    call DATE_AND_TIME(values=this%EndTime)
    
    iinit = this%initTime(8)
    eend =  this%endTime(8)
    
    iinit = iinit/1000
    eend = eend/1000
    
    iinit = iinit + (this%initTime(3)*24*3600)+(this%initTime(5)*3600)+(this%initTime(6)*60)+this%initTime(7)
    eend = eend + (this%endTime(3)*24*3600)+(this%endTime(5)*3600)+(this%endTime(6)*60)+this%endTime(7)
    
    eend = eend - iinit
    eend = eend / 3600
    this%endTime(5) = floor(eend) 
    eend = eend - this%endTime(5)
     
    eend = eend * 60
    this%endTime(6) = floor(eend)
    eend = eend - this%endTime(6)
     
    eend = eend * 60
    this%endTime(7) = floor(eend)
    eend = eend - this%endTime(7)
     
    eend = eend * 1000
    this%endTime(8) = ceiling(eend)
    
    call cpu_time(currentTime)

    this%enlapsetTime =  currentTime - this%startTime

  end subroutine Stopwatch_stop


  subroutine Stopwatch_splitTime()
    implicit none
    real(8) :: currentTime

    call cpu_time(currentTime)

    lowdin_stopwatch%enlapsetTime =  currentTime - lowdin_stopwatch%currentTime
    lowdin_stopwatch%currentTime = currentTime

    write(*,"(F5.2)", advance="no") lowdin_stopwatch%enlapsetTime


  end subroutine Stopwatch_splitTime


  function Stopwatch_getCurretData( this ) result( output )
    implicit none
    type(Stopwatch) :: this
    character(100) :: output

    call FDATE(this%currentDate)
    output=trim(this%currentDate)

  end function Stopwatch_getCurretData

  !>
  !! @brief Muestra informacion del objeto
  !! @param this 
  subroutine Stopwatch_show(this)
    implicit none
    type(Stopwatch) :: this
  end subroutine Stopwatch_show

  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  function Stopwatch_isInstanced( this ) result( output )
    implicit  none
    type(Stopwatch), intent(in) :: this
    logical :: output
    output = this%isInstanced

  end function Stopwatch_isInstanced

end module Stopwatch_
