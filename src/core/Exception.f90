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

!>
!! @brief  Modulo para manejo de errores (Importado de APMO)
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2006-03-10
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-05-15 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapto al estandar de codificacion propuesto.
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin
module Exception_
  use Stopwatch_

  type, public :: Exception
     character(100) :: debugDescription
     character(255) :: description
     integer :: typeMessage
  end type Exception

  integer, public, parameter :: INFORMATION = 1
  integer, public, parameter :: WARNING = 2
  integer, public, parameter :: ERROR = 3

  public :: &
       Exception_constructor, &
       Exception_destructor, &
       Exception_show, &
       Exception_setDebugDescription, &
       Exception_setDescription

contains

  subroutine Exception_constructor( this, typeMessage,  debugDescription, description )
    implicit none
    type(Exception) :: this
    integer :: typeMessage
    character(*), optional :: debugDescription
    character(*), optional :: description

    select case ( typeMessage )

    case( INFORMATION )
       this%typeMessage=INFORMATION

    case( WARNING )
       this%typeMessage=WARNING

    case( ERROR )
       this%typeMessage=ERROR

    end select

    if ( present(debugDescription) ) then

       this%debugDescription = trim(debugDescription)
    else

       this%debugDescription = "Undefined"

    end if

    if ( present(description) ) then

       this%description = trim(description)

    else

       this%description = "Undefined"

    end if

  end subroutine Exception_constructor

  subroutine Exception_destructor( this )
    implicit none
    type(Exception) :: this

    this%description = ""
    this%debugDescription = ""

  end subroutine Exception_destructor

  subroutine Exception_show( this )
    implicit none
    type(Exception) :: this

    select case  ( this%typeMessage )

    case( INFORMATION )

       print *,""
       print *,"+++ INFORMATION +++"
       print *,"   Debug description: ", trim(this%debugDescription)
       print *,"   Description: ", trim(this%description)
       print *,"+++"
       print *,""

    case( WARNING )

       print *,""
       print *,"!!! WARNING !!!"
       print *,"   Debug description: ", trim(this%debugDescription)
       print *,"   Description: ", trim(this%description)
       print *,"!!!"
       print *,""

    case( ERROR )

       print *,""
       print *,"### ERROR ###"
       print *,"   Debug description: ", trim(this%debugDescription)
       print *,"   Description: ", trim(this%description)
       print *,"###"
       print *,""
       call Stopwatch_stop( lowdin_stopwatch )
       write(6, *)
       write(6,"(A16,ES10.2,A4)") "Elapsed Time : ", lowdin_stopwatch%enlapsetTime ," (s)"
       write(6,*) "lowdin execution terminated ABNORMALLY at : ", trim( Stopwatch_getCurretData( lowdin_stopwatch ) )
       call Stopwatch_destructor( lowdin_stopwatch )
       stop

    end select


  end subroutine Exception_show

  subroutine Exception_setDebugDescription( this, debugDescription )
    implicit none
    type(Exception) :: this
    character(*) :: debugDescription

    this%debugDescription = trim(debugDescription)

  end subroutine Exception_setDebugDescription

  subroutine Exception_setDescription( this , description )
    implicit none
    type(Exception) :: this
    character(*) :: description

    this%description = trim(description)

  end subroutine Exception_setDescription


end module Exception_
