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
!! @brief Hessians Module.
!!        This module contains all basic functions of the hessians calculations
!! @author  J.M. Rodas
!! @author  S.A. Gonzalez
!!
!! <b> Creation date : </b> 2009-06-25
!!
!! <b> History: </b>
!!
!!   - <tt> 2009-06-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Basics functions has been created using empirical method
!!   - <tt> 2015-02-24 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0 and prepare the module for new methods
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Hessians_
  use CONTROL_
  use Matrix_
  use Vector_
  use AtomicElement_
  use ParticleManager_
  use MolecularSystem_
  use Exception_

  implicit none

  type, public :: Hessians
     character(20) :: name
     logical :: isInstanced
  end type Hessians

  type(Hessians), target :: Hessians_instance

  public :: &
       Hessians_constructor, &
       Hessians_destructor, &
       Hessians_getEmpirical

  private

contains

  !>
  !! @brief Define class constructor
  !<
  subroutine Hessians_constructor( this )
    implicit none
    type(Hessians) :: this

    this%isInstanced =.true.

  end subroutine Hessians_constructor

  !>
  !! @brief Define class destructor
  !<
  subroutine Hessians_destructor( this)
    implicit none
    type(Hessians) :: this

    this%isInstanced =.false.

  end subroutine Hessians_destructor

  !>
  !! @brief   Retorna una aproximacion a la matriz hessiana, basada en las
  !!      formulas empiricas de Fischer y Almlof para calculo de constantes
  !!      de fuerza. J Phys Chem. 96, 24 1992, 9768-9774
  !!      de las variables independientes
  !<
  function Hessians_getEmpirical( this , system )  result( output )
    implicit none
    type(Hessians) :: this
    type(MolecularSystem) :: system
    type(Matrix) :: output
    type(Vector) ::  cartesianCoordinates
    type(AtomicElement) :: element
    type(Exception) :: ex
    character(10) :: auxSymbol
    real(8) :: covalentRadius
    real(8) :: strengthConstant
    real(8) :: particlesDistance
    real(8) :: componentsOfStengthConstant(3,3)
    real(8) :: origin(3)
    integer :: numberOfCoordinates
    integer :: iteratorOfParticles_1
    integer :: iteratorOfParticles_2
    integer :: beginOfParticle
    integer :: beginOfOtherParticle
    integer :: iteratorOfCenter
    integer :: iteratorOfOtherCenter
    integer :: i
    integer :: j

    i = ParticleManager_getNumberOfCentersOfOptimization()
    i = i*3

    call Matrix_constructor(output, int(i,8), int(i,8) )
    output%values = 0.0_8
    iteratorOfCenter = 0
    iteratorOfOtherCenter = 0

    ! call AtomicElement_constructor( element )

    ! !! Calculo de elementos fuera de la diagonal
    do iteratorOfParticles_1=1, size(ParticleManager_instance)
       if (   ParticleManager_isCenterOfOptimization( iteratorOfParticles_1 ) ) then
          iteratorOfCenter = iteratorOfCenter + 1
          iteratorOfOtherCenter =0

          do iteratorOfParticles_2=1, iteratorOfParticles_1-1
             if (   ParticleManager_isCenterOfOptimization( iteratorOfParticles_2 )  )  then
                iteratorOfOtherCenter = iteratorOfOtherCenter + 1
                auxSymbol=trim( ParticleManager_getSymbol( iteratorOfParticles_1 ) )

                if (  scan( auxSymbol, "_" ) /= 0 ) then
                   auxSymbol = trim(auxSymbol(1: scan( auxSymbol, "_" ) - 1 ) )
                end if

                call AtomicElement_load( element, auxSymbol, 0 )
                covalentRadius = element%covalentRadius
    !             !                       print *,"radio covalente para ",trim(auxSymbol),element%covalentRadius
                auxSymbol=trim(ParticleManager_getSymbol( iteratorOfParticles_2 ))

                if (  scan( auxSymbol, "_" ) /= 0 ) then
                   auxSymbol = trim(auxSymbol(1: scan( auxSymbol, "_" ) - 1 ) )
                end if

                call AtomicElement_load( element, auxSymbol, 0 )
                covalentRadius = covalentRadius + element%covalentRadius
    !             !                       print *,"radio covalente para ",trim(auxSymbol),element%covalentRadius

                origin = ParticleManager_getOrigin( iterator = iteratorOfParticles_1 )
                origin = origin - ParticleManager_getOrigin( iterator = iteratorOfParticles_2 )
                origin = origin * ANGSTROM
                particlesDistance = sum( origin*origin )
                strengthConstant = 0.3601_8 * exp( -1.944_8*( dsqrt(particlesDistance) -  covalentRadius ) )

                componentsOfStengthConstant(1,1) = (origin(1)*origin(1)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(1,2) = (origin(1)*origin(2)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(1,3) = (origin(1)*origin(3)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(2,2) = (origin(2)*origin(2)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(2,3) = (origin(2)*origin(3)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(3,3) = (origin(3)*origin(3)* strengthConstant)/particlesDistance
                componentsOfStengthConstant(2,1) =componentsOfStengthConstant(1,2)
                componentsOfStengthConstant(3,1) =componentsOfStengthConstant(1,3)
                componentsOfStengthConstant(3,2) =componentsOfStengthConstant(2,3)

                beginOfParticle = (iteratorOfCenter - 1) * 3
                beginOfOtherParticle = (iteratorOfOtherCenter - 1) * 3
                do i=1,3
                   do j=1,3

                      output%values( beginOfParticle+i, beginOfParticle +j) = &
                           output%values( beginOfParticle+i, beginOfParticle +j) &
                           + componentsOfStengthConstant(i,j)

                      output%values( beginOfOtherParticle+i, beginOfOtherParticle +j) = &
                           output%values( beginOfOtherParticle+i, beginOfOtherParticle +j) &
                           + componentsOfStengthConstant(i,j)

                      output%values( beginOfParticle+i, beginOfOtherParticle +j) = &
                           output%values( beginOfParticle+i, beginOfOtherParticle +j) &
                           - componentsOfStengthConstant(i,j)

                      output%values( beginOfOtherParticle+i, beginOfParticle +j) = &
                           output%values( beginOfOtherParticle+i, beginOfParticle +j) &
                           - componentsOfStengthConstant(i,j)
                   end do
                end do

             end if
          end do
       end if
    end do

    do i=1, size(output%values, dim=1)
       if(output%values(i,i) < 0.00001_8 ) output%values(i,i) = 0.00001_8
    end do

    ! call AtomicElement_destructor(element)

  end function Hessians_getEmpirical

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine Hessians_exception( typeMessage, description, debugDescription )
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine Hessians_exception


end module Hessians_
