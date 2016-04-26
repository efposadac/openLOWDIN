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
!! @brief Derivative Manager Module.
!!        This module contains all basic functions of the derivatives calculations
!! @author  E. Posada
!! @author  J.M. Rodas
!! 
!! <b> Creation date : </b> 2011-12-14
!!
!! <b> History: </b>
!!
!!   - <tt> 2011-12-14 </tt>: Edwin Fernando Posada ( efposadac@unal.edu.co )
!!        -# Basics functions and functions has been created
!!   - <tt> 2015-02-27 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Rewrite the code to Lowdin v 2.0 and prepare the module for new methods
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module DerivativeManager_
#ifdef intel
  use IFPORT
#endif
  use MolecularSystem_
  use ContractedGaussian_
  use KineticDerivatives_
  use AttractionDerivatives_
  use OverlapDerivatives_
  use RepulsionDerivatives_
  use CosmoCore_
  use Exception_
  use Math_
  use Matrix_
  use String_
  use Exception_
  implicit none

  integer, parameter, public :: OVERLAP_DERIVATIVES       = 7
  integer, parameter, public :: KINETIC_DERIVATIVES       = 8
  integer, parameter, public :: ATTRACTION_DERIVATIVES    = 9
  integer, parameter, public :: MOMENT_DERIVATIVES        = 10
  integer, parameter, public :: MOMENTUM_DERIVATIVES      = 11
  integer, parameter, public :: REPULSION_DERIVATIVES     = 12
  integer, parameter, public :: TWOPARTICLE_REPULSION_DERIVATIVES = 13


  type, public :: DerivativeManager
     character(20) :: name
     logical :: isInstanced
  end type DerivativeManager

  public :: &
       DerivativeManager_constructor, &
       DerivativeManager_destructor, &
                                ! DerivativeManager_show, &
       DerivativeManager_getElement!, &
  ! DerivativeManager_getLabels

  private
contains

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine DerivativeManager_constructor(this)
    implicit none
    type(DerivativeManager) :: this

    this%isInstanced = .true.

  end subroutine DerivativeManager_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine DerivativeManager_destructor(this)
    implicit none
    type(DerivativeManager) :: this

    this%isInstanced = .false.

  end subroutine DerivativeManager_destructor


  !**
  ! @brief Retorna valor de las derivadas sin almacenar su valor en ningun arreglo
  !
  ! @param thisID Corresponde al valor de la instancia. Las instancias soportadas son: OVERLAP_INTEGRALS,
  !               KINETIC_INTEGRALS, ATTRACTION_INTEGRALS, MOMENT_INTEGRALS, MOMENTUM_INTEGRALS,
  !               OVERLAP_DERIVATIVES, KINETIC_DERIVATIVES, ATTRACTION_DERIVATIVES
  !**
  subroutine DerivativeManager_getElement( thisID, deriveVector, surface, i, j, k, l, nameOfSpecie, otherNameOfSpecie, A, B )
    implicit none
    integer :: thisID
    real(8), allocatable :: deriveVector(:)
    integer :: i
    integer :: j
    integer, optional :: k
    integer, optional :: l
    character(*), optional :: nameOfSpecie
    character(*), optional :: otherNameOfSpecie
    integer, optional :: A
    integer, optional :: B
    type(ContractedGaussian), allocatable :: contractions(:)
    type(ContractedGaussian), allocatable :: otherContractions(:)
    type(Exception) :: ex
    integer :: specieID
    integer :: otherSpecieID
    character(30) :: nameOfSpecieSelected
    type(surfaceSegment), intent(in), optional :: surface

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) ) then
       nameOfSpecieSelected= nameOfSpecie
    end if
    specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
    call MolecularSystem_getBasisSet(specieID, contractions)

    if ( present( otherNameOfSpecie ) ) then
       otherSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=otherNameOfSpecie )
    end if

    call MolecularSystem_getBasisSet(specieID, contractions)


    select case (thisID)

    case( KINETIC_DERIVATIVES )

       call KineticDerivatives_getDerive( contractions, i, j,  deriveVector, specieID)

    case( ATTRACTION_DERIVATIVES )

       if(present(surface)) then
          call AttractionDerivatives_getDerive( contractions, i, j,  deriveVector, A, B, specieID, surface)
       else
          call AttractionDerivatives_getDerive( contractions, i, j,  deriveVector, A, B, specieID)
       end if

    case( OVERLAP_DERIVATIVES )

       call OverlapDerivatives_getDerive( contractions, i, j,  deriveVector, specieID)

    case( REPULSION_DERIVATIVES )

       if ( present(k) .and. present(l) ) then
          call RepulsionDerivatives_getDerive( contractions, i, j, k, l, deriveVector, specieID)
       else
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object DerivativeManager in the getElement(i,j) function" )
          call Exception_setDescription( ex, "For repulsion derivatives is necessary four indices (p,q|r,s)" )
          call Exception_show( ex )
       end if

    case( TWOPARTICLE_REPULSION_DERIVATIVES )
       call RepulsionDerivatives_getInterDerive(i, j, k, l, deriveVector, specieID, otherSpecieID)

    case default

       call Exception_constructor( ex , WARNING )
       call Exception_setDebugDescription( ex, "Class object DerivativeManager in the get(i) function" )
       call Exception_setDescription( ex, "This ID hasn't been defined, returning zero value" )
       call Exception_show( ex )

       return

    end select

  end subroutine DerivativeManager_getElement

  ! !<
  ! !! @brief Devuelve los indices necesarios para guardar integrales
  ! !!>
  ! function DerivativeManager_getLabels(specieID, numberOfContractions) result(labelsOfContractions)
  !   implicit none
  !   integer:: specieID
  !   integer:: numberOfContractions
  !   integer:: labelsOfContractions(numberOfContractions)

  !   integer:: auxLabelsOfContractions
  !   integer:: i

  !   auxLabelsOfContractions = 1

  !   do i = 1, numberOfContractions

  !      !!Posicion real de las contracciones con respecto al numero total de contracciones
  !      labelsOfContractions(i) = auxLabelsOfContractions

  !      auxLabelsOfContractions = auxLabelsOfContractions + ParticleManager_instance%particlesPtr(ParticleManager_instance%idsOfContractionsForSpecie(&
  !           specieID)%contractionID(i)%particleID)%basis%contractions( &
  !           ParticleManager_instance%idsOfContractionsForSpecie(specieID)%contractionID(&
  !           i)%contractionIDInParticle)%numCartesianOrbital
  !   end do

  ! end function DerivativeManager_getLabels


  ! !>
  ! !! @brief Muestra informacion del objeto
  ! !!
  ! !! @param this
  ! !<
  ! subroutine DerivativeManager_show(this)
  !   implicit none
  !   type(DerivativeManager) :: this
  ! end subroutine DerivativeManager_show

  ! !!>
  ! !! @brief Indica si el objeto ha sido instanciado o no
  ! !!
  ! !<
  ! function DerivativeManager_isInstanced( this ) result( output )
  !   implicit  none
  !   type(DerivativeManager), intent(in) :: this
  !   logical :: output

  !   output = this%isInstanced

  ! end function DerivativeManager_isInstanced

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine DerivativeManager_exception( typeMessage, description, debugDescription)
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

  end subroutine DerivativeManager_exception

end module DerivativeManager_
