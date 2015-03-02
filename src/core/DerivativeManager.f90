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
  ! use ExternalPotential_
  ! use InterPotential_
  ! use InterPotential_Manager_
  ! use ParticleManager_
  ! use LibintInterface_
  ! use LibintInterface2_
  ! use IndexMap_
  use MolecularSystem_
  use ContractedGaussian_
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
  function DerivativeManager_getElement( thisID, i, j, k, l, nameOfSpecie, otherNameOfSpecie, nuclei, component  ) result ( output )
    implicit none
    integer :: thisID
    integer :: i
    integer :: j
    integer, optional :: k
    integer, optional :: l
    character(*), optional :: nameOfSpecie
    character(*), optional :: otherNameOfSpecie
    integer, optional :: nuclei
    integer, optional :: component
    ! real(8), allocatable :: output(:)
    real(8) :: output

    type(ContractedGaussian), allocatable :: contractions(:)
    type(Exception) :: ex
    integer :: specieID
    integer :: otherSpecieID
    integer :: m
    integer :: numCartesianOrbitalI
    integer :: numCartesianOrbitalJ
    character(30) :: nameOfSpecieSelected
    real(8) :: auxCharge


    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= nameOfSpecie
    specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    call MolecularSystem_getBasisSet(specieID, contractions)

    ! numCartesianOrbitalI = contractions(i)%numCartesianOrbital
    ! numCartesianOrbitalJ = contractions(j)%numCartesianOrbital

    ! if(allocated(output)) deallocate(output)
    ! allocate(output(numCartesianOrbitalI*numCartesianOrbitalJ))

    output = 0.0_8

    ! write(*,"(A,I,A2,I,A4,I,A2,I)") "Numero de orbitales cartesianos i y j: ", i, ": ", numCartesianOrbitalI, " -> ", j, ": ", numCartesianOrbitalJ

    ! select case (thisID)


    ! case( OVERLAP_DERIVATIVES )

    !    output = ContractedGaussian_overlapDerivative( &
    !         ParticleManager_getContractionPtr( specieID,  numberOfContraction = i ), &
    !         ParticleManager_getContractionPtr( specieID,  numberOfContraction = j ), nuclei, component )

    ! case( KINETIC_DERIVATIVES )

    !    !                output = ContractedGaussian_kineticDerivative( &
    !    !                    ParticleManager_getContractionPtr( specieID,  numberOfContraction=i ), &
    !    !                    ParticleManager_getContractionPtr( specieID,  numberOfContraction=j ), nuclei, component ) / &
    !    !                    ParticleManager_getMass( specieID )

    ! case( ATTRACTION_DERIVATIVES )

    !    !                output = 0.0_8
    !    !
    !    !                auxCharge = ParticleManager_getCharge( specieID )
    !    !
    !    !                do m=1, ParticleManager_getNumberOfPuntualParticles()
    !    !
    !    !                    output = output + ContractedGaussian_attractionDerivative( &
    !    !                            ParticleManager_getContractionPtr( specieID,  numberOfContraction=i ), &
    !    !                            ParticleManager_getContractionPtr( specieID,  numberOfContraction=j ), &
    !    !                            ParticleManager_getOriginOfPuntualParticle( m ), nuclei, component, &
    !    !                            ParticleManager_getOwnerOfPuntualParticle( m ) ) * &
    !    !                            auxCharge * ParticleManager_getChargeOfPuntualParticle( m )
    !    !                end do

    ! case( REPULSION_DERIVATIVES )

    !    !                if ( present(k) .and. present(l) ) then
    !    !
    !    !                    specieID = ParticleManager_getSpecieID( nameOfSpecie = trim( nameOfSpecie ) )
    !    !
    !    !                    otherSpecieID = ParticleManager_getSpecieID( nameOfSpecie = trim( otherNameOfSpecie ) )
    !    !
    !    !                    output = ContractedGaussian_repulsionDerivative( &
    !    !                            ParticleManager_getContractionPtr( specieID,  numberOfContraction=i), &
    !    !                            ParticleManager_getContractionPtr( specieID,  numberOfContraction=j), &
    !    !                            ParticleManager_getContractionPtr( otherSpecieID,  numberOfContraction=k), &
    !    !                            ParticleManager_getContractionPtr( otherSpecieID,  numberOfContraction=l), &
    !    !                            nuclei, component )
    !    !
    !    !                    output = output * ( ParticleManager_getCharge( specieID ) &
    !    !                                * ParticleManager_getCharge( specieID ) )

    !    !                else
    !    !
    !    !                    call Exception_constructor( ex , ERROR )
    !    !                    call Exception_setDebugDescription( ex, "Class object IntegralManager in the getElement(i,j) function" )
    !    !                    call Exception_setDescription( ex, "" )
    !    !                    call Exception_show( ex )
    !    !
    !    !                end if

    ! case default

    !    output = 0.0_8

    !    call Exception_constructor( ex , WARNING )
    !    call Exception_setDebugDescription( ex, "Class object IntegralManager in the get(i) function" )
    !    call Exception_setDescription( ex, "This ID hasn't been defined, returning zero value" )
    !    call Exception_show( ex )

    !    return

    ! end select


  end function DerivativeManager_getElement

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
