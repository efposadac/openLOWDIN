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
!! @brief Cuda Interface for Computation of Electron Repulsion Integrals
!!        This module allows to make calculations of ERIs on GPUs
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-10-08
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module CudintInterface_
  use Exception_
  use CudintInterfaceTypes_
  use MolecularSystem_
  use OverlapIntegrals_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use, intrinsic :: iso_c_binding
  implicit none
  
 
  interface

     subroutine cuda_int_intraspecies (&
          numberOfContractions, &
          maxLength, &
          maxNumCartesianOrbital, &
          primNormalizationSize, &
          contractionId, &
          contractionLength, &
          contractionAngularMoment, &
          contractionNumCartesianOrbital, &
          contractionOwner, &
          contractionOrigin, &
          contractionOrbitalExponents, &
          contractionCoefficients, &
          contractionContNormalization, &
          contractionPrimNormalization) bind (C, name = "cuda_int_intraspecies_")
       use, intrinsic :: iso_c_binding
       implicit none
       integer (c_int) :: numberOfContractions
       integer (c_int) :: maxLength
       integer (c_int) :: maxNumCartesianOrbital
       integer (c_int) :: primNormalizationSize
       integer (c_int) :: contractionId(*)
       integer (c_int) :: contractionLength(*)
       integer (c_int) :: contractionAngularMoment(*)
       integer (c_int) :: contractionNumCartesianOrbital(*)
       integer (c_int) :: contractionOwner(*)
       real (c_double) :: contractionOrigin(numberOfContractions,*)
       real (c_double) :: contractionOrbitalExponents(numberOfContractions,*)
       real (c_double) :: contractionCoefficients(numberOfContractions,*)
       real (c_double) :: contractionContNormalization(numberOfContractions,*)
       real (c_double) :: contractionPrimNormalization(*)
     end subroutine cuda_int_intraspecies
    
 
  end interface
  
contains
  

  subroutine CudintInterface_computeIntraSpecies(specieID)
    implicit none
    
    integer, intent(in) :: specieID
    
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: auxCounter

    integer,target :: i, j, k
    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie

    integer, allocatable :: contractionId(:)
    integer, allocatable :: contractionLength(:)
    integer, allocatable :: contractionAngularMoment(:)
    integer, allocatable :: contractionNumCartesianOrbital(:)
    integer, allocatable :: contractionOwner(:)
    real(8), allocatable :: contractionOrigin(:,:)
    real(8), allocatable :: contractionOrbitalExponents(:,:)
    real(8), allocatable :: contractionCoefficients(:,:)
    real(8), allocatable :: contractionContNormalization(:,:)
    real(8), allocatable :: contractionPrimNormalization(:)
    integer :: maxLength
    integer :: maxNumCartesianOrbital
    integer :: primNormalizationSize

  
   !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)

    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

    if(allocated(contractionId)) deallocate(contractionId)
    if(allocated(contractionLength)) deallocate(contractionLength)
    if(allocated(contractionAngularMoment)) deallocate(contractionAngularMoment)
    if(allocated(contractionNumCartesianOrbital)) deallocate(contractionNumCartesianOrbital)
    if(allocated(contractionOwner)) deallocate(contractionOwner)
    if(allocated(contractionOrigin)) deallocate(contractionOrigin)

    allocate(contractionId(numberOfContractions))
    allocate(contractionLength(numberOfContractions))
    allocate(contractionAngularMoment(numberOfContractions))
    allocate(contractionNumCartesianOrbital(numberOfContractions))
    allocate(contractionOwner(numberOfContractions))
    allocate(contractionOrigin(numberOfContractions, 3))

    primNormalizationSize=0
    do i=1, numberOfContractions
       contractionId(i) = contractions(i)%id
       contractionLength(i) = contractions(i)%length
       contractionAngularMoment(i)  = contractions(i)%angularMoment
       contractionNumCartesianOrbital(i) = contractions(i)%numCartesianOrbital
       contractionOwner(i) = contractions(i)%owner
       contractionOrigin(i,1) = contractions(i)%origin(1)
       contractionOrigin(i,2) = contractions(i)%origin(2)
       contractionOrigin(i,3) = contractions(i)%origin(3)
       primNormalizationSize = primNormalizationSize + contractionLength(i)
    end do

    maxLength = 0
    do i=1, numberOfContractions
       maxLength = max(maxLength, contractionLength(i))
    end do
    write(*,*) " Max en Fortran", maxLength
    
    maxNumCartesianOrbital = 0
    do i=1, numberOfContractions
       maxNumCartesianOrbital = max(maxNumCartesianOrbital, contractionNumCartesianOrbital(i))
    end do

    if(allocated(contractionOrbitalExponents)) deallocate(contractionOrbitalExponents)
    if(allocated(contractionCoefficients)) deallocate(contractionCoefficients)
    if(allocated(contractionContNormalization)) deallocate(contractionContNormalization)
    if(allocated(contractionPrimNormalization)) deallocate(contractionPrimNormalization)

    allocate(contractionOrbitalExponents(numberOfContractions, maxLength))
    allocate(contractionCoefficients(numberOfContractions, maxLength))
    allocate(contractionContNormalization(numberOfContractions, maxNumCartesianOrbital))
    allocate(contractionPrimNormalization(primNormalizationSize))

    write(*,*) "En la interfaz"
    auxCounter = 1
    do i=1, numberOfContractions
       do j=1, contractionLength(i)
          contractionOrbitalExponents(i,j) = contractions(i)%orbitalExponents(j)
          write(*,*) contractionOrbitalExponents(i,j)
          contractionPrimNormalization(auxCounter) = contractions(i)%primNormalization(j,1)
          ! write(*,*) contractionPrimNormalization(auxCounter), contractions(i)%primNormalization(j,1)
          contractionCoefficients(i,j) = contractions(i)%contractionCoefficients(j)
          auxCounter = auxCounter + 1
       end do
       if(contractionLength(i)<maxLength) then
          do j=contractionLength(i)+1,maxLength
             contractionOrbitalExponents(i,j) = 0.0_8
             contractionCoefficients(i,j) = 0.0_8
          end do
       end if
       write(*,*) contractionOrbitalExponents(i,j)
    end do

    ! write(*,*) ""
    ! write(*,*) "Constantes de normalizacion"
    do i=1, numberOfContractions
       do j=1, contractionNumCartesianOrbital(i) 
          contractionContNormalization(i,j) = contractions(i)%contNormalization(j)
       end do
       if(contractionNumCartesianOrbital(i)<maxNumCartesianOrbital) then
          do j=contractionNumCartesianOrbital(i)+1,maxNumCartesianOrbital
             contractionContNormalization(i,j) = 0.0_8
          end do
       end if
       ! write(*,*) contractionContNormalization(i,:)
    end do


    call cuda_int_intraspecies(&
         numberOfContractions, &
         maxLength, &
         maxNumCartesianOrbital, &
         primNormalizationSize, &
         contractionId, &
         contractionLength, &
         contractionAngularMoment, &
         contractionNumCartesianOrbital, &
         contractionOwner, &
         contractionOrigin, &
         contractionOrbitalExponents, &
         contractionCoefficients, &
         contractionContNormalization, &
         contractionPrimNormalization &
         )

    
  end subroutine CudintInterface_computeIntraSpecies
  
 !>
  !! @brief  Maneja excepciones de la clase
  subroutine CudintInterface_exception( typeMessage, description, debugDescription)
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

  end subroutine CudintInterface_exception

end module CudintInterface_
