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

  !> @brief the integrals are saved for big records (that reduces the I/O time)
  type, public :: erisStack
     integer*2, allocatable :: a(:)
     integer*2, allocatable :: b(:)
     integer*2, allocatable :: c(:)
     integer*2, allocatable :: d(:)
     real(8), allocatable :: integrals(:)
  end type erisStack
 
  interface

     subroutine cuda_int_intraspecies (&
          numberOfContractions, &
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
          contractionPrimNormalization, &
          contractionIntegrals, &
          contractionIndices) bind (C, name = "cuda_int_intraspecies_")
       use, intrinsic :: iso_c_binding
       implicit none
       integer (c_int) :: numberOfContractions
       integer (c_int) :: maxNumCartesianOrbital
       integer (c_int) :: primNormalizationSize
       integer (c_int) :: contractionId(*)
       integer (c_int) :: contractionLength(*)
       integer (c_int) :: contractionAngularMoment(*)
       integer (c_int) :: contractionNumCartesianOrbital(*)
       integer (c_int) :: contractionOwner(*)
       real (c_double) :: contractionOrigin(*)
       real (c_double) :: contractionOrbitalExponents(*)
       real (c_double) :: contractionCoefficients(*)
       real (c_double) :: contractionContNormalization(numberOfContractions,*)
       real (c_double) :: contractionPrimNormalization(*)
       real (c_double) :: contractionIntegrals(*)
       integer (c_int) :: contractionIndices(*)
     end subroutine cuda_int_intraspecies
    
 
  end interface

  !> @brief Integrals Stack
  type(erisStack), private :: eris
  
contains
  

  subroutine CudintInterface_computeIntraSpecies(specieID)
    implicit none
    
    integer, intent(in) :: specieID
    
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: auxCounter, counter, control

    integer,target :: i, j, k
    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie

    integer, allocatable :: contractionId(:)
    integer, allocatable :: contractionLength(:)
    integer, allocatable :: contractionAngularMoment(:)
    integer, allocatable :: contractionNumCartesianOrbital(:)
    integer, allocatable :: contractionOwner(:)
    real(8), allocatable :: contractionOrigin(:)
    real(8), allocatable :: contractionOrbitalExponents(:)
    real(8), allocatable :: contractionCoefficients(:)
    real(8), allocatable :: contractionContNormalization(:,:)
    real(8), allocatable :: contractionPrimNormalization(:)
    integer :: maxNumCartesianOrbital
    integer :: primNormalizationSize
    integer :: unicIntegrals


    real(8), allocatable :: contractionIntegrals(:)
    integer, allocatable :: contractionIndices(:)
    character(50) :: fileNumber

    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))

    write(fileNumber,*) 1
    fileNumber = trim(adjustl(fileNumber))

   open(UNIT=34,FILE=trim(fileNumber)//trim(MolecularSystem_instance%species(specieID)%name)//".ints", &
        STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='Unformatted')

   !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)

    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

    
    unicIntegrals = ((numberOfContractions*(numberOfContractions+1)/2)+1)*(numberOfContractions*(numberOfContractions+1)/2)/2

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
    allocate(contractionOrigin(numberOfContractions*3))

    primNormalizationSize=0
    do i=1, numberOfContractions
       contractionId(i) = contractions(i)%id
       contractionLength(i) = contractions(i)%length
       contractionAngularMoment(i)  = contractions(i)%angularMoment
       contractionNumCartesianOrbital(i) = contractions(i)%numCartesianOrbital
       contractionOwner(i) = contractions(i)%owner
       contractionOrigin(i*3-2) = contractions(i)%origin(1)
       contractionOrigin(i*3-1) = contractions(i)%origin(2)
       contractionOrigin(i*3) = contractions(i)%origin(3)
       primNormalizationSize = primNormalizationSize + contractionLength(i)
    end do

  
    maxNumCartesianOrbital = 0
    do i=1, numberOfContractions
       maxNumCartesianOrbital = max(maxNumCartesianOrbital, contractionNumCartesianOrbital(i))
    end do

    if(allocated(contractionOrbitalExponents)) deallocate(contractionOrbitalExponents)
    if(allocated(contractionCoefficients)) deallocate(contractionCoefficients)
    if(allocated(contractionContNormalization)) deallocate(contractionContNormalization)
    if(allocated(contractionPrimNormalization)) deallocate(contractionPrimNormalization)

    allocate(contractionOrbitalExponents(primNormalizationSize))
    allocate(contractionCoefficients(primNormalizationSize))
    allocate(contractionContNormalization(numberOfContractions, maxNumCartesianOrbital))
    allocate(contractionPrimNormalization(primNormalizationSize))

    ! write(*,*) "En la interfaz"
    auxCounter = 1
    do i=1, numberOfContractions
       do j=1, contractionLength(i)
          contractionOrbitalExponents(auxCounter) = contractions(i)%orbitalExponents(j)
          contractionPrimNormalization(auxCounter) = contractions(i)%primNormalization(j,1)
          contractionCoefficients(auxCounter) = contractions(i)%contractionCoefficients(j)
          ! write(*,*) contractions(i)%orbitalExponents(j)
          ! write(*,*) contractionPrimNormalization(auxCounter), contractions(i)%primNormalization(j,1)
          auxCounter = auxCounter + 1
       end do
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

    if(allocated(contractionIntegrals)) deallocate(contractionIntegrals)
    allocate(contractionIntegrals(unicIntegrals))

    if(allocated(contractionIndices)) deallocate(contractionIndices)
    allocate(contractionIndices(unicIntegrals*4))


    call cuda_int_intraspecies(&
         numberOfContractions, &
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
         contractionPrimNormalization, &
         contractionIntegrals, &
         contractionIndices &
         )

    auxCounter = 0
    counter = 0
    do i=1, unicIntegrals

       if(abs(contractionIntegrals(i)) > 1.0D-10) then

          auxCounter = auxCounter + 1
          counter = counter + 1

          eris%a(counter) = contractionIndices(i*4 - 3)
          eris%b(counter) = contractionIndices(i*4 - 2)
          eris%c(counter) = contractionIndices(i*4 - 1)
          eris%d(counter) = contractionIndices(i*4)
          eris%integrals(counter) = contractionIntegrals(i)
       end if

       if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

          write(34) &
               eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

          counter = 0

       end if
    end do

    counter = counter + 1 
    eris%a(counter) = -1
    eris%b(counter) = -1
    eris%c(counter) = -1
    eris%d(counter) = -1
    eris%integrals(counter) = 0.0_8

    write(34) &
         eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    close(34)

    write(6,"(A,I12,A,A)") " Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
         trim(MolecularSystem_instance%species(specieID)%name)

    
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
