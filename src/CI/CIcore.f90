 module CIcore_
  use Exception_
  use Matrix_
  use Vector_
  use MolecularSystem_
  use Configuration_
  use MolecularSystem_
  use String_
  use IndexMap_
  use InputCI_
  use omp_lib
  implicit none

 type, public :: CIcore
     logical :: isInstanced
     integer :: numberOfSpecies
     type(matrix) :: hamiltonianMatrix
     type(ivector8) :: auxIndexCIMatrix
     type(matrix) :: eigenVectors
     type(matrix) :: initialEigenVectors
     type(vector8) :: initialEigenValues
     integer(8) :: numberOfConfigurations
     integer :: nproc
     integer :: numberOfQuantumSpecies
     type(ivector) :: numberOfCoreOrbitals
     type(ivector) :: numberOfOccupiedOrbitals
     type(ivector) :: numberOfOrbitals
     type(vector) ::  numberOfSpatialOrbitals2 
     type(vector8) :: eigenvalues
     type(vector) :: groundStateEnergies
     type(vector) :: DDCISDTiming
     type(vector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     type(matrix), allocatable :: fourCenterIntegrals(:,:)
     type(matrix), allocatable :: twoCenterIntegrals(:)
     type(imatrix8), allocatable :: twoIndexArray(:)
     type(imatrix8), allocatable :: fourIndexArray(:)
     type(imatrix), allocatable :: strings(:) !! species, conf, occupations. index for occupied orbitals, e.g. 1 2 5 6
     type(imatrix1), allocatable :: orbitals(:) !! species, conf, occupations. array with 1 for occupied and 0 unoccupied orb, e.g. 1 1 0 0 1 1
     integer, allocatable :: sumstrings(:) !! species
     type(ivector), allocatable :: auxstring(:,:) !! species, occupations
     type(ivector8), allocatable :: numberOfStrings(:) !! species, excitation level, number of strings
     type(ivector8), allocatable :: numberOfStrings2(:) !! species, excitation level, number of strings

     !! species, threads
     type(imatrix), allocatable :: couplingMatrix(:,:)
     type(Vector), allocatable :: couplingMatrixEnergyOne(:,:)
!     type(matrix), allocatable :: couplingMatrixEnergyTwo(:)
     type(ivector), allocatable :: couplingMatrixFactorOne(:,:)
     type(ivector), allocatable :: couplingMatrixOrbOne(:,:)
     type(imatrix), allocatable :: nCouplingOneTwo(:,:)
     type(imatrix), allocatable :: nCouplingSize(:,:)

     type(ivector1), allocatable :: couplingOrderList(:,:)
     type(ivector1), allocatable :: couplingOrderIndex(:,:)

     integer, allocatable :: ciOrderList(:,:)
     integer, allocatable :: auxciOrderList(:)
     integer :: sizeCiOrderList
     integer(8), allocatable :: ciOrderSize1(:,:)
     integer(8), allocatable :: ciOrderSize2(:,:)
     integer(4), allocatable :: allIndexConf(:,:) !! species, total number of configurations

     integer :: ncouplingOrderOne
     integer :: ncouplingOrderTwo
     integer :: ncouplingOrderTwoDiff

     type(imatrix) :: auxConfigurations !! species, configurations for initial hamiltonian
     type(imatrix) :: coreConfigurations !! species, configurations for core SCI space
     type(imatrix) :: targetConfigurations !! species, configurations for target SCI space
     type(imatrix) :: fullConfigurations !! species, configurations for target SCI space
     type(imatrix) :: coreConfigurationsLevel !! species, configurations for CI level of core SCI space
     type(imatrix) :: targetConfigurationsLevel !! species, configurations for CI level target SCI space
     type(imatrix) :: fullConfigurationsLevel !! species, configurations for CI level target SCI space

     type(configuration), allocatable :: configurations(:)
     integer(2), allocatable :: auxconfs(:,:,:) ! nconf, species, occupation
     type (Vector8) :: diagonalHamiltonianMatrix
     type (Vector8) :: diagonalHamiltonianMatrix2
     real(8) :: totalEnergy
     integer, allocatable :: totalNumberOfContractions(:)
     integer, allocatable :: occupationNumber(:)
     integer, allocatable :: recursionVector1(:)
     integer, allocatable :: recursionVector2(:)
     integer, allocatable :: CILevel(:)
     integer, allocatable :: pindexConf(:,:) !! save previous configuration to avoid unneccesary calculations

     integer :: maxCILevel
     type (Matrix) :: initialHamiltonianMatrix
     type (Matrix) :: initialHamiltonianMatrix2
     character(20) :: level
     real(8) :: timeA(7)
     real(8) :: timeB(7)

  end type CIcore

  type, public :: HartreeFock
        real(8) :: totalEnergy
        real(8) :: puntualInteractionEnergy
        type(matrix) :: coefficientsofcombination 
        type(matrix) :: HcoreMatrix 
  end type HartreeFock

  integer, allocatable :: Conf_occupationNumber(:)
  type(HartreeFock) :: HartreeFock_instance
  type(CIcore) :: CIcore_instance

  public :: &
       CIcore_constructor


contains

  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine CIcore_constructor(level)
    implicit none
    character(*) :: level

    integer :: numberOfSpecies
    integer :: i,j,k,l,m,n,p,q,cc,r,s,el, nproc
    integer(8) :: c
    integer :: ma,mb,mc,md,me,pa,pb,pc,pd,pe
    integer :: isLambdaEqual1,lambda,otherlambda
    type(vector) :: occupiedCode
    type(vector) :: unoccupiedCode
    real(8) :: totalEnergy

    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: nameOfSpecie
    integer :: numberOfContractions
    character(50) :: arguments(2)

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])

    CIcore_instance%numberOfQuantumSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies
    CIcore_instance%numberOfSpecies = numberOfSpecies


    do i=1, numberOfSpecies
        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecies( i ) )
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )

        arguments(2) = nameOfSpecie
        arguments(1) = "HCORE"
        HartreeFock_instance%HcoreMatrix  = &
                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "COEFFICIENTS"
        HartreeFock_instance%coefficientsofcombination = &
                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
    end do

    CIcore_instance%isInstanced=.true.
    CIcore_instance%level=level
    CIcore_instance%numberOfConfigurations=0

    call Vector_constructorInteger (CIcore_instance%numberOfCoreOrbitals, numberOfSpecies)
    call Vector_constructorInteger (CIcore_instance%numberOfOccupiedOrbitals, numberOfSpecies)
    call Vector_constructorInteger (CIcore_instance%numberOfOrbitals, numberOfSpecies)
    call Vector_constructor (CIcore_instance%lambda, numberOfSpecies)
    call Vector_constructor (CIcore_instance%numberOfSpatialOrbitals2, numberOfSpecies)

    CIcore_instance%nproc = omp_get_max_threads()

    if  ( allocated (CIcore_instance%strings ) ) &
    deallocate ( CIcore_instance%strings )
    allocate ( CIcore_instance%strings ( numberOfSpecies ) )

    if  ( allocated (CIcore_instance%orbitals ) ) &
    deallocate ( CIcore_instance%orbitals )
    allocate ( CIcore_instance%orbitals ( numberOfSpecies ) )

    if  ( allocated (CIcore_instance%auxstring ) ) &
    deallocate ( CIcore_instance%auxstring )
    allocate ( CIcore_instance%auxstring ( CIcore_instance%nproc, numberOfSpecies ) )

    if  ( allocated (CIcore_instance%couplingMatrix ) ) &
    deallocate ( CIcore_instance%couplingMatrix )
    allocate ( CIcore_instance%couplingMatrix ( numberOfSpecies, CIcore_instance%nproc ) )

    if  ( allocated (CIcore_instance%couplingMatrixEnergyOne ) ) &
    deallocate ( CIcore_instance%couplingMatrixEnergyOne )
    allocate ( CIcore_instance%couplingMatrixEnergyOne ( numberOfSpecies, CIcore_instance%nproc  ) )

    if  ( allocated (CIcore_instance%couplingMatrixFactorOne ) ) &
    deallocate ( CIcore_instance%couplingMatrixFactorOne )
    allocate ( CIcore_instance%couplingMatrixFactorOne ( numberOfSpecies, CIcore_instance%nproc  ) )

    if  ( allocated (CIcore_instance%couplingMatrixOrbOne ) ) &
    deallocate ( CIcore_instance%couplingMatrixOrbOne )
    allocate ( CIcore_instance%couplingMatrixOrbOne ( numberOfSpecies, CIcore_instance%nproc  ) )

    if  ( allocated (CIcore_instance%nCouplingOneTwo ) ) &
    deallocate ( CIcore_instance%nCouplingOneTwo )
    allocate ( CIcore_instance%nCouplingOneTwo ( numberOfSpecies, CIcore_instance%nproc  ) )

    if  ( allocated (CIcore_instance%nCouplingSize ) ) &
    deallocate ( CIcore_instance%nCouplingSize )
    allocate ( CIcore_instance%nCouplingSize ( numberOfSpecies, CIcore_instance%nproc  ) )

    if  ( allocated (CIcore_instance%numberOfStrings ) ) &
    deallocate ( CIcore_instance%numberOfStrings )
    allocate ( CIcore_instance%numberOfStrings ( numberOfSpecies ) )

    if  ( allocated (CIcore_instance%numberOfStrings2 ) ) &
    deallocate ( CIcore_instance%numberOfStrings2 )
    allocate ( CIcore_instance%numberOfStrings2 ( numberOfSpecies ) )

    if  ( allocated (CIcore_instance%sumstrings ) ) &
    deallocate ( CIcore_instance%sumstrings )
    allocate ( CIcore_instance%sumstrings ( numberOfSpecies ) )

    if ( allocated ( CIcore_instance%totalNumberOfContractions ) ) &
    deallocate ( CIcore_instance%totalNumberOfContractions ) 
    allocate ( CIcore_instance%totalNumberOfContractions (numberOfSpecies ) )

    if ( allocated ( CIcore_instance%occupationNumber ) ) &
    deallocate ( CIcore_instance%occupationNumber ) 
    allocate ( CIcore_instance%occupationNumber (numberOfSpecies ) )

    if ( allocated ( CIcore_instance%recursionVector1 ) ) &
    deallocate ( CIcore_instance%recursionVector1 ) 
    allocate ( CIcore_instance%recursionVector1 (numberOfSpecies ) )

    if ( allocated ( CIcore_instance%recursionVector2 ) ) &
    deallocate ( CIcore_instance%recursionVector2 ) 
    allocate ( CIcore_instance%recursionVector2 (numberOfSpecies ) )

    if ( allocated ( CIcore_instance%CILevel) ) &
    deallocate ( CIcore_instance%CILevel ) 
    allocate ( CIcore_instance%CILevel (numberOfSpecies ) )

    if ( allocated ( CIcore_instance%pindexConf) ) &
    deallocate ( CIcore_instance%pindexConf ) 
    allocate ( CIcore_instance%pindexConf (numberOfSpecies, CIcore_instance%nproc ) )

    if ( allocated ( Conf_occupationNumber ) ) &
    deallocate ( Conf_occupationNumber ) 
    allocate ( Conf_occupationNumber (numberOfSpecies ) )


    CIcore_instance%recursionVector1 = 1
    CIcore_instance%recursionVector2 = 0

    CIcore_instance%recursionVector1(numberOfSpecies) = 0
    CIcore_instance%recursionVector2(numberOfSpecies) = 1

    CIcore_instance%pindexConf = 0

    do i=1, numberOfSpecies
       !! We are working in spin orbitals not in spatial orbitals!
       CIcore_instance%lambda%values(i) = MolecularSystem_getLambda( i )
       CIcore_instance%numberOfCoreOrbitals%values(i) = 0
       CIcore_instance%numberOfOccupiedOrbitals%values(i) = int (MolecularSystem_getOcupationNumber( i )* &
                                                                              CIcore_instance%lambda%values(i))
       CIcore_instance%numberOfOrbitals%values(i) = MolecularSystem_getTotalNumberOfContractions( i )* &
                                                                      CIcore_instance%lambda%values(i) 
       CIcore_instance%numberOfSpatialOrbitals2%values(i) = MolecularSystem_getTotalNumberOfContractions( i )
       CIcore_instance%numberOfSpatialOrbitals2%values(i) = &
         CIcore_instance%numberOfSpatialOrbitals2%values(i) *  ( &
         CIcore_instance%numberOfSpatialOrbitals2%values(i) + 1 ) / 2

      
       CIcore_instance%totalNumberOfContractions( i ) = MolecularSystem_getTotalNumberOfContractions( i )
       CIcore_instance%occupationNumber( i ) = int( MolecularSystem_instance%species(i)%ocupationNumber )
       Conf_occupationNumber( i ) =  MolecularSystem_instance%species(i)%ocupationNumber


      !! Take the active space from input
      if ( InputCI_Instance(i)%coreOrbitals /= 0 ) then
       CIcore_instance%numberOfCoreOrbitals%values(i) = InputCI_Instance(i)%coreOrbitals 
      end if
      if ( InputCI_Instance(i)%activeOrbitals /= 0 ) then
        CIcore_instance%numberOfOrbitals%values(i) = InputCI_Instance(i)%activeOrbitals * &
                                    CIcore_instance%lambda%values(i) + &
                                    CIcore_instance%numberOfCoreOrbitals%values(i)
      end if

       !!Uneven occupation number = alpha
       !!Even occupation number = beta     
    end do

    call Configuration_globalConstructor()

    close(wfnUnit)

  end subroutine CIcore_constructor

recursive  function CIcore_gatherConfRecursion(s, numberOfSpecies, indexConf, c, cilevel ) result (os)
    implicit none

    integer(8) :: a,b,c,cc,d
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies
    integer :: os,is
    integer :: size1, size2
    integer(8) :: indexConf(:)
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer :: ssize
    integer :: cilevel(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a
        os = CIcore_gatherConfRecursion( is, numberOfSpecies, indexConf, c, cilevel )
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        CIcore_instance%allIndexConf(:,c) = indexConf

      end do
    end if

  end function CIcore_gatherConfRecursion


  function CIcore_getIndex ( indexConf ) result ( output )
    implicit none
    integer(8) :: indexConf(:)
    integer(8) :: output, ssize
    integer :: i,j, numberOfSpecies

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies
    output = 0 
     !! simplify!!
    do i = 1, numberOfSpecies
      ssize = 1 
      do j = i + 1, numberOfSpecies
        ssize = ssize * CIcore_instance%sumstrings(j)
        !ssize = ssize * sum(CIcore_instance%numberOfStrings(j)%values(1:2))
      end do
      output = output + ( indexConf(i) - 1 ) * ssize
    end do
    output = output + 1

  end function CIcore_getIndex

end module CIcore_

