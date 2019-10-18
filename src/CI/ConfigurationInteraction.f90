!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	  UNIVERSIDAD NACIONAL DE COLOMBIA"
!!	  PROF. ANDRES REYES GROUP"
!!	  http://www.qcc.unal.edu.co"
!!	
!!	  UNIVERSIDAD DE GUADALAJARA"
!!	  PROF. ROBERTO FLORES GROUP"
!!	  http://www.cucei.udg.mx/~robertof"
!!	
!!	AUTHORS
!!		E.F. POSADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		S.A. GONZALEZ. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		F.S. MONCADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		J. ROMERO. UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!	CONTRIBUTORS
!!		N.F.AGUIRRE. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		GABRIEL MERINO. UNIVERSIDAD DE GUANAJUATO
!!   		J.A. CHARRY UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************
                
module ConfigurationInteraction_
  use Exception_
  use Matrix_
  use Vector_
!  use MolecularSystem_
  use Configuration_
  use ReadTransformedIntegrals_
  use MolecularSystem_
  use String_
  use IndexMap_
  use InputCI_
  use omp_lib
  use ArpackInterface_
  implicit none
      
  !>
  !! @brief Configuration Interaction Module, works in spin orbitals
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 07-24-12
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 07-24-12 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# description.
  !!   - <tt> 07-09-16 </tt>: Jorge Charry ( jacharrym@unal.edu.co )
  !!        -# Add CIS, and Fix CISD.
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<

  type, public :: ConfigurationInteraction
     logical :: isInstanced
     integer :: numberOfSpecies
     type(matrix) :: hamiltonianMatrix
     type(ivector8) :: auxIndexCIMatrix
     type(matrix) :: eigenVectors
     type(matrix) :: initialEigenVectors
     type(vector8) :: initialEigenValues
     integer(8) :: numberOfConfigurations
     integer :: nproc
     type(ivector) :: numberOfCoreOrbitals
     type(ivector) :: numberOfOccupiedOrbitals
     type(ivector) :: numberOfOrbitals
     type(vector) ::  numberOfSpatialOrbitals2 
     type(vector8) :: eigenvalues
     type(vector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     type(matrix), allocatable :: fourCenterIntegrals(:,:)
     type(matrix), allocatable :: twoCenterIntegrals(:)
     type(imatrix8), allocatable :: twoIndexArray(:)
     type(imatrix8), allocatable :: fourIndexArray(:)
     type(imatrix), allocatable :: strings(:) !! species, conf, occupations
     type(imatrix1), allocatable :: orbitals(:) !! species, conf, occupations
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
     integer, allocatable :: pindexConf(:,:)
     integer :: maxCILevel
     type (Matrix) :: initialHamiltonianMatrix
     type (Matrix) :: initialHamiltonianMatrix2
     character(20) :: level
     real(8) :: timeA(7)
     real(8) :: timeB(7)

  end type ConfigurationInteraction

  type, public :: HartreeFock
        real(8) :: totalEnergy
        real(8) :: puntualInteractionEnergy
        type(matrix) :: coefficientsofcombination 
        type(matrix) :: HcoreMatrix 
  end type HartreeFock
  
  integer, allocatable :: Conf_occupationNumber(:)
  type(ConfigurationInteraction) :: ConfigurationInteraction_instance
  type(HartreeFock) :: HartreeFock_instance

  public :: &
       ConfigurationInteraction_constructor, &
       ConfigurationInteraction_destructor, &
       ConfigurationInteraction_getTotalEnergy, &
       ConfigurationInteraction_run, &
       ConfigurationInteraction_showEigenVectors, &
       ConfigurationInteraction_densityMatrices, &
       ConfigurationInteraction_show

  private

contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine ConfigurationInteraction_constructor(level)
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

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ConfigurationInteraction_instance%numberOfSpecies = numberOfSpecies


    do i=1, numberOfSpecies
        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
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

    ConfigurationInteraction_instance%isInstanced=.true.
    ConfigurationInteraction_instance%level=level
    ConfigurationInteraction_instance%numberOfConfigurations=0

    call Vector_constructorInteger (ConfigurationInteraction_instance%numberOfCoreOrbitals, numberOfSpecies)
    call Vector_constructorInteger (ConfigurationInteraction_instance%numberOfOccupiedOrbitals, numberOfSpecies)
    call Vector_constructorInteger (ConfigurationInteraction_instance%numberOfOrbitals, numberOfSpecies)
    call Vector_constructor (ConfigurationInteraction_instance%lambda, numberOfSpecies)
    call Vector_constructor (ConfigurationInteraction_instance%numberOfSpatialOrbitals2, numberOfSpecies)

    ConfigurationInteraction_instance%nproc = omp_get_max_threads()

    if  ( allocated (ConfigurationInteraction_instance%strings ) ) &
    deallocate ( ConfigurationInteraction_instance%strings )
    allocate ( ConfigurationInteraction_instance%strings ( numberOfSpecies ) )

    if  ( allocated (ConfigurationInteraction_instance%orbitals ) ) &
    deallocate ( ConfigurationInteraction_instance%orbitals )
    allocate ( ConfigurationInteraction_instance%orbitals ( numberOfSpecies ) )

    if  ( allocated (ConfigurationInteraction_instance%auxstring ) ) &
    deallocate ( ConfigurationInteraction_instance%auxstring )
    allocate ( ConfigurationInteraction_instance%auxstring ( ConfigurationInteraction_instance%nproc, numberOfSpecies ) )

    if  ( allocated (ConfigurationInteraction_instance%couplingMatrix ) ) &
    deallocate ( ConfigurationInteraction_instance%couplingMatrix )
    allocate ( ConfigurationInteraction_instance%couplingMatrix ( numberOfSpecies, ConfigurationInteraction_instance%nproc ) )

    if  ( allocated (ConfigurationInteraction_instance%couplingMatrixEnergyOne ) ) &
    deallocate ( ConfigurationInteraction_instance%couplingMatrixEnergyOne )
    allocate ( ConfigurationInteraction_instance%couplingMatrixEnergyOne ( numberOfSpecies, ConfigurationInteraction_instance%nproc  ) )

    if  ( allocated (ConfigurationInteraction_instance%couplingMatrixFactorOne ) ) &
    deallocate ( ConfigurationInteraction_instance%couplingMatrixFactorOne )
    allocate ( ConfigurationInteraction_instance%couplingMatrixFactorOne ( numberOfSpecies, ConfigurationInteraction_instance%nproc  ) )

    if  ( allocated (ConfigurationInteraction_instance%couplingMatrixOrbOne ) ) &
    deallocate ( ConfigurationInteraction_instance%couplingMatrixOrbOne )
    allocate ( ConfigurationInteraction_instance%couplingMatrixOrbOne ( numberOfSpecies, ConfigurationInteraction_instance%nproc  ) )

    if  ( allocated (ConfigurationInteraction_instance%nCouplingOneTwo ) ) &
    deallocate ( ConfigurationInteraction_instance%nCouplingOneTwo )
    allocate ( ConfigurationInteraction_instance%nCouplingOneTwo ( numberOfSpecies, ConfigurationInteraction_instance%nproc  ) )

    if  ( allocated (ConfigurationInteraction_instance%nCouplingSize ) ) &
    deallocate ( ConfigurationInteraction_instance%nCouplingSize )
    allocate ( ConfigurationInteraction_instance%nCouplingSize ( numberOfSpecies, ConfigurationInteraction_instance%nproc  ) )

    if  ( allocated (ConfigurationInteraction_instance%numberOfStrings ) ) &
    deallocate ( ConfigurationInteraction_instance%numberOfStrings )
    allocate ( ConfigurationInteraction_instance%numberOfStrings ( numberOfSpecies ) )

    if  ( allocated (ConfigurationInteraction_instance%numberOfStrings2 ) ) &
    deallocate ( ConfigurationInteraction_instance%numberOfStrings2 )
    allocate ( ConfigurationInteraction_instance%numberOfStrings2 ( numberOfSpecies ) )

    if  ( allocated (ConfigurationInteraction_instance%sumstrings ) ) &
    deallocate ( ConfigurationInteraction_instance%sumstrings )
    allocate ( ConfigurationInteraction_instance%sumstrings ( numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%totalNumberOfContractions ) ) &
    deallocate ( ConfigurationInteraction_instance%totalNumberOfContractions ) 
    allocate ( ConfigurationInteraction_instance%totalNumberOfContractions (numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%occupationNumber ) ) &
    deallocate ( ConfigurationInteraction_instance%occupationNumber ) 
    allocate ( ConfigurationInteraction_instance%occupationNumber (numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%recursionVector1 ) ) &
    deallocate ( ConfigurationInteraction_instance%recursionVector1 ) 
    allocate ( ConfigurationInteraction_instance%recursionVector1 (numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%recursionVector2 ) ) &
    deallocate ( ConfigurationInteraction_instance%recursionVector2 ) 
    allocate ( ConfigurationInteraction_instance%recursionVector2 (numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%CILevel) ) &
    deallocate ( ConfigurationInteraction_instance%CILevel ) 
    allocate ( ConfigurationInteraction_instance%CILevel (numberOfSpecies ) )

    if ( allocated ( ConfigurationInteraction_instance%pindexConf) ) &
    deallocate ( ConfigurationInteraction_instance%pindexConf ) 
    allocate ( ConfigurationInteraction_instance%pindexConf (numberOfSpecies, ConfigurationInteraction_instance%nproc ) )

    if ( allocated ( Conf_occupationNumber ) ) &
    deallocate ( Conf_occupationNumber ) 
    allocate ( Conf_occupationNumber (numberOfSpecies ) )


    ConfigurationInteraction_instance%recursionVector1 = 1
    ConfigurationInteraction_instance%recursionVector2 = 0

    ConfigurationInteraction_instance%recursionVector1(numberOfSpecies) = 0
    ConfigurationInteraction_instance%recursionVector2(numberOfSpecies) = 1

    ConfigurationInteraction_instance%pindexConf = 0

    do i=1, numberOfSpecies
       !! We are working in spin orbitals not in spatial orbitals!
       ConfigurationInteraction_instance%lambda%values(i) = MolecularSystem_getLambda( i )
       ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) = 0
       ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) = int (MolecularSystem_getOcupationNumber( i )* &
                                                                              ConfigurationInteraction_instance%lambda%values(i))
       ConfigurationInteraction_instance%numberOfOrbitals%values(i) = MolecularSystem_getTotalNumberOfContractions( i )* &
                                                                      ConfigurationInteraction_instance%lambda%values(i) 
       ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) = MolecularSystem_getTotalNumberOfContractions( i )
       ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) = &
         ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) *  ( &
         ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(i) + 1 ) / 2

      
       ConfigurationInteraction_instance%totalNumberOfContractions( i ) = MolecularSystem_getTotalNumberOfContractions( i )
       ConfigurationInteraction_instance%occupationNumber( i ) = int( MolecularSystem_instance%species(i)%ocupationNumber )
       Conf_occupationNumber( i ) =  MolecularSystem_instance%species(i)%ocupationNumber


      !! Take the active space from input
      if ( InputCI_Instance(i)%coreOrbitals /= 0 ) then
       ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) = InputCI_Instance(i)%coreOrbitals 
      end if
      if ( InputCI_Instance(i)%activeOrbitals /= 0 ) then
        ConfigurationInteraction_instance%numberOfOrbitals%values(i) = InputCI_Instance(i)%activeOrbitals * &
                                    ConfigurationInteraction_instance%lambda%values(i) + &
                                    ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i)
      end if

       !!Uneven occupation number = alpha
       !!Even occupation number = beta     
    end do


    call Configuration_globalConstructor()

    close(wfnUnit)

  end subroutine ConfigurationInteraction_constructor

  subroutine ConfigurationInteraction_buildStrings()
    implicit none

    integer(8) :: a,b,c,c1,c2,aa,d
    integer :: ci, oci, cilevel,maxcilevel
    integer :: u,uu,vv, p, nn,z
    integer :: i,j
    integer :: numberOfSpecies, auxnumberOfSpecies,s
    type(ivector) :: order
    integer(8) :: ssize
    real(8) :: timeA, timeB
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( allocated( occupiedCode ) ) deallocate( occupiedCode )
    allocate (occupiedCode ( numberOfSpecies ) )
    if ( allocated( unoccupiedCode ) ) deallocate( unoccupiedCode )
    allocate (unoccupiedCode ( numberOfSpecies ) )

    call Vector_constructorInteger (order, numberOfSpecies, 0 )
    order%values = 0

    s = 0
    do i = 1, numberOfSpecies

      call Vector_constructorInteger8 (ConfigurationInteraction_instance%numberOfStrings(i), &
        int(ConfigurationInteraction_instance%CILevel(i) + 1,8), 0_8)

      ConfigurationInteraction_instance%numberOfStrings(i)%values(1) = 1 !! ground

      write (*,"(A,A)") "  ", MolecularSystem_getNameOfSpecie(i)

      do cilevel = 1,ConfigurationInteraction_instance%CILevel(i) 

        call Vector_constructor (occupiedCode(i), cilevel, real(ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i),8) )
        call Vector_constructor (unoccupiedCode(i), cilevel, 0.0_8)

        unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)  ! it's also a lower bound in a for loop

        if ( cilevel <= ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) ) then

          !! just get the number of strings...
          ci = 0 
          oci = ConfigurationInteraction_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel)

          write (*,"(A,I4,I8)") "    ", cilevel, ConfigurationInteraction_instance%numberOfStrings(i)%values(cilevel+1)

        end if
      end do
      write (*,"(A,I8)") "  Total:", sum(ConfigurationInteraction_instance%numberOfStrings(i)%values)
      write (*,"(A)") ""

      !! allocate the strings arrays
      if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) > 0 ) then
        call Matrix_constructorInteger( ConfigurationInteraction_instance%strings(i), &
          int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i),8), &
          sum(ConfigurationInteraction_instance%numberOfStrings(i)%values), int(0,4))

        call Matrix_constructorInteger1( ConfigurationInteraction_instance%orbitals(i), &
          int(ConfigurationInteraction_instance%numberOfOrbitals%values(i),8), &
          sum(ConfigurationInteraction_instance%numberOfStrings(i)%values), 0_1)

      else
        call Matrix_constructorInteger( ConfigurationInteraction_instance%strings(i), &
         1_8, 1_8, int(0,4))
        call Matrix_constructorInteger1( ConfigurationInteraction_instance%orbitals(i), &
         1_8, 1_8, 0_1)

      end if

      !! zero, build the reference
      call Vector_constructorInteger (order, numberOfSpecies, 0 )

      call Vector_constructor (occupiedCode(i), 1, 0.0_8) !! initialize in zero
      call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)

      c = 0 
      c = c + 1
      call Configuration_constructorB(ConfigurationInteraction_instance%strings(i), ConfigurationInteraction_instance%orbitals(i), &
                                      occupiedCode, unoccupiedCode, i, c, order)
      
      !! now build the strings
      do cilevel = 1,ConfigurationInteraction_instance%CILevel(i) 

        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i) = cilevel

        call Vector_constructor (occupiedCode(i), cilevel, real(ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i),8) )
        call Vector_constructor (unoccupiedCode(i), cilevel, 0.0_8)

        unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)  ! it's also a lower bound in a for loop

        if ( cilevel <= ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) ) then

          !! recursion to build the strings
          ci = 0 
          oci = ConfigurationInteraction_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c)

        end if
      end do

    end do

    !! useful array  
    do i = 1, numberOfSpecies
      ConfigurationInteraction_instance%sumstrings(i) = sum(ConfigurationInteraction_instance%numberOfStrings(i)%values)
    end do

    !! useful array, save the total number of string for a previous CI level. 
    do i = 1, numberOfSpecies
      call Vector_constructorInteger8 (ConfigurationInteraction_instance%numberOfStrings2(i), &
        int(size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) + 1,8), 0_8)

      ssize = 0
      do j = 1, size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) !
        ssize = ssize + ConfigurationInteraction_instance%numberOfStrings(i)%values(j)
        ConfigurationInteraction_instance%numberOfStrings2(i)%values(j+1) = ssize 
      end do
      ConfigurationInteraction_instance%numberOfStrings2(i)%values(1) = 0
    end do
      

  end subroutine ConfigurationInteraction_buildStrings

!! This is just to get the total number of strings...

recursive  function ConfigurationInteraction_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ici, cilevel ) result (oci)
    implicit none

    integer :: i, numberOfSpecies
    integer :: ci, ici, oci, cilevel
    integer :: m, a
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)

    ci = ici + 1

    if ( ci == 1 .and. ci < cilevel ) then ! first
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = ConfigurationInteraction_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel )
        end do
        unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    else if ( ci > 1 .and. ci < cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = ConfigurationInteraction_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel )
        end do
      end do

    else if ( ci == 1 .and. ci == cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          ConfigurationInteraction_instance%numberOfStrings(i)%values(ci+1) = &
            ConfigurationInteraction_instance%numberOfStrings(i)%values(ci+1) + 1
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do

    else !final

      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          ConfigurationInteraction_instance%numberOfStrings(i)%values(ci+1) = &
            ConfigurationInteraction_instance%numberOfStrings(i)%values(ci+1) + 1
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    end if

  end function ConfigurationInteraction_buildStringsRecursion

!! and this is for building the strings
recursive  function ConfigurationInteraction_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, &
                     ici, cilevel, order, c ) result (oci)
    implicit none

    integer :: i, numberOfSpecies
    integer :: ci, ici, oci, cilevel
    integer(8) :: c
    integer :: m, a
    type(ivector) :: order
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)

    ci = ici + 1

    if ( ci == 1 .and. ci < cilevel ) then ! first
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = ConfigurationInteraction_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c )
        end do
        unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    else if ( ci > 1 .and. ci < cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = ConfigurationInteraction_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c )
        end do
      end do

    else if ( ci == 1 .and. ci == cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a

          c = c + 1
          call Configuration_constructorB(ConfigurationInteraction_instance%strings(i), ConfigurationInteraction_instance%orbitals(i), &
                                          occupiedCode, unoccupiedCode, i, c, order)
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do

    else !final

      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(ConfigurationInteraction_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          c = c + 1
          call Configuration_constructorB(ConfigurationInteraction_instance%strings(i), ConfigurationInteraction_instance%orbitals(i), &
                                          occupiedCode, unoccupiedCode, i, c, order)
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    end if


  end function ConfigurationInteraction_buildStringsRecursion2

!! Build the CI table with all combinations of excitations between quantum species.

  subroutine ConfigurationInteraction_buildCIOrderList()
    implicit none

    integer :: c
    integer :: i,j, u,v
    integer :: ci, ii, jj
    integer(8) :: output, auxsize
    integer :: numberOfSpecies, auxnumberOfSpecies,s
    integer(1) :: coupling
    real(8) :: timeA, timeB
    integer :: ncouplingOrderOne
    integer :: ncouplingOrderTwo
    logical :: includecilevel, same
    integer(8) :: ssize, auxssize
    integer, allocatable :: cilevel(:), auxcilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! Allocate size considering all possible combinations, FCI.
    ssize = 1 
    do i = 1, numberOfSpecies
       ssize = ssize * (ConfigurationInteraction_instance%CILevel(i) + 1)
    end do

    allocate ( ConfigurationInteraction_instance%ciOrderList( ssize, numberOfSpecies ) ) 
    allocate ( ConfigurationInteraction_instance%ciOrderSize1( ssize, numberOfSpecies ) ) 
    allocate ( ConfigurationInteraction_instance%ciOrderSize2( ssize, numberOfSpecies ) ) 
    allocate ( ConfigurationInteraction_instance%auxciOrderList( ssize ) ) 

    ConfigurationInteraction_instance%ciOrderList = 0
    ConfigurationInteraction_instance%auxciOrderList = 0

    ConfigurationInteraction_instance%ciOrderSize1 = -1 !! I have reasons... -1 for all species except the last one
    ConfigurationInteraction_instance%ciOrderSize2 = 1 !! and 1 for the last species

    ConfigurationInteraction_instance%sizeCiOrderList = 0

    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( auxciLevel ( numberOfSpecies ) )
    ciLevel = 0
    auxciLevel = 0
    s = 0
    c = 0
    !! Search which combinations of excitations satifies the desired CI level.
    auxnumberOfSpecies = ConfigurationInteraction_buildCIOrderRecursion( s, numberOfSpecies, c, cilevel )


    !! Print list
    write (6,"(T2,A)") "--------------------------"
    write (6,"(T2,A)") "CI level \ Species"
    write (6,"(T2,A)") "--------------------------"
    do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
      do i = 1, numberOfSpecies
        write (6,"(T2,I4)",advance="no") ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(u), i)
      end do
      write (6,"(A)") ""
    end do
    write (6,"(T2,A)") "--------------------------"

    !! Calculates the three required factors in order to get the position of any given configuration.
    !! position = S1 + (indexConf(i,u) - numberOfStrings2(i) -1 )*S2(i,u)  
    !! i: speciesID, u: cilevelID

    !! Factor S1
    ssize = 0
    do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(u), :)

      ssize = 0
      do v = 1,  u-1

        auxcilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(v), :)
        auxnumberOfSpecies = ConfigurationInteraction_getIndexSize(0, ssize, auxcilevel) 

      end do

      ConfigurationInteraction_instance%ciOrderSize1(ConfigurationInteraction_instance%auxciOrderList(u),:) = -1
      ConfigurationInteraction_instance%ciOrderSize1(ConfigurationInteraction_instance%auxciOrderList(u),numberOfSpecies) = ssize !!just the last

    end do

    !! Factor S2
    do i = 1, numberOfSpecies-1
      do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

        cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(u), :)
        ssize = 1
        do j = i+1, numberOfSpecies
          ssize = ssize * ConfigurationInteraction_instance%numberOfStrings(j)%values(cilevel(j)+1)
        end do

        ConfigurationInteraction_instance%ciOrderSize2(ConfigurationInteraction_instance%auxciOrderList(u),i) = ssize

      end do
    end do

    ConfigurationInteraction_instance%ciOrderSize2(:,numberOfSpecies) = 1 

    deallocate ( auxcilevel )
    deallocate ( cilevel )
    
  end subroutine ConfigurationInteraction_buildCIOrderList

    !! Search which combinations of excitations satifies the desired CI level.
recursive  function ConfigurationInteraction_buildCIOrderRecursion( s, numberOfSpecies, c, cilevel ) result (os)
    implicit none

    integer :: u,v,c
    integer :: i, j, ii, jj, nn, k, l
    integer :: s, numberOfSpecies
    integer :: os,is,auxis, auxos
    integer :: cilevel(:)
    integer :: plusOne(3,3) , plusTwo(4,6)

    is = s + 1
    if ( is < numberOfSpecies ) then
      do i = 1, size(ConfigurationInteraction_instance%numberOfStrings(is)%values, dim = 1)
       cilevel(is) = i - 1
       os = ConfigurationInteraction_buildCIOrderRecursion( is, numberOfSpecies, c, cilevel )
      end do
      cilevel(is) = 0
    else 
      do i = 1, size(ConfigurationInteraction_instance%numberOfStrings(is)%values, dim = 1)
       cilevel(is) = i - 1
       c = c + 1

       ConfigurationInteraction_instance%ciOrderList( c, : ) = cilevel(:)
       if ( sum(cilevel) <= ConfigurationInteraction_instance%maxCIlevel ) then
         ConfigurationInteraction_instance%sizeCiOrderList = ConfigurationInteraction_instance%sizeCiOrderList + 1
         ConfigurationInteraction_instance%auxciOrderList(  ConfigurationInteraction_instance%sizeCiOrderList  ) = c
       end if

       if ( trim(ConfigurationInteraction_instance%level) == "CISD+" ) then !!special case. 
         plusOne(:,1) = (/1,1,1/)
         plusOne(:,2) = (/2,0,1/)
         plusOne(:,3) = (/0,2,1/)
       
         do k = 1, 3
           if ( sum(  abs(cilevel(:) - plusOne(:,k)) ) == 0 ) then
           ConfigurationInteraction_instance%sizeCiOrderList = ConfigurationInteraction_instance%sizeCiOrderList + 1
           ConfigurationInteraction_instance%auxciOrderList(  ConfigurationInteraction_instance%sizeCiOrderList  ) = c
           end if
         end do
       
       end if
       
       if ( trim(ConfigurationInteraction_instance%level) == "CISD+2" ) then !!special case. 
         plusTwo(:,1) = (/1,1,1,0/)
         plusTwo(:,2) = (/1,1,0,1/)
         plusTwo(:,3) = (/2,0,1,0/)
         plusTwo(:,4) = (/2,0,0,1/)
         plusTwo(:,5) = (/0,2,1,0/)
         plusTwo(:,6) = (/0,2,0,1/)

         do k = 1, 6
           if ( sum(  abs(cilevel(:) - plusTwo(:,k)) ) == 0 ) then
           ConfigurationInteraction_instance%sizeCiOrderList = ConfigurationInteraction_instance%sizeCiOrderList + 1
           ConfigurationInteraction_instance%auxciOrderList(  ConfigurationInteraction_instance%sizeCiOrderList  ) = c
           end if
         end do
 
       end if

      end do
      cilevel(is) = 0
    end if

  end function ConfigurationInteraction_buildCIOrderRecursion

!! Build a list with all possible combinations of number of different orbitals from all quantum species, coupling (0,1,2)

  subroutine ConfigurationInteraction_buildCouplingOrderList()
    implicit none

    integer(8) :: a,b,c,c1,c2,aa,d
    integer :: u,uu,vv, p, nn,z
    integer :: i
    integer :: numberOfSpecies, auxnumberOfSpecies,s
    integer(1), allocatable :: couplingOrder(:)
    integer(1) :: coupling
    real(8) :: timeA, timeB
    integer :: ncouplingOrderOne
    integer :: ncouplingOrderTwo
    integer :: ssize
    integer, allocatable :: cilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    ssize = 1
    do i = 1, numberOfSpecies
      ssize = ssize * 3 !! ( 0,1,2) different orbitals
    end do

    allocate ( ConfigurationInteraction_instance%couplingOrderList( 3, ssize ) ) !! one, two same, two diff 
    allocate ( ConfigurationInteraction_instance%couplingOrderIndex( 3, ssize ) ) !! one, two same, two diff 

    do a = 1, 3
      do b = 1, ssize
        call Vector_constructorInteger1( ConfigurationInteraction_instance%couplingOrderList(a,b), &
          int( numberOfSpecies,8), int(0,1) )

      end do
    end do

    !! same species
    do b = 1, ssize
      call Vector_constructorInteger1( ConfigurationInteraction_instance%couplingOrderIndex(1,b), 1_8, int(0,1) )
      call Vector_constructorInteger1( ConfigurationInteraction_instance%couplingOrderIndex(2,b), 1_8, int(0,1) )
    end do

    !! diff species
    do b = 1, ssize
      call Vector_constructorInteger1( ConfigurationInteraction_instance%couplingOrderIndex(3,b), 2_8, int(0,1) )
    end do


    allocate ( couplingOrder ( numberOfSpecies )) !! 0, 1, 2
    couplingOrder = 0

    !! call recursion
    s = 0
    ConfigurationInteraction_instance%ncouplingOrderOne = 0
    ConfigurationInteraction_instance%ncouplingOrderTwo = 0
    ConfigurationInteraction_instance%ncouplingOrderTwoDiff = 0

    allocate ( ciLevel ( numberOfSpecies ) )
    ciLevel = 0

    !! get all combinations
    auxnumberOfSpecies = ConfigurationInteraction_buildCouplingOrderRecursion( s, numberOfSpecies, couplingOrder, cilevel )

    !! save the index for species (speciesID) just to avoid a lot of conditionals later!

    do u = 1, ConfigurationInteraction_instance%ncouplingOrderOne
      do i = 1, numberOfSpecies
        if ( ConfigurationInteraction_instance%couplingOrderList(1,u)%values(i) == 1 ) then
          ConfigurationInteraction_instance%couplingOrderIndex(1,u)%values(1) = i
        end if
      end do
    end do

    do u = 1, ConfigurationInteraction_instance%ncouplingOrderTwo
      do i = 1, numberOfSpecies
        if ( ConfigurationInteraction_instance%couplingOrderList(2,u)%values(i) == 2 ) then
          ConfigurationInteraction_instance%couplingOrderIndex(2,u)%values(1) = i
        end if
      end do
    end do

    do u = 1, ConfigurationInteraction_instance%ncouplingOrderTwoDiff
      z = 0 
      do i = 1, numberOfSpecies
        if ( ConfigurationInteraction_instance%couplingOrderList(3,u)%values(i) == 1 ) then
          z = z + 1
          ConfigurationInteraction_instance%couplingOrderIndex(3,u)%values(z) = i
        end if
      end do
    end do


    deallocate ( ciLevel )
    deallocate ( couplingOrder ) 

  end subroutine ConfigurationInteraction_buildCouplingOrderList

!! Get all possible combinations of number of different orbitals from all quantum species.
recursive  function ConfigurationInteraction_buildCouplingOrderRecursion( s, numberOfSpecies, couplingOrder, cilevel ) result (os)
    implicit none

    integer(8) :: a,b,c,d
    integer :: u,v
    integer :: i, j, ii, jj, nn
    integer :: s, numberOfSpecies
    integer :: os,is,auxis, auxos
    integer(1) :: couplingOrder(:)
    logical :: same
    integer :: cilevel(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      if ( sum ( couplingOrder) <= 2 ) then
        do i = 1, 3 - sum ( couplingOrder ) !! 0,1,2
          couplingOrder(is) = i-1
          couplingOrder(is+1:) = 0
          os = ConfigurationInteraction_buildCouplingOrderRecursion( is, numberOfSpecies, couplingOrder, cilevel )
        end do
      end if
    else 
      if ( sum ( couplingOrder) <= 2 ) then
        do i = 1, 3 - sum ( couplingOrder ) !! 0,1,2
          couplingOrder(is) = i-1
          couplingOrder(is+1:) = 0
          os = is
          if ( sum ( couplingOrder ) == 1 ) then

            auxis = 0
            ConfigurationInteraction_instance%ncouplingOrderOne = ConfigurationInteraction_instance%ncouplingOrderOne + 1
            b = ConfigurationInteraction_instance%ncouplingOrderOne
            ConfigurationInteraction_instance%couplingOrderList(1,b)%values = couplingOrder

          else if ( sum ( couplingOrder ) == 2 ) then

            same = .false. 

            do j = 1, numberOfSpecies
              if ( couplingOrder(j) == 2 ) same = .true.
            end do

            if ( same ) then
              auxis = 0
              ConfigurationInteraction_instance%ncouplingOrderTwo = ConfigurationInteraction_instance%ncouplingOrderTwo + 1
              b = ConfigurationInteraction_instance%ncouplingOrderTwo
              ConfigurationInteraction_instance%couplingOrderList(2,b)%values = couplingOrder
            else 
              auxis = 0
              ConfigurationInteraction_instance%ncouplingOrderTwoDiff = ConfigurationInteraction_instance%ncouplingOrderTwoDiff + 1
              b = ConfigurationInteraction_instance%ncouplingOrderTwoDiff
              ConfigurationInteraction_instance%couplingOrderList(3,b)%values = couplingOrder
            end if

          end if
        end do
      end if
    end if

  end function ConfigurationInteraction_buildCouplingOrderRecursion


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine ConfigurationInteraction_destructor()
    implicit none
    integer i,j,m,n,p,q,c
    integer numberOfSpecies
    integer :: isLambdaEqual1

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !!Destroy configurations
    !!Ground State
    if (allocated(ConfigurationInteraction_instance%configurations)) then
      c=1
      call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )
  
      do c=2, ConfigurationInteraction_instance%numberOfConfigurations
         call Configuration_destructor(ConfigurationInteraction_instance%configurations(c) )                
      end do
  
      if (allocated(ConfigurationInteraction_instance%configurations)) deallocate(ConfigurationInteraction_instance%configurations)
    end if

    call Matrix_destructor(ConfigurationInteraction_instance%hamiltonianMatrix)
    call Vector_destructorInteger (ConfigurationInteraction_instance%numberOfOccupiedOrbitals)
    call Vector_destructorInteger (ConfigurationInteraction_instance%numberOfOrbitals)
    call Vector_destructor (ConfigurationInteraction_instance%lambda)

    ConfigurationInteraction_instance%isInstanced=.false.

  end subroutine ConfigurationInteraction_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_show()
    implicit none
    type(ConfigurationInteraction) :: this
    integer :: i
    real(8) :: davidsonCorrection, HFcoefficient, CIcorrection
    integer numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( ConfigurationInteraction_instance%isInstanced ) then

       write(*,"(A)") ""
       write(*,"(A)") " POST HARTREE-FOCK CALCULATION"
       write(*,"(A)") " CONFIGURATION INTERACTION THEORY:"
       write(*,"(A)") "=============================="
       write(*,"(A)") ""
       write (6,"(T8,A30, A5)") "LEVEL = ", ConfigurationInteraction_instance%level
       write (6,"(T8,A30, I8)") "NUMBER OF CONFIGURATIONS = ", ConfigurationInteraction_instance%numberOfConfigurations
       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
        write (6,"(T8,A17,I3,A10, F18.12)") "STATE: ", i, " ENERGY = ", ConfigurationInteraction_instance%eigenvalues%values(i)
       end do
       write(*,"(A)") ""
       CIcorrection = ConfigurationInteraction_instance%eigenvalues%values(1) - &
                HartreeFock_instance%totalEnergy

       write (6,"(T4,A34, F20.12)") "GROUND STATE CORRELATION ENERGY = ", CIcorrection

       if (  ConfigurationInteraction_instance%level == "CISD" ) then
         write(*,"(A)") ""
         write (6,"(T2,A34)") "RENORMALIZED DAVIDSON CORRECTION:"
         write(*,"(A)") ""
         write (6,"(T8,A54)") "E(CISDTQ) \approx E(CISD) + \delta E(Q)               "
         write (6,"(T8,A54)") "\delta E(Q) = (1 - c_0^2) * \delta E(CISD) / c_0^2    "
         write (*,*) ""
         HFcoefficient = ConfigurationInteraction_instance%eigenVectors%values(1,1) 
         davidsonCorrection = ( 1 - HFcoefficient*HFcoefficient) * CIcorrection / (HFcoefficient*HFcoefficient)
  
  
         write (6,"(T8,A19, F20.12)") "HF COEFFICIENT = ", HFcoefficient
         write (6,"(T8,A19, F20.12)") "\delta E(Q) = ", davidsonCorrection
         write (6,"(T8,A19, F20.12)") "E(CISDTQ) ESTIMATE ",  HartreeFock_instance%totalEnergy +&
            CIcorrection + davidsonCorrection
       else 

         write(*,"(A)") ""
         HFcoefficient = ConfigurationInteraction_instance%eigenVectors%values(1,1) 
         write (6,"(T8,A19, F20.12)") "HF COEFFICIENT = ", HFcoefficient

       end if

    else 

    end if

  end subroutine ConfigurationInteraction_show

  subroutine ConfigurationInteraction_showEigenVectors()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v,p
    integer :: ci
    integer :: i, j, ii, jj
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:), auxcilevel(:), dd(:)


    if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "NONE" ) return

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    numberOfConfigurations = ConfigurationInteraction_instance%numberOfConfigurations 

    allocate ( ConfigurationInteraction_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )
    ciLevel = 0
    ConfigurationInteraction_instance%allIndexConf = 0
    indexConf = 0

    !! gather all configurations
    s = 0
    c = 0
    ciLevel = 0

    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = ConfigurationInteraction_gatherConfRecursion( s, numberOfSpecies, indexConf,  c, cilevel )
    end do
    !stop

    deallocate ( ciLevel )

    if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "ORBITALS" ) then
    write (*,*) ""
    write (*, "(T1,A)") "Eigenvectors" 
    write (*,*) ""

    do c = 1, CONTROL_instance%NUMBER_OF_CI_STATES
      write (*, "(T1,A,I4,A,F18.10)") "State: ", c, " Energy: ", ConfigurationInteraction_instance%eigenValues%values(c) 
      write (*, "(T1,A)") "Conf, orbital occupation per species, coefficient"
      write (*,*) ""
      do a = 1, numberOfConfigurations
        if ( abs(ConfigurationInteraction_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
          indexConf(:) = ConfigurationInteraction_instance%allIndexConf(:,a) 

          write (*, "(T1,I8,A1)", advance="no") a, " "
          do i = 1, numberOfSpecies
            do p = 1, ConfigurationInteraction_instance%numberOfOrbitals%values(i)
              write (*, "(I1)", advance="no")  ConfigurationInteraction_instance%orbitals(i)%values(p,indexConf(i)) 
            end do
            write (*, "(A1)", advance="no")  " "
          end do
          write (*, "(F11.8)") ConfigurationInteraction_instance%eigenVectors%values(a,c) 
        end if
      end do
      write (*,*) ""
    end do


    else if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "OCCUPIED" ) then
    write (*,*) ""
    write (*, "(T1,A)") "Eigenvectors" 
    write (*,*) ""

    do c = 1, CONTROL_instance%NUMBER_OF_CI_STATES
      write (*, "(T1,A,I4,A,F18.10)") "State: ", c, " Energy: ", ConfigurationInteraction_instance%eigenValues%values(c) 
      write (*, "(T1,A)") "Conf, occupied orbitals per species, coefficient"
      write (*,*) ""
      do a = 1, numberOfConfigurations
        if ( abs(ConfigurationInteraction_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
          indexConf(:) = ConfigurationInteraction_instance%allIndexConf(:,a) 

          write (*, "(T1,I8,A1)", advance="no") a, " "
          do i = 1, numberOfSpecies
            do p = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
              write (*, "(I3,A1)", advance="no") ConfigurationInteraction_instance%strings(i)%values(p,indexConf(i) ), " "
            end do
            write (*, "(A1)", advance="no")  "|"
          end do
          write (*, "(A,F11.8)") " ", ConfigurationInteraction_instance%eigenVectors%values(a,c) 
        end if
      end do
      write (*,*) ""
    end do

    end if

    deallocate ( indexConf )
    deallocate ( ConfigurationInteraction_instance%allIndexConf )

  end subroutine ConfigurationInteraction_showEigenVectors


  !FELIX IS HERE
  subroutine ConfigurationInteraction_densityMatrices()
    implicit none
    type(ConfigurationInteraction) :: this
    type(Configuration) :: auxthisA, auxthisB
    integer :: i, j, k, l, mu, nu, n
    integer :: factor
    integer :: unit, wfnunit
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state, species, orbital, orbitalA, orbitalB
    character(50) :: file, wfnfile, speciesName, auxstring
    character(50) :: arguments(2)
    type(matrix), allocatable :: coefficients(:), atomicDensityMatrix(:,:), ciDensityMatrix(:,:), auxDensMatrix(:,:)
    integer numberOfSpecies

    type(matrix) :: auxdensityEigenVectors 
    type(matrix) :: densityEigenVectors
    type(vector) :: auxdensityEigenValues
    type(vector) :: densityEigenValues
    integer, allocatable :: cilevel(:), cilevelA(:)
    integer(8) :: numberOfConfigurations, c
    integer(8), allocatable :: indexConf(:)
    type(ivector), allocatable :: stringAinB(:)
    integer :: s, ss, ci, auxnumberOfSpecies
    integer, allocatable :: coupling(:)
    integer :: a, b, AA, BB, bj
    integer :: u, uu, ssize
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(8), allocatable :: jj(:)
    real(8) :: timeDA
    real(8) :: timeDB

  
    ! type(Vector) :: eigenValues
    ! type(Matrix) :: eigenVectors, auxMatrix
    ! real(8) :: sumaPrueba

    !!Iterators: i,j - Configurations .... k,l - molecular orbitals .... mu,nu - atomic orbitals ... n - threads
    if ( ConfigurationInteraction_instance%isInstanced .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0 ) then
       !$  timeDA = omp_get_wtime()

      numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
  
      numberOfConfigurations = ConfigurationInteraction_instance%numberOfConfigurations 
  
      allocate (stringAinB ( numberOfSpecies ))
  
      do i = 1, numberOfSpecies 
        call Vector_constructorInteger (stringAinB(i), ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i), 0)
      end do 
  
      allocate ( ConfigurationInteraction_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
      allocate ( ciLevelA ( numberOfSpecies ) )
      allocate ( ciLevel ( numberOfSpecies ) )
      allocate ( indexConf ( numberOfSpecies ) )
      ciLevelA = 0
      ciLevel = 0
      ConfigurationInteraction_instance%allIndexConf = 0
      indexConf = 0
  
      !! gather all configurations
      s = 0
      c = 0
      ciLevel = 0
  
      do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
  
        cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
        s = 0
        auxnumberOfSpecies = ConfigurationInteraction_gatherConfRecursion( s, numberOfSpecies, indexConf,  c, cilevel )
      end do
      !stop
  
      deallocate ( indexConf )
      allocate ( coupling ( numberOfSpecies ) )


      write (*,*) ""
      write (*,*) "=============================="
      write (*,*) "BUILDING CI DENSITY MATRICES"
      write (*,*) "=============================="
      write (*,*) ""

      allocate( coefficients(numberOfSpecies), atomicDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), &
               ciDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), auxDensMatrix(numberOfSpecies,ConfigurationInteraction_instance%nproc) )

      wfnFile = "lowdin.wfn"
      wfnUnit = 20
      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      !Inicializando las matrices
      do species=1, numberOfSpecies
         speciesName = MolecularSystem_getNameOfSpecie(species)
         
         numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
         ! numberOfOrbitals = ConfigurationInteraction_instance%numberOfOrbitals%values(species)
         numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(species)

        arguments(2) = speciesName
        arguments(1) = "COEFFICIENTS"
        ! print *, "trolo", numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals

        coefficients(species) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        ! print *, "trololo"
        
        do state=1, CONTROL_instance%CI_STATES_TO_PRINT

           call Matrix_constructor ( ciDensityMatrix(species,state) , &
                int(numberOfContractions,8), &
                int(numberOfContractions,8),  0.0_8 )

           do k=1, numberOfOccupiedOrbitals
             ciDensityMatrix(species,state)%values( k, k)=1.0_8
           end do

        end do

        do n=1, ConfigurationInteraction_instance%nproc

           call Matrix_constructor ( auxDensMatrix(species,n) , &
                int(numberOfContractions,8), &
                int(numberOfContractions,8),  0.0_8 )
        end do
      end do
       
      close(wfnUnit)

      allocate ( indexConfA ( numberOfSpecies ) )
      allocate ( indexConfB ( numberOfSpecies ) )
      allocate ( jj ( numberOfSpecies ) )

      indexConfA = 0
      indexConfB = 0
      jj = 0

      !! Building the CI reduced density matrix in the molecular orbital representation in parallel
      ! call Matrix_show (ConfigurationInteraction_instance%eigenVectors)

      !!print *, "        State, Progress"
      
      do state=1, CONTROL_instance%CI_STATES_TO_PRINT

         !$omp parallel & 
         !$omp& firstprivate (stringAinB,indexConfA,indexConfB, jj) &
         !$omp& private(i,j, species, s, numberOfOccupiedOrbitals, k, coupling, orbital, orbitalA, orbitalB, AA, BB, a, b, factor, n, cilevelA, ss, ssize, cilevel, ci, u, uu, bj),&
         !$omp& shared(ConfigurationInteraction_instance, auxDensMatrix )
         n = omp_get_thread_num() + 1
         !$omp do schedule (dynamic) 
         do i=1, ConfigurationInteraction_instance%numberOfConfigurations

            !!if( mod( i , 50000 ) .eq. 0 ) print *, state, floor(real(100*i/ConfigurationInteraction_instance%numberOfConfigurations)), "%"
            !!Filter very small coefficients
            if( abs(ConfigurationInteraction_instance%eigenVectors%values(i,state)) .ge. 1E-10) then

               indexConfA(:) = ConfigurationInteraction_instance%allIndexConf(:,i) 

               !print *, "==", indexConfA , "|", i


               !!Diagonal contributions
               do species=1, numberOfSpecies
                  numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(species)

                  do k=1, numberOfOccupiedOrbitals

                     !!Occupied orbitals
                     auxDensMatrix(species,n)%values(k,k)=auxDensMatrix(species,n)%values(k,k) - ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
                     ! ciDensityMatrix(species,state)%values( k, k) = ciDensityMatrix(species,state)%values( k, k) -  &
                     !      ConfigurationInteraction_instance%eigenVectors%values(i,state)**2

                     !print *, i, j, k, species 
                     !orbital = ConfigurationInteraction_instance%configurations(i)%occupations(k,species) 
                     orbital =  ConfigurationInteraction_instance%strings(species)%values(k,indexConfA(species))
                     !!Unoccupied orbitals

                     auxDensMatrix(species,n)%values(orbital,orbital)=auxDensMatrix(species,n)%values(orbital,orbital) + ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
                     ! ciDensityMatrix(species,state)%values( orbital, orbital)= ciDensityMatrix(species,state)%values( orbital, orbital) + &
                     !      ConfigurationInteraction_instance%eigenVectors%values(i,state)**2

                  end do
               end do

               !!Off Diagonal contributions
               cilevelA = 0
               do ss = 1, numberOfSpecies 
                 stringAinB(ss)%values = 0
                 do k = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(ss)

                   stringAinB(ss)%values(k) = ConfigurationInteraction_instance%orbitals(ss)%values( &
                                             ConfigurationInteraction_instance%strings(ss)%values(k,  ConfigurationInteraction_instance%allIndexConf(ss,1)), indexConfA(ss))
                 end do
                 cilevelA(ss) = configurationinteraction_instance%numberOfOccupiedOrbitals%values(ss) - sum ( stringAinB(ss)%values )
               end do 

               jj = 0
               coupling = 0
               do ss = 1, numberOfSpecies 
                 ssize = 0 

                 indexConfB(:) = indexConfA(:)
                 cilevel = cilevelA

                 do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(ss)%values, dim = 1)
                   cilevel(ss) = ci - 1
                   do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
                     if ( sum(abs(cilevel - &
                          ConfigurationInteraction_instance%ciOrderList( ConfigurationInteraction_instance%auxciOrderList(u), :))) == 0 ) then
                       uu = ConfigurationInteraction_instance%auxciOrderList(u)
                       do bj = 1 + ssize , ConfigurationInteraction_instance%numberOfStrings(ss)%values(ci) + ssize
                         indexConfB(ss) = bj
  
                         do s=1, numberOfSpecies
                           jj(s) = (indexConfB(s) - ConfigurationInteraction_instance%numberOfStrings2(s)%values(cilevel(s)+1) + &
                                    ConfigurationInteraction_instance%ciOrderSize1(uu,s) )* ConfigurationInteraction_instance%ciOrderSize2(uu,s) 
                         end do

                         j = sum(jj)
                         !print *, "  ", indexConfB , "|", j, ConfigurationInteraction_instance%eigenVectors%values(j,state) 
                         if ( j > i ) then
                           if( abs(ConfigurationInteraction_instance%eigenVectors%values(j,state)) .ge. 1E-10) then

                             coupling = 0
                             do s=1, numberOfSpecies
                                stringAinB(s)%values = 0
                                do k = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s)
                                   stringAinB(s)%values(k) = ConfigurationInteraction_instance%orbitals(s)%values( &
                                        ConfigurationInteraction_instance%strings(s)%values(k,indexConfA(s) ), indexConfB(s) ) 
                                end do
                                coupling(s) = configurationinteraction_instance%numberOfOccupiedOrbitals%values(s) - sum ( stringAinB(s)%values )
                             end do
                             if (sum(coupling) == 1) then
    
                               do s = 1, numberOfSpecies
    
                                 if ( coupling(s) == 1) then !!hmm

                                   !print *, "      ", coupling
                                   orbitalA = 0
                                   orbitalB = 0
                                   AA = 0
                                   BB = 0
                                   a = indexConfA(s)
                                   b = indexConfB(s)
    
                                   do k = 1, ConfigurationInteraction_instance%occupationNumber(s) 
                                      if ( ConfigurationInteraction_instance%orbitals(s)%values( &
                                           ConfigurationInteraction_instance%strings(s)%values(k,a),b) == 0 ) then
                                         orbitalA =  ConfigurationInteraction_instance%strings(s)%values(k,a)
                                         AA = k
                                         exit
                                      end if
                                   end do
                                   do k = 1, ConfigurationInteraction_instance%occupationNumber(s) 
                                      if ( ConfigurationInteraction_instance%orbitals(s)%values( &
                                           ConfigurationInteraction_instance%strings(s)%values(k,b),a) == 0 ) then
                                         orbitalB =  ConfigurationInteraction_instance%strings(s)%values(k,b)
                                         BB = k
                                         exit
                                      end if
                                   end do
    
                                   factor = (-1)**(AA-BB)
    
                                   numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s)
    
                                   ! print *, i, j, ConfigurationInteraction_instance%configurations(i)%occupations(:,specie), ConfigurationInteraction_instance%configurations(j)%occupations(:,specie)
                                   ! print *, i, j, auxthisA%occupations(:,specie), auxthisB%occupations(:,specie)
                                   ! print *, i, j, orbitalA, orbitalB, factor*ConfigurationInteraction_instance%eigenVectors%values(i,1)*ConfigurationInteraction_instance%eigenVectors%values(j,1)
    
                                   auxDensMatrix(s,n)%values( orbitalA,orbitalB)= auxDensMatrix(s,n)%values( orbitalA, orbitalB) + &
                                        factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
                                        ConfigurationInteraction_instance%eigenVectors%values(j,state)
                                   auxDensMatrix(s,n)%values( orbitalB,orbitalA)= auxDensMatrix(s,n)%values( orbitalB, orbitalA) + &
                                        factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
                                        ConfigurationInteraction_instance%eigenVectors%values(j,state)
                                  end if
                                end do
                              end if
                           end if
                         end if
                         !! here
                       end do
                       ssize = ssize + ConfigurationInteraction_instance%numberOfStrings(ss)%values(ci)
                       !exit
                     end if

                   end do
                 end do

               end do 

!               do j=i+1, ConfigurationInteraction_instance%numberOfConfigurations
!                  if( abs(ConfigurationInteraction_instance%eigenVectors%values(j,state)) .ge. 1E-12) then

!                     indexConfB(:) = ConfigurationInteraction_instance%allIndexConf(:,j)

!                     coupling = 0
!                     do s=1, numberOfSpecies
!                        stringAinB(s)%values = 0
!                        do k = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s)
!                           stringAinB(s)%values(k) = ConfigurationInteraction_instance%orbitals(s)%values( &
!                                ConfigurationInteraction_instance%strings(s)%values(k,indexConfA(s) ), indexConfB(s) ) 
!                        end do
!                        coupling(s) = configurationinteraction_instance%numberOfOccupiedOrbitals%values(s) - sum ( stringAinB(s)%values )
!                     end do
!
!                     if (sum(coupling) == 1) then
!
!                        do s = 1, numberOfSpecies
!
!                           if ( coupling(s) == 1) then
!                              orbitalA = 0
!                              orbitalB = 0
!                              AA = 0
!                              BB = 0
!                              a = indexConfA(s)
!                              b = indexConfB(s)
!
!                              do k = 1, ConfigurationInteraction_instance%occupationNumber(s) 
!                                 if ( ConfigurationInteraction_instance%orbitals(s)%values( &
!                                      ConfigurationInteraction_instance%strings(s)%values(k,a),b) == 0 ) then
!                                    orbitalA =  ConfigurationInteraction_instance%strings(s)%values(k,a)
!                                    AA = k
!                                    exit
!                                 end if
!                              end do
!                              do k = 1, ConfigurationInteraction_instance%occupationNumber(s) 
!                                 if ( ConfigurationInteraction_instance%orbitals(s)%values( &
!                                      ConfigurationInteraction_instance%strings(s)%values(k,b),a) == 0 ) then
!                                    orbitalB =  ConfigurationInteraction_instance%strings(s)%values(k,b)
!                                    BB = k
!                                    exit
!                                 end if
!                              end do
!
!                              factor = (-1)**(AA-BB)
!
!                              numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s)
!
!                              ! print *, i, j, ConfigurationInteraction_instance%configurations(i)%occupations(:,specie), ConfigurationInteraction_instance%configurations(j)%occupations(:,specie)
!                              ! print *, i, j, auxthisA%occupations(:,specie), auxthisB%occupations(:,specie)
!
!                              ! print *, i, j, orbitalA, orbitalB, factor*ConfigurationInteraction_instance%eigenVectors%values(i,1)*ConfigurationInteraction_instance%eigenVectors%values(j,1)
!
!                              auxDensMatrix(s,n)%values( orbitalA,orbitalB)= auxDensMatrix(s,n)%values( orbitalA, orbitalB) + &
!                                   factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
!                                   ConfigurationInteraction_instance%eigenVectors%values(j,state)
!                              ! ciDensityMatrix(s,state)%values( orbitalA,orbitalB)= ciDensityMatrix(s,state)%values( orbitalA, orbitalB) + &
!                              !      factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
!                              !      ConfigurationInteraction_instance%eigenVectors%values(j,state)
!
!                              auxDensMatrix(s,n)%values( orbitalB,orbitalA)= auxDensMatrix(s,n)%values( orbitalB, orbitalA) + &
!                                   factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
!                                   ConfigurationInteraction_instance%eigenVectors%values(j,state)
!
!                              ! ciDensityMatrix(s,state)%values( orbitalB, orbitalA)= ciDensityMatrix(s,state)%values( orbitalB, orbitalA) + &
!                              !      factor*ConfigurationInteraction_instance%eigenVectors%values(i,state)* &
!                              !      ConfigurationInteraction_instance%eigenVectors%values(j,state)
!
!                           end if
!                        end do
!                     end if
!                  end if
!               end do

            end if
         end do
         !$omp end do nowait
         !$omp end parallel
         
         !! Gather the parallel results
         do species=1, numberOfSpecies
            do n=1, ConfigurationInteraction_instance%nproc
               ciDensityMatrix(species,state)%values = ciDensityMatrix(species,state)%values + auxDensMatrix(species,n)%values
            end do
         end do
      end do


      !! Open file - to write density matrices
     unit = 29
       
     file = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
     open(unit = unit, file=trim(file), status="new", form="formatted")
       
     !! Building the CI reduced density matrix in the atomic orbital representation       
     do species=1, numberOfSpecies
       speciesName = MolecularSystem_getNameOfSpecie(species)
       ! numberOfOrbitals = ConfigurationInteraction_instance%numberOfOrbitals%values(species)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
       numberOfOccupiedOrbitals = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(species)

       do state=1, CONTROL_instance%CI_STATES_TO_PRINT
          
         ! print *, "CI density matrix ", trim(speciesName), state
         ! call Matrix_show ( ciDensityMatrix(species,state))
             
         call Matrix_constructor ( atomicDensityMatrix(species,state) , &
                                   int(numberOfContractions,8), &
                                   int(numberOfContractions,8),  0.0_8 )

         do mu=1, numberOfContractions
           do nu=1, numberOfContractions
             do k=1, numberOfContractions
               do l=k, numberOfContractions

                 if(l .eq. k) then
                   atomicDensityMatrix(species,state)%values(mu,nu) =  &
                     atomicDensityMatrix(species,state)%values(mu,nu) + &
                     ciDensityMatrix(species,state)%values(k,k) *&
                     coefficients(species)%values(mu,k)*coefficients(species)%values(nu,k)

                 else
                   atomicDensityMatrix(species,state)%values(mu,nu) =  &
                     atomicDensityMatrix(species,state)%values(mu,nu) + &
                     ciDensityMatrix(species,state)%values(k,l) *&
                     (coefficients(species)%values(mu,k)*coefficients(species)%values(nu,l) + & 
                     coefficients(species)%values(mu,l)*coefficients(species)%values(nu,k))
                            
                 end if
               end do
             end do
           end do
         end do
             
         !!print *, "atomic density matrix  ", trim(speciesName), state
         !!call Matrix_show ( atomicDensityMatrix(species,state))

         write(auxstring,*) state
         arguments(2) = speciesName
         arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
             
         call Matrix_writeToFile ( atomicDensityMatrix(species,state), unit , arguments=arguments(1:2) )

         end do
       end do


      !! Natural orbitals

      write(*,*) ""
      write(*,*) "=============================="
      write(*,*) " NATURAL ORBITALS: "
      write(*,*) ""

      do state=1, CONTROL_instance%CI_STATES_TO_PRINT

        write(*,*) " STATE: ", state

        do species=1, numberOfSpecies

          write(*,*) ""
          write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(species)%name )
          write(*,*) "-----------------"

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
          speciesName = MolecularSystem_getNameOfSpecie(species)


          call Vector_constructor ( auxdensityEigenValues, &
                                   int(numberOfContractions,4),  0.0_8 )

          call Matrix_constructor ( auxdensityEigenVectors, &
                                   int(numberOfContractions,8), &
                                   int(numberOfContractions,8),  0.0_8 )

          call Vector_constructor ( densityEigenValues, &
                                   int(numberOfContractions,4),  0.0_8 )

          call Matrix_constructor ( densityEigenVectors, &
                                   int(numberOfContractions,8), &
                                   int(numberOfContractions,8),  0.0_8 )

          call Matrix_eigen ( ciDensityMatrix(species,state), auxdensityEigenValues, auxdensityEigenVectors, SYMMETRIC )  

          !! reorder
          do u = 1, numberOfContractions
            densityEigenValues%values(u) =  auxdensityEigenValues%values(numberOfContractions - u + 1)
            densityEigenVectors%values(:,u) = auxdensityEigenVectors%values(:,numberOfContractions - u + 1)
          end do


          !! Transform to atomic basis
          densityEigenVectors%values = matmul( coefficients(species)%values, densityEigenVectors%values )

          call Matrix_show( densityEigenVectors, &
             rowkeys = MolecularSystem_getlabelsofcontractions( species ), &
             columnkeys = string_convertvectorofrealstostring( densityEigenValues ),&
             flags=WITH_BOTH_KEYS)

          write(auxstring,*) state
          arguments(2) = speciesName
          arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 
             
          call Matrix_writeToFile ( atomicDensityMatrix(species,state), unit , arguments=arguments(1:2) )

         !! it's the same
         !!auxdensityEigenVectors%values = 0

         !!do mu=1, numberOfContractions
         !!  do nu=1, numberOfContractions
         !!    do k=1, numberOfContractions
         !!      auxdensityEigenVectors%values(mu,nu) = auxdensityEigenVectors%values(mu,nu) + &
         !!                              densityEigenVectors%values(mu,k) *  densityEigenVectors%values(nu,k)*densityEigenValues%values(k) 
         !!    end do
         !!  end do
         !!end do
         !!print *, "atomic density matrix from natural orbitals"
         !!call Matrix_show ( auxdensityEigenVectors)

        write(*,*) " End of natural orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(species)%name )
        end do
      end do



      write(*,*) ""
      write(*,*) " END OF NATURAL ORBITALS"
      write(*,*) "=============================="
      write(*,*) ""

      close(unit)

      deallocate ( jj )
      deallocate ( indexConfB )
      deallocate ( indexConfA )
      deallocate ( coupling )
      deallocate ( cilevel )
      deallocate ( cilevelA )
      deallocate ( ConfigurationInteraction_instance%allIndexConf )
      deallocate ( stringAinB )

     deallocate( coefficients, atomicDensityMatrix, ciDensityMatrix )

     !$  timeDB = omp_get_wtime()
     !$  write(*,"(A,F10.4,A4)") "** TOTAL Elapsed Time for Building density matrices: ", timeDB - timeDA ," (s)"

     
  end if
       
    ! print *, i, i, orbital, orbital, ConfigurationInteraction_instance%eigenVectors%values(i,1)**2
    
    ! do mu = 1 , numberOfOrbitals
    !    do nu = 1 , numberOfOrbitals
    
    !       densityMatrix%values(mu,nu) =  &
    !            densityMatrix%values(mu,nu) + &
    !            ConfigurationInteraction_instance%eigenVectors%values(i,state)**2 *&
    !            coefficients%values(mu,orbital)*coefficients%values(nu,orbital)
    !    end do
    ! end do
    
    !!off-Diagonal ground state
    
    ! do mu = 1 , numberOfOrbitals
    !    do nu = 1 , numberOfOrbitals
    
    !       densityMatrix%values(mu,nu) =  &
    !            densityMatrix%values(mu,nu) + &
    !            factor *&
    !            ConfigurationInteraction_instance%eigenVectors%values(i,state) *&
    !            ConfigurationInteraction_instance%eigenVectors%values(j,state) *&
    !            (coefficients%values(mu,orbitalA)*coefficients%values(nu,orbitalB) + coefficients%values(mu,orbitalB)*coefficients%values(nu,orbitalA))
    !    end do
    ! end do
    
    ! call Vector_constructor(eigenValues, numberOfOrbitals)
    ! call Matrix_constructor(eigenVectors, int(numberOfOrbitals,8), int(numberOfOrbitals,8))
    ! call Matrix_eigen(ciOccupationMatrix, eigenValues, eigenVectors, SYMMETRIC)
    
    ! print *, "Diagonal sum", sum(eigenValues%values)
    ! call Vector_show(eigenValues)
    
    ! call Matrix_show(eigenVectors)
    ! print *, arguments(1:2)
    ! call Matrix_show ( densityMatrix )
    
    ! call Matrix_constructor ( ciOccupationNumbers , int(numberOfOrbitals,8) , &
    !      int(CONTROL_instance%CI_STATES_TO_PRINT,8),  0.0_8 )
    
    ! do state=1, CONTROL_instance%CI_STATES_TO_PRINT
    !    sumaPrueba=0
    !    do j=1, numberOfOccupiedOrbitals
    !       ciOccupationNumbers%values(j,state) = 1.0
    !    end do
    
    ! ! !Get occupation numbers from each configuration contribution
    
    !    do i=1, ConfigurationInteraction_instance%numberOfConfigurations
    !       do j=1, numberOfOccupiedOrbitals
    
    !          !! Occupied orbitals
    !          ciOccupationNumbers%values( j, state)= ciOccupationNumbers%values( j, state) -  &
    !               ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
    !          !! Unoccupied orbitals
    !          orbital = ConfigurationInteraction_instance%configurations(i)%occupations(j,specie) 
    
    !          ciOccupationNumbers%values( orbital, state)= ciOccupationNumbers%values( orbital, state) + &
    !               ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
    
    !          ! print *, j, orbital, ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
    !          ! sumaPrueba=sumaPrueba+ConfigurationInteraction_instance%eigenVectors%values(i,state)**2
    !       end do
    !       ! end if
    
    !    end do
    
    !    ! print *, "suma", sumaPrueba
    !    !Build a new density matrix (P) in atomic orbitals
    
    !    call Matrix_constructor ( densityMatrix , &
    !         int(numberOfOrbitals,8), &
    !         int(numberOfOrbitals,8),  0.0_8 )
    
    !    wfnFile = "lowdin.wfn"
    !    wfnUnit = 20
    
    !    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    
    !    arguments(2) = speciesName
    !    arguments(1) = "COEFFICIENTS"
    
    !    coefficients = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
    !         columns= int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))
    
    !    close(wfnUnit)
    
    !    do mu = 1 , numberOfOrbitals
    !       do nu = 1 , numberOfOrbitals
    !          do k = 1 , numberOfOrbitals
    
    !             densityMatrix%values(mu,nu) =  &
    !                  densityMatrix%values(mu,nu) + &
    !                  ciOccupationNumbers%values(k, state)**2* &
    !                  coefficients%values(mu,k)*coefficients%values(nu,k)
    !           end do
    !        end do
    !     end do
    
    !     write(auxstring,*) state
    !     arguments(2) = speciesName
    !     arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
    
    !     call Matrix_writeToFile ( densityMatrix, unit , arguments=arguments(1:2) )
    
    !     print *, arguments(1:2)
    !     call Matrix_show ( densityMatrix )
    
    !     call Matrix_destructor(coefficients)          
    !     call Matrix_destructor(densityMatrix)          
    
    
    !  end do
    
    ! !Write occupation numbers to file
    ! write (6,"(T8,A10,A20)") trim(MolecularSystem_getNameOfSpecie(specie)),"OCCUPATIONS:"
    
    ! call Matrix_show ( ciOccupationNumbers )
    
    ! arguments(2) = speciesName
    ! arguments(1) = "OCCUPATIONS"
    
    ! call Matrix_writeToFile ( ciOccupationNumbers, unit , arguments=arguments(1:2) )
    
    ! call Matrix_destructor(ciOccupationNumbers)          
    


  end subroutine ConfigurationInteraction_densityMatrices

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_run()
    implicit none 
    integer :: i,j,m, numberOfSpecies
    real(8), allocatable :: eigenValues(:) 

!    select case ( trim(ConfigurationInteraction_instance%level) )
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    write (*,*) ""
    write (*,*) "==============================================="
    write (*,*) "         BEGIN ", trim(ConfigurationInteraction_instance%level)," CALCULATION"
    write (*,*) "         J. Charry, F. Moncada                "
    write (*,*) "-----------------------------------------------"
    write (*,*) ""

    write (*,"(A32)",advance="no") "Number of orbitals for species: "
    do i = 1, numberOfSpecies-1
      write (*,"(A)",advance="no") trim(MolecularSystem_getNameOfSpecie(i))//", "
    end do
    write (*,"(A)",advance="no") trim(MolecularSystem_getNameOfSpecie(numberOfSpecies))
    write (*,*) ""

    write (*,"(A28)",advance="no") "  occupied orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)", advance="no") ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  virtual orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") int(MolecularSystem_getTotalNumberOfContractions( i )* &
                                                ConfigurationInteraction_instance%lambda%values(i)  - &
                                                ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) )
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  total number of orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") int(MolecularSystem_getTotalNumberOfContractions( i )* &
                       ConfigurationInteraction_instance%lambda%values(i)   )
    end do
    write (*,*) ""


    write (*,"(A28)",advance="no") "  frozen core orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  active occupied orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) - &
                         ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  active virtual orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") ConfigurationInteraction_instance%numberOfOrbitals%values(i) - &
                         ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") " total active orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no")  ConfigurationInteraction_instance%numberOfOrbitals%values(i) - &
                          ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""
    write (*,*) " "

    write (*,*) "Getting transformed integrals..."
    call ConfigurationInteraction_getTransformedIntegrals()
    write (*,*) " "

    !write (*,*) ConfigurationInteraction_instance%fourCenterIntegrals(1,1)%values(171, 1) a bug...
    write (*,*) "Setting CI level..."

    call ConfigurationInteraction_settingCILevel()

   !! write (*,*) "Total number of configurations", ConfigurationInteraction_instance%numberOfConfigurations
    write (*,*) ""
    call Vector_constructor8 ( ConfigurationInteraction_instance%eigenvalues, &
                              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8 )

    select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

    case ("ARPACK")

      write (*,*) "This method was removed"

    case ("JADAMILU")

      write (*,*) "Building Strings..."
      call ConfigurationInteraction_buildStrings()

      write (*,*) "Building CI level table..."
      call ConfigurationInteraction_buildCIOrderList()

      call ConfigurationInteraction_buildCouplingMatrix()
      call ConfigurationInteraction_buildCouplingOrderList()

      write (*,*) "Building diagonal..."
      call ConfigurationInteraction_buildDiagonal()

      write (*,*) "Building initial hamiltonian..."
      call ConfigurationInteraction_buildInitialCIMatrix2()
      !!call ConfigurationInteraction_buildHamiltonianMatrix() This should be  modified to build the CI matrix in memory

      call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
           int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
        call ConfigurationInteraction_loadEigenVector (ConfigurationInteraction_instance%eigenvalues, &
               ConfigurationInteraction_instance%eigenVectors) 
      end if 

      if ( CONTROL_instance%CI_BUILD_FULL_MATRIX ) then
        write (*,*) "Building and saving hamiltonian..."
        call ConfigurationInteraction_buildAndSaveCIMatrix()
      end if

      write(*,*) ""
      write(*,*) "Diagonalizing hamiltonian..."
      write(*,*) "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
      write(*,*) "============================================================="
      write(*,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
      write(*,*) " a software code for computing selected eigenvalues of "
      write(*,*) " large sparse symmetric matrices, "
      write(*,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 
      write(*,*) "============================================================="


      call ConfigurationInteraction_jadamiluInterface(ConfigurationInteraction_instance%numberOfConfigurations, &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), &
           ConfigurationInteraction_instance%eigenvalues, &
           ConfigurationInteraction_instance%eigenVectors )

      if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
        call ConfigurationInteraction_saveEigenVector () 
      end if
    case ("DSYEVX")

      write (*,*) "Building Strings..."
      call ConfigurationInteraction_buildStrings()

      write (*,*) "Building CI level table..."
      call ConfigurationInteraction_buildCIOrderList()

      write (*,*) "Building diagonal..."
      call ConfigurationInteraction_buildDiagonal()

      write (*,*) "Building Hamiltonian..."
      call ConfigurationInteraction_buildHamiltonianMatrix()

      call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
           int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      !! deallocate transformed integrals
      deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
      deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)

      write(*,*) ""
      write(*,*) "Diagonalizing hamiltonian..."
      write(*,*) "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

      call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
           int(1), int(CONTROL_instance%NUMBER_OF_CI_STATES), &  
           eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
           flags = int(SYMMETRIC,4))

!      call Matrix_eigen_select (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!           1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!           flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )


    case ("DSYEVR")

      write (*,*) "Building Strings..."
      call ConfigurationInteraction_buildStrings()

      write (*,*) "Building CI level table..."
      call ConfigurationInteraction_buildCIOrderList()

      write (*,*) "Building diagonal..."
      call ConfigurationInteraction_buildDiagonal()

      write (*,*) "Building Hamiltonian..."
      call ConfigurationInteraction_buildHamiltonianMatrix()

      call Matrix_constructor (ConfigurationInteraction_instance%eigenVectors, &
           int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      !! deallocate transformed integrals
      deallocate(ConfigurationInteraction_instance%twoCenterIntegrals)
      deallocate(ConfigurationInteraction_instance%fourCenterIntegrals)

      call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
           1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
           eigenVectors = ConfigurationInteraction_instance%eigenVectors, &
           flags = SYMMETRIC)

!     call Matrix_eigen_dsyevr (ConfigurationInteraction_instance%hamiltonianMatrix, ConfigurationInteraction_instance%eigenvalues, &
!           1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
!           flags = SYMMETRIC, dm = ConfigurationInteraction_instance%numberOfConfigurations )

    case default

      call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Diagonalization method not implemented")


    end select

    write(*,*) ""
    write(*,*) "-----------------------------------------------"
    write(*,*) "          END ", trim(ConfigurationInteraction_instance%level)," CALCULATION"
    write(*,*) "==============================================="
    write(*,*) ""
         
!    case ( "FCI-oneSpecie" )
!
!       print *, ""
!       print *, ""
!       print *, "==============================================="
!       print *, "|  Full CI for one specie calculation          |"
!       print *, "|  Use fci program to perform the calculation  |"
!       print *, "-----------------------------------------------"
!       print *, ""
!       ! call ConfigurationInteraction_getTransformedIntegrals()
!       !call ConfigurationInteraction_printTransformedIntegralsToFile()
!

  end subroutine ConfigurationInteraction_run

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_settingCILevel()
    implicit none

    integer :: numberOfSpecies
    integer :: i,ii,j,k,l,m,n,p,q,a,b,d,r,s
    integer(8) :: c, cc
    integer :: ma,mb,mc,md,me,pa,pb,pc,pd,pe
    integer :: isLambdaEqual1
    type(ivector) :: order
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)
    integer, allocatable :: auxArray(:,:), auxvector(:),auxvectorA(:)
    integer :: lambda, otherlambda

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( allocated( occupiedCode ) ) deallocate( occupiedCode )
    allocate (occupiedCode ( numberOfSpecies ) )
    if ( allocated( unoccupiedCode ) ) deallocate( unoccupiedCode )
    allocate (unoccupiedCode ( numberOfSpecies ) )

    !1 auxiliary string for omp paralelization
    do n = 1, ConfigurationInteraction_instance%nproc
      do i = 1, numberOfSpecies
        call Vector_constructorInteger( ConfigurationInteraction_instance%auxstring(n,i), &
          int(ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i),4), int(0,4))
      end do  
    end do  

    select case ( trim(ConfigurationInteraction_instance%level) )

    case ( "FCI" )

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
      end do

      ConfigurationInteraction_instance%maxCILevel = sum(ConfigurationInteraction_instance%CILevel)

    case ( "CIS" )

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 1
      end do
      ConfigurationInteraction_instance%maxCILevel = 1

    case ( "CISD" )

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 2
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 2

    case ( "CISD+" )

      if ( .not. numberOfSpecies == 3 ) call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "CISD+ is specific for three quantum species")

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 2
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 2

    case ( "CISD+2" )

      if ( .not. numberOfSpecies == 4 ) call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "CISD+2 is specific for three quantum species")
      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 2
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 2

    case ("CISDT")

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 3
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 3 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 3

    case ("CISDTQ")

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 4
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 4 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 4

    case ("CISDTQQ")

      do i=1, numberOfSpecies
        ConfigurationInteraction_instance%CILevel(i) = 5
        if ( ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) < 5 ) &
          ConfigurationInteraction_instance%CILevel(i) = ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      ConfigurationInteraction_instance%maxCILevel = 5

    case default

       call ConfigurationInteraction_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select


  end subroutine ConfigurationInteraction_settingCILevel

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildCouplingMatrix()
    implicit none

    integer(8) :: a,b,c1,c2
    integer :: u,v,p
    integer :: i,n
    integer :: auxis,auxos
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(1), allocatable :: orbitalsA(:), orbitalsB(:)
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(1), allocatable :: couplingOrder(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    coupling = 0

    !! allocate arrays
    do n = 1, ConfigurationInteraction_instance%nproc
      do i = 1, numberOfSpecies

        call Matrix_constructorInteger ( ConfigurationInteraction_instance%couplingMatrix(i,n), &
          sum(ConfigurationInteraction_instance%numberOfStrings(i)%values), 3_8 , 0)

        call Matrix_constructorInteger(ConfigurationInteraction_instance%nCouplingOneTwo(i,n), &
          3_8, int(size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim=1),8),  0 )
  
        call Matrix_constructorInteger(ConfigurationInteraction_instance%nCouplingSize(i,n), &
          3_8, int(size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim=1) + 1 ,8),  0 )
  
        call Vector_constructor(ConfigurationInteraction_instance%couplingMatrixEnergyOne(i,n), &
          int(sum(ConfigurationInteraction_instance%numberOfStrings(i)%values),4), 0.0_8 )
  
        call Vector_constructorInteger(ConfigurationInteraction_instance%couplingMatrixFactorOne(i,n), &
          int(sum(ConfigurationInteraction_instance%numberOfStrings(i)%values),4), 2 )
  
        call Vector_constructorInteger( ConfigurationInteraction_instance%couplingMatrixOrbOne(i,n), &
          int(sum(ConfigurationInteraction_instance%numberOfStrings(i)%values),4), 0 )

      end do  
    end do  

  end subroutine ConfigurationInteraction_buildCouplingMatrix

  function ConfigurationInteraction_calculateEnergyOneSame( n, ii, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: a, b
    integer :: i,j,s,n, nn,ii
    integer :: l,k,z,kk,ll
    integer :: factor, factor2, auxOcc, AA, BB
    logical(1) :: equalA, equalB
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    integer :: diffOrb(2), otherdiffOrb(2) !! to avoid confusions
    real(8) :: auxCIenergy

    auxCIenergy = 0.0_8
    factor = 1

    !! copy a
    a = thisA(ii)
    b = thisB(ii)

    diffOrb = 0

    do kk = 1, ConfigurationInteraction_instance%occupationNumber(ii) 
      if ( ConfigurationInteraction_instance%orbitals(ii)%values( &
             ConfigurationInteraction_instance%strings(ii)%values(kk,a),b) == 0 ) then
        diffOrb(1) =  ConfigurationInteraction_instance%strings(ii)%values(kk,a)
        AA = kk
        exit
      end if
    end do

    do kk = 1, ConfigurationInteraction_instance%occupationNumber(ii) 
      if ( ConfigurationInteraction_instance%orbitals(ii)%values( &
             ConfigurationInteraction_instance%strings(ii)%values(kk,b),a) == 0 ) then
        diffOrb(2) =  ConfigurationInteraction_instance%strings(ii)%values(kk,b)
        BB = kk
        exit
      end if
    end do

    factor = (-1)**(AA-BB)

    configurationInteraction_instance%couplingMatrixFactorOne(ii,n)%values(b) = factor

    !One particle terms

    auxCIenergy= auxCIenergy +  ConfigurationInteraction_instance%twoCenterIntegrals(ii)%values( diffOrb(1), diffOrb(2) )

    !! save the different orbitals

    auxIndex1= ConfigurationInteraction_instance%twoIndexArray(ii)%values( diffOrb(1), diffOrb(2))
    ConfigurationInteraction_instance%couplingMatrixOrbOne(ii,n)%values(b) = auxIndex1

    do ll=1, ConfigurationInteraction_instance%occupationNumber( ii ) !! the same orbitals pair are excluded by the exchange

      l = ConfigurationInteraction_instance%strings(ii)%values(ll,b) !! or a

      auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(ii)%values( l,l) 
      auxIndex = ConfigurationInteraction_instance%fourIndexArray(ii)%values( auxIndex1, auxIndex2 )

      auxCIenergy = auxCIenergy + ConfigurationInteraction_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)

      auxIndex = ConfigurationInteraction_instance%fourIndexArray(ii)%values( &
                             ConfigurationInteraction_instance%twoIndexArray(ii)%values(diffOrb(1),l), &
                             ConfigurationInteraction_instance%twoIndexArray(ii)%values(l,diffOrb(2)) ) 

      auxCIenergy = auxCIenergy + &
                     MolecularSystem_instance%species(ii)%kappa*ConfigurationInteraction_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)

    end do

    !end if

    auxCIenergy= auxCIenergy * factor

  end function ConfigurationInteraction_calculateEnergyOneSame

  function ConfigurationInteraction_calculateEnergyOneDiff( ii, thisB, nn ) result (auxCIenergy)
    implicit none
    integer(8) :: thisB(:)
    integer(8) :: b
    integer :: i,j,ii, nn
    integer :: l,ll
    integer :: factor
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer :: auxIndex1, auxIndex11, auxIndex
    real(8) :: auxCIenergy

    auxCIenergy = 0.0_8

    b = thisB(ii)

    auxIndex1 = ConfigurationInteraction_instance%couplingMatrixOrbOne(ii,nn)%values(b) 
    factor = ConfigurationInteraction_instance%couplingMatrixFactorOne(ii,nn)%values(b) 

    do j=1, ii - 1 !! avoid ii, same species

      b = thisB(j)

      auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 
      auxIndex11 = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) 

      do ll=1,  ConfigurationInteraction_instance%occupationNumber( j ) 

        l = ConfigurationInteraction_instance%strings(j)%values(ll,b)

        auxIndex = auxIndex11  + ConfigurationInteraction_instance%twoIndexArray(j)%values( l,l) 

        auxCIenergy = auxCIenergy + &
        ConfigurationInteraction_instance%fourCenterIntegrals(ii,j)%values(auxIndex, 1) 

      end do

    end do

    do j= ii + 1, MolecularSystem_instance%numberOfQuantumSpecies!! avoid ii, same species

      b = thisB(j)

      auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

      auxIndex11 = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) 

      do ll=1,  ConfigurationInteraction_instance%occupationNumber( j )

        l = ConfigurationInteraction_instance%strings(j)%values(ll,b)

        auxIndex = auxIndex11  + ConfigurationInteraction_instance%twoIndexArray(j)%values( l,l) 

        auxCIenergy = auxCIenergy + &
        ConfigurationInteraction_instance%fourCenterIntegrals(ii,j)%values(auxIndex, 1) 
      end do

    end do

    auxCIenergy= auxCIenergy * factor

  end function ConfigurationInteraction_calculateEnergyOneDiff


  function ConfigurationInteraction_calculateEnergyTwoSame( ii, a, b ) result (auxCIenergy)
    implicit none
    integer(8) :: a, b
    integer :: ii
    integer :: kk,z
    integer :: factor, AA(2), BB(2)
    integer(8) :: auxIndex
    integer :: diffOrbA(2), diffOrbB(2)  !! to avoid confusions
    real(8) :: auxCIenergy

    !diffOrbA = 0
    !diffOrbB = 0
    z = 0

    do kk = 1, ConfigurationInteraction_instance%occupationNumber(ii) 
      if ( configurationinteraction_instance%orbitals(ii)%values( &
             configurationinteraction_instance%strings(ii)%values(kk,a),b) == 0 ) then
        z = z + 1
        diffOrbA(z) = ConfigurationInteraction_instance%strings(ii)%values(kk,a)
        AA(z) = kk
        if ( z == 2 ) exit
      end if
    end do

    z = 0
    do kk = 1, ConfigurationInteraction_instance%occupationNumber(ii) 
      if ( ConfigurationInteraction_instance%orbitals(ii)%values( &
             ConfigurationInteraction_instance%strings(ii)%values(kk,b),a) == 0 ) then
        z = z + 1
        diffOrbB(z) =  ConfigurationInteraction_instance%strings(ii)%values(kk,b)
        BB(z) = kk
        if ( z == 2 ) exit
      end if
    end do

    factor = (-1)**(AA(1)-BB(1) + AA(2) - BB(2) )
    auxIndex = ConfigurationInteraction_instance%fourIndexArray(ii)%values( &
                ConfigurationInteraction_instance%twoIndexArray(ii)%values(&
                  diffOrbA(1),diffOrbB(1)),&
                ConfigurationInteraction_instance%twoIndexArray(ii)%values(&
                  diffOrbA(2),diffOrbB(2)) )

    auxCIenergy = ConfigurationInteraction_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)

    auxIndex = ConfigurationInteraction_instance%fourIndexArray(ii)%values( &
                 ConfigurationInteraction_instance%twoIndexArray(ii)%values(&
                   diffOrbA(1),diffOrbB(2)),&
                 ConfigurationInteraction_instance%twoIndexArray(ii)%values(&
                   diffOrbA(2),diffOrbB(1)) )
    auxCIenergy = auxCIenergy + &
                MolecularSystem_instance%species(ii)%kappa*ConfigurationInteraction_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)

    auxCIenergy= auxCIenergy * factor

  end function ConfigurationInteraction_calculateEnergyTwoSame

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildDiagonal()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: ci
    integer :: i, j, ii, jj
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:), auxcilevel(:), dd(:)

!$  timeA = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    coupling = 0
    CIenergy = 0
    s = 0
    c = 0
    numberOfConfigurations = 0

    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( auxciLevel ( numberOfSpecies ) )
    allocate ( dd ( numberOfSpecies ) )

    ciLevel = 0
    auxciLevel = 0

    !!auxnumberOfSpecies = ConfigurationInteraction_numberOfConfigurationsRecursion2(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    numberOfConfigurations = 0
    ciLevel = 0

    !! call recursion to get the number of configurations...
    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = ConfigurationInteraction_numberOfConfigurationsRecursion(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    end do

    call Vector_constructor8 ( ConfigurationInteraction_instance%diagonalHamiltonianMatrix2, &
                              numberOfConfigurations, 0.0_8 ) 

    ConfigurationInteraction_instance%numberOfConfigurations = numberOfConfigurations 

    write (*,*) "Number Of Configurations: ", numberOfConfigurations

    allocate ( indexConf ( numberOfSpecies ) )
    indexConf = 0

    !! calculate the diagonal 
    s = 0
    c = 0
    ciLevel = 0

    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      s = 0
      dd = 0

      u = ConfigurationInteraction_instance%auxciOrderList(ci)
      auxnumberOfSpecies = ConfigurationInteraction_buildDiagonalRecursion( s, numberOfSpecies, indexConf,  c, dd, u, cilevel, auxcilevel )
    end do
    !stop

    deallocate ( dd )
    deallocate ( indexConf )
    deallocate ( ciLevel )
    deallocate ( auxciLevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Building diagonal of CI matrix : ", timeB - timeA ," (s)"

    write (*,*) "Reference energy, H_0: ",  ConfigurationInteraction_instance%diagonalHamiltonianMatrix2%values(1)

  end subroutine ConfigurationInteraction_buildDiagonal

recursive  function ConfigurationInteraction_numberOfConfigurationsRecursion(s, numberOfSpecies, c, cilevel) result (os)
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies
    integer :: os,is
    integer :: cilevel(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      i = cilevel(is) + 1
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        os = ConfigurationInteraction_numberOfConfigurationsRecursion( is, numberOfSpecies, c, cilevel )
      end do
    else 
      os = is

      i = cilevel(is) + 1
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        c = c + 1
      end do
    end if

  end function ConfigurationInteraction_numberOfConfigurationsRecursion

recursive  function ConfigurationInteraction_buildDiagonalRecursion(s, numberOfSpecies, indexConf, c, dd, u, cilevel, auxcilevel) result (os)
    implicit none

    integer(8) :: a,b,c,cc,d
    integer :: u,v
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies
    integer :: os,is
    integer :: size1, size2
    integer(8) :: indexConf(:)
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer :: ssize
    integer :: cilevel(:), auxcilevel(:), dd(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      i = cilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a

        dd(is) =(a + ConfigurationInteraction_instance%ciOrderSize1(u,is))* ConfigurationInteraction_instance%ciOrderSize2(u,is) 
        os = ConfigurationInteraction_buildDiagonalRecursion( is, numberOfSpecies, indexConf, c, dd, u, cilevel, auxcilevel )
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        !print *, indexConf
        dd(is) =(a + ConfigurationInteraction_instance%ciOrderSize1(u,is))* ConfigurationInteraction_instance%ciOrderSize2(u,is) 
        d = sum(dd)

        ConfigurationInteraction_instance%diagonalHamiltonianMatrix2%values(c) = &
                              ConfigurationInteraction_calculateEnergyZero ( indexConf )

      end do
    end if

  end function ConfigurationInteraction_buildDiagonalRecursion

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<

  !! Map the indexes of initial CI matrix to the complete matrix.
  subroutine ConfigurationInteraction_getInitialIndexes()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: ci
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:)

!$  timeA = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    s = 0
    c = 0

    call Matrix_constructorInteger ( ConfigurationInteraction_instance%auxConfigurations, int( numberOfSpecies,8), &
          int(CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX,8), 0 )

    !! call recursion

    allocate ( cilevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )

    s = 0
    c = 0
    indexConf = 0
    cilevel = 0

    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = ConfigurationInteraction_getIndexesRecursion( s, numberOfSpecies, indexConf, c, cilevel )
    end do

    deallocate ( indexConf )
    deallocate ( cilevel )

!$  timeB = omp_get_wtime()

!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for getting initial indexes : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_getInitialIndexes

recursive  function ConfigurationInteraction_getIndexesRecursion(s, numberOfSpecies, indexConf, c, cilevel) result (os)
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: i, j, ii, jj
    integer :: s, ss, numberOfSpecies
    integer :: os,is
    integer :: size1, size2
    integer(8) :: indexConf(:)
    integer(1) :: coupling
    integer :: ssize
    integer :: cilevel(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      i = cilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)

      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a
        os = ConfigurationInteraction_getIndexesRecursion( is, numberOfSpecies, indexConf, c, cilevel)
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)

      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        do u = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
          if ( c ==  ConfigurationInteraction_instance%auxIndexCIMatrix%values(u) ) then
            do ss = 1, numberOfSpecies
              ConfigurationInteraction_instance%auxConfigurations%values(ss,u) = indexConf(ss) 
            end do
          end if
        end do
      end do
    end if

  end function ConfigurationInteraction_getIndexesRecursion

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_calculateInitialCIMatrix()
    implicit none

    integer(8) :: a,b,aa,bb
    integer :: u,v
    integer :: i
    integer :: numberOfSpecies
    real(8) :: timeA1, timeB1
    integer(1) :: coupling
    integer(1), allocatable :: orbitalsA(:), orbitalsB(:)
    integer :: initialCIMatrixSize 
    integer :: nproc
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    initialCIMatrixSize = CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 

    allocate ( indexConfA ( numberOfSpecies ) )
    allocate ( indexConfB ( numberOfSpecies ) )

!$    timeA1 = omp_get_wtime()

    do a = 1, initialCIMatrixSize 
      aa = ConfigurationInteraction_instance%auxIndexCIMatrix%values(a)
      do b = a, initialCIMatrixSize 
        bb = ConfigurationInteraction_instance%auxIndexCIMatrix%values(b)
        coupling = 0

        indexConfA = 0
        indexConfB = 0

        do i = 1, numberOfSpecies

          allocate (orbitalsA ( ConfigurationInteraction_instance%numberOfOrbitals%values(i) ))
          allocate (orbitalsB ( ConfigurationInteraction_instance%numberOfOrbitals%values(i) ))
          orbitalsA = 0
          orbitalsB = 0

          indexConfA(i) =  ConfigurationInteraction_instance%auxConfigurations%values(i,a)
          indexConfB(i) =  ConfigurationInteraction_instance%auxConfigurations%values(i,b)

          do u = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
            orbitalsA( ConfigurationInteraction_instance%strings(i)%values(u,indexConfA(i) ) ) = 1
          end do
          do v = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
            orbitalsB( ConfigurationInteraction_instance%strings(i)%values(v,indexConfB(i) ) ) = 1
          end do
          coupling = coupling + &
            ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) - sum ( orbitalsA * orbitalsB ) 

          deallocate (orbitalsA )
          deallocate (orbitalsB )

        end do
        if ( coupling  == 0 ) then
          ConfigurationInteraction_instance%initialHamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_instance%diagonalHamiltonianMatrix2%values(a) 

        else if (  coupling == 1 ) then

          ConfigurationInteraction_instance%initialHamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_calculateEnergyOne ( 1, indexConfA, indexConfB )

        else if ( coupling  == 2 ) then

          ConfigurationInteraction_instance%initialHamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_calculateEnergyTwo ( 1, indexConfA, indexConfB )

        end if
 

      end do


    end do

    deallocate ( indexConfB )
    deallocate ( indexConfA )

!$  timeB1 = omp_get_wtime()
    !! symmetrize
    do a = 1, initialCIMatrixSize 
      do b = a, initialCIMatrixSize 

         ConfigurationInteraction_instance%initialHamiltonianMatrix%values(b,a) = &
            ConfigurationInteraction_instance%initialHamiltonianMatrix%values(a,b) 
      end do
    end do

    !!open(unit=318, file="cimatrix.dat", action = "write", form="formatted")
    !!do a = 1, initialCIMatrixSize 
    !!  do b = 1, initialCIMatrixSize 
    !!     write (318,*) a,b, ConfigurationInteraction_instance%initialHamiltonianMatrix%values(a,b)
    !!  end do
    !!     write (318,*) " "
    !!end do
    !!close(318)
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Calculating initial CI matrix : ", timeB1 - timeA1 ," (s)"

  end subroutine ConfigurationInteraction_calculateInitialCIMatrix 


  subroutine ConfigurationInteraction_buildInitialCIMatrix2()
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    type (Vector8) :: diagonalHamiltonianMatrix
    integer :: a,b,c,aa,bb,i
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer :: initialCIMatrixSize 
    integer :: nproc

    !$ timeA = omp_get_wtime()
    initialCIMatrixSize = CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
    if ( ConfigurationInteraction_instance%numberOfConfigurations <  CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX ) then
      CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX = ConfigurationInteraction_instance%numberOfConfigurations !! assign to an internal variable
    end if

    call Vector_constructorInteger8 ( ConfigurationInteraction_instance%auxIndexCIMatrix, &
                              ConfigurationInteraction_instance%numberOfConfigurations, 0_8 ) !hmm

    do a = 1, ConfigurationInteraction_instance%numberOfConfigurations
      ConfigurationInteraction_instance%auxIndexCIMatrix%values(a)= a
    end do

   !! save the unsorted diagonal Matrix
    call Vector_constructor8 ( ConfigurationInteraction_instance%diagonalHamiltonianMatrix, &
                              ConfigurationInteraction_instance%numberOfConfigurations, 0.0_8 ) 


    ConfigurationInteraction_instance%diagonalHamiltonianMatrix%values = ConfigurationInteraction_instance%diagonalHamiltonianMatrix2%values

   !! To get only the lowest 300 values.
   call Vector_reverseSortElements8( ConfigurationInteraction_instance%diagonalHamiltonianMatrix2, &
          ConfigurationInteraction_instance%auxIndexCIMatrix, int(initialCIMatrixSize,8))

   call Matrix_constructor ( ConfigurationInteraction_instance%initialHamiltonianMatrix, int(initialCIMatrixSize,8) , &
                               int(initialCIMatrixSize,8) , 0.0_8 ) 

    !! get the configurations for the initial hamiltonian matrix
    call ConfigurationInteraction_getInitialIndexes()

    call ConfigurationInteraction_calculateInitialCIMatrix()

    !! diagonalize the initial matrix
    call Vector_constructor8 ( ConfigurationInteraction_instance%initialEigenValues, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    call Matrix_constructor (ConfigurationInteraction_instance%initialEigenVectors, &
          int(initialCIMatrixSize,8), &
          int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

    call Matrix_eigen_select ( ConfigurationInteraction_instance%initialHamiltonianMatrix, &
           ConfigurationInteraction_instance%initialEigenValues, &
           1, int(CONTROL_instance%NUMBER_OF_CI_STATES,4), &  
           eigenVectors = ConfigurationInteraction_instance%initialEigenVectors, &
           flags = int(SYMMETRIC,4))
    
    write(*,*) "Initial eigenValues"
    do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
      write (*,*)  i, ConfigurationInteraction_instance%initialEigenValues%values(i)
    end do

    call Vector_destructor8 ( ConfigurationInteraction_instance%diagonalHamiltonianMatrix2 )

!$    timeB = omp_get_wtime()
!$    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Solving Initial CI : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_buildInitialCIMatrix2


  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildHamiltonianMatrix()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v,p
    integer :: ci
    integer :: i, j, ii, jj
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:)
    integer(8), allocatable :: pindexConf(:,:)
    integer, allocatable :: cilevel(:), auxcilevel(:), dd(:)
    integer(8), allocatable :: indexConfA(:,:)
    integer(8), allocatable :: indexConfB(:,:)
    integer, allocatable :: stringAinB(:)
    integer(1), allocatable :: couplingSpecies(:,:)
    integer :: n,nproc


!$  timeA = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    numberOfConfigurations = ConfigurationInteraction_instance%numberOfConfigurations 

    allocate ( ConfigurationInteraction_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )
    ciLevel = 0
    ConfigurationInteraction_instance%allIndexConf = 0
    indexConf = 0

    !! gather all configurations
    s = 0
    c = 0
    ciLevel = 0

    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 

      cilevel(:) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = ConfigurationInteraction_gatherConfRecursion( s, numberOfSpecies, indexConf,  c, cilevel )
    end do
    !stop

    deallocate ( indexConf )
    deallocate ( ciLevel )

    !! allocate the hamiltonian matrix
    call Matrix_constructor ( ConfigurationInteraction_instance%hamiltonianMatrix, & 
           int(ConfigurationInteraction_instance%numberOfConfigurations,8), &
           int(ConfigurationInteraction_instance%numberOfConfigurations,8), 0.0_8)


    nproc = omp_get_max_threads()
    !! calculate the matrix elements
    allocate ( indexConfA ( numberOfSpecies, nproc ) )
    allocate ( indexConfB ( numberOfSpecies, nproc ) )
    allocate ( pindexConf ( numberOfSpecies, nproc ) )
    allocate ( couplingSpecies ( numberOfSpecies, nproc ) )

    indexConfA = 0
    indexConfB = 0
    pindexConf = 0
    couplingSpecies = 0

!$omp parallel & 
!$omp& private(a,b,coupling,i,p,stringAinB,n),&
!$omp& shared(ConfigurationInteraction_instance, HartreeFock_instance)
    n = omp_get_thread_num() + 1
!$omp do schedule (dynamic) 
    do a = 1, numberOfConfigurations
      indexConfA(:,n) = ConfigurationInteraction_instance%allIndexConf(:,a) 
      do b = a, numberOfConfigurations

        indexConfB(:,n) = ConfigurationInteraction_instance%allIndexConf(:,b) 

        do i = 1, numberOfSpecies
          if ( pindexConf(i,n) /= indexConfB(i,n) ) then
            allocate (stringAinB (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) ))
            stringAinB = 0
            do p = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
              stringAinB(p) = ConfigurationInteraction_instance%orbitals(i)%values( &
                                ConfigurationInteraction_instance%strings(i)%values(p,indexConfA(i,n) ), indexConfB(i,n) ) 
            end do
            couplingSpecies(i,n) = configurationinteraction_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB )
            deallocate (stringAinB )
          end if
        end do
        coupling = sum(couplingSpecies(:,n))

        if ( coupling  == 0 ) then
          ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_instance%diagonalHamiltonianMatrix2%values(a) 
  
        else if (  coupling == 1 ) then
  
          ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_calculateEnergyOne ( n, indexConfA(:,n), indexConfB(:,n) )
  
        else if ( coupling  == 2 ) then
  
          ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) = &
            ConfigurationInteraction_calculateEnergyTwo ( n, indexConfA(:,n), indexConfB(:,n) )
  
        end if

        pindexConf(:,n) = indexConfB(:,n)

      end do
        pindexConf(:,n) = 0
    end do
    !$omp end do nowait
    !$omp end parallel

    deallocate ( pindexConf )
    deallocate ( couplingSpecies )
    deallocate ( indexConfB )
    deallocate ( indexConfA )

    !! symmetrize
    do a = 1, numberOfConfigurations
      do b  = a, numberOfConfigurations
         ConfigurationInteraction_instance%hamiltonianMatrix%values(b,a) = &
            ConfigurationInteraction_instance%hamiltonianMatrix%values(a,b) 
      end do
    end do

    deallocate ( ConfigurationInteraction_instance%allIndexConf )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"

  end subroutine ConfigurationInteraction_buildHamiltonianMatrix

recursive  function ConfigurationInteraction_gatherConfRecursion(s, numberOfSpecies, indexConf, c, cilevel ) result (os)
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
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a
        os = ConfigurationInteraction_gatherConfRecursion( is, numberOfSpecies, indexConf, c, cilevel )
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        ConfigurationInteraction_instance%allIndexConf(:,c) = indexConf

      end do
    end if

  end function ConfigurationInteraction_gatherConfRecursion

recursive  function ConfigurationInteraction_buildMatrixRecursion(nproc, s, indexConf, auxindexConf, cc, c, n, v, w, &
                                                                  cilevel, auxcilevel) result (os)
    implicit none

    integer(8) :: a,c,aa
    integer :: i, n, nn, nproc 
    integer :: s, numberOfSpecies
    integer :: os,is,ss,ssize
    integer(8) :: cc(:)
    integer(8) :: indexConf(:,:)
    integer(8) :: auxindexConf(:,:)
    real(8) :: v(:)
    real(8) :: w(:)
    integer :: cilevel(:,:)
    integer :: auxcilevel(:,:)

    is = s + 1
    !if ( is < numberOfSpecies ) then
    do ss = 1, ConfigurationInteraction_instance%recursionVector1(is) 
      i = cilevel(is,n) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)

      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        indexConf(is,n:) = ssize + a
        os = ConfigurationInteraction_buildMatrixRecursion( nproc, is, indexConf, auxindexConf, cc, c, n, v, w, cilevel, auxcilevel )
      end do
    end do
    !else 
    do ss = 1, ConfigurationInteraction_instance%recursionVector2(is) 
      os = is
      i = cilevel(is,n) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)

      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        c = c + 1

        if ( abs(v(c)) > CONTROL_instance%CI_MATVEC_TOLERANCE ) then
          cc(n) = c 
          indexConf(is,n:) = ssize + a

          auxindexConf = indexConf
          auxcilevel = cilevel

          if ( n == nproc ) then

            !$omp parallel &
            !$omp& private(nn),&
            !$omp& shared(v,w, indexConf, cc, nproc, cilevel) 
            !$omp do schedule (static) 
            do nn = 1, nproc
              call ConfigurationInteraction_buildRow( nn, indexConf(:,nn), cc(nn), w, v(cc(nn)), cilevel(:,nn))
            end do
            !$omp end do nowait
            !$omp end parallel
            n = 0 

            do nn = 1, nproc
              indexConf(:,nn) = indexConf(:,nproc) 
              cilevel(:,nn) = cilevel(:,nproc) 
            end do
          end if 

          n = n + 1

        end if

      end do
    end do
    !end if


  end function ConfigurationInteraction_buildMatrixRecursion

  !! Alternative option to the recursion with the same computational cost... However, it may be helpul some day. 

  function ConfigurationInteraction_buildMatrixRecursion2(nproc, s, indexConf, auxindexConf, cc, c, n, v, w, &
                                                                  cilevel, auxcilevel) result (os)
    implicit none

    integer(8) :: a,c,aa, x
    integer :: i, j, n, nn, nproc, ci
    integer :: s, numberOfSpecies
    integer :: os,is,ss,ssize
    integer(8) :: cc(:)
    integer(8) :: indexConf(:,:)
    integer(8) :: auxindexConf(:,:)
    real(8) :: v(:)
    real(8) :: w(:)
    integer :: cilevel(:,:)
    integer(8) :: totalsize, auxtotalsize
    integer :: auxcilevel(:,:)
    integer, allocatable :: counter(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    
    allocate (counter(numberOfSpecies))
    counter = 0 

    totalsize = 1
    do i = 1 , numberOfSpecies
      totalsize = totalsize * ConfigurationInteraction_instance%numberOfStrings(i)%values(cilevel(i,n) + 1)
    end do

    do i = 1 , numberOfSpecies 
      ci = cilevel(i,n) + 1 
      ssize = ConfigurationInteraction_instance%numberOfStrings2(i)%values(ci)
      indexConf(i,n:) = ssize  + 1
    end do

    indexConf(numberOfSpecies,n:) = indexConf(numberOfSpecies,n:) -1

    do x = 1, totalsize

      indexConf(numberOfSpecies,n:) = indexConf(numberOfSpecies,n:) + 1

      do i = numberOfSpecies, 1 + 1, -1 
        auxtotalsize = 1
        do j = i, numberOfSpecies
          auxtotalsize = auxtotalsize * ConfigurationInteraction_instance%numberOfStrings(j)%values(cilevel(j,n) + 1)
        end do
        if (counter(i) == auxtotalsize) then
          do j = i, numberOfSpecies
            ci = cilevel(j,n) + 1 
            ssize = ConfigurationInteraction_instance%numberOfStrings2(j)%values(ci)
            indexConf(j,n:) = ssize + 1
          end do
          counter(i) = 0
          indexConf(i-1,n:) = indexConf(i-1,n:) + 1
        end if
        counter(i) = counter(i) + 1 

      end do
      !print *, indexConf(:,1)
    end do

    deallocate (counter)

  end function ConfigurationInteraction_buildMatrixRecursion2


  subroutine ConfigurationInteraction_buildRow( nn, indexConfA, c, w, vc, cilevelA)
    implicit none

    integer(8) :: a,b,c,bb,ci,d,cj
    integer :: u,v,uu,vv, p, nn
    integer :: i, j, auxis,auxos,is, ii, aa
    integer :: numberOfSpecies, s
    integer, allocatable :: stringAinB(:)
    integer(4) :: coupling
    integer(4) :: ssize,auxcoupling(3) !! 0,1,2
    integer(8) :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(8), allocatable :: dd(:)
    real(8) :: vc, CIenergy
    real(8) :: w(:)
    integer :: cilevelA(:)
    integer, allocatable :: cilevel(:)


    !ConfigurationInteraction_instance%pindexConf = 0

    !!$ ConfigurationInteraction_instance%timeA(1) = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do i = 1, numberOfSpecies 

      if ( ConfigurationInteraction_instance%pindexConf(i,nn) /= indexConfA(i) ) then

      ConfigurationInteraction_instance%nCouplingOneTwo(i,nn)%values = 0
      auxcoupling = 0

      !allocate (stringBinA (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) ))
      allocate (stringAinB (ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i) ))

      stringAinB = 0
      !stringBinA = 0

      a = indexConfA(i)

    !!$ ConfigurationInteraction_instance%timeA(2) = omp_get_wtime()

      ssize = 0 
      do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1)
        do b = 1 + ssize , ConfigurationInteraction_instance%numberOfStrings(i)%values(ci) + ssize

          !b = ssize + bb
          do p = ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i)+1, &
                 ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)
          !do p = 1, &
          !       ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)

            stringAinB(p) = ConfigurationInteraction_instance%orbitals(i)%values( &
                              ConfigurationInteraction_instance%strings(i)%values(p,a),b) 

            !stringBinA(p) = ConfigurationInteraction_instance%orbitals(i)%values( &
            !                  ConfigurationInteraction_instance%strings(i)%values(p,b),a) 
          end do

          coupling = configurationinteraction_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB ) - &
                     ConfigurationInteraction_instance%numberOfCoreOrbitals%values(i) 

         ! coupling = configurationinteraction_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB ) 

          if ( coupling  <= 2 ) then

            coupling = coupling + 1

            auxcoupling(coupling) = auxcoupling(coupling) + 1 

            ConfigurationInteraction_instance%nCouplingOneTwo(i,nn)%values( coupling, ci) = &
              ConfigurationInteraction_instance%nCouplingOneTwo(i,nn)%values( coupling, ci) + 1

            ConfigurationInteraction_instance%couplingMatrix(i,nn)%values( auxcoupling(coupling), coupling ) = b
          end if

        end do

        ssize = ssize + ConfigurationInteraction_instance%numberOfStrings(i)%values(ci)

        end do

      deallocate (stringAinB)
      !deallocate (stringBinA)
      end if

    end do

    !!$ ConfigurationInteraction_instance%timeB(1) = omp_get_wtime()

    do is = 1, numberOfSpecies
      do i = 1, 3 !! 0,1,2
        ssize = 0
        do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(is)%values, dim = 1) !! 1 is always zero
          ssize = ssize + ConfigurationInteraction_instance%nCouplingOneTwo(is,nn)%values( i,ci ) 
          ConfigurationInteraction_instance%nCouplingSize(is,nn)%values( i,ci+1 ) = ssize
         end do
        ConfigurationInteraction_instance%nCouplingSize(is,nn)%values( i,1 ) = 0 !0?
      end do
   end do


    !!$ ConfigurationInteraction_instance%timeA(2) = omp_get_wtime()
    allocate ( indexConfB ( numberOfSpecies ) )
    allocate ( cilevel ( numberOfSpecies ) )
    allocate ( dd ( numberOfSpecies ) )
    indexConfB = 0

    !!$ ConfigurationInteraction_instance%timeB(2) = omp_get_wtime()
    !!$ ConfigurationInteraction_instance%timeA(3) = omp_get_wtime()

    !!one diff same species
    do i = 1, numberOfSpecies

      if ( ConfigurationInteraction_instance%pindexConf(i,nn) /= indexConfA(i) ) then
      cilevel(:) = 0
      indexConfB = indexConfA

      cilevel = cilevelA

      do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        auxos = ConfigurationInteraction_buildRowRecursionFirstOne( i, indexConfA, indexConfB, nn, cilevel )

      end do      
      end if
    end do      

    !!$ ConfigurationInteraction_instance%timeB(3) = omp_get_wtime()

    !!$ ConfigurationInteraction_instance%timeA(4) = omp_get_wtime()

    !$omp atomic
      w(c) = w(c) + vc*ConfigurationInteraction_instance%diagonalHamiltonianMatrix%values(c) 
    !$omp end atomic

    !!$ ConfigurationInteraction_instance%timeB(4) = omp_get_wtime()

    !!$ ConfigurationInteraction_instance%timeA(5) = omp_get_wtime()
    !! one diff
    do i = 1, numberOfSpecies
      cilevel(:) = 0
      indexConfB = indexConfA

      cilevel = cilevelA

      do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        do u = 1,  configurationinteraction_instance%sizeciorderlist 
          if ( sum(abs(cilevel - &
               configurationinteraction_instance%ciorderlist( configurationinteraction_instance%auxciorderlist(u), :))) == 0 ) then

            uu = configurationinteraction_instance%auxciorderlist(u)
            dd = 0

            auxos = ConfigurationInteraction_buildRowRecursionSecondOne( i, indexConfB, w, vc, dd, nn, cilevel, uu )
            exit

          end if
        end do
      end do      
    end do      

    !!$ ConfigurationInteraction_instance%timeB(5) = omp_get_wtime()
    !!$ ConfigurationInteraction_instance%timeA(6) = omp_get_wtime()

    !! two diff same species
    do i = 1, numberOfSpecies

      cilevel(:) = 0
      indexConfB = indexConfA
      cilevel = cilevelA

      do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
          if ( sum(abs(cilevel - &
               ConfigurationInteraction_instance%ciOrderList( ConfigurationInteraction_instance%auxciOrderList(u), :))) == 0 ) then
            uu = ConfigurationInteraction_instance%auxciOrderList(u)
            dd = 0

            if ( ConfigurationInteraction_instance%pindexConf(i,nn) /= indexConfA(i) ) then
              auxos = ConfigurationInteraction_buildRowRecursionSecondTwoCal( i, indexConfA, indexConfB, w, vc, dd, nn, cilevel, uu )
            else
              auxos = ConfigurationInteraction_buildRowRecursionSecondTwoGet( i, indexConfA, indexConfB, w, vc, dd, nn, cilevel, uu )
            end if

            exit

          end if
        end do
      end do
    end do      

    !!$ ConfigurationInteraction_instance%timeB(6) = omp_get_wtime()
    !!$ ConfigurationInteraction_instance%timeA(7) = omp_get_wtime()

    !! two diff diff species
    do v = 1, ConfigurationInteraction_instance%ncouplingOrderTwoDiff

       i = ConfigurationInteraction_instance%couplingOrderIndex(3,v)%values(1)
       j = ConfigurationInteraction_instance%couplingOrderIndex(3,v)%values(2)

       indexConfB = indexConfA
       cilevel = cilevelA

        do ci = 1,  size(ConfigurationInteraction_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
          cilevel(i) = ci - 1
          do cj = 1,  size(ConfigurationInteraction_instance%numberOfStrings(j)%values, dim = 1) !! 1 is always zero
            cilevel(j) = cj - 1
            do u = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
              if ( sum(abs(cilevel - &
                   ConfigurationInteraction_instance%ciOrderList( ConfigurationInteraction_instance%auxciOrderList(u), :))) == 0 ) then

                uu = ConfigurationInteraction_instance%auxciOrderList(u)
                dd = 0
                auxos = ConfigurationInteraction_buildRowRecursionSecondTwoDiff( i, j, indexConfB, w, vc, dd, nn, cilevel, uu )
                exit
              end if
            end do
          end do
        end do
    end do    

    !!$ ConfigurationInteraction_instance%timeB(7) = omp_get_wtime()

    !!$ print *, "omptime"
    !!$ print *, "1", ConfigurationInteraction_instance%timeB(1) - ConfigurationInteraction_instance%timeA(1)
    !!$ print *, "2", ConfigurationInteraction_instance%timeB(2) - ConfigurationInteraction_instance%timeA(2)
    !!$ print *, "3", ConfigurationInteraction_instance%timeB(3) - ConfigurationInteraction_instance%timeA(3)
    !!$ print *, "4", ConfigurationInteraction_instance%timeB(4) - ConfigurationInteraction_instance%timeA(4)
    !!$ print *, "5", ConfigurationInteraction_instance%timeB(5) - ConfigurationInteraction_instance%timeA(5)
    !!$ print *, "6", ConfigurationInteraction_instance%timeB(6) - ConfigurationInteraction_instance%timeA(6)
    !!$ print *, "7", ConfigurationInteraction_instance%timeB(7) - ConfigurationInteraction_instance%timeA(7)

    ConfigurationInteraction_instance%pindexConf(:,nn) = indexConfA(:)

    deallocate ( dd )
    deallocate ( cilevel )
    deallocate ( indexConfB )

  end subroutine ConfigurationInteraction_buildRow

recursive  function ConfigurationInteraction_buildRowRecursionFirstOne( ii, indexConfA, indexConfB, nn, cilevel ) result (os)
    implicit none

    integer(8) :: a, aa
    integer :: ii, nn, ci
    integer :: os, ssize
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

      ci = cilevel(ii) + 1
      ssize = ConfigurationInteraction_instance%nCouplingSize(ii,nn)%values( 2,ci ) 
      do aa = 1, ConfigurationInteraction_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
        a = ssize + aa

        indexConfB(ii) = ConfigurationInteraction_instance%couplingMatrix(ii,nn)%values(a, 2)
        CIenergy = ConfigurationInteraction_calculateEnergyOneSame ( nn, ii, indexConfA, indexConfB )
        ConfigurationInteraction_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) = CIenergy

      end do

  end function ConfigurationInteraction_buildRowRecursionFirstOne
 
recursive  function ConfigurationInteraction_buildRowRecursionSecondOne( ii, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: ii, nn, ci, u, j
    integer :: ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ci = cilevel(ii) + 1
    ssize = ConfigurationInteraction_instance%nCouplingSize(ii,nn)%values( 2,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - ConfigurationInteraction_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
                 ConfigurationInteraction_instance%ciOrderSize1(u,j) )* ConfigurationInteraction_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, ConfigurationInteraction_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
      a = ssize + aa

      indexConfB(ii) = ConfigurationInteraction_instance%couplingMatrix(ii,nn)%values(a, 2)

      dd(ii) = (indexConfB(ii) - ConfigurationInteraction_instance%numberOfStrings2(ii)%values(ci) + &
                 ConfigurationInteraction_instance%ciOrderSize1(u,ii) )* ConfigurationInteraction_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      CIenergy = ConfigurationInteraction_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      CIenergy = CIenergy + ConfigurationInteraction_calculateEnergyOneDiff ( ii, indexConfB, nn )
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function ConfigurationInteraction_buildRowRecursionSecondOne

  function ConfigurationInteraction_buildRowRecursionSecondTwoCal( ii, indexConfA, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: i, ii, nn, ci, u, j
    integer :: s, ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ci = cilevel(ii) + 1
    ssize = ConfigurationInteraction_instance%nCouplingSize(ii,nn)%values( 3,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - ConfigurationInteraction_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
               ConfigurationInteraction_instance%ciOrderSize1(u,j) )* ConfigurationInteraction_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, ConfigurationInteraction_instance%nCouplingOneTwo(ii,nn)%values( 3,ci ) 
      a = ssize + aa

      indexConfB(ii) = ConfigurationInteraction_instance%couplingMatrix(ii,nn)%values(a, 3)
      dd(ii) = (indexConfB(ii) - ConfigurationInteraction_instance%numberOfStrings2(ii)%values(ci) + &
               ConfigurationInteraction_instance%ciOrderSize1(u,ii) )* ConfigurationInteraction_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      !CIenergy = ConfigurationInteraction_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      CIenergy = ConfigurationInteraction_calculateEnergyTwoSame ( ii, indexConfA(ii), indexConfB(ii) )
      ConfigurationInteraction_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) = CIenergy
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function ConfigurationInteraction_buildRowRecursionSecondTwoCal

  function ConfigurationInteraction_buildRowRecursionSecondTwoGet( ii, indexConfA, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: i, ii, nn, ci, u, j
    integer :: s, ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ci = cilevel(ii) + 1
    ssize = ConfigurationInteraction_instance%nCouplingSize(ii,nn)%values( 3,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - ConfigurationInteraction_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
                 ConfigurationInteraction_instance%ciOrderSize1(u,j) )* ConfigurationInteraction_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, ConfigurationInteraction_instance%nCouplingOneTwo(ii,nn)%values( 3,ci ) 
      a = ssize + aa

      indexConfB(ii) = ConfigurationInteraction_instance%couplingMatrix(ii,nn)%values(a, 3)
      dd(ii) = (indexConfB(ii) - ConfigurationInteraction_instance%numberOfStrings2(ii)%values(ci) + &
               ConfigurationInteraction_instance%ciOrderSize1(u,ii) )* ConfigurationInteraction_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      CIenergy = ConfigurationInteraction_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      !CIenergy = ConfigurationInteraction_calculateEnergyTwoSame ( ii, indexConfA(ii), indexConfB(ii) )
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function ConfigurationInteraction_buildRowRecursionSecondTwoGet

 function ConfigurationInteraction_buildRowRecursionSecondTwoDiff( ii, jj, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: ai,aj,d, aai, aaj
    integer :: ii, nn, ci, u, k, jj, cj
    integer :: ssizei, ssizej
    integer :: bi, bj, factor, factori
    integer :: auxIndex1, auxIndex2, auxIndex
    integer :: os,numberOfSpecies
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ci = cilevel(ii) + 1
    cj = cilevel(jj) + 1
    ssizei = ConfigurationInteraction_instance%nCouplingSize(ii,nn)%values( 2,ci ) 
    ssizej = ConfigurationInteraction_instance%nCouplingSize(jj,nn)%values( 2,cj ) 

    do k = 1, numberOfSpecies
        dd(k) = (indexConfB(k) - ConfigurationInteraction_instance%numberOfStrings2(k)%values(cilevel(k)+1) + &
                ConfigurationInteraction_instance%ciOrderSize1(u,k) )* ConfigurationInteraction_instance%ciOrderSize2(u,k) 
    end do

    do aai = 1, ConfigurationInteraction_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
      ai = ssizei + aai
      indexConfB(ii) = ConfigurationInteraction_instance%couplingMatrix(ii,nn)%values(ai, 2)
      dd(ii) = (indexConfB(ii) - ConfigurationInteraction_instance%numberOfStrings2(ii)%values(ci) + &
                ConfigurationInteraction_instance%ciOrderSize1(u,ii) )* ConfigurationInteraction_instance%ciOrderSize2(u,ii) 

      bi = indexConfB(ii)
      factori = ConfigurationInteraction_instance%couplingMatrixFactorOne(ii,nn)%values(bi) 
      auxIndex1 = ConfigurationInteraction_instance%couplingMatrixOrbOne(ii,nn)%values(bi) 
      auxIndex1 = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(jj) * (auxIndex1 - 1 ) 

      do aaj = 1, ConfigurationInteraction_instance%nCouplingOneTwo(jj,nn)%values( 2,cj ) 
        aj = ssizej + aaj
        indexConfB(jj) = ConfigurationInteraction_instance%couplingMatrix(jj,nn)%values(aj, 2)

        dd(jj) = (indexConfB(jj) - ConfigurationInteraction_instance%numberOfStrings2(jj)%values(cj) + &
                ConfigurationInteraction_instance%ciOrderSize1(u,jj) )* ConfigurationInteraction_instance%ciOrderSize2(u,jj) 

        d = sum(dd)
          !CIenergy = vc*ConfigurationInteraction_calculateEnergyTwoDiff ( ii, jj, indexConfB, nn )

        bj = indexConfB(jj)
        factor = factori * ConfigurationInteraction_instance%couplingMatrixFactorOne(jj,nn)%values(bj) 
        auxIndex2 = ConfigurationInteraction_instance%couplingMatrixOrbOne(jj,nn)%values(bj) 
        auxIndex = auxIndex1 + auxIndex2

        CIenergy = vc * factor *ConfigurationInteraction_instance%fourCenterIntegrals(ii,jj)%values(auxIndex, 1)
        !CIenergy = vc*CIenergy 

        !$omp atomic
        w(d) = w(d) + CIenergy 
        !$omp end atomic
      end do
    end do

  end function ConfigurationInteraction_buildRowRecursionSecondTwoDiff



  function ConfigurationInteraction_getIndex ( indexConf ) result ( output )
    implicit none
    integer(8) :: indexConf(:)
    integer(8) :: output, ssize
    integer :: i,j, numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    output = 0 
     !! simplify!!
    do i = 1, numberOfSpecies
      ssize = 1 
      do j = i + 1, numberOfSpecies
        ssize = ssize * ConfigurationInteraction_instance%sumstrings(j)
        !ssize = ssize * sum(ConfigurationInteraction_instance%numberOfStrings(j)%values(1:2))
      end do
      output = output + ( indexConf(i) - 1 ) * ssize
    end do
    output = output + 1

  end function ConfigurationInteraction_getIndex

recursive  function ConfigurationInteraction_getIndexSize(s, c, auxcilevel) result (os)
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: i, j, ii, jj, ss
    integer :: s, numberOfSpecies
    integer :: os,is,cc, ssize
    integer :: auxcilevel(:)

    is = s + 1
    do ss = 1, ConfigurationInteraction_instance%recursionVector1(is) 
      i = auxcilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      do a = 1, ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
        os = ConfigurationInteraction_getIndexSize( is, c, auxcilevel )
      end do
    end do
    do ss = 1, ConfigurationInteraction_instance%recursionVector2(is) 
      os = is
      i = auxcilevel(is) + 1
      ssize = ConfigurationInteraction_instance%numberOfStrings2(is)%values(i)
      c = c + ConfigurationInteraction_instance%numberOfStrings(is)%values(i)
    end do

  end function ConfigurationInteraction_getIndexSize

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_buildAndSaveCIMatrix()
    implicit none
    type(Configuration) :: auxConfigurationA, auxConfigurationB
    integer(8) :: a,b,c,d,n, nproc,cc
    real(8) :: timeA, timeB
    character(50) :: CIFile
    integer :: CIUnit
    real(8) :: CIenergy
    integer, allocatable :: indexArray(:),auxIndexArray(:)
    real(8), allocatable :: energyArray(:),auxEnergyArray(:)
    integer :: starting, ending, step, maxConfigurations
    character(50) :: fileNumberA, fileNumberB
    integer :: cmax
    integer :: maxStackSize, i, ia, ib, ssize, ci,cj, size1, size2
    integer :: nblocks

    size1 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=1)
    size2 = size(ConfigurationInteraction_instance%configurations(1)%occupations, dim=2) 

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 

    allocate (ConfigurationInteraction_instance%auxconfs (size1,size2, ConfigurationInteraction_instance%numberOfConfigurations ))

    do a=1, ConfigurationInteraction_instance%numberOfConfigurations
        ConfigurationInteraction_instance%auxconfs(:,:,a) = ConfigurationInteraction_instance%configurations(a)%occupations
    end do


    timeA = omp_get_wtime()

    CIFile = "lowdin.ci"
    CIUnit = 4

#ifdef intel
    open(unit=CIUnit, file=trim(CIFile), action = "write", form="unformatted", BUFFERED="YES")
#else
    open(unit=CIUnit, file=trim(CIFile), action = "write", form="unformatted")
#endif
    
    print *, "  OMP Number of threads: " , omp_get_max_threads()
    nproc = omp_get_max_threads()

    !call omp_set_num_threads(omp_get_max_threads())
    !call omp_set_num_threads(nproc)

    !if (allocated(cmax)) deallocate(cmax)
    !allocate(cmax(nproc))
    cmax = 0

    maxConfigurations = ConfigurationInteraction_instance%numberOfConfigurations
    if (allocated(indexArray )) deallocate(indexArray)
    allocate (indexArray(maxConfigurations))
    indexArray = 0
    if (allocated(energyArray )) deallocate(energyArray)
    allocate (energyArray(maxConfigurations))
    energyArray = 0

    do a=1, ConfigurationInteraction_instance%numberOfConfigurations

      !indexArray = 0
      energyArray = 0
      c = 0

!$omp parallel & 
!$omp& private(b,CIenergy),&
!$omp& shared(indexArray,energyArray, HartreeFock_instance),&
!$omp& shared(ConfigurationInteraction_instance) reduction (+:c)
!$omp do schedule(guided)
        do b= a, ConfigurationInteraction_instance%numberOfConfigurations
!          CIenergy = ConfigurationInteraction_calculateCoupling( a, b, size1, size2 )

          if ( abs(CIenergy) > 1E-9 ) then
            c = c +1   
            !indexArray(b) = b
            energyArray(b) = CIenergy
          end if
      end do
!$omp end do nowait
!$omp end parallel

       
      cmax = cmax + c

      write(CIUnit) c
      write(CIUnit) a

      allocate (auxEnergyArray(c))
      allocate (auxIndexArray(c))

      cj = 0
      do ci = a, ConfigurationInteraction_instance%numberOfConfigurations
        !if ( indexArray(ci) > 0 ) then
        if ( abs(energyArray(ci)) > 1E-9 ) then
          cj = cj + 1
          auxIndexArray(cj) =(ci)
          auxEnergyArray(cj) = energyArray(ci)
        end if
      end do
      nblocks = ceiling(real(c) / real(maxStackSize) )

      do i = 1, nblocks - 1
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        write(CIUnit) auxIndexArray(ia:ib)
      end do

      ia = maxStackSize * (nblocks - 1) + 1
      write(CIUnit) auxIndexArray(ia:c)

      deallocate(auxIndexArray)

      do i = 1, nblocks - 1
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        write(CIUnit) auxEnergyArray(ia:ib)
      end do

      ia = maxStackSize * (nblocks - 1) + 1
      write(CIUnit) auxEnergyArray(ia:c)

      deallocate (auxEnergyArray)

    end do

    write(CIUnit) -1

    close(CIUnit)

    deallocate(indexArray)
    deallocate(energyArray)
    
    timeB = omp_get_wtime()
    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Building CI matrix : ", timeB - timeA ," (s)"
    print *, "Nonzero elements", cmax

  end subroutine ConfigurationInteraction_buildAndSaveCIMatrix

  function ConfigurationInteraction_calculateEnergyZero( this ) result (auxCIenergy)
    implicit none

    integer(8) :: this(:)
    integer(8) :: a, b
    integer :: i,j,s
    integer :: l,k,z,kk,ll
    integer :: factor
    integer(2) :: numberOfDiffOrbitals
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    real(8) :: auxCIenergy

    auxCIenergy = 0.0_8

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = this(i)
      do kk=1, ConfigurationInteraction_instance%occupationNumber( i )  !! 1 is from a and 2 from b

        k = ConfigurationInteraction_instance%strings(i)%values(kk,a)

        !One particle terms
        auxCIenergy = auxCIenergy + &
                    ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( k, k )

        !Two particles, same specie
        auxIndex1 = ConfigurationInteraction_instance%twoIndexArray(i)%values(k,k)

        do ll = kk + 1, ConfigurationInteraction_instance%occupationNumber( i )  !! 1 is from a and 2 from b

          l = ConfigurationInteraction_instance%strings(i)%values(ll,a)
          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(i)%values(l,l)
          auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values(auxIndex1,auxIndex2) 

          !Coulomb
          auxCIenergy = auxCIenergy + &
              ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

          !Exchange, depends on spin

          auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                        ConfigurationInteraction_instance%twoIndexArray(i)%values(k,l), &
                        ConfigurationInteraction_instance%twoIndexArray(i)%values(l,k) )

          auxCIenergy = auxCIenergy + &
                  MolecularSystem_instance%species(i)%kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
        end do

        !!Two particles, different species
        do j = i + 1, MolecularSystem_instance%numberOfQuantumSpecies
          b = this(j)
          auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

          do ll = 1, ConfigurationInteraction_instance%occupationNumber( j ) !! 1 is from a and 2 from b
            l = ConfigurationInteraction_instance%strings(j)%values(ll,b)

            auxIndex2= ConfigurationInteraction_instance%twoIndexArray(j)%values(l,l)
            auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

            auxCIenergy = auxCIenergy + &!couplingEnergy
            ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

          end do

        end do

      end do
    end do

    auxCIenergy= auxCIenergy + HartreeFock_instance%puntualInteractionEnergy

  end function ConfigurationInteraction_calculateEnergyZero

  function ConfigurationInteraction_calculateEnergyOne( n, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: a, b
    integer :: i,j,s,n, nn
    integer :: l,k,z,kk,ll
    integer :: factor
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    integer :: diffOrb(2), otherdiffOrb(2) !! to avoid confusions
    real(8) :: auxCIenergy
    integer :: auxOcc

    auxCIenergy = 0.0_8

    factor = 1

    !! copy a
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(i)

      ConfigurationInteraction_instance%auxstring(n,i)%values(:) = ConfigurationInteraction_instance%strings(i)%values(:,a)
    end do

    !! set at maximum coincidence

    do s = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(s)
      b = thisB(s)

      do i = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s) !b
        do j = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s) !a
          if ( ConfigurationInteraction_instance%auxstring(n,s)%values(j) == &
             ConfigurationInteraction_instance%strings(s)%values(i,b) ) then

            auxOcc = ConfigurationInteraction_instance%auxstring(n,s)%values(i) 
            ConfigurationInteraction_instance%auxstring(n,s)%values(i) = ConfigurationInteraction_instance%strings(s)%values(i,b)
            ConfigurationInteraction_instance%auxstring(n,s)%values(j) = auxOcc
            if ( i /= j ) factor = -1*factor
            exit
          end if
        end do
      end do
    end do

    !! calculate
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies

      a = thisA(i)
      b = thisB(i)
      diffOrb = 0

      do kk = 1, ConfigurationInteraction_instance%occupationNumber( i) !! 1 is from a and 2 from b

        if ( ConfigurationInteraction_instance%auxstring(n,i)%values(kk) .ne. &
                 ConfigurationInteraction_instance%strings(i)%values(kk,b) ) then
          diffOrb(1) = ConfigurationInteraction_instance%auxstring(n,i)%values(kk)
          diffOrb(2) = ConfigurationInteraction_instance%strings(i)%values(kk,b)
          exit                   
        end if

      end do
      if (  diffOrb(2) > 0 ) then 

        !One particle terms
        auxCIenergy= auxCIenergy +  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values( &
                           diffOrb(1), diffOrb(2) )

        auxIndex1= ConfigurationInteraction_instance%twoIndexArray(i)%values( & 
                         diffOrb(1), diffOrb(2))

        do ll = 1, ConfigurationInteraction_instance%occupationNumber( i ) !! 1 is from a and 2 from b

          if ( ConfigurationInteraction_instance%auxstring(n,i)%values(ll) .eq. &
                 ConfigurationInteraction_instance%strings(i)%values(ll,b) ) then

            l = ConfigurationInteraction_instance%auxstring(n,i)%values(ll) !! or b

            auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(i)%values( l,l) 

            auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( auxIndex1, auxIndex2 )

            auxCIenergy = auxCIenergy + &
                        ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)


            auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                                ConfigurationInteraction_instance%twoIndexArray(i)%values(diffOrb(1),l), &
                                ConfigurationInteraction_instance%twoIndexArray(i)%values(l,diffOrb(2)) ) 

            auxCIenergy = auxCIenergy + &
                    MolecularSystem_instance%species(i)%kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

          end if
        end do
        if (MolecularSystem_instance%numberOfQuantumSpecies .gt. 1 ) then !.and. spin(1) .eq. spin(2) ) then
          do j=1, MolecularSystem_instance%numberOfQuantumSpecies

            if (i .ne. j) then

              auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 

              do ll=1,  ConfigurationInteraction_instance%occupationNumber( j ) !! 1 is from a and 2 from b
                l = ConfigurationInteraction_instance%auxstring(n,j)%values(ll) !! or b?

                auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values( l,l) 
                auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                auxCIenergy = auxCIenergy + &
                      ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1) 
              end do
            end if
          end do
        end if
      end if
    end do

    auxCIenergy= auxCIenergy * factor


  end function ConfigurationInteraction_calculateEnergyOne


  function ConfigurationInteraction_calculateEnergyTwo( n, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: a, b
    integer :: i,j,s,n
    integer :: l,k,z,kk,ll
    integer :: factor
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    integer :: diffOrb(4), otherdiffOrb(4) !! to avoid confusions
    real(8) :: auxCIenergy
    integer :: auxOcc

    auxCIenergy = 0.0_8
    factor = 1

    !! copy a
    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(i)
      ConfigurationInteraction_instance%auxstring(n,i)%values(:) = ConfigurationInteraction_instance%strings(i)%values(:,a)
    end do

    !! set at maximum coincidence

    do s = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(s)
      b = thisB(s)

      do i = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s) !b
        do j = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(s) !a
          if ( ConfigurationInteraction_instance%auxstring(n,s)%values(j) == &
                 ConfigurationInteraction_instance%strings(s)%values(i,b) ) then

            auxOcc = ConfigurationInteraction_instance%auxstring(n,s)%values(i) 
            ConfigurationInteraction_instance%auxstring(n,s)%values(i) = ConfigurationInteraction_instance%strings(s)%values(i,b)
            ConfigurationInteraction_instance%auxstring(n,s)%values(j) = auxOcc
            if ( i /= j ) factor = -1*factor
            exit
          end if
        end do
      end do
    end do

    !!calculate
    do i=1, MolecularSystem_instance%numberOfQuantumSpecies

      a = thisA(i)
      b = thisB(i)
      diffOrb = 0
      z = 1 
      do k = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(i)

        if ( ConfigurationInteraction_instance%auxstring(n,i)%values(k) .ne. &
                 ConfigurationInteraction_instance%strings(i)%values(k,b) ) then
          diffOrb(z) = ConfigurationInteraction_instance%auxstring(n,i)%values(k) 
          diffOrb(z+2) = ConfigurationInteraction_instance%strings(i)%values(k,b)  
          z = z + 1
          cycle
        end if 
      end do 
      if (  diffOrb(2) > 0 ) then

        !Coulomb
        auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                          ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                             diffOrb(1),diffOrb(3)),&
                          ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                             diffOrb(2),diffOrb(4)) )

         auxCIenergy = auxCIenergy + &
                  ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

         auxIndex = ConfigurationInteraction_instance%fourIndexArray(i)%values( &
                          ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                             diffOrb(1),diffOrb(4)),&
                          ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                             diffOrb(2),diffOrb(3)) )

         auxCIenergy = auxCIenergy + &
               MolecularSystem_instance%species(i)%kappa*ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

      end if
      !! different species
      do j = i + 1, MolecularSystem_instance%numberOfQuantumSpecies
        auxnumberOfOtherSpecieSpatialOrbitals = ConfigurationInteraction_instance%numberOfSpatialOrbitals2%values(j) 
        otherdiffOrb = 0
        a = thisA(j)
        b = thisB(j)

        do k = 1, ConfigurationInteraction_instance%numberOfOccupiedOrbitals%values(j)
          if ( ConfigurationInteraction_instance%auxstring(n,j)%values(k) .ne. &
                ConfigurationInteraction_instance%strings(j)%values(k,b) ) then
            otherdiffOrb(1) = ConfigurationInteraction_instance%auxstring(n,j)%values(k)
            otherdiffOrb(3) = ConfigurationInteraction_instance%strings(j)%values(k,b)
            exit 
          end if 

        end do 

        if ( diffOrb(3) .gt. 0 .and. otherdiffOrb(3) .gt. 0 ) then
          auxIndex1 = ConfigurationInteraction_instance%twoIndexArray(i)%values(&
                                   diffOrb(1),diffOrb(3) )
          auxIndex2 = ConfigurationInteraction_instance%twoIndexArray(j)%values(&
                                   otherdiffOrb(1),otherdiffOrb(3) )
          auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

          auxCIenergy = auxCIenergy + &
                        ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

        end if
      end do
    end do

    auxCIenergy= auxCIenergy * factor

  end function ConfigurationInteraction_calculateEnergyTwo

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine ConfigurationInteraction_getTransformedIntegrals()
    implicit none

    integer :: numberOfSpecies
    integer :: i,j,m,n,mu,nu,a,b
    integer(8) :: c
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecie
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients
    real(8) :: charge
    real(8) :: otherSpecieCharge

    integer :: ssize1, ssize2
    type(Matrix) :: externalPotential

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(ConfigurationInteraction_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))

    allocate(ConfigurationInteraction_instance%twoIndexArray(numberOfSpecies))
    allocate(ConfigurationInteraction_instance%fourIndexArray(numberOfSpecies))

!    print *,""
!    print *,"BEGIN INTEGRALS TRANFORMATION:"
!    print *,"========================================"
!    print *,""
!    print *,"--------------------------------------------------"
!    print *,"    Algorithm Four-index integral tranformation"
!    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!    print *,"  Computer Physics Communications, 2005, 166, 58-65"
!    print *,"--------------------------------------------------"
!    print *,""
!
!    call TransformIntegrals_constructor( repulsionTransformer )

    do i=1, numberOfSpecies
      nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
      specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
      ocupationNumber = MolecularSystem_getOcupationNumber( i )
      numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
      charge=MolecularSystem_getCharge(i)

!        write (6,"(T10,A)")"ONE PARTICLE INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
      call Matrix_constructor (ConfigurationInteraction_instance%twoCenterIntegrals(i), &
        int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

      call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

      !! Open file for wavefunction

      wfnFile = "lowdin.wfn"
      wfnUnit = 20

      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      arguments(2) = MolecularSystem_getNameOfSpecie(i)
      arguments(1) = "COEFFICIENTS"

      coefficients = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

      arguments(1) = "HCORE"

      hcoreMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        !arguments(1) = "FOCK"
        !ConfigurationInteraction_instance%FockMatrix(i) = &
        !  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
        !  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        !arguments(1) = "ORBITALS"
        !call Vector_getFromFile( elementsNum = numberOfContractions, &
        !  unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
        !  output =ConfigurationInteraction_instance%energyofmolecularorbitals(i) )     

        !do m=1,numberOfContractions
        !   ConfigurationInteraction_instance%fockMatrix(i)%values(m,m) = &
        !        ConfigurationInteraction_instance%energyofmolecularorbitals(i)%values(m) 
        !end do

        if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          arguments(1) = "EXTERNAL_POTENTIAL"

          externalPotential = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

          hcoreMatrix%values = hcoreMatrix%values + externalPotential%values
        end if
        !print *, "fock matrix for species", i
        !call matrix_show ( ConfigurationInteraction_instance%fockMatrix(i) )

        do m=1,numberOfContractions
          do n=m, numberOfContractions
             do mu=1, numberOfContractions
                do nu=1, numberOfContractions
                    ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) = &
                        ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) + &
                        coefficients%values(mu,m)* &
                        coefficients%values(nu,n)* &
                        hcoreMatrix%values(mu,nu)
            end do
          end do
        end do
      end do

!! Not implemented yet
!!       if( WaveFunction_HF_instance( specieID )%isThereExternalPotential ) then
!!          do m=1,numberOfContractions
!!             do n=m, numberOfContractions
!!                do mu=1, numberOfContractions
!!                   do nu=1, numberOfContractions
!!                      ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) = &
!!                           ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n) + &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(mu,m)* &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(nu,n) * &
!!                           WaveFunction_HF_instance( specieID )%ExternalPotentialMatrix%values(mu,nu)
!!                   end do
!!                end do
!!             end do
!!          end do
!!       end if

      do m = 1,numberOfContractions
        do n = m, numberOfContractions
          ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(n,m)=&
                  ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n)
        end do
      end do

      call Matrix_constructorInteger8(ConfigurationInteraction_instance%twoIndexArray(i), &
                          int( numberOfContractions,8), int( numberOfContractions,8) , 0_8 )

      c = 0
      do a=1,numberOfContractions
        do b = a, numberOfContractions
          c = c + 1
          ConfigurationInteraction_instance%twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          ConfigurationInteraction_instance%twoIndexArray(i)%values(b,a) = ConfigurationInteraction_instance%twoIndexArray(i)%values(a,b)
        end do 
      end do


      ssize1 = MolecularSystem_getTotalNumberOfContractions( i )
      ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2

      call Matrix_constructorInteger8(ConfigurationInteraction_instance%fourIndexArray(i), &
                          int( ssize1,8), int( ssize1,8) , 0_8 )
      c = 0
      do a = 1, ssize1
        do b = a, ssize1
          c = c + 1
          ConfigurationInteraction_instance%fourIndexArray(i)%values(a,b) = c! IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          ConfigurationInteraction_instance%fourIndexArray(i)%values(b,a) = &
               ConfigurationInteraction_instance%fourIndexArray(i)%values(a,b)
         end do 
       end do


       call ReadTransformedIntegrals_readOneSpecies( specieID, ConfigurationInteraction_instance%fourCenterIntegrals(i,i)   )
       ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values = &
           ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values * charge * charge

       if ( numberOfSpecies > 1 ) then
         do j = 1 , numberOfSpecies
           if ( i .ne. j) then
             nameOfOtherSpecie = trim(  MolecularSystem_getNameOfSpecie( j ) )
             otherSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
             ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
             numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
             otherSpecieCharge = MolecularSystem_getCharge(j)

             call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         ConfigurationInteraction_instance%fourCenterIntegrals(i,j) )
             ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values = &
               ConfigurationInteraction_instance%fourCenterIntegrals(i,j)%values * charge * otherSpeciecharge


           end if
         end do
       end if
     end do
     close (wfnUnit)
     call Matrix_destructor (hcoreMatrix)

  end subroutine ConfigurationInteraction_getTransformedIntegrals

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
!  subroutine ConfigurationInteraction_printTransformedIntegralsToFile()
!    implicit none
!
!!    type(TransformIntegrals) :: repulsionTransformer
!    integer :: numberOfSpecies
!    integer :: i,j,m,n,mu,nu
!    integer :: a,b,r,s,u, auxIndex
!    integer :: z
!    integer :: stats, recNum
!    character(10) :: nameOfSpecie, auxNameOfSpecie
!    character(10) :: nameOfOtherSpecie
!    integer :: ocupationNumber
!    integer :: ocupationNumberOfOtherSpecie
!    integer :: numberOfContractions
!    integer :: numberOfContractionsOfOtherSpecie
!    type(Matrix) :: auxMatrix
!    type(Matrix) :: molecularCouplingMatrix
!    type(Matrix) :: molecularExtPotentialMatrix
!
!    integer :: spin
!
!    real(8) :: totalCoupEnergy
!    real(8) :: fixedPotEnergy
!    real(8) :: fixedIntEnergy
!    real(8) :: KineticEnergy
!    real(8) :: RepulsionEnergy
!    real(8) :: couplingEnergy


!    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
!
!    print *,""
!    print *,"BEGIN INTEGRALS TRANFORMATION:"
!    print *,"========================================"
!    print *,""
!    print *,"--------------------------------------------------"
!    print *,"    Algorithm Four-index integral tranformation"
!    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!    print *,"  Computer Physics Communications, 2005, 166, 58-65"
!    print *,"--------------------------------------------------"
!    print *,""
!
!    totalCoupEnergy = 0.0_8
!    fixedPotEnergy = 0.0_8
!    fixedIntEnergy = 0.0_8
!    KineticEnergy = 0.0_8
!    RepulsionEnergy = 0.0_8
!    couplingEnergy = 0.0_8
!    spin = 0
!
!    do i=1, numberOfSpecies
!        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
!        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
!        spin = MolecularSystem_getMultiplicity(i) - 1
!
!        if(trim(nameOfSpecie) /= "E-BETA" ) then
!
!           if(trim(nameOfSpecie) /= "U-" ) then 
!
!              open(unit=35, file="FCIDUMP-"//trim(nameOfSpecie)//".com", form="formatted", status="replace")
!
!              write(35,"(A)")"gprint basis"
!              write(35,"(A)")"memory 1000 M"
!              write(35,"(A)")"cartesian"
!              write(35,"(A)")"gthresh twoint=1e-12 prefac=1e-14 energy=1e-10 edens=1e-10 zero=1e-12"
!              write(35,"(A)")"basis={"
!              call ConfigurationInteraction_printBasisSetToFile(35)
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"symmetry nosym"
!              write(35,"(A)")"angstrom"
!              write(35,"(A)")"geometry={"
!              call ConfigurationInteraction_printGeometryToFile(35)
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"import 21500.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"jcoup")
!              write(35,"(A)")"import 21510.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"icoup")
!              write(35,"(A)")"import 21520.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"kin")
!              write(35,"(A)")"import 21530.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"coeff")
!
!              if(trim(nameOfSpecie) == "E-ALPHA") then
!
!                 write(35,"(A)")"import 21550.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//"E-BETA"//"."//"coeff")
!
!              end if
!
!              write(35,"(A)")"import 21540.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"dens")
!              !write(35,"(A)")"import 21560.2 "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"pot")
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load Jcoup, SQUARE 21500.2"
!              write(35,"(A)")"load Icoup, SQUARE 21510.2"
!              write(35,"(A)")"load K, SQUARE 21520.2"
!              !write(35,"(A)")"load Pot, SQUARE 21560.2"
!              write(35,"(A)")"add H01, K Icoup Jcoup"! Pot"
!              write(35,"(A)")"save H01, 21511.2 H0"
!              write(35,"(A)")"}"
!
!              if(trim(nameOfSpecie) == "E-ALPHA") then
!                 write(35,"(A)")"{matrop"
!                 write(35,"(A)")"load Ca, SQUARE 21530.2"
!                 write(35,"(A)")"load Cb, SQUARE 21550.2"               
!                 write(35,"(A)")"save Ca, 2100.1 ORBITALS alpha"
!                 write(35,"(A)")"save Cb, 2100.1 ORBITALS beta"
!                 write(35,"(A)")"}"
!              else
!                 write(35,"(A)")"{matrop"
!                 write(35,"(A)")"load C, SQUARE 21530.2"
!                 write(35,"(A)")"save C, 2100.1 ORBITALS"
!                 write(35,"(A)")"}"
!              end if
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load D, SQUARE 21540.2"
!              write(35,"(A)")"save D, 21400.1 DENSITY"
!              write(35,"(A)")"}"
!
!
!              !            write(35,"(A,I3,A,I3,A,I3,A1)")"$FCI NORB=",numberOfContractions, ",NELEC=", MolecularSystem_getNumberOfParticles(i)-spin, ", MS2=", spin,","
!              !
!              !            write(35,"(A)",advance="no") "ORBSYM="
!              !            do z=1, numberOfContractions
!              !                write(35,"(I1,A1)",advance="no") 1,","
!              !            end do
!              !            write(35,"(A)") ""
!              !
!              !            write(35, "(A,I3,A,I9)") "ISYM=",1, ",MEMORY=", 200000000
!              !
!              !            write(35, "(A)") "$"
!              !
!              !            print *, "FOUR CENTER INTEGRALS FOR SPECIE: ", trim(nameOfSpecie)
!              !
!              !            recNum = 0
!              !            do a = 1, numberOfContractions
!              !                n = a
!              !                do b=a, numberOfContractions
!              !                    u = b
!              !                    do r = n, numberOfContractions
!              !                        do s = u, numberOfContractions
!              !
!              !                            auxIndex = IndexMap_tensorR4ToVector( a, b, r, s, numberOfContractions )
!              !                            write(35,"(F20.10,4I3)") ConfigurationInteraction_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1), a, b, r, s
!              !
!              !                        end do
!              !                        u=r+1
!              !                    end do
!              !                end do
!              !            end do
!              !
!              !
!              !            print *, "TWO CENTER TRANSFORMED INTEGRALS FOR SPECIE: ", trim(nameOfSpecie)
!              !
!              !            do m=1,numberOfContractions
!              !                do n=1, m
!              !                    write(35,"(F20.10,4I3)") ConfigurationInteraction_instance%twoCenterIntegrals(i)%values(m,n), m, n, 0, 0
!              !                end do
!              !            end do
!
!              !!Calculating the core energy....
!
!
!
!              totalCoupEnergy = MolecularSystem_instance%totalCouplingEnergy
!              fixedPotEnergy = MolecularSystem_instance%puntualInteractionEnergy
!
!              do j = 1, numberOfSpecies
!
!                 auxNameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
!
!                 if(trim(auxNameOfSpecie) == "E-ALPHA" .or.  trim(auxNameOfSpecie) == "E-BETA" .or.  trim(auxNameOfSpecie) == "e-") cycle
!
!                 fixedIntEnergy = fixedIntEnergy + MolecularSystem_instance%quantumPuntualInteractionEnergy(j)
!                 KineticEnergy = KineticEnergy + MolecularSystem_instance%kineticEnergy(j)
!                 RepulsionEnergy = RepulsionEnergy + MolecularSystem_instance%repulsionEnergy(j)
!                 couplingEnergy = couplingEnergy + MolecularSystem_instance%couplingEnergy(j)
!
!              end do
!
!              !!COMO SEA QUE SE META LA ENERGIA DE CORE
!              !write(35,"(F20.10,4I3)") (couplingEnergy-totalCoupEnergy+fixedPotEnergy+fixedIntEnergy+KineticEnergy+RepulsionEnergy), 0, 0, 0, 0
!              
!              print*, "COREENERGY ", (couplingEnergy-totalCoupEnergy+fixedPotEnergy+fixedIntEnergy+KineticEnergy+RepulsionEnergy)
!
!              write(35,"(A)")"{hf"
!              write(35,"(A)")"maxit 250"
!              write(35,"(A10,I2,A1,A6,I2,A1,A6,I3)")"wf spin=", spin, ",", "charge=",0, ",", "elec=", MolecularSystem_getNumberOfParticles(i)-spin
!              write(35,"(A)")"start 2100.1"
!              write(35,"(A)")"}"
!
!
!              write(35,"(A)")"{fci"
!              write(35,"(A)")"maxit 250"
!              write(35,"(A)")"dm 21400.1, IGNORE_ERROR"
!              write(35,"(A)")"orbit 2100.1, IGNORE_ERROR"
!              write(35,"(A10,I2,A1,A6,I2,A1,A6,I3)")"wf spin=", spin, ",", "charge=",0, ",", "elec=", MolecularSystem_getNumberOfParticles(i)-spin
!              !            write(35,"(A)")"print, orbital=2 integral = 2"
!              !            write(35,"(A)")"CORE"
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"{matrop"
!              write(35,"(A)")"load D, DEN, 21400.1"
!              !	    write(35,"(A)")"print D"
!              write(35,"(A)")"natorb Norb, D"
!              write(35,"(A)")"save Norb, 21570.2"
!              !	    write(35,"(A)")"print Norb"
!              write(35,"(A)")"}"
!
!              write(35,"(A)")"put molden "//trim(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//"molden")//"; orb, 21570.2"
!
!              close(35)
!
!              print*, ""
!
!              stats = system("molpro "//"FCIDUMP-"//trim(nameOfSpecie)//".com ")
!              stats = system("cat "//"FCIDUMP-"//trim(nameOfSpecie)//".out ")
!
!              print*, ""
!
!              print *,"END"
!              
!           end if
!
!         end if
!
!    end do
    
!  end subroutine ConfigurationInteraction_printTransformedIntegralsToFile

!  subroutine ConfigurationInteraction_printGeometryToFile(unit)
!    implicit none
!    integer :: unit
!
!    integer :: i
!    integer :: from, to
!    real(8) :: origin(3)
!    character(50) :: auxString
!
!    
!    do i = 1, MolecularSystem_getTotalNumberOfParticles()
!       
!       origin = MolecularSystem_getOrigin( iterator = i ) * AMSTRONG
!       auxString = trim( MolecularSystem_getNickName( iterator = i ) )
!       
!       if( String_findSubstring( trim( auxString ), "e-") == 1 ) then
!          if( String_findSubstring( trim( auxString ), "BETA") > 1 ) then
!             cycle
!          end if
!            
!          from =String_findSubstring( trim(auxString), "[")
!          to = String_findSubstring( trim(auxString), "]")
!          auxString = auxString(from+1:to-1)
!          
!       else if( String_findSubstring( trim( auxString ), "_") /= 0 ) then
!          cycle
!       end if
!         
!         
!       write (unit,"(A10,3F20.10)") trim( auxString ), origin(1), origin(2), origin(3)
!       
!    end do

!  end subroutine ConfigurationInteraction_printGeometryToFile


!  subroutine ConfigurationInteraction_printBasisSetToFile(unit)
!    implicit none
!
!    integer :: unit
!
!    integer :: i, j
!    character(16) :: auxString
!
!
!    do i =1, MolecularSystem_instance%numberOfQuantumSpecies
!       
!       auxString=trim( Map_getKey( MolecularSystem_instance%speciesID, iterator=i ) )
!       
!       if( String_findSubstring( trim(auxString), "e-") == 1 ) then
!          
!          if( String_findSubstring( trim(auxString), "BETA") > 1 ) then
!             
!             cycle
!             
!          end if
!          
!          
!       end if
!       
!       if(trim(auxString)=="U-") cycle
!
!       do j =1, size(MolecularSystem_instance%particlesPtr)
!
!          if (    trim(MolecularSystem_instance%particlesPtr(j)%symbol) == trim( Map_getKey( MolecularSystem_instance%speciesID, iterator=i ) ) &
!               .and. MolecularSystem_instance%particlesPtr(j)%isQuantum ) then
!             
!             call BasisSet_showInMolproForm( MolecularSystem_instance%particlesPtr(j)%basis, trim(MolecularSystem_instance%particlesPtr(j)%nickname), unit=unit )
!             
!          end if
!          
!       end do
!       
!    end do
    
!  end subroutine ConfigurationInteraction_printBasisSetToFile


  !**
  ! @ Retorna la energia final com correccion Moller-Plesset de orrden dado
  !**
  function ConfigurationInteraction_getTotalEnergy() result(output)
    implicit none
    real(8) :: output

    output = ConfigurationInteraction_instance%totalEnergy

  end function ConfigurationInteraction_getTotalEnergy


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine ConfigurationInteraction_exception( typeMessage, description, debugDescription)
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

  end subroutine ConfigurationInteraction_exception

  subroutine ConfigurationInteraction_saveEigenVector () 
    implicit none
    character(50) :: nameFile
    integer :: unitFile
    integer(8) :: i, ia
    integer :: ib, nonzero
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
    nameFile = "lowdin.civec"
    unitFile = 20

    nonzero = 0
    do i = 1, ConfigurationInteraction_instance%numberOfConfigurations
      if ( abs(ConfigurationInteraction_instance%eigenVectors%values(i,1) ) >= 1E-12 ) nonzero = nonzero + 1
    end do 

    write (*,*) "nonzero", nonzero

    allocate(auxArray(nonzero))
    allocate(auxIndexArray(nonzero))

    ia = 0
    do i = 1, ConfigurationInteraction_instance%numberOfConfigurations
      if ( abs(ConfigurationInteraction_instance%eigenVectors%values(i,1) ) >= 1E-12 ) then 
        ia = ia + 1
        auxIndexArray(ia) = i 
        auxArray(ia) = ConfigurationInteraction_instance%eigenVectors%values(i,1) 
      end if
    end do 

    open(unit=unitFile, file=trim(nameFile), status="replace", form="unformatted")

    write(unitFile) ConfigurationInteraction_instance%eigenValues%values(1)
    write(unitFile) nonzero

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxIndexArray(ia:ib)
    end do
    deallocate(auxIndexArray)

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxArray(ia:ib)
    end do
    deallocate(auxArray)

    close(unitFile)

  end subroutine ConfigurationInteraction_saveEigenVector

  subroutine ConfigurationInteraction_loadEigenVector (eigenValues,eigenVectors) 
    implicit none
    type(Vector8) :: eigenValues
    type(Matrix) :: eigenVectors
    character(50) :: nameFile
    integer :: unitFile
    integer :: i, ia, ib, nonzero
    real(8) :: eigenValue
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
 

    nameFile = "lowdin.civec"
    unitFile = 20


    open(unit=unitFile, file=trim(nameFile), status="old", action="read", form="unformatted")

    readvectors : do
      read (unitFile) eigenValue
      read (unitFile) nonzero
      write (*,*) "eigenValue", eigenValue
      write (*,*) "nonzero", nonzero

      allocate (auxIndexArray(nonzero))
      auxIndexArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxIndexArray(ia:ib)
      end do

      allocate (auxArray(nonzero))
      auxArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxArray(ia:ib)
      end do
      exit readvectors
    end do readvectors

    eigenValues%values(1) = eigenValue
    do i = 1, nonzero
      eigenVectors%values(auxIndexArray(i),1) = auxArray(i)
    end do

    deallocate (auxIndexArray )
    deallocate (auxArray )


    close(unitFile)

  end subroutine ConfigurationInteraction_loadEigenVector

  subroutine av ( nx, v, w)
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    integer(8) nx
    real(8) v(nx)
    real(8) w(nx)
    character(50) :: CIFile
    integer :: CIUnit
    integer, allocatable :: jj(:)
    real(8), allocatable :: CIEnergy(:)
    integer :: nonzero,ii, kk
    integer :: maxStackSize, i, ia, ib

    CIFile = "lowdin.ci"
    CIUnit = 20
    nonzero = 0
    maxStackSize = CONTROL_instance%CI_STACK_SIZE 

    w = 0
#ifdef intel
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted", BUFFERED="YES")
#else
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted")
#endif

    readmatrix : do
      read (CIUnit) nonzero
      if (nonzero > 0 ) then

        read (CIUnit) ii

        if ( allocated(jj)) deallocate (jj)
        allocate (jj(nonzero))
        jj = 0

        if ( allocated(CIEnergy)) deallocate (CIEnergy)
        allocate (CIEnergy(nonzero))
        CIEnergy = 0

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) jj(ia:ib)
    
        end do

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) CIEnergy(ia:ib)
    
        end do

        w(ii) = w(ii) + CIEnergy(1)*v(jj(1)) !! disk
        do kk = 2, nonzero
          !w(ii) = w(ii) + ConfigurationInteraction_calculateCIenergy(ii,jj(kk))*v(jj(kk))  !! direct
          w(ii) = w(ii) + CIEnergy(kk)*v(jj(kk)) !! disk
          w(jj(kk)) = w(jj(kk)) + CIEnergy(kk)*v(ii) !! disk
        end do

      else if ( nonzero == -1 ) then
        exit readmatrix
      end if
    end do readmatrix

!! memory
!    do i = 1, nx
!        w(:) = w(:) + ConfigurationInteraction_instance%hamiltonianMatrix%values(:,i)*v(i)
!    end do 

     close(CIUnit)

    return
  end subroutine av


  subroutine ConfigurationInteraction_jadamiluInterface(n,  maxeig, eigenValues, eigenVectors)
    implicit none
    integer(8) :: maxnev
    real(8) :: CIenergy
    integer(8) :: nproc
    type(Vector8), intent(inout) :: eigenValues
    type(Matrix), intent(inout) :: eigenVectors

!   N: size of the problem
!   MAXEIG: max. number of wanteg eig (NEIG<=MAXEIG)
!   MAXSP: max. value of MADSPACE
    integer(8) :: n, maxeig, MAXSP
    integer(8) :: LX
    real(8), allocatable :: EIGS(:), RES(:), X(:)!, D(:)
!   arguments to pass to the routines
    integer(8) :: NEIG, MADSPACE, ISEARCH, NINIT
    integer(8) :: ICNTL(5)
    integer(8) :: ITER, IPRINT, INFO
    real(8) :: SIGMA, TOL, GAP, MEM, DROPTOL, SHIFT
    integer(8) :: NDX1, NDX2, NDX3
    integer(8) :: IJOB!   some local variables
    integer(8) :: auxSize
    integer(4) :: size1,size2
    integer(8) :: I,J,K,ii,jj,jjj
    integer(4) :: iiter
    logical :: fullMatrix
    
    maxsp = CONTROL_instance%CI_MADSPACE
    !!if ( CONTROL_instance%CI_JACOBI ) then

    LX = N*(3*MAXSP+MAXEIG+1)+4*MAXSP*MAXSP

    if ( allocated ( eigs ) ) deallocate ( eigs )
    allocate ( eigs ( maxeig ) )
    eigs = 0.0_8
    if ( allocated ( res ) ) deallocate ( res )
    allocate ( res ( maxeig ) )
    res = 0.0_8
    if ( allocated ( x ) ) deallocate ( x )
    allocate ( x ( lx ) )
    x = 0.0_8


!    set input variables
!    the matrix is already in the required format

     IPRINT = -6 !     standard report on standard output
     ISEARCH = 1 !    we want the smallest eigenvalues
     NEIG = maxeig !    number of wanted eigenvalues
     !NINIT = 0 !    no initial approximate eigenvectors
     NINIT = NEIG !    initial approximate eigenvectors
     MADSPACE = maxsp !    desired size of the search space
     ITER = 1000*NEIG !    maximum number of iteration steps
     TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual

     NDX1 = 0
     NDX2 = 0
     MEM = 0

!    additional parameters set to default
     ICNTL(1)=0
     ICNTL(2)=0
     ICNTL(3)=0
     ICNTL(4)=0
     ICNTL(5)=1

     IJOB=0

     ! set initial eigenpairs
     if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
       print *, "Loading the eigenvector to the initial guess"
       do j = 1, n 
         X(j) = eigenVectors%values(j,1)
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = eigenValues%values(i)
       end do
     else
       jj = 0 
       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         jj = (i - 1) * n 
         do j = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
          X(jj + ConfigurationInteraction_instance%auxIndexCIMatrix%values(j)) = ConfigurationInteraction_instance%initialEigenVectors%values(j,i)
         end do
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = ConfigurationInteraction_instance%initialEigenValues%values(i)
       end do
     end if

     DROPTOL = 0

     SIGMA = EIGS(1)
     gap = 0 
     SHIFT = EIGS(1)

     do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
       write(6,"(T2,A5,I4,2X,A10,F20.10,2X,A11,F20.10)") "State", i, "Eigenvalue", eigs( i ), "Eigenvector", x((i-1)*n + i)
     end do

     iiter = 0
  
!10     CALL DPJDREVCOM( N, A, JA, IA,EIGS, RES, X, LX, NEIG, &
!                        SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL, &
!                        SHIFT, DROPTOL, MEM, ICNTL, &
!                        IJOB, NDX1, NDX2, IPRINT, INFO, GAP)
10   CALL DPJDREVCOM( N, ConfigurationInteraction_instance%diagonalHamiltonianMatrix%values ,-1_8,-1_8,EIGS, RES, X, LX, NEIG, &
                        SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL, &
                        SHIFT, DROPTOL, MEM, ICNTL, &
                        IJOB, NDX1, NDX2, IPRINT, INFO, GAP)
      if (CONTROL_instance%CI_JACOBI ) then
        fullMatrix = .false.
      else 
        fullMatrix = .true.
      end if
!!    your private matrix-vector multiplication
  
      iiter = iiter +1
      IF (IJOB.EQ.1) THEN
        if ( CONTROL_instance%CI_BUILD_FULL_MATRIX ) then
          call av ( n, x(ndx1), x(ndx2))
        else 
          call matvec2 ( N, X(NDX1), X(NDX2), iiter)
        end if

        GOTO 10
      END IF
  
      !! saving the eigenvalues
      eigenValues%values = EIGS

      !! saving the eigenvectors
      k = 0
      do j = 1, maxeig
         do i = 1, N
          k = k + 1
          eigenVectors%values(i,j) = X(k)
        end do
      end do

!    release internal memory and discard preconditioner
     CALL PJDCLEANUP
     if ( allocated ( x ) ) deallocate ( x )

  end subroutine ConfigurationInteraction_jadamiluInterface

  subroutine matvec2 ( nx, v, w, iter)
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    integer(8) nx
    real(8) v(nx)
    real(8) w(nx)
    real(8) :: CIEnergy
    integer(8) :: nonzero
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: wi
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    !integer(8), allocatable :: indexArray(:)
    logical :: fullMatrix
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b,c
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    integer(8), allocatable :: cc(:) !! ncore
    integer(8), allocatable :: indexConf(:,:) !! ncore, species
    integer(8), allocatable :: auxindexConf(:,:) !! ncore, species
    integer, allocatable :: cilevel(:,:), auxcilevel(:,:)

    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()


    allocate( cc ( nproc ) )
    cc = 0 

    nonzero = 0
    w = 0 
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 

    do i = 1 , nx
       if ( abs(v(i) ) >= tol) nonzero = nonzero + 1
    end do
  
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate ( indexConf ( numberOfSpecies, nproc ) )
    allocate ( auxindexConf ( numberOfSpecies, nproc ) )
    allocate ( cilevel ( numberOfSpecies, nproc ) )
    allocate ( auxcilevel ( numberOfSpecies, nproc ) )

    cilevel = 0
    auxcilevel = 0
    indexConf = 0
    auxindexConf = 0
    !! call recursion
    s = 0
    c = 0
    n = 1
!$  timeA = omp_get_wtime()
    do ci = 1,  ConfigurationInteraction_instance%sizeCiOrderList 
      do nn = n, nproc
        cilevel(:,nn) =  ConfigurationInteraction_instance%ciOrderList(  ConfigurationInteraction_instance%auxciOrderList(ci), :)
      end do
      s = 0 
      auxnumberOfSpecies = ConfigurationInteraction_buildMatrixRecursion(nproc, s, indexConf, auxindexConf,cc, c, n, v, w, &
                             cilevel, auxcilevel )

    end do

    if  ( n > 1 ) then
       do nn = 1, n-1

       call ConfigurationInteraction_buildRow( nn, auxindexConf(:,nn), cc(nn), w, v(cc(nn)), auxcilevel(:,nn))
      end do
    end if
    
    ConfigurationInteraction_instance%pindexConf = 0

!$  timeB = omp_get_wtime()
    deallocate ( cilevel )
    deallocate ( auxindexConf )
    deallocate ( indexConf )
    deallocate ( cc )
!$    write(*,"(A,I2,A,E10.3,A2,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero
!   stop
    

    return

  end subroutine matvec2

end module ConfigurationInteraction_

