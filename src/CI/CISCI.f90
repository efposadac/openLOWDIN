module CISCI_
  use Exception_
  use Matrix_
  use Vector_
  use IndexMap_
  use InputCI_
  use CIcore_
  use CIJadamilu_
  use CIInitial_
  use sort_
  use omp_lib
  use CISort_
  implicit none

  type, public :: CISCI
    !! arrays for storing coefficients
    type(Vector) :: buffer_amplitudeCore
    type(Vector) :: coefficientCore
    !! auxiliary array to trace the original elements in array sorting
    type (ivector8) :: index_amplitudeCore
    !! arrays for storing CI configurations, species, orbitals, vector size
    type (IMatrix1), allocatable :: confCore(:)
    integer(1), allocatable :: confAmplitudeCore_orb(:,:)
    type (IMatrix1), allocatable :: confTarget_orb(:)
    type (IMatrix), allocatable :: confTarget_occ(:)
    !! storing the CI diagonal matrix elements for Jadamilu preconditioner
    !! eigenvalues per SCI iteration
    type (Vector), allocatable :: eigenValues(:) ! eigenvalues per SCI iteration
    real(8), allocatable :: minCoeff(:) ! eigenvalues per SCI iteration
    type(Vector) :: diagonalCore
    type(Vector) :: diagonalTarget
    !! length of SCI search vectors
    integer(8) :: coreSpaceSize
    integer(8) :: targetSpaceSize
    integer(8) :: targetSpaceSize_max
    integer(8) :: buffer_amplitudeCoreSize
    integer, allocatable :: targetSpaceSize_iter(:)
    !! PT2 perturbative energy correction to SCI
    real(8) :: PT2energy
    !! auxiliary variables to map orbitals from vector to array location
    integer :: combinedNumberOfOrbitals
    integer :: combinedNumberOfOccupiedOrbitals
    integer, allocatable :: combinedOrbitalsPositions(:,:) ! lower and upper position, species
    integer, allocatable :: combinedOccupiedOrbitalsPositions(:,:) ! lower and upper position, species
    !! auxiliary array to store the position of the target space for each omp thread 
    integer(8), allocatable :: omp_targetInterval(:,:) ! lower and upper position, n_threads
    integer(8), allocatable :: omp_target_iterator_m(:) ! n_threads
    !! ci level per species renormalization
    integer, allocatable :: maxCIexcitations(:)
    integer, allocatable :: CIorder_list(:,:)
    integer, allocatable :: CIorder_count(:)
    real(8), allocatable :: CIorder_weight(:)
    !! In case of bit-masking representation
    !!store the orbitals for each target configurations, to avoid recomputing them 
    !!type (IVector), allocatable :: targetOrb(:,:) ! species, num of target configurations % num. of orbitals
    type(ivector), allocatable :: canonicalOrder(:)
    !! heat bath sorted two particles contributions
    type(Matrix), allocatable :: heatBathDoubleExcitations(:,:)
    type(IMatrix8), allocatable :: heatBathDoubleExcitations_index(:,:)
    type(IMatrix8), allocatable :: heatBathDoubleExcitations_size(:,:)

  end type CISCI

  type(CISCI) :: CISCI_instance
  integer, parameter :: VARIATIONAL = 1
  integer, parameter :: PERTURBATIVE = 2

contains

  !! show information about the SCI calculation, memory, creators, and libraries used
  subroutine CISCI_show()
    implicit none
    integer :: spi, spj
    integer(8) :: totalSize

    !! estimate memory cost 

    !! take from input
    CISCI_instance%coreSpaceSize = CONTROL_instance%CI_SCI_CORE_SPACE
    CISCI_instance%targetSpaceSize = CONTROL_instance%CI_SCI_TARGET_SPACE !! initial size
    CISCI_instance%targetSpaceSize_max = CONTROL_instance%CI_SCI_TARGET_SPACE * ( CONTROL_instance%CI_SCI_TARGET_GROWTH_FACTOR**(CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS - 1 ))
    !! the first quarter is the real target space, the second quarter is the waiting list if the amplitude could grow more. 
    !! The last half is just a temporary space to avoid sorting a big array for a single addition
    CISCI_instance%buffer_amplitudeCoreSize = CISCI_instance%targetSpaceSize_max * CONTROL_instance%CI_SCI_BUFFER_FACTOR

    totalSize = 0_8
    do spi = 1, CIcore_instance%numberOfSpecies 
      totalSize = totalSize + &
                  ( CISCI_instance%buffer_amplitudeCoreSize * &
                                  ( 8 + 8 + 1*CIcore_instance%numberOfActiveOrbitals%values(spi) ) & ! data type for coeff, index, conf_orb
                  + CISCI_instance%coreSpaceSize * &
                                  ( 8 + 1*CIcore_instance%numberOfActiveOrbitals%values(spi)) & ! coeff, conf
                  + CISCI_instance%targetSpaceSize_max * &
                                  ( 8 + 8 + 2*8 + CIcore_instance%nproc * 8 & !! coeff, diagonal, eigenvectors, W per omp thread
                                    + 1*CIcore_instance%numberOfActiveOrbitals%values(spi) + 4*CIcore_instance%numberOfActiveOrbitals%values(spi) ) )  !! conf_orb, conf_cc

      do spj = spi, CIcore_instance%numberOfSpecies 
        totalSize = totalSize + &
                    size(CIcore_instance%fourCenterIntegrals( spi, spj )%values,1) * 8
      enddo
    enddo

    select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
      case ("DSYEVR")
        totalSize = totalSize + CISCI_instance%targetSpaceSize_max * CISCI_instance%targetSpaceSize_max
      case ("JADAMILU")
        totalSize = totalSize + CISCI_instance%targetSpaceSize_max * CONTROL_instance%CI_MADSPACE
    end select 

    write (6,*) "-----------------------------------------------------------------------"
    write (6,"(T2,A62)") "          SELECTED CONFIGURATION INTERACTION (SCI):          " 
    write (6,"(T2,A62)") "                 Adaptive Sampling CI (ASCI)                 " 
    write (6,"(T2,A62)") "                   Deterministic Algorithm                   " 
    write (6,"(T2,A62)") "                  Based on 10.1063/1.4955109                 "
    write (6,"(T2,A62)") "                        J. Charry                            "
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""
    write (6,*) "  Diagonalizer for target space hamiltonian : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
    write (6,*) "a software code for computing selected eigenvalues of "
    write (6,*) "large sparse symmetric matrices, "
    write (6,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007."
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""
    write (6,*) " Modified sorting algorithm from CENCALC quicksort code "
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) " Code available at https://github.com/dimassuarez/cencalc_quicksort  "
    write (6,*) " E. Suárez, N. Díaz, J. Méndez and D. Suárez. "
    write (6,*) " CENCALC: A Computational Tool for Conformational Entropy Calculations"
    write (6,*) " from Molecular Simulations."
    write (6,*) " J. Comput. Chem. 54, 2031. DOI: 10.1002/jcc.23350 "
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""
    write (6,"(T2,A,F14.3,A3 )") "Estimated memory needed : ", real( totalSize )/(1024**2) , " MB"
    write (6,"(T2,A,F14.3,A3 )") "                          ", real( totalSize )/(1024**3) , " GB"
    write (6,"(T2,A,I10 )") "Length of core (search) space                          :",  CISCI_instance%coreSpaceSize
    write (6,"(T2,A,I10 )") "Length of target (Full-CI subset) space per OMP thread :",  CISCI_instance%targetSpaceSize_max / CIcore_instance%nproc
    write (6,"(T2,A,I10 )") "Length of buffer (auxiliary sort) space per OMP thread :",  CISCI_instance%buffer_amplitudeCoreSize / CIcore_instance%nproc
    write (6,"(T2,A,I10 )") "Length of total target space                           :",  CISCI_instance%targetSpaceSize_max
    write (6,"(T2,A,I10 )") "Length of total buffer space                           :",  CISCI_instance%buffer_amplitudeCoreSize 
    write (6,*) "-----------------------------------------------------------------------"
    
  end subroutine CISCI_show

  !! Allocating arrays 
  subroutine CISCI_constructor( numberOfConfigurations )
    implicit none
    integer(8), intent(out) :: numberOfConfigurations
    integer :: a,b,c,aa,bb,i
    integer :: spi, spj
    real(8) :: CIenergy
    integer :: nproc, n
    integer :: numberOfSpecies
    integer :: m
    integer :: k
    integer :: pi

    numberOfSpecies = CIcore_instance%numberOfSpecies 
    numberOfConfigurations = CISCI_instance%targetSpaceSize !! initial size

    !! auxiliary variables to map orbitals from vector to array location
    CISCI_instance%combinedNumberOfOrbitals = sum(CIcore_instance%numberOfActiveOrbitals%values(:))
    allocate ( CISCI_instance%combinedOrbitalsPositions(2,numberOfSpecies) )
    m = 0
    do spi = 1, numberOfSpecies
      CISCI_instance%combinedOrbitalsPositions(1,spi) = m + 1
      CISCI_instance%combinedOrbitalsPositions(2,spi) = m + CIcore_instance%numberOfActiveOrbitals%values(spi)
      m = m + CIcore_instance%numberOfActiveOrbitals%values(spi)
    enddo 

    CISCI_instance%combinedNumberOfOccupiedOrbitals = sum(CIcore_instance%numberOfOccupiedOrbitals%values(:))
    allocate ( CISCI_instance%combinedOccupiedOrbitalsPositions(2,numberOfSpecies) )
    m = 0
    do spi = 1, numberOfSpecies
      CISCI_instance%combinedOccupiedOrbitalsPositions(1,spi) = m + 1
      CISCI_instance%combinedOccupiedOrbitalsPositions(2,spi) = m + CIcore_instance%numberOfOccupiedOrbitals%values(spi)
      m = m + CIcore_instance%numberOfOccupiedOrbitals%values(spi)
    enddo 

    !! auxiliary arrays to store the position of the target space position for each omp thread within the big arrays
    !! last batch (nproc + 1) is just dummy index when calling to sort the whole buffer
    allocate ( CISCI_instance%omp_targetInterval(2, CIcore_instance%nproc + 1 ) ) !fixed
    allocate ( CISCI_instance%omp_target_iterator_m( CIcore_instance%nproc + 1) ) !variable

    !! arrays for storing coefficients
    call Vector_constructor(CISCI_instance%buffer_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize, 8), 0.0_8)
    call Vector_constructor(CISCI_instance%coefficientCore, int(CISCI_instance%coreSpaceSize, 8), 0.0_8)

    !! auxiliary array to trace the original elements in array sorting
    call Vector_constructorInteger8 ( CISCI_instance%index_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0_8)

    !! arrays for storing CI configurations, species, orbitals, vector size
    if ( allocated ( CISCI_instance%confTarget_orb ) ) deallocate ( CISCI_instance%confTarget_orb )  
    if ( allocated ( CISCI_instance%confTarget_occ ) ) deallocate ( CISCI_instance%confTarget_occ )  

    allocate ( CISCI_instance%confCore ( numberOfSpecies ) ) 
    allocate ( CISCI_instance%confTarget_orb ( numberOfSpecies ) )
    allocate ( CISCI_instance%confTarget_occ ( numberOfSpecies ) ) 
    do spi = 1, numberOfSpecies 
      call Matrix_constructorInteger1 ( CISCI_instance%confCore(spi), int(CIcore_instance%numberOfActiveOrbitals%values(spi),8) , int(CISCI_instance%coreSpaceSize,8), -1_1 )
      call Matrix_constructorInteger1 ( CISCI_instance%confTarget_orb(spi), int(CIcore_instance%numberOfActiveOrbitals%values(spi),8) , int(CISCI_instance%targetSpaceSize,8), -1_1) !! this will be reallocated
      call Matrix_constructorInteger ( CISCI_instance%confTarget_occ(spi), int(CIcore_instance%numberOfOccupiedOrbitals%values(spi),8) , int(CISCI_instance%targetSpaceSize,8), -1_4) !! this will be reallocated
    enddo
    allocate ( CISCI_instance%confAmplitudeCore_orb (  CISCI_instance%combinedNumberOfOrbitals, CISCI_instance%buffer_amplitudeCoreSize ) ) 
    CISCI_instance%confAmplitudeCore_orb = -1_1

    call CISCI_buildCIOrderList(  CISCI_instance%maxCIexcitations, CISCI_instance%CIorder_list, CISCI_instance%CIorder_count, CISCI_instance%CIorder_weight )

    !! this was replaced by a "vectorized" array to avoid using arrays of types inside a recursive function
    !!do spi = 1, numberOfSpecies
    !!allocate ( CISCI_instance%confAmplitudeCore_orb ( numberOfSpecies ) )
      !!call Matrix_constructorInteger1 ( CISCI_instance%confAmplitudeCore_orb(spi), int(CIcore_instance%numberOfActiveOrbitals%values(spi),8) , int(CISCI_instance%buffer_amplitudeCoreSize,8), -1_1) 
    !!enddo

    !! store the orbitals for each target configurations, to avoid recomputing them
    !! this is helpful when using bit-masking approach to avoid transforming from bit to decimal multiple times
    !!allocate ( CISCI_instance%targetOrb ( numberOfSpecies, CISCI_instance%targetSpaceSize ) )
    !!do a = 1, CISCI_instance%targetSpaceSize
    !!  do spi = 1, numberOfSpecies
    !!    call Vector_constructorInteger ( CISCI_instance%targetOrb(spi,a), CIcore_instance%numberOfActiveOrbitals%values(spi), 0 )
    !!  enddo
    !!enddo

    !! storing the CI diagonal matrix elements for Jadamilu preconditioner (moved to reset buffers)
    !! call Vector_constructor ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8)
    call Vector_constructor ( CISCI_instance%diagonalCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)

    !! eigenvalues per SCI iteration
    allocate ( CISCI_instance%eigenValues ( 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS + CONTROL_instance%CI_SCI_REFINEMENT_STEPS ) )
    allocate ( CISCI_instance%minCoeff ( 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS + CONTROL_instance%CI_SCI_REFINEMENT_STEPS ) )
    do k = 1, 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS +  CONTROL_instance%CI_SCI_REFINEMENT_STEPS
      call Vector_constructor ( CISCI_instance%eigenValues(k), int(CONTROL_instance%CI_NUMBER_OF_STATES,8), 0.0_8) !! store the eigenvalues per macro iterations
    enddo

    !! target space size per iteration
    allocate ( CISCI_instance%targetSpaceSize_iter ( 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS + CONTROL_instance%CI_SCI_REFINEMENT_STEPS ) )

    !! initialize buffer arrays: diagonalTarger, omp_targetInterval, omp_target_iterator_m
    !!                           index_amplitudeCore, buffer_amplitudeCore, confAmplitudeCore_orb
    call CISCI_resetBuffer()

    !! initialize sorting subroutines
    call CISort_constructor()

    !! auxiliary orbial vector in canonical order (1,2,3,4...)
    allocate ( CISCI_instance%canonicalOrder ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( CISCI_instance%canonicalOrder(spi), CIcore_instance%numberOfActiveOrbitals%values(spi), 0 )
      do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
        CISCI_instance%canonicalOrder(spi)%values(pi) = pi
      enddo
    enddo

    if ( CONTROL_instance%CI_SELECTIVE_METHOD == "HBCI" ) then
      call CISCI_heatbathIntegralSorting( CISCI_instance%heatBathDoubleExcitations, CISCI_instance%heatBathDoubleExcitations_index, &
                                          CISCI_instance%heatBathDoubleExcitations_size )
    endif

  end subroutine CISCI_constructor

  !! Allocating arrays 
  subroutine CISCI_destructor( )
    implicit none

    deallocate ( CISCI_instance%combinedOrbitalsPositions )
    deallocate ( CISCI_instance%combinedOccupiedOrbitalsPositions )
    deallocate ( CISCI_instance%omp_targetInterval ) !fixed
    deallocate ( CISCI_instance%omp_target_iterator_m ) !variable

    !! arrays for storing coefficients
    call Vector_destructor(CISCI_instance%buffer_amplitudeCore)
    call Vector_destructor(CISCI_instance%coefficientCore)

    !! auxiliary array to trace the original elements in array sorting
    call Vector_destructorInteger8 ( CISCI_instance%index_amplitudeCore)

    !! arrays for storing CI configurations, species, orbitals, vector size
    deallocate ( CISCI_instance%confCore ) 
    deallocate ( CISCI_instance%confAmplitudeCore_orb ) 

    call Vector_destructor ( CISCI_instance%diagonalCore )

    !! eigenvalues per SCI iteration
    deallocate ( CISCI_instance%eigenValues )
    deallocate ( CISCI_instance%minCoeff )

    !! target space size per iteration
    deallocate ( CISCI_instance%targetSpaceSize_iter )
    deallocate ( CISCI_instance%canonicalOrder )

    !! CIorder list
    deallocate ( CISCI_instance%maxCIexcitations )
    deallocate ( CISCI_instance%CIorder_list )
    deallocate ( CISCI_instance%CIorder_count )
    deallocate ( CISCI_instance%CIorder_weight )

    call CISort_destructor()

  end subroutine CISCI_destructor
 
  !! main part
  subroutine CISCI_run( numberOfConfigurations, eigenVectors, initialEnergy, initialStep, unboundReference, computePT2 )
    use sort_
    implicit none
    integer(8), intent(inout) :: numberOfConfigurations
    type(matrix), intent(inout) :: eigenVectors
    real(8), intent(in) :: initialEnergy 
    logical, intent(in):: initialStep, unboundReference, computePT2
    real(8) :: currentEnergy 
    integer(8) :: a, aa, i, j, ii, jj
    integer(8) :: o1, o2
    integer(8) :: m, m1, m2
    integer :: k, finalk ! macro SCI iteration
    integer :: nproc, n
    real(8) :: timeA(20), timeB(20)
    real(8) :: timeAA, timeBB
    real(8) :: timeAS, timeBS
    real(8) :: minValue
    integer :: numberOfSpecies, spi
    logical :: use_guess
    type(matrix) :: hamiltonianMatrix

    numberOfSpecies = CIcore_instance%numberOfSpecies 
    nproc = CIcore_instance%nproc 

    k = 1
    CISCI_instance%eigenValues(k)%values(1) = initialEnergy
    currentEnergy = CISCI_instance%eigenValues(k)%values(1)

    !! initial step
    if ( initialStep ) then
      use_guess = .false. !! usually HF is bad guess
      write (6,*)    ""
      write (6,"(T2,A29 )")    "Starting SCI macro iterations "
      write (6,*)    ""
      call CISCI_initialConfigurations( CISCI_instance%coefficientCore, CISCI_instance%confCore )
    endif

    if ( .not. initialStep ) then
      use_guess = .true.
      write (6,*)    ""
      write (6,"(T2,A29 )")    "Re-starting SCI macro iterations "
      write (6,*)    ""
    endif

    do k = 2, 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS + CONTROL_instance%CI_SCI_REFINEMENT_STEPS
        
!$  timeA(k) = omp_get_wtime()

      CISCI_instance%targetSpaceSize_iter(k) = CISCI_instance%targetSpaceSize

      !! computing the diagonal in the core space, for fast computation of core amplitudes
      call CISCI_buildDiagonal ( CISCI_instance%diagonalCore, CISCI_instance%confCore, CISCI_instance%coreSpaceSize )

      !! search and select
      select case ( CONTROL_instance%CI_SELECTIVE_METHOD )

      case ("ASCI")
        !! calculating the amplitudes in core space. This is the pertubation guess of CI eigenvector
        if ( CIcore_instance%level == "CISD-" .or. unboundReference ) then
          call CISCI_core_amplitudes_cisd ( CISCI_instance%coefficientCore%values, CISCI_instance%confCore, CISCI_instance%coreSpaceSize, currentEnergy )
        endif
        if ( CIcore_instance%level == "FCI" .and. .not. unboundReference ) then
          call CISCI_core_amplitudes ( CISCI_instance%diagonalCore, CISCI_instance%coefficientCore%values, CISCI_instance%confCore, CISCI_instance%coreSpaceSize, currentenergy )
        endif
      case ("HBCI")
          call CISCI_heatBathGenerate (  CISCI_instance%diagonalCore, CISCI_instance%coefficientCore%values,&
          CISCI_instance%confCore, CISCI_instance%coreSpaceSize, currentenergy, VARIATIONAL )
      end select

      !! merge core space with the top amplitude to form a new target space, and save them in saved_conf with coefficients in eigenvectors
      call CISCI_mergeCoreAndTarget( CISCI_instance%confTarget_orb, eigenVectors )

      !! computing the diagonal in the target space, jadamilu requires the diagonal in advance
      call CISCI_buildDiagonal ( CISCI_instance%diagonalTarget, CISCI_instance%confTarget_orb, CISCI_instance%targetSpaceSize )

      !! eigenvalue guess
      CISCI_instance%eigenValues(k)%values(1) = currentEnergy

      !! diagonalize in target space
      select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

      case ("DSYEVR")

        if ( .not. use_guess ) eigenVectors%values = 0.0_8

        !!build full matrix and use lapack...
        call Matrix_constructor ( hamiltonianMatrix, int(CISCI_instance%targetSpaceSize,8), int(CISCI_instance%targetSpaceSize,8), 0.0_8 )
        call CISCI_buildHamiltonian ( hamiltonianMatrix )
        call Matrix_eigen_select ( hamiltonianMatrix, CISCI_instance%eigenValues(k), &
                 int(1), int(CONTROL_instance%CI_NUMBER_OF_STATES), &  
                 eigenVectors =eigenVectors, &
                 flags = int(SYMMETRIC,4))

      case ("JADAMILU")
        call CISCI_jadamiluInterface( int(CISCI_instance%targetSpaceSize,8), &
                 1_8, &
                 CISCI_instance%eigenValues(k), &
                 eigenVectors, timeAA, timeBB, use_guess )

      end select
  
      write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k-1 , " Energy: ", CISCI_instance%eigenValues(k)%values(1)

      !! if the correlation energy is positive, then don't use the guess
      if (  CISCI_instance%eigenValues(k)%values(1) - HartreeFock_instance%totalEnergy < 0 ) use_guess = .true.

      !! reset auxindex arrary, for later use in sorting target coeff. global absolute index
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo

      !! copy confTarget_orb to confTarget_occ ( same order)
      call CISCI_orb2occ()

      !! getting the core absolute largest coefficients
      call CISort_quicksort_vector(  eigenVectors%values(:,1), &
                                     CISCI_instance%index_amplitudeCore%values(1:CISCI_instance%targetSpaceSize), & 
                                     1_8,  int(CISCI_instance%targetSpaceSize,8)  )
    
      !! just some diagnostics, average of last 5 coefficients
      CISCI_instance%minCoeff(k) = sum(eigenVectors%values(CISCI_instance%targetSpaceSize-5:CISCI_instance%targetSpaceSize,1) ) / 5.0

      !! copy confTarget_occ to confTarget_orb (sorted by index_amplitude)
      call CISCI_occ2orb()

      !! copy confTarget_orb to confTarget_occ ( now sortered)
      call CISCI_orb2occ()

!$  timeB(k) = omp_get_wtime()

      finalk = k

      !! convergence criteria. Exit here avoiding matrices reset if: the energy converged or reach max iter, and if at least 3 iterations were achieved  
      if ( abs( CISCI_instance%eigenValues(k)%values(1) - currentEnergy ) < 1.0E-5 .and. k > 2 ) then
        write (6,"(T2,A30)") "Reached SCI Energy Convergence of 1E-5 "
        exit
      end if

      !! don't grow anymore
      if ( k == 1 + CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS + CONTROL_instance%CI_SCI_REFINEMENT_STEPS ) then
        write (6,"(T2,A30)") "Reached Max number of steps "
        exit
      end if

      !! preparation for next iter

      !! storing only the largest coefficients, and rearraing the next eigenvector guess
      do i = 1,  CISCI_instance%coreSpaceSize
        CISCI_instance%coefficientCore%values(i) = eigenVectors%values(i,1)
      enddo

      !! storing the top sorted target conf into the core conf space
      do i = 1,  CISCI_instance%coreSpaceSize
        do spi = 1, numberOfSpecies
          CISCI_instance%confCore(spi)%values(:,i) = CISCI_instance%confTarget_orb(spi)%values(:,i)
        enddo
      enddo

      !! set target space, either grow or refine
      if ( k <= CONTROL_instance%CI_SCI_TARGET_GROWTH_STEPS ) then
        CISCI_instance%targetSpaceSize = CISCI_instance%targetSpaceSize * int (CONTROL_instance%CI_SCI_TARGET_GROWTH_FACTOR )
      else
        CISCI_instance%targetSpaceSize = CISCI_instance%targetSpaceSize
      end if
      numberOfConfigurations = CISCI_instance%targetSpaceSize

      !! growing eigenvector size
      call Matrix_constructor (eigenVectors, &
           int(CISCI_instance%targetSpaceSize, 8), &
           int(CONTROL_instance%CI_NUMBER_OF_STATES,8), 0.0_8)

      !! updating new reference
      currentEnergy = CISCI_instance%eigenValues(k)%values(1) 

      !! reallocate due to increased target size
      do spi = 1, CIcore_instance%numberOfSpecies
        call Matrix_constructorInteger1 ( CISCI_instance%confTarget_orb(spi), int(CIcore_instance%numberOfActiveOrbitals%values(spi),8) , int(CISCI_instance%targetSpaceSize,8), -1_1)
        call Matrix_constructorInteger ( CISCI_instance%confTarget_occ(spi), int(CIcore_instance%numberOfOccupiedOrbitals%values(spi),8) , int(CISCI_instance%targetSpaceSize,8), -1_4)
      enddo

      !! reset iterators for next iter
      call CISCI_resetBuffer()

    enddo !k

    !! summary of the macro iteration 
    write (6,*)    ""
    write (6,"(T2,A107 )")    "                                    Selected CI (SCI)  summary                                            "
    write (6,"(T2,A6,A12,A25,A25,A20,A15,A12)")    "Iter", "TargetSize", "GS Energy", "Corr. Energy", "Delta E.", "Min coeff.","Time(s)"
    do k = 2, finalk
       write (6,"(T2,I6, I12, F25.12, F25.12, F20.12,  ES15.2, F12.2 )") k-1, CISCI_instance%targetSpaceSize_iter(k),  CISCI_instance%eigenValues(k)%values(1),  &
                                                          CISCI_instance%eigenValues(k)%values(1) - HartreeFock_instance%totalEnergy, &
                                                          CISCI_instance%eigenValues(k)%values(1) - CISCI_instance%eigenValues(k-1)%values(1), &
                                                          CISCI_instance%minCoeff(k), &
                                                          timeB(k) - timeA(k)
    enddo !k

    !! save final eigenValues to CIcore instance
    CIcore_instance%eigenValues%values(1) = CISCI_instance%eigenValues(finalk)%values(1)
    write (6,"(T2,A,ES12.4)") "Minimum coefficient in target space: ", CISCI_instance%minCoeff(finalk) 
    write (6,*)    ""

    !! recalculate n conf to prevent the code breaking in case the targetSpace is not fully filled
    numberOfConfigurations = 0
    do a = 1, CISCI_instance%targetSpaceSize
      !if (CISCI_instance%confAmplitudeCore_orb(1,a) == -1_1  ) exit
      if (CISCI_instance%confTarget_orb(1)%values(1,a) == -1_1 ) exit
      numberOfConfigurations = numberOfConfigurations + 1
    enddo

    call CISCI_countSpeciesPairs()

    !! calculating PT2 correction. A pertuberd estimation of configurations not include in the target space
    if ( computePT2 ) then

      !! reset iterators
      call CISCI_resetBuffer()

      m = 1
      !! add the final target configurations at the beginning of the buffer array. in such way, only the non-duplicated connected configurations will be added
      do n = 1, CIcore_instance%nproc
        m1 = CISCI_instance%omp_targetInterval(1, n ) !! position to add 
        m2 = m1 + CISCI_instance%targetSpaceSize / CIcore_instance%nproc - 1 !! number of conf added

        CISCI_instance%omp_target_iterator_m(n) = m2
        CISCI_instance%buffer_amplitudeCore%values(m1:m2) = 1.0E+6 !! big number to ensure this conf won't be discarded after sorting
        do spi = 1, numberOfSpecies

          o1 = CISCI_instance%combinedOrbitalsPositions(1,spi) 
          o2 = CISCI_instance%combinedOrbitalsPositions(2,spi)

          CISCI_instance%confAmplitudeCore_orb( o1:o2, m1:m2 ) = CISCI_instance%confTarget_orb(spi)%values(:, m: m + m2 - m1 )
        enddo
        !! m is the position in confTarget (contiguous), m1 and m2 position in buffer array
        m = m + m2 - m1 + 1
      enddo

      !! computing the diagonal in the target space, for fast computation of core amplitudes
      call CISCI_buildDiagonal ( CISCI_instance%diagonalTarget, CISCI_instance%confTarget_orb, CISCI_instance%targetSpaceSize )

      select case ( CONTROL_instance%CI_SELECTIVE_METHOD )
        case ("ASCI")
        !! recompute amplitudes, but now from target space not core, in this way all connected conf are saved in buffer
        if ( CIcore_instance%level == "CISD-" ) then
          call CISCI_core_amplitudes_cisd (  eigenVectors%values(:,1), CISCI_instance%confTarget_orb, CISCI_instance%targetSpaceSize, currentEnergy )
        endif
        if ( CIcore_instance%level == "FCI" ) then
          call CISCI_core_amplitudes ( CISCI_instance%diagonalTarget, eigenVectors%values(:,1), CISCI_instance%confTarget_orb, CISCI_instance%targetSpaceSize, currentEnergy )
        endif
      case ("HBCI")
          call CISCI_heatBathGenerate ( CISCI_instance%diagonalTarget, eigenVectors%values(:,1), CISCI_instance%confTarget_orb, &
          CISCI_instance%targetSpaceSize, currentEnergy, PERTURBATIVE )
      end select

      !! the real PT2 calculation
      call CISCI_PT2 ( CISCI_instance%targetSpaceSize, CIcore_instance%eigenValues%values(1), CISCI_instance%PT2energy, eigenVectors )
    endif

  end subroutine CISCI_run

  !! reset buffer and other temporary arrays for a next SCI iteration
  subroutine CISCI_resetBuffer()
    implicit none
    integer :: n, m, i
    integer :: spi

    !! storing the CI diagonal matrix elements for Jadamilu preconditioner
    call Vector_constructor ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8)

    !! auxiliary arrays to store the position of the target space for each omp thread within the big arrays 
    !CISCI_instance%omp_targetInterval ! reminder: this one is fixed
    !CISCI_instance%omp_target_iterator_m ! reminder: this one is variable
    m = 0_8
    do n = 1, CIcore_instance%nproc 
      CISCI_instance%omp_targetInterval(1, n ) = m + 1_8
      CISCI_instance%omp_targetInterval(2, n ) = m + CISCI_instance%buffer_amplitudeCoreSize / CIcore_instance%nproc
      m = m + CISCI_instance%buffer_amplitudeCoreSize / CIcore_instance%nproc
      CISCI_instance%omp_target_iterator_m(n) = CISCI_instance%omp_targetInterval(1, n ) - 1
    enddo

    CISCI_instance%omp_targetInterval(1, CIcore_instance%nproc + 1 ) = 1_8
    CISCI_instance%omp_targetInterval(2, CIcore_instance%nproc + 1 ) = CISCI_instance%buffer_amplitudeCoreSize 
    CISCI_instance%omp_target_iterator_m( CIcore_instance%nproc + 1 ) = 0_8

    !! reset auxindex array. relative indexes, this index is relative for each omp thread to simplify internal usage during sorting
    i = 1
    do n = 1, CIcore_instance%nproc 
      do m = 1, CISCI_instance%omp_targetInterval(2, n ) - CISCI_instance%omp_targetInterval(1, n ) + 1  
         CISCI_instance%index_amplitudeCore%values(i) = m
         i = i + 1
      enddo
    enddo

    !! restart amplitudes for next run
    CISCI_instance%buffer_amplitudeCore%values = 0.0_8

    !! restart configurations for next run
    do spi = 1, CIcore_instance%numberOfSpecies
      CISCI_instance%confAmplitudeCore_orb = -1_1
      !CISCI_instance%confTarget_orb(spi)%values = -1_1
    enddo 


  end subroutine CISCI_resetBuffer
  

  !! compute the reference configuration, the HF 
  subroutine CISCI_initialConfigurations ( coefficientCore, confCore )

    implicit none
    type(Vector) :: coefficientCore
    type(IMatrix1) :: confCore(:)
    type(IVector), allocatable :: orbA(:), occA(:), virA(:)
    integer :: spi, spj, numberOfSpecies
    integer(8) :: m 
    real(8) :: indexConf
    integer :: pi, qi
    integer :: oia, via
    integer :: oi1, vi1

    !! Hartree-Fock reference coeff
    m = 1
    coefficientCore%values(m) = 0.50_8

    !! build orbitals references
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( occA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies 
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0) 
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
    
      do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
        orbA(spi)%values(pi) = 1.0
      enddo

      confCore(spi)%values(:,m) = orbA(spi)%values(:)
      !!call CISCI_binaryToDecimal ( orbA(spi)%values, indexConf )
      !!confCore(spi)%values(spi,m) = indexConf
      
    enddo

    !! add all single excited positronic states, useful for unbound HF references 
    if ( CONTROL_instance%CI_UNBOUND_REFERENCE ) then
      coefficientCore%values(m) = 0.10_8 
      singles: do spi = 1, numberOfSpecies 
        if ( trim(  MolecularSystem_getNameOfSpecies( spi ) ) == "E+" ) then

          oia = 0 
          via = 0

          !! build the orbital from the index using the bit mapping
          !!call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )

          !! build auxiliary vectors of occupied and virtuals orbitals
          do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
            if ( orbA(spi)%values(pi) == 1 ) then
              oia = oia + 1
              occA(spi)%values(oia) = pi
            else if ( orbA(spi)%values(pi) == 0 ) then
              via = via + 1
              virA(spi)%values(via) = pi
            end if
          enddo !pi

          !! single excitations
          do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            oi1 = occA(spi)%values(pi)  
            orbA(spi)%values(oi1) = orbA(spi)%values(oi1) - 1 
            do qi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
              vi1 = virA(spi)%values(qi)
              orbA(spi)%values(vi1) = orbA(spi)%values(vi1) + 1
              m = m + 1
              if ( m > CISCI_instance%coreSpaceSize ) exit singles

              !! add the configuration
              coefficientCore%values(m) = 0.10_8

              !! save all species
              do spj = 1, numberOfSpecies 
                confCore(spj)%values(:,m) = orbA(spj)%values(:)
              enddo

              orbA(spi)%values(vi1) = orbA(spi)%values(vi1) - 1
            enddo !qi
            orbA(spi)%values(oi1) = orbA(spi)%values(oi1) - 1
          enddo !pi

        endif ! E+
      enddo singles
    endif !endif unbound

    do spi = 1, numberOfSpecies 
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( orbA(spi)) 
      call Vector_destructorInteger ( virA(spi) )  
    enddo

    deallocate ( occA  )
    deallocate ( virA  )
    deallocate ( orbA  )

  end subroutine CISCI_initialConfigurations

  !! compute the estimated amplitude to form the target space from the core space
  subroutine CISCI_core_amplitudes ( diagonal, coefficientCore, confCore, SCICoreSpaceSize, oldEnergy )

    implicit none
    type(Vector), intent(in) :: diagonal
    real(8), intent(in) :: coefficientCore ( SCICoreSpaceSize )
    type(IMatrix1), intent(in) :: confCore(:)
    integer(8), intent(in) :: SCICoreSpaceSize
    real(8), intent(in) :: oldEnergy
    real(8) :: CIEnergy
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b,c, aa
    integer :: spi, spj, numberOfSpecies
    integer(8), allocatable :: indexConfA(:) !! ncore, species
    real(8) :: diagEnergy, diagEnergy_a
    real(8) :: diagEnergy_ao1, diagEnergy_ao1o2
    real(8) :: shift
    type (ivector), allocatable :: occA(:), occB(:), virA(:), virB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer, allocatable :: CIlevel(:)
    real(8) :: tmpconfCoreConfB
    integer(8) :: pi, qi, ri, si, pj, qj, rj, sj
    integer(8) :: oia, oja, via, vja, aaa
    integer(8) :: oi1, vi1, oi2, vi2, oj2, vj2
    integer :: factor1, factor2, factor2j
    integer :: nonzero

!$  timeA = omp_get_wtime()
    shift = 1E-8 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !! work only with non-zero conf
    nonzero = 0
    do a = 1, SCICoreSpaceSize  
      if ( confCore(1)%values(1,a) == -1_1 .or. abs(coefficientCore(a)) <= 1E-10 ) exit
      nonzero = nonzero + 1
    enddo

    !$omp parallel &
    !$omp& private ( occA, occB, virA, virB, orbA, orbB, CIlevel) &
    !$omp& private ( n, a, oia, via, pi, qi, ri, si, oi1, vi1, oi2, vi2, spi, spj, oj2, vj2, factor1, factor2, factor2j, &
    !$omp&           CIenergy, diagEnergy, diagEnergy_a, diagEnergy_ao1, diagEnergy_ao1o2 ) 

    !! allocating auxiliary arrays (omp) for working with conf and orbitals
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( CIlevel ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    n = omp_get_thread_num() + 1

    !! loop to find all CI configurtions coupled to core space
    !!$omp do schedule (runtime) !with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    do a = 1, nonzero  

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0_8
        via = 0_8

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = confCore(spi)%values(:,a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1_8 ) then
            oia = oia + 1_8
            occA(spi)%values(oia) = pi
          else if ( orbA(spi)%values(pi) == 0_8 ) then
            via = via + 1_8
            virA(spi)%values(via) = pi
          end if
        enddo

        !! copy to conf B, these are variable
        orbB(spi)%values = orbA(spi)%values 
        occB(spi)%values = occA(spi)%values 
        virB(spi)%values = virA(spi)%values 

      enddo
      CIlevel = 0

      !! Use diagonal(a) as the reference energy for the diagonal elements in b
      diagEnergy_a = diagonal%values(a)

      !! building all single sustitutions from configuration A. 
      !! here all configurations pairs are generated in maximum coincidence 
      do spi = 1, numberOfSpecies 

        !! calculate the sign factor for canonical order of the configuration
        factor1 = CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )

        do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1_8

          !! remove energy from the excited orbital
          diagEnergy_ao1 = diagEnergy_a - CISCI_calculateEnergyOne( spi, occA, occA, oi1, oi1 )

          do qi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi) !! occ or core???
            vi1 = virA(spi)%values(qi)
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) + 1_8
            occB(spi)%values(pi) = vi1

            CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

            !! bit mapping from orbital to decimal num
            !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )
            !confCoreConfB(spi)%values = orbB(spi)%values ! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! get spingle sustitutions energy, H_ab C_a
            CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, oi1, vi1  )
            CIenergy = CIenergy * coefficientCore(a) 

            !! calculate the sign factor for canonical order of the configuration
            factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
            CIenergy = CIenergy * factor1 * factor2

            !! add energy of the excited orbital to diagonal, H_bb
            diagEnergy = diagEnergy_ao1 + CISCI_calculateEnergyOne( spi, occB, occB, vi1, vi1 )
            !! alternative way (slower)
            !!diagEnergy = CISCI_calculateEnergyZero( occB )

            !! Amplitude estimation A_b = H_ab C_a / ( H_bb - E_ref )
            CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

            if ( CIenergy /= 0.0_8 ) then
              !! append the amplitude 
              call CISCI_appendAmplitude ( n, CIenergy, orbB )
            endif

            !tmpconfCoreConfB = confCoreConfB(spi) !! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! building all double intraspecies sustitutions from configuration A
            do ri = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
              oi2 = occA(spi)%values(ri)  
              if ( oi1 <= oi2 ) cycle 
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) - 1_8

              !! remove energy from the 2nd excited orbital
              diagEnergy_ao1o2 = diagEnergy_ao1 - CISCI_calculateEnergyOne( spi, occA, occA, oi2, oi2 )
              diagEnergy_ao1o2 = diagEnergy_ao1o2 + CISCI_calculateEnergyTwoSame( spi, occA, occA, oi1, oi2, oi1, oi2 )

              do si = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
                vi2 = virA(spi)%values(si)
                if ( vi1 <= vi2 ) cycle 
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) + 1_8
                occB(spi)%values(ri) = vi2

                !! bit mapping from orbital to decimal num
                !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )

                CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

                !! get double intraspecies sustitutions energy
                CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, oi1, oi2, vi1,  vi2 )
                CIenergy = CIenergy * coefficientCore(a)

                !! calculate the sign factor for canonical order of the configuration
                factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
                CIenergy = CIenergy * factor1 * factor2

                !! add energy of the 1st and 2nd excited orbital to the diagonal H_bb
                diagEnergy = diagEnergy_ao1o2 + CISCI_calculateEnergyOne( spi, occB, occB, vi1, vi1 )
                diagEnergy = diagEnergy + CISCI_calculateEnergyOne( spi, occB, occB, vi2, vi2 )
                diagEnergy = diagEnergy - CISCI_calculateEnergyTwoSame( spi, occB, occB, vi1, vi2, vi1, vi2 )
                !!alternative way (slower)
                !!diagEnergy = CISCI_calculateEnergyZero( occB )

                !! Amplitude estimation A_b = H_ab C_a / ( H_bb - E_ref )
                CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                if ( CIenergy /= 0.0_8 ) then
                  !! append the amplitude 
                  call CISCI_appendAmplitude ( n, CIenergy, orbB )
                endif

                occB(spi)%values(ri) = occA(spi)%values(ri)  
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) - 1_8 ! reset orbital 
              enddo
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) + 1_8 ! reset orbital
            enddo

            CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

            !! building all double interspecies sustitutions from configuration A. maybe build this as a superloop?
            !! get double interspecies sustitutions energy
            do spj = spi + 1, numberOfSpecies 
              do rj = CIcore_instance%numberOfCoreOrbitals%values(spj) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                oj2 = occA(spj)%values(rj)  
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) - 1_8

                !! remove energy from the 2nd excited orbital
                diagEnergy_ao1o2 = diagEnergy_ao1 - CISCI_calculateEnergyOne( spj, occA, occA, oj2, oj2 )
                diagEnergy_ao1o2 = diagEnergy_ao1o2 + CISCI_calculateEnergyTwoDiff( spi, spj, oi1, oj2, oi1, oj2 )

                do sj = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj) - CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                  vj2 = virA(spj)%values(sj)
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) + 1_8
                  occB(spj)%values(rj) = vj2

                  !! bit mapping from orbital to decimal num
                  !call CISCI_binaryToDecimal ( orbB(spj)%values, confCoreConfB(spj) )

                  CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj)+1:) )

                  !! get double interspecies sustitutions energy
                  CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, oi1, oj2, vi1, vj2 )
                  CIenergy = CIenergy * coefficientCore(a) 

                  !! calculate the sign factor for canonical order of the configuration
                  factor2j = CISCI_canonicalOrderFactor( spj, orbB(spj), occB(spj) )
                  CIenergy = CIenergy * factor2 * factor2j

                  !! add energy of the 1st and 2nd excited orbital to the diagonal H_bb
                  diagEnergy = diagEnergy_ao1o2 + CISCI_calculateEnergyOne( spi, occB, occB, vi1, vi1 )
                  diagEnergy = diagEnergy + CISCI_calculateEnergyOne( spj, occB, occB, vj2, vj2 )
                  diagEnergy = diagEnergy - CISCI_calculateEnergyTwoDiff( spi, spj, vi1, vj2, vi1, vj2 )
                  !! alternative way (slower)
                  !!diagEnergy = CISCI_calculateEnergyZero( occB )

                  !! Amplitude estimation A_b = H_ab C_a / ( H_bb - E_ref )
                  CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                  if ( CIenergy /= 0.0_8 ) then
                    !! append the amplitude 
                    call CISCI_appendAmplitude ( n, CIenergy, orbB )
                  endif

                  !! reset the confB
                  occB(spj)%values(rj) = occA(spj)%values(rj)  
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) - 1_8
                enddo ! sj
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) + 1_8
                CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj) + 1_8:) )
              enddo ! rj

              CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj) + 1_8:) )
            enddo !spj

            !! reset the confB
            occB(spi)%values(pi) = occA(spi)%values(pi)  
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) - 1_8
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1_8
       enddo !pi
       CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi) + 1_8:) )

      enddo !spi
      
    enddo !enddo a
    !$omp enddo 

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( virA(spi) )  
      call Vector_destructorInteger ( virB(spi) )  
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( CIlevel )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    !$omp end parallel

    !! ------------------------------
    !! the above code applies the denominator from eq 4 10.1063/1.4955109 on the fly for each configuration. Alternatively, it can be sort by heat-bath ranking first, then final sort by asci ranking

    !! sort and reduce the target arrays among all OMP threads, final run
    call CISCI_sortAmplitude( CIcore_instance%nproc + 1 ) 

!$  timeB = omp_get_wtime()
!$  write(*,"(A,ES10.2,A4)") "** TOTAL Elapsed Time for calculating SCI amplitudes : ", timeB - timeA ," (s)"

  end subroutine CISCI_core_amplitudes

  !! compute the estimated amplitude to form the target space from the core space
  subroutine CISCI_core_amplitudes_cisd ( coefficientCore, confCore, SCICoreSpaceSize, oldEnergy )

    implicit none
    integer(8) SCICoreSpaceSize
    real(8) coefficientCore ( SCICoreSpaceSize )
    type(IMatrix1) :: confCore(:)
    real(8) :: CIEnergy
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b,c, aa
    integer :: spi, spj, numberOfSpecies
    integer(8), allocatable :: indexConfA(:) !! ncore, species
    real(8) :: diagEnergy
    real(8) :: oldEnergy
    real(8) :: shift
    type (ivector), allocatable :: occA(:), occB(:), virA(:), virB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer, allocatable :: CIlevel(:)
    real(8) :: tmpconfCoreConfB
    integer(8) :: pi, qi, ri, si, pj, qj, rj, sj
    integer(8) :: oia, oja, via, vja, aaa
    integer(8) :: oi1, vi1, oi2, vi2, oj2, vj2
    integer :: factor1, factor2, factor2j
    integer :: nonzero

!$  timeA = omp_get_wtime()
    shift = 1E-4 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !! work only with non-zero conf
    nonzero = 0
    do a = 1, SCICoreSpaceSize  
      if ( confCore(1)%values(1,a) == -1_1 ) exit
      nonzero = nonzero + 1
    enddo

    !$omp parallel &
    !$omp& private ( occA, occB, virA, virB, orbA, orbB, CIlevel) &
    !$omp& private ( n, a, oia, via, pi, qi, ri, si, oi1, vi1, oi2, vi2, spi, spj, oj2, vj2, factor1, factor2, factor2j, CIenergy, diagEnergy ) 

    !! allocating auxiliary arrays (omp) for working with conf and orbitals
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( CIlevel ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    n = omp_get_thread_num() + 1

    !! loop to find all CI configurtions coupled to core space
    !!$omp do schedule (runtime) !with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    do a = 1, nonzero  
      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 
        via = 0

        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = confCore(spi)%values(:,a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          else if ( orbA(spi)%values(pi) == 0 ) then
            via = via + 1
            virA(spi)%values(via) = pi
          end if
        enddo

        !! copy to conf B, these are variable
        orbB(spi)%values = orbA(spi)%values 
        occB(spi)%values = occA(spi)%values 
        virB(spi)%values = virA(spi)%values 

      enddo
      CIlevel = 0

      !! building all single sustitutions from configuration A. 
      !! here all configurations pairs are generated in maximum coincidence 
      do spi = 1, numberOfSpecies 

        !! calculate the sign factor for canonical order of the configuration
        factor1 = CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )

        do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1 

          do qi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            vi1 = virA(spi)%values(qi)
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) + 1
            occB(spi)%values(pi) = vi1

            !! get spingle sustitutions energy
            CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, oi1, vi1  )

            !! calculate the sign factor for canonical order of the configuration
            factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )

            CIenergy = CIenergy * factor1 * factor2

            !! bit mapping from orbital to decimal num
            !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )
            !confCoreConfB(spi)%values = orbB(spi)%values ! save the indexconfB to use later in double inter, because double intra will overwritten it 
            CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

            if ( CIenergy /= 0.0_8 .and. sum(CIlevel) <= 2 .and. sum(CIlevel(1:2)) <= 1 ) then
            !if ( CIenergy /= 0.0_8 ) then
              CIenergy = CIenergy * coefficientCore(a) 

              ! alternative approximation
              diagEnergy = CISCI_calculateEnergyZero( occB )
              CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

              !! append the amplitude 
              call CISCI_appendAmplitude ( n, CIenergy, orbB )
            endif

            !tmpconfCoreConfB = confCoreConfB(spi) !! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! building all double intraspecies sustitutions from configuration A
            do ri = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
              oi2 = occA(spi)%values(ri)  
              if ( oi1 <= oi2 ) cycle 
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) - 1 
              do si = 1, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
                vi2 = virA(spi)%values(si)
                if ( vi1 <= vi2 ) cycle 
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) + 1
                occB(spi)%values(ri) = vi2 

                !! get double intraspecies sustitutions energy
                CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, oi1, oi2, vi1,  vi2 )
                CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

                !if ( CIenergy /= 0.0_8 ) then
                if ( CIenergy /= 0.0_8 .and. sum(CIlevel) <= 2 .and. sum(CIlevel(1:2)) <= 1 ) then
                  !! calculate the sign factor for canonical order of the configuration
                  factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
 
                  CIenergy = CIenergy * factor1 * factor2

                  !! bit mapping from orbital to decimal num
                  !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )

                  CIenergy = CIenergy * coefficientCore(a)
                  ! alternative approximation
                  diagEnergy = CISCI_calculateEnergyZero( occB )
                  CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                  !! append the amplitude 
                  call CISCI_appendAmplitude ( n, CIenergy, orbB )

                endif

                occB(spi)%values(ri) = occA(spi)%values(ri)  
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) -1  ! reset orbital 
              enddo
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) + 1 ! reset orbital
            enddo

            CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

            !! building all double interspecies sustitutions from configuration A. maybe build this as a superloop?
            !! get double interspecies sustitutions energy
            do spj = spi + 1, numberOfSpecies 
              !if ( spj == spi ) cycle 
              do rj = CIcore_instance%numberOfCoreOrbitals%values(spj) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                oj2 = occA(spj)%values(rj)  
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) - 1 
                do sj = 1, CIcore_instance%numberOfActiveOrbitals%values(spj) - CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                  vj2 = virA(spj)%values(sj)
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) + 1
                  occB(spj)%values(rj) = vj2 

                  !! get double interspecies sustitutions energy
                  CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, oi1, oj2, vi1, vj2 )
                  CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj)+1:) )

                  !if ( CIenergy /= 0.0_8 ) then
                  if ( CIenergy /= 0.0_8 .and. sum(CIlevel) <= 2 .and. sum(CIlevel(1:2)) <= 1 ) then

                    !! calculate the sign factor for canonical order of the configuration
                    factor2j = CISCI_canonicalOrderFactor( spj, orbB(spj), occB(spj) )
   
                    CIenergy = CIenergy * factor2 * factor2j

                    !! bit mapping from orbital to decimal num
                    !call CISCI_binaryToDecimal ( orbB(spj)%values, confCoreConfB(spj) )

                    CIenergy = CIenergy * coefficientCore(a) 
                    ! alternative approximation
                    diagEnergy = CISCI_calculateEnergyZero( occB )
                    CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                    !! append the amplitude 
                    call CISCI_appendAmplitude ( n, CIenergy, orbB )
                  endif

                  !! reset the confB
                  occB(spj)%values(rj) = occA(spj)%values(rj)  
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) -1 
                enddo ! sj
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) + 1
                CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj)+1:) )
              enddo ! rj

              CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj)+1:) )
            enddo !spj

            !! reset the confB
            occB(spi)%values(pi) = occA(spi)%values(pi)  
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) -1 
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1
        enddo !pi
       CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

      enddo !spi
      
    enddo !enddo a
    !$omp enddo 

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( virA(spi) )  
      call Vector_destructorInteger ( virB(spi) )  
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( CIlevel )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    !$omp end parallel

    !! ------------------------------
    !! the above code applies the denominator from eq 4 10.1063/1.4955109 on the fly for each configuration. Alternatively, it can be sort by heat-bath ranking first, then final sort by asci ranking

    !! sort and reduce the target arrays among all OMP threads, final run
    call CISCI_sortAmplitude( CIcore_instance%nproc + 1 ) 

!$  timeB = omp_get_wtime()
!$  write(*,"(A,ES10.2,A4)") "** TOTAL Elapsed Time for calculating SCI amplitudes : ", timeB - timeA ," (s)"

  end subroutine CISCI_core_amplitudes_cisd


  subroutine CISCI_jadamiluInterface( N, numberOfeigenValues, eigenValues, eigenVectors, timeA, timeB, use_guess )
    implicit none
    external DPJDREVCOM
    integer(8), intent(in) :: numberOfeigenValues
    type(Vector), intent(inout) :: eigenValues
    type(Matrix), intent(inout) :: eigenVectors
    logical, intent(in) :: use_guess
    real(8), intent(out) :: timeA, timeB
    integer(4) :: iiter
    integer(8) :: i, j, k
    !Jadamilu variables
    !N: size of the problem
    !MAXSP: max. value of MADSPACE
    integer(8), intent(in) :: N
    integer(8) :: MAXSP
    integer(8) :: LX
    real(8), allocatable :: EIGS(:), RES(:), X(:)
    !arguments to pass to the routines
    integer(8) :: NEIG, MADSPACE, ISEARCH, NINIT
    integer(8) :: JA(1), IA(1)
    integer(8) :: ICNTL(5)
    integer(8) :: ITER, IPRINT, INFO
    real(8) :: SIGMA, TOL, GAP, MEM, DROPTOL, SHIFT
    integer(8) :: NDX1, NDX2
    integer(8) :: IJOB
    
!$  timeA = omp_get_wtime()
    MAXSP = CONTROL_instance%CI_MADSPACE

    LX = N*(3*MAXSP + numberOfeigenValues + 1) + 4*MAXSP*MAXSP

    if ( allocated ( EIGS ) ) deallocate ( EIGS )
    allocate ( EIGS ( numberOfeigenValues ) )
    EIGS = 0.0_8
    if ( allocated ( RES ) ) deallocate ( RES )
    allocate ( RES ( numberOfeigenValues ) )
    RES = 0.0_8
    if ( allocated ( X ) ) deallocate ( X )
    allocate ( X ( LX ) )
    X = 0.0_8

    !set input variables
    IPRINT = 0 !     standard report on standard output
    ISEARCH = 1 !    we want the smallest eigenvalues
    NEIG = numberOfeigenValues !    number of wanted eigenvalues
    MADSPACE = MAXSP !    desired size of the search space
    ITER = 30*NEIG !    maximum number of iteration steps
    TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual
    TOL = 1e-3 !1.0d-4 !    tolerance for the eigenvector residual, for ASCI this can be higher
    DROPTOL = 1E-3

    NDX1 = 0
    NDX2 = 0
    MEM = 20

    ! additional parameters set to default
    ICNTL(1)=0
    ICNTL(2)=0
    ICNTL(3)=0
    ICNTL(4)=0
    ICNTL(5)=0

    IJOB=0

    JA(1) = -1 
    IA(1) = -1 

    if ( use_guess ) then
      NINIT = NEIG !    initial approximate eigenvectors
      ! set initial eigenpairs
      do j = 1, N
        X(j) = eigenVectors%values(j,1)
      end do
    else
      NINIT = 0 !    no initial approximate eigenvectors
    endif

    do i = 1, CONTROL_instance%CI_NUMBER_OF_STATES
      EIGS(i) = eigenValues%values(i)
    end do


    SIGMA = EIGS(1)
    GAP = 0
    SHIFT = 0

    do i = 1, CONTROL_instance%CI_NUMBER_OF_STATES
      write(6,"(T2,A5,I4,2X,A10,F20.10,2X,A17,F10.6,A5,F10.6)") "State", i, "Eigenvalue", EIGS( i ), "Eigenvector. Max:", X((i-1)*N + i), "Min:", X(i*N )
    end do

    iiter = 0
  
10  CALL DPJDREVCOM( N, CISCI_instance%diagonalTarget%values , JA, IA, EIGS, RES, X, LX, NEIG, &
                       SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL, &
                       SHIFT, DROPTOL, MEM, ICNTL, &
                       IJOB, NDX1, NDX2, IPRINT, INFO, GAP)

    !! the private matrix-vector multiplication
    iiter = iiter +1
    IF (IJOB.EQ.1) THEN
      call CISCI_matvec ( N, X(NDX1), X(NDX2), iiter)
      GOTO 10
    END IF
  
    !! saving the eigenvalues
    eigenValues%values = EIGS

    !! saving the eigenvectors
    k = 0
    do j = 1, numberOfeigenValues
       do i = 1, N
        k = k + 1
        eigenVectors%values(i,j) = X(k)
      end do
    end do

    !! release internal memory and discard preconditioner
    CALL PJDCLEANUP
    if ( allocated ( X ) ) deallocate ( X )
    if ( allocated ( EIGS ) ) deallocate ( EIGS )
    if ( allocated ( RES ) ) deallocate ( RES )

!$  timeB = omp_get_wtime()

  end subroutine CISCI_jadamiluInterface

  subroutine CISCI_matvec ( NX, V, W, iter)
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    integer(8), intent(in) :: NX
    real(8), intent(in) :: V(NX)
    real(8), intent(out) :: W(NX)
    integer(4), intent(in) :: iter
    integer(8) :: a,b,aa,bb
    integer(8) :: nonzero, nonzerow
    real(8) :: tol
    integer :: i, spi, spj
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer(1) :: coupling
    integer(1), allocatable :: couplingS(:)
    integer(8) :: diffOrbi(4)
    integer(8) :: diffOrbj(4)
    integer :: pi
    integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA, factorB
    real(8), allocatable :: W_nproc(:,:)
    integer :: thread_id

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies

    nonzero = 0
    nonzerow = 0
    W = 0.0_8
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 
    do a = 1 , NX
       if ( abs(V(a) ) >= tol) nonzero = nonzero + 1
    end do

!$    timeA= omp_get_wtime()

    !! store matrix-vector product per omp thread
    allocate ( W_nproc ( NX, CIcore_instance%nproc ) )
    W_nproc = 0.0_8

    !! generate occupied orbital representation
    call CISCI_orb2occ()

    !$omp parallel &
    !$omp& private(thread_id, aa, a, spi, oia, orbA, pi, occA, CIenergy, bb, b, oib, orbB, occB, couplingS, coupling, i, diffOrbi, diffOrbj, spj, factorA, factorB )
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    thread_id = omp_get_thread_num() + 1

    !!$omp do schedule (runtime) ! with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    aloop: do aa = 1, NX

      !if ( abs(V(aa) ) <= tol) cycle ! this never happens

      !a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      a = aa ! if index_amplitude is sorted

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,a)
        occA(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:, a) 

      enddo

      bloop: do bb = aa, NX

        !b = CISCI_instance%index_amplitudeCore%values(bb)
        b = bb

        ! getting configuration B
        do spi = 1, numberOfSpecies 

          orbB(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,b)
          occB(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:,b) 

        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - dot_product ( orbA(spi)%values(:), orbB(spi)%values(:) )
        end do
      
        select case ( sum(couplingS) )
  
        !! same conf 
        case (0)
          CIenergy = CISCI_instance%diagonalTarget%values(aa) 
          W_nproc(aa,thread_id) = W_nproc(aa,thread_id) + CIenergy * V(aa)
        !! one orbital different
        case (1)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
          CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, diffOrbi(1), diffOrbi(3)  )

          W_nproc(bb,thread_id) = W_nproc(bb,thread_id) + CIenergy * V(aa) * factorA
          W_nproc(aa,thread_id) = W_nproc(aa,thread_id) + CIenergy * V(bb) * factorA

        case (2)
          select case (maxval(couplingS) )
          !! two orbital different, same species
          case (2)
            do i = 1, numberOfSpecies
              if ( couplingS(i) == 2 ) spi = i
            end do
  
            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
            CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  )
  
            W_nproc(bb,thread_id) = W_nproc(bb,thread_id) + CIenergy * V(aa) * factorA
            W_nproc(aa,thread_id) = W_nproc(aa,thread_id) + CIenergy * V(bb) * factorA

          !! two orbital different, different species
          case (1)
            do i = 1, numberOfSpecies
              if ( couplingS(i) == 1 ) then 
                spi = i
                exit
              end if
            end do
  
            do i = spi+1, numberOfSpecies
              if ( couplingS(i) == 1 ) spj = i
            end do
            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
            diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj)%values, orbB(spj)%values, occA(spj)%values, occB(spj)%values, factorB )
            CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  )
  
            W_nproc(bb,thread_id) = W_nproc(bb,thread_id) + CIenergy * V(aa) * factorA * factorB
            W_nproc(aa,thread_id) = W_nproc(aa,thread_id) + CIenergy * V(bb) * factorA * factorB

          end select ! maxval(couplingS)
        end select ! sum(couplingS)
      end do bloop !b
    end do aloop !a 
    !$omp end do nowait

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) )
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( orbA(spi) )
      call Vector_destructorInteger ( orbB(spi) )
    end do

    deallocate ( couplingS )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    !$omp end parallel

    do thread_id = 1, CIcore_instance%nproc
      W(:) = W(:) + W_nproc(:,thread_id) 
    enddo

    deallocate ( W_nproc )

!$  timeB = omp_get_wtime()

    !! to check how dense is the W vector
    do a = 1 , NX
       if ( abs(W(a) ) >= tol) nonzerow = nonzerow + 1
    end do

!$    write(*,"(A,I2,A,ES10.2,A2,I12,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero, nonzerow
    return

  end subroutine CISCI_matvec

  subroutine CISCI_buildHamiltonian ( hamiltonianMatrix )
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
  
    type(matrix) :: hamiltonianMatrix
    integer(8) :: a,b,aa,bb
    integer(8) :: nonzero, nonzerow
    real(8) :: tol
    integer :: uu,vv
    integer :: i, ii, jj, n, spi, spj
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer(1) :: coupling
    integer(1), allocatable :: couplingS(:)
    integer :: nproc
    integer(8) :: diffOrbi(4)
    integer(8) :: diffOrbj(4)
    integer :: pi
    integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA, factorB

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies

!$    timeA= omp_get_wtime()

    !! generate occupied orbital representation
    call CISCI_orb2occ()

    !$omp parallel &
    !$omp& private(aa, a, spi, oia, orbA, pi, occA, CIenergy, bb, b, oib, orbB, occB, couplingS, coupling, i, ii, diffOrbi, diffOrbj, spj, factorA, factorB ) 
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do
    !!$omp do schedule (runtime) !with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    aloop: do aa = 1, CISCI_instance%targetSpaceSize 

      !a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      a = aa ! if index_amplitude is sorted

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        !! build the orbital from the index using the bit-masking
        !!call CISCI_decimalToBinary ( CISCI_instance%confTarget_orb%values(spi,a), orbA(spi)%values )
        !!orbA(spi)%values = CISCI_instance%targetOrb(spi,a)%values
        orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,a)

        !! build auxiliary vectors of occupied and virtuals orbitals
        !do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
        !  if ( orbA(spi)%values(pi) == 1 ) then
        !    oia = oia + 1
        !    occA(spi)%values(oia) = pi
        !  end if
        !enddo
        occA(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:, a) 

      enddo

      bloop: do bb = aa, CISCI_instance%targetSpaceSize 


        !b = CISCI_instance%index_amplitudeCore%values(bb)
        b = bb
        ! getting configuration B
        do spi = 1, numberOfSpecies 

          !! build the orbital from the index using the bit mapping
          !! call CISCI_decimalToBinary ( CISCI_instance%confTarget_orb%values(spi,b), orbB(spi)%values )
          !! orbB(spi)%values = CISCI_instance%targetOrb(spi,b)%values
          orbB(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,b)
          occB(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:,b) 

        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) 
        end do
      
        !! if any species differs in more than 2 orb, skip
        if ( sum(couplingS) <= 2 ) then
          !do spi = 1, numberOfSpecies 
          !  oib = 0 
          !  !! build auxiliary vectors of occupied and virtuals orbitals
          !  do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          !    if ( orbB(spi)%values(pi) == 1 ) then
          !      oib = oib + 1
          !      occB(spi)%values(oib) = pi
          !    end if
          !  enddo
          !enddo
        else 
          cycle 
        endif

        !! same conf 
        if ( sum(couplingS) == 0 ) then
          CIenergy = CISCI_instance%diagonalTarget%values(aa) 

          hamiltonianMatrix%values(a,a) = CIenergy 
        endif 
        !! one orbital different
        if ( sum(couplingS) == 1 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
          CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, diffOrbi(1), diffOrbi(3)  )

          hamiltonianMatrix%values(a,b) = CIenergy * factorA 
        endif
        !! two orbital different, same species
        if ( sum(couplingS) == 2 .and. maxval(couplingS) == 2 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
          CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  )

          hamiltonianMatrix%values(a,b) = CIenergy * factorA 
        endif
        !! two orbital different, different species
        if ( sum(couplingS) == 2 .and. maxval(couplingS) == 1 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) then 
              spi = i
              exit
            end if
          end do

          do i = spi+1, numberOfSpecies
            if ( couplingS(i) == 1 ) spj = i
          end do
          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
          diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj)%values, orbB(spj)%values, occA(spj)%values, occB(spj)%values, factorB )
          CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  )

          hamiltonianMatrix%values(a,b) = CIenergy * factorA * factorB 
        endif
 
      end do bloop !b
    end do aloop !a 
    !$omp end do nowait

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) )
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( orbA(spi) )
      call Vector_destructorInteger ( orbB(spi) )
    end do

    deallocate ( couplingS )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( orbA  )
    deallocate ( orbB  )
    !$omp end parallel

!$  timeB = omp_get_wtime()

!$    write(*,"(T2,A,ES10.2,A4)") "Time for building CI Hamiltonian ",timeB -timeA, " (S)"

  end subroutine CISCI_buildHamiltonian


  subroutine CISCI_buildDiagonal ( diagonal, configurations_orb, spaceSize )
    implicit none
    type(Vector), intent(inout) :: diagonal 
    type (IMatrix1), allocatable, intent(in) :: configurations_orb(:)
    integer(8), intent(in) :: spaceSize
    integer(8) :: a,b,aa,bb
    integer :: i, ii, jj, n, spi, spj
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer(1) :: coupling
    integer :: nproc
    integer :: pi
    integer :: oia
    type (ivector), allocatable :: occA(:)
    type (ivector), allocatable :: orbA(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies

    allocate ( occA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) 
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

!$  timeA= omp_get_wtime()

    do aa = 1, spaceSize 

      !a = CISCI_instance%index_amplitudeCore%values(aa)
      a = aa
      if ( configurations_orb(1)%values(1,a) == -1_1 ) exit

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        occA(spi)%values = 0
        oia = 0 
        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( CISCI_instance%confTarget_orb%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = configurations_orb(spi)%values(:,a)
        !!call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore_orb%values(spi,a), orbA(spi)%values )

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo

      enddo

      diagonal%values(aa) = CISCI_calculateEnergyZero( occA )

    end do !a 

!$  timeB = omp_get_wtime()

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( orbA(spi) ) 
    end do

    deallocate ( occA  )
    deallocate ( orbA  )

    return

  end subroutine CISCI_buildDiagonal

  function CISCI_calculateEnergyOne( si, occA, occB, a, b ) result ( CIenergy )
    implicit none
    type(ivector), intent(in) :: occA(:), occB(:)
    integer(8), intent(in) :: a, b
    integer, intent(in) :: si
    integer :: sj
    integer(8) :: ab, ll, abll, albl, auxab
    integer(8) :: l, la
    real(8) :: CIenergy

    CIenergy = 0.0_8
  
    CIenergy = CIenergy + CIcore_instance%twoCenterIntegrals(si)%values( a, b )

    ab = CIcore_instance%twoIndexArray(si)%values( a, b )
  
    do la = 1_8, CIcore_instance%occupationNumber( si ) !! the same orbitals pair are excluded by the exchange
  
        l = occA(si)%values(la) ! or b, both are the same
  
        ll = CIcore_instance%twoIndexArray(si)%values( l,l ) 
        abll = CIcore_instance%fourIndexArray(si)%values( ab, ll )
  
        CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals(si,si)%values( abll, 1_8)
  
        albl = CIcore_instance%fourIndexArray(si)%values( &
                               CIcore_instance%twoIndexArray(si)%values( a,l ), &
                               CIcore_instance%twoIndexArray(si)%values( l,b ) ) 
  
        CIenergy = CIenergy + MolecularSystem_instance%species(si)%kappa*CIcore_instance%fourCenterIntegrals(si,si)%values(albl, 1_8)
  
      end do
  
      do sj = 1, si - 1 !! avoid ii, same species

        auxab = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( ab - 1_8 )

        do la = 1_8,  CIcore_instance%occupationNumber( sj ) 

          l = occA(sj)%values(la) ! or b, both are the same

          abll = auxab + CIcore_instance%twoIndexArray(sj)%values( l, l)

          CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals( si, sj )%values(abll, 1_8)

        end do

      end do

      do sj = si + 1, MolecularSystem_instance%numberOfQuantumSpecies !! avoid ii, same species

        auxab = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( ab - 1_8)

        do la = 1_8,  CIcore_instance%occupationNumber( sj )

          l = occA(sj)%values(la) ! or b, both are the same

          abll = auxab  + CIcore_instance%twoIndexArray(sj)%values( l, l) 

          CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals( si, sj )%values(abll, 1_8)

        end do

      end do

  end function CISCI_calculateEnergyOne

  function CISCI_calculateEnergyTwoSame( si, occA, occB, ai, aj, bi, bj ) result ( CIenergy )
    implicit none
    type(ivector), intent(in) :: occA(:), occB(:)
    integer(8), intent(in) :: ai, aj, bi, bj
    integer, intent(in) :: si
    integer(8) :: aibi_ajbj, aibj_ajbi
    real(8) :: CIenergy

    CIenergy = 0.0_8

    aibi_ajbj = CIcore_instance%fourIndexArray(si)%values( &
                    CIcore_instance%twoIndexArray(si)%values( ai, bi ), &
                    CIcore_instance%twoIndexArray(si)%values( aj, bj ) )
  
    CIenergy = CIcore_instance%fourCenterIntegrals(si,si)%values( aibi_ajbj, 1_8)

    aibj_ajbi = CIcore_instance%fourIndexArray(si)%values( &
                    CIcore_instance%twoIndexArray(si)%values( ai, bj ), &
                    CIcore_instance%twoIndexArray(si)%values( aj, bi ) )
 
    CIenergy = CIenergy + MolecularSystem_instance%species(si)%kappa * & 
                          CIcore_instance%fourCenterIntegrals(si,si)%values( aibj_ajbi, 1_8)
  
  end function CISCI_calculateEnergyTwoSame

  function CISCI_calculateEnergyTwoDiff( si, sj, ai, aj, bi, bj ) result ( CIenergy )
    implicit none
    integer(8), intent(in) :: ai, aj, bi, bj
    integer, intent(in) :: si, sj
    integer(8) :: aibi,  aux_aibi, ajbj
    real(8) :: CIenergy

    CIenergy = 0.0_8

    aibi = CIcore_instance%twoIndexArray(si)%values( ai, bi )
    aux_aibi = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( aibi - 1_8 )

    ajbj = CIcore_instance%twoIndexArray(sj)%values( aj, bj )
    CIenergy = CIcore_instance%fourCenterIntegrals( si, sj )%values( aux_aibi +  ajbj, 1_8)

  end function CISCI_calculateEnergyTwoDiff

  function CISCI_calculateEnergyZero( occA ) result (CIenergy)
    implicit none

    type(ivector), intent(in) :: occA(:)
    integer :: si,sj
    integer(8) :: ki,li,lj,k,l,kk,ll,kkll,kllk
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    real(8) :: CIenergy

    CIenergy = 0.0_8

    do si = 1, MolecularSystem_instance%numberOfQuantumSpecies
      do ki = 1_8, CIcore_instance%occupationNumber( si )  !! 1 is from a and 2 from b

        k = occA( si )%values( ki )

        !One particle terms
        CIenergy = CIenergy + CIcore_instance%twoCenterIntegrals( si )%values( k, k )

        !Two particles, same specie
        kk = CIcore_instance%twoIndexArray( si )%values( k, k)

        do li = ki + 1_8, CIcore_instance%occupationNumber( si )  !! 1 is from a and 2 from b

          l = occA( si )%values( li )
          ll = CIcore_instance%twoIndexArray( si )%values( l,l )
          kkll = CIcore_instance%fourIndexArray( si )%values( kk, ll ) 

          !Coulomb
          CIenergy = CIenergy + &
              CIcore_instance%fourCenterIntegrals( si, si )%values( kkll, 1_8)

          !Exchange, depends on spin

          kllk = CIcore_instance%fourIndexArray( si )%values( &
                        CIcore_instance%twoIndexArray( si )%values(k,l), &
                        CIcore_instance%twoIndexArray( si )%values(l,k) )

          CIenergy = CIenergy + &
                  MolecularSystem_instance%species( si )%kappa*CIcore_instance%fourCenterIntegrals( si, si )%values( kllk, 1_8)
        end do

        !!Two particles, different species
        do sj = si + 1, MolecularSystem_instance%numberOfQuantumSpecies

          do lj = 1_8, CIcore_instance%occupationNumber( sj ) !! 1 is from a and 2 from b
            l = occA( sj )%values(lj)

            ll = CIcore_instance%twoIndexArray( sj )%values(l,l)
            kkll = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * (kk - 1_8 ) + ll

            CIenergy = CIenergy + &
            CIcore_instance%fourCenterIntegrals( si, sj )%values( kkll, 1_8)

          end do

        end do

      end do
    end do

    CIenergy = CIenergy + HartreeFock_instance%puntualInteractionEnergy

  end function CISCI_calculateEnergyZero

  subroutine CISCI_PT2 ( SCITargetSpaceSize, refEnergy, energyCorrection, eigenVectors )
    implicit none
    integer(8), intent(in) :: SCITargetSpaceSize
    real(8), intent(in) :: refEnergy
    real(8), intent(inout) :: energyCorrection
    type(matrix), intent(in) :: eigenVectors
    real(8) :: CIEnergy
    integer(8) :: nonzero, nonzeroTarget
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b, aa,bb
    integer :: spi, spj, numberOfSpecies
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer(1), allocatable :: couplingS(:)
    integer :: pi
    integer :: oia, oib
    integer :: factorA, factorB
    real(8) :: diagonal, denominator
    real(8) :: energyIncrement, energyIncrement_corrected, energyCorrection_errorCrumbs, energyCorrection_aux ! For Kahan summ
    integer(1) :: coupling
    integer(8) :: diffOrbi(4)
    integer(8) :: diffOrbj(4)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    energyCorrection = 0.0_8

    nonzeroTarget = 0
    do aa = 1, CISCI_instance%targetSpaceSize
      a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      if (CISCI_instance%confTarget_orb(1)%values(1,a) == -1_1 ) exit
      nonzeroTarget = nonzeroTarget + 1
    enddo

    nonzero = 0
    do aa = CISCI_instance%targetSpaceSize + 1, CISCI_instance%buffer_amplitudeCoreSize
      a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      if (CISCI_instance%confAmplitudeCore_orb(1,a) == -1_1 .or. abs(CISCI_instance%buffer_amplitudeCore%values(aa)) <= 1E-10 ) exit
      nonzero = nonzero + 1
    enddo

    write(6,"(T2,A31)") "Computing SCI-PT2 correction..."
    write(6,"(T2,A26,ES10.2,A5,ES10.2,A9,I10)") "Buffer coefficients. Max: ", &
                                       CISCI_instance%buffer_amplitudeCore%values(CISCI_instance%targetSpaceSize + 1), &
                                       " Min: ", CISCI_instance%buffer_amplitudeCore%values(CISCI_instance%targetSpaceSize + nonzero), &
                                       " Nonzero: ", nonzero
    if ( nonzero == 0 ) then
      write(6,"(T2,A44)") "Buffer is empty, skipping SCI-PT2 correction"
      return
    endif

!$  timeA = omp_get_wtime()

    !! generate occupied orbital representation
    call CISCI_orb2occ()

    !$omp parallel &
    !$omp& private(aa, a, spi, oia, orbA, pi, occA, CIenergy, bb, b, oib, orbB, occB, couplings, coupling, i, ii, &
    !$omp&        diagonal, denominator, diffOrbi, diffOrbj, spj, factorA, factorB, energyIncrement, energyIncrement_corrected, energyCorrection_aux ) &
    !$omp& reduction (+:energyCorrection, energyCorrection_errorCrumbs)

    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) 
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    energyCorrection = 0.0_8
    energyCorrection_errorCrumbs = 0.0_8
    
    !!$omp do schedule (runtime) !with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    aloop: do aa = CISCI_instance%targetSpaceSize + 1,  CISCI_instance%targetSpaceSize + nonzero

      a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      !a = aa ! if index_amplitude is sorted

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        orbA(spi)%values(:) = CISCI_instance%confAmplitudeCore_orb(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo

      enddo

      CIenergy = 0.0_8
      bloop: do bb = 1,  nonzeroTarget

        !b = CISCI_instance%index_amplitudeCore%values(bb)
        b = bb

        ! getting configuration B
        do spi = 1, numberOfSpecies 

          orbB(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,b)
          occB(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:, b) 

        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - dot_product ( orbA(spi)%values(:), orbB(spi)%values(:) ) 
        end do

        select case ( sum(couplingS) )
    
        !! one orbital different
        case (1)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
          CIenergy = CIenergy + CISCI_calculateEnergyOne( spi, occA, occB, diffOrbi(1), diffOrbi(3)  ) * factorA * eigenVectors%values(bb,1)

        !! two orbital different
        case(2)

          select case ( maxval(couplingS) )

          !! two orbital different, same species
          case (2)
            do i = 1, numberOfSpecies
              if ( couplingS(i) == 2 ) spi = i
            end do

            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
            CIenergy = CIenergy + CISCI_calculateEnergyTwoSame( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  ) * factorA * eigenVectors%values(bb,1)
          !! two orbital different, different species
          case (1)
            do i = 1, numberOfSpecies
              if ( couplingS(i) == 1 ) then 
                spi = i
                exit
              end if
            end do
            do i = spi+1, numberOfSpecies
              if ( couplingS(i) == 1 ) spj = i
            end do

            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
            diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj)%values, orbB(spj)%values, occA(spj)%values, occB(spj)%values, factorB )
            CIenergy = CIenergy + CISCI_calculateEnergyTwoDiff( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  ) * factorA * factorB * eigenVectors%values(bb,1)
          end select ! maxval(couplingS)  
        end select ! sum(couplingS) 
 
      end do bloop !b

      !! calculate diagonal term and denominator of Eq5 10.1063/1.4955109
      diagonal = CISCI_calculateEnergyZero( occA )
      denominator = 1 / ( refEnergy - diagonal ) 
      !! PT2 correction, normal way
      energyCorrection = energyCorrection + ( CIenergy**2) * denominator
      !! Kahan summation way. serial or omp static only
      !energyIncrement = (CIenergy**2) * denominator
      !energyIncrement_corrected = energyIncrement - energyCorrection_errorCrumbs
      !energyCorrection_aux = energyCorrection + energyIncrement_corrected
      !energyCorrection_errorCrumbs = ( energyCorrection_aux - energyCorrection ) - energyIncrement_corrected
      !energyCorrection = energyCorrection_aux

    end do aloop !a 
   !$omp end do nowait

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( couplingS )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( orbA  )
    deallocate ( orbB  )
    !$omp end parallel

    !energyCorrection = energyCorrection - energyCorrection_errorCrumbs

!$  timeB = omp_get_wtime()
    write (6,"(T2,A,F25.12,A,ES10.2)") "CI-PT2 energy correction: ", energyCorrection
    !write (6,"(T2,A,F25.12,A,ES10.2)") "CI-PT2 energy correction: ", energyCorrection, " Kahan's error crumbs: ", energyCorrection_errorCrumbs
!$  write(*,"(A,ES10.2)") "Time for CI-PT2 correction: ", timeB -timeA

  end subroutine CISCI_PT2

  subroutine CISCI_saveEigenVector ( eigenVectors )
    implicit none
    type(matrix), intent(in) :: eigenVectors
    character(50) :: nameFile
    integer :: unitFile
    real(8) :: timeA, timeB
    integer(8) :: a
    integer :: spi, spj, numberOfSpecies
    type (ivector), allocatable :: orbA(:), orbReF(:)
    integer, allocatable :: CIlevel(:)

!$  timeA = omp_get_wtime()

    numberOfSpecies = CIcore_instance%numberOfSpecies 

    unitFile = 36 
    nameFile = trim(CONTROL_instance%INPUT_FILE)//"sci"
    open(unit = unitFile, file=trim(nameFile), status="new", form="formatted")

    allocate ( CIlevel ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbRef ( numberOfSpecies ) )

    CIlevel = 0

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbRef(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    !building reference orbitals
    do spi = 1, numberOfSpecies
      orbRef(spi)%values(1:CIcore_instance%numberOfOccupiedOrbitals%values(spi) ) = 1
    enddo 

    do a = 1, CISCI_instance%targetSpaceSize 
      ! getting configuration A
      do spi = 1, numberOfSpecies 
        !! build the orbital from the index using the bit mapping
        orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:, a) 

        CIlevel(spi) = CIcore_instance%numberOfOccupiedOrbitals%values(spi) - sum( orbA(spi)%values * orbRef(spi)%values ) 
      enddo

      write(unitFile,*) a, eigenVectors%values(a,1), CIlevel(:), sum(CIlevel(:))

    end do 

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( orbRef(spi) ) 
      call Vector_destructorInteger ( orbA(spi) ) 
    end do

    deallocate ( orbRef )
    deallocate ( orbA )
    deallocate ( CIlevel )
    close(unitFile)

!$  timeB = omp_get_wtime()
!$  write(*,"(A,ES10.2)") "Time for saving SCI eigenVector: ", timeB -timeA

  end subroutine CISCI_saveEigenVector

  !> Vector is defined in reverse order, e.g. vector (1,1,0), binary : 011, decimal : 3
  subroutine CISCI_binaryToDecimal ( binary, decimalNumber )
    implicit none
    integer, intent(in) :: binary(:)
    real(8), intent(out) :: decimalNumber
    integer :: n

    decimalNumber = 0.0_8

    do n = 1, size(binary) 
      decimalNumber = decimalNumber + binary(n) * ( 2.0_8**(n-1) )
    enddo

  end subroutine CISCI_binaryToDecimal 

  !> Vector is defined in reverse order, e.g. vector (1,1,0), binary : 011, decimal : 3
  subroutine CISCI_decimalToBinary ( decimalNumber, binary )
    implicit none
    real(8), intent(in) :: decimalNumber
    integer, intent(out) :: binary(:)
    real(8) :: auxdecimalNumber
    integer :: n

    binary(:) = 0
    auxdecimalNumber = decimalNumber
    n = 0

    do while ( auxdecimalNumber > 0 ) 
      n = n + 1
      binary(n) = modulo(auxdecimalNumber,2.0_8 )
      auxdecimalNumber = floor( auxdecimalNumber/2.0_8, SELECTED_INT_KIND(16))
    enddo 

  end subroutine CISCI_decimalToBinary 

  !! calculate the sign factor for canonical order of the configuration
  function CISCI_canonicalOrderFactor ( spi, orb, occ ) result (factor) 
    implicit none
    integer, intent(in) :: spi
    type(Ivector), intent(in) :: orb
    type(IVector), intent(in) :: occ
    integer :: factor 
    integer :: n_permu, pi, oi, qi
    integer :: above, below, pos

    ! counter for positions above or below the canonical order
    above = 0
    below = 0  
    qi = 0
    do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
      pos = 0
      qi = qi + orb%values(pi)
      if ( orb%values(pi) == 1 ) then
        do oi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          if ( pi == occ%values(oi) ) then
            pos = oi
            exit
          endif
        enddo
      endif 

      if ( qi * orb%values(pi) - pos > 0 ) above = above + 1 
      if ( qi * orb%values(pi) - pos < 0 ) below = below + 1 
    enddo

    n_permu = max(above,below)
    factor = 1
    if (mod(n_permu, 2) /= 0) factor = -1

  end function CISCI_canonicalOrderFactor 

  !! calculate the sign factor for canonical order of the configuration
  function CISCI_getDiffOrbitals ( spi, orbA, orbB, occA, occB, factor ) result (diffOrb) 
    implicit none
    integer, intent(in) :: spi
    integer, intent(in), contiguous :: orbA(:), orbB(:), occA(:), occB(:)
    integer, intent(out) :: factor
    integer :: diffOrb(4), diffPos(4)
    integer :: pi, z
    integer :: phase_exponent
    integer :: n_occ

    n_occ = CIcore_instance%numberOfOccupiedOrbitals%values(spi)

    diffPos = 0
    diffOrb = 0
    z = 0
    ! different orbital in A
    do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, n_occ
      if ( orbB(occA(pi) ) == 0  ) then
        z = z + 1
        diffOrb(z) = occA(pi)
        diffPos(z) = pi
      endif  
    enddo

    z = 2
    ! different orbital in B
    do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, n_occ
      if ( orbA(occB(pi) ) == 0  ) then
        z = z + 1
        diffOrb(z) = occB(pi)
        diffPos(z) = pi
      endif  
    enddo

    !factor = (-1)**(diffPos(1)-diffPos(3) + diffPos(2) - diffPos(4) )
    phase_exponent = diffPos(1) - diffPos(3) + diffPos(2) - diffPos(4)

    if (iand(phase_exponent, 1) == 0) then
      factor = 1.0
    else
      factor = -1.0
    end if

  end function CISCI_getDiffOrbitals

  subroutine CISCI_appendAmplitude ( n, amplitude, orbB )
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: amplitude
    type (ivector) :: orbB(:)
    integer(8) :: m
    
    real(8) :: diffConf
    logical :: addConf
    integer :: spi, numberOfSpecies

    numberOfSpecies = CIcore_instance%numberOfSpecies 

    !! add the configuration at the end of the array
    CISCI_instance%omp_target_iterator_m(n) = CISCI_instance%omp_target_iterator_m(n) + 1
    m = CISCI_instance%omp_target_iterator_m(n) 

    CISCI_instance%buffer_amplitudeCore%values(m) = amplitude
    do spi = 1, numberOfSpecies
      CISCI_instance%confAmplitudeCore_orb(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), m)  = orbB(spi)%values(:)
    enddo

    !! sort the array when it's full
    if ( m == CISCI_instance%omp_targetInterval(2, n ) ) then
      call CISCI_sortAmplitude( n, auxm = m )
      !! update the iterator after cleaning
      CISCI_instance%omp_target_iterator_m(n) = m 
    endif

  end subroutine CISCI_appendAmplitude

  subroutine CISCI_sortAmplitude( n, auxm )
    implicit none
    integer, intent(in) :: n
    integer(8), intent(inout), optional  :: auxm
    integer(8)  :: m1, m2, m2b
    real(8) :: diffConf
    integer(8) :: halfm
    integer :: i, j, ii, spi
    logical :: addConf
    integer :: nproc

    !! determine the chunk to sort
    if ( n ==  CIcore_instance%nproc + 1 ) then !! sort the full array
      m1 = 1
      m2 = CISCI_instance%buffer_amplitudeCoreSize 
      m2b = CISCI_instance%buffer_amplitudeCoreSize 

      !! reset auxindex array to the full size, global absolute position
      do i = m1, m2
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo

    else
      !! the subarray is full, time to sort, merge and clean 
      m1 = CISCI_instance%omp_targetInterval(1, n ) 
      m2 = CISCI_instance%omp_targetInterval(2, n ) 
      m2b = CISCI_instance%omp_targetInterval(2, n )  !! to avoid zeros
      do i = m1, CISCI_instance%omp_targetInterval(2, n ) 
        if (  CISCI_instance%confAmplitudeCore_orb(1,i) == -1_1 ) then
          m2b = i-1
          exit
        endif
      enddo
    endif 

    !! Sort according to the index of CI configurations per species in order to find duplicates ( N log( N ) )
    call CISort_quicksort_matrix(  CISCI_instance%confAmplitudeCore_orb(:,m1:m2), &
                                     CISCI_instance%index_amplitudeCore%values(m1:m2), & 
                                     1_8, m2b - m1 + 1, n )

    !! Sort amplituted coeff vector according to sorting of index array, keeping both arrays aligned
    call CISort_sortVectorByIndex( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                                   CISCI_instance%index_amplitudeCore%values(m1:m2), &
                                   m2 - m1 + 1 )

    !! merge duplicated configurations (sum amplitudes)
    call CISort_mergeDuplicates ( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                                           CISCI_instance%confAmplitudeCore_orb(:,m1:m2), &
                                           CISCI_instance%index_amplitudeCore%values(m1:m2), &
                                           m2 - m1 + 1 )

    !! reset auxindex arrary, relative position
    do i = 1, m2 - m1 + 1
      CISCI_instance%index_amplitudeCore%values( m1 + i - 1 ) = i
    enddo

    !call MTSort ( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
    !              CISCI_instance%index_amplitudeCore%values(m1:m2), &
    !              m2 - m1 + 1, "D", 1 )

    !! descening sort according to the absolute value of the amplitude
    call CISort_quicksort_vector(  CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                                     CISCI_instance%index_amplitudeCore%values(m1:m2), & 
                                     1_8, m2b - m1 + 1 )
    
    !! organize configurations according to the sorted amplitudes
    call CISort_sortArrayByIndex( CISCI_instance%confAmplitudeCore_orb(:,m1:m2), &
                                  CISCI_instance%index_amplitudeCore%values(m1:m2), &
                                  CISort_instance%combinedNumberOfOrbitals, &
                                  m2 - m1 + 1, n )

    !! reset auxindex arrary, relative positions
    do i = 1, m2 - m1 + 1
      CISCI_instance%index_amplitudeCore%values( m1 + i - 1 ) = i
    enddo

    !! cleaning and finding limit
    if ( present ( auxm ) ) then

      !! get the new end position
      halfm = ((m2 - m1 + 1)/2) + m1 - 1 ! half point of the array
      auxm = halfm ! new limit of the array, new elements will be stored after here

      !! check if there are any zeros in the first half
      do i = m1, halfm +  1
        if ( CISCI_instance%confAmplitudeCore_orb(1,i) == -1_1 ) then
          auxm = i - 1
          exit
        endif
      enddo

      !! discard the last two quarters of tmp_ampltitude for next run, if not keep it fot PT2 corr
      !CISCI_instance%buffer_amplitudeCore%values( halfm + 1 : m2 ) = 0.0_8
      !do spi = 1, CIcore_instance%numberOfSpecies 
      !  CISCI_instance%confAmplitudeCore_orb(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), &
      !                                   halfm + 1 : m2) = -1_1
      ! 
      !enddo
    endif

  end subroutine CISCI_sortAmplitude

  subroutine CISCI_mergeCoreAndTarget( confTarget_orb, eigenVectors )
    implicit none
    type (IMatrix1), intent(inout) :: confTarget_orb(:)
    type (matrix), intent(inout) :: eigenVectors
    integer :: spi
    integer(8) :: i, j, m
    integer :: orb, auxorb
    logical :: is_equal
    real(8) :: timeA, timeB

!$  timeA = omp_get_wtime()

    !! first add the core to the new target space
    m = 0
    do j = 1, CISCI_instance%coreSpaceSize
      if ( CISCI_instance%confCore(1)%values(1,j) == -1_1 ) exit
      m = m + 1
      do spi = 1, CIcore_instance%numberOfSpecies 
        confTarget_orb(spi)%values(:,m) = CISCI_instance%confCore(spi)%values(:,j)
      enddo
      eigenVectors%values(m,1) = CISCI_instance%coefficientCore%values(j)
    enddo

    if ( m == CISCI_instance%targetSpaceSize ) return  ! the targetSpace is full

    !! run over top amplitudes to add to the new target space
    targetSpace: do i = 1, CISCI_instance%targetSpaceSize
      if ( CISCI_instance%confAmplitudeCore_orb(1, i) == -1 ) exit targetSpace
      !! but check if the conf is already included in the core space (previously added)
      is_equal = .false.
      coreSpace : do j = 1, CISCI_instance%coreSpaceSize

        auxorb = 0

        !! compare each orb for all species
        species: do spi = 1, CIcore_instance%numberOfSpecies 

          do orb = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
            auxorb = auxorb + 1
            if (  CISCI_instance%confCore(spi)%values(orb,j) == CISCI_instance%confAmplitudeCore_orb(auxorb, i) ) then
              is_equal = .true.
            else
              is_equal = .false.
              exit species !! this two conf are diff, so go to the next
            endif
          enddo
        enddo species

        if ( is_equal ) exit coreSpace ! this conf is already included

      enddo coreSpace

      !! if the conf is not in core space then add the configuration
      if ( .not. is_equal) then
        m = m + 1
        if ( m > CISCI_instance%targetSpaceSize ) exit targetSpace ! the targetSpace is full

        do spi = 1, CIcore_instance%numberOfSpecies 
          confTarget_orb(spi)%values(:,m) = CISCI_instance%confAmplitudeCore_orb(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), i) 
        enddo
        eigenVectors%values(m,1) = CISCI_instance%buffer_amplitudeCore%values(i)

      endif

    enddo targetSpace

!$  timeB = omp_get_wtime()
!$  write(*,"(A,ES10.2,A4)") "** TOTAL Elapsed Time for merging core and new amplitudes : ", timeB - timeA ," (s)"

  end subroutine CISCI_mergeCoreAndTarget

  !! Generate the occupied configuration representation from the orbital configuration (1 or 0)
  subroutine CISCI_orb2occ()
    implicit none
    integer :: spi, numberOfSpecies
    integer :: pi, oia
    integer :: a, aa
    type (ivector), allocatable :: occA(:)
    type (ivector), allocatable :: orbA(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    do spi = 1, numberOfSpecies
      CISCI_instance%confTarget_occ(spi)%values = -1_4
    enddo

    allocate ( occA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) 
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    do a = 1, CISCI_instance%targetSpaceSize

      if (CISCI_instance%confTarget_orb(1)%values(1,a) == -1_1 ) exit

      do spi = 1, numberOfSpecies 

        oia = 0

        !! Stored code for bit-masking representation
        !! transforming from decimal to binary once
        !!do aa = 1, nonzero
        !!  do spi = 1, numberOfSpecies
        !!    !call CISCI_decimalToBinary ( CISCI_instance%confTarget_orb%values(spi,aa), CISCI_instance%targetOrb(spi,aa)%values )
        !!    CISCI_instance%targetOrb(spi,aa)%values = CISCI_instance%confTarget_orb(spi)%values(:,aa)
        !!  enddo
        !!enddo

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( CISCI_instance%confTarget_orb%values(spi,b), orbB(spi)%values )
        !!orbA(spi)%values = CISCI_instance%targetOrb(spi,a)%values
        orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:, a) 

        occA(spi)%values(:) = 0
        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1_1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo
        !!occB(spi)%values(:) = pack( CISCI_instance%canonicalOrder(spi)%values, orbB(spi)%values(:) == 1 )

        if ( orbA(spi)%values(1) == -1_1 ) occA(spi)%values(:) = -1_4

        CISCI_instance%confTarget_occ(spi)%values(:, a) = occA(spi)%values(:) 

      enddo ! spi, numberOfSpecies
    enddo ! a, buffer_amplitudeCoreSize

    deallocate ( occA )
    deallocate ( orbA )

  end subroutine CISCI_orb2occ

  !! Generate the orbital configuration representation from the occupied orbitals
  subroutine CISCI_occ2orb()
    implicit none
    integer :: spi, numberOfSpecies
    integer :: pi
    integer :: a, aa

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies

    do spi = 1, numberOfSpecies
      CISCI_instance%confTarget_orb(spi)%values = -1_1
    enddo

    do a = 1,  CISCI_instance%targetSpaceSize

      aa = CISCI_instance%index_amplitudeCore%values(a)
      if (CISCI_instance%confTarget_occ(1)%values(1,aa) == -1_4 ) exit

      do spi = 1, numberOfSpecies
        CISCI_instance%confTarget_orb(spi)%values(:,a) = 0_1
        do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          CISCI_instance%confTarget_orb(spi)%values( CISCI_instance%confTarget_occ(spi)%values(pi,aa), a ) = 1_4
        enddo
      enddo
    enddo

  end subroutine CISCI_occ2orb

  !! Generate a list of all possible species combinations for CI
  subroutine CISCI_buildCIOrderList( maxCIexcitations, CIorder_list, CIorder_count, CIorder_weight )
    implicit none
    integer, allocatable, intent(inout) :: maxCIexcitations(:)
    integer, allocatable, intent(inout) :: CIorder_list(:,:)
    integer, allocatable, intent(inout) :: CIorder_count(:)
    real(8), allocatable, intent(inout) :: CIorder_weight(:)
    integer :: spi, numberOfSpecies
    integer :: activeOrbitals
    integer, allocatable :: CIlevel(:)
    integer :: m1, m2
    integer :: i, j
    logical :: done
    integer :: combination, totalCombinations

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies
    allocate ( CIlevel( numberOfSpecies )) 
    allocate ( maxCIexcitations( numberOfSpecies )) 
    CIlevel = 0

    totalCombinations = 1
    do spi = 1, numberOfSpecies 
      maxCIexcitations( spi ) = CIcore_instance%numberOfOccupiedOrbitals%values(spi) - CIcore_instance%numberOfCoreOrbitals%values(spi)  
      totalCombinations = totalCombinations * ( maxCIexcitations( spi ) + 1 )
    enddo

    allocate ( CIorder_list( numberOfSpecies, totalCombinations )) 
    allocate ( CIorder_count( totalCombinations )) 
    allocate ( CIorder_weight( totalCombinations )) 
    CIorder_list = 0
    CIorder_count = 0
    CIorder_weight = 0.0_8

    combination = 0
    done = .false. 
    do while (.not. done)
      combination = combination + 1
      CIorder_list(:, combination) = CIlevel(:)
      !! increment the configuration from right to left 
      do j = numberOfSpecies, 1, -1
        if ( CIlevel(j) < maxCIexcitations(j) ) then 
          CIlevel(j) = CIlevel(j) + 1
          exit  
        else
          CIlevel(j) = 0  ! reset current position and carry over to the left
        end if

        if (j == 1) done = .true.
      end do
    end do

    deallocate ( CIlevel ) 

    write (6, "(T2,A)") "--------------------------"
    write (6, "(T2,A)") "ID | CI level per species "
    write (6, "(T2,A)") "--------------------------"
    do j = 1, size( CIorder_list, dim = 2 )
      write (6, "(T2,I2)", advance="no") j
      write (6, "(A2)", advance="no") " |"
      do spi = 1, numberOfSpecies
        write (6, "(T2,I4)", advance="no") CIorder_list(spi, j)
      end do
      write (6, "(A)") ""
    end do
    write (6, "(T2,A)") "--------------------------"

  end subroutine CISCI_buildCIOrderList

  !! Given a CI excitation level per species, return the position in the list of all CI combinations
  function CISCI_combinationIndex(CIlevel, maxCIexcitations ) result(idx)
    implicit none
    integer, dimension(:), intent(in) :: CIlevel, maxCIexcitations
    integer :: idx
    integer :: numberOfSpecies, i, j, multiplier

    numberOfSpecies = size(CIlevel)
    idx = 1  ! Start with 1 for Fortran 1-based indexing

    do i = 1, numberOfSpecies
       ! Calculate the product of bases for all species to the right of j
       multiplier = 1
       do j = i + 1, numberOfSpecies
          multiplier = multiplier * ( maxCIexcitations(j) + 1 )
       end do
       
       ! Add contribution of current species
       idx = idx + CIlevel(i) * multiplier
    end do

  end function CISCI_combinationIndex

  !! count number of configurations in the buffer space per species subspaces
  subroutine CISCI_countSpeciesPairs()
    implicit none
    integer :: spi, numberOfSpecies
    integer :: activeOrbitals
    integer, allocatable :: CIlevel(:)
    integer :: m1, m2
    integer :: a, c
    logical :: done
    integer :: CIorder_index

    write (6, "(T2,A)") "Counting number of configurations in the buffer space..."

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies
    allocate ( CIlevel( numberOfSpecies )) 
    CIlevel = 0

    do a = 1, CISCI_instance%buffer_amplitudeCoreSize
      if ( CISCI_instance%confAmplitudeCore_orb(1, a) == -1 ) exit 

      !! compare each orb for all species
      do spi = 1, numberOfSpecies 
        m1 = CISCI_instance%combinedOrbitalsPositions(1,spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi)
        m2 = CISCI_instance%combinedOrbitalsPositions(2,spi)
    
        CIlevel(spi) = sum( CISCI_instance%confAmplitudeCore_orb(m1:m2, a) )
      enddo

      CIorder_index = CISCI_combinationIndex(CIlevel, CISCI_instance%maxCIexcitations )
      CISCI_instance%CIorder_count(CIorder_index) = CISCI_instance%CIorder_count(CIorder_index) + 1
      !print *, i, CIlevel, CISCI_combinationIndex(CIlevel, CISCI_instance%maxCIexcitations )

    enddo

    deallocate ( CIlevel ) 

    write (6, "(T2,A)") "--------------------------------------------------"
    write (6, "(T2,A)") "ID | CI_tot | # config | CI level per species     "
    write (6, "(T2,A)") "--------------------------------------------------"
    do c = 1, size( CISCI_instance%CIorder_list, dim = 2 )
      write (6, "(T2,I2,A2)", advance="no") c, " |"
      write (6, "(T2,I6)", advance="no") sum(CISCI_instance%CIorder_list(:, c))
      write (6, "(A2)", advance="no") " |"
      write (6, "(T2,I8)", advance="no") CISCI_instance%CIorder_count(c) 
      write (6, "(A2)", advance="no") " |"
      do spi = 1, numberOfSpecies
        write (6, "(T2,I4)", advance="no") CISCI_instance%CIorder_list(spi, c)
      end do
      write (6, "(A)") ""
    end do
    write (6, "(T2,A)") "--------------------------------------------------"

  end subroutine CISCI_countSpeciesPairs

  !! Sort the two particles contributions according to HeatBath CI method
  !! Section II.A of 10.1021/acs.jctc.6b00407
  !! TODO: check where we need total or active number of orbitals
  subroutine CISCI_heatbathIntegralSorting( heatBathDoubleExcitations, heatBathDoubleExcitations_index, &
    heatBathDoubleExcitations_size )
    implicit none
    type(Matrix), allocatable, intent(out) :: heatBathDoubleExcitations(:,:)
    type(IMatrix8), allocatable, intent(out) :: heatBathDoubleExcitations_index(:,:)
    type(IMatrix8), allocatable, intent(out) :: heatBathDoubleExcitations_size(:,:)
    integer :: spi, spj, numberOfSpecies
    integer(8) :: pi, qi, ri, si
    integer(8) :: pj, qj, rj, sj
    integer(8) :: pq, rs, pr, ps, rq, qs, pqrs, psrq, prqs
    integer(8) :: size_ii, size_ij
    real(8) :: kappa
    logical, allocatable :: occupied(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies

    allocate ( occupied ( maxval(CIcore_instance%numberOfActiveOrbitals%values(:) ) ) )
    occupied = .False. 

    !! matrix allocation and initialization
    allocate ( heatBathDoubleExcitations(numberOfSpecies,numberOfSpecies))
    do spi = 1, numberOfSpecies

      !! triangular with diagonal terms
      size_ii = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * ( CIcore_instance%numberOfActiveOrbitals%values(spi) + 1_8)) / 2.0
      call Matrix_constructor ( heatBathDoubleExcitations(spi,spi), &
                                size_ii, size_ii, 0.0_8 &
                              )
      do spj = spi + 1, numberOfSpecies
        !! full
        size_ij = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * CIcore_instance%numberOfActiveOrbitals%values(spj) )
        call Matrix_constructor ( heatBathDoubleExcitations(spj,spi), &
                                  size_ij, size_ij, 0.0_8 &
                                )
      enddo !spj
    enddo !spi

    allocate ( heatBathDoubleExcitations_index(numberOfSpecies,numberOfSpecies))
    do spi = 1, numberOfSpecies
      !! triangular with diagonal terms
      size_ii = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * ( CIcore_instance%numberOfActiveOrbitals%values(spi) + 1_8)) / 2.0
      call Matrix_constructorInteger8 ( heatBathDoubleExcitations_index(spi,spi), &
                                        size_ii, size_ii, 0_8 &
                                      )
      do spj = spi + 1, numberOfSpecies
        !! full
        size_ij = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * CIcore_instance%numberOfActiveOrbitals%values(spj) )
        call Matrix_constructorInteger8 ( heatBathDoubleExcitations_index(spj,spi), &
                                          size_ij, size_ij, 0_8 &
                                        )
      enddo !spj
    enddo !spi

    allocate ( heatBathDoubleExcitations_size(numberOfSpecies,numberOfSpecies))
    do spi = 1, numberOfSpecies
      !! triangular with diagonal terms
      size_ii = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * ( CIcore_instance%numberOfActiveOrbitals%values(spi) + 1_8)) / 2.0
      call Matrix_constructorInteger8 ( heatBathDoubleExcitations_size(spi,spi), &
                                        size_ii, 2_8, 0_8 &
                                      ) ! Max size for SCI and PT heatbath thresholds
      do spj = spi + 1, numberOfSpecies
        !! full
        size_ij = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * CIcore_instance%numberOfActiveOrbitals%values(spj) )
        call Matrix_constructorInteger8 ( heatBathDoubleExcitations_size(spj,spi), &
                                          size_ij, 2_8, 0_8 &
                                        )  ! Max size for SCI and PT heatbath thresholds
      enddo !spj
    enddo !spi

    !! index initialization (same species), including diag terms
    do spi = 1, numberOfSpecies
      do pi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
        do qi = pi, CIcore_instance%numberOfActiveOrbitals%values(spi)
          pq = CIcore_instance%twoIndexArray(spi)%values(qi,pi)
          do ri = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
            do si = ri, CIcore_instance%numberOfActiveOrbitals%values(spi)
              rs = CIcore_instance%twoIndexArray(spi)%values(si,ri)
              !! missing condition for rs > pq? nope
              heatBathDoubleExcitations_index(spi,spi)%values(rs,pq) = rs
            enddo !si
          enddo !ri
        enddo !qi
      enddo !pi
    enddo !spi

    !! getting contributions to double excitaions (same species), excluding diagonal terms (pp -> rr)
    do spi = 1, numberOfSpecies

      size_ii = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * ( CIcore_instance%numberOfActiveOrbitals%values(spi) + 1_8)) / 2.0

      kappa = MolecularSystem_instance%species(spi)%kappa
      do pi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
        occupied(pi) = .True.

        do qi = pi + 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
          occupied(qi) = .True.

          pq = CIcore_instance%twoIndexArray(spi)%values(qi,pi)
          do ri = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
            if ( occupied(ri) ) cycle
            do si = ri + 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
              if ( occupied(si) ) cycle
              rs = CIcore_instance%twoIndexArray(spi)%values(si,ri)

              ps = CIcore_instance%twoIndexArray(spi)%values(si,pi)
              rq = CIcore_instance%twoIndexArray(spi)%values(qi,ri)

              pqrs = CIcore_instance%fourIndexArray(spi)%values(rs,pq)
              psrq = CIcore_instance%fourIndexArray(spi)%values(rq,ps)

              heatBathDoubleExcitations(spi,spi)%values(rs,pq) = CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs, 1_8)
              heatBathDoubleExcitations(spi,spi)%values(rs,pq) = heatBathDoubleExcitations(spi,spi)%values(rs,pq) + kappa * &
              CIcore_instance%fourCenterIntegrals(spi,spi)%values(psrq, 1_8)

              !print *, ri, si, rs, heatBathDoubleExcitations(spi,spi)%values(rs,pq)

            enddo !si
          enddo !ri

          !! removing self excitation term (pq -> pq ) I think this will never happend with the above loops, but
          heatBathDoubleExcitations(spi,spi)%values(pq,pq) = 0.0_8

          !print *, "========"

          !! sort for all rs
          call CISort_quicksort_vector( heatBathDoubleExcitations(spi,spi)%values(:,pq), &
                                        heatBathDoubleExcitations_index(spi,spi)%values(:,pq), &
                                        1_8, size_ii &
                                      )

          !do rs = 1_8, size_ii
          !  print *, rs, heatBathDoubleExcitations_index(spi,spi)%values(rs,pq), &
          !  heatBathDoubleExcitations(spi,spi)%values(rs,pq), "|", &
          !  IndexMap_vectorToMatrix( heatBathDoubleExcitations_index(spi,spi)%values(rs,pq), CIcore_instance%numberOfActiveOrbitals%values(spi))
          !enddo

          !! finding the size of elements above SCI variational threshold
          do rs = 1_8, size_ii
            if ( abs(heatBathDoubleExcitations(spi,spi)%values(rs,pq)) <= CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(1) ) then
              heatBathDoubleExcitations_size(spi,spi)%values(pq,1) = rs - 1_8
              exit
            endif
          enddo

          !! finding the size of elements above PT2 threshold
          do rs = heatBathDoubleExcitations_size(spi,spi)%values(pq,1) + 1_8, size_ii
            if ( abs(heatBathDoubleExcitations(spi,spi)%values(rs,pq)) <= CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(2) ) then
              heatBathDoubleExcitations_size(spi,spi)%values(pq,2) = rs - 1_8
              exit
            endif
          enddo

          !print *, "min pos", heatBathDoubleExcitations_size(spi,spi)%values(pq,1), heatBathDoubleExcitations_size(spi,spi)%values(pq,2)

          occupied(qi) = .False.
        enddo !qi
        occupied(pi) = .False.
      enddo !pi
    enddo !spi

    !! index initialization (diff species)
    do spi = 1, numberOfSpecies
      do spj = spi + 1, numberOfSpecies
        do pi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
          do qj = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)

            !!index for HBCI
            pq = qj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( pi - 1_8 )

            do ri = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
              do sj = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)

                !!index for HBCI
                rs = sj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( ri - 1_8 )
                heatBathDoubleExcitations_index(spj,spi)%values(rs,pq) = rs

              enddo !sj
            enddo !ri
          enddo !qj
        enddo !pi
      enddo !spj
    enddo !spi

    !! getting contributions to double excitaions (diff species)
    do spi = 1, numberOfSpecies
      do spj = spi + 1, numberOfSpecies

        size_ij = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * CIcore_instance%numberOfActiveOrbitals%values(spj) )
        do pi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
          do qj = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)
            !!index for HBCI
            pq = qj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( pi - 1_8 )
            do ri = pi + 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
              !! integral index
              pr = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( CIcore_instance%twoIndexArray(spi)%values(ri,pi) - 1_8 ) !! aux index
              do sj = qj + 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)

                !!index for HBCI
                rs = sj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( ri - 1_8 )

                !! integral index
                qs = CIcore_instance%twoIndexArray(spj)%values(sj,qj)
                !! integral index
                prqs = pr + qs

                heatBathDoubleExcitations(spj,spi)%values(rs,pq) = CIcore_instance%fourCenterIntegrals(spi,spj)%values(prqs, 1_8)
                !print *, spi, spj, ri, sj, rs, heatBathDoubleExcitations(spj,spi)%values(rs,pq)

              enddo !sj
            enddo !ri

            !print *, "================"
            !! removing self excitation term (pq -> pq)
            heatBathDoubleExcitations(spj,spi)%values(pq,pq) = 0.0_8

            !! sort for all rs
            call CISort_quicksort_vector( heatBathDoubleExcitations(spj,spi)%values(:,pq), &
                                          heatBathDoubleExcitations_index(spj,spi)%values(:,pq), &
                                          1_8, size_ij &
                                        )

            !do ri = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi)
            !  do sj = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)

            !    rs = sj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( ri - 1_8 )
            !    print *, spi, spj, ri, sj, rs, heatBathDoubleExcitations_index(spj,spi)%values(rs,pq), &
            !    heatBathDoubleExcitations(spj,spi)%values(rs,pq), "|", &
            !    ((heatBathDoubleExcitations_index(spj,spi)%values(rs,pq) - 1_8 ) / CIcore_instance%numberOfActiveOrbitals%values(spj)) + 1_8, &
            !    mod(heatBathDoubleExcitations_index(spj,spi)%values(rs,pq) - 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)) + 1_8
            !  enddo ! sj
            !enddo ! ri

            !! finding the size of elements above SCI variational threshold
            do rs = 1_8, size_ij
              if ( abs(heatBathDoubleExcitations(spj,spi)%values(rs,pq)) <= CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(1) ) then
                heatBathDoubleExcitations_size(spj,spi)%values(pq,1) = rs - 1_8
                exit
              endif
            enddo

            !! finding the size of elements above PT2 threshold
            do rs = heatBathDoubleExcitations_size(spj,spi)%values(pq,1) + 1_8, size_ij
              if ( abs(heatBathDoubleExcitations(spj,spi)%values(rs,pq)) <= CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(2) ) then
                heatBathDoubleExcitations_size(spj,spi)%values(pq,2) = rs - 1_8
                exit
              endif
            enddo !rs

            !print *, "min pos", heatBathDoubleExcitations_size(spj,spi)%values(pq,1), heatBathDoubleExcitations_size(spj,spi)%values(pq,2)

          enddo ! qj
        enddo ! pi

      enddo ! spj
    enddo !spi

    deallocate ( occupied )

  end subroutine CISCI_heatbathIntegralSorting

  !! generate the configurations to form the target space from the core space
  !! TODO: check where we need total or active number of orbitals
  subroutine CISCI_heatBathGenerate (  diagonal, coefficientCore, confCore, SCICoreSpaceSize, oldEnergy, mode )

    implicit none
    type(Vector), intent(in) :: diagonal
    real(8), intent(in) :: coefficientCore ( SCICoreSpaceSize )
    type(IMatrix1), intent(in) :: confCore(:)
    integer(8), intent(in) :: SCICoreSpaceSize
    real(8), intent(in) :: oldEnergy
    integer, intent(in) :: mode
    real(8) :: CIEnergy
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b,c, aa
    integer :: spi, spj, numberOfSpecies
    integer(8), allocatable :: indexConfA(:) !! ncore, species
    real(8) :: diagEnergy, diagEnergy_a
    real(8) :: diagEnergy_ao1, diagEnergy_ao1o2
    real(8) :: shift
    type (ivector), allocatable :: occA(:), occB(:), virA(:), virB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer, allocatable :: CIlevel(:)
    real(8) :: tmpconfCoreConfB
    integer(8) :: pi, qi, ri, si, pj, qj, rj, sj
    integer(8) :: pq, rs, rs_pair(2)
    integer(8) :: oia, oja, via, vja, aaa
    integer(8) :: oi1, vi1, oi2, vi2, oj2, vj2
    integer(8) :: size_ii, size_ij
    integer :: factor1, factor2, factor2j
    integer :: nonzero

!$  timeA = omp_get_wtime()
    shift = 1E-8 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !! work only with non-zero conf
    nonzero = 0
    do a = 1, SCICoreSpaceSize  
      if ( confCore(1)%values(1,a) == -1_1 .or. abs(coefficientCore(a)) <= 1E-10 ) exit
      nonzero = nonzero + 1
    enddo

    !$omp parallel &
    !$omp& private ( occA, occB, virA, virB, orbA, orbB, CIlevel) &
    !$omp& private ( n, a, oia, via, pi, qi, ri, si, oi1, vi1, oi2, vi2, spi, spj, oj2, vj2, factor1, factor2, factor2j, &
    !$omp&           CIenergy, diagEnergy, diagEnergy_a, diagEnergy_ao1, diagEnergy_ao1o2 ) 

    !! allocating auxiliary arrays (omp) for working with conf and orbitals
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( CIlevel ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    n = omp_get_thread_num() + 1

    !! loop to find all CI configurtions coupled to core space
    !!$omp do schedule (runtime) !with OMP_SCHEDULE for testing
    !$omp do schedule (dynamic)
    do a = 1, nonzero  

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0_8
        via = 0_8

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = confCore(spi)%values(:,a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1_8 ) then
            oia = oia + 1_8
            occA(spi)%values(oia) = pi
          else if ( orbA(spi)%values(pi) == 0_8 ) then
            via = via + 1_8
            virA(spi)%values(via) = pi
          end if
        enddo

        !! copy to conf B, these are variable
        orbB(spi)%values = orbA(spi)%values 
        occB(spi)%values = occA(spi)%values 
        virB(spi)%values = virA(spi)%values 

      enddo
      CIlevel = 0

      !! building all single sustitutions from configuration A. 
      !! here all configurations pairs are generated in maximum coincidence 
      do spi = 1, numberOfSpecies 

        !! calculate the sign factor for canonical order of the configuration
        factor1 = CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )

        do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1_8

          do qi = 1_8, CIcore_instance%numberOfActiveOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi) !! occ or core???
            vi1 = virA(spi)%values(qi)
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) + 1_8
            occB(spi)%values(pi) = vi1

            CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

            !! bit mapping from orbital to decimal num
            !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )
            !confCoreConfB(spi)%values = orbB(spi)%values ! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! get spingle sustitutions energy, H_ab C_a
            CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, oi1, vi1  )

            diagEnergy = CISCI_calculateEnergyZero( occB )
            CIenergy = abs(CIenergy * coefficientCore(a) ) / ( diagEnergy - oldEnergy + shift)

            !! calculate the sign factor for canonical order of the configuration
            !factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
            !CIenergy = abs(CIenergy * factor1 * factor2)

            if ( CIenergy > CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(mode) ) then
              !! append the amplitude 
              call CISCI_appendAmplitude ( n, CIenergy, orbB )
            endif
            !! reset the confB
            occB(spi)%values(pi) = occA(spi)%values(pi)  
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) - 1_8
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1_8
        enddo !pi
        CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi) + 1_8:) )

        !! generate double excitations (same species)
        do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1_8

          do qi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            oi2 = occA(spi)%values(qi)  
            if ( oi1 <= oi2 ) cycle !?
            orbB(spi)%values(oi2) = orbB(spi)%values(oi2) - 1_8

            pq = CIcore_instance%twoIndexArray(spi)%values(oi2,oi1)

            !do rs = 1_8, CISCI_instance%heatBathDoubleExcitations_size(spi,spi)%values(pq,mode) 
            do rs = 1_8, size_ii

              CIenergy = CISCI_instance%heatBathDoubleExcitations(spi,spi)%values(rs,pq)
              CIenergy = abs( CIenergy * coefficientCore(a) )


              if ( CIenergy > CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(mode) ) then

                !! get the different orbitals (excitation)
                rs_pair = IndexMap_vectorToMatrix( CISCI_instance%heatBathDoubleExcitations_index(spi,spi)%values(rs,pq), CIcore_instance%numberOfActiveOrbitals%values(spi))
                vi1 = rs_pair(1)
                vi2 = rs_pair(2)

                if ( orbB(spi)%values(vi1) == 1_8 ) cycle ! already occupied
                if ( orbB(spi)%values(vi2) == 1_8 ) cycle ! already occupied

                !! generate the configuration
                orbB(spi)%values(vi1) = 1_8
                orbB(spi)%values(vi2) = 1_8

                CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )

                occB(spi)%values(pi) = vi1
                occB(spi)%values(qi) = vi2
                diagEnergy = CISCI_calculateEnergyZero( occB )
                CIenergy = abs(CIenergy * coefficientCore(a) ) / ( diagEnergy - oldEnergy + shift)
                occB(spi)%values(pi) = occA(spi)%values(pi)  
                occB(spi)%values(qi) = occA(spi)%values(qi)  

                !! append the amplitude 
                call CISCI_appendAmplitude ( n, CIenergy, orbB )

                !! restore the configuration
                orbB(spi)%values(vi1) = 0_8
                orbB(spi)%values(vi2) = 0_8
              endif !!
            enddo !! rs
            !! reset the confB
            orbB(spi)%values(oi2) = orbB(spi)%values(oi2) + 1_8
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1_8
        enddo !pi

        !! generate double excitations (diff species)
        do spj = spi + 1, numberOfSpecies

          size_ij = ( CIcore_instance%numberOfActiveOrbitals%values(spi) * CIcore_instance%numberOfActiveOrbitals%values(spj) )
          do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            oi1 = occA(spi)%values(pi)  
            orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1_8
            do qj = CIcore_instance%numberOfCoreOrbitals%values(spj) + 1_8, CIcore_instance%numberOfOccupiedOrbitals%values(spj)
              oj2 = occA(spj)%values(qj)  
              orbB(spj)%values(oj2) = orbB(spj)%values(oj2) - 1_8

              pq = qj + CIcore_instance%numberOfActiveOrbitals%values(spj) * ( pi - 1_8 )

              do rs = 1_8, size_ij
              !do rs = 1_8, CISCI_instance%heatBathDoubleExcitations_size(spj,spi)%values(pq,mode) 
                CIenergy = CISCI_instance%heatBathDoubleExcitations(spj,spi)%values(rs,pq)
                CIenergy = abs(CIenergy * coefficientCore(a) )

                if ( CIenergy > CONTROL_instance%CISCI_HEAT_BATH_THRESHOLD(mode) ) then

                  !! get the different orbitals (excitation)
                  vi1 = ((CISCI_instance%heatBathDoubleExcitations_index(spj,spi)%values(rs,pq) - 1_8 ) / CIcore_instance%numberOfActiveOrbitals%values(spj)) + 1_8
                  vj2 = mod(CISCI_instance%heatBathDoubleExcitations_index(spj,spi)%values(rs,pq) - 1_8, CIcore_instance%numberOfActiveOrbitals%values(spj)) + 1_8

                  if ( orbB(spi)%values(vi1) == 1_8 ) cycle ! already occupied
                  if ( orbB(spj)%values(vj2) == 1_8 ) cycle ! already occupied

                  !! generate the configuration
                  orbB(spi)%values(vi1) = 1_8
                  orbB(spj)%values(vj2) = 1_8

                  CIlevel(spi) = sum(orbB(spi)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spi)+1:) )
                  CIlevel(spj) = sum(orbB(spj)%values(CIcore_instance%numberOfOccupiedOrbitals%values(spj)+1:) )

                  occB(spi)%values(pi) = vi1
                  occB(spj)%values(qj) = vj2
                  diagEnergy = CISCI_calculateEnergyZero( occB )
                  CIenergy = abs(CIenergy * coefficientCore(a) ) / ( diagEnergy - oldEnergy + shift)
                  occB(spi)%values(pi) = occA(spi)%values(pi)  
                  occB(spj)%values(qj) = occA(spj)%values(qj)  

                  !! append the amplitude 
                  call CISCI_appendAmplitude ( n, CIenergy, orbB )

                  !! restore the configuration
                  orbB(spi)%values(vi1) = 0_8
                  orbB(spj)%values(vj2) = 0_8
                endif !!
              enddo !! rs
              !! reset the confB
              orbB(spj)%values(oj2) = orbB(spj)%values(oj2) + 1_8
            enddo !qj
            orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1_8
          enddo !pi
        enddo !spj
      enddo !spi
    enddo !enddo a
    !$omp enddo 

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( virA(spi) )  
      call Vector_destructorInteger ( virB(spi) )  
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( CIlevel )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    !$omp end parallel

    !! ------------------------------
    !! the above code applies the denominator from eq 4 10.1063/1.4955109 on the fly for each configuration. Alternatively, it can be sort by heat-bath ranking first, then final sort by asci ranking

    !! sort and reduce the target arrays among all OMP threads, final run
    call CISCI_sortAmplitude( CIcore_instance%nproc + 1 ) 

!$  timeB = omp_get_wtime()
!$  write(*,"(A,ES10.2,A4)") "** TOTAL Elapsed Time for generating HB-SCI configurations: ", timeB - timeA ," (s)"

  end subroutine CISCI_heatBathGenerate 

end module CISCI_
