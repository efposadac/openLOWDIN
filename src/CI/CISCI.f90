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
    type (Vector8) :: buffer_amplitudeCore
    type (Vector8) :: coefficientCore
    type (Vector8) :: auxcoefficientTarget
    !! auxiliary array to trace the original elements in array sorting
    type (ivector8) :: index_amplitudeCore
    !! arrays for storing CI configurations, species, orbitals, vector size
    type (IMatrix1), allocatable :: confCore(:)
    integer(1), allocatable :: confAmplitudeCore(:,:)
    type (IMatrix1), allocatable :: saved_confTarget(:)
    !! storing the CI diagonal matrix elements for Jadamilu preconditioner
    type (Vector8) :: diagonalCore
    type (Vector8) :: diagonalTarget
    !! eigenvalues per SCI iteration
    type (Vector8) :: eigenValues
    !! length of SCI search vectors
    integer(8) :: coreSpaceSize
    integer(8) :: targetSpaceSize
    integer(8) :: buffer_amplitudeCoreSize
    !! PT2 perturbative energy correction to SCI
    real(8) :: PT2energy
    !! auxiliary variables to map orbitals from vector to array location
    integer :: combinedNumberOfOrbitals
    integer, allocatable :: combinedOrbitalsPositions(:,:) ! lower and upper position, species
    !! auxiliary array to store the position of the target space for each omp thread 
    integer(8), allocatable :: omp_targetInterval(:,:) ! lower and upper position, n_threads
    integer(8), allocatable :: omp_target_iterator_m(:) ! n_threads
    !! Stored code for bit-masking representation
    !!store the orbitals for each target configurations, to avoid recomputing them 
    !!type (IVector), allocatable :: targetOrb(:,:) ! species, num of target configurations % num. of orbitals

  end type CISCI

  type(CISCI) :: CISCI_instance

contains

  !! show information about the SCI calculation, creators, and libraries used
  subroutine CISCI_show()
    implicit none
    integer :: spi, spj
    integer(8) :: totalSize

    !! estimate memory cost 

    !! take from input
    CISCI_instance%coreSpaceSize = CONTROL_instance%CI_SCI_CORE_SPACE
    CISCI_instance%targetSpaceSize = CONTROL_instance%CI_SCI_TARGET_SPACE * CIcore_instance%nproc 
    !! the first quarter is the real target space, the second quarter is the waiting list if the amplitude could grow more. 
    !! The last half is just a temporary space to avoid sorting a big array for a single addition
    CISCI_instance%buffer_amplitudeCoreSize = CISCI_instance%targetSpaceSize * 8 

    totalSize = 0_8
    do spi = 1, CIcore_instance%numberOfSpecies 
      totalSize = totalSize + &
                  ( CISCI_instance%buffer_amplitudeCoreSize * ( 8 + 8 + 1*CIcore_instance%numberOfOrbitals%values(spi) )  + & ! data type for coeff, index, conf
                  CISCI_instance%coreSpaceSize * ( 8 + 8 + 1*CIcore_instance%numberOfOrbitals%values(spi)) + & ! coeff, diagonal conf
                  CISCI_instance%targetSpaceSize * ( 8 + 8 + 2*8 + 1*CIcore_instance%numberOfOrbitals%values(spi) ) ) * & !! coeff, diagonal, eigenvectors, conf
                  CIcore_instance%nproc ! per OMP thread
      do spj = spi, CIcore_instance%numberOfSpecies 
        totalSize = totalSize + &
                    size(CIcore_instance%fourCenterIntegrals( spi, spj )%values,1) * 8
      enddo
    enddo
    print *, totalSize

    write(6,*) "-----------------------------------------------------------------------"
    write (6,"(T2,A62)") "          SELECTED CONFIGURATION INTERACTION (SCI):          " 
    write (6,"(T2,A62)") "                 Adaptive Sampling CI (ASCI)                 " 
    write (6,"(T2,A62)") "                   Deterministic Algorithm                   " 
    write (6,"(T2,A62)") "                  Based on 10.1063/1.4955109                 "
    write (6,"(T2,A62)") "                        J. Charry                            "
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""
    write(6,*) "  Diagonalizer for target space hamiltonian : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
    write(6,*) "a software code for computing selected eigenvalues of "
    write(6,*) "large sparse symmetric matrices, "
    write(6,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""
    write(6,*) " Modified sorting algorithm from CENCALC quicksort code "
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) " Code available at https://github.com/dimassuarez/cencalc_quicksort  "
    write(6,*) " E. Suárez, N. Díaz, J. Méndez and D. Suárez. "
    write(6,*) " CENCALC: A Computational Tool for Conformational Entropy Calculations"
    write(6,*) " from Molecular Simulations."
    write(6,*) " J. Comput. Chem. 54, 2031. DOI: 10.1002/jcc.23350 "
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""
    write (6,"(T2,A,F14.3,A3 )") "Estimated memory needed : ", float( totalSize )/(1024**2) , " MB"
    write (6,"(T2,A,F14.3,A3 )") "                          ", float( totalSize )/(1024**3) , " GB"
    write (6,"(T2,A,I8 )") "Length of core (search) space                          :",  CISCI_instance%coreSpaceSize
    write (6,"(T2,A,I8 )") "Length of target (Full-CI subset) space per OMP thread :",  CISCI_instance%targetSpaceSize / CIcore_instance%nproc
    write (6,"(T2,A,I8 )") "Length of buffer (auxiliary sort) space per OMP thread :",  CISCI_instance%buffer_amplitudeCoreSize / CIcore_instance%nproc
    write (6,"(T2,A,I8 )") "Length of total target space                           :",  CISCI_instance%targetSpaceSize/ CIcore_instance%nproc
    write (6,"(T2,A,I8 )") "Length of total buffer space                           :",  CISCI_instance%buffer_amplitudeCoreSize / CIcore_instance%nproc
    write(6,*) "-----------------------------------------------------------------------"
    

  end subroutine CISCI_show

  !! Allocating arrays 
  subroutine CISCI_constructor( numberOfConfigurations )
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    integer(8), intent(out) :: numberOfConfigurations
    integer :: a,b,c,aa,bb,i, spi
    real(8) :: CIenergy
    integer :: nproc, n
    integer :: numberOfSpecies
    integer :: m

    numberOfSpecies = CIcore_instance%numberOfSpecies 
    numberOfConfigurations = CISCI_instance%targetSpaceSize

    !! auxiliary variables to map orbitals from vector to array location
    CISCI_instance%combinedNumberOfOrbitals = sum(CIcore_instance%numberOfOrbitals%values(:))
    allocate ( CISCI_instance%combinedOrbitalsPositions(2,numberOfSpecies) )
    m = 0
    do spi = 1, numberOfSpecies 
      CISCI_instance%combinedOrbitalsPositions(1,spi) = m + 1
      CISCI_instance%combinedOrbitalsPositions(2,spi) = m + CIcore_instance%numberOfOrbitals%values(spi)
      m = m + CIcore_instance%numberOfOrbitals%values(spi)
    enddo 

    !! auxiliary arrays to store the position of the target space for each omp thread within the big arrays 
    allocate ( CISCI_instance%omp_targetInterval(2, CIcore_instance%nproc + 1 ) ) !fixed
    allocate ( CISCI_instance%omp_target_iterator_m( CIcore_instance%nproc + 1) ) !variable

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

    !! arrays for storing coefficients
    call Vector_constructor8 ( CISCI_instance%buffer_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0.0_8)
    call Vector_constructor8 ( CISCI_instance%coefficientCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Vector_constructor8 ( CISCI_instance%auxcoefficientTarget, int(CISCI_instance%targetSpaceSize,8), 0.0_8) !! meh... just to avoid changing everything from matrix to vector format

    !! auxiliary array to trace the original elements in array sorting
    call Vector_constructorInteger8 ( CISCI_instance%index_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0_8)

    !! arrays for storing CI configurations, species, orbitals, vector size
    allocate ( CISCI_instance%confCore ( numberOfSpecies ) ) 
    !allocate ( CISCI_instance%confAmplitudeCore ( numberOfSpecies ) ) 
    allocate ( CISCI_instance%saved_confTarget ( numberOfSpecies ) ) 
    do spi = 1, numberOfSpecies 
      call Matrix_constructorInteger1 ( CISCI_instance%confCore(spi), int(CIcore_instance%numberOfOrbitals%values(spi),8) , int(CISCI_instance%coreSpaceSize,8), -1_1 )
      ! this was replaced by a "vectorized" array to avoid using arrays of types inside a recursive function
      !call Matrix_constructorInteger1 ( CISCI_instance%confAmplitudeCore(spi), int(CIcore_instance%numberOfOrbitals%values(spi),8) , int(CISCI_instance%buffer_amplitudeCoreSize,8), -1_1) 
      call Matrix_constructorInteger1 ( CISCI_instance%saved_confTarget(spi), int(CIcore_instance%numberOfOrbitals%values(spi),8) , int(CISCI_instance%targetSpaceSize,8), -1_1)
    enddo
    allocate ( CISCI_instance%confAmplitudeCore (  CISCI_instance%combinedNumberOfOrbitals, CISCI_instance%buffer_amplitudeCoreSize ) ) 
    CISCI_instance%confAmplitudeCore = -1_1

    !! storing the CI diagonal matrix elements for Jadamilu preconditioner
    call Vector_constructor8 ( CISCI_instance%diagonalCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Vector_constructor8 ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8) 

    !! eigenvalues per SCI iteration
    call Vector_constructor8 ( CISCI_instance%eigenValues, 20_8, 0.0_8) !! store the eigenvalues per macro iterations

    !! store the orbitals for each target configurations, to avoid recomputing them 
    !! this is helpful when using bit-masking approach to avoid transforming from bit to decimal multiple times
    !!allocate ( CISCI_instance%targetOrb ( numberOfSpecies, CISCI_instance%targetSpaceSize ) )
    !!do a = 1, CISCI_instance%targetSpaceSize
    !!  do spi = 1, numberOfSpecies 
    !!    call Vector_constructorInteger ( CISCI_instance%targetOrb(spi,a), CIcore_instance%numberOfOrbitals%values(spi), 0 ) 
    !!  enddo
    !!enddo

    a = 1
    !! initializing index arrays, this index is relative for each omp thread to simplify internal usage during sorting
    do n = 1, CIcore_instance%nproc 
      do m = 1,  CISCI_instance%omp_targetInterval(2, n ) - CISCI_instance%omp_targetInterval(1, n ) + 1 
         CISCI_instance%index_amplitudeCore%values(a) = m
         a = a + 1
      enddo
    enddo

    !! initialize sorting subroutines
    call CISort_constructor()

  end subroutine CISCI_constructor
 
  !! main part
  subroutine CISCI_run( eigenVectors )
    use sort_
    implicit none
    type(matrix) :: eigenVectors
    integer(8) :: a, aa, i, j, ii, jj, m
    integer :: k ! macro SCI iteration
    integer :: nproc, n
    real(8) :: timeA(20), timeB(20)
    real(8) :: timeAA, timeBB
    real(8) :: timeAS, timeBS
    type(Vector8) :: eigenValuesTarget
    real(8) :: minValue
    real(8) :: currentEnergy 
    integer :: numberOfSpecies, spi

    numberOfSpecies = CIcore_instance%numberOfSpecies 
    nproc = CIcore_instance%nproc 
    currentEnergy = HartreeFock_instance%totalEnergy 
    call Vector_constructor8 ( eigenValuesTarget, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    !! HF determinant coefficient
    CISCI_instance%coefficientCore%values(1) = 1.0_8
    CISCI_instance%eigenValues%values(1) = HartreeFock_instance%totalEnergy 

    write (6,*)    ""
    write (6,"(T2,A29 )")    "Starting SCI macro iterations "
    write (6,*)    ""

    do k = 2, 20
        
!$  timeA(k) = omp_get_wtime()

      !! calculating the amplitudes in core space. This is the pertubation guess of CI eigenvector
      CISCI_instance%buffer_amplitudeCore%values = 0.0_8
      if ( k == 2 ) call CISCI_initialConfigurations(  CISCI_instance%coefficientCore, CISCI_instance%confCore )

      call CISCI_core_amplitudes ( CISCI_instance%coefficientCore%values, CISCI_instance%confCore, CISCI_instance%coreSpaceSize, currentEnergy )

      if ( k == 2 ) CISCI_instance%buffer_amplitudeCore%values(1) = 1.0_8 !! at first iter, this coefficient diverges

      !! reduction 
      call CISCI_sortAmplitude( CIcore_instance%nproc + 1 ) !! final sort after applying the denominator from Epstein-Nesbet estimates, and after fixing HF coeff divergence

      !! saving top-targetSize index configurations sorted from confAmplitudeCore in saved_confTarget.
      do spi = 1, numberOfSpecies
        do i = 1, CISCI_instance%targetSpaceSize
          ii =  CISCI_instance%index_amplitudeCore%values(i) 
          CISCI_instance%saved_confTarget(spi)%values(:,i) = CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), ii) 
        enddo 
      enddo 

      !! computing the diagonal in the target space, jadamilu requires the diagonal in advance
      call CISCI_buildDiagonal ( )

      !! using the amplituded as the initial coeff guess, after that, use the previous diganolized eigenvectors in target space
      do a = 1,  CISCI_instance%targetSpaceSize
        aa = CISCI_instance%index_amplitudeCore%values(a)
        eigenVectors%values(a,1) = CISCI_instance%buffer_amplitudeCore%values(a)
      enddo
      
      !! eigenvalue guess
      eigenValuesTarget%values(1) = currentEnergy 

      !! diagonalize in target space
      call CISCI_jadamiluInterface( int(CISCI_instance%targetSpaceSize,8), &
                 1_8, &
                 eigenValuesTarget, &
                 eigenVectors, timeAA, timeBB )
  
      !! saving the eigenvectors coeff to an aux vector. Only ground state 
      CISCI_instance%auxcoefficientTarget%values(:) = eigenVectors%values(:,1)

      !! convergence criteria
      CISCI_instance%eigenValues%values(k) = eigenValuesTarget%values(1)
      if ( abs( CISCI_instance%eigenValues%values(k) - currentEnergy ) < 1.0E-5 ) then
        write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k , " Energy: ", CISCI_instance%eigenValues%values(k)
        !$  timeB(k) = omp_get_wtime()
        exit
      end if

      !! reset auxindex arrary, for later use in sorting target coeff. global absolute index
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo

      !! getting the core absolute largest coefficients
      call MTSort (  CISCI_instance%auxcoefficientTarget%values, &
           CISCI_instance%index_amplitudeCore%values,  int( CISCI_instance%coreSpaceSize ,8), "D", nproc )

      !! storing only the largest coefficients, and rearraing the next eigenvector guess 
      do i = 1,  CISCI_instance%coreSpaceSize
        CISCI_instance%coefficientCore%values(i) = CISCI_instance%auxcoefficientTarget%values(i)
      enddo

      !! storing sorted confAmplitudeCore to index_core
      do i = 1,  CISCI_instance%coreSpaceSize
        ii = CISCI_instance%index_amplitudeCore%values(i)
        do spi = 1, numberOfSpecies
          CISCI_instance%confCore(spi)%values(:,i) = CISCI_instance%saved_confTarget(spi)%values(:,ii) 
        enddo
      enddo

      !! reset iterators
      do n = 1, CIcore_instance%nproc 
        CISCI_instance%omp_target_iterator_m(n) = CISCI_instance%omp_targetInterval(1, n ) - 1
      enddo

      !! reset auxindex array. relative indexes
      i = 1
      do n = 1, CIcore_instance%nproc 
        do m = 1, CISCI_instance%omp_targetInterval(2, n ) - CISCI_instance%omp_targetInterval(1, n ) + 1  
           CISCI_instance%index_amplitudeCore%values(i) = m
           i = i + 1
        enddo
      enddo

      !! reset auxindex arrary, for later use in sorting target coeff
      !do i = 1, CISCI_instance%buffer_amplitudeCoreSize
      !  print *, "----",  i, CISCI_instance%index_amplitudeCore%values( i ) 
      !enddo

      !! restart amplitudes for next run, except when exiting to do PT2 corr
      CISCI_instance%buffer_amplitudeCore%values = 0.0_8
      do spi = 1, numberOfSpecies
        CISCI_instance%confAmplitudeCore = -1_1
        !CISCI_instance%confAmplitudeCore(spi)%values = -1_1
        CISCI_instance%saved_confTarget(spi)%values = -1_1
      enddo 

      write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k , " Energy: ", CISCI_instance%eigenValues%values(k)

!$  timeB(k) = omp_get_wtime()
      !! updating new reference
      currentEnergy = eigenValuesTarget%values(1)

    enddo !k

    !! summary of the macro iteration 
    write (6,*)    ""
    write (6,"(T2,A95 )")    "                            Selected CI (SCI)  summary                                               "
    write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "
    do k = 2, 20
       write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") k-1,  CISCI_instance%eigenValues%values(k),  &
                                                          CISCI_instance%eigenValues%values(k) - HartreeFock_instance%totalEnergy, &
                                                          CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1), &
                                                          timeB(k) - timeA(k)
      if ( abs( CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1)  ) < 1.0E-5 ) then 
        CIcore_instance%eigenvalues%values(1) = CISCI_instance%eigenValues%values(k) 
        exit
      endif
    enddo !k
    write (6,"(T2,A30)") "SCI Energy Convergence : 1E-5 " 
    write (6,*)    ""

    !! calculating PT2 correction. A pertuberd estimation of configurations not include in the target space
    call CISCI_PT2 ( CISCI_instance%targetSpaceSize, CISCI_instance%eigenValues%values(k), CISCI_instance%PT2energy )

    write (6,"(T2,A,F25.12)") "CI-PT2 energy correction :", CISCI_instance%PT2energy 

  end subroutine CISCI_run

  subroutine CISCI_initialConfigurations ( coefficientCore, confCore )

    implicit none
    type(vector8) :: coefficientCore
    type(IMatrix1) :: confCore(:)
    type(IVector), allocatable :: orbA(:)
    integer :: spi, numberOfSpecies
    integer :: pi
    integer(8) :: m 
    real(8) :: indexConf

    !! finding position in the global arrays
    CISCI_instance%omp_target_iterator_m(1) = CISCI_instance%omp_target_iterator_m(1) + 1 !! appending in the first thread id 
    m = CISCI_instance%omp_target_iterator_m(1) !! appending in the first thread id 

    !! Hartree-Fock reference coeff
    coefficientCore%values(m) = 1.0_8

    !! build orbitals references
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies 
     call Vector_constructorInteger( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0) 
    
      do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
        orbA(spi)%values(pi) = 1.0
      enddo

      confCore(spi)%values(:,m) = orbA(spi)%values(:)
      !call CISCI_binaryToDecimal ( orbA(spi)%values, indexConf )
      !confCore(spi)%values(spi,m) = indexConf
      
    enddo

    deallocate ( orbA  )

    !! append the initial configuration to the amplitude, otherwise this is zero at the first iteration 
    CISCI_instance%buffer_amplitudeCore%values(m) = 1.0_8
    do spi = 1, numberOfSpecies 
      CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), m) = confCore(spi)%values(:,m)
    enddo
    CISCI_instance%index_amplitudeCore%values(m) = 1
    

  end subroutine CISCI_initialConfigurations

  subroutine CISCI_core_amplitudes ( coefficientCore, confCore, SCICoreSpaceSize, oldEnergy )

    implicit none
  
    integer(8) SCICoreSpaceSize
    real(8) coefficientCore ( SCICoreSpaceSize )
    type(IMatrix1) :: confCore(:)
    real(8) :: CIEnergy
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: wi
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
    real(8) :: tmpconfCoreConfB
    integer :: pi, qi, ri, si, pj, qj, rj, sj
    integer :: oia, oja, via, vja, aaa
    integer :: oi1, vi1, oi2, vi2, oj2, vj2
    integer :: factor1, factor2, factor2j
    integer :: nonzero

!$  timeA = omp_get_wtime()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
    shift = 1E-6 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !! work only with non-zero conf
    nonzero = 0
    do a = 1, SCICoreSpaceSize  
      if ( confCore(1)%values(1,a) == -1_1 ) exit
      nonzero = nonzero + 1
    enddo

    !$omp parallel &
    !$omp& private ( occA, occB, virA, virB, orbA, orbB) &
    !$omp& private ( n, a, oia, via, pi, qi, ri, si, oi1, vi1, oi2, vi2, spi, spj, oj2, vj2, factor1, factor2, factor2j, CIenergy ) 
    !$omp& shared ( nonzero, CISCI_instance, coefficientCore, confCore)

    !! allocating auxiliary arrays (omp) for working with conf and orbitals
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    n = omp_get_thread_num() + 1

    !! loop to find all CI configurtions coupled to core space
    !$omp do schedule (static) 
    do a = 1, nonzero !! coreSpace
      
      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 
        via = 0

        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = confCore(spi)%values(:,a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
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


      !! building all single sustitutions from configuration A. 
      !! here all configurations pairs are generated in maximum coincidence 
      do spi = 1, numberOfSpecies 

        !! calculate the sign factor for canonical order of the configuration
        factor1 = CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )

        do pi = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1 

          do qi = 1, CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
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

            if ( CIenergy /= 0.0_8 ) then
              CIenergy = CIenergy * coefficientCore(a) 

              ! alternative approximation
              !diagEnergy = CISCI_calculateEnergyZero( occA )
              !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

              !! append the amplitude 
              call CISCI_appendAmplitude ( n, CIenergy, orbB )
            endif

            !tmpconfCoreConfB = confCoreConfB(spi) !! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! building all double intraspecies sustitutions from configuration A
            do ri = CIcore_instance%numberOfCoreOrbitals%values(spi) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
              oi2 = occA(spi)%values(ri)  
              if ( oi1 <= oi2 ) cycle 
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) - 1 
              do si = 1, CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
                vi2 = virA(spi)%values(si)
                if ( vi1 <= vi2 ) cycle 
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) + 1
                occB(spi)%values(ri) = vi2 

                !! get double intraspecies sustitutions energy
                CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, oi1, oi2, vi1,  vi2 )

                if ( CIenergy /= 0.0_8 ) then
                  !! calculate the sign factor for canonical order of the configuration
                  factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
 
                  CIenergy = CIenergy * factor1 * factor2

                  !! bit mapping from orbital to decimal num
                  !call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )

                  CIenergy = CIenergy * coefficientCore(a)
                  ! alternative approximation
                  !diagEnergy = CISCI_calculateEnergyZero( occA )
                  !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                  !! append the amplitude 
                  call CISCI_appendAmplitude ( n, CIenergy, orbB )

                endif

                occB(spi)%values(ri) = occA(spi)%values(ri)  
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) -1  ! reset orbital 
              enddo
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) + 1 ! reset orbital
            enddo

            !! building all double interspecies sustitutions from configuration A. maybe build this as a superloop?
            !! get double interspecies sustitutions energy
            do spj = spi + 1, numberOfSpecies 
              !if ( spj == spi ) cycle 
              do rj = CIcore_instance%numberOfCoreOrbitals%values(spj) + 1, CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                oj2 = occA(spj)%values(rj)  
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) - 1 
                do sj = 1, CIcore_instance%numberOfOrbitals%values(spj) - CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                  vj2 = virA(spj)%values(sj)
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) + 1
                  occB(spj)%values(rj) = vj2 

                  !! get double interspecies sustitutions energy
                  CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, oi1, oj2, vi1, vj2 )

                  if ( CIenergy /= 0.0_8 ) then

                    !! calculate the sign factor for canonical order of the configuration
                    factor2j = CISCI_canonicalOrderFactor( spj, orbB(spj), occB(spj) )
   
                    CIenergy = CIenergy * factor2 * factor2j

                    !! bit mapping from orbital to decimal num
                    !call CISCI_binaryToDecimal ( orbB(spj)%values, confCoreConfB(spj) )

                    CIenergy = CIenergy * coefficientCore(a) 
                    ! alternative approximation
                    !diagEnergy = CISCI_calculateEnergyZero( occA )
                    !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                    !! append the amplitude 
                    call CISCI_appendAmplitude ( n, CIenergy, orbB )
                  endif

                  !! reset the confB
                  occB(spj)%values(rj) = occA(spj)%values(rj)  
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) -1 
                enddo ! sj
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) + 1
              enddo ! rj

            enddo !spj

            !! reset the confB
            occB(spi)%values(pi) = occA(spi)%values(pi)  
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) -1 
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1
        enddo !pi

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

    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    !$omp end parallel

    !! ------------------------------
    !! apply the denominator from eq 4 10.1063/1.4955109. First sort by heat-bath ranking, then final sort by asci ranking

    !! allocate these again, but now in serial
    allocate ( occA ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    !! reduce the target arrays, final run
    call CISCI_sortAmplitude( CIcore_instance%nproc + 1 ) 

    do aa = 1, CISCI_instance%buffer_amplitudeCoreSize 
      a = CISCI_instance%index_amplitudeCore%values(aa) 

      if (  CISCI_instance%confAmplitudeCore(1, a ) == -1_1 ) exit 

      do spi = 1, numberOfSpecies 
        occA(spi)%values(:) = 0
        virA(spi)%values(:) = 0
        orbA(spi)%values(:) = 0
      enddo

      do spi = 1, numberOfSpecies 

        oia = 0 
        via = 0

        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi, a ),  orbA(spi)%values )
        orbA(spi)%values(:) = CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          else if ( orbA(spi)%values(pi) == 0 ) then
            via = via + 1
            virA(spi)%values(via) = pi
          end if
        enddo

      enddo

      !! calculating the diagonal elements and the denominator
      CIenergy = CISCI_calculateEnergyZero( occA )
      CISCI_instance%buffer_amplitudeCore%values( aa ) = CISCI_instance%buffer_amplitudeCore%values( aa ) / ( CIenergy - oldEnergy + shift)

    enddo    

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( virA(spi) )  
      call Vector_destructorInteger ( orbA(spi) ) 
    end do

    deallocate ( occA  )
    deallocate ( virA  )
    deallocate ( orbA  )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for calculating SCI amplitudes : ", timeB - timeA ," (s)"

  end subroutine CISCI_core_amplitudes

  subroutine CISCI_jadamiluInterface(n,  maxeig, eigenValues, eigenVectors, timeA, timeB)
    implicit none
    external DPJDREVCOM
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
    integer(8) :: JA(1), IA(1)
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
    real(8) :: timeA, timeB
    
!$  timeA = omp_get_wtime()
    maxsp = CONTROL_instance%CI_MADSPACE

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
     IPRINT = 0 !     standard report on standard output
     ISEARCH = 1 !    we want the smallest eigenvalues
     NEIG = maxeig !    number of wanted eigenvalues
     !NINIT = 0 !    no initial approximate eigenvectors
     NINIT = NEIG !    initial approximate eigenvectors
     MADSPACE = maxsp !    desired size of the search space
     ITER = 1000*NEIG !    maximum number of iteration steps
     TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual
     TOL = 1e-4 !1.0d-4 !    tolerance for the eigenvector residual, for ASCI this can be higher

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

     JA(1) = -1 
     IA(1) = -1 

     ! set initial eigenpairs
     do j = 1, n 
       X(j) = eigenVectors%values(j,1)
     end do

     do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
       EIGS(i) = eigenValues%values(i)
     end do

     DROPTOL = 1E-4

     SIGMA = EIGS(1)
     gap = 0 
     SHIFT = 0!EIGS(1)

     do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
       write(6,"(T2,A5,I4,2X,A10,F20.10,2X,A11,F20.10)") "State", i, "Eigenvalue", eigs( i ), "Eigenvector", x((i-1)*n + i)
     end do

     iiter = 0
  
10   CALL DPJDREVCOM( N, CISCI_instance%diagonalTarget%values , JA, IA, EIGS, RES, X, LX, NEIG, &
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
        call CISCI_matvec ( N, X(NDX1), X(NDX2), iiter)
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

!$  timeB = omp_get_wtime()

  end subroutine CISCI_jadamiluInterface

  subroutine CISCI_matvec ( nx, v, w, iter)
  
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
    integer(4) :: iter
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
    integer :: diffOrbi(4)
    integer :: diffOrbj(4)
    integer :: pi
    integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA, factorB

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()

    nonzero = 0
    nonzerow = 0
    w = 0.0_8 
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 
    do a = 1 , nx
       if ( abs(v(a) ) >= tol) nonzero = nonzero + 1
    end do

    !! Stored code for bit-masking representation
    !! transforming from decimal to binary once
    !!do aa = 1, nonzero
    !!  do spi = 1, numberOfSpecies
    !!    !call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,aa), CISCI_instance%targetOrb(spi,aa)%values )
    !!    CISCI_instance%targetOrb(spi,aa)%values = CISCI_instance%saved_confTarget(spi)%values(:,aa)
    !!  enddo
    !!enddo

!$    timeA= omp_get_wtime()

    !$omp parallel &
    !$omp& private(aa, a, spi, oia, orbA, pi, occA, CIenergy, bb, b, oib, orbB, occB, couplingS, coupling, i, ii, diffOrbi, diffOrbj, spj, factorA, factorB ) 
    !$omp& shared ( w, v, numberOfSpecies, CIcore_instance)
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do
    !$omp do schedule (static) 
    aloop: do aa = 1, nonzero

      !a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      a = aa ! if index_amplitude is sorted

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        !! build the orbital from the index using the bit-masking
        !!call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
        !!orbA(spi)%values = CISCI_instance%targetOrb(spi,a)%values
        orbA(spi)%values(:) = CISCI_instance%saved_confTarget(spi)%values(:,a)

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo

      enddo

      bloop: do bb = aa, nonzero

        !b = CISCI_instance%index_amplitudeCore%values(bb)
        b = bb
        ! getting configuration B
        do spi = 1, numberOfSpecies 

          !! build the orbital from the index using the bit mapping
          !! call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,b), orbB(spi)%values )
          !! orbB(spi)%values = CISCI_instance%targetOrb(spi,b)%values
          orbB(spi)%values(:) = CISCI_instance%saved_confTarget(spi)%values(:,b)

        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) 
        end do
      
        !! if any species differs in more than 2 orb, skip
        if ( sum(couplingS) <= 2 ) then
          do spi = 1, numberOfSpecies 
            oib = 0 
            !! build auxiliary vectors of occupied and virtuals orbitals
            do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
              if ( orbB(spi)%values(pi) == 1 ) then
                oib = oib + 1
                occB(spi)%values(oib) = pi
              end if
            enddo
          enddo
        else 
          cycle 
        endif

        !! same conf 
        if ( sum(couplingS) == 0 ) then
          CIenergy = CISCI_instance%diagonalTarget%values(aa) 
          !$omp atomic
          w(aa) = w(aa) + CIenergy * v(aa)
          !$omp end atomic
        endif 
        !! one orbital different
        if ( sum(couplingS) == 1 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CISCI_calculateEnergyOne( spi, occA, occB, diffOrbi(1), diffOrbi(3)  )

          !$omp atomic
          w(bb) = w(bb) + CIenergy * v(aa) * factorA 
          !$omp end atomic
          !$omp atomic
          w(aa) = w(aa) + CIenergy * v(bb) * factorA
          !$omp end atomic
        endif
        !! two orbital different, same species
        if ( sum(couplingS) == 2 .and. maxval(couplingS) == 2 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CISCI_calculateEnergyTwoSame( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  )

          !$omp atomic
          w(bb) = w(bb) + CIenergy * v(aa) * factorA 
          !$omp end atomic
          !$omp atomic
          w(aa) = w(aa) + CIenergy * v(bb) * factorA 
          !$omp end atomic
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
          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj), orbB(spj), occA(spj), occB(spj), factorB )
          CIenergy = CISCI_calculateEnergyTwoDiff( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  )

          !$omp atomic
          w(bb) = w(bb) + CIenergy * v(aa) * factorA * factorB
          !$omp end atomic
          !$omp atomic
          w(aa) = w(aa) + CIenergy * v(bb) * factorA * factorB
          !$omp end atomic
        endif
 
      end do bloop !b
    end do aloop !a 
    !$omp end do nowait

    deallocate ( couplingS )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( orbA  )
    deallocate ( orbB  )
    !$omp end parallel

!$  timeB = omp_get_wtime()

    !! to check how dense is the w vector
    do a = 1 , nx
       if ( abs(w(a) ) >= tol) nonzerow = nonzerow + 1
    end do

!$    write(*,"(A,I2,A,E10.3,A2,I12,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero, nonzerow
    return

  end subroutine CISCI_matvec

  subroutine CISCI_buildDiagonal ( )
    implicit none
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
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()

    allocate ( occA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) 
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

!$    timeA= omp_get_wtime()

    do aa = 1, CISCI_instance%targetSpaceSize 

      !a = CISCI_instance%index_amplitudeCore%values(aa)
      a = aa
      if ( CISCI_instance%saved_confTarget(1)%values(1,a) == -1_1 ) exit
      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 
        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = CISCI_instance%saved_confTarget(spi)%values(:,a)
        !!call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi,a), orbA(spi)%values )

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo

      enddo

      CISCI_instance%diagonalTarget%values(aa) = CISCI_calculateEnergyZero( occA )

    end do !a 

!$  timeB = omp_get_wtime()

    deallocate ( occA  )
    deallocate ( orbA  )

    return

  end subroutine CISCI_buildDiagonal

  function CISCI_calculateEnergyOne( si, occA, occB, a, b ) result ( CIenergy )
    implicit none
    type(ivector), intent(in) :: occA(:), occB(:)
    integer, intent(in) :: a, b
    integer, intent(in) :: si
    integer :: sj
    integer(8) :: ab, ll, abll, albl, auxab
    integer :: l, la
    real(8) :: CIenergy

    CIenergy = 0.0_8
  
    CIenergy = CIenergy + CIcore_instance%twoCenterIntegrals(si)%values( a, b )

    ab = CIcore_instance%twoIndexArray(si)%values( a, b )
  
    do la = 1, CIcore_instance%occupationNumber( si ) !! the same orbitals pair are excluded by the exchange
  
        l = occA(si)%values(la) ! or b, both are the same
  
        ll = CIcore_instance%twoIndexArray(si)%values( l,l ) 
        abll = CIcore_instance%fourIndexArray(si)%values( ab, ll )
  
        CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals(si,si)%values( abll, 1)
  
        albl = CIcore_instance%fourIndexArray(si)%values( &
                               CIcore_instance%twoIndexArray(si)%values( a,l ), &
                               CIcore_instance%twoIndexArray(si)%values( l,b ) ) 
  
        CIenergy = CIenergy + MolecularSystem_instance%species(si)%kappa*CIcore_instance%fourCenterIntegrals(si,si)%values(albl, 1)
  
      end do
  
      do sj = 1, si - 1 !! avoid ii, same species

        auxab = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( ab - 1)

        do la =1,  CIcore_instance%occupationNumber( sj ) 

          l = occA(sj)%values(la) ! or b, both are the same

          abll = auxab  + CIcore_instance%twoIndexArray(sj)%values( l, l) 

          CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals( si, sj )%values(abll, 1) 

        end do

      end do

      do sj = si + 1, MolecularSystem_instance%numberOfQuantumSpecies !! avoid ii, same species

        auxab = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( ab - 1)

        do la = 1,  CIcore_instance%occupationNumber( sj )

          l = occA(sj)%values(la) ! or b, both are the same

          abll = auxab  + CIcore_instance%twoIndexArray(sj)%values( l, l) 

          CIenergy = CIenergy + CIcore_instance%fourCenterIntegrals( si, sj )%values(abll, 1) 

        end do

      end do

  end function CISCI_calculateEnergyOne

  function CISCI_calculateEnergyTwoSame( si, occA, occB, ai, aj, bi, bj ) result ( CIenergy )
    implicit none
    type(ivector), intent(in) :: occA(:), occB(:)
    integer, intent(in) :: ai, aj, bi, bj
    integer, intent(in) :: si
    integer :: sj
    integer(8) :: aibi_ajbj, aibj_ajbi
    integer :: l, la
    real(8) :: CIenergy

    CIenergy = 0.0_8

    aibi_ajbj = CIcore_instance%fourIndexArray(si)%values( &
                    CIcore_instance%twoIndexArray(si)%values( ai, bi ), &
                    CIcore_instance%twoIndexArray(si)%values( aj, bj ) )
  
    CIenergy = CIcore_instance%fourCenterIntegrals(si,si)%values( aibi_ajbj, 1)

    aibj_ajbi = CIcore_instance%fourIndexArray(si)%values( &
                    CIcore_instance%twoIndexArray(si)%values( ai, bj ), &
                    CIcore_instance%twoIndexArray(si)%values( aj, bi ) )
 
    CIenergy = CIenergy + MolecularSystem_instance%species(si)%kappa * & 
                          CIcore_instance%fourCenterIntegrals(si,si)%values( aibj_ajbi, 1)
  
  end function CISCI_calculateEnergyTwoSame

  function CISCI_calculateEnergyTwoDiff( si, sj, ai, aj, bi, bj ) result ( CIenergy )
    implicit none
    integer, intent(in) :: ai, aj, bi, bj
    integer, intent(in) :: si, sj
    integer(8) :: aibi,  aux_aibi, ajbj
    real(8) :: CIenergy

    CIenergy = 0.0_8

    aibi = CIcore_instance%twoIndexArray(si)%values( ai, bi )
    aux_aibi = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( aibi - 1 ) 

    ajbj = CIcore_instance%twoIndexArray(sj)%values( aj, bj )
    CIenergy = CIcore_instance%fourCenterIntegrals( si, sj )%values( aux_aibi +  ajbj, 1)

  end function CISCI_calculateEnergyTwoDiff

  function CISCI_calculateEnergyZero( occA ) result (CIenergy)
    implicit none

    type(ivector), intent(in) :: occA(:)
    integer :: si,sj
    integer :: ki,li,lj,k,l,kk,ll,kkll,kllk
    integer(8) :: auxIndex1, auxIndex2, auxIndex
    real(8) :: CIenergy

    CIenergy = 0.0_8

    do si = 1, MolecularSystem_instance%numberOfQuantumSpecies
      do ki =1, CIcore_instance%occupationNumber( si )  !! 1 is from a and 2 from b

        k = occA( si )%values( ki )

        !One particle terms
        CIenergy = CIenergy + CIcore_instance%twoCenterIntegrals( si )%values( k, k )

        !Two particles, same specie
        kk = CIcore_instance%twoIndexArray( si )%values( k, k)

        do li = ki + 1, CIcore_instance%occupationNumber( si )  !! 1 is from a and 2 from b

          l = occA( si )%values( li )
          ll = CIcore_instance%twoIndexArray( si )%values( l,l )
          kkll = CIcore_instance%fourIndexArray( si )%values( kk, ll ) 

          !Coulomb
          CIenergy = CIenergy + &
              CIcore_instance%fourCenterIntegrals( si, si )%values( kkll, 1)

          !Exchange, depends on spin

          kllk = CIcore_instance%fourIndexArray( si )%values( &
                        CIcore_instance%twoIndexArray( si )%values(k,l), &
                        CIcore_instance%twoIndexArray( si )%values(l,k) )

          CIenergy = CIenergy + &
                  MolecularSystem_instance%species( si )%kappa*CIcore_instance%fourCenterIntegrals( si, si )%values( kllk, 1)
        end do

        !!Two particles, different species
        do sj = si + 1, MolecularSystem_instance%numberOfQuantumSpecies

          do lj = 1, CIcore_instance%occupationNumber( sj ) !! 1 is from a and 2 from b
            l = occA( sj )%values(lj)

            ll = CIcore_instance%twoIndexArray( sj )%values(l,l)
            kkll = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * (kk - 1 ) + ll

            CIenergy = CIenergy + &
            CIcore_instance%fourCenterIntegrals( si, sj )%values( kkll, 1)

          end do

        end do

      end do
    end do

    CIenergy = CIenergy + HartreeFock_instance%puntualInteractionEnergy

  end function CISCI_calculateEnergyZero

  subroutine CISCI_PT2 ( SCITargetSpaceSize, refEnergy, energyCorrection )
  
    implicit none
    integer(8) :: SCITargetSpaceSize
    real(8) :: refEnergy
    real(8) :: energyCorrection
    real(8) :: CIEnergy
    integer(8) :: nonzero
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: wi
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
    integer(1) :: coupling
    integer :: diffOrbi(4)
    integer :: diffOrbj(4)

!$  timeA = omp_get_wtime()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !$omp parallel &
    !$omp& private(aa, a, spi, oia, orbA, pi, occA, CIenergy, bb, b, oib, orbB, occB, couplings, coupling, i, ii, diagonal, denominator, diffOrbi, diffOrbj, spj, factorA, factorB )
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) 
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    energyCorrection = 0.0_8
    
    !$omp& reduction (+:energyCorrection)
    !$omp do schedule (static) 
    aloop: do aa = CISCI_instance%targetSpaceSize + 1, CISCI_instance%buffer_amplitudeCoreSize

      a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      !a = aa ! if index_amplitude is sorted
      if (CISCI_instance%confAmplitudeCore(1,a) == -1_1) cycle ! cycle or exit?
      !if (CISCI_instance%confAmplitudeCore(1)%values(1,a) == -1_1) cycle ! cycle or exit?

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
        orbA(spi)%values(:) = CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), a) 

        !! build auxiliary vectors of occupied and virtuals orbitals
        do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
          if ( orbA(spi)%values(pi) == 1 ) then
            oia = oia + 1
            occA(spi)%values(oia) = pi
          end if
        enddo

      enddo

      CIenergy = 0.0_8
      bloop: do bb = 1,  CISCI_instance%targetSpaceSize

        !b = CISCI_instance%index_amplitudeCore%values(bb)
        b = bb
        !if (CISCI_instance%saved_confTarget(1)%values(1,b) == -1_1) exit ! this is never ocurring

        ! getting configuration B
        do spi = 1, numberOfSpecies 

          oib = 0 

          !! build the orbital from the index using the bit mapping
          !call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,b), orbB(spi)%values )
          orbB(spi)%values(:) = CISCI_instance%saved_confTarget(spi)%values(:,b)

          !! build auxiliary vectors of occupied and virtuals orbitals
          do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
            if ( orbB(spi)%values(pi) == 1 ) then
              oib = oib + 1
              occB(spi)%values(oib) = pi
            end if
          enddo
        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) 
        end do

        !! one orbital different
        if ( sum(couplingS) == 1 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 1 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CIenergy + CISCI_calculateEnergyOne( spi, occA, occB, diffOrbi(1), diffOrbi(3)  ) * factorA * CISCI_instance%auxcoefficientTarget%values(bb)

        endif
        !! two orbital different, same species
        if ( sum(couplingS) == 2 .and. maxval(couplingS) == 2 ) then
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CIenergy + CISCI_calculateEnergyTwoSame( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  ) * factorA * CISCI_instance%auxcoefficientTarget%values(bb)

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

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj), orbB(spj), occA(spj), occB(spj), factorB )
          CIenergy = CIenergy + CISCI_calculateEnergyTwoDiff( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  ) * factorA * factorB * CISCI_instance%auxcoefficientTarget%values(bb)

        endif
 
      end do bloop !b

      !! calculate diagonal term and denominator of Eq5 10.1063/1.4955109
      diagonal = CISCI_calculateEnergyZero( occA )
      denominator = 1 / ( refEnergy - diagonal ) 
      !! PT2 correction
      energyCorrection = energyCorrection + ( CIenergy**2) * denominator

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
!$  write(*,"(A,E10.3)") "Time for CI-PT2 correction: ", timeB -timeA

  end subroutine CISCI_PT2

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
      auxdecimalNumber = floor( auxdecimalNumber/2.0_8, 16)
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
    do pi = 1, CIcore_instance%numberOfOrbitals%values(spi)
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
    factor = (-1)**n_permu

  end function CISCI_canonicalOrderFactor 

  !! calculate the sign factor for canonical order of the configuration
  function CISCI_getDiffOrbitals ( spi, orbA, orbB, occA, occB, factor ) result (diffOrb) 
    implicit none
    integer, intent(in) :: spi
    type(Ivector), intent(in) :: orbA, orbB, occA, occB
    integer, intent(out) :: factor
    integer :: diffOrb(4), diffPos(4)
    integer :: pi, z

    diffPos = 0
    diffOrb = 0
    z = 0
    ! different orbital in A
    do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
      if ( orbB%values(occA%values(pi) ) == 0  ) then
        z = z + 1
        diffOrb(z) = occA%values(pi)
        diffPos(z) = pi
      endif  
    enddo

    z = 2
    ! different orbital in B
    do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
      if ( orbA%values(occB%values(pi) ) == 0  ) then
        z = z + 1
        diffOrb(z) = occB%values(pi)
        diffPos(z) = pi
      endif  
    enddo

    factor = (-1)**(diffPos(1)-diffPos(3) + diffPos(2) - diffPos(4) )

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
      CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), m)  = orbB(spi)%values(:)
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
        if (  CISCI_instance%confAmplitudeCore(1,i) == -1_1 ) then
          m2b = i-1
          exit
        endif
      enddo
    endif 

    !! Sort according to the index of CI configurations per species in order to find duplicates ( N log( N ) )
    call CISort_quicksort(  CISCI_instance%confAmplitudeCore(:,m1:m2), &
                                     CISCI_instance%index_amplitudeCore%values(m1:m2), & 
                                     1_8, m2b - m1 + 1, n )

    !! Sort amplituted coeff vector according to sorting of index array, keeping both arrays aligned
    call CISort_sortVectorByIndex( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                                   CISCI_instance%index_amplitudeCore%values(m1:m2), &
                                   m2 - m1 + 1 )

    !! merge duplicated configurations (sum amplitudes)
    call CISort_mergeDuplicates ( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                                           CISCI_instance%confAmplitudeCore(:,m1:m2), &
                                           CISCI_instance%index_amplitudeCore%values(m1:m2), &
                                           m2 - m1 + 1 )

    !! reset auxindex arrary, relative position
    do i = 1, m2 - m1 + 1
      CISCI_instance%index_amplitudeCore%values( m1 + i - 1 ) = i
    enddo

    call MTSort ( CISCI_instance%buffer_amplitudeCore%values(m1:m2), &
                  CISCI_instance%index_amplitudeCore%values(m1:m2), &
                  m2 - m1 + 1, "D", 1 )

    call CISort_sortArrayByIndex( CISCI_instance%confAmplitudeCore(:,m1:m2), &
                                  CISCI_instance%index_amplitudeCore%values(m1:m2), &
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
        if ( CISCI_instance%buffer_amplitudeCore%values(i) == 0.0_8 ) then
          auxm = i - 1
          exit
        endif
      enddo
      auxm = auxm  ! change from relative to absolute position

      !! discard the last two quarters of tmp_ampltitude for next run, if not keep it fot PT2 corr
      CISCI_instance%buffer_amplitudeCore%values( halfm + 1 : m2 ) = 0.0_8
      do spi = 1, CIcore_instance%numberOfSpecies 
        CISCI_instance%confAmplitudeCore(CISCI_instance%combinedOrbitalsPositions(1,spi) : CISCI_instance%combinedOrbitalsPositions(2,spi), &
                                         halfm + 1 : m2) = -1_1

      enddo
    endif

  end subroutine CISCI_sortAmplitude


end module CISCI_
