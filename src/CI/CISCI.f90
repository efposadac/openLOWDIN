module CISCI_
  use Exception_
  use Matrix_
  use Vector_
  use MolecularSystem_
  use Configuration_
  use MolecularSystem_
  use String_
  use IndexMap_
  use InputCI_
  use CIcore_
  use CIJadamilu_
  use CIInitial_
  use sort_
  use omp_lib
  use IntVectorSort_
  implicit none

  type, public :: CISCI
    type (Vector8) :: amplitudeCore !old code?
    type (Vector8) :: buffer_amplitudeCore
    type (Vector8) :: auxeigenVector
    type (Vector8) :: coefficientCore
    type (Matrix) :: coefficientTarget
    type (Vector8) :: auxcoefficientTarget

    type (Matrix) :: confCore
    type (Matrix) :: confAmplitudeCore
    type (Matrix) :: saved_confTarget

    type (ivector8) :: index_amplitudeCore

    type (Vector8) :: diagonalCore
    type (Vector8) :: diagonalTarget

    type (Vector8) :: eigenValues
!    type (Vector8), allocatable :: orbitalsCore(:,:) ! species, num of core configurations % num. of orbitals
    integer(8) :: coreSpaceSize
    integer(8) :: targetSpaceSize
    integer(8) :: buffer_amplitudeCoreSize
    real(8) :: PT2energy
    !! store the orbitals for each target configurations, to avoid recomputing them 
    type (IVector), allocatable :: targetOrb(:,:) ! species, num of target configurations % num. of orbitals
  end type CISCI

  type(CISCI) :: CISCI_instance

contains

  subroutine CISCI_show()

        write (6,*) ""
        write (6,"(T2,A62)") "          SELECTED CONFIGURATION INTERACTION (SCI):          " 
        write (6,"(T2,A62)") "                 Adaptive Sampling CI (ASCI)                 " 
        write (6,"(T2,A62)") "                   Deterministic Algorithm                   " 
        write (6,"(T2,A62)") "                  Based on 10.1063/1.4955109                 "

        write(6,*) ""
        write(6,*) "  Diagonalizer for target space hamiltonian : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
        write(6,*) "============================================================="
        write(6,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
        write(6,*) "a software code for computing selected eigenvalues of "
        write(6,*) "large sparse symmetric matrices, "
        write(6,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 
        write(6,*) "============================================================="

        write(6,*) ""
        write(6,*) " Modified sorting algorithm from CENCALC quicksort code "
        write(6,*) "============================================================="
        write(6,*) " Code available at https://github.com/dimassuarez/cencalc_quicksort  "
        write(6,*) " E. Suárez, N. Díaz, J. Méndez and D. Suárez. "
        write(6,*) " CENCALC: A Computational Tool for Conformational Entropy Calculations"
        write(6,*) " from Molecular Simulations."
        write(6,*) " J. Comput. Chem. 54, 2031. DOI: 10.1002/jcc.23350 "
        write(6,*) "============================================================="
        write (6,*) ""
        write (6,"(T2,A,F14.5,A3 )") "Estimated memory needed: ", float(CIcore_instance%numberOfConfigurations*3*8)/(1024**3) , " GB"
        write (6,*) ""

  end subroutine CISCI_show

  !! Allocating arrays 
  subroutine CISCI_constructor()
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    integer :: a,b,c,aa,bb,i
    real(8) :: CIenergy
    integer :: nproc
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    CISCI_instance%coreSpaceSize = CONTROL_instance%CI_SCI_CORE_SPACE
    CISCI_instance%targetSpaceSize = CONTROL_instance%CI_SCI_TARGET_SPACE

    !! the first quarter is the real target space, the second quarter is the waiting list if the amplitude could grow more. 
    !! The last half is just a temporary space to avoid sorting a big array for a single addition
    CISCI_instance%buffer_amplitudeCoreSize = CISCI_instance%targetSpaceSize * 4 

    !! copy and destroying diagonal vector... because of the intial matrix subroutine
    call Vector_constructor8 ( CIcore_instance%diagonalHamiltonianMatrix, &
                              CIcore_instance%numberOfConfigurations, 0.0_8 ) 
    CIcore_instance%diagonalHamiltonianMatrix%values = CIcore_instance%diagonalHamiltonianMatrix2%values
    call Vector_destructor8 ( CIcore_instance%diagonalHamiltonianMatrix2 )


    !!This will carry the index changes after sorting configurations
    call Vector_constructorInteger8 ( CIcore_instance%auxIndexCIMatrix, &
                              CIcore_instance%numberOfConfigurations, 0_8 ) !hmm

    do a = 1, CIcore_instance%numberOfConfigurations
      CIcore_instance%auxIndexCIMatrix%values(a)= a
    end do

    !!  auxiliary array to get the index vector to build a configuration. get the configurations for the hamiltonian matrix in the core and target space
    call CISCI_getInitialIndexes2( CIcore_instance%coreConfigurations, CIcore_instance%coreConfigurationsLevel, CISCI_instance%coreSpaceSize )
    call CISCI_getInitialIndexes2( CIcore_instance%targetConfigurations, CIcore_instance%targetConfigurationsLevel, CISCI_instance%targetSpaceSize )

    !! arrays for CISCI
    call Vector_constructor8 ( CISCI_instance%amplitudeCore, CIcore_instance%numberOfConfigurations,  0.0_8)
    call Vector_constructor8 ( CISCI_instance%coefficientCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)

    !! new arrays
    call Matrix_constructor ( CISCI_instance%confCore, int(numberOfSpecies,8), int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Matrix_constructor ( CISCI_instance%confAmplitudeCore, int(numberOfSpecies,8), int(CISCI_instance%buffer_amplitudeCoreSize,8),  0.0_8)
    call Vector_constructor8 ( CISCI_instance%buffer_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0.0_8)
    call Vector_constructorInteger8 ( CISCI_instance%index_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0_8)
    call Vector_constructor8 ( CISCI_instance%diagonalCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Matrix_constructor ( CISCI_instance%saved_confTarget, int(numberOfSpecies,8), int(CISCI_instance%targetSpaceSize,8),  0.0_8)

    call Matrix_constructor ( CISCI_instance%coefficientTarget, int(CISCI_instance%targetSpaceSize,8), &
                               int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)
    call Vector_constructor8 ( CISCI_instance%auxcoefficientTarget, int(CISCI_instance%targetSpaceSize,8), 0.0_8) !! meh... just to avoid changing everything from matrix to vector format
    call Vector_constructor8 ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8) !! Jadamilu requires to store diagonal vector (in target space)
    call Vector_constructor8 ( CISCI_instance%eigenValues, 15_8, 0.0_8) !! store the eigenvalues per macro iterations


    do a = 1, CISCI_instance%buffer_amplitudeCoreSize
      CISCI_instance%index_amplitudeCore%values(a) = a
    enddo

  end subroutine CISCI_constructor

  !! Allocating arrays 
  subroutine CISCI_constructorNew( numberOfConfigurations )
    implicit none

    type(Configuration) :: auxConfigurationA, auxConfigurationB
    integer(8), intent(out) :: numberOfConfigurations
    integer :: a,b,c,aa,bb,i, spi
    real(8) :: CIenergy
    integer :: nproc
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    CISCI_instance%coreSpaceSize = CONTROL_instance%CI_SCI_CORE_SPACE
    CISCI_instance%targetSpaceSize = CONTROL_instance%CI_SCI_TARGET_SPACE
    numberOfConfigurations = CISCI_instance%targetSpaceSize

    !! the first quarter is the real target space, the second quarter is the waiting list if the amplitude could grow more. 
    !! The last half is just a temporary space to avoid sorting a big array for a single addition
    CISCI_instance%buffer_amplitudeCoreSize = CISCI_instance%targetSpaceSize * 8 

    !! arrays for CISCI
    call Vector_constructor8 ( CISCI_instance%coefficientCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)

    !! new arrays
    call Matrix_constructor ( CISCI_instance%confCore, int(numberOfSpecies,8), int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Matrix_constructor ( CISCI_instance%confAmplitudeCore, int(numberOfSpecies,8), int(CISCI_instance%buffer_amplitudeCoreSize,8),  0.0_8)
    call Vector_constructor8 ( CISCI_instance%buffer_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0.0_8)
    call Vector_constructorInteger8 ( CISCI_instance%index_amplitudeCore, int(CISCI_instance%buffer_amplitudeCoreSize,8),  0_8)
    call Vector_constructor8 ( CISCI_instance%diagonalCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Matrix_constructor ( CISCI_instance%saved_confTarget, int(numberOfSpecies,8), int(CISCI_instance%targetSpaceSize,8),  0.0_8)

    call Matrix_constructor ( CISCI_instance%coefficientTarget, int(CISCI_instance%targetSpaceSize,8), &
                               int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)
    call Vector_constructor8 ( CISCI_instance%auxcoefficientTarget, int(CISCI_instance%targetSpaceSize,8), 0.0_8) !! meh... just to avoid changing everything from matrix to vector format
    call Vector_constructor8 ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8) !! Jadamilu requires to store diagonal vector (in target space)
    call Vector_constructor8 ( CISCI_instance%eigenValues, 15_8, 0.0_8) !! store the eigenvalues per macro iterations


    do a = 1, CISCI_instance%buffer_amplitudeCoreSize
      CISCI_instance%index_amplitudeCore%values(a) = a
    enddo

    allocate ( CISCI_instance%targetOrb ( numberOfSpecies, CISCI_instance%targetSpaceSize ) )


    do a = 1, CISCI_instance%targetSpaceSize
      do spi = 1, numberOfSpecies 
        call Vector_constructorInteger ( CISCI_instance%targetOrb(spi,a), CIcore_instance%numberOfOrbitals%values(spi), 0 ) 
      enddo
    enddo

    !! 
    call IntVectorSort_constructor()

  end subroutine CISCI_constructorNew
 
  
  !! main part
  subroutine CISCI_run()
    use sort_
    implicit none
    integer(8) :: i, j, ii, jj, m,a
    integer :: k ! macro SCI iteration
    integer :: nproc
    real(8) :: timeA(15), timeB(15)
    real(8) :: timeAA, timeBB
    real(8) :: timeAS, timeBS
    type(Vector8) :: eigenValuesTarget
    real(8) :: minValue
    real(8) :: currentEnergy 

    nproc = omp_get_max_threads()
    currentEnergy = HartreeFock_instance%totalEnergy 
    call Vector_constructor8 ( eigenValuesTarget, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    !! HF determinant coefficient
    CISCI_instance%coefficientCore%values(1) = 1.0_8
    CISCI_instance%coefficientTarget%values(1,1) = 1.0_8
    CIcore_instance%eigenVectors%values(1,1) = 1.0_8
    CISCI_instance%eigenValues%values(1) = HartreeFock_instance%totalEnergy 


    write (6,*)    ""
    write (6,"(T2,A29 )")    "Starting SCI macro iterations "
    write (6,*)    ""
    do k = 2, 15
        
!$  timeA(k) = omp_get_wtime()

      !! calculating the amplitudes in core space. This is the pertubation guess of CI eigenvector
      CISCI_instance%amplitudeCore%values = 0.0_8

      call CISCI_core_amplitudes( CISCI_instance%amplitudeCore%values, CIcore_instance%numberOfConfigurations, &
                                  CISCI_instance%coefficientCore%values, CISCI_instance%coreSpaceSize, currentEnergy )

      !! setting the HF again... because the HF amplitude diverges 
      if ( k == 2 ) CISCI_instance%amplitudeCore%values(1) = 1.0_8

      !! resetting the original index changes
      do i = 1, CIcore_instance%numberOfConfigurations
        CIcore_instance%auxIndexCIMatrix%values(i)= i
      end do

      !! copy to keep the unsorted for getting the index later
      CIcore_instance%eigenVectors%values(:,1) = CISCI_instance%amplitudeCore%values

      !! getting the target absolute largest coefficients
      call MTSort ( CISCI_instance%amplitudeCore%values, &
            CIcore_instance%auxIndexCIMatrix%values, CIcore_instance%numberOfConfigurations, "D", nproc )

      !! set all target space coefficint to non zero. This is needed for getting the index, and for assigning an initial weight in the diagonalization
      do i = 1,  CISCI_instance%targetSpaceSize
        if ( abs (CISCI_instance%amplitudeCore%values(i) ) <= CONTROL_instance%CI_MATVEC_TOLERANCE ) then
        ii = CIcore_instance%auxIndexCIMatrix%values(i)
        CIcore_instance%eigenVectors%values(ii,1) = 0.0_8
        end if
      end do


      minValue = abs (CISCI_instance%amplitudeCore%values(CISCI_instance%targetSpaceSize) )
      if ( minValue <= 0.0_8 ) minValue = 0.0_8

      !! recover the configurations for the hamiltonian matrix in the target space
      call CISCI_getInitialIndexes3( CIcore_instance%eigenVectors%values(:,1), minValue, &
                                     CIcore_instance%targetConfigurations, CIcore_instance%targetConfigurationsLevel, CISCI_instance%targetSpaceSize )

      !! storing only the largest diagonal elements (for jadamilu)
      do i = 1,  CISCI_instance%targetSpaceSize
        ii = CIcore_instance%auxIndexCIMatrix%values(i)
        !! storing only the largest diagonal elements (for jadamilu)
        CISCI_instance%diagonalTarget%values(i) = CIcore_instance%diagonalHamiltonianMatrix%values(ii)
      enddo

      !! using the amplitued as the initial coeff guess, after that, use the previous diganolized eigenvectors in target space
      if ( k == 2 ) then 
        do i = 1,  CISCI_instance%targetSpaceSize
          CISCI_instance%coefficientTarget%values(i,1) = CISCI_instance%amplitudeCore%values(i)
        enddo
      else 
        do i = 1,  CISCI_instance%targetSpaceSize
          ii = CIcore_instance%auxIndexCIMatrix%values(i)
          CISCI_instance%coefficientTarget%values(i,1) = CIcore_instance%eigenVectors%values(ii,1) 
        enddo
      end if

      !! eigenvalue guess
      eigenValuesTarget%values(1) = currentEnergy 

      !! diagonalize in target space
      call CISCI_jadamiluInterface( int(CISCI_instance%targetSpaceSize,8), &
                 1_8, &
                 eigenValuesTarget, &
                 CISCI_instance%coefficientTarget, timeAA, timeBB )
  
      !! saving the eigenvectors coeff to an aux vector. Only ground state 
      CISCI_instance%auxcoefficientTarget%values(:) = CISCI_instance%coefficientTarget%values(:,1)

      !! updating the full eigenvector with the solution in the target space
      CIcore_instance%eigenVectors%values(:,1) = 0.0_8
      do i = 1, CISCI_instance%targetSpaceSize
        ii = CIcore_instance%auxIndexCIMatrix%values(i)
        CIcore_instance%eigenVectors%values(ii,1) = CISCI_instance%coefficientTarget%values(i,1)
      end do

      !! convergence criteria
      CISCI_instance%eigenValues%values(k) = eigenValuesTarget%values(1)
      if ( abs( CISCI_instance%eigenValues%values(k) - currentEnergy ) < 1.0E-5 ) then
        write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k , " Energy: ", CISCI_instance%eigenValues%values(k)
        !$  timeB(k) = omp_get_wtime()
        exit
      end if

      !! getting the core absolute largest coefficients
      call MTSort (  CISCI_instance%auxcoefficientTarget%values, &
           CIcore_instance%auxIndexCIMatrix%values,  int( CISCI_instance%coreSpaceSize ,8), "D", nproc )

      !! set all coefficient space coefficint to non zero. This is needed for getting the index
      do i = 1,  CISCI_instance%coreSpaceSize
        ii = CIcore_instance%auxIndexCIMatrix%values(i)
        if ( abs (CISCI_instance%auxcoefficientTarget%values(i) ) <= CONTROL_instance%CI_MATVEC_TOLERANCE ) then
        CIcore_instance%eigenVectors%values(ii,1) = CONTROL_instance%CI_MATVEC_TOLERANCE
        end if
      enddo
      minValue = abs (CISCI_instance%auxcoefficientTarget%values(CISCI_instance%coreSpaceSize) )
      if ( minValue <= CONTROL_instance%CI_MATVEC_TOLERANCE ) minValue = CONTROL_instance%CI_MATVEC_TOLERANCE

      !! recover the configurations for the hamiltonian matrix in the core space
      call CISCI_getInitialIndexes3(  CIcore_instance%eigenVectors%values(:,1), minValue, & 
                                      CIcore_instance%coreConfigurations, CIcore_instance%coreConfigurationsLevel, CISCI_instance%coreSpaceSize )
      !! call CISCI_getInitialIndexes( CIcore_instance%fullConfigurations, CIcore_instance%fullConfigurationsLevel , int(CIcore_instance%numberOfConfigurations,4) )

      !! storing only the largest coefficients, and rearraing the next eigenvector guess 
      do i = 1,  CISCI_instance%coreSpaceSize
        CISCI_instance%coefficientCore%values(i) = CISCI_instance%auxcoefficientTarget%values(i)
      enddo

      !! restart amplitudes for next run
      CISCI_instance%amplitudeCore%values = 0.0_8
      write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k , " Energy: ", CISCI_instance%eigenValues%values(k)

!$  timeB(k) = omp_get_wtime()
      !! updating new reference
      currentEnergy = eigenValuesTarget%values(1)

    enddo !k

    !! summary of the macro iteration 
    write (6,*)    ""
    write (6,"(T2,A95 )")    "                            Selected CI (SCI)  summary                                               "
    write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "
    do k = 2, 15
       write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") k-1,  CISCI_instance%eigenValues%values(k),  &
                                                          CISCI_instance%eigenValues%values(k) - HartreeFock_instance%totalEnergy, &
                                                          CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1), &
                                                          timeB(k) - timeA(k)
      if ( abs( CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1)  ) < 1.0E-5 ) then 
        CIcore_instance%eigenvalues%values(1) = CISCI_instance%eigenValues%values(k) 
        exit
      endif
    enddo !k

    !! calculating PT2 correction. A pertuberd estimation of configurations not include in the target space
    call CISCI_PT2 ( CIcore_instance%eigenVectors%values(:,1), CISCI_instance%amplitudeCore%values, &
                     CISCI_instance%targetSpaceSize, CIcore_instance%numberOfConfigurations, &
                     CISCI_instance%eigenValues%values(k), CISCI_instance%PT2energy )

    write (6,*)    ""
    write (6,"(T2,A,F25.12)") "CI-PT2 energy correction :", CISCI_instance%PT2energy 

  end subroutine CISCI_run

  !! main part
  subroutine CISCI_runNew( eigenVectors )
    use sort_
    implicit none
    type(matrix) :: eigenVectors
    integer(8) :: a, aa, i, j, ii, jj, m
    integer :: k ! macro SCI iteration
    integer :: nproc
    real(8) :: timeA(15), timeB(15)
    real(8) :: timeAA, timeBB
    real(8) :: timeAS, timeBS
    type(Vector8) :: eigenValuesTarget
    real(8) :: minValue
    real(8) :: currentEnergy 

    nproc = omp_get_max_threads()
    currentEnergy = HartreeFock_instance%totalEnergy 
    call Vector_constructor8 ( eigenValuesTarget, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    !! HF determinant coefficient
    CISCI_instance%coefficientCore%values(1) = 1.0_8
    !CISCI_instance%coefficientTarget%values(1,1) = 1.0_8
    CISCI_instance%eigenValues%values(1) = HartreeFock_instance%totalEnergy 

    !call CISCI_initialConfigurations(  CISCI_instance%coefficientCore, CISCI_instance%confCore, m )

    write (6,*)    ""
    write (6,"(T2,A29 )")    "Starting SCI macro iterations "
    write (6,*)    ""

    do k = 2, 15
        
!$  timeA(k) = omp_get_wtime()

      !! calculating the amplitudes in core space. This is the pertubation guess of CI eigenvector
      CISCI_instance%buffer_amplitudeCore%values = 0.0_8
      if ( k == 2 ) call CISCI_initialConfigurations(  CISCI_instance%coefficientCore, CISCI_instance%confCore, m )

      call CISCI_core_amplitudesNew ( CISCI_instance%coefficientCore%values, CISCI_instance%confCore, CISCI_instance%coreSpaceSize, currentEnergy, m )
      m = 1

      if ( k == 2 ) CISCI_instance%buffer_amplitudeCore%values(1) = 1.0_8

      call CISCI_sortAmplitude( ) ! do we need this? again?

      !! saving top-targetSize index configurations sorted from confAmplitudeCore in saved_confTarget. I will gather here the omp results, maybe...
      do i = 1, CISCI_instance%targetSpaceSize
        ii =  CISCI_instance%index_amplitudeCore%values(i) 
        CISCI_instance%saved_confTarget%values(:,i) =  CISCI_instance%confAmplitudeCore%values(:,ii) 
      enddo 

      !! computing the diagonal in the target space, jadamilu requires the diagonal in advance
      call CISCI_buildDiagonalNew ( )

      !! using the amplituded as the initial coeff guess, after that, use the previous diganolized eigenvectors in target space
      do a = 1,  CISCI_instance%targetSpaceSize
        aa = CISCI_instance%index_amplitudeCore%values(a)
        eigenVectors%values(a,1) = CISCI_instance%buffer_amplitudeCore%values(a)
      enddo
      
      !! eigenvalue guess
      eigenValuesTarget%values(1) = currentEnergy 

      !! diagonalize in target space
      call CISCI_jadamiluInterfaceNew( int(CISCI_instance%targetSpaceSize,8), &
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

      !! reset auxindex arrary, for later use in sorting target coeff
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
        CISCI_instance%confCore%values(:,i) = CISCI_instance%saved_confTarget%values(:,ii) 
      enddo
      m = CISCI_instance%coreSpaceSize

      !! reset auxindex array
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo

      !! restart amplitudes for next run, except when exiting to do PT2 corr
      CISCI_instance%buffer_amplitudeCore%values = 0.0_8
      CISCI_instance%confAmplitudeCore%values = 0.0_8
      CISCI_instance%saved_confTarget%values = 0.0_8

      write (6,"(T2,A10,I4,A8,F25.12)")    "SCI Iter: ", k , " Energy: ", CISCI_instance%eigenValues%values(k)

!$  timeB(k) = omp_get_wtime()
      !! updating new reference
      currentEnergy = eigenValuesTarget%values(1)

    enddo !k

    !! summary of the macro iteration 
    write (6,*)    ""
    write (6,"(T2,A95 )")    "                            Selected CI (SCI)  summary                                               "
    write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "
    do k = 2, 15
       write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") k-1,  CISCI_instance%eigenValues%values(k),  &
                                                          CISCI_instance%eigenValues%values(k) - HartreeFock_instance%totalEnergy, &
                                                          CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1), &
                                                          timeB(k) - timeA(k)
      if ( abs( CISCI_instance%eigenValues%values(k) - CISCI_instance%eigenValues%values(k-1)  ) < 1.0E-5 ) then 
        CIcore_instance%eigenvalues%values(1) = CISCI_instance%eigenValues%values(k) 
        exit
      endif
    enddo !k

    !! calculating PT2 correction. A pertuberd estimation of configurations not include in the target space
    call CISCI_PT2New ( CISCI_instance%targetSpaceSize, CISCI_instance%eigenValues%values(k), CISCI_instance%PT2energy )

    write (6,*)    ""
    write (6,"(T2,A,F25.12)") "CI-PT2 energy correction :", CISCI_instance%PT2energy 

  end subroutine CISCI_runNew

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  !! Map the indexes of initial CI matrix to the complete matrix.
  subroutine CISCI_getInitialIndexes( auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize )
    implicit none

    type(imatrix) :: auxConfigurationMatrix
    type(imatrix) :: auxConfigurationLevel
    integer :: auxMatrixSize
    integer(8) :: a,b,c
    integer :: u,v
    integer :: ci
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies, auxnumberOfSpecies
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

    call Matrix_constructorInteger ( auxConfigurationMatrix, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )
    call Matrix_constructorInteger ( auxConfigurationLevel, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )

    !! call recursion
    allocate ( cilevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )

    s = 0
    c = 0
    indexConf = 0
    cilevel = 0

    do ci = 1,  CIcore_instance%sizeCiOrderList 
      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = CISCI_getIndexesRecursion( auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize,  s, numberOfSpecies, indexConf, c, cilevel )
    end do

    deallocate ( indexConf )
    deallocate ( cilevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for getting sorted indexes : ", timeB - timeA ," (s)"

  end subroutine CISCI_getInitialIndexes
  subroutine CISCI_getInitialIndexes2( auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize )
    implicit none

    type(imatrix) :: auxConfigurationMatrix
    type(imatrix) :: auxConfigurationLevel
    integer(8) :: auxMatrixSize
    integer(8) :: a,b,c,cc
    integer :: u, u1, u2 
    integer :: ci, aci
    integer :: i, j, ii, jj 
    integer :: s, ss, numberOfSpecies, auxnumberOfSpecies
    integer :: os,is
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(8), allocatable :: indexConf(:,:)
    integer, allocatable :: cilevel(:,:)
    integer, allocatable :: counter(:,:)
    integer :: ssize
    integer(8) :: x, totalsize, auxtotalsize
    integer(8) :: nn, ncore, chunkSize

!$  timeA = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ncore = omp_get_max_threads()

    call Matrix_constructorInteger ( auxConfigurationMatrix, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )
    call Matrix_constructorInteger ( auxConfigurationLevel, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )

    !! call recursion
    allocate ( cilevel ( numberOfSpecies, ncore ) )
    allocate ( indexConf ( numberOfSpecies, ncore ) )
    allocate ( counter (numberOfSpecies, ncore ) ) 

    indexConf = 0
    cilevel = 0
    counter = 0

    chunkSize = ceiling ( float ( auxMatrixSize ) / float (ncore) )

    !$omp parallel &
    !$omp& private( u, u1, u2, nn ),&
    !$omp& shared( auxConfigurationMatrix, auxConfigurationLevel, cilevel, indexConf, counter)
    !$omp do schedule (static) 
    do nn = 1, ncore
      u1 = ( nn - 1)*chunkSize + 1 
      u2 = ( nn ) * chunkSize 
      if ( u2 > auxMatrixSize ) u2 = auxMatrixSize
    do u = u1, u2
!    do u = 1, auxMatrixSize 

      c = 0
      indexConf(:,nn) = 0
      cilevel(:,nn) = 0

      outer: do aci = 1,  CIcore_instance%sizeCiOrderList 
        cilevel(:,nn) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(aci), :)
        counter(:,nn) = 0 

        totalsize = 1
        do i = 1 , numberOfSpecies
          totalsize = totalsize * CIcore_instance%numberOfStrings(i)%values(cilevel(i,nn) + 1)
        end do
  
        do i = 1 , numberOfSpecies 
          ci = cilevel(i,nn) + 1 
          ssize = CIcore_instance%numberOfStrings2(i)%values(ci)
          indexConf(i,nn) = ssize  + 1
        end do
  
        indexConf(numberOfSpecies,nn) = indexConf(numberOfSpecies,nn) -1
        cc = c

        !! delimiting the index
        if (  CIcore_instance%auxIndexCIMatrix%values(u) > c  .and.  CIcore_instance%auxIndexCIMatrix%values(u) <= c + totalSize   ) then

          !! run over the selected window
          do x = cc + 1 , CIcore_instance%auxIndexCIMatrix%values(u) 
            c = c + 1

            indexConf(numberOfSpecies,nn) = indexConf(numberOfSpecies,nn) + 1
  
            do i = numberOfSpecies, 1 + 1, -1 
              auxtotalsize = 1
              do j = i, numberOfSpecies
                auxtotalsize = auxtotalsize * CIcore_instance%numberOfStrings(j)%values(cilevel(j,nn) + 1)
              end do
              if (counter(i,nn) == auxtotalsize) then
                do j = i, numberOfSpecies
                  ci = cilevel(j,nn) + 1 
                  ssize = CIcore_instance%numberOfStrings2(j)%values(ci)
                  indexConf(j,nn) = ssize + 1
                end do
                counter(i,nn) = 0
                indexConf(i-1,nn) = indexConf(i-1,nn) + 1

              end if
              counter(i,nn) = counter(i,nn) + 1
  
            end do
          end do

          !! saving the index 
          do i = 1, numberOfSpecies
            auxConfigurationMatrix%values(i,u) = indexConf(i,nn) 
            auxConfigurationLevel%values(i,u) = cilevel(i,nn) 
          end do

          exit outer

        end if

        c = cc + totalsize

      end do outer
    end do
    end do
      !$omp end do nowait
      !$omp end parallel

    deallocate (counter)
    deallocate ( indexConf )
    deallocate ( cilevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for getting sorted indexes2 : ", timeB - timeA ," (s)"

  end subroutine CISCI_getInitialIndexes2

  subroutine CISCI_getInitialIndexes3( referenceMatrix, minValue,  auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize )
    implicit none

    real(8), intent(in) :: referenceMatrix(:)
    real(8), intent(in) :: minValue
    type(imatrix) :: auxConfigurationMatrix
    type(imatrix) :: sortedauxindexCIMatrix
    type(imatrix) :: auxConfigurationLevel
    integer(8) :: auxMatrixSize
    integer(8) :: a,b,c,cc
    integer :: u, u1, u2 
    integer :: ci, aci
    integer :: i, j, ii, jj 
    integer :: s, ss, numberOfSpecies, auxnumberOfSpecies
    integer :: os,is
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:)
    integer, allocatable :: counter(:)
    integer :: ssize
    integer(8) :: x, totalsize, auxtotalsize
    integer(8) :: nn, ncore, chunkSize

!$  timeA = omp_get_wtime()

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ncore = omp_get_max_threads()

    call Matrix_constructorInteger ( auxConfigurationMatrix, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )
    call Matrix_constructorInteger ( auxConfigurationLevel, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )

    call Matrix_constructorInteger ( sortedauxindexCIMatrix, int( numberOfSpecies,8), &
          int(auxMatrixSize,8), 0 )

    allocate ( cilevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )
    allocate ( counter (numberOfSpecies ) ) 

    indexConf = 0
    cilevel = 0
    counter = 0
    nn = 1

    c = 0
    indexConf = 0
    cilevel = 0

    outer: do aci = 1,  CIcore_instance%sizeCiOrderList 
      cilevel =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(aci), :)
      counter = 0 

      totalsize = 1
      do i = 1 , numberOfSpecies
        totalsize = totalsize * CIcore_instance%numberOfStrings(i)%values(cilevel(i) + 1)
      end do

      do i = 1 , numberOfSpecies 
        ci = cilevel(i) + 1 
        ssize = CIcore_instance%numberOfStrings2(i)%values(ci)
        indexConf(i) = ssize  + 1
      end do

      indexConf(numberOfSpecies) = indexConf(numberOfSpecies) -1
      cc = c

      !! run over the selected window
      do x = 1, totalSize
        c = c + 1

        indexConf(numberOfSpecies) = indexConf(numberOfSpecies) + 1

        do i = numberOfSpecies, 1 + 1, -1 
          auxtotalsize = 1
          do j = i, numberOfSpecies
            auxtotalsize = auxtotalsize * CIcore_instance%numberOfStrings(j)%values(cilevel(j) + 1)
          end do
          if (counter(i) == auxtotalsize) then
            do j = i, numberOfSpecies
              ci = cilevel(j) + 1 
              ssize = CIcore_instance%numberOfStrings2(j)%values(ci)
              indexConf(j) = ssize + 1
            end do
            counter(i) = 0
            indexConf(i-1) = indexConf(i-1) + 1

          end if
          counter(i) = counter(i) + 1

        end do

       !! search only if the configuration was selected in the auxMatrixSize 
       if ( abs(referenceMatrix(c )) >= minValue) then
          do u = 1, auxMatrixSize  
            if (  c == CIcore_instance%auxIndexCIMatrix%values(u) ) then
              !! saving the index 
              do i = 1, numberOfSpecies
                auxConfigurationMatrix%values(i,u) = indexConf(i) 
                auxConfigurationLevel%values(i,u) = cilevel(i) 
              end do
              exit 
            endif
          end do ! u, auxMatrixSize
        end if

      end do ! x totalsize
       
      c = cc + totalsize

    end do outer

    deallocate (counter)
    deallocate ( indexConf )
    deallocate ( cilevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for getting sorted indexes3 : ", timeB - timeA ," (s)"

  end subroutine CISCI_getInitialIndexes3

recursive  function CISCI_getIndexesRecursion(  auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize, s, numberOfSpecies, indexConf, c, cilevel) result (os)
    implicit none

    type(imatrix) :: auxConfigurationMatrix
    type(imatrix) :: auxConfigurationLevel
    integer :: auxMatrixSize
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
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a
        os = CISCI_getIndexesRecursion( auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize, is, numberOfSpecies, indexConf, c, cilevel)

      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        do u = 1, auxMatrixSize 
          if ( c ==  CIcore_instance%auxIndexCIMatrix%values(u) ) then
            do ss = 1, numberOfSpecies
              auxConfigurationMatrix%values(ss,u) = indexConf(ss) 
              auxConfigurationLevel%values(ss,u) = cilevel(ss) !? check...
            end do
            exit 
          end if
        end do
      end do
    end if

  end function CISCI_getIndexesRecursion

  !! Alternative option to the recursion with the same computational cost... However, it may be helpul some day. 

  function CISCI_getIndexes( sortedIndexVector, auxConfigurationMatrix, auxConfigurationLevel, auxMatrixSize,  c, counter, indexConf, cilevel) result (os)
    implicit none

    type(ivector8) :: sortedIndexVector 
    type(imatrix) :: auxConfigurationMatrix
    type(imatrix) :: auxConfigurationLevel
    integer :: auxMatrixSize
    integer(8) :: a,aa, x, u
    integer :: i, j, ci
    integer :: numberOfSpecies
    integer :: os,is,ss,ssize
    integer(8) :: indexConf(:)
    integer :: cilevel(:)
    integer(8) :: c, totalsize, auxtotalsize
    integer :: counter(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    

    totalsize = 1
    do i = 1 , numberOfSpecies
      totalsize = totalsize * CIcore_instance%numberOfStrings(i)%values(cilevel(i) + 1)
    end do

    do i = 1 , numberOfSpecies 
      ci = cilevel(i) + 1 
      ssize = CIcore_instance%numberOfStrings2(i)%values(ci)
      indexConf(i) = ssize  + 1
    end do

    indexConf(numberOfSpecies) = indexConf(numberOfSpecies) -1

    do x = 1, totalsize
      c = c + 1
      indexConf(numberOfSpecies) = indexConf(numberOfSpecies) + 1

      do i = numberOfSpecies, 1 + 1, -1 
        auxtotalsize = 1
        do j = i, numberOfSpecies
          auxtotalsize = auxtotalsize * CIcore_instance%numberOfStrings(j)%values(cilevel(j) + 1)
        end do
        if (counter(i) == auxtotalsize) then
          do j = i, numberOfSpecies
            ci = cilevel(j) + 1 
            ssize = CIcore_instance%numberOfStrings2(j)%values(ci)
            indexConf(j) = ssize + 1
          end do
          counter(i) = 0
          indexConf(i-1) = indexConf(i-1) + 1
        end if
        counter(i) = counter(i) + 1 

      end do

      do u = 1, auxMatrixSize 
          if ( c ==  sortedIndexVector%values(u) ) then
          do i = 1, numberOfSpecies
            auxConfigurationMatrix%values(i,u) = indexConf(i) 
            auxConfigurationLevel%values(i,u) = cilevel(i) 
          end do
          exit 
        end if
      end do

    end do


  end function CISCI_getIndexes

  subroutine CISCI_core_amplitudes ( amplitudeCore, numberOfConfigurations, coefficientCore, SCICoreSpaceSize, oldEnergy )

    implicit none
  
    integer(8) SCICoreSpaceSize
    integer(8) numberOfConfigurations
    real(8) amplitudeCore ( numberOfConfigurations )
    real(8) coefficientCore ( SCICoreSpaceSize )
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
    integer(8) :: a,b,c, aa, bb
    integer :: s, numberOfSpecies
    integer(8), allocatable :: indexConfA(:) !! ncore, species
    integer, allocatable :: cilevel(:)
    real(8) :: diagEnergy
    real(8) :: oldEnergy
    real(8) :: shift

!$  timeA = omp_get_wtime()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
    shift = 1E-6 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( indexConfA ( numberOfSpecies ) )
    allocate ( cilevel ( numberOfSpecies ) )

    cilevel = 0
    indexConfA = 0

    !$omp parallel &
    !$omp& private( n, a, aa, cilevel, indexConfA, diagEnergy ),&
    !$omp& shared( amplitudeCore, coefficientCore, oldEnergy   ) 
    !$omp do schedule (static) 
    do a = 1, SCICoreSpaceSize  
      aa = CIcore_instance%auxIndexCIMatrix%values(a)
      cilevel = CIcore_instance%coreConfigurationsLevel%values(:,a)
      indexConfA = CIcore_instance%coreConfigurations%values(:,a)
      n = OMP_GET_THREAD_NUM() + 1

      !! using jadamilu subroutine to calculate all configurations coupled to configuration a 
      call CIJadamilu_buildRow( n, indexConfA, aa, amplitudeCore, coefficientCore(a), cilevel  )

      !! removing diagonal term
      amplitudeCore(aa) = amplitudeCore(aa) - CIcore_instance%diagonalHamiltonianMatrix%values(aa) * coefficientCore(a)

    end do
    !$omp end do nowait
    !$omp end parallel

    do b = 1, CIcore_instance%numberOfConfigurations
      amplitudeCore(b) = amplitudeCore(b) / ( CIcore_instance%diagonalHamiltonianMatrix%values(b) - oldEnergy + shift  )
    end do

    CIcore_instance%pindexConf = 0
!$  timeB = omp_get_wtime()
    deallocate ( cilevel )
    deallocate ( indexConfA )
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for calculating SCI amplitudes : ", timeB - timeA ," (s)"

  end subroutine CISCI_core_amplitudes 

  subroutine CISCI_initialConfigurations ( coefficientCore, confCore, m )

    implicit none
    type(vector8) :: coefficientCore
    type(Matrix) :: confCore
    type(IVector), allocatable :: orbA(:)
    integer :: spi, numberOfSpecies
    integer :: pi
    integer(8) :: m 
    real(8) :: indexConf

    m = 1
    !! Hartree-Fock reference
    coefficientCore%values(m) = 1.0_8

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies 
     call Vector_constructorInteger( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0) 
    
      do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
        orbA(spi)%values(pi) = 1.0
      enddo

      call CISCI_binaryToDecimal ( orbA(spi)%values, indexConf )

      confCore%values(spi,m) = indexConf
      
    enddo

    deallocate ( orbA  )

    !!! append the initial configuration to the amplitude, otherwise this is zero at the first iteration 
    CISCI_instance%buffer_amplitudeCore%values(m) = 1.0_8
    CISCI_instance%confAmplitudeCore%values(:,m) = confCore%values(:,m)
    CISCI_instance%index_amplitudeCore%values(m) = 1
    

  end subroutine CISCI_initialConfigurations

  subroutine CISCI_core_amplitudesNew ( coefficientCore, confCore, SCICoreSpaceSize, oldEnergy, m )

    implicit none
  
    integer(8) SCICoreSpaceSize
    real(8) coefficientCore ( SCICoreSpaceSize )
    type(matrix) :: confCore
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
    integer, allocatable :: cilevel(:)
    real(8) :: diagEnergy
    real(8) :: oldEnergy
    real(8) :: shift
    type (ivector), allocatable :: occA(:), occB(:), virA(:), virB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    real(8), allocatable :: confCoreConfB(:)
    real(8) :: tmpconfCoreConfB
    integer :: pi, qi, ri, si, pj, qj, rj, sj
    integer :: oia, oja, via, vja, aaa
    integer :: oi1, vi1, oi2, vi2, oj2, vj2
    integer :: factor1, factor2, factor2j
    integer(8) :: m   !! index to run over buffer_amplitudeCore array
    integer :: nonzero

!$  timeA = omp_get_wtime()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
    shift = 1E-6 !! to avoid divergence
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    !m = 1 !!

    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( confCoreConfB ( numberOfSpecies ) )

    nonzero = 0
    do a = 1, SCICoreSpaceSize  
      nonzero = nonzero + 1
      if ( confCore%values(1,a) == 0 ) exit
    enddo

   !!$omp& shared ( m, CISCI_instance%buffer_amplitudeCore, CISCI_instance%confAmplitudeCore, CISCI_instance%index_amplitudeCore)
   ! !$omp parallel & 
   ! !$omp& firstprivate ( occA, occB, virA, virB, orbA, orbB, confCoreConfB) &
   ! !$omp& private ( a, oia, via, pi, qi, ri, si, oi1, vi1, oi2, vi2, spi, spj, oj2, vj2, factor1, factor2, factor2j, CIenergy, tmpconfCoreConfB  ) &
   ! !$omp& shared ( m, CISCI_instance, coefficientCore, confCore )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do
    confCoreConfB = 0.0_8

   ! !$omp do schedule (static) 
    do a = 1, nonzero !! coreSpace
      
      !print *, a
      if ( confCore%values(1,a) == 0 ) exit
      ! print *, "a", a, confCore%values(:,a)
      !print *, "a",  a, OMP_GET_THREAD_NUM()
      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 
        via = 0

        !! build the orbital from the index using the bit mapping
        call CISCI_decimalToBinary ( confCore%values(spi,a), orbA(spi)%values )

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

        !! copy to conf B
        orbB(spi)%values = orbA(spi)%values 
        occB(spi)%values = occA(spi)%values 
        virB(spi)%values = virA(spi)%values 
      enddo

      confCoreConfB = confCore%values(:,a) !! copy decimal representation of confA 

      !! building all single sustitutions from configuration A. 
      !! here all configurations pairs are generated in maximum coincidence 
      do spi = 1, numberOfSpecies 

        !! calculate the sign factor for canonical order of the configuration
        factor1 = CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )

        do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          oi1 = occA(spi)%values(pi)  
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) - 1 

          do qi = 1, CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            vi1 = virA(spi)%values(qi)
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) + 1
            occB(spi)%values(pi) = vi1

            !! get spingle sustitutions energy
            CIenergy = CISCI_calculateEnergyOneNew( spi, occA, occB, oi1, vi1  )

            !! calculate the sign factor for canonical order of the configuration
            factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )

            CIenergy = CIenergy * factor1 * factor2

            !! bit mapping from orbital to decimal num
            call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )

            if ( CIenergy /= 0.0_8 ) then
              CIenergy = CIenergy * coefficientCore(a) 

              !diagEnergy = CISCI_calculateEnergyZero( occA )
              !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

              !! append the amplitude 
              !!$omp critical (single)
              call CISCI_appendAmplitude ( m, CIenergy, confCoreConfB )
              !!$omp end critical (single)
            endif

            tmpconfCoreConfB = confCoreConfB(spi) !! save the indexconfB to use later in double inter, because double intra will overwritten it 

            !! building all double intraspecies sustitutions from configuration A
            do ri = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
              oi2 = occA(spi)%values(ri)  
              if ( oi1 <= oi2 ) cycle 
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) - 1 
              do si = 1, CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi)
                vi2 = virA(spi)%values(si)
                if ( vi1 <= vi2 ) cycle 
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) + 1
                occB(spi)%values(ri) = vi2 

                !! get double intraspecies sustitutions energy
                CIenergy = CISCI_calculateEnergyTwoSameNew( spi, occA, occB, oi1, oi2, vi1,  vi2 )

                if ( CIenergy /= 0.0_8 ) then
                  !! calculate the sign factor for canonical order of the configuration
                  factor2 = CISCI_canonicalOrderFactor( spi, orbB(spi), occB(spi) )
 
                  CIenergy = CIenergy * factor1 * factor2

                  !! bit mapping from orbital to decimal num
                  call CISCI_binaryToDecimal ( orbB(spi)%values, confCoreConfB(spi) )

                  CIenergy = CIenergy * coefficientCore(a)
                  !diagEnergy = CISCI_calculateEnergyZero( occA )
                  !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                  !! append the amplitude 
                  !!$omp critical (intra)
                  call CISCI_appendAmplitude ( m, CIenergy, confCoreConfB )
                  !!$omp end critical (intra)
                endif

                occB(spi)%values(ri) = occA(spi)%values(ri)  
                orbB(spi)%values(vi2) = orbB(spi)%values(vi2) -1 
              enddo
              orbB(spi)%values(oi2) = orbB(spi)%values(oi2) + 1
            enddo

            !reset indexconf i
            confCoreConfB(spi) = tmpconfCoreConfB 

            !! building all double interspecies sustitutions from configuration A. maybe build this as a superloop?
            !! get double interspecies sustitutions energy
            do spj = spi + 1, numberOfSpecies 
              !if ( spj == spi ) cycle 
              do rj = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                oj2 = occA(spj)%values(rj)  
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) - 1 
                do sj = 1, CIcore_instance%numberOfOrbitals%values(spj) - CIcore_instance%numberOfOccupiedOrbitals%values(spj)
                  vj2 = virA(spj)%values(sj)
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) + 1
                  occB(spj)%values(rj) = vj2 

                  !print *, spi, spj, oi1, oj2, vi1, vj2
                  !! get double interspecies sustitutions energy
                  CIenergy = CISCI_calculateEnergyTwoDiffNew( spi, spj, oi1, oj2, vi1, vj2 )

                  if ( CIenergy /= 0.0_8 ) then

                    !! calculate the sign factor for canonical order of the configuration
                    factor2j = CISCI_canonicalOrderFactor( spj, orbB(spj), occB(spj) )
   
                    CIenergy = CIenergy * factor2 * factor2j

                    !! bit mapping from orbital to decimal num
                    call CISCI_binaryToDecimal ( orbB(spj)%values, confCoreConfB(spj) )

                    CIenergy = CIenergy * coefficientCore(a) 
                    !diagEnergy = CISCI_calculateEnergyZero( occA )
                    !CIenergy = CIenergy / ( diagEnergy - oldEnergy + shift)

                    !! append the amplitude 
                    !!$omp critical (inter)
                    call CISCI_appendAmplitude ( m, CIenergy, confCoreConfB )
                    !!$omp end critical (inter)
                  endif

                  !! reset the confB
                  occB(spj)%values(rj) = occA(spj)%values(rj)  
                  orbB(spj)%values(vj2) = orbB(spj)%values(vj2) -1 
                enddo ! sj
                orbB(spj)%values(oj2) = orbB(spj)%values(oj2) + 1
              enddo ! rj

              confCoreConfB(spj) = confCore%values(spj,a) !! reset just spj index when switchcing another species, spi still unchanged
            enddo !spj

            !! reset the confB
            occB(spi)%values(pi) = occA(spi)%values(pi)  
            orbB(spi)%values(vi1) = orbB(spi)%values(vi1) -1 
          enddo !qi
          orbB(spi)%values(oi1) = orbB(spi)%values(oi1) + 1
        enddo !pi

        confCoreConfB = confCore%values(:,a) !! reset when switchcing another species

      enddo !spi
      
    enddo !enddo a
    !!$omp enddo 
    !!$omp end parallel

    deallocate ( confCoreConfB )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( confCoreConfB ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    !! apply the denominator from eq 4 10.1063/1.4955109. First sort by heat-bath ranking, then final sort by asci ranking
 !   print *, "before denominator"
    call CISCI_sortAmplitude( )

    do aa = 1, CISCI_instance%buffer_amplitudeCoreSize !! replace this by m
      a = CISCI_instance%index_amplitudeCore%values(aa) 

     ! print *, "a den",  CISCI_instance%buffer_amplitudeCore%values(a), CISCI_instance%confAmplitudeCore%values(:,a),  CISCI_instance%index_amplitudeCore%values(a) 
       if (  CISCI_instance%confAmplitudeCore%values(1, a ) == 0 ) cycle
      do spi = 1, numberOfSpecies 
        occA(spi)%values(:) = 0
        virA(spi)%values(:) = 0
        orbA(spi)%values(:) = 0
      enddo

      do spi = 1, numberOfSpecies 

        oia = 0 
        via = 0

        !! build the orbital from the index using the bit mapping
        call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi, a ),  orbA(spi)%values )
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
      !print *, aa, CISCI_instance%confAmplitudeCore%values(:, a ), CISCI_instance%buffer_amplitudeCore%values( aa )
      CISCI_instance%buffer_amplitudeCore%values( aa ) = CISCI_instance%buffer_amplitudeCore%values( aa ) / ( CIenergy - oldEnergy + shift)

    enddo    

    !!! append the core space to the target space (i.e. conf from previous step)
!    do a = 1, SCICoreSpaceSize  
!      if ( confCore%values(1,a) == 0 ) cycle
      !! append the amplitude 
!      call CISCI_appendAmplitude ( m,  coefficientCore(a),  confCore%values(:,a) )
!    enddo



!$  timeB = omp_get_wtime()

    deallocate ( confCoreConfB )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for calculating SCI amplitudes New: ", timeB - timeA ," (s)"

  end subroutine CISCI_core_amplitudesNew

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
     TOL = 1e-3 !1.0d-4 !    tolerance for the eigenvector residual, for ASCI this can be higher

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
    integer :: i, ii, jj, n
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    real(8) :: CIenergy2
    integer(1) :: coupling
    integer(1), allocatable :: orbitalsA(:,:), orbitalsB(:,:), couplingS(:)
    integer :: initialCIMatrixSize 
    integer :: nproc
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()

    allocate ( couplingS ( numberOfSpecies ) )
    allocate ( indexConfA ( numberOfSpecies ) )
    allocate ( indexConfB ( numberOfSpecies ) )
    allocate (orbitalsA ( maxval(CIcore_instance%numberOfOrbitals%values(:)), numberOfSpecies))
    allocate (orbitalsB ( maxval(CIcore_instance%numberOfOrbitals%values(:)), numberOfSpecies))

    nonzero = 0
    nonzerow = 0
    w = 0.0_8 
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 
    do a = 1 , nx
       if ( abs(v(a) ) >= tol) nonzero = nonzero + 1
    end do

!$    timeA= omp_get_wtime()

    !$omp parallel &
    !$omp& private( a, aa, b, bb, uu, vv, coupling, couplingS, CIenergy, i, ii, jj, orbitalsA, orbitalsB),&
    !$omp& firstprivate( indexConfA, indexConfB ),&
    !$omp& shared( v, nx, numberOfSpecies, tol )  reduction (+:w)
    !$omp do schedule (dynamic)
    do a = 1, nx

      indexConfA =  CIcore_instance%targetConfigurations%values(:,a)

      orbitalsA = 0
      do i = 1, numberOfSpecies
        do uu = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
          orbitalsA( CIcore_instance%strings(i)%values(uu,indexConfA(i) ), i ) = 1
        end do
      end do

      do b = a, nx

        couplingS = 0
        indexConfB(:) =  CIcore_instance%targetConfigurations%values(:,b)

        orbitalsB = 0
        do i = 1, numberOfSpecies
          do vv = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
            orbitalsB( CIcore_instance%strings(i)%values(vv,indexConfB(i) ), i ) = 1
          end do
        end do

        do i = 1, numberOfSpecies
          couplingS(i) = couplingS(i) + &
            CIcore_instance%numberOfOccupiedOrbitals%values(i) - sum ( orbitalsA(:,i) * orbitalsB(:,i) ) + 1
        end do

        coupling = product(couplingS)

        select case (coupling)

        case(1)
            CIenergy = CISCI_instance%diagonalTarget%values(a)
            w(a) = w(a) + CIEnergy*v(a)

        case(2)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) ii = i
          end do
          CIenergy = CISCI_calculateEnergyOne ( ii, indexConfA, indexConfB )
          w(b) = w(b) + CIenergy * v(a)
          w(a) = w(a) + CIenergy * v(b)

        case(3)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 3 ) ii = i
          end do

          CIenergy = CISCI_calculateEnergyTwoSame ( ii, indexConfA, indexConfB ) 
          w(b) = w(b) + CIenergy * v(a)
          w(a) = w(a) + CIenergy * v(b)

        case(4)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) then 
              ii = i
              exit
            end if
          end do
          do i = ii+1, numberOfSpecies
            if ( couplingS(i) == 2 ) jj = i
          end do
          CIenergy = CISCI_calculateEnergyTwoDiff ( ii, jj, indexConfA, indexConfB ) 
          w(b) = w(b) + CIenergy * v(a)
          w(a) = w(a) + CIenergy * v(b)

        end select
 
      end do !b
    end do !a 
    !$omp end do nowait 
    !$omp end parallel

!$  timeB = omp_get_wtime()
    !! to check how dense is the w vector
    do a = 1 , nx
       if ( abs(w(a) ) >= tol) nonzerow = nonzerow + 1
    end do
    deallocate (orbitalsA )
    deallocate (orbitalsB )
    deallocate ( indexConfB )
    deallocate ( indexConfA )
    deallocate ( couplingS )

!$    write(*,"(A,I2,A,E10.3,A2,I12,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero, nonzerow
    return

  end subroutine CISCI_matvec

  subroutine CISCI_jadamiluInterfaceNew(n,  maxeig, eigenValues, eigenVectors, timeA, timeB)
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
     TOL = 1e-3 !1.0d-4 !    tolerance for the eigenvector residual, for ASCI this can be higher

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
        call CISCI_matvecNew ( N, X(NDX1), X(NDX2), iiter)
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

  end subroutine CISCI_jadamiluInterfaceNew

  subroutine CISCI_matvecNew ( nx, v, w, iter)
  
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
    integer :: initialCIMatrixSize 
    integer :: nproc
    integer :: diffOrbi(4)
    integer :: diffOrbj(4)
    integer :: pi
    integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA, factorB

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()

    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    nonzero = 0
    nonzerow = 0
    w = 0.0_8 
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 
    do a = 1 , nx
       if ( abs(v(a) ) >= tol) nonzero = nonzero + 1
    end do

    !! transforming from decimal to binary once
    do aa = 1, nonzero
      do spi = 1, numberOfSpecies
        call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,aa), CISCI_instance%targetOrb(spi,aa)%values )
      enddo
    enddo

!$    timeA= omp_get_wtime()

    aloop: do aa = 1, nonzero

      !a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      a = aa ! if index_amplitude is sorted

      !if ( CISCI_instance%saved_confTarget%values(1,a) == 0 ) exit !using nonzero now

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        !! build the orbital from the index using the bit mapping
        !call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
        orbA(spi)%values = CISCI_instance%targetOrb(spi,aa)%values

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
        !if ( CISCI_instance%saved_confTarget%values(1,b) == 0 ) exit
        ! getting configuration B
        do spi = 1, numberOfSpecies 

          !oib = 0 

          !! build the orbital from the index using the bit mapping
         ! call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,b), orbB(spi)%values )
          orbB(spi)%values = CISCI_instance%targetOrb(spi,b)%values
          !!call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi,b), orbB(spi)%values )

        enddo

        !! determinate number of diff orbitals
        couplingS = 0
        do spi = 1, numberOfSpecies
          couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                            - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) + 1
          !if ( couplingS(spi) > 2 ) cycle bloop
        end do
      
        !! if any species differs in more than 2 orb, skip
        coupling = product(couplingS)

        if ( coupling <= 4 ) then
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

        select case (coupling)

        case(1)
            CIenergy = CISCI_instance%diagonalTarget%values(aa) 
            w(aa) = w(aa) + CIenergy * v(aa)

        case(2)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CISCI_calculateEnergyOneNew( spi, occA, occB, diffOrbi(1), diffOrbi(3)  )

          w(bb) = w(bb) + CIenergy * v(aa) * factorA 
          w(aa) = w(aa) + CIenergy * v(bb) * factorA

        case(3)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 3 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CISCI_calculateEnergyTwoSameNew( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  )

          w(bb) = w(bb) + CIenergy * v(aa) * factorA 
          w(aa) = w(aa) + CIenergy * v(bb) * factorA 

        case(4)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) then 
              spi = i
              exit
            end if
          end do

          do i = spi+1, numberOfSpecies
            if ( couplingS(i) == 2 ) spj = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj), orbB(spj), occA(spj), occB(spj), factorB )
          CIenergy = CISCI_calculateEnergyTwoDiffNew( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  )

          w(bb) = w(bb) + CIenergy * v(aa) * factorA * factorB
          w(aa) = w(aa) + CIenergy * v(bb) * factorA * factorB

        end select
 
      end do bloop !b
    end do aloop !a 

!$  timeB = omp_get_wtime()

    !! to check how dense is the w vector
    do a = 1 , nx
       if ( abs(w(a) ) >= tol) nonzerow = nonzerow + 1
    end do

    deallocate ( couplingS )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

!$    write(*,"(A,I2,A,E10.3,A2,I12,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero, nonzerow
    return

  end subroutine CISCI_matvecNew

  subroutine CISCI_buildDiagonalNew ( )
    implicit none
    integer(8) :: a,b,aa,bb
    integer :: i, ii, jj, n, spi, spj
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    real(8) :: CIenergy
    integer(1) :: coupling
    integer :: initialCIMatrixSize 
    integer :: nproc
    integer :: pi
    integer :: oia
    type (ivector), allocatable :: occA(:)
    type (ivector), allocatable :: orbA(:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()

    allocate ( occA ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

!$    timeA= omp_get_wtime()

    do aa = 1, CISCI_instance%targetSpaceSize 

      !a = CISCI_instance%index_amplitudeCore%values(aa)
      a = aa
      if ( CISCI_instance%saved_confTarget%values(1,a) == 0 ) exit
      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 
        !! build the orbital from the index using the bit mapping
        call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
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

  end subroutine CISCI_buildDiagonalNew

  function CISCI_calculateEnergyOne( ii, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: a, b
    integer :: i,j,s,n, nn,ii
    integer :: l,k,z,kk,ll
    integer :: factor, factor2, auxOcc, AA, BB
    logical(1) :: equalA, equalB
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex11, auxIndex2, auxIndex
    integer :: diffOrb(2), otherdiffOrb(2) !! to avoid confusions
    real(8) :: auxCIenergy

    auxCIenergy = 0.0_8
    factor = 1
  
    !! copy a
    a = thisA(ii)
    b = thisB(ii)
  
    diffOrb = 0
  
    do kk = 1, CIcore_instance%occupationNumber(ii) 
      if ( CIcore_instance%orbitals(ii)%values( &
             cicore_instance%strings(ii)%values(kk,a),b) == 0 ) then
        diffOrb(1) =  CIcore_instance%strings(ii)%values(kk,a)
        AA = kk
        exit
      end if
    end do
  
    do kk = 1, CIcore_instance%occupationNumber(ii) 
      if ( CIcore_instance%orbitals(ii)%values( &
             CIcore_instance%strings(ii)%values(kk,b),a) == 0 ) then
        diffOrb(2) =  CIcore_instance%strings(ii)%values(kk,b)
        BB = kk
        exit
      end if
    end do
  
    factor = (-1)**(AA-BB)
    !One particle terms
  
    auxCIenergy= auxCIenergy +  CIcore_instance%twoCenterIntegrals(ii)%values( diffOrb(1), diffOrb(2) )
  
    !! save the different orbitals
  
    auxIndex1= CIcore_instance%twoIndexArray(ii)%values( diffOrb(1), diffOrb(2))
  
    do ll=1, CIcore_instance%occupationNumber( ii ) !! the same orbitals pair are excluded by the exchange
  
      l = CIcore_instance%strings(ii)%values(ll,b) !! or a
  
      auxIndex2 = CIcore_instance%twoIndexArray(ii)%values( l,l) 
      auxIndex = CIcore_instance%fourIndexArray(ii)%values( auxIndex1, auxIndex2 )
  
      auxCIenergy = auxCIenergy + CIcore_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)
  
      auxIndex = CIcore_instance%fourIndexArray(ii)%values( &
                             CIcore_instance%twoIndexArray(ii)%values(diffOrb(1),l), &
                             CIcore_instance%twoIndexArray(ii)%values(l,diffOrb(2)) ) 
  
      auxCIenergy = auxCIenergy + &
                     MolecularSystem_instance%species(ii)%kappa*CIcore_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)
  
    end do
  
    !end if
    do j=1, ii - 1 !! avoid ii, same species

      b = thisB(j)

      auxnumberOfOtherSpecieSpatialOrbitals = CIcore_instance%numberOfSpatialOrbitals2%values(j) 
      auxIndex11 = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) 

      do ll=1,  CIcore_instance%occupationNumber( j ) 

        l = CIcore_instance%strings(j)%values(ll,b)

        auxIndex = auxIndex11  + CIcore_instance%twoIndexArray(j)%values( l,l) 

        auxCIenergy = auxCIenergy + &
        CIcore_instance%fourCenterIntegrals(ii,j)%values(auxIndex, 1) 

      end do

    end do

    do j= ii + 1, MolecularSystem_instance%numberOfQuantumSpecies!! avoid ii, same species

      b = thisB(j)

      auxnumberOfOtherSpecieSpatialOrbitals = CIcore_instance%numberOfSpatialOrbitals2%values(j) 

      auxIndex11 = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) 

      do ll=1,  CIcore_instance%occupationNumber( j )

        l = CIcore_instance%strings(j)%values(ll,b)

        auxIndex = auxIndex11  + CIcore_instance%twoIndexArray(j)%values( l,l) 

        auxCIenergy = auxCIenergy + &
        CIcore_instance%fourCenterIntegrals(ii,j)%values(auxIndex, 1) 
      end do

    end do

    auxCIenergy= auxCIenergy * factor

  end function CISCI_calculateEnergyOne

  function CISCI_calculateEnergyOneNew( si, occA, occB, a, b ) result ( CIenergy )
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

  end function CISCI_calculateEnergyOneNew

  function CISCI_calculateEnergyTwoSame( ii, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: a, b
    integer :: ii
    integer :: kk,z
    integer :: factor, AA(2), BB(2)
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: auxIndex
    integer :: diffOrbA(2), diffOrbB(2)  !! to avoid confusions
    real(8) :: auxCIenergy

    a = thisA(ii)
    b = thisB(ii)
    z = 0
    auxCIenergy = 0.0_8

      do kk = 1, CIcore_instance%occupationNumber(ii) 
        if ( CIcore_instance%orbitals(ii)%values( &
               CIcore_instance%strings(ii)%values(kk,a),b) == 0 ) then
          z = z + 1
          diffOrbA(z) = CIcore_instance%strings(ii)%values(kk,a)
          AA(z) = kk
          if ( z == 2 ) exit
        end if
      end do
  
      z = 0
      do kk = 1, CIcore_instance%occupationNumber(ii) 
        if ( CIcore_instance%orbitals(ii)%values( &
               CIcore_instance%strings(ii)%values(kk,b),a) == 0 ) then
          z = z + 1
          diffOrbB(z) =  CIcore_instance%strings(ii)%values(kk,b)
          BB(z) = kk
          if ( z == 2 ) exit
        end if
      end do
  
      factor = (-1)**(AA(1)-BB(1) + AA(2) - BB(2) )
      auxIndex = CIcore_instance%fourIndexArray(ii)%values( &
                  CIcore_instance%twoIndexArray(ii)%values(&
                    diffOrbA(1),diffOrbB(1)),&
                  CIcore_instance%twoIndexArray(ii)%values(&
                    diffOrbA(2),diffOrbB(2)) )
  
      auxCIenergy = CIcore_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)
  
      auxIndex = CIcore_instance%fourIndexArray(ii)%values( &
                   CIcore_instance%twoIndexArray(ii)%values(&
                     diffOrbA(1),diffOrbB(2)),&
                   CIcore_instance%twoIndexArray(ii)%values(&
                     diffOrbA(2),diffOrbB(1)) )
      auxCIenergy = auxCIenergy + &
                  MolecularSystem_instance%species(ii)%kappa*CIcore_instance%fourCenterIntegrals(ii,ii)%values(auxIndex, 1)
  
      auxCIenergy= auxCIenergy * factor

  end function CISCI_calculateEnergyTwoSame
  
  function CISCI_calculateEnergyTwoSameNew( si, occA, occB, ai, aj, bi, bj ) result ( CIenergy )
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
  
  end function CISCI_calculateEnergyTwoSameNew

  function CISCI_calculateEnergyTwoDiff( ii, jj, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: a, b
    integer :: ii, jj
    integer :: kk,z
    integer :: factori, factorj, AA, BB
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: auxIndex, auxIndex1, auxIndex2
    integer :: diffOrb(2)
    real(8) :: auxCIenergy

    a = thisA(ii)
    b = thisB(ii)
  
    diffOrb = 0
  
    do kk = 1, CIcore_instance%occupationNumber(ii) 
      if ( CIcore_instance%orbitals(ii)%values( &
             CIcore_instance%strings(ii)%values(kk,a),b) == 0 ) then
        diffOrb(1) =  CIcore_instance%strings(ii)%values(kk,a)
        AA = kk
        exit
      end if
    end do
  
    do kk = 1, CIcore_instance%occupationNumber(ii) 
      if ( CIcore_instance%orbitals(ii)%values( &
             CIcore_instance%strings(ii)%values(kk,b),a) == 0 ) then
        diffOrb(2) =  CIcore_instance%strings(ii)%values(kk,b)
        BB = kk
        exit
      end if
    end do
  
    factori = (-1)**(AA-BB)
    auxIndex1= CIcore_instance%twoIndexArray(ii)%values( diffOrb(1), diffOrb(2))
    auxIndex1 = CIcore_instance%numberOfSpatialOrbitals2%values(jj) * (auxIndex1 - 1 ) 

    a = thisA(jj)
    b = thisB(jj)
 
    diffOrb = 0
  
    do kk = 1, CIcore_instance%occupationNumber(jj) 
      if ( CIcore_instance%orbitals(jj)%values( &
             CIcore_instance%strings(jj)%values(kk,a),b) == 0 ) then
        diffOrb(1) =  CIcore_instance%strings(jj)%values(kk,a)
        AA = kk
        exit
      end if
    end do
  
    do kk = 1, CIcore_instance%occupationNumber(jj) 
      if ( CIcore_instance%orbitals(jj)%values( &
             CIcore_instance%strings(jj)%values(kk,b),a) == 0 ) then
        diffOrb(2) =  CIcore_instance%strings(jj)%values(kk,b)
        BB = kk
        exit
      end if
    end do
  
    factorj = (-1)**(AA-BB)


    auxIndex2= CIcore_instance%twoIndexArray(jj)%values( diffOrb(1), diffOrb(2))
    auxIndex = auxIndex1 + auxIndex2

    auxCIenergy = factori * factorj *CIcore_instance%fourCenterIntegrals(ii,jj)%values(auxIndex, 1)
 
  end function CISCI_calculateEnergyTwoDiff

  function CISCI_calculateEnergyTwoDiffNew( si, sj, ai, aj, bi, bj ) result ( CIenergy )
    implicit none
    integer, intent(in) :: ai, aj, bi, bj
    integer, intent(in) :: si, sj
    integer(8) :: aibi,  aux_aibi, ajbj
    real(8) :: CIenergy

    CIenergy = 0.0_8

    aibi = CIcore_instance%twoIndexArray(si)%values( ai, bi )
    aux_aibi = CIcore_instance%numberOfSpatialOrbitals2%values( sj ) * ( aibi - 1 ) 

    ajbj = CIcore_instance%twoIndexArray(sj)%values( aj, bj )
    !print *, si, sj, ai, aj, bi, bj, CIenergy 
    CIenergy = CIcore_instance%fourCenterIntegrals( si, sj )%values( aux_aibi +  ajbj, 1)


  end function CISCI_calculateEnergyTwoDiffNew

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

  subroutine CISCI_PT2 ( coefficients, auxenergyCorrection, SCITargetSpaceSize, numberOfConfigurations, refEnergy, energyCorrection )
    implicit none
    integer(8) :: SCITargetSpaceSize
    integer(8) numberOfConfigurations
    real(8) :: coefficients ( numberOfConfigurations )
    real(8) :: auxenergyCorrection ( numberOfConfigurations )
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
    integer(8) :: a,b,c, aa, bb
    integer :: s, numberOfSpecies
    integer(8), allocatable :: indexConfA(:) !! ncore, species
    integer, allocatable :: cilevel(:)
    real(8) :: diagEnergy
    real(8) :: oldEnergy
    real(8) :: shift

    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
!$  timeA = omp_get_wtime()
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( indexConfA ( numberOfSpecies ) )
    allocate ( cilevel ( numberOfSpecies ) )

    cilevel = 0
    indexConfA = 0
    energyCorrection = 0.0_8
    auxenergyCorrection = 0.0_8

    do a = 1, SCITargetSpaceSize  

      aa = CIcore_instance%auxIndexCIMatrix%values(a)
      cilevel = CIcore_instance%targetConfigurationsLevel%values(:,a)
      indexConfA = CIcore_instance%targetConfigurations%values(:,a)
      n = OMP_GET_THREAD_NUM() + 1

      !! using jadamilu subroutine to calculate all configurations coupled to configuration a 
      call CIJadamilu_buildRow( n, indexConfA, aa, auxenergyCorrection, coefficients(aa), cilevel  )

      !! removing the contributions from configurations within the target space
      do b = 1, SCITargetSpaceSize  
        bb = CIcore_instance%auxIndexCIMatrix%values(b)
        auxenergyCorrection(bb) = 0.0_8 !! commented to force test
      enddo

    end do

    do b = 1, numberOfConfigurations
      energyCorrection = energyCorrection + auxenergyCorrection(b) **2 / ( refEnergy - CIcore_instance%diagonalHamiltonianMatrix%values(b)  )
    enddo


    CIcore_instance%pindexConf = 0
!$  timeB = omp_get_wtime()
    deallocate ( cilevel )
    deallocate ( indexConfA )
!$    write(*,"(A,E10.3)") "Time for CI-PT2 correction: ", timeB -timeA

  end subroutine CISCI_PT2

  subroutine CISCI_PT2New ( SCITargetSpaceSize, refEnergy, energyCorrection )
  
    implicit none
    type(matrix) :: confCore !?? or index amplityde?
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
    integer, allocatable :: cilevel(:)
    real(8) :: diagEnergy
    real(8) :: shift
    type (ivector), allocatable :: occA(:), occB(:), virA(:), virB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer(1), allocatable :: couplingS(:)
    real(8) :: indexConf
    real(8), allocatable :: confCoreConfA(:)
    real(8) :: tmpconfCoreConfB
    integer :: pi
    integer :: oia, oib
    integer :: factorA, factorB
    real(8) :: diagonal, denominator
    logical :: found
    integer(8) :: m   !! index to run over buffer_amplitudeCore array
    integer(1) :: coupling
    integer :: diffOrbi(4)
    integer :: diffOrbj(4)

!$  timeA = omp_get_wtime()
    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
!    shift = 1E-6 !! to avoid divergence
    shift = 0.0_8
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( virA ( numberOfSpecies ) )
    allocate ( virB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( confCoreConfA ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( virA(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( virB(spi), CIcore_instance%numberOfOrbitals%values(spi) - CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )  
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfOrbitals%values(spi),  0 ) 
    end do

    energyCorrection = 0.0_8
    
    aloop: do aa = CISCI_instance%targetSpaceSize + 1, CISCI_instance%buffer_amplitudeCoreSize

      a = CISCI_instance%index_amplitudeCore%values(aa) ! if index_amplitude is unsortered
      !a = aa ! if index_amplitude is sorted
      if (CISCI_instance%confAmplitudeCore%values(1,a) == 0) cycle

      ! getting configuration A
      do spi = 1, numberOfSpecies 

        oia = 0 

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,a), orbA(spi)%values )
        call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi,a), orbA(spi)%values )

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

        if (CISCI_instance%saved_confTarget%values(1,b) == 0) exit

        ! getting configuration B
        do spi = 1, numberOfSpecies 

          oib = 0 

          !! build the orbital from the index using the bit mapping
          call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,b), orbB(spi)%values )
          !call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi,b), orbB(spi)%values )

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
                            - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) + 1
          !if (couplingS(spi) > 2 ) cycle bloop
        end do

        coupling = product(couplingS)

        select case (coupling)

        case(2)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CIenergy + CISCI_calculateEnergyOneNew( spi, occA, occB, diffOrbi(1), diffOrbi(3)  ) * factorA * CISCI_instance%auxcoefficientTarget%values(bb)


        case(3)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 3 ) spi = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          CIenergy = CIenergy + CISCI_calculateEnergyTwoSameNew( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  ) * factorA * CISCI_instance%auxcoefficientTarget%values(bb)

        case(4)
          do i = 1, numberOfSpecies
            if ( couplingS(i) == 2 ) then 
              spi = i
              exit
            end if
          end do
          do i = ii+1, numberOfSpecies
            if ( couplingS(i) == 2 ) spj = i
          end do

          diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
          diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj), orbB(spj), occA(spj), occB(spj), factorB )
          CIenergy = CIenergy + CISCI_calculateEnergyTwoDiffNew( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  ) * factorA * factorB * CISCI_instance%auxcoefficientTarget%values(bb)


        end select
 
      end do bloop !b

      !! calculate diagonal term and denominator of Eq5 10.1063/1.4955109
      diagonal = CISCI_calculateEnergyZero( occA )
      denominator = 1 / ( refEnergy - diagonal ) 
      !! PT2 correction
      energyCorrection = energyCorrection + ( CIenergy**2) * denominator
      !print *, a, energyCorrection

    end do aloop !a 

 
!$  timeB = omp_get_wtime()

    deallocate ( couplingS )
    deallocate ( confCoreConfA )
    deallocate ( occA  )
    deallocate ( occB  )
    deallocate ( virA  )
    deallocate ( virB  )
    deallocate ( orbA  )
    deallocate ( orbB  )

!$    write(*,"(A,E10.3)") "Time for CI-PT2 correction: ", timeB -timeA

  end subroutine CISCI_PT2New

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
  function CISCI_maximumCoincidenceFactor ( spi, orbA, orbB ) result (factor) 
    implicit none
    integer, intent(in) :: spi
    type(Ivector), intent(in) :: orbA, orbB
    integer :: factor 
    integer :: n_permu, pi

    n_permu = 0
    do pi = 1, CIcore_instance%numberOfOccupiedOrbitals%values(spi)
      n_permu = n_permu + ( orbA%values(pi) - orbB%values(pi) ) * pi
    end do
    factor = (-1)**n_permu

  end function CISCI_maximumCoincidenceFactor 



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

  subroutine CISCI_appendAmplitude ( m, amplitude, indexConf )
    implicit none
    integer(8), intent(inout) :: m
    real(8), intent(in) :: amplitude
    real(8), intent(in) :: indexConf(:)
    real(8) :: diffConf
    integer(8) :: newm
    integer :: i, j, ii, sp
    logical :: addConf
    integer :: nproc

    !nproc = omp_get_max_threads()
    nproc = 1

    !! add the configuration at the end of the array
    m = m + 1
    CISCI_instance%buffer_amplitudeCore%values(m) = amplitude
    CISCI_instance%confAmplitudeCore%values(:,m) = indexConf(:)
    CISCI_instance%index_amplitudeCore%values(m) = m

     ! print *, "m", m,  CISCI_instance%buffer_amplitudeCore%values(m), CISCI_instance%confAmplitudeCore%values(:,m),  CISCI_instance%index_amplitudeCore%values(m) 

    !! sort the array when it's full
    if ( m == CISCI_instance%buffer_amplitudeCoreSize ) then
      !print *, "sorting"
      call CISCI_sortAmplitude( m )
      !print *, "sorted"
    endif

  end subroutine CISCI_appendAmplitude

  subroutine CISCI_sortAmplitude( auxm )
    implicit none
    integer(8), intent(inout), optional  :: auxm
    integer(8)  :: m
    real(8) :: diffConf
    integer(8) :: newm
    integer :: i, j, ii, sp
    logical :: addConf
    integer :: nproc

    !nproc = omp_get_max_threads()
    nproc = 1
    if ( present ( auxm ) ) then 
      m = auxm
    else 
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        if (  CISCI_instance%confAmplitudeCore%values(1,i) == 0 ) then
          m = i
          exit
        endif
      enddo
    endif

!      print *, "-----------------------------", m

   !   print *, "before finish new sort 1"
   !  do i = 1, m 
   !     print *, i, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i)
   !   enddo

      !! Sort according to the index of CI configurations per species in order to find duplicates ( N log( N ) )
      call IntVectorSort_quicksort(  CISCI_instance%confAmplitudeCore, &
                                     CISCI_instance%index_amplitudeCore, & 
                                     1_8, m )
                                     !1_8,  CISCI_instance%buffer_amplitudeCoreSize )
   !   print *, "after finish new sort 1"
     ! do i = 1, m 
     !   print *, i, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i)
     ! enddo

  
      !! Sort amplituted coeff vector according to sorting of index array, keeping both arrays aligned
      call CISCI_sortVectorByIndex( CISCI_instance%buffer_amplitudeCore, CISCI_instance%index_amplitudeCore, CISCI_instance%buffer_amplitudeCoreSize )

     ! print *, "after sorting the coeff vector"
     ! do i = 1, m 
     !   print *, i, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i)
     ! enddo

    !! merge duplicated configurations (sum amplitudes)
      call IntVectorSort_mergeDuplicates ( CISCI_instance%buffer_amplitudeCore, &
                                           CISCI_instance%confAmplitudeCore, &
                                           CISCI_instance%index_amplitudeCore, &
                                           CISCI_instance%buffer_amplitudeCoreSize )

     ! print *, "after merging duplicates"
     ! do i = 1, m 
     !   print *, i, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i)
     ! enddo

      !! reset auxindex arrary
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo

      call MTSort ( CISCI_instance%buffer_amplitudeCore%values, &
            CISCI_instance%index_amplitudeCore%values,  CISCI_instance%buffer_amplitudeCoreSize, "D", nproc )

     ! print *, "finish sort 1"
     ! do i = 1, m 
     !   print *, i, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
     ! enddo


      call CISCI_sortArrayByIndex( CISCI_instance%confAmplitudeCore, CISCI_instance%index_amplitudeCore, CISCI_instance%buffer_amplitudeCoreSize )

     ! print *, "finish sort index1"
     ! do i = 1, m 
     !  print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
     ! enddo

      !! reset auxindex arrary
      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
        CISCI_instance%index_amplitudeCore%values( i ) = i
      enddo
!
!      !! the next configurations will be added in the last half
!      m = CISCI_instance%buffer_amplitudeCoreSize
!      !m = CISCI_instance%buffer_amplitudeCoreSize/2 
!  !    print *, "after index reset"
!
!      do i = 1, m 
! !      print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
!      enddo
!
!      !! check if the configuration is already in the amplitudes set
!      do i = 1, m
!        !print *, i
!        if ( CISCI_instance%buffer_amplitudeCore%values(i) == 0.0_8 ) cycle
!
!        do j = i + 1, m
!
!          diffConf = 0.0_8
!          do sp = 1, CIcore_instance%numberOfQuantumSpecies 
!            diffConf = diffConf + abs( CISCI_instance%confAmplitudeCore%values(sp, CISCI_instance%index_amplitudeCore%values(i) ) &
!                                    - CISCI_instance%confAmplitudeCore%values(sp, CISCI_instance%index_amplitudeCore%values(j) ) )
!          enddo
!
!          if ( diffConf == 0 ) then
!            !! add and clean
!            CISCI_instance%buffer_amplitudeCore%values(i) = CISCI_instance%buffer_amplitudeCore%values(i) + CISCI_instance%buffer_amplitudeCore%values(j)
!            CISCI_instance%buffer_amplitudeCore%values(j) = 0.0_8
!            CISCI_instance%confAmplitudeCore%values(:,j) = 0.0_8
!            !CISCI_instance%index_amplitudeCore%values(j) = j
!          end if
!        enddo !j
!      enddo !i
!
!      print *, "finish duplicate"
!      do i = 1, m 
!!        print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
!      enddo
!
!      !! sort again removing duplicates
!      call MTSort ( CISCI_instance%buffer_amplitudeCore%values, &
!            CISCI_instance%index_amplitudeCore%values,  CISCI_instance%buffer_amplitudeCoreSize, "D", nproc )
!
!     print *, "finish sort 2"
!      do i = 1, m 
!!        print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
!      enddo
!
!      call CISCI_sortArrayByIndex( CISCI_instance%confAmplitudeCore, CISCI_instance%index_amplitudeCore, CISCI_instance%buffer_amplitudeCoreSize )
!      print *, "finish sort 3"
!
!      do i = 1, m 
!       ! print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
!      enddo
!
!      !! reset auxindex arrary
!      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
!        CISCI_instance%index_amplitudeCore%values( i ) = i
!      enddo

      !! cleaning and finding limit
      if ( present ( auxm ) ) then

        !! get the new end position
        newm = CISCI_instance%buffer_amplitudeCoreSize/2 
        do i = 1, CISCI_instance%buffer_amplitudeCoreSize/2 + 1
          if ( CISCI_instance%buffer_amplitudeCore%values(i) == 0.0_8 ) then
            newm = i - 1
            exit
          endif
        enddo
        auxm = newm

        !! discard the last two quarters of tmp_ampltitude for next run, if not keep it fot PT2 corr
        CISCI_instance%buffer_amplitudeCore%values( CISCI_instance%buffer_amplitudeCoreSize/2 + 1 : CISCI_instance%buffer_amplitudeCoreSize ) = 0.0_8
        CISCI_instance%confAmplitudeCore%values( :, CISCI_instance%buffer_amplitudeCoreSize/2 + 1 : CISCI_instance%buffer_amplitudeCoreSize ) = 0.0_8
      endif

      !print *, "cleaned"
      !do i = 1, m 
      !  print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
      !enddo

      !m = CISCI_instance%buffer_amplitudeCoreSize/2 
   end subroutine CISCI_sortAmplitude



!  subroutine CISCI_appendAmplitudeb ( m, amplitude, indexConf )
!    implicit none
!    integer(8), intent(inout) :: m
!    real(8), intent(in) :: amplitude
!    real(8), intent(in) :: indexConf(:)
!    real(8) :: diffConf
!    integer :: i, ii, sp
!    logical :: addConf
!    integer :: nproc
!
!    nproc = omp_get_max_threads()
!
!    !if ( abs(amplitude) < 1E-16 ) return 
!    addConf = .true.
!
!    !! check if the configuration is already in the amplitudes set
!    do i = 1, m
!      diffConf = 0.0_8
!      do sp = 1, CIcore_instance%numberOfQuantumSpecies 
!        diffConf = diffConf + abs( indexConf(sp) - CISCI_instance%confAmplitudeCore%values(sp, CISCI_instance%index_amplitudeCore%values(i) ) )
!      enddo
!
!      if ( diffConf == 0 ) then
!        CISCI_instance%buffer_amplitudeCore%values(i) = CISCI_instance%buffer_amplitudeCore%values(i) + amplitude
!        addConf = .false.
!        exit 
!      end if
!    enddo
!
!
!    !! add the configuration at the end of the array
!    if ( addConf ) then
!      m = m + 1
!      CISCI_instance%buffer_amplitudeCore%values(m) = amplitude
!      CISCI_instance%confAmplitudeCore%values(:,m) = indexConf(:)
!      CISCI_instance%index_amplitudeCore%values(m) = m
!   !   print *, "m", m,  CISCI_instance%buffer_amplitudeCore%values(m), CISCI_instance%confAmplitudeCore%values(:,m),  CISCI_instance%index_amplitudeCore%values(m) 
!    endif
!
!    !! sort the array when it's full
!    if ( m == CISCI_instance%buffer_amplitudeCoreSize ) then
!
!      call MTSort ( CISCI_instance%buffer_amplitudeCore%values, &
!            CISCI_instance%index_amplitudeCore%values,  CISCI_instance%buffer_amplitudeCoreSize, "D", nproc )
!
!      !! discard the last two quarters of tmp_ampltitude
!      CISCI_instance%buffer_amplitudeCore%values( CISCI_instance%buffer_amplitudeCoreSize/2 + 1 : CISCI_instance%buffer_amplitudeCoreSize ) = 0.0_8
!
!
!      call CISCI_sortArrayByIndex( CISCI_instance%confAmplitudeCore, CISCI_instance%index_amplitudeCore, CISCI_instance%buffer_amplitudeCoreSize )
!
!
!      CISCI_instance%confAmplitudeCore%values( :, CISCI_instance%buffer_amplitudeCoreSize/2 + 1 : CISCI_instance%buffer_amplitudeCoreSize ) = 0.0_8
!
!      !! reset auxindex arrary
!      do i = 1, CISCI_instance%buffer_amplitudeCoreSize
!        CISCI_instance%index_amplitudeCore%values( i ) = i
!      enddo
!
!      !print *, "middle old sorting"
!
!      do i = 1, m 
!     !   print *, CISCI_instance%buffer_amplitudeCore%values(i), CISCI_instance%confAmplitudeCore%values(:,i),  CISCI_instance%index_amplitudeCore%values(i) 
!      enddo
!
!      !! the next configurations will be added in the last half
!      m = CISCI_instance%buffer_amplitudeCoreSize/2 
!
!    endif
!
!  end subroutine CISCI_appendAmplitudeb

  subroutine CISCI_buildPT2Row ( occA, occB, orbA, orbB, couplingS, factor0, CIenergy )
    implicit none
    integer(8) :: b,bb
    integer :: i, ii, jj, n, spi, spj
    integer :: numberOfSpecies
    real(8), intent(out) :: CIenergy
    real(8) :: auxCIenergy
    integer(1) :: coupling
    integer :: diffOrbi(4)
    integer :: diffOrbj(4)
    integer :: pi
    integer :: oib
    integer(1), intent(inout), allocatable :: couplingS(:)
    type (iVector), intent(in), allocatable :: occA(:)
    type (iVector), intent(inout), allocatable :: occB(:)
    type (ivector), intent(in), allocatable :: orbA(:)
    type (ivector), intent(inout), allocatable :: orbB(:)
    integer, intent(in) :: factor0
    integer :: factorA, factorB

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    CIenergy = 0.0_8

    !! calculate the sign factor for canonical order of the configuration A
    factorA = 1
    do spi = 1, numberOfSpecies 
      factorA = factorA * CISCI_canonicalOrderFactor( spi, orbA(spi), occA(spi) )
    enddo 

    do bb = 1, CISCI_instance%targetSpaceSize

      b = CISCI_instance%index_amplitudeCore%values(bb)
      ! getting configuration B
      do spi = 1, numberOfSpecies 

        oib = 0 

        !! build the orbital from the index using the bit mapping
        !!call CISCI_decimalToBinary ( CISCI_instance%saved_confTarget%values(spi,b), orbB(spi)%values )
        call CISCI_decimalToBinary ( CISCI_instance%confAmplitudeCore%values(spi,b), orbB(spi)%values )

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
                          - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) + 1
      end do

      coupling = product(couplingS)

      select case (coupling)

      case(2)
        do i = 1, numberOfSpecies
          if ( couplingS(i) == 2 ) spi = i
        end do

        diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
        CIenergy = CIenergy + CISCI_calculateEnergyOneNew( spi, occA, occB, diffOrbi(1), diffOrbi(3)  ) * factorA * factor0 * CISCI_instance%auxcoefficientTarget%values(bb) 

      case(3)
        do i = 1, numberOfSpecies
          if ( couplingS(i) == 3 ) spi = i
        end do

        diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
        CIenergy = CIenergy + CISCI_calculateEnergyTwoSameNew( spi, occA, occB, diffOrbi(1), diffOrbi(2), diffOrbi(3), diffOrbi(4)  ) * factorA * factor0 * CISCI_instance%auxcoefficientTarget%values(bb) 

      case(4)
        do i = 1, numberOfSpecies
          if ( couplingS(i) == 2 ) then 
            spi = i
            exit
          end if
        end do
        do i = ii+1, numberOfSpecies
          if ( couplingS(i) == 2 ) spj = i
        end do

        diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi), orbB(spi), occA(spi), occB(spi), factorA )
        diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj), orbB(spj), occA(spj), occB(spj), factorB )
        auxCIenergy = CISCI_calculateEnergyTwoDiffNew( spi, spj, diffOrbi(1), diffOrbj(1), diffOrbi(3), diffOrbj(3)  )
        CIenergy = CIenergy + auxCIenergy * factorA * factorB * factor0 * CISCI_instance%auxcoefficientTarget%values(bb) 

      end select

    end do !b

    CIenergy = CIenergy * CIenergy

  end subroutine CISCI_buildPT2Row

  function CISCI_searchConfInTargetSpace( indexConf ) result (found) 
    implicit none
    real(8), intent(in) :: indexConf(:)
    logical :: found
    real(8) :: diffConf
    integer :: i, sp
    integer(8) :: b,bb

    found = .false.
    !! check if the configuration is already in the amplitudes set
    do bb = 1, CISCI_instance%targetSpaceSize
      b = CISCI_instance%index_amplitudeCore%values(bb)
      diffConf = 0.0_8
      do sp = 1, CIcore_instance%numberOfQuantumSpecies 
        diffConf = diffConf + ( indexConf(sp) - CISCI_instance%confAmplitudeCore%values(sp, b ) )
      enddo

      if ( diffConf == 0 ) then
        found = .true.
        exit 
      end if
    enddo

  endfunction CISCI_searchConfInTargetSpace 

  !! Reordering an array based on sorted indices in another array
  !! Algorithm taken from Zecong Hu
  !! https://stackoverflow.com/questions/60917343/

  subroutine CISCI_sortArrayByIndex( indexConf, auxindex_array, vectorSize  )
    implicit none
    type(Matrix), intent(inout) :: indexConf
    type(ivector8), intent(inout) :: auxindex_array
    integer(8) :: vectorSize
    real(8), allocatable :: temp(:)
    integer :: numberOfSpecies 
    integer :: i, x, y

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate ( temp ( numberOfSpecies ) )

    temp(:) = 0.0_8

    do i = 1, vectorSize
      if ( auxindex_array%values(i) == -1 ) cycle 
      temp(:) = indexConf%values(:,i) 
      x = i
      y = auxindex_array%values(i) 
      do while ( y /= i )
        auxindex_array%values(x) = -1
        indexConf%values(:,x) = indexConf%values(:,y)
        x = y
        y = auxindex_array%values(x)
      enddo 
      indexConf%values(:,x) = temp(:)
      auxindex_array%values(x) = -1
    enddo

    deallocate ( temp ) 

  endsubroutine CISCI_sortArrayByIndex


  subroutine CISCI_sortVectorByIndex( targetVector, auxindex_array, vectorSize  )
    implicit none
    type(vector8), intent(inout) :: targetVector
    type(ivector8), intent(inout) :: auxindex_array
    integer(8) :: vectorSize
    real(8) :: temp
    integer :: numberOfSpecies 
    integer :: i, x, y

    temp = 0.0_8

    do i = 1, vectorSize
      if ( auxindex_array%values(i) == -1 ) cycle 
      temp = targetVector%values(i) 
      x = i
      y = auxindex_array%values(i) 
      do while ( y /= i )
        auxindex_array%values(x) = -1
        targetVector%values(x) = targetVector%values(y)
        x = y
        y = auxindex_array%values(x)
      enddo 
      targetVector%values(x) = temp
      auxindex_array%values(x) = -1
    enddo

  endsubroutine CISCI_sortVectorByIndex

end module CISCI_
