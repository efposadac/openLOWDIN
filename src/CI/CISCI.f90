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
  use omp_lib
  implicit none

  type, public :: CISCI
    type (Vector8) :: amplitudeCore
    type (Vector8) :: amplitudeCore2
    type (Vector8) :: coefficientCore
    type (Matrix) :: coefficientTarget
    type (Vector8) :: auxcoefficientTarget
    type (Vector8) :: diagonalTarget
    type (Vector8) :: eigenValues
    integer :: coreSpaceSize
    integer :: targetSpaceSize
    real(8) :: PT2energy
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

    CISCI_instance%coreSpaceSize = CONTROL_instance%CI_SCI_CORE_SPACE
    CISCI_instance%targetSpaceSize = CONTROL_instance%CI_SCI_TARGET_SPACE

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
    !call CISCI_getInitialIndexes( CIcore_instance%coreConfigurations, CIcore_instance%coreConfigurationsLevel, CISCI_instance%coreSpaceSize )
    !call CISCI_getInitialIndexes( CIcore_instance%targetConfigurations, CIcore_instance%targetConfigurationsLevel, CISCI_instance%targetSpaceSize )

    call CISCI_getInitialIndexes2( CIcore_instance%coreConfigurations, CIcore_instance%coreConfigurationsLevel, CISCI_instance%coreSpaceSize )
    call CISCI_getInitialIndexes2( CIcore_instance%targetConfigurations, CIcore_instance%targetConfigurationsLevel, CISCI_instance%targetSpaceSize )



    !! arrays for CISCI
    call Vector_constructor8 ( CISCI_instance%amplitudeCore, CIcore_instance%numberOfConfigurations,  0.0_8)
    call Vector_constructor8 ( CISCI_instance%coefficientCore, int(CISCI_instance%coreSpaceSize,8),  0.0_8)
    call Matrix_constructor ( CISCI_instance%coefficientTarget, int(CISCI_instance%targetSpaceSize,8), &
                               int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)
    call Vector_constructor8 ( CISCI_instance%auxcoefficientTarget, int(CISCI_instance%targetSpaceSize,8), 0.0_8) !! meh... just to avoid changing everything from matrix to vector format
    call Vector_constructor8 ( CISCI_instance%diagonalTarget, int(CISCI_instance%targetSpaceSize,8),  0.0_8) !! Jadamilu requires to store diagonal vector (in target space)
    call Vector_constructor8 ( CISCI_instance%eigenValues, 15_8, 0.0_8) !! store the eigenvalues per macro iterations

  end subroutine CISCI_constructor
  
  !! main part
  subroutine CISCI_run()
    implicit none
    integer(8) :: i, j, ii, jj
    integer :: k ! macro SCI iteration
    real(8) :: timeA(15), timeB(15)
    real(8) :: timeAA, timeBB
    real(8) :: timeAS, timeBS
    type(Vector8) :: eigenValuesTarget
    real(8) :: currentEnergy 

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
      call CISCI_core_amplitudes ( CISCI_instance%amplitudeCore%values, CIcore_instance%numberOfConfigurations, &
                                  CISCI_instance%coefficientCore%values, CISCI_instance%coreSpaceSize, currentEnergy )

      !! setting the HF again... because the HF amplitude diverges 
      if ( k == 2 ) CISCI_instance%amplitudeCore%values(1) = 1.0_8

      !! resetting the original index changes
      do i = 1, CIcore_instance%numberOfConfigurations
        CIcore_instance%auxIndexCIMatrix%values(i)= i
      end do

      !! getting the target absolute largest coefficients
      call Vector_sortElementsAbsolute8( CISCI_instance%amplitudeCore, &
            CIcore_instance%auxIndexCIMatrix, int( CISCI_instance%targetSpaceSize ,8), CONTROL_instance%CI_MATVEC_TOLERANCE )

     ! do i = 1,  CISCI_instance%targetSpaceSize
     !   print *, i, CISCI_instance%amplitudeCore%values(i), CIcore_instance%auxIndexCIMatrix%values(i)
     ! end do
     ! stop 
      !! recover the configurations for the hamiltonian matrix in the target space
      call CISCI_getInitialIndexes2( CIcore_instance%targetConfigurations, CIcore_instance%targetConfigurationsLevel, CISCI_instance%targetSpaceSize )

      !! storing only the largest diagonal elements (for jadamilu)
      do i = 1,  CISCI_instance%targetSpaceSize
        ii = CIcore_instance%auxIndexCIMatrix%values(i)
        !! storing only the largest diagonal elements (for jadamilu)
        CISCI_instance%diagonalTarget%values(i) = CIcore_instance%diagonalHamiltonianMatrix%values(ii)
      enddo

      !! using the amplitued as the initial coeff guess, after that, use the previous diganolized eigenvectors in target space
      if ( k == 2 ) then 
        do i = 1,  CISCI_instance%targetSpaceSize
          ii = CIcore_instance%auxIndexCIMatrix%values(i)
          CISCI_instance%coefficientTarget%values(i,1) = CISCI_instance%amplitudeCore%values(ii)
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
      call Vector_sortElementsAbsolute8( CISCI_instance%auxcoefficientTarget, &
           CIcore_instance%auxIndexCIMatrix, int( CISCI_instance%coreSpaceSize ,8), CONTROL_instance%CI_MATVEC_TOLERANCE )

      !! recover the configurations for the hamiltonian matrix in the core space
      call CISCI_getInitialIndexes2( CIcore_instance%coreConfigurations, CIcore_instance%coreConfigurationsLevel, CISCI_instance%coreSpaceSize )
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
    integer :: auxMatrixSize
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
      if ( u2 > auxMatrixSiez ) u2 = auxMatrixSize
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

!      print *, CIcore_instance%numberOfStrings(is)%values(i),  CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        do u = 1, auxMatrixSize 
          if ( c ==  CIcore_instance%auxIndexCIMatrix%values(u) ) then
!            print *, c, indexConf(:)
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
        if ( j2 > n ) j2 = n
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

!      print *, c, indexConf
      !print *, indexConf(:,1)
    end do


  end function CISCI_getIndexes

!  function CISCI_buildMatrixRecursion2( c, indexConf, cilevel) result (os)
!    implicit none
!
!    integer(8) :: a,aa, x
!    integer :: i, j, ci
!    integer :: numberOfSpecies
!    integer :: os,is,ss,ssize
!    integer(8) :: indexConf(:)
!    integer :: cilevel(:)
!    integer(8) :: c, totalsize, auxtotalsize
!    integer, allocatable :: counter(:)
!
!    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
!    
!    allocate (counter(numberOfSpecies))
!    counter = 0 
!
!    totalsize = 1
!    do i = 1 , numberOfSpecies
!      totalsize = totalsize * CIcore_instance%numberOfStrings(i)%values(cilevel(i) + 1)
!    end do
!
!    do i = 1 , numberOfSpecies 
!      ci = cilevel(i) + 1 
!      ssize = CIcore_instance%numberOfStrings2(i)%values(ci)
!      indexConf(i) = ssize  + 1
!    end do
!
!    indexConf(numberOfSpecies) = indexConf(numberOfSpecies) -1
!
!    do x = 1, totalsize
!      c = c + 1
!      indexConf(numberOfSpecies) = indexConf(numberOfSpecies) + 1
!
!      do i = numberOfSpecies, 1 + 1, -1 
!        auxtotalsize = 1
!        do j = i, numberOfSpecies
!          auxtotalsize = auxtotalsize * CIcore_instance%numberOfStrings(j)%values(cilevel(j) + 1)
!        end do
!        if (counter(i) == auxtotalsize) then
!          do j = i, numberOfSpecies
!            ci = cilevel(j) + 1 
!            ssize = CIcore_instance%numberOfStrings2(j)%values(ci)
!            indexConf(j) = ssize + 1
!          end do
!          counter(i) = 0
!          indexConf(i-1) = indexConf(i-1) + 1
!        end if
!        counter(i) = counter(i) + 1 
!
!      end do
!        print *, c, indexConf
!      !print *, indexConf(:,1)
!    end do
!
!    deallocate (counter)
!
!  end function CISCI_buildMatrixRecursion2




  subroutine CISCI_core_amplitudes ( amplitudeCore, numberOfConfigurations, coefficientCore, SCICoreSpaceSize, oldEnergy )

    implicit none
  
    integer(4) SCICoreSpaceSize
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
    integer(8) :: a,b,c, aa
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
!    shift = 0.0_8
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
    !stop
    deallocate (orbitalsA )
    deallocate (orbitalsB )
    deallocate ( indexConfB )
    deallocate ( indexConfA )
    deallocate ( couplingS )

!$    write(*,"(A,I2,A,E10.3,A2,I12,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero, nonzerow
    return

  end subroutine CISCI_matvec

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
    !diffOrbA = 0
    !diffOrbB = 0
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

  subroutine CISCI_PT2 ( coefficients, auxenergyCorrection, SCITargetSpaceSize, numberOfConfigurations, refEnergy, energyCorrection )
  
  !*******************************************************************************
  !! AV computes w <- A * V where A is a discretized Laplacian.
  !  Parameters:
  !    Input, integer NX, the length of the vectors.
  !    Input, real V(NX), the vector to be operated on by A.
  !    Output, real W(NX), the result of A*V.
  !
    implicit none
    integer :: SCITargetSpaceSize
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

    !$omp parallel &
    !$omp& private( n, a, aa, cilevel, indexConfA ),&
    !$omp& shared( auxenergyCorrection, coefficients, refEnergy   ) 
    !$omp do schedule (static) 
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
        auxenergyCorrection(bb) = 0.0_8
      enddo

    end do
    !$omp end do nowait
    !$omp end parallel

    do b = 1, numberOfConfigurations
      energyCorrection = energyCorrection + auxenergyCorrection(b) **2 / ( refEnergy - CIcore_instance%diagonalHamiltonianMatrix%values(b)  )
    enddo


    CIcore_instance%pindexConf = 0
!$  timeB = omp_get_wtime()
    deallocate ( cilevel )
    deallocate ( indexConfA )
!$    write(*,"(A,E10.3)") "Time for CI-PT2 correction: ", timeB -timeA

  end subroutine CISCI_PT2

end module CISCI_
