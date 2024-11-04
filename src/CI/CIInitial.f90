module CIInitial_
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
  use CIcore_

contains


  subroutine CIInitial_buildInitialCIMatrix2()
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
    if ( CIcore_instance%numberOfConfigurations <  CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX ) then
      CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX = CIcore_instance%numberOfConfigurations !! assign to an internal variable
    end if

    call Vector_constructorInteger8 ( CIcore_instance%auxIndexCIMatrix, &
                              CIcore_instance%numberOfConfigurations, 0_8 ) !hmm

    do a = 1, CIcore_instance%numberOfConfigurations
      CIcore_instance%auxIndexCIMatrix%values(a)= a
    end do

   !! save the unsorted diagonal Matrix
    call Vector_constructor8 ( CIcore_instance%diagonalHamiltonianMatrix, &
                              CIcore_instance%numberOfConfigurations, 0.0_8 ) 


    CIcore_instance%diagonalHamiltonianMatrix%values = CIcore_instance%diagonalHamiltonianMatrix2%values

   !! To get only the lowest 300 values.
   call Vector_reverseSortElements8( CIcore_instance%diagonalHamiltonianMatrix2, &
          CIcore_instance%auxIndexCIMatrix, int(initialCIMatrixSize,8))

   call Matrix_constructor ( CIcore_instance%initialHamiltonianMatrix, int(initialCIMatrixSize,8) , &
                               int(initialCIMatrixSize,8) , 0.0_8 ) 

    !! get the configurations for the initial hamiltonian matrix
    call CIInitial_getInitialIndexes()

    call CIInitial_calculateInitialCIMatrix()

    !! diagonalize the initial matrix
    call Vector_constructor8 ( CIcore_instance%initialEigenValues, int(CONTROL_instance%NUMBER_OF_CI_STATES,8),  0.0_8)

    call Matrix_constructor (CIcore_instance%initialEigenVectors, &
          int(initialCIMatrixSize,8), &
          int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

    call Matrix_eigen_select ( CIcore_instance%initialHamiltonianMatrix, &
           CIcore_instance%initialEigenValues, &
           1, int(CONTROL_instance%NUMBER_OF_CI_STATES,4), &  
           eigenVectors = CIcore_instance%initialEigenVectors, &
           flags = int(SYMMETRIC,4))
    
    write(*,*) "Initial eigenValues"
    do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
      write (*,*)  i, CIcore_instance%initialEigenValues%values(i)
    end do

    call Vector_destructor8 ( CIcore_instance%diagonalHamiltonianMatrix2 )

!$    timeB = omp_get_wtime()
!$    write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time for Solving Initial CI : ", timeB - timeA ," (s)"

  end subroutine CIInitial_buildInitialCIMatrix2

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  !! Map the indexes of initial CI matrix to the complete matrix.
  subroutine CIInitial_getInitialIndexes()
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

    call Matrix_constructorInteger ( CIcore_instance%auxConfigurations, int( numberOfSpecies,8), &
          int(CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX,8), 0 )

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
      auxnumberOfSpecies = CIInitial_getIndexesRecursion( s, numberOfSpecies, indexConf, c, cilevel )
    end do

    deallocate ( indexConf )
    deallocate ( cilevel )

!$  timeB = omp_get_wtime()

!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for getting initial indexes : ", timeB - timeA ," (s)"

  end subroutine CIInitial_getInitialIndexes


recursive  function CIInitial_getIndexesRecursion(s, numberOfSpecies, indexConf, c, cilevel) result (os)
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
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a
        os = CIInitial_getIndexesRecursion( is, numberOfSpecies, indexConf, c, cilevel)
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        do u = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
          if ( c ==  CIcore_instance%auxIndexCIMatrix%values(u) ) then
            do ss = 1, numberOfSpecies
              CIcore_instance%auxConfigurations%values(ss,u) = indexConf(ss) 
            end do
          end if
        end do
      end do
    end if

  end function CIInitial_getIndexesRecursion

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CIInitial_calculateInitialCIMatrix()
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
      aa = CIcore_instance%auxIndexCIMatrix%values(a)
      do b = a, initialCIMatrixSize 
        bb = CIcore_instance%auxIndexCIMatrix%values(b)
        coupling = 0

        indexConfA = 0
        indexConfB = 0

        do i = 1, numberOfSpecies

          allocate (orbitalsA ( CIcore_instance%numberOfOrbitals%values(i) ))
          allocate (orbitalsB ( CIcore_instance%numberOfOrbitals%values(i) ))
          orbitalsA = 0
          orbitalsB = 0

          indexConfA(i) =  CIcore_instance%auxConfigurations%values(i,a)
          indexConfB(i) =  CIcore_instance%auxConfigurations%values(i,b)

          do u = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
            orbitalsA( CIcore_instance%strings(i)%values(u,indexConfA(i) ) ) = 1
          end do
          do v = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
            orbitalsB( CIcore_instance%strings(i)%values(v,indexConfB(i) ) ) = 1
          end do
          coupling = coupling + &
            CIcore_instance%numberOfOccupiedOrbitals%values(i) - sum ( orbitalsA * orbitalsB ) 

          deallocate (orbitalsA )
          deallocate (orbitalsB )

        end do
        if ( coupling  == 0 ) then
          CIcore_instance%initialHamiltonianMatrix%values(a,b) = &
            CIcore_instance%diagonalHamiltonianMatrix2%values(a) 

        else if (  coupling == 1 ) then

          CIcore_instance%initialHamiltonianMatrix%values(a,b) = &
            CIInitial_calculateEnergyOne ( 1, indexConfA, indexConfB )

        else if ( coupling  == 2 ) then

          CIcore_instance%initialHamiltonianMatrix%values(a,b) = &
            CIInitial_calculateEnergyTwo ( 1, indexConfA, indexConfB )

        end if
 

      end do


    end do

    deallocate ( indexConfB )
    deallocate ( indexConfA )

!$  timeB1 = omp_get_wtime()
    !! symmetrize
    do a = 1, initialCIMatrixSize 
      do b = a, initialCIMatrixSize 

         CIcore_instance%initialHamiltonianMatrix%values(b,a) = &
            CIcore_instance%initialHamiltonianMatrix%values(a,b) 
      end do
    end do

    !!open(unit=318, file="cimatrix.dat", action = "write", form="formatted")
    !!do a = 1, initialCIMatrixSize 
    !!  do b = 1, initialCIMatrixSize 
    !!     write (318,*) a,b, CIcore_instance%initialHamiltonianMatrix%values(a,b)
    !!  end do
    !!     write (318,*) " "
    !!end do
    !!close(318)
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Calculating initial CI matrix : ", timeB1 - timeA1 ," (s)"

  end subroutine CIInitial_calculateInitialCIMatrix 

  function CIInitial_calculateEnergyOne( n, thisA, thisB ) result (auxCIenergy)
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

      CIcore_instance%auxstring(n,i)%values(:) = CIcore_instance%strings(i)%values(:,a)
    end do

    !! set at maximum coincidence

    do s = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(s)
      b = thisB(s)

      do i = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s) !b
        do j = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s) !a
          if ( CIcore_instance%auxstring(n,s)%values(j) == &
             CIcore_instance%strings(s)%values(i,b) ) then

            auxOcc = CIcore_instance%auxstring(n,s)%values(i) 
            CIcore_instance%auxstring(n,s)%values(i) = CIcore_instance%strings(s)%values(i,b)
            CIcore_instance%auxstring(n,s)%values(j) = auxOcc
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

      do kk = 1, CIcore_instance%occupationNumber( i) !! 1 is from a and 2 from b

        if ( CIcore_instance%auxstring(n,i)%values(kk) .ne. &
                 CIcore_instance%strings(i)%values(kk,b) ) then
          diffOrb(1) = CIcore_instance%auxstring(n,i)%values(kk)
          diffOrb(2) = CIcore_instance%strings(i)%values(kk,b)
          exit                   
        end if

      end do
      if (  diffOrb(2) > 0 ) then 

        !One particle terms
        auxCIenergy= auxCIenergy +  CIcore_instance%twoCenterIntegrals(i)%values( &
                           diffOrb(1), diffOrb(2) )

        auxIndex1= CIcore_instance%twoIndexArray(i)%values( & 
                         diffOrb(1), diffOrb(2))

        do ll = 1, CIcore_instance%occupationNumber( i ) !! 1 is from a and 2 from b

          if ( CIcore_instance%auxstring(n,i)%values(ll) .eq. &
                 CIcore_instance%strings(i)%values(ll,b) ) then

            l = CIcore_instance%auxstring(n,i)%values(ll) !! or b

            auxIndex2 = CIcore_instance%twoIndexArray(i)%values( l,l) 

            auxIndex = CIcore_instance%fourIndexArray(i)%values( auxIndex1, auxIndex2 )

            auxCIenergy = auxCIenergy + &
                        CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)


            auxIndex = CIcore_instance%fourIndexArray(i)%values( &
                                CIcore_instance%twoIndexArray(i)%values(diffOrb(1),l), &
                                CIcore_instance%twoIndexArray(i)%values(l,diffOrb(2)) ) 

            auxCIenergy = auxCIenergy + &
                    MolecularSystem_instance%species(i)%kappa*CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

          end if
        end do
        if (MolecularSystem_instance%numberOfQuantumSpecies .gt. 1 ) then !.and. spin(1) .eq. spin(2) ) then
          do j=1, MolecularSystem_instance%numberOfQuantumSpecies

            if (i .ne. j) then

              auxnumberOfOtherSpecieSpatialOrbitals = CIcore_instance%numberOfSpatialOrbitals2%values(j) 

              do ll=1,  CIcore_instance%occupationNumber( j ) !! 1 is from a and 2 from b
                l = CIcore_instance%auxstring(n,j)%values(ll) !! or b?

                auxIndex2 = CIcore_instance%twoIndexArray(j)%values( l,l) 
                auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

                auxCIenergy = auxCIenergy + &
                      CIcore_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1) 
              end do
            end if
          end do
        end if
      end if
    end do

    auxCIenergy= auxCIenergy * factor


  end function CIInitial_calculateEnergyOne


  function CIInitial_calculateEnergyTwo( n, thisA, thisB ) result (auxCIenergy)
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
      CIcore_instance%auxstring(n,i)%values(:) = CIcore_instance%strings(i)%values(:,a)
    end do

    !! set at maximum coincidence

    do s = 1, MolecularSystem_instance%numberOfQuantumSpecies
      a = thisA(s)
      b = thisB(s)

      do i = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s) !b
        do j = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s) !a
          if ( CIcore_instance%auxstring(n,s)%values(j) == &
                 CIcore_instance%strings(s)%values(i,b) ) then

            auxOcc = CIcore_instance%auxstring(n,s)%values(i) 
            CIcore_instance%auxstring(n,s)%values(i) = CIcore_instance%strings(s)%values(i,b)
            CIcore_instance%auxstring(n,s)%values(j) = auxOcc
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
      do k = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)

        if ( CIcore_instance%auxstring(n,i)%values(k) .ne. &
                 CIcore_instance%strings(i)%values(k,b) ) then
          diffOrb(z) = CIcore_instance%auxstring(n,i)%values(k) 
          diffOrb(z+2) = CIcore_instance%strings(i)%values(k,b)  
          z = z + 1
          cycle
        end if 
      end do 
      if (  diffOrb(2) > 0 ) then

        !Coulomb
        auxIndex = CIcore_instance%fourIndexArray(i)%values( &
                          CIcore_instance%twoIndexArray(i)%values(&
                             diffOrb(1),diffOrb(3)),&
                          CIcore_instance%twoIndexArray(i)%values(&
                             diffOrb(2),diffOrb(4)) )

         auxCIenergy = auxCIenergy + &
                  CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

         auxIndex = CIcore_instance%fourIndexArray(i)%values( &
                          CIcore_instance%twoIndexArray(i)%values(&
                             diffOrb(1),diffOrb(4)),&
                          CIcore_instance%twoIndexArray(i)%values(&
                             diffOrb(2),diffOrb(3)) )

         auxCIenergy = auxCIenergy + &
               MolecularSystem_instance%species(i)%kappa*CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

      end if
      !! different species
      do j = i + 1, MolecularSystem_instance%numberOfQuantumSpecies
        auxnumberOfOtherSpecieSpatialOrbitals = CIcore_instance%numberOfSpatialOrbitals2%values(j) 
        otherdiffOrb = 0
        a = thisA(j)
        b = thisB(j)

        do k = 1, CIcore_instance%numberOfOccupiedOrbitals%values(j)
          if ( CIcore_instance%auxstring(n,j)%values(k) .ne. &
                CIcore_instance%strings(j)%values(k,b) ) then
            otherdiffOrb(1) = CIcore_instance%auxstring(n,j)%values(k)
            otherdiffOrb(3) = CIcore_instance%strings(j)%values(k,b)
            exit 
          end if 

        end do 

        if ( diffOrb(3) .gt. 0 .and. otherdiffOrb(3) .gt. 0 ) then
          auxIndex1 = CIcore_instance%twoIndexArray(i)%values(&
                                   diffOrb(1),diffOrb(3) )
          auxIndex2 = CIcore_instance%twoIndexArray(j)%values(&
                                   otherdiffOrb(1),otherdiffOrb(3) )
          auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

          auxCIenergy = auxCIenergy + &
                        CIcore_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)
  
        end if
      end do
    end do

    auxCIenergy= auxCIenergy * factor

  end function CIInitial_calculateEnergyTwo

end module CIInitial_
