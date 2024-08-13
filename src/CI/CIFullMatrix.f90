 module CIFullMatrix_
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

!>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CIFullMatrix_buildHamiltonianMatrix()
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

    numberOfConfigurations = CIcore_instance%numberOfConfigurations 

    allocate ( CIcore_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )
    ciLevel = 0
    CIcore_instance%allIndexConf = 0
    indexConf = 0

    !! gather all configurations
    s = 0
    c = 0
    ciLevel = 0

    do ci = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = CIcore_gatherConfRecursion( s, numberOfSpecies, indexConf,  c, cilevel )
    end do
    !stop

    deallocate ( indexConf )
    deallocate ( ciLevel )

    !! allocate the hamiltonian matrix
    call Matrix_constructor ( CIcore_instance%hamiltonianMatrix, & 
           int(CIcore_instance%numberOfConfigurations,8), &
           int(CIcore_instance%numberOfConfigurations,8), 0.0_8)


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
!$omp& shared(CIcore_instance, HartreeFock_instance)
    n = omp_get_thread_num() + 1
!$omp do schedule (dynamic) 
    do a = 1, numberOfConfigurations
      indexConfA(:,n) = CIcore_instance%allIndexConf(:,a) 
      do b = a, numberOfConfigurations

        indexConfB(:,n) = CIcore_instance%allIndexConf(:,b) 

        do i = 1, numberOfSpecies
          if ( pindexConf(i,n) /= indexConfB(i,n) ) then
            allocate (stringAinB (CIcore_instance%numberOfOccupiedOrbitals%values(i) ))
            stringAinB = 0
            do p = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
              stringAinB(p) = CIcore_instance%orbitals(i)%values( &
                                CIcore_instance%strings(i)%values(p,indexConfA(i,n) ), indexConfB(i,n) ) 
            end do
            couplingSpecies(i,n) = CIcore_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB )
            deallocate (stringAinB )
          end if
        end do
        coupling = sum(couplingSpecies(:,n))

        if ( coupling  == 0 ) then
          CIcore_instance%hamiltonianMatrix%values(a,b) = &
            CIcore_instance%diagonalHamiltonianMatrix2%values(a) 
  
        else if (  coupling == 1 ) then
  
          CIcore_instance%hamiltonianMatrix%values(a,b) = &
            CIFullMatrix_calculateEnergyOne ( n, indexConfA(:,n), indexConfB(:,n) )
  
        else if ( coupling  == 2 ) then
  
          CIcore_instance%hamiltonianMatrix%values(a,b) = &
            CIFullMatrix_calculateEnergyTwo ( n, indexConfA(:,n), indexConfB(:,n) )
  
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
         CIcore_instance%hamiltonianMatrix%values(b,a) = &
            CIcore_instance%hamiltonianMatrix%values(a,b) 
      end do
    end do

    deallocate ( CIcore_instance%allIndexConf )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"

  end subroutine CIFullMatrix_buildHamiltonianMatrix

  function CIFullMatrix_calculateEnergyOne( n, thisA, thisB ) result (auxCIenergy)
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


  end function CIFullMatrix_calculateEnergyOne

  function CIFullMatrix_calculateEnergyTwo( n, thisA, thisB ) result (auxCIenergy)
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

  end function CIFullMatrix_calculateEnergyTwo


 end module CIFullMatrix_
