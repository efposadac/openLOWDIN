module CIDiag_
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
  subroutine CIDiag_buildDiagonal()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: ci
    integer :: nproc, n, nn
    integer :: i, j, ii, jj
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:,:)
    integer(8), allocatable :: cc(:)
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

    ciLevel = 0
    auxciLevel = 0

    !!auxnumberOfSpecies = CIcore_numberOfConfigurationsRecursion2(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    numberOfConfigurations = 0
    ciLevel = 0

    !! call recursion to get the number of configurations...
    do ci = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = CIDiag_numberOfConfigurationsRecursion(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    end do

    call Vector_constructor8 ( CIcore_instance%diagonalHamiltonianMatrix2, &
                              numberOfConfigurations, 0.0_8 ) 

    CIcore_instance%numberOfConfigurations = numberOfConfigurations 

    write (*,*) "Number Of Configurations: ", numberOfConfigurations

    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()
    allocate ( indexConf ( numberOfSpecies, nproc ) )
    allocate ( cc ( nproc ) )
    indexConf = 0

    !! calculate the diagonal 
    s = 0
    c = 0
    ciLevel = 0
    n = 1

    do ci = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0

      auxnumberOfSpecies = CIDiag_buildDiagonalRecursion( s, numberOfSpecies, nproc, indexConf,  c, cc, n, cilevel )
    end do

    !! computing the final batch < ncore
    if  ( n > 1 ) then
      do nn = 1, n-1
        CIcore_instance%diagonalHamiltonianMatrix2%values(cc(nn)) = CIDiag_calculateEnergyZero ( indexConf(:,nn) )
      end do
    end if

    !stop

    deallocate ( cc )
    deallocate ( indexConf )
    deallocate ( ciLevel )
    deallocate ( auxciLevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Building diagonal of CI matrix : ", timeB - timeA ," (s)"

    write (*,*) "Reference energy, H_0: ",  CIcore_instance%diagonalHamiltonianMatrix2%values(1)

  end subroutine CIDiag_buildDiagonal

recursive  function CIDiag_numberOfConfigurationsRecursion(s, numberOfSpecies, c, cilevel) result (os)
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
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        os = CIDiag_numberOfConfigurationsRecursion( is, numberOfSpecies, c, cilevel )
      end do
    else 
      os = is

      i = cilevel(is) + 1
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
      end do
    end if

  end function CIDiag_numberOfConfigurationsRecursion


recursive  function CIDiag_buildDiagonalRecursion(s, numberOfSpecies, nproc, indexConf, c, cc, n, cilevel ) result (os)
    implicit none

    integer(8), intent(inout) :: c
    integer(8), intent(inout) :: cc(:)
    integer(8), intent(inout) :: indexConf(:,:)
    integer, intent(in) :: nproc
    integer, intent(in) :: cilevel(:)
    integer(8) :: a
    integer :: i, j, ii, jj 
    integer :: s, numberOfSpecies
    integer :: os,is
    integer :: n, nn
    integer :: ssize

    is = s + 1
    if ( is < numberOfSpecies ) then
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is,n:) = ssize + a

        !dd(is) =(a + CIcore_instance%ciOrderSize1(u,is))* CIcore_instance%ciOrderSize2(u,is) 
        os = CIDiag_buildDiagonalRecursion( is, numberOfSpecies, nproc, indexConf, c, cc, n, cilevel )
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        cc(n) = c
        indexConf(is,n:) = ssize + a
        !dd(is) =(a + CIcore_instance%ciOrderSize1(u,is))* CIcore_instance%ciOrderSize2(u,is) 
        !d = sum(dd)

        !! saving to a stock to send in parallel
        if ( n == nproc ) then
          !$omp parallel &
          !$omp& private(nn),&
          !$omp& shared( indexConf, nproc) 
          !$omp do schedule (static) 
          do nn = 1, nproc
            CIcore_instance%diagonalHamiltonianMatrix2%values(cc(nn)) = CIDiag_calculateEnergyZero ( indexConf(:,nn) )
          end do
          !$omp end do nowait
          !$omp end parallel

          !!reset the stock
          n = 0 
          do nn = 1, nproc
            indexConf(:,nn) = indexConf(:,nproc) 
            !cilevel(:,nn) = cilevel(:,nproc) 
          end do

        end if

        n = n + 1

      end do
    end if

  end function CIDiag_buildDiagonalRecursion

  function CIDiag_calculateEnergyZero( this ) result (auxCIenergy)
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
      do kk=1, CIcore_instance%occupationNumber( i )  !! 1 is from a and 2 from b

        k = CIcore_instance%strings(i)%values(kk,a)

        !One particle terms
        auxCIenergy = auxCIenergy + &
                    CIcore_instance%twoCenterIntegrals(i)%values( k, k )

        !Two particles, same specie
        auxIndex1 = CIcore_instance%twoIndexArray(i)%values(k,k)

        do ll = kk + 1, CIcore_instance%occupationNumber( i )  !! 1 is from a and 2 from b

          l = CIcore_instance%strings(i)%values(ll,a)
          auxIndex2 = CIcore_instance%twoIndexArray(i)%values(l,l)
          auxIndex = CIcore_instance%fourIndexArray(i)%values(auxIndex1,auxIndex2) 

          !Coulomb
          auxCIenergy = auxCIenergy + &
              CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)

          !Exchange, depends on spin

          auxIndex = CIcore_instance%fourIndexArray(i)%values( &
                        CIcore_instance%twoIndexArray(i)%values(k,l), &
                        CIcore_instance%twoIndexArray(i)%values(l,k) )

          auxCIenergy = auxCIenergy + &
                  MolecularSystem_instance%species(i)%kappa*CIcore_instance%fourCenterIntegrals(i,i)%values(auxIndex, 1)
        end do

        !!Two particles, different species
        do j = i + 1, MolecularSystem_instance%numberOfQuantumSpecies
          b = this(j)
          auxnumberOfOtherSpecieSpatialOrbitals = CIcore_instance%numberOfSpatialOrbitals2%values(j) 

          do ll = 1, CIcore_instance%occupationNumber( j ) !! 1 is from a and 2 from b
            l = CIcore_instance%strings(j)%values(ll,b)

            auxIndex2= CIcore_instance%twoIndexArray(j)%values(l,l)
            auxIndex = auxnumberOfOtherSpecieSpatialOrbitals * (auxIndex1 - 1 ) + auxIndex2

            auxCIenergy = auxCIenergy + &!couplingEnergy
            CIcore_instance%fourCenterIntegrals(i,j)%values(auxIndex, 1)

          end do

        end do

      end do
    end do

    auxCIenergy= auxCIenergy + HartreeFock_instance%puntualInteractionEnergy

  end function CIDiag_calculateEnergyZero


end module CIDiag_
