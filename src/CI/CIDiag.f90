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
  subroutine CIcore_buildDiagonal()
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

    !!auxnumberOfSpecies = CIcore_numberOfConfigurationsRecursion2(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    numberOfConfigurations = 0
    ciLevel = 0

    !! call recursion to get the number of configurations...
    do ci = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0
      auxnumberOfSpecies = CIcore_numberOfConfigurationsRecursion(s, numberOfSpecies,  numberOfConfigurations, ciLevel) 

    end do

    call Vector_constructor8 ( CIcore_instance%diagonalHamiltonianMatrix2, &
                              numberOfConfigurations, 0.0_8 ) 

    CIcore_instance%numberOfConfigurations = numberOfConfigurations 

    write (*,*) "Number Of Configurations: ", numberOfConfigurations

    allocate ( indexConf ( numberOfSpecies ) )
    indexConf = 0

    !! calculate the diagonal 
    s = 0
    c = 0
    ciLevel = 0

    do ci = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      s = 0
      dd = 0

      u = CIcore_instance%auxciOrderList(ci)
      auxnumberOfSpecies = CIcore_buildDiagonalRecursion( s, numberOfSpecies, indexConf,  c, dd, u, cilevel, auxcilevel )
    end do
    !stop

    deallocate ( dd )
    deallocate ( indexConf )
    deallocate ( ciLevel )
    deallocate ( auxciLevel )

!$  timeB = omp_get_wtime()
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for Building diagonal of CI matrix : ", timeB - timeA ," (s)"

    write (*,*) "Reference energy, H_0: ",  CIcore_instance%diagonalHamiltonianMatrix2%values(1)

  end subroutine CIcore_buildDiagonal

recursive  function CIcore_numberOfConfigurationsRecursion(s, numberOfSpecies, c, cilevel) result (os)
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
        os = CIcore_numberOfConfigurationsRecursion( is, numberOfSpecies, c, cilevel )
      end do
    else 
      os = is

      i = cilevel(is) + 1
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
      end do
    end if

  end function CIcore_numberOfConfigurationsRecursion


recursive  function CIcore_buildDiagonalRecursion(s, numberOfSpecies, indexConf, c, dd, u, cilevel, auxcilevel) result (os)
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
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is) = ssize + a

        dd(is) =(a + CIcore_instance%ciOrderSize1(u,is))* CIcore_instance%ciOrderSize2(u,is) 
        os = CIcore_buildDiagonalRecursion( is, numberOfSpecies, indexConf, c, dd, u, cilevel, auxcilevel )
      end do
    else 
      os = is
      i = cilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1
        indexConf(is) = ssize + a
        !print *, indexConf
        dd(is) =(a + CIcore_instance%ciOrderSize1(u,is))* CIcore_instance%ciOrderSize2(u,is) 
        d = sum(dd)

        CIcore_instance%diagonalHamiltonianMatrix2%values(c) = &
                              CIcore_calculateEnergyZero ( indexConf )

      end do
    end if

  end function CIcore_buildDiagonalRecursion

  function CIcore_calculateEnergyZero( this ) result (auxCIenergy)
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

  end function CIcore_calculateEnergyZero


end module CIDiag_
