module CIStrings_
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

  subroutine CIStrings_buildStrings()
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

      call Vector_constructorInteger8 (CIcore_instance%numberOfStrings(i), &
        int(CIcore_instance%CILevel(i) + 1,8), 0_8)

      CIcore_instance%numberOfStrings(i)%values(1) = 1 !! ground

      write (*,"(A,A)") "  ", MolecularSystem_getSymbolOfSpecies(i)

      do cilevel = 1,CIcore_instance%CILevel(i) 

        call Vector_constructor (occupiedCode(i), cilevel, real(CIcore_instance%numberOfCoreOrbitals%values(i),8) )
        call Vector_constructor (unoccupiedCode(i), cilevel, 0.0_8)

        unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i)  ! it's also a lower bound in a for loop

        if ( cilevel <= CIcore_instance%numberOfOccupiedOrbitals%values(i) ) then

          !! just get the number of strings...
          ci = 0 
          oci = CIStrings_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel)

          write (*,"(A,I4,I12)") "    ", cilevel, CIcore_instance%numberOfStrings(i)%values(cilevel+1)

        end if
      end do
      write (*,"(A,I8)") "  Total:", sum(CIcore_instance%numberOfStrings(i)%values)
      write (*,"(A)") ""

      !! allocate the strings arrays
      if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) > 0 ) then
        call Matrix_constructorInteger( CIcore_instance%strings(i), &
          int(CIcore_instance%numberOfOccupiedOrbitals%values(i),8), &
          sum(CIcore_instance%numberOfStrings(i)%values), int(0,4))

        call Matrix_constructorInteger1( CIcore_instance%orbitals(i), &
          int(CIcore_instance%numberOfOrbitals%values(i),8), &
          sum(CIcore_instance%numberOfStrings(i)%values), 0_1)

      else
        call Matrix_constructorInteger( CIcore_instance%strings(i), &
         1_8, 1_8, int(0,4))
        call Matrix_constructorInteger1( CIcore_instance%orbitals(i), &
         1_8, 1_8, 0_1)

      end if

      !! zero, build the reference
      call Vector_constructorInteger (order, numberOfSpecies, 0 )

      call Vector_constructor (occupiedCode(i), 1, 0.0_8) !! initialize in zero
      call Vector_constructor (unoccupiedCode(i), 1, 0.0_8)

      c = 0 
      c = c + 1
      call Configuration_constructorB(CIcore_instance%strings(i), CIcore_instance%orbitals(i), &
                                      occupiedCode, unoccupiedCode, i, c, order)
      
      !! now build the strings
      do cilevel = 1,CIcore_instance%CILevel(i) 

        call Vector_constructorInteger (order, numberOfSpecies, 0 )
        order%values(i) = cilevel

        call Vector_constructor (occupiedCode(i), cilevel, real(CIcore_instance%numberOfCoreOrbitals%values(i),8) )
        call Vector_constructor (unoccupiedCode(i), cilevel, 0.0_8)

        unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i)  ! it's also a lower bound in a for loop

        if ( cilevel <= CIcore_instance%numberOfOccupiedOrbitals%values(i) ) then

          !! recursion to build the strings
          ci = 0 
          oci = CIStrings_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c)

        end if
      end do

    end do

    !! useful array  
    do i = 1, numberOfSpecies
      CIcore_instance%sumstrings(i) = sum(CIcore_instance%numberOfStrings(i)%values)
    end do

    !! useful array, save the total number of string for a previous CI level. 
    do i = 1, numberOfSpecies
      call Vector_constructorInteger8 (CIcore_instance%numberOfStrings2(i), &
        int(size(CIcore_instance%numberOfStrings(i)%values, dim = 1) + 1,8), 0_8)

      ssize = 0
      do j = 1, size(CIcore_instance%numberOfStrings(i)%values, dim = 1) !
        ssize = ssize + CIcore_instance%numberOfStrings(i)%values(j)
        CIcore_instance%numberOfStrings2(i)%values(j+1) = ssize 
      end do
      CIcore_instance%numberOfStrings2(i)%values(1) = 0
    end do
      

  end subroutine CIStrings_buildStrings

!! This is just to get the total number of strings...
recursive  function CIStrings_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ici, cilevel ) result (oci)
    implicit none

    integer :: i, numberOfSpecies
    integer :: ci, ici, oci, cilevel
    integer :: m, a
    type(vector), allocatable :: occupiedCode(:)
    type(vector), allocatable :: unoccupiedCode(:)

    ci = ici + 1

    if ( ci == 1 .and. ci < cilevel ) then ! first
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = CIStrings_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel )
        end do
        unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    else if ( ci > 1 .and. ci < cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = CIStrings_buildStringsRecursion( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel )
        end do
      end do

    else if ( ci == 1 .and. ci == cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          CIcore_instance%numberOfStrings(i)%values(ci+1) = &
            CIcore_instance%numberOfStrings(i)%values(ci+1) + 1
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do

    else !final

      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          CIcore_instance%numberOfStrings(i)%values(ci+1) = &
            CIcore_instance%numberOfStrings(i)%values(ci+1) + 1
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    end if

  end function CIStrings_buildStringsRecursion

!! and this is for building the strings
recursive  function CIStrings_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, &
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
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = CIStrings_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c )
        end do
        unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    else if ( ci > 1 .and. ci < cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          oci = CIStrings_buildStringsRecursion2( i, numberOfSpecies, occupiedCode, unoccupiedCode, ci, cilevel, order, c )
        end do
      end do

    else if ( ci == 1 .and. ci == cilevel ) then ! mid
      do m = int(occupiedCode(i)%values(ci)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a

          c = c + 1
          call Configuration_constructorB(CIcore_instance%strings(i), CIcore_instance%orbitals(i), &
                                          occupiedCode, unoccupiedCode, i, c, order)
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do

    else !final

      do m = int(occupiedCode(i)%values(ci-1)) + 1,  int(CIcore_instance%numberOfOccupiedOrbitals%values(i))
        do a = int(unoccupiedCode(i)%values(ci-1)) + 1, int(CIcore_instance%numberOfOrbitals%values(i) )
          occupiedCode(i)%values(ci) = m
          unoccupiedCode(i)%values(ci) = a
          c = c + 1
          call Configuration_constructorB(CIcore_instance%strings(i), CIcore_instance%orbitals(i), &
                                          occupiedCode, unoccupiedCode, i, c, order)
        end do
        if ( ci == 1 ) unoccupiedCode(i)%values = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
    end if


  end function CIStrings_buildStringsRecursion2

end module CIStrings_
