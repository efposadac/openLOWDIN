module CIOrder_
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
  subroutine CIOrder_settingCILevel()
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
    do n = 1, CIcore_instance%nproc
      do i = 1, numberOfSpecies
        call Vector_constructorInteger( CIcore_instance%auxstring(n,i), &
          int(CIcore_instance%numberOfOccupiedOrbitals%values(i),4), int(0,4))
      end do  
    end do  

    select case ( trim(CIcore_instance%level) )

    case ( "FCI" )

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i)
      end do

      CIcore_instance%maxCILevel = sum(CIcore_instance%CILevel)

    case ( "SCI" ) !! same as FCI

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i)
      end do

      CIcore_instance%maxCILevel = sum(CIcore_instance%CILevel)

    case ( "CIS" )

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 1
      end do
      CIcore_instance%maxCILevel = 1

    case ( "CISD" )

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 2
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 2

    case ( "CISD+" )

      if ( .not. numberOfSpecies == 3 ) call CIOrder_exception( ERROR, "CIOrder setting CI level ", "CISD+ is specific for three quantum species")

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 2
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 2

    case ( "CISD+2" )

      if ( .not. numberOfSpecies == 4 ) call CIOrder_exception( ERROR, "CIOrder setting CI level", "CISD+2 is specific for three quantum species")
      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 2
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 2 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 2

    case ("CISDT")

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 3
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 3 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 3

    case ("CISDTQ")

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 4
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 4 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 4

    case ("CISDTQQ")

      do i=1, numberOfSpecies
        CIcore_instance%CILevel(i) = 5
        if ( CIcore_instance%numberOfOccupiedOrbitals%values(i) < 5 ) &
          CIcore_instance%CILevel(i) = CIcore_instance%numberOfOccupiedOrbitals%values(i) 
      end do
      CIcore_instance%maxCILevel = 5

    case default

       call CIOrder_exception( ERROR, "Configuration interactor constructor", "Correction level not implemented")

    end select

    if ( CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT == "CISD" .and. trim(CIcore_instance%level)  /= "CISD" ) then

       call CIOrder_exception( ERROR, "Configuration interactor constructor", "DDCISD shift are only valid for CISD level!")

    end if


  end subroutine CIOrder_settingCILevel




!! Build the CI table with all combinations of excitations between quantum species.
  subroutine CIOrder_buildCIOrderList()
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
       ssize = ssize * (CIcore_instance%CILevel(i) + 1)
    end do

    allocate ( CIcore_instance%ciOrderList( ssize, numberOfSpecies ) ) 
    allocate ( CIcore_instance%ciOrderSize1( ssize, numberOfSpecies ) ) 
    allocate ( CIcore_instance%ciOrderSize2( ssize, numberOfSpecies ) ) 
    allocate ( CIcore_instance%auxciOrderList( ssize ) ) 

    CIcore_instance%ciOrderList = 0
    CIcore_instance%auxciOrderList = 0

    CIcore_instance%ciOrderSize1 = -1 !! I have reasons... -1 for all species except the last one
    CIcore_instance%ciOrderSize2 = 1 !! and 1 for the last species

    CIcore_instance%sizeCiOrderList = 0

    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( auxciLevel ( numberOfSpecies ) )
    ciLevel = 0
    auxciLevel = 0
    s = 0
    c = 0
    !! Search which combinations of excitations satifies the desired CI level.
    auxnumberOfSpecies = CIOrder_buildCIOrderRecursion( s, numberOfSpecies, c, cilevel )


    !! Print list
    write (6,"(T2,A)") "--------------------------"
    write (6,"(T2,A)") "CI level \ Species"
    write (6,"(T2,A)") "--------------------------"
    do u = 1,  CIcore_instance%sizeCiOrderList 
      do i = 1, numberOfSpecies
        write (6,"(T2,I4)",advance="no") CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(u), i)
      end do
      write (6,"(A)") ""
    end do
    write (6,"(T2,A)") "--------------------------"

    !! Calculates the three required factors in order to get the position of any given configuration.
    !! position = S1 + (indexConf(i,u) - numberOfStrings2(i) -1 )*S2(i,u)  
    !! i: speciesID, u: cilevelID

    !! Factor S1
    ssize = 0
    do u = 1,  CIcore_instance%sizeCiOrderList 

      cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(u), :)

      ssize = 0
      do v = 1,  u-1

        auxcilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(v), :)
        auxnumberOfSpecies = CIOrder_getIndexSize(0, ssize, auxcilevel) 

      end do

      CIcore_instance%ciOrderSize1(CIcore_instance%auxciOrderList(u),:) = -1
      CIcore_instance%ciOrderSize1(CIcore_instance%auxciOrderList(u),numberOfSpecies) = ssize !!just the last

    end do

    !! Factor S2
    do i = 1, numberOfSpecies-1
      do u = 1,  CIcore_instance%sizeCiOrderList 

        cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(u), :)
        ssize = 1
        do j = i+1, numberOfSpecies
          ssize = ssize * CIcore_instance%numberOfStrings(j)%values(cilevel(j)+1)
        end do

        CIcore_instance%ciOrderSize2(CIcore_instance%auxciOrderList(u),i) = ssize

      end do
    end do

    CIcore_instance%ciOrderSize2(:,numberOfSpecies) = 1 

    deallocate ( auxcilevel )
    deallocate ( cilevel )
    
  end subroutine CIOrder_buildCIOrderList

    !! Search which combinations of excitations satifies the desired CI level.
recursive  function CIOrder_buildCIOrderRecursion( s, numberOfSpecies, c, cilevel ) result (os)
    implicit none

    integer :: u,v,c
    integer :: i, j, ii, jj, nn, k, l
    integer :: s, numberOfSpecies
    integer :: os,is,auxis, auxos
    integer :: cilevel(:)
    integer :: plusOne(3,3) , plusTwo(4,6)

    is = s + 1
    if ( is < numberOfSpecies ) then
      do i = 1, size(CIcore_instance%numberOfStrings(is)%values, dim = 1)
       cilevel(is) = i - 1
       os = CIOrder_buildCIOrderRecursion( is, numberOfSpecies, c, cilevel )
      end do
      cilevel(is) = 0
    else 
      do i = 1, size(CIcore_instance%numberOfStrings(is)%values, dim = 1)
       cilevel(is) = i - 1
       c = c + 1

       CIcore_instance%ciOrderList( c, : ) = cilevel(:)
       if ( sum(cilevel) <= CIcore_instance%maxCIlevel ) then
         CIcore_instance%sizeCiOrderList = CIcore_instance%sizeCiOrderList + 1
         CIcore_instance%auxciOrderList(  CIcore_instance%sizeCiOrderList  ) = c
       end if

       if ( trim(CIcore_instance%level) == "CISD+" ) then !!special case. 
         plusOne(:,1) = (/1,1,1/)
         plusOne(:,2) = (/2,0,1/)
         plusOne(:,3) = (/0,2,1/)
       
         do k = 1, 3
           if ( sum(  abs(cilevel(:) - plusOne(:,k)) ) == 0 ) then
           CIcore_instance%sizeCiOrderList = CIcore_instance%sizeCiOrderList + 1
           CIcore_instance%auxciOrderList(  CIcore_instance%sizeCiOrderList  ) = c
           end if
         end do
       
       end if
       
       if ( trim(CIcore_instance%level) == "CISD+2" ) then !!special case. 
         plusTwo(:,1) = (/1,1,1,0/)
         plusTwo(:,2) = (/1,1,0,1/)
         plusTwo(:,3) = (/2,0,1,0/)
         plusTwo(:,4) = (/2,0,0,1/)
         plusTwo(:,5) = (/0,2,1,0/)
         plusTwo(:,6) = (/0,2,0,1/)

         do k = 1, 6
           if ( sum(  abs(cilevel(:) - plusTwo(:,k)) ) == 0 ) then
           CIcore_instance%sizeCiOrderList = CIcore_instance%sizeCiOrderList + 1
           CIcore_instance%auxciOrderList(  CIcore_instance%sizeCiOrderList  ) = c
           end if
         end do
 
       end if

      end do
      cilevel(is) = 0
    end if

  end function CIOrder_buildCIOrderRecursion

recursive  function CIOrder_getIndexSize(s, c, auxcilevel) result (os)
    implicit none

    integer(8) :: a,b,c
    integer :: u,v
    integer :: i, j, ii, jj, ss
    integer :: s, numberOfSpecies
    integer :: os,is,cc, ssize
    integer :: auxcilevel(:)

    is = s + 1
    do ss = 1, CIcore_instance%recursionVector1(is) 
      i = auxcilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        os = CIOrder_getIndexSize( is, c, auxcilevel )
      end do
    end do
    do ss = 1, CIcore_instance%recursionVector2(is) 
      os = is
      i = auxcilevel(is) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)
      c = c + CIcore_instance%numberOfStrings(is)%values(i)
    end do

  end function CIOrder_getIndexSize

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine CIOrder_exception( typeMessage, description, debugDescription)
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

  end subroutine CIOrder_exception


end module CIOrder_
