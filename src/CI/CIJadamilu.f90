 module CIJadamilu_
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
  subroutine CIJadamilu_buildCouplingMatrix()
    implicit none

    integer(8) :: a,b,c1,c2
    integer :: u,v,p
    integer :: i,n
    integer :: auxis,auxos
    integer :: numberOfSpecies
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(1), allocatable :: orbitalsA(:), orbitalsB(:)
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(1), allocatable :: couplingOrder(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    coupling = 0

    !! allocate arrays
    do n = 1, CIcore_instance%nproc
      do i = 1, numberOfSpecies

        call Matrix_constructorInteger ( CIcore_instance%couplingMatrix(i,n), &
          sum(CIcore_instance%numberOfStrings(i)%values), 3_8 , 0)

        call Matrix_constructorInteger(CIcore_instance%nCouplingOneTwo(i,n), &
          3_8, int(size(CIcore_instance%numberOfStrings(i)%values, dim=1),8),  0 )
  
        call Matrix_constructorInteger(CIcore_instance%nCouplingSize(i,n), &
          3_8, int(size(CIcore_instance%numberOfStrings(i)%values, dim=1) + 1 ,8),  0 )
  
        call Vector_constructor(CIcore_instance%couplingMatrixEnergyOne(i,n), &
          int(sum(CIcore_instance%numberOfStrings(i)%values),4), 0.0_8 )
  
        call Vector_constructorInteger(CIcore_instance%couplingMatrixFactorOne(i,n), &
          int(sum(CIcore_instance%numberOfStrings(i)%values),4), 2 )
  
        call Vector_constructorInteger( CIcore_instance%couplingMatrixOrbOne(i,n), &
          int(sum(CIcore_instance%numberOfStrings(i)%values),4), 0 )

      end do  
    end do  

  end subroutine CIJadamilu_buildCouplingMatrix

!! Build a list with all possible combinations of number of different orbitals from all quantum species, coupling (0,1,2)
  subroutine CIJadamilu_buildCouplingOrderList()
    implicit none

    integer(8) :: a,b,c,c1,c2,aa,d
    integer :: u,uu,vv, p, nn,z
    integer :: i
    integer :: numberOfSpecies, auxnumberOfSpecies,s
    integer(1), allocatable :: couplingOrder(:)
    integer(1) :: coupling
    real(8) :: timeA, timeB
    integer :: ncouplingOrderOne
    integer :: ncouplingOrderTwo
    integer :: ssize
    integer, allocatable :: cilevel(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    ssize = 1
    do i = 1, numberOfSpecies
      ssize = ssize * 3 !! ( 0,1,2) different orbitals
    end do

    allocate ( CIcore_instance%couplingOrderList( 3, ssize ) ) !! one, two same, two diff 
    allocate ( CIcore_instance%couplingOrderIndex( 3, ssize ) ) !! one, two same, two diff 

    do a = 1, 3
      do b = 1, ssize
        call Vector_constructorInteger1( CIcore_instance%couplingOrderList(a,b), &
          int( numberOfSpecies,8), int(0,1) )

      end do
    end do

    !! same species
    do b = 1, ssize
      call Vector_constructorInteger1( CIcore_instance%couplingOrderIndex(1,b), 1_8, int(0,1) )
      call Vector_constructorInteger1( CIcore_instance%couplingOrderIndex(2,b), 1_8, int(0,1) )
    end do

    !! diff species
    do b = 1, ssize
      call Vector_constructorInteger1( CIcore_instance%couplingOrderIndex(3,b), 2_8, int(0,1) )
    end do


    allocate ( couplingOrder ( numberOfSpecies )) !! 0, 1, 2
    couplingOrder = 0

    !! call recursion
    s = 0
    CIcore_instance%ncouplingOrderOne = 0
    CIcore_instance%ncouplingOrderTwo = 0
    CIcore_instance%ncouplingOrderTwoDiff = 0

    allocate ( ciLevel ( numberOfSpecies ) )
    ciLevel = 0

    !! get all combinations
    auxnumberOfSpecies = CIJadamilu_buildCouplingOrderRecursion( s, numberOfSpecies, couplingOrder, cilevel )

    !! save the index for species (speciesID) just to avoid a lot of conditionals later!

    do u = 1, CIcore_instance%ncouplingOrderOne
      do i = 1, numberOfSpecies
        if ( CIcore_instance%couplingOrderList(1,u)%values(i) == 1 ) then
          CIcore_instance%couplingOrderIndex(1,u)%values(1) = i
        end if
      end do
    end do

    do u = 1, CIcore_instance%ncouplingOrderTwo
      do i = 1, numberOfSpecies
        if ( CIcore_instance%couplingOrderList(2,u)%values(i) == 2 ) then
          CIcore_instance%couplingOrderIndex(2,u)%values(1) = i
        end if
      end do
    end do

    do u = 1, CIcore_instance%ncouplingOrderTwoDiff
      z = 0 
      do i = 1, numberOfSpecies
        if ( CIcore_instance%couplingOrderList(3,u)%values(i) == 1 ) then
          z = z + 1
          CIcore_instance%couplingOrderIndex(3,u)%values(z) = i
        end if
      end do
    end do


    deallocate ( ciLevel )
    deallocate ( couplingOrder ) 

  end subroutine CIJadamilu_buildCouplingOrderList


!! Get all possible combinations of number of different orbitals from all quantum species.
recursive  function CIJadamilu_buildCouplingOrderRecursion( s, numberOfSpecies, couplingOrder, cilevel ) result (os)
    implicit none

    integer(8) :: a,b,c,d
    integer :: u,v
    integer :: i, j, ii, jj, nn
    integer :: s, numberOfSpecies
    integer :: os,is,auxis, auxos
    integer(1) :: couplingOrder(:)
    logical :: same
    integer :: cilevel(:)

    is = s + 1
    if ( is < numberOfSpecies ) then
      if ( sum ( couplingOrder) <= 2 ) then
        do i = 1, 3 - sum ( couplingOrder ) !! 0,1,2
          couplingOrder(is) = i-1
          couplingOrder(is+1:) = 0
          os = CIJadamilu_buildCouplingOrderRecursion( is, numberOfSpecies, couplingOrder, cilevel )
        end do
      end if
    else 
      if ( sum ( couplingOrder) <= 2 ) then
        do i = 1, 3 - sum ( couplingOrder ) !! 0,1,2
          couplingOrder(is) = i-1
          couplingOrder(is+1:) = 0
          os = is
          if ( sum ( couplingOrder ) == 1 ) then

            auxis = 0
            CIcore_instance%ncouplingOrderOne = CIcore_instance%ncouplingOrderOne + 1
            b = CIcore_instance%ncouplingOrderOne
            CIcore_instance%couplingOrderList(1,b)%values = couplingOrder

          else if ( sum ( couplingOrder ) == 2 ) then

            same = .false. 

            do j = 1, numberOfSpecies
              if ( couplingOrder(j) == 2 ) same = .true.
            end do

            if ( same ) then
              auxis = 0
              CIcore_instance%ncouplingOrderTwo = CIcore_instance%ncouplingOrderTwo + 1
              b = CIcore_instance%ncouplingOrderTwo
              CIcore_instance%couplingOrderList(2,b)%values = couplingOrder
            else 
              auxis = 0
              CIcore_instance%ncouplingOrderTwoDiff = CIcore_instance%ncouplingOrderTwoDiff + 1
              b = CIcore_instance%ncouplingOrderTwoDiff
              CIcore_instance%couplingOrderList(3,b)%values = couplingOrder
            end if

          end if
        end do
      end if
    end if

  end function CIJadamilu_buildCouplingOrderRecursion

  subroutine CIJadamilu_jadamiluInterface(n,  maxeig, eigenValues, eigenVectors, timeA, timeB)
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
    !!if ( CONTROL_instance%CI_JACOBI ) then

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
!    the matrix is already in the required format

     IPRINT = 0 !     standard report on standard output
     ISEARCH = 1 !    we want the smallest eigenvalues
     NEIG = maxeig !    number of wanted eigenvalues
     !NINIT = 0 !    no initial approximate eigenvectors
     NINIT = NEIG !    initial approximate eigenvectors
     MADSPACE = maxsp !    desired size of the search space
     ITER = 1000*NEIG !    maximum number of iteration steps
     TOL = CONTROL_instance%CI_CONVERGENCE !1.0d-4 !    tolerance for the eigenvector residual

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
     if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
       print *, "Loading the eigenvector to the initial guess"
       do j = 1, n 
         X(j) = eigenVectors%values(j,1)
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = eigenValues%values(i)
       end do
     else
       jj = 0 
       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         jj = (i - 1) * n 
         do j = 1, CONTROL_instance%CI_SIZE_OF_GUESS_MATRIX 
          X(jj + CIcore_instance%auxIndexCIMatrix%values(j)) = CIcore_instance%initialEigenVectors%values(j,i)
         end do
       end do

       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
         EIGS(i) = CIcore_instance%initialEigenValues%values(i)
       end do
     end if

     DROPTOL = 0

     SIGMA = EIGS(1)
     gap = 0 
     SHIFT = EIGS(1)

     do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
       write(6,"(T2,A5,I4,2X,A10,F20.10,2X,A11,F20.10)") "State", i, "Eigenvalue", eigs( i ), "Eigenvector", x((i-1)*n + i)
     end do

     iiter = 0
  
!10     CALL DPJDREVCOM( N, A, JA, IA,EIGS, RES, X, LX, NEIG, &
!                        SIGMA, ISEARCH, NINIT, MADSPACE, ITER, TOL, &
!                        SHIFT, DROPTOL, MEM, ICNTL, &
!                        IJOB, NDX1, NDX2, IPRINT, INFO, GAP)
10   CALL DPJDREVCOM( N, CIcore_instance%diagonalHamiltonianMatrix%values , JA, IA, EIGS, RES, X, LX, NEIG, &
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
        if ( CONTROL_instance%CI_BUILD_FULL_MATRIX ) then
          call av ( n, x(ndx1), x(ndx2))
        else 
          call matvec2 ( N, X(NDX1), X(NDX2), iiter)
        end if

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

  end subroutine CIJadamilu_jadamiluInterface

  subroutine matvec2 ( nx, v, w, iter)
  
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
    real(8) :: CIEnergy
    integer(8) :: nonzero
    integer(8) :: i, j, ia, ib, ii, jj, iii, jjj
    integer(4) :: nproc, n, nn
    real(8) :: wi
    real(8) :: timeA, timeB
    real(8) :: tol
    integer(4) :: iter, size1, size2
    !integer(8), allocatable :: indexArray(:)
    logical :: fullMatrix
    integer :: ci
    integer :: auxSize
    integer(8) :: a,b,c
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    integer(8), allocatable :: cc(:) !! ncore
    integer(8), allocatable :: indexConf(:,:) !! ncore, species
    integer(8), allocatable :: auxindexConf(:,:) !! ncore, species
    integer, allocatable :: cilevel(:,:), auxcilevel(:,:)

    call omp_set_num_threads(omp_get_max_threads())
    nproc = omp_get_max_threads()


    allocate( cc ( nproc ) )
    cc = 0 

    nonzero = 0
    w = 0 
    tol = CONTROL_instance%CI_MATVEC_TOLERANCE 

    do i = 1 , nx
       if ( abs(v(i) ) >= tol) nonzero = nonzero + 1
    end do
  
    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    allocate ( indexConf ( numberOfSpecies, nproc ) )
    allocate ( auxindexConf ( numberOfSpecies, nproc ) )
    allocate ( cilevel ( numberOfSpecies, nproc ) )
    allocate ( auxcilevel ( numberOfSpecies, nproc ) )

    cilevel = 0
    auxcilevel = 0
    indexConf = 0
    auxindexConf = 0
    !! call recursion
    s = 0
    c = 0
    n = 1
!$  timeA = omp_get_wtime()
    do ci = 1,  CIcore_instance%sizeCiOrderList 
      do nn = n, nproc
        cilevel(:,nn) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
      end do
      s = 0 
      auxnumberOfSpecies = CIJadamilu_buildMatrixRecursion(nproc, s, indexConf, auxindexConf,cc, c, n, v, w, &
                             cilevel, auxcilevel )

    end do

    if  ( n > 1 ) then
       do nn = 1, n-1

       call CIJadamilu_buildRow( nn, auxindexConf(:,nn), cc(nn), w, v(cc(nn)), auxcilevel(:,nn))
      end do
    end if
    
    CIcore_instance%pindexConf = 0

!$  timeB = omp_get_wtime()
    deallocate ( cilevel )
    deallocate ( auxindexConf )
    deallocate ( indexConf )
    deallocate ( cc )
!$    write(*,"(A,I2,A,E10.3,A2,I12)") "  ", iter, "  ", timeB -timeA ,"  ", nonzero
!   stop
    return

  end subroutine matvec2

  subroutine av ( nx, v, w)
  
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
    character(50) :: CIFile
    integer :: CIUnit
    integer, allocatable :: jj(:)
    real(8), allocatable :: CIEnergy(:)
    integer :: nonzero,ii, kk
    integer :: maxStackSize, i, ia, ib

    CIFile = "lowdin.ci"
    CIUnit = 20
    nonzero = 0
    maxStackSize = CONTROL_instance%CI_STACK_SIZE 

    w = 0
#ifdef intel
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted", BUFFERED="YES")
#else
    open(unit=CIUnit, file=trim(CIFile), action = "read", form="unformatted")
#endif

    readmatrix : do
      read (CIUnit) nonzero
      if (nonzero > 0 ) then

        read (CIUnit) ii

        if ( allocated(jj)) deallocate (jj)
        allocate (jj(nonzero))
        jj = 0

        if ( allocated(CIEnergy)) deallocate (CIEnergy)
        allocate (CIEnergy(nonzero))
        CIEnergy = 0

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) jj(ia:ib)
    
        end do

        do i = 1, ceiling(real(nonZero) / real(maxStackSize) )

          ib = maxStackSize * i  
          ia = ib - maxStackSize + 1
          if ( ib >  nonZero ) ib = nonZero
          read (CIUnit) CIEnergy(ia:ib)
    
        end do

        w(ii) = w(ii) + CIEnergy(1)*v(jj(1)) !! disk
        do kk = 2, nonzero
          !w(ii) = w(ii) + CIcore_calculateCIenergy(ii,jj(kk))*v(jj(kk))  !! direct
          w(ii) = w(ii) + CIEnergy(kk)*v(jj(kk)) !! disk
          w(jj(kk)) = w(jj(kk)) + CIEnergy(kk)*v(ii) !! disk
        end do

      else if ( nonzero == -1 ) then
        exit readmatrix
      end if
    end do readmatrix

!! memory
!    do i = 1, nx
!        w(:) = w(:) + CIcore_instance%hamiltonianMatrix%values(:,i)*v(i)
!    end do 

     close(CIUnit)

    return
  end subroutine av

recursive  function CIJadamilu_buildMatrixRecursion(nproc, s, indexConf, auxindexConf, cc, c, n, v, w, &
                                                                  cilevel, auxcilevel) result (os)
    implicit none

    integer(8) :: a,c,aa
    integer :: i, n, nn, nproc 
    integer :: s, numberOfSpecies
    integer :: os,is,ss,ssize
    integer(8) :: cc(:)
    integer(8) :: indexConf(:,:)
    integer(8) :: auxindexConf(:,:)
    real(8) :: v(:)
    real(8) :: w(:)
    integer :: cilevel(:,:)
    integer :: auxcilevel(:,:)

    is = s + 1
    !if ( is < numberOfSpecies ) then
    do ss = 1, CIcore_instance%recursionVector1(is) 
      i = cilevel(is,n) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        indexConf(is,n:) = ssize + a
        os = CIJadamilu_buildMatrixRecursion( nproc, is, indexConf, auxindexConf, cc, c, n, v, w, cilevel, auxcilevel )
      end do
    end do
    !else 
    do ss = 1, CIcore_instance%recursionVector2(is) 
      os = is
      i = cilevel(is,n) + 1
      ssize = CIcore_instance%numberOfStrings2(is)%values(i)

      do a = 1, CIcore_instance%numberOfStrings(is)%values(i)
        c = c + 1

        if ( abs(v(c)) > CONTROL_instance%CI_MATVEC_TOLERANCE ) then
          cc(n) = c 
          indexConf(is,n:) = ssize + a

          auxindexConf = indexConf
          auxcilevel = cilevel

          if ( n == nproc ) then

            !$omp parallel &
            !$omp& private(nn),&
            !$omp& shared(v,w, indexConf, cc, nproc, cilevel) 
            !$omp do schedule (static) 
            do nn = 1, nproc
              call CIJadamilu_buildRow( nn, indexConf(:,nn), cc(nn), w, v(cc(nn)), cilevel(:,nn))
            end do
            !$omp end do nowait
            !$omp end parallel
            n = 0 

            do nn = 1, nproc
              indexConf(:,nn) = indexConf(:,nproc) 
              cilevel(:,nn) = cilevel(:,nproc) 
            end do
          end if 

          n = n + 1

        end if

      end do
    end do
    !end if


  end function CIJadamilu_buildMatrixRecursion

  !! Alternative option to the recursion with the same computational cost... However, it may be helpul some day. 

  function CIJadamilu_buildMatrixRecursion2(nproc, s, indexConf, auxindexConf, cc, c, n, v, w, &
                                                                  cilevel, auxcilevel) result (os)
    implicit none

    integer(8) :: a,c,aa, x
    integer :: i, j, n, nn, nproc, ci
    integer :: s, numberOfSpecies
    integer :: os,is,ss,ssize
    integer(8) :: cc(:)
    integer(8) :: indexConf(:,:)
    integer(8) :: auxindexConf(:,:)
    real(8) :: v(:)
    real(8) :: w(:)
    integer :: cilevel(:,:)
    integer(8) :: totalsize, auxtotalsize
    integer :: auxcilevel(:,:)
    integer, allocatable :: counter(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    
    allocate (counter(numberOfSpecies))
    counter = 0 

    totalsize = 1
    do i = 1 , numberOfSpecies
      totalsize = totalsize * CIcore_instance%numberOfStrings(i)%values(cilevel(i,n) + 1)
    end do

    do i = 1 , numberOfSpecies 
      ci = cilevel(i,n) + 1 
      ssize = CIcore_instance%numberOfStrings2(i)%values(ci)
      indexConf(i,n:) = ssize  + 1
    end do

    indexConf(numberOfSpecies,n:) = indexConf(numberOfSpecies,n:) -1

    do x = 1, totalsize

      indexConf(numberOfSpecies,n:) = indexConf(numberOfSpecies,n:) + 1

      do i = numberOfSpecies, 1 + 1, -1 
        auxtotalsize = 1
        do j = i, numberOfSpecies
          auxtotalsize = auxtotalsize * CIcore_instance%numberOfStrings(j)%values(cilevel(j,n) + 1)
        end do
        if (counter(i) == auxtotalsize) then
          do j = i, numberOfSpecies
            ci = cilevel(j,n) + 1 
            ssize = CIcore_instance%numberOfStrings2(j)%values(ci)
            indexConf(j,n:) = ssize + 1
          end do
          counter(i) = 0
          indexConf(i-1,n:) = indexConf(i-1,n:) + 1
        end if
        counter(i) = counter(i) + 1 

      end do
      !print *, indexConf(:,1)
    end do

    deallocate (counter)

  end function CIJadamilu_buildMatrixRecursion2

  subroutine CIJadamilu_buildRow( nn, indexConfA, c, w, vc, cilevelA)
    implicit none

    integer(8) :: a,b,c,bb,ci,d,cj
    integer :: u,v,uu,vv, p, nn
    integer :: i, j, auxis,auxos,is, ii, aa
    integer :: numberOfSpecies, s
    integer, allocatable :: stringAinB(:)
    integer(4) :: coupling
    integer(4) :: ssize,auxcoupling(3) !! 0,1,2
    integer(8) :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(8), allocatable :: dd(:)
    real(8) :: vc, CIenergy
    real(8) :: w(:)
    integer :: cilevelA(:)
    integer, allocatable :: cilevel(:)


    !CIcore_instance%pindexConf = 0

    !!$ CIcore_instance%timeA(1) = omp_get_wtime()

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 

    do i = 1, numberOfSpecies 

      if ( CIcore_instance%pindexConf(i,nn) /= indexConfA(i) ) then

      CIcore_instance%nCouplingOneTwo(i,nn)%values = 0
      auxcoupling = 0

      !allocate (stringBinA (CIcore_instance%numberOfOccupiedOrbitals%values(i) ))
      allocate (stringAinB (CIcore_instance%numberOfOccupiedOrbitals%values(i) ))

      stringAinB = 0
      !stringBinA = 0

      a = indexConfA(i)

    !!$ CIcore_instance%timeA(2) = omp_get_wtime()

      ssize = 0 
      do ci = 1,  size(CIcore_instance%numberOfStrings(i)%values, dim = 1)
        do b = 1 + ssize , CIcore_instance%numberOfStrings(i)%values(ci) + ssize

          !b = ssize + bb
          do p = CIcore_instance%numberOfCoreOrbitals%values(i)+1, &
                 CIcore_instance%numberOfOccupiedOrbitals%values(i)
          !do p = 1, &
          !       CIcore_instance%numberOfOccupiedOrbitals%values(i)

            stringAinB(p) = CIcore_instance%orbitals(i)%values( &
                              CIcore_instance%strings(i)%values(p,a),b) 

            !stringBinA(p) = CIcore_instance%orbitals(i)%values( &
            !                  CIcore_instance%strings(i)%values(p,b),a) 
          end do

          coupling = CIcore_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB ) - &
                     CIcore_instance%numberOfCoreOrbitals%values(i) 

         ! coupling = CIcore_instance%numberOfOccupiedOrbitals%values(i) - sum ( stringAinB ) 

          if ( coupling  <= 2 ) then

            coupling = coupling + 1

            auxcoupling(coupling) = auxcoupling(coupling) + 1 

            CIcore_instance%nCouplingOneTwo(i,nn)%values( coupling, ci) = &
              CIcore_instance%nCouplingOneTwo(i,nn)%values( coupling, ci) + 1

            CIcore_instance%couplingMatrix(i,nn)%values( auxcoupling(coupling), coupling ) = b
          end if

        end do

        ssize = ssize + CIcore_instance%numberOfStrings(i)%values(ci)

        end do

      deallocate (stringAinB)
      !deallocate (stringBinA)
      end if

    end do

    !!$ CIcore_instance%timeB(1) = omp_get_wtime()

    do is = 1, numberOfSpecies
      do i = 1, 3 !! 0,1,2
        ssize = 0
        do ci = 1,  size(CIcore_instance%numberOfStrings(is)%values, dim = 1) !! 1 is always zero
          ssize = ssize + CIcore_instance%nCouplingOneTwo(is,nn)%values( i,ci ) 
          CIcore_instance%nCouplingSize(is,nn)%values( i,ci+1 ) = ssize
         end do
        CIcore_instance%nCouplingSize(is,nn)%values( i,1 ) = 0 !0?
      end do
   end do


    !!$ CIcore_instance%timeA(2) = omp_get_wtime()
    allocate ( indexConfB ( numberOfSpecies ) )
    allocate ( cilevel ( numberOfSpecies ) )
    allocate ( dd ( numberOfSpecies ) )
    indexConfB = 0

    !!$ CIcore_instance%timeB(2) = omp_get_wtime()
    !!$ CIcore_instance%timeA(3) = omp_get_wtime()

    !!one diff same species
    do i = 1, numberOfSpecies

      if ( CIcore_instance%pindexConf(i,nn) /= indexConfA(i) ) then
      cilevel(:) = 0
      indexConfB = indexConfA

      cilevel = cilevelA

      do ci = 1,  size(CIcore_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        auxos = CIJadamilu_buildRowRecursionFirstOne( i, indexConfA, indexConfB, nn, cilevel )

      end do      
      end if
    end do      

    !!$ CIcore_instance%timeB(3) = omp_get_wtime()

    !!$ CIcore_instance%timeA(4) = omp_get_wtime()

    !$omp atomic
      w(c) = w(c) + vc*CIcore_instance%diagonalHamiltonianMatrix%values(c) 
    !$omp end atomic

    !!$ CIcore_instance%timeB(4) = omp_get_wtime()

    !!$ CIcore_instance%timeA(5) = omp_get_wtime()
    !! one diff
    do i = 1, numberOfSpecies
      cilevel(:) = 0
      indexConfB = indexConfA

      cilevel = cilevelA

      do ci = 1,  size(CIcore_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        do u = 1,  CIcore_instance%sizeciorderlist 
          if ( sum(abs(cilevel - &
               CIcore_instance%ciorderlist( CIcore_instance%auxciorderlist(u), :))) == 0 ) then

            uu = CIcore_instance%auxciorderlist(u)
            dd = 0

            auxos = CIJadamilu_buildRowRecursionSecondOne( i, indexConfB, w, vc, dd, nn, cilevel, uu )
            exit

          end if
        end do
      end do      
    end do      

    !!$ CIcore_instance%timeB(5) = omp_get_wtime()
    !!$ CIcore_instance%timeA(6) = omp_get_wtime()

    !! two diff same species
    do i = 1, numberOfSpecies

      cilevel(:) = 0
      indexConfB = indexConfA
      cilevel = cilevelA

      do ci = 1,  size(CIcore_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
        cilevel(i) = ci - 1

        do u = 1,  CIcore_instance%sizeCiOrderList 
          if ( sum(abs(cilevel - &
               CIcore_instance%ciOrderList( CIcore_instance%auxciOrderList(u), :))) == 0 ) then
            uu = CIcore_instance%auxciOrderList(u)
            dd = 0

            if ( CIcore_instance%pindexConf(i,nn) /= indexConfA(i) ) then
              auxos = CIJadamilu_buildRowRecursionSecondTwoCal( i, indexConfA, indexConfB, w, vc, dd, nn, cilevel, uu )
            else
              auxos = CIJadamilu_buildRowRecursionSecondTwoGet( i, indexConfA, indexConfB, w, vc, dd, nn, cilevel, uu )
            end if

            exit

          end if
        end do
      end do
    end do      

    !!$ CIcore_instance%timeB(6) = omp_get_wtime()
    !!$ CIcore_instance%timeA(7) = omp_get_wtime()

    !! two diff diff species
    do v = 1, CIcore_instance%ncouplingOrderTwoDiff

       i = CIcore_instance%couplingOrderIndex(3,v)%values(1)
       j = CIcore_instance%couplingOrderIndex(3,v)%values(2)

       indexConfB = indexConfA
       cilevel = cilevelA

        do ci = 1,  size(CIcore_instance%numberOfStrings(i)%values, dim = 1) !! 1 is always zero
          cilevel(i) = ci - 1
          do cj = 1,  size(CIcore_instance%numberOfStrings(j)%values, dim = 1) !! 1 is always zero
            cilevel(j) = cj - 1
            do u = 1,  CIcore_instance%sizeCiOrderList 
              if ( sum(abs(cilevel - &
                   CIcore_instance%ciOrderList( CIcore_instance%auxciOrderList(u), :))) == 0 ) then

                uu = CIcore_instance%auxciOrderList(u)
                dd = 0
                auxos = CIJadamilu_buildRowRecursionSecondTwoDiff( i, j, indexConfB, w, vc, dd, nn, cilevel, uu )
                exit
              end if
            end do
          end do
        end do
    end do    

    !!$ CIcore_instance%timeB(7) = omp_get_wtime()

    !!$ print *, "omptime"
    !!$ print *, "1", CIcore_instance%timeB(1) - CIcore_instance%timeA(1)
    !!$ print *, "2", CIcore_instance%timeB(2) - CIcore_instance%timeA(2)
    !!$ print *, "3", CIcore_instance%timeB(3) - CIcore_instance%timeA(3)
    !!$ print *, "4", CIcore_instance%timeB(4) - CIcore_instance%timeA(4)
    !!$ print *, "5", CIcore_instance%timeB(5) - CIcore_instance%timeA(5)
    !!$ print *, "6", CIcore_instance%timeB(6) - CIcore_instance%timeA(6)
    !!$ print *, "7", CIcore_instance%timeB(7) - CIcore_instance%timeA(7)

    CIcore_instance%pindexConf(:,nn) = indexConfA(:)

    deallocate ( dd )
    deallocate ( cilevel )
    deallocate ( indexConfB )

  end subroutine CIJadamilu_buildRow

recursive  function CIJadamilu_buildRowRecursionFirstOne( ii, indexConfA, indexConfB, nn, cilevel ) result (os)
    implicit none

    integer(8) :: a, aa
    integer :: ii, nn, ci
    integer :: os, ssize
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

      ci = cilevel(ii) + 1
      ssize = CIcore_instance%nCouplingSize(ii,nn)%values( 2,ci ) 
      do aa = 1, CIcore_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
        a = ssize + aa

        indexConfB(ii) = CIcore_instance%couplingMatrix(ii,nn)%values(a, 2)
        CIenergy = CIJadamilu_calculateEnergyOneSame ( nn, ii, indexConfA, indexConfB )
        CIcore_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) = CIenergy

      end do

  end function CIJadamilu_buildRowRecursionFirstOne
 
recursive  function CIJadamilu_buildRowRecursionSecondOne( ii, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: ii, nn, ci, u, j
    integer :: ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    ci = cilevel(ii) + 1
    ssize = CIcore_instance%nCouplingSize(ii,nn)%values( 2,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - CIcore_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
                 CIcore_instance%ciOrderSize1(u,j) )* CIcore_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, CIcore_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
      a = ssize + aa

      indexConfB(ii) = CIcore_instance%couplingMatrix(ii,nn)%values(a, 2)

      dd(ii) = (indexConfB(ii) - CIcore_instance%numberOfStrings2(ii)%values(ci) + &
                 CIcore_instance%ciOrderSize1(u,ii) )* CIcore_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      CIenergy = CIcore_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      CIenergy = CIenergy + CIJadamilu_calculateEnergyOneDiff ( ii, indexConfB, nn )
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function CIJadamilu_buildRowRecursionSecondOne


  function CIJadamilu_buildRowRecursionSecondTwoCal( ii, indexConfA, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: i, ii, nn, ci, u, j
    integer :: s, ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    ci = cilevel(ii) + 1
    ssize = CIcore_instance%nCouplingSize(ii,nn)%values( 3,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - CIcore_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
               CIcore_instance%ciOrderSize1(u,j) )* CIcore_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, CIcore_instance%nCouplingOneTwo(ii,nn)%values( 3,ci ) 
      a = ssize + aa

      indexConfB(ii) = CIcore_instance%couplingMatrix(ii,nn)%values(a, 3)
      dd(ii) = (indexConfB(ii) - CIcore_instance%numberOfStrings2(ii)%values(ci) + &
               CIcore_instance%ciOrderSize1(u,ii) )* CIcore_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      !CIenergy = CIcore_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      CIenergy = CIJadamilu_calculateEnergyTwoSame ( ii, indexConfA(ii), indexConfB(ii) )
      CIcore_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) = CIenergy
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function CIJadamilu_buildRowRecursionSecondTwoCal

  function CIJadamilu_buildRowRecursionSecondTwoGet( ii, indexConfA, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer(8) :: a,d, aa
    integer :: i, ii, nn, ci, u, j
    integer :: s, ssize
    integer :: os,numberOfSpecies
    integer(8) :: indexConfA(:)
    integer(8) :: indexConfB(:)
    integer(8) :: dd(:)
    real(8) :: vc
    real(8) :: w(:)
    real(8) :: CIenergy
    integer :: cilevel(:)

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    ci = cilevel(ii) + 1
    ssize = CIcore_instance%nCouplingSize(ii,nn)%values( 3,ci ) 

    do j = 1, numberOfSpecies
      dd(j) = (indexConfB(j) - CIcore_instance%numberOfStrings2(j)%values(cilevel(j)+1) + &
                 CIcore_instance%ciOrderSize1(u,j) )* CIcore_instance%ciOrderSize2(u,j) 
    end do

    do aa = 1, CIcore_instance%nCouplingOneTwo(ii,nn)%values( 3,ci ) 
      a = ssize + aa

      indexConfB(ii) = CIcore_instance%couplingMatrix(ii,nn)%values(a, 3)
      dd(ii) = (indexConfB(ii) - CIcore_instance%numberOfStrings2(ii)%values(ci) + &
               CIcore_instance%ciOrderSize1(u,ii) )* CIcore_instance%ciOrderSize2(u,ii) 

      d = sum(dd)

      CIenergy = CIcore_instance%couplingMatrixEnergyOne(ii,nn)%values(indexConfB(ii)) 
      !CIenergy = CIcore_calculateEnergyTwoSame ( ii, indexConfA(ii), indexConfB(ii) )
      CIenergy = vc*CIenergy 

      !$omp atomic
      w(d) = w(d) + CIenergy 
      !$omp end atomic
    end do

  end function CIJadamilu_buildRowRecursionSecondTwoGet

 function CIJadamilu_buildRowRecursionSecondTwoDiff( ii, jj, indexConfB, w, vc, dd, nn, cilevel, u ) result (os)
    implicit none

    integer, intent(in) :: ii, nn, u, jj
    integer, intent(in) :: cilevel(:)
    integer(8), intent(out) :: dd(:)
    real(8), intent(in) :: vc
    integer(8), intent(inout) :: indexConfB(:)
    real(8), intent(inout) :: w(:)
    integer(8) :: ai,aj,d, aai, aaj
    integer :: ci, k, cj
    integer(8) :: ssizei, ssizej
    integer(8) :: dd_i_shift, dd_j_shift
    integer :: bi, bj, factor, factori
    integer :: auxIndex1, auxIndex2, auxIndex
    integer :: os,numberOfSpecies
    real(8) :: CIenergy

    numberOfSpecies = CIcore_instance%numberOfQuantumSpecies 
    ci = cilevel(ii) + 1
    cj = cilevel(jj) + 1
    ssizei = CIcore_instance%nCouplingSize(ii,nn)%values( 2,ci ) 
    ssizej = CIcore_instance%nCouplingSize(jj,nn)%values( 2,cj ) 

    do k = 1, numberOfSpecies
        dd(k) = (indexConfB(k) - CIcore_instance%numberOfStrings2(k)%values(cilevel(k)+1) + &
                CIcore_instance%ciOrderSize1(u,k) )* CIcore_instance%ciOrderSize2(u,k) 
    end do

    dd_i_shift  = - CIcore_instance%numberOfStrings2(ii)%values(ci) + &
                CIcore_instance%ciOrderSize1(u,ii) 

    dd_j_shift =  - CIcore_instance%numberOfStrings2(jj)%values(cj) + &
                CIcore_instance%ciOrderSize1(u,jj)

    do aai = 1, CIcore_instance%nCouplingOneTwo(ii,nn)%values( 2,ci ) 
      ai = ssizei + aai
      indexConfB(ii) = CIcore_instance%couplingMatrix(ii,nn)%values(ai, 2)
      dd(ii) = (indexConfB(ii) + dd_i_shift )* CIcore_instance%ciOrderSize2(u,ii) 

      bi = indexConfB(ii)
      factori = CIcore_instance%couplingMatrixFactorOne(ii,nn)%values(bi) 
      auxIndex1 = CIcore_instance%couplingMatrixOrbOne(ii,nn)%values(bi) 
      auxIndex1 = CIcore_instance%numberOfSpatialOrbitals2%values(jj) * (auxIndex1 - 1 ) 

      !$omp simd
      do aaj = 1, CIcore_instance%nCouplingOneTwo(jj,nn)%values( 2,cj ) 
        aj = ssizej + aaj
        indexConfB(jj) = CIcore_instance%couplingMatrix(jj,nn)%values(aj, 2)

        dd(jj) = (indexConfB(jj) + dd_j_shift )* CIcore_instance%ciOrderSize2(u,jj) 

        d = sum(dd)
          !CIenergy = vc*CIcore_calculateEnergyTwoDiff ( ii, jj, indexConfB, nn )

        bj = indexConfB(jj)
        factor = factori * CIcore_instance%couplingMatrixFactorOne(jj,nn)%values(bj) 
        auxIndex2 = CIcore_instance%couplingMatrixOrbOne(jj,nn)%values(bj) 
        auxIndex = auxIndex1 + auxIndex2

        CIenergy = vc * factor *CIcore_instance%fourCenterIntegrals(ii,jj)%values(auxIndex, 1)
        !CIenergy = vc*CIenergy 

        !$omp atomic
        w(d) = w(d) + CIenergy 
        !$omp end atomic
      end do
    end do

  end function CIJadamilu_buildRowRecursionSecondTwoDiff


  function CIJadamilu_calculateEnergyOneSame( n, ii, thisA, thisB ) result (auxCIenergy)
    implicit none
    integer(8) :: thisA(:), thisB(:)
    integer(8) :: a, b
    integer :: i,j,s,n, nn,ii
    integer :: l,k,z,kk,ll
    integer :: factor, factor2, auxOcc, AA, BB
    logical(1) :: equalA, equalB
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer(8) :: auxIndex1, auxIndex2, auxIndex
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

    CIcore_instance%couplingMatrixFactorOne(ii,n)%values(b) = factor

    !One particle terms

    auxCIenergy= auxCIenergy +  CIcore_instance%twoCenterIntegrals(ii)%values( diffOrb(1), diffOrb(2) )

    !! save the different orbitals

    auxIndex1= CIcore_instance%twoIndexArray(ii)%values( diffOrb(1), diffOrb(2))
    CIcore_instance%couplingMatrixOrbOne(ii,n)%values(b) = auxIndex1

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

    auxCIenergy= auxCIenergy * factor

  end function CIJadamilu_calculateEnergyOneSame

  function CIJadamilu_calculateEnergyOneDiff( ii, thisB, nn ) result (auxCIenergy)
    implicit none
    integer(8) :: thisB(:)
    integer(8) :: b
    integer :: i,j,ii, nn
    integer :: l,ll
    integer :: factor
    integer :: auxnumberOfOtherSpecieSpatialOrbitals
    integer :: auxIndex1, auxIndex11, auxIndex
    real(8) :: auxCIenergy

    auxCIenergy = 0.0_8

    b = thisB(ii)

    auxIndex1 = CIcore_instance%couplingMatrixOrbOne(ii,nn)%values(b) 
    factor = CIcore_instance%couplingMatrixFactorOne(ii,nn)%values(b) 

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

  end function CIJadamilu_calculateEnergyOneDiff


  function CIJadamilu_calculateEnergyTwoSame( ii, a, b ) result (auxCIenergy)
    implicit none
    integer(8) :: a, b
    integer :: ii
    integer :: kk,z
    integer :: factor, AA(2), BB(2)
    integer(8) :: auxIndex
    integer :: diffOrbA(2), diffOrbB(2)  !! to avoid confusions
    real(8) :: auxCIenergy

    !diffOrbA = 0
    !diffOrbB = 0
    z = 0

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

  end function CIJadamilu_calculateEnergyTwoSame

end module CIJadamilu_
