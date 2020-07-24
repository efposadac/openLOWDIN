!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Clase encargada de realizar transformacion de integrales atomicas a  moleculares
!!
!!  Esta clase reliza la transformacion de integrales de orbitales atomicos a orbitales moleculares,
!!  creando una interface al algoritmo de   Yamamoto, Shigeyoshi; Nagashima, Umpei.
!!  Computer Physics Communications, 2005, 166, 58-65
!!
!! @author Sergio Gonzalez
!!
!! <b> Fecha de creacion : </b> 2009-07-07
!!   - <tt> 2009-07-07 </tt>: Sergio Gonzalez ( sagonzalez@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)     
!!        -# Adapts to Lowdin 2               
!!   - <tt> 2014-08-26 </tt>: Jorge Charry (jacharrym@unal.edu.co)     
!!        -# Adapts this module to works indepently from MP2 program
!<
module TransformIntegralsE_
  use MolecularSystem_
  use InputCI_
  use LapackInterface_
  use InputManager_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
!  use omp_lib
  implicit none

  type, public :: TransformIntegralsE
     character(30) :: name
     character(255) :: fileForCoefficients
     character(255) :: fileForIntegrals
     character(255) :: prefixOfFile
     integer :: numberOfContractions
     integer :: otherNumberOfContractions
     integer :: bias
     integer :: specieID
     integer :: otherSpecieID
     integer :: unidOfOutputForCoefficients
     integer :: unidOfOutputForIntegrals
     integer :: nproc
     integer :: integralStackSize

     integer :: p_l, p_u
     integer :: q_l, q_u
     integer :: r_l, r_u
     integer :: s_l, s_u

     character(50) :: partialTransform
     
  end type TransformIntegralsE

!  interface
!      subroutine TransformIntegralsE_setmem( ssize, twoParticlesIntegrals ) &
!        bind(C, name"c_TransformIntegralsE_setmem")
!    use, intrinsic :: iso_c_binding
!    implicit none 
!    integer :: ssize
!    real(8), allocatable :: twoParticlesIntegrals(:)
!    integer(8) :: ssize, ssize2, ssize4
!
!    ssize2 = (ssize * (ssize + 1))/2 
!    ssize4 = (ssize2 * (ssize2 + 1))/2 
!
!    allocate (twoParticlesIntegrals ( ssize4 ) )
!
!    twoParticlesIntegrals = 0
!    end subroutine TransformIntegralsE_setmem
!  end interface





  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE = 0
  integer, parameter :: TWO_SPECIES = 1
  !! }

  public :: &
       TransformIntegralsE_constructor, &
       TransformIntegralsE_destructor, &
       TransformIntegralsE_show, &
       TransformIntegralsE_atomicToMolecularOfOneSpecie, &
       TransformIntegralsE_atomicToMolecularOfTwoSpecies
  !       TransformIntegralsE_readIntegralsTransformed

  private

contains


  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsE_constructor(this,partial)
    implicit none
    type(TransformIntegralsE) :: this
    character(*) :: partial

    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"

    this%partialTransform=trim(partial)



  end subroutine TransformIntegralsE_constructor

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsE_destructor(this)
    implicit none
    type(TransformIntegralsE) :: this

  end subroutine TransformIntegralsE_destructor

  !>
  !! @brief show
  !<
  subroutine TransformIntegralsE_show()
    implicit none

    print *,"--------------------------------------------------"
    print *,"   Two-half Four-index integral tranformation     "
    print *,"--------------------------------------------------"

  end subroutine TransformIntegralsE_show

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsE_atomicToMolecularOfOneSpecie( this, coefficientsOfAtomicOrbitals, &
       molecularIntegrals, specieID, nameOfSpecie )
    implicit none
    type(TransformIntegralsE) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID
    character(*) :: nameOfSpecie
    integer :: integralStackSize

    integer :: status
    integer(8) :: ssize,ssize2, ssize2ij, ssize2kl, buffersize

    !!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
!    real(8)  :: twoParticlesIntegrals(1073741824)
    real(8) , allocatable :: twoParticlesIntegrals(:)
    real(8), allocatable :: twoParticlesIntegrals2(:,:)
    !!    integer(kind=8), allocatable :: indexTwoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: intIJ(:,:)
    integer(8), allocatable :: indexIJ(:,:)
    integer(8), allocatable :: nIJ(:)
    integer(8), allocatable :: totalnIJ(:)
    integer(8), allocatable :: totalnIJ2(:)

    real(8), allocatable :: intPQ(:,:)
    integer(8), allocatable :: indexPQ(:,:)
    integer(8), allocatable :: nPQ(:)
    integer(8), allocatable :: totalnPQ2(:)

    integer(8), allocatable :: lastnIJ(:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempBC(:,:)
    real(8), allocatable :: tempCBC(:,:)
    real(8), allocatable :: auxtempB(:,:)
    integer(8), allocatable :: xy(:,:)
    integer(8), allocatable :: xypair(:,:)
    integer(8), allocatable :: ioff(:)
    integer(8), allocatable :: ioff2(:)
    integer(8), allocatable :: ioff3(:)
    real(8), allocatable :: twoint(:)
    real(8), allocatable :: auxtwoint(:)
    integer(8), allocatable :: twoindex(:)
    real(8), allocatable :: tempC(:)
    real(8), allocatable :: coeffMatrix(:,:)
    integer(8), allocatable :: coeffMap(:,:)
    integer(8), allocatable :: coeffSize(:)
    integer(8) :: nonZeroIntegrals
    type(matrix) :: densityMatrix
    real(8), allocatable :: tmpArray(:,:)
    !$ real(8) :: timeA(10), timeB(10)

    integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: pp(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: qq(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer(8) :: ab, cd, ee, abcd, pq, rs, ij, kl, ij2, kl2

    real(8), allocatable :: auxIntegrals(:)
    integer(8), allocatable :: auxij(:), auxkl(:)

    integer :: p, q, r, s, mu, nu, auxnu, lambda, sigma, m, n, u, x, i, j, k, l
    integer :: topp
    integer(8), allocatable :: ijmap(:), klmap(:)
    integer(8) nstacks, stacktop
    integer(8) :: index, index2, pqrs, posxy, mm, ij_l, ij_u, kl_l, kl_u


    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: unittmp
    integer :: unittmp2
    integer(8) :: filesize
    logical :: disk
    real(8) :: coulomb
    real(8) :: exchange

    disk = .false.
    !disk = .true.

!     call system(" lowdin-ints.x TWO_PARTICLE_R12")

    if ( disk ) then

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    this%specieID = specieID

    call TransformIntegralsE_checkMOIntegralType(specieID, this)

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    ssize = this%numberOfContractions

    !!call TransformIntegralsE_setmem( ssize, twoParticlesIntegrals )

!    call matrix_constructor(densityMatrix, ssize, ssize, 0.0_8)
!    densityMatrix%values = 1
    ssize2 = (ssize * (ssize + 1))/2 
    !ssize2ab = (this%s_u - this%s_l + 1) * (this%r_u - this%r_l + 1) + ((this%s_u - this%s_l + 1)* (this%s_u - this%s_l + 2)/2)
    !ssize2ab = ssize2

    allocate (xy ( ssize , ssize ) )
    allocate (xypair ( 2 , ssize2 ) )
    allocate (ioff ( ssize2 ) )
    allocate (ioff3 ( ssize2 ) )

    xy = 0
    xypair = 0
    ioff = 0
    ioff3 = 0

    m = 0
    do p = 1,  ssize
      do q = p,  ssize
        m = m + 1
        xy(p,q) = m        
        xy(q,p) = m        
        !print *, m, p, q
        xypair(1,m) = p
        xypair(2,m) = q
      end do
    end do

    m = 0
    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
      end do
    end do

    ssize2ij = m
    !print *, ssize2ab

    allocate( ijmap (m)) 
    ijmap = 0
    m = 0

    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
        ijmap(m) = xy(i,j)
      end do
    end do

    m = 0
    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
      end do
    end do

    ssize2kl = m
    !print *, ssize2ab

    allocate( klmap (m)) 
    klmap = 0
    m = 0

    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
        klmap(m) = xy(k,l)
      end do
    end do


    ioff(1) = 0
    do pq = 2, ssize2 
      ioff(pq) = ioff(pq-1) + ssize2 - pq + 1 
    end do

    do pq = 1, ssize2  !!
      ioff3(pq) = (pq-1)*ssize2*8*2  !! two 8-bits values 

    end do

    buffersize = CONTROL_instance%IT_BUFFERSIZE
    write(*,"(T4,A36,I8)") "Buffer Size: ", buffersize

    allocate (intPQ ( buffersize, ssize2 ))  !!! -> intIJ 
    intPQ = 0

    allocate (indexPQ ( buffersize, ssize2 )) !!! indexIJ
    indexPQ = 0

    allocate (nPQ ( ssize2 ) ) !!! nIJ
    nPQ = 0

    allocate (totalnPQ2 ( ssize2 ) )  !!!
    totalnPQ2 = 0

    !!

    allocate (intIJ ( buffersize, ssize2ij ))  !!! -> intIJ 
    intIJ = 0

    allocate (indexIJ ( buffersize, ssize2ij )) !!! indexIJ
    indexIJ = 0

    allocate (nIJ ( ssize2ij ) ) !!! nIJ
    nIJ = 0

    !allocate (totalnIJ ( ssize2 ) )
    !totalnIJ = 0

    allocate (totalnIJ2 ( ssize2ij ) )  !!!
    totalnIJ2 = 0

!$  timeA(1) = omp_get_wtime()
    !! Read integrals

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    unittmp = 2000
    open( unit=unittmp,FILE=trim(nameOfSpecie)//"it.tmp", status='replace',access='stream', form='Unformatted')

    if ( trim(nameOfSpecie) == "E-BETA" ) then
       open( UNIT=unitid,FILE=trim(fileid)//trim("E-ALPHA")//".ints", status='old',access='stream', form='Unformatted')
    else 
       open( unit=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//".ints", status='old',access='stream', form='Unformatted')
    end if

    rewind(unitid)
    flush(unitid)

    !! get size
    inquire(unit=unitid, size=filesize)
    filesize= filesize/24
    filesize= filesize/CONTROL_instance%INTEGRAL_STACK_SIZE

    !allocate (twoParticlesIntegrals ( nonzeroIntegrals ) )
    !twoParticlesIntegrals = 0.0_8

    do p = 1, filesize -1

      read(UNIT=unitid, iostat=status) pp,qq,rr,ss, shellIntegrals

      do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

          pq = xy(pp(i),qq(i))
          rs = xy(rr(i),ss(i))

          nPQ(pq) = nPQ(pq) + 1
          intPQ(nPQ(pq),pq) = shellIntegrals(i)
          indexPQ(nPQ(pq),pq) = rs

          if ( nPQ(pq) == buffersize ) then

            write(unittmp,pos=ioff3(pq)+totalnPQ2(pq)*16+1) indexPQ(:,pq), intPQ(:,pq)
    
            totalnPQ2(pq) = totalnPQ2(pq) + nPQ(pq)
            indexPQ(:,pq) = 0
            intPQ(:,pq) = 0
            nPQ(pq) = 0
          end if

          nPQ(rs) = nPQ(rs) + 1
          intPQ(nPQ(rs),rs) = shellIntegrals(i)
          indexPQ(nPQ(rs),rs) = pq

          if ( nPQ(rs) == buffersize ) then

            write(unittmp,pos=ioff3(rs)+totalnPQ2(rs)*16+1) indexPQ(:,rs), intPQ(:,rs)
    
            totalnPQ2(rs) = totalnPQ2(rs) + nPQ(rs)
            indexPQ(:,rs) = 0
            intPQ(:,rs) = 0
            nPQ(rs) = 0
          end if




          !if ( pq >= rs ) then 
          !  index2 = ioff(rs) + pq
          !else 
          !  index2 = ioff(pq) + rs
          !end if
          !twoParticlesIntegrals(index2) = shellIntegrals(i)

      end do
    end do 

    !! read last stack
    read(UNIT=unitid, iostat=status) pp, qq, rr, ss, shellIntegrals

    do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
      if( pp(i) == -1 ) exit

      pq = xy(pp(i),qq(i))
      rs = xy(rr(i),ss(i))

          nPQ(pq) = nPQ(pq) + 1
          intPQ(nPQ(pq),pq) = shellIntegrals(i)
          indexPQ(nPQ(pq),pq) = rs

          if ( nPQ(pq) == buffersize ) then

            write(unittmp,pos=ioff3(pq)+totalnPQ2(pq)*16+1) indexPQ(:,pq), intPQ(:,pq)
    
            totalnPQ2(pq) = totalnPQ2(pq) + nPQ(pq)
            indexPQ(:,pq) = 0
            intPQ(:,pq) = 0
            nPQ(pq) = 0
          end if

          nPQ(rs) = nPQ(rs) + 1
          intPQ(nPQ(rs),rs) = shellIntegrals(i)
          indexPQ(nPQ(rs),rs) = pq

          if ( nPQ(rs) == buffersize ) then

            write(unittmp,pos=ioff3(rs)+totalnPQ2(rs)*16+1) indexPQ(:,rs), intPQ(:,rs)
    
            totalnPQ2(rs) = totalnPQ2(rs) + nPQ(rs)
            indexPQ(:,rs) = 0
            intPQ(:,rs) = 0
            nPQ(rs) = 0
          end if

    end do 


    do i = 1, ssize2
      if ( nPQ(i) > 0 ) then
        write(unittmp,pos=ioff3(i)+totalnPQ2(i)*16+1) indexPQ(1:nPQ(i),i), intPQ(1:nPQ(i),i)
!        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*8+ssize2*8+1) intIJ(:,i)
      end if
    end do

    totalnPQ2 = totalnPQ2 + nPQ



    close (unitid)
    close (unittmp)

    !$OMP END PARALLEL

!$  timeB(1) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Yoshimine's sorting time(s): ", timeB(1) -timeA(1) 

    !! First half-transformation
!$  timeA(2) = omp_get_wtime()

    allocate (twoindex ( buffersize ) )
    twoindex = 0

    allocate (auxtwoint ( buffersize ) )
    auxtwoint = 0

    allocate (twoint ( ssize2 ) )
    twoint = 0

    allocate (tempB ( ssize, ssize ) )
    tempB = 0

    allocate (tempBC ( ssize, ssize ) )
    tempBC = 0

    unittmp2 = 2200
    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"it2.tmp", status='replace',access='stream', form='Unformatted')

    open( unit=unittmp,FILE=trim(nameOfSpecie)//"it.tmp", status='old',access='stream', form='Unformatted')

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !! read seque
    !do p = 1, this%numberOfContractions
    !  do q = 1, p
    do pq = 1, ssize2
   
      tempB = 0

      if ( totalnPQ2(pq) == 0 ) cycle
      twoindex = 0
      twoint = 0
      auxtwoint = 0
      !!read(unittmp, pos=ioff2(pq)+8 ) twoint
      nstacks = ceiling(real(totalnPQ2(pq)/real(buffersize))) !! store this

      do n = 0, nstacks-1 -1 
        read(unittmp, pos=ioff3(pq)+n*buffersize*16+1 ) twoindex(1:buffersize), auxtwoint(1:buffersize) 
        do i = 1, buffersize
          rs = twoindex(i)
          r = xypair(1,rs)
          s = xypair(2,rs)
          tempB(s,r) = auxtwoint(i) 
          tempB(r,s) = tempB(s,r)

!          twoint(twoindex(i)) = auxtwoint(i)
        end do
      end do

      stacktop = totalnPQ2(pq) - (nstacks -1) * buffersize
      !print *, totalnIJ2(pq), nstacks, stacktop
      read(unittmp, pos=ioff3(pq)+(nstacks-1)*buffersize*16+1 ) twoindex(1:stacktop), auxtwoint(1:stacktop) 
      do i = 1, stacktop
!        twoint(twoindex(i)) = auxtwoint(i)
          rs = twoindex(i)
          r = xypair(1,rs)
          s = xypair(2,rs)
          tempB(s,r) = auxtwoint(i) 
          tempB(r,s) = tempB(s,r)

      end do

!      do rs = 1, pq
!        !index2 = ioff(rs) + pq
!        r = xypair(1,rs)
!        s = xypair(2,rs)
!
!        !tempB(s,r) = twoParticlesIntegrals(index2) 
!        tempB(s,r) = twoint(index2) 
!        tempB(r,s) = tempB(s,r)
!      end do 
!
!      do rs = pq+1, ssize2
!        index2 = ioff(pq) + rs
!        r = xypair(1,rs)
!        s = xypair(2,rs)
!
!        !tempB(s,r) = twoParticlesIntegrals(index2) 
!        tempB(s,r) = twoint(index2) 
!        tempB(r,s) = tempB(s,r)
!      end do 

      !! matmul
      !!tempB = matmul(matmul(transpose( coefficientsOfAtomicOrbitals%values ),tempB),  coefficientsOfAtomicOrbitals%values )

      !! dgemm
      !!tempBC = 0
      !!tempCBC = 0
      !!call dgemm ("T","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, coefficientsOfAtomicOrbitals%values, int(ssize,4), tempB, int(ssize,4), 0.0_8, tempBC, int(ssize,4) )
      !!call dgemm ("N","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, tempBc,  int(ssize,4),  coefficientsOfAtomicOrbitals%values,  int(ssize,4), 0.0_8, tempCBC,  int(ssize,4))
      !!tempB = tempCBC

      !! d_ij = sum_mu b_imu * c_muj = sum_mu b_imu * sum_nu a_munu b_nuj

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (j, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do j = this%q_l,  this%q_u
        do mu = 1, ssize
!          do nu = 1, ssize
!            tempBC(mu,j) = tempBC(mu,j) + tempB(nu,mu) * coefficientsOfAtomicOrbitals%values(nu,j)
!          end do

          tempBC(mu,j) = sum ( tempB(:,mu) * coefficientsOfAtomicOrbitals%values(:,j) )

        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0
        
!        do i = this%r_l, this%r_u
!          do j = this%s_l,  i!this%s_u
!            if ( j > this%s_u ) cycle
      do ij = 1, ssize2ij
      !do ij = 1, ssize2 !! -> ssize2ab
      !do ij = ij_l, ij_u
        ij2 = ijmap(ij)
        j = xypair(1,ij2)
        i = xypair(2,ij2)

!        do mu = 1, this%numberOfContractions
!          tempB(i,j) = tempB(i,j) + coefficientsOfAtomicOrbitals%values(mu,i)*tempBC(mu,j)
!        end do 

        tempB(i,j) = sum ( coefficientsOfAtomicOrbitals%values(:,i)*tempBC(:,j) )

        !rs = xy(j,i) 
        if ( abs(tempB(i,j)) > 1E-10 ) then

          nIJ(ij) = nIJ(ij) + 1
          intIJ(nIJ(ij),ij) = tempB(i,j)
          indexIJ(nIJ(ij),ij) = pq

          if ( nIJ(ij) == buffersize ) then

            write(unittmp2,pos=ioff3(ij)+totalnIJ2(ij)*16+1) indexIJ(:,ij), intIJ(:,ij)
    
            totalnIJ2(ij) = totalnIJ2(ij) + nIJ(ij)
            indexIJ(:,ij) = 0
            intIJ(:,ij) = 0
            nIJ(ij) = 0
          end if

        end if

      end do !! ij
    end do  !!pq

    do i = 1, ssize2ij
      if ( nIJ(i) > 0 ) then
        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*16+1) indexIJ(1:nIJ(i),i), intIJ(1:nIJ(i),i)
!        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*8+ssize2*8+1) intIJ(:,i)
      end if
    end do

    totalnIJ2 = totalnIJ2 + nIJ

    close (unittmp)
    !deallocate (twoParticlesIntegrals)

!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "First Half transformation time(s): ", timeB(2) -timeA(2) 

    !! Second half-transformation
!$  timeA(3) = omp_get_wtime()

    allocate (auxIntegrals(integralStackSize) )
    allocate (auxij(integralStackSize) )
    allocate (auxkl(integralStackSize) )

    auxIntegrals = 0.0_8
    auxij = 0_8
    auxkl = 0_8

    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"it2.tmp", status='old',access='stream', form='Unformatted')

    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    mm = 0

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !kl_l = (this%s_u - this%s_l + 1 ) 
    !kl_u = (this%r_u - this%r_l + 1 )*(this%s_u - this%s_l + 1 )  + (this%s_u - this%s_l + 1 )*(this%s_u - this%s_l + 2 )/2

    !do p = this%p_l, this%p_u
    !  do q = this%q_l, p
    !    if ( q > this%q_u ) cycle
    !do ij = 1, ssize2
    do ij = 1, ssize2ij
      ij2 = ijmap(ij)
    !do ij = ij_l, ij_u
      !pq = xy(p,q)
      !read(unittmp2, rec=ioff2(pq)+1) twoint(1:ssize2)

      if ( totalnIJ2(ij) == 0 ) cycle
      twoindex = 0
      twoint = 0
      auxtwoint = 0
      !!read(unittmp, pos=ioff2(pq)+8 ) twoint
      nstacks = ceiling(real(totalnIJ2(ij)/real(buffersize))) !! store this

      do n = 0, nstacks-1 -1 
        read(unittmp2, pos=ioff3(ij)+n*buffersize*16+1 ) twoindex(1:buffersize), auxtwoint(1:buffersize) 
        do i = 1, buffersize
          twoint(twoindex(i)) = auxtwoint(i)
        end do
      end do

      stacktop = totalnIJ2(ij) - (nstacks -1) * buffersize
      !print *, totalnIJ2(pq), nstacks, stacktop
      read(unittmp2, pos=ioff3(ij)+(nstacks-1)*buffersize*16+1 ) twoindex(1:stacktop), auxtwoint(1:stacktop) 
      do i = 1, stacktop
        twoint(twoindex(i)) = auxtwoint(i)
      end do

      !do i = 1, totalnIJ2(pq)
      !  twoint(twoindex(i)) = auxtwoint(i)
      !end do

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (k, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do k = this%s_l,  this%s_u
        do mu = 1, this%numberOfContractions
          do nu = 1, this%numberOfContractions
            rs = xy(mu,nu)
            tempBC(mu,k) = tempBC(mu,k) +  twoint(rs)* coefficientsOfAtomicOrbitals%values(nu,k)
          end do
        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0

      !do i = this%r_l, this%r_u
      !  do j = this%s_l,  i!this%s_u
      !    if ( j > this%s_u ) cycle
      !do kl = 1, ssize2
      do kl = 1, ssize2kl
      !do kl = kl_l, kl_u
        kl2 = klmap(kl)
        l = xypair(1,kl2)
        k = xypair(2,kl2)

        !do mu = 1, this%numberOfContractions
        !  tempB(k,l) = tempB(k,l) + coefficientsOfAtomicOrbitals%values(mu,k)*tempBC(mu,l)
        !end do 
        tempB(k,l) = sum( coefficientsOfAtomicOrbitals%values(:,k)*tempBC(:,l) )

        !rs = xy(i,j) 
        if ( abs(tempB(k,l)) > 1E-10 ) then

          m = m + 1
          auxIntegrals(m) = tempB(k,l)
          auxij(m) = ij2
          auxkl(m) = kl2

          if (m == integralStackSize ) then
            write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
            mm = mm + m
            m = 0
            auxIntegrals = 0
            auxij = 0
            auxkl = 0
          end if
        end if

      end do !! rs
    end do !! pq

    close (unittmp2)
    mm = mm + m 
    m = m + 1
    auxij(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

!$  timeB(3) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Second half transformation time(s): ", timeB(3) -timeA(3) 

    write (*,"(T4,A36,I12)") "Non-zero transformed integrals2: ", mm




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    else !! memory

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    this%specieID = specieID

    call TransformIntegralsE_checkMOIntegralType(specieID, this)

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    ssize = this%numberOfContractions

    call TransformIntegralsE_setmem( ssize, twoParticlesIntegrals )

!    call matrix_constructor(densityMatrix, ssize, ssize, 0.0_8)
!    densityMatrix%values = 1
    ssize2 = (ssize * (ssize + 1))/2 
    !ssize2ab = (this%s_u - this%s_l + 1) * (this%r_u - this%r_l + 1) + ((this%s_u - this%s_l + 1)* (this%s_u - this%s_l + 2)/2)
    !ssize2ab = ssize2

    allocate (xy ( ssize , ssize ) )
    allocate (xypair ( 2 , ssize2 ) )
    allocate (ioff ( ssize2 ) )

    xy = 0
    xypair = 0
    ioff = 0

    m = 0
    do p = 1,  ssize
      do q = p,  ssize
        m = m + 1
        xy(p,q) = m        
        xy(q,p) = m        
        !print *, m, p, q
        xypair(1,m) = p
        xypair(2,m) = q
      end do
    end do

    m = 0
    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
      end do
    end do

    ssize2ij = m
    !print *, ssize2ab

    allocate( ijmap (m)) 
    ijmap = 0
    m = 0

    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
        ijmap(m) = xy(i,j)
      end do
    end do

    m = 0
    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
      end do
    end do

    ssize2kl = m
    !print *, ssize2ab

    allocate( klmap (m)) 
    klmap = 0
    m = 0

    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
        klmap(m) = xy(k,l)
      end do
    end do


    ioff(1) = 0
    do pq = 2, ssize2 
      ioff(pq) = ioff(pq-1) + ssize2 - pq + 1 
    end do

    allocate (ioff3 ( ssize2ij ) )
    ioff3 = 0

    do pq = 1, ssize2ij  !!
      ioff3(pq) = (pq-1)*ssize2*8*2  !! two 8-bits values 
    end do

    buffersize = CONTROL_instance%IT_BUFFERSIZE
    write(*,"(T4,A36,I8)") "Buffer Size: ", buffersize

!$  timeA(1) = omp_get_wtime()
    !! Read integrals

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    unittmp = 2000
    !open( unit=unittmp,FILE=trim(nameOfSpecie)//"it.tmp", status='replace',access='stream', form='Unformatted')

    if ( trim(nameOfSpecie) == "E-BETA" ) then
       open( UNIT=unitid,FILE=trim(fileid)//trim("E-ALPHA")//".ints", status='old',access='stream', form='Unformatted')
    else 
       open( unit=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//".ints", status='old',access='stream', form='Unformatted')
    end if

    rewind(unitid)
    flush(unitid)

    !! get size
    inquire(unit=unitid, size=filesize)
    filesize= filesize/24
    filesize= filesize/CONTROL_instance%INTEGRAL_STACK_SIZE

    !allocate (twoParticlesIntegrals ( nonzeroIntegrals ) )
    !twoParticlesIntegrals = 0.0_8

    do p = 1, filesize -1

      read(UNIT=unitid, iostat=status) pp,qq,rr,ss, shellIntegrals

      do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

          pq = xy(pp(i),qq(i))
          rs = xy(rr(i),ss(i))

          if ( pq >= rs ) then 
            index2 = ioff(rs) + pq
          else 
            index2 = ioff(pq) + rs
          end if

          twoParticlesIntegrals(index2) = shellIntegrals(i)
      end do
    end do 

    !! read last stack
    read(UNIT=unitid, iostat=status) pp, qq, rr, ss, shellIntegrals

    do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
      if( pp(i) == -1 ) exit

      pq = xy(pp(i),qq(i))
      rs = xy(rr(i),ss(i))

      if ( pq >= rs ) then 
        index2 = ioff(rs) + pq
      else 
        index2 = ioff(pq) + rs
      end if

      twoParticlesIntegrals(index2) = shellIntegrals(i)

    end do 

    close (unitid)

    !$OMP END PARALLEL

!$  timeB(1) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Loading AO integrals time(s): ", timeB(1) -timeA(1) 

    !! First half-transformation
!$  timeA(2) = omp_get_wtime()

    allocate (intIJ ( buffersize, ssize2ij ))  !!! -> intIJ 
    intIJ = 0

    allocate (indexIJ ( buffersize, ssize2ij )) !!! indexIJ
    indexIJ = 0

    allocate (nIJ ( ssize2ij ) ) !!! nIJ
    nIJ = 0

    !allocate (totalnIJ ( ssize2 ) )
    !totalnIJ = 0

    allocate (totalnIJ2 ( ssize2ij ) )  !!!
    totalnIJ2 = 0

    allocate (twoindex ( buffersize ) )
    twoindex = 0

    allocate (auxtwoint ( buffersize ) )
    auxtwoint = 0

    allocate (twoint ( ssize2 ) )
    twoint = 0

    allocate (tempB ( ssize, ssize ) )
    tempB = 0

    allocate (tempBC ( ssize, ssize ) )
    tempBC = 0

    unittmp2 = 2200
    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"it2.tmp", status='replace',access='stream', form='Unformatted')

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !! read seque
    !do p = 1, this%numberOfContractions
    !  do q = 1, p
    do pq = 1, ssize2
   
      tempB = 0

      do rs = 1, pq
        index2 = ioff(rs) + pq
        r = xypair(1,rs)
        s = xypair(2,rs)

        tempB(s,r) = twoParticlesIntegrals(index2) 
        tempB(r,s) = tempB(s,r)
      end do 

      do rs = pq+1, ssize2
        index2 = ioff(pq) + rs
        r = xypair(1,rs)
        s = xypair(2,rs)

        tempB(s,r) = twoParticlesIntegrals(index2) 
        tempB(r,s) = tempB(s,r)
      end do 

      !! matmul
      !!tempB = matmul(matmul(transpose( coefficientsOfAtomicOrbitals%values ),tempB),  coefficientsOfAtomicOrbitals%values )

      !! dgemm
      !!tempBC = 0
      !!tempCBC = 0
      !!call dgemm ("T","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, coefficientsOfAtomicOrbitals%values, int(ssize,4), tempB, int(ssize,4), 0.0_8, tempBC, int(ssize,4) )
      !!call dgemm ("N","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, tempBc,  int(ssize,4),  coefficientsOfAtomicOrbitals%values,  int(ssize,4), 0.0_8, tempCBC,  int(ssize,4))
      !!tempB = tempCBC

      !! d_ij = sum_mu b_imu * c_muj = sum_mu b_imu * sum_nu a_munu b_nuj

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (j, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do j = this%q_l,  this%q_u
        do mu = 1, ssize
!          do nu = 1, ssize
!            tempBC(mu,j) = tempBC(mu,j) + tempB(nu,mu) * coefficientsOfAtomicOrbitals%values(nu,j)
!          end do

          tempBC(mu,j) = sum ( tempB(:,mu) * coefficientsOfAtomicOrbitals%values(:,j) )

        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0
        
!        do i = this%r_l, this%r_u
!          do j = this%s_l,  i!this%s_u
!            if ( j > this%s_u ) cycle
      do ij = 1, ssize2ij
      !do ij = 1, ssize2 !! -> ssize2ab
      !do ij = ij_l, ij_u
        ij2 = ijmap(ij)
        j = xypair(1,ij2)
        i = xypair(2,ij2)

!        do mu = 1, this%numberOfContractions
!          tempB(i,j) = tempB(i,j) + coefficientsOfAtomicOrbitals%values(mu,i)*tempBC(mu,j)
!        end do 

        tempB(i,j) = sum ( coefficientsOfAtomicOrbitals%values(:,i)*tempBC(:,j) )

        !rs = xy(j,i) 
        if ( abs(tempB(i,j)) > 1E-10 ) then

          nIJ(ij) = nIJ(ij) + 1
          intIJ(nIJ(ij),ij) = tempB(i,j)
          indexIJ(nIJ(ij),ij) = pq

          if ( nIJ(ij) == buffersize ) then

            write(unittmp2,pos=ioff3(ij)+totalnIJ2(ij)*16+1) indexIJ(:,ij), intIJ(:,ij)
    
            totalnIJ2(ij) = totalnIJ2(ij) + nIJ(ij)
            indexIJ(:,ij) = 0
            intIJ(:,ij) = 0
            nIJ(ij) = 0
          end if

        end if

      end do !! ij
    end do  !!pq

    do i = 1, ssize2ij
      if ( nIJ(i) > 0 ) then
        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*16+1) indexIJ(1:nIJ(i),i), intIJ(1:nIJ(i),i)
!        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*8+ssize2*8+1) intIJ(:,i)
      end if
    end do

    totalnIJ2 = totalnIJ2 + nIJ

    close (unittmp)
    deallocate (twoParticlesIntegrals)

!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "First Half transformation time(s): ", timeB(2) -timeA(2) 

    !! Second half-transformation
!$  timeA(3) = omp_get_wtime()

    allocate (auxIntegrals(integralStackSize) )
    allocate (auxij(integralStackSize) )
    allocate (auxkl(integralStackSize) )

    auxIntegrals = 0.0_8
    auxij = 0_8
    auxkl = 0_8

    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"it2.tmp", status='old',access='stream', form='Unformatted')

    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    mm = 0

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !kl_l = (this%s_u - this%s_l + 1 ) 
    !kl_u = (this%r_u - this%r_l + 1 )*(this%s_u - this%s_l + 1 )  + (this%s_u - this%s_l + 1 )*(this%s_u - this%s_l + 2 )/2

    !do p = this%p_l, this%p_u
    !  do q = this%q_l, p
    !    if ( q > this%q_u ) cycle
    !do ij = 1, ssize2
    do ij = 1, ssize2ij
      ij2 = ijmap(ij)
    !do ij = ij_l, ij_u
      !pq = xy(p,q)
      !read(unittmp2, rec=ioff2(pq)+1) twoint(1:ssize2)

      if ( totalnIJ2(ij) == 0 ) cycle
      twoindex = 0
      twoint = 0
      auxtwoint = 0
      !!read(unittmp, pos=ioff2(pq)+8 ) twoint
      nstacks = ceiling(real(totalnIJ2(ij)/real(buffersize))) !! store this

      do n = 0, nstacks-1 -1 
        read(unittmp2, pos=ioff3(ij)+n*buffersize*16+1 ) twoindex(1:buffersize), auxtwoint(1:buffersize) 
        do i = 1, buffersize
          twoint(twoindex(i)) = auxtwoint(i)
        end do
      end do

      stacktop = totalnIJ2(ij) - (nstacks -1) * buffersize
      !print *, totalnIJ2(pq), nstacks, stacktop
      read(unittmp2, pos=ioff3(ij)+(nstacks-1)*buffersize*16+1 ) twoindex(1:stacktop), auxtwoint(1:stacktop) 
      do i = 1, stacktop
        twoint(twoindex(i)) = auxtwoint(i)
      end do

      !do i = 1, totalnIJ2(pq)
      !  twoint(twoindex(i)) = auxtwoint(i)
      !end do

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (k, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do k = this%s_l,  this%s_u
        do mu = 1, this%numberOfContractions
          do nu = 1, this%numberOfContractions
            rs = xy(mu,nu)
            tempBC(mu,k) = tempBC(mu,k) +  twoint(rs)* coefficientsOfAtomicOrbitals%values(nu,k)
          end do
        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0

      !do i = this%r_l, this%r_u
      !  do j = this%s_l,  i!this%s_u
      !    if ( j > this%s_u ) cycle
      !do kl = 1, ssize2
      do kl = 1, ssize2kl
      !do kl = kl_l, kl_u
        kl2 = klmap(kl)
        l = xypair(1,kl2)
        k = xypair(2,kl2)

        !do mu = 1, this%numberOfContractions
        !  tempB(k,l) = tempB(k,l) + coefficientsOfAtomicOrbitals%values(mu,k)*tempBC(mu,l)
        !end do 
        tempB(k,l) = sum( coefficientsOfAtomicOrbitals%values(:,k)*tempBC(:,l) )

        !rs = xy(i,j) 
        if ( abs(tempB(k,l)) > 1E-10 ) then

          m = m + 1
          auxIntegrals(m) = tempB(k,l)
          auxij(m) = ij2
          auxkl(m) = kl2

          if (m == integralStackSize ) then
            write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
            mm = mm + m
            m = 0
            auxIntegrals = 0
            auxij = 0
            auxkl = 0
          end if
        end if

      end do !! rs
    end do !! pq

    close (unittmp2)
    mm = mm + m 
    m = m + 1
    auxij(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

!$  timeB(3) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Second half transformation time(s): ", timeB(3) -timeA(3) 

    write (*,"(T4,A36,I12)") "Non-zero transformed integrals: ", mm
    
    end if

  end subroutine TransformIntegralsE_atomicToMolecularOfOneSpecie

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !!    a integrales moleculares.
  !<
  !!  subroutine TransformIntegralsE_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
  !!       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
  subroutine TransformIntegralsE_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, &
       otherSpecieID, nameOfOtherSpecie )
    implicit none
    type(TransformIntegralsE) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID, otherSpecieID
    character(*) :: nameOfSpecie, nameOfOtherSpecie
    integer :: nproc
    integer :: integralStackSize

    integer :: status
    integer(8) :: ssizea,ssizeb, ssize2a, ssize2b, ssize2ij, ssize2kl, buffersize

    !!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
!    real(8)  :: twoParticlesIntegrals(1073741824)
    real(8) , allocatable :: twoParticlesIntegrals(:)
    !!    integer(kind=8), allocatable :: indexTwoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: intIJ(:,:)
    integer(8), allocatable :: indexIJ(:,:)
    integer(8), allocatable :: nIJ(:)
    integer(8), allocatable :: totalnIJ(:)
    integer(8), allocatable :: totalnIJ2(:)

    real(8), allocatable :: intPQ(:,:)
    integer(8), allocatable :: indexPQ(:,:)
    integer(8), allocatable :: nPQ(:)
    integer(8), allocatable :: totalnPQ2(:)

    integer(8), allocatable :: lastnIJ(:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempBC(:,:)
    real(8), allocatable :: tempCBC(:,:)
    real(8), allocatable :: auxtempB(:,:)
    integer(8), allocatable :: xya(:,:), xyb(:,:)
    integer(8), allocatable :: xypaira(:,:), xypairb(:,:)
    integer(8), allocatable :: ioffa(:), ioffb(:)
    integer(8), allocatable :: ioff3(:)
    real(8), allocatable :: twoint(:)
    real(8), allocatable :: auxtwoint(:)
    integer(8), allocatable :: twoindex(:)
    real(8), allocatable :: tempC(:)
    integer(8) :: nonZeroIntegrals
    real(8), allocatable :: tmpArray(:,:)
    !$ real(8) :: timeA(10), timeB(10)

    integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: pp(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: qq(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer(8) :: ab, cd, ee, abcd, pq, rs, ij, kl, ij2, kl2

    real(8), allocatable :: auxIntegrals(:)
    integer(8), allocatable :: auxij(:), auxkl(:)

    integer :: p, q, r, s, mu, nu, auxnu, lambda, sigma, m, n, u, x, i, j, k, l
    integer :: topp
    integer(8), allocatable :: ijmap(:), klmap(:)
    integer(8) nstacks, stacktop
    integer(8) :: index, index2, pqrs, posxy, mm, ij_l, ij_u, kl_l, kl_u

    !! OpenMP related variables
    integer :: otherSsize
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: unittmp
    integer :: unittmp2
    integer(8) :: filesize

    ! Reads the number of cores

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    ! call TransformIntegralsE_getNumberOfNonZeroCouplingIntegrals( specieID, otherSpecieID, nproc, nonZeroIntegrals )

    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"

    call TransformIntegralsE_checkInterMOIntegralType(specieID, otherSpecieID, this)

    !if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    !allocate (twoParticlesIntegrals ( nonZeroIntegrals ) )
    !twoParticlesIntegrals = 0

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%otherNumberOfContractions=size(otherCoefficientsOfAtomicOrbitals%values,dim=1)
    this%specieID = specieID

    ssizea = this%numberOfContractions
    ssizeb = this%otherNumberOfContractions

    call TransformIntegralsE_setmemInter( ssizea, ssizeb, twoParticlesIntegrals )

    ssize2a = (ssizea * (ssizea + 1))/2 
    ssize2b = (ssizeb * (ssizeb + 1))/2 

    allocate (xya ( ssizea , ssizea ) )
    allocate (xyb ( ssizeb , ssizeb ) )
    allocate (xypaira ( 2 , ssize2a ) )
    allocate (xypairb ( 2 , ssize2b ) )
    allocate (ioffa ( ssize2a ) )
    allocate (ioffb ( ssize2b ) )

    xya = 0
    xyb = 0
    xypaira = 0
    xypairb = 0
    ioffa = 0
    ioffb = 0

    m = 0
    do p = 1,  ssizea
      do q = p,  ssizea
        m = m + 1
        xya(p,q) = m        
        xya(q,p) = m        
        xypaira(1,m) = p
        xypaira(2,m) = q
      end do
    end do

    m = 0
    do p = 1,  ssizeb
      do q = p,  ssizeb
        m = m + 1
        xyb(p,q) = m        
        xyb(q,p) = m        
        xypairb(1,m) = p
        xypairb(2,m) = q
      end do
    end do

    m = 0
    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
      end do
    end do

    ssize2ij = m

    allocate( ijmap (m)) 
    ijmap = 0
    m = 0

    do i = this%p_l, this%p_u
      do j = this%q_l, i
        if ( j > this%q_u) cycle
        m = m + 1
        ijmap(m) = xya(i,j)
      end do
    end do

    m = 0
    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
      end do
    end do

    ssize2kl = m

    allocate( klmap (m)) 
    klmap = 0
    m = 0

    do k = this%r_l, this%r_u
      do l = this%s_l, k
        if ( l > this%s_u) cycle
        m = m + 1
        klmap(m) = xyb(k,l)
      end do
    end do


    ioffa(1) = 0
    do pq = 2, ssize2a
      ioffa(pq) = ioffa(pq-1) + ssize2a - pq + 1 
    end do

    ioffb(1) = 0
    do pq = 2, ssize2b 
      ioffb(pq) = ioffb(pq-1) + ssize2b - pq + 1 
    end do

    allocate (ioff3 ( ssize2a ) )
    ioff3 = 0
    do pq = 1, ssize2a  !!
!    do pq = 1, ssize2ij
      ioff3(pq) = (pq-1)*ssize2b*8*2  !! two 8-bits values 
    end do

    buffersize = CONTROL_instance%IT_BUFFERSIZE
    write(*,"(T4,A36,I8)") "Buffer Size: ", buffersize

!$  timeA(1) = omp_get_wtime()
 

    !! Read integrals

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))


    unittmp = 2000
    !! open file for integrals
    open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
         STATUS='OLD', ACCESS='stream', FORM='Unformatted')

    !! get size
    inquire(unit=unitid, size=filesize)
    filesize= filesize/24
    filesize= filesize/CONTROL_instance%INTEGRAL_STACK_SIZE

    !allocate (twoParticlesIntegrals ( nonzeroIntegrals ) )
    !twoParticlesIntegrals = 0.0_8

    do p = 1, filesize -1

      read(UNIT=unitid, iostat=status) pp,qq,rr,ss, shellIntegrals

      do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

          pq = xya(pp(i),qq(i))
          rs = xyb(rr(i),ss(i))

          !if ( pq >= rs ) then 
          !  index2 = ioff(rs) + pq
          !else 
          !  index2 = ioff(pq) + rs
          !end if

          index2 = (rs-1)*ssize2a + pq
          twoParticlesIntegrals(index2) = shellIntegrals(i)

      end do
    end do 

    !! read last stack
    read(UNIT=unitid, iostat=status) pp, qq, rr, ss, shellIntegrals

    do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
      if( pp(i) == -1 ) exit

      pq = xya(pp(i),qq(i))
      rs = xyb(rr(i),ss(i))

      !if ( pq >= rs ) then 
      !  index2 = ioff(rs) + pq
      !else 
      !  index2 = ioff(pq) + rs
      !end if
      index2 = (rs-1)*ssize2a + pq
      twoParticlesIntegrals(index2) = shellIntegrals(i)

    end do 

    close (unitid)

    !$OMP END PARALLEL

!$  timeB(1) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Loading AO integrals time(s): ", timeB(1) -timeA(1) 

    !! First half-transformation
!$  timeA(2) = omp_get_wtime()

    allocate (intIJ ( buffersize, ssize2ij ))  !!! -> intIJ 
    intIJ = 0

    allocate (indexIJ ( buffersize, ssize2ij )) !!! indexIJ
    indexIJ = 0

    allocate (nIJ ( ssize2ij ) ) !!! nIJ
    nIJ = 0

    !allocate (totalnIJ ( ssize2 ) )
    !totalnIJ = 0

    allocate (totalnIJ2 ( ssize2ij ) )  !!!
    totalnIJ2 = 0

    allocate (twoindex ( buffersize ) )
    twoindex = 0

    allocate (auxtwoint ( buffersize ) )
    auxtwoint = 0

    allocate (twoint ( ssize2b ) ) !!!!!!
    twoint = 0

    allocate (tempB ( ssizea, ssizea ) )
    tempB = 0

    allocate (tempBC ( ssizea, ssizea ) )
    tempBC = 0

    unittmp2 = 2200
    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".it2.tmp", status='replace',access='stream', form='Unformatted')

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !! read seque
    !do p = 1, this%numberOfContractions
    !  do q = 1, p
    do pq = 1, ssize2b
   
      tempB = 0

      do rs = 1, ssize2a
        index2 = (pq-1)*ssize2a + rs
        r = xypaira(1,rs)
        s = xypaira(2,rs)

        tempB(s,r) = twoParticlesIntegrals(index2) 
        tempB(r,s) = tempB(s,r)
      end do 

      !! matmul
      !!tempB = matmul(matmul(transpose( coefficientsOfAtomicOrbitals%values ),tempB),  coefficientsOfAtomicOrbitals%values )

      !! dgemm
      !!tempBC = 0
      !!tempCBC = 0
      !!call dgemm ("T","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, coefficientsOfAtomicOrbitals%values, int(ssize,4), tempB, int(ssize,4), 0.0_8, tempBC, int(ssize,4) )
      !!call dgemm ("N","N", int(ssize,4), int(ssize,4), int(ssize,4), 1.0_8, tempBc,  int(ssize,4),  coefficientsOfAtomicOrbitals%values,  int(ssize,4), 0.0_8, tempCBC,  int(ssize,4))
      !!tempB = tempCBC

      !! d_ij = sum_mu b_imu * c_muj = sum_mu b_imu * sum_nu a_munu b_nuj

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (j, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do j = this%q_l,  this%q_u
        do mu = 1, ssizea
!          do nu = 1, ssize
!            tempBC(mu,j) = tempBC(mu,j) + tempB(nu,mu) * coefficientsOfAtomicOrbitals%values(nu,j)
!          end do

          tempBC(mu,j) = sum ( tempB(:,mu) * coefficientsOfAtomicOrbitals%values(:,j) )

        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0
        
      do ij = 1, ssize2ij
      !do ij = 1, ssize2 !! -> ssize2ab
      !do ij = ij_l, ij_u
        ij2 = ijmap(ij)
        j = xypaira(1,ij2)
        i = xypaira(2,ij2)

!        do mu = 1, this%numberOfContractions
!          tempB(i,j) = tempB(i,j) + coefficientsOfAtomicOrbitals%values(mu,i)*tempBC(mu,j)
!        end do 

        tempB(i,j) = sum ( coefficientsOfAtomicOrbitals%values(:,i)*tempBC(:,j) )

        !rs = xy(j,i) 
        if ( abs(tempB(i,j)) > 1E-10 ) then

          nIJ(ij) = nIJ(ij) + 1
          intIJ(nIJ(ij),ij) = tempB(i,j)
          indexIJ(nIJ(ij),ij) = pq

          if ( nIJ(ij) == buffersize ) then
            write(unittmp2,pos=ioff3(ij)+totalnIJ2(ij)*16+1) indexIJ(:,ij), intIJ(:,ij)
    
            totalnIJ2(ij) = totalnIJ2(ij) + nIJ(ij)
            indexIJ(:,ij) = 0
            intIJ(:,ij) = 0
            nIJ(ij) = 0
          end if

        end if

      end do !! ij
    end do  !!pq

    do i = 1, ssize2ij
      if ( nIJ(i) > 0 ) then
        write(unittmp2,pos=ioff3(i)+totalnIJ2(i)*16+1) indexIJ(1:nIJ(i),i), intIJ(1:nIJ(i),i)
      end if
    end do

    totalnIJ2 = totalnIJ2 + nIJ

    close (unittmp)
    deallocate (twoParticlesIntegrals)

!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "First Half transformation time(s): ", timeB(2) -timeA(2) 

    !! Second half-transformation
!$  timeA(3) = omp_get_wtime()

    deallocate (tempB )
    deallocate (tempBC )
    allocate (tempB ( ssizeb, ssizeb ) )
    tempB = 0

    allocate (tempBC ( ssizeb, ssizeb ) )
    tempBC = 0



    allocate (auxIntegrals(integralStackSize) )
    allocate (auxij(integralStackSize) )
    allocate (auxkl(integralStackSize) )

    auxIntegrals = 0.0_8
    auxij = 0_8
    auxkl = 0_8

    open( unit=unittmp2,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".it2.tmp", status='old',access='stream', form='Unformatted')

    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    mm = 0

    !ij_l = (this%q_u - this%q_l + 1 ) 
    !ij_u = (this%p_u - this%p_l + 1 )*(this%q_u - this%q_l + 1 )  + (this%q_u - this%q_l + 1 )*(this%q_u - this%q_l + 2 )/2

    !kl_l = (this%s_u - this%s_l + 1 ) 
    !kl_u = (this%r_u - this%r_l + 1 )*(this%s_u - this%s_l + 1 )  + (this%s_u - this%s_l + 1 )*(this%s_u - this%s_l + 2 )/2

    !do p = this%p_l, this%p_u
    !  do q = this%q_l, p
    !    if ( q > this%q_u ) cycle
    !do ij = 1, ssize2
    do ij = 1, ssize2ij
      ij2 = ijmap(ij)

      if ( totalnIJ2(ij) == 0 ) cycle
      twoindex = 0
      twoint = 0
      auxtwoint = 0
      nstacks = ceiling(real(totalnIJ2(ij)/real(buffersize))) !! store this

      do n = 0, nstacks-1 -1 
        read(unittmp2, pos=ioff3(ij)+n*buffersize*16+1 ) twoindex(1:buffersize), auxtwoint(1:buffersize) 
        do i = 1, buffersize
          twoint(twoindex(i)) = auxtwoint(i)
        end do
      end do

      stacktop = totalnIJ2(ij) - (nstacks -1) * buffersize
      read(unittmp2, pos=ioff3(ij)+(nstacks-1)*buffersize*16+1 ) twoindex(1:stacktop), auxtwoint(1:stacktop) 
      do i = 1, stacktop
        twoint(twoindex(i)) = auxtwoint(i)
      end do

      tempBC = 0
      !!$omp parallel & 
      !!$omp& private (k, mu, nu)  shared (ssize, tempBC, coefficientsOfAtomicOrbitals)
      !!$omp do schedule (dynamic)
      do k = this%s_l,  this%s_u
        do mu = 1, ssizeb
          do nu = 1, ssizeb
            rs = xyb(mu,nu)
            tempBC(mu,k) = tempBC(mu,k) +  twoint(rs)* otherCoefficientsOfAtomicOrbitals%values(nu,k)
          end do
        end do
      end do
      !!$omp end do 
      !!$omp end parallel

      tempB = 0

      !do i = this%r_l, this%r_u
      !  do j = this%s_l,  i!this%s_u
      !    if ( j > this%s_u ) cycle
      !do kl = 1, ssize2
      do kl = 1, ssize2kl
      !do kl = kl_l, kl_u
        kl2 = klmap(kl)
        l = xypairb(1,kl2)
        k = xypairb(2,kl2)

        !do mu = 1, this%numberOfContractions
        !  tempB(k,l) = tempB(k,l) + coefficientsOfAtomicOrbitals%values(mu,k)*tempBC(mu,l)
        !end do 
        tempB(k,l) = sum( otherCoefficientsOfAtomicOrbitals%values(:,k)*tempBC(:,l) )

        if ( abs(tempB(k,l)) > 1E-10 ) then

          m = m + 1
          auxIntegrals(m) = tempB(k,l)
          auxij(m) = ij2
          auxkl(m) = kl2

          if (m == integralStackSize ) then
            write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
            mm = mm + m
            m = 0
            auxIntegrals = 0
            auxij = 0
            auxkl = 0
          end if
        end if

      end do !! rs
    end do !! pq

    close (unittmp2)
    mm = mm + m 
    m = m + 1
    auxij(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) auxij, auxkl, auxIntegrals
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

!$  timeB(3) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Second half transformation time(s): ", timeB(3) -timeA(3) 

    write (*,"(T4,A36,I12)") "Non-zero transformed integrals: ", mm
    
  end subroutine TransformIntegralsE_atomicToMolecularOfTwoSpecies



  subroutine TransformIntegralsE_setSizeOfInterIntegralsArray ( numberOfContractions, otherNumberOfContractions, otherSsize8, &
       twoParticlesIntegrals)
    implicit none 
    integer :: numberOfContractions, otherNumberOfContractions
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer(kind=8) :: ssize8
    integer :: otherSsize8 !! Beyond 360 cartesian funtions

    ssize8 = numberOfContractions
    ssize8 = (ssize8 * (ssize8 + 1))/2 

    otherSsize8 = otherNumberOfContractions
    otherSsize8 = (otherSsize8 * (otherSsize8 + 1))/2 

    ssize8 = ( ssize8 * otherSsize8 )

    if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    allocate (twoParticlesIntegrals ( ssize8 ) )

    twoParticlesIntegrals = 0

  end subroutine TransformIntegralsE_setSizeOfInterIntegralsArray


  subroutine TransformIntegralsE_buildArrayAInter( integralArray, i, ssize , otherSsize, auxOtherSsize, auxtempA)
    implicit none 
    integer, intent(in) :: ssize, otherSsize, auxOtherSsize
    integer, intent(in) :: i
    real(8), intent(in) :: integralArray(:)
    real(8) :: auxtempA(ssize,otherSsize,otherSsize)
    integer :: j,k,l,ij, kl


    !$OMP PARALLEL DO private(j,k,l,ij,kl) shared(ssize,otherSsize,i,auxOtherSsize,integralArray)
    do j = 1, ssize
       ij = int(IndexMap_tensorR2ToVectorB( int(i,8), int(j,8), int(ssize,8)), 4)
       kl = 0
       do k = 1, otherSsize
          do l = k, otherSsize
             kl = kl + 1
             !         auxm = indexArray(auxIndex)
             auxtempA(j,k,l) = integralArray( ( int(ij,8) - 1 ) * int(auxOtherSsize,8) + int(kl,8) )
          end do
       end do
    end do
    !$OMP END PARALLEL DO  

    do l = 1, otherSsize
       do k = l+1, otherSsize
          auxtempA(:,k,l) = auxtempA(:,l,k)
       end do
    end do

  end subroutine TransformIntegralsE_buildArrayAInter

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsE_checkMOIntegralType(speciesID, this)
    implicit none
    integer :: speciesID
    type(TransformIntegralsE) :: this
    integer :: totalOccupation 
    integer :: coreOrbitals
    integer :: totalActiveOrbitals

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )

    !! Take the active space from input
    coreOrbitals=0
    if ( InputCI_Instance(speciesID)%coreOrbitals /= 0 ) coreOrbitals=InputCI_Instance(speciesID)%coreOrbitals 

    totalActiveOrbitals=MolecularSystem_getTotalNumberOfContractions( speciesID )
    if ( InputCI_Instance(speciesID)%activeOrbitals /= 0 ) totalActiveOrbitals=InputCI_Instance(speciesID)%activeOrbitals 

    !! Boundary orbitals. Default
    this%p_l = coreOrbitals+1
    this%p_u = totalActiveOrbitals
    this%q_l = coreOrbitals+1
    this%q_u = totalActiveOrbitals
    this%r_l = coreOrbitals+1
    this%r_u = totalActiveOrbitals
    this%s_l = coreOrbitals+1
    this%s_u = totalActiveOrbitals
    
    if ( trim(this%partialTransform)=="ALLACTIVE") then
       this%p_l = 1
       this%p_u = totalActiveOrbitals
       this%q_l = 1
       this%q_u = totalActiveOrbitals
       this%r_l = 1
       this%r_u = totalActiveOrbitals
       this%s_l = 1
       this%s_u = totalActiveOrbitals !this%numberOfContractions
    end if

    !! only the (ia|jb) integrals will be transformed
    if ( trim(this%partialTransform)=="MP2"  ) then

       this%p_l = totalOccupation + 1
       this%p_u = totalActiveOrbitals
       this%q_l = coreOrbitals+1
       this%q_u = totalOccupation
       this%r_l = totalOccupation + 1
       this%r_u = totalActiveOrbitals
       this%s_l = coreOrbitals+1
       this%s_u = totalOccupation

    end if

    !! only the (ip|aq) integrals will be transformed
    if ( trim(this%partialTransform)=="PT2"  ) then
    
      if ( CONTROL_instance%IONIZE_MO == 0 ) then
          !! all
          this%q_l = coreOrbitals+1!totalOccupation     !! HOMO 
          this%q_u = totalOccupation+1 !! LUMO
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation 
          this%r_l = totalOccupation + 1
          this%r_u = totalActiveOrbitals

          !! half symmetric
          !this%p_l = 1 
          !this%p_u = totalOccupation + 1
          !this%q_l = 1
          !this%q_u = totalNumberOfContractions

          !this%r_l = 1
          !this%r_u = totalOccupation 
          !this%s_l = totalOccupation + 1
          !this%s_u = totalNumberOfContractions

          !! fully symmetric
          !this%q_l = 1 
          !this%q_u = totaloccupation  
          !this%p_l = 1
          !this%p_u = totalNumberOfContractions

          !this%s_l = 1
          !this%s_u = totalNumberOfContractions
          !this%r_l =  1
          !this%r_u = totalNumberOfContractions

      else

        if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
          this%q_l = CONTROL_instance%IONIZE_MO  
          this%q_u = CONTROL_instance%IONIZE_MO 
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation !totalNumberOfContractions
          !this%s_u = totalOccupation !totalNumberOfContractions
          this%r_l = coreOrbitals+1
          this%r_u = totalActiveOrbitals

        else 

          this%q_l = CONTROL_instance%IONIZE_MO  
          this%q_u = CONTROL_instance%IONIZE_MO 
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation 
          this%r_l = totalOccupation + 1
          this%r_u = totalActiveOrbitals

        end if
      end if

    end if

    !for a simultaneous PT2 - MP2 calculation
    if ( trim(this%partialTransform)=="MP2-PT2" ) then
    
      if ( CONTROL_instance%IONIZE_MO == 0 ) then
          !! all
          this%q_l = coreOrbitals+1!totalOccupation     !! HOMO 
          this%q_u = totalOccupation + 1 !! LUMO
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation 
          this%r_l = totalOccupation + 1
          this%r_u = totalActiveOrbitals

      else

        if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
          this%q_l = coreOrbitals+1
          this%q_u = max(CONTROL_instance%IONIZE_MO,totalOccupation)
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation !totalNumberOfContractions
          !this%s_u = totalOccupation !totalNumberOfContractions
          this%r_l = coreOrbitals+1
          this%r_u = totalActiveOrbitals

        else 

          this%q_l = coreOrbitals+1
          this%q_u = max(CONTROL_instance%IONIZE_MO,totalOccupation)
          this%p_l = coreOrbitals+1
          this%p_u = totalActiveOrbitals

          this%s_l = coreOrbitals+1
          this%s_u = totalOccupation 
          this%r_l = totalOccupation + 1
          this%r_u = totalActiveOrbitals

        end if
      end if
   end if
   
   write(*,"(T15,A)") "Transformation boundaries "
   write(*,"(T15,A10,A6,A6)") "orbital","lower", "upper"
   write(*,"(T20,A5,I6,I6)") "p", this%p_l, this%p_u
   write(*,"(T20,A5,I6,I6)") "q", this%q_l, this%q_u
   write(*,"(T20,A5,I6,I6)") "r", this%r_l, this%r_u
   write(*,"(T20A5,I6,I6)") "s", this%s_l, this%s_u
   print *, ""

   
 end subroutine TransformIntegralsE_checkMOIntegralType


  subroutine TransformIntegralsE_checkInterMOIntegralType(speciesID, otherSpeciesID, this)
    implicit none
    integer :: speciesID, otherSpeciesID
    type(TransformIntegralsE) :: this
    integer :: totalOccupation, otherTotalOccupation
    integer :: coreOrbitals, otherCoreOrbitals
    integer :: totalActiveOrbitals, otherTotalActiveOrbitals
    logical :: ionizeA, ionizeB
    integer :: s
    character(10) :: nameOfSpecies
    character(10) :: nameOfOtherSpecies

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )
    otherTotalOccupation = MolecularSystem_getOcupationNumber( otherSpeciesID )

    !! Take the active space from input
    coreOrbitals=0
    otherCoreOrbitals=0
    if ( InputCI_Instance(speciesID)%coreOrbitals /= 0 ) coreOrbitals=InputCI_Instance(speciesID)%coreOrbitals
    if ( InputCI_Instance(otherSpeciesID)%coreOrbitals /= 0 ) otherCoreOrbitals=InputCI_Instance(otherSpeciesID)%coreOrbitals

    totalActiveOrbitals=MolecularSystem_getTotalNumberOfContractions( speciesID )
    otherTotalActiveOrbitals=MolecularSystem_getTotalNumberOfContractions( otherSpeciesID )
    if ( InputCI_Instance(speciesID)%activeOrbitals /= 0 ) totalActiveOrbitals=InputCI_Instance(speciesID)%activeOrbitals 
    if ( InputCI_Instance(otherSpeciesID)%activeOrbitals /= 0 ) otherTotalActiveOrbitals=InputCI_Instance(otherSpeciesID)%activeOrbitals 

    !! Boundary orbitals. Default
    this%p_l = coreOrbitals+1
    this%p_u = totalActiveOrbitals
    this%q_l = coreOrbitals+1
    this%q_u = totalActiveOrbitals
    this%r_l = otherCoreOrbitals+1
    this%r_u = otherTotalActiveOrbitals
    this%s_l = otherCoreOrbitals+1
    this%s_u = otherTotalActiveOrbitals

    if ( trim(this%partialTransform) .eq. "ALLACTIVE") then
       this%p_l = 1
       this%p_u = totalActiveOrbitals!this%numberOfContractions
       this%q_l = 1
       this%q_u = totalActiveOrbitals!this%numberOfContractions
       this%r_l = 1
       this%r_u = otherTotalActiveOrbitals!this%otherNumberOfContractions
       this%s_l = 1
       this%s_u = otherTotalActiveOrbitals!this%otherNumberOfContractions
    end if

    !! only the (ia|jb) integrals will be transformed
    if ( trim(this%partialTransform) .eq. "MP2"  ) then
       this%p_l = totalOccupation + 1
       this%p_u = totalActiveOrbitals
       this%q_l = coreOrbitals+1
       this%q_u = totalOccupation
       this%r_l = otherTotalOccupation + 1
       this%r_u = otherTotalActiveOrbitals
       this%s_l = otherCoreOrbitals+1
       this%s_u = otherTotalOccupation
    end if

    !! only the (ip|IP) integrals will be transformed.
    if ( trim(this%partialTransform) .eq. "PT2" ) then

       this%q_l = coreOrbitals+1
       this%q_u = totalOccupation + 1 !! occ + lumo (default)
       this%p_l = coreOrbitals+1
       this%p_u = totalActiveOrbitals
       this%s_l = coreOrbitals+1
       this%s_u = otherTotalOccupation + 1 !!occ + lumo
       this%r_l = coreOrbitals+1
       this%r_u = otherTotalActiveOrbitals

      if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE" ) then

        ionizeA = .false.
        ionizeB = .false.

         nameOfSpecies= trim(  MolecularSystem_getNameOfSpecie( speciesID ) )
         nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( otherSpeciesID ) )

         do s = 1, size(CONTROL_instance%IONIZE_SPECIE )
           if ( nameOfSpecies == trim(CONTROL_instance%IONIZE_SPECIE(s)) ) then
             ionizeA = .true. 
           end if
           if ( nameOfOtherSpecies == trim(CONTROL_instance%IONIZE_SPECIE(s)) ) then
             ionizeB = .true. 
           end if
         end do

        if ( CONTROL_instance%IONIZE_MO == 0 ) then

          if ( ionizeA .and. ionizeB ) then

            this%q_l = coreOrbitals+1
            this%q_u = totalOccupation + 1 !! occ + lumo (default)
            this%p_l = coreOrbitals+1
            this%p_u = totalActiveOrbitals
            this%s_l = otherCoreOrbitals+1
            this%s_u = otherTotalOccupation + 1 !!occ + lumo
            this%r_l = otherCoreOrbitals+1
            this%r_u = otherTotalActiveOrbitals

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%q_l = coreOrbitals+1
            this%q_u = totalOccupation + 1 !! occ + lumo (default)
            this%p_l = coreOrbitals+1
            this%p_u = totalActiveOrbitals
  
            this%s_l = otherCoreOrbitals+1
            this%s_u = othertotaloccupation 
            this%r_l = othertotaloccupation + 1
            this%r_u = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%q_l = coreOrbitals+1
            this%q_u = totalOccupation
            this%p_l = totalOccupation + 1
            this%p_u = totalActiveOrbitals
            this%s_l = otherCoreOrbitals+1
            this%s_u = otherTotalOccupation + 1
            this%r_l = otherCoreOrbitals+1
            this%r_u = otherTotalActiveOrbitals

          end if

        else if ( CONTROL_instance%IONIZE_MO /= 0 ) then !!occ and vir..

          if ( ionizeA .and. ionizeB ) then
            if ( CONTROL_instance%IONIZE_MO <= totalOccupation .and. CONTROL_instance%IONIZE_MO <= othertotalOccupation ) then
              this%q_l = coreOrbitals+1
              this%q_u = totalOccupation
              this%p_l = coreOrbitals+1
              this%p_u = totalActiveOrbitals
              this%s_l = otherCoreOrbitals+1
              this%s_u = otherTotalOccupation
              this%r_l = otherCoreOrbitals+1
              this%r_u = otherTotalActiveOrbitals

            else if ( CONTROL_instance%IONIZE_MO > totalOccupation .and. CONTROL_instance%IONIZE_MO > othertotalOccupation ) then

              this%q_l = coreOrbitals+1
              this%q_u = totalActiveOrbitals!CONTROL_instance%IONIZE_MO 
              this%p_l = coreOrbitals+1
              this%p_u = totalActiveOrbitals
              this%s_l = otherCoreOrbitals+1
              this%s_u = otherTotalActiveOrbitals!CONTROL_instance%IONIZE_MO 
              this%r_l = otherCoreOrbitals+1
              this%r_u = otherTotalActiveOrbitals

            end if

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%q_l = CONTROL_instance%IONIZE_MO
            this%q_u = CONTROL_instance%IONIZE_MO


        if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
            this%q_l = coreOrbitals+1
            this%q_u = totalActiveOrbitals
        end if
            this%p_l = coreOrbitals+1
            this%p_u = totalActiveOrbitals
  
            this%s_l = otherCoreOrbitals+1
            this%s_u = othertotaloccupation 
            this%r_l = othertotaloccupation + 1
            this%r_u = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%q_l = coreOrbitals+1
            this%q_u = totalOccupation
            this%p_l = totalOccupation + 1
            this%p_u = totalActiveOrbitals

            !this%s_l = CONTROL_instance%IONIZE_MO !...
            !this%s_u = CONTROL_instance%IONIZE_MO !...
            this%s_l = otherCoreOrbitals+1
            this%s_u = otherTotalActiveOrbitals

            this%r_l = otherCoreOrbitals+1
            this%r_u = otherTotalActiveOrbitals

          end if



        end if

          !if ( CONTROL_instance%IONIZE_MO /= 0 ) then
      end if

   end if

   !for a simultaneous PT2 - MP2 calculation
   if ( trim(this%partialTransform)=="MP2-PT2" ) then


      this%q_l = coreOrbitals+1
      this%q_u = totalOccupation + 1 !! occ + lumo (default)
      this%p_l = coreOrbitals+1
      this%p_u = totalActiveOrbitals
      this%s_l = otherCoreOrbitals+1
      this%s_u = otherTotalOccupation + 1 !!occ + lumo
      this%r_l = otherCoreOrbitals+1
      this%r_u = totalActiveOrbitals

      if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE" ) then

         ionizeA = .false.
         ionizeB = .false.

         nameOfSpecies= trim(  MolecularSystem_getNameOfSpecie( speciesID ) )
         nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( otherSpeciesID ) )

         do s = 1, size(CONTROL_instance%IONIZE_SPECIE )
            if ( nameOfSpecies == trim(CONTROL_instance%IONIZE_SPECIE(s)) ) then
               ionizeA = .true. 
            end if
            if ( nameOfOtherSpecies == trim(CONTROL_instance%IONIZE_SPECIE(s)) ) then
               ionizeB = .true. 
            end if
         end do

         if ( CONTROL_instance%IONIZE_MO == 0 ) then

            if ( ionizeA .and. ionizeB ) then

               this%q_l = coreOrbitals+1
               this%q_u = totalOccupation + 1 !! occ + lumo (default)
               this%p_l = coreOrbitals+1
               this%p_u = totalActiveOrbitals
               this%s_l = otherCoreOrbitals+1
               this%s_u = otherTotalOccupation + 1 !!occ + lumo
               this%r_l = otherCoreOrbitals+1
               this%r_u = otherTotalActiveOrbitals

            else if ( ionizeA .and. .not. ionizeB ) then

               this%q_l = coreOrbitals+1
               this%q_u = totalOccupation + 1 !! occ + lumo (default)
               this%p_l = coreOrbitals+1
               this%p_u = totalActiveOrbitals

               this%s_l = otherCoreOrbitals+1
               this%s_u = othertotaloccupation 
               this%r_l = othertotaloccupation + 1
               this%r_u = otherTotalActiveOrbitals

            else if ( .not. ionizeA .and. ionizeB ) then

               this%q_l = coreOrbitals+1
               this%q_u = totalOccupation
               this%p_l = totalOccupation + 1
               this%p_u = totalActiveOrbitals
               this%s_l = otherCoreOrbitals+1
               this%s_u = otherTotalOccupation + 1
               this%r_l = otherCoreOrbitals+1
               this%r_u = otherTotalActiveOrbitals

            end if

         else if ( CONTROL_instance%IONIZE_MO /= 0 ) then !!occ and vir..

            if ( ionizeA .and. ionizeB ) then
               if ( CONTROL_instance%IONIZE_MO <= totalOccupation .and. CONTROL_instance%IONIZE_MO <= othertotalOccupation ) then
                  this%q_l = coreOrbitals+1
                  this%q_u = totalOccupation
                  this%p_l = coreOrbitals+1
                  this%p_u = totalActiveOrbitals
                  this%s_l = otherCoreOrbitals+1
                  this%s_u = otherTotalOccupation
                  this%r_l = otherCoreOrbitals+1
                  this%r_u = otherTotalActiveOrbitals

               else if ( CONTROL_instance%IONIZE_MO > totalOccupation .and. CONTROL_instance%IONIZE_MO > othertotalOccupation ) then

                  this%q_l = coreOrbitals+1
                  this%q_u = totalActiveOrbitals!CONTROL_instance%IONIZE_MO 
                  this%p_l = coreOrbitals+1
                  this%p_u = totalActiveOrbitals
                  this%s_l = otherCoreOrbitals+1
                  this%s_u = otherTotalActiveOrbitals!CONTROL_instance%IONIZE_MO 
                  this%r_l = otherCoreOrbitals+1
                  this%r_u = otherTotalActiveOrbitals

               end if

            else if ( ionizeA .and. .not. ionizeB ) then

               this%q_l = coreOrbitals+1
               this%q_u = max(CONTROL_instance%IONIZE_MO,totalOccupation)

               if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
                  this%q_l = coreOrbitals+1
                  this%q_u = totalActiveOrbitals
               end if
               this%p_l = coreOrbitals+1
               this%p_u = totalActiveOrbitals

               this%s_l = otherCoreOrbitals+1
               this%s_u = othertotaloccupation 
               this%r_l = othertotaloccupation + 1
               this%r_u = otherTotalActiveOrbitals

            else if ( .not. ionizeA .and. ionizeB ) then

               this%q_l = coreOrbitals+1
               this%q_u = totalOccupation
               this%p_l = totalOccupation + 1
               this%p_u = totalActiveOrbitals

               !this%s_l = CONTROL_instance%IONIZE_MO !...
               !this%s_u = CONTROL_instance%IONIZE_MO !...
               this%s_l = otherCoreOrbitals+1
               this%s_u = otherTotalActiveOrbitals

               this%r_l = otherCoreOrbitals+1
               this%r_u = otherTotalActiveOrbitals

            end if



         end if
         
         !if ( CONTROL_instance%IONIZE_MO /= 0 ) then
      end if


   end if
   
   write(*,"(T15,A)") "Transformation boundaries "
   write(*,"(T15,A10,A6,A6)") "orbital","lower", "upper"
   write(*,"(T20,A5,I6,I6)") "p", this%p_l, this%p_u
   write(*,"(T20,A5,I6,I6)") "q", this%q_l, this%q_u
   write(*,"(T20,A5,I6,I6)") "r", this%r_l, this%r_u
   write(*,"(T20,A5,I6,I6)") "s", this%s_l, this%s_u
   print *, ""

 end subroutine TransformIntegralsE_checkInterMOIntegralType

  subroutine TransformIntegralsE_setmem( ssize, twoParticlesIntegrals )
    implicit none 
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer(8) :: ssize, ssize2, ssize4

    ssize2 = (ssize * (ssize + 1))/2 
    ssize4 = (ssize2 * (ssize2 + 1))/2 

    allocate (twoParticlesIntegrals ( ssize4 ) )
!    allocate (twoParticlesIntegrals ( 1073741824 ) )
    !allocate (twoParticlesIntegrals ( 2147483648_8 ) )
    !ssize4 = ceiling(log(real(ssize4))/log(2.0))
    !ssize4 = 2**(ssize4)
    !print *, ssize4
    !allocate (twoParticlesIntegrals ( ssize4 ) )


    twoParticlesIntegrals = 0

  end subroutine TransformIntegralsE_setmem

  subroutine TransformIntegralsE_setmemInter( ssizea, ssizeb, twoParticlesIntegrals )
    implicit none 
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer(8) :: ssizea, ssizeb, ssize2a, ssize2b, ssize4

    ssize2a = (ssizea * (ssizea + 1))/2_8
    ssize2b = (ssizeb * (ssizeb + 1))/2_8
    ssize4 = (ssize2a * ssize2b) 

    allocate (twoParticlesIntegrals ( ssize4 ) )

    twoParticlesIntegrals = 0

  end subroutine TransformIntegralsE_setmemInter



  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine TransformIntegralsE_exception( typeMessage, description, debugDescription)
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

  end subroutine TransformIntegralsE_exception

end module TransformIntegralsE_
