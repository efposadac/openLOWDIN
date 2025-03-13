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
module TransformIntegralsC_
  use MolecularSystem_
  use InputCI_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
  use DirectIntegralManager_
  use omp_lib
  implicit none

  type, public :: TransformIntegralsC
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

     integer :: p_lowerOrbital, p_upperOrbital
     integer :: q_lowerOrbital, q_upperOrbital
     integer :: r_lowerOrbital, r_upperOrbital
     integer :: s_lowerOrbital, s_upperOrbital

     integer :: p_l, p_u
     integer :: q_l, q_u
     integer :: r_l, r_u
     integer :: s_l, s_u

     character(50) :: partialTransform



  end type TransformIntegralsC

  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE = 0
  integer, parameter :: TWO_SPECIES = 1
  !! }

  public :: &
       TransformIntegralsC_constructor, &
       TransformIntegralsC_destructor, &
       TransformIntegralsC_show, &
       TransformIntegralsC_atomicToMolecularOfOneSpecie, &
       TransformIntegralsC_atomicToMolecularOfTwoSpecies
  !       TransformIntegralsC_readIntegralsTransformed

  private

contains


  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsC_constructor(this,partial)
    implicit none
    type(TransformIntegralsC) :: this
    character(*) :: partial
    
    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"

    this%partialTransform=trim(partial)

    if ( trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then
       if (allocated(Libint2Instance)) call DirectIntegralManager_destructor(Libint2Instance)
       if (.not. allocated(Libint2Instance)) allocate(Libint2Instance(size(molecularSystem_instance%species)))
       call DirectIntegralManager_constructor(Libint2Instance, molecularSystem_instance)
    end if
    
  end subroutine TransformIntegralsC_constructor

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsC_destructor(this)
    implicit none
    type(TransformIntegralsC) :: this

  end subroutine TransformIntegralsC_destructor

  !>
  !! @brief show
  !<
  subroutine TransformIntegralsC_show()
    implicit none

    print *,"--------------------------------------------------"
    print *,"   4N^5 Algorithm Four-index integral tranformation"
    print *,"--------------------------------------------------"

  end subroutine TransformIntegralsC_show

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsC_atomicToMolecularOfOneSpecie( this, densityMatrix, coefficientsOfAtomicOrbitals, &
       speciesID, nameOfSpecie  )
    implicit none
    type(TransformIntegralsC) :: this
    type(Matrix) :: densityMatrix
    type(Matrix) :: coefficientsOfAtomicOrbitals
    integer :: speciesID
    character(*) :: nameOfSpecie
    integer :: integralStackSize

    integer :: status
    integer(8) :: ssize,ssize2

    !!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
!    real(8)  :: twoParticlesIntegrals(1073741824)
    real(8) , allocatable :: twoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable, target :: auxtempA(:,:,:)
    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)
    integer(8), allocatable :: xy(:,:)
    integer(8), allocatable :: ioff(:)
    integer(8) :: nonZeroIntegrals
    !$ real(8) :: timeA(10), timeB(10)

    integer :: pp(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: qq(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: auxIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer(8) :: ab, cd, ee, abcd, pq, rs, ij, kl, ij2, kl2

    integer :: p, q, r, s, mu, nu, auxnu, lambda, sigma, m, n, u, x, i, j, k, l
    integer(8) nstacks, stacktop
    integer(8) :: index, index2, pqrs, posxy, mm

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer(8) :: filesize
    logical :: symmetric

    if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    this%specieID = speciesID
    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)

    ssize = this%numberOfContractions
    call TransformIntegralsC_checkMOIntegralType(speciesID, this, symmetric)

    !call TransformIntegralsE_setmem( ssize, twoParticlesIntegrals )
    !call TransformIntegralsC_setSizeOfIndexArray( ssize, twoParticlesIntegrals)
    call TransformIntegralsC_setSizeOfIntegralsArray( ssize, twoParticlesIntegrals)


    ssize2 = (ssize * (ssize + 1))/2 

    allocate (xy ( ssize , ssize ) )
    allocate (ioff ( ssize2 ) )

    xy = 0
    ioff = 0

    m = 0
    do p = 1,  ssize
      do q = p,  ssize
        m = m + 1
        xy(p,q) = m        
        xy(q,p) = m        
      end do
    end do

    ioff(1) = 0
    do pq = 2, ssize2 
      ioff(pq) = ioff(pq-1) + ssize2 - pq + 1 
    end do

    !$  timeA(1) = omp_get_wtime()
    !! Read integrals

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

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

      !    print *, pp(i), qq(i), rr(i), ss(i), shellIntegrals(i)
      twoParticlesIntegrals(index2) = shellIntegrals(i)

    end do 

    close (unitid)

    !$OMP END PARALLEL

!$  timeB(1) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Loading AO integrals time(s): ", timeB(1) -timeA(1) 

    !! First half-transformation
!$  timeA(2) = omp_get_wtime()

    auxIntegrals = 0.0_8
    pp = 0
    qq = 0
    rr = 0
    ss = 0

    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    mm = 0

    call omp_set_num_threads(omp_get_max_threads())

    !! begin transformation
    mm = 0
    !$OMP PARALLEL &
    !$omp& private(p,q,r,s,j,n,u,k,l,ij,kl,tempA, tempB, tempC, index2, mu, nu,lambda, sigma, auxTransformedTwoParticlesIntegral) &
    !$omp& shared(pp,qq,rr,ss,m,auxIntegrals) reduction(+:mm) 
    if ( allocated (tempC)) deallocate (tempC )
    allocate (tempC ( ssize ) )
    tempC = 0

    if ( allocated (tempB)) deallocate (tempB )
    allocate (tempB ( ssize, ssize ) )
    tempB = 0

    if ( allocated (tempA)) deallocate (tempA )
    allocate (tempA ( ssize, ssize, ssize ) )
    tempA = 0

    !$omp do schedule (dynamic)
    do p = this%p_l, this%p_u
       n = p
       tempA = 0
       auxTransformedTwoParticlesIntegral = 0

       !! First quarter
       do mu = 1, this%numberOfContractions

          do j = 1, ssize
             ij = xy(j,mu)
             kl = 0
             do k = 1, ssize
                do l = k, ssize
                   kl = kl + 1

                   if ( ij >= kl ) then 
                     index2 = ioff(kl) + ij
                   else 
                     index2 = ioff(ij) + kl
                   end if

                   tempA(l,k,j) = tempA(l,k,j) + twoParticlesIntegrals(index2)*coefficientsOfAtomicOrbitals%values( mu, p ) 
                   tempA(k,l,j) = tempA(l,k,j) 
                end do
             end do
          end do

       end do

       do q = this%q_l, this%q_u
       !do q = p, this%q_u
          u = q
          tempB = 0

          !if ( q < this%q_l ) cycle
          if ( q < p .and. symmetric ) cycle
          !! second quarter
          do nu = 1, this%numberOfContractions
             !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
             tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                  tempA(:,:,nu)
          end do

          !do r = n, this%r_u
          do r = this%r_l, this%r_u

             tempC = 0

             !if ( r <  this%r_l  ) cycle
             if ( r < n  .and. symmetric ) cycle

             !! third quarter
             do lambda = 1, this%numberOfContractions
                !if ( abs(coefficientsOfAtomicOrbitals%values( lambda, r )) < 1E-10 ) cycle
                tempC(:) = tempC(:) + coefficientsOfAtomicOrbitals%values( lambda, r )* &
                     !tempB(lambda,:)
                     tempB(:,lambda)
             end do

             do s = this%s_l, this%s_u
             !do s = u, this%s_u
                auxTransformedTwoParticlesIntegral = 0

                !if ( s < this%s_l ) cycle
                if ( s < r .and. symmetric ) cycle
                !! fourth quarter
                do sigma = 1, this%numberOfContractions
                   auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                        coefficientsOfAtomicOrbitals%values( sigma, s )* &
                        tempC(sigma)

                end do

                if ( abs(auxTransformedTwoParticlesIntegral ) > 1E-10 ) then
                  !$omp critical
                  m = m + 1
                  auxIntegrals(m) = auxTransformedTwoParticlesIntegral
                  pp(m) = p
                  qq(m) = q
                  rr(m) = r
                  ss(m) = s

                  ! print *, p, q, r, s, auxTransformedTwoParticlesIntegral
                  
                  if (m == integralStackSize ) then
                    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
                    mm = mm + m
                    m = 0
                    auxIntegrals = 0
                    pp = 0
                    qq = 0
                    rr = 0
                    ss = 0
                  end if
                  !$omp end critical
                end if

             end do
             u = r + 1
          end do
       end do
    end do

    !$omp end do 

    !$omp critical
    mm = mm + m 
    m = m + 1
    pp(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
    !$omp end critical
 
    !$omp end parallel



!    close (unittmp)
    deallocate (twoParticlesIntegrals)

!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Integral transformation time(s): ", timeB(2) -timeA(2) 
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

    write (*,"(T4,A36,I12)") "Non-zero transformed integrals: ", mm

    !!!!!!!!!!!!!!!!!!!!!!!!
    else !! DIRECT

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    this%specieID = speciesID

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    ssize = this%numberOfContractions
    
    call TransformIntegralsC_checkMOIntegralType(speciesID, this, symmetric)

    !call TransformIntegralsE_setmem( ssize, twoParticlesIntegrals )
    !call TransformIntegralsC_setSizeOfIndexArray( ssize, twoParticlesIntegrals)
    !!call TransformIntegralsC_setSizeOfIntegralsArray( ssize, twoParticlesIntegrals)

    ssize2 = (ssize * (ssize + 1))/2 

    allocate (xy ( ssize , ssize ) )
    allocate (ioff ( ssize2 ) )

    xy = 0
    ioff = 0

    m = 0
    do p = 1,  ssize
      do q = p,  ssize
        m = m + 1
        xy(p,q) = m        
        xy(q,p) = m        
      end do
    end do

    ioff(1) = 0
    do pq = 2, ssize2 
      ioff(pq) = ioff(pq-1) + ssize2 - pq + 1 
    end do

    !! First half-transformation
!$  timeA(2) = omp_get_wtime()

    allocate (tempB ( ssize, ssize ) )
    tempB = 0

    if ( allocated (tempC)) deallocate (tempC )
    allocate (tempC ( ssize ) )

    tempC = 0

    !if ( allocated (tempA)) deallocate (tempA )
    !allocate (tempA ( ssize, ssize, ssize ) )
    !tempA = 0

    auxIntegrals = 0.0_8
    pp = 0
    qq = 0
    rr = 0
    ss = 0

    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    mm = 0

    !call omp_set_num_threads(omp_get_max_threads())
    !print *, "fortran", coefficientsOfAtomicOrbitals%values(1,2)
    !print *, "fortran", coefficientsOfAtomicOrbitals%values(2,1)

    !! begin transformation
    mm = 0
    do p = this%p_l, this%p_u
       n = p
       !tempA = 0
       auxTransformedTwoParticlesIntegral = 0

       !! First quarter
       !do mu = 1, this%numberOfContractions

          call DirectIntegralManager_getDirectIntraRepulsionFirstQuarter(&
               speciesID, &
               trim(CONTROL_instance%INTEGRAL_SCHEME), &
               densityMatrix, & 
               coefficientsOfAtomicOrbitals, &
               auxtempA, p-1 )

          !do j = 1, ssize
          !   !ij = int(IndexMap_tensorR2ToVectorB( int(mu,8), int(j,8), int(ssize,8)), 4)
          !   ij = xy(j,mu)
          !   kl = 0
          !   do k = 1, ssize
          !      do l = k, ssize
          !         kl = kl + 1
          !                  !auxIndex = IndexMap_tensorR4ToVectorB(int(i,8),int(j,8),int(k,8),int(l,8),int(ssize,8))
          !         !index2 =(IndexMap_tensorR2ToVectorB( int(ij,8), int(kl,8), int(ssize2,8)))
          !         if ( ij >= kl ) then 
          !           index2 = ioff(kl) + ij
          !         else 
          !           index2 = ioff(ij) + kl
          !         end if
          !      
          !         tempA(l,k,j) = tempA(l,k,j) + twoParticlesIntegrals(index2)*coefficientsOfAtomicOrbitals%values( mu, p ) 
          !         tempA(k,l,j) = tempA(l,k,j) 
          !      end do
          !   end do
          !end do

          !tempA(:,:,:) = tempA(:,:,:) + coefficientsOfAtomicOrbitals%values( mu, p ) * & 
          !     auxtempA(:,:,:)

       !end do

       do q = p, this%q_u
          u = q
          tempB = 0

          if ( q < this%q_l ) cycle
          !! second quarter
          do nu = 1, this%numberOfContractions
             !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
             tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                  auxtempA(:,:,nu)
          end do

          do r = n, this%r_u

             tempC = 0

             !if ( r <  this%r_l  ) cycle

             !! third quarter
             do lambda = 1, this%numberOfContractions
                !if ( abs(coefficientsOfAtomicOrbitals%values( lambda, r )) < 1E-10 ) cycle
                tempC(:) = tempC(:) + coefficientsOfAtomicOrbitals%values( lambda, r )* &
                     !tempB(lambda,:)
                     tempB(:,lambda)
             end do

             do s = u, this%s_u
                auxTransformedTwoParticlesIntegral = 0

                if ( s < this%s_l ) cycle
                !! fourth quarter
                do sigma = 1, this%numberOfContractions
                   auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                        coefficientsOfAtomicOrbitals%values( sigma, s )* &
                        tempC(sigma)

                end do

                if ( abs(auxTransformedTwoParticlesIntegral ) > 1E-10 ) then
                  m = m + 1
                  auxIntegrals(m) = auxTransformedTwoParticlesIntegral
                  pp(m) = p
                  qq(m) = q
                  rr(m) = r
                  ss(m) = s
        
                  if (m == integralStackSize ) then
                    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
                    mm = mm + m
                    m = 0
                    auxIntegrals = 0
                    pp = 0
                    qq = 0
                    rr = 0
                    ss = 0
                  end if
                end if

             end do
             u = r + 1
          end do
       end do
    end do

    mm = mm + m 
    m = m + 1
    pp(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals

!    close (unittmp)
!    deallocate (twoParticlesIntegrals)

!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Integral transformation time(s): ", timeB(2) -timeA(2) 
    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

    write (*,"(T4,A36,I12)") "Non-zero transformed integrals: ", mm

    end if
    
  end subroutine TransformIntegralsC_atomicToMolecularOfOneSpecie

  subroutine TransformIntegralsC_setSizeOfIntegralsArray ( numberOfContractions, twoParticlesIntegrals)
    implicit none 
    integer(8) :: numberOfContractions
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer :: ssize 
    integer(8) :: ssize8 !! Beyond 360 cartesian funtions

       ssize8 = numberOfContractions
       ssize8 = (ssize8 * (ssize8 + 1))/2 
       ssize8 = (ssize8 * (ssize8 + 1))/2 

       if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
       allocate (twoParticlesIntegrals ( ssize8 ) )

       twoParticlesIntegrals = 0

  end subroutine TransformIntegralsC_setSizeOfIntegralsArray

  subroutine TransformIntegralsC_buildArrayA( integralArray, i, ssize , ssize2, auxtempA)
    implicit none 
    integer, intent(in) :: ssize,ssize2
    integer, intent(in) :: i
    real(8), intent(in) :: integralArray(:)
    real(8) :: auxtempA(ssize,ssize,ssize)
    integer :: j,k,l,ij, kl

    ij = 0

    !$OMP PARALLEL DO private(j,k,l,ij,kl) shared(ssize,ssize2,i,integralArray)
    do j = 1, ssize
       ij = int(IndexMap_tensorR2ToVectorB( int(i,8), int(j,8), int(ssize,8)), 4)
       kl = 0
       do k = 1, ssize
          do l = k, ssize
             kl = kl + 1
             !         auxIndex = IndexMap_tensorR4ToVectorB(int(i,8),int(j,8),int(k,8),int(l,8),int(ssize,8))
             auxtempA(j,k,l) = integralArray(IndexMap_tensorR2ToVectorB( int(ij,8), int(kl,8), int(ssize2,8)))
          end do
       end do
    end do
    !$OMP END PARALLEL DO  

    !$OMP PARALLEL DO private(k,l) shared(ssize)
    do l = 1, ssize
       do k = l+1, ssize
          auxtempA(:,k,l) = auxtempA(:,l,k)
       end do
    end do
    !$OMP END PARALLEL DO  

  end subroutine TransformIntegralsC_buildArrayA


  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !!    a integrales moleculares.
  !<
  !!  subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
  !!       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
  subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies( this, densityMatrix, coefficientsOfAtomicOrbitals, &
       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, &
       otherSpecieID, nameOfOtherSpecie )
    implicit none
    type(TransformIntegralsC) :: this
    type(Matrix) :: densityMatrix
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID, otherSpecieID
    character(*) :: nameOfSpecie, nameOfOtherSpecie
    character(50) :: fileName
    integer :: integralStackSize

    integer :: status
    integer :: nonZeroIntegrals

    !!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8), allocatable :: twoParticlesIntegrals(:)
    !!integer, allocatable :: indexTwoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: auxtempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer(8), allocatable :: xya(:,:), xyb(:,:)
    integer(8), allocatable :: ioffa(:), ioffb(:)
    integer(8) :: ssizea,ssizeb, ssize2a, ssize2b 
    !$ real(8) :: timeA(10), timeB(10)

    integer :: pp(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: qq(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: auxIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m, mm, i, j, k, l
    integer(8) :: ab, cd, ee, abcd, pq, rs, ij, kl, ij2, kl2
    integer(8) :: index2

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer(8) :: index
    integer(8) :: filesize
    logical :: symmetric

    ! Reads the number of cores

    if ( .not. trim(String_getUppercase(CONTROL_instance%INTEGRAL_STORAGE)) == "DIRECT" ) then

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    ! call TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals( specieID, otherSpecieID, nproc, nonZeroIntegrals )

    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%otherNumberOfContractions=size(otherCoefficientsOfAtomicOrbitals%values,dim=1)
    
    call TransformIntegralsC_checkInterMOIntegralType(specieID, otherSpecieID, this, symmetric)

    !if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    !allocate (twoParticlesIntegrals ( nonZeroIntegrals ) )
    !twoParticlesIntegrals = 0

    !! Setting size of index array
    call TransformIntegralsC_setSizeOfInterIntegralsArray ( this%numberOfContractions, this%otherNumberOfContractions, &
         twoParticlesIntegrals)

    this%specieID = specieID

    !! change the order to ge the AO integrals
    ssizea = this%numberOfContractions
    ssizeb = this%otherNumberOfContractions

    ssize2a = (ssizea * (ssizea + 1))/2 
    ssize2b = (ssizeb * (ssizeb + 1))/2 

    allocate (xya ( ssizea , ssizea ) )
    allocate (xyb ( ssizeb , ssizeb ) )

    xya = 0
    xyb = 0

    m = 0
    do p = 1,  ssizea
      do q = p,  ssizea
        m = m + 1
        xya(p,q) = m        
        xya(q,p) = m        
      end do
    end do

    m = 0
    do p = 1,  ssizeb
      do q = p,  ssizeb
        m = m + 1
        xyb(p,q) = m        
        xyb(q,p) = m        
      end do
    end do


    if ( specieID < otherSpecieID ) then

       if ( trim(nameOfSpecie) == "E-ALPHA" .and. trim(nameOfOtherSpecie) == "E-BETA" ) then
          fileName = "E-ALPHA.E-BETA"
       else if( trim(nameOfOtherSpecie) == "E-BETA" ) then
          fileName = trim(nameOfSpecie)//".E-ALPHA"
       else if(trim(nameOfSpecie) == "E-BETA") then
          fileName = "E-ALPHA."//trim(nameOfOtherSpecie)
       else 
          fileName = trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
       end if

    !! Read integrals
!$  timeA(1) = omp_get_wtime()
    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    !! open file for integrals

      open(UNIT=unitid,FILE=trim(fileid)//trim(fileName)//".ints", &
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

      index2 = (rs-1)*ssize2a + pq
      twoParticlesIntegrals(index2) = shellIntegrals(i)

    end do 

    close (unitid)

    !$OMP END PARALLEL

    else

       if ( trim(nameOfOtherSpecie) == "E-ALPHA" .and. trim(nameOfSpecie) == "E-BETA" ) then
          fileName = "E-ALPHA.E-BETA"
       else if( trim(nameOfOtherSpecie) == "E-BETA" ) then
          fileName = "E-ALPHA."//trim(nameOfSpecie)
       else if(trim(nameOfSpecie) == "E-BETA") then
          fileName = trim(nameOfOtherSpecie)//".E-ALPHA"
       else 
          fileName = trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)
       end if

    !! Read integrals
!$  timeA(1) = omp_get_wtime()
    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, pp, qq, rr, ss, p, shellIntegrals, i, index2, filesize, pq, rs)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    !! open file for integrals

      open(UNIT=unitid,FILE=trim(fileid)//trim(fileName)//".ints", &
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

          pq = xyb(pp(i),qq(i))
          rs = xya(rr(i),ss(i))

          index2 = (pq-1)*ssize2a + rs
          twoParticlesIntegrals(index2) = shellIntegrals(i)

      end do
    end do 

    !! read last stack
    read(UNIT=unitid, iostat=status) pp, qq, rr, ss, shellIntegrals

    do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
      if( pp(i) == -1 ) exit

      pq = xyb(pp(i),qq(i))
      rs = xya(rr(i),ss(i))

      index2 = (pq-1)*ssize2a + rs
      twoParticlesIntegrals(index2) = shellIntegrals(i)

    end do 

    close (unitid)

    !$OMP END PARALLEL


    end if




!$  timeB(1) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Loading AO integrals time(s): ", timeB(1) -timeA(1) 

!$  timeA(2) = omp_get_wtime()

    deallocate (xya )
    deallocate (xyb )

    ssizea = this%numberOfContractions
    ssizeb = this%otherNumberOfContractions

    ssize2a = (ssizea * (ssizea + 1))/2 
    ssize2b = (ssizeb * (ssizeb + 1))/2 

    allocate (xya ( ssizea , ssizea ) )
    allocate (xyb ( ssizeb , ssizeb ) )

    xya = 0
    xyb = 0

    m = 0
    do p = 1,  ssizea
      do q = p,  ssizea
        m = m + 1
        xya(p,q) = m        
        xya(q,p) = m        
      end do
    end do

    m = 0
    do p = 1,  ssizeb
      do q = p,  ssizeb
        m = m + 1
        xyb(p,q) = m        
        xyb(q,p) = m        
      end do
    end do



    !! allocate some auxiliary arrays
    if ( allocated (tempA)) deallocate (tempA )
    allocate (tempA ( ssizeb , &
         ssizeb, &
         ssizea ) )
    tempA = 0

    if ( allocated (auxtempA)) deallocate (auxtempA )
    allocate (auxtempA ( this%numberOfContractions , &
         this%otherNumberOfContractions, &
         this%otherNumberOfContractions ) )
    auxtempA = 0


    if ( allocated (tempB)) deallocate (tempB )
    allocate (tempB ( ssizeb, ssizeb ) )
    tempB = 0

    if ( allocated (tempC)) deallocate (tempC )
    allocate (tempC ( ssizeb ) )
    tempC = 0

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    call omp_set_num_threads(omp_get_max_threads())

    !! begin transformation
    m = 0
    mm = 0
    do p = this%p_lowerOrbital, this%p_upperOrbital

       tempA = 0
       auxTransformedTwoParticlesIntegral = 0

       !! First quarter
       do mu = 1, this%numberOfContractions
          !if ( abs(coefficientsOfAtomicOrbitals%values( mu, p )) < 1E-10 ) cycle
          !! auxtemp is the twoparticlesintegrals reduced to a three dimensional array
          !call TransformIntegralsC_buildArrayAInter( twoParticlesIntegrals, mu, &
          !     this%numberOfContractions, this%otherNumberOfContractions, int(ssize2b,4), &
          !     auxtempA )

          do j = 1, ssizea
            !ij = int(IndexMap_tensorR2ToVectorB( int(i,8), int(j,8), int(ssize,8)), 4)
            ij = xya(j,mu)
            kl = 0
            do k = 1, ssizeb
              do l = k, ssizeb
                kl = kl + 1
                !         auxm = indexArray(auxIndex)
                index2 = (kl-1)*ssize2a + ij
                !auxtempA(j,k,l) = integralArray( ( int(ij,8) - 1 ) * int(auxOtherSsize,8) + int(kl,8) )
                tempA(l,k,j) = tempA(l,k,j) + twoParticlesIntegrals(index2) * coefficientsOfAtomicOrbitals%values( mu, p ) 
                tempA(k,l,j) = tempA(l,k,j)
              end do
            end do
          end do

          !tempA(:,:,:) = tempA(:,:,:) + coefficientsOfAtomicOrbitals%values( mu, p ) * & 
          !     auxtempA(:,:,:)

       end do

       !do q = p, this%q_upperOrbital
       do q = this%q_lowerOrbital, this%q_upperOrbital
          tempB = 0

          !if ( q < this%q_lowerOrbital ) cycle
          if ( q < p .and. symmetric ) cycle
          !! second quarter
          do nu = 1, this%numberOfContractions
             !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle

             tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                  tempA(:,:,nu)
          end do

          do r = this%r_lowerOrbital , this%r_upperOrbital

             tempC = 0

             !!if ( r >  this%upperOccupiedOrbital  ) cycle

             !! third quarter
             do lambda = 1, this%otherNumberOfContractions

                tempC(:) = tempC(:) + otherCoefficientsOfAtomicOrbitals%values( lambda, r )* &
                     tempB(:,lambda)

             end do
             do s = this%s_lowerOrbital, this%s_upperOrbital
                auxTransformedTwoParticlesIntegral = 0

                !if ( s < this%s_lowerOrbital ) cycle
                if ( s < r .and. symmetric ) cycle
                !! fourth quarter
                do sigma = 1, this%otherNumberOfContractions
                   auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                        otherCoefficientsOfAtomicOrbitals%values( sigma, s )* &
                        tempC(sigma)

                end do
                !write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral
                !mm = mm + 1
                if ( abs(auxTransformedTwoParticlesIntegral ) > 1E-10 ) then
                  !!$omp critical
                  m = m + 1
                  auxIntegrals(m) = auxTransformedTwoParticlesIntegral
                  pp(m) = p
                  qq(m) = q
                  rr(m) = r
                  ss(m) = s

                  if (m == integralStackSize ) then
                    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
          
                    mm = mm + m
                    m = 0
                    auxIntegrals = 0
                    pp = 0
                    qq = 0
                    rr = 0
                    ss = 0
                  end if
                  !!$omp end critical
                end if



             end do
          end do
       end do
    end do

    !!$omp critical
    mm = mm + m 
    m = m + 1
    pp(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
    !!$omp end critical
 
!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Integral transformation time(s): ", timeB(2) -timeA(2) 
    print *, "Non zero transformed coupling integrals: ", mm

    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

    else !! DIRECT

    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    ! call TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals( specieID, otherSpecieID, nproc, nonZeroIntegrals )

    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%otherNumberOfContractions=size(otherCoefficientsOfAtomicOrbitals%values,dim=1)

    call TransformIntegralsC_checkInterMOIntegralType(specieID, otherSpecieID, this, symmetric)

    !if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    !allocate (twoParticlesIntegrals ( nonZeroIntegrals ) )
    !twoParticlesIntegrals = 0


    !! Setting size of index array
    !!call TransformIntegralsC_setSizeOfInterIntegralsArray ( this%numberOfContractions, this%otherNumberOfContractions, &
    !!     twoParticlesIntegrals)

    this%specieID = specieID

    !! change the order to ge the AO integrals
    ssizea = this%numberOfContractions
    ssizeb = this%otherNumberOfContractions

    ssize2a = (ssizea * (ssizea + 1))/2 
    ssize2b = (ssizeb * (ssizeb + 1))/2 

    allocate (xya ( ssizea , ssizea ) )
    allocate (xyb ( ssizeb , ssizeb ) )

    xya = 0
    xyb = 0

    m = 0
    do p = 1,  ssizea
      do q = p,  ssizea
        m = m + 1
        xya(p,q) = m        
        xya(q,p) = m        
      end do
    end do

    m = 0
    do p = 1,  ssizeb
      do q = p,  ssizeb
        m = m + 1
        xyb(p,q) = m        
        xyb(q,p) = m        
      end do
    end do

!$  timeA(2) = omp_get_wtime()

    !! allocate some auxiliary arrays
    !if ( allocated (auxtempA)) deallocate (auxtempA )
    !allocate (auxtempA ( this%numberOfContractions , &
    !     this%otherNumberOfContractions, &
    !     this%otherNumberOfContractions ) )
    !auxtempA = 0


    if ( allocated (tempB)) deallocate (tempB )
    allocate (tempB ( ssizeb, ssizeb ) )
    tempB = 0

    if ( allocated (tempC)) deallocate (tempC )
    allocate (tempC ( ssizeb ) )
    tempC = 0

    auxIntegrals = 0.0_8
    pp = 0
    qq = 0
    rr = 0
    ss = 0

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    !call omp_set_num_threads(omp_get_max_threads())

    !! begin transformation
    m = 0
    mm = 0
    do p = this%p_lowerOrbital, this%p_upperOrbital

       auxTransformedTwoParticlesIntegral = 0

       !! First quarter
       !do mu = 1, this%numberOfContractions

          call DirectIntegralManager_getDirectInterRepulsionFirstQuarter(&
               specieID, otherSpecieID, &
               trim(CONTROL_instance%INTEGRAL_SCHEME), &
               densityMatrix, & 
               coefficientsOfAtomicOrbitals, &
               auxtempA, p-1 )

       !   do j = 1, ssizea
       !     ij = xya(j,mu)
       !     kl = 0
       !     do k = 1, ssizeb
       !       do l = k, ssizeb
       !         kl = kl + 1
       !         index2 = (kl-1)*ssize2a + ij
       !         tempA(l,k,j) = tempA(l,k,j) + twoParticlesIntegrals(index2) * coefficientsOfAtomicOrbitals%values( mu, p ) 
       !         tempA(k,l,j) = tempA(l,k,j)
       !       end do
       !     end do
       !   end do

       !end do

       !do q = p, this%q_upperOrbital
       do q = this%q_lowerOrbital, this%q_upperOrbital
          tempB = 0

          !if ( q < this%q_lowerOrbital ) cycle
          if ( q < p .and. symmetric ) cycle
          !! second quarter
          do nu = 1, this%numberOfContractions
             !if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle

             tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                  auxtempA(:,:,nu)
          end do

          do r = this%r_lowerOrbital , this%r_upperOrbital

             tempC = 0

             !!if ( r >  this%upperOccupiedOrbital  ) cycle

             !! third quarter
             do lambda = 1, this%otherNumberOfContractions

                tempC(:) = tempC(:) + otherCoefficientsOfAtomicOrbitals%values( lambda, r )* &
                     tempB(:,lambda)

             end do
             do s = this%s_lowerOrbital, this%s_upperOrbital
                auxTransformedTwoParticlesIntegral = 0

                !if ( s < this%s_lowerOrbital ) cycle
                if ( s < r .and. symmetric ) cycle
                !! fourth quarter
                do sigma = 1, this%otherNumberOfContractions
                   auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                        otherCoefficientsOfAtomicOrbitals%values( sigma, s )* &
                        tempC(sigma)

                end do
                !write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral
                !mm = mm + 1
                if ( abs(auxTransformedTwoParticlesIntegral ) > 1E-10 ) then
                  !!$omp critical
                  m = m + 1
                  auxIntegrals(m) = auxTransformedTwoParticlesIntegral
                  pp(m) = p
                  qq(m) = q
                  rr(m) = r
                  ss(m) = s
                  if (m == integralStackSize ) then
                    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
          
                    mm = mm + m
                    m = 0
                    auxIntegrals = 0
                    pp = 0
                    qq = 0
                    rr = 0
                    ss = 0
                  end if
                  !!$omp end critical
                end if



             end do
          end do
       end do
    end do

    !!$omp critical
    mm = mm + m 
    m = m + 1
    pp(m) = -1_8

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) pp, qq, rr, ss, auxIntegrals
    !!$omp end critical
 
!$  timeB(2) = omp_get_wtime()
!$  write(*,"(T4,A36,E10.3)") "Integral transformation time(s): ", timeB(2) -timeA(2) 
    print *, "Non zero transformed coupling integrals: ", mm

    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)

    end if

  end subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies



  subroutine TransformIntegralsC_setSizeOfInterIntegralsArray ( numberOfContractions, otherNumberOfContractions, &
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

  end subroutine TransformIntegralsC_setSizeOfInterIntegralsArray


  subroutine TransformIntegralsC_buildArrayAInter( integralArray, i, ssize , otherSsize, auxOtherSsize, auxtempA)
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

  end subroutine TransformIntegralsC_buildArrayAInter

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsC_checkMOIntegralType(speciesID, this, symmetric)
    implicit none
    integer :: speciesID
    type(TransformIntegralsC) :: this
    integer :: totalOccupation 
    integer :: coreOrbitals
    integer :: totalActiveOrbitals
    logical :: symmetric

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )

    !! Take the active space from input
    coreOrbitals=0
    if ( InputCI_Instance(speciesID)%coreOrbitals /= 0 ) coreOrbitals=InputCI_Instance(speciesID)%coreOrbitals 

    totalActiveOrbitals=this%numberOfContractions
    if ( InputCI_Instance(speciesID)%activeOrbitals /= 0 ) totalActiveOrbitals=InputCI_Instance(speciesID)%activeOrbitals 

    symmetric = .true. 

    !! Boundary orbitals. Default
    this%p_l = coreOrbitals+1
    this%p_u = totalActiveOrbitals
    this%q_l = coreOrbitals+1
    this%q_u = totalActiveOrbitals
    this%r_l = coreOrbitals+1
    this%r_u = totalActiveOrbitals
    this%s_l = coreOrbitals+1
    this%s_u = totalActiveOrbitals

    if ( trim(this%partialTransform)=="ALL") then
       this%p_l = 1
       this%p_u = this%numberOfContractions
       this%q_l = 1
       this%q_u = this%numberOfContractions
       this%r_l = 1
       this%r_u = this%numberOfContractions
       this%s_l = 1
       this%s_u = this%numberOfContractions !this%numberOfContractions
    end if

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

       this%p_l = coreOrbitals+1
       this%p_u = totalOccupation 
       this%q_l = totalOccupation+1 !+coreOrbitals+1
       this%q_u = totalActiveOrbitals

       this%r_l = coreOrbitals+1
       this%r_u = totalOccupation 
       this%s_l = totalOccupation+1 !+coreOrbitals+1
       this%s_u = totalActiveOrbitals

    end if

    !! only the (ip|aq) integrals will be transformed
    if ( trim(this%partialTransform)=="PT2"  ) then

       symmetric = .false. 
       if ( CONTROL_instance%IONIZE_MO(1) == 0 ) then
          !! all
          this%p_l = totalOccupation     !! HOMO 
          this%p_u = totalOccupation + 1 !! LUMO
          this%q_l = coreOrbitals+1
          this%q_u = totalActiveOrbitals

          this%r_l = coreOrbitals+1
          this%r_u = totalOccupation 
          this%s_l = totalOccupation + 1
          this%s_u = totalActiveOrbitals

          !! half symmetric
          !this%p_l = 1 
          !this%p_u = totalOccupation + 1
          !this%q_l = 1
          !this%q_u = totalActiveOrbitals

          !this%r_l = 1
          !this%r_u = totalOccupation 
          !this%s_l = totalOccupation + 1
          !this%s_u = totalActiveOrbitals

          !! fully symmetric
          !this%p_l = 1 
          !this%p_u = totaloccupation  
          !this%q_l = 1
          !this%q_u = totalActiveOrbitals

          !this%r_l = 1
          !this%r_u = totalActiveOrbitals
          !this%s_l =  1
          !this%s_u = totalActiveOrbitals

       else

          if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
             this%p_l = CONTROL_instance%IONIZE_MO(1)  
             this%p_u = CONTROL_instance%IONIZE_MO(1) 
             this%q_l = coreOrbitals+1
             this%q_u = totalActiveOrbitals

             this%r_l = coreOrbitals+1
             this%r_u = totalOccupation !totalActiveOrbitals
             this%s_l = coreOrbitals+1
             this%s_u = totalActiveOrbitals
          else 

             this%p_l = CONTROL_instance%IONIZE_MO(1)  
             this%p_u = CONTROL_instance%IONIZE_MO(1) 
             this%q_l = coreOrbitals+1
             this%q_u = totalActiveOrbitals

             this%r_l = coreOrbitals+1
             this%r_u = totalOccupation 
             this%s_l = totalOccupation + 1
             this%s_u = totalActiveOrbitals

          end if
       end if

    end if

    !for a simultaneous PT2 - MP2 calculation
    if ( trim(this%partialTransform)=="MP2-PT2" ) then

       symmetric = .false. 
       if ( CONTROL_instance%IONIZE_MO(1) == 0 ) then

          this%p_l = coreOrbitals+1
          this%p_u = totalOccupation + 1
          this%q_l = coreOrbitals+1
          this%q_u = totalActiveOrbitals

          this%r_l = coreOrbitals+1
          this%r_u = totalOccupation 
          this%s_l = totalOccupation + 1
          this%s_u = totalActiveOrbitals
       else

          if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
             this%p_l = coreOrbitals+1
             this%p_u = max(CONTROL_instance%IONIZE_MO(1),totalOccupation)
             this%q_l = coreOrbitals+1
             this%q_u = totalActiveOrbitals

             this%r_l = coreOrbitals+1
             this%r_u = totalOccupation !totalActiveOrbitals
             this%s_l = coreOrbitals+1
             this%s_u = totalActiveOrbitals
          else 

             this%p_l = coreOrbitals+1  
             this%p_u = max(CONTROL_instance%IONIZE_MO(1),totalOccupation)
             this%q_l = coreOrbitals+1
             this%q_u = totalActiveOrbitals

             this%r_l = coreOrbitals+1
             this%r_u = totalOccupation 
             this%s_l = totalOccupation + 1
             this%s_u = totalActiveOrbitals

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
   
  end subroutine TransformIntegralsC_checkMOIntegralType


  subroutine TransformIntegralsC_checkInterMOIntegralType(speciesID, otherSpeciesID, this, symmetric)
    implicit none
    integer :: speciesID, otherSpeciesID
    type(TransformIntegralsC) :: this
    integer :: totalOccupation, otherTotalOccupation
    integer :: coreOrbitals, otherCoreOrbitals
    integer :: totalActiveOrbitals, otherTotalActiveOrbitals
    logical :: ionizeA, ionizeB
    integer :: s
    character(10) :: nameOfSpecies
    character(10) :: nameOfOtherSpecies
    logical :: symmetric

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )
    otherTotalOccupation = MolecularSystem_getOcupationNumber( otherSpeciesID )

    !! Take the active space from input
    coreOrbitals=0
    otherCoreOrbitals=0
    if ( InputCI_Instance(speciesID)%coreOrbitals /= 0 ) coreOrbitals=InputCI_Instance(speciesID)%coreOrbitals
    if ( InputCI_Instance(otherSpeciesID)%coreOrbitals /= 0 ) otherCoreOrbitals=InputCI_Instance(otherSpeciesID)%coreOrbitals

    totalActiveOrbitals=this%numberOfContractions
    otherTotalActiveOrbitals=this%otherNumberOfContractions
    if ( InputCI_Instance(speciesID)%activeOrbitals /= 0 ) totalActiveOrbitals=InputCI_Instance(speciesID)%activeOrbitals 
    if ( InputCI_Instance(otherSpeciesID)%activeOrbitals /= 0 ) otherTotalActiveOrbitals=InputCI_Instance(otherSpeciesID)%activeOrbitals 
    
    
    symmetric = .true.
    !! Boundary orbitals. Default
    this%p_lowerOrbital = coreOrbitals+1
    this%p_upperOrbital = totalActiveOrbitals
    this%q_lowerOrbital = coreOrbitals+1
    this%q_upperOrbital = totalActiveOrbitals
    this%r_lowerOrbital = otherCoreOrbitals+1
    this%r_upperOrbital = otherTotalActiveOrbitals
    this%s_lowerOrbital = otherCoreOrbitals+1
    this%s_upperOrbital = otherTotalActiveOrbitals

    if ( trim(this%partialTransform) .eq. "ALL") then
       this%p_lowerOrbital = 1
       this%p_upperOrbital = this%numberOfContractions
       this%q_lowerOrbital = 1
       this%q_upperOrbital = this%numberOfContractions
       this%r_lowerOrbital = 1
       this%r_upperOrbital = this%otherNumberOfContractions
       this%s_lowerOrbital = 1
       this%s_upperOrbital = this%otherNumberOfContractions
    end if

    if ( trim(this%partialTransform) .eq. "ALLACTIVE") then
       this%p_lowerOrbital = 1
       this%p_upperOrbital = totalActiveOrbitals!this%numberOfContractions
       this%q_lowerOrbital = 1
       this%q_upperOrbital = totalActiveOrbitals!this%numberOfContractions
       this%r_lowerOrbital = 1
       this%r_upperOrbital = otherTotalActiveOrbitals!this%otherNumberOfContractions
       this%s_lowerOrbital = 1
       this%s_upperOrbital = otherTotalActiveOrbitals!this%otherNumberOfContractions
    end if
    
    !! only the (ia|jb) integrals will be transformed
    if ( trim(this%partialTransform) .eq. "MP2"  ) then

       this%p_lowerOrbital = coreOrbitals+1
       this%p_upperOrbital = totalOccupation
       this%q_lowerOrbital = totalOccupation + 1
       this%q_upperOrbital = totalActiveOrbitals
       this%r_lowerOrbital = otherCoreOrbitals+1
       this%r_upperOrbital = otherTotalOccupation
       this%s_lowerOrbital = otherTotalOccupation + 1
       this%s_upperOrbital = otherTotalActiveOrbitals

    end if

    !! only the (ip|IP) integrals will be transformed.
    if ( trim(this%partialTransform) .eq. "PT2" ) then

       this%p_lowerOrbital = coreOrbitals+1
       this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
       this%q_lowerOrbital = coreOrbitals+1
       this%q_upperOrbital = totalActiveOrbitals
       this%r_lowerOrbital = otherCoreOrbitals+1
       this%r_upperOrbital = otherTotalOccupation + 1 !!occ + lumo
       this%s_lowerOrbital = otherCoreOrbitals+1
       this%s_upperOrbital = otherTotalActiveOrbitals
       symmetric = .true.

      if (CONTROL_instance%IONIZE_SPECIES(1) /= "NONE" ) then

        symmetric = .false.
        ionizeA = .false.
        ionizeB = .false.

         nameOfSpecies= trim(  MolecularSystem_getNameOfSpecies( speciesID ) )
         nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecies( otherSpeciesID ) )

         do s = 1, size(CONTROL_instance%IONIZE_SPECIES )
           if ( nameOfSpecies == trim(CONTROL_instance%IONIZE_SPECIES(s)) ) then
             ionizeA = .true. 
           end if
           if ( nameOfOtherSpecies == trim(CONTROL_instance%IONIZE_SPECIES(s)) ) then
             ionizeB = .true. 
           end if
         end do

        if ( CONTROL_instance%IONIZE_MO(1) == 0 ) then

          if ( ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = otherTotalOccupation + 1 !!occ + lumo
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
  
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = othertotaloccupation 
            this%s_lowerOrbital = othertotaloccupation + 1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation
            this%q_lowerOrbital = totalOccupation + 1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = otherTotalOccupation + 1
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          end if

        else if ( CONTROL_instance%IONIZE_MO(1) /= 0 ) then !!occ and vir..

          if ( ionizeA .and. ionizeB ) then
            if ( CONTROL_instance%IONIZE_MO(1) <= totalOccupation .and. CONTROL_instance%IONIZE_MO(1) <= othertotalOccupation ) then
              this%p_lowerOrbital = coreOrbitals+1
              this%p_upperOrbital = totalOccupation
              this%q_lowerOrbital = coreOrbitals+1
              this%q_upperOrbital = totalActiveOrbitals
              this%r_lowerOrbital = otherCoreOrbitals+1
              this%r_upperOrbital = otherTotalOccupation
              this%s_lowerOrbital = otherCoreOrbitals+1
              this%s_upperOrbital = otherTotalActiveOrbitals

            else if ( CONTROL_instance%IONIZE_MO(1) > totalOccupation .and. CONTROL_instance%IONIZE_MO(1) > othertotalOccupation ) then

              this%p_lowerOrbital = coreOrbitals+1
              this%p_upperOrbital = CONTROL_instance%IONIZE_MO(1)
              this%q_lowerOrbital = coreOrbitals+1
              this%q_upperOrbital = totalActiveOrbitals
              this%r_lowerOrbital = otherCoreOrbitals+1
              this%r_upperOrbital = CONTROL_instance%IONIZE_MO(1) 
              this%s_lowerOrbital = otherCoreOrbitals+1
              this%s_upperOrbital = otherTotalActiveOrbitals

            end if

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%p_lowerOrbital = CONTROL_instance%IONIZE_MO(1)
            this%p_upperOrbital = CONTROL_instance%IONIZE_MO(1)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
  
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = othertotaloccupation 
            this%s_lowerOrbital = othertotaloccupation + 1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation
            this%q_lowerOrbital = totalOccupation + 1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = CONTROL_instance%IONIZE_MO(1)
            this%r_upperOrbital = CONTROL_instance%IONIZE_MO(1)
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          end if

        end if

          !if ( CONTROL_instance%IONIZE_MO /= 0 ) then
      end if

   end if

   !for a simultaneous PT2 - MP2 calculation
   if ( trim(this%partialTransform) .eq. "MP2-PT2"  ) then

       this%p_lowerOrbital = coreOrbitals+1
       this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
       this%q_lowerOrbital = coreOrbitals+1
       this%q_upperOrbital = totalActiveOrbitals
       this%r_lowerOrbital = otherCoreOrbitals+1
       this%r_upperOrbital = otherTotalOccupation + 1 !!occ + lumo
       this%s_lowerOrbital = otherCoreOrbitals+1
       this%s_upperOrbital = otherTotalActiveOrbitals
       symmetric = .true.

      if (CONTROL_instance%IONIZE_SPECIES(1) /= "NONE" ) then

        symmetric = .false.
        ionizeA = .false.
        ionizeB = .false.

         nameOfSpecies= trim(  MolecularSystem_getNameOfSpecies( speciesID ) )
         nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecies( otherSpeciesID ) )

         do s = 1, size(CONTROL_instance%IONIZE_SPECIES )
           if ( nameOfSpecies == trim(CONTROL_instance%IONIZE_SPECIES(s)) ) then
             ionizeA = .true. 
           end if
           if ( nameOfOtherSpecies == trim(CONTROL_instance%IONIZE_SPECIES(s)) ) then
             ionizeB = .true. 
           end if
         end do

        if ( CONTROL_instance%IONIZE_MO(1) == 0 ) then

          if ( ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = otherTotalOccupation + 1 !!occ + lumo
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation + 1 !! occ + lumo (default)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
  
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = othertotaloccupation 
            this%s_lowerOrbital = othertotaloccupation + 1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation
            this%q_lowerOrbital = totalOccupation + 1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = otherTotalOccupation + 1
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          end if

        else if ( CONTROL_instance%IONIZE_MO(1) /= 0 ) then !!occ and vir..

          if ( ionizeA .and. ionizeB ) then
            if ( CONTROL_instance%IONIZE_MO(1) <= totalOccupation .and. CONTROL_instance%IONIZE_MO(1) <= othertotalOccupation ) then
              this%p_lowerOrbital = coreOrbitals+1
              this%p_upperOrbital = totalOccupation
              this%q_lowerOrbital = coreOrbitals+1
              this%q_upperOrbital = totalActiveOrbitals
              this%r_lowerOrbital = otherCoreOrbitals+1
              this%r_upperOrbital = otherTotalOccupation
              this%s_lowerOrbital = otherCoreOrbitals+1
              this%s_upperOrbital = otherTotalActiveOrbitals

            else if ( CONTROL_instance%IONIZE_MO(1) > totalOccupation .and. CONTROL_instance%IONIZE_MO(1) > othertotalOccupation ) then

              this%p_lowerOrbital = coreOrbitals+1
              this%p_upperOrbital = CONTROL_instance%IONIZE_MO(1)
              this%q_lowerOrbital = coreOrbitals+1
              this%q_upperOrbital = totalActiveOrbitals
              this%r_lowerOrbital = otherCoreOrbitals+1
              this%r_upperOrbital =  CONTROL_instance%IONIZE_MO(1) 
              this%s_lowerOrbital = otherCoreOrbitals+1
              this%s_upperOrbital = otherTotalActiveOrbitals

            end if

          else if ( ionizeA .and. .not. ionizeB ) then
        
            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = max(CONTROL_instance%IONIZE_MO(1),totalOccupation)
            this%q_lowerOrbital = coreOrbitals+1
            this%q_upperOrbital = totalActiveOrbitals
  
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = othertotaloccupation 
            this%s_lowerOrbital = othertotaloccupation + 1
            this%s_upperOrbital = otherTotalActiveOrbitals

          else if ( .not. ionizeA .and. ionizeB ) then

            this%p_lowerOrbital = coreOrbitals+1
            this%p_upperOrbital = totalOccupation
            this%q_lowerOrbital = totalOccupation + 1
            this%q_upperOrbital = totalActiveOrbitals
            this%r_lowerOrbital = otherCoreOrbitals+1
            this%r_upperOrbital = max(CONTROL_instance%IONIZE_MO(1),otherTotalOccupation)
            this%s_lowerOrbital = otherCoreOrbitals+1
            this%s_upperOrbital = otherTotalActiveOrbitals

          end if

        end if

          !if ( CONTROL_instance%IONIZE_MO /= 0 ) then
      end if

   end if

    write(*,"(T15,A)") "Transformation boundaries "
    write(*,"(T15,A10,A6,A6)") "orbital","lower", "upper"
    write(*,"(T20,A5,I6,I6)") "p", this%p_lowerOrbital, this%p_upperOrbital
    write(*,"(T20,A5,I6,I6)") "q", this%q_lowerOrbital, this%q_upperOrbital
    write(*,"(T20,A5,I6,I6)") "r", this%r_lowerOrbital, this%r_upperOrbital
    write(*,"(T20,A5,I6,I6)") "s", this%s_lowerOrbital, this%s_upperOrbital
    print *, ""

   
  end subroutine TransformIntegralsC_checkInterMOIntegralType


  subroutine TransformIntegralsC_getNumberOfNonZeroRepulsionIntegrals( specieID, nproc, nonZeroIntegrals )
    implicit none
    integer :: specieID, nproc, auxNonZeroIntegrals
    integer :: nonZeroIntegrals
    integer :: ifile, unit
    character(50) :: sfile
    character(30) :: nameOfSpecie

    nonZeroIntegrals = 0

    do ifile = 1, nproc

       write(sfile,*) ifile
       sfile = trim(adjustl(sfile))
       unit = ifile+50

       nameOfSpecie = MolecularSystem_getNameOfSpecies( specieID )          

       if ( trim(nameOfSpecie) == "E-BETA" ) nameOfSpecie =""//trim("E-ALPHA")

       !! open file (order, integral(shell))
       open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".nints", status='old',access='sequential', form='Unformatted')

       read (unit) auxNonZeroIntegrals
       nonZeroIntegrals = nonZeroIntegrals + auxNonZeroIntegrals
       close (unit)
    end do

  end subroutine TransformIntegralsC_getNumberOfNonZeroRepulsionIntegrals

  subroutine TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals( i, j,  nproc, nonZeroIntegrals )
    implicit none
    integer :: i, j, nproc
    integer :: nonZeroIntegrals, auxNonZeroIntegrals
    integer :: ifile, unit
    character(50) :: sfile
    character(30) :: nameOfSpecie, nameOfOtherSpecie

    nonZeroIntegrals = 0

    do ifile = 1, nproc

       write(sfile,*) ifile
       sfile = trim(adjustl(sfile))
       unit = ifile+50

       nameOfSpecie = MolecularSystem_getNameOfSpecies( i )          
       nameOfOtherSpecie = MolecularSystem_getNameOfSpecies( j )          


       open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".nints", status='old',access='sequential', form='Unformatted')

       read (unit) auxNonZeroIntegrals
       nonZeroIntegrals = nonZeroIntegrals + auxNonZeroIntegrals
       close (unit)
    end do

    
  end subroutine TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine TransformIntegralsC_exception( typeMessage, description, debugDescription)
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

  end subroutine TransformIntegralsC_exception

end module TransformIntegralsC_
