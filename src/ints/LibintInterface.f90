!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Libint 1.04 interface
!!
!! @author E. F. Posada
!!
!! <b> Creation data : </b> 10-07-10
!!
!! <b> History change: </b>
!!
!!   - <tt> 10-07-10 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# Modulo que sirve como interfaz para las librerias LIBINT
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
!!
module LibintInterface_
  use Exception_
  use LibintInterfaceTypes_
  use MolecularSystem_
  use OverlapIntegrals_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use, intrinsic :: iso_c_binding
  implicit none
  
#define contr(n,m) contractions(n)%contractions(m)
  
  !> @brief Type for libint/c++ library
  type, public :: LibintInterface
     character(5) :: job
     logical :: isInstanced
     integer :: maxAngularMoment
     integer :: numberOfPrimitives
     integer :: libintStorage
     type(lib_int) :: libint
  end type LibintInterface

  !> @brief the integrals are saved for big records (that reduces the I/O time)
  type, public :: erisStack
     integer*2, allocatable :: a(:)
     integer*2, allocatable :: b(:)
     integer*2, allocatable :: c(:)
     integer*2, allocatable :: d(:)
     real(8), allocatable :: integrals(:)
  end type erisStack

  !> @brief type for convenence
  type, private :: auxBasis
     type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie
  end type auxBasis
  
  !>
  !! pointers to array functions on  libint.a, libderiv.a y libr12.a
  type(c_funptr), dimension(0:4,0:4,0:4,0:4), bind(c) :: build_eri
  
  interface
     
     !>
     !!Interfaz a libint.a
     subroutine LibintInterface_initLibIntBase() bind(C, name="init_libint_base")
       implicit none
       
     end subroutine LibintInterface_initLibIntBase
     
     function LibintInterface_initLibInt(libInt, maxAngMoment, numberOfPrimitives) bind(C, name="init_libint")
       use LibintInterfaceTypes_
       use, intrinsic :: iso_c_binding
       implicit none

       integer(kind=c_int) :: LibintInterface_initLibInt
       type(lib_int) :: libInt
       integer(kind=c_int), value :: maxAngMoment
       integer(kind=c_int), value :: numberOfPrimitives
       
     end function LibintInterface_initLibInt
     
     function LibintInterface_buildLibInt(libInt, numberOfPrimitives) bind(C)
       use LibintInterfaceTypes_
       use, intrinsic :: iso_c_binding
       implicit none
       
       type(c_ptr) :: LibintInterface_buildLibint
       type(lib_int) :: libInt
       integer(kind=c_int), value :: numberOfPrimitives
       
     end function LibintInterface_buildLibInt
     
     
     subroutine LibintInterface_freeLibInt(libInt) bind(C, name="free_libint")
       use LibintInterfaceTypes_
       use, intrinsic :: iso_c_binding
       implicit none
       
       type(lib_int) :: libInt
       
     end subroutine LibintInterface_freeLibInt
     
  end interface
  
  !> @brief Singleton lock
  type(LibintInterface), public :: LibintInterface_instance

  !> @brief Integrals Stack
  type(erisStack), private :: eris

contains
  
  !>
  !! @brief Constructor por omision, se llama una sola vez!
  !! @author E. F. Posada, 2010
  !! @param job, trabajo a realizar  
  !! @version 2.0
  subroutine LibintInterface_constructor( maxAngMoment, numberOfPrimitives, job )
    implicit none
    
    integer, intent(in) :: maxAngMoment
    integer, intent(in) :: numberOfPrimitives
    character(*) :: job
    
    integer(kind=c_int) :: libIntStorage_c
    integer(kind=c_int) :: maxAngMoment_c
    integer(kind=c_int) :: numberOfPrimitives_c

    if (.not. LibintInterface_isInstanced()) then

       LibintInterface_instance%job = trim(job)
       LibintInterface_instance%maxAngularMoment = maxAngMoment
       LibintInterface_instance%numberOfPrimitives = numberOfPrimitives


       select case(trim(LibintInterface_instance%job))

       case("ERIS")

          call LibintInterface_initLibIntBase()
          maxAngMoment_c = LibintInterface_instance%maxAngularMoment
          numberOfPrimitives_c = LibintInterface_instance%numberOfPrimitives
          libIntStorage_c = LibintInterface_initLibInt(LibintInterface_instance%libint, maxAngMoment_c, numberOfPrimitives_c)
          LibintInterface_instance%libintStorage = libIntStorage_c
          
       end select

       allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))

       LibintInterface_instance%isInstanced = .true.
       
    else
       
       call LibintInterface_exception( ERROR, "in libint interface constructor function",&
            "you must destroy this object before to use it again")
       
    end if

  end subroutine LibintInterface_constructor

  !>
  !! @brief Destructor por omision para ser llamado una sola vez
  !! @author E. F. Posada, 2010
  !! @version 2.0
  subroutine LibintInterface_destructor()
    implicit none
    
    if (LibintInterface_isInstanced()) then
       
       select case(trim(LibintInterface_instance%job))
          
       case("ERIS")
          
          !call LibintInterface_freeLibInt(LibintInterface_instance%libint) !! it does not work
          
       end select
       
       deallocate(eris%a, eris%b, eris%c, eris%d, eris%integrals)
       
       LibintInterface_instance%isInstanced = .false.
       
    else
       
       call LibintInterface_exception( ERROR, "in libint interface destructor function", "you must instantiate this object before to destroying it")
       
    end if
    
  end subroutine LibintInterface_destructor

  !>
  !! @brief Muestra informacion del objeto (una sola vez)
  !! @author E. F. Posada, 2010
  !! @version 2.0
  subroutine LibintInterface_show()
    implicit none

    write(*,  "(A)")  " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
    write(*, "(A)")   " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2011   "
    write(*, "(A)")   " ----------------------------------------------------------------------"
    write(*, "(A)" )  " LIBINT parameters: "
    !write(*, "(A, T44,A)")  " work ", trim(LibintInterface_instance%job)
    !write(*, "(A, T44,A)")  " Storage ", "DISK"
    write(*, "(A, T43,I5)") " Stack size", CONTROL_instance%INTEGRAL_STACK_SIZE
    !write(*, "(A, T43,I5)") " maxAngularMoment", LibintInterface_instance%maxAngularMoment
    write(*, "(A, T38,I10)") " number of primitives", LibintInterface_instance%numberOfPrimitives
    write(*, "(A, T38,I10/)")" Memory required (in words)", LibintInterface_instance%libintStorage

  end subroutine LibintInterface_show


  !>
  !! @brief calculate eris using libint library for all basis set (intra-specie)
  !! @author E. F. Posada, 2010
  !! @version 2.0
  !! @info Tested
  subroutine LibintInterface_computeIntraSpecies(specieID, job, starting, ending, process)
    implicit none
    
    character(*), intent(in) :: job
    integer, intent(in) :: specieID
    integer(8), intent(in) :: starting
    integer(8), intent(in) :: ending
    integer, intent(in) :: process
    
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    integer :: numberOfPrimitives
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: arraySize !< number of cartesian orbitals for maxAngularMoment
    integer :: n,u,m !< auxiliary itetators
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT)
    integer :: a, b, r, s !< not permuted iterators (original)
    integer :: pa, pb, pr, ps !< labels index
    integer :: apa, apb, apr, aps !< labels index
    integer :: ii, jj, kk, ll !< cartesian iterators for primitives and contractions
    integer :: aux, order !<auxiliary index
    integer :: arraySsize(1)
    integer :: counter, auxCounter

    integer(8) :: control

    integer,target :: i, j, k, l !< contraction length iterators
    integer,pointer :: pi, pj, pk, pl !< pointer to contraction length iterators
    integer,pointer :: poi, poj, pok, pol !< pointer to contraction length iterators

    integer, allocatable :: labelsOfContractions(:) !< cartesian position of contractions in all basis set

    real(8), dimension(:), pointer :: integralsPtr !< pointer to C array of integrals
    real(8), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: auxIntegrals(:) !<array with permuted integrals aux!
    real(8), allocatable :: integralsValue (:) !<array with permuted integrals
    real(8), allocatable :: incompletGamma(:) !<array with incomplete gamma integrals

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !<geometric values that appear in gaussian product
    real(8) :: zeta, eta, rho !< exponents... that appear in gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !<geometric quantities that appear in gaussian product
    real(8) :: incompletGammaArgument    
    
    character(50) :: fileNumber
    integer(8) :: ssize

    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie
    type(prim_data), target :: primitiveQuartet !<Primquartet object needed by LIBINT
    type(c_ptr) :: resultPc !< array of integrals from C (LIBINT)
    procedure(LibintInterface_buildLibInt), pointer :: pBuild !<procedure to calculate eris on LIBINT
    
    integer :: contractionNumberdebug, primitiveCounterdebug, contractionNumberperOrbital, totalIntegralswithP

    write(fileNumber,*) process
    fileNumber = trim(adjustl(fileNumber))
    
    !! open file for integrals
    open(UNIT=34,FILE=trim(fileNumber)//trim(MolecularSystem_instance%species(specieID)%name)//".ints", &
         STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='Unformatted')

    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID)



    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

    ssize = (numberOfContractions * (numberOfContractions + 1))/2
    ssize = (ssize * (ssize + 1))/2

    !! Libint constructor (solo una vez)
    if( LibintInterface_isInstanced() ) then
       call LibintInterface_destructor()
    end if
    
    if( .not. LibintInterface_isInstanced() ) then
       call LibintInterface_constructor( maxAngularMoment, numberOfPrimitives, trim(job))
       !! DEBUG
       !call LibintInterface_show()
    end if

    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))

    !! Debug
    !write(*,*) "SSize: ", ssize
    !write(*,*) "Max Angular Moment: ", maxAngularMoment
    !write(*,*) "Number of Primitives: ", numberOfPrimitives
    !write(*,*) "Number of Contractions: ", numberOfContractions
    !write(*,*) "Total number of Contractions: ", numberOfContractions

    !!Real labels for contractions
    aux = 1
    do i = 1, numberOfContractions
       !!position for cartesian contractions
       labelsOfContractions(i) = aux
       ! write(*,*) "Label of Contractions: ", labelsOfContractions(i)
       aux = aux + contractions(i)%numCartesianOrbital          
    end do

    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2

    ! write(*,*) "Array size: ", arraySize

    if(allocated(incompletGamma)) deallocate(incompletGamma)
    if(allocated(auxIntegrals)) deallocate(auxIntegrals)
    if(allocated(integralsValue)) deallocate(integralsValue)
    
    allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize))
    allocate(integralsValue(arraySize* arraySize* arraySize* arraySize))
    allocate(incompletGamma(0:MaxAngularMoment*4))

    counter = 0
    auxCounter = 0    
    control = 0
    contractionNumberdebug = 0
    totalIntegralswithP = 0
   ! write(*,*) "Starting and ending", starting, ending
    !!Start Calculating integrals for each shell
    ! write(*,*) "Number of Contractions Libint: ", numberOfContractions
    do a = 1, numberOfContractions
       n = a
       do b = a, numberOfContractions
          u = b
          do r = n , numberOfContractions
             do s = u,  numberOfContractions
                                
                contractionNumberdebug = contractionNumberdebug + 1
                control = control + 1                 

                if( control >= starting ) then                    
                   
                   if (control > ending) then
                      
                      counter = counter + 1
                      eris%a(counter) = -1
                      eris%b(counter) = -1
                      eris%c(counter) = -1
                      eris%d(counter) = -1
                      eris%integrals(counter) = 0.0_8
                      
                      write(34) &
                           eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                           eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                           eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                           eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                           eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   
                      write(6,"(A,I12,A,A)") "**Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
                           trim(MolecularSystem_instance%species(specieID)%name)
                      
                      ! close(69)
                      close(34)

                      return
                      
                   end if
                
                   ! Calcula el momento angular total
                   sumAngularMoment =  contractions(a)%angularMoment + &
                        contractions(b)%angularMoment + &
                        contractions(r)%angularMoment + &
                        contractions(s)%angularMoment
                   
                   !!Calcula el tamano del arreglo de integrales para la capa (ab|rs)
                   arraySize = contractions(a)%numCartesianOrbital * &
                        contractions(b)%numCartesianOrbital * &
                        contractions(r)%numCartesianOrbital * &
                        contractions(s)%numCartesianOrbital

                   !!For (ab|rs)  ---> RESTRICTION a>b && r>s && r+s > a+b
                   aux = 0
                   order = 0
                   
                   !!permuted index
                   aa = a
                   bb = b
                   rr = r
                   ss = s
                   
                   !!pointer to permuted index under a not permuted loop
                   pi => i
                   pj => j
                   pk => k
                   pl => l
                   
                   !!pointer to not permuted index under a permuted loop
                   poi => i
                   poj => j
                   pok => k
                   pol => l
                   
                   ! write(*,*) contractions(a)%angularMoment, contractions(b)%angularMoment, contractions(r)%angularMoment, contractions(s)%angularMoment
                   
                   if (contractions(a)%angularMoment < contractions(b)%angularMoment) then
                      
                      aa = b
                      bb = a
                      
                      pi => j
                      pj => i
                      
                      poi => j
                      poj => i
                      
                      order = order + 1
                      
                   end if
                   
                   if (contractions(r)%angularMoment < contractions(s)%angularMoment) then
                      
                      rr = s
                      ss = r
                      
                      pk => l
                      pl => k
                      
                      pok => l
                      pol => k
                      
                      order = order + 3
                      
                   end if
                   
                   if((contractions(a)%angularMoment + contractions(b)%angularMoment) > &
                        (contractions(r)%angularMoment + contractions(s)%angularMoment)) then
                      
                      aux = aa
                      aa = rr
                      rr = aux
                      
                      aux = bb
                      bb = ss
                      ss = aux
                      
                      select case(order)
                      case(0)
                         pi => k
                         pj => l
                         pk => i
                         pl => j
                         
                         poi => k
                         poj => l
                         pok => i
                         pol => j
                         
                      case(1)
                         pi => k
                         pj => l
                         pk => j
                         pl => i
                         
                         poi => l
                         poj => k
                         pok => i
                         pol => j
                         
                      case(3)
                         pi => l
                         pj => k
                         pk => i
                         pl => j
                         
                         poi => k
                         poj => l
                         pok => j
                         pol => i
                         
                      case(4)
                         pi => l
                         pj => k
                         pk => j
                         pl => i
                         
                         poi => l
                         poj => k
                         pok => j
                         pol => i
                         
                      end select
                      
                   end if
                   
                   !!************************************
                   !! Calculate iteratively primitives
                   !!
                   
                   !!Distancias AB, CD
                   AB = contractions(aa)%origin - contractions(bb)%origin
                   CD = contractions(rr)%origin - contractions(ss)%origin
                   
                   AB2 = dot_product(AB, AB)
                   CD2 = dot_product(CD, CD)
                   
                   !!Asigna valores a la estrucutra Libint
                   LibintInterface_instance%libint%AB = AB
                   LibintInterface_instance%libint%CD = CD
                   
                   !!start :)
                   integralsValue(1:arraySize) = 0.0_8
                   
                   ! write(*,*)"Pointer printing contraida"
                   ! write(*,*) a,b,r,s
                   ! write(*,*) aa,bb,rr,ss
                   
                   ! write(*,*) "**************************************************************************"
                   ! write(*,"(A20,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,4(I3))") "Contraida numero: ", contractionNumberdebug, &
                   !      " (", a, ",", b, "|", r, ",", s, ") ", &
                   !      contractions(s)%length, &
                   !      contractions(r)%length, &
                   !      contractions(b)%length, &
                   !      contractions(a)%length
                   ! write(*,*) "**************************************************************************"

                   primitiveCounterdebug = 0
                   do l = 1, contractions(s)%length
                      do k = 1, contractions(r)%length
                         do j = 1, contractions(b)%length
                            do i = 1, contractions(a)%length
                               
                               primitiveCounterdebug = primitiveCounterdebug + 1
                               !!LIBINT PRIMQUARTET
                               
                               !!Exponentes
                               zeta = contractions(aa)%orbitalExponents(pi) + &
                                    contractions(bb)%orbitalExponents(pj)
                               
                               eta =  contractions(rr)%orbitalExponents(pk) + &
                                    contractions(ss)%orbitalExponents(pl)
                               
                               rho  = (zeta * eta) / (zeta + eta) !Exponente reducido ABCD
                               
                               !!prim_data.U
                               P  = (( contractions(aa)%orbitalExponents(pi) * contractions(aa)%origin ) + &
                                    ( contractions(bb)%orbitalExponents(pj) * contractions(bb)%origin )) / zeta
                               
                               Q  = (( contractions(rr)%orbitalExponents(pk) * contractions(rr)%origin ) + &
                                    ( contractions(ss)%orbitalExponents(pl) * contractions(ss)%origin )) / eta
                               
                               W  = ((zeta * P) + (eta * Q)) / (zeta + eta)
                               
                               PQ = P - Q
                               
                               primitiveQuartet%U(1:3,1)= (P - contractions(aa)%origin)
                               primitiveQuartet%U(1:3,3)= (Q - contractions(rr)%origin)
                               primitiveQuartet%U(1:3,5)= (W - P)
                               primitiveQuartet%U(1:3,6)= (W - Q)
                               
                               !!Distancias ABCD(PQ2)
                               PQ2 = dot_product(PQ, PQ)
                               
                               !!Evalua el argumento de la funcion gamma incompleta
                               incompletGammaArgument = rho*PQ2
                               
                               !!Overlap Factor
                               s12 = ((Math_PI/zeta)**1.5_8) * exp(-(( contractions(aa)%orbitalExponents(pi) * &
                                    contractions(bb)%orbitalExponents(pj)) / zeta) * AB2)
                               
                               s34 = ((Math_PI/ eta)**1.5_8) * exp(-(( contractions(rr)%orbitalExponents(pk) * &
                                    contractions(ss)%orbitalExponents(pl)) /  eta) * CD2)
                               
                               s1234 = sqrt(rho/Math_PI) * s12 * s34
                               
                               !!Screening
                               if(abs(s1234) < 1.0D-12) cycle
                               
                               call Math_fgamma0(sumAngularMoment,incompletGammaArgument,incompletGamma(0:sumAngularMoment))
                               
                               
                               !!prim_data.F
                               primitiveQuartet%F(1:sumAngularMoment+1) = 2.0_8 * incompletGamma(0:sumAngularMoment) * s1234

                               !!Datos restantes para prim.U
                               primitiveQuartet%oo2z = (0.5_8 / zeta)
                               primitiveQuartet%oo2n = (0.5_8 / eta)
                               primitiveQuartet%oo2zn = (0.5_8 / (zeta + eta))
                               primitiveQuartet%poz = (rho/zeta)
                               primitiveQuartet%pon = (rho/eta)
                               primitiveQuartet%oo2p = (0.5_8 / rho)
                               
                               if(arraySize == 1) then
                                  ! write(*,*)"Pointer printing (ss|ss)"
                                  ! write(*,*)pi,pj,pk,pl
                                  auxIntegrals(1) = primitiveQuartet%F(1)
                                  ! write(*,*) auxIntegrals(1), s1234,  incompletGamma(0:sumAngularMoment), CD2

                                  ! write(*,"(A20,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,<sumAngularMoment+1>(F16.12), A3, <sumAngularMoment+1>(F16.12), A3, F16.12)") "Primitiva numero: ", &
                                  !      primitiveCounterdebug, &
                                  !      " (", pi, ",", pj, "|", pk, ",", pl, ") ", &
                                  !      incompletGamma(0:sumAngularMoment), &
                                  !      " | ", &
                                  !      primitiveQuartet%F(1:sumAngularMoment+1), &
                                  !      " | ", &
                                  !      auxIntegrals(1)

                               else
                                  write(*,*)"Pointer printing"
                                  write(*,*)aa,bb,rr,ss
                                  arraySsize(1) = arraySize
                                  
                                  LibintInterface_instance%libint%PrimQuartet = c_loc(primitiveQuartet)
                                  
                                  call c_f_procpointer(build_eri( contractions(ss)%angularMoment , &
                                       contractions(rr)%angularMoment , &
                                       contractions(bb)%angularMoment , &
                                       contractions(aa)%angularMoment), pBuild)
                                  
                                  resultPc = pBuild(LibintInterface_instance%libint,1)

                                  call c_f_pointer(resultPc, temporalPtr, arraySsize)
                                  
                                  integralsPtr => temporalPtr
                                  
                                  auxIntegrals(1:arraySize) = integralsPtr(1:arraySize)
                               end if !!done by contractions

                                                                                      
                               !!Normalize by primitive
                               m = 0
                               do ii = 1, contractions(aa)%numCartesianOrbital
                                  do jj = 1, contractions(bb)%numCartesianOrbital
                                     do kk = 1, contractions(rr)%numCartesianOrbital
                                        do ll = 1, contractions(ss)%numCartesianOrbital
                                           m = m + 1
                                           !if(aa == 1 .and. bb == 1 .and. rr == 3 .and. ss == 1) then
                                           !   write(*,*) "Integral sin Norm, K= ", kk
                                           !   write(*,*) "( ", aa, ", ", bb, " | ", rr, ", ", ss, " )"
                                           !   write(*,*) "( ", pi, ", ", pj, " | ", pk, ", ", pl, " )"
                                           !   write(*,*) auxIntegrals(m)
                                           !   write(*,*) contractions(rr)%primNormalization(pk,kk)
                                           !end if
                                           auxIntegrals(m) = auxIntegrals(m) &
                                                * contractions(aa)%primNormalization(pi,ii) &
                                                * contractions(bb)%primNormalization(pj,jj) &
                                                * contractions(rr)%primNormalization(pk,kk) &
                                                * contractions(ss)%primNormalization(pl,ll)
                                           ! if(aa == 1 .and. bb == 1 .and. rr == 3 .and. ss == 3) then
                                           !    write(*,*) "Integral Normalizada ( ", aa, ", ", bb, " | ", rr, ", ", ss, " )"
                                           !    write(*,*) "( ", pi, ", ", pj, " | ", pk, ", ", pl, " )", auxIntegrals(m)
                                           !    write(*,*) contractions(rr)%primNormalization(pk,kk)
                                           ! end if
                                        end do
                                     end do
                                  end do
                               end do !! done by cartesian of contractions                               
                               
                               auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
                                    * contractions(aa)%contractionCoefficients(pi) &
                                    * contractions(bb)%contractionCoefficients(pj) &
                                    * contractions(rr)%contractionCoefficients(pk) &
                                    * contractions(ss)%contractionCoefficients(pl)                                                             
                               ! write(*,*) "(",aa,",",bb,"|",rr,",",ss,")", ": ", &
                               !      "(",pi,",",pj,"|",pk,",",pl,")", ": ", &
                               !      contractions(aa)%contractionCoefficients(pi), &
                               !      contractions(bb)%contractionCoefficients(pj), &
                               !      contractions(rr)%contractionCoefficients(pk), &
                               !      contractions(ss)%contractionCoefficients(pl)
                               ! write(*,*) "( ", pi, ", ", pj, " | ", pk, ", ", pl, " )", auxIntegrals(1:arraySize)

                               integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)
                               
                            end do
                         end do
                      end do
                   end do !!done integral by shell                
                   

                   !!normalize by shell
                   m = 0
                   contractionNumberperOrbital = 0
                   ! write(*,*) "Con numero cartesian Orbital"
                   do ii = 1,  contractions(aa)%numCartesianOrbital
                      do jj = 1,  contractions(bb)%numCartesianOrbital
                         do kk = 1,  contractions(rr)%numCartesianOrbital
                            do ll = 1, contractions(ss)%numCartesianOrbital
                               m = m + 1
                               contractionNumberperOrbital = contractionNumberperOrbital + 1
                               integralsValue(m) = integralsValue(m) &
                                    * contractions(aa)%contNormalization(ii) &
                                    * contractions(bb)%contNormalization(jj) &
                                    * contractions(rr)%contNormalization(kk) &
                                    * contractions(ss)%contNormalization(ll)                            
                               ! if(aa == 1 .and. bb == 1 .and. rr == 3 .and. ss == 3) then
                               !    write(*,*) "Integral Contraida ( ", aa, ", ", bb, " | ", rr, ", ", ss, " )", integralsValue(m)
                               ! end if
                            end do
                         end do
                      end do
                   end do !! done by shell                                         
                   totalIntegralswithP = totalIntegralswithP + contractionNumberperOrbital
                   ! write(*,*) "Number of integrals : ", contractionNumberperOrbital, totalIntegralswithP
                   !!write to disk
                   m = 0
                   do i = 1, contractions(aa)%numCartesianOrbital
                      do j = 1, contractions(bb)%numCartesianOrbital
                         do k = 1, contractions(rr)%numCartesianOrbital
                            do l = 1, contractions(ss)%numCartesianOrbital
                               
                               m = m + 1

                               !! index not permuted
                               pa=labelsOfContractions(a)+poi-1
                               pb=labelsOfContractions(b)+poj-1
                               pr=labelsOfContractions(r)+pok-1
                               ps=labelsOfContractions(s)+pol-1
                               
                               apa=pa
                               apb=pb
                               apr=pr
                               aps=ps
                               
                               if( pa <= pb .and. pr <= ps .and. (pa*1000)+pb >= (pr*1000)+ps) then
                                  
                                  aux = pa
                                  pa = pr
                                  pr = aux
                                  
                                  aux = pb
                                  pb = ps
                                  ps = aux
                                  
                               end if
                               
                               if( pa <= pb .and. pr <= ps ) then
                                  
                                  if (  (apa /= pa .or. apb/=pb .or. apr/=pr .or. aps/=ps) .and. ( b==s ) ) then
                                     
                                     !! Descarted! (they are repeated)
                                  else

                                     write(*,"(A20,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2, A3, F16.12)") &
                                          "Contraida numero: ", &
                                          contractionNumberdebug, &
                                          " (", pa, ",", pb, "|", pr, ",", ps, ") ", &
                                          " | ", &
                                          integralsValue(m)
                                                                         
                                     if(abs(integralsValue(m)) > 1.0D-10) then

                                        counter = counter + 1
                                        auxCounter = auxCounter + 1
                                        
                                        eris%a(counter) = pa
                                        eris%b(counter) = pb
                                        eris%c(counter) = pr
                                        eris%d(counter) = ps
                                        eris%integrals(counter) = integralsValue(m)
                                     end if
                                     
                                     if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

                                        write(34) &
                                             eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                             eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                             eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                             eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                             eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                                        counter = 0
                                        
                                     end if
                                     
                                  end if
                               end if

                            end do
                         end do
                      end do
                   end do !! Primitives loop
                   
                end if

                ! write(*,*) "Info libintB"
                ! write(*,*) a, b, r, s, sumAngularMoment
                ! write(*,*) aa,bb,rr,ss, sumAngularMoment
             end do
             u=r+1
          end do
       end do
    end do !! shells loop

    counter = counter + 1
    eris%a(counter) = -1
    eris%b(counter) = -1
    eris%c(counter) = -1
    eris%d(counter) = -1
    eris%integrals(counter) = 0.0_8
       
    write(34) &
         eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
    
    ! close(69)
    ! write(*,*) &
    !      "Contraida numero: ", &
    !      contractionNumberdebug, &
    !      " (", eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), "," , &
    !      eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), "|", &
    !      eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), ",", &
    !      eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), ") ", &
    !      " | ", &
    !      eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    close(34)

    write(6,"(A,I12,A,A)") " Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
         trim(MolecularSystem_instance%species(specieID)%name)
    
  end subroutine LibintInterface_computeIntraSpecies
  
  !>
  !! @brief calculate eris using libint library for all basis set (inter-species)
  !! @author E. F. Posada, 2010
  !! @version 2.0
  subroutine LibintInterface_computeInterSpecies(specieID, otherSpecieID, job, isInterSpecies, isCouplingA, isCouplingB)
    implicit none

    integer,target :: specieID
    integer,target :: otherSpecieID
    character(*), intent(in) :: job
    logical, optional :: isInterSpecies
    logical, optional :: isCouplingA
    logical, optional :: isCouplingB
    
    logical :: interSpecies
    logical :: couplingA
    logical :: couplingB
    
    integer :: auxSpecieID
    integer :: auxOtherSpecieID    
    integer :: totalNumberOfContractions
    integer :: otherTotalNumberOfContractions
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: numberOfPrimitives
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: arraySize !! number of cartesian orbitals for maxAngularMoment
    integer :: n,u,m !! auxiliary itetators
    integer :: aa, bb, rr, ss !! permuted iterators (LIBINT)
    integer :: a, b, r, s !! not permuted iterators (original)
    integer*2 :: pa, pb, pr, ps !! labels index
    integer :: ii, jj, kk, ll !! cartesian iterators for contractions and contractions
    integer :: aux, order !!auxiliary index
    integer :: arraySsize(1)
    integer :: sizeTotal
    integer :: auxIndex
    integer :: counter, auxCounter
    
    integer,target :: i, j, k, l !! contraction length iterators
    integer,pointer :: pi, pj, pk, pl !! pointer to contraction length iterators
    integer,pointer :: poi, poj, pok, pol !! pointer to contraction length iterators
    integer, pointer :: pSpecieID, pOtherSpecieID !! pointer to species ID
    integer, allocatable :: labelsOfContractions(:) !! cartesian position of contractions in all basis set
    integer, allocatable :: otherLabelsOfContractions(:) !! cartesian position of contractions in all basis set

    real(8), dimension(:), pointer :: integralsPtr !! pointer to C array of integrals
    real(8), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: auxIntegrals(:) !!array with permuted integrals aux!
    real(8), allocatable :: integralsValue (:) !!array with permuted integrals
    real(8), allocatable :: incompletGamma(:) !!array with incomplete gamma integrals
    
    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !!geometric values that appear in gaussian product
    real(8) :: zeta, eta, rho !! exponents... that appear in gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !!geometric quantities that appear in gaussian product
    real(8) :: incompletGammaArgument
    
    type(auxBasis) :: contractions(2)
    type(prim_data), target :: primitiveQuartet !!Primquartet object needed by LIBINT
    type(c_ptr) :: resultPc !! array of integrals from C (LIBINT)
    procedure(LibintInterface_buildLibInt), pointer :: pBuild !!procedure to calculate eris on LIBINT
    
    interSpecies = .true.
    couplingA = .false.
    couplingB = .false.
    
    if(present(isInterSpecies)) interSpecies = isInterSpecies
    if(present(isCouplingA)) couplingA = isCouplingA
    if(present(isCouplingB)) couplingB = isCouplingB
    
    !! open file for integrals
    open(UNIT=34,FILE=trim(MolecularSystem_instance%species(specieID)%name)//"."//trim(MolecularSystem_instance%species(otherSpecieID)%name)//".ints", &
         STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='Unformatted')
    
    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions(1)%contractions)
    call MolecularSystem_getBasisSet(otherSpecieID, contractions(2)%contractions)
    
    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions(1)%contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

    !! Get number of shells and cartesian contractions (other specie)
    otherNumberOfContractions = size(contractions(2)%contractions)
    otherTotalNumberOfContractions = MolecularSystem_instance%species(otherSpecieID)%basisSetSize
    
    !! Libint constructor (solo una vez)
    maxAngularMoment = max(MolecularSystem_getMaxAngularMoment(specieID), MolecularSystem_getMaxAngularMoment(otherSpecieID))
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID) + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
    
    if( .not. LibintInterface_isInstanced() ) then
       call LibintInterface_constructor( maxAngularMoment, numberOfPrimitives, trim(job))
       call LibintInterface_show()
    end if
    
    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2
    
    sizeTotal = ((totalNumberOfContractions *(totalNumberOfContractions + 1 ))/2) * ((otherTotalNumberOfContractions *(otherTotalNumberOfContractions + 1 ))/2)
    
    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))
    
    !!Real labels for contractions
    aux = 1
    do i = 1, numberOfContractions
       labelsOfContractions(i) = aux
       aux = aux + contractions(1)%contractions(i)%numCartesianOrbital
    end do

    !! Get contractions labels for integrals index (other specie)
    if (allocated(otherLabelsOfContractions)) deallocate(otherLabelsOfContractions)
    allocate(otherLabelsOfContractions(otherNumberOfContractions))

    !!Real labels for contractions (other specie)
    aux = 1
    do i = 1, otherNumberOfContractions
       otherLabelsOfContractions(i) = aux
       aux = aux + contractions(2)%contractions(i)%numCartesianOrbital
    end do
    
    !! Allocating some space
    if(allocated(incompletGamma)) deallocate(incompletGamma)
    if(allocated(auxIntegrals)) deallocate(auxIntegrals)
    if(allocated(integralsValue)) deallocate(integralsValue)

    allocate(incompletGamma(0:MaxAngularMoment*4))    
    allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize))
    allocate(integralsValue(arraySize* arraySize* arraySize* arraySize))

    
    counter = 0
    auxCounter = 0
    
    auxSpecieID = specieID
    auxOtherSpecieID = otherSpecieID

    specieID = 1
    otherSpecieID = 2
    
    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       do b = a, numberOfContractions
          do r = 1 , otherNumberOfContractions
             do s = r,  otherNumberOfContractions
                
                !!Calcula el momento angular total
                sumAngularMoment =  contractions(1)%contractions(a)%angularMoment + &
                     contractions(1)%contractions(b)%angularMoment + &
                     contractions(2)%contractions(r)%angularMoment + &
                     contractions(2)%contractions(s)%angularMoment
                
                !! Calcula el tamano del arreglo de integrales para la capa (ab|rs)
                arraySize = contractions(1)%contractions(a)%numCartesianOrbital * &
                     contractions(1)%contractions(b)%numCartesianOrbital * &
                     contractions(2)%contractions(r)%numCartesianOrbital * &
                     contractions(2)%contractions(s)%numCartesianOrbital
                
                !! For (ab|rs)  ---> RESTRICTION a>b && r>s && r+s > a+b
                aux = 0
                order = 0

                !! permuted index
                aa = a
                bb = b
                rr = r
                ss = s
                
                !!pointer to permuted index under a not permuted loop
                pi => i
                pj => j
                pk => k
                pl => l
                
                !!pointer to not permuted index under a permuted loop
                poi => i
                poj => j
                pok => k
                pol => l

                !!Pointer to specie ID
                pSpecieID => specieID
                pOtherSpecieID => otherSpecieID
                
                if (contractions(1)%contractions(a)%angularMoment < contractions(1)%contractions(b)%angularMoment) then

                   aa = b
                   bb = a
                   
                   pi => j
                   pj => i

                   poi => j
                   poj => i

                   order = order + 1
                end if
                
                if (contractions(2)%contractions(r)%angularMoment < contractions(2)%contractions(s)%angularMoment) then
                   
                   rr = s
                   ss = r

                   pk => l
                   pl => k

                   pok => l
                   pol => k

                   order = order + 3

                end if

                if((contractions(1)%contractions(a)%angularMoment + contractions(1)%contractions(b)%angularMoment) > &
                     (contractions(2)%contractions(r)%angularMoment + contractions(2)%contractions(s)%angularMoment)) then
                   
                   aux = aa
                   aa = rr
                   rr = aux

                   aux = bb
                   bb = ss
                   ss = aux

                   pSpecieID => otherSpecieID
                   pOtherSpecieID => specieID

                   select case(order)
                   case(0)
                      pi => k
                      pj => l
                      pk => i
                      pl => j

                      poi => k
                      poj => l
                      pok => i
                      pol => j

                   case(1)
                      pi => k
                      pj => l
                      pk => j
                      pl => i

                      poi => l
                      poj => k
                      pok => i
                      pol => j

                   case(3)
                      pi => l
                      pj => k
                      pk => i
                      pl => j

                      poi => k
                      poj => l
                      pok => j
                      pol => i

                   case(4)
                      pi => l
                      pj => k
                      pk => j
                      pl => i

                      poi => l
                      poj => k
                      pok => j
                      pol => i

                   end select

                   order = order + 5

                end if

                !!************************************
                !! Calculate iteratively contractions
                !!

                !!Distancias AB, CD
                AB = contr(pSpecieID,aa)%origin - contr(pSpecieID,bb)%origin
                CD = contr(pOtherSpecieID,rr)%origin - contr(pOtherSpecieID,ss)%origin

                AB2 = dot_product(AB, AB)
                CD2 = dot_product(CD, CD)

                !! Asigna valores a la estrucutra Libint
                LibintInterface_instance%libint%AB = AB
                LibintInterface_instance%libint%CD = CD

                !!start :)
                integralsValue(1:arraySize) = 0.0_8

                !! not-permuted loop
                do l = 1, contractions(otherSpecieID)%contractions(s)%length
                   do k = 1, contractions(otherSpecieID)%contractions(r)%length
                      do j = 1, contractions(specieID)%contractions(b)%length
                         do i = 1, contractions(specieID)%contractions(a)%length

                            !!LIBINT PRIMQUARTET

                            !!Exponentes
                            
                            zeta = contr(pSpecieID,aa)%orbitalExponents(pi) + &
                                 contr(pSpecieID,bb)%orbitalExponents(pj)

                            eta =  contr(pOtherSpecieID,rr)%orbitalExponents(pk) + &
                                 contr(pOtherSpecieID,ss)%orbitalExponents(pl)

                            rho  = (zeta * eta) / (zeta + eta) !Exponente reducido ABCD

                            !!prim_data.U
                            P  = (( contr(pSpecieID,aa)%orbitalExponents(pi) * contr(pSpecieID,aa)%origin ) + &
                                 ( contr(pSpecieID,bb)%orbitalExponents(pj) * contr(pSpecieID,bb)%origin )) / zeta

                            Q  = (( contr(pOtherSpecieID,rr)%orbitalExponents(pk) * contr(pOtherSpecieID,rr)%origin ) + &
                                 ( contr(pOtherSpecieID,ss)%orbitalExponents(pl) * contr(pOtherSpecieID,ss)%origin )) / eta
                            
                            W  = ((zeta * P) + (eta * Q)) / (zeta + eta)
                            
                            PQ = P - Q
                            
                            primitiveQuartet%U(1:3,1)= (P - contr(pSpecieID,aa)%origin)
                            primitiveQuartet%U(1:3,3)= (Q - contr(pOtherSpecieID,rr)%origin)
                            primitiveQuartet%U(1:3,5)= (W - P)
                            primitiveQuartet%U(1:3,6)= (W - Q)

                            !!Distancias ABCD(PQ2)
                            PQ2 = dot_product(PQ, PQ)

                            !!Evalua el argumento de la funcion gamma incompleta
                            incompletGammaArgument = rho*PQ2

                            !!Overlap Factor
                            s12 = ((Math_PI/zeta)**1.5_8) * exp(-(( contr(pSpecieID,aa)%orbitalExponents(pi) * &
                                 contr(pSpecieID,bb)%orbitalExponents(pj)) / zeta) * AB2)

                            s34 = ((Math_PI/ eta)**1.5_8) * exp(-(( contr(pOtherSpecieID,rr)%orbitalExponents(pk) * &
                                 contr(pOtherSpecieID,ss)%orbitalExponents(pl)) /  eta) * CD2)

                            s1234 = sqrt(rho/Math_PI) * s12 * s34

                            call Math_fgamma0(sumAngularMoment,incompletGammaArgument,incompletGamma(0:sumAngularMoment))

                            !!prim_data.F
                            primitiveQuartet%F(1:sumAngularMoment+1) = 2.0_8 * incompletGamma(0:sumAngularMoment) * s1234

                            !!Datos restantes para prim.U
                            primitiveQuartet%oo2z = (0.5_8 / zeta)
                            primitiveQuartet%oo2n = (0.5_8 / eta)
                            primitiveQuartet%oo2zn = (0.5_8 / (zeta + eta))
                            primitiveQuartet%poz = (rho/zeta)
                            primitiveQuartet%pon = (rho/eta)
                            primitiveQuartet%oo2p = (0.5_8 / rho)

                            if(arraySize == 1) then
                               
                               auxIntegrals(1) = primitiveQuartet%F(1)

                            else


                               arraySsize(1) = arraySize
                               
                               LibintInterface_instance%libint%PrimQuartet = c_loc(primitiveQuartet)

                               !! calculate integrals (finally)                               
                               call c_f_procpointer(build_eri( contr(pOtherSpecieID,ss)%angularMoment , &
                                    contr(pOtherSpecieID,rr)%angularMoment , &
                                    contr(pSpecieID,bb)%angularMoment , &
                                    contr(pSpecieID,aa)%angularMoment), pBuild)
                               
                               resultPc = pBuild(LibintInterface_instance%libint,1)

                               !! Get numbers from memory!
                               call c_f_pointer(resultPc, temporalPtr, arraySsize)
                               
                               integralsPtr => temporalPtr
                               
                               !! Copy to fortran casting
                               auxIntegrals(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy

                            end if !!done by contractions

                            !!Normalize by primitive
                            m = 0
                            do ii = 1, contr(pSpecieID,aa)%numCartesianOrbital
                               do jj = 1, contr(pSpecieID,bb)%numCartesianOrbital
                                  do kk = 1, contr(pOtherSpecieID,rr)%numCartesianOrbital
                                     do ll = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital
                                        m = m + 1
                                        auxIntegrals(m) = auxIntegrals(m) &
                                             * contr(pSpecieID,aa)%primNormalization(pi,ii) &
                                             * contr(pSpecieID,bb)%primNormalization(pj,jj) &
                                             * contr(pOtherSpecieID,rr)%primNormalization(pk,kk) &
                                             * contr(pOtherSpecieID,ss)%primNormalization(pl,ll)
                                     end do
                                  end do
                               end do
                            end do !! done by primitives
                            
                            auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
                                 * contr(pSpecieID,aa)%contractionCoefficients(pi) &
                                 * contr(pSpecieID,bb)%contractionCoefficients(pj) &
                                 * contr(pOtherSpecieID,rr)%contractionCoefficients(pk) &
                                 * contr(pOtherSpecieID,ss)%contractionCoefficients(pl)

                            integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)

                         end do
                      end do
                   end do
                end do !!done Integral of pair of shells

                !!normalize by contraction
                m = 0
                do ii = 1,  contr(pSpecieID,aa)%numCartesianOrbital
                   do jj = 1,  contr(pSpecieID,bb)%numCartesianOrbital
                      do kk = 1,  contr(pOtherSpecieID,rr)%numCartesianOrbital
                         do ll = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital
                            m = m + 1

                            integralsValue(m) = integralsValue(m) &
                                 * contr(pSpecieID,aa)%contNormalization(ii) &
                                 * contr(pSpecieID,bb)%contNormalization(jj) &
                                 * contr(pOtherSpecieID,rr)%contNormalization(kk) &
                                 * contr(pOtherSpecieID,ss)%contNormalization(ll)
                         end do
                      end do
                   end do
                end do !! done by cartesian of contractions
                
                !!write to disk
                m = 0
                do i = 1, contr(pSpecieID,aa)%numCartesianOrbital
                   do j = 1, contr(pSpecieID,bb)%numCartesianOrbital
                      do k = 1, contr(pOtherSpecieID,rr)%numCartesianOrbital
                         do l = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital

                            m = m + 1

                            !! index not permuted
                            pa=labelsOfContractions(a)+poi-1
                            pb=labelsOfContractions(b)+poj-1
                            pr=otherLabelsOfContractions(r)+pok-1
                            ps=otherLabelsOfContractions(s)+pol-1

                            if( pa <= pb .and. pr <= ps ) then

                               if(abs(integralsValue(m)) > 1.0D-10) then
                                  
                                  counter = counter + 1
                                  auxCounter = auxCounter + 1
                                  
                                  eris%a(counter) = pa
                                  eris%b(counter) = pb
                                  eris%c(counter) = pr
                                  eris%d(counter) = ps
                                  eris%integrals(counter) = integralsValue(m)

                               end if

                               
                               if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

                                  write(34) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
                                  
                                  counter = 0

                               end if
                            end if !! Stack control

                         end do
                      end do
                   end do
                end do !! Done write to disk
                
             end do
             u=r+1
          end do
       end do
    end do !! done by basis set

    eris%a(counter+1) = -1
    eris%b(counter+1) = -1
    eris%c(counter+1) = -1
    eris%d(counter+1) = -1
    eris%integrals(counter+1) = 0.0_8	       
    
    write(34) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    specieID = auxSpecieID
    otherSpecieID = auxOtherSpecieID
    
    write(*,"(A,I12,A,A)") " Stored ", &
         auxCounter, &
         " non-zero repulsion integrals between species: ", &
         trim(MolecularSystem_instance%species(specieID)%name)//" / "//&
         trim(MolecularSystem_instance%species(otherSpecieID)%name)

    close(34)
    
  end subroutine LibintInterface_computeInterSpecies

  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  function LibintInterface_isInstanced( ) result( output )
    implicit  none

    logical :: output
    
    output = LibintInterface_instance%isInstanced
    
  end function LibintInterface_isInstanced  

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine LibintInterface_exception( typeMessage, description, debugDescription)
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

  end subroutine LibintInterface_exception

end module LibintInterface_
