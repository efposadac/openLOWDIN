!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!      http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!      http://www.cucei.udg.mx/~robertof
!!
!!      Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief G12 integrals calculation using Libint 2 library
!!
!! @author E. F. Posada
!!
!! <b> Creation data : </b> 20-10-2016
!!
!! <b> History change: </b>
!!
!!   - <tt> 28-04-16 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# Creation of the module.

module G12Integrals_
  use iso_c_binding
  use Libint2Interface_
  use CONTROL_
  use MolecularSystem_
  use ContractedGaussian_
  use InterPotential_
  implicit none

#define contr(n,m) contractions(n)%contractions(m)

  type, public :: G12Integrals
     character(5) :: job
     logical :: isInstanced
     integer :: maxAngularMoment
     integer :: numberOfPrimitives
     integer :: libintStorage
     type(lib_int) :: libint
     type(libint2) :: libintG12
  end type G12Integrals

  type, public :: erisStack
     integer, allocatable :: a(:)
     integer, allocatable :: b(:)
     integer, allocatable :: c(:)
     integer, allocatable :: d(:)
     real(8), allocatable :: integrals(:)
  end type erisStack

  !> @brief type for convenience
  type, private :: auxBasis
     type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie
  end type auxBasis

  interface

     subroutine LibintInterface_libint2StaticInit() bind(C, name="libint2_static_init")
       implicit none

     end subroutine LibintInterface_libint2StaticInit

     function LibintInterface_buildLibInt(libInt, numberOfPrimitives) bind(C)
       use Libint2Interface_
       use, intrinsic :: iso_c_binding
       implicit none

       type(c_ptr) :: LibintInterface_buildLibint
       type(lib_int) :: libInt
       integer(kind=c_int), value :: numberOfPrimitives

     end function LibintInterface_buildLibInt

     subroutine LibintInterface_libint2StaticCleanup() bind(C, name="libint2_static_cleanup")
       implicit none

     end subroutine LibintInterface_libint2StaticCleanup

     subroutine LibintInterface_libint2InitR12kg12(libIntG12, maxAngMoment, buf) bind(C, name="libint2_init_eri")
       use Libint2Interface_
       use, intrinsic :: iso_c_binding
       implicit none

       type(libint2) :: libIntG12
       integer(kind=c_int), value :: maxAngMoment
       real(kind=c_double) :: buf

     end subroutine LibintInterface_libint2InitR12kg12

     function LibintInterface_libint2NeedMemoryR12kg12(maxAngMoment) bind(C, name="libint2_need_memory_eri")
       use, intrinsic :: iso_c_binding
       implicit none

       integer(kind=c_int) :: LibintInterface_libint2NeedMemoryR12kg12
       integer(kind=c_int), value :: maxAngMoment

     end function LibintInterface_libint2NeedMemoryR12kg12

     subroutine LibintInterface_libint2CleanupR12kg12(LibIntG12) bind(C, name="libint2_cleanup_eri")
       use Libint2Interface_
       use, intrinsic :: iso_c_binding
       implicit none

       type(libint2) :: libIntG12

     end subroutine LibintInterface_libint2CleanupR12kg12

  end interface

  type(G12Integrals), public, target :: G12Integrals_instance
  type(erisStack), private :: eris

contains

  !>
  !! @brief Constructor by default
  !! @param job
  subroutine G12Integrals_constructor( maxAngMoment, numberOfPrimitives )
    implicit none

    integer, intent(in) :: maxAngMoment
    integer, intent(in) :: numberOfPrimitives

    integer(kind=c_int) :: libIntStorage_c
    integer(kind=c_int) :: maxAngMoment_c
    integer(kind=c_int) :: numberOfPrimitives_c

    if (.not. G12Integrals_isInstanced()) then

       G12Integrals_instance%maxAngularMoment = maxAngMoment
       G12Integrals_instance%numberOfPrimitives = numberOfPrimitives

       call LibintInterface_libint2StaticInit()
       maxAngMoment_c = G12Integrals_instance%maxAngularMoment
       numberOfPrimitives_c = G12Integrals_instance%numberOfPrimitives

       call LibintInterface_libint2InitR12kg12(G12Integrals_instance%libintG12, maxAngMoment_c, 0.0_8)
       libIntStorage_c = LibintInterface_libint2NeedMemoryR12kg12(maxAngMoment_c)
       G12Integrals_instance%libintStorage = libIntStorage_c

       G12Integrals_instance%isInstanced = .true.

    end if

  end subroutine G12Integrals_constructor



  !<
  !! @brief calculate G12 eris using libint library for all basis set (intra-specie)
  subroutine G12Integrals_diskIntraSpecie(specieID)
    implicit none

    integer, intent(in) :: specieID

    character(150) :: nameOfSpecie
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    integer :: numberOfPrimitives
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: arraySize !! number of Cartesian orbitals for maxAngularMoment
    integer :: n,u,m !! auxiliary iterators
    integer :: aa, bb, rr, ss !! permuted iterators (LIBINT)
    integer :: a, b, r, s !! not permuted iterators (original)
    integer :: pa, pb, pr, ps !! labels index
    integer :: apa, apb, apr, aps !! labels index

    integer :: ii, jj, kk, ll !! cartesian iterators for primitives and contractions
    integer :: aux, order !!auxiliary index
    integer :: arraySsize(1)
    integer :: sizeTotal !!For Large systems

    integer :: counter, auxCounter

    integer,target :: i, j, k, l !! contraction length iterators
    integer,pointer :: pi, pj, pk, pl !! pointer to contraction length iterators
    integer,pointer :: poi, poj, pok, pol !! pointer to contraction length iterators
    integer, allocatable :: labelsOfContractions(:) !! Cartesian position of contractions in all basis set

    integer :: potID, potSize, potLength

    real(8), dimension(:), pointer :: integralsPtr !! pointer to C array of integrals
    real(8), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: auxIntegrals(:) !!array with permuted integrals aux!
    real(8), allocatable :: auxIntegralsB(:) !!array with permuted integrals aux!
    real(8), allocatable :: integralsValue (:) !!array with permuted integrals
    real(8), allocatable :: incompletGamma(:) !!array with incomplete gamma integrals

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3), AC(3), BD(3)!!geometric values that appear in Gaussian product
    real(8) :: zeta, eta, rho !! exponents... that appear in Gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !!geometric quantities that appear in Gaussian product
    real(8) :: incompletGammaArgument, ssss(21)
    real(8) :: prefactor

    real(8) :: startTime, endTime

    type(libint2), pointer :: G12_ptr
    type(ContractedGaussian), pointer :: contractionG12
    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie
    type(libint2), target :: primitiveQuartet !!Prim-quartet object needed by LIBINT
    type(erisStack) :: buffer

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: i1234

    procedure(LibintInterface_buildLibInt), pointer :: pBuild !!procedure to calculate eris on LIBINT

    G12_ptr => G12Integrals_instance%libintG12

    nameOfSpecie = trim(MolecularSystem_getNameOfSpecie(specieID))

    call cpu_time(startTime)

    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID)

    !! Get number of shells and Cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize
    
    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))

    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))

    !!Real labels for contractions
    aux = 1
    do i = 1, numberOfContractions
       !!position for Cartesian contractions
       labelsOfContractions(i) = aux
       aux = aux + contractions(i)%numCartesianOrbital          
    end do

    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2

    sizeTotal = (totalNumberOfContractions *(totalNumberOfContractions + 1 ))/2
    sizeTotal = (sizeTotal *(sizeTotal + 1))/2

    !Get potential ID
    do i=1, InterPotential_instance%ssize
       if ( trim( MolecularSystem_instance%species(specieID)%symbol) == trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. trim( MolecularSystem_instance%species(specieID)%symbol) == trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie))) ) then
          potID=i
          exit
       end if
    end do

    !$OMP PARALLEL default(private), shared(contractions, maxAngularMoment, numberOfPrimitives, numberOfContractions), &
    !$OMP& shared(totalNumberOfContractions, labelsOfContractions, arraySize, sizeTotal, InterPotential_instance, nameOfSpecie, potID), &
    !$OMP& shared(G12Integrals_instance, CONTROL_instance, MolecularSystem_instance, specieID, auxCounter)

    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 400 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    !! open file for integrals
    if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(specieID)%isElectron) then
      open( UNIT=unitid,FILE=trim(fileid)//"E-ALPHA.ints", status='unknown', access='stream', form='unformatted')
    else
      open( UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//".ints", status='unknown', access='stream', form='unformatted')
    end if
    
    if( .not. G12Integrals_isInstanced()) then
      call G12Integrals_Constructor(maxAngularMoment, numberOfPrimitives)
    end if

    i1234 = -1
    G12_ptr => G12Integrals_instance%libIntG12

    if( .not. G12Integrals_isInstanced() ) then
       call G12Integrals_Constructor( maxAngularMoment, numberOfPrimitives)
    end if

    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    if(allocated(incompletGamma)) deallocate(incompletGamma)
    if(allocated(auxIntegrals)) deallocate(auxIntegrals)
    if(allocated(auxIntegralsB)) deallocate(auxIntegralsB)
    if(allocated(integralsValue)) deallocate(integralsValue)

    allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize), &
         auxIntegralsB(arraySize* arraySize* arraySize * arraySize), &
         integralsValue(arraySize* arraySize* arraySize* arraySize), &
         incompletGamma(0:MaxAngularMoment*4))

    allocate (buffer%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            buffer%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            buffer%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            buffer%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            buffer%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))


    counter = 0
    auxCounter = 0

    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       n = a
       do b = a, numberOfContractions
          u = b
          do r = n , numberOfContractions
             do s = u,  numberOfContractions
              
              i1234 = i1234 + 1
              if ( mod(i1234, nthreads) /= threadid) cycle

                !$OMP CRITICAL

                !! Calculates total angular moment
                sumAngularMoment =  contractions(a)%angularMoment + &
                     contractions(b)%angularMoment + &
                     contractions(r)%angularMoment + &
                     contractions(s)%angularMoment

                !! Calculates length of array for quartet (ab|rs)
                arraySize = contractions(a)%numCartesianOrbital * &
                     contractions(b)%numCartesianOrbital * &
                     contractions(r)%numCartesianOrbital * &
                     contractions(s)%numCartesianOrbital

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

                !! Distance AB, CD
                AB = contractions(aa)%origin - contractions(bb)%origin
                CD = contractions(rr)%origin - contractions(ss)%origin

                AB2 = dot_product(AB, AB)
                CD2 = dot_product(CD, CD)

                !! Setting up Libint structure
                G12Integrals_instance%libint%AB = AB
                G12Integrals_instance%libint%CD = CD

                !!start :)                                
                integralsValue(1:arraySize) = 0.0_8

                do l = 1, contractions(s)%length
                   do k = 1, contractions(r)%length
                      do j = 1, contractions(b)%length
                         do i = 1, contractions(a)%length

                            auxIntegrals = 0

                            do potSize=1, size(InterPotential_instance%Potentials(potID)%gaussianComponents)
                               contractionG12 => InterPotential_instance%Potentials(potID)%gaussianComponents(potSize)

                               do potLength = 1, contractionG12%length
                                  !!LIBINT2 PRIMQUARTET

                                  zeta = contractions(aa)%orbitalExponents(pi) + contractions(bb)%orbitalExponents(pj)
                                  eta = contractions(rr)%orbitalExponents(pk) + contractions(ss)%orbitalExponents(pl)
                                  rho  = (zeta * eta) / (zeta + eta) !Reduced exponent ABCD

                                  !! Exponents
                                  G12_ptr%zeta_A = contractions(aa)%orbitalExponents(pi)
                                  G12_ptr%zeta_B = contractions(bb)%orbitalExponents(pj)
                                  G12_ptr%zeta_C = contractions(rr)%orbitalExponents(pk)
                                  G12_ptr%zeta_D = contractions(ss)%orbitalExponents(pl)

                                  !! Squared exponents
                                  G12_ptr%zeta_A_2 = contractions(aa)%orbitalExponents(pi)**2
                                  G12_ptr%zeta_B_2 = contractions(bb)%orbitalExponents(pj)**2
                                  G12_ptr%zeta_C_2 = contractions(rr)%orbitalExponents(pk)**2
                                  G12_ptr%zeta_D_2 = contractions(ss)%orbitalExponents(pl)**2

                                  !! Appear in OS RR for ERIs

                                  !! One over 2.0*zeta
                                  G12_ptr%oo2z = (0.5_8 / zeta)
                                  !! One over 2.0*eta
                                  G12_ptr%oo2e = (0.5_8 / eta)
                                  !! One over 2.0*(zeta+eta)
                                  G12_ptr%oo2ze = (0.5_8 / (zeta + eta))
                                  !! rho over zeta
                                  G12_ptr%roz = (rho/zeta)
                                  !! rho over eta
                                  G12_ptr%roe = (rho/eta)

                                  !! Appear in standard OS RR for ERI and almost all other recurrence relations

                                  G12_ptr%AB_x = contractions(aa)%origin(1) - contractions(bb)%origin(1)
                                  G12_ptr%AB_y = contractions(aa)%origin(2) - contractions(bb)%origin(2)
                                  G12_ptr%AB_z = contractions(aa)%origin(3) - contractions(bb)%origin(3)

                                  G12_ptr%CD_x = contractions(rr)%origin(1) - contractions(ss)%origin(1)
                                  G12_ptr%CD_y = contractions(rr)%origin(2) - contractions(ss)%origin(2)
                                  G12_ptr%CD_z = contractions(rr)%origin(3) - contractions(ss)%origin(3)

                                  P  = ((contractions(aa)%orbitalExponents(pi) * contractions(aa)%origin) + &
                                       (contractions(bb)%orbitalExponents(pj) * contractions(bb)%origin)) / zeta
                                  Q  = ((contractions(rr)%orbitalExponents(pk) * contractions(rr)%origin) + &
                                       (contractions(ss)%orbitalExponents(pl) * contractions(ss)%origin)) / eta

                                  W  = ((zeta * P) + (eta * Q)) / (zeta + eta)

                                  G12_ptr%WP_x = W(1) - P(1)
                                  G12_ptr%WP_y = W(2) - P(2)
                                  G12_ptr%WP_z = W(3) - P(3)

                                  G12_ptr%WQ_x = W(1) - Q(1)
                                  G12_ptr%WQ_y = W(2) - Q(2)
                                  G12_ptr%WQ_z = W(3) - Q(3)

                                  G12_ptr%PA_x = P(1) - contractions(aa)%origin(1)
                                  G12_ptr%PA_y = P(2) - contractions(aa)%origin(2)
                                  G12_ptr%PA_z = P(3) - contractions(aa)%origin(3)

                                  G12_ptr%QC_x = Q(1) - contractions(rr)%origin(1)
                                  G12_ptr%QC_y = Q(2) - contractions(rr)%origin(2)
                                  G12_ptr%QC_z = Q(3) - contractions(rr)%origin(3)

                                  AB = contractions(aa)%origin - contractions(bb)%origin
                                  CD = contractions(rr)%origin - contractions(ss)%origin
                                  AC = contractions(aa)%origin - contractions(rr)%origin
                                  BD = contractions(bb)%origin - contractions(ss)%origin
                                  PQ = P - Q

                                  AB2 = dot_product(AB, AB)
                                  CD2 = dot_product(CD, CD)
                                  PQ2 = dot_product(PQ, PQ)

                                  !! Gamma incomplete argument
                                  incompletGammaArgument = rho*PQ2

                                  !!Overlap Factor
                                  s12 = ((Math_PI/zeta)**1.5_8) * exp(-((contractions(aa)%orbitalExponents(pi) * contractions(bb)%orbitalExponents(pj)) / zeta) * AB2)
                                  s34 = ((Math_PI/ eta)**1.5_8) * exp(-((contractions(rr)%orbitalExponents(pk) * contractions(ss)%orbitalExponents(pl)) /  eta) * CD2)
                                  s1234 = sqrt(rho/Math_PI) * s12 * s34

                                  call Math_fgamma0(sumAngularMoment,incompletGammaArgument,incompletGamma(0:sumAngularMoment))

                                  ssss = 0.0_8
                                  ssss(1:sumAngularMoment+1) = 2.0_8 * incompletGamma(0:sumAngularMoment) * s1234

                                  G12_ptr%LIBINT_T_SS_EREP_SS0 = ssss(1)
                                  G12_ptr%LIBINT_T_SS_EREP_SS1 = ssss(2)
                                  G12_ptr%LIBINT_T_SS_EREP_SS2 = ssss(3)
                                  G12_ptr%LIBINT_T_SS_EREP_SS3 = ssss(4)
                                  G12_ptr%LIBINT_T_SS_EREP_SS4 = ssss(5)
                                  G12_ptr%LIBINT_T_SS_EREP_SS5 = ssss(6)
                                  G12_ptr%LIBINT_T_SS_EREP_SS6 = ssss(7)
                                  G12_ptr%LIBINT_T_SS_EREP_SS7 = ssss(8)
                                  G12_ptr%LIBINT_T_SS_EREP_SS8 = ssss(9)
                                  G12_ptr%LIBINT_T_SS_EREP_SS9 = ssss(10)
                                  G12_ptr%LIBINT_T_SS_EREP_SS10 = ssss(11)
                                  G12_ptr%LIBINT_T_SS_EREP_SS11 = ssss(12)
                                  G12_ptr%LIBINT_T_SS_EREP_SS12 = ssss(13)
                                  G12_ptr%LIBINT_T_SS_EREP_SS13 = ssss(14)
                                  G12_ptr%LIBINT_T_SS_EREP_SS14 = ssss(15)
                                  G12_ptr%LIBINT_T_SS_EREP_SS15 = ssss(16)
                                  G12_ptr%LIBINT_T_SS_EREP_SS16 = ssss(17)
                                  G12_ptr%LIBINT_T_SS_EREP_SS17 = ssss(18)
                                  G12_ptr%LIBINT_T_SS_EREP_SS18 = ssss(19)
                                  G12_ptr%LIBINT_T_SS_EREP_SS19 = ssss(20)
                                  G12_ptr%LIBINT_T_SS_EREP_SS20 = ssss(21)

                                  !! Prefactor for recurrence relations from Weber and Daul, Comp. Phys. Comm. 158, 1 (2004).

                                  G12_ptr%LIBINT_T_SS_K0G12_SS_0 = &
                                       (((Math_PI**2) / (zeta * eta))**1.5_8) * ((rho/(rho+contractionG12%orbitalExponents(potLength)))**1.5_8) * &
                                       (exp(-((contractions(aa)%orbitalExponents(pi) * contractions(bb)%orbitalExponents(pj))/(zeta))*AB2)) * &
                                       (exp(-((contractions(rr)%orbitalExponents(pk) * contractions(ss)%orbitalExponents(pl))/(eta))*CD2)) * &
                                       (exp(-((rho * contractionG12%orbitalExponents(potLength))/(rho + contractionG12%orbitalExponents(potLength)))*PQ2))


                                  prefactor = 1.0_8/((zeta*eta) + (contractionG12%orbitalExponents(potLength)*(zeta+eta)))


                                  !! WD2004, Eq. 30, prefactor in front of (a0|k|c0) ---  -2 (ζb(ζq + γ )ABi + γ ζd CDi + γ ζq ACi)

                                  G12_ptr%R12kG12_pfac0_0_x = -2.0_8*((contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(ss)%orbitalExponents(pl)*CD(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(1)))

                                  G12_ptr%R12kG12_pfac0_0_y = -2.0_8*((contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(2) ) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(ss)%orbitalExponents(pl)*CD(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(2)))

                                  G12_ptr%R12kG12_pfac0_0_z = -2.0_8*((contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(ss)%orbitalExponents(pl)*CD(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(3)))

                                  G12_ptr%R12kG12_pfac0_1_x = -2.0_8*((contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(bb)%orbitalExponents(pj)*AB(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(1)))

                                  G12_ptr%R12kG12_pfac0_1_y = -2.0_8*((contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(bb)%orbitalExponents(pj)*AB(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(2)))

                                  G12_ptr%R12kG12_pfac0_1_z = -2.0_8*((contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(bb)%orbitalExponents(pj)*AB(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(3)))

                                  !! WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0)) -- +(ζq + γ )
                                  G12_ptr%R12kG12_pfac1_0 = (eta+contractionG12%orbitalExponents(potLength)) * (prefactor/2.0_8)
                                  G12_ptr%R12kG12_pfac1_1 = (zeta+contractionG12%orbitalExponents(potLength)) * (prefactor/2.0_8)

                                  !! WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0) -- gamma G12
                                  G12_ptr%R12kG12_pfac2 = contractionG12%orbitalExponents(potLength) * (prefactor/2.0_8)

                                  !! WD2004, Eq. 30, prefactor in front of curly brackets (excludes k) -- + k ζ
                                  G12_ptr%R12kG12_pfac3_0 = eta
                                  G12_ptr%R12kG12_pfac3_1 = zeta

                                  primitiveQuartet = G12Integrals_instance%libintG12

                                  if(arraySize == 1) then

                                     auxIntegralsB(1) = primitiveQuartet%LIBINT_T_SS_K0G12_SS_0

                                  else

                                     arraySsize(1) = arraySize

                                     allocate(temporalPtr(arraySize))
                                     temporalPtr = 0.0_8

                                     call LibintInterface_buildG12(contractions(ss)%angularMoment , &
                                          contractions(rr)%angularMoment , &
                                          contractions(bb)%angularMoment , &
                                          contractions(aa)%angularMoment , &
                                          maxAngularMoment, &
                                          arraySize, primitiveQuartet, temporalPtr)

                                     integralsPtr => temporalPtr

                                     auxIntegralsB(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy

                                  end if !!done by primitives

                                  auxIntegrals = auxIntegrals + auxIntegralsB &
                                       * InterPotential_instance%Potentials(potID)%gaussianComponents(potSize)%contractionCoefficients(potLength)


                               end do !! G12 length
                            end do!! done for potential


                            !!Normalize by primitive
                            m = 0
                            do ii = 1, contractions(aa)%numCartesianOrbital
                               do jj = 1, contractions(bb)%numCartesianOrbital
                                  do kk = 1, contractions(rr)%numCartesianOrbital
                                     do ll = 1, contractions(ss)%numCartesianOrbital
                                        m = m + 1
                                        auxIntegrals(m) = auxIntegrals(m) &
                                             * contractions(aa)%primNormalization(pi,ii) &
                                             * contractions(bb)%primNormalization(pj,jj) &
                                             * contractions(rr)%primNormalization(pk,kk) &
                                             * contractions(ss)%primNormalization(pl,ll)
                                     end do
                                  end do
                               end do
                            end do !! done by cartesian of contractions 

                            auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
                                 * contractions(aa)%contractionCoefficients(pi) &
                                 * contractions(bb)%contractionCoefficients(pj) &
                                 * contractions(rr)%contractionCoefficients(pk) &
                                 * contractions(ss)%contractionCoefficients(pl) !&

                            integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)

                         end do !l
                      end do !k
                   end do !j
                end do !i   !! done by contractions

                !!normalize by shell
                m = 0
                do ii = 1,  contractions(aa)%numCartesianOrbital
                   do jj = 1,  contractions(bb)%numCartesianOrbital
                      do kk = 1,  contractions(rr)%numCartesianOrbital
                         do ll = 1, contractions(ss)%numCartesianOrbital
                            m = m + 1
                            integralsValue(m) = integralsValue(m) &
                                 * contractions(aa)%contNormalization(ii) &
                                 * contractions(bb)%contNormalization(jj) &
                                 * contractions(rr)%contNormalization(kk) &
                                 * contractions(ss)%contNormalization(ll)
                         end do
                      end do
                   end do
                end do !! done by shell

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

                                  !! (they are repeated)

                               else 

                                  if(abs(integralsValue(m)) > 1.0D-10) then

                                     !print *, m,integralsValue(m)
                                     counter = counter + 1
                                     
                                     !$OMP ATOMIC
                                     auxCounter = auxCounter + 1

                                     buffer%a(counter) = pa
                                     buffer%b(counter) = pb
                                     buffer%c(counter) = pr
                                     buffer%d(counter) = ps
                                     buffer%integrals(counter) = integralsValue(m)

                                  end if

                                  if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

                                     write(unitid) &
                                          buffer%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                          buffer%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                          buffer%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                          buffer%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                          buffer%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
                                     counter = 0

                                  end if

                               end if
                            end if

                         end do
                      end do
                   end do
                end do

                !$OMP END CRITICAL

             end do
             u=r+1
          end do
       end do
    end do !! done by basis set

    buffer%a(counter+1) = -1
    buffer%b(counter+1) = -1
    buffer%c(counter+1) = -1
    buffer%d(counter+1) = -1
    buffer%integrals(counter+1) = 0.0_8


    write(unitid) &
         buffer%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         buffer%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         buffer%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         buffer%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         buffer%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    close(unitid)

    deallocate(buffer%a, buffer%b, buffer%c, buffer%d, buffer%integrals)

    !$OMP END PARALLEL

    call cpu_time(endTime)

    write(6,"(A,I10,A,F12.4,A)") "*****Time for ", auxCounter, "   "//trim(nameOfSpecie)//" non-zero G12 intra-species integrals ", (endTime - startTime), "(s)"

  end subroutine G12Integrals_DiskIntraSpecie

  !<
  !! calculate eris using libint library for all basis set (inter-specie)
  subroutine G12Integrals_G12diskInterSpecie(nameOfSpecie, otherNameOfSpecie,  specieID, otherSpecieID, isInterSpecies, isCouplingA, isCouplingB)

    implicit none

    integer,target :: specieID
    integer,target :: otherSpecieID
    character(*) :: nameOfSpecie
    character(*) :: otherNameOfSpecie
    logical, optional :: isInterSpecies
    logical, optional :: isCouplingA
    logical, optional :: isCouplingB

    logical :: interSpecies
    logical :: couplingA
    logical :: couplingB

    integer :: totalNumberOfContractions
    integer :: otherTotalNumberOfContractions
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: numberOfPrimitives
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: arraySize !! number of cartesian orbitals for maxAngularMoment
    integer :: u,m !! auxiliary iterators
    integer :: aa, bb, rr, ss !! permuted iterators (LIBINT)
    integer :: a, b, r, s !! not permuted iterators (original)
    integer :: pa, pb, pr, ps !! labels index
    integer :: ii, jj, kk, ll !! cartesian iterators for primitives and contractions
    integer :: aux, order !!auxiliary index
    integer :: arraySsize(1)
    integer :: sizeTotal
    integer :: counter, auxCounter

    integer :: potID, potLength, potSize

    integer,target :: i, j, k, l !! contraction length iterators
    integer,pointer :: pi, pj, pk, pl !! pointer to contraction length iterators
    integer,pointer :: poi, poj, pok, pol !! pointer to contraction length iterators
    integer, pointer :: pSpecieID, pOtherSpecieID !! pointer to species ID
    integer, allocatable :: labelsOfContractions(:) !! cartesian position of contractions in all basis set
    integer, allocatable :: otherLabelsOfContractions(:) !! cartesian position of contractions in all basis set
    !integer, allocatable :: buffer(:) !! avoid rep integrals
    !integer*1, allocatable :: buffer(:) !! avoid rep integrals

    real(8), dimension(:), pointer :: integralsPtr !! pointer to C array of integrals
    real(8), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: auxIntegrals(:) !!array with permuted integrals aux!
    real(8), allocatable :: auxIntegralsB(:) !!array with permuted integrals
    real(8), allocatable :: integralsValue (:) !!array with permuted integrals
    real(8), allocatable :: incompletGamma(:) !!array with incomplete gamma integrals

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3), AC(3), BD(3) !!geometric values that appear in Gaussian product
    real(8) :: zeta, eta, rho !! exponents... that appear in Gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !!geometric quantities that appear in Gaussian product
    real(8) :: incompletGammaArgument, ssss(21)
    real(8) :: prefactor

    real(8) :: startTime, endTime

    type(libint2), pointer :: G12_ptr
    type(ContractedGaussian), pointer :: contractionG12
    type(auxBasis) :: contractions(2)
    type(libint2), target :: primitiveQuartet !!Prim-quartet object needed by LIBINT
    type(erisStack) :: eris

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: i1234

    procedure(LibintInterface_buildLibInt), pointer :: pBuild !!procedure to calculate eris on LIBINT

    G12_ptr => G12Integrals_instance%libintG12

    call cpu_time(startTime)
    
    interSpecies = .true.
    couplingA = .false.
    couplingB = .false.

    if(present(isInterSpecies)) interSpecies = isInterSpecies
    if(present(isCouplingA)) couplingA = isCouplingA
    if(present(isCouplingB)) couplingB = isCouplingB

    !Get potential ID
    do i=1, InterPotential_instance%ssize
       if ( (trim(MolecularSystem_instance%species(specieID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. &
            trim(MolecularSystem_instance%species(otherSpecieID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie)) ) &
            ) .or. &
            (trim( MolecularSystem_instance%species(otherSpecieID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%specie))) .and. &
            trim( MolecularSystem_instance%species(specieID)%symbol) == &
            trim(String_getUpperCase(trim(InterPotential_instance%potentials(i)%otherSpecie)) ) &
            ) &
            ) then
          potID=i
          exit
       end if
    end do

    !$OMP PARALLEL default(private), shared(contractions, maxAngularMoment, numberOfPrimitives, numberOfContractions, otherNumberOfContractions, totalNumberOfContractions, otherTotalNumberOfContractions, labelsOfContractions, otherLabelsOfContractions, arraySize, sizeTotal, InterPotential_instance, nameOfSpecie, otherNameOfSpecie, potID, G12Integrals_instance, CONTROL_instance, MolecularSystem_instance, specieID, otherSpecieID, auxCounter)
    
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 400 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))
    
    !! open file for integrals
    open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecie)//"."//trim(otherNameOfSpecie)//".ints", &
         status='unknown', access='stream', form='unformatted')

    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions(1)%contractions)
    call MolecularSystem_getBasisSet(otherSpecieID, contractions(2)%contractions)

    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions(1)%contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

    !! Get number of shells and cartesian contractions (other specie)
    otherNumberOfContractions = size(contractions(2)%contractions)
    otherTotalNumberOfContractions = MolecularSystem_instance%species(otherSpecieID)%basisSetSize

    !! Libint constructor (just one time)
    maxAngularMoment = max(MolecularSystem_getMaxAngularMoment(specieID), MolecularSystem_getMaxAngularMoment(otherSpecieID))
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID) + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

    i1234 = -1
    G12_ptr => G12Integrals_instance%libIntG12

    if( .not. G12Integrals_isInstanced() ) then
       call G12Integrals_Constructor( maxAngularMoment, numberOfPrimitives)
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
    if(allocated(auxIntegralsB)) deallocate(auxIntegralsB)
    !if(allocated(buffer)) deallocate(buffer)

    allocate(integralsValue(arraySize* arraySize* arraySize* arraySize), &
         auxIntegrals(arraySize* arraySize* arraySize * arraySize), &
         auxIntegralsB(arraySize* arraySize* arraySize* arraySize), &
         incompletGamma(0:MaxAngularMoment*4))    

    counter = 0
    auxCounter = 0

    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
            eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))

    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       do b = a, numberOfContractions
          do r = 1 , otherNumberOfContractions
             do s = r,  otherNumberOfContractions

              i1234 = i1234 + 1
              if ( mod(i1234, nthreads) /= threadid) cycle

              !$OMP CRITICAL

              !! Total angular moment for this shell
                sumAngularMoment =  contractions(1)%contractions(a)%angularMoment + &
                     contractions(1)%contractions(b)%angularMoment + &
                     contractions(2)%contractions(r)%angularMoment + &
                     contractions(2)%contractions(s)%angularMoment

                !! size of array for quartet (ab|rs)
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
                !! Calculate iteratively primitives
                !!

                !!Distances AB, CD
                AB = contractions(pSpecieID)%contractions(aa)%origin - contractions(pSpecieID)%contractions(bb)%origin
                CD = contractions(pSpecieID)%contractions(rr)%origin - contractions(pSpecieID)%contractions(ss)%origin

                AB2 = dot_product(AB, AB)
                CD2 = dot_product(CD, CD)


                !! Setting up Libint structure
                G12Integrals_instance%libint%AB = AB
                G12Integrals_instance%libint%CD = CD

                !!start :)
                integralsValue(1:arraySize) = 0.0_8

                !! not-permuted loop
                
                do l = 1, contractions(otherSpecieID)%contractions(s)%length
                   do k = 1, contractions(otherSpecieID)%contractions(r)%length
                      do j = 1, contractions(specieID)%contractions(b)%length
                         do i = 1, contractions(specieID)%contractions(a)%length

                            auxIntegrals=0

                            do potSize=1, size(InterPotential_instance%Potentials(potID)%gaussianComponents)

                               !auxIntegralsB(1:arraySize) = 0.0_8
                               contractionG12 => InterPotential_instance%Potentials(potID)%gaussianComponents(potSize)

                               do potLength = 1, InterPotential_instance%Potentials(potID)%gaussianComponents(potSize)%length

                                  !! not-permuted loop

                                  !!LIBINT2 PRIMQUARTET

                                  zeta = contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi) + contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj) 
                                  eta = contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk) + contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)
                                  rho  = (zeta * eta) / (zeta + eta) !Reduced exponent ABCD

                                  !! Exponents
                                  G12_ptr%zeta_A = contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi)
                                  G12_ptr%zeta_B = contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)
                                  G12_ptr%zeta_C = contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk)
                                  G12_ptr%zeta_D = contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)

                                  !! Squared exponents
                                  G12_ptr%zeta_A_2 = contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi)**2
                                  G12_ptr%zeta_B_2 = contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)**2
                                  G12_ptr%zeta_C_2 = contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk)**2
                                  G12_ptr%zeta_D_2 = contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)**2

                                  !! Appear in OS RR for ERIs

                                  !! One over 2.0*zeta
                                  G12_ptr%oo2z = (0.5_8 / zeta)
                                  !! One over 2.0*eta
                                  G12_ptr%oo2e = (0.5_8 / eta)
                                  !! One over 2.0*(zeta+eta)
                                  G12_ptr%oo2ze = (0.5_8 / (zeta + eta))
                                  !! rho over zeta
                                  G12_ptr%roz = (rho/zeta)
                                  !! rho over eta
                                  G12_ptr%roe = (rho/eta)

                                  !! Appear in standard OS RR for ERI and almost all other recurrence relations

                                  G12_ptr%AB_x = contractions(pSpecieID)%contractions(aa)%origin(1) - contractions(pSpecieID)%contractions(bb)%origin(1)
                                  G12_ptr%AB_y = contractions(pSpecieID)%contractions(aa)%origin(2) - contractions(pSpecieID)%contractions(bb)%origin(2)
                                  G12_ptr%AB_z = contractions(pSpecieID)%contractions(aa)%origin(3) - contractions(pSpecieID)%contractions(bb)%origin(3)

                                  G12_ptr%CD_x = contractions(pOtherSpecieID)%contractions(rr)%origin(1) - contractions(pOtherSpecieID)%contractions(ss)%origin(1)
                                  G12_ptr%CD_y = contractions(pOtherSpecieID)%contractions(rr)%origin(2) - contractions(pOtherSpecieID)%contractions(ss)%origin(2)
                                  G12_ptr%CD_z = contractions(pOtherSpecieID)%contractions(rr)%origin(3) - contractions(pOtherSpecieID)%contractions(ss)%origin(3)

                                  P  = ((contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi) * contractions(pSpecieID)%contractions(aa)%origin) + &
                                       (contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj) * contractions(pSpecieID)%contractions(bb)%origin)) / zeta
                                  Q  = ((contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk) * contractions(pOtherSpecieID)%contractions(rr)%origin) + &
                                       (contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl) * contractions(pOtherSpecieID)%contractions(ss)%origin)) / eta

                                  W  = ((zeta * P) + (eta * Q)) / (zeta + eta)

                                  G12_ptr%WP_x = W(1) - P(1)
                                  G12_ptr%WP_y = W(2) - P(2)
                                  G12_ptr%WP_z = W(3) - P(3)

                                  G12_ptr%WQ_x = W(1) - Q(1)
                                  G12_ptr%WQ_y = W(2) - Q(2)
                                  G12_ptr%WQ_z = W(3) - Q(3)

                                  G12_ptr%PA_x = P(1) - contractions(pSpecieID)%contractions(aa)%origin(1)
                                  G12_ptr%PA_y = P(2) - contractions(pSpecieID)%contractions(aa)%origin(2)
                                  G12_ptr%PA_z = P(3) - contractions(pSpecieID)%contractions(aa)%origin(3)

                                  G12_ptr%QC_x = Q(1) - contractions(pOtherSpecieID)%contractions(rr)%origin(1)
                                  G12_ptr%QC_y = Q(2) - contractions(pOtherSpecieID)%contractions(rr)%origin(2)
                                  G12_ptr%QC_z = Q(3) - contractions(pOtherSpecieID)%contractions(rr)%origin(3)

                                  AB = contractions(pSpecieID)%contractions(aa)%origin - contractions(pSpecieID)%contractions(bb)%origin
                                  CD = contractions(pOtherSpecieID)%contractions(rr)%origin - contractions(pOtherSpecieID)%contractions(ss)%origin
                                  AC = contractions(pSpecieID)%contractions(aa)%origin - contractions(pOtherSpecieID)%contractions(rr)%origin
                                  BD = contractions(pSpecieID)%contractions(bb)%origin - contractions(pOtherSpecieID)%contractions(ss)%origin
                                  PQ = P - Q

                                  AB2 = dot_product(AB, AB)
                                  CD2 = dot_product(CD, CD)
                                  PQ2 = dot_product(PQ, PQ)

                                  !! Gamma incomplete argument
                                  incompletGammaArgument = rho*PQ2

                                  !!Overlap Factor
                                  s12 = ((Math_PI/zeta)**1.5_8) * exp(-((contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi) * contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)) / zeta) * AB2)
                                  s34 = ((Math_PI/ eta)**1.5_8) * exp(-((contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk) * contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)) /  eta) * CD2)
                                  s1234 = sqrt(rho/Math_PI) * s12 * s34

                                  call Math_fgamma0(sumAngularMoment,incompletGammaArgument,incompletGamma(0:sumAngularMoment))

                                  ssss = 0.0_8
                                  ssss(1:sumAngularMoment+1) = 2.0_8 * incompletGamma(0:sumAngularMoment) * s1234

                                  G12_ptr%LIBINT_T_SS_EREP_SS0 = ssss(1)
                                  G12_ptr%LIBINT_T_SS_EREP_SS1 = ssss(2)
                                  G12_ptr%LIBINT_T_SS_EREP_SS2 = ssss(3)
                                  G12_ptr%LIBINT_T_SS_EREP_SS3 = ssss(4)
                                  G12_ptr%LIBINT_T_SS_EREP_SS4 = ssss(5)
                                  G12_ptr%LIBINT_T_SS_EREP_SS5 = ssss(6)
                                  G12_ptr%LIBINT_T_SS_EREP_SS6 = ssss(7)
                                  G12_ptr%LIBINT_T_SS_EREP_SS7 = ssss(8)
                                  G12_ptr%LIBINT_T_SS_EREP_SS8 = ssss(9)
                                  G12_ptr%LIBINT_T_SS_EREP_SS9 = ssss(10)
                                  G12_ptr%LIBINT_T_SS_EREP_SS10 = ssss(11)
                                  G12_ptr%LIBINT_T_SS_EREP_SS11 = ssss(12)
                                  G12_ptr%LIBINT_T_SS_EREP_SS12 = ssss(13)
                                  G12_ptr%LIBINT_T_SS_EREP_SS13 = ssss(14)
                                  G12_ptr%LIBINT_T_SS_EREP_SS14 = ssss(15)
                                  G12_ptr%LIBINT_T_SS_EREP_SS15 = ssss(16)
                                  G12_ptr%LIBINT_T_SS_EREP_SS16 = ssss(17)
                                  G12_ptr%LIBINT_T_SS_EREP_SS17 = ssss(18)
                                  G12_ptr%LIBINT_T_SS_EREP_SS18 = ssss(19)
                                  G12_ptr%LIBINT_T_SS_EREP_SS19 = ssss(20)
                                  G12_ptr%LIBINT_T_SS_EREP_SS20 = ssss(21)

                                  !! Prefactor for recurrence relations from Weber and Daul, Comp. Phys. Comm. 158, 1 (2004).

                                  G12_ptr%LIBINT_T_SS_K0G12_SS_0 = &
                                       (((Math_PI**2) / (zeta * eta))**1.5_8) * ((rho/(rho+contractionG12%orbitalExponents(potLength)))**1.5_8) * &
                                       (exp(-((contractions(pSpecieID)%contractions(aa)%orbitalExponents(pi) * contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj))/(zeta))*AB2)) * &
                                       (exp(-((contractions(pOtherSpecieID)%contractions(rr)%orbitalExponents(pk) * contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl))/(eta))*CD2)) * &
                                       (exp(-((rho * contractionG12%orbitalExponents(potLength))/(rho + contractionG12%orbitalExponents(potLength)))*PQ2))


                                  prefactor = 1.0_8/((zeta*eta) + (contractionG12%orbitalExponents(potLength)*(zeta+eta)))


                                  !! WD2004, Eq. 30, prefactor in front of (a0|k|c0) ---  -2 (ζb(ζq + γ )ABi + γ ζd CDi + γ ζq ACi)

                                  G12_ptr%R12kG12_pfac0_0_x = -2.0_8*((contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*CD(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(1)))

                                  G12_ptr%R12kG12_pfac0_0_y = -2.0_8*((contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(2) ) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*CD(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(2)))

                                  G12_ptr%R12kG12_pfac0_0_z = -2.0_8*((contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*(zeta + contractionG12%orbitalExponents(potLength)) * AB(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*CD(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*zeta*AC(3)))

                                  G12_ptr%R12kG12_pfac0_1_x = -2.0_8*((contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*AB(1)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(1)))

                                  G12_ptr%R12kG12_pfac0_1_y = -2.0_8*((contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*AB(2)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(2)))

                                  G12_ptr%R12kG12_pfac0_1_z = -2.0_8*((contractions(pOtherSpecieID)%contractions(ss)%orbitalExponents(pl)*(eta + contractionG12%orbitalExponents(potLength)) * CD(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*contractions(pSpecieID)%contractions(bb)%orbitalExponents(pj)*AB(3)) &
                                       + (contractionG12%orbitalExponents(potLength)*eta*BD(3)))

                                  !! WD2004, Eq. 30, prefactor in front of (a-1 0|k|c0)) -- +(ζq + γ )
                                  G12_ptr%R12kG12_pfac1_0 = (eta+contractionG12%orbitalExponents(potLength)) * (prefactor/2.0_8)
                                  G12_ptr%R12kG12_pfac1_1 = (zeta+contractionG12%orbitalExponents(potLength)) * (prefactor/2.0_8)

                                  !! WD2004, Eq. 30, prefactor in front of (a0|k|c-1 0) -- gamma G12
                                  G12_ptr%R12kG12_pfac2 = contractionG12%orbitalExponents(potLength) * (prefactor/2.0_8)

                                  !! WD2004, Eq. 30, prefactor in front of curly brackets (excludes k) -- + k ζ
                                  G12_ptr%R12kG12_pfac3_0 = eta
                                  G12_ptr%R12kG12_pfac3_1 = zeta

                                  primitiveQuartet = G12Integrals_instance%libintG12

                                  if(arraySize == 1) then

                                     auxIntegralsB(1) = primitiveQuartet%LIBINT_T_SS_K0G12_SS_0

                                  else

                                     arraySsize(1) = arraySize

                                     allocate(temporalPtr(arraySize))
                                     temporalPtr = 0.0_8

                                     call LibintInterface_buildG12(&
                                          contractions(pOtherSpecieID)%contractions(ss)%angularMoment , &
                                          contractions(pOtherSpecieID)%contractions(rr)%angularMoment , &
                                          contractions(pSpecieID)%contractions(bb)%angularMoment , &
                                          contractions(pSpecieID)%contractions(aa)%angularMoment , &
                                          maxAngularMoment, &
                                          arraySize, primitiveQuartet, temporalPtr)

                                     integralsPtr => temporalPtr

                                     auxIntegralsB(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy

                                  end if !!done by primitives

                                  auxIntegrals = auxIntegrals + auxIntegralsB &
                                       * InterPotential_instance%Potentials(potID)%gaussianComponents(potSize)%contractionCoefficients(potLength)  
                                  
                               end do !! G12 length
                            end do!! done for potential

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
                end do !!done by contractions

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

                                  !$OMP ATOMIC

                                  auxCounter = auxCounter + 1

                                  eris%a(counter) = pa
                                  eris%b(counter) = pb
                                  eris%c(counter) = pr
                                  eris%d(counter) = ps
                                  eris%integrals(counter) = integralsValue(m)

                               end if


                               if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

                                  write(unitid) &
                                       eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                                  counter = 0

                               end if
                            end if !! Stack control

                         end do
                      end do
                   end do
                end do !! Done write to disk

                !$OMP END CRITICAL

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

    write(unitid) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE),&
         eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    close(unitid)

    deallocate(eris%a, eris%b, eris%c, eris%d, eris%integrals)
    !deallocate(buffer)

    !$OMP END PARALLEL

    call cpu_time(endTime)

    write(6,"(A,I10,A,F12.2,A)") "*****Time for ", auxCounter, "   "//trim(nameOfSpecie)//"/"//trim(otherNameOfSpecie)//" non-zero G12 Intra-species integrals ", (endTime - startTime), "(s)"

  end subroutine G12Integrals_G12diskInterSpecie

  !!>
  !! @brief Returns whether the object has been instanced or not.
  function G12Integrals_isInstanced( ) result( output )
    implicit  none

    logical :: output

    output = G12Integrals_instance%isInstanced

  end function G12Integrals_isInstanced

end module G12Integrals_
