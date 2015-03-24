!!******************************************************************************
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
!! @brief Libint 1.04 interface for Repulsion Derivative Module.
!!        This module contains all basic functions of the repulsion derivatives calculations
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2015-03-17
!!
!! <b> History: </b>
!!
!!   - <tt> 2015-03-17 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs,
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module RepulsionDerivatives_
  use Exception_
  use MolecularSystem_
  use RepulsionDerivativesTypes_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use, intrinsic :: iso_c_binding
  implicit none
  
! #define contr(n,m) contractions(n)%contractions(m)
  
  !> @brief Type for libint/c++ library
  type, public :: RepulsionDerivatives
     logical :: isInstanced
     integer :: maxAngularMoment
     integer :: numberOfPrimitives
     integer :: numberOfCartesians
     integer :: libintStorage
     integer :: libderivStorage
     type(lib_int) :: libint
     type(lib_deriv) :: libderiv
  end type RepulsionDerivatives

  !>
  !! pointers to array functions on  libint.a, libderiv.a y libr12.a
  type(c_funptr), dimension(0:3,0:3,0:3,0:3), bind(c) :: build_deriv1_eri
  
  interface

     !>
     !!Interfaz a libint.a
     subroutine RepulsionDerivatives_initLibIntBase() bind(C, name="init_libint_base")
       implicit none
       
     end subroutine RepulsionDerivatives_initLibIntBase

     function RepulsionDerivatives_initLibInt(libInt, maxAngMoment, numberOfPrimitives) bind(C, name="init_libint")
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none

       integer(kind=c_int) :: RepulsionDerivatives_initLibInt
       type(lib_int) :: libInt
       integer(kind=c_int), value :: maxAngMoment
       integer(kind=c_int), value :: numberOfPrimitives

     end function RepulsionDerivatives_initLibInt

     subroutine RepulsionDerivatives_freeLibInt(libInt) bind(C, name="free_libint")
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none
       
       type(lib_int) :: libInt
       
     end subroutine RepulsionDerivatives_freeLibInt
       
     !>
     !!Interfaz a libderiv.a
     subroutine RepulsionDerivatives_initLibDerivBase() bind(C, name="init_libderiv_base")
       implicit none
       
     end subroutine RepulsionDerivatives_initLibDerivBase

     function RepulsionDerivatives_libDerivStorage(&
          maxAngMoment, &
          numberOfPrimitives, &
          numberOfCartesians) bind(C, name="libderiv1_storage_required")
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none

       integer(kind=c_int) :: RepulsionDerivatives_libDerivStorage
       integer(kind=c_int), value :: maxAngMoment
       integer(kind=c_int), value :: numberOfPrimitives
       integer(kind=c_int), value :: numberOfCartesians
     end function RepulsionDerivatives_libDerivStorage
     
     function RepulsionDerivatives_initLibDeriv1(&
          libderiv, &
          maxAngMoment, &
          numberOfPrimitives, &
          numberOfCartesians) bind(C, name="init_libderiv1")
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none

       integer(kind=c_int) :: RepulsionDerivatives_initLibDeriv1
       type(lib_deriv) :: libderiv
       integer(kind=c_int), value :: maxAngMoment
       integer(kind=c_int), value :: numberOfPrimitives
       integer(kind=c_int), value :: numberOfCartesians
     end function RepulsionDerivatives_initLibDeriv1
     
     subroutine RepulsionDerivatives_buildLibDeriv(libderiv, numberOfPrimitives) bind(C)
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none
       
       ! type(c_ptr) :: RepulsionDerivatives_buildLibDeriv
       type(lib_deriv) :: libderiv
       integer(kind=c_int), value :: numberOfPrimitives
       
     end subroutine RepulsionDerivatives_buildLibDeriv
     
     subroutine RepulsionDerivatives_freeLibDeriv(libderiv) bind(C, name="free_libderiv")
       use RepulsionDerivativesTypes_
       use, intrinsic :: iso_c_binding
       implicit none
       
       type(lib_deriv) :: libderiv
       
     end subroutine RepulsionDerivatives_freeLibDeriv
     
  end interface
  
  !> @brief Singleton lock
  type(RepulsionDerivatives), public :: RepulsionDerivatives_instance

!   !> @brief Integrals Stack
!   type(erisStack), private :: eris

contains
  
  !>
  ! @brief Constructor por omision, se llama una sola vez!
  ! @author J.M. Rodas 2015
  ! @version 1.0
  subroutine RepulsionDerivatives_constructor(maxAngMoment, numberOfPrimitives, numberOfCartesians)
    implicit none
    integer, intent(in) :: maxAngMoment
    integer, intent(in) :: numberOfPrimitives
    integer, intent(in) :: numberOfCartesians

    integer :: storage
    integer(kind=c_int) :: libIntStorage_c
    integer(kind=c_int) :: libDerivStorage_c
    integer(kind=c_int) :: maxAngMoment_c
    integer(kind=c_int) :: numberOfPrimitives_c
    integer(kind=c_int) :: numberOfCartesians_c

    if (.not. RepulsionDerivatives_isInstanced()) then
       
       RepulsionDerivatives_instance%maxAngularMoment = maxAngMoment
       RepulsionDerivatives_instance%numberOfPrimitives = numberOfPrimitives
       RepulsionDerivatives_instance%numberOfCartesians = numberOfCartesians

       call RepulsionDerivatives_initLibIntBase()
       call RepulsionDerivatives_initLibDerivBase()

       maxAngMoment_c = RepulsionDerivatives_instance%maxAngularMoment
       numberOfPrimitives_c = RepulsionDerivatives_instance%numberOfPrimitives
       numberOfCartesians_c = maxAngMoment_c*maxAngMoment_c*maxAngMoment_c*maxAngMoment_c

       libDerivStorage_c = RepulsionDerivatives_initLibDeriv1(&
            RepulsionDerivatives_instance%libderiv, &
            maxAngMoment_c, &
            numberOfPrimitives_c, &
            numberOfCartesians_c)

       print*, "LibDerivStorage: ", libDerivStorage_c

       RepulsionDerivatives_instance%libintStorage = libIntStorage_c
       RepulsionDerivatives_instance%libderivStorage = libDerivStorage_c
          
       RepulsionDerivatives_instance%isInstanced = .true.
       
    else
       
       call RepulsionDerivatives_exception( ERROR, "in libint interface constructor function",&
            "you must destroy this object before to use it again")
       
    end if

  end subroutine RepulsionDerivatives_constructor

  !>
  !! @brief Destructor por omision para ser llamado una sola vez
  !! @author J.M. Rodas
  !! @version 1.0
  subroutine RepulsionDerivatives_destructor()
    implicit none
    
    if (RepulsionDerivatives_isInstanced()) then
       ! call RepulsionDerivatives_freeLibInt(RepulsionDerivatives_instance%libint)
       ! call RepulsionDerivatives_freeLibDeriv(RepulsionDerivatives_instance%libderiv) !! it does not work
          
       RepulsionDerivatives_instance%isInstanced = .false.
       
    else
       
       call RepulsionDerivatives_exception( ERROR, "in libint interface destructor function", "you must instantiate this object before to destroying it")
       
    end if
    
  end subroutine RepulsionDerivatives_destructor

  ! !>
  ! !! @brief Muestra informacion del objeto (una sola vez)
  ! !! @author E. F. Posada, 2010
  ! !! @version 2.0
  ! subroutine RepulsionDerivatives_show()
  !   implicit none

  !   write(*,  "(A)")  " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
  !   write(*, "(A)")   " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2011   "
  !   write(*, "(A)")   " ----------------------------------------------------------------------"
  !   write(*, "(A)" )  " LIBINT parameters: "
  !   !write(*, "(A, T44,A)")  " work ", trim(RepulsionDerivatives_instance%job)
  !   !write(*, "(A, T44,A)")  " Storage ", "DISK"
  !   write(*, "(A, T43,I5)") " Stack size", CONTROL_instance%INTEGRAL_STACK_SIZE
  !   !write(*, "(A, T43,I5)") " maxAngularMoment", RepulsionDerivatives_instance%maxAngularMoment
  !   write(*, "(A, T38,I10)") " number of primitives", RepulsionDerivatives_instance%numberOfPrimitives
  !   write(*, "(A, T38,I10/)")" Memory required (in words)", RepulsionDerivatives_instance%libintStorage

  ! end subroutine RepulsionDerivatives_show


  !>
  !! @brief calculate eris using libint library for all basis set (intra-specie)
  !! @author E. F. Posada, 2010
  !! @version 2.0
  !! @info Tested
  subroutine RepulsionDerivatives_getDerive(this, a, b, r, s, deriveValue, specieID)
    implicit none
    type(ContractedGaussian), intent(in):: this(:)
    integer, intent(in) :: a, b, r, s, specieID
    real(8), allocatable :: deriveValue(:)
    
    integer :: maxAngularMoment !x
    integer :: sumAngularMoment !x
    integer :: numberOfPrimitives !x
    integer :: maxNumberOfPrimitives !x
    integer :: maxNumberOfCartesians !x
    integer :: maxNPrimSize, maxNCartSize !x
    logical :: p12, p34, p13p24 !x
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT) !x
    integer :: am1, am2, am3, am4
    integer :: pa, pb, pr, ps !< labels index !x
    real(8) :: auxExponentA, auxCoefficientA, auxContConstantA, auxPrimConstantA, c1
    real(8) :: auxExponentB, auxCoefficientB, auxContConstantB, auxPrimConstantB, c2
    real(8) :: auxExponentR, auxCoefficientR, auxContConstantR, auxPrimConstantR, c3
    real(8) :: auxExponentS, auxCoefficientS, auxContConstantS, auxPrimConstantS, c4
    ! integer :: apa, apb, apr, aps !< labels index
    ! integer :: ii, jj, kk, ll !< cartesian iterators for primitives and contractions
    integer :: aux, order !<auxiliary index !x
    integer :: arraySsize(1)
    integer :: counter, auxCounter

    integer(8) :: control

    integer,target :: i, j, k, l !< contraction length iterators
    integer,pointer :: pi, pj, pk, pl !< pointer to contraction length iterators !x
    integer,pointer :: poi, poj, pok, pol !< pointer to contraction length iterators !x

    real(8), dimension(:), pointer :: integralsPtr !< pointer to C array of integrals
    real(c_double), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: auxIntegrals(:) !<array with permuted integrals aux!
    real(8), allocatable :: integralsValue (:) !<array with permuted integrals
    real(8), allocatable :: incompletGamma(:) !<array with incomplete gamma integrals !x

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !<geometric values that appear in gaussian product !x
    real(8) :: zeta, ooz, oo2z, nu, oon, oo2n, oo2zn, rho !< exponents... that appear in gaussian product !x
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !<geometric quantities that appear in gaussian product !x
    real(8) :: incompletGammaArgument !x   
    
    character(50) :: fileNumber
    integer(8) :: ssize

    type(prim_data), target :: primitiveQuartet !<Primquartet object needed by LIBINT
    type(c_ptr) :: resultPc !< array of integrals from C (LIBINT)
    procedure(RepulsionDerivatives_buildLibDeriv), pointer :: pBuild !<procedure to calculate eris on LIBINT !x
    
    integer :: contractionNumberdebug, primitiveCounterdebug, contractionNumberperOrbital, totalIntegralswithP

    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    maxNumberOfPrimitives = MolecularSystem_getMaxNumberofPrimitives(specieID)
    maxNumberOfCartesians = MolecularSystem_getMaxNumberofCartesians(specieID)
    maxNPrimSize = maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives
    maxNCartSize = maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians


    arraySsize(1) = this(a)%numCartesianOrbital * &
         this(b)%numCartesianOrbital * &
         this(r)%numCartesianOrbital * &
         this(s)%numCartesianOrbital

    sumAngularMoment = this(a)%angularMoment + this(b)%angularMoment + this(r)%angularMoment + this(s)%angularMoment

    ! if(allocated(incompletGamma)) deallocate(incompletGamma)
    ! allocate(incompletGamma(0:maxAngularMoment+1))
    if(allocated(incompletGamma)) deallocate(incompletGamma)
    allocate(incompletGamma(0:sumAngularMoment+1))

    ! if(allocated(primitiveQuartet)) deallocate(primitiveQuartet)
    ! allocate(primitiveQuartet(maxNumberOfPrimitives))

    ! Libderiv constructor (solo una vez)
    if( RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_destructor()
    end if
    
    if( .not. RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_constructor(sumAngularMoment,maxNPrimSize,maxNCartSize)
       !! DEBUG
       !call RepulsionDerivatives_show()
    end if

    if(allocated(deriveValue)) deallocate(deriveValue)
    allocate(deriveValue(0:6*9))

    deriveValue = 0.0_8

    order = 0

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

    p12 = .false.
    p34 = .false.
    p13p24 = .false.
    
    if(this(a)%angularMoment .lt. this(b)%angularMoment) then
       aa = b
       bb = a

       pi => j
       pj => i

       poi => j
       poj => i

       order = order + 1
       p12 = .true.
    end if

    if(this(r)%angularMoment .lt. this(s)%angularMoment) then
       rr = s
       ss = r

       pk => l
       pl => k

       pok => l
       pol => k

       order = order + 3
       p34 = .true.
    end if

    if((this(a)%angularMoment + this(b)%angularMoment) .gt. (this(r)%angularMoment + this(s)%angularMoment)) then
       aux = aa
       aa = rr
       rr = aux

       aux = bb
       bb = ss
       ss = aux
       p13p24 = .true.
       
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

    write(*,"(A)") "Contraida:"
    ! write(*,"(3(A,I))") "maxam: ",  maxAngularMoment, " maxnprim: ", maxNPrimSize, " maxncart: ", maxNCartSize
    !write(*,"(A,I,A,I,A,I,A,I,A)") "(",aa,",",bb,"|",rr,",",ss,")"
    !! Si los offsets son necesarios hay que ponerlos aqui
    am1 = this(aa)%angularMoment
    am2 = this(bb)%angularMoment
    am3 = this(rr)%angularMoment
    am4 = this(ss)%angularMoment

   
    !!Distancias AB, CD
    AB = this(aa)%origin - this(bb)%origin
    CD = this(rr)%origin - this(ss)%origin

    AB2 = dot_product(AB, AB)
    CD2 = dot_product(CD, CD)    

    ! Asigna valores a la estrucutra Libderiv
    
    RepulsionDerivatives_instance%libderiv%AB = AB
    RepulsionDerivatives_instance%libderiv%CD = CD
    ! write(*,"(A)") "----------------------------------------"
    ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
    ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A6,I1,A1,I1,A1,I1,A1,I1,A1)") "(",am1,",",am2,"|",am3,",",am4,") -> (",am4,",",am3,"|",am2,",",am1,")"
    ! write(*,"(A)") "-----------------------------------------"
    ! write(*,"(A,3(F17.12))") "AB: ", AB(:)
    ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%AB(:)
    ! write(*,"(A)") "-----------------------------------------"

    numberOfPrimitives = 1

    do pa=1, this(aa)%length
       do pb=1, this(bb)%length
          do pr=1, this(rr)%length
             do ps=1, this(ss)%length
                auxExponentS = this(ss)%orbitalExponents(ps)
                auxCoefficientS = this(ss)%contractionCoefficients(ps)
                auxPrimConstantS = this(ss)%primNormalization(ps,1)
                auxContConstantS = this(ss)%contNormalization(1)
                c4 = auxCoefficientS*auxPrimConstantS*auxContConstantS

                auxExponentR = this(rr)%orbitalExponents(pr)
                auxCoefficientR = this(rr)%contractionCoefficients(pr)
                auxPrimConstantR = this(rr)%primNormalization(pr,1)
                auxContConstantR = this(rr)%contNormalization(1)
                c3 = auxCoefficientR*auxPrimConstantR*auxContConstantR

                auxExponentB = this(bb)%orbitalExponents(pb)
                auxCoefficientB = this(bb)%contractionCoefficients(pb)
                auxPrimConstantB = this(bb)%primNormalization(pb,1)
                auxContConstantB = this(bb)%contNormalization(1)
                c2 = auxCoefficientB*auxPrimConstantB*auxContConstantB


                auxExponentA = this(aa)%orbitalExponents(pa)
                auxCoefficientA = this(aa)%contractionCoefficients(pa)
                auxPrimConstantA = this(aa)%primNormalization(pa,1)
                auxContConstantA = this(aa)%contNormalization(1)
                c1 = auxCoefficientA*auxPrimConstantA*auxContConstantA

                zeta = auxExponentA + auxExponentB
                ooz = 1.0_8/zeta
                oo2z = 1.0_8/(2.0_8*zeta)

                P = (auxExponentA*this(aa)%origin + auxExponentB*this(bb)%origin)*ooz

                primitiveQuartet%U(1:3,1)= P - this(aa)%origin
                primitiveQuartet%U(1:3,2)= P - this(bb)%origin

                s12 = ((Math_PI*ooz)**1.5_8) * exp(-auxExponentA*auxExponentB*ooz*AB2)!*c1*c2

                nu = auxExponentR + auxExponentS
                oon = 1.0_8/nu
                oo2n = 1.0_8/(2.0_8*nu)
                oo2zn = 1.0_8/(2.0_8*(zeta+nu))
                rho = (zeta*nu)/(zeta+nu)

                Q = (auxExponentR*this(rr)%origin + auxExponentS*this(ss)%origin)*oon
                primitiveQuartet%U(1:3,3)= Q - this(rr)%origin
                primitiveQuartet%U(1:3,4)= Q - this(ss)%origin

                PQ = P - Q
                PQ2 = dot_product(PQ, PQ)

                W  = ((zeta * P) + (nu * Q)) / (zeta + nu)

                primitiveQuartet%U(1:3,5)= (W - P)
                primitiveQuartet%U(1:3,6)= (W - Q)
                               
                primitiveQuartet%oo2z = oo2z
                primitiveQuartet%oo2n = oo2n
                primitiveQuartet%oo2zn = oo2zn
                primitiveQuartet%poz = rho * ooz
                primitiveQuartet%pon = rho * oon
                primitiveQuartet%twozeta_a = 2.0 * auxExponentA
                primitiveQuartet%twozeta_b = 2.0 * auxExponentB
                primitiveQuartet%twozeta_c = 2.0 * auxExponentR
                primitiveQuartet%twozeta_d = 2.0 * auxExponentS

                incompletGammaArgument = rho*PQ2
                ! call Math_fgamma0(maxAngularMoment+1,incompletGammaArgument,incompletGamma(0:maxAngularMoment+1))
                call Math_fgamma0(sumAngularMoment+1,incompletGammaArgument,incompletGamma(0:sumAngularMoment+1))

                s34 = ((Math_PI*oon)**1.5_8) * exp(-auxExponentR*auxExponentS*oon*CD2)!*c3*c4

                s1234 = 2.0_8*sqrt(rho/Math_PI) * s12 * s34

                do i=1, sumAngularMoment+1
                   primitiveQuartet%F(i) = incompletGamma(i-1)*s1234
                end do

                ! do i=1, maxAngularMoment+1
                !    primitiveQuartet%F(i) = incompletGamma(i-1)*s1234
                ! end do
                ! arraySsize(1) = arraySize

                RepulsionDerivatives_instance%libderiv%PrimQuartet = c_loc(primitiveQuartet) !ok

                write(*,"(A)") "-----------------------------------------"
                write(*,"(A,I8)") "Array Size: ", arraySsize(1)
                write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",pa,",",pb,"|",pr,",",ps,")"
                write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",am1,",",am2,"|",am3,",",am4,")"
                write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2z
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2n
                ! write(*,"(A)") "-----------------------------------------"
                ! RepulsionDerivatives_instance%libderiv%PrimQuartet(numberOfPrimitives) = c_loc(primitiveQuartet)
                call c_f_procpointer(build_deriv1_eri(am4,am3,am2,am1), pBuild)

                write(*,"(A)") " Pase el llamado al c_f_procpointer"

                call pBuild(RepulsionDerivatives_instance%libderiv,1)

                write(*,"(A)") " Pase el llamado al pBuild"

                write(*,"(A)") "-----------------------------------------"
                do k=1,12
                   if(k==4 .or. k==5 .or. k==6) cycle
                    resultPc = RepulsionDerivatives_instance%libderiv%ABCD(k)
                    call c_f_pointer(resultPc, temporalPtr, arraySsize)
                   
                   do i=1,arraySsize(1)
                      write(*,"(2I, f17.12)") k, i, temporalPtr(i)
                      ! work_forces(i,k) = tmp_data(i)
                   end do
                end do
                write(*,"(A)") "-----------------------------------------"
                !! Get numbers from memory!
                ! arraySsize(1) = 3_8

                
                ! integralsPtr => temporalPtr
                
                ! call c_f_pointer(resultPc, temporalPtr, 1)

                ! integralsPtr => temporalPtr
                
                numberOfPrimitives = numberOfPrimitives + 1
             end do
          end do
       end do
    end do

    ! call c_f_procpointer(build_deriv1_eri(am1,am2,am3,am4), pBuild)

    ! write(*,"(A)") " Pase el llamado al c_f_procpointer"

    ! call pBuild(RepulsionDerivatives_instance%libderiv,numberOfPrimitives-1)

    ! write(*,"(A)") " Pase el llamado al pBuild"
    ! write(*,"(A)") "-----------------------------------------"
    ! do k=1,12
    !    if(k==4 .or. k==5 .or. k==6) cycle
    !    resultPc = RepulsionDerivatives_instance%libderiv%ABCD(k)
    !    ! write(*,"(f17.12)") RepulsionDerivatives_instance%libderiv%ABCD(k)
    !    call c_f_pointer(resultPc, temporalPtr, arraySsize)
    !    do i=1,arraySsize(1)
    !       integralsPtr => temporalPtr
    !       write(*,"(f17.12)") integralsPtr(i)
    !       ! work_forces(i,k) = tmp_data(i)
    !    end do
    ! end do

    ! call c_f_procpointer(build_deriv1_eri(am1,am2,am3,am4), pBuild)
    
    ! call pBuild(RepulsionDerivatives_instance%libderiv,(numberOfPrimitives-1))

    ! build_deriv1_eri[am1][am2][am3][am4](&libderiv_, nprim);
    ! call build_deriv1_eri(am1,am2,am3,am4)
    ! do i=1, 166
    !    write(*,"(f17.12)") RepulsionDerivatives_instance%libderiv%ABCD(i)
    ! end do    
    ! write(*,"(A)") "Antes de llamar a libint"
!(RepulsionDerivatives_instance%libderiv, numberOfPrimitives-1))

    ! write(*,"(A)") "-----------------------------------------"
    ! do i=1, 166
    !    write(*,"(f17.12)") RepulsionDerivatives_instance%libderiv%ABCD(i)
    ! end do
    write(*,"(A)") "-----------------------------------------"
    ! write(*,"(A)") "Despues de llamar a libint"

  end subroutine RepulsionDerivatives_getDerive
  
  !>
  !! @brief calculate eris using libint library for all basis set (inter-species)
  !! @author E. F. Posada, 2010
  !! @version 2.0
  ! subroutine RepulsionDerivatives_computeInterSpecies(specieID, otherSpecieID, job, isInterSpecies, isCouplingA, isCouplingB)
  !   implicit none

  !   integer,target :: specieID
  !   integer,target :: otherSpecieID
  !   character(*), intent(in) :: job
  !   logical, optional :: isInterSpecies
  !   logical, optional :: isCouplingA
  !   logical, optional :: isCouplingB
    
  !   logical :: interSpecies
  !   logical :: couplingA
  !   logical :: couplingB
    
  !   integer :: auxSpecieID
  !   integer :: auxOtherSpecieID    
  !   integer :: totalNumberOfContractions
  !   integer :: otherTotalNumberOfContractions
  !   integer :: numberOfContractions
  !   integer :: otherNumberOfContractions
  !   integer :: numberOfPrimitives
  !   integer :: maxAngularMoment
  !   integer :: sumAngularMoment
  !   integer :: arraySize !! number of cartesian orbitals for maxAngularMoment
  !   integer :: n,u,m !! auxiliary itetators
  !   integer :: aa, bb, rr, ss !! permuted iterators (LIBINT)
  !   integer :: a, b, r, s !! not permuted iterators (original)
  !   integer*2 :: pa, pb, pr, ps !! labels index
  !   integer :: ii, jj, kk, ll !! cartesian iterators for contractions and contractions
  !   integer :: aux, order !!auxiliary index
  !   integer :: arraySsize(1)
  !   integer :: sizeTotal
  !   integer :: auxIndex
  !   integer :: counter, auxCounter
    
  !   integer,target :: i, j, k, l !! contraction length iterators
  !   integer,pointer :: pi, pj, pk, pl !! pointer to contraction length iterators
  !   integer,pointer :: poi, poj, pok, pol !! pointer to contraction length iterators
  !   integer, pointer :: pSpecieID, pOtherSpecieID !! pointer to species ID
  !   integer, allocatable :: labelsOfContractions(:) !! cartesian position of contractions in all basis set
  !   integer, allocatable :: otherLabelsOfContractions(:) !! cartesian position of contractions in all basis set

  !   real(8), dimension(:), pointer :: integralsPtr !! pointer to C array of integrals
  !   real(8), dimension(:), pointer :: temporalPtr

  !   real(8), allocatable :: auxIntegrals(:) !!array with permuted integrals aux!
  !   real(8), allocatable :: integralsValue (:) !!array with permuted integrals
  !   real(8), allocatable :: incompletGamma(:) !!array with incomplete gamma integrals
    
  !   real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !!geometric values that appear in gaussian product
  !   real(8) :: zeta, eta, rho !! exponents... that appear in gaussian product
  !   real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !!geometric quantities that appear in gaussian product
  !   real(8) :: incompletGammaArgument
    
  !   type(auxBasis) :: contractions(2)
  !   type(prim_data), target :: primitiveQuartet !!Primquartet object needed by LIBINT
  !   type(c_ptr) :: resultPc !! array of integrals from C (LIBINT)
  !   procedure(RepulsionDerivatives_buildLibInt), pointer :: pBuild !!procedure to calculate eris on LIBINT
    
  !   interSpecies = .true.
  !   couplingA = .false.
  !   couplingB = .false.
    
  !   if(present(isInterSpecies)) interSpecies = isInterSpecies
  !   if(present(isCouplingA)) couplingA = isCouplingA
  !   if(present(isCouplingB)) couplingB = isCouplingB
    
  !   !! open file for integrals
  !   open(UNIT=34,FILE=trim(MolecularSystem_instance%species(specieID)%name)//"."//trim(MolecularSystem_instance%species(otherSpecieID)%name)//".ints", &
  !        STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='Unformatted')
    
  !   !! Get basisSet
  !   call MolecularSystem_getBasisSet(specieID, contractions(1)%contractions)
  !   call MolecularSystem_getBasisSet(otherSpecieID, contractions(2)%contractions)
    
  !   !! Get number of shells and cartesian contractions
  !   numberOfContractions = size(contractions(1)%contractions)
  !   totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize

  !   !! Get number of shells and cartesian contractions (other specie)
  !   otherNumberOfContractions = size(contractions(2)%contractions)
  !   otherTotalNumberOfContractions = MolecularSystem_instance%species(otherSpecieID)%basisSetSize
    
  !   !! Libint constructor (solo una vez)
  !   maxAngularMoment = max(MolecularSystem_getMaxAngularMoment(specieID), MolecularSystem_getMaxAngularMoment(otherSpecieID))
  !   numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID) + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
    
  !   if( .not. RepulsionDerivatives_isInstanced() ) then
  !      call RepulsionDerivatives_constructor( maxAngularMoment, numberOfPrimitives, trim(job))
  !      call RepulsionDerivatives_show()
  !   end if
    
  !   !! allocating space for integrals just one time (do not put it inside do loop!!!)
  !   arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2
    
  !   sizeTotal = ((totalNumberOfContractions *(totalNumberOfContractions + 1 ))/2) * ((otherTotalNumberOfContractions *(otherTotalNumberOfContractions + 1 ))/2)
    
  !   !! Get contractions labels for integrals index
  !   if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
  !   allocate(labelsOfContractions(numberOfContractions))
    
  !   !!Real labels for contractions
  !   aux = 1
  !   do i = 1, numberOfContractions
  !      labelsOfContractions(i) = aux
  !      aux = aux + contractions(1)%contractions(i)%numCartesianOrbital
  !   end do

  !   !! Get contractions labels for integrals index (other specie)
  !   if (allocated(otherLabelsOfContractions)) deallocate(otherLabelsOfContractions)
  !   allocate(otherLabelsOfContractions(otherNumberOfContractions))

  !   !!Real labels for contractions (other specie)
  !   aux = 1
  !   do i = 1, otherNumberOfContractions
  !      otherLabelsOfContractions(i) = aux
  !      aux = aux + contractions(2)%contractions(i)%numCartesianOrbital
  !   end do
    
  !   !! Allocating some space
  !   if(allocated(incompletGamma)) deallocate(incompletGamma)
  !   if(allocated(auxIntegrals)) deallocate(auxIntegrals)
  !   if(allocated(integralsValue)) deallocate(integralsValue)

  !   allocate(incompletGamma(0:MaxAngularMoment*4))    
  !   allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize))
  !   allocate(integralsValue(arraySize* arraySize* arraySize* arraySize))

    
  !   counter = 0
  !   auxCounter = 0
    
  !   auxSpecieID = specieID
  !   auxOtherSpecieID = otherSpecieID

  !   specieID = 1
  !   otherSpecieID = 2
    
  !   !!Start Calculating integrals for each shell
  !   do a = 1, numberOfContractions
  !      do b = a, numberOfContractions
  !         do r = 1 , otherNumberOfContractions
  !            do s = r,  otherNumberOfContractions
                
  !               !!Calcula el momento angular total
  !               sumAngularMoment =  contractions(1)%contractions(a)%angularMoment + &
  !                    contractions(1)%contractions(b)%angularMoment + &
  !                    contractions(2)%contractions(r)%angularMoment + &
  !                    contractions(2)%contractions(s)%angularMoment
                
  !               !! Calcula el tamano del arreglo de integrales para la capa (ab|rs)
  !               arraySize = contractions(1)%contractions(a)%numCartesianOrbital * &
  !                    contractions(1)%contractions(b)%numCartesianOrbital * &
  !                    contractions(2)%contractions(r)%numCartesianOrbital * &
  !                    contractions(2)%contractions(s)%numCartesianOrbital
                
  !               !! For (ab|rs)  ---> RESTRICTION a>b && r>s && r+s > a+b
  !               aux = 0
  !               order = 0

  !               !! permuted index
  !               aa = a
  !               bb = b
  !               rr = r
  !               ss = s
                
  !               !!pointer to permuted index under a not permuted loop
  !               pi => i
  !               pj => j
  !               pk => k
  !               pl => l
                
  !               !!pointer to not permuted index under a permuted loop
  !               poi => i
  !               poj => j
  !               pok => k
  !               pol => l

  !               !!Pointer to specie ID
  !               pSpecieID => specieID
  !               pOtherSpecieID => otherSpecieID
                
  !               if (contractions(1)%contractions(a)%angularMoment < contractions(1)%contractions(b)%angularMoment) then

  !                  aa = b
  !                  bb = a
                   
  !                  pi => j
  !                  pj => i

  !                  poi => j
  !                  poj => i

  !                  order = order + 1
  !               end if
                
  !               if (contractions(2)%contractions(r)%angularMoment < contractions(2)%contractions(s)%angularMoment) then
                   
  !                  rr = s
  !                  ss = r

  !                  pk => l
  !                  pl => k

  !                  pok => l
  !                  pol => k

  !                  order = order + 3

  !               end if

  !               if((contractions(1)%contractions(a)%angularMoment + contractions(1)%contractions(b)%angularMoment) > &
  !                    (contractions(2)%contractions(r)%angularMoment + contractions(2)%contractions(s)%angularMoment)) then
                   
  !                  aux = aa
  !                  aa = rr
  !                  rr = aux

  !                  aux = bb
  !                  bb = ss
  !                  ss = aux

  !                  pSpecieID => otherSpecieID
  !                  pOtherSpecieID => specieID

  !                  select case(order)
  !                  case(0)
  !                     pi => k
  !                     pj => l
  !                     pk => i
  !                     pl => j

  !                     poi => k
  !                     poj => l
  !                     pok => i
  !                     pol => j

  !                  case(1)
  !                     pi => k
  !                     pj => l
  !                     pk => j
  !                     pl => i

  !                     poi => l
  !                     poj => k
  !                     pok => i
  !                     pol => j

  !                  case(3)
  !                     pi => l
  !                     pj => k
  !                     pk => i
  !                     pl => j

  !                     poi => k
  !                     poj => l
  !                     pok => j
  !                     pol => i

  !                  case(4)
  !                     pi => l
  !                     pj => k
  !                     pk => j
  !                     pl => i

  !                     poi => l
  !                     poj => k
  !                     pok => j
  !                     pol => i

  !                  end select

  !                  order = order + 5

  !               end if

  !               !!************************************
  !               !! Calculate iteratively contractions
  !               !!

  !               !!Distancias AB, CD
  !               AB = contr(pSpecieID,aa)%origin - contr(pSpecieID,bb)%origin
  !               CD = contr(pOtherSpecieID,rr)%origin - contr(pOtherSpecieID,ss)%origin

  !               AB2 = dot_product(AB, AB)
  !               CD2 = dot_product(CD, CD)

  !               !! Asigna valores a la estrucutra Libint
  !               RepulsionDerivatives_instance%libint%AB = AB
  !               RepulsionDerivatives_instance%libint%CD = CD

  !               !!start :)
  !               integralsValue(1:arraySize) = 0.0_8

  !               !! not-permuted loop
  !               do l = 1, contractions(otherSpecieID)%contractions(s)%length
  !                  do k = 1, contractions(otherSpecieID)%contractions(r)%length
  !                     do j = 1, contractions(specieID)%contractions(b)%length
  !                        do i = 1, contractions(specieID)%contractions(a)%length

  !                           !!LIBINT PRIMQUARTET

  !                           !!Exponentes
                            
  !                           zeta = contr(pSpecieID,aa)%orbitalExponents(pi) + &
  !                                contr(pSpecieID,bb)%orbitalExponents(pj)

  !                           eta =  contr(pOtherSpecieID,rr)%orbitalExponents(pk) + &
  !                                contr(pOtherSpecieID,ss)%orbitalExponents(pl)

  !                           rho  = (zeta * eta) / (zeta + eta) !Exponente reducido ABCD

  !                           !!prim_data.U
  !                           P  = (( contr(pSpecieID,aa)%orbitalExponents(pi) * contr(pSpecieID,aa)%origin ) + &
  !                                ( contr(pSpecieID,bb)%orbitalExponents(pj) * contr(pSpecieID,bb)%origin )) / zeta

  !                           Q  = (( contr(pOtherSpecieID,rr)%orbitalExponents(pk) * contr(pOtherSpecieID,rr)%origin ) + &
  !                                ( contr(pOtherSpecieID,ss)%orbitalExponents(pl) * contr(pOtherSpecieID,ss)%origin )) / eta
                            
  !                           W  = ((zeta * P) + (eta * Q)) / (zeta + eta)
                            
  !                           PQ = P - Q
                            
  !                           primitiveQuartet%U(1:3,1)= (P - contr(pSpecieID,aa)%origin)
  !                           primitiveQuartet%U(1:3,3)= (Q - contr(pOtherSpecieID,rr)%origin)
  !                           primitiveQuartet%U(1:3,5)= (W - P)
  !                           primitiveQuartet%U(1:3,6)= (W - Q)

  !                           !!Distancias ABCD(PQ2)
  !                           PQ2 = dot_product(PQ, PQ)

  !                           !!Evalua el argumento de la funcion gamma incompleta
  !                           incompletGammaArgument = rho*PQ2

  !                           !!Overlap Factor
  !                           s12 = ((Math_PI/zeta)**1.5_8) * exp(-(( contr(pSpecieID,aa)%orbitalExponents(pi) * &
  !                                contr(pSpecieID,bb)%orbitalExponents(pj)) / zeta) * AB2)

  !                           s34 = ((Math_PI/ eta)**1.5_8) * exp(-(( contr(pOtherSpecieID,rr)%orbitalExponents(pk) * &
  !                                contr(pOtherSpecieID,ss)%orbitalExponents(pl)) /  eta) * CD2)

  !                           s1234 = sqrt(rho/Math_PI) * s12 * s34

  !                           call Math_fgamma0(sumAngularMoment,incompletGammaArgument,incompletGamma(0:sumAngularMoment))

  !                           !!prim_data.F
  !                           primitiveQuartet%F(1:sumAngularMoment+1) = 2.0_8 * incompletGamma(0:sumAngularMoment) * s1234

  !                           !!Datos restantes para prim.U
  !                           primitiveQuartet%oo2z = (0.5_8 / zeta)
  !                           primitiveQuartet%oo2n = (0.5_8 / eta)
  !                           primitiveQuartet%oo2zn = (0.5_8 / (zeta + eta))
  !                           primitiveQuartet%poz = (rho/zeta)
  !                           primitiveQuartet%pon = (rho/eta)
  !                           primitiveQuartet%oo2p = (0.5_8 / rho)

  !                           if(arraySize == 1) then
                               
  !                              auxIntegrals(1) = primitiveQuartet%F(1)

  !                           else


  !                              arraySsize(1) = arraySize
                               
  !                              RepulsionDerivatives_instance%libint%PrimQuartet = c_loc(primitiveQuartet)

  !                              !! calculate integrals (finally)                               
  !                              call c_f_procpointer(build_eri( contr(pOtherSpecieID,ss)%angularMoment , &
  !                                   contr(pOtherSpecieID,rr)%angularMoment , &
  !                                   contr(pSpecieID,bb)%angularMoment , &
  !                                   contr(pSpecieID,aa)%angularMoment), pBuild)
                               
  !                              resultPc = pBuild(RepulsionDerivatives_instance%libint,1)

  !                              !! Get numbers from memory!
  !                              call c_f_pointer(resultPc, temporalPtr, arraySsize)
                               
  !                              integralsPtr => temporalPtr
                               
  !                              !! Copy to fortran casting
  !                              auxIntegrals(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy

  !                           end if !!done by contractions

  !                           !!Normalize by primitive
  !                           m = 0
  !                           do ii = 1, contr(pSpecieID,aa)%numCartesianOrbital
  !                              do jj = 1, contr(pSpecieID,bb)%numCartesianOrbital
  !                                 do kk = 1, contr(pOtherSpecieID,rr)%numCartesianOrbital
  !                                    do ll = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital
  !                                       m = m + 1
  !                                       auxIntegrals(m) = auxIntegrals(m) &
  !                                            * contr(pSpecieID,aa)%primNormalization(pi,ii) &
  !                                            * contr(pSpecieID,bb)%primNormalization(pj,jj) &
  !                                            * contr(pOtherSpecieID,rr)%primNormalization(pk,kk) &
  !                                            * contr(pOtherSpecieID,ss)%primNormalization(pl,ll)
  !                                    end do
  !                                 end do
  !                              end do
  !                           end do !! done by primitives
                            
  !                           auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
  !                                * contr(pSpecieID,aa)%contractionCoefficients(pi) &
  !                                * contr(pSpecieID,bb)%contractionCoefficients(pj) &
  !                                * contr(pOtherSpecieID,rr)%contractionCoefficients(pk) &
  !                                * contr(pOtherSpecieID,ss)%contractionCoefficients(pl)

  !                           integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)

  !                        end do
  !                     end do
  !                  end do
  !               end do !!done Integral of pair of shells

  !               !!normalize by contraction
  !               m = 0
  !               do ii = 1,  contr(pSpecieID,aa)%numCartesianOrbital
  !                  do jj = 1,  contr(pSpecieID,bb)%numCartesianOrbital
  !                     do kk = 1,  contr(pOtherSpecieID,rr)%numCartesianOrbital
  !                        do ll = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital
  !                           m = m + 1

  !                           integralsValue(m) = integralsValue(m) &
  !                                * contr(pSpecieID,aa)%contNormalization(ii) &
  !                                * contr(pSpecieID,bb)%contNormalization(jj) &
  !                                * contr(pOtherSpecieID,rr)%contNormalization(kk) &
  !                                * contr(pOtherSpecieID,ss)%contNormalization(ll)
  !                        end do
  !                     end do
  !                  end do
  !               end do !! done by cartesian of contractions
                
  !               !!write to disk
  !               m = 0
  !               do i = 1, contr(pSpecieID,aa)%numCartesianOrbital
  !                  do j = 1, contr(pSpecieID,bb)%numCartesianOrbital
  !                     do k = 1, contr(pOtherSpecieID,rr)%numCartesianOrbital
  !                        do l = 1, contr(pOtherSpecieID,ss)%numCartesianOrbital

  !                           m = m + 1

  !                           !! index not permuted
  !                           pa=labelsOfContractions(a)+poi-1
  !                           pb=labelsOfContractions(b)+poj-1
  !                           pr=otherLabelsOfContractions(r)+pok-1
  !                           ps=otherLabelsOfContractions(s)+pol-1

  !                           if( pa <= pb .and. pr <= ps ) then

  !                              if(abs(integralsValue(m)) > 1.0D-10) then
                                  
  !                                 counter = counter + 1
  !                                 auxCounter = auxCounter + 1
                                  
  !                                 eris%a(counter) = pa
  !                                 eris%b(counter) = pb
  !                                 eris%c(counter) = pr
  !                                 eris%d(counter) = ps
  !                                 eris%integrals(counter) = integralsValue(m)

  !                              end if

                               
  !                              if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

  !                                 write(34) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                                      eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                                      eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
                                  
  !                                 counter = 0

  !                              end if
  !                           end if !! Stack control

  !                        end do
  !                     end do
  !                  end do
  !               end do !! Done write to disk
                
  !            end do
  !            u=r+1
  !         end do
  !      end do
  !   end do !! done by basis set

  !   eris%a(counter+1) = -1
  !   eris%b(counter+1) = -1
  !   eris%c(counter+1) = -1
  !   eris%d(counter+1) = -1
  !   eris%integrals(counter+1) = 0.0_8	       
    
  !   write(34) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !        eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !   specieID = auxSpecieID
  !   otherSpecieID = auxOtherSpecieID
    
  !   if(CONTROL_instance%LAST_STEP) then
  !      write(*,"(A,I12,A,A)") " Stored ", &
  !           auxCounter, &
  !           " non-zero repulsion integrals between species: ", &
  !           trim(MolecularSystem_instance%species(specieID)%name)//" / "//&
  !           trim(MolecularSystem_instance%species(otherSpecieID)%name)
  !   end if

  !   close(34)
    
  ! end subroutine RepulsionDerivatives_computeInterSpecies

  !>
  !! @brief Indica si el objeto ha sido instanciado o no
  
  function RepulsionDerivatives_isInstanced( ) result( output )
    implicit  none

    logical :: output
    
    output = RepulsionDerivatives_instance%isInstanced
    
  end function RepulsionDerivatives_isInstanced  

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine RepulsionDerivatives_exception( typeMessage, description, debugDescription)
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

  end subroutine RepulsionDerivatives_exception

end module RepulsionDerivatives_
