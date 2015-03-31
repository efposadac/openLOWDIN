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

#define contr(n,m) contractions(n)%contractions(m)
  
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

  !> @brief type for convenence
  type, private :: auxBasis
     type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie
  end type auxBasis

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


  !>
  !! @brief calculate derivative eris using libint library for all basis set (intra-specie)
  !! @author J.M. Rodas 2015
  !! @version 2.0
  !! @info Tested
  subroutine RepulsionDerivatives_getDerive(this, a, b, r, s, deriveValue, specieID)
    implicit none
    type(ContractedGaussian), intent(in):: this(:)
    integer, intent(in) :: a, b, r, s, specieID
    real(8), allocatable :: deriveValue(:)
    
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: numberOfPrimitives
    integer :: maxNumberOfPrimitives 
    integer :: maxNumberOfCartesians 
    integer :: maxNPrimSize, maxNCartSize
    logical :: p12, p34, p13p24
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT) 
    integer :: am1, am2, am3, am4
    integer :: pa, pb, pr, ps !< labels index
    real(8) :: auxExponentA, auxCoefficientA, auxContConstantA, auxPrimConstantA, c1
    real(8) :: auxExponentB, auxCoefficientB, auxContConstantB, auxPrimConstantB, c2
    real(8) :: auxExponentR, auxCoefficientR, auxContConstantR, auxPrimConstantR, c3
    real(8) :: auxExponentS, auxCoefficientS, auxContConstantS, auxPrimConstantS, c4
    real(8), allocatable :: workForces(:,:), work(:,:)
    integer :: aux !<auxiliary index
    integer :: arraySsize(1)
    integer :: arraySize
    integer :: counter, auxCounter
    integer(8) :: control
    integer,target :: i, j, k, l !< contraction length iterators
    integer :: nc1, nc2, nc3, nc4
    real(8), dimension(:), pointer :: derivativesPtr !< pointer to C array of derivatives
    real(c_double), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: incompletGamma(:) !<array with incomplete gamma integrals

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !<geometric values that appear in gaussian product
    real(8) :: zeta, ooz, oo2z, nu, oon, oo2n, oo2zn, rho !< exponents... that appear in gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !<geometric quantities that appear in gaussian product
    real(8) :: incompletGammaArgument
    integer :: bufferOffsets(0:3)
    type(prim_data), target :: primitiveQuartet !<Primquartet object needed by LIBINT
    type(c_ptr) :: resultPc !< array of integrals from C (LIBINT)
    procedure(RepulsionDerivatives_buildLibDeriv), pointer :: pBuild !<procedure to calculate eris on LIBINT !x
    
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    maxNumberOfPrimitives = MolecularSystem_getMaxNumberofPrimitives(specieID)
    maxNumberOfCartesians = MolecularSystem_getMaxNumberofCartesians(specieID)
    maxNPrimSize = maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives
    maxNCartSize = maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians

    arraySsize(1) = this(a)%numCartesianOrbital * &
         this(b)%numCartesianOrbital * &
         this(r)%numCartesianOrbital * &
         this(s)%numCartesianOrbital


    arraySize = arraySsize(1)
    sumAngularMoment = this(a)%angularMoment + this(b)%angularMoment + this(r)%angularMoment + this(s)%angularMoment

    if(allocated(workForces)) deallocate(workForces)
    allocate(workForces(arraySize,12))

    if(allocated(work)) deallocate(work)
    allocate(work(arraySize,12))

    work = 0.0_8

    if(allocated(incompletGamma)) deallocate(incompletGamma)
    allocate(incompletGamma(0:maxAngularMoment+1))
    ! if(allocated(incompletGamma)) deallocate(incompletGamma)
    ! allocate(incompletGamma(0:sumAngularMoment+1))

    ! Libderiv constructor (solo una vez)
    if( RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_destructor()
    end if
    
    if( .not. RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_constructor(sumAngularMoment,maxNPrimSize,maxNCartSize)
       !! DEBUG
       !call RepulsionDerivatives_show()
    end if

    aa = a
    bb = b
    rr = r
    ss = s

    p12 = .false.
    p34 = .false.
    p13p24 = .false.
    
    if(this(a)%angularMoment .lt. this(b)%angularMoment) then
       aa = b
       bb = a

       p12 = .true.
    end if

    if(this(r)%angularMoment .lt. this(s)%angularMoment) then
       rr = s
       ss = r

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
    end if

    if(p12) then
        if(p34) then
            if(p13p24) then
                ! (AB|CD) -> (DC|BA)
                bufferOffsets(0) = 9 
                bufferOffsets(1) = 6
                bufferOffsets(2) = 3 
                bufferOffsets(3) = 0
            else
                ! (AB|CD) -> (BA|DC)
                bufferOffsets(0) = 3 
                bufferOffsets(1) = 0
                bufferOffsets(2) = 9 
                bufferOffsets(3) = 6
            end if
        else
            if(p13p24) then
                ! (AB|CD) -> (CD|BA)
                bufferOffsets(0) = 6 
                bufferOffsets(1) = 9
                bufferOffsets(2) = 3 
                bufferOffsets(3) = 0
            else
                ! (AB|CD) -> (BA|CD)
                bufferOffsets(0) = 3 
                bufferOffsets(1) = 0
                bufferOffsets(2) = 6 
                bufferOffsets(3) = 9
            end if
        end if
    else
        if(p34) then
            if(p13p24) then
                ! (AB|CD) -> (DC|AB)
                bufferOffsets(0) = 9 
                bufferOffsets(1) = 6
                bufferOffsets(2) = 0 
                bufferOffsets(3) = 3
            else
                ! (AB|CD) -> (AB|DC)
                bufferOffsets(0) = 0 
                bufferOffsets(1) = 3
                bufferOffsets(2) = 9 
                bufferOffsets(3) = 6
            end if
        else
            if(p13p24) then
                ! (AB|CD) -> (CD|AB)
                bufferOffsets(0) = 6 
                bufferOffsets(1) = 9
                bufferOffsets(2) = 0 
                bufferOffsets(3) = 3
            else
                ! (AB|CD) -> (AB|CD)
                bufferOffsets(0) = 0 
                bufferOffsets(1) = 3
                bufferOffsets(2) = 6 
                bufferOffsets(3) = 9
            end if
        end if
    end if

    ! write(*,"(A)") "Contraida:"
    ! write(*,"(3(A,I))") "maxam: ",  maxAngularMoment, " maxnprim: ", maxNPrimSize, " maxncart: ", maxNCartSize
    !write(*,"(A,I,A,I,A,I,A,I,A)") "(",aa,",",bb,"|",rr,",",ss,")"

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

                s12 = ((Math_PI*ooz)**1.5_8) * exp(-auxExponentA*auxExponentB*ooz*AB2)*c1*c2

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
                call Math_fgamma0(maxAngularMoment+1,incompletGammaArgument,incompletGamma(0:maxAngularMoment+1))
                ! call Math_fgamma0(sumAngularMoment+1,incompletGammaArgument,incompletGamma(0:sumAngularMoment+1))

                s34 = ((Math_PI*oon)**1.5_8) * exp(-auxExponentR*auxExponentS*oon*CD2)*c3*c4

                s1234 = 2.0_8*sqrt(rho/Math_PI) * s12 * s34

                ! do i=1, sumAngularMoment+1
                !    primitiveQuartet%F(i) = incompletGamma(i-1)*s1234
                ! end do

                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                ! do i=1,3
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,1)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,2)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,3)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,4)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,5)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,6)
                ! end do

                do i=0, maxAngularMoment+1
                   primitiveQuartet%F(i+1) = incompletGamma(i)*s1234
                   ! write(*,"(A1,I1,A,f17.12)") "F", i, " : ", primitiveQuartet%F(i+1)
                end do
                ! write(*,"(A)") "-----------------------------------------"
                ! arraySsize(1) = arraySize

                RepulsionDerivatives_instance%libderiv%PrimQuartet = c_loc(primitiveQuartet) !ok

                ! write(*,"(A)") "-----------------------------------------"
                ! ! write(*,"(A,I8)") "Array Size: ", arraySsize(1)
                ! ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",pa,",",pb,"|",pr,",",ps,")"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                ! write(*,"(A,f17.12)") "oo2z : ", primitiveQuartet%oo2z
                ! write(*,"(A,f17.12)") "oo2n : ", primitiveQuartet%oo2n 
                ! write(*,"(A,f17.12)") "oo2zn : ", primitiveQuartet%oo2zn
                ! write(*,"(A,f17.12)") "poz : ", primitiveQuartet%poz
                ! write(*,"(A,f17.12)") "pon : ", primitiveQuartet%pon
                ! write(*,"(A,f22.8)") "a : ", primitiveQuartet%twozeta_a 
                ! write(*,"(A,f22.8)") "b : ", primitiveQuartet%twozeta_b 
                ! write(*,"(A,f22.8)") "c : ", primitiveQuartet%twozeta_c 
                ! write(*,"(A,f22.8)") "d : ", primitiveQuartet%twozeta_d 
                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2z
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2n
                ! write(*,"(A)") "-----------------------------------------"
                ! RepulsionDerivatives_instance%libderiv%PrimQuartet(numberOfPrimitives) = c_loc(primitiveQuartet)
                call c_f_procpointer(build_deriv1_eri(am4,am3,am2,am1), pBuild)

                call pBuild(RepulsionDerivatives_instance%libderiv,1)


                workForces = 0.0_8
                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                do k=1,12
                   if(k==4 .or. k==5 .or. k==6) cycle
                    resultPc = RepulsionDerivatives_instance%libderiv%ABCD(k)
                    call c_f_pointer(resultPc, temporalPtr, arraySsize)
                   
                    derivativesPtr => temporalPtr
                    
                    do i=1,arraySsize(1)
                       workForces(i,k) = derivativesPtr(i)
                       ! write(*,"(2I, 2f17.12)") k, i, derivativesPtr(i), workForces(i,k)
                    end do
                 end do

                 do k=4,6
                    do i=1,arraySsize(1)
                       workForces(i,k) = 0.0_8
                    end do
                 end do
                 
                
                 do k=1,12
                    do i=1, arraySsize(1)
                       work(i,k) = work(i,k) + workForces(i,k)
                    end do
                 end do
                 ! write(*,"(A)") "-----------------------------------------"
                 numberOfPrimitives = numberOfPrimitives + 1
             end do
          end do
       end do
    end do

    call RepulsionDerivatives_reordering1(bufferOffsets, work, deriveValue, arraySsize(1))

    if (p12 .or. p34 .or. p13p24) then
       nc1 =  this(aa)%numCartesianOrbital
       nc2 =  this(bb)%numCartesianOrbital
       nc3 =  this(rr)%numCartesianOrbital
       nc4 =  this(ss)%numCartesianOrbital
       call RepulsionDerivatives_permute(nc1,nc2,nc3,nc4, p12, p34, p13p24, arraySsize(1), deriveValue)
    end if


    ! if (p12 .or. p34 .or. p13p24) then
    !    call RepulsionDerivatives_permute(this, aa, bb, rr, ss, p12, p34, p13p24, arraySsize(1), deriveValue)
    ! end if

    ! write(*,"(A)") "-----------------------------------------"
    ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",a,",",b,"|",r,",",s,")"
    ! do k=0,9*arraySsize(1) - 1
    !    write(*,"(I2,f17.12)") k, deriveValue(k)
    ! end do
    ! ! do k=1,12
    ! !    do i=1,arraySsize(1)
    ! !       write(*,"(I2,I2,f17.12)") k, 1, work(i,k)
    ! !    end do
    ! ! end do
    ! write(*,"(A)") "-----------------------------------------"

  end subroutine RepulsionDerivatives_getDerive

  !>
  !! @brief calculate derivative eris using libint library for all basis set (inter-specie)
  !! @author J.M. Rodas 2015
  !! @version 2.0
  !! @info Tested
  subroutine RepulsionDerivatives_getInterDerive(a, b, r, s, deriveValue, specieID, otherSpecieID)
    implicit none
    integer, intent(in) :: a, b, r, s
    integer :: specieID
    integer :: otherSpecieID
    real(8), allocatable :: deriveValue(:)
    
    type(auxBasis) :: contractions(2)
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: numberOfPrimitives
    integer :: maxNumberOfPrimitives 
    integer :: maxNumberOfCartesians 
    integer :: maxNPrimSize, maxNCartSize
    logical :: p12, p34, p13p24
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT) 
    integer :: am1, am2, am3, am4
    integer :: pa, pb, pr, ps !< labels index
    real(8) :: auxExponentA, auxCoefficientA, auxContConstantA, auxPrimConstantA, c1
    real(8) :: auxExponentB, auxCoefficientB, auxContConstantB, auxPrimConstantB, c2
    real(8) :: auxExponentR, auxCoefficientR, auxContConstantR, auxPrimConstantR, c3
    real(8) :: auxExponentS, auxCoefficientS, auxContConstantS, auxPrimConstantS, c4
    real(8), allocatable :: workForces(:,:), work(:,:)
    integer :: aux !<auxiliary index
    integer :: arraySsize(1)
    integer :: arraySize
    integer :: counter, auxCounter
    integer(8) :: control
    integer,target :: i, j, k, l !< contraction length iterators
    integer :: nc1, nc2, nc3, nc4
    integer :: pSpecieID, pOtherSpecieID !! pointer to species ID

    real(8), dimension(:), pointer :: derivativesPtr !< pointer to C array of derivatives
    real(c_double), dimension(:), pointer :: temporalPtr

    real(8), allocatable :: incompletGamma(:) !<array with incomplete gamma integrals

    real(8) :: P(3), Q(3), W(3), AB(3), CD(3), PQ(3) !<geometric values that appear in gaussian product
    real(8) :: zeta, ooz, oo2z, nu, oon, oo2n, oo2zn, rho !< exponents... that appear in gaussian product
    real(8) :: s1234, s12, s34, AB2, CD2, PQ2 !<geometric quantities that appear in gaussian product
    real(8) :: incompletGammaArgument
    integer :: bufferOffsets(0:3)
    type(prim_data), target :: primitiveQuartet !<Primquartet object needed by LIBINT
    type(c_ptr) :: resultPc !< array of integrals from C (LIBINT)
    procedure(RepulsionDerivatives_buildLibDeriv), pointer :: pBuild !<procedure to calculate eris on LIBINT !x
    
    maxAngularMoment = max(MolecularSystem_getMaxAngularMoment(specieID), MolecularSystem_getMaxAngularMoment(otherSpecieID))
    maxNumberOfPrimitives = max(MolecularSystem_getMaxNumberofPrimitives(specieID), MolecularSystem_getMaxNumberofPrimitives(otherSpecieID))
    maxNumberOfCartesians = max(MolecularSystem_getMaxNumberofCartesians(specieID), MolecularSystem_getMaxNumberofCartesians(otherSpecieID))
    maxNPrimSize = maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives*maxNumberOfPrimitives
    maxNCartSize = maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians*maxNumberOfCartesians

    call MolecularSystem_getBasisSet(specieID, contractions(1)%contractions)
    call MolecularSystem_getBasisSet(otherSpecieID, contractions(2)%contractions)

    arraySsize(1) =contractions(1)%contractions(a)%numCartesianOrbital * &
         contractions(1)%contractions(b)%numCartesianOrbital * &
         contractions(2)%contractions(r)%numCartesianOrbital * &
         contractions(2)%contractions(s)%numCartesianOrbital

    arraySize = arraySsize(1)
    
    sumAngularMoment =  contractions(1)%contractions(a)%angularMoment + &
         contractions(1)%contractions(b)%angularMoment + &
         contractions(2)%contractions(r)%angularMoment + &
         contractions(2)%contractions(s)%angularMoment

    if(allocated(workForces)) deallocate(workForces)
    allocate(workForces(arraySize,12))

    if(allocated(work)) deallocate(work)
    allocate(work(arraySize,12))

    work = 0.0_8

    if(allocated(incompletGamma)) deallocate(incompletGamma)
    allocate(incompletGamma(0:maxAngularMoment+1))
    ! if(allocated(incompletGamma)) deallocate(incompletGamma)
    ! allocate(incompletGamma(0:sumAngularMoment+1))

    ! Libderiv constructor (solo una vez)
    if( RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_destructor()
    end if
    
    if( .not. RepulsionDerivatives_isInstanced() ) then
       call RepulsionDerivatives_constructor(sumAngularMoment,maxNPrimSize,maxNCartSize)
       !! DEBUG
       !call RepulsionDerivatives_show()
    end if

    aa = a
    bb = b
    rr = r
    ss = s

    p12 = .false.
    p34 = .false.
    p13p24 = .false.

    pSpecieID = 1
    pOtherSpecieID = 2
    
    if(contractions(1)%contractions(a)%angularMoment .lt. contractions(1)%contractions(b)%angularMoment) then
       aa = b
       bb = a

       p12 = .true.
    end if

    if(contractions(2)%contractions(r)%angularMoment .lt. contractions(2)%contractions(s)%angularMoment) then
       rr = s
       ss = r

       p34 = .true.
    end if

    if((contractions(1)%contractions(a)%angularMoment + contractions(1)%contractions(b)%angularMoment) .gt. &
         (contractions(2)%contractions(r)%angularMoment + contractions(2)%contractions(s)%angularMoment)) then
       aux = aa
       aa = rr
       rr = aux

       aux = bb
       bb = ss
       ss = aux


       pSpecieID = 2
       pOtherSpecieID = 1
       p13p24 = .true.
    end if

    if(p12) then
        if(p34) then
            if(p13p24) then
                ! (AB|CD) -> (DC|BA)
                bufferOffsets(0) = 9 
                bufferOffsets(1) = 6
                bufferOffsets(2) = 3 
                bufferOffsets(3) = 0
            else
                ! (AB|CD) -> (BA|DC)
                bufferOffsets(0) = 3 
                bufferOffsets(1) = 0
                bufferOffsets(2) = 9 
                bufferOffsets(3) = 6
            end if
        else
            if(p13p24) then
                ! (AB|CD) -> (CD|BA)
                bufferOffsets(0) = 6 
                bufferOffsets(1) = 9
                bufferOffsets(2) = 3 
                bufferOffsets(3) = 0
            else
                ! (AB|CD) -> (BA|CD)
                bufferOffsets(0) = 3 
                bufferOffsets(1) = 0
                bufferOffsets(2) = 6 
                bufferOffsets(3) = 9
            end if
        end if
    else
        if(p34) then
            if(p13p24) then
                ! (AB|CD) -> (DC|AB)
                bufferOffsets(0) = 9 
                bufferOffsets(1) = 6
                bufferOffsets(2) = 0 
                bufferOffsets(3) = 3
            else
                ! (AB|CD) -> (AB|DC)
                bufferOffsets(0) = 0 
                bufferOffsets(1) = 3
                bufferOffsets(2) = 9 
                bufferOffsets(3) = 6
            end if
        else
            if(p13p24) then
                ! (AB|CD) -> (CD|AB)
                bufferOffsets(0) = 6 
                bufferOffsets(1) = 9
                bufferOffsets(2) = 0 
                bufferOffsets(3) = 3
            else
                ! (AB|CD) -> (AB|CD)
                bufferOffsets(0) = 0 
                bufferOffsets(1) = 3
                bufferOffsets(2) = 6 
                bufferOffsets(3) = 9
            end if
        end if
    end if

    ! write(*,"(A)") "Contraida:"
    ! write(*,"(3(A,I))") "maxam: ",  maxAngularMoment, " maxnprim: ", maxNPrimSize, " maxncart: ", maxNCartSize
    !write(*,"(A,I,A,I,A,I,A,I,A)") "(",aa,",",bb,"|",rr,",",ss,")"

    am1 = contr(pSpecieID,aa)%angularMoment
    am2 = contr(pSpecieID,bb)%angularMoment
    am3 = contr(pOtherSpecieID,rr)%angularMoment
    am4 = contr(pOtherSpecieID,ss)%angularMoment
   
    !!Distancias AB, CD
    AB = contr(pSpecieID,aa)%origin - contr(pSpecieID,bb)%origin
    CD = contr(pOtherSpecieID,rr)%origin - contr(pOtherSpecieID,ss)%origin

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

    do pa=1, contr(pSpecieID,aa)%length
       do pb=1, contr(pSpecieID,bb)%length
          do pr=1, contr(pOtherSpecieID,rr)%length
             do ps=1, contr(pOtherSpecieID,ss)%length
                auxExponentS = contr(pOtherSpecieID,ss)%orbitalExponents(ps)
                auxCoefficientS = contr(pOtherSpecieID,ss)%contractionCoefficients(ps)
                auxPrimConstantS = contr(pOtherSpecieID,ss)%primNormalization(ps,1)
                auxContConstantS = contr(pOtherSpecieID,ss)%contNormalization(1)
                c4 = auxCoefficientS*auxPrimConstantS*auxContConstantS

                auxExponentR = contr(pOtherSpecieID,rr)%orbitalExponents(pr)
                auxCoefficientR = contr(pOtherSpecieID,rr)%contractionCoefficients(pr)
                auxPrimConstantR = contr(pOtherSpecieID,rr)%primNormalization(pr,1)
                auxContConstantR = contr(pOtherSpecieID,rr)%contNormalization(1)
                c3 = auxCoefficientR*auxPrimConstantR*auxContConstantR

                auxExponentB = contr(pSpecieID,bb)%orbitalExponents(pb)
                auxCoefficientB = contr(pSpecieID,bb)%contractionCoefficients(pb)
                auxPrimConstantB = contr(pSpecieID,bb)%primNormalization(pb,1)
                auxContConstantB = contr(pSpecieID,bb)%contNormalization(1)
                c2 = auxCoefficientB*auxPrimConstantB*auxContConstantB


                auxExponentA = contr(pSpecieID,aa)%orbitalExponents(pa)
                auxCoefficientA = contr(pSpecieID,aa)%contractionCoefficients(pa)
                auxPrimConstantA = contr(pSpecieID,aa)%primNormalization(pa,1)
                auxContConstantA = contr(pSpecieID,aa)%contNormalization(1)
                c1 = auxCoefficientA*auxPrimConstantA*auxContConstantA

                zeta = auxExponentA + auxExponentB
                ooz = 1.0_8/zeta
                oo2z = 1.0_8/(2.0_8*zeta)

                P = (auxExponentA*contr(pSpecieID,aa)%origin + auxExponentB*contr(pSpecieID,bb)%origin)*ooz

                primitiveQuartet%U(1:3,1)= P - contr(pSpecieID,aa)%origin
                primitiveQuartet%U(1:3,2)= P - contr(pSpecieID,bb)%origin

                s12 = ((Math_PI*ooz)**1.5_8) * exp(-auxExponentA*auxExponentB*ooz*AB2)*c1*c2

                nu = auxExponentR + auxExponentS
                oon = 1.0_8/nu
                oo2n = 1.0_8/(2.0_8*nu)
                oo2zn = 1.0_8/(2.0_8*(zeta+nu))
                rho = (zeta*nu)/(zeta+nu)

                Q = (auxExponentR*contr(pOtherSpecieID,rr)%origin + auxExponentS*contr(pOtherSpecieID,ss)%origin)*oon
                primitiveQuartet%U(1:3,3)= Q - contr(pOtherSpecieID,rr)%origin
                primitiveQuartet%U(1:3,4)= Q - contr(pOtherSpecieID,ss)%origin

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
                call Math_fgamma0(maxAngularMoment+1,incompletGammaArgument,incompletGamma(0:maxAngularMoment+1))
                ! call Math_fgamma0(sumAngularMoment+1,incompletGammaArgument,incompletGamma(0:sumAngularMoment+1))

                s34 = ((Math_PI*oon)**1.5_8) * exp(-auxExponentR*auxExponentS*oon*CD2)*c3*c4

                s1234 = 2.0_8*sqrt(rho/Math_PI) * s12 * s34

                ! do i=1, sumAngularMoment+1
                !    primitiveQuartet%F(i) = incompletGamma(i-1)*s1234
                ! end do

                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                ! do i=1,3
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,1)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,2)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,3)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,4)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,5)
                !    write(*,"(A1,I1,A,f17.12)") "U", i, " , ", primitiveQuartet%U(i,6)
                ! end do

                do i=0, maxAngularMoment+1
                   primitiveQuartet%F(i+1) = incompletGamma(i)*s1234
                   ! write(*,"(A1,I1,A,f17.12)") "F", i, " : ", primitiveQuartet%F(i+1)
                end do
                ! write(*,"(A)") "-----------------------------------------"
                ! arraySsize(1) = arraySize

                RepulsionDerivatives_instance%libderiv%PrimQuartet = c_loc(primitiveQuartet) !ok

                ! write(*,"(A)") "-----------------------------------------"
                ! ! write(*,"(A,I8)") "Array Size: ", arraySsize(1)
                ! ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",pa,",",pb,"|",pr,",",ps,")"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                ! write(*,"(A,f17.12)") "oo2z : ", primitiveQuartet%oo2z
                ! write(*,"(A,f17.12)") "oo2n : ", primitiveQuartet%oo2n 
                ! write(*,"(A,f17.12)") "oo2zn : ", primitiveQuartet%oo2zn
                ! write(*,"(A,f17.12)") "poz : ", primitiveQuartet%poz
                ! write(*,"(A,f17.12)") "pon : ", primitiveQuartet%pon
                ! write(*,"(A,f22.8)") "a : ", primitiveQuartet%twozeta_a 
                ! write(*,"(A,f22.8)") "b : ", primitiveQuartet%twozeta_b 
                ! write(*,"(A,f22.8)") "c : ", primitiveQuartet%twozeta_c 
                ! write(*,"(A,f22.8)") "d : ", primitiveQuartet%twozeta_d 
                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2z
                ! write(*,"(A,3(F17.12))") "AB deriv: ", RepulsionDerivatives_instance%libderiv%PrimQuartet%oo2n
                ! write(*,"(A)") "-----------------------------------------"
                ! RepulsionDerivatives_instance%libderiv%PrimQuartet(numberOfPrimitives) = c_loc(primitiveQuartet)
                call c_f_procpointer(build_deriv1_eri(am4,am3,am2,am1), pBuild)

                call pBuild(RepulsionDerivatives_instance%libderiv,1)


                workForces = 0.0_8
                ! write(*,"(A)") "-----------------------------------------"
                ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",aa,",",bb,"|",rr,",",ss,")"
                do k=1,12
                   if(k==4 .or. k==5 .or. k==6) cycle
                    resultPc = RepulsionDerivatives_instance%libderiv%ABCD(k)
                    call c_f_pointer(resultPc, temporalPtr, arraySsize)
                   
                    derivativesPtr => temporalPtr
                    
                    do i=1,arraySsize(1)
                       workForces(i,k) = derivativesPtr(i)
                       ! write(*,"(2I, 2f17.12)") k, i, derivativesPtr(i), workForces(i,k)
                    end do
                 end do

                 do k=4,6
                    do i=1,arraySsize(1)
                       workForces(i,k) = 0.0_8
                    end do
                 end do
                 
                
                 do k=1,12
                    do i=1, arraySsize(1)
                       work(i,k) = work(i,k) + workForces(i,k)
                    end do
                 end do
                 ! write(*,"(A)") "-----------------------------------------"
                 numberOfPrimitives = numberOfPrimitives + 1
             end do
          end do
       end do
    end do

    call RepulsionDerivatives_reordering1(bufferOffsets, work, deriveValue, arraySsize(1))

    if (p12 .or. p34 .or. p13p24) then
       nc1 =  contr(pSpecieID,aa)%numCartesianOrbital
       nc2 =  contr(pSpecieID,bb)%numCartesianOrbital
       nc3 =  contr(pOtherSpecieID,rr)%numCartesianOrbital
       nc4 =  contr(pOtherSpecieID,ss)%numCartesianOrbital
       call RepulsionDerivatives_permute(nc1,nc2,nc3,nc4, p12, p34, p13p24, arraySsize(1), deriveValue)
    end if

    ! write(*,"(A)") "-----------------------------------------"
    ! write(*,"(A1,I1,A1,I1,A1,I1,A1,I1,A1)") "(",a,",",b,"|",r,",",s,")"
    ! do k=0,9*arraySsize(1) - 1
    !    write(*,"(I2,f17.12)") k, deriveValue(k)
    ! end do
    ! ! do k=1,12
    ! !    do i=1,arraySsize(1)
    ! !       write(*,"(I2,I2,f17.12)") k, 1, work(i,k)
    ! !    end do
    ! ! end do
    ! write(*,"(A)") "-----------------------------------------"

  end subroutine RepulsionDerivatives_getInterDerive
  

  subroutine RepulsionDerivatives_reordering1(bufferOffsets, prederivatives, deriveValue, size)
    implicit none
    integer, intent(in) :: bufferOffsets(0:3)
    real(8), allocatable, intent(in) :: prederivatives(:,:)
    real(8), allocatable, intent(inout) :: deriveValue(:)
    integer, intent(in) :: size
    integer :: A, B, C, D
    integer :: i


    if(allocated(deriveValue)) deallocate(deriveValue)
    allocate(deriveValue(0:(size*9 - 1)))

    deriveValue = 0.0_8

    do i=0, 3
       if(bufferOffsets(i) == 0) then
          A = 3*i + 1
       end if
       if(bufferOffsets(i) == 3) then 
          B = 3*i + 1
       end if
       if(bufferOffsets(i) == 6) then 
          C = 3*i + 1
       end if
       if(bufferOffsets(i) == 9) then 
          D = 3*i + 1
       end if
    end do
    
    if(bufferOffsets(1)==0) then
       do i = 0, size - 1
          deriveValue(i+0*size) = deriveValue(i+0*size) - prederivatives(i+1,B+0) - prederivatives(i+1,C+0) - prederivatives(i+1,D+0) !Ax
          deriveValue(i+1*size) = deriveValue(i+1*size) - prederivatives(i+1,B+1) - prederivatives(i+1,C+1) - prederivatives(i+1,D+1) !Ay
          deriveValue(i+2*size) = deriveValue(i+2*size) - prederivatives(i+1,B+2) - prederivatives(i+1,C+2) - prederivatives(i+1,D+2) !Az
          deriveValue(i+3*size) = deriveValue(i+3*size) + prederivatives(i+1,C+0) !Cx
          deriveValue(i+4*size) = deriveValue(i+4*size) + prederivatives(i+1,C+1) !Cy
          deriveValue(i+5*size) = deriveValue(i+5*size) + prederivatives(i+1,C+2) !Cz
          deriveValue(i+6*size) = deriveValue(i+6*size) + prederivatives(i+1,D+0) !Dx
          deriveValue(i+7*size) = deriveValue(i+7*size) + prederivatives(i+1,D+1) !Dy
          deriveValue(i+8*size) = deriveValue(i+8*size) + prederivatives(i+1,D+2) !Dz
       end do
    else if(bufferOffsets(1)==6) then
       do i = 0, size - 1
          deriveValue(i+0*size) = deriveValue(i+0*size) + prederivatives(i+1,A+0) !Ax
          deriveValue(i+1*size) = deriveValue(i+1*size) + prederivatives(i+1,A+1) !Ay
          deriveValue(i+2*size) = deriveValue(i+2*size) + prederivatives(i+1,A+2) !Az
          deriveValue(i+3*size) = deriveValue(i+3*size) - prederivatives(i+1,A+0) - prederivatives(i+1,B+0) - prederivatives(i+1,D+0) !Cx
          deriveValue(i+4*size) = deriveValue(i+4*size) - prederivatives(i+1,A+1) - prederivatives(i+1,B+1) - prederivatives(i+1,D+1) !Cy
          deriveValue(i+5*size) = deriveValue(i+5*size) - prederivatives(i+1,A+2) - prederivatives(i+1,B+2) - prederivatives(i+1,D+2) !Cz
          deriveValue(i+6*size) = deriveValue(i+6*size) + prederivatives(i+1,D+0) !Dx
          deriveValue(i+7*size) = deriveValue(i+7*size) + prederivatives(i+1,D+1) !Dy
          deriveValue(i+8*size) = deriveValue(i+8*size) + prederivatives(i+1,D+2) !Dz
       end do
    else if(bufferOffsets(1)==9) then
       do i = 0, size - 1
          deriveValue(i+0*size) = deriveValue(i+0*size) + prederivatives(i+1,A+0) !Ax
          deriveValue(i+1*size) = deriveValue(i+1*size) + prederivatives(i+1,A+1) !Ay
          deriveValue(i+2*size) = deriveValue(i+2*size) + prederivatives(i+1,A+2) !Az
          deriveValue(i+3*size) = deriveValue(i+3*size) + prederivatives(i+1,C+0) !Cx
          deriveValue(i+4*size) = deriveValue(i+4*size) + prederivatives(i+1,C+1) !Cy
          deriveValue(i+5*size) = deriveValue(i+5*size) + prederivatives(i+1,C+2) !Cz
          deriveValue(i+6*size) = deriveValue(i+6*size) - prederivatives(i+1,A+0) - prederivatives(i+1,B+0) - prederivatives(i+1,C+0) !Dx
          deriveValue(i+7*size) = deriveValue(i+7*size) - prederivatives(i+1,A+1) - prederivatives(i+1,B+1) - prederivatives(i+1,C+1) !Dy
          deriveValue(i+8*size) = deriveValue(i+8*size) - prederivatives(i+1,A+2) - prederivatives(i+1,B+2) - prederivatives(i+1,C+2) !Dz
       end do
    else
       do i = 0, size - 1
          deriveValue(i+0*size) = deriveValue(i+0*size) + prederivatives(i+1,A+0) !Ax
          deriveValue(i+1*size) = deriveValue(i+1*size) + prederivatives(i+1,A+1) !Ay
          deriveValue(i+2*size) = deriveValue(i+2*size) + prederivatives(i+1,A+2) !Az
          deriveValue(i+3*size) = deriveValue(i+3*size) + prederivatives(i+1,C+0) !Cx
          deriveValue(i+4*size) = deriveValue(i+4*size) + prederivatives(i+1,C+1) !Cy
          deriveValue(i+5*size) = deriveValue(i+5*size) + prederivatives(i+1,C+2) !Cz
          deriveValue(i+6*size) = deriveValue(i+6*size) + prederivatives(i+1,D+0) !Dx
          deriveValue(i+7*size) = deriveValue(i+7*size) + prederivatives(i+1,D+1) !Dy
          deriveValue(i+8*size) = deriveValue(i+8*size) + prederivatives(i+1,D+2) !Dz
       end do
    end if
  end subroutine RepulsionDerivatives_reordering1

  subroutine RepulsionDerivatives_permute(nbf1, nbf2, nbf3, nbf4, p12, p34, p13p24, size, deriveValue)
    implicit none
    ! type(ContractedGaussian), intent(in):: this(:)
    ! integer, intent(in) :: a, b, r, s
    integer, intent(in) :: nbf1, nbf2, nbf3, nbf4
    logical, intent(in) :: p12, p34, p13p24
    integer, intent(in) :: size
    real(8), allocatable, intent(inout) :: deriveValue(:)
    integer :: bf1, bf2, bf3, bf4, derivIter
    integer :: f1, f2, f3, f4
    integer :: newPtr, auxIter
    real(8), allocatable :: auxDerivatives(:)


    if(allocated(auxDerivatives)) deallocate(auxDerivatives)
    allocate(auxDerivatives(0:(size*9 - 1)))

    auxDerivatives = deriveValue
    deriveValue = 0.0_8

    ! nbf1 =  this(a)%numCartesianOrbital
    ! nbf2 =  this(b)%numCartesianOrbital
    ! nbf3 =  this(r)%numCartesianOrbital
    ! nbf4 =  this(s)%numCartesianOrbital

    do derivIter=0, 8
       if (.not.p13p24) then
          if (p12) then
             if (p34) then
                f1=nbf2
                f2=nbf1
                f3=nbf4
                f4=nbf3
                auxIter = 0
                ! (a,b|r,s) -> (b,a|s,r)
                do bf2=0, f2-1
                   do bf1=0,f1-1
                      do bf4=0, f4-1
                         do bf3=0, f3-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             else
                f1=nbf2
                f2=nbf1
                f3=nbf3
                f4=nbf4
                auxIter = 0
                ! (a,b|r,s) -> (b,a|r,s)
                do bf2=0, f2-1
                   do bf1=0, f1-1
                      do bf3=0, f3-1
                         do bf4=0, f4-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             end if
          else
             f1=nbf1
             f2=nbf2
             f3=nbf4
             f4=nbf3
             auxIter = 0
             ! (a,b|r,s) -> (a,b|s,r)
             do bf1=0, f1-1
                do bf2=0, f2-1
                   do bf4=0, f4-1
                      do bf3=0, f3-1
                         newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                         deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                         auxIter = auxIter + 1
                      end do
                   end do
                end do
             end do
          end if
       else
          if (p12) then
             if (p34) then
                f1=nbf4
                f2=nbf3
                f3=nbf2
                f4=nbf1
                auxIter = 0
                ! (a,b|r,s) -> (s,r|b,a)
                do bf4=0, f4-1
                   do bf3=0,f3-1
                      do bf2=0, f2-1
                         do bf1=0, f1-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             else
                f1=nbf4
                f2=nbf3
                f3=nbf1
                f4=nbf2
                auxIter = 0
                ! (a,b|r,s) -> (s,r|a,b)
                do bf3=0, f3-1
                   do bf4=0, f4-1
                      do bf2=0, f2-1
                         do bf1=0, f1-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             end if
          else
             if (p34) then
                f1=nbf3
                f2=nbf4
                f3=nbf2
                f4=nbf1
                auxIter = 0
                ! (a,b|r,s) -> (r,s|b,a)
                do bf4=0, f4-1
                   do bf3=0, f3-1
                      do bf1=0, f1-1
                         do bf2=0, f2-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             else
                f1=nbf3
                f2=nbf4
                f3=nbf1
                f4=nbf2
                auxIter = 0
                ! (a,b|r,s) -> (r,s|a,b)
                do bf3=0, f3-1
                   do bf4=0, f4-1
                      do bf1=0, f1-1
                         do bf2=0, f2-1
                            newPtr = ((bf1*f2 + bf2)*f3 + bf3)*f4 + bf4
                            deriveValue(newPtr+derivIter*size) = auxDerivatives(auxIter+derivIter*size)
                            auxIter = auxIter + 1
                         end do
                      end do
                   end do
                end do
             end if
          end if
       end if
    end do

  end subroutine RepulsionDerivatives_permute

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
