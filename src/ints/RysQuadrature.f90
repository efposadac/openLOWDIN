!!******************************************************************************
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
!! @brief RysQuad interface
!!
!! @author E. F. Posada
!!
!! <b> Creation data : </b> 11-06-2013
!!
!! <b> History change: </b>
!!
!!   - <tt> 11-06-13 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# This module interfaces the rysquad module form Ruben Guerrero.
!!
module RysQuadrature_
  use Exception_
  use MolecularSystem_
  use OverlapIntegrals_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use Matrix_
  use, intrinsic :: iso_c_binding
  implicit none

  public :: &
  RysQuadrature_computeIntraSpecies, &
  RysQuadrature_directIntraSpecies
  
#define contr(n,m) contractions(n)%contractions(m)
  
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
  
  interface
     
      !>
      !!Interface to coulomb_repulsion
!!      function RysQuadrature_CoulRep(&
!!           xa, ya, za, norma, la, alphaa, &
!!           xb, yb, zb, normb, lb, alphab, &
!!           xc, yc, zc, normc, lc, alphac, &
!!           xd, yd, zd, normd, ld, alphad) bind(C, name="Coul_Rep")
!!        use, intrinsic :: iso_c_binding
!!        implicit none
!!      
!!        real(kind=c_double) :: RysQuadrature_CoulRep
!!      
!!        real(kind=c_double), value :: xa, ya, za
!!        real(kind=c_double), value :: xb, yb, zb
!!        real(kind=c_double), value :: xc, yc, zc
!!        real(kind=c_double), value :: xd, yd, zd
!!      
!!        real(kind=c_double), value :: norma
!!        real(kind=c_double), value :: normb
!!        real(kind=c_double), value :: normc
!!        real(kind=c_double), value :: normd
!!      
!!        real(kind=c_double), value :: alphaa
!!        real(kind=c_double), value :: alphab
!!        real(kind=c_double), value :: alphac
!!        real(kind=c_double), value :: alphad
!!      
!!        integer(kind=c_int), value :: la
!!        integer(kind=c_int), value :: lb
!!        integer(kind=c_int), value :: lc
!!        integer(kind=c_int), value :: ld
!!      
!!      end function RysQuadrature_CoulRep

  end interface
  
  !> @brief Integrals Stack
  type(erisStack), private :: eris
  
contains
  
  !>
  !! @brief calculate eris using quadrature of Rys module
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @info Tested! up to "d" angular momentum
  subroutine RysQuadrature_computeIntraSpecies(specieID, job, starting, ending, process)
    implicit none

    character(*), intent(in) :: job
    character(50) :: wfnFile
    integer, intent(in) :: specieID
    integer(8), intent(in) :: starting
    integer(8), intent(in) :: ending
    integer, intent(in) :: process

    character(50) :: arguments(20)
    
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
    integer :: wfnUnit

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
    type(Matrix) :: coefficients
    
    external coul_rep


    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))
    
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
    
    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))
    
    !!Real labels for contractions
    aux = 1
    do i = 1, numberOfContractions
       !!position for cartesian contractions
       labelsOfContractions(i) = aux
       aux = aux + contractions(i)%numCartesianOrbital          
    end do
    
    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2
    
    if(allocated(auxIntegrals)) deallocate(auxIntegrals)
    if(allocated(integralsValue)) deallocate(integralsValue)
    
    allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize))
    allocate(integralsValue(arraySize* arraySize* arraySize* arraySize))
    allocate(integralsPtr(arraySize* arraySize* arraySize* arraySize))    

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    
    counter = 0
    auxCounter = 0    
    control = 0
    integralsPtr = 0.0_8
    auxIntegrals = 0.0_8
    integralsValue = 0.0_8
    
    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       n = a
       do b = a, numberOfContractions
          u = b
          do r = n , numberOfContractions
             do s = u,  numberOfContractions
                                
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
                      
                      close(34)

                      return
                      
                   end if
                
                   !!Calcula el momento angular total
                   sumAngularMoment =  contractions(a)%angularMoment + &
                        contractions(b)%angularMoment + &
                        contractions(r)%angularMoment + &
                        contractions(s)%angularMoment
                   
                   !!Calcula el tamano del arreglo de integrales para la capa (ab|rs)
                   arraySize = contractions(a)%numCartesianOrbital * &
                        contractions(b)%numCartesianOrbital * &
                        contractions(r)%numCartesianOrbital * &
                        contractions(s)%numCartesianOrbital

                   
                   !!************************************
                   !! Calculate iteratively primitives
                   !!

                   !!start :)
                   integralsValue(1:arraySize) = 0.0_8
                   
                   do l = 1, contractions(s)%length
                      do k = 1, contractions(r)%length
                         do j = 1, contractions(b)%length
                            do i = 1, contractions(a)%length                               
                               
                               call coul_rep(&
                                    contractions(a)%origin(1), contractions(a)%origin(2), contractions(a)%origin(3), &
                                    contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&
                                    contractions(b)%origin(1), contractions(b)%origin(2), contractions(b)%origin(3), &
                                    contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&
                                    contractions(r)%origin(1), contractions(r)%origin(2), contractions(r)%origin(3), &
                                    contractions(r)%angularMoment, contractions(r)%orbitalExponents(k), labelsOfContractions(r),&
                                    contractions(s)%origin(1), contractions(s)%origin(2), contractions(s)%origin(3), &
                                    contractions(s)%angularMoment, contractions(s)%orbitalExponents(l), labelsOfContractions(s),&
                                    b, s, integralsPtr(1:arraySize))                               
                               
                               
                               auxIntegrals(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy                                                                    
                               
                               m = 0
                               do ii = 1, contractions(a)%numCartesianOrbital
                                  do jj = 1, contractions(b)%numCartesianOrbital
                                     do kk = 1, contractions(r)%numCartesianOrbital
                                        do ll = 1, contractions(s)%numCartesianOrbital
                                           m = m + 1
                                           auxIntegrals(m) = auxIntegrals(m) &
                                                * contractions(a)%primNormalization(i,ii) &
                                                * contractions(b)%primNormalization(j,jj) &
                                                * contractions(r)%primNormalization(k,kk) &
                                                * contractions(s)%primNormalization(l,ll)
                                        end do
                                     end do
                                  end do
                               end do !! done by cartesian of contractions                               
                               
                               auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
                                    * contractions(a)%contractionCoefficients(i) &
                                    * contractions(b)%contractionCoefficients(j) &
                                    * contractions(r)%contractionCoefficients(k) &
                                    * contractions(s)%contractionCoefficients(l)                                                             
                               
                               integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)
                               
                            end do
                         end do
                      end do
                   end do !!done integral by shell                                   
                   
                   !!normalize by shell
                   m = 0
                   do ii = 1,  contractions(a)%numCartesianOrbital
                      do jj = 1,  contractions(b)%numCartesianOrbital
                         do kk = 1,  contractions(r)%numCartesianOrbital
                            do ll = 1, contractions(s)%numCartesianOrbital
                               m = m + 1
                               integralsValue(m) = integralsValue(m) &
                                    * contractions(a)%contNormalization(ii) &
                                    * contractions(b)%contNormalization(jj) &
                                    * contractions(r)%contNormalization(kk) &
                                    * contractions(s)%contNormalization(ll)                            
                            end do
                         end do
                      end do
                   end do !! done by shell
                   
                   !!write to disk
                   m = 0
                   do i = 1, contractions(a)%numCartesianOrbital
                      do j = 1, contractions(b)%numCartesianOrbital
                         do k = 1, contractions(r)%numCartesianOrbital
                            do l = 1, contractions(s)%numCartesianOrbital
                               
                               m = m + 1
                               
                               !! index not permuted
                               pa=labelsOfContractions(a)+i-1
                               pb=labelsOfContractions(b)+j-1
                               pr=labelsOfContractions(r)+k-1
                               ps=labelsOfContractions(s)+l-1
                               
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
    
    close(34)

    write(6,"(A,I12,A,A)") " Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
         trim(MolecularSystem_instance%species(specieID)%name)

    deallocate (eris%a)
    deallocate (eris%b)
    deallocate (eris%c)
    deallocate (eris%d)
    deallocate (eris%integrals)
    
  end subroutine RysQuadrature_computeIntraSpecies

  !>
  !! @brief calculate eris using quadrature of Rys module
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @info Tested! up to "d" angular momentum
  recursive subroutine RysQuadrature_directIntraSpecies(specieID, job, starting, ending, process, &
                 densityMatrix, twoParticlesMatrix, factor)
    implicit none

    character(*), intent(in) :: job
    character(50) :: wfnFile
    integer, intent(in) :: specieID
    integer(8), intent(in) :: starting
    integer(8), intent(in) :: ending
    integer, intent(in) :: process

    type(matrix), intent(inout):: twoParticlesMatrix
    type(matrix), intent(in) :: densityMatrix
    real(8), intent(in) :: factor
    real(8) :: coulomb,exchange

    character(50) :: arguments(20)
    
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
    integer :: wfnUnit

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
    type(Matrix) :: coefficients
    
!!    external coul_rep

!!    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))
    
    write(fileNumber,*) process
    fileNumber = trim(adjustl(fileNumber))
    
    !! open file for integrals
!!    open(UNIT=34,FILE=trim(fileNumber)//trim(MolecularSystem_instance%species(specieID)%name)//".ints", &
!!         STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='Unformatted')
    
    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID)
    
    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize
    
    ssize = (numberOfContractions * (numberOfContractions + 1))/2
    ssize = (ssize * (ssize + 1))/2
    
    !! Get contractions labels for integrals index
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)
    allocate(labelsOfContractions(numberOfContractions))
    
    !!Real labels for contractions
    aux = 1
    do i = 1, numberOfContractions
       !!position for cartesian contractions
       labelsOfContractions(i) = aux
       aux = aux + contractions(i)%numCartesianOrbital          
    end do
    
    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2
    
    if(allocated(auxIntegrals)) deallocate(auxIntegrals)
    if(allocated(integralsValue)) deallocate(integralsValue)
    
    allocate(auxIntegrals(arraySize* arraySize* arraySize * arraySize))
    allocate(integralsValue(arraySize* arraySize* arraySize* arraySize))
    allocate(integralsPtr(arraySize* arraySize* arraySize* arraySize))    

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    
    counter = 0
    auxCounter = 0    
    control = 0
    integralsPtr = 0.0_8
    auxIntegrals = 0.0_8
    integralsValue = 0.0_8
    
    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       n = a
       do b = a, numberOfContractions
          u = b
          do r = n , numberOfContractions
             do s = u,  numberOfContractions
                                
                control = control + 1                 
                if( control >= starting ) then                    
                   
                   if (control > ending) then
!!                      
!!                      counter = counter + 1
!!                      eris%a(counter) = -1
!!                      eris%b(counter) = -1
!!                      eris%c(counter) = -1
!!                      eris%d(counter) = -1
!!                      eris%integrals(counter) = 0.0_8
!!                      
!!                      write(34) &
!!                           eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                           eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                           eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                           eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                           eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
!!                   
!!                      write(6,"(A,I12,A,A)") "**Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
!!                           trim(MolecularSystem_instance%species(specieID)%name)
!!                      
!!                      close(34)
!!
                      return
!!                      
                   end if
                
                   !!Calcula el momento angular total
                   sumAngularMoment =  contractions(a)%angularMoment + &
                        contractions(b)%angularMoment + &
                        contractions(r)%angularMoment + &
                        contractions(s)%angularMoment
                   
                   !!Calcula el tamano del arreglo de integrales para la capa (ab|rs)
                   arraySize = contractions(a)%numCartesianOrbital * &
                        contractions(b)%numCartesianOrbital * &
                        contractions(r)%numCartesianOrbital * &
                        contractions(s)%numCartesianOrbital

                   
                   !!************************************
                   !! Calculate iteratively primitives
                   !!

                   !!start :)
                   integralsValue(1:arraySize) = 0.0_8
                   
                   do l = 1, contractions(s)%length
                      do k = 1, contractions(r)%length
                         do j = 1, contractions(b)%length
                            do i = 1, contractions(a)%length                               
                               
                               call coul_rep(&
                                    contractions(a)%origin(1), contractions(a)%origin(2), contractions(a)%origin(3), &
                                    contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&
                                    contractions(b)%origin(1), contractions(b)%origin(2), contractions(b)%origin(3), &
                                    contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&
                                    contractions(r)%origin(1), contractions(r)%origin(2), contractions(r)%origin(3), &
                                    contractions(r)%angularMoment, contractions(r)%orbitalExponents(k), labelsOfContractions(r),&
                                    contractions(s)%origin(1), contractions(s)%origin(2), contractions(s)%origin(3), &
                                    contractions(s)%angularMoment, contractions(s)%orbitalExponents(l), labelsOfContractions(s),&
                                    b, s, integralsPtr(1:arraySize))                               
                               
                               
                               auxIntegrals(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy                                                                   
!                               auxIntegrals(1:arraySize) = 1
                               
                               m = 0
                               do ii = 1, contractions(a)%numCartesianOrbital
                                  do jj = 1, contractions(b)%numCartesianOrbital
                                     do kk = 1, contractions(r)%numCartesianOrbital
                                        do ll = 1, contractions(s)%numCartesianOrbital
                                           m = m + 1
                                           auxIntegrals(m) = auxIntegrals(m) &
                                                * contractions(a)%primNormalization(i,ii) &
                                                * contractions(b)%primNormalization(j,jj) &
                                                * contractions(r)%primNormalization(k,kk) &
                                                * contractions(s)%primNormalization(l,ll)
                                        end do
                                     end do
                                  end do
                               end do !! done by cartesian of contractions                               
                               
                               auxIntegrals(1:arraySize) = auxIntegrals(1:arraySize) &
                                    * contractions(a)%contractionCoefficients(i) &
                                    * contractions(b)%contractionCoefficients(j) &
                                    * contractions(r)%contractionCoefficients(k) &
                                    * contractions(s)%contractionCoefficients(l)                                                             
                               
                               integralsValue(1:arraySize) = integralsValue(1:arraySize) + auxIntegrals(1:arraySize)
                               
                            end do
                         end do
                      end do
                   end do !!done integral by shell                                   
                   
                   !!normalize by shell
                   m = 0
                   do ii = 1,  contractions(a)%numCartesianOrbital
                      do jj = 1,  contractions(b)%numCartesianOrbital
                         do kk = 1,  contractions(r)%numCartesianOrbital
                            do ll = 1, contractions(s)%numCartesianOrbital
                               m = m + 1
                               integralsValue(m) = integralsValue(m) &
                                    * contractions(a)%contNormalization(ii) &
                                    * contractions(b)%contNormalization(jj) &
                                    * contractions(r)%contNormalization(kk) &
                                    * contractions(s)%contNormalization(ll)                            
                            end do
                         end do
                      end do
                   end do !! done by shell
                   
                   !!write to disk
                   m = 0
                   do i = 1, contractions(a)%numCartesianOrbital
                      do j = 1, contractions(b)%numCartesianOrbital
                         do k = 1, contractions(r)%numCartesianOrbital
                            do l = 1, contractions(s)%numCartesianOrbital
                               
                               m = m + 1
                               
                               !! index not permuted
                               pa=labelsOfContractions(a)+i-1
                               pb=labelsOfContractions(b)+j-1
                               pr=labelsOfContractions(r)+k-1
                               ps=labelsOfContractions(s)+l-1
                               
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
                                     
!                                    if(abs(integralsValue(m)) > 1.0D-10) then
                                        
                                        counter = counter + 1
                                        auxCounter = auxCounter + 1
!                                        integralsValue(m) = 1
                                        
!!                                        eris%a(counter) = pa
!!                                        eris%b(counter) = pb
!!                                        eris%c(counter) = pr
!!                                        eris%d(counter) = ps
!!                                        eris%integrals(counter) = integralsValue(m)                                                                                
                                        
!!                                     end if
                                     
!!                                     if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then
!!                                        
!!                                        write(34) &
!!                                             eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                                             eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                                             eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                                             eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!                                             eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
!!                                        
!!                                        counter = 0
!!                                        
!!                                     end if

                                  coulomb = densityMatrix%values(pr,ps)*integralsValue(m)

                                  !!*****************************************************************************
                                  !! 	Adiciona aportes debidos al operador de coulomb
                                  !!*****
                                  if( pa == pr .and. pb == ps ) then

                                     twoParticlesMatrix%values(pa,pb) = &
                                          twoParticlesMatrix%values(pa,pb) + coulomb

                                     if( pr /= ps ) then

                                        twoParticlesMatrix%values(pa,pb) = &
                                             twoParticlesMatrix%values(pa,pb) + coulomb

                                     end if

                                  else

                                     twoParticlesMatrix%values(pa,pb) = &
                                          twoParticlesMatrix%values(pa,pb) + coulomb

                                     if( pr /= ps ) then

                                        twoParticlesMatrix%values(pa,pb) = &
                                             twoParticlesMatrix%values(pa,pb) + coulomb

                                     end if

                                     coulomb = densityMatrix%values(pa,pb)*integralsValue(m)

                                     twoParticlesMatrix%values(pr,ps) = &
                                          twoParticlesMatrix%values(pr,ps) + coulomb

                                     if ( pa /= pb ) then

                                        twoParticlesMatrix%values( pr, ps ) = &
                                             twoParticlesMatrix%values( pr, ps ) + coulomb

                                     end if

                                  end if

                                  !!
                                  !!*****************************************************************************

                                  !!*****************************************************************************
                                  !! 	Adicionaa aportess debidoss al operadorr de intercambio
                                  !!*****
                                  if( pr /= ps ) then

                                     exchange =densityMatrix%values(pb,ps)*integralsValue(m)* factor

                                     twoParticlesMatrix%values( pa, pr ) = &
                                          twoParticlesMatrix%values( pa, pr ) + exchange

                                     if( pa == pr .and. pb /= ps ) then

                                        twoParticlesMatrix%values( pa, pr ) = &
                                             twoParticlesMatrix%values( pa, pr ) + exchange

                                     end if

                                  end if

                                  if ( pa /= pb ) then

                                     exchange = densityMatrix%values(pa,pr)*integralsValue(m) * factor

                                     if( pb > ps ) then

                                        twoParticlesMatrix%values( ps, pb ) = &
                                             twoParticlesMatrix%values( ps, pb) + exchange

                                     else

                                        twoParticlesMatrix%values( pb, ps ) = &
                                             twoParticlesMatrix%values( pb, ps ) + exchange

                                        if( pb==ps .and. pa /= pr ) then

                                           twoParticlesMatrix%values( pb, ps ) = &
                                                twoParticlesMatrix%values( pb, ps ) + exchange

                                        end if

                                     end if

                                     if ( pr /= ps ) then

                                        exchange=densityMatrix%values(pa,ps)*integralsValue(m) * factor

                                        if( pb <= pr ) then

                                           twoParticlesMatrix%values( pb, pr ) = &
                                                twoParticlesMatrix%values( pb, pr ) + exchange

                                           if( pb == pr ) then

                                              twoParticlesMatrix%values( pb, pr ) = &
                                                   twoParticlesMatrix%values( pb, pr ) + exchange

                                           end if

                                        else

                                           twoParticlesMatrix%values( pr, pb ) = &
                                                twoParticlesMatrix%values( pr, pb) + exchange

                                           if( pa == pr .and. ps == pb ) goto 30

                                        end if

                                     end if

                                  end if


                                  exchange = densityMatrix%values(pb,pr)*integralsValue(m) * factor

                                  twoParticlesMatrix%values( pa, ps ) = &
                                       twoParticlesMatrix%values( pa, ps ) + exchange

30                                continue

                                  !!
                                  !!*****************************************************************************
                                                     
!                                  end if
                                  end if
                               end if

                            end do
                         end do
                      end do
                   end do !! Primitives loop

                   
                end if !starting
                
             end do
             u=r+1
          end do
       end do
    end do !! shells loop

    counter = counter + 1
!!    eris%a(counter) = -1
!!    eris%b(counter) = -1
!!    eris%c(counter) = -1
!!    eris%d(counter) = -1
!!    eris%integrals(counter) = 0.0_8
!!       
!!    write(34) &
!!         eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
!!         eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
!!    
!!    close(34)
!!
!!    write(6,"(A,I12,A,A)") " Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
!!         trim(MolecularSystem_instance%species(specieID)%name)
!!
!!    deallocate (eris%a)
!!    deallocate (eris%b)
!!    deallocate (eris%c)
!!    deallocate (eris%d)
!!    deallocate (eris%integrals)
!!    
!!    write(6,"(A,I12,A,A)") " Stored ", auxCounter, " non-zero repulsion integrals of species: ", &
!!         trim(MolecularSystem_instance%species(specieID)%name)
        
    
  end subroutine RysQuadrature_directIntraSpecies

  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine RysQuadrature_exception( typeMessage, description, debugDescription)
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

  end subroutine RysQuadrature_exception

end module RysQuadrature_
