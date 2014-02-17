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
!! @authors E. F. Posada,
!!
!! <b> Creation data : </b> 11-06-2013
!!
!! <b> History change: </b>
!!
!!   - <tt> 11-06-13 </tt>:  Edwin Posada ( efposadac@unal.edu.co )
!!        -# This module interfaces the rysquad module form Ruben Guerrero.
!!   -<tt>  </tt0
!!
module RysQuadrature_
  use Exception_
  use LibintInterfaceTypes_
  use MolecularSystem_
  use OverlapIntegrals_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use Matrix_
  use AttractionIntegrals_
  use, intrinsic :: iso_c_binding
  implicit none
  
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
     
     ! !>
     ! !!Interface to coulomb_repulsion
     ! function RysQuadrature_CoulRep(&
     !      xa, ya, za, norma, la, alphaa, &
     !      xb, yb, zb, normb, lb, alphab, &
     !      xc, yc, zc, normc, lc, alphac, &
     !      xd, yd, zd, normd, ld, alphad) bind(C, name="Coul_Rep")
     !   use, intrinsic :: iso_c_binding
     !   implicit none
       
     !   real(kind=c_double) :: RysQuadrature_CoulRep
       
     !   real(kind=c_double), value :: xa, ya, za
     !   real(kind=c_double), value :: xb, yb, zb
     !   real(kind=c_double), value :: xc, yc, zc
     !   real(kind=c_double), value :: xd, yd, zd
       
     !   real(kind=c_double), value :: norma
     !   real(kind=c_double), value :: normb
     !   real(kind=c_double), value :: normc
     !   real(kind=c_double), value :: normd
       
     !   real(kind=c_double), value :: alphaa
     !   real(kind=c_double), value :: alphab
     !   real(kind=c_double), value :: alphac
     !   real(kind=c_double), value :: alphad
       
     !   integer(kind=c_int), value :: la
     !   integer(kind=c_int), value :: lb
     !   integer(kind=c_int), value :: lc
     !   integer(kind=c_int), value :: ld
       
     ! end function RysQuadrature_CoulRep

  end interface
  
  !> @brief Integrals Stack
  type(erisStack), private :: eris
  
contains
  
  !>
  !! @brief calculate eris using quadrature of Rys module
  !! @authors E. F. Posada, 2013
  !! @version 1.0
  !! @info Tested! up to "d" angular momentum
  subroutine RysQuadrature_computeIntraSpecies(specieID, job, starting, ending, process)
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
    integer*2 :: pa, pb, pr, ps !< labels index
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
    external coul_rep_3d
    external coul_rep_2d
    external coul_rep_1d

    external ecg_coul_rep_3d

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

                   select case(CONTROL_instance%DIMENSIONALITY)
                   case(3)
                      do l = 1, contractions(s)%length
                         do k = 1, contractions(r)%length
                            do j = 1, contractions(b)%length
                               do i = 1, contractions(a)%length                               
                               
                                     ! no calcular "por ahora"
                                     ! call ecg_coul_rep_3d(&
                                     !   contractions(a)%origin(1), contractions(a)%origin(2), contractions(a)%origin(3), &
                                     !   contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&
                                     !   contractions(b)%origin(1), contractions(b)%origin(2), contractions(b)%origin(3), &
                                     !   contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&
                                     !   contractions(r)%origin(1), contractions(r)%origin(2), contractions(r)%origin(3), &
                                     !   contractions(r)%angularMoment, contractions(r)%orbitalExponents(k), labelsOfContractions(r),&
                                     !   contractions(s)%origin(1), contractions(s)%origin(2), contractions(s)%origin(3), &
                                     !   contractions(s)%angularMoment, contractions(s)%orbitalExponents(l), labelsOfContractions(s),&
                                     !   b, s, integralsPtr(1:arraySize))

                                     call coul_rep_3d(&
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
                   case(2)
                      do l = 1, contractions(s)%length
                         do k = 1, contractions(r)%length
                            do j = 1, contractions(b)%length
                               do i = 1, contractions(a)%length                               

                                  call coul_rep_2d(&
                                       contractions(a)%origin(1), contractions(a)%origin(2),  &
                                       contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&
                                       contractions(b)%origin(1), contractions(b)%origin(2), &
                                       contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&
                                       contractions(r)%origin(1), contractions(r)%origin(2), &
                                       contractions(r)%angularMoment, contractions(r)%orbitalExponents(k), labelsOfContractions(r),&
                                       contractions(s)%origin(1), contractions(s)%origin(2),  &
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
                   case(1)
                      do l = 1, contractions(s)%length
                         do k = 1, contractions(r)%length
                            do j = 1, contractions(b)%length
                               do i = 1, contractions(a)%length                               

                                  call coul_rep_1d(&
                                       contractions(a)%origin(1), &
                                       contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&
                                       contractions(b)%origin(1), &
                                       contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&
                                       contractions(r)%origin(1), &
                                       contractions(r)%angularMoment, contractions(r)%orbitalExponents(k), labelsOfContractions(r),&
                                       contractions(s)%origin(1), &
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
                      
                   
                   end select
                   
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

  
  subroutine RysQuadrature_attractionECG (rrx,rry,rrz,rram,rrExp, &
       rx,ry,rz,aam,aExp,sam,sExp, &
       sx,sy,sz,bam,bExp,ram,rExp, &
       ssx,ssy,ssz,ssam,ssExp,integral) 

    implicit none

    real(8), intent(in) :: rrx,rry,rrz
    real(8), intent(in) :: rx,ry,rz
    real(8), intent(in) :: sx,sy,sz
    real(8), intent(in) :: ssx,ssy,ssz
    integer, intent(in) :: rram,sam,ram,ssam
    integer, intent(in) :: aam,bam
    real(8), intent(in) :: rrExp,sExp,rExp,ssExp
    real(8), intent(in) :: aExp,bExp
    real(8), intent(out) :: integral(:)

    integer :: counter
    real(8) :: preFactor,A,B,expFactorX, expFactorY, expFactorZ, erfFactorX, erfFactorY, erfFactorZ,erfFactor

    A = aExp + bExp
    B = rrExp + sExp + rExp + ssExp

    expFactorX = sx**2 *rrExp*sExp + sx**2*rExp*sExp + ssx**2 *rrExp*ssExp + ssx**2*rExp*ssExp &
         + sx**2 * sExp*ssExp - 2*sx*ssx*sExp*ssExp + ssx**2* sExp*ssExp &
         + rx**2 * rExp * ( rrExp + sExp + ssExp ) + rrx**2 * rrExp * (rExp + sExp + ssExp ) &
         - 2*rx* rExp * ( sx*sExp + ssx*ssExp ) - 2*rrx* rrExp * ( rx*rExp + sx*sExp + ssx*ssExp) 

    expFactorY = sy**2 *rrExp*sExp + sy**2*rExp*sExp + ssy**2 *rrExp*ssExp + ssy**2*rExp*ssExp &
         + sy**2 * sExp*ssExp - 2*sy*ssy*sExp*ssExp + ssy**2* sExp*ssExp &
         + ry**2 * rExp * ( rrExp + sExp + ssExp ) + rry**2 * rrExp * (rExp + sExp + ssExp ) &
         - 2*ry* rExp * ( sy*sExp + ssy*ssExp ) - 2*rry* rrExp * ( ry*rExp + sy*sExp + ssy*ssExp) 

    expFactorY = sz**2 *rrExp*sExp + sz**2*rExp*sExp + ssz**2 *rrExp*ssExp + ssz**2*rExp*ssExp &
         + sz**2 * sExp*ssExp - 2*sz*ssz*sExp*ssExp + ssz**2* sExp*ssExp &
         + rz**2 * rExp * ( rrExp + sExp + ssExp ) + rrz**2 * rrExp * (rExp + sExp + ssExp ) &
         - 2*rz* rExp * ( sz*sExp + ssz*ssExp ) - 2*rrz* rrExp * ( rz*rExp + sz*sExp + ssz*ssExp) 

    preFactor = exp( -(expFactorX + expFactorY + expFactorZ )/(B) )
    preFactor = preFactor * 4.0* Math_PI**(3.0/2.0) / (A**(3.0/2.0) * B**(3.0/2.0) ) 

    erfFactor = ( A*B  ) / ( A+B )

    integral = Math_PI * preFactor * ( 1 + erf (erfFactor) )

    !print *, "integral", integral

  end subroutine RysQuadrature_attractionECG
  
  subroutine RysQuadrature_attractionECG2N2E ( &
       rx,ry,rz,aam,aExp,sam,sExp, &
       sx,sy,sz,bam,bExp,ram,rExp, integral)

    implicit none

    real(8), intent(in) :: rx,ry,rz
    real(8), intent(in) :: sx,sy,sz
    integer, intent(in) :: sam,ram
    integer, intent(in) :: aam,bam
    real(8), intent(in) :: sExp,rExp
    real(8), intent(in) :: aExp,bExp
    real(8), intent(out) :: integral(:)

    integer :: counter
    real(8) :: preFactor,A,B,expFactorX, expFactorY, expFactorZ, erfFactorX, erfFactorY, erfFactorZ

    A = aExp + bExp
    B = sExp + rExp

    expFactorX = (sx - rx)**2*rExp*sExp / ( rExp + sExp ) 
    expFactorY = (sy - ry)**2*rExp*sExp / ( rExp + sExp ) 
    expFactorZ = (sz - rz)**2*rExp*sExp / ( rExp + sExp ) 


    preFactor = exp( -(expFactorX + expFactorY + expFactorZ ) )
    preFactor = preFactor * Math_PI**3 / (2.0**3 * A**(3.0/2.0) * B**(3.0/2.0) ) 

    erfFactorX = ( rx*rExp + sx*sExp ) / ( B**(1.0/2.0) )
    erfFactorY = ( ry*rExp + sy*sExp ) / ( B**(1.0/2.0) )
    erfFactorZ = ( rz*rExp + sz*sExp ) / ( B**(1.0/2.0) )

    integral = preFactor * ( 1 + erf (erfFactorX) )*( 1 + erf (erfFactorY) )*( 1 + erf (erfFactorZ) )

  end subroutine RysQuadrature_attractionECG2N2E
  
  !>
  !! @brief calculate eris using quadrature of Rys module
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @info Tested! up to "d" angular momentum
  subroutine  RysQuadrature_computeInterSpecies(specieID, otherSpecieID, job )

    implicit none

    character(*), intent(in) :: job
    integer, intent(in) :: specieID
    integer, intent(in) :: otherSpecieID

    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: totalNumberOfContractions
    integer :: otherTotalNumberOfContractions
    integer :: numberOfPrimitives
    integer :: otherNumberOfPrimitives
    integer :: maxAngularMoment
    integer :: arraySize !< number of cartesian orbitals for maxAngularMoment
    integer :: u,m !< auxiliary itetators for this specie
    
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT)
    integer :: a, b, r, s !< not permuted iterators (original) for this specie
    
    integer*2 :: pa, pb, pr, ps !< labels index
    
    integer :: ii, jj, kk, ll, nn, oo !< cartesian iterators for primitives and contractions
    integer :: aux, order !<auxiliary index
    
    integer :: counter, auxCounter

    integer :: i, j, k, l, n, o !< contraction length iterators for this specie
    
    integer, allocatable :: labelsOfContractions(:) !< cartesian position of contractions in all basis set
    integer, allocatable :: otherLabelsOfContractions(:) !! cartesian position of contractions in all basis set    

    real(8), allocatable :: integralsPtr(:) !< pointer to C array of integrals

    real(8), allocatable :: auxIntegralsECG(:) !<array with permuted integrals aux!
    real(8), allocatable :: integralsValueECG (:) !<array with permuted integrals
    
    real(8) :: auxIntegralValue

    character(50) :: fileNumber
    integer(8) :: ssize

    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for one of the species
    type(ContractedGaussian), allocatable :: otherContractions(:) !< Basis set for other specie
    
    external coul_rep_3d
    external coul_rep_2d
    external coul_rep_1d

    external ecg_coul_rep_3d

    external ecg_coul_attract_3d

    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))

    !! open file for integrals
    open(UNIT=34,FILE=trim(MolecularSystem_instance%species(specieID)%name)//"."//trim(MolecularSystem_instance%species(otherSpecieID)%name)//".ints", &
         STATUS='REPLACE', ACCESS='SEQUENTIAL', FORM='Unformatted')

    !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    call MolecularSystem_getBasisSet(otherSpecieID, otherContractions)

    maxAngularMoment = max(MolecularSystem_getMaxAngularMoment(specieID), MolecularSystem_getMaxAngularMoment(otherSpecieID))
    numberOfPrimitives = MolecularSystem_getTotalNumberOfContractions(specieID) + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

    !! Get number of shells and cartesian contractions                                                  
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(specieID)%basisSetSize
    
    !! Get number of shells and cartesian contractions (other specie)                                                                                
    otherNumberOfContractions = size(otherContractions)
    otherTotalNumberOfContractions = MolecularSystem_instance%species(otherSpecieID)%basisSetSize

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
    
    !! Get contractions labels for integrals index (other specie)                                                                                     
    if (allocated(otherLabelsOfContractions)) deallocate(otherLabelsOfContractions)
    allocate(otherLabelsOfContractions(otherNumberOfContractions))
    
    !!Real labels for contractions (other specie)                                                                                                     
    aux = 1
    do i = 1, otherNumberOfContractions
       !!position for cartesian contractions
       otherLabelsOfContractions(i) = aux
       aux = aux + otherContractions(i)%numCartesianOrbital
    end do

    !! allocating space for integrals just one time (do not put it inside do loop!!!)
    arraySize = ((maxAngularMoment + 1)*(maxAngularMoment + 2))/2
        
    if(allocated(auxIntegralsECG)) deallocate(auxIntegralsECG)    
    if(allocated(integralsValueECG)) deallocate(integralsValueECG)
    if(allocated(integralsPtr)) deallocate(integralsPtr)

    allocate(auxIntegralsECG(arraySize* arraySize* arraySize * arraySize))
    allocate(integralsValueECG(arraySize* arraySize* arraySize* arraySize))
    allocate(integralsPtr(arraySize* arraySize* arraySize* arraySize))

    counter = 0
    auxCounter = 0    
    integralsValueECG = 0.0_8
    
    !!Start Calculating integrals for each shell
    do a = 1, numberOfContractions
       do b = a, numberOfContractions
          do r = 1 , othernumberofcontractions
             do s = r,  othernumberofcontractions
                
                !!Calcula el tamano del arreglo de integrales para la capa (ab|rs)
                arraySize = contractions(a)%numCartesianOrbital * &
                     contractions(b)%numCartesianOrbital * &
                     otherContractions(r)%numCartesianOrbital  * &
                     otherContractions(s)%numCartesianOrbital
                
                !!start :)
                integralsValueECG(1:arraySize) = 0.0_8
                
                !! not-permuted loop
                do l = 1, otherContractions(s)%length
                   do k = 1, otherContractions(r)%length
                      do j = 1, contractions(b)%length
                         do i = 1, contractions(a)%length
                            
                            auxIntegralsECG(1:arraySize) = 0.0_8
                            integralsPtr(1:arraySize) = 0.0_8
                            
                            call coul_rep_3d(&
                                 contractions(a)%origin(1), contractions(a)%origin(2), contractions(a)%origin(3), &
                                 contractions(a)%angularMoment, contractions(a)%orbitalExponents(i), labelsOfContractions(a),&                                 
                                 contractions(b)%origin(1), contractions(b)%origin(2), contractions(b)%origin(3), &
                                 contractions(b)%angularMoment, contractions(b)%orbitalExponents(j), labelsOfContractions(b),&                                 
                                 othercontractions(r)%origin(1), othercontractions(r)%origin(2), othercontractions(r)%origin(3), &
                                 othercontractions(r)%angularMoment, othercontractions(r)%orbitalExponents(k), otherlabelsOfContractions(r),&                                 
                                 othercontractions(s)%origin(1), othercontractions(s)%origin(2), othercontractions(s)%origin(3), &
                                 othercontractions(s)%angularMoment, othercontractions(s)%orbitalExponents(l), otherlabelsOfContractions(s),&                                 
                                 b, s, integralsPtr(1:arraySize))           
                            
                            auxIntegralsECG(1:arraySize) = integralsPtr(1:arraySize) !!it is to slow with pointer...! so.. copy
                            
                            !write(*, "(4I3, F12.8)") l, k, j, i, auxIntegralsECG(1:arraySize)
                            
                            m = 0
                            do ii = 1, contractions(a)%numCartesianOrbital
                               do jj = 1, contractions(b)%numCartesianOrbital
                                  do kk = 1, otherContractions(r)%numCartesianOrbital
                                     do ll = 1, otherContractions(s)%numCartesianOrbital
                                        m = m + 1
                                        
                                        auxIntegralsECG(m) = auxIntegralsECG(m) &
                                             * contractions(a)%primNormalization(i,ii) &
                                             * othercontractions(r)%primNormalization(k,kk) &
                                             * Contractions(b)%primNormalization(j,jj) &
                                             * otherContractions(s)%primNormalization(l,ll) 
                                        
                                     end do
                                  end do
                               end do
                            end do !! done by primitives
                            
                            auxIntegralsECG(1:arraySize) = auxIntegralsECG(1:arraySize) &
                                 * contractions(a)%contractionCoefficients(i) &
                                 * othercontractions(r)%contractionCoefficients(k) &
                                 * Contractions(b)%contractionCoefficients(j) &
                                 * otherContractions(s)%contractionCoefficients(l)     
                            
                            
                            integralsValueECG(1:arraySize) = integralsValueECG(1:arraySize) + auxIntegralsECG(1:arraySize)
                            
                         end do
                      end do
                   end do
                end do !!done Integral primitives
                
                m = 0
                do ii = 1, contractions(a)%numCartesianOrbital
                   do jj = 1, contractions(b)%numCartesianOrbital
                      do kk = 1, otherContractions(r)%numCartesianOrbital
                         do ll = 1, otherContractions(s)%numCartesianOrbital
                            m = m + 1
                            
                            integralsValueECG(m) = integralsValueECG(m) &
                                 * contractions(a)%contNormalization(ii) &
                                 * othercontractions(r)%contNormalization(kk) &
                                 * Contractions(b)%contNormalization(jj) &
                                 * otherContractions(s)%contNormalization(ll) 
                            
                         end do
                      end do
                   end do
                end do !! done by contractions
                
                !!write to disk
                m = 0
                do i = 1, contractions(a)%numCartesianOrbital
                   do j = 1, contractions(b)%numCartesianOrbital
                      do k = 1, otherContractions(r)%numCartesianOrbital
                         do l = 1, otherContractions(s)%numCartesianOrbital
                            
                            m = m + 1
                            auxIntegralValue = integralsValueECG(m)
                            
                            !! index not permuted
                            pa=labelsOfContractions(a)+i-1
                            pb=labelsOfContractions(b)+j-1
                            pr=otherLabelsOfContractions(r)+k-1
                            ps=otherLabelsOfContractions(s)+l-1
                            
                            !write(*, "(4I3,F12.8)") pa, pb, pr, ps, auxIntegralValue
                            
                            if( pa <= pb .and. pr <= ps ) then

                               if(abs(auxIntegralValue) > 1.0D-10) then

                                  counter = counter + 1
                                  auxCounter = auxCounter + 1

                                  eris%a(counter) = pa
                                  eris%b(counter) = pb
                                  eris%c(counter) = pr
                                  eris%d(counter) = ps
                                  eris%integrals(counter) = auxIntegralValue

                               end if


                               if( counter == CONTROL_instance%INTEGRAL_STACK_SIZE ) then

                                  write(34) eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                                       eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                                  counter = 0
                                  
                               end if
                            end if

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

    write(34) &
         eris%a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

    write(*,"(A,I12,A,A)") " Stored ", auxCounter, &
         " non-zero ECG Coulomb integrals between species: ", &
         trim(MolecularSystem_instance%species(specieID)%name)//" / "//&
         trim(MolecularSystem_instance%species(otherSpecieID)%name)

    close(34)

    deallocate (eris%a)
    deallocate (eris%b)
    deallocate (eris%c)
    deallocate (eris%d)
    deallocate (eris%integrals)        

  end subroutine RysQuadrature_computeInterSpecies
  
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
