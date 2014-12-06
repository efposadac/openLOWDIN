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
!! @brief Cuda Interface for Computation of Electron Repulsion Integrals
!!        This module allows to make calculations of ERIs on GPUs
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-10-08
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module CudintInterface_
  use Exception_
  use MolecularSystem_
  use MatrixInteger_
  use OverlapIntegrals_
  use ContractedGaussian_
  use Math_
  use CONTROL_
  use, intrinsic :: iso_c_binding
  implicit none

  !> @brief the integrals are saved for big records (that reduces the I/O time)
  type, public :: erisStack
     integer*2, allocatable :: a(:)
     integer*2, allocatable :: b(:)
     integer*2, allocatable :: c(:)
     integer*2, allocatable :: d(:)
     real(8), allocatable :: integrals(:)
  end type erisStack
 
  interface

     subroutine cuda_int_intraspecies (&
          numberOfContractions, &
          totalContIntegrals, &
          totalPrimitives, &
          maxNumCartesianOrbital, &
          primNormalizationSize, &
          contractionId, &
          contractionLength, &
          contractionAngularMoment, &
          contractionNumCartesianOrbital, &
          contractionOwner, &
          contractionOrigin, &
          contractionOrbitalExponents, &
          contractionCoefficients, &
          contractionContNormalization, &
          contractionPrimNormalization, &
          contractionIntegrals, &
          contractionIndices, &
          primitiveIndices, &
          numberOfPPUC, &
          labelsOfContractions) bind (C, name = "cuda_int_intraspecies_")
       use, intrinsic :: iso_c_binding
       implicit none
       integer (c_int) :: numberOfContractions
       integer (c_int) :: totalContIntegrals
       integer (c_int) :: totalPrimitives
       integer (c_int) :: maxNumCartesianOrbital
       integer (c_int) :: primNormalizationSize
       integer (c_int) :: contractionId(*)
       integer (c_int) :: contractionLength(*)
       integer (c_int) :: contractionAngularMoment(*)
       integer (c_int) :: contractionNumCartesianOrbital(*)
       integer (c_int) :: contractionOwner(*)
       real (c_double) :: contractionOrigin(*)
       real (c_double) :: contractionOrbitalExponents(*)
       real (c_double) :: contractionCoefficients(*)
       real (c_double) :: contractionContNormalization(numberOfContractions,*)
       real (c_double) :: contractionPrimNormalization(*)
       real (c_double) :: contractionIntegrals(*)
       integer (c_int) :: contractionIndices(*)
       integer (c_int) :: primitiveIndices(*)
       integer (c_int) :: numberOfPPUC(*)
       integer (c_int) :: labelsOfContractions(*)
     end subroutine cuda_int_intraspecies
    
 
  end interface

  !> @brief Integrals Stack
  type(erisStack), private :: eris
  
contains
  

  subroutine CudintInterface_computeIntraSpecies(specieID)
    implicit none
    
    integer, intent(in) :: specieID
    
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    integer :: maxAngularMoment
    integer :: sumAngularMoment
    integer :: counter, control, aux, order
    integer :: auxCounter, auxCounter2, auxCounter3


    integer :: n,u,m !< auxiliary itetators
    integer :: aa, bb, rr, ss !< permuted iterators (LIBINT style)
    integer :: a, b, r, s !< not permuted iterators (original)
    integer :: pa, pb, pr, ps !< labels index
    integer :: apa, apb, apr, aps
    integer,target :: ii, jj, kk, ll !< cartesian iterators for primitives and contractions
    integer,target :: i, j, k, l !< contraction length iterators
    integer,pointer :: pi, pj, pk, pl !< pointer to contraction length iterators
    integer,pointer :: poi, poj, pok, pol !< pointer to contraction length iterators
    type(ContractedGaussian), allocatable :: contractions(:) !< Basis set for specie

    integer, allocatable :: contractionId(:)
    integer, allocatable :: labelsOfContractions(:)
    integer, allocatable :: contractionLength(:)
    integer, allocatable :: contractionAngularMoment(:)
    integer, allocatable :: contractionNumCartesianOrbital(:)
    integer, allocatable :: contractionOwner(:)
    real(8), allocatable :: contractionOrigin(:)
    real(8), allocatable :: contractionOrbitalExponents(:)
    real(8), allocatable :: contractionCoefficients(:)
    real(8), allocatable :: contractionContNormalization(:,:)
    real(8), allocatable :: contractionPrimNormalization(:)
    integer :: maxNumCartesianOrbital
    integer :: primNormalizationSize
    integer :: unicIntegrals

    integer :: totalContIntegrals, totalPrimitives
    type(MatrixInteger) :: auxContIndices
    type(MatrixInteger) :: auxPrimIndices
    type(MatrixInteger) :: defContIndices
    type(MatrixInteger) :: auxNumberOfPPUC

    real(8), allocatable :: contractionIntegrals(:)
    integer, allocatable :: contractionIndices(:), primitiveIndices(:), numberOfPPUC(:)
    character(50) :: fileNumber

    allocate (eris%a(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%b(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%c(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%d(CONTROL_instance%INTEGRAL_STACK_SIZE), &
         eris%integrals(CONTROL_instance%INTEGRAL_STACK_SIZE))

    write(fileNumber,*) 1
    fileNumber = trim(adjustl(fileNumber))

   open(UNIT=34,FILE=trim(fileNumber)//trim(MolecularSystem_instance%species(specieID)%name)//".ints", &
        STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='Unformatted')

   !! Get basisSet
    call MolecularSystem_getBasisSet(specieID, contractions)
    
    maxAngularMoment = MolecularSystem_getMaxAngularMoment(specieID)

    !! Get number of shells and cartesian contractions
    numberOfContractions = size(contractions)
    totalNumberOfContractions = MolecularSystem_instance%species(SpecieID)%basisSetSize
   
    unicIntegrals = ((numberOfContractions*(numberOfContractions+1)/2)+1)*(numberOfContractions*(numberOfContractions+1)/2)/2

    if(allocated(contractionId)) deallocate(contractionId)
    if(allocated(contractionLength)) deallocate(contractionLength)
    if(allocated(contractionAngularMoment)) deallocate(contractionAngularMoment)
    if(allocated(contractionNumCartesianOrbital)) deallocate(contractionNumCartesianOrbital)
    if(allocated(contractionOwner)) deallocate(contractionOwner)
    if(allocated(contractionOrigin)) deallocate(contractionOrigin)
    if (allocated(labelsOfContractions)) deallocate(labelsOfContractions)

    allocate(contractionId(numberOfContractions))
    allocate(contractionLength(numberOfContractions))
    allocate(contractionAngularMoment(numberOfContractions))
    allocate(contractionNumCartesianOrbital(numberOfContractions))
    allocate(contractionOwner(numberOfContractions))
    allocate(contractionOrigin(numberOfContractions*3))
    allocate(labelsOfContractions(numberOfContractions))

    primNormalizationSize=0
    aux = 1
    do i=1, numberOfContractions
       contractionId(i) = contractions(i)%id
       contractionLength(i) = contractions(i)%length
       contractionAngularMoment(i)  = contractions(i)%angularMoment
       contractionNumCartesianOrbital(i) = contractions(i)%numCartesianOrbital
       contractionOwner(i) = contractions(i)%owner
       contractionOrigin(i*3-2) = contractions(i)%origin(1)
       contractionOrigin(i*3-1) = contractions(i)%origin(2)
       contractionOrigin(i*3) = contractions(i)%origin(3)
       primNormalizationSize = primNormalizationSize + contractionLength(i)
       !!position for cartesian contractions
       labelsOfContractions(i) = aux
       ! write(*,*) "labels: ", labelsOfContractions(i)
       aux = aux + contractions(i)%numCartesianOrbital          
    end do
 
    maxNumCartesianOrbital = 0
    do i=1, numberOfContractions
       maxNumCartesianOrbital = max(maxNumCartesianOrbital, contractionNumCartesianOrbital(i))
    end do

    if(allocated(contractionOrbitalExponents)) deallocate(contractionOrbitalExponents)
    if(allocated(contractionCoefficients)) deallocate(contractionCoefficients)
    if(allocated(contractionContNormalization)) deallocate(contractionContNormalization)
    if(allocated(contractionPrimNormalization)) deallocate(contractionPrimNormalization)

    allocate(contractionOrbitalExponents(primNormalizationSize))
    allocate(contractionCoefficients(primNormalizationSize))
    allocate(contractionContNormalization(numberOfContractions, maxNumCartesianOrbital))
    allocate(contractionPrimNormalization(primNormalizationSize))

    auxCounter = 1
    do i=1, numberOfContractions
       do j=1, contractionLength(i)
          contractionOrbitalExponents(auxCounter) = contractions(i)%orbitalExponents(j)
          contractionPrimNormalization(auxCounter) = contractions(i)%primNormalization(j,1)
          contractionCoefficients(auxCounter) = contractions(i)%contractionCoefficients(j)
          ! write(*,*) contractions(i)%orbitalExponents(j)
          ! write(*,*) contractionPrimNormalization(auxCounter), contractions(i)%primNormalization(j,1)
          auxCounter = auxCounter + 1
       end do
    end do

    ! write(*,*) ""
    ! write(*,*) "Constantes de normalizacion"
    do i=1, numberOfContractions
       do j=1, contractionNumCartesianOrbital(i) 
          contractionContNormalization(i,j) = contractions(i)%contNormalization(j)
       end do
       if(contractionNumCartesianOrbital(i)<maxNumCartesianOrbital) then
          do j=contractionNumCartesianOrbital(i)+1,maxNumCartesianOrbital
             contractionContNormalization(i,j) = 0.0_8
          end do
       end if
       ! write(*,*) contractionContNormalization(i,:)
    end do

    call MatrixInteger_constructor(auxContIndices, 1, 4)
    call MatrixInteger_constructor(defContIndices, 1, 4)
    call MatrixInteger_constructor(auxPrimIndices, 1, 9)
    call MatrixInteger_constructor(auxNumberOfPPUC, 1, 3)

    totalContIntegrals = 0
    totalPrimitives = 0
    auxCounter2 = 0
    auxCounter3 = 0
    do a = 1, numberOfContractions
       n = a
       do b = a, numberOfContractions
          u = b
          do r = n , numberOfContractions
             do s = u,  numberOfContractions

                !!For (ab|rs)  ---> RESTRICTION a>b && r>s && r+s > a+b
                aux = 0
                order = 0

                !!permuted index
                aa = a
                bb = b
                rr = r
                ss = s

                !!pointer to permuted index under a not permuted loop
                pi => ii
                pj => jj
                pk => kk
                pl => ll

                !!pointer to not permuted index under a permuted loop
                poi => i
                poj => j
                pok => k
                pol => l

                ! write(*,*) contractions(a)%angularMoment, contractions(b)%angularMoment, contractions(r)%angularMoment, contractions(s)%angularMoment

                if (contractions(a)%angularMoment < contractions(b)%angularMoment) then

                   aa = b
                   bb = a

                   pi => jj
                   pj => ii

                   poi => j
                   poj => i

                   order = order + 1

                end if

                if (contractions(r)%angularMoment < contractions(s)%angularMoment) then

                   rr = s
                   ss = r

                   pk => ll
                   pl => kk

                   pok => l
                   pol => k

                   order = order + 3

                end if

                if((contractions(a)%angularMoment + contractions(b)%angularMoment) < &
                     (contractions(r)%angularMoment + contractions(s)%angularMoment)) then

                   aux = aa
                   aa = rr
                   rr = aux

                   aux = bb
                   bb = ss
                   ss = aux

                   select case(order)
                   case(0)
                      pi => kk
                      pj => ll
                      pk => ii
                      pl => jj

                      poi => k
                      poj => l
                      pok => i
                      pol => j

                   case(1)
                      pi => kk
                      pj => ll
                      pk => jj
                      pl => ii

                      poi => l
                      poj => k
                      pok => i
                      pol => j

                   case(3)
                      pi => ll
                      pj => kk
                      pk => ii
                      pl => jj

                      poi => k
                      poj => l
                      pok => j
                      pol => i

                   case(4)
                      pi => ll
                      pj => kk
                      pk => jj
                      pl => ii

                      poi => l
                      poj => k
                      pok => j
                      pol => i

                   end select

                end if

                do i = 1, contractions(aa)%numCartesianOrbital
                   do j = 1, contractions(bb)%numCartesianOrbital
                      do k = 1, contractions(rr)%numCartesianOrbital
                         do l = 1, contractions(ss)%numCartesianOrbital

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
                                  auxCounter=0
                                  auxCounter=contractions(aa)%length*contractions(bb)%length*contractions(rr)%length*contractions(ss)%length
                                  auxCounter2 = auxCounter3
                                  auxCounter3 = auxCounter3 + auxCounter
                                  if(totalContIntegrals==0) then
                                     auxContIndices%values(1,1) = aa
                                     auxContIndices%values(1,2) = bb
                                     auxContIndices%values(1,3) = rr
                                     auxContIndices%values(1,4) = ss
                                     defContIndices%values(1,1) = pa
                                     defContIndices%values(1,2) = pb
                                     defContIndices%values(1,3) = pr
                                     defContIndices%values(1,4) = ps
                                     auxNumberOfPPUC%values(1,1) = auxCounter
                                     auxNumberOfPPUC%values(1,2) = auxCounter2
                                     auxNumberOfPPUC%values(1,3) = auxCounter3
                                  else
                                     call MatrixInteger_addRow(auxContIndices)
                                     call MatrixInteger_addRow(defContIndices)
                                     call MatrixInteger_addRow(auxNumberOfPPUC)
                                     m = 0
                                     m = totalContIntegrals + 1
                                     auxContIndices%values(m,1) = aa
                                     auxContIndices%values(m,2) = bb
                                     auxContIndices%values(m,3) = rr
                                     auxContIndices%values(m,4) = ss
                                     defContIndices%values(m,1) = pa
                                     defContIndices%values(m,2) = pb
                                     defContIndices%values(m,3) = pr
                                     defContIndices%values(m,4) = ps
                                     auxNumberOfPPUC%values(m,1) = auxCounter
                                     auxNumberOfPPUC%values(m,2) = auxCounter2
                                     auxNumberOfPPUC%values(m,3) = auxCounter3
                                  end if

                                  write(*,"(A20,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,A6,I1,A1,I1,A1,I1,A1,I1,A3,A2,I1,A1,I1,A1,I1,A1,I1,A)") &
                                       "Contraida numero: ", &
                                       totalContIntegrals+1, " (", &
                                       auxContIndices%values(totalContIntegrals+1,1), ",", &
                                       auxContIndices%values(totalContIntegrals+1,2), "|", &
                                       auxContIndices%values(totalContIntegrals+1,3), ",", &
                                       auxContIndices%values(totalContIntegrals+1,4), ") ", " ->  (", &
                                       defContIndices%values(totalContIntegrals+1,1), ",", &
                                       defContIndices%values(totalContIntegrals+1,2), "|", &
                                       defContIndices%values(totalContIntegrals+1,3), ",", &
                                       defContIndices%values(totalContIntegrals+1,4), ") :", &
                                       " [", i, ",", j, "|", k, ",", l, "]"

                                  do ii = 1, contractions(aa)%length
                                     do jj = 1, contractions(bb)%length
                                        do kk = 1, contractions(rr)%length
                                           do ll = 1, contractions(ss)%length
                                              if(totalPrimitives == 0) then
                                                 auxPrimIndices%values(1,1) = totalContIntegrals
                                                 auxPrimIndices%values(1,2) = pi
                                                 auxPrimIndices%values(1,3) = pj
                                                 auxPrimIndices%values(1,4) = pk
                                                 auxPrimIndices%values(1,5) = pl
                                                 auxPrimIndices%values(1,6) = i
                                                 auxPrimIndices%values(1,7) = j
                                                 auxPrimIndices%values(1,8) = k
                                                 auxPrimIndices%values(1,9) = l
                                              else
                                                 call MatrixInteger_addRow(auxPrimIndices)
                                                 m = 0
                                                 m = totalPrimitives + 1
                                                 auxPrimIndices%values(m,1) = totalContIntegrals
                                                 auxPrimIndices%values(m,2) = pi
                                                 auxPrimIndices%values(m,3) = pj
                                                 auxPrimIndices%values(m,4) = pk
                                                 auxPrimIndices%values(m,5) = pl
                                                 auxPrimIndices%values(m,6) = i
                                                 auxPrimIndices%values(m,7) = j
                                                 auxPrimIndices%values(m,8) = k
                                                 auxPrimIndices%values(m,9) = l
                                              end if
                                              totalPrimitives = totalPrimitives + 1

                                              ! write(*,"(A30,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,A6,I1,A1,I1,A1,I1,A1,I1,A2)") &
                                              !      "Primitiva numero: ", &
                                              !      totalPrimitives, " (", &
                                              !      auxPrimIndices%values(totalPrimitives,2), ",", &
                                              !      auxPrimIndices%values(totalPrimitives,3), "|", &
                                              !      auxPrimIndices%values(totalPrimitives,4), ",", &
                                              !      auxPrimIndices%values(totalPrimitives,5), ") ", " ->  [", &
                                              !      auxPrimIndices%values(totalPrimitives,6), ",", &
                                              !      auxPrimIndices%values(totalPrimitives,7), "|", &
                                              !      auxPrimIndices%values(totalPrimitives,8), ",", &
                                              !      auxPrimIndices%values(totalPrimitives,9), "] "
                                           end do
                                        end do
                                     end do
                                  end do
                                  totalContIntegrals = totalContIntegrals + 1
                               end if
                            end if

                         end do
                      end do
                   end do
                end do !! Primitives loop
             end do
             u=r+1
          end do
       end do
    end do !! shells loop



    if(allocated(contractionIntegrals)) deallocate(contractionIntegrals)
    if(allocated(contractionIndices)) deallocate(contractionIndices)
    if(allocated(primitiveIndices)) deallocate(primitiveIndices)
    if(allocated(numberOfPPUC)) deallocate(numberOfPPUC)

    allocate(contractionIntegrals(totalContIntegrals))
    allocate(contractionIndices(totalContIntegrals*4))
    allocate(primitiveIndices(totalPrimitives*9))
    allocate(numberOfPPUC(totalContIntegrals*3))

    do i = 1, totalContIntegrals
       contractionIndices(i*4-3) = auxContIndices%values(i,1)
       contractionIndices(i*4-2) = auxContIndices%values(i,2)
       contractionIndices(i*4-1) = auxContIndices%values(i,3)
       contractionIndices(i*4) = auxContIndices%values(i,4)
       numberOfPPUC(i*3-2) = auxNumberOfPPUC%values(i,1)
       numberOfPPUC(i*3-1) = auxNumberOfPPUC%values(i,2)
       numberOfPPUC(i*3) = auxNumberOfPPUC%values(i,3)
       ! write(*,"(A20,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,A6,I1,A1,I1,A1,I1,A1,I1,A)") &
       !      "Contraida numero: ", &
       !      i, " (", &
       !      contractionIndices(i*4-3), ",", &
       !      contractionIndices(i*4-2), "|", &
       !      contractionIndices(i*4-1), ",", &
       !      contractionIndices(i*4), ") ", " ->  (", &
       !      defContIndices%values(i,1), ",", &
       !      defContIndices%values(i,2), "|", &
       !      defContIndices%values(i,3), ",", &
       !      defContIndices%values(i,4), ")"
    end do

    do i = 1, totalPrimitives
       primitiveIndices(i*9-8) = auxPrimIndices%values(i,1)
       primitiveIndices(i*9-7) = auxPrimIndices%values(i,2)
       primitiveIndices(i*9-6) = auxPrimIndices%values(i,3)
       primitiveIndices(i*9-5) = auxPrimIndices%values(i,4)
       primitiveIndices(i*9-4) = auxPrimIndices%values(i,5)
       primitiveIndices(i*9-3) = auxPrimIndices%values(i,6)
       primitiveIndices(i*9-2) = auxPrimIndices%values(i,7)
       primitiveIndices(i*9-1) = auxPrimIndices%values(i,8)
       primitiveIndices(i*9) = auxPrimIndices%values(i,9)
       ! write(*,"(A30,I3,A2,I1,A1,I1,A1,I1,A1,I1,A2,A6,I1,A1,I1,A1,I1,A1,I1,A2)") &
       !      "Primitiva numero: ", &
       !      i, " (", &
       !      primitiveIndices(i*9-7), ",", &
       !      primitiveIndices(i*9-6), "|", &
       !      primitiveIndices(i*9-5), ",", &
       !      primitiveIndices(i*9-4), ") ", " ->  [", &
       !      primitiveIndices(i*9-3), ",", &
       !      primitiveIndices(i*9-2), "|", &
       !      primitiveIndices(i*9-1), ",", &
       !      primitiveIndices(i*9), "] "
    end do


    call cuda_int_intraspecies(&
         numberOfContractions, &
         totalContIntegrals, &
         totalPrimitives, &
         maxNumCartesianOrbital, &
         primNormalizationSize, &
         contractionId, &
         contractionLength, &
         contractionAngularMoment, &
         contractionNumCartesianOrbital, &
         contractionOwner, &
         contractionOrigin, &
         contractionOrbitalExponents, &
         contractionCoefficients, &
         contractionContNormalization, &
         contractionPrimNormalization, &
         contractionIntegrals, &
         contractionIndices, &
         primitiveIndices, &
         numberOfPPUC, &
         labelsOfContractions &
         )

    auxCounter = 0
    counter = 0
    ! write(*,*) "Contraida"
    do i=1, unicIntegrals

      
       if(abs(contractionIntegrals(i)) > 1.0D-10) then

          auxCounter = auxCounter + 1
          counter = counter + 1

          eris%a(counter) = defContIndices%values(i,1)
          eris%b(counter) = defContIndices%values(i,2)
          eris%c(counter) = defContIndices%values(i,3)
          eris%d(counter) = defContIndices%values(i,4)
          eris%integrals(counter) = contractionIntegrals(i)
          ! write(*,*) "(",eris%a(counter),",",eris%b(counter),"|",eris%c(counter),",",eris%d(counter),")", ": ", contractionIntegrals(i)
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
    end do

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

    
  end subroutine CudintInterface_computeIntraSpecies
  
 !>
  !! @brief  Maneja excepciones de la clase
  subroutine CudintInterface_exception( typeMessage, description, debugDescription)
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

  end subroutine CudintInterface_exception

end module CudintInterface_
