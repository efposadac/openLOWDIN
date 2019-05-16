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
module TransformIntegralsA_
  use MolecularSystem_
  use InputManager_
  use IntegralManager_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
  use omp_lib
  implicit none

  type, public :: TransformIntegralsA
     character(30) :: name
     character(255) :: fileForCoefficients
     character(255) :: fileForIntegrals
     character(255) :: prefixOfFile
     integer :: numberOfContractions
     integer :: bias
     integer :: specieID
     integer :: otherSpecieID
     integer :: unidOfOutputForCoefficients
     integer :: unidOfOutputForIntegrals
     integer :: nproc
     integer :: integralStackSize

  end type TransformIntegralsA

  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE = 0
  integer, parameter :: TWO_SPECIES = 1
  !! }

  public :: &
       TransformIntegralsA_constructor, &
       TransformIntegralsA_destructor, &
       TransformIntegralsA_show, &
       TransformIntegralsA_atomicToMolecularOfOneSpecie, &
       TransformIntegralsA_atomicToMolecularOfTwoSpecies
  !       TransformIntegralsA_readIntegralsTransformed

  private

  interface

     !**                                                                                                                                                                                                
     ! Realiza proceso de minimizacion restringida o no de una funcion arbitraria                                                                                                                       
     ! tomado de:    Zhu, C.; Lu, P. and Nocedal, J. Department of Electrical Engineering                                                                                                               
     !               and Computer Science. Northwestern University. 1996                                                                                                                                
     !               L-BFGS-B FORTRAN SUBROUTINES FOR LARGE-SCALE BOUND CONSTRAINED                                                                                                                     
     !               OPTIMIZATION                                                                                                                                                                       
     !**                                                                                                                                                                                                
     subroutine  setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa, task, iprint,  csave, lsave, isave, dsave)
       character*60    ::  task, csave
       logical ::   lsave(4)
       integer :: n, m, iprint, nbd(*), iwa(*), isave(44)
       real(8) :: f, factr, pgtol, x(*), l(*), u(*), g(*), wa(*), dsave(29)
     end subroutine setulb

     !**                                                                                                                                                                                                
     ! Realiza transformacion de cuatro indices de integrales en OA a integrales                                                                                                                        
     ! en orbitales moleculares.                                                                                                                                                                        
     ! tomado de:    Yamamoto, Shigeyoshi; Nagashima, U.                                                                                                                                                
     !               Four-index integral tranformation exploiting symmetry.                                                                                                                             
     !               Computer Physics Communications, 2005, 166, 58-65                                                                                                                                  
     !**                                                                                                                                                                                                
     subroutine  fourIndexTransformation(numberOfContractions, otherNumberOfContractions, nameFile, nproc, integralStackSize )
       integer :: numberOfContractions
       integer :: otherNumberOfContractions
       integer :: nproc
       integer :: integralStackSize
       character(*) :: nameFile
     end subroutine fourIndexTransformation

  end interface


contains


  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsA_constructor(this)
    implicit none
    type(TransformIntegralsA) :: this

    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"

  end subroutine TransformIntegralsA_constructor


  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsA_destructor(this)
    implicit none
    type(TransformIntegralsA) :: this

  end subroutine TransformIntegralsA_destructor

  !>
  !! @brief show
  !<
  subroutine TransformIntegralsA_show()
    implicit none

    print *,"--------------------------------------------------"
    print *,"    Algorithm Four-index integral tranformation"
    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
    print *,"  Computer Physics Communications, 2005, 166, 58-65"
    print *,"--------------------------------------------------"
    print *,""

  end subroutine TransformIntegralsA_show


  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsA_atomicToMolecularOfOneSpecie( this, coefficientsOfAtomicOrbitals, &
       molecularIntegrals, specieID, nameOfSpecie  )
    implicit none
    type(TransformIntegralsA) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID
    character(*) :: nameOfSpecie
    integer :: nproc
    integer :: integralStackSize
    integer :: errorNum
    real(8) :: initialTime
    real(8) :: finalTime

    ! Reads the number of cores
    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    if ( .not.CONTROL_instance%OPTIMIZE ) then
       call cpu_time(initialTime)
    end if

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%specieID = specieID

    call TransformIntegralsA_writeCoefficients(this, coefficientsOfAtomicOrbitals)

    !! Inicia proceso de transformacion
    call fourIndexTransformation( this%numberOfContractions, 0_4,  trim(this%prefixOfFile ), nproc, integralStackSize )

    !    !! Lee  de disco las integrales tranformadas
    !    call TransformIntegralsA_readIntegralsTransformed( this, molecularIntegrals, ONE_SPECIE )

    !    !! Remueve archivos empleados en proceso de transformacion
    !    call system("rm "// trim(this%prefixOfFile)//"*.dat "// trim(this%prefixOfFile) // "*.values "  )

    !    if ( .not.CONTROL_instance%OPTIMIZE ) then
    !       call cpu_time(finalTime)
    !       write (6,"(T15,A30,ES10.2,A4)") "cpu-time  for transformation:  ", finalTime-initialTime ," (s)"
    !       print *,""
    !    end if

  end subroutine TransformIntegralsA_atomicToMolecularOfOneSpecie

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsA_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
    implicit none
    type(TransformIntegralsA) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID
    character(*) :: nameOfSpecie
    integer :: otherSpecieID
    character(*) :: nameOfOtherSpecie
    integer :: nproc
    integer :: integralStackSize
    integer :: errorNum
    real(8) :: initialTime
    real(8) :: finalTime

    if ( .not.CONTROL_instance%OPTIMIZE ) then
       call cpu_time(initialTime)
    end if

    ! Reads the number of cores
    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE


    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"


    this%numberOfContractions = size(coefficientsOfAtomicOrbitals%values, dim=1)+size(otherCoefficientsOfAtomicOrbitals%values, dim=1)

    this%bias = size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%specieID = specieID
    this%otherSpecieID = otherSpecieID


    call TransformIntegralsA_writeCoefficients( this, coefficientsOfAtomicOrbitals, otherCoefficientsOfAtomicOrbitals )

    !! Inicia proceso de transformacion
    !! this%numberOfContractions = Total number of contractions, it is the sum of contractions beetwen specieID and otherSpecieID
    call fourIndexTransformation( this%numberOfContractions, size(coefficientsOfAtomicOrbitals%values,dim=1), trim(this%prefixOfFile), 0_4, integralStackSize )

    ! Lee  de disco las integrales tranformadas
    !    call TransformIntegralsA_readIntegralsTransformed( this, molecularIntegrals, TWO_SPECIES )

    !! Remueve archivos empleados en proceso de transformacion
    !    call system("rm "// trim(this%prefixOfFile)//"*.dat "// trim(this%prefixOfFile) // "*.values "  )

    !    if ( .not.CONTROL_instance%OPTIMIZE ) then
    !       call cpu_time(finalTime)
    !       write (6,"(T15,A30,ES10.2,A4)") "cpu-time  for transformation:  ", finalTime-initialTime ," (s)"
    !       print *,""
    !    end if

  end subroutine TransformIntegralsA_atomicToMolecularOfTwoSpecies

  !>
  !! @brief Escribe los coefficientes de combinacion para los orbitales atomicos.
  !!    El almacenamiento requiere guardar columnas completas una tras de otra
  !<
  subroutine TransformIntegralsA_writeCoefficients( this, coefficients, otherCoefficients )
    implicit none
    type(TransformIntegralsA) :: this
    type(Matrix) :: coefficients
    type(Matrix), optional :: otherCoefficients
    integer :: a
    integer :: b

    open( UNIT=this%unidOfOutputForCoefficients,FILE=trim(this%fileforcoefficients),STATUS='REPLACE', &
         ACCESS='SEQUENTIAL', FORM='FORMATTED' )

    if ( .not.present(otherCoefficients) ) then

       this%numberOfContractions=size(coefficients%values,dim=1)

       do a=1, this%numberOfContractions
          do b=1,this%numberOfContractions

             write(this%unidOfOutputForCoefficients,*) a,b,coefficients%values(b,a)

          end do
       end do

    else

       this%numberOfContractions=size(coefficients%values,dim=1)+size(otherCoefficients%values,dim=1)
       this%bias = size(coefficients%values,dim=1)

       !*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! Escribe en disco los coeficientes de combinacion  para un par de especies,  haciendo
       !! un "this%bias" de uno de los conjuntos sobre el otro
       !!
       do a=1, this%numberOfContractions
          do b=1,this%numberOfContractions

             if ( ( a <= this%bias ) .and. ( b <= this%bias ) ) then
                write(this%unidOfOutputForCoefficients,*) a, b, coefficients%values( b, a )

             else &
                  if ( ( a > this%bias ) .and. ( b > this%bias ) ) then
                write(this%unidOfOutputForCoefficients,*) a, b, otherCoefficients%values( b-this%bias, a-this%bias )

             else
                write(this%unidOfOutputForCoefficients,*) a,b,0.0_8

             end if

          end do
       end do

    end if

    close( UNIT=this%unidOfOutputForCoefficients )


  end subroutine TransformIntegralsA_writeCoefficients

  !! This subroutine was moved to lowdinCore  
  !
  !  subroutine TransformIntegralsA_readIntegralsTransformed(this, matrixContainer, typeOfIntegrals )
  !    implicit none
  !    type(TransformIntegralsA) :: this
  !    type(Matrix) :: matrixContainer
  !    integer :: typeOfIntegrals
  !
  !    integer(8) :: numberOfIntegrals
  !    integer(8) :: auxIndex
  !    real(8),dimension(791) :: integralValue
  !    real(8) :: auxValue
  !    integer :: iter
  !    integer :: errorValue
  !    integer :: bufferA
  !    integer :: bufferB
  !    integer :: bufferSize
  !    integer :: lowerIndices(2), upperIndeces(2), counter(2)
  !    integer,dimension(4,791) :: indexBuffer
  !
  !
  !    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
  !    open(unit=this%unidOfOutputForIntegrals, file=trim(this%prefixOfFile)//"moint.dat", &
  !         status='old',access='sequential', form='unformatted' )
  !
  !    select case( typeOfIntegrals)
  !
  !    case(ONE_SPECIE)
  !
  !       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)
  !
  !       numberOfIntegrals   =  int( ( (  this%numberOfContractions * (  this%numberOfContractions + 1.0_8 ) / 4.0_8 ) * &
  !            ( (  this%numberOfContractions * (  this%numberOfContractions + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )
  !
  !       call Matrix_constructor( matrixContainer, numberOfIntegrals, 1_8, 0.0_8 )
  !       matrixContainer%values = 0.0_8
  !
  !       do
  !          read(UNIT=this%unidOfOutputForIntegrals,IOSTAT=errorValue) bufferA,bufferB,integralValue,indexBuffer
  !
  !          bufferSize = iabs( bufferA)
  !          if ( bufferA /= 0 ) then
  !             do iter=1,bufferSize
  !
  !                auxIndex = IndexMap_tensorR4ToVector(indexBuffer(1,iter),indexBuffer(2,iter), &
  !                     indexBuffer(3,iter), indexBuffer(4,iter), this%numberOfContractions )
  !
  !                matrixContainer%values( auxIndex, 1 ) = integralValue(iter)
  !
  !             end do
  !          end if
  !
  !          if ( bufferA <= 0 ) exit
  !
  !       end do
  !
  !
  !    case(TWO_SPECIES)
  !
  !
  !       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)
  !
  !       numberOfIntegrals = ( this%bias    *  ( ( this%bias + 1.0_8) / 2.0_8 ) ) * &
  !            ( (this%numberOfContractions-this%bias) * ( ( (this%numberOfContractions-this%bias) + 1.0_8 ) / 2.0_8 ) )
  !
  !       call Matrix_constructor( matrixContainer, numberOfIntegrals, 1_8, 0.0_8 )
  !
  !       matrixContainer%values = 0.0_8
  !
  !       do
  !          read(UNIT=this%unidOfOutputForIntegrals,IOSTAT=errorValue) bufferA,bufferB,integralValue,indexBuffer
  !          bufferSize = iabs( bufferA)
  !          if ( bufferA /= 0 ) then
  !             do iter=1,bufferSize
  !                counter=1
  !                if ( indexBuffer(1,iter ) > this%bias )  then
  !                   indexBuffer(1,iter)=indexBuffer(1,iter)-this%bias
  !                   upperIndeces( counter(1) ) = indexBuffer(1,iter)
  !                   counter(1)=2
  !                else
  !                   lowerIndices( counter(2) ) = indexBuffer(1,iter)
  !                   counter(2)=2
  !                end if
  !
  !                if ( indexBuffer(2,iter ) > this%bias ) then
  !                   indexBuffer(2,iter)=indexBuffer(2,iter)-this%bias
  !                   upperIndeces( counter(1) ) = indexBuffer(2,iter)
  !                   counter(1)=2
  !                else
  !                   lowerIndices( counter(2) ) = indexBuffer(2,iter)
  !                   counter(2)=2
  !                end if
  !
  !                if ( indexBuffer(3,iter ) > this%bias ) then
  !                   indexBuffer(3,iter)=indexBuffer(3,iter)-this%bias
  !                   upperIndeces( counter(1) ) = indexBuffer(3,iter)
  !                   counter(1)=2
  !                else
  !                   lowerIndices( counter(2) ) = indexBuffer(3,iter)
  !                   counter(2)=2
  !                end if
  !
  !                if ( indexBuffer(4,iter ) > this%bias ) then
  !                   indexBuffer(4,iter)=indexBuffer(4,iter)-this%bias
  !                   upperIndeces( counter(1) ) = indexBuffer(4,iter)
  !                   counter(1)=2
  !                else
  !                   lowerIndices( counter(2) ) = indexBuffer(4,iter)
  !                   counter(2)=2
  !                end if
  !
  !
  !
  !                auxIndex = IndexMap_tensorR4ToVector(lowerIndices(1),lowerIndices(2),upperIndeces(1),upperIndeces(2), &
  !                     this%bias, this%numberOfContractions - this%bias )
  !
  !                matrixContainer%values(auxIndex,1)=integralValue(iter)
  !
  !             end do
  !          end if
  !
  !          if ( bufferA <= 0 ) exit
  !
  !       end do
  !
  !    end select
  !
  !    close(this%unidOfOutputForIntegrals)
  !
  !  end subroutine TransformIntegralsA_readIntegralsTransformed


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine TransformIntegralsA_exception( typeMessage, description, debugDescription)
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

  end subroutine TransformIntegralsA_exception
  
end module TransformIntegralsA_
