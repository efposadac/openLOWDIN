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
!! @brief Clase encargada de realizar transformacion de integrales atomicas a  moleculares
!!
!! 	Esta clase reliza la transformacion de integrales de orbitales atomicos a orbitales moleculares,
!!	creando una interface al algoritmo de   Yamamoto, Shigeyoshi; Nagashima, Umpei.
!!	Computer Physics Communications, 2005, 166, 58-65
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
  use InputManager_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
  use omp_lib
  implicit none

  type, public :: TransformIntegralsC
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
     
     integer :: lowerOccupiedOrbital
     integer :: upperOccupiedOrbital
     integer :: lowerVirtualOrbital
     integer :: upperVirtualOrbital

  end type TransformIntegralsC

  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE			= 0
  integer, parameter :: TWO_SPECIES			= 1
  !! }

  public :: &
       TransformIntegralsC_constructor, &
       TransformIntegralsC_destructor, &
       TransformIntegralsC_show, &
       TransformIntegralsC_atomicToMolecularOfOneSpecie!, &
!       TransformIntegralsC_atomicToMolecularOfTwoSpecies
!       TransformIntegralsC_readIntegralsTransformed

  private

contains

  
  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsC_constructor(this)
    implicit none
    type(TransformIntegralsC) :: this

    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"



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

     print *,""
     print *,"BEGIN FOUR-INDEX INTEGRALS TRANFORMATION:"
     print *,"========================================"
     print *,""
     print *,"--------------------------------------------------"
     print *,"   4N^5 Algorithm Four-index integral tranformation"
     print *,"--------------------------------------------------"
     print *,""

  end subroutine TransformIntegralsC_show
  
  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !! 		a integrales moleculares.
  !<
  subroutine TransformIntegralsC_atomicToMolecularOfOneSpecie( this, coefficientsOfAtomicOrbitals, &
       molecularIntegrals, specieID, nameOfSpecie  )
    implicit none
    type(TransformIntegralsC) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID
    character(*) :: nameOfSpecie
    integer :: nproc
    integer :: integralStackSize
    real(8) :: initialTime
    real(8) :: finalTime

    integer :: ifile, i
    integer :: unit
    character(50) :: sfile
    integer :: status

    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer*2 :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m, n, u, nn, uu, qq, rr, ss

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

    call TransformIntegralsC_checkMOIntegralType(specieID, this)

    if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    allocate (twoParticlesIntegrals ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )

    twoParticlesIntegrals = 0

    do ifile = 1, nproc

      write(sfile,*) ifile
      sfile = trim(adjustl(sfile))
      unit = ifile+50

      open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".ints", status='old',access='sequential', form='Unformatted')

      loadintegrals : do
  
        read(UNIT=unit, iostat=status) aa(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
              bb(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
              cc(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
              dd(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
              shellIntegrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
    
    
         do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
    
            if( aa(i) == -1 ) exit loadintegrals

            twoParticlesIntegrals(aa(i),bb(i),cc(i),dd(i)) = shellIntegrals(i)
    
         end do

       end do loadintegrals
       close (unit)

    end do 

   !! symmetrize 
    do mu = 1, this%numberOfContractions
      do nu = 1, this%numberOfContractions
        do lambda = 1, this%numberOfContractions
          do sigma = 1, this%numberOfContractions
            if (abs(twoParticlesIntegrals(mu,nu,lambda,sigma)) < 1E-10 ) cycle
            twoParticlesIntegrals(nu,mu,lambda,sigma) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
            twoParticlesIntegrals(mu,nu,sigma,lambda) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
            twoParticlesIntegrals(lambda,sigma,mu,nu) = twoParticlesIntegrals(mu,nu,lambda,sigma) 

          end do
        end do  
      end do  
    end do  

    if ( allocated (tempA)) deallocate (tempA )
    allocate (tempA ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )

    tempA = 0

    if ( allocated (tempB)) deallocate (tempB )
    allocate (tempB ( this%numberOfContractions , &
                                      this%numberOfContractions ) )
    tempB = 0

    if ( allocated (tempC)) deallocate (tempC )
    allocate (tempC ( this%numberOfContractions ) )

    tempC = 0

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )
     
    m = 0
    do p = 1, this%numberOfContractions
      n = p
      tempA = 0
      auxTransformedTwoParticlesIntegral = 0

      if ( p >  this%upperOccupiedOrbital  ) cycle

      !! First quarter
      do mu = 1, this%numberOfContractions
        if ( abs(coefficientsOfAtomicOrbitals%values( mu, p )) < 1E-10 ) cycle
        tempA(:,:,:) = tempA(:,:,:) + coefficientsOfAtomicOrbitals%values( mu, p )* &
                                        twoParticlesIntegrals(mu,:,:,: ) 
      end do

      do q = p, this%numberOfContractions
        u = q
        tempB = 0

        if ( q < this%lowerVirtualOrbital ) cycle
        !! second quarter
        do nu = 1, this%numberOfContractions
          if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
          tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                                        tempA(nu,:,:)
        end do

        do r = n, this%numberOfContractions
           tempC = 0

           if ( r >  this%upperOccupiedOrbital  ) cycle

           !! third quarter
           do lambda = 1, this%numberOfContractions
             tempC(:) = tempC(:) + coefficientsOfAtomicOrbitals%values( lambda, r )* &
                                            tempB(lambda,:)

           end do
           do s = u, this%numberOfContractions
             auxTransformedTwoParticlesIntegral = 0

             if ( s < this%lowerVirtualOrbital ) cycle
             !! fourth quarter
             do sigma = 1, this%numberOfContractions
               auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                                                         coefficientsOfAtomicOrbitals%values( sigma, s )* &
                                                         tempC(sigma)

             end do

             write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral

           end do
           u = r + 1
         end do
       end do
     end do

     write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0  

  end subroutine TransformIntegralsC_atomicToMolecularOfOneSpecie





  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !! 		a integrales moleculares.
  !<
!!  subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
!!       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
!!    implicit none
!!    type(TransformIntegralsC) :: this
!!    type(Matrix) :: coefficientsOfAtomicOrbitals
!!    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
!!    type(Matrix) :: molecularIntegrals
!!    integer :: specieID
!!    character(*) :: nameOfSpecie
!!    integer :: otherSpecieID
!!    character(*) :: nameOfOtherSpecie
!!    integer :: nproc
!!    integer :: integralStackSize
!!    integer :: errorNum
!!    real(8) :: initialTime
!!    real(8) :: finalTime
!!    
!!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!!       call cpu_time(initialTime)
!!    end if
!!
!!    ! Reads the number of cores
!!    nproc = CONTROL_instance%NUMBER_OF_CORES
!!    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE
!!
!!    
!!    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
!!    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"
!!
!!
!!    this%numberOfContractions = size(coefficientsOfAtomicOrbitals%values, dim=1)+size(otherCoefficientsOfAtomicOrbitals%values, dim=1)
!!
!!    this%bias = size(coefficientsOfAtomicOrbitals%values,dim=1)
!!    this%specieID = specieID
!!    this%otherSpecieID = otherSpecieID
!!
!!
!!    call TransformIntegralsC_writeCoefficients( this, coefficientsOfAtomicOrbitals, otherCoefficientsOfAtomicOrbitals )
!!
!!    !! Inicia proceso de transformacion
!!    !! this%numberOfContractions = Total number of contractions, it is the sum of contractions beetwen specieID and otherSpecieID
!!    call fourIndexTransformation( this%numberOfContractions, size(coefficientsOfAtomicOrbitals%values,dim=1), trim(this%prefixOfFile), 0_4, integralStackSize )
!!
!!    ! Lee  de disco las integrales tranformadas
!!!    call TransformIntegralsC_readIntegralsTransformed( this, molecularIntegrals, TWO_SPECIES )
!!
!!    !! Remueve archivos empleados en proceso de transformacion
!!!    call system("rm "// trim(this%prefixOfFile)//"*.dat "// trim(this%prefixOfFile) // "*.values "  )
!!
!!!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!!!       call cpu_time(finalTime)
!!!       write (6,"(T15,A30,ES10.2,A4)") "cpu-time  for transformation:  ", finalTime-initialTime ," (s)"
!!!       print *,""
!!!    end if
!!
!!  end subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies
!!
!!  !>
!!  !! @brief Escribe los coefficientes de combinacion para los orbitales atomicos.
!!  !! 		El almacenamiento requiere guardar columnas completas una tras de otra
!!  !<
!!  subroutine TransformIntegralsC_writeCoefficients( this, coefficients, otherCoefficients )
!!    implicit none
!!    type(TransformIntegralsC) :: this
!!    type(Matrix) :: coefficients
!!    type(Matrix), optional :: otherCoefficients
!!    integer :: a
!!    integer :: b
!!
!!    open( UNIT=this%unidOfOutputForCoefficients,FILE=trim(this%fileforcoefficients),STATUS='REPLACE', &
!!         ACCESS='SEQUENTIAL', FORM='FORMATTED' )
!!
!!    if ( .not.present(otherCoefficients) ) then
!!
!!       this%numberOfContractions=size(coefficients%values,dim=1)
!!
!!       do a=1, this%numberOfContractions
!!          do b=1,this%numberOfContractions
!!
!!             write(this%unidOfOutputForCoefficients,*) a,b,coefficients%values(b,a)
!!
!!          end do
!!       end do
!!
!!    else
!!
!!       this%numberOfContractions=size(coefficients%values,dim=1)+size(otherCoefficients%values,dim=1)
!!       this%bias = size(coefficients%values,dim=1)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       !! Escribe en disco los coeficientes de combinacion  para un par de especies,  haciendo
!!       !! un "this%bias" de uno de los conjuntos sobre el otro
!!       !!
!!       do a=1, this%numberOfContractions
!!          do b=1,this%numberOfContractions
!!
!!             if ( ( a <= this%bias ) .and. ( b <= this%bias ) ) then
!!                write(this%unidOfOutputForCoefficients,*) a, b, coefficients%values( b, a )
!!
!!             else &
!!                  if ( ( a > this%bias ) .and. ( b > this%bias ) ) then
!!                write(this%unidOfOutputForCoefficients,*) a, b, otherCoefficients%values( b-this%bias, a-this%bias )
!!
!!             else
!!                write(this%unidOfOutputForCoefficients,*) a,b,0.0_8
!!
!!             end if
!!
!!          end do
!!       end do
!!
!!    end if
!!
!!    close( UNIT=this%unidOfOutputForCoefficients )
!!
!!
!!  end subroutine TransformIntegralsC_writeCoefficients

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsC_checkMOIntegralType(speciesID, this)
    implicit none
    integer :: speciesID
    type(TransformIntegralsC) :: this
    integer :: totalOccupation 
    integer :: totalNumberOfContractions

    totalOccupation = MolecularSystem_getOcupationNumber( speciesID )
    totalNumberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID)

    !! All orbitals. Default
    this%lowerOccupiedOrbital = 1 
    this%upperOccupiedOrbital = totalNumberOfContractions
    this%lowerVirtualOrbital = 1
    this%upperVirtualOrbital = totalNumberOfContractions

    !! only the (ia|jb) integrals will be transformed
    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION == 2 .or. &
         ( CONTROL_instance%PT_ORDER == 2 .and.  CONTROL_instance%IONIZE_MO <= totalOCcupation ) ) then
      this%lowerOccupiedOrbital = 1 
      this%upperOccupiedOrbital = totalOccupation
      this%lowerVirtualOrbital = totalOccupation + 1
      this%upperVirtualOrbital = totalNumberOfContractions
    end if


  end subroutine TransformIntegralsC_checkMOIntegralType
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
