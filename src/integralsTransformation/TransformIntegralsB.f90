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
module TransformIntegralsB_
  use MolecularSystem_
  use InputManager_
  use ParticleManager_
  use Matrix_
  use IndexMap_
  use Exception_
  use omp_lib
  implicit none

  type, public :: TransformIntegralsB
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

     integer :: lowerOccupiedOrbital
     integer :: upperOccupiedOrbital
     integer :: lowerVirtualOrbital
     integer :: upperVirtualOrbital

  end type TransformIntegralsB

  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE = 0
  integer, parameter :: TWO_SPECIES = 1
  !! }

  public :: &
       TransformIntegralsB_constructor, &
       TransformIntegralsB_destructor, &
       TransformIntegralsB_show, &
       TransformIntegralsB_atomicToMolecularOfOneSpecie, &
       TransformIntegralsB_atomicToMolecularOfTwoSpecies
  !       TransformIntegralsB_readIntegralsTransformed

  private

contains


  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsB_constructor(this)
    implicit none
    type(TransformIntegralsB) :: this

    this%unidOfOutputForCoefficients = CONTROL_instance%UNIT_FOR_MOLECULAR_ORBITALS_FILE
    this%unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE
    this%fileForIntegrals = trim(CONTROL_INSTANCE%INPUT_FILE)//".ints"



  end subroutine TransformIntegralsB_constructor

  !>
  !! @brief Contructor de la clase
  !<
  subroutine TransformIntegralsB_destructor(this)
    implicit none
    type(TransformIntegralsB) :: this

  end subroutine TransformIntegralsB_destructor

  !>
  !! @brief show
  !<
  subroutine TransformIntegralsB_show()
    implicit none

    print *,""
    print *,"BEGIN FOUR-INDEX INTEGRALS TRANFORMATION:"
    print *,"========================================"
    print *,""
    print *,"--------------------------------------------------"
    print *,"    N^8 Algorithm Four-index integral tranformation"
    print *,"--------------------------------------------------"
    print *,""

  end subroutine TransformIntegralsB_show


  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de la misma especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsB_atomicToMolecularOfOneSpecie( this, coefficientsOfAtomicOrbitals, &
       molecularIntegrals, specieID, nameOfSpecie  )
    implicit none
    type(TransformIntegralsB) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID
    character(*) :: nameOfSpecie
    integer :: nproc
    integer :: integralStackSize

    integer :: ifile, i
    integer :: unit
    character(50) :: sfile
    integer :: status

    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8)  auxTransformedTwoParticlesIntegral

    integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m, n, u

    ! Reads the number of cores
    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%specieID = specieID

    this%lowerOccupiedOrbital = 1 
    this%upperOccupiedOrbital = MolecularSystem_getOcupationNumber( specieID )
    this%lowerVirtualOrbital = MolecularSystem_getOcupationNumber( specieID ) + 1
    this%upperVirtualOrbital = this%numberOfContractions

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

       open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".ints", status='old',access='stream', form='Unformatted')

       loadintegrals : do

          read(UNIT=unit, iostat=status) aa(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               bb(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               rr(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               ss(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
               shellIntegrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)


          do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

             if( aa(i) == -1 ) exit loadintegrals

             twoParticlesIntegrals(aa(i),bb(i),rr(i),ss(i)) = shellIntegrals(i)

          end do

       end do loadintegrals
       close (unit)

    end do

    !! symmetrize 
    do mu = 1, this%numberOfContractions
       do nu = 1, this%numberOfContractions
          do lambda = 1, this%numberOfContractions
             do sigma = 1, this%numberOfContractions
                twoParticlesIntegrals(nu,mu,lambda,sigma) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
                twoParticlesIntegrals(mu,nu,sigma,lambda) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
                twoParticlesIntegrals(lambda,sigma,mu,nu) = twoParticlesIntegrals(mu,nu,lambda,sigma) 

             end do
          end do
       end do
    end do

    !!    print *, "this 0", this%lowerOccupiedOrbital
    !!    print *, "this 0", this%upperOccupiedOrbital
    !!    print *, "this 0", this%lowerVirtualOrbital 
    !!    print *, "this 0", this%upperVirtualOrbital 

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    do p = 1, this%numberOfContractions
       n = p
       do q = p, this%numberOfContractions
          u = q 
          do r = n, this%numberOfContractions
             do s = u, this%numberOfContractions

                if ( q >= this%lowerVirtualOrbital .and. s >= this%lowerVirtualOrbital .and. &
                     p <= this%upperOccupiedOrbital .and. r <= this%upperOccupiedOrbital ) then

                   auxTransformedTwoParticlesIntegral = 0
                   do mu = 1, this%numberOfContractions
                      do nu = 1, this%numberOfContractions
                         do lambda = 1, this%numberOfContractions
                            do sigma = 1, this%numberOfContractions

                               auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                                    coefficientsOfAtomicOrbitals%values( mu, p )* &
                                    coefficientsOfAtomicOrbitals%values( nu, q )* &
                                    coefficientsOfAtomicOrbitals%values( lambda, r )* &
                                    coefficientsOfAtomicOrbitals%values( sigma, s )* &
                                    twoParticlesIntegrals(mu, nu, lambda, sigma) 

                            end do
                         end do
                      end do
                   end do

                   write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral

                end if

             end do
             u = r + 1
          end do
       end do
    end do

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0  

  end subroutine TransformIntegralsB_atomicToMolecularOfOneSpecie

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !!    a integrales moleculares.
  !<
  subroutine TransformIntegralsB_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
    implicit none
    type(TransformIntegralsB) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID, otherSpecieID
    character(*) :: nameOfSpecie, nameOfOtherSpecie
    integer :: nproc
    integer :: integralStackSize

    integer :: i

    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8)  auxTransformedTwoParticlesIntegral

    integer :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m

    ! Reads the number of cores

    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    this%prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//"mo.values"

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    this%otherNumberOfContractions=size(otherCoefficientsOfAtomicOrbitals%values,dim=1)

    this%specieID = specieID


    if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    allocate (twoParticlesIntegrals (  this%numberOfContractions, &
         this%numberOfContractions, &
         this%otherNumberOfContractions , &
         this%otherNumberOfContractions )  )

    twoParticlesIntegrals = 0

    m = 0

    !! Read integrals

    !! open file for integrals
    open(UNIT=34,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

    loadintegrals : do

       read(34)   aa(1:CONTROL_instance%INTEGRAL_STACK_SIZE), bb(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            cc(1:CONTROL_instance%INTEGRAL_STACK_SIZE), dd(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            shellIntegrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

       do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

          if (aa(i) == -1) exit loadintegrals

          m = m + 1
          twoParticlesIntegrals(aa(i),bb(i),cc(i),dd(i)) = shellIntegrals(i)

       end do

    end do loadintegrals

    close (34)

    !! symmetrize 
    do mu = 1, this%numberOfContractions
       do nu = 1, this%numberOfContractions
          do lambda = 1, this%otherNumberOfContractions
             do sigma = 1, this%otherNumberOfContractions
                twoParticlesIntegrals(nu,mu,lambda,sigma) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
                twoParticlesIntegrals(mu,nu,sigma,lambda) = twoParticlesIntegrals(mu,nu,lambda,sigma) 
                !            twoParticlesIntegrals(lambda,sigma,mu,nu) = twoParticlesIntegrals(mu,nu,lambda,sigma) 

             end do
          end do
       end do
    end do

    !!    print *, "this 0", this%lowerOccupiedOrbital
    !!    print *, "this 0", this%upperOccupiedOrbital
    !!    print *, "this 0", this%lowerVirtualOrbital 
    !!    print *, "this 0", this%upperVirtualOrbital 

    !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
    open(unit=CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE, file=trim(this%prefixOfFile)//"moint.dat", &
         status='replace',access='sequential', form='unformatted' )

    m = 0
    do p = 1, this%numberOfContractions
       do q = p, this%numberOfContractions
          do r = 1, this%otherNumberOfContractions
             do s = r, this%otherNumberOfContractions

                !            if ( q >= this%lowerVirtualOrbital .and. s >= this%lowerVirtualOrbital .and. &
                !                 p <= this%upperOccupiedOrbital .and. r <= this%upperOccupiedOrbital ) then

                auxTransformedTwoParticlesIntegral = 0
                do mu = 1, this%numberOfContractions
                   do nu = 1, this%numberOfContractions
                      do lambda = 1, this%otherNumberOfContractions
                         do sigma = 1, this%otherNumberOfContractions

                            auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                                 coefficientsOfAtomicOrbitals%values( mu, p )* &
                                 coefficientsOfAtomicOrbitals%values( nu, q )* &
                                 otherCoefficientsOfAtomicOrbitals%values( lambda, r )* &
                                 otherCoefficientsOfAtomicOrbitals%values( sigma, s )* &
                                 twoParticlesIntegrals(mu, nu, lambda, sigma) 

                         end do
                      end do
                   end do
                end do

                write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral

                !end if

             end do
          end do
       end do
    end do

    write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0  



    close(CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE)


  end subroutine TransformIntegralsB_atomicToMolecularOfTwoSpecies

  !>
  !! @brief Escribe los coefficientes de combinacion para los orbitales atomicos.
  !!    El almacenamiento requiere guardar columnas completas una tras de otra
  !<
  subroutine TransformIntegralsB_writeCoefficients( this, coefficients, otherCoefficients )
    implicit none
    type(TransformIntegralsB) :: this
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


  end subroutine TransformIntegralsB_writeCoefficients

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine TransformIntegralsB_exception( typeMessage, description, debugDescription)
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

  end subroutine TransformIntegralsB_exception

end module TransformIntegralsB_