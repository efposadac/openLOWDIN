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
     
     integer :: p_lowerOrbital, p_upperOrbital
     integer :: q_lowerOrbital, q_upperOrbital
     integer :: r_lowerOrbital, r_upperOrbital
     integer :: s_lowerOrbital, s_upperOrbital

  end type TransformIntegralsC

  !! TypeOfIntegrals {
  integer, parameter :: ONE_SPECIE			= 0
  integer, parameter :: TWO_SPECIES			= 1
  !! }

  public :: &
       TransformIntegralsC_constructor, &
       TransformIntegralsC_destructor, &
       TransformIntegralsC_show, &
       TransformIntegralsC_atomicToMolecularOfOneSpecie, &
       TransformIntegralsC_atomicToMolecularOfTwoSpecies
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
    integer :: nonZeroIntegrals

!!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer, allocatable :: indexTwoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: auxtempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer*2 :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m, n, u, mm
    integer :: ssize

    ! Reads the number of cores

    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    call TransformIntegralsC_getNumberOfNonZeroRepulsionIntegrals( specieID, nproc, nonZeroIntegrals )
 
    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"

    if ( .not.CONTROL_instance%OPTIMIZE ) then
       call cpu_time(initialTime)
    end if

    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    ssize = this%numberOfContractions
    ssize = (ssize * (ssize + 1))/2 + ssize
    ssize = (ssize * (ssize + 1))/2 + ssize

    this%specieID = specieID

    call TransformIntegralsC_checkMOIntegralType(specieID, this)

    if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    allocate (twoParticlesIntegrals ( nonZeroIntegrals ) )

    twoParticlesIntegrals = 0
    if ( allocated (indexTwoParticlesIntegrals)) deallocate (indexTwoParticlesIntegrals )
    allocate (indexTwoParticlesIntegrals ( ssize ) )

    indexTwoParticlesIntegrals = 0

    m = 0

    !! Read integrals
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
            m = m + 1
            twoParticlesIntegrals(m) = shellIntegrals(i)
            indexTwoParticlesIntegrals(IndexMap_tensorR4ToVectorC(int(aa(i),4),int(bb(i),4),int(cc(i),4),int(dd(i),4), &
                                       this%numberOfContractions )) = m

         end do

       end do loadintegrals

       close (unit)
    end do 

    !! allocate some auxiliary arrays
    if ( allocated (tempA)) deallocate (tempA )
    allocate (tempA ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )
    tempA = 0

    if ( allocated (auxtempA)) deallocate (auxtempA )
    allocate (auxtempA ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )
    auxtempA = 0


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
     
    call omp_set_num_threads(omp_get_max_threads())

    !! begin transformation
    mm = 0
    do p = this%p_lowerOrbital, this%p_upperOrbital

      n = p
      tempA = 0
      auxTransformedTwoParticlesIntegral = 0

      !! First quarter
      do mu = 1, this%numberOfContractions
        if ( abs(coefficientsOfAtomicOrbitals%values( mu, p )) < 1E-10 ) cycle
        !! auxtemp is the twoparticlesintegrals reduced to a three dimensional array
        call TransformIntegralsC_buildArrayA( twoParticlesIntegrals, mu, indexTwoParticlesIntegrals, &
                                                this%numberOfContractions, auxtempA )
        tempA(:,:,:) = tempA(:,:,:) + coefficientsOfAtomicOrbitals%values( mu, p ) * & 
                                      auxtempA(:,:,:)

      end do

      do q = p, this%q_upperOrbital
        u = q
        tempB = 0

        if ( q < this%q_lowerOrbital ) cycle
        !! second quarter
        do nu = 1, this%numberOfContractions
          if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
          tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                                        tempA(nu,:,:)
        end do

        do r = n, this%r_upperOrbital

           tempC = 0

!           if ( r >  this%upperOccupiedOrbital  ) cycle

           !! third quarter
           do lambda = 1, this%numberOfContractions
             tempC(:) = tempC(:) + coefficientsOfAtomicOrbitals%values( lambda, r )* &
                                            tempB(lambda,:)

           end do
           do s = u, this%s_upperOrbital
             auxTransformedTwoParticlesIntegral = 0

             if ( s < this%s_lowerOrbital ) cycle
             !! fourth quarter
             do sigma = 1, this%numberOfContractions
               auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                                                         coefficientsOfAtomicOrbitals%values( sigma, s )* &
                                                         tempC(sigma)

             end do

             write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral
             mm = mm + 1

           end do
           u = r + 1
         end do
       end do
     end do

     write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0  
     print *, "non zero transformed integrals", mm

  end subroutine TransformIntegralsC_atomicToMolecularOfOneSpecie

  subroutine TransformIntegralsC_buildArrayA( integralArray, i, indexArray, ssize , auxtempA)
  implicit none 
  integer, intent(in) :: ssize
  integer, intent(in) :: i
  real(8), intent(in) :: integralArray(:)
  real(8) :: auxtempA(ssize,ssize,ssize)
  integer, intent(in) :: indexArray(:)
  integer :: auxIndex
  integer :: j,k,l,auxm, n, u
  

  !$OMP PARALLEL DO private(j,k,l,u,auxIndex,auxm) shared(ssize,i,n,indexArray,integralArray)
  do j = 1, ssize
    do k = 1, ssize
       do l = k, ssize
         auxIndex = IndexMap_tensorR4ToVectorC(i,j,k,l,ssize )
         auxm = indexArray(auxIndex)
         auxtempA(j,k,l) = integralArray(auxm)
       end do 
     end do 
   end do 
  !$OMP END PARALLEL DO  

    do l = 1, ssize
       do k = l+1, ssize
         auxtempA(:,k,l) = auxtempA(:,l,k)
       end do 
     end do 

  end subroutine TransformIntegralsC_buildArrayA

  !>
  !! @brief Transforma integrales de repulsion atomicas entre particulas de diferente especie
  !! 		a integrales moleculares.
  !<
!!  subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
!!       otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )
  subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies( this, coefficientsOfAtomicOrbitals, &
                          otherCoefficientsOfAtomicOrbitals, molecularIntegrals, specieID, nameOfSpecie, &
                          otherSpecieID, nameOfOtherSpecie )
    implicit none
    type(TransformIntegralsC) :: this
    type(Matrix) :: coefficientsOfAtomicOrbitals
    type(Matrix) :: otherCoefficientsOfAtomicOrbitals
    type(Matrix) :: molecularIntegrals
    integer :: specieID, otherSpecieID
    character(*) :: nameOfSpecie, nameOfOtherSpecie
    integer :: nproc
    integer :: integralStackSize
    real(8) :: initialTime
    real(8) :: finalTime

    integer :: ifile, i
    integer :: unit
    character(50) :: sfile
    integer :: status
    integer :: nonZeroIntegrals

!!    real(8), allocatable :: twoParticlesIntegrals(:,:,:,:)
    real(8), allocatable :: twoParticlesIntegrals(:)
    integer, allocatable :: indexTwoParticlesIntegrals(:)
    real(8)  auxTransformedTwoParticlesIntegral

    real(8), allocatable :: tempA(:,:,:)
    real(8), allocatable :: auxtempA(:,:,:)
    real(8), allocatable :: tempB(:,:)
    real(8), allocatable :: tempC(:)

    integer*2 :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: cc(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: dd(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    integer :: p, q, r, s, mu, nu, lambda, sigma, m, n, u, mm
    integer :: ssize

    ! Reads the number of cores

    nproc = CONTROL_instance%NUMBER_OF_CORES
    integralStackSize = CONTROL_instance%INTEGRAL_STACK_SIZE

    call TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals( specieID, otherSpecieID, nonZeroIntegrals )
    print *, "non zero ints", nonZeroIntegrals
 
    this%prefixOfFile =""//trim(nameOfSpecie)
    this%fileForCoefficients =""//trim(nameOfSpecie)//"mo.values"


    this%numberOfContractions=size(coefficientsOfAtomicOrbitals%values,dim=1)
    ssize = this%numberOfContractions
    ssize = (ssize * (ssize + 1))/2 + ssize
    ssize = (ssize * (ssize + 1))/2 + ssize

    this%specieID = specieID

    call TransformIntegralsC_checkMOIntegralType(specieID, this)

    if ( allocated (twoParticlesIntegrals)) deallocate (twoParticlesIntegrals )
    allocate (twoParticlesIntegrals ( nonZeroIntegrals ) )

    twoParticlesIntegrals = 0
    if ( allocated (indexTwoParticlesIntegrals)) deallocate (indexTwoParticlesIntegrals )
    allocate (indexTwoParticlesIntegrals ( ssize ) )

    indexTwoParticlesIntegrals = 0

    m = 0

    !! Read integrals
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
            m = m + 1
            twoParticlesIntegrals(m) = shellIntegrals(i)
            indexTwoParticlesIntegrals(IndexMap_tensorR4ToVectorC(int(aa(i),4),int(bb(i),4),int(cc(i),4),int(dd(i),4), &
                                       this%numberOfContractions )) = m

         end do

       end do loadintegrals

       close (unit)
    end do 

    !! allocate some auxiliary arrays
    if ( allocated (tempA)) deallocate (tempA )
    allocate (tempA ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )
    tempA = 0

    if ( allocated (auxtempA)) deallocate (auxtempA )
    allocate (auxtempA ( this%numberOfContractions , &
                                      this%numberOfContractions, &
                                      this%numberOfContractions ) )
    auxtempA = 0


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
     
    call omp_set_num_threads(omp_get_max_threads())

    !! begin transformation
    mm = 0
    do p = this%p_lowerOrbital, this%p_upperOrbital

      n = p
      tempA = 0
      auxTransformedTwoParticlesIntegral = 0

      !! First quarter
      do mu = 1, this%numberOfContractions
        if ( abs(coefficientsOfAtomicOrbitals%values( mu, p )) < 1E-10 ) cycle
        !! auxtemp is the twoparticlesintegrals reduced to a three dimensional array
        call TransformIntegralsC_buildArrayA( twoParticlesIntegrals, mu, indexTwoParticlesIntegrals, &
                                                this%numberOfContractions, auxtempA )
        tempA(:,:,:) = tempA(:,:,:) + coefficientsOfAtomicOrbitals%values( mu, p ) * & 
                                      auxtempA(:,:,:)

      end do

      do q = p, this%q_upperOrbital
        u = q
        tempB = 0

        if ( q < this%q_lowerOrbital ) cycle
        !! second quarter
        do nu = 1, this%numberOfContractions
          if ( abs(coefficientsOfAtomicOrbitals%values( nu, q )) < 1E-10 ) cycle
          tempB(:,:) = tempB(:,:) + coefficientsOfAtomicOrbitals%values( nu, q )* &
                                        tempA(nu,:,:)
        end do

        do r = n, this%r_upperOrbital

           tempC = 0

!           if ( r >  this%upperOccupiedOrbital  ) cycle

           !! third quarter
           do lambda = 1, this%numberOfContractions
             tempC(:) = tempC(:) + coefficientsOfAtomicOrbitals%values( lambda, r )* &
                                            tempB(lambda,:)

           end do
           do s = u, this%s_upperOrbital
             auxTransformedTwoParticlesIntegral = 0

             if ( s < this%s_lowerOrbital ) cycle
             !! fourth quarter
             do sigma = 1, this%numberOfContractions
               auxTransformedTwoParticlesIntegral = auxTransformedTwoParticlesIntegral + &
                                                         coefficientsOfAtomicOrbitals%values( sigma, s )* &
                                                         tempC(sigma)

             end do

             write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) p,q,r,s, auxTransformedTwoParticlesIntegral
             mm = mm + 1

           end do
           u = r + 1
         end do
       end do
     end do

     write (CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE) -1,0,0,0, 0  
     print *, "non zero transformed integrals", mm




  end subroutine TransformIntegralsC_atomicToMolecularOfTwoSpecies

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
    this%p_lowerOrbital = 1
    this%p_upperOrbital = totalNumberOfContractions
    this%q_lowerOrbital = 1
    this%q_upperOrbital = totalNumberOfContractions
    this%r_lowerOrbital = 1
    this%r_upperOrbital = totalNumberOfContractions
    this%s_lowerOrbital = 1
    this%s_upperOrbital = totalNumberOfContractions


    !! only the (ia|jb) integrals will be transformed
    if ( CONTROL_instance%MOLLER_PLESSET_CORRECTION == 2 .or. &
         ( CONTROL_instance%PT_ORDER == 2 .and.  CONTROL_instance%IONIZE_MO <= totalOCcupation ) ) then

      this%p_lowerOrbital = 1
      this%p_upperOrbital = totalOccupation
      this%q_lowerOrbital = totalOccupation + 1
      this%q_upperOrbital = totalNumberOfContractions
      this%r_lowerOrbital = 1
      this%r_upperOrbital = totalOccupation
      this%s_lowerOrbital = totalOccupation + 1
      this%s_upperOrbital = totalNumberOfContractions

    end if


  end subroutine TransformIntegralsC_checkMOIntegralType

  subroutine TransformIntegralsC_getNumberOfNonZeroRepulsionIntegrals( specieID, nproc, nonZeroIntegrals )
    implicit none
    integer :: specieID, nproc, auxNonZeroIntegrals
    integer :: nonZeroIntegrals
    integer :: ifile, unit
    character(50) :: sfile
    character(30) :: nameOfSpecie

    nonZeroIntegrals = 0

    do ifile = 1, nproc

       write(sfile,*) ifile
       sfile = trim(adjustl(sfile))
       unit = ifile+50

       nameOfSpecie = MolecularSystem_getNameOfSpecie( specieID )          

       !! open file (order, integral(shell))
       open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".nints", status='old',access='sequential', form='Unformatted')

       read (unit) auxNonZeroIntegrals
       nonZeroIntegrals = nonZeroIntegrals + auxNonZeroIntegrals
       close (unit)
     end do

  end subroutine TransformIntegralsC_getNumberOfNonZeroRepulsionIntegrals

  subroutine TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals( i, j, nonZeroIntegrals )
    implicit none
    integer :: i, j, auxNonZeroIntegrals
    integer :: nonZeroIntegrals
    integer :: unit
    character(30) :: nameOfSpecie, nameOfOtherSpecie

    nonZeroIntegrals = 0

    nameOfSpecie = MolecularSystem_getNameOfSpecie( i )          
    nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( j )          
    !! open file for integrals
    open(UNIT=54,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".nints", &
                 STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

    read (54) nonZeroIntegrals
    close (54)

  end subroutine TransformIntegralsC_getNumberOfNonZeroCouplingIntegrals


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
