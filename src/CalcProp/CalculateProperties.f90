
!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!  Prof. G. MERINO's Lab. Universidad de Guanajuato
!!    http://quimera.ugto.mx/qtc/gmerino.html
!!
!!  Authors:
!!    E. F. Posada (efposadac@unal.edu.co)
!!
!!  Contributors:
!!
!!    Todos los derechos reservados, 2011
!!
!!******************************************************************************

module CalculateProperties_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use Units_
  use Exception_
  use WaveFunction_
  use ContractedGaussian_
  implicit none

  !>
  !!
  !!  Este modulo define una seudoclase para calculo de propiedades derivadas de
  !! la funcion de onda como cargas, dipolos, polarizabilidades, etc.
  !!
  !! @author Sergio A. Gonzalez Monico
  !!
  !! <b> Fecha de creacion : </b> 2007-09-18
  !!
  !! <b> Historial de modificaciones: </b>
  !!
  !!   - <tt> 2007-09-18 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
  !!        -# Creacion de modulo y metodos basicos.
  !!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
  !!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
  !!   - <tt> 2011-11-23 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# Adds numerical integration properties, ADPT calculations and brings population analyses 
  !!   - <tt> 2014-01-23 </tt>: Matheus Rodriguez ( matrodriguezalv@unal.edu.co )
        !!        -# Reescribe y adapta el modulo de Calculate properties en Lowdin2
  !<

  type, public :: CalculateProperties
        character(30) :: name
       ! type(Matrix) :: contributionsOfdipoleMoment
        type(Matrix) :: expectedPositions
        type(Vector) :: expectedR2
       ! type(Matrix) :: polarizabilityTensor
       ! type(Matrix) :: hyperPolarizabilityTensor(3)
       ! type(Matrix) :: interparticleDistances
       ! type(Matrix) :: interparticleDistancesErrors
        type(Matrix) :: interparticleOverlap
       ! type(Vector) :: volume
       ! type(Vector) :: cumulativeDensity
       ! type(Cube), allocatable :: densityCube(:)
       ! type(Cube), allocatable :: orbitalCube(:)
        type(Matrix) :: negativeFukui
        type(Matrix) :: positiveFukui
        type(Matrix) :: overlapMatrix   !!! JORGE
        type(Matrix) :: densityMatrix
  end type


       integer, parameter, public :: MULLIKEN  =  1
       integer, parameter, public :: LOWDIN    =  2



  !private :: &
    !CalculateProperties_getDipoleOfPuntualCharges
    ! CalculateProperties_getDipoleOfQuantumSpecie

  public :: &
!    CalculateProperties_constructor, &
!    CalculateProperties_destructor, &
!    CalculateProperties_dipole, &
     CalculateProperties_expectedPosition, &
!     CalculateProperties_expectedR2, &
!    CalculateProperties_polarizability, &
!    CalculateProperties_showContributionsToElectrostaticMoment, &
     CalculateProperties_showExpectedPositions, &
!     CalculateProperties_showExpectedR2, &
!    CalculateProperties_showPolarizabilityTensor, &
!    CalculateProperties_interparticleDistance,  &
!    CalculateProperties_interparticleOverlap, &
!    CalculateProperties_distanceToPoint, &
!    CalculateProperties_buildDensityCubesLimits, &
!    CalculateProperties_buildDensityCubes, &
!    CalculateProperties_volumes, &
     CalculateProperties_showPopulationAnalyses, &
     CalculateProperties_getPopulation
!    CalculateProperties_getPartialCharges, &
!    CalculateProperties_showIonizationPotentials, &
!    CalculateProperties_showCharges
!    CalculateProperties_showVolumes, &
!    CalculateProperties_getFukuiAt
    
contains



  subroutine CalculateProperties_showPopulationAnalyses()
    implicit none
               ! type (CalculateProperties) :: this ! por medio de este this accedo a todo lo que este en la estructura o type
                                                   ! calculate properties 

    real(8) :: total
    character(10) :: specieName
    integer :: i,specieID
    logical :: showPopulations 


    specieID=1
    !! Recorre las especies buscando electrones
                                        
    search_specie: do i = 1, MolecularSystem_getNumberOfQuantumSpecies()
      specieName=""
      specieName = trim(MolecularSystem_getNameOfSpecie(i))

      if( scan(trim(specieName),"E")==1 ) then
        if( scan(trim(specieName),"-")>1 ) then
          showPopulations=.true.
          specieID=i
          exit search_specie
        end if
      else
        showPopulations=.false.
      end if

    end do search_specie

    if( showPopulations ) then
      !!Obtiene Poblaciones de Mulliken
      print *,""
      print *," POPULATION ANALYSES: "
      print *,"===================="
      print *,""
      print *, " Mulliken Population: "
      print *,"---------------------"
      print *,""
      call Vector_show( CalculateProperties_getPopulation(MULLIKEN, total,trim(specieName)),&
        flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )

      write (6,"(T25,A10)") "__________"
      write (6,"(T10,A15,F10.6)") "Total = ", total
      print *,""
      print *,"...end of Mulliken Population"
      print *,""
      print *, " Lowdin Population:"
      print *,"---------------------"
      print *,""
      call Vector_show( CalculateProperties_getPopulation( LOWDIN, total,trim(specieName)),&
        flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )
      write (6,"(T25,A10)") "__________"
      write (6,"(T10,A15,F10.6)") "Total = ", total
      print *,""
      print *,"...end of Lowdin Population"
      print *,""
      print *,"END POPULATION ANALYSES "
      print *,""
    end if

  end subroutine CalculateProperties_showPopulationAnalyses


 !<
 !! @brief Retorna la poblacion de Mulliken o Lowdin del sistema molecular
 !>
 function CalculateProperties_getPopulation( typeOfPopulation, totalSum, nameOfSpecie, fukuiType )  result( output )
   implicit none
  ! type (CalculateProperties) :: this 
   integer :: typeOfPopulation
   real(8), optional, intent(out) :: totalSum
   character(*),optional  :: nameOfSpecie
   character(*), optional :: fukuiType
   type(Vector) :: output
   
   type(Matrix) :: densityMatrix
   type(Matrix) :: overlapMatrix
   type(Matrix) :: auxMatrix
   type(Matrix) :: auxMatrixB
   character(10) :: auxNameOfSpecie
   integer :: numberOfcontractions
   integer :: speciesID
   integer :: i
   character(50) :: wfnFile
   integer :: wfnUnit
   character(50) :: arguments(20)
   character(50) ::  integralsFile
   integer ::  integralsUnit


  integralsFile = "lowdin.opints"
  integralsUnit = 30


  wfnFile = "lowdin.wfn"
  wfnUnit = 20


  !! Open file for wavefunction                                                                                                                                             
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

  open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted") 


  auxNameOfSpecie="E-" 
  if (present( nameOfSpecie ) )  then
     auxNameOfSpecie = trim(nameOfSpecie)
  end if

  !    if ( MolecularSystem_isSet() ) then
  if ( .not. present( fukuiType) .or. (present(fukuiType) .and. trim(auxNameofSpecie) .eq. "E-") ) then
     speciesID =MolecularSystem_getSpecieID (  nameOfSpecie = trim(auxNameOfSpecie) )
     numberOfcontractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )
     call Matrix_constructor( auxMatrix, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )
     call Vector_constructor( output, numberOfcontractions   )


     arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
     arguments(1) = "DENSITY"
     !  WaveFunction_instance(speciesID)%densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

!     print *, "CalculateProperties_getPopulation 0"

     densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
!!! Solo para esta subrutina

     arguments(1) = "OVERLAP"
 

     !!! Abrir el archivo lowdin.opints

     !  WaveFunction_instance(speciesID)%overlapMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))


!     print *, "CalculateProperties_getPopulation 1"

     overlapMatrix = Matrix_getFromFile(unit=integralsUnit, rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
  !!! General


!     print *, "CalculateProperties_getPopulation 2"

     select case( typeOfPopulation )

     case( MULLIKEN )
        ! Reto por que cambia bastante en lowdin 2
        ! Leer el modulo wave function de HF en lowdin 2
        ! Leer todo el directorio HF de lowdin 2    
    !!    overlapMatrix = WaveFunction_instance( speciesID )%overlapMatrix
        !  if (trim(fukuiType)=="positive" .and. trim(auxNameOfSpecie) == "E-") then
        !    densityMatrix = this%positiveFukui
        !  else if  (trim(fukuiType)=="negative" .and. trim(auxNameOfSpecie) == "E-") then
        !     densityMatrix = this%negativeFukui
        ! else              
    !!    densityMatrix =  WaveFunction_instance( speciesID )%densityMatrix      
        ! end if
        auxMatrix%values = matmul(densityMatrix%values, overlapMatrix%values )
        
     case (LOWDIN)
        
     !!   overlapMatrix = WaveFunction_instance( speciesID )%overlapMatrix
     !!   densityMatrix =  WaveFunction_instance( speciesID )%densityMatrix          
        auxMatrix%values = matmul(densityMatrix%values, overlapMatrix%values )
        
        auxMatrix = Matrix_pow( overlapMatrix, 0.5_8 )
        auxMatrixB = auxMatrix
        auxMatrix%values = matmul( matmul( auxMatrixB%values , densityMatrix%values), auxMatrixB%values )
        
        call Matrix_destructor(auxMatrixB)
        
     case default
        
     end select

     
     
     do i=1, numberOfcontractions
        output%values(i) = auxMatrix%values(i,i)
     end do
                     
     !print*,"auxMatrix%values", auxMatrix%values      
     if ( present( totalSum ) ) totalSum = sum(output%values)

     call Matrix_destructor(overlapMatrix)
     call Matrix_destructor(auxMatrix)
     call Matrix_destructor(auxMatrixB)
     call Matrix_destructor(densityMatrix)
     
  end if
                       
!else

  ! call CalculateProperties_exception(ERROR, "You should set the molecular system before use this function", &
  !      "Class object CalculateProperties in the getPopulation function" )
   
   close(wfnUnit)
   close(integralsUnit)
!end if

  end function CalculateProperties_getPopulation

!  !<
!  !! @brief  Calculates the expected position for each quantum specie
!  !>
!  subroutine CalculateProperties_expectedR2( this )
!    implicit none
!    type(CalculateProperties) :: this
!    type(Matrix) :: densityMatrix
!    type(Matrix) :: R2Matrix
!    type(ContractedGaussian):: dxx
!    type(ContractedGaussian):: dyy
!    type(ContractedGaussian):: dzz
!    type(ExternalPotential) :: R2Operator(1)
!    real(8) :: expo(1)
!    real(8) :: coefficient(1)
!    real(8) :: orig(3)
!    integer(8) ::  angMom
!    integer(8) ::  angMomIndex(3)
!    character(30) :: nameOfSpecieSelected
!    integer :: i
!                integer :: numberOfContractions
!    integer(8) :: numberOfSpecies
!
!    numberOfSpecies = Particle_Manager_getNumberOfQuantumSpecies()
!    call Vector_constructor(this%expectedR2,int(numberOfSpecies,4))
!
!    !! Preparing the R2 operator (Treated as a potential)
!
!    expo(1)=0.0
!    coefficient(1)=1.0
!    orig=0.0
!    angMom=2
!    angMomIndex(1)=2.0
!    angMomIndex(2)=0.0
!    angMomIndex(3)=0.0
!
!    call ContractedGaussian_constructor( dxx , orbitalsExponents=expo , &
!    oefficients=coefficient , origin=orig , angularMoment=angMom, angularMomentIndex=angMomIndex, noNormalize=.true. )
!
!    angMomIndex(1)=0.0
!    angMomIndex(2)=2.0
!    angMomIndex(3)=0.0
!
!    call ContractedGaussian_constructor( dyy , orbitalsExponents=expo , &
!    oefficients=coefficient , origin=orig , angularMoment=angMom, angularMomentIndex=angMomIndex, noNormalize=.true. )
!
!    angMomIndex(1)=0.0
!    angMomIndex(2)=0.0
!    angMomIndex(3)=2.0
!
!    call ContractedGaussian_constructor( dzz , orbitalsExponents=expo , &
!    contractionCoefficients=coefficient , origin=orig , angularMoment=angMom, angularMomentIndex=angMomIndex, noNormalize=.true. )
!
!
!    do i=1, numberOfSpecies
!      nameOfSpecieSelected = trim( Particle_Manager_getNameOfSpecie( i ) )
!      numberOfContractions = Particle_Manager_getTotalNumberOfContractions( i )
!      call Matrix_constructor (densityMatrix, int(numberOfContractions,8), int(numberOfContractions,8))
!      densityMatrix = MolecularSystem_getDensityMatrix( trim(nameOfSpecieSelected) )
!      call ExternalPotential_constructor(R2Operator(1), "R2", nameOfSpecieSelected)
!      R2Operator(1)%numOfComponents=3
!      allocate(R2Operator(1)%gaussianComponents(R2Operator(1)%numOfComponents))
!      R2Operator(1)%gaussianComponents(1)=dxx
!      R2Operator(1)%gaussianComponents(2)=dyy
!      R2Operator(1)%gaussianComponents(3)=dzz
!      R2Matrix=IntegralManager_getInteractionWithPotentialMatrix(R2Operator, sspecieID=i )
!      this%expectedR2%values(i)= sum( densityMatrix%values * R2Matrix%values )
!      call Matrix_destructor( densityMatrix )
!      call Matrix_destructor( R2Matrix )
!      call ExternalPotential_destructor( R2Operator(1) )
!    end do
!
!  end subroutine CalculateProperties_expectedR2
!
!  !<
!  !! @brief Muestra las contrinuciones al dipolo de cada especie
!  !>
!  subroutine CalculateProperties_showExpectedR2(this)
!    implicit none
!    type(CalculateProperties) :: this
!    character(30) :: nameOfSpecieSelected
!    integer :: i,j
!    integer :: numberOfSpecies
!    real(8) :: output(3)
!
!                if( externalPotential_Manager_instance%isInstanced ) then
!                   numberOfSpecies = Particle_Manager_getNumberOfQuantumSpecies()
!                   print *,""
!                   print *," EXPECTED <R^2> OF QUANTUM SPECIES:"
!                   print *,"======================"
!                   print *,""
!                   print *,"IN BOHR^2"
!                   print *,"------"
!                   print *,""
!                   write (6,"(T19,A9)") "<R^2>"
!                   do i=1, numberOfSpecies
!                      write (6,"(T5,A15,F9.4)") trim(Particle_Manager_getNameOfSpecie( i )), (this%expectedR2%values(i))
!                   end do
!                   print *,""
!                   print *,"END EXPECTED <R^2>"
!                   print *,""
!                end if
!
!  end subroutine CalculateProperties_showExpectedR2

  subroutine CalculateProperties_expectedPosition( this )
    implicit none
    type(CalculateProperties) :: this
    type(Matrix) :: densityMatrix
    type(Matrix) :: momentMatrix
    character(30) :: nameOfSpecieSelected
    integer :: i
    integer :: numberOfSpecies
    integer :: unit
    integer :: totalNumberOfContractions
    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20)
    character(50) ::  integralsFile
    integer ::  integralsUnit
 
    integralsFile = "lowdin.opints"
    integralsUnit = 30
  
    wfnFile = "lowdin.wfn"
    wfnUnit = 20
  
    !! Open file for wavefunction                                                                                                                                             

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    call Matrix_constructor(this%expectedPositions,int(numberOfSpecies,8),3_8)


    do i=1, numberOfSpecies

      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
      open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted") 

      !! Get number of shells and number of cartesian contractions
      totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )          

      arguments(2) = trim(MolecularSystem_getNameOfSpecie(i))

      !! Load density Matrix
      arguments(1) = "DENSITY"    
      densityMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=wfnUnit, binary=.true., arguments=arguments(1:2))

      !! Load moment Matrix
      arguments(1) = "MOMENTX"    
      momentMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=integralsUnit, binary=.true., arguments=arguments(1:2))
      !! Calcula integrales asociadas al momento electrico dipolar en la direccion x
      this%expectedPositions%values(i,1)=sum( densityMatrix%values * momentMatrix%values ) * 0.52917720859

      arguments(1) = "MOMENTY"    
      momentMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=integralsUnit, binary=.true., arguments=arguments(1:2))
      !! Calcula integrales asociadas al momento electrico dipolar en la direccion y
      this%expectedPositions%values(i,2) =sum( densityMatrix%values * momentMatrix%values ) * 0.52917720859

      arguments(1) = "MOMENTZ"    
      momentMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=integralsUnit, binary=.true., arguments=arguments(1:2))
      !! Calcula integrales asociadas al momento electrico dipolar en la direccion z
      this%expectedPositions%values(i,3) =sum( densityMatrix%values * momentMatrix%values ) * 0.52917720859

      close(wfnUnit)
      close(integralsUnit)

      call Matrix_destructor(momentMatrix)
      call Matrix_destructor(densityMatrix)

   end do

  end subroutine CalculateProperties_expectedPosition

  !<
  !! @brief Muestra las contrinuciones al dipolo de cada especie
  !>
  subroutine CalculateProperties_showExpectedPositions(this)
    implicit none
    type(CalculateProperties) :: this
    character(30) :: nameOfSpecieSelected
    integer :: i,j
    integer :: numberOfSpecies
    real(8) :: output(3)
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
                print *,""
    print *," EXPECTED POSITIONS OF QUANTUM SPECIES:"
    print *,"======================"
    print *,""
    print *,"POSITIONS IN ANGSTROMS"
    print *,"------"
    print *,""
    write (6,"(T19,4A9)") "<x>","<y>", "<z>", ""
    do i=1, numberOfSpecies
       write (6,"(T5,A15,3F9.4)") trim(MolecularSystem_getNameOfSpecie( i )), (this%expectedPositions%values(i,j),j=1,3)
     end do
    print *,""
    print *,"END EXPECTED POSITIONS"
    print *,""
  end subroutine CalculateProperties_showExpectedPositions






  subroutine CalculateProperties_exception( typeMessage, description, debugDescription)
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
  
  end subroutine CalculateProperties_exception

end module CalculateProperties_
