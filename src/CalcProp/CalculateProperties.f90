
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
  !!   - <tt> 2017-01-19 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# Trae de vuelta el calculo de dipolo, incluye la matriz de densidad CI y reordena el modulo
  !<

  type, public :: CalculateProperties
     type(Matrix), allocatable :: overlapMatrix(:)   !!! JORGE
     type(Matrix), allocatable :: densityMatrix(:)
     type(Matrix), allocatable :: momentMatrices(:,:)
     ! character(30) :: name
     ! type(Matrix) :: contributionsOfdipoleMoment
     ! type(Matrix) :: expectedPositions
     ! type(Vector) :: expectedR2
     ! type(Matrix) :: polarizabilityTensor
     ! type(Matrix) :: hyperPolarizabilityTensor(3)
     ! type(Matrix) :: interparticleDistances
     ! type(Matrix) :: interparticleDistancesErrors
     ! type(Matrix) :: interparticleOverlap
     ! type(Vector) :: volume
     ! type(Vector) :: cumulativeDensity
     ! type(Cube), allocatable :: densityCube(:)
     ! type(Cube), allocatable :: orbitalCube(:)
     ! type(Matrix) :: negativeFukui
     ! type(Matrix) :: positiveFukui
  end type CalculateProperties


  integer, parameter, public :: MULLIKEN  =  1
  integer, parameter, public :: LOWDIN    =  2



  !private :: &

  public :: &
       CalculateProperties_constructor, &
       CalculateProperties_destructor, &
       CalculateProperties_showExpectedPositions, &
       CalculateProperties_getExpectedPosition, &
       CalculateProperties_showPopulationAnalyses, &
       CalculateProperties_getPopulation, &
       CalculateProperties_showContributionsToElectrostaticMoment, &
       CalculateProperties_getDipoleOfPuntualCharges, &
       CalculateProperties_getDipoleOfQuantumSpecie
  !     CalculateProperties_expectedR2, &
  !    CalculateProperties_polarizability, &
  !     CalculateProperties_showExpectedR2, &
  !    CalculateProperties_showPolarizabilityTensor, &
  !    CalculateProperties_interparticleDistance,  &
  !    CalculateProperties_interparticleOverlap, &
  !    CalculateProperties_distanceToPoint, &
  !    CalculateProperties_buildDensityCubesLimits, &
  !    CalculateProperties_buildDensityCubes, &
  !    CalculateProperties_volumes, &
  !    CalculateProperties_getPartialCharges, &
  !    CalculateProperties_showIonizationPotentials, &
  !    CalculateProperties_showCharges
  !    CalculateProperties_showVolumes, &
  !    CalculateProperties_getFukuiAt

contains

  !<
  !! @brief Constructor para la clase
  !>
  subroutine CalculateProperties_constructor( this )
    implicit none
    type(CalculateProperties) :: this
    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20)
    character(50) ::  integralsFile
    integer ::  integralsUnit
    character(50) :: occupationsFile, auxstring
    integer :: occupationsUnit
    integer :: numberOfSpecies, speciesID, numberOfContractions

    integralsFile = "lowdin.opints"
    integralsUnit = 30

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate(this%overlapMatrix(numberOfSpecies))
    allocate(this%densityMatrix(numberOfSpecies))
    allocate(this%momentMatrices(numberOfSpecies,3))

    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted") 

    do speciesID=1, numberOfSpecies
       numberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )

       ! Check if there are CI density matrices and read those or the HF matrix
       if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE"  ) then
          print *, "We are calculating properties for ", trim(MolecularSystem_getNameOfSpecie(speciesID)), &
               " in the CI ground state"

          occupationsUnit = 29
          occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"CIOccupations.occ"

          open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")

          auxstring="1" !ground state
          arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
          arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 

          this%densityMatrix(speciesID)= Matrix_getFromFile(unit=occupationsUnit, rows= int(numberOfcontractions,4), &
               columns= int(numberOfcontractions,4), binary=.false., arguments=arguments(1:2))

          close(occupationsUnit)     

       else

          print *, "We are calculating properties for ", trim(MolecularSystem_getNameOfSpecie(speciesID)), &
               " in the HF/KS ground state"

          !! Read density matrix
          arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
          arguments(1) = "DENSITY"

          this%densityMatrix(speciesID) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
               columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

       end if

       ! Overlap matrix
       arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
       arguments(1) = "OVERLAP"

       this%overlapMatrix(speciesID) = Matrix_getFromFile(unit=integralsUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

       !! Load moment Matrices
       arguments(1) = "MOMENTX"    
       this%momentMatrices(speciesID,1) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
            unit=integralsUnit, binary=.true., arguments=arguments(1:2))

       arguments(1) = "MOMENTY"    
       this%momentMatrices(speciesID,2) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
            unit=integralsUnit, binary=.true., arguments=arguments(1:2))

       arguments(1) = "MOMENTZ"    
       this%momentMatrices(speciesID,3) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
            unit=integralsUnit, binary=.true., arguments=arguments(1:2))

    end do

  end subroutine CalculateProperties_constructor

  !<
  !! @brief Destructor para la clase
  !>
  subroutine CalculateProperties_destructor( this )
    implicit none
    type(CalculateProperties) :: this
    integer :: numberOfSpecies, speciesID

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    do speciesID=1, numberOfSpecies
       call Matrix_destructor(this%densityMatrix(speciesID) )
       call Matrix_destructor(this%overlapMatrix(speciesID) )
       call Matrix_destructor(this%momentMatrices(speciesID,1) )
       call Matrix_destructor(this%momentMatrices(speciesID,2) )
       call Matrix_destructor(this%momentMatrices(speciesID,3) )
    end do

    deallocate(this%overlapMatrix)
    deallocate(this%densityMatrix)
    deallocate(this%momentMatrices)


  end subroutine CalculateProperties_destructor

  subroutine CalculateProperties_showPopulationAnalyses(this)
    implicit none
    type (CalculateProperties) :: this ! por medio de este this accedo a todo lo que este en la estructura o type
    ! calculate properties 

    real(8) :: total
    character(10) :: specieName
    integer :: specieID

    !Felix: Vamos a hacer el analisis de poblaciones para todas las especies

    do specieID = 1, MolecularSystem_getNumberOfQuantumSpecies()
       specieName = trim(MolecularSystem_getNameOfSpecie( specieID ))

       !!Obtiene Poblaciones de Mulliken
       print *,""
       print *," POPULATION ANALYSES: "
       print *,"===================="
       print *,""
       print *, " Mulliken Population: for ", specieName
       print *,"---------------------"
       print *,""
       call Vector_show( CalculateProperties_getPopulation(this, "MULLIKEN", specieID, total), &
            flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )

       write (6,"(T25,A10)") "__________"
       write (6,"(T10,A15,F10.6)") "Total = ", total
       print *,""
       print *,"...end of Mulliken Population"
       print *,""
       print *, " Lowdin Population: for ", specieName
       print *,"---------------------"
       print *,""
       call Vector_show( CalculateProperties_getPopulation( this, "LOWDIN", specieID, total),&
            flags = VERTICAL+WITH_KEYS, keys=MolecularSystem_getlabelsofcontractions( specieID ) )
       write (6,"(T25,A10)") "__________"
       write (6,"(T10,A15,F10.6)") "Total = ", total
       print *,""
       print *,"...end of Lowdin Population"
       print *,""
       print *,"END POPULATION ANALYSES "
       print *,""
    end do

    ! specieID=1
    !! Antes: Recorre las especies buscando electrones

    ! search_specie: do i = 1, MolecularSystem_getNumberOfQuantumSpecies()
    !   specieName=""
    !   specieName = trim(MolecularSystem_getNameOfSpecie(i))

    !   if( scan(trim(specieName),"E")==1 ) then
    !     if( scan(trim(specieName),"-")>1 ) then
    !       showPopulations=.true.
    !       specieID=i
    !       exit search_specie
    !     end if
    !   else
    !     showPopulations=.false.
    !   end if

    ! end do search_specie

  end subroutine CalculateProperties_showPopulationAnalyses


  !<
  !! @brief Retorna la poblacion de Mulliken o Lowdin del sistema molecular
  !>
  function CalculateProperties_getPopulation( this, typeOfPopulation, specieID, totalSum, fukuiType )  result( output )
    implicit none
    type (CalculateProperties) :: this 
    character(*) :: typeOfPopulation
    integer  :: specieID
    real(8), optional, intent(out) :: totalSum
    character(*), optional :: fukuiType
    type(Vector) :: output

    type(Matrix) :: auxMatrix
    type(Matrix) :: auxMatrixB
    integer :: numberOfcontractions
    integer :: i

    numberOfcontractions=MolecularSystem_getTotalNumberOfContractions (specieID )

    call Matrix_constructor( auxMatrix, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )
    call Vector_constructor( output, numberOfcontractions   )

    select case( typeOfPopulation )

    case("MULLIKEN")
       auxMatrix%values = matmul(this%densityMatrix(specieID)%values, this%overlapMatrix(specieID)%values )

    case ("LOWDIN")

       auxMatrix%values = matmul(this%densityMatrix(specieID)%values, this%overlapMatrix(specieID)%values )

       auxMatrix = Matrix_pow( this%overlapMatrix(specieID), 0.5_8 )
       auxMatrixB = auxMatrix
       auxMatrix%values = matmul( matmul( auxMatrixB%values , this%densityMatrix(specieID)%values), auxMatrixB%values )

    case default

    end select

    do i=1, numberOfcontractions
       output%values(i) = auxMatrix%values(i,i)
    end do

    !print*,"auxMatrix%values", auxMatrix%values      
    if ( present( totalSum ) ) totalSum = sum(output%values)

    call Matrix_destructor(auxMatrix)
    call Matrix_destructor(auxMatrixB)

  end function CalculateProperties_getPopulation

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
       write (6,"(T5,A15,3F9.4)") trim(MolecularSystem_getNameOfSpecie( i )), CalculateProperties_getExpectedPosition(this, i)
    end do
    print *,""
    print *,"END EXPECTED POSITIONS"
    print *,""
  end subroutine CalculateProperties_showExpectedPositions

  function CalculateProperties_getExpectedPosition( this , specieID) result(output)
    implicit none
    type(CalculateProperties) :: this
    integer :: specieID
    real(8) :: output(3)

    !! Open file for wavefunction                                                                                                          
    output=0.0_8

    output(1)=sum( this%densityMatrix(specieID)%values * this%momentMatrices(specieID,1)%values ) * 0.52917720859
    output(2)=sum( this%densityMatrix(specieID)%values * this%momentMatrices(specieID,2)%values ) * 0.52917720859
    output(3)=sum( this%densityMatrix(specieID)%values * this%momentMatrices(specieID,3)%values ) * 0.52917720859

  end function CalculateProperties_getExpectedPosition

  !<
  !! @brief Muestra las contrinuciones al dipolo de cada especie
  !>
  subroutine CalculateProperties_showContributionsToElectrostaticMoment(this)
    implicit none
    type(CalculateProperties) :: this

    integer :: i, numberOfSpecies
    real(8), allocatable :: dipole(:,:)
    real(8) :: totalDipole(3)

    totalDipole=0.0_8
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

    allocate(dipole(numberOfSpecies+1,3))

    print *,""
    print *," ELECTROSTATIC MOMENTS:"
    print *,"======================"
    print *,""
    print *,"DIPOLE: (DEBYE)"
    print *,"------"
    print *,""
    write (6,"(T19,4A9)") "<Dx>","<Dy>", "<Dz>"," |D|"

    do i=1, numberOfSpecies
       dipole(i,:)=CalculateProperties_getDipoleOfQuantumSpecie(this, i)*2.54174619
       totalDipole(:)=totalDipole(:)+dipole(i,:)
       write (6,"(T5,A15,3F9.4)") trim(MolecularSystem_getNameOfSpecie( i )), dipole(i,:)
    end do

    dipole(numberOfSpecies+1,:)=CalculateProperties_getDipoleOfPuntualCharges()*2.54174619
    totalDipole(:)=totalDipole(:)+dipole(numberOfSpecies+1,:)
    write (6,"(T5,A15,3F9.4)") "Point charges: ", dipole(numberOfSpecies+1,:)

    write (6,"(T22,A28)") "___________________________________"

    write (6,"(T5,A15,3F9.4, F9.4)") "Total ", totalDipole(:), sqrt(sum(totalDipole(:)**2.0 ) )

    print *,""
    print *,"END ELECTROSTATIC MOMENTS"
    print *,""

    deallocate(dipole)

  end subroutine CalculateProperties_showContributionsToElectrostaticMoment

  ! !<
  ! !! @brief Calcula el aporte al dipolo de las cargas puntuales presentes
  ! !>
  function CalculateProperties_getDipoleOfPuntualCharges() result( output )
    implicit none
    real(8) :: output(3)
    integer :: i

    output = 0.0_8

    
    do i=1, size( MolecularSystem_instance%pointCharges )      
       output(:) = output(:) + MolecularSystem_instance%pointCharges(i)%origin(:) * MolecularSystem_instance%pointCharges(i)%charge
    end do

    
  end function CalculateProperties_getDipoleOfPuntualCharges


  !<
  !! @brief Calcula el aporte al dipolo debido a particulas no fijas
  !>
  function CalculateProperties_getDipoleOfQuantumSpecie( this, i ) result( output )
    implicit none
    type(CalculateProperties) :: this
    integer :: i !specieID
    real(8) :: output(3)

    output(1) =sum( this%densityMatrix(i)%values * this%momentMatrices(i,1)%values )
    output(2) =sum( this%densityMatrix(i)%values * this%momentMatrices(i,2)%values )
    output(3) =sum( this%densityMatrix(i)%values * this%momentMatrices(i,3)%values )

    output = output * MolecularSystem_getCharge( i )

  end function CalculateProperties_getDipoleOfQuantumSpecie






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

