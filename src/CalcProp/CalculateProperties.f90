
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
  use ContractedGaussian_
  use DirectIntegralManager_
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
       CalculateProperties_getDipoleOfQuantumSpecies, &
       CalculateProperties_getExpectedR2, &
       CalculateProperties_showRMSradius
  !    CalculateProperties_polarizability, &
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
  subroutine CalculateProperties_constructor( this,fileName )
    implicit none
    type(CalculateProperties) :: this
    character(*) :: fileName

    character(50) :: wfnFile
    integer :: wfnUnit
    character(50) :: arguments(20)
    character(50) ::  integralsFile
    integer ::  integralsUnit
    character(50) :: occupationsFile, auxstring
    integer :: occupationsUnit
    integer :: numberOfSpecies, speciesID, numberOfContractions
    logical :: existFile, existCIFile

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate(this%densityMatrix(numberOfSpecies))
    allocate(this%overlapMatrix(numberOfSpecies))
    allocate(this%momentMatrices(numberOfSpecies,9))

    !! Open file for HF wavefunction   
    wfnFile = trim(fileName)//".wfn"
    wfnUnit = 20
    inquire(FILE = wfnFile, EXIST = existFile )
    !! Open file for CI wavefunction   
    occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    occupationsUnit = 29
    inquire(FILE = occupationsFile, EXIST = existCIFile )
    if ( existCIFile ) then
       open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")
       do speciesID=1, numberOfSpecies
          numberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )
          print *, "We are calculating properties for ", trim(MolecularSystem_getSymbolOfSpecies(speciesID)), &
               " in the CI ground state"
          auxstring="1" !ground state
          arguments(2) = MolecularSystem_getNameOfSpecies(speciesID)
          arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
          this%densityMatrix(speciesID)= Matrix_getFromFile(unit=occupationsUnit, rows= int(numberOfcontractions,4), &
               columns= int(numberOfcontractions,4), binary=.false., arguments=arguments(1:2))
       end do
       close(occupationsUnit)     
    else if( existFile) then
       open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
       do speciesID=1, numberOfSpecies
          numberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )
          print *, "We are calculating properties for ", trim(MolecularSystem_getSymbolOfSpecies(speciesID)), &
               " in the HF/KS ground state"
          arguments(2) = MolecularSystem_getNameOfSpecies(speciesID)
          arguments(1) = "DENSITY"
          this%densityMatrix(speciesID) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
               columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
       end do
       close(wfnUnit)     
    else
       call Exception_stopError("I did not find the file "//trim(wfnFile)//" or "//trim(occupationsFile), "in CalculatateProperties constructor")
       ! Check if there are CI density matrices and read those or the HF matrix
    end if

    integralsFile = "lowdin.opints"
    integralsUnit = 30
    inquire(FILE = integralsFile, EXIST = existFile )
    if( existFile) then
       open(unit=integralsUnit, file=trim(integralsFile), status="old", form="unformatted") 
       do speciesID=1, numberOfSpecies
          numberOfContractions =  MolecularSystem_getTotalNumberOfContractions (speciesID )
          ! Overlap matrix
          arguments(2) = MolecularSystem_getNameOfSpecies(speciesID)
          arguments(1) = "OVERLAP"
          this%overlapMatrix(speciesID) = Matrix_getFromFile(unit=integralsUnit, rows= int(numberOfContractions,4), &
               columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
          !! Load moment Matrices
          arguments(1) = "MOMENTX0"    
          this%momentMatrices(speciesID,1) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTY0"    
          this%momentMatrices(speciesID,2) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTZ0"    
          this%momentMatrices(speciesID,3) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          !! Load moment Matrices
          arguments(1) = "MOMENTXX"    
          this%momentMatrices(speciesID,4) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTYY"    
          this%momentMatrices(speciesID,5) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTZZ"    
          this%momentMatrices(speciesID,6) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          !! Load moment Matrices
          arguments(1) = "MOMENTXY"    
          this%momentMatrices(speciesID,7) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTXZ"    
          this%momentMatrices(speciesID,8) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))

          arguments(1) = "MOMENTYZ"    
          this%momentMatrices(speciesID,9) = Matrix_getFromFile(rows=numberOfContractions, columns=numberOfContractions, &
               unit=integralsUnit, binary=.true., arguments=arguments(1:2))
       end do
       close(integralsUnit)
    else
       do speciesID=1, numberOfSpecies
          call DirectIntegralManager_getOverlapIntegrals(molecularSystem_instance,speciesID,this%overlapMatrix(speciesID))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,1,this%momentMatrices(speciesID,1))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,2,this%momentMatrices(speciesID,2))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,3))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,4))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,5))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,6))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,7))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,8))
          call DirectIntegralManager_getMomentIntegrals(molecularSystem_instance,speciesID,3,this%momentMatrices(speciesID,9))
       end do
    end if

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
       call Matrix_destructor(this%momentMatrices(speciesID,4) )
       call Matrix_destructor(this%momentMatrices(speciesID,5) )
       call Matrix_destructor(this%momentMatrices(speciesID,6) )
       call Matrix_destructor(this%momentMatrices(speciesID,7) )
       call Matrix_destructor(this%momentMatrices(speciesID,8) )
       call Matrix_destructor(this%momentMatrices(speciesID,9) )
    end do

    deallocate(this%overlapMatrix)
    deallocate(this%densityMatrix)
    deallocate(this%momentMatrices)


  end subroutine CalculateProperties_destructor

  subroutine CalculateProperties_showPopulationAnalyses(this)
    implicit none
    type (CalculateProperties) :: this ! por medio de este this accedo a todo lo que este en la estructura o type
    ! calculate properties 
    real(8) :: total(2), atomSum, atomSpin
    character(10) :: speciesName
    character(10) :: speciesNickname
    character(10) :: analysis(2)
    character(10) :: header(2)
    integer :: speciesID, i, j, k, l, type
    type(Matrix) :: populations

    analysis(1)="MULLIKEN"
    analysis(2)="LOWDIN"

    !Felix: Vamos a hacer el analisis de poblaciones para todas las especies

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

       do type= 1, size(analysis)

          speciesName = trim(MolecularSystem_getNameOfSpecies( speciesID ))

          if(trim(speciesName) .eq. "E-ALPHA") then
             speciesNickname="E-"
             header(1)="total"
             header(2)="spin"
          else if(trim(speciesName) .eq. "E-BETA") then
             cycle
          else
             speciesNickname=trim(MolecularSystem_getSymbolOfSpecies( speciesID ))
             header(1)=""
             header(2)=""
          end if

          populations = CalculateProperties_getPopulation(this, trim(analysis(type) ), speciesID, total)

          !!Obtiene Poblaciones de Mulliken
          print *,""
          print *," POPULATION ANALYSES: "
          print *,"===================="
          print *,""
          print *, trim(analysis(type)), " Population: for ", speciesNickname
          print *,"---------------------"
          print *,""
          call Matrix_show( populations, &
               rowkeys = MolecularSystem_getlabelsofcontractions( speciesID ), &
               columnkeys = header,&
               flags=WITH_BOTH_KEYS)
          if(trim(speciesName) .eq. "E-ALPHA") then
             write (*,"(T28,A25)") "__________     __________"
             write (*,"(T13,A10,F15.6,F15.6)") "Total = ", total(1), total(2)
          else
             write (*,"(T28,A10)") "__________"
             write (*,"(T13,A10,F15.6)") "Total = ", total(1)
          end if
          print *,""
          print *, " Atomic ", trim(analysis(type) ), " Population: for ", speciesNickname
          l=0
          do i=1, size(MolecularSystem_instance%species(speciesID)%particles)
             atomSum=0
             atomSpin=0
             do j=1, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)
                do k=1, MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital
                   l=l+1
                   atomSum=atomSum+populations%values(l,1)
                   if(trim(speciesName) .eq. "E-ALPHA") atomSpin=atomSpin+populations%values(l,2)

                end do
             end do

             if(trim(speciesName) .eq. "E-ALPHA") then
                write(*,"(T13,A10,F15.6,F15.6)") ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner), atomSum, atomSpin
             else
                write(*,"(T13,A10,F15.6)") ParticleManager_getSymbol(MolecularSystem_instance%species(speciesID)%particles(i)%owner), atomSum
             end if
          end do

          print *,""
          print *,"...end of ", trim(analysis(type) )," Population"
          print *,""

          print *,""
          print *,"END POPULATION ANALYSES "
          print *,""
       end do
    end do

    ! specieID=1
    !! Antes: Recorre las especies buscando electrones

    ! search_specie: do i = 1, MolecularSystem_getNumberOfQuantumSpecies()
    !   speciesName=""
    !   speciesName = trim(MolecularSystem_getNameOfSpecies(i))

    !   if( scan(trim(speciesName),"E")==1 ) then
    !     if( scan(trim(speciesName),"-")>1 ) then
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
  function CalculateProperties_getPopulation( this, typeOfPopulation, speciesID, totalSum, fukuiType )  result( output )
    implicit none
    type (CalculateProperties) :: this 
    character(*) :: typeOfPopulation
    integer  :: speciesID
    real(8), optional, intent(out) :: totalSum(2)
    character(*), optional :: fukuiType
    type(Matrix) :: output

    type(Matrix) :: auxMatrix
    type(Matrix) :: auxMatrixB
    type(Matrix) :: otherAuxMatrix
    type(Matrix) :: otherAuxMatrixB
    integer :: numberOfcontractions
    integer :: i
    integer  :: otherSpeciesID
    character(10) :: speciesName

    numberOfcontractions=MolecularSystem_getTotalNumberOfContractions (speciesID )

    call Matrix_constructor( auxMatrix, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )

    speciesName=trim(MolecularSystem_getNameOfSpecies( speciesID ))
    if(trim(speciesName) .eq. "E-ALPHA") then
       call Matrix_constructor( output, int( numberOfcontractions, 8), 2_8 )
       otherSpeciesID=speciesID+1
       call Matrix_constructor( otherAuxMatrix, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )
       call Matrix_constructor( otherAuxMatrixB, int( numberOfcontractions, 8), int( numberOfcontractions, 8) )
    else
       call Matrix_constructor( output, int( numberOfcontractions, 8), 1_8 )
    end if

    select case( typeOfPopulation )

    case("MULLIKEN")

       auxMatrix=Matrix_product_dgemm(this%densityMatrix(speciesID), this%overlapMatrix(speciesID))
       if(trim(speciesName) .eq. "E-ALPHA") otherAuxMatrix = Matrix_product_dgemm(this%densityMatrix(otherSpeciesID), this%overlapMatrix(otherSpeciesID) )

    case ("LOWDIN")

       auxMatrix = Matrix_product_dgemm(this%densityMatrix(speciesID), this%overlapMatrix(speciesID) )

       auxMatrix = Matrix_pow( this%overlapMatrix(speciesID), 0.5_8, method="SVD" )
       auxMatrixB = auxMatrix
       auxMatrix = Matrix_product_dgemm( Matrix_product_dgemm( auxMatrixB , this%densityMatrix(speciesID)), auxMatrixB )

       if(trim(speciesName) .eq. "E-ALPHA") then
          otherAuxMatrix = Matrix_product_dgemm(this%densityMatrix(otherSpeciesID), this%overlapMatrix(otherSpeciesID) )

          otherAuxMatrix = Matrix_pow( this%overlapMatrix(otherSpeciesID), 0.5_8, method="SVD"  )
          otherAuxMatrixB = otherAuxMatrix
          otherAuxMatrix = Matrix_product_dgemm( Matrix_product_dgemm( otherAuxMatrixB , this%densityMatrix(otherSpeciesID)), otherAuxMatrixB )
       end if

    case default

    end select

    if(trim(speciesName) .eq. "E-ALPHA") then
       do i=1, numberOfcontractions
          output%values(i,1) = auxMatrix%values(i,i)+otherAuxMatrix%values(i,i)
          output%values(i,2) = auxMatrix%values(i,i)-otherAuxMatrix%values(i,i)
       end do

    else
       do i=1, numberOfcontractions
          output%values(i,1) = auxMatrix%values(i,i)
       end do
    end if
    !print*,"auxMatrix%values", auxMatrix%values      
    if( present( totalSum ) ) totalSum(1) = sum(output%values(:,1))
    if( present( totalSum ) .and. trim(speciesName) .eq. "E-ALPHA") totalSum(2) = sum(output%values(:,2))

    call Matrix_destructor(auxMatrix)
    call Matrix_destructor(auxMatrixB)
    call Matrix_destructor(otherAuxMatrix)
    call Matrix_destructor(otherAuxMatrixB)

  end function CalculateProperties_getPopulation

  !<
  !! @brief Muestra las contrinuciones al dipolo de cada especie
  !>
  subroutine CalculateProperties_showExpectedPositions(this)
    implicit none
    type(CalculateProperties) :: this
    integer :: i, numberOfSpecies
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    print *,""
    print *," EXPECTED POSITIONS OF QUANTUM SPECIES:"
    print *,"======================"
    print *,""
    print *,"POSITIONS IN ANGSTROMS"
    print *,"------"
    print *,""
    write (*,"(T19,4A9)") "<x>","<y>", "<z>", ""
    do i=1, numberOfSpecies
       write (*,"(T5,A15,3F9.4)") trim(MolecularSystem_getSymbolOfSpecies( i )), CalculateProperties_getExpectedPosition(this, i)* ANGSTROM
    end do
    if(trim(CONTROL_instance%UNITS) .eq. "BOHR") then
       print *,""
       print *,"POSITIONS IN BOHR"
       print *,"------"
       print *,""
       write (*,"(T19,4A9)") "<x>","<y>", "<z>", ""
       do i=1, numberOfSpecies
          write (*,"(T5,A15,3F9.4)") trim(MolecularSystem_getSymbolOfSpecies( i )), CalculateProperties_getExpectedPosition(this, i)
       end do
    end if
    print *,""
    print *,"END EXPECTED POSITIONS"
    print *,""
  end subroutine CalculateProperties_showExpectedPositions

  function CalculateProperties_getExpectedPosition( this , speciesID) result(output)
    implicit none
    type(CalculateProperties) :: this
    integer :: speciesID
    real(8) :: output(3)

    !! Open file for wavefunction                                                                                                          
    output=0.0_8

    output(1)=sum( this%densityMatrix(speciesID)%values * this%momentMatrices(speciesID,1)%values ) 
    output(2)=sum( this%densityMatrix(speciesID)%values * this%momentMatrices(speciesID,2)%values ) 
    output(3)=sum( this%densityMatrix(speciesID)%values * this%momentMatrices(speciesID,3)%values ) 

  end function CalculateProperties_getExpectedPosition

  !<
  !! @brief Muestra las contrinuciones al dipolo de cada especie
  !>
  subroutine CalculateProperties_showContributionsToElectrostaticMoment(this)
    implicit none
    type(CalculateProperties) :: this

    integer :: i, numberOfSpecies
    real(8), allocatable :: dipole(:,:)
    real(8), allocatable :: quadrupole(:,:)
    real(8) :: totalDipole(3)
    real(8) :: totalQuadrupole(6)

    totalDipole=0.0_8
    totalQuadrupole=0.0_8
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

    allocate(dipole(numberOfSpecies+1,3))
    allocate(quadrupole(numberOfSpecies+1,6))

    print *,""
    print *," ELECTROSTATIC MOMENTS:"
    print *,"======================"
    print *,""
    print *,""
    print *,"DIPOLE: (A.U.)"
    print *,"------"
    print *,""
    write (*,"(T19,4A13)") "<Dx>","<Dy>", "<Dz>"," |D|"

    do i=1, numberOfSpecies
       dipole(i,:)=CalculateProperties_getDipoleOfQuantumSpecies(this, i)
       totalDipole(:)=totalDipole(:)+dipole(i,:)
       write (*,"(T5,A15,3F13.8)") trim(MolecularSystem_getSymbolOfSpecies( i )), dipole(i,:)
    end do
    dipole(numberOfSpecies+1,:)=CalculateProperties_getDipoleOfPuntualCharges()
    totalDipole(:)=totalDipole(:)+dipole(numberOfSpecies+1,:)
    write (*,"(T5,A15,3F13.8)") "Point charges: ", dipole(numberOfSpecies+1,:)
    write (*,"(T22,A28)") "___________________________________"
    write (*,"(T5,A15,3F13.8, F13.8)") "Total Dipole:", totalDipole(:), sqrt(sum(totalDipole(:)**2.0 ) )

    print *,""
    print *,"DIPOLE: (DEBYE)"
    print *,"------"
    print *,""
    write (*,"(T19,4A13)") "<Dx>","<Dy>", "<Dz>"," |D|"

    totalDipole=0.0_8
    do i=1, numberOfSpecies
       dipole(i,:)=CalculateProperties_getDipoleOfQuantumSpecies(this, i)*DEBYE
       totalDipole(:)=totalDipole(:)+dipole(i,:)
       write (*,"(T5,A15,3F13.8)") trim(MolecularSystem_getSymbolOfSpecies( i )), dipole(i,:)
    end do

    dipole(numberOfSpecies+1,:)=CalculateProperties_getDipoleOfPuntualCharges()*DEBYE
    totalDipole(:)=totalDipole(:)+dipole(numberOfSpecies+1,:)
    write (*,"(T5,A15,3F13.8)") "Point charges: ", dipole(numberOfSpecies+1,:)

    write (*,"(T22,A28)") "___________________________________"

    write (*,"(T5,A15,3F13.8, F13.8)") "Total Dipole:", totalDipole(:), sqrt(sum(totalDipole(:)**2.0 ) )


    print *,""
    print *,"QUADRUPOLE NON-TRACELESS: (DEBYE ANGS)"
    print *,"------"
    print *,""
    write (*,"(T19,6A13)") "<xx>","<yy>", "<zz>", "<xy>","<xz>","<yz>"

    do i=1, numberOfSpecies
       quadrupole(i,:)=CalculateProperties_getQuadrupoleOfQuantumSpecies(this, i)*DEBYE*ANGSTROM
       totalQuadrupole(:)=totalQuadrupole(:)+quadrupole(i,:)
       write (*,"(T5,A15,6F14.8)") trim(MolecularSystem_getSymbolOfSpecies( i )), quadrupole(i,:)
    end do

    quadrupole(numberOfSpecies+1,:)=CalculateProperties_getQuadrupoleOfPuntualCharges()*DEBYE*ANGSTROM
    totalquadrupole(:)=totalquadrupole(:)+quadrupole(numberOfSpecies+1,:)
    write (*,"(T5,A15,6F14.8)") "Point charges: ", quadrupole(numberOfSpecies+1,:)

    write (*,"(T2,A18,6F14.8)") "Total Quadrupole:", totalQuadrupole(:) 

    write (*,"(T22,A28)") "___________________________________"

    !write (*,"(T5,A15,3F13.8, F13.8)") "Total Dipole:", totalDipole(:), sqrt(sum(totalDipole(:)**2.0 ) )


    print *,""
    print *,"END ELECTROSTATIC MOMENTS"
    print *,""

    deallocate(dipole)
    deallocate(quadrupole)

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


  ! !<
  ! !! @brief Calcula el aporte al dipolo de las cargas puntuales presentes
  ! !>
  function CalculateProperties_getQuadrupoleOfPuntualCharges() result( output )
    implicit none
    real(8) :: output(6)
    integer :: i

    output = 0.0_8
    do i=1, size( MolecularSystem_instance%pointCharges )      
       output(1) = output(1) + MolecularSystem_instance%pointCharges(i)%origin(1)* MolecularSystem_instance%pointCharges(i)%origin(1)* MolecularSystem_instance%pointCharges(i)%charge
       output(2) = output(2) + MolecularSystem_instance%pointCharges(i)%origin(2)* MolecularSystem_instance%pointCharges(i)%origin(2)* MolecularSystem_instance%pointCharges(i)%charge
       output(3) = output(3) + MolecularSystem_instance%pointCharges(i)%origin(3)* MolecularSystem_instance%pointCharges(i)%origin(3)* MolecularSystem_instance%pointCharges(i)%charge
       output(4) = output(4) + MolecularSystem_instance%pointCharges(i)%origin(1)* MolecularSystem_instance%pointCharges(i)%origin(2)* MolecularSystem_instance%pointCharges(i)%charge
       output(5) = output(5) + MolecularSystem_instance%pointCharges(i)%origin(1)* MolecularSystem_instance%pointCharges(i)%origin(3)* MolecularSystem_instance%pointCharges(i)%charge
       output(6) = output(6) + MolecularSystem_instance%pointCharges(i)%origin(2)* MolecularSystem_instance%pointCharges(i)%origin(3)* MolecularSystem_instance%pointCharges(i)%charge
    end do
  end function CalculateProperties_getQuadrupoleOfPuntualCharges


  !<
  !! @brief calcula el aporte al dipolo debido a particulas no fijas
  !>
  function CalculateProperties_getDipoleOfQuantumSpecies( this, i ) result( output )
    implicit none
    type(calculateproperties) :: this
    integer :: i !specieid
    real(8) :: output(3)

    output = CalculateProperties_getExpectedPosition(this, i) * molecularsystem_getcharge( i )

  end function CalculateProperties_getDipoleOfQuantumSpecies


  !<
  !! @brief calcula el aporte al dipolo debido a particulas no fijas
  !>
  function CalculateProperties_getQuadrupoleOfQuantumSpecies( this, i ) result( output )
    implicit none
    type(calculateproperties) :: this
    integer :: i !specieid
    real(8) :: output(6)

    output(1) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,4)%values )
    output(2) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,5)%values )
    output(3) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,6)%values )
    output(4) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,7)%values )
    output(5) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,8)%values )
    output(6) =sum( this%densitymatrix(i)%values * this%momentmatrices(i,9)%values )

    output = output * molecularsystem_getcharge( i )

  end function CalculateProperties_getQuadrupoleOfQuantumSpecies

  !<
  !! @brief  Calculates the expected position^2 for each quantum species
  !! @author F.M. mar-2025  
  !>
  function CalculateProperties_getExpectedR2(this,speciesID) result(output)
    implicit none
    type(CalculateProperties) :: this
    integer :: speciesID
    real(8) :: output

    real(8) :: orig(3)
    type(Matrix) :: harmonicIntegrals

    orig=CalculateProperties_getExpectedPosition( this , speciesID)
    call DirectIntegralManager_getHarmonicIntegrals(MolecularSystem_instance,speciesID,orig,harmonicIntegrals)
    output= sum( this%densitymatrix(speciesID)%values*harmonicIntegrals%values )

  end function CalculateProperties_getExpectedR2

  !<
  !! @brief  Shows the square root of the expected position^2 for each quantum species
  !! @author F.M. mar-2025  
  !>
  subroutine CalculateProperties_showRMSradius(this)
    implicit none
    type(CalculateProperties) :: this
    real(8) :: R2
    integer :: numberOfSpecies   
    integer :: speciesID

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    print *,""
    print *," EXPECTED RMS RADIUS OF QUANTUM SPECIES:"
    print *,"======================"
    print *,"Relative to the position expectation value"
    print *,"Computed as sqrt(<R^2 - <R>>)"
    if(trim(CONTROL_instance%UNITS) .eq. "BOHR") then
       write (*,"(T20,A)") "IN BOHR"
       write (*,"(T20,A)") "------"
       do speciesID=1, numberOfSpecies
          R2=CalculateProperties_getExpectedR2(this,speciesID)
          if(sqrt(R2).gt.1E-3 .and. sqrt(R2).lt.1E5) then 
             write (*,"(T5,A15,F12.6)") trim(MolecularSystem_getSymbolOfSpecies(speciesID)), sqrt(R2)
          else
             write (*,"(T5,A15,ES12.5)") trim(MolecularSystem_getSymbolOfSpecies(speciesID)), sqrt(R2)
          end if
       end do
    end if
    print *,""
    write (*,"(T20,A)") "IN ANGSTROM"
    write (*,"(T20,A)") "-----------"
    do speciesID=1, numberOfSpecies
       R2=CalculateProperties_getExpectedR2(this,speciesID)
       if(sqrt(R2).gt.1E-3 .and. sqrt(R2).lt.1E5) then 
          write (*,"(T5,A15,F12.6)") trim(MolecularSystem_getSymbolOfSpecies(speciesID)), sqrt(R2)*ANGSTROM
       else
          write (*,"(T5,A15,ES12.5)") trim(MolecularSystem_getSymbolOfSpecies(speciesID)), sqrt(R2)*ANGSTROM
       end if
    end do
    print *,""
    print *,"END EXPECTED RMS RADIUS"
    print *,""

  end subroutine CalculateProperties_showRMSradius

end module CalculateProperties_

