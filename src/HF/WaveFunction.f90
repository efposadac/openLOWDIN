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
!! @brief Modulo para definicion de funciones de onda RHF
!!  Este modulo define funciones de onda para diferentes especies cuanticas dentro del formalismo OMNE
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-30
!! <b> Historial de modificaciones: </b>
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y metodos
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el modulo para su inclusion en Lowdin
module WaveFunction_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use MolecularSystem_
  implicit none


  !< enum Matrix_type {
  integer, parameter :: CANONICAL_ORTHOGONALIZATION	= 1
  integer, parameter :: SYMMETRIC_ORTHOGONALIZATION	= 2
  !< }

  !< enum type of orbital graph {
  integer, parameter, public :: ORBITAL_ALONE = 1
  integer, parameter, public :: ORBITAL_WITH_POTENTIAL = 2
  !< }

  type, public :: WaveFunction
     
     !!**************************************************************
     !! Matrices asociadas al calculo de funciones de onda RHF
     !!
     type(Matrix) :: overlapMatrix
     type(Matrix) :: transformationMatrix
     type(Matrix) :: kineticMatrix
     type(Matrix) :: puntualInteractionMatrix
     type(Matrix) :: HcoreMatrix
     type(Matrix) :: densityMatrix
     type(Matrix) :: twoParticlesMatrix
     type(Matrix) :: couplingMatrix
     type(Matrix) :: externalPotentialMatrix
     type(Matrix) :: coefficientsofcombination
     type(vector) :: energyofmolecularorbital
     
     !!**************************************************************


     !!**************************************************************
     !! Atributos asociados a los valores de energia al final de
     !! la funcion de onda RHF
     real(8) :: totalEnergyForSpecie
     real(8) :: independentSpecieEnergy
     real(8) :: kineticEnergy
     real(8) :: puntualInteractionEnergy
     real(8) :: independentParticleEnergy
     real(8) :: repulsionEnergy
     real(8) :: couplingEnergy
     real(8) :: externalPotentialEnergy

     !!**************************************************************

  end type WaveFunction
  
  type(WaveFunction), public, allocatable, target :: WaveFunction_instance(:)
  
  public

  
contains


  !>
  !! @brief Define el constructor para la clase
  subroutine WaveFunction_constructor()
    implicit none
    
    integer(8) :: numberOfContractions
    integer :: speciesID

    allocate( WaveFunction_instance( MolecularSystem_instance%numberOfQuantumSpecies ) )

    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

       WaveFunction_instance( speciesID )%totalEnergyForSpecie = 0.0_8
       WaveFunction_instance( speciesID )%independentSpecieEnergy =0.0_8
       WaveFunction_instance( speciesID )%kineticEnergy = 0.0_8
       WaveFunction_instance( speciesID )%puntualInteractionEnergy = 0.0_8
       WaveFunction_instance( speciesID )%independentParticleEnergy = 0.0_8
       WaveFunction_instance( speciesID )%repulsionEnergy = 0.0_8
       WaveFunction_instance( speciesID )%externalPotentialEnergy = 0.0_8
       WaveFunction_instance( speciesID )%couplingEnergy = 0.0_8
       
       call Matrix_constructor( WaveFunction_instance(speciesID)%externalPotentialMatrix, numberOfContractions, numberOfContractions )
       
    end do

  end subroutine WaveFunction_constructor

  !>
  !! @brief Contruye la matrix de overlap.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_buildOverlapMatrix(file, speciesID)
    implicit none
    
    character(*), intent(in) :: file
    integer, intent(in) :: speciesID
    
    integer :: unit
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)
    
    arguments(1) = "OVERLAP"
    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
    
    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    
    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          
    
    WaveFunction_instance( speciesID )%overlapMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)
    
    close(34)
    
    !! DEBUG
    ! print *,"Matriz de overlap: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show(WaveFunction_instance( speciesID )%overlapMatrix)
    
  end subroutine WaveFunction_buildOverlapMatrix
  
   
  !>
  !! @brief Contruye la matrix de de transformacion.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_buildTransformationMatrix(file, speciesID, typeOfOrthogonalization )
    implicit none
    
    character(*), intent(in) :: file
    integer, intent(in) :: speciesID    
    integer, optional, intent(in) :: typeOfOrthogonalization
    
    type(Matrix) :: eigenVectors
    type(Vector) :: eigenValues
    integer(8) :: numberOfContractions
    integer :: i, j
    
    !! Numero de contracciones "totales"
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
    
    if ( .not. allocated( WaveFunction_instance( speciesID )%overlapMatrix%values ) ) then
          
       call WaveFunction_buildOverlapMatrix( trim(file), speciesID )
       
    end if
    
    if ( .not. allocated(WaveFunction_instance( speciesID )%transformationMatrix%values)) then
          
       call Matrix_constructor( WaveFunction_instance( speciesID )%transformationMatrix, &
            int(numberOfContractions,8), int(numberOfContractions,8), Math_NaN )
          
    end if
       
    if ( numberOfContractions > 1) then
          
       call Vector_constructor( eigenValues, int(numberOfContractions) )
          
       call Matrix_constructor( eigenVectors, numberOfContractions, numberOfContractions)
          
       !!****************************************************************
       !! diagonaliza la matriz de overlap obteniendo una matriz unitaria
       !!          
       call Matrix_eigen( WaveFunction_instance( speciesID )%overlapMatrix, eigenValues, eigenVectors, SYMMETRIC  )
       
       do i = 1 , numberOfContractions
          do j = 1 , numberOfContractions
                
             WaveFunction_instance( speciesID )%transformationMatrix%values(i,j) = &
                  eigenVectors%values(i,j)/sqrt( eigenValues%values(j) )
             
          end do
       end do
       !!
       !!****************************************************************
       
       !!****************************************************************
       !! Calcula matriz de transformacion
       !!
       select case (typeOfOrthogonalization)
          
       !! Ortogonalizacion canonica
       case (CANONICAL_ORTHOGONALIZATION)
          
          WaveFunction_instance( speciesID )%transformationMatrix%values = &
               WaveFunction_instance( speciesID )%transformationMatrix%values
          
       !!Ortogonalizacion simetrica
       case (SYMMETRIC_ORTHOGONALIZATION)
             
          WaveFunction_instance( speciesID )%transformationMatrix%values  = &
               matmul(WaveFunction_instance( speciesID )%transformationMatrix%values, transpose(eigenVectors%values))
          
       case default
             
          WaveFunction_instance( speciesID )%transformationMatrix%values  = &
               matmul(WaveFunction_instance( speciesID )%transformationMatrix%values, transpose(eigenVectors%values))
          
       end select
       
       call Vector_destructor( eigenValues )
       call Matrix_destructor( eigenVectors )
       !!
       !!****************************************************************
       
    else
       
       WaveFunction_instance( speciesID )%transformationMatrix%values= 1.0_8
       
    end if
    
    !! Debug
    ! print *,"Matriz de transformacion: ", trim(MolecularSystem_instance%species(speciesID)%symbol)
    ! call Matrix_show( WaveFunction_instance( speciesID )%transformationMatrix )
       
  end subroutine WaveFunction_buildTransformationMatrix

  !>
  !! @brief Contruye la matrix de particula independiente.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  subroutine WaveFunction_HCoreMatrix(file, speciesID)
    implicit none
    
    character(*), intent(in) :: file
    integer, intent(in) :: speciesID
    
    integer :: unit
    integer :: k, l, r, s
    integer :: ParticleID, ParticleID_2
    integer :: contractionID, contractionID_2
    integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    integer :: owner, owner_2
    integer :: auxCharge
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    real(8) :: aux

    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")
    
    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))    

    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          
    
    !! Load Kinetic Matrix
    arguments(1) = "KINETIC"    

    WaveFunction_instance( speciesID )%kineticMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)

    !! Incluiding mass effect       
    if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
       WaveFunction_instance( speciesID )%kineticMatrix%values =  &
            WaveFunction_instance( speciesID )%kineticMatrix%values * &
            ( 1.0_8/MolecularSystem_getMass( speciesID ) -1.0_8 / ParticleManager_getTotalMass() )
       
    else
       WaveFunction_instance( speciesID )%kineticMatrix%values = &
            WaveFunction_instance( speciesID )%kineticMatrix%values / &
            MolecularSystem_getMass( speciesID )
    end if
    
    !! Finite Nuclear Mass Correction
    if ( CONTROL_instance%FINITE_MASS_CORRECTION ) then
       k=1
       do particleID = 1, size(MolecularSystem_instance%species(speciesID)%particles)
          do contractionID = 1, size(MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction)
             
             numberOfCartesiansOrbitals = MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction(contractionID)%numCartesianOrbital
             owner = MolecularSystem_instance%species(speciesID)%particles(particleID)%basis%contraction(contractionID)%owner
             
             do s = 1, numberOfCartesiansOrbitals
                l=k
                
                do particleID_2 = particleID, size(MolecularSystem_instance%species(speciesID)%particles)
                   do contractionID_2 = contractionID, size(MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction)
                      
                      numberOfCartesiansOrbitals_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%numCartesianOrbital
                      owner_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%owner
             
                      do r = 1, numberOfCartesiansOrbitals_2
                         
                         if ( owner .eq. owner_2) then

                            WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)=&
                                 WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)*&
                                 ( 1 + MolecularSystem_getMass( speciesID ) / MolecularSystem_getOwnerMass(owner) )
                            
                            WaveFunction_instance( speciesID )%kineticMatrix%values(l,k) = WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)

                         end if
                         
                      end do
                      l=l+1
                      
                   end do
                end do
                k=k+1

             end do

          end do
       end do

    end if
   
    !! DEBUG
    !print *,"Matriz de energia cinetica: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !call Matrix_show( WaveFunction_instance(speciesID)%kineticMatrix )

    !! Load N-Q- Attraction  Matrix
    arguments(1) = "ATTRACTION"

    WaveFunction_instance( speciesID )%puntualInteractionMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    

    !! Incluiding charge effect
    auxCharge = MolecularSystem_getCharge( speciesID )
    
    WaveFunction_instance( speciesID )%puntualInteractionMatrix%values = &
         WaveFunction_instance( speciesID )%puntualInteractionMatrix%values * (-auxCharge)

    !! DEBUG
    !print *,"Matriz de interaccion n-e: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !call Matrix_show( WaveFunction_instance(speciesID)%puntualInteractionMatrix )
    
    close(34)    
    
    !! Build Hcore Matrix
    if ( .not.allocated(WaveFunction_instance( speciesID )%HcoreMatrix%values ) ) then
       
       call Matrix_constructor( WaveFunction_instance( speciesID )%HcoreMatrix, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
       
    end if
    
    WaveFunction_instance(speciesID)%HCoreMatrix%values = &
         WaveFunction_instance(speciesID)%kineticMatrix%values + &
         WaveFunction_instance(speciesID)%puntualInteractionMatrix%values
         
    !! DEBUG
    !print *,"Matriz de hcore: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !call Matrix_show( WaveFunction_instance( speciesID )%HcoreMatrix )
    
  end subroutine WaveFunction_HCoreMatrix

  !>
  !! @brief Ajusta la matriz de densidad para una especie espcificada
  subroutine WaveFunction_setDensityMatrix( densityMatrix, speciesID )
    implicit none
    
    type(Matrix), intent(in) :: densityMatrix
    integer :: speciesID
    
    integer :: totalNumberOfContractions
    
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)
    
    if( .not. allocated(WaveFunction_instance( speciesID )%densityMatrix%values )) then
       call Matrix_constructor( WaveFunction_instance( speciesID )%densityMatrix, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
    end if
    
    call Matrix_copyConstructor( WaveFunction_instance( speciesID )%densityMatrix, densityMatrix )

    !! Debug
    ! print*, "Matriz de densidad inicial ", MolecularSystem_getNameOfSpecie(speciesID)
    ! call Matrix_show(WaveFunction_instance( speciesID )%densityMatrix)
    
  end subroutine WaveFunction_setDensityMatrix
  


  !>
  !! @brief Calcula las componentes de energia para la especie especificada
  !! @warning Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  subroutine WaveFunction_obtainEnergyComponents(specieID )
    implicit none
    
    integer :: specieID

    !! Calcula la energia de repulsion
    WaveFunction_instance( specieID )%repulsionEnergy = 0.5_8 * &
         sum( transpose( WaveFunction_instance( specieID )%densityMatrix%values ) * &
         WaveFunction_instance( specieID )%twoParticlesMatrix%values )
    
    !! Calcula energia de particula independiente
    WaveFunction_instance( specieID )%independentParticleEnergy = &
         sum( transpose( WaveFunction_instance(specieID)%densityMatrix%values ) * &
         WaveFunction_instance( specieID )%hcoreMatrix%values )

    !! Calcula energia cinetica para la especie dada
    WaveFunction_instance( specieID )%kineticEnergy = &
         sum( transpose(WaveFunction_instance( specieID )%densityMatrix%values) * &
         WaveFunction_instance( specieID )%kineticMatrix%values )

    !! Calcula energia de potencial externo para la especie dada
    WaveFunction_instance( specieID )%externalPotentialEnergy = &
         sum( transpose(WaveFunction_instance( specieID )%densityMatrix%values) * &
         WaveFunction_instance( specieID )%externalPotentialMatrix%values )
    
    !! Calcula energia de interaccion entre particulas puntuales y cuanticas
    WaveFunction_instance( specieID )%puntualInteractionEnergy =  &
         WaveFunction_instance( specieID )%independentParticleEnergy - &
         WaveFunction_instance( specieID )%kineticEnergy

    !! Calula enegia de especie independiente (  sin considerar el termino de acoplamiento )
    WaveFunction_instance( specieID )%independentSpecieEnergy = &
         sum( transpose(WaveFunction_instance( specieID )%densityMatrix%values) * &
         (  ( WaveFunction_instance( specieID )%hcoreMatrix%values ) + &
         0.5_8 * WaveFunction_instance( specieID )%twoParticlesMatrix%values + &
         WaveFunction_instance( specieID )%externalPotentialMatrix%values))
    
    
    WaveFunction_instance( specieID )%independentSpecieEnergy = &
         WaveFunction_instance( specieID )%independentSpecieEnergy 

    !! Calcula energia de acoplamiento en caso de mas de una especie presente
    WaveFunction_instance( specieID )%couplingEnergy = &
         sum( transpose( WaveFunction_instance( specieID )%densityMatrix%values ) * &
         WaveFunction_instance( specieID )%couplingMatrix%values )

    !! Total energy for species
    WaveFunction_instance( specieID )%totalEnergyForSpecie = &
         WaveFunction_instance( specieID )%independentSpecieEnergy +  &
         WaveFunction_instance( specieID )%couplingEnergy


    ! print *, "__________________ ENERGY COMPONENTS _______________________"
    ! print *, "	Specie                       ", MolecularSystem_getNameOfSpecie( specieID )
    ! print *, "	Total Energy                =", WaveFunction_instance( specieID )%totalEnergyForSpecie
    ! print *, "	Indepent Specie Energy      =", WaveFunction_instance( specieID )%independentSpecieEnergy
    ! print *, "	Kinetic Energy              =",WaveFunction_instance( specieID )%kineticEnergy
    ! print *, "	Puntual Interaction Energy  =",WaveFunction_instance( specieID )%puntualInteractionEnergy
    ! print *, "	Independent Particle Energy =",WaveFunction_instance( specieID )%independentParticleEnergy
    ! print *, "	Repultion Energy            =",WaveFunction_instance( specieID )%repulsionEnergy
    ! print *, "	Coupling Energy             =", WaveFunction_instance( specieID )%couplingEnergy
    ! print *, "____________________________________________________________"

  end subroutine WaveFunction_obtainEnergyComponents

  !>
  !! @brief retorna la matrix de particula independiente.
  !! @param nameOfSpecie nombre de la especie seleccionada.
  function WaveFunction_getHcoreMatrix( nameOfSpecie ) result(output)
    implicit none
    character(*), optional :: nameOfSpecie
    type(Matrix) :: output

    integer :: specieID
    character(30) :: nameOfSpecieSelected
    
    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )
    
    specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    if ( .not. allocated(WaveFunction_instance( specieID )%HcoreMatrix%values ) ) then
       
       call WaveFunction_exception(ERROR, "You need to build the Hcore matrix before to use it.", "At HF program, at WaveFunction_getHcoreMatrix function.")
       
    end if

    call Matrix_copyConstructor( output, WaveFunction_instance(specieID)%HcoreMatrix )

  end function WaveFunction_getHcoreMatrix

  !>
  !! @brief Retorna la matrix de transformationMatrix
  function WaveFunction_getTransformationMatrix( nameOfSpecie ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    type(Matrix) ::  output

    character(30) :: nameOfSpecieSelected
    integer :: specieID

    nameOfSpecieSelected = "e-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
    
    if ( .not.allocated(WaveFunction_instance(specieID)%transformationMatrix%values) ) then

       call WaveFunction_exception(ERROR, "You need build the transformation Matrix before to use it.", "At HF program, at WaveFunction_getTransformationMatrix function.")

    end if

    call Matrix_copyConstructor( output, WaveFunction_instance(specieID)%transformationMatrix )

  end function WaveFunction_getTransformationMatrix

!   !<
!   !! @brief Contruye una matriz de interaccion con un potencial externo
!   !!
!   !! @param nameOfSpecie nombre de la especie seleccionada.
!   !>
!   subroutine WaveFunction_buildExternalPotentialMatrix( nameOfSpecie )
!     implicit none
!     character(*), optional :: nameOfSpecie

!     character(30) :: nameOfSpecieSelected
!     integer :: specieID

!     nameOfSpecieSelected = "e-"
!     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

!     specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

!     if ( nameOfspecie /= "e-BETA" ) then

! 	if( WaveFunction_instance(specieID)%isThereExternalPotential ) then
! 		WaveFunction_instance(specieID)%externalPotentialMatrix%values = 0.0_8


! 	if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )	then	!! Numerical integration
! 		if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
! 			WaveFunction_instance(specieID)%externalPotentialMatrix = &
! 				IntegralManager_getNumericalInteractionWithPotentialMatrix( &
! 				ExternalPotential_Manager_instance%externalsPots, specieID, integralName="external" )

! 		else 		!! From xml file
! 			WaveFunction_instance(specieID)%externalPotentialMatrix = &
! 				IntegralManager_getNumericalPotentialMatrixFromXml( &
! 				ExternalPotential_Manager_instance%externalsPots, specieID, integralName="external" )
! 		end if
! 	else		!! Analytical Integration	

! 		WaveFunction_instance(specieID)%externalPotentialMatrix = &
! 		IntegralManager_getInteractionWithPotentialMatrix( &
! 		ExternalPotential_Manager_instance%externalsPots, specieID, "external" )

!           end if

!        end if

!    else !! Use the same matrix for e-beta and e-alpha

!        WaveFunction_instance(specieID)%externalPotentialMatrix = &
!             WaveFunction_instance( MolecularSystem_getSpecieID( nameOfSpecie="e-ALPHA" ))%externalPotentialMatrix

!     end if

!     			print *,"EXTERNAL POTENTIAL MATRIX FOR: ", nameOfSpecie
!     			call Matrix_show(WaveFunction_instance(specieID)%externalPotentialMatrix)

!   end subroutine WaveFunction_buildExternalPotentialMatrix


!   function WaveFunction_getValueForOrbitalAt( nameOfSpecie, orbitalNum, coordinate ) result(output)
!     implicit none
!     character(*), optional, intent(in) :: nameOfSpecie
!     integer :: orbitalNum
!     real(8) :: coordinate(3)
!     real(8) :: output

!     integer :: specieID
!     character(30) :: nameOfSpecieSelected
!     integer :: numberOfContractions
!     integer :: totalNumberOfContractions
!     integer :: particleID
!     integer :: contractionID
!     integer :: i, j
!     real(8), allocatable :: auxVal(:)


!     nameOfSpecieSelected = "e-"
!     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

!     specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

!     numberOfContractions = MolecularSystem_getNumberOfContractions( specieID )

!     output=0.0_8
!     do i=1,numberOfContractions

!        particleID = MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%particleID
!        contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(specieID)%contractionID(i)%contractionIDInParticle

!        totalNumberOfContractions = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital

!        if( allocated(auxVal)) deallocate(auxVal)
!        allocate(auxVal(totalNumberOfContractions))

!        auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( specieID,  numberOfContraction=i ), coordinate )

!        do j = 1, totalNumberOfContractions

!           output = output + auxVal(j) * WaveFunction_instance( specieID )%waveFunctionCoefficients%values(j,orbitalNum)

!        end do

!     end do


!   end function WaveFunction_getValueForOrbitalAt
!   !
!   !
!   subroutine WaveFunction_draw2DOrbital( nameOfSpecie, orbitalNum, flags )
!     implicit none
!     character(*), optional, intent(in) :: nameOfSpecie
!     integer :: orbitalNum
!     integer :: flags


!     character(30) :: nameOfSpecieSelected
!     character(50) :: fileName
!     character(50) :: xRange
!     integer :: specieID
!     integer :: numberOfContractions
!     integer :: j
!     integer :: i
!     integer :: numOfGraphs
!     integer :: auxInitOrbitalNum
!     integer :: auxLastOrbitalNum

!     ! 	nameOfSpecieSelected = "e-"
!     ! 	if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

!     ! 	specieID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

!     ! 	numberOfContractions = MolecularSystem_getTotalNumberOfContractions( specieID )
!     ! 	fileName=trim(CONTROL_instance%INPUT_FILE)//'orbital.'//trim(String_convertIntegerToString(orbitalNum))//'.'//trim(nameOfSpecieSelected)

!     ! 	select case( flags )

!     ! 		case(ORBITAL_ALONE)

!     ! 			open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
!     ! 			do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
!     ! 				CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
!     ! 			write (5,"(F20.10,F20.10)") j*CONTROL_instance%STEP_OF_GRAPHS, &
!     ! 			WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected,&
!     ! 			orbitalNum, [0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) !&
!     ! !			!!! + (WaveFunction_instance( specieID )%molecularOrbitalsEnergy%values(orbitalNum) * CM_NEG1)
!     ! 			end do
!     ! 			close(5)
!     ! 			call InputOutput_make2DGraph(trim(fileName),&
!     ! 				"Nuclear Wave Function",&
!     ! 				"r / Bohr",&
!     ! 				"U / a.u.", &
!     ! 				y_format="%.2e")

!     ! 		case(ORBITAL_WITH_POTENTIAL)

!     ! 			auxInitOrbitalNum=orbitalNum
!     ! 			auxLastOrbitalNum=orbitalNum
!     ! 			numOfGraphs=2
!     ! 			if(orbitalNum==0) then
!     ! 				auxInitOrbitalNum=1
!     ! 				auxLastOrbitalNum=numberOfContractions
!     ! 				numOfGraphs=numberOfContractions+1
!     ! 			end if
!     ! 		open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
!     ! 		do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
!     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
!     ! 			write (5,"(2ES20.10$)") &
!     ! 			j*CONTROL_instance%STEP_OF_GRAPHS, &
!     ! 			ExternalPotential_getPotential(ExternalPotential_Manager_instance%externalsPots(1),&
!     ! 			[j*CONTROL_instance%STEP_OF_GRAPHS,0.0_8,0.0_8])*CM_NEG1
!     ! 			do i=auxInitOrbitalNum,auxLastOrbitalNum
!     ! 				write (5,"(ES20.10$)") &
!     ! 				CONTROL_instance%WAVE_FUNCTION_SCALE&
!     ! 				*WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected, i,&
!     ! 					[0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) &
!     ! 				+ (WaveFunction_instance( specieID )%molecularOrbitalsEnergy%values(i) * CM_NEG1)
!     ! 			end do
!     ! 			write (5,"(A)") ""
!     ! 		end do
!     ! 		close(5)

!     ! 		xRange=trim(adjustl(String_convertRealToString(real(&
!     ! 			-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
!     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))//':'//trim(adjustl(String_convertRealToString(real(&
!     ! 			CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS&
!     ! 			*CONTROL_instance%STEP_OF_GRAPHS,8))))

!     ! 			call InputOutput_make2DGraph(trim(fileName),&
!     ! 			"Nuclear Wave Function in potential ",&
!     ! 			"r / Bohr",&
!     ! 			"U / cm-1",&
!     ! 			y_format="%.2e",numOfGraphs=numOfGraphs,x_range=trim(xRange))

!     ! 	end select

!   end subroutine WaveFunction_draw2DOrbital

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine WaveFunction_exception( typeMessage, description, debugDescription)
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
  end subroutine WaveFunction_exception
  
end module WaveFunction_
