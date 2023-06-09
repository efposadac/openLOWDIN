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
  use CosmoCore_
  implicit none


  !< enum Matrix_type {
  integer, parameter :: CANONICAL_ORTHOGONALIZATION = 1
  integer, parameter :: SYMMETRIC_ORTHOGONALIZATION = 2
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
     type(Matrix) :: exchangeCorrelationMatrix
     type(Matrix) :: externalPotentialMatrix
     type(Matrix) :: coefficientsofcombination
     type(vector) :: energyofmolecularorbital
     !! Cosmo Things
     type(Matrix) :: cosmo1
     type(Matrix) :: cosmo2
     type(Matrix) :: cosmo4
     type(Matrix) :: cosmoCoupling
     type(Matrix) :: electricField(9)
     type(Matrix) :: harmonic
     real(8) :: cosmoCharge
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
     real(8) :: exchangeCorrelationEnergy
     !! Cosmo Things
     real(8) :: cosmoEnergy
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
       WaveFunction_instance( speciesID )%exchangeCorrelationEnergy = 0.0_8

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
    if (  CONTROL_instance%DEBUG_SCFS) then
     print *,"Matriz de overlap: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
     call Matrix_show(WaveFunction_instance( speciesID )%overlapMatrix)
    end if

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

       !!do i = 1 , numberOfContractions
       !!   print *, eigenvalues%values(i) 
       !!end do

       do i = 1 , numberOfContractions
         do j = 1 , numberOfContractions
           if ( abs(eigenValues%values(j)) >= CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
             WaveFunction_instance( speciesID )%transformationMatrix%values(i,j) = &
                            eigenVectors%values(i,j)/sqrt( eigenvalues%values(j) )
           else
             WaveFunction_instance( speciesID )%transformationMatrix%values(i,j) = 0
           end if
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
    real(8) :: auxCharge
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(40) :: arguments(2)
    type(matrix) :: auxmatrix

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

                do particleID_2 = 1, size(MolecularSystem_instance%species(speciesID)%particles)
                   do contractionID_2 = 1, size(MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction)

                      numberOfCartesiansOrbitals_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%numCartesianOrbital
                      owner_2 = MolecularSystem_instance%species(speciesID)%particles(particleID_2)%basis%contraction(contractionID_2)%owner

                      do r = 1, numberOfCartesiansOrbitals_2

                         if ( owner .eq. owner_2) then
                            WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)=&
                                 WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)*&
                                 ( 1 + MolecularSystem_getMass( speciesID ) / MolecularSystem_instance%species(speciesID)%particles(particleID)%mass  )

                            WaveFunction_instance( speciesID )%kineticMatrix%values(l,k)=&
                                 WaveFunction_instance( speciesID )%kineticMatrix%values(k,l)
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
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Matriz de energia cinetica: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
        call Matrix_show( WaveFunction_instance(speciesID)%kineticMatrix )
    end if

    !! Load N-Q- Attraction  Matrix
    arguments(1) = "ATTRACTION"

    WaveFunction_instance( speciesID )%puntualInteractionMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    

    !! Incluiding charge effect
    auxcharge = MolecularSystem_getCharge( speciesID )

    WaveFunction_instance( speciesID )%puntualInteractionMatrix%values = &
         WaveFunction_instance( speciesID )%puntualInteractionMatrix%values * (-auxCharge)

    !! DEBUG
    ! print *,"Matriz de interaccion n-quantum: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%puntualInteractionMatrix )

    !! Build Hcore Matrix
    if ( .not.allocated(WaveFunction_instance( speciesID )%HcoreMatrix%values ) ) then

       call Matrix_constructor( WaveFunction_instance( speciesID )%HcoreMatrix, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )

    end if

    WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%kineticMatrix%values + &
        WaveFunction_instance(speciesID)%puntualInteractionMatrix%values

    !! Add electric field F_i < \mu | e_i | \nu >
!    if ( sum(abs(CONTROL_instance%ELECTRIC_FIELD )) .ne. 0 ) then
    if ( CONTROL_instance%MULTIPOLE_ORDER /= 0 ) then

      write (*,"(T2,A15,I2)") "MULTIPOLE ORDER R**X: ", CONTROL_instance%MULTIPOLE_ORDER
      write (*,"(T2,A15,3F12.8)") "ELECTRIC FIELD:", CONTROL_instance%ELECTRIC_FIELD

      select case ( CONTROL_instance%MULTIPOLE_ORDER )

      case (1)

        arguments(1) = "MOMENTX0"
        WaveFunction_instance(speciesID)%electricField(1) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                             columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTY0"
        WaveFunction_instance(speciesID)%electricField(2) = Matrix_getFromFile(rows=totalNumberOfContractions, & 
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTZ0"
        WaveFunction_instance(speciesID)%electricField(3) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
  
        WaveFunction_instance(speciesID)%HCoreMatrix%values = &
          WaveFunction_instance(speciesID)%HCoreMatrix%values + &
                  auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(speciesID)%electricField(1)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(speciesID)%electricField(2)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(speciesID)%electricField(3)%values )

      case (2)

        arguments(1) = "MOMENTXX"
        WaveFunction_instance(speciesID)%electricField(4) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                             columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTYY"
        WaveFunction_instance(speciesID)%electricField(5) = Matrix_getFromFile(rows=totalNumberOfContractions, & 
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTZZ"
        WaveFunction_instance(speciesID)%electricField(6) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
  
        arguments(1) = "MOMENTXY"
        WaveFunction_instance(speciesID)%electricField(7) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                             columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTXZ"
        WaveFunction_instance(speciesID)%electricField(8) = Matrix_getFromFile(rows=totalNumberOfContractions, & 
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
        arguments(1) = "MOMENTYZ"
        WaveFunction_instance(speciesID)%electricField(9) = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                              columns=totalNumberOfContractions, &
                                                              unit=unit, binary=.true., arguments=arguments)    
  
  
        WaveFunction_instance(speciesID)%HCoreMatrix%values = &
          WaveFunction_instance(speciesID)%HCoreMatrix%values + &
                  (1.0/3.0)*auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(speciesID)%electricField(4)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(speciesID)%electricField(5)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(speciesID)%electricField(6)%values + &
           CONTROL_instance%ELECTRIC_FIELD(4)*WaveFunction_instance(speciesID)%electricField(7)%values + &
           CONTROL_instance%ELECTRIC_FIELD(5)*WaveFunction_instance(speciesID)%electricField(8)%values + &
           CONTROL_instance%ELECTRIC_FIELD(6)*WaveFunction_instance(speciesID)%electricField(9)%values )

      case default 

        call WaveFunction_exception ( ERROR, "WaveFunction_HCoreMatrix", "Multipole order not implemented")

      end select
    end if

     if ( CONTROL_instance%HARMONIC_CONSTANT /= 0.0_8 ) then 

      write (*,"(T2,A15,3F12.8)") "HARMONIC CONSTANT:", CONTROL_instance%HARMONIC_CONSTANT

      arguments(1) = "HARMONIC"
      WaveFunction_instance(speciesID)%harmonic = Matrix_getFromFile(rows=totalNumberOfContractions, &
                                                            columns=totalNumberOfContractions, &
                                                            unit=unit, binary=.true., arguments=arguments)   

      WaveFunction_instance(speciesID)%HCoreMatrix%values = &
        WaveFunction_instance(speciesID)%HCoreMatrix%values + &
        CONTROL_instance%HARMONIC_CONSTANT * WaveFunction_instance(speciesID)%harmonic%values

     end if

    close(34)    
         !WaveFunction_instance(speciesID)%externalPotentialMatrix%values 

    !! DEBUG
    !!   print *,"Matriz de hcore: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    !!   call Matrix_show( WaveFunction_instance( speciesID )%HcoreMatrix )

  end subroutine WaveFunction_HCoreMatrix

  !>
  !! @brief Ajusta la matriz de densidad para una especie espcificada
  subroutine WaveFunction_setDensityMatrix( densityMatrix, speciesID )
    implicit none

    type(Matrix), intent(in) :: densityMatrix
    integer :: speciesID

    integer :: totalNumberOfContractions
    character(50) :: densFile, labels(2)
    integer :: densUnit

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

    if( .not. allocated(WaveFunction_instance( speciesID )%densityMatrix%values )) then
       call Matrix_constructor( WaveFunction_instance( speciesID )%densityMatrix, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
    end if

    call Matrix_copyConstructor( WaveFunction_instance( speciesID )%densityMatrix, densityMatrix )

    !! Debug
    ! print*, "Matriz de densidad inicial ", MolecularSystem_getNameOfSpecie(speciesID)
    ! call Matrix_show(WaveFunction_instance( speciesID )%densityMatrix)

    !!Save this matrix for DFT calculations, because reasons
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       
       densUnit = 78
       densFile = trim(CONTROL_instance%INPUT_FILE)//trim(MolecularSystem_getNameOfSpecie(speciesID))//".densmatrix"
       open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")

       labels(1) = "DENSITY-MATRIX"
       labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
     
       call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )

       close (78)
    end if

  end subroutine WaveFunction_setDensityMatrix


  !>
  !! @brief Ajusta la matriz de densidad para una especie espcificada
  subroutine WaveFunction_setCoefficientsMatrix(coefficientsMatrix, speciesID )
    implicit none

    type(Matrix), intent(in) :: coefficientsMatrix
    integer :: speciesID

    integer :: totalNumberOfContractions

    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

    if( .not. allocated(WaveFunction_instance( speciesID )%coefficientsofcombination%values )) then
       call Matrix_constructor( WaveFunction_instance( speciesID )%coefficientsofcombination, &
            int(totalNumberOfContractions,8), int(totalNumberOfContractions,8), Math_NaN )
    end if

    call Matrix_copyConstructor( WaveFunction_instance( speciesID )%coefficientsofcombination, coefficientsMatrix )

    !! Debug
    ! print*, "Matriz de coefficients inicial ", MolecularSystem_getNameOfSpecie(speciesID)
    ! call Matrix_show(WaveFunction_instance( speciesID )%coefficientsofcombination)

  end subroutine WaveFunction_setCoefficientsMatrix



  !>
  !! @brief Calcula las componentes de energia para la especie especificada
  !! @warning Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  subroutine WaveFunction_obtainEnergyComponents(specieID )
    implicit none

    integer :: specieID
    real(8) :: auxCharge

    auxcharge = MolecularSystem_getCharge( specieID )

    !! Remove the electric field matrix to calculate the energy components
    if ( CONTROL_instance%MULTIPOLE_ORDER /= 0 ) then

      select case ( CONTROL_instance%MULTIPOLE_ORDER )

      case(1)
        WaveFunction_instance(specieID)%HCoreMatrix%values = &
          WaveFunction_instance(specieID)%HCoreMatrix%values - &
                  auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(specieID)%electricField(1)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(specieID)%electricField(2)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(specieID)%electricField(3)%values )

      case(2)

        WaveFunction_instance(specieID)%HCoreMatrix%values = &
          WaveFunction_instance(specieID)%HCoreMatrix%values - &
                  (1.0/3.0)*auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(specieID)%electricField(4)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(specieID)%electricField(5)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(specieID)%electricField(6)%values + &
           CONTROL_instance%ELECTRIC_FIELD(4)*WaveFunction_instance(specieID)%electricField(7)%values + &
           CONTROL_instance%ELECTRIC_FIELD(5)*WaveFunction_instance(specieID)%electricField(8)%values + &
           CONTROL_instance%ELECTRIC_FIELD(6)*WaveFunction_instance(specieID)%electricField(9)%values )

      end select

    end if




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
         WaveFunction_instance( specieID )%couplingEnergy + &
         WaveFunction_instance( specieID )%exchangeCorrelationEnergy


    !! Calcula la energia COSMO	

    if(CONTROL_instance%COSMO)then

       write(*,*)"COSMO energy contributions"

       write(*,*)"Especie = ",trim(MolecularSystem_instance%species(specieid)%name)

       WaveFunction_instance( specieID )%cosmoEnergy =  &
            0.5_8* (sum( transpose( wavefunction_instance( specieid )%densitymatrix%values ) * &
            wavefunction_instance( specieid )%cosmo1%values )) +0.5_8 *  &
            (sum( transpose( WaveFunction_instance( specieID )%densityMatrix%values ) * &
            WaveFunction_instance( specieID )%cosmo4%values )) + &
            0.5_8*(sum( transpose( WaveFunction_instance( specieID )%densityMatrix%values ) * &
            WaveFunction_instance( specieID )%cosmo2%values ) + &
            sum( transpose( WaveFunction_instance( specieID )%densityMatrix%values ) * &
            WaveFunction_instance( specieID )%cosmoCoupling%values )) 

       write(*,*)"COSMO energy 1 = ",0.5_8*(sum( transpose( wavefunction_instance( specieid )%densitymatrix%values ) * wavefunction_instance( specieid )%cosmo1%values )) 
       write(*,*)"COSMO energy 2 = ",0.5_8*(sum( transpose( wavefunction_instance( specieid )%densitymatrix%values ) * wavefunction_instance( specieid )%cosmo2%values )) 
       write(*,*)"COSMO energy 4 = ",0.5_8*(sum( transpose( wavefunction_instance( specieid )%densitymatrix%values ) * wavefunction_instance( specieid )%cosmo4%values )) 
       write(*,*)"COSMO coupling energy = ",0.5_8*(sum( transpose( wavefunction_instance( specieid )%densitymatrix%values ) * wavefunction_instance( specieid )%cosmoCoupling%values )) 
       write(*,*)"COSMO total energy = ",WaveFunction_instance( specieID )%cosmoEnergy


    end if

    !! Put back the electric field matrix to the Hcore matrix
    if ( CONTROL_instance%MULTIPOLE_ORDER /= 0 ) then

      select case ( CONTROL_instance%MULTIPOLE_ORDER )

      case(1)
        WaveFunction_instance(specieID)%HCoreMatrix%values = &
          WaveFunction_instance(specieID)%HCoreMatrix%values + &
                  auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(specieID)%electricField(1)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(specieID)%electricField(2)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(specieID)%electricField(3)%values )

      case(2)

        WaveFunction_instance(specieID)%HCoreMatrix%values = &
          WaveFunction_instance(specieID)%HCoreMatrix%values + &
                  (1.0/3.0)*auxcharge * &
          (CONTROL_instance%ELECTRIC_FIELD(1)*WaveFunction_instance(specieID)%electricField(4)%values + &
           CONTROL_instance%ELECTRIC_FIELD(2)*WaveFunction_instance(specieID)%electricField(5)%values + &
           CONTROL_instance%ELECTRIC_FIELD(3)*WaveFunction_instance(specieID)%electricField(6)%values + &
           CONTROL_instance%ELECTRIC_FIELD(4)*WaveFunction_instance(specieID)%electricField(7)%values + &
           CONTROL_instance%ELECTRIC_FIELD(5)*WaveFunction_instance(specieID)%electricField(8)%values + &
           CONTROL_instance%ELECTRIC_FIELD(6)*WaveFunction_instance(specieID)%electricField(9)%values )

      end select

    end if




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
    !
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

  !!cosmo matrices construction

  subroutine WaveFunction_cosmoHCoreMatrix(file, speciesID)
    implicit none

    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    ! integer :: k, l, r, s
    ! integer :: ParticleID, ParticleID_2
    ! integer :: contractionID, contractionID_2
    ! integer :: numberOfCartesiansOrbitals, numberOfCartesiansOrbitals_2
    ! integer :: owner, owner_2
    ! integer :: auxCharge
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(10) :: arguments(2)

    !! Open file
    unit = 44
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))  


    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          


    !! Load electron potential vs clasical charges cosmo matrix
    arguments(1) = "COSMO1"    
    WaveFunction_instance( speciesID )%cosmo1 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)



    ! DEBUG
    ! print *,"Matriz cosmo1: ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%cosmo1 )

    !! Load clasical potential vs quantum charges cosmo matrix
    arguments(1) = "COSMO4"


    WaveFunction_instance( speciesID )%cosmo4 = Matrix_getFromFile(rows=totalNumberOfContractions, columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments)    


    !! DEBUG
    ! print *,"Matriz cosmo 4 ", trim(MolecularSystem_getNameOfSpecie(speciesID))
    ! call Matrix_show( WaveFunction_instance(speciesID)%cosmo4 )

    close(44)    


  end subroutine WaveFunction_cosmoHcoreMatrix

  !<
  !! @brief Contruye una matriz de interaccion con un potencial externo
  !!
  !! @param nameOfSpecie nombre de la especie seleccionada.
  !>
  subroutine WaveFunction_buildExternalPotentialMatrix( file, speciesID )
    implicit none
    character(30) :: nameOfSpecieSelected
    character(*), intent(in) :: file
    integer, intent(in) :: speciesID

    integer :: unit
    integer :: numberOfContractions
    integer :: totalNumberOfContractions
    character(50) :: arguments(2)

    arguments(1) = "EXTERNAL_POTENTIAL"
    arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))

    !! Open file
    unit = 34
    open(unit = unit, file=trim(file), status="old", form="unformatted")

    WaveFunction_instance(speciesID)%externalPotentialMatrix%values = 0.0_8

    !! Get number of shells and number of cartesian contractions
    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )          
    WaveFunction_instance( speciesID )%externalPotentialMatrix = Matrix_getFromFile(rows=totalNumberOfContractions, &
         columns=totalNumberOfContractions, &
         unit=unit, binary=.true., arguments=arguments(1:2))
    close(34)

      !!if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )  then  !! Numerical integration
      !!if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
      !!  WaveFunction_instance(specieID)%externalPotentialMatrix = &
      !!    IntegralManager_getNumericalInteractionWithPotentialMatrix( &
      !!    ExternalPotential_Manager_instance%externalsPots, specieID, integralName="external" )

      !!else     !! From xml file
      !!  WaveFunction_instance(specieID)%externalPotentialMatrix = &
      !!    IntegralManager_getNumericalPotentialMatrixFromXml( &
      !!    ExternalPotential_Manager_instance%externalsPots, specieID, integralName="external" )
      !!end if
      !!else    !! Analytical Integration  

      !!end if

    !! Debug
    if (  CONTROL_instance%DEBUG_SCFS) then
      print *,"EXTERNAL POTENTIAL MATRIX FOR: ", arguments(2)
      call Matrix_show(WaveFunction_instance(speciesID)%externalPotentialMatrix)
    end if

  end subroutine WaveFunction_buildExternalPotentialMatrix


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
