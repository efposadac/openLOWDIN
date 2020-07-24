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
!! @brief Modulo para realizacion de iteraciones SCF sobre diferentes especies cuanticas
!!
!!  Este modulo define una seudoclase para realizacion de ciclos SCF sobre todas las
!! especies cuanticas.
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2008-08-30
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y metodos
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe el modulo para su inclusion en Lowdin
module MultiSCF_
  use Exception_
  use Stopwatch_
  use CONTROL_
  use Matrix_
  use List_  
  use MolecularSystem_
  use WaveFunction_
  use SingleSCF_

  implicit none    

  !< enum Convergece_status_type {
  ! integer, parameter :: SCF_GLOBAL_CONVERGENCE_FAILED = 0
  ! integer, parameter :: SCF_GLOBAL_CONVERGENCE_CONTINUE = 1
  ! integer, parameter :: SCF_GLOBAL_CONVERGENCE_SUCCESS = 2
  !< }

  !< For effective Fock matrix for restricted open-shell SCF
  logical :: isROHF
  !>

  type, public :: MultiSCF

     type(List) :: energyOMNE
     character(50) :: name
     integer :: numberOfIterations
     integer :: status
     real(8), allocatable :: singleEnergyTolerance(:)
     real(8), allocatable :: singleDensityTolerance(:)
     integer, allocatable :: singleMaxIterations(:)

     !! Global energies
     real(8) :: totalEnergy
     real(8) :: totalCouplingEnergy
     real(8) :: electronicRepulsionEnergy

     !! Global Density Standard Deviation
     real(8) :: totalDensityMatrixStandardDeviation
     
     !! Cosmo
     real(8) :: cosmo3Energy


  end type MultiSCF

  type(MultiSCF), public, target :: MultiSCF_instance

  public :: &
       MultiSCF_constructor, &
       MultiSCF_destructor, &
       MultiSCF_getNumberOfIterations, &
       MultiSCF_getLastEnergy, &
       MultiSCF_iterate, &
       MultiSCF_restart, &
       MultiSCF_reset

contains


  !>
  !! @brief Define el constructor para la clase
  subroutine MultiSCF_constructor(iterationScheme)
    implicit none
    integer :: iterationScheme
    integer :: i
    
    isROHF = .false.

    select case(iterationScheme)
    case(0)
       MultiSCF_instance%name="NONELECTRONIC_FULLY_CONVERGED_BY_ELECTRONIC_ITERATION"
    case(1)
       MultiSCF_instance%name="ELECTRONIC_FULLY_CONVERGED_BY_NONELECTRONIC_ITERATION"
    case(2)
       MultiSCF_instance%name="SPECIES_FULLY_CONVERGED_INDIVIDUALLY"
    case(3)
       MultiSCF_instance%name="SIMULTANEOUS"
    case(4)
       call MultiSCF_exception( ERROR, "selected an iteration scheme not implemented (>3)", "at SCF program, MultiSCF_constructor")
    end select
       
    call List_constructor( MultiSCF_instance%energyOMNE,"ENERGY", CONTROL_instance%LISTS_SIZE)
    MultiSCF_instance%numberOfIterations = 0
    MultiSCF_instance%status = 0

    allocate(MultiSCF_instance%singleEnergyTolerance(MolecularSystem_getNumberOfQuantumSpecies()),&
         MultiSCF_instance%singleDensityTolerance(MolecularSystem_getNumberOfQuantumSpecies()),&
         MultiSCF_instance%singleMaxIterations(MolecularSystem_getNumberOfQuantumSpecies()))

    do i = 1, MolecularSystem_getNumberOfQuantumSpecies()
       if(MolecularSystem_instance%species(i)%isElectron ) then
          MultiSCF_instance%singleEnergyTolerance(i)=CONTROL_instance%ELECTRONIC_ENERGY_TOLERANCE
          MultiSCF_instance%singleDensityTolerance(i)=CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
       else
          MultiSCF_instance%singleEnergyTolerance(i)=CONTROL_instance%NONELECTRONIC_ENERGY_TOLERANCE
          MultiSCF_instance%singleDensityTolerance(i)=CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE
       end if

       !The maximum number of iterations for each species in selected according to the SCF scheme chosen
       select case(iterationScheme)

       case( 0 )          
          ! we perform single species iterations for nonelectrons
          if(MolecularSystem_instance%species(i)%isElectron ) then
             MultiSCF_instance%singleMaxIterations(i)=1
          else
             MultiSCF_instance%singleMaxIterations(i)=CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
          end if
          
       case( 1 )
          ! we perform single species iterations for nelectrons
          if(MolecularSystem_instance%species(i)%isElectron ) then
             MultiSCF_instance%singleMaxIterations(i)=CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
          else
             MultiSCF_instance%singleMaxIterations(i)=1
          end if


       case( 2 )
          ! we perform single species  for all species
          if(MolecularSystem_instance%species(i)%isElectron ) then
             MultiSCF_instance%singleMaxIterations(i)=CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS
          else
             MultiSCF_instance%singleMaxIterations(i)=CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS
          end if
          
       case ( 3 )
          ! we do not perform single species SCF
          if(MolecularSystem_instance%species(i)%isElectron ) then
             MultiSCF_instance%singleMaxIterations(i)=1
          else
             MultiSCF_instance%singleMaxIterations(i)=1
          end if
       end select

    end do
    
    MultiSCF_instance%totalEnergy = 0.0_8
    MultiSCF_instance%cosmo3Energy = 0.0_8
    MultiSCF_instance%totalCouplingEnergy = 0.0_8
    MultiSCF_instance%electronicRepulsionEnergy = 0.0_8

  end subroutine MultiSCF_constructor

  !>
  !! @brief Define el destructor para clase
  subroutine MultiSCF_destructor()
    implicit none

    call List_destructor( MultiSCF_instance%energyOMNE )

  end subroutine MultiSCF_destructor

  !>
  !! @brief Retorna la energia obtenida en la ultima iteracion iter-especies
  function MultiSCF_getNumberOfIterations() result( output )
    implicit none

    integer :: output

    output= MultiSCF_instance%numberOfIterations

  end function MultiSCF_getNumberOfIterations

  !>
  !! @brief Retorna la energia obtenida en la ultima iteracion iter-especies
  function MultiSCF_getLastEnergy( ) result( output )
    implicit none

    real(8) :: output

    call List_end( MultiSCF_instance%energyOMNE )
    output= List_current( MultiSCF_instance%energyOMNE )

  end function MultiSCF_getLastEnergy

  !>
  !! @brief Realiza esquema de iteracion SCF para todas las especies cuanticas presentes
  !! @param
  !!
  subroutine MultiSCF_iterate(iterationScheme)
    implicit none
    integer, intent(in) :: iterationScheme

    integer :: i
    integer :: numberOfSpecies
    character(30) :: nameOfSpecies
    integer :: speciesID
    character(50) :: densFile
    real(8) :: singleEnergyTolerance, singleDensityTolerance
    integer :: singleMaxIterations, singleIterator
    real(8) :: oldEnergy
    real(8) :: deltaEnergy

    
    
    !!!We start with an update of the global energy and matrices
    
    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    densFile="lowdin.densmatrix"

    !!DFT calculations are performed only in global iterations
    !!Save density matrices to file for DFT calculations
    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       call WaveFunction_writeDensityMatricesToFile(trim(densFile))
       call system("lowdin-DFT.x SCF_DFT "//trim(densFile))
    end if

    !Build Fock matrices and calculate energy with the initial densities
    !Coupling Matrix is only updated in global SCF cycles
    do i = 1, numberOfSpecies
       nameOfSpecies = MolecularSystem_getNameOfSpecie(i)

       call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies))

       call WaveFunction_buildCouplingMatrix( trim(nameOfSpecies))

       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          call WaveFunction_readExchangeCorrelationMatrix( trim(nameOfSpecies))
       end if

       if (CONTROL_instance%COSMO) then
          call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecies))
          call WaveFunction_buildCosmoCoupling( trim(nameOfSpecies) )
       end if

       call WaveFunction_buildFockMatrix( trim(nameOfSpecies) )

       call WaveFunction_obtainTotalEnergyForSpecie( nameOfSpecies )

    end do

    !Get total energy
    call WaveFunction_obtainTotalEnergy(&
         MultiSCF_instance%totalEnergy, &
         MultiSCF_instance%totalCouplingEnergy, &
         MultiSCF_instance%electronicRepulsionEnergy, &
         MultiSCF_instance%cosmo3Energy)

    call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)             
    MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1

!!!Now we procede to update each species density matrices according to the iteration scheme selected
    MultiSCF_instance%totalDensityMatrixStandardDeviation=0.0
    do i = 1, numberOfSpecies
       nameOfSpecies = MolecularSystem_getNameOfSpecie(i)
       oldEnergy=WaveFunction_instance(i)%totalEnergyForSpecie
       deltaEnergy=1.0E16_8
       singleIterator=0

       !The maximum number of iterations for each species in selected according to the SCF scheme chosen
       
       do while( (abs(deltaEnergy) .gt. MultiSCF_instance%singleEnergyTolerance(i) .or. &
            Matrix_standardDeviation( WaveFunction_instance(i)%beforeDensityMatrix, WaveFunction_instance(i)%densityMatrix ) .gt. MultiSCF_instance%singleDensityTolerance(i)) .and. &
            singleIterator .lt. MultiSCF_instance%singleMaxIterations(i) )

          singleIterator=singleIterator+1
          !! Obtains the new eigenvectors
          call SingleSCF_iterate( trim(nameOfSpecies))

          !! Determina la desviacion estandar de los elementos de la matriz de densidad
          call Matrix_copyConstructor( WaveFunction_instance(i)%beforeDensityMatrix, WaveFunction_instance(i)%densityMatrix )
          call WaveFunction_buildDensityMatrix( trim(nameOfSpecies) )
          call List_push_back( WaveFunction_instance(i)%standardDesviationOfDensityMatrixElements, &
               Matrix_standardDeviation( WaveFunction_instance(i)%beforeDensityMatrix, WaveFunction_instance(i)%densityMatrix ) )

          !! Calcula la energy de la especie con la nueva densidad
          call WaveFunction_obtainTotalEnergyForSpecie( trim(nameOfSpecies) )
          call List_push_back( WaveFunction_instance(i)%energySCF, WaveFunction_instance(i)%totalEnergyForSpecie )
          call List_push_back( WaveFunction_instance(i)%diisError, Convergence_getDiisError( WaveFunction_instance(i)%convergenceMethod) )
          
          if(MultiSCF_instance%singleMaxIterations(i).gt.1) then
             !! Updates two particle matrix
             call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecies))

             if (CONTROL_instance%COSMO) then
                call WaveFunction_buildCosmo2Matrix(trim(nameOfSpecies))
             end if
             !!Builds new fock Matrix
             call WaveFunction_buildFockMatrix( trim(nameOfSpecies))

             !!Checks energy convergence
             deltaEnergy = oldEnergy -WaveFunction_instance(i)%totalEnergyForSpecie
             oldEnergy = WaveFunction_instance(i)%totalEnergyForSpecie

             !!Prints iteration results
             if (  CONTROL_instance%DEBUG_SCFS) then
                write (*,"(A10,I5,F20.12,F20.12,F20.12)") trim(nameOfSpecies), singleIterator , &
                     WaveFunction_instance(i)%totalEnergyForSpecie, deltaEnergy, &
                     List_current(WaveFunction_instance(i)%standardDesviationOfDensityMatrixElements)
             end if

             !!Prints convergence messages
             if(singleIterator .ge. MultiSCF_instance%singleMaxIterations(i) ) &
                  write (*,"(T35,A10,A30,I4,A)") trim(nameOfSpecies), " Max. subcycles reached(", MultiSCF_instance%singleMaxIterations(i),")"
             
             if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "DENSITY" .and. &
                  Matrix_standardDeviation( WaveFunction_instance(i)%beforeDensityMatrix, WaveFunction_instance(i)%densityMatrix ) .lt. MultiSCF_instance%singleDensityTolerance(i)) & 
                  write (*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Density converged in", singleIterator ," subcycles"
             
             if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "ENERGY" .and. &
                  abs(deltaEnergy) .lt. MultiSCF_instance%singleEnergyTolerance(i) )&
                  write (*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Energy converged in", singleIterator ," subcycles"
             
             if(trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) .eq. "BOTH" .and. & 
                  Matrix_standardDeviation( WaveFunction_instance(i)%beforeDensityMatrix, WaveFunction_instance(i)%densityMatrix ) .lt. MultiSCF_instance%singleDensityTolerance(i) .and. & 
                  abs(deltaEnergy) .lt. MultiSCF_instance%singleEnergyTolerance(i) )&
                  write (*,"(T35,A10,A30,I4,A)")  trim(nameOfSpecies), " Energy-density converged in", singleIterator ," subcycles"
          end if

       end do
       MultiSCF_instance%totalDensityMatrixStandardDeviation=MultiSCF_instance%totalDensityMatrixStandardDeviation+&
            sqrt(List_current(WaveFunction_instance(i)%standardDesviationOfDensityMatrixElements)**2)
    end do

  end subroutine MultiSCF_iterate

  ! >
  ! @brief Prueba si la energia del la ultima iteracion a sufrido un cambio por debajo de cierta tolerancia
  function MultiSCF_checkConvergence(  tolerace ) result( output )
    implicit none
    real(8), intent(in) :: tolerace
    integer :: output

    character(30) :: nameOfSpecie
    real(8) :: deltaEnergy
    real(8) :: finalEnergy
    real(8) :: toleraceOfSpecie
    type(Exception) :: ex
    integer :: speciesID
    logical :: auxVar

    ! if ( MultiSCF_instance%numberOfIterations > CONTROL_instance%SCF_GLOBAL_MAX_ITERATIONS ) then

    !    output = SCF_GLOBAL_CONVERGENCE_SUCCESS
    !    call List_end( MultiSCF_instance%energyOMNE )
    !    finalEnergy= List_current( MultiSCF_instance%energyOMNE )

    !    !! Obtiene el valor  de energia anterior a la ultima iteracion
    !    call List_iterate( MultiSCF_instance%energyOMNE, -1 )

    !    !! Obtiene el cambio de energia en las ultimas dos iteraciones
    !    deltaEnergy = finalEnergy - List_current( MultiSCF_instance%energyOMNE )
    !    print *,""
    !    write (6,"(A20,A20)") "   Current Energy   ", "   Current Change   "
    !    print *,"-----------------------------------------------"
    !    write (6,"(F20.10,F20.10)") finalEnergy, deltaEnergy

    !    print *,  "The number of Iterations was exceded, the convergence had failed"

    ! else

    !    if ( MultiSCF_instance%numberOfIterations > 1 ) then

    !       auxVar=.true.
    !       do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

    !          nameOfSpecie = MolecularSystem_getNameOfSpecie(speciesID)

    !          toleraceOfSpecie = MultiSCF_instance%electronicTolerance

    !          if (.not. MolecularSystem_instance%species(speciesID)%isElectron ) then
    !             toleraceOfSpecie = MultiSCF_instance%nonelectronicTolerance
    !          end if

    !          if ( SingleSCF_testDensityMatrixChange( nameOfSpecie, &
    !               toleraceOfSpecie ) == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then
    !             auxVar = auxVar .and. .true.
    !          else
    !             auxVar = auxVar .and. .false.
    !          end if

    !       end do

    !       call List_end( MultiSCF_instance%energyOMNE )
    !       deltaEnergy= List_current( MultiSCF_instance%energyOMNE )

    !       !! Obtiene el valor  de energia anterior a la ultima iteracion
    !       call List_iterate( MultiSCF_instance%energyOMNE, -1 )

    !       !! Obtiene el cambio de energia en las ultimas dos iteraciones
    !       deltaEnergy = deltaEnergy - List_current( MultiSCF_instance%energyOMNE )

    !       if( ( ( abs( deltaEnergy ) < tolerace ) .and. auxVar ) .or. abs(deltaEnergy) < CONTROL_instance%TOTAL_ENERGY_TOLERANCE ) then
    !          output = SCF_GLOBAL_CONVERGENCE_SUCCESS
    !       else
    !          output = SCF_GLOBAL_CONVERGENCE_CONTINUE
    !       end if

    !    else

    !       output = SCF_GLOBAL_CONVERGENCE_CONTINUE

    !    end if

    ! end if

  end function MultiSCF_checkConvergence


  
  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine MultiSCF_restart()
    implicit none

    MultiSCF_instance%numberOfIterations = 0

    call List_clear( MultiSCF_instance%energyOMNE )

  end subroutine MultiSCF_restart

  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine MultiSCF_reset()
    implicit none

    integer :: speciesIterator

    MultiSCF_instance%numberOfIterations = 0

    call List_clear( MultiSCF_instance%energyOMNE )
    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

    do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies()
       call SingleSCF_reset(speciesIterator)
    end do

  end subroutine MultiSCF_reset

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine MultiSCF_exception( typeMessage, description, debugDescription)
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

  end subroutine MultiSCF_exception

end module MultiSCF_

  ! !>
  ! !! @brief Realiza la convergencia completa de la parte no electronica por cada iteracion
  ! !!		SCF de la parte electronica
  ! !! @todo falta implementacion completa
  ! subroutine MultiSCF_iterateNonElectronicFullyByElectronicIteration()
  !   implicit none

  !   integer :: iteratorOfSpecie
  !   integer :: iteratorOfElectronicSpecie
  !   character(30) :: nameOfSpecie
  !   character(30) :: nameOfElectronicSpecie
  !   character(30) :: nameOfInitialSpecie
  !   real(8) :: tolerace
  !   integer :: statusSystem
  !   integer :: densUnit
  !   character(50) :: densFile
  !   character(30) :: labels(2)
  !   integer :: speciesID

  !   MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !   iteratorOfSpecie = 1
  !   nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)    

  !   if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !      species_loop : do

  !         if(iteratorOfSpecie <= MolecularSystem_getNumberOfQuantumSpecies()) then          

  !            nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)   

  !            if( .not. MolecularSystem_instance%species(iteratorOfSpecie)%isElectron) then

  !               iteratorOfSpecie = iteratorOfSpecie + 1 
  !               cycle species_loop

  !            else 

  !               non_electronic_loop : do iteratorOfElectronicSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

  !                  if( .not. MolecularSystem_instance%species(iteratorOfElectronicSpecie)%isElectron ) then

  !                     MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !                     nameOfElectronicSpecie = MolecularSystem_getNameOfSpecie(iteratorOfElectronicSpecie)                   
  !                     tolerace = MultiSCF_instance%electronicTolerance

  !                     do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !                          (SingleSCF_getNumberOfIterations(iteratorOfElectronicSpecie) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

  !                        call WaveFunction_buildTwoParticlesMatrix( trim(nameOfElectronicSpecie))

  !                        if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then

  !                           !!Save density matrices to file for DFT calculations
  !                           densUnit = 78
  !                           densFile = trim(CONTROL_instance%INPUT_FILE)//"densmatrix"
  !                           open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
  !                           labels(1) = "DENSITY-MATRIX"
  !                           do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
  !                              labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
  !                              call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
  !                           end do
  !                           close (densUnit)

  !                           statusSystem = system ("lowdin-DFT.x SCF_DFT "//trim(densFile))
  !                           call WaveFunction_readExchangeCorrelationMatrix( trim(nameOfSpecie))

  !                        end if

                         
  !                        if (CONTROL_instance%COSMO) then
  !                           call  WaveFunction_buildCosmo2Matrix( trim(nameOfElectronicSpecie))
  !                           if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
  !                              call WaveFunction_buildCosmoCoupling( trim(nameOfElectronicSpecie) )
  !                           end if
  !                        end if

  !                        !! At first iteration is not included the coupling operator.
  !                        ! if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
  !                        call WaveFunction_buildCouplingMatrix( trim(nameOfElectronicSpecie) )
  !                        ! end if

  !                        !! Build Fock Matrix
  !                        call WaveFunction_buildFockMatrix( trim(nameOfElectronicSpecie) )

  !                        !! Perform SCF iteration for Single species
  !                        call SingleSCF_iterate( trim(nameOfElectronicSpecie) )

  !                        !! Test energy or density change (as requested)
  !                        if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

  !                           MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, tolerace )

  !                           MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, &
  !                                MultiSCF_instance%electronicTolerance )

  !                        else

  !                           MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfElectronicSpecie, &
  !                                MultiSCF_instance%electronicTolerance )

  !                        end if

  !                     end do !! end convergence for electrons

  !                  end if

  !               end do non_electronic_loop

  !               !! Realiza iteracion SCF para una especie cuantica particular
  !               call SingleSCF_iterate( trim( nameOfSpecie ))

  !               if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

  !                  MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
  !                       MultiSCF_instance%nonelectronicTolerance )
  !               else

  !                  MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
  !                       MultiSCF_instance%nonelectronicTolerance )

  !               end if

  !               if ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !                    (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS ) ) then

  !                  cycle species_loop

  !               else if ( MultiSCF_instance%status == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then

  !                  iteratorOfSpecie = iteratorOfSpecie + 1 
  !                  cycle species_loop                      

  !               end if

  !            end if

  !         else

  !            !! Finaliza un ciclo SCF completo para para todas la especies presentes
  !            call WaveFunction_obtainTotalEnergy(&
  !                 MultiSCF_instance%totalEnergy, &
  !                 MultiSCF_instance%totalCouplingEnergy, &
  !                 MultiSCF_instance%electronicRepulsionEnergy, &
  !                 MultiSCF_instance%cosmo3Energy )

  !            call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
  !            MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1              
  !            exit species_loop

  !         end if

  !      end do species_loop

  !   else

  !      if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
  !         call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
  !      end if

  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy )

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

  !   end if

  ! end subroutine MultiSCF_iterateNonElectronicFullyByElectronicIteration

  ! !>
  ! !! @brief Realiza la convergencia completa de la parte electronica por cada iteracion
  ! !!	    SCF de la parte no electronica
  ! subroutine MultiSCF_iterateElectronicFullyByNonElectronicIteration( )
  !   implicit none
  !   integer :: iteratorOfSpecie
  !   integer :: iteratorOfElectronicSpecie
  !   character(30) :: nameOfSpecie
  !   character(30) :: nameOfElectronicSpecie
  !   character(30) :: nameOfInitialSpecie
  !   real(8) :: tolerace
  !   integer :: densUnit
  !   character(50) :: densFile
  !   character(30) :: labels(2)


  !   MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !   iteratorOfSpecie = 1
  !   nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)    

  !   if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !      species_loop : do

  !         if(iteratorOfSpecie <= MolecularSystem_getNumberOfQuantumSpecies()) then          

  !            nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)   

  !            if(MolecularSystem_instance%species(iteratorOfSpecie)%isElectron) then

  !               iteratorOfSpecie = iteratorOfSpecie + 1 
  !               cycle species_loop

  !            else 

  !               electronic_loop : do iteratorOfElectronicSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

  !                  if(MolecularSystem_instance%species(iteratorOfElectronicSpecie)%isElectron) then

  !                     MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !                     nameOfElectronicSpecie = MolecularSystem_getNameOfSpecie(iteratorOfElectronicSpecie)                   
  !                     tolerace = MultiSCF_instance%electronicTolerance

  !                     do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !                          (SingleSCF_getNumberOfIterations(iteratorOfElectronicSpecie) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

  !                        call WaveFunction_buildTwoParticlesMatrix( trim(nameOfElectronicSpecie) )

  !                        !! At first iteration is not included the coupling operator.
  !                         if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
  !                           call WaveFunction_buildCouplingMatrix( trim(nameOfElectronicSpecie) )
  !                         end if

  !                        if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
  !                           call WaveFunction_readExchangeCorrelationMatrix( trim(nameOfSpecie))
  !                        end if
                         
  !                        ! call Matrix_show(wavefunction_instance(iteratorOfSpecie)%couplingMatrix)


  !                        if (CONTROL_instance%COSMO) then
  !                           call  WaveFunction_buildCosmo2Matrix( trim(nameOfElectronicSpecie))
  !                           if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
  !                              call WaveFunction_buildCosmoCoupling( trim(nameOfElectronicSpecie) )
  !                           end if
  !                        end if



  !                        !! Build Fock Matrix
  !                        call WaveFunction_buildFockMatrix( trim(nameOfElectronicSpecie) )

  !                        !! Perform SCF iteration for Single species
  !                        call SingleSCF_iterate( trim(nameOfElectronicSpecie) )

  !                        !! Test energy or density change (as requested)
  !                        if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

  !                           MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, tolerace )

  !                           MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, &
  !                                MultiSCF_instance%electronicTolerance )

  !                        else

  !                           MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfElectronicSpecie, &
  !                                MultiSCF_instance%electronicTolerance )

  !                        end if

  !                     end do !! end convergence for electrons

  !                  end if

  !               end do electronic_loop

  !               !! Realiza iteracion SCF para una especie cuantica particular
  !               call SingleSCF_iterate( trim( nameOfSpecie ))

  !               if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

  !                  MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
  !                       MultiSCF_instance%nonelectronicTolerance )
  !               else

  !                  MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
  !                       MultiSCF_instance%nonelectronicTolerance )

  !               end if

  !               if ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !                    (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS ) ) then

  !                  cycle species_loop

  !               else if ( MultiSCF_instance%status == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then

  !                  iteratorOfSpecie = iteratorOfSpecie + 1 
  !                  cycle species_loop                      

  !               end if

  !            end if

  !         else

  !            !! Finaliza un ciclo SCF completo para para todas la especies presentes
  !            call WaveFunction_obtainTotalEnergy(&
  !                 MultiSCF_instance%totalEnergy, &
  !                 MultiSCF_instance%totalCouplingEnergy, &
  !                 MultiSCF_instance%electronicRepulsionEnergy, &
  !                 MultiSCF_instance%cosmo3Energy)

  !            call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
  !            MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1              
  !            exit species_loop

  !         end if

  !      end do species_loop

  !   else

  !      if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
  !         call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
  !      end if

  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy)

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

  !   end if

  ! end subroutine MultiSCF_iterateElectronicFullyByNonElectronicIteration

  ! !>
  ! !! @brief Realiza la convergencia completa de cada especie de manera individual
  ! subroutine MultiSCF_iterateSpecieFullyConvergedIndividually()
  !   implicit none
  !   integer :: iteratorOfSpecie
  !   character(30) :: nameOfSpecie
  !   character(30) :: nameOfInitialSpecie
  !   real(8) :: tolerace

  !   MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !   iteratorOfSpecie = 1
  !   nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)

  !   if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !      do iteratorOfSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

  !         nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)

  !         if ( MolecularSystem_instance%species(iteratorOfSpecie)%isElectron ) then
  !            tolerace = MultiSCF_instance%electronicTolerance
  !         else
  !            tolerace = MultiSCF_instance%nonelectronicTolerance
  !         end if

  !         do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !              (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

  !            !             print*, "Multi-SCF", nameOfSpecie
  !            !             call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

  !            if (CONTROL_instance%COSMO) then
  !               call  WaveFunction_buildCosmo2Matrix( trim(nameOfSpecie))
  !               if(SingleSCF_getNumberOfIterations( iteratorOfSpecie ) > 0) then
  !                  call WaveFunction_buildCosmoCoupling( trim(nameOfSpecie) )
  !               end if
  !            end if

  !            !! At first iteration is not included the coupling operator.
  !            if(SingleSCF_getNumberOfIterations( iteratorOfSpecie ) > 0) then
  !               call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie) )
  !            end if

  !            !! Build Fock Matrix
  !            call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )

  !            !! Perform SCF iteration for Single species

  !            call SingleSCF_iterate( trim(nameOfSpecie))

  !            !! Test energy or density change (as requested)
  !            if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "DENSITY" ) then

  !               MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, tolerace )

  !               MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
  !                    MultiSCF_instance%electronicTolerance )

  !            else

  !               MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
  !                    MultiSCF_instance%electronicTolerance )

  !            end if

  !         end do

  !         MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

  !      end do

  !      !! Finaliza un ciclo SCF completo para para todas la especies presentes
  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy)

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
  !      MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1       

  !   else

  !      if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
  !         call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
  !      end if

  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy)

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

  !   end if

  ! end subroutine MultiSCF_iterateSpecieFullyConvergedIndividually


  ! subroutine MultiSCF_iterateSimultaneous()
  !   implicit none

  !   integer :: i
  !   integer :: numberOfSpecies
  !   character(30) :: nameOfSpecie
  !   character(30) :: nameOfInitialSpecie
  !   real(8) :: startTime, endTime
  !   real(8) :: time1,time2
  !   logical :: auxValue
  !   integer :: statusSystem
  !   integer :: speciesID
  !   integer :: densUnit
  !   character(50) :: densFile
  !   character(30) :: labels(2)

  !   MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE
  !   numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

  !   if( numberOfSpecies > 1 ) then
       
  !      if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
  !         !!Save density matrices to file for DFT calculations
  !         densUnit = 78
  !         densFile = trim(CONTROL_instance%INPUT_FILE)//"densmatrix"
  !         open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
  !         labels(1) = "DENSITY-MATRIX"
  !         do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
  !            labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
  !            call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
  !         end do
  !         close (densUnit)
  !         call system ("lowdin-DFT.x SCF_DFT "//trim(densFile))
  !      end if

  !      do i = 1, numberOfSpecies
          
  !         nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
  !         ! call SingleSCF_iterate( trim(nameOfSpecie), actualizeDensityMatrix=.false.)
  !      end do

  !      !! Get effective Fock matrix for restricted open-shell SCF
  !      ! if( isROHF ) then
  !      ! call MultiSCF_instance%mergeFockMatrix()
  !      ! end if

  !      do i=1, numberOfSpecies

  !         nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
  !         ! call SingleSCF_actualizeDensityMatrix( trim(nameOfSpecie) )

  !      end do

  !      !Call DFT_actualizeExchangePotentialContributions FELIX
       
  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy)

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)             
  !      MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1
  !      return

  !   else
       
  !      !!Especifica el procedimiento a seguir si lo electrones han sido congelados
  !      if ( .not. CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
  !         call MultiSCF_iterateUniqueSpecie( speciesID = numberOfSpecies )
  !      end if

  !      call WaveFunction_obtainTotalEnergy(&
  !           MultiSCF_instance%totalEnergy, &
  !           MultiSCF_instance%totalCouplingEnergy, &
  !           MultiSCF_instance%electronicRepulsionEnergy, &
  !           MultiSCF_instance%cosmo3Energy)

  !      call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)             

  !   end if

  ! end subroutine MultiSCF_iterateSimultaneous

  ! !>
  ! !! @brief Realiza la convergencia completa de la unica especie presente
  ! subroutine MultiSCF_iterateUniqueSpecie( speciesID )
  !   implicit none

  !   integer :: speciesID
 
  !   character(30) :: nameOfSpecie    
  !   real(8) :: tolerace
  !   real(8) :: diisError
  !   character :: typeConvergence
  !   integer :: statusSystem
  !   integer :: densUnit
  !   character(50) :: densFile
  !   character(30) :: labels(2)
    
  !   nameOfSpecie = MolecularSystem_getNameOfSpecie(speciesID)

  !   if ( MolecularSystem_instance%species(speciesID)%isElectron ) then
  !      tolerace = MultiSCF_instance%electronicTolerance
  !   else
  !      tolerace = MultiSCF_instance%nonelectronicTolerance
  !   end if

  !   !! Realiza ciclo de convergencia SCF para la unica especie cuantica presente

  !   if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
  !      print *,""
  !      print *,"Begin SCF calculation by: ",trim(nameOfSpecie)
  !      print *,"-------------------------"
  !      print *,""
  !      print *,"-----------------------------------------------------------------"
  !      write (6,"(A10,A12,A25,A20)") "Iteration", "Energy", " Density Change","         DIIS Error "
  !      print *,"-----------------------------------------------------------------"
  !   end if

  !   !! Build an initial two particles matrix, which it will be recalculated in SingleSCF_iterate
  !   call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie))

    
  !   ! write(*,*)"entre al unique specie"
  !   do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
  !        ( SingleSCF_getNumberOfIterations(speciesID) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

  !      !!      This is not necessary, the two particles matrix can be the one calculated in the previous iteration. 
  !      !!       call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

  !      if (CONTROL_instance%COSMO) then
  !         call  WaveFunction_buildCosmo2Matrix( trim(nameOfSpecie))
  !         call WaveFunction_buildCosmoCoupling( trim(nameOfSpecie) )
  !      end if

  !      if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
  !         !!Save density matrices to file for DFT calculations
  !         densUnit = 78
  !         densFile = trim(CONTROL_instance%INPUT_FILE)//"densmatrix"
  !         open(unit = densUnit, file=trim(densFile), status="replace", form="unformatted")
  !         labels(1) = "DENSITY-MATRIX"
  !         do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
  !            labels(2) = MolecularSystem_getNameOfSpecie(speciesID)
  !            call Matrix_writeToFile(WaveFunction_instance(speciesID)%densityMatrix, unit=densUnit, binary=.true., arguments = labels )
  !         end do
  !         close (densUnit)
  !         statusSystem = system ("lowdin-DFT.x SCF_DFT "//trim(densFile))
  !         call WaveFunction_readExchangeCorrelationMatrix( trim(nameOfSpecie))
  !      end if

       
  !      call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie) )
  !      call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )

  !      call SingleSCF_iterate( nameOfSpecie )

  !      diisError = SingleSCF_getDiisError(speciesID)

  !      if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then

  !         typeConvergence = " "
  !         if ( diisError > CONTROL_instance%DIIS_SWITCH_THRESHOLD ) typeConvergence = "*"

  !         if (abs(diisError) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
  !            write (6,"(I5,F20.10,F20.10,A20,A1)") SingleSCF_getNumberOfIterations( speciesID ), &
  !                 SingleSCF_getLastEnergy( speciesID ), &
  !                 SingleSCF_getStandardDeviationOfDensityMatrix( speciesID ) , "         --         ",typeConvergence
  !         else
  !            write (6,"(I5,F20.10,F20.10,F20.10,A1)") SingleSCF_getNumberOfIterations( speciesID), &
  !                 SingleSCF_getLastEnergy( speciesID ), &
  !                 SingleSCF_getStandardDeviationOfDensityMatrix( speciesID ), diisError,typeConvergence
  !         end if

  !      end if

  !      MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, tolerace )

  !   end do

  !   if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then

  !      call SingleSCF_showIteratiosStatus( MultiSCF_instance%status, nameOfSpecie )

  !      write(*,*) "... end SCF calculation"

  !   end if

  ! end subroutine MultiSCF_iterateUniqueSpecie

  !>
  !! @brief Ajusta el criterio de parada que se emplea en los ciclos SCF de
  !!		entre particulas de la misma especie.
  ! subroutine MultiSCF_setStopingThreshold( electronicThreshold, nonElectronicThreshold )
  !   implicit none
  !   real(8), intent(in) :: electronicThreshold
  !   real(8), intent(in) :: nonElectronicThreshold


  !   MultiSCF_instance%electronicTolerance = electronicThreshold
  !   MultiSCF_instance%nonelectronicTolerance = nonElectronicThreshold


  ! end subroutine MultiSCF_setStopingThreshold
