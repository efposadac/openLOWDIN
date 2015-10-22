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
  integer, parameter :: SCF_GLOBAL_CONVERGENCE_FAILED	= 0
  integer, parameter :: SCF_GLOBAL_CONVERGENCE_CONTINUE = 1
  integer, parameter :: SCF_GLOBAL_CONVERGENCE_SUCCESS = 2
  !< }

  !< enum Iteration_scheme_type {
  integer, parameter :: SCHEME_NONELECRONIC_FULLY_CONVERGED_BY_ELECTRONIC_ITERATION = 0
  integer, parameter :: SCHEME_ELECTRONIC_FULLY_CONVERGED_BY_NONELECRONIC_ITERATION = 1
  integer, parameter :: SCHEME_SPECIE_FULLY_CONVERGED_INDIVIDIALLY = 2
  integer, parameter :: SCHEME_SIMULTANEOUS = 3
  !< }

  !< For effective Fock matrix for restricted open-shell SCF
  logical :: isROHF
  !>

  type, public :: MultiSCF

     type(List) :: energyOMNE
     character(30) :: name
     integer :: numberOfIterations
     integer :: status
     integer :: nproc
     real(8) :: electronicTolerance
     real(8) :: nonelectronicTolerance

     !! Global energies
     real(8) :: totalEnergy
     real(8) :: totalCouplingEnergy
     real(8) :: electronicRepulsionEnergy

		 !! Cosmo
     real(8) :: cosmo3Energy


  end type MultiSCF

  type(MultiSCF), public, target :: MultiSCF_instance

  private :: &
       MultiSCF_iterateNonElectronicFullyByElectronicIteration, &
       MultiSCF_iterateElectronicFullyByNonElectronicIteration, &
       MultiSCF_iterateUniqueSpecie, &
       MultiSCF_iterateSpecieFullyConvergedIndividually

  public :: &
       MultiSCF_constructor, &
       MultiSCF_destructor, &
       MultiSCF_setStopingThreshold, &
       MultiSCF_getNumberOfIterations, &
       MultiSCF_getLastEnergy, &
       MultiSCF_iterate, &
       MultiSCF_restart, &
       MultiSCF_reset, &
       MultiSCF_testEnergyChange

contains


  !>
  !! @brief Define el constructor para la clase
  subroutine MultiSCF_constructor(nproc)
    implicit none

    integer :: nproc

    isROHF = .false.

    call List_constructor( MultiSCF_instance%energyOMNE,"ENERGY", CONTROL_instance%LISTS_SIZE)
    MultiSCF_instance%numberOfIterations = 0
    MultiSCF_instance%status = 0
    MultiSCF_instance%nproc = nproc

    if ( CONTROL_instance%OPTIMIZE .and.  .not.CONTROL_instance%MINIMIZATION_WITH_SINGLE_POINT ) then

       MultiSCF_instance%electronicTolerance = CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE
       MultiSCF_instance%nonelectronicTolerance = CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE

    else

       if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "DENSITY" ) then

          MultiSCF_instance%electronicTolerance = CONTROL_instance%ELECTRONIC_DENSITY_MATRIX_TOLERANCE
          MultiSCF_instance%nonelectronicTolerance = CONTROL_instance%NONELECTRONIC_DENSITY_MATRIX_TOLERANCE

       else

          MultiSCF_instance%electronicTolerance = CONTROL_instance%SCF_ELECTRONIC_ENERGY_TOLERANCE
          MultiSCF_instance%nonelectronicTolerance = CONTROL_instance%SCF_NONELECTRONIC_ENERGY_TOLERANCE

       end if

    end if

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
  !! @brief Ajusta el criterio de parada que se emplea en los ciclos SCF de
  !!		entre particulas de la misma especie.
  subroutine MultiSCF_setStopingThreshold( electronicThreshold, nonElectronicThreshold )
    implicit none
    real(8), intent(in) :: electronicThreshold
    real(8), intent(in) :: nonElectronicThreshold


    MultiSCF_instance%electronicTolerance = electronicThreshold
    MultiSCF_instance%nonelectronicTolerance = nonElectronicThreshold


  end subroutine MultiSCF_setStopingThreshold

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
  !! @param nameOfSpecie nombre de la especie seleccionada.
  !! @todo Adicionar condiciones de existencia de electrones en caso 1 y 2
  !!		y plantear caso analogo iterando sobre la particula mas ligera
  subroutine MultiSCF_iterate( iterationScheme, ROHF)
    implicit none
    integer, intent(in) :: iterationScheme
    logical, optional :: ROHF

    if(present(ROHF)) isROHF = ROHF

    select  case( iterationScheme )

    case( SCHEME_NONELECRONIC_FULLY_CONVERGED_BY_ELECTRONIC_ITERATION )

       !! done single - multi
       call MultiSCF_iterateNonElectronicFullyByElectronicIteration()

    case( SCHEME_ELECTRONIC_FULLY_CONVERGED_BY_NONELECRONIC_ITERATION )

       !! done single - multi
       call MultiSCF_iterateElectronicFullyByNonElectronicIteration()

    case( SCHEME_SPECIE_FULLY_CONVERGED_INDIVIDIALLY )

       !! done single - multi
       call MultiSCF_iterateSpecieFullyConvergedIndividually()

    case( SCHEME_SIMULTANEOUS )

       !! done single - multi
       call MultiSCF_iterateSimultaneous()

    case default
       call MultiSCF_iterateElectronicFullyByNonElectronicIteration()

    end select

  end subroutine MultiSCF_iterate

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
  !! @brief Realiza la convergencia completa de la parte no electronica por cada iteracion
  !!		SCF de la parte electronica
  !! @todo falta implementacion completa
  subroutine MultiSCF_iterateNonElectronicFullyByElectronicIteration()
    implicit none

    integer :: iteratorOfSpecie
    integer :: iteratorOfElectronicSpecie
    character(30) :: nameOfSpecie
    character(30) :: nameOfElectronicSpecie
    character(30) :: nameOfInitialSpecie
    real(8) :: tolerace

    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

    iteratorOfSpecie = 1
    nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)    

    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

       species_loop : do

          if(iteratorOfSpecie <= MolecularSystem_getNumberOfQuantumSpecies()) then          

             nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)   

             if( .not. MolecularSystem_instance%species(iteratorOfSpecie)%isElectron) then

                iteratorOfSpecie = iteratorOfSpecie + 1 
                cycle species_loop

             else 

                non_electronic_loop : do iteratorOfElectronicSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

                   if( .not. MolecularSystem_instance%species(iteratorOfElectronicSpecie)%isElectron ) then

                      MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

                      nameOfElectronicSpecie = MolecularSystem_getNameOfSpecie(iteratorOfElectronicSpecie)                   
                      tolerace = MultiSCF_instance%electronicTolerance

                      do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
                           (SingleSCF_getNumberOfIterations(iteratorOfElectronicSpecie) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

                         call WaveFunction_buildTwoParticlesMatrix( trim(nameOfElectronicSpecie), nproc = MultiSCF_instance%nproc )

                         if (CONTROL_instance%COSMO) then
                            call  WaveFunction_buildCosmo2Matrix( trim(nameOfElectronicSpecie))
                            if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
                               call WaveFunction_buildCosmoCoupling( trim(nameOfElectronicSpecie) )
                            end if
                         end if

                         !! At first iteration is not included the coupling operator.
                         if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
                            call WaveFunction_buildCouplingMatrix( trim(nameOfElectronicSpecie) )
                         end if

                         !! Build Fock Matrix
                         call WaveFunction_buildFockMatrix( trim(nameOfElectronicSpecie) )
                         waveFunction_instance(iteratorOfElectronicSpecie)%wasBuiltFockMatrix = .true.

                         !! Perform SCF iteration for Single species
                         call SingleSCF_iterate( trim(nameOfElectronicSpecie), nproc = MultiSCF_instance%nproc )

                         !! Test energy or density change (as requested)
                         if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

                            MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, tolerace )

                            MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, &
                                 MultiSCF_instance%electronicTolerance )

                         else

                            MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfElectronicSpecie, &
                                 MultiSCF_instance%electronicTolerance )

                         end if

                      end do !! end convergence for electrons

                   end if

                end do non_electronic_loop

                !! Realiza iteracion SCF para una especie cuantica particular
                call SingleSCF_iterate( trim( nameOfSpecie ), nproc = MultiSCF_instance%nproc )

                if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

                   MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
                        MultiSCF_instance%nonelectronicTolerance )
                else

                   MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
                        MultiSCF_instance%nonelectronicTolerance )

                end if

                if ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
                     (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS ) ) then

                   cycle species_loop

                else if ( MultiSCF_instance%status == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then

                   iteratorOfSpecie = iteratorOfSpecie + 1 
                   cycle species_loop                      

                end if

             end if

          else

             !! Finaliza un ciclo SCF completo para para todas la especies presentes
             call WaveFunction_obtainTotalEnergy(&
                  MultiSCF_instance%totalEnergy, &
                  MultiSCF_instance%totalCouplingEnergy, &
                  MultiSCF_instance%electronicRepulsionEnergy, &
                  MultiSCF_instance%cosmo3Energy )

             call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
             MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1              
             exit species_loop

          end if

       end do species_loop

    else

       if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
          call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
       end if

       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
            MultiSCF_instance%cosmo3Energy )

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

    end if

  end subroutine MultiSCF_iterateNonElectronicFullyByElectronicIteration

  !>
  !! @brief Realiza la convergencia completa de la parte electronica por cada iteracion
  !!	    SCF de la parte no electronica
  subroutine MultiSCF_iterateElectronicFullyByNonElectronicIteration( )
    implicit none
    integer :: iteratorOfSpecie
    integer :: iteratorOfElectronicSpecie
    character(30) :: nameOfSpecie
    character(30) :: nameOfElectronicSpecie
    character(30) :: nameOfInitialSpecie
    real(8) :: tolerace


    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

    iteratorOfSpecie = 1
    nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)    

    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

       species_loop : do

          if(iteratorOfSpecie <= MolecularSystem_getNumberOfQuantumSpecies()) then          

             nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)   

             if(MolecularSystem_instance%species(iteratorOfSpecie)%isElectron) then

                iteratorOfSpecie = iteratorOfSpecie + 1 
                cycle species_loop

             else 

                electronic_loop : do iteratorOfElectronicSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

                   if(MolecularSystem_instance%species(iteratorOfElectronicSpecie)%isElectron) then

                      MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

                      nameOfElectronicSpecie = MolecularSystem_getNameOfSpecie(iteratorOfElectronicSpecie)                   
                      tolerace = MultiSCF_instance%electronicTolerance

                      do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
                           (SingleSCF_getNumberOfIterations(iteratorOfElectronicSpecie) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

                         call WaveFunction_buildTwoParticlesMatrix( trim(nameOfElectronicSpecie), nproc = MultiSCF_instance%nproc )

                         !! At first iteration is not included the coupling operator.
                         if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
                            call WaveFunction_buildCouplingMatrix( trim(nameOfElectronicSpecie) )
                         end if

                         ! call Matrix_show(wavefunction_instance(iteratorOfSpecie)%couplingMatrix)


                         if (CONTROL_instance%COSMO) then
                            call  WaveFunction_buildCosmo2Matrix( trim(nameOfElectronicSpecie))
                            if(SingleSCF_getNumberOfIterations( iteratorOfElectronicSpecie ) > 0) then
                               call WaveFunction_buildCosmoCoupling( trim(nameOfElectronicSpecie) )
                            end if
                         end if



                         !! Build Fock Matrix
                         call WaveFunction_buildFockMatrix( trim(nameOfElectronicSpecie) )
                         waveFunction_instance(iteratorOfElectronicSpecie)%wasBuiltFockMatrix = .true.

                         !! Perform SCF iteration for Single species
                         call SingleSCF_iterate( trim(nameOfElectronicSpecie), nproc = MultiSCF_instance%nproc )

                         !! Test energy or density change (as requested)
                         if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

                            MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, tolerace )

                            MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfElectronicSpecie, &
                                 MultiSCF_instance%electronicTolerance )

                         else

                            MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfElectronicSpecie, &
                                 MultiSCF_instance%electronicTolerance )

                         end if

                      end do !! end convergence for electrons

                   end if

                end do electronic_loop

                !! Realiza iteracion SCF para una especie cuantica particular
                call SingleSCF_iterate( trim( nameOfSpecie ), nproc = MultiSCF_instance%nproc )

                if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "density" ) then

                   MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
                        MultiSCF_instance%nonelectronicTolerance )
                else

                   MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
                        MultiSCF_instance%nonelectronicTolerance )

                end if

                if ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
                     (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_NONELECTRONIC_MAX_ITERATIONS ) ) then

                   cycle species_loop

                else if ( MultiSCF_instance%status == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then

                   iteratorOfSpecie = iteratorOfSpecie + 1 
                   cycle species_loop                      

                end if

             end if

          else

             !! Finaliza un ciclo SCF completo para para todas la especies presentes
             call WaveFunction_obtainTotalEnergy(&
                  MultiSCF_instance%totalEnergy, &
                  MultiSCF_instance%totalCouplingEnergy, &
                  MultiSCF_instance%electronicRepulsionEnergy, &
									MultiSCF_instance%cosmo3Energy)

             call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
             MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1              
             exit species_loop

          end if

       end do species_loop

    else

       if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
          call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
       end if

       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
						MultiSCF_instance%cosmo3Energy)

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

    end if

  end subroutine MultiSCF_iterateElectronicFullyByNonElectronicIteration

  !>
  !! @brief Realiza la convergencia completa de cada especie de manera individual
  subroutine MultiSCF_iterateSpecieFullyConvergedIndividually()
    implicit none
    integer :: iteratorOfSpecie
    character(30) :: nameOfSpecie
    character(30) :: nameOfInitialSpecie
    real(8) :: tolerace

    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

    write(*,*)"entre a SpecieFullyConvergedIndividually"

    iteratorOfSpecie = 1
    nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)

    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

       do iteratorOfSpecie = 1, MolecularSystem_getNumberOfQuantumSpecies()

          nameOfSpecie = MolecularSystem_getNameOfSpecie(iteratorOfSpecie)

          if ( MolecularSystem_instance%species(iteratorOfSpecie)%isElectron ) then
             tolerace = MultiSCF_instance%electronicTolerance
          else
             tolerace = MultiSCF_instance%nonelectronicTolerance
          end if

          do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
               (SingleSCF_getNumberOfIterations( iteratorOfSpecie ) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

             call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

             if (CONTROL_instance%COSMO) then
                call  WaveFunction_buildCosmo2Matrix( trim(nameOfSpecie))
	                if(SingleSCF_getNumberOfIterations( iteratorOfSpecie ) > 0) then
                call WaveFunction_buildCosmoCoupling( trim(nameOfSpecie) )
                end if
             end if

             !! At first iteration is not included the coupling operator.
             if(SingleSCF_getNumberOfIterations( iteratorOfSpecie ) > 0) then
                call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie) )
             end if

             !! Build Fock Matrix
             call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )
             waveFunction_instance( iteratorOfSpecie )%wasBuiltFockMatrix = .true.

             !! Perform SCF iteration for Single species

             call SingleSCF_iterate( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

             !! Test energy or density change (as requested)
             if ( trim(CONTROL_instance%SCF_CONVERGENCE_CRITERIUM) == "DENSITY" ) then

                MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, tolerace )

                MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, &
                     MultiSCF_instance%electronicTolerance )

             else

                MultiSCF_instance%status =  SingleSCF_testEnergyChange( nameOfSpecie, &
                     MultiSCF_instance%electronicTolerance )

             end if

          end do

          MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE

       end do

       !! Finaliza un ciclo SCF completo para para todas la especies presentes
       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
						MultiSCF_instance%cosmo3Energy)

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)
       MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1       

    else

       if (.not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
          call MultiSCF_iterateUniqueSpecie(iteratorOfSpecie )
       end if

       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
						MultiSCF_instance%cosmo3Energy)

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)

    end if

  end subroutine MultiSCF_iterateSpecieFullyConvergedIndividually


  subroutine MultiSCF_iterateSimultaneous()
    implicit none

    integer :: i
    integer :: numberOfSpecies
    character(30) :: nameOfSpecie
    character(30) :: nameOfInitialSpecie
    real(8) :: startTime, endTime
    real(8) :: time1,time2
    logical :: auxValue

    MultiSCF_instance%status =  SCF_INTRASPECIES_CONVERGENCE_CONTINUE
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if( numberOfSpecies > 1 ) then

       do i = 1, numberOfSpecies

          nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
          call SingleSCF_iterate( trim(nameOfSpecie), actualizeDensityMatrix=.false., nproc = MultiSCF_instance%nproc )
       end do

       !! Get effective Fock matrix for restricted open-shell SCF
       ! if( isROHF ) then
       ! call MultiSCF_instance%mergeFockMatrix()
       ! end if

       do i=1, numberOfSpecies

          nameOfSpecie = MolecularSystem_getNameOfSpecie(i)
          call SingleSCF_actualizeDensityMatrix( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

       end do

       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
						MultiSCF_instance%cosmo3Energy)

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)             
       MultiSCF_instance%numberOfIterations = MultiSCF_instance%numberOfIterations + 1
       return

    else

       !!Especifica el procedimiento a seguir si lo electrones han sido congelados
       if ( .not. CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS) then
          call MultiSCF_iterateUniqueSpecie( speciesID = numberOfSpecies )
       end if

       call WaveFunction_obtainTotalEnergy(&
            MultiSCF_instance%totalEnergy, &
            MultiSCF_instance%totalCouplingEnergy, &
            MultiSCF_instance%electronicRepulsionEnergy, &
						MultiSCF_instance%cosmo3Energy)

       call List_push_back( MultiSCF_instance%energyOMNE, MultiSCF_instance%totalEnergy)             

    end if

  end subroutine MultiSCF_iterateSimultaneous

  !>
  !! @brief Realiza la convergencia completa de la unica especie presente
  subroutine MultiSCF_iterateUniqueSpecie( speciesID )
    implicit none

    integer :: speciesID

    character(30) :: nameOfSpecie    
    real(8) :: tolerace
    real(8) :: diisError
    character :: typeConvergence

    nameOfSpecie = MolecularSystem_getNameOfSpecie(speciesID)

    if ( MolecularSystem_instance%species(speciesID)%isElectron ) then
       tolerace = MultiSCF_instance%electronicTolerance
    else
       tolerace = MultiSCF_instance%nonelectronicTolerance
    end if

    !! Realiza ciclo de convergencia SCF para la unica especie cuantica presente

    if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then
       print *,""
       print *,"Begin SCF calculation by: ",trim(nameOfSpecie)
       print *,"-------------------------"
       print *,""
       print *,"-----------------------------------------------------------------"
       write (6,"(A10,A12,A25,A20)") "Iteration", "Energy", " Density Change","         DIIS Error "
       print *,"-----------------------------------------------------------------"
    end if

    ! write(*,*)"entre al unique specie"
    do while ( ( MultiSCF_instance%status ==  SCF_INTRASPECIES_CONVERGENCE_CONTINUE ) .and. &
         ( SingleSCF_getNumberOfIterations(speciesID) <= CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) )

       call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), nproc = MultiSCF_instance%nproc )

       if (CONTROL_instance%COSMO) then
          call  WaveFunction_buildCosmo2Matrix( trim(nameOfSpecie))
          call WaveFunction_buildCosmoCoupling( trim(nameOfSpecie) )
       end if

       call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie) )
       call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )

       waveFunction_instance(speciesID)%wasBuiltFockMatrix = .true.

       call SingleSCF_iterate( nameOfSpecie, nproc = MultiSCF_instance%nproc )

       diisError = SingleSCF_getDiisError(speciesID)

       if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then

          typeConvergence = " "
          if ( diisError > CONTROL_instance%DIIS_SWITCH_THRESHOLD ) typeConvergence = "*"

          if (abs(diisError) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
             write (6,"(I5,F20.10,F20.10,A20,A1)") SingleSCF_getNumberOfIterations( speciesID ), &
                  SingleSCF_getLastEnergy( speciesID ), &
                  SingleSCF_getStandardDeviationOfDensityMatrix( speciesID ) , "         --         ",typeConvergence
          else
             write (6,"(I5,F20.10,F20.10,F20.10,A1)") SingleSCF_getNumberOfIterations( speciesID), &
                  SingleSCF_getLastEnergy( speciesID ), &
                  SingleSCF_getStandardDeviationOfDensityMatrix( speciesID ), diisError,typeConvergence
          end if

       end if

       MultiSCF_instance%status = SingleSCF_testDensityMatrixChange( nameOfSpecie, tolerace )

    end do

    if ( .not.CONTROL_instance%OPTIMIZE .or. CONTROL_instance%DEBUG_SCFS ) then

       call SingleSCF_showIteratiosStatus( MultiSCF_instance%status, nameOfSpecie )

       write(*,*) "... end SCF calculation"

    end if

  end subroutine MultiSCF_iterateUniqueSpecie

  !>
  !! @brief Prueba si la energia del la ultima iteracion a sufrido un cambio por debajo de cierta tolerancia
  function MultiSCF_testEnergyChange(  tolerace ) result( output )
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

    if ( MultiSCF_instance%numberOfIterations > CONTROL_instance%SCF_GLOBAL_MAXIMUM_ITERATIONS ) then

       output =	SCF_GLOBAL_CONVERGENCE_SUCCESS
       call List_end( MultiSCF_instance%energyOMNE )
       finalEnergy= List_current( MultiSCF_instance%energyOMNE )

       !! Obtiene el valor  de energia anterior a la ultima iteracion
       call List_iterate( MultiSCF_instance%energyOMNE, -1 )

       !! Obtiene el cambio de energia en las ultimas dos iteraciones
       deltaEnergy = finalEnergy - List_current( MultiSCF_instance%energyOMNE )
       print *,""
       write (6,"(A20,A20)") "   Current Energy   ", "   Current Change   "
       print *,"-----------------------------------------------"
       write (6,"(F20.10,F20.10)") finalEnergy, deltaEnergy

       print *,  "The number of Iterations was exceded, the convergence had failed"

    else

       if ( MultiSCF_instance%numberOfIterations > 1 ) then

          auxVar=.true.
          do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()

             nameOfSpecie = MolecularSystem_getNameOfSpecie(speciesID)

             toleraceOfSpecie = MultiSCF_instance%electronicTolerance

             if (.not. MolecularSystem_instance%species(speciesID)%isElectron ) then
                toleraceOfSpecie = MultiSCF_instance%nonelectronicTolerance
             end if

             if ( SingleSCF_testDensityMatrixChange( nameOfSpecie, &
                  toleraceOfSpecie ) == SCF_INTRASPECIES_CONVERGENCE_SUCCESS ) then
                auxVar = auxVar .and. .true.
             else
                auxVar = auxVar .and. .false.
             end if

          end do

          call List_end( MultiSCF_instance%energyOMNE )
          deltaEnergy= List_current( MultiSCF_instance%energyOMNE )

          !! Obtiene el valor  de energia anterior a la ultima iteracion
          call List_iterate( MultiSCF_instance%energyOMNE, -1 )

          !! Obtiene el cambio de energia en las ultimas dos iteraciones
          deltaEnergy = deltaEnergy - List_current( MultiSCF_instance%energyOMNE )

          if( ( ( abs( deltaEnergy ) < tolerace ) .and. auxVar ) .or. abs(deltaEnergy) < CONTROL_instance%STRONG_ENERGY_TOLERANCE ) then
             output =	SCF_GLOBAL_CONVERGENCE_SUCCESS
          else
             output =	SCF_GLOBAL_CONVERGENCE_CONTINUE
          end if

       else

          output = SCF_GLOBAL_CONVERGENCE_CONTINUE

       end if

    end if

  end function MultiSCF_testEnergyChange

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
