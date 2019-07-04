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
!! @brief Modulo para diagonalizacion de matrices de especies cuanticas
!!  Este modulo define una seudoclase para realizacion de ciclos SCF de
!! especies cuanticas. (una sola especie)
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
!!
module SingleSCF_
  use Exception_
  use CONTROL_
  use Matrix_
  use Vector_
  use List_
  use Convergence_
  use WaveFunction_
  use MolecularSystem_
  implicit none

  !< enum Matrix_type {
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_FAILED =  0
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_CONTINUE =  1
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_SUCCESS =  2
  !< }

  public :: &
       SingleSCF_showIteratiosStatus, &
       SingleSCF_setNumberOfIterations,&
       SingleSCF_getNumberOfIterations,&
       SingleSCF_getLastEnergy, &
       SingleSCF_getStandardDeviationOfDensityMatrix,&
       SingleSCF_getDiisError, &
       SingleSCF_iterate, &
       SingleSCF_restart, &
       SingleSCF_reset, &
       SingleSCF_actualizeDensityMatrix, &
       SingleSCF_testEnergyChange, &
       SingleSCF_testDensityMatrixChange

contains

  !>
  !! @brief Retorna el numero de iteraciones realizadas para una especie especificada
  subroutine SingleSCF_showIteratiosStatus( status, nameOfSpecie )
    implicit none
    character(*), intent(in) :: nameOfSpecie
    integer :: status

    select case(status)

    case(SCF_INTRASPECIES_CONVERGENCE_FAILED)

       print *,""
       print *,"SCF STATUS: CONVERGENCE FAILED BY: ",trim(nameOfSpecie)
       print *,""

    case(SCF_INTRASPECIES_CONVERGENCE_CONTINUE)

       print *,""
       print *,"SCF STATUS: CONVERGENCE CONTINUE BY: ",trim(nameOfSpecie)
       print *,""

    case(SCF_INTRASPECIES_CONVERGENCE_SUCCESS)

       print *,""
       print *,"SCF STATUS: CONVERGENCE SUCCESS BY: ",trim(nameOfSpecie)
       print *,""

    end select

  end subroutine SingleSCF_showIteratiosStatus

  !>
  !! @brief set the iterations number
  subroutine SingleSCF_setNumberOfIterations( numberOfIterations )
    implicit none
    integer, intent(in) :: numberOfIterations

    WaveFunction_instance%numberOfIterations = numberOfIterations

  end subroutine SingleSCF_setNumberOfIterations

  !>
  !! @brief Retorna el numero de iteraciones realizadas para una especie especificada
  function SingleSCF_getNumberOfIterations( speciesID ) result( output )
    implicit none

    integer :: speciesID
    integer :: output

    output = WaveFunction_instance(speciesID)%numberOfIterations

  end function SingleSCF_getNumberOfIterations

  !>
  !! @brief Retorna la energia calculada en la ultima iteracion
  function SingleSCF_getLastEnergy( speciesID ) result( output )
    implicit none

    integer :: speciesID
    real(8) :: output

    !! Obtiene el valor de energia de la ultima iteracion
    call List_end( WaveFunction_instance(speciesID)%energySCF )
    output= List_current( WaveFunction_instance(speciesID)%energySCF )

  end function SingleSCF_getLastEnergy

  !>
  !! @brief Retorna de la densviacion estandar de la matrix de densidad en la ultima iteracion
  function SingleSCF_getStandardDeviationOfDensityMatrix( speciesID ) result( output )
    implicit none

    integer :: speciesID
    real(8) :: output

    call List_end( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )
    output= List_current( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )

  end function SingleSCF_getStandardDeviationOfDensityMatrix

  !>
  !! @brief Retorna de la densviacion estandar de la matrix de densidad en la ultima iteracion
  function SingleSCF_getDiisError(speciesID ) result( output )
    implicit none
    real(8) :: output

    integer :: speciesID

    call List_end( WaveFunction_instance(speciesID)%diisError )
    output= List_current( WaveFunction_instance(speciesID)%diisError )

  end function SingleSCF_getDiisError

  !>
  !! @brief Realiza iteracion SCF
  !! @param nameOfSpecie nombre de la especie seleccionada.
  !! @warning Se debe garantizar la actualizacion de la matriz de densidad  y calculos de esta derivados
  !! cuando se emplee este metodo con el parametro actualizeDensityMatrix=.false.
  subroutine SingleSCF_iterate( nameOfSpecie, actualizeDensityMatrix )
    implicit none

    character(*), optional :: nameOfSpecie
    logical,optional :: actualizeDensityMatrix

    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: previousWavefunctionCoefficients, matchingMatrix, matchingMatrix2
    type(Vector) :: deltaVector 
    type(Matrix) :: auxOverlapMatrix

    character(30) :: nameOfSpecieSelected
    integer(8) :: total3
    integer(8) :: numberOfContractions
    real(8) :: threshold, hold
    real(8) :: levelShiftingFactor, deltaEnergy    
    integer :: speciesID
    integer :: ocupados, occupatedOrbitals, prodigals
    integer :: totales
    integer :: search, trial
    integer :: ii, jj, astrayOrbitals, startIteration
    integer :: i
    logical :: existFile
    logical :: internalActualizeDensityMatrix 

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

!    wfnFile = trim(CONTROL_instance%INPUT_FILE)//"lowdin.vec"
    wfnUnit = 30

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    internalActualizeDensityMatrix = .true.

    if ( present(  actualizeDensityMatrix ) ) then
       internalActualizeDensityMatrix =  actualizeDensityMatrix
    end if

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

    numberOfContractions = MolecularSystem_getTotalnumberOfContractions( speciesID )

    !!**********************************************************************************************
    !! If fock matrix has not been created! creates the fock matrix
    !!
    if ( .not. WaveFunction_instance(speciesID)%wasBuiltFockMatrix ) then

       ! call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecie), nproc )
       ! call WaveFunction_buildCouplingMatrix( trim(nameOfSpecie) )       


       if (CONTROL_instance%COSMO) then
          call WaveFunction_buildCosmo2Matrix( trim(nameOfSpecie) )
          call WaveFunction_buildCosmoCoupling( trim(nameOfSpecie) )
       end if


       call WaveFunction_buildFockMatrix( trim(nameOfSpecie) )

       WaveFunction_instance(speciesID)%wasBuiltFockMatrix = .true.


    end if
    !!
    !!**********************************************************************************************

    if ( .not.CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS ) then

       call SingleSCF_convergenceMethod( speciesID )

       call Matrix_copyConstructor( fockMatrixTransformed, WaveFunction_instance( speciesID )%fockMatrix )

       
       !!**********************************************************************************************
       !! Level Shifting Convergence Method       
       !!
       if ( CONTROL_instance%ACTIVATE_LEVEL_SHIFTING .eqv. .true. ) then

          if ( MolecularSystem_instance%species(speciesID)%isElectron ) then
             levelShiftingFactor=CONTROL_instance%ELECTRONIC_LEVEL_SHIFTING
          else
             levelShiftingFactor=CONTROL_instance%NONELECTRONIC_LEVEL_SHIFTING
          end if

          !The level shifting should be turned off when we are close to convergence, should be a parameter in CONTROL 
          if ( WaveFunction_instance(speciesID)%numberOfIterations .gt. 1 ) then
             !! Obtiene el valor de energia de la ultima iteracion
             call List_end( WaveFunction_instance(speciesID)%energySCF )
             deltaEnergy= List_current( WaveFunction_instance(speciesID)%energySCF )
             
             !! Obtiene el valor  de energia anterior a la ultima iteracion
             call List_iterate( WaveFunction_instance(speciesID)%energySCF, -1)
             
             !! Obtiene el cambio de energia en las ultimas dos iteraciones
             deltaEnergy = abs(deltaEnergy - List_current( WaveFunction_instance(speciesID)%energySCF ))
             
             ! if( levelShiftingFactor .gt. 0.0_8 .and. deltaEnergy .lt. 1E-6 ) then
             !    levelShiftingFactor = 0.0_8
             !    print *, "turned off level shifting for ",  trim(nameOfSpecie) 
             ! end if
          end if
          
          fockMatrixTransformed%values = &
               matmul( matmul( transpose( WaveFunction_instance( speciesID )%waveFunctionCoefficients%values ) , &
               fockMatrixTransformed%values), WaveFunction_instance( speciesID )%waveFunctionCoefficients%values )

          totales = Matrix_getNumberOfRows(fockMatrixTransformed)
          ocupados = MolecularSystem_getOcupationNumber( speciesID )

          do i= ocupados + 1, totales
             fockMatrixTransformed%values(i,i) = levelShiftingFactor + fockMatrixTransformed%values(i,i)
          end do

          fockMatrixTransformed%values = &
               matmul( matmul(WaveFunction_instance( speciesID )%waveFunctionCoefficients%values, &
               fockMatrixTransformed%values), transpose(WaveFunction_instance( speciesID )%waveFunctionCoefficients%values ))

          call Matrix_copyConstructor(auxOverlapMatrix,WaveFunction_instance( speciesID )%OverlapMatrix )

          fockMatrixTransformed%values = &
               matmul( matmul(auxOverlapMatrix%values, &
               fockMatrixTransformed%values), transpose(auxOverlapMatrix%values ))                    

       end if

       if (CONTROL_instance%MO_FRACTION_OCCUPATION < 1.0_8 .or. CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF ) then

          call Matrix_copyConstructor( previousWavefunctionCoefficients, WaveFunction_instance( speciesID )%waveFunctionCoefficients )

       end if       
       !!
       !! End of Level Shifting Convergence Routine       
       !!**********************************************************************************************

       ! print *, "Coefficients before for ", speciesID
       ! call Matrix_show(WaveFunction_instance( speciesID )%waveFunctionCoefficients)
       
       !!**********************************************************************************************
       !! Iteration begins
       !!
       fockMatrixTransformed%values = &
            matmul( matmul( transpose( WaveFunction_instance(speciesID)%transformationMatrix%values ) , &
            fockMatrixTransformed%values), WaveFunction_instance(speciesID)%transformationMatrix%values )

       !! Calcula valores y vectores propios de matriz de Fock transformada.
       call Matrix_eigen( fockMatrixTransformed, WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, &
            WaveFunction_instance(speciesID)%waveFunctionCoefficients, SYMMETRIC )

       !! Calcula los  vectores propios para matriz de Fock       
       WaveFunction_instance(speciesID)%waveFunctionCoefficients%values = &
            matmul( WaveFunction_instance(speciesID)%transformationMatrix%values, &
            WaveFunction_instance(speciesID)%waveFunctionCoefficients%values )
       !!
       !! Interation ends
       !!**********************************************************************************************

       ! print *, "Coefficients after for ", speciesID
       ! call Matrix_show(WaveFunction_instance( speciesID )%waveFunctionCoefficients)

       
       !!**********************************************************************************************
       !! Orbital exchange for TOM calculation
       !!
       if (CONTROL_instance%MO_FRACTION_OCCUPATION < 1.0_8 .or. CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF ) then

          startIteration=5

          threshold=0.8_8

          call Matrix_copyConstructor(auxOverlapMatrix,WaveFunction_instance(speciesID)%overlapMatrix)

          occupatedOrbitals = MolecularSystem_getOcupationNumber( speciesID )
          total3 = int(occupatedOrbitals,8)

          call Matrix_constructor (matchingMatrix, total3, total3)

          matchingMatrix%values = &
               matmul( matmul( transpose( previousWavefunctionCoefficients%values ) , &
               auxOverlapMatrix%values), WaveFunction_instance(speciesID)%waveFunctionCoefficients%values )

          astrayOrbitals=0 !!! number of orbitals that must be reordered

          do ii= 1, occupatedOrbitals

             !DEBUG
             ! print *,"Antes ","Orbital: ",ii," ",abs(matchingMatrix%values(ii,ii))," energy ",WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(ii)

             if ( abs(matchingMatrix%values(ii,ii)) < threshold ) then
                astrayOrbitals=astrayOrbitals+1
             end if

          end do

          !DEBUG
          ! print *,"Number of astrayOrbitals",astrayOrbitals            
          ! print *,"Initial MM:"
          ! call Matrix_show (matchingMatrix)

          call Vector_constructor (deltaVector,astrayOrbitals)

          jj=0
          do ii= 1, occupatedOrbitals
             if ( abs(matchingMatrix%values(ii,ii)) < threshold ) then
                jj=jj+1
                deltaVector%values(jj)=ii
             end if
          end do

          prodigals = astrayOrbitals

          if(astrayOrbitals .ge. 1 .and. WaveFunction_instance(speciesID)%numberOfIterations > startIteration) then
             ! print *, "Switching orbitals..."

             do ii=1,astrayOrbitals

                call Matrix_copyConstructor( auxiliaryMatrix, WaveFunction_instance(speciesID)%waveFunctionCoefficients )            
                call Matrix_copyConstructor( matchingMatrix2, matchingMatrix )            

                search=deltaVector%values(ii)            
                ! print *,"search: ",search

                do jj=1, numberOfContractions

                   !trial=deltaVector%values(jj)
                   trial=jj
                   ! print *,"trial: ",trial

                   ! print *,"valor: ",abs(matchingMatrix%values(search,trial))

                   if ( abs(matchingMatrix%values(search,trial)) > threshold ) then

                      WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,search)=auxiliaryMatrix%values(:,trial)
                      WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,trial)=auxiliaryMatrix%values(:,search)

                      hold=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(search)
                      WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(search)=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(trial)
                      WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(trial)=hold

                      matchingMatrix%values(:,trial)=matchingMatrix2%values(:,search)
                      matchingMatrix%values(:,search)=matchingMatrix2%values(:,trial)

                      prodigals = prodigals - 1
                      
                      print *, "Switching orbital... ",search," with ",trial, " for ", trim(nameOfSpecie)

                      exit

                   end if

                end do

                ! if (prodigals<2) then
                !    exit
                ! end if

             end do

          end if

          ! print *, "Coefficients after switching for ", speciesID
          ! call Matrix_show(WaveFunction_instance( speciesID )%waveFunctionCoefficients)

          
          call Matrix_destructor(matchingMatrix)
          call Vector_destructor(deltaVector)           

       end if
       !!**********************************************************************************************
       !! Orbital exchange for TOM calculation
       !!
       if (CONTROL_instance%MO_FRACTION_OCCUPATION < 1.0_8 .or. CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF ) then

          startIteration=5

          threshold=0.8_8

          call Matrix_copyConstructor(auxOverlapMatrix,WaveFunction_instance(speciesID)%overlapMatrix)

          occupatedOrbitals = MolecularSystem_getOcupationNumber( speciesID )
          total3 = int(occupatedOrbitals,8)

          call Matrix_constructor (matchingMatrix, total3, total3)

          matchingMatrix%values = &
               matmul( matmul( transpose( previousWavefunctionCoefficients%values ) , &
               auxOverlapMatrix%values), WaveFunction_instance(speciesID)%waveFunctionCoefficients%values )

          astrayOrbitals=0 !!! number of orbitals that must be reordered

          do ii= 1, occupatedOrbitals

             !DEBUG
             ! print *,"Antes ","Orbital: ",ii," ",abs(matchingMatrix%values(ii,ii))," energy ",WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(ii)

             if ( abs(matchingMatrix%values(ii,ii)) < threshold ) then
                astrayOrbitals=astrayOrbitals+1
             end if

          end do

          !DEBUG
          ! print *,"Number of astrayOrbitals",astrayOrbitals            
          ! print *,"Initial MM:"
          ! call Matrix_show (matchingMatrix)

          call Vector_constructor (deltaVector,astrayOrbitals)

          jj=0
          do ii= 1, occupatedOrbitals
             if ( abs(matchingMatrix%values(ii,ii)) < threshold ) then
                jj=jj+1
                deltaVector%values(jj)=ii
             end if
          end do

          prodigals = astrayOrbitals

          if(astrayOrbitals .ge. 1 .and. WaveFunction_instance(speciesID)%numberOfIterations > startIteration) then
             ! print *, "Switching orbitals..."

             do ii=1,astrayOrbitals

                call Matrix_copyConstructor( auxiliaryMatrix, WaveFunction_instance(speciesID)%waveFunctionCoefficients )            
                call Matrix_copyConstructor( matchingMatrix2, matchingMatrix )            

                search=deltaVector%values(ii)            
                ! print *,"search: ",search

                do jj=1, numberOfContractions

                   !trial=deltaVector%values(jj)
                   trial=jj
                   ! print *,"trial: ",trial

                   ! print *,"valor: ",abs(matchingMatrix%values(search,trial))

                   if ( abs(matchingMatrix%values(search,trial)) > threshold ) then

                      WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,search)=auxiliaryMatrix%values(:,trial)
                      WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,trial)=auxiliaryMatrix%values(:,search)

                      hold=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(search)
                      WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(search)=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(trial)
                      WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(trial)=hold

                      matchingMatrix%values(:,trial)=matchingMatrix2%values(:,search)
                      matchingMatrix%values(:,search)=matchingMatrix2%values(:,trial)

                      prodigals = prodigals - 1
                      
                      print *, "Switching orbital... ",search," with ",trial, " for ", trim(nameOfSpecie)

                      exit

                   end if

                end do

                ! if (prodigals<2) then
                !    exit
                ! end if

             end do

          end if

          ! print *, "Coefficients after switching for ", speciesID
          ! call Matrix_show(WaveFunction_instance( speciesID )%waveFunctionCoefficients)

          
          call Matrix_destructor(matchingMatrix)
          call Vector_destructor(deltaVector)           

       end if
       !!
       !!**********************************************************************************************

       !! If NO SCF cicle is desired, read the coefficients from the ".vec" file again
       if ( (CONTROL_instance%NO_SCF .and. CONTROL_instance%READ_COEFFICIENTS) .or. &
            (SingleSCF_getNumberOfIterations(speciesID) == 0 .and. CONTROL_instance%READ_COEFFICIENTS) ) then

          arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
          arguments(1) = "COEFFICIENTS"

          wfnFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
          inquire(FILE = wfnFile, EXIST = existFile )

          if ( existFile) then
             open(unit=wfnUnit, file=trim(wfnFile), status="old", form="formatted")

             WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(unit=wfnUnit, &
                  rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.false.,  & 
                  arguments=arguments(1:2))

             close(wfnUnit)

          else 
             wfnFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
             inquire(FILE = wfnFile, EXIST = existFile )

             if ( existFile) then
                open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

                WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(unit=wfnUnit, &
                     rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., & 
                     arguments=arguments(1:2))

                close(wfnUnit)

             else
                call  SingleSCF_exception( ERROR, "I did not find any .vec coefficients file", "At SCF program, at SingleSCF_Iterate")
             end if

          end if
       end if

       !! If NO SCF cicle is desired, read the coefficients from the ".vec" file again
       if ( (CONTROL_instance%NO_SCF .and. CONTROL_instance%READ_EIGENVALUES) .or. &
            (SingleSCF_getNumberOfIterations(speciesID) == 0 .and. CONTROL_instance%READ_EIGENVALUES) ) then

          arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
          arguments(1) = "ORBITALS"

          wfnFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
          inquire(FILE = wfnFile, EXIST = existFile )

          if ( existFile) then
             open(unit=wfnUnit, file=trim(wfnFile), status="old", form="formatted")

                call Vector_getFromFile(unit=wfnUnit, &
                          output = WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, &
                          elementsNum= int(numberOfContractions,4), binary=.true., & 
                          arguments=arguments(1:2))

             close(wfnUnit)

          else 
             wfnFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
             inquire(FILE = wfnFile, EXIST = existFile )

             if ( existFile) then
                open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

                call Vector_getFromFile(unit=wfnUnit, &
                          output = WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, &
                          elementsNum= int(numberOfContractions,4), binary=.false., & 
                          arguments=arguments(1:2))

                close(wfnUnit)

             else
                call  SingleSCF_exception( ERROR, "I did not find any .vec coefficients file", "At SCF program, at SingleSCF_Iterate")
             end if

          end if

      end if

       !! Not implemented yet
       !! If NO SCF cicle is desired, read the coefficients from the ".vec" file again
       if ( CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS) then
          if ( CONTROL_instance%READ_COEFFICIENTS ) then
             if ( .not. MolecularSystem_instance%species(speciesID)%isElectron ) then

                inquire(FILE = trim(nameOfSpecieSelected)//".vec", EXIST = existFile )

                if ( existFile) then

                   WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(int(numberOfContractions,4), int(numberOfContractions,4), &
                        file=trim(nameOfSpecie)//".vec", binary = .false.)

                end if
             end if
          end if
       end if

       if ( CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS) then
          if ( CONTROL_instance%READ_COEFFICIENTS ) then
             if ( MolecularSystem_instance%species(speciesID)%isElectron ) then

                inquire(FILE = trim(nameOfSpecieSelected)//".vec", EXIST = existFile )

                if ( existFile) then

                   WaveFunction_instance(speciesID)%waveFunctionCoefficients = Matrix_getFromFile(int(numberOfContractions,4), int(numberOfContractions,4), &
                        file=trim(nameOfSpecie)//".vec", binary = .false.)

                end if
             end if
          end if
       end if





       !! Exchanging orbitals just for calculation excited states
       if( nameOfSpecie == trim(CONTROL_instance%EXCITE_SPECIE) ) then

          call Matrix_constructor (auxiliaryMatrix, numberOfContractions, numberOfContractions)
          call Matrix_copyConstructor( auxiliaryMatrix, WaveFunction_instance(speciesID)%waveFunctionCoefficients )
          WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,2)=auxiliaryMatrix%values(:,1)
          WaveFunction_instance(speciesID)%waveFunctionCoefficients%values(:,1)=auxiliaryMatrix%values(:,2)
          call Matrix_destructor(auxiliaryMatrix)
          hold=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(1)
          WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(1)=WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(2)
          WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values(2)=hold

       end if

       if (  internalActualizeDensityMatrix ) then

          !! Determina la desviacion estandar de los elementos de la matriz de densidad
          call Matrix_copyConstructor( WaveFunction_instance(speciesID)%beforeDensityMatrix, WaveFunction_instance(speciesID)%densityMatrix )
          call WaveFunction_builtDensityMatrix( trim(nameOfSpecieSelected) )
          call List_push_back( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements, &
               Matrix_standardDeviation( WaveFunction_instance(speciesID)%beforeDensityMatrix, WaveFunction_instance(speciesID)%densityMatrix ) )

          !! Calcula energia total para la especie especificada
          call WaveFunction_obtainTotalEnergyForSpecie( trim(nameOfSpecieSelected) )
          call List_push_back( WaveFunction_instance(speciesID)%energySCF, WaveFunction_instance(speciesID)%totalEnergyForSpecie )
          call List_push_back( WaveFunction_instance(speciesID)%diisError, Convergence_getDiisError( WaveFunction_instance(speciesID)%convergenceMethod) )

       end if
       
       !! Actualiza el contador para el numero de iteraciones SCF de la especie actual
       WaveFunction_instance(speciesID)%numberOfIterations =  WaveFunction_instance(speciesID)%numberOfIterations + 1
       call Matrix_destructor(fockMatrixTransformed)

    else

       !! Warning:
       !! This part has not been tested

       call WaveFunction_builtDensityMatrix(trim(nameOfSpecieSelected) )

       !! Calcula energia total para la especie especificada
       call WaveFunction_obtainTotalEnergyForSpecie( trim(nameOfSpecieSelected) )

       call List_push_back( WaveFunction_instance(speciesID)%energySCF, WaveFunction_instance(speciesID)%totalEnergyForSpecie )

    end if

    WaveFunction_instance(speciesID)%wasBuiltFockMatrix = .false.

  end subroutine SingleSCF_iterate
  

  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine SingleSCF_restart(speciesID)
    implicit none

    integer :: speciesID

    WaveFunction_instance(speciesID)%numberOfIterations = 0

    call List_clear( WaveFunction_instance(speciesID)%energySCF )
    call List_clear( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )
    call List_clear( WaveFunction_instance(speciesID)%diisError )
    call Convergence_reset()

  end subroutine SingleSCF_restart

  !>
  !! @brief Reinicializa el modulo
  subroutine SingleSCF_reset( speciesID )
    implicit none

    integer :: speciesID

    WaveFunction_instance(speciesID)%numberOfIterations = 0
    WaveFunction_instance(speciesID)%wasBuiltFockMatrix =.false.

    call List_clear( WaveFunction_instance(speciesID)%energySCF )
    call List_clear( WaveFunction_instance(speciesID)%diisError )
    call List_clear( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )

    call Convergence_destructor( WaveFunction_instance(speciesID)%convergenceMethod )
    call Convergence_constructor(WaveFunction_instance(speciesID)%convergenceMethod, &
         WaveFunction_instance(speciesID)%name,CONTROL_instance%CONVERGENCE_METHOD)

    call Convergence_reset()

  end subroutine SingleSCF_reset


  !>
  !! @brief Updates density matrix after one SCF iteration
  subroutine SingleSCF_actualizeDensityMatrix( nameOfSpecie )
    implicit none
    character(*), optional :: nameOfSpecie

    integer :: speciesID
    character(30) :: nameOfSpecieSelected

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

    !! Determina la desviacion estandar de los elementos de la matriz de densidad
    call Matrix_copyConstructor( WaveFunction_instance(speciesID)%beforeDensityMatrix, WaveFunction_instance(speciesID)%densityMatrix )


    call WaveFunction_builtDensityMatrix( nameOfSpecieSelected )

    call List_push_back( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements, &
         Matrix_standardDeviation( WaveFunction_instance(speciesID)%beforeDensityMatrix, WaveFunction_instance(speciesID)%densityMatrix) )

    !! Calcula energia total para la especie especificada
    call WaveFunction_obtainTotalEnergyForSpecie( nameOfSpecieSelected )

    call List_push_back( WaveFunction_instance(speciesID)%energySCF, WaveFunction_instance(speciesID)%totalEnergyForSpecie )

    call List_push_back( WaveFunction_instance(speciesID)%diisError, Convergence_getDiisError( WaveFunction_instance(speciesID)%convergenceMethod) )

  end subroutine SingleSCF_actualizeDensityMatrix

  !>
  !! @brief Prueba si la energia del la ultima iteracion a sufrido un cambio por debajo de cierta tolerancia
  function SingleSCF_testEnergyChange( nameOfSpecie, tolerace ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    real(8), intent(in) :: tolerace
    integer :: output

    integer :: speciesID
    character(30) :: nameOfSpecieSelected
    real(8) :: deltaEnergy

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID(trim(nameOfSpecieSelected))

    if ( SingleSCF_getNumberOfIterations( speciesID ) &
         >CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) then

       output = SCF_INTRASPECIES_CONVERGENCE_FAILED

       call SingleSCF_exception(ERROR, "SCF_INTRASPECIES_CONVERGENCE_FAILED BY: "//trim(nameOfSpecieSelected)//" ITERATIONS EXCEEDED", &
            "Class object SCF in the testEnergyChange function" )

    else

       if ( SingleSCF_getNumberOfIterations( speciesID ) > 1 ) then

          !! Obtiene el valor de energia de la ultima iteracion
          call List_end( WaveFunction_instance(speciesID)%energySCF )
          deltaEnergy= List_current( WaveFunction_instance(speciesID)%energySCF )

          !! Obtiene el valor  de energia anterior a la ultima iteracion
          call List_iterate( WaveFunction_instance(speciesID)%energySCF, -1)

          !! Obtiene el cambio de energia en las ultimas dos iteraciones
          deltaEnergy = deltaEnergy - List_current( WaveFunction_instance(speciesID)%energySCF )

          if( abs( deltaEnergy -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) < tolerace  ) then

             output = SCF_INTRASPECIES_CONVERGENCE_SUCCESS

          else

             output = SCF_INTRASPECIES_CONVERGENCE_CONTINUE


          end if

       else

          output = SCF_INTRASPECIES_CONVERGENCE_CONTINUE

       end if

    end if

  end function SingleSCF_testEnergyChange

  !>
  !! @brief Prueba si la matriz de densidad ha sufrido un cambio respecto a una tolerancia dada
  !!
  function SingleSCF_testDensityMatrixChange( nameOfSpecie, tolerace ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    real(8), intent(in) :: tolerace
    integer :: output

    integer :: speciesID
    character(30) :: nameOfSpecieSelected
    real(8) :: deltaDensityMatrix

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID(trim(nameOfSpecieSelected))

    if ( SingleSCF_getNumberOfIterations( speciesID ) &
         >CONTROL_instance%SCF_ELECTRONIC_MAX_ITERATIONS ) then

       output = SCF_INTRASPECIES_CONVERGENCE_FAILED

       call SingleSCF_exception(ERROR, "SCF_INTRASPECIES_CONVERGENCE_FAILED BY: "//trim(nameOfSpecieSelected)//" ITERATIONS EXCEEDED", &
            "Class object SCF in the testDensityMatrixChange function")

    else

       if ( SingleSCF_getNumberOfIterations( speciesID ) > 1 ) then

          !! Obtiene la desviacion de la ultima iteracion
          call List_end( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )
          deltaDensityMatrix= List_current( WaveFunction_instance(speciesID)%standartDesviationOfDensityMatrixElements )

          if( abs( deltaDensityMatrix -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) > tolerace ) then

             output = SCF_INTRASPECIES_CONVERGENCE_CONTINUE

          else

             output = SCF_INTRASPECIES_CONVERGENCE_SUCCESS

          end if

       else

          output = SCF_INTRASPECIES_CONVERGENCE_CONTINUE

       end if

    end if

  end function SingleSCF_testDensityMatrixChange

  !>
  !!  @brief Aplica el metodo de convegencia SCF especificado para la matriz de Fock
  subroutine SingleSCF_convergenceMethod( speciesID )
    implicit none
    integer, intent(in) :: speciesID

    !!**************************************************************************************************
    !! Emplea metodo de convergecia SCF
    if (  WaveFunction_instance(speciesID)%numberOfIterations  > 1  )then

       if ( Convergence_isInstanced( WaveFunction_instance(speciesID)%convergenceMethod ) ) then

          if (CONTROL_instance%CONVERGENCE_METHOD == SCF_CONVERGENCE_DAMPING ) then

             call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
                  WaveFunction_instance(speciesID)%fockMatrix, WaveFunction_instance(speciesID)%densityMatrix, &
                  WaveFunction_instance(speciesID)%OverlapMatrix, &
                  methodType=SCF_CONVERGENCE_DAMPING, &
                  coefficientMatrix=WaveFunction_instance(speciesID)%waveFunctionCoefficients, speciesID=speciesID )

          else if (CONTROL_instance%CONVERGENCE_METHOD == SCF_CONVERGENCE_DIIS ) then

             call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
                  WaveFunction_instance(speciesID)%fockMatrix, WaveFunction_instance(speciesID)%densityMatrix, &
                  WaveFunction_instance(speciesID)%overlapMatrix, &
                  methodType=SCF_CONVERGENCE_DIIS, &
                  coefficientMatrix=WaveFunction_instance(speciesID)%waveFunctionCoefficients, speciesID=speciesID )

          else if (CONTROL_instance%CONVERGENCE_METHOD == SCF_CONVERGENCE_MIXED ) then

             call Convergence_setMethod( WaveFunction_instance(speciesID)%convergenceMethod, &
                  WaveFunction_instance(speciesID)%fockMatrix, WaveFunction_instance(speciesID)%densityMatrix, &
                  WaveFunction_instance(speciesID)%overlapMatrix, &
                  methodType=SCF_CONVERGENCE_MIXED, &
                  coefficientMatrix=WaveFunction_instance(speciesID)%waveFunctionCoefficients)

          end if


          call Convergence_run( WaveFunction_instance(speciesID)%convergenceMethod )


       end if

    else

       call Convergence_setInitialDensityMatrix( WaveFunction_instance(speciesID)%convergenceMethod, &
            WaveFunction_instance(speciesID)%densityMatrix )

       call Convergence_setInitialFockMatrix( WaveFunction_instance(speciesID)%convergenceMethod, &
            WaveFunction_instance(speciesID)%fockMatrix )

    end if
    !!
    !!**************************************************************************************************

  end subroutine SingleSCF_convergenceMethod

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine SingleSCF_exception( typeMessage, description, debugDescription)
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

  end subroutine SingleSCF_exception

end module SingleSCF_

