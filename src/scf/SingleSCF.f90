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
!!   - <tt> 2023-03-15 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
!!        -# Reescribe el modulo para simplificarlo
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
  use OrbitalLocalizer_
  implicit none

  !< enum Matrix_type {
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_FAILED =  0
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_CONTINUE =  1
  integer, parameter, public :: SCF_INTRASPECIES_CONVERGENCE_SUCCESS =  2
  !< }

  public :: &
       SingleSCF_showIterationStatus, &
       SingleSCF_getNumberOfIterations,&
       SingleSCF_getLastEnergy, &
       SingleSCF_getStandardDeviationOfDensityMatrix,&
       SingleSCF_getDiisError, &
       SingleSCF_iterate, &
       SingleSCF_restart, &
       SingleSCF_reset

contains

  !>
  !! @brief Retorna el numero de iteraciones realizadas para una especie especificada
  subroutine SingleSCF_showIterationStatus(wfObject, status )
    implicit none
    type(WaveFunction) :: wfObject
    integer :: status

    select case(status)

    case(SCF_INTRASPECIES_CONVERGENCE_FAILED)

       print *,""
       print *,"SCF STATUS: CONVERGENCE FAILED BY: ",trim(wfObject%name)
       print *,""

    case(SCF_INTRASPECIES_CONVERGENCE_CONTINUE)

       print *,""
       print *,"SCF STATUS: CONVERGENCE CONTINUE BY: ",trim(wfObject%name)
       print *,""

    case(SCF_INTRASPECIES_CONVERGENCE_SUCCESS)

       print *,""
       print *,"SCF STATUS: CONVERGENCE SUCCESS BY: ",trim(wfObject%name)
       print *,""

    end select

  end subroutine SingleSCF_showIterationStatus


  !>
  !! @brief Retorna el numero de iteraciones realizadas para una especie especificada
  function SingleSCF_getNumberOfIterations(wfObject) result( output )
    implicit none
    Type(WaveFunction) :: wfObject
    integer :: output
    output = wfObject%numberOfIterations
  end function SingleSCF_getNumberOfIterations

  !>
  !! @brief Retorna la energia calculada en la ultima iteracion
  function SingleSCF_getLastEnergy(wfObject) result( output )
    implicit none
    Type(WaveFunction) :: wfObject

    real(8) :: output

    !! Obtiene el valor de energia de la ultima iteracion
    call List_end(wfObject%energySCF )
    output= List_current(wfObject%energySCF )

  end function SingleSCF_getLastEnergy

  !>
  !! @brief Retorna de la densviacion estandar de la matrix de densidad en la ultima iteracion
  function SingleSCF_getStandardDeviationOfDensityMatrix(wfObject) result( output )
    implicit none
    Type(WaveFunction) :: wfObject

    real(8) :: output

    call List_end(wfObject%standardDesviationOfDensityMatrixElements )
    output= List_current(wfObject%standardDesviationOfDensityMatrixElements )

  end function SingleSCF_getStandardDeviationOfDensityMatrix

  !>
  !! @brief Retorna de la densviacion estandar de la matrix de densidad en la ultima iteracion
  function SingleSCF_getDiisError(wfObject) result( output )
    implicit none
    Type(WaveFunction) :: wfObject
    real(8) :: output

    call List_end(wfObject%diisError )
    output= List_current(wfObject%diisError )

  end function SingleSCF_getDiisError

  !>
  !! @brief Realiza iteracion SCF
  !! @param wfObject es el objeto de la clase funcion de onda con las matrices del SCF
  subroutine SingleSCF_iterate(wfObject)
    implicit none
    type(WaveFunction) :: wfObject

    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: previousWavefunctionCoefficients, matchingMatrix, matchingMatrix2
    type(Vector) :: deltaVector 
    type(Matrix) :: auxOverlapMatrix

    integer(8) :: total3
    integer(8) :: numberOfContractions
    real(8) :: threshold, hold
    real(8) :: levelShiftingFactor, normCheck    
    integer :: ocupados, occupatedOrbitals, prodigals
    integer :: totales
    integer :: search, trial
    integer :: ii, jj, astrayOrbitals, startIteration
    integer :: i,j, mu, nu, index
    logical :: existFile

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    !    wfnFile = trim(CONTROL_instance%INPUT_FILE)//"lowdin.vec"
    wfnUnit = 30
    numberOfContractions = MolecularSystem_getTotalnumberOfContractions(wfObject%species )

    !!**********************************************************************************************

    !!Does this do anything?
    if ( CONTROL_instance%ELECTRONIC_WAVEFUNCTION_ANALYSIS ) then
       !! Warning:
       !! This part has not been tested
       call WaveFunction_buildDensityMatrix(wfObject)
       !! Calcula energia total para la especie especificada
       call WaveFunction_obtainTotalEnergyForSpecies(wfObject)
       call List_push_back(wfObject%energySCF, wfObject%totalEnergyForSpecies )
       return       
    end if

    !!**************************************************************************************************
    !! Modify the fock matrix according to the convergence method chosen
    call Convergence_setMethod( wfObject%convergenceMethod, &
         wfObject%fockMatrix, &
         wfObject%densityMatrix, &
         wfObject%overlapMatrix, &
         CONTROL_instance%CONVERGENCE_METHOD, &
         wfObject%waveFunctionCoefficients, &
         wfObject%species)
    call Convergence_run(wfObject%convergenceMethod )
    !!
    !!**************************************************************************************************

    call Matrix_copyConstructor( fockMatrixTransformed, wfObject%fockMatrix )

    if (CONTROL_instance%MO_FRACTION_OCCUPATION < 1.0_8 .or. CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF ) then
       call Matrix_copyConstructor( previousWavefunctionCoefficients, wfObject%waveFunctionCoefficients )
    end if

    !!**********************************************************************************************
    !! Iteration begins
    !!
    fockMatrixTransformed%values = &
         matmul( matmul( transpose(wfObject%transformationMatrix%values ) , &
         fockMatrixTransformed%values), wfObject%transformationMatrix%values )

    !! Calcula valores y vectores propios de matriz de Fock transformada.
    call Matrix_eigen( fockMatrixTransformed, wfObject%molecularOrbitalsEnergy, &
         wfObject%waveFunctionCoefficients, SYMMETRIC )

    !! Calcula los  vectores propios para matriz de Fock       
    wfObject%waveFunctionCoefficients%values = &
         matmul(wfObject%transformationMatrix%values, &
         wfObject%waveFunctionCoefficients%values )

    ! print *, "Coefficients before for ", wfObject%species
    ! call Matrix_show(wfObject%waveFunctionCoefficients)

    !! Filtra los orbitales eliminados por el umbral de solapamiento
    ! print *, "threshold", CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
    i=0
    do index = 1 , numberOfContractions
       i=i+1
       normCheck=0.0
       if ( abs(wfObject%molecularOrbitalsEnergy%values(i)) .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD ) then
          do mu = 1 , numberOfContractions
             do nu = 1 , numberOfContractions
                normCheck=normCheck+wfObject%waveFunctionCoefficients%values(mu,i)*&
                     wfObject%waveFunctionCoefficients%values(nu,i)*&
                     wfObject%overlapMatrix%values(mu,nu)
             end do
          end do
          ! print *, "eigenvalue", i, wfObject%molecularOrbitalsEnergy%values(i), "normCheck", normCheck

          if ( normCheck .lt. CONTROL_instance%OVERLAP_EIGEN_THRESHOLD) then
             ! Shift orbital coefficients to the end of the matrix and Make energy a very large number
             do j = i , numberOfContractions-1
                wfObject%molecularOrbitalsEnergy%values(j)=wfObject%molecularOrbitalsEnergy%values(j+1)
                wfObject%waveFunctionCoefficients%values(:,j) = wfObject%waveFunctionCoefficients%values(:,j+1)
             end do
             wfObject%molecularOrbitalsEnergy%values(numberOfContractions)=1/CONTROL_instance%OVERLAP_EIGEN_THRESHOLD
             wfObject%waveFunctionCoefficients%values(:,numberOfContractions)=0.0
             i=i-1
          end if

       end if
    end do

    if ( CONTROL_instance%ACTIVATE_LEVEL_SHIFTING) &
         call Convergence_removeLevelShifting(wfObject%convergenceMethod,wfObject%molecularOrbitalsEnergy)   
    !!
    !! Interation ends
    !!**********************************************************************************************

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *, "Coefficients after SCF cycle for ", wfObject%species
       call Matrix_show(wfObject%waveFunctionCoefficients)
    end if

    if (CONTROL_instance%MO_FRACTION_OCCUPATION < 1.0_8 .or. CONTROL_instance%EXCHANGE_ORBITALS_IN_SCF ) &
         call SingleSCF_orbitalExchange(wfObject,previousWavefunctionCoefficients)

    if ( CONTROL_instance%NO_SCF .or. CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS .or. CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS) &
         call SingleSCF_readFrozenOrbitals(wfObject)

    !! Exchanging orbitals just for calculation excited states
    if(wfObject%name == trim(CONTROL_instance%EXCITE_SPECIE) ) then

       call Matrix_constructor (auxiliaryMatrix, numberOfContractions, numberOfContractions)
       call Matrix_copyConstructor( auxiliaryMatrix, wfObject%waveFunctionCoefficients )
       wfObject%waveFunctionCoefficients%values(:,2)=auxiliaryMatrix%values(:,1)
       wfObject%waveFunctionCoefficients%values(:,1)=auxiliaryMatrix%values(:,2)
       call Matrix_destructor(auxiliaryMatrix)
       hold=wfObject%molecularOrbitalsEnergy%values(1)
       wfObject%molecularOrbitalsEnergy%values(1)=wfObject%molecularOrbitalsEnergy%values(2)
       wfObject%molecularOrbitalsEnergy%values(2)=hold

    end if


    !! Actualiza el contador para el numero de iteraciones SCF de la especie actual
    wfObject%numberOfIterations =  wfObject%numberOfIterations + 1
    call Matrix_destructor(fockMatrixTransformed)

    ! wfObject%wasBuiltFockMatrix = .false.

  end subroutine SingleSCF_iterate


  !>
  !! @brief Reinicia el proceso de iteracion SCF para la especie especificada
  subroutine SingleSCF_restart(wfObject)
    implicit none
    type(WaveFunction) :: wfObject

    wfObject%numberOfIterations = 0

    call List_clear(wfObject%energySCF )
    call List_clear(wfObject%standardDesviationOfDensityMatrixElements )
    call List_clear(wfObject%diisError )
    call Convergence_reset()

  end subroutine SingleSCF_restart

  !>
  !! @brief Reinicializa el modulo
  subroutine SingleSCF_reset(wfObject)
    implicit none
    type(WaveFunction) :: wfObject

    wfObject%numberOfIterations = 0

    call List_clear(wfObject%energySCF )
    call List_clear(wfObject%diisError )
    call List_clear(wfObject%standardDesviationOfDensityMatrixElements )

    call Convergence_destructor(wfObject%convergenceMethod)
    call Convergence_constructor(wfObject%convergenceMethod, &
         wfObject%name,CONTROL_instance%CONVERGENCE_METHOD)

    call Convergence_reset()

  end subroutine SingleSCF_reset


  !>
  !! @brief Updates density matrix after one SCF iteration
  ! subroutine SingleSCF_actualizeDensityMatrix(wfObject%name )
  !   implicit none
  !   character(*), optional :: wfObject%name

  !   integer :: wfObject%species
  !   character(30) :: wfObject%name

  !   wfObject%name = "E-"
  !   if ( present(wfObject%name ) )  wfObject%name= trim(wfObject%name )

  !   wfObject%species = MolecularSystem_getSpecieID(wfObject%name=trim(wfObject%name ) )

  !   !! Determina la desviacion estandar de los elementos de la matriz de densidad
  !   call Matrix_copyConstructor(wfObject%beforeDensityMatrix, wfObject%densityMatrix )

  !   call WaveFunction_builtDensityMatrix(wfObject%name )

  !   call List_push_back(wfObject%standardDesviationOfDensityMatrixElements, &
  !        Matrix_standardDeviation(wfObject%beforeDensityMatrix, wfObject%densityMatrix) )

  !   !! Calcula energia total para la especie especificada
  !   call WaveFunction_obtainTotalEnergyForSpecies(wfObject%name )

  !   call List_push_back(wfObject%energySCF, wfObject%totalEnergyForSpecies )

  !   call List_push_back(wfObject%diisError, Convergence_getDiisError(wfObject%convergenceMethod) )

  ! end subroutine SingleSCF_actualizeDensityMatrix

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
  
  !! @brief  Orbital exchange for TOM calculation
  !!
  subroutine SingleSCF_orbitalExchange(wfObject,previousWavefunctionCoefficients)
    type(WaveFunction) :: wfObject
    type(Matrix) :: previousWavefunctionCoefficients
    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: matchingMatrix, matchingMatrix2
    type(Vector) :: deltaVector 
    type(Matrix) :: auxOverlapMatrix

    integer(8) :: total3
    integer(8) :: numberOfContractions
    real(8) :: threshold, hold
    real(8) :: levelShiftingFactor, normCheck    
    integer :: ocupados, occupatedOrbitals, prodigals
    integer :: totales
    integer :: search, trial
    integer :: ii, jj, astrayOrbitals, startIteration
    integer :: i,j, mu, nu, index
    logical :: existFile

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    startIteration=1
    threshold=0.8_8

    call Matrix_copyConstructor(auxOverlapMatrix,wfObject%overlapMatrix)

    occupatedOrbitals = MolecularSystem_getOcupationNumber(wfObject%species )
    total3 = int(occupatedOrbitals,8)

    call Matrix_constructor (matchingMatrix, total3, total3)

    matchingMatrix%values = &
         matmul( matmul( transpose( previousWavefunctionCoefficients%values ) , &
         auxOverlapMatrix%values), wfObject%waveFunctionCoefficients%values )

    astrayOrbitals=0 !!! number of orbitals that must be reordered

    do ii= 1, occupatedOrbitals

       !DEBUG
       ! print *,"Antes ","Orbital: ",ii," ",abs(matchingMatrix%values(ii,ii))," energy ",wfObject%molecularOrbitalsEnergy%values(ii)

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

    if(astrayOrbitals .ge. 1 .and. wfObject%numberOfIterations > startIteration) then
       ! print *, "Switching orbitals..."

       do ii=1,astrayOrbitals

          call Matrix_copyConstructor( auxiliaryMatrix, wfObject%waveFunctionCoefficients )            
          call Matrix_copyConstructor( matchingMatrix2, matchingMatrix )            

          search=deltaVector%values(ii)            
          ! print *,"search: ",search

          do jj=1, numberOfContractions

             !trial=deltaVector%values(jj)
             trial=jj
             ! print *,"trial: ",trial

             ! print *,"valor: ",abs(matchingMatrix%values(search,trial))

             if ( abs(matchingMatrix%values(search,trial)) > threshold ) then

                wfObject%waveFunctionCoefficients%values(:,search)=auxiliaryMatrix%values(:,trial)
                wfObject%waveFunctionCoefficients%values(:,trial)=auxiliaryMatrix%values(:,search)

                hold=wfObject%molecularOrbitalsEnergy%values(search)
                wfObject%molecularOrbitalsEnergy%values(search)=wfObject%molecularOrbitalsEnergy%values(trial)
                wfObject%molecularOrbitalsEnergy%values(trial)=hold

                matchingMatrix%values(:,trial)=matchingMatrix2%values(:,search)
                matchingMatrix%values(:,search)=matchingMatrix2%values(:,trial)

                prodigals = prodigals - 1

                print *, "Switching orbital... ",search," with ",trial, " for ", trim(wfObject%name)

                exit

             end if

          end do

          ! if (prodigals<2) then
          !    exit
          ! end if

       end do

       if (  CONTROL_instance%DEBUG_SCFS) then
          print *, "Coefficients after switching for ", wfObject%species
          call Matrix_show(wfObject%waveFunctionCoefficients)
       end if

       call Matrix_destructor(matchingMatrix)
       call Vector_destructor(deltaVector)           
    end if
  end subroutine SingleSCF_OrbitalExchange

  !!**********************************************************************************************
  !! @brief If NO SCF cicle for a species is desired, read the coefficients from the ".vec" file again
  subroutine SingleSCF_readFrozenOrbitals(wfObject)
    type(WaveFunction) :: wfObject
    
    type(Matrix) :: fockMatrixTransformed
    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: previousWavefunctionCoefficients, matchingMatrix, matchingMatrix2
    type(Vector) :: deltaVector 
    type(Matrix) :: auxOverlapMatrix

    integer(8) :: total3
    integer(8) :: numberOfContractions
    real(8) :: threshold, hold
    real(8) :: levelShiftingFactor, normCheck    
    integer :: ocupados, occupatedOrbitals, prodigals
    integer :: totales
    integer :: search, trial
    integer :: ii, jj, astrayOrbitals, startIteration
    integer :: i,j, mu, nu, index
    logical :: existFile

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    wfnUnit = 30
    numberOfContractions = MolecularSystem_getTotalnumberOfContractions(wfObject%species )
    !! NO SCF cicle for electrons or non-electrons
    if ( CONTROL_instance%FREEZE_ELECTRONIC_ORBITALS .and. .not. MolecularSystem_instance%species(wfObject%species)%isElectron) return
    if ( CONTROL_instance%FREEZE_NON_ELECTRONIC_ORBITALS .and. MolecularSystem_instance%species(wfObject%species)%isElectron ) return

    !! Read coefficients from various possible files
    if (CONTROL_instance%READ_FCHK) then

       call Matrix_constructor (auxiliaryMatrix, numberOfContractions, numberOfContractions)
       call MolecularSystem_readFchk( trim(CONTROL_instance%INPUT_FILE)//trim(wfObject%name)//".fchk", &
            wfObject%waveFunctionCoefficients, auxiliaryMatrix, wfObject%name )
       call Matrix_destructor(auxiliaryMatrix)

    else if (CONTROL_instance%READ_COEFFICIENTS) then
       arguments(2) = MolecularSystem_getNameOfSpecie(wfObject%species)
       arguments(1) = "COEFFICIENTS"

       wfnFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
       inquire(FILE = wfnFile, EXIST = existFile )

       if ( existFile) then
          open(unit=wfnUnit, file=trim(wfnFile), status="old", form="formatted")

          wfObject%waveFunctionCoefficients = Matrix_getFromFile(unit=wfnUnit, &
               rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.false.,  & 
               arguments=arguments(1:2))

          close(wfnUnit)

       else 
          wfnFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
          inquire(FILE = wfnFile, EXIST = existFile )

          if ( existFile) then
             open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

             wfObject%waveFunctionCoefficients = Matrix_getFromFile(unit=wfnUnit, &
                  rows= int(numberOfContractions,4), columns= int(numberOfContractions,4), binary=.true., & 
                  arguments=arguments(1:2))

             close(wfnUnit)

          else
             call  SingleSCF_exception( ERROR, "I did not find any .vec coefficients file", "At SCF program, at SingleSCF_Iterate")
          end if

       end if

    else
       call  SingleSCF_exception( ERROR, "I did not find any coefficients file for the noSCF procedure", "At SCF program, at SingleSCF_Iterate")
    end if
    
    !! Calculate orbital energies from density matrix contributions
    wfObject%molecularOrbitalsEnergy%values=0.0
    !molecular orbital energy
    do i=1, numberOfContractions
       do mu=1, numberOfContractions
          do nu=1, numberOfContractions
             wfObject%molecularOrbitalsEnergy%values(i)=&
                  wfObject%molecularOrbitalsEnergy%values(i)+&
                  wfObject%fockMatrix%values(mu,nu)*&
                  wfObject%waveFunctionCoefficients%values(mu,i)*&
                  wfObject%waveFunctionCoefficients%values(nu,i)
          end do
       end do
    end do

  end subroutine SingleSCF_readFrozenOrbitals
end module SingleSCF_

