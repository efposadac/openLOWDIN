!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	PROF. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module handles all matrices for SCF program.
!! @author E. F. Posada, 20130
!! @warning This module is differs from the Wavefunction.f90 located in HF program.
module WaveFunction_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use List_
  use Convergence_
  use MolecularSystem_
  use CosmoCore_

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

     character(30) :: name

     !!**************************************************************
     !! Matrices requeridas y alteradas en la realizacion del ciclo SCF
     !!
     type(Matrix) :: overlapMatrix     
     type(Matrix) :: fockMatrix     
     type(Matrix) :: densityMatrix
     type(Matrix) :: hcoreMatrix
     type(Matrix) :: twoParticlesMatrix
     type(Matrix) :: couplingMatrix
     type(Matrix) :: interParticleCorrMatrix
     type(Matrix) :: externalPotentialMatrix
     type(Matrix) :: beforeDensityMatrix
     type(Matrix) :: transformationMatrix
     type(Matrix) :: waveFunctionCoefficients
     type(Vector) :: molecularOrbitalsEnergy     

     !! Cosmo Things

     type(Matrix) :: cosmo1
     type(Matrix) :: cosmo2
     type(Matrix) :: cosmo4
     type(Matrix) :: cosmoCoupling
     real(8) :: cosmoChargeValue

     !!
     !!**************************************************************

     logical :: wasBuiltFockMatrix

     !!**************************************************************
     !!  Variables y objetos asociados al metodo SCF
     !!
     integer :: numberOfIterations
     type(List) :: energySCF
     type(List) :: standartDesviationOfDensityMatrixElements
     type(List) :: diisError
     type(Convergence) :: convergenceMethod

     !!**************************************************************
     !! Variable por conveniencia
     real(8) :: totalEnergyForSpecie
     real(8) :: independentSpecieEnergy
     real(8) :: nuclearElectronicCorrelationEnergy

  end type WaveFunction

  type(WaveFunction), public, allocatable :: WaveFunction_instance(:)

contains

  !>
  !! @brief Define el constructor para la clase
  subroutine WaveFunction_constructor( wfnUnit )
    implicit none

    integer, intent(in) :: wfnUnit

    integer :: speciesID    
    integer(8) :: numberOfContractions
    character(50) :: labels(2)

    !! Allocate memory.
    allocate(WaveFunction_instance(MolecularSystem_instance%numberOfQuantumSpecies))

    !! Allocate memory for specie in system and load some matrices.
    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       labels = ""
       labels(2) = trim(MolecularSystem_getNameOfSpecie(speciesID))
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

       !! Parametros Asociados con el SCF
       call List_constructor( WaveFunction_instance( speciesID )%energySCF,"energy",CONTROL_instance%LISTS_SIZE )
       call List_constructor( WaveFunction_instance( speciesID )%diisError,"diisError",CONTROL_instance%LISTS_SIZE )
       call List_constructor( WaveFunction_instance( speciesID )%standartDesviationOfDensityMatrixElements, "densitySD",CONTROL_instance%LISTS_SIZE )

       !! Instancia un objeto para manejo de aceleracion y convergencia del metodo SCF
       call Convergence_constructor(WaveFunction_instance( speciesID )%convergenceMethod, &
            WaveFunction_instance( speciesID )%name,CONTROL_instance%CONVERGENCE_METHOD)

       !! Set defaults
       WaveFunction_instance( speciesID )%totalEnergyForSpecie = 0.0_8
       WaveFunction_instance( speciesID )%independentSpecieEnergy =0.0_8
       WaveFunction_instance( speciesID )%nuclearElectronicCorrelationEnergy = 0.0_8
       WaveFunction_instance( speciesID )%numberOfIterations = 0 

       !! Cosmo things
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo1, numberOfContractions, numberOfContractions, 0.0_8 )     
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo4,numberOfContractions, numberOfContractions, 0.0_8 )

       !! Load integrals form lowdin.wfn
       labels(1) = "OVERLAP"
       WaveFunction_instance(speciesID)%overlapMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=labels)

       labels(1) = "HCORE"
       WaveFunction_instance(speciesID)%HcoreMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=labels)

       labels(1) = "DENSITY"
       WaveFunction_instance(speciesID)%densityMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=labels)

       labels(1) = "TRANSFORMATION"
       WaveFunction_instance(speciesID)%transformationMatrix = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=labels)

       !! Cosmo things

       if (CONTROL_instance%COSMO) then


          labels(1) = "COSMO1"
          WaveFunction_instance(speciesID)%cosmo1 = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
               columns= int(numberOfContractions,4), binary=.true., arguments=labels)
          labels(1) = "COSMO4"
          WaveFunction_instance(speciesID)%cosmo4 = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
               columns= int(numberOfContractions,4), binary=.true., arguments=labels)

          WaveFunction_instance(speciesID)%cosmoChargeValue=0.0_8

       end if


       if (  CONTROL_instance%DEBUG_SCFS) then
          write(*,*) "Matriz cosmo1 "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%cosmo1)
          write(*,*) "Matriz cosmo4 "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%cosmo4)
          write(*,*) "Matriz de Overlap "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%overlapMatrix)
          write(*,*) "Matriz de Hcore "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%hcoreMatrix)
          write(*,*) "Matriz de Densidad "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%densityMatrix)
          write(*,*) "Matriz de Transformacion "//trim(MolecularSystem_getNameOfSpecie(speciesID))
          call Matrix_show(WaveFunction_instance(speciesID)%transformationMatrix)          
       end if

       !! Build some matrices
       call Matrix_constructor( WaveFunction_instance(speciesID)%fockMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%twoParticlesMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%couplingMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%interParticleCorrMatrix, numberOfContractions, numberOfContractions, 0.0_8 )       
       call Matrix_constructor( WaveFunction_instance(speciesID)%externalPotentialMatrix, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%waveFunctionCoefficients,numberOfContractions, numberOfContractions, 0.0_8 )
       call Vector_constructor( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy, int(numberOfContractions) )

       !!cosmo things
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmo2, numberOfContractions, numberOfContractions, 0.0_8 )
       call Matrix_constructor( WaveFunction_instance(speciesID)%cosmoCoupling, numberOfContractions, numberOfContractions, 0.0_8 )

       WaveFunction_instance(speciesID)%wasBuiltFockMatrix = .false.

    end do

  end subroutine WaveFunction_constructor

  !>
  !! @brief Builds two-particles matrix.
  subroutine WaveFunction_buildTwoParticlesMatrix( nameOfSpecie, nproc )
    implicit none

    character(*), optional :: nameOfSpecie
    integer :: nproc

    character(30) :: nameOfSpecieSelected

    real(8) :: integralValue
    real(8) :: coulomb,exchange
    real(8) :: factor

    integer, target :: numberOfContractions
    integer, target :: a, b, r, s
    integer :: speciesID
    integer :: n, u, v, i
    integer :: status

    integer*2 :: aa(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: bb(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: rr(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: ss(CONTROL_instance%INTEGRAL_STACK_SIZE)

    real(8) :: shellIntegrals(CONTROL_instance%INTEGRAL_STACK_SIZE)

    logical :: hadAdded

    integer :: ifile
    integer :: unit
    character(50) :: sfile

    real(8), allocatable :: tmpArray(:,:)

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )


    !! This matrix is only calculated if there are more than one particle for speciesID or if the user want to calculate it.
    if ( MolecularSystem_getNumberOfParticles( speciesID ) > 1 .or. CONTROL_instance%BUILD_TWO_PARTICLES_MATRIX_FOR_ONE_PARTICLE ) then

       wavefunction_instance(speciesID)%twoParticlesMatrix%values = 0.0_8
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
       factor = MolecularSystem_getFactorOfInterchangeIntegrals( speciesID )       

       !$OMP PARALLEL private(ifile,sfile,unit,aa,bb,rr,ss,shellIntegrals,i,coulomb,exchange, tmpArray), shared(wavefunction_instance)
       !$OMP DO 
       do ifile = 1, nproc

          write(sfile,*) ifile
          sfile = trim(adjustl(sfile))
          unit = ifile+50

          !! open file (order, integral(shell))
          if(CONTROL_instance%IS_OPEN_SHELL .and. MolecularSystem_instance%species(speciesID)%isElectron) then

             open( UNIT=unit,FILE=trim(sfile)//"E-ALPHA.ints", status='old',access='sequential', form='Unformatted')

          else

             open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".ints", status='old',access='sequential', form='Unformatted')

          end if

          if(allocated(tmpArray))deallocate(tmpArray)
          allocate(tmpArray(numberOfContractions,numberOfContractions))
          tmpArray = 0.0_8

          loadintegrals : do

             read(UNIT=unit, iostat=status) aa(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  bb(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  rr(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  ss(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                  shellIntegrals(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

             if(status == -1 ) then
                print*, "end of file! file: ",trim(sfile)//"E-ALPHA.ints"
                exit loadintegrals
             end if

             do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                if( aa(i) == -1 ) exit loadintegrals

                coulomb = wavefunction_instance(speciesID)%densityMatrix%values(rr(i),ss(i))*shellIntegrals(i)

                !!*****************************************************************************
                !! Adds coulomb operator contributions
                if( aa(i) == rr(i) .and. bb(i) == ss(i) ) then

                   tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   if( rr(i) /= ss(i) ) then

                      tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   end if

                else

                   tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   if( rr(i) /= ss(i) ) then

                      tmpArray(aa(i),bb(i)) = tmpArray(aa(i),bb(i)) + coulomb

                   end if

                   coulomb = wavefunction_instance(speciesID)%densityMatrix%values(aa(i),bb(i))*shellIntegrals(i)

                   tmpArray(rr(i),ss(i)) = tmpArray(rr(i),ss(i)) + coulomb

                   if ( aa(i) /= bb(i) ) then

                      tmpArray( rr(i), ss(i) ) = tmpArray( rr(i), ss(i) ) + coulomb

                   end if

                end if

                !!
                !!*****************************************************************************

                !!*****************************************************************************
                !! Adds exchange operator contributions
                if( rr(i) /= ss(i) ) then

                   exchange =wavefunction_instance(speciesID)%densityMatrix%values(bb(i),ss(i))*shellIntegrals(i)* factor

                   tmpArray( aa(i), rr(i) ) = tmpArray( aa(i), rr(i) ) + exchange

                   if( aa(i) == rr(i) .and. bb(i) /= ss(i) ) then

                      tmpArray( aa(i), rr(i) ) = tmpArray( aa(i), rr(i) ) + exchange

                   end if

                end if

                if ( aa(i) /= bb(i) ) then

                   exchange = wavefunction_instance(speciesID)%densityMatrix%values(aa(i),rr(i))*shellIntegrals(i) * factor

                   if( bb(i) > ss(i) ) then

                      tmpArray( ss(i), bb(i) ) = tmpArray( ss(i), bb(i)) + exchange

                   else

                      tmpArray( bb(i), ss(i) ) = tmpArray( bb(i), ss(i) ) + exchange

                      if( bb(i)==ss(i) .and. aa(i) /= rr(i) ) then

                         tmpArray( bb(i), ss(i) ) = tmpArray( bb(i), ss(i) ) + exchange

                      end if

                   end if

                   if ( rr(i) /= ss(i) ) then

                      exchange = wavefunction_instance(speciesID)%densityMatrix%values(aa(i),ss(i))*shellIntegrals(i) * factor

                      if( bb(i) <= rr(i) ) then

                         tmpArray( bb(i), rr(i) ) = tmpArray( bb(i), rr(i) ) + exchange

                         if( bb(i) == rr(i) ) then

                            tmpArray( bb(i), rr(i) ) = tmpArray( bb(i), rr(i) ) + exchange

                         end if

                      else

                         tmpArray( rr(i), bb(i) ) = tmpArray( rr(i), bb(i)) + exchange

                         if( aa(i) == rr(i) .and. ss(i) == bb(i) ) goto 30

                      end if

                   end if

                end if

                exchange = wavefunction_instance(speciesID)%densityMatrix%values(bb(i),rr(i))*shellIntegrals(i) * factor

                tmpArray( aa(i), ss(i) ) = tmpArray( aa(i), ss(i) ) + exchange

30              continue

                !!
                !!*****************************************************************************

             end do
          end do loadintegrals

          close(unit)

          do u = 1, numberOfContractions
             do v = 1, numberOfContractions
                !$OMP ATOMIC
                wavefunction_instance(speciesID)%twoParticlesMatrix%values(u,v) = &
                     wavefunction_instance(speciesID)%twoParticlesMatrix%values(u,v) + tmpArray(u,v) 
             end do
          end do

       end do
       !$OMP END DO
       !$OMP END PARALLEL

       do u = 1 , numberOfContractions
          do v = u , numberOfContractions

             wavefunction_instance(speciesID)%twoParticlesMatrix%values(v,u) = wavefunction_instance(speciesID)%twoParticlesMatrix%values(u,v)

          end do
       end do

       wavefunction_instance(speciesID)%twoParticlesMatrix%values = wavefunction_instance(speciesID)%twoParticlesMatrix%values * ( MolecularSystem_getCharge(speciesID=speciesID ) )**2.0_8

    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "two particle matrix for: ", trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%twoParticlesMatrix)
    end if

  end subroutine WaveFunction_buildTwoParticlesMatrix

  !>
  !! @brief Builds the coupling matrix.
  subroutine WaveFunction_buildCouplingMatrix( nameOfSpecie )
    implicit none

    character(*), optional :: nameOfSpecie

    character(30) :: nameOfSpecieSelected
    character(30) :: nameOfOtherSpecie
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: currentSpecieID
    integer :: otherSpecieID
    integer :: speciesIterator
    integer :: ssize
    integer :: i, j, u
    real(8), allocatable :: auxMatrix(:,:)
    real(8) :: coulomb

    integer*2 :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)    

    nameOfSpecieSelected = "E-"    
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    currentSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpecieID)

    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

       ssize = size(wavefunction_instance(currentSpecieID)%couplingMatrix%values,dim=1)

       allocate(auxMatrix(ssize, ssize))
       auxMatrix=0.0_8                
       wavefunction_instance(currentSpecieID)%couplingMatrix%values = 0.0_8

       do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies()

          otherSpecieID = speciesIterator
          nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )          
          OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpecieID /= currentSpecieID ) then


             if( currentSpecieID > otherSpecieID) then       
                auxMatrix = 0.0_8

                !! open file for integrals
                open(UNIT=34,FILE=trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)//".ints", &
                     STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

                readIntegrals1 : do

                   read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE                      

                      if (a(u) == -1) exit readIntegrals1

                      coulomb = wavefunction_instance(otherSpecieID)%densityMatrix%values(a(u),b(u))*integral(u)

                      auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

                      if( a(u) /= b(u) ) auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

                   end do

                end do readIntegrals1

                close(34)

                auxMatrix= auxMatrix * MolecularSystem_getCharge(currentSpecieID ) * MolecularSystem_getCharge( otherSpecieID )
                wavefunction_instance(currentSpecieID)%couplingMatrix%values = wavefunction_instance(currentSpecieID)%couplingMatrix%values + auxMatrix

             else
                auxMatrix=0.0_8

                !! open file for integrals
                open(UNIT=34,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
                     STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

                readIntegrals2 : do

                   read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                      if (a(u) == -1) exit readIntegrals2

                      coulomb = wavefunction_instance(otherSpecieID)%densityMatrix%values(r(u),s(u))*integral(u)

                      auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

                      if( r(u) /= s(u) ) auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

                   end do

                end do readIntegrals2

                close(34)

                auxMatrix= auxMatrix * MolecularSystem_getCharge(currentSpecieID ) * MolecularSystem_getCharge( otherSpecieID )
                wavefunction_instance(currentSpecieID)%couplingMatrix%values = wavefunction_instance(currentSpecieID)%couplingMatrix%values + auxMatrix

             end if

          end if

       end do

       deallocate(auxMatrix)

       !! Simetrize
       do i = 1 , ssize
          do j = i , ssize
             wavefunction_instance(currentSpecieID)%couplingMatrix%values(j,i) = wavefunction_instance(currentSpecieID)%couplingMatrix%values(i,j)
          end do
       end do

    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "Matriz de acoplamiento: ", trim(nameOfSpecieSelected)
       call Matrix_show( wavefunction_instance(currentSpecieID)%couplingMatrix )
    end if

  end subroutine WaveFunction_buildCouplingMatrix

  !>
  !! @brief Builds fock Matrix
  subroutine WaveFunction_buildFockMatrix( nameOfSpecie )       
    implicit none

    character(*), optional :: nameOfSpecie

    character(30) :: nameOfSpecieSelected
    integer :: speciesID
    type(Matrix)::cosmoContribution

    nameOfSpecieSelected = "E-"    
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%hcoreMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 1 (hcore): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    !! cosmo fock matrix

    !!full coupling
    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
         0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
         wavefunction_instance(speciesID)%cosmo4%values)+ &
         wavefunction_instance(speciesID)%cosmo2%values + &
         wavefunction_instance(speciesID)%cosmoCoupling%values 

    !!half coupling
    ! wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
    !      0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
    !      wavefunction_instance(speciesID)%cosmo4%values)+ &
    !      wavefunction_instance(speciesID)%cosmo2%values +0.5_8*( &
    !      wavefunction_instance(speciesID)%cosmoCoupling%values) 

    !!without coupling
    ! wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + &
    !      0.5_8*(wavefunction_instance(speciesID)%cosmo1%values + &
    !      wavefunction_instance(speciesID)%cosmo4%values)+ &
    !      wavefunction_instance(speciesID)%cosmo2%values


    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + wavefunction_instance(speciesID)%twoParticlesMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 2 (+ two particles): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    wavefunction_instance(speciesID)%fockMatrix%values = wavefunction_instance(speciesID)%fockMatrix%values + wavefunction_instance(speciesID)%couplingMatrix%values

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK 3 (+ coupling): "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"MATRIZ DE FOCK: "//trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%fockMatrix)
    end if

    WaveFunction_instance(speciesID)%wasBuiltFockMatrix = .true.

  end subroutine WaveFunction_buildFockMatrix

  !>
  !! @brief Calcula la matriz de densidad para una especie especificada
  subroutine WaveFunction_builtDensityMatrix( nameOfSpecie )
    implicit none

    character(*), optional :: nameOfSpecie

    character(30) :: nameOfSpecieSelected
    integer :: orderMatrix
    integer :: speciesID
    integer :: ocupationNumber
    integer :: i
    integer :: j
    integer :: k

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( trim(nameOfSpecieSelected) )

    orderMatrix = size( wavefunction_instance(speciesID)%densityMatrix%values, DIM = 1 )

    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )

    wavefunction_instance(speciesID)%densityMatrix%values = 0.0_8

    !! Segment for fractional occupations: 1
    if (CONTROL_instance%IONIZE_MO /= 0 .and. trim(nameOfSpecieSelected) == trim(CONTROL_instance%IONIZE_SPECIE(1)) ) then
       wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO) = &
            wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO)*sqrt(CONTROL_instance%MO_FRACTION_OCCUPATION)

    end if

    do i = 1 , orderMatrix
       do j = 1 , orderMatrix
          do k = 1 , ocupationNumber

             wavefunction_instance(speciesID)%densityMatrix%values(i,j) =  &
                  wavefunction_instance(speciesID)%densityMatrix%values( i,j ) + &
                  ( wavefunction_instance(speciesID)%waveFunctionCoefficients%values(i,k) &
                  * wavefunction_instance(speciesID)%waveFunctionCoefficients%values(j,k) )
          end do
       end do
    end do

    wavefunction_instance(speciesID)%densityMatrix%values =  MolecularSystem_getEta( speciesID )  * wavefunction_instance(speciesID)%densityMatrix%values

    !! Segment for fractional occupations: 1
    if (CONTROL_instance%IONIZE_MO /= 0 .and. trim(nameOfSpecieSelected) == trim(CONTROL_instance%IONIZE_SPECIE(1))) then
       wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO) = &
            wavefunction_instance(speciesID)%waveFunctionCoefficients%values(:,CONTROL_instance%IONIZE_MO)/sqrt(CONTROL_instance%MO_FRACTION_OCCUPATION)
    end if

    !!DEBUG
    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Density Matrix ", trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%densityMatrix)
    end if

  end subroutine WaveFunction_builtDensityMatrix

  !>
  !! @brief Calculates total energy for one species
  subroutine WaveFunction_obtainTotalEnergyForSpecie( nameOfSpecie, nproc )

    implicit none

    character(*), optional :: nameOfSpecie
    integer :: nproc

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    !! Recalcula matriz de dos particulas (G) Con la nueva matriz de densidad
    call WaveFunction_buildTwoParticlesMatrix( trim(nameOfSpecieSelected) , nproc)

    if( .not. allocated(wavefunction_instance(speciesID)%externalPotentialMatrix%values) ) then

       wavefunction_instance(speciesID)%totalEnergyForSpecie = &
            sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
            *  (( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
            + 0.5_8 *wavefunction_instance(speciesID)%twoParticlesMatrix%values &
            + wavefunction_instance(speciesID)%couplingMatrix%values)) &
            + wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy


    else if(CONTROL_instance%COSMO)then

       wavefunction_instance(speciesID)%totalEnergyForSpecie = &
            sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
            *  (( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
            + 0.5_8 *wavefunction_instance(speciesID)%twoParticlesMatrix%values &
            + wavefunction_instance(speciesID)%couplingMatrix%values)) &
            + wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy


       wavefunction_instance( speciesID )%totalEnergyForSpecie =wavefunction_instance( speciesID )%totalEnergyForSpecie + 0.5_8 * &
            (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo1%values )+ &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo2%values ) + &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmo4%values ) + &
            sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
            wavefunction_instance( speciesID )%cosmoCoupling%values ))

       ! wavefunction_instance( speciesID )%totalEnergyForSpecie =wavefunction_instance( speciesID )%totalEnergyForSpecie + 0.5_8 * &
       !      (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
       !      wavefunction_instance( speciesID )%cosmo1%values )+ &
       !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
       !      wavefunction_instance( speciesID )%cosmo2%values ) + &
       !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
       !      wavefunction_instance( speciesID )%cosmo4%values )) 

    else

       wavefunction_instance(speciesID)%totalEnergyForSpecie = &
            sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
            *  (  ( wavefunction_instance(speciesID)%hcoreMatrix%values ) &
            + 0.5_8 *wavefunction_instance(speciesID)%twoParticlesMatrix%values &
            + wavefunction_instance(speciesID)%couplingMatrix%values &
            + wavefunction_instance(speciesID)%externalPotentialMatrix%values ))&
            + wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy

    end if

    if (  CONTROL_instance%DEBUG_SCFS) then
       print *,"Total energy for "// trim(nameOfSpecieSelected) //"= ", wavefunction_instance(speciesID)%totalEnergyForSpecie
    end if

  end subroutine WaveFunction_obtainTotalEnergyForSpecie

  !>
  !! @brief Calcula la energia total para el sistema estudiado  
  subroutine WaveFunction_obtainTotalEnergy( totalEnergy, totalCouplingEnergy,  electronicRepulsionEnergy, cosmo3Energy)
    implicit none

    real(8) :: totalEnergy
    real(8) :: totalCouplingEnergy
    real(8) :: electronicRepulsionEnergy

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    !! cosmo

    type(surfaceSegment) :: surface_aux2
    real(8) :: cosmo3Energy

    totalEnergy = 0.0_8
    cosmo3Energy = 0.0_8

    do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies

       !! Calula enegia de especie independiente ( sin considerar el termino de acoplamiento )
       if( .not. allocated( WaveFunction_instance( speciesID )%externalPotentialMatrix%values ) ) then

          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values))

       else if(CONTROL_instance%COSMO)then


          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values))


          wavefunction_instance( speciesID )%independentSpecieEnergy =wavefunction_instance( speciesID )%independentSpecieEnergy + 0.5_8 * &
               (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo1%values )+ &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo2%values ) +  &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmo4%values ) + &
               sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
               wavefunction_instance( speciesID )%cosmoCoupling%values))

          ! wavefunction_instance( speciesID )%independentSpecieEnergy =wavefunction_instance( speciesID )%independentSpecieEnergy + 0.5_8 * &
          !      (sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo1%values )+ &
          !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo2%values ) +  &
          !      sum( transpose( WaveFunction_instance( speciesID )%densityMatrix%values ) * &
          !      wavefunction_instance( speciesID )%cosmo4%values)) 


       else

          WaveFunction_instance( speciesID )%independentSpecieEnergy = &
               sum(  transpose(WaveFunction_instance( speciesID )%densityMatrix%values) &
               *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
               + 0.5_8 * WaveFunction_instance( speciesID )%twoParticlesMatrix%values &
               + WaveFunction_instance( speciesID )%externalPotentialMatrix%values))

       end if

       WaveFunction_instance( speciesID )%independentSpecieEnergy = &
            WaveFunction_instance( speciesID )%independentSpecieEnergy + &
            WaveFunction_instance( speciesID )%nuclearElectronicCorrelationEnergy

       totalEnergy = totalEnergy + WaveFunction_instance( speciesID )%independentSpecieEnergy

    end do

    !! Adicionado energia de interaccion entre particulas puntuales
    totalEnergy = totalEnergy + MolecularSystem_getPointChargesEnergy()


    !! cosmo potential nuclei-charges nuclei

    if(CONTROL_instance%COSMO)then
       call CosmoCore_lines(surface_aux2)
       call CosmoCore_filler(surface_aux2)

       call CosmoCore_nucleiPotentialNucleiCharges(surface_aux2,cosmo3Energy)
       totalEnergy=totalEnergy+cosmo3Energy

    end if

    !! Adicionar  energia de acoplamiento
    totalCouplingEnergy = WaveFunction_getTotalCouplingEnergy()

    !! Adds inter-electron species coupling energy
    electronicRepulsionEnergy = WaveFunction_getAlphaBetaRepulsion()

    !! Total Energy
    totalEnergy = totalEnergy +  totalCouplingEnergy + electronicRepulsionEnergy 

    ! print *,"Total energy = ", totalEnergy

  end subroutine WaveFunction_obtainTotalEnergy

  !>
  !! @brief calcula la energia total de acoplamiento para una especie especificada
  function WaveFunction_getTotalCouplingEnergy() result( output )
    implicit none
    real(8) :: output

    character(30) :: nameOfSpecie
    character(30) :: nameOfOtherSpecie
    real(8) :: auxValue
    real(8) :: auxRepulsion
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    integer :: numberOfTotalContractions
    integer :: numberOfTotalContractionsOfOtherSpecie
    integer :: speciesID
    integer :: otherSpecieID
    integer :: outFile
    integer*2 :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: u, v
    integer :: k, l, m
    integer :: arrayNumber

    output = 0.0_8

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
       do otherSpecieID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies()

          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpecieID /= speciesID ) then

             nameOfSpecie = MolecularSystem_getNameOfSpecie( speciesID )
             numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
             numberOfTotalContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

             nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )
             numberOfContractionsOfOtherSpecie = MolecularSystem_getNumberOfContractions( otherSpecieID )
             numberOfTotalContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( otherSpecieID )

             !Restringe la suma de terminos repulsivos electronicos
             if(trim(nameOfSpecie)=="E-ALPHA" .and. trim(nameOfOtherSpecie)=="E-BETA") cycle

             auxValue = 0.0_8
             m = 0

             !! open file for integrals
             open(UNIT=34,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
                  STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

             auxValue=0.0_8

             readIntegrals : do

                read(34) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                     r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                     integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                   if (a(u) == -1) exit readIntegrals

                   m = m + 1

                   auxValue = auxValue +&
                        (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                        * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                        *  integral(u))

                   if(b(u) /= a(u)) then

                      m = m + 1

                      auxValue = auxValue +&
                           (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                           * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                           *  integral(u))
                   end if

                   if(s(u) /= r(u)) then

                      m = m + 1

                      auxValue = auxValue +&
                           (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                           * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                           *  integral(u))
                   end if

                   if(b(u) /= a(u) .and. s(u) /= r(u)) then

                      m = m + 1

                      auxValue = auxValue +&
                           (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                           * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                           *  integral(u))
                   end if


                end do

             end do readIntegrals

             auxValue = auxValue *  MolecularSystem_getCharge( speciesID=speciesID ) &
                  * MolecularSystem_getCharge( speciesID=otherSpecieID )

             output = output + auxValue

             close(34)

          end if
       end do
    end do

  end function WaveFunction_getTotalCouplingEnergy

  !>
  !! @brief calcula la energia total de acoplamiento para una especie especificada
  function WaveFunction_getAlphaBetaRepulsion() result( output )
    implicit none

    real(8) :: output

    character(30) :: nameOfSpecie
    character(30) :: nameOfOtherSpecie
    real(8) :: auxValue
    real(8) :: auxRepulsion
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    integer :: numberOfTotalContractions
    integer :: numberOfTotalContractionsOfOtherSpecie
    integer :: speciesID
    integer :: otherSpecieID
    integer :: outFile
    integer*2 :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer*2 :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: u, v
    integer :: k, l, m
    integer :: arrayNumber

    output =0.0_8
    auxValue = 0.0_8

    do speciesID = 1, MolecularSystem_getNumberOfQuantumSpecies()
       do otherSpecieID = speciesID+1, MolecularSystem_getNumberOfQuantumSpecies()

          !! Restringe suma de terminos repulsivos de la misma especie.
          if ( otherSpecieID /= speciesID ) then

             nameOfSpecie = MolecularSystem_getNameOfSpecie( speciesID )
             numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
             numberOfTotalContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )

             nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )
             numberOfContractionsOfOtherSpecie = MolecularSystem_getNumberOfContractions( otherSpecieID )
             numberOfTotalContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( otherSpecieID )

             !Restringe la suma a solo electrones
             if(trim(nameOfSpecie)=="E-ALPHA" .and. trim(nameOfOtherSpecie)=="E-BETA") then

                auxValue = 0.0_8
                m = 0             
                !! open file for integrals
                open(UNIT=34,FILE=trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
                     STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

                readIntegrals : do

                   read(34) a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                        integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

                   do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

                      if (a(u) == -1) exit readIntegrals

                      m = m + 1

                      auxValue = auxValue +&
                           (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                           * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                           *  integral(u))

                      if(b(u) /= a(u)) then

                         m = m + 1

                         auxValue = auxValue +&
                              (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                              * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                              *  integral(u))
                      end if

                      if(s(u) /= r(u)) then

                         m = m + 1

                         auxValue = auxValue +&
                              (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                              * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                              *  integral(u))
                      end if

                      if(b(u) /= a(u) .and. s(u) /= r(u)) then

                         m = m + 1

                         auxValue = auxValue +&
                              (  wavefunction_instance(speciesID)%densityMatrix%values(b(u),a(u)) &
                              * WaveFunction_instance( otherSpecieID)%densityMatrix%values(r(u),s(u)) &
                              *  integral(u))
                      end if

                   end do

                end do readIntegrals

                auxValue = auxValue *  MolecularSystem_getCharge( speciesID=speciesID ) &
                     * MolecularSystem_getCharge( speciesID=otherSpecieID )

                output = output + auxValue

             end if
          end if

          close(34)                

       end do
    end do

  end function WaveFunction_getAlphaBetaRepulsion


  !   !>
  !   !! @brief   Contruye la matrix de acoplamineto para la especie especificada (C)
  !   !!      la matrix resultante sera tenida en cuenta en la construccion de la matriz de Fock
  !   !!
  !   !>
  !   subroutine WaveFunction_buildCouplingMatrixElectronFree( nameOfSpecie, output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(matrix) :: output

  !     character(30) :: nameOfSpecieSelected
  !     character(30) :: nameOfOtherSpecie
  !     integer :: numberOfContractions
  !     integer :: otherNumberOfContractions
  !     integer :: currentSpecieID
  !     integer :: otherSpecieID
  !     integer :: speciesIterator
  !     integer(8) :: ssize
  !     integer :: i, j
  !     integer :: k, l
  !     integer :: arrayNumber
  !     real(8), allocatable :: auxMatrix(:,:)
  !     real(8) :: coulomb

  !     integer*2 :: a(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer*2 :: b(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer*2 :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer*2 :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
  !     integer :: u

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     currentSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
  !     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpecieID)

  !     if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !        ssize = size(WaveFunction_instance( currentSpecieID )%wavefunction_instance(speciesID)%couplingMatrix%values,dim=1)
  !        allocate(auxMatrix(ssize, ssize))
  !        call Matrix_constructor(output, ssize, ssize)
  !        output%values=0.0_8

  !        do speciesIterator = MolecularSystem_beginSpecie(), MolecularSystem_endSpecie()

  !           otherSpecieID = MolecularSystem_getSpecieID( iteratorOfSpecie = speciesIterator )
  !           nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )

  !           if(trim(nameOfSpecieSelected) == "e-ALPHA" .and. trim(nameOfOtherSpecie) == "e-BETA" ) cycle
  !           if(trim(nameOfSpecieSelected) == "e-BETA" .and. trim(nameOfOtherSpecie) == "e-ALPHA" ) cycle

  !           OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

  !           !! Restringe suma de terminos repulsivos de la misma especie.
  !           if ( otherSpecieID /= currentSpecieID ) then

  !              !! ALL IN MEMORY
  !              if( .not. IntegralManager_instance%toDisk ) then

  !                 if( currentSpecieID > otherSpecieID) then

  !                    call IntegralManager_interspecieRepulsionIntegral (otherSpecieID, currentSpecieID , isInterSpecies=.true., arrayNumber=arrayNumber)

  !                    auxMatrix=0.0_8
  !                    u = 0

  !                    do i=1,otherNumberOfContractions
  !                       do j=i, otherNumberOfContractions
  !                          do k=1,numberOfContractions
  !                             do l=k,numberOfContractions

  !                                u = u + 1

  !                                coulomb = WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values(i,j) * &
  !                                     IntegralManager_instance%interSpecieRepulsionIntegrals(arrayNumber)%values(u)

  !                                auxMatrix(k,l) = auxMatrix(k,l) + coulomb

  !                                if( i /= j ) auxMatrix(k,l) = auxMatrix(k,l) + coulomb

  !                             end do
  !                          end do
  !                       end do
  !                    end do

  !                 else

  !                    auxMatrix=0.0_8
  !                    u = 0

  !                    call IntegralManager_interspecieRepulsionIntegral (currentSpecieID, otherSpecieID, isInterSpecies=.true., arrayNumber=arrayNumber)

  !                    do i=1,numberOfContractions
  !                       do j=i, numberOfContractions
  !                          do k=1,otherNumberOfContractions
  !                             do l=k,otherNumberOfContractions

  !                                u = u + 1

  !                                coulomb = WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values(k,l)* &
  !                                     IntegralManager_instance%interSpecieRepulsionIntegrals(arrayNumber)%values(u)

  !                                auxMatrix(i,j) = auxMatrix(i,j) + coulomb

  !                                if( k /= l ) auxMatrix(i,j) = auxMatrix(i,j) + coulomb

  !                             end do
  !                          end do
  !                       end do
  !                    end do

  !                 end if

  !                 !!ALL IN DISK
  !              else

  !                 if( currentSpecieID > otherSpecieID) then

  !                    call IntegralManager_interspecieRepulsionIntegral (otherSpecieID, currentSpecieID , isInterSpecies=.true.)

  !                    auxMatrix=0.0_8

  !                    !! open file for integrals
  !                    open(UNIT=34,FILE=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)//".ints", &
  !                         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

  !                    do
  !                       read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !                       do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

  !                          if (a(u) == -1) goto 30

  !                          coulomb = WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values(a(u),b(u))*integral(u)

  !                          auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

  !                          if( a(u) /= b(u) ) auxMatrix(r(u),s(u)) = auxMatrix(r(u),s(u)) + coulomb

  !                       end do

  !                    end do

  ! 30                 continue

  !                    close(34)

  !                 else

  !                    call IntegralManager_interspecieRepulsionIntegral (currentSpecieID, otherSpecieID, isInterSpecies=.true.)

  !                    auxMatrix=0.0_8

  !                    !! open file for integrals
  !                    open(UNIT=34,FILE=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)//".ints", &
  !                         STATUS='OLD', ACCESS='SEQUENTIAL', FORM='Unformatted')

  !                    do
  !                       read(34)   a(1:CONTROL_instance%INTEGRAL_STACK_SIZE), b(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
  !                            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

  !                       do u = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

  !                          if (a(u) == -1) goto 20

  !                          coulomb = WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values(r(u),s(u))*integral(u)

  !                          auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

  !                          if( r(u) /= s(u) ) auxMatrix(a(u),b(u)) = auxMatrix(a(u),b(u)) + coulomb

  !                       end do

  !                    end do

  ! 20                 continue

  !                    close(34)

  !                 end if

  !              end if

  !              auxMatrix= auxMatrix * MolecularSystem_getCharge( speciesID=currentSpecieID ) &
  !                   * MolecularSystem_getCharge( speciesID=otherSpecieID )

  !              output%values= &
  !                   output%values + auxMatrix



  !           end if

  !        end do

  !        deallocate(auxMatrix)

  !        !! Simetriza la matriz de Acoplamineto
  !        do i = 1 , ssize
  !           do j = i , ssize
  !              output%values(j,i) = &
  !                   output%values(i,j)
  !           end do
  !        end do

  !        WaveFunction_instance( currentSpecieID )%addCouplingMatrix =.true.
  !        WaveFunction_instance( currentSpecieID )%wasBuiltFockMatrix = .false.

  !     end if

  !     !       print *,"Matriz de acoplamiento: ", trim(nameOfSpecieSelected)
  !     !       call Matrix_show( output )

  !   end subroutine WaveFunction_buildCouplingMatrixElectronFree

  !   !! Add nuclear-electron correlation with ADFT (this could become useful)
  !   subroutine WaveFunction_buildInterParticleCorrMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     character(30) :: nameOfOtherSpecie
  !     integer :: numberOfContractions
  !     integer :: otherNumberOfContractions
  !     integer :: currentSpecieID
  !     integer :: otherSpecieID
  !     integer :: speciesIterator
  !     real(8) :: coulomb
  !     real(8) :: startTime, endTime
  !     real(8) :: correlationEnergy

  !     call cpu_time(startTime)

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     currentSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpecieID)

  !     WaveFunction_instance( currentSpecieID )%addInterParticleCorrMatrix =.false.

  !     if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then

  !        otherSpecieID = MolecularSystem_getSpecieID( nameOfSpecie = "e-" ) ! Only for electron and nuclei
  !        nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )
  !        OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

  !        !! Restringe suma de terminos repulsivos de la misma especie.
  !        if ( otherSpecieID /= currentSpecieID ) then
  !           if ( nameofOtherSpecie .ne. nameOfSpecieSelected ) then

  !              if ( .not.allocated(WaveFunction_instance( otherSpecieID)%wavefunction_instance(speciesID)%densityMatrix%values) ) then


  !                 call WaveFunction_exception(ERROR, "Class object WaveFunction_RHF in the builtCouplingMatrix(" &
  !                      // trim(nameOfOtherSpecie) //") function", &
  !                      "Density matrix for "// trim(nameOfOtherSpecie) //" specie, hasn't been defined." )

  !              end if

  !              if ( CONTROL_instance%CALL_DFT ) then

  !                 call bld_aux_ks_c_mat(6, &
  !                      WaveFunction_instance( otherSpecieID )%bridge%system, &
  !                      WaveFunction_instance( otherSpecieID )%bridge%system, &
  !                      WaveFunction_instance( currentSpecieID )%bridge%system, &
  !                      WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( currentSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values, &
  !                      WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      WaveFunction_instance( currentSpecieID )%wavefunction_instance(speciesID)%interParticleCorrMatrix%values, &
  !                      otherNumberOfContractions,NumberOfContractions, &
  !                      correlationEnergy,.false.)

  ! 		WaveFunction_instance(currentSpecieID)%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy=correlationEnergy
  ! 		      WaveFunction_instance( otherSpecieID )%addInterParticleCorrMatrix =.true.
  ! 		      WaveFunction_instance( currentSpecieID )%addInterParticleCorrMatrix =.true.
  ! 		      WaveFunction_instance( currentSpecieID )%wasBuiltFockMatrix = .false.

  !              else 
  !                 WaveFunction_instance(currentSpecieID)%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy= 0.0
  !              end if
  !           end if
  !        end if


  !     end if

  !   end subroutine WaveFunction_buildInterParticleCorrMatrix

  !   !<
  !   !! @brief Contruye una matriz de interaccion con un potencial externo
  !   !!
  !   !! @param nameOfSpecie nombre de la especie seleccionada.
  !   !>
  !   subroutine WaveFunction_buildExternalPotentialMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( nameOfspecie /= "e-BETA" ) then

  ! 	if( WaveFunction_instance(speciesID)%isThereExternalPotential ) then
  ! 		WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix%values = 0.0_8


  ! 	if ( CONTROL_instance%NUMERICAL_INTEGRATION_FOR_EXTERNAL_POTENTIAL )	then	!! Numerical integration
  ! 		if ( trim(ExternalPotential_Manager_instance%externalsPots(1)%name) == "none" ) then
  ! 			WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 				IntegralManager_getNumericalInteractionWithPotentialMatrix( &
  ! 				ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )

  ! 		else 		!! From xml file
  ! 			WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 				IntegralManager_getNumericalPotentialMatrixFromXml( &
  ! 				ExternalPotential_Manager_instance%externalsPots, speciesID, integralName="external" )
  ! 		end if
  ! 	else		!! Analytical Integration	

  ! 		WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  ! 		IntegralManager_getInteractionWithPotentialMatrix( &
  ! 		ExternalPotential_Manager_instance%externalsPots, speciesID, "external" )

  !           end if

  !        end if

  !    else !! Use the same matrix for e-beta and e-alpha

  !        WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix = &
  !             WaveFunction_instance( MolecularSystem_getSpecieID( nameOfSpecie="e-ALPHA" ))%wavefunction_instance(speciesID)%externalPotentialMatrix

  !     end if

  !     			print *,"EXTERNAL POTENTIAL MATRIX FOR: ", nameOfSpecie
  !     			call Matrix_show(WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix)

  !   end subroutine WaveFunction_buildExternalPotentialMatrix





  !   subroutine WaveFunction_buildPuntualParticleMatrix( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if(.not. associated(WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr)) then
  !        call IntegralManager_buildMatrix( ATTRACTION_INTEGRALS, trim(nameOfSpecieSelected  ) )

  !        WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr => 	&
  !             IntegralManager_getMatrixPtr(ATTRACTION_INTEGRALS, nameOfSpecieSelected)

  !     end if

  !     WaveFunction_instance(speciesID)%puntualParticleMatrix%values = &
  !          WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr%values

  !   end subroutine WaveFunction_buildPuntualParticleMatrix





  !   !**
  !   ! @brief Calcula las componentes de energia para la especie especificada
  !   !
  !   ! @warning 	Debe garantizarse el llamdo de esta funcion solo si previamente a llamado a
  !   !			obtainTotalEnergyForSpecie
  !   !**
  !   subroutine WaveFunction_obtainEnergyComponents( nameOfSpecie )
  !     implicit none
  !     character(*), optional :: nameOfSpecie

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID
  !     integer :: otherSpecieID
  !     type(Matrix) :: auxMatrix


  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     if( trim(nameOfSpecieSelected) == "e-ALPHA") then

  !        otherSpecieID = MolecularSystem_getSpecieID(nameOfSpecie="e-BETA")

  !     else if (trim(nameOfSpecieSelected) == "e-BETA") then

  !        otherSpecieID = MolecularSystem_getSpecieID(nameOfSpecie="e-ALPHA")

  !     end if

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     ! 		if( nameOfSpecieSelected == "e-") then
  !     ! 		    auxMatrix=IntegralManager_getInteractionBetweenQuantumClassicMtx(speciesID,3)
  !     ! 		    print *,"ATRACION H-e:",sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values )&
  !     ! 			 * auxMatrix%values )
  !     ! 		end if

  !     !! Calcula energia de repulsion para la especie dada                
  !     !                if ( trim(nameOfSpecieSelected) == "e-ALPHA" .or. trim(nameOfSpecieSelected) == "e-BETA") then
  !     !
  !     !                   WaveFunction_instance( speciesID )%repulsionEnergy= (0.5_8 * sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !     !                        * wavefunction_instance(speciesID)%twoParticlesMatrix%values )) + &
  !     !                        (0.5_8 * sum( transpose( WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%densityMatrix%values ) &
  !     !                        * WaveFunction_instance( otherSpecieID )%wavefunction_instance(speciesID)%twoParticlesMatrix%values ))

  !     !                else 

  !     WaveFunction_instance( speciesID )%repulsionEnergy= 0.5_8 * sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !          * wavefunction_instance(speciesID)%twoParticlesMatrix%values )
  !     !                end if

  !     !! Calcula energia de particula independiente
  !     WaveFunction_instance( speciesID )%independentParticleEnergy = sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values )&
  !          * WaveFunction_instance( speciesID )%hcoreMatrix%values )

  !     !! Calcula energia cinetica para la especie dada
  !     WaveFunction_instance( speciesID )%kineticEnergy = sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values)  &
  !          * WaveFunction_instance( speciesID )%kineticMatrixValuesPtr%values )

  !     !! Calcula energia de potencial externo para la especie dada
  !     if( WaveFunction_instance( speciesID )%isThereExternalPotential ) &
  !          WaveFunction_instance( speciesID )%externalPotentialEnergy = sum( transpose(wavefunction_instance(speciesID)%densityMatrix%values)  &
  !          * WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%externalPotentialMatrix%values )

  !     !! Calcula energia de interaccion entre particulas puntuales y cuanticas
  !     WaveFunction_instance( speciesID )%puntualInteractionEnergy =  WaveFunction_instance( speciesID )%independentParticleEnergy &
  !          - WaveFunction_instance( speciesID )%kineticEnergy

  !     !! Calula enegia de especie independiente (  sin considerar el termino de acoplamiento)
  !     if( .not.WaveFunction_instance( speciesID )%isThereExternalPotential ) then

  !        WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !             sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !             *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
  !             + 0.5_8 * wavefunction_instance(speciesID)%twoParticlesMatrix%values))
  !     else

  !        WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !             sum(  transpose(wavefunction_instance(speciesID)%densityMatrix%values) &
  !             *  (  ( WaveFunction_instance( speciesID )%hcoreMatrix%values ) &
  !             + 0.5_8 * wavefunction_instance(speciesID)%twoParticlesMatrix%values &
  !             + WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%externalPotentialMatrix%values))

  !     end if

  !     WaveFunction_instance( speciesID )%independentSpecieEnergy = &
  !          WaveFunction_instance( speciesID )%independentSpecieEnergy + &
  !          WaveFunction_instance( SpecieID )%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy




  !     !! Calcula energia de acoplamiento en caso de mas de una especie presente
  !     if ( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) &
  !          WaveFunction_instance( speciesID )%couplingEnergy = sum( transpose( wavefunction_instance(speciesID)%densityMatrix%values ) &
  !          * WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%couplingMatrix%values )

  !     !print *, "__________________ ENERGY COMPONENTS _______________________"
  !     !print *, "	Specie                       ", MolecularSystem_getNameOfSpecie( speciesID )
  !     !print *, "	Total Energy                =", WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%totalEnergyForSpecie
  !     !print *, "	Indepent Specie Energy      =", WaveFunction_instance( speciesID )%independentSpecieEnergy
  !     !print *, "	Kinetic Energy              =",WaveFunction_instance( speciesID )%kineticEnergy
  !     !print *, "	Puntual Interaction Energy  =",WaveFunction_instance( speciesID )%puntualInteractionEnergy
  !     !print *, "	Independent Particle Energy =",WaveFunction_instance( speciesID )%independentParticleEnergy
  !     !print *, "	N-E Correlation Energy      =",WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%nuclearElectronicCorrelationEnergy
  !     !print *, "	Repultion Energy            =",WaveFunction_instance( speciesID )%repulsionEnergy
  !     !print *, "	Coupling Energy             =", WaveFunction_instance( speciesID )%couplingEnergy
  !     !print *, "____________________________________________________________"

  !   end subroutine WaveFunction_obtainEnergyComponents

  !   !**
  !   ! @brief indica el objeto ha sido instanciado
  !   !
  !   !**
  !   function WaveFunction_isInstanced() result(output)
  !     implicit none
  !     logical :: output

  !     if ( allocated( WaveFunction_instance ) ) then
  !        output = .true.
  !     else
  !        output = .false.
  !     end if

  !   end function WaveFunction_isInstanced


  !   !**
  !   ! @brief Retorna la matrix de Overlap
  !   !
  !   ! @param nameOfSpecie nombre de la especie seleccionada.
  !   !**
  !   function WaveFunction_getOverlapMatrix( nameOfSpecie, sspeciesID ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     integer, optional :: sspeciesID
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     if ( present(sspeciesID) ) then

  !        speciesID = sspeciesID

  !     else

  !        nameOfSpecieSelected = "e-"

  !        if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !        speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     end if

  !     if ( .not.associated(WaveFunction_instance(speciesID)%overlapMatrixValuesPtr) ) then

  !        call WaveFunction_buildOverlapMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%overlapMatrixValuesPtr )

  !   end function WaveFunction_getOverlapMatrix





  !   !**
  !   ! @brief Retorna la matrix de dos particulas para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getTwoParticlesMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( wavefunction_instance(speciesID)%twoParticlesMatrix%values) ) then

  !        call WaveFunction_buildTwoParticlesMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, wavefunction_instance(speciesID)%twoParticlesMatrix )

  !   end function WaveFunction_getTwoParticlesMatrix

  !   !**
  !   ! @brief Retorna la matrix de acoplamiento para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getCouplingMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     if(trim(nameOfSpecieSelected) == "e-ALPHA" .or. trim(nameOfSpecieSelected) == "e-BETA" ) then
  !        call WaveFunction_buildCouplingMatrixElectronFree( nameOfSpecieSelected, output )
  !        return
  !     end if

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%couplingMatrix%values) ) then

  !        call WaveFunction_buildCouplingMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%couplingMatrix )

  !   end function WaveFunction_getCouplingMatrix

  !   !**
  !   ! @brief Retorna la matrix de interaccion con un potencial externo
  !   !
  !   !**
  !   function WaveFunction_getExternalPotentialMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix%values) ) then

  !        call WaveFunction_buildExternalPotentialMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%externalPotentialMatrix )

  !   end function WaveFunction_getExternalPotentialMatrix


  !   !**
  !   ! @brief Contruye la matrix de Fock para la especie especificada
  !   !
  !   !**
  !   function WaveFunction_getFockMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( .not.allocated( WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%fockMatrix%values) ) then

  !        call WaveFunction_buildFockMatrix( nameOfSpecieSelected )

  !     end if

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%wavefunction_instance(speciesID)%fockMatrix )

  !   end function WaveFunction_getFockMatrix


  !   !**
  !   ! @brief Retorna la matriz  la ultima matriz de densidad calculada para una especie especificada
  !   !
  !   !**
  !   function WaveFunction_getDensityMatrix( nameOfSpecie, sspeciesID ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     integer, optional :: sspeciesID
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     if ( present(sspeciesID) ) then

  !        speciesID = sspeciesID

  !     else

  !        nameOfSpecieSelected = "e-"

  !        if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !        speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     end if

  !     if ( allocated( wavefunction_instance(speciesID)%densityMatrix%values) ) then

  !        call Matrix_copyConstructor( output, wavefunction_instance(speciesID)%densityMatrix )

  !     end if

  !   end function WaveFunction_getDensityMatrix


  !   function WaveFunction_getKineticMatrix( nameOfSpecie ) result(output)
  !     implicit none

  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"

  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     if ( allocated( WaveFunction_instance(speciesID)%kineticMatrix%values) ) then

  !        call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%kineticMatrix)

  !     end if

  !   end function WaveFunction_getKineticMatrix


  !   !**
  !   ! @brief Retorna la matrix de acoplamiento para la especie especificada
  !   !        respecto a las particula puntuales
  !   !
  !   !**
  !   function WaveFunction_getPuntualParticleMatrix( nameOfSpecie ) result( output )
  !     implicit none
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) ::  output

  !     character(30) :: nameOfSpecieSelected
  !     integer :: speciesID

  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

  !     call WaveFunction_buildPuntualParticleMatrix( nameOfSpecieSelected )

  !     call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%puntualParticleMatrix )

  !   end function WaveFunction_getPuntualParticleMatrix




  !**
  ! @brief Retorna la matriz  de coeficientes de combinacion
  !
  !**
  function WaveFunction_getWaveFunctionCoefficients( nameOfSpecie ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    type(Matrix) ::  output

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    if ( allocated( WaveFunction_instance(speciesID)%waveFunctionCoefficients%values) ) then

       call Matrix_copyConstructor( output, WaveFunction_instance(speciesID)%waveFunctionCoefficients )

    end if

  end function WaveFunction_getWaveFunctionCoefficients

  !**
  ! @brief Retorna valores propios del sistema molecular
  !
  !**
  function WaveFunction_getMolecularOrbitalsEnergy( nameOfSpecie ) result( output )
    implicit none
    character(*), optional :: nameOfSpecie
    type(Vector) ::  output

    character(30) :: nameOfSpecieSelected
    integer :: speciesID

    nameOfSpecieSelected = "E-"
    if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

    speciesID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )

    if ( allocated( WaveFunction_instance(speciesID)%molecularOrbitalsEnergy%values) ) then

       call Vector_copyConstructor( output, WaveFunction_instance(speciesID)%molecularOrbitalsEnergy )

    end if

  end function WaveFunction_getMolecularOrbitalsEnergy


  !   function WaveFunction_getValueForOrbitalAt( nameOfSpecie, orbitalNum, coordinate ) result(output)
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     integer :: orbitalNum
  !     real(8) :: coordinate(3)
  !     real(8) :: output

  !     integer :: speciesID
  !     character(30) :: nameOfSpecieSelected
  !     integer :: numberOfContractions
  !     integer :: totalNumberOfContractions
  !     integer :: particleID
  !     integer :: contractionID
  !     integer :: i, j
  !     real(8), allocatable :: auxVal(:)


  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )

  !     output=0.0_8
  !     do i=1,numberOfContractions

  !        particleID = MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%particleID
  !        contractionID=MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(i)%contractionIDInParticle

  !        totalNumberOfContractions = MolecularSystem_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital

  !        if( allocated(auxVal)) deallocate(auxVal)
  !        allocate(auxVal(totalNumberOfContractions))

  !        auxVal = ContractedGaussian_getValueAt(MolecularSystem_getContractionPtr( speciesID,  numberOfContraction=i ), coordinate )

  !        do j = 1, totalNumberOfContractions

  !           output = output + auxVal(j) * wavefunction_instance(speciesID)%waveFunctionCoefficients%values(j,orbitalNum)

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
  !     integer :: speciesID
  !     integer :: numberOfContractions
  !     integer :: j
  !     integer :: i
  !     integer :: numOfGraphs
  !     integer :: auxInitOrbitalNum
  !     integer :: auxLastOrbitalNum

  !     ! 	nameOfSpecieSelected = "e-"
  !     ! 	if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     ! 	speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )

  !     ! 	numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
  !     ! 	fileName=trim(CONTROL_instance%INPUT_FILE)//'orbital.'//trim(String_convertIntegerToString(orbitalNum))//'.'//trim(nameOfSpecieSelected)

  !     ! 	select case( flags )

  !     ! 		case(ORBITAL_ALONE)

  !     ! 			open ( 5,FILE=trim(fileName)//".dat", STATUS='REPLACE',ACTION='WRITE')
  !     ! 			do j=-CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,&
  !     ! 				CONTROL_instance%MAXIMUM_RANGE_OF_GRAPHS,1
  !     ! 			write (5,"(F20.10,F20.10)") j*CONTROL_instance%STEP_OF_GRAPHS, &
  !     ! 			WaveFunction_getValueForOrbitalAt( nameOfSpecieSelected,&
  !     ! 			orbitalNum, [0.0_8,0.0_8,j*CONTROL_instance%STEP_OF_GRAPHS] ) !&
  !     ! !			!!! + (WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values(orbitalNum) * CM_NEG1)
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
  !     ! 				+ (WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values(i) * CM_NEG1)
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




  !   !**
  !   ! @brief Resetea los atributos de clase
  !   !**
  !   subroutine WaveFunction_reset()
  !     implicit none


  !     integer :: speciesIterator
  !     integer :: speciesID

  !     do speciesIterator = MolecularSystem_beginSpecie(), MolecularSystem_endSpecie()

  !        speciesID = MolecularSystem_getSpecieID( iteratorOfSpecie=speciesIterator )
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%totalEnergyForSpecie = 0.0_8
  !        WaveFunction_instance( speciesID )%independentSpecieEnergy =0.0_8
  !        WaveFunction_instance( speciesID )%kineticEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%puntualInteractionEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%independentParticleEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%repulsionEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%couplingEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%externalPotentialEnergy = 0.0_8
  !        WaveFunction_instance( speciesID )%addTwoParticlesMatrix = .false.
  !        WaveFunction_instance( speciesID )%addCouplingMatrix = .false.
  !        WaveFunction_instance( speciesID )%addInterParticleCorrMatrix = .false.
  !        WaveFunction_instance( speciesID )%wasBuiltFockMatrix = .false.

  !        wavefunction_instance(speciesID)%waveFunctionCoefficients%values = 0.0_8
  !        WaveFunction_instance( speciesID )%molecularOrbitalsEnergy%values = 0.0_8
  !        WaveFunction_instance( speciesID )%hcoreMatrix%values = 0.0_8
  !        wavefunction_instance(speciesID)%densityMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%transformationMatrix%values = 0.0_8
  !        wavefunction_instance(speciesID)%twoParticlesMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%couplingMatrix%values = 0.0_8
  !        WaveFunction_instance( speciesID )%wavefunction_instance(speciesID)%fockMatrix%values = 0.0_8

  !        if ( associated( WaveFunction_instance(speciesID)%kineticMatrixValuesPtr ) )  &
  !             WaveFunction_instance(speciesID)%kineticMatrixValuesPtr => null()

  !        if ( associated( WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr )) &
  !             WaveFunction_instance(speciesID)%puntualInteractionMatrixValuesPtr => null()

  !        if ( associated( WaveFunction_instance(speciesID)%overlapMatrixValuesPtr )) &
  !             WaveFunction_instance(speciesID)%overlapMatrixValuesPtr => null()


  !     end do

  !   end subroutine WaveFunction_reset


  !   !**
  !   ! @brief Asigna los coeficientes de la función de onda
  !   !**
  !   subroutine WaveFunction_setWaveFunctionCoefficients(wavefunction_instance(speciesID)%waveFunctionCoefficients, nameOfSpecie)
  !     implicit none
  !     character(*), optional, intent(in) :: nameOfSpecie
  !     character(30) :: nameOfSpecieSelected
  !     type(Matrix), intent(in) :: wavefunction_instance(speciesID)%waveFunctionCoefficients
  !     integer :: speciesID
  !     nameOfSpecieSelected = "e-"
  !     if ( present( nameOfSpecie ) )  nameOfSpecieSelected= trim( nameOfSpecie )

  !     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecieSelected ) )


  !     call Matrix_copyConstructor( wavefunction_instance(speciesID)%waveFunctionCoefficients, wavefunction_instance(speciesID)%waveFunctionCoefficients)

  !   end subroutine WaveFunction_setWaveFunctionCoefficients


  !   !<
  !   !! Write out matrices in MOLPRO'S format
  !   !>
  !   subroutine WaveFunction_writeMatrices(nameOfSpecie)
  !     character(*), optional :: nameOfSpecie
  !     type(Matrix) :: A
  !     type(Matrix) :: P, Pnew
  !     type(Matrix) :: C, Cnew
  !     type(Matrix) :: S, Snew
  !     type(Matrix) :: K, Knew
  !     type(Matrix) :: Pot, PotNew
  !     type(Matrix) :: Jcoup, JcoupNew
  !     type(Matrix) :: Icoup, IcoupNew
  !     type(Exception) :: ex

  !     character(20) :: fileName
  !     character(50) :: nameOfSpecieSelected

  !     integer :: numberOfContractions
  !     integer :: totalNumberOfContractions
  !     integer :: speciesID
  !     integer :: st(1), pt(3), dt(6), ft(10), gt(15)
  !     integer :: angularMoment,numCartesianOrbital
  !     integer, allocatable :: trans(:)
  !     integer :: i, j
  !     integer :: u, v
  !     real(8) :: Epn, Ecn
  !     real(8) :: aVal, bVal

  !     !!Esta rutina no se usa mas... resucitoo....!!!
  !     !return

  !     nameOfSpecieSelected="e-"
  !     if(present(nameOfSpecie)) nameOfSpecieSelected=trim(nameOfSpecie)

  !     fileName = CONTROL_instance%INPUT_FILE

  !     !!Id de la especie seleccionada (por defecto e-)
  !     speciesID = MolecularSystem_getSpecieID(nameOfSpecie = trim(nameOfSpecieSelected))

  !     !!Numero total de contracciones
  !     numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)
  !     totalNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(speciesID)

  !     !!Tamano del arreglo de nuevas etiquetas
  !     if(allocated(trans)) deallocate(trans)
  !     allocate(trans(totalNumberOfContractions))

  !     !! Reglas de transformacion de indices
  !     st(1)    = 1
  !     pt(1:3)  = [1, 2, 3]
  !     dt(1:6)  = [1, 4, 5, 2, 6, 3]
  !     ft(1:10) = [1, 4, 5, 6, 10, 8, 2, 7, 9, 3]
  !     gt(1:15) = [1, 4, 5, 10, 13, 11, 6, 14, 15, 8, 2, 7, 12, 9, 3]

  !     write(*,*) "RULES FOR INDEX TRANSFORMATION"
  !     write(*,*) "================================"
  !     write(*,*) ""
  !     write(*,"(' s =',I3)") st(1)
  !     write(*,"(' p =',3I3)") ( pt(i), i=1,3 )
  !     write(*,"(' d =',6I3)") ( dt(i), i=1,6 )
  !     write(*,"(' f =',10I3)") ( ft(i), i=1,10 )
  !     write(*,"(' g =',15I3)") ( gt(i), i=1,15 )
  !     write(*,*) ""

  !     trans = 0
  !     u = 1
  !     v = 0

  !     do i = 1, numberOfContractions
  !        angularMoment = MolecularSystem_instance%particlesPtr(MolecularSystem_instance%idsOfContractionsForSpecie(&
  !             speciesID)%contractionID(i)%particleID)%basis%contractions( &
  !             MolecularSystem_instance%idsOfContractionsForSpecie(speciesID)%contractionID(&
  !             i)%contractionIDInParticle)%angularMoment

  !        select case(CONTROL_instance%DIMENSIONALITY)
  !        case(3)
  !           numCartesianOrbital = ((angularMoment + 1_8)*(angularMoment + 2_8))/2_8
  !        case(2)
  !           numCartesianOrbital = ((angularMoment + 1_8))
  !        case(1)
  !           numCartesianOrbital = 1 
  !        case default
  !           call Exception_constructor( ex , WARNING )
  !           call Exception_setDebugDescription( ex, "Class object WaveFunction in the writeMatrices function" )
  !           call Exception_setDescription( ex, "this Dimensionality is not implemented (D>3)" )
  !           call Exception_show( ex )

  !           numCartesianOrbital = 0
  !        end select

  !        select case(angularMoment)
  !        case(0)
  !           trans(u) = st(1) + v
  !           u = u + 1

  !        case(1)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = pt(j) + v
  !              u = u + 1
  !           end do

  !        case(2)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = dt(j) + v
  !              u = u + 1
  !           end do
  !        case(3)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = ft(j) + v
  !              u = u + 1
  !           end do
  !        case(4)
  !           do j = 1, numCartesianOrbital
  !              trans(u) = gt(j) + v
  !              u = u + 1
  !           end do

  !        case default
  !           call Exception_constructor( ex , WARNING )
  !           call Exception_setDebugDescription( ex, "Class object WaveFunction in the writeMatrices function" )
  !           call Exception_setDescription( ex, "this angular moment is not implemented (l>4)" )
  !           call Exception_show( ex )

  !           return

  !        end select
  !        v = u - 1
  !     end do
  !     i = 1
  !     !write(*,"(' trans =',<totalNumberOfContractions>I3)") ( trans(i), i=1,totalNumberOfContractions )

  !     !Epn = MolecularSystem_instance%repulsionEnergy(2)
  !     !Ecn = MolecularSystem_instance%kineticEnergy(2)

  !     !write(*,*) "Epn = ", Epn
  !     !write(*,*) "Ecn = ", Ecn

  !     P = WaveFunction_getDensityMatrix( trim(nameOfSpecieSelected) )
  !     C = WaveFunction_getWaveFunctionCoefficients( trim(nameOfSpecieSelected) )
  !     S = WaveFunction_getOverlapMatrix( trim(nameOfSpecieSelected) )
  !     K = WaveFunction_getKineticMatrix( trim(nameOfSpecieSelected) )
  !     Pot= WaveFunction_getExternalPotentialMatrix(trim(nameOfSpecieSelected)) 
  !     Jcoup = WaveFunction_getCouplingMatrix( trim(nameOfSpecieSelected) )
  !     Icoup = WaveFunction_getPuntualParticleMatrix( trim(nameOfSpecieSelected) )

  !     call Matrix_constructor( Pnew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Cnew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Snew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( Knew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( JcoupNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )
  !     call Matrix_constructor( IcoupNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )    
  !     call Matrix_constructor( PotNew, int(totalNumberOfContractions, 8), int(totalNumberOfContractions, 8) )

  !     do i=1,totalNumberOfContractions
  !        do j=1,totalNumberOfContractions
  !           Pnew%values( trans(i), trans(j) ) = P%values( i, j )
  !           Cnew%values( trans(i), trans(j) ) = C%values( i, j )
  !           Snew%values( trans(i), trans(j) ) = S%values( i, j )
  !           Knew%values( trans(i), trans(j) ) = K%values( i, j )
  !           JcoupNew%values( trans(i), trans(j) ) = Jcoup%values( i, j )
  !           IcoupNew%values( trans(i), trans(j) ) = Icoup%values( i, j )
  !           PotNew%values( trans(i), trans(j) ) = Pot%values( i, j )
  !        end do
  !     end do

  !     ! 		do i=1,totalNumberOfContractions
  !     ! 			aVal = Pnew.values( 3, i )
  !     ! 			bVal = Pnew.values( 6, i )
  !     ! 			Pnew.values( 6, i ) = aVal
  !     ! 			Pnew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = Cnew.values( 3, i )
  !     ! 			bVal = Cnew.values( 6, i )
  !     ! 			Cnew.values( 6, i ) = aVal
  !     ! 			Cnew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = Snew.values( 3, i )
  !     ! 			bVal = Snew.values( 6, i )
  !     ! 			Snew.values( 6, i ) = aVal
  !     ! 			Snew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = JcoupNew.values( 3, i )
  !     ! 			bVal = JcoupNew.values( 6, i )
  !     ! 			JcoupNew.values( 6, i ) = aVal
  !     ! 			JcoupNew.values( 3, i ) = bVal
  !     !
  !     ! 			aVal = IcoupNew.values( 3, i )
  !     ! 			bVal = IcoupNew.values( 6, i )
  !     ! 			IcoupNew.values( 6, i ) = aVal
  !     ! 			IcoupNew.values( 3, i ) = bVal
  !     ! 		end do


  !     write(*,"(3A)", advance="no") " Saving coefficients matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"coeff"), " ) ... "
  !     open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"coeff")), action='write', form='unformatted' )
  !     write(20) int(size(Cnew%values), 8)
  !     write(20) Cnew%values
  !     close(20)
  !     write(*,*) "OK"

  !     if(trim(nameOfSpecieSelected) /= "e-BETA") then 

  !        write(*,"(3A)", advance="no") " Saving potential matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"pot"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"pot")), action='write', form='unformatted')
  !        write(20) int(size(PotNew%values), 8)
  !        write(20) PotNew%values
  !        close(20)
  !        write(*,*) "OK"


  !        write(*,"(3A)", advance="no") " Saving density matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"dens"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"dens")), action='write', form='unformatted' )
  !        write(20) int(size(Pnew%values), 8)
  !        write(20) Pnew%values
  !        close(20)
  !        write(*,*) "OK"

  !        !		write(*,*) "Jcoup ="
  !        !		call Matrix_show( JcoupNew )

  !        !		write(*,*) "Icoup ="
  !        !		call Matrix_show( IcoupNew )

  !        write(*,"(3A)", advance="no") " Saving overlap matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"over"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"over")), action='write', form='unformatted' )
  !        write(20) int(size(Snew%values), 8)
  !        write(20) Snew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving kinetic matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"kin"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"kin")), action='write', form='unformatted' )
  !        write(20) int(size(Snew%values), 8)
  !        write(20) Knew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving coupling matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"jcoup"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"jcoup")), action='write', form='unformatted' )
  !        write(20) int(size(JcoupNew%values), 8)
  !        write(20) JcoupNew%values
  !        close(20)
  !        write(*,*) "OK"

  !        write(*,"(3A)", advance="no") " Saving fixed interaction matrix ( ", trim(trim(fileName)//trim(nameOfSpecieSelected)//"."//"icoup"), " ) ... "
  !        open( 20, file=trim(String_getLowercase(trim(fileName)//trim(nameOfSpecieSelected)//"."//"icoup")), action='write', form='unformatted' )
  !        write(20) int(size(IcoupNew%values), 8)
  !        write(20) IcoupNew%values
  !        close(20)
  !        write(*,*) "OK"

  !     end if

  !   end subroutine WaveFunction_writeMatrices

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

  !! build cosmo 2 matrix 

  subroutine WaveFunction_buildCosmo2Matrix(nameOfSpecie)

    character(*), optional :: nameOfSpecie
    type(species) :: specieSelected
    character(30) :: nameOfSpecieSelected
    integer :: speciesID


    integer, allocatable :: labels(:)
    ! real(8), allocatable :: cosmo_int(:)
    real(8), allocatable :: ints_mat_aux(:,:)
    real(8), allocatable :: cosmo2_aux(:,:)

    real(8), allocatable :: qe(:)
    real(8) :: cosmo_int


    integer :: g,i,ii,h,hh,j,jj,k,l,m,o,p
    integer :: iii,jjj,hhh,gg,ll,pp,oo

    integer:: auxLabelsOfContractions
    integer:: a, b, c


    nameOfSpecieSelected = "E-"
    if (present(nameOfSpecie))  nameOfSpecieSelected= trim(nameOfSpecie)
    speciesID = MolecularSystem_getSpecieID(nameOfSpecie=trim(nameOfSpecieSelected))
    specieSelected=MolecularSystem_instance%species(speciesID)


    open(unit=110, file=trim(nameOfSpecieSelected)//"_qq.inn", status='old', form="unformatted")
    read(110)m


    ! if(allocated(cosmo_int)) deallocate(cosmo_int)
    ! allocate(cosmo_int(m))

    ! close(unit=110)


    if(allocated(labels)) deallocate(labels)
    allocate(labels(MolecularSystem_instance%species(speciesID)%basisSetSize))

    if(allocated(ints_mat_aux)) deallocate(ints_mat_aux)
    allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID)))


    if(allocated(cosmo2_aux)) deallocate(cosmo2_aux)
    allocate(cosmo2_aux(MolecularSystem_getTotalNumberOfContractions(speciesID), MolecularSystem_getTotalNumberOfContractions(speciesID)))


    auxLabelsOfContractions = 1

    c = 0
    do a = 1, size(specieSelected%particles)
       do b = 1, size(specieSelected%particles(a)%basis%contraction)

          c = c + 1

          !!position for cartesian contractions

          labels(c) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital


       end do
    end do


    ! call Matrix_show(wavefunction_instance(speciesID)%densityMatrix)

    m = 0

    ii = 0
    do g = 1, size(MolecularSystem_instance%species(speciesID)%particles)
       do h = 1, size(MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction)

          hh = h
          ii = ii + 1
          jj = ii - 1

          do i = g, size(MolecularSystem_instance%species(speciesID)%particles)
             do j = hh, size(MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction)

                jj = jj + 1

                !!saving integrals on Matrix
                do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(speciesID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(speciesID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      iii=0
                      do gg = 1, size(MolecularSystem_instance%species(speciesID)%particles)
                         do ll = 1, size(MolecularSystem_instance%species(speciesID)%particles(gg)%basis%contraction)

                            hhh = ll
                            iii = iii + 1
                            jjj = iii - 1
                            do p = gg, size(MolecularSystem_instance%species(speciesID)%particles)
                               do o = hhh, size(MolecularSystem_instance%species(speciesID)%particles(p)%basis%contraction)
                                  jjj = jjj + 1

                                  !!saving integrals on Matrix
                                  do pp = labels(iii), labels(iii) + (MolecularSystem_instance%species(speciesID)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                     do oo = labels(jjj), labels(jjj) + (MolecularSystem_instance%species(speciesID)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
                                        m = m + 1

                                        read(110)cosmo_int

                                        ints_mat_aux(pp, oo) =(wavefunction_instance(speciesID)%densityMatrix%values(pp,oo))* cosmo_int
                                        ints_mat_aux(oo, pp) = ints_mat_aux(pp, oo)

                                     end do
                                  end do
                               end do
                               hhh = 1
                            end do

                         end do
                      end do
                      cosmo2_aux(k,l)=0.0_8
                      do pp=1,size(ints_mat_aux,DIM=1)
                         do oo=1,size(ints_mat_aux,DIM=1)
                            cosmo2_aux(k,l)=cosmo2_aux(k,l)+ints_mat_aux(pp,oo)
                            wavefunction_instance(speciesID)%cosmo2%values(k,l)=cosmo2_aux(k,l)
                            wavefunction_instance(speciesID)%cosmo2%values(l,k)=wavefunction_instance(speciesID)%cosmo2%values(k,l)
                         end do
                      end do
                   end do
                end do
             end do
             hh = 1
          end do
       end do
    end do

    close(unit=110)

    if (  CONTROL_instance%DEBUG_SCFS) then
       write(*,*) "COSMO 2 matrix for: ", trim(nameOfSpecieSelected)
       call Matrix_show(wavefunction_instance(speciesID)%cosmo2)
    end if



  end subroutine WaveFunction_buildCosmo2Matrix


  subroutine WaveFunction_buildCosmoCoupling(nameOfSpecie)

    character(*), optional :: nameOfSpecie
    type(species) :: specieSelected
    type(species) :: otherSpecieSelected
    character(30) :: nameOfSpecieSelected
    character(30) :: nameOfOtherSpecie

    integer, allocatable :: labels(:)
    integer, allocatable :: otherLabels(:)
    ! real(8), allocatable :: cosmo_int(:)
    real(8), allocatable :: ints_mat_aux(:,:)
    real(8), allocatable :: cosmoCoup_aux(:,:)

    real(8), allocatable :: auxMatrix(:,:)

    real(8), allocatable :: qe(:)
    real(8) :: cosmo_int


    integer :: currentSpecieID
    integer :: otherSpecieID
    integer :: numberOfContractions
    integer :: otherNumberOfContractions
    integer :: speciesIterator

    integer :: g,i,ii,h,hh,j,jj,k,l,m,o,p
    integer :: iii,jjj,hhh,gg,ll,pp,oo

    integer:: auxLabelsOfContractions
    integer:: otherAuxLabelsOfContractions
    integer:: a, b, c


    nameOfSpecieSelected = "E-"
    if (present(nameOfSpecie))  nameOfSpecieSelected= trim(nameOfSpecie)

    currentSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecieSelected )
    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(currentSpecieID)
    specieSelected=MolecularSystem_instance%species(currentSpecieID)

    if(allocated(labels)) deallocate(labels)
    allocate(labels(MolecularSystem_instance%species(currentSpecieID)%basisSetSize))

    wavefunction_instance(currentSpecieID)%cosmoCoupling%values(:,:)=0.0_8

    auxLabelsOfContractions = 1



    c = 0
    do a = 1, size(specieSelected%particles)
       do b = 1, size(specieSelected%particles(a)%basis%contraction)

          c = c + 1

          !!position for cartesian contractions

          labels(c) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital


       end do
    end do


    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then
             
			wavefunction_instance(currentSpecieID)%cosmoCoupling%values = 0.0_8


       do speciesIterator = 1, MolecularSystem_getNumberOfQuantumSpecies()

          otherSpecieID = speciesIterator
          nameOfOtherSpecie = MolecularSystem_getNameOfSpecie( otherSpecieID )          
          OtherNumberOfContractions = MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          otherSpecieSelected=MolecularSystem_instance%species(otherSpecieID)

          if ( otherSpecieID /= currentSpecieID ) then

             ! write(*,*)"hola other and current", otherSpecieID,currentSpecieID 

      !      wavefunction_instance(currentSpecieID)%cosmoCoupling%values = 0.0_8

             open(unit=110, file=trim(nameOfOtherSpecie)//trim(nameOfSpecieSelected)//"_qq.cup", status='old', form="unformatted")
             ! open(unit=110, file=trim(nameOfSpecieSelected)//trim(nameOfOtherSpecie)//"_qq.cup", status='old', form="unformatted")
             read(110)m

             ! if(allocated(cosmo_int)) deallocate(cosmo_int)
             ! allocate(cosmo_int(m))


             if(allocated(otherLabels)) deallocate(otherLabels)
             allocate(otherLabels(MolecularSystem_instance%species(otherSpecieID)%basisSetSize))

             otherAuxLabelsOfContractions=1

             c = 0

             do a = 1, size(otherSpecieSelected%particles)
                do b = 1, size(otherSpecieSelected%particles(a)%basis%contraction)

                   c = c + 1

                   !!position for cartesian contractions

                   otherlabels(c) = otherAuxLabelsOfContractions
                   otherAuxLabelsOfContractions = otherAuxLabelsOfContractions + otherSpecieSelected%particles(a)%basis%contraction(b)%numCartesianOrbital

                end do
             end do


             if(allocated(ints_mat_aux)) deallocate(ints_mat_aux)
             allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(otherSpecieID), MolecularSystem_getTotalNumberOfContractions(otherSpecieID)))

             ints_mat_aux=0.0_8                


             if(allocated(cosmoCoup_aux)) deallocate(cosmoCoup_aux)
             allocate(cosmoCoup_aux(MolecularSystem_getTotalNumberOfContractions(currentSpecieID), MolecularSystem_getTotalNumberOfContractions(currentSpecieID)))


             m = 0

             ii = 0
             do g = 1, size(MolecularSystem_instance%species(currentSpecieID)%particles)
                do h = 1, size(MolecularSystem_instance%species(currentSpecieID)%particles(g)%basis%contraction)

                   hh = h
                   ii = ii + 1
                   jj = ii - 1

                   do i = g, size(MolecularSystem_instance%species(currentSpecieID)%particles)
                      do j = hh, size(MolecularSystem_instance%species(currentSpecieID)%particles(i)%basis%contraction)

                         jj = jj + 1

                         !!saving integrals on Matrix
                         do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(currentSpecieID)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                            do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(currentSpecieID)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                               iii=0
                               do gg = 1, size(MolecularSystem_instance%species(otherSpecieID)%particles)
                                  do ll = 1, size(MolecularSystem_instance%species(otherSpecieID)%particles(gg)%basis%contraction)

                                     hhh = ll
                                     iii = iii + 1
                                     jjj = iii - 1

                                     do p = gg, size(MolecularSystem_instance%species(otherSpecieID)%particles)
                                        do o = hhh, size(MolecularSystem_instance%species(otherSpecieID)%particles(p)%basis%contraction)
                                           jjj = jjj + 1

                                           !!saving integrals on Matrix
                                           do pp = otherlabels(iii), otherlabels(iii) + (MolecularSystem_instance%species(otherSpecieID)%particles(gg)%basis%contraction(ll)%numCartesianOrbital - 1)
                                              do oo = otherlabels(jjj), otherlabels(jjj) + (MolecularSystem_instance%species(otherSpecieID)%particles(p)%basis%contraction(o)%numCartesianOrbital - 1)
                                                 m = m + 1

                                                 ! write(*,*)"m,cosmo_int(m),P_element,pp,oo",m,cosmo_int(m),wavefunction_instance(otherSpecieID)%densityMatrix%values(pp,oo),pp,oo
                                                 read(110)cosmo_int
                                                 ints_mat_aux(pp, oo) =(wavefunction_instance(otherSpecieID)%densityMatrix%values(pp,oo))* cosmo_int
                                                 ints_mat_aux(oo, pp) = ints_mat_aux(pp, oo)

                                              end do
                                           end do

                                        end do
                                        hhh = 1
                                     end do

                                  end do
                               end do
															 ! write(*,*)"m ", m
                               cosmoCoup_aux(k,l)=0.0_8
                               do pp=1,size(ints_mat_aux,DIM=1)
                                  do oo=1,size(ints_mat_aux,DIM=1)
                                     cosmoCoup_aux(k,l)=cosmoCoup_aux(k,l)+ints_mat_aux(pp,oo)
                                  end do
                               end do
                            end do
                         end do
                      end do
                      hh = 1
                   end do
                end do
             end do
                               do k=1,size(cosmoCoup_aux,DIM=1)
                                  do l=k,size(cosmoCoup_aux,DIM=1)
																		wavefunction_instance(currentSpecieID)%cosmoCoupling%values(k,l)=cosmoCoup_aux(k,l)+wavefunction_instance(currentSpecieID)%cosmoCoupling%values(k,l)
																		wavefunction_instance(currentSpecieID)%cosmoCoupling%values(l,k)=wavefunction_instance(currentSpecieID)%cosmoCoupling%values(k,l)
                                  end do
                               end do



             !! debug

             if (  CONTROL_instance%DEBUG_SCFS) then

                write(*,*)"cosmo Coupling = "//trim(nameofSpecieSelected)

                call Matrix_show(wavefunction_instance(currentSpecieID)%cosmoCoupling)
								
								write(*,*)"cosmo density matrix used = "//trim(nameOfOtherSpecie)

                call Matrix_show(wavefunction_instance(otherSpecieID)%densityMatrix)

             end if
             close(unit=110)
          end if
       end do


    end if



  end subroutine WaveFunction_buildCosmoCoupling



end module WaveFunction_
