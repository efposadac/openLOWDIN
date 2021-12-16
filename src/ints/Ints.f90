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

!> @brief This program calculates one-particle integrals for contracted gaussian functions representations
!! this integrals are necessary for any Hartree-Fock scheme that use gaussian functions as initial wave function.
!! @author E. F. Posada (efposadac@unal.edu.co)
!! @version 1.0
!! <b> Fecha de creacion : </b> 2013-02-28
!!
!! <b> Historial de modificaciones: </b>
!!   - <tt> 2013-02-28 </tt>: E. F. Posada ( efposadac@unal.edu.co )
!!        -# Creacion del programa, depuracion y pruebas exahustivas
Program Ints
  use CONTROL_
  use MolecularSystem_
  use EnergyGradients_
  use IntegralManager_
  use String_
  use Stopwatch_
  use CosmoCore_
  use ParticleManager_
  use Libint2Interface_
  use G12Integrals_

  implicit none

  character(50) :: job
  integer :: speciesID, i, j

  type(surfaceSegment) :: surface_aux

  !!Cosmo test
  ! real(8) :: x,y,z
  ! integer :: j
  !Cosmo test


  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.bas format
  call MolecularSystem_loadFromFile( "LOWDIN.BAS" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  select case(trim(job))

  case("ONE_PARTICLE")

     if(CONTROL_instance%LAST_STEP) then
        write(*,"(A)")"----------------------------------------------------------------------"
        write(*,"(A)")"** PROGRAM INTS                          Author: E. F. Posada, 2013   "
        write(*,"(A)")"----------------------------------------------------------------------"
        write(*,"(A)") "INFO: RUNNING IN "//trim(job)//" MODE."
        write(*,"(A)")" "
     end if

     !!Open file to store integrals
     open(unit=30, file="lowdin.opints", status="unknown", form="unformatted")

     !!write global info on output
     write(30) size(MolecularSystem_instance%species)


     ! !!Calculate overlap integrals
     ! call Libint2Interface_compute1BodyInts(1)
     call IntegralManager_getOverlapIntegrals()
     !call IntegralManager_getThreeCenterIntegrals()

     ! !!Calculate kinetic integrals
     ! call Libint2Interface_compute1BodyInts(2)
     call IntegralManager_getKineticIntegrals()

     if ( CONTROL_instance%REMOVE_TRANSLATIONAL_CONTAMINATION ) then
       call IntegralManager_getFirstDerivativeIntegrals()
     end if

     if ( CONTROL_instance%HARMONIC_CONSTANT /= 0.0_8 ) then 
       call IntegralManager_getHarmonicIntegrals()
     end if

     ! !!Calculate attraction integrals
     ! call Libint2Interface_compute1BodyInts(3)
     call IntegralManager_getAttractionIntegrals()

     ! !!Calculate moment integrals
     call IntegralManager_getMomentIntegrals()
      
     !! Calculate integrals with external potential
     if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
       call IntegralManager_getThreeCenterIntegrals()
       !call IntegralManager_getThreeCenterIntegralsByProduct()
     end if
     !stop time
     call Stopwatch_stop(lowdin_stopwatch)

     if(CONTROL_instance%LAST_STEP) then     
        write(*, *) ""
        write(*,"(A,F10.3,A4)") "** TOTAL CPU Time INTS : ", lowdin_stopwatch%enlapsetTime ," (s)"
        write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time INTS : ", lowdin_stopwatch%elapsetWTime ," (s)"
        write(*, *) ""
     end if
     close(30)

  case("COSMO")

     call CosmoCore_lines(surface_aux)
     call CosmoCore_filler(surface_aux)

     !!Open file to store integrals
     open(unit=40, file="cosmo.opints", status="unknown", form="unformatted")

     !!write global info on output
     write(40) size(MolecularSystem_instance%species)

     !!Calculate cosmo integrals and charges
     call IntegralManager_getAttractionIntegrals(surface_aux)

     !stop time
     call Stopwatch_stop(lowdin_stopwatch)

     if(CONTROL_instance%LAST_STEP) then
        write(*, *) ""
        write(*,"(A,F10.3,A4)") "** TOTAL CPU Time Cosmo-INTS : ", lowdin_stopwatch%enlapsetTime ," (s)"
        write(*,"(A,F10.3,A4)") "** TOTAL Elapsed Time Cosmo-INTS : ", lowdin_stopwatch%elapsetWTime ," (s)"
        write(*, *) ""
     end if
     close(40)

  case("TWO_PARTICLE_R12")

     if(CONTROL_instance%LAST_STEP) then
        write(*, "(A)") " "
        write(*, "(A)") " TWO-BODY INTEGRAL SETUP: "
        write(*, "(A)") "------------------------- "
        write(*, "(A, A6)") " Storage: ", trim(String_getUppercase( CONTROL_instance%INTEGRAL_STORAGE ))
        write(*, '(A, A6)') " Scheme: ", trim(String_getUppercase(trim(CONTROL_instance%INTEGRAL_SCHEME)))
        write(*, '(A, I6)') " Stack size: ", CONTROL_instance%INTEGRAL_STACK_SIZE
        write(*, "(A)") " "

        select case (trim(String_getUppercase(trim(CONTROL_instance%INTEGRAL_SCHEME))))

        case("RYS")
           write(*, "(A)") " RYS QUADRATURE SCHEME                 " 
           write(*, "(A)") " LOWDIN-RYS Implementation V. 1.0   Guerrero R. D. ; Posada E. F. 2013 "
           write(*, "(A)") " ----------------------------------------------------------------------"

        case("LIBINT")
           write(*, "(A)") " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
           write(*, "(A)") " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2016   "
           write(*, "(A)") " ----------------------------------------------------------------------"

        case("CUDINT")
           write(*, "(A)") " CUDA ERI Integrals Calculations has been implemented based on:         " 
           write(*, "(A)") " Ufimtsev, I. S.; Martinez, T. J.; JCTC 2008, 4, 222           " 
           write(*, "(A)") " LOWDIN-CUDINT Implementation V. 1.0:  "
           write(*, "(A)") " Rodas, J. M.; Hernandez, R.; Zapata, A.; Galindo, J. F.; Reyes A. 2014   "
           write(*, "(A)") " ----------------------------------------------------------------------"

        case default
           write(*, "(A)") " LIBINT library, Fermann, J. T.; Valeev, F. L. 2010                   " 
           write(*, "(A)") " LOWDIN-LIBINT Implementation V. 2.1  Posada E. F. ; Reyes A. 2016   "
           write(*, "(A)") " ----------------------------------------------------------------------"

        end select
     end if


     !! intra-species two-boy integration
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
        !!Calculate attraction integrals (intra-species)
        call IntegralManager_getIntraRepulsionIntegrals(trim(MolecularSystem_getNameOfSpecie(speciesID)), &
             trim(CONTROL_instance%INTEGRAL_SCHEME))
     end do

     !stop time
     if(CONTROL_instance%LAST_STEP) then
        write(*,"(/A)", advance="no") "*** TOTAL CPU time intra-species integrals : "
        call Stopwatch_splitTime()
        write(*,"(A4)") " (s)"

        write(*,"(A)", advance="no") "*** TOTAL elapsed time intra-species integrals : "
        call Stopwatch_splitWTime()
        write(*,"(A4/)") " (s)"

     end if

     !! inter-species two-boy integration
     if(Molecularsystem_instance%numberOfQuantumSpecies > 1) then

        !!Calculate attraction integrals (inter-species)
        call IntegralManager_getInterRepulsionIntegrals(trim(CONTROL_instance%INTEGRAL_SCHEME))

        !stop time
        call Stopwatch_stop(lowdin_stopwatch)

        if(CONTROL_instance%LAST_STEP) then
           write(*,"(/A,F10.3,A4)") "*** TOTAL CPU time inter-species integrals : ", lowdin_stopwatch%enlapsetTime ," (s)"
           write(*,"(A,F10.3,A4/)") "*** TOTAL elapsed time inter-species integrals : ", lowdin_stopwatch%elapsetWTime ," (s)"
        end if
     end if

  case("GET_GRADIENTS")
     call EnergyGradients_constructor()
     call EnergyGradients_getAnalyticDerivative()

  case("TWO_PARTICLE_G12")

     !! intra-species G12 integration
     do speciesID = 1, MolecularSystem_instance%numberOfQuantumSpecies
        !!Calculate repulsion integrals (intra-species)

        ! call G12Integrals_diskIntraSpecie(speciesID)

        call Libint2Interface_computeG12Intraspecies_disk(speciesID)


     end do

     !stop time
     if(CONTROL_instance%LAST_STEP) then
        write(*,"(/A)", advance="no") "*** TOTAL CPU time  G12 intra-species integrals : "
        call Stopwatch_splitTime()
        write(*,"(A4/)") " (s)"

        write(*,"(/A)", advance="no") "*** TOTAL elapsed time  G12 intra-species integrals : "
        call Stopwatch_splitWTime()
        write(*,"(A4/)") " (s)"

     end if

     !! inter-species two-boy integration
     if(Molecularsystem_instance%numberOfQuantumSpecies > 1) then

        !!Calculate attraction integrals (inter-species)
        do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
          do j = i+1, MolecularSystem_instance%numberOfQuantumSpecies

             call Libint2Interface_computeG12Interspecies_disk(i, j)
        
             ! call G12Integrals_G12diskInterSpecie(trim(MolecularSystem_getNameOfSpecie(i)), &
             !  trim(MolecularSystem_getNameOfSpecie(j)), i, j)

          end do
        end do

        !stop time
        call Stopwatch_stop(lowdin_stopwatch)

        if(CONTROL_instance%LAST_STEP) then
           write(*,"(/A,F10.3,A4/)") "*** TOTAL CPU time G12 inter-species integrals : ", lowdin_stopwatch%enlapsetTime ," (s)"
           write(*,"(/A,F10.3,A4/)") "*** TOTAL elapsed time G12 inter-species integrals : ", lowdin_stopwatch%elapsetWTime ," (s)"
        end if
     end if

  case default

     write(*,*) "USAGE: lowdin-ints.x job "
     write(*,*) "Where job can be: "
     write(*,*) "  ONE_PARTICLE"
     write(*,*) "  TWO_PARTICLE_R12"
     write(*,*) "  GET_GRADIENTS"
     write(*,*) "  TWO_PARTICLE_G12"
     stop "ERROR"

  end select

end Program Ints
