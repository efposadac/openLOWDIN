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
  use MolecularSystem_
  use IntegralManager_
  use String_
  use Stopwatch_
  implicit none
  
  character(50) :: job
  character(50) :: nprocs
  character(50) :: proc
  character(50) :: speciesName

  integer(8) :: nprocess
  integer(8) :: process

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
  
  select case(trim(job))

  case("ONE_PARTICLE")
     
     write(*,"(A)")"----------------------------------------------------------------------"
     write(*,"(A)")"** PROGRAM INTS                          Author: E. F. Posada, 2013   "
     write(*,"(A)")"----------------------------------------------------------------------"
     write(*,"(A)") "INFO: RUNNING IN "//trim(job)//" MODE."
     write(*,"(A)")" "

     !!Open file to store integrals
     open(unit=30, file="lowdin.opints", status="unknown", form="unformatted")
     
     !!write global info on output
     write(30) size(MolecularSystem_instance%species)
     
     !!Calculate overlap integrals
     call IntegralManager_getOverlapIntegrals()

     !!Calculate kinetic integrals
     call IntegralManager_getKineticIntegrals()
     
     !!Calculate attraction integrals
     call IntegralManager_getAttractionIntegrals()
     
     !!Calculate moment integrals
     call IntegralManager_getMomentIntegrals()
     
     !stop time
     call Stopwatch_stop(lowdin_stopwatch)
     
     write(*, *) ""
     write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time INTS : ", lowdin_stopwatch%enlapsetTime ," (s)"
     write(*, *) ""
     close(30)
     
  case("TWO_PARTICLE_R12_INTRA")
     
     nprocs = ""
     proc = ""
     speciesName=""
     
     call get_command_argument(2,value=nprocs)
     call get_command_argument(3,value=proc)
     call get_command_argument(4,value=speciesName)
     
     nprocess = 1
     process = 1

     if(trim(nprocs)/= "")read(nprocs,*) nprocess
     if(trim(proc)/= "")read(proc,*) process     
     
     if(trim(speciesName) == "") speciesName="E-"
     
     !!Calculate attraction integrals (intra-species)
     call IntegralManager_getIntraRepulsionIntegrals(nprocess, process, speciesName, trim(CONTROL_instance%INTEGRAL_SCHEME))
     
     !stop time
     call Stopwatch_stop(lowdin_stopwatch)
     
      write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time INTS : ", lowdin_stopwatch%enlapsetTime ," (s)"

   case("TWO_PARTICLE_R12_INTER")
      
      if(Molecularsystem_instance%numberOfQuantumSpecies > 1) then
         
         !!Calculate attraction integrals (inter-species)
         call IntegralManager_getInterRepulsionIntegrals(trim(CONTROL_instance%INTEGRAL_SCHEME))
         
         !stop time
         call Stopwatch_stop(lowdin_stopwatch)
     
         write(*, *) ""
         write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Time INTS : ", lowdin_stopwatch%enlapsetTime ," (s)"
         write(*, *) ""
         
     end if
     
  case("TWO_PARTICLE_F12")
     
     stop "NOT IMPLEMENTED"
     
  case default
     
     write(*,*) "USAGE: lowdin-ints.x job "
     write(*,*) "Where job can be: "
     write(*,*) "  ONE_PARTICLE"
     write(*,*) "  TWO_PARTICLE_R12_INTRA"
     write(*,*) "  TWO_PARTICLE_R12_INTER"
     write(*,*) "  TWO_PARTICLE_F12"
     stop "ERROR"
     
  end select
  
end Program Ints
