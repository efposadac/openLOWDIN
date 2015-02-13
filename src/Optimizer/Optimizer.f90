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
!! @brief Molecular Mechanics program.
!!        This module minimize the energy
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-09-30
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-09-30 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions using Conjugate Gradients has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program Optimizer
  use CONTROL_
  use MolecularSystem_
  use String_
  use Exception_
  implicit none

       character(50) :: job
       integer :: i, j,k 
       integer :: numberOfParticles
       integer :: particleID

  job = ""  
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)
  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.BAS" )
  
  ! print *, "numberofquatumspecies", MolecularSystem_instance%numberOfQuantumSpecies
  ! do i = 1, MolecularSystem_instance%numberOfQuantumSpecies

  !    print *, "numberofparticles", size(MolecularSystem_instance%species(i)%particles) 
  !    do j = 1, size(MolecularSystem_instance%species(i)%particles)
  !       do k = 1, size(MolecularSystem_instance%species(i)%particles(j)%basis%contraction)

  !          print *, MolecularSystem_instance%species(i)%particles(j)%basis%contraction(k)%origin
  !       end do
  !    end do
  ! end do

  call system("lowdin-MolecularMechanics.x CONTROL_instance%FORCE_FIELD")

!!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
     !       do j = 1, size(MolecularSystem_instance%species(i)%particles)
     MolecularSystem_instance%species(i)%particles(1)%basis%origin=MolecularSystem_instance%species(i)%particles(1)%basis%origin-0.1
     do k = 1, size(MolecularSystem_instance%species(i)%particles(1)%basis%contraction)
        MolecularSystem_instance%species(i)%particles(1)%basis%contraction(k)%origin = &
             MolecularSystem_instance%species(i)%particles(1)%basis%contraction(k)%origin-0.1
        ! print *, MolecularSystem_instance%species(i)%particles(1)%basis%contraction(k)%origin
     end do
     !       end do
  end do

    numberOfParticles = size(MolecularSystem_instance%allParticles)
    particleID = (numberOfParticles/2) + 1
!    do i=1, size(MolecularSystem_instance%allParticles)
       MolecularSystem_instance%allParticles(1)%particlePtr%origin = &
            MolecularSystem_instance%allParticles(1)%particlePtr%origin-0.1
       MolecularSystem_instance%allParticles(particleID)%particlePtr%origin = &
            MolecularSystem_instance%allParticles(particleID)%particlePtr%origin-0.1
!    end do

  call MolecularSystem_saveToFile()

  call system("lowdin-MolecularMechanics.x CONTROL_instance%FORCE_FIELD")

  !!stop time
  call Stopwatch_stop(lowdin_stopwatch)
  
  write(*, *) ""
  write(*,"(A,F10.3,A4)") "** TOTAL Enlapsed Optimization Time : ", lowdin_stopwatch%enlapsetTime ," (s)"
  write(*, *) ""
  close(30)

end program Optimizer
