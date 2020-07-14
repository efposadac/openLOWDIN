
!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

!>
!!
!! This program organizes the modules necessary to calcutate Wavefuncion plots and
!! to write the HF wavefunction to the molden format, implemented in Lowdin1
!!
!! @author Jorge Charry ( jacharrym@unaledu.co ), Mauricio Rodas (jmrodasr@unal.edu.co)
!!
!! <b> Fecha de creacion : </b> 2014-01-31
!!
!! <b> Historial de modificaciones: </b>
!!
!<

program Output_
  use MolecularSystem_
  use Matrix_
  use InputOutput_
  use OutputBuilder_
  implicit none

  character(50) :: job
  integer :: numberOfOutputs, i
  type(OutputBuilder), allocatable :: outputs(:)
  
  job = ""
  call get_command_argument(1,value=job)  
  job = trim(String_getUppercase(job))

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  if(job.eq."FCHK") then
     
     allocate(outputs(1) )

     call OutputBuilder_constructor( outputs(1), 1, &
          "fchkFile", "ALL")
     
     call OutputBuilder_buildOutput(outputs(1))
     call OutputBuilder_show(outputs(1))
     
  else
     read(job,"(I10)") numberOfOutputs

       call InputOutput_constructor( numberOfOutputs )
       call InputOutput_load( )

       allocate(outputs(numberOfOutputs) )
  
       do i=1, numberOfOutputs
          call OutputBuilder_constructor( outputs(i), i, &
               InputOutput_Instance(i)%type, &
               InputOutput_Instance(i)%species, & 
               InputOutput_Instance(i)%state, &
               InputOutput_Instance(i)%orbital, &
               InputOutput_Instance(i)%dimensions, &
               InputOutput_Instance(i)%cubeSize, &
               InputOutput_Instance(i)%point1, & 
               InputOutput_Instance(i)%point2, &
               InputOutput_Instance(i)%point3  )

     call OutputBuilder_buildOutput(outputs(i))
     call OutputBuilder_show(outputs(i))

  end do

  

end program Output_





