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
!! @brief COSMO and COSMO-APMO program.
!!        This module allows to make calculations in the COSMO-APMO framework
!! @author D. Gonzalez.
!!
!! <b> Creation date : </b> 2014-21-08
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos  para calculos con solvente implicito
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program Cosmo
  use CosmoTools_
  use CONTROL_
  use MolecularSystem_
  use Matrix_
  use String_
	use CosmoCore_
  implicit none 

  
	integer(8) :: n
	integer :: np

	
	type(Matrix) :: cmatin
	type(Matrix) :: qc
	type(Matrix) :: qq

	n=4

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

	!cmatin es el dummy de cmatinv
	call CosmoTools_constructor(surfaceSegment_instance,cmatin)	
	
	!call CosmoTools_quantum(n,np)

	call CosmoTools_clasical(surfaceSegment_instance,n,cmatin,qc)

	!write(*,*)quantum charges
	
	call CosmoTools_quantum(surfaceSegment_instance,n,cmatin,qq)

end program Cosmo

