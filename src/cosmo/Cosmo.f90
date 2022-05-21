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
!! <b> Creation date : </b> 2014-08-21
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-08-21 </tt>: Danilo Gonzalez F. ( dagonzalezfo@unal.edu.co )
!!        -# Creacion de modulo y procedimientos  para calculos con solvente implicito
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program Cosmo
  use CONTROL_
  use MolecularSystem_
  use Matrix_
  use String_
  use CosmoCore_

  implicit none 

  integer(8) :: n

  type(Matrix) :: cmatin
  type(Matrix) :: qc
  type(Matrix) :: qq


  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )

  !cmatin es el dummy de cmatinv
  call CosmoCore_constructor(surfaceSegment_instance,cmatin)	

  n=MolecularSystem_instance%numberOfParticles

  call CosmoCore_clasical(surfaceSegment_instance,n,cmatin,qc)

  call system(" lowdin-ints.x COSMO ")

end program Cosmo

