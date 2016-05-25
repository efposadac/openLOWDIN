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
!! @brief Integrals transformation program.
!!        This module allows to make calculations in the APMO framework
!!
!! @warning This programs only works linked to lowdincore library,
!!          provided by LOWDIN quantum chemistry package
!!
program IntegralTransformation
  use, intrinsic ::  iso_c_binding
  use CONTROL_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use ReadIntegrals_
  use Interface_
  implicit none

  integer :: i, j, k, nao, sze
  integer :: wfnUnit
  integer :: numberOfQuantumSpecies

  character(50) :: nameOfSpecies
  character(50) :: wfnFile
  character(50) :: arguments(2)

  type(Matrix) :: coefficients
  type(Vector) :: orb
  type(c_ptr) :: coeff_ptr, ints_ptr

  logical :: transformThisSpecies
  
  real(8), allocatable, target :: ints(:)
  real(8), allocatable, target :: coeff(:, :)

  print*, "Starting..."

  wfnFile = "lowdin.wfn"
  wfnUnit = 20

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


  !!*******************************************************************************************
  !! BEGIN
  open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 
  rewind(wfnUnit)

  numberOfQuantumSpecies = MolecularSystem_getNumberOfQuantumSpecies()

  do i=1, numberOfQuantumSpecies

    nameOfSpecies = MolecularSystem_getNameOfSpecie(i)

    if ( .not.CONTROL_instance%OPTIMIZE .and. transformThisSpecies) then

      write (6,"(T2,A)")"Integrals transformation for: "//trim(nameOfSpecies)

    end if

    !! Reading the coefficients
    nao = MolecularSystem_getTotalNumberOfContractions(i)

    arguments(2) = nameOfSpecies
    arguments(1) = "COEFFICIENTS"

    coefficients = Matrix_getFromFile(unit=wfnUnit, rows=int(nao, 4), &
                               columns=int(nao, 4), binary=.true., &
                               arguments=arguments(1:2))

    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(nao, nao))

    do j = 1, nao
      do k = 1, nao
        coeff(j, k) = coefficients%values(j, k)
      end do
    end do

    arguments(1) = "ORBITALS"

    call Vector_getFromFile(elementsNum=nao, unit=wfnUnit, binary=.true., &
                            arguments=arguments(1:2), output=orb)

    !! Read Integrals
    sze = nao * (nao + 1) / 2
    sze = sze * (sze + 1) / 2

    if(allocated(ints)) deallocate(ints)
    allocate(ints(sze))

    call ReadIntegrals_intraSpecies(trim(nameOfSpecies), ints, CONTROL_instance%NUMBER_OF_CORES)

    ! Calling C function
    coeff_ptr = c_loc(coeff(1, 1))
    ints_ptr = c_loc(ints(1))

    call c_test(coeff_ptr, ints_ptr, nao)

  end do
  

end program IntegralTransformation

