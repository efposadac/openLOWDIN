module DFTBCore_
  use CONTROL_
  use Exception_
  use MolecularSystem_

contains
  subroutine DFTBCore_ini()
    if ( CONTROL_instance%dftbplus .eqv. .true. ) then
       print *, " Con DFTB plus"
    else
       print *, "Sin DFTB plus"
    end if

    call MolecularSystem_showCartesianMatrix()
  end subroutine DFTBCore_ini
end module DFTBCore_


