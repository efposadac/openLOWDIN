module DFTBCore_
  use CONTROL_
  use Exception_
  use MolecularSystem_

contains
  !>
  !! @brief Construye el archivo de entrada DFTB+
  subroutine DFTBCore_ini()

       print *, " Con DFTB plus"


    call MolecularSystem_showCartesianMatrix()
  end subroutine DFTBCore_ini
end module DFTBCore_


