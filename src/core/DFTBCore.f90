module DFTBCore_
  use CONTROL_
  use Exception_
  use MolecularSystem_
  use Particle_
  use ParticleManager_

contains
  !>
  !! @brief Construye el archivo de entrada DFTB+
  subroutine DFTBCore_ini()

    call MolecularSystem_showHSDFormat()
    
  end subroutine DFTBCore_ini
end module DFTBCore_


