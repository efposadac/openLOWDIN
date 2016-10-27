module CCSD_
  use MolecularSystem_
  use CoupledCluster_
  implicit none

contains
    
  subroutine CCSD_constructor()
      implicit none
      print*, "INFORMATION IN CCSD_constructor() HF_energy: ", CoupledCluster_instance%HF_energy
      call CoupledCluster_pairing_function(1)
      
  end subroutine CCSD_constructor
end module CCSD_