!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!    http://www.qcc.unal.edu.co/
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Coupled Cluster Singles Doubles module
!!        This module contains the iteration of equations of Coupled Cluster Singles Doubles over APMO approach (CCSD-APMO).
!! @author  Carlos Andres Ortiz Mahecha (CAOM) (caraortizmah@unal.edu.co)
!!
!! <b> Creation date : </b> 2016-10-26
!!
!! <b> History: </b>
!!
!!   - <tt> 2016-10-26 </tt>: (CAOM) ( caraortizmah@unal.edu.co )
!!        -# Development of CCSD-APMO module:
!!                This Program calls xxx modules that iterate the intermediates for amplitude equations of CCSD-APMO.
!!   - <tt> data </tt>:  
!!
!!
!! @warning <em>  All characters and events in this module -- even those based on real source code -- are entirely fictional. </br>
!!                All celebrity lines are impersonated.....poorly. </br> 
!!                The following module contains corase language and due to it's cintent should not be viewed by anyone. </em>
!!
!!
!!
module CCSD_
  use MolecularSystem_
  use CoupledCluster_
  implicit none

  type, public :: CCSD
      
      
      real(8) :: sum
      real(8), allocatable :: Dai(:,:)

      logical :: isInstanced

  end type CCSD

  type(CCSD), public :: CCSD_instance


contains

  !>
  ! @brief Constructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSD_constructor()
      implicit none

      CCSD_instance%isInstanced = .true.

      ! allocate all that you can...
      ! First denominator in T2
      if (allocated(CCSD_instance%Dai)) deallocate (CCSD_instance%Dai)
      allocate(CCSD_instance%Dai(CoupledCluster_instance%noc,CoupledCluster_instance%noc))
      CCSD_instance%Dai(:,:) = 0.0_8
      
  end subroutine CCSD_constructor

  !>
  ! @brief Destructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSD_destructor()
      implicit none

      CCSD_instance%isInstanced = .false.
      
  end subroutine CCSD_destructor

  subroutine CCSD_init()
      implicit none

      integer :: a, i

      do a=CoupledCluster_instance%nop+1, CoupledCluster_instance%noc
         do i=1, CoupledCluster_instance%nop
            CCSD_instance%Dai(a,i) = CoupledCluster_instance%HF_fs%values(i,i) - CoupledCluster_instance%HF_fs%values(a,a)
            write(*,*) a,i,CCSD_instance%Dai(a,i)
         end do
      end do
      
  end subroutine CCSD_init

  subroutine CCSD_run()
      implicit none
      
      print*, "CCSD_init(): "
      call CCSD_init()
      call CCSD_show()
      
  end subroutine CCSD_run

  subroutine CCSD_show()
      implicit none

      print*, "INFORMATION IN CCSD_constructor() HF_energy: ", CoupledCluster_instance%HF_energy
      print*, "INFORMATION IN CCSD_constructor() MP2_energy: ", CoupledCluster_instance%MP2_EnergyCorr
      CCSD_instance%sum = CoupledCluster_instance%HF_energy + CoupledCluster_instance%MP2_EnergyCorr
      print*, "INFORMATION IN CCSD_constructor() Total_energy: ", CCSD_instance%sum
      call CoupledCluster_pairing_function(1)

  end subroutine CCSD_show
end module CCSD_