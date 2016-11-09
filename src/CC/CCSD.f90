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
      
      
      ! real(8), allocatable :: Dai(:,:)
      real(8), allocatable :: Tssame(:,:)
      real(8), allocatable :: Tdsame(:,:,:,:)
      real(8), allocatable :: tau(:,:,:,:)
      real(8), allocatable :: ttau(:,:,:,:)
      real(8) :: sum

      logical :: isInstanced

  end type CCSD

  type, public :: CCSDiter
      

      real(8), allocatable :: Dai(:,:)

  end type CCSDiter

  type(CCSD), public :: CCSD_instance
  type(CCSDiter), public :: CCSDinit


contains

  !>
  ! @brief Constructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSD_constructor()
      implicit none

      integer noc, nocs, nop, nops

      noc = CoupledCluster_instance%noc
      nocs = CoupledCluster_instance%nocs
      nop = CoupledCluster_instance%nop
      nops = CoupledCluster_instance%nops

      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine
      call CoupledCluster_pairing_function(1,2)

      ! Denominator in T1 D^{a}_{i}
      if (allocated(CCSDinit%Dai)) deallocate (CCSDinit%Dai)
      allocate(CCSDinit%Dai(CoupledCluster_instance%noc,CoupledCluster_instance%noc))
      CCSDinit%Dai(:,:) = 0.0_8

      ! t^{a}_{i} amplitude for single excitation
      if (allocated(CCSD_instance%Tssame)) deallocate(CCSD_instance%Tssame)
      allocate(CCSD_instance%Tssame(noc-nop,nop)) ! 01f
      CCSD_instance%Tssame(:,:) = 0.0_8

      ! t^{ab}_{ij} amplitude for double excitation for same species
      if (allocated(CCSD_instance%Tdsame)) deallocate(CCSD_instance%Tdsame)
      allocate(CCSD_instance%Tdsame(noc-nop,noc-nop,nop,nop)) ! 01f
      CCSD_instance%Tdsame(:,:,:,:) = 0.0_8

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tilde{\tau}
      if (allocated(CCSD_instance%ttau)) deallocate (CCSD_instance%ttau)
      allocate(CCSD_instance%ttau(noc-nop,noc-nop,nop,nop))
      CCSD_instance%ttau(:,:,:,:) = 0.0_8

      ! \tau
      if (allocated(CCSD_instance%tau)) deallocate (CCSD_instance%tau)
      allocate(CCSD_instance%tau(noc-nop,noc-nop,nop,nop))
      CCSD_instance%tau(:,:,:,:) = 0.0_8


      CCSD_instance%isInstanced = .true.
      
  end subroutine CCSD_constructor

  !>
  ! @brief Destructor of the class
  ! @author CAOM
  subroutine CCSD_destructor()
      implicit none

      ! if (allocated(CCSD_instance%Dai)) deallocate (CCSD_instance%Dai)
      if (allocated(CCSD_instance%Tssame)) deallocate (CCSD_instance%Tssame)
      ! if (allocated(CCSD_instance%Tdsame)) deallocate (CCSD_instance%Tdsame)
      ! if (allocated(CCSD_instance%ttau)) deallocate (CCSD_instance%ttau)
      ! if (allocated(CCSD_instance%tau)) deallocate (CCSD_instance%tau)

      CCSD_instance%isInstanced = .false.
      
  end subroutine CCSD_destructor

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information
  ! @author CAOM
  subroutine CCSD_init()
      implicit none

      integer noc, nocs, nop, nops, ax
      integer :: a, b, i, j
      real(8) :: ccsdE=0.0_8
      real(8) :: convergence = 1.0_8
      real(8) :: prev_ccsdE
      integer :: speciesId=1

      noc = CoupledCluster_instance%noc
      nocs = CoupledCluster_instance%nocs
      nop = CoupledCluster_instance%nop
      nops = CoupledCluster_instance%nops

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tau^{ab}_{ij} = t^{ab}_{ij} + \frac{1}{2}(t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j})
      ! \tau^{ab}_{ij} = t^{ab}_{ij} + t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j}

      print*, "before loop"
      do a=nop+1, noc
        do b=nop+1, noc
          do i=1, nop
            do j=1, nop

              CCSD_instance%Tdsame(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                +( (spints(speciesId)%valuesp(i,j,a,b))/( Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(speciesId)%HF_fs%values(j,j) &
                  -Allspecies(speciesId)%HF_fs%values(a,a)-Allspecies(speciesId)%HF_fs%values(b,b) ) ) 
              
              !under construction        
              
              CCSD_instance%ttau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) + 0.5*( CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) &
                -CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j) )

              CCSD_instance%tau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) !+ CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) &
                ! -CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j)

              write(*,*) CCSD_instance%Tdsame(a,b,i,j), "Tdsame"
            end do
          end do
        end do
      end do

      print*, "CCSD_init():"
      ! Denominator D^{a}_{i}
      do a=nop+1, noc
         do i=1, nop
            CCSDinit%Dai(a,i) = Allspecies(speciesId)%HF_fs%values(i,i) - Allspecies(speciesId)%HF_fs%values(a,a)
            write(*,*) a,i,CCSDinit%Dai(a,i)
         end do
      end do


      ! call Vector_destructor (Allspecies(speciesId)%HF_ff)
      ! call Matrix_destructor (Allspecies(speciesId)%HF_fs)
      ! Loop to obtain T1 and T2 intermediates values

      ! This could be a subroutine: 

      ! do while (convergence >= 1.0D-8)
          
      !   prev_ccsdE = ccsdE


      !   !intermediates loop

      !   convergence = abs( ccsdE - prev_ccsdE )

      !   write (*,*) speciesID, "Species"
      !   write (*,*) ccsdE, "CCSD Energy " 

      ! end do


  end subroutine CCSD_init

  subroutine CCSD_run()
      implicit none
      
      call CCSD_init()
      call CCSD_show()
      print*, "CCSD_show()"
      
  end subroutine CCSD_run

  subroutine CCSD_show()
      implicit none

      print*, "INFORMATION IN CCSD_constructor() HF_energy: ", CoupledCluster_instance%HF_energy
      print*, "INFORMATION IN CCSD_constructor() MP2_energy: ", CoupledCluster_instance%MP2_EnergyCorr
      CCSD_instance%sum = CoupledCluster_instance%HF_energy + CoupledCluster_instance%MP2_EnergyCorr
      print*, "INFORMATION IN CCSD_constructor() Total_energy: ", CCSD_instance%sum
      print*, "INFORMATION IN CCSD_constructor() Td: ", CCSD_instance%Tdsame(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() Dai: ", CCSDinit%Dai(1,3)
      ! print*, "INFORMATION IN CCSD_constructor() ttau: ",CCSD_instance%ttau(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() tau: ",CCSD_instance%tau(1,1,1,1)

      !if (allocated(CCSD_instance%tau)) deallocate (CCSD_instance%tau)

      !call CCSD_destructor()
      print*, "CCSD_show()x2"

  end subroutine CCSD_show
end module CCSD_