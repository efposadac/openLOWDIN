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
module CCSDPT_
  use MolecularSystem_
  use CoupledCluster_
  use CCSD_
  use Tensor_
  use omp_lib
  implicit none

  type, public :: CCSD_T
      
      
      integer :: max
      integer :: min
      integer :: cont, aux_cont
      integer :: num_i, e_cont
      real(8) :: suma
      real(8), allocatable :: convergence_same(:)
      real(8), allocatable :: convergence_same_amp(:)
      real(8), allocatable :: convergence_diff(:)
      real(8), allocatable :: convergence_diff_amp(:)

      logical :: isInstanced

  end type CCSD_T

  type, public :: CCSD_Titer
      

      real(8), allocatable :: Dai(:,:)
      real(8), allocatable :: Fac(:,:)
      real(8), allocatable :: Fki(:,:)
      real(8), allocatable :: Fkc_aba(:,:)
      real(8), allocatable :: Fkc_aa(:,:)
      real(8), allocatable :: Fkca_ab(:,:)
      real(8), allocatable :: Faca(:,:)
      real(8), allocatable :: Fkia(:,:)
      real(8), allocatable :: Fbcb(:,:)
      real(8), allocatable :: Fkjb(:,:)
      real(8), allocatable :: Wklij(:,:,:,:)
      real(8), allocatable :: Wabcd(:,:,:,:)
      real(8), allocatable :: Wkbcj(:,:,:,:)
      real(8), allocatable :: Wkkcc(:,:,:,:)
      real(8), allocatable :: Waka(:,:,:,:)
      real(8), allocatable :: Wcia(:,:,:,:)
      real(8), allocatable :: Wbkb(:,:,:,:)
      real(8), allocatable :: Wcjb(:,:,:,:)
      real(8), allocatable :: Wakic_a(:,:,:,:)
      real(8), allocatable :: Wbkjc_b(:,:,:,:)
      real(8), allocatable :: Wklcd_b(:,:,:,:)
      real(8), allocatable :: Wklcd_a(:,:,:,:)
      real(8), allocatable :: Wakic(:,:,:,:)
      real(8), allocatable :: Wbkjc(:,:,:,:)
      real(8), allocatable :: Tai(:,:)
      real(8), allocatable :: Tai_AB(:,:)
      real(8), allocatable :: Taibj_int(:,:,:,:)
      real(8), allocatable :: Tabij(:,:,:,:)
      real(8), allocatable :: Tabij_AB(:,:,:,:)

  end type CCSD_Titer

  type(CCSD_T), public :: CCSDPT_instance
  type(CCSD_Titer), public, allocatable :: CCSDPTinit(:)
  type(CCSD_Titer), public, allocatable :: CCSDPT_loop(:)
  type(CCSD_Titer), public, allocatable :: CCSDPT_T1T2(:)
  type(CCSD_Titer), public, allocatable :: CCSDPT_inter(:)


contains

  !>
  ! @brief Constructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSDPT_constructor(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nop
      
      ! call CoupledCluster_pairing_function(1,2)
      
      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSDPT_constructor: noc=", noc, "nop=", nop
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! ! t^{ab}_{ij} amplitude for double excitation for same species
      if (allocated(Allspecies(speciesId)%Tt4same)) deallocate(Allspecies(speciesId)%Tt4same)
      allocate(Allspecies(speciesId)%Tt4same(noc-nop,noc-nop,noc-nop,nop,nop,nop))
      Allspecies(speciesId)%Tt4same(:,:,:,:,:,:) = 0.0_8

      if (allocated(Allspecies(speciesId)%Tt5same)) deallocate(Allspecies(speciesId)%Tt5same)
      allocate(Allspecies(speciesId)%Tt5same(noc-nop,noc-nop,noc-nop,nop,nop,nop))
      Allspecies(speciesId)%Tt5same(:,:,:,:,:,:) = 0.0_8

      CCSDPT_instance%isInstanced = .true.
  end subroutine CCSDPT_constructor

  !>
  ! @brief Destructor of the class
  ! @author CAOM
  subroutine CCSDPT_destructor()
      implicit none

      ! if (allocated(CCSDPT_instance%Dai)) deallocate (CCSDPT_instance%Dai)
      ! if (allocated(Allspecies(speciesId)%Tssame)) deallocate (Allspecies(speciesId)%Tssame)
      ! if (allocated(Allspecies(speciesId)%Tdsame)) deallocate (Allspecies(speciesId)%Tdsame)
      ! if (allocated(Allspecies(speciesId)%ttau)) deallocate (Allspecies(speciesId)%ttau)
      ! if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)

      CCSDPT_instance%isInstanced = .false.
  end subroutine CCSDPT_destructor

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information
  ! @author CAOM
  subroutine CCSDPT_init(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nop

      integer :: a, b, c, d, i, j, k, l

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tilde{\tau}^{ab}_{ij} = t^{ab}_{ij} + \frac{1}{2}(t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j})
      ! \tau^{ab}_{ij} = t^{ab}_{ij} + t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j}

      ! print*, "before loop"
      ! print*, Allspecies(speciesId)%Tdsame(1,1,1,1)
      ! print*, Allspecies(speciesId)%Tssame(1,1)
      ! print*, Allspecies(speciesId)%ttau(1,1,1,1)
      ! print*, Allspecies(speciesId)%tau(1,1,1,1)
      ! ! Allspecies(speciesId)%Tdsame(1,1,1,1) = Allspecies(speciesId)%tau(1,1,1,1)
      ! print*, "before loop"
      if (nop>=3) then
        do a=nop+1, noc
          do b=nop+1, noc
            do c=nop+1, noc
              do i=1, nop
                do j=1, nop
                  do k=1, nop


                    do d=nop+1, noc

                      Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) &
                        + (spints(speciesId)%valuesp(d,i,b,c)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,j,k)) &
                          - (spints(speciesId)%valuesp(d,i,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,j,k)) &
                            - (spints(speciesId)%valuesp(d,i,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,j,k)) &
                              - (spints(speciesId)%valuesp(d,j,b,c)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,k,i)) &
                                + (spints(speciesId)%valuesp(d,j,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,k,i)) &
                                  + (spints(speciesId)%valuesp(d,j,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,k,i)) &
                                    - (spints(speciesId)%valuesp(d,k,b,c)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,i,j)) &
                                      + (spints(speciesId)%valuesp(d,k,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,i,j)) &
                                        + (spints(speciesId)%valuesp(d,k,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,i,j))

                    end do
                    do l=1, nop

                      Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) &
                        - (spints(speciesId)%valuesp(l,a,j,k)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,i,l)) &
                          + (spints(speciesId)%valuesp(l,b,j,k)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,i,l)) &
                            + (spints(speciesId)%valuesp(l,c,j,k)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,l)) &
                              + (spints(speciesId)%valuesp(l,a,k,i)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,j,l)) &
                                - (spints(speciesId)%valuesp(l,b,k,i)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,j,l)) &
                                  - (spints(speciesId)%valuesp(l,c,k,i)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,l)) &
                                    + (spints(speciesId)%valuesp(l,a,i,j)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,k,l)) &
                                      - (spints(speciesId)%valuesp(l,b,i,j)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,k,l)) &
                                        - (spints(speciesId)%valuesp(l,c,i,j)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,k,l))

                      ! print*, "Allspecies(speciesId)%Tt4same: ", spints(speciesId)%valuesp(a,b,i,j), Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j), &
                      !   spints(speciesId)%valuesp(a,b,i,j)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j)
                    end do

                  Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k) &
                    + (spints(speciesId)%valuesp(j,k,b,c)*Allspecies(speciesId)%Tssame(a-nop,i)) &
                      - (spints(speciesId)%valuesp(j,k,c,a)*Allspecies(speciesId)%Tssame(b-nop,i)) &
                        - (spints(speciesId)%valuesp(j,k,a,b)*Allspecies(speciesId)%Tssame(c-nop,i)) &
                          - (spints(speciesId)%valuesp(k,i,b,c)*Allspecies(speciesId)%Tssame(a-nop,j)) &
                            + (spints(speciesId)%valuesp(k,i,c,a)*Allspecies(speciesId)%Tssame(b-nop,j)) &
                              + (spints(speciesId)%valuesp(k,i,a,b)*Allspecies(speciesId)%Tssame(c-nop,j)) &
                                - (spints(speciesId)%valuesp(i,j,b,c)*Allspecies(speciesId)%Tssame(a-nop,k)) &
                                  + (spints(speciesId)%valuesp(i,j,c,a)*Allspecies(speciesId)%Tssame(b-nop,k)) &
                                    + (spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tssame(c-nop,k))


                  !***
                  Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)/ &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                      + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))

                  Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k)/ &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                      + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c)) 

                  ! print*, "Dabcijk: ", Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                  !     + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                  !       - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c)  
                  ! print*, "tt4same: ", Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)

                  if ((Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                      + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))==0) then
                    Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)=0.0_8
                    Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k)=0.0_8
                  end if
              
                    ! write(*,*) Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j), "Tdsame"
                    ! , Allspecies(speciesId)%HF_fs%values(a,a), "HF_fs%values", spints(speciesId)%valuesp(i,j,a,b), "spints"

                  end do
                end do
              end do
            end do
          end do
        end do
      end if
      
      ! print*, "CCSDPT_init(", speciesId, ")"
      ! ! Denominator D^{a}_{i}
      ! do a=nop+1, noc
      !   do i=1, nop
      !     CCSDinit(speciesId)%Dai(a,i) = Allspecies(speciesId)%HF_fs%values(i,i) - Allspecies(speciesId)%HF_fs%values(a,a)
      !     ! write(*,*) a,i,CCSDinit(speciesId)%Dai(a,i)
      !   end do
      ! end do
      ! do a=nop+1, noc
      !   do i=1, nop
      !     write(*,*) a,i,CCSDinit(speciesId)%Dai(a,i)
      !     print*, "HF_fs: ", Allspecies(speciesId)%HF_fs%values(i,i), Allspecies(speciesId)%HF_fs%values(a,a)
      !   end do
      ! end do

  end subroutine CCSDPT_init

  !>
  ! @brief Build a intermediates that will be used in Coupled Cluster loop
  ! @author CAOM
  ! subroutine CCSD_loop_constructor(speciesId)
  !     implicit none

  !     integer, intent(in) :: speciesId
  !     integer noc, nop
      
  !     noc = Allspecies(speciesId)%noc
  !     ! nocs = CoupledCluster_instance%nocs
  !     nop = Allspecies(speciesId)%nop
  !     ! nops = CoupledCluster_instance%nops
      
  !     write(*, *) " CCSD_loop_constructor: noc=", noc, " nop=", nop

  !     !
  !     if (allocated(CCSDloop(speciesId)%Fac)) deallocate (CCSDloop(speciesId)%Fac)
  !     allocate(CCSDloop(speciesId)%Fac(noc-nop,noc-nop))
  !     CCSDloop(speciesId)%Fac=0.0_8

  !     !
  !     if (allocated(CCSDloop(speciesId)%Fki)) deallocate (CCSDloop(speciesId)%Fki)
  !     allocate(CCSDloop(speciesId)%Fki(nop,nop))
  !     CCSDloop(speciesId)%Fki=0.0_8
      
  !     !
  !     if (allocated(CCSDloop(speciesId)%Fkc_aa)) deallocate (CCSDloop(speciesId)%Fkc_aa)
  !     allocate(CCSDloop(speciesId)%Fkc_aa(nop,noc-nop))
  !     CCSDloop(speciesId)%Fkc_aa=0.0_8

  !     !
  !     if (allocated(CCSDloop(speciesId)%Wklij)) deallocate (CCSDloop(speciesId)%Wklij)
  !     allocate(CCSDloop(speciesId)%Wklij(nop,nop,nop,nop))
  !     CCSDloop(speciesId)%Wklij=0.0_8

  !     !
  !     if (allocated(CCSDloop(speciesId)%Wabcd)) deallocate (CCSDloop(speciesId)%Wabcd)
  !     allocate(CCSDloop(speciesId)%Wabcd(noc-nop,noc-nop,noc-nop,noc-nop))
  !     CCSDloop(speciesId)%Wabcd=0.0_8

  !     !
  !     if (allocated(CCSDloop(speciesId)%Wkbcj)) deallocate (CCSDloop(speciesId)%Wkbcj)
  !     allocate(CCSDloop(speciesId)%Wkbcj(nop,noc-nop,noc-nop,nop))
  !     CCSDloop(speciesId)%Wkbcj=0.0_8
      
  !     ! stop"ME"
  ! end subroutine CCSD_loop_constructor

  !>
  ! @brief Build T1 and T2 amplitude equations that will be information of intermediates
  ! @author CAOM
  ! subroutine CCSD_T1T2_constructor(speciesId)
  !     implicit none

  !     integer, intent(in) :: speciesId
  !     integer noc, nop

  !     noc = Allspecies(speciesId)%noc
  !     ! nocs = CoupledCluster_instance%nocs
  !     nop = Allspecies(speciesId)%nop
  !     ! nops = CoupledCluster_instance%nops
  !     ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_constructor: noc=", noc, "nop=", nop     
  !     print*, " CCSD_T1T2_constructor(: ", nop, noc, speciesId
      
  !     !
  !     if (allocated(CCSDT1T2(speciesId)%Tai)) deallocate (CCSDT1T2(speciesId)%Tai)
  !     allocate(CCSDT1T2(speciesId)%Tai(noc-nop,nop))
  !     CCSDT1T2(speciesId)%Tai=0.0_8

  !     !
  !     if (allocated(CCSDT1T2(speciesId)%Tai_AB)) deallocate (CCSDT1T2(speciesId)%Tai_AB)
  !     allocate(CCSDT1T2(speciesId)%Tai_AB(noc-nop,nop))
  !     CCSDT1T2(speciesId)%Tai_AB=0.0_8

  !     !
  !     if (allocated(CCSDT1T2(speciesId)%Taibj_int)) deallocate (CCSDT1T2(speciesId)%Taibj_int)
  !     allocate(CCSDT1T2(speciesId)%Taibj_int(noc-nop,noc-nop,nop,nop))
  !     CCSDT1T2(speciesId)%Taibj_int=0.0_8
      
  !     !
  !     if (allocated(CCSDT1T2(speciesId)%Tabij)) deallocate (CCSDT1T2(speciesId)%Tabij)
  !     allocate(CCSDT1T2(speciesId)%Tabij(noc-nop,noc-nop,nop,nop))
  !     CCSDT1T2(speciesId)%Tabij=0.0_8
  ! end subroutine CCSD_T1T2_constructor

  !>
  ! @brief Calculate F intermediates for intra-species
  ! @author CAOM
  ! subroutine F_onespecies_intermediates(speciesId)
  !     implicit none

  !     integer, intent(in) :: speciesId

  !     integer :: noc, nop
  !     integer :: a, e, i, f, m, n

  !     noc = Allspecies(speciesId)%noc
  !     nop = Allspecies(speciesId)%nop

  !     print*,"noc: ", noc, "nop: ", nop
  !     print*, "speciesId: ", speciesId
  !     ! if (speciesId>2) stop "inside loop"

  !     ! CCSDloop(speciesId)%Fac
  !     do a=nop+1, noc
  !       do e=nop+1, noc

  !         CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) & 
  !           +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e)
            
  !         do m=1, nop
  !           CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
  !             + (-0.5*Allspecies(speciesId)%HF_fs%values(m,e)*Allspecies(speciesId)%Tssame(a-nop,m))
            
  !           if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
  !             do f=nop+1, noc
  !               CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
  !                 + Allspecies(speciesId)%Tssame(f-nop,m)*spints(speciesId)%valuesp(m,a,f,e)
  !               do n=1, nop
  !                 CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
  !                   + (-0.5*Allspecies(speciesId)%ttau(a-nop,f-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f))
  !                 ! write(*,*) a,e,CCSDloop(speciesId)%Fac(a,e)
  !               end do
  !             end do
  !           end if
  !         end do
  !       end do
  !     end do

  !     ! CCSDloop(speciesId)%Fki(m,i)
  !     do m=1, nop
  !       do i=1, nop
             
  !         CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
  !           + (1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i)
            
  !         do e=nop+1, noc
  !           CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
  !              + 0.5*Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%HF_fs%values(m,e)
            
  !           if (nop>=2) then ! kind of interaction just for two or more particles of the principal species  
  !             do n=1, nop
  !               CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
  !                 + Allspecies(speciesId)%Tssame(e-nop,n)*spints(speciesId)%valuesp(m,n,i,e)
  !               do f=nop+1, noc
  !                 CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
  !                   + 0.5*Allspecies(speciesId)%ttau(e-nop,f-nop,i,n)*spints(speciesId)%valuesp(m,n,e,f)
  !                 ! write(*,*) a,e,CCSDloop(speciesId)%Fki(a,e)
  !               end do
  !             end do
  !           end if
  !         end do
  !       end do
  !     end do

  !     ! CCSDloop(speciesId)%Fkc_aa
  !     do m=1, nop
  !       do e=nop+1, noc

  !         CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
  !           + Allspecies(speciesId)%HF_fs%values(m,e)
  !         if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
  !           do n=1, nop
  !             do f=nop+1, noc
  !               CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
  !                 + Allspecies(speciesId)%Tssame(f-nop,n)*spints(speciesId)%valuesp(m,n,e,f)
  !               ! write(*,*) a,e,CCSDloop(speciesId)%Fkc_aa(a,e)
  !             end do
  !           end do
  !         end if
  !       end do
  !     end do
  ! end subroutine F_onespecies_intermediates


  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSDPT_same_species(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer :: noc, nop
      integer :: num_species
      
      integer :: a, b, c, i, j, k
      real(8) :: prev_ccsdE
      real(8) :: prev_ampl
      real(8) :: tmp_ccsd_t4E
      real(8) :: tmp_ccsd_ts5E
      real(8) :: ccsdE=0.0_8
      real(8) :: ampl=0.0_8
      real(8) :: convergence = 1.0_8
      real(8) :: conver_amplitude = 1.0_8
      real(8) :: auxtdsame = 0.0_8
      real(8) :: auxtssame = 0.0_8

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (conver_amplitude /= 1.D-8) conver_amplitude = 1.0_8
      if (ccsdE /= 0.0_8) ccsdE = 0.0_8
      if (ampl /= 0.0_8) ampl = 0.0_8

      !Initialization of private variables from public variables
      noc = Allspecies(speciesId)%noc
      nop = Allspecies(speciesId)%nop
      num_species = CoupledCluster_instance%num_species
      ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
      print*, "T1T2_constructor", convergence, noc, nop, speciesId

      ! times_i = CoupledCluster_instance%times_intersp
      ! max = CCSD_instance%max
      ! min = CCSD_instance%min

      !   !**change position of do while
      !   !do i=1, num_species
      !     !all intra CCSD
      !   !end do
      !   !do i=min, max
      !   ! all inter-species
      !   !end do
      ! prev_ccsdE = e_ccsd
      ! prev_ampl = v_ampl

      ! print*, "speciesId CCSD_loop: ", speciesId

      ! call CCSD_T1T2_constructor(speciesId)
      ! call CCSD_loop_constructor(speciesId)
      ! ! stop"FER CCSD_same_species"
      ! !intermediates loop for:
      ! !singles excitations
      ! print*, "F_onespecies_intermediates(): "
      ! call F_onespecies_intermediates(speciesId)
      ! !doubles excitations
      ! if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
      !   print*, "W_onespecies_intermediates(): "
      !   call W_onespecies_intermediates(speciesId)
      ! end if

      ! Resolve CCSD equation of energy

      tmp_ccsd_t4E=0.0_8
      tmp_ccsd_ts5E=0.0_8

      ! for same species
        do a=nop+1, noc
          do b=nop+1, noc
            do c=nop+1, noc
              do i=1, nop
                do j=1, nop
                  do k=1, nop
                    tmp_ccsd_t4E = tmp_ccsd_t4E + ((0.0277777777)*Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)* &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                        + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))* &
                            Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k))
                    ! print*, "tmp_ccsd_t4E: ", tmp_ccsd_t4E, Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k), &
                    !   (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                    !     + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                    !       - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))

                    tmp_ccsd_ts5E = tmp_ccsd_ts5E + (0.25*Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)* &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                        + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))* &
                            Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k))

                  end do
                end do
              end do
            end do
          end do
        end do

      print *, "ccsd_t4E intra", tmp_ccsd_t4E
      print *, "ccsd_ts5E intra", tmp_ccsd_ts5E

      ! !change in values for the intra-species loop
      ! convergence = abs( ccsdE - prev_ccsdE )
      ! ! print*, "Otherspecies Tssame: ", Allspecies(speciesId)%Tssame
      ! write (*,*) ccsdE, "CCSD Energy same species", prev_ccsdE, "previous Energy" 
      ! write (*,*) convergence, "Convergence " 

      ! call CCSD_T1(speciesId)
      ! if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
      !   call CCSD_T2(speciesId)
      ! end if

      ! if (convergence > 100) then 
      !   stop "Error: There are not convergence. The differences between energies is more than 100 eV"
      ! end if

      CoupledCluster_instance%CCSD_T4_E_intra(speciesId) = tmp_ccsd_t4E
      CoupledCluster_instance%CCSD_TS5_E_intra(speciesId) = tmp_ccsd_ts5E
  end subroutine CCSDPT_same_species
  
  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  ! subroutine CCSD_diff_species(speciesId, OtherspeciesId, e_ccsd, v_ampl_int)
  !     implicit none

  !     integer, intent(in) :: speciesId
  !     integer, intent(in) :: OtherspeciesId
  !     real(8), intent(in) :: e_ccsd
  !     real(8), intent(in) :: v_ampl_int

  !     integer :: noc, nocs, nop, nops
  !     integer :: num_species
      
  !     integer :: max, e_cont
  !     integer :: a, i
  !     integer :: aa, ii
  !     integer :: p, q, r, s
  !     integer :: num_inter
  !     real(8) :: prev_ccsdE_int
  !     real(8) :: prev_ampl
  !     real(8) :: tmp_ccsdE_int
  !     real(8) :: convergence = 1.0_8
  !     real(8) :: conver_amplitude = 1.0_8
  !     real(8) :: ccsdE_int=0.0_8
  !     real(8) :: ampl=0.0_8
  !     real(8) :: auxtdsame = 0.0_8
  !     real(8) :: auxtssame = 0.0_8

  !     if (convergence /= 1.0D-8) convergence = 1.0_8
  !     if (conver_amplitude /= 1.D-8) conver_amplitude = 1.0_8
  !     if (ccsdE_int /= 0.0_8) ccsdE_int = 0.0_8
  !     if (ampl /= 0.0_8) ampl = 0.0_8

  !     !Initialization of private variables from public variables
  !     noc = Allspecies(speciesId)%noc
  !     nocs = Allspecies(OtherspeciesId)%noc
  !     nop = Allspecies(speciesId)%nop
  !     nops = Allspecies(OtherspeciesId)%nop
  !     num_species = CoupledCluster_instance%num_species
  !     ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
  !     print*, "T1T2_constructor", convergence, noc, nop, speciesId

  !     prev_ccsdE_int = e_ccsd
  !     prev_ampl = v_ampl_int

  !     print*, "speciesId: ", speciesId

  !     num_inter = CCSD_instance%num_i
  !     e_cont = CCSD_instance%e_cont

  !     ! call CCSD_constructor_inter(min, jj)
  !     print*, "OtherspeciesId: ", OtherspeciesId
  !     ! if (speciesId>OtherspeciesId) print*, "before CCSD_T2AB_constructor: ", Allinterspecies(speciesId)%Tdsame
  !     call CCSD_T2AB_constructor(speciesId, OtherspeciesId, num_inter)
  !     ! if (speciesId>OtherspeciesId) print*, "before CCSD_loop_constructor_inter: ", Allinterspecies(speciesId)%Tdsame
  !     call CCSD_loop_constructor_inter(speciesId, OtherspeciesId, num_inter)

  !     !intermediates loop for:
  !     !If there are interspecies?
  !     print*, "ciclo: OtherspeciesId: ", OtherspeciesId, "speciesId: ", speciesId
  !     ! if (speciesId>OtherspeciesId) print*, "before F_twospecies: ", Allinterspecies(speciesId)%Tdsame
  !     call F_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
  !     ! if (speciesId>OtherspeciesId) print*, "before W_twospecies: ", Allinterspecies(speciesId)%Tdsame
  !     call W_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
  !     ! if (speciesId>OtherspeciesId) print*, "before F_T2_AB: ", Allinterspecies(speciesId)%Tdsame
  !     call F_T2_AB(speciesId, OtherspeciesId, num_inter)
  !     ! if (speciesId>OtherspeciesId) print*, "between F_T2_AB and W_T2_AB: ", Allinterspecies(speciesId)%Tdsame
  !     call W_T2_AB(speciesId, OtherspeciesId, num_inter)
        
  !     ! Resolve CCSD equation of energy

  !     tmp_ccsdE_int=0.0_8
  !     print*, "energy CCSD-APMO: ", speciesId+1, max
  !     print*, "e_cont: ", e_cont
  !     auxtdsame = 0 
  !     auxtssame = 0 
  !     do i=1, nop
  !       do a=nop+1, noc
  !         do ii=1, nops
  !           do aa=nops+1, nocs

  !             if (speciesId<OtherspeciesId) then
  !               p=i
  !               q=ii
  !               r=a
  !               s=aa
  !             else
  !               p=ii
  !               q=i
  !               r=aa
  !               s=a
  !             end if

  !             tmp_ccsdE_int = tmp_ccsdE_int + (0.5*spintm(e_cont)%valuesp(p,q,r,s)* &
  !                 Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) ) &
  !                   + (0.5*spintm(e_cont)%valuesp(p,q,r,s)*Allinterspecies(speciesId)%Tssame(a-nop,i)* &
  !                     Allinterspecies(OtherspeciesId)%Tssame(aa-nops,ii))
  !             auxtdsame = auxtdsame + (0.5*spintm(e_cont)%valuesp(p,q,r,s)* Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) )
  !             auxtssame = auxtssame + (0.5*spintm(e_cont)%valuesp(p,q,r,s)*Allinterspecies(speciesId)%Tssame(a-nop,i)*Allinterspecies(OtherspeciesId)%Tssame(aa-nops,ii))

  !             ! print*, "spintm(e_cont): ", spintm(e_cont)%valuesp(p,q,r,s)
  !             ! if (speciesId>OtherspeciesId) print*, "CCSD diff Allinterspecies(speciesId)%Tdsame: ", &
  !             !   Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii)
  !             ! print*, "Allspecies(speciesId)%Tssame: ", Allspecies(speciesId)%Tssame(a-nop,i)
  !             ! print*, "Allspecies(OtherspeciesId)%Tssame: ", Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)
  !           end do
  !         end do
  !       end do
  !     end do
  !     print *, "speciesId", speciesId
  !     print *, "tdsame E inter", auxtdsame
  !     print *, "tssame E inter", auxtssame
  !     print *, "Tdsame: ",Allspecies(speciesId)%Tdsame(1,1,1,1)! + sum(Allspecies(speciesId)%Tdsame,dim=2) +&
  !     print *, "Tdsame inter: ", Allinterspecies(speciesId)%Tdsame(1,1,1,1)! + sum(Allinterspecies(speciesId)%Tdsame,dim=2) +&
  !     print *, "Tdsame diff: ", Allinterspecies(speciesId)%Tddiff(1,1,1,1)! + sum(Allinterspecies(speciesId)%Tddiff,dim=2) +&
  !     print *, "intau: ", Allinterspecies(speciesId)%intau(1,1,1,1)
  !     print *, "tau: ", Allinterspecies(speciesId)%tau(1,1,1,1)
  !     print *, "Tssame: ", Allspecies(speciesId)%Tssame(1,1) !+ sum(Allspecies(speciesId)%Tssame,dim=2)
  !     print *, "Fac: ", CCSDloop(speciesId)%Fac(1,1)
  !     print *, "Fac T1T2: ", CCSDT1T2(speciesId)%Fac(1,1)
  !     print *, "Fki: ", CCSDloop(speciesId)%Fki(1,1)
  !     print *, "Fki: T1T2", CCSDT1T2(speciesId)%Fki(1,1)
  !     print *, "Fkc_aa: ", CCSDloop(speciesId)%Fkc_aa(1,1)
  !     print *, "Fkc_aa: T1T2", CCSDT1T2(speciesId)%Fkc_aa(1,1)
  !     print *, "Fkc_aba: ", CCSDinter(speciesId)%Fkc_aba(1,1)
  !     print *, "Fkca_ab: ", CCSDinter(OtherspeciesId)%Fkca_ab(1,1)
  !     print *, "Faca: ", CCSDinter(speciesId)%Faca(1,1)
  !     print *, "Fkia: ", CCSDinter(speciesId)%Fkia(1,1)
  !     print *, "Fbcb: ", CCSDinter(OtherspeciesId)%Fbcb(1,1)
  !     print *, "Fkjb: ", CCSDinter(OtherspeciesId)%Fkjb(1,1)
  !     print *, "Waka: ", CCSDinter(speciesId)%Waka(1,1,1,1)
  !     print *, "Wcia: ", CCSDinter(speciesId)%Wcia(1,1,1,1)
  !     print *, "Wbkb: ", CCSDinter(speciesId)%Wbkb(1,1,1,1)
  !     print *, "Wcjb: ", CCSDinter(speciesId)%Wcjb(1,1,1,1)
  !     print *, "Wakic: ", CCSDinter(speciesId)%Wakic(1,1,1,1)
  !     print *, "Wbkjc: ", CCSDinter(OtherspeciesId)%Wbkjc(1,1,1,1)
  !     print *, "spintm: ", spintm(e_cont)%valuesp(1,1,1,1)


  !     print*, "energy CCSD-APMO"

  !     ccsdE_int = tmp_ccsdE_int

  !     !change in values for the inter-species loop
  !     convergence = abs( ccsdE_int - prev_ccsdE_int )
    
  !     ampl = sum(Allinterspecies(speciesId)%Tdsame)
  !     conver_amplitude = abs(ampl-prev_ampl)

  !     write (*,*) ccsdE_int, "CCSD Energy inter-species ", prev_ccsdE_int, "previous Energy inter-species" 
  !     write (*,*) convergence, "Convergence inter-species" 

  !     ! ! T1 and T2 equation energies for interspecies
  !     ! print*, "before CCSD_T1_inter()", speciesId, OtherspeciesId
  !     ! call CCSD_T1_inter(speciesId, OtherspeciesId, num_inter)
  !     ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T1 Allinterspecies(speciesId)%Tdsame: ", &
  !     ! !   Allinterspecies(speciesId)%Tdsame
  !     ! print*, "before CCSD_T2_inter()", speciesId, OtherspeciesId
  !     ! call CCSD_T2_inter(speciesId, OtherspeciesId, num_inter)
  !     ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2 Allinterspecies(speciesId)%Tdsame: ", &
  !     ! !   Allinterspecies(speciesId)%Tdsame
  !     ! print*, "CCSD_T2_AB()"
  !     ! call CCSD_T2_AB(speciesId, OtherspeciesId, num_inter)
  !     ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2AB Allinterspecies(speciesId)%Tdsame: ", &
  !     ! !   Allinterspecies(speciesId)%Tdsame
  !     ! print*, "num_inter: ", num_inter
  
  !     if (convergence > 100) then 
  !       stop "Error: There are not convergence. The differences between energies is more than 100 eV"
  !     end if


  !     CoupledCluster_instance%CCSD_E_inter(num_inter) = ccsdE_int
  !     CoupledCluster_instance%CCSD_A_inter(speciesId) = ampl
  !     CCSD_instance%convergence_diff(num_inter) = convergence
  !     CCSD_instance%convergence_diff_amp(speciesId) = conver_amplitude
      
  !     num_inter = num_inter + 1!final
  !     CCSD_instance%num_i = num_inter
  !     CCSD_instance%cont = CCSD_instance%cont + 1!final

  !     e_cont = e_cont + 1!final
  !     CCSD_instance%e_cont = e_cont
  ! end subroutine CCSD_diff_species


  !>
  ! @brief Calculate T1 energy equations for intra-species
  ! @author CAOM
  ! subroutine CCSD_T1(speciesId)
  !     implicit none

  !     integer, intent(in) :: speciesId

  !     integer noc, nop
  !     integer :: num_species
  !     integer :: times_i
  !     integer :: nthreads, threadid
  !     integer :: a, e, f, i, m, n
  !     real(8), allocatable :: tai_tmp_1(:)
  !     real(8), allocatable :: tai_tmp_2(:)
  !     real(8), allocatable :: tai_tmp_3(:)
  !     real(8), allocatable :: tai_tmp_4(:)
  !     real(8), allocatable :: tai_tmp_5(:)

  !     noc = Allspecies(speciesId)%noc
  !     ! nocs = CoupledCluster_instance%nocs
  !     nop = Allspecies(speciesId)%nop
  !     ! nops = CoupledCluster_instance%nops
  !     num_species = CoupledCluster_instance%num_species
  !     times_i = CoupledCluster_instance%times_intersp

  !     !
  !     if (allocated(tai_tmp_1)) deallocate(tai_tmp_1)
  !     if (allocated(tai_tmp_2)) deallocate(tai_tmp_2)
  !     if (allocated(tai_tmp_3)) deallocate(tai_tmp_3)
  !     if (allocated(tai_tmp_4)) deallocate(tai_tmp_4)
  !     if (allocated(tai_tmp_5)) deallocate(tai_tmp_5)

  !     !$OMP PARALLEL 
  !       nthreads = OMP_GET_NUM_THREADS()
  !     !$OMP END PARALLEL
  !     allocate(tai_tmp_1(0:nthreads))
  !     allocate(tai_tmp_2(0:nthreads))
  !     allocate(tai_tmp_3(0:nthreads))
  !     allocate(tai_tmp_4(0:nthreads))
  !     allocate(tai_tmp_5(0:nthreads))


  !     !Basic parallelization
  !     !$OMP PARALLEL shared(tai_tmp_1,tai_tmp_2,tai_tmp_3,tai_tmp_4,tai_tmp_5,speciesId,nop,noc)
  !     !$OMP& private(a,i,e,m,f,n,nthreads,threadid) collapse(2)
  !     ! T^{a}_{i}D^{a}_{i} = ...
  !     nthreads = OMP_GET_NUM_THREADS()
  !     threadid =  OMP_GET_THREAD_NUM()

  !     tai_tmp_1(threadid) = 0.0_8
      
  !     do a=nop+1, noc
  !       do i=1, nop
  !         CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !           + Allspecies(speciesId)%HF_fs%values(i,a) ! eq1 1

  !         do e=nop+1, noc
  !           CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !           ! tai_tmp_1(threadid) = tai_tmp_1(threadid) &
  !             + Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop(speciesId)%Fac(a-nop,e-nop) !eq 4 1

  !           ! tai_tmp_4_buff(threadid)
  !         end do
  !         if (threadid == 0) then
  !           CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + sum(tai_tmp_1,dim=1)
  !         end if
  !         !$OMP BARRIER
          
  !         tai_tmp_1(threadid) = 0.0_8
  !         tai_tmp_2(threadid) = 0.0_8
  !         tai_tmp_3(threadid) = 0.0_8
  !         tai_tmp_4(threadid) = 0.0_8
  !         tai_tmp_5(threadid) = 0.0_8                

  !         do m=1, nop
  !           CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !           ! tai_tmp_1(threadid) = tai_tmp_1(threadid) &
  !             + (-Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop(speciesId)%Fki(m,i)) !eq 5 2
          
  !           if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
              
  !             !$OMP CRITICAL
  !             do e=nop+1, noc
  !               CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !               ! tai_tmp_2(threadid) = tai_tmp_2(threadid) &
  !                 + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop(speciesId)%Fkc_aa(m,e-nop) !eq 6 3

  !               do f=nop+1, noc
  !                 CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !                 ! tai_tmp_3(threadid) = tai_tmp_3(threadid) &
  !                   + (-0.5*Allspecies(speciesId)%Tdsame(e-nop,f-nop,i,m)*spints(speciesId)%valuesp(m,a,e,f)) !eq 3 7
  !               end do

  !               do n=1, nop
  !                 CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !                 ! tai_tmp_4(threadid) = tai_tmp_4(threadid) &
  !                   + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,n)*spints(speciesId)%valuesp(n,m,e,i)) !eq 2 6
  !               end do
                

  !               ! if (threadid == 0) then
  !               !   tai_tmp_2(threadid) = tai_tmp_2(threadid) & + sum(tai_tmp_1,dim=1)
  !               ! end if

  !               ! tai_tmp_1(e-nop) = tai_tmp_1(e-nop) + sum(tai_tmp_2,dim=1) + sum(tai_tmp_3,dim=1)
  !             end do
  !             !$OMP END CRITICAL


  !           end if

  !         end do

  !         if (threadid == 0) then
  !           CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + sum(tai_tmp_1,dim=1)
  !             !+ sum(tai_tmp_2,dim=1) + sum(tai_tmp_3,dim=1) + sum(tai_tmp_4,dim=1)
  !         end if

  !         !$OMP CRITICAL
  !         ! CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + sum(tai_tmp_5,dim=1)
  !         if ((nop>=2)) then !  .and. (times_i==0) kind of interaction just for two or more particles of the principal species 
  !           !if times_i/=0 then the product below will be calculated in Fkc_aba in subroutine CCSD_T1_inter()
  !           do m=1,nop !n
  !             ! tai_tmp_2(threadid) = 0.0_8
  !             do e=nop+1, noc !f
  !               CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
  !               ! tai_tmp_2(threadid) = tai_tmp_2(threadid) &
  !                 + (-Allspecies(speciesId)%Tssame(e-nop,m)*spints(speciesId)%valuesp(m,a,i,e)) !eq 7 5 !n,a,i,f
  !             end do
  !             ! tai_tmp_3(m) = tai_tmp_3(m) + sum(tai_tmp_1,dim=1)
  !           end do
  !           ! CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + sum(tai_tmp_3,dim=1)
  !         end if
  !         !$OMP END CRITICAL

  !         ! if (threadid == 0) then
  !         !   CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + sum(tai_tmp_2,dim=1)
  !         ! end if

  !           CCSDT1T2(speciesId)%Tai_AB(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
  !         ! print*, "Tai: ", CCSDT1T2(speciesId)%Tai(a-nop,i), "Tai_AB: ", CCSDT1T2(speciesId)%Tai_AB(a-nop,i)
  !           CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit(speciesId)%Dai(a,i)
            
  !           if(CCSDinit(speciesId)%Dai(a,i)==0) CCSDT1T2(speciesId)%Tai(a-nop,i) = 0.0_8
            
  !           Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
  !         ! end if
  !         ! print*, "Tai: ", CCSDT1T2(speciesId)%Tai(a-nop,i)
  !         ! print*, "Dai: ", CCSDinit(speciesId)%Dai(a,i)
  !         ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
  !       end do
  !     end do
  !     ! print*, "CCSD_T1"
  !     ! print*, "Tai: ", CCSDT1T2(speciesId)%Tai
  !     ! print*, "Dai: ", CCSDinit(speciesId)%Dai
  !     !$OMP END PARALLEL
  ! end subroutine CCSD_T1

  !>
  ! @brief Manager of Coupled Cluster maths
  ! @author CAOM
  subroutine CCSDPT_run()
      implicit none

      integer :: n_intersp(10), i_counterID(10)
      integer :: num_species, counterID, finalID
      integer :: num_inter, num_i
      integer :: i, j
      integer :: times_i
      ! integer :: stoped=0
      integer :: aux_cont=0
      integer :: n_sp=0
      ! integer :: max, min
      ! real(8) :: intra=0
      ! real(8) :: convergence = 1.0_8
      ! real(8) :: convergence_int = 1.0_8
      ! real(8), allocatable :: e_same_ccd(:)
      ! real(8), allocatable :: a_same_ccd(:)
      ! real(8), allocatable :: e_diff_ccd(:)
      ! real(8), allocatable :: a_diff_ccd(:)

      call CCSD_run()
      
      !Initialization of variables
      num_species = CoupledCluster_instance%num_species
      ! finalID = CoupledCluster_instance%finalID
      counterID = CoupledCluster_instance%counterID
      num_inter = CoupledCluster_instance%num_intersp
      times_i = CoupledCluster_instance%times_intersp
      ! CCSD_instance%max=0
      ! CCSD_instance%min=0

      ! if (convergence /= 1.0D-8) convergence = 1.0_8
      ! if (convergence_int /= 1.0D-8) convergence_int = 1.0_8

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(CCSDT1T2)) deallocate(CCSDT1T2)
      ! allocate(CCSDT1T2(num_species))

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(CCSDinit)) deallocate(CCSDinit)
      ! allocate(CCSDinit(num_species))

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(CCSDloop)) deallocate(CCSDloop)
      ! allocate(CCSDloop(num_species))

      ! ! allocate array for results of inter-species in CCSD
      ! if (allocated(CCSDinter)) deallocate(CCSDinter)
      ! allocate(CCSDinter(num_inter*2))

      ! ! allocate array for results of intra-species in CCSD
      ! if (allocated(CCSD_instance%convergence_same)) deallocate(CCSD_instance%convergence_same)
      ! allocate(CCSD_instance%convergence_same(num_species))

      ! ! allocate array for results of intra-species in CCSD
      ! if (allocated(CCSD_instance%convergence_same_amp)) deallocate(CCSD_instance%convergence_same_amp)
      ! allocate(CCSD_instance%convergence_same_amp(num_species))

      ! ! allocate array for results of intra-species in CCSD
      ! if (allocated(CCSD_instance%convergence_diff_amp)) deallocate(CCSD_instance%convergence_diff_amp)
      ! allocate(CCSD_instance%convergence_diff_amp(num_inter*2))

      ! ! allocate array for results of intra-species in CCSD
      ! if (allocated(CCSD_instance%convergence_diff)) deallocate(CCSD_instance%convergence_diff)
      ! allocate(CCSD_instance%convergence_diff(num_inter))

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(e_same_ccd)) deallocate(e_same_ccd)
      ! allocate(e_same_ccd(num_species))
      ! e_same_ccd=0.0_8

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(a_same_ccd)) deallocate(a_same_ccd)
      ! allocate(a_same_ccd(num_species))
      ! a_same_ccd=0.0_8

      ! if (allocated(a_diff_ccd)) deallocate(a_diff_ccd)
      ! allocate(a_diff_ccd(num_inter))
      ! a_same_ccd=0.0_8

      ! ! allocate array for many results by species in CCSD
      ! if (allocated(e_diff_ccd)) deallocate(e_diff_ccd)
      ! allocate(e_diff_ccd(num_inter))
      ! e_diff_ccd=0.0_8


      do i=counterID, num_species
        print*, "num_species: CCSD(T)_constructor(): ", num_species, i
        call CCSDPT_constructor(i)
      end do

      do i=counterID, num_species
        print*, "num_species: CCSD(T)_run(): ", num_species, i, finalID
        call CCSDPT_init(i)
        print*, "CCSD(T)_init(i): ", i
      end do

        
      ! if (times_i>0) then
          
      !   !search the principal species (speciesId) in all the options
      !   do j=1, times_i
      !     print*, "times_i", times_i
      !     print*, "counterID: ", counterID, "finalID: ", finalID
      !     !public to private            
      !     i_counterID(j) = CoupledCluster_instance%i_counterID(j)
      !     n_intersp(j) = CoupledCluster_instance%n_intersp(j)
      !     !Find the appropriate species in all the options
      !     ! if (i_counterID(j)==counterID) then
      !     !public variables used in the loop
      !     CCSD_instance%max=n_intersp(j)
      !     CCSD_instance%min=i_counterID(j)
      !     !for the order in call in spintm matrix
      !     aux_cont = n_sp + 1
      !     CCSD_instance%aux_cont = aux_cont
      !     !make all combinations with speciesId
      !     n_sp = n_sp + 1
      !     CCSD_instance%cont = n_sp

      !     call CCSD_constructor_inter(i_counterID(j), n_intersp(j))!, num_inter)
      !     call CCSD_constructor_inter(n_intersp(j), i_counterID(j))!, num_inter)
      !     call CCSD_init_inter(i_counterID(j), n_intersp(j))!, num_inter)
      !     call CCSD_init_inter(n_intersp(j), i_counterID(j))!, num_inter)

      !     num_inter = num_inter + 1
      !   end do

      ! end if

      do i=counterID, num_species
        call CCSDPT_same_species(i)
          ! e_same_ccd(i) = CoupledCluster_instance%CCSD_E_intra(i)
          ! a_same_ccd(i) = CoupledCluster_instance%CCSD_a_intra(i)
      end do
        ! if (times_i>0) then
        !   do i=1, times_i
        !     CCSD_instance%num_i = 1
        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont
        !     max = CCSD_instance%max
        !     min = CCSD_instance%min

        !       num_i = CCSD_instance%num_i
        !       print*, "num_i: ", num_i 
        !       print*, "e_cont: ", CCSD_instance%e_cont
        !       call CCSD_diff_species(i_counterID(i),n_intersp(i),e_diff_ccd(num_i),a_diff_ccd(num_i))
        !       e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
        !       a_diff_ccd(num_i) = CoupledCluster_instance%CCSD_a_inter(num_i)
           
        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont

        !       num_i = num_i+1
        !       print*, "num_i: ", num_i
        !       print*, "e_cont: ", CCSD_instance%e_cont
        !       call CCSD_diff_species(n_intersp(i),i_counterID(i),e_diff_ccd(num_i),a_diff_ccd(num_i))
        !       e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
        !       a_diff_ccd(num_i) = CoupledCluster_instance%CCSD_a_inter(num_i)

        !   end do

        !   do i=1, times_i
        !     CCSD_instance%num_i = 1
        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont
        !     max = CCSD_instance%max
        !     min = CCSD_instance%min
        !     num_i = CCSD_instance%num_i
        !     call CCSD_diff_species_energy(i_counterID(i),n_intersp(i))
        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont
        !     call CCSD_diff_species_energy(n_intersp(i),i_counterID(i))
        !   end do

        !   do i=1, times_i
        !     CCSD_instance%num_i = 1
        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont
        !     max = CCSD_instance%max
        !     min = CCSD_instance%min
        !     num_i = CCSD_instance%num_i
        !     print*, "CCSD_T2_AB()"
        !     call CCSD_T2_AB(i_counterID(i),n_intersp(i), num_i)

        !     CCSD_instance%cont = CCSD_instance%aux_cont
        !     CCSD_instance%e_cont = CCSD_instance%aux_cont
        !     call CCSD_T2_AB(n_intersp(i),i_counterID(i), num_i)            
        !   end do
          
        ! end if

        ! print*, "CoupledCluster_instance%CCSD_E_intra: ", CoupledCluster_instance%CCSD_E_intra
        ! print*, "e_same_ccd: ", CCSD_instance%convergence_same
        ! print*, "CoupledCluster_instance%CCSD_A_intra: ", CoupledCluster_instance%CCSD_A_intra
        ! print*, "CoupledCluster_instance%CCSD_E_inter: ", CoupledCluster_instance%CCSD_E_inter
        ! print*, "e_diff_ccd: ", CCSD_instance%convergence_diff
        ! print*, "CoupledCluster_instance%CCSD_A_inter: ", CoupledCluster_instance%CCSD_A_inter


        ! convergence = CCSD_instance%convergence_same(1) !test

        ! For amplitudes: CCSD_instance%convergence_same_amp and CCSD_instance%convergence_same_amp

        ! For energy: CCSD_instance%convergence_same_amp and CCSD_instance%convergence_diff_amp
        ! if (stoped>=3) then
        !   if (convergence >= 1.0D-6) then
        !     convergence = CCSD_instance%convergence_same_amp(1) !test
        !   else 
        !     convergence = &!sum(CCSD_instance%convergence_same,dim=1) + & 
        !        sum(CCSD_instance%convergence_diff_amp,dim=1) !original
        !   end if
        ! end if

        ! print*, "Temporary convergence: ", convergence
        ! print*, "Temporary total ccsd-apmo energy: ", sum(e_same_ccd, dim=1) + &
        !   sum(e_diff_ccd, dim=1)


      do i=counterID, num_species
        call CCSDPT_show(i)
        ! intra = intra + CoupledCluster_instance%CCSD_E_intra(i)
      end do

      ! CoupledCluster_instance%CCSD_once_Energy = sum(CoupledCluster_instance%CCSD_E_intra,dim=1)
      ! CoupledCluster_instance%CCSD_twice_Energy = sum(CoupledCluster_instance%CCSD_E_inter,dim=1)
      ! CoupledCluster_instance%CCSD_total_Energy = CoupledCluster_instance%CCSD_once_Energy + &
      !   CoupledCluster_instance%CCSD_twice_Energy

      ! ! CoupledCluster_instance%CCSD_ones_Energy = intra

      ! print*, "Total correction same species CCSD energy: ", CoupledCluster_instance%CCSD_once_Energy
      ! print*, "Total correction different species CCSD energy: ", CoupledCluster_instance%CCSD_twice_Energy
      ! print*, "Total correction CCSD energy: ", CoupledCluster_instance%CCSD_total_Energy
      ! print*, "Total CCSD energy: ", CoupledCluster_instance%CCSD_total_Energy &
      !   + CoupledCluster_instance%HF_energy
      ! print*, "CCSDPT_show()"  
  end subroutine CCSDPT_run

  !>
  ! @brief Show results of Coupled Cluster
  !        call any energy that you want
  ! @author CAOM
  subroutine CCSDPT_show(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      print*, "INFORMATION OF CCSD(T) METHOD"
      print*, "speciesId: **", speciesId, "**"
      print*, "E_T[4] Energy of species " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "E_TS[5] Energy of different species " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "CCSD + E_T[4] = CCSD[T] correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "CCSD + E_TS[5] = CCSD(T) correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId) + CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "Total CCSD[T] energy: ", CoupledCluster_instance%HF_energy + CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "Total CCSD(T) energy: ", CoupledCluster_instance%HF_energy + CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId) + CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)

  end subroutine CCSDPT_show

end module CCSDPT_