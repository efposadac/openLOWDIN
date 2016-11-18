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
  use Tensor_
  implicit none

  type, public :: CCSD
      
      
      ! real(8), allocatable :: Tssame(:,:)
      ! real(8), allocatable :: Tdsame(:,:,:,:)
      ! real(8), allocatable :: tau(:,:,:,:)
      ! real(8), allocatable :: ttau(:,:,:,:)
      real(8) :: sum

      logical :: isInstanced

  end type CCSD

  type, public :: CCSDiter
      

      real(8), allocatable :: Dai(:,:)
      real(8), allocatable :: Fac(:,:)
      real(8), allocatable :: Fki(:,:)
      real(8), allocatable :: Fkc_aba(:,:)
      real(8), allocatable :: Fkc_aa(:,:)
      real(8), allocatable :: Wklij(:,:,:,:)
      real(8), allocatable :: Wabcd(:,:,:,:)
      real(8), allocatable :: Wkbcj(:,:,:,:)
      real(8), allocatable :: Tai(:,:)
      real(8), allocatable :: Tabij(:,:,:,:)

  end type CCSDiter

  type(CCSD), public :: CCSD_instance
  type(CCSDiter), public :: CCSDinit
  type(CCSDiter), public :: CCSDloop
  type(CCSDiter), public, allocatable :: CCSDT1T2(:)


contains

  !>
  ! @brief Constructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSD_constructor(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nop
      
      ! call CoupledCluster_pairing_function(1,2)
      
      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSD_constructor: noc=", noc, "nop=", nop
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! Denominator in T1 D^{a}_{i}
      if (allocated(CCSDinit%Dai)) deallocate (CCSDinit%Dai)
      allocate(CCSDinit%Dai(noc,noc))
      CCSDinit%Dai(:,:) = 0.0_8

      ! t^{a}_{i} amplitude for single excitation
      if (allocated(Allspecies(speciesId)%Tssame)) deallocate(Allspecies(speciesId)%Tssame)
      allocate(Allspecies(speciesId)%Tssame(noc-nop,nop))
      Allspecies(speciesId)%Tssame(:,:) = 0.0_8

      ! t^{ab}_{ij} amplitude for double excitation for same species
      if (allocated(Allspecies(speciesId)%Tdsame)) deallocate(Allspecies(speciesId)%Tdsame)
      allocate(Allspecies(speciesId)%Tdsame(noc-nop,noc-nop,nop,nop))
      Allspecies(speciesId)%Tdsame(:,:,:,:) = 0.0_8

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tilde{\tau}
      if (allocated(Allspecies(speciesId)%ttau)) deallocate (Allspecies(speciesId)%ttau)
      allocate(Allspecies(speciesId)%ttau(noc-nop,noc-nop,nop,nop))
      Allspecies(speciesId)%ttau(:,:,:,:) = 0.0_8

      ! \tau
      if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)
      allocate(Allspecies(speciesId)%tau(noc-nop,noc-nop,nop,nop))
      Allspecies(speciesId)%tau(:,:,:,:) = 0.0_8


      CCSD_instance%isInstanced = .true.
      
  end subroutine CCSD_constructor

  !>
  ! @brief Destructor of the class
  ! @author CAOM
  subroutine CCSD_destructor()
      implicit none

      ! if (allocated(CCSD_instance%Dai)) deallocate (CCSD_instance%Dai)
      ! if (allocated(Allspecies(speciesId)%Tssame)) deallocate (Allspecies(speciesId)%Tssame)
      ! if (allocated(Allspecies(speciesId)%Tdsame)) deallocate (Allspecies(speciesId)%Tdsame)
      ! if (allocated(Allspecies(speciesId)%ttau)) deallocate (Allspecies(speciesId)%ttau)
      ! if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)

      CCSD_instance%isInstanced = .false.
      
  end subroutine CCSD_destructor

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information
  ! @author CAOM
  subroutine CCSD_init(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nocs, nop, nops

      integer :: a, b, i, j

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tilde{\tau}^{ab}_{ij} = t^{ab}_{ij} + \frac{1}{2}(t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j})
      ! \tau^{ab}_{ij} = t^{ab}_{ij} + t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j}

      print*, "before loop"
      print*, Allspecies(speciesId)%Tdsame(1,1,1,1)
      print*, Allspecies(speciesId)%Tssame(1,1)
      print*, Allspecies(speciesId)%ttau(1,1,1,1)
      print*, Allspecies(speciesId)%tau(1,1,1,1)
      ! Allspecies(speciesId)%Tdsame(1,1,1,1) = Allspecies(speciesId)%tau(1,1,1,1)
      print*, "before loop"
      do a=nop+1, noc
        do b=nop+1, noc
          do i=1, nop
            do j=1, nop


              Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                +( (spints(speciesId)%valuesp(i,j,a,b))/( Allspecies(speciesId)%HF_fs%values(i,i)+ &
                  Allspecies(speciesId)%HF_fs%values(j,j) -Allspecies(speciesId)%HF_fs%values(a,a)- &
                    Allspecies(speciesId)%HF_fs%values(b,b) ) ) 
              
              Allspecies(speciesId)%ttau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + 0.5*( Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                  -Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j) )

              Allspecies(speciesId)%tau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                  -Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j)

              ! write(*,*) Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j), "Tdsame"
              ! , Allspecies(speciesId)%HF_fs%values(a,a), "HF_fs%values", spints(speciesId)%valuesp(i,j,a,b), "spints"
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

  end subroutine CCSD_init

  !>
  ! @brief Build a intermediates that will be used in Coupled Cluster loop
  ! @author CAOM
  subroutine CCSD_loop_constructor(speciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer noc, nocs, nop, nops
      
      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      
      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop_constructor: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

      !
      if (allocated(CCSDloop%Fac)) deallocate (CCSDloop%Fac)
      allocate(CCSDloop%Fac(noc-nop,noc-nop))
      CCSDloop%Fac=0.0_8

      !
      if (allocated(CCSDloop%Fki)) deallocate (CCSDloop%Fki)
      allocate(CCSDloop%Fki(nop,nop))
      CCSDloop%Fki=0.0_8
      
      !
      if (allocated(CCSDloop%Fkc_aa)) deallocate (CCSDloop%Fkc_aa)
      allocate(CCSDloop%Fkc_aa(nop,noc-nop))
      CCSDloop%Fkc_aa=0.0_8

      !
      if (allocated(CCSDloop%Wklij)) deallocate (CCSDloop%Wklij)
      allocate(CCSDloop%Wklij(nop,nop,nop,nop))
      CCSDloop%Wklij=0.0_8

      !
      if (allocated(CCSDloop%Wabcd)) deallocate (CCSDloop%Wabcd)
      allocate(CCSDloop%Wabcd(noc-nop,noc-nop,noc-nop,noc-nop))
      CCSDloop%Wabcd=0.0_8

      !
      if (allocated(CCSDloop%Wkbcj)) deallocate (CCSDloop%Wkbcj)
      allocate(CCSDloop%Wkbcj(nop,noc-nop,noc-nop,nop))
      CCSDloop%Wkbcj=0.0_8
      
  end subroutine CCSD_loop_constructor

  !>
  ! @brief Build T1 and T2 amplitude equations that will be information of intermediates
  ! @author CAOM
  subroutine CCSD_T1T2_constructor(speciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer noc, nocs, nop, nops

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_constructor: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops     
      
      !
      if (allocated(CCSDT1T2(speciesId)%Tai)) deallocate (CCSDT1T2(speciesId)%Tai)
      allocate(CCSDT1T2(speciesId)%Tai(noc-nop,nop))
      CCSDT1T2(speciesId)%Tai=0.0_8
      
      !
      if (allocated(CCSDT1T2(speciesId)%Tabij)) deallocate (CCSDT1T2(speciesId)%Tabij)
      allocate(CCSDT1T2(speciesId)%Tabij(noc-nop,noc-nop,nop,nop))
      CCSDT1T2(speciesId)%Tabij=0.0_8

  end subroutine CCSD_T1T2_constructor

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information
  !         for inter-species
  ! @author CAOM
  subroutine CCSD_constructor_inter(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId

      integer noc, nop, nocs, nops
      
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSD_constructor: noc=", noc, "nop=", nop, "nocs=", nocs, "nops=", nops
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! lowercase = alpha species | uppercase = beta species

      ! t^{aB}_{iJ} amplitude for double excitation for different species
      if (allocated(Allinterspecies(speciesId)%Tdsame)) deallocate(Allinterspecies(speciesId)%Tdsame)
      allocate(Allinterspecies(speciesId)%Tdsame(noc-nop,nocs-nops,nop,nops))
      Allinterspecies(speciesId)%Tdsame(:,:,:,:) = 0.0_8

      ! Effective two-particle excitation operators \tilde{\tau} and \tau for different species:

      ! \tilde{\tau}
      if (allocated(Allinterspecies(speciesId)%ttau)) deallocate (Allinterspecies(speciesId)%ttau)
      allocate(Allinterspecies(speciesId)%ttau(noc-nop,noc-nop,nop,nop))
      Allinterspecies(speciesId)%ttau(:,:,:,:) = 0.0_8

      ! \tau
      if (allocated(Allinterspecies(speciesId)%tau)) deallocate (Allinterspecies(speciesId)%tau)
      allocate(Allinterspecies(speciesId)%tau(noc-nop,noc-nop,nop,nop))
      Allinterspecies(speciesId)%tau(:,:,:,:) = 0.0_8

  end subroutine CCSD_constructor_inter

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information for inter-species
  ! @author CAOM
  subroutine CCSD_init_inter(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId

      integer noc, nocs, nop, nops

      integer :: a, b, i, j

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      ! lowercase = alpha species | uppercase = beta species

      ! Effective two-particle excitation operators \tilde{\tau} and \tau for inter-species:

      ! \tilde{\tau}^{aB}_{iJ} = t^{aB}_{iJ} + \frac{1}{2}(t^{a}_{i}t^{B}_{J} - t^{B}_{i}t^{a}_{J})
      ! \tau^{aB}_{iJ} = t^{aB}_{iJ} + t^{a}_{i}t^{B}_{J} - t^{B}_{i}t^{a}_{J}

      print*, "before loop inter"
      print*, Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      print*, Allinterspecies(speciesId)%ttau(1,1,1,1)
      print*, Allinterspecies(speciesId)%tau(1,1,1,1)
      print*, "before loop inter"

      do a=nop+1, noc
        do b=nops+1, nocs
          do i=1, nop
            do j=1, nops


              Allinterspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) = Allinterspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                +( (spintm(speciesId)%valuesp(i,j,a,b))/( Allinterspecies(speciesId)%HF_fs%values(i,i)+ &
                  Allinterspecies(speciesId)%HF_fs%values(j,j) -Allinterspecies(speciesId)%HF_fs%values(a,a)- &
                    Allinterspecies(speciesId)%HF_fs%values(b,b) ) ) 
              
              Allinterspecies(speciesId)%ttau(a-nop,b-nop,i,j) = Allinterspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + 0.5*( Allinterspecies(speciesId)%Tssame(a-nop,i)*Allinterspecies(speciesId)%Tssame(b-nop,j) &
                  -Allinterspecies(speciesId)%Tssame(b-nop,i)*Allinterspecies(speciesId)%Tssame(a-nop,j) )

              Allinterspecies(speciesId)%tau(a-nop,b-nop,i,j) = Allinterspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + Allinterspecies(speciesId)%Tssame(a-nop,i)*Allinterspecies(speciesId)%Tssame(b-nop,j) &
                  -Allinterspecies(speciesId)%Tssame(b-nop,i)*Allinterspecies(speciesId)%Tssame(a-nop,j)

 ! (auxspinints(ii,iii,aa,aaa)/(Fs%values(ii,ii)+otherFs%values(iii,iii)-Fs%values(aa,aa)-otherFs%values(aaa,aaa)))
 !                         auxtaus(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) + 0.5*(CoupledClusterold_instance%Tstest(aa-nop,ii)*otherTs(aaa-nops,iii)) ! Eliminated crossed t1 terms
 !                         auxtau(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) + CoupledClusterold_instance%Tstest(aa-nop,ii)*otherTs(aaa-nops,iii) ! same here


              ! write(*,*) Allinterspecies(speciesId)%Tdsame(a-nop,b-nop,i,j), "Tdsame"
              ! , Allinterspecies(speciesId)%HF_fs%values(a,a), "HF_fs%values", spintm(speciesId)%valuesp(i,j,a,b), "spintm"
            end do
          end do
        end do
      end do
      print*, "CCSD_init_inter():"
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

  end subroutine CCSD_init_inter

  !>
  ! @brief Calculate F intermediates for intra-species
  ! @author CAOM
  subroutine F_onespecies_intermediates(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer :: noc, nop
      integer :: a, b, e, i, j, f, m, n

      noc = Allspecies(speciesId)%noc
      nop = Allspecies(speciesId)%nop

      print*,"noc: ", noc, "nop: ", nop
      print*, "speciesId: ", speciesId
      if (speciesId>2) stop "inside loop"

      ! CCSDloop%Fac
      do a=nop+1, noc
        do e=nop+1, noc

          CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) & 
            +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e)
            
          do m=1, nop
            CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
              + (-0.5*Allspecies(speciesId)%HF_fs%values(m,e)*Allspecies(speciesId)%Tssame(a-nop,m))
            do f=nop+1, noc
              CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
                + Allspecies(speciesId)%Tssame(f-nop,m)*spints(speciesId)%valuesp(m,a,f,e)
              do n=1, nop
                CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
                  + (-0.5*Allspecies(speciesId)%ttau(a-nop,f-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f))
                ! write(*,*) a,e,CCSDloop%Fac(a,e)
                ! write(*,*) a-nop,f-nop,m,n,spints(speciesId)%valuesp(m,n,e,f)*(-0.5*Allspecies(speciesId)%ttau(a-nop,f-nop,m,n))
              end do
            end do
          end do
        end do
      end do

      ! CCSDloop%Fki(m,i)
      do m=1, nop
        do i=1, nop
             
          CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
            + (1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i)
            
          do e=nop+1, noc
            CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
               + 0.5*Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%HF_fs%values(m,e)
            do n=1, nop
              CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
                + Allspecies(speciesId)%Tssame(e-nop,n)*spints(speciesId)%valuesp(m,n,i,e)
              do f=nop+1, noc
                CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
                  + 0.5*Allspecies(speciesId)%ttau(e-nop,f-nop,i,n)*spints(speciesId)%valuesp(m,n,e,f)
                ! write(*,*) a,e,CCSDloop%Fki(a,e)
                ! write(*,*) a,e,e-nop,f-nop,0.5*Allspecies(speciesId)%ttau(e-nop,f-nop,i,n)*spints(speciesId)%valuesp(m,n,e,f)
              end do
            end do
          end do
        end do
      end do

      ! CCSDloop%Fkc_aa
      do m=1, nop
        do e=nop+1, noc

          CCSDloop%Fkc_aa(m,e-nop) = CCSDloop%Fkc_aa(m,e-nop) &
            + Allspecies(speciesId)%HF_fs%values(m,e)
          do n=1, nop
            do f=nop+1, noc
              CCSDloop%Fkc_aa(m,e-nop) = CCSDloop%Fkc_aa(m,e-nop) &
                + Allspecies(speciesId)%Tssame(f-nop,n)*spints(speciesId)%valuesp(m,n,e,f)
              ! write(*,*) a,e,CCSDloop%Fkc_aa(a,e)
            end do
          end do
        end do
      end do

  end subroutine

  !>
  ! @brief Calculate W intermediates for intra-species
  ! @author CAOM
  subroutine W_onespecies_intermediates(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer :: noc, nop
      integer :: a, b, e, i, j, f, m, n

      noc = Allspecies(speciesId)%noc
      nop = Allspecies(speciesId)%nop

      ! CCSDloop%Wklij
      do m=1, nop
        do n=1, nop
          do i=1, nop
            do j=1, nop

              CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                + spints(speciesId)%valuesp(m,n,i,j)
              do e=nop+1, noc
                CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                  + (Allspecies(speciesId)%Tssame(e-nop,j)*spints(speciesId)%valuesp(m,n,i,e) &
                      -Allspecies(speciesId)%Tssame(e-nop,i)*spints(speciesId)%valuesp(m,n,j,e))
                do f=nop+1, noc
                  CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                    + 0.25*Allspecies(speciesId)%tau(e-nop,f-nop,i,j)*spints(speciesId)%valuesp(m,n,e,f)
                  ! write(*,*) m,n,i,j,CCSDloop%Wklij(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do

      !CCSDloop%Wabcd
      do a=nop+1, noc
        do b=nop+1, noc
          do e=nop+1, noc
            do f=nop+1, noc

              CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                + spints(speciesId)%valuesp(a,b,e,f)
              do m=1, nop
                CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                  + (-Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(a,m,e,f) &
                      +Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(b,m,e,f))
                do n=1, nop
                  CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                    + 0.25*Allspecies(speciesId)%tau(a-nop,b-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f)
                    ! write(*,*) m,n,i,j,CCSDloop%Wabcd(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do

      !CCSDloop%Wkbcj
      do m=1, nop
        do b=nop+1, noc
          do e=nop+1, noc
            do j=1, nop
              CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                + spints(speciesId)%valuesp(m,b,e,j)
              do f=nop+1, noc
                CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                 + Allspecies(speciesId)%Tssame(f-nop,j)*spints(speciesId)%valuesp(m,b,e,f)
              end do
              do n=1, nop
                CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                  - Allspecies(speciesId)%Tssame(b-nop,n)*spints(speciesId)%valuesp(m,n,e,j)
                do f=nop+1, noc
                  CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                    - ((0.5*Allspecies(speciesId)%Tdsame(f-nop,b-nop,j,n) &
                        + Allspecies(speciesId)%Tssame(f-nop,j)*Allspecies(speciesId)%Tssame(b-nop,n))*spints(speciesId)%valuesp(m,n,e,f))
                  ! write(*,*) m,n,i,j,CCSDloop%Wkbcj(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do

  end subroutine

  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSD_loop(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nocs, nop, nops
      integer num_species
      
      integer :: a, b, e, i, j, f, m, n
      real(8) :: ccsdE=0.0_8
      real(8) :: convergence = 1.0_8
      real(8) :: prev_ccsdE
      real(8) :: tmp_ccsdE

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (ccsdE /= 0.0_8) ccsdE = 0.0_8

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      num_species = CoupledCluster_instance%num_species
      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

      print*, "T1T2_constructor", convergence
      
      do while (convergence >= 1.0D-8)

        prev_ccsdE = ccsdE

        print*, "speciesId CCSD_loop: ", speciesId

        call CCSD_T1T2_constructor(speciesId)
        call CCSD_loop_constructor(speciesId)

        !intermediates loop for:
        !singles excitations
        call F_onespecies_intermediates(speciesId)
        !doubles excitations
        call W_onespecies_intermediates(speciesId)

 
        if (speciesId>1) then
          do j=1, num_species
            if (i /= j) then
              call CCSD_constructor_inter(speciesId, j)
              !call F_inter()
              !call W_inter()
            end if
          end do
        end if


        ! Resolve CCSD equation of energy

        tmp_ccsdE=0.0_8
        do i=1, nop
          do a=nop+1, noc
            tmp_ccsdE = tmp_ccsdE + Allspecies(speciesId)%HF_fs%values(i,a)*Allspecies(speciesId)%Tssame(a-nop,i)
            do j=1, nop
              do b=nop+1, noc
                tmp_ccsdE = tmp_ccsdE + (0.25*spints(speciesId)%valuesp(i,j,a,b) *Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                  + 0.5*spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j))
              end do
            end do
          end do
        end do
        ccsdE = tmp_ccsdE

        if (speciesId==2) print*, "tmp_ccsdE: ", tmp_ccsdE

        !change in values for the loop
        convergence = abs( ccsdE - prev_ccsdE )

        write (*,*) ccsdE, "CCSD Energy ", prev_ccsdE, "previous Energy" 
        write (*,*) convergence, "Convergence " 
        ! if ((convergence > 10) .and. speciesId>1) stop "test"
        ! Resolve T1 and T2 amplitude equations
        call CCSD_T1T2(speciesId)

        if (convergence > 100) then 
          stop "Error: There are not convergence. The differences between energies is more than 100 eV"
        end if

      end do

      CoupledCluster_instance%CCSD_E_intra(speciesId) = ccsdE
      
  end subroutine CCSD_loop

  !>
  ! @brief Calculate T1 and T2 energy equations for intra-species
  ! @author CAOM
  subroutine CCSD_T1T2(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nocs, nop, nops
      integer :: a, b, e, f, i, j, m, n

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      ! T^{a}_{i}D^{a}_{i} = ...
      do a=nop+1, noc
        do i=1, nop
          CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) + Allspecies(speciesId)%HF_fs%values(i,a)
          do e=nop+1, noc
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop%Fac(a-nop,e-nop)
          end do
          do m=1, nop
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + (-Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop%Fki(m,i))
            do e=nop+1, noc
              CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop%Fkc_aa(m,e-nop)
              do f=nop+1, noc
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + (-0.5*Allspecies(speciesId)%Tdsame(e-nop,f-nop,i,m)*spints(speciesId)%valuesp(m,a,e,f))
              end do
              do n=1, nop
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,n)*spints(speciesId)%valuesp(n,m,e,i))
              end do
            end do
          end do
          do n=1,nop
            do f=nop+1, noc
              CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                + (-Allspecies(speciesId)%Tssame(f-nop,n)*spints(speciesId)%valuesp(n,a,i,f))
            end do
          end do
          CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit%Dai(a,i)
          Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
          ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
        end do
      end do


      ! T^{ab}_{ij}D^{ab}_{ij} = ...
      do a=nop+1, noc
         do b=nop+1, noc
            do i=1, nop
               do j=1, nop
                  CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                    + spints(speciesId)%valuesp(i,j,a,b) !A
                  ! 1er ciclo
                  do e=nop+1, noc
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,j)*CCSDloop%Fac(b-nop,e-nop) &
                        -Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,j)*CCSDloop%Fac(a-nop,e-nop)) !B
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (Allspecies(speciesId)%Tssame(e-nop,i)*spints(speciesId)%valuesp(a,b,e,j) &
                        -Allspecies(speciesId)%Tssame(e-nop,j)*spints(speciesId)%valuesp(a,b,e,i)) !G
                     do f=nop+1, noc
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*Allspecies(speciesId)%tau(e-nop,f-nop,i,j)*CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) !D
                     end do
                     do m=1, nop
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,j)*Allspecies(speciesId)%Tssame(b-nop,m)*CCSDloop%Fkc_aa(m,e-nop) &
                              +0.5*Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,j)*Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop%Fkc_aa(m,e-nop)) !B'
                     end do
                  end do
                  ! 2do ciclo
                  do m=1, nop
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (-Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,m)*CCSDloop%Fki(m,j) &
                          +Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,m)*CCSDloop%Fki(m,i)) !C
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (-Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,i,j) &
                          +Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,i,j)) !H
                     do n=1, nop
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*Allspecies(speciesId)%tau(a-nop,b-nop,m,n)*CCSDloop%Wklij(m,n,i,j) !E
                     end do
                     do e=nop+1, noc
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,m)*Allspecies(speciesId)%Tssame(e-nop,j)*CCSDloop%Fkc_aa(m,e-nop) &
                              +0.5*Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,m)*Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop%Fkc_aa(m,e-nop)) !C'
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                            - Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,e,j) &
                              -Allspecies(speciesId)%Tdsame(a-nop,e-nop,j,m)*CCSDloop%Wkbcj(m,b-nop,e-nop,i) &
                                + Allspecies(speciesId)%Tssame(e-nop,j)*Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,e,i) &
                                  -Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,m)*CCSDloop%Wkbcj(m,a-nop,e-nop,j) &
                                    - Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,e,j) &
                                      + Allspecies(speciesId)%Tdsame(b-nop,e-nop,j,m)*CCSDloop%Wkbcj(m,a-nop,e-nop,i) &
                                       - Allspecies(speciesId)%Tssame(e-nop,j)*Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,e,i) !F
                     end do
                  end do
                  ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
                  CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                    /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(speciesId)%HF_fs%values(j,j)-Allspecies(speciesId)%HF_fs%values(a,a) &
                       -Allspecies(speciesId)%HF_fs%values(b,b))
                  
                  Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j)

                  Allspecies(speciesId)%ttau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                    + 0.5*(Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j))
                  Allspecies(speciesId)%tau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                    + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j)
               ! write(*,*) a,b,i,j,Allspecies(speciesId)%Tdsame(a,b,i,j),CCSDT1T2(speciesId)%Tabij(a,b,i,j)
               end do
            end do
         end do
      end do
      
  end subroutine CCSD_T1T2

  !>
  ! @brief Manager of Coupled Cluster maths
  ! @author CAOM
  subroutine CCSD_run()
      implicit none

      integer :: num_species
      integer :: num_inter
      integer :: i
      real(8) :: intra=0
      
      num_species = CoupledCluster_instance%num_species
      num_inter = CoupledCluster_instance%num_intersp

      ! allocate array for many results by species in CCSD
      if (allocated(CCSDT1T2)) deallocate(CCSDT1T2)
      allocate(CCSDT1T2(num_species))

      do i=1, num_species
        
        call CCSD_constructor(i)
        print*, "num_species: CCSD_run(): ", num_species, i
        call CCSD_init(i)
        print*, "num_species: CCSD_loop(): ", num_species, i
        call CCSD_loop(i)
        call CCSD_show(i)

        intra = intra + CoupledCluster_instance%CCSD_E_intra(i)

      end do

      CoupledCluster_instance%CCSD_ones_Energy = intra

      print*, "Total CCSD energy: ", CoupledCluster_instance%CCSD_ones_Energy &
        + CoupledCluster_instance%HF_energy
      ! call CCSD_show()
      print*, "CCSD_show()"
      
  end subroutine CCSD_run

  !>
  ! @brief Show results of Coupled Cluster
  !        call any energy that you want
  ! @author CAOM
  subroutine CCSD_show(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      print*, "INFORMATION IN CCSD_constructor() HF_energy: ", CoupledCluster_instance%HF_energy
      print*, "INFORMATION IN CCSD_constructor() MP2_energy: ", CoupledCluster_instance%MP2_EnergyCorr
      CCSD_instance%sum = CoupledCluster_instance%HF_energy + CoupledCluster_instance%MP2_EnergyCorr
      print*, "INFORMATION IN CCSD_constructor() TotalMP2_Energy: ", CCSD_instance%sum
      print*, "INFORMATION IN CCSD_constructor() Td: ", Allspecies(speciesId)%Tdsame(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() Dai: ", CCSDinit%Dai(1,3)
      print*, "INFORMATION IN CCSD_constructor() ttau: ",Allspecies(speciesId)%ttau(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() tau: ",Allspecies(speciesId)%tau(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() ccsdE: ", CoupledCluster_instance%CCSD_E_intra(1)
      print*, "INFORMATION IN Tensor Tensor_index2: ", IndexMap_tensorR2ToVector(3,4,2)
      

      !if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)

      !call CCSD_destructor()
      print*, "CCSD_show()x2"

  end subroutine CCSD_show
end module CCSD_