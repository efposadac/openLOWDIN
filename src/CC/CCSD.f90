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
      
      integer :: max
      integer :: min
      integer :: cont, aux_cont
      integer :: num_i, e_cont
      real(8) :: suma
      real(8), allocatable :: convergence_same(:)
      real(8), allocatable :: convergence_diff(:)

      logical :: isInstanced

  end type CCSD

  type, public :: CCSDiter
      

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
      real(8), allocatable :: Tabij(:,:,:,:)
      real(8), allocatable :: Tabij_AB(:,:,:,:)

  end type CCSDiter

  type(CCSD), public :: CCSD_instance
  type(CCSDiter), public, allocatable :: CCSDinit(:)
  type(CCSDiter), public, allocatable :: CCSDloop(:)
  type(CCSDiter), public, allocatable :: CCSDT1T2(:)
  type(CCSDiter), public, allocatable :: CCSDinter(:)


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
      if (allocated(CCSDinit(speciesId)%Dai)) deallocate (CCSDinit(speciesId)%Dai)
      allocate(CCSDinit(speciesId)%Dai(noc,noc))
      CCSDinit(speciesId)%Dai(:,:) = 0.0_8

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
      if (nop>=2) then
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
      end if
      print*, "CCSD_init(", speciesId, ")"
      ! Denominator D^{a}_{i}
      do a=nop+1, noc
         do i=1, nop
            CCSDinit(speciesId)%Dai(a,i) = Allspecies(speciesId)%HF_fs%values(i,i) - Allspecies(speciesId)%HF_fs%values(a,a)
            write(*,*) a,i,CCSDinit(speciesId)%Dai(a,i)
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
      
      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop_constructor: noc=", noc, "nop=", nop

      !
      if (allocated(CCSDloop(speciesId)%Fac)) deallocate (CCSDloop(speciesId)%Fac)
      allocate(CCSDloop(speciesId)%Fac(noc-nop,noc-nop))
      CCSDloop(speciesId)%Fac=0.0_8

      !
      if (allocated(CCSDloop(speciesId)%Fki)) deallocate (CCSDloop(speciesId)%Fki)
      allocate(CCSDloop(speciesId)%Fki(nop,nop))
      CCSDloop(speciesId)%Fki=0.0_8
      
      !
      if (allocated(CCSDloop(speciesId)%Fkc_aa)) deallocate (CCSDloop(speciesId)%Fkc_aa)
      allocate(CCSDloop(speciesId)%Fkc_aa(nop,noc-nop))
      CCSDloop(speciesId)%Fkc_aa=0.0_8

      !
      if (allocated(CCSDloop(speciesId)%Wklij)) deallocate (CCSDloop(speciesId)%Wklij)
      allocate(CCSDloop(speciesId)%Wklij(nop,nop,nop,nop))
      CCSDloop(speciesId)%Wklij=0.0_8

      !
      if (allocated(CCSDloop(speciesId)%Wabcd)) deallocate (CCSDloop(speciesId)%Wabcd)
      allocate(CCSDloop(speciesId)%Wabcd(noc-nop,noc-nop,noc-nop,noc-nop))
      CCSDloop(speciesId)%Wabcd=0.0_8

      !
      if (allocated(CCSDloop(speciesId)%Wkbcj)) deallocate (CCSDloop(speciesId)%Wkbcj)
      allocate(CCSDloop(speciesId)%Wkbcj(nop,noc-nop,noc-nop,nop))
      CCSDloop(speciesId)%Wkbcj=0.0_8
      
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
      ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_constructor: noc=", noc, "nop=", nop     
      
      !
      if (allocated(CCSDT1T2(speciesId)%Tai)) deallocate (CCSDT1T2(speciesId)%Tai)
      allocate(CCSDT1T2(speciesId)%Tai(noc-nop,nop))
      CCSDT1T2(speciesId)%Tai=0.0_8
      
      !
      if (allocated(CCSDT1T2(speciesId)%Tabij)) deallocate (CCSDT1T2(speciesId)%Tabij)
      allocate(CCSDT1T2(speciesId)%Tabij(noc-nop,noc-nop,nop,nop))
      CCSDT1T2(speciesId)%Tabij=0.0_8
      print*, " CCSD_T1T2_constructor(: ", nop, noc, speciesId

  end subroutine CCSD_T1T2_constructor

  !>
  ! @brief Build T1 and T2 amplitude equations that will be information of intermediates
  ! @author CAOM
  subroutine CCSD_T2AB_constructor(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter
      integer noc, nocs, nop, nops

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T2AB_constructor: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops     

      !
      if (allocated(CCSDT1T2(speciesId)%Tabij_AB)) deallocate (CCSDT1T2(speciesId)%Tabij_AB)
      allocate(CCSDT1T2(speciesId)%Tabij_AB(noc-nop,nocs-nops,nop,nops))
      CCSDT1T2(speciesId)%Tabij_AB=0.0_8
      print*, "nops 34 "

  end subroutine CCSD_T2AB_constructor

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information
  !         for inter-species
  ! @author CAOM
  subroutine CCSD_constructor_inter(speciesId, OtherspeciesId)!, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      ! integer, intent(in) :: num_inter

      integer noc, nop, nocs, nops
      
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSD_constructor_inter: noc=", noc, "nop=", nop, "nocs=", nocs, "nops=", nops
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! lowercase = alpha species | uppercase = beta species

      ! t^{aB}_{iJ} amplitude for double excitation for different species
      if (allocated(Allinterspecies(speciesId)%Tdsame)) deallocate(Allinterspecies(speciesId)%Tdsame)
      allocate(Allinterspecies(speciesId)%Tdsame(noc-nop,nocs-nops,nop,nops))
      Allinterspecies(speciesId)%Tdsame(:,:,:,:) = 0.0_8

      ! Effective two-particle excitation operators \check{\tau}, \ddot{\tau} and \tilde{\tau} 
      !   for different species:

      !There are two kind of \check{\tau}}:

      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
        ! \check{\tau} alpha
        if (allocated(Allinterspecies(speciesId)%chtau_a)) deallocate (Allinterspecies(speciesId)%chtau_a)
        allocate(Allinterspecies(speciesId)%chtau_a(noc-nop,noc-nop,nocs-nops,nop,nop,nops))
        Allinterspecies(speciesId)%chtau_a(:,:,:,:,:,:) = 0.0_8
      end if

      if (nops>=2) then ! kind of interaction just for two or more particles of the another species
        ! \check{\tau} beta
        if (allocated(Allinterspecies(OtherspeciesId)%chtau_b)) deallocate (Allinterspecies(OtherspeciesId)%chtau_b)
        allocate(Allinterspecies(OtherspeciesId)%chtau_b(nocs-nops,nocs-nops,noc-nop,nops,nops,nop))
        Allinterspecies(OtherspeciesId)%chtau_b(:,:,:,:,:,:) = 0.0_8
      end if

      ! \ddot{\tau}
      if (allocated(Allinterspecies(speciesId)%tau)) deallocate (Allinterspecies(speciesId)%tau)
      allocate(Allinterspecies(speciesId)%tau(noc-nop,noc-nop,nop,nop))
      Allinterspecies(speciesId)%tau(:,:,:,:) = 0.0_8

      ! \tilde{\tau}
      if (allocated(Allinterspecies(speciesId)%intau)) deallocate (Allinterspecies(speciesId)%intau)
      allocate(Allinterspecies(speciesId)%intau(noc-nop,noc-nop,nop,nop))
      Allinterspecies(speciesId)%intau(:,:,:,:) = 0.0_8

      print*, "fin"

  end subroutine CCSD_constructor_inter

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information for inter-species
  ! @author CAOM
  subroutine CCSD_init_inter(speciesId, OtherspeciesId)!, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      ! integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops, n_sp, num_species

      integer :: a, b, c, i, ii, j, k
      integer :: aa, bb, jj, cc, kk
      integer :: p, q, r, s

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      ! n_sp = Tix2(speciesId, OtherspeciesId, num_species)
      ! number of transformed integrals matrix for speciesId and OtherspeciesId
      n_sp =CCSD_instance%cont

      ! lowercase = alpha species | uppercase = beta species

      ! Effective two-particle excitation operators \ddot{\tau}, \check{\tau} and \tilde{\tau} for inter-species:

      ! \check{\tau} _{ijk}^{abc}=t_{i}^{a}t_{j}^{b}t_{k}^{c}+\frac{1}{2}\left(t_{i}^{a}t_{jk}^{bc} &
      !   + t_{j}^{b}t_{ik}^{ac} + t_{k}^{c}t_{ji}^{ba} \right)

      ! \ddot{\tau} _{iJ}^{aB}=t_{iJ}^{aB} - t_{i}^{a}t_{J}^{B}

      ! \tilde{\tau} _{iJ}^{aB}=t_{iJ}^{aB} - 0.5*t_{i}^{a}t_{J}^{B}

      print*, "before loop inter"
      print*, Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_a(1,1,1,1,1,1)
      !print*, Allinterspecies(OtherspeciesId)%chtau_b(1,1,1,1,1,1)
      print*, Allinterspecies(speciesId)%tau(1,1,1,1)
      print*, Allinterspecies(speciesId)%intau(1,1,1,1)

      do a=nop+1, noc
        do bb=nops+1, nocs
          do i=1, nop
            do jj=1, nops

              if (speciesId<OtherspeciesId) then
                p=i
                q=jj
                r=a
                s=bb
              else
                p=jj
                q=i
                r=bb
                s=a
              end if
              Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                +( (spintm(n_sp)%valuesp(p,q,r,s))/( Allspecies(speciesId)%HF_fs%values(i,i)+ &
                  Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) -Allspecies(speciesId)%HF_fs%values(a,a)- &
                    Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) ) ) 

              ! \ddot{\tau} 
              Allinterspecies(speciesId)%tau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)

              ! \tilde{\tau} 
              Allinterspecies(speciesId)%intau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - 0.5*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)

              if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                do b=nop+1, noc
                  do cc=nops+1, nocs
                    do j=1, nop
                      do kk=1, nops
                        ! \check{\tau} triple excitation holy shit!!
                        Allinterspecies(speciesId)%chtau_a(a-nop,b-nop,cc-nops,i,j,kk) = Allspecies(speciesId)%Tssame(a-nop,i)* &
                          Allspecies(speciesId)%Tssame(b-nop,j)*Allspecies(OtherspeciesId)%Tssame(cc-nops,kk) &
                          + 0.5*( Allspecies(speciesId)%Tssame(a-nop,i)*Allinterspecies(speciesId)%Tdsame(b-nop,cc-nops,j,kk) &
                            + Allspecies(speciesId)%Tssame(b-nop,j)*Allinterspecies(speciesId)%Tdsame(a-nop,cc-nops,i,kk) &
                              + Allspecies(OtherspeciesId)%Tssame(cc-nops,kk)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) )
                      end do
                    end do
                  end do
                end do
              end if

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do aa=nops+1, nocs
                  do c=nop+1, noc
                    do ii=1, nops
                      do k=1, nop
                        !  \check{\tau} triple excitation holy shit!!
                        Allinterspecies(OtherspeciesId)%chtau_b(aa-nops,bb-nops,c-nop,ii,jj,k) = Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)* &
                          Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(speciesId)%Tssame(c-nop,k) &
                          + 0.5*( Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)*Allinterspecies(speciesId)%Tdsame(bb-nops,c-nop,jj,k) &
                            + Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allinterspecies(speciesId)%Tdsame(aa-nops,c-nop,ii,k) &
                              + Allspecies(speciesId)%Tssame(c-nop,k)*Allspecies(OtherspeciesId)%Tdsame(aa-nops,bb-nops,ii,jj) )
                      end do
                    end do
                  end do
                end do
              end if

              ! write(*,*) Allinterspecies(speciesId)%Tdsame(a-nop,bb-nop,i,jj), "Tdsame"
              ! , Allinterspecies(speciesId)%HF_fs%values(a,a), "HF_fs%values", spintm(n_sp)%valuesp(i,jj,a,bb), "spintm"
            end do
          end do
        end do
        print*, "nops 1 - 5"
      end do

      print*, Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_a(1,1,1,1,1,1)
      !print*, Allinterspecies(OtherspeciesId)%chtau_b(1,1,1,1,1,1)
      print*, Allinterspecies(speciesId)%tau(1,1,1,1)
      print*, Allinterspecies(speciesId)%intau(1,1,1,1)
      print*, "before before inter" 

  end subroutine CCSD_init_inter

  !>
  ! @brief Build a intermediates that will be used in Coupled Cluster loop
  ! @author CAOM
  subroutine CCSD_loop_constructor_inter(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops, ierr
      
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      
      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop_constructor_inter: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

      !
      if (allocated(CCSDinter(speciesId)%Fkc_aba)) deallocate (CCSDinter(speciesId)%Fkc_aba)
      allocate(CCSDinter(speciesId)%Fkc_aba(nop,noc-nop))
      CCSDinter(speciesId)%Fkc_aba=0.0_8
      
      !Same species just initialized for interspecies
      if (allocated(CCSDinter(OtherspeciesId)%Fkca_ab)) deallocate (CCSDinter(OtherspeciesId)%Fkca_ab)
      allocate(CCSDinter(OtherspeciesId)%Fkca_ab(nops,nocs-nops))
      CCSDinter(OtherspeciesId)%Fkca_ab=0.0_8

      !
      if (allocated(CCSDinter(speciesId)%Faca)) deallocate (CCSDinter(speciesId)%Faca)
      allocate(CCSDinter(speciesId)%Faca(noc-nop,noc-nop))
      CCSDinter(speciesId)%Faca=0.0_8

      !
      if (allocated(CCSDinter(speciesId)%Fkia)) deallocate (CCSDinter(speciesId)%Fkia)
      allocate(CCSDinter(speciesId)%Fkia(nop,nop))
      CCSDinter(speciesId)%Fkia=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Fbcb)) deallocate (CCSDinter(OtherspeciesId)%Fbcb)
      allocate(CCSDinter(OtherspeciesId)%Fbcb(nocs-nops,nocs-nops))
      CCSDinter(OtherspeciesId)%Fbcb=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Fkjb)) deallocate (CCSDinter(OtherspeciesId)%Fkjb)
      allocate(CCSDinter(OtherspeciesId)%Fkjb(nops,nops))
      CCSDinter(OtherspeciesId)%Fkjb=0.0_8

      !
      if (allocated(CCSDinter(speciesId)%Wkkcc)) deallocate (CCSDinter(speciesId)%Wkkcc)
      allocate(CCSDinter(speciesId)%Wkkcc(noc-nop,nops,nop,nocs-nops))
      CCSDinter(speciesId)%Wkkcc=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Waka)) deallocate (CCSDinter(speciesId)%Waka)
      allocate(CCSDinter(speciesId)%Waka(nop,nocs-nops,nop,nops))
      CCSDinter(speciesId)%Waka=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wcia)) deallocate (CCSDinter(speciesId)%Wcia)
      allocate(CCSDinter(speciesId)%Wcia(noc-nop,nocs-nops,noc-nop,nops))
      CCSDinter(speciesId)%Wcia=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wbkb)) deallocate (CCSDinter(speciesId)%Wbkb)
      allocate(CCSDinter(speciesId)%Wbkb(noc-nop,nops,nop,nops))
      CCSDinter(speciesId)%Wbkb=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wcjb)) deallocate (CCSDinter(speciesId)%Wcjb)
      allocate(CCSDinter(speciesId)%Wcjb(noc-nop,nocs-nops,nop,nocs-nops))
      CCSDinter(speciesId)%Wcjb=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wakic_a)) deallocate (CCSDinter(speciesId)%Wakic_a)
      allocate(CCSDinter(speciesId)%Wakic_a(noc-nop,nop,nop,noc-nop))
      CCSDinter(speciesId)%Wakic_a=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Wbkjc_b)) deallocate (CCSDinter(OtherspeciesId)%Wbkjc_b)
      allocate(CCSDinter(OtherspeciesId)%Wbkjc_b(nocs-nops,nops,nops,nocs-nops))
      CCSDinter(OtherspeciesId)%Wbkjc_b=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Wklcd_b)) deallocate (CCSDinter(OtherspeciesId)%Wklcd_b)
      allocate(CCSDinter(OtherspeciesId)%Wklcd_b(nops,nops,nocs-nops,nocs-nops))
      CCSDinter(OtherspeciesId)%Wklcd_b=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wklcd_a)) deallocate (CCSDinter(speciesId)%Wklcd_a)
      allocate(CCSDinter(speciesId)%Wklcd_a(nop,nop,noc-nop,noc-nop))!,stat = ierr)
      CCSDinter(speciesId)%Wklcd_a=0.0_8

      print*, "nops 14"
      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wakic)) deallocate (CCSDinter(speciesId)%Wakic)
      allocate(CCSDinter(speciesId)%Wakic(noc-nop,nops,nop,nocs-nops))
      CCSDinter(speciesId)%Wakic=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Wbkjc)) deallocate (CCSDinter(OtherspeciesId)%Wbkjc)
      allocate(CCSDinter(OtherspeciesId)%Wbkjc(nop,nocs-nops,noc-nop,nops))
      CCSDinter(OtherspeciesId)%Wbkjc=0.0_8

  end subroutine CCSD_loop_constructor_inter

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
      ! if (speciesId>2) stop "inside loop"

      ! CCSDloop(speciesId)%Fac
      do a=nop+1, noc
        do e=nop+1, noc

          CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) & 
            +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e)
            
          do m=1, nop
            CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
              + (-0.5*Allspecies(speciesId)%HF_fs%values(m,e)*Allspecies(speciesId)%Tssame(a-nop,m))
            
            if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
              do f=nop+1, noc
                CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                  + Allspecies(speciesId)%Tssame(f-nop,m)*spints(speciesId)%valuesp(m,a,f,e)
                do n=1, nop
                  CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                    + (-0.5*Allspecies(speciesId)%ttau(a-nop,f-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f))
                  ! write(*,*) a,e,CCSDloop(speciesId)%Fac(a,e)
                end do
              end do
            end if
          end do
        end do
      end do

      ! CCSDloop(speciesId)%Fki(m,i)
      do m=1, nop
        do i=1, nop
             
          CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
            + (1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i)
            
          do e=nop+1, noc
            CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
               + 0.5*Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%HF_fs%values(m,e)
            
            if (nop>=2) then ! kind of interaction just for two or more particles of the principal species  
              do n=1, nop
                CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                  + Allspecies(speciesId)%Tssame(e-nop,n)*spints(speciesId)%valuesp(m,n,i,e)
                do f=nop+1, noc
                  CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                    + 0.5*Allspecies(speciesId)%ttau(e-nop,f-nop,i,n)*spints(speciesId)%valuesp(m,n,e,f)
                  ! write(*,*) a,e,CCSDloop(speciesId)%Fki(a,e)
                end do
              end do
            end if
          end do
        end do
      end do

      ! CCSDloop(speciesId)%Fkc_aa
      do m=1, nop
        do e=nop+1, noc

          CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
            + Allspecies(speciesId)%HF_fs%values(m,e)
          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do n=1, nop
              do f=nop+1, noc
                CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
                  + Allspecies(speciesId)%Tssame(f-nop,n)*spints(speciesId)%valuesp(m,n,e,f)
                ! write(*,*) a,e,CCSDloop(speciesId)%Fkc_aa(a,e)
              end do
            end do
          end if
        end do
      end do

  end subroutine F_onespecies_intermediates

  ! @brief Calculate F intermediates for intra-species
  ! @author CAOM
  subroutine F_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      ! integer :: cont
      integer :: noc, nop, nocs, nops
      integer :: n_sp, num_species
      integer :: a, aa, b, e, ee, i, j, f
      integer :: ff, m, mm, n, nn, ii

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      ! n_sp = Tix2(speciesId, OtherspeciesId, num_species)
      !Tix2 function return the number of inter-specie matrix depending on two species assigned
      n_sp = CCSD_instance%cont

      print*,"noc: ", noc, "nop: ", nop, "nocs: ", nocs, "nops: ", nops

      ! CCSDloop(speciesId)%Fac 1
      do a=nop+1, noc
        do e=nop+1, noc
          do mm=1, nops
            do ee=nops+1, nocs

              CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                + ( 0.25*Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)* &
                    spintm(n_sp)%valuesp(a,mm,e,ee)) ! check this

              do m=1, nop

                CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                  - ( 0.25*Allinterspecies(speciesId)%intau(a-nop,ee-nops,m,mm)* &
                      spintm(n_sp)%valuesp(m,mm,e,ee))
                  ! write(*,*) a,e,CCSDloop(speciesId)%Fac(a,e)
              end do
            end do
          end do

        end do
      end do
      print*, "nops 6"

      ! CCSDloop(speciesId)%Fki(m,i) 2
      do m=1, nop
        do i=1, nop
             
          do mm=1, nops
            do ee=nops+1, nocs

              CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                 + ( 0.25*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)* &
                    spintm(n_sp)%valuesp(m,mm,i,ee))
            
              do e=nop+1, noc

                CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                  + ( 0.25*Allinterspecies(n_sp)%intau(e-nop,ee-nops,i,mm)* &
                      spintm(n_sp)%valuesp(m,mm,e,ee))
                ! write(*,*) a,e,CCSDloop(speciesId)%Fki(m,i)
              end do
            end do
          end do

        end do
      end do
      print*, "nops 7"

      ! CCSDloop(speciesId)%Fkc_aa(m,e-nop) 3
      do m=1, nop
        do e=nop+1, noc

          do mm=1, nops
            do ee=nops+1, nocs
              CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
                + ( 0.25*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)* &
                    spintm(n_sp)%valuesp(m,mm,e,ee))
              ! write(*,*) a,e,CCSDloop(speciesId)%Fkc_aa(m,e-nop)
            end do
          end do

        end do
      end do
      print*, "nops 8"

      ! ! CCSDinter(speciesId)%Fkc_aba(m,e-nop) 4
      do m=1, nop
        do e=nop+1, noc
          do aa=nops+1, nocs
            do ii=1, nops

              CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                - (0.5*spintm(n_sp)%valuesp(m,aa,e,ii))

              do ee=1, nops

                CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                  -(0.25* Allspecies(OtherspeciesId)%Tssame(ee-nops,ii)*spintm(n_sp)%valuesp(m,aa,e,ee))

                if (nops>=2) then ! kind of interaction just for two or more particles of the another species

                  do mm=1, nops

                    CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                      -(0.25* Allspecies(OtherspeciesId)%Tdsame(aa-nops,ee-nops,ii,mm)* &
                          spintm(n_sp)%valuesp(m,mm,e,ee))

                    CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                      +(0.25* Allspecies(OtherspeciesId)%Tssame(aa-nops,mm)* &
                          Allspecies(OtherspeciesId)%Tssame(ee-nops,ii) &
                          *spintm(n_sp)%valuesp(m,mm,e,ee))
                  end do
                end if
              end do

              do mm=1, nops

                CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                  +(0.25* Allspecies(OtherspeciesId)%Tssame(aa-nops,mm)*spintm(n_sp)%valuesp(m,mm,e,ii))
              end do

            end do
          end do

          do a=nop+1, noc
            do i=1, nop
              do mm=1, nops
                do ee=nops+1, nocs

                  CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                    -(0.125* Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* & 
                        spintm(n_sp)%valuesp(m,mm,e,ee))
                end do
              end do
            end do
          end do
          ! write(*,*) a,e,CCSDinter(num_inter)%Fkc_aba(m,e-nop)
        end do
      end do
      print*, "nops 9"

      ! CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) 5
      do mm=1, nops
        do ee=nops+1, nocs

          if(nops>=2) then ! kind of interaction just for two or more particles of the another species

            CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
              + ( 0.5*Allspecies(OtherspeciesId)%HF_fs%values(mm,ee))

            do nn=1, nops
              do ff=nops+1, nocs

                CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
                  + ( 0.5*Allspecies(OtherspeciesId)%Tssame(ff-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,ee,ff))
              end do
            end do
          end if

          do m=1, nop
            do e=nop+1, noc

              CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
                + ( 0.125*Allspecies(speciesId)%Tssame(e-nop,m)* &
                    spintm(n_sp)%valuesp(m,mm,e,ee))
              ! write(*,*) a,e,CCSDloop(speciesId)%Fkca_ab(mm,ee-nops)
            end do
          end do

        end do
      end do
      print*, "nops 10"

  end subroutine F_twospecies_intermediates

  ! @brief Calculate F intermediates for intra-species
  ! @author CAOM
  subroutine F_T2_AB(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer :: noc, nop, nocs, nops
      integer :: num_intersp, num_species
      integer :: a, aa, b, bb, e, ee, i, j
      integer :: f, ff, m, mm, n, nn, ii, jj

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      !Tix2 function return the number of inter-specie matrix depending on two species assigned
      ! num_intersp = Tix2(speciesId, OtherspeciesId, num_species)

      print*,"noc: ", noc, "nop: ", nop, "nocs: ", nocs, "nops: ", nops

      ! CCSDinter(speciesId)%Faca
      do a=nop+1, noc
        do e=nop+1, noc

          CCSDinter(speciesId)%Faca(a-nop,e-nop) = CCSDinter(speciesId)%Faca(a-nop,e-nop) &
            +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e)

          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do m=1, nop
              do f=nop+1, noc
            
                CCSDinter(speciesId)%Faca(a-nop,e-nop) = CCSDinter(speciesId)%Faca(a-nop,e-nop) &
                  - ( Allspecies(speciesId)%Tssame(a-nop,m)* &
                      spints(speciesId)%valuesp(m,a,e,f))
              end do
            end do
          end if
        ! write(*,*) a,e,CCSDinter(speciesId)%Faca(a-nop,e-nop)
        end do
      end do
      print*, "nops 15"

      ! CCSDinter(speciesId)%Fkia
      do m=1, nop
        do i=1, nop

          CCSDinter(speciesId)%Fkia(m,i) = CCSDinter(speciesId)%Fkia(m,i) &
            +(1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i)

          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do n=1, nop
              do e=nop+1, noc
 
                CCSDinter(speciesId)%Fkia(m,i) = CCSDinter(speciesId)%Fkia(m,i) &
                  - ( Allspecies(speciesId)%Tssame(e-nop,n)* &
                      spints(speciesId)%valuesp(m,n,i,e))
              end do
            end do
          end if
        ! write(*,*) a,e,CCSDinter(speciesId)%Fkia(m,i)
        end do
      end do
      print*, "nops 16"

      ! CCSDinter(OtherspeciesId)%Fbcb
      do bb=nops+1, nocs
        do ee=nops+1, nocs

          CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) = CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) &
            +(1 - logic2dbl(bb==ee))*Allspecies(OtherspeciesId)%HF_fs%values(bb,ee)

          if (nops>=2) then ! kind of interaction just for two or more particles of the another species
            do mm=1, nops
              do ff=nops+1, nocs
            
                CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) = CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) &
                  - ( Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      spints(OtherspeciesId)%valuesp(mm,bb,ee,ff))
              end do
            end do
          end if
        ! write(*,*) bb,ee,CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops)
        end do
      end do
      print*, "nops 17"

      ! CCSDinter(OtherspeciesId)%Fkjb
      do mm=1, nops
        do jj=1, nops

          CCSDinter(OtherspeciesId)%Fkjb(mm,jj) = CCSDinter(OtherspeciesId)%Fkjb(mm,jj) &
            +(1 - logic2dbl(mm==jj))*Allspecies(OtherspeciesId)%HF_fs%values(mm,jj)

          if (nops>=2) then ! kind of interaction just for two or more particles of the another species          
            do nn=1, nops
              do ee=nops+1, nocs

                CCSDinter(OtherspeciesId)%Fkjb(mm,jj) = CCSDinter(OtherspeciesId)%Fkjb(mm,jj) &
                  - ( Allspecies(OtherspeciesId)%Tssame(ee-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,jj,ee))
              end do
            end do
          end if
        ! write(*,*) a,e,CCSDinter(OtherspeciesId)%Fkjb(mm,jj)
        end do
      end do
      print*, "nops 18"

  end subroutine F_T2_AB

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

      !Baisc parallelization
      !$OMP PARALLEL
      !$OMP DO
      ! CCSDloop(speciesId)%Wklij
      do m=1, nop
        do n=1, nop
          do i=1, nop
            do j=1, nop

              CCSDloop(speciesId)%Wklij(m,n,i,j) = CCSDloop(speciesId)%Wklij(m,n,i,j) &
                + spints(speciesId)%valuesp(m,n,i,j)
              do e=nop+1, noc
                CCSDloop(speciesId)%Wklij(m,n,i,j) = CCSDloop(speciesId)%Wklij(m,n,i,j) &
                  + (Allspecies(speciesId)%Tssame(e-nop,j)*spints(speciesId)%valuesp(m,n,i,e) &
                      -Allspecies(speciesId)%Tssame(e-nop,i)*spints(speciesId)%valuesp(m,n,j,e))
                do f=nop+1, noc
                  CCSDloop(speciesId)%Wklij(m,n,i,j) = CCSDloop(speciesId)%Wklij(m,n,i,j) &
                    + 0.25*Allspecies(speciesId)%tau(e-nop,f-nop,i,j)*spints(speciesId)%valuesp(m,n,e,f)
                  ! write(*,*) m,n,i,j,CCSDloop(speciesId)%Wklij(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do
      !$OMP END DO

      !$OMP DO
      !CCSDloop(speciesId)%Wabcd
      do a=nop+1, noc
        do b=nop+1, noc
          do e=nop+1, noc
            do f=nop+1, noc

              CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                + spints(speciesId)%valuesp(a,b,e,f)
              do m=1, nop
                CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                  + (-Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(a,m,e,f) &
                      +Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(b,m,e,f))
                do n=1, nop
                  CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                    + 0.25*Allspecies(speciesId)%tau(a-nop,b-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f)
                    ! write(*,*) m,n,i,j,CCSDloop(speciesId)%Wabcd(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do
      !$OMP END DO

      !$OMP DO
      !CCSDloop(speciesId)%Wkbcj
      do m=1, nop
        do b=nop+1, noc
          do e=nop+1, noc
            do j=1, nop
              CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) = CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) &
                + spints(speciesId)%valuesp(m,b,e,j)
              do f=nop+1, noc
                CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) = CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) &
                 + Allspecies(speciesId)%Tssame(f-nop,j)*spints(speciesId)%valuesp(m,b,e,f)
              end do
              do n=1, nop
                CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) = CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) &
                  - Allspecies(speciesId)%Tssame(b-nop,n)*spints(speciesId)%valuesp(m,n,e,j)
                do f=nop+1, noc
                  CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) = CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) &
                    - ((0.5*Allspecies(speciesId)%Tdsame(f-nop,b-nop,j,n) &
                        + Allspecies(speciesId)%Tssame(f-nop,j)*Allspecies(speciesId)%Tssame(b-nop,n))* &
                          spints(speciesId)%valuesp(m,n,e,f))
                  ! write(*,*) m,n,i,j,CCSDloop(speciesId)%Wkbcj(m,n,i,j)
                end do
              end do
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

  end subroutine W_onespecies_intermediates

  !>
  ! @brief Calculate W intermediates for intra-species
  ! @author CAOM
  subroutine W_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer :: noc, nop, nocs, nops
      integer :: n_sp, num_species
      integer :: a, aa, b, e, ee, i, j, f
      integer :: ff, m, mm, n, nn, ii, bb, jj

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      ! num_intersp = Tix2(speciesId, OtherspeciesId, num_species)
      !Tix2 function return the number of inter-specie matrix depending on two species assigned
      n_sp = CCSD_instance%cont

      print*,"noc: ", noc, "nop: ", nop, "nocs: ", nocs, "nops: ", nops

      ! CCSDinter(speciesiD)%Wkkcc 6
      do b=nop+1, noc
        do j=1, nop
          do mm=1, nops
            do ee=nops+1, nocs

              CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                + (0.5*spintm(n_sp)%valuesp(b,mm,j,ee))
              
              do e=nop+1, noc
                do m=1, nop

                  CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                    - (0.75*(Allspecies(speciesId)%Tssame(b-nop,m)*spintm(n_sp)%valuesp(m,mm,j,ee) &
                      - Allspecies(speciesId)%Tssame(e-nop,j)*spintm(n_sp)%valuesp(b,mm,e,ee)))
                  if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
                    CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                      - ((0.5*Allspecies(speciesId)%Tdsame(e-nop,b-nop,j,m)) &
                        + (Allspecies(speciesId)%Tssame(e-nop,j)*Allspecies(speciesId)%Tssame(b-nop,m)))* &
                          spintm(n_sp)%valuesp(m,mm,e,ee)
                  end if
                end do
              end do

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species

                do ff=nops+1, nocs
                  do nn=1, nops
                  
                    CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                      + (0.25*Allinterspecies(OtherspeciesId)%Tdsame(b,ff-nops,j,nn)* &
                          spints(OtherspeciesId)%valuesp(mm,nn,ee,ff))
                  end do
                end do
              end if
              ! write(*,*) m,n,i,j,CCSDloop(speciesId)%Wklij(m,n,i,j)
            end do
          end do
        end do
      end do
      print*, "nops 11"

  end subroutine W_twospecies_intermediates

  !>
  ! @brief Calculate W intermediates for intra-species
  ! @author CAOM
  subroutine W_T2_AB(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer :: noc, nop, nocs, nops
      integer :: n_sp, num_species
      integer :: a, aa, b, bb, e, ee, i, j
      integer :: f, ff, m, mm, n, nn, ii, jj

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species
      n_sp = CCSD_instance%cont

      !Tix2 function return the number of inter-specie matrix depending on two species assigned
      ! num_intersp = Tix2(speciesId, OtherspeciesId, num_species)

      print*,"noc: ", noc, "nop: ", nop, "nocs: ", nocs, "nops: ", nops

      ! CCSDinter(speciesId)%Waka
      do m=1, nop
        do bb=1+nops, nocs
          do i=1, nop
            do jj=1, nops

              CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                + spintm(n_sp)%valuesp(m,bb,i,jj)

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops
                    CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                      + (0.5*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        spintm(n_sp)%valuesp(m,mm,i,ee))
                  end do
                end do
              end if

              do e=nop+1, noc
                CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                  + (0.5*Allspecies(speciesId)%Tssame(e-nop,i)* &
                      spintm(n_sp)%valuesp(m,bb,e,jj))
              end do
              ! write(*,*) m,bb-nops,i,jj,CCSDinter(speciesId)%Waka(m,bb-nops,i,jj)
            end do
          end do
        end do
      end do
      print*, "nops 19"

      ! CCSDinter(speciesId)%Wcia
      do a=nop+1, noc
        do bb=1+nops, nocs
          do e=nop+1, noc
            do jj=1, nops

              CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                + spintm(n_sp)%valuesp(a,bb,e,jj)

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops
                    CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                      + (0.5*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        spintm(n_sp)%valuesp(a,mm,e,ee))
                  end do
                end do
              end if

              do m=1, nop
                CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                  + (0.5*Allspecies(speciesId)%Tssame(a-nop,m)* &
                       spintm(n_sp)%valuesp(m,bb,e,jj))
              end do
              ! write(*,*) a-nop,bb-nops,e-nop,jj,CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj)
            end do
          end do
        end do
      end do
      print*, "nops 20"

      ! CCSDinter(speciesId)%Wbkb
      do a=nop+1, noc
        do mm=1, nops
          do i=1, nop
            do jj=1, nops

              CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                + spintm(n_sp)%valuesp(a,mm,i,jj)

              do e=nop+1, noc
                do m=1, nop
                  CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                    + (0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                      spintm(n_sp)%valuesp(m,mm,e,jj))
                end do
              end do

              do ee=nops+1, nocs
                CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                      spintm(n_sp)%valuesp(a,mm,i,ee))
              end do
              ! write(*,*) a-nop,mm,i,jj,CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj)
            end do
          end do
        end do
      end do
      print*, "nops 21"

      ! CCSDinter(speciesId)%Wcjb
      do a=nop+1, noc
        do bb=1+nops, nocs
          do i=1, nop
            do ee=nops+1, nocs

              CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                + spintm(n_sp)%valuesp(a,bb,i,ee)

              do e=nop+1, noc
                do m=1, nop
                  CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                    + (0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                      spintm(n_sp)%valuesp(m,mm,e,ee))
                end do
              end do

              do mm=1, nops
                CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      spintm(n_sp)%valuesp(a,mm,i,ee))
              end do
              ! write(*,*) a-nop,bb-nops,i,ee-nops,CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops)
            end do
          end do
        end do
      end do
      print*, "nops 22"

      ! CCSDinter(speciesId)%Wakic_a
      do a=nop+1, noc
        do m=1, nop
          do i=1, nop
            do e=nop+1, noc

              CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                  + (Allspecies(speciesId)%HF_fs%values(m,e)* &
                      Allspecies(speciesId)%Tssame(a-nop,i))

              if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
                CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                  + spints(speciesId)%valuesp(a,m,i,e)

                do n=1, nop
                  CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (Allspecies(speciesId)%Tssame(a-nop,n)* &
                      spints(speciesId)%valuesp(m,n,i,e))
                end do

                do f=nop+1, noc
                  CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (Allspecies(speciesId)%Tssame(f-nop,i)* &
                      spints(speciesId)%valuesp(m,a,e,f))
                end do
              end if

              do ee=nops+1, nocs
                do mm=1, nops
                  CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (0.125*Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                        spintm(n_sp)%valuesp(m,mm,e,ee))
                end do
              end do
              ! write(*,*) a-nop,m,i,e-nop,CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop)
            end do
          end do
        end do
      end do
      print*, "nops 23"

      ! CCSDinter(speciesId)%Wbkjc_b 
      do bb=nops+1, nocs
        do mm=1, nops
          do jj=1, nops
            do ee=nops+1, nocs

              CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                  + (Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)* &
                      Allspecies(OtherspeciesId)%Tssame(bb-nops,jj))

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                  + spints(OtherspeciesId)%valuesp(bb,mm,jj,ee)

                do nn=1, nops
                  CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (Allspecies(OtherspeciesId)%Tssame(bb-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,jj,ee))
                end do

                do ff=nops+1, nocs
                  CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (Allspecies(OtherspeciesId)%Tssame(ff-nops,jj)* &
                      spints(OtherspeciesId)%valuesp(mm,bb,ee,ff))
                end do
              end if

              do e=nop+1, noc
                do m=1, nop
                  CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (0.125*Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)* &
                        spintm(n_sp)%valuesp(m,mm,e,ee))
                end do
              end do
              ! write(*,*) bb-nops,mm,jj,ee-nops,CCSDinter(speciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops)
            end do
          end do
        end do
      end do
      print*, "nops 24"

      ! CCSDinter(speciesId)%Wklcd_a
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        do m=1, nop
          do n=1, nop
            do e=nop+1, noc
              do f=nop+1, noc

                do bb=nops+1, nocs
                  do jj=1, nops
                    CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop) = CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop) &
                      + spints(speciesId)%valuesp(m,n,e,f)* &
                        (Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj) &
                          + (Allspecies(speciesId)%Tssame(e-nop,m)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)))
                  end do
                end do
                ! write(*,*) m,n,e-nop,f-nop,CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop)
              end do
            end do
          end do
        end do
      end if
      print*, "nops 25"

      ! CCSDinter(OtherspeciesId)%Wklcd_b
      if (nops>=2) then ! kind of interaction just for two or more particles of the another species
        do mm=1, nops
          do nn=1, nops
            do ee=nops+1, nocs
              do ff=nops+1, nocs

                do a=nop+1, noc
                  do i=1, nop
                    CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops) = CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops) &
                      + spints(OtherspeciesId)%valuesp(mm,nn,ee,ff)* &
                        (Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm) &
                          + (Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)))
                  end do
                end do
                ! write(*,*) mm,nn,ee-nops,ff-nops,CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops)
              end do
            end do
          end do
        end do
      end if
      print*, "nops 26"

      ! CCSDinter(speciesId)%Wakic
      do a=nop+1, noc
        do mm=1, nops
          do i=1, nop
            do ee=nops+1, nocs

              CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                + spintm(n_sp)%valuesp(a,mm,i,ee) &
                  + (Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)* &
                    Allspecies(speciesId)%Tssame(a-nop,i))

              do m=1, nop
                CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                  - (0.5*Allspecies(speciesId)%Tssame(a-nop,m)* &
                    spintm(n_sp)%valuesp(m,mm,i,ee))
              end do

              do e=nop+1, noc
                CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                  + (0.5*Allspecies(speciesId)%Tssame(e-nop,i)* &
                    spintm(n_sp)%valuesp(a,mm,e,ee))
              end do

              if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
                do e=nop+1, noc
                  do m=1, nop
                    CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                      + ( Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(e-nop,m) &
                        + (0.5*Allspecies(speciesId)%tau(a-nop,e-nop,i,m)) &
                          - (0.25*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,i)))* &
                           spintm(n_sp)%valuesp(m,mm,e,ee)
                  end do
                end do
              end if
              ! write(*,*) a-nop,mm,i,ee-nops,CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops)
            end do
          end do
        end do
      end do
      print*, "nops 27"

      ! CCSDinter(speciesId)%Wbkjc
      do m=1, nop
        do bb=nops+1, nocs
          do e=nop+1, noc
            do jj=1, nops

              CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                + spintm(n_sp)%valuesp(m,bb,e,jj) &
                  + (Allspecies(speciesId)%HF_fs%values(m,e)* &
                    Allspecies(OtherspeciesId)%Tssame(bb-nops,jj))

              do mm=1, nops
                CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                  - (0.5*Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                    spintm(n_sp)%valuesp(m,mm,e,jj))
              end do

              do ee=nops+1, nocs
                CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                    spintm(n_sp)%valuesp(m,bb,e,ee))
              end do

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops
                    CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                      +((Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)) &
                        + (0.5*Allspecies(OtherspeciesId)%tau(bb-nops,ee-nops,jj,mm)) &
                          - (0.25*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,mm,jj)) )* &
                          spintm(n_sp)%valuesp(m,mm,e,ee)
                  end do
                end do
              end if
              ! write(*,*) m,bb-nops,e-nop,jj,CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj)
            end do
          end do
        end do
      end do
      print*, "nops 28"

  end subroutine W_T2_AB

  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSD_loop(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, optional, intent(in) :: OtherspeciesId

      integer :: i_counterID(10)
      integer :: n_intersp(10)
      integer :: times_i
      integer :: noc, nocs, nop, nops
      integer :: num_species
      integer :: n_sp=0
      integer :: OtspId=0
      
      integer :: max, min, e_cont
      integer :: a, b, e, i, j,jj
      integer :: aa, ii, f, m, n
      integer :: num_inter!=1
      real(8) :: prev_ccsdE
      real(8) :: tmp_ccsdE
      real(8) :: prev_ccsdE_int
      real(8) :: tmp_ccsdE_int
      real(8) :: ccsdE=0.0_8
      real(8) :: convergence = 1.0_8
      real(8) :: convergence_intra = 1.0_8
      real(8) :: convergence_inter = 1.0_8
      real(8) :: ccsdE_int=0.0_8

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (ccsdE /= 0.0_8) ccsdE = 0.0_8
      if (convergence_intra /= 1.0D-8) convergence_intra = 1.0_8
      if (convergence_inter /= 1.0D-8) convergence_inter = 1.0_8
      if (ccsdE_int /= 0.0_8) ccsdE_int = 0.0_8
      if (present(OtherspeciesId)) OtspId = OtherspeciesId

      !Initialization of private variables from public variables
      noc = Allspecies(speciesId)%noc
      ! nocs = Allspecies(OtspId)%noc
      nop = Allspecies(speciesId)%nop
      ! nops = Allspecies(OtspId)%nop
      num_species = CoupledCluster_instance%num_species
      ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
      print*, "T1T2_constructor", convergence, noc, nop, speciesId

      times_i = CoupledCluster_instance%times_intersp
      max = CCSD_instance%max
      min = CCSD_instance%min

      ! n_sp = Tix2(speciesId, OtspId, num_species)
      ! num_inter = num_inter + 1
      
      do while (convergence >= 1.0D-8)
        !**change position of do while
        !do i=1, num_species
          !all intra CCSD
        !end do
        !do i=min, max
        ! all inter-species
        !end do
        prev_ccsdE = ccsdE

        prev_ccsdE_int = ccsdE_int

        print*, "speciesId CCSD_loop: ", speciesId

        call CCSD_T1T2_constructor(speciesId)
        call CCSD_loop_constructor(speciesId)

        !if there are inter-species
        if ((max/=0) .and. (min/=0)) then
          num_inter = 0
          CCSD_instance%cont = CCSD_instance%aux_cont
          do jj=min+1, max
          
            ! call CCSD_constructor_inter(min, jj)
            print*, "jj: ", jj
            call CCSD_T2AB_constructor(min, jj, num_inter)
            call CCSD_loop_constructor_inter(min, jj, num_inter)

            num_inter = num_inter + 1
            CCSD_instance%cont = CCSD_instance%cont + 1                
          end do
        end if

        !intermediates loop for:
        !singles excitations
        print*, "F_onespecies_intermediates(): "
        call F_onespecies_intermediates(speciesId)
        !doubles excitations
        print*, "W_onespecies_intermediates(): "
        if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
          call W_onespecies_intermediates(speciesId)
        end if

        !If there are interspecies?
        if ((max/=0) .and. (min/=0)) then
          num_inter = 0
          CCSD_instance%cont = CCSD_instance%aux_cont
          do jj=min+1, max
                print*, "ciclo: jj: ", jj
                call F_twospecies_intermediates(min, jj, num_inter)
                call W_twospecies_intermediates(min, jj, num_inter)
                call F_T2_AB(min, jj, num_inter)
                call W_T2_AB(min, jj, num_inter)
        
            num_inter = num_inter + 1
            CCSD_instance%cont = CCSD_instance%cont + 1                
          end do
        end if

        ! Resolve CCSD equation of energy

        tmp_ccsdE=0.0_8

        ! for same species
        do i=1, nop
          do a=nop+1, noc
            tmp_ccsdE = tmp_ccsdE + Allspecies(speciesId)%HF_fs%values(i,a)*Allspecies(speciesId)%Tssame(a-nop,i)
            if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
              do j=1, nop
                do b=nop+1, noc
                  tmp_ccsdE = tmp_ccsdE + (0.25*spints(speciesId)%valuesp(i,j,a,b)* &
                      Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                        + 0.5*spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tssame(a-nop,i)* &
                          Allspecies(speciesId)%Tssame(b-nop,j))
                end do
              end do
            end if
          end do
        end do

        ccsdE = tmp_ccsdE

        if (times_i>0) then
          tmp_ccsdE_int=0.0_8
          ! for different species OtspId
          e_cont = CCSD_instance%aux_cont
          print*, "energy CCSD-APMO: ", min+1, max
          do jj=min+1, max

            nops = Allspecies(jj)%nop
            nocs = Allspecies(jj)%noc
            do i=1, nop
              do a=nop+1, noc
                do ii=1, nops
                  do aa=nops+1, nocs
                    tmp_ccsdE_int = tmp_ccsdE_int + (0.25*spintm(e_cont)%valuesp(i,ii,a,aa)* &
                        Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) ) &
                          + (0.5*spintm(e_cont)%valuesp(i,ii,a,aa)*Allspecies(speciesId)%Tssame(a-nop,i)* &
                            Allspecies(jj)%Tssame(aa-nops,ii))
                  end do
                end do
              end do
            end do
            e_cont = e_cont + 1
            print*, "energy CCSD-APMO"
          end do
        end if

        ccsdE_int = tmp_ccsdE_int

        !change in values for the intra-species loop
        convergence_intra = abs( ccsdE - prev_ccsdE )
        !change in values for the inter-species loop
        convergence_inter = abs( ccsdE_int - prev_ccsdE_int )
        ! total convergence
        if (times_i>0) then
          convergence = abs( (ccsdE+ccsdE_int) - (prev_ccsdE+prev_ccsdE_int) )
        else
          convergence = convergence_intra
        end if

        write (*,*) ccsdE, "CCSD Energy ", prev_ccsdE, "previous Energy" 
        write (*,*) convergence, "Convergence " 
        ! if ((convergence > 10) .and. speciesId>1) stop "test"
        ! Resolve T1 and T2 amplitude equations
        if (times_i>0) then
          do jj=min, max
            call CCSD_T1(jj)
          end do
        else
          call CCSD_T1(speciesId)
        end if
        call CCSD_T2(speciesId)

        ! T1 and T2 equation energies for interspecies
        if (times_i>0) then
          num_inter = 0
          CCSD_instance%cont = CCSD_instance%aux_cont
          do jj=min+1, max
              num_inter = num_inter + 1
              print*, "before CCSD_T1_inter()", min, jj
              call CCSD_T1_inter(min, jj, num_inter)
              print*, "before CCSD_T2_inter()", min, jj
              call CCSD_T2_inter(min, jj, num_inter)
              print*, "CCSD_T2_AB()"
              call CCSD_T2_AB(min, jj, num_inter)
              CCSD_instance%cont = CCSD_instance%cont + 1
          end do
        end if

        if (convergence > 100) then 
          stop "Error: There are not convergence. The differences between energies is more than 100 eV"
        end if

      end do

      CoupledCluster_instance%CCSD_E_intra(speciesId) = ccsdE
      
  end subroutine CCSD_loop

  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSD_same_species(speciesId, e_ccsd)
      implicit none

      integer, intent(in) :: speciesId
      real(8), intent(in) :: e_ccsd

      integer :: i_counterID(10)
      integer :: n_intersp(10)
      integer :: times_i
      integer :: noc, nocs, nop, nops
      integer :: num_species
      integer :: n_sp=0
      
      integer :: a, b, e, i, j,jj
      integer :: aa, ii, f, m, n
      integer :: num_inter!=1
      real(8) :: prev_ccsdE
      real(8) :: tmp_ccsdE
      real(8) :: prev_ccsdE_int
      real(8) :: tmp_ccsdE_int
      real(8) :: ccsdE=0.0_8
      real(8) :: convergence = 1.0_8

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (ccsdE /= 0.0_8) ccsdE = 0.0_8

      !Initialization of private variables from public variables
      noc = Allspecies(speciesId)%noc
      nop = Allspecies(speciesId)%nop
      num_species = CoupledCluster_instance%num_species
      ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
      print*, "T1T2_constructor", convergence, noc, nop, speciesId

      ! times_i = CoupledCluster_instance%times_intersp
      ! max = CCSD_instance%max
      ! min = CCSD_instance%min

        !**change position of do while
        !do i=1, num_species
          !all intra CCSD
        !end do
        !do i=min, max
        ! all inter-species
        !end do
      prev_ccsdE = e_ccsd

      print*, "speciesId CCSD_loop: ", speciesId

      call CCSD_T1T2_constructor(speciesId)
      call CCSD_loop_constructor(speciesId)

      !intermediates loop for:
      !singles excitations
      print*, "F_onespecies_intermediates(): "
      call F_onespecies_intermediates(speciesId)
      !doubles excitations
      print*, "W_onespecies_intermediates(): "
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        call W_onespecies_intermediates(speciesId)
      end if

      ! Resolve CCSD equation of energy

      tmp_ccsdE=0.0_8

      ! for same species
      do i=1, nop
        do a=nop+1, noc
          tmp_ccsdE = tmp_ccsdE + Allspecies(speciesId)%HF_fs%values(i,a)*Allspecies(speciesId)%Tssame(a-nop,i)
          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do j=1, nop
              do b=nop+1, noc
                tmp_ccsdE = tmp_ccsdE + (0.25*spints(speciesId)%valuesp(i,j,a,b)* &
                    Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                      + 0.5*spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tssame(a-nop,i)* &
                        Allspecies(speciesId)%Tssame(b-nop,j))
              end do
            end do
          end if
        end do
      end do

      ccsdE = tmp_ccsdE

      !change in values for the intra-species loop
      convergence = abs( ccsdE - prev_ccsdE )

      write (*,*) ccsdE, "CCSD Energy same species", prev_ccsdE, "previous Energy" 
      write (*,*) convergence, "Convergence " 

      call CCSD_T1(speciesId)
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        call CCSD_T2(speciesId)
      end if

      if (convergence > 100) then 
        stop "Error: There are not convergence. The differences between energies is more than 100 eV"
      end if

      CoupledCluster_instance%CCSD_E_intra(speciesId) = ccsdE
      CCSD_instance%convergence_same(speciesId) = convergence
      
  end subroutine CCSD_same_species
  
  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSD_diff_species(speciesId, OtherspeciesId, e_ccsd)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      real(8), intent(in) :: e_ccsd

      integer :: i_counterID(10)
      integer :: n_intersp(10)
      integer :: times_i
      integer :: noc, nocs, nop, nops
      integer :: num_species
      integer :: n_sp=0
      
      integer :: max, min, e_cont
      integer :: a, b, e, i, j,jj
      integer :: aa, ii, f, m, n
      integer :: num_inter
      real(8) :: prev_ccsdE_int
      real(8) :: tmp_ccsdE_int
      real(8) :: convergence = 1.0_8
      real(8) :: ccsdE_int=0.0_8

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (ccsdE_int /= 0.0_8) ccsdE_int = 0.0_8

      !Initialization of private variables from public variables
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      num_species = CoupledCluster_instance%num_species
      ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
      print*, "T1T2_constructor", convergence, noc, nop, speciesId

      prev_ccsdE_int = e_ccsd

      print*, "speciesId CCSD_loop: ", speciesId

      num_inter = CCSD_instance%num_i
      e_cont = CCSD_instance%e_cont

      ! call CCSD_constructor_inter(min, jj)
      print*, "jj: ", OtherspeciesId
      call CCSD_T2AB_constructor(speciesId, OtherspeciesId, num_inter)
      call CCSD_loop_constructor_inter(speciesId, OtherspeciesId, num_inter)

      !intermediates loop for:
      !If there are interspecies?
      print*, "ciclo: OtherspeciesId: ", OtherspeciesId
      call F_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      call W_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      call F_T2_AB(speciesId, OtherspeciesId, num_inter)
      call W_T2_AB(speciesId, OtherspeciesId, num_inter)
        
      ! Resolve CCSD equation of energy

      tmp_ccsdE_int=0.0_8
      print*, "energy CCSD-APMO: ", speciesId+1, max

      do i=1, nop
        do a=nop+1, noc
          do ii=1, nops
            do aa=nops+1, nocs
              tmp_ccsdE_int = tmp_ccsdE_int + (0.25*spintm(e_cont)%valuesp(i,ii,a,aa)* &
                  Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) ) &
                    + (0.5*spintm(e_cont)%valuesp(i,ii,a,aa)*Allspecies(speciesId)%Tssame(a-nop,i)* &
                      Allspecies(OtherspeciesId)%Tssame(aa-nops,ii))
            end do
          end do
        end do
      end do

      print*, "energy CCSD-APMO"

      ccsdE_int = tmp_ccsdE_int

      !change in values for the inter-species loop
      convergence = abs( ccsdE_int - prev_ccsdE_int )

      write (*,*) ccsdE_int, "CCSD Energy inter-species ", prev_ccsdE_int, "previous Energy inter-species" 
      write (*,*) convergence, "Convergence inter-species" 

      ! T1 and T2 equation energies for interspecies
      print*, "before CCSD_T1_inter()", speciesId, OtherspeciesId
      call CCSD_T1_inter(speciesId, OtherspeciesId, num_inter)
      print*, "before CCSD_T2_inter()", speciesId, OtherspeciesId
      call CCSD_T2_inter(speciesId, OtherspeciesId, num_inter)
      print*, "CCSD_T2_AB()"
      call CCSD_T2_AB(speciesId, OtherspeciesId, num_inter)
  
      CCSD_instance%e_cont = e_cont
      CCSD_instance%num_i = num_inter

      if (convergence > 100) then 
        stop "Error: There are not convergence. The differences between energies is more than 100 eV"
      end if

      CoupledCluster_instance%CCSD_E_inter(num_inter) = ccsdE_int
      CCSD_instance%convergence_diff(num_inter) = convergence
      
      num_inter = num_inter + 1!final
      CCSD_instance%cont = CCSD_instance%cont + 1!final
      e_cont = e_cont + 1!final

  end subroutine CCSD_diff_species

  !>
  ! @brief Calculate T1 energy equations for intra-species
  ! @author CAOM
  subroutine CCSD_T1(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nocs, nop, nops
      integer :: num_species
      integer :: times_i
      integer :: a, b, e, f, i, j, m, n

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      num_species = CoupledCluster_instance%num_species
      times_i = CoupledCluster_instance%times_intersp
      
      !Basic parallelization
      !!$OMP PARALLEL
      !!$OMP DO
      ! T^{a}_{i}D^{a}_{i} = ...
      do a=nop+1, noc
        do i=1, nop
          CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
            + Allspecies(speciesId)%HF_fs%values(i,a)
          do e=nop+1, noc
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop(speciesId)%Fac(a-nop,e-nop)
          end do
          do m=1, nop
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + (-Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop(speciesId)%Fki(m,i))
            if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
              do e=nop+1, noc
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop(speciesId)%Fkc_aa(m,e-nop)
                do f=nop+1, noc
                  CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                    + (-0.5*Allspecies(speciesId)%Tdsame(e-nop,f-nop,i,m)*spints(speciesId)%valuesp(m,a,e,f))
                end do
                do n=1, nop
                  CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                    + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,n)*spints(speciesId)%valuesp(n,m,e,i))
                end do
              end do
            end if
          end do
          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do n=1,nop
              do f=nop+1, noc
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + (-Allspecies(speciesId)%Tssame(f-nop,n)*spints(speciesId)%valuesp(n,a,i,f))
              end do
            end do
          end if
          if (times_i==0) then
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit(speciesId)%Dai(a,i)
            Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
          end if
          ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
        end do
      end do
      print*, "CCSD_T1"
      !!$OMP END DO
      !!$OMP END PARALLEL

  end subroutine CCSD_T1

  !>
  ! @brief Calculate T2 energy equations for intra-species
  ! @author CAOM
  subroutine CCSD_T2(speciesId)
      implicit none

      integer, intent(in) :: speciesId

      integer noc, nocs, nop, nops
      integer :: num_species
      integer :: times_i
      integer :: a, b, e, f, i, j, m, n

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      num_species = CoupledCluster_instance%num_species
      times_i = CoupledCluster_instance%times_intersp
      print*, "CCSD_T2"
      
      !Basic parallelization
      !!$OMP PARALLEL
      !!$OMP DO 
      ! T^{ab}_{ij}D^{ab}_{ij} = ...
      ! if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
      do a=nop+1, noc
         do b=nop+1, noc
            do i=1, nop
               do j=1, nop
                  CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                    + spints(speciesId)%valuesp(i,j,a,b) !A
                  ! 1er ciclo
                  do e=nop+1, noc
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,j)*CCSDloop(speciesId)%Fac(b-nop,e-nop) &
                        -Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,j)*CCSDloop(speciesId)%Fac(a-nop,e-nop)) !B
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (Allspecies(speciesId)%Tssame(e-nop,i)*spints(speciesId)%valuesp(a,b,e,j) &
                        -Allspecies(speciesId)%Tssame(e-nop,j)*spints(speciesId)%valuesp(a,b,e,i)) !G
                     do f=nop+1, noc
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*Allspecies(speciesId)%tau(e-nop,f-nop,i,j)*CCSDloop(speciesId)%Wabcd(a-nop,b-nop,e-nop,f-nop) !D
                     end do
                     do m=1, nop
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,j)* &
                              Allspecies(speciesId)%Tssame(b-nop,m)*CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
                                +0.5*Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,j)* &
                                  Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop(speciesId)%Fkc_aa(m,e-nop)) !B'
                     end do
                  end do
                  ! 2do ciclo
                  do m=1, nop
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (-Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,m)*CCSDloop(speciesId)%Fki(m,j) &
                          +Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,m)*CCSDloop(speciesId)%Fki(m,i)) !C
                     CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      + (-Allspecies(speciesId)%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,i,j) &
                          +Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,i,j)) !H
                     do n=1, nop
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*Allspecies(speciesId)%tau(a-nop,b-nop,m,n)*CCSDloop(speciesId)%Wklij(m,n,i,j) !E
                     end do
                     do e=nop+1, noc
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,m)* &
                              Allspecies(speciesId)%Tssame(e-nop,j)*CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
                                +0.5*Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,m)* &
                                  Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop(speciesId)%Fkc_aa(m,e-nop)) !C'
                        CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                          + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,j) &
                            - Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%Tssame(a-nop,m)* &
                              spints(speciesId)%valuesp(m,b,e,j) &
                                -Allspecies(speciesId)%Tdsame(a-nop,e-nop,j,m)*CCSDloop(speciesId)%Wkbcj(m,b-nop,e-nop,i) &
                                  + Allspecies(speciesId)%Tssame(e-nop,j)*Allspecies(speciesId)%Tssame(a-nop,m)* &
                                    spints(speciesId)%valuesp(m,b,e,i) &
                                    -Allspecies(speciesId)%Tdsame(b-nop,e-nop,i,m)*CCSDloop(speciesId)%Wkbcj(m,a-nop,e-nop,j) &
                                      - Allspecies(speciesId)%Tssame(e-nop,i)*Allspecies(speciesId)%Tssame(b-nop,m)* &
                                        spints(speciesId)%valuesp(m,a,e,j) &
                                        + Allspecies(speciesId)%Tdsame(b-nop,e-nop,j,m)* &
                                          CCSDloop(speciesId)%Wkbcj(m,a-nop,e-nop,i)- Allspecies(speciesId)%Tssame(e-nop,j)* &
                                          Allspecies(speciesId)%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,e,i) !F
                     end do
                  end do
                  ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
                  if (times_i==0) then
                    CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                      /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(speciesId)%HF_fs%values(j,j) &
                        -Allspecies(speciesId)%HF_fs%values(a,a) -Allspecies(speciesId)%HF_fs%values(b,b))
                  
                    Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j)

                    Allspecies(speciesId)%ttau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                      + 0.5*(Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                        - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j))
                    Allspecies(speciesId)%tau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                      + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                        - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j)
                  end if
                  ! write(*,*) a,b,i,j,Allspecies(speciesId)%Tdsame(a,b,i,j),CCSDT1T2(speciesId)%Tabij(a,b,i,j)
               end do
            end do
         end do
      end do
      ! end if
      !!$OMP END DO
      !!$OMP END PARALLEL

  end subroutine CCSD_T2

  !>
  ! @brief Calculate T1 energy equations for inter-species
  ! @author CAOM
  subroutine CCSD_T1_inter(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops
      integer :: a, b, e, ee, f, i
      integer :: j, m, mm, n, aa, ii
      integer :: n_sp

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      n_sp = CCSD_instance%cont

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_inte: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops
      ! T^{a}_{i}D^{a}_{i} = ...
      do a=nop+1, noc
        do i=1, nop

          do ee=nops+1, nocs
            do mm=1, nops

              CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                + (Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                    CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops))
              do m=1, nop
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  - Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,m,mm)* &
                      (0.25*spintm(n_sp)%valuesp(m,mm,i,ee))

                do e=nop+1, noc
                  CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                    - Allspecies(speciesId)%Tssame(e-nop,m)*CCSDinter(speciesId)%Fkc_aba(m,e-nop)
                end do
              end do

              do e=nop+1, noc
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + Allinterspecies(speciesId)%Tdsame(e-nop,ee-nops,i,mm)* &
                      (0.25*spintm(n_sp)%valuesp(a,mm,e,ee))
              end do

            end do
          end do

          CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit(speciesId)%Dai(a,i)
          Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
          ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
        end do
      end do
      print*, "nops 12"
      
  end subroutine CCSD_T1_inter

  !>
  ! @brief Calculate T2 energy equations for inter-species
  ! @author CAOM
  subroutine CCSD_T2_inter(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops
      integer :: a, b, e, ee, f, i
      integer :: j, m, mm, n, aa, ii
      integer :: n_sp

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      n_sp = CCSD_instance%cont

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_inte: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

      ! T^{ab}_{ij}D^{ab}_{ij} = ...
      do a=nop+1, noc
        do b=nop+1, noc
          do i=1, nop
            do j=1, nop

              do ee=nops+1, nocs
                do mm=1, nops
                  CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                   + 0.5*((Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                      CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops)) &
                    - (Allinterspecies(speciesId)%Tdsame(b-nop,ee-nops,i,mm)* &
                        CCSDinter(speciesId)%Wkkcc(a-nop,mm,j,ee-nops)) &
                      - (Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,j,mm)* &
                          CCSDinter(speciesId)%Wkkcc(b-nop,mm,i,ee-nops)) &
                        + (Allinterspecies(speciesId)%Tdsame(b-nop,ee-nops,j,mm)* &
                          CCSDinter(speciesId)%Wkkcc(a-nop,mm,i,ee-nops)) )
                end do
              end do
              ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
              CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j) &
                /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(speciesId)%HF_fs%values(j,j) &
                  -Allspecies(speciesId)%HF_fs%values(a,a) -Allspecies(speciesId)%HF_fs%values(b,b))
                  
              Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) = CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j)

              Allspecies(speciesId)%ttau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + 0.5*(Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                  - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j))
              Allspecies(speciesId)%tau(a-nop,b-nop,i,j) = Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) &
                + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(b-nop,j) &
                  - Allspecies(speciesId)%Tssame(b-nop,i)*Allspecies(speciesId)%Tssame(a-nop,j)
              !case for inter-species... is it necessary?
  
              ! write(*,*) a,b,i,j,Allspecies(speciesId)%Tdsame(a,b,i,j),CCSDT1T2(speciesId)%Tabij(a,b,i,j)
            end do
          end do
        end do
      end do
      print*, "nops 13"
      
  end subroutine CCSD_T2_inter

  !>
  ! @brief Calculate T1 and T2 energy equations for intra-species
  ! @author CAOM
  subroutine CCSD_T2_AB(speciesId, OtherspeciesId, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops
      integer num_species, n_sp
      integer :: a, b, bb, e, ee, f, ff
      integer :: i, j, jj, m, mm, n, nn
      integer :: aa, ii, c, cc, k, kk

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      ! number of transformed integrals matrix for speciesId and OtherspeciesId
      ! n_sp = Tix2(speciesId, OtherspeciesId, num_species)
      n_sp = CCSD_instance%cont

      ! T^{aB}_{iJ}D^{aB}_{iJ} = ...
      do a=nop+1, noc
        do bb=nops+1, nocs
          do i=1, nop
            do jj=1, nops

              do e=nop+1, noc
                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 + 0.5*(Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,i,jj)* &
                    CCSDinter(speciesId)%Faca(a-nop,e-nop))

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  + 0.5*(Allspecies(speciesId)%Tssame(e-nop,i)* &
                      CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj))

                do mm=1, nops
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    - 0.5*(Allinterspecies(speciesId)%tau(e-nop,bb-nops,i,mm)* &
                        spintm(n_sp)%valuesp(a,bb,e,jj))
                end do

                do ee=nops+1, nocs
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    + 0.5*(Allinterspecies(speciesId)%tau(e-nop,ee-nops,i,jj)* &
                        spintm(n_sp)%valuesp(a,bb,e,ee))
                end do
              end do
              !print*, "nops 29"

              do ee=nops+1, nocs
                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 + 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,jj)* &
                    CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops))

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  + 0.5*(Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                      CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops))

                do m=1, nop
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    - 0.5*(Allinterspecies(speciesId)%tau(a-nop,ee-nops,m,jj)* &
                        spintm(n_sp)%valuesp(m,bb,i,ee))
                end do
              end do
              !print*, "nops 30"

              do m=1, nop
                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 - 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,m,jj)* &
                    CCSDinter(speciesId)%Fkia(m,i))

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  - 0.5*(Allspecies(speciesId)%Tssame(a-nop,m)* &
                      CCSDinter(speciesId)%Waka(m,bb-nops,i,jj))

                do e=nop+1, noc
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                   + 0.5*(Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)* &
                      CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop))

                  if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                  
                    CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                      + 0.5*(Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                          CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj))

                    do f=nop+1, noc
                      do n=1, nop
                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%Tdsame(a-nop,f-nop,i,n)* &
                            CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop))

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%tau(a-nop,e-nop,i,m)* &
                            Allinterspecies(speciesId)%Tdsame(f-nop,bb-nops,n,jj)* &
                              spints(speciesId)%valuesp(m,n,e,f))
                      end do
                    end do
                  end if

                  do ee=nops+1, nocs
                    !\check{\tau} alpha
                    if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_a(a-nop,e-nop,ee-nops,m,i,jj)* &
                            spintm(n_sp)%valuesp(m,bb,e,ee))
                    end if

                    do mm=1, nops
                      if ((nop>=2) .and. (nops>=2)) then ! kind of interaction for two or more particles for each species
                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%ttau(a-nop,e-nop,i,m)* &
                            Allspecies(OtherspeciesId)%ttau(bb-nops,ee-nops,jj,mm)* &
                              spintm(n_sp)%valuesp(m,mm,e,ee))
                      end if

                    end do
                  end do

                  if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                    do mm=1, nops
                      !\check{\tau} alpha
                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_a(a-nop,e-nop,bb-nops,i,m,mm)* &
                            spintm(n_sp)%valuesp(m,mm,e,jj))
                    end do
                  end if

                end do
              end do
              !print*, "nops 31"

              do mm=1, nops
                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 - 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,mm)* &
                    CCSDinter(OtherspeciesId)%Fkjb(mm,jj))

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  - 0.5*(Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj))

                do ee=nops+1, nocs
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                   + 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                      CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops))

                  if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                    CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                     + 0.5*(Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops))

                    do ff=nops+1, nocs
                      do nn=1, nocs
                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(OtherspeciesId)%Tdsame(bb-nops,ff-nops,jj,nn)* &
                            CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops))

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(OtherspeciesId)%tau(bb-nops,ff-nops,jj,nn)* &
                            Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                              spints(OtherspeciesId)%valuesp(mm,nn,ee,ff))
                      end do
                    end do
                  end if

                  if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                    do e=nop+1, nocs
                    !\check{\tau} beta
                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(OtherspeciesId)%chtau_b(bb-nops,ee-nops,e-nop,mm,jj,i)* &
                            spintm(n_sp)%valuesp(a,mm,e,ee))
                    end do

                    do m=1, nop
                    !\check{\tau} beta
                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(OtherspeciesId)%chtau_b(bb-nops,ee-nops,a-nop,jj,mm,m)* &
                            spintm(n_sp)%valuesp(m,mm,i,ee))
                    end do
                  end if

                end do
                
                do m=1, nop
                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    + 0.5*(Allinterspecies(speciesId)%tau(a-nop,bb-nops,m,mm)* &
                        spintm(n_sp)%valuesp(m,mm,i,jj))
                end do
              end do
              !print*,"nops 32"

              ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
              CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                  -Allspecies(speciesId)%HF_fs%values(a,a) -Allspecies(OtherspeciesId)%HF_fs%values(bb,bb))
                  
              Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj)

              ! \ddot{\tau} 
              Allinterspecies(speciesId)%tau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)
              ! \tilde{\tau} 
              Allinterspecies(speciesId)%intau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - 0.5*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)

              if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                do b=nop+1, noc
                  do cc=nops+1, nocs
                    do j=1, nop
                      do kk=1, nops
                        !  \check{\tau} triple excitation holy shit!!
                        Allinterspecies(speciesId)%chtau_a(a-nop,b-nop,cc-nops,i,j,kk) = Allspecies(speciesId)%Tssame(a-nop,i)* &
                          Allspecies(speciesId)%Tssame(b-nop,j)*Allspecies(OtherspeciesId)%Tssame(cc-nops,kk) &
                          + 0.5*( Allspecies(speciesId)%Tssame(a-nop,i)*Allinterspecies(speciesId)%Tdsame(b-nop,cc-nops,j,kk) &
                            + Allspecies(speciesId)%Tssame(b-nop,j)*Allinterspecies(speciesId)%Tdsame(a-nop,cc-nops,i,kk) &
                              + Allspecies(OtherspeciesId)%Tssame(cc-nops,kk)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) )
                      end do
                    end do
                  end do
                end do
              end if

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do aa=nops+1, nocs
                  do c=nop+1, noc
                    do ii=1, nops
                      do k=1, nop
                        !  \check{\tau} triple excitation holy shit!!
                        Allinterspecies(speciesId)%chtau_b(aa-nops,bb-nops,c-nop,ii,jj,k) = Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)* &
                          Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(speciesId)%Tssame(c-nop,k) &
                          + 0.5*( Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)*Allinterspecies(speciesId)%Tdsame(c-nop,bb-nops,k,jj) &
                            + Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allinterspecies(speciesId)%Tdsame(c-nop,aa-nops,k,ii) &
                              + Allspecies(speciesId)%Tssame(c-nop,k)*Allspecies(OtherspeciesId)%Tdsame(aa-nops,bb-nops,ii,jj) )
                      end do
                    end do
                  end do
                end do
              end if

              ! write(*,*) a,b,i,j,Allspecies(speciesId)%Tdsame(a,b,i,j),CCSDT1T2(speciesId)%Tabij_AB(a,b,i,j)
            end do
          end do
        end do
      end do
      print*, "nops 33"
      
  end subroutine CCSD_T2_AB

  !>
  ! @brief Manager of Coupled Cluster maths
  ! @author CAOM
  subroutine CCSD_run()
      implicit none

      integer :: n_intersp(10), i_counterID(10)
      integer :: num_species, counterID, finalID
      integer :: num_inter, num_i
      integer :: i, j, jj
      integer :: times_i
      integer :: aux_cont=0
      integer :: e_cont=0
      integer :: n_sp=0
      integer :: max, min
      real(8) :: intra=0
      real(8) :: convergence = 1.0_8
      real(8) :: convergence_int = 1.0_8
      real(8), allocatable :: e_same_ccd(:)
      real(8), allocatable :: e_diff_ccd(:)
      
      !Initialization of variables
      num_species = CoupledCluster_instance%num_species
      finalID = CoupledCluster_instance%finalID
      counterID = CoupledCluster_instance%counterID
      num_inter = CoupledCluster_instance%num_intersp
      times_i = CoupledCluster_instance%times_intersp
      CCSD_instance%max=0
      CCSD_instance%min=0

      if (convergence /= 1.0D-8) convergence = 1.0_8
      if (convergence_int /= 1.0D-8) convergence_int = 1.0_8

      ! allocate array for many results by species in CCSD
      if (allocated(CCSDT1T2)) deallocate(CCSDT1T2)
      allocate(CCSDT1T2(num_species))

      ! allocate array for many results by species in CCSD
      if (allocated(CCSDinit)) deallocate(CCSDinit)
      allocate(CCSDinit(num_species))

      ! allocate array for many results by species in CCSD
      if (allocated(CCSDloop)) deallocate(CCSDloop)
      allocate(CCSDloop(num_species))

      ! allocate array for results of inter-species in CCSD
      if (allocated(CCSDinter)) deallocate(CCSDinter)
      allocate(CCSDinter(num_inter))

      ! allocate array for results of intra-species in CCSD
      if (allocated(CCSD_instance%convergence_same)) deallocate(CCSD_instance%convergence_same)
      allocate(CCSD_instance%convergence_same(num_species))

      ! allocate array for results of intra-species in CCSD
      if (allocated(CCSD_instance%convergence_diff)) deallocate(CCSD_instance%convergence_diff)
      allocate(CCSD_instance%convergence_diff(num_inter))

      ! allocate array for many results by species in CCSD
      if (allocated(e_same_ccd)) deallocate(e_same_ccd)
      allocate(e_same_ccd(num_species))

      ! allocate array for many results by species in CCSD
      if (allocated(e_diff_ccd)) deallocate(e_diff_ccd)
      allocate(e_diff_ccd(num_inter))



      do i=counterID, num_species
        print*, "num_species: CCSD_constructor(): ", num_species, i
        call CCSD_constructor(i)
      end do

      do i=counterID, num_species!finalID
        print*, "num_species: CCSD_run(): ", num_species, i, finalID
        call CCSD_init(i)
        print*, "CCSD_init(i): ", i
      end do
        
      do i=counterID, finalID!num_species
        
        !If there are interspecies?
        if (times_i>0) then
          
          !search the principal species (speciesId) in all the options
          do j=1, times_i
            print*, "times_i", times_i
            !public to private            
            i_counterID(j) = CoupledCluster_instance%i_counterID(j)
            n_intersp(j) = CoupledCluster_instance%n_intersp(j)
            !Find the appropriate species in all the options
            if (i_counterID(j)==counterID) then
              !public variables used in the loop
              CCSD_instance%max=n_intersp(j)
              CCSD_instance%min=i_counterID(j)
              !for the order in call in spintm matrix
              aux_cont = n_sp + 1
              CCSD_instance%aux_cont = aux_cont
              !make all combinations with speciesId
              do jj=i_counterID(j)+1, n_intersp(j)
                n_sp = n_sp + 1
                CCSD_instance%cont = n_sp

                call CCSD_constructor_inter(i_counterID(j), jj)!, num_inter)
                call CCSD_init_inter(i_counterID(j), jj)!, num_inter)
                call CCSD_constructor_inter(jj, i_counterID(j))!, num_inter)
                call CCSD_init_inter(jj, i_counterID(j))!, num_inter)

                num_inter = num_inter + 1
              end do
            else
              !locating 
              n_sp = n_sp + (n_intersp(j) - i_counterID(j))
              CCSD_instance%cont = n_sp
            end if
          end do

          ! print*, "CCSD_loop(i, i+1): ", i
          ! call CCSD_loop(i, i+1) !fix i+1

        ! else

        !   call CCSD_loop(i)
        !   print*, "CCSD_loop(i): ", i
        
        end if
      end do

      do while (convergence >= 1.0D-8)

        do i=counterID, num_species
          e_same_ccd(i) = CoupledCluster_instance%CCSD_E_intra(i)
          call CCSD_same_species(i,e_same_ccd(i))
        end do

        if (times_i>0) then
          do i=counterID, finalID
            CCSD_instance%num_i = 1
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            max = CCSD_instance%max
            min = CCSD_instance%min
            num_i = CCSD_instance%num_i
            do jj=min+1, max
              e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
              call CCSD_diff_species(i,jj,e_diff_ccd(num_i))
              ! convergence_int = CCSD_instance%convergence_diff(num_i)
            end do
            do jj=min+1, max
              e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
              call CCSD_diff_species(jj,i,e_diff_ccd(num_i))
              ! convergence_int = CCSD_instance%convergence_diff(num_i)
            end do
          end do
        end if
        convergence = sum(CCSD_instance%convergence_same,dim=1) + & 
          sum(CCSD_instance%convergence_diff,dim=1)

      end do

      do i=counterID, finalID
        call CCSD_show(i)
        intra = intra + CoupledCluster_instance%CCSD_E_intra(i)
      end do

      CoupledCluster_instance%CCSD_once_Energy = sum(CoupledCluster_instance%CCSD_E_intra,dim=1)
      CoupledCluster_instance%CCSD_twice_Energy = sum(CoupledCluster_instance%CCSD_E_inter,dim=1)
      CoupledCluster_instance%CCSD_total_Energy = CoupledCluster_instance%CCSD_once_Energy + &
        CoupledCluster_instance%CCSD_twice_Energy

      ! CoupledCluster_instance%CCSD_ones_Energy = intra

      print*, "Total correction same species CCSD energy: ", CoupledCluster_instance%CCSD_once_Energy
      print*, "Total correction different species CCSD energy: ", CoupledCluster_instance%CCSD_twice_Energy
      print*, "Total CCSD energy: ", CoupledCluster_instance%CCSD_total_Energy &
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
      CCSD_instance%suma = CoupledCluster_instance%HF_energy + CoupledCluster_instance%MP2_EnergyCorr
      print*, "INFORMATION IN CCSD_constructor() TotalMP2_Energy: ", CCSD_instance%suma
      print*, "INFORMATION IN CCSD_constructor() Td: ", Allspecies(speciesId)%Tdsame(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() Dai: ", CCSDinit(speciesId)%Dai(1,3)
      print*, "INFORMATION IN CCSD_constructor() ttau: ",Allspecies(speciesId)%ttau(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() tau: ",Allspecies(speciesId)%tau(1,1,1,1)
      print*, "INFORMATION IN CCSD_constructor() ccsdE: ", CoupledCluster_instance%CCSD_E_intra(1)
      print*, "INFORMATION IN Tensor Tensor_index2: ", Tix2(2,1,3)
      

      !if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)

      !call CCSD_destructor()
      print*, "CCSD_show()x2"

  end subroutine CCSD_show
end module CCSD_
