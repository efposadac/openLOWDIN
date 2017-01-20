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
      real(8), allocatable :: convergence_same_amp(:)
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
          ! write(*,*) a,i,CCSDinit(speciesId)%Dai(a,i)
        end do
      end do
      do a=nop+1, noc
        do i=1, nop
          write(*,*) a,i,CCSDinit(speciesId)%Dai(a,i)
          print*, "HF_fs: ", Allspecies(speciesId)%HF_fs%values(i,i), Allspecies(speciesId)%HF_fs%values(a,a)
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
      integer noc, nocs, nop, nops, status
      
      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      
      write(*, *) " CCSD_loop_constructor: noc=", noc, " nop=", nop

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
      
      ! stop"ME"
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
      print*, " CCSD_T1T2_constructor(: ", nop, noc, speciesId
      
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
      ! print*, "Allinterspecies(speciesId)%Tdsame: (w:x:y:z): ", noc-nop,nocs-nops,nop,nops

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
        if (allocated(Allinterspecies(speciesId)%chtau_b)) deallocate (Allinterspecies(speciesId)%chtau_b)
        allocate(Allinterspecies(speciesId)%chtau_b(nocs-nops,nocs-nops,noc-nop,nops,nops,nop))
        Allinterspecies(speciesId)%chtau_b(:,:,:,:,:,:) = 0.0_8
      end if

      ! \ddot{\tau}
      if (allocated(Allinterspecies(speciesId)%tau)) deallocate (Allinterspecies(speciesId)%tau)
      allocate(Allinterspecies(speciesId)%tau(noc-nop,nocs-nops,nop,nops))
      Allinterspecies(speciesId)%tau(:,:,:,:) = 0.0_8

      ! \tilde{\tau}
      if (allocated(Allinterspecies(speciesId)%intau)) deallocate (Allinterspecies(speciesId)%intau)
      allocate(Allinterspecies(speciesId)%intau(noc-nop,nocs-nops,nop,nops))
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

      print*, "before loop inter: n_sp", n_sp, "nop: ", nop, "noc: ", noc, "nops: ", nops, "nocs: ", nocs
      print*, Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_a(1,1,1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_b(1,1,1,1,1,1)
      print*, Allinterspecies(speciesId)%tau(1,1,1,1)
      print*, Allinterspecies(speciesId)%intau(1,1,1,1)

      do a=nop+1, noc
        do i=1, nop
          do bb=nops+1, nocs
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
              ! print*, "spintm(n_sp)%valuesp(p,q,r,s): ", spintm(n_sp)%valuesp(p,q,r,s)
              ! print*, "Fii: ", Allspecies(speciesId)%HF_fs%values(i,i), "Fjjjj: ", Allspecies(OtherspeciesId)%HF_fs%values(jj,jj), &
              !   "Faa: ", -Allspecies(speciesId)%HF_fs%values(a,a), "Fbbbb: ", -Allspecies(OtherspeciesId)%HF_fs%values(bb,bb)
              
              ! print*, "2E-.2E+: ", n_sp
              Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                +( (spintm(n_sp)%valuesp(p,q,r,s)) /( Allspecies(speciesId)%HF_fs%values(i,i)+ &
                  Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) -Allspecies(speciesId)%HF_fs%values(a,a)- &
                    Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) ) )

              ! if (speciesId>OtherspeciesId) print*, "CCSD_init_inter Allinterspecies(speciesId)%Tdsame: ", &
              !   Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj)

              ! \ddot{\tau} 
              Allinterspecies(speciesId)%tau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)

              ! \tilde{\tau} 
              Allinterspecies(speciesId)%intau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                + 0.5*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj) !REVISAR

            end do
          end do

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
              do bb=nops+1, nocs
                do ii=1, nops
                  do jj=1, nops
                    !  \check{\tau} triple excitation holy shit!!
                    Allinterspecies(speciesId)%chtau_b(aa-nops,bb-nops,a-nop,ii,jj,i) = Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)* &
                      Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(speciesId)%Tssame(a-nop,i) &
                      + 0.5*( Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)*Allinterspecies(OtherspeciesId)%Tdsame(bb-nops,a-nop,jj,i) &
                        + Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allinterspecies(OtherspeciesId)%Tdsame(aa-nops,a-nop,ii,i) &
                          + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tdsame(aa-nops,bb-nops,ii,jj) )
                  end do
                end do
              end do
            end do
          end if

              ! write(*,*) Allinterspecies(speciesId)%Tdsame(a-nop,bb-nop,i,jj), "Tdsame"
              ! , Allinterspecies(speciesId)%HF_fs%values(a,a), "HF_fs%values", spintm(n_sp)%valuesp(i,jj,a,bb), "spintm"
        end do
        ! print*, "nops 1 - 5: nops+1 y nocs: ", bb, nocs, "nop+1 y noc: ", a, noc
        ! print*, "nops 1 - 5: i y nop: ", i, nop, "jj nops: ", jj, nops
        ! print*, Allinterspecies(speciesId)%Tdsame!(a-nop,bb-nops,i,jj)
      end do

      print*, Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_a(1,1,1,1,1,1)
      !print*, Allinterspecies(speciesId)%chtau_b(1,1,1,1,1,1)
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

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(speciesId)%Wakic)) deallocate (CCSDinter(speciesId)%Wakic)
      allocate(CCSDinter(speciesId)%Wakic(noc-nop,nops,nop,nocs-nops))
      CCSDinter(speciesId)%Wakic=0.0_8

      !Initial guess: Information of another species T2: alpha-beta
      if (allocated(CCSDinter(OtherspeciesId)%Wbkjc)) deallocate (CCSDinter(OtherspeciesId)%Wbkjc)
      allocate(CCSDinter(OtherspeciesId)%Wbkjc(nop,nocs-nops,noc-nop,nops))
      CCSDinter(OtherspeciesId)%Wbkjc=0.0_8
      print*, "nops 14"

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
      integer :: p, q, r, s

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

              if (speciesId<OtherspeciesId) then
                p=a
                q=mm
                r=e
                s=ee
              else
                p=mm
                q=a
                r=ee
                s=e
              end if
              CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                + ( 0.25*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)* & 
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 4-5
              ! print*, "test OtherspeciesId: p=",p,", q=",q,", r=",r,", s=", s, " n_sp=", n_sp
              ! Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)

              do m=1, nop

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=e
                end if

                CCSDloop(speciesId)%Fac(a-nop,e-nop) = CCSDloop(speciesId)%Fac(a-nop,e-nop) &
                  - ( 0.25*Allinterspecies(speciesId)%intau(a-nop,ee-nops,m,mm)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 4-6
                  ! write(*,*) a,e,CCSDloop(speciesId)%Fac(a,e)
              end do
            end do
          end do

        end do
      end do
      ! print*, "nops 6"

      ! CCSDloop(speciesId)%Fki(m,i) 2
      do m=1, nop
        do i=1, nop
             
          do mm=1, nops
            do ee=nops+1, nocs

              if (speciesId<OtherspeciesId) then
                p=m
                q=mm
                r=i
                s=ee
              else
                p=mm
                q=m
                r=ee
                s=i
              end if

              CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                 + ( 0.25*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 5-5
            
              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=e
                end if

                CCSDloop(speciesId)%Fki(m,i) = CCSDloop(speciesId)%Fki(m,i) &
                  + ( 0.25*Allinterspecies(speciesId)%intau(e-nop,ee-nops,i,mm)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 5-6
                ! write(*,*) a,e,CCSDloop(speciesId)%Fki(m,i)
              end do
            end do
          end do

        end do
      end do
      ! print*, "nops 7"

      ! CCSDloop(speciesId)%Fkc_aa(m,e-nop) 3
      do m=1, nop
        do e=nop+1, noc

          do mm=1, nops
            do ee=nops+1, nocs

              if (speciesId<OtherspeciesId) then
                p=m
                q=mm
                r=e
                s=ee
              else
                p=mm
                q=m
                r=ee
                s=e
              end if

              CCSDloop(speciesId)%Fkc_aa(m,e-nop) = CCSDloop(speciesId)%Fkc_aa(m,e-nop) &
                + ( 0.25*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)* &
                    spintm(n_sp)%valuesp(p,q,r,s))
              ! write(*,*) a,e,CCSDloop(speciesId)%Fkc_aa(m,e-nop)
            end do
          end do

        end do
      end do
      ! print*, "nops 8"

      ! ! CCSDinter(speciesId)%Fkc_aba(m,e-nop) 4
      do m=1, nop
        do e=nop+1, noc
          do aa=nops+1, nocs
            do ii=1, nops

              if (speciesId<OtherspeciesId) then
                p=m
                q=aa
                r=e
                s=ii
              else
                p=aa
                q=m
                r=ii
                s=e
              end if

              !REVISAR 01/10/2017
              ! CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
              !   - (0.5*spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-2

              ! print*, "%Fkc_aba: ", CCSDinter(speciesId)%Fkc_aba(m,e-nop)

              do ee=nops+1, nocs

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=aa
                  r=e
                  s=ee
                else
                  p=aa
                  q=m
                  r=ee
                  s=e
                end if

                CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                  -(0.25* Allspecies(OtherspeciesId)%Tssame(ee-nops,ii)*spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-3

                if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                  
                  do mm=1, nops

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=e
                    end if

                    CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                      -(0.25* Allspecies(OtherspeciesId)%Tdsame(aa-nops,ee-nops,ii,mm)* &
                          spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-5

                    CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                      +(0.25* Allspecies(OtherspeciesId)%Tssame(aa-nops,mm)* &
                          Allspecies(OtherspeciesId)%Tssame(ee-nops,ii) &
                          *spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-7
                  end do
                end if
              end do

              do mm=1, nops

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=e
                  s=ii
                else
                  p=mm
                  q=m
                  r=ii
                  s=e
                end if

                CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                  +(0.25* Allspecies(OtherspeciesId)%Tssame(aa-nops,mm)*spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-4
              end do

            end do
          end do

          do a=nop+1, noc
            do i=1, nop

              !REVISAR
              CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                + spints(speciesId)%valuesp(m,a,i,e) !eq 7-1
              
              do mm=1, nops
                do ee=nops+1, nocs

                  if (speciesId<OtherspeciesId) then
                    p=m
                    q=mm
                    r=e
                    s=ee
                  else
                    p=mm
                    q=m
                    r=ee
                    s=e
                  end if

                  CCSDinter(speciesId)%Fkc_aba(m,e-nop) = CCSDinter(speciesId)%Fkc_aba(m,e-nop) &
                    -(0.125* Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* & 
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 7-6
                end do
              end do
            end do
          end do

          ! write(*,*) a,e,CCSDinter(num_inter)%Fkc_aba(m,e-nop)
        end do
      end do
      ! print*, "nops 9: "

      ! CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) 5
      do mm=1, nops
        do ee=nops+1, nocs

          if(nops>=2) then ! kind of interaction just for two or more particles of the another species

            CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
              + ( 0.5*Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)) !eq 8-1

            do nn=1, nops
              do ff=nops+1, nocs

                CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
                  + ( 0.5*Allspecies(OtherspeciesId)%Tssame(ff-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,ee,ff)) !eq 8-2
              end do
            end do
          end if

          do m=1, nop
            do e=nop+1, noc

              if (speciesId<OtherspeciesId) then
                p=m
                q=mm
                r=e
                s=ee
              else
                p=mm
                q=m
                r=ee
                s=e
              end if

              CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) = CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops) &
                + ( 0.125*Allspecies(speciesId)%Tssame(e-nop,m)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 8-3
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
            +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e) !eq 5-1

          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do m=1, nop
              do f=nop+1, noc
            
                CCSDinter(speciesId)%Faca(a-nop,e-nop) = CCSDinter(speciesId)%Faca(a-nop,e-nop) &
                  - ( Allspecies(speciesId)%Tssame(a-nop,m)* &
                      spints(speciesId)%valuesp(m,a,e,f)) !eq 5-2
              end do
            end do
          end if
        ! write(*,*) a,e,CCSDinter(speciesId)%Faca(a-nop,e-nop)
        end do
      end do
      ! print*, "nops 15"

      ! CCSDinter(speciesId)%Fkia
      do m=1, nop
        do i=1, nop

          CCSDinter(speciesId)%Fkia(m,i) = CCSDinter(speciesId)%Fkia(m,i) &
            +(1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i) !eq 7-1

          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
            do n=1, nop
              do e=nop+1, noc
 
                CCSDinter(speciesId)%Fkia(m,i) = CCSDinter(speciesId)%Fkia(m,i) &
                  - ( Allspecies(speciesId)%Tssame(e-nop,n)* &
                      spints(speciesId)%valuesp(m,n,i,e)) !eq 7-2
              end do
            end do
          end if
          ! write(*,*) a,e,CCSDinter(speciesId)%Fkia(m,i)
        end do
      end do
      ! print*, "nops 16"

      ! CCSDinter(OtherspeciesId)%Fbcb
      do bb=nops+1, nocs
        do ee=nops+1, nocs

          CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) = CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) &
            +(1 - logic2dbl(bb==ee))*Allspecies(OtherspeciesId)%HF_fs%values(bb,ee) !eq 6-1

          if (nops>=2) then ! kind of interaction just for two or more particles of the another species
            do mm=1, nops
              do ff=nops+1, nocs
            
                CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) = CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops) &
                  - ( Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      spints(OtherspeciesId)%valuesp(mm,bb,ee,ff)) !eq 6-2
              end do
            end do
          end if
        ! write(*,*) bb,ee,CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops)
        end do
      end do
      ! print*, "nops 17"

      ! CCSDinter(OtherspeciesId)%Fkjb
      do mm=1, nops
        do jj=1, nops

          CCSDinter(OtherspeciesId)%Fkjb(mm,jj) = CCSDinter(OtherspeciesId)%Fkjb(mm,jj) &
            +(1 - logic2dbl(mm==jj))*Allspecies(OtherspeciesId)%HF_fs%values(mm,jj) !eq 8-1

          if (nops>=2) then ! kind of interaction just for two or more particles of the another species          
            do nn=1, nops
              do ee=nops+1, nocs

                CCSDinter(OtherspeciesId)%Fkjb(mm,jj) = CCSDinter(OtherspeciesId)%Fkjb(mm,jj) &
                  - ( Allspecies(OtherspeciesId)%Tssame(ee-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,jj,ee)) !eq 8-2
              end do
            end do
          end if
        ! write(*,*) a,e,CCSDinter(OtherspeciesId)%Fkjb(mm,jj)
        end do
      end do
      print*, "nops 18"
      ! print*, "F_two_species Allinterspecies(speciesId)%Tdsame", Allinterspecies(OtherspeciesId)%Tdsame(ff-nops,b-nop,nn,j)

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
      !!$OMP PARALLEL
      !!$OMP DO
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
      !!$OMP END DO

      !!$OMP DO
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
      !!$OMP END DO

      !!$OMP DO
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
      !!$OMP END DO
      !!$OMP END PARALLEL

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
      integer :: p, q, r, s, pp, qq, rr, ss

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

              if (speciesId<OtherspeciesId) then
                p=b
                q=mm
                r=j
                s=ee
              else
                p=mm
                q=b
                r=ee
                s=j
              end if

              CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                + (0.5*spintm(n_sp)%valuesp(p,q,r,s)) !eq 9-1
              
              do m=1, nop

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=j
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=j
                end if

                CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                  - 0.75*(Allspecies(speciesId)%Tssame(b-nop,m)*spintm(n_sp)%valuesp(p,q,r,s)) !eq 9-2
              end do

              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=b
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=b
                  r=ee
                  s=e
                end if

                CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                  + (0.75*Allspecies(speciesId)%Tssame(e-nop,j)*spintm(n_sp)%valuesp(p,q,r,s)) !eq 9-3
              end do

            end do
          end do
                  
          if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
          
            do m=1, nop
              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=e
                end if

                CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                  - ((0.5*Allspecies(speciesId)%Tdsame(e-nop,b-nop,j,m)) & !eq 9-4
                    + (Allspecies(speciesId)%Tssame(e-nop,j)*Allspecies(speciesId)%Tssame(b-nop,m)))* & !eq 9-5
                      spintm(n_sp)%valuesp(p,q,r,s)
              end do
            end do
          end if

          if (nops>=2) then ! kind of interaction just for two or more particles of the another species

            do ff=nops+1, nocs
              do nn=1, nops
                  
                CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) = CCSDinter(speciesId)%Wkkcc(b-nop,mm,j,ee-nops) &
                  + (0.25*Allinterspecies(OtherspeciesId)%Tdsame(ff-nops,b-nop,nn,j)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,ee,ff)) !eq 9-6
              end do
            end do
                ! if (speciesId<OtherspeciesId) print*, "W 2 Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(OtherspeciesId)%Tdsame(ff-nops,b-nop,nn,j)
          end if
              ! write(*,*) m,n,i,j,CCSDloop(speciesId)%Wklij(m,n,i,j)
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
      integer :: p, q, r, s

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
        do bb=nops+1, nocs
          do i=1, nop
            do jj=1, nops

              if (speciesId<OtherspeciesId) then
                p=m
                q=bb
                r=i
                s=jj
              else
                p=bb
                q=m
                r=jj
                s=i
              end if

              CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                + spintm(n_sp)%valuesp(p,q,r,s) !eq 1-1

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=i
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=i
                    end if

                    CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                      + (0.5*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 1-2
                  end do
                end do
              end if

              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=bb
                  r=e
                  s=jj
                else
                  p=bb
                  q=m
                  r=jj
                  s=e
                end if

                CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) = CCSDinter(speciesId)%Waka(m,bb-nops,i,jj) &
                  + (0.5*Allspecies(speciesId)%Tssame(e-nop,i)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 1-3
              end do
              ! write(*,*) m,bb-nops,i,jj,CCSDinter(speciesId)%Waka(m,bb-nops,i,jj)
            end do
          end do
        end do
      end do
      ! print*, "nops 19"

      ! CCSDinter(speciesId)%Wcia
      do a=nop+1, noc
        do bb=1+nops, nocs
          do e=nop+1, noc
            do jj=1, nops

              if (speciesId<OtherspeciesId) then
                p=a
                q=bb
                r=e
                s=jj
              else
                p=bb
                q=a
                r=jj
                s=e
              end if

              CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                + spintm(n_sp)%valuesp(p,q,r,s) !eq 2-1

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops

                    if (speciesId<OtherspeciesId) then
                      p=a
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=a
                      r=ee
                      s=e
                    end if

                    CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                      + (0.5*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 2-2
                  end do
                end do
              end if

              do m=1, nop

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=bb
                  r=e
                  s=jj
                else
                  p=bb
                  q=m
                  r=jj
                  s=e
                end if

                CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) = CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj) &
                  + (0.5*Allspecies(speciesId)%Tssame(a-nop,m)* &
                       spintm(n_sp)%valuesp(p,q,r,s)) !eq 2-3
              end do
              ! write(*,*) a-nop,bb-nops,e-nop,jj,CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj)
            end do
          end do
        end do
      end do
      ! print*, "nops 20"

      ! CCSDinter(speciesId)%Wbkb
      do a=nop+1, noc
        do mm=1, nops
          do i=1, nop
            do jj=1, nops

              if (speciesId<OtherspeciesId) then
                p=a
                q=mm
                r=i
                s=jj
              else
                p=mm
                q=a
                r=jj
                s=i
              end if

              CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                + spintm(n_sp)%valuesp(p,q,r,s) !eq 3-1

              do e=nop+1, noc
                do m=1, nop

                  if (speciesId<OtherspeciesId) then
                    p=m
                    q=mm
                    r=e
                    s=jj
                  else
                    p=mm
                    q=m
                    r=jj
                    s=e
                  end if

                  CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                    + (0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 3-2
                end do
              end do

              do ee=nops+1, nocs
              
                if (speciesId<OtherspeciesId) then
                  p=a
                  q=mm
                  r=i
                  s=ee
                else
                  p=mm
                  q=a
                  r=ee
                  s=i
                end if

                CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) = CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 3-3
              end do
              ! write(*,*) a-nop,mm,i,jj,CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj)
            end do
          end do
        end do
      end do
      ! print*, "nops 21"

      ! CCSDinter(speciesId)%Wcjb
      do a=nop+1, noc
        do bb=nops+1, nocs
          do i=1, nop
            do ee=nops+1, nocs

              if (speciesId<OtherspeciesId) then
                p=a
                q=bb
                r=i
                s=ee
              else
                p=bb
                q=a
                r=ee
                s=i
              end if

              CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                + spintm(n_sp)%valuesp(p,q,r,s) !eq 4-1

              do e=nop+1, noc
                do m=1, nop

                  if (speciesId<OtherspeciesId) then
                    p=m
                    q=bb
                    r=e
                    s=ee
                  else
                    p=bb
                    q=m
                    r=ee
                    s=e
                  end if

                  CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                    + (0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 4-2
                end do
              end do

              do mm=1, nops

              if (speciesId<OtherspeciesId) then
                p=a
                q=mm
                r=i
                s=ee
              else
                p=mm
                q=a
                r=ee
                s=i
              end if

                CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) = CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      spintm(n_sp)%valuesp(p,q,r,s)) !eq 4-3
              end do
              ! write(*,*) a-nop,bb-nops,i,ee-nops,CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops)
            end do
          end do
        end do
      end do
      ! print*, "nops 22"

      ! CCSDinter(speciesId)%Wakic_a
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        do a=nop+1, noc
          do m=1, nop
            do i=1, nop
              do e=nop+1, noc

                CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (Allspecies(speciesId)%HF_fs%values(m,e)* &
                        Allspecies(speciesId)%Tssame(a-nop,i)) !eq 9-3

                CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                  + spints(speciesId)%valuesp(a,m,i,e) !eq 9-1

                do n=1, nop
                  CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (Allspecies(speciesId)%Tssame(a-nop,n)* &
                      spints(speciesId)%valuesp(m,n,i,e)) !eq 9-4
                end do

                do f=nop+1, noc
                  CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                    + (Allspecies(speciesId)%Tssame(f-nop,i)* &
                      spints(speciesId)%valuesp(m,a,e,f)) !eq 9-5
                end do

                do ee=nops+1, nocs
                  do mm=1, nops

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=e
                    end if

                    CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) = CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop) &
                      + (0.125*Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                          spintm(n_sp)%valuesp(p,q,r,s)) !eq 9-2
                  end do
                end do
                ! if (speciesId>OtherspeciesId) print*, "W T2 Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)
                ! write(*,*) a-nop,m,i,e-nop,CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop)
              end do
            end do
          end do
        end do
      end if
      ! print*, "nops 23"

      ! CCSDinter(speciesId)%Wbkjc_b 
      if (nops>=2) then ! kind of interaction just for two or more particles of the another species
        do bb=nops+1, nocs
          do mm=1, nops
            do jj=1, nops
              do ee=nops+1, nocs

                CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)* &
                        Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)) !eq 10-3

                CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                  + spints(OtherspeciesId)%valuesp(bb,mm,jj,ee) !eq 10-1

                do nn=1, nops
                  CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (Allspecies(OtherspeciesId)%Tssame(bb-nops,nn)* &
                      spints(OtherspeciesId)%valuesp(mm,nn,jj,ee)) !eq 10-4
                end do

                do ff=nops+1, nocs
                  CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                    + (Allspecies(OtherspeciesId)%Tssame(ff-nops,jj)* &
                      spints(OtherspeciesId)%valuesp(mm,bb,ee,ff)) !eq 10-5
                end do

                do e=nop+1, noc
                  do m=1, nop

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=e
                    end if

                    CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) = CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops) &
                      + (0.125*Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)* &
                          spintm(n_sp)%valuesp(p,q,r,s)) !eq 10-2
                    ! if (speciesId>OtherspeciesId) print*, "W T2 Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)
                  end do
                end do
                ! write(*,*) bb-nops,mm,jj,ee-nops,CCSDinter(speciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops)
              end do
            end do
          end do
        end do
      end if
      ! print*, "nops 24"

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
                        (Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj) & !eq 11-1
                          + (Allspecies(speciesId)%Tssame(e-nop,m)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj))) !eq 11-2
                  end do
                end do
                ! if (speciesId>OtherspeciesId) print*, "W T2 Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)
                ! write(*,*) m,n,e-nop,f-nop,CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop)
              end do
            end do
          end do
        end do
      end if
      ! print*, "nops 25"

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
                        (Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm) & !eq 12-1
                          + (Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm))) !eq 12-2
                  end do
                end do
                ! if (speciesId>OtherspeciesId) print*, "Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)
                ! write(*,*) mm,nn,ee-nops,ff-nops,CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops)
              end do
            end do
          end do
        end do
      end if
      ! print*, "nops 26"

      ! CCSDinter(speciesId)%Wakic
      do a=nop+1, noc
        do mm=1, nops
          do i=1, nop
            do ee=nops+1, nocs

              if (speciesId<OtherspeciesId) then
                p=a
                q=mm
                r=i
                s=ee
              else
                p=mm
                q=a
                r=ee
                s=i
              end if

              CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                + spintm(n_sp)%valuesp(p,q,r,s) & !eq 14-1
                  + (Allspecies(OtherspeciesId)%HF_fs%values(mm,ee)* &
                    Allspecies(speciesId)%Tssame(a-nop,i)) !eq 14-2

              do m=1, nop

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=i
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=i
                end if

                CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                  - (0.5*Allspecies(speciesId)%Tssame(a-nop,m)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 14-3
              end do

              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=a
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=a
                  r=ee
                  s=e
                end if

                CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                  + (0.5*Allspecies(speciesId)%Tssame(e-nop,i)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 14-4
              end do

              if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
                do e=nop+1, noc
                  do m=1, nop

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=e
                    end if

                    CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) = CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops) &
                      + ( Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(speciesId)%Tssame(e-nop,m) & !eq 14-5
                        + (0.5*Allspecies(speciesId)%tau(a-nop,e-nop,i,m)) & !eq 14-6
                          - (0.25*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,i)))* & !eq 14-7
                           spintm(n_sp)%valuesp(p,q,r,s)
                  end do
                end do
              end if
              ! write(*,*) a-nop,mm,i,ee-nops,CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops)
            end do
          end do
        end do
      end do
      ! print*, "nops 27"

      ! CCSDinter(speciesId)%Wbkjc
      do m=1, nop
        do bb=nops+1, nocs
          do e=nop+1, noc
            do jj=1, nops

              if (speciesId<OtherspeciesId) then
                p=m
                q=bb
                r=e
                s=jj
              else
                p=bb
                q=m
                r=jj
                s=e
              end if

              CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                + spintm(n_sp)%valuesp(p,q,r,s) & !eq 13-1
                  + (Allspecies(speciesId)%HF_fs%values(m,e)* &
                    Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)) !eq 13-2

              do mm=1, nops

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=e
                  s=jj
                else
                  p=mm
                  q=m
                  r=jj
                  s=e
                end if

                CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                  - (0.5*Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 13-3
              end do

              do ee=nops+1, nocs

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=bb
                  r=e
                  s=ee
                else
                  p=bb
                  q=m
                  r=ee
                  s=e
                end if

                CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                  + (0.5*Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                    spintm(n_sp)%valuesp(p,q,r,s)) !eq 13-4
              end do

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ee=nops+1, nocs
                  do mm=1, nops

                    if (speciesId<OtherspeciesId) then
                      p=m
                      q=mm
                      r=e
                      s=ee
                    else
                      p=mm
                      q=m
                      r=ee
                      s=e
                    end if

                    CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) = CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj) &
                      +((Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(OtherspeciesId)%Tssame(ee-nops,mm)) & !eq 13-5
                        + (0.5*Allspecies(OtherspeciesId)%tau(bb-nops,ee-nops,jj,mm)) & !eq 13-6
                          - (0.25*Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,mm,jj)) )* & !eq 13-7
                          spintm(n_sp)%valuesp(p,q,r,s)
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
  subroutine CCSD_same_species(speciesId, e_ccsd, v_ampl)
      implicit none

      integer, intent(in) :: speciesId
      real(8), intent(in) :: e_ccsd
      real(8), intent(in) :: v_ampl

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
      real(8) :: prev_ampl
      real(8) :: tmp_ccsdE
      real(8) :: prev_ccsdE_int
      real(8) :: tmp_ccsdE_int
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

        !**change position of do while
        !do i=1, num_species
          !all intra CCSD
        !end do
        !do i=min, max
        ! all inter-species
        !end do
      prev_ccsdE = e_ccsd
      prev_ampl = v_ampl

      print*, "speciesId CCSD_loop: ", speciesId

      call CCSD_T1T2_constructor(speciesId)
      call CCSD_loop_constructor(speciesId)
      ! stop"FER CCSD_same_species"
      !intermediates loop for:
      !singles excitations
      print*, "F_onespecies_intermediates(): "
      call F_onespecies_intermediates(speciesId)
      !doubles excitations
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        print*, "W_onespecies_intermediates(): "
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
                 auxtdsame = auxtdsame + 0.25*spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,j) 
                 auxtssame = auxtssame + 0.5*spints(speciesId)%valuesp(i,j,a,b)*Allspecies(speciesId)%Tssame(a-nop,i)* Allspecies(speciesId)%Tssame(b-nop,j)

              end do
            end do
          end if
        end do
      end do
      print *, "tdsame intra", auxtdsame
      print *, "tssame intra", auxtssame
      ccsdE = tmp_ccsdE

      ampl = sum(Allspecies(speciesId)%Tdsame) + sum(Allspecies(speciesId)%Tssame)

      conver_amplitude = abs(ampl-prev_ampl)

      !change in values for the intra-species loop
      convergence = abs( ccsdE - prev_ccsdE )
      print*, "Otherspecies Tssame: ", Allspecies(speciesId)%Tssame
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
      CoupledCluster_instance%CCSD_A_intra(speciesId) = ampl
      CCSD_instance%convergence_same(speciesId) = convergence
      CCSD_instance%convergence_same_amp(speciesId) = conver_amplitude
      
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
      integer :: p, q, r, s
      integer :: num_inter
      real(8) :: prev_ccsdE_int
      real(8) :: tmp_ccsdE_int
      real(8) :: convergence = 1.0_8
      real(8) :: ccsdE_int=0.0_8
      real(8) :: auxtdsame = 0.0_8
      real(8) :: auxtssame = 0.0_8

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

      print*, "speciesId: ", speciesId

      num_inter = CCSD_instance%num_i
      e_cont = CCSD_instance%e_cont

      ! call CCSD_constructor_inter(min, jj)
      print*, "OtherspeciesId: ", OtherspeciesId
      ! if (speciesId>OtherspeciesId) print*, "before CCSD_T2AB_constructor: ", Allinterspecies(speciesId)%Tdsame
      call CCSD_T2AB_constructor(speciesId, OtherspeciesId, num_inter)
      ! if (speciesId>OtherspeciesId) print*, "before CCSD_loop_constructor_inter: ", Allinterspecies(speciesId)%Tdsame
      call CCSD_loop_constructor_inter(speciesId, OtherspeciesId, num_inter)

      !intermediates loop for:
      !If there are interspecies?
      print*, "ciclo: OtherspeciesId: ", OtherspeciesId, "speciesId: ", speciesId
      ! if (speciesId>OtherspeciesId) print*, "before F_twospecies: ", Allinterspecies(speciesId)%Tdsame
      call F_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      ! if (speciesId>OtherspeciesId) print*, "before W_twospecies: ", Allinterspecies(speciesId)%Tdsame
      call W_twospecies_intermediates(speciesId, OtherspeciesId, num_inter)
      ! if (speciesId>OtherspeciesId) print*, "before F_T2_AB: ", Allinterspecies(speciesId)%Tdsame
      call F_T2_AB(speciesId, OtherspeciesId, num_inter)
      ! if (speciesId>OtherspeciesId) print*, "between F_T2_AB and W_T2_AB: ", Allinterspecies(speciesId)%Tdsame
      call W_T2_AB(speciesId, OtherspeciesId, num_inter)
        
      ! Resolve CCSD equation of energy

      tmp_ccsdE_int=0.0_8
      print*, "energy CCSD-APMO: ", speciesId+1, max
      print*, "e_cont: ", e_cont
      auxtdsame = 0 
      auxtssame = 0 
      do i=1, nop
        do a=nop+1, noc
          do ii=1, nops
            do aa=nops+1, nocs

              if (speciesId<OtherspeciesId) then
                p=i
                q=ii
                r=a
                s=aa
              else
                p=ii
                q=i
                r=aa
                s=a
              end if

              tmp_ccsdE_int = tmp_ccsdE_int + (0.5*spintm(e_cont)%valuesp(p,q,r,s)* &
                  Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) ) &
                    + (0.5*spintm(e_cont)%valuesp(p,q,r,s)*Allspecies(speciesId)%Tssame(a-nop,i)* &
                      Allspecies(OtherspeciesId)%Tssame(aa-nops,ii))
              auxtdsame = auxtdsame + (0.5*spintm(e_cont)%valuesp(p,q,r,s)* Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) )
              auxtssame = auxtssame + (0.5*spintm(e_cont)%valuesp(p,q,r,s)*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(aa-nops,ii))

              ! print*, "spintm(e_cont): ", spintm(e_cont)%valuesp(p,q,r,s)
              ! if (speciesId>OtherspeciesId) print*, "CCSD diff Allinterspecies(speciesId)%Tdsame: ", &
              !   Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii)
              ! print*, "Allspecies(speciesId)%Tssame: ", Allspecies(speciesId)%Tssame(a-nop,i)
              ! print*, "Allspecies(OtherspeciesId)%Tssame: ", Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)
            end do
          end do
        end do
      end do
      print *, "speciesId", speciesId
      print *, "tdsame E inter", auxtdsame
      print *, "tssame E inter", auxtssame
      print *, "Tdsame: ", Allinterspecies(speciesId)%Tdsame(1,1,1,1)
      print *, "intau: ", Allinterspecies(speciesId)%intau(1,1,1,1)
      print *, "tau: ", Allinterspecies(speciesId)%tau(1,1,1,1)
      print *, "Tssame: ", Allspecies(speciesId)%Tssame(1,1)
      print *, "Fac: ", CCSDloop(speciesId)%Fac(1,1)
      print *, "Fki: ", CCSDloop(speciesId)%Fki(1,1)
      print *, "Fkc_aa: ", CCSDloop(speciesId)%Fkc_aa(1,1)
      print *, "Fkc_aba: ", CCSDinter(speciesId)%Fkc_aba(1,1)
      print *, "Fkca_ab: ", CCSDinter(OtherspeciesId)%Fkca_ab(1,1)
      print *, "Faca: ", CCSDinter(speciesId)%Faca(1,1)
      print *, "Fkia: ", CCSDinter(speciesId)%Fkia(1,1)
      print *, "Fbcb: ", CCSDinter(OtherspeciesId)%Fbcb(1,1)
      print *, "Fkjb: ", CCSDinter(OtherspeciesId)%Fkjb(1,1)
      print *, "Waka: ", CCSDinter(speciesId)%Waka(1,1,1,1)
      print *, "Wcia: ", CCSDinter(speciesId)%Wcia(1,1,1,1)
      print *, "Wbkb: ", CCSDinter(speciesId)%Wbkb(1,1,1,1)
      print *, "Wcjb: ", CCSDinter(speciesId)%Wcjb(1,1,1,1)
      print *, "Wakic: ", CCSDinter(speciesId)%Wakic(1,1,1,1)
      print *, "Wbkjc: ", CCSDinter(OtherspeciesId)%Wbkjc(1,1,1,1)
      print *, "spintm: ", spintm(e_cont)%valuesp(1,1,1,1)


      print*, "energy CCSD-APMO"

      ccsdE_int = tmp_ccsdE_int

      !change in values for the inter-species loop
      convergence = abs( ccsdE_int - prev_ccsdE_int )

      write (*,*) ccsdE_int, "CCSD Energy inter-species ", prev_ccsdE_int, "previous Energy inter-species" 
      write (*,*) convergence, "Convergence inter-species" 

      ! ! T1 and T2 equation energies for interspecies
      ! print*, "before CCSD_T1_inter()", speciesId, OtherspeciesId
      ! call CCSD_T1_inter(speciesId, OtherspeciesId, num_inter)
      ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T1 Allinterspecies(speciesId)%Tdsame: ", &
      ! !   Allinterspecies(speciesId)%Tdsame
      ! print*, "before CCSD_T2_inter()", speciesId, OtherspeciesId
      ! call CCSD_T2_inter(speciesId, OtherspeciesId, num_inter)
      ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2 Allinterspecies(speciesId)%Tdsame: ", &
      ! !   Allinterspecies(speciesId)%Tdsame
      ! print*, "CCSD_T2_AB()"
      ! call CCSD_T2_AB(speciesId, OtherspeciesId, num_inter)
      ! ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2AB Allinterspecies(speciesId)%Tdsame: ", &
      ! !   Allinterspecies(speciesId)%Tdsame
      ! print*, "num_inter: ", num_inter
  
      if (convergence > 100) then 
        stop "Error: There are not convergence. The differences between energies is more than 100 eV"
      end if


      CoupledCluster_instance%CCSD_E_inter(num_inter) = ccsdE_int
      CCSD_instance%convergence_diff(num_inter) = convergence
      
      num_inter = num_inter + 1!final
      CCSD_instance%num_i = num_inter
      CCSD_instance%cont = CCSD_instance%cont + 1!final

      e_cont = e_cont + 1!final
      CCSD_instance%e_cont = e_cont

  end subroutine CCSD_diff_species

  !>
  ! @brief Make convergence of amplitude and energy equations for Coupled Cluster
  ! @author CAOM
  subroutine CCSD_diff_species_energy(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId

      integer :: num_inter
      integer :: nop, noc

      noc = Allspecies(speciesId)%noc
      nop = Allspecies(speciesId)%nop

      num_inter = CCSD_instance%num_i

      ! T1 and T2 equation energies for interspecies
      print*, "before CCSD_T1_inter()", speciesId, OtherspeciesId
      call CCSD_T1_inter(speciesId, OtherspeciesId, num_inter)
      ! if (speciesId>OtherspeciesId) print*, "after CCSD_T1 Allinterspecies(speciesId)%Tdsame: ", &
      !   Allinterspecies(speciesId)%Tdsame
      print*, "before CCSD_T2_inter()", speciesId, OtherspeciesId
      if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
        call CCSD_T2_inter(speciesId, OtherspeciesId, num_inter)
      end if
      ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2 Allinterspecies(speciesId)%Tdsame: ", &
      !   Allinterspecies(speciesId)%Tdsame

  end subroutine CCSD_diff_species_energy


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
            + Allspecies(speciesId)%HF_fs%values(i,a) ! eq1 1
          do e=nop+1, noc
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + Allspecies(speciesId)%Tssame(e-nop,i)*CCSDloop(speciesId)%Fac(a-nop,e-nop) !eq 4 1
          end do
          do m=1, nop
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
              + (-Allspecies(speciesId)%Tssame(a-nop,m)*CCSDloop(speciesId)%Fki(m,i)) !eq 5 2
            if (nop>=2) then ! kind of interaction just for two or more particles of the principal species 
              do e=nop+1, noc
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)*CCSDloop(speciesId)%Fkc_aa(m,e-nop) !eq 6 3
                do f=nop+1, noc
                  CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                    + (-0.5*Allspecies(speciesId)%Tdsame(e-nop,f-nop,i,m)*spints(speciesId)%valuesp(m,a,e,f)) !eq 3 7
                end do
                do n=1, nop
                  CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                    + (-0.5*Allspecies(speciesId)%Tdsame(a-nop,e-nop,m,n)*spints(speciesId)%valuesp(n,m,e,i)) !eq 2 6
                end do
              end do
            end if
          end do
          if ((nop>=2)) then !  .and. (times_i==0) kind of interaction just for two or more particles of the principal species 
            !if times_i/=0 then the product below will be calculated in Fkc_aba in subroutine CCSD_T1_inter()
            do m=1,nop !n
              do e=nop+1, noc !f
                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + (-Allspecies(speciesId)%Tssame(e-nop,m)*spints(speciesId)%valuesp(m,a,i,e)) !eq 7 5 !n,a,i,f
              end do
            end do
          end if
          if (times_i==0) then
            CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit(speciesId)%Dai(a,i)
            Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
          end if
          print*, "Tai: ", CCSDT1T2(speciesId)%Tai(a-nop,i)
          print*, "Dai: ", CCSDinit(speciesId)%Dai(a,i)
          ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
        end do
      end do
      print*, "CCSD_T1"
      ! print*, "Tai: ", CCSDT1T2(speciesId)%Tai
      ! print*, "Dai: ", CCSDinit(speciesId)%Dai
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

                  ! print*, "Tabij: ", CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j)
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
      integer :: p, q, r, s
      integer :: n_sp

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      n_sp = CCSD_instance%cont

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1_inter: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops
      ! T^{a}_{i}D^{a}_{i} = ...
      do a=nop+1, noc
        do i=1, nop

          do e=nop+1, noc
            do m=1, nop

              CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                - Allspecies(speciesId)%Tssame(e-nop,m)*CCSDinter(speciesId)%Fkc_aba(m,e-nop) !eq 7
            end do
          end do          

          do ee=nops+1, nocs
            do mm=1, nops

              CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                + (Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                    CCSDinter(OtherspeciesId)%Fkca_ab(mm,ee-nops)) !eq 8

              ! if (speciesId>OtherspeciesId) print*, "Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,mm)
              do m=1, nop

                if (speciesId<OtherspeciesId) then
                  p=m
                  q=mm
                  r=i
                  s=ee
                else
                  p=mm
                  q=m
                  r=ee
                  s=i
                end if

                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  - Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,m,mm)* &
                      (0.25*spintm(n_sp)%valuesp(p,q,r,s)) !eq 9

              end do

              do e=nop+1, noc

                if (speciesId<OtherspeciesId) then
                  p=a
                  q=mm
                  r=e
                  s=ee
                else
                  p=mm
                  q=a
                  r=ee
                  s=e
                end if

                CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i) &
                  + Allinterspecies(speciesId)%Tdsame(e-nop,ee-nops,i,mm)* &
                      (0.25*spintm(n_sp)%valuesp(p,q,r,s)) !eq 10
              end do

            end do
          end do

          CCSDT1T2(speciesId)%Tai(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)/CCSDinit(speciesId)%Dai(a,i)
          Allspecies(speciesId)%Tssame(a-nop,i) = CCSDT1T2(speciesId)%Tai(a-nop,i)
          print*, "Tai: ", CCSDT1T2(speciesId)%Tai(a-nop,i)
          print*, "Dai: ", CCSDinit(speciesId)%Dai(a,i)
          ! write(*,*) a,i,Allspecies(speciesId)%Tssame(a,i),CCSDT1T2(speciesId)%Tai(a,i)
        end do
      end do
      print*, "nops 12"
      ! print*, "Tai: ", CCSDT1T2(speciesId)%Tai
      ! print*, "Dai: ", CCSDinit(speciesId)%Dai
      
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

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T2_inter: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

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
                          CCSDinter(speciesId)%Wkkcc(a-nop,mm,i,ee-nops)) ) !eq 9

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

              ! print*, "Tabij: ", CCSDT1T2(speciesId)%Tabij(a-nop,b-nop,i,j)
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
      integer :: p, q, r, s

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species
      ! if (speciesId>OtherspeciesId) print*, "inside CCSD_T2AB Allinterspecies(speciesId)%Tdsame: ", &
        ! Allinterspecies(speciesId)%Tdsame

      ! number of transformed integrals matrix for speciesId and OtherspeciesId
      ! n_sp = Tix2(speciesId, OtherspeciesId, num_species)
      n_sp = CCSD_instance%cont
      ! T^{aB}_{iJ}D^{aB}_{iJ} = ...
      do a=nop+1, noc
        do i=1, nop
          do bb=nops+1, nocs
            do jj=1, nops

              do e=nop+1, noc
                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 + 0.5*(Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,i,jj)* &
                    CCSDinter(speciesId)%Faca(a-nop,e-nop)) !eq 5

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  + 0.5*(Allspecies(speciesId)%Tssame(e-nop,i)* &
                      CCSDinter(speciesId)%Wcia(a-nop,bb-nops,e-nop,jj)) !eq 2

                do mm=1, nops

                  if (speciesId<OtherspeciesId) then
                    p=a
                    q=bb
                    r=e
                    s=jj
                  else
                    p=bb
                    q=a
                    r=jj
                    s=e
                  end if

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    - 0.5*(Allinterspecies(speciesId)%tau(e-nop,bb-nops,i,mm)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 23
                end do

                do ee=nops+1, nocs

                  if (speciesId<OtherspeciesId) then
                    p=a
                    q=bb
                    r=e
                    s=ee
                  else
                    p=bb
                    q=a
                    r=ee
                    s=e
                  end if

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    + 0.5*(Allinterspecies(speciesId)%tau(e-nop,ee-nops,i,jj)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 25
                end do
              end do
              !print*, "nops 29"

              do ee=nops+1, nocs

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 + 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,ee-nops,i,jj)* &
                    CCSDinter(OtherspeciesId)%Fbcb(bb-nops,ee-nops)) !eq 6

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  + 0.5*(Allspecies(OtherspeciesId)%Tssame(ee-nops,jj)* &
                      CCSDinter(speciesId)%Wcjb(a-nop,bb-nops,i,ee-nops)) !eq 4

                do m=1, nop

                  if (speciesId<OtherspeciesId) then
                    p=m
                    q=bb
                    r=i
                    s=ee
                  else
                    p=bb
                    q=m
                    r=ee
                    s=i
                  end if

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    - 0.5*(Allinterspecies(speciesId)%tau(a-nop,ee-nops,m,jj)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 24
                end do
              end do
              !print*, "nops 30"

              do m=1, nop

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 - 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,m,jj)* &
                    CCSDinter(speciesId)%Fkia(m,i)) !eq 7

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  - 0.5*(Allspecies(speciesId)%Tssame(a-nop,m)* &
                      CCSDinter(speciesId)%Waka(m,bb-nops,i,jj)) !eq 1

                do e=nop+1, noc

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                   + 0.5*(Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,m,jj)* &
                      CCSDinter(speciesId)%Wakic_a(a-nop,m,i,e-nop)) !eq 9

                  if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                  
                    CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                      + 0.5*(Allspecies(speciesId)%Tdsame(a-nop,e-nop,i,m)* &
                          CCSDinter(OtherspeciesId)%Wbkjc(m,bb-nops,e-nop,jj)) !eq 13

                    do f=nop+1, noc
                      do n=1, nop

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%Tdsame(a-nop,f-nop,i,n)* &
                            CCSDinter(speciesId)%Wklcd_a(m,n,e-nop,f-nop)) !eq 11                

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%tau(a-nop,e-nop,i,m)* &
                            Allinterspecies(speciesId)%Tdsame(f-nop,bb-nops,n,jj)* &
                              spints(speciesId)%valuesp(m,n,e,f)) !eq 16
                      end do
                    end do
                  end if

                  do ee=nops+1, nocs
                    !\check{\tau} alpha
                    if (nop>=2) then ! kind of interaction just for two or more particles of the principal species

                      if (speciesId<OtherspeciesId) then
                        p=m
                        q=bb
                        r=e
                        s=ee
                      else
                        p=bb
                        q=m
                        r=ee
                        s=e
                      end if

                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_a(a-nop,e-nop,ee-nops,m,i,jj)* &
                            spintm(n_sp)%valuesp(p,q,r,s)) !eq 20
                    end if

                    do mm=1, nops
                      if ((nop>=2) .and. (nops>=2)) then ! kind of interaction for two or more particles for each species

                        if (speciesId<OtherspeciesId) then
                          p=m
                          q=mm
                          r=e
                          s=ee
                        else
                          p=mm
                          q=m
                          r=ee
                          s=e
                        end if

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(speciesId)%ttau(a-nop,e-nop,i,m)* &
                            Allspecies(OtherspeciesId)%ttau(bb-nops,ee-nops,jj,mm)* &
                              spintm(n_sp)%valuesp(p,q,r,s)) !eq 15
                      end if

                    end do
                  end do

                  if (nop>=2) then ! kind of interaction just for two or more particles of the principal species
                    do mm=1, nops
                      !\check{\tau} alpha

                      if (speciesId<OtherspeciesId) then
                        p=m
                        q=mm
                        r=e
                        s=jj
                      else
                        p=mm
                        q=m
                        r=jj
                        s=e
                      end if

                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_a(a-nop,e-nop,bb-nops,i,m,mm)* &
                            spintm(n_sp)%valuesp(p,q,r,s)) !eq 21
                    end do
                  end if

                end do
              end do
              !print*, "nops 31"

              if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                do ff=nops+1, nocs
                  do nn=1, nops

                    CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                     + 0.5*(Allspecies(OtherspeciesId)%Tdsame(bb-nops,ff-nops,jj,nn)* &
                        CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops)) !eq 12
                  end do
                end do
              end if


              do mm=1, nops

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                 - 0.5*(Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,mm)* & !OtherspeciesId to speciesId
                    CCSDinter(OtherspeciesId)%Fkjb(mm,jj)) !eq 8

                CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                  - 0.5*(Allspecies(OtherspeciesId)%Tssame(bb-nops,mm)* &
                      CCSDinter(speciesId)%Wbkb(a-nop,mm,i,jj)) !eq 3

                do ee=nops+1, nocs

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                   + 0.5*(Allinterspecies(OtherspeciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                      CCSDinter(OtherspeciesId)%Wbkjc_b(bb-nops,mm,jj,ee-nops)) !eq 10

                  if (nops>=2) then ! kind of interaction just for two or more particles of the another species

                    CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                     + 0.5*(Allspecies(OtherspeciesId)%Tdsame(bb-nops,ee-nops,jj,mm)* &
                        CCSDinter(speciesId)%Wakic(a-nop,mm,i,ee-nops)) !eq 14

                    do ff=nops+1, nocs
                      do nn=1, nops

                        ! CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        !  + 0.5*(Allspecies(OtherspeciesId)%Tdsame(bb-nops,ff-nops,jj,nn)* &
                        !     CCSDinter(OtherspeciesId)%Wklcd_b(mm,nn,ee-nops,ff-nops))                      !eq 12

                        CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                         + 0.5*(Allspecies(OtherspeciesId)%tau(bb-nops,ff-nops,jj,nn)* &
                            Allinterspecies(OtherspeciesId)%Tdsame(a-nop,ee-nops,i,mm)* &
                              spints(OtherspeciesId)%valuesp(mm,nn,ee,ff)) !eq 17

                      end do
                    end do
                  end if

                  if (nops>=2) then ! kind of interaction just for two or more particles of the another species
                    do e=nop+1, noc
                    !\check{\tau} beta

                      if (speciesId<OtherspeciesId) then
                        p=a
                        q=mm
                        r=e
                        s=ee
                      else
                        p=mm
                        q=a
                        r=ee
                        s=e
                      end if

                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_b(bb-nops,ee-nops,e-nop,mm,jj,i)* &
                            spintm(n_sp)%valuesp(p,q,r,s)) !eq 19
                    end do

                    do m=1, nop
                    !\check{\tau} beta

                      if (speciesId<OtherspeciesId) then
                        p=m
                        q=mm
                        r=i
                        s=ee
                      else
                        p=mm
                        q=m
                        r=ee
                        s=i
                      end if

                      CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                        - 0.5*(Allinterspecies(speciesId)%chtau_b(bb-nops,ee-nops,a-nop,jj,mm,m)* &
                            spintm(n_sp)%valuesp(p,q,r,s)) !eq 18
                    end do
                  end if

                end do
                
                do m=1, nop

                  if (speciesId<OtherspeciesId) then
                    p=m
                    q=mm
                    r=i
                    s=jj
                  else
                    p=mm
                    q=m
                    r=jj
                    s=i
                  end if

                  CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                    + 0.5*(Allinterspecies(speciesId)%tau(a-nop,bb-nops,m,mm)* &
                        spintm(n_sp)%valuesp(p,q,r,s)) !eq 22
                end do
              end do
              !print*,"nops 32"
              
              ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
              CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj) &
                /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                  -Allspecies(speciesId)%HF_fs%values(a,a) -Allspecies(OtherspeciesId)%HF_fs%values(bb,bb))
                  
              ! if (speciesId>OtherspeciesId) print*, "inside CCSD_T2AB ...%Tabij_AB", CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj)
              ! if (speciesId>OtherspeciesId) print*, "inside CCSD_T2AB ...%Tdsame", Allinterspecies(speciesId)%Tdsame(e-nop,bb-nops,i,jj)
              Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) = CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj)
              ! print*, "Allinterspecies(speciesId)%Tdsame: ", Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj)

              ! \ddot{\tau} 
              Allinterspecies(speciesId)%tau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                - Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)
              
              ! \tilde{\tau} 
              Allinterspecies(speciesId)%intau(a-nop,bb-nops,i,jj) = Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                + 0.5*Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)
              
              ! print*, "Tabij_AB: ", CCSDT1T2(speciesId)%Tabij_AB(a-nop,bb-nops,i,jj)

            end do
          end do

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
              do bb=nops+1, nocs
                do ii=1, nops
                  do jj=1, nops
                    !  \check{\tau} triple excitation holy shit!!
                    Allinterspecies(speciesId)%chtau_b(aa-nops,bb-nops,a-nop,ii,jj,i) = Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)* &
                      Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allspecies(speciesId)%Tssame(a-nop,i) &
                      + 0.5*( Allspecies(OtherspeciesId)%Tssame(aa-nops,ii)*Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,jj) &
                        + Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)*Allinterspecies(speciesId)%Tdsame(a-nop,aa-nops,i,ii) &
                          + Allspecies(speciesId)%Tssame(a-nop,i)*Allspecies(OtherspeciesId)%Tdsame(aa-nops,bb-nops,ii,jj) )
                  end do
                end do
              end do
            end do
          end if

              ! write(*,*) a,b,i,j,Allspecies(speciesId)%Tdsame(a,b,i,j),CCSDT1T2(speciesId)%Tabij_AB(a,b,i,j)
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
      integer :: stoped=0
      integer :: aux_cont=0
      integer :: e_cont=0
      integer :: n_sp=0
      integer :: max, min
      real(8) :: intra=0
      real(8) :: convergence = 1.0_8
      real(8) :: convergence_int = 1.0_8
      real(8), allocatable :: e_same_ccd(:)
      real(8), allocatable :: a_same_ccd(:)
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
      if (allocated(CCSD_instance%convergence_same_amp)) deallocate(CCSD_instance%convergence_same_amp)
      allocate(CCSD_instance%convergence_same_amp(num_species))

      ! allocate array for results of intra-species in CCSD
      if (allocated(CCSD_instance%convergence_diff)) deallocate(CCSD_instance%convergence_diff)
      allocate(CCSD_instance%convergence_diff(num_inter))

      ! allocate array for many results by species in CCSD
      if (allocated(e_same_ccd)) deallocate(e_same_ccd)
      allocate(e_same_ccd(num_species))
      e_same_ccd=0.0_8

      ! allocate array for many results by species in CCSD
      if (allocated(a_same_ccd)) deallocate(a_same_ccd)
      allocate(a_same_ccd(num_species))
      a_same_ccd=0.0_8


      ! allocate array for many results by species in CCSD
      if (allocated(e_diff_ccd)) deallocate(e_diff_ccd)
      allocate(e_diff_ccd(num_inter))
      e_diff_ccd=0.0_8

      print*,"times_i: ppp: ", times_i

      do i=counterID, num_species
        print*, "num_species: CCSD_constructor(): ", num_species, i
        call CCSD_constructor(i)
      end do

      do i=counterID, num_species!finalID
        print*, "num_species: CCSD_run(): ", num_species, i, finalID
        call CCSD_init(i)
        print*, "CCSD_init(i): ", i
      end do

      ! do i=counterID, finalID!num_species
        
        !If there are interspecies? REVISAR 01/10/2017
        if (times_i>0) then
          
          !search the principal species (speciesId) in all the options
          do j=1, times_i
            print*, "times_i", times_i
            print*, "counterID: ", counterID, "finalID: ", finalID
            !public to private            
            i_counterID(j) = CoupledCluster_instance%i_counterID(j)
            n_intersp(j) = CoupledCluster_instance%n_intersp(j)
            !Find the appropriate species in all the options
            ! if (i_counterID(j)==counterID) then
              !public variables used in the loop
              CCSD_instance%max=n_intersp(j)
              CCSD_instance%min=i_counterID(j)
              !for the order in call in spintm matrix
              aux_cont = n_sp + 1
              CCSD_instance%aux_cont = aux_cont
              !make all combinations with speciesId
              ! do jj=i_counterID(j)+1, n_intersp(j)
                n_sp = n_sp + 1
                CCSD_instance%cont = n_sp

                call CCSD_constructor_inter(i_counterID(j), n_intersp(j))!, num_inter)
                call CCSD_constructor_inter(n_intersp(j), i_counterID(j))!, num_inter)
                call CCSD_init_inter(i_counterID(j), n_intersp(j))!, num_inter)
                call CCSD_init_inter(n_intersp(j), i_counterID(j))!, num_inter)

                ! print*, "after CCSD_init_inter: ", Allinterspecies(jj)%Tdsame
                num_inter = num_inter + 1
              ! end do
            ! else
            !   !locating 
            !   n_sp = n_sp + (n_intersp(j) - i_counterID(j))
            !   CCSD_instance%cont = n_sp
            ! end if
          end do

          ! print*, "CCSD_loop(i, i+1): ", i
          ! call CCSD_loop(i, i+1) !fix i+1

        ! else

        !   call CCSD_loop(i)
        !   print*, "CCSD_loop(i): ", i
        
        end if
      ! end do

      ! do i=counterID, num_species

      !   do while (convergence >= 1.0D-8)
        
      !     call CCSD_same_species(i,e_same_ccd(i))
      !     e_same_ccd(i) = CoupledCluster_instance%CCSD_E_intra(i)

      !     convergence = CCSD_instance%convergence_same(i)
      !   end do
      !   convergence = 1.0_8
      ! end do


      ! if (times_i>0) then
      do while (convergence >= 1.0D-8)

        ! if ((CCSD_instance%min+1)>counterID) print*, "before times_i>0: ", Allinterspecies(CCSD_instance%min+1)%Tdsame
        stoped = stoped +1
        ! if (stoped>1) then
        !   print*, "CCSDT1T2(speciesId)%Tai(a-nop,i): ", Allspecies(1)%Tssame
        !   print*, "CCSDT1T2(OtherspeciesId)%Tai(a-nop,i): ", Allspecies(2)%Tssame
        ! end if

        do i=counterID, num_species
          call CCSD_same_species(i,e_same_ccd(i),a_same_ccd(i))
          e_same_ccd(i) = CoupledCluster_instance%CCSD_E_intra(i)
          a_same_ccd(i) = CoupledCluster_instance%CCSD_a_intra(i)
        end do
        if (times_i>0) then
          print*, "test PPPPP"
          do i=1, times_i!counterID, finalID
            CCSD_instance%num_i = 1
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            max = CCSD_instance%max
            min = CCSD_instance%min
            ! if ((min+1)>i) print*, "before do: ", Allinterspecies(min+1)%Tdsame
            ! do jj=min+1, max
              num_i = CCSD_instance%num_i
              print*, "num_i: ", num_i 
              print*, "e_cont: ", CCSD_instance%e_cont
              call CCSD_diff_species(i_counterID(i),n_intersp(i),e_diff_ccd(num_i))
              e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
              print*, "CCSDT1T2(speciesId)%Tai(a-nop,i): ", CCSDT1T2(i_counterID(i))%Tai
              ! convergence_int = CCSD_instance%convergence_diff(num_i)
            ! end do
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            ! do jj=min+1, max
              ! if (jj>i) print*, "before CCSD_diff_species: ", Allinterspecies(jj)%Tdsame
              num_i = num_i+1
              print*, "num_i: ", num_i
              print*, "e_cont: ", CCSD_instance%e_cont
              call CCSD_diff_species(n_intersp(i),i_counterID(i),e_diff_ccd(num_i))
              e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)

              ! if (jj>i) print*, "after CCSD_diff_species: ", Allinterspecies(jj)%Tdsame
              ! convergence_int = CCSD_instance%convergence_diff(num_i)
            ! end do
          end do

          do i=1, times_i
            CCSD_instance%num_i = 1
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            max = CCSD_instance%max
            min = CCSD_instance%min
            num_i = CCSD_instance%num_i
            call CCSD_diff_species_energy(i_counterID(i),n_intersp(i))
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            call CCSD_diff_species_energy(n_intersp(i),i_counterID(i))
          end do

          do i=1, times_i
            CCSD_instance%num_i = 1
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            max = CCSD_instance%max
            min = CCSD_instance%min
            num_i = CCSD_instance%num_i
            print*, "CCSD_T2_AB()"
            call CCSD_T2_AB(i_counterID(i),n_intersp(i), num_i)
            ! if (speciesId>OtherspeciesId) print*, "after CCSD_T2AB Allinterspecies(speciesId)%Tdsame: ", &
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            call CCSD_T2_AB(n_intersp(i),i_counterID(i), num_i)            
          end do
          
        end if

        print*, "CoupledCluster_instance%CCSD_E_intra: ", CoupledCluster_instance%CCSD_E_intra
        print*, "e_same_ccd: ", CCSD_instance%convergence_same
        print*, "CoupledCluster_instance%CCSD_A_intra: ", CoupledCluster_instance%CCSD_A_intra
        print*, "CoupledCluster_instance%CCSD_E_inter: ", CoupledCluster_instance%CCSD_E_inter
        print*, "e_diff_ccd: ", CCSD_instance%convergence_diff


        ! if (CCSD_instance%convergence_same(1) >= 1.0D-8) then
          ! convergence = CCSD_instance%convergence_same(1)
        ! else if (CCSD_instance%convergence_same(2) >= 1.0D-8) then
        !   convergence = CCSD_instance%convergence_same(2)
        ! else
         convergence = sum(CCSD_instance%convergence_same,dim=1) + & 
           sum(CCSD_instance%convergence_diff,dim=1)
        ! end if
        !convergence = CCSD_instance%convergence_same(1)!,dim=1)
        ! if (stoped>2) then
        !   convergence = CCSD_instance%convergence_same_amp(1)
        ! end if

        print*, "Temporary convergence: ", convergence
        print*, "Temporary total ccsd-apmo energy: ", sum(e_same_ccd, dim=1) + &
          sum(e_diff_ccd, dim=1)
        ! if (stoped>=3) stop "nops 1 - 5"

      end do
      ! end if

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
      print*, "Total correction CCSD energy: ", CoupledCluster_instance%CCSD_total_Energy
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
      print*, "INFORMATION IN CCSD_constructor() MP2_Energy for species: ", speciesId, ": ", CoupledCluster_instance%MP2_ECorr%values(speciesId)
      print*, "INFORMATION IN CCSD_constructor() CCSD_Energy of species " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_E_intra(speciesId)
      print*, "INFORMATION IN CCSD_constructor() CCSD_Energy of different species " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_E_inter(speciesId)
      ! print*, "INFORMATION IN CCSD_constructor() Td: ", sum(Allspecies(speciesId)%Tdsame,dim=1)
      ! print*, "INFORMATION IN CCSD_constructor() Dai: ", sum(CCSDinit(speciesId)%Dai,dim=1)
      ! print*, "INFORMATION IN CCSD_constructor() ttau: ", sum(Allspecies(speciesId)%ttau,dim=1)
      ! print*, "INFORMATION IN CCSD_constructor() tau: ",sum(Allspecies(speciesId)%tau,dim=1)
      print*, "INFORMATION IN CCSD_constructor() ccsdE: ", CoupledCluster_instance%CCSD_E_intra(speciesId)
      print*, "INFORMATION IN Tensor Tensor_index2: ", Tix2(2,1,3)
      

      !if (allocated(Allspecies(speciesId)%tau)) deallocate (Allspecies(speciesId)%tau)

      !call CCSD_destructor()
      print*, "CCSD_show()x2"

  end subroutine CCSD_show
end module CCSD_
