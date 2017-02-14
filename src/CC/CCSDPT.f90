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
  !         for inter-species
  ! @author CAOM
  subroutine CCSDPT_constructor_inter(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      ! integer, intent(in) :: num_inter

      integer noc, nop, nocs, nops
      
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSDPT_constructor_inter: noc=", noc, "nop=", nop, "nocs=", nocs, "nops=", nops
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! lowercase = alpha species | uppercase = beta species

      ! ! t^{ab}_{ij} amplitude for double excitation for same species
      if (allocated(Allinterspecies(speciesId)%Tt4same_abb)) deallocate(Allinterspecies(speciesId)%Tt4same_abb)
      allocate(Allinterspecies(speciesId)%Tt4same_abb(noc-nop,nocs-nops,nocs-nops,nop,nops,nops))
      Allinterspecies(speciesId)%Tt4same_abb(:,:,:,:,:,:) = 0.0_8

      if (allocated(Allinterspecies(speciesId)%Tt4same_aab)) deallocate(Allinterspecies(speciesId)%Tt4same_aab)
      allocate(Allinterspecies(speciesId)%Tt4same_aab(noc-nop,noc-nop,nocs-nops,nop,nop,nops))
      Allinterspecies(speciesId)%Tt4same_aab(:,:,:,:,:,:) = 0.0_8

      if (allocated(Allinterspecies(speciesId)%Tt5same_abb)) deallocate(Allinterspecies(speciesId)%Tt5same_abb)
      allocate(Allinterspecies(speciesId)%Tt5same_abb(noc-nop,nocs-nops,nocs-nops,nop,nops,nops))
      Allinterspecies(speciesId)%Tt5same_abb(:,:,:,:,:,:) = 0.0_8

      if (allocated(Allinterspecies(speciesId)%Tt5same_aab)) deallocate(Allinterspecies(speciesId)%Tt5same_aab)
      allocate(Allinterspecies(speciesId)%Tt5same_aab(noc-nop,noc-nop,nocs-nops,nop,nop,nops))
      Allinterspecies(speciesId)%Tt5same_aab(:,:,:,:,:,:) = 0.0_8

      print*, "fin"

  end subroutine CCSDPT_constructor_inter

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
                          + (spints(speciesId)%valuesp(d,i,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,j,k)) &
                            + (spints(speciesId)%valuesp(d,i,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,j,k)) &
                              + (spints(speciesId)%valuesp(d,j,b,c)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,k,i)) &
                                + (spints(speciesId)%valuesp(d,j,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,k,i)) &
                                  + (spints(speciesId)%valuesp(d,j,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,k,i)) &
                                    + (spints(speciesId)%valuesp(d,k,b,c)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,i,j)) &
                                      + (spints(speciesId)%valuesp(d,k,c,a)*Allspecies(speciesId)%Tdsame(b-nop,d-nop,i,j)) &
                                        + (spints(speciesId)%valuesp(d,k,a,b)*Allspecies(speciesId)%Tdsame(c-nop,d-nop,i,j))

                    end do
                    do l=1, nop

                      Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k) &
                        - (spints(speciesId)%valuesp(l,a,j,k)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,i,l)) &
                          - (spints(speciesId)%valuesp(l,b,j,k)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,i,l)) &
                            - (spints(speciesId)%valuesp(l,c,j,k)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,l)) &
                              - (spints(speciesId)%valuesp(l,a,k,i)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,j,l)) &
                                - (spints(speciesId)%valuesp(l,b,k,i)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,j,l)) &
                                  - (spints(speciesId)%valuesp(l,c,k,i)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,j,l)) &
                                    - (spints(speciesId)%valuesp(l,a,i,j)*Allspecies(speciesId)%Tdsame(b-nop,c-nop,k,l)) &
                                      - (spints(speciesId)%valuesp(l,b,i,j)*Allspecies(speciesId)%Tdsame(c-nop,a-nop,k,l)) &
                                        - (spints(speciesId)%valuesp(l,c,i,j)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,k,l))

                    end do

                  Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k) = Allspecies(speciesId)%Tt5same(a-nop,b-nop,c-nop,i,j,k) &
                    + (spints(speciesId)%valuesp(j,k,b,c)*Allspecies(speciesId)%Tssame(a-nop,i)) &
                      + (spints(speciesId)%valuesp(j,k,c,a)*Allspecies(speciesId)%Tssame(b-nop,i)) &
                        + (spints(speciesId)%valuesp(j,k,a,b)*Allspecies(speciesId)%Tssame(c-nop,i)) &
                          + (spints(speciesId)%valuesp(k,i,b,c)*Allspecies(speciesId)%Tssame(a-nop,j)) &
                            + (spints(speciesId)%valuesp(k,i,c,a)*Allspecies(speciesId)%Tssame(b-nop,j)) &
                              + (spints(speciesId)%valuesp(k,i,a,b)*Allspecies(speciesId)%Tssame(c-nop,j)) &
                                + (spints(speciesId)%valuesp(i,j,b,c)*Allspecies(speciesId)%Tssame(a-nop,k)) &
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
      
  end subroutine CCSDPT_init

  !>
  ! @brief Build a amplitudes and Denominators guesses from MP2 information for inter-species
  ! @author CAOM
  subroutine CCSDPT_init_inter(speciesId, OtherspeciesId)!, num_inter)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId
      ! integer, intent(in) :: num_inter

      integer noc, nocs, nop, nops, n_sp, num_species

      integer :: a, b, i, ii, j, d, l
      integer :: aa, bb, jj, cc, kk, dd, ll
      integer :: p, q, r, s, pp, qq, rr, ss
      integer :: p1, q1, r1, s1, p2, q2, r2, s2

      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop

      num_species = CoupledCluster_instance%num_species

      ! n_sp = Tix2(speciesId, OtherspeciesId, num_species)
      ! number of transformed integrals matrix for speciesId and OtherspeciesId
      n_sp =CCSD_instance%cont

      ! lowercase = alpha species | uppercase = beta species

      if (nops>=2) then
        do a=nop+1, noc
          do bb=nops+1, nocs
            do cc=nops+1, nocs
              do i=1, nop
                do jj=1, nops
                  do kk=1, nops

                    ! if(OtherspeciesId<speciesId) print*, "CCSDTTTTTTTT"


                    do dd=nops+1, nocs

                      if (speciesId<OtherspeciesId) then
                        p=i
                        q=dd
                        r=a
                        s=cc
                        pp=i
                        qq=dd
                        rr=a
                        ss=bb
                      else
                        p=dd
                        q=i
                        r=cc
                        s=a
                        pp=dd
                        qq=i
                        rr=bb
                        ss=a
                      end if

                      Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) = &
                        Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) &
                          + (spints(OtherspeciesId)%valuesp(dd,kk,bb,cc)*Allinterspecies(speciesId)%Tdsame(a-nop,dd-nops,i,jj)) &
                            + (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(OtherspeciesId)%Tdsame(bb-nops,dd-nops,jj,kk)) &
                              + (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(OtherspeciesId)%Tdsame(cc-nops,dd-nops,kk,jj))

                    end do
                    do ll=1, nops


                      if (speciesId<OtherspeciesId) then
                        p=a
                        q=ll
                        r=i
                        s=kk
                        pp=a
                        qq=ll
                        rr=i
                        ss=jj
                      else
                        p=ll
                        q=a
                        r=kk
                        s=i
                        pp=ll
                        qq=a
                        rr=jj
                        ss=i
                      end if

                      Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) = &
                        Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) &
                          - (spints(OtherspeciesId)%valuesp(ll,cc,jj,kk)*Allinterspecies(speciesId)%Tdsame(a-nop,bb-nops,i,ll)) &
                            - (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(OtherspeciesId)%Tdsame(bb-nops,cc-nops,jj,ll)) &
                              - (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(OtherspeciesId)%Tdsame(bb-nops,cc-nops,ll,kk))

                    end do

                  if (speciesId<OtherspeciesId) then
                    p=i
                    q=kk
                    r=a
                    s=cc
                    pp=i
                    qq=jj
                    rr=a
                    ss=cc
                    p1=i
                    q1=kk
                    r1=a
                    s1=bb
                    p2=i
                    q2=jj
                    r2=a
                    s2=bb                    
                  else
                    p=kk
                    q=i
                    r=cc
                    s=a
                    pp=jj
                    qq=i
                    rr=cc
                    ss=a
                    p1=kk
                    q1=i
                    r1=bb
                    s1=a
                    p2=jj
                    q2=i
                    r2=bb
                    s2=a                    
                  end if

                  Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) = &
                    Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) &
                      + (spints(OtherspeciesId)%valuesp(jj,kk,bb,cc)*Allspecies(speciesId)%Tssame(a-nop,i)) &
                        + (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(OtherspeciesId)%Tssame(bb-nops,jj)) &
                          - (spintm(n_sp)%valuesp(p1,q1,r1,s1)*Allspecies(OtherspeciesId)%Tssame(cc-nops,jj)) &
                            - (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(OtherspeciesId)%Tssame(bb-nops,kk)) &
                              + (spintm(n_sp)%valuesp(p2,q2,r2,s2)*Allspecies(OtherspeciesId)%Tssame(cc-nops,kk))


                  !***
                  Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) = &
                    Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)/ &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))

                  Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk) = &
                    Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)/ &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))

                  
                  if ((Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))==0) then
                    Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)=0.0_8
                    Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)=0.0_8
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


      if (nop>=2) then
        do a=nop+1, noc
          do b=nop+1, noc
            do cc=nops+1, nocs
              do i=1, nop
                do j=1, nop
                  do kk=1, nops


                    do d=nops+1, nocs

                      if (speciesId<OtherspeciesId) then
                        p=d
                        q=kk
                        r=b
                        s=cc
                        pp=d
                        qq=kk
                        rr=a
                        ss=cc
                      else
                        p=kk
                        q=d
                        r=cc
                        s=b
                        pp=kk
                        qq=d
                        rr=cc
                        ss=a
                      end if

                      Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk) = &
                        Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk) &
                          + (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(speciesId)%Tdsame(a-nop,d-nop,i,j)) &
                            + (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(speciesId)%Tdsame(d-nop,b-nop,i,j)) &
                              + (spints(speciesId)%valuesp(d,j,a,b)*Allinterspecies(speciesId)%Tdsame(d-nop,cc-nops,i,kk))

                    end do
                    do l=1, nops

                      if (speciesId<OtherspeciesId) then
                        p=l
                        q=cc
                        r=j
                        s=kk
                        pp=l
                        qq=cc
                        rr=i
                        ss=kk
                      else
                        p=cc
                        q=l
                        r=kk
                        s=j
                        pp=cc
                        qq=l
                        rr=kk
                        ss=i
                      end if

                      Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk) = &
                        Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk) &
                          - (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,i,l)) &
                            - (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(speciesId)%Tdsame(a-nop,b-nop,l,j)) &
                              - (spints(speciesId)%valuesp(l,b,i,j)*Allinterspecies(speciesId)%Tdsame(a-nop,cc-nops,l,kk))

                    end do

                  if (speciesId<OtherspeciesId) then
                    p=j
                    q=kk
                    r=b
                    s=cc
                    pp=j
                    qq=kk
                    rr=a
                    ss=cc
                    p1=i
                    q1=kk
                    r1=b
                    s1=cc
                    p2=i
                    q2=kk
                    r2=a
                    s2=cc                    
                  else
                    p=kk
                    q=j
                    r=cc
                    s=b
                    pp=kk
                    qq=j
                    rr=cc
                    ss=a
                    p1=kk
                    q1=i
                    r1=cc
                    s1=b
                    p2=kk
                    q2=i
                    r2=cc
                    s2=a                    
                  end if

                  Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk) = &
                    Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk) &
                      + (spintm(n_sp)%valuesp(p,q,r,s)*Allspecies(speciesId)%Tssame(a-nop,i)) &
                        - (spintm(n_sp)%valuesp(pp,qq,rr,ss)*Allspecies(speciesId)%Tssame(b-nop,i)) &
                          - (spintm(n_sp)%valuesp(p1,q1,r1,s1)*Allspecies(speciesId)%Tssame(a-nop,j)) &
                            + (spintm(n_sp)%valuesp(p2,q2,r2,s2)*Allspecies(speciesId)%Tssame(b-nop,j)) &
                              + (spints(speciesId)%valuesp(i,j,a,b)*Allspecies(OtherspeciesId)%Tssame(cc-nops,kk))


                  !***
                  Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk) = &
                    Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk)/ &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))

                  Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk) = &
                    Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk)/ &
                      (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))

                  if ((Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                        + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                          - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))==0) then
                    Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk)=0.0_8
                    Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk)=0.0_8
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

      print*, "before before inter" 

  end subroutine CCSDPT_init_inter

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

                  tmp_ccsd_ts5E = tmp_ccsd_ts5E + ((0.0277777777)*Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k)* &
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
  subroutine CCSDPT_diff_species(speciesId, OtherspeciesId)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: OtherspeciesId

      integer :: noc, nocs, nop, nops
      integer :: a, aa, b, bb, c, cc, d, dd
      integer :: i, ii, j, jj, k, kk, l, ll
      real(8) :: tmp_ccsd_t4E_inter_abb
      real(8) :: tmp_ccsd_t4E_inter_aab
      real(8) :: tmp_ccsd_ts5E_inter_abb
      real(8) :: tmp_ccsd_ts5E_inter_aab



      ! !Initialization of private variables from public variables
      noc = Allspecies(speciesId)%noc
      nocs = Allspecies(OtherspeciesId)%noc
      nop = Allspecies(speciesId)%nop
      nops = Allspecies(OtherspeciesId)%nop
      ! num_species = CoupledCluster_instance%num_species
      ! ! write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nop=", nop
      ! print*, "T1T2_constructor", convergence, noc, nop, speciesId

      ! prev_ccsdE_int = e_ccsd
      ! prev_ampl = v_ampl_int

      ! print*, "speciesId: ", speciesId

      ! num_inter = CCSD_instance%num_i
      ! e_cont = CCSD_instance%e_cont

       
      ! ! Resolve CCSD equation of energy

      tmp_ccsd_t4E_inter_abb=0.0_8
      tmp_ccsd_ts5E_inter_abb=0.0_8
      do a=nop+1, noc
        do bb=nops+1, nocs
          do cc=nops+1, nocs
            do i=1, nop
              do jj=1, nops
                do kk=1, nops
                  tmp_ccsd_t4E_inter_abb = tmp_ccsd_t4E_inter_abb + ((0.08333333333)*Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)* &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                      + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))* &
                            Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk))
                  ! print*, "tmp_ccsd_t4E_inter_abb: ", tmp_ccsd_t4E_inter_abb, Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k), &
                  !   (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                  !     + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                  !       - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))

                  tmp_ccsd_ts5E_inter_abb = tmp_ccsd_ts5E_inter_abb + ((0.08333333333)*Allinterspecies(speciesId)%Tt4same_abb(a-nop,bb-nops,cc-nops,i,jj,kk)* &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(OtherspeciesId)%HF_fs%values(jj,jj) &
                      + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(OtherspeciesId)%HF_fs%values(bb,bb) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))* &
                          Allinterspecies(speciesId)%Tt5same_abb(a-nop,bb-nops,cc-nops,i,jj,kk))

                end do
              end do
            end do
          end do
        end do
      end do

      print *, "ccsd_t4E inter_aab", tmp_ccsd_t4E_inter_aab
      print *, "ccsd_ts5E inter_aab", tmp_ccsd_ts5E_inter_aab

      tmp_ccsd_t4E_inter_aab=0.0_8
      tmp_ccsd_ts5E_inter_aab=0.0_8
      do a=nop+1, noc
        do b=nop+1, noc
          do cc=nops+1, nocs
            do i=1, nop
              do j=1, nop
                do kk=1, nops
                  tmp_ccsd_t4E_inter_aab = tmp_ccsd_t4E_inter_aab + ((0.08333333333)*Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk)* &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                      + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))* &
                            Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk))
                  ! print*, "tmp_ccsd_t4E_inter_abb: ", tmp_ccsd_t4E_inter_abb, Allspecies(speciesId)%Tt4same(a-nop,b-nop,c-nop,i,j,k), &
                  !   (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                  !     + Allspecies(speciesId)%HF_fs%values(k,k) - Allspecies(speciesId)%HF_fs%values(a,a) &
                  !       - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(speciesId)%HF_fs%values(c,c))

                  tmp_ccsd_ts5E_inter_aab = tmp_ccsd_ts5E_inter_aab + ((0.08333333333)*Allinterspecies(speciesId)%Tt4same_aab(a-nop,b-nop,cc-nops,i,j,kk)* &
                    (Allspecies(speciesId)%HF_fs%values(i,i) + Allspecies(speciesId)%HF_fs%values(j,j) &
                      + Allspecies(OtherspeciesId)%HF_fs%values(kk,kk) - Allspecies(speciesId)%HF_fs%values(a,a) &
                        - Allspecies(speciesId)%HF_fs%values(b,b) - Allspecies(OtherspeciesId)%HF_fs%values(cc,cc))* &
                          Allinterspecies(speciesId)%Tt5same_aab(a-nop,b-nop,cc-nops,i,j,kk))

                end do
              end do
            end do
          end do
        end do
      end do

      print *, "ccsd_t4E inter", tmp_ccsd_t4E_inter_aab
      print *, "ccsd_ts5E inter", tmp_ccsd_ts5E_inter_aab

      CoupledCluster_instance%CCSD_T4_E_inter_aab(speciesId) = (0.5)*tmp_ccsd_t4E_inter_aab
      CoupledCluster_instance%CCSD_TS5_E_inter_aab(speciesId) = (0.5)*tmp_ccsd_ts5E_inter_aab

      CoupledCluster_instance%CCSD_T4_E_inter_abb(speciesId) = (0.5)*tmp_ccsd_t4E_inter_abb
      CoupledCluster_instance%CCSD_TS5_E_inter_abb(speciesId) = (0.5)*tmp_ccsd_ts5E_inter_abb
      !0.5 is due to interaction alpha-beta and beta-alpha
  end subroutine CCSDPT_diff_species

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

      do i=counterID, num_species
        print*, "num_species: CCSD(T)_constructor(): ", num_species, i
        call CCSDPT_constructor(i)
      end do

      do i=counterID, num_species
        print*, "num_species: CCSD(T)_run(): ", num_species, i, finalID
        call CCSDPT_init(i)
        print*, "CCSD(T)_init(i): ", i
      end do

        
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
          n_sp = n_sp + 1
          CCSD_instance%cont = n_sp

          call CCSDPT_constructor_inter(i_counterID(j), n_intersp(j))!, num_inter)
          call CCSDPT_constructor_inter(n_intersp(j), i_counterID(j))!, num_inter)
          call CCSDPT_init_inter(i_counterID(j), n_intersp(j))!, num_inter)
          call CCSDPT_init_inter(n_intersp(j), i_counterID(j))!, num_inter)

          num_inter = num_inter + 1
        end do

      end if

      do i=counterID, num_species
        call CCSDPT_same_species(i)
          ! e_same_ccd(i) = CoupledCluster_instance%CCSD_E_intra(i)
          ! a_same_ccd(i) = CoupledCluster_instance%CCSD_a_intra(i)
      end do
        if (times_i>0) then
          do i=1, times_i
            CCSD_instance%num_i = 1
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont
            ! max = CCSD_instance%max
            ! min = CCSD_instance%min

            num_i = CCSD_instance%num_i
            print*, "num_i: ", num_i 
            print*, "e_cont: ", CCSD_instance%e_cont
            call CCSDPT_diff_species(i_counterID(i),n_intersp(i))
            ! e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
            ! a_diff_ccd(num_i) = CoupledCluster_instance%CCSD_a_inter(num_i)
           
            CCSD_instance%cont = CCSD_instance%aux_cont
            CCSD_instance%e_cont = CCSD_instance%aux_cont

            num_i = num_i+1
            print*, "num_i: ", num_i
            print*, "e_cont: ", CCSD_instance%e_cont
            call CCSDPT_diff_species(n_intersp(i),i_counterID(i))
            ! e_diff_ccd(num_i) = CoupledCluster_instance%CCSD_E_inter(num_i)
            ! a_diff_ccd(num_i) = CoupledCluster_instance%CCSD_a_inter(num_i)

          end do
          
        end if


      do i=counterID, num_species
        call CCSDPT_show(i)
        ! intra = intra + CoupledCluster_instance%CCSD_E_intra(i)
      end do

      CoupledCluster_instance%CCSDPT_once_Energy = sum(CoupledCluster_instance%CCSD_T4_E_intra,dim=1) + &
        sum(CoupledCluster_instance%CCSD_TS5_E_intra,dim=1)
      CoupledCluster_instance%CCSDPT_twice_Energy = CoupledCluster_instance%CCSDPT_once_Energy +&
        sum(CoupledCluster_instance%CCSD_T4_E_inter_aab,dim=1) + sum(CoupledCluster_instance%CCSD_TS5_E_inter_aab,dim=1) + &
          sum(CoupledCluster_instance%CCSD_T4_E_inter_abb,dim=1) + sum(CoupledCluster_instance%CCSD_TS5_E_inter_abb,dim=1)
      CoupledCluster_instance%CCSD_T_total_Energy = CoupledCluster_instance%CCSD_total_Energy + CoupledCluster_instance%CCSDPT_once_Energy + &
        CoupledCluster_instance%CCSDPT_twice_Energy
      ! CoupledCluster_instance%CCSD_once_Energy = sum(CoupledCluster_instance%CCSD_E_intra,dim=1)
      ! CoupledCluster_instance%CCSD_twice_Energy = sum(CoupledCluster_instance%CCSD_E_inter,dim=1)
      ! CoupledCluster_instance%CCSD_total_Energy = CoupledCluster_instance%CCSD_once_Energy + &
      !   CoupledCluster_instance%CCSD_twice_Energy

      ! ! CoupledCluster_instance%CCSD_ones_Energy = intra

      print*, "Total [T] correction same species energy: ", sum(CoupledCluster_instance%CCSD_T4_E_intra,dim=1) 
      print*, "Total (T) correction same species energy: ", CoupledCluster_instance%CCSDPT_once_Energy 
  	  print*, "Total [T] correction different species energy: ", sum(CoupledCluster_instance%CCSD_T4_E_inter_aab,dim=1) +&
  	  	sum(CoupledCluster_instance%CCSD_T4_E_inter_abb,dim=1) + sum(CoupledCluster_instance%CCSD_T4_E_intra,dim=1)
      print*, "Total (T) correction different species energy: ", CoupledCluster_instance%CCSDPT_twice_Energy
      print*, "Total CCSD(T) correction energy: ", CoupledCluster_instance%CCSD_T_total_Energy
      print*, "Total CCSD(T) energy: ", CoupledCluster_instance%CCSD_T_total_Energy &
        + CoupledCluster_instance%HF_energy
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
      print*, "E_TS[5] Energy of species " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "CCSD + E_T[4] = CCSD[T] correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "CCSD + E_TS[5] = CCSD(T) correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId) + CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "CCSD + E_T[4] = CCSD[T] correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "CCSD + E_TS[5] = CCSD(T) correction energy: ", CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId) + CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "Total CCSD[T] energy: ", CoupledCluster_instance%HF_energy + CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId)
      print*, "Total CCSD(T) energy: ", CoupledCluster_instance%HF_energy + CoupledCluster_instance%CCSD_E_intra(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_intra(speciesId) + CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)
      print*, "INFORMATION OF CCSD(T)-APMO METHOD"
      print*, "E_T[4]-APMO alpha-beta-beta Energy " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_T4_E_inter_abb(speciesId)
      print*, "E_TS[5]-APMO alpha-beta-beta Energy " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_TS5_E_inter_abb(speciesId)
      print*, "E_T[4]-APMO alpha-alpha-beta Energy" , speciesId, ": ", &
        CoupledCluster_instance%CCSD_T4_E_inter_aab(speciesId)
      print*, "E_TS[5]-APMO alpha-alpha-beta Energy " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_TS5_E_inter_aab(speciesId)
      print*, "E_T[4] APMO correction Energy" , speciesId, ": ", CoupledCluster_instance%CCSD_T4_E_inter_aab(speciesId) + &
        CoupledCluster_instance%CCSD_T4_E_inter_abb(speciesId)
      print*, "E_TS[5] APMO correction Energy " , speciesId, ": ", CoupledCluster_instance%CCSD_TS5_E_inter_aab(speciesId) + &
        CoupledCluster_instance%CCSD_TS5_E_inter_abb(speciesId)

      print*, "CCSD-APMO + E_T[4] + E_T[4]-APMO = CCSD[T]-APMO correction energy inter-species: " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_E_inter(speciesId) + (2*CoupledCluster_instance%CCSD_T4_E_intra(speciesId)) + &
          CoupledCluster_instance%CCSD_T4_E_inter_aab(speciesId) + CoupledCluster_instance%CCSD_T4_E_inter_abb(speciesId)
      print*, "CCSD[T]-APMO + E_TS[5] + E_TS[5]-APMO = CCSD(T)-APMO correction energy inter-species: " , speciesId, ": ", &
        CoupledCluster_instance%CCSD_E_inter(speciesId) + (2*CoupledCluster_instance%CCSD_T4_E_intra(speciesId)) + &
          (2*CoupledCluster_instance%CCSD_TS5_E_intra(speciesId)) + CoupledCluster_instance%CCSD_T4_E_inter_aab(speciesId) + &
            CoupledCluster_instance%CCSD_T4_E_inter_abb(speciesId) + CoupledCluster_instance%CCSD_TS5_E_inter_aab(speciesId) + &
              CoupledCluster_instance%CCSD_TS5_E_inter_abb(speciesId)




  end subroutine CCSDPT_show

end module CCSDPT_