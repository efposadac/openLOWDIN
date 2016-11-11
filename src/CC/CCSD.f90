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
      
      
      real(8), allocatable :: Tssame(:,:)
      real(8), allocatable :: Tdsame(:,:,:,:)
      real(8), allocatable :: tau(:,:,:,:)
      real(8), allocatable :: ttau(:,:,:,:)
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
  type(CCSDiter), public :: CCSDT1T2


contains

  !>
  ! @brief Constructor of the class
  ! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CCSD_constructor()
      implicit none

      integer noc, nocs, nop, nops
      integer :: speciesId=1
      
      ! call CoupledCluster_pairing_function(1,2)
      
      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      write(*, "(A,I4,A,I4,A,I4,A,I4) ") "CCSD_constructor: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops
      ! allocate all that you can...
      ! All transformed integrals are loaded using the previous subroutine

      ! Denominator in T1 D^{a}_{i}
      if (allocated(CCSDinit%Dai)) deallocate (CCSDinit%Dai)
      allocate(CCSDinit%Dai(noc,noc))
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
      ! if (allocated(CCSD_instance%Tssame)) deallocate (CCSD_instance%Tssame)
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

      integer noc, nocs, nop, nops
      integer :: a, b, i, j
      integer :: speciesId=1

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      ! Effective two-particle excitation operators \tilde{\tau} and \tau:

      ! \tilde{\tau}^{ab}_{ij} = t^{ab}_{ij} + \frac{1}{2}(t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j})
      ! \tau^{ab}_{ij} = t^{ab}_{ij} + t^{a}_{i}t^{b}_{j} - t^{b}_{i}t^{a}_{j}

      print*, "before loop"
      print*, CCSD_instance%Tdsame(1,1,1,1)
      print*, CCSD_instance%Tssame(1,1)
      print*, CCSD_instance%ttau(1,1,1,1)
      print*, CCSD_instance%tau(1,1,1,1)
      ! CCSD_instance%Tdsame(1,1,1,1) = CCSD_instance%tau(1,1,1,1)
      print*, "before loop"
      do a=nop+1, noc
        do b=nop+1, noc
          do i=1, nop
            do j=1, nop


              CCSD_instance%Tdsame(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                +( (spints(speciesId)%valuesp(i,j,a,b))/( Allspecies(speciesId)%HF_fs%values(i,i)+ &
                  Allspecies(speciesId)%HF_fs%values(j,j) -Allspecies(speciesId)%HF_fs%values(a,a)- &
                    Allspecies(speciesId)%HF_fs%values(b,b) ) ) 
              
      !         !under construction        
              CCSD_instance%ttau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                + 0.5*( CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) &
                  -CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j) )

              CCSD_instance%tau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                + CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) &
                  -CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j)

      !         ! write(*,*) CCSD_instance%Tdsame(a,b,i,j), "Tdsame"
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

  subroutine CCSD_loop_constructor()
      implicit none

      integer noc, nocs, nop, nops
      integer :: speciesId=1
      
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

  subroutine CCSD_T1T2_constructor()
      implicit none

      integer noc, nocs, nop, nops
      integer :: speciesId=1

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_T1T2_constructor: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops     

      !
      if (allocated(CCSDT1T2%Tai)) deallocate (CCSDT1T2%Tai)
      allocate(CCSDT1T2%Tai(noc-nop,nop))
      CCSDT1T2%Tai=0.0_8
      
      !
      if (allocated(CCSDT1T2%Tabij)) deallocate (CCSDT1T2%Tabij)
      allocate(CCSDT1T2%Tabij(noc-nop,noc-nop,nop,nop))
      CCSDT1T2%Tabij=0.0_8

  end subroutine CCSD_T1T2_constructor

  subroutine CCSD_loop()
      implicit none

      integer noc, nocs, nop, nops
      
      integer :: a, b, e, i, j, f, m, n
      real(8) :: ccsdE=0.0_8
      real(8) :: convergence = 1.0D-8
      real(8) :: prev_ccsdE
      real(8) :: tmp_ccsdE
      integer :: speciesId=1

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops
      write(*, "(A,I4,A,I4,A,I4,A,I4)") "CCSD_loop: noc=", noc, "nocs=", nocs, "nop=", nop, "nops=", nops

      
      do while (convergence >= 1.0D-8)

        prev_ccsdE = ccsdE

        call CCSD_loop_constructor()
        call CCSD_T1T2_constructor()
          
        !intermediates loop

        ! CCSDloop%Fac
        do a=nop+1, noc
          do e=nop+1, noc

            CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
              +(1 - logic2dbl(a==e))*Allspecies(speciesId)%HF_fs%values(a,e)
            
            do m=1, nop
              CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
                + (-0.5*Allspecies(speciesId)%HF_fs%values(m,e)*CCSD_instance%Tssame(a-nop,m))
              do f=nop+1, noc
                CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
                  + CCSD_instance%Tssame(f-nop,m)*spints(speciesId)%valuesp(m,a,f,e)
                do n=1, nop
                  CCSDloop%Fac(a-nop,e-nop) = CCSDloop%Fac(a-nop,e-nop) &
                    + (-0.5*CCSD_instance%ttau(a-nop,f-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f))
                  ! write(*,*) a,e,CCSDloop%Fac(a,e)
                end do
              end do
            end do
          end do
        end do
        ! write(*,*) CCSD_instance%ttau

        ! CCSDloop%Fki
        do m=1, nop
          do i=1, nop
             
            CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
              + (1 - logic2dbl(m==i))*Allspecies(speciesId)%HF_fs%values(m,i)
            
            do e=nop+1, noc
              CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
                 + 0.5*CCSD_instance%Tssame(e-nop,i)*Allspecies(speciesId)%HF_fs%values(m,e)
              do n=1, nop
                CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
                  + CCSD_instance%Tssame(e-nop,n)*spints(speciesId)%valuesp(m,n,i,e)
                do f=nop+1, noc
                  CCSDloop%Fki(m,i) = CCSDloop%Fki(m,i) &
                    + 0.5*CCSD_instance%ttau(e-nop,f-nop,i,n)*spints(speciesId)%valuesp(m,n,e,f)
                  ! write(*,*) a,e,CCSDloop%Fki(a,e)
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
                  + CCSD_instance%Tssame(f-nop,n)*spints(speciesId)%valuesp(m,n,e,f)
                ! write(*,*) a,e,CCSDloop%Fkc_aa(a,e)
              end do
            end do
          end do
        end do

        !CCSDloop%Wklij
        do m=1, nop
          do n=1, nop
            do i=1, nop
              do j=1, nop

                CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                  + spints(speciesId)%valuesp(m,n,i,j)
                do e=nop+1, noc
                  CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                    + (CCSD_instance%Tssame(e-nop,j)*spints(speciesId)%valuesp(m,n,i,e) &
                        -CCSD_instance%Tssame(e-nop,i)*spints(speciesId)%valuesp(m,n,j,e))
                  do f=nop+1, noc
                    CCSDloop%Wklij(m,n,i,j) = CCSDloop%Wklij(m,n,i,j) &
                      + 0.25*CCSD_instance%tau(e-nop,f-nop,i,j)*spints(speciesId)%valuesp(m,n,e,f)
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
                    + (-CCSD_instance%Tssame(b-nop,m)*spints(speciesId)%valuesp(a,m,e,f) &
                        +CCSD_instance%Tssame(a-nop,m)*spints(speciesId)%valuesp(b,m,e,f))
                  do n=1, nop
                    CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) = CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) &
                      + 0.25*CCSD_instance%tau(a-nop,b-nop,m,n)*spints(speciesId)%valuesp(m,n,e,f)
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
                   + CCSD_instance%Tssame(f-nop,j)*spints(speciesId)%valuesp(m,b,e,f)
                end do
                do n=1, nop
                  CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                    - CCSD_instance%Tssame(b-nop,n)*spints(speciesId)%valuesp(m,n,e,j)
                  do f=nop+1, noc
                    CCSDloop%Wkbcj(m,b-nop,e-nop,j) = CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                      - ((0.5*CCSD_instance%Tdsame(f-nop,b-nop,j,n) &
                          + CCSD_instance%Tssame(f-nop,j)*CCSD_instance%Tssame(b-nop,n))*spints(speciesId)%valuesp(m,n,e,f))
                    ! write(*,*) m,n,i,j,CCSDloop%Wkbcj(m,n,i,j)
                  end do
                end do
              end do
            end do
          end do
        end do

        ! CCSD Energy
        tmp_ccsdE=0.0_8
        do i=1, nop
          do a=nop+1, noc
            tmp_ccsdE = tmp_ccsdE + Allspecies(speciesId)%HF_fs%values(i,a)*CCSD_instance%Tssame(a-nop,i)
            do j=1, nop
              do b=nop+1, noc
                tmp_ccsdE = tmp_ccsdE + (0.25*spints(speciesId)%valuesp(i,j,a,b) *CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                  + 0.5*spints(speciesId)%valuesp(i,j,a,b)*CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j))
              end do
            end do
          end do
        end do
        ccsdE = tmp_ccsdE

        convergence = abs( ccsdE - prev_ccsdE )

        write (*,*) ccsdE, "CCSD Energy ", prev_ccsdE, "previous Energy" 
        write (*,*) convergence, "Convergence " 

        call CCSD_T1T2()

        if (convergence > 100) then 
          stop "test"
        end if

      end do
      
  end subroutine CCSD_loop

  subroutine CCSD_T1T2()
      implicit none

      integer noc, nocs, nop, nops
      integer :: a, b, e, f, i, j, m, n
      integer :: speciesId=1

      noc = Allspecies(speciesId)%noc
      ! nocs = CoupledCluster_instance%nocs
      nop = Allspecies(speciesId)%nop
      ! nops = CoupledCluster_instance%nops

      ! T^{a}_{i}D^{a}_{i} = ...
      do a=nop+1, noc
        do i=1, nop
          CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) + Allspecies(speciesId)%HF_fs%values(i,a)
          do e=nop+1, noc
            CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
              + CCSD_instance%Tssame(e-nop,i)*CCSDloop%Fac(a-nop,e-nop)
          end do
          do m=1, nop
            CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
              + (-CCSD_instance%Tssame(a-nop,m)*CCSDloop%Fki(m,i))
            do e=nop+1, noc
              CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
                + CCSD_instance%Tdsame(a-nop,e-nop,i,m)*CCSDloop%Fkc_aa(m,e-nop)
              do f=nop+1, noc
                CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
                  + (-0.5*CCSD_instance%Tdsame(e-nop,f-nop,i,m)*spints(speciesId)%valuesp(m,a,e,f))
              end do
              do n=1, nop
                CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
                  + (-0.5*CCSD_instance%Tdsame(a-nop,e-nop,m,n)*spints(speciesId)%valuesp(n,m,e,i))
              end do
            end do
          end do
          do n=1,nop
            do f=nop+1, noc
              CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i) &
                + (-CCSD_instance%Tssame(f-nop,n)*spints(speciesId)%valuesp(n,a,i,f))
            end do
          end do
          CCSDT1T2%Tai(a-nop,i) = CCSDT1T2%Tai(a-nop,i)/CCSDinit%Dai(a,i)
          CCSD_instance%Tssame(a-nop,i) = CCSDT1T2%Tai(a-nop,i)
          ! write(*,*) a,i,CCSD_instance%Tssame(a,i),CCSDT1T2%Tai(a,i)
        end do
      end do

      ! T^{ab}_{ij}D^{ab}_{ij} = ...
      do a=nop+1, noc
         do b=nop+1, noc
            do i=1, nop
               do j=1, nop
                  CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                    + spints(speciesId)%valuesp(i,j,a,b) !A
                  ! 1er ciclo
                  do e=nop+1, noc
                     CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                      + (CCSD_instance%Tdsame(a-nop,e-nop,i,j)*CCSDloop%Fac(b-nop,e-nop) &
                        -CCSD_instance%Tdsame(b-nop,e-nop,i,j)*CCSDloop%Fac(a-nop,e-nop)) !B
                     CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                      + (CCSD_instance%Tssame(e-nop,i)*spints(speciesId)%valuesp(a,b,e,j) &
                        -CCSD_instance%Tssame(e-nop,j)*spints(speciesId)%valuesp(a,b,e,i)) !G
                     do f=nop+1, noc
                        CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*CCSD_instance%tau(e-nop,f-nop,i,j)*CCSDloop%Wabcd(a-nop,b-nop,e-nop,f-nop) !D
                     end do
                     do m=1, nop
                        CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*CCSD_instance%Tdsame(a-nop,e-nop,i,j)*CCSD_instance%Tssame(b-nop,m)*CCSDloop%Fkc_aa(m,e-nop) &
                              +0.5*CCSD_instance%Tdsame(b-nop,e-nop,i,j)*CCSD_instance%Tssame(a-nop,m)*CCSDloop%Fkc_aa(m,e-nop)) !B'
                     end do
                  end do
                  ! 2do ciclo
                  do m=1, nop
                     CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                      + (-CCSD_instance%Tdsame(a-nop,b-nop,i,m)*CCSDloop%Fki(m,j) &
                          +CCSD_instance%Tdsame(a-nop,b-nop,j,m)*CCSDloop%Fki(m,i)) !C
                     CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                      + (-CCSD_instance%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,i,j) &
                          +CCSD_instance%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,i,j)) !H
                     do n=1, nop
                        CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                          + 0.5*CCSD_instance%tau(a-nop,b-nop,m,n)*CCSDloop%Wklij(m,n,i,j) !E
                     end do
                     do e=nop+1, noc
                        CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                          + (-0.5*CCSD_instance%Tdsame(a-nop,b-nop,i,m)*CCSD_instance%Tssame(e-nop,j)*CCSDloop%Fkc_aa(m,e-nop) &
                              +0.5*CCSD_instance%Tdsame(a-nop,b-nop,j,m)*CCSD_instance%Tssame(e-nop,i)*CCSDloop%Fkc_aa(m,e-nop)) !C'
                        CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                          + CCSD_instance%Tdsame(a-nop,e-nop,i,m)*CCSDloop%Wkbcj(m,b-nop,e-nop,j) &
                            - CCSD_instance%Tssame(e-nop,i)*CCSD_instance%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,e,j) &
                              -CCSD_instance%Tdsame(a-nop,e-nop,j,m)*CCSDloop%Wkbcj(m,b-nop,e-nop,i) &
                                + CCSD_instance%Tssame(e-nop,j)*CCSD_instance%Tssame(a-nop,m)*spints(speciesId)%valuesp(m,b,e,i) &
                                  -CCSD_instance%Tdsame(b-nop,e-nop,i,m)*CCSDloop%Wkbcj(m,a-nop,e-nop,j) &
                                    - CCSD_instance%Tssame(e-nop,i)*CCSD_instance%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,e,j) &
                                      + CCSD_instance%Tdsame(b-nop,e-nop,j,m)*CCSDloop%Wkbcj(m,a-nop,e-nop,i) &
                                       - CCSD_instance%Tssame(e-nop,j)*CCSD_instance%Tssame(b-nop,m)*spints(speciesId)%valuesp(m,a,e,i) !F
                     end do
                  end do
                  ! Make denominator array D^{ab}_{ij} = F_{ii}+F_{jj}-F_{a,a}-F_{b,b}
                  CCSDT1T2%Tabij(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j) &
                    /(Allspecies(speciesId)%HF_fs%values(i,i)+Allspecies(speciesId)%HF_fs%values(j,j)-Allspecies(speciesId)%HF_fs%values(a,a) &
                       -Allspecies(speciesId)%HF_fs%values(b,b))
                  
                  CCSD_instance%Tdsame(a-nop,b-nop,i,j) = CCSDT1T2%Tabij(a-nop,b-nop,i,j)

                  CCSD_instance%ttau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                    + 0.5*(CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) - CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j))
                  CCSD_instance%tau(a-nop,b-nop,i,j) = CCSD_instance%Tdsame(a-nop,b-nop,i,j) &
                    + CCSD_instance%Tssame(a-nop,i)*CCSD_instance%Tssame(b-nop,j) - CCSD_instance%Tssame(b-nop,i)*CCSD_instance%Tssame(a-nop,j)
               ! write(*,*) a,b,i,j,CCSD_instance%Tdsame(a,b,i,j),CCSDT1T2%Tabij(a,b,i,j)
               end do
            end do
         end do
      end do
      
  end subroutine CCSD_T1T2

  subroutine CCSD_run()
      implicit none
      
      call CCSD_init()
      call CCSD_loop()
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