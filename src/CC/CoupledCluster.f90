!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!    http://www.qcc.unal.edu.co/
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Coupled Cluster module
!!        This module initializes all information necessary to make a calculation of Coupled Cluster over APMO approach (CC-APMO).
!! @author  Carlos Andres Ortiz Mahecha (CAOM) (caraortizmah@unal.edu.co)
!!
!! <b> Creation date : </b> 2016-10-26
!!
!! <b> History: </b>
!!
!!   - <tt> 2016-10-26 </tt>: (CAOM) ( caraortizmah@unal.edu.co )
!!        -# Development of Coupled Cluster (CC) module:
!!                This Program calls three modules that load information about of equations of CC-APMO.
!!   - <tt> data </tt>:  
!!
!!
!! @warning <em>  All characters and events in this module -- even those based on real source code -- are entirely fictional. </br>
!!                All celebrity lines are impersonated.....poorly. </br> 
!!                The following module contains corase language and due to it's cintent should not be viewed by anyone. </em>
!!
!!
!!
module CoupledCluster_
  use Vector_
  use Matrix_
  use MolecularSystem_
  use MPFunctions_
  use Tensor_
  use String_
  implicit none

  type, public :: CoupledCluster
      
      integer :: noc, nocs, nop, nops
      integer(8) :: num_species, counterID, finalID
      integer(8) :: num_intersp
      integer :: spintm_m
      integer :: times_intersp
      integer :: n_intersp(10), i_counterID(10), info_ispecies(20)
      real(8) :: CCSD_once_Energy
      real(8) :: CCSD_twice_Energy
      real(8) :: CCSD_total_Energy
      real(8) :: HF_energy
      real(8) :: HF_puntualInteractionEnergy
      real(8) :: MP2_EnergyCorr
      real(8), allocatable :: CCSD_E_intra(:)
      real(8), allocatable :: CCSD_E_inter(:)
      type(Vector) :: HF_orbitals
      type(matrix) :: HF_orbitals_dmatrix

      real(8), allocatable :: Tssame(:,:)
      real(8), allocatable :: Tdsame(:,:,:,:)
      real(8), allocatable :: tau(:,:,:,:)
      real(8), allocatable :: ttau(:,:,:,:)
      real(8), allocatable :: intau(:,:,:,:)
      real(8), allocatable :: chtau_a(:,:,:,:,:,:)
      real(8), allocatable :: chtau_b(:,:,:,:,:,:)

      type(matrix) :: HF_fs
      type(Tensor) :: MP2_axVc1sp
      type(Tensor) :: MP2_axVc2sp
      type(Vector) :: HF_ff
      type(vector) :: CCSDCorr
      type(vector) :: ECorr2dOr
      type(vector) :: ECpCorr2dOr
      type(vector) :: MP2_ECpCorr
      type(vector) :: MP2_ECorr
      logical :: isInstanced

  end type CoupledCluster

  ! Tensor is used to auxiliry matrices from integrals transformed
  type, public :: TensorCC

      real(8), allocatable :: valuesp(:,:,:,:)
      logical :: isInstanced
  end type TensorCC

  
  type(CoupledCluster), public :: CoupledCluster_instance
  ! for same species
  type(CoupledCluster), allocatable :: Allspecies(:)
  ! for different species
  type(CoupledCluster), public, allocatable :: Allinterspecies(:)
  type(TensorCC), public :: TensorCC_instance
  !values in a array of arrays of single species <ab||ij> a,b,i,j are alpha species
  type(TensorCC), allocatable :: spints(:)
  !values of multiple species <ab|ij> a,i are alpha while b,j are beta species
  type(TensorCC), allocatable :: spintm(:)

  
  character(50), private :: wfnFile = "lowdin.wfn"
  integer, private :: wfnUnit = 20

contains
  
  !This function is temporal
  double precision function logic2dbl(a)
     logical, intent(in) :: a

      if (a) then
         logic2dbl = 1.d0
      else
         logic2dbl = 0.d0
      end if
  end function logic2dbl

  !>
  !! @brief Constructor of the class
  !! @author Carlos Andres Ortiz-Mahecha (CAOM) 
  subroutine CoupledCluster_constructor()
    implicit none
    
    !type(Vector) :: this
    CoupledCluster_instance%isInstanced = .true.
    ! print *, "cc coonstructor"

    ! call ReadTransformedIntegrals_readTwoSpecies( 1, 2, this)

    ! print *, "end cc coonstructor"

  end subroutine CoupledCluster_constructor

  !>
  !! @brief Destructor of the class
  !! @author CAOM
  subroutine CoupledCluster_destructor()
    implicit none

    integer(8) :: i, num_species

    num_species = MolecularSystem_getNumberOfQuantumSpecies()

    ! if (allocated(spints(1)%valuesp)) deallocate(spints(1)%valuesp)
    !if (allocated(spints)) deallocate(spints)
    ! do i=1, num_species

    !   !call Matrix_destructor(Allspecies(i)%HF_fs)
    !   !call Vector_destructor(Allspecies(i)%HF_ff)
    !   if (allocated(spints(i)%valuesp)) deallocate(spints(i))!!print*,"deallocate spints"!
    
    ! end do

    ! if (allocated(spints)) !print*, "CoupledCluster_destructor"!deallocate(spints)
    
    ! if (allocated(Allspecies)) deallocate(Allspecies)

    CoupledCluster_instance%isInstanced = .false.

  end subroutine CoupledCluster_destructor

  !>
  !! @brief Initialize Coupled Cluster Theory for SD and D cases
  !! @author CAOM
  subroutine CoupledCluster_init(cc_full,same_species,inter_species)
    implicit none
    logical, intent(in) :: cc_full
    character(50), intent(in) :: same_species(10)
    character(50), intent(in) :: inter_species(10)

    integer :: i, ii, output
    integer :: i_s=0
    integer :: i_i=1
    integer :: counter_i
    integer :: times_i=0
    integer :: min, max, split
    integer :: s_species(10), i_species(20), maxm(10), minm(10)
    integer :: num_species, counterID
    logical :: default_s=.true.
    logical :: default_i=.true.
    logical :: deft_s=.true.
    logical :: deft_i=.true.
    logical :: excp=.false.
    logical :: order_1=.false.
    logical :: order_2=.false.
    character(50) :: inter_sp(20)
    
    !INPUT CONTROL Coupled Cluster
    !ccFull=
    if (.not.cc_full) then
      !print*, "cc_full False: ", cc_full
      do i=1, size(same_species)
        !ccSameSpecies
        if (same_species(i)/= "NONE") then
          i_s=i_s+1
          print*, "same_species: ", same_species(i)
          default_s=.false.
        end if
        do ii = 1, MolecularSystem_instance%numberOfQuantumSpecies
          if( trim(MolecularSystem_instance%species(ii)%symbol) == trim(same_species(i))) then
             output = ii
          end if
        end do
        s_species(i) = MolecularSystem_getSpecieID(MolecularSystem_instance%species(output)%name)

        if (s_species(i)/=0) deft_s=.false.

      end do

      ! print*, "s_species: ", s_species

      if (.not.default_s) then
        if (.not.deft_s) then
          print*, "not default_s: ", MolecularSystem_instance%species(1)%symbol, same_species(1), s_species(1)
          max = maxval(s_species,mask=s_species.lt.10)
          min = minval(s_species,mask=s_species.gt.0)
        else
          print*, "default_s"
          min = 1
          max = MolecularSystem_getNumberOfQuantumSpecies()
        end if
        print*, "min and max value for same species: ", min, max

      else

        min = 1
        max = MolecularSystem_getNumberOfQuantumSpecies()

      end if

      !All subroutines in CoupledCluster_init depends on the next variables (1/2)
      counterID = min
      num_species = max
      !**

      do i=1, size(inter_species)
        deft_i=.true.
        !ccInterSpecies
        if (inter_species(i)/= "NONE") then
          ! print*, "inter_species: ", inter_species(i)
          split = index(inter_species(i),".")
          if (split==0) excp=.true.
          inter_sp(i_i) = inter_species(i)(1:split-1)
          inter_sp(i_i+1) = inter_species(i)(split+1:)
          i_i=i_i+2
            
          default_i=.false.
          i_species((i*2)-1) = MolecularSystem_getSpecieID(inter_sp((i*2)-1))
          i_species(i*2) = MolecularSystem_getSpecieID(inter_sp(i*2))

          print*, "inter_species: ", inter_sp((i*2)-1)
          print*, "inter_species: ", inter_sp(i*2)

        else

          i_species((i*2)-1) = MolecularSystem_getSpecieID(inter_species(i))
          i_species(i*2) = MolecularSystem_getSpecieID(inter_species(i))

          ! print*, "inter_species: ", inter_species(i)

        end if

        ! print*, "inter_species ID: ", i_species((i*2)-1)
        ! print*, "inter_species ID: ", i_species(i*2)

      
        if (i_species((i*2)-1)/=0) deft_i=.false.

        if (.not.default_i) then
          if (.not.deft_i) then
            !if the principal interaction species is wrong (equal to 0) then ignore the order
            if (i_species((i*2)-1)==0) then
              minm(i) = 0
              maxm(i) = 0
            
            else if (i_species((i*2)-1)==i_species(i*2)) then
              minm(i) = 0
              maxm(i) = 0
              order_1 = .true.

            else

            !number of correct orders in inter-species input control options
            times_i=times_i+1
            
              if(i_species((i*2)-1)>i_species(i*2)) then
                maxm(times_i) = i_species((i*2)-1)
                minm(times_i) = i_species(i*2)
              
                ! if the seconday interaction species are wrong (equal to 0) then 
                ! principal interaction species will be with all quantum species
                if (minm(times_i)==0) then
                  maxm(times_i) = MolecularSystem_getNumberOfQuantumSpecies()
                  minm(times_i) = i_species((i*2)-1)
                  order_2 = .true.
                end if
              
              else

                minm(times_i) = i_species((i*2)-1)
                maxm(times_i) = i_species(i*2)
              
              end if
            
            
            end if
            print*, "min and max value for inter species: ", minm(i), maxm(i)
            ! print*, "times_i: ", times_i

          else

            maxm(i) = 0
            minm(i) = 0

          end if

        else

          maxm(i) = 0
          minm(i) = 0

        end if

        if (excp) print*, "Error: Check the input CONTROL, ccInterSpecies has a wrong parameter"
        if (order_1) print*, "Default ccInterSpecies:  order ", i, ": speciesID (", minm(times_i), &
        ") will interact with all quantum species in this input"
        if (order_2) then
          print*, "Default ccInterSpecies: Order ", i, " was rejected: speciesID (", i_species((i*2)-1), &
            ") has not interaction with speciesID (", i_species(i*2),")."
          print*, "In inter-species APMO-CC there is not interaction between quantum species of the same type"
        end if
      
      end do
      
      !All subroutines in CoupledCluster_init depends on the next variables (2/2)
      counter_i=1
      !times_i already has allocation
      CoupledCluster_instance%i_counterID = minm
      CoupledCluster_instance%n_intersp = maxm
      ! information for CCSD_run()
      CoupledCluster_instance%info_ispecies = i_species
      CoupledCluster_instance%times_intersp = times_i
      !**

    else

      !All subroutines in CoupledCluster_init depends on the next variables (default)
      !For same species
      counterID = 1
      num_species = MolecularSystem_getNumberOfQuantumSpecies()
      !For different species
      counter_i = 0!counterID + 1
      times_i = 0
      CoupledCluster_instance%i_counterID = 0
      CoupledCluster_instance%n_intersp = 0
      ! information for CCSD_run()
      CoupledCluster_instance%info_ispecies = 0
      CoupledCluster_instance%times_intersp = times_i
      !**

    end if

    ! stop "fin"

    ! CoupledCluster_loadWaveFunction() and CoupledCluster_MP2() work independently of 
    ! input control information of Coupled Cluster

    ! Load matrices
    ! Build diagonal matrix
    call CoupledCluster_loadWaveFunction()

    ! Calculate MP2 and load energies
    call CoupledCluster_MP2()

    ! Depends on input control information of Coupled Cluster
    ! Load matrices from transformed integrals
    call CoupledCluster_load_PFunction(counterID, num_species)

  end subroutine CoupledCluster_init

  !>
  !! @brief Load HF wave-function matrices
  !!        Build in a diagonal matrix with the eigenvectors
  !! @author CAOM
  subroutine CoupledCluster_loadWaveFunction()
      implicit none

      ! integer, intent(in) :: counterID
      ! integer, intent(in) :: num_species
      
      character(50) :: arguments(2)
      integer :: speciesId, lambda
      integer(8) :: x, i
      integer(8) :: nao
      integer :: num_species
      integer :: num_intersp

      num_species = MolecularSystem_getNumberOfQuantumSpecies()
      CoupledCluster_instance%num_species = num_species

      num_intersp = num_species*(num_species-1)
      CoupledCluster_instance%num_intersp = num_intersp

      !! Open file for wave-function
      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      ! Load TotalEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_energy, &
        arguments=["TOTALENERGY"])

      ! Load PuntualInteractionEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_puntualInteractionEnergy, &
        arguments=["PUNTUALINTERACTIONENERGY"])

      ! allocate array for many results by full number of species
      if (allocated(Allspecies)) deallocate(Allspecies)
      allocate(Allspecies(num_species))

      ! allocate array for many results by full number of interactions of interspecies
      ! if num_species=1 Allinterspecies can not be used
      if (allocated(Allinterspecies)) deallocate(Allinterspecies)
      allocate(Allinterspecies(num_species))

      ! All species
      do speciesId = 1, num_species ! number of species is the size of CoupledCluster_instance

        nao = MolecularSystem_getTotalNumberOfContractions(speciesId)

        arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesId))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = int(nao,4), unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = CoupledCluster_instance%HF_orbitals )

        ! Build amplitudes from diagonal matrix HF for all species
        call Matrix_diagonalConstructor (CoupledCluster_instance%HF_orbitals_dmatrix, &
            CoupledCluster_instance%HF_orbitals)

        ! FF vector

        lambda = MolecularSystem_getLambda(speciesId)
        ! call Vector_constructor (CoupledCluster_instance%HF_ff,nao*2,0.0_8)
        call Vector_constructor (Allspecies(speciesId)%HF_ff,nao*lambda,0.0_8)

        do x=1, nao
          do i=1, lambda
            Allspecies(speciesId)%HF_ff%values((x-1)*lambda+i) = CoupledCluster_instance%HF_orbitals%values(x)
          end do
        end do

        ! FS matrix
        call Matrix_diagonalConstructor (Allspecies(speciesId)%HF_fs, Allspecies(speciesId)%HF_ff)

        print*, "HF_fs%values"
        call matrix_show( Allspecies(speciesId)%HF_fs)!%values)


        ! Load initial guess for T2 from MP2 correction
        call Vector_constructor(CoupledCluster_instance%CCSDCorr, num_species) !CCSD
        call Vector_constructor(CoupledCluster_instance%ECorr2dOr, num_species) !MP2

      end do

      close(wfnUnit)
      
  end subroutine CoupledCluster_loadWaveFunction

  !>
  ! @brief load: MP2 energy
  !              Two Vectors that contain all of energies by species of correction and coupling energies
  ! @author CAOM
  subroutine CoupledCluster_MP2()
      implicit none

      integer :: num_species

      ! integer, intent(in) :: counterID
      ! integer, intent(in) :: num_species

      integer :: i, m, j
          
      call MollerPlesset_constructor( 2 )
      call MollerPlesset_run()
      ! call MollerPlesset_secondOrderCorrection()
      call MollerPlesset_show()
      ! allocation of MP2 energy to CoupledCluster class
      CoupledCluster_instance%MP2_EnergyCorr= MollerPlesset_instance%secondOrderCorrection

      num_species = CoupledCluster_instance%num_species

      ! allocation of energy contributions of MP2 correction energy and coupling energy

      call Vector_constructor(CoupledCluster_instance%MP2_ECorr, num_species)
      call Vector_constructor(CoupledCluster_instance%MP2_ECpCorr, &
        num_species * (num_species -1) / 2 )

      m = 0
      do i=1, num_species
     
        CoupledCluster_instance%MP2_ECorr%values(i) = MollerPlesset_instance%energyCorrectionOfSecondOrder%values(i)
        
        if (num_species > 1) then

          do j = i + 1, num_species
            
            m = m + 1
            CoupledCluster_instance%MP2_ECpCorr%values(m) = MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder%values(m)

          end do

        end if
      
      end do
                
  end subroutine CoupledCluster_MP2

  subroutine CoupledCluster_load_PFunction(counterID, num_species)
      implicit none

      integer, intent(in) :: counterID
      integer, intent(in) :: num_species
      
      integer :: times_i
      integer :: i_counterID(10)
      integer :: n_intersp(10)

      integer :: i, j
      integer :: cont=0
      integer :: m=0

      CoupledCluster_instance%finalID = num_species
      CoupledCluster_instance%counterID = counterID
      times_i = CoupledCluster_instance%times_intersp
      i_counterID = CoupledCluster_instance%i_counterID
      n_intersp = CoupledCluster_instance%n_intersp

      ! print*, "times_i", times_i
      if (times_i>0) then
        ! print*, "times>0"
        do i=1, times_i
          do j=i+1, times_i
            if (i_counterID(i) == i_counterID(j)) then
              if (n_intersp(i) < n_intersp(j)) then
                n_intersp(i) = n_intersp(j)
              end if
              n_intersp(j) = 0
              i_counterID(i) = 0
            end if
          end do
          cont = cont + (n_intersp(i) - i_counterID(i))
          ! print*, "n_intersp(i) - i_counterID(i)", n_intersp(i), i_counterID(i)
        end do
        ! print*, "cont:", cont
      else
        cont = 0
      end if

      !for the intra-species energies in CC
      if (allocated(CoupledCluster_instance%CCSD_E_intra)) deallocate(CoupledCluster_instance%CCSD_E_intra)
      allocate(CoupledCluster_instance%CCSD_E_intra(num_species))

      !for the inter-species energies in CC
      if (allocated(CoupledCluster_instance%CCSD_E_inter)) deallocate(CoupledCluster_instance%CCSD_E_inter)
      allocate(CoupledCluster_instance%CCSD_E_inter(cont)) 

      !for the one-species matrices from transformed integrals for CC    
      if (allocated(spints)) deallocate(spints)
      allocate(spints(num_species-(counterID-1)))

      !for the two-species matrices from transformed integrals for CC
      if (allocated(spintm)) deallocate(spintm)
      allocate(spintm(cont)) ! nc = (n*(n-1))/2
      ! nc is a number of posible combinations if there are more than one species (n>1)

      do i=counterID, num_species
        print*, "load_PF. i: ", i, " num_species: ", num_species, "combinations: ", (num_species*(num_species-1))/2
        call CoupledCluster_pairing_function(i, num_species)
      end do

      if (times_i>0) then
        do i=1, times_i
          i_counterID(i) = CoupledCluster_instance%i_counterID(i)
          n_intersp(i) = CoupledCluster_instance%n_intersp(i)
          print*, "load_PF. i: ", i, " times_i: ", times_i, "Total inter-species matrices: ", cont
          call CoupledCluster_pairing_function_interspecies(i_counterID(i), n_intersp(i), m)
          ! public counter for spintm matrix
          m = CoupledCluster_instance%spintm_m
        end do
        !At final, m should be equal to cont
      end if
      
  end subroutine CoupledCluster_load_PFunction

  !>
  ! @brief load: Auxiliary matrices from transformed integrals of one species
  !        build: Coulomb and Exchange integrals of one species
  ! @author CAOM
  subroutine CoupledCluster_pairing_function(speciesId, num_species)
      implicit none
      type(Tensor), pointer :: axVc1sp
      integer, intent(in) :: speciesId
      integer, intent(in) :: num_species

      integer :: p, q, r, s
      integer :: noc, nop, lambda
      real(8) :: v_a, v_b, xv_a, xv_b

      
      Allspecies(speciesId)%nop = MolecularSystem_getNumberOfParticles(speciesID)
      nop = Allspecies(speciesId)%nop

      Allspecies(speciesId)%noc = MolecularSystem_getTotalNumberOfContractions(speciesId)
      ! This information is necessary in other modules: 
      ! number of contraction to number of molecular orbitals
      lambda = MolecularSystem_getLambda(speciesId)
      Allspecies(speciesId)%noc = Allspecies(speciesId)%noc*lambda
      ! For simplicity here
      noc = Allspecies(speciesId)%noc
      ! print*, "in one species. noc: ", noc, " nop: ", nop

      if (allocated(spints(speciesId)%valuesp)) deallocate(spints(speciesId)%valuesp)
      allocate(spints(speciesId)%valuesp(noc,noc,noc,noc))

      spints(speciesId)%valuesp(:,:,:,:)=0.0_8

      ! Read transformed integrals from file
      ! call ReadTransformedIntegrals_readOneSpecies( speciesID, CoupledCluster_instance%MP2_axVc1sp)
      call Tensor_constructor(CoupledCluster_instance%MP2_axVc1sp, speciesID, isMolecular=.true.)
      print*, "end Tensor_constructor one species", speciesId, num_species

      ! pairing function
      ! same species
      ! print*,"same species"
      do p=1, noc
        do q=1, noc
          do r=1, noc
            do s=1, noc

              v_a = Tensor_getValue(CoupledCluster_instance%MP2_axVc1sp, (p+1)/lambda,(r+1)/lambda,(q+1)/lambda,(s+1)/lambda, noc/lambda)
              v_b = Tensor_getValue(CoupledCluster_instance%MP2_axVc1sp, (p+1)/lambda,(s+1)/lambda,(q+1)/lambda,(r+1)/lambda, noc/lambda)
            
              xv_a = v_a * logic2dbl(mod(p,lambda) == mod(r,lambda)) * logic2dbl(mod(q,lambda) == mod(s,lambda))
              xv_b = v_b * logic2dbl(mod(p,lambda) == mod(s,lambda)) * logic2dbl(mod(q,lambda) == mod(r,lambda))
              ! spints
              spints(speciesId)%valuesp(p,q,r,s) = xv_a - xv_b
              ! print*, "spints speciesId=2"
              ! write (*,*) spints(speciesId)%valuesp(p,q,r,s)
            end do
          end do
        end do
      end do
      ! print*, "spints(speciesId) complete"

      call Tensor_destructor(CoupledCluster_instance%MP2_axVc1sp)
      
      print*, "fin CoupledCluster_pairing_function"

  end subroutine CoupledCluster_pairing_function

    !>
  ! @brief load: Auxiliary matrices from transformed integrals of two species
  !        build: Coulomb integrals of two species
  !               There are not Exchange integrals of two species
  ! WARNING:  Do not confuse the order in input variables in ...load_PFunction to
  !             CoupledCluster_pairing_function_interspecies(...)
  !           Names of variables could be confused due to the option default:
  !           In default opcion speciesId=counter_i and num_species=times_i
  !           If the user chooses the inter-species calculation speciesId/=counter_i 
  !              and num_species/=times_i
  !
  ! @author CAOM
  subroutine CoupledCluster_pairing_function_interspecies(speciesId, aux_max, m)
      implicit none

      integer, intent(in) :: speciesId
      integer, intent(in) :: aux_max
      integer, intent(inout) :: m

      integer :: i
      integer :: p, q, r, s
      integer :: noc, nocs, nop
      integer :: lambda1, lambda2
      real(8) :: v_a, xv_a

      
      Allspecies(speciesId)%nop = MolecularSystem_getNumberOfParticles(speciesID)
      nop = Allspecies(speciesId)%nop

      Allspecies(speciesId)%noc = MolecularSystem_getTotalNumberOfContractions(speciesId)
      ! This information is necessary in other modules: 
      ! number of contraction to number of molecular orbitals
      lambda1 = MolecularSystem_getLambda(speciesId)
      Allspecies(speciesId)%noc = Allspecies(speciesId)%noc*lambda1
      ! For simplicity here
      noc = Allspecies(speciesId)%noc

      ! Read transformed integrals from file
      ! ! If there are two or more different species

      do i = speciesId+1, aux_max

        m = m + 1
        nocs=0
        ! print*, "inside interspecies loop. i=speciesId+1: ", i
        Allspecies(i)%nop = MolecularSystem_getNumberOfParticles(i)
        ! print*, "Allspecies(i)%nop: ", Allspecies(i)%nop

        Allspecies(i)%noc = MolecularSystem_getTotalNumberOfContractions(i)
        ! This information is necesarry in other modules: 
        ! number of contraction to number of molecular orbitals
        lambda2 = MolecularSystem_getLambda(i)
        Allspecies(i)%noc = Allspecies(i)%noc*lambda2
        ! For simplicity here
        nocs = Allspecies(i)%noc
        ! print*, "i of allspecies: ",i
        ! print*, "noc: ", noc, " nocs: ",nocs
        if (allocated(spintm(m)%valuesp)) deallocate(spintm(m)%valuesp)
        allocate(spintm(m)%valuesp(noc,nocs,noc,nocs))
        spintm(m)%valuesp(:,:,:,:)=0.0_8

        ! Read transformed integrals from file
        ! call ReadTransformedIntegrals_readTwoSpecies( speciesId, OtherspeciesId, CoupledCluster_instance%MP2_axVc2sp)
        call Tensor_constructor(CoupledCluster_instance%MP2_axVc2sp, speciesID, otherSpeciesID=i, &
          isMolecular=.true.)
        print*, "end Tensor_constructor two species"
  
        !different species
        ! print*,"different species"
        do p=1, noc
          do q=1, nocs
            do r=1, noc
              do s=1, nocs
                v_a = Tensor_getValue(CoupledCluster_instance%MP2_axVc2sp, (p+1)/lambda1,(r+1)/lambda1,(q+1)/lambda2,(s+1)/lambda2, &
                    noc/lambda1,nocs/lambda2) !! Coulomb integrals
          
                xv_a = v_a * logic2dbl(mod(p,lambda1) == mod(r,lambda1)) * logic2dbl(mod(q,lambda2) == mod(s,lambda2))
                spintm(m)%valuesp(p,q,r,s) = xv_a
                ! write (*,*) spintm(m)%valuesp(p,q,r,s)    
              end do
            end do
          end do
        end do

        call Tensor_destructor(CoupledCluster_instance%MP2_axVc2sp)

      end do

      CoupledCluster_instance%spintm_m = m

      call Tensor_destructor(CoupledCluster_instance%MP2_axVc1sp)
      
      print*, "fin CoupledCluster_pairing_function_interspecies"
      ! stop "test"

  end subroutine CoupledCluster_pairing_function_interspecies

  ! ! Factorial function just used to know the combination if there are more than one species
  ! recursive function f(n) result(output)
  !     implicit none
  !     integer :: n
  !     integer :: output

  !     if (n<2) then 
  !       output=1
  !     else
  !       output = n*f(n-1)
  !     end if

  ! end function f

end module CoupledCluster_
