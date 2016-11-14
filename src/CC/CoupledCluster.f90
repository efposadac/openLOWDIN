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
      integer(8) :: num_species
      real(8) :: CCSD_ones_Energy
      real(8) :: HF_energy
      real(8) :: HF_puntualInteractionEnergy
      real(8) :: MP2_EnergyCorr
      real(8), allocatable :: CCSD_E_intra(:)
      real(8), allocatable :: CCSD_E_inter(:)
      type(Vector) :: HF_orbitals
      type(matrix) :: HF_orbitals_dmatrix

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
  type(CoupledCluster), allocatable :: Allspecies(:)
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
  subroutine CoupledCluster_init()
    implicit none
    
    ! Load matrices
    ! Build diagonal matrix
    call CoupledCluster_loadWaveFunction()

    ! Calculate MP2 and load energies
    call CoupledCluster_MP2()

    ! Load matrices from transformed integrals
    call CoupledCluster_load_PFunction()

  end subroutine CoupledCluster_init

  !>
  !! @brief Load HF wave-function matrices
  !!        Build in a diagonal matrix with the eigenvectors
  !! @author CAOM
  subroutine CoupledCluster_loadWaveFunction()
      implicit none
      
      character(50) :: arguments(2)
      integer :: speciesId
      integer(8) :: x, i
      integer(8) :: nao
      integer :: num_species

      num_species = MolecularSystem_getNumberOfQuantumSpecies()
      CoupledCluster_instance%num_species = num_species

      !! Open file for wave-function
      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      ! Load TotalEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_energy, &
        arguments=["TOTALENERGY"])

      ! Load PuntualInteractionEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_puntualInteractionEnergy, &
        arguments=["PUNTUALINTERACTIONENERGY"])

      ! allocate array for many results by species
      if (allocated(Allspecies)) deallocate(Allspecies)
      allocate(Allspecies(num_species))

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

        ! call Vector_constructor (CoupledCluster_instance%HF_ff,nao*2,0.0_8)
        call Vector_constructor (Allspecies(speciesId)%HF_ff,nao*2,0.0_8)

        do x=1, nao
          do i=1, 2
            Allspecies(speciesId)%HF_ff%values((x-1)*2+i) = CoupledCluster_instance%HF_orbitals%values(x)
          end do
        end do

        ! FS matrix
        call Matrix_diagonalConstructor (Allspecies(speciesId)%HF_fs, Allspecies(speciesId)%HF_ff)

        ! !print*,"fs%new"
        ! !write(*,*) Allspecies(speciesId)%HF_fs%values


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
      integer :: i, m, j
          
      call MollerPlesset_constructor( 2 )
      call MollerPlesset_run()
      ! call MollerPlesset_secondOrderCorrection()
      call MollerPlesset_show()
      ! allocation of MP2 energy to CoupledCluster class
      CoupledCluster_instance%MP2_EnergyCorr= MollerPlesset_instance%secondOrderCorrection

      ! allocation of energy contributions of MP2 correction energy and coupling energy

      call Vector_constructor(CoupledCluster_instance%MP2_ECorr, MollerPlesset_instance%numberOfSpecies)
      call Vector_constructor(CoupledCluster_instance%MP2_ECpCorr, &
        MollerPlesset_instance%numberOfSpecies * (MollerPlesset_instance%numberOfSpecies -1) / 2 )

      m = 0
      do i=1, MollerPlesset_instance%numberOfSpecies
     
        CoupledCluster_instance%MP2_ECorr%values(i) = MollerPlesset_instance%energyCorrectionOfSecondOrder%values(i)
        
        if (MollerPlesset_instance%numberOfSpecies > 1) then

          do j = i + 1, MollerPlesset_instance%numberOfSpecies
            
            m = m + 1
            CoupledCluster_instance%MP2_ECpCorr%values(m) = MollerPlesset_instance%energyOfCouplingCorrectionOfSecondOrder%values(m)

          end do

        end if
      
      end do
                
  end subroutine CoupledCluster_MP2

  subroutine CoupledCluster_load_PFunction()
      implicit none
      integer :: num_species
      integer :: i

      num_species = CoupledCluster_instance%num_species

      !for the intra-species energies in CC
      if (allocated(CoupledCluster_instance%CCSD_E_intra)) deallocate(CoupledCluster_instance%CCSD_E_intra)
      allocate(CoupledCluster_instance%CCSD_E_intra(num_species))

      !for the inter-species energies in CC
      if (allocated(CoupledCluster_instance%CCSD_E_inter)) deallocate(CoupledCluster_instance%CCSD_E_inter)
      allocate(CoupledCluster_instance%CCSD_E_inter(f(num_species)/(2*f(num_species-2)))) 

      !for the one-species matrices from transformed integrals for CC    
      if (allocated(spints)) deallocate(spints)
      allocate(spints(num_species))

      !for the two-species matrices from transformed integrals for CC
      if (allocated(spintm)) deallocate(spintm)
      allocate(spintm(f(num_species)/(2*f(num_species-2)))) ! nc = n!/(2!*(n-2)!)
      ! nc is a number of posible combinations if there are more than one species (n>1)

      do i=1, num_species
        ! print*, "load_PF. i: ", i, " num_species: ", num_species, "combinations: ", f(num_species)/(2*f(num_species-2))
        call CoupledCluster_pairing_function(i, num_species)
      end do
      
  end subroutine CoupledCluster_load_PFunction


  !>
  ! @brief load: Auxiliary matrices from transformed integrals of one and two species
  !        build: Coulomb integrals of one and two species
  !               Exchange integrals of one species
  ! @author CAOM
  subroutine CoupledCluster_pairing_function(speciesId, num_species)
      implicit none
      type(Tensor), pointer :: axVc1sp
      integer, intent(in) :: speciesId
      integer, intent(in) :: num_species

      integer :: i
      integer :: m=0
      integer :: p, q, r, s
      integer :: noc, nocs, nop
      real(8) :: v_a, v_b, xv_a, xv_b

      
      Allspecies(speciesId)%nop = MolecularSystem_getNumberOfParticles(speciesID)
      nop = Allspecies(speciesId)%nop

      Allspecies(speciesId)%noc = MolecularSystem_getTotalNumberOfContractions(speciesId)
      ! This information is necessary in other modules: 
      ! number of contraction to number of molecular orbitals
      Allspecies(speciesId)%noc = Allspecies(speciesId)%noc*2
      ! For simplicity here
      noc = Allspecies(speciesId)%noc
      ! print*, "in one species. noc: ", noc, " nop: ", nop

      if (allocated(spints(speciesId)%valuesp)) deallocate(spints(speciesId)%valuesp)
      allocate(spints(speciesId)%valuesp(noc,noc,noc,noc))

      spints(speciesId)%valuesp(:,:,:,:)=0.0_8

      ! Read transformed integrals from file
      ! call ReadTransformedIntegrals_readOneSpecies( speciesID, CoupledCluster_instance%MP2_axVc1sp)
      call Tensor_constructor(CoupledCluster_instance%MP2_axVc1sp, speciesID, isMolecular=.true.)
      print*, "end Tensor_constructor one species"

      ! pairing function
      ! same species
      do p=1, noc
        do q=1, noc
          do r=1, noc
            do s=1, noc
              v_a = Tensor_getValue(CoupledCluster_instance%MP2_axVc1sp, (p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2)
              v_b = Tensor_getValue(CoupledCluster_instance%MP2_axVc1sp, (p+1)/2,(s+1)/2,(q+1)/2,(r+1)/2)
            
              xv_a = v_a * logic2dbl(mod(p,2) == mod(r,2)) * logic2dbl(mod(q,2) == mod(s,2))
              xv_b = v_b * logic2dbl(mod(p,2) == mod(s,2)) * logic2dbl(mod(q,2) == mod(r,2))
              ! spints
              spints(speciesId)%valuesp(p,q,r,s) = xv_a - xv_b
              ! write (*,*) spints(speciesId)%valuesp(p,q,r,s)
            end do
          end do
        end do
      end do
      ! print*, "spints(speciesId) complete"

      ! ! If there are two or more different species

      if ( num_species>1 ) then

        do i = speciesId + 1, num_species

          m = m + 1
          ! print*, "inside interspecies loop. i=speciesId+1: ", i
          Allspecies(i)%nop = MolecularSystem_getNumberOfParticles(i)
          ! print*, "Allspecies(i)%nop: ", Allspecies(i)%nop

          Allspecies(i)%noc = MolecularSystem_getTotalNumberOfContractions(i)
          ! This information is necesarry in other modules: 
          ! number of contraction to number of molecular orbitals
          Allspecies(i)%noc = Allspecies(i)%noc*2
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
          !print*,"different species"
          do p=1, noc
            do q=1, nocs
              do r=1, noc
                do s=1, nocs
                  v_a = Tensor_getValue(CoupledCluster_instance%MP2_axVc2sp, (p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,nocs/2) !! Coulomb integrals
          
                  xv_a = v_a * logic2dbl(mod(p,2) == mod(r,2)) * logic2dbl(mod(q,2) == mod(s,2))
                  spintm(m)%valuesp(p,q,r,s) = xv_a
                  write (*,*) spintm(m)%valuesp(p,q,r,s)    
                end do
              end do
            end do
          end do

          call Tensor_destructor(CoupledCluster_instance%MP2_axVc2sp)

        end do
          
      end if

      call Tensor_destructor(CoupledCluster_instance%MP2_axVc1sp)
      
      print*, "fin CoupledCluster_pairing_function"

  end subroutine CoupledCluster_pairing_function

  ! Factorial function just used to know the combination if there are more than one species
  recursive function f(n) result(output)
      implicit none
      integer :: n
      integer :: output

      if (n<2) then 
        output=1
      else
        output = n*f(n-1)
      end if

  end function f

end module CoupledCluster_
