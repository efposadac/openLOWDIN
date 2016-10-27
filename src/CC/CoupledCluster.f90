module CoupledCluster_
  use Vector_
  use Matrix_
  use MolecularSystem_
  use MPFunctions_
  implicit none

  type, public :: CoupledCluster
      
      real(8) :: HF_energy
      real(8) :: HF_puntualInteractionEnergy
      real(8) :: MP2_EnergyCorr
      type(Vector) :: HF_orbitals
      type(matrix) :: HF_orbitals_dmatrix

      type(matrix) :: HF_fs
      type(Matrix) :: MP2_axMx1sp
      type(Matrix), pointer :: MP2_axMx2sp
      type(Vector) :: HF_ff
      type(vector) :: CCSDCorr
      type(vector) :: ECorr2dOr
      type(vector) :: ECpCorr2dOr
      type(vector) :: MP2_ECpCorr
      type(vector) :: MP2_ECorr
      logical :: isInstanced

  end type CoupledCluster

  ! Tensor is used to auxiliar matrices from integrals transformed
  type, public :: Tensor

      real(8), allocatable :: valuesp(:,:,:,:) !values of one species
  end type

  type(CoupledCluster), public :: CoupledCluster_instance

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
    !integer :: nspecies

    CoupledCluster_instance%isInstanced = .true.


  end subroutine CoupledCluster_constructor

  !>
  !! @brief Destructor of the class
  !! @author CAOM
  subroutine CoupledCluster_destructor()
    implicit none

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
    ! Point auxiliar matriz from MPFuntions.f90

    call CoupledCluster_MP2()

         
    !call CoupledCluster_pairing_function()

  end subroutine CoupledCluster_init

  !>
  !! @brief Load HF wave-function matrices
  !!        Build in a diagonal matrix with the eigenvectors
  !! @author Carlos Ortiz
  subroutine CoupledCluster_loadWaveFunction()
      implicit none
      
      character(50) :: arguments(2)
      integer :: speciesId
      integer(8) :: x, i
      integer :: nao, num_species

      num_species = MolecularSystem_getNumberOfQuantumSpecies()

      !! Open file for wave-function
      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      ! Load TotalEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_energy, &
        arguments=["TOTALENERGY"])

      ! Load PuntualInteractionEnergy
      call Vector_getFromFile(unit=wfnUnit, binary=.true., value=CoupledCluster_instance%HF_puntualInteractionEnergy, &
        arguments=["PUNTUALINTERACTIONENERGY"])

      ! One species
      do speciesId = 1, num_species ! number of species is the size of CoupledCluster_instance

        nao = MolecularSystem_getTotalNumberOfContractions(speciesId)

        arguments(2) = trim(MolecularSystem_getNameOfSpecie(speciesId))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = nao, unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = CoupledCluster_instance%HF_orbitals )

        ! Build amplitudes from diagonal matrix HF for all species
        call Matrix_diagonalConstructor (CoupledCluster_instance%HF_orbitals_dmatrix, &
            CoupledCluster_instance%HF_orbitals)

        ! FF vector
        call Vector_constructor (CoupledCluster_instance%HF_ff,nao*2,0.0_8)

        do x=1, nao
          do i=1, 2
            CoupledCluster_instance%HF_ff%values((x-1)*2+i) = CoupledCluster_instance%HF_orbitals%values(x)
          end do
        end do

        ! FS matrix
        call Matrix_diagonalConstructor (CoupledCluster_instance%HF_fs, CoupledCluster_instance%HF_ff)

        ! Load initial guess for T2 from MP2 correction
        call Vector_constructor(CoupledCluster_instance%CCSDCorr, num_species) !CCSD
        call Vector_constructor(CoupledCluster_instance%ECorr2dOr, num_species) !MP2

      end do
      
  end subroutine CoupledCluster_loadWaveFunction

  !>
  ! @brief load: MP2 energy
  !              Two Vectors that contain all of energies by species of correction and coupling energies
  !              Point auxiliar matrix of transformed integrals of one and two species from MPFunctions.f90 to here
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
      
      ! Point auxiliar matrix of integrals transformed from MollerPlesset class to CoupledCluster class:
      
      !CoupledCluster_instance%MP2_axMx1sp => MollerPlesset_instance%auxMatrix1sp
      !CoupledCluster_instance%MP2_axMx2sp => MollerPlesset_instance%auxMatrix2sp
            
  end subroutine CoupledCluster_MP2

  subroutine CoupledCluster_pairing_function(speciesId)
      implicit none
      integer, intent(in) :: speciesId
      integer :: OtherspeciesId, i
      integer :: p, q, r, s
      integer :: num_species
      integer :: noc, nocs
      real(8) :: v1, v2, v_a, v_b, xv_a, xv_b
      character(30) :: nameOtherSpId
      ! real(8), allocatable ::

      !array of matrices for auxiliar matrices from integrals transformed

      type(Tensor), allocatable :: spints(:)

      num_species = MolecularSystem_getNumberOfQuantumSpecies()

      if (allocated(spints)) deallocate(spints)
      allocate(spints(num_species))

      ! do speciesId = 1, num_species

      noc = MolecularSystem_getTotalNumberOfContractions(speciesId)
      ! number of contraction to nomber of molecular orbitals
      noc = noc*2

      if (allocated(spints(speciesId)%valuesp)) deallocate(spints(speciesId)%valuesp)
      allocate(spints(speciesId)%valuesp(noc,noc,noc,noc))

      spints(speciesId)%valuesp(:,:,:,:)=0.0_8

      ! Read transformed integrals from file
      call ReadTransformedIntegrals_readOneSpecies( speciesID, CoupledCluster_instance%MP2_axMx1sp)

      !! pairing function
      !same species
      do p=1, noc
        do q=1, noc
          do r=1, noc
            do s=1, noc
              v1 = IndexMap_tensorR4ToVector((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,noc/2) !! integrales de Coulomb
              v_a= CoupledCluster_instance%MP2_axMx1sp%values(v1, 1)
              v2 = IndexMap_tensorR4ToVector((p+1)/2,(s+1)/2,(q+1)/2,(r+1)/2,noc/2) !! integrales de intercambio
              v_b= CoupledCluster_instance%MP2_axMx1sp%values(v2, 1)
              xv_a = v_a * logic2dbl(mod(p,2) == mod(r,2)) * logic2dbl(mod(q,2) == mod(s,2))
              xv_b = v_b * logic2dbl(mod(p,2) == mod(s,2)) * logic2dbl(mod(q,2) == mod(r,2))
              spints(speciesId)%valuesp(p,q,r,s) = xv_a - xv_b
              write (*,*) spints(speciesId)%valuesp(p,q,r,s)
            end do
          end do
        end do
      end do
      
  end subroutine CoupledCluster_pairing_function


end module CoupledCluster_
