!!******************************************************************************
!!
!!		HERE WE GO !!
!!
!!******************************************************************************

module CoupledCluster_
  use Exception_
  use Matrix_
  use Vector_
  use Coupled_
  use ReadTransformedIntegrals_
  use MolecularSystem_
  use String_
  use IndexMap_
  implicit none

  !>
  !! @brief Coupled Interaction Module, works in spin orbitals
  !!
  !! @author Alejandro
  !!
  !! Optimization of CCSD module
  !! CCSD Interespecies
  !! CCSD(T) Code 
  !! @author Carlos (CAOM)
  !!
  !! <b> Creation data : </b> 2015
  !!
  !! <b> History change: </b>
  !! 
  !! <li> Optimization of CC: January 2016</li>
  !!
  !<
  type, public :: CoupledCluster
     logical :: isInstanced
     integer :: numberOfSpecies
     type(matrix) :: hamiltonianMatrix
!     integer :: numberOfCouplings
     type(vector) :: numberOfOccupiedOrbitals
     type(vector) :: numberOfOrbitals
     type(Vector) :: energyCorrectionOfSecondOrder
     type(Vector) :: mp2Coupling
     type(vector) :: eigenvalues
     type(vector) :: lambda !!Number of particles per orbital, module only works for 1 or 2 particles per orbital
     real(8) :: mp2Correction
     real(8) :: ccsdCorrection
     real(8) :: ccsdTwoParticlesCorrection
     real(8) :: secondOrderCorrection
!     real(8) :: coupledClusterSDCorrection
     type(vector) :: coupledClusterSDCorrection
     type(matrix), allocatable :: fourCenterIntegrals(:,:)
     type(matrix), allocatable :: twoCenterIntegrals(:)
!     type(coupled), allocatable :: couplings(:)
     real(8) :: totalEnergy

     character(20) :: level


     real(8), allocatable :: Tstest(:,:),Tdtest(:,:,:,:),Wmbejtest(:,:,:,:)

  end type CoupledCluster




  type, public :: HartreeFock
        real(8) :: totalEnergy
        real(8) :: puntualInteractionEnergy
        type(matrix) :: coefficientsofcombination 
        type(matrix) :: HcoreMatrix 
  end type HartreeFock
  
  type(CoupledCluster) :: CoupledCluster_instance
  type(HartreeFock) :: HartreeFock_instance

  public :: &
       CoupledCluster_constructor, &
       CoupledCluster_destructor, &
!       CoupledCluster_getTotalEnergy, &
       CoupledCluster_run, &
       CoupledCluster_show

  private		
contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine CoupledCluster_constructor(level)
    implicit none
    character(*) :: level

!    integer :: numberOfSpecies
!    integer :: i,j,m,n,p,q,c
!    integer :: isLambdaEqual1
!    type(vector) :: order
!    type(vector) :: occupiedCode
!    type(vector) :: unoccupiedCode
!    real(8) :: totalEnergy
!
!    character(50) :: wfnFile
!    integer :: wfnUnit
!    character(50) :: nameOfSpecie
!    integer :: numberOfContractions
!    character(50) :: arguments(2)
!
!    wfnFile = "lowdin.wfn"
!    wfnUnit = 20
!
!    !! Open file for wavefunction
!    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
!
!    !! Load results...
!    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%totalEnergy, &
!         arguments=["TOTALENERGY"])
!    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%puntualInteractionEnergy, &
!         arguments=["PUNTUALINTERACTIONENERGY"])
!
!    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
!    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(i)
!
!    do i=1, numberOfSpecies
!        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
!        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
!
!        arguments(2) = nameOfSpecie
!        arguments(1) = "HCORE"
!        HartreeFock_instance%HcoreMatrix  = &
!                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
!                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
!
!        arguments(1) = "COEFFICIENTS"
!        HartreeFock_instance%coefficientsofcombination = &
!                  Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
!                  columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
!
!!         arguments(1) = "ORBITALS"
!!         call Vector_getFromFile( elementsNum = numberOfContractions, &
!!              unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
!!              output = eigenValues )
!
!
!
!    end do
!
    integer :: numberOfContractionss,nopp,numberOfParticless,nocc,i
!!*
!	 do i=1, MolecularSystem_getNumberOfQuantumSpecies()
!!*
!    numberOfContractionss = MolecularSystem_getTotalNumberOfContractions( i )
!!*
!     numberOfParticless = MolecularSystem_getNumberOfParticles(i)
!!*
!    nopp=numberOfParticless
!    nocc=numberOfContractionss
!!*
!!*
!!*
!    if (allocated(CoupledCluster_instance%Tstest)) deallocate(CoupledCluster_instance%Tstest)
!    allocate(CoupledCluster_instance%Tstest(nocc-nopp,nopp)) ! 01f
!!*
!!*
!!*    if (allocated(CoupledCluster_instance%Ts)) deallocate(CoupledCluster_instance%Ts)
!!*    allocate(CoupledCluster_instance%Ts(nocc-nopp,nopp)) ! 01f
!!*
 !   CoupledCluster_instance%Tstest = 136
!!*
!    end do


    CoupledCluster_instance%isInstanced=.true.
    CoupledCluster_instance%level=level

! Not necessary... maybe    
!    CoupledCluster_instance%numberOfCouplings=0

!    call Vector_constructor (CoupledCluster_instance%numberOfOccupiedOrbitals, numberOfSpecies)
!    call Vector_constructor (CoupledCluster_instance%numberOfOrbitals, numberOfSpecies)
!    call Vector_constructor (CoupledCluster_instance%lambda, numberOfSpecies)

!    do i=1, numberOfSpecies
       !! We are working in spin orbitals not in spatial orbitals!
!       CoupledCluster_instance%lambda%values(i) = MolecularSystem_getLambda( i )
!       CoupledCluster_instance%numberOfOccupiedOrbitals%values(i)=MolecularSystem_getOcupationNumber( i )* CoupledCluster_instance%lambda%values(i)
!       CoupledCluster_instance%numberOfOrbitals%values(i)=MolecularSystem_getTotalNumberOfContractions( i )* CoupledCluster_instance%lambda%values(i)
       !!Uneven occupation number = alpha
       !!Even occupation number = beta     
!    end do

!    close(wfnUnit)

  end subroutine CoupledCluster_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine CoupledCluster_destructor()
    implicit none

    CoupledCluster_instance%isInstanced=.false.

  end subroutine CoupledCluster_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CoupledCluster_show()
    implicit none
    type(CoupledCluster) :: this
    integer :: i,j,k,l
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
 
    if ( CoupledCluster_instance%isInstanced ) then

       print *,""
       print *," POST HARTREE-FOCK CALCULATION"
       print *," COUPLED CLUSTER THEORY:"
       print *,"=============================="
       print *,""
       write (6,"(T10,A30, A5)") "LEVEL = ", CoupledCluster_instance%level
       write (6,"(T10,A30, F20.12)") "HF ENERGY = ", HartreeFock_instance%totalEnergy
       write (6,"(T10,A30, F20.12)") "MP2 CORR. ENERGY = ", CoupledCluster_instance%secondOrderCorrection
       write (6,"(T10,A30, F20.12)") "CCSD CORR. ENERGY = ", CoupledCluster_instance%ccsdCorrection
       write (6,"(T10,A30, F20.12)") "============================================================"
       write (6,"(T10,A30, F20.12)") "Total Energy (HF+CCSD) = ", HartreeFock_instance%totalEnergy+CoupledCluster_instance%ccsdCorrection
    

	if ( numberOfSpecies > 1 ) then

       print *,""
       print *,""
       print *," COUPLED CLUSTER THEORY FOR DIFF. SPECIES"
       print *,"========================================="
       print *,""
       print *,""
       print *,""

                    k=0
		    l=0
                    do i=1, numberOfSpecies
                         do j=i+1, numberOfSpecies
                                   k=k+1
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "MP2 Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
					"} = ", CoupledCluster_instance%mp2Coupling%values(k)
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCSD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( i ) ), &
					"} = ",  CoupledCluster_instance%coupledClusterSDCorrection%values( i )
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCSD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( j ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
					"} = ",  CoupledCluster_instance%coupledClusterSDCorrection%values( j )
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCSD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
					"} = ", CoupledCluster_instance%ccsdTwoParticlesCorrection
                         end do
                    end do

       print *,""
       print *,""
				   write (6,'(T10,A30,F20.8)') "========================================================="
				   write (6,'(T10,A30,F20.8)') "MP2 Total Energy = ", HartreeFock_instance%totalEnergy+CoupledCluster_instance%secondOrderCorrection+CoupledCluster_instance%mp2Coupling%values 
				   write (6,'(T10,A30,F20.8)') "CCSD Total Energy = ", HartreeFock_instance%totalEnergy+CoupledCluster_instance%ccsdCorrection+CoupledCluster_instance%ccsdTwoParticlesCorrection 

	end if

    else 

    end if

  end subroutine CoupledCluster_show

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CoupledCluster_run()
    implicit none 
    integer :: m
    real(8), allocatable :: eigenValues(:) 

    select case ( trim(CoupledCluster_instance%level) )

    case ( "CCSD" )

       print *, ""
       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN CCSD CALCULATION           |"
       print *, "-----------------------------------------------"
       print *, ""

	print *, "  ________________________________________________ "
	print *, " |Initial Guess... From MP2 Calculation           |"
	print *, " |------------------------------------------------|"
	print *, " |Iteration of Coupled Cluster intermediates . . .|"
	print *, " |------------------------------------------------|"
	print *, " |Computing CCSD Energy...                        |"
	print *, " |________________________________________________|"
	print *, ""
	print *, "Patience is not the ability to wait, but the ability to keep a good attitude while waiting... Or just destroy the computer."

	call CoupledCluster_iterateIntermediates_SameSpecies()
	
  print *, "mi variable copartida Ts",  CoupledCluster_instance%Tstest(1,1) 
  print *, "mi variable copartida Td",  CoupledCluster_instance%Tdtest(1,1,1,1) 
  print *, "mi variable copartida Wmbej",  CoupledCluster_instance%Wmbejtest(1,1,1,1) 

   call CoupledCluster_iterateIntermediates_DiffSpecies()
 
  print *, "mi variable copartida Ts",  CoupledCluster_instance%Tstest(1,1) *2
  print *, "mi variable copartida Td",  CoupledCluster_instance%Tdtest(1,1,1,1) *2
  print *, "mi variable copartida Wmbej",  CoupledCluster_instance%Wmbejtest(1,1,1,1) *2

       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END CCSD CALCULATION           |"
       print *, "==============================================="
       print *, ""

    case default

       call CoupledCluster_exception( ERROR, "Coupled interactor constructor", "Correction level not implemented")

    end select


  end subroutine CoupledCluster_run

 !>
  !! @brief Calculation of the intermediates
  !!
  !<
 subroutine CoupledCluster_iterateIntermediates_SameSpecies()
   implicit none

   integer :: numberOfSpecies
   integer :: a,b,c,d
   integer :: aa,bb,cc,dd
   integer :: i,j,k,l,p,q,r,s,h,t
   integer :: ii,jj,kk,ll,pp,qq,rr,ss,hh,tt
   integer(8) :: x,y,z,auxIndex,auxIndex2
   integer :: e,f,m,n,iii,jjj,aaa,bbb
   integer :: ee,mm,nn,fff,oo
   integer :: speciesID
   integer :: otherSpeciesID
   character(10) :: nameOfSpecie
   character(10) :: nameOfOtherSpecie
   integer :: electronsID
   integer :: numberOfParticles
   integer :: numberOfOtherSpecieParticles
   integer :: ocupationNumber
   integer :: ocupationNumberOfOtherSpecie
   integer :: numberOfContractions
   integer :: numberOfContractionsOfOtherSpecie
   integer(8) :: numberOfOrbitals
   integer(8) :: numberOfSpatialOrbitals
   integer :: numberOfOtherSpecieOrbitals
   integer :: numberOfOtherSpecieSpatialOrbitals
   integer, allocatable :: spin(:)
   integer, allocatable :: spatialOrbital(:)
   type(Vector) :: eigenValues, ff
   type(Vector) :: eigenValuesOfOtherSpecie, otherff
   type(Vector) :: coupledClusterValue
   type(Matrix) :: auxMatrix!   type(TransformIntegrals) :: repulsionTransformer
   real(8) :: lambda
   real(8) :: lambdaOfOtherSpecie
   real(8) :: kappa !positive or negative exchange
   real(8) :: charge
   real(8) :: TwoParticlesEnergy, ECCSD, DECC, OLDCC, auxECCSD, ECCSD1, ECCSD2
   real(8) :: otherSpecieCharge
   real(8) :: independentEnergyCorrection
   real(8) :: mp2CouplingCorrection
   real(8) :: auxVal,auxVal1,auxVal2,value1,value2
   real(8) :: auxVal_A
   real(8) :: auxVal_B
   type(Matrix) :: eigenVec, eigenVec1, Fs
   type(Matrix) :: eigenVecOtherSpecie, eigenVecOtherSpecie1, otherFs
!   real(8), allocatable :: ff(:)
! 01 febrero 2016
   real(8), allocatable :: Fae(:,:),Fmi(:,:),Fme(:,:),Dai(:,:),TsNew(:,:)
   real(8), allocatable :: Ts(:,:)
   real(8), allocatable :: auxTs(:,:),auxFae(:,:),auxFmi(:,:),auxFme(:,:),auxDai(:,:),auxTsNew(:,:)
   real(8), allocatable :: otherTs(:,:),otherFae(:,:),otherFmi(:,:),otherFme(:,:),otherDai(:,:),otherTsNew(:,:)
   real(8), allocatable :: Td(:,:,:,:), spinints(:,:,:,:), Dabij(:,:,:,:), taus(:,:,:,:), tau(:,:,:,:), Wmnij(:,:,:,:), Wabef(:,:,:,:), Wmbej(:,:,:,:), TdNew(:,:,:,:)
   real(8), allocatable :: otherTd(:,:,:,:), otherspinints(:,:,:,:), otherDabij(:,:,:,:), othertaus(:,:,:,:), auxspinints(:,:,:,:), auxTd(:,:,:,:)
   real(8), allocatable :: othertau(:,:,:,:), otherWmnij(:,:,:,:), otherWabef(:,:,:,:), otherWmbej(:,:,:,:), otherTdNew(:,:,:,:)
   real(8), allocatable :: auxWmnij(:,:,:,:), auxWabef(:,:,:,:), auxWmbej(:,:,:,:)
   real(8), allocatable :: auxTdNew(:,:,:,:), auxtau(:,:,:,:), auxtaus(:,:,:,:), auxDabij(:,:,:,:)
   type(Matrix), allocatable :: auxMatrix1(:,:)
   character(50) :: wfnFile
   character(50) :: arguments(2)
   integer :: wfnUnit
   !! 21 enero 2016
   integer :: noc
   integer :: nop
   integer :: nocs
   integer :: nops
!! 27 de enero 2016
   integer :: kro

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


!!******************************************************************************************************************************
!! Begin CCSD calculation same specie. 
!!******************************************************************************************************************************



!	speciesID=1 !!! FOR NOW, only electrons 
    do speciesID=1, numberOfSpecies
	numberOfParticles = MolecularSystem_getNumberOfParticles(speciesID)

    if ( numberOfParticles > 1 ) then

	lambda=MolecularSystem_getLambda(speciesID) !Particles per orbital
        kappa=MolecularSystem_getKappa(speciesID) !exchange sign
        charge=MolecularSystem_getCharge(speciesID)
    	ocupationNumber = MolecularSystem_getOcupationNumber(speciesID)
	numberOfContractions = MolecularSystem_getTotalNumberOfContractions( speciesID )
        numberOfOrbitals = numberOfContractions*lambda
        numberOfSpatialOrbitals = (numberOfOrbitals/lambda)*2 !Twice the dimension of spatial orbitals

!
         arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)

    	!! Read transformed integrals from file
        call ReadTransformedIntegrals_readOneSpecies( speciesID, auxMatrix)

         arguments(1) = "ORBITALS"
         call Vector_getFromFile( elementsNum = numberOfContractions, &
              unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
              output = eigenValues )

	call Matrix_diagonalConstructor (eigenVec, eigenValues) ! put MO energies in diagonal array (eigenVec) 
	call Vector_constructor (ff,numberOfContractions*2,0.0_8)

	a=0
	b=0
	do x=1, numberOfContractions
		a=a+1
		do i=1, 2
			b=b+1
			ff%values(b) = eigenValues%values(a)
		end do
	end do

	call Matrix_diagonalConstructor (Fs, ff)

	
	! Inital guess for T2, initial guess for the cluster amplitudes are the Moller-Plesset first-order perturbed wave function.
	! t_i^a=0 ; t_{ij}^{ab}= <ij||ab>/(E_i + E_j - E_a - E_b)

	call Vector_constructor( CoupledCluster_instance%coupledClusterSDCorrection, numberOfSpecies)
	call Vector_constructor( CoupledCluster_instance%energyCorrectionOfSecondOrder, numberOfSpecies)
   
   !! Aqui empieza MP2

	do a=1, ocupationNumber
		do b=1, ocupationNumber
			do i=ocupationNumber+1, numberOfContractions
				do j=i, numberOfContractions
					auxIndex = IndexMap_tensorR4ToVector(a,i,b,j, numberOfContractions)
					auxVal_A= auxMatrix%values(auxIndex, 1)

						if (  dabs( auxVal_A)  > 1.0E-10_8 ) then

                                                        if ( j>i ) then

                                                                if (a==b) then

                                                                        if( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                                                                                independentEnergyCorrection = independentEnergyCorrection + 2.0_8 *  auxVal_A**2.0  &
                                                                                * ( lambda  -  1.0_8 ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                                                - eigenValues%values(i) - eigenValues%values(j) )

                                                                        end if

                                                                else

                                                                        auxIndex = IndexMap_tensorR4ToVector(i, b, j, a, numberOfContractions )
                                                                        auxVal_B= auxMatrix%values(auxIndex, 1)

                                                                        independentEnergyCorrection = independentEnergyCorrection + 2.0_8 *  auxVal_A  &
                                                                        * ( lambda * auxVal_A  - auxVal_B ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                                        - eigenValues%values(i) - eigenValues%values(j) )

                                                                end if

                                                        else if ( a==b .and. i==j ) then

                                                                if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

                                                                        independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0_8  &
                                                                                * ( lambda - 1.0_8 ) / ( 2.0_8*( eigenValues%values(a)-eigenValues%values(i)))

                                                                end if

                                                        else

                                                                if ( abs( lambda  -  1.0_8 ) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
                                                                        independentEnergyCorrection = independentEnergyCorrection +  auxVal_A**2.0  &
                                                                                * ( lambda  - 1.0_8 ) / ( eigenValues%values(a) + eigenValues%values(b) &
                                                                                - eigenValues%values(i) - eigenValues%values(j) )
                                                                end if

                                                        end if
                                                end if


                                        end do
                                end do
                        end do

		end do

                CoupledCluster_instance%energyCorrectionOfSecondOrder%values(speciesID) = independentEnergyCorrection &
                        * ( ( MolecularSystem_getCharge( speciesID ) )**4.0_8 )

		!! Suma las correcciones de energia para especies independientes
		CoupledCluster_instance%secondOrderCorrection = sum( CoupledCluster_instance%energyCorrectionOfSecondOrder%values ) !MP2 Corr. Energy


!!! From scratch...

        print *, "numberOfContractions", numberOfContractions
	numberOfContractions=numberOfContractions*2
        print *, "numberOfContractions", numberOfContractions

!! 21 de enero 2016
  noc=numberOfContractions
  nop=numberOfParticles


        if (allocated(spinints)) deallocate (spinints)
        allocate(spinints(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        spinints(:,:,:,:) = 0.0_8

!! pasa de 4 indices a 1 (pairing function)
        do p=1, numberOfContractions
                do q=1, numberOfContractions
                        do r=1, numberOfContractions
                                do s=1, numberOfContractions
                                        value1 = IndexMap_tensorR4ToVector((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,numberOfContractions/2) !! integrales de Coulomb
                                        auxVal_A= auxMatrix%values(value1, 1)
                                        value2 = IndexMap_tensorR4ToVector((p+1)/2,(s+1)/2,(q+1)/2,(r+1)/2,numberOfContractions/2) !! integrales de intercambio
                                        auxVal_B= auxMatrix%values(value2, 1)
                                        auxVal1 = auxVal_A * (mod(p,2) == mod(r,2)) * (mod(q,2) == mod(s,2))
                                        auxVal2 = auxVal_B * (mod(p,2) == mod(s,2)) * (mod(q,2) == mod(r,2))
                                        spinints(p,q,r,s) = auxVal1 - auxVal2 !! p+1 o p-1? Revisar !! ecuacion 1
!					write (*,*) spinints(p,q,r,s)
                                end do
                        end do
                end do
        end do


!!!! Initial guesses T1 and T2



!variables compartidas

    if (allocated(CoupledCluster_instance%Tstest)) deallocate(CoupledCluster_instance%Tstest)
    allocate(CoupledCluster_instance%Tstest(noc-nop,nop)) ! 01f
    CoupledCluster_instance%Tstest(:,:) = 0.0_8


    if (allocated(CoupledCluster_instance%Tdtest)) deallocate(CoupledCluster_instance%Tdtest)
    allocate(CoupledCluster_instance%Tdtest(noc-nop,noc-nop,nop,nop)) ! 01f
    CoupledCluster_instance%Tdtest(:,:,:,:) = 0.0_8


    if (allocated(CoupledCluster_instance%Wmbejtest)) deallocate(CoupledCluster_instance%Wmbejtest)
    allocate(CoupledCluster_instance%Wmbejtest(nop,noc-nop,noc-nop,nop)) ! 01f
    CoupledCluster_instance%Wmbejtest(:,:,:,:) = 0.0_8


!   print *, "55555mi variable copartida ",  CoupledCluster_instance%Tstest(1,1) 

! if (allocated(Ts)) deallocate (Ts)
! allocate(Ts(noc-nop,nop))
! Ts(:,:) = 0.0_8
 

        !! Td solo corre sobre virtual,virtual,ocupado,ocupado
        if (allocated(Td)) deallocate (Td)
 !! 22 de enero 2016
 !!      allocate(Td(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(Td(noc-nop,noc-nop,nop,nop))
        Td(:,:,:,:) = 0.0_8

        if (allocated(taus)) deallocate (taus)
!! 22 enero 2016 
!!      allocate(taus(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(taus(noc-nop,noc-nop,nop,nop))
        taus(:,:,:,:) = 0.0_8

        if (allocated(tau)) deallocate (tau)
!! 22 enero 2016
!!      allocate(tau(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(tau(noc-nop,noc-nop,nop,nop))
        tau(:,:,:,:) = 0.0_8

        do a=numberOfParticles+1, numberOfContractions
                do b=numberOfParticles+1, numberOfContractions
                        do i=1, numberOfParticles
                                do j=1, numberOfParticles

! 01 febrero 2016

!CoupledCluster_instance%Tstest(:,:) = 0.0_8

                                        CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) = CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) + (spinints(i,j,a,b)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)))
                                        taus(a-nop,b-nop,i,j) = CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) + 0.5*(CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j))
                                        tau(a-nop,b-nop,i,j) = CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) + CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j)
!! 22 de enero 2016
                                !        Td(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + (spinints(i,j,a,b)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)))
                                !        taus(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + 0.5*(Ts(a-nop,i)*Ts(b-nop,j) - Ts(b-nop,i)*Ts(a-nop,j))
                                !        tau(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + Ts(a-nop,i)*Ts(b-nop,j) - Ts(b-nop,i)*Ts(a-nop,j)
!			                  write(*,*) Td(a,b,i,j), taus(a,b,i,j), tau(a,b,i,j)
                                end do
                        end do
                end do
        end do

!!!! End Initial Guesses

!!!! Make denominator arrays Dai, Dabij

        !! Equation 12 from Stanton

        if (allocated(Dai)) deallocate (Dai)
        allocate(Dai(numberOfContractions,numberOfContractions))
        Dai(:,:) = 0.0_8

        do a=numberOfParticles+1, numberOfContractions
                do i=1, numberOfParticles
                        Dai(a,i) = Fs%values(i,i) - Fs%values(a,a)
!			write(*,*) a,i,Dai(a,i)
                end do
        end do
        

!!!! End Denominators

!!!! Main Loop



  ECCSD = 0.0_8
  DECC = 1.0_8 
!  OLDCC = 0 

!   do while (OLDCC < 5)
  do while (DECC >= 1.0D-8)
!	OLDCC = OLDCC + 1
	OLDCC = ECCSD


	if (allocated(Fae)) deallocate (Fae)
!! 22 de enero 2016
!!      allocate(Fae(numberOfContractions,numberOfContractions))
        allocate(Fae(noc-nop,noc-nop))
   Fae=0.0_8
	if (allocated(Fmi)) deallocate (Fmi)
!! 21-enero-2016
!!      allocate(Fmi(numberOfContractions,numberOfContractions))
        allocate(Fmi(nop,nop))
	Fmi=0.0_8
	if (allocated(Fme)) deallocate (Fme)
!!        allocate(Fme(numberOfContractions,numberOfContractions))
          allocate(Fme(nop,noc-nop))
   Fme=0.0_8

	if (allocated(Wmnij)) deallocate (Wmnij)
!! 21-enero-2016
!!        allocate(Wmnij(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(Wmnij(nop,nop,nop,nop))
	Wmnij=0.0_8
	if (allocated(Wabef)) deallocate (Wabef)
!! 21-enero-2016
!!        allocate(Wabef(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(Wabef(noc-nop,noc-nop,noc-nop,noc-nop))
	Wabef=0.0_8
	if (allocated(Wmbej)) deallocate (Wmbej)
        allocate(Wmbej(nop,noc-nop,noc-nop,nop))
	Wmbej=0.0_8

	if (allocated(TsNew)) deallocate (TsNew)
!!	allocate(TsNew(numberOfContractions,numberOfContractions))
   allocate(TsNew(noc-nop,nop))
	TsNew=0.0_8
	if (allocated(TdNew)) deallocate (TdNew)
        allocate(TdNew(noc-nop,noc-nop,nop,nop))
	TdNew=0.0_8



!!!! Intermediates

	!! Equation 3
! 27 de enero 2016 
! eq OK

	do a=numberOfParticles+1, numberOfContractions
		do e=numberOfParticles+1, numberOfContractions
!! 27 de enero 2016
! eq ok 
!!       Fae(a-nop,e-nop) = (1 + (a==e))*Fs%values(a,e)
         if (a==e) then 
            kro = 1 
         else 
            kro = 0 
         end if
         
         Fae(a-nop,e-nop) = (1 - kro)*Fs%values(a,e)
			do m=1, numberOfParticles
				Fae(a-nop,e-nop) = Fae(a-nop,e-nop) + (-0.5*Fs%values(m,e)*CoupledCluster_instance%Tstest(a-nop,m))
				do f=numberOfParticles+1, numberOfContractions
					Fae(a-nop,e-nop) = Fae(a-nop,e-nop) + CoupledCluster_instance%Tstest(f-nop,m)*spinints(m,a,f,e)
					do n=1, numberOfParticles
						Fae(a-nop,e-nop) = Fae(a-nop,e-nop) + (-0.5*taus(a-nop,f-nop,m,n)*spinints(m,n,e,f))
!						write(*,*) a,e,Fae(a,e)
					end do
				end do
			end do
		end do
	end do

!	write(*,*) numberOfContractions
!	write(*,*) Fae

	!! Equation 4

	do m=1, numberOfParticles
		do i=1, numberOfParticles
!! 27 de enero 2016
! eq ok
!!       Fmi(m,i) = (1 + (m==i))*Fs%values(m,i)
         if (m==i) then 
            kro = 1 
         else 
            kro = 0 
         end if
         
         Fmi(m,i) = (1 - kro)*Fs%values(m,i)
			do e=numberOfParticles+1, numberOfContractions
				Fmi(m,i) = Fmi(m,i) + 0.5*CoupledCluster_instance%Tstest(e-nop,i)*Fs%values(m,e)
				do n=1, numberOfParticles
					Fmi(m,i) = Fmi(m,i) + CoupledCluster_instance%Tstest(e-nop,n)*spinints(m,n,i,e)
					do f=numberOfParticles+1, numberOfContractions
						Fmi(m,i) = Fmi(m,i) + 0.5*taus(e-nop,f-nop,i,n)*spinints(m,n,e,f)
					end do
				end do
			end do
		end do
	end do


	!! Equation 5
! eq ok
	do m=1, numberOfParticles
		do e=numberOfParticles+1, numberOfContractions
			Fme(m,e-nop) = Fs%values(m,e)
			do n=1, numberOfParticles
				do f=numberOfParticles+1, numberOfContractions
				Fme(m,e-nop) = Fme(m,e-nop) + CoupledCluster_instance%Tstest(f-nop,n)*spinints(m,n,e,f)
				end do
			end do
		end do
	end do

	!! Equation 6
!! eq ok
	do m=1, numberOfParticles
		do n=1, numberOfParticles
			do i=1, numberOfParticles
				do j=1, numberOfParticles
					Wmnij(m,n,i,j) = spinints(m,n,i,j)
					do e=numberOfParticles+1, numberOfContractions
						Wmnij(m,n,i,j) = Wmnij(m,n,i,j) + (CoupledCluster_instance%Tstest(e-nop,j)*spinints(m,n,i,e)-CoupledCluster_instance%Tstest(e-nop,i)*spinints(m,n,j,e))
						do f=numberOfParticles+1, numberOfContractions
							Wmnij(m,n,i,j) = Wmnij(m,n,i,j) + 0.25*tau(e-nop,f-nop,i,j)*spinints(m,n,e,f)
						end do
					end do
				end do
			end do
		end do
	end do

	!! Equation 7

! eq ok
	do a=numberOfParticles+1, numberOfContractions
		do b=numberOfParticles+1, numberOfContractions
			do e=numberOfParticles+1, numberOfContractions
				do f=numberOfParticles+1, numberOfContractions
					Wabef(a-nop,b-nop,e-nop,f-nop) = spinints(a,b,e,f)
					do m=1, numberOfParticles
						Wabef(a-nop,b-nop,e-nop,f-nop) = Wabef(a-nop,b-nop,e-nop,f-nop) + (-CoupledCluster_instance%Tstest(b-nop,m)*spinints(a,m,e,f)+CoupledCluster_instance%Tstest(a-nop,m)*spinints(b,m,e,f))
						do n=1, numberOfParticles
							Wabef(a-nop,b-nop,e-nop,f-nop) = Wabef(a-nop,b-nop,e-nop,f-nop) + 0.25*tau(a-nop,b-nop,m,n)*spinints(m,n,e,f)
						end do
					end do
				end do
			end do
		end do
	end do


!! Equation 8
! eq ok

	do m=1, numberOfParticles
		do b=numberOfParticles+1, numberOfContractions
			do e=numberOfParticles+1, numberOfContractions
				do j=1, numberOfParticles
					CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) = spinints(m,b,e,j)
					do f=numberOfParticles+1, numberOfContractions
						CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) = CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) + CoupledCluster_instance%Tstest(f-nop,j)*spinints(m,b,e,f)
					end do
					do n=1, numberOfParticles
!! 26 de enero 2016
                   CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) = CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) + (-CoupledCluster_instance%Tstest(b-nop,n)*spinints(m,n,e,j))
						do f=numberOfParticles+1, numberOfContractions
							CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) = CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) + (-(0.5*CoupledCluster_instance%Tdtest(f-nop,b-nop,j,n)+CoupledCluster_instance%Tstest(f-nop,j)*CoupledCluster_instance%Tstest(b-nop,n))*spinints(m,n,e,f))
						end do
					end do
				end do
			end do
		end do
	end do

!!!! End Intermediates


auxECCSD = 0.0_8

	do i=1, numberOfParticles
		do a=numberOfParticles+1, numberOfContractions
			auxECCSD = auxECCSD + Fs%values(i,a)*CoupledCluster_instance%Tstest(a-nop,i)
			do j=1, numberOfParticles
				do b=numberOfParticles+1, numberOfContractions
					auxECCSD = auxECCSD + (0.25*spinints(i,j,a,b)*CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j)+0.5*spinints(i,j,a,b)*CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j))
		!			write(*,*) auxECCSD
				end do
			end do
		end do
	end do
  ECCSD = auxECCSD
	
  DECC = abs( ECCSD - OLDCC )

   write (*,*) speciesID 
   write (*,*) ECCSD 
!  write(*,*) ECCSD, OLDCC, DECC

!!!! Let's make T1 and T2 (the new ones) 


    !! Equation 1
	 !eq ok

	do a=numberOfParticles+1, numberOfContractions
		do i=1, numberOfParticles
			TsNew(a-nop,i) = Fs%values(i,a)
			do e=numberOfParticles+1, numberOfContractions
				TsNew(a-nop,i) = TsNew(a-nop,i) + CoupledCluster_instance%Tstest(e-nop,i)*Fae(a-nop,e-nop)
			end do
			do m=1, numberOfParticles
				TsNew(a-nop,i) = TsNew(a-nop,i) + (-CoupledCluster_instance%Tstest(a-nop,m)*Fmi(m,i))
				do e=numberOfParticles+1, numberOfContractions
					TsNew(a-nop,i) = TsNew(a-nop,i) + CoupledCluster_instance%Tdtest(a-nop,e-nop,i,m)*Fme(m,e-nop)
					do f=numberOfParticles+1, numberOfContractions
						TsNew(a-nop,i) = TsNew(a-nop,i) + (-0.5*CoupledCluster_instance%Tdtest(e-nop,f-nop,i,m)*spinints(m,a,e,f))
					end do
					do n=1, numberOfParticles
						TsNew(a-nop,i) = TsNew(a-nop,i) + (-0.5*CoupledCluster_instance%Tdtest(a-nop,e-nop,m,n)*spinints(n,m,e,i))
					end do
				end do
			end do
			do n=1,numberOfParticles
				do f=numberOfParticles+1, numberOfContractions
					TsNew(a-nop,i) = TsNew(a-nop,i) + (-CoupledCluster_instance%Tstest(f-nop,n)*spinints(n,a,i,f))
				end do
			end do
			TsNew(a-nop,i) = TsNew(a-nop,i)/Dai(a,i)
			CoupledCluster_instance%Tstest(a-nop,i) = TsNew(a-nop,i)
!			write(*,*) a,i,Ts(a,i),TsNew(a,i)
		end do
	end do

!write (*,*) Ts(:,:)
			

	!! Equation 2

	
	do a=numberOfParticles+1, numberOfContractions
		do b=numberOfParticles+1, numberOfContractions
			do i=1, numberOfParticles
				do j=1, numberOfParticles

					TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + spinints(i,j,a,b) !A
	
   ! 1er ciclo
               do e=numberOfParticles+1, numberOfContractions
	
                  TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (CoupledCluster_instance%Tdtest(a-nop,e-nop,i,j)*Fae(b-nop,e-nop)-CoupledCluster_instance%Tdtest(b-nop,e-nop,i,j)*Fae(a-nop,e-nop)) !B

    					TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (CoupledCluster_instance%Tstest(e-nop,i)*spinints(a,b,e,j)-CoupledCluster_instance%Tstest(e-nop,j)*spinints(a,b,e,i)) !G

						do f=numberOfParticles+1, numberOfContractions
							TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + 0.5*tau(e-nop,f-nop,i,j)*Wabef(a-nop,b-nop,e-nop,f-nop) !D
						end do

						do m=1, numberOfParticles
!! 27 de enero 2016
!                      TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-0.5*Td(a-nop,e-nop,i,j)*Ts(b-nop,m)*Fme(m,e-nop)+0.5*Td(a-nop,e-nop,i,j)*Ts(a-nop,m)*Fme(m,e-nop)) !B'
                      TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-0.5*CoupledCluster_instance%Tdtest(a-nop,e-nop,i,j)*CoupledCluster_instance%Tstest(b-nop,m)*Fme(m,e-nop)+0.5*CoupledCluster_instance%Tdtest(b-nop,e-nop,i,j)*CoupledCluster_instance%Tstest(a-nop,m)*Fme(m,e-nop)) !B'

                  end do
					end do
	! 2do ciclo
               do m=1, numberOfParticles
						TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-CoupledCluster_instance%Tdtest(a-nop,b-nop,i,m)*Fmi(m,j)+CoupledCluster_instance%Tdtest(a-nop,b-nop,j,m)*Fmi(m,i)) !C

						TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-CoupledCluster_instance%Tstest(a-nop,m)*spinints(m,b,i,j)+CoupledCluster_instance%Tstest(b-nop,m)*spinints(m,a,i,j)) !H

						do n=1, numberOfParticles
							TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + 0.5*tau(a-nop,b-nop,m,n)*Wmnij(m,n,i,j) !E
						end do

						do e=numberOfParticles+1, numberOfContractions
!! 27 de enero 2016
!                    TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-0.5*Td(a-nop,b-nop,i,m)*Ts(e-nop,j)*Fme(m,e-nop)+0.5*Td(a-nop,b-nop,i,m)*Ts(e-nop,i)*Fme(m,e-nop)) !C'
                    TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-0.5*CoupledCluster_instance%Tdtest(a-nop,b-nop,i,m)*CoupledCluster_instance%Tstest(e-nop,j)*Fme(m,e-nop)+0.5*CoupledCluster_instance%Tdtest(a-nop,b-nop,j,m)*CoupledCluster_instance%Tstest(e-nop,i)*Fme(m,e-nop)) !C'

                  TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + CoupledCluster_instance%Tdtest(a-nop,e-nop,i,m)*CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,j) - CoupledCluster_instance%Tstest(e-nop,i)*CoupledCluster_instance%Tstest(a-nop,m)*spinints(m,b,e,j) + -CoupledCluster_instance%Tdtest(a-nop,e-nop,j,m)*CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,i) + CoupledCluster_instance%Tstest(e-nop,j)*CoupledCluster_instance%Tstest(a-nop,m)*spinints(m,b,e,i) + -CoupledCluster_instance%Tdtest(b-nop,e-nop,i,m)*CoupledCluster_instance%Wmbejtest(m,a-nop,e-nop,j) - CoupledCluster_instance%Tstest(e-nop,i)*CoupledCluster_instance%Tstest(b-nop,m)*spinints(m,a,e,j) + CoupledCluster_instance%Tdtest(b-nop,e-nop,j,m)*CoupledCluster_instance%Wmbejtest(m,a-nop,e-nop,i) - CoupledCluster_instance%Tstest(e-nop,j)*CoupledCluster_instance%Tstest(b-nop,m)*spinints(m,a,e,i) !F

                  end do
					end do


!!* Eq 13 Stanton

!!	22 enero 2016			  
!!!! Make denominator array Dabij     |   Dabij(a-nop,b-nop,i,j) = Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)
               
               TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b))
               CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j)
                                        taus(a-nop,b-nop,i,j) = CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) + 0.5*(CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j))
                                        tau(a-nop,b-nop,i,j) = CoupledCluster_instance%Tdtest(a-nop,b-nop,i,j) + CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j)
				end do
			end do
		end do
	end do


!!!! End New T1 and T2

	CoupledCluster_instance%coupledClusterSDCorrection%values(speciesID) = ECCSD ! CCSD Correlation Energy

 end do !! Main Loop



	end if !! If numberOfParticles > 1	


 end do !! Species


	call Vector_show (CoupledCluster_instance%coupledClusterSDCorrection)
	CoupledCluster_instance%ccsdCorrection = sum( CoupledCluster_instance%coupledClusterSDCorrection%values ) !CCSD Corr. Energy
 
! 28 de enero 2016

!	if (allocated(Ts)) deallocate (Ts)
        if (allocated(Td)) deallocate (Td)
        if (allocated(taus)) deallocate (taus)
        if (allocated(tau)) deallocate (tau)
	if (allocated(Fae)) deallocate (Fae)
	if (allocated(Fmi)) deallocate (Fmi)
	if (allocated(Fme)) deallocate (Fme)
	if (allocated(Wmnij)) deallocate (Wmnij)
	if (allocated(Wabef)) deallocate (Wabef)
	if (allocated(Wmbej)) deallocate (Wmbej)
	if (allocated(TsNew)) deallocate (TsNew)
	if (allocated(TdNew)) deallocate (TdNew)
	call Vector_destructor (ff)
	call Matrix_destructor (Fs)
	call Matrix_destructor (auxMatrix)

!!*********************************************************************************************************************************
!! End CCSD calculation same specie!
!!*********************************************************************************************************************************

   close(wfnUnit)

 end subroutine CoupledCluster_iterateIntermediates_SameSpecies


!!*********************************************************************************************************************************
!! Begin CCSD calculation different species!
!!*********************************************************************************************************************************



 subroutine CoupledCluster_iterateIntermediates_DiffSpecies()
   implicit none

   integer :: numberOfSpecies
   integer :: a,b,c,d
   integer :: aa,bb,cc,dd
   integer :: i,j,k,l,p,q,r,s,h,t
   integer :: ii,jj,kk,ll,pp,qq,rr,ss,hh,tt
   integer(8) :: x,y,z,auxIndex,auxIndex2
   integer :: e,f,m,mmm,n,iii,jjj,aaa,bbb
   integer :: ee,mm,nn,fff,oo
   integer :: speciesID
   integer :: otherSpeciesID
   character(10) :: nameOfSpecie
   character(10) :: nameOfOtherSpecie
   integer :: electronsID
   integer :: numberOfParticles
   integer :: numberOfOtherSpecieParticles
   integer :: ocupationNumber
   integer :: ocupationNumberOfOtherSpecie
   integer :: numberOfContractions
   integer :: numberOfContractionsOfOtherSpecie
   integer(8) :: numberOfOrbitals
   integer(8) :: numberOfSpatialOrbitals
   integer :: numberOfOtherSpecieOrbitals
   integer :: numberOfOtherSpecieSpatialOrbitals
   integer, allocatable :: spin(:)
   integer, allocatable :: spatialOrbital(:)
   type(Vector) :: eigenValues, ff
   type(Vector) :: eigenValuesOfOtherSpecie, otherff
   type(Vector) :: coupledClusterValue
   type(Matrix) :: auxMatrix!   type(TransformIntegrals) :: repulsionTransformer
   real(8) :: lambda
   real(8) :: lambdaOfOtherSpecie
   real(8) :: kappa !positive or negative exchange
   real(8) :: charge
   real(8) :: TwoParticlesEnergy, ECCSD, DECC, OLDCC, auxECCSD, ECCSD1, ECCSD2
   real(8) :: otherSpecieCharge
   real(8) :: independentEnergyCorrection
   real(8) :: mp2CouplingCorrection
   real(8) :: auxVal,auxVal1,auxVal2,value1,value2
   real(8) :: auxVal_A
   real(8) :: auxVal_B
   type(Matrix) :: eigenVec, eigenVec1, Fs
   type(Matrix) :: eigenVecOtherSpecie, eigenVecOtherSpecie1, otherFs
!   real(8), allocatable :: ff(:)
   real(8), allocatable :: Ts(:,:),Fae(:,:),Fmi(:,:),Fme(:,:),Dai(:,:),TsNew(:,:)
   real(8), allocatable :: auxTs(:,:),auxFae(:,:),auxFmi(:,:),auxFme(:,:),auxDai(:,:),auxTsNew(:,:)
   real(8), allocatable :: otherTs(:,:),otherFae(:,:),otherFmi(:,:),otherFme(:,:),otherDai(:,:),otherTsNew(:,:)
   real(8), allocatable :: Td(:,:,:,:), spinints(:,:,:,:), Dabij(:,:,:,:), taus(:,:,:,:), tau(:,:,:,:), Wmnij(:,:,:,:), Wabef(:,:,:,:), Wmbej(:,:,:,:), TdNew(:,:,:,:)
   real(8), allocatable :: otherTd(:,:,:,:), otherspinints(:,:,:,:), otherDabij(:,:,:,:), othertaus(:,:,:,:), auxspinints(:,:,:,:), auxTd(:,:,:,:)
   real(8), allocatable :: othertau(:,:,:,:), otherWmnij(:,:,:,:), otherWabef(:,:,:,:), otherWmbej(:,:,:,:), otherTdNew(:,:,:,:)
   real(8), allocatable :: auxWmnij(:,:,:,:), auxWabef(:,:,:,:), auxWmbej(:,:,:,:)
   real(8), allocatable :: auxTdNew(:,:,:,:), auxtau(:,:,:,:), auxtaus(:,:,:,:), auxDabij(:,:,:,:)
   type(Matrix), allocatable :: auxMatrix1(:,:)
   character(50) :: wfnFile
   character(50) :: arguments(2)
   integer :: wfnUnit
   !! 21 enero 2016
   integer :: noc
   integer :: nop
   integer :: nocs
   integer :: nops
!! 27 de enero 2016
   integer :: kro

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! Load results...
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFock_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


!!**      !!!! Initial guesses T1 and T2
!!**      
              if (allocated(Ts)) deallocate (Ts)
                allocate(Ts(noc-nop,nop))
              Ts(:,:) = 0.0_8
!!**              !! Td solo corre sobre virtual,virtual,ocupado,ocupado
      	if (allocated(TsNew)) deallocate (TsNew)
!!**      !!	allocate(TsNew(numberOfContractions,numberOfContractions))
         allocate(TsNew(noc-nop,nop))
      	TsNew=0.0_8


        !! Td solo corre sobre virtual,virtual,ocupado,ocupado
        if (allocated(Td)) deallocate (Td)
 !! 22 de enero 2016
 !!      allocate(Td(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(Td(noc-nop,noc-nop,nop,nop))
        Td(:,:,:,:) = 0.0_8

        if (allocated(taus)) deallocate (taus)
!! 22 enero 2016 
!!      allocate(taus(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(taus(noc-nop,noc-nop,nop,nop))
        taus(:,:,:,:) = 0.0_8

        if (allocated(tau)) deallocate (tau)
!! 22 enero 2016
!!      allocate(tau(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
        allocate(tau(noc-nop,noc-nop,nop,nop))
        tau(:,:,:,:) = 0.0_8

!  29 de enero 2016

	if (allocated(Wmbej)) deallocate (Wmbej)
        allocate(Wmbej(nop,noc-nop,noc-nop,nop))
	Wmbej=0.0_8


!CoupledCluster_instance%Tstest(:,:) = 5*5


 !  print *, "55555mi variable copartida ",  CoupledCluster_instance%Tstest(1,1) 
!01 de febrero 2016

!        do a=numberOfParticles+1, numberOfContractions
!                do b=numberOfParticles+1, numberOfContractions
!                        do i=1, numberOfParticles
!                                do j=1, numberOfParticles
 !!                                       Td(a,b,i,j) = Td(a,b,i,j) + (spinints(i,j,a,b)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)))
 !!                                       taus(a,b,i,j) = Td(a,b,i,j) + 0.5*(Ts(a,i)*Ts(b,j) - Ts(b,i)*Ts(a,j))
 !!                                       tau(a,b,i,j) = Td(a,b,i,j) + Ts(a,i)*Ts(b,j) - Ts(b,i)*Ts(a,j)

!                                        Td(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + (spinints(i,j,a,b)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)))
                 !! 22 de enero 2016     
!                                        taus(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + 0.5*(CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j))
!                                        tau(a-nop,b-nop,i,j) = Td(a-nop,b-nop,i,j) + CoupledCluster_instance%Tstest(a-nop,i)*CoupledCluster_instance%Tstest(b-nop,j) - CoupledCluster_instance%Tstest(b-nop,i)*CoupledCluster_instance%Tstest(a-nop,j)
!			                  write(*,*) Td(a,b,i,j), taus(a,b,i,j), tau(a,b,i,j)
!                                end do
!                        end do
!                end do
!        end do

!!!! End Initial Guesses




!!**      !!!! Let's make T1 and T2 (the new ones) 
!!**      
!!**      
!!**          !! Equation 1 Ts intra
!!**      	 !eq ok
!!**     
!!**     Serà usado nuevamente al combinar las excitaciones entre diferentes especies
 
! 01 febrero 2016

!!**         do a=numberOfParticles+1, numberOfContractions
!!**      		do i=1, numberOfParticles
!!**      			TsNew(a-nop,i) = Fs%values(i,a)
!!**      			do e=numberOfParticles+1, numberOfContractions
!!**      				TsNew(a-nop,i) = TsNew(a-nop,i) + CoupledCluster_instance%Tstest(e-nop,i)*Fae(a-nop,e-nop)
!!**      			end do
!!**      			do m=1, numberOfParticles
!!**      				TsNew(a-nop,i) = TsNew(a-nop,i) + (-CoupledCluster_instance%Tstest(a-nop,m)*Fmi(m,i))
!!**      				do e=numberOfParticles+1, numberOfContractions
!!**      					TsNew(a-nop,i) = TsNew(a-nop,i) + Td(a-nop,e-nop,i,m)*Fme(m,e-nop)
!!**      					do f=numberOfParticles+1, numberOfContractions
!!**      					TsNew(a-nop,i) = TsNew(a-nop,i) + (-0.5*Td(e-nop,f-nop,i,m)*spinints(m,a,e,f))
!!**      					end do
!!**      					do n=1, numberOfParticles
!!**      						TsNew(a-nop,i) = TsNew(a-nop,i) + (-0.5*Td(a-nop,e-nop,m,n)*spinints(n,m,e,i))
!!**      					end do
!!**      				end do
!!**      			end do
!!**      			do n=1,numberOfParticles
!!**      				do f=numberOfParticles+1, numberOfContractions
!!**      					TsNew(a-nop,i) = TsNew(a-nop,i) + (-CoupledCluster_instance%Tstest(f-nop,n)*spinints(n,a,i,f))
!!**      				end do
!!**      			end do
!!**      			TsNew(a-nop,i) = TsNew(a-nop,i)/Dai(a,i)
!!**      			Ts(a-nop,i) = TsNew(a-nop,i)
!!**      !			write(*,*) a,i,Ts(a,i),TsNew(a,i)
!!**      		end do
!!**      	end do
      


      !write (*,*) Ts(:,:)


if ( numberOfSpecies > 1 ) then

mp2CouplingCorrection = 0.0_8
mmm = 0

print *, " "
print *, " "
print *, " |------------------------------------------------|"
print *, " |Calculating CCSD corr. for different species....|"
print *, " |------------------------------------------------|"



!   CoupledCluster_instance%Tstest = 136 * 100

!   print *, "mi variable copartida *100 ",  CoupledCluster_instance%Tstest(1,1) 

f = numberOfSpecies * ( numberOfSpecies-1 ) / 2

call Vector_constructor( CoupledCluster_instance%mp2Coupling, f)
call Vector_constructor( coupledClusterValue, numberOfSpecies)

  do i=1, numberOfSpecies

     numberOfContractions = MolecularSystem_getTotalNumberOfContractions(i)
     arguments(2) = trim(MolecularSystem_getNameOfSpecie(i))

     arguments(1) = "COEFFICIENTS"
     eigenVec = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
             columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

     arguments(1) = "ORBITALS"
     call Vector_getFromFile( elementsNum = numberOfContractions, &
          unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
          output = eigenValues )

     nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
     speciesID =MolecularSystem_getSpecieID( nameOfSpecie=trim(nameOfSpecie) )
     ocupationNumber = MolecularSystem_getOcupationNumber( i )
     numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
     lambda = MolecularSystem_getlambda( i )

     do j = i + 1, numberOfSpecies
        mmm = mmm + 1


        numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )

        arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))

        arguments(1) = "COEFFICIENTS"
        eigenVecOtherSpecie = &
             Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractionsOfOtherSpecie,4), &
             columns= int(numberOfContractionsOfOtherSpecie,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "ORBITALS"
        call Vector_getFromFile( elementsNum = numberOfContractionsofOtherSpecie, &
             unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
             output = eigenValuesOfOtherSpecie )


        numberOfParticles = MolecularSystem_getNumberOfParticles(i)
        numberOfOtherSpecieParticles = MolecularSystem_getNumberOfParticles(j)
        nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
        otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
        ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )

        lambdaOfOtherSpecie = MolecularSystem_instance%species(j)%lambda

	!! Read transformed integrals from file
        call ReadTransformedIntegrals_readTwoSpecies( speciesID, otherSpeciesID, auxMatrix)
	 auxMatrix%values = auxMatrix%values * MolecularSystem_getCharge( speciesID ) * MolecularSystem_getCharge( otherSpeciesID )


	do a=1, ocupationNumber
           do p=1,ocupationNumberOfOtherSpecie
              do r=ocupationNumber+1, numberOfContractions
                 do t=ocupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie

                    auxIndex = IndexMap_tensorR4ToVector(a,r,p,t, numberOfContractions, &
                         numberOfContractionsOfOtherSpecie )

                    mp2CouplingCorrection = mp2CouplingCorrection + ( ( auxMatrix%values(auxIndex,1) )**2.0_8 ) / (eigenValues%values(a) + eigenValuesOfOtherSpecie%values(p)-eigenValues%values(r)-eigenValuesOfOtherSpecie%values(t) )
                 end do
              end do
           end do
        end do

        CoupledCluster_instance%mp2Coupling%values(mmm)= ( lambda * lambdaOfOtherSpecie * mp2CouplingCorrection )

!end do !!! Number Of Species
!end do !!! 

!!! From Scratch

	call Matrix_diagonalConstructor (eigenVec1, eigenValues) ! put MO energies in diagonal array (eigenVec) 
        call Vector_constructor (ff,numberOfContractions*2,0.0_8)

        a=0
        b=0
        do x=1, numberOfContractions
                a=a+1
                do h=1, 2
                        b=b+1
                        ff%values(b) = eigenValues%values(a)
                end do
        end do

        call Matrix_diagonalConstructor (Fs, ff)

	call Matrix_diagonalConstructor (eigenVecOtherSpecie1, eigenValuesOfOtherSpecie) ! put MO energies in diagonal array (eigenVec) 
        call Vector_constructor (otherff,numberOfContractionsOfOtherSpecie*2,0.0_8)

        a=0
        b=0
        do x=1, numberOfContractionsOfOtherSpecie
                a=a+1
                do h=1, 2
                        b=b+1
                        otherff%values(b) = eigenValuesOfOtherSpecie%values(a)
                end do
        end do

        call Matrix_diagonalConstructor (otherFs, otherff)

	numberOfContractions=numberOfContractions*2

        noc=numberOfContractions
        nop=numberOfParticles
	
   numberOfContractionsOfOtherSpecie=numberOfContractionsOfOtherSpecie*2

!! 25 de enero 2016

        nocs=numberOfContractionsOfOtherSpecie 
        nops=numberOfOtherSpecieParticles
!
!

!01 de febrero 2016
!!**        if (allocated(spinints)) deallocate (spinints)
!!**        allocate(spinints(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!!**        spinints(:,:,:,:) = 0.0_8
!!**
!!**!! pasa de 4 indices a 1 (pairing function)
!!**        do p=1, numberOfContractions
!!**                do q=1, numberOfContractions
!!**                        do r=1, numberOfContractions
!!**                                do s=1, numberOfContractions
!!**                                        value1 = IndexMap_tensorR4ToVector((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,numberOfContractions/2) !! integrales de Coulomb
!!**                                        auxVal_A= auxMatrix%values(value1, 1)
!!**                                        value2 = IndexMap_tensorR4ToVector((p+1)/2,(s+1)/2,(q+1)/2,(r+1)/2,numberOfContractions/2) !! integrales de intercambio
!!**                                        auxVal_B= auxMatrix%values(value2, 1)
!!**                                        auxVal1 = auxVal_A * (mod(p,2) == mod(r,2)) * (mod(q,2) == mod(s,2))
!!**                                        auxVal2 = auxVal_B * (mod(p,2) == mod(s,2)) * (mod(q,2) == mod(r,2))
!!**                                        spinints(p,q,r,s) = auxVal1 - auxVal2 !! p+1 o p-1? Revisar !! ecuacion 1
!!**!					write (*,*) spinints(p,q,r,s)
!!**                                end do
!!**                        end do
!!**                end do
!!**        end do






        if (allocated(auxspinints)) deallocate (auxspinints)
        allocate(auxspinints(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))


	!!!! TEST for combined spinints
!	call ReadTransformedIntegrals_readTwoSpecies( speciesID, otherSpeciesID, auxspinints)
        do p=1, numberOfContractions
                do q=1, numberOfContractionsOfOtherSpecie
                        do r=1, numberOfContractions
                                do s=1, numberOfContractionsOfOtherSpecie
                                        value1 = IndexMap_tensorR4ToVector((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,numberOfContractions/2, numberOfContractionsOfOtherSpecie/2)
                                        auxVal_A= auxMatrix%values(value1, 1)
                                        auxVal1 = auxVal_A * (mod(p,2) == mod(r,2)) * (mod(q,2) == mod(s,2))
                                        auxspinints(p,q,r,s) = auxVal1
                                end do
                        end do
                end do
        end do

!!!! Initial guesses T1 and T2

        if (allocated(Ts)) deallocate (Ts)
!! 26 de enero de 2016
!!      allocate(Ts(numberOfContractions,numberOfContractions))
        allocate(Ts(noc-nop,nop))
        Ts(:,:) = 0.0_8

        if (allocated(otherTs)) deallocate (otherTs)
!! 26 enerdo de 2016
!!      allocate(otherTs(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
        allocate(otherTs(nocs-nops,nops))
        otherTs(:,:) = 0.0_8
        
!	if (allocated(auxTs)) deallocate (auxTs)
!        allocate(auxTs(numberOfContractions,numberOfContractionsOfOtherSpecie))
!        auxTs(:,:) = 0.0_8

!        if (allocated(otherTd)) deallocate (otherTd)
!        allocate(otherTd(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
!        otherTd(:,:,:,:) = 0.0_8

        if (allocated(auxTd)) deallocate (auxTd)
!! 25 de enero de 2016
!!      allocate(auxTd(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxTd(noc-nop,nocs-nops,nop,nops))
        auxTd(:,:,:,:) = 0.0_8

!        if (allocated(othertaus)) deallocate (othertaus)
!        allocate(othertaus(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
!        othertaus(:,:,:,:) = 0.0_8

!        if (allocated(othertau)) deallocate (othertau)
!        allocate(othertau(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
!        othertau(:,:,:,:) = 0.0_8
        
        if (allocated(auxtaus)) deallocate (auxtaus)
!! 25 de enero 2016
!!      allocate(auxtaus(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxtaus(noc-nop,nocs-nops,nop,nops))
        auxtaus(:,:,:,:) = 0.0_8

        if (allocated(auxtau)) deallocate (auxtau)
        allocate(auxtau(noc-nop,nocs-nops,nop,nops))
        auxtau(:,:,:,:) = 0.0_8

!        do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                        do ii=1, numberOfOtherSpecieParticles
!                                do jj=1, numberOfOtherSpecieParticles
!                                        otherTd(a,b,ii,jj) = otherTd(a,b,ii,jj) + (otherspinints(ii,jj,a,b)/(otherFs%values(ii,ii)+otherFs%values(jj,jj)-otherFs%values(a,a)-otherFs%values(b,b)))
!                                        othertaus(a,b,ii,jj) = otherTd(a,b,ii,jj) + 0.5*(otherTs(a,ii)*otherTs(b,jj) - otherTs(b,ii)*otherTs(a,jj))
!                                        othertau(a,b,ii,jj) = otherTd(a,b,ii,jj) + otherTs(a,ii)*otherTs(b,jj) - otherTs(b,ii)*otherTs(a,jj)
!			                  write(*,*) a,b,ii,jj,otherTd(a,b,i,j)
!                                end do
!                        end do
!                end do
!        end do






	!!! TEST for combined Td auxTd

        do aa=numberOfParticles+1, numberOfContractions
                do aaa=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
                        do ii=1, numberOfParticles
                                do iii=1, numberOfOtherSpecieParticles
                                        auxTd(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) + (auxspinints(ii,iii,aa,aaa)/(Fs%values(ii,ii)+otherFs%values(iii,iii)-Fs%values(aa,aa)-otherFs%values(aaa,aaa)))
                                        auxtaus(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) + 0.5*(CoupledCluster_instance%Tstest(aa-nop,ii)*otherTs(aaa-nops,iii)) ! Eliminated crossed t1 terms
                                        auxtau(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) + CoupledCluster_instance%Tstest(aa-nop,ii)*otherTs(aaa-nops,iii) ! same here
!			                  write(*,*) auxTd(aa,aaa,ii,iii)
                                end do
                        end do
                end do
        end do


!!!! End Initial Guesses

!!!! Make denominator arrays Dai, Dabij
 
         !! Equation 12 from Stanton
 
!         if (allocated(Dai)) deallocate (Dai)
!         allocate(Dai(numberOfContractions,numberOfContractions))
!         Dai(:,:) = 0.0_8
! 
!         do a=numberOfParticles+1, numberOfContractions
!                 do ii=1, numberOfParticles
!                         Dai(a,ii) = Fs%values(ii,ii) - Fs%values(a,a)
! !                       write(*,*) a,i,Dai(a,i)
!                 end do
!         end do
!
!         if (allocated(otherDai)) deallocate (otherDai)
!         allocate(otherDai(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
!         otherDai(:,:) = 0.0_8
! 
!         do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                 do ii=1, numberOfOtherSpecieParticles
!                         otherDai(a,ii) = otherFs%values(ii,ii) - otherFs%values(a,a)
! !                       write(*,*) a,i,Dai(a,i)
!                 end do
!         end do
! 
!         if (allocated(auxDai)) deallocate (auxDai)
!         allocate(auxDai(numberOfContractions,numberOfContractionsOfOtherSpecie))
!         auxDai(:,:) = 0.0_8
! 
!         do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                 do ii=1, numberOfOtherSpecieParticles
!                         otherDai(a,ii) = otherFs%values(ii,ii) - otherFs%values(a,a)
! !                       write(*,*) a,i,Dai(a,i)
!                 end do
!         end do

	!! Equation 13 from Stanton
 
!        if (allocated(Dabij)) deallocate (Dabij)
!        allocate(Dabij(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!        Dabij(:,:,:,:) = 0.0_8
!
!
!        do a=numberOfParticles+1, numberOfContractions
!                do b=numberOfParticles+1, numberOfContractions
!                        do ii=1, numberOfParticles
!                                do jj=1, numberOfParticles
!                                        Dabij(a,b,ii,jj) = Fs%values(ii,ii)+Fs%values(jj,jj)-Fs%values(a,a)-Fs%values(b,b)
!!					write(*,*) a,b,i,j,Dabij(a,b,i,j)
!                                end do
!                        end do
!                end do
!        end do
!
! 
!         if (allocated(otherDabij)) deallocate (otherDabij)
!         allocate(otherDabij(numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie,numberOfContractionsOfOtherSpecie))
!         otherDabij(:,:,:,:) = 0.0_8
! 
!         do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                 do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                         do ii=1, numberOfOtherSpecieParticles
!                                 do jj=1, numberOfOtherSpecieParticles
!                                         otherDabij(a,b,ii,jj) = otherFs%values(ii,ii)+otherFs%values(jj,jj)-otherFs%values(a,a)-otherFs%values(b,b)
!!                                         write(*,*) a,b,ii,jj,Dabij(a,b,i,j)
!                                 end do
!                         end do
!                 end do
!         end do
!
!         if (allocated(auxDabij)) deallocate (auxDabij)
!! 25 de enero 2016
!!       allocate(auxDabij(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
!!         allocate(auxDabij(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
!!         auxDabij(:,:,:,:) = 0.0_8

!!         do aa=numberOfParticles+1, numberOfContractions
!!                 do aaa=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!!                         do ii=1, numberOfParticles
!!                                 do iii=1, numberOfOtherSpecieParticles
!!                                         auxDabij(aa,aaa,ii,iii) = Fs%values(ii,ii)+otherFs%values(iii,iii)-Fs%values(aa,aa)-otherFs%values(aaa,aaa)
!                                         write(*,*) a,b,ii,jj,Dabij(a,b,i,j)
!!                                 end do
!!                         end do
!!                 end do
!!         end do

 !!!! End Denominators


!!! Loop

  auxECCSD = 0.0_8
  DECC = 1.0_8
  h = 0
  do while (DECC >= 1.0D-8)
!  do while (h<20)
!	h=h+1
        OLDCC = auxECCSD

!	write(*,*) OLDCC


!!!!! Intermediates for different particles

!	if (allocated(Fae)) deallocate (Fae)
!        allocate(Fae(numberOfContractions,numberOfContractions))
!        Fae=0.0_8
!        if (allocated(Fmi)) deallocate (Fmi)
!        allocate(Fmi(numberOfContractions,numberOfContractions))
!        Fmi=0.0_8
!        if (allocated(Fme)) deallocate (Fme)
!        allocate(Fme(numberOfContractions,numberOfContractions))
!        Fme=0.0_8
!
!        if (allocated(Wmnij)) deallocate (Wmnij)
!        allocate(Wmnij(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!        Wmnij=0.0_8
!        if (allocated(Wabef)) deallocate (Wabef)
!        allocate(Wabef(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!        Wabef=0.0_8
!        if (allocated(Wmbej)) deallocate (Wmbej)
!        allocate(Wmbej(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!        Wmbej=0.0_8
!
!        if (allocated(TsNew)) deallocate (TsNew)
!        allocate(TsNew(numberOfContractions,numberOfContractions))
!        TsNew=0.0_8
!        if (allocated(TdNew)) deallocate (TdNew)
!        allocate(TdNew(numberOfContractions,numberOfContractions,numberOfContractions,numberOfContractions))
!        TdNew=0.0_8
!
!!
        if (allocated(auxFae)) deallocate (auxFae)
        allocate(auxFae(noc,nocs))
        auxFae=0.0_8
!        if (allocated(auxFmi)) deallocate (auxFmi)
!        allocate(auxFmi(numberOfContractions,numberOfContractionsOfOtherSpecie))
!        auxFmi=0.0_8
        if (allocated(otherFme)) deallocate (otherFme)
        allocate(otherFme(nops,nocs-nops))
        otherFme=0.0_8

        if (allocated(auxWmnij)) deallocate (auxWmnij)
!! 25 de enero 2016
!!        allocate(auxWmnij(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxWmnij(nop,nocs,nop,nocs))
        auxWmnij=0.0_8
        if (allocated(auxWabef)) deallocate (auxWabef)
!! 25 de enero 2016
!!      allocate(auxWabef(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxWabef(noc-nop,nocs-nops,noc-nop,nocs-nops))
        auxWabef=0.0_8
        if (allocated(auxWmbej)) deallocate (auxWmbej)
!! 25 de enero 2016
!!      allocate(auxWmbej(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxWmbej(nop,nocs-nops,noc-nop,nops))
        auxWmbej=0.0_8

!        if (allocated(auxTsNew)) deallocate (auxTsNew)
!        allocate(auxTsNew(numberOfContractions,numberOfContractionsOfOtherSpecie))
!        auxTsNew=0.0_8
        if (allocated(auxTdNew)) deallocate (auxTdNew)
!! 25 de enero 2016
!!      allocate(auxTdNew(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))
        allocate(auxTdNew(noc-nop,nocs-nops,nop,nops))
        auxTdNew=0.0_8


	!! Equation 3

! 02 de febrero 2016 

	do a=nops+1, nocs
		do e=nops+1, nocs
 
!Termino de misma especie ya está descrito en intraespecie 
!!         if (B==E) then 
!!            kro = 1 
!!         else 
!!            kro = 0 
!!         end if
!!         
!!         auxFae(B-nops,E-nops) = (1 - kro)*otherFs%values(A,E)
		
         do m=1, numberOfParticles
!Termino de misma especie ya está descrito en intraespecie
!				auxFae(B-nops,E-nops) = auxFae(B-nops,E-nops) + (-0.5*otherFs%values(M,E)*CoupledCluster_instance%Tstest(B-nops,M))
				do f=numberOfParticles+1, numberOfContractions
					auxFae(a-nops,e-nops) = auxFae(a-nops,e-nops) + CoupledCluster_instance%Tstest(f-nop,m)*auxspinints(m,a,f,e)
            end do
         end do
         do m=1, nops
            do n=1, numberOfParticles
               do f=numberOfParticles+1, numberOfContractions
                 ! Tau puede presentarse de dos formas:  t_(Mn)´(Bf) = t_(nM)´(fB) 
						auxFae(a-nops,e-nops) = auxFae(a-nops,e-nops) + (-0.5*auxtaus(f-nop,a-nops,n,m)*auxspinints(n,m,f,e))
!						write(*,*) a,e,Fae(a,e)
					end do
				end do
			end do

		end do
	end do


  print *, "Variable F(BE)",  auxFae(1,1) 



!	write(*,*) auxFae
!
!        do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!             do e=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                     auxFae(a,e) = auxFae(a,e) + (1 + (a==e))*otherFs%values(a,e)
!                     do m=1, numberOfOtherSpecieParticles
!                             auxFae(a,e) = auxFae(a,e) + (-0.5*otherFs%values(m,e)*otherTs(a,m))
!                             do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                                     auxFae(a,e) = auxFae(a,e) + otherTs(f,m)*otherspinints(m,a,f,e)
!                                     do n=1, numberOfOtherSpecieParticles
!                                             auxFae(a,e) = auxFae(a,e) + (-0.5*othertaus(a,f,m,n)*otherspinints(m,n,e,f))
!!                                             write(*,*) auxFae(a,e)
!                                     end do
!                             end do
!                     end do
!             end do
!        end do

	!! Equation 4


	
!
!	do m=1, numberOfParticles
!		do ii=1, numberOfParticles
!			do mm=1, numberOfOtherSpecieParticles
!				do iii=1, numberOfOtherSpecieParticles
!					auxFmi(m,ii) = (1 + (m==ii))*Fs%values(m,ii)
!					auxFmi(m,ii) = auxFmi(m,ii) + (1 + (mm==iii))*otherFs%values(mm,iii)
!						do e=numberOfParticles+1, numberOfContractions
!						do ee=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!							auxFmi(m,ii) = auxFmi(m,ii) + 0.5*Ts(e,ii)*Fs%values(m,e)
!							auxFmi(m,ii) = auxFmi(m,ii) + 0.5*otherTs(ee,iii)*otherFs%values(mm,ee)
!								do n=1, numberOfParticles
!								do nn=1, numberOfOtherSpecieParticles
!									auxFmi(m,ii) = auxFmi(m,ii) + Ts(e,n)*spinints(m,n,ii,e)
!									auxFmi(m,ii) = auxFmi(m,ii) + otherTs(ee,nn)*otherspinints(mm,nn,iii,ee)
!										do f=numberOfParticles+1, numberOfContractions
!										do fff=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!											auxFmi(m,ii) = auxFmi(m,ii) + 0.5*taus(e,f,ii,n)*spinints(m,n,e,f)
!											auxFmi(m,ii) = auxFmi(m,ii) + 0.5*othertaus(ee,fff,iii,nn)*otherspinints(mm,nn,ee,fff)
!										end do
!										end do
!								end do
!								end do
!						end do
!						end do					
!				end do
!			end do
!		end do
!	end do
	!write(*,*) auxFmi
!
!	do m=1, numberOfOtherSpecieParticles
!             do ii=1, numberOfOtherSpecieParticles
!                     otherFmi(m,ii) = (1 + (m==ii))*otherFs%values(m,ii)
!                     do e=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                             otherFmi(m,ii) = otherFmi(m,ii) + 0.5*otherTs(e,ii)*otherFs%values(m,e)
!                             do n=1, numberOfOtherSpecieParticles
!                                     otherFmi(m,ii) = otherFmi(m,ii) + otherTs(e,n)*otherspinints(m,n,ii,e)
!                                     do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                                             otherFmi(m,ii) = otherFmi(m,ii) + 0.5*taus(e,f,ii,n)*otherspinints(m,n,e,f)
!                                     end do  
!                             end do  
!                     end do  
!             end do  
!	end do 
!
	!! Equation 5


	!! Equation 5
! eq ok
	do m=1, nops
		do e=nops+1, nocs
!   Esta parte de la diagonalización ya está descrita en intraespecie  
!			otherFme(m,e-nops) = otherFs%values(m,e)
			do n=1, numberOfParticles
				do f=numberOfParticles+1, numberOfContractions
				otherFme(m,e-nops) = otherFme(m,e-nops) + CoupledCluster_instance%Tstest(f-nop,n)*auxspinints(n,m,f,e)
				end do
			end do
		end do
	end do



  print *, "Variable F(ME)",  otherFme(1,1) 


!	do m=1, numberOfParticles
!		do e=numberOfParticles+1, numberOfContractions
!			Fme(m,e) = Fs%values(m,e)
!			do n=1, numberOfParticles
!				do f=numberOfParticles+1, numberOfContractions
!				Fme(m,e) = Fme(m,e) + Ts(f,n)*spinints(m,n,e,f)
!				end do
!			end do
!		end do
!	end do
!
!	do m=1, numberOfOtherSpecieParticles
!             do e=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                     otherFme(m,e) = otherFs%values(m,e)
!                     do n=1, numberOfOtherSpecieParticles
!                             do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                             otherFme(m,e) = otherFme(m,e) + otherTs(f,n)*otherspinints(m,n,e,f)
!                             end do          
!                     end do                  
!             end do                                  
!	end do     
!
	!! Equation 6

! Antes

!*   do m=1, numberOfParticles
!*		do n=1, numberOfOtherSpecieParticles
!*			do ii=1, numberOfParticles
!*				do jj=1, numberOfOtherSpecieParticles
!*					auxWmnij(m,n,ii,jj) = auxspinints(m,n,ii,jj)
!*					do e=numberOfParticles+1, numberOfContractions
!*      !        		auxWmnij(m,n,ii,jj) = auxWmnij(m,n,ii,jj) + (Ts(e,jj)*auxspinints(m,n,ii,e)-Ts(e,ii)*spinints(m,n,jj,e))
!*						do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!*							auxWmnij(m,n,ii,jj) = auxWmnij(m,n,ii,jj) + 0.25*auxtau(e-nop,f-nops,ii,jj)*auxspinints(m,n,e,f)
!*						end do
!*					end do
!*				end do
!*			end do
!*		end do
!*	end do




!! Eq 6
!!****test Ahora


	do m=1, numberOfParticles
		do n=1, numberOfOtherSpecieParticles
			do ii=1, numberOfParticles
				do jj=1, numberOfOtherSpecieParticles
					auxWmnij(m,n,ii,jj) = auxspinints(m,n,ii,jj)
					do e=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
             		auxWmnij(m,n,ii,jj) = auxWmnij(m,n,ii,jj) + otherTs(e-nops,jj)*auxspinints(m,n,ii,e) !! termino anulado -otherTs(e-nops,ii)*auxspinints(m,n,jj,e))
					end do
               do e=numberOfParticles+1, numberOfContractions
                  do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
   						auxWmnij(m,n,ii,jj) = auxWmnij(m,n,ii,jj) + 0.25*auxtau(e-nop,f-nops,ii,jj)*auxspinints(m,n,e,f)
						end do
					end do
				end do
			end do
		end do
	end do


!!***test





!	!! Equation 7
!
!	do a=numberOfParticles+1, numberOfContractions
!		do b=numberOfParticles+1, numberOfContractions
!			do e=numberOfParticles+1, numberOfContractions
!				do f=numberOfParticles+1, numberOfContractions
!					Wabef(a,b,e,f) = spinints(a,b,e,f)
!					do m=1, numberOfParticles
!						Wabef(a,b,e,f) = Wabef(a,b,e,f) + (-Ts(b,m)*spinints(a,m,e,f)+Ts(a,m)*spinints(b,m,e,f))
!						do n=1, numberOfParticles
!							Wabef(a,b,e,f) = Wabef(a,b,e,f) + 0.25*tau(a,b,m,n)*spinints(m,n,e,f)
!						end do
!					end do
!				end do
!			end do
!		end do
!	end do

	do a=numberOfParticles+1, numberOfContractions
		do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
			do e=numberOfParticles+1, numberOfContractions
				do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
					auxWabef(a-nop,b-nops,e-nop,f-nops) = auxspinints(a,b,e,f)
	        		do m=1, numberOfOtherSpecieParticles
             		auxWabef(a-nop,b-nops,e-nop,f-nops) = auxWabef(a-nop,b-nops,e-nop,f-nops) + otherTs(b-nops,m)*auxspinints(a,m,e,f) !! termino anulado -otherTs(a-nop,m)*auxspinints(b,m,e,f))
					end do
               do m=1, numberOfParticles
						do n=1, numberOfOtherSpecieParticles
							auxWabef(a-nop,b-nops,e-nop,f-nops) = auxWabef(a-nop,b-nops,e-nop,f-nops) + 0.25*auxtau(a-nop,b-nops,m,n)*auxspinints(m,n,e,f)
						end do
					end do
				end do
			end do
		end do
	end do


!	!! Equation 8
! eq ok

	do m=1, numberOfParticles
		do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
			do e=numberOfParticles+1, numberOfContractions
				do jj=1, numberOfOtherSpecieParticles
					auxWmbej(m,b-nops,e-nop,jj) = auxspinints(m,b,e,jj)
					do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
						auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + otherTs(f-nops,jj)*auxspinints(m,b,e,f)
					end do
					do n=1, numberOfOtherSpecieParticles !numberOfParticles (anterior)
!! 26 de enero 2016
!! error  auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + (-Ts(e-nop,n)*auxspinints(m,n,e,f)) Ts deberia ser otherTs u auxspinints deberia llamar otros contadores
						auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + (-otherTs(b-nops,n)*auxspinints(m,n,e,jj)) !doble jj ahora n es otherparticle
!						auxWmbej(m,b,e,jj) = auxWmbej(m,b,e,jj) + (-Ts(b,n)*auxspinints(m,n,e,f))
						do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!! 26 de enero 2016 cambio en el contador de la matriz del termino Ts(b,n). Siempre ha sido virtuales - ocupados en electronico, pero aparece de interespecie.
!!							auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + -(0.5*auxTd(f-nop,b-nops,jj,n)+otherTs(f-nops,jj)*Ts(b,n))*auxspinints(m,n,e,f) ! 0.5 before auxTd
                     auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + -(0.5*auxTd(f-nops,b-nops,jj,n)+otherTs(f-nops,jj)*otherTs(b-nops,n))*auxspinints(m,n,e,f) ! 0.5 before auxTd
						end do
					end do
				end do
			end do
		end do
	end do 

!!!! End Intermediates for different particles

 !!!! Begin Energy Calculation

	auxECCSD = 0.0_8
 
!         do ii=1, numberOfParticles
!                 do a=numberOfParticles+1, numberOfContractions
!                         ECCSD1 = ECCSD1 + Fs%values(ii,a)*Ts(a,ii)
!                         do jj=1, numberOfParticles
!                                 do b=numberOfParticles+1, numberOfContractions
!                                         ECCSD1 = ECCSD1 + (0.25*spinints(ii,jj,a,b)*Td(a,b,ii,jj)+0.5*spinints(ii,jj,a,b)*Ts(a,ii)*Ts(b,jj))
!!					write (*,*) "CC electron", ECCSD1
!                                 end do
!                         end do
!                 end do
!         end do
!
!
!
!         do ii=1, numberOfOtherSpecieParticles
!                 do a=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                         ECCSD2 = ECCSD2 + otherFs%values(ii,a)*otherTs(a,ii)
!                         do jj=1, numberOfOtherSpecieParticles
!                                 do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                                         ECCSD2 = ECCSD2 + (0.25*otherspinints(ii,jj,a,b)*otherTd(a,b,ii,jj)+0.5*otherspinints(ii,jj,a,b)*otherTs(a,ii)*otherTs(b,jj))
!!					write (*,*) "CC positron", ECCSD2
!                                 end do
!                         end do
!                 end do
!         end do
!
	!!! TEST for CCSD combined

!        	do ii=1, numberOfParticles
!                	do aa=numberOfParticles+1, numberOfContractions
!                        auxECCSD = auxECCSD + Fs%values(ii,aa)*auxTs(aa,ii)
!			end do
!		end do
!        	do ii=1, numberOfOtherSpecieParticles
!                	do aa=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
!                        auxECCSD = auxECCSD + otherFs%values(ii,aa)*auxTs(aa,ii)
!			end do
!		end do


         do ii=1, numberOfParticles
                 do aa=numberOfParticles+1, numberOfContractions
                         do iii=1, numberOfOtherSpecieParticles
                                 do aaa=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
                                         auxECCSD = auxECCSD + (auxspinints(ii,iii,aa,aaa)*auxTd(aa-nop,aaa-nops,ii,iii)+auxspinints(ii,iii,aa,aaa)*CoupledCluster_instance%Tstest(aa-nop,ii)*otherTs(aaa-nops,iii)) ! 0.25 y 0.5
                                 end do
                         end do
                 end do
         end do

write(*,*) i,j,auxECCSD

!write (*,*) "CC inter", auxECCSD

  DECC = abs( auxECCSD - OLDCC )
!
!	write (*,*) "CC electron", ECCSD1
!	write (*,*) "CC other", ECCSD2
!	write (*,*) "CC e-/other", auxECCSD
!
!!!! End Energy Calculation


!!!! Let's make T1 and T2 (the new ones) 

	!! Equation 1
		
!	do a=numberOfParticles+1, numberOfContractions
!		do ii=1, numberOfParticles
!			TsNew(a,ii) = Fs%values(ii,a)
!			do e=numberOfParticles+1, numberOfContractions
!				TsNew(a,ii) = TsNew(a,ii) + Ts(e,ii)*Fae(a,e)
!			end do
!			do m=1, numberOfParticles
!				TsNew(a,ii) = TsNew(a,ii) + (-Ts(a,m)*Fmi(m,ii))
!				do e=numberOfParticles+1, numberOfContractions
!					TsNew(a,ii) = TsNew(a,ii) + Td(a,e,ii,m)*Fme(m,e)
!					do f=numberOfParticles+1, numberOfContractions
!						TsNew(a,ii) = TsNew(a,ii) + (-0.5*Td(e,f,ii,m)*spinints(m,a,e,f))
!					end do
!					do n=1, numberOfParticles
!						TsNew(a,ii) = TsNew(a,ii) + (-0.5*Td(a,e,m,n)*spinints(n,m,e,ii))
!					end do
!				end do
!			end do
!			do n=1,numberOfParticles
!				do f=numberOfParticles+1, numberOfContractions
!					TsNew(a,ii) = TsNew(a,ii) + (-Ts(f,n)*spinints(n,a,ii,f))
!				end do
!			end do
!			TsNew(a,ii) = TsNew(a,ii)/Dai(a,ii)
!			Ts(a,ii) = TsNew(a,ii)
!!			write(*,*) a,i,Ts(a,i),TsNew(a,i)
!		end do
!	end do
!
!!write (*,*) Ts(:,:)
!			
	!! Equation 2

	
	do a=numberOfParticles+1, numberOfContractions
		do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
			do ii=1, numberOfParticles
				do jj=1, numberOfOtherSpecieParticles
					auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + auxspinints(ii,jj,a,b) ! termino 1
!					do e=numberOfParticles+1, numberOfContractions
!						TdNew(a,b,ii,jj) = TdNew(a,b,ii,jj) + (Td(a,e,ii,jj)*Fae(b,e)-Td(b,e,ii,jj)*Fae(a,e))
!						do m=1, numberOfParticles
!							TdNew(a,b,ii,jj) = TdNew(a,b,ii,jj) + (-0.5*Td(a,e,ii,jj)*Ts(b,m)*Fme(m,e)+0.5*Td(a,e,ii,jj)*Ts(a,m)*Fme(m,e)) 
!						end do
!					end do
!					do m=1, numberOfParticles
!						TdNew(a,b,ii,jj) = TdNew(a,b,ii,jj) + (-Td(a,b,ii,m)*Fmi(m,jj)+Td(a,b,jj,m)*Fmi(m,ii))
!						do e=numberOfParticles+1, numberOfContractions
!							TdNew(a,b,ii,jj) = TdNew(a,b,ii,jj) + (-0.5*Td(a,b,ii,m)*Ts(e,j)*Fme(m,e)+0.5*Td(a,b,ii,m)*Ts(e,ii)*Fme(m,e))
!						end do
!					end do
					do e=numberOfParticles+1, numberOfContractions
!						auxTdNew(a,b,ii,jj) = auxTdNew(a,b,ii,jj) + (Ts(e,ii)*spinints(a,b,e,jj)-Ts(e,jj)*spinints(a,b,e,ii))
						do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
							auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + 0.5*auxtau(e-nop,f-nops,ii,jj)*auxWabef(a-nop,b-nops,e-nop,f-nops) ! termino 5
						end do
					end do
					do m=1, numberOfParticles

                  auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) - CoupledCluster_instance%Tstest(a-nop,m)*auxspinints(m,b,ii,jj) ! termino 10

!						auxTdNew(a,b,ii,jj) = auxTdNew(a,b,ii,jj) + (-Ts(a,m)*spinints(m,b,ii,jj)+Ts(b,m)*spinints(m,a,ii,jj))
!						do e=numberOfParticles+1, numberOfContractions
!							auxTdNew(a,b,ii,jj) = auxTdNew(a,b,ii,jj) + Td(a,e,ii,m)*Wmbej(m,b,e,jj) - Ts(e,ii)*Ts(a,m)*spinints(m,b,e,jj) + -Td(a,e,jj,m)*Wmbej(m,b,e,ii) + Ts(e,jj)*Ts(a,m)*spinints(m,b,e,ii) + -Td(b,e,ii,m)*Wmbej(m,a,e,jj) - Ts(e,ii)*Ts(b,m)*spinints(m,a,e,jj) + Td(b,e,jj,m)*Wmbej(m,a,e,ii) - Ts(e,jj)*Ts(b,m)*spinints(m,a,e,ii)
!						end do
						do n=1, numberOfOtherSpecieParticles
							auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + 0.5*auxtau(a-nop,b-nops,m,n)*auxWmnij(m,n,ii,jj)  ! termino 4
						end do
					end do
    
      !! 02 febrero


               do e=nops+1, nocs
                  auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + auxTd(a-nop,e-nops,ii,jj)*auxFae(b-nops,e-nops)    ! parte 1 termino 2
                  do m=1, nops
                  auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + (-0.5*(otherTs(b-nops,m)*otherFme(m,e-nops)))     ! parte 2 termino 2
                  
                  end do
               end do




       !! 29 de enero 2016
               do m=1, numberOfParticles
                  do e=numberOfParticles+1, numberOfContractions
       ! Correlacion e-e- ---> e+ (posible?) 
                     auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + CoupledCluster_instance%Tdtest(a-nop,e-nop,ii,m)*auxWmbej(m,b-nops,e-nop,jj)  ! termino 6 revisar


                     auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) - CoupledCluster_instance%Tstest(e-nop,ii)*CoupledCluster_instance%Tstest(a-nop,m)*auxspinints(m,b,e,jj) ! termino 7 

!auxTd(b-nops,e-nop,ii,m) para continuar con las dimensiones de la matrix auxTd y por reglas de simetria el valor de Td para la proxima ecuacion puede ser: auxTd(e-nop,b-nops,m,ii)

!                    auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + auxTd(b-nops,e-nop,jj,m)*CoupledCluster_instance%Wmbejtest(m,a-nop,e-nop,ii) ! termino 8 revisar
                  end do
               end do

               do e=numberOfParticles+1, numberOfContractions
                  auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + CoupledCluster_instance%Tstest(e-nop,ii)*auxspinints(m,b,e,jj)  ! termino 9
               end do


!! Make a auxDabij denominator | It's the same Dabij denominator that is used in interspecies part
!!                                         auxDabij(aa,aaa,ii,iii) = Fs%values(ii,ii)+otherFs%values(iii,iii)-Fs%values(aa,aa)-otherFs%values(aaa,aaa)

!!             auxTdNew(a-nop,b,ii,jj) = auxTdNew(a-nop,b,ii,jj)/auxDabij(a,b,ii,jj)
               
               auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj)/(Fs%values(ii,ii)+otherFs%values(jj,jj)-Fs%values(a,a)-otherFs%values(b,b))
					auxTd(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj)
!                                        auxtaus(a,b,ii,jj) = auxTd(a,b,ii,jj) + 0.5*(Ts(a,ii)*otherTs(b,jj) - Ts(b,ii)*Ts(a,jj))
!                                        auxtau(a,b,ii,jj) = auxTd(a,b,ii,jj) + Ts(a,ii)*Ts(b,jj) - Ts(b,ii)*Ts(a,jj)
                                        auxtaus(a-nop,b-nops,ii,jj) = auxTd(a-nop,b-nops,ii,jj) + 0.5*(CoupledCluster_instance%Tstest(a-nop,ii)*otherTs(b-nops,jj))
                                        auxtau(a-nop,b-nops,ii,jj) = auxTd(a-nop,b-nops,ii,jj) + CoupledCluster_instance%Tstest(a-nop,ii)*otherTs(b-nops,jj)


				end do
			end do
		end do
	end do





!! Equation 8 electronica
! eq ok

!!**  	do m=1, numberOfParticles
!!**  		do b=numberOfParticles+1, numberOfContractions
!!**  			do e=numberOfParticles+1, numberOfContractions
!!**  				do jj=1, numberOfParticles
!!**  			!		CoupledCluster_instance%Wmbejtest(m,b-nop,e-nop,jj) = spinints(m,b,e,jj)
!!**  					do f=numberOfParticles+1, numberOfContractions
!!**  						Wmbej(m,b-nop,e-nop,jj) = Wmbej(m,b-nop,e-nop,jj) + Ts(f-nop,jj)*spinints(m,b,e,f)
!!**  					end do
!!**  					do n=1, numberOfParticles
!!**  !! 26 de enero 2016
!!**                     Wmbej(m,b-nop,e-nop,jj) = Wmbej(m,b-nop,e-nop,jj) + (-Ts(b-nop,n)*spinints(m,n,e,jj))
!!**  						do f=numberOfParticles+1, numberOfContractions
!!**  							Wmbej(m,b-nop,e-nop,jj) = Wmbej(m,b-nop,e-nop,jj) + (-(0.5*Td(f-nop,b-nop,jj,n)+Ts(f-nop,jj)*Ts(b-nop,n))*spinints(m,n,e,f))
!!**  						end do
!!**  					end do
!!**  				end do
!!**  			end do
!!**  		end do
!!**  	end do









end do !! Loop

	coupledClusterValue%values(mmm) = auxECCSD

      end do !! Number Of
  end do     !! Species
!
	call Matrix_destructor(auxMatrix)
!	CoupledCluster_instance%mp2Correction = sum ( CoupledCluster_instance%mp2Coupling%values )
	CoupledCluster_instance%ccsdTwoParticlesCorrection = auxECCSD
	call Vector_show(coupledClusterValue)

end if

!!*********************************************************************************************************************************
!! End CCSD Calculation different species
!!*********************************************************************************************************************************

   close(wfnUnit)



 end subroutine CoupledCluster_iterateIntermediates_DiffSpecies





 
 !**
 ! @ Returns final energy (with CCSD correction)
 !**

! function CoupledCluster_getTotalEnergy() result(output)
!	implicit none
!	real(8) :: output
!
!	output = CoupledCluster_instance%totalEnergy
!
! end function CoupledCluster_getTotalEnergy

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CoupledCluster_getTransformedIntegrals()
    implicit none

!    type(TransformIntegrals) :: repulsionTransformer
    integer :: numberOfSpecies
    integer :: i,j,m,n,mu,nu
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecie
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    type(Matrix) :: auxMatrix
    type(Matrix) :: molecularCouplingMatrix
    type(Matrix) :: molecularExtPotentialMatrix
    type(Matrix) :: couplingMatrix
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(CoupledCluster_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(CoupledCluster_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))

!    print *,""
!    print *,"BEGIN INTEGRALS TRANFORMATION:"
!    print *,"========================================"
!    print *,""
!    print *,"--------------------------------------------------"
!    print *,"    Algorithm Four-index integral tranformation"
!    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!    print *,"  Computer Physics Communications, 2005, 166, 58-65"
!    print *,"--------------------------------------------------"
!    print *,""
!
!    call TransformIntegrals_constructor( repulsionTransformer )

    do i=1, numberOfSpecies
        nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
        specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
        ocupationNumber = MolecularSystem_getOcupationNumber( i )
        numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )

!        write (6,"(T10,A)")"ONE PARTICLE INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
        call Matrix_constructor (CoupledCluster_instance%twoCenterIntegrals(i), int(numberOfContractions,8), int(numberOfContractions,8) )
        call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8))
        call Matrix_constructor (couplingMatrix,int(numberOfContractions,8), int(numberOfContractions,8))

        hcoreMatrix  = HartreeFock_instance%HcoreMatrix 

!        hcoreMatrix%values = WaveFunction_HF_instance( specieID )%IndependentParticleMatrix%values

        !! Open file for wavefunction

        wfnFile = "lowdin.wfn"
        wfnUnit = 20

        open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

        arguments(2) = MolecularSystem_getNameOfSpecie(i)
        arguments(1) = "COEFFICIENTS"

        coefficients = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        arguments(1) = "HCORE"

        hcoreMatrix = &
          Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
          columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

        do m=1,numberOfContractions
          do n=m, numberOfContractions
             do mu=1, numberOfContractions
                do nu=1, numberOfContractions
                    CoupledCluster_instance%twoCenterIntegrals(i)%values(m,n) = &
                        CoupledCluster_instance%twoCenterIntegrals(i)%values(m,n) + &
                        coefficients%values(mu,m)* &
                        coefficients%values(nu,n)* &
                        hcoreMatrix%values(mu,nu)
                end do
             end do
          end do
       end do

!! Not implemented yet
!!       if( WaveFunction_HF_instance( specieID )%isThereExternalPotential ) then
!!          do m=1,numberOfContractions
!!             do n=m, numberOfContractions
!!                do mu=1, numberOfContractions
!!                   do nu=1, numberOfContractions
!!                      CoupledCluster_instance%twoCenterIntegrals(i)%values(m,n) = &
!!                           CoupledCluster_instance%twoCenterIntegrals(i)%values(m,n) + &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(mu,m)* &
!!                           WaveFunction_HF_instance( specieID )%waveFunctionCoefficients%values(nu,n) * &
!!                           WaveFunction_HF_instance( specieID )%ExternalPotentialMatrix%values(mu,nu)
!!                   end do
!!                end do
!!             end do
!!          end do
!!       end if

       do m=1,numberOfContractions
          do n=m, numberOfContractions
             CoupledCluster_instance%twoCenterIntegrals(i)%values(n,m)=&
                  CoupledCluster_instance%twoCenterIntegrals(i)%values(m,n)
          end do
       end do
          
       ! print *, "Independent Particle"
       ! call Matrix_show ( CoupledCluster_instance%twoCenterIntegrals(i) )

!       write (6,"(T10,A)")"TWO PARTICLES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
!       print *,""

!       call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer, MolecularSystem_getEigenvectors(specieID), &
!            CoupledCluster_instance%fourCenterIntegrals(i,i), specieID, trim(nameOfSpecie) )

        call ReadTransformedIntegrals_readOneSpecies( specieID, CoupledCluster_instance%fourCenterIntegrals(i,i)   )

       if ( numberOfSpecies > 1 ) then
          do j = 1 , numberOfSpecies
             if ( i .ne. j) then
                nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
                otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
                ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
                numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )

!                write (6,"(T10,A)") "INTER-SPECIES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)//"/"//trim(nameOfOtherSpecie)
!                print *,""

!                call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!                     MolecularSystem_getEigenVectors(i), MolecularSystem_getEigenVectors(j), &
!                     CoupledCluster_instance%fourCenterIntegrals(i,j), specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )

                 call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         CoupledCluster_instance%fourCenterIntegrals(i,j) )

             end if
          end do
       end if
    end do
!    print *,"END INTEGRALS TRANFORMATION:"
    close (wfnUnit)
	call Matrix_destructor (hcoreMatrix)
        call Matrix_destructor (couplingMatrix)

  end subroutine CoupledCluster_getTransformedIntegrals

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine CoupledCluster_exception( typeMessage, description, debugDescription)
    implicit none
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription

    type(Exception) :: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine CoupledCluster_exception

end module CoupledCluster_
