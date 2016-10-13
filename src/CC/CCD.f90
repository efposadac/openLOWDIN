!!******************************************************************************
!!
!!                       CCD-APMO
!!
!!******************************************************************************

module CCD_
  use Exception_
  use Matrix_
  use Vector_
  !use Coupled_
  use ReadTransformedIntegrals_
  use MolecularSystem_
  use String_
  use IndexMap_
  implicit none

  !>
  !! @brief Coupled Interaction Module, works in spin orbitals
  !!
  !! @author Carlos Andrés Ortiz Mahecha (CAOM)
  !!
  !! CCD Module:
  !! Intraespecie
  !! Interespecie
  !! Based on CCSD Code (Alejandro @author)
  !!
  !! <b> Creation data : 05 February</b> 2016
  !!
  !! <b> History change: </b>
  !! 
  !!
  !<
  type, public :: CCD
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

  end type CCD

  type, public :: HartreeFockD
        real(8) :: totalEnergy
        real(8) :: puntualInteractionEnergy
        type(matrix) :: coefficientsofcombination 
        type(matrix) :: HcoreMatrix 
  end type HartreeFockD
  
  type(CCD) :: CCD_instance
  type(HartreeFockD) :: HartreeFockD_instance

  public :: &
       CCD_constructor, &
       CCD_destructor, &
!       CCD_getTotalEnergy, &
       CCD_run, &
       CCD_show, &
       CCD_exception

  private   
contains

  double precision function logic2dbl(a)
     logical, intent(in) :: a

      if (a) then
         logic2dbl = 1.d0
      else
         logic2dbl = 0.d0
      end if
  end function logic2dbl
  
  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  !<
  subroutine CCD_constructor(level)
    implicit none
    character(*) :: level

!
    integer :: numberOfContractionss,nopp,numberOfParticless,nocc,i

    CCD_instance%isInstanced=.true.
    CCD_instance%level=level

  end subroutine CCD_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine CCD_destructor()
    implicit none

    CCD_instance%isInstanced=.false.

  end subroutine CCD_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CCD_show()
    implicit none
    type(CCD) :: this
    integer :: i,j,k,l
    integer :: numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
 
    if ( CCD_instance%isInstanced ) then

       print *,""
       print *," POST HARTREE-FOCK CALCULATION"
       print *," COUPLED CLUSTER THEORY:"
       print *,"=============================="
       print *,""
       write (6,"(T10,A30, A5)") "LEVEL = ", CCD_instance%level
       write (6,"(T10,A30, F20.12)") "HF ENERGY = ", HartreeFockD_instance%totalEnergy
       write (6,"(T10,A30, F20.12)") "MP2 CORR. ENERGY = ", CCD_instance%secondOrderCorrection
       write (6,"(T10,A30, F20.12)") "CCD CORR. ENERGY = ", CCD_instance%ccsdCorrection
       write (6,"(T10,A30, F20.12)") "============================================================"
       write (6,"(T10,A30, F20.12)") "Total Energy (HF+CCD) = ", HartreeFockD_instance%totalEnergy+CCD_instance%ccsdCorrection
    

  if ( numberOfSpecies > 1 ) then

       print *,""
       print *,""
       print *," COUPLED CLUSTER DOUBLES THEORY FOR DIFF. SPECIES"
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
          "} = ", CCD_instance%mp2Coupling%values(k)
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( i ) ), &
          "} = ",  CCD_instance%coupledClusterSDCorrection%values( i )
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( j ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
          "} = ",  CCD_instance%coupledClusterSDCorrection%values( j )
                                   write (*,'(T10,A30,A8,A4,ES16.8)') "CCD Corr.{", &
                                             trim(  MolecularSystem_getNameOfSpecie( i ) ) // "/" // trim(  MolecularSystem_getNameOfSpecie( j ) ), &
          "} = ", CCD_instance%ccsdTwoParticlesCorrection
                         end do
                    end do

       print *,""
       print *,""
           write (6,'(T10,A30,F20.8)') "========================================================="
           write (6,'(T10,A30,F20.8)') "MP2 Total Energy = ", HartreeFockD_instance%totalEnergy+CCD_instance%secondOrderCorrection+CCD_instance%mp2Coupling%values 
           write (6,'(T10,A30,F20.8)') "CCD Total Energy = ", HartreeFockD_instance%totalEnergy+CCD_instance%ccsdCorrection+CCD_instance%ccsdTwoParticlesCorrection 

  end if

    else 

    end if

  end subroutine CCD_show

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CCD_run()
    implicit none 
    integer :: m
    real(8), allocatable :: eigenValues(:) 

       print *, ""
       print *, ""
       print *, "==============================================="
       print *, "|            BEGIN CCD CALCULATION            |"
       print *, "-----------------------------------------------"
       print *, ""

  print *, "  ________________________________________________ "
  print *, " |Initial Guess... From MP2 Calculation           |"
  print *, " |------------------------------------------------|"
  print *, " |Iteration of Coupled Cluster intermediates . . .|"
  print *, " |------------------------------------------------|"
  print *, " |Computing CCD Energy...                         |"
  print *, " |________________________________________________|"
  print *, ""
  print *, "Patience is not the ability to wait, but the ability to keep a good attitude while waiting... Or just destroy the computer."

  call CCD_iterateIntermediates_SameSpecies()
  
  print *, "mi variable copartida Td",  CCD_instance%Tdtest(1,1,1,1) 

   call CCD_iterateIntermediates_DiffSpecies()
 
  print *, "mi variable copartida Td",  CCD_instance%Tdtest(1,1,1,1) *2

       print *,""
       print *, "-----------------------------------------------"
       print *, "|              END CCD CALCULATION            |"
       print *, "==============================================="
       print *, ""
         

  end subroutine CCD_run

 !>
  !! @brief Calculation of the intermediates
  !!
  !<
 subroutine CCD_iterateIntermediates_SameSpecies()
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
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFockD_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFockD_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


!!******************************************************************************************************************************
!! Begin CCD calculation same specie. 
!!******************************************************************************************************************************



! speciesID=1 !!! FOR NOW, only electrons 
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

  call Vector_constructor( CCD_instance%coupledClusterSDCorrection, numberOfSpecies)
  call Vector_constructor( CCD_instance%energyCorrectionOfSecondOrder, numberOfSpecies)
   
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

                CCD_instance%energyCorrectionOfSecondOrder%values(speciesID) = independentEnergyCorrection &
                        * ( ( MolecularSystem_getCharge( speciesID ) )**4.0_8 )

    !! Suma las correcciones de energia para especies independientes
    CCD_instance%secondOrderCorrection = sum( CCD_instance%energyCorrectionOfSecondOrder%values ) !MP2 Corr. Energy


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
                                        auxVal1 = auxVal_A * logic2dbl(mod(p,2) == mod(r,2)) * logic2dbl(mod(q,2) == mod(s,2))
                                        auxVal2 = auxVal_B * logic2dbl(mod(p,2) == mod(s,2)) * logic2dbl(mod(q,2) == mod(r,2))
                                        spinints(p,q,r,s) = auxVal1 - auxVal2 !! p+1 o p-1? Revisar !! ecuacion 1
!         write (*,*) spinints(p,q,r,s)
                                end do
                        end do
                end do
        end do


!!!! Initial guesses T1 and T2



!variables compartidas

    if (allocated(CCD_instance%Tstest)) deallocate(CCD_instance%Tstest)
    allocate(CCD_instance%Tstest(noc-nop,nop)) ! 01f
    CCD_instance%Tstest(:,:) = 0.0_8


    if (allocated(CCD_instance%Tdtest)) deallocate(CCD_instance%Tdtest)
    allocate(CCD_instance%Tdtest(noc-nop,noc-nop,nop,nop)) ! 01f
    CCD_instance%Tdtest(:,:,:,:) = 0.0_8


    if (allocated(CCD_instance%Wmbejtest)) deallocate(CCD_instance%Wmbejtest)
    allocate(CCD_instance%Wmbejtest(nop,noc-nop,noc-nop,nop)) ! 01f
    CCD_instance%Wmbejtest(:,:,:,:) = 0.0_8


 

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


                                        CCD_instance%Tdtest(a-nop,b-nop,i,j) = CCD_instance%Tdtest(a-nop,b-nop,i,j) + (spinints(i,j,a,b)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)))
                                        taus(a-nop,b-nop,i,j) = CCD_instance%Tdtest(a-nop,b-nop,i,j)
                                        tau(a-nop,b-nop,i,j) = CCD_instance%Tdtest(a-nop,b-nop,i,j) 
                                        
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
!     write(*,*) a,i,Dai(a,i)
                end do
        end do
        

!!!! End Denominators

!!!! Main Loop



  ECCSD = 0.0_8
  DECC = 1.0_8 
!  OLDCC = 0 

!   do while (OLDCC < 5)
  do while (DECC >= 1.0D-8)
! OLDCC = OLDCC + 1
  OLDCC = ECCSD


  if (allocated(Fae)) deallocate (Fae)
        allocate(Fae(noc-nop,noc-nop))
   Fae=0.0_8
  if (allocated(Fmi)) deallocate (Fmi)
        allocate(Fmi(nop,nop))
  Fmi=0.0_8
  if (allocated(Fme)) deallocate (Fme)
          allocate(Fme(nop,noc-nop))
   Fme=0.0_8

  if (allocated(Wmnij)) deallocate (Wmnij)
        allocate(Wmnij(nop,nop,nop,nop))
  Wmnij=0.0_8
  if (allocated(Wabef)) deallocate (Wabef)
        allocate(Wabef(noc-nop,noc-nop,noc-nop,noc-nop))
  Wabef=0.0_8
  if (allocated(Wmbej)) deallocate (Wmbej)
        allocate(Wmbej(nop,noc-nop,noc-nop,nop))
  Wmbej=0.0_8

  if (allocated(TsNew)) deallocate (TsNew)
   allocate(TsNew(noc-nop,nop))
  TsNew=0.0_8
  if (allocated(TdNew)) deallocate (TdNew)
        allocate(TdNew(noc-nop,noc-nop,nop,nop))
  TdNew=0.0_8



!!!! Intermediates

  !! Equation 3
! eq OK

  do a=numberOfParticles+1, numberOfContractions
    do e=numberOfParticles+1, numberOfContractions
         if (a==e) then 
            kro = 1 
         else 
            kro = 0 
         end if
         
         Fae(a-nop,e-nop) = (1 - kro)*Fs%values(a,e)
      do m=1, numberOfParticles
        do f=numberOfParticles+1, numberOfContractions
          do n=1, numberOfParticles
            Fae(a-nop,e-nop) = Fae(a-nop,e-nop) + (-0.5*taus(a-nop,f-nop,m,n)*spinints(m,n,e,f))
          end do
        end do
      end do

    end do
  end do


  !! Equation 4

  do m=1, numberOfParticles
    do i=1, numberOfParticles
         if (m==i) then 
            kro = 1 
         else 
            kro = 0 
         end if
         
         Fmi(m,i) = (1 - kro)*Fs%values(m,i)
      do e=numberOfParticles+1, numberOfContractions
        do n=1, numberOfParticles
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
          CCD_instance%Wmbejtest(m,b-nop,e-nop,j) = spinints(m,b,e,j)
          do n=1, numberOfParticles
            do f=numberOfParticles+1, numberOfContractions
              CCD_instance%Wmbejtest(m,b-nop,e-nop,j) = CCD_instance%Wmbejtest(m,b-nop,e-nop,j) + -(0.5*CCD_instance%Tdtest(f-nop,b-nop,j,n)*spinints(m,n,e,f))
                end do
          end do
        end do
      end do
    end do
  end do

!!!! End Intermediates


auxECCSD = 0.0_8


    do a=numberOfParticles+1, numberOfContractions
        do i=1, numberOfParticles
          do j=i+1, numberOfParticles
            do b=a+1, numberOfContractions
              auxECCSD = auxECCSD + spinints(i,j,a,b)*CCD_instance%Tdtest(a-nop,b-nop,i,j) 
            end do
          end do
       end do
     end do

  ECCSD = auxECCSD
  
  DECC = abs( ECCSD - OLDCC )

   write (*,*) speciesID 
   write (*,*) ECCSD 


  !! Equation 2

  
  do a=numberOfParticles+1, numberOfContractions
    do b=numberOfParticles+1, numberOfContractions
      do i=1, numberOfParticles
        do j=1, numberOfParticles

          TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + spinints(i,j,a,b) !A
  
   ! 1er ciclo
               do e=numberOfParticles+1, numberOfContractions
  
                  TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (CCD_instance%Tdtest(a-nop,e-nop,i,j)*Fae(b-nop,e-nop)-CCD_instance%Tdtest(b-nop,e-nop,i,j)*Fae(a-nop,e-nop)) !B

            do f=numberOfParticles+1, numberOfContractions
              TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + 0.5*tau(e-nop,f-nop,i,j)*Wabef(a-nop,b-nop,e-nop,f-nop) !D
            end do

          end do
  ! 2do ciclo
               do m=1, numberOfParticles
            TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + (-CCD_instance%Tdtest(a-nop,b-nop,i,m)*Fmi(m,j)+CCD_instance%Tdtest(a-nop,b-nop,j,m)*Fmi(m,i)) !C


            do n=1, numberOfParticles
              TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + 0.5*tau(a-nop,b-nop,m,n)*Wmnij(m,n,i,j) !E
            end do

            do e=numberOfParticles+1, numberOfContractions

                  TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + CCD_instance%Tdtest(a-nop,e-nop,i,m)*CCD_instance%Wmbejtest(m,b-nop,e-nop,j) - CCD_instance%Tdtest(a-nop,e-nop,j,m)*CCD_instance%Wmbejtest(m,b-nop,e-nop,i) - CCD_instance%Tdtest(b-nop,e-nop,i,m)*CCD_instance%Wmbejtest(m,a-nop,e-nop,j) + CCD_instance%Tdtest(b-nop,e-nop,j,m)*CCD_instance%Wmbejtest(m,a-nop,e-nop,i)  !F

                  end do
          end do


!Nueva ecuacion
!!***               do m=1, numberOfParticles
!!***                  do n=1, numberOfParticles
!!***                     do e=numberOfParticles+1, numberOfContractions
!!***                        do f=numberOfParticles+1, numberOfContractions
!!***                           TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j) + spinints(i,j,a,b)*(0.25*taus(e-nop,f-nop,i,j)*taus(a-nop,b-nop,m,n) - 0.5*(taus(a-nop,e-nop,i,j)*taus(b-nop,f-nop,m,n) + taus(b-nop,f-nop,i,j)*taus(a-nop,e-nop,m,n)) - 0.5*(taus(a-nop,b-nop,i,m)*taus(e-nop,f-nop,j,n) + taus(e-nop,f-nop,i,m)) + taus(a-nop,e-nop,i,m)*taus(b-nop,f-nop,j,n) + taus(b-nop,f-nop,i,m)*taus(a-nop,e-nop,j,n))
!!***                        end do
!!***                     end do
!!***                  end do
!!***              end do



!!* Eq 13 Stanton

!!  22 enero 2016       
!!!! Make denominator array Dabij     |   Dabij(a-nop,b-nop,i,j) = Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b)
               
               TdNew(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j)/(Fs%values(i,i)+Fs%values(j,j)-Fs%values(a,a)-Fs%values(b,b))
               CCD_instance%Tdtest(a-nop,b-nop,i,j) = TdNew(a-nop,b-nop,i,j)
                                        taus(a-nop,b-nop,i,j) = CCD_instance%Tdtest(a-nop,b-nop,i,j)
                                        tau(a-nop,b-nop,i,j) = CCD_instance%Tdtest(a-nop,b-nop,i,j) 
        end do
      end do
    end do
  end do


!!!! End New T1 and T2

  CCD_instance%coupledClusterSDCorrection%values(speciesID) = ECCSD ! CCD Correlation Energy

 end do !! Main Loop



  end if !! If numberOfParticles > 1  


 end do !! Species


  call Vector_show (CCD_instance%coupledClusterSDCorrection)
  CCD_instance%ccsdCorrection = sum( CCD_instance%coupledClusterSDCorrection%values ) !CCD Corr. Energy
 
! 28 de enero 2016

! if (allocated(Ts)) deallocate (Ts)
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
!! End CCD calculation same specie!
!!*********************************************************************************************************************************

   close(wfnUnit)

 end subroutine CCD_iterateIntermediates_SameSpecies


!!*********************************************************************************************************************************
!! Begin CCD calculation different species!
!!*********************************************************************************************************************************



 subroutine CCD_iterateIntermediates_DiffSpecies()
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
   real(8) :: auxVal_A,auxVal_AA,auxVal_AA1
   real(8) :: auxVal_B,auxVal_BB,auxVal_BB1
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
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFockD_instance%totalEnergy, &
         arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=HartreeFockD_instance%puntualInteractionEnergy, &
         arguments=["PUNTUALINTERACTIONENERGY"])


    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()


if ( numberOfSpecies > 1 ) then

mp2CouplingCorrection = 0.0_8
mmm = 0

print *, " "
print *, " "
print *, " |------------------------------------------------|"
print *, " |Calculating CCD corr. for different species.....|"
print *, " |------------------------------------------------|"



!   CCD_instance%Tstest = 136 * 100

!   print *, "mi variable copartida *100 ",  CCD_instance%Tstest(1,1) 

f = numberOfSpecies * ( numberOfSpecies-1 ) / 2

call Vector_constructor( CCD_instance%mp2Coupling, f)
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

        CCD_instance%mp2Coupling%values(mmm)= ( lambda * lambdaOfOtherSpecie * mp2CouplingCorrection )

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


        if (allocated(auxspinints)) deallocate (auxspinints)
        allocate(auxspinints(numberOfContractions,numberOfContractionsOfOtherSpecie,numberOfContractions,numberOfContractionsOfOtherSpecie))


  !!!! TEST for combined spinints
! call ReadTransformedIntegrals_readTwoSpecies( speciesID, otherSpeciesID, auxspinints)
        do p=1, numberOfContractions
                do q=1, numberOfContractionsOfOtherSpecie
                        do r=1, numberOfContractions
                                do s=1, numberOfContractionsOfOtherSpecie
                                        value1 = IndexMap_tensorR4ToVector((p+1)/2,(r+1)/2,(q+1)/2,(s+1)/2,numberOfContractions/2, numberOfContractionsOfOtherSpecie/2)
                                        auxVal_A= auxMatrix%values(value1, 1)
                                        auxVal1 = auxVal_A * logic2dbl(mod(p,2) == mod(r,2)) * logic2dbl(mod(q,2) == mod(s,2))
                                        auxspinints(p,q,r,s) = auxVal1
                                end do
                        end do
                end do
        end do

!!!! Initial guess T2

        if (allocated(auxTd)) deallocate (auxTd)
        allocate(auxTd(noc-nop,nocs-nops,nop,nops))
        auxTd(:,:,:,:) = 0.0_8

        if (allocated(auxtaus)) deallocate (auxtaus)
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
!                       write(*,*) a,b,ii,jj,otherTd(a,b,i,j)
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
                                        auxtaus(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii)
                                        auxtau(aa-nop,aaa-nops,ii,iii) = auxTd(aa-nop,aaa-nops,ii,iii) 
                                end do
                        end do
                end do
        end do


!!!! End Initial Guess

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
!!          write(*,*) a,b,i,j,Dabij(a,b,i,j)
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
! h=h+1
        OLDCC = auxECCSD


       if (allocated(otherFae)) deallocate (otherFae)
        allocate(otherFae(noc-nop,noc-nop))
       otherFae=0.0_8
       if (allocated(otherFmi)) deallocate (otherFmi)
        allocate(otherFmi(nop,nop))
       otherFmi=0.0_8
!      if (allocated(otherFme)) deallocate (otherFme)
!       allocate(otherFme(nop,nop-noc))
!      otherFme=0.0_8
       if (allocated(Wmbej)) deallocate (Wmbej)
        allocate(Wmbej(nop,noc-nop,noc-nop,nop))
       Wmbej=0.0_8
        if (allocated(auxWmnij)) deallocate (auxWmnij)
        allocate(auxWmnij(nop,nocs,nop,nocs))
        auxWmnij=0.0_8
        if (allocated(auxWabef)) deallocate (auxWabef)
        allocate(auxWabef(noc-nop,nocs-nops,noc-nop,nocs-nops))
        auxWabef=0.0_8
        if (allocated(auxWmbej)) deallocate (auxWmbej)
        allocate(auxWmbej(nop,nocs-nops,noc-nop,nops))
        auxWmbej=0.0_8

        if (allocated(auxTdNew)) deallocate (auxTdNew)
        allocate(auxTdNew(noc-nop,nocs-nops,nop,nops))
        auxTdNew=0.0_8



!!!! Intermediates

  !! Equation 3
! eq OK

  do a=numberOfParticles+1, numberOfContractions
    do e=numberOfParticles+1, numberOfContractions
! Terminos ya incluidos en intraespecie
!!**         if (a==e) then 
!!**            kro = 1 
!!**         else 
!!**            kro = 0 
!!**         end if
!!**         
!!**         Fae(a-nop,e-nop) = (1 - kro)*Fs%values(a,e)
      do m=1, numberOfParticles
        do n=1, nops
               do f=nops+1, nocs
            otherFae(a-nop,e-nop) = otherFae(a-nop,e-nop) + (-0.5*auxtaus(a-nop,f-nops,m,n)*auxspinints(m,n,e,f))
          end do
        end do
      end do

    end do
  end do


  !! Equation 4

  do m=1, numberOfParticles
    do ii=1, numberOfParticles
! Terminos ya incluidos en intraespecie
!         if (m==i) then 
!            kro = 1 
!         else 
!            kro = 0 
!         end if
!         
!         Fmi(m,i) = (1 - kro)*Fs%values(m,i)
         do n=1, nops
             do e=numberOfParticles+1, numberOfContractions
             do f=nops+1, nocs
            otherFmi(m,ii) = otherFmi(m,ii) + 0.5*auxtaus(e-nop,f-nops,ii,n)*auxspinints(m,n,e,f)
           end do
         end do
       end do

    end do
  end do


! !! Equation 5
!! eq ok
! Terminos ya incluidos en intraespecie
! do m=1, numberOfParticles
!   do e=numberOfParticles+1, numberOfContractions
!     otherFme(m,e-nop) = otherFs%values(m,e)
!   end do
! end do



!! Eq 6
! eq ok

  do m=1, numberOfParticles
    do n=1, numberOfOtherSpecieParticles
      do ii=1, numberOfParticles
        do jj=1, numberOfOtherSpecieParticles
          auxWmnij(m,n,ii,jj) = auxspinints(m,n,ii,jj)
               do e=numberOfParticles+1, numberOfContractions
                  do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
              auxWmnij(m,n,ii,jj) = auxWmnij(m,n,ii,jj) + 0.25*auxtau(e-nop,f-nops,ii,jj)*auxspinints(m,n,e,f)
            end do
          end do
        end do
      end do
    end do
  end do


  !! Equation 7

! eq ok


  do a=numberOfParticles+1, numberOfContractions
    do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
      do e=numberOfParticles+1, numberOfContractions
        do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
          auxWabef(a-nop,b-nops,e-nop,f-nops) = auxspinints(a,b,e,f)
          do m=1, numberOfParticles
            do n=1, numberOfOtherSpecieParticles
              auxWabef(a-nop,b-nops,e-nop,f-nops) = auxWabef(a-nop,b-nops,e-nop,f-nops) + 0.25*auxtau(a-nop,b-nops,m,n)*auxspinints(m,n,e,f)
            end do
          end do
        end do
      end do
    end do
  end do


!! Equation 8
! eq ok

  do m=1, numberOfParticles
    do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
      do e=numberOfParticles+1, numberOfContractions
        do jj=1, numberOfOtherSpecieParticles
          auxWmbej(m,b-nops,e-nop,jj) = auxspinints(m,b,e,jj)
          do n=1, numberOfOtherSpecieParticles 
            do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
                     auxWmbej(m,b-nops,e-nop,jj) = auxWmbej(m,b-nops,e-nop,jj) + -(0.5*auxTd(f-nops,b-nops,jj,n)*auxspinints(m,n,e,f)) ! Revisar (Podria comentarse?)
            end do
          end do
        end do
      end do
    end do
  end do 

!!!! End Intermediates for different particles

  auxECCSD = 0.0_8
 
  
            do ii=1, numberOfParticles
              do iii=1, numberOfOtherSpecieParticles
                 do aa=numberOfParticles+1, numberOfContractions
                   do aaa=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
                         auxECCSD = auxECCSD + auxspinints(ii,iii,aa,aaa)*auxTd(aa-nop,aaa-nops,ii,iii)
                   end do
                 end do
               end do
           end do


write(*,*) i,j,auxECCSD

!write (*,*) "CC inter", auxECCSD

  DECC = abs( auxECCSD - OLDCC )
!
! write (*,*) "CC electron", ECCSD1
! write (*,*) "CC other", ECCSD2
! write (*,*) "CC e-/other", auxECCSD
!
!!!! End Energy Calculation


!!!! Let's make T1 and T2 (the new ones) 

  !! Equation 2

  
  do a=numberOfParticles+1, numberOfContractions
    do b=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
      do ii=1, numberOfParticles
        do jj=1, numberOfOtherSpecieParticles
          auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + auxspinints(ii,jj,a,b) ! termino 1
  
               do e=numberOfParticles+1, numberOfContractions
            do f=numberOfOtherSpecieParticles+1, numberOfContractionsOfOtherSpecie
              auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + 0.5*auxtau(e-nop,f-nops,ii,jj)*auxWabef(a-nop,b-nops,e-nop,f-nops) ! termino 5
            end do
          end do
          do m=1, numberOfParticles
            do n=1, numberOfOtherSpecieParticles
              auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + 0.5*auxtau(a-nop,b-nops,m,n)*auxWmnij(m,n,ii,jj)  ! termino 4
            end do
          end do
               do m=1, numberOfParticles
                  do e=numberOfParticles+1, numberOfContractions
       ! Correlacion e-e- ---> e+ (posible?) 
                     auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + CCD_instance%Tdtest(a-nop,e-nop,ii,m)*auxWmbej(m,b-nops,e-nop,jj)  ! termino 6 revisado

                     auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj) + auxTd(b-nops,e-nop,jj,m)*auxWmbej(m,b-nops,e-nop,jj) ! termino 8 revisado
                  end do
               end do

!! Make a auxDabij denominator | It's the same Dabij denominator that is used in interspecies part
!!                                         auxDabij(aa,aaa,ii,iii) = Fs%values(ii,ii)+otherFs%values(iii,iii)-Fs%values(aa,aa)-otherFs%values(aaa,aaa)

!!             auxTdNew(a-nop,b,ii,jj) = auxTdNew(a-nop,b,ii,jj)/auxDabij(a,b,ii,jj)
               
               auxTdNew(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj)/(Fs%values(ii,ii)+otherFs%values(jj,jj)-Fs%values(a,a)-otherFs%values(b,b))
          auxTd(a-nop,b-nops,ii,jj) = auxTdNew(a-nop,b-nops,ii,jj)
                                        auxtaus(a-nop,b-nops,ii,jj) = auxTd(a-nop,b-nops,ii,jj)
                                        auxtau(a-nop,b-nops,ii,jj) = auxTd(a-nop,b-nops,ii,jj)

        end do
      end do
    end do
  end do

end do !! Loop

  coupledClusterValue%values(mmm) = auxECCSD

      end do !! Number Of
  end do     !! Species
!
  call Matrix_destructor(auxMatrix)
! CCD_instance%mp2Correction = sum ( CCD_instance%mp2Coupling%values )
  CCD_instance%ccsdTwoParticlesCorrection = auxECCSD
  call Vector_show(coupledClusterValue)

end if

!!*********************************************************************************************************************************
!! End CCSD Calculation different species
!!*********************************************************************************************************************************

   close(wfnUnit)



 end subroutine CCD_iterateIntermediates_DiffSpecies





 
 !**
 ! @ Returns final energy (with CCSD correction)
 !**

! function CCD_getTotalEnergy() result(output)
! implicit none
! real(8) :: output
!
! output = CCD_instance%totalEnergy
!
! end function CCD_getTotalEnergy

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CCD_getTransformedIntegrals()
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
    allocate(CCD_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(CCD_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))

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
        call Matrix_constructor (CCD_instance%twoCenterIntegrals(i), int(numberOfContractions,8), int(numberOfContractions,8) )
        call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8))
        call Matrix_constructor (couplingMatrix,int(numberOfContractions,8), int(numberOfContractions,8))

        hcoreMatrix  = HartreeFockD_instance%HcoreMatrix 

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
                    CCD_instance%twoCenterIntegrals(i)%values(m,n) = &
                        CCD_instance%twoCenterIntegrals(i)%values(m,n) + &
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
!!                      CCD_instance%twoCenterIntegrals(i)%values(m,n) = &
!!                           CCD_instance%twoCenterIntegrals(i)%values(m,n) + &
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
             CCD_instance%twoCenterIntegrals(i)%values(n,m)=&
                  CCD_instance%twoCenterIntegrals(i)%values(m,n)
          end do
       end do
          
       ! print *, "Independent Particle"
       ! call Matrix_show ( CCD_instance%twoCenterIntegrals(i) )

!       write (6,"(T10,A)")"TWO PARTICLES INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
!       print *,""

!       call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer, MolecularSystem_getEigenvectors(specieID), &
!            CCD_instance%fourCenterIntegrals(i,i), specieID, trim(nameOfSpecie) )

        call ReadTransformedIntegrals_readOneSpecies( specieID, CCD_instance%fourCenterIntegrals(i,i)   )

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
!                     CCD_instance%fourCenterIntegrals(i,j), specieID, nameOfSpecie, otherSpecieID, nameOfOtherSpecie )

                 call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         CCD_instance%fourCenterIntegrals(i,j) )

             end if
          end do
       end if
    end do
!    print *,"END INTEGRALS TRANFORMATION:"
    close (wfnUnit)
  call Matrix_destructor (hcoreMatrix)
        call Matrix_destructor (couplingMatrix)

  end subroutine CCD_getTransformedIntegrals

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine CCD_exception( typeMessage, description, debugDescription)
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

  end subroutine CCD_exception

end module CCD_
