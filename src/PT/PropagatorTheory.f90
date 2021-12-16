!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!
!!		J. Romero (jromerof@unal.edu.co)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module PropagatorTheory_
#ifdef intel
	use IFPORT
#endif
	use MolecularSystem_
        use InputCI_
!	use IntegralManager_
!	use GenericInterface_
!        use PInterface_ 
	use Exception_
        use Matrix_
	use Vector_
        use ReadTransformedIntegrals_
	use IndexMap_
  use omp_lib
!	use TransformIntegrals_
!	use TransformIntegrals2_
	implicit NONE

	!>
	!! @brief Implementation of propagator theory
	!!
	!!
	!<

	!< enum PropagatorTheory_correctionFlags {
	integer, parameter :: FIRST_ORDER = 1
	integer, parameter :: SECOND_ORDER = 2
	integer, parameter :: THIRD_ORDER = 3
	!< }

	type, private :: PropagatorTheory

		character(50) :: name
		integer :: orderOfCorrection
		integer :: numberOfSpecies
		integer :: occupationBoundary
		integer :: virtualBoundary

		!! Matrices to store the energy corrections
		type(Matrix) :: energyCorrectionsOfSecondOrder
		type(Matrix),allocatable :: secondOrderCorrections(:)
		type(Matrix),allocatable :: thirdOrderCorrections(:)
                type(Matrix) :: energyCorrections
                type(IMatrix8), allocatable :: xy(:)
                type(IVector8), allocatable :: ioff(:)
                integer(8), allocatable :: ssize2(:)
		logical :: isInstanced
                logical :: externalSCS

	end type PropagatorTheory

	type(PropagatorTheory), private, target :: PropagatorTheory_instance
	
	private :: &
		PropagatorTheory_secondOrderCorrection, &
                ! PropagatorTheory_nonDiagonalSecondOrderCorrection, &
                ! PropagatorTheory_nonDiagonalSecondOrderTDACorrection, &
!		PropagatorTheory_thirdOrderCorrection, &  !! Commented 30th August 2014
!        	PropagatorTheory_thirdOrderCorrection2, & !! Commented 30th August 2014
!		PropagatorTheory_thirdOrderCorrection3, & !! Commented 30th August 2014
!		PropagatorTheory_thirdOrderCorrection4, & !! Commented 30th August 2014
		PropagatorTheory_thirdOrderCorrection5
	
	public :: &
		PropagatorTheory_constructor, &
		PropagatorTheory_destructor, &
		PropagatorTheory_show, &
		PropagatorTheory_run

contains
	  
  !**
  ! Defines the class' constructor
  !
  !**
  subroutine PropagatorTheory_constructor( orderOfCorrection )
    implicit NONE
    integer, intent(in) :: orderOfCorrection
    integer :: maxi
    
    integer :: i
    type(Exception) :: ex
    
    if( .not.PropagatorTheory_instance%isInstanced ) then
       
       PropagatorTheory_instance%orderOfCorrection = orderOfCorrection
       PropagatorTheory_instance%numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
       ! Define the maximal number of orbitals whose corrections will be calculated
       maxi= MolecularSystem_getOcupationNumber( 1 )
       if (PropagatorTheory_instance%numberOfSpecies.gt.1) then
          do i=2, PropagatorTheory_instance%numberOfSpecies
             maxi=max(maxi,MolecularSystem_getOcupationNumber( i ))                        
          end do
       end if
       
       if ( PropagatorTheory_instance%orderOfCorrection == 2 ) then
          
          call Matrix_constructor( PropagatorTheory_instance%energyCorrectionsOfSecondOrder, int(PropagatorTheory_instance%numberOfSpecies,8), &
               4*int(maxi,8))

       end if
       
       if ( PropagatorTheory_instance%orderOfCorrection >= 4 ) then
          
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object PropagatorTheory in the constructor() function" )
          call Exception_setDescription( ex, "This order correction hasn't been implemented" )
          call Exception_show( ex )
          
       end if
       
       PropagatorTheory_instance%isInstanced =.true.

       call PropagatorTheory_setIndexMap
       
    end if
    
  end subroutine PropagatorTheory_constructor

  subroutine PropagatorTheory_setIndexMap( )
    implicit NONE
    integer :: s
    integer(8) :: m, p, q, pq, ssize, ssize2

    allocate ( PropagatorTheory_instance%xy (PropagatorTheory_instance%numberOfSpecies ))
    allocate ( PropagatorTheory_instance%ioff (PropagatorTheory_instance%numberOfSpecies ))
    allocate ( PropagatorTheory_instance%ssize2 (PropagatorTheory_instance%numberOfSpecies ))

    do s = 1, PropagatorTheory_instance%numberOfSpecies 

      ssize = MolecularSystem_getTotalNumberOfContractions( s )
      call Matrix_constructorInteger8( PropagatorTheory_instance%xy(s), ssize, ssize, 0_8 )

      m = 0
      do p = 1, ssize
        do q = p, ssize
          m = m + 1
          PropagatorTheory_instance%xy(s)%values(p,q) = m
          PropagatorTheory_instance%xy(s)%values(q,p) = m
        end do
      end do
    
      ssize2 = ssize * ( ssize + 1 ) / 2
      PropagatorTheory_instance%ssize2(s) = ssize2
      call Vector_constructorInteger8( PropagatorTheory_instance%ioff(s), ssize2, 0_8 )

      PropagatorTheory_instance%ioff(s)%values(1) = 0 
      do pq = 2, ssize2 
        PropagatorTheory_instance%ioff(s)%values(pq) = PropagatorTheory_instance%ioff(s)%values(pq-1) + ssize2 - pq + 1 
      end do
    end do 

  end subroutine PropagatorTheory_setIndexMap

  function PropagatorTheory_IndexMapAA( i, j, k, l, s ) result ( ijkl )
    implicit NONE
    integer :: s
    integer(4) :: i, j, k, l
    integer(8) :: ijkl, ij, kl

    ij = PropagatorTheory_instance%xy(s)%values(i,j)
    kl = PropagatorTheory_instance%xy(s)%values(k,l)

    if ( ij >= kl ) then
      ijkl = PropagatorTheory_instance%ioff(s)%values(kl) + ij
    else
      ijkl = PropagatorTheory_instance%ioff(s)%values(ij) + kl
    end if

  end function PropagatorTheory_IndexMapAA

  function PropagatorTheory_IndexMapAB( i, j, k, l, sa, sb ) result ( ijkl )
    implicit NONE
    integer :: sa, sb
    integer(4) :: i, j, k, l
    integer(8) :: ijkl, ij, kl

    ij = PropagatorTheory_instance%xy(sa)%values(i,j)
    kl = PropagatorTheory_instance%xy(sb)%values(k,l)

    !if ( sa <= sb ) then
    !ijkl = (kl-1) * PropagatorTheory_instance%ssize2(sa) + ij
    !print *, "B1", ijkl
    !else  
    ijkl = (ij-1) * PropagatorTheory_instance%ssize2(sb) + kl
    !print *, "B2", ijkl
    !end if

  end function PropagatorTheory_IndexMapAB
  
  !**
  ! Defines the class' destructor
  !
  !**
  subroutine PropagatorTheory_destructor()
    implicit NONE
    
    PropagatorTheory_instance%isInstanced =.false.
    
  end subroutine PropagatorTheory_destructor
  
  !**
  ! @brief Prints final results of propagation theory
  !**
  subroutine PropagatorTheory_show()
    implicit NONE
    
    integer :: i, j, p, q, m, n, z
    integer :: species1ID
    integer :: species2ID
    integer :: occupationNumber
    real(8) :: orbital
    type(Vector) :: eigenValues
    character(10) :: nameOfSpecies
    
    if ( PropagatorTheory_instance%isInstanced )  then

       if (PropagatorTheory_instance%orderOfCorrection==2) then
          
          print *,""
          print *," POST HARTREE-FOCK CALCULATION"
          print *," PROPAGATOR THEORY:"
          print *,"=============================="
          print *,""
          write (6,"(T10,A50)") "PROPAGATOR FORMALISM FOR SEVERAL FERMIONS SPECIES "
          
          if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
             write (6,"(T10,A23,I5,A23)") "ORDER OF CORRECTION = ",PropagatorTheory_instance%orderOfCorrection," + TRANSITION OPERATOR"
             write (6,"(T10,A40)") "The following articles must be cited:"
             write (6,"(T10,A40)") "-------------------------------------"
             write (6,"(T10,A40)") "J. Chem. Phys. 127, 134106 (2007) "
             write (6,"(T10,A40)") "J. Chem. Phys. 137, 074105 (2012) "
             write (6,"(T10,A40)") "J. Chem. Phys. 138, 194108 (2013) "
             write (6,"(T10,A40)") "CPL coming soon "
             write (6,"(T10,A40)") "-------------------------------------"
          else if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
             write (6,"(T10,A23,I5)") "ORDER OF CORRECTION = ",PropagatorTheory_instance%orderOfCorrection
             write (6,"(T10,A40)") "The following articles must be cited:"
             write (6,"(T10,A40)") "-------------------------------------"
             write (6,"(T10,A40)") "J. Chem. Phys. 137, 074105 (2012) "
             write (6,"(T10,A40)") "J. Chem. Phys. 138, 194108 (2013) "
             write (6,"(T10,A40)") "CPL coming soon "
             write (6,"(T10,A40)") "-------------------------------------"
          else
             write (6,"(T10,A23,I5)") "ORDER OF CORRECTION = ",PropagatorTheory_instance%orderOfCorrection
             write (6,"(T10,A40)") "The following articles must be cited:"
             write (6,"(T10,A40)") "-------------------------------------"
             write (6,"(T10,A40)") "J. Chem. Phys. 137, 074105 (2012) "
             write (6,"(T10,A40)") "J. Chem. Phys. 138, 194108 (2013) "
             write (6,"(T10,A40)") "CPL coming soon "
             write (6,"(T10,A40)") "-------------------------------------"
          end if

          if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE") then
             species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE(1) )
             do z = 1, size(CONTROL_instance%IONIZE_SPECIE )
               if (CONTROL_instance%IONIZE_SPECIE(z) /= "NONE" ) then
               species2ID= MolecularSystem_getSpecieID(CONTROL_instance%IONIZE_SPECIE(z))
               end if
             end do 
          else
             species1ID=1
             species2ID=PropagatorTheory_instance%numberOfSpecies
          end if

          i = 0

          do q= species1ID, species2ID
             
             i = i + 1
             
             n=size(PropagatorTheory_instance%secondOrderCorrections(i)%values,DIM=1)

             nameOfSpecies=trim(MolecularSystem_getNameOfSpecie( q ))

             write (6,"(T10,A8,A10)") "SPECIES: ",nameOfSpecies

             if (nameOfSpecies=="E-ALPHA".or.nameOfSpecies=="E-BETA") then
                

                if ( PropagatorTheory_instance%externalSCS .eqv. .false. ) then

                  write ( 6,'(T10,A90)') "-------------------------------------------------------------------------------------------------"
                  write ( 6,'(T10,A12,A12,A12,A12,A12,A12,A12,A12)') " Orbital ","  KT (eV) "," EP2 (eV) ","  P.S  "," SCS-EP2(eV)"&
                       ,"  P.S  "," SOS-EP2(eV)","  P.S  "
                  write ( 6,'(T10,A90)') "-------------------------------------------------------------------------------------------------"
                  
                  do j=1,n
                     write (*,'(T10,A4,I4,A4,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6)') "    ",&
                          int(PropagatorTheory_instance%secondOrderCorrections(i)%values(j,1)),&
                          "    ",PropagatorTheory_instance%secondOrderCorrections(i)%values(j,2), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,3), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,4), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,5), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,6), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,7), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,8)
                  end do
                  write ( 6,'(T10,A90)') "----------------------------------------------------------------------------------------------"

                else 
                  write ( 6,'(T10,A110)') "---------------------------------------------------------------------------------------------------------------------"
                  write ( 6,'(T10,A12,A12,A12,A12,A12,A12,A12,A12,A12,A12)') " Orbital ","  KT (eV) "," EP2 (eV) ","  P.S  "," SCS-EP2(eV)"&
                       ,"  P.S  "," SOS-EP2(eV)","  P.S  ","*SCS-EP2(eV)","  P.S  "


                  write ( 6,'(T10,A110)') "---------------------------------------------------------------------------------------------------------------------"
                  
                  do j=1,n
                     write (*,'(T10,A4,I4,A4,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6)') "    ",&
                          int(PropagatorTheory_instance%secondOrderCorrections(i)%values(j,1)),&
                          "    ",PropagatorTheory_instance%secondOrderCorrections(i)%values(j,2), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,3), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,4), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,5), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,6), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,7), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,8), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,9), &
                          PropagatorTheory_instance%secondOrderCorrections(i)%values(j,10)

                  end do

                  write ( 6,'(T10,A110)') "---------------------------------------------------------------------------------------------------------------------"
                  write ( 6,'(T10,A12,A11,F12.6,A12,F12.6)') "*SCS-EP2    ","Factor OS: ",  CONTROL_instance%PT_FACTOR_OS, &
                         "Factor SS: ",  CONTROL_instance%PT_FACTOR_SS 
                end if
                
             else
                   
                write ( 6,'(T10,A50)') "-----------------------------------------------------------"
                write ( 6,'(T10,A12,A12,A12,A12)') " Orbital ","  KT (eV) ","  EP2 (eV)","  P.S  "
                write ( 6,'(T10,A50)') "-----------------------------------------------------------"
                
                do j=1,n
                   write (*,'(T10,A4,I4,A4,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6,F12.6)') "    ",&
                        int(PropagatorTheory_instance%secondOrderCorrections(i)%values(j,1)),&
                        "    ",PropagatorTheory_instance%secondOrderCorrections(i)%values(j,2), &
                        PropagatorTheory_instance%secondOrderCorrections(i)%values(j,3), &
                        PropagatorTheory_instance%secondOrderCorrections(i)%values(j,4)
                end do
                write ( 6,'(T10,A50)') "-----------------------------------------------------------"

             end if

          end do

       else if (PropagatorTheory_instance%orderOfCorrection==3) then
          
          print *,""
          print *," POST HARTREE-FOCK CALCULATION"
          print *," PROPAGATOR THEORY:"
          print *,"=============================="
          print *,""
          write (6,"(T10,A50)") "PROPAGATOR FORMALISM FOR SEVERAL FERMIONS SPECIES "
          
          write (6,"(T10,A23,I5,A23)") "ORDER OF CORRECTION = ",PropagatorTheory_instance%orderOfCorrection," PARTIAL AND COMPLETE"
          write (*,"(T10,A40)") "The following articles must be cited:"
          write (*,"(T10,A64)") "---------------------------------------------------------------"
          write (*,"(T10,A60)") "Partial third order: J. Chem. Phys. 104, 1008 (1996)"
          write (*,"(T10,A60)") "LOWDIN implementation: J. Chem. Phys. 137, 074105 (2012)"
          write (*,"(T10,A60)") "LOWDIN implementation: J. Chem. Phys. 138, 194108 (2013)"
          write (*,"(T10,A64)") "---------------------------------------------------------------"

          if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE") then
             species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE(1) )
             do z = 1, size(CONTROL_instance%IONIZE_SPECIE )
               if (CONTROL_instance%IONIZE_SPECIE(z) /= "NONE" ) then
               species2ID= MolecularSystem_getSpecieID(CONTROL_instance%IONIZE_SPECIE(z))
               end if
             end do 
          else
             species1ID=1
             species2ID=PropagatorTheory_instance%numberOfSpecies
          end if

!          i = 0
!
!          do q= species1ID, species2ID
!
!             i = i + 1
!
!             n=size(PropagatorTheory_instance%thirdOrderCorrections(i)%values,DIM=1)
!             write (6,"(T10,A8,A10)")"SPECIE: ",trim(MolecularSystem_getNameOfSpecie( q ))
!             write ( 6,'(T10,A85)') "--------------------------------------------------------------------------------------------"
!             write ( 6,'(T10,A10,A10,A10,A10,A10,A10,A10,A10)') " Orbital ","  KT (eV) ","  EP2 (eV)","  P.S  ","  P3 (eV)"&
!                  ,"  P.S  "," OVGF (eV)","  P.S  "
!             write ( 6,'(T10,A85)') "--------------------------------------------------------------------------------------------"
!
!             do j=1,n
!                write (*,'(T10,A4,I4,A4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)') "    ",&
!                     int(PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,1)),&
!                     "    ",PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,2), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,3), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,4), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,5), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,6), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,7), &
!                     PropagatorTheory_instance%thirdOrderCorrections(i)%values(j,8)
!             end do
!             write ( 6,'(T10,A85)') "-------------------------------------------------------------------------------------------"
!
!          end do

       end if

    end if
    
  end subroutine PropagatorTheory_show
  
  !**
  ! @brief Evaluate numerators and denominators of the self-energy expression and evaluate the corrected 
  ! koopmans energy at the second order level, include modifications for spin scaled second order and
  ! transition operator method
  !**

  subroutine PropagatorTheory_secondOrderCorrection()
    implicit NONE
    
    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
    integer :: i, j, k ! counters for species
    integer :: z, zz
    integer :: m, n, o, p, q, ni, nc, limit, id1, id2, id3 ! auxiliar counters
    integer :: speciesAID, speciesBID, speciesCID
    integer :: species1ID, species2ID
    integer :: electronsID
    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
    integer :: activeOrbitalsOfSpeciesC
    integer(8) :: vectorSize1, vectorSize2, vectorSize3 !!! Sizes for diagrams
    integer(8) :: auxIndex
    integer(4) :: errorNum
    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!    type(TransformIntegrals) :: repulsionTransformer
    type(Matrix),allocatable :: auxMatrix2(:), selfEnergy2hp(:), selfEnergy2ph(:), selfEnergyhp(:)
    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
    real(8), allocatable :: factorOS(:), factorSS(:)
    real(8) :: E2hp, E2ph, Ehp, dE2hp, dE2ph, dEhp, TE2hp, TE2ph, orx, prx, Tprx, Torx
    real(8) :: lastOmega, newOmega, residual, koopmans, selfEnergy, selfEnergyDerivative 
    real(8) :: selfEnergySS, selfEnergyDerivativeSS, selfEnergyOS, selfEnergyDerivativeOS 
    real(8) :: a1, a2, b, c, d, poleStrenght
    logical :: paso1, paso2, paso3, paso4
    ! *******************************************************************************************
    ! Determinate the numerators and denominators of the second Oder propapator 
    character(50) :: wfnFile
    character(50) :: arguments(2)
    integer :: wfnUnit

    ! if ( .not. CONTROL_instance%LOCALIZE_ORBITALS) then
       wfnFile = "lowdin.wfn"
    ! else
    !    wfnFile = "lowdin-subsystemA.wfn"
    ! end if

    wfnUnit = 20

   open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 
   rewind(wfnUnit)

 
    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
    
    print *,"*******************************************************************"
    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"

    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified

    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies))

    !!! Defining for which species the correction will be applied
    
!    if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE") then
!       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE(1) )
!       species2ID= species1ID
!       m=1
!    else
!       species1ID=1
!       species2ID=PropagatorTheory_instance%numberOfSpecies
!       m = species2ID
!    end if

      if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE") then
             species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE(1) )
             zz = 0 
             do z = 1, size(CONTROL_instance%IONIZE_SPECIE )
               if (CONTROL_instance%IONIZE_SPECIE(z) /= "NONE" ) then
               zz = zz + 1
               species2ID= MolecularSystem_getSpecieID(CONTROL_instance%IONIZE_SPECIE(z))
               end if
             end do 
             m = zz 
          else
             species1ID=1
             species2ID=PropagatorTheory_instance%numberOfSpecies
             m = species2ID
          end if


    if (allocated(PropagatorTheory_instance%secondOrderCorrections)) deallocate(PropagatorTheory_instance%secondOrderCorrections)
    allocate(PropagatorTheory_instance%secondOrderCorrections(m))

    ! Factors for spin scaling
    
    if ( CONTROL_instance%PT_FACTOR_SS == 0 .and. CONTROL_instance%PT_FACTOR_OS == 0 ) then
        PropagatorTheory_instance%externalSCS = .false. 

        if (allocated(factorOS)) deallocate (factorOS)
        allocate(factorOS(3))
        if (allocated(factorSS)) deallocate (factorSS)
        allocate(factorSS(3))
    else 
        PropagatorTheory_instance%externalSCS = .true. 

        if (allocated(factorOS)) deallocate (factorOS)
        allocate(factorOS(4))
        if (allocated(factorSS)) deallocate (factorSS)
        allocate(factorSS(4))

        factorOS(4)=CONTROL_instance%PT_FACTOR_OS
        factorSS(4)=CONTROL_instance%PT_FACTOR_SS 

    end if
    ! Regular calculation
    factorOS(1)=1.0_8
    factorSS(1)=1.0_8
    ! SCS
    factorOS(2)=1.2_8
    factorSS(2)=1.0_8/3.0_8
    ! SOS
    factorOS(3)=1.3_8
    factorSS(3)=0.0_8
    
    !!! Start loop for species

    q = 0

    do i = species1ID , species2ID

       q = q + 1

       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )

       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i ) 
       if ( InputCI_Instance(i)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesA = InputCI_Instance(i)%activeOrbitals

       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA

       !!! Defining the number of orbitals !!! Insert a parameter for the else option

       if (CONTROL_instance%PT_TRANSITION_OPERATOR.or.CONTROL_instance%PT_JUST_ONE_ORBITAL) then
          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
          n = 1
       !else if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
          !PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
          !PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
          !n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
       else if ( CONTROL_instance%IONIZE_MO /= 0) then
          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
          n = 1
       else
          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
          if(occupationNumberOfSpeciesA .eq. 0) then
             PropagatorTheory_instance%occupationBoundary = 1
          else
             PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
          end if
          n = 2
       end if

       if (nameOfSpeciesA=="E-ALPHA".or.nameOfSpeciesA=="E-BETA") then
          
          call Matrix_constructor(PropagatorTheory_instance%secondOrderCorrections(q), int(n,8), 10_8, 0.0_8)

       else 

          call Matrix_constructor(PropagatorTheory_instance%secondOrderCorrections(q), int(n,8), 4_8, 0.0_8)

       end if

       ! Occupations

       call Vector_constructor(occupationsOfSpeciesA,occupationNumberOfSpeciesA,1.0_8)

       if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
          
          occupationsOfSpeciesA%values(CONTROL_instance%IONIZE_MO)=CONTROL_instance%MO_FRACTION_OCCUPATION
       
       end if

       ! Storing transformed integrals !!!! We need a more efficient algorithm for this
       
!       call TransformIntegrals_constructor( repulsionTransformer )

        arguments(2) = MolecularSystem_getNameOfSpecie(i)

        arguments(1) = "ORBITALS"

        call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( i ), &
               unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
               output =  eigenValuesOfSpeciesA )     
       
       do p = 1 , PropagatorTheory_instance%numberOfSpecies
          
             activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( p )
             if ( InputCI_Instance(p)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(p)%activeOrbitals

             arguments(2) = trim(MolecularSystem_getNameOfSpecie(p))

             arguments(1) = "ORBITALS"
             call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( p ), &
                     unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                     output =  eigenValuesOfSpeciesB  )    

          if (p==i) then

!             call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p), p, trim(nameOfSpeciesA) )

          !! Read transformed integrals from file
             call ReadTransformedIntegrals_readOneSpecies( p, auxMatrix2(p) )
             
             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( p ) &
                  * MolecularSystem_getCharge( p )
             
          else
             
              !! Read transformed integrals from file
              call ReadTransformedIntegrals_readTwoSpecies( i, p, auxMatrix2(p) )

!             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!                  MolecularSystem_getEigenVectors(i), MolecularSystem_getEigenVectors(p), &
!                  auxMatrix2(p), i, nameOfSpeciesA, p, nameOfSpeciesB )
             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( i ) &
                  * MolecularSystem_getCharge( p )
          end if
          
       end do
       
       !**************************************************************************
       !	Storing of denominators and numerators in the corresponding vectors
       !****

       m =0
       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	

          m=m+1          

          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
          
          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))

          if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
             
             if (allocated(selfEnergyhp)) deallocate(selfEnergyhp)
             allocate(selfEnergyhp(PropagatorTheory_instance%numberOfSpecies))
             
          end if
          
          do j = 1 , PropagatorTheory_instance%numberOfSpecies             

             activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
             if ( InputCI_Instance(j)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(j)%activeOrbitals

             arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))

             arguments(1) = "ORBITALS"
             call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( j ), &
                     unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                     output =  eigenValuesOfSpeciesB  )     


             
             if (j==i) then ! Intraspecies factors
                
                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA

                call Matrix_constructor(selfEnergy2ph(j), 2_8, vectorSize1, 0.0_8)
                
                if (occupationNumberOfSpeciesA>1) then
                   
                   call Matrix_constructor(selfEnergy2hp(j), 6_8, vectorSize2, 0.0_8)
                   !! 2hp(1,2), prx(3,4), orx(5,6) (numerator and denominator)

                end if

                if (CONTROL_instance%PT_TRANSITION_OPERATOR) then

                   vectorSize3 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
                   call Matrix_constructor(selfEnergyhp(j), 2_8, vectorSize3, 0.0_8)

                end if

                id1 = 0
                id2 = 0

                ! factor 2ph
                
                do ia = 1 , occupationNumberOfSpeciesA
                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         
                         !!auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ia, ba, i )
                         auxValue_A= auxMatrix2(j)%values(auxIndex, 1)

                         !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ia, aa, i )
                         auxValue_B= auxMatrix2(j)%values(auxIndex, 1)

                         id1 = id1 + 1
     
                         selfEnergy2ph(j)%values(1,id1) = auxValue_A*(lambdaOfSpeciesA*auxValue_A - auxValue_B)
                         
                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
                              - eigenValuesOfSpeciesA%values(ba)
                         
                      end do
                   end do
                end do
                
                if (occupationNumberOfSpeciesA > 1) then

                   ! factor 2hp

                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ia = 1 , occupationNumberOfSpeciesA
                         do ja = 1 , occupationNumberOfSpeciesA
                            
                            id2 = id2 + 1
                            
                            !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
                            auxIndex = PropagatorTheory_IndexMapAA(pa, ia, ja, aa, i )
                            auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
                            !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, aa, activeOrbitalsOfSpeciesA )
                            auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ia, aa, i )
                            auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
                            
                            selfEnergy2hp(j)%values(1,id2) = occupationsOfSpeciesA%values(ia)*occupationsOfSpeciesA%values(ja)*&
                                 auxValue_A*(lambdaOfSpeciesA*auxValue_A - auxValue_B)

                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
                                 - eigenValuesOfSpeciesA%values(ja) 

                            if ( pa == ia ) then
                                 selfEnergy2hp(j)%values(5,id2) = selfEnergy2hp(j)%values(1,id2) 
                                 selfEnergy2hp(j)%values(6,id2) = selfEnergy2hp(j)%values(2,id2) 
                            else 
                                 selfEnergy2hp(j)%values(3,id2) = selfEnergy2hp(j)%values(1,id2) 
                                 selfEnergy2hp(j)%values(4,id2) = selfEnergy2hp(j)%values(2,id2) 
                            end if

                         end do
                      end do
                   end do

                end if

                if (CONTROL_instance%PT_TRANSITION_OPERATOR) then

                   id3 = 0

                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ia = 1 , occupationNumberOfSpeciesA
                         
                         id3 = id3 + 1
                         
                         !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, pa, ia, aa, i )
                         auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
                         !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ia, pa, i )
                         auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
                          
                         selfEnergyhp(j)%values(1,id3) = ( 1 - occupationsOfSpeciesA%values(pa)) * &
                            auxValue_A*(lambdaOfSpeciesA*auxValue_A - auxValue_B)
                         
                         selfEnergyhp(j)%values(2,id3) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(pa) &
                              - eigenValuesOfSpeciesA%values(aa)
                                                  
                      end do
                   end do
                   
                end if
                
             else ! interspecies

                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
                chargeOfSpeciesB = MolecularSystem_getCharge( j )
!                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
                occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
                lambdaOfSpeciesB = MolecularSystem_getLambda( j )
                virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB

                call Vector_constructor(occupationsOfSpeciesB,occupationNumberOfSpeciesB,1.0_8)

                vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
                vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB

                call Matrix_constructor(selfEnergy2ph(j), 2_8, vectorSize1, 0.0_8)
                call Matrix_constructor(selfEnergy2hp(j), 6_8, vectorSize2, 0.0_8)

                if (CONTROL_instance%PT_TRANSITION_OPERATOR) then

                   vectorSize3 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesB
                   call Matrix_constructor(selfEnergyhp(j), 2_8, vectorSize3, 0.0_8)

                end if

                id1 = 0
                id2 = 0

                ! diagram A
                
                do ib = 1 , occupationNumberOfSpeciesB
                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                         
                        !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                        !print *, "A", auxIndex
                         auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ib, ab, i, j )
                         auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
                         
                         id1 = id1 + 1
                         
                         selfEnergy2ph(j)%values(1,id1) = lambdaOfSpeciesA*lambdaOfSpeciesB*((auxValue_A)**2.0_8)
                         
                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
                              - eigenValuesOfSpeciesB%values(ab)
                         
                      end do
                   end do
                end do
                
                ! diagram B
                
                do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                   do ia = 1 , occupationNumberOfSpeciesA
                      do ib = 1 , occupationNumberOfSpeciesB
                         
                         id2 = id2 + 1
                         
                        !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                         auxIndex = PropagatorTheory_IndexMapAB(pa, ia, ib, ab, i, j )
                         auxValue_A = auxMatrix2(j)%values(auxIndex, 1)

                         selfEnergy2hp(j)%values(1,id2) = occupationsOfSpeciesA%values(ia)*&
                              lambdaOfSpeciesA*lambdaOfSpeciesB*((auxValue_A)**2.0_8)
                         
                         selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
                              - eigenValuesOfSpeciesB%values(ib)

                         if ( pa == ia ) then
                              selfEnergy2hp(j)%values(5,id2) = selfEnergy2hp(j)%values(1,id2) 
                              selfEnergy2hp(j)%values(6,id2) = selfEnergy2hp(j)%values(2,id2) 
                         else 
                              selfEnergy2hp(j)%values(3,id2) = selfEnergy2hp(j)%values(1,id2) 
                              selfEnergy2hp(j)%values(4,id2) = selfEnergy2hp(j)%values(2,id2) 
                         end if


                      end do
                   end do
                end do

             if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
                
                id3 = 0
                
                do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                   do ib = 1 , occupationNumberOfSpeciesB
                      
                      id3 = id3 + 1
                      
                      !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                      auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ib, ab, i, j )
                      auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
                      
                      selfEnergyhp(j)%values(1,id3) = ( 1- occupationsOfSpeciesA%values(pa))*&
                        lambdaOfSpeciesA*lambdaOfSpeciesB*((auxValue_A)**2.0_8)
                      selfEnergyhp(j)%values(2,id3) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(pa) &
                           - eigenValuesOfSpeciesB%values(ab)
                      
                   end do
                end do
                
              end if
            end if

          end do

          ! Koopmans
          koopmans = eigenValuesOfSpeciesA%values(pa)
          PropagatorTheory_instance%secondOrderCorrections(q)%values(m,1)=real(pa,8)
          PropagatorTheory_instance%secondOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans

          print *,"----------------------------------------------------------------"
          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%secondOrderCorrections(q)%values(m,1)),&
               " of species: ",nameOfSpeciesA
          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%secondOrderCorrections(q)%values(m,2)

          ! Selecting value of o, for spin-component-scaled calculations

          if (nameOfSpeciesA=="E-ALPHA".or.nameOfSpeciesA=="E-BETA") then
             o = size(factorSS)
          else
             o = 1
          end if
          
          ! Second order calculations with different spin-component-scaling factors
          do n = 1, o 
             print *, "   SCS :" , n  
             ! Initial guess
             newOmega = koopmans
             lastOmega = 0.0_8
             
             ni = 0
             limit = 50
             residual = 1.0_8
             
             ! Calculation of second order pole
             
             do while ((residual>0.0001_8).or.(limit.lt.ni))
                
                ni = ni + 1
                
                lastOmega = newOmega
                selfEnergy = lastOmega - koopmans
                selfEnergyDerivative = 1.0_8

                do j = 1 , PropagatorTheory_instance%numberOfSpecies             

                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )                   

                   E2hp = 0.0_8
                   E2ph= 0.0_8
                   Ehp = 0.0_8
                   dE2hp = 0.0_8
                   dE2ph= 0.0_8
                   dEhp = 0.0_8

                   do id1 = 1, size(selfEnergy2ph(j)%values,DIM=2)
                      
                      b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                      
                      E2ph = E2ph + selfEnergy2ph(j)%values(1,id1)/b
                      dE2ph = dE2ph + selfEnergy2ph(j)%values(1,id1)/(b**2.0_8)
                      
                   end do


                   if (occupationNumberOfSpeciesA==1.and.i==j) goto 20

                   do id2 = 1, size(selfEnergy2hp(j)%values,DIM=2)
                      
                      b = selfEnergy2hp(j)%values(2,id2) + lastOmega
                      
                      E2hp = E2hp + selfEnergy2hp(j)%values(1,id2)/b
                      dE2hp = dE2hp + selfEnergy2hp(j)%values(1,id2)/(b**2.0_8)
                      
                   end do

                   !!print *, "Specie ij", i,j, "E2hp", E2hp, "E2ph", E2ph
                   
20                 continue
                   
                   if (CONTROL_instance%PT_TRANSITION_OPERATOR) then

                      do id3 = 1, size(selfEnergyhp(j)%values,DIM=2)
                         
                         b = selfEnergyhp(j)%values(2,id3) + lastOmega
                         
                         Ehp = Ehp + selfEnergyhp(j)%values(1,id3)/b
                         dEhp = dEhp + selfEnergyhp(j)%values(1,id3)/(b**2.0_8)
                         
                      end do
                      
                   end if

                   paso1=(nameOfSpeciesA=="E-ALPHA".or.nameOfSpeciesA=="E-BETA")
                   paso2=(j==i)
                   paso3=(nameOfSpeciesA=="E-ALPHA".and.nameOfSpeciesB=="E-BETA").or.(nameOfSpeciesA=="E-BETA".and.nameOfSpeciesB=="E-ALPHA")

                   if (paso1.and.paso2) then

                      selfEnergy = selfEnergy - factorSS(n)*( E2ph + E2hp + Ehp )
                      
                      selfEnergyDerivative = selfEnergyDerivative + factorSS(n)*( dE2ph + dE2hp + dEhp )
                      
                   else if (paso3) then

                      selfEnergy = selfEnergy - factorOS(n)*( E2ph + E2hp + Ehp )

                      selfEnergyDerivative = selfEnergyDerivative + factorOS(n)*( dE2ph + dE2hp + dEhp )
                      
                   else

                      selfEnergy = selfEnergy - ( E2ph + E2hp + Ehp )

                      selfEnergyDerivative = selfEnergyDerivative + ( dE2ph + dE2hp + dEhp )

                   end if

                end do
                
                ! NR iteration
                newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
                
                residual = abs(newOmega-lastOmega)
                
                write (*,"(T5,A11,I2,A10,F8.5,A10,F8.5)") "Iteration: ",ni," newOmega: ",newOmega," residual: ",residual
                
             end do ! end while
             
             poleStrenght = 1.0_8/selfEnergyDerivative
              
             !! P2 decomposition
             print *, ""
             write (*,"(T5,A43,I2,A13,A12)") "P2 decomposition (in eV) for spin-orbital:", &
                   int(PropagatorTheory_instance%secondOrderCorrections(q)%values(m,1)),&
                   " of species a: ",nameOfSpeciesA

             write (*, "(T6,A50)") "--------------------------------------------------"
             write (*, "(T6,A14,A11,A11,A12,A13)"),"Species b     ", "    PRX    ", "    ORX    ", "E_2ph (PRM) ", "\Sigma_{ab}^2"
             write (*, "(T6,A50)") "--------------------------------------------------"
             TE2hp = 0.0_8
             TE2ph = 0.0_8
             Tprx = 0.0_8
             Torx = 0.0_8

             do j = 1 , PropagatorTheory_instance%numberOfSpecies             

                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )                   

                E2hp = 0.0_8
                E2ph= 0.0_8
                Ehp = 0.0_8
                prx = 0.0_8
                orx = 0.0_8

                do id1 = 1, size(selfEnergy2ph(j)%values,DIM=2)
                   
                   b = selfEnergy2ph(j)%values(2,id1) + newOmega
                   E2ph = E2ph + selfEnergy2ph(j)%values(1,id1)/b
                   
                end do


                if (occupationNumberOfSpeciesA==1.and.i==j) goto 30

                do id2 = 1, size(selfEnergy2hp(j)%values,DIM=2)

                   b = selfEnergy2hp(j)%values(2,id2) + newOmega
                   E2hp = E2hp + selfEnergy2hp(j)%values(1,id2)/b
                   
                   b = selfEnergy2hp(j)%values(4,id2) + newOmega
                   prx = prx +  selfEnergy2hp(j)%values(3,id2)/b

                   b = selfEnergy2hp(j)%values(6,id2) + newOmega
                   orx = orx +  selfEnergy2hp(j)%values(5,id2)/b

                end do

                if ((nameOfSpeciesA=="E-ALPHA".and.nameOfSpeciesB=="E-BETA") &
                        .or.(nameOfSpeciesA=="E-BETA".and.nameOfSpeciesB=="E-ALPHA")) then
                        E2ph = factorOS(n) * E2ph
                        E2hp = factorOS(n) * E2hp
                        prx = factorOS(n) * prx
                        orx = factorOS(n) * orx
                end if
                if ((nameOfSpeciesA=="E-ALPHA".and.nameOfSpeciesB=="E-ALPHA") &
                        .or.(nameOfSpeciesA=="E-BETA".and.nameOfSpeciesB=="E-BETA")) then
                        E2ph = factorSS(n) * E2ph
                        E2hp = factorSS(n) * E2hp
                        prx = factorSS(n) * prx
                        orx = factorSS(n) * orx
                end if

                
30              continue
                
                if (CONTROL_instance%PT_TRANSITION_OPERATOR) then

                   do id3 = 1, size(selfEnergyhp(j)%values,DIM=2)
                      
                      b = selfEnergyhp(j)%values(2,id3) + newOmega
                      
                      Ehp = Ehp + selfEnergyhp(j)%values(1,id3)/b
                      
                   end do
                   
                end if

                write (*,"(T6,A10,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5)"), nameOfSpeciesB, prx*27.211396_8, orx*27.211396_8, &
                        E2ph*27.211396_8,(E2hp+E2ph)*27.211396_8

                TE2hp = TE2hp + E2hp
                TE2ph = TE2ph + E2ph
                Torx = Torx + orx
                Tprx = Tprx + prx

             end do ! end do species j

             write (*, "(T6,A50)") "--------------------------------------------------"

             !! Total 
             write (*,"(T6,A10,2X,F10.5,2X,F10.5,2X,F10.5,2X,F10.5)"), "Sum for b " , Tprx*27.211396_8,Torx*27.211396_8, &
                         TE2ph*27.211396_8,(TE2hp+TE2ph)*27.211396_8

             write (*, "(T6,A50)") "--------------------------------------------------"
             write (*, *) ""

             ! Storing corrections
             
             PropagatorTheory_instance%secondOrderCorrections(q)%values(m,2*n+1)=27.211396_8 * newOmega
             PropagatorTheory_instance%secondOrderCorrections(q)%values(m,2*n+2)=poleStrenght
             
             write (*,"(T5,A10,F8.5,A10,F8.5)") " FactorOS: ",factorOS(n)," FactorSS: ",factorSS(n)
             write (*,"(T5,A30,F8.4,A7,I2,A12)") " Optimized second order pole: ",&
                  PropagatorTheory_instance%secondOrderCorrections(q)%values(m,2*n+1),&
                  " after ",ni," iterations."
             write (*,"(T5,A17,F8.4,A15,F7.4)") "Correction(eV): ",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
             write (*, *) ""

          end do          

       end do

       ! call Matrix_destructor(auxMatrix2(:))          
!       call TransformIntegrals_destructor( repulsionTransformer )
       
    end do
       
    !!
    !!************************************************************************************************
    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
    print *,"***************************************************************"
    close(wfnUnit)

  end subroutine PropagatorTheory_secondOrderCorrection  

!
!  !**
!  ! @brief Evaluate partial third order, EP3 and OVGF poles, Correction to OVGF for considering a correction
!  ! for each kind of correlation  ! v2
!  !**

  subroutine PropagatorTheory_thirdOrderCorrection5()
    implicit NONE
    
    type(Exception) :: ex
    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
    integer :: i, j, k ! counters for species
    integer :: m, n, o, p, q, r, s, ni, nc, limit, id1, id2 ! auxiliar counters
    integer :: speciesAID, speciesBID, speciesCID
    integer :: species1ID, species2ID
    integer :: electronsID
    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
    integer :: activeOrbitalsOfSpeciesC
    integer(8) :: vectorSize1, vectorSize2 !!! Sizes for diagrams
    integer(8) :: auxIndex, electrons(2)
    integer(4) :: errorNum
    character(10) :: thirdOrderMethods(6)
    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!    type(TransformIntegrals) :: repulsionTransformer
    type(Matrix),allocatable :: auxMatrix2(:,:), selfEnergy2hp(:), selfEnergy2ph(:)
    type(Matrix),allocatable :: secondOrderDensities(:)
    type(Matrix) :: diagram_A, diagram_B, auxMatrix, auxMatrix3
    type(Matrix) :: partialMO1, partialMO2    
    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
    real(8) :: valueOfW, valueOfU, valueOfdU, sub2, subW, subU, subd2, subdW, subdU
    real(8) :: lastOmega, newOmega, residual, threshold, selfEnergy, selfEnergyDerivative, koopmans 
    real(8) :: a1, a2, b, c, d, poleStrenght, partialValue, partialValue2, initialValue
    real(8) :: fW, fI, thirdOrderResults(2,6)
    real(8) :: value1, value2, value3, value4
    real(8),allocatable :: s2hp(:), s2ph(:), W2hp(:), W2ph(:), U2hp(:), U2ph(:), constantSelfEnergy(:,:), factors(:,:,:)
    real(8),allocatable :: W22hp(:), W22ph(:), U22hp(:), U22ph(:), factors2(:,:,:)
    logical :: paso1, paso2, paso3
    character(50) :: wfnFile
    character(50) :: arguments(2)
    integer :: wfnUnit
    integer :: oo,oarray(6), ooarray(6), maxoo

    ! if ( .not. CONTROL_instance%LOCALIZE_ORBITALS) then
       wfnFile = "lowdin.wfn"
    ! else
    !    wfnFile = "lowdin-subsystemA.wfn"
    ! end if
    wfnUnit = 20

    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted") 
    rewind(wfnUnit)

    ! *******************************************************************************************
    ! Determinate the numerators and denominators of the second Oder propapator 
    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
    
    print *,"*******************************************************************"
    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"
    print *, ""

    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified

    !!! Determining who are electrons

    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies,PropagatorTheory_instance%numberOfSpecies))

    if (allocated(secondOrderDensities)) deallocate(secondOrderDensities)
    allocate(secondOrderDensities(PropagatorTheory_instance%numberOfSpecies))

    ! Containers for different types of corrections

    if (allocated(s2hp)) deallocate(s2hp)
    allocate(s2hp(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(s2ph)) deallocate(s2ph)
    allocate(s2ph(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(W2hp)) deallocate(W2hp)
    allocate(W2hp(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(W2ph)) deallocate(W2ph)
    allocate(W2ph(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(U2hp)) deallocate(U2hp)
    allocate(U2hp(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(U2ph)) deallocate(U2ph)
    allocate(U2ph(PropagatorTheory_instance%numberOfSpecies))

    ! terms for saving a-a-b terms

    if (allocated(W22hp)) deallocate(W22hp)
    allocate(W22hp(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(W22ph)) deallocate(W22ph)
    allocate(W22ph(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(U22hp)) deallocate(U22hp)
    allocate(U22hp(PropagatorTheory_instance%numberOfSpecies))

    if (allocated(U22ph)) deallocate(U22ph)
    allocate(U22ph(PropagatorTheory_instance%numberOfSpecies))

    !!!

    if (allocated(constantSelfEnergy)) deallocate(constantSelfEnergy)
    allocate(constantSelfEnergy(PropagatorTheory_instance%numberOfSpecies,PropagatorTheory_instance%numberOfSpecies))

    if (allocated(factors)) deallocate(factors)
    allocate(factors(PropagatorTheory_instance%numberOfSpecies,3,6))

    ! new terms
    
    if (allocated(factors2)) deallocate(factors2)
    allocate(factors2(PropagatorTheory_instance%numberOfSpecies,3,6))

    ! Defining for which species the correction will be applied
    
    if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE") then
       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE(1) )
       species2ID= species1ID
       m=1
    else
       species1ID=1
       species2ID=PropagatorTheory_instance%numberOfSpecies
       m = species2ID
    end if

    if (allocated(PropagatorTheory_instance%thirdOrderCorrections)) deallocate(PropagatorTheory_instance%thirdOrderCorrections)
    allocate(PropagatorTheory_instance%thirdOrderCorrections(m))

    ooarray = 0

    do oo = 1, 7

      select case ( CONTROL_instance%PT_P3_METHOD(oo) ) 

      case ("ALL")
        ooarray(1) = 1
        ooarray(2) = 2
        ooarray(3) = 3
        ooarray(4) = 4
        ooarray(5) = 5
        ooarray(6) = 6
      case ( "P3" )
        ooarray(1) = 1
      case ( "EP3" )
        ooarray(2) = 2
      case( "OVGF-A" )
        ooarray(2) = 2
        ooarray(3) = 3
      case( "OVGF-B" )
        ooarray(2) = 2
        ooarray(4) = 4
      case ( "OVGF-C" )
        ooarray(2) = 2
        ooarray(5) = 5
      case ( "OVGF" )
        ooarray(2) = 2
        ooarray(3) = 3
        ooarray(4) = 4
        ooarray(5) = 5
      case ( "REN-P3" )
        ooarray(1) = 1
        ooarray(6) = 6
      case ( "TROLOLO" )
        print *, "Trololo"
        print *, "Eduard Khil"
        print *, "https://www.youtube.com/watch?v=oavMtUWDBTM"
        print *, "Lyrics:"
        print *, "Ahhhh, ya ya yaaah,"
        print *, "ya ya yah yah ya yaaah."
        print *, "Oh oh oh oh oooh, oh ya yah,"
        print *, "ya ya yah yah ya yah."
        print *, ""
        print *, "Ye ye ye ye ye, ye ye yeh, ye ye yeh."
        print *, "Oh oh oh oh oooh."
        print *, "Ye ye ye ye ye, ye ye yeh, ye ye yeh."
        print *, "Oh oh oh oh oooh, lololol."
        print *, "Oh oh oooh oooh, la lah."
        print *, "Na na na na nah na na nah na na nah na na nah na na nah."
        print *, "Na na na na nan na na nan, na na nah,"
        print *, "na na na na nah."
        print *, ""
        print *, "Na na na na naaaaah, na na naaaah..."
        print *, "Na na nah nah na na."
        print *, "Lololololoooool,"
        print *, "la la lah."
        print *, "La la lah lah la lah."
        print *, ""
        print *, "Oh oh oh oh oh, oh oh oh, oh oh oh."
        print *, "Oh oh oh oh oh."
        print *, "Oh oh oh oh oh, oh oh oh, oh oh oh,"
        print *, "Lololololol!"
        print *, ""
        print *, "Ah-eeeeeee,"
        print *, "ee-ee-eeeh!"
        print *, "La la lah lah la lah."
        print *, "Oh oh oh oh oooh,"
        print *, "bop a-da da da dah da da dah."
        print *, "Da da dah dah da dah."
        print *, ""
        print *, "Lolololo lol, lololol, lololol."
        print *, "La la la la lah."
        print *, "Trololololol, lololol, lololol,"
        print *, "Oh ha ha ha oh!"
        print *, "Oh ha ha ha oh!"
        print *, "Oh ha ha ha oh!"
        print *, "Oh ha ha ha oh!"
        print *, ""
        print *, "Lolololololol,"
        print *, "lolololololol,"
        print *, "lolololololol,"
        print *, "lolololol!"
        print *, ""
        print *, "Laah la la lah,"
        print *, "la la lah lah la lah."
        print *, "Lololololol, la la lah,"
        print *, "la la lah lah la lah."
        print *, ""
        print *, "Lololololol, lololol, lololol,"
        print *, "oh oh oh oh oh."
        print *, "Lololololol, lololol, lololol,"
        print *, "oh oh oh oh oooooooh!"
        print *, ""

      case ( "NONE" )
      case default 

          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object PropagatorTheory in PropagatorTheory_thirdOrderCorrection5 function" )
          call Exception_setDescription( ex, "This correction hasn't been implemented: "//trim(CONTROL_instance%PT_P3_METHOD(oo))) 
          call Exception_show( ex )     

      end select 
    end do

    oarray = 0
    maxoo = 0
    o = 0
    do oo = 1, 6
      if (ooarray(oo) /= 0 ) then
        o = o + 1
        oarray(o) = ooarray(oo)
      end if
    end do
    maxoo = o


    ! Storing transformed integrals !!!! We need a more efficient algorithm to do this

    r = 0 ! counter for locating who are electrons

    electrons(1)=0
    electrons(2)=0

!    call TransformIntegrals_constructor( repulsionTransformer )
    
    do p = 1 , PropagatorTheory_instance%numberOfSpecies

       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( p ) )

       if (nameOfSpeciesA=="E-ALPHA".or.nameOfSpeciesA=="E-BETA") then

          r = r+1

          electrons(r)=p

       end if

       do n = p , PropagatorTheory_instance%numberOfSpecies
          
          if (n==p) then
             
!JC            call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p,p), p, trim(nameOfSpeciesA) )

          !! Read transformed integrals from file
             call ReadTransformedIntegrals_readOneSpecies( p, auxMatrix2(p,p) )
             
             auxMatrix2(p,p)%values = auxMatrix2(p,p)%values * MolecularSystem_getCharge( p ) &
                  * MolecularSystem_getCharge( p )
             
          else
             
             nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( n ) )
             
!JC             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!                  MolecularSystem_getEigenVectors(p), MolecularSystem_getEigenVectors(n), &
!                  auxMatrix2(p,n), p, nameOfSpeciesA, n, nameOfSpeciesB )

              !! Read transformed integrals from file
              call ReadTransformedIntegrals_readTwoSpecies( p, n, auxMatrix2(p,n) )
             
             auxMatrix2(p,n)%values = auxMatrix2(p,n)%values * MolecularSystem_getCharge( n ) &
                  * MolecularSystem_getCharge( p )
             
          end if
          
       end do

    end do

    ! Start loop for species
    
    q = 0
    
    do i = species1ID , species2ID
       
       q = q + 1
       
       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!JC       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )
       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i ) 
       if ( InputCI_Instance(i)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesA = InputCI_Instance(i)%activeOrbitals
       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA

       arguments(2) = trim(MolecularSystem_getNameOfSpecie(i))

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( i ), &
                 unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                 output =  eigenValuesOfSpeciesA  )    

       ! paso
       
       paso1=(nameOfSpeciesA=="E-ALPHA".or.nameOfSpeciesA=="E-BETA")

       ! Defining the number of orbitals !!! Insert a parameter for the else option
       
       if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
          n = 1
       else if (CONTROL_instance%IONIZE_SPECIE(1) /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
          PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
          n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
       else
          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
          if(occupationNumberOfSpeciesA .eq. 0) then
             PropagatorTheory_instance%occupationBoundary = 1
          else
             PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
          end if
          n = 2
       end if

       call Matrix_constructor(PropagatorTheory_instance%thirdOrderCorrections(q), int(n,8), 8_8, 0.0_8)

       !**************************************************************************
       !	Storing of denominators and numerators in the corresponding vectors
       !****
       
       m =0
       
       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	

          m=m+1          
          
          ! calculation of constant self energy
          if ( ooarray(2) == 2 ) then !! EP3 based methods
          
            constantSelfEnergy(:,:) = 0.0_8
            
            do p = 1 , PropagatorTheory_instance%numberOfSpecies
               
               if (p==i) then
  
  !JC               print *,"entro al if"                
                  ! alpha-alpha-alpha
  
                  do ia = 1 , occupationNumberOfSpeciesA
                     do ja = 1 , occupationNumberOfSpeciesA
                        
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, ia, ja, i )
                        auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
                        !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ia, pa, i )
                        auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8
                        
                        do ka = 1 , occupationNumberOfSpeciesA
                           do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                              do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ka, ba, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, aa, ka, ba, i )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ka, aa, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ba, ka, aa, i )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ja, aa, ka, ba, i )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ja, ba, ka, aa, i )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue &
                                      - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
                                      +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))&
                                      *( eigenValuesOfSpeciesA%values(ja)&
                                      + eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,i) = constantSelfEnergy(i,i) + partialValue*(auxValue_E-auxValue_F)                       
  
                     end do
                  end do
                  
                  do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                        
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, aa, ba, i )
                        auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
                        !auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, ba, aa, pa, i )

                        auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8                      
                        
                        do ia = 1 , occupationNumberOfSpeciesA
                           do ja = 1 , occupationNumberOfSpeciesA
                              do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ja, ca, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, aa, ja, ca, i )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, aa, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ca, ja, aa, i )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ba, ja, ca, i )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ca, ja, ba, i )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue &
                                      + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
                                      +eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca))&
                                      *( eigenValuesOfSpeciesA%values(ia)&
                                      + eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(ba)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,i) = constantSelfEnergy(i,i) + partialValue*(auxValue_E-auxValue_F) 
                        
                     end do
                  end do
                  
                  do ia = 1 , occupationNumberOfSpeciesA
                     do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                        
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, ia, aa, i )
                        auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
                        !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ia, pa, i )

                        auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8
                        
                        do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                           do ja = 1 , occupationNumberOfSpeciesA
                              do ka = 1 , occupationNumberOfSpeciesA
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ja, ba, ka, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ja, ba, ka, i )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ka, ba, ja, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ka, ba, ja, i )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ja, aa, ka, ba, i )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ja, ba, ka, aa, i )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
                                      +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))
                                 
                              end do
                           end do
                        end do
                        
                        do ja = 1 , occupationNumberOfSpeciesA
                           do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                              do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ba, aa, ca, ja, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ba, aa, ca, ja, i )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, aa, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ba, ja, ca, aa, i )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ba, ja, ca, i )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAA(ia, ca, ja, ba, i )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
                                      +eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,i) = constantSelfEnergy(i,i) &
                             + (auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
                        
                     end do
                  end do
  
  !JC               print *,"constant sigma a-a-a:", constantSelfEnergy(i,i)
  
               else
  
  !JC               print *,"entro al else"                
  
                  nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( p ) )
                  chargeOfSpeciesB = MolecularSystem_getCharge( p )
  !JC                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( p )
                  occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( p )
                  activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( p )
                  if ( InputCI_Instance(p)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(p)%activeOrbitals
                  lambdaOfSpeciesB = MolecularSystem_getLambda( p )
                  virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
  
                  arguments(2) = trim(MolecularSystem_getNameOfSpecie(p))
                  arguments(1) = "ORBITALS"
                  call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( p ), &
                       unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                       output =  eigenValuesOfSpeciesB  )     
  
                  ! alpha-beta-beta 
  
                  do ib = 1 , occupationNumberOfSpeciesB
                     do jb = 1 , occupationNumberOfSpeciesB
                        
                        if (p>i) then
                           
                           !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                           auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ib, jb, i, p )

                           auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
                           
                        else
  
                           !auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                           auxIndex = PropagatorTheory_IndexMapAB(ib, jb, pa, pa, p, i )

                           auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
  
                        end if
  
                        partialValue = 0.0_8
                        
                        do kb = 1 , occupationNumberOfSpeciesB
                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                              do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, ab, kb, bb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, ab, kb, bb, p )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, bb, kb, ab, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, bb, kb, ab, p )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(jb, ab, kb, bb, p )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(jb, bb, kb, ab, p )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue &
                                      - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
                                      +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))&
                                      *( eigenValuesOfSpeciesB%values(jb)&
                                      + eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(p,p) = constantSelfEnergy(p,p) + partialValue*auxValue_E
                        
                     end do
                  end do
                  
                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                     do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
  
                        if (p>i) then
  
                           !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                           auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ab, bb, i, p)

                           auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
                        
                        else
  
                           !auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                           auxIndex = PropagatorTheory_IndexMapAB(ab, bb, pa, pa, p, i)

                           auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
  
                        end if
                        partialValue = 0.0_8                      
                        
                        do ib = 1 , occupationNumberOfSpeciesB
                           do jb = 1 , occupationNumberOfSpeciesB
                              do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, ab, jb, cb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, ab, jb, cb, p )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, ab, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, cb, jb, ab, p )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, bb, jb, cb, p )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, cb, jb, bb, p )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue &
                                      + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
                                      +eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(cb))&
                                      *( eigenValuesOfSpeciesB%values(ib)&
                                      + eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(cb) - eigenValuesOfSpeciesB%values(bb)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(p,p) = constantSelfEnergy(p,p) + partialValue*auxValue_E 
                        
                     end do
                  end do
                  
                  do ib = 1 , occupationNumberOfSpeciesB
                     do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                        
                        if (p>i) then
                           
                           !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                           auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ib, ab, i, p  )

                           auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
                        
                        else
                           
                           !auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                           auxIndex = PropagatorTheory_IndexMapAB(ib, ab, pa, pa, p, i )

                           auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
  
                        end if
  
                        partialValue = 0.0_8
                        
                        do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                           do jb = 1 , occupationNumberOfSpeciesB
                              do kb = 1 , occupationNumberOfSpeciesB
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, kb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, jb, bb, kb, p )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, kb, bb, jb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, kb, bb, jb, p )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(jb, ab, kb, bb, p )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(jb, bb, kb, ab, p )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
                                      +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))
                                 
                              end do
                           end do
                        end do
                        
                        do jb = 1 , occupationNumberOfSpeciesB
                           do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                              do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(bb, ab, cb, jb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(bb, ab, cb, jb, p )
                                 auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(bb, jb, cb, ab, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(bb, jb, cb, ab, p )
                                 auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, bb, jb, cb, p )
                                 auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAA(ib, cb, jb, bb, p )

                                 auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
                                 
                                 partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
                                      +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesB%values(cb))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(p,p) = constantSelfEnergy(p,p) &
                             + auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
                        
                     end do
                  end do
  
  !JC                print *,"constant sigma after a-b-b:", constantSelfEnergy(p,p)
  
                  ! alpha-alpha-beta
  
                  do ia = 1 , occupationNumberOfSpeciesA
                     do ja = 1 , occupationNumberOfSpeciesA
  
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, ia, ja, i )
                        auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
                        !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ia, pa, i )

                        auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8
  
                        do ib = 1 , occupationNumberOfSpeciesB
                           do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                              do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
  
                                 if (p>i) then
                                    
                                    !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ib, ab, i, p )
                                    auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ja, aa, ib, ab, i, p )

                                    auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
  
                                 else
  
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, aa, p, i )
                                    auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ja, aa, p, i )

                                    auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
  
                                 end if
  
                                 partialValue = partialValue &
                                      - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
                                      +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
                                      *( eigenValuesOfSpeciesA%values(ja)&
                                      + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,p) = constantSelfEnergy(i,p) + partialValue*(auxValue_E-auxValue_F) 
                        
                     end do
                  end do
                  
                  do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                        
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, aa, ba, i )
                        auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
                       ! auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, ba, aa, pa, i )

                        auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8
                        
                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                           do ia = 1 , occupationNumberOfSpeciesA
                              do ib = 1 , occupationNumberOfSpeciesB
  
                                 if (p>i) then                               
                                    
                                    !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ib, ab, i, p )
                                    auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ia, ba, ib, ab, i, p )

                                    auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
  
                                 else
  
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, aa, p, i )
                                    auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, ba, p, i )

                                    auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
  
                                 end if
  
                                 partialValue = partialValue &
                                      + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
                                      + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
                                      *( eigenValuesOfSpeciesA%values(ia)&
                                      + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab)))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,p) = constantSelfEnergy(i,p) + partialValue*(auxValue_E-auxValue_F) 
                        
                     end do
                  end do
                  
                  do ia = 1 , occupationNumberOfSpeciesA
                     do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                        
                        !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, pa, ia, aa, i )
                        auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
                        !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
                        auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ia, pa, i )

                        auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
                        
                        partialValue = 0.0_8
                        
                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                           do ja = 1 , occupationNumberOfSpeciesA
                              do ib =  1 , occupationNumberOfSpeciesB
  
                                 if (p>i) then                                                              
                                    
                                    !auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ia, ja, ib, ab, i, p )
                                    auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ja, aa, ib, ab, i, p )

                                    auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
  
                                 else
  
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ja, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, ja, p, i )
                                    auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ja, aa, p, i )

                                    auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
  
                                 end if
  
                                 partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ja)&
                                      +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))
                                 
                              end do
                           end do
                        end do
                        
                        do ib = 1 , occupationNumberOfSpeciesB
                           do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                              do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                 
                                 if (p>i) then                                                              
  
                                    !auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(aa, ba, ib, ab, i, p )
                                    auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                    auxIndex = PropagatorTheory_IndexMapAB(ia, ba, ib, ab, i, p )

                                    auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
  
                                 else
  
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, aa, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, aa, ba, p, i )
                                    auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
                                    !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                    auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, ba, p, i )

                                    auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
  
  
                                 end if
  
                                 partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ia)&
                                      +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab))
                                 
                              end do
                           end do
                        end do
                        
                        constantSelfEnergy(i,p) = constantSelfEnergy(i,p) &
                             + 2.0_8*(auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
                        
                     end do
                  end do
  
  !JC                print *,"constant sigma after a-a-b:", constantSelfEnergy(i,p)
  
                  ! Three species
  
                  ! alpha-beta-alpha
  
                  ! alpha-beta-gamma
  
                  do r = 1, PropagatorTheory_instance%numberOfSpecies
  
                     if ( r /= p) then
  
  !JC                      print *,"entro a r diferente de p"
  
                        nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( r ) )
                        chargeOfSpeciesC = MolecularSystem_getCharge( r )
  !                      eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( r )
                        occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( r )
                        activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( r )
                        if ( InputCI_Instance(r)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesC = InputCI_Instance(r)%activeOrbitals
                        lambdaOfSpeciesC = MolecularSystem_getLambda( r )
                        virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
  
                        arguments(2) = trim(MolecularSystem_getNameOfSpecie(r))
                        arguments(1) = "ORBITALS"
                        call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( r ), &
                             unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                             output =  eigenValuesOfSpeciesC  )     
                        
                        do ib = 1 , occupationNumberOfSpeciesB
                           do jb = 1 , occupationNumberOfSpeciesB
  
                              if (p>i) then                                                              
  
                                 !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ib, jb, i, p )

                                 auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
  
                              else
  
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAB(ib, jb, pa, pa, p, i )

                                 auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
  
                              end if
  
                              partialValue = 0.0_8
                              
                              do ic = 1 , occupationNumberOfSpeciesC
                                 do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                    do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
  
                                       if (r>p) then                                                                                                   
                                          !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ic, ac, p, r )
                                          auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(jb, ab, ic, ac, p, r )

                                          auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
  
                                       else
  
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, ab, r, p )
                                          auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, jb, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, jb, ab, r, p )

                                          auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
  
                                       end if
  
                                       partialValue = partialValue &
                                            - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
                                            +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
                                            *( eigenValuesOfSpeciesB%values(jb)&
                                            + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac)))
                                       
                                    end do
                                 end do
                              end do
                              
                              constantSelfEnergy(p,r) = constantSelfEnergy(p,r) + partialValue*auxValue_E
                              
                           end do
                        end do
                        
                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                           do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
  
                              if (p>i) then
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ab, bb, i, p )

                                 auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
  
                              else
  
                                 !auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAB(ab, bb, pa, pa, p, i )

                                 auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
  
                              end if
  
                                 partialValue = 0.0_8
                              
                              do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
                                 do ib = 1 , occupationNumberOfSpeciesB
                                    do ic = 1 , occupationNumberOfSpeciesC
        
                                       if (r>p) then                                                                                                                                  

                                          !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ic, ac, p, r )
                                          auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ib, bb, ic, ac, p, r )

                                          auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
  
                                       else
  
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, ab, r, p )
                                          auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, bb, r, p )

                                          auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
  
                                       end if
  
                                       partialValue = partialValue &
                                            + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
                                            + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
                                            *( eigenValuesOfSpeciesB%values(ib)&
                                            + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac)))
                                       
                                    end do
                                 end do
                              end do
                              
                              constantSelfEnergy(p,r) = constantSelfEnergy(p,r) + partialValue*auxValue_E
                              
                           end do
                        end do
                        
                        do ib = 1 , occupationNumberOfSpeciesB
                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
  
                              if (p>i) then
                                 
                                 !auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                 auxIndex = PropagatorTheory_IndexMapAB(pa, pa, ib, ab, i, p )

                                 auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
  
                              else
  
                                 !auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                 auxIndex = PropagatorTheory_IndexMapAB(ib, ab, pa, pa, p, i  )

                                 auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
                                 
                              end if
  
                              partialValue = 0.0_8
                              
                              do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
                                 do jb = 1 , occupationNumberOfSpeciesB
                                    do ic =  1 , occupationNumberOfSpeciesC
  
                                       if (r>p) then
                                                                                  
                                          !auxIndex = IndexMap_tensorR4ToVector(ib, jb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ib, jb, ic, ac, p, r )
                                          auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(jb, ab, ic, ac, p, r )

                                          auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
  
                                       else
  
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, jb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, jb, r, p )
                                          auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, jb, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, jb, ab, r, p )

                                          auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
  
                                       end if
                                       
                                       partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(jb)&
                                            +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))
                                       
                                    end do
                                 end do
                              end do
                              
                              do ic = 1 , occupationNumberOfSpeciesC
                                 do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                    do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
  
                                       if (r>p) then
                                          
                                          !auxIndex = IndexMap_tensorR4ToVector(ab, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ab, bb, ic, ac, p, r )
                                          auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
                                          
                                          !auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                          auxIndex = PropagatorTheory_IndexMapAB(ib, bb, ic, ac, p, r )

                                          auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
  
                                       else
                                          
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ab, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ab, bb, r, p )
                                          auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
                                          
                                          !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                          auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, bb, r, p )

                                          auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
                                          
                                       end if
  
                                       partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ib)&
                                            +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac))
                                       
                                    end do
                                 end do
                              end do
                              
                              constantSelfEnergy(p,r) = constantSelfEnergy(p,r) &
                                   + 2.0_8*auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
                              
                           end do
                        end do
  
                     end if
  
  !JC                   print *,"constant sigma after a-b-a:", constantSelfEnergy(p,r)
  
                  end do
                  
               end if
                  
            end do

          end if 
          ! end of constant self-energy calculation
          
          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
          
          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))

          do j = 1 , PropagatorTheory_instance%numberOfSpecies             
             
             if (j==i) then ! Intraspecies factors
                
                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
                
                call Matrix_constructor(selfEnergy2ph(j), 3_8, vectorSize1, 0.0_8)
                call Matrix_constructor(selfEnergy2hp(j), 3_8, vectorSize2, 0.0_8)
                
                id1 = 0
                id2 = 0
                
                ! factor 2ph
                
                do ia = 1 , occupationNumberOfSpeciesA
                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         
                         !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ia, ba, i )
                         auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
                         !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
                         auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ia, aa, i )

                         auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
                         
                         id1 = id1 + 1
                         
                         selfEnergy2ph(j)%values(1,id1) = auxValue_A - auxValue_B
                         
                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
                              - eigenValuesOfSpeciesA%values(ba)
                         
                         valueOfW = 0.0_8
                         
                         do ja = 1, occupationNumberOfSpeciesA
                            do ka = 1, occupationNumberOfSpeciesA
                               
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, ka, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ia, ka, i )
                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ka, ia, ja, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ka, ia, ja, i )
                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, aa, ka, ba, i )
                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, ba, ka, aa, i )

                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                               
                               valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
                               
                            end do
                         end do
                         
                         do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ja = 1, occupationNumberOfSpeciesA
                               
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, aa, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ca, aa, i )
                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ca, ja, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ca, ja, i )
                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, ba, ia, ca, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, ba, ia, ca, i )
                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, ba, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, ca, ia, ba, i )
                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                               
                               valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
                                    - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
                               
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, ba, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ca, ba, i )
                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ca, ja, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ca, ja, i )
                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ia, ca, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, aa, ia, ca, i )
                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, aa, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(ja, ca, ia, aa, i )
                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                               
                               valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca) )
                               
                            end do
                         end do
                         
                         selfEnergy2ph(j)%values(3,id1) = valueOfW
                         
                      end do
                   end do
                end do
                
                if (occupationNumberOfSpeciesA > 1) then
                   
                   ! factor 2hp
                   
                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                      do ia = 1 , occupationNumberOfSpeciesA
                         do ja = 1 , occupationNumberOfSpeciesA
                            
                            id2 = id2 + 1
                            
                            !auxIndex = IndexMap_tensorR4ToVector(pa, ia, aa, ja, activeOrbitalsOfSpeciesA )
                            auxIndex = PropagatorTheory_IndexMapAA(pa, ia, aa, ja, i )
                            auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
                            !auxIndex = IndexMap_tensorR4ToVector(pa, ja, aa, ia, activeOrbitalsOfSpeciesA )
                            auxIndex = PropagatorTheory_IndexMapAA(pa, ja, aa, ia, i )
                            auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
                            
                            selfEnergy2hp(j)%values(1,id2) = auxValue_A - auxValue_B
                            
                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
                                 - eigenValuesOfSpeciesA%values(ja) 
                               
                               valueOfW = 0.0_8
                               
                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                     
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, ca, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ba, aa, ca, i )
                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ca, aa, ba, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ca, aa, ba, i )
                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ia, ca, ja, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ia, ca, ja, i )
                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, ia, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ja, ca, ia, i )
                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                                     
                                     valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
                                     
                                  end do
                               end do
                               
                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                  do ka = 1, occupationNumberOfSpeciesA
                                     
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ia, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ka, ia, i )
                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ka, ba, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ia, ka, ba, i )
                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ja, aa, ka, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ja, aa, ka, i )
                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ja, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ka, aa, ja, i )
                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                                     
                                     valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                          /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
                                     
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ja, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ka, ja, i )
                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ka, ba, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ka, ba, i )
                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ia, aa, ka, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ia, aa, ka, i )
                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ia, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAA(ba, ka, aa, ia, i )
                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
                                     
                                     valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ka) &
                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
                                     
                                  end do
                               end do
                               
                               selfEnergy2hp(j)%values(3,id2) = valueOfW
                               
                            end do
                         end do
                      end do
                      
                   end if
                   
                else ! interspecies
                   
                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
                   if ( InputCI_Instance(j)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(j)%activeOrbitals
                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB

                   arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))
                   arguments(1) = "ORBITALS"
                   call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( j ), &
                         unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                         output =  eigenValuesOfSpeciesB  )
                   
                   ! paso
                   
                   paso2=(nameOfSpeciesA=="E-ALPHA".and.nameOfSpeciesB=="E-BETA").or.&
                        (nameOfSpeciesA=="E-BETA".and.nameOfSpeciesB=="E-ALPHA")

                   vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
                   vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB
                   
                   call Matrix_constructor(selfEnergy2ph(j), 3_8, vectorSize1, 0.0_8)
                   call Matrix_constructor(selfEnergy2hp(j), 3_8, vectorSize2, 0.0_8)
                   
                   id1 = 0
                   id2 = 0
                   
                   ! diagram A
                   
                   do ib = 1 , occupationNumberOfSpeciesB
                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                            if (j>i) then
                               
                               !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                               auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ib, ab, i, j )
                               auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)

                            else
                               
                               !auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAB(ib, ab, pa, aa, j, i )
                               auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                               
                            end if

                            id1 = id1 + 1
                            
                            selfEnergy2ph(j)%values(1,id1) = auxValue_A
                            
                            selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
                                 - eigenValuesOfSpeciesB%values(ab)
                            
                            valueOfW = 0.0_8
                            
                            do ia = 1, occupationNumberOfSpeciesA
                               do jb = 1, occupationNumberOfSpeciesB

                                  if (j>i) then                                  

                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, jb, &
                                          !activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(pa, ia, ib, jb, i, j )
                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
                                          !activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(ia, aa, jb, ab, i, j )
                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                  else

                                     !auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, ia, &
                                          !activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(ib, jb, pa, ia, j, i )
                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(jb, ab, ia, aa, &
                                          !activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(jb, ab, ia, aa, j, i )
                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                  end if

                                  valueOfW = valueOfW + (auxValue_A * auxValue_B)&
                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
                                  
                               end do
                            end do
                            
                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                               do ia = 1 , occupationNumberOfSpeciesA

                                  if (j>i) then                                                                    
                                     
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, ab, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(pa, ia, bb, ab, i, j )
                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib,  bb, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ib,  bb, i, j )
                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                  else

                                     !auxIndex = IndexMap_tensorR4ToVector(bb, ab, pa, ia, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(bb, ab, pa, ia, j, i )
                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ib, bb, ia, aa, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(ib, bb, ia, aa, j, i )
                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                  end if

                                  valueOfW = valueOfW - (auxValue_A * auxValue_B)&
                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
                                  
                               end do
                            end do
                            
                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                               do ia = 1 , occupationNumberOfSpeciesA                                  

                                  !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, aa, activeOrbitalsOfSpeciesA )
                                  auxIndex = PropagatorTheory_IndexMapAA(pa, ia, ba, aa, i )
                                  auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
                                  !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ba, ia, activeOrbitalsOfSpeciesA )
                                  auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ba, ia, i )
                                  auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
                                  if (j>i) then
                                     !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(ia, ba, ib, ab, i, j )
                                     auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                  else
                                     !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ia, ba, j, i )
                                     auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
                                  end if
                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
                                  
                               end do
                            end do
                            
                            do jb = 1 , occupationNumberOfSpeciesB
                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                  
                                  !auxIndex = IndexMap_tensorR4ToVector(jb, ab, ib, bb, activeOrbitalsOfSpeciesB )
                                  auxIndex = PropagatorTheory_IndexMapAA(jb, ab, ib, bb, j )
                                  auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
                                  !auxIndex = IndexMap_tensorR4ToVector(jb, bb, ib, ab, activeOrbitalsOfSpeciesB )
                                  auxIndex = PropagatorTheory_IndexMapAA(jb, bb, ib, ab, j )
                                  auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
                                  if (j>i) then
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, aa, jb, bb, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(pa, aa, jb, bb, i, j)
                                     auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                  else
                                     !auxIndex = IndexMap_tensorR4ToVector(jb, bb, pa, aa, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(jb, bb, pa, aa, j, i)
                                     auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
                                  end if
                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
                                       /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
                                       - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
                                  
                               end do
                            end do

                            selfEnergy2ph(j)%values(3,id1) = valueOfW
                            
                         end do
                      end do
                   end do
                   
                   ! diagram B
                   
                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                      do ia = 1 , occupationNumberOfSpeciesA
                         do ib = 1 , occupationNumberOfSpeciesB
                            
                            id2 = id2 + 1
                            
                            if (j>i) then                                                                    

                               !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, &
                               !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                               auxIndex = PropagatorTheory_IndexMapAB(pa, ia, ab, ib, i, j)
                               auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)

                            else

                               !auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ia, &
                               !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAB(ab, ib, pa, ia, j, i)
                               auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)

                            end if

                            selfEnergy2hp(j)%values(1,id2) = auxValue_A
                            
                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
                                 - eigenValuesOfSpeciesB%values(ib)
                            
                            valueOfW = 0.0_8

                            do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                                  if (j>i) then                                                                                                                                           
                                     !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, bb, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ab, bb, i, j)
                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, bb, &
                                     !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ib, bb, i, j)
                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                  else

                                     !auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, aa, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(ab, bb, pa, aa, j, i)
                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(ib, bb, ia, aa, &
                                     !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                     auxIndex = PropagatorTheory_IndexMapAB(ib, bb, ia, aa, j, i)
                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                  end if

                               valueOfW = valueOfW + (auxValue_A * auxValue_B)&
                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
                               
                            end do
                         end do
                         
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do jb = 1 , occupationNumberOfSpeciesB
                               
                               if (j>i) then                                                                                                                                           
                                  !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, jb, &
                                  !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                  auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ib, jb, i, j )
                                  auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                  !auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
                                  !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                  auxIndex = PropagatorTheory_IndexMapAB(ia, aa, jb, ab, i, j )
                                  auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                               else

                                  !auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, aa, &
                                  !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                  auxIndex = PropagatorTheory_IndexMapAB(ib, jb, pa, aa, j, i )
                                  auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                  !auxIndex = IndexMap_tensorR4ToVector(jb, ab, ia, aa, &
                                  !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                  auxIndex = PropagatorTheory_IndexMapAB(jb, ab, ia, aa, j, i )
                                  auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                               end if

                               valueOfW = valueOfW - (auxValue_A * auxValue_B)&
                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
                               
                            end do
                         end do

                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ja = 1 , occupationNumberOfSpeciesA

                               !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ia, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ja, ia, i )
                               auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAA(pa, ia, ja, aa, i )
                               auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
                               if (j>i) then
                                  !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, &
                                  !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                  auxIndex = PropagatorTheory_IndexMapAB(ja, aa, ib, ab, i, j )
                                  auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                               else
                                  !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, &
                                  !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                  auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ja, aa, j, i )
                                  auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
                               end if
                               valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
                               
                            end do
                         end do

                         do jb = 1 , occupationNumberOfSpeciesB
                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                               
                               !auxIndex = IndexMap_tensorR4ToVector(ab, ib, bb, jb, activeOrbitalsOfSpeciesB )
                               auxIndex = PropagatorTheory_IndexMapAA(ab, ib, bb, jb, j )
                               auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
                               !auxIndex = IndexMap_tensorR4ToVector(ab, jb, bb, ib, activeOrbitalsOfSpeciesB )
                               auxIndex = PropagatorTheory_IndexMapAA(ab, jb, bb, ib, j )
                               auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
                               if (j>i) then
                               !auxIndex = IndexMap_tensorR4ToVector(pa, ia, jb, bb,&
                               !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                               auxIndex = PropagatorTheory_IndexMapAB(pa, ia, jb, bb, i, j )
                               auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                               else
                               !auxIndex = IndexMap_tensorR4ToVector(jb, bb, pa, ia,&
                               !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                               auxIndex = PropagatorTheory_IndexMapAB(jb, bb, pa, ia, j, i )
                               auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
                               end if
                               valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
                                    /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
                                    - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
                               
                            end do
                         end do
                         
                         selfEnergy2hp(j)%values(3,id2) = valueOfW
                         
                      end do
                   end do
                end do
                
             end if
             
          end do

          ! Initial guess
          koopmans = eigenValuesOfSpeciesA%values(pa)
          newOmega = koopmans
          lastOmega = 0.0_8
          
          ni = 0
          limit = 50
          residual = 1.0_8

          ! Calculation of second order pole


          write (*,"(T2,A10,A13,A15)") "Iteration ","  New Omega ","    Residual "
          do while ((residual>0.001_8).or.(limit.lt.ni))
             
             ni = ni + 1
             
             lastOmega = newOmega
             selfEnergy = lastOmega - koopmans
             selfEnergyDerivative = 1.0_8
             
             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
                
                if (j==i) then ! Intraspecies term
                   
                   id1=0
                   id2=0
                   
                   do ia = 1 , occupationNumberOfSpeciesA
                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            
                            id1 = id1 + 1
                            
                            a1 = selfEnergy2ph(j)%values(1,id1)
                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                            
                            selfEnergy = selfEnergy - 0.5_8*a1*a1/b
                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
                            
                         end do
                      end do
                   end do
                   
                   if (occupationNumberOfSpeciesA > 1) then
                      
                      ! factor 2hp
                      
                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         do ia = 1 , occupationNumberOfSpeciesA
                            do ja = 1 , occupationNumberOfSpeciesA
                               
                               id2 = id2 + 1
                               
                               a1 = selfEnergy2hp(j)%values(1,id2)
                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
                               
                               selfEnergy = selfEnergy - 0.5_8*a1*a1/b
                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
                                                           
                            end do
                         end do
                      end do
                      
                   end if
                   
                else ! Interspecies term

                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
                   if ( InputCI_Instance(j)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(j)%activeOrbitals
                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB

                   arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))
                   arguments(1) = "ORBITALS"
                   call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( j ), &
                        unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                        output =  eigenValuesOfSpeciesB  )     

                   id1 = 0
                   id2 = 0

                   ! diagram A
                   
                   do ib = 1 , occupationNumberOfSpeciesB
                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                            
                            id1 = id1 + 1

                            a1 = selfEnergy2ph(j)%values(1,id1)
                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                            
                            selfEnergy = selfEnergy - a1*a1/b
                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
                            
                         end do
                      end do
                   end do
                   
                   ! diagram B
                   
                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                      do ia = 1 , occupationNumberOfSpeciesA
                         do ib = 1 , occupationNumberOfSpeciesB
                            
                            id2 = id2 + 1

                            a1 = selfEnergy2hp(j)%values(1,id2)
                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
                            
                            selfEnergy = selfEnergy -a1*a1/b
                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
                            
                         end do
                      end do
                   end do
                   
                end if
                
             end do
             
             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
             
             residual = abs(newOmega-lastOmega)
             
!             print *,"Iteration ",ni,"New Omega ",newOmega,"Residual ",residual
             write (*,"(T5,I2,A5,E12.5E2,A5,E12.5E2)") ,ni,"     ",newOmega,"     ",residual
          end do ! while

          poleStrenght = 1.0_8/selfEnergyDerivative

          ! Storing corrections

!          print *,"M y Q:",m,q

          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)=real(pa,8)
          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans
          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3)=27.211396_8 * newOmega
          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)=poleStrenght

          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)),&
               " of species: ",nameOfSpeciesA
          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
          write (*,"(T5,A29,F8.4,A7,I2,A12)") "Optimized second order pole: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
               " after ",ni," iterations."
          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
          print *,"----------------------------------------------------------------"

          ! Calculation of third order poles

          ! factors for different algorythms

          factors(:,:,:) = 0.0_8
          factors2(:,:,:) = 0.0_8
          thirdOrderResults(:,:) = 0.0_8

          ! o=1 P3
          thirdOrderMethods(1)="P3"
          ! o=2 EP3
          thirdOrderMethods(2)="EP3"
          ! o=3 OVGF version A
          thirdOrderMethods(3)="OVGF-A"
          ! o=4 OVGF version B
          thirdOrderMethods(4)="OVGF-B"
          ! o=5 OVGF version C
          thirdOrderMethods(5)="OVGF-C"
          ! o=6 
          thirdOrderMethods(6)="REN-P3"


          
          !do o = 1 , 6 ! Options for third order and renormalized third order
          do oo = 1, maxoo
             o = oarray(oo)

             ! Initial guess             
             koopmans = eigenValuesOfSpeciesA%values(pa)
             newOmega = koopmans
             ! if (o>2.and.o/=6) newOmega = thirdOrderResults(1,2)/27.211396_8
             ! if (o==6) newOmega = thirdOrderResults(1,1)/27.211396_8
             lastOmega = 0.0_8
             
             ni = 0
             limit = 15
             residual = 1.0_8

             ! factor for W
             
             fW = 2.0_8
             fI = 1.0_8
             if (o==1 .or. o==6) fW=1.0_8
             if (o==1 .or. o==6) fI=0.0_8
             threshold=0.001_8/27.211396_8
             ! if (o==1 .or. o==2) threshold=0.00001_8
             
             ! NR procedure Hola

             write (*,"(T2,A10,A13,A15)") "Iteration ","  New Omega ","    Residual "
             do while ((residual>threshold))
                
                ni = ni + 1
                
                lastOmega = newOmega
                selfEnergy = lastOmega - koopmans
                selfEnergyDerivative = 1.0_8
                s2hp(:) = 0.0_8
                s2ph(:) = 0.0_8
                W2hp(:) = 0.0_8
                W2ph(:) = 0.0_8
                U2hp(:) = 0.0_8
                U2ph(:) = 0.0_8
                ! new terms
                W22hp(:) = 0.0_8
                W22ph(:) = 0.0_8
                U22hp(:) = 0.0_8
                U22ph(:) = 0.0_8
                
                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
                   
                   if (j==i) then ! Intraspecies term
                      
                      id1=0
                      id2=0

                      sub2 = 0.0_8
                      subW = 0.0_8 
                      subU = 0.0_8 
                      subd2 = 0.0_8
                      subdW = 0.0_8
                      subdU = 0.0_8                      

                      ! 2ph terms a-a-a
                      if ( (.not.paso1).or.(o/=1.and.o/=6) ) then

                      id1=0
                      do ia = 1 , occupationNumberOfSpeciesA
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                               
                               id1 = id1 + 1

                               valueOfU = 0.0_8
                               valueOfdU = 0.0_8

                               !if ( (.not.paso1).or.(o/=1.and.o/=6) ) then
                                  
                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                     do da = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ca, ia, da, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ca, ia, da, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, da, ia, ca, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, da, ia, ca, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ca, aa, da, ba, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ca, aa, da, ba, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ca, ba, da, aa, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ca, ba, da, aa, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(ia) &
                                             - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(da)
                                        
                                        valueOfU = valueOfU + 0.5_8*a2/c
                                        valueOfdU = valueOfdU - 0.5_8*a2/(c**2.0_8)
                                        
                                     end do
                                  end do
                                  
                                  do ja = 1 , occupationNumberOfSpeciesA
                                     do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ja, ca, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ba, ja, ca, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, ba, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ca, ja, ba, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, aa, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ia, ja, ca, aa, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ca, ja, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ia, aa, ca, ja, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca)
                                        
                                        valueOfU = valueOfU + a2/c
                                        valueOfdU = valueOfdU - a2/(c**2.0_8)
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ca, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, aa, ja, ca, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, aa, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ca, ja, aa, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, ba, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ia, ja, ca, ba, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ca, ja, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ia, ba, ca, ja, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca)
                                        
                                        valueOfU = valueOfU - a2/c
                                        valueofdU = valueOfdU + a2/(c**2.0_8)
                                        
                                     end do
                                  end do                                  

                               !end if

                               a1 = selfEnergy2ph(j)%values(1,id1)
                               a2 = selfEnergy2ph(j)%values(3,id1)
                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                               
                               sub2 = sub2 + (a1**2.0_8)/b
                               subW = subW + (a1*a2)/b
                               subU = subU + (a1*valueOfU)/b
                               
                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                               subdW = subdW + (a1*a2)/(b**2.0_8)
                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                                                        
                            end do
                         end do
                      end do
                      end if

                      if (paso1.and.(o==1.or.o==6)) then
  
                      id1=0
                      do ia = 1 , occupationNumberOfSpeciesA
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                               
                               id1 = id1 + 1

                               valueOfdU = 0.0_8

                               a1 = selfEnergy2ph(j)%values(1,id1)
                               a2 = selfEnergy2ph(j)%values(3,id1)
                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                               
                               sub2 = sub2 + (a1**2.0_8)/b
                               
                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                           
                            end do
                         end do
                      end do

                      subdW = 0.0_8
                      subdU = 0.0_8                                                       
                      subW = 0.0_8
                      subU = 0.0_8 

                      end if

                      !if (paso1.and.(o==1.or.o==6)) subW=0.0_8
                      !if (paso1.and.(o==1.or.o==6)) subdW=0.0_8                      

                      s2ph(j) = s2ph(j) + 0.5_8*sub2
                      W2ph(j) = W2ph(j) + 0.5_8*(fW*subW)/(1.0_8-factors(i,1,o))                       
                      U2ph(j) = U2ph(j) + 0.5_8*(subU)/(1.0_8-factors(i,1,o))                       

                      selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(i,1,o)) )                       
                      
                      selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(i,1,o)) )                     
                      ! Diagram 2hp a-a-a
                        
                      if (occupationNumberOfSpeciesA > 1) then

                         sub2 = 0.0_8
                         subW = 0.0_8 
                         subU = 0.0_8 
                         subd2 = 0.0_8
                         subdW = 0.0_8
                         subdU = 0.0_8
                         
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ia = 1 , occupationNumberOfSpeciesA
                               do ja = 1 , occupationNumberOfSpeciesA
                                  
                                  id2 = id2 + 1
                                  
                                  valueOfU = 0.0_8
                                  valueOfdU = 0.0_8

                                  do ka = 1 , occupationNumberOfSpeciesA
                                     do la = 1 , occupationNumberOfSpeciesA
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ka, aa, la, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ka, aa, la, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, la, aa, ka, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, la, aa, ka, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ka, ia, la, ja, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ka, ia, la, ja, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ka, ja, la, ia, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(ka, ja, la, ia, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(aa) &
                                             - eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(la)
                                        
                                        valueOfU = valueOfU - 0.5_8*a2/c
                                        valueOfdU = valueOfdU + 0.5_8*a2/(c**2.0_8)
                                        
                                     end do
                                  end do
                                  
                                  do ka = 1 , occupationNumberOfSpeciesA
                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ba, ka, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ja, ba, ka, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ja, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ka, ba, ja, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ia, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(aa, ba, ka, ia, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(aa, ia, ka, ba, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(aa, ia, ka, ba, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
                                             - eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ka)
                                        
                                        valueOfU = valueOfU - a2/c
                                        valueOfdU = valueOfdU + a2/(c**2.0_8)
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, ka, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ia, ba, ka, i )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ia, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(pa, ka, ba, ia, i )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ja, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(aa, ba, ka, ja, i )
                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(aa, ja, ka, ba, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAA(aa, ja, ka, ba, i )
                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
                                             - eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ka)
                                        
                                        valueOfU = valueOfU + a2/c
                                        valueofdU = valueOfdU - a2/(c**2.0_8)
                                        
                                     end do
                                  end do
                                                                                                   
                                  a1 = selfEnergy2hp(j)%values(1,id2)
                                  a2 = selfEnergy2hp(j)%values(3,id2)
                                  b = selfEnergy2hp(j)%values(2,id2) + lastOmega
                                  
                                  sub2 = sub2 + (a1**2.0_8)/b
                                  subW = subW + (a1*a2)/b
                                  subU = subU + (a1*valueOfU)/b
                                  
                                  subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                                  subdW = subdW + (a1*a2)/(b**2.0_8)
                                  subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                                  
                               end do
                            end do
                         end do
                         
                         s2hp(i) = s2hp(i) + 0.5_8*sub2
                         W2hp(i) = W2hp(i) + 0.5_8*(fW*subW)/(1.0_8-factors(i,2,o))                       
                         U2hp(i) = U2hp(i) + 0.5_8*(subU)/(1.0_8-factors(i,2,o))                       
                         
                         selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(i,2,o)) )                       
                         
                         selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(i,2,o)) )                                                                                           
                      end if

                      ! terms a-a-b
                         
                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
                         
                         if (k .ne. i)  then
                            
                            id1=0
                            id2=0

                            nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
                            chargeOfSpeciesB = MolecularSystem_getCharge( k )
!                            eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
                            occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
                            activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
                            if ( InputCI_Instance(k)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(k)%activeOrbitals
                            lambdaOfSpeciesB = MolecularSystem_getLambda( k )
                            virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB

                            arguments(2) = trim(MolecularSystem_getNameOfSpecie(k))
                            arguments(1) = "ORBITALS"
                            call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( k ), &
                                     unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                                     output =  eigenValuesOfSpeciesB  )     

                            ! 2ph a-a-b

                            if ( (.not.paso1).or.(o/=1.and.o/=6)) then                            

!                               print *,"ENTRA A LA PARTE 2PH A-A-B CON A:",i,"Y B:",k

                               subU = 0.0_8 
                               subdU = 0.0_8                      
                               subW = 0.0_8
                               subdW = 0.0_8
                            
                               do ia = 1 , occupationNumberOfSpeciesA
                                  do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                        
                                        id1 = id1 + 1

                                        valueOfU = 0.0_8
                                        valueOfdU = 0.0_8
                                        valueOfW = 0.0_8
                                                                                
                                        do ib = 1 , occupationNumberOfSpeciesB
                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                                              if (k>i) then
                                                 !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib,&
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(pa, ba, ab, ib, i, k )
                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib,&
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ab, ib, i, k )
                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
                                              else
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ba,&
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, pa, ba, k, i )
                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, aa, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, ia, aa, k, i )
                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
                                              end if

                                              a2 = auxValue_A*auxValue_B
                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)

                                              valueOfU = valueOfU - 2.0_8*a2/c
                                              valueOfW = valueOfW - a2/d
                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
                                              
                                              if (k>i) then
                                                 !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ab, ib, i, k )
                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ia, ba, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ia, ba, ab, ib, i, k )
                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
                                              else
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, aa, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, pa, aa, k, i )
                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, ba, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, ia, ba, k, i )
                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
                                              end if
                                                 
                                              a2 = auxValue_A*auxValue_B
                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)

                                              valueOfU = valueOfU + 2.0_8*a2/c
                                              valueOfW = valueOfW + a2/d
                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
                                              
                                           end do
                                        end do
                                        
                                        a1 = selfEnergy2ph(j)%values(1,id1)
                                        a2 = selfEnergy2ph(j)%values(3,id1)
                                        b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                                        
                                        subU = subU + (a1*valueOfU)/b                                        
                                        subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                                        
                                        subW = subW + (a1*valueOfW)/b
                                        subdW = subdW + (a1*valueOfW)/(b**2.0_8)

                                        ! print *,"a-a-b 2ph subU:",subU
                                        ! print *,"a-a-b 2ph subU:",subW
                                        
                                     end do
                                  end do
                               end do

                               U22ph(k) = U22ph(k) + 0.5_8*(subU)/(1.0_8-factors2(k,1,o))                    
                               W22ph(k) = W22ph(k) + 0.5_8*(fW*subW)/(1.0_8-factors2(k,1,o))                          

                               selfEnergy = selfEnergy - 0.5_8*(subU+fW*subW)/(1.0_8-factors2(k,1,o))                       
                               
                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*(subdU+fW*subdW)/(1.0_8-factors2(k,1,o))                     

                            end if
                            
                            ! 2hp a-a-b

                            if (occupationNumberOfSpeciesA > 1) then

                               subU = 0.0_8 
                               subdU = 0.0_8     
                               subW = 0.0_8
                               subdW = 0.0_8
                               
                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                  do ia = 1 , occupationNumberOfSpeciesA
                                     do ja = 1 , occupationNumberOfSpeciesA
                                        
                                        id2 = id2 + 1
                                        
                                        valueOfU = 0.0_8
                                        valueOfdU = 0.0_8
                                        valueOfW = 0.0_8
                                        
                                        do ib = 1 , occupationNumberOfSpeciesB
                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                                              if (k>i) then
                                                 
                                                 !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(pa, ja, ab, ib, i, k )
                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ia, aa, ab, ib, i, k)
                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 
                                              else

                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ja, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, pa, ja, k, i )
                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, aa, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, ia, aa, k, i )
                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)

                                              end if
                                                    
                                              a2 = auxValue_A*auxValue_B
                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
                                              
                                              valueOfU = valueOfU + 2.0_8*a2/c
                                              valueOfW = valueOfW - a2/d
                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)

                                              if (k>i) then                                              
                                                 
                                                 !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(pa, ia, ab, ib, i, k )
                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ja, aa, ab, ib, &
                                                 !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ja, aa, ab, ib, i, k )
                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
                                                 
                                              else

                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ia, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, pa, ia, k, i )
                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                                 !auxIndex = IndexMap_tensorR4ToVector(ab, ib, ja, aa, &
                                                 !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                                 auxIndex = PropagatorTheory_IndexMapAB(ab, ib, ja, aa, k, i )
                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)

                                              end if

                                              a2 = auxValue_A*auxValue_B
                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ia)
                                              d = eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
                                              
                                              valueOfU = valueOfU - 2.0_8*a2/c
                                              valueOfW = valueOfW + a2/d
                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
                                              
                                           end do
                                        end do
                                        
                                        a1 = selfEnergy2hp(j)%values(1,id2)
                                        a2 = selfEnergy2hp(j)%values(3,id2)
                                        b = selfEnergy2hp(j)%values(2,id2) + lastOmega

                                        subU = subU + (a1*valueOfU)/b
                                        subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)

                                        subW = subW + (a1*valueOfW)/b
                                        subdW = subdW + (a1*valueOfW)/(b**2.0_8)

                                        ! print *,"a-a-b 2hp subU:",subU
                                        ! print *,"a-a-b 2hp subU:",subW
                                                 
                                     end do
                                  end do
                               end do
                               
                               U22hp(k) = U22hp(k) + 0.5_8*(subU)/(1.0_8-factors2(k,2,o))                       
                               W22hp(k) = W22hp(k) + 0.5_8*(fW*subW)/(1.0_8-factors2(k,2,o))                          
                               
                               selfEnergy = selfEnergy - 0.5_8*(subU+fW*subW)/(1.0_8-factors2(k,2,o))                       
                               
                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*(subdU+fW*subdW)/(1.0_8-factors2(k,2,o))                     
                               
                            end if

                         end if
                            
                      end do
                     
                      selfEnergy = selfEnergy - fI*constantSelfEnergy(i,i)/(1.0_8-factors(i,3,o))                
                      
                   else ! Interspecies term
                      
                      nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
                      chargeOfSpeciesB = MolecularSystem_getCharge( j )
!                      eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
                      occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
                      activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
                      if ( InputCI_Instance(j)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesB = InputCI_Instance(j)%activeOrbitals
                      lambdaOfSpeciesB = MolecularSystem_getLambda( j )
                      virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB

                      arguments(2) = trim(MolecularSystem_getNameOfSpecie(j))
                      arguments(1) = "ORBITALS"
                      call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( j ), &
                           unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                           output =  eigenValuesOfSpeciesB  )   
                      
                      paso2=(nameOfSpeciesA=="E-ALPHA".and.nameOfSpeciesB=="E-BETA").or.&
                           (nameOfSpeciesA=="E-BETA".and.nameOfSpeciesB=="E-ALPHA")
                      
                      id1 = 0
                      id2 = 0
                      
                      ! Diagram 2ph interspecies

                      sub2 = 0.0_8
                      subW = 0.0_8 
                      subU = 0.0_8 
                      subd2 = 0.0_8
                      subdW = 0.0_8
                      subdU = 0.0_8
                      
                      !if ( (.not.paso2) .or. (ooarray(2) ==2) ) then
                      if ( (.not.paso2).or.(o/=1.and.o/=6) ) then
                      id1 = 0
                      do ib = 1 , occupationNumberOfSpeciesB
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                               
                               id1 = id1 + 1
                               
                               valueOfU = 0.0_8
                               valueOfdU = 0.0_8
                               
                               !if ( (.not.paso2).or.(o/=1.and.o/=6) ) then
                                  
                                  do jb = 1 , occupationNumberOfSpeciesB
                                     
                                     do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                                        if (j>i) then
                                           
                                           !auxIndex = IndexMap_tensorR4ToVector(pa, aa, bb, jb, &
                                           !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                           auxIndex = PropagatorTheory_IndexMapAB(pa, aa, bb, jb, i, j )
                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)

                                        else

                                           !auxIndex = IndexMap_tensorR4ToVector(bb, jb, pa, aa, &
                                           !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                           auxIndex = PropagatorTheory_IndexMapAB(bb, jb, pa, aa, j, i )
                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)

                                        end if

                                        !auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, ab, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAA(ib, jb, bb, ab, j )
                                        auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ib, ab, bb, jb, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAA(ib, ab, bb, jb, j )
                                        auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
                                        
                                        a2 = (auxValue_A)*(auxValue_B - auxValue_C)
                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(aa)
                                        
                                        valueOfU = valueOfU - a2/c
                                        valueofdU = valueOfdU + a2/(c**2.0_8)
                                        
                                     end do
                                  end do
                                  
                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                        
                                        if (j>i) then
                                           
                                           !auxIndex = IndexMap_tensorR4ToVector(ba, aa, bb, ab, &
                                           !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                           auxIndex = PropagatorTheory_IndexMapAB(ba, aa, bb, ab, i, j )
                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                           !auxIndex = IndexMap_tensorR4ToVector(pa, ba, ib, bb, &
                                           !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                           auxIndex = PropagatorTheory_IndexMapAB(pa, ba, ib, bb, i, j )
                                           auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                        else

                                           !auxIndex = IndexMap_tensorR4ToVector(bb, ab, ba, aa, &
                                           !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                           auxIndex = PropagatorTheory_IndexMapAB(bb, ab, ba, aa, j, i )
                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                           !auxIndex = IndexMap_tensorR4ToVector(ib, bb, pa, ba, &
                                           !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                           auxIndex = PropagatorTheory_IndexMapAB(ib, bb, pa, ba, j, i )
                                           auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                        end if

                                        a2 = auxValue_A*auxValue_B
                                        c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(ba)
                                        
                                        valueOfU = valueOfU + (a2/c)
                                        valueofdU = valueOfdU - (a2/(c**2.0_8))
                                        
                                     end do
                                  end do
                                  
                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                     do jb = 1 , occupationNumberOfSpeciesB

                                        if (j>i) then                                        
                                           
                                           !auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, jb, &
                                           !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                           auxIndex = PropagatorTheory_IndexMapAB(aa, ba, ib, jb, i, j )
                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                           !auxIndex = IndexMap_tensorR4ToVector(pa, ba, jb, ab, &
                                           !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                           auxIndex = PropagatorTheory_IndexMapAB(pa, ba, jb, ab, i, j )
                                           auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                        else

                                           !auxIndex = IndexMap_tensorR4ToVector(ib, jb, aa, ba, &
                                           !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                           auxIndex = PropagatorTheory_IndexMapAB(ib, jb, aa, ba, j, i )
                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                           !auxIndex = IndexMap_tensorR4ToVector(jb, ab, pa, ba, &
                                           !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                           auxIndex = PropagatorTheory_IndexMapAB(jb, ab, pa, ba, j, i )
                                           auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                        end if

                                        a2 = auxValue_A*auxValue_B
                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
                                             - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
                                        
                                        valueOfU = valueOfU - a2/c
                                        valueofdU = valueOfdU + a2/(c**2.0_8)
                                        
                                     end do
                                  end do
                                  
                               !end if
                               
                               a1 = selfEnergy2ph(j)%values(1,id1)
                               a2 = selfEnergy2ph(j)%values(3,id1)
                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                               
                               sub2 = sub2 + (a1**2.0_8)/b
                               subW = subW + (a1*a2)/b
                               subU = subU + (a1*valueOfU)/b
                               
                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                               subdW = subdW + (a1*a2)/(b**2.0_8)
                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                               
                            end do
                         end do
                      end do
                      end if

                      if ( (paso2).and. ( o ==1 .or. o==6) ) then
                      !if ( (paso2) .or. (ooarray(1) ==1) ) then
                      id1 = 0
                      do ib = 1 , occupationNumberOfSpeciesB
                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                               
                               id1 = id1 + 1

                               valueOfdU = 0.0_8
                               
                               a1 = selfEnergy2ph(j)%values(1,id1)
                               a2 = selfEnergy2ph(j)%values(3,id1)
                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                               
                               sub2 = sub2 + (a1**2.0_8)/b
                               
                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                               
                            end do
                         end do
                      end do
                      subW = 0.0_8
                      subU = 0.0_8
                      subdW = 0.0_8
                      subdU = 0.0_8
                      
                      end if

                      !if (paso2.and.(o==1.or.o==6)) subW=0.0_8
                      !if (paso2.and.(o==1.or.o==6)) subdW=0.0_8                      

                      s2ph(j) = s2ph(j) + sub2
                      W2ph(j) = W2ph(j) + (fW*subW)/(1.0_8-factors(j,1,o))                       
                      U2ph(j) = U2ph(j) + (subU)/(1.0_8-factors(j,1,o))                       
                      
                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(j,1,o)) )                       
                      
                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(j,1,o)) )                                      
                      ! Diagram 2hp interspecies

                      sub2 = 0.0_8
                      subW = 0.0_8 
                      subU = 0.0_8 
                      subd2 = 0.0_8
                      subdW = 0.0_8
                      subdU = 0.0_8
                      
                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                         do ia = 1 , occupationNumberOfSpeciesA
                            do ib = 1 , occupationNumberOfSpeciesB
                               
                               id2 = id2 + 1
                               
                               valueOfU = 0.0_8
                               valueOfdU = 0.0_8
                               
                               do jb = 1 , occupationNumberOfSpeciesB
                                  
                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB

                                     if (j>i) then
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, jb,&
                                        !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAB(pa, ia, bb, jb, i, j )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)

                                     else

                                        !auxIndex = IndexMap_tensorR4ToVector(bb, jb, pa, ia, &
                                        !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAB(bb, jb, pa, ia, j, i )
                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                        
                                     end if

                                     !auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAA(bb, ab, ib, jb, j )
                                     auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
                                     !auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
                                     auxIndex = PropagatorTheory_IndexMapAA(bb, jb, ib, ab, j )
                                     auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
                                     
                                     a2 = (auxValue_A)*(auxValue_B - auxValue_C)
                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ia)
                                     
                                     valueOfU = valueOfU + a2/c
                                     valueofdU = valueOfdU - a2/(c**2.0_8)
                                     
                                  end do
                               end do
                               
                               do jb = 1 , occupationNumberOfSpeciesB
                                  do ja = 1 , occupationNumberOfSpeciesA
                                      
                                     if (j>i) then
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, jb,&
                                        !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAB(ia, ja, ib, jb, i, j )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, jb,&
                                        !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAB(pa, ja, ab, jb, i, j )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                     else

                                        !auxIndex = IndexMap_tensorR4ToVector(ib, jb, ia, ja, &
                                        !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAB(ib, jb, ia, ja, j, i )
                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(ab, jb, pa, ja,&
                                        !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAB(ab, jb, pa, ja, j, i )
                                        auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                     end if

                                     a2 = auxValue_A*auxValue_B
                                     c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ja)
                                     
                                     valueOfU = valueOfU - a2/c
                                     valueofdU = valueOfdU + a2/(c**2.0_8)
                                     
                                  end do
                               end do
                               
                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                  do ja = 1 , occupationNumberOfSpeciesA

                                     if (j>i) then                                     
                                        
                                        !auxIndex = IndexMap_tensorR4ToVector(ia, ja, bb, ab,&
                                        !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAB(ia, ja, bb, ab, i, j )
                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(pa, ja, bb, ib,&
                                        !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
                                        auxIndex = PropagatorTheory_IndexMapAB(pa, ja, bb, ib, i, j )
                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)

                                     else

                                        !auxIndex = IndexMap_tensorR4ToVector(bb, ab, ia, ja,&
                                        !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAB(bb, ab, ia, ja, j, i )
                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
                                        !auxIndex = IndexMap_tensorR4ToVector(bb, ib, pa, ja,&
                                        !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
                                        auxIndex = PropagatorTheory_IndexMapAB(bb, ib, pa, ja, j, i )
                                        auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)

                                     end if

                                     a2 = auxValue_A*auxValue_B
                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
                                          - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
                                     
                                     valueOfU = valueOfU + a2/c
                                     valueofdU = valueOfdU - a2/(c**2.0_8)
                                     
                                  end do
                               end do
                               
                               a1 = selfEnergy2hp(j)%values(1,id2)
                               a2 = selfEnergy2hp(j)%values(3,id2)
                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
                               
                               sub2 = sub2 + (a1**2.0_8)/b
                               subW = subW + (a1*a2)/b
                               subU = subU + (a1*valueOfU)/b
                               
                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
                               subdW = subdW + (a1*a2)/(b**2.0_8)
                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                               
                            end do
                         end do
                      end do
                      
                      s2hp(j) = s2hp(j) + sub2
                      W2hp(j) = W2hp(j) + (fW*subW)/(1.0_8-factors(j,2,o))                       
                      U2hp(j) = U2hp(j) + (subU)/(1.0_8-factors(j,2,o))                       
                      
                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(j,2,o)) )                       
                      
                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(j,2,o)) )                                      
                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
                         
                         id1=0
                         id2=0
                         
                         if (k.ne.i .and. k.ne.j)  then
                            
!                            print *,"ENTRO AL TERMINO DE TRES PARTICULAS:",i,j,k
                            
                            nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( k ) )
                            chargeOfSpeciesC = MolecularSystem_getCharge( k )
!                            eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( k )
                            occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( k )
                            activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( k )
                            if ( InputCI_Instance(k)%activeOrbitals /= 0 ) activeOrbitalsOfSpeciesC = InputCI_Instance(k)%activeOrbitals
                            lambdaOfSpeciesC = MolecularSystem_getLambda( k )
                            virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC

                            arguments(2) = trim(MolecularSystem_getNameOfSpecie(k))
                            arguments(1) = "ORBITALS"
                            call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions( k ), &
                                    unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                                    output =  eigenValuesOfSpeciesC  )    
                            
                            paso3=(nameOfSpeciesB=="E-ALPHA".and.nameOfSpeciesC=="E-BETA").or.&
                                 (nameOfSpeciesB=="E-BETA".and.nameOfSpeciesC=="E-ALPHA")
                        
                            subU=0.0_8
                            subdU=0.0_8                            
                            
                            do ib = 1 , occupationNumberOfSpeciesB
                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
                                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                                     
                                     id1 = id1 + 1
                                     
                                     valueOfU=0.0_8
                                     valueOfdU=0.0_8
                                     
                                     do ic = 1 , occupationNumberOfSpeciesC
                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
                                           
                                           if (k>i) then
                                              !auxIndex = IndexMap_tensorR4ToVector(pa, aa, ic, ac,&
                                              !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
                                              auxIndex = PropagatorTheory_IndexMapAB(pa, aa, ic, ac, i, k )
                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)                                              
                                           else
                                              !auxIndex = IndexMap_tensorR4ToVector(ic, ac, pa, aa,&
                                              !     activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesA )
                                              auxIndex = PropagatorTheory_IndexMapAB(ic, ac, pa, aa, k, i )
                                              auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                           end if
                                           if (k>j) then
                                              !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac,&
                                              !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                              auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ic, ac, j, k )
                                              auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
                                           else
                                              !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab,&
                                              !     activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                              auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, ab, k , j)
                                              auxValue_B = auxMatrix2(k,j)%values(auxIndex, 1)
                                           end if
                                              
                                           a2 = auxValue_A*auxValue_B
                                           c = lastOmega + eigenValuesOfSpeciesC%values(ic) &
                                                - eigenValuesOfSpeciesC%values(ac) - eigenValuesOfSpeciesA%values(aa)
                                           
                                           
                                           valueOfU = valueOfU - a2/c
                                           valueofdU = valueOfdU + a2/(c**2.0_8)
                                           
                                        end do
                                     end do
                                     
                                     a1 = selfEnergy2ph(j)%values(1,id1)
                                     b = selfEnergy2ph(j)%values(2,id1) + lastOmega
                                     
                                     subU = subU + (a1*valueOfU)/b
                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                                     
                                  end do
                               end do
                            end do

                            if (paso3) then
                               
                               U2ph(j)=U2ph(j)+subU
                               selfEnergy = selfEnergy - (subU)/(1.0_8-factors(j,1,o))
                               selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(j,1,o))

                            else

                               selfEnergy = selfEnergy - (subU)
                               selfEnergyDerivative = selfEnergyDerivative + subdU

                            end if

                            subU=0.0_8
                            subdU=0.0_8
                            
                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
                               do ia = 1 , occupationNumberOfSpeciesA
                                  do ib = 1 , occupationNumberOfSpeciesB
                                     
                                     id2 = id2 + 1
                                     
                                     valueOfU=0.0_8
                                     valueOfdU=0.0_8
                                     
                                     do ic = 1 , occupationNumberOfSpeciesC
                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
                                           
                                           if (k>i) then
                                              !auxIndex = IndexMap_tensorR4ToVector(pa, ia, ic, ac,&
                                              !     activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
                                              auxIndex = PropagatorTheory_IndexMapAB(pa, ia, ic, ac, i, k )
                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
                                           else
                                              !auxIndex = IndexMap_tensorR4ToVector(ic, ac, pa, ia,&
                                              !     activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesA )
                                              auxIndex = PropagatorTheory_IndexMapAB(ic, ac, pa, ia, k, i )
                                              auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
                                           end if
                                           if (k>j) then
                                              !auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac,&
                                              !     activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
                                              auxIndex = PropagatorTheory_IndexMapAB(ib, ab, ic, ac, j, k )
                                              auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
                                           else
                                              !auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab,&
                                              !     activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
                                              auxIndex = PropagatorTheory_IndexMapAB(ic, ac, ib, ab, k, j )
                                              auxValue_B = auxMatrix2(k,j)%values(auxIndex, 1)
                                           end if

                                           a2 = auxValue_A*auxValue_B
                                           c = lastOmega + eigenValuesOfSpeciesC%values(ac) &
                                                - eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesA%values(ia)
                                           
                                           valueOfU = valueOfU + a2/c
                                           valueofdU = valueOfdU - a2/(c**2.0_8)
                                           
                                        end do
                                     end do
                                     
                                     a1 = selfEnergy2hp(j)%values(1,id2)
                                     b = selfEnergy2hp(j)%values(2,id2) + lastOmega

                                     subU = subU + (a1*valueOfU)/b
                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
                                     
                                  end do
                               end do
                            end do

                            if (paso3) then

                               U2hp(j)=U2hp(j)+subU
                               selfEnergy = selfEnergy - (subU)/(1.0_8-factors(j,2,o))
                               selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(j,2,o))
                               selfEnergy = selfEnergy - fI*constantSelfEnergy(j,k)/(1.0_8-factors(j,3,o))

                            else

                               selfEnergy = selfEnergy - (subU)
                               selfEnergyDerivative = selfEnergyDerivative + subdU
                               selfEnergy = selfEnergy - fI*constantSelfEnergy(j,k)

                            end if
                            
                         end if
                         
                      end do

                      selfEnergy = selfEnergy - fI*constantSelfEnergy(j,j)/(1.0_8-factors(j,3,o))                
                      selfEnergy = selfEnergy - fI*constantSelfEnergy(i,j)/(1.0_8-factors(j,3,o))                
                      selfEnergy = selfEnergy - fI*constantSelfEnergy(j,i)/(1.0_8-factors(j,3,o))                
                                            
                   end if
                   
                end do
                
                newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
                
                residual = abs(newOmega-lastOmega)
                
                write (*,"(T5,I2,A5,E12.5E2,A5,E12.5E2)") ,ni,"     ",newOmega,"     ",residual
             end do ! while

             if (o==1) then

                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
                   
                   if (j/=electrons(1).and.j/=electrons(2)) then

                      value1=s2ph(j)
                      value2=W2ph(j)
                      value3=s2hp(j)
                      value4=W2hp(j)
                   
                      ! Renormalized P3
                      if (value1/=0.0_8) factors(j,1,6) = value2/value1
                      if (value3/=0.0_8) factors(j,2,6) = value4/value3
                      factors(j,3,6) = 0.0_8

                      value1=s2ph(i)
                      value2=W22ph(j)
                      value3=s2hp(i)
                      value4=W22hp(j)
                   
                      ! Renormalized P3
                      if (value1/=0.0_8) factors2(j,1,6) = value2/value1
                      if (value3/=0.0_8) factors2(j,2,6) = value4/value3
                      factors2(j,3,6) = 0.0_8
                      
                   end if
                   
                end do

                r = electrons(1)
                s = electrons(2)
!                print *,"values of r and s:",r,s
!                print *,"value of i:",i

                if (r/=0.and.s/=0) then
                   
                   if (r==i.or.s==i) then

!                      print *,"entro a la parte cuchi-cuchi"

                      value1=s2ph(r)+s2ph(s)
                      value2=W2ph(r)+W2ph(s)+W22ph(r)+W22ph(s)
                      value3=s2hp(r)+s2hp(s)
                      value4=W2hp(r)+W2hp(s)+W22hp(r)+W22hp(s)
                      
                      ! Renormalized P3
                      
                      if (value1/=0.0_8) factors(r,1,6)=value2/value1 
                      if (value3/=0.0_8) factors(r,2,6)=value4/value3 
                      factors(r,3,6) = 0.0_8
                      factors(s,1,6) = factors(r,1,6)
                      factors(s,2,6) = factors(r,2,6)
                      factors(s,3,6) = factors(r,3,6)                      
                      factors2(r,1,6) = factors(r,1,6)
                      factors2(r,2,6) = factors(r,2,6)
                      factors2(r,3,6) = factors(r,3,6)                      
                      factors2(s,1,6) = factors(r,1,6)
                      factors2(s,2,6) = factors(r,2,6)
                      factors2(s,3,6) = factors(r,3,6)                      

                   else

!                      print *,"entro a la OTRA parte cuchi-cuchi"

                      value1=s2ph(r)+s2ph(s)
                      value2=W2ph(r)+W2ph(s)
                      value3=s2hp(r)+s2hp(s)
                      value4=W2hp(r)+W2hp(s)
                      
                      ! Renormalized P3
                      
                      if (value1/=0.0_8) factors(r,1,6)=value2/value1 
                      if (value3/=0.0_8) factors(r,2,6)=value4/value3 
                      factors(r,3,6) = 0.0_8
                      factors(s,1,6) = factors(r,1,6)
                      factors(s,2,6) = factors(r,2,6)
                      factors(s,3,6) = factors(r,3,6)

                      value1=s2ph(i)
                      value2=W22ph(r)+W22ph(s)
                      value3=s2hp(i)
                      value4=W22hp(r)+W22hp(s)
                      
                      ! Renormalized P3
                      
                      if (value1/=0.0_8) factors2(r,1,6)=value2/value1 
                      if (value3/=0.0_8) factors2(r,2,6)=value4/value3 
                      factors2(r,3,6) = 0.0_8
                      factors2(s,1,6) = factors2(r,1,6)
                      factors2(s,2,6) = factors2(r,2,6)
                      factors2(s,3,6) = factors2(r,3,6)

                   end if

                end if

             end if
             
             if (o==2) then
                
                do k = 1 , PropagatorTheory_instance%numberOfSpecies             
                   
                      if (k/=electrons(1).and.k/=electrons(2)) then

                         ! OVGF version A

                         value2=W2hp(k)+W2ph(k)
                         value1=s2hp(k)+s2ph(k)

                         if (value1/=0.0_8) factors(k,1,3) = value2/value1
                         factors(k,2,3) = factors(k,1,3)
                         factors(k,3,3) = factors(k,1,3)

                         value2=W22hp(k)+W22ph(k)
                         value1=s2hp(i)+s2ph(i)

                         if (value1/=0.0_8) factors2(k,1,3) = value2/value1
                         factors2(k,2,3) = factors2(k,1,3)
                         factors2(k,3,3) = factors2(k,1,3)
                         
                         ! OVGF version B

                         if (s2ph(k)/=0.0_8) factors(k,1,4) = W2ph(k)/s2ph(k)
                         if (s2hp(k)/=0.0_8) factors(k,2,4) = W2hp(k)/s2hp(k)
                         factors(k,3,4) = 0.0_8

                         if (s2ph(i)/=0.0_8) factors2(k,1,4) = W22ph(k)/s2ph(i)
                         if (s2hp(i)/=0.0_8) factors2(k,2,4) = W22hp(k)/s2hp(i)
                         factors2(k,3,4) = 0.0_8
                         
                         ! OVGF version C

                         value1=W2hp(k)+W2ph(k)+U2hp(k)+U2ph(k)
                         value2=factors(k,2,4)*(U2hp(k)+W2hp(k))+factors(k,1,4)*(U2ph(k)+W2ph(k))

                         if (value1/=0.0_8) factors(k,1,5) = value2/value1
                         factors(k,2,5) = factors(k,1,5)
                         factors(k,3,5) = factors(k,1,5)

                         value1=W22hp(k)+W22ph(k)+U22hp(k)+U22ph(k)
                         value2=factors2(k,2,4)*(U22hp(k)+W22hp(k))+factors2(k,1,4)*(U22ph(k)+W22ph(k))

                         if (value1/=0.0_8) factors2(k,1,5) = value2/value1
                         factors2(k,2,5) = factors2(k,1,5)
                         factors2(k,3,5) = factors2(k,1,5)

                      end if

                   end do

                   r = electrons(1)
                   s = electrons(2)
!                   print *,"values of r and s:",r,s

                   if (r/=0.and.s/=0) then

                      if (r==i.or.s==i) then

!                         print *,"entro a la parte cuchi-cuchi"

                         ! OVGF version A
                         value1=s2hp(r)+s2ph(r)+s2hp(s)+s2ph(s)
                         value2=W2hp(r)+W2ph(r)+W2hp(s)+W2ph(s)+W22hp(r)+W22ph(r)+W22hp(s)+W22ph(s)
                         
                         if (value1/=0.0_8) factors(r,1,3) = value2/value1
                         factors(r,2,3) = factors(r,1,3)
                         factors(r,3,3) = factors(r,1,3)
                         factors(s,1,3) = factors(r,1,3)
                         factors(s,2,3) = factors(r,1,3)
                         factors(s,3,3) = factors(r,1,3)
                         factors2(r,1,3) = factors(r,1,3)
                         factors2(r,2,3) = factors(r,1,3)
                         factors2(r,3,3) = factors(r,1,3)
                         factors2(s,1,3) = factors(r,1,3)
                         factors2(s,2,3) = factors(r,1,3)
                         factors2(s,3,3) = factors(r,1,3)

                         ! OVGF version B
                         
                         value1=s2ph(r)+s2ph(s)
                         value2=W2ph(r)+W2ph(s)+W22ph(r)+W22ph(s)
                         value3=s2hp(r)+s2hp(s)
                         value4=W2hp(r)+W2hp(s)+W22hp(r)+W22hp(s)
                         
                         if (value1/=0.0_8) factors(r,1,4)=value2/value1 
                         if (value3/=0.0_8) factors(r,2,4) =value4/value3 
                         factors(r,3,4) = 0.0_8
                         factors(s,1,4) = factors(r,1,4)
                         factors(s,2,4) = factors(r,2,4)
                         factors(s,3,4) = factors(r,3,4)
                         factors2(r,1,4) = factors(r,1,4)
                         factors2(r,2,4) = factors(r,2,4)
                         factors2(r,3,4) = factors(r,3,4)
                         factors2(s,1,4) = factors(r,1,4)
                         factors2(s,2,4) = factors(r,2,4)
                         factors2(s,3,4) = factors(r,3,4)

                         ! OVGF version C
                         value2=factors(r,2,4)*(U2hp(r)+W2hp(r)+U2hp(s)+W2hp(s)+U22hp(r)+W22hp(r)+U22hp(s)+W22hp(s))+&
                              factors(r,1,4)*(U2ph(r)+W2ph(r)+U2ph(s)+W2ph(s)+U22ph(r)+W22ph(r)+U22ph(s)+W22ph(s))
                         value1=W2hp(r)+W2ph(r)+U2hp(r)+U2ph(r)+W2hp(s)+W2ph(s)+U2hp(s)+U2ph(s)+&
                              U22ph(r)+W22ph(r)+U22ph(s)+W22ph(s)+U22hp(r)+W22hp(r)+U22hp(s)+W22hp(s)
                         
                         if (value1/=0.0_8) factors(r,1,5) = value2/value1
                         factors(r,2,5) = factors(r,1,5)
                         factors(r,3,5) = factors(r,1,5)
                         factors(s,1,5) = factors(r,1,5)
                         factors(s,2,5) = factors(r,2,5)
                         factors(s,3,5) = factors(r,3,5)
                         factors2(r,1,5) = factors(r,1,5)
                         factors2(r,2,5) = factors(r,2,5)
                         factors2(r,3,5) = factors(r,3,5)
                         factors2(s,1,5) = factors(r,1,5)
                         factors2(s,2,5) = factors(r,2,5)
                         factors2(s,3,5) = factors(r,3,5)

                      else

                         print *,"entro a la OTRA parte cuchi-cuchi"

                         ! OVGF version A
                         value1=s2hp(r)+s2ph(r)+s2hp(s)+s2ph(s)
                         value2=W2hp(r)+W2ph(r)+W2hp(s)+W2ph(s)
                         
                         if (value1/=0.0_8) factors(r,1,3) = value2/value1
                         factors(r,2,3) = factors(r,1,3)
                         factors(r,3,3) = factors(r,1,3)
                         factors(s,1,3) = factors(r,1,3)
                         factors(s,2,3) = factors(r,1,3)
                         factors(s,3,3) = factors(r,1,3)
                         
                         !
                         
                         value1=s2hp(i)+s2ph(i)
                         value2=W22hp(r)+W22ph(r)+W22hp(s)+W22ph(s)
                         
                         if (value1/=0.0_8) factors2(r,1,3) = value2/value1
                         factors2(r,2,3) = factors2(r,1,3)
                         factors2(r,3,3) = factors2(r,1,3)
                         factors2(s,1,3) = factors2(r,1,3)
                         factors2(s,2,3) = factors2(r,1,3)
                         factors2(s,3,3) = factors2(r,1,3)
                         
                         ! OVGF version B
                         
                         value1=s2ph(r)+s2ph(s)
                         value2=W2ph(r)+W2ph(s)
                         value3=s2hp(r)+s2hp(s)
                         value4=W2hp(r)+W2hp(s)
                         
                         if (value1/=0.0_8) factors(r,1,4)=value2/value1 
                         if (value3/=0.0_8) factors(r,2,4) =value4/value3 
                         factors(r,3,4) = 0.0_8
                         factors(s,1,4) = factors(r,1,4)
                         factors(s,2,4) = factors(r,2,4)
                         factors(s,3,4) = factors(r,3,4)
                         
                         !
                         
                         value1=s2ph(i)
                         value2=W22ph(r)+W22ph(s)
                         value3=s2hp(i)
                         value4=W22hp(r)+W22hp(s)
                         
                         if (value1/=0.0_8) factors2(r,1,4)=value2/value1 
                         if (value3/=0.0_8) factors2(r,2,4) =value4/value3 
                         factors2(r,3,4) = 0.0_8
                         factors2(s,1,4) = factors2(r,1,4)
                         factors2(s,2,4) = factors2(r,2,4)
                         factors2(s,3,4) = factors2(r,3,4)
                         
                         ! OVGF version C
                         value2=factors(r,2,4)*(U2hp(r)+W2hp(r)+U2hp(s)+W2hp(s))+&
                              factors(r,1,4)*(U2ph(r)+W2ph(r)+U2ph(s)+W2ph(s))
                         value1=W2hp(r)+W2ph(r)+U2hp(r)+U2ph(r)+W2hp(s)+W2ph(s)+U2hp(s)+U2ph(s)
                         
                         if (value1/=0.0_8) factors(r,1,5) = value2/value1
                         factors(r,2,5) = factors(r,1,5)
                         factors(r,3,5) = factors(r,1,5)
                         factors(s,1,5) = factors(r,1,5)
                         factors(s,2,5) = factors(r,2,5)
                         factors(s,3,5) = factors(r,3,5)
                         
                         !
                         
                         value2=factors2(r,2,4)*(U22hp(r)+W22hp(r)+U22hp(s)+W22hp(s))+&
                              factors2(r,1,4)*(U22ph(r)+W22ph(r)+U22ph(s)+W22ph(s))
                         value1=W22hp(r)+W22ph(r)+U22hp(r)+U22ph(r)+W22hp(s)+W22ph(s)+U22hp(s)+U22ph(s)
                         
                         if (value1/=0.0_8) factors2(r,1,5) = value2/value1
                         factors2(r,2,5) = factors2(r,1,5)
                         factors2(r,3,5) = factors2(r,1,5)
                         factors2(s,1,5) = factors2(r,1,5)
                         factors2(s,2,5) = factors2(r,2,5)
                         factors2(s,3,5) = factors2(r,3,5)

                      end if

                   end if

             end if
             
             poleStrenght = 1.0_8/(selfEnergyDerivative)
             thirdOrderResults(1,o) = 27.211396_8 * newOmega
             thirdOrderResults(2,o) = poleStrenght

!             print *,"value of o:",o
!             ! print *,"FACTORS:"
!             ! print *,factors(:,:,:)
             print *, "Constant self-energy: ", constantSelfEnergy*27.211396_8
             print *, "2hp(2):               ",s2hp(:)
             print *, "2ph(2):               ",s2ph(:)
             print *, "W 2hp(3):             ",W2hp(:)
             print *, "W 2ph(3):             ",W2ph(:)
             print *, "U 2hp(3):             ",U2hp(:)
             print *, "U 2ph(3):             ",U2ph(:)
             print *, "W2 2hp(3):            ",W22hp(:)
             print *, "W2 2ph(3):            ",W22ph(:)
             print *, "U2 2hp(3):            ",U22hp(:)
             print *, "U2 2ph(3):            ",U22ph(:)
             print *, "Factor 2hp:           ",factors(:,2,o)
             print *, "Factor 2ph:           ",factors(:,1,o)
             print *, "Factor2 2hp:          ",factors2(:,2,o)
             print *, "Factor2 2ph:         ",factors2(:,1,o)
             print *, ""
             write (*,"(T5,A10,A10,A6,F8.4,A7,I2,A12)") "Optimized ",thirdOrderMethods(o),"pole: ",newOmega*27.211396_8," after ",ni," iterations."
             write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
             print *,"----------------------------------------------------------------"
                          
          end do ! options for third order

          ! printing results for one spin-orbital

          write (*,"(T5,A55,I2,A13,A8)") "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL:",&
               int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1))," OF SPECIES:",nameOfSpeciesA 
          write (*, "(T5,A45)") "--------------------------------------------" 
          write (*, "(T10,A12,A12,A12)") " Method ","BE (eV)","Pole S."
          write (*, "(T5,A45)") "--------------------------------------------"
          write (*,"(T10,A12,F12.6)") "KT        ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
          write (*,"(T10,A12,F12.6,F12.6)") "EP2       ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
               PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)
          do o=1,6

             write (*,"(T10,A12,F12.6,F12.6)") thirdOrderMethods(o),thirdOrderResults(1,o),thirdOrderResults(2,o)
             
          end do
          write (*, "(T5,A45)") "--------------------------------------------"

          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,5)=27.211396_8 * newOmega
          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,6)=poleStrenght
                              
          ! call Matrix_destructor(auxMatrix2(:))          
!          call TransformIntegrals_destructor( repulsionTransformer )
          
       end do
       
    end do

    do i = 1, PropagatorTheory_instance%numberOfSpecies

!JC ?       print *,"Second order densities for species:",i
!           print *,secondOrderDensities(i)%values(:,:)

    end do
    
    !!
    !!************************************************************************************************
    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
    print *,"***************************************************************"

    close(wfnUnit)

  end subroutine PropagatorTheory_thirdOrderCorrection5
  
! 
!  !**
!  ! @brief Sum the energy corrections in the total propagator results
!  !**
  subroutine PropagatorTheory_run()
    implicit NONE
    
    type(Exception) :: ex
    integer :: i, j, nproc


    nproc = CONTROL_instance%NUMBER_OF_CORES
    
    if ( PropagatorTheory_instance%isInstanced ) then

       i = PropagatorTheory_instance%orderOfCorrection

       select case( i )
          
       case(  SECOND_ORDER )
          
          call PropagatorTheory_secondOrderCorrection()
          ! call PropagatorTheory_nonDiagonalSecondOrderTDACorrection()
          call Matrix_copyConstructor(PropagatorTheory_instance%energyCorrections , &
               PropagatorTheory_instance%energyCorrectionsOfSecondOrder)

       case ( THIRD_ORDER )

          call PropagatorTheory_thirdOrderCorrection5()
   
       case default
          
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object PropagatorTheory in run() function" )
          call Exception_setDescription( ex, "This order correction hasn't been implemented" )
          call Exception_show( ex )
          
       end select
       
    else
       
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object PropagatorTheory in run() function" )
       call Exception_setDescription( ex, "You should to instance PropagatorTheory module before use this function" )
       call Exception_show( ex )
       
    end if
    
    
  end subroutine PropagatorTheory_run
  
  
  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine PropagatorTheory_exception( typeMessage, description, debugDescription)
    implicit NONE
    integer :: typeMessage
    character(*) :: description
    character(*) :: debugDescription
    
    type(Exception) :: ex
    
    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )
    
  end subroutine PropagatorTheory_exception
  
end module PropagatorTheory_


  ! !**
  ! ! @brief Evaluate numerators and denominators of the self-energy expression and evaluate the corrected 
  ! ! koopmans energy at the second order level
  ! !**
  
  ! subroutine PropagatorTheory_secondOrderCorrection()
  !   implicit NONE
    
  !   integer :: a ! Index for orbital
  !   integer :: b ! Index for orbital
  !   integer :: r ! Index for orbital
  !   integer :: s ! Index for orbital
  !   integer :: p ! Index for orbital
  !   integer :: q ! Index for orbital
  !   integer :: t ! Index for orbital
  !   integer :: i ! Index to count species
  !   integer :: j ! Index to count species
  !   integer :: k ! Index to count orbitals
  !   integer :: id ! Index to count numerators and denominator
  !   integer :: id1 ! Index to count numerators and denominator
  !   integer :: id2 ! Index to count numerators and denominator
  !   integer :: id3 ! Index to count numerators and denominator
  !   integer :: ab ! Index to count numerators and denominators
  !   integer :: ac ! Index to count numerators and denominators
  !   integer :: m
  !   integer :: n
  !   integer :: u
  !   integer :: specieID
  !   integer :: specie1ID
  !   integer :: specie2ID
  !   integer :: otherSpecieID
  !   integer :: electronsID
  !   integer :: occupationNumber
  !   integer :: occupationNumberOfOtherSpecie
  !   integer :: vectorSize
  !   integer :: vectorSize2
  !   integer :: ocupation1
  !   integer :: numberOfContractions
  !   integer :: numberOfContractionsOfOtherSpecie
  !   integer(4) :: errorNum
  !   integer(8) :: auxIndex
  !   character(10) :: nameOfSpecie
  !   character(10) :: nameOfOtherSpecie
  !   type(Vector) :: eigenValues
  !   type(Vector) :: eigenValues1
  !   type(Vector) :: eigenValuesOfOtherSpecie
  !   type(Matrix) :: auxMatrix
  !   type(Matrix),allocatable :: auxMatrix2(:)
  !   type(Vector) :: auxNumeratorsVector  
  !   type(Vector) :: auxDenominatorsVector  
  !   type(Vector) :: auxNumeratorsVector2  
  !   type(Vector) :: auxDenominatorsVector2  
  !   type(Vector) :: auxNumeratorsVector3  
  !   type(Vector) :: auxDenominatorsVector3  
  !   type(Vector) :: observador
  !   type(Vector) :: occupations
  !   type(Vector) :: occupationsOfOtherSpecie
  !   real(8) :: lambda
  !   real(8) :: charge, chargeSpecie, chargeOtherSpecie
  !   real(8) :: lambdaOfOtherSpecie
  !   real(8) :: independentEnergyCorrection
  !   real(8) :: couplingEnergyCorrection
  !   real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
  !   real(8) :: auxVal_A, auxVal_B, auxVal_C, auxVal_D 
  !   real(8) :: auxVal_E, auxVal_F, auxVal_G, auxVal_H
  !   real(8) :: preliminary
  !   real(8) :: range
  !   real(8) :: correctionWithoutInterParticlesTerm 
  !   real(8) :: optimizedOmega 
  !   real(8) :: selfEnergyInterSpecie 
  !   real(8) :: selfEnergyIntraSpecie 
    
  !   !!*******************************************************************************************
  !   !! Determinate the numerators and denominators of the second Oder propapator 
  !   !!
  !   ! if ( .not.CONTROL_instance%OPTIMIZE ) then
  !   !    print *,"===================================================="
  !   !    print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
  !   !    print *,"===================================================="
  !   !    print *,"    Algorithm Four-index integral tranformation"
  !   !    print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
  !   !    print *,"  Computer Physics Communications, 2005, 166, 58-65 "
  !   !    print *,"--------------------------------------------------"
  !   !    print *,""
       
  !   ! end if
        
  !   print *,"******************************************************************"
  !   print *,"BEGINNING OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR CALCULATIONS"

  !   electronsID = MolecularSystem_getSpecieID(  nameOfSpecie="e-" )

  !   if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
       
  !      if (allocated(auxMatrix2)) deallocate(auxMatrix2)
  !      allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies))
       
  !   end if

  !   if (CONTROL_instance%PT_TRANSITION_OPERATOR.or.CONTROL_instance%PT_JUST_ONE_ORBITAL) then
  !      specie1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
  !      specie2ID= specie1ID
  !   else
  !      specie1ID=1
  !      specie2ID=PropagatorTheory_instance%numberOfSpecies
  !   end if

  !   do i = specie1ID , specie2ID
       
  !      nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
       
  !      specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
  !      charge = MolecularSystem_getCharge( specieID )
  !      eigenValues = MolecularSystem_getEigenValues(i)
  !      occupationNumber = MolecularSystem_getOcupationNumber( i )
  !      numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
  !      lambda = MolecularSystem_getLambda( i )

  !      if (CONTROL_instance%PT_TRANSITION_OPERATOR.or.CONTROL_instance%PT_JUST_ONE_ORBITAL) then
  !         PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
  !         PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
  !      else
  !         PropagatorTheory_instance%virtualBoundary= occupationNumber
  !         PropagatorTheory_instance%occupationBoundary = 1
  !      end if
       
  !      vectorSize=0
       
  !      do j = 1 , PropagatorTheory_instance%numberOfSpecies
          
  !         occupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
  !         numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
  !         vectorSize=vectorSize+(numberOfContractions-occupationNumber)* & 
  !              (numberOfContractionsOfOtherSpecie-occupationNumberOfOtherSpecie)*occupationNumberOfOtherSpecie &
  !              + (occupationNumberOfOtherSpecie)*(occupationNumber)* &
  !              (numberOfContractionsOfOtherSpecie-occupationNumberOfOtherSpecie)
  !         if (nameOfSpecie == CONTROL_instance%IONIZE_SPECIE) then
  !            vectorSize=vectorSize+(numberOfContractionsOfOtherSpecie-occupationNumberOfOtherSpecie)*occupationNumberOfOtherSpecie
  !         end if
  !      end do
       
  !      if (CONTROL_instance%PT_TRANSITION_OPERATOR) then
  !         call Vector_constructor(occupations,occupationNumber)
  !         do k=1, occupationNumber
  !            if (k==CONTROL_instance%IONIZE_MO) then
  !               occupations%values(k)=CONTROL_instance%MO_FRACTION_OCCUPATION
  !            else
  !               occupations%values(k)=1.0_8
  !            end if
  !         end do
  !      end if
       
  !      !!**************************************************************************
  !      !!	Storing of denominators and numerators in the corresponding vectors
  !      !!****

  !      m=0 !orbital counter
  !      vectorSize2=0

  !      do p=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	
  !         id=0
  !         m=m+1
          
  !         call Vector_constructor(auxDenominatorsVector,vectorSize)
  !         call Vector_constructor(auxNumeratorsVector,vectorSize)
          
  !         do a=1, occupationNumber
  !            do b=1, occupationNumber
  !               do r=occupationNumber+1, numberOfContractions
  !                  id=id+1

  !                  auxVal_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, a, r, b )
  !                  auxVal_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, b, a, r )
  !                  auxVal_C = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, a, b, b )
  !                  auxVal_D = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, b, b, a )
  !                  auxVal_E = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, a, r, r )
  !                  auxVal_F = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, r, r, a )
  !                  auxVal_G = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, b, b, r, r )
  !                  auxVal_H = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, b, r, r, b )
  !                  ! print *, p, a, b, r, " auxVal_A:", auxVal_A, "auxVal_B:", auxVal_B

  !                  if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfSpecie == CONTROL_instance%IONIZE_SPECIE) then

  !                     auxNumeratorsVector%values(id)= ( occupations%values(a)* &
  !                          occupations%values(b))*(auxVal_A  &
  !                          * ( lambda * auxVal_A  - auxVal_B ) * (charge**4.0_8 ))
  !                     auxDenominatorsVector%values(id)=  ( eigenValues%values(r) &
  !                          - eigenValues%values(a) - eigenValues%values(b)  &
  !                    + (0.5_8 *  auxVal_C * ( lambda * auxVal_C  - auxVal_D ) * (charge**4.0_8 )) &
  !                    - (auxVal_E * ( lambda * auxVal_E  - auxVal_F ) * (charge**4.0_8 )) &
  !                    - (auxVal_G * ( lambda * auxVal_G  - auxVal_H ) * (charge**4.0_8 )) )

             
  !                  else
                      
  !                     auxNumeratorsVector%values(id)= (auxVal_A  &
  !                          * ( lambda * auxVal_A  - auxVal_B ) * (charge**4.0_8 ))
  !                     auxDenominatorsVector%values(id)=  ( eigenValues%values(r) &
  !                          - eigenValues%values(a) - eigenValues%values(b) &
  !                    + (0.5_8 *  auxVal_C * ( lambda * auxVal_C  - auxVal_D ) * (charge**4.0_8 )) &
  !                    - (auxVal_E * ( lambda * auxVal_E  - auxVal_F ) * (charge**4.0_8 )) &
  !                    - (auxVal_G * ( lambda * auxVal_G  - auxVal_H ) * (charge**4.0_8 )) )
                           

  !                  end if
                   
  !               end do
  !            end do
  !         end do
  !         do a=1, occupationNumber
  !            do r=occupationNumber+1, numberOfContractions
  !               do s=occupationNumber+1, numberOfContractions
  !                  id=id+1

  !                  auxVal_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, r, a, s )
  !                  auxVal_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, s, r, a )
  !                  auxVal_C = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, a, r, r )
  !                  auxVal_D = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, r, r, a )
  !                  auxVal_E = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, a, s, s )
  !                  auxVal_F = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, s, s, a )
  !                  auxVal_G = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, r, r, s, s )
  !                  auxVal_H = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, r, s, s, r )
  !                  ! print *, p, a, r, s, " auxVal_A:", auxVal_A, "auxVal_B:", auxVal_B
                   
  !                  auxNumeratorsVector%values(id)= auxVal_A  &
  !                       * ( lambda * auxVal_A  - auxVal_B ) * ( (charge)**4.0_8 )
  !                  auxDenominatorsVector%values(id)=  ( eigenValues%values(a) &
  !                       - eigenValues%values(r) - eigenValues%values(s) &
  !                    + (0.5_8 *  auxVal_C * ( lambda * auxVal_C  - auxVal_D ) * (charge**4.0_8 )) &
  !                    + (auxVal_E * ( lambda * auxVal_E  - auxVal_F ) * (charge**4.0_8 )) &
  !                    - (auxVal_G * ( lambda * auxVal_G  - auxVal_H ) * (charge**4.0_8 )) )

                   
  !               end do
  !            end do
  !         end do
          
  !         if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfSpecie == CONTROL_instance%IONIZE_SPECIE .and. p==CONTROL_instance%IONIZE_MO) then
  !            do a=1, occupationNumber
  !               do s=occupationNumber+1, numberOfContractions
  !                  id=id+1

  !                  auxVal_A = TransformIntegrals2_atomicToMolecularOfOneSpecieTOP2A( nameOfSpecie, p, p, a, s )
  !                  auxVal_B = TransformIntegrals2_atomicToMolecularOfOneSpecieTOP2B( nameOfSpecie, p, a, s, p )
  !                  auxVal_C = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, p, a, a )
  !                  auxVal_D = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, a, a, p )
  !                  auxVal_E = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, p, s, s )
  !                  auxVal_F = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, p, s, s, p )
  !                  auxVal_G = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, a, s, s )
  !                  auxVal_H  = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecie, a, s, s, a )
  !                  ! print *, p, a, p, s, " auxVal_A:", auxVal_A, "auxVal_B:", auxVal_B

  !                  auxNumeratorsVector%values(id)= (1.0_8-occupations%values(p))* auxVal_A  &
  !                       * ( lambda * auxVal_A  - auxVal_B ) * ( (charge)**4.0_8 )
  !                  auxDenominatorsVector%values(id)=  ( eigenValues%values(a) &
  !                       - eigenValues%values(p) - eigenValues%values(s) &
  !                    + (1.0_8-occupations%values(p))*(-0.5_8*(auxVal_C * ( lambda * auxVal_C  - auxVal_D ) &
  !                  * (charge**4.0_8 )) + (auxVal_E * ( lambda * auxVal_E  - auxVal_F ) * (charge**4.0_8 ))) &
  !                    + ( auxVal_G * ( lambda * auxVal_G  - auxVal_H ) * (charge**4.0_8 )) )
                   
  !               end do
  !            end do
  !         end if
          
  !         preliminary=eigenValues%values(p)
  !         id1=id
  !         vectorSize2=id
          
  !         if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
             
  !            do j = 1 , PropagatorTheory_instance%numberOfSpecies
  !               if (j.ne.i) then
                   
  !                  nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( j ) )
  !                  otherSpecieID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
  !                  eigenValuesOfOtherSpecie = MolecularSystem_getEigenValues(j)
  !                  occupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
  !                  numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
  !                  lambdaOfOtherSpecie = MolecularSystem_getLambda( j )

  !                  if ( (nameOfSpecie=="e-ALPHA" .and. nameOfOtherSpecie=="e-BETA") .or. (nameOfSpecie=="e-BETA" .and. nameOfOtherSpecie=="e-ALPHA") ) id2=id
 
  !                  if (nameOfOtherSpecie == CONTROL_instance%IONIZE_SPECIE) then
  !                     call Vector_constructor(occupationsOfOtherSpecie,occupationNumberOfOtherSpecie)
  !                     do k=1, occupationNumberOfOtherSpecie
  !                        if (k==CONTROL_instance%IONIZE_MO) then
  !                           occupationsOfOtherSpecie%values(k)=CONTROL_instance%MO_FRACTION_OCCUPATION
  !                        else
  !                           occupationsOfOtherSpecie%values(k)=1.0_8
  !                        end if
  !                     end do
  !                  end if

  !                  chargeSpecie = MolecularSystem_getCharge( specieID )
  !                  chargeOtherSpecie = MolecularSystem_getCharge( otherSpecieID )
                   
  !                  do a=1, occupationNumberOfOtherSpecie
  !                     do r=occupationNumber+1, numberOfContractions
  !                        do s=occupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie
  !                           id=id+1

  !                           auxVal = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, p, r, a, s )
  !                           auxVal_1 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, a, a, r, r )
  !                           auxVal_2 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, a, a, s, s )
  !                           auxVal_3 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, r, r, s, s )
  !                           auxVal = auxVal * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_1 = auxVal_1 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_2 = auxVal_2 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_3 = auxVal_3 * chargeSpecie * chargeOtherSpecie 
                            
  !                           ! print *, p, a, r, s, " auxVal", auxVal

  !                           if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfOtherSpecie == CONTROL_instance%IONIZE_SPECIE) then 
  !                              auxNumeratorsVector%values(id)= (occupationsOfOtherSpecie%values(a))*( (auxVal)**2.0_8 ) &
  !                                   * lambda * lambdaOfOtherSpecie
  !                              auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(a) &
  !                                   -eigenValues%values(r)-eigenValuesOfOtherSpecie%values(s) + auxVal_1 + auxVal_2 - 0.5_8*auxVal_3  )
                                                              
  !                           else
                               
  !                              auxNumeratorsVector%values(id)= ( (auxVal)**2.0_8 ) &
  !                                   * lambda * lambdaOfOtherSpecie
  !                              auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(a) &
  !                                   -eigenValues%values(r)-eigenValuesOfOtherSpecie%values(s) + auxVal_1 + auxVal_2 - 0.5_8*auxVal_3 )
                                                              
  !                           end if
                            
  !                        end do
  !                     end do
  !                  end do
                   
  !                  do r=occupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie
  !                     do a=1, occupationNumber
  !                        do b=1, occupationNumberOfOtherSpecie 
  !                           id=id+1

  !                           auxVal =  TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, p, a, r, b )
  !                           auxVal_1 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, a, a, b, b )
  !                           auxVal_2 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, a, a, r, r )
  !                           auxVal_3 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, b, b, r, r )
  !                           auxVal = auxVal * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_1 = auxVal_1 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_2 = auxVal_2 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_3 = auxVal_3 * chargeSpecie * chargeOtherSpecie 

  !                           ! print *, p, r, a, b, " auxVal", auxVal

  !                          if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfSpecie == CONTROL_instance%IONIZE_SPECIE) then   

  !                              auxNumeratorsVector%values(id)= (occupations%values(a))*(( auxVal )**2.0_8 ) &
  !                                * lambda * lambdaOfOtherSpecie
  !                              auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(r) &
  !                                   -eigenValues%values(a)-eigenValuesOfOtherSpecie%values(b) + 0.5_8*auxVal_1 - auxVal_2 - auxVal_3  )

  !                           else if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfOtherSpecie == CONTROL_instance%IONIZE_SPECIE) then
                               
  !                              auxNumeratorsVector%values(id)= (occupationsOfOtherSpecie%values(b))*(( auxVal )**2.0_8 ) &
  !                                * lambda * lambdaOfOtherSpecie
  !                              auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(r) &
  !                                   -eigenValues%values(a)-eigenValuesOfOtherSpecie%values(b) + 0.5_8*auxVal_1 - auxVal_2 - auxVal_3)

  !                           else 

  !                              auxNumeratorsVector%values(id)= (( auxVal)**2.0_8 ) &
  !                                * lambda * lambdaOfOtherSpecie
  !                              auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(r) &
  !                                   -eigenValues%values(a)-eigenValuesOfOtherSpecie%values(b) + 0.5_8*auxVal_1 - auxVal_2 - auxVal_3 )
                               
  !                           end if
                            
  !                        end do
  !                     end do
  !                  end do
                   
  !                  if (CONTROL_instance%PT_TRANSITION_OPERATOR .and. nameOfSpecie == CONTROL_instance%IONIZE_SPECIE .and. p==CONTROL_instance%IONIZE_MO) then
  !                     do a=1, occupationNumberOfOtherSpecie
  !                        do s=occupationNumberOfOtherSpecie+1, numberOfContractionsOfOtherSpecie
                            
  !                           id=id+1

  !                           auxVal = TransformIntegrals2_atomicToMolecularOfTwoSpeciesTOP2 ( nameOfSpecie , nameOfOtherSpecie, p, p, a, s )
  !                           auxVal_1 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, p, p, a, a )
  !                           auxVal_2 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, p, p, r, r )
  !                           auxVal_3 = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecie , nameOfOtherSpecie, a, a, r, r )
  !                           auxVal = auxVal * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_1 = auxVal_1 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_2 = auxVal_2 * chargeSpecie * chargeOtherSpecie 
  !                           auxVal_3 = auxVal_3 * chargeSpecie * chargeOtherSpecie 
  !                           ! print *, p, p, a, s, " auxVal", auxVal

  !                           auxNumeratorsVector%values(id)= (1.0_8-occupations%values(p))*( ( auxVal )**2.0_8 ) &
  !                                * lambda * lambdaOfOtherSpecie
  !                           auxDenominatorsVector%values(id)= ( eigenValuesOfOtherSpecie%values(a) &
  !                                -eigenValues%values(p)-eigenValuesOfOtherSpecie%values(s) &
  !                                + (1.0_8-occupations%values(p))*(-0.5_8*auxVal_1 + auxVal_2) + auxVal_3 )
                            
  !                        end do
  !                     end do
  !                  end if

  !                  if ( (nameOfSpecie=="e-ALPHA" .and. nameOfOtherSpecie=="e-BETA") .or. (nameOfSpecie=="e-BETA" .and. nameOfOtherSpecie=="e-ALPHA") ) then
                      
  !                     id3=id
  !                     vectorSize2=vectorSize2+(id3-id2)

  !                  end if
                   
  !               end if

  !            end do
  !         end if
          
  !         call Vector_constructor(auxDenominatorsVector2,vectorSize2)
  !         call Vector_constructor(auxNumeratorsVector2,vectorSize2)

  !         auxDenominatorsVector2%values(1:vectorSize2)=auxDenominatorsVector%values(1:vectorSize2)
  !         auxNumeratorsVector2%values(1:vectorSize2)=auxNumeratorsVector%values(1:vectorSize2)
          
  !         ! if ( (nameOfSpecie=="e-ALPHA" .and. nameOfOtherSpecie=="e-BETA") .or. (nameOfSpecie=="e-BETA" .and. nameOfOtherSpecie=="e-ALPHA") ) then
             
  !         !    auxDenominatorsVector2%values((id1+1):vectorSize2)=auxDenominatorsVector%values(id2+1:id3)
  !         !    auxNumeratorsVector2%values((id1+1):vectorSize2)=auxNumeratorsVector%values(id2+1:id3)
             
  !         ! end if
          
  !         ! print *,"Numerator a-a",auxNumeratorsVector2%values
  !         ! print *,"Numerator a-b",auxNumeratorsVector3%values
          
  !         ! Calculation of optimized omega

  !         print *,"Iteration scheme employed:",CONTROL_instance%PT_ITERATION_SCHEME,"NR Default"
  !         print *,""
          
  !         optimizedOmega=PropagatorTheory_evaluateSelfEnergy( auxNumeratorsVector, &
  !              auxDenominatorsVector, id, eigenValues%values(p), eigenValues%values(p))
          
  !         print *,"size1",vectorSize,"size2",vectorSize2,"id1",id1,"id2",id2,"id3",id3
          
  !         selfEnergyIntraSpecie=PropagatorTheory_getSelfEnergy(auxNumeratorsVector2, auxDenominatorsVector2, vectorSize2, optimizedOmega)
          
  !         selfEnergyInterSpecie=optimizedOmega-selfEnergyIntraSpecie-eigenValues%values(p)
          
  !         print *,"IntraSpecie term:",selfEnergyIntraSpecie,"(a.u)"
  !         print *,"InterSpecie term:",selfEnergyInterSpecie,"(a.u)"
          
  !         ! Storing of corrections

  !         PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,4*m)= &
  !              27.211396_8 * optimizedOmega
          
  !         PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-1))= &
  !              PropagatorTheory_getPolarStrength( auxNumeratorsVector, &
  !              auxDenominatorsVector, id, optimizedOmega)

  !         PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-2))= 27.211396_8 * selfEnergyIntraSpecie
  !         PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-3))= 27.211396_8 * selfEnergyInterSpecie
          
  !         if ( CONTROL_instance%PT_SELF_ENERGY_SCAN ) then    
             
  !            range=(CONTROL_instance%PT_SELF_ENERGY_RANGE)/27.211396_8
             
  !            call PropagatorTheory_scanSelfEnergy( auxNumeratorsVector, &
  !                 auxDenominatorsVector, id, p, nameOfSpecie, & 
  !                 (range+eigenValues%values(p)), (-range+eigenValues%values(p)), &
  !                 eigenValues%values(p)) 
             
  !         end if
          
  !         call Vector_destructor(auxDenominatorsVector)
  !         call Vector_destructor(auxNumeratorsVector)		 		  
  !         call Vector_destructor(auxNumeratorsVector2)
  !         call Vector_destructor(auxDenominatorsVector2)            
          
  !      end do
  !      call Matrix_destructor(auxMatrix)          
  !   end do
  !   !!
  !   !!************************************************************************************************
  !   print *,"END OF CALCULATION OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR"
  !   print *,"**************************************************************"
  ! end subroutine PropagatorTheory_secondOrderCorrection

  !!! This procedure calculate the superoperator Hamiltonian Matrix of Second order and the eigenvectors and eigenvalues corresponding
  !!! to the requested ionization energies

  ! subroutine PropagatorTheory_nonDiagonalSecondOrderCorrection()
  !   implicit NONE
    
  !   integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha species
  !   integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta species
  !   integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha species
  !   integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta species
  !   integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha species
  !   integer :: pb, qb, rb, sb ! Indices for general orbitals of beta species
  !   integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
  !   integer :: i, j ! counters for species
  !   integer :: m, n, o ! auxiliar counters
  !   integer :: speciesID
  !   integer :: otherSpeciesID
  !   integer :: electronsID
  !   integer :: occupationNumberOfSpecies, virtualNumberOfSpecies
  !   integer :: occupationNumberOfOtherSpecies, virtualNumberOfOtherSpecies
  !   integer :: numberOfContractionsOfSpecies
  !   integer :: numberOfContractionsOfOtherSpecies
  !   integer :: vectorSize, vectorSizeaHf
  !   integer(4) :: errorNum
  !   integer(8) :: auxIndex, HamiltonianSize
  !   character(10) :: nameOfSpecies
  !   character(10) :: nameOfOtherSpecies
  !   type(Vector) :: occupationsOfOtherSpecies
  !   type(Vector) :: eigenValuesOfSpecies
  !   type(Vector) :: eigenValuesOfOtherSpecies
  !   real(8) :: lambdaOfSpecies, lambdaOfOtherSpecies 
  !   real(8) :: chargeOfSpecies, chargeOfOtherSpecies
  !   real(8) :: auxValue, auxValue_A, auxValue_B, auxValue_1, auxValue_2
  !   real(8) :: lastEigenvalue, newEigenvalue, norm
  !   type(Vector) :: superEigenvalues, subEigenvalues, auxVector, W, P 
  !   type(Vector) :: lastEigenvector, newEigenvector, residual, direction    
  !   type(Matrix) :: aHa, fHf, subEigenvectors
  !   type(Matrix) :: aHf, fHa, superHamiltonian, superEigenvectors
  !   type(Matrix) :: subHamiltonian, zeroH, inverseZeroH
  !   integer :: position
  !   real(8) :: u, t, v, remainder, eW, eHo, difference

  !   if ( .not.CONTROL_instance%OPTIMIZE ) then
  !      print *,"================================================================================="
  !      print *,"  If you are running propagator calculations the following papers must be cited: "
  !      print *,"================================================================================="
  !      print *," J. Chem. Phys. 137, 074105, (2012). J. Chem. Phys. 108, 1008, (1998).           "
  !      print *," J. Chem. Phys. 99, 6716, (1993).                                                "
  !      print *," Computational Chemistry: Reviews of Current Trends. Vol 2. Pag 1-61. 1997       "
  !      print *,"================================================================================="
  !   end if
    
  !   print *,"******************************************************************"
  !   print *,"BEGINNING OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR CALCULATIONS"

  !   speciesID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )

  !   nameOfSpecies= trim(  MolecularSystem_getNameOfSpecie( speciesID ) )
       
  !   chargeOfSpecies = MolecularSystem_getCharge( speciesID )
    
  !   eigenValuesOfSpecies = MolecularSystem_getEigenValues( speciesID )

  !   occupationNumberOfSpecies = MolecularSystem_getOcupationNumber( speciesID )
       
  !   numberOfContractionsOfSpecies = MolecularSystem_getTotalNumberOfContractions( speciesID )
    
  !   virtualNumberOfSpecies = numberOfContractionsOfSpecies - occupationNumberOfSpecies 
 
  !   lambdaOfSpecies = MolecularSystem_getLambda( speciesID )

  !   if ( CONTROL_instance%PT_JUST_ONE_ORBITAL) then
  !      PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
  !      PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
  !   else
  !      PropagatorTheory_instance%virtualBoundary= occupationNumberOfSpecies
  !      PropagatorTheory_instance%occupationBoundary = 1
  !   end if

  !   !!! Construction of the superoperator Hamiltonian matrix

  !   vectorSizeaHf=0 !!! length of the (a|fH) block

  !   do j = 1 , PropagatorTheory_instance%numberOfSpecies

  !         occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
  !         numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
  !         virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
          
  !      if ( j == speciesID ) then !!! Here, m*(n-1)*n/2
          
  !         vectorSizeaHf = vectorSizeaHf +  0.5*( (virtualNumberOfSpecies - 1)*occupationNumberOfSpecies*& 
  !              virtualNumberOfSpecies &
  !              +  (occupationNumberOfSpecies - 1) * occupationNumberOfSpecies*&
  !              virtualNumberOfOtherSpecies )
          
  !      else
          
  !         vectorSizeaHf = vectorSizeaHf + (virtualNumberOfOtherSpecies*occupationNumberOfOtherSpecies*& 
  !              virtualNumberOfSpecies) &
  !              +  (occupationNumberOfOtherSpecies* occupationNumberOfSpecies*&
  !              virtualNumberOfOtherSpecies)
          
  !      end if
       
  !   end do

  !   print *,"vectorSizeaHf:",vectorSizeaHf

  !   !!! Construction of the superoperator Hamiltonian

  !   HamiltonianSize = numberOfContractionsOfSpecies + vectorSizeaHf

  !   print *,"HamiltonianSize:",HamiltonianSize

  !   call Matrix_constructor(superHamiltonian, HamiltonianSize, HamiltonianSize)

  !   !!! Section (a|Ha) !!! Diagonal terms of this block are stored into aHa matrix

  !   do pa = 1, numberOfContractionsOfSpecies

  !      superHamiltonian%values(pa,pa) = eigenValuesOfSpecies%values(pa)

  !   end do

  !   ! aHa = MolecularSystem_getEigenValues( speciesID )

  !   !!! Section (f|Hf) !!! This part is the bottleneck for other orders, but for the second order ... well it is too easy
      
  !   !!!	Storing (f|Hf) terms

  !   idfHf = numberOfContractionsOfSpecies

  !   do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !      do ia = 1 , occupationNumberOfSpecies
  !         do ja = ia + 1 , occupationNumberOfSpecies
             
  !            idfHf = idfHf + 1
  !            print*,"idfHf:",idfHf
   
  !            superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(ia) + eigenValuesOfSpecies%values(ja) &
  !                                 - eigenValuesOfSpecies%values(aa) )

  !         end do
  !      end do
  !   end do

  !   do ia = 1 , occupationNumberOfSpecies
  !      do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !         do ba = aa + 1 , numberOfContractionsOfSpecies

  !            idfHf = idfHf + 1
  !            print*,"idfHf:",idfHf
             
  !            superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(aa) + eigenValuesOfSpecies%values(ba) &
  !                                 - eigenValuesOfSpecies%values(ia) )

  !         end do
  !      end do
  !   end do

  !   if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
             
  !      do j = 1 , PropagatorTheory_instance%numberOfSpecies
  !         if (j.ne.speciesID) then
             
  !            nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( j ) )
  !            otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
  !            eigenValuesOfOtherSpecies = MolecularSystem_getEigenValues(j)
  !            occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
  !            numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
  !            lambdaOfOtherSpecies = MolecularSystem_getLambda( j )
  !            virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies

  !            do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
  !               do ia = 1 , occupationNumberOfSpecies 
  !                 do ib = 1 , occupationNumberOfOtherSpecies

  !                     idfHf = idfHf + 1
  !                     print*,"idfHf:",idfHf
                      
  !                     superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(ia) + eigenValuesOfOtherSpecies%values(ib) &
  !                          - eigenValuesOfOtherSpecies%values(ab) )
                      
  !                  end do
  !               end do
  !            end do
             
  !            do ib = 1 , occupationNumberOfOtherSpecies
  !               do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !                  do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
                      
  !                     idfHf = idfHf + 1
  !                     print*,"idfHf:",idfHf

  !                     superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(aa) + eigenValuesOfOtherSpecies%values(ab) &
  !                          - eigenValuesOfOtherSpecies%values(ib) )
                      
  !                  end do
  !               end do
  !            end do
             
  !         end if
  !      end do
       
  !   end if

  !   !!! Building the (a|Hf) blocks
    
  !   ! call Matrix_constructor(aHf, vectorSizefHf, numberOfContractionsOfSpecies)

  !   do pa = 1 , numberOfContractionsOfSpecies

  !      idaHf = numberOfContractionsOfSpecies
       
  !      do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !         do ia = 1 , occupationNumberOfSpecies
  !            do ja = ia + 1 , occupationNumberOfSpecies
                
  !               idaHf = idaHf + 1

  !               auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ia, aa, ja )
  !               auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ja, ia, aa )
  !               ! print *,"auxValue_A:",auxValue_A,"auxValue_B:",auxValue_B
  !               ! print*,"idaHf:",idaHf

  !               superHamiltonian%values(idaHf, pa) = (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B)
  !               superHamiltonian%values(pa, idaHf) = superHamiltonian%values(idaHf, pa)

  !               print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)
                
  !            end do
  !         end do
  !      end do
       
  !      do ia = 1 , occupationNumberOfSpecies
  !         do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !            do ba = aa + 1 , numberOfContractionsOfSpecies
                
  !               idaHf = idaHf + 1

  !               auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, aa, ia, ba )
  !               auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ba, aa, ia )
  !               ! print *,"auxValue_A:",auxValue_A,"auxValue_B:",auxValue_B
  !               ! print*,"idaHf:",idaHf

  !               superHamiltonian%values( idaHf, pa) = (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B)
  !               superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

  !               print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)

  !            end do
  !         end do
  !      end do

  !      if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
          
  !         do j = 1 , PropagatorTheory_instance%numberOfSpecies
  !            if (j.ne.speciesID) then
                
  !               nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( j ) )
  !               otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
  !               eigenValuesOfOtherSpecies = MolecularSystem_getEigenValues(j)
  !               occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
  !               numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
  !               lambdaOfOtherSpecies = MolecularSystem_getLambda( j )
  !               virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
  !               chargeOfOtherSpecies = MolecularSystem_getCharge( j )
                
  !               do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
  !                  do ia = 1 , occupationNumberOfSpecies
  !                     do ib = 1 , occupationNumberOfOtherSpecies
                         
  !                        idaHf = idaHf + 1

  !                        auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, pa, ia, ab, ib )
  !                        ! print *,"auxValue:",auxValue
  !                        ! print*,"idaHf:",idaHf

  !                        superHamiltonian%values( idaHf, pa) = auxValue * chargeOfSpecies* chargeOfOtherSpecies !!! OJO lambdas
  !                        superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

  !                        print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)

  !                     end do
  !                  end do
  !               end do
                
  !               do ib = 1 , occupationNumberOfOtherSpecies
  !                  do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
  !                     do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
                         
  !                        idaHf = idaHf + 1

  !                        auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, pa, aa, ib, ab )     
  !                        ! print *,"auxValue:",auxValue
  !                        ! print*,"idaHf:",idaHf

  !                        superHamiltonian%values( idaHf, pa) = auxValue * chargeOfSpecies* chargeOfOtherSpecies !!! OJO lambdas
  !                        superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

  !                        print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)                         

  !                     end do
  !                  end do
  !               end do
                
  !            end if
  !         end do
          
  !      end if
     
  !   end do

  !   print *,"idaHf:",idaHf,"idfHf:",idfHf,"HamiltonianSize:",HamiltonianSize

  !   !!! HERE THE CALCULATION OF OPTIMAL ENERGIES START

  !   !!! DAVIDSON

  !   call Vector_constructor( lastEigenvector, int(HamiltonianSize) )
  !   call Vector_constructor( newEigenvector, int(HamiltonianSize) )
  !   call Vector_constructor( residual, int(HamiltonianSize) )
  !   call Vector_constructor( P, int(HamiltonianSize) )
  !   call Vector_constructor( direction, int(HamiltonianSize) )
  !   call Vector_constructor( auxVector, int(HamiltonianSize) )
  !   call Vector_constructor( subEigenvalues, 2 )
  !   call Matrix_constructor( subHamiltonian, 2, 2 )
  !   call Matrix_constructor( subEigenvectors, 2, 2 )
  !   call Vector_constructor( superEigenvalues, int(HamiltonianSize) )
  !   call Matrix_constructor( superEigenvectors, HamiltonianSize, HamiltonianSize )

  !   !!! Initialization of the algorithm

  !   do pa = 1, HamiltonianSize
       
  !      if ( pa==CONTROL_instance%IONIZE_MO ) then
          
  !         lastEigenvector%values(pa) = 1.0_8       
          
  !      else

  !         lastEigenvector%values(pa) = 0.0_8      
    
  !      end if

  !   end do

  !   lastEigenvalue = eigenValuesOfSpecies%values(CONTROL_instance%IONIZE_MO)

  !   print *,lastEigenvector%values

  !   !!! DAVIDSON ALGORYTHM STARTS

  !   ! !!! loop starts here, on n, C must be normalized at this point
  !   n = 0
  !   difference = 1.0_8

  !   do while ( difference > 0.00001_8 )

  !      do pa = 1, HamiltonianSize
          
  !         superHamiltonian%values(pa,pa) = superHamiltonian%values(pa,pa) - lastEigenvalue 
  !         print *,"diagonal element of H:",pa,"value:",superHamiltonian%values(pa,pa)
          
  !      end do

  !      n = n+1
  !      !!! Substraction of En1
              
  !      !!! Calculating residual

  !      residual%values = matmul(superHamiltonian%values,lastEigenvector%values)

  !      !!!! Calculating new direction
       
  !      do pa = 1, numberOfContractionsOfSpecies
          
  !         !!! Be careful in the first iterations

  !         if ( abs(superHamiltonian%values(pa,pa)) < 0.0000001_8  ) then

  !            P%values(pa) = 0.0_8
  !            print *,"ENTRO:",pa

  !         else

  !            P%values(pa) = residual%values(pa)/superHamiltonian%values(pa,pa) 

  !         end if
          
  !      end do

  !      do pa = numberOfContractionsOfSpecies+1, HamiltonianSize
          
  !            P%values(pa) = -residual%values(pa)/lastEigenvalue
          
  !      end do

  !      ! Orthogonalize again C

  !      direction%values = P%values - (dot_product(P%values,lastEigenvector%values)/dot_product(lastEigenvector%values,lastEigenvector%values))*lastEigenvector%values

  !      !!! Normalization of the new direction
       
  !      norm = Vector_norm( direction )
       
  !      direction%values = direction%values/norm

  !      print *,"dot product of C and d:",dot_product(direction%values,lastEigenvector%values)

  !   !!! Constructing superHamiltonian in the basis of C and d

  !      do pa = 1, HamiltonianSize
          
  !         superHamiltonian%values(pa,pa) = superHamiltonian%values(pa,pa) + lastEigenvalue 
          
  !      end do

  !      subHamiltonian%values(1,1) = lastEigenvalue
  
  !      subHamiltonian%values(1,2) = dot_product( lastEigenvector%values, matmul(superHamiltonian%values,direction%values) )
       
  !      subHamiltonian%values(2,1) = dot_product( direction%values, matmul(superHamiltonian%values, lastEigenvector%values) )
       
  !      subHamiltonian%values(2,2) = dot_product( direction%values, matmul(superHamiltonian%values, direction%values) )

  !      call Matrix_eigen( subHamiltonian, subEigenvalues, subEigenvectors, SYMMETRIC )

  !      if (abs(lastEigenvalue-subEigenvalues%values(1))<abs(lastEigenvalue-subEigenvalues%values(2))) then

  !         newEigenvalue=subEigenvalues%values(1)
  !         norm=1.0_8/sqrt( dot_product( subEigenvectors%values(:,1), subEigenvectors%values(:,1) ) )
  !         newEigenvector%values=norm*subEigenvectors%values(1,1)*lastEigenvector%values + &
  !                               norm*subEigenvectors%values(2,1)*direction%values

  !      else

  !         newEigenvalue=subEigenvalues%values(2)
  !         norm=1.0_8/sqrt( dot_product( subEigenvectors%values(:,2), subEigenvectors%values(:,2) ) )
  !         newEigenvector%values=norm*subEigenvectors%values(1,2)*lastEigenvector%values + &
  !                               norm*subEigenvectors%values(2,2)*direction%values

  !      end if

  !      lastEigenvector%values = newEigenvector%values/(Vector_norm(newEigenvector))

  !      difference = abs(lastEigenvalue - newEigenvalue)

  !      lastEigenvalue = newEigenvalue

  !      print *,"----------------------------------"
  !      print *,"RESULTS OF ITERATION NUMBER:",n
  !      print *,"this the subHamiltonian:"
  !      call Matrix_show(subHamiltonian)
  !      print*,"Eigenvector"
  !      call Vector_show(lastEigenvector)
  !      print *,"eigenvalue:",lastEigenvalue
  !      print*,"Residual"
  !      call Vector_show(residual)
  !      print*,"Direction"
  !      call Vector_show(direction)
  !      print*,"subHamiltonian matrix eigenvectors and eigenvalues"
  !      call Matrix_show(subEigenvectors)
  !      call Vector_show(subEigenvalues)
  !      print *,"----------------------------------"
       
  !   !!! loop finishs here
  !   end do

  !   ! Storing of corrections
    
  !   PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(speciesID,4*CONTROL_instance%IONIZE_MO)= &
  !        27.211396_8 * lastEigenvalue

  !   !!! DAVIDSON ALGORYTHM ENDS

  !   call Vector_destructor( lastEigenvector )
  !   call Vector_destructor( newEigenvector )
  !   call Vector_destructor( residual )
  !   call Vector_destructor( direction )
  !   call Vector_destructor( P )
  !   call Vector_destructor( auxVector )
  !   call Vector_destructor( subEigenvalues )
  !   call Matrix_destructor( subHamiltonian )
  !   call Matrix_destructor( subEigenvectors )
  !   call Vector_destructor( superEigenvalues )
  !   call Matrix_destructor( superEigenvectors )

  !   !!************************************************************************************************
  !   print *,"END OF CALCULATION OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR"
  !   print *,"**************************************************************"
  ! end subroutine PropagatorTheory_nonDiagonalSecondOrderCorrection


!   subroutine PropagatorTheory_nonDiagonalSecondOrderTDACorrection()
!     implicit NONE
    
!     integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha species
!     integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta species
!     integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma species
!     integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha species
!     integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta species
!     integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma species
!     integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha species
!     integer :: pb, qb, rb, sb ! Indices for general orbitals of beta species
!     integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma species
!     integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
!     integer :: i, j, k ! counters for species
!     integer :: m, n, o, p, q, r, s, t, u ! auxiliar counters
!     integer :: i2hp, i2ph
!     integer :: speciesID
!     integer :: otherSpeciesID, otherSpeciesID1, otherSpeciesID2
!     integer :: electronsID
!     integer :: occupationNumberOfSpecies, virtualNumberOfSpecies
!     integer :: occupationNumberOfOtherSpecies, virtualNumberOfOtherSpecies
!     integer :: occupationNumberOfOtherSpecies1, virtualNumberOfOtherSpecies1
!     integer :: occupationNumberOfOtherSpecies2, virtualNumberOfOtherSpecies2
!     integer :: numberOfContractionsOfSpecies
!     integer :: numberOfContractionsOfOtherSpecies
!     integer :: numberOfContractionsOfOtherSpecies1, numberOfContractionsOfOtherSpecies2
!     integer :: vectorSize, vectorSizeaHf
!     integer(4) :: errorNum
!     integer(8) :: auxIndex, HamiltonianSize, size2ph, size2hp
!     character(10) :: nameOfSpecies
!     character(10) :: nameOfOtherSpecies, nameOfOtherSpecies1, nameOfOtherSpecies2
!     type(Vector) :: occupationsOfOtherSpecies, occupationsOfOtherSpecies1, occupationsOfOtherSpecies2
!     type(Vector) :: eigenValuesOfSpecies
!     type(Vector) :: eigenValuesOfOtherSpecies, eigenValuesOfOtherSpecies1, eigenValuesOfOtherSpecies2
!     real(8) :: lambdaOfSpecies, lambdaOfOtherSpecies, lambdaOfOtherSpecies1, lambdaOfOtherSpecies2 
!     real(8) :: chargeOfSpecies, chargeOfOtherSpecies, chargeOfOtherSpecies1, chargeOfOtherSpecies2
!     real(8) :: auxValue, auxValue_A, auxValue_B, auxValue_1, auxValue_2
!     real(8) :: lastEigenvalue, newEigenvalue, norm
!     type(Vector) :: superEigenvalues, subEigenvalues, auxVector, W, projection, sizes
!     type(Vector) :: lastEigenvector, newEigenvector, residual, direction
!     type(Matrix),allocatable :: hph(:)      
!     type(Matrix),allocatable :: php(:)  
!     type(Matrix) :: aHa, subEigenvectors
!     type(Matrix) :: aHf, fHa, superHamiltonian, superEigenvectors
!     type(Matrix) :: subHamiltonian, zeroH, inverseZeroH
!     integer :: position
!     real(8) :: v, remainder, eW, eHo, difference, dot, im1, im2, im3
!     logical :: davidson

!     if ( .not.CONTROL_instance%OPTIMIZE ) then
!        print *,"================================================================================="
!        print *,"  If you are running propagator calculations the following papers must be cited: "
!        print *,"================================================================================="
!        print *," J. Chem. Phys. 137, 074105, (2012). J. Chem. Phys. 108, 1008, (1998).           "
!        print *," J. Chem. Phys. 99, 6716, (1993).                                                "
!        print *," Computational Chemistry: Reviews of Current Trends. Vol 2. Pag 1-61. 1997       "
!        print *,"================================================================================="
!     end if
    
!     print *,"******************************************************************"
!     print *,"BEGINNING OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR CALCULATIONS"

!     speciesID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
!     nameOfSpecies= trim(  MolecularSystem_getNameOfSpecie( speciesID ) )       
!     chargeOfSpecies = MolecularSystem_getCharge( speciesID )
!     eigenValuesOfSpecies = MolecularSystem_getEigenValues( speciesID )
!     occupationNumberOfSpecies = MolecularSystem_getOcupationNumber( speciesID )
!     numberOfContractionsOfSpecies = MolecularSystem_getTotalNumberOfContractions( speciesID )
!     virtualNumberOfSpecies = numberOfContractionsOfSpecies - occupationNumberOfSpecies 
!     lambdaOfSpecies = MolecularSystem_getLambda( speciesID )

!     !!! Construction of the 2hp and 2ph indices arranges (h-alpha h-p beta and p-alpha hp-beta)

!     if ( allocated(hph) ) deallocate(hph)
!     allocate(hph(PropagatorTheory_instance%numberOfSpecies))

!     if (allocated(php)) deallocate(php)
!     allocate(php(PropagatorTheory_instance%numberOfSpecies))

!     call Vector_constructor(sizes, int(PropagatorTheory_instance%numberOfSpecies))

!     vectorSizeaHf=0 !!! length of the (a|fH) block
!     size2hp=0 !!! length of a f(2hp) block
!     size2ph=0 !!! length of a f(2ph) block

!     !!! Constructing the 2hp and 2ph matrices

!     do j = 1 , PropagatorTheory_instance%numberOfSpecies

!           occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
!           numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
!           virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
          
!        if ( j == speciesID ) then !!! Here, m*(n-1)*n/2
          
!           size2ph = 0.5*(virtualNumberOfSpecies - 1)*occupationNumberOfSpecies*& 
!                virtualNumberOfSpecies 
!           size2hp = 0.5*(occupationNumberOfSpecies - 1)*occupationNumberOfSpecies*&
!                virtualNumberOfSpecies

!           if ( allocated(hph(j)%values ) ) deallocate(hph(j)%values)
!           call Matrix_constructor( hph(j), 3, size2hp)

!           if ( allocated(php(j)%values ) ) deallocate(php(j)%values)
!           call Matrix_constructor( php(j), 3, size2ph)

!           print *,"sizes:",size2hp,"",size2ph

!           sizes%values(j) = size2hp + size2ph
!           vectorSizeaHf = vectorSizeaHf + size2hp + size2ph

!        else
          
!           size2ph = virtualNumberOfOtherSpecies*occupationNumberOfOtherSpecies*& 
!                virtualNumberOfSpecies
!           size2hp = occupationNumberOfOtherSpecies* occupationNumberOfSpecies*&
!                virtualNumberOfOtherSpecies

!           print *,"sizes:",size2hp,"",size2ph

!           if ( allocated(hph(j)%values ) ) deallocate(hph(j)%values)
!           call Matrix_constructor( hph(j), 3, size2hp)

!           if ( allocated(php(j)%values ) ) deallocate(php(j)%values)
!           call Matrix_constructor( php(j), 3, size2ph)

!           sizes%values(j) = size2hp + size2ph
!           vectorSizeaHf = vectorSizeaHf + size2hp + size2ph
          
!        end if
       
!     end do

!     print *,"vectorSizeaHf:",vectorSizeaHf

!     !!! Storing 2hp (alpha alpha alpha) operators

!     m=0

!     do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!        do ia = 1 , occupationNumberOfSpecies
!           do ja = ia + 1 , occupationNumberOfSpecies
             
!              m = m + 1
!              hph(speciesID)%values(1,m)=aa
!              hph(speciesID)%values(2,m)=ia
!              hph(speciesID)%values(3,m)=ja

!           end do
!        end do
!     end do

!     !!! Storing 2ph operators (alpha alpha alpha)

!     m=0

!     do ia = 1 , occupationNumberOfSpecies
!        do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!           do ba = aa + 1 , numberOfContractionsOfSpecies

!              m = m + 1
!              php(speciesID)%values(1,m)=ia
!              php(speciesID)%values(2,m)=aa
!              php(speciesID)%values(3,m)=ba

!           end do
!        end do
!     end do

!     if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
             
!        do j = 1 , PropagatorTheory_instance%numberOfSpecies
!           if (j.ne.speciesID) then
             
!              occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
!              numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )

!              !!! Storing 2hp (beta alpha beta) operators

!              m = 0

!              do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
!                 do ia = 1 , occupationNumberOfSpecies 
!                   do ib = 1 , occupationNumberOfOtherSpecies

!                      m = m + 1
!                      hph(j)%values(1,m)=ab
!                      hph(j)%values(2,m)=ia
!                      hph(j)%values(3,m)=ib
                      
!                    end do
!                 end do
!              end do
             
!              !!! Storing 2ph (beta alpha beta) operators

!              m = 0             

!              do ib = 1 , occupationNumberOfOtherSpecies
!                 do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!                    do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
                      
!                      m = m + 1
!                      php(j)%values(1,m)=ib
!                      php(j)%values(2,m)=aa
!                      php(j)%values(3,m)=ab
                      
!                    end do
!                 end do
!              end do
             
!           end if
!        end do
       
!     end if

!     do j = 1, PropagatorTheory_instance%numberOfSpecies
    
!        print *,"hph of species:",j
!        call Matrix_show(hph(j))
!        print *,"php of species:",j
!        call Matrix_show(php(j))
       
!     end do

!    !!! CONSTRUCTION OF THE SUPEROPERATOR HAMILTONIAN BEGINS

!     HamiltonianSize = numberOfContractionsOfSpecies + vectorSizeaHf

!     print *,"HamiltonianSize:",HamiltonianSize

!     call Matrix_constructor(superHamiltonian, HamiltonianSize, HamiltonianSize)

!     !!!	BUILDING (a|Ha) BLOCK !!! I wish all the blocks were like this

!     do pa = 1, numberOfContractionsOfSpecies

!        superHamiltonian%values(pa,pa) = eigenValuesOfSpecies%values(pa)

!     end do

!     !!!	BUILDING (a|Hf) BLOCKS
    
!     do pa = 1 , numberOfContractionsOfSpecies !!! Run over the alpha holes and particles

!        idaHf = numberOfContractionsOfSpecies

!        !!! any (alpha) 2hp (alpha alpha alpha)

!        do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!           do ia = 1 , occupationNumberOfSpecies
!              do ja = ia + 1 , occupationNumberOfSpecies
                
!                 idaHf = idaHf + 1

!                 auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ia, aa, ja )
!                 auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ja, ia, aa )
!                 ! print *,"auxValue_A:",auxValue_A,"auxValue_B:",auxValue_B
!                 ! print*,"idaHf:",idaHf

!                 superHamiltonian%values(idaHf, pa) = (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B)
!                 superHamiltonian%values(pa, idaHf) = superHamiltonian%values(idaHf, pa)

!                 print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)
                
!              end do
!           end do
!        end do

!        !!! any (alpha) 2ph (alpha alpha alpha)
       
!        do ia = 1 , occupationNumberOfSpecies
!           do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!              do ba = aa + 1 , numberOfContractionsOfSpecies
                
!                 idaHf = idaHf + 1

!                 auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, aa, ia, ba )
!                 auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, pa, ba, aa, ia )
!                 ! print *,"auxValue_A:",auxValue_A,"auxValue_B:",auxValue_B
!                 ! print*,"idaHf:",idaHf

!                 superHamiltonian%values( idaHf, pa) = (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B)
!                 superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

!                 print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)

!              end do
!           end do
!        end do

!        if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
          
!           do j = 1 , PropagatorTheory_instance%numberOfSpecies
!              if (j.ne.speciesID) then
                
!                 nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( j ) )
!                 otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
!                 eigenValuesOfOtherSpecies = MolecularSystem_getEigenValues(j)
!                 occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
!                 numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
!                 lambdaOfOtherSpecies = MolecularSystem_getLambda( j )
!                 virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
!                 chargeOfOtherSpecies = MolecularSystem_getCharge( j )

!                !!! any (alpha) 2hp (beta alpha beta)                

!                 do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
!                    do ia = 1 , occupationNumberOfSpecies
!                       do ib = 1 , occupationNumberOfOtherSpecies
                         
!                          idaHf = idaHf + 1

!                          auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, pa, ia, ab, ib )
!                          ! print *,"auxValue:",auxValue
!                          ! print*,"idaHf:",idaHf

!                          superHamiltonian%values( idaHf, pa) = auxValue * chargeOfSpecies* chargeOfOtherSpecies !!! OJO lambdas
!                          superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

!                          print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)

!                       end do
!                    end do
!                 end do

!                !!! any (alpha) 2ph (beta alpha beta)                                

!                 do ib = 1 , occupationNumberOfOtherSpecies
!                    do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!                       do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
                         
!                          idaHf = idaHf + 1

!                          auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, pa, aa, ib, ab )     
!                          ! print *,"auxValue:",auxValue
!                          ! print*,"idaHf:",idaHf

!                          superHamiltonian%values( idaHf, pa) = auxValue * chargeOfSpecies* chargeOfOtherSpecies !!! OJO lambdas
!                          superHamiltonian%values( pa, idaHf) = superHamiltonian%values( idaHf, pa)

!                          print *,"idaHf:",idaHf,"pa:",pa,"value:",superHamiltonian%values( idaHf, pa)                         

!                       end do
!                    end do
!                 end do
                
!              end if
!           end do
          
!        end if
     
!     end do

!     !!!	BUILDING (f|Hf) BLOCK !!! THIS IS MADNEEEEESSSSSSSSSS

!     idfHf = numberOfContractionsOfSpecies !!! VERY IMPORTANT INDEX

!     !!! Intraspecies diagonal (f|Hf) blocks

!     !!! 2hp (alpha alpha alpha) 2hp (alpha alpha alpha)

!     o = Matrix_getNumberOfColumns(hph(speciesID)) + numberOfContractionsOfSpecies

!     do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!        do ia = 1 , occupationNumberOfSpecies
!           do ja = ia + 1 , occupationNumberOfSpecies
             
!              idfHf = idfHf + 1
!              print*,"idfHf:",idfHf
             
!              superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(ia) + eigenValuesOfSpecies%values(ja) &
!                                   - eigenValuesOfSpecies%values(aa) )

!              do m = idfHf, o
                
!                 n = m - numberOfContractionsOfSpecies

!                 if (aa==hph(speciesID)%values(1,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(hph(speciesID)%values(2,n)), ja, int(hph(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(hph(speciesID)%values(3,n)), int(hph(speciesID)%values(2,n)), ja )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ia==hph(speciesID)%values(2,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(1,n)), ja, int(hph(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(3,n)),int( hph(speciesID)%values(1,n)), ja )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ia==hph(speciesID)%values(3,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(1,n)), ja, int(hph(speciesID)%values(2,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(2,n)), int(hph(speciesID)%values(1,n)), ja )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ja==hph(speciesID)%values(2,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(1,n)), ia, int(hph(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(3,n)), int(hph(speciesID)%values(1,n)), ia )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ja==hph(speciesID)%values(3,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(1,n)), ia, int(hph(speciesID)%values(2,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ba, int(hph(speciesID)%values(2,n)), int(hph(speciesID)%values(1,n)), ia )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 superHamiltonian%values(m, idfHf) = superHamiltonian%values(idfHf, m)

!              end do

!           end do
!        end do
!     end do

!     !!! 2ph (alpha alpha alpha) 2ph (alpha alpha alpha)

!     p = o
!     o = sizes%values(speciesID) + numberOfContractionsOfSpecies

!     do ia = 1 , occupationNumberOfSpecies
!        do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!           do ba = aa + 1 , numberOfContractionsOfSpecies

!              idfHf = idfHf + 1
!              print*,"idfHf:",idfHf
             
!              superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(aa) + eigenValuesOfSpecies%values(ba) &
!                                   - eigenValuesOfSpecies%values(ia) )

!              do m = idfHf, o

!                 n = m - p

!                 if (ia==php(speciesID)%values(1,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(2,n)), ba, int(php(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(3,n)), int(php(speciesID)%values(2,n)), ba )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (aa==php(speciesID)%values(2,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(php(speciesID)%values(1,n)), ba, int(php(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(php(speciesID)%values(3,n)), int(php(speciesID)%values(1,n)), ba )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (aa==php(speciesID)%values(3,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(php(speciesID)%values(1,n)), ba, int(php(speciesID)%values(2,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, ia, int(php(speciesID)%values(2,n)), int(php(speciesID)%values(1,n)), ba )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ba==php(speciesID)%values(2,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(1,n)), ia, int(php(speciesID)%values(3,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(3,n)), int(php(speciesID)%values(1,n)), ia )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 if (ba==php(speciesID)%values(3,n)) then
                   
!                    auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(1,n)), ia, int(php(speciesID)%values(2,n)) )
!                    auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfSpecies, aa, int(php(speciesID)%values(2,n)), int(php(speciesID)%values(1,n)), ia )

!                    superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas

!                 end if

!                 superHamiltonian%values(m, idfHf) = superHamiltonian%values(idfHf, m)

!              end do

!           end do
!        end do
!     end do

!     !!! Interspecies (f|Hf) diagonal blocks

!     if ( PropagatorTheory_instance%numberOfSpecies.gt.1 ) then
             
!        do j = 1 , PropagatorTheory_instance%numberOfSpecies
!           if (j.ne.speciesID) then
             
!              nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( j ) )
!              otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
!              eigenValuesOfOtherSpecies = MolecularSystem_getEigenValues(j)
!              occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
!              numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
!              lambdaOfOtherSpecies = MolecularSystem_getLambda( j )
!              virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
!              chargeOfOtherSpecies = MolecularSystem_getCharge( j )

!              p = o
!              o = o + Matrix_getNumberOfColumns(hph(j))

!              print *,"p:",p,"o:",o

!             !!! 2hp (beta alpha beta) 2hp (beta alpha beta)             

!              do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
!                 do ia = 1 , occupationNumberOfSpecies 
!                   do ib = 1 , occupationNumberOfOtherSpecies

!                       idfHf = idfHf + 1
!                       i2hp = i2hp + 1
!                       print*,"idfHf:",idfHf
                      
!                       superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(ia) + eigenValuesOfOtherSpecies%values(ib) &
!                            - eigenValuesOfOtherSpecies%values(ab) )

!                       do m = idfHf, o
                         
!                          n = m - p
                         
!                          if (ab==hph(j)%values(1,n)) then

!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, ia, int(hph(j)%values(2,n)), ib, int(hph(j)%values(3,n)) )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                            
!                          end if
                         
!                          if (ia==hph(j)%values(2,n)) then
                            
!                             auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfOtherSpecies, ab, int(hph(j)%values(1,n)), ib, int(hph(j)%values(3,n)) )
!                             auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfOtherSpecies, ab, int(hph(j)%values(3,n)), int(hph(j)%values(1,n)), ib )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + (chargeOfOtherSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas
                            
!                          end if
                         
!                          if (ib==hph(j)%values(3,n)) then
                            
!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, ia, int(hph(j)%values(2,n)), ab, int(hph(j)%values(1,n)) )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                            
!                          end if
                         
!                          superHamiltonian%values(m, idfHf) = superHamiltonian%values(idfHf, m)
                         
!                       end do
                      
!                    end do
!                 end do
!              end do

!             !!! 2ph (beta alpha beta) 2ph (beta alpha beta)                          

!              p = o
!              o = o + Matrix_getNumberOfColumns(php(j))

!              print *,"p:",p,"o:",o

!              do ib = 1 , occupationNumberOfOtherSpecies
!                 do aa = occupationNumberOfSpecies+1 , numberOfContractionsOfSpecies
!                    do ab = occupationNumberOfOtherSpecies+1 , numberOfContractionsOfOtherSpecies
                      
!                       idfHf = idfHf + 1
!                       i2ph = i2ph + 1
!                       print*,"idfHf:",idfHf

!                       superHamiltonian%values(idfHf,idfHf) = ( eigenValuesOfSpecies%values(aa) + eigenValuesOfOtherSpecies%values(ab) &
!                            - eigenValuesOfOtherSpecies%values(ib) )

!                       do m = idfHf, o

!                          n = m - p

!                          if (ib==php(j)%values(1,n)) then

!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, aa, int(php(j)%values(2,n)), ab, int(php(j)%values(3,n)) )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) + ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                            
!                          end if
                         
!                          if (aa==php(j)%values(2,n)) then
                            
!                             auxValue_A = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfOtherSpecies, ib, int(php(j)%values(1,n)), ab, int(php(j)%values(3,n)) )
!                             auxValue_B = TransformIntegrals2_atomicToMolecularOfOneSpecieP2( nameOfOtherSpecies, ib, int(php(j)%values(3,n)), int(php(j)%values(1,n)), ab )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - (chargeOfOtherSpecies**2.0)*(auxValue_A - auxValue_B) !!! OJO lambdas
                            
!                          end if
                         
!                          if (ab==php(j)%values(3,n)) then
                            
!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, aa, int(php(j)%values(2,n)), ib, int(php(j)%values(1,n)) )
                            
!                             superHamiltonian%values(idfHf, m) = superHamiltonian%values(idfHf, m) - ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                            
!                          end if
                         
!                          superHamiltonian%values(m, idfHf) = superHamiltonian%values(idfHf, m)
                         
!                       end do
                      
!                    end do
!                 end do
!              end do
             
!           end if
!        end do

!        print *,"SUPERHAMILTONIAN"

!        call Matrix_show(superHamiltonian)

!        ! print *,"AQUI PAROOOOO"
!        ! stop 

!        print *,"counters after diagonal (f|Hf) blocks creation m:",m,"o:",o,"idfHf:",idfHf       
       
!        !!! (f|Hf) NON-DIAGONAL BLOCKS
       
!        !!! (alpha alpha alpha) vs (beta alpha beta) operators 
       
!        q = numberOfContractionsOfSpecies + Matrix_getNumberOfColumns(hph(speciesID)) !!! for alpha
!        r = numberOfContractionsOfSpecies + sizes%values(speciesID) !!! for alpha
       
!        s = numberOfContractionsOfSpecies + sizes%values(speciesID)  !!! for beta
       
!        do j = 1 , PropagatorTheory_instance%numberOfSpecies
!           if (j.ne.speciesID) then
             
!              nameOfOtherSpecies= trim(  MolecularSystem_getNameOfSpecie( j ) )
!              otherSpeciesID =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
!              eigenValuesOfOtherSpecies = MolecularSystem_getEigenValues(j)
!              occupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
!              numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
!              lambdaOfOtherSpecies = MolecularSystem_getLambda( j )
!              virtualNumberOfOtherSpecies = numberOfContractionsOfOtherSpecies - occupationNumberOfOtherSpecies
!              chargeOfOtherSpecies = MolecularSystem_getCharge( j )
             
!              t = s +  Matrix_getNumberOfColumns(hph(j)) !!! for beta
             
!              !!! 2hp (alpha alpha alpha) 2hp (beta alpha beta) 
             
!              do m = numberOfContractionsOfSpecies + 1, q !!! Indices for the position in the H matrix
!                 do n = s + 1, t                           !!! Indices for the position in the H matrix
                   
!                    o = m - numberOfContractionsOfSpecies !!! Indices for locating the corresponding 2ph (a-a-a) operator
!                    p = n - s                              !!! Indices for locating the corresponding 2ph (b-a-b) operator
                   
!                    if (hph(speciesID)%values(2,o) == hph(j)%values(2,p)) then
                      
!                       auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, int(hph(speciesID)%values(1,o)) , int(hph(speciesID)%values(3,o)),  int(hph(j)%values(1,p)), int(hph(j)%values(3,p)) )
                      
!                       superHamiltonian%values(m, n) = superHamiltonian%values(m, n) + ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                      
!                    end if
                   
!                    if (hph(speciesID)%values(3,o) == hph(j)%values(2,p)) then
                      
!                       auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, int(hph(speciesID)%values(1,o)) , int(hph(speciesID)%values(2,o)),  int(hph(j)%values(1,p)), int(hph(j)%values(3,p)) )
                      
!                       superHamiltonian%values(idfHf, m) = superHamiltonian%values(m, n) - ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                      
                      
!                    end if
                   
!                    superHamiltonian%values(n, m) = superHamiltonian%values(m, n)
                   
!                 end do
!              end do
             
!              s = t !!! for beta
!              t = s +  Matrix_getNumberOfColumns(php(j)) !!! for beta
             
!             !!! 2ph (alpha alpha alpha) 2ph (beta alpha beta) 
             
!              do m = q + 1, r     !!! Indices for the position in the H matrix
!                 do n = s + 1, t  !!! Indices for the position in the H matrix
                   
!                    o = m - q  !!! Indices for locating the corresponding 2ph (a-a-a) operator
!                    p = n - s  !!! Indices for locating the corresponding 2ph (b-a-b) operator
                   
!                    if (php(speciesID)%values(2,o) == php(j)%values(2,p)) then
                      
!                       auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, int(php(speciesID)%values(1,o)) , int(php(speciesID)%values(3,o)),  int(php(j)%values(1,p)), int(php(j)%values(3,p)) )
                      
!                       superHamiltonian%values(m, n) = superHamiltonian%values(m, n) - ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                      
!                    end if
                   
!                    if (php(speciesID)%values(3,o) == php(j)%values(2,p)) then
                      
!                       auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfSpecies , nameOfOtherSpecies, int(php(speciesID)%values(1,o)) , int(php(speciesID)%values(2,o)),  int(php(j)%values(1,p)), int(php(j)%values(3,p)) )
                      
!                       superHamiltonian%values(idfHf, m) = superHamiltonian%values(m, n) + ( chargeOfSpecies*chargeOfOtherSpecies* auxValue ) !!! OJO lambdas
                      
                      
!                    end if
                   
!                    superHamiltonian%values(n, m) = superHamiltonian%values(m, n)
                   
!                 end do
!              end do
             
!           end if
!        end do
       
!     end if
    
!     !!! (beta alpha beta) vs (gamma alpha gamma) operators 

!     if ( PropagatorTheory_instance%numberOfSpecies.gt.2 ) then
       
!        do i = 1 , PropagatorTheory_instance%numberOfSpecies - 1
!           if (i.ne.speciesID) then
             
!              nameOfOtherSpecies1= trim(  MolecularSystem_getNameOfSpecie( i ) )
!              otherSpeciesID1 =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
!              eigenValuesOfOtherSpecies1 = MolecularSystem_getEigenValues(j)
!              occupationNumberOfOtherSpecies1 = MolecularSystem_getOcupationNumber( i )
!              numberOfContractionsOfOtherSpecies1 = MolecularSystem_getTotalNumberOfContractions( i )
!              lambdaOfOtherSpecies1 = MolecularSystem_getLambda( i )
!              virtualNumberOfOtherSpecies1 = numberOfContractionsOfOtherSpecies1 - occupationNumberOfOtherSpecies1
!              chargeOfOtherSpecies1 = MolecularSystem_getCharge( i )
             
!              u = r
!              q = r + Matrix_getNumberOfColumns(hph(i)) !!! for beta
!              r = r + sizes%values(i) !!! for beta
             
!              s = u + sizes%values(i)  !!! for gamma
             
!              do j = i+1 , PropagatorTheory_instance%numberOfSpecies
!                 if (j.ne.speciesID) then
                   
!                    nameOfOtherSpecies2= trim(  MolecularSystem_getNameOfSpecie( j ) )
!                    otherSpeciesID2 =MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecies )
!                    eigenValuesOfOtherSpecies2 = MolecularSystem_getEigenValues(j)
!                    occupationNumberOfOtherSpecies2 = MolecularSystem_getOcupationNumber( j )
!                    numberOfContractionsOfOtherSpecies2 = MolecularSystem_getTotalNumberOfContractions( j )
!                    lambdaOfOtherSpecies2 = MolecularSystem_getLambda( j )
!                    virtualNumberOfOtherSpecies2 = numberOfContractionsOfOtherSpecies2 - occupationNumberOfOtherSpecies2
!                    chargeOfOtherSpecies2 = MolecularSystem_getCharge( j )

!                    !!! 2hp (beta alpha beta) 2hp (gamma alpha gamma) 
                   
!                    t = s +  Matrix_getNumberOfColumns(hph(j)) !!! for gamma
                   
!                    do m = u + 1, q !!! Indices for the position in the H matrix
!                       do n = s + 1, t                           !!! Indices for the position in the H matrix
                         
!                          o = m - u  !!! Indices for locating the corresponding 2hp (b-a-b) operator
!                          p = n - s  !!! Indices for locating the corresponding 2hp (c-a-c) operator
                         
!                          if (hph(i)%values(2,o) == hph(j)%values(2,p)) then
                            
!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfOtherSpecies1 , nameOfOtherSpecies2, int(hph(i)%values(1,o)) , int(hph(i)%values(3,o)),  int(hph(j)%values(1,p)), int(hph(j)%values(3,p)) )
                            
!                             superHamiltonian%values(m, n) = chargeOfOtherSpecies1*chargeOfOtherSpecies2*auxValue !!! OJO lambdas
                            
!                             superHamiltonian%values(n, m) = superHamiltonian%values(m, n)                         
                            
!                          end if
                         
!                       end do
!                    end do
 
!                    !!! 2ph (beta alpha beta) 2ph (gamma alpha gamma)                    
                   
!                    s = t !!! for gamma
!                    t = s +  Matrix_getNumberOfColumns(php(j)) !!! for gamma
                   
!                    do m = q + 1, r     !!! Indices for the position in the H matrix
!                       do n = s + 1, t  !!! Indices for the position in the H matrix
                         
!                          o = m - q  !!! Indices for locating the corresponding 2ph (b-a-b) operator
!                          p = n - s  !!! Indices for locating the corresponding 2ph (c-a-c) operator
                         
!                          if (php(i)%values(2,o) == php(j)%values(2,p)) then
                            
!                             auxValue = TransformIntegrals2_atomicToMolecularOfTwoSpeciesP2 ( nameOfOtherSpecies1 , nameOfOtherSpecies2, int(php(i)%values(1,o)) , int(php(i)%values(3,o)),  int(php(j)%values(1,p)), int(php(j)%values(3,p)) )
                            
!                             superHamiltonian%values(m, n) = chargeOfOtherSpecies1*chargeOfOtherSpecies2*auxValue !!! OJO lambdas
                            
!                             superHamiltonian%values(n, m) = superHamiltonian%values(m, n)                         
                            
!                          end if
                         
!                       end do
!                    end do
                   
!                 end if
!              end do
             
!           end if
!        end do
!     end if
    
! !!! CONSTRUCTION OF THE SUPEROPERATOR HAMILTONIAN ENDS 
    
!     print *,"idaHf:",idaHf,"idfHf:",idfHf,"HamiltonianSize:",HamiltonianSize

!     ! if ( CONTROL_instance%PT_JUST_ONE_ORBITAL) then
!     !    PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
!     !    PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
!     ! else
!     !    PropagatorTheory_instance%virtualBoundary= occupationNumberOfSpecies
!     !    PropagatorTheory_instance%occupationBoundary = 1
!     ! end if

!     !!! HERE THE CALCULATION OF OPTIMAL ENERGIES START

!     !!! DAVIDSON

!     call Vector_constructor( lastEigenvector, int(HamiltonianSize) )
!     call Vector_constructor( newEigenvector, int(HamiltonianSize) )

!     !!! Initialization of the algorithm

!     lastEigenvector%values = 0.0_8

!     lastEigenvector%values(CONTROL_instance%IONIZE_MO) = 1.0_8       

!     lastEigenvalue = eigenValuesOfSpecies%values(CONTROL_instance%IONIZE_MO)

!     print *,"Initial lastEigenvector values:"
!     call Vector_show(lastEigenvector)

!     !!! all the eigenvalues

!     call Vector_constructor(superEigenvalues, int(HamiltonianSize))
!     call Matrix_constructor(superEigenvectors, HamiltonianSize, HamiltonianSize)

!     call Matrix_eigen(superHamiltonian, superEigenvalues, superEigenvectors, SYMMETRIC)

!     print *,"koopmans:",eigenvaluesOfSpecies%values(CONTROL_instance%IONIZE_MO)

!     ! call Vector_show(superEigenvalues)
!     ! print *,"superEigenvalue:",superEigenvalues%values(CONTROL_instance%IONIZE_MO)

!     do k = 1, HamiltonianSize
!        dot=dot_product(superEigenvectors%values(:,k),lastEigenvector%values)
!        if (abs(dot)>0.8) then
!           print *,"Eigenvector:",k,"Eigenvalue:",superEigenvalues%values(k)
!           print *,superEigenvectors%values(:,k)
!        end if
!     end do

!     call Vector_Destructor(superEigenvalues)
!     call Matrix_Destructor(superEigenvectors)

!     davidson=.false.

!     if (davidson==.true.) then

!        call Vector_constructor( residual, int(HamiltonianSize) )
!        call Vector_constructor( projection, int(HamiltonianSize) )
!        call Vector_constructor( direction, int(HamiltonianSize) )
!        call Vector_constructor( auxVector, int(HamiltonianSize) )
!        call Vector_constructor( subEigenvalues, 2 )
!        call Matrix_constructor( subHamiltonian, 2, 2 )
!        call Matrix_constructor( subEigenvectors, 2, 2 )
       
!        !!! DAVIDSON ALGORYTHM STARTS
       
!        !!! loop starts here, on n, C must be normalized at this point
       
!        n = 0
!        difference = 1.0_8
       
!        do while ( difference > 0.00001_8 )
          
!           n = n+1
          
!           do pa = 1, HamiltonianSize
             
!              superHamiltonian%values(pa,pa) = superHamiltonian%values(pa,pa) - lastEigenvalue 
!              print *,"diagonal element of H:",pa,"value:",superHamiltonian%values(pa,pa)
             
!           end do
          
! !!! Substraction of En1
          
! !!! Calculating residual
          
!           residual%values = matmul(superHamiltonian%values,lastEigenvector%values)
          
! !!!! Calculating new direction
          
!           ! print *,"residual:"
          
!           ! call Vector_show(residual)
          
!           do pa = 1, numberOfContractionsOfSpecies
             
! !!! Be careful in the first iterations
             
!              if ( abs(superHamiltonian%values(pa,pa)) < 0.0000001_8  ) then
                
!                 projection%values(pa) = 0.0_8
!                 print *,"ENTRO:",pa
                
!              else
                
!                 projection%values(pa) = residual%values(pa)/superHamiltonian%values(pa,pa) 
                
!              end if
             
!           end do
          
!           do pa = numberOfContractionsOfSpecies+1, HamiltonianSize
             
!              projection%values(pa) = -residual%values(pa)/lastEigenvalue
             
!           end do
          
!           ! print *,"projection:"
          
!           ! call Vector_show(projection)
          
!           ! Orthogonalize again C
          
!           direction%values = projection%values - (dot_product(projection%values,lastEigenvector%values)/dot_product(lastEigenvector%values,lastEigenvector%values))*lastEigenvector%values
          
! !!! Normalization of the new direction
          
!           norm = Vector_norm( direction )
          
!           direction%values = direction%values/norm
          
!           ! print *,"direction:"
          
!           ! call Vector_show(direction)
          
!           print *,"dot product of C and d:",dot_product(direction%values,lastEigenvector%values)
          
! !!! Constructing superHamiltonian in the basis of C and d
          
!           do pa = 1, HamiltonianSize
             
!              superHamiltonian%values(pa,pa) = superHamiltonian%values(pa,pa) + lastEigenvalue 
             
!           end do
          
!           subHamiltonian%values(1,1) = lastEigenvalue
          
!           subHamiltonian%values(1,2) = dot_product( lastEigenvector%values, matmul(superHamiltonian%values,direction%values) )
          
!           subHamiltonian%values(2,1) = dot_product( direction%values, matmul(superHamiltonian%values, lastEigenvector%values) )
          
!           subHamiltonian%values(2,2) = dot_product( direction%values, matmul(superHamiltonian%values, direction%values) )
          
!           call Matrix_eigen( subHamiltonian, subEigenvalues, subEigenvectors, SYMMETRIC )
          
!           if (abs(lastEigenvalue-subEigenvalues%values(1))<abs(lastEigenvalue-subEigenvalues%values(2))) then
             
!              newEigenvalue=subEigenvalues%values(1)
!              norm=1.0_8/sqrt( dot_product( subEigenvectors%values(:,1), subEigenvectors%values(:,1) ) )
!              newEigenvector%values=norm*subEigenvectors%values(1,1)*lastEigenvector%values + &
!                   norm*subEigenvectors%values(2,1)*direction%values
             
!           else
             
!              newEigenvalue=subEigenvalues%values(2)
!              norm=1.0_8/sqrt( dot_product( subEigenvectors%values(:,2), subEigenvectors%values(:,2) ) )
!              newEigenvector%values=norm*subEigenvectors%values(1,2)*lastEigenvector%values + &
!                   norm*subEigenvectors%values(2,2)*direction%values
             
!           end if
          
!           lastEigenvector%values = newEigenvector%values/(Vector_norm(newEigenvector))
          
!           difference = abs(newEigenvalue - lastEigenvalue)
          
!           lastEigenvalue = newEigenvalue
          
!           print *,"----------------------------------"
!           print *,"RESULTS OF ITERATION NUMBER:",n
!           print *,"this the subHamiltonian:"
!           call Matrix_show(subHamiltonian)
!           print*,"Eigenvector"
!           call Vector_show(lastEigenvector)
!           print *,"eigenvalue:",lastEigenvalue
!           ! print*,"Residual"
!           ! call Vector_show(residual)
!           ! print*,"Direction"
!           ! call Vector_show(direction)
!           print*,"subHamiltonian matrix eigenvectors and eigenvalues"
!           call Matrix_show(subEigenvectors)
!           call Vector_show(subEigenvalues)
!           print *,"----------------------------------"

!           call Vector_destructor( residual )
!           call Vector_destructor( direction )
!           call Vector_destructor( projection )
!           call Vector_destructor( auxVector )
!           call Vector_destructor( subEigenvalues )
!           call Matrix_destructor( subHamiltonian )
!           call Matrix_destructor( subEigenvectors )
          
! !!! loop finishs here
!        end do
! !!! DAVIDSON ALGORYTHM ENDS    

!     else

! !!! INVERSE-ITERATION METHOD
       
!        call Vector_constructor( residual, int(HamiltonianSize) )
!        call Vector_constructor( projection, int(HamiltonianSize) )
!        call Vector_constructor( direction, int(HamiltonianSize) )
!        call Vector_constructor( auxVector, int(HamiltonianSize) )
!        call Matrix_constructor( zeroH, HamiltonianSize, HamiltonianSize, 0.0_8 )
       
!        n=0

!        do while(n<5)

!           n=n+1
          
! !!! Getting new eigenvalue
          
!           newEigenvalue = dot_product(lastEigenvector%values,matmul(superHamiltonian%values,lastEigenvector%values))
          
! !!! Building V
          
!           do ia = 1, numberOfContractionsOfSpecies
             
!              superHamiltonian%values(ia,ia) = superHamiltonian%values(ia,ia) - eigenValuesOfSpecies%values(ia)
             
!           end do
          
! !!! Building inverse matrix (EI-H_0)^-1
          
!           do ia = 1, numberOfContractionsOfSpecies
             
!              if (abs(newEigenvalue - eigenValuesOfSpecies%values(ia))<0.000001_8) then
                
!                 zeroH%values(ia,ia)=0.0_8
!                 position=ia
                
!              else
                
!                 zeroH%values(ia,ia) = 1.0_8/(newEigenvalue - eigenValuesOfSpecies%values(ia))
                
!              end if
             
!           end do
          
!           do ia = numberOfContractionsOfSpecies + 1, HamiltonianSize
             
!              zeroH%values(ia,ia) = 1.0_8/newEigenvalue
             
!           end do
          
! !!! Building W=projection
          
!           projection%values = matmul(superHamiltonian%values,lastEigenvector%values)
          
!           im1=dot_product(lastEigenvector%values,matmul(zeroH%values,lastEigenvector%values))
!           im2=dot_product(lastEigenvector%values,matmul(zeroH%values,projection%values))
!           im3=(im2-1.0_8)/im1
          
!           print *,"im1:",im1,"im2:",im2,"im3:",im3
          
!           newEigenvector%values=matmul(zeroH%values,(projection%values-im3*lastEigenvector%values))
          
!           norm=Vector_norm(newEigenvector)
          
!           lastEigenvector%values=newEigenvector%values/norm
          
! !!! Restoring superHamiltonian
          
!           do ia = 1, numberOfContractionsOfSpecies
             
!              superHamiltonian%values(ia,ia) = superHamiltonian%values(ia,ia) + eigenValuesOfSpecies%values(ia)
             
!           end do
          
!           print *,"----------------------------------"
!           print *,"RESULTS OF ITERATION NUMBER:",n
!           print*,"Eigenvector"
!           call Vector_show(lastEigenvector)
!           print *,"eigenvalue:",lastEigenvalue
!           print *,"----------------------------------"

!        end do

!        call Vector_destructor( residual )
!        call Vector_destructor( direction )
!        call Vector_destructor( projection )
!        call Vector_destructor( auxVector )
!        call Matrix_destructor( zeroH )       

!     end if

!     call Vector_destructor( lastEigenvector )
!     call Vector_destructor( newEigenvector )

    
!     ! Storing of corrections
    
!     PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(speciesID,4*CONTROL_instance%IONIZE_MO)= &
!          27.211396_8 * lastEigenvalue
    
!     ! PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-1))= &
!     !      PropagatorTheory_getPolarStrength( auxNumeratorsVector, &
!     !      auxDenominatorsVector, id, optimizedOmega)
    
!     ! PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-2))= 27.211396_8 * selfEnergyIntraSpecie
!     ! PropagatorTheory_instance%energyCorrectionsOfSecondOrder%values(i,(4*m-3))= 27.211396_8 * selfEnergyInterSpecie

!     !!************************************************************************************************
!     print *,"END OF CALCULATION OF SECOND ORDER ELECTRON-NUCLEAR PROPAGATOR"
!     print *,"**************************************************************"
!   end subroutine PropagatorTheory_nonDiagonalSecondOrderTDACorrection

  !**
  ! @brief Evaluate partial third order, EP3 and OVGF poles
  ! ... well EP3 and OVGF soon
  !**

!  subroutine PropagatorTheory_thirdOrderCorrection()
!    implicit NONE
!    
!    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
!    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
!    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
!    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
!    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
!    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
!    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
!    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
!    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
!    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
!    integer :: i, j, k ! counters for species
!    integer :: m, n, o, p, q, ni, nc, limit, id1, id2 ! auxiliar counters
!    integer :: speciesAID, speciesBID, speciesCID
!    integer :: species1ID, species2ID
!    integer :: electronsID
!    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
!    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
!    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
!    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
!    integer :: activeOrbitalsOfSpeciesC
!    integer(8) :: vectorSize1, vectorSize2 !!! Sizes for diagrams
!    integer(8) :: auxIndex
!    integer(4) :: errorNum
!    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
!    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
!    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
!    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
!    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!!    type(TransformIntegrals) :: repulsionTransformer
!    type(Matrix),allocatable :: auxMatrix2(:), selfEnergy2hp(:), selfEnergy2ph(:)
!    type(Matrix) :: diagram_A, diagram_B, auxMatrix, auxMatrix3
!    type(Matrix) :: partialMO1, partialMO2    
!    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
!    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
!    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
!    real(8) :: factorA_numerator, factorB_numerator, factorA_denominator, factorB_denominator
!    real(8) :: valueOfW, valueOfU, valueOfdU
!    real(8) :: lastOmega, newOmega, residual, selfEnergy, selfEnergyDerivative, koopmans 
!    real(8) :: a1, a2, b, c, d, poleStrenght, secondOrderStrength, thirdOrderStrength
!    real(8) :: numerator1, denominator1, numerator2, denominator2
!    logical :: paso1, paso2
!    ! *******************************************************************************************
!    ! Determinate the numerators and denominators of the second Oder propapator 
!    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
!    
!    print *,"*******************************************************************"
!    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"
!
!    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified
!
!    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
!    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies))
!
!    !!! Defining for which species the correction will be applied
!    
!    if (CONTROL_instance%IONIZE_SPECIE /= "NONE") then
!       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
!       species2ID= species1ID
!       m=1
!    else
!       species1ID=1
!       species2ID=PropagatorTheory_instance%numberOfSpecies
!       m = species2ID
!    end if
!
!    if (allocated(PropagatorTheory_instance%thirdOrderCorrections)) deallocate(PropagatorTheory_instance%thirdOrderCorrections)
!    allocate(PropagatorTheory_instance%thirdOrderCorrections(m))
!
!    !!! Start loop for species
!
!    q = 0
!
!    do i = species1ID , species2ID
!       
!       q = q + 1
!       
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
!       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!!       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
!       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )
!       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i )
!       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
!       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA
!
!       ! paso
!
!       paso1=(nameOfSpeciesA=="e-ALPHA".or.nameOfSpeciesA=="e-BETA")
!
!       !!! Defining the number of orbitals !!! Insert a parameter for the else option
!
!       if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
!          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
!          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
!          n = 1
!       else if (CONTROL_instance%IONIZE_SPECIE /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
!          n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
!       else
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
!          n = 2
!       end if
!
!       call Matrix_constructor(PropagatorTheory_instance%thirdOrderCorrections(q), int(n,8), 8, 0.0_8)
!       ! Storing transformed integrals !!!! We need a more efficient algorithm for this
!       
!!       call TransformIntegrals_constructor( repulsionTransformer )
!       
!       do p = 1 , PropagatorTheory_instance%numberOfSpecies
!          
!          if (p==i) then
!             
!!             call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p), p, trim(nameOfSpeciesA) )
!
!             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( p ) &
!                  * MolecularSystem_getCharge( p )
!             
!          else
!             
!             nameOfSpeciesB= trim(  MolecularSystem_getNameOfSpecie( p ) )
!             
!!             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                  MolecularSystem_getEigenVectors(i), MolecularSystem_getEigenVectors(p), &
!!                  auxMatrix2(p), i, nameOfSpeciesA, p, nameOfSpeciesB )
!
!             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( i ) &
!                  * MolecularSystem_getCharge( p )
!             
!          end if
!          
!       end do
!       
!       !**************************************************************************
!       !	Storing of denominators and numerators in the corresponding vectors
!       !****
!
!       m =0
!
!       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	
!
!          m=m+1          
!
!          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
!          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
!          
!          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
!          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))
!          
!          do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!             
!             if (j==i) then ! Intraspecies factors
!
!                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
!                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
!
!                call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!
!                id1 = 0
!                id2 = 0
!
!                ! factor 2ph
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
!                         auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                         auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
!                         auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                         
!                         id1 = id1 + 1
!                         
!                         selfEnergy2ph(j)%values(1,id1) = auxValue_A - auxValue_B
!                         
!                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
!                              - eigenValuesOfSpeciesA%values(ba)
!                         
!                         valueOfW = 0.0_8
!
!                         if (.not.paso1) then
!                            
!                            do ja = 1, occupationNumberOfSpeciesA
!                               do ka = 1, occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, ka, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ka, ia, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                  
!                               end do
!                            end do
!                            
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ja = 1, occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, ba, ia, ca, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, aa, ia, ca, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca) )
!                                  
!                               end do
!                            end do
!                            
!                            do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                               
!                               if (k .ne. i)  then
!                                  
!                                  nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                  chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                  eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                  occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                  activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                  lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                  virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                  
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                        
!                                        valueOfW = valueOfW + 1.0_8*(auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                        
!                                        valueOfW = valueOfW - 1.0_8*(auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                        
!                                     end do
!                                  end do
!                                  
!                               end if
!                               
!                            end do
!
!                         end if
!
!                         selfEnergy2ph(j)%values(3,id1) = valueOfW
!
!                      end do
!                   end do
!                end do
!                
!                if (occupationNumberOfSpeciesA > 1) then
!
!                   ! factor 2hp
!                   
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            
!                            id2 = id2 + 1
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ia, aa, ja, activeOrbitalsOfSpeciesA )
!                            auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ja, aa, ia, activeOrbitalsOfSpeciesA )
!                            auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                            
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A - auxValue_B
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesA%values(ja) 
!
!                            valueOfW = 0.0_8
!
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, ca, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ca, aa, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ia, ca, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!
!                               end do
!                            end do
!
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ka = 1, occupationNumberOfSpeciesA
!
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ja, aa, ka, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!
!                                  valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                  auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                  auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ia, aa, ka, activeOrbitalsOfSpeciesA )
!                                  auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ka) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!
!                               end do
!                            end do
!
!                            do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                               
!                               if (k .ne. i)  then
!
!                                  nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                  chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                  eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                  occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                  activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                  lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                  virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!
!                                        valueOfW = valueOfW + 1.0_8*(auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                                                                
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!
!                                        valueOfW = valueOfW - 1.0_8*(auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                                                                
!                                     end do
!                                  end do
!
!                               end if
!                               
!                            end do
!
!                            selfEnergy2hp(j)%values(3,id2) = valueOfW
!                            
!                         end do
!                      end do
!                   end do
!
!                end if
!                
!             else ! interspecies
!
!                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                ! paso
!
!                paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                     (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!            
!                vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
!                vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB
!
!!                call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                     MolecularSystem_getEigenvectors(j), auxMatrix, j, trim(nameOfSpeciesB) )
!                
!                auxMatrix%values = auxMatrix%values * ( chargeOfSpeciesB**2.0_8 )
!                         
!                call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!
!                id1 = 0
!                id2 = 0
!
!                ! diagram A
!                
!                do ib = 1 , occupationNumberOfSpeciesB
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                         auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                         
!                         id1 = id1 + 1
!                         
!                         selfEnergy2ph(j)%values(1,id1) = auxValue_A
!                         
!                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
!                              - eigenValuesOfSpeciesB%values(ab)
!                         
!                         valueOfW = 0.0_8
!
!                         if (.not.paso2) then
!                            
!                            do ia = 1, occupationNumberOfSpeciesA
!                               do jb = 1, occupationNumberOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, jb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib,  bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!                            
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_A = auxMatrix2(i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ba, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_B = auxMatrix2(i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, bb, ib, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, jb, bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!                            
!                         end if
!                         
!                         selfEnergy2ph(j)%values(3,id1) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!                ! diagram B
!                
!                do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         
!                         id2 = id2 + 1
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                         auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                         
!                         selfEnergy2hp(j)%values(1,id2) = auxValue_A
!                         
!                         selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
!                              - eigenValuesOfSpeciesB%values(ib)
!                         
!                         valueOfW = 0.0_8
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do jb = 1 , occupationNumberOfSpeciesB
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ia, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B = auxMatrix2(i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ab, ib, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ab, jb, bb, ib, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, jb, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!
!                         selfEnergy2hp(j)%values(3,id2) = valueOfW
!
!                      end do
!                   end do
!                end do
!                                
!             end if
!
!          end do
!
!          ! Initial guess
!          secondOrderStrength = 0.0_8
!          koopmans = eigenValuesOfSpeciesA%values(pa)
!          newOmega = koopmans
!          lastOmega = 0.0_8
!          
!          ni = 0
!          limit = 50
!          residual = 1.0_8
!
!          ! Calculation of second order pole
!
!          do while ((residual>0.0001_8).or.(limit.lt.ni))
!             
!             ni = ni + 1
!             
!             lastOmega = newOmega
!             selfEnergy = lastOmega - koopmans
!             selfEnergyDerivative = 1.0_8
!             
!             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                
!                if (j==i) then ! Intraspecies term
!                   
!                   id1=0
!                   id2=0
!                   
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            
!                            id1 = id1 + 1
!                            
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!
!                            if (ni==1) secondOrderStrength = secondOrderStrength - 0.5_8*a1*a1/b
!                            
!                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!
!                               if (ni==1) secondOrderStrength = secondOrderStrength - 0.5_8*a1*a1/b
!
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                                                           
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! Interspecies term
!
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                   id1 = 0
!                   id2 = 0
!
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            id1 = id1 + 1
!
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*a1/b
!
!                            if (ni==1) secondOrderStrength = secondOrderStrength - 0.5_8*a1*a1/b
!
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!
!                            a1 = selfEnergy2hp(j)%values(1,id2)
!                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                            
!                            selfEnergy = selfEnergy -a1*a1/b
!
!                            if (ni==1) secondOrderStrength = secondOrderStrength - 0.5_8*a1*a1/b
!
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                end if
!                
!             end do
!             
!             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!             
!             residual = abs(newOmega-lastOmega)
!             
!             print *,"iteration",ni,"newOmega",newOmega,"residual",residual
!
!          end do ! while
!
!          poleStrenght = 1.0_8/selfEnergyDerivative
!
!          ! Storing corrections
!
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)=real(pa,8)
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3)=27.211396_8 * newOmega
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)=poleStrenght
!
!          print *,"----------------------------------------------------------------"
!          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)),&
!               " of species: ",nameOfSpeciesA
!          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T5,A29,F8.4,A7,I2,A12)") "Optimized second order pole: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               " after ",ni," iterations."
!          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!
!          ! Initial guess
!          koopmans = eigenValuesOfSpeciesA%values(pa)
!          newOmega = koopmans
!          lastOmega = 0.0_8
!          
!          ni = 0
!          limit = 50
!          residual = 1.0_8
!
!          ! Calculation of P3 pole
!
!          do while ((residual>0.0001_8).or.(limit.lt.ni))
!             
!             ni = ni + 1
!             
!             secondOrderStrength = 0.0_8
!             thirdOrderStrength = 0.0_8
!             lastOmega = newOmega
!             selfEnergy = lastOmega - koopmans
!             selfEnergyDerivative = 1.0_8
!             
!             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                
!                if (j==i) then ! Intraspecies term
!                   
!                   id1=0
!                   id2=0
!                   
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            
!                            id1 = id1 + 1
!                            
!                            valueOfU = 0.0_8
!                            valueOfdU = 0.0_8
!
!                            if (.not.paso1) then
!                               
!                               do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do da = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, ia, da, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, da, ia, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ca, aa, da, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ca, ba, da, aa, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(ia) &
!                                          - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(da)
!                                     
!                                     valueOfU = valueOfU + 0.5_8*a2/c
!                                     valueOfdU = valueOfdU - 0.5_8*a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueOfdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                  
!                                  if (k .ne. i)  then
!                                     
!                                     nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                     chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                     eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                     occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                     activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                     lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                     virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesB
!                                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                           
!                                           valueOfU = valueOfU + 2.0_8*a2/c
!                                           valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ia, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                           
!                                           valueOfU = valueOfU - 2.0_8*a2/c
!                                           valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                  end if
!                                  
!                               end do
!
!                            end if
!                                                        
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            a2 = selfEnergy2ph(j)%values(3,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - 0.5_8*a1*(a1+a2+valueOfU)/b
!
!                            secondOrderStrength = secondOrderStrength - 0.5_8*a1*(a1)/b
!                            thirdOrderStrength = thirdOrderStrength - 0.5_8*a1*(a2+valueOfU)/b
!
!                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( a1*(a1+a2+valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               do ka = 1 , occupationNumberOfSpeciesA
!                                  do la = 1 , occupationNumberOfSpeciesA
!
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ka, aa, la, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, la, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ka, ia, la, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ka, ja, la, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(aa) &
!                                          - eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(la)
!                                     
!                                     valueOfU = valueOfU - 0.5_8*a2/c
!                                     valueOfdU = valueOfdU + 0.5_8*a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!
!                               do ka = 1 , occupationNumberOfSpeciesA
!                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                          - eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ka)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueOfdU = valueOfdU + a2/(c**2.0_8)
!
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!
!                                     a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                     c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                          - eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ka)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!
!                                  end do
!                               end do
!
!                               do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                  
!                                  if (k .ne. i)  then
!                                     
!                                     nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                     chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                     eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                     occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                     activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                     lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                     virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesB
!                                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                           
!                                           valueOfU = valueOfU + 2.0_8*a2/c
!                                           valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
! 
!                                          auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ja, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ia)
!                                           
!                                           valueOfU = valueOfU - 2.0_8*a2/c
!                                           valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                  end if
!                                  
!                               end do
!
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               a2 = selfEnergy2hp(j)%values(3,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!
!                               selfEnergy = selfEnergy - 0.5_8*a1*(a1+a2+valueOfU)/b
!
!                               secondOrderStrength = secondOrderStrength - 0.5_8*a1*(a1)/b
!                               thirdOrderStrength = thirdOrderStrength - 0.5_8*a1*(a2+valueOfU)/b
!
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( a1*(a1+a2+valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!                            
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! Interspecies term
!
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                   paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                        (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                   
!!                   call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                        MolecularSystem_getEigenvectors(j), auxMatrix, j, trim(nameOfSpeciesB) )
!                   
!                   auxMatrix%values = auxMatrix%values * ( chargeOfSpeciesB**2.0_8 )
!                   
!                   id1 = 0
!                   id2 = 0
!
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            id1 = id1 + 1
!
!                            valueOfU = 0.0_8
!                            valueOfdU = 0.0_8
!
!                            if (.not.paso2) then
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, aa, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!
!                                     a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                          - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(aa)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, aa, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, aa, ab, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                          - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(ba)
!                                     
!                                     ! valueOfU = valueOfU + (a2/c)*(1.0_8/(1.0_8-auxValue_C))
!                                     ! valueofdU = valueOfdU - (a2/(c**2.0_8))*(1.0_8/(1.0_8-auxValue_C))
!                                     valueOfU = valueOfU + (a2/c)
!                                     valueofdU = valueOfdU - (a2/(c**2.0_8))
!                                     
!                                  end do
!                               end do
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do jb = 1 , occupationNumberOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                          - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!
!                            end if
!                            
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            a2 = selfEnergy2ph(j)%values(3,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*(a1+a2+valueOfU)/b 
!
!                            secondOrderStrength = secondOrderStrength - 0.5_8*a1*(a1)/b
!                            thirdOrderStrength = thirdOrderStrength - 0.5_8*a1*(a2+valueOfU)/b
!
!                            selfEnergyDerivative = selfEnergyDerivative + ( a1*(a1+a2+valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!
!                            valueOfU = 0.0_8
!                            valueOfdU = 0.0_8
!
!                            do jb = 1 , occupationNumberOfSpeciesB
!
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                                                    
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!
!                                  auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_C= auxMatrix%values(auxIndex, 1)
!                               
!                                  a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                  c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                       - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ia)
!                                  
!                                  valueOfU = valueOfU + a2/c
!                                  valueofdU = valueOfdU - a2/(c**2.0_8)
!                                  
!                               end do
!                            end do
!                            
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  a2 = auxValue_A*auxValue_B
!                                  c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                       - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ja)
!                                  
!                                  valueOfU = valueOfU - a2/c
!                                  valueofdU = valueOfdU + a2/(c**2.0_8)
!                                  
!                               end do
!                            end do
!                            
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ja, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ja, bb, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  a2 = auxValue_A*auxValue_B
!                                  c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                       - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                  
!                                  valueOfU = valueOfU + a2/c
!                                  valueofdU = valueOfdU - a2/(c**2.0_8)
!                                  
!                               end do
!                            end do
!
!                            a1 = selfEnergy2hp(j)%values(1,id2)
!                            a2 = selfEnergy2hp(j)%values(3,id2)
!                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*(a1+a2+valueOfU)/b 
!
!                            secondOrderStrength = secondOrderStrength - 0.5_8*a1*(a1)/b
!                            thirdOrderStrength = thirdOrderStrength - 0.5_8*a1*(a2+valueOfU)/b
!
!                            selfEnergyDerivative = selfEnergyDerivative + ( a1*(a1+a2+valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!                            
!                         end do
!                      end do
!                   end do
!
!                   do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!
!                      id1=0
!                      id2=0
!
!                      if (k.ne.i .and. k.ne.j)  then
!                         
!                         print *,"entro al manolito",k
!                         
!                         nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                         chargeOfSpeciesC = MolecularSystem_getCharge( k )
!!                         eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( k )
!                         occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( k )
!                         activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( k )
!                         lambdaOfSpeciesC = MolecularSystem_getLambda( k )
!                         virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                         
!!                         call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                              MolecularSystem_getEigenVectors(j), MolecularSystem_getEigenVectors(k), &
!!                              auxMatrix3, j, nameOfSpeciesB, k, nameOfSpeciesC )
!                         
!                         auxMatrix3%values = auxMatrix3%values * (chargeOfSpeciesB*chargeOfSpeciesC)
!                         
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                  id1 = id1 + 1
!
!                                  valueOfU=0.0_8
!                                  valueOfdU=0.0_8
!
!                                  do ic = 1 , occupationNumberOfSpeciesC
!                                     do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                               
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                        auxValue_B = auxMatrix3%values(auxIndex, 1)
!                                        
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesC%values(ic) &
!                                             - eigenValuesOfSpeciesC%values(ac) - eigenValuesOfSpeciesA%values(aa)
!                                        
!                                       
!                                        valueOfU = valueOfU + a2/c
!                                        valueofdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!
!                                  a1 = selfEnergy2ph(j)%values(1,id1)
!                                  a2 = selfEnergy2ph(j)%values(3,id1)
!                                  b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                                  
!                                  selfEnergy = selfEnergy - a1*(valueOfU)/b 
!
!                                  thirdOrderStrength = thirdOrderStrength - a1*(valueOfU)/b 
!
!                                  selfEnergyDerivative = selfEnergyDerivative + ( a1*(valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!                                  
!                               end do
!                            end do
!                         end do
!
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            do ia = 1 , occupationNumberOfSpeciesA
!                               do ib = 1 , occupationNumberOfSpeciesB
!                                  
!                                  id2 = id2 + 1
!                         
!                                  valueOfU=0.0_8
!                                  valueOfdU=0.0_8
!
!                                  do ic = 1 , occupationNumberOfSpeciesC
!                                     do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                        auxValue_B = auxMatrix3%values(auxIndex, 1)
!                                        
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesC%values(ac) &
!                                             - eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesA%values(ia)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!
!                                  a1 = selfEnergy2hp(j)%values(1,id2)
!                                  a2 = selfEnergy2hp(j)%values(3,id2)
!                                  b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                  
!                                  selfEnergy = selfEnergy - a1*(valueOfU)/b 
!
!                                  thirdOrderStrength = thirdOrderStrength - a1*(valueOfU)/b 
!
!                                  selfEnergyDerivative = selfEnergyDerivative + ( a1*(valueOfU)/(b**2.0_8) - a1*valueOfdU/b )
!
!                               end do
!                            end do
!                         end do
!                         
!                      end if
!                      
!                   end do
!                   
!                end if
!                
!             end do
!             
!             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!             
!             residual = abs(newOmega-lastOmega)
!             
!             print *,"ni",ni,"newOmega",newOmega,"residual",residual
!
!          end do ! while
!
!          poleStrenght = 1.0_8/(selfEnergyDerivative)
!
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,5)=27.211396_8 * newOmega
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,6)=poleStrenght
!
!          write (*,"(T5,A26,F8.4,A7,I2,A12)") "Optimized P3 order pole: ",newOmega*27.211396_8," after ",ni," iterations."
!          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!          print *,"----------------------------------------------------------------"
!
!       end do
!
!       ! call Matrix_destructor(auxMatrix2(:))          
!!       call TransformIntegrals_destructor( repulsionTransformer )
!       
!    end do
!       
!    !!
!    !!************************************************************************************************
!    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
!    print *,"***************************************************************"
!  end subroutine PropagatorTheory_thirdOrderCorrection
!
!  !**
!  ! @brief Evaluate partial third order, EP3 and OVGF poles
!  ! version 2
!  !**
!
!  subroutine PropagatorTheory_thirdOrderCorrection2()
!    implicit NONE
!    
!    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
!    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
!    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
!    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
!    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
!    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
!    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
!    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
!    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
!    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
!    integer :: i, j, k ! counters for species
!    integer :: m, n, o, p, q, ni, nc, limit, id1, id2 ! auxiliar counters
!    integer :: speciesAID, speciesBID, speciesCID
!    integer :: species1ID, species2ID
!    integer :: electronsID
!    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
!    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
!    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
!    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
!    integer :: activeOrbitalsOfSpeciesC
!    integer(8) :: vectorSize1, vectorSize2 !!! Sizes for diagrams
!    integer(8) :: auxIndex
!    integer(4) :: errorNum
!    character(10) :: thirdOrderMethods(5)
!    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
!    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
!    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
!    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
!    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!!    type(TransformIntegrals) :: repulsionTransformer
!    type(Matrix),allocatable :: auxMatrix2(:), selfEnergy2hp(:), selfEnergy2ph(:)
!    type(Matrix),allocatable :: secondOrderDensities(:)
!    type(Matrix) :: diagram_A, diagram_B, auxMatrix, auxMatrix3
!    type(Matrix) :: partialMO1, partialMO2    
!    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
!    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
!    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
!    real(8) :: valueOfW, valueOfU, valueOfdU, sub2, subW, subU, subd2, subdW, subdU
!    real(8) :: lastOmega, newOmega, residual, threshold, selfEnergy, selfEnergyDerivative, koopmans 
!    real(8) :: a1, a2, b, c, d, poleStrenght, partialValue, initialValue
!    real(8) :: factors(3,5), fW, fI, constantSelfEnergy,thirdOrderResults(2,5)
!    real(8) :: s2hp, s2ph, W2hp, W2ph, U2hp, U2ph
!    logical :: paso1, paso2
!    ! *******************************************************************************************
!    ! Determinate the numerators and denominators of the second Oder propapator 
!    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
!    
!    print *,"*******************************************************************"
!    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"
!
!    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified
!
!    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
!    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(secondOrderDensities)) deallocate(secondOrderDensities)
!    allocate(secondOrderDensities(PropagatorTheory_instance%numberOfSpecies))
!
!    !!! Defining for which species the correction will be applied
!    
!    if (CONTROL_instance%IONIZE_SPECIE /= "NONE") then
!       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
!       species2ID= species1ID
!       m=1
!    else
!       species1ID=1
!       species2ID=PropagatorTheory_instance%numberOfSpecies
!       m = species2ID
!    end if
!
!    if (allocated(PropagatorTheory_instance%thirdOrderCorrections)) deallocate(PropagatorTheory_instance%thirdOrderCorrections)
!    allocate(PropagatorTheory_instance%thirdOrderCorrections(m))
!
!    !!! Start loop for species
!    
!    q = 0
!    
!    do i = species1ID , species2ID
!       
!       q = q + 1
!       
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
!       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!!       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
!       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )
!       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i )
!       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
!       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA
!       
!       ! paso
!       
!       paso1=(nameOfSpeciesA=="e-ALPHA".or.nameOfSpeciesA=="e-BETA")
!       
!       ! Defining the number of orbitals !!! Insert a parameter for the else option
!       
!       if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
!          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
!          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
!          n = 1
!       else if (CONTROL_instance%IONIZE_SPECIE /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
!          n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
!       else
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
!          n = 2
!       end if
!
!       call Matrix_constructor(PropagatorTheory_instance%thirdOrderCorrections(q), int(n,8), 8, 0.0_8)
!       
!       ! Storing transformed integrals !!!! We need a more efficient algorithm for this
!       
!!       call TransformIntegrals_constructor( repulsionTransformer )
!       
!       do p = 1 , PropagatorTheory_instance%numberOfSpecies
!          
!          if (p==i) then
!             
!!             call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p), p, trim(nameOfSpeciesA) )
!             
!             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( p ) &
!                  * MolecularSystem_getCharge( p )
!             
!          else
!             
!             nameOfSpeciesB= trim(  MolecularSystem_getNameOfSpecie( p ) )
!             
!!             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                  MolecularSystem_getEigenVectors(i), MolecularSystem_getEigenVectors(p), &
!!                  auxMatrix2(p), i, nameOfSpeciesA, p, nameOfSpeciesB )
!             
!             auxMatrix2(p)%values = auxMatrix2(p)%values * MolecularSystem_getCharge( i ) &
!                  * MolecularSystem_getCharge( p )
!             
!          end if
!          
!       end do
!
!       if (q==1) then
!          
!          ! second order densities 1 (for alpha)
!          print *,"entro a densities 1"
!          
!          call Matrix_constructor(secondOrderDensities(i), int(activeOrbitalsOfSpeciesA,8),&
!               int(activeOrbitalsOfSpeciesA,8) , 0.0_8)
!          
!          do p = 1 , PropagatorTheory_instance%numberOfSpecies
!             
!             if (p==i) then
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!                      
!                      do ka = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               secondOrderDensities(i)%values(ia,ja) = secondOrderDensities(i)%values(ia,ja) &
!                                    - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               secondOrderDensities(i)%values(aa,ba) = secondOrderDensities(i)%values(aa,ba) &
!                                    + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      partialValue = 0.0_8
!                      
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ka = 1 , occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))
!                               
!                            end do
!                         end do
!                      end do
!                       
!                      do ja = 1 , occupationNumberOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ba, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      secondOrderDensities(i)%values(ia,aa) = secondOrderDensities(i)%values(ia,aa) &
!                           + partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!                
!             else
!                
!                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( p ) )
!                chargeOfSpeciesB = MolecularSystem_getCharge( p )
!!                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( p )
!                occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( p )
!                activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( p )
!                lambdaOfSpeciesB = MolecularSystem_getLambda( p )
!                virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               secondOrderDensities(i)%values(ia,ja) = secondOrderDensities(i)%values(ia,ja) &
!                                    - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end do
!                end do
!                
!                do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                   do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               secondOrderDensities(i)%values(aa,ba) = secondOrderDensities(i)%values(aa,ba) &
!                                    + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      partialValue = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      secondOrderDensities(i)%values(ia,aa) = secondOrderDensities(i)%values(ia,aa) &
!                           + 2.0_8*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!                
!             end if
!             
!          end do
!       
!          ! second order densities 1 (for alpha)
!          
!       end if
!       
!       !**************************************************************************
!       !	Storing of denominators and numerators in the corresponding vectors
!       !****
!       
!       m =0
!       
!       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	
!          
!          m=m+1          
!          
!          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
!          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
!          
!          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
!          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))
!          
!          do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!             
!             if (j==i) then ! Intraspecies factors
!                
!                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
!                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
!                
!                call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                
!                id1 = 0
!                id2 = 0
!                
!                ! factor 2ph
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
!                         auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                         auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
!                         auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                         
!                         id1 = id1 + 1
!                         
!                         selfEnergy2ph(j)%values(1,id1) = auxValue_A - auxValue_B
!                         
!                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
!                              - eigenValuesOfSpeciesA%values(ba)
!                         
!                         valueOfW = 0.0_8
!                         
!                         do ja = 1, occupationNumberOfSpeciesA
!                            do ka = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ka, ia, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                               
!                            end do
!                         end do
!                         
!                         do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                            end do
!                         end do
!                         
!                         do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                            
!                            if (k .ne. i)  then
!                               
!                               nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                  chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                  eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                  occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                  activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                  lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                  virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                  
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                        
!                                        valueOfW = valueOfW + (auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                        
!                                        valueOfW = valueOfW - (auxValue_A*auxValue_B)&
!                                             /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                        
!                                     end do
!                                  end do
!                                  
!                               end if
!                               
!                            end do
!                            
!                            selfEnergy2ph(j)%values(3,id1) = valueOfW
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, aa, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, aa, ia, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               selfEnergy2hp(j)%values(1,id2) = auxValue_A - auxValue_B
!                               
!                               selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
!                                    - eigenValuesOfSpeciesA%values(ja) 
!                               
!                               valueOfW = 0.0_8
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, aa, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, ca, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                                     
!                                  end do
!                               end do
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ka = 1, occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                  end do
!                               end do
!                               
!                               do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                  
!                                  if (k .ne. i)  then
!                                     
!                                     nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                     chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                     eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                     occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                     activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                     lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                     virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesB
!                                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           valueOfW = valueOfW + (auxValue_A*auxValue_B)&
!                                                /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                           
!                                           valueOfW = valueOfW - (auxValue_A*auxValue_B)&
!                                                /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                           
!                                        end do
!                                     end do
!                                     
!                                  end if
!                                  
!                               end do
!                               
!                               selfEnergy2hp(j)%values(3,id2) = valueOfW
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! interspecies
!                   
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                   
!                   ! paso
!                   
!                   paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                        (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                   
!                   vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   
!!                   call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                        MolecularSystem_getEigenvectors(j), auxMatrix, j, trim(nameOfSpeciesB) )
!                   
!                   auxMatrix%values = auxMatrix%values * ( chargeOfSpeciesB**2.0_8 )
!
!                   if (q==1.and.m==1) then
!                      
!                      print *,"entro a densities 2"
!                      ! Second order densities 2
!                      
!                      call Matrix_constructor(secondOrderDensities(j), int(activeOrbitalsOfSpeciesB,8),&
!                           int(activeOrbitalsOfSpeciesB,8) , 0.0_8)
!                      
!                      do ia = 1 , occupationNumberOfSpeciesB
!                         do ja = 1 , occupationNumberOfSpeciesB
!                            
!                            do ka = 1 , occupationNumberOfSpeciesB
!                               do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ka, ba, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ka, aa, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesB )
!                                     auxValue_D= auxMatrix%values(auxIndex, 1)
!                                     
!                                     secondOrderDensities(j)%values(ia,ja) = secondOrderDensities(j)%values(ia,ja) &
!                                          - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ia)&
!                                          +eigenValuesOfSpeciesB%values(ka) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesB%values(ba))&
!                                          *( eigenValuesOfSpeciesB%values(ja)&
!                                          + eigenValuesOfSpeciesB%values(ka) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesB%values(ba)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ib = 1 , occupationNumberOfSpeciesA
!                               do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ab = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     secondOrderDensities(j)%values(ia,ja) = secondOrderDensities(j)%values(ia,ja) &
!                                          - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ia)&
!                                          +eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesA%values(ab))&
!                                          *( eigenValuesOfSpeciesB%values(ja)&
!                                          + eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesA%values(ab)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                         end do
!                      end do
!                      
!                      do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            do ia = 1 , occupationNumberOfSpeciesB
!                               do ja = 1 , occupationNumberOfSpeciesB
!                                  do ca = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ja, ca, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, aa, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesB )
!                                     auxValue_D= auxMatrix%values(auxIndex, 1)
!                                     
!                                     secondOrderDensities(j)%values(aa,ba) = secondOrderDensities(j)%values(aa,ba) &
!                                          + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ia)&
!                                          +eigenValuesOfSpeciesB%values(ja) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesB%values(ca))&
!                                          *( eigenValuesOfSpeciesB%values(ia)&
!                                          + eigenValuesOfSpeciesB%values(ja) - eigenValuesOfSpeciesB%values(ca) - eigenValuesOfSpeciesB%values(ba)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ab = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ia = 1 , occupationNumberOfSpeciesB
!                                  do ib = 1 , occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     secondOrderDensities(j)%values(aa,ba) = secondOrderDensities(j)%values(aa,ba) &
!                                          +  (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ia)&
!                                          + eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesA%values(ab))&
!                                          *( eigenValuesOfSpeciesB%values(ia)&
!                                          + eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(ba) - eigenValuesOfSpeciesA%values(ab)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                         end do
!                      end do
!                      
!                      do ia = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            partialValue = 0.0_8
!                            
!                            do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ja = 1 , occupationNumberOfSpeciesB
!                                  do ka = 1 , occupationNumberOfSpeciesB
!                                       
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, ba, ka, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ka, ba, ja, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesB )
!                                     auxValue_D= auxMatrix%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(ja)&
!                                          +eigenValuesOfSpeciesB%values(ka) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesB%values(ba))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ab = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ja = 1 , occupationNumberOfSpeciesB
!                                  do ib = 1 , occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ja, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue - 2.0_8*(auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ja)&
!                                          +eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesA%values(ab))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ja = 1 , occupationNumberOfSpeciesB
!                               do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ca = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, aa, ca, ja, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, aa, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesB )
!                                     auxValue_D= auxMatrix%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(ja)&
!                                          +eigenValuesOfSpeciesB%values(ia) - eigenValuesOfSpeciesB%values(ba) - eigenValuesOfSpeciesB%values(ca))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ib = 1 , occupationNumberOfSpeciesA
!                               do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ab = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ba, aa, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A= auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue + 2.0_8*(auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ia)&
!                                          +eigenValuesOfSpeciesA%values(ib) - eigenValuesOfSpeciesB%values(ba) - eigenValuesOfSpeciesA%values(ab))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            secondOrderDensities(j)%values(ia,aa) = secondOrderDensities(j)%values(ia,aa) &
!                                 + partialValue/(eigenValuesOfSpeciesB%values(ia) - eigenValuesOfSpeciesB%values(aa))
!                            
!                         end do
!                      end do
!                      
!                      ! Second order densities 2
!
!                   end if
!
!                   call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                   call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                   
!                   id1 = 0
!                   id2 = 0
!                   
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                            
!                            id1 = id1 + 1
!                            
!                            selfEnergy2ph(j)%values(1,id1) = auxValue_A
!                            
!                            selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
!                                 - eigenValuesOfSpeciesB%values(ab)
!                            
!                            valueOfW = 0.0_8
!                            
!                            do ia = 1, occupationNumberOfSpeciesA
!                               do jb = 1, occupationNumberOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, jb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib,  bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!                            
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_A = auxMatrix2(i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ba, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_B = auxMatrix2(i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, bb, ib, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, jb, bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!
!                            selfEnergy2ph(j)%values(3,id1) = valueOfW
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                            
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesB%values(ib)
!                            
!                            valueOfW = 0.0_8
!
!                            do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                  
!                               valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ia, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B = auxMatrix2(i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ab, ib, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ab, jb, bb, ib, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, jb, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         selfEnergy2hp(j)%values(3,id2) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!             end if
!             
!          end do
!
!          ! Initial guess
!          koopmans = eigenValuesOfSpeciesA%values(pa)
!          newOmega = koopmans
!          lastOmega = 0.0_8
!          
!          ni = 0
!          limit = 50
!          residual = 1.0_8
!
!          ! Calculation of second order pole
!
!          do while ((residual>0.001_8).or.(limit.lt.ni))
!             
!             ni = ni + 1
!             
!             lastOmega = newOmega
!             selfEnergy = lastOmega - koopmans
!             selfEnergyDerivative = 1.0_8
!             
!             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                
!                if (j==i) then ! Intraspecies term
!                   
!                   id1=0
!                   id2=0
!                   
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            
!                            id1 = id1 + 1
!                            
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                                                           
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! Interspecies term
!
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                   id1 = 0
!                   id2 = 0
!
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            id1 = id1 + 1
!
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!
!                            a1 = selfEnergy2hp(j)%values(1,id2)
!                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                            
!                            selfEnergy = selfEnergy -a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                end if
!                
!             end do
!             
!             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!             
!             residual = abs(newOmega-lastOmega)
!             
!             print *,"iteration",ni,"newOmega",newOmega,"residual",residual
!
!          end do ! while
!
!          poleStrenght = 1.0_8/selfEnergyDerivative
!
!          ! Storing corrections
!
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)=real(pa,8)
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3)=27.211396_8 * newOmega
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)=poleStrenght
!
!          print *,"----------------------------------------------------------------"
!          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)),&
!               " of species: ",nameOfSpeciesA
!          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T5,A29,F8.4,A7,I2,A12)") "Optimized second order pole: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               " after ",ni," iterations."
!          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!
!          ! Calculation of third order poles
!
!          ! factors for different algorythms
!
!          factors(:,:) = 0.0_8
!          thirdOrderResults(:,:) = 0.0_8
!          constantSelfEnergy = 0.0_8
!
!          ! o=1 P3
!          thirdOrderMethods(1)="P3"
!          ! o=2 EP3
!          thirdOrderMethods(2)="EP3"
!          ! o=3 OVGF version A
!          thirdOrderMethods(3)="OVGF A"
!          ! o=4 OVGF version B
!          thirdOrderMethods(4)="OVGF B"
!          ! o=5 OVGF version C
!          thirdOrderMethods(5)="OVGF C"
!
!          do o = 1 , 5 ! Options for third order
!
!             ! Initial guess             
!             koopmans = eigenValuesOfSpeciesA%values(pa)
!             newOmega = koopmans
!             lastOmega = 0.0_8
!             
!             ni = 0
!             limit = 15
!             residual = 1.0_8
!
!             ! factor for W
!             
!             fW = 2.0_8
!             fI = 1.0_8
!             if (o==1) fW=1.0_8
!             if (o==1) fI=0.0_8
!             threshold=0.001_8
!             if (o==2) threshold=0.0005_8
!             
!             ! NR procedure
!
!             do while ((residual>threshold))
!                
!                ni = ni + 1
!                
!                lastOmega = newOmega
!                selfEnergy = lastOmega - koopmans
!                selfEnergyDerivative = 1.0_8
!                s2hp = 0.0_8
!                s2ph = 0.0_8
!                W2hp = 0.0_8
!                W2ph = 0.0_8
!                U2hp = 0.0_8
!                U2ph = 0.0_8
!                
!                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                   
!                   if (j==i) then ! Intraspecies term
!                      
!                      id1=0
!                      id2=0
!
!                      ! Diagram 2ph intraspecies
!
!                      ! if (o==1 .and. ni==2) then
!                         
!                      !    do ra = 1, activeOrbitalsOfSpeciesA
!                      !       do sa = 1, activeOrbitalsOfSpeciesA
!                               
!                      !          auxIndex = IndexMap_tensorR4ToVector(pa, pa, ra, sa, activeOrbitalsOfSpeciesA )
!                      !          auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!
!                      !          auxIndex = IndexMap_tensorR4ToVector(pa, sa, ra, pa, activeOrbitalsOfSpeciesA )
!                      !          auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                      !          constantSelfEnergy = constantSelfEnergy + (auxValue_A-auxValue_B)*secondOrderDensities(j)%values(ra,sa)
!                               
!                      !       end do
!                      !    end do
!                         
!                      ! end if
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               id1 = id1 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!
!                               if ( (.not.paso1).or.(o/=1)) then
!                                  
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do da = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ia, da, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, da, ia, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, aa, da, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, ba, da, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ia) &
!                                             - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(da)
!                                        
!                                        valueOfU = valueOfU + 0.5_8*a2/c
!                                        valueOfdU = valueOfdU - 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueOfdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                     
!                                     if (k .ne. i)  then
!                                        
!                                        nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                        chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                        eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                        occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                        activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                        lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                        virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                        
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                              
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                              
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                     end if
!                                     
!                                  end do
!                                  
!                               end if
!                               
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      if (paso1.and.(o==1)) subW=0.0_8
!                      if (paso1.and.(o==1)) subdW=0.0_8                      
!
!                      s2ph = s2ph + 0.5_8*sub2
!                      W2ph = w2ph + 0.5_8*(fW*subW)/(1.0_8-factors(1,o))                       
!                      U2ph = U2ph + 0.5_8*(subU)/(1.0_8-factors(1,o))                       
!
!                      selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(1,o)) )                       
!
!                      selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(1,o)) )                     
!                      ! Diagram 2hp intraspecies
!                        
!                      if (occupationNumberOfSpeciesA > 1) then
!                         
!                         sub2 = 0.0_8
!                         subW = 0.0_8 
!                         subU = 0.0_8 
!                         subd2 = 0.0_8
!                         subdW = 0.0_8
!                         subdU = 0.0_8
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ia = 1 , occupationNumberOfSpeciesA
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  
!                                  id2 = id2 + 1
!                                  
!                                  valueOfU = 0.0_8
!                                  valueOfdU = 0.0_8
!                                  
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do la = 1 , occupationNumberOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, aa, la, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, la, aa, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ia, la, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ja, la, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(aa) &
!                                             - eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(la)
!                                        
!                                        valueOfU = valueOfU - 0.5_8*a2/c
!                                        valueOfdU = valueOfdU + 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueOfdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueofdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                     
!                                     if (k .ne. i)  then
!                                        
!                                        nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                        chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                        eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                        occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                        activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                        lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                        virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                        
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                              
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ja, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ia)
!                                              
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                     end if
!                                     
!                                  end do
!                                  
!                                  a1 = selfEnergy2hp(j)%values(1,id2)
!                                  a2 = selfEnergy2hp(j)%values(3,id2)
!                                  b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                  
!                                  sub2 = sub2 + (a1**2.0_8)/b
!                                  subW = subW + (a1*a2)/b
!                                  subU = subU + (a1*valueOfU)/b
!                                  
!                                  subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                                  subdW = subdW + (a1*a2)/(b**2.0_8)
!                                  subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                  
!                               end do
!                            end do
!                         end do
!                         
!                         s2hp = s2hp + 0.5_8*sub2
!                         W2hp = w2hp + 0.5_8*(fW*subW)/(1.0_8-factors(2,o))                       
!                         U2hp = U2hp + 0.5_8*(subU)/(1.0_8-factors(2,o))                       
!                         
!                         selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(2,o)) )                       
!                         
!                         selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(2,o)) )                     
!                      end if
!                      
!                   else ! Interspecies term
!                      
!                      nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                      chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                      eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                      occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                      activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                      lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                      virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                      
!                      paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                           (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                      
!!                      call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                           MolecularSystem_getEigenvectors(j), auxMatrix, j, trim(nameOfSpeciesB) )
!                      
!                      auxMatrix%values = auxMatrix%values * ( chargeOfSpeciesB**2.0_8 )
!                      
!                      id1 = 0
!                      id2 = 0
!                      
!                      if (o==1 .and. ni==2) then
!
!                         do rb = 1, activeOrbitalsOfSpeciesB
!                            do sb = 1, activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, pa, rb, sb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                               
!                               constantSelfEnergy = constantSelfEnergy + auxValue_A*secondOrderDensities(j)%values(rb,sb)
!                               
!                            end do
!                         end do
!                      
!                      end if
!                                            
!                      ! Diagram 2ph interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               id1 = id1 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               if ( (.not.paso2).or.(o/=1) ) then
!                                  
!                                  do jb = 1 , occupationNumberOfSpeciesB
!                                     
!                                     do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                        auxValue_C= auxMatrix%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(aa)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ba, aa, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, aa, ab, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_C = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU + (a2/c)
!                                        valueofdU = valueOfdU - (a2/(c**2.0_8))
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do jb = 1 , occupationNumberOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                        
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                               end if
!                               
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      if (paso2.and.(o==1)) subW=0.0_8
!                      if (paso2.and.(o==1)) subdW=0.0_8                      
!
!                      s2ph = s2ph + sub2
!                      W2ph = w2ph + (fW*subW)/(1.0_8-factors(1,o))                       
!                      U2ph = U2ph + (subU)/(1.0_8-factors(1,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(1,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(1,o)) )                  
!                      
!                      ! Diagram 2hp interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               id2 = id2 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ia)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, bb, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               a2 = selfEnergy2hp(j)%values(3,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      s2hp = s2hp + sub2
!                      W2hp = w2hp + (fW*subW)/(1.0_8-factors(2,o))                       
!                      U2hp = U2hp + (subU)/(1.0_8-factors(2,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(2,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(2,o)) )                  
!                      
!                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                         
!                         id1=0
!                         id2=0
!                         
!                         if (k.ne.i .and. k.ne.j)  then
!                            
!                            print *,"entro al manolito",k
!                            
!                            nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                            chargeOfSpeciesC = MolecularSystem_getCharge( k )
!!                            eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( k )
!                            occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( k )
!                            activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( k )
!                            lambdaOfSpeciesC = MolecularSystem_getLambda( k )
!                            virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                            
!!                           call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                                 MolecularSystem_getEigenVectors(j), MolecularSystem_getEigenVectors(k), &
!!                                 auxMatrix3, j, nameOfSpeciesB, k, nameOfSpeciesC )
!
!                            auxMatrix3%values = auxMatrix3%values * (chargeOfSpeciesB*chargeOfSpeciesC)                            
!
!                            if (q==1 .and. o==1 .and. ni==1) then
!                               
!                               ! Second order densities 3
!
!                               print *,"entro a densities 3"
! 
!                               do ia = 1 , occupationNumberOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesB
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesC
!                                        do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_A= auxMatrix3%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B= auxMatrix3%values(auxIndex, 1)
!                                              
!                                              secondOrderDensities(j)%values(ia,ja) = secondOrderDensities(j)%values(ia,ja) &
!                                                   - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ia)&
!                                                   +eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesC%values(ab))&
!                                                   *( eigenValuesOfSpeciesB%values(ja)&
!                                                   + eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesC%values(ab)))
!                                              
!                                           end do
!                                        end do
!                                     end do
!                                     
!                                  end do
!                               end do
!                               
!                               do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     do ab = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                        do ia = 1 , occupationNumberOfSpeciesB
!                                           do ib = 1 , occupationNumberOfSpeciesC
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_A= auxMatrix3%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B= auxMatrix3%values(auxIndex, 1)
!                                              
!                                              secondOrderDensities(j)%values(aa,ba) = secondOrderDensities(j)%values(aa,ba) &
!                                                   + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ia)&
!                                                   + eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesC%values(ab))&
!                                                   *( eigenValuesOfSpeciesB%values(ia)&
!                                                   + eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(ba) - eigenValuesOfSpeciesC%values(ab)))
!                                              
!                                           end do
!                                        end do
!                                     end do
!                                     
!                                  end do
!                               end do
!                               
!                               do ia = 1 , occupationNumberOfSpeciesB
!                                  do aa = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     partialValue = 0.0_8
!                                     
!                                     do ab = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                        do ja = 1 , occupationNumberOfSpeciesB
!                                           do ib = 1 , occupationNumberOfSpeciesC
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_A= auxMatrix3%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B= auxMatrix3%values(auxIndex, 1)
!                                              
!                                              partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ja)&
!                                                   +eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(aa) - eigenValuesOfSpeciesC%values(ab))
!                                              
!                                           end do
!                                        end do
!                                     end do
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesC
!                                        do ba = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(ba, aa, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_A= auxMatrix3%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B= auxMatrix3%values(auxIndex, 1)
!                                              
!                                              partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ia)&
!                                                   +eigenValuesOfSpeciesC%values(ib) - eigenValuesOfSpeciesB%values(ba) - eigenValuesOfSpeciesC%values(ab))
!                                              
!                                           end do
!                                        end do
!                                     end do
!                                     
!                                     secondOrderDensities(j)%values(ia,aa) = secondOrderDensities(j)%values(ia,aa) &
!                                          + 2.0_8*partialValue/(eigenValuesOfSpeciesB%values(ia) - eigenValuesOfSpeciesB%values(aa))
!                                  
!                                  end do
!                               end do
!                               
!                               ! end of second order densities 3
!
!                            end if
!                            
!                            subU=0.0_8
!                            subdU=0.0_8                            
!                            
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     id1 = id1 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, aa, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_B = auxMatrix3%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ic) &
!                                                - eigenValuesOfSpeciesC%values(ac) - eigenValuesOfSpeciesA%values(aa)
!                                           
!                                           
!                                           valueOfU = valueOfU + a2/c
!                                           valueofdU = valueOfdU - a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2ph(j)%values(1,id1)
!                                     b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                                     
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            U2ph = U2ph + (subU)/(1.0_8-factors(1,o))                                                   
!                            
!                            selfEnergy = selfEnergy - (subU)/(1.0_8-factors(1,o))                       
!                            
!                            selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(1,o))                                           
!                            
!                            subU=0.0_8
!                            subdU=0.0_8
!                            
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     
!                                     id2 = id2 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ia, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_A = auxMatrix2(k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_B = auxMatrix3%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ac) &
!                                                - eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesA%values(ia)
!                                           
!                                           valueOfU = valueOfU - a2/c
!                                           valueofdU = valueOfdU + a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2hp(j)%values(1,id2)
!                                     b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                     
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!
!                            U2hp = U2hp + (subU)/(1.0_8-factors(2,o))    
!                            
!                            selfEnergy = selfEnergy - (subU)/(1.0_8-factors(2,o))                       
!                            
!                            selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(2,o))                   
!                            
!                         end if
!                         
!                      end do
!                      
!                   end if
!                   
!                end do
!
!                selfEnergy = selfEnergy - fI*constantSelfEnergy/(1.0_8-factors(3,o))
!                
!                newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!                
!                residual = abs(newOmega-lastOmega)
!                
!                print *,"ni",ni,"newOmega",newOmega,"residual",residual
!                
!             end do ! while
!
!             if (o==2) then
!                
!                ! OVGF version A
!                factors(1,3) = (W2hp+W2ph)/(s2hp+s2ph)
!                factors(2,3) = factors(1,3)
!                factors(3,3) = factors(1,3)
!
!                ! OVGF version B
!                factors(1,4) = W2ph/s2ph
!                factors(2,4) = W2hp/s2hp
!                factors(3,4) = 0.0_8
!
!                ! OVGF version C
!                factors(1,5) = (factors(2,4)*(U2hp+W2hp)+factors(1,4)*(U2ph+W2ph))/(W2hp+W2ph+U2hp+U2ph)
!                factors(2,5) = factors(1,5)
!                factors(3,5) = factors(1,5)
!
!             end if
!
!             poleStrenght = 1.0_8/(selfEnergyDerivative)
!             thirdOrderResults(1,o) = 27.211396_8 * newOmega
!             thirdOrderResults(2,o) = poleStrenght
!
!             print *,"value of o:",o
!             print *,"constant self-energy:", constantSelfEnergy
!             print *,"2hp(2):",s2hp*27.211396_8,"2ph(2):",s2ph*27.211396_8
!             print *,"W 2hp(3):",W2hp*27.211396_8,"W 2ph(3):",W2ph*27.211396_8
!             print *,"U 2hp(3):",U2hp*27.211396_8,"U 2ph(3):",U2ph*27.211396_8
!             print *,"factor 2hp:",factors(2,o),"factor 2ph:",factors(1,o)
!             write (*,"(T5,A10,A10,A6,F8.4,A7,I2,A12)") "Optimized ",thirdOrderMethods(o),"pole: ",newOmega*27.211396_8," after ",ni," iterations."
!             write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!             print *,"----------------------------------------------------------------"
!             
!             
!          end do ! options of third order
!
!          ! printing results for one spin-orbital
!
!          write (*,"(T5,A55,I2,A13,A8)") "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL:",&
!               int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1))," OF SPECIES:",nameOfSpeciesA
!          write (*, "(T5,A45)") "--------------------------------------------"
!          write (*, "(T10,A10,A10,A10)") " Method ","BE (eV)","Pole S."
!          write (*, "(T5,A45)") "--------------------------------------------"
!          write (*,"(T10,A10,F10.3)") "KT        ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T10,A10,F10.3,F10.4)") "EP2       ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)
!          do o=1,5
!
!             write (*,"(T10,A10,F10.3,F10.4)") thirdOrderMethods(o),thirdOrderResults(1,o),thirdOrderResults(2,o)
!             
!          end do
!          write (*, "(T5,A45)") "--------------------------------------------"
!
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,5)=27.211396_8 * newOmega
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,6)=poleStrenght
!                              
!          ! call Matrix_destructor(auxMatrix2(:))          
!!          call TransformIntegrals_destructor( repulsionTransformer )
!          
!       end do
!       
!    end do
!
!    do i = 1, PropagatorTheory_instance%numberOfSpecies
!
!       print *,"Second order densities for species:",i
!       print *,secondOrderDensities(i)%values(:,:)
!
!    end do
!    
!    !!
!    !!************************************************************************************************
!    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
!    print *,"***************************************************************"
!  end subroutine PropagatorTheory_thirdOrderCorrection2
!
!  !**
!  ! @brief Evaluate partial third order, EP3 and OVGF poles
!  ! version 3
!  !**
!
!  subroutine PropagatorTheory_thirdOrderCorrection3()
!    implicit NONE
!    
!    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
!    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
!    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
!    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
!    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
!    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
!    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
!    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
!    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
!    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
!    integer :: i, j, k ! counters for species
!    integer :: m, n, o, p, q, r, ni, nc, limit, id1, id2 ! auxiliar counters
!    integer :: speciesAID, speciesBID, speciesCID
!    integer :: species1ID, species2ID
!    integer :: electronsID
!    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
!    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
!    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
!    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
!    integer :: activeOrbitalsOfSpeciesC
!    integer(8) :: vectorSize1, vectorSize2 !!! Sizes for diagrams
!    integer(8) :: auxIndex
!    integer(4) :: errorNum
!    character(10) :: thirdOrderMethods(6)
!    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
!    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
!    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
!    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
!    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!!    type(TransformIntegrals) :: repulsionTransformer
!    type(Matrix),allocatable :: auxMatrix2(:,:), selfEnergy2hp(:), selfEnergy2ph(:)
!    type(Matrix),allocatable :: secondOrderDensities(:)
!    type(Matrix) :: diagram_A, diagram_B, auxMatrix, auxMatrix3
!    type(Matrix) :: partialMO1, partialMO2    
!    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
!    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
!    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
!    real(8) :: valueOfW, valueOfU, valueOfdU, sub2, subW, subU, subd2, subdW, subdU
!    real(8) :: lastOmega, newOmega, residual, threshold, selfEnergy, selfEnergyDerivative, koopmans 
!    real(8) :: a1, a2, b, c, d, poleStrenght, partialValue, partialValue2, initialValue
!    real(8) :: factors(3,6), fW, fI, constantSelfEnergy,thirdOrderResults(2,5)
!    real(8) :: s2hp, s2ph, W2hp, W2ph, U2hp, U2ph
!    logical :: paso1, paso2
!    ! *******************************************************************************************
!    ! Determinate the numerators and denominators of the second Oder propapator 
!    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
!    
!    print *,"*******************************************************************"
!    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"
!
!    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified
!
!    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
!    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies,PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(secondOrderDensities)) deallocate(secondOrderDensities)
!    allocate(secondOrderDensities(PropagatorTheory_instance%numberOfSpecies))
!
!    !!! Defining for which species the correction will be applied
!    
!    if (CONTROL_instance%IONIZE_SPECIE /= "NONE") then
!       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
!       species2ID= species1ID
!       m=1
!    else
!       species1ID=1
!       species2ID=PropagatorTheory_instance%numberOfSpecies
!       m = species2ID
!    end if
!
!    if (allocated(PropagatorTheory_instance%thirdOrderCorrections)) deallocate(PropagatorTheory_instance%thirdOrderCorrections)
!    allocate(PropagatorTheory_instance%thirdOrderCorrections(m))
!
!    ! Storing transformed integrals !!!! We need a more efficient algorithm to do this
!    
!!    call TransformIntegrals_constructor( repulsionTransformer )
!    
!    do p = 1 , PropagatorTheory_instance%numberOfSpecies
!
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( p ) )
!
!       do n = 1 , PropagatorTheory_instance%numberOfSpecies
!          
!          if (n==p) then
!             
!!             call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p,p), p, trim(nameOfSpeciesA) )
!             
!             auxMatrix2(p,p)%values = auxMatrix2(p,p)%values * MolecularSystem_getCharge( p ) &
!                  * MolecularSystem_getCharge( p )
!             
!          else
!             
!             nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( n ) )
!             
!!             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                  MolecularSystem_getEigenVectors(p), MolecularSystem_getEigenVectors(n), &
!!                  auxMatrix2(p,n), p, nameOfSpeciesA, n, nameOfSpeciesB )
!             
!             auxMatrix2(p,n)%values = auxMatrix2(p,n)%values * MolecularSystem_getCharge( n ) &
!                  * MolecularSystem_getCharge( p )
!             
!          end if
!          
!       end do
!
!    end do
!
!    ! Start loop for species
!    
!    q = 0
!    
!    do i = species1ID , species2ID
!       
!       q = q + 1
!       
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
!       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!!       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
!       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )
!       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i )
!       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
!       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA
!       
!       ! paso
!       
!       paso1=(nameOfSpeciesA=="e-ALPHA".or.nameOfSpeciesA=="e-BETA")
!       
!       ! Defining the number of orbitals !!! Insert a parameter for the else option
!       
!       if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
!          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
!          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
!          n = 1
!       else if (CONTROL_instance%IONIZE_SPECIE /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
!          n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
!       else
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
!          n = 2
!       end if
!
!       call Matrix_constructor(PropagatorTheory_instance%thirdOrderCorrections(q), int(n,8), 8, 0.0_8)
!
!       !**************************************************************************
!       !	Storing of denominators and numerators in the corresponding vectors
!       !****
!       
!       m =0
!       
!       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	
!
!          m=m+1          
!          
!          ! calculation of constant self energy
!          
!          constantSelfEnergy = 0.0_8
!          
!          do p = 1 , PropagatorTheory_instance%numberOfSpecies
!             
!             if (p==i) then
!
!                print *,"entro al if"                
!
!                ! alpha-alpha-alpha
!
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ka = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*(auxValue_E-auxValue_F)                       
!
!                   end do
!                end do
!                
!                do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                   do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8                      
!                      
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ka = 1 , occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do ja = 1 , occupationNumberOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ba, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy &
!                           + (auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma after a-a-a:", constantSelfEnergy
!
!             else
!
!                print *,"entro al else"                
!
!                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( p ) )
!                chargeOfSpeciesB = MolecularSystem_getCharge( p )
!!                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( p )
!                occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( p )
!                activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( p )
!                lambdaOfSpeciesB = MolecularSystem_getLambda( p )
!                virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                
!                ! alpha-beta-beta 
!
!                do ib = 1 , occupationNumberOfSpeciesB
!                   do jb = 1 , occupationNumberOfSpeciesB
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                      auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do kb = 1 , occupationNumberOfSpeciesB
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
!                                    +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))&
!                                    *( eigenValuesOfSpeciesB%values(jb)&
!                                    + eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*auxValue_E
!                      
!                   end do
!                end do
!                
!                do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                   do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                      auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8                      
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
!                                    +eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(cb))&
!                                    *( eigenValuesOfSpeciesB%values(ib)&
!                                    + eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(cb) - eigenValuesOfSpeciesB%values(bb)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*auxValue_E 
!                      
!                   end do
!                end do
!                
!                do ib = 1 , occupationNumberOfSpeciesB
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                      auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do kb = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, kb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, kb, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
!                                    +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do jb = 1 , occupationNumberOfSpeciesB
!                         do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(bb, ab, cb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(bb, jb, cb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesB%values(cb))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy &
!                           + auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma after a-b-b:", constantSelfEnergy
!
!                ! alpha-alpha-beta
!
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                   do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ib =  1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy = constantSelfEnergy &
!                           + 2.0_8*(auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma after a-a-b:", constantSelfEnergy
!
!                ! Three species
!
!                ! alpha-beta-alpha
!
!                ! alpha-beta-gamma
!
!                do r = 1, PropagatorTheory_instance%numberOfSpecies
!
!                   if ( r /= p) then
!
!                      print *,"entro a r diferente de p"
!
!                      nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( r ) )
!                      chargeOfSpeciesC = MolecularSystem_getCharge( r )
!!                      eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( r )
!                      occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( r )
!                      activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( r )
!                      lambdaOfSpeciesC = MolecularSystem_getLambda( r )
!                      virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                            
!                            partialValue = 0.0_8
!                            
!                            do ic = 1 , occupationNumberOfSpeciesC
!                               do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue &
!                                          - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
!                                          *( eigenValuesOfSpeciesB%values(jb)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy = constantSelfEnergy + partialValue*auxValue_E
!                            
!                         end do
!                      end do
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                            
!                            partialValue = 0.0_8
!                            
!                            do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                               do ib = 1 , occupationNumberOfSpeciesB
!                                  do ic = 1 , occupationNumberOfSpeciesC
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue &
!                                          + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
!                                          *( eigenValuesOfSpeciesB%values(ib)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy = constantSelfEnergy + partialValue*auxValue_E
!                            
!                         end do
!                      end do
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                            
!                            partialValue = 0.0_8
!                            
!                            do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  do ic =  1 , occupationNumberOfSpeciesC
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, jb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(jb)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ic = 1 , occupationNumberOfSpeciesC
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ab, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                     auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!                                     
!                                     partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ib)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy = constantSelfEnergy &
!                                 + 2.0_8*auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
!                            
!                         end do
!                      end do
!
!                   end if
!
!                end do
!
!                print *,"constant sigma after a-b-a:", constantSelfEnergy
!                
!             end if
!                
!          end do
!
!          ! end of constant self-energy calculation
!          
!          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
!          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
!          
!          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
!          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))
!
!          do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!             
!             if (j==i) then ! Intraspecies factors
!                
!                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
!                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
!                
!                call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                
!                id1 = 0
!                id2 = 0
!                
!                ! factor 2ph
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
!                         auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
!                         auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
!                         auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
!                         
!                         id1 = id1 + 1
!                         
!                         selfEnergy2ph(j)%values(1,id1) = auxValue_A - auxValue_B
!                         
!                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
!                              - eigenValuesOfSpeciesA%values(ba)
!                         
!                         valueOfW = 0.0_8
!                         
!                         do ja = 1, occupationNumberOfSpeciesA
!                            do ka = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ka, ia, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                               
!                            end do
!                         end do
!                         
!                         do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                            end do
!                         end do
!                         
!                         do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                            
!                            if (k .ne. i)  then
!                               
!                               nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                               chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                               eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                               occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                               activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                               lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                               virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                               
!                               do ib = 1 , occupationNumberOfSpeciesB
!                                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + (auxValue_A*auxValue_B)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW - (auxValue_A*auxValue_B)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                     
!                                  end do
!                               end do
!                               
!                            end if
!                            
!                         end do
!                         
!                         selfEnergy2ph(j)%values(3,id1) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!                if (occupationNumberOfSpeciesA > 1) then
!                   
!                   ! factor 2hp
!                   
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            
!                            id2 = id2 + 1
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ia, aa, ja, activeOrbitalsOfSpeciesA )
!                            auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ja, aa, ia, activeOrbitalsOfSpeciesA )
!                            auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
!                            
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A - auxValue_B
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesA%values(ja) 
!                               
!                               valueOfW = 0.0_8
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, aa, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, ca, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                                     
!                                  end do
!                               end do
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ka = 1, occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                  end do
!                               end do
!                               
!                               do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                  
!                                  if (k .ne. i)  then
!                                     
!                                     nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                     chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                     eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                     occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                     activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                     lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                     virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                     
!                                     do ib = 1 , occupationNumberOfSpeciesB
!                                        do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           
!                                           valueOfW = valueOfW + (auxValue_A*auxValue_B)&
!                                                /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           
!                                           valueOfW = valueOfW - (auxValue_A*auxValue_B)&
!                                                /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                                - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                           
!                                        end do
!                                     end do
!                                     
!                                  end if
!                                  
!                               end do
!                               
!                               selfEnergy2hp(j)%values(3,id2) = valueOfW
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! interspecies
!                   
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                   
!                   ! paso
!                   
!                   paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                        (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                   
!                   vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   
!                   call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                   call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                   
!                   id1 = 0
!                   id2 = 0
!                   
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                            
!                            id1 = id1 + 1
!                            
!                            selfEnergy2ph(j)%values(1,id1) = auxValue_A
!                            
!                            selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
!                                 - eigenValuesOfSpeciesB%values(ab)
!                            
!                            valueOfW = 0.0_8
!                            
!                            do ia = 1, occupationNumberOfSpeciesA
!                               do jb = 1, occupationNumberOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, jb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib,  bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!                            
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ba, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, ab, ib, bb, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, bb, ib, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, jb, bb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!
!                            selfEnergy2ph(j)%values(3,id1) = valueOfW
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                            auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                            
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesB%values(ib)
!                            
!                            valueOfW = 0.0_8
!
!                            do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  
!                               valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ia, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ab, ib, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ab, jb, bb, ib, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, jb, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         selfEnergy2hp(j)%values(3,id2) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!             end if
!             
!          end do
!
!          ! Initial guess
!          koopmans = eigenValuesOfSpeciesA%values(pa)
!          newOmega = koopmans
!          lastOmega = 0.0_8
!          
!          ni = 0
!          limit = 50
!          residual = 1.0_8
!
!          ! Calculation of second order pole
!
!          do while ((residual>0.001_8).or.(limit.lt.ni))
!             
!             ni = ni + 1
!             
!             lastOmega = newOmega
!             selfEnergy = lastOmega - koopmans
!             selfEnergyDerivative = 1.0_8
!             
!             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                
!                if (j==i) then ! Intraspecies term
!                   
!                   id1=0
!                   id2=0
!                   
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            
!                            id1 = id1 + 1
!                            
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                                                           
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! Interspecies term
!
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                   id1 = 0
!                   id2 = 0
!
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            id1 = id1 + 1
!
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!
!                            a1 = selfEnergy2hp(j)%values(1,id2)
!                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                            
!                            selfEnergy = selfEnergy -a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                end if
!                
!             end do
!             
!             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!             
!             residual = abs(newOmega-lastOmega)
!             
!             print *,"iteration",ni,"newOmega",newOmega,"residual",residual
!
!          end do ! while
!
!          poleStrenght = 1.0_8/selfEnergyDerivative
!
!          ! Storing corrections
!
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)=real(pa,8)
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3)=27.211396_8 * newOmega
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)=poleStrenght
!
!          print *,"----------------------------------------------------------------"
!          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)),&
!               " of species: ",nameOfSpeciesA
!          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T5,A29,F8.4,A7,I2,A12)") "Optimized second order pole: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               " after ",ni," iterations."
!          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!
!          ! Calculation of third order poles
!
!          ! factors for different algorythms
!
!          factors(:,:) = 0.0_8
!          thirdOrderResults(:,:) = 0.0_8
!
!          ! o=1 P3
!          thirdOrderMethods(1)="P3"
!          ! o=2 EP3
!          thirdOrderMethods(2)="EP3"
!          ! o=3 OVGF version A
!          thirdOrderMethods(3)="OVGF A"
!          ! o=4 OVGF version B
!          thirdOrderMethods(4)="OVGF B"
!          ! o=5 OVGF version C
!          thirdOrderMethods(5)="OVGF C"
!          ! o=6 
!          thirdOrderMethods(6)="REN-P3"
!          
!          do o = 1 , 6 ! Options for third order
!
!             ! Initial guess             
!             koopmans = eigenValuesOfSpeciesA%values(pa)
!             newOmega = koopmans
!             lastOmega = 0.0_8
!             
!             ni = 0
!             limit = 15
!             residual = 1.0_8
!
!             ! factor for W
!             
!             fW = 2.0_8
!             fI = 1.0_8
!             if (o==1 .or. o==6) fW=1.0_8
!             if (o==1 .or. o==6) fI=0.0_8
!             threshold=0.001_8
!             if (o==1 .or. o==2) threshold=0.00005_8
!             
!             ! NR procedure
!
!             do while ((residual>threshold))
!                
!                ni = ni + 1
!                
!                lastOmega = newOmega
!                selfEnergy = lastOmega - koopmans
!                selfEnergyDerivative = 1.0_8
!                s2hp = 0.0_8
!                s2ph = 0.0_8
!                W2hp = 0.0_8
!                W2ph = 0.0_8
!                U2hp = 0.0_8
!                U2ph = 0.0_8
!                
!                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                   
!                   if (j==i) then ! Intraspecies term
!                      
!                      id1=0
!                      id2=0
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!
!                      ! 2ph terms
!
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               id1 = id1 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!
!                               if ( (.not.paso1).or.(o/=1)) then
!                                  
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do da = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ia, da, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, da, ia, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, aa, da, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, ba, da, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ia) &
!                                             - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(da)
!                                        
!                                        valueOfU = valueOfU + 0.5_8*a2/c
!                                        valueOfdU = valueOfdU - 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueOfdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                     
!                                     if (k .ne. i)  then
!                                        
!                                        nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                        chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                        eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                        occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                        activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                        lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                        virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                        
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                              
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, ba, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                              
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                     end if
!                                     
!                                  end do
!                                  
!                               end if
!                               
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!
!                      if (paso1.and.(o==1.or.o==6)) subW=0.0_8
!                      if (paso1.and.(o==1.or.o==6)) subdW=0.0_8                      
!
!                      s2ph = s2ph + 0.5_8*sub2
!                      W2ph = w2ph + 0.5_8*(fW*subW)/(1.0_8-factors(1,o))                       
!                      U2ph = U2ph + 0.5_8*(subU)/(1.0_8-factors(1,o))                       
!
!                      selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(1,o)) )                       
!
!                      selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(1,o)) )                     
!
!                      ! Diagram 2hp intraspecies
!                        
!                      if (occupationNumberOfSpeciesA > 1) then
!                         
!                         sub2 = 0.0_8
!                         subW = 0.0_8 
!                         subU = 0.0_8 
!                         subd2 = 0.0_8
!                         subdW = 0.0_8
!                         subdU = 0.0_8
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ia = 1 , occupationNumberOfSpeciesA
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  
!                                  id2 = id2 + 1
!                                  
!                                  valueOfU = 0.0_8
!                                  valueOfdU = 0.0_8
!                                  
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do la = 1 , occupationNumberOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, aa, la, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, la, aa, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ia, la, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ja, la, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(aa) &
!                                             - eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(la)
!                                        
!                                        valueOfU = valueOfU - 0.5_8*a2/c
!                                        valueOfdU = valueOfdU + 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueOfdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueofdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                                     
!                                     if (k .ne. i)  then
!                                        
!                                        nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                                        chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                                        eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                                        occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                                        activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                                        lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                                        virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                                        
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                              
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                              
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              auxIndex = IndexMap_tensorR4ToVector(ja, aa, ab, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                              
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ia)
!                                              
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                     end if
!                                     
!                                  end do
!                                  
!                                  a1 = selfEnergy2hp(j)%values(1,id2)
!                                  a2 = selfEnergy2hp(j)%values(3,id2)
!                                  b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                  
!                                  sub2 = sub2 + (a1**2.0_8)/b
!                                  subW = subW + (a1*a2)/b
!                                  subU = subU + (a1*valueOfU)/b
!                                  
!                                  subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                                  subdW = subdW + (a1*a2)/(b**2.0_8)
!                                  subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                  
!                               end do
!                            end do
!                         end do
!                         
!                         s2hp = s2hp + 0.5_8*sub2
!                         W2hp = w2hp + 0.5_8*(fW*subW)/(1.0_8-factors(2,o))                       
!                         U2hp = U2hp + 0.5_8*(subU)/(1.0_8-factors(2,o))                       
!                         
!                         selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(2,o)) )                       
!                         
!                         selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(2,o)) )                     
!                      
!                      end if
!                      
!                   else ! Interspecies term
!                      
!                      nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                      chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                      eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                      occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                      activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                      lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                      virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                      
!                      paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                           (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                      
!                      id1 = 0
!                      id2 = 0
!                      
!                      ! Diagram 2ph interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               id1 = id1 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               if ( (.not.paso2).or.(o/=1) ) then
!                                  
!                                  do jb = 1 , occupationNumberOfSpeciesB
!                                     
!                                     do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, ab, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, bb, jb, activeOrbitalsOfSpeciesB )
!                                        auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(aa)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ba, aa, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ib, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                                                                
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU + (a2/c)
!                                        valueofdU = valueOfdU - (a2/(c**2.0_8))
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do jb = 1 , occupationNumberOfSpeciesB
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, jb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                               end if
!                               
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      if (paso2.and.(o==1.or.o==6)) subW=0.0_8
!                      if (paso2.and.(o==1.or.o==6)) subdW=0.0_8                      
!
!                      s2ph = s2ph + sub2
!                      W2ph = w2ph + (fW*subW)/(1.0_8-factors(1,o))                       
!                      U2ph = U2ph + (subU)/(1.0_8-factors(1,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(1,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(1,o)) )                  
!                      
!                      ! Diagram 2hp interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               id2 = id2 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ia)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                      
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ja, bb, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, bb, ib, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               a2 = selfEnergy2hp(j)%values(3,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      s2hp = s2hp + sub2
!                      W2hp = w2hp + (fW*subW)/(1.0_8-factors(2,o))                       
!                      U2hp = U2hp + (subU)/(1.0_8-factors(2,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(2,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(2,o)) )                 
!                      
!                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                         
!                         id1=0
!                         id2=0
!                         
!                         if (k.ne.i .and. k.ne.j)  then
!                            
!                            print *,"entro al manolito",k
!                            
!                            nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                            chargeOfSpeciesC = MolecularSystem_getCharge( k )
!!                            eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( k )
!                            occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( k )
!                            activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( k )
!                            lambdaOfSpeciesC = MolecularSystem_getLambda( k )
!                            virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                                                        
!                            subU=0.0_8
!                            subdU=0.0_8                            
!                            
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     id1 = id1 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, aa, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ic) &
!                                                - eigenValuesOfSpeciesC%values(ac) - eigenValuesOfSpeciesA%values(aa)
!                                           
!                                           
!                                           valueOfU = valueOfU - a2/c
!                                           valueofdU = valueOfdU + a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2ph(j)%values(1,id1)
!                                     b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                                     
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            U2ph = U2ph + (subU)/(1.0_8-factors(1,o))                                                   
!                            
!                            selfEnergy = selfEnergy - (subU)/(1.0_8-factors(1,o))                       
!                            
!                            selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(1,o))                                                                       
!                            subU=0.0_8
!                            subdU=0.0_8
!                            
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     
!                                     id2 = id2 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ia, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                           auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
!                                           
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ac) &
!                                                - eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesA%values(ia)
!                                           
!                                           valueOfU = valueOfU + a2/c
!                                           valueofdU = valueOfdU - a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2hp(j)%values(1,id2)
!                                     b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                     
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!
!                            U2hp = U2hp + (subU)/(1.0_8-factors(2,o))    
!                            
!                            selfEnergy = selfEnergy - (subU)/(1.0_8-factors(2,o))                       
!                            
!                            selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(2,o))                   
!                            
!                         end if
!                         
!                      end do
!                      
!                   end if
!                   
!                end do
!
!                selfEnergy = selfEnergy - fI*constantSelfEnergy/(1.0_8-factors(3,o))
!                ! selfEnergy = selfEnergy - fI*constantSelfEnergy
!                
!                newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!                
!                residual = abs(newOmega-lastOmega)
!                
!                print *,"ni",ni,"newOmega",newOmega,"residual",residual
!                
!             end do ! while
!             
!             if (o==1) then
!
!                ! Renormalized P3
!                factors(1,6) = W2ph/s2ph
!                factors(2,6) = W2hp/s2hp
!                factors(3,6) = 0.0_8
!                                
!             end if
!             
!             if (o==2) then
!                
!                ! OVGF version A
!                factors(1,3) = (W2hp+W2ph)/(s2hp+s2ph)
!                factors(2,3) = factors(1,3)
!                factors(3,3) = factors(1,3)
!
!                ! OVGF version B
!                factors(1,4) = W2ph/s2ph
!                factors(2,4) = W2hp/s2hp
!                factors(3,4) = 0.0_8
!
!                ! OVGF version C
!                factors(1,5) = (factors(2,4)*(U2hp+W2hp)+factors(1,4)*(U2ph+W2ph))/(W2hp+W2ph+U2hp+U2ph)
!                factors(2,5) = factors(1,5)
!                factors(3,5) = factors(1,5)
!
!             end if
!
!             poleStrenght = 1.0_8/(selfEnergyDerivative)
!             thirdOrderResults(1,o) = 27.211396_8 * newOmega
!             thirdOrderResults(2,o) = poleStrenght
!
!             print *,"value of o:",o
!             print *,"constant self-energy:", constantSelfEnergy*27.211396_8
!             print *,"2hp(2):",s2hp*27.211396_8,"2ph(2):",s2ph*27.211396_8
!             print *,"W 2hp(3):",W2hp*27.211396_8,"W 2ph(3):",W2ph*27.211396_8
!             print *,"U 2hp(3):",U2hp*27.211396_8,"U 2ph(3):",U2ph*27.211396_8
!             print *,"factor 2hp:",factors(2,o),"factor 2ph:",factors(1,o)
!             write (*,"(T5,A10,A10,A6,F8.4,A7,I2,A12)") "Optimized ",thirdOrderMethods(o),"pole: ",newOmega*27.211396_8," after ",ni," iterations."
!             write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!             print *,"----------------------------------------------------------------"
!             
!          end do ! options of third order
!          
!          ! printing results for one spin-orbital
!
!          write (*,"(T5,A55,I2,A13,A8)") "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL:",&
!               int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1))," OF SPECIES:",nameOfSpeciesA
!          write (*, "(T5,A45)") "--------------------------------------------"
!          write (*, "(T10,A10,A10,A10)") " Method ","BE (eV)","Pole S."
!          write (*, "(T5,A45)") "--------------------------------------------"
!          write (*,"(T10,A10,F10.3)") "KT        ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T10,A10,F10.3,F10.4)") "EP2       ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)
!          do o=1,6
!
!             write (*,"(T10,A10,F10.3,F10.4)") thirdOrderMethods(o),thirdOrderResults(1,o),thirdOrderResults(2,o)
!             
!          end do
!          write (*, "(T5,A45)") "--------------------------------------------"
!
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,5)=27.211396_8 * newOmega
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,6)=poleStrenght
!                              
!          ! call Matrix_destructor(auxMatrix2(:))          
!!          call TransformIntegrals_destructor( repulsionTransformer )
!          
!       end do
!       
!    end do
!
!    do i = 1, PropagatorTheory_instance%numberOfSpecies
!
!       print *,"Second order densities for species:",i
!       print *,secondOrderDensities(i)%values(:,:)
!
!    end do
!    
!    !!
!    !!************************************************************************************************
!    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
!    print *,"***************************************************************"
!  end subroutine PropagatorTheory_thirdOrderCorrection3
!
!  !**
!  ! @brief Evaluate partial third order, EP3 and OVGF poles, Correction to OVGF for considering a correction
!  ! for each kind of correlation  ! v1
!  !**
!
!  subroutine PropagatorTheory_thirdOrderCorrection4()
!    implicit NONE
!    
!    integer :: ia, ja, ka, la ! Indices for occupied orbitals of alpha (A) species
!    integer :: ib, jb, kb, lb ! Indices for occupied orbitals of beta (B) species
!    integer :: ic, jc, kc, lc ! Indices for occupied orbitals of gamma (C) species
!    integer :: aa, ba, ca, da ! Indices for virtual orbitals of alpha (A) species
!    integer :: ab, bb, cb, db ! Indices for virtual orbitals of beta (B) species
!    integer :: ac, bc, cc, dc ! Indices for virtual orbitals of gamma (C) species
!    integer :: pa, qa, ra, sa ! Indices for general orbitals of alpha (A) species
!    integer :: pb, qb, rb, sb ! Indices for general orbitals of beta (B) species
!    integer :: pc, qc, rc, sc ! Indices for general orbitals of gamma (C) species
!    integer :: idfHf, idaHf ! Counters for elements in fHf and aHf blocks
!    integer :: i, j, k ! counters for species
!    integer :: m, n, o, p, q, r, s, ni, nc, limit, id1, id2 ! auxiliar counters
!    integer :: speciesAID, speciesBID, speciesCID
!    integer :: species1ID, species2ID
!    integer :: electronsID
!    integer :: occupationNumberOfSpeciesA, virtualNumberOfSpeciesA
!    integer :: occupationNumberOfSpeciesB, virtualNumberOfSpeciesB
!    integer :: occupationNumberOfSpeciesC, virtualNumberOfSpeciesC
!    integer :: activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB
!    integer :: activeOrbitalsOfSpeciesC
!    integer(8) :: vectorSize1, vectorSize2 !!! Sizes for diagrams
!    integer(8) :: auxIndex, electrons(2)
!    integer(4) :: errorNum
!    character(10) :: thirdOrderMethods(6)
!    character(10) :: nameOfSpeciesA, nameOfSpeciesB, nameOfSpeciesC
!    type(Vector) :: occupationsOfSpeciesA, occupationsOfSpeciesB, occupationsOfSpeciesC
!    type(Vector) :: eigenValuesOfSpeciesA, eigenValuesOfSpeciesB, eigenValuesOfSpeciesC
!    real(8) :: lambdaOfSpeciesA,  lambdaOfSpeciesB,  lambdaOfSpeciesC 
!    real(8) :: chargeOfSpeciesA, chargeOfSpeciesB, chargeOfSpeciesC
!!    type(TransformIntegrals) :: repulsionTransformer
!    type(Matrix),allocatable :: auxMatrix2(:,:), selfEnergy2hp(:), selfEnergy2ph(:)
!    type(Matrix),allocatable :: secondOrderDensities(:)
!    type(Matrix) :: diagram_A, diagram_B, auxMatrix, auxMatrix3
!    type(Matrix) :: partialMO1, partialMO2    
!    real(8) :: auxVal, auxVal_1, auxVal_2, auxVal_3
!    real(8) :: auxValue_A, auxValue_B, auxValue_C, auxValue_D 
!    real(8) :: auxValue_E, auxValue_F, auxValue_G, auxValue_H
!    real(8) :: valueOfW, valueOfU, valueOfdU, sub2, subW, subU, subd2, subdW, subdU
!    real(8) :: lastOmega, newOmega, residual, threshold, selfEnergy, selfEnergyDerivative, koopmans 
!    real(8) :: a1, a2, b, c, d, poleStrenght, partialValue, partialValue2, initialValue
!    real(8) :: fW, fI, thirdOrderResults(2,5)
!    real(8) :: value1, value2, value3, value4
!    real(8),allocatable :: s2hp(:), s2ph(:), W2hp(:), W2ph(:), U2hp(:), U2ph(:), constantSelfEnergy(:,:), factors(:,:,:)
!    logical :: paso1, paso2, paso3
!    ! *******************************************************************************************
!    ! Determinate the numerators and denominators of the second Oder propapator 
!    
!    if ( .not.CONTROL_instance%OPTIMIZE ) then
!       print *,"===================================================="
!       print *,"      BEGIN FOUR-INDEX INTEGRALS TRANSFORMATION:    "
!       print *,"===================================================="
!       print *,"    Algorithm Four-index integral tranformation"
!       print *,"      Yamamoto, Shigeyoshi; Nagashima, Umpei. "
!       print *,"  Computer Physics Communications, 2005, 166, 58-65 "
!       print *,"--------------------------------------------------"
!       print *,""
!       
!    end if
!    
!    print *,"*******************************************************************"
!    print *,"BEGINNING OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS:"
!
!    !!! Allocating matrix for transformed integrals !!! The algorithm for Integral transformation should be modified
!
!    !!! Determining who are electrons
!
!    if (allocated(auxMatrix2)) deallocate(auxMatrix2)
!    allocate(auxMatrix2(PropagatorTheory_instance%numberOfSpecies,PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(secondOrderDensities)) deallocate(secondOrderDensities)
!    allocate(secondOrderDensities(PropagatorTheory_instance%numberOfSpecies))
!
!    ! Containers for different types of corrections
!
!    if (allocated(s2hp)) deallocate(s2hp)
!    allocate(s2hp(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(s2ph)) deallocate(s2ph)
!    allocate(s2ph(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(W2hp)) deallocate(W2hp)
!    allocate(W2hp(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(W2ph)) deallocate(W2ph)
!    allocate(W2ph(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(U2hp)) deallocate(U2hp)
!    allocate(U2hp(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(U2ph)) deallocate(U2ph)
!    allocate(U2ph(PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(constantSelfEnergy)) deallocate(constantSelfEnergy)
!    allocate(constantSelfEnergy(PropagatorTheory_instance%numberOfSpecies,PropagatorTheory_instance%numberOfSpecies))
!
!    if (allocated(factors)) deallocate(factors)
!    allocate(factors(PropagatorTheory_instance%numberOfSpecies,3,6))
!
!    ! Defining for which species the correction will be applied
!    
!    if (CONTROL_instance%IONIZE_SPECIE /= "NONE") then
!       species1ID = MolecularSystem_getSpecieID( nameOfSpecie=CONTROL_instance%IONIZE_SPECIE )
!       species2ID= species1ID
!       m=1
!    else
!       species1ID=1
!       species2ID=PropagatorTheory_instance%numberOfSpecies
!       m = species2ID
!    end if
!
!    if (allocated(PropagatorTheory_instance%thirdOrderCorrections)) deallocate(PropagatorTheory_instance%thirdOrderCorrections)
!    allocate(PropagatorTheory_instance%thirdOrderCorrections(m))
!
!    ! Storing transformed integrals !!!! We need a more efficient algorithm to do this
!
!    r = 0 ! counter for locating who are electrons
!
!    electrons(1)=0
!    electrons(2)=0
!
!!    call TransformIntegrals_constructor( repulsionTransformer )
!    
!    do p = 1 , PropagatorTheory_instance%numberOfSpecies
!
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( p ) )
!
!       if (nameOfSpeciesA=="e-ALPHA".or.nameOfSpeciesA=="e-BETA") then
!
!          r = r+1
!
!          electrons(r)=p
!
!       end if
!
!       do n = p , PropagatorTheory_instance%numberOfSpecies
!          
!          if (n==p) then
!             
!!             call TransformIntegrals_atomicToMolecularOfOneSpecie( repulsionTransformer,&
!!                  MolecularSystem_getEigenvectors(p), auxMatrix2(p,p), p, trim(nameOfSpeciesA) )
!             
!             auxMatrix2(p,p)%values = auxMatrix2(p,p)%values * MolecularSystem_getCharge( p ) &
!                  * MolecularSystem_getCharge( p )
!             
!          else
!             
!             nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( n ) )
!             
!!             call TransformIntegrals_atomicToMolecularOfTwoSpecies( repulsionTransformer, &
!!                  MolecularSystem_getEigenVectors(p), MolecularSystem_getEigenVectors(n), &
!!                  auxMatrix2(p,n), p, nameOfSpeciesA, n, nameOfSpeciesB )
!             
!             auxMatrix2(p,n)%values = auxMatrix2(p,n)%values * MolecularSystem_getCharge( n ) &
!                  * MolecularSystem_getCharge( p )
!             
!          end if
!          
!       end do
!
!    end do
!
!    ! Start loop for species
!    
!    q = 0
!    
!    do i = species1ID , species2ID
!       
!       q = q + 1
!       
!       nameOfSpeciesA = trim(  MolecularSystem_getNameOfSpecie( i ) )
!       chargeOfSpeciesA = MolecularSystem_getCharge( i )
!!       eigenValuesOfSpeciesA = MolecularSystem_getEigenValues( i )
!       occupationNumberOfSpeciesA = MolecularSystem_getOcupationNumber( i )
!       activeOrbitalsOfSpeciesA = MolecularSystem_getTotalNumberOfContractions( i )
!       lambdaOfSpeciesA = MolecularSystem_getLambda( i )
!       virtualNumberOfSpeciesA = activeOrbitalsOfSpeciesA - occupationNumberOfSpeciesA
!       
!       ! paso
!       
!       paso1=(nameOfSpeciesA=="e-ALPHA".or.nameOfSpeciesA=="e-BETA")
!
!       ! Defining the number of orbitals !!! Insert a parameter for the else option
!       
!       if (CONTROL_instance%PT_JUST_ONE_ORBITAL) then
!          PropagatorTheory_instance%virtualBoundary=CONTROL_instance%IONIZE_MO
!          PropagatorTheory_instance%occupationBoundary=CONTROL_instance%IONIZE_MO
!          n = 1
!       else if (CONTROL_instance%IONIZE_SPECIE /= "NONE".and.CONTROL_instance%IONIZE_MO /= 0) then
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = CONTROL_instance%IONIZE_MO
!          n = PropagatorTheory_instance%virtualBoundary-PropagatorTheory_instance%occupationBoundary+1
!       else
!          PropagatorTheory_instance%virtualBoundary = occupationNumberOfSpeciesA + 1
!          PropagatorTheory_instance%occupationBoundary = occupationNumberOfSpeciesA
!          n = 2
!       end if
!
!       call Matrix_constructor(PropagatorTheory_instance%thirdOrderCorrections(q), int(n,8), 8, 0.0_8)
!
!       !**************************************************************************
!       !	Storing of denominators and numerators in the corresponding vectors
!       !****
!       
!       m =0
!       
!       do pa=PropagatorTheory_instance%occupationBoundary, PropagatorTheory_instance%virtualBoundary	
!
!          m=m+1          
!          
!          ! calculation of constant self energy
!          
!          constantSelfEnergy(:,:) = 0.0_8
!          
!          do p = 1 , PropagatorTheory_instance%numberOfSpecies
!             
!             if (p==i) then
!
!                print *,"entro al if"                
!
!                ! alpha-alpha-alpha
!
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ka = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,i) = constantSelfEnergy(i,i) + partialValue*(auxValue_E-auxValue_F)                       
!
!                   end do
!                end do
!                
!                do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                   do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8                      
!                      
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(ba)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,i) = constantSelfEnergy(i,i) + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(p,p)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(p,p)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ka = 1 , occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do ja = 1 , occupationNumberOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ba, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ia, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,i) = constantSelfEnergy(i,i) &
!                           + (auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma a-a-a:", constantSelfEnergy(i,i)
!
!             else
!
!                print *,"entro al else"                
!
!                nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( p ) )
!                chargeOfSpeciesB = MolecularSystem_getCharge( p )
!!                eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( p )
!                occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( p )
!                activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( p )
!                lambdaOfSpeciesB = MolecularSystem_getLambda( p )
!                virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                
!                ! alpha-beta-beta 
!
!                do ib = 1 , occupationNumberOfSpeciesB
!                   do jb = 1 , occupationNumberOfSpeciesB
!                      
!                      if (p>i) then
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                         auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                         
!                      else
!
!                         auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                         auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                      end if
!
!                      partialValue = 0.0_8
!                      
!                      do kb = 1 , occupationNumberOfSpeciesB
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    - 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
!                                    +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))&
!                                    *( eigenValuesOfSpeciesB%values(jb)&
!                                    + eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(p,p) = constantSelfEnergy(p,p) + partialValue*auxValue_E
!                      
!                   end do
!                end do
!                
!                do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                   do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                      if (p>i) then
!
!                         auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                         auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                      
!                      else
!
!                         auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                         auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                      end if
!                      partialValue = 0.0_8                      
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue &
!                                    + 0.5_8*(auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/(( eigenValuesOfSpeciesB%values(ib)&
!                                    +eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(cb))&
!                                    *( eigenValuesOfSpeciesB%values(ib)&
!                                    + eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesB%values(cb) - eigenValuesOfSpeciesB%values(bb)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(p,p) = constantSelfEnergy(p,p) + partialValue*auxValue_E 
!                      
!                   end do
!                end do
!                
!                do ib = 1 , occupationNumberOfSpeciesB
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      
!                      if (p>i) then
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                         auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!                      
!                      else
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                         auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                      end if
!
!                      partialValue = 0.0_8
!                      
!                      do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do kb = 1 , occupationNumberOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, kb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, kb, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(jb, ab, kb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(jb, bb, kb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue - (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
!                                    +eigenValuesOfSpeciesB%values(kb) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do jb = 1 , occupationNumberOfSpeciesB
!                         do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            do cb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(bb, ab, cb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(bb, jb, cb, ab, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, bb, jb, cb, activeOrbitalsOfSpeciesB )
!                               auxValue_C= auxMatrix2(p,p)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ib, cb, jb, bb, activeOrbitalsOfSpeciesB )
!                               auxValue_D= auxMatrix2(p,p)%values(auxIndex, 1)
!                               
!                               partialValue = partialValue + (auxValue_A-auxValue_B)*(auxValue_C-auxValue_D)/( eigenValuesOfSpeciesB%values(jb)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesB%values(cb))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(p,p) = constantSelfEnergy(p,p) &
!                           + auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma after a-b-b:", constantSelfEnergy(p,p)
!
!                ! alpha-alpha-beta
!
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do ja = 1 , occupationNumberOfSpeciesA
!
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, ja, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                               if (p>i) then
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                               else
!
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                               end if
!
!                               partialValue = partialValue &
!                                    - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ja)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,p) = constantSelfEnergy(i,p) + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                   do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, aa, ba, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!
!                               if (p>i) then                               
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                               else
!
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                               end if
!
!                               partialValue = partialValue &
!                                    + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))&
!                                    *( eigenValuesOfSpeciesA%values(ia)&
!                                    + eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab)))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,p) = constantSelfEnergy(i,p) + partialValue*(auxValue_E-auxValue_F) 
!                      
!                   end do
!                end do
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      
!                      auxIndex = IndexMap_tensorR4ToVector(pa, pa, ia, aa, activeOrbitalsOfSpeciesA )
!                      auxValue_E= auxMatrix2(i,i)%values(auxIndex, 1)
!                      auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, pa, activeOrbitalsOfSpeciesA )
!                      auxValue_F= auxMatrix2(i,i)%values(auxIndex, 1)
!                      
!                      partialValue = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            do ib =  1 , occupationNumberOfSpeciesB
!
!                               if (p>i) then                                                              
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                               else
!
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ja, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                               end if
!
!                               partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ja)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               if (p>i) then                                                              
!
!                                  auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(i,p)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                               else
!
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, aa, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_A= auxMatrix2(p,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_B= auxMatrix2(p,i)%values(auxIndex, 1)
!
!
!                               end if
!
!                               partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesA%values(ia)&
!                                    +eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab))
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      constantSelfEnergy(i,p) = constantSelfEnergy(i,p) &
!                           + 2.0_8*(auxValue_E-auxValue_F)*partialValue/(eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa))
!                      
!                   end do
!                end do
!
!                print *,"constant sigma after a-a-b:", constantSelfEnergy(i,p)
!
!                ! Three species
!
!                ! alpha-beta-alpha
!
!                ! alpha-beta-gamma
!
!                do r = 1, PropagatorTheory_instance%numberOfSpecies
!
!                   if ( r /= p) then
!
!                      print *,"entro a r diferente de p"
!
!                      nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( r ) )
!                      chargeOfSpeciesC = MolecularSystem_getCharge( r )
!!                      eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( r )
!                      occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( r )
!                      activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( r )
!                      lambdaOfSpeciesC = MolecularSystem_getLambda( r )
!                      virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do jb = 1 , occupationNumberOfSpeciesB
!
!                            if (p>i) then                                                              
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, jb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                            else
!
!                               auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                            end if
!
!                            partialValue = 0.0_8
!                            
!                            do ic = 1 , occupationNumberOfSpeciesC
!                               do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!
!                                     if (r>p) then                                                                                                   
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, jb, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
!
!                                     end if
!
!                                     partialValue = partialValue &
!                                          - (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
!                                          *( eigenValuesOfSpeciesB%values(jb)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy(p,r) = constantSelfEnergy(p,r) + partialValue*auxValue_E
!                            
!                         end do
!                      end do
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                            if (p>i) then
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, pa, ab, bb, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                            else
!
!                               auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!
!                            end if
!
!                               partialValue = 0.0_8
!                            
!                            do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                               do ib = 1 , occupationNumberOfSpeciesB
!                                  do ic = 1 , occupationNumberOfSpeciesC
!      
!                                     if (r>p) then                                                                                                                                  
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
!
!                                     end if
!
!                                     partialValue = partialValue &
!                                          + (auxValue_A*auxValue_B)/(( eigenValuesOfSpeciesB%values(ib)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))&
!                                          *( eigenValuesOfSpeciesB%values(ib)&
!                                          + eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac)))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy(p,r) = constantSelfEnergy(p,r) + partialValue*auxValue_E
!                            
!                         end do
!                      end do
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                            if (p>i) then
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, pa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_E= auxMatrix2(i,p)%values(auxIndex, 1)
!
!                            else
!
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, pa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(p,i)%values(auxIndex, 1)
!                               
!                            end if
!
!                            partialValue = 0.0_8
!                            
!                            do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  do ic =  1 , occupationNumberOfSpeciesC
!
!                                     if (r>p) then
!                                                                                
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, jb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(jb, ab, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, jb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, jb, ab, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
!
!                                     end if
!                                     
!                                     partialValue = partialValue - (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(jb)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesC%values(ac))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            do ic = 1 , occupationNumberOfSpeciesC
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!
!                                     if (r>p) then
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ab, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_A= auxMatrix2(p,r)%values(auxIndex, 1)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, bb, ic, ac, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                        auxValue_B= auxMatrix2(p,r)%values(auxIndex, 1)
!
!                                     else
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ab, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_A= auxMatrix2(r,p)%values(auxIndex, 1)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, bb, activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(r,p)%values(auxIndex, 1)
!                                        
!                                     end if
!
!                                     partialValue = partialValue + (auxValue_A*auxValue_B)/( eigenValuesOfSpeciesB%values(ib)&
!                                          +eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesC%values(ac))
!                                     
!                                  end do
!                               end do
!                            end do
!                            
!                            constantSelfEnergy(p,r) = constantSelfEnergy(p,r) &
!                                 + 2.0_8*auxValue_E*partialValue/(eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesB%values(ab))
!                            
!                         end do
!                      end do
!
!                   end if
!
!                   print *,"constant sigma after a-b-a:", constantSelfEnergy(p,r)
!
!                end do
!                
!             end if
!                
!          end do
!
!          ! end of constant self-energy calculation
!          
!          if (allocated(selfEnergy2hp)) deallocate(selfEnergy2hp)
!          allocate(selfEnergy2hp(PropagatorTheory_instance%numberOfSpecies))
!          
!          if (allocated(selfEnergy2ph)) deallocate(selfEnergy2ph)
!          allocate(selfEnergy2ph(PropagatorTheory_instance%numberOfSpecies))
!
!          do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!             
!             if (j==i) then ! Intraspecies factors
!                
!                vectorSize1 = occupationNumberOfSpeciesA * virtualNumberOfSpeciesA * virtualNumberOfSpeciesA
!                vectorSize2 = occupationNumberOfSpeciesA * occupationNumberOfSpeciesA * virtualNumberOfSpeciesA
!                
!                call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                
!                id1 = 0
!                id2 = 0
!                
!                ! factor 2ph
!                
!                do ia = 1 , occupationNumberOfSpeciesA
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         
!                         auxIndex = IndexMap_tensorR4ToVector(pa, aa, ia, ba, activeOrbitalsOfSpeciesA )
!                         auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
!                         auxIndex = IndexMap_tensorR4ToVector(pa, ba, ia, aa, activeOrbitalsOfSpeciesA )
!                         auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
!                         
!                         id1 = id1 + 1
!                         
!                         selfEnergy2ph(j)%values(1,id1) = auxValue_A - auxValue_B
!                         
!                         selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(aa) &
!                              - eigenValuesOfSpeciesA%values(ba)
!                         
!                         valueOfW = 0.0_8
!                         
!                         do ja = 1, occupationNumberOfSpeciesA
!                            do ka = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ia, ka, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ka, ia, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ka, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ka, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                               
!                            end do
!                         end do
!                         
!                         do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1, occupationNumberOfSpeciesA
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ba, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                               auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                               auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, aa, ia, ca, activeOrbitalsOfSpeciesA )
!                               auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ja, ca, ia, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                               
!                               valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca) )
!                               
!                            end do
!                         end do
!                         
!                         selfEnergy2ph(j)%values(3,id1) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!                if (occupationNumberOfSpeciesA > 1) then
!                   
!                   ! factor 2hp
!                   
!                   do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ja = 1 , occupationNumberOfSpeciesA
!                            
!                            id2 = id2 + 1
!                            
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ia, aa, ja, activeOrbitalsOfSpeciesA )
!                            auxValue_A= auxMatrix2(i,j)%values(auxIndex, 1)
!                            auxIndex = IndexMap_tensorR4ToVector(pa, ja, aa, ia, activeOrbitalsOfSpeciesA )
!                            auxValue_B= auxMatrix2(i,j)%values(auxIndex, 1)
!                            
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A - auxValue_B
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesA%values(ja) 
!                               
!                               valueOfW = 0.0_8
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, aa, ca, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ca, aa, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, ca, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, ca, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + 0.5_8*(auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ja) &
!                                          - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca) )
!                                     
!                                  end do
!                               end do
!                               
!                               do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ka = 1, occupationNumberOfSpeciesA
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ja, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW + (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                     auxValue_C= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                     auxValue_D= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ia, aa, ka, activeOrbitalsOfSpeciesA )
!                                     auxValue_E= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ba, ka, aa, ia, activeOrbitalsOfSpeciesA )
!                                     auxValue_F= auxMatrix2(i,j)%values(auxIndex, 1)
!                                     
!                                     valueOfW = valueOfW - (auxValue_C - auxValue_D)*(auxValue_E - auxValue_F)&
!                                          /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesA%values(ka) &
!                                          - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ba) )
!                                     
!                                  end do
!                               end do
!                               
!                               selfEnergy2hp(j)%values(3,id2) = valueOfW
!                               
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! interspecies
!                   
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                   
!                   ! paso
!                   
!                   paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                        (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!
!                   vectorSize1 = occupationNumberOfSpeciesB * virtualNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   vectorSize2 = occupationNumberOfSpeciesB * occupationNumberOfSpeciesA * virtualNumberOfSpeciesB
!                   
!                   call Matrix_constructor(selfEnergy2ph(j), 3, vectorSize1, 0.0_8)
!                   call Matrix_constructor(selfEnergy2hp(j), 3, vectorSize2, 0.0_8)
!                   
!                   id1 = 0
!                   id2 = 0
!                   
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                            if (j>i) then
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, ab, activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                            else
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ib, ab, pa, aa, activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                               
!                            end if
!
!                            id1 = id1 + 1
!                            
!                            selfEnergy2ph(j)%values(1,id1) = auxValue_A
!                            
!                            selfEnergy2ph(j)%values(2,id1) = eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(aa) &
!                                 - eigenValuesOfSpeciesB%values(ab)
!                            
!                            valueOfW = 0.0_8
!                            
!                            do ia = 1, occupationNumberOfSpeciesA
!                               do jb = 1, occupationNumberOfSpeciesB
!
!                                  if (j>i) then                                  
!
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, ib, jb, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                  else
!
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, ia, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(jb, ab, ia, aa, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                  end if
!
!                                  valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!
!                                  if (j>i) then                                                                    
!                                     
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, ab, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib,  bb, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                  else
!
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, ab, pa, ia, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, bb, ia, aa, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                  end if
!
!                                  valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!                            
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do ia = 1 , occupationNumberOfSpeciesA                                  
!
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, aa, activeOrbitalsOfSpeciesA )
!                                  auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ba, ia, activeOrbitalsOfSpeciesA )
!                                  auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
!                                  if (j>i) then
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, ba, ib, ab, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  else
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, ab, ia, ba, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
!                                  end if
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                       - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesB%values(ab) )
!                                  
!                               end do
!                            end do
!                            
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, ab, ib, bb, activeOrbitalsOfSpeciesB )
!                                  auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, bb, ib, ab, activeOrbitalsOfSpeciesB )
!                                  auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                  if (j>i) then
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, aa, jb, bb, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  else
!                                     auxIndex = IndexMap_tensorR4ToVector(jb, bb, pa, aa, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
!                                  end if
!                                  valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                       /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                       - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                                  
!                               end do
!                            end do
!
!                            selfEnergy2ph(j)%values(3,id1) = valueOfW
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!                            
!                            if (j>i) then                                                                    
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, &
!                                    activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                            else
!
!                               auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ia, &
!                                    activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                            end if
!
!                            selfEnergy2hp(j)%values(1,id2) = auxValue_A
!                            
!                            selfEnergy2hp(j)%values(2,id2) = eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ia) &
!                                 - eigenValuesOfSpeciesB%values(ib)
!                            
!                            valueOfW = 0.0_8
!
!                            do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                  if (j>i) then                                                                                                                                           
!                                     auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, bb, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ia, aa, ib, bb, &
!                                          activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                     auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                  else
!
!                                     auxIndex = IndexMap_tensorR4ToVector(ab, bb, pa, aa, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(ib, bb, ia, aa, &
!                                          activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                     auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                  end if
!
!                               valueOfW = valueOfW + (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do jb = 1 , occupationNumberOfSpeciesB
!                               
!                               if (j>i) then                                                                                                                                           
!                                  auxIndex = IndexMap_tensorR4ToVector(pa, aa, ib, jb, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(ia, aa, jb, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                               else
!
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, jb, pa, aa, &
!                                       activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                  auxIndex = IndexMap_tensorR4ToVector(jb, ab, ia, aa, &
!                                       activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                               end if
!
!                               valueOfW = valueOfW - (auxValue_A * auxValue_B)&
!                                    /( eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!
!                               auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ia, activeOrbitalsOfSpeciesA )
!                               auxValue_A = auxMatrix2(i,i)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, ja, aa, activeOrbitalsOfSpeciesA )
!                               auxValue_B = auxMatrix2(i,i)%values(auxIndex, 1)
!                               if (j>i) then
!                                  auxIndex = IndexMap_tensorR4ToVector(ja, aa, ib, ab, &
!                                       activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                  auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                               else
!                                  auxIndex = IndexMap_tensorR4ToVector(ib, ab, ja, aa, &
!                                       activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                  auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
!                               end if
!                               valueOfW = valueOfW - (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                    - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesB%values(ab) )
!                               
!                            end do
!                         end do
!
!                         do jb = 1 , occupationNumberOfSpeciesB
!                            do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               auxIndex = IndexMap_tensorR4ToVector(ab, ib, bb, jb, activeOrbitalsOfSpeciesB )
!                               auxValue_A= auxMatrix2(j,j)%values(auxIndex, 1)
!                               auxIndex = IndexMap_tensorR4ToVector(ab, jb, bb, ib, activeOrbitalsOfSpeciesB )
!                               auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                               if (j>i) then
!                               auxIndex = IndexMap_tensorR4ToVector(pa, ia, jb, bb,&
!                                    activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                               auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                               else
!                               auxIndex = IndexMap_tensorR4ToVector(jb, bb, pa, ia,&
!                                    activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                               auxValue_C = auxMatrix2(j,i)%values(auxIndex, 1)
!                               end if
!                               valueOfW = valueOfW + (auxValue_A - auxValue_B)*(auxValue_C)&
!                                    /( eigenValuesOfSpeciesB%values(ib) + eigenValuesOfSpeciesB%values(jb) &
!                                    - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesB%values(bb) )
!                               
!                            end do
!                         end do
!                         
!                         selfEnergy2hp(j)%values(3,id2) = valueOfW
!                         
!                      end do
!                   end do
!                end do
!                
!             end if
!             
!          end do
!
!          ! Initial guess
!          koopmans = eigenValuesOfSpeciesA%values(pa)
!          newOmega = koopmans
!          lastOmega = 0.0_8
!          
!          ni = 0
!          limit = 50
!          residual = 1.0_8
!
!          ! Calculation of second order pole
!
!          do while ((residual>0.001_8).or.(limit.lt.ni))
!             
!             ni = ni + 1
!             
!             lastOmega = newOmega
!             selfEnergy = lastOmega - koopmans
!             selfEnergyDerivative = 1.0_8
!             
!             do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                
!                if (j==i) then ! Intraspecies term
!                   
!                   id1=0
!                   id2=0
!                   
!                   do ia = 1 , occupationNumberOfSpeciesA
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            
!                            id1 = id1 + 1
!                            
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   if (occupationNumberOfSpeciesA > 1) then
!                      
!                      ! factor 2hp
!                      
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ja = 1 , occupationNumberOfSpeciesA
!                               
!                               id2 = id2 + 1
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               selfEnergy = selfEnergy - 0.5_8*a1*a1/b
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*a1*a1/(b**2.0_8)
!                                                           
!                            end do
!                         end do
!                      end do
!                      
!                   end if
!                   
!                else ! Interspecies term
!
!                   nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                   chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                   eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                   occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                   activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                   lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                   virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!
!                   id1 = 0
!                   id2 = 0
!
!                   ! diagram A
!                   
!                   do ib = 1 , occupationNumberOfSpeciesB
!                      do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                         do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                            
!                            id1 = id1 + 1
!
!                            a1 = selfEnergy2ph(j)%values(1,id1)
!                            b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                            
!                            selfEnergy = selfEnergy - a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                   ! diagram B
!                   
!                   do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do ib = 1 , occupationNumberOfSpeciesB
!                            
!                            id2 = id2 + 1
!
!                            a1 = selfEnergy2hp(j)%values(1,id2)
!                            b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                            
!                            selfEnergy = selfEnergy -a1*a1/b
!                            selfEnergyDerivative = selfEnergyDerivative + a1*a1/(b**2.0_8)
!                            
!                         end do
!                      end do
!                   end do
!                   
!                end if
!                
!             end do
!             
!             newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!             
!             residual = abs(newOmega-lastOmega)
!             
!             print *,"iteration",ni,"newOmega",newOmega,"residual",residual
!
!          end do ! while
!
!          poleStrenght = 1.0_8/selfEnergyDerivative
!
!          ! Storing corrections
!
!          print *,"M y Q:",m,q
!
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)=real(pa,8)
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)=27.211396_8 * koopmans
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3)=27.211396_8 * newOmega
!          PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)=poleStrenght
!
!          print *,"----------------------------------------------------------------"
!          write (*,"(T5,A25,I2,A13,A8)") "Results for spin-orbital:",int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1)),&
!               " of species: ",nameOfSpeciesA
!          write (*,"(T5,A17,F8.4)") "Koopmans' value: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T5,A29,F8.4,A7,I2,A12)") "Optimized second order pole: ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               " after ",ni," iterations."
!          write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!
!          ! Calculation of third order poles
!
!          ! factors for different algorythms
!
!          factors(:,:,:) = 0.0_8
!          thirdOrderResults(:,:) = 0.0_8
!
!          ! o=1 P3
!          thirdOrderMethods(1)="P3"
!          ! o=2 EP3
!          thirdOrderMethods(2)="EP3"
!          ! o=3 OVGF version A
!          thirdOrderMethods(3)="OVGF A"
!          ! o=4 OVGF version B
!          thirdOrderMethods(4)="OVGF B"
!          ! o=5 OVGF version C
!          thirdOrderMethods(5)="OVGF C"
!          ! o=6 
!          thirdOrderMethods(6)="REN-P3"
!          
!          do o = 1 , 6 ! Options for third order and renormalized third order
!
!             ! Initial guess             
!             koopmans = eigenValuesOfSpeciesA%values(pa)
!             newOmega = koopmans
!             lastOmega = 0.0_8
!             
!             ni = 0
!             limit = 15
!             residual = 1.0_8
!
!             ! factor for W
!             
!             fW = 2.0_8
!             fI = 1.0_8
!             if (o==1 .or. o==6) fW=1.0_8
!             if (o==1 .or. o==6) fI=0.0_8
!             threshold=0.00005_8
!             if (o==1 .or. o==2) threshold=0.00001_8
!             
!             ! NR procedure
!
!             do while ((residual>threshold))
!                
!                ni = ni + 1
!                
!                lastOmega = newOmega
!                selfEnergy = lastOmega - koopmans
!                selfEnergyDerivative = 1.0_8
!                s2hp(:) = 0.0_8
!                s2ph(:) = 0.0_8
!                W2hp(:) = 0.0_8
!                W2ph(:) = 0.0_8
!                U2hp(:) = 0.0_8
!                U2ph(:) = 0.0_8
!                
!                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                   
!                   if (j==i) then ! Intraspecies term
!                      
!                      id1=0
!                      id2=0
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8                      
!
!                      ! 2ph terms a-a-a
!
!                      do ia = 1 , occupationNumberOfSpeciesA
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                               
!                               id1 = id1 + 1
!
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!
!                               if ( (.not.paso1).or.(o/=1)) then
!                                  
!                                  do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do da = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ia, da, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, da, ia, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, aa, da, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ca, ba, da, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ia) &
!                                             - eigenValuesOfSpeciesA%values(ca) - eigenValuesOfSpeciesA%values(da)
!                                        
!                                        valueOfU = valueOfU + 0.5_8*a2/c
!                                        valueOfdU = valueOfdU - 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                     do ca = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ba, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, aa, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(ba) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueOfdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, aa, ja, ca, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ca, ja, aa, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ca, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ba, ca, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ja) &
!                                             - eigenValuesOfSpeciesA%values(aa) - eigenValuesOfSpeciesA%values(ca)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do                                  
!
!                               end if
!
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                                        
!                            end do
!                         end do
!                      end do
!
!                      if (paso1.and.(o==1.or.o==6)) subW=0.0_8
!                      if (paso1.and.(o==1.or.o==6)) subdW=0.0_8                      
!
!                      s2ph(j) = s2ph(j) + 0.5_8*sub2
!                      W2ph(j) = W2ph(j) + 0.5_8*(fW*subW)/(1.0_8-factors(i,1,o))                       
!                      U2ph(j) = U2ph(j) + 0.5_8*(subU)/(1.0_8-factors(i,1,o))                       
!
!                      selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(i,1,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(i,1,o)) )                     
!                      ! Diagram 2hp a-a-a
!                        
!                      if (occupationNumberOfSpeciesA > 1) then
!
!                         sub2 = 0.0_8
!                         subW = 0.0_8 
!                         subU = 0.0_8 
!                         subd2 = 0.0_8
!                         subdW = 0.0_8
!                         subdU = 0.0_8
!                         
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ia = 1 , occupationNumberOfSpeciesA
!                               do ja = 1 , occupationNumberOfSpeciesA
!                                  
!                                  id2 = id2 + 1
!                                  
!                                  valueOfU = 0.0_8
!                                  valueOfdU = 0.0_8
!
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do la = 1 , occupationNumberOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, aa, la, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, la, aa, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ia, la, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ka, ja, la, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(aa) &
!                                             - eigenValuesOfSpeciesA%values(ka) - eigenValuesOfSpeciesA%values(la)
!                                        
!                                        valueOfU = valueOfU - 0.5_8*a2/c
!                                        valueOfdU = valueOfdU + 0.5_8*a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ka = 1 , occupationNumberOfSpeciesA
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ia, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ja) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueOfdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, ba, ka, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ka, ba, ia, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ba, ka, ja, activeOrbitalsOfSpeciesA )
!                                        auxValue_C = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(aa, ja, ka, ba, activeOrbitalsOfSpeciesA )
!                                        auxValue_D = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A - auxValue_B)*(auxValue_C - auxValue_D)
!                                        c = lastOmega + eigenValuesOfSpeciesA%values(ba) &
!                                             - eigenValuesOfSpeciesA%values(ia) - eigenValuesOfSpeciesA%values(ka)
!                                        
!                                        valueOfU = valueOfU + a2/c
!                                        valueofdU = valueOfdU - a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                                                                                   
!                                  a1 = selfEnergy2hp(j)%values(1,id2)
!                                  a2 = selfEnergy2hp(j)%values(3,id2)
!                                  b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                                  
!                                  sub2 = sub2 + (a1**2.0_8)/b
!                                  subW = subW + (a1*a2)/b
!                                  subU = subU + (a1*valueOfU)/b
!                                  
!                                  subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                                  subdW = subdW + (a1*a2)/(b**2.0_8)
!                                  subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                  
!                               end do
!                            end do
!                         end do
!                         
!                         s2hp(i) = s2hp(i) + 0.5_8*sub2
!                         W2hp(i) = W2hp(i) + 0.5_8*(fW*subW)/(1.0_8-factors(i,2,o))                       
!                         U2hp(i) = U2hp(i) + 0.5_8*(subU)/(1.0_8-factors(i,2,o))                       
!                         
!                         selfEnergy = selfEnergy - 0.5_8*( sub2+ (fW*subW+subU)/(1.0_8-factors(i,2,o)) )                       
!                         
!                         selfEnergyDerivative = selfEnergyDerivative + 0.5_8*( subd2+ (fW*subdW+subdU)/(1.0_8-factors(i,2,o)) )                                                                                           
!                      end if
!
!                      ! terms a-a-b
!                         
!                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                         
!                         if (k .ne. i)  then
!                            
!                            id1=0
!                            id2=0
!
!                            nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                            chargeOfSpeciesB = MolecularSystem_getCharge( k )
!!                            eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( k )
!                            occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( k )
!                            activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( k )
!                            lambdaOfSpeciesB = MolecularSystem_getLambda( k )
!                            virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                            
!                            ! 2ph a-a-b
!
!                            if ( (.not.paso1).or.(o/=1.and.o/=6)) then                            
!
!                               print *,"ENTRA A LA PARTE 2PH A-A-B CON A:",i,"Y B:",k
!
!                               subU = 0.0_8 
!                               subdU = 0.0_8                      
!                               subW = 0.0_8
!                               subdW = 0.0_8
!                            
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        id1 = id1 + 1
!
!                                        valueOfU = 0.0_8
!                                        valueOfdU = 0.0_8
!                                        valueOfW = 0.0_8
!                                                                                
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                              if (k>i) then
!                                                 
!                                                 auxIndex = IndexMap_tensorR4ToVector(pa, ba, ab, ib,&
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib,&
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!
!                                              else
!
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ba,&
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, aa, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
!
!                                              end if
!
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueOfW = valueOfW - a2/d
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                              if (k>i) then
!                                                 
!                                                 auxIndex = IndexMap_tensorR4ToVector(pa, aa, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ia, ba, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!
!                                              else
!
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, aa, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, ba, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
!
!                                              end if
!                                                 
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueOfW = valueOfW + a2/d
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                        a1 = selfEnergy2ph(j)%values(1,id1)
!                                        a2 = selfEnergy2ph(j)%values(3,id1)
!                                        b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                                        
!                                        subU = subU + (a1*valueOfU)/b                                        
!                                        subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                        
!                                        subW = subW + (a1*valueOfW)/b
!                                        subdW = subdW + (a1*valueOfW)/(b**2.0_8)
!                                        
!                                     end do
!                                  end do
!                               end do
!                               
!                               U2ph(k) = U2ph(k) + 0.5_8*(subU)/(1.0_8-factors(k,1,o))                    
!                               W2ph(k) = W2ph(k) + 0.5_8*(fW*subW)/(1.0_8-factors(k,1,o))                          
!
!                               selfEnergy = selfEnergy - 0.5_8*(subU+fW*subW)/(1.0_8-factors(k,1,o))                       
!                               
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*(subdU+fW*subdW)/(1.0_8-factors(k,1,o))                     
!
!                            end if
!                            
!                            ! 2hp a-a-b
!
!                            if (occupationNumberOfSpeciesA > 1) then
!
!                               subU = 0.0_8 
!                               subdU = 0.0_8     
!                               subW = 0.0_8
!                               subdW = 0.0_8
!                               
!                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ia = 1 , occupationNumberOfSpeciesA
!                                     do ja = 1 , occupationNumberOfSpeciesA
!                                        
!                                        id2 = id2 + 1
!                                        
!                                        valueOfU = 0.0_8
!                                        valueOfdU = 0.0_8
!                                        valueOfW = 0.0_8
!                                        
!                                        do ib = 1 , occupationNumberOfSpeciesB
!                                           do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                              if (k>i) then
!                                                 
!                                                 auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ia, aa, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 
!                                              else
!
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ja, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, ia, aa, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
!
!                                              end if
!                                                    
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                              d = eigenValuesOfSpeciesA%values(ia) + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                              
!                                              valueOfU = valueOfU + 2.0_8*a2/c
!                                              valueOfW = valueOfW - a2/d
!                                              valueofdU = valueOfdU - 2.0_8*a2/(c**2.0_8)
!
!                                              if (k>i) then                                              
!                                                 
!                                                 auxIndex = IndexMap_tensorR4ToVector(pa, ia, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ja, aa, ab, ib, &
!                                                      activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                                 auxValue_B = auxMatrix2(i,k)%values(auxIndex, 1)
!                                                 
!                                              else
!
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, pa, ia, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                                 auxIndex = IndexMap_tensorR4ToVector(ab, ib, ja, aa, &
!                                                      activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                                 auxValue_B = auxMatrix2(k,i)%values(auxIndex, 1)
!
!                                              end if
!
!                                              a2 = auxValue_A*auxValue_B
!                                              c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                                   - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ia)
!                                              d = eigenValuesOfSpeciesA%values(ja) + eigenValuesOfSpeciesB%values(ib) &
!                                                   - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(aa)
!                                              
!                                              valueOfU = valueOfU - 2.0_8*a2/c
!                                              valueOfW = valueOfW + a2/d
!                                              valueofdU = valueOfdU + 2.0_8*a2/(c**2.0_8)
!                                              
!                                           end do
!                                        end do
!                                        
!                                        a1 = selfEnergy2hp(j)%values(1,id2)
!                                        a2 = selfEnergy2hp(j)%values(3,id2)
!                                        b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!
!                                        subU = subU + (a1*valueOfU)/b
!                                        subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!
!                                        subW = subW + (a1*valueOfW)/b
!                                        subdW = subdW + (a1*valueOfW)/(b**2.0_8)
!                                        
!                                     end do
!                                  end do
!                               end do
!                               
!                               U2hp(k) = U2hp(k) + 0.5_8*(subU)/(1.0_8-factors(k,2,o))                       
!                               W2hp(k) = W2hp(k) + 0.5_8*(fW*subW)/(1.0_8-factors(k,2,o))                          
!                               
!                               selfEnergy = selfEnergy - 0.5_8*(subU+fW*subW)/(1.0_8-factors(k,2,o))                       
!                               
!                               selfEnergyDerivative = selfEnergyDerivative + 0.5_8*(subdU+fW*subdW)/(1.0_8-factors(k,2,o))                     
!                               
!                            end if
!
!                         end if
!                            
!                      end do
!                     
!                      selfEnergy = selfEnergy - fI*constantSelfEnergy(i,i)/(1.0_8-factors(i,3,o))                
!                      
!                   else ! Interspecies term
!                      
!                      nameOfSpeciesB = trim(  MolecularSystem_getNameOfSpecie( j ) )
!                      chargeOfSpeciesB = MolecularSystem_getCharge( j )
!!                      eigenValuesOfSpeciesB = MolecularSystem_getEigenValues( j )
!                      occupationNumberOfSpeciesB = MolecularSystem_getOcupationNumber( j )
!                      activeOrbitalsOfSpeciesB = MolecularSystem_getTotalNumberOfContractions( j )
!                      lambdaOfSpeciesB = MolecularSystem_getLambda( j )
!                      virtualNumberOfSpeciesB = activeOrbitalsOfSpeciesB - occupationNumberOfSpeciesB
!                      
!                      paso2=(nameOfSpeciesA=="e-ALPHA".and.nameOfSpeciesB=="e-BETA").or.&
!                           (nameOfSpeciesA=="e-BETA".and.nameOfSpeciesB=="e-ALPHA")
!                      
!                      id1 = 0
!                      id2 = 0
!                      
!                      ! Diagram 2ph interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ib = 1 , occupationNumberOfSpeciesB
!                         do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               
!                               id1 = id1 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               if ( (.not.paso2).or.(o/=1.and.o/=6) ) then
!                                  
!                                  do jb = 1 , occupationNumberOfSpeciesB
!                                     
!                                     do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                        if (j>i) then
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, aa, bb, jb, &
!                                                activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                        else
!
!                                           auxIndex = IndexMap_tensorR4ToVector(bb, jb, pa, aa, &
!                                                activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                        end if
!
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, jb, bb, ab, activeOrbitalsOfSpeciesB )
!                                        auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, ab, bb, jb, activeOrbitalsOfSpeciesB )
!                                        auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
!                                        
!                                        a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(aa)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                        
!                                        if (j>i) then
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(ba, aa, bb, ab, &
!                                                activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ba, ib, bb, &
!                                                activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                        else
!
!                                           auxIndex = IndexMap_tensorR4ToVector(bb, ab, ba, aa, &
!                                                activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, bb, pa, ba, &
!                                                activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                           auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                        end if
!
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(ib) &
!                                             - eigenValuesOfSpeciesB%values(bb) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU + (a2/c)
!                                        valueofdU = valueOfdU - (a2/(c**2.0_8))
!                                        
!                                     end do
!                                  end do
!                                  
!                                  do ba = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                     do jb = 1 , occupationNumberOfSpeciesB
!
!                                        if (j>i) then                                        
!                                           
!                                           auxIndex = IndexMap_tensorR4ToVector(aa, ba, ib, jb, &
!                                                activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(pa, ba, jb, ab, &
!                                                activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                           auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                        else
!
!                                           auxIndex = IndexMap_tensorR4ToVector(ib, jb, aa, ba, &
!                                                activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                           auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                           auxIndex = IndexMap_tensorR4ToVector(jb, ab, pa, ba, &
!                                                activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                           auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                        end if
!
!                                        a2 = auxValue_A*auxValue_B
!                                        c = lastOmega + eigenValuesOfSpeciesB%values(jb) &
!                                             - eigenValuesOfSpeciesB%values(ab) - eigenValuesOfSpeciesA%values(ba)
!                                        
!                                        valueOfU = valueOfU - a2/c
!                                        valueofdU = valueOfdU + a2/(c**2.0_8)
!                                        
!                                     end do
!                                  end do
!                                  
!                               end if
!                               
!                               a1 = selfEnergy2ph(j)%values(1,id1)
!                               a2 = selfEnergy2ph(j)%values(3,id1)
!                               b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      if (paso2.and.(o==1.or.o==6)) subW=0.0_8
!                      if (paso2.and.(o==1.or.o==6)) subdW=0.0_8                      
!
!                      s2ph(j) = s2ph(j) + sub2
!                      W2ph(j) = W2ph(j) + (fW*subW)/(1.0_8-factors(j,1,o))                       
!                      U2ph(j) = U2ph(j) + (subU)/(1.0_8-factors(j,1,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(j,1,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(j,1,o)) )                                      
!                      ! Diagram 2hp interspecies
!
!                      sub2 = 0.0_8
!                      subW = 0.0_8 
!                      subU = 0.0_8 
!                      subd2 = 0.0_8
!                      subdW = 0.0_8
!                      subdU = 0.0_8
!                      
!                      do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                         do ia = 1 , occupationNumberOfSpeciesA
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               
!                               id2 = id2 + 1
!                               
!                               valueOfU = 0.0_8
!                               valueOfdU = 0.0_8
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  
!                                  do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!
!                                     if (j>i) then
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ia, bb, jb,&
!                                             activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(bb, jb, pa, ia, &
!                                             activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                        
!                                     end if
!
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, ab, ib, jb, activeOrbitalsOfSpeciesB )
!                                     auxValue_B= auxMatrix2(j,j)%values(auxIndex, 1)
!                                     auxIndex = IndexMap_tensorR4ToVector(bb, jb, ib, ab, activeOrbitalsOfSpeciesB )
!                                     auxValue_C= auxMatrix2(j,j)%values(auxIndex, 1)
!                                     
!                                     a2 = (auxValue_A)*(auxValue_B - auxValue_C)
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ia)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do jb = 1 , occupationNumberOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!                                      
!                                     if (j>i) then
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, ib, jb,&
!                                             activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, ab, jb,&
!                                             activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(ib, jb, ia, ja, &
!                                             activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(ab, jb, pa, ja,&
!                                             activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                     end if
!
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(ab) &
!                                          - eigenValuesOfSpeciesB%values(jb) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU - a2/c
!                                     valueofdU = valueOfdU + a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               do bb = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                  do ja = 1 , occupationNumberOfSpeciesA
!
!                                     if (j>i) then                                     
!                                        
!                                        auxIndex = IndexMap_tensorR4ToVector(ia, ja, bb, ab,&
!                                             activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_A = auxMatrix2(i,j)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(pa, ja, bb, ib,&
!                                             activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesB )
!                                        auxValue_B = auxMatrix2(i,j)%values(auxIndex, 1)
!
!                                     else
!
!                                        auxIndex = IndexMap_tensorR4ToVector(bb, ab, ia, ja,&
!                                             activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                        auxValue_A = auxMatrix2(j,i)%values(auxIndex, 1)
!                                        auxIndex = IndexMap_tensorR4ToVector(bb, ib, pa, ja,&
!                                             activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesA )
!                                        auxValue_B = auxMatrix2(j,i)%values(auxIndex, 1)
!
!                                     end if
!
!                                     a2 = auxValue_A*auxValue_B
!                                     c = lastOmega + eigenValuesOfSpeciesB%values(bb) &
!                                          - eigenValuesOfSpeciesB%values(ib) - eigenValuesOfSpeciesA%values(ja)
!                                     
!                                     valueOfU = valueOfU + a2/c
!                                     valueofdU = valueOfdU - a2/(c**2.0_8)
!                                     
!                                  end do
!                               end do
!                               
!                               a1 = selfEnergy2hp(j)%values(1,id2)
!                               a2 = selfEnergy2hp(j)%values(3,id2)
!                               b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!                               
!                               sub2 = sub2 + (a1**2.0_8)/b
!                               subW = subW + (a1*a2)/b
!                               subU = subU + (a1*valueOfU)/b
!                               
!                               subd2 = subd2 + (a1**2.0_8)/(b**2.0_8)
!                               subdW = subdW + (a1*a2)/(b**2.0_8)
!                               subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                               
!                            end do
!                         end do
!                      end do
!                      
!                      s2hp(j) = s2hp(j) + sub2
!                      W2hp(j) = W2hp(j) + (fW*subW)/(1.0_8-factors(j,2,o))                       
!                      U2hp(j) = U2hp(j) + (subU)/(1.0_8-factors(j,2,o))                       
!                      
!                      selfEnergy = selfEnergy - ( sub2+ (fW*subW+subU)/(1.0_8-factors(j,2,o)) )                       
!                      
!                      selfEnergyDerivative = selfEnergyDerivative + ( subd2+ (fW*subdW+subdU)/(1.0_8-factors(j,2,o)) )                                      
!                      do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                         
!                         id1=0
!                         id2=0
!                         
!                         if (k.ne.i .and. k.ne.j)  then
!                            
!                            print *,"ENTRO AL TERMINO DE TRES PARTICULAS:",i,j,k
!                            
!                            nameOfSpeciesC = trim(  MolecularSystem_getNameOfSpecie( k ) )
!                            chargeOfSpeciesC = MolecularSystem_getCharge( k )
!!                            eigenValuesOfSpeciesC = MolecularSystem_getEigenValues( k )
!                            occupationNumberOfSpeciesC = MolecularSystem_getOcupationNumber( k )
!                            activeOrbitalsOfSpeciesC = MolecularSystem_getTotalNumberOfContractions( k )
!                            lambdaOfSpeciesC = MolecularSystem_getLambda( k )
!                            virtualNumberOfSpeciesC = activeOrbitalsOfSpeciesC - occupationNumberOfSpeciesC
!                            
!                            paso3=(nameOfSpeciesB=="e-ALPHA".and.nameOfSpeciesC=="e-BETA").or.&
!                                 (nameOfSpeciesB=="e-BETA".and.nameOfSpeciesC=="e-ALPHA")
!                        
!                            subU=0.0_8
!                            subdU=0.0_8                            
!                            
!                            do ib = 1 , occupationNumberOfSpeciesB
!                               do aa = occupationNumberOfSpeciesA+1 , activeOrbitalsOfSpeciesA
!                                  do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                                     
!                                     id1 = id1 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           if (k>i) then
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, aa, ic, ac,&
!                                                   activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)                                              
!                                           else
!                                              auxIndex = IndexMap_tensorR4ToVector(ic, ac, pa, aa,&
!                                                   activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesA )
!                                              auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                           end if
!                                           if (k>j) then
!                                              auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac,&
!                                                   activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
!                                           else
!                                              auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab,&
!                                                   activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k,j)%values(auxIndex, 1)
!                                           end if
!                                              
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ic) &
!                                                - eigenValuesOfSpeciesC%values(ac) - eigenValuesOfSpeciesA%values(aa)
!                                           
!                                           
!                                           valueOfU = valueOfU - a2/c
!                                           valueofdU = valueOfdU + a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2ph(j)%values(1,id1)
!                                     b = selfEnergy2ph(j)%values(2,id1) + lastOmega
!                                     
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!
!                            if (paso3) then
!                               
!                               U2ph(j)=U2ph(j)+subU
!                               selfEnergy = selfEnergy - (subU)/(1.0_8-factors(j,1,o))
!                               selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(j,1,o))
!
!                            else
!
!                               selfEnergy = selfEnergy - (subU)
!                               selfEnergyDerivative = selfEnergyDerivative + subdU
!
!                            end if
!
!                            subU=0.0_8
!                            subdU=0.0_8
!                            
!                            do ab = occupationNumberOfSpeciesB+1 , activeOrbitalsOfSpeciesB
!                               do ia = 1 , occupationNumberOfSpeciesA
!                                  do ib = 1 , occupationNumberOfSpeciesB
!                                     
!                                     id2 = id2 + 1
!                                     
!                                     valueOfU=0.0_8
!                                     valueOfdU=0.0_8
!                                     
!                                     do ic = 1 , occupationNumberOfSpeciesC
!                                        do ac = occupationNumberOfSpeciesC+1 , activeOrbitalsOfSpeciesC
!                                           
!                                           if (k>i) then
!                                              auxIndex = IndexMap_tensorR4ToVector(pa, ia, ic, ac,&
!                                                   activeOrbitalsOfSpeciesA, activeOrbitalsOfSpeciesC )
!                                              auxValue_A = auxMatrix2(i,k)%values(auxIndex, 1)
!                                           else
!                                              auxIndex = IndexMap_tensorR4ToVector(ic, ac, pa, ia,&
!                                                   activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesA )
!                                              auxValue_A = auxMatrix2(k,i)%values(auxIndex, 1)
!                                           end if
!                                           if (k>j) then
!                                              auxIndex = IndexMap_tensorR4ToVector(ib, ab, ic, ac,&
!                                                   activeOrbitalsOfSpeciesB, activeOrbitalsOfSpeciesC )
!                                              auxValue_B = auxMatrix2(j,k)%values(auxIndex, 1)
!                                           else
!                                              auxIndex = IndexMap_tensorR4ToVector(ic, ac, ib, ab,&
!                                                   activeOrbitalsOfSpeciesC, activeOrbitalsOfSpeciesB )
!                                              auxValue_B = auxMatrix2(k,j)%values(auxIndex, 1)
!                                           end if
!
!                                           a2 = auxValue_A*auxValue_B
!                                           c = lastOmega + eigenValuesOfSpeciesC%values(ac) &
!                                                - eigenValuesOfSpeciesC%values(ic) - eigenValuesOfSpeciesA%values(ia)
!                                           
!                                           valueOfU = valueOfU + a2/c
!                                           valueofdU = valueOfdU - a2/(c**2.0_8)
!                                           
!                                        end do
!                                     end do
!                                     
!                                     a1 = selfEnergy2hp(j)%values(1,id2)
!                                     b = selfEnergy2hp(j)%values(2,id2) + lastOmega
!
!                                     subU = subU + (a1*valueOfU)/b
!                                     subdU = subdU + a1*(valueOfU/(b**2.0_8) - valueOfdU/b)
!                                     
!                                  end do
!                               end do
!                            end do
!
!                            if (paso3) then
!
!                               U2hp(j)=U2hp(j)+subU
!                               selfEnergy = selfEnergy - (subU)/(1.0_8-factors(j,2,o))
!                               selfEnergyDerivative = selfEnergyDerivative + subdU/(1.0_8-factors(j,2,o))
!                               selfEnergy = selfEnergy - fI*constantSelfEnergy(j,k)/(1.0_8-factors(j,3,o))
!
!                            else
!
!                               selfEnergy = selfEnergy - (subU)
!                               selfEnergyDerivative = selfEnergyDerivative + subdU
!                               selfEnergy = selfEnergy - fI*constantSelfEnergy(j,k)
!
!                            end if
!                            
!                         end if
!                         
!                      end do
!
!                      selfEnergy = selfEnergy - fI*constantSelfEnergy(j,j)/(1.0_8-factors(j,3,o))                
!                      selfEnergy = selfEnergy - fI*constantSelfEnergy(i,j)/(1.0_8-factors(j,3,o))                
!                      selfEnergy = selfEnergy - fI*constantSelfEnergy(j,i)/(1.0_8-factors(j,3,o))                
!                                            
!                   end if
!                   
!                end do
!
!                newOmega = lastOmega - (selfEnergy/selfEnergyDerivative)
!                
!                residual = abs(newOmega-lastOmega)
!                
!                print *,"ni",ni,"newOmega",newOmega,"residual",residual
!                
!             end do ! while
!
!             if (o==1) then
!
!                do j = 1 , PropagatorTheory_instance%numberOfSpecies             
!                   
!                   if (j/=electrons(1).and.j/=electrons(2)) then
!
!                      value1=s2ph(j)
!                      value2=W2ph(j)
!                      value3=s2hp(j)
!                      value4=W2hp(j)
!                   
!                      ! Renormalized P3
!                      if (value1/=0.0_8) factors(j,1,6) = value2/value1
!                      if (value3/=0.0_8) factors(j,2,6) = value4/value3
!                      factors(j,3,6) = 0.0_8
!                      
!                   end if
!                   
!                end do
!
!                r = electrons(1)
!                s = electrons(2)
!                print *,"values of r and s:",r,s
!
!                if (r/=0.and.s/=0) then
!
!                   value1=s2ph(r)+s2ph(s)
!                   value2=W2ph(r)+W2ph(s)
!                   value3=s2hp(r)+s2hp(s)
!                   value4=W2hp(r)+W2hp(s)
!
!                   ! Renormalized P3
!
!                   if (value1/=0.0_8) factors(r,1,6)=value2/value1 
!                   if (value3/=0.0_8) factors(r,2,6)=value4/value3 
!                   factors(r,3,6) = 0.0_8
!                   factors(s,1,6) = factors(r,1,6)
!                   factors(s,2,6) = factors(r,2,6)
!                   factors(s,3,6) = factors(r,3,6)
!                   
!                end if
!
!             end if
!             
!             if (o==2) then
!                
!                do k = 1 , PropagatorTheory_instance%numberOfSpecies             
!                   
!                      if (k/=electrons(1).and.k/=electrons(2)) then
!
!                         ! OVGF version A
!
!                         value2=W2hp(k)+W2ph(k)
!                         value1=s2hp(k)+s2ph(k)
!
!                         if (value1/=0.0_8) factors(k,1,3) = value2/value1
!                         factors(k,2,3) = factors(k,1,3)
!                         factors(k,3,3) = factors(k,1,3)
!                         
!                         ! OVGF version B
!
!                         if (s2ph(k)/=0.0_8) factors(k,1,4) = W2ph(k)/s2ph(k)
!                         if (s2hp(k)/=0.0_8) factors(k,2,4) = W2hp(k)/s2hp(k)
!                         factors(k,3,4) = 0.0_8
!                         
!                         ! OVGF version C
!
!                         value1=W2hp(k)+W2ph(k)+U2hp(k)+U2ph(k)
!                         value2=factors(k,2,4)*(U2hp(k)+W2hp(k))+factors(k,1,4)*(U2ph(k)+W2ph(k))
!
!                         if (value1/=0.0_8) factors(k,1,5) = value2/value1
!                         factors(k,2,5) = factors(k,1,5)
!                         factors(k,3,5) = factors(k,1,5)
!
!                      end if
!
!                   end do
!
!                   r = electrons(1)
!                   s = electrons(2)
!                   print *,"values of r and s:",r,s
!
!                   if (r/=0.and.s/=0) then
!
!                      ! OVGF version A
!                      value1=s2hp(r)+s2ph(r)+s2hp(s)+s2ph(s)
!                      value2=W2hp(r)+W2ph(r)+W2hp(s)+W2ph(s)
!
!                      if (value1/=0.0_8) factors(r,1,3) = value2/value1
!                      factors(r,2,3) = factors(r,1,3)
!                      factors(r,3,3) = factors(r,1,3)
!                      factors(s,1,3) = factors(r,1,3)
!                      factors(s,2,3) = factors(r,1,3)
!                      factors(s,3,3) = factors(r,1,3)
!                      
!                      ! OVGF version B
!
!                      value1=s2ph(r)+s2ph(s)
!                      value2=W2ph(r)+W2ph(s)
!                      value3=s2hp(r)+s2hp(s)
!                      value4=W2hp(r)+W2hp(s)
!
!                      if (value1/=0.0_8) factors(r,1,4)=value2/value1 
!                      if (value3/=0.0_8) factors(r,2,4) =value4/value3 
!                      factors(r,3,4) = 0.0_8
!                      factors(s,1,4) = factors(r,1,4)
!                      factors(s,2,4) = factors(r,2,4)
!                      factors(s,3,4) = factors(r,3,4)
!                      
!                      ! OVGF version C
!                      value2=factors(r,2,4)*(U2hp(r)+W2hp(r)+U2hp(s)+W2hp(s))+&
!                           factors(r,1,4)*(U2ph(r)+W2ph(r)+U2ph(s)+W2ph(s))
!                      value1=W2hp(r)+W2ph(r)+U2hp(r)+U2ph(r)+W2hp(s)+W2ph(s)+U2hp(s)+U2ph(s)
!
!                      if (value1/=0.0_8) factors(r,1,5) = value2/value1
!                      factors(r,2,5) = factors(r,1,5)
!                      factors(r,3,5) = factors(r,1,5)
!                      factors(s,1,5) = factors(r,1,5)
!                      factors(s,2,5) = factors(r,2,5)
!                      factors(s,3,5) = factors(r,3,5)
!                      
!                   end if
!
!             end if
!             
!             poleStrenght = 1.0_8/(selfEnergyDerivative)
!             thirdOrderResults(1,o) = 27.211396_8 * newOmega
!             thirdOrderResults(2,o) = poleStrenght
!
!             print *,"value of o:",o
!             ! print *,"FACTORS:"
!             ! print *,factors(:,:,:)
!             print *,"constant self-energy:", constantSelfEnergy*27.211396_8
!             print *,"2hp(2):",s2hp(:)
!             print *,"2ph(2):",s2ph(:)
!             print *,"W 2hp(3):",W2hp(:)
!             print *,"W 2ph(3):",W2ph(:)
!             print *,"U 2hp(3):",U2hp(:)
!             print *,"U 2ph(3):",U2ph(:)
!             print *,"factor 2hp:",factors(:,2,o)
!             print *,"factor 2ph:",factors(:,1,o)
!             write (*,"(T5,A10,A10,A6,F8.4,A7,I2,A12)") "Optimized ",thirdOrderMethods(o),"pole: ",newOmega*27.211396_8," after ",ni," iterations."
!             write (*,"(T5,A11,F8.4,A15,F7.4)") "Correction:",(newOmega-koopmans)*27.211396_8," Pole strength:",poleStrenght
!             print *,"----------------------------------------------------------------"
!                          
!          end do ! options for third order
!
!          ! printing results for one spin-orbital
!
!          write (*,"(T5,A55,I2,A13,A8)") "SUMMARY OF PROPAGATOR RESULTS FOR THE SPIN-ORBITAL:",&
!               int(PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,1))," OF SPECIES:",nameOfSpeciesA 
!          write (*, "(T5,A45)") "--------------------------------------------" 
!          write (*, "(T10,A10,A10,A10)") " Method ","BE (eV)","Pole S."
!          write (*, "(T5,A45)") "--------------------------------------------"
!          write (*,"(T10,A10,F10.3)") "KT        ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,2)
!          write (*,"(T10,A10,F10.3,F10.4)") "EP2       ",PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,3),&
!               PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,4)
!          do o=1,6
!
!             write (*,"(T10,A10,F10.3,F10.4)") thirdOrderMethods(o),thirdOrderResults(1,o),thirdOrderResults(2,o)
!             
!          end do
!          write (*, "(T5,A45)") "--------------------------------------------"
!
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,5)=27.211396_8 * newOmega
!          ! PropagatorTheory_instance%thirdOrderCorrections(q)%values(m,6)=poleStrenght
!                              
!          ! call Matrix_destructor(auxMatrix2(:))          
!!          call TransformIntegrals_destructor( repulsionTransformer )
!          
!       end do
!       
!    end do
!
!    do i = 1, PropagatorTheory_instance%numberOfSpecies
!
!       print *,"Second order densities for species:",i
!       print *,secondOrderDensities(i)%values(:,:)
!
!    end do
!    
!    !!
!    !!************************************************************************************************
!    print *,"END OF GENERALIZED ANY-PARTICLE PROPAGATOR CALCULATIONS"
!    print *,"***************************************************************"
!  end subroutine PropagatorTheory_thirdOrderCorrection4
