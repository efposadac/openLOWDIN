!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	  UNIVERSIDAD NACIONAL DE COLOMBIA"
!!	  PROF. ANDRES REYES GROUP"
!!	  http://www.qcc.unal.edu.co"
!!	
!!	  UNIVERSIDAD DE GUADALAJARA"
!!	  PROF. ROBERTO FLORES GROUP"
!!	  http://www.cucei.udg.mx/~robertof"
!!	
!!	AUTHORS
!!		E.F. POSADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		S.A. GONZALEZ. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		F.S. MONCADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		J. ROMERO. UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!	CONTRIBUTORS
!!		N.F.AGUIRRE. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		GABRIEL MERINO. UNIVERSIDAD DE GUANAJUATO
!!   		J.A. CHARRY UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************
                
module CImod_
  use Exception_
  use Matrix_
  use Vector_
  use MolecularSystem_
  use Configuration_
  use ReadTransformedIntegrals_
  use String_
  use IndexMap_
  use InputCI_
  use omp_lib
  use JadamiluInterface_
  use CIcore_
  use CIDiag_
  use CIFullMatrix_
  use CIInitial_
  use CISCI_
  use CIJadamilu_
  use CIOrder_
  use CIStrings_
  use CIMCSCF_

  ! use ArpackInterface_
  implicit none
      
  !>
  !! @brief Configuration Interaction Module, works in spin orbitals
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 07-24-12
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 07-24-12 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# description.
  !!   - <tt> 07-09-16 </tt>: Jorge Charry ( jacharrym@unal.edu.co )
  !!        -# Add CIS, and Fix CISD.
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<

  public :: &
!       CIcore_constructor, &
       CImod_destructor, &
       CImod_getTotalEnergy, &
       CImod_run, &
       CImod_showEigenVectors, &
       CImod_show

  private

contains

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CImod_run()
    implicit none 
    type(MCSCF) :: MCSCF_instance
    integer :: i, k, numberOfSpecies
    integer :: a, ms
    real(8) :: timeA, timeB
    real(8) :: ecorr

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    ms = CONTROL_instance%CI_MADSPACE

    !! printing header
    write (*,*) ""
    write (*,*) "         BEGIN ", trim(CIcore_instance%level)," CALCULATION"
    write (*,*) "         J. Charry, F. Moncada                "
    write (*,*) "-----------------------------------------------"
    write (*,*) ""

    !! printing active space
    write (*,"(A32)",advance="no") "Number of orbitals for species: "
    do i = 1, numberOfSpecies-1
      write (*,"(A)",advance="no") trim(MolecularSystem_getSymbolOfSpecies(i))//", "
    end do
    write (*,"(A)",advance="no") trim(MolecularSystem_getSymbolOfSpecies(numberOfSpecies))
    write (*,*) ""

    write (*,"(A28)",advance="no") "  occupied orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)", advance="no") CIcore_instance%numberOfOccupiedOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  virtual orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") int(MolecularSystem_getTotalNumberOfContractions( i )* &
                                                CIcore_instance%lambda%values(i)  - &
                                                CIcore_instance%numberOfOccupiedOrbitals%values(i) )
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  total number of orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") int(MolecularSystem_getTotalNumberOfContractions( i )* &
                       CIcore_instance%lambda%values(i)   )
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  frozen core orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") CIcore_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  active occupied orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") CIcore_instance%numberOfOccupiedOrbitals%values(i) - &
                         CIcore_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") "  active virtual orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no") CIcore_instance%numberOfActiveOrbitals%values(i) - &
                         CIcore_instance%numberOfOccupiedOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") " total active orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no")  CIcore_instance%numberOfActiveOrbitals%values(i) - &
                          CIcore_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""
    write (*,*) " "

    !! printing header of the diagonalizers
    select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

      case ("JADAMILU")
        write(*,*) ""
        write(*,*) "  Diagonalizer : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
        write(*,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
        write(*,*) " a software code for computing selected eigenvalues of "
        write(*,*) " large sparse symmetric matrices, "
        write(*,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 

      case ("DSYEVX")
        write(*,*) ""
        write(*,*) "  Diagonalizer : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
        write(*,*) "LAPACK (Linear Algebra Package), standard software library for numerical linear algebra  "
        write(*,*) "https://netlib.org/lapack/"
        write(*,*) "DSYEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices"
        write (6,*) ""

    case ("DSYEVR")
        write(*,*) ""
        write(*,*) "  Diagonalizer : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
        write(*,*) "LAPACK (Linear Algebra Package), standard software library for numerical linear algebra  "
        write(*,*) "https://netlib.org/lapack/"
        write(*,*) "DSYEVR computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices"
        write (6,*) ""

    case default

      call CImod_exception( ERROR, "CImod run", "Diagonalization method not implemented")

    end select

    !! setting the requested CI level
    write (*, *) "Setting CI level..."
    call CIOrder_settingCILevel()

    !! just a rough estimate of the number of configurations in FCI. Just for info
    call CIOrder_estimate_FCI_numberOfConf()

    !! -------------------------------- Standard CI -------------------------------------
    if ( CONTROL_instance%CI_SELECTIVE_METHOD == "NONE" ) then


      write (*,*) "Building Strings..."
      call CIStrings_buildStrings()

      write (*,*) "Building CI level table..."
      call CIOrder_buildCIOrderList()

      if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) == "JADAMILU"  ) then
        !! additional arrays for matrix-vector diagonalizer
        call CIJadamilu_buildCouplingMatrix()
        call CIJadamilu_buildCouplingOrderList()
      endif

      !! getting the transformed AO to MO integrals, and transforming the one-particle integrals
      write (*, *) "Getting transformed integrals..."
      call CImod_getTransformedIntegrals()
      write (*,*) ""

      write (*, *) "Building diagonal..." !! and get number of configurations
      call CIDiag_buildDiagonal()

      if ( trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)) == "JADAMILU"  ) then
        !! initial guess of CI solutions for matrix-vector diagonalizer
        write (*,*) "Building initial hamiltonian..."
        call CIInitial_buildInitialCIMatrix2()
      endif

      !! allocating eigenValues array
      select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

      case ("JADAMILU")
        call Vector_constructor(CIcore_instance%eigenValues, &
                                int(CONTROL_instance%CI_NUMBER_OF_STATES, 8), 0.0_8)
      case ("DSYEVX")
        call Vector_constructor(CIcore_instance%eigenValues, &
                                int(CIcore_instance%numberOfConfigurations, 8), 0.0_8)

      case ("DSYEVR")
        call Vector_constructor(CIcore_instance%eigenValues, &
                                int(CIcore_instance%numberOfConfigurations, 8), 0.0_8)

      case default
        call CImod_exception(ERROR, "CImod run", "Diagonalization method not implemented")

      end select

      !! allocating eigenVector array
      call Matrix_constructor (CIcore_instance%eigenVectors, &
           int(CIcore_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%CI_NUMBER_OF_STATES,8), 0.0_8)

      if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
        call CImod_loadEigenVector (CIcore_instance%eigenvalues, &
               CIcore_instance%eigenVectors) 
      end if 

      !! diagonal correction. See 10.1016/j.chemphys.2007.07.001
      if ( CONTROL_instance%CI_DRESSING_SHIFT == "CISD") then

        call Vector_constructor(CIcore_instance%groundStateEnergies, 30_8, 0.0_8)
        call Vector_constructor(CIcore_instance%DDCISDTiming, 30_8, 0.0_8)
  
        write (6,*) ""
        write (6,"(T2,A50, A12)") "          ITERATIVE DIAGONAL DRESSED CISD SHIFT:   " , CONTROL_instance%CI_DRESSING_SHIFT
        write (6,"(T2,A62)")     "               ( Size-extensive correction)                   "
        write (6,"(T2,A62)")     " Based on 10.1016/j.chemphys.2007.07.001 and 10.1063/5.0182498"
        write (6,*) ""
  
        ecorr = 0.0_8
  
        do i = 2, 31
  
          !! add the diagonal shift
          do a = 2, CIcore_instance%numberOfConfigurations 
            CIcore_instance%diagonalHamiltonianMatrix%values(a) = CIcore_instance%diagonalHamiltonianMatrix%values(a) + ecorr
          end do
  
          select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
  
          case ("JADAMILU")

            write (6, *) ""
            write (6, "(T2,A,F14.5,A3 )") "Estimated memory needed: ", &
              real(CIcore_instance%numberOfConfigurations*(2 + (3*ms + CONTROL_instance%CI_NUMBER_OF_STATES + 1) + 4*ms*ms)*8)/(1024**3), " GB"
            write (6, *) ""

            call CIJadamilu_jadamiluInterface(CIcore_instance%numberOfConfigurations, &
               int(CONTROL_instance%CI_NUMBER_OF_STATES,8), &
                                              CIcore_instance%eigenValues, &
               CIcore_instance%eigenVectors, timeA, timeB)
  
            !! restore the original diagonal
            do a = 2, CIcore_instance%numberOfConfigurations 
              CIcore_instance%diagonalHamiltonianMatrix%values(a) = CIcore_instance%diagonalHamiltonianMatrix%values(a) - ecorr
            end do
  
          case ("DSYEVX")
            write (6, "(T2,A,F14.5,A3 )") "Estimated memory needed: ", &
              real((CIcore_instance%numberOfConfigurations**2 + 3)*8)/(1024**3), " GB"
            write (6, *) ""

            call CIFullMatrix_buildHamiltonianMatrix(timeA, timeB)

            call CIFullMatrix_buildHamiltonianMatrix( timeA, timeB)
    
            !! adding the shift
            do a = 2, CIcore_instance%numberOfConfigurations 
              CIcore_instance%hamiltonianMatrix%values(a,a) = CIcore_instance%hamiltonianMatrix%values(a,a) + ecorr
            end do
  
            call Matrix_eigen_select(CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenValues, &
               int(1), int(CONTROL_instance%CI_NUMBER_OF_STATES), &  
               eigenVectors = CIcore_instance%eigenVectors, &
               flags = int(SYMMETRIC,4))
  
          case ("DSYEVR")

            write (6, "(T2,A,F14.5,A3 )") "Estimated memory needed: ", &
              real((CIcore_instance%numberOfConfigurations**2 + 3)*8)/(1024**3), " GB"
            write (6, *) ""

            call CIFullMatrix_buildHamiltonianMatrix( timeA, timeB)
    
            do a = 2, CIcore_instance%numberOfConfigurations 
              CIcore_instance%hamiltonianMatrix%values(a,a) = CIcore_instance%hamiltonianMatrix%values(a,a) + ecorr
            end do
  
            call Matrix_eigen_dsyevr(CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenValues, &
                                     1_4, int(CONTROL_instance%CI_NUMBER_OF_STATES, 4), &
                   eigenVectors = CIcore_instance%eigenVectors, &
                   flags = SYMMETRIC)
  
          end select

          CIcore_instance%DDCISDTiming%values(i) = timeB - timeA
          CIcore_instance%groundStateEnergies%values(i) = CIcore_instance%eigenValues%values(1)

          !! current correlation energy
          if ( i == 2 ) then
            ecorr = CIcore_instance%groundStateEnergies%values(i) - HartreeFock_instance%totalEnergy
          else
            !! damped correlation energy
            ecorr = (1.0_8 - 0.5_8 ) * ( CIcore_instance%groundStateEnergies%values(i) - HartreeFock_instance%totalEnergy ) + &
                    ( 0.5_8 ) * ( CIcore_instance%groundStateEnergies%values(i-1) - HartreeFock_instance%totalEnergy )
          endif

          write (6, "(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i - 1, CIcore_instance%groundStateEnergies%values(i), &
                                       ecorr, &
                                       (CIcore_instance%groundStateEnergies%values(i - 1) - CIcore_instance%groundStateEnergies%values(i)), &
                                       timeB - timeA

          !! Restart ci matrix diagonalization from previous eigenvectors
          CONTROL_instance%CI_LOAD_EIGENVECTOR = .True.
  
          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit
  
        end do !! loop iterative CI
  
        !! printing the results of iterative CI
        write (6,*) ""
        write (6,"(T2,A42 )")    "  ITERATIVE DIAGONAL DRESSED CONVERGENCE  "
        write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "
        do i = 2, 31

          !! current correlation energy
          if ( i == 2 ) then
            ecorr = CIcore_instance%groundStateEnergies%values(i) - HartreeFock_instance%totalEnergy
          else
            !! damped correlation energy
            ecorr = (1.0_8 - 0.5_8 ) * ( CIcore_instance%groundStateEnergies%values(i) - HartreeFock_instance%totalEnergy ) + &
                    ( 0.5_8 ) * ( CIcore_instance%groundStateEnergies%values(i-1) - HartreeFock_instance%totalEnergy )
          endif

          write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i-1, CIcore_instance%groundStateEnergies%values(i), &
                                       ecorr, &
                                       (CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i)), &
                                       CIcore_instance%DDCISDTiming%values(i)

          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit
        end do
  
        if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
          call CImod_saveEigenVector () 
        end if
  
      !!-----------------------------------------------
      else !! no diagonal correction
  
         select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
  
        case ("JADAMILU")
  
          call CIJadamilu_jadamiluInterface(CIcore_instance%numberOfConfigurations, &
               int(CONTROL_instance%CI_NUMBER_OF_STATES,8), &
                                            CIcore_instance%eigenValues, &
               CIcore_instance%eigenVectors, timeA, timeB )
  
          if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
            call CImod_saveEigenVector () 
          end if
  
        case ("DSYEVX")
  
          call CIFullMatrix_buildHamiltonianMatrix(timeA, timeB)
          !$ write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"
    
          call Matrix_eigen_select(CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenValues, &
                 int(1), int(CONTROL_instance%CI_NUMBER_OF_STATES), &  
                 eigenVectors = CIcore_instance%eigenVectors, &
                 flags = int(SYMMETRIC,4))
  
        case ("DSYEVR")
  
          call CIFullMatrix_buildHamiltonianMatrix(timeA, timeB)
            !$ write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"
    
          call Matrix_eigen_dsyevr(CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenValues, &
                                   1_4, int(CONTROL_instance%CI_NUMBER_OF_STATES, 4), &
               eigenVectors = CIcore_instance%eigenVectors, &
               flags = SYMMETRIC)
  
         end select
  
       endif !! standard CI methods 
  
    !! -------------------------------- SCI -------------------------------------
    else !if ( CONTROL_instance%CI_SELECTIVE_METHOD == "SCI" ) then

      !! MCSCF
      if ( CONTROL_instance%CI_MCSCF ) then

        call CIMCSCF_show()
        call CIMCSCF_constructor( MCSCF_instance )

          do k = 1, CONTROL_instance%CI_MCSCF_MAX_ITER 
          !! getting the transformed AO to MO integrals, and transforming the one-particle integrals
          write (*, *) "Getting transformed integrals..."
          call CImod_getTransformedIntegrals()
          write (*,*) ""

          call CISCI_show()

          write (*,*) "allocating arrays for sci ..."
          call CISCI_constructor( CIcore_instance%numberOfConfigurations, &
                                  CIcore_instance%eigenValues, &
                                  CIcore_instance%eigenVectors )

          call CISCI_run( CIcore_instance%numberOfConfigurations, CIcore_instance%eigenVectors, &
                          initialEnergy = HartreeFock_instance%totalEnergy, &
                          initialStep = .true., unboundReference = .false., computePT2 = .false. )

          call CISCI_destructor()

          call CIMCSCF_compute( MCSCF_instance )

          if ( abs ( MCSCF_instance%energyChange%values( k + 1 ) ) <= 1E-5 ) exit

        enddo ! MCSCF iter

        call CIMCSCF_summary( MCSCF_instance )

        call CIMCSCF_destructor( MCSCF_instance )

      endif ! MCSCF macro

      !! single SCI ( or final SCI after MCSCF)

      !! getting the transformed AO to MO integrals, and transforming the one-particle integrals
      write (*, *) "Getting transformed integrals..."
      call CImod_getTransformedIntegrals()
      write (*,*) ""

      call CISCI_show()

      write (*,*) "Allocating arrays for SCI ..."
      call CISCI_constructor( CIcore_instance%numberOfConfigurations, &
                                  CIcore_instance%eigenValues, &
                                  CIcore_instance%eigenVectors )

      if ( CONTROL_instance%CI_UNBOUND_REFERENCE ) then
        call CISCI_run( CIcore_instance%numberOfConfigurations, CIcore_instance%eigenVectors, &
                        initialEnergy = HartreeFock_instance%totalEnergy, &
                        initialStep = .true., unboundReference = .true., computePT2 = .false. ) ! do a cisd- first
        call CISCI_run( CIcore_instance%numberOfConfigurations, CIcore_instance%eigenVectors, &
                        initialEnergy = CIcore_instance%eigenValues%values(1), &
                        initialStep = .false., unboundReference = .false., computePT2 = CONTROL_instance%CI_SCI_PT2_CORRECTION ) ! fci
      else  
        call CISCI_run( CIcore_instance%numberOfConfigurations, CIcore_instance%eigenVectors, &
                        initialEnergy = HartreeFock_instance%totalEnergy, &
                        initialStep = .true., unboundReference = .false., computePT2 = CONTROL_instance%CI_SCI_PT2_CORRECTION ) 
      endif

      call CISCI_destructor()

      call CISCI_saveEigenVector ( CIcore_instance%eigenVectors )

    end if !! SCI or not SCI

    write(6,*) ""
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) "          END ", trim(CIcore_instance%level)," CALCULATION"
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""
         
!    case ( "FCI-oneSpecie" )
!
!       print *, ""
!       print *, ""
!       print *, "==============================================="
!       print *, "|  Full CI for one specie calculation          |"
!       print *, "|  Use fci program to perform the calculation  |"
!       print *, "-----------------------------------------------"
!       print *, ""
!       ! call CIcore_getTransformedIntegrals()
!       !call CIcore_printTransformedIntegralsToFile()

  end subroutine CImod_run

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CImod_getTransformedIntegrals()
    implicit none

    integer :: numberOfSpecies
    integer :: i, j, m, n, mu, nu
    integer(8) :: a, b, c
    integer :: speciesID
    integer :: otherSpeciesID
    character(10) :: nameOfSpecies
    character(10) :: nameOfOtherSpecies
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecies
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecies
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients
    real(8) :: charge
    real(8) :: otherSpeciesCharge
    integer :: ssize1
    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if (allocated(CIcore_instance%twoCenterIntegrals)) deallocate (CIcore_instance%twoCenterIntegrals)
    if (allocated(CIcore_instance%fourCenterIntegrals)) deallocate (CIcore_instance%fourCenterIntegrals)
    if (allocated(CIcore_instance%twoIndexArray)) deallocate (CIcore_instance%twoIndexArray)
    if (allocated(CIcore_instance%fourIndexArray)) deallocate (CIcore_instance%fourIndexArray)

    allocate(CIcore_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(CIcore_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))

    allocate(CIcore_instance%twoIndexArray(numberOfSpecies))
    allocate(CIcore_instance%fourIndexArray(numberOfSpecies))

    do i=1, numberOfSpecies
      nameOfSpecies= trim(  MolecularSystem_getNameOfSpecies( i ) )
      speciesID = MolecularSystem_getSpeciesID( nameOfSpecies=nameOfSpecies )
      ocupationNumber = MolecularSystem_getOcupationNumber( i )
      numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
      charge=MolecularSystem_getCharge(i)

      call Matrix_constructor (CIcore_instance%twoCenterIntegrals(i), &
        int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

      call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

      !! Open file for wavefunction

      wfnFile = "lowdin.wfn"
      wfnUnit = 20

      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      arguments(2) = MolecularSystem_getNameOfSpecies(i)
      !if ( i == 2 ) arguments(2) =   MolecularSystem_getNameOfSpecies(1)
      !if ( i == 3 ) arguments(2) =   MolecularSystem_getNameOfSpecies(1)
      arguments(1) = "COEFFICIENTS"

      coefficients = &
        Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                           columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))

      arguments(1) = "HCORE"

      hcoreMatrix = &
        Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                           columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))

      !! transform two center integrals (one body operators)
        do m=1,numberOfContractions
          do n=m, numberOfContractions
             do mu=1, numberOfContractions
                do nu=1, numberOfContractions
                    CIcore_instance%twoCenterIntegrals(i)%values(m,n) = &
                        CIcore_instance%twoCenterIntegrals(i)%values(m,n) + &
                        coefficients%values(mu,m)* &
                        coefficients%values(nu,n)* &
                        hcoreMatrix%values(mu,nu)
            end do
          end do
        end do
      end do

      !! symmetrization
      do m = 1,numberOfContractions
        do n = m, numberOfContractions
          CIcore_instance%twoCenterIntegrals(i)%values(n,m)=&
                  CIcore_instance%twoCenterIntegrals(i)%values(m,n)
        end do
      end do

      !! auxilary 2-index array
      call Matrix_constructorInteger8(CIcore_instance%twoIndexArray(i), &
                          int( numberOfContractions,8), int( numberOfContractions,8) , 0_8 )

      c = 0_8
      do a = 1, numberOfContractions
        do b = a, numberOfContractions
          c = c + 1_8
          CIcore_instance%twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          CIcore_instance%twoIndexArray(i)%values(b,a) = CIcore_instance%twoIndexArray(i)%values(a,b)
        end do 
      end do

      !! auxilary 4-index array
      ssize1 = MolecularSystem_getTotalNumberOfContractions( i )
      ssize1 = ( ssize1 * ( ssize1 + 1_8 ) ) / 2_8

      call Matrix_constructorInteger8(CIcore_instance%fourIndexArray(i), &
                          int( ssize1,8), int( ssize1,8) , 0_8 )
      c = 0_8
      do a = 1, ssize1
        do b = a, ssize1
          c = c + 1_8
          CIcore_instance%fourIndexArray(i)%values(a,b) = c! IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          CIcore_instance%fourIndexArray(i)%values(b,a) = &
               CIcore_instance%fourIndexArray(i)%values(a,b)
         end do 
       end do

       call ReadTransformedIntegrals_readOneSpecies( speciesID, CIcore_instance%fourCenterIntegrals(i,i)   )
       CIcore_instance%fourCenterIntegrals(i,i)%values = &
           CIcore_instance%fourCenterIntegrals(i,i)%values * charge * charge

       if ( numberOfSpecies > 1 ) then
         do j = 1 , numberOfSpecies
           if ( i .ne. j) then
             nameOfOtherSpecies = trim(  MolecularSystem_getNameOfSpecies( j ) )
             otherSpeciesID = MolecularSystem_getSpeciesID( nameOfSpecies=nameOfOtherSpecies )
             ocupationNumberOfOtherSpecies = MolecularSystem_getOcupationNumber( j )
             numberOfContractionsOfOtherSpecies = MolecularSystem_getTotalNumberOfContractions( j )
             otherSpeciesCharge = MolecularSystem_getCharge(j)

             call ReadTransformedIntegrals_readTwoSpecies( speciesID, otherSpeciesID, &
                         CIcore_instance%fourCenterIntegrals(i,j) )
             CIcore_instance%fourCenterIntegrals(i,j)%values = &
               CIcore_instance%fourCenterIntegrals(i,j)%values * charge * otherSpeciescharge

           end if
         end do
       end if
     end do
     close (wfnUnit)
     call Matrix_destructor (hcoreMatrix)
    call Matrix_destructor(coefficients)

  end subroutine CImod_getTransformedIntegrals

  !**
  ! @ Retorna la energia final com correccion Moller-Plesset de orrden dado
  !**
  function CImod_getTotalEnergy() result(output)
    implicit none
    real(8) :: output

    output = CIcore_instance%totalEnergy

  end function CImod_getTotalEnergy

  subroutine CImod_saveEigenVector () 
    implicit none
    character(50) :: nameFile
    integer :: unitFile
    integer(8) :: i, ia
    integer :: ib, nonzero
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
    nameFile = "lowdin.civec"
    unitFile = 20

    nonzero = 0
    do i = 1, CIcore_instance%numberOfConfigurations
      if ( abs(CIcore_instance%eigenVectors%values(i,1) ) >= 1E-12 ) nonzero = nonzero + 1
    end do 

    write (*,*) "nonzero", nonzero

    allocate(auxArray(nonzero))
    allocate(auxIndexArray(nonzero))

    ia = 0
    do i = 1, CIcore_instance%numberOfConfigurations
      if ( abs(CIcore_instance%eigenVectors%values(i,1) ) >= 1E-12 ) then 
        ia = ia + 1
        auxIndexArray(ia) = i 
        auxArray(ia) = CIcore_instance%eigenVectors%values(i,1) 
      end if
    end do 

    open(unit=unitFile, file=trim(nameFile), status="replace", form="unformatted")

    write(unitFile) CIcore_instance%eigenValues%values(1)
    write(unitFile) nonzero

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxIndexArray(ia:ib)
    end do
    deallocate(auxIndexArray)

    do i = 1, ceiling(real(nonzero) / real(maxStackSize) )
      ib = maxStackSize * i  
      ia = ib - maxStackSize + 1
      if ( ib > nonzero ) ib = nonzero
      write(unitFile) auxArray(ia:ib)
    end do
    deallocate(auxArray)

    close(unitFile)

  end subroutine CImod_saveEigenVector

  subroutine CImod_loadEigenVector (eigenValues,eigenVectors) 
    implicit none
    type(Vector) :: eigenValues
    type(Matrix) :: eigenVectors
    character(50) :: nameFile
    integer :: unitFile
    integer :: i, ia, ib, nonzero
    real(8) :: eigenValue
    integer, allocatable :: auxIndexArray(:)
    real(8), allocatable :: auxArray(:)
    integer :: maxStackSize

    maxStackSize = CONTROL_instance%CI_STACK_SIZE 
 
    nameFile = "lowdin.civec"
    unitFile = 20

    open(unit=unitFile, file=trim(nameFile), status="old", action="read", form="unformatted")

    readvectors : do
      read (unitFile) eigenValue
      read (unitFile) nonzero
      write (*,*) "eigenValue", eigenValue
      write (*,*) "nonzero", nonzero

      allocate (auxIndexArray(nonzero))
      auxIndexArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxIndexArray(ia:ib)
      end do

      allocate (auxArray(nonzero))
      auxArray = 0

      do i = 1, ceiling(real(nonZero) / real(maxStackSize) )
        ib = maxStackSize * i  
        ia = ib - maxStackSize + 1
        if ( ib >  nonZero ) ib = nonZero
       read (unitFile) auxArray(ia:ib)
      end do
      exit readvectors
    end do readvectors

    eigenValues%values(1) = eigenValue
    do i = 1, nonzero
      eigenVectors%values(auxIndexArray(i),1) = auxArray(i)
    end do

    deallocate (auxIndexArray )
    deallocate (auxArray )

    close(unitFile)

  end subroutine CImod_loadEigenVector

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CImod_show()
    implicit none
    integer :: i
    integer(8) :: a
    real(8) :: davidsonCorrection, HFcoefficient, CIcorrection, MR
    integer numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( CIcore_instance%isInstanced ) then

      CIcorrection = CIcore_instance%eigenValues%values(1) - HartreeFock_instance%totalEnergy

      write(*,"(A)") " SUMMARY OF                       "
      write(*,"(A)") " POST HARTREE-FOCK CALCULATION    "
      write(*,"(A)") " CONFIGURATION INTERACTION THEORY:"
      write(6,*) "-----------------------------------------------------------------------"
      write(*,"(A)") ""
      write (6,"(T8,A30, A5)") "LEVEL = ", CIcore_instance%level
      write (6,"(T8,A30, I8)") "NUMBER OF CONFIGURATIONS = ", CIcore_instance%numberOfConfigurations
      write (6,"(T4,A34, F25.12)") "HF ENERGY = ", HartreeFock_instance%totalEnergy
      write (6,"(T4,A34, F25.12)") "GROUND STATE CORRELATION ENERGY = ", CIcorrection
      do i = 1, CONTROL_instance%CI_NUMBER_OF_STATES
        write (6, "(T18,A7,I3,A10, F25.12)") "STATE: ", i, " ENERGY = ", CIcore_instance%eigenValues%values(i)
      end do
      write(*,"(A)") ""

      if (  CIcore_instance%level == "CISD" ) then
        write(*,"(A)") ""
        write (6,"(T2,A34)") "RENORMALIZED DAVIDSON CORRECTION:"
        write(*,"(A)") ""
        write (6,"(T8,A54)") "E(CISDTQ) \approx E(CISD) + \delta E(Q)               "
        write (6,"(T8,A54)") "\delta E(Q) = (1 - c_0^2) * \delta E(CISD) / c_0^2    "
        write (*,*) ""
        HFcoefficient = CIcore_instance%eigenVectors%values(1,1) 
        davidsonCorrection = ( 1 - HFcoefficient*HFcoefficient) * CIcorrection / (HFcoefficient*HFcoefficient)
  
        write (6,"(T8,A19, F25.12)") "HF COEFFICIENT = ", HFcoefficient
        write (6,"(T8,A19, F25.12)") "\delta E(Q) = ", davidsonCorrection
        write (6,"(T8,A19, F25.12)") "E(CISDTQ) ESTIMATE ",  HartreeFock_instance%totalEnergy +&
           CIcorrection + davidsonCorrection
      endif

      if (  CONTROL_instance%CI_SELECTIVE_METHOD /= "NONE" ) then

        write(*,"(A)") ""
        write (6,"(T2,A34)") "EPSTEIN-NESBET PT2 CORRECTION:"
        write(*,"(A)") ""
        write (6,"(T8,A19, F25.12)") "E_PT2 :", CISCI_instance%PT2energy 
        write (6, "(T8,A19, F25.12)") "E_SCI + E_PT2 :", CIcore_instance%eigenValues%values(1) + CISCI_instance%PT2energy
      endif

      MR = 0.0_8
      do a = 1, CIcore_instance%numberOfConfigurations
        MR = MR +  abs( CIcore_instance%eigenVectors%values(a,1) )**2  - abs( CIcore_instance%eigenVectors%values(a,1) )**4 
      enddo

      write(*,"(A)") ""
      write(*,"(A)") "MULTI-REFERENCE CHARACTER ANALYSIS:"
      HFcoefficient = CIcore_instance%eigenVectors%values(1,1) 
      write(*,"(A)") ""
      write (6,"(T8,A19, F25.12)") "   HF COEFFICIENT = ", HFcoefficient
      write (6,"(T8,A19, F25.12)") "HF COEFFICIENT**2 = ", HFcoefficient**2
      write (6,"(T8,A19, F25.12)") "               MR = ", MR

      write(6,*) "-----------------------------------------------------------------------"

    end if

  end subroutine CImod_show

  subroutine CImod_showEigenVectors()
    implicit none

    integer(8) :: a,c
    integer :: p
    integer :: ci
    integer :: i
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer(8) :: numberOfConfigurations
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:)

    if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "NONE" ) return

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    numberOfConfigurations = CIcore_instance%numberOfConfigurations 

    write (*,*) ""
    write (*, "(T1,A)") "CI EIGENVECTORS" 
    write (*,*) ""
    write (*, "(T1,A,ES8.1)") "Printing coefficients larger than:", CONTROL_instance%CI_PRINT_THRESHOLD 
    write (*,*) ""

    if ( CONTROL_instance%CI_SELECTIVE_METHOD == "NONE" ) then

      allocate ( CIcore_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
      allocate ( indexConf ( numberOfSpecies ) )
      allocate (ciLevel(numberOfSpecies))
      ciLevel = 0
      CIcore_instance%allIndexConf = 0
      indexConf = 0
  
      !! gather all configurations
      s = 0
      c = 0
      ciLevel = 0
  
      do ci = 1,  CIcore_instance%sizeCiOrderList 
  
        cilevel(:) =  CIcore_instance%ciOrderList(  CIcore_instance%auxciOrderList(ci), :)
        s = 0
        auxnumberOfSpecies = CIcore_gatherConfRecursion( s, numberOfSpecies, indexConf,  c, cilevel )
      end do
  
      deallocate ( ciLevel )

      if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "ORBITALS" ) then
 
        do c = 1, CONTROL_instance%CI_NUMBER_OF_STATES
          write (*, "(T1,A,I4,A,F25.12)") "State: ", c, " Energy: ", CIcore_instance%eigenValues%values(c) 
          write (*, "(T1,A)") "Conf, orbital occupation per species, coefficient"
          write (*,*) ""
          do a = 1, numberOfConfigurations
            if ( abs(CIcore_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
              indexConf(:) = CIcore_instance%allIndexConf(:,a) 
  
              write (*, "(T1,I8,A1)", advance="no") a, " "
              do i = 1, numberOfSpecies
                do p = 1, CIcore_instance%numberOfActiveOrbitals%values(i)
                  write (*, "(I1)", advance="no")  CIcore_instance%orbitals(i)%values(p,indexConf(i)) 
                end do
                write (*, "(A1)", advance="no")  " "
              end do
              write (*, "(F11.8)") CIcore_instance%eigenVectors%values(a,c) 
            end if
          end do
          write (*,*) ""
        end do
  
      else if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "OCCUPIED" ) then
  
        do c = 1, CONTROL_instance%CI_NUMBER_OF_STATES
          write (*, "(T1,A,I4,A,F25.12)") "State: ", c, " Energy: ", CIcore_instance%eigenValues%values(c) 
          write (*, "(T1,A)") "Conf, occupied orbitals per species, coefficient"
          write (*,*) ""
          do a = 1, numberOfConfigurations
            if ( abs(CIcore_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
              indexConf(:) = CIcore_instance%allIndexConf(:,a) 
  
              write (*, "(T1,I8,A1)", advance="no") a, " "
              do i = 1, numberOfSpecies
                do p = 1, CIcore_instance%numberOfOccupiedOrbitals%values(i)
                  write (*, "(I3,A1)", advance="no") CIcore_instance%strings(i)%values(p,indexConf(i) ), " "
                end do
                write (*, "(A1)", advance="no")  "|"
              end do
              write (*, "(A1,F11.8)") " ", CIcore_instance%eigenVectors%values(a,c) 
            end if
          end do
          write (*,*) ""
        end do
  
      end if
  
      deallocate ( indexConf )
      deallocate ( CIcore_instance%allIndexConf )

    else !if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL == "SCI" ) then

      if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "ORBITALS" ) then
  
        do c = 1, CONTROL_instance%CI_NUMBER_OF_STATES
          write (*, "(T1,A,I4,A,F25.12)") "State: ", c, " Energy: ", CIcore_instance%eigenValues%values(c) 
          write (*, "(T1,A)") "Conf, orbital occupation per species, coefficient"
          write (*,*) ""
          do a = 1, numberOfConfigurations
            if ( abs(CIcore_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
              write (*, "(T1,I8,A1)", advance="no") a, " "
              do i = 1, numberOfSpecies
                do p = 1, CIcore_instance%numberOfActiveOrbitals%values(i)
                  write (*, "(I1)", advance="no") CISCI_instance%confTarget_orb(i)%values(p,a)
                                                  !CISCI_instance%targetOrb(i,a)%values(p)
                end do
                write (*, "(A1)", advance="no")  " "
              end do
              write (*, "(F11.8)") CIcore_instance%eigenVectors%values(a,c) 
            end if
          end do
          write (*,*) ""
        end do
  
      else if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "OCCUPIED" ) then

        do c = 1, CONTROL_instance%CI_NUMBER_OF_STATES
          write (*, "(T1,A,I4,A,F25.12)") "State: ", c, " Energy: ", CIcore_instance%eigenValues%values(c) 
          write (*, "(T1,A)") "Conf, occupied orbitals per species, coefficient"
          write (*,*) ""
          do a = 1, numberOfConfigurations
            if ( abs(CIcore_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
  
              write (*, "(T1,I8,A1)", advance="no") a, " "
              do i = 1, numberOfSpecies
                do p = 1, CIcore_instance%numberOfActiveOrbitals%values(i)
                  if ( CISCI_instance%confTarget_orb(i)%values(p,a)  == 1 ) then
                    !CISCI_instance%targetOrb(i,a)%values(p)
                    write (*, "(I3,A1)", advance="no") p, " "
                  endif
                end do
                write (*, "(A1)", advance="no")  "|"
              end do
              write (*, "(A1,F11.8)") " ", CIcore_instance%eigenVectors%values(a,c) 
            end if
          end do
          write (*,*) ""
        end do
  
      end if

    endif

  end subroutine CImod_showEigenVectors

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine CImod_exception( typeMessage, description, debugDescription)
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

  end subroutine CImod_exception

  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine CImod_destructor()
    implicit none
    integer c
    integer numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !!Destroy configurations
    !!Ground State
    if (allocated(CIcore_instance%configurations)) then
      c=1
      call Configuration_destructor(CIcore_instance%configurations(c) )
  
      do c=2, CIcore_instance%numberOfConfigurations
         call Configuration_destructor(CIcore_instance%configurations(c) )                
      end do
  
      if (allocated(CIcore_instance%configurations)) deallocate(CIcore_instance%configurations)
    end if

    call Vector_destructor(CIcore_instance%diagonalHamiltonianMatrix2)
    call Matrix_destructor(CIcore_instance%hamiltonianMatrix)
    call Vector_destructorInteger8(CIcore_instance%numberOfOccupiedOrbitals)
    call Vector_destructorInteger8(CIcore_instance%numberOfActiveOrbitals)
    call Vector_destructor (CIcore_instance%lambda)

    call Matrix_destructor(CIcore_instance%eigenVectors)
    call Vector_destructor(CIcore_instance%eigenValues)

    if (allocated(CIcore_instance%twoCenterIntegrals)) deallocate (CIcore_instance%twoCenterIntegrals)
    if (allocated(CIcore_instance%fourCenterIntegrals)) deallocate (CIcore_instance%fourCenterIntegrals)
    if (allocated(CIcore_instance%twoIndexArray)) deallocate (CIcore_instance%twoIndexArray)
    if (allocated(CIcore_instance%fourIndexArray)) deallocate (CIcore_instance%fourIndexArray)

    CIcore_instance%isInstanced=.false.

  end subroutine CImod_destructor

end module CImod_

