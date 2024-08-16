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
  use MolecularSystem_
  use String_
  use IndexMap_
  use InputCI_
  use omp_lib
  use JadamiluInterface_
  use CIcore_
  use CIDiag_
  use CIFullMatrix_
  use CIInitial_
  use CIJadamilu_
  use CIOrder_
  use CIStrings_

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
       CImod_densityMatrices, &
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
    integer :: i,j,m, numberOfSpecies
    integer :: a
    real(8) :: timeA, timeB
    real(8), allocatable :: eigenValues(:) 
    real(8) :: ecorr

!    select case ( trim(CIcore_instance%level) )
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    write (*,*) ""
    write (*,*) "==============================================="
    write (*,*) "         BEGIN ", trim(CIcore_instance%level)," CALCULATION"
    write (*,*) "         J. Charry, F. Moncada                "
    write (*,*) "-----------------------------------------------"
    write (*,*) ""

    write (*,"(A32)",advance="no") "Number of orbitals for species: "
    do i = 1, numberOfSpecies-1
      write (*,"(A)",advance="no") trim(MolecularSystem_getNameOfSpecie(i))//", "
    end do
    write (*,"(A)",advance="no") trim(MolecularSystem_getNameOfSpecie(numberOfSpecies))
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
      write (*,"(I5)",advance="no") CIcore_instance%numberOfOrbitals%values(i) - &
                         CIcore_instance%numberOfOccupiedOrbitals%values(i) 
    end do
    write (*,*) ""

    write (*,"(A28)",advance="no") " total active orbitals: "
    do i = 1, numberOfSpecies
      write (*,"(I5)",advance="no")  CIcore_instance%numberOfOrbitals%values(i) - &
                          CIcore_instance%numberOfCoreOrbitals%values(i) 
    end do
    write (*,*) ""
    write (*,*) " "

    write (*,*) "Getting transformed integrals..."
    call CImod_getTransformedIntegrals()
    write (*,*) " "

    !write (*,*) CIcore_instance%fourCenterIntegrals(1,1)%values(171, 1) a bug...
    write (*,*) "Setting CI level..."

    call CIOrder_settingCILevel()

   !! write (*,*) "Total number of configurations", CIcore_instance%numberOfConfigurations
    write (*,*) ""
    call Vector_constructor8 ( CIcore_instance%eigenvalues, &
                              int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8 )

    select case (trim(String_getUppercase(CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

    ! case ("ARPACK")

    !   write (*,*) "This method was removed"

    case ("JADAMILU")

      write (*,*) "Building Strings..."
      call CIStrings_buildStrings()

      write (*,*) "Building CI level table..."
      call CIOrder_buildCIOrderList()

      call CIJadamilu_buildCouplingMatrix()
      call CIJadamilu_buildCouplingOrderList()

      write (*,*) "Building diagonal..."
      call CIDiag_buildDiagonal()

      write (*,*) "Building initial hamiltonian..."
      call CIInitial_buildInitialCIMatrix2()

      call Matrix_constructor (CIcore_instance%eigenVectors, &
           int(CIcore_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      if ( CONTROL_instance%CI_LOAD_EIGENVECTOR ) then 
        call CImod_loadEigenVector (CIcore_instance%eigenvalues, &
               CIcore_instance%eigenVectors) 
      end if 

      write(*,*) ""
      write(*,*) "Diagonalizing hamiltonian..."
      write(*,*) "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))
      write(*,*) "============================================================="
      write(*,*) "M. BOLLHÖFER AND Y. NOTAY, JADAMILU:"
      write(*,*) " a software code for computing selected eigenvalues of "
      write(*,*) " large sparse symmetric matrices, "
      write(*,*) "Computer Physics Communications, vol. 177, pp. 951-964, 2007." 
      write(*,*) "============================================================="

      !! diagonal correction. See 10.1016/j.chemphys.2007.07.001
      if ( CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT == "CISD") then

      call Vector_constructor  (  CIcore_instance%groundStateEnergies, 30, 0.0_8)
      call Vector_constructor  (  CIcore_instance%DDCISDTiming, 30, 0.0_8)

        write (6,*) ""
        write (6,"(T2,A50, A12)") "          ITERATIVE DIAGONAL DRESSED CISD SHIFT:   " , CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT
        write (6,"(T2,A62)")     "               ( Size-extensive correction)                   "
        write (6,"(T2,A62)")     " Based on 10.1016/j.chemphys.2007.07.001 and 10.1063/5.0182498"
        write (6,*) ""

        ecorr = 0.0_8

        do i = 2, 31
  
          !! add the diagonal shift
          do a = 2, CIcore_instance%numberOfConfigurations 
            CIcore_instance%diagonalHamiltonianMatrix%values(a) = CIcore_instance%diagonalHamiltonianMatrix%values(a) + ecorr
          end do

          call CIJadamilu_jadamiluInterface(CIcore_instance%numberOfConfigurations, &
             int(CONTROL_instance%NUMBER_OF_CI_STATES,8), &
             CIcore_instance%eigenvalues, &
             CIcore_instance%eigenVectors, timeA, timeB)

          !! restore the original diagonal
          do a = 2, CIcore_instance%numberOfConfigurations 
            CIcore_instance%diagonalHamiltonianMatrix%values(a) = CIcore_instance%diagonalHamiltonianMatrix%values(a) - ecorr
          end do

          ecorr = CIcore_instance%eigenvalues%values(1)   - HartreeFock_instance%totalEnergy
          CIcore_instance%groundStateEnergies%values(i) = CIcore_instance%eigenvalues%values(1) 
          CIcore_instance%DDCISDTiming%values(i) = timeB - timeA

          write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i-1, CIcore_instance%groundStateEnergies%values(i), ecorr, (CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i)) , timeB - timeA

          !! Restart ci matrix diagonalization from previous eigenvectors
          CONTROL_instance%CI_LOAD_EIGENVECTOR = .True.

          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit

        end do 


        write (6,*) ""
        write (6,"(T2,A42 )")    "  ITERATIVE DIAGONAL DRESSED CONVERGENCE  "
        write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "
        do i = 2, 31
          write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i-1, CIcore_instance%groundStateEnergies%values(i), ecorr, (CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i)) , CIcore_instance%DDCISDTiming%values(i)
          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit
        end do

        if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
          call CImod_saveEigenVector () 
        end if

      else !! standard CI, no diagonal correction

        call CIJadamilu_jadamiluInterface(CIcore_instance%numberOfConfigurations, &
             int(CONTROL_instance%NUMBER_OF_CI_STATES,8), &
             CIcore_instance%eigenvalues, &
             CIcore_instance%eigenVectors, timeA, timeB )
  
        if ( CONTROL_instance%CI_SAVE_EIGENVECTOR ) then 
          call CImod_saveEigenVector () 
        end if

      end if

    case ("DSYEVX")

      write (*,*) "Building Strings..."
      call CIStrings_buildStrings()

      write (*,*) "Building CI level table..."
      call CIOrder_buildCIOrderList()

      write (*,*) "Building diagonal..."
      call CIDiag_buildDiagonal()

      call Matrix_constructor (CIcore_instance%eigenVectors, &
           int(CIcore_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      write(*,*) ""
      write(*,*) "Diagonalizing hamiltonian..."
      write(*,*) "  Using : ", trim(String_getUppercase((CONTROL_instance%CI_DIAGONALIZATION_METHOD)))

      !! diagonal correction. See 10.1016/j.chemphys.2007.07.001
      if ( CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT == "CISD") then

      call Vector_constructor  (  CIcore_instance%groundStateEnergies, 30, 0.0_8)

        write (6,*) ""
        write (6,"(T2,A50, A12)") "          ITERATIVE DIAGONAL DRESSED CISD SHIFT:   " , CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT
        write (6,"(T2,A62)")     "               ( Size-extensive correction)                   "
        write (6,"(T2,A62)")     " Based on 10.1016/j.chemphys.2007.07.001 and 10.1063/5.0182498"
        write (6,*) ""
        write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "

        ecorr = 0.0_8

        do i = 2, 31
  
          call CIFullMatrix_buildHamiltonianMatrix( timeA, timeB)
  
          do a = 2, CIcore_instance%numberOfConfigurations 
            CIcore_instance%hamiltonianMatrix%values(a,a) = CIcore_instance%hamiltonianMatrix%values(a,a) + ecorr
          end do

  
          call Matrix_eigen_select (CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenvalues, &
             int(1), int(CONTROL_instance%NUMBER_OF_CI_STATES), &  
             eigenVectors = CIcore_instance%eigenVectors, &
             flags = int(SYMMETRIC,4))
  
          ecorr = CIcore_instance%eigenvalues%values(1)   - HartreeFock_instance%totalEnergy
          CIcore_instance%groundStateEnergies%values(i) = CIcore_instance%eigenvalues%values(1) 

          write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i-1, CIcore_instance%groundStateEnergies%values(i), ecorr, (CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i)) , timeB - timeA

          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit

        end do 

      else !! standard CI, no diagonal correction

        call CIFullMatrix_buildHamiltonianMatrix(timeA, timeB)
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"
  
        call Matrix_eigen_select (CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenvalues, &
             int(1), int(CONTROL_instance%NUMBER_OF_CI_STATES), &  
             eigenVectors = CIcore_instance%eigenVectors, &
             flags = int(SYMMETRIC,4))

      end if

      !! deallocate transformed integrals
      deallocate(CIcore_instance%twoCenterIntegrals)
      deallocate(CIcore_instance%fourCenterIntegrals)

    case ("DSYEVR")

      write (*,*) "Building Strings..."
      call CIStrings_buildStrings()

      write (*,*) "Building CI level table..."
      call CIOrder_buildCIOrderList()

      write (*,*) "Building diagonal..."
      call CIDiag_buildDiagonal()

      call Matrix_constructor (CIcore_instance%eigenVectors, &
           int(CIcore_instance%numberOfConfigurations,8), &
           int(CONTROL_instance%NUMBER_OF_CI_STATES,8), 0.0_8)

      !! diagonal correction. See 10.1016/j.chemphys.2007.07.001
      if ( CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT == "CISD") then

      call Vector_constructor  (  CIcore_instance%groundStateEnergies, 30, 0.0_8)

        write (6,*) ""
        write (6,"(T2,A50, A12)") "          ITERATIVE DIAGONAL DRESSED CISD SHIFT:   " , CONTROL_instance%CI_DIAGONAL_DRESSED_SHIFT
        write (6,"(T2,A62)")     "               ( Size-extensive correction)                   "
        write (6,"(T2,A62)")     " Based on 10.1016/j.chemphys.2007.07.001 and 10.1063/5.0182498"
        write (6,*) ""
        write (6,"(T2,A95 )")    "Iter      Ground-State Energy       Correlation Energy           Energy Diff.          Time(s) "

        ecorr = 0.0_8

        do i = 2, 31
  
          call CIFullMatrix_buildHamiltonianMatrix( timeA, timeB)
  
          do a = 2, CIcore_instance%numberOfConfigurations 
            CIcore_instance%hamiltonianMatrix%values(a,a) = CIcore_instance%hamiltonianMatrix%values(a,a) + ecorr
          end do

         call Matrix_eigen_dsyevr (CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenvalues, &
               1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
               eigenVectors = CIcore_instance%eigenVectors, &
               flags = SYMMETRIC)
  
          ecorr = CIcore_instance%eigenvalues%values(1)   - HartreeFock_instance%totalEnergy
          CIcore_instance%groundStateEnergies%values(i) = CIcore_instance%eigenvalues%values(1) 

          write (6,"(T2,I2, F25.12, F25.12, F25.12, F16.4 )") i-1, CIcore_instance%groundStateEnergies%values(i), ecorr, (CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i)) , timeB - timeA

          if ( abs( CIcore_instance%groundStateEnergies%values(i-1) - CIcore_instance%groundStateEnergies%values(i) ) <= 1e-6) exit

        end do 

      else !! standard CI, no diagonal correction

        call CIFullMatrix_buildHamiltonianMatrix(timeA, timeB)
!$  write(*,"(A,E10.3,A4)") "** TOTAL Elapsed Time for building Hamiltonian Matrix : ", timeB - timeA ," (s)"
  
        call Matrix_eigen_dsyevr (CIcore_instance%hamiltonianMatrix, CIcore_instance%eigenvalues, &
             1, CONTROL_instance%NUMBER_OF_CI_STATES, &  
             eigenVectors = CIcore_instance%eigenVectors, &
             flags = SYMMETRIC)

      end if

      !! deallocate transformed integrals
      deallocate(CIcore_instance%twoCenterIntegrals)
      deallocate(CIcore_instance%fourCenterIntegrals)

    case default

      call CImod_exception( ERROR, "CImod run", "Diagonalization method not implemented")


    end select

    write(*,*) ""
    write(*,*) "-----------------------------------------------"
    write(*,*) "          END ", trim(CIcore_instance%level)," CALCULATION"
    write(*,*) "==============================================="
    write(*,*) ""
         
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
!

  end subroutine CImod_run

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine CImod_getTransformedIntegrals()
    implicit none

    integer :: numberOfSpecies
    integer :: i,j,m,n,mu,nu,a,b
    integer(8) :: c
    integer :: specieID
    integer :: otherSpecieID
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    integer :: ocupationNumber
    integer :: ocupationNumberOfOtherSpecie
    integer :: numberOfContractions
    integer :: numberOfContractionsOfOtherSpecie
    type(Matrix) :: hcoreMatrix
    type(Matrix) :: coefficients
    real(8) :: charge
    real(8) :: otherSpecieCharge

    integer :: ssize1, ssize2
    type(Matrix) :: externalPotential

    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    allocate(CIcore_instance%twoCenterIntegrals(numberOfSpecies))
    allocate(CIcore_instance%fourCenterIntegrals(numberOfSpecies,numberOfSpecies))

    allocate(CIcore_instance%twoIndexArray(numberOfSpecies))
    allocate(CIcore_instance%fourIndexArray(numberOfSpecies))

    do i=1, numberOfSpecies
      nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( i ) )
      specieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfSpecie )
      ocupationNumber = MolecularSystem_getOcupationNumber( i )
      numberOfContractions = MolecularSystem_getTotalNumberOfContractions( i )
      charge=MolecularSystem_getCharge(i)

!        write (6,"(T10,A)")"ONE PARTICLE INTEGRALS TRANSFORMATION FOR: "//trim(nameOfSpecie)
      call Matrix_constructor (CIcore_instance%twoCenterIntegrals(i), &
        int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8 )

      call Matrix_constructor (hcoreMatrix,int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

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

      c = 0
      do a=1,numberOfContractions
        do b = a, numberOfContractions
          c = c + 1
          CIcore_instance%twoIndexArray(i)%values(a,b) = c !IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          CIcore_instance%twoIndexArray(i)%values(b,a) = CIcore_instance%twoIndexArray(i)%values(a,b)
        end do 
      end do


      !! auxilary 4-index array
      ssize1 = MolecularSystem_getTotalNumberOfContractions( i )
      ssize1 = ( ssize1 * ( ssize1 + 1 ) ) / 2

      call Matrix_constructorInteger8(CIcore_instance%fourIndexArray(i), &
                          int( ssize1,8), int( ssize1,8) , 0_8 )
      c = 0
      do a = 1, ssize1
        do b = a, ssize1
          c = c + 1
          CIcore_instance%fourIndexArray(i)%values(a,b) = c! IndexMap_tensorR2ToVectorC( a, b, numberOfContractions )
          CIcore_instance%fourIndexArray(i)%values(b,a) = &
               CIcore_instance%fourIndexArray(i)%values(a,b)
         end do 
       end do


       call ReadTransformedIntegrals_readOneSpecies( specieID, CIcore_instance%fourCenterIntegrals(i,i)   )
       CIcore_instance%fourCenterIntegrals(i,i)%values = &
           CIcore_instance%fourCenterIntegrals(i,i)%values * charge * charge

       if ( numberOfSpecies > 1 ) then
         do j = 1 , numberOfSpecies
           if ( i .ne. j) then
             nameOfOtherSpecie = trim(  MolecularSystem_getNameOfSpecie( j ) )
             otherSpecieID = MolecularSystem_getSpecieID( nameOfSpecie=nameOfOtherSpecie )
             ocupationNumberOfOtherSpecie = MolecularSystem_getOcupationNumber( j )
             numberOfContractionsOfOtherSpecie = MolecularSystem_getTotalNumberOfContractions( j )
             otherSpecieCharge = MolecularSystem_getCharge(j)

             call ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, &
                         CIcore_instance%fourCenterIntegrals(i,j) )
             CIcore_instance%fourCenterIntegrals(i,j)%values = &
               CIcore_instance%fourCenterIntegrals(i,j)%values * charge * otherSpeciecharge


           end if
         end do
       end if
     end do
     close (wfnUnit)
     call Matrix_destructor (hcoreMatrix)

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
    type(Vector8) :: eigenValues
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
    type(CIcore) :: this
    integer :: i
    real(8) :: davidsonCorrection, HFcoefficient, CIcorrection
    integer numberOfSpecies

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    if ( CIcore_instance%isInstanced ) then

       write(*,"(A)") ""
       write(*,"(A)") " POST HARTREE-FOCK CALCULATION"
       write(*,"(A)") " CONFIGURATION INTERACTION THEORY:"
       write(*,"(A)") "=============================="
       write(*,"(A)") ""
       write (6,"(T8,A30, A5)") "LEVEL = ", CIcore_instance%level
       write (6,"(T8,A30, I8)") "NUMBER OF CONFIGURATIONS = ", CIcore_instance%numberOfConfigurations
       do i = 1, CONTROL_instance%NUMBER_OF_CI_STATES
        write (6,"(T8,A17,I3,A10, F25.12)") "STATE: ", i, " ENERGY = ", CIcore_instance%eigenvalues%values(i)
       end do
       write(*,"(A)") ""
       CIcorrection = CIcore_instance%eigenvalues%values(1) - &
                HartreeFock_instance%totalEnergy

       write (6,"(T4,A34, F25.12)") "GROUND STATE CORRELATION ENERGY = ", CIcorrection

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
       else 

         write(*,"(A)") ""
         HFcoefficient = CIcore_instance%eigenVectors%values(1,1) 
         write (6,"(T8,A19, F25.12)") "HF COEFFICIENT = ", HFcoefficient

       end if

    else 

    end if

  end subroutine CImod_show

  subroutine CImod_showEigenVectors()
    implicit none

    integer(8) :: a,b,c
    integer :: u,v,p
    integer :: ci
    integer :: i, j, ii, jj
    integer :: s, numberOfSpecies, auxnumberOfSpecies
    integer :: size1, size2
    real(8) :: timeA, timeB
    integer(1) :: coupling
    integer(8) :: numberOfConfigurations
    real(8) :: CIenergy
    integer(8), allocatable :: indexConf(:)
    integer, allocatable :: cilevel(:), auxcilevel(:), dd(:)


    if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "NONE" ) return

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    numberOfConfigurations = CIcore_instance%numberOfConfigurations 

    allocate ( CIcore_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
    allocate ( ciLevel ( numberOfSpecies ) )
    allocate ( indexConf ( numberOfSpecies ) )
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
    !stop

    deallocate ( ciLevel )

    if ( CONTROL_instance%CI_PRINT_EIGENVECTORS_FORMAT == "ORBITALS" ) then
    write (*,*) ""
    write (*, "(T1,A)") "Eigenvectors" 
    write (*,*) ""

    do c = 1, CONTROL_instance%NUMBER_OF_CI_STATES
      write (*, "(T1,A,I4,A,F25.12)") "State: ", c, " Energy: ", CIcore_instance%eigenValues%values(c) 
      write (*, "(T1,A)") "Conf, orbital occupation per species, coefficient"
      write (*,*) ""
      do a = 1, numberOfConfigurations
        if ( abs(CIcore_instance%eigenVectors%values(a,c)) > CONTROL_instance%CI_PRINT_THRESHOLD ) then  
          indexConf(:) = CIcore_instance%allIndexConf(:,a) 

          write (*, "(T1,I8,A1)", advance="no") a, " "
          do i = 1, numberOfSpecies
            do p = 1, CIcore_instance%numberOfOrbitals%values(i)
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
    write (*,*) ""
    write (*, "(T1,A)") "Eigenvectors" 
    write (*,*) ""

    do c = 1, CONTROL_instance%NUMBER_OF_CI_STATES
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
          write (*, "(A,F11.8)") " ", CIcore_instance%eigenVectors%values(a,c) 
        end if
      end do
      write (*,*) ""
    end do

    end if

    deallocate ( indexConf )
    deallocate ( CIcore_instance%allIndexConf )

  end subroutine CImod_showEigenVectors


  !FELIX IS HERE
  subroutine CImod_densityMatrices()
    implicit none
    type(CIcore) :: this
    type(Configuration) :: auxthisA, auxthisB
    integer :: i, j, k, l, mu, nu, n
    integer :: factor
    integer :: unit, wfnunit
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state, species, orbital, orbitalA, orbitalB
    character(50) :: file, wfnfile, speciesName, auxstring
    character(50) :: arguments(2)
    type(matrix), allocatable :: coefficients(:), atomicDensityMatrix(:,:), ciDensityMatrix(:,:), auxDensMatrix(:,:)
    type(matrix), allocatable :: kineticMatrix(:), attractionMatrix(:), externalPotMatrix(:)
    integer numberOfSpecies

    type(matrix) :: auxdensityEigenVectors 
    type(matrix) :: densityEigenVectors
    type(vector) :: auxdensityEigenValues
    type(vector) :: densityEigenValues
    integer, allocatable :: cilevel(:), cilevelA(:)
    integer(8) :: numberOfConfigurations, c
    integer(8), allocatable :: indexConf(:)
    type(ivector), allocatable :: stringAinB(:)
    integer :: s, ss, ci, auxnumberOfSpecies
    integer, allocatable :: coupling(:)
    integer :: a, b, AA, BB, bj
    integer :: u, uu, ssize
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(8), allocatable :: jj(:)
    real(8) :: timeDA
    real(8) :: timeDB

  
    ! type(Vector) :: eigenValues
    ! type(Matrix) :: eigenVectors, auxMatrix
    ! real(8) :: sumaPrueba

    !!Iterators: i,j - Configurations .... k,l - molecular orbitals .... mu,nu - atomic orbitals ... n - threads
    if ( CIcore_instance%isInstanced .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0 ) then
       !$  timeDA = omp_get_wtime()

      numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
  
      numberOfConfigurations = CIcore_instance%numberOfConfigurations 
  
      allocate (stringAinB ( numberOfSpecies ))
  
      do i = 1, numberOfSpecies 
        call Vector_constructorInteger (stringAinB(i), CIcore_instance%numberOfOccupiedOrbitals%values(i), 0)
      end do 
  
      allocate ( CIcore_instance%allIndexConf( numberOfSpecies, numberOfConfigurations ) )
      allocate ( ciLevelA ( numberOfSpecies ) )
      allocate ( ciLevel ( numberOfSpecies ) )
      allocate ( indexConf ( numberOfSpecies ) )
      ciLevelA = 0
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
      !stop
  
      deallocate ( indexConf )
      allocate ( coupling ( numberOfSpecies ) )


      write (*,*) ""
      write (*,*) "=============================="
      write (*,*) "BUILDING CI DENSITY MATRICES"
      write (*,*) "=============================="
      write (*,*) ""

      allocate( coefficients(numberOfSpecies), &
           kineticMatrix(numberOfSpecies), &
           attractionMatrix(numberOfSpecies), &
           externalPotMatrix(numberOfSpecies), &
           atomicDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), &
           ciDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), &
           auxDensMatrix(numberOfSpecies,CIcore_instance%nproc) )

      wfnFile = "lowdin.wfn"
      wfnUnit = 20
      open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

      !Inicializando las matrices
      do species=1, numberOfSpecies
         speciesName = MolecularSystem_getNameOfSpecie(species)
         
         numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
         ! numberOfOrbitals = CIcore_instance%numberOfOrbitals%values(species)
         numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(species)

         arguments(2) = speciesName
         ! print *, "trolo", numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals

         arguments(1) = "COEFFICIENTS"
         coefficients(species) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
              columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

         arguments(1) = "KINETIC"
         kineticMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
              columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
         
         arguments(1) = "ATTRACTION"
         attractionMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
              columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

         arguments(1) = "EXTERNAL_POTENTIAL"
         if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
              externalPotMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
              columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))
         ! print *, "trololo"
        
         do state=1, CONTROL_instance%CI_STATES_TO_PRINT

            call Matrix_constructor ( ciDensityMatrix(species,state) , &
                 int(numberOfContractions,8), &
                 int(numberOfContractions,8),  0.0_8 )

            do k=1, numberOfOccupiedOrbitals
               ciDensityMatrix(species,state)%values( k, k)=1.0_8
            end do

         end do

         do n=1, CIcore_instance%nproc

            call Matrix_constructor ( auxDensMatrix(species,n) , &
                 int(numberOfContractions,8), &
                 int(numberOfContractions,8),  0.0_8 )
         end do
      end do
       
      close(wfnUnit)

      allocate ( indexConfA ( numberOfSpecies ) )
      allocate ( indexConfB ( numberOfSpecies ) )
      allocate ( jj ( numberOfSpecies ) )

      indexConfA = 0
      indexConfB = 0
      jj = 0

      !! Building the CI reduced density matrix in the molecular orbital representation in parallel
      ! call Matrix_show (CIcore_instance%eigenVectors)

      !!print *, "        State, Progress"
      
      do state=1, CONTROL_instance%CI_STATES_TO_PRINT

         !$omp parallel & 
         !$omp& firstprivate (stringAinB,indexConfA,indexConfB, jj) &
         !$omp& private(i,j, species, s, numberOfOccupiedOrbitals, k, coupling, orbital, orbitalA, orbitalB, AA, BB, a, b, factor, n, cilevelA, ss, ssize, cilevel, ci, u, uu, bj),&
         !$omp& shared(CIcore_instance, auxDensMatrix )
         n = omp_get_thread_num() + 1
         !$omp do schedule (dynamic) 
         do i=1, CIcore_instance%numberOfConfigurations

            !!if( mod( i , 50000 ) .eq. 0 ) print *, state, floor(real(100*i/CIcore_instance%numberOfConfigurations)), "%"
            !!Filter very small coefficients
            if( abs(CIcore_instance%eigenVectors%values(i,state)) .ge. 1E-10) then

               indexConfA(:) = CIcore_instance%allIndexConf(:,i) 

               !print *, "==", indexConfA , "|", i


               !!Diagonal contributions
               do species=1, numberOfSpecies
                  numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(species)

                  do k=1, numberOfOccupiedOrbitals

                     !!Occupied orbitals
                     auxDensMatrix(species,n)%values(k,k)=auxDensMatrix(species,n)%values(k,k) - CIcore_instance%eigenVectors%values(i,state)**2
                     ! ciDensityMatrix(species,state)%values( k, k) = ciDensityMatrix(species,state)%values( k, k) -  &
                     !      CIcore_instance%eigenVectors%values(i,state)**2

                     !print *, i, j, k, species 
                     !orbital = CIcore_instance%configurations(i)%occupations(k,species) 
                     orbital =  CIcore_instance%strings(species)%values(k,indexConfA(species))
                     !!Unoccupied orbitals

                     auxDensMatrix(species,n)%values(orbital,orbital)=auxDensMatrix(species,n)%values(orbital,orbital) + CIcore_instance%eigenVectors%values(i,state)**2
                     ! ciDensityMatrix(species,state)%values( orbital, orbital)= ciDensityMatrix(species,state)%values( orbital, orbital) + &
                     !      CIcore_instance%eigenVectors%values(i,state)**2

                  end do
               end do

               !!Off Diagonal contributions
               cilevelA = 0
               do ss = 1, numberOfSpecies 
                 stringAinB(ss)%values = 0
                 do k = 1, CIcore_instance%numberOfOccupiedOrbitals%values(ss)

                   stringAinB(ss)%values(k) = CIcore_instance%orbitals(ss)%values( &
                                             CIcore_instance%strings(ss)%values(k,  CIcore_instance%allIndexConf(ss,1)), indexConfA(ss))
                 end do
                 cilevelA(ss) = CIcore_instance%numberOfOccupiedOrbitals%values(ss) - sum ( stringAinB(ss)%values )
               end do 

               jj = 0
               coupling = 0
               do ss = 1, numberOfSpecies 
                 ssize = 0 

                 indexConfB(:) = indexConfA(:)
                 cilevel = cilevelA

                 do ci = 1,  size(CIcore_instance%numberOfStrings(ss)%values, dim = 1)
                   cilevel(ss) = ci - 1
                   do u = 1,  CIcore_instance%sizeCiOrderList 
                     if ( sum(abs(cilevel - &
                          CIcore_instance%ciOrderList( CIcore_instance%auxciOrderList(u), :))) == 0 ) then
                       uu = CIcore_instance%auxciOrderList(u)
                       do bj = 1 + ssize , CIcore_instance%numberOfStrings(ss)%values(ci) + ssize
                         indexConfB(ss) = bj
  
                         do s=1, numberOfSpecies
                           jj(s) = (indexConfB(s) - CIcore_instance%numberOfStrings2(s)%values(cilevel(s)+1) + &
                                    CIcore_instance%ciOrderSize1(uu,s) )* CIcore_instance%ciOrderSize2(uu,s) 
                         end do

                         j = sum(jj)
                         !print *, "  ", indexConfB , "|", j, CIcore_instance%eigenVectors%values(j,state) 
                         if ( j > i ) then
                           if( abs(CIcore_instance%eigenVectors%values(j,state)) .ge. 1E-10) then

                             coupling = 0
                             do s=1, numberOfSpecies
                                stringAinB(s)%values = 0
                                do k = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s)
                                   stringAinB(s)%values(k) = CIcore_instance%orbitals(s)%values( &
                                        CIcore_instance%strings(s)%values(k,indexConfA(s) ), indexConfB(s) ) 
                                end do
                                coupling(s) = CIcore_instance%numberOfOccupiedOrbitals%values(s) - sum ( stringAinB(s)%values )
                             end do
                             if (sum(coupling) == 1) then
    
                               do s = 1, numberOfSpecies
    
                                 if ( coupling(s) == 1) then !!hmm

                                   !print *, "      ", coupling
                                   orbitalA = 0
                                   orbitalB = 0
                                   AA = 0
                                   BB = 0
                                   a = indexConfA(s)
                                   b = indexConfB(s)
    
                                   do k = 1, CIcore_instance%occupationNumber(s) 
                                      if ( CIcore_instance%orbitals(s)%values( &
                                           CIcore_instance%strings(s)%values(k,a),b) == 0 ) then
                                         orbitalA =  CIcore_instance%strings(s)%values(k,a)
                                         AA = k
                                         exit
                                      end if
                                   end do
                                   do k = 1, CIcore_instance%occupationNumber(s) 
                                      if ( CIcore_instance%orbitals(s)%values( &
                                           CIcore_instance%strings(s)%values(k,b),a) == 0 ) then
                                         orbitalB =  CIcore_instance%strings(s)%values(k,b)
                                         BB = k
                                         exit
                                      end if
                                   end do
    
                                   factor = (-1)**(AA-BB)
    
                                   numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(s)
    
                                   ! print *, i, j, CIcore_instance%configurations(i)%occupations(:,specie), CIcore_instance%configurations(j)%occupations(:,specie)
                                   ! print *, i, j, auxthisA%occupations(:,specie), auxthisB%occupations(:,specie)
                                   ! print *, i, j, orbitalA, orbitalB, factor*CIcore_instance%eigenVectors%values(i,1)*CIcore_instance%eigenVectors%values(j,1)
    
                                   auxDensMatrix(s,n)%values( orbitalA,orbitalB)= auxDensMatrix(s,n)%values( orbitalA, orbitalB) + &
                                        factor*CIcore_instance%eigenVectors%values(i,state)* &
                                        CIcore_instance%eigenVectors%values(j,state)
                                   auxDensMatrix(s,n)%values( orbitalB,orbitalA)= auxDensMatrix(s,n)%values( orbitalB, orbitalA) + &
                                        factor*CIcore_instance%eigenVectors%values(i,state)* &
                                        CIcore_instance%eigenVectors%values(j,state)
                                  end if
                                end do
                              end if
                           end if
                         end if
                         !! here
                       end do
                       ssize = ssize + CIcore_instance%numberOfStrings(ss)%values(ci)
                       !exit
                     end if

                   end do
                 end do

               end do 

!               do j=i+1, CIcore_instance%numberOfConfigurations
!                  if( abs(CIcore_instance%eigenVectors%values(j,state)) .ge. 1E-12) then

!                     indexConfB(:) = CIcore_instance%allIndexConf(:,j)

!                     coupling = 0
!                     do s=1, numberOfSpecies
!                        stringAinB(s)%values = 0
!                        do k = 1, CIcore_instance%numberOfOccupiedOrbitals%values(s)
!                           stringAinB(s)%values(k) = CIcore_instance%orbitals(s)%values( &
!                                CIcore_instance%strings(s)%values(k,indexConfA(s) ), indexConfB(s) ) 
!                        end do
!                        coupling(s) = CIcore_instance%numberOfOccupiedOrbitals%values(s) - sum ( stringAinB(s)%values )
!                     end do
!
!                     if (sum(coupling) == 1) then
!
!                        do s = 1, numberOfSpecies
!
!                           if ( coupling(s) == 1) then
!                              orbitalA = 0
!                              orbitalB = 0
!                              AA = 0
!                              BB = 0
!                              a = indexConfA(s)
!                              b = indexConfB(s)
!
!                              do k = 1, CIcore_instance%occupationNumber(s) 
!                                 if ( CIcore_instance%orbitals(s)%values( &
!                                      CIcore_instance%strings(s)%values(k,a),b) == 0 ) then
!                                    orbitalA =  CIcore_instance%strings(s)%values(k,a)
!                                    AA = k
!                                    exit
!                                 end if
!                              end do
!                              do k = 1, CIcore_instance%occupationNumber(s) 
!                                 if ( CIcore_instance%orbitals(s)%values( &
!                                      CIcore_instance%strings(s)%values(k,b),a) == 0 ) then
!                                    orbitalB =  CIcore_instance%strings(s)%values(k,b)
!                                    BB = k
!                                    exit
!                                 end if
!                              end do
!
!                              factor = (-1)**(AA-BB)
!
!                              numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(s)
!
!                              ! print *, i, j, CIcore_instance%configurations(i)%occupations(:,specie), CIcore_instance%configurations(j)%occupations(:,specie)
!                              ! print *, i, j, auxthisA%occupations(:,specie), auxthisB%occupations(:,specie)
!
!                              ! print *, i, j, orbitalA, orbitalB, factor*CIcore_instance%eigenVectors%values(i,1)*CIcore_instance%eigenVectors%values(j,1)
!
!                              auxDensMatrix(s,n)%values( orbitalA,orbitalB)= auxDensMatrix(s,n)%values( orbitalA, orbitalB) + &
!                                   factor*CIcore_instance%eigenVectors%values(i,state)* &
!                                   CIcore_instance%eigenVectors%values(j,state)
!                              ! ciDensityMatrix(s,state)%values( orbitalA,orbitalB)= ciDensityMatrix(s,state)%values( orbitalA, orbitalB) + &
!                              !      factor*CIcore_instance%eigenVectors%values(i,state)* &
!                              !      CIcore_instance%eigenVectors%values(j,state)
!
!                              auxDensMatrix(s,n)%values( orbitalB,orbitalA)= auxDensMatrix(s,n)%values( orbitalB, orbitalA) + &
!                                   factor*CIcore_instance%eigenVectors%values(i,state)* &
!                                   CIcore_instance%eigenVectors%values(j,state)
!
!                              ! ciDensityMatrix(s,state)%values( orbitalB, orbitalA)= ciDensityMatrix(s,state)%values( orbitalB, orbitalA) + &
!                              !      factor*CIcore_instance%eigenVectors%values(i,state)* &
!                              !      CIcore_instance%eigenVectors%values(j,state)
!
!                           end if
!                        end do
!                     end if
!                  end if
!               end do

            end if
         end do
         !$omp end do nowait
         !$omp end parallel
         
         !! Gather the parallel results
         do species=1, numberOfSpecies
            do n=1, CIcore_instance%nproc
               ciDensityMatrix(species,state)%values = ciDensityMatrix(species,state)%values + auxDensMatrix(species,n)%values
               auxDensMatrix(species,n)%values=0.0
            end do
         end do

      end do


      !! Open file - to write density matrices
     unit = 29
       
     file = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
     open(unit = unit, file=trim(file), status="new", form="formatted")
       
     !! Building the CI reduced density matrix in the atomic orbital representation       
     do species=1, numberOfSpecies
       speciesName = MolecularSystem_getNameOfSpecie(species)
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )

       do state=1, CONTROL_instance%CI_STATES_TO_PRINT
          
         ! print *, "CI density matrix ", trim(speciesName), state
         ! call Matrix_show ( ciDensityMatrix(species,state))
             
         call Matrix_constructor ( atomicDensityMatrix(species,state) , &
                                   int(numberOfContractions,8), &
                                   int(numberOfContractions,8),  0.0_8 )

         do mu=1, numberOfContractions
            do nu=1, numberOfContractions
               do k=1, numberOfContractions
                  atomicDensityMatrix(species,state)%values(mu,nu) =  &
                       atomicDensityMatrix(species,state)%values(mu,nu) + &
                       ciDensityMatrix(species,state)%values(k,k) *&
                       coefficients(species)%values(mu,k)*coefficients(species)%values(nu,k)

                  do l=k+1, numberOfContractions

                     atomicDensityMatrix(species,state)%values(mu,nu) =  &
                          atomicDensityMatrix(species,state)%values(mu,nu) + &
                          ciDensityMatrix(species,state)%values(k,l) *&
                          (coefficients(species)%values(mu,k)*coefficients(species)%values(nu,l) + & 
                          coefficients(species)%values(mu,l)*coefficients(species)%values(nu,k))

                  end do
               end do
            end do
         end do
       
         ! print *, "atomic density matrix  ", trim(speciesName), state
         ! call Matrix_show ( atomicDensityMatrix(species,state))

         write(auxstring,*) state
         arguments(2) = speciesName
         arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
             
         call Matrix_writeToFile ( atomicDensityMatrix(species,state), unit , arguments=arguments(1:2) )

         end do
       end do

       write(*,*) ""
       write(*,*) "==============================="
       write(*,*) " ONE BODY ENERGY CONTRIBUTIONS:"
       write(*,*) ""
       do state=1, CONTROL_instance%CI_STATES_TO_PRINT
          write(*,*) " STATE: ", state
          do species=1, molecularSystem_instance%numberOfQuantumSpecies
             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%name ) // &
                  " Kinetic energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*kineticMatrix(species)%values)
             write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%name ) // &
                  "/Fixed interact. energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*attractionMatrix(species)%values)
             if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                  write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%name) // &
                  " Ext Pot energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*externalPotMatrix(species)%values)
             print *, ""
          end do
          print *, ""
       end do
 
      !! Natural orbitals

       if (CONTROL_instance%CI_NATURAL_ORBITALS) then

          write(*,*) ""
          write(*,*) "=============================="
          write(*,*) " NATURAL ORBITALS: "
          write(*,*) ""

          do state=1, CONTROL_instance%CI_STATES_TO_PRINT

             write(*,*) " STATE: ", state

             do species=1, numberOfSpecies

                write(*,*) ""
                write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(species)%name )
                write(*,*) "-----------------"

                numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
                speciesName = MolecularSystem_getNameOfSpecie(species)


                call Vector_constructor ( auxdensityEigenValues, &
                     int(numberOfContractions,4),  0.0_8 )

                call Matrix_constructor ( auxdensityEigenVectors, &
                     int(numberOfContractions,8), &
                     int(numberOfContractions,8),  0.0_8 )

                call Vector_constructor ( densityEigenValues, &
                     int(numberOfContractions,4),  0.0_8 )

                call Matrix_constructor ( densityEigenVectors, &
                     int(numberOfContractions,8), &
                     int(numberOfContractions,8),  0.0_8 )

                call Matrix_eigen ( ciDensityMatrix(species,state), auxdensityEigenValues, auxdensityEigenVectors, SYMMETRIC )  

                ! reorder and count significant occupations
                k=0
                do u = 1, numberOfContractions
                   densityEigenValues%values(u) =  auxdensityEigenValues%values(numberOfContractions - u + 1)
                   densityEigenVectors%values(:,u) = auxdensityEigenVectors%values(:,numberOfContractions - u + 1)
                   if(densityEigenValues%values(u) .ge. 5.0E-5 ) k=k+1
                end do

                !! Transform to atomic basis
                densityEigenVectors%values = matmul( coefficients(species)%values, densityEigenVectors%values )

                ! Print eigenvectors with occupation larger than 5.0E-5
                call Matrix_constructor(auxdensityEigenVectors,int(numberOfContractions,8),int(k,8),0.0_8)
                do u=1, numberOfContractions
                   do j=1, k
                      auxdensityEigenVectors%values(u,j)=densityEigenVectors%values(u,j)
                   end do
                end do
                call Matrix_show( auxdensityEigenVectors, &
                     rowkeys = MolecularSystem_getlabelsofcontractions( species ), &
                     columnkeys = string_convertvectorofrealstostring( densityEigenValues ),&
                     flags=WITH_BOTH_KEYS)

                write(auxstring,*) state
                arguments(2) = speciesName
                arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 

                call Matrix_writeToFile ( densityEigenVectors, unit , arguments=arguments(1:2) )
                arguments(1) = "OCCUPATIONS"//trim(adjustl(auxstring))

                call Vector_writeToFile( densityEigenValues, unit, arguments=arguments(1:2) )
                !! it's the same
                !!auxdensityEigenVectors%values = 0

                !!do mu=1, numberOfContractions
                !!  do nu=1, numberOfContractions
                !!    do k=1, numberOfContractions
                !!      auxdensityEigenVectors%values(mu,nu) = auxdensityEigenVectors%values(mu,nu) + &
                !!                              densityEigenVectors%values(mu,k) *  densityEigenVectors%values(nu,k)*densityEigenValues%values(k) 
                !!    end do
                !!  end do
                !!end do
                !!print *, "atomic density matrix from natural orbitals"
                !!call Matrix_show ( auxdensityEigenVectors)
                write(*,"(A10,A10,A40,F17.12)") "sum of ", trim(speciesName) , "natural orbital occupations", sum(densityEigenValues%values)

                write(*,*) " End of natural orbitals in state: ", state, " for: ", trim(speciesName)
             end do
          end do



          write(*,*) ""
          write(*,*) " END OF NATURAL ORBITALS"
          write(*,*) "=============================="
          write(*,*) ""

       end if
   
      close(unit)

      deallocate ( jj )
      deallocate ( indexConfB )
      deallocate ( indexConfA )
      deallocate ( coupling )
      deallocate ( cilevel )
      deallocate ( cilevelA )
      deallocate ( CIcore_instance%allIndexConf )
      deallocate ( stringAinB )

     deallocate( coefficients, atomicDensityMatrix, ciDensityMatrix )

     !$  timeDB = omp_get_wtime()
     !$  write(*,"(A,F10.4,A4)") "** TOTAL Elapsed Time for Building density matrices: ", timeDB - timeDA ," (s)"

     
  end if
       
    ! print *, i, i, orbital, orbital, CIcore_instance%eigenVectors%values(i,1)**2
    
    ! do mu = 1 , numberOfOrbitals
    !    do nu = 1 , numberOfOrbitals
    
    !       densityMatrix%values(mu,nu) =  &
    !            densityMatrix%values(mu,nu) + &
    !            CIcore_instance%eigenVectors%values(i,state)**2 *&
    !            coefficients%values(mu,orbital)*coefficients%values(nu,orbital)
    !    end do
    ! end do
    
    !!off-Diagonal ground state
    
    ! do mu = 1 , numberOfOrbitals
    !    do nu = 1 , numberOfOrbitals
    
    !       densityMatrix%values(mu,nu) =  &
    !            densityMatrix%values(mu,nu) + &
    !            factor *&
    !            CIcore_instance%eigenVectors%values(i,state) *&
    !            CIcore_instance%eigenVectors%values(j,state) *&
    !            (coefficients%values(mu,orbitalA)*coefficients%values(nu,orbitalB) + coefficients%values(mu,orbitalB)*coefficients%values(nu,orbitalA))
    !    end do
    ! end do
    
    ! call Vector_constructor(eigenValues, numberOfOrbitals)
    ! call Matrix_constructor(eigenVectors, int(numberOfOrbitals,8), int(numberOfOrbitals,8))
    ! call Matrix_eigen(ciOccupationMatrix, eigenValues, eigenVectors, SYMMETRIC)
    
    ! print *, "Diagonal sum", sum(eigenValues%values)
    ! call Vector_show(eigenValues)
    
    ! call Matrix_show(eigenVectors)
    ! print *, arguments(1:2)
    ! call Matrix_show ( densityMatrix )
    
    ! call Matrix_constructor ( ciOccupationNumbers , int(numberOfOrbitals,8) , &
    !      int(CONTROL_instance%CI_STATES_TO_PRINT,8),  0.0_8 )
    
    ! do state=1, CONTROL_instance%CI_STATES_TO_PRINT
    !    sumaPrueba=0
    !    do j=1, numberOfOccupiedOrbitals
    !       ciOccupationNumbers%values(j,state) = 1.0
    !    end do
    
    ! ! !Get occupation numbers from each configuration contribution
    
    !    do i=1, CIcore_instance%numberOfConfigurations
    !       do j=1, numberOfOccupiedOrbitals
    
    !          !! Occupied orbitals
    !          ciOccupationNumbers%values( j, state)= ciOccupationNumbers%values( j, state) -  &
    !               CIcore_instance%eigenVectors%values(i,state)**2
    !          !! Unoccupied orbitals
    !          orbital = CIcore_instance%configurations(i)%occupations(j,specie) 
    
    !          ciOccupationNumbers%values( orbital, state)= ciOccupationNumbers%values( orbital, state) + &
    !               CIcore_instance%eigenVectors%values(i,state)**2
    
    !          ! print *, j, orbital, CIcore_instance%eigenVectors%values(i,state)**2
    !          ! sumaPrueba=sumaPrueba+CIcore_instance%eigenVectors%values(i,state)**2
    !       end do
    !       ! end if
    
    !    end do
    
    !    ! print *, "suma", sumaPrueba
    !    !Build a new density matrix (P) in atomic orbitals
    
    !    call Matrix_constructor ( densityMatrix , &
    !         int(numberOfOrbitals,8), &
    !         int(numberOfOrbitals,8),  0.0_8 )
    
    !    wfnFile = "lowdin.wfn"
    !    wfnUnit = 20
    
    !    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    
    !    arguments(2) = speciesName
    !    arguments(1) = "COEFFICIENTS"
    
    !    coefficients = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
    !         columns= int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))
    
    !    close(wfnUnit)
    
    !    do mu = 1 , numberOfOrbitals
    !       do nu = 1 , numberOfOrbitals
    !          do k = 1 , numberOfOrbitals
    
    !             densityMatrix%values(mu,nu) =  &
    !                  densityMatrix%values(mu,nu) + &
    !                  ciOccupationNumbers%values(k, state)**2* &
    !                  coefficients%values(mu,k)*coefficients%values(nu,k)
    !           end do
    !        end do
    !     end do
    
    !     write(auxstring,*) state
    !     arguments(2) = speciesName
    !     arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
    
    !     call Matrix_writeToFile ( densityMatrix, unit , arguments=arguments(1:2) )
    
    !     print *, arguments(1:2)
    !     call Matrix_show ( densityMatrix )
    
    !     call Matrix_destructor(coefficients)          
    !     call Matrix_destructor(densityMatrix)          
    
    
    !  end do
    
    ! !Write occupation numbers to file
    ! write (6,"(T8,A10,A20)") trim(MolecularSystem_getNameOfSpecie(specie)),"OCCUPATIONS:"
    
    ! call Matrix_show ( ciOccupationNumbers )
    
    ! arguments(2) = speciesName
    ! arguments(1) = "OCCUPATIONS"
    
    ! call Matrix_writeToFile ( ciOccupationNumbers, unit , arguments=arguments(1:2) )
    
    ! call Matrix_destructor(ciOccupationNumbers)          
    
  end subroutine CImod_densityMatrices

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
    integer i,j,m,n,p,q,c
    integer numberOfSpecies
    integer :: isLambdaEqual1

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

    call Matrix_destructor(CIcore_instance%hamiltonianMatrix)
    call Vector_destructorInteger (CIcore_instance%numberOfOccupiedOrbitals)
    call Vector_destructorInteger (CIcore_instance%numberOfOrbitals)
    call Vector_destructor (CIcore_instance%lambda)

    CIcore_instance%isInstanced=.false.

  end subroutine CImod_destructor


end module CImod_


