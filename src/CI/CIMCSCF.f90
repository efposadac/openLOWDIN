module CIMCSCF_
  use Matrix_
  use Vector_
  use IndexMap_
  use CIcore_
  use CIdensity_
  use IntegralsTransformation_, only : IntegralsTransformation_main

  implicit none

  type, public :: MCSCF
    type(vector) :: energy
    type(vector) :: energyChange
    type(matrix) :: maxGradient
    type(matrix) :: totalGradient
    integer :: iter
  endtype MCSCF 

  public :: &
    CIMCSCF_compute, &
    CIMCSCF_show, &
    CIMCSCF_summary

  private :: &
    CIMCSCF_energy, &
    CIMCSCF_gradient, &
    CIMCSCF_hessian, &
    CIMCSCF_newtonRaphson, &
    CIMCSCF_buildUnitaryMatrix, &
    CIMCSCF_rotateCoefficients 

contains
  subroutine CIMCSCF_show () 
    implicit none
    integer(8) :: totalSize
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: numberOfContractions_j
    integer(8) :: sizeAA, sizeBB, sizeAAAA, sizeAABB

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    totalSize = 0_8

    do spi = 1, CIcore_instance%numberOfSpecies 
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
    
      sizeAA = ( numberOfContractions_i * ( numberOfContractions_i + 1 ) ) / 2
      sizeAAAA = ( sizeAA * sizeAA )

      totalSize = totalSize + sizeAA * 8 * ( &
                                              + 1 & ! 1-RDM
                                              + 2 & ! AO2MO one-body, indexmap
                                              + 2 & ! Generalized Fock + Gradient
                                              + 3 & ! Hessian + W_xyxy + W_yxxy
                                              + 4 & ! Rotations, unitary, old, new coeff
                                              + 6 & ! hcore, kinetic, attraction, ext, overlap, density
                                            ) 

      totalSize = totalSize + sizeAAAA * 8 * ( &
                                              + 2 & ! AO2MO two body, indexmap
                                              + 1 & ! 2RDM 
                                             )

      do spj = spi + 1 , CIcore_instance%numberOfSpecies 
        numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
        sizeBB = ( numberOfContractions_j * ( numberOfContractions_j + 1 ) ) / 2
        sizeAABB = ( sizeAA * sizeBB )

        totalSize = totalSize + sizeAABB * 8 * ( &
                                                + 1 & ! AO2MO 
                                                + 1 & ! 2RDM
                                               ) 
      enddo ! spj

      !! SCI conf
      totalSize = totalSize + CISCI_instance%targetSpaceSize_max * &
                                  ( 8 + 2*8 & !! coeff, eigenvectors, W per omp thread
                                    + 1*CIcore_instance%numberOfActiveOrbitals%values(spi) + 4*CIcore_instance%numberOfActiveOrbitals%values(spi) )  !! conf_orb, conf_cc

    enddo ! spi 

    write (6,*) ""
    write (6,*) "-----------------------------------------------------------------------"
    write (6,"(T2,A62)") "       Multi-Configurational Self-Consistent Field (MCSCF)       " 
    write (6,"(T2,A62)") "  Based on J. Chem. Phys. 74, 2384 (1981); doi: 10.1063/1.441359 "
    write (6,"(T2,A62)") "                          J. Charry                              "
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""
    write (6,"(T2,A,F14.3,A3 )") "Estimated memory needed       :", real( totalSize )/(1024**2) , " MB"
    write (6,"(T2,A,F14.3,A3 )") "                               ", real( totalSize )/(1024**3) , " GB"
    write (6,"(T2,A,I8 )")       "Maximum number of iterations  :", CONTROL_instance%CI_MCSCF_MAX_ITER 
    write (6,"(T2,A,F14.3 )")    "Netwon-Raphson damping factor :", CONTROL_instance%CI_MCSCF_DAMPING_FACTOR_NR
    write (6,"(T2,A,F14.3 )")    "Convergence energy criteria   :", 1E-5  
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""
    
  end subroutine CIMCSCF_show

  subroutine CIMCSCF_summary( MCSCF_instance )
    implicit none
    type(MCSCF), intent(in) :: MCSCF_instance
    integer :: k

    write (6,*) ""
    write (6,*) "CONVERGED MCSCF!"
    write (6,*) "____________________________________________________________________"
    write (6,*) "| Iter |    Energy |    Δ Energy |  Max. Gradient |  Tot. Gradient |"
    write (6,*) "____________________________________________________________________"
    do k = 1, MCSCF_instance%iter 
      write (6,"(T2,I6,F14.6,ES14.3,A3,ES14.3,A3,ES14.3)") k, MCSCF_instance%energy%values(k), MCSCF_instance%energyChange%values(k), &
                  "  ", maxval(MCSCF_instance%maxGradient%values(k,:)), "   ", sum(MCSCF_instance%totalGradient%values(k,:))
    enddo
    write (6,*) "____________________________________________________________________"
    write (6,"(A,F14.8)")  "FINAL MCSCF Energy        =", MCSCF_instance%energy%values( MCSCF_instance%iter )
    write (6,"(A,ES14.3)") "FINAL MCSCF Delta Energy  =", MCSCF_instance%energyChange%values( MCSCF_instance%iter )
    write (6,"(A,ES14.3)") "FINAL MCSCF Max. Gradient =", maxval(MCSCF_instance%maxGradient%values( MCSCF_instance%iter, :))
    write (6,"(A,ES14.3)") "FINAL MCSCF Tot. Gradient =", sum(MCSCF_instance%totalGradient%values( MCSCF_instance%iter, :)) 

    write (6,*) ""
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) "          END MCSCF CALCULATION"
    write (6,*) "-----------------------------------------------------------------------"
    write (6,*) ""

  end subroutine CIMCSCF_summary 

  subroutine CIMCSCF_compute( MCSCF_instance )
    implicit none
    type(MCSCF), intent(inout) :: MCSCF_instance
    type(matrix), allocatable :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable :: gradient(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable :: rotations(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: state
    integer(8) :: numberOfElements
    real(8) :: timeDA, timeDB

    !$ timeDA = omp_get_wtime()

    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) "BUILDING MCSCF MATRICES "
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""

    MCSCF_instance%iter = MCSCF_instance%iter + 1

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    !! constructor

    !! matrix allocation 1RDM
    allocate( CI1RDM(numberOfSpecies, CONTROL_instance%CI_NUMBER_OF_STATES) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
       numberOfOccupiedOrbitals_i = CIcore_instance%numberOfOccupiedOrbitals%values( spi )

       do state = 1, CONTROL_instance%CI_NUMBER_OF_STATES
          call Matrix_constructor ( CI1RDM( spi, state) , &
               int( numberOfContractions_i, 8), &
               int( numberOfContractions_i, 8), 0.0_8 )
       enddo ! state
    enddo ! spi

    !! matrix allocation 2RDM
    allocate( CI2RDM( numberOfSpecies, numberOfSpecies ) )

    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      numberOfOccupiedOrbitals_i = CIcore_instance%numberOfOccupiedOrbitals%values( spi )

      !! not the same as integrals!
      numberOfElements = int( ( (numberOfContractions_i * numberOfContractions_i )* &
                                (numberOfContractions_i * numberOfContractions_i + 1) / 2.0 ), 8)

      call Vector_constructor(CI2RDM(spi,spi), numberOfElements, 0.0_8)
      CI2RDM(spi,spi)%values = 0.0_8

      do spj = spi + 1, numberOfSpecies
        numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
        numberOfOccupiedOrbitals_j = CIcore_instance%numberOfOccupiedOrbitals%values( spj )

        numberOfContractions = numberOfContractions_i + numberOfContractions_j

        numberOfElements = int( (numberOfContractions_i * numberOfContractions_i )* &
                                (numberOfContractions_j * numberOfContractions_j ), 8)


        call Vector_constructor( CI2RDM(spi,spj), numberOfElements, 0.0_8)

      enddo ! spj
    enddo ! spi

    write (6,*) "Building 1-RDM ..."
    call CIdensity_1RDM_SCI( CI1RDM )

    write (6,*) "Building 2-RDM ..."
    call CIdensity_2RDM_SCI( CI2RDM )

    write (6,*) "Computing the MCSCF energy from 1- and 2-RDM ..."
    call CIMCSCF_energy( MCSCF_instance, CI1RDM, CI2RDM  )

    write (6,*) "Building MCSCF gradient ..."
    call CIMCSCF_gradient( MCSCF_instance, CI1RDM, CI2RDM, fock, gradient )

    write (6,*) "Building MCSCF hessian (diagonal) ..."
    call CIMCSCF_hessian( CI1RDM, CI2RDM, fock, hessian, gradient )

    !! Newtown-Rapshon
    write (6,*) "Performing MCSCF Newton Raphson minimization ..."
    call CIMCSCF_newtonRaphson ( gradient, hessian, rotations )

    !! build unitary transformation
    write (6,*) "Building MCSCF unitary transformation matrix U ..."
    call CIMCSCF_buildUnitaryMatrix ( rotations, UnitaryMatrix ) 

    !! transform coefficients and save them in a file
    call CIMCSCF_rotateCoefficients ( unitaryMatrix, CI1RDM )

    !! Transform integrals
    call IntegralsTransformation_main()

    !! destructor
    do spi = 1, numberOfSpecies
      do state = 1, CONTROL_instance%CI_NUMBER_OF_STATES
        call Matrix_destructor( CI1RDM(spi, state))
      enddo
      do spj = spi, numberOfSpecies
        call Vector_destructor( CI2RDM(spi, spj))
      enddo
      call Matrix_destructor( rotations(spi) )
      call Matrix_destructor( hessian(spi) )
      call Matrix_destructor( gradient(spi) )
      call Matrix_destructor( fock(spi) )
    enddo

    deallocate( rotations )
    deallocate( hessian )
    deallocate( gradient )
    deallocate( fock )
    deallocate( CI2RDM )
    deallocate( CI1RDM )

    !$  timeDB = omp_get_wtime()
    !$  write(*,"(A,F10.4,A4)") "** TOTAL Elapsed Time for Building MCSCF matrices: ", timeDB - timeDA ," (s)"

  end subroutine CIMCSCF_compute

  subroutine CIMCSCF_energy( MCSCF_instance, CI1RDM, CI2RDM  )
    implicit none
    type(MCSCF), intent(inout) :: MCSCF_instance
    type(matrix), allocatable, intent(inout) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(inout) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    integer :: p,q,r,s, pqrs, pq, rs, pq_aux
    integer :: pqrs_rdm, pq_rdm, rs_rdm
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: numberOfContractions_j
    real(8) :: energy_one, energy_two_aa, energy_two_ab, energy_total
    real(8) :: n_pairs

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    energy_one = 0.0_8
    energy_two_aa = 0.0_8
    energy_two_ab = 0.0_8
    energy_total = 0.0_8
    n_pairs = 0.0_8

    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      !! one body
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          pq = CIcore_instance%twoIndexArray(spi)%values(p,q)
          pq_rdm = p + ( q - 1) * numberOfContractions_i 

          energy_one = energy_one + CIcore_instance%twoCenterIntegrals(spi)%values(p,q) * &
                            CI1RDM(spi,1)%values(p,q)
                            
          !! two body intra
          do r = 1, numberOfContractions_i
            do s = 1, numberOfContractions_i
              rs = CIcore_instance%twoIndexArray(spi)%values(r,s )
              rs_rdm = r + ( s - 1) * numberOfContractions_i 

              pqrs = CIcore_instance%fourIndexArray(spi)%values(pq,rs)
              pqrs_rdm = CIcore_two2one ( pq_rdm, rs_rdm )

              energy_two_aa = energy_two_aa + 0.5_8 * CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs,1) * &
                            CI2RDM(spi,spi)%values(pqrs_rdm)
            enddo ! s
          enddo ! r
        enddo ! q 
      enddo ! p
    enddo ! spi

    !! two body inter
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do spj = spi + 1, numberOfSpecies
      numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
        do p = 1, numberOfContractions_i
          do q = 1, numberOfContractions_i
            pq = CIcore_instance%twoIndexArray(spi)%values(p,q)
            pq_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( pq - 1_8 )

            pq_rdm = p + ( q - 1) * numberOfContractions_i 
            do r = 1, numberOfContractions_j
              do s = 1, numberOfContractions_j
                rs = CIcore_instance%twoIndexArray(spj)%values(r,s )
                pqrs = pq_aux + rs

                rs_rdm = r + ( s - 1) * numberOfContractions_j 
                pqrs_rdm = pq_rdm + ( rs_rdm - 1 ) * numberOfContractions_i * numberOfContractions_i
                energy_two_ab = energy_two_ab + CIcore_instance%fourCenterIntegrals(spi,spj)%values(pqrs,1) * &
                            CI2RDM(spi,spj)%values(pqrs_rdm)
              enddo ! s
            enddo ! r
          enddo ! q 
        enddo ! p

      enddo ! spj
    enddo ! spi

    energy_total = energy_one + energy_two_aa + energy_two_ab + HartreeFock_instance%puntualInteractionEnergy
 
    write (6,*) ""
    write (6,*) "MCSCF initial energy components: "
    write (6,"(T2,A37,F25.12)") "MCSCF Initial fixed potential =       ", HartreeFock_instance%puntualInteractionEnergy
    write (6,"(T2,A37,F25.12)") "MCSCF Initial one-body energy =       ", energy_one
    write (6,"(T2,A37,F25.12)") "MCSCF Initial two-body intra energy = ", energy_two_aa
    write (6,"(T2,A37,F25.12)") "MCSCF Initial two-body inter energy = ", energy_two_ab
    write (6,"(T2,A37,F25.12)") "MCSCF Initial total energy =          ", energy_total
    write (6,*) ""
    MCSCF_instance%energy%values(MCSCF_instance%iter) = energy_total 
    MCSCF_instance%energyChange%values(MCSCF_instance%iter) = MCSCF_instance%energy%values(MCSCF_instance%iter) - &
                                                                MCSCF_instance%energy%values(MCSCF_instance%iter - 1) 

  end subroutine CIMCSCF_energy

  subroutine CIMCSCF_gradient( MCSCF_instance, CI1RDM, CI2RDM, fock, gradient )
    implicit none
    type(MCSCF), intent(inout) :: MCSCF_instance
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(in) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: gradient(:) ! species % numcontractions, numcontractions
    integer :: p,q,r,s,t, pr, pr_rdm, qr, qr_aux, st, st_rdm, qrst, prst_rdm
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: numberOfContractions_j

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( fock(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( fock( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    allocate( gradient(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( gradient( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    !! generalized Fock matrix
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )

      !! one body 
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i

          do r = 1, numberOfContractions_i

            fock(spi)%values( p, q ) = fock(spi)%values( p, q ) + CIcore_instance%twoCenterIntegrals(spi)%values(q,r) * &
                                                                      CI1RDM(spi,1)%values(p,r)

          enddo
        enddo
      enddo

      !! two body intra
      !$omp parallel &
      !$omp& private ( p, q, r, pr, pr_rdm, qr, s, t, st, st_rdm, qrst, prst_rdm, & 
      !$omp            spj, numberOfContractions_j, qr_aux ) &  
      !$omp& shared ( fock ) 
      !$omp do schedule ( dynamic ) collapse(2)
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i

          do r = 1, numberOfContractions_i

            pr = CIcore_instance%twoIndexArray(spi)%values(p,r)
            pr_rdm = p + ( r - 1) * numberOfContractions_i 

            qr = CIcore_instance%twoIndexArray(spi)%values( q, r )

           do s = 1, numberOfContractions_i
            do t = 1, numberOfContractions_i

                st = CIcore_instance%twoIndexArray(spi)%values(s,t)
                st_rdm = s + ( t - 1) * numberOfContractions_i 

                qrst = CIcore_instance%fourIndexArray(spi)%values(qr,st)
                prst_rdm = CIcore_two2one ( pr_rdm, st_rdm )

                fock(spi)%values( p, q ) = fock(spi)%values( p, q ) + CIcore_instance%fourCenterIntegrals(spi,spi)%values(qrst,1) * &
                                                                      CI2RDM(spi,spi)%values(prst_rdm)
              enddo ! t
            enddo ! s
          enddo ! r

          !! two body inter BA
          do spj = 1, spi - 1

            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )

            do r = 1, numberOfContractions_i
              qr = CIcore_instance%twoIndexArray(spi)%values( q, r )
              pr_rdm = p + ( r - 1) * numberOfContractions_i 

              qr_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( qr - 1_8 )
              do s = 1, numberOfContractions_j
                do t = 1, numberOfContractions_j
                  st = CIcore_instance%twoIndexArray(spj)%values( s, t )
                  qrst = qr_aux + st

                  st_rdm = s + ( t - 1) * numberOfContractions_j 
                  prst_rdm = st_rdm + ( pr_rdm - 1 ) * numberOfContractions_j * numberOfContractions_j  

                  fock(spi)%values( p, q ) = fock(spi)%values( p, q ) + CIcore_instance%fourCenterIntegrals(spi,spj)%values(qrst,1) * &
                                                                      CI2RDM(spj,spi)%values(prst_rdm)

                enddo ! t
              enddo ! s
            enddo ! r
          enddo ! spj

          !! two body inter AB
          do spj = spi + 1, numberOfSpecies
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )

            do r = 1, numberOfContractions_i

              qr = CIcore_instance%twoIndexArray(spi)%values( q, r )
              pr_rdm = p + ( r - 1) * numberOfContractions_i 

              qr_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( qr - 1_8 )

              do s = 1, numberOfContractions_j
                do t = 1, numberOfContractions_j
                  st = CIcore_instance%twoIndexArray(spj)%values( s, t )
                  qrst = qr_aux + st

                  st_rdm = s + ( t - 1) * numberOfContractions_j 
                  prst_rdm = pr_rdm + ( st_rdm - 1 ) * numberOfContractions_i * numberOfContractions_i

                  fock(spi)%values( p, q ) = fock(spi)%values( p, q ) + &
                                              CIcore_instance%fourCenterIntegrals(spi,spj)%values(qrst,1) * &
                                              CI2RDM(spi,spj)%values(prst_rdm)

                enddo ! t
              enddo ! s
            enddo ! r
          enddo ! spj
        enddo ! q 
      enddo ! p
      !$omp enddo
      !$omp end parallel

    enddo ! spi

    !do spi = 1, numberOfSpecies
    !  print *, "fock for spi", spi
    !  call Matrix_show (fock(spi))
    !enddo ! spi

    !! gradient 
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          gradient(spi)%values(p,q) = gradient(spi)%values(p,q) + 2.0_8 * ( fock(spi)%values(p,q) - fock(spi)%values(q,p) )
        enddo ! q 
      enddo ! p

      MCSCF_instance%maxGradient%values(MCSCF_instance%iter,spi) = maxval(gradient(spi)%values)
      MCSCF_instance%totalGradient%values(MCSCF_instance%iter,spi) = sum(abs(gradient(spi)%values))
      !call Matrix_show (gradient(spi))

    enddo ! spi

  end subroutine CIMCSCF_gradient

  subroutine CIMCSCF_hessian( CI1RDM, CI2RDM, fock, hessian, gradient )
    implicit none
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(in) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable, intent(in) :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable, intent(in) :: gradient(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable :: W_xyxy(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable :: W_yxxy(:) ! species % numcontractions, numcontractions (diagonal)
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: numberOfContractions_j 
    integer :: x, y, m, n
    integer :: yy, yy_aux, mn, ym, yn, yymn, ymyn
    integer :: xx_rdm, mn_rdm, xn_rdm, mx_rdm, xm_rdm, xxmn_rdm, xnmx_rdm, xnxm_rdm
    integer :: xy, xy_aux, xm, xymn, xmyn
    integer :: xy_rdm, my_rdm, ym_rdm, xymn_rdm, xnmy_rdm, xnym_rdm
    integer :: p, q

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( hessian(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( hessian( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    allocate( W_xyxy(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( W_xyxy( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo
    allocate( W_yxxy(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( W_yxxy( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    ! W_xyxy
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do x = 1, numberOfContractions_i

        xx_rdm = ( x - 1) * numberOfContractions_i + x
        do y = 1, numberOfContractions_i

          W_xyxy(spi)%values(x,y) = W_xyxy(spi)%values(x,y) + 2.0_8 * CI1RDM(spi,1)%values(x,x) * &
                                                              CIcore_instance%twoCenterIntegrals(spi)%values(y,y)
                                                                
          yy = CIcore_instance%twoIndexArray(spi)%values(y,y)
                                                                     
          do m = 1, numberOfContractions_i

            ym = CIcore_instance%twoIndexArray(spi)%values(y,m)
            mx_rdm = m + ( x - 1) * numberOfContractions_i 
            xm_rdm = x + ( m - 1) * numberOfContractions_i 

            do n = 1, numberOfContractions_i

              mn_rdm = m + ( n - 1) * numberOfContractions_i 

              xxmn_rdm = CIcore_two2one ( xx_rdm, mn_rdm )

              mn = CIcore_instance%twoIndexArray(spi)%values(m,n)
              yymn = CIcore_instance%fourIndexArray(spi)%values(yy,mn)

              W_xyxy(spi)%values( x, y ) = W_xyxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spi,spi)%values(xxmn_rdm) * &
                                                            CIcore_instance%fourCenterIntegrals(spi,spi)%values(yymn,1)

              yn = CIcore_instance%twoIndexArray(spi)%values(y,n)
              ymyn = CIcore_instance%fourIndexArray(spi)%values(ym,yn)

              xn_rdm = x + ( n - 1) * numberOfContractions_i
              xnmx_rdm = CIcore_two2one ( xn_rdm, mx_rdm )
              xnxm_rdm = CIcore_two2one ( xn_rdm, xm_rdm )

              W_xyxy(spi)%values( x, y ) = W_xyxy(spi)%values( x, y ) + 2.0_8 * &
                                              ( CI2RDM(spi,spi)%values(xnmx_rdm) + CI2RDM(spi,spi)%values(xnxm_rdm) ) * &
                                              CIcore_instance%fourCenterIntegrals(spi,spi)%values(ymyn,1)
            enddo ! n
          enddo ! m

          do spj = 1, spi - 1

            yy_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( yy - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do m = 1, numberOfContractions_j
              do n = 1, numberOfContractions_j

                mn = CIcore_instance%twoIndexArray(spj)%values( m, n )
                yymn = yy_aux + mn

                mn_rdm = m + ( n - 1) * numberOfContractions_j 
                xxmn_rdm = mn_rdm + ( xx_rdm - 1 ) * numberOfContractions_j * numberOfContractions_j 
                W_xyxy(spi)%values( x, y ) = W_xyxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spj,spi)%values(xxmn_rdm) * &
                                                      CIcore_instance%fourCenterIntegrals(spi,spj)%values(yymn,1) 
                                                                   
              enddo ! n 
            enddo ! m
          enddo ! spj 

          do spj = spi + 1, numberOfSpecies

            yy_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( yy - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do m = 1, numberOfContractions_j
              do n = 1, numberOfContractions_j

                mn = CIcore_instance%twoIndexArray(spj)%values( m, n )
                yymn = yy_aux + mn

                mn_rdm = m + ( n - 1) * numberOfContractions_j 
                xxmn_rdm = xx_rdm + ( mn_rdm - 1 ) * numberOfContractions_i * numberOfContractions_i  
                W_xyxy(spi)%values( x, y ) = W_xyxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spi,spj)%values(xxmn_rdm) * &
                                                  CIcore_instance%fourCenterIntegrals(spi,spj)%values(yymn,1)
              enddo ! n 
            enddo ! m
          enddo ! spj 

        enddo ! y

        W_xyxy(spi)%values( x, x ) = W_xyxy(spi)%values( x, x ) + fock(spi)%values( x, x ) + fock(spi)%values( x, x )
      enddo ! x

    enddo ! spi

    ! W_yxxy
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do x = 1, numberOfContractions_i

        do y = 1, numberOfContractions_i

          W_yxxy(spi)%values(x,y) = W_yxxy(spi)%values(x,y) + 2.0_8 * CI1RDM(spi,1)%values(y,x) * &
                                                            CIcore_instance%twoCenterIntegrals(spi)%values(x,y) 
                                                                 
          xy_rdm = x + ( y - 1) * numberOfContractions_i 
          xy = CIcore_instance%twoIndexArray(spi)%values(x,y)
                                                                     
          do m = 1, numberOfContractions_i

            my_rdm = m + ( y - 1) * numberOfContractions_i 
            ym_rdm = y + ( m - 1) * numberOfContractions_i

            do n = 1, numberOfContractions_i

              mn = CIcore_instance%twoIndexArray(spi)%values(m,n)
              mn_rdm = m + ( n - 1) * numberOfContractions_i 
              xymn_rdm = CIcore_two2one ( xy_rdm, mn_rdm )

              mn = CIcore_instance%twoIndexArray(spi)%values(m,n)
              xymn = CIcore_instance%fourIndexArray(spi)%values(xy,mn)

              W_yxxy(spi)%values( x, y ) = W_yxxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spi,spi)%values(xymn_rdm) * &
                                             CIcore_instance%fourCenterIntegrals(spi,spi)%values(xymn,1)
                                                            

              xm = CIcore_instance%twoIndexArray(spi)%values(x,m)
              yn = CIcore_instance%twoIndexArray(spi)%values(y,n)
              xmyn = CIcore_instance%fourIndexArray(spi)%values(xm,yn)

              xn_rdm = x + ( n - 1) * numberOfContractions_i
              xnmy_rdm = CIcore_two2one ( xn_rdm, my_rdm )
              xnym_rdm = CIcore_two2one ( xn_rdm, ym_rdm )

              W_yxxy(spi)%values( x, y ) = W_yxxy(spi)%values( x, y ) + 2.0_8 * &
                                           ( CI2RDM(spi,spi)%values(xnmy_rdm) + CI2RDM(spi,spi)%values(xnym_rdm) ) * &
                                             CIcore_instance%fourCenterIntegrals(spi,spi)%values(xmyn,1) 
            enddo ! n
          enddo ! m

          do spj = 1, spi - 1

            xy_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( xy - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do m = 1, numberOfContractions_j
              do n = 1, numberOfContractions_j

                mn = CIcore_instance%twoIndexArray(spj)%values( m, n )
                xymn = xy_aux + mn

                mn_rdm = m + ( n - 1) * numberOfContractions_j
                xymn_rdm = mn_rdm + ( xy_rdm - 1 ) * numberOfContractions_j * numberOfContractions_j
                W_yxxy(spi)%values( x, y ) = W_yxxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spj,spi)%values(xymn_rdm) * &
                                              CIcore_instance%fourCenterIntegrals(spi,spj)%values(xymn,1) 
              enddo ! n 
            enddo ! m
          enddo ! spj 

          do spj = spi + 1, numberOfSpecies

            xy_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( xy - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do m = 1, numberOfContractions_j
              do n = 1, numberOfContractions_j
                mn = CIcore_instance%twoIndexArray(spj)%values( m, n )
                xymn = xy_aux + mn
                mn_rdm = m + ( n - 1) * numberOfContractions_j 
                xymn_rdm = xy_rdm + ( mn_rdm - 1 ) * numberOfContractions_i * numberOfContractions_i 
                W_yxxy(spi)%values( x, y ) = W_yxxy(spi)%values( x, y ) + 2.0_8 * CI2RDM(spi,spj)%values(xymn_rdm) * &
                                            CIcore_instance%fourCenterIntegrals(spi,spj)%values(xymn,1) 
              enddo ! n 
            enddo ! m
          enddo ! spj 

        enddo ! y

        do y = 1, numberOfContractions_i
          W_yxxy(spi)%values( x, y ) = W_yxxy(spi)%values( x, y ) + 2.0_8 * fock(spi)%values( y, y )
        enddo ! y
      enddo ! x
    enddo ! spi

    !! building the hessian
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )


      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          hessian(spi)%values( p, q ) = hessian(spi)%values( p, q ) + W_xyxy(spi)%values( p, q ) &
                                                                    - W_yxxy(spi)%values( p, q ) &
                                                                    - W_yxxy(spi)%values( q, p ) &
                                                                    + W_xyxy(spi)%values( q, p ) 
        enddo ! q
      enddo ! p 
    enddo ! spi

    !do spi = 1, numberOfSpecies
    !  print *, "xyxy"
    !  call Matrix_show(W_xyxy(spi))
    !  print *, "yxxy"
    !  call Matrix_show(W_yxxy(spi))
    !  print *, "xyyx"
    !  call Matrix_show(W_xyyx(spi))
    !  print *, "yxyx"
    !  call Matrix_show(W_yxyx(spi))
    !enddo ! spi

    deallocate( W_yxxy )
    deallocate( W_xyxy )

    !do spi = 1, numberOfSpecies
    !  print *, "hessian for spi", spi
    !  call Matrix_show (hessian(spi))
    !enddo ! spi

  end subroutine CIMCSCF_hessian

  subroutine CIMCSCF_newtonRaphson( gradient,  hessian, rotations )
    implicit none
    type(matrix), allocatable, intent(in) :: gradient(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(in) :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable, intent(inout) :: rotations(:) ! species % numcontractions, numcontractions 
    integer :: spi, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: p, q
    real(8) :: epsilon
    real(8) :: update
    real(8) :: dampingFactor

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    dampingFactor = CONTROL_instance%CI_MCSCF_DAMPING_FACTOR_NR
    epsilon = 1.0E-6

    allocate( rotations(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( rotations( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    !! rotations = - gradient / hessian
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )

      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          update = - gradient(spi)%values(p,q) * dampingFactor / ( hessian(spi)%values(p,q) + epsilon) 
          !if ( abs(update) > 0.50_8 ) update = sign(0.50_8, update ) 
          if ( abs(update) > 0.50_8 ) write (*,*) "Warning! large rotations detected"
          rotations(spi)%values(p,q) = update
        enddo ! q
      enddo ! p 
    enddo ! spi

    !! antisymmetrizing the orbital rotations
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do p = 1, numberOfContractions_i
        do q = p + 1, numberOfContractions_i
          rotations(spi)%values(q,p) = - rotations(spi)%values(p,q)
        enddo ! q
        rotations(spi)%values(p,p) = 0.0_8
      enddo ! p
    enddo ! spi

    !do spi = 1, numberOfSpecies
    !  print *, "rotations for spi", spi
    !  call Matrix_show (rotations(spi))
    !enddo ! spi

  end subroutine CIMCSCF_newtonRaphson 

  subroutine CIMCSCF_buildUnitaryMatrix ( rotations, UnitaryMatrix ) 
   implicit none
    type(matrix), allocatable, intent(inout) :: rotations(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable, intent(inout) :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
    integer :: spi, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: p, q, r
    real(8), allocatable :: A(:,:)  
    real(8), allocatable :: B(:,:)  
    real(8), allocatable :: work(:)  
    integer, allocatable :: IPIV(:)
    integer :: info
    integer :: lwork
    real(8) :: auxvalue

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( unitaryMatrix(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( unitaryMatrix( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    !! Caley transformation ( exponential of a matrix)
    do spi = 1, numberOfSpecies

      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      lwork = numberOfContractions_i 

      !! building A = I - rotations and B = I + rotations. Keep B in unitary matrix U
      allocate ( A ( numberOfContractions_i, numberOfContractions_i ) )
      allocate ( B ( numberOfContractions_i, numberOfContractions_i ) )
      allocate ( work ( lwork ) )
      allocate ( IPIV ( lwork ) )

      A = 0.0_8
      B = 0.0_8
      work = 0.0_8
      IPIV = 0

      do p = 1, numberOfContractions_i 
        A(p, p ) = 1.0_8
        B(p, p ) = 1.0_8
      enddo ! p 

      do q = 1, numberOfContractions_i
        do p = 1, numberOfContractions_i 
          A(p, q) = A(p, q) - 1.0_8 * rotations(spi)%values(p, q)
          B(p, q) = B(p, q) + 1.0_8 * rotations(spi)%values(p, q)
        enddo ! p 
      enddo ! q

      !! LAPACK LU Factorization of A
      !! DGETRF computes: A = P * L * U
      call dgetrf( numberOfContractions_i, &
                   numberOfContractions_i, &
                   B, &
                   numberOfContractions_i, &
                   IPIV, &
                   info)

      if (info /= 0) then
          write(*,*) 'Error in DGETRF during Cayley Transform! INFO = ', info
      end if

      lwork = -1
      call dgetri(numberOfContractions_i, &     ! N
                  B, &                          ! A
                  numberOfContractions_i, &     ! LDA
                  IPIV, &                       ! IPIV
                  work, &                       ! work 
                  lwork, &
                  info)

      lwork = work(1)
      deallocate ( work )
      allocate ( work ( lwork ) )
      work = 0.0_8

      call dgetri(numberOfContractions_i, &     ! N
                  B                     , &     ! A
                  numberOfContractions_i, &     ! LDA
                  IPIV, &                       ! IPIV
                  work, &                       ! work 
                  lwork, &
                  info)

      !!U = A * B 
      do q = 1, numberOfContractions_i
        do p = 1, numberOfContractions_i 
          auxvalue = 0.0_8
          do r = 1, numberOfContractions_i 
            auxvalue = auxvalue + A(p,r) * B(r,q) 
          enddo
          unitaryMatrix(spi)%values(p,q) = auxvalue
        enddo
      enddo

      if (info /= 0) then
          write(*,*) 'Error in DGETRS during Cayley Transform! INFO = ', info
      end if

      deallocate ( work )
      deallocate ( IPIV )
      deallocate ( A )
      deallocate ( B )

      !print *, "U"
      !call Matrix_show ( unitaryMatrix(spi) )

    enddo ! spi

  end subroutine CIMCSCF_buildUnitaryMatrix 

  subroutine CIMCSCF_rotateCoefficients ( UnitaryMatrix, CI1RDM ) 
   implicit none
    type(matrix), allocatable, intent(in) :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    integer :: spi, numberOfSpecies
    integer :: numberOfContractions_i
    integer :: p, q, r
    type(matrix), allocatable :: coefficients_old(:)
    type(matrix), allocatable :: coefficients_new(:)
    integer :: wfnunit
    character(50) :: speciesName
    character(50) :: wfnfile
    character(100) :: arguments(2)
    type(Matrix), allocatable :: hcoreMatrix(:)
    type(Matrix), allocatable :: densityMatrix(:), overlapMatrix(:)
    type(matrix), allocatable :: kineticMatrix(:), attractionMatrix(:), externalPotMatrix(:)
    real(8) :: auxvalue

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( coefficients_old(numberOfSpecies) )
    allocate( coefficients_new(numberOfSpecies) )
    allocate( hcoreMatrix(numberOfSpecies) )
    allocate( kineticMatrix(numberOfSpecies) )
    allocate( attractionMatrix(numberOfSpecies) )
    allocate( externalPotMatrix(numberOfSpecies) )
    allocate( densityMatrix(numberOfSpecies) )
    allocate( overlapMatrix(numberOfSpecies) )

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! open to get old coefficient Matrix 
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! matrix allocation
    do spi = 1, numberOfSpecies

      speciesName = MolecularSystem_getNameOfSpecies( spi )
       
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )

      arguments(2) = speciesName

      !! read hcore to save it again in the new file
      arguments(1) = "HCORE"
      hcoreMatrix( spi ) = Matrix_getFromFile( unit = wfnUnit, rows = int(numberOfContractions_i, 8), &
                                                 columns = int(numberOfContractions_i, 8), &
                                                 binary = .true., arguments = arguments(1:2) )

      arguments(1) = "COEFFICIENTS"
      coefficients_old( spi ) = Matrix_getFromFile( unit = wfnUnit, rows = int(numberOfContractions_i, 8), &
                                                 columns = int(numberOfContractions_i, 8), &
                                                 binary = .true., arguments = arguments(1:2) )

      arguments(1) = "KINETIC"
      kineticMatrix( spi ) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions_i, 8), &
                                                  columns=int(numberOfContractions_i, 8), binary=.true., arguments=arguments(1:2))
       
      arguments(1) = "ATTRACTION"
      attractionMatrix( spi ) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions_i, 8), &
                                                     columns=int(numberOfContractions_i, 8), binary=.true., arguments=arguments(1:2))
      arguments(1) = "EXTERNAL-POTENTIAL"
      if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
      externalPotMatrix( spi ) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions_i, 8), &
                                                       columns=int(numberOfContractions_i, 8), binary=.true., arguments=arguments(1:2))

      arguments(1) = "OVERLAP"
      overlapMatrix( spi ) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions_i, 8), &
                                                     columns=int(numberOfContractions_i, 8), binary=.true., arguments=arguments(1:2))

      call Matrix_constructor( coefficients_new(spi), int(numberOfContractions_i, 8), int(numberOfContractions_i, 8), 0.0_8 )

      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i 
          auxvalue = 0.0_8
          do r = 1, numberOfContractions_i 
            auxvalue = auxvalue + coefficients_old(spi)%values(p,r) * unitaryMatrix(spi)%values(r,q) 
          enddo
          coefficients_new(spi)%values(p,q) = auxvalue
        enddo
      enddo
      call Matrix_constructor( densityMatrix(spi), int(numberOfContractions_i, 8), int(numberOfContractions_i, 8), 0.0_8 )

      !print *, " new coefficients ", spi
      !call Matrix_show ( coefficients_new(spi) )

      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i 
          auxvalue = 0.0_8
          do r = 1, CIcore_instance%numberOfOccupiedOrbitals%values( spi )
            auxvalue = auxvalue + coefficients_new(spi)%values(p,r) * coefficients_new(spi)%values(q,r)
          enddo
          densityMatrix(spi)%values(p,q) = auxvalue
        enddo
      enddo

    enddo ! spi

    close (wfnUnit)

    !! open again to save new coefficient Matrix 
    open(unit=wfnUnit, file=trim(wfnFile), status="replace", form="unformatted")

    do spi = 1, numberOfSpecies

      speciesName = MolecularSystem_getNameOfSpecies( spi )
       
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      arguments(2) = speciesName

      arguments(1) = "COEFFICIENTS"
      call Matrix_writeToFile ( coefficients_new(spi), wfnUnit, arguments=arguments(1:2), binary=.true. )

      arguments(1) = "DENSITY"
      call Matrix_writeToFile ( CI1RDM(spi,1), wfnUnit , arguments=arguments(1:2), binary=.true. )

      arguments(1) = "HCORE"
      call Matrix_writeToFile ( hcoreMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )

      arguments(1) = "KINETIC"
      call Matrix_writeToFile ( kineticMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )
       
      arguments(1) = "ATTRACTION"
      call Matrix_writeToFile ( attractionMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )

      arguments(1) = "OVERLAP"
      call Matrix_writeToFile ( overlapMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )

      arguments(1) = "DENSITY"
      call Matrix_writeToFile ( densityMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )

      arguments(1) = "EXTERNAL-POTENTIAL"
      if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
      call Matrix_writeToFile ( externalPotMatrix(spi), wfnUnit , arguments=arguments(1:2), binary=.true. )

      call Matrix_destructor ( externalPotMatrix(spi) ) 
      call Matrix_destructor ( attractionMatrix(spi) ) 
      call Matrix_destructor ( kineticMatrix(spi) ) 
      call Matrix_destructor ( hcoreMatrix(spi) ) 

    enddo ! spi

    close (wfnUnit)
    deallocate ( externalPotMatrix ) 
    deallocate ( attractionMatrix ) 
    deallocate ( kineticMatrix ) 
    deallocate ( hcoreMatrix ) 
    deallocate ( coefficients_old ) 
    deallocate ( coefficients_new ) 
    deallocate ( densityMatrix )
    deallocate ( overlapMatrix )

  end subroutine CIMCSCF_rotateCoefficients

end module CIMCSCF_

