module CIMCSCF_
  use Matrix_
  use Vector_
  use IndexMap_
  use CIcore_
  use CIdensity_

  implicit none

  public :: &
    CIMCSCF_compute

contains

  subroutine CIMCSCF_compute()
    implicit none

    type(matrix), allocatable :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions, numberOfOccupiedOrbitals
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: i, ii, iiii, ii_aux, k, kk, iikk
    integer :: p,q,r,s, pqrs, pq, rs
    integer :: state
    integer(8) :: numberOfElements
    real(8) :: timeDA, timeDB

    !$ timeDA = omp_get_wtime()

    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) "BUILDING MCSCF MATRICES "
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""
  
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
          !! initializing with HF occupancy, although this is not really neccesary, but RDM subroutines asssume this 
          do i = 1, numberOfOccupiedOrbitals_i
            CI1RDM( spi, state )%values( i, i ) = 1.0_8
          enddo ! i 
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
    call CIMCSCF_energy( CI1RDM, CI2RDM  )

    !call CIMCSCF_gradient()
    !call CIMCSCF_hessian()

    !! Newtown-Rapshon

    !! build unitary transformation

    !! transform coefficients

    !! Transform integrals

    !! Get transformed Integrals

    !! destructor
    do spi = 1, numberOfSpecies
      do state = 1, CONTROL_instance%CI_NUMBER_OF_STATES
        call Matrix_destructor( CI1RDM(spi, state))
      enddo
      do spj = spi, numberOfSpecies
        call Vector_destructor( CI2RDM(spi, spj))
      enddo
    enddo
    deallocate( CI1RDM )
    deallocate( CI2RDM )

    !$  timeDB = omp_get_wtime()
    !$  write(*,"(A,F10.4,A4)") "** TOTAL Elapsed Time for Building MCSCF matrices: ", timeDB - timeDA ," (s)"

  end subroutine CIMCSCF_compute

  subroutine CIMCSCF_energy( CI1RDM, CI2RDM  )
    implicit none
    type(matrix), allocatable :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    integer :: p,q,r,s, pqrs, pq, rs, pq_aux, pqpq
    integer :: pqrs_rdm, pq_rdm, rs_rdm
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions, numberOfOccupiedOrbitals
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    real(8) :: energy_one, energy_two_aa, energy_two_ab, energy_total
    real(8) :: n_pairs

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    energy_one = 0.0_8
    energy_two_aa = 0.0_8
    energy_two_ab = 0.0_8
    energy_total = 0.0_8
    n_pairs = 0.0_8

    !open(unit=1018, file="2rdm", status="replace", form="formatted")

    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          pq = CIcore_instance%twoIndexArray(spi)%values(p,q)
          pq_rdm = ( p - 1) * numberOfContractions_i + q

          energy_one = energy_one + CIcore_instance%twoCenterIntegrals(spi)%values(p,q) * &
                            CI1RDM(spi,1)%values(p,q)

          do r = 1, numberOfContractions_i
            do s = 1, numberOfContractions_i
              rs = CIcore_instance%twoIndexArray(spi)%values(r,s )
              rs_rdm = ( r - 1) * numberOfContractions_i + s

              pqrs = CIcore_instance%fourIndexArray(spi)%values(pq,rs)
              pqrs_rdm = CIcore_two2one ( pq_rdm, rs_rdm )

              energy_two_aa = energy_two_aa + 0.5_8 * CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs,1) * &
                            CI2RDM(spi,spi)%values(pqrs_rdm)
              !write (1018,"(I2,I2,I2,I2,I6,F12.8,F12.8 )" ) p,q,r,s, pqrs_rdm, CI2RDM(spi,spi)%values(pqrs_rdm), CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs,1)

            enddo ! s
          enddo ! r

          pqpq = CIcore_two2one ( pq_rdm, pq_rdm )

          do spj = spi + 1, numberOfSpecies

            pq_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( pq - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do r = 1, numberOfContractions_j
              do s = 1, numberOfContractions_j
                rs = CIcore_instance%twoIndexArray(spj)%values(r,s )
                pqrs = pq_aux + rs

                rs_rdm = ( r - 1) * numberOfContractions_j + s
                pqrs_rdm = ( pq_rdm - 1 ) * numberOfContractions_j * numberOfContractions_j + rs_rdm 
                energy_two_ab = energy_two_ab + CIcore_instance%fourCenterIntegrals(spi,spj)%values(pqrs,1) * &
                            CI2RDM(spi,spj)%values(pqrs_rdm)

              enddo ! s
            enddo ! r


          enddo ! spj

        enddo ! q 
      enddo ! p
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

    !close(1018)

  end subroutine CIMCSCF_energy

  subroutine CIMCSCF_gradient()
    implicit none

  end subroutine CIMCSCF_gradient

  subroutine CIMCSCF_hessian()
    implicit none

  end subroutine CIMCSCF_hessian

end module CIMCSCF_

