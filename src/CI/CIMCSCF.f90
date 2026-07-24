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

      !! initializing with HF occupancy, although this is not really neccesary, but RDM subroutines asssume this 
      !do i = 1, numberOfOccupiedOrbitals_i
      !  ii = CIcore_instance%twoIndexArray(spi)%values(i,i)
      !  iiii = CIcore_instance%fourIndexArray(spi)%values(ii,ii)
      !  CI2RDM( spi, spi )%values( iiii ) = 1.0_8
      !enddo ! i

      do spj = spi + 1, numberOfSpecies
        numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
        numberOfOccupiedOrbitals_j = CIcore_instance%numberOfOccupiedOrbitals%values( spj )

        numberOfContractions = numberOfContractions_i + numberOfContractions_j

        numberOfElements = (numberOfContractions_i*((numberOfContractions_i + 1.0_8)/2.0_8))* &
                           (numberOfContractions_j*((numberOfContractions_j + 1.0_8)/2.0_8))

        call Vector_constructor( CI2RDM(spi,spj), numberOfElements, 0.0_8)

        !!! initializing with HF occupancy, although this is not really neccesary, but RDM subroutines asssume this 
        !do i = 1, numberOfOccupiedOrbitals_i
        !  ii = CIcore_instance%twoIndexArray(spi)%values(i,i)
        !  ii_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( ii - 1_8 )
        !  do k = 1, numberOfOccupiedOrbitals_j
        !    kk = CIcore_instance%twoIndexArray(spj)%values(k,k)
        !    iikk = ii_aux + kk 
        !    CI2RDM( spi, spj )%values( iikk ) = 1.0_8
        !  enddo ! k 
        !enddo ! i 
      enddo ! spj
    enddo ! spi

    call CIdensity_1RDM_SCI( CI1RDM )

    call CIdensity_2RDM_SCI( CI2RDM )

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
              !print *, p,q,r,s, pq_rdm, rs_rdm, pqrs_rdm, CI2RDM(spi,spi)%values(pqrs_rdm), CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs,1)

            enddo ! s
          enddo ! r

          pqpq = CIcore_two2one ( pq_rdm, pq_rdm )
          n_pairs = n_pairs + CI2RDM(spi,spi)%values(pqpq)

          do spj = spi + 1, numberOfSpecies

            pq_aux = CIcore_instance%numberOfSpatialOrbitals2%values( spj ) * ( pq - 1_8 )
            numberOfContractions_j = MolecularSystem_getTotalNumberOfContractions( spj )
            do r = 1, numberOfContractions_j
              do s = 1, numberOfContractions_j
                rs = CIcore_instance%twoIndexArray(spj)%values(r,s )
                pqrs = pq_aux + rs

                energy_two_ab = energy_two_ab + CIcore_instance%fourCenterIntegrals(spi,spj)%values(pqrs,1) * &
                            CI2RDM(spi,spj)%values(pqrs)

              enddo ! s
            enddo ! r

          pqpq = pq_aux + pq
          n_pairs = n_pairs + CI2RDM(spi,spj)%values(pqpq)

          enddo ! spj

        enddo ! q 
      enddo ! p
    enddo ! spi

    energy_total = energy_one + energy_two_aa + energy_two_ab + HartreeFock_instance%puntualInteractionEnergy

    print *, "MCSCF Energy point c",  HartreeFock_instance%puntualInteractionEnergy
    print *, "MCSCF Energy one    ", energy_one
    print *, "MCSCF Energy two aa ", energy_two_aa
    print *, "MCSCF Energy two ab ", energy_two_ab
    print *, "MCSCF Energy total  ", energy_total
    print *, "MCSCF n_pairs       ", n_pairs

  end subroutine CIMCSCF_energy

  subroutine CIMCSCF_gradient()
    implicit none

  end subroutine CIMCSCF_gradient

  subroutine CIMCSCF_hessian()
    implicit none

  end subroutine CIMCSCF_hessian

end module CIMCSCF_

