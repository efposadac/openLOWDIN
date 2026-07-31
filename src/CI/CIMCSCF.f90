module CIMCSCF_
  use Matrix_
  use Vector_
  use IndexMap_
  use CIcore_
  use CIdensity_
  use IntegralsTransformation_, only : IntegralsTransformation_main

  implicit none

  public :: &
    CIMCSCF_compute

contains

  subroutine CIMCSCF_compute()
    implicit none

    type(matrix), allocatable :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable :: gradient(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable :: rotations(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
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

    !do spi = 1, numberOfSpecies
    !  print *, "1rdm for", spi
    !  call Matrix_show (CI1RDM(spi,1))
    !enddo 

    write (6,*) "Building 2-RDM ..."
    call CIdensity_2RDM_SCI( CI2RDM )

    !CI2RDM(1,3)%values(:) = CI2RDM(1,2)%values(:)
    !CI2RDM(2,3)%values(:) = CI2RDM(1,2)%values(:)

    write (6,*) "Computing the MCSCF energy from 1- and 2-RDM ..."
    call CIMCSCF_energy( CI1RDM, CI2RDM  )

    write (6,*) "Building MCSCF gradient ..."
    call CIMCSCF_gradient( CI1RDM, CI2RDM, fock, gradient )

    write (6,*) "Building MCSCF hessian (diagonal) ..."
    call CIMCSCF_hessian( CI1RDM, CI2RDM, fock, hessian )

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

  subroutine CIMCSCF_energy( CI1RDM, CI2RDM  )
    implicit none
    type(matrix), allocatable, intent(inout) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(inout) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    integer :: p,q,r,s, pqrs, pq, rs, pq_aux
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

    open(unit=1018, file="2rdm", status="replace", form="formatted")
    open(unit=1019, file="2rdm.inter", status="replace", form="formatted")

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
              write (1018,"(I2,I2,I2,I2,I6,F12.8,F12.8 )" ) p,q,r,s, pqrs_rdm, CI2RDM(spi,spi)%values(pqrs_rdm), CIcore_instance%fourCenterIntegrals(spi,spi)%values(pqrs,1)

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

                write (1019,"(I2,I2,I2,I2,I6,F12.8,F12.8 )" ) p,q,r,s, pqrs_rdm, CI2RDM(spi,spj)%values(pqrs_rdm), &
                CIcore_instance%fourCenterIntegrals(spi,spj)%values(pqrs,1)

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

    close(1018)
    close(1019)

  end subroutine CIMCSCF_energy

  subroutine CIMCSCF_gradient( CI1RDM, CI2RDM, fock, gradient )
    implicit none
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(in) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: gradient(:) ! species % numcontractions, numcontractions
    integer :: p,q,r,s,t, pr, pr_rdm, qr, qr_aux, st, st_rdm, qrst, prst_rdm, st_aux
    integer :: pqpq_rdm, pq_rdm 
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j

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
                  !! or stpr_rdm

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

      print *, "gradient for spi", spi
      call Matrix_show (gradient(spi))

    enddo ! spi


  end subroutine CIMCSCF_gradient

  subroutine CIMCSCF_hessian( CI1RDM, CI2RDM, fock, hessian )
    implicit none
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    type(vector), allocatable, intent(in) :: CI2RDM(:,:) ! species, species % numcontractions, numcontractions, numcontractions, numcontractions
    type(matrix), allocatable, intent(in) :: fock(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(inout) :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable :: W_I(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable :: W_II(:) ! species % numcontractions, numcontractions (diagonal)
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: x, y, xy, xy_aux, m, n, yy, yy_aux, mn, xm, ym, ym_aux, yn, xmyn, xymn, ymyn, yymn
    integer :: xy_rdm, xx_rdm, mn_rdm, xn_rdm, mx_rdm, my_rdm, xm_rdm 
    integer :: xnmx_rdm, xnmy_rdm, xnxm_rdm, xnym_rdm, xxmn_rdm, xymn_rdm, ym_rdm
    integer :: p, q

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( hessian(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( hessian( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    allocate( W_I(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( W_I( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo
    allocate( W_II(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( W_II( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    ! W_I
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do x = 1, numberOfContractions_i

        xx_rdm = ( x - 1) * numberOfContractions_i + x
        do y = 1, numberOfContractions_i

          W_I(spi)%values(x,y) = W_I(spi)%values(x,y) + 2.0_8 * CIcore_instance%twoCenterIntegrals(spi)%values(y,y) * &
                                                                CI1RDM(spi,1)%values(x,x)

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

              W_I(spi)%values( x, y ) = W_I(spi)%values( x, y ) + 2.0_8 * &
                                                                  CIcore_instance%fourCenterIntegrals(spi,spi)%values(yymn,1) * &
                                                                  CI2RDM(spi,spi)%values(xxmn_rdm)

              yn = CIcore_instance%twoIndexArray(spi)%values(y,n)
              ymyn = CIcore_instance%fourIndexArray(spi)%values(ym,yn)

              xn_rdm = x + ( n - 1) * numberOfContractions_i
              xnmx_rdm = CIcore_two2one ( xn_rdm, mx_rdm )
              xnxm_rdm = CIcore_two2one ( xn_rdm, xm_rdm )

              W_I(spi)%values( x, y ) = W_I(spi)%values( x, y ) + 2.0_8 * &
                                                        CIcore_instance%fourCenterIntegrals(spi,spi)%values(ymyn,1) * &
                                                        ( CI2RDM(spi,spi)%values(xnmx_rdm) + CI2RDM(spi,spi)%values(xnxm_rdm) )
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
                W_I(spi)%values( x, y ) = W_I(spi)%values( x, y ) + 2.0_8 * &
                                                                   CIcore_instance%fourCenterIntegrals(spi,spj)%values(yymn,1) * &
                                                                   CI2RDM(spj,spi)%values(xxmn_rdm)

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
                W_I(spi)%values( x, y ) = W_I(spi)%values( x, y ) + 2.0_8 * &
                                                                   CIcore_instance%fourCenterIntegrals(spi,spj)%values(yymn,1) * &
                                                                   CI2RDM(spi,spj)%values(xxmn_rdm)
              enddo ! n 
            enddo ! m
          enddo ! spj 



        enddo ! q

        W_I(spi)%values( x, x ) = W_I(spi)%values( x, x ) + 2.0_8 * fock(spi)%values( x, x )

      enddo ! x

    enddo ! spi

    ! W_II
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do x = 1, numberOfContractions_i

        do y = 1, numberOfContractions_i

          W_II(spi)%values(y,x) = W_I(spi)%values(y,x) + 2.0_8 * CIcore_instance%twoCenterIntegrals(spi)%values(x,y) * &
                                                                 CI1RDM(spi,1)%values(y,x)

          xy_rdm = x + ( y - 1) * numberOfContractions_i 
          xy = CIcore_instance%twoIndexArray(spi)%values(x,y)
                                                                     
          do m = 1, numberOfContractions_i

            xm = CIcore_instance%twoIndexArray(spi)%values(x,m)
            my_rdm = m + ( y - 1) * numberOfContractions_i 
            ym_rdm = y + ( m - 1) * numberOfContractions_i

            do n = 1, numberOfContractions_i

              mn_rdm = m + ( n - 1) * numberOfContractions_i 
              xymn_rdm = CIcore_two2one ( xy_rdm, mn_rdm )

              mn = CIcore_instance%twoIndexArray(spi)%values(m,n)
              xymn = CIcore_instance%fourIndexArray(spi)%values(xy,mn)

              W_II(spi)%values( y, x ) = W_II(spi)%values( y, x ) + 2.0_8 * &
                                                            CIcore_instance%fourCenterIntegrals(spi,spi)%values(xymn,1) * &
                                                            CI2RDM(spi,spi)%values(xymn_rdm)

              yn = CIcore_instance%twoIndexArray(spi)%values(y,n)
              xmyn = CIcore_instance%fourIndexArray(spi)%values(xm,yn)

              xn_rdm = x + ( n - 1) * numberOfContractions_i
              xnmy_rdm = CIcore_two2one ( xn_rdm, my_rdm )
              xnym_rdm = CIcore_two2one ( xn_rdm, ym_rdm )

              W_II(spi)%values( y, x ) = W_II(spi)%values( y, x ) + 2.0_8 * &
                                                             CIcore_instance%fourCenterIntegrals(spi,spi)%values(xmyn,1) * &
                                                           ( CI2RDM(spi,spi)%values(xnmy_rdm) + CI2RDM(spi,spi)%values(xnym_rdm) )
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
                W_II(spi)%values( y, x ) = W_II(spi)%values( y, x ) + 2.0_8 * &
                                                                   CIcore_instance%fourCenterIntegrals(spj,spi)%values(xymn,1) * &
                                                                   CI2RDM(spj,spi)%values(xymn_rdm)
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
                W_II(spi)%values( y, x ) = W_II(spi)%values( y, x ) + 2.0_8 * &
                                                                   CIcore_instance%fourCenterIntegrals(spi,spj)%values(xymn,1) * &
                                                                   CI2RDM(spi,spj)%values(xymn_rdm)
              enddo ! n 
            enddo ! m
          enddo ! spj 

        enddo ! y

        W_II(spi)%values( x, x ) = W_II(spi)%values( x, x ) + 2.0_8 * fock(spi)%values( x, x )

      enddo ! x
    enddo ! spi

    !! building the hessian
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do p = 1, numberOfContractions_i
        do q = 1, numberOfContractions_i
          hessian(spi)%values( p, q ) = hessian(spi)%values( p, q ) + W_I(spi)%values( p, q ) &
                                                                    - W_II(spi)%values( q, p ) &
                                                                    - W_II(spi)%values( p, q ) &
                                                                    + W_I(spi)%values( q, p ) 
        enddo ! q
      enddo ! p 
  
      !print *, "hessian for spi", spi
      !call Matrix_show (hessian(spi))

    enddo ! spi

    deallocate( W_II )
    deallocate( W_I )

  end subroutine CIMCSCF_hessian

  subroutine CIMCSCF_newtonRaphson( gradient,  hessian, rotations )
    implicit none
    type(matrix), allocatable, intent(in) :: gradient(:) ! species % numcontractions, numcontractions
    type(matrix), allocatable, intent(in) :: hessian(:) ! species % numcontractions, numcontractions (diagonal)
    type(matrix), allocatable, intent(inout) :: rotations(:) ! species % numcontractions, numcontractions 
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: p, q
    real(8) :: epsilon

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    epsilon = 1.0E-2

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
          rotations(spi)%values(p,q) = rotations(spi)%values(p,q) - gradient(spi)%values(p,q) / &
                                                                    ( hessian(spi)%values(p,q) + epsilon) 
        enddo ! q
      enddo ! p 
  
      print *, "rotations for spi", spi
      call Matrix_show (rotations(spi))

    enddo ! spi

  end subroutine CIMCSCF_newtonRaphson 

  subroutine CIMCSCF_buildUnitaryMatrix ( rotations, UnitaryMatrix ) 
   implicit none
    type(matrix), allocatable, intent(inout) :: rotations(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable, intent(inout) :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: p, q
    real(8), allocatable :: A(:,:)  
    integer, allocatable :: IPIV(:)
    integer :: info

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( unitaryMatrix(numberOfSpecies) )
    do spi = 1, numberOfSpecies
       numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      call Matrix_constructor ( unitaryMatrix( spi ), int( numberOfContractions_i, 8), int( numberOfContractions_i, 8), 0.0_8 ) 
    enddo

    !! antisymmetrizing the orbital rotations
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )
      do p = 1, numberOfContractions_i
        do q = p + 1, numberOfContractions_i
          rotations(spi)%values(q,p) = - rotations(spi)%values(p,q)
        enddo ! q
      enddo ! p
    enddo ! spi

    !! Caley transformation ( exponential of a matrix)
    do spi = 1, numberOfSpecies
      numberOfContractions_i = MolecularSystem_getTotalNumberOfContractions( spi )

      !! building A = I - 0.5*rotations and B = I + 0.5*rotations. Keep B in unitary matrix U
      allocate ( A ( numberOfContractions_i, numberOfContractions_i ) )
      allocate ( IPIV ( numberOfContractions_i ) )
      A = 0.0_8
      IPIV = 0

      do p = 1, numberOfContractions_i 
        A(p, p ) = 1.0_8
        unitaryMatrix(spi)%values(p, p) = 1.0_8
      enddo ! p 

      do q = 1, numberOfContractions_i
        do p = 1, numberOfContractions_i 
          A(p, q) = A(p, q) - 0.5_8 * rotations(spi)%values(p, q)
          unitaryMatrix(spi)%values(p, q) = unitaryMatrix(spi)%values(p, q) + 0.5_8 * rotations(spi)%values(p, q)
        enddo ! p 
      enddo ! q

      !! LAPACK LU Factorization of A
      !! DGETRF computes: A = P * L * U
      call dgetrf( numberOfContractions_i, &
                   numberOfContractions_i, &
                   A, &
                   numberOfContractions_i, &
                   IPIV, &
                   info)

      if (info /= 0) then
          write(*,*) 'Error in DGETRF during Cayley Transform! INFO = ', info
      end if

      !! LAPACK Linear Solver
      !! DGETRS solves A * U = B (where U holds B initially, and is overwritten with the solution U)
      call dgetrs('N', & !! 'N' specifies no transpose on A.
                  numberOfContractions_i, &
                  numberOfContractions_i, &
                  A, & 
                  numberOfContractions_i, &
                  IPIV, &
                  unitaryMatrix(spi)%values, & 
                  numberOfContractions_i, &
                  info)
      if (info /= 0) then
          write(*,*) 'Error in DGETRS during Cayley Transform! INFO = ', info
      end if

      deallocate ( IPIV )
      deallocate ( A )

      !print *, "U"
      !call Matrix_show ( unitaryMatrix(spi) )

    enddo ! spi

  end subroutine CIMCSCF_buildUnitaryMatrix 

  subroutine CIMCSCF_rotateCoefficients ( UnitaryMatrix, CI1RDM ) 
   implicit none
    type(matrix), allocatable, intent(in) :: unitaryMatrix(:) ! species % numcontractions, numcontractions 
    type(matrix), allocatable, intent(in) :: CI1RDM(:,:) ! species, state % numcontractions, numcontractions
    integer :: spi, spj, numberOfSpecies
    integer :: numberOfContractions_i, numberOfOccupiedOrbitals_i
    integer :: numberOfContractions_j, numberOfOccupiedOrbitals_j
    integer :: p, q, r
    type(matrix), allocatable :: coefficients_old(:)
    type(matrix), allocatable :: coefficients_new(:)
    integer :: unit 
    integer :: wfnunit
    character(50) :: file, speciesName, auxstring
    character(50) :: wfnfile
    character(100) :: arguments(2)
    type(Matrix), allocatable :: hcoreMatrix(:)
    type(matrix), allocatable :: kineticMatrix(:), attractionMatrix(:), externalPotMatrix(:)
    real(8) :: auxvalue

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( coefficients_old(numberOfSpecies) )
    allocate( coefficients_new(numberOfSpecies) )
    allocate( hcoreMatrix(numberOfSpecies) )
    allocate( kineticMatrix(numberOfSpecies) )
    allocate( attractionMatrix(numberOfSpecies) )
    allocate( externalPotMatrix(numberOfSpecies) )

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

      print *, " old coefficients ", spi
      call Matrix_show( coefficients_old(spi) )

      call Matrix_constructor( coefficients_new(spi), int(numberOfContractions_i, 8), int(numberOfContractions_i, 8), 0.0_8 )

      do q = 1, numberOfContractions_i
        do p = 1, numberOfContractions_i 
          auxvalue = 0.0_8
          do r = 1, numberOfContractions_i 
            auxvalue = auxvalue + coefficients_old(spi)%values(p,r) * unitaryMatrix(spi)%values(r,q) 
          enddo
          coefficients_new(spi)%values(p,q) = auxvalue
          !coefficients_new(spi)%values(p,q) = coefficients_old(spi)%values(p,q) 
        enddo
      enddo

      print *, " new coefficients ", spi
      call Matrix_show ( coefficients_new(spi) )
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

  end subroutine CIMCSCF_rotateCoefficients

end module CIMCSCF_

