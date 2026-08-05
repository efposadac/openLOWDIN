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
 
module CIdensity_
  use Exception_
  use Matrix_
  use Vector_
  use MolecularSystem_
  use IndexMap_
  use InputCI_
  use omp_lib
  use CIcore_
  use CISCI_
 
  implicit none

  public :: &
       CIdensity_compute, &
       CIdensity_1RDM_SCI, &
       CIdensity_2RDM_SCI
  private :: &
       CIdensity_1RDM_CI, &
       CIdensity_energyTerms, &
       CIdensity_naturalOrbitals

contains

  subroutine CIdensity_compute()
    implicit none
    type(matrix), allocatable :: ciDensityMatrix(:,:)
    type(matrix), allocatable :: coefficients(:)
    integer :: state, species, numberOfSpecies
    integer :: k, n
    integer :: numberOfContractions, numberOfOccupiedOrbitals
    real(8) :: timeDA, timeDB

    !$  timeDA = omp_get_wtime()

    if ( .not. CIcore_instance%isInstanced ) return

    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) "BUILDING CI DENSITY MATRICES"
    write(6,*) "-----------------------------------------------------------------------"
    write(6,*) ""

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    !! constructor
    allocate( ciDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), &
              coefficients(numberOfSpecies) )

    !! matrix allocation
    do species=1, numberOfSpecies
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
       numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(species)

       do state = 1, CONTROL_instance%CI_STATES_TO_PRINT
          call Matrix_constructor ( ciDensityMatrix(species,state) , &
               int(numberOfContractions,8), &
               int(numberOfContractions,8),  0.0_8 )
          do k = 1, numberOfOccupiedOrbitals
            ciDensityMatrix(species,state)%values( k, k)=1.0_8
          end do
       end do

    end do

    !! computing all density related stuff
    if ( CONTROL_instance%CI_ONE_REDUCED_DENSITY_MATRIX ) then 
      if ( CONTROL_instance%CI_SELECTIVE_METHOD == "NONE" ) then
        call CIdensity_1RDM_CI( ciDensityMatrix )
      else
        call CIdensity_1RDM_SCI( ciDensityMatrix )
      endif
    endif
    !if ( CONTROL_instance%CI_TWO_REDUCED_DENSITY_MATRIX ) call CIdensity_2RDM_SCI( ) 
    if ( CONTROL_instance%CI_STATES_TO_PRINT > 0 ) call CIdensity_energyTerms( coefficients, ciDensityMatrix )
    if ( CONTROL_instance%CI_STATES_TO_PRINT > 0 .and. CONTROL_instance%CI_NATURAL_ORBITALS) call CIdensity_naturalOrbitals( coefficients, ciDensityMatrix )

    !! destructor
    do species = 1, numberOfSpecies
      do state = 1, CONTROL_instance%CI_STATES_TO_PRINT
        call Matrix_destructor(ciDensityMatrix(species, state))
      end do
      call Matrix_destructor(coefficients(species))
    end do
    deallocate( ciDensityMatrix, coefficients )

    !$  timeDB = omp_get_wtime()
    !$  write(*,"(A,F10.4,A4)") "** TOTAL Elapsed Time for Building density matrices: ", timeDB - timeDA ," (s)"

  end subroutine CIdensity_compute

  subroutine CIdensity_1RDM_CI( ciDensityMatrix )
    implicit none
    integer :: i, j, k, l, mu, nu, n
    integer :: factor
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state, species, orbital, orbitalA, orbitalB
    type(matrix), allocatable, intent(inout) :: ciDensityMatrix(:,:)
    type(matrix), allocatable :: auxDensMatrix(:,:)
    integer :: numberOfSpecies
    integer(8) :: numberOfConfigurations, c
    integer, allocatable :: cilevel(:), cilevelA(:)
    integer(8), allocatable :: indexConf(:)
    type(ivector), allocatable :: stringAinB(:)
    integer :: s, ss, ci, auxnumberOfSpecies
    integer, allocatable :: coupling(:)
    integer :: a, b, AA, BB, bj
    integer :: u, uu, ssize
    integer(8), allocatable :: indexConfA(:)
    integer(8), allocatable :: indexConfB(:)
    integer(8), allocatable :: jj(:)

    !!Iterators: i,j - Configurations .... k,l - molecular orbitals .... mu,nu - atomic orbitals ... n - threads

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    numberOfConfigurations = CIcore_instance%numberOfConfigurations 

    allocate ( auxDensMatrix(numberOfSpecies,CIcore_instance%nproc))

    !! matrix allocation
    do species=1, numberOfSpecies
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )

       do n = 1, CIcore_instance%nproc
          call Matrix_constructor ( auxDensMatrix(species,n) , &
               int(numberOfContractions,8), &
               int(numberOfContractions,8),  0.0_8 )
       end do
    end do

  
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
  
    deallocate ( indexConf )
    allocate ( coupling ( numberOfSpecies ) )
    allocate ( indexConfA ( numberOfSpecies ) )
    allocate ( indexConfB ( numberOfSpecies ) )
    allocate ( jj ( numberOfSpecies ) )

    indexConfA = 0
    indexConfB = 0
    jj = 0

    !! Building the CI reduced density matrix in the molecular orbital representation in parallel
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

            !!Diagonal contributions
            do species=1, numberOfSpecies
               numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(species)

               do k=1, numberOfOccupiedOrbitals

                  !!Occupied orbitals
                  auxDensMatrix(species,n)%values(k,k)=auxDensMatrix(species,n)%values(k,k) - CIcore_instance%eigenVectors%values(i,state)**2
                  orbital =  CIcore_instance%strings(species)%values(k,indexConfA(species))
                  !!Unoccupied orbitals

                  auxDensMatrix(species,n)%values(orbital,orbital)=auxDensMatrix(species,n)%values(orbital,orbital) + CIcore_instance%eigenVectors%values(i,state)**2

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
                    end do
                    ssize = ssize + CIcore_instance%numberOfStrings(ss)%values(ci)
                  end if

                end do
              end do

            end do 

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
              
    end do !! number of CI states

    deallocate ( jj )
    deallocate ( indexConfB )
    deallocate ( indexConfA )
    deallocate ( coupling )
    deallocate ( cilevel )
    deallocate ( cilevelA )
    deallocate ( CIcore_instance%allIndexConf )
    deallocate ( stringAinB )

    do species = 1, numberOfSpecies
      do n = 1, CIcore_instance%nproc
        call Matrix_destructor(auxDensMatrix(species, n))
      end do
    end do
    deallocate( auxDensMatrix )

  end subroutine CIdensity_1RDM_CI

  subroutine CIdensity_1RDM_SCI( ciDensityMatrix )
    implicit none
    integer :: i, j, k, l, mu, nu, n
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state, species, orbital, orbitalA, orbitalB
    type(matrix), allocatable, intent(inout) :: ciDensityMatrix(:,:)
    integer :: numberOfSpecies
    integer(8) :: numberOfConfigurations, a, b, c
    !! Auxiliary variables for SCI
    integer(1), allocatable :: couplingS(:)
    integer :: spi
    integer :: pi
    integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA
    integer :: diffOrbi(4)

    !!Iterators: i,j - Configurations .... k,l - molecular orbitals .... mu,nu - atomic orbitals ... n - threads

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
  
    numberOfConfigurations = CIcore_instance%numberOfConfigurations 

    !CIcore_instance%eigenVectors%values(1,1) = 1.0_8 
  
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )
    
    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    !! Building the CI reduced density matrix in the molecular orbital representation in parallel
    do state=1, CONTROL_instance%CI_STATES_TO_PRINT

      !do a = 1, 1
      do a = 1, CIcore_instance%numberOfConfigurations
        n = 1

        do spi = 1, numberOfSpecies 
          oia = 0 

          !!orbA(spi)%values = CISCI_instance%targetOrb(spi,a)%values
          orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,a)

          !! build auxiliary vectors of occupied and virtuals orbitals
          do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
            if ( orbA(spi)%values(pi) == 1 ) then
              oia = oia + 1
              occA(spi)%values(oia) = pi
            end if
          enddo
     
        enddo

        !!Diagonal contributions
        do spi = 1, numberOfSpecies
          numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(spi)

          do k = 1, numberOfOccupiedOrbitals

            !!Occupied orbitals
            !ciDensityMatrix(spi,state)%values(k,k) = ciDensityMatrix(spi,state)%values(k,k) - CIcore_instance%eigenVectors%values(a,state)**2
            orbital = occA(spi)%values(k) 

            !!Unoccupied orbitals
            ciDensityMatrix(spi,state)%values(orbital,orbital) = ciDensityMatrix(spi,state)%values(orbital,orbital) + CIcore_instance%eigenVectors%values(a,state)**2
    
           end do
         end do

        !!Off Diagonal contributions
        !do b = a + 1, 1
        do b = a + 1, CICore_instance%numberOfConfigurations 

          do spi = 1, numberOfSpecies 
            !orbB(spi)%values = CISCI_instance%targetOrb(spi,b)%values
            orbB(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,b)
          enddo

          !! determinate number of diff orbitals
          couplingS = 0
          do spi = 1, numberOfSpecies
            couplingS(spi) = couplingS(spi) + CIcore_instance%numberOfOccupiedOrbitals%values(spi) &
                              - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) 
          end do
    
          !! just single particle diff 
          if ( sum(couplingS) == 1 ) then

            do spi = 1, numberOfSpecies 
              oib = 0 
              !! build auxiliary vectors of occupied and virtuals orbitals
              do pi = 1, CIcore_instance%numberOfActiveOrbitals%values(spi)
                if ( orbB(spi)%values(pi) == 1 ) then
                  oib = oib + 1
                  occB(spi)%values(oib) = pi
                end if
              enddo
            enddo

            do i = 1, numberOfSpecies
                if ( couplingS(i) == 1 ) spi = i
            end do

            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )

            ciDensityMatrix(spi,state)%values( diffOrbi(1), diffOrbi(3) ) = ciDensityMatrix(spi,state)%values( diffOrbi(1), diffOrbi(3) ) + &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(a,state) * &
                                                            CIcore_instance%eigenVectors%values(b,state)
           cidensitymatrix(spi,state)%values( difforbi(3), difforbi(1) ) = ciDensityMatrix(spi,state)%values( diffOrbi(3), diffOrbi(1) ) + &
                                                            factorA * &
                                                            CIcore_instance%eigenVectors%values(a,state) * &
                                                            CIcore_instance%eigenVectors%values(b,state)

           
          endif !! coupling 
        enddo !! b
      end do !! a

    end do !! number of CI states

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( occA )
    deallocate ( occB )
    deallocate ( orbA )
    deallocate ( orbB )
    deallocate ( couplingS )
       
  end subroutine CIdensity_1RDM_SCI

  subroutine CIdensity_2RDM_SCI( ciDensityMatrix )
    implicit none
    type(vector), allocatable, intent(inout) :: ciDensityMatrix(:,:)
    integer :: mu, nu !! particles
    integer :: p, pp, pppp, q, r, rr, pprr, pr, rp, prrp, rrpp, rppr !! orbitals general
    integer :: pp_aux, ij_aux, ji_aux, kl_aux, lk_aux
    integer :: i, j, k, l
    integer :: iq, qj, ij, ppij, iqqj
    integer :: jq, qi, ji, ppji, jqqi, ijpp, jipp
    integer :: kl, ijkl, klij
    integer :: il, kj, ilkj, kjil
    integer :: lk, lkji
    integer :: jilk
    integer :: jk, li, lijk
    integer(8) :: II, JJ !! configurations
    integer(8) :: numberOfConfigurations
    integer :: numberOfOccupiedOrbitals_spi, numberOfOccupiedOrbitals_spj
    integer :: numberOfOrbitals_spi, numberOfOrbitals_spj
    integer :: state
    integer :: species, numberOfSpecies
    integer :: spi, spj
    !! Auxiliary variables for SCI
    integer(1), allocatable :: couplingS(:)
    !integer :: oia, oib
    type (ivector), allocatable :: occA(:), occB(:)
    type (ivector), allocatable :: orbA(:), orbB(:)
    integer :: factorA, factorB
    integer :: diffOrbi(4), diffOrbj(4)

    !!Iterators: i,j - Configurations .... k,l - molecular orbitals .... mu,nu - atomic orbitals ... 

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
  
    numberOfConfigurations = CIcore_instance%numberOfConfigurations 
  
    allocate ( occA ( numberOfSpecies ) )
    allocate ( occB ( numberOfSpecies ) )
    allocate ( orbA ( numberOfSpecies ) )
    allocate ( orbB ( numberOfSpecies ) )
    allocate ( couplingS ( numberOfSpecies ) )
    
    do spi = 1, numberOfSpecies
      call Vector_constructorInteger ( occA(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 ) ! use core here? yes
      call Vector_constructorInteger ( occB(spi), CIcore_instance%numberOfOccupiedOrbitals%values(spi), 0 )
      call Vector_constructorInteger ( orbA(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
      call Vector_constructorInteger ( orbB(spi), CIcore_instance%numberOfActiveOrbitals%values(spi),  0 ) 
    end do

    !CIcore_instance%eigenVectors%values(1,1) = 1.0_8 

    !! Building the CI reduced density matrix in the molecular orbital representation in parallel
    do state = 1, CONTROL_instance%CI_NUMBER_OF_STATES
      !do II = 1, 1!CIcore_instance%numberOfConfigurations
      do II = 1, CIcore_instance%numberOfConfigurations

        do spi = 1, numberOfSpecies 
          orbA(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,II)
          occA(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:,II) 
        enddo

        !!Diagonal contributions
        do spi = 1, numberOfSpecies
          numberOfOccupiedOrbitals_spi = CIcore_instance%numberOfOccupiedOrbitals%values(spi)
          numberOfOrbitals_spi = CIcore_instance%numberOfOrbitals%values(spi)

          !! alpha-alpha
          do mu = 1, numberOfOccupiedOrbitals_spi
            p = occA(spi)%values(mu) 
            pp = p + ( p - 1 )*numberOfOrbitals_spi  

              do nu = mu, numberOfOccupiedOrbitals_spi
              r = occA(spi)%values(nu) 
              rr = r + ( r - 1 )*numberOfOrbitals_spi 

              pprr = CIcore_two2one ( pp, rr )
              ciDensityMatrix(spi,spi)%values(pprr) = ciDensityMatrix(spi,spi)%values(pprr) + &
                                                        CIcore_instance%eigenVectors%values(II,state)**2
              rp = r + ( p - 1 )*numberOfOrbitals_spi  
              pr = p + ( r - 1 )*numberOfOrbitals_spi 
              prrp = CIcore_two2one ( pr, rp )
                ciDensityMatrix(spi,spi)%values(prrp) = ciDensityMatrix(spi,spi)%values(prrp) - &
                                                         CIcore_instance%eigenVectors%values(II,state)**2
            end do !nu
          end do ! mu

          !! alpha-beta 
          do spj = spi + 1, numberOfSpecies
            numberOfOccupiedOrbitals_spj = CIcore_instance%numberOfOccupiedOrbitals%values(spj)
            numberOfOrbitals_spj = CIcore_instance%numberOfOrbitals%values(spj)
            do mu = 1, numberOfOccupiedOrbitals_spi
              p = occA(spi)%values(mu) 
              pp = p + ( p - 1 )*numberOfOrbitals_spi 

              do nu = 1, numberOfOccupiedOrbitals_spj
                r = occA(spj)%values(nu) 
                rr = r + ( r - 1 )*numberOfOrbitals_spj  

                pprr = pp + ( rr - 1 ) * numberOfOrbitals_spi * numberOfOrbitals_spi 
                ciDensityMatrix(spi,spj)%values(pprr) = ciDensityMatrix(spi,spj)%values(pprr) + &
                                                          CIcore_instance%eigenVectors%values(II,state)**2
              end do !nu
            end do ! mu
          end do ! spj

        enddo !spi

        !Off Diagonal contributions
        !do JJ = II + 1, 1!CICore_instance%numberOfConfigurations 
        do JJ = II + 1, CICore_instance%numberOfConfigurations 

          !print *, II, JJ
          do spi = 1, numberOfSpecies 
            orbB(spi)%values(:) = CISCI_instance%confTarget_orb(spi)%values(:,JJ)
            occB(spi)%values(:) = CISCI_instance%confTarget_occ(spi)%values(:,JJ) 
          enddo

          !! determinate number of diff orbitals
          couplingS = 0
          do spi = 1, numberOfSpecies
            numberOfOccupiedOrbitals_spi = CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            couplingS(spi) = couplingS(spi) + numberOfOccupiedOrbitals_spi &
                              - sum ( orbA(spi)%values(:) * orbB(spi)%values(:) ) 
          end do
      
          select case ( sum(couplingS) )
      
          !! one orbital different
          case (1) !! this is adding noise to alpha alpha repulsion
            do species = 1, numberOfSpecies
              if ( couplingS(species) == 1 ) spi = species
            end do

            numberOfOccupiedOrbitals_spi = CIcore_instance%numberOfOccupiedOrbitals%values(spi)
            numberOfOrbitals_spi = CIcore_instance%numberOfOrbitals%values(spi)
  
            diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )

            i = diffOrbi(1)
            j = diffOrbi(3)
            ij = i + ( j - 1 )*numberOfOrbitals_spi 
            ji = j + ( i - 1 )*numberOfOrbitals_spi 

            do mu = 1, numberOfOccupiedOrbitals_spi
              p = occA(spi)%values(mu) 
              pp = p + ( p - 1 )*numberOfOrbitals_spi 
              ppij = CIcore_two2one ( pp, ij )
  
              ciDensityMatrix(spi,spi)%values( ppij ) = ciDensityMatrix(spi,spi)%values( ppij ) + &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              q = occA(spi)%values(mu) 
              iq = i + ( q - 1 )*numberOfOrbitals_spi
              qj = q + ( j - 1 )*numberOfOrbitals_spi 
              iqqj = CIcore_two2one ( iq, qj )

              ciDensityMatrix(spi,spi)%values( iqqj ) = ciDensityMatrix(spi,spi)%values( iqqj ) - &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              ppji = CIcore_two2one ( pp, ji )
  
              ciDensityMatrix(spi,spi)%values( ppji ) = ciDensityMatrix(spi,spi)%values( ppji ) + &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              jq = j + ( q - 1 )*numberOfOrbitals_spi
              qi = q + ( i - 1 )*numberOfOrbitals_spi
              jqqi = CIcore_two2one ( jq, qi )

              ciDensityMatrix(spi,spi)%values( jqqi ) = ciDensityMatrix(spi,spi)%values( jqqi ) - &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

            enddo

            !! beta-beta, diff in alpha
            do spj = 1, spi - 1

              numberOfOccupiedOrbitals_spj = CIcore_instance%numberOfOccupiedOrbitals%values(spj)
              numberOfOrbitals_spj = CIcore_instance%numberOfOrbitals%values(spj)
              do mu = 1, numberOfOccupiedOrbitals_spj
                p = occA(spj)%values(mu)
                pp = p + ( p - 1 )*numberOfOrbitals_spj

                ppij = pp + ( ij - 1 ) * numberOfOrbitals_spj * numberOfOrbitals_spj

                ciDensityMatrix(spj,spi)%values( ppij ) = ciDensityMatrix(spj,spi)%values( ppij ) + &
                                                              factorA * & 
                                                              CIcore_instance%eigenVectors%values(II,state) * &
                                                              CIcore_instance%eigenVectors%values(JJ,state)

                ppji = pp + ( ji - 1 ) * numberOfOrbitals_spj * numberOfOrbitals_spj
                ciDensityMatrix(spj,spi)%values( ppji ) = ciDensityMatrix(spj,spi)%values( ppji ) + &
                                                              factorA * & 
                                                              CIcore_instance%eigenVectors%values(II,state) * &
                                                              CIcore_instance%eigenVectors%values(JJ,state)
              enddo !mu
            enddo !spj 

            !! alpha-beta, diff in alpha
            do spj = spi + 1, numberOfSpecies

              numberOfOccupiedOrbitals_spj = CIcore_instance%numberOfOccupiedOrbitals%values(spj)
              numberOfOrbitals_spj = CIcore_instance%numberOfOrbitals%values(spj)
              do mu = 1, numberOfOccupiedOrbitals_spj
                p = occA(spj)%values(mu)
                pp = p + ( p - 1 )*numberOfOrbitals_spj 
                ijpp = ij + ( pp - 1 ) * numberOfOrbitals_spi * numberOfOrbitals_spi 

                ciDensityMatrix(spi,spj)%values( ijpp ) = ciDensityMatrix(spi,spj)%values( ijpp ) + &
                                                              factorA * & 
                                                              CIcore_instance%eigenVectors%values(II,state) * &
                                                              CIcore_instance%eigenVectors%values(JJ,state)

                jipp = ji + ( pp - 1 ) * numberOfOrbitals_spi * numberOfOrbitals_spi 

                ciDensityMatrix(spi,spj)%values( jipp ) = ciDensityMatrix(spi,spj)%values( jipp ) + &
                                                              factorA * & 
                                                              CIcore_instance%eigenVectors%values(II,state) * &
                                                              CIcore_instance%eigenVectors%values(JJ,state)
              enddo !mu
            enddo !spj 

          !! two orbital different
          case(2)
  
            select case ( maxval(couplingS) )
            !! two orbital different, same species
            case (2)
              do species = 1, numberOfSpecies
                if ( couplingS(species) == 2 ) spi = species
              end do
  
              diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, occA(spi)%values, occB(spi)%values, factorA )
              numberOfOrbitals_spi = CIcore_instance%numberOfOrbitals%values(spi)

              i = diffOrbi(1) ! 1 diff orb in a
              j = diffOrbi(3) ! 1 diff orb in b
              k = diffOrbi(2) ! 2 diff orb in a
              l = diffOrbi(4) ! 2 diff orb in b

              ij = i + ( j - 1 )*numberOfOrbitals_spi 
              kl = k + ( l - 1 )*numberOfOrbitals_spi 

              ijkl = CIcore_two2one ( ij, kl )
  
              ciDensityMatrix(spi,spi)%values( ijkl ) = ciDensityMatrix(spi,spi)%values( ijkl ) + &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              il = i + ( l - 1 )*numberOfOrbitals_spi 
              kj = k + ( j - 1 )*numberOfOrbitals_spi 

              ilkj = CIcore_two2one ( il, kj )
  
              ciDensityMatrix(spi,spi)%values( ilkj ) = ciDensityMatrix(spi,spi)%values( ilkj ) - &
                                                            factorA * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              ji = j + ( i - 1 )*numberOfOrbitals_spi                                                        
              lk = l + ( k - 1 )*numberOfOrbitals_spi                                                         
                                                                                                             
              jilk= CIcore_two2one ( ji, lk )                                                                
                                                                                                             
              ciDensityMatrix(spi,spi)%values( jilk ) = ciDensityMatrix(spi,spi)%values( jilk ) + &          
                                                            factorA * &                                      
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)    
                                                                                                             
              li = l + ( i - 1 )*numberOfOrbitals_spi                                                       
              jk = j + ( k - 1 )*numberOfOrbitals_spi                                                        
                                                                                                             
              lijk = CIcore_two2one ( li, jk )                                                               
                                                                                                             
              ciDensityMatrix(spi,spi)%values( lijk ) = ciDensityMatrix(spi,spi)%values( lijk ) - &          
                                                            factorA * &                                      
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)    

            !! two orbital different, different species
            case (1)
              do species = 1, numberOfSpecies
                if ( couplingS(species) == 1 ) then 
                  spi = species
                  exit
                end if
              end do
              do species = spi + 1, numberOfSpecies
                if ( couplingS(species) == 1 ) spj = species
              end do
  
              diffOrbi = CISCI_getDiffOrbitals ( spi, orbA(spi)%values, orbB(spi)%values, &
                                                 occA(spi)%values, occB(spi)%values, factorA )
              diffOrbj = CISCI_getDiffOrbitals ( spj, orbA(spj)%values, orbB(spj)%values, &
                                                 occA(spj)%values, occB(spj)%values, factorB )

              numberOfOrbitals_spi = CIcore_instance%numberOfOrbitals%values(spi)
              numberOfOrbitals_spj = CIcore_instance%numberOfOrbitals%values(spj)

              i = diffOrbi(1) ! 1 diff orb in a spi
              j = diffOrbi(3) ! 1 diff orb in b spi
              k = diffOrbj(1) ! 1 diff orb in a spj
              l = diffOrbj(3) ! 1 diff orb in b spj

              ij = i + ( j - 1 )*numberOfOrbitals_spi 
              kl = k + ( l - 1 )*numberOfOrbitals_spj 
              ijkl = ij + ( kl - 1 ) * numberOfOrbitals_spi * numberOfOrbitals_spi 

              ciDensityMatrix(spi,spj)%values( ijkl ) = ciDensityMatrix(spi,spj)%values( ijkl ) + &
                                                            factorA * factorB * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

              ji = j + ( i - 1 )*numberOfOrbitals_spi 
              lk = l + ( k - 1 )*numberOfOrbitals_spj 

              jilk = ji + ( lk - 1 ) * numberOfOrbitals_spi * numberOfOrbitals_spi

              ciDensityMatrix(spi,spj)%values( jilk ) = ciDensityMatrix(spi,spj)%values( jilk ) + &
                                                            factorA * factorB * & 
                                                            CIcore_instance%eigenVectors%values(II,state) * &
                                                            CIcore_instance%eigenVectors%values(JJ,state)

            end select ! maxval(couplingS)  

          end select ! sum(couplingS) 

        enddo !! JJ
      end do !! II
    end do !! number of CI states

    do spi = 1, numberOfSpecies
      call Vector_destructorInteger ( occA(spi) ) 
      call Vector_destructorInteger ( occB(spi) )
      call Vector_destructorInteger ( orbA(spi) ) 
      call Vector_destructorInteger ( orbB(spi) ) 
    end do

    deallocate ( occA )
    deallocate ( occB )
    deallocate ( orbA )
    deallocate ( orbB )
    deallocate ( couplingS )

  end subroutine CIdensity_2RDM_SCI

  subroutine CIdensity_energyTerms(  coefficients, ciDensityMatrix )
    implicit none

    integer :: unit 
    integer :: wfnunit
    character(50) :: file, speciesName, auxstring
    character(50) :: wfnfile
    character(100) :: arguments(2)
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state
    integer :: orbital, orbitalA, orbitalB
    integer :: species, numberOfSpecies
    integer :: mu, nu
    integer :: k, l
    type(matrix), allocatable, intent(inout) :: coefficients(:)
    type(matrix), allocatable :: atomicDensityMatrix(:,:)
    type(matrix), allocatable :: kineticMatrix(:), attractionMatrix(:), externalPotMatrix(:)
    type(matrix), allocatable, intent(in) :: ciDensityMatrix(:,:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate( atomicDensityMatrix(numberOfSpecies,CONTROL_instance%CI_STATES_TO_PRINT), &
              kineticMatrix(numberOfSpecies), &
              attractionMatrix(numberOfSpecies), &
              externalPotMatrix(numberOfSpecies) )

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    !! matrix allocation
    do species = 1, numberOfSpecies
       speciesName = MolecularSystem_getNameOfSpecies(species)
       
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
       numberOfOccupiedOrbitals = CIcore_instance%numberOfOccupiedOrbitals%values(species)

       arguments(2) = speciesName

       arguments(1) = "COEFFICIENTS"
       coefficients(species) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                                                 columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))

       arguments(1) = "KINETIC"
       kineticMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                                                  columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))
       
       arguments(1) = "ATTRACTION"
       attractionMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                                                     columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))
       arguments(1) = "EXTERNAL-POTENTIAL"
       if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
       externalPotMatrix(species) = Matrix_getFromFile(unit=wfnUnit, rows=int(numberOfContractions, 8), &
                                                        columns=int(numberOfContractions, 8), binary=.true., arguments=arguments(1:2))
    end do
     
    close(wfnUnit)

    !! Open file - to write density matrices
    unit = 29
    file = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    open(unit = unit, file=trim(file), status="unknown", form="formatted")
       
    !! Building the CI reduced density matrix in the atomic orbital representation       
    do species = 1, numberOfSpecies
      speciesName = MolecularSystem_getNameOfSpecies(species)
      numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )

      do state = 1, CONTROL_instance%CI_STATES_TO_PRINT
         
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
      
        write(auxstring,*) state
        arguments(2) = speciesName
        arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 
            
        call Matrix_writeToFile ( atomicDensityMatrix(species,state), unit , arguments=arguments(1:2) )

        end do !state
      end do !species

      write(6,*) "-----------------------------------------------------------------------"
      write(*,*) " ONE BODY ENERGY CONTRIBUTIONS:"
      write(*,*) ""
      do state=1, CONTROL_instance%CI_STATES_TO_PRINT
         write(*,*) " STATE: ", state
         do species=1, molecularSystem_instance%numberOfQuantumSpecies
            write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%symbol ) // &
                 " Kinetic energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*kineticMatrix(species)%values)
            write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%symbol ) // &
                 "/Fixed interact. energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*attractionMatrix(species)%values)
            if( CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) &
                 write(*,"(A38,F25.12)") trim( MolecularSystem_instance%species(species)%symbol) // &
                 " Ext Pot energy = ", sum(transpose(atomicDensityMatrix(species,state)%values)*externalPotMatrix(species)%values)
            print *, ""
         end do
         print *, ""
      end do ! state

     close(unit)

      do species = 1, numberOfSpecies
        do state = 1, CONTROL_instance%CI_STATES_TO_PRINT
          call Matrix_destructor(atomicDensityMatrix(species, state))
        end do

        call Matrix_destructor(kineticMatrix(species))
        call Matrix_destructor(attractionMatrix(species))
        call Matrix_destructor(externalPotMatrix(species))

      end do

      deallocate( kineticMatrix, attractionMatrix, externalPotMatrix, atomicDensityMatrix )

  end subroutine CIdensity_energyTerms

  subroutine CIdensity_naturalOrbitals( coefficients, ciDensityMatrix )
    implicit none
    integer :: unit 
    character(50) :: file, speciesName, auxstring
    character(100) :: arguments(2)
    integer :: numberOfOrbitals, numberOfContractions, numberOfOccupiedOrbitals
    integer :: state
    integer :: orbital, orbitalA, orbitalB
    integer :: species, numberOfSpecies
    integer :: k, j, u
    type(vector) :: auxdensityEigenValues
    type(vector) :: densityEigenValues
    type(matrix) :: auxdensityEigenVectors 
    type(matrix) :: densityEigenVectors
    type(matrix), allocatable, intent(in) :: coefficients(:)
    type(matrix), allocatable, intent(in) :: ciDensityMatrix(:,:)

    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()
    !! Open file - to write density matrices
    unit = 30
    file = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    open(unit = unit, file=trim(file), status="unknown", position="append", form="formatted")

    write(6,*) "-----------------------------------------------------------------------"
    write(*,*) " NATURAL ORBITALS: "
    write(*,*) ""

    do state = 1, CONTROL_instance%CI_STATES_TO_PRINT

      write(*,*) " STATE: ", state

      do species=1, numberOfSpecies

         write(*,*) ""
         write(*,*) " Natural Orbitals in state: ", state, " for: ", trim( MolecularSystem_instance%species(species)%symbol )
         write(*,*) "-----------------"

         numberOfContractions = MolecularSystem_getTotalNumberOfContractions( species )
         speciesName = MolecularSystem_getNameOfSpecies(species)

         call Vector_constructor ( auxdensityEigenValues, &
                             int(numberOfContractions, 8), 0.0_8)

         call Matrix_constructor ( auxdensityEigenVectors, &
              int(numberOfContractions,8), &
              int(numberOfContractions,8),  0.0_8 )

         call Vector_constructor ( densityEigenValues, &
                             int(numberOfContractions, 8), 0.0_8)

         call Matrix_constructor ( densityEigenVectors, &
              int(numberOfContractions,8), &
              int(numberOfContractions,8),  0.0_8 )

         call Matrix_eigen ( ciDensityMatrix(species,state), auxdensityEigenValues, auxdensityEigenVectors, SYMMETRIC )  

         ! reorder and count significant occupations
         k = 0
         do u = 1, numberOfContractions
            densityEigenValues%values(u) =  auxdensityEigenValues%values(numberOfContractions - u + 1)
            densityEigenVectors%values(:,u) = auxdensityEigenVectors%values(:,numberOfContractions - u + 1)
            if(densityEigenValues%values(u) .ge. 5.0E-5 ) k=k+1
         end do

         !! Transform to atomic basis
         densityEigenVectors%values = matmul( coefficients(species)%values, densityEigenVectors%values )

         ! Print eigenvectors with occupation larger than 5.0E-5
         call Matrix_constructor(auxdensityEigenVectors,int(numberOfContractions,8),int(k,8),0.0_8)
         do u = 1, numberOfContractions
            do j = 1, k
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

         !! it's the same as
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

         write(*,*) ""
         write(*,*) " Natural orbital occupation for: ", trim( MolecularSystem_instance%species(species)%symbol )
    
         write(*,*) ""
         do u = 1, numberOfContractions                           
           write(*,"(T2,I4,F17.12)") u, densityEigenValues%values(u)
         enddo
         write(*,*) ""

         write(*,"(A10,A10,A40,F17.12)") "sum of ", trim(MolecularSystem_instance%species(species)%symbol) , "natural orbital occupations", sum(densityEigenValues%values)
         write(*,*) ""

         write(*,*) " End of natural orbitals in state: ", state, " for: ", trim(MolecularSystem_instance%species(species)%symbol)

        call Vector_destructor(auxdensityEigenValues)
        call Matrix_destructor(auxdensityEigenVectors)
        call Vector_destructor(densityEigenValues)
        call Matrix_destructor(densityEigenVectors)

       enddo ! species
    enddo ! state

    write(*,*) ""
    write(*,*) " END OF NATURAL ORBITALS"
    write(6,*) "-----------------------------------------------------------------------"
    write(*,*) ""

    close(unit)

  end subroutine CIdensity_naturalOrbitals

end module CIdensity_
