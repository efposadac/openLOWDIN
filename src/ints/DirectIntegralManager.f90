!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Module to handle direct integral calculations
!! @author Jorge Charry
!! @version 1.0
!! <b> Fecha de creacion : </b> 2015-12-28
!!
!! <b> Historial de modificaciones: </b> Merge DirectIntegralManagers from different programs Felix 2022-4
!!
module DirectIntegralManager_
  use CONTROL_  
  use MolecularSystem_
  use IntegralManager_
  use OverlapIntegrals_
  use ThreeCOverlapIntegrals_
  use AttractionIntegrals_
  use KineticIntegrals_
  use Libint2Interface_
  use RysQuadrature_
  use Matrix_
  use Stopwatch_
  use ExternalPotential_
  !# use RysQInts_  !! Please do not remove this line

  implicit none

  public :: &
       DirectIntegralManager_getDirectIntraRepulsionMatrix, &
       DirectIntegralManager_getDirectInterRepulsionMatrix, &
       DirectIntegralManager_getDirectIntraRepulsionFirstQuarter, &
       DirectIntegralManager_getDirectInterRepulsionFirstQuarter, &
       DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix, &
       DirectIntegralManager_destructor, &
       DirectIntegralManager_getOverlapIntegrals, &
       DirectIntegralManager_getKineticIntegrals, &
       DirectIntegralManager_getAttractionIntegrals, &
       DirectIntegralManager_getExternalPotentialIntegrals
contains

  !> 
  !! @brief Calculate Intra-species repulsion integrals directly and use it to obtain two particles matrix
  !! @author J. A. Charry, 2015
  !! @version 1.0
  !! @par History
  !!    
  recursive subroutine DirectIntegralManager_getDirectIntraRepulsionMatrix(speciesID, scheme, &
       densityMatrix, twoParticlesMatrix, factor )
    implicit none

    integer :: speciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: twoParticlesMatrix(:,:)
    real(8) :: factor

    ! integer :: numberOfContractions
    ! integer(8) :: integralsByProcess
    ! integer(8) :: nprocess
    ! integer(8) :: process
    ! integer(8) :: starting
    ! integer(8) :: ending

    real(8), allocatable, target :: density(:,:)
    integer(8) :: ssize

    ssize = size(densityMatrix%values, DIM=1)
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    !     numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)

    !     ssize = (numberOfContractions * (numberOfContractions + 1))/2
    !     ssize = (ssize * (ssize + 1))/2

    !     integralsByProcess = ceiling( real(ssize,8)/real(nprocess,8) )

    !     ending = process * integralsByProcess
    !     starting = ending - integralsByProcess + 1

    !     if( starting > ssize ) return

    !     if( ending > ssize ) ending = ssize

    !! Calculate integrals
    select case (trim(String_getUppercase(trim(scheme))))

       !     case("RYS")
       !        call RysQuadrature_directIntraSpecies( speciesID, "ERIS", starting, ending, int( process ) , &
       !               densityMatrix, & 
       !               twoParticlesMatrix, factor)
    case("LIBINT")
       call Libint2Interface_compute2BodyIntraspecies_direct(speciesID, density, twoParticlesMatrix, factor )

       !     ! case("CUDINT")
       !     !    call CudintInterface_computeIntraSpecies(speciesID)
    case default
       call Libint2Interface_compute2BodyIntraspecies_direct(speciesID, density, twoParticlesMatrix, factor )
    end select

    deallocate(density)
  end subroutine DirectIntegralManager_getDirectIntraRepulsionMatrix


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly and use it to obtain coupling matrix
  !! @author E. F. Posada 2016
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionMatrix(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, couplingMatrix )
    implicit none
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    real(8), allocatable, target :: couplingMatrix(:,:)

    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    select case (trim(String_getUppercase(trim(scheme))))

       !case("RYS")
       ! Not implemented

    case("LIBINT")
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix)
    case default
       call Libint2Interface_compute2BodyInterspecies_direct(speciesID, otherSpeciesID, density, couplingMatrix)
    end select

    deallocate(density)

  end subroutine DirectIntegralManager_getDirectInterRepulsionMatrix

  !> 
  !! @brief Calculate Intra-species repulsion integrals directly and transform the first index from atomic to molecular orbitals
  !! @author J. A. Charry, 2015
  !! @version 1.0
  !! @par History
  !!    
  recursive subroutine DirectIntegralManager_getDirectIntraRepulsionFirstQuarter(speciesID, scheme, &
       densityMatrix, coeffMatrix, matrixA, p )
    implicit none

    integer :: speciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: coeffMatrix
    real(8), allocatable, target :: matrixA(:,:,:)
    integer :: p

    integer :: numberOfContractions

    real(8), allocatable, target :: coefficients(:,:)
    real(8), allocatable, target :: density(:,:)
    integer :: ssize

    ssize = size(coeffMatrix%values, DIM=1)
    allocate(coefficients(ssize, ssize))
    coefficients = coeffMatrix%values

    allocate(density(ssize, ssize))
    density = densityMatrix%values

    select case (trim(String_getUppercase(trim(scheme))))

       !     case("RYS")
       !        call RysQuadrature_directIntraSpecies( speciesID, "ERIS", starting, ending, int( process ) , &
       !               densityMatrix, & 
       !               twoParticlesMatrix, factor)
    case("LIBINT")
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p )

       !     ! case("CUDINT")
       !     !    call CudintInterface_computeIntraSpecies(speciesID)
    case default
       call Libint2Interface_compute2BodyIntraspecies_direct_IT(speciesID, density, coefficients, matrixA, p )
    end select

    deallocate(coefficients,density)

  end subroutine DirectIntegralManager_getDirectIntraRepulsionFirstQuarter


  !> 
  !! @brief Calculate Inter-species repulsion integrals directly and transform the first index from atomic to molecular orbitals
  !! @author E. F. Posada 2016
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_getDirectInterRepulsionFirstQuarter(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, coeffMatrix, couplingMatrix, p )
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: densityMatrix
    type(matrix) :: coeffMatrix
    real(8), allocatable, target :: couplingMatrix(:,:,:)

    real(8), allocatable, target :: coefficients(:,:)
    real(8), allocatable, target :: density(:,:)
    integer :: ssize
    integer :: p

    ssize = size(coeffMatrix%values, DIM=1)

    allocate(coefficients(ssize, ssize))
    coefficients = coeffMatrix%values
    
    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    select case (trim(String_getUppercase(trim(scheme))))

       !case("RYS")
       ! Not implemented

    case("LIBINT")
       call Libint2Interface_compute2BodyInterspecies_direct_IT(speciesID, otherSpeciesID, density, coefficients, couplingMatrix, p)
    case default
       call Libint2Interface_compute2BodyInterspecies_direct_IT(speciesID, otherSpeciesID, density, coefficients, couplingMatrix, p)
    end select

    deallocate(coefficients,density)
    
  end subroutine DirectIntegralManager_getDirectInterRepulsionFirstQuarter

  subroutine DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix(speciesID, OtherSpeciesID, scheme, &
       densityMatrix, otherdensityMatrix, coupling )
    integer :: speciesID
    integer :: otherSpeciesID
    character(*) :: scheme
    type(matrix) :: couplingMatrix
    type(matrix) :: densityMatrix
    type(matrix) :: otherdensityMatrix
    real(8), allocatable, target :: coupling(:,:)
    real(8), allocatable, target :: density(:,:)
    real(8), allocatable, target :: otherdensity(:,:)
    integer :: ssize

    ssize = size(densityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    if(allocated(density))deallocate(density)
    allocate(density(ssize, ssize))
    density = densityMatrix%values

    ssize = size(otherdensityMatrix%values, DIM=1)
    ! print*, "DIRECT, SIZE DENS:", ssize
    if(allocated(otherdensity))deallocate(otherdensity)
    allocate(otherdensity(ssize, ssize))
    otherdensity = otherdensityMatrix%values


    select case (trim(String_getUppercase(trim(scheme))))

       !case("RYS")
       ! Not implemented

    case("LIBINT")
       call Libint2Interface_compute2BodyAlphaBeta_direct(speciesID, otherSpeciesID, density, otherdensity, coupling)
    case default
       call Libint2Interface_compute2BodyAlphaBeta_direct(speciesID, otherSpeciesID, density, otherdensity, coupling)
    end select


    deallocate(density)

  end subroutine DirectIntegralManager_getDirectAlphaBetaRepulsionMatrix

  !> 
  !! @brief Destroy libint interface objects, to update the molecular system  
  !! @author F. Moncada 2022
  !! @version 1.0
  !! @par History
  !!    
  subroutine DirectIntegralManager_destructor()
    implicit none
    integer speciesID

    do speciesID=1, size(Libint2Instance(:))
       call Libint2Interface_destructor(Libint2Instance(speciesID))
    end do
  end subroutine DirectIntegralManager_destructor

  !> 
  !! @brief Calculate overlap integrals and return them as matrices
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  subroutine DirectIntegralManager_getOverlapIntegrals(output)
    implicit none
    type(Matrix), intent(out) :: output(:)

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    ! real(8) :: maxOverlap
    character(100) :: colNum

    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       ! maxOverlap=0.0

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = f), MolecularSystem_getTotalNumberOfContractions(specieID = f)))

       integralsMatrix = 0.0_8

       ii = 0
       do g = 1, size(MolecularSystem_instance%species(f)%particles)
          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)

             hh = h

             ii = ii + 1
             jj = ii - 1

             do i = g, size(MolecularSystem_instance%species(f)%particles)
                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)

                   jj = jj + 1

                   !! allocating memory Integrals for shell
                   if(allocated(integralValue)) deallocate(integralValue)
                   allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                   integralValue = 0.0_8

                   !! Calculating integrals for shell
                   call OverlapIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue)

                   !! saving integrals on Matrix
                   m = 0
                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                         m = m + 1

                         integralsMatrix(k, l) = integralValue(m)
                         integralsMatrix(l, k) = integralsMatrix(k, l)

                         ! if( (integralValue(m) .gt. maxOverlap) .and. (l .ne. k) ) maxOverlap=integralValue(m)

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do

       !! Write integrals to file (unit 30)
       ! if(CONTROL_instance%LAST_STEP) then
       ! write(*,"(A, A ,A,I6)")" Number of Overlap integrals for species ", &
       ! trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       ! end if

       ! write (*, "(A,A,F15.6)"), "Maximum overlap value for ", trim( MolecularSystem_getNameOfSpecie(f)), maxOverlap 
       call Matrix_constructor(output(f),int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), &
            int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), 0.0_8 )

       output(f)%values=integralsMatrix

       !!Depuration block
       ! print*, "Overlap Matrix for species: ", f

       ! do  k = 1, ceiling((size( integralsMatrix, DIM=2 )/5.0_8))

       !    l = 5 * (k-1)+1
       !    h = 5 * k
       !    g = 5

       !    if(h > size( integralsMatrix, DIM=2 )) then
       !       g = 5 - h + size( integralsMatrix, DIM=2 )
       !       h = size( integralsMatrix, DIM=2 )
       !    end if

       !    write(colNum,*) g

       !    write (*,"(5X,"//trim(colNum)//"F15.6)")((integralsMatrix(i,j),j=l,h),i=1,size(integralsMatrix,DIM=1))

       !    print*, ""
       !    print*, ""

       ! end do

    end do !done!    

  end subroutine DirectIntegralManager_getOverlapIntegrals

  !> 
  !! @brief Calculate kinetic energy integrals and return them as a matrix
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  subroutine DirectIntegralManager_getKineticIntegrals(output)
    implicit none
    type(Matrix), intent(out) :: output(:)

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)

    !!Kinetic Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = f), MolecularSystem_getTotalNumberOfContractions(specieID = f)))

       ii = 0
       do g = 1, size(MolecularSystem_instance%species(f)%particles)
          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)

             hh = h

             ii = ii + 1
             jj = ii - 1

             do i = g, size(MolecularSystem_instance%species(f)%particles)
                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)

                   jj = jj + 1

                   !! allocating memory Integrals for shell
                   if(allocated(integralValue)) deallocate(integralValue)
                   allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                   !!Calculating integrals for shell
                   call KineticIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue)

                   !!saving integrals on Matrix
                   m = 0
                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                         m = m + 1

                         integralsMatrix(k, l) = integralValue(m)
                         integralsMatrix(l, k) = integralsMatrix(k, l)

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do

       !!Write integrals to file (unit 30)
       ! if(CONTROL_instance%LAST_STEP) then
       !    ! write(*,"(A, A ,A,I6)")" Number of Kinetic integrals for species ", &
       !    !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       ! end if

       call Matrix_constructor(output(f),int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), &
            int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), 0.0_8 )

       output(f)%values=integralsMatrix

       ! !!Depuration block
       ! print*, "Kinetic Matrix for specie: ", f

       ! do  k = 1, ceiling((size( integralsMatrix, DIM=2 )/5.0_8))

       !    l = 5 * (k-1)+1
       !    h = 5 * k
       !    g = 5

       !    if(h > size( integralsMatrix, DIM=2 )) then
       !       g = 5 - h + size( integralsMatrix, DIM=2 )
       !       h = size( integralsMatrix, DIM=2 )
       !    end if

       !    write(colNum,*) g

       !    write (*,"(5X,"//trim(colNum)//"F15.6)")((integralsMatrix(i,j),j=l,h),i=1,size(integralsMatrix,DIM=1))

       !    print*, ""
       !    print*, ""

       ! end do

    end do !done! 

  end subroutine DirectIntegralManager_getKineticIntegrals


  !> 
  !! @brief Calculate point charge - quantum particle attraction integrals, return them as a matrix
  !! @author E. F. Posada, 2013 - F. Moncada, 2022
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: reads point charge information from lowdin.bas file
  !!			- 2014.16.09: modify subroutines to calcule cosmo monoelectronic integrals Danilo
  subroutine DirectIntegralManager_getAttractionIntegrals(output) !, surface
    implicit none
    type(Matrix), intent(out) :: output(:)
    ! type(surfaceSegment),intent(in), optional :: surface
    !
    integer :: f, g, h, i, c
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer :: numberOfPointCharges
    integer :: numberOfClasicalCharges
    ! integer :: numberOfSurfaceSegments
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralBuffer(:)
    real(8), allocatable :: integralsMatrix(:,:)
    type(pointCharge), allocatable :: point(:)
    ! character(20) :: colNum
    character(100) :: job

    ! Variables para calcular las integrales monoelectronicas para cosmo
    ! integer :: times
    ! integer :: a
    ! logical :: isCosmo
    ! real(8), allocatable :: integralValueCosmo(:,:)
    ! real(8), allocatable ::  cosmoV(:)
    ! real(8), allocatable ::  qCharges(:)
    ! ! real(8), allocatable ::  qChargesAc(:,:)
    ! character(100) :: cosmoIntegralFile
    ! character(100) :: cosmoQuantumChargeFile
    ! character(100) :: cosmoClasicalChargeFile
    ! type(Matrix) :: cmatin
    ! real(8) ::  cosmoPotInt
    ! integer :: b,d,total_aux

    ! integer, allocatable ::  totals(:)

    ! a=0

    ! if(present(surface)) then
    !    job= "COSMO"
    !    ! en lugar de pasar todo el arreglo de cargas puntuales, debo pasar cada
    !    ! una de las cargas para así no alterar tanto el código
    !    isCosmo=.true.
    !    ! write(*,'(A)')"remplazando los valores de point charges por surface"
    !    numberOfPointCharges=surface%sizeSurface
    !    numberOfClasicalCharges = MolecularSystem_instance%numberOfPointCharges
    !    ! write (*,*) "cargas puntuales cosmo", numberOfPointCharges
    !    if(allocated(point)) deallocate(point)
    !    allocate(point(1:numberOfPointCharges))
    !    point%charge = 0.0_8
    !    ! write(*,*) "remplazadas por estas"

    !    if(allocated(totals)) deallocate(totals)
    !    allocate(totals(size(MolecularSystem_instance%species)))

    !    write(40) job

    !    call Matrix_constructor(cmatin, int(surface%sizeSurface,8), int(surface%sizeSurface,8))

    !    call CosmoCore_cmat(surface,cmatin)

    !    ! call Matrix_Show(cmatin)
    !    ! do i=1,120
    !    ! 	do j=i,120
    !    ! 		cmatin%values(j,i)=cmatin%values(i,j)
    !    ! 	end do
    !    ! end do

    !    !!do sobre las especies

    !    do f = 1, size(MolecularSystem_instance%species)
    !       write(40) MolecularSystem_instance%species(f)%name

    !       total_aux=0

    !       cosmoIntegralFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".opints"
    !       cosmoQuantumChargeFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".charges"

    !       open(unit=70, file=trim(cosmoIntegralFile), status="unknown",form="unformatted")
    !       open(unit=80, file=trim(cosmoQuantumChargeFile), status="unknown",form="unformatted")

    !       if(allocated(labels)) deallocate(labels)
    !       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
    !       labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

    !       if(allocated(integralBuffer)) deallocate(integralBuffer)
    !       allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )


    !       integralBuffer = 0.0_8
    !       write(80)numberOfPointCharges

    !       ii = 0
    !       !! do sobre las particulas
    !       do g = 1, size(MolecularSystem_instance%species(f)%particles)
    !          !! do sobre las contracciones 
    !          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)
    !             hh = h
    !             ii = ii + 1
    !             jj = ii - 1
    !             !! do sobre las particulas diferentes a g de la misma especie
    !             do i = g, size(MolecularSystem_instance%species(f)%particles)
    !                !! do sobre las bases de las particulas diferentes a g de la misma especie
    !                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)
    !                   jj = jj + 1
    !                   if(allocated(integralValue)) deallocate(integralValue)
    !                   allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
    !                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))

    !                   if(allocated(integralValueCosmo)) deallocate(integralValueCosmo)
    !                   allocate(integralValueCosmo((MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
    !                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital), numberOfPointCharges))

    !                   integralValueCosmo = 0
    !                   b=MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital
    !                   d=MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital
    !                   totals(f)=b*d

    !                   total_aux=total_aux+totals(f)

    !                   !Calculating integrals for shell
    !                   do c = 1, numberOfPointCharges
    !                      ! do sobre las cargas puntuales 
    !                      point(1)%charge = 1.0 
    !                      point(1)%x  =surface%xs(c)
    !                      point(1)%y  =surface%ys(c)
    !                      point(1)%z  =surface%zs(c)

    !                      !Calculating integrals for shell
    !                      call AttractionIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
    !                           MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), point, 1, integralValue)
    !                      m=0

    !                      do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
    !                         do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
    !                            m = m + 1

    !                            integralValueCosmo(m,c)=integralValue(m)*(-MolecularSystem_getCharge(f))

    !                         end do
    !                      end do

    !                      ! todas las integrales para un c específico 
    !                   end do
    !                   ! aqui es donde se hace el calculo para cada carga


    !                   if(allocated(qCharges)) deallocate(qCharges)
    !                   allocate(qCharges(numberOfPointCharges))


    !                   m=0
    !                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
    !                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
    !                         m = m + 1
    !                         if(allocated(cosmoV)) deallocate(cosmoV)
    !                         allocate(cosmoV(numberOfPointCharges))
    !                         cosmoV = 0 
    !                         cosmoV(:)=integralValueCosmo(m,:)

    !                         call CosmoCore_q_builder(cmatin, cosmoV, numberOfPointCharges, qCharges,f)

    !                         write(70)integralValueCosmo(m,:)
    !                         write(80)qCharges
    !                      end do
    !                   end do

    !                end do
    !                !! end do de bases
    !                hh = 1
    !             end do
    !             !!end do particulas
    !          end do
    !          ! end do bases
    !       end do
    !       !! end do particulas

    !       close(80)
    !       close(70)

    !       !!quantum
    !       totals(f)=total_aux
    !       ! ####################################
    !       ! Nuevo, ojo
    !       ! write(80)numberOfPointCharges
    !       ! write(80)totals(f)

    !       ! ####################################
    !       call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoQuantumChargeFile,numberOfPointCharges,totals(f),totals(f),f,f)


    !       !!clasical vs clasical

    !       cosmoClasicalChargeFile="cosmo.clasical"

    !       job="COSMO1"

    !       write(40) job
    !       write(40) MolecularSystem_instance%species(f)%name

    !       call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoClasicalChargeFile,numberOfPointCharges,1,totals(f),f,f,labels)

    !       !clasical vs quantum


    !       job="COSMO4"
    !       write(40) job
    !       write(40) MolecularSystem_instance%species(f)%name

    !       call CosmoCore_nucleiPotentialQuantumCharges(surface,cosmoQuantumChargeFile,totals(f),labels,f)

    !    end do
    !    !! end do especies


    !    if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then


    !       do f = 1, size(MolecularSystem_instance%species)
    !          do g= 1, size(MolecularSystem_instance%species)

    !             !! do all possible combinations of potentials and charges.

    !             if ( f /= g ) then

    !                cosmoQuantumChargeFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".charges"
    !                cosmoIntegralFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( g ) )//".opints"

    !                call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoQuantumChargeFile,numberOfPointCharges,totals(f),totals(g),f,g)

    !             end if
    !          end do
    !       end do

    !    end if


    ! else

    numberOfPointCharges = MolecularSystem_instance%numberOfPointCharges

    !! Allocating memory for point charges objects
    if(allocated(point)) deallocate(point)
    allocate(point(0:numberOfPointCharges - 1))


    do f = 0, numberOfPointCharges - 1
       point(f)%charge = MolecularSystem_instance%pointCharges(f+1)%charge
       point(f)%x  = MolecularSystem_instance%pointCharges(f+1)%origin(1)
       point(f)%y  = MolecularSystem_instance%pointCharges(f+1)%origin(2)
       point(f)%z  = MolecularSystem_instance%pointCharges(f+1)%origin(3)
    end do

    !!Attraction Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = f), MolecularSystem_getTotalNumberOfContractions(specieID = f)))

       ii = 0
       do g = 1, size(MolecularSystem_instance%species(f)%particles)
          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)

             hh = h
             ii = ii + 1
             jj = ii - 1

             do i = g, size(MolecularSystem_instance%species(f)%particles)
                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)

                   jj = jj + 1

                   !! allocating memory Integrals for shell

                   if(allocated(integralValue)) deallocate(integralValue)
                   allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))


                   !!Calculating integrals for shell
                   call AttractionIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), point, numberOfPointCharges, integralValue)

                   !!saving integrals on Matrix
                   m = 0
                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                         m = m + 1


                         ! write(*,*)"lowdin integrals:f, m,k,l, integral value",f,m,k,l,integralValue(m)
                         integralsMatrix(k, l) = integralValue(m)
                         integralsMatrix(l, k) = integralsMatrix(k, l)

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do


       !!Write integrals to file (unit 30)
       ! if(CONTROL_instance%LAST_STEP) then
       ! write(*,"(A, A ,A,I6)")" Number of Nuclear integrals for species ", &
       !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       ! end if
       call Matrix_constructor(output(f),int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), &
            int(MolecularSystem_getTotalNumberOfContractions(specieID = f),8), 0.0_8 )

       output(f)%values=integralsMatrix


       !!Depuration block
       !        print*, "Attraction  Matrix for specie: ", f

       !        do  k = 1, ceiling((size( integralsMatrix, DIM=2 )/5.0_8))

       !           l = 5 * (k-1)+1
       !           h = 5 * k
       !           g = 5

       !           if(h > size( integralsMatrix, DIM=2 )) then
       !              g = 5 - h + size( integralsMatrix, DIM=2 )
       !              h = size( integralsMatrix, DIM=2 )
       !           end if

       !           write(colNum,*) g

       !           write (*,"(5X,"//trim(colNum)//"F15.6)")((integralsMatrix(i,j),j=l,h),i=1,size(integralsMatrix,DIM=1))

       !           print*, ""
       !           print*, ""

       !        end do

    end do !done! 
    ! end if

  end subroutine DirectIntegralManager_getAttractionIntegrals

  ! !>
  ! !! @brief Return real labels for integrals
  ! !! @autor E. F. Posada, 2011
  ! !! @version 1.0
  ! function DirectIntegralManager_getLabels(specieSelected) result(labelsOfContractions)
  !   implicit none

  !   type(species) :: specieSelected
  !   integer:: labelsOfContractions(specieSelected%basisSetSize)

  !   integer:: auxLabelsOfContractions
  !   integer:: i, j, k

  !   auxLabelsOfContractions = 1

  !   ! write(*,*)"labels data from integral manager"

  !   k = 0
  !   ! write(*,*)size(specieSelected%particles)
  !   do i = 1, size(specieSelected%particles)

  !      do j = 1, size(specieSelected%particles(i)%basis%contraction)

  !         k = k + 1

  !         !!position for cartesian contractions
  !         labelsOfContractions(k) = auxLabelsOfContractions
  !         auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital
  !         ! write(*,*)specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital

  !      end do
  !   end do


  ! end function DirectIntegralManager_getLabels

  subroutine DirectIntegralManager_getExternalPotentialIntegrals(output)
    implicit none
    type(Matrix), intent(out) :: output(:)
    integer :: f, g, h, i
    integer :: j, k, l, m, r
    integer :: o, p, potID
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:), auxintegralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(100) :: job, speciesName
    type (ContractedGaussian) :: auxContractionA, auxContractionB, auxContraction
    type (ExternalPotential) :: potential
    

    job = "EXTERNAL_POTENTIAL"
    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       potID = 0

       do i= 1, ExternalPotential_instance%ssize
         !if( trim(potential(i)%specie)==trim(interactNameSelected) ) then ! This does not work for UHF
         ! if ( String_findSubstring(trim( MolecularSystem_instance%species(f)%name  ), &
         !      trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%specie)))) == 1 ) then

         if ( trim( MolecularSystem_instance%species(f)%symbol) == trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%specie))) ) then
           potID=i
           exit
         end if
       end do

       ! write(30) job
       speciesName = MolecularSystem_instance%species(f)%name
       ! write(30) speciesName

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       call Matrix_constructor(output(f), int(MolecularSystem_getTotalNumberOfContractions(f),8), &
            int(MolecularSystem_getTotalNumberOfContractions(f),8), 0.0_8)

       ii = 0
       do g = 1, size(MolecularSystem_instance%species(f)%particles)
          do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)

             hh = h

             ii = ii + 1
             jj = ii - 1

             do i = g, size(MolecularSystem_instance%species(f)%particles)
                do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)

                   jj = jj + 1

                   !! allocating memory Integrals for shell
                   if(allocated(integralValue)) deallocate(integralValue)
                   allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                   if(allocated(auxintegralValue)) deallocate(auxintegralValue)
                   allocate(auxintegralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))


                   integralValue = 0.0_8
                   do p = 1, ExternalPotential_instance%potentials(potID)%numOfComponents
                     auxintegralValue = 0.0_8

                      
                     do r = 1, ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%numcartesianOrbital
                       ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%primNormalization( &
                       1:ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%length,r) = 1

                     ExternalPotential_instance%potentials(potID)%gaussianComponents(p)%contNormalization(r) = 1
                    end do

                     !! Calculating integrals for shell
                     call ThreeCOverlapIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                           MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), &
                           ExternalPotential_instance%potentials(potID)%gaussianComponents(p), auxintegralValue)
                    integralValue = integralValue + auxintegralValue

                   end do !! potential
                   !! saving integrals on Matrix
                   m = 0
                   do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                      do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                         m = m + 1

                         output(f)%values(k, l) = integralValue(m)
                         output(f)%values(l, k) = output(f)%values(k, l)

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do
    end do !done!    

  end subroutine DirectIntegralManager_getExternalPotentialIntegrals

  
end module DirectIntegralManager_
