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
!! @brief Module to handle integral calculations
!! @author Edwin Fernando Posada
!! @version 1.0
!! <b> Fecha de creacion : </b> 2011-02-15
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 02-15-11 </tt>:  E. F. Posada ( efposadac@unal.edu.co )
!!        -# Creacion del modulo y metodos basado en APMO para su inclusion en Lowdin
!!
module IntegralManager_
  use CONTROL_
  use String_
  use MolecularSystem_
  use ContractedGaussian_
  use OverlapIntegrals_
  use ThreeCOverlapIntegrals_
  use AttractionIntegrals_
  use MomentIntegrals_
  use KineticIntegrals_
  use FirstDerivativeIntegrals_
  use HarmonicIntegrals_
  use Libint2Interface_
  ! use CudintInterface_
  use RysQuadrature_
  use Matrix_
  use CosmoCore_
  use Stopwatch_
  use ExternalPotential_


  implicit none

  public :: &
       IntegralManager_getOverlapIntegrals, &
       IntegralManager_getKineticIntegrals, &
       IntegralManager_getFirstDerivativeIntegrals, &
       IntegralManager_getAttractionIntegrals, &
       IntegralManager_getMomentIntegrals, &
       IntegralManager_getInterRepulsionIntegrals, &
       IntegralManager_getIntraRepulsionIntegrals
  private :: &
       IntegralManager_getLabels

contains

  !> 
  !! @brief Calculate overlap integrals and write it on file as a matrix
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_getOverlapIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    ! real(8) :: maxOverlap
    character(100) :: job, colNum

    job = "OVERLAP"

    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       ! maxOverlap=0.0
       
       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of Overlap integrals for species ", &
          ! trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if

       ! write (*, "(A,A,F15.6)"), "Maximum overlap value for ", trim( MolecularSystem_getNameOfSpecie(f)), maxOverlap 
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

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

  end subroutine IntegralManager_getOverlapIntegrals

  !> 
  !! @brief Calculate kinetic energy integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_getKineticIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(100) :: job

    job = "KINETIC"

    !!Kinetic Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of Kinetic integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

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

  end subroutine IntegralManager_getKineticIntegrals

  subroutine IntegralManager_getFirstDerivativeIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(100) :: job
    integer :: ijob

    job = "FIRSTDX"
    ijob = 0 

    !!First derivative Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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

                   !!Calculating integrals for shell
                   call FirstDerivativeIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue, ijob)

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

       write(*,"(A, A ,A,I6)")" Number of First derivative dx integrals for species ", &
                              trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2

       !!Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of First derivative integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done! 

    job = "FIRSTDY"
    ijob = 1

    !!First derivative Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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

                   !!Calculating integrals for shell
                   call FirstDerivativeIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue, ijob)

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

       write(*,"(A, A ,A,I6)")" Number of First derivative dy integrals for species ", &
                              trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2

       !!Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of First derivative integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done! 

    job = "FIRSTDZ"
    ijob = 2 

    !!First derivative Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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

                   !!Calculating integrals for shell
                   call FirstDerivativeIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue, ijob)

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

       write(*,"(A, A ,A,I6)")" Number of First derivative dz integrals for species ", &
                              trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2

       !!Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of First derivative integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done! 

  end subroutine IntegralManager_getFirstDerivativeIntegrals

  subroutine IntegralManager_getHarmonicIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(100) :: job
    integer :: ijob

    job = "HARMONIC"
    ijob = 0 

    !!First derivative Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

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

                   !!Calculating integrals for shell
                   call HarmonicIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
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

       write(*,"(A, A ,A,I6)")" Number of Harmonic Oscillator integrals for species ", &
                              trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2

       !!Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of First derivative integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done! 

  end subroutine IntegralManager_getHarmonicIntegrals


  !> 
  !! @brief Calculate point charge - quantum particle attraction integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: reads point charge information from lowdin.bas file
  !!			- 2014.16.09: modify subroutines to calcule cosmo monoelectronic integrals Danilo
  subroutine IntegralManager_getAttractionIntegrals(surface)
    implicit none

    integer :: f, g, h, i, c
    integer :: j, k, l, m
    integer :: ii, jj, hh, ff
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
    type(surfaceSegment),intent(in), optional :: surface
    ! integer :: times
    integer :: a
    logical :: isCosmo
    real(8), allocatable :: integralValueCosmo(:,:)
    real(8), allocatable ::  cosmoV(:)
    real(8), allocatable ::  qCharges(:)
    ! real(8), allocatable ::  qChargesAc(:,:)
    character(100) :: cosmoIntegralFile
    character(100) :: cosmoQuantumChargeFile
    character(100) :: cosmoClasicalChargeFile


    type(Matrix) :: cmatin
    ! real(8) ::  cosmoPotInt

    integer :: b,d,total_aux

    integer, allocatable ::  totals(:)


    a=0

    job = "ATTRACTION"

    if(present(surface)) then
       job= "COSMO"
       ! en lugar de pasar todo el arreglo de cargas puntuales, debo pasar cada
       ! una de las cargas para así no alterar tanto el código
       isCosmo=.true.
       ! write(*,'(A)')"remplazando los valores de point charges por surface"
       numberOfPointCharges=surface%sizeSurface
       numberOfClasicalCharges = MolecularSystem_instance%numberOfPointCharges
       ! write (*,*) "cargas puntuales cosmo", numberOfPointCharges
       if(allocated(point)) deallocate(point)
       allocate(point(1:numberOfPointCharges))
       point%charge = 0.0_8
       ! write(*,*) "remplazadas por estas"

       if(allocated(totals)) deallocate(totals)
       allocate(totals(size(MolecularSystem_instance%species)))

       write(40) job

       call Matrix_constructor(cmatin, int(surface%sizeSurface,8), int(surface%sizeSurface,8))

       call CosmoCore_cmat(surface,cmatin)

       ! call Matrix_Show(cmatin)
       ! do i=1,120
       ! 	do j=i,120
       ! 		cmatin%values(j,i)=cmatin%values(i,j)
       ! 	end do
       ! end do

       !!do sobre las especies

       do f = 1, size(MolecularSystem_instance%species)
          write(40) MolecularSystem_instance%species(f)%name

          total_aux=0

          cosmoIntegralFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".opints"
          cosmoQuantumChargeFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".charges"

          open(unit=70, file=trim(cosmoIntegralFile), status="unknown",form="unformatted")
          open(unit=80, file=trim(cosmoQuantumChargeFile), status="unknown",form="unformatted")

          if(allocated(labels)) deallocate(labels)
          allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
          labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

          if(allocated(integralBuffer)) deallocate(integralBuffer)
          allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )


          integralBuffer = 0.0_8
          write(80)numberOfPointCharges

          ii = 0
          !! do sobre las particulas
          do g = 1, size(MolecularSystem_instance%species(f)%particles)
             !! do sobre las contracciones 
             do h = 1, size(MolecularSystem_instance%species(f)%particles(g)%basis%contraction)
                hh = h
                ii = ii + 1
                jj = ii - 1
                !! do sobre las particulas diferentes a g de la misma especie
                do i = g, size(MolecularSystem_instance%species(f)%particles)
                   !! do sobre las bases de las particulas diferentes a g de la misma especie
                   do j = hh, size(MolecularSystem_instance%species(f)%particles(i)%basis%contraction)
                      jj = jj + 1
                      if(allocated(integralValue)) deallocate(integralValue)
                      allocate(integralValue(MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                           MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital))

                      if(allocated(integralValueCosmo)) deallocate(integralValueCosmo)
                      allocate(integralValueCosmo((MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital * &
                           MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital), numberOfPointCharges))

                      integralValueCosmo = 0
                      b=MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital
                      d=MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital
                      totals(f)=b*d

                      total_aux=total_aux+totals(f)

                      !Calculating integrals for shell
                      do c = 1, numberOfPointCharges
                         ! do sobre las cargas puntuales 
                         point(1)%charge = 1.0 
                         point(1)%x  =surface%xs(c)
                         point(1)%y  =surface%ys(c)
                         point(1)%z  =surface%zs(c)

                         !Calculating integrals for shell
                         call AttractionIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                              MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), point, 1, integralValue)
                         m=0

                         do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                            do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                               m = m + 1

                               integralValueCosmo(m,c)=integralValue(m)*(-MolecularSystem_getCharge(f))

                            end do
                         end do

                         ! todas las integrales para un c específico 
                      end do
                      ! aqui es donde se hace el calculo para cada carga


                      if(allocated(qCharges)) deallocate(qCharges)
                      allocate(qCharges(numberOfPointCharges))


                      m=0
                      do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                         do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                            m = m + 1
                            if(allocated(cosmoV)) deallocate(cosmoV)
                            allocate(cosmoV(numberOfPointCharges))
                            cosmoV = 0 
                            cosmoV(:)=integralValueCosmo(m,:)

                            call CosmoCore_q_builder(cmatin, cosmoV, numberOfPointCharges, qCharges,f)

                            write(70)integralValueCosmo(m,:)
                            write(80)qCharges
                         end do
                      end do

                   end do
                   !! end do de bases
                   hh = 1
                end do
                !!end do particulas
             end do
             ! end do bases
          end do
          !! end do particulas

          close(80)
          close(70)

          !!quantum
          totals(f)=total_aux
          ! ####################################
          ! Nuevo, ojo
          ! write(80)numberOfPointCharges
          ! write(80)totals(f)

          ! ####################################
          call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoQuantumChargeFile,numberOfPointCharges,totals(f),totals(f),f,f)


          !!clasical vs clasical

          cosmoClasicalChargeFile="cosmo.clasical"

          job="COSMO1"

          write(40) job
          write(40) MolecularSystem_instance%species(f)%name

          call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoClasicalChargeFile,numberOfPointCharges,1,totals(f),f,f,labels)
          !clasical vs quantum


          job="COSMO4"
          write(40) job
          write(40) MolecularSystem_instance%species(f)%name

          call CosmoCore_nucleiPotentialQuantumCharges(surface,cosmoQuantumChargeFile,totals(f),labels,f)

       end do
       !! end do especies


       if( MolecularSystem_getNumberOfQuantumSpecies() > 1 ) then


          do f = 1, size(MolecularSystem_instance%species)
             do g= 1, size(MolecularSystem_instance%species)

                !! do all possible combinations of potentials and charges.

                if ( f /= g ) then

                   cosmoQuantumChargeFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( f ) )//".charges"
                   cosmoIntegralFile="cosmo"//trim( MolecularSystem_getNameOfSpecie( g ) )//".opints"

                   call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoQuantumChargeFile,numberOfPointCharges,totals(f),totals(g),f,g)

                end if
             end do
          end do

       end if


    else

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

          write(30) job
          write(30) MolecularSystem_instance%species(f)%name

         !do ff = 0, numberOfPointCharges - 1
         !   point(ff)%charge = - MolecularSystem_instance%pointCharges(ff+1)%charge * MolecularSystem_getCharge( f )
         !end do

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
          if(CONTROL_instance%LAST_STEP) then
             ! write(*,"(A, A ,A,I6)")" Number of Nuclear integrals for species ", &
             !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
          end if
          write(30) int(size(integralsMatrix),8)
          write(30) integralsMatrix



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
    end if

  end subroutine IntegralManager_getAttractionIntegrals

  !> 
  !! @brief Calculate moment integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_getMomentIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer :: component
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralBuffer(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(10) :: coordinate(3)
    character(100) :: job

    job = "MOMENT"

    coordinate = ["X", "Y", "Z"]

    !!Moment Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralBuffer)) deallocate(integralBuffer)
       allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )
       integralBuffer = 0.0_8

       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = f), MolecularSystem_getTotalNumberOfContractions(specieID = f)))

       do component = 1, 3 !! components x, y, z

          write(30) trim(job)//trim(coordinate(component))
          write(30) MolecularSystem_instance%species(f)%name
          !write(30) coordinate(component)

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
                      call MomentIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                           MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), [0.0_8, 0.0_8, 0.0_8], component, integralValue)

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
          write(30) int(size(integralsMatrix),8)          
          write(30) integralsMatrix

          ! !!Depuration block
          ! print*, "Moment Matrix for specie: ", f, " and component ", component

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

       end do !! Done by component

       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of Moment integrals for species ", &
          !    trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if

    end do !done! 

  end subroutine IntegralManager_getMomentIntegrals

  !> 
  !! @brief Calculate Intra-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  !!      - 2016.06.13: Use Libint V 2.x
  subroutine IntegralManager_getIntraRepulsionIntegrals(nameOfSpecies, scheme)
    implicit none

    character(*) :: nameOfSpecies
    character(*) :: scheme

    integer :: speciesID
    integer :: numberOfContractions

    !! Skip integrals calculation two times for electrons alpha and beta    
    if(CONTROL_instance%IS_OPEN_SHELL .and. ( trim(nameOfSpecies) == "E-BETA" )) return

    speciesID = MolecularSystem_getSpecieID(trim(nameOfSpecies))
    numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)

    if ( trim(String_getUppercase( CONTROL_instance%INTEGRAL_STORAGE )) == "DIRECT") return 

    !! Calculate integrals (stored on disk)           
    select case (trim(String_getUppercase(trim(scheme))))

    case("RYS")

      !! Check critical OMP
       call RysQuadrature_computeIntraSpecies( speciesID )

    case("LIBINT")

       call Libint2Interface_compute2BodyIntraspecies_disk(speciesID)

    case default

       call Libint2Interface_compute2BodyIntraspecies_disk(speciesID)

    end select

    ! call IntegralManager_SaveNumberOfNonZeroIntegrals(speciesID, auxCounter, process ) 

  end subroutine IntegralManager_getIntraRepulsionIntegrals


  !> 
  !! @brief Calculate Inter-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  !!      - 2016.06.13: Use Libint V 2.x
  subroutine IntegralManager_getInterRepulsionIntegrals(scheme)
    implicit none
    character(*) :: scheme    
    integer :: i, j
    integer :: auxCounter = 0

    if ( trim(String_getUppercase( CONTROL_instance%INTEGRAL_STORAGE )) == "DIRECT") return 

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       do j = i+1, MolecularSystem_instance%numberOfQuantumSpecies

          !! Calculate integrals (stored on disk)       
          select case (trim(String_getUppercase(trim(scheme))))

          case("LIBINT")

             call Libint2Interface_compute2BodyInterspecies_disk(i, j)
             ! call LibintInterface_computeInterSpecies( i, j, "ERIS", auxCounter )
             ! case("CUDINT")
             !    call CudintInterface_computeInterSpecies( i, j, "ERIS" )

          case default

             call Libint2Interface_compute2BodyInterspecies_disk(i, j)
             ! call LibintInterface_computeInterSpecies( i, j, "ERIS", auxCounter )

          end select

          call IntegralManager_SaveNumberOfNonZeroCouplingIntegrals(i, j, auxCounter) 

       end do
    end do

    close(30) 

  end subroutine IntegralManager_getInterRepulsionIntegrals

  !> 
  !! @brief Calculate three center integrals and write it on file as a matrix
  !! @version 1.0
  subroutine IntegralManager_getThreeCenterIntegrals()
    implicit none

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

       write(30) job
       speciesName = MolecularSystem_instance%species(f)%name
       write(30) speciesName

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

                         integralsMatrix(k, l) = integralValue(m)
                         integralsMatrix(l, k) = integralsMatrix(k, l)

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do

       !! Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
           write(*,"(A, A ,A,I6)")" Number of 3C Overlap integrals for species ", &
           trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done!    

  end subroutine IntegralManager_getThreeCenterIntegrals

  !> 
  !! @brief Calculate three center integrals and write it on file as a matrix
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_getThreeCenterIntegralsByProduct()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
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

       write(30) job
       speciesName = MolecularSystem_instance%species(f)%name
       write(30) speciesName

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralsMatrix)) deallocate(integralsMatrix)
       allocate(integralsMatrix(MolecularSystem_getTotalNumberOfContractions(specieID = f), MolecularSystem_getTotalNumberOfContractions(specieID = f)))

       integralsMatrix = 0.0_8

       do i= 1, ExternalPotential_instance%ssize
         !if( trim(potential(i)%specie)==trim(interactNameSelected) ) then ! This does not work for UHF
         ! if ( String_findSubstring(trim( MolecularSystem_instance%species(f)%name  ), &
         !      trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%specie)))) == 1 ) then
         if ( trim( MolecularSystem_instance%species(f)%symbol) == trim(String_getUpperCase(trim(ExternalPotential_instance%potentials(i)%specie))) ) then
           potID=i
           exit
         end if
       end do

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

                   call ContractedGaussian_copyConstructor ( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                        auxContractionA )
                   call ContractedGaussian_copyConstructor ( MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), &
                        auxContractionB )

                   do p = 1, ExternalPotential_instance%potentials(potID)%numOfComponents

                     auxintegralValue = 0.0_8

                     !! Calculating integrals for shell
  
                     call ContractedGaussian_product(auxContractionA, &
                        ExternalPotential_instance%potentials(potID)%gaussianComponents(p), auxContraction)

                     call OverlapIntegrals_computeShell( auxContraction, auxContractionB, auxintegralValue)
                     !call OverlapIntegrals_computeShell( auxContractionA, auxContractionB, auxintegralValue)
                     !call OverlapIntegrals_computeShell ( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                     !  MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), auxintegralValue )

                    integralValue = integralValue + auxintegralValue

                   end do !! potential
                   !! saving integrals on Matrix
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

       !! Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          ! write(*,"(A, A ,A,I6)")" Number of Overlap integrals for species ", &
          ! trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done!    

  end subroutine IntegralManager_getThreeCenterIntegralsByProduct

  !>


  !>
  !! @brief Return real labels for integrals
  !! @autor E. F. Posada, 2011
  !! @version 1.0
  function IntegralManager_getLabels(specieSelected) result(labelsOfContractions)
    implicit none

    type(species) :: specieSelected
    integer:: labelsOfContractions(specieSelected%basisSetSize)

    integer:: auxLabelsOfContractions
    integer:: i, j, k

    auxLabelsOfContractions = 1

    ! write(*,*)"labels data from integral manager"

    k = 0
    ! write(*,*)size(specieSelected%particles)
    do i = 1, size(specieSelected%particles)

       do j = 1, size(specieSelected%particles(i)%basis%contraction)

          k = k + 1

          !!position for cartesian contractions
          labelsOfContractions(k) = auxLabelsOfContractions
          auxLabelsOfContractions = auxLabelsOfContractions + specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital
          ! write(*,*)specieSelected%particles(i)%basis%contraction(j)%numCartesianOrbital

       end do
    end do


  end function IntegralManager_getLabels

  subroutine IntegralManager_SaveNumberOfNonZeroIntegrals(speciesID, auxCounter, process) 
    implicit none

    integer :: speciesID
    integer :: auxCounter
    integer(8) :: process
    character(50) :: fileNumber

    write(fileNumber,*) process
    fileNumber = trim(adjustl(fileNumber))

    open(UNIT=49,FILE=trim(fileNumber)//trim(MolecularSystem_instance%species(speciesID)%name)//".nints", &
         STATUS='replace', ACCESS='SEQUENTIAL', FORM='Unformatted')
    
    write (49) auxCounter 

    close(49)

  end subroutine IntegralManager_SaveNumberOfNonZeroIntegrals

  subroutine IntegralManager_SaveNumberOfNonZeroCouplingIntegrals( i, j, auxCounter) 
    implicit none

    integer :: i, j 
    integer :: auxCounter

    !! open file for integrals
    open(UNIT=59,FILE=trim(MolecularSystem_instance%species(i)%name)//"."//trim(MolecularSystem_instance%species(j)%name)//".nints", &
         STATUS='replace', ACCESS='SEQUENTIAL', FORM='Unformatted')

    write (59) auxCounter 
    close(59)


  end subroutine IntegralManager_SaveNumberOfNonZeroCouplingIntegrals


end module IntegralManager_
