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
!!   - <tt> 10-03-23 </tt>:  F. Moncada ( fsmoncadaa@unal.edu.co )
!!        -# Separa el modulo - mueve rutinas a DirectIntegralManager - quedan rutinas que escriben a disco
!!
module IntegralManager_
  use CONTROL_
  use String_
  use MolecularSystem_
  use ContractedGaussian_
  use DirectIntegralManager_
  use Libint2Interface_
  ! use CudintInterface_
  use RysQuadrature_
  use Matrix_
  use CosmoCore_
  use Stopwatch_
  use ExternalPotential_
  use HarmonicIntegrals_
  use FirstDerivativeIntegrals_

  implicit none

  public :: &
       IntegralManager_writeOverlapIntegrals, &
       IntegralManager_writeKineticIntegrals, &
       IntegralManager_writeAttractionIntegrals, &
       IntegralManager_writeMomentIntegrals, &
       IntegralManager_writeInterRepulsionIntegrals, &
       IntegralManager_writeIntraRepulsionIntegrals, &
       IntegralManager_writeHarmonicIntegrals
  ! private :: &

contains

  !> 
  !! @brief Calculate overlap integrals and write it on file as a matrix
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_writeOverlapIntegrals()
    implicit none

    integer :: f
    Type(Matrix) :: integralsMatrix
    character(100) :: arguments(2) 

    arguments(1) = "OVERLAP"

    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       arguments(2) = MolecularSystem_instance%species(f)%name

       call DirectIntegralManager_getOverlapIntegrals(MolecularSystem_instance,f,integralsMatrix)

       call Matrix_writeToFile(integralsMatrix, unit=30, binary=.true., arguments = arguments )  

       !! Write integrals to file (unit 30)
       if(CONTROL_instance%LAST_STEP) then
          write(*,"(A, A ,A,I6)")" Number of Overlap integrals for species ", &
               trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix%values,DIM=1)**2
       end if

       !!Depuration block
       ! print*, "Overlap Matrix for species: ", f
       !Call Matrix_show(integralsMatrix)
       
    end do !done!    

  end subroutine IntegralManager_writeOverlapIntegrals

  !> 
  !! @brief Calculate kinetic energy integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_writeKineticIntegrals()
    implicit none

    integer :: f
    Type(Matrix) :: integralsMatrix
    character(100) :: arguments(2) 

    arguments(1) = "KINETIC"

    !!Kinetic Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       arguments(2) = MolecularSystem_instance%species(f)%name

       call DirectIntegralManager_getKineticIntegrals(MolecularSystem_instance,f,integralsMatrix)
       !!Write integrals to file (unit 30)
       call Matrix_writeToFile(integralsMatrix, unit=30, binary=.true., arguments = arguments )
       
       if(CONTROL_instance%LAST_STEP) then
          write(*,"(A, A ,A,I6)")" Number of Kinetic integrals for species ", &
               trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix%values,DIM=1)**2
       end if

       !!Depuration block
       ! print*, "Kinetic Matrix for specie: ", f
       !Call Matrix_show(integralsMatrix)

    end do !done! 

  end subroutine IntegralManager_writeKineticIntegrals

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
       labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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
       labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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
       labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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

  subroutine IntegralManager_writeHarmonicIntegrals()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    real(8) :: origin(3)
    character(100) :: job
    integer :: ijob

    job = "HARMONIC"
    ijob = 0 

    !!First derivative Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

      if ( MolecularSystem_getOmega(f) /= 0.0_8 ) then
        origin = MolecularSystem_getQDOcenter( f ) 

         write(30) job
         write(30) MolecularSystem_instance%species(f)%name

         if(allocated(labels)) deallocate(labels)
         allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
         labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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
                          MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), integralValue, origin)

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
      end if
    end do !done! 

  end subroutine IntegralManager_writeHarmonicIntegrals


  !> 
  !! @brief Calculate point charge - quantum particle attraction integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: reads point charge information from lowdin.bas file
  !!			- 2014.16.09: modify subroutines to calcule cosmo monoelectronic integrals Danilo
  subroutine IntegralManager_writeAttractionIntegrals(surface)
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
    type(pointCharge), allocatable :: point(:)
    type(Matrix) :: integralsMatrix
    ! character(20) :: colNum
    character(100) :: job, arguments(2)

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
       !  do j=i,120
       !  	cmatin%values(j,i)=cmatin%values(i,j)
       !  end do
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
          labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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
                         point(1)%qdoCenterOf = "NONE"

                         !Calculating integrals for shell
                         call AttractionIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                              MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), point, 1, integralValue, f, "NONE")
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
       arguments(1) = job

       !!Attraction Integrals for all species
       do f = 1, size(MolecularSystem_instance%species)

          arguments(2) = MolecularSystem_instance%species(f)%name

          call DirectIntegralManager_getAttractionIntegrals(MolecularSystem_instance,f,integralsMatrix)         

          !!Write integrals to file (unit 30)
          call Matrix_writeToFile(integralsMatrix, unit=30, binary=.true., arguments = arguments )  
          
          if(CONTROL_instance%LAST_STEP) then
             write(*,"(A, A ,A,I6)")" Number of Nuclear integrals for species ", &
                  trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix%values,DIM=1)**2
          end if
          
       !!Depuration block
       !        print*, "Attraction  Matrix for specie: ", f
       !Call Matrix_show(integralsMatrix)

       end do !done! 
    end if

  end subroutine IntegralManager_writeAttractionIntegrals

  !> 
  !! @brief Calculate moment integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_writeMomentIntegrals()
    implicit none

    integer :: f
    integer :: component
    Type(Matrix) :: integralsMatrix
    character(10) :: coordinate(9)
    character(100) :: arguments(2)

    coordinate = ["X0", "Y0", "Z0","XX","YY","ZZ","XY","XZ","YZ"]

    !!Moment Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)
       arguments(2) = MolecularSystem_instance%species(f)%name

       do component = 1, 9 !! components x, y, z, XX, YY, ZZ, XY, XZ, YZ

          arguments(1) = "MOMENT"//trim(coordinate(component))

          call DirectIntegralManager_getMomentIntegrals(MolecularSystem_instance,f,component,integralsMatrix)

          !!Write integrals to file (unit 30)
          call Matrix_writeToFile(integralsMatrix, unit=30, binary=.true., arguments = arguments )  
          
          !!Depuration block
          ! print*, "Moment Matrix for specie: ", f, " and component ", component
          !Call Matrix_show(integralsMatrix)
       end do !! Done by component
       
       if(CONTROL_instance%LAST_STEP) then
          write(*,"(A, A ,A,I6)")" Number of Moment integrals for species ", &
               trim(MolecularSystem_instance%species(f)%name), ": ", 3*size(integralsMatrix%values,DIM=1)**2
       end if

    end do !done! 

  end subroutine IntegralManager_writeMomentIntegrals

  !> 
  !! @brief Calculate Intra-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  !!      - 2016.06.13: Use Libint V 2.x
  subroutine IntegralManager_writeIntraRepulsionIntegrals(nameOfSpecies, scheme)
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

  end subroutine IntegralManager_writeIntraRepulsionIntegrals


  !> 
  !! @brief Calculate Inter-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  !!      - 2016.06.13: Use Libint V 2.x
  subroutine IntegralManager_writeInterRepulsionIntegrals(scheme)
    implicit none
    character(*) :: scheme    
    integer :: i, j
    integer :: auxCounter = 0

    !! Skip integrals calculation two times for electrons alpha and beta    
    ! if(CONTROL_instance%IS_OPEN_SHELL .and. ( trim(nameOfSpecies) == "E-BETA" )) return

    if ( trim(String_getUppercase( CONTROL_instance%INTEGRAL_STORAGE )) == "DIRECT") return 

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies

       if(CONTROL_instance%IS_OPEN_SHELL .and. trim(MolecularSystem_getNameOfSpecie(i)) == "E-BETA" ) cycle

       do j = i+1, MolecularSystem_instance%numberOfQuantumSpecies

          if(trim(MolecularSystem_getNameOfSpecie(j)) == "E-BETA" .and. .not. trim(MolecularSystem_getNameOfSpecie(i)) == "E-ALPHA" ) cycle

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

  end subroutine IntegralManager_writeInterRepulsionIntegrals

  !> 
  !! @brief Calculate three center integrals and write it on file as a matrix
  !! @version 1.0
  subroutine IntegralManager_writeThreeCenterIntegrals()
    implicit none

    integer :: f
    Type(Matrix) :: integralsMatrix
    character(100) :: arguments(2)

    arguments(1) = "EXTERNAL_POTENTIAL"
    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       arguments(2) = MolecularSystem_instance%species(f)%name

       call DirectIntegralManager_getExternalPotentialIntegrals(MolecularSystem_instance,f,integralsMatrix)

       !! Write integrals to file (unit 30)
       call Matrix_writeToFile(integralsMatrix, unit=30, binary=.true., arguments = arguments )  

       if(CONTROL_instance%LAST_STEP) then
          write(*,"(A, A ,A,I6)")" Number of 3C Overlap integrals for species ", &
               trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix%values,DIM=1)**2
       end if

    end do !done!    

  end subroutine IntegralManager_writeThreeCenterIntegrals

  !> 
  !! @brief Calculate three center integrals and write it on file as a matrix
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine IntegralManager_writeThreeCenterIntegralsByProduct()
    implicit none

    integer :: f, g, h, i
    integer :: j, k, l, m
    integer :: p, potID
    integer :: ii, jj, hh
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:), auxintegralValue(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(100) :: job, speciesName
    type (ContractedGaussian) :: auxContractionA, auxContractionB, auxContraction

    job = "EXTERNAL_POTENTIAL"
    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       speciesName = MolecularSystem_instance%species(f)%name
       write(30) speciesName

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = DirectIntegralManager_getLabels(MolecularSystem_instance%species(f))

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
          write(*,"(A, A ,A,I6)")" Number of Overlap integrals for species ", &
               trim(MolecularSystem_instance%species(f)%name), ": ", size(integralsMatrix,DIM=1)**2
       end if
       write(30) int(size(integralsMatrix),8)
       write(30) integralsMatrix

    end do !done!    

  end subroutine IntegralManager_writeThreeCenterIntegralsByProduct


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
