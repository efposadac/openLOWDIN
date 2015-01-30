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
  use String_
  use MolecularSystem_
  use ContractedGaussian_
  use OverlapIntegrals_
  use AttractionIntegrals_
  use MomentIntegrals_
  use KineticIntegrals_
  use LibintInterface_
  ! use CudintInterface_
  use RysQuadrature_
  use Matrix_
  use CosmoCore_
  use Stopwatch_


  implicit none

  public :: &
       IntegralManager_getOverlapIntegrals, &
       IntegralManager_getKineticIntegrals, &
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
    integer :: numCartesianOrbitalI
    integer :: numCartesianOrbitalJ
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralBuffer(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(20) :: colNum
    character(100) :: job

    job = "OVERLAP"

    !!Overlap Integrals for all species    
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralBuffer)) deallocate(integralBuffer)
       allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )
       integralBuffer = 0.0_8

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

                      end do
                   end do

                end do
                hh = 1
             end do

          end do
       end do

       !! Write integrals to file (unit 30)
       write(*,"(A,I6,A,A,A)")" Stored ", size(integralsMatrix,DIM=1)**2," overlap integrals of specie ",trim(MolecularSystem_instance%species(f)%name),&
            " in file lowdin.opints"
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
    real(8), allocatable :: integralBuffer(:)
    real(8), allocatable :: integralsMatrix(:,:)
    character(20) :: colNum
    character(100) :: job

    job = "KINETIC"

    !!Kinetic Integrals for all species
    do f = 1, size(MolecularSystem_instance%species)

       write(30) job
       write(30) MolecularSystem_instance%species(f)%name

       if(allocated(labels)) deallocate(labels)
       allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
       labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

       if(allocated(integralBuffer)) deallocate(integralBuffer)
       allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )
       integralBuffer = 0.0_8

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
       write(*,"(A,I6,A,A,A)")" Stored ", size(integralsMatrix,DIM=1)**2," kinetic integrals of specie ",trim(MolecularSystem_instance%species(f)%name),&
            " in file lowdin.opints"
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


  !> 
  !! @brief Calculate point charge - quantum particle attraction integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: reads point charge information from lowdin.bas file
  !!			- 2014.16.09: modify subroutines to calcule cosmo monoelectronic integrals Danilo
  subroutine	IntegralManager_getAttractionIntegrals(surface)
    implicit none

    integer :: f, g, h, i, c
    integer :: j, k, l, m
    integer :: ii, jj, hh
    integer :: numberOfPointCharges
    integer :: numberOfClasicalCharges
    integer :: numberOfSurfaceSegments
    integer, allocatable :: labels(:)
    real(8), allocatable :: integralValue(:)
    real(8), allocatable :: integralBuffer(:)
    real(8), allocatable :: integralsMatrix(:,:)
    type(pointCharge), allocatable :: point(:)
    character(20) :: colNum
    character(100) :: job

    ! Variables para calcular las integrales monoelectronicas para cosmo
    type(surfaceSegment),intent(in), optional :: surface
    integer :: times
    integer :: n,a
    logical :: isCosmo
    real(8), allocatable :: integralValueCosmo(:,:)
    real(8), allocatable ::  cosmoV(:)
    real(8), allocatable ::  qCharges(:)
    real(8), allocatable ::  qChargesAc(:,:)
    character(100) :: cosmoIntegralFile
    character(100) :: cosmoQuantumChargeFile
    character(100) :: cosmoClasicalChargeFile


    type(Matrix) :: cmatin
    real(8) ::  cosmoPotInt

    integer :: b,d,totals,total_aux

    a=0


    !

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
       ! write(*,*) "remplazadas por estas"

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

          ! write(*,*)"species", f
          ! write(70) job
          ! write(70) MolecularSystem_instance%species(f)%name

          if(allocated(labels)) deallocate(labels)
          allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
          labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

          if(allocated(integralBuffer)) deallocate(integralBuffer)
          allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )


          integralBuffer = 0.0_8


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

                      b=MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital
                      d=MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital
                      totals=b*d

                      total_aux=total_aux+totals


                      do c = 1, numberOfPointCharges
                         ! do sobre las cargas puntuales 
                         point(1)%charge = 1.0
                         point(1)%x  =surface%xs(c)
                         point(1)%y  =surface%ys(c)
                         point(1)%z  =surface%zs(c)
												 if (c==1) then
													write(*,*) point(c)%x,point(c)%y,point(c)%z
												end if
                         !Calculating integrals for shell

                         call AttractionIntegrals_computeShell( MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h), &
                              MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j), point, numberOfPointCharges, integralValue)
                         m=0

                         do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                            do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                               m = m + 1
                               integralValueCosmo(m,c)=integralValue(m)
                               !! debug
															 
															 write(*,*)"cosmo integrals: m,k,l",m,k,l
                               write(*,*)integralValueCosmo(m,c)
                               
															 write(37,'(F10.5)')integralValueCosmo(m,c)
                            end do
                         end do

                         ! todas las integrales para un c específico 

                         ! aqui es donde se hace el calculo para cada carga

                      end do
                      !! end do de cargas puntuales

                      if(allocated(qCharges)) deallocate(qCharges)
                      allocate(qCharges(numberOfPointCharges))


                      m=0
                      do k = labels(ii), labels(ii) + (MolecularSystem_instance%species(f)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                         do l = labels(jj), labels(jj) + (MolecularSystem_instance%species(f)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                            m = m + 1
                            if(allocated(cosmoV)) deallocate(cosmoV)
                            allocate(cosmoV(numberOfPointCharges))
                            cosmoV(:)=integralValueCosmo(m,:)
                            call CosmoCore_q_builder(cmatin, cosmoV, numberOfPointCharges, qCharges)

														! write(*,*)"cosmo integrals: m,k,l",m,k,l
                            ! write(*,*)integralValueCosmo(m,:)
														! write(*,*)"cosmo charges"
                            ! write(*,*)qCharges
														
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
          totals=total_aux
          call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoQuantumChargeFile,numberOfPointCharges,totals,totals)

          !!clasical vs clasical

          cosmoClasicalChargeFile="cosmo.clasical"

					job="COSMO1"

					write(40) job
          write(40) MolecularSystem_instance%species(f)%name

          call CosmoCore_q_int_builder(cosmoIntegralFile,cosmoClasicalChargeFile,numberOfPointCharges,1,totals,labels,f)

          !clasical vs quantum
          
					
					job="COSMO4"
					write(40) job
          write(40) MolecularSystem_instance%species(f)%name

          call CosmoCore_nucleiPotentialQuantumCharges(surface,cosmoQuantumChargeFile,totals,labels,f)

       end do
       !! end do especies

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

          if(allocated(labels)) deallocate(labels)
          allocate(labels(MolecularSystem_instance%species(f)%basisSetSize))
          labels = IntegralManager_getLabels(MolecularSystem_instance%species(f))

          if(allocated(integralBuffer)) deallocate(integralBuffer)
          allocate(integralBuffer((MolecularSystem_instance%species(f)%basisSetSize * (MolecularSystem_instance%species(f)%basisSetSize + 1)) / 2 ) )
          integralBuffer = 0.0_8

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

                            ! write(*,*)"lowdin integrals: m,k,l",m,k,l

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
          write(*,"(A,I6,A,A,A)")" Stored ",size(integralsMatrix,DIM=1)**2," attraction integrals of specie ",trim(MolecularSystem_instance%species(f)%name),&
               " in file lowdin.opints"
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
    character(20) :: colNum
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

          write(30) job
          write(30) MolecularSystem_instance%species(f)%name
          write(30) coordinate(component)

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

       write(*,"(A,I6,A,A,A)")" Stored ",(size(integralsMatrix,DIM=1)**2)*3,&
            " Moment integrals of specie ",trim(MolecularSystem_instance%species(f)%name),&
            " in file lowdin.opints"

    end do !done! 

  end subroutine IntegralManager_getMomentIntegrals

  !> 
  !! @brief Calculate Intra-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  subroutine IntegralManager_getIntraRepulsionIntegrals(nprocess, process, nameOfSpecies, scheme)
    implicit none

    integer(8) :: nprocess
    integer(8) :: process
    character(*) :: nameOfSpecies
    character(*) :: scheme

    integer :: speciesID
    integer :: numberOfContractions
    integer(8) :: integralsByProcess
    integer(8) :: ssize
    integer(8) :: starting
    integer(8) :: ending

    !! Skip integrals calculation two times for electrons alpha and beta    
    if(CONTROL_instance%IS_OPEN_SHELL .and. ( trim(nameOfSpecies) == "E-BETA" )) return

    speciesID = MolecularSystem_getSpecieID(trim(nameOfSpecies))
    numberOfContractions = MolecularSystem_getNumberOfContractions(speciesID)


    ssize = (numberOfContractions * (numberOfContractions + 1))/2
    ssize = (ssize * (ssize + 1))/2


    integralsByProcess = ceiling( real(ssize,8)/real(nprocess,8) )

    ending = process * integralsByProcess
    starting = ending - integralsByProcess + 1

    if( starting > ssize ) return

    if( ending > ssize ) ending = ssize

    !! Calculate integrals (stored on disk)           
    select case (trim(String_getUppercase(trim(scheme))))
    
		case("RYS")
       call RysQuadrature_computeIntraSpecies( speciesID, "ERIS", starting, ending, int(process) )
    case("LIBINT")
       call LibintInterface_computeIntraSpecies( speciesID, "ERIS", starting, ending, int(process) )
    ! case("CUDINT")
    !    call CudintInterface_computeIntraSpecies(speciesID)
    case default
       call LibintInterface_computeIntraSpecies( speciesID, "ERIS", starting, ending, int(process) )
    end select

  end subroutine IntegralManager_getIntraRepulsionIntegrals

  !> 
  !! @brief Calculate Inter-species repulsion integrals
  !! @author E. F. Posada, 2013
  !! @version 1.0
  !! @par History
  !!      - 2013.03.05: Use Libint V 1.1.4
  subroutine IntegralManager_getInterRepulsionIntegrals(scheme)
    implicit none
    character(*) :: scheme    
    integer :: i, j

    do i = 1, MolecularSystem_instance%numberOfQuantumSpecies
       do j = i+1, MolecularSystem_instance%numberOfQuantumSpecies

          !! Calculate integrals (stored on disk)       
          select case (trim(String_getUppercase(trim(scheme))))
          case("LIBINT")
             call LibintInterface_computeInterSpecies( i, j, "ERIS" )
          ! case("CUDINT")
          !    call CudintInterface_computeInterSpecies( i, j, "ERIS" )
          case default
             call LibintInterface_computeInterSpecies( i, j, "ERIS" )
          end select

       end do       
    end do

  end subroutine IntegralManager_getInterRepulsionIntegrals

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
		
		write(*,*)"labels data from integral manager"

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

end module IntegralManager_
