!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	PROF. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief This module creates grids for integration of the exchange correlation energy, potential or kernel. Partially based on R. Flores-Moreno Parakata's modules
!! @author F. Moncada, 2017
module Grid_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use String_
  use Exception_
  use Stopwatch_
  use Lebedev_

  implicit none

  type, public :: Grid

     type(MolecularSystem), pointer :: molSys
     
     character(30) :: nameOfSpecies
     integer :: totalSize
     type(Matrix) :: points !! x,y,z,weight
     type(Matrix), allocatable :: orbitalsWithGradient(:) !!x,y,z per orbital
     type(Vector) :: potential
     type(Vector) :: gradientPotential(3)
     type(Vector) :: density
     type(Vector) :: densityGradient(3)

  end type Grid

  public :: &
       Grid_constructor, &
       Grid_copyPoints,&
       Grid_buildAtomic, &
       Grid_radialCutoff, &
       Grid_weightCutoff, &
       Grid_radialQuadrature, &
       Grid_EulerMaclaurinQuadrature

contains

  !>
  !! @brief Builds a grid for each species - Different sizes are possible, all points in memory
  ! Felix Moncada, 2017
  ! Roberto Flores-Moreno, 2009
  subroutine Grid_constructor( this, speciesID, type, molSys )
    implicit none
    type(Grid) :: this
    integer :: speciesID
    character(*) :: type
    type(MolecularSystem), target :: molSys
    
    type(Matrix) :: atomicGrid, molecularGrid
    integer :: numberOfSpecies, numberOfCenters
    integer :: radialSize, angularSize, numberOfShells, initialGridSize, molecularGridSize
    integer :: particleID, particleID2, particleID3, point, i
    real(8) :: cutoff, sum, r, w, mu
    real(8), allocatable :: origins(:,:), distance(:),factor(:)
    integer, allocatable :: atomicGridSize(:)

    this%molSys=>molSys

    this%nameOfSpecies=trim(MolecularSystem_getNameOfSpecies(speciesID,this%molSys))
    
    if (trim(type) .eq. "INITIAL") then
       radialSize=CONTROL_instance%GRID_RADIAL_POINTS
       angularSize=CONTROL_instance%GRID_ANGULAR_POINTS
       numberOfShells=CONTROL_instance%GRID_NUMBER_OF_SHELLS
    elseif (trim(type) .eq. "FINAL") then
       radialSize=CONTROL_instance%FINAL_GRID_RADIAL_POINTS
       angularSize=CONTROL_instance%FINAL_GRID_ANGULAR_POINTS
       numberOfShells=CONTROL_instance%FINAL_GRID_NUMBER_OF_SHELLS
    else
       STOP "What grid were you trying to build?"
    end if
    
    if(CONTROL_instance%PRINT_LEVEL .gt. 0) &
         write (*,"(A,I4,A,I2,A,A)") " Building an atomic grid with", radialSize, " radial points in ", numberOfShells, " shells for ", trim(this%nameOfSpecies)
    
    numberOfCenters=size(this%molSys%species(speciesID)%particles)
    allocate(origins(numberOfCenters,3), atomicGridSize(numberOfCenters))

    !Get Atomic Grid
    call Grid_buildAtomic(radialSize , angularSize, numberOfShells, atomicGrid, initialGridSize )
    ! call Matrix_show(atomicGrid)

    !Get origins and Build the complete Grid for the species
    !We are screening the points with delocalized orbital values lower than 1E-6
    molecularGridSize=0
    do particleID = 1, numberOfCenters
       call Grid_radialCutoff( atomicGrid, initialGridSize, speciesID, particleID, atomicGridSize(particleID), this%molSys)
       molecularGridSize=molecularGridSize + atomicGridSize(particleID)
    end do
    
    call Matrix_constructor( molecularGrid, int(molecularGridSize,8), int(4,8), 0.0_8 )            

    i=1
    do particleID = 1, numberOfCenters
       origins(particleID,1:3) = this%molSys%species(speciesID)%particles(particleID)%origin(1:3)
       do point = 1, atomicGridSize(particleID)
          molecularGrid%values(i,1)=atomicGrid%values(point,1)+origins(particleID,1)
          molecularGrid%values(i,2)=atomicGrid%values(point,2)+origins(particleID,2)
          molecularGrid%values(i,3)=atomicGrid%values(point,3)+origins(particleID,3)
          molecularGrid%values(i,4)=atomicGrid%values(point,4)
          i=i+1
       end do
    end do
    
    ! call Matrix_show(molecularGrid)

    !Calculate adecuate weights with Becke's
    allocate(distance(numberOfCenters),factor(numberOfCenters))

    i=1
    do particleID3 = 1, numberOfCenters
       do point = 1, atomicGridSize(particleID3)
          ! call Grid_Becke(numberOfCenters)

          do particleID=1, numberOfCenters
             distance(particleID)= sqrt( (origins(particleID,1)-molecularGrid%values(i,1))**2 +&
                  (origins(particleID,2)-molecularGrid%values(i,2))**2 +&
                  (origins(particleID,3)-molecularGrid%values(i,3))**2)
             factor(particleID)=1.0
          end do

          do particleID=2, numberOfCenters
             do particleID2=1, particleID-1

                r=sqrt( (origins(particleID,1)-origins(particleID2,1) )**2 +&
                     (origins(particleID,2)-origins(particleID2,2) )**2 +&
                     (origins(particleID,3)-origins(particleID2,3) )**2 )
                w=(distance(particleID)-distance(particleID2))/r
                mu=Grid_Beckefun(Grid_Beckefun(Grid_Beckefun(w)))

                factor(particleID)=(1.0-mu)*factor(particleID)
                factor(particleID2)=(1.0+mu)*factor(particleID2)

             end do
          end do

          sum = 0.0
          do particleID=1, numberOfCenters
             sum=sum+factor(particleID)
          end do

          molecularGrid%values(i,4)=molecularGrid%values(i,4)*factor(particleID3)/sum
          i=i+1

       end do
    end do

    !We are screening points with weight lower than 1E-16
    this%totalSize=0
    call Grid_weightCutoff( molecularGrid, molecularGridSize, this%points, this%totalSize )

    ! print *, "Grid for species", speciesID 
    ! call Matrix_show(this%points)


    if(CONTROL_instance%PRINT_LEVEL .gt. 0) then
       write(*,"(A,ES9.3,A,ES9.3,A)") "Screening delocalized orbital(<", CONTROL_instance%ELECTRON_DENSITY_THRESHOLD,&
            ") and low weight(<",CONTROL_instance%GRID_WEIGHT_THRESHOLD,") points ..."
       print *, "Final molecular grid size for: ", trim(this%nameOfSpecies), this%totalSize ," points"
       print *, " "
    end if
    
    call Matrix_destructor( atomicGrid)
    call Matrix_destructor( molecularGrid)
    deallocate(origins, distance,factor, atomicGridSize)



  end subroutine Grid_constructor

  !>
  !! @brief Copies a grid object otherThis->this
  ! Felix Moncada, 2025
  subroutine Grid_copyPoints( this, otherThis )
    implicit none
    type(Grid) :: this
    type(Grid) :: otherThis

    this%totalSize=otherThis%totalSize
    call Matrix_copyConstructor( this%points, otherThis%points )

  end subroutine Grid_copyPoints
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine Grid_exception( typeMessage, description, debugDescription)
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

  end subroutine Grid_exception
  
  subroutine Grid_radialCutoff(atomicGrid, gridSize, speciesID, particleID, relevantPoints, molSys)
    ! Gets the radial point where basis sets take negligible values
    ! Felix Moncada, 2017
    implicit none
    type(matrix) :: atomicGrid
    integer :: gridSize, speciesID, particleID   
    integer :: relevantPoints !output
    type(molecularSystem) :: molSys

    integer :: numberOfContractions
    integer :: numberOfPrimitives

    integer :: point, mu, i
    real(8) :: minExp, normC, orbitalCutoff, radialCutoff

    ! Search for the lowest exponent in the atomic basis set
    minExp=1E12
    numberOfContractions = size(molSys%species(speciesID)%particles(particleID)%basis%contraction)
    do mu = 1, numberOfContractions
       numberOfPrimitives = size(molSys%species(speciesID)%particles(particleID)%basis%contraction(mu)%orbitalExponents)
       do i = 1, numberOfPrimitives
          if (molSys%species(speciesID)%particles(particleID)%basis%contraction(mu)%orbitalExponents(i) .lt. minExp) then
             minExp=molSys%species(speciesID)%particles(particleID)%basis%contraction(mu)%orbitalExponents(i)
          end if
       end do
    end do
    
    ! Search for the points where the GTF with the lowest exponent in the atomic basis set has a value larger than the cutoff
    ! We assume that it is an s-function for simplicity

    normC=(2*minExp/Math_PI)**(3/4)
    orbitalCutoff=CONTROL_instance%ELECTRON_DENSITY_THRESHOLD
    radialCutoff=sqrt(1/minExp*log(normC/orbitalCutoff))

    relevantPoints=0
    do point = 1, gridSize
       if ( sqrt(atomicGrid%values(point,1)**2 + atomicGrid%values(point,2)**2 + atomicGrid%values(point,3)**2 ) .lt. radialCutoff ) then
          relevantPoints=relevantPoints+1
       end if
    end do
    
    ! print *, "radialCutoff for", speciesID, particleID, radialCutoff, relevantPoints
  
  end subroutine Grid_radialCutoff

  subroutine Grid_weightCutoff(molecularGrid, molecularGridSize, finalGrid, finalGridSize)
    ! Removes the points where weights take negligible values
    ! Felix Moncada, 2017
    implicit none
    type(matrix) :: molecularGrid, finalGrid !molecular=input, final=output
    integer :: molecularGridSize, finalGridSize
    integer :: speciesID, particleID
    integer :: weightCutOff

    integer :: point, rpoint

    !Should be a control parameter
    weightCutoff=CONTROL_instance%GRID_WEIGHT_THRESHOLD
    
    finalGridSize=0
    do point = 1, molecularGridSize
       if ( molecularGrid%values(point,4) .gt. weightCutoff ) then
          finalGridSize=finalGridSize+1
       end if
    end do

    call Matrix_constructor( finalGrid, int(finalGridSize, 8) , int(4,8) , 0.0_8 )

    rpoint=0
    do point = 1, molecularGridSize
       if ( molecularGrid%values(point,4) .gt. weightCutoff ) then
          rpoint=rpoint+1
          finalGrid%values(rpoint,1:4)=molecularGrid%values(point,1:4)
       end if
    end do

    ! print *, "radialCutoff for", speciesID, particleID, radialCutoff, relevantPoints
    
  end subroutine Grid_weightCutoff


  subroutine Grid_buildAtomic(nrad,nang, numberOfShells, agrid,initialSize)
    ! Generate an atomic grid
    ! Felix Moncada, 2017
    ! Roberto Flores-Moreno, May 2009
    implicit none
    integer :: nrad, nang !! grid sizes, Inputs
    integer :: numberOfShells !! Input, to build grids with a different sizes for regions near the nuclei
    type(Matrix) :: agrid !! points and weights, Output
    integer :: initialSize !! total size of the built grid, Output
    ! real(8) :: cutoff
    real(8), allocatable :: r(:),wr(:)
    real(8), allocatable :: x(:,:),y(:,:),z(:,:),wa(:,:)

    integer :: iang,irad, n, i, j
    real(8) :: solida, switch
    integer :: lebedevNumbers(10) !(32)
    integer :: maxLebedev
    integer, allocatable :: shellAng(:), shellRadialChange(:)

    
    ! Get angular points, for numberOfShells different distributions
    ! lebedevNumbers(1:32)=(/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810/)

    lebedevNumbers(1:10)=(/6,14,26,50,110,194,302,434,590,770/)

    maxLebedev=0
    do i=1, size(lebedevNumbers)
       if( lebedevNumbers(i) .eq. nang) maxLebedev=i
    end do

    if( maxLebedev .eq. 0) STOP "The number of angular points chosen is not supported!, choose from 6,14,26,50,110,194,302,434,590,770"

    
    allocate(x(nang,numberOfShells),y(nang,numberOfShells),z(nang,numberOfShells),wa(nang,numberOfShells))
    allocate(shellAng(numberOfShells), shellRadialChange(numberOfShells+1))

    ! Get radial points
    allocate(r(nrad),wr(nrad))
    call Grid_radialQuadrature(r,wr,nrad)

    ! Depending on the number of shells we build our angular grid
    select case(numberOfShells)         
    case(1)
       !The whole space with the same angular grid
       shellAng(1)=lebedevNumbers(maxLebedev) 
       call Lebedev_angularGrid(x(1:shellAng(1),1),y(1:shellAng(1),1),z(1:shellAng(1),1),wa(1:shellAng(1),1),shellAng(1))

       shellRadialChange(1)=0
       shellRadialChange(2)=nrad

       
    case(2)
       !Near the nuclei few points that far from it
       
       shellAng(1)=lebedevNumbers(maxLebedev-1)
       call Lebedev_angularGrid(x(1:shellAng(1),1),y(1:shellAng(1),1),z(1:shellAng(1),1),wa(1:shellAng(1),1),shellAng(1))
       shellAng(2)=lebedevNumbers(maxLebedev) 
       call Lebedev_angularGrid(x(1:shellAng(2),2),y(1:shellAng(2),2),z(1:shellAng(2),2),wa(1:shellAng(2),2),shellAng(2))

       !We divide at 1 a.u.
       shellRadialChange(1)=0
       shellRadialChange(numberOfShells+1)=nrad

       switch=1.0_8
       do irad=nrad, 1, -1
          if( log10(r(irad)) .lt. log10(switch) ) then
             shellRadialChange(2)=irad
             exit
          end if
       end do
       
    case default
       !The most important (valence) part of the density lies between 1 a.u. and 10 a.u. 
       !The largest angular grid goes in that region

       do i=1, numberOfShells-1
          shellAng(i)=lebedevNumbers(maxLebedev-numberOfShells+1+i) !! Inner shells need less points than outer shells
          call Lebedev_angularGrid(x(1:shellAng(i),i),y(1:shellAng(i),i),z(1:shellAng(i),i),wa(1:shellAng(i),i),shellAng(i))
       end do

       shellAng(numberOfShells)=lebedevNumbers(maxLebedev-1) !!the outermost shell is not that important
       call Lebedev_angularGrid(x(1:shellAng(numberOfShells),numberOfShells),y(1:shellAng(numberOfShells),numberOfShells),&
            z(1:shellAng(numberOfShells),numberOfShells), wa(1:shellAng(numberOfShells),numberOfShells), shellAng(numberOfShells))

       !The last change occurs at 10 a.u., below that, logaritmically
       shellRadialChange(1)=0
       shellRadialChange(numberOfShells+1)=nrad

       switch=10.0_8
       i=numberOfShells
       do irad=nrad, 1, -1
          if( log10(r(irad)) .lt. log10(switch) ) then
             shellRadialChange(i)=irad
             switch=switch/10
             if(i .eq. 2 ) then
                exit
             else
                i = i-1
             end if
          end if
       end do
       
    end select
    

    !We are distributing evenly the points - this has not been tested

    ! do i=2, numberOfShells
       
    !       if( irad .lt. shellRadialChange(i-1)+nrad/numberOfShells) shellRadialChange(i)=irad !I need a function to define this, the first one should be zero
    ! end do

    ! print *, "radial points", "angular grid"
    
    ! do i=1, numberOfShells
    !    do j= shellRadialChange(i)+1 , shellRadialChange(i+1)
    !       print *, j, r(j), shellAng(i)
    !    end do
    ! end do
    
    initialSize=0
    do i=1, numberOfShells
       if(CONTROL_instance%PRINT_LEVEL .gt. 0) &
            write (*,"(A,F10.4,A,F10.4,A,I4,A)")  " Between radii ", r(shellRadialChange(i)+1), " a.u. and ", r(shellRadialChange(i+1)), " a.u. with", shellAng(i), " angular points"
       do j= shellRadialChange(i)+1 , shellRadialChange(i+1)

          do iang=1,shellAng(i)
             initialSize=initialSize+1
          end do
       end do
    end do

    ! print *, "Lebedev Grids Used", shellAng
    if(CONTROL_instance%PRINT_LEVEL .gt. 0) &
         print *, "Number of points in the atomic grid", initialSize
    
    call Matrix_constructor( agrid, int(initialSize,8), int(4,8), 0.0_8 )
    ! Generate angular distributions
    n = 0
    irad=0
    do i=1, numberOfShells
       do j= shellRadialChange(i)+1 , shellRadialChange(i+1)
          irad=irad + 1
          solida = 4.0*Math_PI*r(irad)**2 !! 4*pi*r^2
          do iang=1,shellAng(i)
             n = n + 1
             agrid%values(n,1) = x(iang,i)*r(irad) 
             agrid%values(n,2) = y(iang,i)*r(irad) 
             agrid%values(n,3) = z(iang,i)*r(irad) 
             agrid%values(n,4) = wa(iang,i)*wr(irad)*solida
             ! print *, i, j, n, r(irad), x(iang,i), y(iang,i), z(iang,i)
             ! print *, i, j, n, wr(irad), wa(iang,i), agrid%values(n,4)
          end do
       end do
    end do

    ! call Matrix_show(agrid)
    
    deallocate(r,wr,x,y,z,wa,shellAng)

  end subroutine Grid_buildAtomic

  subroutine Grid_radialQuadrature(r,w,n)
    ! Radial quadrature
    ! Roberto Flores-Moreno, May 2009
    implicit none
    integer n
    real(8) r(*),w(*)

    integer allocs
    real(8),allocatable :: rp(:),wp(:)

    allocate(rp(n),wp(n),stat=allocs)
    if (allocs.gt.0) stop 'allocation failed'

    ! Transformed Gauss-Chebyshev quadratures
    !rfm call tskgc(r,w,n)
    !w(:n) = (w(:n)/(1.0-r(:n)))/log(2.0)
    !r(:n) = 1.0 - log(1.0-r(:n))/log(2.0)

    ! Euler-Maclaurin
    call Grid_EulerMaclaurinQuadrature(rp,wp,n)
    r(:n) = 0.31*(rp(:n)/(1.0 - rp(:n)))**2
    w(:n) = 2.0*wp(:n)*r(:n)/(rp(:n)-rp(:n)**2)

    deallocate(rp,wp,stat=allocs)
    if (allocs.gt.0) stop 'deallocation failed'
  end subroutine Grid_radialQuadrature

  subroutine Grid_EulerMaclaurinQuadrature(x,w,n)
    ! Abscissas and weights for the Euler-Maclaurin quadrature.
    ! Robert Flores-Moreno, Jun 2010
    implicit none
    integer n
    real(8) x(*),w(*)

    integer i

    do i=1,n
       x(i) = float(i)/float(n+1)
       w(i) = 1.0/float(n+1)
    end do
  end subroutine Grid_EulerMaclaurinQuadrature
  
  real(8) function Grid_Beckefun(mu)
    ! Becke iteration function
    ! Roberto Flores-Moreno, Jun 2010
    implicit none
    real(8) mu
    Grid_Beckefun = 0.5*mu*(3.0 - mu**2)
  end function Grid_Beckefun

end module Grid_


  
  ! subroutine set_grid(m,g,ng,atom)
  !   ! Build molecular grid
  !   ! Roberto Flores-Moreno, May 2009, Jun 2010
  !   implicit none
  !   integer ng,atom
  !   real(8) g(4,*)
  !   type(pmolecule) m

  !   integer ig

  !   ! Get atomic grid
  !   do ig=1,ng
  !      g(1:3,ig) = rg(1:3,ig) + m%atom(atom)%pos(1:3)
  !      g(4,ig) = rg(4,ig)
  !   end do

  !   ! Becke weights
  !   call becke(m,atom,g,ng)
  ! end subroutine set_grid
  
  ! subroutine Grid_Becke(m,atom,g,n)
  !   ! Build and apply Becke weights
  !   ! Roberto Flores-Moreno, Jun 2010
  !   !
  !   ! Lit.: A.D. Becke, J. Chem. Phys. 88, 2547 (1988)
  !   !
  !   implicit none
  !   integer atom,n
  !   real(8) g(4,*)
  !   type(pmolecule) m

  !   integer i,iatom,jatom
  !   real(8) s,mu
  !   real(8) r(pmaxatom),p(pmaxatom)

  !   do i=1,n
  !      s = 0.0
  !      do iatom=1,m%natom
  !         r(iatom) = sqrt((g(1,i)-m%atom(iatom)%pos(1))**2 +              &
  !              (g(2,i)-m%atom(iatom)%pos(2))**2 +              &
  !              (g(3,i)-m%atom(iatom)%pos(3))**2)
  !         p(iatom) = 1.0
  !      end do
  !      do iatom=2,m%natom
  !         do jatom=1,iatom-1
  !            mu = sqrt((m%atom(jatom)%pos(1)-m%atom(iatom)%pos(1))**2 +    &
  !                 (m%atom(jatom)%pos(2)-m%atom(iatom)%pos(2))**2 +    &
  !                 (m%atom(jatom)%pos(3)-m%atom(iatom)%pos(3))**2)
  !            mu = (r(iatom)-r(jatom))/mu
  !            mu = beckefun(beckefun(beckefun(mu)))
  !            p(iatom) = (1.0 - mu)*p(iatom)
  !            p(jatom) = (1.0 + mu)*p(jatom)
  !         end do
  !      end do
  !      do iatom=1,m%natom
  !         s = s + p(iatom)
  !      end do
  !      g(4,i) = g(4,i)*p(atom)/s
  !   end do
  ! end subroutine Grid_Becke
  
!     subroutine dft_initialize(m)
!   ! Initialize dft calculation
!   ! Roberto Flores-Moreno, Oct 2010
!   implicit none
!     integer :: iatom, ng
!     type(pmolecule) :: m
!     ! real(8) :: startTime, endTime
!     ! call cpu_time(startTime)
   
 	
!     call bld_grid(rg,ng,nradp,nangp)

! 	    do iatom=1,m%natom
! 	       call set_grid(m,g,ng,iatom)
! !	       call writeGrid(g,ngp,iatom,"write",m%specie)
! 	    end do

!     call screenBasis(m)

!     ! call cpu_time(endTime)
!     ! write(6,"(A,F15.3,A)") "*****Time for grid construction + screen Basis", (endTime - startTime), "(s)"

! end subroutine dft_initialize


  
! !!!!!!!!!!!!!!!!!!!!!
! !!! I/O 
!   subroutine writeGrid(grid,ng,iatom,option,nameOfSpecie)
!     ! Write/read atomic grids for DFT
!     ! Roberto Flores-Moreno, Oct 2008
!     implicit none

!     character*(*) option
!     integer :: iatom
!     integer :: ng
!     character(20) :: nameOfSpecie
!     real(8) :: grid(*)

!     integer i,j,tape

!     ! Determine record number

!     ! Open file
!     call open_binary_file(trim(nameOfSpecie)//'.grid',tape,4*ng)

!     ! Write data
!     if (option.eq.'write') then
!        write(tape,rec=iatom) (grid(i), i=1, 4*ng)

!        ! Read data
!     else if (option.eq.'read') then
!        read(tape,rec=iatom)  (grid(i), i=1, 4*ng)
!     end if

!     ! Close file
!     call close_file(tape)
!   end subroutine writeGrid

