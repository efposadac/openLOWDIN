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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module ChargesEQeq_
  use Vertex_
  use Vector_
  use MolecularSystem_
  use Matrix_
  use Exception_
  implicit none


  public :: &
       ChargesEQeq_getCharges

contains

  subroutine ChargesEQeq_getCharges(partialCharges, vertices)
    implicit none
    type(Vertex), intent(in) :: vertices
    real(8), allocatable, intent(in out) :: partialCharges(:)
    integer :: i, j
    real(8), allocatable :: electronegativities(:) !! Electronegativities = X = (IP + EA)/2
    real(8), allocatable :: hardness(:) !! Hardnes for atom J = IP - EA
    real(8) :: ionizationPotential
    real(8) :: electronAffinity
    real(8), allocatable :: B(:) !! electronegativities diferences
    type(Matrix) :: A
    integer :: ionizationPosition
    integer :: chargeCenter
    integer :: totalCharge
    integer(8) :: numberOfCenters
    real(8) :: Ja, Jb !! dummies for idempotentials
    integer :: info

    numberOfCenters = vertices%numberOfVertices
    totalCharge = MolecularSystem_instance%charge
    ! write(*,"(T20,A,I)") "Carga Total: ", totalCharge
    
    allocate( partialCharges( vertices%numberOfVertices ) )
    allocate( electronegativities( vertices%numberOfVertices ) )
    allocate( hardness( vertices%numberOfVertices ) )

    do i=1,vertices%numberOfVertices
       partialCharges(i) = 0.0
    end do

    chargeCenter = 0 !! In the future we can change this for any charges

    do i=1,vertices%numberOfVertices
       if(trim(vertices%type(i)) == "H_" .or. trim(vertices%type(i)) == "H_b") then !! Correction for Hydrogen
          ionizationPotential = 13.598
          electronAffinity = -2.0
          electronegativities(i) = 0.5*(ionizationPotential + electronAffinity)
          hardness(i) = ionizationPotential - electronAffinity
       else
          ionizationPosition = chargeCenter + 1
          ionizationPotential = vertices%ionizationPotential(i)%values(1,ionizationPosition+1)
          electronAffinity = vertices%ionizationPotential(i)%values(1,ionizationPosition)
          electronegativities(i) = 0.5*(ionizationPotential + electronAffinity)
          hardness(i) = ionizationPotential - electronAffinity
          electronegativities(i) = electronegativities(i) - chargeCenter*hardness(i) 
       end if
    end do

    !! We need to solve the system Ax = b, but, first we need to obtain A and b matrix
    allocate( B( vertices%numberOfVertices ) )
    B(1) = totalCharge

    do i=2,vertices%numberOfVertices
       B(i) = electronegativities(i) - electronegativities(i-1)
    end do
     
    call Matrix_constructor( A, numberOfCenters, numberOfCenters)
    do i=1,vertices%numberOfVertices
       A%values(1,i) = 1.0
    end do

    do i=2,vertices%numberOfVertices
       do j=1,vertices%numberOfVertices
          call ChargesEQeq_getIdempotential(i-1,j,hardness,vertices,Ja)
          call ChargesEQeq_getIdempotential(i,j,hardness,vertices,Jb)          
          A%values(i,j) = Ja - Jb
       end do
    end do

    ! write(*,"(T20,A)") "B original"
    ! do i=1,vertices%numberOfVertices
    !    write(*,"(T20,F12.5)") B(i)
    ! end do

    ! write(*,"(T20,A)") "A values"
    ! do i=1,vertices%numberOfVertices
    !    write(*,"(T20,<vertices%numberOfVertices>F12.5)") A%values(i,:)
    ! end do

    ! partialCharges = Matrix_solveLinearEquation(A, B)

    call Matrix_linear(vertices%numberOfVertices, 1, A, vertices%numberOfVertices, B, vertices%numberOfVertices, partialCharges, info)

    
    write(*,"(T20,A,I)") "Info: ", info
    ! do i=1,vertices%numberOfVertices
    !    write(*,"(T20,F12.5)") partialCharges(i)
    ! end do

  end subroutine ChargesEQeq_getCharges

  subroutine ChargesEQeq_getIdempotential(i,j,hardness,vertices,idempotential) 
    implicit none
    integer, intent(in) :: i, j
    real(8), allocatable, intent(in) :: hardness(:)
    type(Vertex), intent(in) :: vertices
    real(8), intent(out) :: idempotential
    real(8) :: separationOfCenters, Rab !! distance between atoms
    real(8) :: Jij !! partial idempotential (dummy)
    real(8) :: a
    real(8) :: orbitalOverlap
    real(8) :: K !!vacuum permittivity
    real(8) :: lambda !!Coulomb scaling parameter

    K = 14.4
    lambda = 1.2

    if(i==j) then
       idempotential = hardness(i)
    else
       separationOfCenters = sum( ( vertices%cartesianMatrix%values(i,:) - vertices%cartesianMatrix%values(j,:))**2.0 )
       Rab = (sqrt(separationOfCenters)) * AMSTRONG
       Jij = dsqrt(hardness(i)*hardness(j))
       a = Jij/K
       orbitalOverlap = (exp(-(a*a*Rab*Rab)))*(2*a - a*a*Rab - 1/Rab)
       idempotential = lambda*(K/2)*((1/Rab) + orbitalOverlap)
    end if

  end subroutine ChargesEQeq_getIdempotential

end module ChargesEQeq_
