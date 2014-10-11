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
!! @brief Molecular Mechanics program.
!!        This module calculates the partial charges  system using the EQeq(Extended Charge Equilibration) approach,
!! this calculation is optional and must be activated in the input like this:
!! <BLOCKQUOTE>
!! CONTROL \n
!! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;electrostaticMM = T \n
!! END CONTROL \n
!! </BLOCKQUOTE>
!! @note For reference see: \n
!! <b>Original Charge Equilibration approach (QEq):</b> \n
!! \n
!! Rappe, A.K., Goddard III, W.A., <b>Charge Equilibration for Molecular Dynamics Simulations</b>,
!! J. Phys. Chem., 95, 3358--3363, 1991 \n
!! \n
!! <b>Extended Charge Equilibration approach (EQeq)</b> \n
!! \n
!! Wilmer, C.E., Kim, K.C., Snurr, R.Q., <b>An Extended Charge Equilibration Method</b>,
!! J. Phys. Chem. Lett, 3, 2506--2511, 2012 
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
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

  !>
  !! @brief This routine calculates the partial charges using the EQeq approach
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] partialCharges REAL(8) ARRAY with the partial charges of the system
  !! @param [in] vertices Class with the information of the vertices
  !! @note For reference see: \n
  !! <b>Extended Charge Equilibration approach (EQeq)</b> \n
  !! \n
  !! Wilmer, C.E., Kim, K.C., Snurr, R.Q., <b>An Extended Charge Equilibration Method</b>,
  !! J. Phys. Chem. Lett, 3, 2506--2511, 2012 \n
  !! \n
  !! This routine needs the Idempotencials and electronegativities Matrices and resolves the system: \n
  !! \f{eqnarray*}{
  !!\left(\begin{bmatrix}
  !!       J_{11} & \cdots & J_{1N} \\[0.3em]
  !!       \vdots & \ddots & \cdots \\[0.3em]
  !!   J_{(N-1)1} & \cdots & J_{(N-1)N}
  !!     \end{bmatrix}-\begin{bmatrix}
  !!       J_{21} & \cdots & J_{2N} \\[0.3em]
  !!       \vdots & \ddots & \cdots \\[0.3em]
  !!       J_{N1} & \cdots & J_{NN}
  !!     \end{bmatrix}\right)\begin{bmatrix}
  !!       Q_{1} \\[0.3em]
  !!       \vdots \\[0.3em]
  !!       Q_{N}
  !!     \end{bmatrix}=\begin{bmatrix}
  !!       \chi_{2}-\chi_{1} \\[0.3em]
  !!       \vdots \\[0.3em]
  !!       \chi_{N}-\chi_{N-1}
  !!     \end{bmatrix}
  !! \f}
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

    
    ! write(*,"(T20,A,I)") "Info: ", info
    ! do i=1,vertices%numberOfVertices
    !    write(*,"(T20,F12.5)") partialCharges(i)
    ! end do

  end subroutine ChargesEQeq_getCharges

  !>
  !! @brief This routine calculates the idempotentials between two atoms of the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] i INTEGER first atom to evaluate
  !! @param [in] j INTEGER second atom to evaluate
  !! @param [in] vertices Class with the information of the vertices
  !! @return [out] idempotential REAL(8) return the idempotential of the atoms evaluated
  !! @note For reference see: \n
  !! <b>Extended Charge Equilibration approach (EQeq)</b> \n
  !! \n
  !! Wilmer, C.E., Kim, K.C., Snurr, R.Q., <b>An Extended Charge Equilibration Method</b>,
  !! J. Phys. Chem. Lett, 3, 2506--2511, 2012 \n
  !! \n
  !! \f[
  !! J_{ij} = \lambda\left(\frac{K}{2}\right)\left[\frac{1}{R_{ij}} + E_{0}(R_{ij})\right] 
  !! \f]
  !! \f[
  !! E_{0}(R_{ab}) = e^{-\left(\frac{J_{ab}R_{ab}}{K}\right)^2}\left(\frac{2J_{ab}}{K}-\frac{J_{ab}^{2}R_{ab}}{K^{2}}-\frac{1}{R_{ab}}\right)
  !! \f]
  !! \f[
  !! J_{ab} = \sqrt{J_{a}J_{b}}
  !! \f]
  !! where: \n
  !! - \f$J_{ij}\f$ is the idempotential
  !! - \f$R_{ij} = R_{ab}\f$ is distance between atoms
  !! - \f$E_{0}(R_{ij})\f$ is the orbital energy term
  !! - \f$J_{ab}\f$ is geometric mean of the chemical hardness
  !! - \f$J_{a}\f$ and \f$J_{b}\f$ are the hardnees of atoms
  !! - \f$K = 14.4 \f$ 
  !! - \f$\lambda = 1.2 \f$ 
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
