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
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

!>
!! @brief Modulo para calculo de integrales de energia cinetica.
!!
!! Este modulo define una seudoclase cuyos metodos devuelven la integral de
!! energia cinetica  para una particula descrita mediante una distribucion
!! gausiana de sin normalizar. La integral es de la forma:
!!
!! \f[(\bf{a} \parallel T \parallel \bf{b}) = \int_{TE} {{\varphi(\bf{r};
!!  {\zeta}_a, \bf{n_a,R_A})}{(-\frac{1}{2}) {\nabla}^2} {\varphi(\bf{r};
!!  {\zeta}_b , \bf{n_b,R_B})}}\,dr\f]
!!
!! Este tipo de integral corresponde a una integral de dos centros calalada de
!! forma cerrada.
!!
!! Donde:
!!
!! <table>
!! <tr> <td> \f$ \zeta \f$ : <td> <dfn> exponente orbital. </dfn>
!! <tr> <td> \f$ r \f$ : <td>  <dfn> coordenas espaciales de la funcion. </dfn>
!! <tr> <td> \f$ n_a , n_b \f$ : <td> <dfn> indice de momento angular.</dfn>
!! <tr> <td> \f$ R_A , R_B \f$ : <td> <dfn>origen de las gausianas </dfn>
!! </table>
!!
!! La integral de energia cinetica para la particula dada se calcula de acuerdo
!! metodo recursivo propuesto por Obara-Sayka, cuya expresion general para la
!! integral es:
!!
!! \f[({\bf{a + 1_i}} \parallel T \parallel {\bf{b}}) = \f]
!! \f[ (P_i -A_i) ({\bf{a}} \parallel T \parallel {\bf{b}} ) + \frac{1}
!! {2 \zeta} N_i(\bf{a}) ({\bf{a-1_i}} \parallel T \parallel {\bf{b}} )
!! + \frac{1}{2 \zeta} N_i(\bf{b}) + ({\bf{a}} \parallel T \parallel
!! {\bf{b-1_i}}) + 2  \xi \left[({\bf{a - 1_i}} \parallel {\bf{b}})
!! \frac{1}{2 {\zeta}_a} N_i(\bf{a}) ({\bf{a - 1_i}} \parallel {\bf{b}} )
!!  \right]\f]
!!
!! Los parametros <b> P </b> , \f$ \xi \f$ y \f$ \zeta\f$ de la expresion provienen del
!! producto de dos funciones gaussianas.
!!
!! @author Edwin Posada
!!
!! <b> Fecha de creacion : </b> 2010-03-08
!!
!! <b> Historial de modificaciones: </b>
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin
!!
module KineticIntegrals_
  use Exception_
  use ContractedGaussian_
  use Math_
  implicit none

  
  public :: &
       kineticIntegrals_computeShell, &
       kineticIntegrals_computePrimitive

  private :: &
       KineticIntegrals_obaraSaikaRecursion

contains
  
  !>
  !! Calculates kinetic integral between two contractions (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2013.02.04: E.F.Posada: change for use in opints
  !! @return  output: kinetic integral of a shell (all combinations)
  !! @version 1.0
  subroutine kineticIntegrals_computeShell(contractedGaussianA, contractedGaussianB, integral)
    implicit none
    
    type(ContractedGaussian), intent(in) :: contractedGaussianA, contractedGaussianB
    real(8), intent(inout) :: integral(contractedGaussianA%numCartesianOrbital * contractedGaussianB%numCartesianOrbital)

    integer ::  am1(0:3)
    integer ::  am2(0:3)
    integer ::  nprim1
    integer ::  nprim2
    real(8) ::  A(0:3)
    real(8) ::  B(0:3)
    real(8) ::  exp1(0:contractedGaussianA%length)
    real(8) ::  exp2(0:contractedGaussianB%length)
    real(8) ::  coef1(0:contractedGaussianA%length)
    real(8) ::  coef2(0:contractedGaussianB%length)
    real(8) ::  nor1(0:contractedGaussianA%length)
    real(8) ::  nor2(0:contractedGaussianB%length)
    real(8) :: auxIntegral
    integer, allocatable :: angularMomentIndexA(:,:)
    integer, allocatable :: angularMomentIndexB(:,:)
    integer ::  i, m, p, q
    
    if(allocated(angularMomentIndexA)) deallocate(angularMomentIndexA)
    if(allocated(angularMomentIndexB)) deallocate(angularMomentIndexB)

    allocate(angularMomentIndexA(3, contractedGaussianA%numCartesianOrbital))
    allocate(angularMomentIndexB(3, contractedGaussianB%numCartesianOrbital))

    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexA, contractedGaussianA)
    call contractedGaussian_getAllAngularMomentIndex(angularMomentIndexB, contractedGaussianB)

    nprim1 = contractedGaussianA%length
    A(0) = contractedGaussianA%origin(1)
    A(1) = contractedGaussianA%origin(2)
    A(2) = contractedGaussianA%origin(3)
    coef1(0:nprim1-1) =  contractedGaussianA%contractionCoefficients(1:nprim1)


    nprim2 = contractedGaussianB%length
    B(0) = contractedGaussianB%origin(1)
    B(1) = contractedGaussianB%origin(2)
    B(2) = contractedGaussianB%origin(3)
    coef2(0:nprim2-1) =  contractedGaussianB%contractionCoefficients(1:nprim2)
    
    m = 0

    do p = 1, contractedGaussianA%numcartesianOrbital
       do q = 1, contractedGaussianB%numcartesianOrbital

          m = m + 1

          exp1(0:nprim1-1) = contractedGaussianA%orbitalExponents(1:nprim1)
          nor1(0:nprim1-1) = contractedGaussianA%primNormalization(1:nprim1,p)

          exp2(0:nprim2-1) = contractedGaussianB%orbitalExponents(1:nprim2)
          nor2(0:nprim2-1) = contractedGaussianB%primNormalization(1:nprim2,q)
             
          am1 = 0
          am2 = 0

          am1(0:2) = angularMomentIndexA(1:3, p)
          am2(0:2) = angularMomentIndexB(1:3, q)


          call KineticIntegrals_computePrimitive(am1, am2, nprim1, nprim2, A, B, exp1, exp2, coef1, coef2, nor1, nor2, auxIntegral)

          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q)
          
          integral(m) = auxIntegral

       end do
    end do
    
  end subroutine KineticIntegrals_computeShell
  
  !>
  !! @brief Evalua integrales overlap para cualquier momento angular
  !! @author Edwin Posada, 2010
  !! @return devuelve los valores de integrales de atraccion (output)
  !! @version 1.0
  subroutine kineticIntegrals_computePrimitive(angularMomentIndexA, angularMomentIndexB, lengthA, lengthB, A, B, orbitalExponentsA, orbitalExponentsB, contractionCoefficientsA, contractionCoefficientsB, normalizationConstantA, normalizationConstantB, integralValue)
    implicit none

    integer, intent(in) :: angularMomentIndexA(0:3), angularMomentIndexB(0:3)
    integer, intent(in) :: lengthA, lengthB
    real(8), intent(in) :: A(0:3), B(0:3)
    real(8), intent(in) :: orbitalExponentsA(0:lengthA), orbitalExponentsB(0:lengthB)
    real(8), intent(in) :: contractionCoefficientsA(0:lengthA), contractionCoefficientsB(0:lengthB)
    real(8), intent(in) :: normalizationConstantA(0:lengthA), normalizationConstantB(0:lengthB)
    real(8), intent(out) :: integralValue

    real(8), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: AB2
    real(8) :: auxExponentA, auxCoefficientA, auxConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxConstantB
    real(8) :: zeta, zetaInv
    real(8) :: PA(0:3), PB(0:3), P(0:3)
    real(8) :: commonPreFactor
    real(8) :: x0, y0, z0
    real(8) :: I1, I2, I3, I4
    real(8) :: Ix, Iy, Iz

    integer :: angularMomentA, angularMomentB
    integer :: maxAngularMoment
    integer :: p1, p2
    integer :: ao12
    integer :: ii, jj, kk, ll
    integer :: l1, m1, n1
    integer :: l2, m2, n2

    integralValue = 0.0_8

    angularMomentA = sum(angularMomentIndexA)
    angularMomentB = sum(angularMomentIndexB)

    maxAngularMoment = max(angularMomentA, angularMomentB) + 1

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2), y(0:maxAngularMoment+2, 0:maxAngularMoment+2), z(0:maxAngularMoment+2, 0:maxAngularMoment+2))

    AB2 = 0.0_8
    AB2 = AB2 + (A(0) - B(0)) * (A(0) - B(0))
    AB2 = AB2 + (A(1) - B(1)) * (A(1) - B(1))
    AB2 = AB2 + (A(2) - B(2)) * (A(2) - B(2))

    do p1=0, lengthA - 1
       auxExponentA = orbitalExponentsA(p1)
       auxCoefficientA = contractionCoefficientsA(p1)
       auxConstantA = normalizationConstantA(p1)
       do p2=0, lengthB - 1
          auxExponentB = orbitalExponentsB(p2)
          auxCoefficientB = contractionCoefficientsB(p2)
          auxConstantB = normalizationConstantB(p2)
          zeta = auxExponentA + auxExponentB
          zetaInv = 1.0/zeta

          P(0) = (auxExponentA*A(0) + auxExponentB*B(0))*zetaInv
          P(1) = (auxExponentA*A(1) + auxExponentB*B(1))*zetaInv
          P(2) = (auxExponentA*A(2) + auxExponentB*B(2))*zetaInv
          PA(0) = P(0) - A(0)
          PA(1) = P(1) - A(1)
          PA(2) = P(2) - A(2)
          PB(0) = P(0) - B(0)
          PB(1) = P(1) - B(1)
          PB(2) = P(2) - B(2)

          commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*zetaInv) * sqrt(Math_PI*zetaInv) * Math_PI * zetaInv * auxCoefficientA * auxCoefficientB * auxConstantA * auxConstantB

          !! recursion
          call KineticIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, zeta, angularMomentA+2, angularMomentB+2)

          if (angularMomentIndexA(0) == 0 .or. angularMomentIndexB(0) == 0) then
             I1 = 0.0_8
          else
             I1 = x(angularMomentIndexA(0)-1,angularMomentIndexB(0)-1) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          I2 = x(angularMomentIndexA(0)+1,angularMomentIndexB(0)+1) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor

          if (angularMomentIndexB(0) == 0) then
             I3 = 0.0_8
          else
             I3 = x(angularMomentIndexA(0)+1,angularMomentIndexB(0)-1) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          if (angularMomentIndexA(0) == 0) then
             I4 = 0.0_8
          else
             I4 = x(angularMomentIndexA(0)-1,angularMomentIndexB(0)+1) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          Ix = 0.5 * angularMomentIndexA(0) * angularMomentIndexB(0) * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * angularMomentIndexB(0) * I3 - angularMomentIndexA(0) * auxExponentB * I4

          if (angularMomentIndexA(1) == 0 .or. angularMomentIndexB(1) == 0) then
             I1 = 0.0_8
          else
             I1 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1)-1,angularMomentIndexB(1)-1) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          I2 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1)+1,angularMomentIndexB(1)+1) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor

          if (angularMomentIndexB(1) == 0) then
             I3 = 0.0_8
          else
             I3 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1)+1,angularMomentIndexB(1)-1) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          if (angularMomentIndexA(1) == 0) then
             I4 = 0.0_8
          else
             I4 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1)-1,angularMomentIndexB(1)+1) * z(angularMomentIndexA(2),angularMomentIndexB(2)) * commonPreFactor
          end if

          Iy = 0.5 * angularMomentIndexA(1) * angularMomentIndexB(1) * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * angularMomentIndexB(1) * I3 - angularMomentIndexA(1) * auxExponentB * I4

          if (angularMomentIndexA(2) == 0 .or. angularMomentIndexB(2) == 0) then
             I1 = 0.0_8
          else
             I1 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2)-1,angularMomentIndexB(2)-1) * commonPreFactor
          end if

          I2 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2)+1,angularMomentIndexB(2)+1) * commonPreFactor

          if (angularMomentIndexB(2) == 0) then
             I3 = 0.0_8
          else
             I3 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2)+1,angularMomentIndexB(2)-1) * commonPreFactor
          end if

          if (angularMomentIndexA(2) == 0) then
             I4 = 0.0_8
          else
             I4 = x(angularMomentIndexA(0),angularMomentIndexB(0)) * y(angularMomentIndexA(1),angularMomentIndexB(1)) * z(angularMomentIndexA(2)-1,angularMomentIndexB(2)+1) * commonPreFactor
          end if

          Iz = 0.5 * angularMomentIndexA(2) * angularMomentIndexB(2) * I1 + 2.0 * auxExponentA * auxExponentB * I2 - auxExponentA * angularMomentIndexB(2) * I3 - angularMomentIndexA(2) * auxExponentB * I4

          integralValue = integralValue + (Ix + Iy + Iz)

       end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)

  end subroutine kineticIntegrals_computePrimitive

  !>
  !! @brief Implementacion de la recursion propuesta por Obara - Saika para el caso de integrales de energia cinetica
  !! @author E. F. Posada, 2010
  !! @version 1.0
  subroutine KineticIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, zeta, angularMomentIndexA, angularMomentIndexB)
    implicit none

    real(8), intent(inout), allocatable :: x(:,:), y(:,:), z(:,:)
    real(8), intent(in) :: PA(0:3), PB(0:3)
    real(8), intent(in) :: zeta
    integer, intent(in) :: angularMomentIndexA, angularMomentIndexB

    real(8) :: twoZetaInv
    integer :: i, j, k

    twoZetaInv = 1_8/(2_8*zeta)

    x(0,0) = 1.0_8
    y(0,0) = 1.0_8
    z(0,0) = 1.0_8

    !! Upward recursion in j for i=0
    x(0,1) = PB(0)
    y(0,1) = PB(1)
    z(0,1) = PB(2)

    do j=1, angularMomentIndexB -1
       x(0,j+1) = PB(0)*x(0,j)
       y(0,j+1) = PB(1)*y(0,j)
       z(0,j+1) = PB(2)*z(0,j)
       x(0,j+1) = x(0,j+1) + j*twoZetaInv*x(0,j-1)
       y(0,j+1) = y(0,j+1) + j*twoZetaInv*y(0,j-1)
       z(0,j+1) = z(0,j+1) + j*twoZetaInv*z(0,j-1)
    end do

    !! Upward recursion in i for all j's
    x(1,0) = PA(0)
    y(1,0) = PA(1)
    z(1,0) = PA(2)

    do j=1, angularMomentIndexB
       x(1,j) = PA(0)*x(0,j)
       y(1,j) = PA(1)*y(0,j)
       z(1,j) = PA(2)*z(0,j)
       x(1,j) = x(1,j) + j*twoZetaInv*x(0,j-1)
       y(1,j) = y(1,j) + j*twoZetaInv*y(0,j-1)
       z(1,j) = z(1,j) + j*twoZetaInv*z(0,j-1)
    end do

    do i=1, angularMomentIndexA - 1
       x(i+1,0) = PA(0)*x(i,0)
       y(i+1,0) = PA(1)*y(i,0)
       z(i+1,0) = PA(2)*z(i,0)
       x(i+1,0) = x(i+1,0) + i*twoZetaInv*x(i-1,0)
       y(i+1,0) = y(i+1,0) + i*twoZetaInv*y(i-1,0)
       z(i+1,0) = z(i+1,0) + i*twoZetaInv*z(i-1,0)
       do j=1, angularMomentIndexB
          x(i+1,j) = PA(0)*x(i,j)
          y(i+1,j) = PA(1)*y(i,j)
          z(i+1,j) = PA(2)*z(i,j)
          x(i+1,j) = x(i+1,j) + i*twoZetaInv*x(i-1,j)
          y(i+1,j) = y(i+1,j) + i*twoZetaInv*y(i-1,j)
          z(i+1,j) = z(i+1,j) + i*twoZetaInv*z(i-1,j)
          x(i+1,j) = x(i+1,j) + j*twoZetaInv*x(i,j-1)
          y(i+1,j) = y(i+1,j) + j*twoZetaInv*y(i,j-1)
          z(i+1,j) = z(i+1,j) + j*twoZetaInv*z(i,j-1)
       end do
    end do

  end subroutine KineticIntegrals_obaraSaikaRecursion

  !>
  !! @brief  Maneja excepciones de la clase
  !! @author Sergio Gonzalez
  subroutine KineticIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine KineticIntegrals_exception

end module KineticIntegrals_
