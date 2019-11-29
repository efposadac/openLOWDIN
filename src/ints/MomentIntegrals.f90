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
!! @brief Modulo para calculo de integrales de momento
!!
!! Este modulo contiene los algoritmos necesarios para la evaluacion de integrales de
!! momento entre pares de funciones gaussianas primitivas (PrimitiveGaussian_),
!! sin normalizar.
!!
!! Las integrales son de la forma:
!!
!! \f[ (\bf{a} \mid {\mathbf{R}(\mu)} \mid \mathbf{b}) = \int_{TE} {{\varphi(\bf{r};{\zeta}_a,
!! \bf{a,A})} {{\mathbf{R}(\mu)}} {\varphi(\bf{r};{\zeta}_b,\mathbf{b,B})}},dr \f]
!!
!! Donde:
!!
!! <table>
!! <tr> <td> \f$ \zeta \f$ : <td> <dfn> exponente orbital. </dfn>
!! <tr> <td> <b> r </b> : <td> <dfn> coordenas espaciales de la funcion. </dfn>
!! <tr> <td> <b> n </b> : <td> <dfn> indice de momento angular. </dfn>
!! <tr> <td> <b> R </b> : <td> <dfn> origen de la funcion gaussiana cartesiana. </dfn>
!! </table>
!!
!! @author E. F. Posada
!!
!! <b> Fecha de creacion : </b> 2010-03-11
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
!!
module MomentIntegrals_
  use Exception_
  use Math_
  use ContractedGaussian_
  implicit none

  public ::  &
       MomentIntegrals_computeShell, &
       MomentIntegrals_computePrimitive

  private :: &
       MomentIntegrals_obaraSaikaRecursion

contains

  !>
  !! @brief Calculates moment integral between two contractions with respect to one component: x, y or z (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2011.02.04: E.F.Posada: chage for usage on opints
  !! @return  output: overlap integral of a shell (all combinations)
  !! @version 1.0
  subroutine MomentIntegrals_computeShell(contractedGaussianA, contractedGaussianB, originRC, component, integral)
    implicit none
    
    type(ContractedGaussian), intent(in) :: contractedGaussianA, contractedGaussianB
    real(8) :: originRC(3)
    integer :: component
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

          call MomentIntegrals_computePrimitive(am1, am2, nprim1, nprim2, A, B, exp1, exp2, coef1, coef2, nor1, nor2, &
               originRC, component, auxIntegral)

          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q)
          
          integral(m) = auxIntegral

       end do
    end do

  end subroutine MomentIntegrals_computeShell

  
  !>
  !!@brief Evalua integrales overlap para cualquier momento angular
  !!@author Edwin Posada, 2010
  !!@return devuelve los valores de integrales de overlap para una capa o individual
  !!@version 1.0
  subroutine MomentIntegrals_computePrimitive(angularMomentIndexA, angularMomentIndexB, lengthA, lengthB, A, B, &
       orbitalExponentsA, orbitalExponentsB, &
       contractionCoefficientsA, contractionCoefficientsB, &
       normalizationConstantsA, normalizationConstantsB, &
       originRC, component, integralValue)
    implicit none

    integer, intent(in) :: angularMomentIndexA(0:3), angularMomentIndexB(0:3)
    integer, intent(in) :: lengthA, lengthB
    real(8), intent(in) :: A(0:3), B(0:3)
    real(8), intent(in) :: orbitalExponentsA(0:lengthA), orbitalExponentsB(0:lengthB)
    real(8), intent(in) :: contractionCoefficientsA(0:lengthA), contractionCoefficientsB(0:lengthB)
    real(8), intent(in) :: normalizationConstantsA(0:lengthA), normalizationConstantsB(0:lengthB)
    real(8), intent(in) :: originRC(3)
    integer, intent(in) :: component
    real(8), intent(out):: integralValue

    real(8), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: AB2
    real(8) :: auxExponentA, auxCoefficientA, auxConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxConstantB
    real(8) :: gamma, gammaInv
    real(8) :: PA(0:3), PB(0:3), P(0:3)
    real(8) :: commonPreFactor
    real(8) :: x00, y00, z00
    real(8) :: x01, y01, z01
    real(8) :: x02, y02, z02
    real(8) :: x10, y10, z10
    real(8) :: x11, y11, z11

    integer :: angularMomentA, angularMomentB
    integer :: maxAngularMoment
    integer :: p1, p2 !< iteradores

    angularMomentA = sum(angularMomentIndexA)
    angularMomentB = sum(angularMomentIndexB)

    integralValue = 0.0_8

    maxAngularMoment = max(angularMomentA, angularMomentB) + 1

    allocate(x(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(y(0:maxAngularMoment+2, 0:maxAngularMoment+2))
    allocate(z(0:maxAngularMoment+2, 0:maxAngularMoment+2))

    x = 0.0_8
    y = 0.0_8
    z = 0.0_8

    AB2 = 0.0_8
    AB2 = AB2 + (A(0) - B(0)) * (A(0) - B(0))
    AB2 = AB2 + (A(1) - B(1)) * (A(1) - B(1))
    AB2 = AB2 + (A(2) - B(2)) * (A(2) - B(2))

    do p1=0, lengthA - 1
       auxExponentA = orbitalExponentsA(p1)
       auxCoefficientA = contractionCoefficientsA(p1)
       auxConstantA = normalizationConstantsA(p1)
       do p2=0, lengthB - 1
          auxExponentB = orbitalExponentsB(p2)
          auxCoefficientB = contractionCoefficientsB(p2)
          auxConstantB = normalizationConstantsB(p2)
          gamma = auxExponentA + auxExponentB
          gammaInv = 1.0/gamma

          P(0) = (auxExponentA*A(0) + auxExponentB*B(0))*gammaInv
          P(1) = (auxExponentA*A(1) + auxExponentB*B(1))*gammaInv
          P(2) = (auxExponentA*A(2) + auxExponentB*B(2))*gammaInv
          PA(0) = P(0) - A(0)
          PA(1) = P(1) - A(1)
          PA(2) = P(2) - A(2)
          PB(0) = P(0) - B(0)
          PB(1) = P(1) - B(1)
          PB(2) = P(2) - B(2)

          commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*gammaInv) &
               * sqrt(Math_PI*gammaInv) * Math_PI * gammaInv &
               * auxCoefficientA * auxCoefficientB * auxConstantA * auxConstantB

          !! recursion
          call MomentIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+2, angularMomentB+2)

          x00 = x(angularMomentIndexA(0),angularMomentIndexB(0))
          y00 = y(angularMomentIndexA(1),angularMomentIndexB(1))
          z00 = z(angularMomentIndexA(2),angularMomentIndexB(2))

          x01 = x(angularMomentIndexA(0),angularMomentIndexB(0)+1)
          y01 = y(angularMomentIndexA(1),angularMomentIndexB(1)+1)
          z01 = z(angularMomentIndexA(2),angularMomentIndexB(2)+1)

          x02 = x(angularMomentIndexA(0),angularMomentIndexB(0)+2)
          y02 = y(angularMomentIndexA(1),angularMomentIndexB(1)+2)
          z02 = z(angularMomentIndexA(2),angularMomentIndexB(2)+2)

          x10 = x(angularMomentIndexA(0)+1,angularMomentIndexB(0))
          y10 = y(angularMomentIndexA(1)+1,angularMomentIndexB(1))
          z10 = z(angularMomentIndexA(2)+1,angularMomentIndexB(2))

          x11 = x(angularMomentIndexA(0)+1,angularMomentIndexB(0)+1)
          y11 = y(angularMomentIndexA(1)+1,angularMomentIndexB(1)+1)
          z11 = z(angularMomentIndexA(2)+1,angularMomentIndexB(2)+1)

          select case (component)

          !! Dipole
          case(1) !X
             integralValue = integralValue + (commonPreFactor*(x01+x00*(B(0)-originRC(1)))*y00*z00)
          case(2) !Y
             integralValue = integralValue + (commonPreFactor*x00*(y01+y00*(B(1)-originRC(2)))*z00)
          case(3) !Z
             integralValue = integralValue + (commonPreFactor*x00*y00*(z01+z00*(B(2)-originRC(3))))

          !! Quadrupole
          case(4) !XX
             integralValue = integralValue + (1.0/1.0)*(commonPreFactor*(x02+x00*(B(0)-originRC(1)))*y00*z00)
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*x00*(y02+y00*(B(1)-originRC(2)))*z00)
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*x00*y00*(z02+z00*(B(2)-originRC(3))))
          case(5) !YY
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*(x02+x00*(B(0)-originRC(1)))*y00*z00)
             integralValue = integralValue + (1.0/1.0)*(commonPreFactor*x00*(y02+y00*(B(1)-originRC(2)))*z00)
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*x00*y00*(z02+z00*(B(2)-originRC(3))))
          case(6) !ZZ
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*(x02+x00*(B(0)-originRC(1)))*y00*z00)
             integralValue = integralValue - (1.0/2.0)*(commonPreFactor*x00*(y02+y00*(B(1)-originRC(2)))*z00)
             integralValue = integralValue + (1.0/1.0)*(commonPreFactor*x00*y00*(z02+z00*(B(2)-originRC(3))))
          case(7) !XY
             integralValue = integralValue + (3.0/2.0)*(commonPreFactor*(x01+x00*(B(0)-originRC(1)))*(y01+y00*(B(1)-originRC(2)))*z00)
          case(8) !XZ
             integralValue = integralValue + (3.0/2.0)*(commonPreFactor*(x01+x00*(B(0)-originRC(1)))*y00*(z01+z00*(B(2)-originRC(3))))
          case(9) !YZ
             integralValue = integralValue + (3.0/2.0)*(commonPreFactor*x00*(y01+y00*(B(1)-originRC(2)))*(z01+z00*(B(2)-originRC(3))))
          end select

       end do
    end do

    deallocate(x)
    deallocate(y)
    deallocate(z)

  end subroutine MomentIntegrals_computePrimitive

  !>
  !! @brief Implementacion de la recursion propuesta por Obara - Saika para el caso de integrales de momento
  subroutine MomentIntegrals_obaraSaikaRecursion(x, y, z, PA, PB, zeta, angularMomentIndexA, angularMomentIndexB)
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

  end subroutine MomentIntegrals_obaraSaikaRecursion

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine MomentIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine MomentIntegrals_exception

end module MomentIntegrals_
