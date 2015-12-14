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
!! @brief   Modulo para calculo de integrales de atraccion repulsion con cargas puntuales.
!!
!! Este modulo define una seudoclase cuyos metodos devuelven la integral de atraccion
!! o repulsion entre una  particula cuantica descrita mediante una distribucion
!! gausiana sin normalizar  y una carga puntual localizada en \f$R_c\f$. El calculo de la
!! integral  no considera la carga de la particula puntual.
!!
!! \f[ (\bf{a} \mid A(0) \mid \bf{b}) = \int_{TE} {{\varphi(\bf{r};
!!  {\bf{\zeta}}_a, \bf{n_a,R_A})} {\frac{1}{\mid \bf{r-R_C} \mid
!!  }} {\varphi(\bf{r};{\zeta}_b,\bf{n_b,R_B})}},dr\f]
!!
!! Donde:
!!
!! <table>
!! <tr> <td> \f$ \zeta \f$ : <td> <dfn> Exponente orbital. </dfn>
!! <tr> <td> <b> r  </b> : <td> <dfn> Coordenas espaciales de la funcion. </dfn>
!! <tr> <td> \f$ n_a ,n_b \f$ : <td> <dfn> Indice de momento angular. </dfn>
!! <tr> <td> \f$ R_A , R_B \f$ : <td> <dfn> Origen de la funcion gaussiana cartesiana.</dfn>
!! <tr> <td> \f$ R_C  \f$ : <td> <dfn> origen de la particula puntual </dfn>
!! </table>
!!
!! Este tipo de integral corresponde a una integral de dos centros, calculada por
!! aproximacion numerica, utilizando la integral auxiliar (ver Math_):
!!
!!    \f[ F_m(U)= \int_{0}^{1} {t^{2m}e^{-Tt^2}}\,dt \f]
!!
!! La integral de repulsion-atraccion con cargas puntuales se calcula de
!! acuerdo metodo recursivo propuesto por Obara-Sayka. La expresion general
!! de la integral es:
!!
!! \f[({\bf{a + 1_i}} \parallel A(0) \parallel {\bf{b}})^{(m)} = \f]
!! \f[ (P_i -A_i) ({\bf{a}} \parallel A(0) \parallel {\bf{b}} )^{(m)}
!!  - (P_i -C_i)({\bf{a}} \parallel A(0) \parallel {\bf{b}} )^{(m+1)} \f]
!! \f[ + \frac{1}{2 \zeta} N_i(\bf{a}) \left\{ ({\bf{a-1_i}} \parallel A(0)
!! \parallel {\bf{b}})^{(m)} - ({\bf{a-1_i}} \parallel A(0) \parallel
!! {\bf{b}})^{(m+1)}  \right\}\f]
!! \f[ + \frac{1}{2 \zeta} N_i(\bf{b}) \left\{ ({\bf{a}} \parallel A(0)
!! \parallel {\bf{b-1_i}})^{(m)} - ({\bf{a}} \parallel A(0) \parallel
!! {\bf{b-1_i}})^{(m+1)}  \right\}\f]
!!
!! Donde (m) es un entero no negativo que hace referencia al orden de la
!! integral dentro de la recursion. Los parametros <b> P </b> y \f$ \zeta \f$ de la
!! expresion provienen del producto de dos funciones gaussianas.
!!
!! @author Fernando Posada (efposadac@unal.edu.co)
!!
!! <b> Fecha de creacion : </b> 20010-03-10
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-05-15 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapta al estandar de codificacion propuesto.
!!   - <tt> 2007-05-15 </tt>: Fernando Posada. ( efposadac@unal.edu.co )
!!        -# Reescribe el módulo.
!!   - <tt> 2010-06-05 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Implementa metodos para el calculo de integrales con L > 2
!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
module AttractionIntegrals_
  use Exception_
  use Math_
  use ContractedGaussian_
  implicit none
  

  !>
  !! Puntual particle atributes
  type, public :: pointCharge
     real(8) :: x
     real(8) :: y
     real(8) :: z
     real(8) :: charge
  end type pointCharge

  public ::  &
       AttractionIntegrals_computeShell, &
       AttractionIntegrals_computePrimitive

  private :: &
       AttractionIntegrals_obaraSaikaRecursion

contains

  !>
  !! Calculates point charge - quantum specie integral between two contractions (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2013.02.04: E.F.Posada: change for use in opints
  !! @return  output: attraction integral of a shell (all combinations)
  !! @version 1.0
  subroutine AttractionIntegrals_computeShell(contractedGaussianA, contractedGaussianB, point, npoints, integral)
    implicit none
    
    type(ContractedGaussian), intent(in) :: contractedGaussianA, contractedGaussianB
    type(pointCharge), intent(in), allocatable :: point(:)

    integer, intent(in) :: npoints
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
    
    integral = 0.0_8
    auxIntegral = 0.0_8

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

          call AttractionIntegrals_computePrimitive(am1, am2, nprim1, nprim2, npoints, A, B, exp1, exp2, coef1, coef2, nor1, nor2, point, auxintegral)


          auxIntegral = auxIntegral * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q)

          integral(m) = auxIntegral

       end do
    end do

  end subroutine AttractionIntegrals_computeShell

  !>
  !! @brief Evaluates attraction - repulsion integrals for any angular momentum
  !! @author Edwin Posada, 2010
  !! @return integral values for all shell (all possibles combinations if angular momentum index) (output)
  subroutine AttractionIntegrals_computePrimitive(angularMomentindexA, angularMomentindexB, lengthA, lengthB, &
       numberOfPointCharges, A, B, &
       orbitalExponentsA, orbitalExponentsB, &
       contractionCoefficientsA, contractionCoefficientsB, &
       normalizationConstantsA, normalizationConstantsB, &
       pointCharges, integralValue )
    implicit none
		
    integer, intent(in) :: angularMomentindexA(0:3), angularMomentindexB(0:3)
    integer, intent(in) :: lengthA, lengthB
    integer, intent(in) :: numberOfPointCharges
    real(8), intent(in) :: A(0:3), B(0:3)
    real(8), intent(in) :: orbitalExponentsA(0:lengthA) ,orbitalExponentsB(0:lengthB)
    real(8), intent(in) :: contractionCoefficientsA(0:lengthA), contractionCoefficientsB(0:lengthB)
    real(8), intent(in) :: normalizationConstantsA(0:lengthA), normalizationConstantsB(0:lengthB)
    type(pointCharge), intent(in) :: pointCharges(0:numberOfPointCharges-1)
    real(8), intent(inout) :: integralValue
		

    real(8), allocatable :: AI0(:,:,:)
    real(8) :: PA(0:3), PB(0:3), PC(0:3), P(0:3)
    real(8) :: auxExponentA, auxCoefficentA, auxConstantA
    real(8) :: auxExponentB, auxCoefficentB, auxConstantB
    real(8) :: commonPreFactor
    real(8) :: AB2
    real(8) :: zeta, zetaInv
    integer :: angularMomentA, angularMomentB
    integer :: izm, iym, ixm
    integer :: jzm, jym, jxm
    integer :: maxIndex
    integer :: p1, p2
    integer :: atom
    integer :: maxAngularMoment, sumAngularMoment
    integer :: indexI
    integer :: indexJ
    integer :: i, j

    integralValue = 0.0_8

    angularMomentA = sum(angularMomentindexA)
    angularMomentB = sum(angularMomentindexB)

    maxAngularMoment = max(angularMomentA, angularMomentB) + 1

    maxIndex = (maxAngularMoment-1)*maxAngularMoment*maxAngularMoment+1

    if(allocated(AI0))deallocate(AI0)
    allocate(AI0(0:maxIndex, 0:maxIndex, 0:2*maxAngularMoment+1))

    AI0 = 0.0_8

    AB2 = 0.0_8
    AB2 = AB2 + (A(0) - B(0)) * (A(0) - B(0))
    AB2 = AB2 + (A(1) - B(1)) * (A(1) - B(1))
    AB2 = AB2 + (A(2) - B(2)) * (A(2) - B(2))

    izm = 1
    iym = angularMomentA+1
    ixm = iym*iym

    jzm = 1
    jym = angularMomentB+1
    jxm = jym*jym

    do p1=0, lengthA-1
       auxExponentA = orbitalExponentsA(p1)
       auxCoefficentA = contractionCoefficientsA(p1)
       auxConstantA = normalizationConstantsA(p1)
       do p2=0, lengthB -1
          auxExponentB = orbitalExponentsB(p2)
          auxCoefficentB = contractionCoefficientsB(p2)
          auxConstantB = normalizationConstantsB(p2)
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

          commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*zetaInv) * sqrt(Math_PI*zetaInv) * Math_PI * zetaInv * auxCoefficentA * auxCoefficentB * auxConstantA * auxConstantB
          ! write(*,*)"fragmentos y length a y b",numberOfPointCharges,lengthA,lengthB

          do atom = 0, numberOfPointCharges - 1

             PC(0) = P(0) - pointCharges(atom)%x
             PC(1) = P(1) - pointCharges(atom)%y
             PC(2) = P(2) - pointCharges(atom)%z

             sumAngularMoment = angularMomentA + angularMomentB + 1

             call AttractionIntegrals_obaraSaikaRecursion(AI0,PA,PB,PC,zeta,sumAngularMoment,angularMomentA,angularMomentB)

             indexI = angularMomentindexA(2)*izm + angularMomentindexA(1)*iym + angularMomentindexA(0)*ixm

             indexJ = angularMomentindexB(2)*jzm + angularMomentindexB(1)*jym + angularMomentindexB(0)*jxm

             integralValue = integralValue - AI0(indexI,indexJ,0) * pointCharges(atom)%charge * commonPreFactor

          end do
          ! write(*,*) "se ha llamado obara-saika ",atom," veces"
       end do
    end do
		! write(*,*)"finaliza_computePrimitives"

  end subroutine AttractionIntegrals_computePrimitive

  !> @brief Obara-Saika recursion for nucleo-electron attraction integrals. Supports all angular momentum numbers
  !! @author E. F. Posada, 2010
  !! @version 1.0  
  subroutine AttractionIntegrals_obaraSaikaRecursion(AI0, PA, PB, PC, zeta, sumAngularMoment, angularMomentA, angularMomentB)
    implicit none
    real(8), intent(inout), allocatable :: AI0(:,:,:)
    real(8), intent(in) :: PA(0:3)
    real(8), intent(in) :: PB(0:3)
    real(8), intent(in) :: PC(0:3)
    real(8), intent(in) :: zeta
    integer, intent(in) :: sumAngularMoment
    integer, intent(in) :: angularMomentA
    integer, intent(in) :: angularMomentB

    real(8), dimension(0: sumAngularMoment) :: F
    real(8) :: twoZetaInv
    real(8) :: tmp
    real(8) :: u
    integer :: a, b, m
    integer :: izm, iym, ixm
    integer :: jzm, jym, jxm
    integer :: ix,iy,iz,jx,jy,jz
    integer :: iind,jind

    izm = 1
    iym = angularMomentA + 1
    ixm = iym * iym
    jzm = 1
    jym = angularMomentB + 1
    jxm = jym * jym
    twoZetaInv = 1 / (2 * zeta)
    tmp = sqrt(zeta)*(2.0_8/Math_SQRT_PI)
    u = zeta*(PC(0) * PC(0) + PC(1) * PC(1) + PC(2) * PC(2))

    call Math_fgamma0(sumAngularMoment,u,F)

    do m = 0, sumAngularMoment
       AI0(0, 0, m) = tmp * F(m)
    end do

    !! Upward recursion in j with i=0
    do b = 1, angularMomentB
       do jx = 0, b
          do jy=0, b - jx
             jz = b - jx - jy
             jind = jx * jxm + jy * jym + jz * jzm
             if (jz > 0) then
                do m=0, sumAngularMoment - b	!! Electrostatic potential integrals
                   AI0(0,jind,m) = PB(2)*AI0(0,jind-jzm,m) - PC(2)*AI0(0, jind-jzm, m+1)
                end do

                if (jz > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jz-1)*(AI0(0,jind-2*jzm,m) - AI0(0,jind-2*jzm,m+1))
                   end do
                end if

             else if (jy > 0) then
                do m=0, sumAngularMoment-b
                   AI0(0,jind,m) = PB(1)*AI0(0,jind-jym,m) -	PC(1)*AI0(0,jind-jym,m+1)
                end do

                if (jy > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jy-1)*(AI0(0,jind-2*jym,m) - AI0(0,jind-2*jym,m+1))
                   end do
                end if

             else if (jx > 0) then
                do m=0, sumAngularMoment-b
                   AI0(0,jind,m) = PB(0)*AI0(0,jind-jxm,m) -	PC(0)*AI0(0,jind-jxm,m+1)
                end do

                if (jx > 1) then
                   do m=0, sumAngularMoment-b
                      AI0(0,jind,m) = AI0(0,jind,m) + twoZetaInv*(jx-1)*(AI0(0,jind-2*jxm,m) - AI0(0,jind-2*jxm,m+1))
                   end do
                end if

             else
                call AttractionIntegrals_exception( ERROR, "AttractionIntegrals in obaraSaikaRecursion Function", &
                     "There's some error in the obaraSaika algorithm")
             end if
          end do
       end do
    end do

    !! The following fragment cannot be vectorized easily, I guess :-)
    !! Upward recursion in i with all possible j's
    do b=0, angularMomentB
       do jx=0, b
          do jy=0, b-jx
             jz = b-jx-jy
             jind = jx*jxm + jy*jym + jz*jzm
             do a=1, angularMomentA
                do ix=0, a
                   do iy=0, a-ix
                      iz = a-ix-iy
                      iind = ix*ixm + iy*iym + iz*izm
                      if (iz > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(2)*AI0(iind-izm,jind,m) - PC(2)*AI0(iind-izm,jind,m+1)
                         end do

                         if (iz > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(iz-1)*(AI0(iind-2*izm,jind,m) - AI0(iind-2*izm,jind,m+1))
                            end do
                         end if

                         if (jz > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jz*(AI0(iind-izm,jind-jzm,m) - AI0(iind-izm,jind-jzm,m+1))
                            end do
                         end if

                      else if (iy > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(1)*AI0(iind-iym,jind,m) - PC(1)*AI0(iind-iym,jind,m+1)
                         end do

                         if (iy > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(iy-1)*(AI0(iind-2*iym,jind,m) - AI0(iind-2*iym,jind,m+1))
                            end do
                         end if

                         if (jy > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jy*(AI0(iind-iym,jind-jym,m) - AI0(iind-iym,jind-jym,m+1))
                            end do
                         end if

                      else if (ix > 0) then
                         do m=0, sumAngularMoment-a-b
                            AI0(iind,jind,m) = PA(0)*AI0(iind-ixm,jind,m) - PC(0)*AI0(iind-ixm,jind,m+1)
                         end do

                         if (ix > 1) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*(ix-1)*(AI0(iind-2*ixm,jind,m) - AI0(iind-2*ixm,jind,m+1))
                            end do
                         end if

                         if (jx > 0) then
                            do m=0, sumAngularMoment-a-b
                               AI0(iind,jind,m) = AI0(iind,jind,m) + twoZetaInv*jx*(AI0(iind-ixm,jind-jxm,m) - AI0(iind-ixm,jind-jxm,m+1))
                            end do
                         end if
                      else
                         call AttractionIntegrals_exception( ERROR, "AttractionIntegrals in obaraSaikaRecursion Function", &
                              "There's some error in the obaraSaika algorithm")
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

  end subroutine AttractionIntegrals_obaraSaikaRecursion

  !>
  !! @brief  Handle exceptions of this module
  subroutine AttractionIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine AttractionIntegrals_exception

end module AttractionIntegrals_
