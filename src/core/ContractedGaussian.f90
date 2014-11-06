!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief   Modulo para calculo de integrales entre contracciones de funciones gausianas
!!
!! Este modulo define una seudoclase para la construccion de contracciones de funciones
!! gausianas primitivas. Ademas incluye los metodos necesarios para el calculos de
!! integrales moleculares entre las contracciones definidas. Las funciones contraidas
!! tienen la siguiente forma:
!!
!! \f[ \chi(\bf{r-R_i}) = \sum_{i=1}^{L} {C_i{\phi (\bf{r;n_i},{\zeta}_i ,\bf{R_i})}} \f]
!!
!! y la integral de una particula asociada:
!!
!! \f[  \int_{TE}{ {\chi}^{*} \hat{O} \chi } \f]
!! Donde:
!!
!! <table>
!! <tr> <td> \f$ \chi \f$ : <td> <dfn> gausiana contraida. </dfn>
!! <tr> <td> <b> L </b> : <td><dfn> Longitud de la contraccion. </dfn>
!! <tr> <td> \f$ C_i \f$ :<td> <dfn> Coeficiente de contraccion de la i-enesima primitiva.</dfn>
!! <tr> <td> \f$ R_i \f$ :<td> <dfn> Origen de la  i-enesima gausiana primitiva.</dfn>
!! <tr> <td> \f$ n_i \f$ :<td> <dfn> Indice de momento anngular de la  i-enesima primitiva.</dfn>
!! <tr> <td> \f$ {\zeta}_i \f$ : <td> <dfn> Exponente orbital de la  i-enesima primitiva.</dfn>
!! </table>
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2007-02-06
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-07-12 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapta al estandar de codificacion propuesto.
!!   - <tt> 2010-09-25 </tt>: Edwin F. Posada C. (efposadac@unal.edu.co)
!!	  -# Cambia de indices de momento angular a capas, amplía hasta momento angular l.
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
!!
module ContractedGaussian_
  use Math_
  use CONTROL_
  use Exception_
  use String_
  implicit none
  
  type :: contractedGaussian
     integer :: id
     integer :: length
     integer :: angularMoment
     integer :: numCartesianOrbital
     integer :: owner
     real(8) :: origin(3)
     real(8) , allocatable :: orbitalExponents(:)
     real(8) , allocatable :: contractionCoefficients(:)
     real(8) , allocatable :: contNormalization(:)
     real(8) , allocatable :: primNormalization(:,:)
  end type contractedGaussian
  

  public :: &
       ContractedGaussian_saveToFile, &
       ContractedGaussian_loadFromFile, &
       ContractedGaussian_showInCompactForm, &
       ContractedGaussian_showInSimpleForm, 	&
       ContractedGaussian_normalizePrimitive, &
       ContractedGaussian_normalizeContraction, &
       ContractedGaussian_overlapIntegral, &
       ContractedGaussian_getShellCode, &
       ContractedGaussian_getAllAngularMomentIndex
  
  private :: &
       ContractedGaussian_obaraSaikaRecursion

contains

  !>
  !! @brief Saves the ContractedGaussian structure to file.
  !! @param this contracted gaussian
  !! @author E. F. Posada, 2013
  subroutine ContractedGaussian_saveToFile(this, unit)
    implicit none
    
    type(ContractedGaussian) :: this
    integer :: unit

    write(unit,*) this%id
    write(unit,*) this%length
    write(unit,*) this%angularMoment
    write(unit,*) this%numCartesianOrbital
    write(unit,*) this%owner
    write(unit,*) this%origin
    write(unit,*) this%orbitalExponents
    write(unit,*) this%contractionCoefficients
    write(unit,*) this%contNormalization
    write(unit,*) this%primNormalization

  end subroutine ContractedGaussian_saveToFile

  !>
  !! @brief Loads the ContractedGaussian structure from file.
  !! @param this contracted gaussian
  !! @author E. F. Posada, 2013
  subroutine ContractedGaussian_loadFromFile(this, unit)
    implicit none
    
    type(ContractedGaussian) :: this
    integer :: unit
    
    read(unit,*) this%id
    read(unit,*) this%length
    read(unit,*) this%angularMoment
    read(unit,*) this%numCartesianOrbital
    read(unit,*) this%owner
    read(unit,*) this%origin

    allocate(this%orbitalExponents(this%length))
    allocate(this%contractionCoefficients(this%length))
    
    read(unit,*) this%orbitalExponents
    read(unit,*) this%contractionCoefficients

    allocate(this%contNormalization(this%numCartesianOrbital))
    allocate(this%primNormalization(this%length, this%length*this%numCartesianOrbital))
                 
    read(unit,*) this%contNormalization
    read(unit,*) this%primNormalization

  end subroutine ContractedGaussian_loadFromFile


  !>
  !! Muestra en pantalla el valor de los atributos asociados a la gausiana
  !! contraida solicitada, En caso de indicarse una de sus gausianas primitivas,
  !! muestra informacion sobre dicha primitiva.
  subroutine ContractedGaussian_showInCompactForm( this )
    implicit none
    
    type(ContractedGaussian) , intent(in) :: this

    integer :: i, j
    character(9) :: shellCode(this%numCartesianOrbital)

    shellCode=ContractedGaussian_getShellCode(this)

    do i=1,this%length

       write (6,"(T10,I5,A9,A1,A6,F20.8,F20.8)") i,"        ", trim(shellCode(1)),"      ",&
            this%orbitalExponents(i) ,this%contractionCoefficients(i)
    end do
    
    !! TEST
    ! print*, "contraction normalization"
    ! do i=1, this%numCartesianOrbital
    !    print*, this%contNormalization(i)
    ! end do
    ! print*, "primitive normalization"
    ! do i = 1, this%length
    !    do j = 1, this%length*this%numCartesianOrbital
    !       print*, this%primNormalization(i, j)
    !    end do
    ! end do

  end subroutine ContractedGaussian_showInCompactForm

  subroutine ContractedGaussian_showInSimpleForm( this, unidOfOutput)
    implicit none
    type(ContractedGaussian) , intent(in) :: this
    integer :: unidOfOutput

    integer :: i
    character(9) :: shellCode(this%numCartesianOrbital)

    shellCode=ContractedGaussian_getShellCode(this)

    write (unidOfOutput,"(A1,I3,F5.2)") trim(String_getLowercase(shellCode(1)(1:1))),this%length,1.00

    do i=1,this%length
       write (unidOfOutput,"(ES19.10,ES19.10)") this%orbitalExponents(i) ,this%contractionCoefficients(i)
    end do

  end subroutine ContractedGaussian_showInSimpleForm


    
  !>
  !! @brief Normalization constant for a primitive gaussian
  !! @author S. A. Gonzalez
  !! @par History
  !!      -2011.02.04: E.F.Posada: extends for shells
  !! @version 1.1
  subroutine ContractedGaussian_normalizePrimitive( this )
    implicit none
    type(contractedGaussian) , intent( inout ) :: this
    
    integer, allocatable :: angularMomentIndex(:,:)
    integer :: i, j, m
    
    if(allocated(angularMomentIndex)) deallocate(angularMomentIndex)
    allocate(angularMomentIndex(3, this%numCartesianOrbital))
    
    call ContractedGaussian_getAllAngularMomentIndex( angularMomentIndex, this)
    
    
    ! write(*,*) "Normalization"
    do i = 1, this%length
       m = 1
       do j = 1, this%numCartesianOrbital
          this%primNormalization(i,m) =( ( 2.0_8*this%orbitalExponents(i)/Math_PI )**0.75_8 ) &
               / sqrt( &
               Math_factorial( 2_8 * angularMomentIndex(1, j) - 1_8,2 )&
               * Math_factorial( 2_8 * angularMomentIndex(2, j) - 1_8,2 )&
               * Math_factorial( 2_8 * angularMomentIndex(3, j) - 1_8,2 )/&
               ((4.0_8*this%orbitalExponents(i))**this%angularMoment))
          m = m + 1
       end do
          ! write(*,*) this%primNormalization(i,:)
    end do
    
  end subroutine ContractedGaussian_normalizePrimitive

  
  !>
  !! @brief Normalization constant for a contraction (shell)
  !! @author S. A. Gonzalez
  !! @par History
  !!      -2011.02.04: E.F.Posada: extends for shells
  !! @version 1.1  
  subroutine ContractedGaussian_normalizeContraction( this )
    implicit none

    type(ContractedGaussian) , intent(inout) :: this
    real(8) :: integralValue(this%numCartesianOrbital * this%numCartesianOrbital)
    integer :: i, j, m, n

    this%contNormalization = 1.0_8
    
    call ContractedGaussian_overlapIntegral( this , this, integralValue )

    !print*,  "integralValue", integralValue
    
    m = 0
    do i=1, this%numCartesianOrbital
       do j = 1, this%numCartesianOrbital
          m = m + 1
          if (j == i) then
             this%contNormalization(i) = 1.0_8 / sqrt( integralValue(m))
          end if
       end do
    end do
    
  end subroutine ContractedGaussian_normalizeContraction
  
  !>
  !! @brief Calculates overlap integral between two contractions (shell)
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2011.02.04: E.F.Posada: chage for usage on opints
  !! @return  output: overlap integral of a shell (all combinations)
  !! @version 1.0
  subroutine ContractedGaussian_overlapIntegral(contractedGaussianA, contractedGaussianB, integral)
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
    integer, allocatable :: angularMomentIndexA(:,:)
    integer, allocatable :: angularMomentIndexB(:,:)
    integer ::  i, m, p, q
    
    real(8), allocatable ::  x(:,:), y(:,:), z(:,:)
    real(8) :: AB2
    real(8) :: auxExponentA, auxCoefficientA, auxConstantA
    real(8) :: auxExponentB, auxCoefficientB, auxConstantB
    real(8) :: gamma, gammaInv
    real(8) :: PA(0:3), PB(0:3), P0(0:3)
    real(8) :: commonPreFactor
    real(8) :: x0, y0, z0

    integer :: angularMomentA, angularMomentB
    integer :: maxAngularMoment
    integer :: p1, p2 !< iteradores

    real(8) :: integralValue
    
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
          
          !! Start calculating integrals
          angularMomentA = sum(am1)
          angularMomentB = sum(am2)
          
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
          
          do p1=0, nprim1 - 1
             auxExponentA = exp1(p1)
             auxCoefficientA = coef1(p1)
             auxConstantA = nor1(p1)
             do p2=0, nprim2 - 1
                auxExponentB = exp2(p2)
                auxCoefficientB = coef2(p2)
                auxConstantB = nor2(p2)
                gamma = auxExponentA + auxExponentB
                gammaInv = 1.0/gamma
                
                P0(0) = (auxExponentA*A(0) + auxExponentB*B(0))*gammaInv
                P0(1) = (auxExponentA*A(1) + auxExponentB*B(1))*gammaInv
                P0(2) = (auxExponentA*A(2) + auxExponentB*B(2))*gammaInv
                PA(0) = P0(0) - A(0)
                PA(1) = P0(1) - A(1)
                PA(2) = P0(2) - A(2)
                PB(0) = P0(0) - B(0)
                PB(1) = P0(1) - B(1)
                PB(2) = P0(2) - B(2)
                
                commonPreFactor = exp(-auxExponentA*auxExponentB*AB2*gammaInv) &
                     * sqrt(Math_PI*gammaInv) * Math_PI * gammaInv &
                     * auxCoefficientA * auxCoefficientB * auxConstantA * auxConstantB
                
                !! recursion
                call ContractedGaussian_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMomentA+2, angularMomentB+2)
                
                x0 = x(am1(0),am2(0))
                y0 = y(am1(1),am2(1))
                z0 = z(am1(2),am2(2))
                
                !! Calculating integrals for primitives
                integralValue = integralValue + commonPreFactor*x0*y0*z0

             end do
          end do
          
          deallocate(x)
          deallocate(y)
          deallocate(z)

          !! Integral for shell
          integralValue = integralValue * contractedGaussianA%contNormalization(p) &
               * contractedGaussianB%contNormalization(q)

          integral(m) = integralValue

       end do
    end do

  end subroutine ContractedGaussian_overlapIntegral

  !>
  !!@brief Implementation of recursion proposed by Obara-Saika for overlap integrals.
  !!@author Edwin Posada, 2010
  !!@return x, y, z : recursion matrix
  !!@param PA, PB : reduced origin for gaussian A and B
  !!@param gamma : reduced exponent
  !!@see Gaussian product: if you want to know what reduced exponent and origin is.
  subroutine ContractedGaussian_obaraSaikaRecursion(x, y, z, PA, PB, gamma, angularMoment1, angularMoment2)
    implicit none

    real(8), intent(inout), allocatable :: x(:,:), y(:,:), z(:,:)
    real(8), intent(in) :: PA(0:3), PB(0:3)
    real(8), intent(in) :: gamma
    integer, intent(in) :: angularMoment1, angularMoment2

    real(8) :: pp
    integer :: i, j, k

    pp = 1/(2*gamma)

    x(0,0) = 1.0_8
    y(0,0) = 1.0_8
    z(0,0) = 1.0_8

    !! Upward recursion in j for i=0
    x(0,1) = PB(0)
    y(0,1) = PB(1)
    z(0,1) = PB(2)

    do j=1, angularMoment2 -1
       x(0,j+1) = PB(0)*x(0,j)
       y(0,j+1) = PB(1)*y(0,j)
       z(0,j+1) = PB(2)*z(0,j)
       x(0,j+1) = x(0,j+1) + j*pp*x(0,j-1)
       y(0,j+1) = y(0,j+1) + j*pp*y(0,j-1)
       z(0,j+1) = z(0,j+1) + j*pp*z(0,j-1)
    end do

    !! Upward recursion in i for all j
    x(1,0) = PA(0)
    y(1,0) = PA(1)
    z(1,0) = PA(2)

    do j=1, angularMoment2
       x(1,j) = PA(0)*x(0,j)
       y(1,j) = PA(1)*y(0,j)
       z(1,j) = PA(2)*z(0,j)
       x(1,j) = x(1,j) + j*pp*x(0,j-1)
       y(1,j) = y(1,j) + j*pp*y(0,j-1)
       z(1,j) = z(1,j) + j*pp*z(0,j-1)
    end do

    do i=1, angularMoment1 - 1
       x(i+1,0) = PA(0)*x(i,0)
       y(i+1,0) = PA(1)*y(i,0)
       z(i+1,0) = PA(2)*z(i,0)
       x(i+1,0) = x(i+1,0) + i*pp*x(i-1,0)
       y(i+1,0) = y(i+1,0) + i*pp*y(i-1,0)
       z(i+1,0) = z(i+1,0) + i*pp*z(i-1,0)
       do j=1, angularMoment2
          x(i+1,j) = PA(0)*x(i,j)
          y(i+1,j) = PA(1)*y(i,j)
          z(i+1,j) = PA(2)*z(i,j)
          x(i+1,j) = x(i+1,j) + i*pp*x(i-1,j)
          y(i+1,j) = y(i+1,j) + i*pp*y(i-1,j)
          z(i+1,j) = z(i+1,j) + i*pp*z(i-1,j)
          x(i+1,j) = x(i+1,j) + j*pp*x(i,j-1)
          y(i+1,j) = y(i+1,j) + j*pp*y(i,j-1)
          z(i+1,j) = z(i+1,j) + j*pp*z(i,j-1)
       end do
    end do

  end subroutine ContractedGaussian_obaraSaikaRecursion

  !>
  !! @brief 	Retorna el codigo de capa y la direccion espacial
  !!		de la contraccion especificada
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2011.02.04: E.F.Posada: creation.
  !! @version 1.0
  function ContractedGaussian_getShellCode( this ) result( output )
    implicit none
    type(ContractedGaussian) , intent(in) :: this
    character(9) :: output (this%numCartesianOrbital)

    type(Exception) :: ex

    character(1) :: indexCode(0:this%angularMoment)	!< Codigo para solo un indice de momento angular
    character(1) :: shellCode(0:8)			!< Codigo para una capa dada
    character(1) :: coordCode(3)			!< Codigo de las coordenadas cartesianas
    integer :: nx, ny, nz	 			!< Indices de momento angular
    integer :: i, j, m, u, v				!< Iteradores

    if ( this%angularMoment <= 8 ) then
       
       shellCode(0:8) = ["S", "P", "D", "F", "G", "H", "I", "J", "L"]
       coordCode(1:3) = ["x", "y", "z"]
       
       indexCode(0) = trim(shellCode(this%angularMoment))
       
       m = 0
       
       select case(CONTROL_instance%DIMENSIONALITY)
          
       case(3)	
          
          do i = 0 , this%angularMoment
             nx = this%angularMoment - i
             do j = 0 , i
                ny = i - j
                nz = j
                m = m + 1
                !! nx
                u = 0
                do v = 1, nx
                   u = u + 1
                   indexCode(u) = trim(coordCode(1))
                end do
                !! ny
                do v = 1, ny
                   u = u + 1
                   indexCode(u) = trim(coordCode(2))
                end do
                !! nz
                do v = 1, nz
                   u = u + 1
                   indexCode(u) = trim(coordCode(3))
                end do
                output(m) = trim(indexCode(1)(0:this%angularMoment))
             end do
          end do
          
       case(2)
          
          do i = 0 , this%angularMoment
             nx = this%angularMoment - i
             ny = i 
             nz = 0
             m = m + 1
             !! nx
             u = 0
             do v = 1, nx
                u = u + 1
                indexCode(u) = trim(coordCode(1))
             end do
             !! ny
             do v = 1, ny
                u = u + 1
                indexCode(u) = trim(coordCode(2))
             end do
             !! nz
             do v = 1, nz
                u = u + 1
                indexCode(u) = trim(coordCode(3))
             end do
             output(m) = trim(indexCode(1)(0:this%angularMoment))
          end do
          
       case(1)
          
          nx = this%angularMoment 
          ny = 0
          nz = 0
          
          m = m + 1
          !! nx
          u = 0
          do v = 1, nx
             u = u + 1
             indexCode(u) = trim(coordCode(1))
          end do
          !! ny
          do v = 1, ny
             u = u + 1
             indexCode(u) = trim(coordCode(2))
          end do
          !! nz
          do v = 1, nz
             u = u + 1
             indexCode(u) = trim(coordCode(3))
          end do
          output(m) = trim(indexCode(1)(0:this%angularMoment))
          
       case default
          
          call ContractedGaussian_exception( ERROR, "Class object ContratedGaussian in the getShellCode function",&
               "This Dimensionality is not avaliable") 
          
       end select
       
    else
       
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object ContractedGaussian in the getShellCode function" )
       call Exception_setDescription( ex, "This angular moment  isn't implemented" )
       call Exception_show( ex )
       
    end if
    
  end function ContractedGaussian_getShellCode

    !>                                                                      
  !! @brief return all angularMomentIndex for a contraction               
  !! @author E. F. Posada, efposadac@unal.edu.co
  !! @par History
  !!      -2011.02.04: E.F.Posada: creation.
  !! @version 1.0
  subroutine ContractedGaussian_getAllAngularMomentIndex(output, this, angularMoment)
    implicit none
    type(contractedGaussian), optional :: this
    integer, optional :: angularMoment
    integer, allocatable, intent(inout) :: output(:,:)
    
    integer :: totalAngularMoment
    integer :: counter
    integer :: x, y, z
    integer :: i, j
    
    if(present(this)) totalAngularMoment = this%angularMoment
    if(present(angularMoment)) totalAngularMoment = angularMoment
    
    counter = 1
    
    do i = 0 , totalAngularMoment
       x = totalAngularMoment - i
       do j = 0 , i
          y = i - j
          z = j          
          output(1:3, counter) = [x, y, z]
          counter = counter + 1
       end do
    end do
    
  end subroutine  ContractedGaussian_getAllAngularMomentIndex

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine ContractedGaussian_exception( typeMessage, description, debugDescription)
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
    
  end subroutine ContractedGaussian_exception
  
end module ContractedGaussian_
