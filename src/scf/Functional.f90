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
!! @brief This module manages the orbital and density represented in the DFT grids. Partially based on R. Flores-Moreno Parakata's modules
!! @author F. Moncada, 2017
module Functional_
  use Matrix_
  use Exception_
  use String_
  use MolecularSystem_
  use Grid_
  implicit none

  type, public :: Functional
     character(50) :: type
     character(30) :: species1
     character(30) :: species2

  end type Functional

  type(Functional), public, allocatable :: Functionals(:)

  public :: &
       Functional_constructor, &
       Functional_show, &
       padevwn, &
       dpadevwn, &
       ecvwn, &
       vcvwn, &
       exdirac, &
       vxdirac, &
       srs, &
       pi
  ! CalculateWaveFunction_getDensityAt, &
  !   	CalculateWaveFunction_getOrbitalValueAt
  !		CalculateWaveFunction_getFukuiFunctionAt

  private
  real(8) vwna1,vwnb1,vwnc1,vwnx01
  real(8) vwna2,vwnb2,vwnc2,vwnx02
  real(8) vwna3,vwnb3,vwnc3,vwnx03
  parameter (vwna1=0.0621814,   vwna2=0.0310907,   vwna3=-0.0337737,  &
             vwnb1=3.7274400,   vwnb2=7.0604200,   vwnb3=1.1310710,   &
             vwnc1=12.9352000,  vwnc2=18.0578000,  vwnc3=13.0045000,  &
             vwnx01=-0.1049800, vwnx02=-0.3250000, vwnx03=-0.0047584)

   
contains

  subroutine Functional_constructor(this, type, speciesID, otherSpeciesID)
    implicit none
    type(Functional) :: this 
    character(*) :: type
    integer :: speciesID
    integer :: otherSpeciesID

    this%type=type
    this%species1=MolecularSystem_getNameOfSpecie(speciesID)
    this%species2=MolecularSystem_getNameOfSpecie(otherSpeciesID)

    if( (this%species1 == "E-" .or. this%species1 == "E-ALPHA" .or. this%species1 == "E-BETA") .and. this%type=="exchange" ) then
       this%type=this%type//CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL
    else if( (this%species1 == "E-" .or. this%species1 == "E-ALPHA" .or. this%species1 == "E-BETA") .and. this%type=="correlation" ) then
       this%type=this%type//CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL
    else if( (this%species1 == "E-" .or. this%species1 == "E-ALPHA" .or. this%species1 == "E-BETA") .and. this%species1 .ne. this%species2)  then
       this%type=this%type//CONTROL_instance%ELECTRON_NUCLEAR_CORRELATION_FUNCTIONAL
    else
       this%type="NONE"
    end if

  end subroutine Functional_constructor

  subroutine Functional_show(this)
    implicit none
    type(Functional) :: this 

    print *, this%species1, this%species2, this%type

  end subroutine Functional_show

  subroutine exdirac(rho,ex,n)
    ! Evaluates Dirac exchange functional
    ! Roberto Flores-Moreno, May 2009
    implicit none
    integer n
    real(8) rho(*),ex(*)

    real(8) factor

    factor = 4.0/3.0
    ex(1:n) = (2.0*rho(1:n))**factor
    factor = -(0.75*(3.0/pi())**(1.0/3.0))
    ex(1:n) = factor*ex(1:n)

  end subroutine exdirac

  subroutine vxdirac(rho,vx,n)
    ! Evaluates Dirac exchange potential
    ! Roberto Flores-Moreno, May 2009
    implicit none
    integer n
    real(8) rho(*),vx(*)

    vx(1:n) = -((3.0/pi())**(1.0/3.0))*(2.0*rho(1:n))**(1.0/3.0)
  end subroutine vxdirac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! VWN Correlation functional
  !
  ! S. Vosko, L. Wilk, M. Nusair,
  ! Can. J. Phys. 58, 1200 (1980)

  real(8) function padevwn(x,a,b,c,x0)
    implicit none
    real(8) a,b,c,x,x0

    real(8) den,q,x0s,xs

    q = sqrt(4.0*c-b*b)
    xs = x*x
    x0s = x0*x0
    den = (x - x0)
    den = den*den
    padevwn = a*(log(xs/den) - ((x0s + c)*log((xs+b*x+c)/den) +         &
         2.0*b*(x0s-c)*atan(q/(2.0*x+b))/q)/(x0s+b*x0+c))
  end function padevwn

  real(8) function dpadevwn(x,a,b,c,x0)
    implicit none
    real(8) a,b,c,x,x0

    real(8) xs
    ! DPADE(x) = -(x/6) p'(x)
    xs = x*x
    dpadevwn = (1.0/3.0)*a*((1.0+b/(x-x0))*(xs/(xs+b*x+c))-1.0)
  end function dpadevwn

  subroutine ecvwn(rhoa,rhob,ec,n)
    ! Evaluates VWN correlation energy functional
    ! Roberto Flores-Moreno, Jun 2010
    implicit none
    integer n
    real(8) rhoa(*),rhob(*),ec(*)

    integer i
    real(8) rhot,s,ep,ef,ez,sp,vwnf2z,vwns,zs,zs4

    vwns = 2.0*(2.0**(1.0/3.0) - 1.0)
    vwnf2z = 4.0/(9.0*(2.0**(1.0/3.0) - 1.0))

    do i=1,n
       rhot = rhoa(i) + rhob(i)
       if (rhot .lt. 1E-12) cycle      
       s = srs(rhot)
       ep = padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
       zs = (rhoa(i) - rhob(i))/rhot
       ec(i) = ep
       if (zs.ne.0.0) then
          sp = ((1.0+zs)**(4.0/3.0)+(1.0-zs)**(4.0/3.0) - 2.0)/vwns
          ef = padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
          ez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)
          zs4 = zs**4
          ec(i)= ec(i) + sp*((ef-ep)*zs4 + ez*(1.0-zs4)/vwnf2z)
       end if
       ec(i) = 0.5*ec(i)*rhot
    end do
  end subroutine ecvwn


  subroutine vcvwn(rhoa,rhob,vca,vcb,n)
    ! Evaluates VWN correlation potential
    ! Roberto Flores-Moreno, Jun 2010
    implicit none
    integer n
    real(8) rhoa(*),rhob(*),vca(*),vcb(*)

    integer i
    real(8) def,defpz,dep,dez,dsp,ef,efpz,ep,ez,rhot
    real(8) s,sp,vcfp,vcpol,vwnf2z,vwns,zs,zs3

    vwns = 2.0*(2.0**(1.0/3.0) - 1.0)
    vwnf2z = 4.0/(9.0*(2.0**(1.0/3.0) - 1.0))

    do i=1,n
       rhot = rhoa(i) + rhob(i)
       if (rhot .lt. 1E-12) cycle
       s = srs(rhot)
       ep = padevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
       dep = dpadevwn(s,vwna1,vwnb1,vwnc1,vwnx01)
       vca(i) = ep + dep
       vcb(i) = vca(i)
       zs = (rhoa(i) - rhob(i))/rhot
       if (zs.ne.0.0) then
          sp = ((1.0+zs)**(4.0/3.0)+(1.0-zs)**(4.0/3.0) - 2.0)/vwns
          dsp = (4.0/3.0)*((1.0+zs)**(1.0/3.0)-(1.0-zs)**(1.0/3.0))/vwns
          zs3 = zs**3
          ef = padevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
          ez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
          efpz = zs3*(ef - ep - ez)
          def = dpadevwn(s,vwna2,vwnb2,vwnc2,vwnx02)
          dez = padevwn(s,vwna3,vwnb3,vwnc3,vwnx03)/vwnf2z
          defpz = zs3*(def - dep - dez)
          vcfp = sp*(zs*(efpz + defpz) + ez + dez)
          vcpol = dsp*(zs*efpz + ez) + 4.0*sp*efpz
          vca(i) = vca(i) + vcfp + (1.0 - zs)*vcpol
          vcb(i) = vcb(i) + vcfp - (1.0 + zs)*vcpol
       end if
       vca(i) = 0.5*vca(i)
       vcb(i) = 0.5*vcb(i)
    end do
  end subroutine vcvwn

  real(8) function srs(rhot)
    implicit none
    real(8) rhot
    srs = sqrt((0.75/(pi()*rhot))**(1.0/3.0))
  end function srs

  real(8) function pi()
    implicit none
    pi = 2.0*acos(0.0)
  end function pi

  
  
end module Functional_

