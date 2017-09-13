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
!! @brief This module manages the orbital and density represented in the DFT grids. Includes calls to libxc library. Partially based on R. Flores-Moreno Parakata's modules
!! @author F. Moncada, 2017
module Functional_
  use Matrix_
  use Exception_
  use String_
  use MolecularSystem_
  use Grid_ 
  use LibxcInterface_
  implicit none

  type, public :: Functional
     character(30) :: species1
     character(30) :: species2
     character(50) :: name
     real(8) :: exactExchangeFraction
     TYPE(xc_f03_func_t) :: xc1 
     TYPE(xc_f03_func_info_t) :: info1
     TYPE(xc_f03_func_t) :: xc2
     TYPE(xc_f03_func_info_t) :: info2
  end type Functional

  type(Functional), public, allocatable :: Functionals(:)

  public :: &
       Functional_constructor, &
       Functional_show, &
       Functional_getIndex, &
       Functional_getExchangeFraction, &
       Functional_libxcEvaluate, &
       Functional_LDAEvaluate, &
       Functional_CSEvaluate, &
       padevwn, &
       dpadevwn, &
       ecvwn, &
       vcvwn, &
       exdirac, &
       vxdirac, &
       srs, &
       pi

  private
  real(8) vwna1,vwnb1,vwnc1,vwnx01
  real(8) vwna2,vwnb2,vwnc2,vwnx02
  real(8) vwna3,vwnb3,vwnc3,vwnx03
  parameter (vwna1=0.0621814,   vwna2=0.0310907,   vwna3=-0.0337737,  &
             vwnb1=3.7274400,   vwnb2=7.0604200,   vwnb3=1.1310710,   &
             vwnc1=12.9352000,  vwnc2=18.0578000,  vwnc3=13.0045000,  &
             vwnx01=-0.1049800, vwnx02=-0.3250000, vwnx03=-0.0047584)

   
contains

  subroutine Functional_constructor(this, speciesID, otherSpeciesID)
    implicit none
    type(Functional) :: this 
    integer :: speciesID
    integer :: otherSpeciesID

    this%name="NONE"
    this%exactExchangeFraction=1.0_8
    this%species1=MolecularSystem_getNameOfSpecie(speciesID)
    this%species2=MolecularSystem_getNameOfSpecie(otherSpeciesID)

    if( CONTROL_instance%CALL_LIBXC)  then

       if(this%species1 == "E-" .and. this%species2 == "E-")  then
          if( CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL .ne. "NONE" ) then

             call xc_f03_func_init(this%xc1,&
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)), XC_UNPOLARIZED)
             this%info1 = xc_f03_func_get_info(this%xc1)

          else

             call xc_f03_func_init(this%xc1, &
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL)), XC_UNPOLARIZED)
             this%info1 = xc_f03_func_get_info(this%xc1)

             call xc_f03_func_init(this%xc2, &
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL)), XC_UNPOLARIZED)

             this%info2 = xc_f03_func_get_info(this%xc2)

          end if

          this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

          
       elseif ( (this%species1 == "E-ALPHA" .and. this%species2 == "E-ALPHA") .or. &
            (this%species1 == "E-ALPHA" .and. this%species2 == "E-BETA") .or. &
            (this%species1 == "E-BETA" .and. this%species2 == "E-BETA") .or. &
            (this%species1 == "E-BETA" .and. this%species2 == "E-ALPHA") ) then

          if( CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL .ne. "NONE" ) then

             call xc_f03_func_init(this%xc1,&
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)), XC_POLARIZED)
             this%info1 = xc_f03_func_get_info(this%xc1)

          else

             call xc_f03_func_init(this%xc1, &
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL)), XC_POLARIZED)
             this%info1 = xc_f03_func_get_info(this%xc1)

             call xc_f03_func_init(this%xc2, &
                  xc_f03_functional_get_number( "XC_"//trim(CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL)), XC_POLARIZED)

             this%info2 = xc_f03_func_get_info(this%xc2)

          end if

          this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)
          
       end if
       
    else

       if((this%species1 == "E-" .and. this%species2 == "E-") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-ALPHA") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-ALPHA") ) then

          this%name="x:"//trim(CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL)//"-c:"//trim(CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL)

          this%exactExchangeFraction=0.0_8 !We only have LDA outside libxc

       end if
       
    end if
       
    !Provisional, only correlation between electron and other species (nuclei)
    if( (this%species1 == "E-" .or. this%species1 == "E-ALPHA" .or. this%species1 == "E-BETA") .and. &
         (this%species2 .ne.  "E-" .and. this%species2 .ne. "E-ALPHA" .and. this%species2 .ne. "E-BETA") )   then
       
       this%name="c:"//trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
       
    end if
    
  end subroutine Functional_constructor

  subroutine Functional_show(this)
    implicit none
    type(Functional) :: this 


    if ((this%species1 == "E-" .and. this%species2 == "E-") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-ALPHA") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-ALPHA") ) then

       if( CONTROL_instance%CALL_LIBXC)  then

          if( CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL .ne. "NONE" ) then

             print *, trim(this%species1), trim(this%species2), "xc: ", xc_f03_func_info_get_name(this%info1)

          else
             print *, trim(this%species1), trim(this%species2), "x: ", xc_f03_func_info_get_name(this%info1)
             print *, trim(this%species1), trim(this%species2), "c: ", xc_f03_func_info_get_name(this%info2)
          end if

       else
             print *, trim(this%species1), trim(this%species2), this%name

       end if
    else 
       print *, trim(this%species1), trim(this%species2), this%name
       
    end if
    
  end subroutine Functional_show

    function Functional_getIndex(speciesID, otherSpeciesID) result( output)
    implicit none
    integer :: speciesID
    integer, optional :: otherSpeciesID
    integer :: output
    character(50) :: nameOfSpecies, otherNameOfSpecies
    integer i

    nameOfSpecies = MolecularSystem_getNameOfSpecie( speciesID )
    if( present(otherSpeciesID) )  then
       otherNameOfSpecies = MolecularSystem_getNameOfSpecie( otherSpeciesID )
    else
       otherNameOfSpecies = nameOfSpecies
    end if
    
    do i=1, size(Functionals(:))
       if (Functionals(i)%species1 == nameOfSpecies .and. Functionals(i)%species1 == otherNameOfSpecies) then
          output=i
          return
       end if
    end do
    
  end function Functional_getIndex

  function Functional_getExchangeFraction(speciesID) result( output)
    implicit none
    integer :: speciesID
    real(8) :: output
    integer :: index

    index=Functional_getIndex(speciesID)

    output=Functionals(index)%exactExchangeFraction

  end function Functional_getExchangeFraction

  subroutine Functional_libxcEvaluate(this, n,rho,sigma, exc, vxc, vsigma)
    ! Call LIBXC to evaluate electronic exchange correlation functionals
    ! Roberto Flores-Moreno, May 2009
    implicit none
    Type(Functional) :: this
    integer :: n !!gridSize
    real(8) :: rho(*), sigma(*) !!Density and gradient - input
    real(8) :: exc(*), vxc(*) !! Energy density and potential - output   
    real(8) :: vsigma(*)  !! Energy derivative with respect to the gradient - output   
    real(8), allocatable :: e_exchange(:),e_correlation(:)
    real(8), allocatable :: v_exchange(:),v_correlation(:)
    real(8), allocatable :: vs_exchange(:),vs_correlation(:)
    integer :: i, vmajor, vminor, vmicro, func_id = 1
    TYPE(xc_f03_func_t) :: xc_func
    TYPE(xc_f03_func_info_t) :: xc_info

    allocate( e_exchange(n), e_correlation(n), v_exchange(n), v_correlation(n), vs_exchange(n), vs_correlation(n) )
    
    e_exchange=0.0_8
    e_correlation=0.0_8
    v_exchange=0.0_8
    v_correlation=0.0_8
    vs_exchange=0.0_8
    vs_correlation=0.0_8

    select case ( xc_f03_func_info_get_family(this%info1))
    case(XC_FAMILY_LDA)
       call xc_f03_lda_exc_vxc(this%xc1, n, rho, e_exchange, v_exchange)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
       call xc_f03_gga_exc_vxc(this%xc1, n, rho, sigma, e_exchange, v_exchange, vs_exchange)
    end select

    if ( CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL .eq. "NONE") then

       select case ( xc_f03_func_info_get_family(this%info2))
       case(XC_FAMILY_LDA)
          call xc_f03_lda_exc_vxc(this%xc2, n, rho, e_correlation, v_correlation)
       case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc(this%xc2, n, rho, sigma, e_correlation, v_correlation, vs_correlation )
       end select
    
    end if
    
    exc(1:n)=e_exchange(1:n)+e_correlation(1:n)
    vxc(1:n)=v_exchange(1:n)+v_correlation(1:n)
    vsigma(1:n)=vs_exchange(1:n)+vs_correlation(1:n)
    
    ! print *, "rho, energy density, energy density*rho, potential"
    ! do i = 1, n
    !    write(*,"(F10.6,1X,3F9.6)") rho(i), exc(i), exc(i)*rho(i), vxc(i)
    ! end do

    deallocate( e_exchange, e_correlation, v_exchange, v_correlation, vs_exchange, vs_correlation )
   
  end subroutine Functional_libxcEvaluate

  subroutine Functional_CSEvaluate(n,rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Hamess-Schiffer's Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2017
    implicit none
    integer :: n !!gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a,b,c
    real(8), allocatable :: denominator(:)
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    a=2.35
    b=2.4
    c=6.6

    allocate(denominator(n))
    
    denominator(1:n)=a-b*sqrt(rhoE(1:n)*rhoN(1:n))+c*rhoE(1:n)*rhoN(1:n)
    
!!!Energy
    ec(1:n)= -rhoE(1:n)*rhoN(1:n)/denominator(1:n)
    
!!!Potential  
    vcE(1:n)=vcE(1:n) -rhoN(1:n)/denominator(1:n) + (c*rhoE(1:n)*rhoN(1:n)**2 - b/2*sqrt(rhoE(1:n))*rhoN(1:n)**(3/2))/denominator(1:n)**2

    vcN(1:n)=vcN(1:n) -rhoE(1:n)/denominator(1:n) + (c*rhoN(1:n)*rhoE(1:n)**2 - b/2*sqrt(rhoN(1:n))*rhoE(1:n)**(3/2))/denominator(1:n)**2

    deallocate(denominator)

    ! print *, "rhoE, rhoN, energy density, potentialE, potentialN"
    ! do i = 1, n
    !    write(*,"(5F9.6)") rhoE(i), rhoN(i), ec(i), vcE(i), vcN(i)
    ! end do

    
  end subroutine Functional_CSEvaluate

  
  subroutine Functional_LDAEvaluate(n,rhoA, rhoB, exc, vxcA, vxcB )
    ! Evaluates Dirac exchange and VWN correlation functionals
    ! Roberto Flores-Moreno, May 2009
    implicit none
    integer :: n !!gridSize
    real(8) :: rhoA(*), rhoB(*) !!Alpha and beta Densities - input
    real(8) :: exc(*), vxcA(*) !! Energy density and potential - output   
    real(8) , optional :: vxcB(*) !! Energy density and potential - output   

    real(8), allocatable :: e_exchange(:),e_correlation(:)
    real(8), allocatable :: v_exchange(:),va_correlation(:),vb_correlation(:)
    
    allocate( e_exchange(n), e_correlation(n), v_exchange(n), va_correlation(n), vb_correlation(n) )

    e_exchange=0.0_8
    e_correlation=0.0_8
    v_exchange=0.0_8
    va_correlation=0.0_8
    vb_correlation=0.0_8

    
!!!Energy
    call exdirac( (rhoA(1:n) + rhoB(1:n) )/2 , e_exchange, n)         
    call ecvwn( rhoA, rhoB, e_correlation, n)
    exc(1:n)=e_exchange(1:n)+e_correlation(1:n)

!!!Potential  
    call vxdirac( (rhoA(1:n) + rhoB(1:n) )/2 , v_exchange, n)         
    call vcvwn( rhoA, rhoB , va_correlation, vb_correlation, n)
    vxcA(1:n)= v_exchange(1:n) + va_correlation(1:n)
    if (present(vxcB)) vxcB(1:n) = v_exchange(1:n) + vb_correlation(1:n)
       
    deallocate( e_exchange, e_correlation, v_exchange, va_correlation, vb_correlation )
    
  end subroutine Functional_LDAEvaluate

  
  
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

