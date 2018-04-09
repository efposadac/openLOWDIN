!!*****************************************************************************
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
  use LibxcInterface_
  implicit none

  type, public :: Functional
     character(30) :: species1
     character(30) :: species2
     character(50) :: name
     character(50) :: exchangeName
     character(50) :: correlationName
     integer :: shell !closedShell or openShell using libxc values
     real(8) :: exactExchangeFraction
     TYPE(xc_f03_func_t) :: xc1 
     TYPE(xc_f03_func_info_t) :: info1
     TYPE(xc_f03_func_t) :: xc2
     TYPE(xc_f03_func_info_t) :: info2
  end type Functional

  type(Functional), public, allocatable :: Functionals(:)

  public :: &
       Functional_createFunctionals, &
       Functional_constructor, &
       Functional_show, &
       Functional_getIndex, &
       Functional_getExchangeFraction, &
       Functional_libxcEvaluate, &
       Functional_LDAEvaluate, &
       Functional_CSEvaluate, &
       Functional_myCSEvaluate, &
       Functional_PSNEvaluate, &
       Functional_lowLimitEvaluate, &
       padevwn, &
       dpadevwn, &
       ecvwn, &
       vcvwn, &
       exdirac, &
       vxdirac, &
       srs, &
       pi, &
       eap_homogeneus,&
       dr_eap_homogeneus

  private
  real(8) vwna1,vwnb1,vwnc1,vwnx01
  real(8) vwna2,vwnb2,vwnc2,vwnx02
  real(8) vwna3,vwnb3,vwnc3,vwnx03
  parameter (vwna1=0.0621814,   vwna2=0.0310907,   vwna3=-0.0337737,  &
             vwnb1=3.7274400,   vwnb2=7.0604200,   vwnb3=1.1310710,   &
             vwnc1=12.9352000,  vwnc2=18.0578000,  vwnc3=13.0045000,  &
             vwnx01=-0.1049800, vwnx02=-0.3250000, vwnx03=-0.0047584)

   
contains

  subroutine Functional_createFunctionals()
    implicit none

    integer :: numberOfSpecies
    integer :: speciesID, otherSpeciesID, i
    character(50) :: labels(2), dftFile
    integer :: dftUnit
    numberOfSpecies = MolecularSystem_getNumberOfQuantumSpecies()

    allocate(Functionals(numberOfSpecies+numberOfSpecies*(numberOfSpecies-1)/2))

    i=1
    do speciesID=1, numberOfSpecies
       call Functional_constructor(Functionals(i), speciesID, speciesID)
       i=i+1
    end do

    do speciesID=1, numberOfSpecies-1
       do otherSpeciesID=speciesID+1, numberOfSpecies  
          call Functional_constructor(Functionals(i), speciesID, otherSpeciesID)
          i=i+1
       end do
    end do

  end subroutine Functional_createFunctionals

  
  subroutine Functional_constructor(this, speciesID, otherSpeciesID)
    implicit none
    type(Functional) :: this 
    integer :: speciesID
    integer :: otherSpeciesID
    
    this%name="NONE"
    this%correlationName="NONE"
    this%exchangeName="NONE"
    this%species1=MolecularSystem_getNameOfSpecie(speciesID)
    this%species2=MolecularSystem_getNameOfSpecie(otherSpeciesID)
    this%exactExchangeFraction=1.0_8

    if((this%species1 == "E-" .and. this%species2 == "E-") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-ALPHA") .or. &
         (this%species1 == "E-ALPHA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-BETA") .or. &
         (this%species1 == "E-BETA" .and. this%species2 == "E-ALPHA") ) then

       if(this%species1 == "E-") then
          this%shell=XC_UNPOLARIZED !closed shell
       else
          this%shell=XC_POLARIZED !open shell
       end if
       
       this%exactExchangeFraction=0.0_8
       
       if( CONTROL_instance%CALL_LIBXC)  then

          select case(trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL))         

          case("NONE")

             if (CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL .ne. "NONE") then
                
                this%name="exchange:"//CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL
                this%exchangeName="XC_"//CONTROL_instance%ELECTRON_EXCHANGE_FUNCTIONAL
                
                call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
                this%info1 = xc_f03_func_get_info(this%xc1)
                          
                if( xc_f03_func_info_get_family(this%info1) .eq. XC_FAMILY_HYB_GGA ) this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)
             
             end if
             
             if (CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL .ne. "NONE") then
                
                this%name=this%name//"-correlation:"//CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL
                this%correlationName="XC_"//CONTROL_instance%ELECTRON_CORRELATION_FUNCTIONAL
             
                call xc_f03_func_init(this%xc2, xc_f03_functional_get_number( this%correlationName), this%shell)
                this%info2 = xc_f03_func_get_info(this%xc2)
             
             end if
                          
          case("LDA")
             this%name="exchange:Slater-correlation:VWN5"
             this%exchangeName="XC_LDA_X"
             this%correlationName="XC_LDA_C_VwN"

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             call xc_f03_func_init(this%xc2, xc_f03_functional_get_number( this%correlationName), this%shell)
             this%info2 = xc_f03_func_get_info(this%xc2)

             
          case("PBE")
             this%name="exchange:PBE-correlation:PBE"
             this%exchangeName="XC_GGA_X_PBE"
             this%correlationName="XC_GGA_C_PBE"

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             call xc_f03_func_init(this%xc2, xc_f03_functional_get_number( this%correlationName),this%shell)
             this%info2 = xc_f03_func_get_info(this%xc2)

          case("BLYP")
             this%name="exchange:B88-correlation:LYP"
             this%exchangeName="XC_GGA_X_B88"
             this%correlationName="XC_GGA_C_LYP"

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             call xc_f03_func_init(this%xc2, xc_f03_functional_get_number( this%correlationName), this%shell)
             this%info2 = xc_f03_func_get_info(this%xc2)

          case("PBE0")
             this%name="exchange-correlation:PBE0"
             this%exchangeName="XC_HYB_GGA_XC_PBEH"
             this%correlationName="NONE" !!Both are combined
             
             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

          case("B3LYP")
             this%name="exchange-correlation:B3LYP"
             this%exchangeName="XC_HYB_GGA_XC_B3LYP5"
             this%correlationName="NONE" !!Both are combined

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

          case default

             this%name="exchange-correlation:"//trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)
             this%exchangeName=trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)
             this%correlationName=trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

             if( xc_f03_func_info_get_family(this%info1) .eq. XC_FAMILY_HYB_GGA) this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

             ! stop "ERROR: Please select an electronic functional to call LIBXC (LDA, PBE, BLYP, PBE0, B3LYP, or choose one from their web page)"
             
          end select
          
          ! print *, "sere yo maestro", this%species1, this%exactExchangeFraction
          
       else

          select case(trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL))         
             
          case("NONE")
             
             this%name="NONE"
             this%exactExchangeFraction=1.0_8

          case("LDA")

             this%name="exchange:Slater-correlation:VWN5"
             this%exchangeName="Slater"
             this%correlationName="VwN"
             this%exactExchangeFraction=0.0_8

          case default

             stop "ERROR: Please use LDA or call LIBXC to evaluate the selected electronic exchange correlation functional "

          end select

       end if

    !Provisional, only correlation between electron and other species (nuclei)
    elseif( (this%species1 .eq. "E-" .or.  this%species1 .eq. "E-ALPHA" .or.  this%species1 .eq. "E-BETA") .and. &
            (this%species2 .ne. "E-" .and. this%species2 .ne. "E-ALPHA" .and. this%species2 .ne. "E-BETA") )   then
       
       this%name="correlation:"//trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
       this%correlationName=trim(CONTROL_instance%NUCLEAR_ELECTRON_CORRELATION_FUNCTIONAL)
       this%exactExchangeFraction=1.0_8  !should be irrelevant

    else

       this%name="NONE"
       this%exactExchangeFraction=1.0_8

       
    end if

  end subroutine Functional_constructor

  subroutine Functional_show()
    implicit none
    type(Functional) :: this 
    integer :: i

    if( CONTROL_instance%CALL_LIBXC) then
       print *, "--------------------------------------------------------------------------------------"
       print *, "LIBXC library, Fermann, Miguel A. L. Marques, Micael J. T. Oliveira, and Tobias Burnus"
       print *, "Comput. Phys. Commun. 183, 2272 (2012) OAI: arXiv:1203.1739"
       print *, "LOWDIN-LIBXC Implementation V. 1.0  Moncada F. ; Reyes A. 2017"
    end if
    print *, ""
    print *, "--------------------------------------------------------------------------------------"
    print *, "|-------------------------------Functionals summary ---------------------------------|"
    print *, "--------------------------------------------------------------------------------------"
    print *, ""

    do i=1, size(Functionals)
       this=Functionals(i)

       if ((this%species1 == "E-" .and. this%species2 == "E-") .or. &
            (this%species1 == "E-ALPHA" .and. this%species2 == "E-ALPHA") .or. &
            (this%species1 == "E-ALPHA" .and. this%species2 == "E-BETA") .or. &
            (this%species1 == "E-BETA" .and. this%species2 == "E-BETA") .or. &
            (this%species1 == "E-BETA" .and. this%species2 == "E-ALPHA") ) then

          if( CONTROL_instance%CALL_LIBXC)  then

             if( this%correlationName .ne. "NONE" ) then

                write(*, "(T5,A10,A10,A5,A12,A)") trim(this%species1), trim(this%species2), "","exchange:", xc_f03_func_info_get_name(this%info1)
                print *, "family", xc_f03_func_info_get_family(this%info1), "shell", this%shell

                write(*, "(T5,A10,A10,A5,A12,A)") trim(this%species1), trim(this%species2), "","correlation:", xc_f03_func_info_get_name(this%info2)
                print *, "family", xc_f03_func_info_get_family(this%info2), "shell", this%shell

             else

                write(*, "(T5,A10,A10,A5,A21,A)") trim(this%species1), trim(this%species2), "", "exchange-correlation:", xc_f03_func_info_get_name(this%info1)

                ! print *, "family", xc_f03_func_info_get_family(this%info1), "shell", this%shell

             end if

          else
             write(*, "(T5,A10,A10,A5,A)") trim(this%species1), trim(this%species2), "",this%name

          end if
       else 
          write(*, "(T5,A10,A10,A5,A)") trim(this%species1), trim(this%species2), "",this%name

       end if

    end do

    print *, ""
    print *, "--------------------------------------------------------------------------------------"

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
       if (Functionals(i)%species1 == nameOfSpecies .and. Functionals(i)%species2 == otherNameOfSpecies) then
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
    ! Felix Moncada, Sep. 2017
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
    
    if ( this%shell .eq. XC_UNPOLARIZED) then
       allocate( e_exchange(n), e_correlation(n), v_exchange(n), v_correlation(n), vs_exchange(n), vs_correlation(n) )
    else
       allocate( e_exchange(n), e_correlation(n), v_exchange(2*n), v_correlation(2*n), vs_exchange(3*n), vs_correlation(3*n) )
    end if
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

    if( this%correlationName .ne. "NONE" ) then

       select case ( xc_f03_func_info_get_family(this%info2))
       case(XC_FAMILY_LDA)
          call xc_f03_lda_exc_vxc(this%xc2, n, rho, e_correlation, v_correlation)
       case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc_vxc(this%xc2, n, rho, sigma, e_correlation, v_correlation, vs_correlation )
       end select
    
    end if
    

    if ( this%shell .eq. XC_UNPOLARIZED) then
       exc(1:n)=e_exchange(1:n)+e_correlation(1:n)
       vxc(1:n)=v_exchange(1:n)+v_correlation(1:n)
       vsigma(1:n)=vs_exchange(1:n)+vs_correlation(1:n)
    else
       exc(1:n)=e_exchange(1:n)+e_correlation(1:n)
       vxc(1:2*n)=v_exchange(1:2*n)+v_correlation(1:2*n)
       vsigma(1:3*n)=vs_exchange(1:3*n)+vs_correlation(1:3*n)
    end if

    
    ! print *, "rho, energy density, energy density*rho, potential"
    ! do i = 1, n
    !    write(*,"(F10.6,1X,3F9.6)") rho(i), exc(i), exc(i)*rho(i), vxc(i)
    ! end do

    deallocate( e_exchange, e_correlation, v_exchange, v_correlation, vs_exchange, vs_correlation )
   
  end subroutine Functional_libxcEvaluate

  subroutine Functional_CSEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Hamess-Schiffer's Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2017
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a,b,c, q
    real(8), allocatable :: denominator(:)
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:epc17-2" ) then
       a=2.35
       b=2.4
       c=6.6
    else if(this%name .eq. "correlation:epc17-1" ) then
       a=2.35
       b=2.4
       c=3.2
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if
       
    allocate(denominator(n))

    denominator(1:n)=a-b*sqrt(rhoE(1:n)*rhoN(1:n))+c*rhoE(1:n)*rhoN(1:n)
    ! denominator(1:n)=a

    
!!!Energy
    ec(1:n)= -rhoE(1:n)*rhoN(1:n)/denominator(1:n)
    
!!!Potential  
    ! vcE(1:n)=vcE(1:n) -rhoN(1:n)/denominator(1:n)! + (c*rhoE(1:n)*rhoN(1:n)**2 - b/2*sqrt(rhoE(1:n))*rhoN(1:n)**(3/2))/denominator(1:n)**2
    ! vcN(1:n)=vcN(1:n) -rhoE(1:n)/denominator(1:n)! + (c*rhoN(1:n)*rhoE(1:n)**2 - b/2*sqrt(rhoN(1:n))*rhoE(1:n)**(3/2))/denominator(1:n)**2

    vcE(1:n)=vcE(1:n) + (b*sqrt(rhoE(1:n))*rhoN(1:n)**(3/2)-2*a*rhoN(1:n))/denominator(1:n)**2/2
    vcN(1:n)=vcN(1:n) + (b*sqrt(rhoN(1:n))*rhoE(1:n)**(3/2)-2*a*rhoE(1:n))/denominator(1:n)**2/2
    
    ! do i = 1, n
    !    if(denominator(i) .lt. 1E-6 .or. rhoE(i) .lt. 1E-6 .or. rhoN(i) .lt. 1E-6) then
    !       ec(i)=0.0
    !       vcE(i)=0.0
    !       vcN(i)=0.0
    !    end if
    ! end do


    ! print *, "i, rhoE, rhoN, denominator, energy density, potentialE, potentialN"
    ! do i = 1, n
    !    write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i),  ec(i), vcE(i), vcN(i)
    ! end do

    deallocate(denominator)
    
  end subroutine Functional_CSEvaluate

  subroutine Functional_myCSEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Hamess-Schiffer's Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2017
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a0,a1,a2,a3,a4,a5,a6,q
    real(8) :: beta, dbetaDE, dbetaDN
    real(8) :: denominator, ddenominatorDbeta
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:CS-myfit" ) then
       q=2.5
       a0=0.0011335821964296108
       a1=-0.027186701045457354
       a2=0.2559301045397849
       a3=-0.4955389585548209
       a4=0.3247814704559331
       a5=-0.05045899525215592
       a6=0.002906302490376316
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if
        
    do i = 1, n
       
       ! beta=q*(rhoE(i)*rhoN(i))**(1.0/6.0)
       ! dbetaDE=(1.0/6.0)*q*rhoN(i)**(1.0/6.0)/rhoE(i)**(5.0/6.0)
       ! dbetaDN=(1.0/6.0)*q*rhoE(i)**(1.0/6.0)/rhoN(i)**(5.0/6.0)

       beta=q*(rhoE(i)**(1.0/3.0)+rhoN(i)**(1.0/3.0))
       dbetaDE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
       dbetaDN=(1.0/3.0)*q*rhoN(i)**(-2.0/3.0)
       
       denominator=a0+a1*beta+a2*beta**2+a3*beta**3+a4*beta**4+a5*beta**5+a6*beta**6
       ddenominatorDbeta=a1+2*a2*beta+3*a3*beta**2+4*a4*beta**3+5*a5*beta**4+6*a6*beta**5
       
!!!Energy
       ec(i)= -rhoE(i)*rhoN(i)/denominator

       
!!!Potential
       
       vcE(i)=vcE(i) !- rhoN(i)/denominator*(1-rhoE(i)*dbetaDE*ddenominatorDbeta/denominator)
       vcN(i)=vcN(i) !- rhoE(i)/denominator*(1-rhoN(i)*dbetaDN*ddenominatorDbeta/denominator)

       ! if(denominator .lt. 1E-6) then
       !    ec(i)=0.0
       ! end if

       if(rhoE(i) .lt. 1E-6) then
          vcE(i)=0.0
       end if
       
       ! if(rhoE(i) .lt. denominator) then
       !    vcN(i)=vcN(i) - rhoE(i)/denominator*(1-rhoN(i)*dbetaDN*ddenominatorDbeta/denominator)
       ! else
       !    vcN(i)=vcN(i) - denominator/rhoE(i)*(1-rhoN(i)*dbetaDN*ddenominatorDbeta/denominator)
       ! end if

       if(rhoN(i) .lt. 1E-6) then
          vcN(i)=0.0
       end if

       
       ! write(*,"(I0.1,7F16.6)") i, rhoN(i), rhoE(i), rhoE(i)/denominator, -(1-rhoN(i)*dbetaDN*ddenominatorDbeta/denominator) ,vcN(i)  

    end do

  end subroutine Functional_myCSEvaluate

  
  subroutine Functional_PSNEvaluate( this, mass, n, rhoE, rhoP, ec, vcE, vcP )
    ! Evaluates Hamess-Schiffer's Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2017
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoP(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcP(*) !! Potentials - output   

    real(8) :: Aa,Ba,Ca,Bb,Cb,Cc,rse,rsp,drho_rse,drho_rsp
    real(8), allocatable :: denominator(:)
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:psn" ) then
       !*2 to convert from Rydbergs to a.u.
        Aa=69.7029*2.0_8
        Ba=-107.4927*2.0_8
        Bb=141.8458*2.0_8
        Ca=23.7182*2.0_8
        Cb=-33.6472*2.0_8
        Cc=5.21152*2.0_8
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if
           
    allocate(denominator(n))

    ! print *, "i, rhoE, rhoN, denominator, energy density, potentialE, potentialN"
    do i = 1, n

       rse= (3.0/(4.0*Math_PI*rhoE(i)))**(1.0/3.0)
       rsp= (3.0/(4.0*Math_PI*rhoP(i)))**(1.0/3.0)

       drho_rse= -0.206783/rhoE(i)**(4.0/3.0)
       drho_rsp= -0.206783/rhoP(i)**(4.0/3.0)
       !(1.0/3.0)*(3.0/(4.0*Math_PI))**(1.0/3.0)

       !Energy
       
       !This should be a parameter in CONTROL
       if( rse .ge. 8.0 .and. rsp .ge. 8.0) then
          denominator(i)= 4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          ec(i)=ec(i) + 1/denominator(i)

          
       else
          denominator(i)= Aa+Ba*(rse+rsp)+Ca*(rse**2+rsp**2) &
               +Bb*rse*rsp+Cb*(rse**2 *rsp + rse*rsp**2) + Cc*rse**2 *rsp**2 &
               +4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)
          ec(i)=ec(i) + 1/denominator(i)

       end if

       !Potential
       
       if( rse .ge. 8.0 .and. rsp .ge. 8.0) then

          vcE(i)=vcE(i)
          vcP(i)=vcP(i)
          ! denominator(i)= 4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)
          
          ! if( rse .ge. 50) then            
          !    vcE(i)=vcE(i) - drho_rse/denominator(i)**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2)
          ! else
          !    vcE(i)=vcE(i) - drho_rse/denominator(i)**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2+4*Math_PI*rse**2/eap_homogeneus(rsp))
          ! end if

          ! if( rsp .ge. 50) then            
          !   vcP(i)=vcP(i) - drho_rsp/denominator(i)**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2)
          ! else
          !    vcP(i)=vcP(i) - drho_rsp/denominator(i)**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2+4*Math_PI*rsp**2/eap_homogeneus(rse))
          ! end if
          
       else

          denominator(i)= Aa+Ba*(rse+rsp)+Ca*(rse**2+rsp**2) &
               +Bb*rse*rsp+Cb*(rse**2 *rsp + rse*rsp**2) + Cc*rse**2 *rsp**2 &
               +4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          vcE(i)=vcE(i) - drho_rse/denominator(i)**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2 &
               +4*Math_PI*rse**2/eap_homogeneus(rsp)&
               +Ba+2*Ca*rse+Bb*rsp+Cb*(2*rse*rsp + rsp**2)+2*Cc*rse*rsp**2)
          vcP(i)=vcP(i) - drho_rsp/denominator(i)**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2 &
               +4*Math_PI*rsp**2/eap_homogeneus(rse)&
               +Ba+2*Ca*rsp+Bb*rse+Cb*(rse**2 + 2*rse*rsp)+2*Cc*rse**2*rsp)
          
       end if
       
     !print *, i, rse, rsp , ec(i), vcE(i), vcP(i) !, vcE(i), - drho_rse/denominator(i)**2*(4*Math_PI*rse**2/eap_homogeneus(rsp)), rse**2, eap_homogeneus(rsp)
     !drho_rse, denominator(i)**2, dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2, 4*Math_PI*rse**2/eap_homogeneus(rsp)
       
       ! if(denominator(i) .ge. 0.0) then
       !    print *, "no jodas, en serio?"
       !    print *, i, rse, rsp, eap_homogeneus(rse), eap_homogeneus(rsp), denominator(i)
       ! end if
                 
    end do
    
    deallocate(denominator)
    
  end subroutine Functional_PSNEvaluate


  subroutine Functional_lowLimitEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates E/sqrt(PePn)
    ! Felix Moncada, 2018
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: energyDensity, b
    integer :: i

    print *, this%name
    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:lowlimit" ) then
       energyDensity=-0.5_8*mass/(mass+1.0_8)
       b=-0.5_8
    else
       ! STOP  "The nuclear electron functional chosen is not implemented"
    end if
        
    do i = 1, n

       ec(i)=rhoE(i)*rhoN(i)*energyDensity/(sqrt(rhoE(i)*rhoN(i)) + b*rhoE(i)*rhoN(i)**(3.0/2.0))

       ! vcE(i)=0.5_8*energyDensity*sqrt(rhoN(i))/sqrt(rhoE(i))
       ! vcN(i)=0.5_8*energyDensity*sqrt(rhoE(i))/sqrt(rhoN(i))
       
       ! write(*,"(I0.1,5F16.6)") i, rhoN(i), rhoE(i), ec(i), vcE(i), vcN(i)

    end do

  end subroutine Functional_lowLimitEvaluate


  
  
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


  real(8) function eap_homogeneus(rs)
    implicit none
    real(8) rs

    !!This expressions are in rydbergs
    if(rs .lt. 0.302) eap_homogeneus=-1.56/sqrt(rs)+(0.051*log(rs)-0.081)*log(rs)+1.14
    if(rs .ge. 0.302 .and. rs .lt. 0.56) eap_homogeneus=-0.92305-0.05459/rs**2.0
    if(rs .ge. 0.25 .and. rs .lt. 8.0) eap_homogeneus=-13.15111/(rs+2.5)**2.0 + 2.8655/(rs+2.5) - 0.6298
    if(rs .ge. 8.0) eap_homogeneus=-10250.57860/rs**6.0  + 44.50466/rs**3.0-0.524

    !!changing to a.u.
    eap_homogeneus=eap_homogeneus/2.0_8
  end function eap_homogeneus

  real(8) function dr_eap_homogeneus(rs)
    implicit none
    real(8) rs

    !!This expressions are in rydbergs
    if(rs .lt. 0.302) dr_eap_homogeneus=0.78/rs**(3.0/2.0)+(0.102*log(rs)-0.081)/rs
    if(rs .ge. 0.302 .and. rs .lt. 0.56) dr_eap_homogeneus=0.10918/rs**3.0
    if(rs .ge. 0.25 .and. rs .lt. 8.0) dr_eap_homogeneus=26.30222/(rs+2.5)**3.0 - 2.8655/(rs+2.5)**2.0
    if(rs .ge. 8.0) dr_eap_homogeneus=61503.4716/rs**7.0  -133.51399/rs**4.0

    !!changing to a.u.
    dr_eap_homogeneus=dr_eap_homogeneus/2.0_8
    
  end function dr_eap_homogeneus

  
end module Functional_

