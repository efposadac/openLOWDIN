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
     real(8) :: mass1
     real(8) :: mass2
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
       Functional_EPCEvaluate, &
       Functional_IKNEvaluate, &
       Functional_MLCSEvaluate, &
       Functional_MLCSAEvaluate, &
       Functional_MLCSANEvaluate, &
       Functional_myCSEvaluate, &
       Functional_expCSEvaluate, &
       Functional_expCSGGAEvaluate, &
       Functional_PSNEvaluate, &
       Functional_PSNAPEvaluate, &
       Functional_lowLimitEvaluate, &
       Functional_getBeta, &
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
    this%mass1=MolecularSystem_getMass(speciesID)
    this%mass2=MolecularSystem_getMass(otherSpeciesID)

    
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
                          
                !if( xc_f03_func_info_get_family(this%info1) .eq. XC_FAMILY_HYB_GGA )
                this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)
             
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
             this%exchangeName="XC_"//trim(CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL)
             this%correlationName="NONE"

             call xc_f03_func_init(this%xc1, xc_f03_functional_get_number( this%exchangeName), this%shell)
             this%info1 = xc_f03_func_get_info(this%xc1)

             this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

             !if( xc_f03_func_info_get_family(this%info1) .eq. XC_FAMILY_HYB_GGA)
             this%exactExchangeFraction=xc_f03_hyb_exx_coef(this%xc1)

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
    real(8) :: p,qe,qn,qen,q2en,q3en,Eab,Ea2b,Eab2,a0,q0,q2,q4
    
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
                ! print *, "family", xc_f03_func_info_get_family(this%info1), "shell", this%shell

                write(*, "(T5,A10,A10,A5,A12,A)") trim(this%species1), trim(this%species2), "","correlation:", xc_f03_func_info_get_name(this%info2)
                ! print *, "family", xc_f03_func_info_get_family(this%info2), "shell", this%shell

             else

                write(*, "(T5,A10,A10,A5,A21,A)") trim(this%species1), trim(this%species2), "", "exchange-correlation:", xc_f03_func_info_get_name(this%info1)

                ! print *, "family", xc_f03_func_info_get_family(this%info1), "shell", this%shell

             end if

          else
             write(*, "(T5,A10,A10,A5,A)") trim(this%species1), trim(this%species2), "",this%name

          end if
       else 
          write(*, "(T5,A10,A10,A5,A)") trim(this%species1), trim(this%species2), "",this%name

          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "rhoE3") print *, "Using as correlation length: beta=q*rhoE^(1/3)"
          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "rhoE6rhoN6") print *, "Using as correlation length: beta=q*rhoE^(1/6)*rhoN^(1/6)"
          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN") print *, "Using as correlation length: beta=1/(q*rhoE^(-1/3)+r*rhoN^(-1))"
          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN3As") print *, "Using as correlation length: beta=qe*rhoE^(1/3)+qn*rhoN^(1/3)"
          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN3rhoEN6") print *, "Using as correlation length: beta=q*rhoE^(1/3)+p*rhoN^(1/3)+r*rhoE^(1/6)*rhoN^(1/6)"
          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "newBeta") print *, "beta=(qe*rhoE(i)+qn*rhoN(i)+q2en*(rhoE(i)-rhoN(i))**2/(rhoE(i)+rhoN(i))+q3en*(rhoE(i)-rhoN(i))**3/(rhoE(i)+rhoN(i))**2)**(1.0/3.0)"

          if(this%name .ne. "NONE" .and. CONTROL_instance%BETA_FUNCTION .eq. "newnewBeta") print *, "beta=(q0*(rhoE(i)+rhoN(i))+q2*(rhoE(i)-rhoN(i))**2/(rhoE(i)+rhoN(i))+q3*(rhoE(i)-rhoN(i))**4/(rhoE(i)+rhoN(i))**3)**(1.0/3.0)"

          if(CONTROL_instance%BETA_FUNCTION .eq. "newBeta") then

             if(this%mass2 .gt. 2.0) then !hydrogen
                print *, "electron-hydrogen correlation parameters"
                a0=4.5839773752240566113
                Ea2b=0.527444
                Eab2=0.597139 
             else !positron
                print *, "electron-positron correlation parameters"
                a0=2.2919886876120283056
                Ea2b=0.262005
                Eab2=0.262005
             end if

             p=1.0

             if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) then
                Ea2b=CONTROL_instance%DUMMY_REAL(1)
                Eab2=CONTROL_instance%DUMMY_REAL(2)
                p=CONTROL_instance%DUMMY_REAL(3)
             end if
                          
             print *, "Ea2b=", Ea2b
             print *, "Eab2=", Eab2
             qe=a0*(11.0/8.0/Ea2b-3.0/4.0/Eab2)
             qn=a0*(11.0/8.0/Eab2-3.0/4.0/Ea2b)
             q2en=a0*3.0/16.0*(1.0/Ea2b+1.0/Eab2)
             q3en=a0*9.0/16.0*(1.0/Eab2-1.0/Ea2b)
             print *, "qe=", qe
             print *, "qn=", qn
             print *, "q2en=", q2en
             print *, "q3en=", q3en

             print *, "p", CONTROL_instance%DUMMY_REAL(3)

          else if(CONTROL_instance%BETA_FUNCTION .eq. "newnewBeta") then

             if(this%mass2 .gt. 2.0) then !hydrogen
                STOP "this beta function only works for electron-positron"
             else !positron
                a0=2.2919886876120283056
                Eab=0.25
                Eab2=0.2620050702329801
             end if


             if (this%name .eq. "correlation:expCS-GGA-noA") a0=4.5839773752240566113                 

             p=1.0

             if(CONTROL_instance%BETA_PARAMETER_A .ne. 0 .or. CONTROL_instance%BETA_PARAMETER_B .ne. 0 .or. CONTROL_instance%BETA_PARAMETER_C .ne. 0) then
                Eab=CONTROL_instance%BETA_PARAMETER_A
                Eab2=CONTROL_instance%BETA_PARAMETER_B
                p=CONTROL_instance%BETA_PARAMETER_C
             end if
                          
             print *, "Eab=", Eab
             print *, "Eab2=", Eab2
             q0=a0/2/Eab
             q2=a0*(-5/Eab+53/Eab2/8)
             q4=a0*(9/Eab/2-45/Eab2/8)
             print *, "q0=", q0
             print *, "q2=", q2
             print *, "q4=", q4

             print *, "p", CONTROL_instance%DUMMY_REAL(3)


          else 
             if(this%name .ne. "NONE" .and. (CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) ) then
                print *, "q", CONTROL_instance%DUMMY_REAL(1)
                print *, "p", CONTROL_instance%DUMMY_REAL(2)
                print *, "r", CONTROL_instance%DUMMY_REAL(3)
             end if

          end if

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
       call xc_f03_lda_exc_vxc(this%xc1, int(n,8), rho, e_exchange, v_exchange)
    case(XC_FAMILY_GGA)
       call xc_f03_gga_exc_vxc(this%xc1, int(n,8), rho, sigma, e_exchange, v_exchange, vs_exchange)
    end select

    if( this%correlationName .ne. "NONE" ) then

       select case ( xc_f03_func_info_get_family(this%info2))
       case(XC_FAMILY_LDA)
          call xc_f03_lda_exc_vxc(this%xc2, int(n,8), rho, e_correlation, v_correlation)
       case(XC_FAMILY_GGA)
          call xc_f03_gga_exc_vxc(this%xc2, int(n,8), rho, sigma, e_correlation, v_correlation, vs_correlation )
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

  subroutine Functional_EPCEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
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
    real(8) :: denominator, densityThreshold
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD !TODO: add to other functionals

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

    !$omp parallel private(denominator)
    !$omp do schedule (dynamic)
    do i = 1, n

       denominator=a-b*sqrt(rhoE(i)*rhoN(i))+c*rhoE(i)*rhoN(i)

!!!Energy density
       ! ec(i)= -rhoE(1:n)*rhoN(1:n)/denominator(1:n)
       ec(i)= -rhoN(i)/denominator

!!!Potential  

       if( rhoE(i)+rhoN(i) .gt. densityThreshold ) then !
          vcE(i)= (rhoE(i)*rhoN(i)*(c*rhoN(i)-b*rhoN(i)/(2*sqrt(rhoE(i)*rhoN(i))))-rhoN(i)*denominator)/denominator**2
          vcN(i)= (rhoN(i)*rhoE(i)*(c*rhoE(i)-b*rhoE(i)/(2*sqrt(rhoE(i)*rhoN(i))))-rhoE(i)*denominator)/denominator**2
       else
          vcE(i)=0.0
          vcN(i)=0.0
       end if
       
    end do
    !$omp end do 
    !$omp end parallel
    
    ! print *, "i, rhoE, rhoN, denominator, energy density, potentialE, potentialN"
    ! do i = 1, n
    !    write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i),  ec(i), vcE(i), vcN(i)
    ! end do

  end subroutine Functional_EPCEvaluate

  subroutine Functional_IKNEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates YUTAKA IMAMURA, HIROYOSHI KIRYU, HIROMI NAKAI Colle Salvetti nuclear electron correlation functional J Comput Chem 29: 735–740, 2008
    ! Only works for Hydrogen - can be extended
    ! Felix Moncada, 2019
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a1,a2,a3,a4,q
    real(8) :: beta, numerator, denominator, F, dbetaE, dbetaN, dFdbeta
    integer :: i, sign

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    a1=14.48724087490633
    a2=-19.119574638807872
    a3=1.2732395447351628
    a4=-2.256758334191025
       
    if(this%name .eq. "correlation:ikn-nsf" ) then
       if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .and. CONTROL_instance%DUMMY_REAL(2) .ne. 0) then
          q=CONTROL_instance%DUMMY_REAL(1)
       else
          q=4.971
       end if
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    do i = 1, n

       if(rhoE(i) .gt. 1E-6 .and. rhoN(i) .gt. 1E-6 ) then !
!!!Energy
          beta=q*rhoE(i)**(1.0/3.0)

          if(beta.gt. 0.1) then
             ! beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
             numerator=a1+a2*beta
             denominator=2*Math_PI*beta**4*exp(a3/beta**2+a4/beta)

             F=numerator/denominator

             ec(i)=rhoN(i)*F
             ! ec(i)=rhoE(i)*rhoN(i)*F

!!!Potential
             dFdbeta=a2/denominator-4*numerator/denominator/beta+(2*a3/beta**3+a4/beta**2)*numerator/denominator

             dbetaE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
             dbetaN=0.0

             vcE(i)=rhoN(i)*(F+rhoE(i)*dbetaE*dFdbeta)
             vcN(i)=rhoE(i)*(F+rhoN(i)*dbetaN*dFdbeta)

             ! write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i), ec(i), vcE(i), vcN(i)
             ! write(*,"(I0.1,5F16.6)") 1000, beta,  F, dGdbeta, dFdbeta, dbetaE

          end if
       end if
    end do

    ! STOP
    
  end subroutine Functional_IKNEvaluate
 
  subroutine Functional_MLCSEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Mejia-de la Lande Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2018
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a1,a2,r1,r2,p,q
    real(8) :: beta, attractiveExp, repulsiveExp, F, G, dbetaE, dbetaN, dFdbeta, dGdbeta
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:mlcs-fit" ) then

       a1=2.0943951023931953
       a2=2.849869774919022
       r1=0.7607437033951782
       r2=0.433750315318239
       if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .and. CONTROL_instance%DUMMY_REAL(2) .ne. 0) then
          q=CONTROL_instance%DUMMY_REAL(1)
          p=CONTROL_instance%DUMMY_REAL(2)
       else
          q=1.2
          p=0.5
       end if
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    do i = 1, n

       if(rhoE(i) .gt. 1E-6 .and. rhoN(i) .gt. 1E-6 ) then !
          !!!Energy
          beta=q*rhoE(i)**(1.0/3.0)
          ! beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
          attractiveExp=exp(-a2*(beta**2))
          repulsiveExp=exp(-r2/beta)
          F=a1*attractiveExp/beta**2-r1*repulsiveExp/beta**3
          ! ec(i)= -p*rhoE(i)*rhoN(i)*F
          ec(i)= -p*rhoN(i)*F

          !!!Potential
          dFdbeta= -2*a1*(a2*beta**2+1)*attractiveExp/beta**3 - r1*(r2-3*beta)*repulsiveExp/beta**5
          
          !dGdbeta=-2*a1*a2*beta*attractiveExp - r1/beta**2*repulsiveExp + r1*r2/beta**3*repulsiveExp
          !dFdbeta=p*(-2*G/beta**3 + dGdbeta/beta**2)
          
          ! dbetaE=(1.0/6.0)*q*rhoE(i)**(-5.0/6.0)*rhoN(i)**(1.0/6.0)
          ! dbetaN=(1.0/6.0)*q*rhoN(i)**(-5.0/6.0)*rhoE(i)**(1.0/6.0)
          dbetaE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
          dbetaN=0.0
          
          vcE(i)=-p*rhoN(i)*(F+rhoE(i)*dbetaE*dFdbeta)
          vcN(i)=-p*rhoE(i)*(F+rhoN(i)*dbetaN*dFdbeta)
          
          ! write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i), ec(i), vcE(i), vcN(i)
          ! write(*,"(I0.1,5F16.6)") 1000, beta,  F, dGdbeta, dFdbeta, dbetaE

       end if
    end do


    ! STOP
    
  end subroutine Functional_MLCSEvaluate

  subroutine Functional_MLCSAEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Mejia-de la Lande Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2018
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a1,a2,r1,r2,p,q
    real(8) :: beta, attractiveExp, repulsiveExp, F, G, dbetaE, dbetaN, dFdbeta, dGdbeta
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:mlcs-a" ) then

       a1=2.0943951023931953
       a2=2.849869774919022
       r1=0.0
       r2=0.433750315318239
       if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .and. CONTROL_instance%DUMMY_REAL(2) .ne. 0) then
          q=CONTROL_instance%DUMMY_REAL(1)
          p=CONTROL_instance%DUMMY_REAL(2)
       else
          q=1.2
          p=0.5
       end if
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    do i = 1, n

       if(rhoE(i) .gt. 1E-6 .and. rhoN(i) .gt. 1E-6 ) then !
          !!!Energy
          beta=q*rhoE(i)**(1.0/3.0)
          ! beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
          attractiveExp=exp(-a2*(beta**2))
          repulsiveExp=exp(-r2/beta)
          F=a1*attractiveExp/beta**2-r1*repulsiveExp/beta**3
          ! ec(i)= -p*rhoE(i)*rhoN(i)*F
          ec(i)= -p*rhoN(i)*F

          !!!Potential
          dFdbeta= -2*a1*(a2*beta**2+1)*attractiveExp/beta**3 - r1*(r2-3*beta)*repulsiveExp/beta**5
          
          !dGdbeta=-2*a1*a2*beta*attractiveExp - r1/beta**2*repulsiveExp + r1*r2/beta**3*repulsiveExp
          !dFdbeta=p*(-2*G/beta**3 + dGdbeta/beta**2)
          
          ! dbetaE=(1.0/6.0)*q*rhoE(i)**(-5.0/6.0)*rhoN(i)**(1.0/6.0)
          ! dbetaN=(1.0/6.0)*q*rhoN(i)**(-5.0/6.0)*rhoE(i)**(1.0/6.0)
          dbetaE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
          dbetaN=0.0
          
          vcE(i)=-p*rhoN(i)*(F+rhoE(i)*dbetaE*dFdbeta)
          vcN(i)=-p*rhoE(i)*(F+rhoN(i)*dbetaN*dFdbeta)
          
          ! write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i), ec(i), vcE(i), vcN(i)
          ! write(*,"(I0.1,5F16.6)") 1000, beta,  F, dGdbeta, dFdbeta, dbetaE

       end if
    end do


    ! STOP
    
  end subroutine Functional_MLCSAEvaluate

  subroutine Functional_MLCSANEvaluate( this, mass, n, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Mejia-de la Lande Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2018
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a1,a2,r1,r2,p,q
    real(8) :: beta, attractiveExp, repulsiveExp, F, G, dbetaE, dbetaN, dFdbeta, dGdbeta
    real(8) :: v_exchange(n),va_correlation(n),vb_correlation(n)
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:mlcs-an" ) then

       a1=2.0943951023931953
       a2=2.849869774919022
       r1=0.0
       r2=0.433750315318239
       if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .and. CONTROL_instance%DUMMY_REAL(2) .ne. 0) then
          q=CONTROL_instance%DUMMY_REAL(1)
          p=CONTROL_instance%DUMMY_REAL(2)
       else
          q=1.2
          p=0.5
       end if
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    do i = 1, n

       if(rhoE(i) .gt. 1E-6 .and. rhoN(i) .gt. 1E-6 ) then !
          !!!Energy
          ! beta=q*rhoE(i)**(1.0/3.0)
          beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
          attractiveExp=exp(-a2*(beta**2))
          repulsiveExp=exp(-r2/beta)
          F=a1*attractiveExp/beta**2-r1*repulsiveExp/beta**3
          ! ec(i)= -p*rhoE(i)*rhoN(i)*F
          ec(i)= -p*rhoN(i)*F

          !!!Potential
          dFdbeta= -2*a1*(a2*beta**2+1)*attractiveExp/beta**3 - r1*(r2-3*beta)*repulsiveExp/beta**5
          
          !dGdbeta=-2*a1*a2*beta*attractiveExp - r1/beta**2*repulsiveExp + r1*r2/beta**3*repulsiveExp
          !dFdbeta=p*(-2*G/beta**3 + dGdbeta/beta**2)
          
          dbetaE=(1.0/6.0)*q*rhoE(i)**(-5.0/6.0)*rhoN(i)**(1.0/6.0)
          dbetaN=(1.0/6.0)*q*rhoN(i)**(-5.0/6.0)*rhoE(i)**(1.0/6.0)
          ! dbetaE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
          ! dbetaN=0.0
          
          ! vcE(i)=-p*rhoN(i)*(F+rhoE(i)*dbetaE*dFdbeta)
          ! vcN(i)=-p*rhoE(i)*(F+rhoN(i)*dbetaN*dFdbeta)
          vcE(i)=-p*rhoN(i)*rhoE(i)*(F+dbetaE*dFdbeta)
          vcN(i)=-p*rhoE(i)*rhoN(i)*(F+dbetaN*dFdbeta)
          
          ! write(*,"(I0.1,5F16.6)") i, rhoE(i), rhoN(i), ec(i), vcE(i), vcN(i)
          ! write(*,"(I0.1,5F16.6)") 1000, beta,  F, dGdbeta, dFdbeta, dbetaE

       end if
    end do


    ! STOP
    
  end subroutine Functional_MLCSANEvaluate
  
  
  subroutine Functional_myCSEvaluate( this, mass, npoints, rhoE, rhoN, ec, vcE, vcN )
    ! Evaluates Hamess-Schiffer's Colle Salvetti nuclear electron correlation functional
    ! Felix Moncada, 2019
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: npoints !!nuclear gridSize
    real(8) :: rhoE(*), rhoN(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcN(*) !! Potentials - output   

    real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,az,an,b0,b1,b2,b3,b4,bz,bn,p,q,r
    real(8) :: ke,kn    
    real(8) :: beta, dbetaE, dbetaN, F, dFdbeta, aPoly,aExp,bPoly,bExp, daPolydbeta, daExpdbeta, dbPolydbeta, dbExpdbeta
    integer :: i,n

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:CS-myfit" ) then
       a0=20.799624419985    
       a1=-19.251078153683    
       a2=68.341133121032    
       a3=-349.312459728153  
       a4=621.759104204928   
       a5=0.0
       a6=0.0
       a7=0.0
       a8=0.0
       az=7.262349396307      
       an=9.101297095659     
       b0=3.042892080859     
       b1=1.606594302396     
       b2=38.431539423963    
       b3=-110.273758102273  
       b4=130.813870751403   
       bz=4.357404128407      
       bn=12.297852188465    
    else if(this%name .eq. "correlation:Imamura-myfit" ) then
       a0=-1.810905109364
       a1=-1.167135966502
       a2=12.330825453317
       a3=-143.958134459753
       a4=773.724912824595
       a5=-2246.847075516435
       a6=3681.632763571440
       a7=-3248.547769822670
       a8=1216.285109989069
       az=1.736460214413
       an=3.614606996422
       b0=3.042974813585
       b1=2.952784475691
       b2=30.539736874367
       b3=-95.876667515705
       b4=157.531741188944
       bz=0.195361337190
       bn=17.465905244385
       ! a0=-1.810905109364
       ! a1=-1.052850613233
       ! a2=10.008633190911
       ! a3=-128.659319481101
       ! a4=743.178299971580
       ! a5=-2257.693608998759
       ! a6=3731.656228439009
       ! a7=-3205.916666708522
       ! a8=1147.270968266223
       ! az=1.665373819617
       ! an=3.624419432661
       ! b0=3.042974813585
       ! b1=3.248066178957
       ! b2=26.755007188150
       ! b3=-82.805230932937
       ! b4=144.122496204253
       ! bz=0.306222273912
       ! bn=21.383069425991
    else if(this%name .eq. "correlation:Mejia-myfit" ) then
       a0=-2.094395102394
       a1=-1.025867935018
       a2=-0.538148600474
       a3=11.072663674859
       a4=-11.728894706789
       a5=0
       a6=0
       a7=0
       a8=0
       az=5.485962579901
       an=6.089408324396
       b0=1.521487406791
       b1=0.075813099062
       b2=-0.153966977767
       b3=-3.962056400173
       b4=2.453311545282
       bz=2.864883246299
       bn=9.451138253283
    else if(this%name .eq. "correlation:MejiaA-myfit" ) then
       a0=-2.094395102393
       a1=1.919828725183
       a2=9.666337157333
       a3=-15.283848155658
       a4=6.031393373580
       a5=0
       a6=0
       a7=0
       a8=0
       az=0.911004697711
       an=6.773374962808
       b0=0.760743703395
       b1=0.055672495857
       b2=-0.291021401973
       b3=-0.206323289017
       b4=0
       bz=1.807726544299
       bn=5.580898227664
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .and. CONTROL_instance%DUMMY_REAL(2) .ne. 0) then
       q=CONTROL_instance%DUMMY_REAL(1)
       p=CONTROL_instance%DUMMY_REAL(2)
       r=CONTROL_instance%DUMMY_REAL(3)
    else
       q=1.0
       p=1.0
       r=1.0
    end if
    
    !$omp parallel & 
    !$omp& private(i,beta, dbetaE, dbetaN, F, dFdbeta, aPoly,aExp,bPoly,bExp, daPolydbeta, daExpdbeta, dbPolydbeta, dbExpdbeta),&
    !$omp& shared(mass, n, rhoE, rhoN, ec, vcE, vcN)
    n = omp_get_thread_num() + 1
    !$omp do schedule (dynamic) 
    
    do i = 1, npoints
       if(CONTROL_instance%BETA_FUNCTION .eq. "rhoE3") then
          beta=q*rhoE(i)**(1.0/3.0)
          dbetaE=(1.0/3.0)*q*rhoE(i)**(-2.0/3.0)
          dbetaN=0.0
       else if(CONTROL_instance%BETA_FUNCTION .eq. "rhoE6rhoN6") then
          beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
          dbetaE=(1.0/6.0)*q*rhoE(i)**(-5.0/6.0)*rhoN(i)**(1.0/6.0)
          dbetaN=(1.0/6.0)*q*rhoN(i)**(-5.0/6.0)*rhoE(i)**(1.0/6.0)
       else if(CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN3As") then
          ke=3*(3.0/(4.0*Math_PI))**(4.0/3.0)
          kn=8.0/3.0*sqrt(2.0/3.0)
          beta=q/(1/(ke*rhoE(i)**(1.0/3.0))+1/(kn*rhoN(i)**(1.0/3.0)))
          dbetaE=q/(3*ke*rhoE(i)**(4.0/3.0)*(1/(ke*rhoE(i)**(1.0/3.0))+1/(kn*rhoN(i)**(1.0/3.0)))**2.0)
          dbetaN=q/(3*kn*rhoN(i)**(4.0/3.0)*(1/(ke*rhoE(i)**(1.0/3.0))+1/(kn*rhoN(i)**(1.0/3.0)))**2.0)
       else if(CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN") then
          beta=1/(q*rhoE(i)**(-1.0/3.0)+r*rhoN(i)**(-1.0))
          dbetaE=q/(3*rhoE(i)**(4.0/3.0)*(q*rhoE(i)**(-1.0/3.0)+r*rhoN(i)**(-1.0))**2.0)
          dbetaN=r/(rhoN(i)**(2.0)*(q*rhoE(i)**(-1.0/3.0)+r*rhoN(i)**(-1.0))**2.0)
       else
          beta=q*rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)
          dbetaE=(1.0/6.0)*q*rhoE(i)**(-5.0/6.0)*rhoN(i)**(1.0/6.0)
          dbetaN=(1.0/6.0)*q*rhoN(i)**(-5.0/6.0)*rhoE(i)**(1.0/6.0)
       end if

       if(beta .gt. 0.01 ) then !
!       if(rhoE(i) .gt. 1E-6 .and. rhoN(i) .gt. 1E-6 ) then !
          !!!Energy
          aPoly=a0+a1*beta+a2*beta**2+a3*beta**3+a4*beta**4+a5*beta**5+a6*beta**6+a7*beta**7+a8*beta**8
          bPoly=b0+b1/beta+b2/beta**2+b3/beta**3+b4/beta**4
          aExp=exp(-az*beta**an)
          bExp=(1-exp(-bz*beta**bn))
          
          F=aPoly*aExp/beta**2+bPoly*bExp/beta**3
          ! ec(i)= -p*rhoE(i)*rhoN(i)*F
          ec(i)= -p*rhoN(i)*F

          !!!Potential
          daPolydbeta=a1+2*a2*beta+3*a3*beta**2+4*a4*beta**3+5*a5*beta**4+6*a6*beta**5+7*a7*beta**6+8*a8*beta**7
          daExpdbeta=-az*an*beta**(an-1)*exp(-az*beta**an)
          dbPolydbeta=-b1/beta**2-2*b2/beta**3-3*b3/beta**4-4*b4/beta**5
          dbExpdbeta=bz*bn*beta**(bn-1)*exp(-bz*beta**bn)
          
          dFdbeta=-2*aPoly*aExp/beta**3+(aPoly*daExpdbeta+daPolydbeta*aExp)/beta**2-3*bPoly*bExp/beta**4+(bPoly*dbExpdbeta+dbPolydbeta*bExp)/beta**3
          
          vcE(i)=-p*rhoN(i)*(F+rhoE(i)*dbetaE*dFdbeta)
          if(rhoN(i) .gt. 0.01 .or. mass .lt. 10.0) then !This is a cut off to achieve convergence in the SCF cycle
             vcN(i)=-p*rhoE(i)*(F+rhoN(i)*dbetaN*dFdbeta)
          end if
          ! write(*,"(I0.1,6F16.6)") i, beta, dFdbeta, dbetaE, dbetaN, vcE(i), vcN(i)

       end if
    end do
    !$omp end do nowait
    !$omp end parallel

  end subroutine Functional_myCSEvaluate


    subroutine Functional_expCSEvaluate( this, mass, npoints, rhoE, rhoP, ec, vcE, vcP )
    ! Evaluates my Exponential Jastrow factor functional
    ! Felix Moncada, 2019
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: npoints !!nuclear gridSize
    real(8) :: rhoE(*), rhoP(*) !! electron and positive particle Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcP(*) !! Potentials - output   

    real(8) :: a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,p,qe,qn,qen,q2en,q3en,Ea2b,Eab2,Eab,q0,q2,q4
    real(8) :: beta, dbetaE, dbetaP, F, dFdbeta, aPoly,aExp,bPoly,bExp, daPolydbeta, daExpdbeta, dbPolydbeta, dbExpdbeta
    real(8) :: d2BdE2, d2BdP2, d2BdEP !! dummys - not required
    real(8) :: deltaQ
    real :: time1, time2
    integer :: i,n

    real(8) :: densityThreshold
    
    densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD !TODO: add to other functionals
    
    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:expCS-A" ) then
       if(mass .gt. 2.0) then !nuclear- adiabatic
          a4=1.0
          a3=330.096328569496899372
          a2=126.389946461302141674
          a1=37.8099831062836374483
          a0=4.5839773752240566113
          b4=0.6572515786440476772
          b3=216.7870644040942790
          b2=39.66058188679972412
          b1=10.22481338790345212
          b0=1.0
       else  !positron-adiabatic
          ! a4=1.0
          ! a3=1303.3072101399865635
          ! a2=249.92864475592566026
          ! a1=37.282895111835410721
          ! a0=2.2919886876120283056
          ! b4=1.3145031572880953546
          ! b3=1712.8945795597543538
          ! b2=156.85114166883440203
          ! b1=20.216312587194126413
          ! b0=1.0
          a4=1.0
          a3=1151.5935487650133
          a2=223.8321298117117
          a1=32.8929178106386
          a0=2.2919886876120283056
          b4=1.3145031572880953546
          b3=1513.5873924109194
          b2=141.52688227410547
          b1=18.267727642586607
          b0=1.0
       end if
    else if(this%name .eq. "correlation:expCS-noA" ) then
       if(mass .gt. 2.0) then !nuclear- no adiabatic
          a4=1.0
          a3=259.74684041327130047
          a2=118.92498372366549632
          a1=40.237074092688264982
          a0=9.1679547504481132225
          b4=0.3286257893220238386
          b3=85.2776888283076625
          b2=12.8341294979177963
          b1=5.45967946556013866
          b0=1.0
       else  !positron-no adiabatic
          a4=1.0
          a3=1043.1272348010363592
          a2=238.60672146759363273
          a1=40.579167476742475919
          a0=4.5839773752240566113
          b4=0.6572515786440476772
          b3=685.58251729720346024
          b2=51.522539301737416155
          b1=10.995665380524416865
          b0=1.0
       end if          
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if

    p=1.0
    if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 ) then
       p=CONTROL_instance%DUMMY_REAL(1) 
    end if
    
    
    time1=omp_get_wtime()
    !$omp parallel private(beta, dbetaE, dbetaP, d2BdE2, d2BdP2, d2BdEP, F, dFdbeta, aPoly,aExp,bPoly,bExp, daPolydbeta, daExpdbeta, dbPolydbeta, dbExpdbeta)
    !$omp do schedule (dynamic)
    do i = 1, npoints

       call Functional_getBeta( rhoE(i), rhoP(i), mass, a0, beta, dbetaE, dbetaP, d2BdE2, d2BdP2, d2BdEP)

       if( rhoE(i) .gt. densityThreshold .and. rhoP(i) .gt. densityThreshold ) then !
          !!!Energy
          aPoly=a0+a1*beta+a2*beta**2+a3*beta**3+a4*beta**4
          bPoly=b0*beta**3+b1*beta**4+b2*beta**5+b3*beta**6+b4*beta**7
          F=aPoly/bPoly

         !!!Potential
          daPolydbeta=a1+2*a2*beta+3*a3*beta**2+4*a4*beta**3
          dbPolydbeta=3.0*b0*beta**2+4.0*b1*beta**3+5.0*b2*beta**4+6.0*b3*beta**5+7.0*b4*beta**6
          dFdbeta=(daPolydbeta*bPoly-aPoly*dbPolydbeta)/bPoly**2

       else if( rhoE(i) .gt. densityThreshold .or. rhoP(i) .gt. densityThreshold) then !
          F=a0/b0*beta**(-3)
          dFdbeta=-3*a0/b0*beta**(-4)

       else
          F=0.0
          dFdbeta=0.0
          
       end if
       
       ! ec(i)= -p*rhoE(i)*rhoN(i)*F
       ec(i)= -p*rhoP(i)*F
       vcE(i)=-p*rhoP(i)*(F+rhoE(i)*dbetaE*dFdbeta)
       vcP(i)=-p*rhoE(i)*(F+rhoP(i)*dbetaP*dFdbeta)

       ! write(*,"(I0.1,9F20.10)") i, rhoE(i), rhoP(i), beta, dbetaE, dbetaP, F, ec(i)*rhoE(i), vcE(i), vcP(i) 
       ! write(*,"(I0.1,5F16.6)") i, rhoE(i)/rhoN(i), rhoE(i)**(1.0/6.0)*rhoN(i)**(1.0/6.0)/rhoN(i)**(1.0/3.0), beta/rhoN(i)**(1.0/3.0), F*beta**3.0, ec(i)/rhoN(i)
       ! write(*,"(I0.1,6E20.10)") i, beta, dFdbeta, dbetaE, dbetaN, vcE(i), vcN(i)
       ! write(*,"(I0.1,5F16.6)") i, aPoly, bPoly, daPolydbeta, dbPolydbeta, dFdbeta
    end do
    !$omp end do 
    !$omp end parallel

    time2=omp_get_wtime()
    ! write(*,"(A,F10.3,A4)") "**expCSEvaluate:", time2-time1 ," (s)"


  end subroutine Functional_expCSEvaluate


  subroutine Functional_expCSGGAEvaluate( this, mass, npoints, electronDensity, electronGradient, positronDensity, positronGradient, &
       ec, vcE, vcgE, vcP, vcgP )
    ! Evaluates my Exponential Jastrow factor functional with Gradient terms
    ! Felix Moncada, 2020
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: npoints !!nuclear gridSize
    type(Vector) :: electronDensity, positronDensity !! electron and nuclear Densities - input
    type(Vector) :: electronGradient(3), positronGradient(3) !! electron and nuclear gradient - input
    type(Vector) :: ec !! Energy density - output, per electron density 
    type(Vector) :: vcE, vcP !! Potentials - output   
    type(Vector) :: vcgE(3), vcgP(3) !! Gradient Potentials - output   

    real(8) :: densityThreshold
    
    real(8) :: a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,p,Ea2b,Eab2,Eab,q0,q2,q4,g1,g2
    real(8) :: rhoE,rhoP,rhoTot,rhoDif,sigmaEE,sigmaPP,sigmaEP
    real(8) :: beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP
    real(8) :: aPoly,bPoly, daPolydB, dbPolydB, d2aPolydB2, d2bPolydB2
    real(8) :: F, dFdB, d2FdB2, dFdE, dFdP
    real(8) :: G, dGdB, d2GdB2, dGdE, dGdP, d2GdE2, d2GdP2, d2GdEP
    real :: time1, time2
    integer :: i, dir

    densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD !TODO: add to other functionals

    if(this%name .eq. "correlation:expCS-GGA-noA" ) then
       !positron-no adiabatic
       a4=1.0
       a3=1043.1272348010363592
       a2=238.60672146759363273
       a1=40.579167476742475919
       a0=4.5839773752240566113
       b4=0.6572515786440476772
       b3=685.58251729720346024
       b2=51.522539301737416155
       b1=10.995665380524416865
       b0=1.0
    else !positron-adiabatic
       a4=1.0
       a3=1151.5935487650133
       a2=223.8321298117117
       a1=32.8929178106386
       a0=2.2919886876120283056
       b4=1.3145031572880953546
       b3=1513.5873924109194
       b2=141.52688227410547
       b1=18.267727642586607
       b0=1.0
    end if
    
    if(mass .gt. 2.0) STOP "the expCSGGA functional only works for positron-electron correlation at the moment"

    p=1.0
    g1=1.0
    g2=1.0
    
    if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) then
       g1=CONTROL_instance%DUMMY_REAL(1) !coefficient 
       g2=CONTROL_instance%DUMMY_REAL(2) !exponent
       p=CONTROL_instance%DUMMY_REAL(3)
    end if

    
    ! time1=omp_get_wtime()
    !$omp parallel private( rhoE,rhoP,rhoTot,rhoDif,sigmaEE,sigmaPP,sigmaEP, &
    !$omp& beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP, &
    !$omp& aPoly,bPoly, daPolydB, dbPolydB, d2aPolydB2, d2bPolydB2, &
    !$omp& F, dFdB, d2FdB2, dFdE, dFdP, &
    !$omp& G, dGdB, d2GdB2, dGdE, dGdP, d2GdE2, d2GdP2, d2GdEP, &
    !$omp& dir)
    !$omp do schedule (dynamic)
    do i = 1, npoints
       rhoE=electronDensity%values(i)
       rhoP=positronDensity%values(i)
       sigmaEE=electronGradient(1)%values(i)**2+electronGradient(2)%values(i)**2+electronGradient(3)%values(i)**2
       sigmaPP=positronGradient(1)%values(i)**2+positronGradient(2)%values(i)**2+positronGradient(3)%values(i)**2
       sigmaEP=electronGradient(1)%values(i)*positronGradient(1)%values(i)+&
            electronGradient(2)%values(i)*positronGradient(2)%values(i)+&
            electronGradient(3)%values(i)*positronGradient(3)%values(i)
       
       rhoTot=rhoE+rhoP
       rhoDif=rhoE-rhoP

       call Functional_getBeta( rhoE, rhoP, mass, a0, &
            beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP)
       
!!!LDA terms
       if( rhoE .gt. densityThreshold .and. rhoP .gt. densityThreshold ) then !
          aPoly=a0+a1*beta+a2*beta**2+a3*beta**3+a4*beta**4
          bPoly=b0*beta**3+b1*beta**4+b2*beta**5+b3*beta**6+b4*beta**7
          daPolydB=a1+2.0*a2*beta+3.0*a3*beta**2+4.0*a4*beta**3
          dbPolydB=3.0*b0*beta**2+4.0*b1*beta**3+5.0*b2*beta**4+6.0*b3*beta**5+7.0*b4*beta**6
          d2aPolydB2=2.0*a2+6.0*a3*beta+12.0*a4*beta**2
          d2bPolydB2=6.0*b0*beta+12.0*b1*beta**2+20.0*b2*beta**3+30.0*b3*beta**4+42.0*b4*beta**5

          F=aPoly/bPoly
          dFdB=(daPolydB*bPoly-aPoly*dbPolydB)/bPoly**2
          d2FdB2=( 2.0*aPoly*dbPolydB**2 &
               +bPoly**2*d2aPolydB2&
               -bPoly*(2*daPolydB*dbPolydB&
               +aPoly*d2bPolydB2)&
               )/bPoly**3
          
       else if( (rhoE .gt. densityThreshold) .or. (rhoP .gt. densityThreshold)  ) then 
          F=a0/b0*beta**(-3)
          dFdB=-3.0*a0/b0*beta**(-4)
          d2FdB2=12.0*a0/b0*beta**(-5)
       else
          F=0.0
          dFdB=0.0
          d2FdB2=0.0
       end if

!!!GGA Terms
       if( beta**3 .gt. densityThreshold ) then
          G=g1*exp(-g2/beta)*F/beta**2
          dGdB=g1*exp(-g2/beta)*(F*(g2-2.0*beta) + dFdB*beta**2)/beta**4
          d2GdB2=g1*exp(-g2/beta)*(F*(6.0*beta**2-6*beta*g2+g2**2) + beta**2*(dFdB*(2.0*g2-4.0*beta) + d2FdB2*beta**2))/beta**6
       else
          G=0.0
          dGdB=0.0
          d2GdB2=0.0
       end if


!!!! F and B derivatives with respect to the densities       
       dFdE=dBdE*dFdB
       dFdP=dBdP*dFdB
       dGdE=dBdE*dGdB
       dGdP=dBdP*dGdB
       d2GdE2=dGdB*d2BdE2+d2GdB2*dBdE**2
       d2GdP2=dGdB*d2BdP2+d2GdB2*dBdP**2
       d2GdEP=dGdB*d2BdEP+d2GdB2*dBdE*dBdP

       
       !The multiplication with the electronic density is performed again in Grid Manager
       ec%values(i)= -p*(rhoE*rhoP*F&
            -sigmaEP*G&
            -1.0/4.0*sigmaEE*rhoP*dGdE&
            -1.0/4.0*sigmaEP*rhoP*dGdP&
            -1.0/4.0*sigmaPP*rhoE*dGdP&
            -1.0/4.0*sigmaEP*rhoE*dGdE)/rhoE
       
       vcE%values(i)=-p*(rhoP*(F+rhoE*dFdE)&
            -1.0/4.0*rhoE*sigmaEP*d2GdE2&
            -1.0/4.0*rhoE*sigmaPP*d2GdEP&
            -1.0/4.0*rhoP*sigmaEE*d2GdE2&
            -1.0/4.0*rhoP*sigmaEP*d2GdEP&
            -5.0/4.0*sigmaEP*dGdE&
            -1.0/4.0*sigmaPP*dGdP)
       
       vcP%values(i)=-p*(rhoE*(F+rhoP*dFdP)&
            -1.0/4.0*rhoE*sigmaEP*d2GdEP&
            -1.0/4.0*rhoE*sigmaPP*d2GdP2&
            -1.0/4.0*rhoP*sigmaEE*d2GdEP&
            -1.0/4.0*rhoP*sigmaEP*d2GdP2&
            -5.0/4.0*sigmaEP*dGdP&
            -1.0/4.0*sigmaEE*dGdE)

       !Gradien potential is returned
       do dir=1,3
          vcgE(dir)%values(i)=-p*(-1.0/2.0*rhoP*dGdE*electronGradient(dir)%values(i)&
               -(G+1.0/4.0*rhoE*dGdE+1.0/4.0*rhoE*dGdP)*positronGradient(dir)%values(i))

          vcgP(dir)%values(i)=-p*(-1.0/2.0*rhoE*dGdP*positronGradient(dir)%values(i)&
               -(G+1.0/4.0*rhoP*dGdP+1.0/4.0*rhoE*dGdE)*electronGradient(dir)%values(i))
       end do
       
       ! write(*,"(I0.1,6E20.10)") i, vcgE(1)%values(i), vcgE(2)%values(i), vcgE(3)%values(i), vcgP(1)%values(i), vcgP(2)%values(i), vcgP(3)%values(i)
       ! write(*,"(I0.1,7E20.10)") i, rhoE, rhoP, beta, F, ec%values(i)*rhoE, vcE%values(i), vcP%values(i)
       ! write(*,"(3E20.10)") sigmaEE, sigmaPP, sigmaEP
       ! write(*,"(5E20.10)") dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP
       ! write(*,"(4E20.10)") dFdB, d2FdB2, dGdB, d2GdB2
       ! write(*,"(5E20.10)") dGdE, dGdP, d2GdE2, d2GdP2, d2GdEP

       
    !    ! write(*,"(I0.1,5F16.6)") i, rhoE(i)/rhoP(i), rhoE(i)**(1.0/6.0)*rhoP(i)**(1.0/6.0)/rhoP(i)**(1.0/3.0), beta/rhoP(i)**(1.0/3.0), F*beta**3.0, ec(i)/rhoP(i)
    !    ! write(*,"(I0.1,6E20.10)") i, beta, dFdB, dBdE, dBdP, vcE(i), vcN(i)
    !    ! write(*,"(I0.1,5F16.6)") i, aPoly, bPoly, daPolydB, dbPolydB, dFdB
    end do
    !$omp end do 
    !$omp end parallel

    ! time2=omp_get_wtime()
    ! ! write(*,"(A,F10.3,A4)") "**expCSEvaluate:", time2-time1 ," (s)"


  end subroutine Functional_expCSGGAEvaluate

  subroutine Functional_getBeta( rhoE, rhoP, positiveMass, functionalLimitConstant, &
       beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP)
    ! Evaluates beta function for electron-positive particle correlation energy from low density energy limits
    ! Felix Moncada, 2020
    implicit none
    real(8) :: rhoE, rhoP !! electron and positive particle Densities - input
    real(8) :: positiveMass !!positive particle mass
    real(8) :: functionalLimitConstant !!kf: Assuming that at low densities F=kf/beta^3
    real(8) :: beta, dBdE, dBdP, d2BdE2, d2BdP2, d2BdEP !! - outputs

    real(8) :: Eab, Eab2, Ea2b
    real(8) :: q0,q1,q2,q3,q4    
    real(8) :: rhoTot,rhoDif
    real(8) :: cutOff
    
    real(8) :: densityThreshold

    densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD !TODO: add to other functionals
    rhoTot=rhoE+rhoP
    rhoDif=rhoE-rhoP

    select case(CONTROL_instance%BETA_FUNCTION)

    case("newnewBeta")
       if(CONTROL_instance%BETA_PARAMETER_A .ne. 0.0 .and. CONTROL_instance%BETA_PARAMETER_B .ne. 0.0) then
          Eab=CONTROL_instance%BETA_PARAMETER_A
          Eab2=CONTROL_instance%BETA_PARAMETER_B
       else
          if( positiveMass .gt. 2.0) then
             Eab=0.5
             Eab2=0.53
          else
             Eab=0.25
             Eab2=0.2620050702329801             
          end if
       end if

       q0=functionalLimitConstant/2/Eab
       q2=functionalLimitConstant*(-5/Eab+53/Eab2/8)
       q4=functionalLimitConstant*(9/Eab/2-45/Eab2/8)

!!! Beta terms
       if( rhoTot .gt. densityThreshold ) then !

          beta=(q0*rhoTot&
               +q2*rhoDif**2/rhoTot&
               +q4*rhoDif**4/rhoTot**3&
               )**(1.0/3.0)

          dBdE=(q0&
               -q2*rhoDif**2/rhoTot**2&
               +2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               +4.0*q4*rhoDif**3/rhoTot**3&
               )/(3.0*beta**2)

          dBdP=(q0&
               -q2*rhoDif**2/rhoTot**2&
               -2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               -4.0*q4*rhoDif**3/rhoTot**3&
               )/(3.0*beta**2)

          d2BdE2=-2.0*(q0&
               -q2*rhoDif**2/rhoTot**2&
               +2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               +4.0*q4*rhoDif**3/rhoTot**3&
               )**2/(9.0*beta**5)&
               +(2.0*q2/rhoTot&
               -4.0*q2*rhoDif/rhoTot**2&
               +(2.0*q2+12.0*q4)*rhoDif**2/rhoTot**3&
               -24.0*q4*rhoDif**3/rhoTot**4&
               +12.0*q4*rhoDif**4/rhoTot**5&
               )/(3.0*beta**2)

          d2BdP2=-2.0*(q0&
               -q2*rhoDif**2/rhoTot**2&
               -2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               -4.0*q4*rhoDif**3/rhoTot**3&
               )**2/(9.0*beta**5)&
               +(2.0*q2/rhoTot&
               +4.0*q2*rhoDif/rhoTot**2&
               +(2.0*q2+12.0*q4)*rhoDif**2/rhoTot**3&
               +24.0*q4*rhoDif**3/rhoTot**4&
               +12.0*q4*rhoDif**4/rhoTot**5&
               )/(3.0*beta**2)

          d2BdEP=-(2.0*(q0&
               -q2*rhoDif**2/rhoTot**2&
               +2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               +4.0*q4*rhoDif**3/rhoTot**3&
               )*(q0&
               -q2*rhoDif**2/rhoTot**2&
               -2.0*q2*rhoDif/rhoTot&
               -3.0*q4*rhoDif**4/rhoTot**4&
               -4.0*q4*rhoDif**3/rhoTot**3&               
               ))/(9.0*beta**5)&
               +(-2.0*q2/rhoTot&
               +(2.0*q2-12.0*q4)*rhoDif**2/rhoTot**3&
               +12.0*q4*rhoDif**4/rhoTot**5&
               )/(3.0*beta**2)

       else if( rhoE .gt. densityThreshold ) then !
          beta=((q0+q2+q4)*rhoE)**(1.0/3.0)
          dBdE=1.0/3.0*(q0+q2+q4)**(1.0/3.0)*rhoE**(-2.0/3.0)
          dBdP=0.0
          d2BdE2=-2.0/9.0*(q0+q2+q4)**(1.0/3.0)*rhoE**(-5.0/3.0)
          d2BdP2=0.0
          d2BdEP=0.0
       else if( rhoP .gt. densityThreshold ) then !
          beta=((q0+q2+q4)*rhoP)**(1.0/3.0)
          dBdE=0.0
          dBdP=1.0/3.0*(q0+q2+q4)**(1.0/3.0)*rhoP**(-2.0/3.0)
          d2BdE2=0.0
          d2BdP2=-2.0/9.0*(q0+q2+q4)**(1.0/3.0)*rhoP**(-5.0/3.0)
          d2BdEP=0.0
       else
          beta=0.0
          dBdE=0.0
          dBdP=0.0
          d2BdE2=0.0
          d2BdP2=0.0
          d2BdEP=0.0
       end if

    case("PsBeta")
       if(CONTROL_instance%BETA_PARAMETER_A .ne. 0.0 .and. CONTROL_instance%BETA_PARAMETER_B .ne. 0.0) then
          Eab=CONTROL_instance%BETA_PARAMETER_A
          cutOff=CONTROL_instance%BETA_PARAMETER_B
       else
          if( positiveMass .gt. 2.0) then
             Eab=-0.5
          else
             Eab=-0.25
          end if
          cutOff=0.1
       end if

       q0=-functionalLimitConstant/2/Eab
       ! print *, q0
       
       if( rhoTot .gt. densityThreshold) then !
       ! if( rhoTot .gt. densityThreshold .and. Sqrt(rhoDif**2) .gt. densityThreshold ) then !
          
          ! beta=(q0*(rhoTot+Sqrt(rhoDif**2)))**(1.0/3.0)

          ! dBdE=q0*(1.0+rhoDif/Sqrt(rhoDif**2))/&
          !      (3.0*beta**2)

          ! dBdP=q0*(1.0-rhoDif/Sqrt(rhoDif**2))/&
          !      (3.0*beta**2)

          ! d2BdE2=-4.0*q0**2*(1.0+rhoDif/Sqrt(rhoDif**2))/&
          !      (9.0*beta**5)     
          
          ! d2BdP2=-4.0*q0**2*(1.0-rhoDif/Sqrt(rhoDif**2))/&
          !      (9.0*beta**5)     

          beta=(q0*(rhoTot+rhoDif*tanh(rhoDif/rhoTot/cutOff)))**(1.0/3.0)

          dBdE=q0*(2.0*rhoDif*rhoP*(1.0/cosh(rhoDif/rhoTot/cutOff))**2&
               +rhoTot**2*cutOff*(1.0+tanh(rhoDif/rhoTot/cutOff)))/&
               (3.0*beta**2*rhoTot**2*cutOff)

          dBdP=-q0*(2.0*rhoDif*rhoE*(1.0/cosh(rhoDif/rhoTot/cutOff))**2&
               +rhoTot**2*cutOff*(-1.0+tanh(rhoDif/rhoTot/cutOff)))/&
               (3.0*beta**2*rhoTot**2*cutOff)

          d2BdE2=0.0
          ! -4.0*q0**2*(1.0+rhoDif/Sqrt(rhoDif**2))/&
          !      (9.0*beta**5)     
          
          d2BdP2=0.0
          ! -4.0*q0**2*(1.0-rhoDif/Sqrt(rhoDif**2))/&
          !      (9.0*beta**5)     

          ! beta=(q0*rhoTot*(1.0+Sqrt(rhoDif**2)/rhoTot))**(1.0/3.0)

          ! dBdE=q0*2.0/(1.0+exp(-logisticExp*rhoDif))/&
          !      (3.0*beta**2)

          ! dBdP=q0*2.0/(1.0+exp(logisticExp*rhoDif))/&
          !      (3.0*beta**2)

          ! d2BdE2=-4.0*q0**2*2.0/(1.0+exp(-logisticExp*rhoDif))/&
          !      (9.0*beta**5)     
          
          ! d2BdP2=-4.0*q0**2*2.0/(1.0+exp(logisticExp*rhoDif))/&
          !      (9.0*beta**5)     
          
          ! d2BdEP=0.0
       ! else if( rhoTot .gt. densityThreshold .and. Sqrt(rhoDif**2) .lt. densityThreshold ) then !
       !    beta=(q0*rhoTot)**(1.0/3.0)
       !    dBdE=0.0
       !    dBdP=0.0
       !    d2BdE2=0.0
       !    d2BdP2=0.0
       !    d2BdEP=0.0

       else if( rhoE .gt. densityThreshold ) then !
          beta=(2.0*q0*rhoE)**(1.0/3.0)
          dBdE=1.0/3.0*(2.0*q0)**(1.0/3.0)*rhoE**(-2.0/3.0)
          dBdP=0.0
          d2BdE2=-2.0/9.0*(2.0*q0)**(1.0/3.0)*rhoE**(-5.0/3.0)
          d2BdP2=0.0
          d2BdEP=0.0
       else if( rhoP .gt. densityThreshold ) then !
          beta=(2.0*q0*rhoP)**(1.0/3.0)
          dBdE=0.0
          dBdP=1.0/3.0*(2.0*q0)**(1.0/3.0)*rhoP**(-2.0/3.0)
          d2BdE2=0.0
          d2BdP2=-2.0/9.0*(2.0*q0)**(1.0/3.0)*rhoP**(-5.0/3.0)
          d2BdEP=0.0
       else
          beta=0.0
          dBdE=0.0
          dBdP=0.0
          d2BdE2=0.0
          d2BdP2=0.0
          d2BdEP=0.0
       end if

       ! print *, rhoE, rhoP, Sqrt(rhoDif**2), beta, dBdE, dBdP
       
    case("PsBetaMax")
       if(CONTROL_instance%BETA_PARAMETER_A .ne. 0.0 ) then
          Eab=CONTROL_instance%BETA_PARAMETER_A
       else
          if( positiveMass .gt. 2.0) then
             Eab=-0.5
          else
             Eab=-0.25
          end if
       end if

       q0=(-functionalLimitConstant/Eab)**(1.0/3.0)
       ! print *, q0

       if( rhoTot .gt. densityThreshold) then !

          if( abs(rhoDif) .lt. densityThreshold) then !could be .lt. cutoff

             beta=q0*(rhoTot/2.0)**(1.0/3.0)
             dBdE=q0**3/6.0/beta**2
             dBdP=q0**3/6.0/beta**2
             d2BdE2=-q0**6/18.0/beta**5
             d2BdP2=-q0**6/18.0/beta**5
             d2BdEP=-q0**6/18.0/beta**5

          elseif (rhoE .gt. rhoP) then

             beta=q0*rhoE**(1.0/3.0)
             dBdE=q0**3/3.0/beta**2
             dBdP=0.0
             d2BdE2=-2.0*q0**6/9.0/beta**5
             d2BdP2=0.0
             d2BdEP=0.0

          elseif (rhoP .gt. rhoE) then

             beta=q0*rhoP**(1.0/3.0)
             dBdE=0.0
             dBdP=q0**3/3.0/beta**2
             d2BdE2=0.0
             d2BdP2=-2.0*q0**6/9.0/beta**5
             d2BdEP=0.0
          end if

       else
          beta=0.0
          dBdE=0.0
          dBdP=0.0
          d2BdE2=0.0
          d2BdP2=0.0
          d2BdEP=0.0
       end if

       ! print *, rhoE, rhoP, beta
    case default

    end select

    ! if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) then
    !    qe=CONTROL_instance%DUMMY_REAL(1)
    !    qn=CONTROL_instance%DUMMY_REAL(2)
    !    if(CONTROL_instance%BETA_FUNCTION .eq. "rhoE3rhoN3rhoEN6") then
    !       qen=CONTROL_instance%DUMMY_REAL(3)
    !       p=1
    !    else
    !       p=CONTROL_instance%DUMMY_REAL(3)
    !    end if
    ! else
    !    qe=1.1487573585337585
    !    qn=1.1487573585337585
    !    p=1.0
    ! end if

    ! if(CONTROL_instance%BETA_FUNCTION .eq. "newBeta") then

    !    if(mass .gt. 2.0) then !hydrogen
    !       Ea2b=0.527444
    !       Eab2=0.597139 
    !    else !positron
    !       Ea2b=0.262005
    !       Eab2=0.262005
    !    end if
    !    p=1.0
       
    !    if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) then
    !       Ea2b=CONTROL_instance%DUMMY_REAL(1)
    !       Eab2=CONTROL_instance%DUMMY_REAL(2)
    !       p=CONTROL_instance%DUMMY_REAL(3)
    !    end if

    !    qe=a0*(11.0/8.0/Ea2b-3.0/4.0/Eab2)
    !    qn=a0*(11.0/8.0/Eab2-3.0/4.0/Ea2b)
    !    q2en=a0*3.0/16.0*(1.0/Ea2b+1.0/Eab2)
    !    q3en=a0*9.0/16.0*(1.0/Eab2-1.0/Ea2b)
       
    ! else if(CONTROL_instance%BETA_FUNCTION .eq. "newnewBeta") then

    !    if(this%mass2 .gt. 2.0) then !hydrogen
    !       STOP "this beta function only works for electron-positron"
    !    else !positron
    !       Eab=0.25
    !       Eab2=0.2620050702329801
    !    end if

    !    p=1.0

    !    if(CONTROL_instance%DUMMY_REAL(1) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(2) .ne. 0 .or. CONTROL_instance%DUMMY_REAL(3) .ne. 0) then
    !       Eab=CONTROL_instance%DUMMY_REAL(1)
    !       Eab2=CONTROL_instance%DUMMY_REAL(2)
    !       p=CONTROL_instance%DUMMY_REAL(3)
    !    end if

    !    q0=a0/2/Eab
    !    q2=a0*(-5/Eab+53/Eab2/8)
    !    q4=a0*(9/Eab/2-45/Eab2/8)
    ! end if

  end subroutine Functional_getBeta
  
  subroutine Functional_PSNEvaluate( this, mass, n, rhoE, rhoP, ec, vcE, vcP )
    ! Evaluates Puska-Seitsonen-Nieminen electron-positron correlation functional
    ! Felix Moncada, 2017
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoP(*) !! electron and nuclear Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcP(*) !! Potentials - output   

    real(8) :: Aa,Ba,Ca,Bb,Cb,Cc,rse,rsp,drho_rse,drho_rsp
    real(8) :: denominator
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
    
    ! densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD    

    ! print *, "i, rhoE, rhoN, denominator, energy density, potentialE, potentialN"
    do i = 1, n

       rse= (3.0/(4.0*Math_PI*rhoE(i)))**(1.0/3.0)
       rsp= (3.0/(4.0*Math_PI*rhoP(i)))**(1.0/3.0)
       
       drho_rse= -1/(rhoE(i)**(4.0/3.0)*6.0**(2.0/3.0)*Math_PI**(1.0/3.0))
       drho_rsp= -1/(rhoP(i)**(4.0/3.0)*6.0**(2.0/3.0)*Math_PI**(1.0/3.0))
       !(1.0/3.0)*(3.0/(4.0*Math_PI))**(1.0/3.0)

       !Energy and potential

       !This should be a parameter in CONTROL
       if( rse .lt. 8.0 .and. rsp .lt. 8.0) then
          denominator= Aa+Ba*(rse+rsp)+Ca*(rse**2+rsp**2) &
               +Bb*rse*rsp+Cb*(rse**2 *rsp + rse*rsp**2) + Cc*rse**2 *rsp**2 &
               +4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          !Filter positive values that may be problematic
          if(denominator .lt. 0.0 ) then
             ec(i)=ec(i) + 1/denominator

             vcE(i)=vcE(i) - drho_rse/denominator**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2 &
                  +4*Math_PI*rse**2/eap_homogeneus(rsp)&
                  +Ba+2*Ca*rse+Bb*rsp+Cb*(2*rse*rsp + rsp**2)+2*Cc*rse*rsp**2)
             vcP(i)=vcP(i) - drho_rsp/denominator**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2 &
                  +4*Math_PI*rsp**2/eap_homogeneus(rse)&
                  +Ba+2*Ca*rsp+Bb*rse+Cb*(rse**2 + 2*rse*rsp)+2*Cc*rse**2*rsp)
          end if
       else if( rse .lt. 20 .or. rsp .lt. 20) then            
          denominator= 4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          if(denominator .lt. 0.0 ) then
             ec(i)=ec(i) + 1/denominator
             vcE(i)=vcE(i) - drho_rse/denominator**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2 &
                  +4*Math_PI*rse**2/eap_homogeneus(rsp))
             vcP(i)=vcP(i) - drho_rsp/denominator**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2 &
                  +4*Math_PI*rsp**2/eap_homogeneus(rse))                  
          end if
       end if

       !We multiply again by rhoE in grid manager
       if(rhoE(i) .gt. 0.0) then
          ec(i)=ec(i)/rhoE(i)
       else
          ec(i)=0.0
       end if

       ! print *, i, rhoE(i), rhoP(i), rse, rsp, ec(i)*rhoE(i), vcE(i), vcP(i)
       
       !Potential
       
       ! denominator= 4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)
                    
       ! print *, i, rse, rsp , ec(i), vcE(i), vcP(i) !, vcE(i), - drho_rse/denominator**2*(4*Math_PI*rse**2/eap_homogeneus(rsp)), rse**2, eap_homogeneus(rsp)
     !drho_rse, denominator**2, dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2, 4*Math_PI*rse**2/eap_homogeneus(rsp)
       
       ! if(denominator .ge. 0.0) then
       !    print *, "no jodas, en serio?"
       !    print *, i, rse, rsp, eap_homogeneus(rse), eap_homogeneus(rsp), denominator
       ! end if
                 
    end do
    
  end subroutine Functional_PSNEvaluate

  subroutine Functional_PSNAPEvaluate( this, mass, n, rhoE, rhoP, ec, vcE, vcP )
    ! Evaluates Puska-Seitsonen-Nieminen electron-positron correlation functional for high re and rp values
    ! Evaluates Arponen-Pajanne single particle limit, as interpolated by Boronsky and Nieminen of the electron-positron correlation for low  re or rp values
    ! Includes a smooth switch function at a r cutoff
    ! Felix Moncada, 2022
    implicit none
    type(Functional):: this !!type of functional
    real(8) :: mass !!nuclear mass
    integer :: n !!nuclear gridSize
    real(8) :: rhoE(*), rhoP(*) !! electron and positron Densities - input
    real(8) :: ec(*) !! Energy density - output
    real(8) :: vcE(*), vcP(*) !! Potentials - output   

    real(8) :: Aa,Ba,Ca,Bb,Cb,Cc,rse,rsp,drho_rse,drho_rsp
    real(8) :: rcut, xcut !switch parameters
    real(8) :: f, dfde, dfdp !switch results
    real(8) :: denominatorPSN, EPSN, vcEPSN, vcPPSN !PSN results     
    real(8) :: denominatorAP, EAP, vcEAP, vcPAP !AP results     
    real(8) :: densityThreshold
    integer :: i

    !!The idea is that the parameters are a functional of the nuclear mass and charge
    if(this%name .eq. "correlation:psnap" ) then
       !*2 to convert from Rydbergs to a.u.
        Aa=69.7029*2.0_8
        Ba=-107.4927*2.0_8
        Bb=141.8458*2.0_8
        Ca=23.7182*2.0_8
        Cb=-33.6472*2.0_8
        Cc=5.21152*2.0_8
        rcut=8.0
        xcut=6.0
    else
       print *, this%name
       STOP "The nuclear electron functional chosen is not implemented"
    end if
    
    densityThreshold=CONTROL_instance%NUCLEAR_ELECTRON_DENSITY_THRESHOLD    

    ! print *, "i, rhoE, rhoN, denominator, energy density, potentialE, potentialN"
    do i = 1, n

       EPSN=0.0
       vcEPSN=0.0
       vcPPSN=0.0
       
       if(rhoE(i) .gt. densityThreshold .and. rhoP(i) .gt. densityThreshold) then
          rse= (3.0/(4.0*Math_PI*rhoE(i)))**(1.0/3.0)
          rsp= (3.0/(4.0*Math_PI*rhoP(i)))**(1.0/3.0)

          drho_rse= -1/(rhoE(i)**(4.0/3.0)*6.0**(2.0/3.0)*Math_PI**(1.0/3.0))
          drho_rsp= -1/(rhoP(i)**(4.0/3.0)*6.0**(2.0/3.0)*Math_PI**(1.0/3.0))

          !switch function
          f=exp(-(rse**xcut+rsp**xcut)/rcut**xcut)
          dfde =drho_rse*( -(xcut*(rse/rcut)**xcut)/rse)*f
          dfdp =drho_rsp*( -(xcut*(rsp/rcut)**xcut)/rsp)*f

          !Energy and potential

          !PSN part
          !This should be a parameter in CONTROL
          denominatorPSN= Aa+Ba*(rse+rsp)+Ca*(rse**2+rsp**2) &
               +Bb*rse*rsp+Cb*(rse**2 *rsp + rse*rsp**2) + Cc*rse**2 *rsp**2 &
               +4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          !Filter positive values that may be problematic
          if(denominatorPSN .lt. 0.0 ) then
             EPSN=1/denominatorPSN

             vcEPSN=-drho_rse/denominatorPSN**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2 &
                  +4*Math_PI*rse**2/eap_homogeneus(rsp)&
                  +Ba+2*Ca*rse+Bb*rsp+Cb*(2*rse*rsp + rsp**2)+2*Cc*rse*rsp**2)
             vcPPSN=-drho_rsp/denominatorPSN**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2 &
                  +4*Math_PI*rsp**2/eap_homogeneus(rse)&
                  +Ba+2*Ca*rsp+Bb*rse+Cb*(rse**2 + 2*rse*rsp)+2*Cc*rse**2*rsp)
          end if

          !AP part
          denominatorAP= 4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse) + 4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)

          if(denominatorAP .lt. 0.0 ) then
             EAP=1/denominatorAP
             vcEAP= - drho_rse/denominatorAP**2*(-dr_eap_homogeneus(rse)*4.0/3.0*Math_PI*rsp**3/eap_homogeneus(rse)**2 &
                  +4*Math_PI*rse**2/eap_homogeneus(rsp))
             vcPAP= - drho_rsp/denominatorAP**2*(-dr_eap_homogeneus(rsp)*4.0/3.0*Math_PI*rse**3/eap_homogeneus(rsp)**2 &
                  +4*Math_PI*rsp**2/eap_homogeneus(rse))                  
          end if

          ec(i)=ec(i) + EPSN*f + EAP*(1-f)
          vcE(i)=vcE(i) + vcEPSN*f + vcEAP*(1-f) + EPSN*dfde - EAP*dfde
          vcP(i)=vcP(i) + vcPPSN*f + vcPAP*(1-f) + EPSN*dfdp - EAP*dfdp

       end if
       !We multiply again by rhoE in grid manager
       if(rhoE(i) .gt. 0.0) then
          ec(i)=ec(i)/rhoE(i)
       else
          ec(i)=0.0
       end if

       ! print *, i, rhoE(i), rhoP(i), rse, rsp, ec(i)*rhoE(i), vcE(i), vcP(i)
       
    end do
    
  end subroutine Functional_PSNAPEvaluate
  

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

    factor = 1.0/3.0
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
    ! Evaluates Boronsky and Nieminen positron-electron correlation single particle limit energy
    ! Felix Moncada, 2017
    implicit none
    real(8) rs

    !!This expressions are in rydbergs
    if(rs .lt. 0.302) eap_homogeneus=-1.56/sqrt(rs)+(0.051*log(rs)-0.081)*log(rs)+1.14
    if(rs .ge. 0.302 .and. rs .lt. 0.56) eap_homogeneus=-0.92305-0.05459/rs**2.0
    if(rs .ge. 0.56 .and. rs .lt. 8.0) eap_homogeneus=-13.15111/(rs+2.5)**2.0 + 2.8655/(rs+2.5) - 0.6298
    if(rs .ge. 8.0) eap_homogeneus=-10250.57860/rs**6.0  + 44.50466/rs**3.0-0.524

    !!changing to a.u.
    eap_homogeneus=eap_homogeneus/2.0_8
  end function eap_homogeneus

  real(8) function dr_eap_homogeneus(rs)
    ! Evaluates Boronsky and Nieminen positron-electron correlation single particle limit potential
    ! Felix Moncada, 2017
    implicit none
    real(8) rs

    !!This expressions are in rydbergs
    if(rs .lt. 0.302) dr_eap_homogeneus=0.78/rs**(3.0/2.0)+(0.102*log(rs)-0.081)/rs
    if(rs .ge. 0.302 .and. rs .lt. 0.56) dr_eap_homogeneus=0.10918/rs**3.0
    if(rs .ge. 0.56 .and. rs .lt. 8.0) dr_eap_homogeneus=26.30222/(rs+2.5)**3.0 - 2.8655/(rs+2.5)**2.0
    if(rs .ge. 8.0) dr_eap_homogeneus=61503.4716/rs**7.0  -133.51399/rs**4.0

    !!changing to a.u.
    dr_eap_homogeneus=dr_eap_homogeneus/2.0_8
    
  end function dr_eap_homogeneus

  
end module Functional_

