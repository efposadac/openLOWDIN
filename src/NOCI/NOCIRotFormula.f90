!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                
!!	
!!	this program has been developed under direction of:
!!
!!	  UNIVERSIDAD NACIONAL DE COLOMBIA"
!!	  PROF. ANDRES REYES GROUP"
!!	  http://www.qcc.unal.edu.co"
!!	
!!	  UNIVERSIDAD DE GUADALAJARA"
!!	  PROF. ROBERTO FLORES GROUP"
!!	  http://www.cucei.udg.mx/~robertof"
!!	
!!	AUTHORS
!!		E.F. POSADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		S.A. GONZALEZ. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		F.S. MONCADA. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		J. ROMERO. UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!	CONTRIBUTORS
!!		N.F.AGUIRRE. UNIVERSIDAD NACIONAL DE COLOMBIA
!!   		GABRIEL MERINO. UNIVERSIDAD DE GUANAJUATO
!!   		J.A. CHARRY UNIVERSIDAD NACIONAL DE COLOMBIA
!!
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

module NOCIRotFormula_
  use NOCIBuild_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use Math_
  use DirectIntegralManager_
  use Libint2Interface_
  use omp_lib
  implicit none

  !>
  !! @brief non Orthogonal Configuration Interaction Module. APMO implementation of Skone et al 2005 10.1063/1.2039727
  !!
  !! @author Felix
  !!
  !! <b> Creation data : </b> 02-22
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 02-22 </tt>: Felix Moncada ( fsmoncadaa@unal.edu.co )
  !!        -# Creation of the module.
  !!
  !<

  public :: &
       NOCIRotFormula_compute

  private

contains

  !>
  !! @brief Computes the rotational CI energies from previously computed overlap and hamiltonian non orthogonal CI elements
  !!
  !! @param this
  !<
  subroutine NOCIRotFormula_compute(this)
    implicit none
    type(NonOrthogonalCI) :: this
    type(Vector) :: angles, signs
    type(Matrix) :: weights
    integer :: i,state,sysI,npoints,nstates,speciesID,otherSpeciesID
    real(8) :: overlapIntegral, auxEnergyIntegral, auxEnergy, e0, sc

    ! real(8) :: timeMerging, timePrescreen, timeOverlap, timeTwoIntegrals
    ! real(8) :: timeA
    ! real(8) :: timeB

    if((this%transformationType).ne."ROTATION_AROUND_Z") then
       print *, "The Rotational Configuration Interaction formula for the rotational states energy is only available for molecular systems rotated around the z-axis"
       print *, "Please set rotationalScanGridAroundZ=N in the input and restart the calculation"
    end if

    nstates=min(CONTROL_instance%NUMBER_OF_CI_STATES,this%numberOfDisplacedSystems)
    if(nstates .lt. 2) nstates=2
    npoints=this%numberOfIndividualTransformations
    
    call Vector_constructor(angles,npoints,0.0_8)    
    call Matrix_constructor(weights,int(npoints,8),int(nstates,8),1.0_8)
    call Vector_constructor(signs,this%numberOfDisplacedSystems,1.0_8)    
    call Vector_constructor(this%statesEigenvalues, this%numberOfDisplacedSystems, 0.0_8)

    do i=1,npoints
       angles%values(i)=(i-1)*CONTROL_instance%ROTATION_AROUND_Z_STEP*Math_PI/180
    end do

    if(this%molecularSystems(1)%numberOfPointCharges .gt. 1) then
       print *, "Using 1D formula: cos(m gamma) as weights, with trapezoid integration rule"
       weights%values(1,:)=0.5_8
       do state=1,nstates
          weights%values(2:npoints,state)=cos((state-1)*angles%values(2:npoints))
       end do
       weights%values(npoints,:)=0.5_8*weights%values(npoints,:)
    else
       print *, "Using 2D formula: sin(gamma) P_l(cos(gamma)) as weights, with trapezoid integration rule"
       ! weights%values(1,:)=0.5_8
       call Math_p_polynomial_value (npoints , nstates-1, cos(angles%values), weights%values)
       call flush()
       do i=1,npoints
          weights%values(i,:)=sin(angles%values(i))*weights%values(i,:)
       end do
       weights%values(npoints,:)=0.5_8*weights%values(npoints,:)
    end if

    write(*,"(A)") ""
    write(*,"(A)") " MIXED HARTREE-FOCK CALCULATION"
    write(*,"(A)") " NON ORTHOGONAL CONFIGURATION INTERACTION"
    write(*,"(A)") " ROTATIONAL CI FORMULA UNSCALED ENERGIES: "
    write(*,"(A)") "========================================="
    write(*,"(A)") ""

    do state=1,nstates
       overlapIntegral=0

       do sysI=1,this%numberOfDisplacedSystems
          signs%values(sysI)=this%configurationOverlapMatrix%values(1,sysI)/abs(this%configurationOverlapMatrix%values(1,sysI))
          overlapIntegral=overlapIntegral+signs%values(sysI)*this%configurationOverlapMatrix%values(1,sysI)*weights%values(sysI,state)
       end do

       auxEnergyIntegral=0
       do sysI=1,this%numberOfDisplacedSystems
          auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationHamiltonianMatrix%values(1,sysI)*weights%values(sysI,state)
       end do
       this%statesEigenvalues%values(state)=auxEnergyIntegral/overlapIntegral
       ! print *, "state", state, "overlapIntegral", overlapIntegral, "energyIntegral", auxEnergyIntegral, "energy", auxEnergyIntegral/overlapIntegral
       write (*,"(A)")  ""
       write (*,"(T9,A17,I3,A10, F25.12)") "STATE: ", state, " ENERGY = ", this%statesEigenvalues%values(state)

       write (*,"(A38)")  "Components: "
       write(*,"(A38,F25.12)") " Point charges energy = ", MolecularSystem_getPointChargesEnergy(this%molecularSystems(1))

       do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
          auxEnergyIntegral=0
          do sysI=1,this%numberOfDisplacedSystems
             auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationKineticMatrix(speciesID)%values(1,sysI)*weights%values(sysI,state)
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Kinetic energy = ", auxEnergyIntegral/overlapIntegral
       end do

       do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
          auxEnergyIntegral=0
          do sysI=1,this%numberOfDisplacedSystems
             auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationPuntualMatrix(speciesID)%values(1,sysI)*weights%values(sysI,state)
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Puntual energy = ", auxEnergyIntegral/overlapIntegral
       end do

       if(CONTROL_instance%IS_THERE_EXTERNAL_POTENTIAL) then
          do speciesID = 1, this%molecularSystems(1)%numberOfQuantumSpecies                
             auxEnergyIntegral=0
             do sysI=1,this%numberOfDisplacedSystems
                auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationExternalMatrix(speciesID)%values(1,sysI)*weights%values(sysI,state)
             end do
             write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                  " External energy = ", auxEnergyIntegral/overlapIntegral
          end do
       end if

       do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
          auxEnergyIntegral=0
          do sysI=1,this%numberOfDisplacedSystems
             auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationHartreeMatrix(speciesID,speciesID)%values(1,sysI)*weights%values(sysI,state)
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               "/"//trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Hartree energy = ", auxEnergyIntegral/overlapIntegral
       end do

       do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
          auxEnergyIntegral=0
          do sysI=1,this%numberOfDisplacedSystems
             auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationExchangeMatrix(speciesID)%values(1,sysI)*weights%values(sysI,state)
          end do
          write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
               " Exchange energy = ", auxEnergyIntegral/overlapIntegral

       end do

       do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies-1
          do otherSpeciesID=speciesID+1, this%molecularSystems(1)%numberOfQuantumSpecies
             auxEnergyIntegral=0
             do sysI=1,this%numberOfDisplacedSystems
                auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationHartreeMatrix(speciesID,otherSpeciesID)%values(1,sysI)*weights%values(sysI,state)
             end do
             write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                  "/"//trim( this%molecularSystems(1)%species(otherSpeciesID)%name ) // &
                  " Hartree energy = ", auxEnergyIntegral/overlapIntegral
          end do
       end do
       if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
          do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
             do otherSpeciesID=speciesID, this%molecularSystems(1)%numberOfQuantumSpecies
                auxEnergyIntegral=0
                do sysI=1,this%numberOfDisplacedSystems
                   auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(1,sysI)*weights%values(sysI,state)
                end do
                write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                     "/"//trim( this%molecularSystems(1)%species(otherSpeciesID)%name ) // &
                     " DFTcorrelation energy = ", auxEnergyIntegral/overlapIntegral
             end do
          end do
       end if
       
       write(*,"(A)") ""

    end do

    if ( CONTROL_instance%METHOD .eq. "RKS" .or. CONTROL_instance%METHOD .eq. "UKS" ) then
       write(*,"(A)") ""
       write(*,"(A)") " MIXED HARTREE-FOCK CALCULATION"
       write(*,"(A)") " NON ORTHOGONAL CONFIGURATION INTERACTION"
       write(*,"(A)") " ROTATIONAL CI FORMULA SCALED ENERGIES: "
       write(*,"(A)") "========================================="
       write(*,"(A)") ""

       print *, "Using a sigmoid function, e0+(1+e0)exp(-(1-|S|)^4/sc^4), to scale down the interspecies correlation functional energy"
       print *, "All the other energy contributions remain equal"
       print *, ""

       if(CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0 .ne. 0.0_8 .or. CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_SC .ne. 0.0_8 ) then
          sc=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_SC
          e0=CONTROL_instance%EMPIRICAL_OVERLAP_PARAMETER_E0
          write (*,'(A63)') "Employing sigmoid parameters provided in the input"
          write (*,'(A48,F15.8)') "sc =", sc
          write (*,'(A48,F15.8)') "e0 =", e0
          print *, ""
       else
          call NOCIRotFormula_getScalingParameters(this,sc,e0)
       end if
       
       do state=1,nstates
          overlapIntegral=0
          auxEnergy=0
          do sysI=1,this%numberOfDisplacedSystems
             overlapIntegral=overlapIntegral+signs%values(sysI)*this%configurationOverlapMatrix%values(1,sysI)*weights%values(sysI,state)
          end do

          do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
             do otherSpeciesID=speciesID+1, this%molecularSystems(1)%numberOfQuantumSpecies
                auxEnergyIntegral=0
                do sysI=1,this%numberOfDisplacedSystems
                   auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(1,sysI)*weights%values(sysI,state)*&
                   ((e0+(1-e0)*exp(-(1-abs(this%configurationOverlapMatrix%values(1,sysI)))**4/sc**4))-1)
                        
                end do
                write(*,"(A38,F25.12)") trim( this%molecularSystems(1)%species(speciesID)%name ) // &
                     "/"//trim( this%molecularSystems(1)%species(otherSpeciesID)%name ) // &
                     " DFTcorrection energy = ", auxEnergyIntegral/overlapIntegral
                auxEnergy=auxEnergy+auxEnergyIntegral/overlapIntegral
             end do
          end do
          
          auxEnergyIntegral=0
          do sysI=1,this%numberOfDisplacedSystems
             auxEnergyIntegral=auxEnergyIntegral+signs%values(sysI)*this%configurationHamiltonianMatrix%values(1,sysI)*weights%values(sysI,state)
          end do
          write (*,"(T9,A10,I3,A17, F25.12)") "STATE: ", state, "SCALED ENERGY = ", auxEnergyIntegral/overlapIntegral+auxEnergy
          write (*,"(A)")  ""
       end do

    end if
    
  end subroutine NOCIRotFormula_compute

  !>
  !! @brief Get the empirical scaling parameters from the multilinear regression 
  !!
  !! @param this
  !<
  subroutine NOCIRotFormula_getScalingParameters(this,sc,e0)
    implicit none
    type(NonOrthogonalCI) :: this
    real(8), intent(out) :: sc,e0

    type(Matrix) :: momentMatrices(1:3), densityMatrix
    integer :: i,j,k,speciesID,otherSpeciesID
    real(8) :: b(0:6,1:2) !1 for sc, 2 for e0
    real(8) :: x(1:6)

!Regression parameters
! Intercept
! $x_1$:	dim
! $x_2$:	DFT $T_p$
! $x_3$:	$\langle r_p \rangle$
! $x_4$:	$-E^c_{ep}$
! $x_5$:	E$_{0\ \mathrm{RoDFT}}-$E$_\mathrm{DFT}$
! $x_6$:	$\Delta E_{0-1 \mathrm{RoDFT}}$
    b(0,1)=0.3881
    b(1,1)=-0.0802
    b(2,1)=17.22
    b(3,1)=0.00948
    b(4,1)=-7.41
    b(5,1)=16.21
    b(6,1)=122.3

    b(0,2)=-0.6042
    b(1,2)=0.5374
    b(2,2)=-34.29
    b(3,2)=-0.02577
    b(4,2)=37.33
    b(5,2)=-65.85
    b(6,2)=-256.5

    x(:)=0.0
    
    x(1)=2
    if(this%molecularSystems(1)%numberOfPointCharges .gt. 1) x(1)=1

    do speciesID=1, this%molecularSystems(1)%numberOfQuantumSpecies
       !find proton kinetic energy
       if(this%molecularSystems(1)%species(speciesID)%mass .ge. 1500_8 .and. this%molecularSystems(1)%species(speciesID)%mass .lt. 2500_8 .and. this%molecularSystems(1)%species(speciesID)%charge .eq. 1.0_8) then
          x(2)=this%configurationKineticMatrix(speciesID)%values(1,1)

          call DirectIntegralManager_getMomentIntegrals(this%molecularSystems(1),speciesID,1,momentMatrices(1))
          call DirectIntegralManager_getMomentIntegrals(this%molecularSystems(1),speciesID,2,momentMatrices(2))
          call DirectIntegralManager_getMomentIntegrals(this%molecularSystems(1),speciesID,3,momentMatrices(3))

          call Matrix_constructor(densityMatrix,int(size( this%HFCoefficients(1,speciesID)%values, DIM = 1),8),int(size( this%HFCoefficients(1,speciesID)%values, DIM = 1),8),0.0_8)
          
          do i = 1 , size( this%HFCoefficients(1,speciesID)%values, DIM = 1 )
             do j = 1 , size( this%HFCoefficients(1,speciesID)%values, DIM = 1 )
                do k = 1 , MolecularSystem_getOcupationNumber(speciesID,this%molecularSystems(1))

                   densityMatrix%values(i,j) =  &
                        densityMatrix%values( i,j ) + &
                        ( this%HFCoefficients(1,speciesID)%values(i,k) &
                        * this%HFCoefficients(1,speciesID)%values(j,k) )
                end do
             end do
          end do
          densityMatrix%values =  MolecularSystem_getEta(speciesID,this%molecularSystems(1))  * densityMatrix%values
          
          if(this%molecularSystems(1)%numberOfPointCharges .gt. 1) then
             x(3)=sqrt(sum( densityMatrix%values * momentMatrices(1)%values ) **2 +&
                  sum( densityMatrix%values * momentMatrices(2)%values ) **2)
             
          else
             x(3)=sqrt(sum( densityMatrix%values * momentMatrices(1)%values ) **2 +&
                  sum( densityMatrix%values * momentMatrices(2)%values ) **2 +&
                  sum( densityMatrix%values * momentMatrices(3)%values ) **2)
          end if

          do otherSpeciesID=1,speciesID-1
             x(4)=x(4)+this%configurationDFTcorrelationMatrix(otherSpeciesID,speciesID)%values(1,1)
          end do
          do otherSpeciesID=speciesID+1,this%molecularSystems(1)%numberOfQuantumSpecies
             x(4)=x(4)+this%configurationDFTcorrelationMatrix(speciesID,otherSpeciesID)%values(1,1)
          end do
          x(4)=-x(4)
          exit
       end if
    end do

    x(5)=this%configurationHamiltonianMatrix%values(1,1)-this%statesEigenvalues%values(1)
    x(6)=this%statesEigenvalues%values(2)-this%statesEigenvalues%values(1)
    
    sc=b(0,1)
    e0=b(0,2)

    do i=1,6
       sc=sc+b(i,1)*x(i)
       e0=e0+b(i,2)*x(i)
    end do

    if(sc<1.0E-8) sc=1.0E-8
    if(e0<0.0) e0=0.0
    if(e0>1.0) e0=1.0

    write (*,'(A63)') "The sigmoid parameters"
    write (*,'(A48,F15.8)') "e0 =", e0
    write (*,'(A48,F15.8)') "sc =", sc
    write (*,'(A63)') "were obtained from the regression parameters:"
    write (*,'(A48,I6)') "rotational dimensions =", int(x(1))
    write (*,'(A48,F15.8)') "proton kinetic energy =", x(2)
    write (*,'(A48,F15.8)') "proton rotation radius =", x(3)
    write (*,'(A48,F15.8)') "-proton/electron correlation energy =", x(4)
    write (*,'(A48,F15.8)') "-Unscaled rotational ground state correction =", x(5)
    write (*,'(A48,F15.8)') "Unscaled rotational first transition energy =", x(6)
    write (*,'(A63)') ""
    
  end subroutine NOCIRotFormula_getScalingParameters
  
end module NOCIRotFormula_

