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
!!        This module calculates the total energies with the Universal Force Field (UFF)
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions using Universal Force Field has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module EnergyUFF_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Graph_
  use Exception_
  implicit none

  type , public :: EnergyUFF
     
     real(8) :: totalStretchingEnergy
     real(8) :: totalStretchingEnergyKJ
     real(8) :: totalBendingEnergy
     real(8) :: totalBendingEnergyKJ
     real(8) :: totalTorsionEnergy
     real(8) :: totalTorsionEnergyKJ
     real(8) :: totalVDWEnergy
     real(8) :: totalVDWEnergyKJ
     real(8) :: totalInversionEnergy
     real(8) :: totalInversionEnergyKJ
     real(8) :: totalElectrostaticEnergy
     real(8) :: totalElectrostaticEnergyKJ
     real(8) :: totalEnergy
     real(8) :: totalEnergyKJ

  end type EnergyUFF
  
  public :: &
       EnergyUFF_run
       ! EnergyUFF_show

  !>Singleton
  type(EnergyUFF), public, target :: EnergyUFF_instance

contains
    
  !>
  !! @brief This routine calculates the total energies with the UFF
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class with all information about the Graph (System)
  !! @param [in] electrostatic LOGICAL evaluates if the user requires Electrostatic Energy
  !! @note The total energies are calculated using the potential energy of UFF, Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[ 
  !! E = E_{R} + E_{\theta} + E_{\phi} + E_{\omega} + E_{vdw} + E_{el}
  !! \f]
  !! where: \n
  !! - \f$E_{R}\f$ is the stretching energy
  !! - \f$E_{\theta}\f$ is the bending energy
  !! - \f$E_{\phi}\f$ is the torsion energy
  !! - \f$E_{\omega}\f$ is the inversion energy
  !! - \f$E_{vdw}\f$ is the Van der Waals energy
  !! - \f$E_{el}\f$ is the electrostatic energy (this energy is optional)
  subroutine EnergyUFF_run( this, electrostatic )
    implicit none
    type(Graph) :: this
    logical, intent(in) :: electrostatic
    integer :: i


    EnergyUFF_instance%totalStretchingEnergy = 0.0
    EnergyUFF_instance%totalBendingEnergy = 0.0
    EnergyUFF_instance%totalTorsionEnergy = 0.0
    EnergyUFF_instance%totalVDWEnergy = 0.0
    EnergyUFF_instance%totalInversionEnergy = 0.0
    EnergyUFF_instance%totalElectrostaticEnergy = 0.0

    do i=1, this%edges%numberOfEdges
       EnergyUFF_instance%totalStretchingEnergy = EnergyUFF_instance%totalStretchingEnergy + this%edges%stretchingEnergy(i)
    end do

    do i=1, this%angles%numberOfAngles
       EnergyUFF_instance%totalBendingEnergy = EnergyUFF_instance%totalBendingEnergy + this%angles%bendingEnergy(i)
    end do

    if(this%torsions%hasTorsion) then
       do i=1, this%torsions%numberOfTorsions
          EnergyUFF_instance%totalTorsionEnergy = EnergyUFF_instance%totalTorsionEnergy + this%torsions%torsionEnergy(i)
       end do
    end if

    if(this%vdwaals%VDW) then
       do i=1, this%vdwaals%numberOfVDWaals
          EnergyUFF_instance%totalVDWEnergy = EnergyUFF_instance%totalVDWEnergy + this%vdwaals%VDWEnergy(i)
       end do
    end if

    if(electrostatic) then
       if(this%electrostatic%isElectrostatic) then
          do i=1, this%electrostatic%numberOfElectrostatics
             EnergyUFF_instance%totalElectrostaticEnergy = EnergyUFF_instance%totalElectrostaticEnergy + this%electrostatic%electrostaticEnergy(i)
          end do
       end if
    end if

    if(this%inversions%hasInversions) then
       do i=1, this%inversions%numberOfInversions
          EnergyUFF_instance%totalInversionEnergy = EnergyUFF_instance%totalInversionEnergy + this%inversions%inversionEnergy(i)
       end do
    end if


    EnergyUFF_instance%totalStretchingEnergyKJ = EnergyUFF_instance%totalStretchingEnergy*4.1868
    EnergyUFF_instance%totalBendingEnergyKJ = EnergyUFF_instance%totalBendingEnergy*4.1868
    EnergyUFF_instance%totalTorsionEnergyKJ = EnergyUFF_instance%totalTorsionEnergy*4.1868
    EnergyUFF_instance%totalVDWEnergyKJ = EnergyUFF_instance%totalVDWEnergy*4.1868
    EnergyUFF_instance%totalInversionEnergyKJ = EnergyUFF_instance%totalInversionEnergy*4.1868
    EnergyUFF_instance%totalElectrostaticEnergyKJ = EnergyUFF_instance%totalElectrostaticEnergy*4.1868

    EnergyUFF_instance%totalEnergy = EnergyUFF_instance%totalStretchingEnergy + EnergyUFF_instance%totalBendingEnergy + &
         EnergyUFF_instance%totalTorsionEnergy + EnergyUFF_instance%totalVDWEnergy + EnergyUFF_instance%totalInversionEnergy + &
         EnergyUFF_instance%totalElectrostaticEnergy
    EnergyUFF_instance%totalEnergyKJ = EnergyUFF_instance%totalEnergy*4.1868
   
  end subroutine EnergyUFF_run

  !>
  !! @brief This routine calculates the total energies with the UFF
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-09-30
  !! @param [in] this Class with all information about the Energies (System)
  !! @param [in] graphs Class with all information about the Graph (System)
  !! @param [in] electrostatic LOGICAL evaluates if the user requires Electrostatic Energy  
  subroutine EnergyUFF_show(this, graphs, electrostatic)
    implicit none
    type(EnergyUFF) :: this
    type(Graph), intent(in) :: graphs
    logical, intent(in) :: electrostatic

    write(*,"(T5,A)") ""
    write(*,"(T5,A)") ""
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T12,A)") "Summary of total energies"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T7,A,T25,A,T39,A)") "Type", "kcal/mol", "kJ/mol"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Stretching", this%totalStretchingEnergy, this%totalStretchingEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Bending", this%totalBendingEnergy, this%totalBendingEnergyKJ
    if(graphs%torsions%hasTorsion) then
       write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Torsional", this%totalTorsionEnergy, this%totalTorsionEnergyKJ
    end if
    if(graphs%vdwaals%VDW) then
       write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Van der Waals", this%totalVDWEnergy, this%totalVDWEnergyKJ
    end if
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Out of Plane", this%totalInversionEnergy, this%totalInversionEnergyKJ
    if(electrostatic) then
       if(graphs%electrostatic%isElectrostatic) then
          write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Electrostatic", this%totalElectrostaticEnergy, this%totalElectrostaticEnergyKJ
       end if
    end if
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "TOTAL ENERGY", this%totalEnergy, this%totalEnergyKJ
    write(*,"(T5,A)") "------------------------------------------"

  end subroutine EnergyUFF_show

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine EnergyUFF_exception( typeMessage, description, debugDescription)
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
    
  end subroutine EnergyUFF_exception
  
end module EnergyUFF_
