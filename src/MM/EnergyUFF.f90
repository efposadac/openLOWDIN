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
  
  
  public :: &
       EnergyUFF_run
       ! EnergyUFF_show

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
    integer :: i


    totalStretchingEnergy = 0.0
    totalBendingEnergy = 0.0
    totalTorsionEnergy = 0.0
    totalVDWEnergy = 0.0
    totalInversionEnergy = 0.0
    totalElectrostaticEnergy = 0.0

    do i=1, this%edges%numberOfEdges
       totalStretchingEnergy = totalStretchingEnergy + this%edges%stretchingEnergy(i)
    end do

    do i=1, this%angles%numberOfAngles
       totalBendingEnergy = totalBendingEnergy + this%angles%bendingEnergy(i)
    end do

    if(this%torsions%hasTorsion) then
       do i=1, this%torsions%numberOfTorsions
          totalTorsionEnergy = totalTorsionEnergy + this%torsions%torsionEnergy(i)
       end do
    end if

    if(this%vdwaals%VDW) then
       do i=1, this%vdwaals%numberOfVDWaals
          totalVDWEnergy = totalVDWEnergy + this%vdwaals%VDWEnergy(i)
       end do
    end if

    if(electrostatic) then
       if(this%electrostatic%isElectrostatic) then
          do i=1, this%electrostatic%numberOfElectrostatics
             totalElectrostaticEnergy = totalElectrostaticEnergy + this%electrostatic%electrostaticEnergy(i)
          end do
       end if
    end if

    if(this%inversions%hasInversions) then
       do i=1, this%inversions%numberOfInversions
          totalInversionEnergy = totalInversionEnergy + this%inversions%inversionEnergy(i)
       end do
    end if


    totalStretchingEnergyKJ = totalStretchingEnergy*4.1868
    totalBendingEnergyKJ = totalBendingEnergy*4.1868
    totalTorsionEnergyKJ = totalTorsionEnergy*4.1868
    totalVDWEnergyKJ = totalVDWEnergy*4.1868
    totalInversionEnergyKJ = totalInversionEnergy*4.1868
    totalElectrostaticEnergyKJ = totalElectrostaticEnergy*4.1868

    totalEnergy = totalStretchingEnergy + totalBendingEnergy + totalTorsionEnergy + totalVDWEnergy + totalInversionEnergy + totalElectrostaticEnergy
    totalEnergyKJ = totalEnergy*4.1868

    write(*,"(T5,A)") ""
    write(*,"(T5,A)") ""
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T12,A)") "Summary of total energies"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T7,A,T25,A,T39,A)") "Type", "kcal/mol", "kJ/mol"
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Stretching", totalStretchingEnergy, totalStretchingEnergyKJ
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Bending", totalBendingEnergy, totalBendingEnergyKJ
    if(this%torsions%hasTorsion) then
       write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Torsional", totalTorsionEnergy, totalTorsionEnergyKJ
    end if
    if(this%vdwaals%VDW) then
       write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Van der Waals", totalVDWEnergy, totalVDWEnergyKJ
    end if
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Out of Plane", totalInversionEnergy, totalInversionEnergyKJ
    if(electrostatic) then
       if(this%electrostatic%isElectrostatic) then
          write(*,"(T5,A,T20,F12.5,T34,F12.5)") "Electrostatic", totalElectrostaticEnergy, totalElectrostaticEnergyKJ
       end if
    end if
    write(*,"(T5,A)") "------------------------------------------"
    write(*,"(T5,A,T20,F12.5,T34,F12.5)") "TOTAL", totalEnergy, totalEnergyKJ
    write(*,"(T5,A)") "------------------------------------------"

    
  end subroutine EnergyUFF_run
  
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
