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
!!        This module call the graph constructor and energy module
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
module MMFunctions_
  use CONTROL_
  use ParticleManager_
  use Matrix_
  use Graph_
  use EnergyUFF_
  use Exception_
  implicit none


	type :: MolecularMechanics
    
       	character(50) :: ffmethod
        logical :: electrostaticEnergy
        logical :: isInstanced

	end type MolecularMechanics

       type(MolecularMechanics), target :: MolecularMechanics_instance

       public :: &
            MolecularMechanics_constructor, &
            MolecularMechanics_destructor, &
            MolecularMechanics_run, &
            MolecularMechanics_show

contains
  !>
  !! @brief Defines the class constructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine MolecularMechanics_constructor()
    implicit none
  
    MolecularMechanics_instance%isInstanced =.true.

  end subroutine MolecularMechanics_constructor

  !>
  !! @brief Defines the class destructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine MolecularMechanics_destructor()
    implicit none

    MolecularMechanics_instance%isInstanced =.false.

  end subroutine MolecularMechanics_destructor

  !>
  !! @brief Initialize all Molecular Mechanics calculation
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] ffmethod CHARACTER Force Field selected by the user, for now only UFF has been implemented
  !! @param [in] electrostaticEnergy LOGICAL evaluates if the user requires Electrostatic Energy
  !! @param [in] printAllMM LOGICAL evaluate if the user requires aditional print information
  subroutine MolecularMechanics_run( ffmethod, electrostaticEnergy, printAllMM )
    implicit none
    character(50), intent(in) :: ffmethod
    logical, intent(in) :: electrostaticEnergy
    logical, intent(in) :: printAllMM
    type(Exception) :: ex
    integer :: i, j

    write(*,"(A)")"----------------------------------------------------------------------"
    write(*,"(A)")"** PROGRAM: MM (Molecular Mechanics).      Author: J.M Rodas  "
    write(*,"(A)")"----------------------------------------------------------------------"

    write(*,"(A)") "INFO: RUNNING WITH "//trim(ffmethod)//" FORCE FIELD."
    write(*,"(A)")" "

    MolecularMechanics_instance%ffmethod = ffmethod
    MolecularMechanics_instance%electrostaticEnergy = electrostaticEnergy
    if ( MolecularMechanics_instance%isInstanced ) then
       !! If force field = UFF initialize the graph
       if ( MolecularMechanics_instance%ffmethod == "UFF" ) then
          call Graph_initialize(MolecularMechanics_instance%ffmethod,MolecularMechanics_instance%electrostaticEnergy)
          !! Print all results with UFF
          if (printAllMM) then
             call MolecularMechanics_show(Graph_instance)
          end if
          !! Calculate total energies with UFF
          call EnergyUFF_run(Graph_instance,MolecularMechanics_instance%electrostaticEnergy)
          call EnergyUFF_show(EnergyUFF_instance, Graph_instance, MolecularMechanics_instance%electrostaticEnergy)
       else
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object MolecularMechanics in run() function" )
          call Exception_setDescription( ex, "This Force Field hasn't been implemented" )
          call Exception_show( ex )
       end if
    else
       call Exception_constructor( ex , ERROR )
       call Exception_setDebugDescription( ex, "Class object MolecularMechanics in run() function" )
       call Exception_setDescription( ex, "You should to instance MolecularMechanics module before use this function" )
       call Exception_show( ex )
    end if

  end subroutine MolecularMechanics_run

  !>
  !! @brief Show all information about Molecular Mechanics calculation,
  !! If the user needs aditional information about the Molecular Mechanics calculations, it must be activated in the input like this:
  !! <BLOCKQUOTE>
  !! CONTROL \n
  !! &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;printMM = T \n
  !! END CONTROL \n
  !! </BLOCKQUOTE> 
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-09-23
  !! @param [in] this Class with the information of the molecular graph
  subroutine MolecularMechanics_show(this)
    implicit none
    type(Graph), intent(in) :: this
    integer :: atomAIdx, AtomBIdx, AtomCIdx, AtomDIdx
    character(10) :: atomA, AtomB, AtomC, AtomD
    integer :: i

          write(*,"(T5,A)") ""
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          write(*,"(T30,A)") "INITIAL GEOMETRY: ANGSTROM"
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          write (*,"(T5,A5,T14,A,T20,A,T30,A,T45,A,T60,A,T75,A)") "Idx", &
               "Atom", "Type", "Charge(Z)", &
               "<x>","<y>","<z>"
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          do i=1,this%vertex%numberOfVertices
                write(*,"(T5,I5,T15,A,T20,A,T28,F8.2,T40,F10.5,T55,F10.5,T70,F10.5)") i, &
                     trim(this%vertex%symbol(i)), &
                     trim( this%vertex%type(i) ), &
                     this%vertex%charges(i), &
                     this%vertex%cartesianMatrix%values(i,1), &
                     this%vertex%cartesianMatrix%values(i,2), &
                     this%vertex%cartesianMatrix%values(i,3)
          end do
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"

          
          write(*,"(T5,A)") ""
          write(*,"(T5,A)") ""
          write(*,"(T45,A)") "STRETCHING ENERGY"
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          write(*,"(T15,A)") "Bond"
          write(*,"(T5,A,T35,A,T50,A,T65,A,T80,A,T100,A)") "----------------------------", "Bond order", "Bond length", &
               "Ideal length", "Force constant", "Energy"
          write(*,"(T5,A5,T15,A,T25,A,T50,A,T66,A,T80,A,T99,A)") "Idx", "atom A", &
               "atom B", "(Amstrong)", "(Amstrong)", "(kcal/mol*A^2)", "(kJ/mol)"
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          do i=1,this%edges%numberOfEdges
             atomAIdx=this%edges%connectionMatrix%values(i,1)
             Write( atomA, '(i10)' ) atomAIdx
             atomA = adjustl(trim(atomA))
             atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
             atomBIdx=this%edges%connectionMatrix%values(i,2)
             Write( atomB, '(i10)' ) atomBIdx
             atomB = adjustl(trim(atomB))
             atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
             write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T50,F8.5,T65,F8.5,T80,F12.5,T95,F12.5)") i, atomA, atomB, &
                  this%edges%bondOrder(i), &
                  this%edges%distance(i), &
                  this%edges%idealDistance(i), &
                  this%edges%forceConstant(i), &
                  this%edges%stretchingEnergyKJ(i)
          end do
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          write(*,"(T5,A)") ""

          write(*,"(T5,A)") ""
          write(*,"(T5,A)") ""
          write(*,"(T47,A)") "BENDING ENERGY"
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          write(*,"(T22,A)") "Angle"
          write(*,"(T5,A,T48,A,T60,A,T75,A,T95,A)") "--------------------------------------", &
               "Angle", &
               "Ideal Angle", "Force constant", "Energy"
          write(*,"(T5,A5,T15,A,T25,A,T35,A,T46,A,T61,A,T75,A,T94,A)") "Idx", "atom A", &
               "atom B", "atom C", "(Degrees)", "(Degrees)", "(kcal/mol*rad^2)", "(kJ/mol)"
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          do i=1,this%angles%numberOfAngles
             atomAIdx=this%angles%connectionMatrix%values(i,1)
             Write( atomA, '(i10)' ) atomAIdx
             atomA = adjustl(trim(atomA))
             atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
             atomBIdx=this%angles%connectionMatrix%values(i,2)
             Write( atomB, '(i10)' ) atomBIdx
             atomB = adjustl(trim(atomB))
             atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
             atomCIdx=this%angles%connectionMatrix%values(i,3)
             Write( atomC, '(i10)' ) atomCIdx
             atomC = adjustl(trim(atomC))
             atomC=trim(this%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
             write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,F10.5,T60,F10.5,T75,F12.5,T91,F12.5)") i, atomA, &
                  atomB, &
                  atomC, &
                  this%angles%theta(i), &
                  this%angles%idealTheta(i), &
                  this%angles%forceConstant(i), &
                  this%angles%bendingEnergyKJ(i)
          end do
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          write(*,"(T5,A)") ""

          if(this%torsions%hasTorsion) then
             write(*,"(T5,A)") ""
             write(*,"(T5,A)") ""
             write(*,"(T47,A)") "TORSIONAL ENERGY"
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             write(*,"(T22,A)") "Dihedral Angle"
             write(*,"(T5,A,T58,A,T70,A,T85,A,T105,A)") "------------------------------------------------", &
                  "Angle", &
                  "Ideal Angle", "Force constant", "Energy"
             write(*,"(T5,A5,T15,A,T25,A,T35,A,T45,A,T56,A,T71,A,T85,A,T104,A)") "Idx", "atom A", &
                  "atom B", "atom C", "atom D", "(Degrees)", "(Degrees)", "(kcal/mol*rad^2)", "(kJ/mol)"
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             do i=1,this%torsions%numberOfTorsions
                atomAIdx=this%torsions%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=this%torsions%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                atomCIdx=this%torsions%connectionMatrix%values(i,3)
                Write( atomC, '(i10)' ) atomCIdx
                atomC = adjustl(trim(atomC))
                atomC=trim(this%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
                atomDIdx=this%torsions%connectionMatrix%values(i,4)
                Write( atomD, '(i10)' ) atomDIdx
                atomD = adjustl(trim(atomD))
                atomD=trim(this%vertex%symbol(atomDIdx))//"("//trim(atomD)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,A,T55,F10.5,T70,F10.5,T85,F12.5,T100,F12.5)") i, atomA, &
                     atomB, &
                     atomC, &
                     atomD, &
                     this%torsions%phi(i), &
                     this%torsions%idealPhi(i), &
                     this%torsions%rotationalBarrier(i), &
                     this%torsions%torsionEnergyKJ(i)
             end do
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"
             write(*,"(T5,A)") ""
          end if

          if(this%inversions%hasInversions) then
             write(*,"(T5,A)") ""
             write(*,"(T5,A)") ""
             write(*,"(T32,A)") "INVERSION ENERGY (OUT OF PLANE)"
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"
             write(*,"(T22,A)") "Improper Angle"
             write(*,"(T5,A,T58,A,T70,A,T90,A)") "------------------------------------------------", &
                  "Angle", &
                  "Force constant", "Energy"
             write(*,"(T5,A5,T15,A,T25,A,T35,A,T45,A,T56,A,T70,A,T89,A)") "Idx", "atom A", &
                  "atom B", "atom C", "atom D", "(Degrees)", "(kcal/mol)", "(kJ/mol)"
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"
             do i=1,this%inversions%numberOfInversions
                atomAIdx=this%inversions%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=this%inversions%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                atomCIdx=this%inversions%connectionMatrix%values(i,3)
                Write( atomC, '(i10)' ) atomCIdx
                atomC = adjustl(trim(atomC))
                atomC=trim(this%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
                atomDIdx=this%inversions%connectionMatrix%values(i,4)
                Write( atomD, '(i10)' ) atomDIdx
                atomD = adjustl(trim(atomD))
                atomD=trim(this%vertex%symbol(atomDIdx))//"("//trim(atomD)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,A,T55,F10.5,T70,F12.5,T85,F12.5)") i, atomA, &
                     atomB, &
                     atomC, &
                     atomD, &
                     this%inversions%omega(i), &
                     this%inversions%forceConstant(i), &
                     this%inversions%inversionEnergyKJ(i)
             end do
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"    
             write(*,"(T5,A)") ""
          end if

          if(this%vdwaals%VDW) then
             write(*,"(T5,A)") ""
             write(*,"(T5,A)") ""
             write(*,"(T45,A)") "VAN DER WAALS ENERGY"
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             write(*,"(T15,A)") "VDW distance"
             write(*,"(T5,A,T36,A,T49,A,T66,A,T85,A)") "----------------------------", "Distance", &
                  "Ideal distance", "Well Depth", "Energy"
             write(*,"(T5,A5,T15,A,T25,A,T35,A,T51,A,T66,A,T84,A)") "Idx", "atom A", &
                  "atom B", "(Amstrong)", "(Amstrong)", "(kcal/mol)", "(kJ/mol)"
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             do i=1,this%vdwaals%numberOfVDWaals
                atomAIdx=this%vdwaals%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=this%vdwaals%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T50,F8.5,T62,F12.5,T80,F12.5)") i, atomA, atomB, &
                     this%vdwaals%distance(i) , &
                     this%vdwaals%idealDistance(i), &
                     this%vdwaals%wellDepth(i), &
                     this%vdwaals%VDWEnergyKJ(i)
             end do
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             write(*,"(T5,A)") ""
          end if

          if(MolecularMechanics_instance%electrostaticEnergy) then
             if(this%electrostatic%isElectrostatic) then
                write(*,"(T5,A)") ""
                write(*,"(T5,A)") ""
                write(*,"(T22,A)") "ELECTROSTATIC ENERGY"
                write(*,"(T5,A)") "--------------------------------------------------------"
                write(*,"(T15,A)") "Charges center"
                write(*,"(T5,A,T36,A,T53,A)") "----------------------------", "Distance", "Energy"
                write(*,"(T5,A5,T15,A,T25,A,T35,A,T52,A)") "Idx", "atom A", &
                     "atom B", "(Amstrong)", "(kJ/mol)"
                write(*,"(T5,A)") "--------------------------------------------------------"
                do i=1,this%electrostatic%numberOfElectrostatics
                   atomAIdx=this%electrostatic%connectionMatrix%values(i,1)
                   Write( atomA, '(i10)' ) atomAIdx
                   atomA = adjustl(trim(atomA))
                   atomA=trim(this%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                   atomBIdx=this%electrostatic%connectionMatrix%values(i,2)
                   Write( atomB, '(i10)' ) atomBIdx
                   atomB = adjustl(trim(atomB))
                   atomB=trim(this%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                   write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T47,F12.5)") i, atomA, atomB, &
                        this%electrostatic%distance(i) , &
                        this%electrostatic%electrostaticEnergyKJ(i)
                end do
                write(*,"(T5,A)") "--------------------------------------------------------"
                write(*,"(T5,A)") ""
             end if

             write(*,"(T5,A)") ""
             write(*,"(T5,A)") "-----------------------------"
             write(*,"(T5,A)") "PARTIAL CHARGES: EQeq Method"
             write(*,"(T5,A)") "-----------------------------"
             write (*,"(T5,A5,T14,A,T20,A)") "Idx", "Atom", "Charge(Z)"
             write(*,"(T5,A)") "-----------------------------"
             do i=1,this%vertex%numberOfVertices
                write(*,"(T5,I5,T15,A,T20,F8.5)") i, &
                     trim(this%vertex%symbol(i)), &
                     this%electrostatic%partialCharge(i)
             end do
             write(*,"(T5,A)") "-----------------------------"
          end if


  end subroutine MolecularMechanics_show

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  subroutine MolecularMechanics_exception( typeMessage, description, debugDescription)
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

  end subroutine MolecularMechanics_exception

  
end module MMFunctions_
