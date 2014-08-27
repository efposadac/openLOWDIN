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
!! @brief Moller-Plesset and APMO-Moller-Plesset program.
!!        This module allows to make calculations in the APMO-Moller-Plesset framework
!! @author  J.M. Rodas, E. F. Posada and S. A. Gonzalez.
!!
!! <b> Creation date : </b> 2013-10-03
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos basicos para correccion de segundo orden
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin 1
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
!!        -# Rewrite the module as a program and adapts to Lowdin 2
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module MMFunctions_
  use CONTROL_
  use Graph_
  use EnergyUFF_
  use Exception_
  implicit none


	type :: MolecularMechanics
    
       	character(50) :: ffmethod
        logical :: isInstanced

	end type MolecularMechanics

       type(MolecularMechanics), target :: MolecularMechanics_instance

       public :: &
            MolecularMechanics_constructor, &
            MolecularMechanics_destructor, &
            MolecularMechanics_run

contains
	!**
	! Define el constructor para la clase
	!
	!**
  subroutine MolecularMechanics_constructor()
    implicit none
  
    MolecularMechanics_instance%isInstanced =.true.

  end subroutine MolecularMechanics_constructor

	!**
	! Define el destructor para clase
	!
	!**
  subroutine MolecularMechanics_destructor()
    implicit none

    MolecularMechanics_instance%isInstanced =.false.

  end subroutine MolecularMechanics_destructor

  subroutine MolecularMechanics_run( ffmethod )
    implicit none
    character(50), intent(in) :: ffmethod
    type(Exception) :: ex
!! Parametros para impresion borrar luego
    integer :: atomAIdx, AtomBIdx, AtomCIdx, AtomDIdx
    character(10) :: atomA, AtomB, AtomC, AtomD
    integer :: i

    
    MolecularMechanics_instance%ffmethod = ffmethod

    if ( MolecularMechanics_instance%isInstanced ) then
    
       if ( MolecularMechanics_instance%ffmethod == "UFF" ) then
          call Graph_initialize(MolecularMechanics_instance%ffmethod)


          write(*,"(T5,A)") ""
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          write(*,"(T30,A)") "INITIAL GEOMETRY: AMSTRONG"
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          write (*,"(T5,A5,T14,A,T20,A,T30,A,T45,A,T60,A,T75,A)") "Idx", &
               "Atom", "Type", "Charge(Z)", &
               "<x>","<y>","<z>"
          write(*,"(T5,A)") "-----------------------------------------------------------------------------"
          do i=1,Graph_instance%vertex%numberOfVertices
                write(*,"(T5,I5,T15,A,T20,A,T28,F8.2,T40,F10.5,T55,F10.5,T70,F10.5)") i, &
                     trim(Graph_instance%vertex%symbol(i)), &
                     trim( Graph_instance%vertex%type(i) ), &
                     Graph_instance%vertex%charges(i), &
                     Graph_instance%vertex%cartesianMatrix%values(i,1), &
                     Graph_instance%vertex%cartesianMatrix%values(i,2), &
                     Graph_instance%vertex%cartesianMatrix%values(i,3)
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
          do i=1,Graph_instance%edges%numberOfEdges
             atomAIdx=Graph_instance%edges%connectionMatrix%values(i,1)
             Write( atomA, '(i10)' ) atomAIdx
             atomA = adjustl(trim(atomA))
             atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
             atomBIdx=Graph_instance%edges%connectionMatrix%values(i,2)
             Write( atomB, '(i10)' ) atomBIdx
             atomB = adjustl(trim(atomB))
             atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
             write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T50,F8.5,T65,F8.5,T80,F12.5,T95,F12.5)") i, atomA, atomB, &
                  Graph_instance%edges%bondOrder(i), &
                  Graph_instance%edges%distance(i), &
                  Graph_instance%edges%idealDistance(i), &
                  Graph_instance%edges%forceConstant(i), &
                  ! Graph_instance%edges%stretchingEnergy(i), &
                  Graph_instance%edges%stretchingEnergyKJ(i)
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
          do i=1,Graph_instance%angles%numberOfAngles
             atomAIdx=Graph_instance%angles%connectionMatrix%values(i,1)
             Write( atomA, '(i10)' ) atomAIdx
             atomA = adjustl(trim(atomA))
             atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
             atomBIdx=Graph_instance%angles%connectionMatrix%values(i,2)
             Write( atomB, '(i10)' ) atomBIdx
             atomB = adjustl(trim(atomB))
             atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
             atomCIdx=Graph_instance%angles%connectionMatrix%values(i,3)
             Write( atomC, '(i10)' ) atomCIdx
             atomC = adjustl(trim(atomC))
             atomC=trim(Graph_instance%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
             write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,F10.5,T60,F10.5,T75,F12.5,T91,F12.5)") i, atomA, &
                  atomB, &
                  atomC, &
                  Graph_instance%angles%theta(i), &
                  Graph_instance%angles%idealTheta(i), &
                  Graph_instance%angles%forceConstant(i), &
                  Graph_instance%angles%bendingEnergyKJ(i)
          end do
          write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
          write(*,"(T5,A)") ""

          if(Graph_instance%torsions%hasTorsion) then
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
             do i=1,Graph_instance%torsions%numberOfTorsions
                atomAIdx=Graph_instance%torsions%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=Graph_instance%torsions%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                atomCIdx=Graph_instance%torsions%connectionMatrix%values(i,3)
                Write( atomC, '(i10)' ) atomCIdx
                atomC = adjustl(trim(atomC))
                atomC=trim(Graph_instance%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
                atomDIdx=Graph_instance%torsions%connectionMatrix%values(i,4)
                Write( atomD, '(i10)' ) atomDIdx
                atomD = adjustl(trim(atomD))
                atomD=trim(Graph_instance%vertex%symbol(atomDIdx))//"("//trim(atomD)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,A,T55,F10.5,T70,F10.5,T85,F12.5,T100,F12.5)") i, atomA, &
                     atomB, &
                     atomC, &
                     atomD, &
                     Graph_instance%torsions%phi(i), &
                     Graph_instance%torsions%idealPhi(i), &
                     Graph_instance%torsions%rotationalBarrier(i), &
                     Graph_instance%torsions%torsionEnergyKJ(i)
             end do
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"
             write(*,"(T5,A)") ""
          end if

          if(Graph_instance%inversions%hasInversions) then
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
             do i=1,Graph_instance%inversions%numberOfInversions
                atomAIdx=Graph_instance%inversions%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=Graph_instance%inversions%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                atomCIdx=Graph_instance%inversions%connectionMatrix%values(i,3)
                Write( atomC, '(i10)' ) atomCIdx
                atomC = adjustl(trim(atomC))
                atomC=trim(Graph_instance%vertex%symbol(atomCIdx))//"("//trim(atomC)//")"
                atomDIdx=Graph_instance%inversions%connectionMatrix%values(i,4)
                Write( atomD, '(i10)' ) atomDIdx
                atomD = adjustl(trim(atomD))
                atomD=trim(Graph_instance%vertex%symbol(atomDIdx))//"("//trim(atomD)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,A,T45,A,T55,F10.5,T70,F12.5,T85,F12.5)") i, atomA, &
                     atomB, &
                     atomC, &
                     atomD, &
                     Graph_instance%inversions%omega(i), &
                     Graph_instance%inversions%forceConstant(i), &
                     Graph_instance%inversions%inversionEnergyKJ(i)
             end do
             write(*,"(T5,A)") "----------------------------------------------------------------------------------------------"    
             write(*,"(T5,A)") ""
          end if

          if(Graph_instance%vdwaals%VDW) then
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
             do i=1,Graph_instance%vdwaals%numberOfVDWaals
                atomAIdx=Graph_instance%vdwaals%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=Graph_instance%vdwaals%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T50,F8.5,T62,F12.5,T80,F12.5)") i, atomA, atomB, &
                     Graph_instance%vdwaals%distance(i) , &
                     Graph_instance%vdwaals%idealDistance(i), &
                     Graph_instance%vdwaals%wellDepth(i), &
                     Graph_instance%vdwaals%VDWEnergyKJ(i)
             end do
             write(*,"(T5,A)") "-------------------------------------------------------------------------------------------------------"
             write(*,"(T5,A)") ""
          end if

          if(Graph_instance%electrostatic%isElectrostatic) then
             write(*,"(T5,A)") ""
             write(*,"(T5,A)") ""
             write(*,"(T22,A)") "ELECTROSTATIC ENERGY"
             write(*,"(T5,A)") "--------------------------------------------------------"
             write(*,"(T15,A)") "Charges center"
             write(*,"(T5,A,T36,A,T53,A)") "----------------------------", "Distance", "Energy"
             write(*,"(T5,A5,T15,A,T25,A,T35,A,T52,A)") "Idx", "atom A", &
                  "atom B", "(Amstrong)", "(kJ/mol)"
             write(*,"(T5,A)") "--------------------------------------------------------"
             do i=1,Graph_instance%electrostatic%numberOfElectrostatics
                atomAIdx=Graph_instance%electrostatic%connectionMatrix%values(i,1)
                Write( atomA, '(i10)' ) atomAIdx
                atomA = adjustl(trim(atomA))
                atomA=trim(Graph_instance%vertex%symbol(atomAIdx))//"("//trim(atomA)//")"
                atomBIdx=Graph_instance%electrostatic%connectionMatrix%values(i,2)
                Write( atomB, '(i10)' ) atomBIdx
                atomB = adjustl(trim(atomB))
                atomB=trim(Graph_instance%vertex%symbol(atomBIdx))//"("//trim(atomB)//")"
                write(*,"(T5,I5,T15,A,T25,A,T35,F8.5,T47,F12.5)") i, atomA, atomB, &
                     Graph_instance%electrostatic%distance(i) , &
                     Graph_instance%electrostatic%electrostaticEnergyKJ(i)
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
          do i=1,Graph_instance%vertex%numberOfVertices
                write(*,"(T5,I5,T15,A,T20,F8.5)") i, &
                     trim(Graph_instance%vertex%symbol(i)), &
                     Graph_instance%electrostatic%partialCharge(i)
          end do
          write(*,"(T5,A)") "-----------------------------"


          call EnergyUFF_run(Graph_instance)

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
