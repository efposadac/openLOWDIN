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
module AtomTypeUFF_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Rings_
  use MMCommons_
  use MatrixInteger_
  use Vector_
  use AromaticityFinder_
  use RingFinder_
  use Exception_
  implicit none


       public :: &
            AtomTypeUFF_run

contains

  subroutine AtomTypeUFF_run(ffAtomType, ring)
    implicit none
    type(Rings), intent(in) :: ring
    integer :: i   
    integer :: numberOfCenterofOptimization
    character(10), allocatable :: labelOfCenters(:)
    real(8), allocatable :: chargesOfCenters(:)
    character(10), allocatable :: ffAtomType(:)
    integer :: connectivity
    real(8) :: angleAverage
    real(8) :: SP2SP3AngleCutoff
    real(8) :: SPSP2AngleCutoff
    real(8) :: AngleCutoff1
    real(8) :: AngleCutoff2
    logical :: isAromatic
    logical :: isOrganometallic
    integer :: numberOfColumns
    type(Exception) :: ex


    SP2SP3AngleCutoff = 115.00000000
    SPSP2AngleCutoff = 160.00000000
    AngleCutoff1 = 100.00000000
    AngleCutoff2 = 105.00000000

    call MMCommons_constructor( MolecularSystem_instance )

    numberOfCenterofOptimization = ParticleManager_getNumberOfCentersOfOptimization()
  
    allocate( labelOfCenters( numberOfCenterofOptimization ) )
    labelOfCenters = ParticleManager_getLabelsOfCentersOfOptimization()

    allocate( chargesOfCenters( numberOfCenterofOptimization ) )
    chargesOfCenters = ParticleManager_getChargesOfCentersOfOptimization()

    allocate( ffAtomType( numberOfCenterofOptimization ) )

    do i=1, size(labelOfCenters)
!!******************************************************************************
!! Se evaluan los Hidrogenos
!!******************************************************************************
       if( trim( labelOfCenters(i) ) == "H" ) then
          ffAtomType(i) = "H_"
!!******************************************************************************
!! Se evaluan los carbonos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "C" ) then

          !! Se chequea la conectividad del carbono
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          !! Imprime el valor de conectividad, borrar luego
          ! write (*,"(T20,A,I)") "Conectividad: ", connectivity
          !! Si connectivity >= 4 se asume hibridacion sp3 (C_3)
          if ( connectivity >= 4 ) then
             ffAtomType(i) = "C_3"
          
          !! Si connectivity == 3 hay tres posibles carbonos: hibridacion sp3: (C_3); hibridacion sp2: (C_2); 
          !! hibridacion sp2 y aromatico: (C_R). 
          !! Por eso es necesario ademas chequear aromaticidad y angulos de enlace
          else if ( connectivity == 3 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Imprime el angulo promedio, borrar despues
             ! write (*,"(T20,A,F16.5)") "Angulo: ", angleAverage
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (C_3)
             if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "C_3"
             !! Falta programar la aromaticidad   
             else if ( ring%hasRings ) then
                isAromatic = AromaticityFinder_isAromatic( ring, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "C_R"
                else
                   ffAtomType(i) = "C_2"
                end if
             else
                ffAtomType(i) = "C_2"
             end if

          else if ( connectivity == 2 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (C_3)
             if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "C_3"
             !! Si el angulo es menor a 160 entonces se asume hibridacion sp2 (C_2)
             else if ( angleAverage < SPSP2AngleCutoff ) then
                ffAtomType(i) = "C_2"
             else
                ffAtomType(i) = "C_1"
             end if
!!!! Falta implementar las otras opciones para conectividad == 1
          else
             ffAtomType(i) = "C_1"
          end if
!!******************************************************************************
!! Se evaluan los Nitrogenos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "N" ) then
          !! Se chequea la conectividad del nitrogeno
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )

          !! Si connectivity >= 4 se asume hibridacion sp3 (N_3)
          if ( connectivity >= 4 ) then
             ffAtomType(i) = "N_3"
          !! Si connectivity == 3 hay tres posibles nitrogenos: hibridacion sp3: (N_3); hibridacion sp2: (N_2); 
          !! hibridacion sp2 y aromatico: (N_R). 
          !! Por eso es necesario ademas chequear aromaticidad y angulos de enlace
          else if ( connectivity == 3 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Falta programar la aromaticidad   
             if ( ring%hasRings ) then
                isAromatic = AromaticityFinder_isAromatic( ring, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "N_R"
                else
                   ffAtomType(i) = "N_2"
                end if
             else if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "N_3"
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (N_3)
             else
                ffAtomType(i) = "N_2"
             end if
          else if ( connectivity == 2 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Falta programar la aromaticidad   
             if ( ring%hasRings ) then
                isAromatic = AromaticityFinder_isAromatic( ring, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "N_R"
                else if ( angleAverage < SPSP2AngleCutoff ) then
                   ffAtomType(i) = "N_2"
                else
                   ffAtomType(i) = "N_1"
                end if
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (N_3)
             else if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "N_3"
             !! Si el angulo es menor a 160 entonces se asume hibridacion sp2 (C_2)
             else if ( angleAverage < SPSP2AngleCutoff ) then
                ffAtomType(i) = "N_2"
             else
                ffAtomType(i) = "N_1"
             end if
!!!! Falta implementar las otras opciones para conectividad == 1
          else
             ffAtomType(i) = "N_1"
          end if
!!******************************************************************************
!! Se evaluan los Oxigenos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "O" ) then
          !! Se chequea la conectividad del oxigeno
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity >= 2 ) then
                !! Falta programar la aromaticidad   
             if ( ring%hasRings ) then
                isAromatic = AromaticityFinder_isAromatic( ring, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "O_R"
                else
                   ffAtomType(i) = "O_3"
                end if
             else
                ffAtomType(i) = "O_3"
             end if
          else
             ffAtomType(i) = "O_2"
          end if
!!******************************************************************************
!! Se evaluan los Fosforos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "P" ) then
          !! Se chequea la conectividad del oxigeno
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity >= 4 ) then
             isOrganometallic = MMCommons_isOrganometallic( MolecularSystem_instance, i, connectivity, labelOfCenters )
             if (isOrganometallic) then
                ffAtomType(i) = "P_3+q"
             else
                ffAtomType(i) = "P_3+5"
             end if
          else
             ffAtomType(i) = "P_3+3"
          end if
!!******************************************************************************
!! Se evaluan los Azufres
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "S" ) then
          !! Se chequea la conectividad del azufre
          angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i ) 
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "S_3+6"
          else if ( connectivity == 2 ) then
             if ( ring%hasRings ) then
                isAromatic = AromaticityFinder_isAromatic( ring, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "S_R"
                else if ( angleAverage <= AngleCutoff1 ) then
                   ffAtomType(i) = "S_3+2"
                else if ( angleAverage <= AngleCutoff2 ) then
                   ffAtomType(i) = "S_3+4"
                else if ( angleAverage <= SP2SP3AngleCutoff ) then
                   ffAtomType(i) = "S_3+6"
                else
                   ffAtomType(i) = "S_2"
                end if
             end if
          else if ( angleAverage <= AngleCutoff1 ) then
             ffAtomType(i) = "S_3+2"
          else if ( angleAverage <= AngleCutoff2 ) then
             ffAtomType(i) = "S_3+4"
          else if ( angleAverage <= SP2SP3AngleCutoff ) then
             ffAtomType(i) = "S_3+6"
          else
             ffAtomType(i) = "S_2"
          end if
!!******************************************************************************
!! Se evalua el grupo 17 (Halogenos)
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "F" ) then
          ffAtomType(i) = "F_"
       else if( trim( labelOfCenters(i) ) == "CL" ) then
          ffAtomType(i) = "Cl"
       else if( trim( labelOfCenters(i) ) == "BR" ) then
          ffAtomType(i) = "Br"
       else if( trim( labelOfCenters(i) ) == "I" ) then
          ffAtomType(i) = "I_"
       else if( trim( labelOfCenters(i) ) == "AT" ) then
          ffAtomType(i) = "At"
!!******************************************************************************
!! Se evalua el grupo 1
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "LI" ) then
          ffAtomType(i) = "Li"
       else if( trim( labelOfCenters(i) ) == "NA" ) then
          ffAtomType(i) = "Na"
       else if( trim( labelOfCenters(i) ) == "K" ) then
          ffAtomType(i) = "K_"
       else if( trim( labelOfCenters(i) ) == "RB" ) then
          ffAtomType(i) = "Rb"
       else if( trim( labelOfCenters(i) ) == "CS" ) then
          ffAtomType(i) = "Cs"
       else if( trim( labelOfCenters(i) ) == "FR" ) then
          ffAtomType(i) = "Fr"
!!******************************************************************************
!! Se evalua el grupo 2
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "BE" ) then
          ffAtomType(i) = "Be3+2"
       else if( trim( labelOfCenters(i) ) == "MG" ) then
          ffAtomType(i) = "Mg3+2"
       else if( trim( labelOfCenters(i) ) == "CA" ) then
          ffAtomType(i) = "Ca6+2"
       else if( trim( labelOfCenters(i) ) == "SR" ) then
          ffAtomType(i) = "Sr6+2"
       else if( trim( labelOfCenters(i) ) == "BA" ) then
          ffAtomType(i) = "Ba6+2"
       else if( trim( labelOfCenters(i) ) == "RA" ) then
          ffAtomType(i) = "Ra6+2"
!!******************************************************************************
!! Se evalua el grupo 13
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "B" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity >= 4 ) then
             ffAtomType(i) = "B_3"
          else
             ffAtomType(i) = "B_2"
          end if
       else if( trim( labelOfCenters(i) ) == "AL" ) then
          ffAtomType(i) = "Al3"
       else if( trim( labelOfCenters(i) ) == "GA" ) then
          ffAtomType(i) = "Ga3+3"
       else if( trim( labelOfCenters(i) ) == "IN" ) then
          ffAtomType(i) = "In3+3"
       else if( trim( labelOfCenters(i) ) == "TL" ) then
          ffAtomType(i) = "Tl3+3"
!!******************************************************************************
!! Se evalua el grupo 14
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "SI" ) then
          ffAtomType(i) = "Si3"
       else if( trim( labelOfCenters(i) ) == "GE" ) then
          ffAtomType(i) = "Ge3"
       else if( trim( labelOfCenters(i) ) == "SN" ) then
          ffAtomType(i) = "Sn3"
       else if( trim( labelOfCenters(i) ) == "PB" ) then
          ffAtomType(i) = "Pb3"
!!******************************************************************************
!! Se evalua el grupo 15
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "AS" ) then
          ffAtomType(i) = "As3+3"
       else if( trim( labelOfCenters(i) ) == "SB" ) then
          ffAtomType(i) = "Sb3+3"
       else if( trim( labelOfCenters(i) ) == "BI" ) then
          ffAtomType(i) = "Bi3+3"
!!******************************************************************************
!! Se evalua el grupo 16
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "SE" ) then
          ffAtomType(i) = "Se3+2"
       else if( trim( labelOfCenters(i) ) == "TE" ) then
          ffAtomType(i) = "Te3+2"
       else if( trim( labelOfCenters(i) ) == "PO" ) then
          ffAtomType(i) = "Po3+2"
!!******************************************************************************
!! Se evalua el grupo 18 (Gases Nobles)
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "HE" ) then
          ffAtomType(i) = "He4+4"
       else if( trim( labelOfCenters(i) ) == "NE" ) then
          ffAtomType(i) = "Ne4+4"
       else if( trim( labelOfCenters(i) ) == "AR" ) then
          ffAtomType(i) = "Ar4+4"
       else if( trim( labelOfCenters(i) ) == "KR" ) then
          ffAtomType(i) = "Kr4+4"
       else if( trim( labelOfCenters(i) ) == "XE" ) then
          ffAtomType(i) = "Xe4+4"
       else if( trim( labelOfCenters(i) ) == "RN" ) then
          ffAtomType(i) = "Rn4+4"
!!******************************************************************************
!! Se evalua la primera serie de transicion 
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "SC" ) then
          ffAtomType(i) = "Sc3+3"
       else if( trim( labelOfCenters(i) ) == "TI" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "Ti6+4"
          else
             ffAtomType(i) = "Ti3+4"
          end if
       else if( trim( labelOfCenters(i) ) == "V" ) then
          ffAtomType(i) = "V_3+5"
       else if( trim( labelOfCenters(i) ) == "CR" ) then
          ffAtomType(i) = "Cr6+3"
       else if( trim( labelOfCenters(i) ) == "MN" ) then
          ffAtomType(i) = "Mn6+2"
       else if( trim( labelOfCenters(i) ) == "FE" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "Fe6+2"
          else
             ffAtomType(i) = "Fe3+2"
          end if
       else if( trim( labelOfCenters(i) ) == "CO" ) then
          ffAtomType(i) = "Co6+3"
       else if( trim( labelOfCenters(i) ) == "NI" ) then
          ffAtomType(i) = "Ni4+2"
       else if( trim( labelOfCenters(i) ) == "CU" ) then
          ffAtomType(i) = "Cu3+1"
       else if( trim( labelOfCenters(i) ) == "ZN" ) then
          ffAtomType(i) = "Zn3+2"
!!******************************************************************************
!! Se evalua la segunda serie de transicion 
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "Y" ) then
          ffAtomType(i) = "Y_3+3"
       else if( trim( labelOfCenters(i) ) == "ZR" ) then
          ffAtomType(i) = "Zr3+4"
       else if( trim( labelOfCenters(i) ) == "NB" ) then
          ffAtomType(i) = "Nb3+5"
       else if( trim( labelOfCenters(i) ) == "MO" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "Mo6+6"
          else
             ffAtomType(i) = "Mo3+6"
          end if
       else if( trim( labelOfCenters(i) ) == "TC" ) then
          ffAtomType(i) = "Tc6+5"
       else if( trim( labelOfCenters(i) ) == "RU" ) then
          ffAtomType(i) = "Ru6+2"
       else if( trim( labelOfCenters(i) ) == "RH" ) then
          ffAtomType(i) = "Rh6+3"
       else if( trim( labelOfCenters(i) ) == "PD" ) then
          ffAtomType(i) = "Pd4+2"
       else if( trim( labelOfCenters(i) ) == "AG" ) then
          ffAtomType(i) = "Ag1+1"
       else if( trim( labelOfCenters(i) ) == "CD" ) then
          ffAtomType(i) = "Cd3+2"
!!******************************************************************************
!! Se evalua la tercera serie de transicion y lantanidos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "LA" ) then
          ffAtomType(i) = "La3+3"
       else if( trim( labelOfCenters(i) ) == "CE" ) then
          ffAtomType(i) = "Ce6+3"
       else if( trim( labelOfCenters(i) ) == "PR" ) then
          ffAtomType(i) = "Pr6+3"
       else if( trim( labelOfCenters(i) ) == "ND" ) then
          ffAtomType(i) = "Nd6+3"
       else if( trim( labelOfCenters(i) ) == "PM" ) then
          ffAtomType(i) = "Pm6+3"
       else if( trim( labelOfCenters(i) ) == "SM" ) then
          ffAtomType(i) = "Sm6+3"
       else if( trim( labelOfCenters(i) ) == "EU" ) then
          ffAtomType(i) = "Eu6+3"
       else if( trim( labelOfCenters(i) ) == "GD" ) then
          ffAtomType(i) = "Gd6+3"
       else if( trim( labelOfCenters(i) ) == "TB" ) then
          ffAtomType(i) = "Tb6+3"
       else if( trim( labelOfCenters(i) ) == "DY" ) then
          ffAtomType(i) = "Dy6+3"
       else if( trim( labelOfCenters(i) ) == "HO" ) then
          ffAtomType(i) = "Ho6+3"
       else if( trim( labelOfCenters(i) ) == "ER" ) then
          ffAtomType(i) = "Er6+3"
       else if( trim( labelOfCenters(i) ) == "TM" ) then
          ffAtomType(i) = "Tm6+3"
       else if( trim( labelOfCenters(i) ) == "YB" ) then
          ffAtomType(i) = "Yb6+3"
       else if( trim( labelOfCenters(i) ) == "LU" ) then
          ffAtomType(i) = "Lu6+3"
       else if( trim( labelOfCenters(i) ) == "HF" ) then
          ffAtomType(i) = "Hf3+4"
       else if( trim( labelOfCenters(i) ) == "TA" ) then
          ffAtomType(i) = "Ta3+5"
       else if( trim( labelOfCenters(i) ) == "W" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "W_6+6"
          else
             ffAtomType(i) = "W_3+6"
          end if
       else if( trim( labelOfCenters(i) ) == "RE" ) then
          connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
          if ( connectivity == 6 ) then
             ffAtomType(i) = "Re6+5"
          else
             ffAtomType(i) = "Re3+7"
          end if
       else if( trim( labelOfCenters(i) ) == "OS" ) then
          ffAtomType(i) = "Os6+6"
       else if( trim( labelOfCenters(i) ) == "IR" ) then
          ffAtomType(i) = "Ir6+3"
       else if( trim( labelOfCenters(i) ) == "PT" ) then
          ffAtomType(i) = "Pt4+2"
       else if( trim( labelOfCenters(i) ) == "AU" ) then
          ffAtomType(i) = "Au4+3"
       else if( trim( labelOfCenters(i) ) == "HG" ) then
          ffAtomType(i) = "Hg1+2"
!!******************************************************************************
!! Se evalua los Actinidos
!!******************************************************************************
       else if( trim( labelOfCenters(i) ) == "AC" ) then
          ffAtomType(i) = "Ac6+3"
       else if( trim( labelOfCenters(i) ) == "TH" ) then
          ffAtomType(i) = "Th6+4"
       else if( trim( labelOfCenters(i) ) == "PA" ) then
          ffAtomType(i) = "Pa6+4"
       else if( trim( labelOfCenters(i) ) == "U" ) then
          ffAtomType(i) = "U_6+4"
       else if( trim( labelOfCenters(i) ) == "NP" ) then
          ffAtomType(i) = "Np6+4"
       else if( trim( labelOfCenters(i) ) == "PU" ) then
          ffAtomType(i) = "Pu6+4"
       else if( trim( labelOfCenters(i) ) == "AM" ) then
          ffAtomType(i) = "Am6+4"
       else if( trim( labelOfCenters(i) ) == "CM" ) then
          ffAtomType(i) = "Cm6+3"
       else if( trim( labelOfCenters(i) ) == "BK" ) then
          ffAtomType(i) = "Bk6+3"
       else if( trim( labelOfCenters(i) ) == "CF" ) then
          ffAtomType(i) = "Cf6+3"
       else if( trim( labelOfCenters(i) ) == "ES" ) then
          ffAtomType(i) = "Es6+3"
       else if( trim( labelOfCenters(i) ) == "FM" ) then
          ffAtomType(i) = "Fm6+3"
       else if( trim( labelOfCenters(i) ) == "MD" ) then
          ffAtomType(i) = "Md6+3"
       else if( trim( labelOfCenters(i) ) == "NO" ) then
          ffAtomType(i) = "No6+3"
       else if( trim( labelOfCenters(i) ) == "LR" ) then
          ffAtomType(i) = "Lr6+3"
!!******************************************************************************
       else
          call Exception_constructor( ex , ERROR )
          call Exception_setDebugDescription( ex, "Class object AtomTypeUFF in run() function" )
          call Exception_setDescription( ex, "This Atom type hasn't been implemented in UFF Force Field" )
          call Exception_show( ex )
       end if
    end do
    
  end subroutine AtomTypeUFF_run

  subroutine AtomTypeUFF_exception( typeMessage, description, debugDescription)
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

  end subroutine AtomTypeUFF_exception

end module AtomTypeUFF_
