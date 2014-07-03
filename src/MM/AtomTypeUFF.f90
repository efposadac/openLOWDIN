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
  use MMCommons_
  use MatrixInteger_
  use Vector_
  use AromaticityFinder_
  use RingFinder_
  use Exception_
  implicit none


       !  type :: AtomTypeUFF
    
       ! 	character(50) :: ffmethod
       !  logical :: isInstanced

       !  end type AtomTypeUFF

       ! type(AtomTypeUFF), target :: AtomTypeUFF_instance
!       character(50) :: job
       ! private :: &
            ! AtomTypeUFF_ffmethod

       public :: &
            ! AtomTypeUFF_constructor, &
            ! AtomTypeUFF_destructor, &
            ! AtomTypeUFF_show, &
            AtomTypeUFF_run
            ! AtomTypeUFF_getTotalEnergy, &
            ! AtomTypeUFF_getEnergyCorrection, &
            ! AtomTypeUFF_getSpecieCorrection


contains
	!**
	! Define el constructor para la clase
	!
	!**
  ! subroutine AtomTypeUFF_constructor()
  !   implicit none
  
  !   AtomTypeUFF_instance%isInstanced =.true.

  ! end subroutine AtomTypeUFF_constructor

	!**
	! Define el destructor para clase
	!
	!**
  ! subroutine AtomTypeUFF_destructor()
  !   implicit none

  !   AtomTypeUFF_instance%isInstanced =.false.

  ! end subroutine AtomTypeUFF_destructor

  subroutine AtomTypeUFF_run()
    implicit none
    integer :: i   
    integer :: numberOfCenterofOptimization
    character(10), allocatable :: labelOfCenters(:)
    character(10), allocatable :: ffAtomType(:)
    integer :: connectivity
    real(8) :: angleAverage
    real(8) :: SP2SP3AngleCutoff
    real(8) :: SPSP2AngleCutoff
    logical :: isAromatic
    integer :: numberOfEdges
    integer :: cyclomaticNumber
    type(MatrixInteger), allocatable :: edges(:)
    type(MatrixInteger) :: connectivityMatrix
    type(Vector) :: bonds
    type(MatrixInteger), allocatable :: rings(:)
    integer :: numberOfRings
    integer :: numberOfColumns

    SP2SP3AngleCutoff = 115.00000000
    SPSP2AngleCutoff = 160.00000000

    call MMCommons_constructor( MolecularSystem_instance )

    numberOfCenterofOptimization = ParticleManager_getNumberOfCentersOfOptimization()
    
!!******************************************************************************
!! Se calcula el cyclomatic number el cual es aquivalente al numero de anillos
!! L. Matyska, J. Comp. Chem. 9(5), 455 (1988)
!! si cyclomaticNumber = 0 no hay anillos 
!!******************************************************************************
    numberOfEdges=size(MolecularSystem_instance%intCoordinates%distanceBondValue%values)
    cyclomaticNumber = numberOfEdges - numberOfCenterofOptimization + 2
    
    allocate( labelOfCenters( numberOfCenterofOptimization ) )
    labelOfCenters = ParticleManager_getLabelsOfCentersOfOptimization()

    allocate( ffAtomType( numberOfCenterofOptimization ) )

    if (cyclomaticNumber>=1) then
       call MMCommons_pruningGraph( MolecularSystem_instance, numberOfCenterofOptimization, edges, connectivityMatrix, bonds )
       call RingFinder_getRings( edges, connectivityMatrix, cyclomaticNumber, rings )
       numberOfRings=size(rings)
       write (*,"(T10,A)") " Rings "
       write (*,"(T10,A)") "--------------------------------------------"
       do i=1,numberOfRings
          numberOfColumns = size(rings(i)%values)
          write (*,"(T10,<numberOfColumns>I)") rings(i)%values(1,:)
       end do
    end if

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

          !! Si connectivity >= 4 se asume hibridacion sp3 (C_3)
          if ( connectivity >= 4 ) then
             ffAtomType(i) = "C_3"
          
          !! Si connectivity == 3 hay tres posibles carbonos: hibridacion sp3: (C_3); hibridacion sp2: (C_2); 
          !! hibridacion sp2 y aromatico: (C_R). 
          !! Por eso es necesario ademas chequear aromaticidad y angulos de enlace
          else if ( connectivity == 3 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (C_3)
             if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "C_3"
             !! Falta programar la aromaticidad   
             else if ( cyclomaticNumber >= 1 ) then
                isAromatic = AromaticityFinder_isAromatic( rings, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "C_R"
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
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (N_3)
             if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "N_3"
             !! Falta programar la aromaticidad   
             else if ( cyclomaticNumber >= 1 ) then
                isAromatic = AromaticityFinder_isAromatic( rings, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "N_R"
                end if
             else
                ffAtomType(i) = "N_2"
             end if
          else if ( connectivity == 2 ) then
             !! Se chequea el angulo promedio de enlace
             angleAverage = MMCommons_getAngleAverage( MolecularSystem_instance, i )
             !! Si el angulo es menor a 115 entonces se asume hibridacion sp3 (N_3)
             if ( angleAverage < SP2SP3AngleCutoff ) then
                ffAtomType(i) = "N_3"
                !! Falta programar la aromaticidad   
             else if ( cyclomaticNumber >= 1 ) then
                isAromatic = AromaticityFinder_isAromatic( rings, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "N_R"
                end if
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
             if ( cyclomaticNumber >= 1 ) then
                isAromatic = AromaticityFinder_isAromatic( rings, i )
                if ( isAromatic ) then
                   ffAtomType(i) = "O_R"
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
       ! else if( trim( labelOfCenters(i) ) == "P" ) then
       !    !! Se chequea la conectividad del oxigeno
       !    connectivity = MMCommons_getConnectivity( MolecularSystem_instance, i )
       !    if ( connectivity >= 4 ) then


       !    else
       !       ffAtomType(i) = "O_2"
       !    end if
!!******************************************************************************
       else
          ffAtomType(i) = labelOfCenters(i)
       end if
    end do

    do i=1, size(ffAtomType)
       write (*,"(T10,I,A,A,A,A)") i, ": ", trim( labelOfCenters(i) ), " ----> ", trim( ffAtomType(i) )
    end do
    
  end subroutine AtomTypeUFF_run

! subroutine AtomTypeUFF_exception( typeMessage, description, debugDescription)
!   implicit none
!   integer :: typeMessage


!   character(*) :: description
!   character(*) :: debugDescription

!   type(Exception) :: ex

!   call Exception_constructor( ex , typeMessage )
!   call Exception_setDebugDescription( ex, debugDescription )
!   call Exception_setDescription( ex, description )
!   call Exception_show( ex )
!   call Exception_destructor( ex )

! end subroutine AtomTypeUFF_exception

end module AtomTypeUFF_
