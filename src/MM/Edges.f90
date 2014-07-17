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
module Edges_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use MMCommons_
  use MatrixInteger_
  use Vector_
  use Exception_
  implicit none

  type , public :: Edges

     integer :: numberOfEdges
     type(MatrixInteger) :: connectionMatrix
     type(Vector) :: distance
     type(Vector) :: bondOrder

  end type Edges


       public :: &
            Edges_constructor, &
            Edges_getBondOrders


contains

  subroutine Edges_constructor( this )
    implicit none
    type(Edges), intent(in out) :: this
    integer :: i, j
    real(8), allocatable :: orders(:)

    
    call MMCommons_constructor( MolecularSystem_instance )
    this%numberOfEdges = size(MolecularSystem_instance%intCoordinates%distanceBondValue%values)

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfEdges, 2 )
    call Vector_constructor( this%distance, this%numberOfEdges )
    call Vector_constructor( this%bondOrder, this%numberOfEdges )

    call Edges_getBondOrders(orders)

    do i=1,this%numberOfEdges
       do j=1,2
          this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForBonds%values(i,j)
       end do
       this%distance%values(i) = MolecularSystem_instance%intCoordinates%distanceBondValue%values(i) * AMSTRONG
       this%bondOrder%values(i) = orders(i)
    end do

  end subroutine Edges_constructor


  subroutine Edges_getBondOrders(bondOrders)
    implicit none
    character(10), allocatable :: labelOfCenters(:)
    integer :: numberOfCenters
    integer :: numberOfEdges
    integer :: connectivity
    type(MatrixInteger) :: connectivityMatrix
    real(8), allocatable :: valences(:)
    ! type(Vector) :: bonds
    ! type(Exception) :: ex
    integer :: i, j, row
    type(MatrixInteger), allocatable :: auxBonds (:)
    real(8), allocatable, intent(in out) :: bondOrders(:)
    integer, allocatable :: neighbor(:)
    integer, allocatable :: edgesRow(:)
    type(Vector) :: bonds
    real(8) :: aromaticBondCutoff
!! Enlaces de carbono
    real(8) :: singleBondCutoff
    real(8) :: doubleBondCutoff
    real(8) :: singleOBondCutoff
    real(8) :: singleNBondCutoff
    real(8) :: doubleNBondCutoff
    real(8) :: singleSBondCutoff
    real(8) :: singleTeBondCutoff
    real(8) :: singleAsBondCutoff
    real(8) :: singlePBondCutoff
!! Enlaces de Nitrogeno
    real(8) :: singleONBondCutoff
    real(8) :: singleNNBondCutoff
    real(8) :: singleSNBondCutoff
    real(8) :: singlePNBondCutoff
    real(8) :: singleSeNBondCutoff
    real(8) :: singleAsNBondCutoff
    real(8) :: doubleNNBondCutoff
!! Enlaces de Oxigeno
    real(8) :: singleOOBondCutoff
    real(8) :: singleSOBondCutoff
    real(8) :: singlePOBondCutoff
    real(8) :: singleAsOBondCutoff
!! Enlaces de Fosforo
    real(8) :: singlePPBondCutoff
    real(8) :: singleSPBondCutoff
!! Enlaces de Azufre
    real(8) :: singleSSBondCutoff
    real(8) :: singleAsSBondCutoff
!! Enlaces de Arsenico
    real(8) :: singleAsAsBondCutoff

!!--------------------------------------------------------------
!! Enlaces de carbono
!!--------------------------------------------------------------
    singleBondCutoff = 1.4200000
    doubleBondCutoff = 1.2200000
    singleOBondCutoff = 1.300000
    singleNBondCutoff = 1.380000
    doubleNBondCutoff = 1.200000
    singleSBondCutoff = 1.700000
    singlePBondCutoff = 1.720000
    singleTeBondCutoff = 2.100000
    singleAsBondCutoff = 1.900000
    aromaticBondCutoff = 1.3600000
!!--------------------------------------------------------------
!! Enlaces de Nitrogeno
!!--------------------------------------------------------------
    singleONBondCutoff = 1.250000
    singleNNBondCutoff = 1.320000
    singleSNBondCutoff = 1.580000
    singlePNBondCutoff = 1.620000
    singleSeNBondCutoff = 1.800000
    singleAsNBondCutoff = 1.845000
    doubleNNBondCutoff = 1.100000
!!--------------------------------------------------------------
!! Enlaces de Oxigeno
!!--------------------------------------------------------------
    singleOOBondCutoff = 1.230000
    singleSOBondCutoff = 1.540000
    singlePOBondCutoff = 1.520000
    singleAsOBondCutoff = 1.680000
!!--------------------------------------------------------------
!! Enlaces de Fosforo
!!--------------------------------------------------------------
    singlePPBondCutoff = 2.060000
    singleSPBondCutoff = 1.950000
!!--------------------------------------------------------------
!! Enlaces de  Azufre
!!--------------------------------------------------------------
    singleSSBondCutoff = 1.900000
    singleAsSBondCutoff = 2.150000
!!--------------------------------------------------------------
!! Enlaces de  Arsenico
!!--------------------------------------------------------------
    singleAsAsBondCutoff = 2.390000



    call MMCommons_constructor( MolecularSystem_instance )

    numberOfCenters = ParticleManager_getNumberOfCentersOfOptimization()
    numberOfEdges=size(MolecularSystem_instance%intCoordinates%distanceBondValue%values)

    call MatrixInteger_constructor( connectivityMatrix, numberOfCenters, 2 )
    call MMCommons_getConnectivityMatrix( MolecularSystem_instance, numberOfCenters, connectivityMatrix )

    allocate( labelOfCenters( numberOfCenters ) )
    labelOfCenters = ParticleManager_getLabelsOfCentersOfOptimization()

    call MMCommons_getValences(labelOfCenters, numberOfCenters, valences)

    ! write(*,"(T20,A)") ""
    ! write(*,"(T20,A)") "Informacion en Bond orders"
    ! do i=1,numberOfCenters
    !    write(*,"(T20,I5,2x,A1,I5,2x,F8.5)") connectivityMatrix%values(i,1), labelOfCenters(i), connectivityMatrix%values(i,2), valences(i)
    ! end do
    ! write(*,"(T20,A)") ""

    allocate( auxBonds( numberOfEdges )  )
    do i=1,numberOfEdges
       call MatrixInteger_constructor( auxBonds(i), 1, 2 )
    end do

    call Vector_constructor( bonds, numberOfEdges )
    do i=1,numberOfEdges
       do j=1,2
          bonds%values(i) = MolecularSystem_instance%intCoordinates%distanceBondValue%values(i) * AMSTRONG
       end do
    end do

    allocate( bondOrders( numberOfEdges ) )
    do i=1,numberOfEdges
       do j=1,2
          auxBonds(i)%values(1,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForBonds%values(i,j)
       end do
       bondOrders(i) = 0
    end do

    do i=1,numberOfCenters
!**********************************************************************************************************
!! Enlaces con hidrogeno primero
!**********************************************************************************************************
       if( trim( labelOfCenters(i) ) == "H" .AND. connectivityMatrix%values(i,2) > 0 ) then
          call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
          call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
          do j=1,size(edgesRow)
             row=edgesRow(j)
             bondOrders(row) = 1
             connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
             connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
             valences(i) = valences(i) - 1
             valences(neighbor(j)) = valences(neighbor(j)) - 1
             auxBonds(row)%values(1,1) = 0
             auxBonds(row)%values(1,2) = 0
          end do
!**********************************************************************************************************
!! Enlaces de carbono
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "C" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Carbono solo con enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Carbono insaturado o un carboanion o carbocation
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if(bonds%values(row)<singleBondCutoff .AND. &
                        bonds%values(row)>aromaticBondCutoff ) then
                      bondOrders(row) = 1.5
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1.5
                      valences(neighbor(j)) = valences(neighbor(j)) - 1.5
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if(bonds%values(row)<doubleBondCutoff) then
                      bondOrders(row) = 3
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 3
                      valences(neighbor(j)) = valences(neighbor(j)) - 3
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singleOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if(bonds%values(row)<doubleNBondCutoff) then
                      bondOrders(row) = 3
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 3
                      valences(neighbor(j)) = valences(neighbor(j)) - 3
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleSBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "TE" ) then
                   if( bonds%values(row)>singleTeBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "P" ) then
                   if( bonds%values(row)>singlePBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "AS" ) then
                   if( bonds%values(row)>singleAsBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de nitrogeno
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "N" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Nitrogeno con conectividad == valencia
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Nitrogeno con conectividades menores a 5 (la mayoria de los casos)
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if(bonds%values(row)<doubleNBondCutoff) then
                      bondOrders(row) = 3
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 3
                      valences(neighbor(j)) = valences(neighbor(j)) - 3
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singleONBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleNNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if(bonds%values(row)<doubleNNBondCutoff) then
                      bondOrders(row) = 3
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 3
                      valences(neighbor(j)) = valences(neighbor(j)) - 3
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleSNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "SE" ) then
                   if( bonds%values(row)>singleSeNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "P" ) then
                   if( bonds%values(row)>singlePNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "AS" ) then
                   if( bonds%values(row)>singleAsNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Oxigeno
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "O" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Oxigeno con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Oxigeno con conectividades menores a 2
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singleOOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleONBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleSOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "P" ) then
                   if( bonds%values(row)>singlePOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "AS" ) then
                   if( bonds%values(row)>singleAsOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Fosforo
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "P" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Fosforo con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Fosforo con conectividades menores a a la valencias posibles insaturaciones
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singlePBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singlePOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singlePNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "P" ) then
                   if( bonds%values(row)>singlePPBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleSPBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Azufre
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "S" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Azufre con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Azufre con conectividades menores a a la valencias posibles insaturaciones
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleSBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singleSOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleSNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "P" ) then
                   if( bonds%values(row)>singleSPBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleSSBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "AS" ) then
                   if( bonds%values(row)>singleAsSBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Arsenico
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "AS" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Arsenico con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Arsenico con conectividades menores a a la valencias posibles insaturaciones
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleAsBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "O" ) then
                   if( bonds%values(row)>singleAsOBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleAsNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "S" ) then
                   if( bonds%values(row)>singleAsSBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else if( trim( labelOfCenters(neighbor(j))) == "AS" ) then
                   if( bonds%values(row)>singleAsAsBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Selenio
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "SE" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Selenio con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Selenio con conectividades menores a a la valencias posibles insaturaciones
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "N" ) then
                   if( bonds%values(row)>singleSeNBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces de Telurio
!**********************************************************************************************************
       else if( trim( labelOfCenters(i) ) == "TE" .AND. connectivityMatrix%values(i,2) > 0 ) then
          !************************************************************************************************
          !! Telurio con conectividad >= valencia, enlaces sencillos
          !************************************************************************************************
          if( connectivityMatrix%values(i,2) >= valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                bondOrders(row) = 1
                connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                valences(i) = valences(i) - 1
                valences(neighbor(j)) = valences(neighbor(j)) - 1
                auxBonds(row)%values(1,1) = 0
                auxBonds(row)%values(1,2) = 0
             end do
          !************************************************************************************************
          !! Telurio con conectividades menores a a la valencias posibles insaturaciones
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if( bonds%values(row)>singleTeBondCutoff ) then
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 2
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 2
                      valences(neighbor(j)) = valences(neighbor(j)) - 2
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   end if
                else
                   bondOrders(row) = 1
                   connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                   connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                   valences(i) = valences(i) - 1
                   valences(neighbor(j)) = valences(neighbor(j)) - 1
                   auxBonds(row)%values(1,1) = 0
                   auxBonds(row)%values(1,2) = 0
                end if
             end do
          end if
!**********************************************************************************************************
!! Enlaces del resto de los elementos
!**********************************************************************************************************
       else
          call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
          call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
          do j=1,size(edgesRow)
             row=edgesRow(j)
             bondOrders(row) = 1
             connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
             connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
             valences(i) = valences(i) - 1
             valences(neighbor(j)) = valences(neighbor(j)) - 1
             auxBonds(row)%values(1,1) = 0
             auxBonds(row)%values(1,2) = 0
          end do
       end if
    end do

  end subroutine Edges_getBondOrders

  subroutine Edges_exception( typeMessage, description, debugDescription)
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

  end subroutine Edges_exception

end module Edges_
