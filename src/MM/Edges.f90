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
!!        This module creates a class with the information of the edges in the system
!! @author  J.M. Rodas
!!
!! <b> Creation date : </b> 2014-06-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2014-06-02 </tt>: Jose Mauricio Rodas R. ( jmrodasr@unal.edu.co )
!!        -# Basics functions has been created
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
module Edges_
  use CONTROL_
  use MolecularSystem_
  use ParticleManager_
  use Matrix_
  use Vertex_
  use Rings_
  use MMCommons_
  use MatrixInteger_
  use Vector_
  use Exception_
  implicit none

  type , public :: Edges

     integer :: numberOfEdges
     type(MatrixInteger) :: connectionMatrix
     real(8), allocatable :: distance(:)
     real(8), allocatable :: bondOrder(:)
     real(8), allocatable :: idealDistance(:)
     real(8), allocatable :: forceConstant(:)
     real(8), allocatable :: stretchingEnergy(:) !! Kcal/mol
     real(8), allocatable :: stretchingEnergyKJ(:) !! KJ/mol

  end type Edges


       public :: &
            Edges_constructor, &
            Edges_getIdealDistance, &
            Edges_getDistance, &
            Edges_getBondOrders, &
            Edges_getOrder


contains

  !>
  !! @brief Defines the class constructor
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the edges
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] ring Class with the information of the rings
  !! @param numberOfEdges INTEGER number of edges in the system
  !! @param connectionMatrix INTEGER ARRAY with the information about the vertices in a bond
  !! @param distance REAL ARRAY with the bond distances(Angstroms) of the system
  !! @param bondOrder REAL ARRAY with the bond orders of the system
  !! @param idealDistance REAL ARRAY with the ideal bond distances(Angstroms) of the system
  !! @param forceConstant REAL ARRAY with the force constants of the system
  !! @param stretchingEnergy REAL ARRAY with the stretching energies (kcal/mol) of the system
  !! @param stretchingEnergyKJ REAL ARRAY with the stretching energies (kJ/mol) of the system
  subroutine Edges_constructor( this, vertices, ring )
    implicit none
    type(Edges), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    type(Rings), intent(in) :: ring
    integer :: i, j
    real(8), allocatable :: orders(:)

    call MMCommons_constructor( MolecularSystem_instance )
    
    this%numberOfEdges = size(MolecularSystem_instance%intCoordinates%distanceBondValue%values)

    call MatrixInteger_constructor( this%connectionMatrix, this%numberOfEdges, 2 )
    allocate( this%distance( this%numberOfEdges ) )
    allocate( this%bondOrder( this%numberOfEdges ) )

    call Edges_getBondOrders(orders, vertices, ring)

    do i=1,this%numberOfEdges
       do j=1,2
          this%connectionMatrix%values(i,j) = MolecularSystem_instance%intCoordinates%connectionMatrixForBonds%values(i,j)
       end do
       this%distance(i) = MolecularSystem_instance%intCoordinates%distanceBondValue%values(i) * ANGSTROM
       this%bondOrder(i) = orders(i)
    end do

    call Edges_getIdealDistance(this, vertices)

    call Edges_getForceConstants(this, vertices)

    call Edges_getStretchingEnergies(this)

  end subroutine Edges_constructor

  !>
  !! @brief This routine calculates the ideal distance of a bond using UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the edges
  !! @param [in] vertices Class with the information of the vertices
  !! @note The ideal distance (Natural Bond Length) is calculated using the equations 2-4 in Rappe et. al. paper (1992)
  !! with the correction proposed by Marcus G. Martin and implemented on TOWHEE (http://towhee.sourceforge.net/forcefields/uff.html) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[ 
  !! r_{ij} = r_{i} + r_{j} + r_{BO} - r_{EN}
  !! \f]
  !! \f[ 
  !! r_{BO} = -\lambda(r_{i} + r_{j})\ln(n)\:\:\:\:\:\lambda = 0.1332
  !! \f]
  !! \f[ 
  !! r_{EN} = \frac{r_{i}r_{j}(\sqrt{X_{i}}-\sqrt{X_{j}})^{2}}{X_{i}r_{i}+X_{j}r_{j}}
  !! \f]
  !! where: \n
  !! - \f$r_i\f$ and \f$r_j\f$ are the valence radius, parameters in the UFF
  !! - \f$X_i\f$ and \f$X_j\f$ are the electronegativities parameters in the UFF
  !! - \f$n\f$ is the bond order
  subroutine Edges_getIdealDistance(this, vertices)
    implicit none
    type(Edges), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    real(8) :: rbo !! Pauling-Type bond order correction
    real(8) :: ren !! Electronegativity correction
    integer :: i
    real(8) :: ri, rj, xi, xj
    real(8) :: lambda

    lambda = 0.1332

    allocate( this%idealDistance( this%numberOfEdges ) )

    do i=1, this%numberOfEdges
       ri=vertices%bondValence(this%connectionMatrix%values(i,1))
       rj=vertices%bondValence(this%connectionMatrix%values(i,2))
       xi=vertices%electronegativityGMP(this%connectionMatrix%values(i,1))
       xj=vertices%electronegativityGMP(this%connectionMatrix%values(i,2))
       !! equation 3
       rbo = -lambda*(ri+rj)*log(this%bondOrder(i)) 
       !! equation 4
       ren = ri*rj*(((sqrt(xi)-sqrt(xj))**2)/((xi*ri)+(xj*rj)))
       !! equation 2
       this%idealDistance(i) = ri + rj +rbo - ren
    end do

  end subroutine Edges_getIdealDistance

  !>
  !! @brief This function evaluates the ideal distance between two atoms
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class with the information of the edges
  !! @param [in] atomA INTEGER first atom to evaluate
  !! @param [in] atomB INTEGER second atom to evaluate
  !! @return [out] output REAL(8) ideal distance between atomA and atomB 
  function Edges_getDistance(this, atomA, atomB) result(output)
    implicit none
    type(Edges), intent(in) :: this
    integer, intent(in) :: atomA
    integer, intent(in) :: atomB
    real(8) :: output
    integer :: i,j

    do i=1, this%numberOfEdges
       do j=1,2
          if(this%connectionMatrix%values(i,1) == atomA .AND. this%connectionMatrix%values(i,2) == atomB) then
             output = this%idealDistance(i)
          else if(this%connectionMatrix%values(i,1) == atomB .AND. this%connectionMatrix%values(i,2) == atomA) then
             output = this%idealDistance(i)
          end if
       end do
    end do

  end function Edges_getDistance

  !>
  !! @brief This routine calculates the force constants using UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the edges
  !! @param [in] vertices Class with the information of the vertices
  !! @note The Force constants are calculated using the equation 6 in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[ 
  !! k_{ij} = 664.12\frac{Z_{i}^{*}Z_{j}^{*}}{r_{ij}^{3}}
  !! \f]
  !! where: \n
  !! - \f$Z_{i}^{*}\f$ and \f$Z_{j}^{*}\f$ are the effective atomic charges, parameters in the UFF
  !! - \f$r_{ij}\f$ is the ideal distance
  subroutine Edges_getForceConstants(this, vertices)
    implicit none
    type(Edges), intent(in out) :: this
    type(Vertex), intent(in) :: vertices
    integer :: i
    real(8) :: Zi, Zj !! Effective charges
    real(8) :: lambda

    lambda = 664.12

    allocate( this%forceConstant( this%numberOfEdges ) )

    do i=1, this%numberOfEdges
       Zi=vertices%effectiveCharge(this%connectionMatrix%values(i,1))
       Zj=vertices%effectiveCharge(this%connectionMatrix%values(i,2))
       !! equation 6
       this%forceConstant(i) = lambda*((Zi*Zj)/(this%idealDistance(i)**3))
    end do

  end subroutine Edges_getForceConstants

  !>
  !! @brief This routine calculates the Stretching energies using UFF parameters
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this Class with the information of the edges
  !! @note The Stretching energies are calculated using the equation 1a in Rappe et. al. paper (1992) \n
  !! A.K. Rappe, C.J. Casewit, K.S. Colwell, W.A. Goddard III, W.M. Skiff. 
  !! <b>UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular 
  !! Dynamics Simulations</b>. J. Am. Chem. Soc. 114, 10024-10035, 1992 \n
  !! \f[ 
  !! E_{R} = \frac{1}{2}k_{ij}(r-r_{ij})^{2}
  !! \f]
  !! where: \n
  !! - \f$k_{ij}\f$ is the force constant
  !! - \f$r\f$ is the bond distance
  !! - \f$r_{ij}\f$ is the ideal distance
  subroutine Edges_getStretchingEnergies(this)
    implicit none
    type(Edges), intent(in out) :: this
    integer :: i

    allocate( this%stretchingEnergy( this%numberOfEdges ) )
    allocate( this%stretchingEnergyKJ( this%numberOfEdges ) )

    do i=1, this%numberOfEdges
       !! equation 6
       this%stretchingEnergy(i) = 0.5*this%forceConstant(i)*((this%distance(i)-this%idealDistance(i))**2)
       this%stretchingEnergyKJ(i) = this%stretchingEnergy(i)*4.1868
    end do

  end subroutine Edges_getStretchingEnergies

  !>
  !! @brief This routine calculates the bond order using the connectivities, valences, distances, angles and aromaticity information
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] bondOrders REAL(8) ARRAY with the information of the bond orders
  !! @param [in] vertices Class with the information of the vertices
  !! @param [in] ring Class with the information of the rings
  !! @note 1 = single bond \n
  !! 2 = double bond \n
  !! 3 = triple bond \n
  !! 1.5 = aromatic bond \n
  !! 1.41 = amide C-N bond
  subroutine Edges_getBondOrders(bondOrders, vertices, ring)
    implicit none
    type(Vertex), intent(in) :: vertices
    type(Rings), intent(in) :: ring
    character(10), allocatable :: labelOfCenters(:)
    integer :: numberOfCenters
    integer :: numberOfEdges
    integer :: connectivity
    type(MatrixInteger) :: connectivityMatrix
    real(8), allocatable :: valences(:)
    logical :: isNeighborAromaticRing
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

    isNeighborAromaticRing = .false.

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
          bonds%values(i) = MolecularSystem_instance%intCoordinates%distanceBondValue%values(i) * ANGSTROM
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
                if(trim( vertices%type(i) ) == "C_R") then
                   isNeighborAromaticRing = Rings_isNeighborAromaticRing(ring, i, neighbor(j)) !!! ojo volver aqui
                   if(isNeighborAromaticRing) then
                      bondOrders(row) = 1.5
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1.0
                      valences(neighbor(j)) = valences(neighbor(j)) - 1.0
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else
                      bondOrders(row) = 1
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1
                      valences(neighbor(j)) = valences(neighbor(j)) - 1
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
          !************************************************************************************************
          !! Carbono insaturado o un carboanion o carbocation
          !************************************************************************************************
          else if( connectivityMatrix%values(i,2) < valences(i) ) then
             call MMCommons_searchEdgesRow( auxBonds, numberOfEdges, i, edgesRow )
             call MMCommons_searchNeighbor( auxBonds, numberOfEdges, size(edgesRow), i, neighbor )
             do j=1,size(edgesRow)
                row=edgesRow(j)
                if( trim( labelOfCenters(neighbor(j))) == "C" ) then
                   if(trim( vertices%type(i) ) == "C_R") then
                      bondOrders(row) = 1.5
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1.5
                      valences(neighbor(j)) = valences(neighbor(j)) - 1.5
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if( bonds%values(row)>singleBondCutoff ) then
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
                   if(trim( vertices%type(neighbor(j)) ) == "N_R") then
                      bondOrders(row) = 1.5
                      connectivityMatrix%values(i,2) = connectivityMatrix%values(i,2) - 1
                      connectivityMatrix%values(neighbor(j),2) = connectivityMatrix%values(neighbor(j),2) - 1
                      valences(i) = valences(i) - 1.5
                      valences(neighbor(j)) = valences(neighbor(j)) - 1.5
                      auxBonds(row)%values(1,1) = 0
                      auxBonds(row)%values(1,2) = 0
                   else if( bonds%values(row)>singleNBondCutoff ) then
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

  !>
  !! @brief This function evaluates the bond order between two atoms
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this Class with the information of the edges
  !! @param [in] atomA INTEGER first atom to evaluate
  !! @param [in] atomB INTEGER second atom to evaluate
  !! @return [out] output REAL(8) bond order between atomA and atomB 
  function Edges_getOrder(this, atomA, atomB) result(output)
    implicit none
    type(Edges), intent(in) :: this
    integer, intent(in) :: atomA, atomB
    real(8) :: output
    integer :: i, j

    do i=i,this%numberOfEdges
       do j=1,2
          if( this%connectionMatrix%values(i,1) == atomA .AND. this%connectionMatrix%values(i,2) == atomB ) then
             output = this%bondOrder(i)
          else if( this%connectionMatrix%values(i,1) == atomB .AND. this%connectionMatrix%values(i,2) == atomA ) then
             output = this%bondOrder(i)
          end if
       end do
    end do

  end function Edges_getOrder

  !>
  !! @brief Defines the class exception
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
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
