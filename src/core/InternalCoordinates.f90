!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!	Prof. G. MERINO's Lab. Universidad de Guanajuato
!!		http://quimera.ugto.mx/qtc/gmerino.html
!!
!!	Authors:
!!		E. F. Posada (efposadac@unal.edu.co)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
!!
!!******************************************************************************

  !>
  !! @brief Clase encargada de obtener coordenadas primitivas para proceso de
  !!             minimizacion multidimencional (Optimizacion de geometria)
  !!
  !! Esta clase manipula toda la informacion necesaria del sistema molecular para
  !! para obtener un conjunto de coordenadas primitivas que puedas\n ser empleadas
  !! eficientemente en un metodo de minimizacion multidimensional, para realizar
  !! la optimizacion geometrica de nucleos tratados dentro de la aproximacion de
  !! Born-Oppenheimer y centros de funciones gausianas asociadas a nucleos dentro
  !! del formalismo de la teoria del orbital nuclear y electronico.
  !!
  !! El metodo para la optimizacion ha sido tomado de:
  !!	 Reveles and Koster, J Comp Chem, 25, 9, 2004, p1109-1116.
  !!
  !! @author Sergio Gonzalez
  !!
  !! <b> Fecha de creacion : </b> 2009-04-05
  !!   - <tt> 2009-04-05 </tt>: Sergio Gonzalez ( sagonzalez@unal.edu.co )
  !!        -# Creacion del archivo y las funciones basicas
  !!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
  !!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
  !!   - <tt> 2014-05-14 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)
  !!        -# Reescribe y adapta el módulo para su inclusion en Lowdin2
  !!
  !<

module InternalCoordinates_
  use ParticleManager_
  use Matrix_
  use MatrixInteger_
  use Vector_
  use Map_
  use Math_
  use List_
  use ListInteger_
  use AtomicElement_
  use Exception_
  implicit none

  type , public :: InternalCoordinates

     character(30) :: name
     type(Matrix) :: cartesianCoordinates
     type(Matrix) :: wilsonMatrix
     type(Matrix) :: symmetricGMatrix
     type(Matrix) :: nonRedundantEigenvectors !< Contiene los vectores propios no redundates del espacio de coordenadas original
     type(Matrix) :: BWilsonMatrix
     type(Matrix) :: inverseOfBWilsonMatrix

     !>
     !!  Matrices para definicion de conectividad quimica
     !<
     type(MatrixInteger) :: connectionMatrixForBonds
     type(Vector) :: distanceBondValue
     type(MatrixInteger) :: connectionMatrixForAngles
     type(Vector) :: angleOfBondValue
     type(MatrixInteger) :: connectionMatrixForDihedrals
     type(Vector) :: dihedralsAngleValue

     !>
     !!   Almacena informacion sobre numero de parametros
     !<
     integer :: numberOfBonds
     integer :: numberOfAnglesOfBond
     integer :: numberOfDihedrals
     integer :: numberOfPrimitivesCoordinates
     integer :: numberOfCenterOfOptimization

     type(Map) :: covalentRadius!<  Mapa empleado por conveniencia

  end type InternalCoordinates

  public :: &
       InternalCoordinates_constructor, &
       InternalCoordinates_destructor, &
       InternalCoordinates_show, &
       InternalCoordinates_obtainCoordinates,&
       InternalCoordinates_getStartHessian, &
       InternalCoordinates_getBonds, &
       InternalCoordinates_getAnglesOfBond

  private		
contains


  !>
  !! Define el constructor para la clase
  !!
  !<
  subroutine InternalCoordinates_constructor( this )
    implicit none
    type(InternalCoordinates) :: this

    this%numberOfCenterOfOptimization = ParticleManager_getNumberOfCentersOfOptimization()
    this%cartesianCoordinates = ParticleManager_getCartesianMatrixOfCentersOfOptimization()

  end subroutine InternalCoordinates_constructor

  !>
  !! @brief Define el destructor para clase
  !!
  !! @param thisPointer Funcion base
  !<
  subroutine InternalCoordinates_destructor( this )
    implicit none
    type(InternalCoordinates) :: this

    call Matrix_destructor( this%cartesianCoordinates )
    call Matrix_destructor( this%cartesianCoordinates)
    call Matrix_destructor( this%wilsonMatrix)
    call Matrix_destructor( this%symmetricGMatrix)
    call Matrix_destructor( this%nonRedundantEigenvectors)
    call Matrix_destructor( this%BWilsonMatrix )
    call Matrix_destructor( this%inverseOfBWilsonMatrix )
    call MatrixInteger_destructor( this%connectionMatrixForBonds)
    call MatrixInteger_destructor( this%connectionMatrixForAngles)
    call MatrixInteger_destructor( this%connectionMatrixForDihedrals)
    call Vector_destructor( this%distanceBondValue)
    call Vector_destructor( this%angleOfBondValue)
    call Vector_destructor( this%dihedralsAngleValue)
    call Map_destructor( this%covalentRadius)

  end subroutine InternalCoordinates_destructor

  !>
  !! @brief Muestra los atributos de todas las part\'iculas en el Administrador de part\'iculas
  !!
  !! @todo Falta adicionar procedimiento para que muestre una sola part\'icula
  !<
  subroutine InternalCoordinates_show( this  )
    implicit none
    type(InternalCoordinates) :: this

  end subroutine InternalCoordinates_show

  !>
  !! @brief Selecciona enlances interatomicos de acuerdo a criterio de conectividad  \f$ R_{AB} < Fs (r_{WA}+r_{WB}) \f$
  !!
  !! Donde:
  !!
  !! <table>
  !! <tr> <td> \f[Fs:\f]   <td> <dfn> Factor de Conectividad. </dfn>
  !! <tr> <td> \f[r_{WX}:\f]  <td>  <dfn> Radio de Van der Waals para atomo X. </dfn>
  !! <tr> <td> \f[R_{AB}:\f] <td> <dfn> distancia entre atomo A y B..</dfn>
  !! </table>
  !! @msc
  !!  Bonds;
  !!    Bonds=>Bonds [label="getBonds(BEGIN)"];
  !!    Bonds=>Bonds [label="defineConnectivity()"];
  !!    Bonds=>Bonds [label="calculateAtomicDistance(AB)"];
  !!    Bonds=>Bonds [label="selectBond(if R_AB< Threshold)"];
  !!    Bonds=>Bonds [label="getBonds(END)"];
  !!    --- [label="Return matrix of bonds "];

  !<
  function InternalCoordinates_getBonds( this, distanceFactor ) result( output )
    implicit none
    type( InternalCoordinates ) :: this
    real(8), intent(in) :: distanceFactor
    type(MatrixInteger) :: output

    type(List) :: bondDistance
    type(ListInteger) :: currentAtom
    type(ListInteger) :: otherAtom
    character(10), allocatable :: labelOfCenters(:)
    real(8) :: connectivityCriteria
    real(8) :: separationOfCenters
    real(8) :: covalentRadiusValues(2)
    integer :: i
    integer :: j
    integer :: auxIndex

    call List_constructor( bondDistance, ssize=-1 )
    call ListInteger_constructor( currentAtom, ssize=-1 )
    call ListInteger_constructor( otherAtom, ssize=-1 )
    call Map_constructor(this%covalentRadius)

    allocate( labelOfCenters( this%numberOfCenterOfOptimization ) )
    labelOfCenters = ParticleManager_getLabelsOfCentersOfOptimization()
    
    do i=1,this%numberOfCenterOfOptimization
       do j=i+1, this%numberOfCenterOfOptimization

          separationOfCenters=dsqrt( sum( ( this%cartesianCoordinates%values(i,:) - this%cartesianCoordinates%values(j,:))**2.0) )

          auxIndex = Map_find( this%covalentRadius, labelOfCenters(i) )
          if( auxIndex == 0 ) then
             covalentRadiusValues(1) = AtomicElement_getCovalentRadius( labelOfCenters(i) )
             call Map_insert( this%covalentRadius, labelOfCenters(i), covalentRadiusValues(1) )
          else
             covalentRadiusValues(1) = this%covalentRadius%value( auxIndex )
          end if

          auxIndex = Map_find( this%covalentRadius, labelOfCenters(j) )
          if ( auxIndex == 0 ) then
             covalentRadiusValues(2) = AtomicElement_getCovalentRadius( labelOfCenters(j) )
             call Map_insert ( this%covalentRadius, labelOfCenters(j), covalentRadiusValues(2) )
          else
             covalentRadiusValues(2) = this%covalentRadius%value( auxIndex )
          end if
         
          connectivityCriteria =  (	distanceFactor * sum( covalentRadiusValues ) ) /  AMSTRONG

          if ( separationOfCenters < connectivityCriteria ) then

             call ListInteger_push_back(currentAtom, i)
             call ListInteger_push_back(otherAtom, j)
             call List_push_back( bondDistance, separationOfCenters )

          end if

       end do
    end do

    call MatrixInteger_constructor( output, ListInteger_size(currentAtom), 2,0 )
    call Vector_constructor(this%distanceBondValue,ListInteger_size(currentAtom) )
    output%values(:,1) = currentAtom%data(:)
    output%values(:,2) = otherAtom%data(:)
    this%distanceBondValue%values = bondDistance%data

    this%numberOfBonds = size(output%values, dim=1)

    deallocate( labelOfCenters )
    call List_destructor( bondDistance )
    call ListInteger_destructor( currentAtom )
    call ListInteger_destructor( otherAtom )
    !!
    !!********************************************************

  end function InternalCoordinates_getBonds

  !>
  !! @brief Selecciona los angulos de enlace formado entre pares de enlaces definidos
  !! que tienen un atomo en comun
  !! v.g A-B y B-C o A-B y C-B  -> A-B-C
  !!
  !! @msc
  !!  Bonds,Angles;
  !!    Bonds=>Bonds [label="getBonds()"];
  !!    Bonds>>Angles[label="connectionMatrixForBonds[...]"];
  !!    Angles=>Angles [ label = "getAngles(BEGIN)"];
  !!    Angles=>Angles [ label = "defineConnectivity()"];
  !!    Angles=>Angles [ label = "calculateAngles(ABC)"];
  !!    Angles=>Angles [label="selectAngle(if arc(ABC) < Threshold)"];
  !!    Angles=>Angles [ label = "getAngles(END)"];
  !!    --- [label="Return matrix of angles of bond "];
  !!  @endmsc
  !<
  function InternalCoordinates_getAnglesOfBond( this, angleThreshold ) result( output )
    implicit none
    type( InternalCoordinates ) :: this
    real(8), intent(in) :: angleThreshold
    type(MatrixInteger) :: output

    type(Vector) :: vectorA
    type(Vector) :: vectorB
    real(8) :: angleValue
    type(List) :: anglesOfBond
    type(ListInteger) :: currentAtom
    type(ListInteger) :: otherAtom
    type(ListInteger) :: atomForDefineAngle
    integer :: connectivity(3)
    integer :: otherConnectivity(3)
    integer :: i
    integer :: j
    integer :: k

    if ( allocated(this%connectionMatrixForBonds%values) ) then

       if ( size(this%connectionMatrixForBonds%values) >=2 ) then
          call ListInteger_constructor( currentAtom,ssize=-1 )
          call ListInteger_constructor( otherAtom, ssize=-1 )
          call ListInteger_constructor( atomForDefineAngle, ssize=-1 )
          call List_constructor( anglesOfBond, ssize=-1 )

          call Vector_constructor( vectorA, 3 )
          call Vector_constructor( vectorB, 3 )

          this%numberOfBonds = size(this%connectionMatrixForBonds%values, dim=1)

          do i=1,this%numberOfBonds
             do j=i+1, this%numberOfBonds

                connectivity = 0

                if ( this%connectionMatrixForBonds%values(i,1) == this%connectionMatrixForBonds%values(j,1) ) then

                   connectivity = [this%connectionMatrixForBonds%values(i,2), this%connectionMatrixForBonds%values(i,1),&
                        this%connectionMatrixForBonds%values(j,2) ]

                else if ( this%connectionMatrixForBonds%values(i,1) == this%connectionMatrixForBonds%values(j,2) ) then

                   connectivity = [this%connectionMatrixForBonds%values(i,2), this%connectionMatrixForBonds%values(i,1), &
                        this%connectionMatrixForBonds%values(j,1) ]

                else if ( this%connectionMatrixForBonds%values(i,2) == this%connectionMatrixForBonds%values(j,1) ) then
                   connectivity = [this%connectionMatrixForBonds%values(i,1), this%connectionMatrixForBonds%values(i,2), &
                        this%connectionMatrixForBonds%values(j,2) ]

                else if ( this%connectionMatrixForBonds%values(i,2) == this%connectionMatrixForBonds%values(j,2)  ) then

                   connectivity = [this%connectionMatrixForBonds%values(i,1), this%connectionMatrixForBonds%values(i,2), &
                        this%connectionMatrixForBonds%values(j,1) ]

                end if


                if ( connectivity(1) /= 0 ) then

                   vectorA%values = this%cartesianCoordinates%values( connectivity(1), : ) &
                        - this%cartesianCoordinates%values( connectivity(2), : )

                   vectorB%values = this%cartesianCoordinates%values( connectivity(3), : ) &
                        - this%cartesianCoordinates%values( connectivity(2), : )

                   angleValue =acos( dot_product( vectorA%values,vectorB%values) &
                        / ( sqrt( dot_product( vectorA%values, vectorA%values ) )  &
                        * sqrt( dot_product( vectorB%values, vectorB%values ) ) ) ) * DEGREES

                   ! if ( angleValue <= angleThreshold ) then
                      call ListInteger_push_back( currentAtom, connectivity(1) )
                      call ListInteger_push_back( otherAtom, connectivity(2) )
                      call ListInteger_push_back( atomForDefineAngle, connectivity(3) )
                      call List_push_back( anglesOfBond, angleValue )
                   ! end if

                end if

             end do
          end do

          call MatrixInteger_constructor( output, ListInteger_size(currentAtom), 3,0 )
          call Vector_constructor( this%angleOfBondValue, ListInteger_size(currentAtom) )
          output%values(:,1) = currentAtom%data(:)
          output%values(:,2) = otherAtom%data(:)
          output%values(:,3) = atomForDefineAngle%data(:)
          this%angleOfBondValue%values = anglesOfBond%data

          !!***********************************************************************
          !! Remueve angulos de anlace equivalentes
          !!
          i=1
          do while( i < size(output%values, dim=1) )

             connectivity= output%values(i,:)

             j=i+1
             do while ( j < size(output%values, dim=1) )

                if ( 	connectivity(2) == output%values(j,1) .and. connectivity(3) == output%values(j,2) .and. &
                     connectivity(1) == output%values(j,3)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%angleOfBondValue,j)

                else &
                     if ( 	connectivity(3) == output%values(j,1) .and. connectivity(1) == output%values(j,2) .and. &
                     connectivity(2) == output%values(j,3)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%angleOfBondValue,j)

                else &
                     if ( 	connectivity(1) == output%values(j,1) .and. connectivity(3) == output%values(j,2) .and. &
                     connectivity(2) == output%values(j,3)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%angleOfBondValue,j)

                else &
                     if ( 	connectivity(3) == output%values(j,1) .and. connectivity(2) == output%values(j,2) .and. &
                     connectivity(1) == output%values(j,3)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%angleOfBondValue,j)

                else &
                     if ( 	connectivity(2) == output%values(j,1) .and. connectivity(1) == output%values(j,2) .and. &
                     connectivity(3) == output%values(j,3)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%angleOfBondValue,j)

                else
                   j=j+1
                end if

             end do
             i=i+1
          end do

          !!
          !!***********************************************************************

          this%numberOfAnglesOfBond = size(output%values, dim=1)

          call ListInteger_destructor( currentAtom )
          call ListInteger_destructor( otherAtom )
          call ListInteger_destructor( atomForDefineAngle )
          call List_destructor( anglesOfBond )
          call Vector_destructor( vectorA )
          call Vector_destructor( vectorB )

       end if

    else

       call InternalCoordinates_exception( ERROR, "The chemical bonds haven't been defined",&
            "Class object InternalCoordinates in the getAnglesOfBond() function")

    end if


  end function InternalCoordinates_getAnglesOfBond

  !>
  !! @brief Selecciona los angulos dihedros formados entre los planos de los angulos primitivos
  !! que tienen dos atomos en comun
  !! v.g A-B-C  y  B-C-D o A-B-C y D-A-B  -> A-B-C-D
  !! @msc
  !! Bonds, Angles, Dihedrals;
  !!    Bonds=>Bonds [label="getBonds()"];
  !!    Bonds>>Angles[label="connectionMatrixForBonds[...]"];
  !!    Angles=>Angles [ label = "getAngles()"];
  !!    Angles>>Dihedrals [ label = "connectionMatrixForAngles[...]"];
  !!    Dihedrals=>Dihedrals [ label = "getDihedralsAngles(BEGIN)"];
  !!    Dihedrals=>Dihedrals [ label = "defineConnectivity()"];
  !!    Dihedrals=>Dihedrals [ label = "calculateAngles(ABCD)"];
  !!    Dihedrals=>Dihedrals [label="selectAngle(if arc(ABCD) < Threshold)"];
  !!    Dihedrals=>Dihedrals [ label = "getDihedralsAngles(END)"];
  !!    --- [label="Return matrix of dihedrals angles"];
  !!  @endmsc
  !<
  function InternalCoordinates_getDihedralAngles( this, angleThreshold ) result( output )
    implicit none
    type( InternalCoordinates ) :: this
    real(8), intent(in) :: angleThreshold
    type(MatrixInteger) :: output

    type(Vector) :: vectorA
    type(Vector) :: vectorB
    type(Vector) :: auxVector
    real(8) :: dihedralAngle
    real(8) :: angleSign
    type(List) :: dihedralAngles
    type(ListInteger) :: currentAtom
    type(ListInteger) :: otherAtom
    type(ListInteger) :: atomForDefineAngle
    type(ListInteger) :: atomForDefineDihedral
    integer :: connectivity(4)
    integer :: i
    integer :: j

    if ( allocated(this%connectionMatrixForAngles%values)   ) then

       if (size(this%angleOfBondValue%values)>=2) then

          call ListInteger_constructor( currentAtom,ssize=-1 )
          call ListInteger_constructor( otherAtom, ssize=-1 )
          call ListInteger_constructor( atomForDefineAngle, ssize=-1 )
          call ListInteger_constructor( atomForDefineDihedral, ssize=-1)
          call List_constructor( dihedralAngles, ssize=-1 )

          call Vector_constructor( vectorA, 3 )
          call Vector_constructor( vectorB, 3 )
          call Vector_constructor( auxVector, 3)

          this%numberOfAnglesOfBond = size(this%connectionMatrixForAngles%values, dim=1)

          do i=1,this%numberOfAnglesOfBond
             do j=i+1, this%numberOfAnglesOfBond

                connectivity = 0

                if (  this%connectionMatrixForAngles%values(i,2) == this%connectionMatrixForAngles%values(j,1) ) then

                   if ( this%connectionMatrixForAngles%values(i,1) == this%connectionMatrixForAngles%values(j,2) ) then

                      connectivity = [this%connectionMatrixForAngles%values(i,3), this%connectionMatrixForAngles%values(i,2), &
                           this%connectionMatrixForAngles%values(i,1), this%connectionMatrixForAngles%values(j,3) ]

                   else if ( this%connectionMatrixForAngles%values(i,3) == this%connectionMatrixForAngles%values(j,2) ) then

                      connectivity = [this%connectionMatrixForAngles%values(i,1), this%connectionMatrixForAngles%values(i,2), &
                           this%connectionMatrixForAngles%values(i,3), this%connectionMatrixForAngles%values(j,3) ]

                   end if

                else if ( this%connectionMatrixForAngles%values(i,2) == this%connectionMatrixForAngles%values(j,3) ) then

                   if ( this%connectionMatrixForAngles%values(i,1) == this%connectionMatrixForAngles%values(j,2) ) then

                      connectivity = [this%connectionMatrixForAngles%values(j,1), this%connectionMatrixForAngles%values(j,2),&
                           this%connectionMatrixForAngles%values(j,3), this%connectionMatrixForAngles%values(i,3) ]

                   else if ( this%connectionMatrixForAngles%values(i,3) == this%connectionMatrixForAngles%values(j,2) ) then

                      connectivity = [this%connectionMatrixForAngles%values(j,1), this%connectionMatrixForAngles%values(j,2), &
                           this%connectionMatrixForAngles%values(j,3), this%connectionMatrixForAngles%values(i,1) ]

                   end if

                end if

                if ( connectivity(1) /= 0 ) then

                   vectorA%values = this%cartesianCoordinates%values( connectivity(1),: ) &
                        - this%cartesianCoordinates%values( connectivity(2),: )

                   vectorB%values = this%cartesianCoordinates%values( connectivity(3),: ) &
                        - this%cartesianCoordinates%values( connectivity(2),: )

                   auxVector = Vector_cross( vectorA, vectorB )

                   vectorA%values = this%cartesianCoordinates%values( connectivity(4),: ) &
                        - this%cartesianCoordinates%values( connectivity(2),: )

                   vectorB = Vector_cross( vectorA, vectorB )
                   vectorA = auxVector

                   !!
                   !! Calcula el signo del  angulo diedro proyectando el vector V_zy sobre el segundo angulo
                   !! de enlace (j) sobre la normal del plano definido por el primer angulo de enlace (i).
                   !!
                   auxVector%values = this%cartesianCoordinates%values( connectivity(3),: ) &
                        - this%cartesianCoordinates%values( connectivity(4),: )

                   angleSign = dsign(1.0_8,dot_product(auxVector%values,vectorA%values ) / Vector_norm(vectorA) )

                   dihedralAngle = angleSign * acos( dot_product( vectorA%values,vectorB%values) &
                        / ( sqrt( dot_product( vectorA%values, vectorA%values ) )  &
                        * sqrt( dot_product( vectorB%values, vectorB%values ) ) ) ) * DEGREES

                   if ( dihedralAngle < angleThreshold ) then

                      call ListInteger_push_back( currentAtom, connectivity(1) )
                      call ListInteger_push_back( otherAtom, connectivity(2) )
                      call ListInteger_push_back( atomForDefineAngle, connectivity(3) )
                      call ListInteger_push_back( atomForDefineDihedral, connectivity(4) )
                      call List_push_back( dihedralAngles, dihedralAngle )

                   end if

                end if

             end do
          end do

          call MatrixInteger_constructor( output, ListInteger_size(currentAtom), 4,0 )
          call Vector_constructor( this%dihedralsAngleValue, ListInteger_size(currentAtom) )
          output%values(:,1) = currentAtom%data(:)
          output%values(:,2) = otherAtom%data(:)
          output%values(:,3) = atomForDefineAngle%data(:)
          output%values(:,4) = atomForDefineDihedral%data(:)
          this%dihedralsAngleValue%values = dihedralAngles%data

          !!***********************************************************************
          !! Remueve angulos diedros equivalentes
          !!
          i=1
          do while( i < size(output%values, dim=1) )

             connectivity= output%values(i,:)

             j=i+1
             do while ( j < size(output%values, dim=1) )

                if ( 	connectivity(2) == output%values(j,1) .and. connectivity(3) == output%values(j,2) .and. &
                     connectivity(4) == output%values(j,3) .and. connectivity(1) == output%values(j,4)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else &
                     if (	connectivity(3) == output%values(j,1) .and. connectivity(4) == output%values(j,2) .and. &
                     connectivity(1) == output%values(j,3) .and. connectivity(2) == output%values(j,4) ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else &
                     if (	connectivity(4) == output%values(j,1) .and. connectivity(1) == output%values(j,2) .and. &
                     connectivity(2) == output%values(j,3) .and. connectivity(3) == output%values(j,4) ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else &
                     if (	connectivity(2) == output%values(j,1) .and. connectivity(1) == output%values(j,2) .and. &
                     connectivity(4) == output%values(j,3) .and. connectivity(3) == output%values(j,4) ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)


                else &
                     if ( 	connectivity(3) == output%values(j,1) .and. connectivity(2) == output%values(j,2) .and. &
                     connectivity(1) == output%values(j,3) .and. connectivity(4) == output%values(j,4)  ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else &
                     if (	connectivity(4) == output%values(j,1) .and. connectivity(3) == output%values(j,2) .and. &
                     connectivity(2) == output%values(j,3) .and. connectivity(1) == output%values(j,4) ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else &
                     if (	connectivity(1) == output%values(j,1) .and. connectivity(4) == output%values(j,2) .and. &
                     connectivity(3) == output%values(j,3) .and. connectivity(2) == output%values(j,4) ) then

                   call MatrixInteger_removeRow(output, j)
                   call Vector_removeElement(this%dihedralsAngleValue,j)

                else
                   j=j+1
                end if

             end do
             i=i+1
          end do
          !!
          !!***********************************************************************

          this%numberOfDihedrals = size(output%values, dim=1)

          call ListInteger_destructor( currentAtom )
          call ListInteger_destructor( otherAtom )
          call ListInteger_destructor( atomForDefineAngle )
          call ListInteger_destructor( atomForDefineDihedral )
          call List_destructor( dihedralAngles )
          call Vector_destructor( vectorA )
          call Vector_destructor( vectorB )
          call Vector_destructor( auxVector )

       end if

    else

       call InternalCoordinates_exception( ERROR, "The chemical bonds haven't been defined",&
            "Class object InternalCoordinates in the getDihedralsAngles() function")

    end if

  end function InternalCoordinates_getDihedralAngles

  !>
  !! @brief Retorna el numero de coordenadas primitivas que han sido seleccionadas
  !<
  function InternalCoordinates_getNumberOfCoordinates(this) result(output)
    implicit none
    type(InternalCoordinates) :: this
    integer :: output

    output = 0
    this%numberOfBonds=0
    this%numberOfAnglesOfBond=0
    this%numberOfDihedrals=0

    if( allocated( this%connectionMatrixForBonds%values ) ) then

       this%numberOfBonds = size( this%connectionMatrixForBonds%values, dim=1 )
       output = this%numberOfBonds

       if( allocated( this%connectionMatrixForAngles%values ) ) then

          this%numberOfAnglesOfBond = size( this%connectionMatrixForAngles%values, dim=1 )
          output = output + this%numberOfAnglesOfBond

          if( allocated( this%connectionMatrixForDihedrals%values ) ) then
             this%numberOfDihedrals = size( this%connectionMatrixForDihedrals%values, dim=1 )
             output = output + this%numberOfDihedrals
          end if

       end if

    end if

  end function InternalCoordinates_getNumberOfCoordinates

  !>
  !! @brief Retorna una matriz auxiliar M, a traves de la cual se remueven los grados
  !!		de libertad externos del sistema molecular, de la transformacion de coordenadas
  !!		del gradiente y de la hessiana.
  !!
  !!
  !! La matriz auxiliar tiene la forma:
  !!\f[
  !! M =
  !!  \left(\begin{array}{ccccccccccc}
  !!   0&0& \cdots &&&&&&&&0 \\  \vdots&&&&&&&&&& \\  &&0&&&&&&&&\\  &&&0&&&&&&&\\  &&&&0&&&&&& \\  &&&&&1&
  !!	&&&& \\  &&&&&&1&&&& \\  &&&&&&&0&&& \\  &&&&&&&&1&& \\  &&&&&&&&&\ddots& \\  0&&&&&&&&&&1
  !! \end{array} \right)
  !! \f]
  !!
  !<
  function InternalCoordinates_getAuxiliaryMatrix( this ) result(output)
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    integer(8) :: numberOfVariables

    numberOfVariables = this%numberOfCenterOfOptimization*3_8


    call Matrix_constructor( output, numberOfVariables, numberOfVariables )
    call Matrix_setIdentity( output )

    !**
    !! Fija arbitrariamente las  coordenadas de los tres primeros atomos del sistemas molecular
    !! con el fin de remover los grados de libertad externos. El primer atomo se ubica en el origen del sistema
    !! de coordenadas; el segundo atomo se ubica  sobre el eje  z y el tercero en el plano xz. Para ello se
    !! se asigna el valor de cero a los siguientes elementos diagonales de la matriz auxiliar: M_11, M_22,
    !! M_33, M_44, M_55, M_88.
    !!
    !! Ver: Reveles and Koster, J Comp Chem, 25, 9, 2004, p1109-1116.
    !**
    output%values(1,1) = 0.0_8
    output%values(2,2) = 0.0_8
    output%values(3,3) = 0.0_8
    output%values(4,4) = 0.0_8
    output%values(5,5) = 0.0_8
    output%values(8,8) = 0.0_8

  end function InternalCoordinates_getAuxiliaryMatrix

  !>
  !! @brief Retorna la matriz \f$ U\f$ de vectores propios no-redundantes de  la matriz G
  !!
  !! La matriz G es  una matriz simetrica factorizable de la forma:
  !! \f[
  !!   {\bf G = PMP^T=(UR)}
  !!  \left(\begin{array}{cc}
  !!   \Lambda&0 \\  0&0
  !! \end{array} \right)
  !!\f]
  !! donde \f$ U \f$ esta asociada a los vectores propios correspondientes a valores propios \f$ \lambda_i > 0\f$.
  !! \f$ P \f$  y \f$ M  \f$, son la matriz de Wilson y la matriz auxiliar respectivamete
  !<
  function InternalCoordinates_getNonRedundantEigenvectors( this ) result( output)
    implicit none
    type( InternalCoordinates ):: this
    type( Matrix ) :: output

    type(Vector) :: eigenValues
    integer :: numberOfNonredundantCoordinates
    integer :: aux

    if (  allocated(this%symmetricGMatrix%values) ) then

       numberOfNonredundantCoordinates = size(this%symmetricGMatrix%values,dim=1)
       call Vector_constructor( eigenValues, numberOfNonredundantCoordinates )
       call Matrix_constructor( output, int(numberOfNonredundantCoordinates,8), int(numberOfNonredundantCoordinates,8) )

       call Matrix_eigen( this%symmetricGMatrix, eigenValues, output, SYMMETRIC  )

       !!Elimina los vectores propios redundantes de la matriz de valores propios
       aux = minloc( eigenValues%values,dim=1 )
       do  while ( eigenValues%values(aux) < 1.0D-7)
          call Matrix_removeColumn(output, aux)
          call Vector_removeElement(eigenValues, aux)
          aux = minloc( eigenValues%values, dim=1 )
       end do

       call Vector_destructor( eigenValues )
    else

       call InternalCoordinates_exception( ERROR, "Symmetric G matrix hasn't been defined",&
            "Class object InternalCoordinates in the getNonRedundatEigenvectors() function")

    end if

  end function InternalCoordinates_getNonRedundantEigenvectors


  !>
  !! @brief calcula los pesos de la coordenadas primitivas
  !!
  !<
  function InternalCoordinates_getWeightOfCoordinate( this ) result( output )
    implicit none
    type( InternalCoordinates ) :: this
    type(Vector) :: output


    integer  :: i
    integer :: ssize

    if  ( allocated( this%nonRedundantEigenvectors%values) ) then

       ssize = size(this%nonRedundantEigenvectors%values, dim=1 )
       call Vector_constructor( output, ssize )

       output%values= 0.0_8

       do i=1, ssize
          output%values(i)= sum( ( this%nonRedundantEigenvectors%values(i,:) )**2.0 )
       end do

    end if

  end function InternalCoordinates_getWeightOfCoordinate



  !>
  !! @brief Ajusta la matriz de conectividad para los enlaces quimicos seleccionados
  !<
  subroutine InternalCoordinates_setBonds( this, bonds )
    implicit none
    type(InternalCoordinates) :: this
    type(MatrixInteger) :: bonds

    if( size(bonds%values, dim=2) == 2 ) then

       this%connectionMatrixForBonds = bonds

    else
       call InternalCoordinates_exception( ERROR, "The size of matrices don't match",&
            "Class object InternalCoordinates in the setBonds() function")
    end if

  end subroutine InternalCoordinates_setBonds

  !>
  !! @brief Ajusta la matriz de conectividad para los angulos de enlace quimicos seleccionados
  !<
  subroutine InternalCoordinates_setAnglesOfBond( this, anglesOfBond )
    implicit none
    type(InternalCoordinates) :: this
    type(MatrixInteger) :: anglesOfBond

    if( size(anglesOfBond%values, dim=2) == 3 ) then

       this%connectionMatrixForAngles = anglesOfBond

    else
       call InternalCoordinates_exception( ERROR, "The size of matrices don't match",&
            "Class object InternalCoordinates in the setAnglesOfBonds() function")
    end if

  end subroutine InternalCoordinates_setAnglesOfBond


  !>
  !! @brief Ajusta la matriz de conectividad para los angulos de enlace quimicos seleccionados
  !<
  subroutine InternalCoordinates_setDihedraslAngles( this, dihedralAngles )
    implicit none
    type(InternalCoordinates) :: this
    type(MatrixInteger) :: dihedralAngles

    if( size(dihedralAngles%values, dim=2) == 4 ) then

       this%connectionMatrixForDihedrals = dihedralAngles

    else
       call InternalCoordinates_exception( ERROR, "The size of matrices don't match",&
            "Class object InternalCoordinates in the setAnglesOfBonds() function")
    end if

  end subroutine InternalCoordinates_setDihedraslAngles

  !>
  !! @brief Selecciona todos los enlaces, angulos y torsiones segun criterio de conectividad
  !<
  subroutine InternalCoordinates_obtainCoordinates( this )
    implicit none
    type(InternalCoordinates) :: this

    integer :: i


    !! Selecciona enlaces de acuerdo a criterio de conectividad
    this%connectionMatrixForBonds = InternalCoordinates_getBonds( this, CONTROL_instance%BOND_DISTANCE_FACTOR )

    ! if ( size(this%connectionMatrixForBonds%values,dim=1) >1 ) then

    !    write (6,"(T20,A24)") "CHEMICAL BONDS: AMSTRONG"
    !    write (6,"(T20,A23)") "======================="
    !    write (6,*) ""

    !    do i=1,size(this%distanceBondValue%values)
    !       write (6,"(T20,I5,A1,2I5,F15.8)") i,":",this%connectionMatrixForBonds%values(i,:), this%distanceBondValue%values(i) &
    !            * AMSTRONG
    !    end do
    ! end if

    !! Selecciona angulos de enlace de acuerdo a enlaces quimicos definidos
    if ( this%numberOfBonds >= 2) then

       this%connectionMatrixForAngles = InternalCoordinates_getAnglesOfBond(this, CONTROL_instance%BOND_ANGLE_THRESHOLD )

       ! if ( size(this%connectionMatrixForAngles%values,dim=1)>=1 .and. &
       !      sum(this%connectionMatrixForAngles%values(1,:))>0 .and. &
       !      sum(abs(this%angleOfBondValue%values)) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
       !    write (6,*) ""
       !    write (6,"(T20,A26)") "ANGLES OF VALENCE: DEGREES"
       !    write (6,"(T20,A25)") "========================="
       !    write (6,*) ""


       !    do i=1,size(this%angleOfBondValue%values)
       !       write (6,"(T20,I5,A1,3I5,F15.8)") i,":",this%connectionMatrixForAngles%values(i,:), this%angleOfBondValue%values(i)
       !    end do
       ! end if
    end if

    !! Selecciona angulos diedros de acuerdo a agulos de enlaces definidos
    if ( this%numberOfAnglesOfBond >= 2)  then

       this%connectionMatrixForDihedrals = &
            InternalCoordinates_getDihedralAngles(this, CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD )
       ! if ( size(this%connectionMatrixForDihedrals%values,dim=1)>=1 .and. &
       !      sum(this%connectionMatrixForDihedrals%values(1,:))>0 .and. &
       !      sum(abs(this%dihedralsAngleValue%values)) > CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then 
       !    write (6,*) ""
       !    write (6,"(T20,A27)") "ANGLES OF TORSION: DEGREES"
       !    write (6,"(T20,A26)") "=========================="
       !    write (6,*) ""


       !    do i=1,size(this%dihedralsAngleValue%values)
       !       write (6,"(T20,I5,A1,4I5,F15.8)") i,":",this%connectionMatrixForDihedrals%values(i,:), &
       !            this%dihedralsAngleValue%values(i)
       !    end do
       ! end if
    end if

  end subroutine InternalCoordinates_obtainCoordinates

  subroutine InternalCoordinates_primitiveCoordinateSelection( this )
    implicit none
    type(InternalCoordinates) :: this

    type(Matrix) :: auxiliaryMatrix
    !		type(Matrix) :: temp
    !		type(Matrix) :: auxiliaryMatrix
    !		type(Matrix) :: auxiliaryMatrix
    type(Vector) :: weightOfCoordinate
    integer :: i
    integer :: j
    integer :: k
    integer :: aux

    do k=1,3

       select case (k)

       case(1)
          this%connectionMatrixForBonds = InternalCoordinates_getBonds( this, CONTROL_instance%BOND_DISTANCE_FACTOR )

       case(2)

          if ( this%numberOfBonds >= 2) &
               this%connectionMatrixForAngles = InternalCoordinates_getAnglesOfBond(this, CONTROL_instance%BOND_ANGLE_THRESHOLD )
      
       case(3)

          if ( this%numberOfAnglesOfBond >= 2) &
               this%connectionMatrixForDihedrals = &
               InternalCoordinates_getDihedralAngles(this, CONTROL_instance%DIHEDRAL_ANGLE_THRESHOLD )
       end select


       !! Contruye la matriz de Wilson para los enlaces quimicos definidos
       this%wilsonMatrix = InternalCoordinates_getWilsonMatrix( this )

       call updateWeights()

       aux = minloc(weightOfCoordinate%values, dim=1)

       do while (  ( weightOfCoordinate%values(aux) < 0.3)  .and. &
            InternalCoordinates_getNumberOfCoordinates(this) > 2*k*this%numberOfCenterOfOptimization )

          call Vector_removeElement(weightOfCoordinate, aux)
          call InternalCoordinates_removeCoordinate(this, aux)
          call updateWeights()
          aux = minloc(weightOfCoordinate%values, dim=1)

       end do

    end do

    print *,"Numero de enlaces: ",this%numberOfBonds
    print *,"Numero de angulos de enlace: ",this%numberOfAnglesOfBond
    print *,"Numero de angulos diedros: ",this%numberOfDihedrals

    call Vector_destructor(weightOfCoordinate)
    call Matrix_destructor( auxiliaryMatrix )

  contains

    subroutine updateWeights()
      implicit none

      !! Obtiene la matriz auxiliar para eliminar  los grados de libertad externos.
      auxiliaryMatrix = InternalCoordinates_getAuxiliaryMatrix(this)

      !!
      !! Construye una  matriz G simetrica, a traves de la cual se realiza la
      !! separacion de coordenadas redundantes y no-redundantes
      !!
      aux=InternalCoordinates_getNumberOfCoordinates(this)
      call Matrix_constructor(this%symmetricGMatrix, int(aux,8),int(aux,8) )

      this%symmetricGMatrix%values = matmul( this%wilsonMatrix%values, &
           matmul( auxiliaryMatrix%values, transpose( this%wilsonMatrix%values ) ) )

      !! Se diagonaliza matriz G para obtener vectores propios redundantes R y no redundantes U
      this%nonRedundantEigenvectors= InternalCoordinates_getNonRedundantEigenvectors( this )

      !! Obtiene los pesos de la coordenadas primitivas
      weightOfCoordinate=InternalCoordinates_getWeightOfCoordinate(this)

    end subroutine updateWeights


  end subroutine InternalCoordinates_primitiveCoordinateSelection

  !>
  !! @brief Remueve la coordena primitiva indicada como argumento, si la matrix
  !! 		de wilson ya ha sido calculada tambien elimina la fila correspondiente.
  !!
  !! @param indexOfCoordinate Indice de la coordena primitiva que se elimina
  !<
  subroutine InternalCoordinates_removeCoordinate(this, indexOfCoordinate  )
    implicit none
    type(InternalCoordinates) :: this
    integer indexOfCoordinate

    integer :: numberOfPrimitivesCoordinates
    integer :: auxIndex

    numberOfPrimitivesCoordinates = InternalCoordinates_getNumberOfCoordinates( this )

    if ( indexOfCoordinate <= this%numberOfBonds) then

       this%numberOfBonds = this%numberOfBonds - 1
       this%numberOfPrimitivesCoordinates = this%numberOfPrimitivesCoordinates - 1

       call MatrixInteger_removeRow(this%connectionMatrixForBonds, indexOfCoordinate )
       call Vector_removeElement(this%distanceBondValue, indexOfCoordinate)

    else &
         if ( indexOfCoordinate <= this%numberOfBonds + this%numberOfAnglesOfBond ) then

       this%numberOfAnglesOfBond = this%numberOfAnglesOfBond - 1
       this%numberOfPrimitivesCoordinates = this%numberOfPrimitivesCoordinates - 1
       auxIndex= indexOfCoordinate - this%numberOfBonds

       call MatrixInteger_removeRow(this%connectionMatrixForAngles, auxIndex )
       call Vector_removeElement(this%angleOfBondValue, auxIndex )

    else &
         if ( indexOfCoordinate <= numberOfPrimitivesCoordinates ) then

       this%numberOfDihedrals = this%numberOfDihedrals - 1
       this%numberOfPrimitivesCoordinates = this%numberOfPrimitivesCoordinates - 1
       auxIndex= indexOfCoordinate - this%numberOfBonds  - this%numberOfAnglesOfBond

       call MatrixInteger_removeRow(this%connectionMatrixForDihedrals, auxIndex )
       call Vector_removeElement(this%dihedralsAngleValue, auxIndex )

    end if

  end subroutine InternalCoordinates_removeCoordinate



  !>
  !! @brief 	Retorna una aproximacion de los elementos diagonales  de la matriz hessiana, basada en las
  !!		formulas empiricas de Fischer y Almlof para calculo de constantes
  !!		de fuerza. J Phys Chem. 96, 24 1992, 9768-9774
  !!
  !<
  function InternalCoordinates_getStartHessian( this )  result( output )
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    character(10), allocatable :: labelOfCenters(:)
    character(10) :: symbolOfAtom
    character(10) :: symbolOfOtherAtom
    real(8) :: covalentRadiusValues(2)
    real(8) :: vectorAB(3)
    real(8) :: vectorCD(3)
    real(8) :: distanceAB
    real(8) :: distanceCD
    integer :: auxIndex
    integer :: valence
    integer :: indexOfTerminalAtom
    integer :: i
    integer :: j

    allocate( labelOfCenters( this%numberOfCenterOfOptimization ) )
    labelOfCenters = ParticleManager_getLabelsOfCentersOfOptimization()

    call Matrix_constructor(output, int(this%numberOfPrimitivesCoordinates,8), int(this%numberOfPrimitivesCoordinates,8), 0.0_8 )

    !! Calculo de elementos fuera de la diagonal
    do i=1, this%numberOfPrimitivesCoordinates

       !! Calculo de H_streching
       if ( i <= this%numberOfBonds) then


          symbolOfAtom=trim( labelOfCenters( this%connectionMatrixForBonds%values( i, 1 ) ) )
          covalentRadiusValues(1) = Map_getValue( this%covalentRadius, symbolOfAtom )
          symbolOfOtherAtom=trim( labelOfCenters( this%connectionMatrixForBonds%values( i, 2 ) ) )
          covalentRadiusValues(2) = Map_getValue( this%covalentRadius, symbolOfOtherAtom )

          output%values(i,i) = 0.3601_8 * exp( -1.944_8*( ( this%distanceBondValue%values(i)* AMSTRONG ) &
               -  sum( covalentRadiusValues ) ) )


          !! Calculo de H_bend
       else &
            if ( i <= this%numberOfBonds + this%numberOfAnglesOfBond ) then

          auxIndex= i - this%numberOfBonds

          !!
          !! Calculo de distancias entre atomos terminales y central
          !!
          vectorAB=	this%cartesianCoordinates%values( this%connectionMatrixForAngles%values(auxIndex,1),: ) - &
               this%cartesianCoordinates%values( this%connectionMatrixForAngles%values(auxIndex,2),: )
          distanceAB = dsqrt(sum(vectorAB**2.0) ) * AMSTRONG

          vectorCD=	this%cartesianCoordinates%values( this%connectionMatrixForAngles%values(auxIndex,3),: ) - &
               this%cartesianCoordinates%values( this%connectionMatrixForAngles%values(auxIndex,2),: )
          distanceCD = dsqrt(sum(vectorCD**2.0) ) * AMSTRONG

          !!
          !! Calculo de sumas de radios covalentes entre atomos terminales y central
          !!
          symbolOfAtom=trim( labelOfCenters( this%connectionMatrixForAngles%values( auxIndex, 1 ) ) )
          covalentRadiusValues(1) = Map_getValue( this%covalentRadius, symbolOfAtom )
          symbolOfOtherAtom=trim( labelOfCenters( this%connectionMatrixForAngles%values( auxIndex, 2 ) ) )
          covalentRadiusValues(2) = Map_getValue( this%covalentRadius, symbolOfOtherAtom )
          covalentRadiusValues(1)= covalentRadiusValues(1) + covalentRadiusValues(2)
          symbolOfAtom=trim( labelOfCenters( this%connectionMatrixForAngles%values( auxIndex, 3 ) ) )
          covalentRadiusValues(2)=covalentRadiusValues(2)+ Map_getValue( this%covalentRadius, symbolOfAtom )

          output%values(i,i) = 0.089_8 + (0.11_8/( ( covalentRadiusValues(1) * covalentRadiusValues(2) )**(-0.42_8) ) ) &
               * exp( -0.44_8 * ( distanceAB + distanceCD - sum(covalentRadiusValues) ) )

          !! Calculo de H_torsional
       else &
            if ( i <= this%numberOfPrimitivesCoordinates ) then

          auxIndex= i - this%numberOfBonds - this%numberOfAnglesOfBond

          valence = -2
          indexOfTerminalAtom=this%connectionMatrixForDihedrals%values(auxIndex,1)
          do j=1, this%numberOfBonds
             if ( this%connectionMatrixForBonds%values(j,1)==indexOfTerminalAtom &
                  .or. this%connectionMatrixForBonds%values(j,2)==indexOfTerminalAtom ) valence= valence+1
          end do
          indexOfTerminalAtom=this%connectionMatrixForDihedrals%values(auxIndex,4)
          do j=1, this%numberOfBonds
             if ( this%connectionMatrixForBonds%values(j,1)==indexOfTerminalAtom &
                  .or. this%connectionMatrixForBonds%values(j,2)==indexOfTerminalAtom ) valence= valence+1
          end do

          !!
          !! Calculo de distancias entre atomos terminales
          !!
          vectorAB=	this%cartesianCoordinates%values( this%connectionMatrixForDihedrals%values(auxIndex,1),: ) - &
               this%cartesianCoordinates%values( this%connectionMatrixForDihedrals%values(auxIndex,4),: )
          distanceAB = dsqrt(sum(vectorAB**2.0) ) * AMSTRONG

          !!
          !! Calculo de sumas de radios covalentes entre atomos terminales y central
          !!
          symbolOfAtom=trim( labelOfCenters( this%connectionMatrixForDihedrals%values( auxIndex, 1 ) ) )
          covalentRadiusValues(1) = Map_getValue( this%covalentRadius, symbolOfAtom )
          symbolOfOtherAtom=trim( labelOfCenters( this%connectionMatrixForDihedrals%values( auxIndex, 4 ) ) )
          covalentRadiusValues(2) = Map_getValue( this%covalentRadius, symbolOfOtherAtom )

          output%values(i,i)= 0.0015_8 + ( (14.0_8 * valence**0.57_8 )/( ( distanceAB * ( sum(covalentRadiusValues ) ) )**4.0 ) ) &
               * exp( -2.85_8 * (distanceAB- sum(covalentRadiusValues ) ) )

       end if

    end do

    deallocate(labelOfCenters)

  end function InternalCoordinates_getStartHessian

  !>
  !! @brief Retorna la matriz B de Wilson de dimension \f$m\f$x\f$3N, m=3N-6(5), N:\f$ numero de centros de optimizacion.
  !!
  !!\f[
  !!	\mathbf{B}=\mathbf{U}^T\mathbf{P}
  !! \f]
  !!
  !! donde \f$ \mathbf{U} \f$ es la matrix de vectores no redundantes, tras la diagonalizacion de la matriz G.
  !<
  function InternalCoordinates_getBWilsonMatrix( this ) result(output)
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    call Matrix_constructor( output , int(size( this%nonRedundantEigenvectors%values, dim=2),8) , &
         int( this%numberOfCenterOfOptimization*3,8 ) )

    output%values= matmul( transpose(this%nonRedundantEigenvectors%values), this%wilsonMatrix%values )

  end function InternalCoordinates_getBWilsonMatrix


  !>
  !! @brief Retorna la inversa  generalizada de la transpuesta de la matriz B de Wilson.
  !!
  !!\f[
  !!	\mathbf{\~B^-}=(\mathbf{BM\~B})^{-1}\mathbf{BM}
  !! \f]
  !!
  !<
  function InternalCoordinates_getInverseOfTransposeBWilsonMatrix( this ) result(output)
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: otherAuxiliaryMatrix
    type(Matrix) :: eigenVectors

    auxiliaryMatrix = InternalCoordinates_getAuxiliaryMatrix(this)

    call Matrix_constructor( output , int(size( this%nonRedundantEigenvectors%values, dim=2),8) , &
         int( this%numberOfCenterOfOptimization*3,8 ) )

    call Matrix_constructor( otherAuxiliaryMatrix , int(size( this%nonRedundantEigenvectors%values, dim=2),8) , &
         int(size( this%nonRedundantEigenvectors%values, dim=2),8) )

    otherAuxiliaryMatrix%values = matmul( this%BWilsonMatrix%values, matmul( auxiliaryMatrix%values, &
         transpose(this%BWilsonMatrix%values) ) )

    eigenVectors = Matrix_inverse( otherAuxiliaryMatrix )
    output%values =matmul( eigenVectors%values, matmul( this%BWilsonMatrix%values, auxiliaryMatrix%values ) )

    call Matrix_destructor( auxiliaryMatrix )
    call Matrix_destructor( eigenVectors )

  end function InternalCoordinates_getInverseOfTransposeBWilsonMatrix


  !>
  !! @brief Retorna la inversa  de la matriz B de Wilson.
  !!
  !!\f[
  !!	\mathbf{\~B}=(\mathbf{\~BM'B})^{-1}\mathbf{\~BM}
  !! \f]
  !!
  !<
  function InternalCoordinates_getInverseOfBWilsonMatrix( this ) result(output)
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    type(Matrix) :: auxiliaryMatrix
    type(Matrix) :: otherAuxiliaryMatrix
    type(Matrix) :: eigenVectors

    auxiliaryMatrix = InternalCoordinates_getAuxiliaryMatrix(this)

    call Matrix_constructor( output , int(size( this%nonRedundantEigenvectors%values, dim=2),8) , &
         int( this%numberOfCenterOfOptimization*3,8 ) )

    call Matrix_constructor( otherAuxiliaryMatrix , int(size( this%nonRedundantEigenvectors%values, dim=2),8) , &
         int(size( this%nonRedundantEigenvectors%values, dim=2),8) )

    otherAuxiliaryMatrix%values = matmul( this%BWilsonMatrix%values, matmul( auxiliaryMatrix%values, &
         transpose(this%BWilsonMatrix%values) ) )

    eigenVectors = Matrix_inverse( otherAuxiliaryMatrix )
    output%values =matmul( eigenVectors%values, matmul( this%BWilsonMatrix%values, auxiliaryMatrix%values ) )

    call Matrix_destructor( auxiliaryMatrix )
    call Matrix_destructor( eigenVectors )

  end function InternalCoordinates_getInverseOfBWilsonMatrix


  !>
  !! @brief Retorna la matriz de Wilson para las m coordenadas primitivas seleccionadas \f$  P_{ij}=\frac{\partial p_i}{\partial x_j} \f$
  !!
  !!\f[
  !! P =
  !!  \left(\begin{array}{ccc}
  !!  \frac{\partial p_1}{\partial x_1} & \cdots & \frac{\partial p_1}{\partial x_n} \\ \vdots &  \ddots &  \\  \frac{\partial p_m}{\partial x_1} &  &
  !!  \frac{\partial p_m}{\partial x_n}
  !! \end{array} \right)
  !! \f]
  !!
  !! @todo Simplificar las expresiones para c\'alculo de derivadas de \'angulos
  !!             diedros respecto a coordenadas cartesianas
  !! @todo Verificar que el signo de la derivada de angulos  diedros es correcto
  !!
  !! @warning Las expresiones para derivadas de angulos diedros no han sido simplificadas.
  !<
  function InternalCoordinates_getWilsonMatrix( this ) result(output)
    implicit none
    type(InternalCoordinates) :: this
    type(Matrix) :: output

    integer :: i
    integer :: j
    integer :: k
    integer :: atom
    integer :: component
    integer :: aux
    real(8) :: originNormOfAtom
    real(8) :: originNormOfOtherAtom
    real(8) :: dotProduct
    real(8) :: ssign
    type(Vector) :: originAtomOne
    type(Vector) :: originAtomTwo
    type(Vector) :: originAtomThree
    type(Vector) :: originAtomFour

    this%numberOfPrimitivesCoordinates = InternalCoordinates_getNumberOfCoordinates(this)

    call Matrix_constructor( output, int(this%numberOfPrimitivesCoordinates,8), &
         int(this%numberOfCenterOfOptimization*3,8), 0.0_8 )

    !!
    !! Calcula los elementos  de  la matrix de wilson para enlaces quimicos definidos.
    !! P_ij = dp_i/ dx_j
    !!
    if ( allocated(this%connectionMatrixForBonds%values) ) then

       do i=1, this%numberOfBonds
          j=0
          do atom =1, this%numberOfCenterOfOptimization
             do component =1,3
                j=j+1

                aux = Math_kroneckerDelta(  this%connectionMatrixForBonds%values(i,1), atom ) &
                     - Math_kroneckerDelta(  this%connectionMatrixForBonds%values(i,2), atom )

                if ( abs(aux) > 0 ) then

                   output%values(i,j) = aux * &
                        ( this%cartesianCoordinates%values( this%connectionMatrixForBonds%values(i,1), component ) &
                        - this%cartesianCoordinates%values( this%connectionMatrixForBonds%values(i,2), component ) ) &
                        / this%distanceBondValue%values(i)

                end if

             end do
          end do
       end do

    end if

    !!
    !! Calcula los elementos  de  la matrix de wilson para angulos de enlace%
    !! P_ij = dp_i/ dx_j, donde p_i = \theta
    !!
    if ( allocated(this%connectionMatrixForAngles%values) ) then

       call Vector_constructor( originAtomOne, 3)
       call Vector_constructor( originAtomTwo, 3)

       do i = 1, this%numberOfAnglesOfBond
          j=0
          k= this%numberOfBonds+i

          originAtomOne%values = this%cartesianCoordinates%values(this%connectionMatrixForAngles%values(i,1), : ) &
               - this%cartesianCoordinates%values(this%connectionMatrixForAngles%values(i,2), : )

          originNormOfAtom = Vector_norm( originAtomOne )

          originAtomTwo%values = 	this%cartesianCoordinates%values(this%connectionMatrixForAngles%values(i,3), : ) &
               - this%cartesianCoordinates%values(this%connectionMatrixForAngles%values(i,2), : )

          originNormOfOtherAtom = Vector_norm( originAtomTwo )

          dotProduct = dot_product(originAtomOne%values, originAtomTwo%values)


          do atom =1, this%numberOfCenterOfOptimization
             do component =1,3
                j=j+1

                if ( Math_kroneckerDelta(  this%connectionMatrixForAngles%values(i,1), atom ) > 0 )  then

                   output%values(k,j) =  ( ( originAtomOne%values(component) * dotProduct ) &
                        / ( ( originNormOfAtom**3.0_8 ) * originNormOfOtherAtom )  &
                        - originAtomTwo%values(component) / (originNormOfAtom*originNormOfOtherAtom) ) &
                        / dsqrt( 1.0_8 - (dotProduct**2.0_8) /( ( originNormOfAtom * originNormOfOtherAtom )**2.0_8 ) )

                else &
                     if ( Math_kroneckerDelta(  this%connectionMatrixForAngles%values(i,3), atom ) > 0 )  then

                   output%values(k,j) = ( ( originAtomTwo%values(component) * dotProduct ) &
                        / ( originNormOfAtom * ( originNormOfOtherAtom**3.0_8 )  )  &
                        - originAtomOne%values(component) / ( originNormOfAtom*originNormOfOtherAtom ) ) &
                        / dsqrt( 1.0_8 - ( dotProduct**2.0_8 ) /( ( originNormOfAtom * originNormOfOtherAtom ) ** 2.0_8 ) )

                else &
                     if  ( Math_kroneckerDelta(  this%connectionMatrixForAngles%values(i,2), atom ) > 0 ) then

                   output%values(k,j) = &
                        ( ( originAtomOne%values(component) + originAtomTwo%values(component) ) &
                        / ( originNormOfAtom * originNormOfOtherAtom  ) &
                        - ( ( originAtomTwo%values(component)*originNormOfAtom**2.0_8 &
                        + originAtomOne%values( component)  * originNormOfOtherAtom**2.0_8 ) &
                        * dotProduct ) &
                        / (  (originNormOfAtom * originNormOfOtherAtom )**3.0_8  ) ) &
                        / dsqrt( 1.0_8 - ( dotProduct**2.0_8) / ( ( originNormOfAtom * originNormOfOtherAtom )**2.0_8 ) )

                end if

             end do

          end do
       end do

       call Vector_destructor( originAtomOne )
       call Vector_destructor( originAtomTwo )

    end if


    !!
    !! Calculo de elementos  de  la matrix de wilson para angulos diedros.
    !! P_ij = dp_i/ dx_j, donde p_i = \gamma
    !!
    if ( allocated(this%connectionMatrixForDihedrals%values) ) then

       call Vector_constructor( originAtomOne, 3)
       call Vector_constructor( originAtomTwo, 3)
       call Vector_constructor( originAtomThree, 3)
       call Vector_constructor( originAtomFour, 3)

       do i = 1, this%numberOfDihedrals
          j=0
          k= this%numberOfBonds + this%numberOfAnglesOfBond + i
          originAtomOne%values = 	this%cartesianCoordinates%values(this%connectionMatrixForDihedrals%values(i,1), : )
          originAtomTwo%values = 	this%cartesianCoordinates%values(this%connectionMatrixForDihedrals%values(i,2), : )
          originAtomThree%values = 	this%cartesianCoordinates%values(this%connectionMatrixForDihedrals%values(i,3), : )
          originAtomFour%values = 	this%cartesianCoordinates%values(this%connectionMatrixForDihedrals%values(i,4), : )
          ssign = sign( 1.0_8, this%dihedralsAngleValue%values(i) )

          do atom =1, this%numberOfCenterOfOptimization
             do component =1,3
                j=j+1

                if ( Math_kroneckerDelta(  this%connectionMatrixForDihedrals%values(i,1), atom ) > 0 )  then

                   if( Math_kroneckerDelta(component,1) > 0 ) then


                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(originAtomThree%values(2)- &
                           originAtomTwo%values(2))+((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(-originAtomThree%values(3)+originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(originAtomThree%values(2)- &
                           originAtomTwo%values(2))+2.0_8*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomOne%values(3)-originAtomTwo%values(3))- &
                           (originAtomOne%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(-originAtomThree%values(3)+originAtomTwo%values(3))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   else &
                        if( Math_kroneckerDelta(component,2) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((((-originAtomThree%values(1)+originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-originAtomThree%values(1)+originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))


                   else &
                        if( Math_kroneckerDelta(component,3) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((-((2.0_8*(originAtomThree%values(1)-originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(-originAtomThree%values(2)+originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8)+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))+(-originAtomThree%values(2)+originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   end if


                else &
                     if ( Math_kroneckerDelta(  this%connectionMatrixForDihedrals%values(i,2), atom ) > 0 )  then


                   if( Math_kroneckerDelta(component,1) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((-((2.0_8*(originAtomOne%values(2)-originAtomThree%values(2))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-originAtomOne%values(3)+originAtomThree%values(3))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8)+ &
                           ((originAtomOne%values(2)-originAtomThree%values(2))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (originAtomFour%values(2)-originAtomThree%values(2))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (-originAtomOne%values(3)+originAtomThree%values(3))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-originAtomFour%values(3)+originAtomThree%values(3))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           ((2.0_8*(originAtomFour%values(2)-originAtomThree%values(2))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-originAtomFour%values(3)+originAtomThree%values(3))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))


                   else &
                        if( Math_kroneckerDelta(component,2) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((-((2.0_8*(-originAtomOne%values(1)+originAtomThree%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(originAtomOne%values(3)-originAtomThree%values(3))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8)+ &
                           ((-originAtomOne%values(1)+originAtomThree%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (-originAtomFour%values(1)+originAtomThree%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (originAtomOne%values(3)-originAtomThree%values(3))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (originAtomFour%values(3)-originAtomThree%values(3))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           ((2.0_8*(-originAtomFour%values(1)+originAtomThree%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(originAtomFour%values(3)-originAtomThree%values(3))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   else &
                        if( Math_kroneckerDelta(component,3) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((-((2.0_8*(originAtomOne%values(1)-originAtomThree%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(-originAtomOne%values(2)+originAtomThree%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8)+ &
                           ((originAtomOne%values(1)-originAtomThree%values(1))*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(1)-originAtomThree%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-originAtomOne%values(2)+originAtomThree%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-originAtomFour%values(2)+originAtomThree%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           ((2.0_8*(originAtomFour%values(1)-originAtomThree%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(-originAtomFour%values(2)+originAtomThree%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   end if

                else &
                     if  ( Math_kroneckerDelta(  this%connectionMatrixForDihedrals%values(i,3), atom ) > 0 ) then

                   if( Math_kroneckerDelta(component,1) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -( &
                           (-(((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(-originAtomFour%values(2)+ &
                           originAtomTwo%values(2))+2.0_8*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(originAtomFour%values(3)-originAtomTwo%values(3))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))+ &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(-originAtomFour%values(2)+ &
                           originAtomTwo%values(2))+(-((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(2)-originAtomTwo%values(2)))+ &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(2)- &
                           originAtomTwo%values(2)))*(-originAtomOne%values(2)+originAtomTwo%values(2))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)-originAtomTwo%values(3))- &
                           (originAtomOne%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(-originAtomOne%values(2)+ &
                           originAtomTwo%values(2))+2.0_8*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomOne%values(3)-originAtomTwo%values(3))- &
                           (originAtomOne%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(originAtomOne%values(3)-originAtomTwo%values(3))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   else &
                        if( Math_kroneckerDelta(component,2) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -( &
                           (-(((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(-originAtomFour%values(3)+ &
                           originAtomTwo%values(3))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))+ &
                           ((originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(-originAtomFour%values(3)+ &
                           originAtomTwo%values(3))+(-((originAtomThree%values(2)-originAtomTwo%values(2))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3)))+ &
                           (originAtomFour%values(2)-originAtomTwo%values(2))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(-originAtomOne%values(3)+originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(-originAtomOne%values(3)+ &
                           originAtomTwo%values(3))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))


                   else &
                        if( Math_kroneckerDelta(component,3) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((-((2.0_8*(-originAtomOne%values(1)+originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8)+ &
                           ((-originAtomOne%values(1)+originAtomTwo%values(1))*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))+(-originAtomFour%values(1)+originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           ((2.0_8*(-originAtomFour%values(1)+originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))


                   end if

                else &
                     if ( Math_kroneckerDelta(  this%connectionMatrixForDihedrals%values(i,4), atom ) > 0 )  then

                   if( Math_kroneckerDelta(component,1) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -( &
                           (-(((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(originAtomThree%values(2)- &
                           originAtomTwo%values(2))+2.0_8*((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomFour%values(3)-originAtomTwo%values(3))- &
                           (originAtomFour%values(1)-originAtomTwo%values(1))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))*(-originAtomThree%values(3)+originAtomTwo%values(3))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))+ &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))*(originAtomThree%values(2)- &
                           originAtomTwo%values(2))+((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomOne%values(3)-originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(-originAtomThree%values(3)+ &
                           originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   else &
                        if( Math_kroneckerDelta(component,2) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -( &
                           (-(((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           (2.0_8*(-originAtomThree%values(1)+originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           2.0_8*(-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))+ &
                           ((-originAtomThree%values(1)+originAtomTwo%values(1))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))*(originAtomThree%values(3)- &
                           originAtomTwo%values(3)))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   else &
                        if( Math_kroneckerDelta(component,3) > 0 ) then

                      !< @warning La siguiente expesion  es provisional, debe ser simplificada
                      output%values(k,j)= &
                           -((((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))* &
                           (originAtomOne%values(3)-originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-originAtomThree%values(2)+originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))/ &
                           (dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))- &
                           ((2.0_8*(originAtomThree%values(1)-originAtomTwo%values(1))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           2.0_8*(-originAtomThree%values(2)+originAtomTwo%values(2))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))))/ &
                           (2.0_8*((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)**1.5_8* &
                           dsqrt((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)))/ &
                           dsqrt( &
                           1-((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))* &
                           (-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))* &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3))))**2.0_8/ &
                           (((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(2)- &
                           originAtomTwo%values(2)))+(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3))-(originAtomFour%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomFour%values(3)- &
                           originAtomTwo%values(3)))+(originAtomFour%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8)* &
                           ((-((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(2)- &
                           originAtomTwo%values(2)))+(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(2)-originAtomTwo%values(2)))**2.0_8+ &
                           ((originAtomThree%values(1)-originAtomTwo%values(1))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3))-(originAtomOne%values(1)-originAtomTwo%values(1))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8+ &
                           (-((originAtomThree%values(2)-originAtomTwo%values(2))*(originAtomOne%values(3)- &
                           originAtomTwo%values(3)))+(originAtomOne%values(2)-originAtomTwo%values(2))* &
                           (originAtomThree%values(3)-originAtomTwo%values(3)))**2.0_8))))

                   end if

                end if

                output%values(k,j)=ssign*output%values(k,j)

             end do

          end do
       end do


       call Vector_destructor( originAtomOne )
       call Vector_destructor( originAtomTwo )
       call Vector_destructor( originAtomThree )
       call Vector_destructor( originAtomFour )

    end if

  end function InternalCoordinates_getWilsonMatrix


  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine InternalCoordinates_exception( typeMessage, description, debugDescription)
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

  end subroutine InternalCoordinates_exception

end module InternalCoordinates_
