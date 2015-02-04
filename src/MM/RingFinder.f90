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
!!        This module search rings in the system using the <b> Hanser's algorithm </b>
!! @note Hanser, Th.; Jauffret, Ph.; Kaufmann, G., 
!!        <b>A New Algorithm for Exhaustive Ring Perception in a Molecular Graph</b>,
!!        J. Chem. Inf. Comput. Sci., 36, 1146--1152, 1996
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
module RingFinder_
  use MatrixInteger_
  use Vector_
  use MMCommons_
  use ListInteger_
  use Exception_
  implicit none

  public :: &
       RingFinder_getRings

contains

  !>
  !! @brief This routine searchs for the rings in the system
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this INTEGER ARRAY with the information of the edges in the pruned graph 
  !! @param [in] connectivityMatrix INTEGER ARRAY with the connectivity in the pruned graph 
  !! @param [in] numberOfRings INTEGER number of rings in the system 
  !! @return [out] rings INTEGER ARRAY with the rings found
  !! @see mmcommons_::mmcommons_pruninggraph
  !! @see rings_::rings_constructor
  subroutine RingFinder_getRings( this, connectivityMatrix, numberOfRings, rings)
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    type(MatrixInteger), intent(in) :: connectivityMatrix
    integer, intent(in) :: numberOfRings
    type(MatrixInteger), intent(out), allocatable :: rings(:)
    ! type(MatrixInteger), intent(out), allocatable :: rings3D(:)
    type(MatrixInteger), allocatable :: auxEdges(:)
    integer :: connectivitySize
    integer :: edgesSize
    integer :: i, ringRow, s
    integer :: j, k, numberOfColumns, numberOfRows
    ! integer :: l, m
    type(MatrixInteger) :: vertices
    logical :: isCycle
    logical :: isPath
    integer :: totalRings
    ! type(ListInteger) :: aux3DRings
    ! type(ListInteger) :: VertexRepeated
    ! integer :: size3Dring
    ! integer :: sizeRing
    ! integer :: numberOfRepeated
    type(MatrixInteger), allocatable :: edgesMatrix(:)

    !! Imprime valores de conectividad luego de cortar, borrar luego
    edgesSize = size(this)
    connectivitySize = size(connectivityMatrix%values)/2
    ! write (*,"(T10,A)") " "
    ! write (*,"(T10,A)") " Connectivity Matrix in Pruned Graph"
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! write (*,"(T20,2A)") "Idx ", "  | Connectivity"
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! do i=1, connectivitySize
    !    write (*,"(T10,I,I)") connectivityMatrix%values(i,1)
    ! end do
    ! !!******************************************************************************

    ! !! Imprime valores de ejes luego de cortar, borrar luego
    ! write (*,"(T10,A)") " Edges and Bonds in Pruned Graph " 
    ! do i=1,edgesSize
    !    write (6,"(T20,I5,A1,2I5)") i,": ", this%values(i,:) !! borrar luego, imprime ejes
    ! end do
    ! !!******************************************************************************

    call MatrixInteger_constructor( vertices, connectivitySize, 1 )

    do i=1, connectivitySize
       vertices%values(i,1) = connectivityMatrix%values(i,1)
    end do

    !!    Imprime valores de vertices del pruned graph, borrar luego
    ! write (*,"(T10,A)") " "
    ! write (*,"(T10,A)") " Atoms in Pruned Graph"
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! do i=1, connectivitySize
    !    write (*,"(T10,I)") vertices%values(i,:)
    ! end do
    !!    ******************************************************************************
    !!    Imprime los ejes, borrar luego
    ! write (*,"(T10,A)") " "
    ! write (*,"(T10,A)") " Edges in Pruned Graph"
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! write (*,"(T10,A,I)") " Tamanio =>", edgesSize
    ! do i=1,edgesSize
    !    write (*,"(T10,2I)") this(i)%values(1,:)
    ! end do

    allocate( auxEdges( edgesSize ))
    do i=1,edgesSize
       numberOfColumns = size(this(i)%values)
       call MatrixInteger_constructor( auxEdges(i), 1, numberOfColumns )
    end do

    do i=1,edgesSize
       numberOfColumns = size(this(i)%values)
       do j=1,numberOfColumns
          auxEdges(i)%values(1,j)=this(i)%values(1,j)
       end do
    end do

    allocate( rings( numberOfRings ))
    ringRow=1
    do i=1,size(vertices%values)
       if(size(auxEdges)>0) then
          call RingFinder_removeVertex(auxEdges,vertices%values(i,1))
          numberOfRows=size(auxEdges)
          ! write(*,"(T20,A,I)") "Numero de ejes luego del splice: ", numberOfRows
          j=1
          do while (j<=numberOfRows)
             ! write(*,"(T20,A,I)") "Evaluando la fila ", j
             ! write(*,"(T20,A,I,A)") "Quedan ", numberOfRows, " ejes."
             numberOfColumns = size(auxEdges(j)%values)
             isPath=RingFinder_isPath(auxEdges,j)
             if(.not.isPath) then
                ! write(*,"(T20,A)") "Evaluando el path "
                call MMCommons_removeEdge( auxEdges, j, numberOfRows)
                j=j-1
                numberOfRows=numberOfRows-1
             else if(numberOfColumns>=9) then
                ! write(*,"(T20,A)") "Evaluando el tamanio del eje "
                call MMCommons_removeEdge( auxEdges, j, numberOfRows)
                j=j-1
                numberOfRows=numberOfRows-1
                numberOfColumns=0
             else
                ! write(*,"(T20,A)") "Evaluando si hay ciclos "
                numberOfColumns = size(auxEdges(j)%values)
                ! write(*,"(T20,A,I,A)") "eje con ", numberOfColumns, " vertices"
                ! write(*,"(T20,<numberOfColumns>I)") auxEdges(j)%values(1,:)
                isCycle=RingFinder_isCycle(auxEdges,j)
                ! write(*,"(T20,A)") "ya evalue si hay ciclos "
                if(isCycle) then
                   ! write(*,"(T20,A,I)") " Encontre el anillo numero: ", ringRow
                   numberOfColumns = size(auxEdges(j)%values)
                   call MatrixInteger_constructor( rings(ringRow), 1, numberOfColumns )
                   do k=1,numberOfColumns
                      rings(ringRow)%values(1,k)=auxEdges(j)%values(1,k)
                   end do
                   call MMCommons_removeEdge( auxEdges, j, numberOfRows)
                   ringRow=ringRow+1
                   numberOfRows=numberOfRows-1
                   j=j-1
                end if
                ! write(*,"(T20,A)") " Salgo del IF de anillos "
             end if
             j=j+1
             ! write(*,"(T20,A,I)") "Sigue la fila ", j
             ! write(*,"(T20,A,I,A)") "Sobreviven ", numberOfRows, " ejes."
             ! write (*,"(T10,A)") " Anillos hallados "
             ! write (*,"(T10,A)") "--------------------------------------------"
             ! do s=1,ringRow-1
             !    numberOfColumns = size(rings(s)%values)
             !    write (*,"(T10,<numberOfColumns>I)") rings(s)%values(1,:)
             ! end do
          end do
          ! !! Imprime los anillos encontrados, borrar luego
          ! write (*,"(T10,A)") " Ejes luego de evaluar paths, tamanio y anillos "
          ! write (*,"(T10,A)") "--------------------------------------------"
          ! do s=1,numberOfRows
          !    numberOfColumns = size(auxEdges(s)%values)
          !    write (*,"(T10,<numberOfColumns>I)") auxEdges(s)%values(1,:)
          ! end do
!!!!**********************************************
       end if
    end do

    totalRings = ringRow-1 
    ! if(totalRings==numberOfRings) then
    !    size3Dring = size(rings(1)%values) - 1
    !    call ListInteger_constructor( aux3DRings, ssize=-1 )
    !    call ListInteger_constructor( VertexRepeated, ssize=-1 )
    !    do i=1,size3Dring
    !       call ListInteger_push_back(aux3DRings, rings(1)%values(i))
    !    end do
    !    do i=2,totalRings
    !       sizeRing = size(rings(i)%values) - 1
    !       do j=1,size3Dring
    !          do k=1,sizeRing
    !             if(aux3DRings%data(j)==rings(i)%values(k)) then
    !                call ListInteger_push_back(VertexRepeated, rings(i)%values(k))
    !             endif
    !          end do
    !       end do
    !       numberOfRepeated = ListInteger_size(VertexRepeated)
    !       if(numberOfRepeated >= 2) then
    !          do l=1,sizeRing
    !             do m=1,numberOfRepeated
    !                if(rings(i)%values(l) /= VertexRepeated%data(m)) then
    !                   !! ojo seguir aqui
    !                end if
    !             end do
    !          end do
    !       end if
    !    end do
    if(totalRings<numberOfRings) then
       do i=ringRow,numberOfRings
          call MMCommons_removeEdge( rings, i, numberOfRings)
       end do
    end if

    ! write(*,"(T20,A,I)") " Anillos hallados: ", ringRow-1
    ! write(*,"(T20,A,I)") " Anillos esperados: ", numberOfRings

    !! Imprime los anillos encontrados, borrar luego
    ! write (*,"(T10,A)") " Rings "
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! do i=1,numberOfRings
    !    numberOfColumns = size(rings(i)%values)
    !    write (*,"(T10,<numberOfColumns>I)") rings(i)%values(1,:)
    ! end do       
    
  end subroutine RingFinder_getRings

  !>
  !! @brief This routine removes a vertex and actualize the path
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in,out] this INTEGER ARRAY with the information of the edges in the pruned graph 
  !! @param [in] vertex INTEGER vertex that will be remove
  subroutine RingFinder_removeVertex(this,vertex)
    implicit none
    type(MatrixInteger), allocatable :: this(:)  
    integer, intent(in) :: vertex
    integer :: edgesInputSize
    integer, allocatable :: edgesRow(:)
    integer :: numberOfRows, numberOfColumns
    integer :: vertexEdgesSize
    integer :: i, j, k, l, m, row
    type(MatrixInteger), allocatable :: vertexEdges(:)
    type(MatrixInteger), allocatable :: edgesOutput(:)
    integer :: edgesOutputSize
    integer, allocatable :: auxEdge(:)
    integer, allocatable :: auxEdge2(:)
    ! type(MatrixInteger), allocatable :: auxEdge3(:)
    integer :: numberOfVertex
    integer :: numberOfVertex2
    integer :: numberOfVertex3
    integer :: invert
    integer :: invert2
    integer :: n
    integer :: adittionalEdges, edgesSpliceSize, edgesOutputRow


    edgesInputSize = size(this)
    
    call MMCommons_searchEdgesRow( this, edgesInputSize, vertex, edgesRow )
    
    vertexEdgesSize = size(edgesRow)

    allocate( vertexEdges( vertexEdgesSize ))
    do i=1,vertexEdgesSize
       row=edgesRow(i)
       numberOfColumns = size(this(row)%values)
       call MatrixInteger_constructor( vertexEdges(i), 1, numberOfColumns )
    end do

    do i=1,vertexEdgesSize
       row=edgesRow(i)
       numberOfColumns = size(this(row)%values)
       do j=1,numberOfColumns
          vertexEdges(i)%values(1,j) = this(row)%values(1,j)
       end do
    end do
    
    ! write(*,"(T20,A,I)") "Vertice a borrar: ", vertex
    if (vertexEdgesSize>=2) then

       edgesSpliceSize = (((vertexEdgesSize-1)*(vertexEdgesSize-1))+(vertexEdgesSize-1))/2
       adittionalEdges = edgesSpliceSize - vertexEdgesSize
       edgesOutputSize = edgesInputSize + adittionalEdges
       ! write(*,"(T20,I,A)") vertexEdgesSize, " ejes "
       allocate( edgesOutput( edgesOutputSize ))

       ! do i=1,vertexEdgesSize
       !    numberOfColumns = size(vertexEdges(i)%values)
       !    write(*,"(T20,<numberOfColumns>I)") vertexEdges(i)%values(1,:)
       ! end do
       
       n=vertexEdgesSize-1
       edgesOutputRow = 1
       do i=1,n
          do j=i+1,vertexEdgesSize
             ! write(*,"(T20,A,I,A,I)") " ejes ", i, " y ", j
             numberOfVertex = size(vertexEdges(i)%values)-1
             numberOfVertex2 = size(vertexEdges(j)%values)
             numberOfVertex3 = numberOfVertex+numberOfVertex2
             allocate( auxEdge( numberOfVertex3 )  )
             
             invert=numberOfVertex+1
             invert2=numberOfVertex2
             
             if(vertexEdges(i)%values(1,1)==vertex) then
                do k=1,numberOfVertex
                   auxEdge(k)=vertexEdges(i)%values(1,invert)
                   invert=invert-1
                end do
                if(vertexEdges(j)%values(1,1)==vertex) then
                   do l=1,numberOfVertex2
                      auxEdge(l+numberOfVertex)=vertexEdges(j)%values(1,l)
                   end do
                else
                   do l=1,numberOfVertex2
                      auxEdge(l+numberOfVertex)=vertexEdges(j)%values(1,invert2)
                      invert2=invert2-1
                   end do
                end if
             else
                do k=1,numberOfVertex
                   auxEdge(k)=vertexEdges(i)%values(1,k)
                end do
                if(vertexEdges(j)%values(1,1)==vertex) then
                   do l=1,numberOfVertex2
                      auxEdge(l+numberOfVertex)=vertexEdges(j)%values(1,l)
                   end do
                else
                   do l=1,numberOfVertex2
                      auxEdge(l+numberOfVertex)=vertexEdges(j)%values(1,invert2)
                      invert2=invert2-1
                   end do
                end if
             end if

             numberOfColumns = size(auxEdge)
             call MatrixInteger_constructor( edgesOutput(edgesOutputRow), 1, numberOfColumns )
             do m=1,numberOfColumns
                edgesOutput(edgesOutputRow)%values(1,m)=auxEdge(m)
             end do
             edgesOutputRow = edgesOutputRow + 1
             deallocate(auxEdge)
          end do
       end do
       
       ! write(*,"(T20,A)") " Ejes luego del splice "       
       ! do i=1,edgesSpliceSize
       !    numberOfColumns = size(edgesOutput(i)%values)
       !    write(*,"(T20,<numberOfColumns>I)") edgesOutput(i)%values(1,:)
       ! end do

       numberOfRows = size(edgesRow)
       do while (numberOfRows >= 1)
          k=edgesRow(1)
          call MMCommons_removeEdge( this, k, edgesInputSize)
          edgesInputSize = size(this)
          deallocate( edgesRow  )
          call MMCommons_searchEdgesRow( this, edgesInputSize, vertex, edgesRow )
          numberOfRows = size(edgesRow)
       end do

       do i=edgesSpliceSize+1,edgesOutputSize
          numberOfColumns = size(this(i-edgesSpliceSize)%values)
          call MatrixInteger_constructor( edgesOutput(i), 1, numberOfColumns )
          do j=1,numberOfColumns
             edgesOutput(i)%values(1,j)=this(i-edgesSpliceSize)%values(1,j)
          end do
       end do
       
!!! Imprime ejes despues de splice, borrar luego
       ! write (*,"(T20,A)") ""
       ! write (*,"(T20,A)") "Ejes completos luego del splice"
       ! write (*,"(T20,A)") "-----------------------------------------------------------"
       ! do i=1,edgesOutputSize
       !    numberOfColumns = size(edgesOutput(i)%values)
       !    write (*,"(T20,<numberOfColumns>I)") edgesOutput(i)%values(1,:)
       ! end do
!!!  
       deallocate( this )
       allocate( this( edgesOutputSize ) )
       do i=1,edgesOutputSize
          numberOfColumns = size(edgesOutput(i)%values)
          call MatrixInteger_constructor( this(i), 1, numberOfColumns )
       end do

       do i=1,edgesOutputSize
          numberOfColumns = size(edgesOutput(i)%values)
          do j=1,numberOfColumns
             this(i)%values(1,j) = edgesOutput(i)%values(1,j)
          end do
       end do

       deallocate( edgesOutput ) 

       
    else if (numberOfRows==1) then
       do while (numberOfRows >= 1)
          k=edgesRow(1)
          call MMCommons_removeEdge( this, k, edgesInputSize)
          edgesInputSize = size(this)
          deallocate( edgesRow  )
          call MMCommons_searchEdgesRow( this, edgesInputSize, vertex, edgesRow )
          numberOfRows = size(edgesRow)
       end do
    end if
    deallocate( edgesRow )

  end subroutine RingFinder_removeVertex

  !>
  !! @brief This function evaluates if a path is a cycle (a-b-c-d-e-f-a)
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @param [in] this INTEGER ARRAY with the information of the edges in the pruned graph 
  !! @param [in] row INTEGER path that will be evaluate
  !! @return [out] output LOGICAL if the path is a cycle returns .true.
  function RingFinder_isCycle(this,row) result(output)
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    integer, intent(in) :: row
    logical :: output
    integer :: numberOfColumns

    output=.false.
    numberOfColumns = size(this(row)%values)

    ! write(*,"(T20,A)") "Dentro de ciclos "
    ! write(*,"(T20,A,I,A)") "eje con ", numberOfColumns, " vertices"
    ! write(*,"(T20,<numberOfColumns>I)") this(row)%values(1,:)

    if ( this(row)%values(1,1) == this(row)%values(1,numberOfColumns)) then
       output = .true.
    end if

  end function RingFinder_isCycle

  !>
  !! @brief This function evaluates if a path a true path
  !! @author J.M. Rodas
  !! <b> Creation date : </b> 2014-06-02
  !! @note if the path is the type [a-b-c-d-e-c-f] is not a path
  !! @param [in] this INTEGER ARRAY with the information of the edges in the pruned graph 
  !! @param [in] row INTEGER path that will be evaluate
  !! @return [out] output LOGICAL if it is a real path returns .true.
  function RingFinder_isPath( this , row) result(output)
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    integer, intent(in) :: row
    logical :: output
    integer :: numberOfColumns
    integer :: i,j

    output=.true.
    numberOfColumns = size(this(row)%values)-1
    do i=2,numberOfColumns-1
       do j=i+1,numberOfColumns
          if (this(row)%values(1,i)==this(row)%values(1,j)) then
             output=.false.
          end if
       end do
    end do

  end function RingFinder_isPath

end module RingFinder_
