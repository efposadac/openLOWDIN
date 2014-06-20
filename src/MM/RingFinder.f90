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
module RingFinder_
  use MatrixInteger_
  use Vector_
  use MMCommons_
  use Exception_
  implicit none


  public :: &
       RingFinder_getRings

contains

  subroutine RingFinder_getRings( this, connectivityMatrix, rings )
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    type(MatrixInteger), intent(in) :: connectivityMatrix
    type(MatrixInteger), intent(out), allocatable :: rings(:)
    integer :: connectivitySize
    integer :: edgesSize
    integer :: i
    integer :: j
    type(MatrixInteger) :: vertices
    ! type(MatrixInteger), allocatable :: edgesMatrix(:)

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

    !! Imprime valores de vertices del pruned graph, borrar luego
    write (*,"(T10,A)") " "
    write (*,"(T10,A)") " Atoms in Pruned Graph"
    write (*,"(T10,A)") "--------------------------------------------"
    do i=1, connectivitySize
       write (*,"(T10,I)") vertices%values(i,:)
    end do
    !!******************************************************************************

    ! allocate( edgesMatrix( edgesSize )  )
    ! do i=1,edgesSize
    !    call MatrixInteger_constructor( edgesMatrix(i), 1, 2 )
    ! end do

    write (*,"(T10,A)") " "
    write (*,"(T10,A)") " Edges in Pruned Graph"
    write (*,"(T10,A)") "--------------------------------------------"
    write (*,"(T10,A,I)") " Tamanio =>", edgesSize
    do i=1,edgesSize
       ! do j=1,2
       !    edgesMatrix(i)%values(1,j) = this%values(i,j)
       ! end do
       write (*,"(T10,2I)") this(i)%values(1,:)
    end do


    do i=1,size(vertices%values)
       call RingFinder_removeVertex(this,i)
    end do

    ! call MMCommons_removeEdge( edgesMatrix, 3, edgesSize )

    ! edgesSize = size(edgesMatrix)

    ! write (*,"(T10,A)") " "
    ! write (*,"(T10,A)") " Edges in Pruned Graph"
    ! write (*,"(T10,A)") "--------------------------------------------"
    ! write (*,"(T10,A,I)") " Tamanio =>", edgesSize
    ! do i=1,edgesSize
    !    write (*,"(T10,2I)") edgesMatrix(i)%values(1,:)
    ! end do       

  end subroutine RingFinder_getRings

  subroutine RingFinder_removeVertex(this,vertex)
    implicit none
    type(MatrixInteger), allocatable :: this(:)  
    integer, intent(in) :: vertex
    integer :: edgesSize
    integer, allocatable :: edgesRow(:)
    integer :: numberOfRows, numberOfColumns
    integer :: i, j, row
    type(MatrixInteger), allocatable :: oldEdges(:)
    type(MatrixInteger), allocatable :: auxRings(:)
    logical :: isCycle

    edgesSize = size(this)
        
    call MMCommons_searchEdgesRow( this, edgesSize, vertex, edgesRow )
    
    write (*,"(T20,A)") "Ejes para splice"
    numberOfRows = size(edgesRow)



    allocate( oldEdges( numberOfRows ))
    do i=1,numberOfRows
       row=edgesRow(i)
       numberOfColumns = size(this(row)%values)
       call MatrixInteger_constructor( oldEdges(i), 1, numberOfColumns )
    end do

    do i=1, numberOfRows
       row=edgesRow(i)
       numberOfColumns = size(this(row)%values)
       do j=1,numberOfColumns
          oldEdges(i)%values(1,j) = this(row)%values(1,j)
       end do
          write (*,"(T20,2I)") oldEdges(i)%values(1,:)
    end do

    do i=1,numberOfRows
       isCycle=RingFinder_isCycle(oldEdges,i)
       if(isCycle) then
          write (*,"(T20,A)") "Es un anillo"
       end if
    end do

  end subroutine RingFinder_removeVertex

  function RingFinder_isCycle(this,row) result(output)
    implicit none
    type(MatrixInteger), allocatable :: this(:)
    integer, intent(in) :: row
    logical :: output
    integer :: numberOfColumns

    output=.false.
    numberOfColumns = size(this(row)%values)
    if ( this(row)%values(1,1) == this(row)%values(1,numberOfColumns)) then
       output = .true.
    end if

  end function RingFinder_isCycle

end module RingFinder_
