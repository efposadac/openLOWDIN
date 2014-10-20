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

module IndexMap_
  use Exception_
	implicit none

	!>
	! @brief  Clase estatica que contiene las funciones para el manejo de índices de contracciones
	!!
	!! @author Sergio A. Gonzalez Monico
	!!
	!! <b> Fecha de creacion : </b> 2008-08-30
	!!
	!! <b> Historial de modificaciones: </b>
	!!
	!!   - <tt> 2008-08-30 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
	!!        -# Creacion de modulo y metodos
	!!   - <tt> 2010-10-28 </tt>: Edwin Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Implementación de métodos
	!!   - <tt> 2011-02-13 </tt>: Fernando Posada ( efposadac@unal.edu.co )
	!!        -# Reescribe y adapta el módulo para su inclusion en Lowdin
	!!
	!<

public  &
		IndexMap_vectorToMatrix, &
		IndexMap_vectorToTensorR4, &
		IndexMap_matrixToTensorR4, &
		IndexMap_matrixToVector, &
		IndexMap_TensorR4ToVector, &
		IndexMap_TensorR4ToMatrix, &
		IndexMap_TransformIndexPair, &
		IndexMap_tensorR4ToVector2, &
		IndexMap_tensorR4ToVector3

contains

	!<
	!! @brief Transforma un indice i en una dupla i'.j'
	!! @author Edwin Posada, 2010
	!>
	function IndexMap_vectorToMatrix( out, cont ) result ( output )
		implicit none
		integer(8) :: cont, out, long
		integer    :: output(2)
		real	   :: tmp1
		integer    :: tmp2, tmp3

		!! calculamos ii
		long = (cont * (cont + 1)) / 2
		tmp1 = (((long- out) * 8) + 1)
		tmp2 = ((SQRT(tmp1)) - 1) / 2
		output(1) = cont - tmp2

		!!calculamos jj
		tmp3 = (tmp2 * (tmp2 + 1)) / 2
		output(2) = cont - (long - tmp3 - out)

	end function IndexMap_vectorToMatrix

	!<
	!! @brief Transforma un indice i en una cuadrupla i'.j',k',l'
	!! @author Edwin Posada, 2010
	!>
	function IndexMap_vectorToTensorR4( out , cont) result ( output )
		implicit none
		integer(8) :: cont, out, long
		integer	   :: output(4)
		integer(8) :: tmp1(2), tmp2(2), tmp3(2)

		long = (int(cont,8) * (int(cont,8) + 1)) / 2

		tmp1 = IndexMap_vectorToMatrix( out, long )

		!!calculo i, j
		tmp2 = IndexMap_vectorToMatrix( tmp1(1), cont )

		!!Calculo k, l
		tmp3 = IndexMap_vectorToMatrix( tmp1(2), cont )

		output(1:2) = tmp2(1:2)
		output(3:4) = tmp3(1:2)

	end function IndexMap_vectorToTensorR4

	!<
	!! @brief Transforma una dupla de indices i,j en una cuadrupladupla i'.j',k',l'
	!!
	!! @todo Falta por implementar
	!>
	function IndexMap_matrixToTensorR4( i, j ) result ( output )
		implicit none
		integer(8) :: i
		integer(8) :: j
		integer :: output(4)

		output = 0
		
	end function IndexMap_matrixToTensorR4

	!<
	!! @brief Transforma una dupla de indices i,j en un unico indice i'
	!!
	!! @todo Falta por implementar
	!>
	function IndexMap_matrixToVector( i, j ) result ( output )
		implicit none
		integer :: i
		integer :: j
		integer :: output

		output = 0

	end function IndexMap_matrixToVector

	!<
	!! @brief Transforma cuatro indices i,j,k,l para un tensor de rango cuatro en un unico indice
	!!		asociado a un vector(procedimiento para intra -  especies)
	!>
	function IndexMap_tensorR4ToVector( i, j, k, l, basisSizeA, basisSizeB ) result ( output )
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer, intent(in) :: k
		integer, intent(in) :: l
		integer, optional :: basisSizeA
		integer, optional :: basisSizeB

		integer(8) :: output

		integer(8) :: aux_i
		integer(8) :: aux_j
		integer(8) :: aux_k
		integer(8) :: aux_l
		integer(8) :: auxIndex

		!!************************************************************
		!! Orderna los indices de entrada segun i<j,k<l y i<=k
		!!
		if ( i > j ) then

			aux_i = j
			aux_j = i
		else
			aux_i = i
			aux_j = j

		end if

		if ( k > l ) then
			aux_k = l
			aux_l = k

		else

			aux_k = k
			aux_l = l

		end if

		!!************************************************************

		if ( .not.present( basisSizeB ) ) then

			if ( ( aux_i > aux_k ) .or.  (  ( aux_j > aux_l) .and. ( aux_i == aux_k ) ) ) then
				auxIndex = aux_i
				aux_i = aux_k
				aux_k = auxIndex
				auxIndex = aux_j
				aux_j = aux_l
				aux_l = auxIndex
			end if

			auxIndex = ( basisSizeA * ( basisSizeA + 1_8) ) / 2_8

			!! Calcula el entero asociado a los indices i,j,k,l
			output = IndexMap_transformIndexPair( 	IndexMap_transformIndexPair(aux_i, aux_j, int(basisSizeA, 8) ), &
										 	IndexMap_transformIndexPair(aux_k, aux_l, int(basisSizeA, 8) ), auxIndex )

		else

			auxIndex = ( basisSizeB * ( basisSizeB + 1_8) ) / 2_8
			output = auxIndex * ( IndexMap_transformIndexPair( aux_i, aux_j, int(basisSizeA, 8) ) - 1_8) &
				   + IndexMap_transformIndexPair( aux_k, aux_l, int(basisSizeB, 8) )

		end if

	end function IndexMap_TensorR4ToVector

	!<
	!! @brief Transforma cuatro indices i,j,k,l para un tensor de rango cuatro en un unico indice
	!!		asociado a un vector(procedimiento para inter -  especies)
	!! @author Edwin Posada, 2010
	!>
	function IndexMap_tensorR4ToVector2( i, j, k, l, basisSizeA, basisSizeB, order) result ( output )
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer, intent(in) :: k
		integer, intent(in) :: l
		integer, intent(in) :: order
		integer :: basisSizeA !! numero total de contracciones particula A
		integer :: basisSizeB !! numero total de contracciones particula B
		integer(8) :: output

		integer(8) :: aux_i, ii
		integer(8) :: aux_j, jj
		integer(8) :: aux_k, kk
		integer(8) :: aux_l, ll
		integer(8) :: auxA
		integer(8) :: auxB

		output = 0

		!!************************************************************
		!! Orderna los indices de entrada segun order (ver LibintInterface)
		!!
		select case (order)
			case(0)
				!!(SS|SS)
				ii = i
				jj = j
				kk = k
				ll = l
			case(1)
				!!(AB|CD)
				ii = i
				jj = j
				kk = k
				ll = l
			case(2)
				!!(BA|CD)
				ii = j
				jj = i
				kk = k
				ll = l
			case(3)
				!!(AB|DC)
				ii = i
				jj = j
				kk = l
				ll = k
			case(4)
				!!(BA|DC)
				ii = j
				jj = i
				kk = l
				ll = k
			case(5)
				!!(CD|AB)
				ii = k
				jj = l
				kk = i
				ll = j
			case(6)
				!!(DC|AB)
				ii = l
				jj = k
				kk = i
				ll = j
			case(7)
				!!(CD|BA)
				ii = l
				jj = k
				kk = i
				ll = j
			case(8)
				!!(DC|BA)
				ii = l
				jj = k
				kk = j
				ll = i

		end select

		if ( ii > jj ) then

			aux_i = jj
			aux_j = ii
		else
			aux_i = ii
			aux_j = jj

		end if

		if ( kk > ll ) then
			aux_k = ll
			aux_l = kk

		else

			aux_k = kk
			aux_l = ll

		end if

		auxA = aux_j - aux_i + ( ( ( 2_8 * basisSizeA * ( aux_i - 1_8 )) - ( aux_i ** 2_8) + (3_8 * aux_i) ) / 2_8 )
        auxB = aux_l - aux_k + ( ( ( 2_8 * basisSizeB * ( aux_k - 1_8 )) - ( aux_k ** 2_8) + (3_8 * aux_k) ) / 2_8 )

        output = ((basisSizeB*(basisSizeB + 1))/2) * auxA - (((basisSizeB*(basisSizeB + 1))/2) - auxB)

	end function IndexMap_TensorR4ToVector2

	!!LIBINTINTERFFACE 2 USA ESTA
	function IndexMap_tensorR4ToVector22( i, j, k, l, basisSizeA, basisSizeB) result ( output )
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer, intent(in) :: k
		integer, intent(in) :: l
		integer :: basisSizeA !! numero total de contracciones particula A
		integer :: basisSizeB !! numero total de contracciones particula B
		integer(8) :: output

		integer(8) :: aux_i, ii
		integer(8) :: aux_j, jj
		integer(8) :: aux_k, kk
		integer(8) :: aux_l, ll
		integer(8) :: auxA
		integer(8) :: auxB

		output = 0

        ii = i
        jj = j
        kk = k
        ll = l

		if ( ii > jj ) then

			aux_i = jj
			aux_j = ii
		else
			aux_i = ii
			aux_j = jj

		end if

		if ( kk > ll ) then
			aux_k = ll
			aux_l = kk

		else

			aux_k = kk
			aux_l = ll

		end if

		auxA = aux_j - aux_i + ( ( ( 2_8 * basisSizeA * ( aux_i - 1_8 )) - ( aux_i ** 2_8) + (3_8 * aux_i) ) / 2_8 )
        auxB = aux_l - aux_k + ( ( ( 2_8 * basisSizeB * ( aux_k - 1_8 )) - ( aux_k ** 2_8) + (3_8 * aux_k) ) / 2_8 )

        output = ((basisSizeB*(basisSizeB + 1))/2) * auxA - (((basisSizeB*(basisSizeB + 1))/2) - auxB)

	end function IndexMap_TensorR4ToVector22


	!<
	!! @brief Transforma cuatro indices i,j,k,l para un tensor de rango cuatro en un unico indice
	!!		asociado a un vector(procedimiento para coupling energy)
	!! @author Edwin Posada
	!>
	function IndexMap_tensorR4ToVector3( i, j, k, l, basisSizeA, basisSizeB, order ) result ( output )
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer, intent(in) :: k
		integer, intent(in) :: l
		integer, intent(in) :: basisSizeA
		integer, intent(in) :: basisSizeB
		integer, intent(in) :: order

		integer(8) :: output

		integer(8) :: aux_i, ii
		integer(8) :: aux_j, jj
		integer(8) :: aux_k, kk
		integer(8) :: aux_l, ll
		integer(8) :: auxIndex

		!!************************************************************
		!! Orderna los indices de entrada segun i<j,k<l y i<=k
		!!
		select case (order)
			case(0)
				!!(SS|SS)
				ii = i
				jj = j
				kk = k
				ll = l
			case(1)
				!!(AB|CD)
				ii = i
				jj = j
				kk = k
				ll = l
			case(2)
				!!(BA|CD)
				ii = j
				jj = i
				kk = k
				ll = l
			case(3)
				!!(AB|DC)
				ii = i
				jj = j
				kk = l
				ll = k
			case(4)
				!!(BA|DC)
				ii = j
				jj = i
				kk = l
				ll = k
			case(5)
				!!(CD|AB)
				ii = k
				jj = l
				kk = i
				ll = j
			case(6)
				!!(DC|AB)
				ii = l
				jj = k
				kk = i
				ll = j
			case(7)
				!!(CD|BA)
				ii = l
				jj = k
				kk = i
				ll = j
			case(8)
				!!(DC|BA)
				ii = l
				jj = k
				kk = j
				ll = i

		end select

		if ( ii > jj ) then

			aux_i = jj
			aux_j = ii
		else
			aux_i = ii
			aux_j = jj

		end if

		if ( kk > ll ) then
			aux_k = ll
			aux_l = kk

		else

			aux_k = kk
			aux_l = ll

		end if

		!!************************************************************

		auxIndex = ( basisSizeB * ( basisSizeB + 1_8) ) / 2_8
		output = auxIndex * ( IndexMap_transformIndexPair( aux_i, aux_j, int(basisSizeA, 8) ) - 1_8) &
			   + IndexMap_transformIndexPair( aux_k, aux_l, int(basisSizeB, 8) )

	end function IndexMap_TensorR4ToVector3


	function IndexMap_tensorR4ToVector33( i, j, k, l, basisSizeA, basisSizeB ) result ( output )
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: j
		integer, intent(in) :: k
		integer, intent(in) :: l
		integer, intent(in) :: basisSizeA
		integer, intent(in) :: basisSizeB

		integer(8) :: output

		integer(8) :: aux_i, ii
		integer(8) :: aux_j, jj
		integer(8) :: aux_k, kk
		integer(8) :: aux_l, ll
		integer(8) :: auxIndex

        ii = i
        jj = j
        kk = k
        ll = l

		if ( ii > jj ) then

			aux_i = jj
			aux_j = ii
		else
			aux_i = ii
			aux_j = jj

		end if

		if ( kk > ll ) then
			aux_k = ll
			aux_l = kk

		else

			aux_k = kk
			aux_l = ll

		end if

		!!************************************************************

		auxIndex = ( basisSizeB * ( basisSizeB + 1_8) ) / 2_8
		output = auxIndex * ( IndexMap_transformIndexPair( aux_i, aux_j, int(basisSizeA, 8) ) - 1_8) &
			   + IndexMap_transformIndexPair( aux_k, aux_l, int(basisSizeB, 8) )

	end function IndexMap_TensorR4ToVector33

	!<
	!! @brief Transforma una cuadrupla de indices i,j,k,l en una dupla i'.j'
	!!
	!! @todo Falta por implementar
	!>
	function IndexMap_TensorR4ToMatrix( i, j, k, l ) result ( output )
		implicit none
		integer :: i
		integer :: j
		integer :: k
		integer :: l
		integer(8) :: output(2)

		output = 0

	end function IndexMap_TensorR4ToMatrix

	!<
	!! @brief Transforma un par de indices i,j en un unico indice
	!>
	function IndexMap_transformIndexPair( ii,jj,maximunValueOfIndex ) result( output )
		implicit none
		integer(8) , intent(in) :: ii
		integer(8) , intent(in) :: jj
		integer(8) , intent(in) :: maximunValueOfIndex
		integer(8) :: output
		
		integer(8) :: i
		integer(8) :: j

		if ( ii > jj ) then

			i = jj
			j = ii
		else
			i = ii
			j = jj

		end if


		output = j - i + ( ( ( 2_8 * maximunValueOfIndex * ( i - 1_8 )) - ( i ** 2_8) + (3_8 * i) ) / 2_8 )

	end function IndexMap_transformIndexPair

	!>
	!! @brief  Maneja excepciones de la clase
	!<
	subroutine IndexMap_exception( typeMessage, description, debugDescription)
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
	
	end subroutine IndexMap_exception

end module IndexMap_
