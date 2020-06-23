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
!! @brief  Clase estatica que contiene los procedimientos para construir
!!	   matrices de densidad iniciales relevantes en procesos tipo SCF
!! @author S. A. Gonzalez
!! <b> Creation data : </b> 02-16-11
!! <b> History change: </b>
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el modulo para su inclusion en Lowdin
module DensityMatrixSCFGuess_
  use MolecularSystem_
  use Matrix_
  use Vector_
  use WaveFunction_
  use String_
  use Exception_
  implicit none
  
  
  public  &
       DensityMatrixSCFGuess_getGuess
  
  
  private
  
  
contains
  
  !>
  !! @brief Obtiene la matriz de densidad inicial
  function DensityMatrixSCFGuess_getGuess( densityType, speciesID ) result( output )
    implicit none
    character(*), intent(in) :: densityType
    integer, intent(in) :: speciesID
    type(Matrix) :: output
    
    type(Matrix) :: orbitals
    character(30) :: nameOfSpecies
    character(50) :: wfnFile
    integer(8) :: orderOfMatrix
    logical :: existPlain, existBinnary
    
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
    nameOfSpecies = MolecularSystem_instance%species(speciesID)%name
    
    call Matrix_constructor(output, orderOfMatrix, orderOfMatrix, 0.0_8  )

    !!Verifica el archivo que contiene los coeficientes para una especie dada
    if ( CONTROL_instance%READ_FCHK ) then


       call Matrix_constructor(orbitals, orderOfMatrix, orderOfMatrix, 0.0_8  )
       call MolecularSystem_readFchk( orbitals, output, nameOfSpecies )
       call Matrix_destructor(orbitals)
       
    else if ( CONTROL_instance%READ_COEFFICIENTS ) then

       wfnFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
       inquire(FILE = wfnFile, EXIST = existBinnary )

       if ( existBinnary ) then

       call DensityMatrixSCFGuess_read( output, speciesID, existBinnary )

       else

          wfnFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
          inquire(FILE = wfnFile, EXIST = existPlain )
          
          if ( existPlain ) then

             call DensityMatrixSCFGuess_read( output, speciesID, existBinnary )

          else

             call DensityMatrixSCFGuess_exception( ERROR, "I did not find any .vec coefficients file", "At HF program, at DensityMatrixSCFGuess_getGuess")
          
          end if

       
       end if
    

    else
       
       write(*, '(A13, A6, A28, A10)') &
          "INFO: Usign ", trim( String_getUppercase( densityType ) ), " density guess for species: ", trim(nameOfSpecies)
       
       select case( trim( String_getUppercase( densityType ) ) )
          
       case( "ONES" )
          call DensityMatrixSCFGuess_ones( output, speciesID )
          
       case( "HCORE" )

          call DensityMatrixSCFGuess_hcore( output, speciesID )
          
       case( "HUCKEL" )
          call DensityMatrixSCFGuess_huckel( output, speciesID )
          
       case default
          call DensityMatrixSCFGuess_hcore( output, speciesID )
          
       end select
       
    end if
    
  end function DensityMatrixSCFGuess_getGuess
  
  !>
  !! @brief Realiza un calculo Huckel extendido utilizando la base
  !! 		MINI de Huzinaga, proyectando esta sobre la base que se
  !! 		este utilizando
  subroutine DensityMatrixSCFGuess_huckel( densityMatrix, speciesID )
    implicit none
    type(Matrix), intent(inout) :: densityMatrix
    integer, intent(in) :: speciesID
    
    integer :: ocupationNumber
    integer(8) :: orderOfMatrix
    integer :: i
    
    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )

    if ( .not.allocated(densityMatrix%values) ) then
       
       call Matrix_constructor(  densityMatrix , orderOfMatrix, orderOfMatrix, 0.0_8 )
       
    else
       
       densityMatrix%values=0.0_8
       
    end if
    
    do i=1, ocupationNumber
       
       densityMatrix%values(i,i) =1.0_8
       
    end do
    
  end subroutine DensityMatrixSCFGuess_huckel
  
  !>
  !! @brief Diagonaliza el hamiltoniano monoelectronico para obtener los
  !! 			orbitales iniciales de partida
  subroutine DensityMatrixSCFGuess_hcore( hcoreDensityMatrix, speciesID )
    implicit none
    type(Matrix), intent(inout) :: hcoreDensityMatrix
    integer, intent(in) :: speciesID
    
    type(Matrix) :: auxGuess
    type(Matrix) :: auxTransformation
    type(Matrix) :: eigenVectors
    type(Vector) :: eigenValues
    character(30) :: nameOfSpecie
    integer :: ocupationNumber
    integer(8) :: orderOfMatrix
    integer :: i
    integer :: j
    integer :: k
    
    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
 
    nameOfSpecie = MolecularSystem_instance%species(speciesID)%name

    if ( .not.allocated(hcoreDensityMatrix%values) ) then
     
       call Matrix_constructor(  hcoreDensityMatrix , orderOfMatrix, orderOfMatrix, 0.0_8 )
       
    end if
    
    call Matrix_constructor(auxGuess, orderOfMatrix, orderOfMatrix )
    call Matrix_constructor(auxTransformation, orderOfMatrix, orderOfMatrix )
    call Matrix_constructor(eigenVectors, orderOfMatrix, orderOfMatrix )
    call Vector_constructor(eigenValues, int( orderOfMatrix ) )
    
    auxGuess=WaveFunction_getHcoreMatrix( trim(nameOfSpecie ) )
    auxTransformation = WaveFunction_getTransformationMatrix( trim(nameOfSpecie ) )
        
    auxGuess%values = matmul( matmul( transpose(auxTransformation%values ) , auxGuess%values ) , auxTransformation%values )
    
    call Matrix_eigen( auxGuess, eigenValues, eigenVectors, SYMMETRIC )
    
    auxGuess%values = matmul( auxTransformation%values , eigenVectors%values)
    
    hcoreDensityMatrix%values = 0.0_8
    
    do i = 1 , orderOfMatrix
       do j = 1 , orderOfMatrix
          do k = 1 , ocupationNumber
             
             hcoreDensityMatrix%values(i,j) = hcoreDensityMatrix%values( i,j ) + ( auxGuess%values(i,k) * auxGuess%values(j,k) )
             
          end do
       end do
    end do
    
    if ( CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX ) then

      hcoreDensityMatrix%values(ocupationNumber,:) = 0.1*hcoreDensityMatrix%values(ocupationNumber,:)*hcoreDensityMatrix%values(ocupationNumber+1,:)

    end if


    hcoreDensityMatrix%values =  MolecularSystem_getEta( speciesID ) * hcoreDensityMatrix%values
    
    call Matrix_destructor( auxGuess )
    call Matrix_destructor( auxTransformation )
    call Matrix_destructor( eigenVectors )
    call Vector_destructor( eigenValues )
    
  end subroutine DensityMatrixSCFGuess_hcore
  
  !>
  !! @brief Retorna una matriz con unos en la diagonal, hasta el numero
  !! 			de ocupacion
  subroutine DensityMatrixSCFGuess_ones( densityMatrix, speciesID )
    implicit none
    
    type(Matrix), intent(inout) :: densityMatrix
    integer, intent(in) :: speciesID
    
    integer :: ocupationNumber
    integer(8) :: orderOfMatrix
    integer :: i
    
    ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
    
    if ( .not.allocated(densityMatrix%values) ) then
       
       call Matrix_constructor(  densityMatrix , orderOfMatrix, orderOfMatrix, 0.0_8 )
       
    else
       
       densityMatrix%values=0.0_8
       
    end if
    
    do i=1, ocupationNumber
       
       densityMatrix%values(i,i) =1.0_8
       
    end do
  
  end subroutine DensityMatrixSCFGuess_ones

  !>
  !! @brief Genera la matriz de densidad incial a partir de los coeficientes de una especie dada
  subroutine DensityMatrixSCFGuess_read( densityMatrix, speciesID, binnary )
    implicit none

    type(Matrix), intent(inout) :: densityMatrix
    integer, intent(in) :: speciesID
    logical, intent(in) :: binnary

    type(Matrix) :: vectors
    character(30) :: nameOfSpecie
    integer :: orderOfMatrix
    integer(8) :: numberOfMatrixElements
    integer :: numberOfContractions
    integer :: ocupationNumber
    integer :: i, j, k
    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit

    wfnUnit = 30


    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
    nameOfSpecie = MolecularSystem_instance%species(speciesID)%name

    numberOfContractions = MolecularSystem_getNumberOfContractions( speciesID )
    numberOfMatrixElements = int(orderOfMatrix, 8) ** 2_8

    ocupationNumber = MolecularSystem_instance%species(speciesID)%ocupationNumber

    !    vectors = Matrix_getFromFile(orderOfMatrix, orderOfMatrix, &
    !         file=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//".vec", binary = .false.)

    arguments(2) = MolecularSystem_getNameOfSpecie(speciesID)
    arguments(1) = "COEFFICIENTS"

    call Matrix_constructor( vectors, int(orderOfMatrix,8), int(orderOfMatrix,8), 0.0_8 )

    if ( binnary ) then

       wfnFile = trim(CONTROL_instance%INPUT_FILE)//"vec"

       !! Open file for wavefunction
       open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

       vectors = Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2))

       close(wfnUnit)

    else 

       wfnFile = trim(CONTROL_instance%INPUT_FILE)//"plainvec"

       !! Open file for wavefunction
       open(unit=wfnUnit, file=trim(wfnFile), status="old", form="formatted")

       vectors = Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
            columns= int(orderOfMatrix,4), binary=.false., arguments=arguments(1:2))
       close(wfnUnit)

    end if

    call matrix_constructor( densitymatrix, int(orderofmatrix,8), int(orderofmatrix,8), 0.0_8 )

    do i = 1 , orderOfMatrix
       do j = 1 , orderOfMatrix
          do k = 1 , ocupationNumber
             densityMatrix%values(i,j) = &
                  densityMatrix%values( i,j ) &
                  + ( vectors%values(i,k) * vectors%values(j,k) )
          end do
       end do
    end do

    densityMatrix%values =  MolecularSystem_getEta( speciesID ) * densityMatrix%values

    print *, "Combination coefficients for ", trim(nameOfSpecie), " were read from ", trim(wfnFile)
    
    call Matrix_destructor(vectors)

  end subroutine DensityMatrixSCFGuess_read
    
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine DensityMatrixSCFGuess_exception( typeMessage, description, debugDescription)
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
    
  end subroutine DensityMatrixSCFGuess_exception
  
end module DensityMatrixSCFGuess_
