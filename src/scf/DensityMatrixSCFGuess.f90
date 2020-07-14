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
  function DensityMatrixSCFGuess_getGuess( speciesID ) result( output )
    implicit none
    integer, intent(in) :: speciesID
    type(Matrix) :: output
    
    type(Matrix) :: orbitals
    character(30) :: nameOfSpecies
    integer(8) :: orderOfMatrix, occupationNumber
    logical :: existPlain, existBinnary, readSuccess
    character(50) :: guessType
    character(50) :: wfnFile
    character(50) :: arguments(20)
    integer :: wfnUnit
    integer :: i,j,k
    
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
    nameOfSpecies = MolecularSystem_instance%species(speciesID)%name
    occupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    readSuccess=.false.
    call Matrix_constructor(output, orderOfMatrix, orderOfMatrix, 0.0_8  )
    call Matrix_constructor(orbitals, int(orderOfMatrix,8), int(orderOfMatrix,8), 0.0_8 )

    arguments(2) = nameOfSpecies
    arguments(1) = "COEFFICIENTS"
    
    !!Verifica el archivo que contiene los coeficientes para una especie dada
    if ( CONTROL_instance%READ_FCHK ) then

       call MolecularSystem_readFchk(trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".fchk", orbitals, output, nameOfSpecies )
       call Matrix_destructor(orbitals)
       return

    else if ( CONTROL_instance%READ_COEFFICIENTS ) then

       wfnFile=trim(CONTROL_instance%INPUT_FILE)//"plainvec"
       inquire(FILE = wfnFile, EXIST = existPlain )

       if ( existPlain ) then

          open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
          orbitals = Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
               columns= int(orderOfMatrix,4), binary=.false., arguments=arguments(1:2),failContinue=.true.)
          close(wfnUnit)
          readSuccess=.true.
       else
          wfnFile=trim(CONTROL_instance%INPUT_FILE)//"vec"
          inquire(FILE = wfnFile, EXIST = existBinnary )

          if ( existBinnary ) then

             open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
             orbitals = Matrix_getFromFile(unit=wfnUnit, rows= int(orderOfMatrix,4), &
                  columns= int(orderOfMatrix,4), binary=.true., arguments=arguments(1:2),failContinue=.true.)
             close(wfnUnit)             
             readSuccess=.true.
          end if
       end if
    end if

    !check if the orbitals were read correctly
    if(.not. allocated(orbitals%values) ) readSuccess=.false.

    if(readSuccess) print *, "Combination coefficients for ", trim(nameOfSpecies), " were read from ", trim(wfnFile)

    
    if(.not. readSuccess) then
       call Matrix_constructor(orbitals, orderOfMatrix, orderOfMatrix, 0.0_8  )

       if ( MolecularSystem_instance%species(speciesID)%isElectron ) then
          guessType=CONTROL_instance%SCF_ELECTRONIC_TYPE_GUESS
       else
          guessType=CONTROL_instance%SCF_NONELECTRONIC_TYPE_GUESS
       end if

       write(*, '(A13, A6, A28, A10)') &
            "INFO: Usign ", trim(guessType), " density guess for species: ", trim(nameOfSpecies)

       select case( trim( String_getUppercase( guessType ) ) )

       case( "ONES" )

          do i=1, occupationNumber
             output%values(i,i) =1.0_8
          end do
          return

       case( "HCORE" )
          call DensityMatrixSCFGuess_hcore( orbitals, speciesID )

          !case( "HUCKEL" )
          ! call DensityMatrixSCFGuess_huckel( output, speciesID )

       case default
          call DensityMatrixSCFGuess_exception( ERROR, "the selected guess method for "//nameOfSpecies//" is not implemented", "at program SCF module DensityMatrixSCFGuess")

       end select

    end if

    do i = 1 , orderOfMatrix
       do j = 1 , orderOfMatrix
          do k = 1 , occupationNumber
             output%values(i,j) = output%values( i,j ) + orbitals%values(i,k) * orbitals%values(j,k)
          end do
       end do
    end do
    output%values=output%values*MolecularSystem_getEta( speciesID )
      
    if ( CONTROL_instance%BUILD_MIXED_DENSITY_MATRIX ) then
       output%values(occupationNumber,:) = 0.1*output%values(occupationNumber,:)*output%values(occupationNumber+1,:)
    end if


  call Matrix_destructor(orbitals)

    
  end function DensityMatrixSCFGuess_getGuess
  
  !>
  !! @brief Diagonaliza el hamiltoniano monoelectronico para obtener los
  !! 			orbitales iniciales de partida
  subroutine DensityMatrixSCFGuess_hcore( eigenVectors, speciesID )
    implicit none
    type(Matrix), intent(inout) :: eigenVectors
    integer, intent(in) :: speciesID
    
    type(Matrix) :: hcore
    type(Matrix) :: hcoreTransformed
    type(Matrix) :: transformation
    type(Vector) :: eigenValues

    integer(8) :: orderOfMatrix
    
    orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
 
    if ( .not.allocated(eigenVectors%values) ) then
       call Matrix_constructor(eigenVectors, orderOfMatrix, orderOfMatrix )
    end if
    
    call Matrix_constructor(hcore, orderOfMatrix, orderOfMatrix )
    call Matrix_constructor(hcoreTransformed, orderOfMatrix, orderOfMatrix )
    call Matrix_constructor(transformation, orderOfMatrix, orderOfMatrix )
    call Vector_constructor(eigenValues, int( orderOfMatrix ) )
    
    hcore=WaveFunction_instance(speciesID)%hCoreMatrix
    transformation = WaveFunction_instance(speciesID)%transformationMatrix
        
    hcoreTransformed%values = matmul( matmul( transpose(transformation%values ) , hcore%values ) , transformation%values )

    call Matrix_eigen( hcoreTransformed, eigenValues, eigenVectors, SYMMETRIC )

    eigenVectors%values = matmul( transformation%values , eigenVectors%values )
    
    call Matrix_destructor( hcore )
    call Matrix_destructor( hcoreTransformed )
    call Matrix_destructor( transformation )
    call Vector_destructor( eigenValues )
    
  end subroutine DensityMatrixSCFGuess_hcore
      
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

  
  !>
  !! @brief Realiza un calculo Huckel extendido utilizando la base
  !! 		MINI de Huzinaga, proyectando esta sobre la base que se
  !! 		este utilizando
  ! subroutine DensityMatrixSCFGuess_huckel( densityMatrix, speciesID )
  !   implicit none
  !   type(Matrix), intent(inout) :: densityMatrix
  !   integer, intent(in) :: speciesID
    
    ! integer :: ocupationNumber
    ! integer(8) :: orderOfMatrix
    ! integer :: i
    
    ! ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
    ! orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )

    ! if ( .not.allocated(densityMatrix%values) ) then
       
    !    call Matrix_constructor(  densityMatrix , orderOfMatrix, orderOfMatrix, 0.0_8 )
       
    ! else
       
    !    densityMatrix%values=0.0_8
       
    ! end if
    
    ! do i=1, ocupationNumber
       
    !    densityMatrix%values(i,i) =1.0_8
       
    ! end do
    
  ! end subroutine DensityMatrixSCFGuess_huckel

  !>
  !! @brief Retorna una matriz con unos en la diagonal, hasta el numero
  !! 			de ocupacion
  ! subroutine DensityMatrixSCFGuess_ones( densityMatrix, speciesID )
  !   implicit none
    
  !   type(Matrix), intent(inout) :: densityMatrix
  !   integer, intent(in) :: speciesID
    
  !   integer :: ocupationNumber
  !   integer(8) :: orderOfMatrix
  !   integer :: i
    
  !   ocupationNumber = MolecularSystem_getOcupationNumber( speciesID )
  !   orderOfMatrix = MolecularSystem_getTotalnumberOfContractions( speciesID )
    
  !   if ( .not.allocated(densityMatrix%values) ) then
       
  !      call Matrix_constructor(  densityMatrix , orderOfMatrix, orderOfMatrix, 0.0_8 )
       
  !   else
       
  !      densityMatrix%values=0.0_8
       
  !   end if
    
  !   do i=1, ocupationNumber
       
  !      densityMatrix%values(i,i) =1.0_8
       
  !   end do
  
  ! end subroutine DensityMatrixSCFGuess_ones

  !>
  !! @brief Genera la matriz de densidad incial a partir de los coeficientes de una especie dada
  ! subroutine DensityMatrixSCFGuess_read( densityMatrix, speciesID, binnary )
  !   implicit none

  !   type(Matrix), intent(inout) :: densityMatrix
  !   integer, intent(in) :: speciesID
  !   logical, intent(in) :: binnary

  !   type(Matrix) :: vectors
  !   character(30) :: nameOfSpecie
  !   integer :: orderOfMatrix
  !   integer(8) :: numberOfMatrixElements
  !   integer :: ocupationNumber
  !   integer :: i, j, k

  !   wfnUnit = 30


  !   numberOfMatrixElements = int(orderOfMatrix, 8) ** 2_8

  !   ocupationNumber = MolecularSystem_instance%species(speciesID)%ocupationNumber

  !   !    vectors = Matrix_getFromFile(orderOfMatrix, orderOfMatrix, &
  !   !         file=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecie)//".vec", binary = .false.)


  !   call matrix_constructor( densitymatrix, int(orderofmatrix,8), int(orderofmatrix,8), 0.0_8 )

  !   densityMatrix%values =  MolecularSystem_getEta( speciesID ) * densityMatrix%values

    
  !   call Matrix_destructor(vectors)

  ! end subroutine DensityMatrixSCFGuess_read
