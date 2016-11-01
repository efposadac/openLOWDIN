!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Clase encargada de realizar transformacion de integrales atomicas a  moleculares
!!
!!  Esta clase reliza la transformacion de integrales de orbitales atomicos a orbitales moleculares,
!!  creando una interface al algoritmo de   Yamamoto, Shigeyoshi; Nagashima, Umpei.
!!  Computer Physics Communications, 2005, 166, 58-65
!!
!! @author Sergio Gonzalez
!!
!! <b> Fecha de creacion : </b> 2009-07-07
!!   - <tt> 2009-07-07 </tt>: Sergio Gonzalez ( sagonzalez@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2011-02-15 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin
!!   - <tt> 2013-10-03 </tt>: Jose Mauricio Rodas (jmrodasr@unal.edu.co)     
!!        -# Adapts to Lowdin 2               
!!   - <tt> 2014-08-26 </tt>: Jorge Charry (jacharrym@unal.edu.co)     
!!        -# Adapts this module to works indepently from MP2 program
!<
module ReadTransformedIntegrals_
  use Matrix_
  use Vector_
  use IndexMap_
  use Exception_
  use MolecularSystem_
  implicit none

  public :: &
       ReadTransformedIntegrals_readOneSpecies, &
       ReadTransformedIntegrals_readTwoSpecies

  private

contains

  !! This subroutine was moved to lowdinCore  

  subroutine ReadTransformedIntegrals_readOneSpecies( specieID, matrixContainer )
    implicit none
    type(Vector) :: matrixContainer

    integer(8) :: numberOfIntegrals
    integer(8) :: auxIndex
    real(8),dimension(791) :: integralValue
    integer :: iter
    integer :: errorValue
    integer :: bufferA
    integer :: bufferB
    integer :: bufferSize
    integer,dimension(4,791) :: indexBuffer
    integer :: specieID
    integer :: numberOfContractions
    character(10) :: nameOfSpecie
    character(255) :: prefixOfFile
    integer :: unidOfOutputForIntegrals

    integer :: p, q, r, s
    real(8) :: auxIntegralValue

    numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
    nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )

    prefixOfFile =""//trim(nameOfSpecie)

    !    if ( nameOfSpecie == "E-BETA" ) then
    !        prefixOfFile =""//trim("E-ALPHA")
    !    end if

    unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

    select case( CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD )

    case ( "A" ) 

       !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
       open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
            status='old',access='sequential', form='unformatted' )

       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

       numberOfIntegrals = int( ( (  numberOfContractions * ( numberOfContractions + 1.0_8 ) / 4.0_8 ) * &
            ( (  numberOfContractions * (  numberOfContractions + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )

       call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )
       matrixContainer%values = 0.0_8

       do
          read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) bufferA,bufferB,integralValue,indexBuffer

          bufferSize = iabs( bufferA)
          !*!!  improve this! remove if
          if ( bufferA /= 0 ) then
             do iter=1,bufferSize

                auxIndex = IndexMap_tensorR4ToVector(indexBuffer(1,iter),indexBuffer(2,iter), &
                     indexBuffer(3,iter), indexBuffer(4,iter), numberOfContractions )

                matrixContainer%values( auxIndex ) = integralValue(iter)

             end do
          end if

          if ( bufferA <= 0 ) exit

       end do

       close(unidOfOutputForIntegrals)

    case ( "B" ) 

       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

       numberOfIntegrals = int( ( (  numberOfContractions * ( numberOfContractions + 1.0_8 ) / 4.0_8 ) * &
            ( (  numberOfContractions * (  numberOfContractions + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )

       call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )
       matrixContainer%values = 0.0_8


       !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
       open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
            status='old',access='sequential', form='unformatted' )


       do
          read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

          if ( p <= 0 ) exit

          auxIndex = IndexMap_tensorR4ToVector( p, q, r, s, numberOfContractions )
          matrixContainer%values( auxIndex ) = auxIntegralValue

       end do

       close(unidOfOutputForIntegrals)

    case ( "C" ) 

       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

       numberOfIntegrals = int( ( (  numberOfContractions * ( numberOfContractions + 1.0_8 ) / 4.0_8 ) * &
            ( (  numberOfContractions * (  numberOfContractions + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )

       call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )
       matrixContainer%values = 0.0_8

       !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
       open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
            status='old',access='sequential', form='unformatted' )


       do
          read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

          if ( p <= 0 ) exit

          auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(numberOfContractions,8 ))
          matrixContainer%values( auxIndex ) = auxIntegralValue

       end do

       close(unidOfOutputForIntegrals)

    case ( "D" ) 

       if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

       numberOfIntegrals = int( ( (  numberOfContractions * ( numberOfContractions + 1.0_8 ) / 4.0_8 ) * &
            ( (  numberOfContractions * (  numberOfContractions + 1.0_8) / 2.0_8 ) + 1.0_8) ), 8 )

       call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )
       matrixContainer%values = 0.0_8

       !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
       open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
            status='old',access='sequential', form='unformatted' )


       do
          read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

          if ( p <= 0 ) exit

          auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(numberOfContractions,8 ))
          matrixContainer%values( auxIndex ) = auxIntegralValue

       end do

       close(unidOfOutputForIntegrals)

    end select

  end subroutine ReadTransformedIntegrals_readOneSpecies

  subroutine ReadTransformedIntegrals_readTwoSpecies( specieID, otherSpecieID, matrixContainer )
    implicit none
    type(Vector) :: matrixContainer

    integer(8) :: numberOfIntegrals
    integer(8) :: auxIndex
    real(8),dimension(791) :: integralValue
    integer :: iter
    integer :: errorValue
    integer :: bufferA
    integer :: bufferB
    integer :: bufferSize
    integer :: lowerIndices(2), upperIndeces(2), counter(2)
    integer,dimension(4,791) :: indexBuffer

    integer :: specieID
    integer :: otherSpecieID
    integer :: numberOfContractions, bias
    character(10) :: nameOfSpecie
    character(10) :: nameOfOtherSpecie
    character(255) :: prefixOfFile
    integer :: unidOfOutputForIntegrals

    integer :: p, q, r, s
    real(8) :: auxIntegralValue

    select case( CONTROL_instance%INTEGRALS_TRANSFORMATION_METHOD )

    case ( "A" ) 

       if ( otherSpecieID > SpecieID ) then

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )

          prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)

          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE


          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) bufferA,bufferB,integralValue,indexBuffer
             bufferSize = iabs( bufferA)
             if ( bufferA /= 0 ) then

                do iter=1,bufferSize

                   counter=1
                   if ( indexBuffer(1,iter ) > bias )  then
                      indexBuffer(1,iter)=indexBuffer(1,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(1,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(1,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(2,iter ) > bias ) then
                      indexBuffer(2,iter)=indexBuffer(2,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(2,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(2,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(3,iter ) > bias ) then
                      indexBuffer(3,iter)=indexBuffer(3,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(3,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(3,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(4,iter ) > bias ) then
                      indexBuffer(4,iter)=indexBuffer(4,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(4,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(4,iter)
                      counter(2)=2
                   end if

                   auxIndex = IndexMap_tensorR4ToVector(lowerIndices(1),lowerIndices(2),upperIndeces(1),upperIndeces(2), &
                        bias, numberOfContractions - bias )

                   matrixContainer%values(auxIndex)=integralValue(iter)

                end do
             end if

             if ( bufferA <= 0 ) exit

          end do

          close(unidOfOutputForIntegrals)


       else ! ( otherSpecieID < specieID ) 

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(otherSpecieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )

          prefixOfFile =""//trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)

          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) bufferA,bufferB,integralValue,indexBuffer
             bufferSize = iabs( bufferA)
             if ( bufferA /= 0 ) then

                do iter=1,bufferSize

                   counter=1
                   if ( indexBuffer(1,iter ) > bias )  then
                      indexBuffer(1,iter)=indexBuffer(1,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(1,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(1,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(2,iter ) > bias ) then
                      indexBuffer(2,iter)=indexBuffer(2,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(2,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(2,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(3,iter ) > bias ) then
                      indexBuffer(3,iter)=indexBuffer(3,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(3,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(3,iter)
                      counter(2)=2
                   end if

                   if ( indexBuffer(4,iter ) > bias ) then
                      indexBuffer(4,iter)=indexBuffer(4,iter)-bias
                      upperIndeces( counter(1) ) = indexBuffer(4,iter)
                      counter(1)=2
                   else
                      lowerIndices( counter(2) ) = indexBuffer(4,iter)
                      counter(2)=2
                   end if

                   auxIndex = IndexMap_tensorR4ToVector(upperIndeces(1),upperIndeces(2),lowerIndices(1),lowerIndices(2), &
                        numberOfContractions - bias, bias )

                   matrixContainer%values(auxIndex)=integralValue(iter)

                end do
             end if

             if ( bufferA <= 0 ) exit

          end do

          close(unidOfOutputForIntegrals)



       end if

    case ( "B" ) 

       if ( otherSpecieID > SpecieID ) then


          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )

          prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)

          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

             if ( p <= 0 ) exit

             auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(bias,8), &
                  int(numberOfContractions - bias,8)  )
             matrixContainer%values( auxIndex ) = auxIntegralValue

          end do

          close(unidOfOutputForIntegrals)

       end if


    case ( "C" ) 

       if ( otherSpecieID > SpecieID ) then
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )
          prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)


          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

             if ( p <= 0 ) exit

             auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(bias,8),  &
                  int(numberOfContractions - bias,8)  )
             matrixContainer%values( auxIndex ) = auxIntegralValue

          end do

          close(unidOfOutputForIntegrals)

       else 
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )
          prefixOfFile =""//trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)

          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) r, s, p, q, auxIntegralValue

             if ( p <= 0 ) exit

             auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(bias,8),  &
                  int(numberOfContractions - bias,8)  )
             matrixContainer%values( auxIndex ) = auxIntegralValue

          end do

          close(unidOfOutputForIntegrals)



       end if

    case ( "D" ) 

       if ( otherSpecieID > SpecieID ) then
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )
          prefixOfFile =""//trim(nameOfSpecie)//"."//trim(nameOfOtherSpecie)


          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) p, q, r, s, auxIntegralValue

             if ( p <= 0 ) exit

             auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(bias,8),  &
                  int(numberOfContractions - bias,8)  )
             matrixContainer%values( auxIndex ) = auxIntegralValue

          end do

          close(unidOfOutputForIntegrals)

       else 
          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID) &
               + MolecularSystem_getTotalNumberOfContractions(otherSpecieID)
          bias = MolecularSystem_getTotalNumberOfContractions(specieID)

          nameOfSpecie= trim(  MolecularSystem_getNameOfSpecie( specieID ) )
          nameOfOtherSpecie= trim(  MolecularSystem_getNameOfSpecie( otherSpecieID ) )
          prefixOfFile =""//trim(nameOfOtherSpecie)//"."//trim(nameOfSpecie)

          unidOfOutputForIntegrals = CONTROL_instance%UNIT_FOR_MP2_INTEGRALS_FILE

          !! Accesa el archivo binario con las integrales en terminos de orbitales moleculares
          open(unit=unidOfOutputForIntegrals, file=trim(prefixOfFile)//"moint.dat", &
               status='old',access='sequential', form='unformatted' )

          if ( allocated(matrixContainer%values ) ) deallocate(matrixContainer%values)

          numberOfIntegrals = ( bias    *  ( ( bias + 1.0_8) / 2.0_8 ) ) * &
               ( (numberOfContractions-bias) * ( ( (numberOfContractions-bias) + 1.0_8 ) / 2.0_8 ) )

          call Vector_constructor( matrixContainer, numberOfIntegrals, 0.0_8 )

          matrixContainer%values = 0.0_8

          do
             read(UNIT=unidOfOutputForIntegrals,IOSTAT=errorValue) r, s, p, q, auxIntegralValue

             if ( p <= 0 ) exit

             auxIndex = IndexMap_tensorR4ToVectorB( int(p,8), int(q,8), int(r,8), int(s,8), int(bias,8),  &
                  int(numberOfContractions - bias,8)  )
             matrixContainer%values( auxIndex ) = auxIntegralValue

          end do

          close(unidOfOutputForIntegrals)

       end if


    end select

  end subroutine ReadTransformedIntegrals_readTwoSpecies

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine ReadTransformedIntegrals_exception( typeMessage, description, debugDescription)
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

  end subroutine ReadTransformedIntegrals_exception

end module ReadTransformedIntegrals_
