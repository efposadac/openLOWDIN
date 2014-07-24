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
!! @brief Modulo para administracion de particulas.
!!  Este modulo crea y manipula grupos de particulas cuanticas o puntuales
!! @author Sergio A. Gonzalez Monico
!! <b> Fecha de creacion : </b> 2008-08-14
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-07-20 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulos y metodos basicos
!!   - <tt> 2011-02-14 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el modulo para su inclusion en Lowdin, inplementa open-shell case.
!!   - <tt> 2013-04-30 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Eliminates the XML dependence
module ParticleManager_
  use Exception_
  use Units_
  use Matrix_
  use Vector_
  use Particle_
  implicit none
    
  type ParticleManager
     
     type(particle), pointer :: particlePtr
     
  end type ParticleManager

  
  !< enum ParticleManager__showFlags {
  integer, parameter :: LABELS_NUMERATED = 1
  !< }
  
  public :: &
        ParticleManager_show, &
!        ParticleManager_getValuesOfFreeCoordinates, &
!                ParticleManager_getPositionOfCenterOfOptimizacion, &
!        ParticleManager_getNumberOfFreeCoordinates, &
!        ParticleManager_getNumberOfCoordinates, &
        ParticleManager_getNumberOfCentersOfOptimization, &
        ParticleManager_getCartesianMatrixOfCentersOfOptimization, &
!        ParticleManager_getDistanceMatrix, &
!        ParticleManager_getCenterOfOptimization, &
!        ParticleManager_getSpecieID, &
!        ParticleManager_getNickName, &
!        ParticleManager_getSymbol, &
!        ParticleManager_getParticlePtr, &
!        ParticleManager_getContractionPtr, &
!        ParticleManager_getFactorOfInterchangeIntegrals, &
!        ParticleManager_getOwnerOfPuntualParticle, &
!        ParticleManager_getOwnerCenter, &
!        ParticleManager_getNameOfPuntualParticle, &
!        ParticleManager_getSymbolOfPuntualParticle, &
!        ParticleManager_getOcupationNumber,&
!        ParticleManager_getMultiplicity,&
!        ParticleManager_getEta, &
!        ParticleManager_getLambda,&
!        ParticleManager_getKappa, &
!        ParticleManager_getParticlesFraction,&
!        ParticleManager_getMass, &
!        ParticleManager_getMassInPosition, &
!        ParticleManager_getTotalMass, &
!        ParticleManager_getCenterOfMass, &
!        ParticleManager_getCharge, &
!        ParticleManager_getLabelsOfContractions, &
        ParticleManager_getLabelsOfCentersOfOptimization, &
        ParticleManager_getChargesOfCentersOfOptimization, &
!        ParticleManager_isQuantum, &
!        ParticleManager_isCenterOfOptimization, &
!        ParticleManager_isComponentFixed, &
!        ParticleManager_iterateSpecie, &
!        ParticleManager_beginSpecie, &
!        ParticleManager_endSpecie, &
!        ParticleManager_rewindSpecies, &
        ParticleManager_setOrigin, &
        ParticleManager_setParticlesPositions, &
        ParticleManager_setOwner
!        ParticleManager_puntualParticlesEnergy,&
!        ParticleManager_changeOriginOfSystem, &
!        ParticleManager_searchSpecie

   type(ParticleManager), pointer :: ParticleManager_instance(:)
   
contains
  
  !>
  !! @brief Muestra los atributos de todas las particulas en el Administrador de particulas
  !! @todo Falta adicionar procedimineto para que muestre una sola particula
  subroutine ParticleManager_show( nameOfParticle )
    implicit none

    character(*), optional ::nameOfParticle
    
    integer :: i, j
    integer :: cosa
    
    if ( present( nameOfParticle ) ) then
       
       !! Aqui adicionar codigo para una sola particula
       
    else
       
       write(*,*) "================================"
       write(*,*) " LOADED PARTICLES IN THE SYSTEM "
       write(*,*) "================================"
       
       do i=1, size(ParticleManager_instance)

          call Particle_show( ParticleManager_instance(i)%particlePtr )
             
       end do
       
    end if
       
  end subroutine ParticleManager_show
  
  
  !>
  !! @brief Ajusta el owner de la particula especificada
  subroutine ParticleManager_setOwner()
    implicit none

    integer :: i
    integer counter
    integer :: particleID
    
    do particleID = 1, size(ParticleManager_instance)
       do i = 1, particleID - 1

          if( abs( ParticleManager_instance(i)%particlePtr%origin(1) - ParticleManager_instance(particleID)%particlePtr%origin(1) ) &
               < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
               abs( ParticleManager_instance(i)%particlePtr%origin(2) - ParticleManager_instance(particleID)%particlePtr%origin(2) ) &
               < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
               abs( ParticleManager_instance(i)%particlePtr%origin(3) - ParticleManager_instance(particleID)%particlePtr%origin(3) ) &
               < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
             
             !!
             !! Selecciona quien debe ser  la particula propietaria, cuando encuentra dos
             !! especies sobre el mismo origen. Se opta por la mas pesada como la propietaria
             !! y la liviana como hija
             !!
             if ( ParticleManager_instance(particleID)%particlePtr%mass  > ParticleManager_instance(i)%particlePtr%mass )	then
                
                
                call Particle_setOwner(ParticleManager_instance(i)%particlePtr, owner = particleID )
                ParticleManager_instance(particleID)%particlePtr%isCenterOfOptimization =.true.
                ParticleManager_instance(i)%particlePtr%isCenterOfOptimization =.false.
                call Particle_setChild(ParticleManager_instance(particleID)%particlePtr, i)
                call Particle_removeChilds( ParticleManager_instance(i)%particlePtr )
                
             else
                
                call Particle_setOwner( ParticleManager_instance(particleID)%particlePtr, &
                     owner= ParticleManager_instance(i)%particlePtr%owner )
                call Particle_removeChilds( ParticleManager_instance(particleID)%particlePtr )
                ParticleManager_instance(particleID)%particlePtr%isCenterOfOptimization =.false.
                ParticleManager_instance(i)%particlePtr%isCenterOfOptimization =.true.
                
             end if
             
          end if
          
       end do
    end do

    
  end subroutine ParticleManager_setOwner

  !>
  !! @brief Ajusta el origen de la particula especificada, o el origen de una
  !! unica  contraccion o gaussina especificadas.
  !! @param this Particula cuantica o carga puntual.
  !! @param origin Origen de la particula o contraccion especificada.
  !! @param contractionNumber Numero de contraccion dentro de una lista.
  !! @param gaussianNumber Numero de gausiana dentro de una lista.
  recursive subroutine ParticleManager_setOrigin( this, origin )
    implicit none
    type(Particle)  :: this
    real(8)  :: origin(3)
    
    integer :: i
    
    this%origin = origin
    
    if (this%isQuantum .eqv. .true.) then

       do i = 1, this%basis%length
          this%basis%contraction(i)%origin = origin
       end do
       
    end if
    
    if ( this%isCenterOfOptimization .and. allocated(this%childs) ) then
       
       do i=1, size(this%childs)
          call ParticleManager_setOrigin( ParticleManager_instance(this%childs(i))%particlePtr, this%origin )
       end do
       
    end if
    
  end subroutine ParticleManager_setOrigin

  !>
  !! @brief calcula la energia total para una especie especificada
  !! @todo Este algoritmo debe ser optimizado
  subroutine ParticleManager_setParticlesPositions( cartesianVector )
    implicit none
    type(Vector), intent(in),optional :: cartesianVector

    integer :: i
    integer :: j
    integer :: k
    integer :: m
    real(8) :: origin(3)
    
    k = 0
    
    do i=1, size(ParticleManager_instance)

       if( k < size(cartesianVector%values) ) then

          if ( ParticleManager_instance(i)%particlePtr%isCenterOfOptimization ) then
             
             origin = ParticleManager_instance(i)%particlePtr%origin
             
             do j=1,3
                k=k+1
                if( .not.ParticleManager_instance(i)%particlePtr%fixComponent(j) )  origin(j)=cartesianVector%values( k )
             end do

             call ParticleManager_setOrigin( ParticleManager_instance(i)%particlePtr, origin )

             if ( allocated(ParticleManager_instance(i)%particlePtr%childs) ) then

                do j=1, size( ParticleManager_instance(i)%particlePtr%childs )
                   m=ParticleManager_instance(i)%particlePtr%childs(j)
                   call ParticleManager_setOrigin( ParticleManager_instance(m)%particlePtr, origin )
                end do

             end if
             
          end if

       else
          
          return
          
       end if

    end do
    
  end subroutine ParticleManager_setParticlesPositions


  !>
  !! @brief gets the number of centers of optimization.
  function ParticleManager_getNumberOfCentersOfOptimization( ) result( output)
    implicit none
    
    integer :: output
    
    integer :: i

    output = 0
    
    do i=1, size(ParticleManager_instance)
       if ( ParticleManager_instance(i)%particlePtr%isCenterOfOptimization ) output =output + 1
    end do
       
  end function ParticleManager_getNumberOfCentersOfOptimization

  !>
  !! @brief Retorna la masa total del sistema en las unidades solicitadas
  function ParticleManager_getTotalMass( unid ) result( output )
    implicit none
    
    character(*), optional :: unid
    real(8) :: output

    integer :: i

    output = 0.0_8
    
    do i=1, size( ParticleManager_instance)
       output = output + ParticleManager_instance(i)%particlePtr%mass * ParticleManager_instance(i)%particlePtr%internalSize
    end do

    if ( present(unid) ) then
       
       select case( trim(unid) )
          
       case ("AU")
          
       case("SI")
          
          output = output * kg
          
       case("AMU")
          
          output = output * AMU

       case default
          
       end select
       
    end if
    
  end function ParticleManager_getTotalMass

  !>
  !! @brief Retorna el valor de masa localizada en un punto dado del espacio.
  function ParticleManager_getMassInPosition( position, unid) result( output )
    implicit none
    real(8),intent(in) :: position(3)
    
    character(*), optional :: unid
    real(8) :: output
    
    integer :: i

    output = 0.0_8
    
    do i=1,size( ParticleManager_instance )

       if( abs( position(1) - ParticleManager_instance(i)%particlePtr%origin(1) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
            abs( position(2) - ParticleManager_instance(i)%particlePtr%origin(2) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
            abs( position(3) - ParticleManager_instance(i)%particlePtr%origin(3) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then

          output = output + ParticleManager_instance(i)%particlePtr%mass * ParticleManager_instance(i)%particlePtr%internalSize

       end if

    end do

    if ( present(unid) ) then

       select case( trim(unid) )

       case ("AU")

       case("SI")

          output = output * kg

       case("AMU")

          output = output * AMU

       case default

       end select

    end if


  end function ParticleManager_getMassInPosition

  
  !> 
  !! @brief Returns a matrix of the position of centers for optimization.
  function ParticleManager_getCartesianMatrixOfCentersOfOptimization() result( output)
    
    implicit none
    type(Matrix) :: output
    integer :: ssize
    integer :: i
    integer :: j
    
    ssize = ParticleManager_getNumberOfCentersOfOptimization()
    
    call Matrix_constructor( output, int(ssize,8), 3_8 )
    
    j =0

    do i = 1, size(ParticleManager_instance)
       
       if ( ParticleManager_instance(i)%particlePtr%isCenterOfOptimization ) then
          
          j = j + 1          
          output%values(j,:) = ParticleManager_instance(i)%particlePtr%origin

       end if
       
    end do

  end function ParticleManager_getCartesianMatrixOfCentersOfOptimization

  !> 
  !! @brief Returns the distance matrix
  function ParticleManager_getDistanceMatrix() result( output )
    implicit none
    type(Matrix) :: output
    
    type(Matrix) :: auxMatrix
    integer :: i
    integer :: j
    integer :: ssize
    
    auxMatrix = ParticleManager_getCartesianMatrixOfCentersOfOptimization()

    ssize = size(auxMatrix%values,dim=1)

    call Matrix_constructor(output, int(ssize,8), int(ssize,8) )
    
    do i=1,ssize
       do j=1, ssize
          output%values(i,j) = sqrt( sum( ( auxMatrix%values(i,:) - auxMatrix%values(j,:) )**2) )
       end do
    end do
    
    call Matrix_destructor( auxMatrix )
    
  end function ParticleManager_getDistanceMatrix
  
  !>
  !! @brief Retorna las etiquetas asocoadas a los centros de optimizacion
  function ParticleManager_getLabelsOfCentersOfOptimization(flags) result( output )
    implicit none
    character(10),allocatable :: output(:)
    integer, optional :: flags
    
    character(20) :: number
    integer :: numberOfCenters
    integer :: internalFlags
    integer :: i
    integer :: j
    
    internalFlags=0
    if( present(flags) ) internalFlags=flags

    numberOfCenters = ParticleManager_getNumberOfCentersOfOptimization()
    
    if ( allocated( output ) ) deallocate( output )
    allocate( output( numberOfCenters ) )

    j = 0
    
    do i=1, size(ParticleManager_instance )
       
       if ( ParticleManager_instance(i)%particlePtr%isCenterOfOptimization ) then
          
          j = j + 1
          
          select case(internalFlags)
             
          case( LABELS_NUMERATED )
             
             write(number,*) j
             number = adjustl(trim(number))
             output(j) = trim( ParticleManager_instance(i)%particlePtr%nickname )//"("//trim(number)//")"
             
          case default
             
             output(j) = trim( ParticleManager_instance(i)%particlePtr%nickname )
             
          end select
          
       end if
       
    end do
    
  end function ParticleManager_getLabelsOfCentersOfOptimization

!   !<
!   !! @brief Indica si la particula es un centro de optimizacion o no
!   !!
!   !>
!   function ParticleManager_isCenterOfOptimization( iterator ) result( output )
!     implicit none
!     integer, intent(in) :: iterator
!     logical :: output

!     output = ParticleManager_instance%particlesPtr(iterator)%isCenterOfOptimization

!   end function ParticleManager_isCenterOfOptimization



!   !<
!   !! @brief Indica si la componete cartesiana debe modificarse durante la optimzacion de geometria
!   !!
!   !>
!   function ParticleManager_isComponentFixed( iterator, component ) result( output )
!     implicit none
!     integer, intent(in) :: iterator
!     integer, intent(in) :: component
!     logical :: output

!     output = ParticleManager_instance%particlesPtr(iterator)%fixComponent(component)

!   end function ParticleManager_isComponentFixed


  
!   !<
!   !! @brief Retorna el nombre de especie cuantica
!   !!
!   !>
!   function ParticleManager_iterateSpecie() result( output )
!     implicit none
!     character(30) :: output

!     ParticleManager_instance%iteratorOfSpecie =ParticleManager_instance%iteratorOfSpecie+ 1

!     if ( ParticleManager_instance%iteratorOfSpecie > ParticleManager_instance%numberOfQuantumSpecies ) &
!          ParticleManager_instance%iteratorOfSpecie = 1

!     if ( ParticleManager_instance%isInstanced ) then

!        output = trim( Map_getKey( ParticleManager_instance%speciesID,iterator=ParticleManager_instance%iteratorOfSpecie ) )
!     else

!        call ParticleManager_exception( ERROR, "Error in interateSpecies algorithm", &
!             "Class object ParticleManager in the iterateSpecie function")

!     end if

!   end function ParticleManager_iterateSpecie

!   !<
!   !! @brief Retorna el nombre de especie cuantica
!   !!
!   !>
!   subroutine ParticleManager_rewindSpecies()
!     implicit none


!     if ( ParticleManager_instance%isInstanced ) then

!        ParticleManager_instance%iteratorOfSpecie = 0

!     end if

!   end subroutine ParticleManager_rewindSpecies

!   !<
!   !! @brief Retorna un ID asociado a la especie especificada
!   !!
!   !>
!   function ParticleManager_getNameOfSpecie( specieID ) result( output )
!     implicit none
!     integer, intent(in) :: specieID
!     character(30) :: output

!     if ( ParticleManager_instance%isInstanced ) then

!        output = Map_getKey( ParticleManager_instance%speciesID, real(specieID,8) )

!     else
!        call ParticleManager_exception( ERROR, "You should instance the ParticleManager before use this function", &
!             "Class object ParticleManager in the getNameOfSpecie function")
!     end if

!   end function ParticleManager_getNameOfSpecie

!   !<
!   !!   @brief Retorna un iterador a laprimera especie en el sistema molecular
!   !!
!   !>
!   function ParticleManager_beginSpecie() result( output )
!     implicit none
!     integer :: output
!     if ( ParticleManager_instance%isInstanced ) then

!        output = Map_begin(ParticleManager_instance%speciesID )

!     else
!        call ParticleManager_exception( ERROR, "You should instance the ParticleManager before use this function", &
!             "Class object ParticleManager in the beginSpecie function")
!     end if

!   end function ParticleManager_beginSpecie

!   !**
!   !   @brief Retorna un iterador a la ultima especie en el sistema molecular
!   !
!   !**
!   function ParticleManager_endSpecie() result( output )
!     implicit none
!     integer :: output
!     if (ParticleManager_instance%isInstanced ) then

!        output = Map_end( ParticleManager_instance%speciesID )

!     else
!        call ParticleManager_exception( ERROR, "You should instance the ParticleManager before use this function", &
!             "Class object ParticleManager in the endSpecie function")
!     end if

!   end function ParticleManager_endSpecie

!   !<
!   !! Retorna un puntero la contraccion iesima del tipo de particula dado
!   !!
!   !! @todo Adicionar condicionales para verificar limites de especie o contraccion, antes de devolverla
!   !>
!   function ParticleManager_getContractionPtr( specieID, numberOfContraction ) result( output )
!     implicit none
!     integer, intent(in) :: specieID
!     integer, intent(in) :: numberOfContraction

!     type(ContractedGaussian), pointer :: output


!     integer :: particleID
!     integer :: contractionID

!     particleID = ParticleManager_instance%idsOfContractionsForSpecie(specieID)%contractionID(numberOfContraction)%particleID
!     contractionID=ParticleManager_instance%idsOfContractionsForSpecie(specieID)%contractionID(numberOfContraction)%contractionIDInParticle

!     output => ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID)

!   end function ParticleManager_getContractionPtr

!   !<
!   !! @brief Retorna el origen de la particula puntual especificada
!   !!
!   !>
!   function ParticleManager_getOriginOfPuntualParticle( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     real(8) :: output(3)

!     output = ParticleManager_instance%particlesPtr( ParticleManager_instance%idsOfPuntualParticles%particleID( specieID ) )%origin

!   end function ParticleManager_getOriginOfPuntualParticle

!   !<
!   !! @brief Retorna el origen de la particula especificada
!   !!
!   !>
!   function ParticleManager_getOrigin( iterator) result( output )
!     implicit none
!     integer, intent(in) :: iterator
!     real(8) :: output(3)


!     output = ParticleManager_instance%particlesPtr( iterator )%origin

!   end function ParticleManager_getOrigin

!   !<
!   !! @brief Retorna el nombre de la particula especificada
!   !!
!   !>
!   function ParticleManager_getName( iterator) result( output )
!     implicit none
!     integer, intent(in) :: iterator

!     character(30) :: output

!     output = trim( ParticleManager_instance%particlesPtr( iterator )%name )


!   end function ParticleManager_getName

!   !<
!   !! @brief Retorna el nickname de la particula especificada
!   !!
!   !>
!   function ParticleManager_getNickName( iterator) result( output )
!     implicit none
!     integer, intent(in) :: iterator

!     character(30) :: output

!     output = trim(ParticleManager_instance%particlesPtr(iterator)%nickname)


!   end function ParticleManager_getNickName

!   !<
!   !! @brief Retorna el simbolo de la particula especificada
!   !!
!   !>
!   function ParticleManager_getSymbol( iterator) result( output )
!     implicit none
!     integer, intent(in) :: iterator

!     character(30) :: output

!     output = trim( ParticleManager_instance%particlesPtr( iterator )%symbol )

!   end function ParticleManager_getSymbol

!   !<
!   !! @brief Retorna un puntero a la particula especificada .
!   !!
!   !>
!   function ParticleManager_getParticlePtr( iindex) result( output )
!     implicit none
!     integer, intent(in) :: iindex

!     type(Particle), pointer :: output

!     output => ParticleManager_instance%particlesPtr( iindex )


!   end function ParticleManager_getParticlePtr



!   !<
!   !! @brief Retorna un vector con los valores de las coordenadas libres del sistema molecular
!   !!
!   !! @todo Este algoritmo dede ser optimizado
!   !>
!   function ParticleManager_getValuesOfFreeCoordinates() result( output )
!     implicit none
!     type(Vector) :: output


!     real(8), allocatable :: auxVector(:)
!     integer, allocatable :: auxCentersOfOptimization(:,:)
!     integer :: numberOfParticles
!     integer :: i
!     integer :: j
!     integer :: k

!     numberOfParticles = size(ParticleManager_instance%particlesPtr)
!     allocate(auxVector(numberOfParticles*3) )
!     allocate(auxCentersOfOptimization(numberOfParticles*3,2) )

!     k=0
!     !! Determina el numero de coordenadas libres durante la optimizacion
!     do i=1, numberOfParticles

!        if (  ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) then
!           do j=1, 3
!              if ( .not.ParticleManager_instance%particlesPtr(i)%fixComponent(j) ) then
!                 k = k +1
!                 auxVector( k) = ParticleManager_instance%particlesPtr(i)%origin(j)
!                 auxCentersOfOptimization(k,1) = ParticleManager_instance%particlesPtr(i)%owner
!                 auxCentersOfOptimization(k,2) = j
!              end if
!           end do
!        end if

!     end do

!     call Vector_constructor(output, k)
!     output%values = auxVector(1:k)
!     allocate( ParticleManager_instance%centersOfOptimization(k,2) )
!     ParticleManager_instance%centersOfOptimization = auxCentersOfOptimization(1:k,:)

!     deallocate(auxCentersOfOptimization)
!     deallocate(auxVector)

!   end function ParticleManager_getValuesOfFreeCoordinates

!   !<
!   !! @brief Retorna un vector con los valores de las coordenadas libres del sistema molecular
!   !!
!   !! @todo Este algoritmo dede ser optimizado
!   !>
!   function ParticleManager_getPositionOfCenterOfOptimizacion() result( output )
!     implicit none
!     type(Vector) :: output


!     real(8), allocatable :: auxVector(:)
!     integer, allocatable :: auxCentersOfOptimization(:,:)
!     integer :: numberOfParticles
!     integer :: i
!     integer :: j
!     integer :: k

!     numberOfParticles = size(ParticleManager_instance%particlesPtr)
!     allocate(auxVector(numberOfParticles*3) )
!     allocate(auxCentersOfOptimization(numberOfParticles*3,2) )

!     k=0
!     !! Determina el numero de coordenadas libres durante la optimizacion
!     do i=1, numberOfParticles

!        if (  ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) then
!           do j=1, 3
!              k = k +1
!              auxVector( k) = ParticleManager_instance%particlesPtr(i)%origin(j)
!              auxCentersOfOptimization(k,1) = ParticleManager_instance%particlesPtr(i)%owner
!              auxCentersOfOptimization(k,2) = j
!           end do
!        end if

!     end do

!     call Vector_constructor(output, k)
!     output%values = auxVector(1:k)
!     if( .not.allocated( ParticleManager_instance%centersOfOptimization)) allocate( ParticleManager_instance%centersOfOptimization(k,2) )
!     ParticleManager_instance%centersOfOptimization = auxCentersOfOptimization(1:k,:)

!     deallocate(auxCentersOfOptimization)
!     deallocate(auxVector)

!   end function ParticleManager_getPositionOfCenterOfOptimizacion



!   function ParticleManager_getNumberOfCoordinates() result( output )
!     implicit none
!     integer :: output


!     integer :: i
!     integer :: j

!     output = 0

!     do i=1, size(ParticleManager_instance%particlesPtr )

!        if ( ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) then

!           output =output + 1

!        end if

!     end do

!     output = output*3

!   end function ParticleManager_getNumberOfCoordinates


!   function ParticleManager_getNumberOfFreeCoordinates() result( output )
!     implicit none
!     integer :: output
!     integer :: i

!     output = 0

!     do i=1, size(ParticleManager_instance%particlesPtr )

!        if ( ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) then

!           if (.not.ParticleManager_instance%particlesPtr(i)%fixComponent(1) ) output =output + 1
!           if (.not.ParticleManager_instance%particlesPtr(i)%fixComponent(2) ) output =output + 1
!           if (.not.ParticleManager_instance%particlesPtr(i)%fixComponent(3) ) output =output + 1

!        end if

!     end do

!   end function ParticleManager_getNumberOfFreeCoordinates


!   function ParticleManager_getNumberOfCentersOfOptimization( fragmentNumber ) result( output)

!     implicit none
!     integer :: output
!     integer, optional, intent(in) :: fragmentNumber
!     integer :: i

!     output = 0
!     !
!     if( present(fragmentNumber)) then
!        do i=1, size(ParticleManager_instance%particlesPtr)
!           if ( ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization .and. &
!                ParticleManager_instance%particlesPtr(i)%fragmentNumber == fragmentNumber) output =output + 1
!        end do
!     else
!        do i=1, size(ParticleManager_instance%particlesPtr)
!           if ( ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) output =output + 1
!        end do

!     end if


!   end function ParticleManager_getNumberOfCentersOfOptimization




!   function ParticleManager_getCenterOfOptimization( coordinate ) result( output )
!     implicit none
!     integer, intent(in) :: coordinate
!     integer :: output(2)

!     output(1) = ParticleManager_instance%centersOfOptimization(coordinate,1)
!     output(2) = ParticleManager_instance%centersOfOptimization(coordinate,2)

!   end  function ParticleManager_getCenterOfOptimization


!   !<
!   !! @brief Retorna la carga de la particula puntual especificada
!   !!
!   !>
!   function ParticleManager_getChargeOfPuntualParticle( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     real(8) :: output

!     output = ParticleManager_instance%particlesPtr( ParticleManager_instance%idsOfPuntualParticles%particleID( specieID ) )%charge

!   end function ParticleManager_getChargeOfPuntualParticle

!   !<
!   !! @brief Retorna el indice del porpietario de la particula puntual especificada
!   !!
!   !>
!   function ParticleManager_getOwnerOfPuntualParticle( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     integer :: output

!     output = ParticleManager_instance%particlesPtr( ParticleManager_instance%idsOfPuntualParticles%particleID( specieID ) )%owner

!   end function ParticleManager_getOwnerOfPuntualParticle


!   !<
!   !! @brief Retorna el numero de fragmento en el sistema molecular
!   !!
!   !>
!   function ParticleManager_getNumberOfFragments() result( output )
!     implicit none
!     integer :: output

!     output = ParticleManager_instance%numberOfFragments

!   end function ParticleManager_getNumberOfFragments


!   !<
!   !! @brief Retorna el indice del porpietario de la particula puntual especificada
!   !!
!   !>
!   function ParticleManager_getOwnerCenter( specieID) result( output )
!     implicit none
!     integer, intent(in), optional :: specieID

!     integer :: output

!     output = ParticleManager_instance%particlesPtr( specieID )%owner

!   end function ParticleManager_getOwnerCenter

!   !<
!   !! @brief Retorna el nombre de la particula puntual especificada
!   !!
!   !>
!   function ParticleManager_getNameOfPuntualParticle( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     character(30) :: output

!     output = ParticleManager_instance%particlesPtr( ParticleManager_instance%idsOfPuntualParticles%particleID(specieID) )%name

!   end function ParticleManager_getNameOfPuntualParticle

!   function ParticleManager_getSymbolOfPuntualParticle( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     character(30) :: output

!     output = ParticleManager_instance%particlesPtr( ParticleManager_instance%idsOfPuntualParticles%particleID(specieID) )%symbol

!   end function ParticleManager_getSymbolOfPuntualParticle


!   !<
!   !! @brief 	Retorna el factor multiplicativo para las integrals de intercambio en la construccion
!   !!		de la matrix de particula independiente (G)
!   !>
!   function ParticleManager_getFactorOfInterchangeIntegrals( specieID) result( output )
!     implicit none
!     integer :: specieID

!     real(8) :: output

!     output = Map_getValue( ParticleManager_instance%kappa, iterator=specieID) &
!          / Map_getValue( ParticleManager_instance%eta, iterator=specieID)

!   end function ParticleManager_getFactorOfInterchangeIntegrals


!   !<
!   !! @brief 	Retorna el numero de ocupacion para la especie solicitada
!   !>
!   function ParticleManager_getOcupationNumber( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     integer :: output

!     output = Map_getValue( ParticleManager_instance%ocupationNumber, iterator=specieID)

!   end function ParticleManager_getOcupationNumber


!   !<
!   !! @brief 	Retorna la multiplicidad para la especie solicitada
!   !>
!   function ParticleManager_getMultiplicity( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID
!     integer :: output


!     output = Map_getValue( ParticleManager_instance%multiplicity, iterator=specieID )

!   end function ParticleManager_getMultiplicity


!   !<
!   !! @brief 	Retorna la constante de acoplamiento "eta" para la particula solicitada
!   !>
!   function ParticleManager_getEta( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     real(8) :: output

!     output = real( Map_getValue( ParticleManager_instance%eta, iterator=specieID), 8 )

!   end function ParticleManager_getEta

!   !<
!   !! @brief 	Retorna la constante de acoplamiento "lambda" para la particula solicitada
!   !>
!   function ParticleManager_getLambda( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     real(8) :: output

!     output = real( Map_getValue( ParticleManager_instance%lambda, iterator=specieID), 8 )

!   end function ParticleManager_getLambda

!   !<
!   !! @brief 	Retorna la constante de acoplamiento "kappa" para la particula solicitada
!   !>
!   function ParticleManager_getKappa( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     real(8) :: output

!     output = real( Map_getValue( ParticleManager_instance%kappa, iterator=specieID), 8 )

!   end function ParticleManager_getKappa

!   function ParticleManager_getParticlesFraction( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID
!     real(8) :: output


!     output = real( Map_getValue( ParticleManager_instance%particlesFraction, iterator=specieID), 8 )

!   end function ParticleManager_getParticlesFraction


!   !<
!   !! @brief 	Retorna la masa de la especie solicitada
!   !>
!   function ParticleManager_getMass( specieID) result( output )
!     implicit none
!     integer :: specieID
!     real(8) :: output


!     output = real( Map_getValue( ParticleManager_instance%mass, iterator=specieID), 8 )

!   end function ParticleManager_getMass


!   !<
!   !! @brief Retorna la masa de la especie solicitada
!   !>
!   function ParticleManager_getTotalMass( unid) result( output )
!     implicit none
!     real(8) :: output
!     character(*), optional :: unid


!     integer :: i

!     output = 0.0_8
!     do i=1, size( ParticleManager_instance%particlesPtr)
!        output = output + ParticleManager_instance%particlesPtr(i)%mass*ParticleManager_instance%particlesPtr(i)%internalSize
!     end do

!     if ( present(unid) ) then

!        select case( trim(unid) )

!        case ("AU")

!        case("SI")

!           output = output * kg

!        case("AMU")

!           output = output * AMU

!        case default

!        end select

!     end if

!   end function ParticleManager_getTotalMass


!   !<
!   !! @brief 	Retorna el centro de masa para las particulas en el administrados
!   !!
!   !! @warning Se asume que la masa de toda parrticula esta localizada donde se centre
!   !!		    su funcion base o la particula fija.
!   !>
!   function ParticleManager_getCenterOfMass( ) result( output )
!     implicit none

!     real(8) :: output(3)

!     integer :: i

!     output = 0.0_8

!   end function ParticleManager_getCenterOfMass


!   !<
!   !! @brief 	Retorna la carga de la especie solicitada
!   !>
!   function ParticleManager_getCharge( specieID, iterator) result( output )
!     implicit none
!     integer, optional :: specieID

!     integer, optional :: iterator
!     real(8) :: output

!     if (present(iterator) ) then
!        output = ParticleManager_instance%particlesPtr(iterator)%charge
!     else
!        output = real( Map_getValue( ParticleManager_instance%charge, iterator=specieID), 8 )
!     end if

!   end function ParticleManager_getCharge

!   !<
!   !! @brief 	Retorna las etiquetas de las contracciones gaussianas asociadas al
!   !!		momento angular de la especie es especificada
!   !>
!   function ParticleManager_getLabelsOfContractions( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID
!     character(19),allocatable :: output(:)

!     integer :: numberOfContractions
!     integer :: totalNumberOfContractions
!     character(30) :: symbolOfSpecie
!     character(9), allocatable :: shellCode(:)
!     integer :: i
!     integer :: j
!     integer :: k
!     integer :: m
!     integer :: particleID
!     integer :: contractionID

!     m = 0
!     do i=1, ParticleManager_instance%numberOfQuantumSpecies

!        if ( specieID ==  i ) then

!           numberOfContractions = int( Map_getValue(ParticleManager_instance%numberOfContractions, iterator = i ) )
!           totalNumberOfContractions = ParticleManager_getTotalNumberOfContractions(i)

!           if (allocated(output)) deallocate(output)
!           allocate( output(totalNumberOfContractions) )
!           output=""

!           do j=1, numberOfContractions

!              particleID = ParticleManager_instance%idsOfContractionsForSpecie(i)%contractionID(j)%particleID
!              contractionID = ParticleManager_instance%idsOfContractionsForSpecie(i)%contractionID(j)%contractionIDInParticle

!              if(allocated(shellCode)) deallocate(shellCode)
!              allocate(shellCode(ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital))

!              shellCode = ContractedGaussian_getShellCode(ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID))
!              do k = 1, ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID)%numCartesianOrbital				
!                 m = m + 1

!                 write (output(m),"(I5,A2,A6,A2,A4)") m,"  ", &
!                      trim(ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID)%name), "  ", &
!                      trim(shellCode(k))//" "


!              end do
!           end do

!        end if

!     end do

!   end function ParticleManager_getLabelsOfContractions


  !>                                                                                                                  !! @brief Retorna las cargas asociadas a los centros de optimizacion. 24-02-2014 
  !>

   function ParticleManager_getChargesOfCentersOfOptimization( flags ) result( output )
    implicit none
    real(8), allocatable :: output(:)
    integer, optional :: flags

    character(20) :: number
    integer :: numberOfCenters
    integer :: internalFlags
    integer :: i
    integer :: j

    internalFlags=0
    if( present(flags) ) internalFlags=flags

    numberOfCenters = ParticleManager_getNumberOfCentersOfOptimization()

    if ( allocated( output ) ) deallocate( output )
    allocate( output( numberOfCenters ) )

    j = 0

    do i=1, size(ParticleManager_instance )

       if ( ParticleManager_instance(i)%particlePtr%isCenterOfOptimization ) then

          j = j + 1

          select case(internalFlags)

          case( LABELS_NUMERATED )

             write(number,*) j
             number = adjustl(trim(number))
             output(j) = ParticleManager_instance(i)%particlePtr%totalCharge

          case default

             output(j) = ParticleManager_instance(i)%particlePtr%totalCharge
          end select

       end if

    end do


   end function ParticleManager_getChargesOfCentersOfOptimization

!   !<
!   !! @brief 	Indica si es una particula fija
!   !>
!   function ParticleManager_isQuantum( specieID) result( output )
!     implicit none
!     integer, intent(in) :: specieID

!     logical :: output

!     output = ParticleManager_instance%particlesPtr(specieID)%isQuantum

!   end function ParticleManager_isQuantum

!   !<
!   !! @brief calcula la energia total para una especie especificada
!   !!
!   !>
!   function ParticleManager_puntualParticlesEnergy() result( output )
!     implicit none
!     real(8) :: output

!     integer :: i
!     integer :: j
!     real(8) :: deltaOrigin(3)
!     real(8) :: tmp

!     output =0.0_8

!     do i=1, size( ParticleManager_instance%particlesPtr )

!        if ( .not.ParticleManager_instance%particlesPtr(i)%isQuantum ) then

!           do j = i + 1 , size( ParticleManager_instance%particlesPtr )

!              if ( .not.ParticleManager_instance%particlesPtr(j)%isQuantum ) then

!                 deltaOrigin = 	ParticleManager_instance%particlesPtr(i)%origin &
!                      - ParticleManager_instance%particlesPtr(j)%origin

!                 output=output + ( ( ParticleManager_instance%particlesPtr(i)%charge &
!                      * ParticleManager_instance%particlesPtr(j)%charge )&
!                      / sqrt( sum( deltaOrigin**2.0_8 ) ) )

!              end if
!           end do

!        end if

!     end do

!   end function ParticleManager_puntualParticlesEnergy





!   subroutine ParticleManager_fixAllNucleous()
!     implicit none
!     integer :: i

!     i = 0
!     do i=1, size(ParticleManager_instance%particlesPtr )

!        if ( ParticleManager_instance%particlesPtr(i)%isCenterOfOptimization ) then

!           ParticleManager_instance%particlesPtr(i)%isQuantum = .false.

!        end if

!     end do

!     ParticleManager_instance%numberOfPuntualParticles=0

!     do i=1, size(ParticleManager_instance%particlesPtr )

!        if ( .not.ParticleManager_instance%particlesPtr(i)%isQuantum ) then

!           ParticleManager_instance%numberOfPuntualParticles = ParticleManager_instance%numberOfPuntualParticles +1
!        end if

!     end do

!   end subroutine ParticleManager_fixAllNucleous

!   !<
!   !! @brief metodo con propositos de depuracion
!   !!	permite evaluar el cambio en coeficientes de contraccion durante tiempo de ejecucion
!   !!	se pueden incluir otras verificaciones
!   !>
!   subroutine ParticleManager_verifyCoefficients( message, threshold )
!     implicit none
!     character(*) :: message
!     real(8) :: threshold


!     integer ::i
!     integer ::j
!     integer :: particleID
!     integer :: contractionID



!     print *,""
!     print *,"Begin Verification"
!     print *,"---------------------------------------------------------------------------------"
!     print *,trim(message)


!     do i=1,ParticleManager_instance%numberOfQuantumSpecies

!        do j=1,size( ParticleManager_instance%idsOfContractionsForSpecie(i)%contractionID )

!           particleID = ParticleManager_instance%idsOfContractionsForSpecie(i)%contractionID(j)%particleID
!           contractionID = ParticleManager_instance%idsOfContractionsForSpecie(i)%contractionID(j)%contractionIDInParticle

!           if ( sum(ParticleManager_instance%particlesPtr(particleID)%basis%contractions(contractionID)%contractionCoefficients) &
!                >threshold ) then

!              call ParticleManager_exception( ERROR, "Coefficients has been changed", &
!                   "Class object ParticleManager in the verifyCoefficients function")

!           end if

!        end do

!     end do
!     print *,""
!     print *,"---------------------------------------------------------------------------------"
!     print *,""

!   end subroutine ParticleManager_verifyCoefficients

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine ParticleManager_exception( typeMessage, description, debugDescription)
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

  end subroutine ParticleManager_exception

end module ParticleManager_
