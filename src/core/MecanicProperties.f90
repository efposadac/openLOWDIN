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
!! @brief Clase encargada de calcular propiedades mecanicas del sistema molecular
!!		como momentos de inercia, matrices de constantes de fuerza, etc.
!!
!! @author Sergio Gonzalez
!!
!! <b> Fecha de creacion : </b> 2009-05-26
!!   - <tt> 2009-05-26 </tt>: Sergio Gonzalez ( sagonzalez@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!   - <tt> 2013-04-28 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Update to the new LOWDIN standar
module MecanicProperties_
  use Exception_
  use Species_
  use Matrix_
  use Vector_
  use AtomicElement_
  use LinearVectorialSpace_
  use ParticleManager_

  type , public :: MecanicProperties
     
     character(30) :: name
     integer :: numberOfPoints
     integer :: numberOfVibrationalVars
     integer :: numberOfRotationalAndTranslationalVars
     real(8), allocatable :: puntualMasses(:)
     real(8) :: molecularMass
     real(8) :: centerOfMass(3)
     type(Matrix) :: cartesianCoordinates
     type(Matrix) :: molecularInertiaTensor
     type(Matrix) :: forceConstansProjector
     type(Matrix) :: transformationMatrix
     
  end type MecanicProperties
  
  public :: &
       MecanicProperties_constructor, &
       MecanicProperties_destructor, &
       MecanicProperties_isInstanced
  
  private :: &
       MecanicProperties_exception

contains
  
  
  !>
  !! Define el constructor para la clase
  subroutine MecanicProperties_constructor( this )
    implicit none
    type(MecanicProperties) :: this
    
    integer :: i
    
    this%numberOfPoints = ParticleManager_getNumberOfCentersOfOptimization()    
    allocate(  this%puntualMasses( this%numberOfPoints ) )
    this%molecularMass = ParticleManager_getTotalMass( unid="AMU")
    this%cartesianCoordinates = ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    
    do i=1, this%numberOfPoints
       this%puntualMasses(i) = ParticleManager_getMassInPosition( this%cartesianCoordinates%values(i,:), unid="AMU" )
    end do
    
  end subroutine MecanicProperties_constructor
  
  !>
  !! @brief Define el destructor para clase
  !! @param thisPointer Funcion base
  subroutine MecanicProperties_destructor( this )
    implicit none
    type(MecanicProperties) :: this
    
    call Matrix_destructor( this%molecularInertiaTensor )
    call Matrix_destructor( this%cartesianCoordinates )
    if(allocated(this%puntualMasses) ) deallocate( this%puntualMasses )
    
  end subroutine MecanicProperties_destructor
  
  !>
  !! @brief  Calcula y retorna el centro de masa del sistema
  subroutine MecanicProperties_getCenterOfMass( this )
    implicit none
    type(MecanicProperties) :: this
    
    integer :: i
    
    this%cartesianCoordinates = ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    
    do i=1, this%numberOfPoints
       this%puntualMasses(i) = ParticleManager_getMassInPosition( this%cartesianCoordinates%values(i,:), unid="AMU" )
    end do
    
    !! Calculo el centro de masa molecular
    this%centerOfMass(1) = sum(this%puntualMasses(:)*this%cartesianCoordinates%values(:,1))
    this%centerOfMass(2) = sum(this%puntualMasses(:)*this%cartesianCoordinates%values(:,2))
    this%centerOfMass(3) = sum(this%puntualMasses(:)*this%cartesianCoordinates%values(:,3))
    this%centerOfMass =  this%centerOfMass / this%molecularMass
    
  end subroutine  MecanicProperties_getCenterOfMass

  !>
  !! @brief  Calcula y retorna el tensor de inercia molecular
  function MecanicProperties_getMolecularInertiaTensor( this ) result(output)
    implicit  none
    type(MecanicProperties) :: this
    type(Matrix) :: output
    
    integer :: i
    real(8) ::  geometryInCM(3) !! Vector con posicion desplazada al centro de masa
    real(8) :: ssign
    type(Matrix) :: inertiaMatrix
    type(Vector) :: eigenValues
    
    call Matrix_constructor(output, 3_8, 3_8, 0.0_8)
    call Matrix_constructor(inertiaMatrix, 3_8, 3_8, 0.0_8)
    call Vector_constructor(eigenValues,3,0.0_8)
    
    call MecanicProperties_getCenterOfMass(this)
    
    do i=1, this%numberOfPoints
       !! Traslada el vector de posicion al centro de masa molecular
       geometryInCM= this%cartesianCoordinates%values(i,:)-this%centerOfMass
       
       !!*************
       !! Calcula el aporte de la i-esima componente al sistema de ecuaciones
       !!*****
       inertiaMatrix%values(1,1) = inertiaMatrix%values(1,1) + this%puntualMasses(i) * ( geometryInCM(2)**2 + geometryInCM(3)**2)
       inertiaMatrix%values(2,2) = inertiaMatrix%values(2,2) + this%puntualMasses(i) * ( geometryInCM(1)**2 + geometryInCM(3)**2)
       inertiaMatrix%values(3,3) = inertiaMatrix%values(3,3) + this%puntualMasses(i) * ( geometryInCM(1)**2 + geometryInCM(2)**2)
       inertiaMatrix%values(1,2) = inertiaMatrix%values(1,2) - this%puntualMasses(i) * ( geometryInCM(1) * geometryInCM(2) )
       inertiaMatrix%values(2,3) = inertiaMatrix%values(2,3) - this%puntualMasses(i) * ( geometryInCM(2) * geometryInCM(3) )
       inertiaMatrix%values(1,3) = inertiaMatrix%values(1,3) - this%puntualMasses(i) * ( geometryInCM(1) * geometryInCM(3) )
    end do
    inertiaMatrix%values(2,1) =inertiaMatrix%values(1,2)
    inertiaMatrix%values(3,2) =inertiaMatrix%values(2,3)
    inertiaMatrix%values(3,1) =inertiaMatrix%values(1,3)
    
    !! Calcula la componentes del vector de inercia molecular
    call Matrix_eigen( inertiaMatrix, eigenValues, output, SYMMETRIC )
    
    !! Verifica el signo del determinante
    ssign= output%values(1,1)*(output%values(2,2)*output%values(3,3)-output%values(3,2)*output%values(2,3)) &
         -  output%values(1,2)*(output%values(2,1)*output%values(3,3)-output%values(3,1)*output%values(2,3)) &
         + output%values(1,3)*(output%values(2,1)*output%values(3,2)-output%values(3,1)*output%values(2,2) )
    
    !! Preserva la handedness del tensor de inercia
    if (  ssign < 0.0 ) then
       output%values(1,2)=-output%values(1,2)
       output%values(2,2)=-output%values(2,2)
       output%values(3,2)=-output%values(3,2)
    endif
    
    !! Verifica la validez del vector de inercia molecular calculado
    if ( abs( eigenValues%values(1) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
         abs( eigenValues%values(2) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. &
         abs( eigenValues%values(3) ) < CONTROL_instance%DOUBLE_ZERO_THRESHOLD ) then
       
       call MecanicProperties_exception( ERROR, "The inertia tensor is invalid, The diagonal elements are zero !! ", &
            "Class object MecanicProperties in getMolecularInertiaTensor() " )
       
    end if
    
    call Matrix_destructor( inertiaMatrix )
    call Vector_destructor(eigenValues)
    
  end function MecanicProperties_getMolecularInertiaTensor


  
!   function MecanicProperties_getTotallyAsymmetricTensor(this) result(output)
!     implicit none
!     type(MecanicProperties) :: this
!     real(8) :: output(3,3,3)
    
!     output=0.0
!     output(2,3,1) =  1.0
!     output(3,2,1) = -1.0
!     output(1,3,2) = -1.0
!     output(3,1,2) =  1.0
!     output(1,2,3) =  1.0
!     output(2,1,3) = -1.0
    
!   end function MecanicProperties_getTotallyAsymmetricTensor
  
!   !>
!   !! @brief  Calcula y retorna el centro de masa unitario del sistema
!   function MecanicProperties_getUnitaryCenterOfMass( this ) result( output )
!     implicit none
!     type(MecanicProperties) :: this
!     real(8) :: output(3)
    
!     this%cartesianCoordinates = ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    
!     output(1)=sum(this%cartesianCoordinates%values(:,1))
!     output(2)=sum(this%cartesianCoordinates%values(:,2))
!     output(3)=sum(this%cartesianCoordinates%values(:,3))
!     output=output/this%numberOfPoints
    
!   end function MecanicProperties_getUnitaryCenterOfMass
  
  
  
!   function MecanicProperties_getUnitaryMolecularInertiaTensor( this, cartesianCoordinatesInCM,info ) result( output )
!     implicit  none
!     type(MecanicProperties) :: this
!     type(Matrix) :: output
!     type(Vector), optional :: cartesianCoordinatesInCM
!     integer,optional :: info
    
!     real(8) :: unitaryCenterOfMass(3)
!     type(vector) :: coordinates
!     integer :: infoProcess
!     integer :: i
!     integer :: j
    
!     call Vector_constructor(coordinates,3*this%numberOfPoints)
    
!     !!**************************************************************************
!     !!   Calcula el centro de masa unitario y desplaza la coordenadas
!     !!*************
!     unitaryCenterOfMass=MecanicProperties_getUnitaryCenterOfMass(this)
!     do i=1,this%numberOfPoints
!        coordinates%values(3*i-2:3*i) = this%cartesianCoordinates%values(i,:) - unitaryCenterOfMass
!     end do
    
!     !! Si se desea retorna una matriz de coordenadas centras en el centro de masa
!     if( present(cartesianCoordinatesInCM) ) cartesianCoordinatesInCM=coordinates
    
!     !!
!     !!**************************************************************************
    
!     !!**************************************************************************
!     !!  Determina el tesor de inercial para una estructira de masa unitarias
!     !!************
!     call Matrix_constructor(output,3_8,3_8,0.0_8)
!     do i=1,this%numberOfPoints
!        j=3*(i-1)+1
!        output%values(1,1)=output%values(1,1)+coordinates%values(j+1)**2+coordinates%values(j+2)**2
!        output%values(1,2)=output%values(1,2)-coordinates%values(j)*coordinates%values(j+1)
!        output%values(1,3)=output%values(1,3)-coordinates%values(j)*coordinates%values(j+2)
!        output%values(2,2)=output%values(2,2)+coordinates%values(j)**2+coordinates%values(j+2)**2
!        output%values(2,3)=output%values(2,3)-coordinates%values(j+1)*coordinates%values(j+2)
!        output%values(3,3)=output%values(3,3)+coordinates%values(j)**2+coordinates%values(j+1)**2
!     end do
!     output%values(2,1)=output%values(1,2)
!     output%values(3,1)=output%values(1,3)
!     output%values(3,2)=output%values(2,3)
    
!     output = Matrix_inverse( output, info=infoProcess,printFormatFlags=WITHOUT_MESSAGES )
!     call Vector_destructor(coordinates)
!     if (present(info)) info=infoProcess
    
!     !!
!     !!**************************************************************************
    
!   end function MecanicProperties_getUnitaryMolecularInertiaTensor
  
!   !>
!   !! @brief Calcula y retorna una proyector de constantes de fuerza para
!   !!		realizar traslaciones y rotaciones  infinetesimales de la molecula
!   function MecanicProperties_getForceConstantsProjector( this ) result( output )
!     implicit none
!     type(MecanicProperties) :: this
!     type(Matrix) :: output
    
!     integer :: i
!     integer :: j
!     integer :: index_x
!     integer :: index_y
!     integer :: index_z
!     integer :: aux
!     logical :: isNull
!     real(8) :: sqrtMass
!     real(8) :: coordinatesProyected(3)
!     real(8) :: geometryInCM(3)
!     real(8) :: squareNorm
!     type(LinearVectorialSpace) :: spaceOfForceConstants
    
!     call Matrix_constructor( output, int(this%numberOfPoints*3,8), int(this%numberOfPoints*3,8), 0.0_8 )
!     this%molecularInertiaTensor = MecanicProperties_getMolecularInertiaTensor( this )
    
!     do i=1, this%numberOfPoints
       
!        index_x = 3*i - 2
!        index_y = 3*i - 1
!        index_z = 3*i
       
!        sqrtMass= sqrt( this%puntualMasses(i) )
!        geometryInCM=this%cartesianCoordinates%values(i,:)-this%centerOfMass(:)
       
!        !!
!        !! Proyecta las coordenadas cartesianas sobre el tensor de inercia molecular
!        !!
!        coordinatesProyected(1)=dot_product(geometryInCM,this%molecularInertiaTensor%values(1,:) )
!        coordinatesProyected(2)=dot_product(geometryInCM,this%molecularInertiaTensor%values(2,:) )
!        coordinatesProyected(3)=dot_product(geometryInCM,this%molecularInertiaTensor%values(3,:) )
       
!        output%values(index_x,1) = sqrtMass
!        output%values(index_y,2) = sqrtMass
!        output%values(index_z,3) = sqrtMass
       
!        output%values(index_x,4) = 	(coordinatesProyected(2)*this%molecularInertiaTensor%values(1,3) &
!             - coordinatesProyected(3)*this%molecularInertiaTensor%values(1,2) )/sqrtMass
!        output%values(index_y,4) = 	(coordinatesProyected(2)*this%molecularInertiaTensor%values(2,3) &
!             - coordinatesProyected(3)*this%molecularInertiaTensor%values(2,2) )/sqrtMass
!        output%values(index_z,4) = 	(coordinatesProyected(2)*this%molecularInertiaTensor%values(3,3) &
!             - coordinatesProyected(3)*this%molecularInertiaTensor%values(3,2) )/sqrtMass
       
!        output%values(index_x,5) = 	(coordinatesProyected(3)*this%molecularInertiaTensor%values(1,1) &
!             - coordinatesProyected(1)*this%molecularInertiaTensor%values(1,3) )/sqrtMass
!        output%values(index_y,5) = 	(coordinatesProyected(3)*this%molecularInertiaTensor%values(2,1) &
!             - coordinatesProyected(1)*this%molecularInertiaTensor%values(2,3) )/sqrtMass
!        output%values(index_z,5) = 	(coordinatesProyected(3)*this%molecularInertiaTensor%values(3,1) &
!             - coordinatesProyected(1)*this%molecularInertiaTensor%values(3,3) )/sqrtMass
       
!        output%values(index_x,6) =	(coordinatesProyected(1)*this%molecularInertiaTensor%values(1,2) &
!             - coordinatesProyected(2)*this%molecularInertiaTensor%values(1,1) )/sqrtMass
!        output%values(index_y,6) = 	(coordinatesProyected(1)*this%molecularInertiaTensor%values(2,2) &
!             - coordinatesProyected(2)*this%molecularInertiaTensor%values(2,1) )/sqrtMass
!        output%values(index_z,6) =	(coordinatesProyected(1)*this%molecularInertiaTensor%values(3,2) &
!             - coordinatesProyected(2)*this%molecularInertiaTensor%values(3,1) )/sqrtMass
!     end do
    
!     !!***********************************************************************
!     !! Verifica que los 6 vectores correspcondan verdaderamente a modos normales
!     !! rotacionales y traslacionales
!     !!
!     this%numberOfRotationalAndTranslationalVars = 0
!     isNull=.false.
!     aux = 0
    
!     do i=1,6
       
!        squareNorm=dot_product( output%values(:,i),output%values(:,i) )
       
!        if ( squareNorm > 1.0D-6) then
          
!           output%values(:,i) = output%values(:,i) / sqrt( squareNorm )
!           this%numberOfRotationalAndTranslationalVars = this%numberOfRotationalAndTranslationalVars + 1
          
!           if ( isNull ) then
!              output%values(:,i-aux) = output%values(:,i)
!              output%values(:,i) = 0.0_8
!           end if
!        else
          
!           isNull = .true.
!           aux = aux+1
          
!        end if
       
!     end do
!     !!
!     !!***********************************************************************
!     this%numberOfVibrationalVars = 3*this%numberOfPoints - this%numberOfRotationalAndTranslationalVars
!     !!Adiciona una serie de vectores asociados a vibraciones con el fin
!     !! de completar la matriz de transformacion
!     j=1
!     do i=this%numberOfRotationalAndTranslationalVars+1,3*this%numberOfPoints
!        output%values(j,i)=1.0_8
!        j=j+1
!     end do
    
!     !! Construye un espacio vectorial lineal con los vectores previemente generados
!     call LinearVectorialSpace_constructor( spaceOfForceConstants, mmatrix= output%values )
    
!     !!**********
!     !! 	Proyecta los vectores asociados a grados de libertad vibracionales sobre los
!     !! 	asociados a grados de libertad rotacionales y traslacionales. - Los primeros
!     !!	N vectores (numberOfRotationalAndTranslationalVars ),  se asocian a los grados
!     !!	de libertad rotacionales y translacionales el resto se asocian a los grados de libertad
!     !!	vibracionales -
!     !!**
!     call LinearVectorialSpace_orthogonalize( spaceOfForceConstants )
!     output = LinearVectorialSpace_getElements(spaceOfForceConstants)
!     call LinearVectorialSpace_destructor( spaceOfForceConstants )
    
!     !!
!     !! Reordena el proyector colocando los vectores asociados al vibraciones al principio y los
!     !! asociados a rotaciones y traslaciones al final
!     !!
!     call Matrix_swapBlockOfColumns( output, [this%numberOfRotationalAndTranslationalVars+1, 3*this%numberOfPoints ] )
    
!   end function MecanicProperties_getForceConstantsProjector
  
!   !>
!   !! @brief Calcula y retorna el proyector de Handy
!   function MecanicProperties_getHandyProjector( this,infoProcess ) result( output )
!     implicit none
!     type(MecanicProperties) :: this
!     integer,optional,intent(inout) :: infoProcess
!     type(Matrix) :: output
    
!     type(vector) :: coordinates
!     type(Matrix) :: molecularInertiaTensor
!     real(8) :: totallyAsymmetricTensor(3,3,3)
!     real(8) :: L_ilk
!     integer  :: i,j,IP,INDX,KNDX,JP,JNDX,LNDX,IC,II,KK,JEND,JC,JJ,LL,IA,IB,JA,JB
    
!     totallyAsymmetricTensor=MecanicProperties_getTotallyAsymmetricTensor(this)
    
!     call Vector_constructor(coordinates,3*this%numberOfPoints)
!     call Matrix_constructor(output,3_8*this%numberOfPoints,3_8*this%numberOfPoints,0.0_8)
!     call Matrix_constructor(molecularInertiaTensor,3_8,3_8,0.0_8)
    
!     !! Obtiene tensor de inercia molecular unitaria y la matriz de coodenadas en el centro de masa
!     molecularInertiaTensor=MecanicProperties_getUnitaryMolecularInertiaTensor(this,coordinates,info=infoProcess)
    
!     !!**************************************************************************
!     !! Calcula el proyector de Miller-Handy-Adams
!     !!************
!     do IP=1,this%numberOfPoints
!        INDX=3*(IP-1)
!        KNDX=MAX(3*(IP-1),6*(IP-1)-3*this%numberOfPoints)
!        do JP=1,IP
!           JNDX=3*(JP-1)
!           LNDX=MAX(3*(JP-1),6*(JP-1)-3*this%numberOfPoints)
!           do IC=1,3
!              II=INDX+IC
!              KK=KNDX+IC
!              JEND=3
!              if(JP == IP) JEND=IC
!              do JC=1,JEND
!                 JJ=JNDX+JC
!                 LL=LNDX+JC
!                 L_ilk=0.0
!                 do IA=1,3
!                    do IB=1,3
!                       if(abs(totallyAsymmetricTensor(IA,IB,IC)) ==1.0 ) then
!                          do JA=1,3
!                             do JB=1,3
!                                if ( abs(totallyAsymmetricTensor(JA,JB,JC)) == 1.0)  then
!                                   L_ilk = L_ilk +totallyAsymmetricTensor(IA,IB,IC) &
!                                        * totallyAsymmetricTensor(JA,JB,JC) &
!                                        * molecularInertiaTensor%values(IA,JA) &
!                                        * coordinates%values(INDX+IB) &
!                                        * coordinates%values(JNDX+JB)
!                                end if
!                             end do
!                          end do
!                       end if
!                    end do
!                 end do
!                 output%values(KK,LL)=L_ilk
!                 if (IC.EQ.JC) output%values(KK,LL)=output%values(KK,LL)  + 1.0 /this%numberOfPoints
!              end do
!           end do
!        end do
!     end do
!     !!
!     !!**************************************************************************
    
!     !! Calcula el termino I-P
!     do i=1,3*this%numberOfPoints
!        do j=1,i
!           output%values(i,j)=-output%values(i,j)
!           if(i == j) output%values(i,j) = 1.0 + output%values(i,j)
!           if( abs(output%values(i,j) ) < 1.0D-8) output%values(i,j)= 0.0
!           output%values(j,i) = output%values(i,j)
!        end do
!     end do
    
!     call Vector_destructor(coordinates)
!     call Matrix_destructor(molecularInertiaTensor)
    
!   end function MecanicProperties_getHandyProjector
  
!   !>
!   !!  @brief Calcula el core de la matriz de transformacion para eliminar  grados de libertad externos
!   !!		del gradiente y de la hessiana.
!   function MecanicProperties_getTransformationMatrix( this ) result(output)
!     implicit none
!     type(MecanicProperties) :: this
!     type(Matrix) :: output
    
!     call Matrix_constructor( output, int(this%numberOfPoints*3,8), int(this%numberOfPoints*3,8), 0.0_8 )
!     this%forceConstansProjector = MecanicProperties_getForceConstantsProjector(this)
    
!     output%values= -1.0_8 * matmul( this%forceConstansProjector%values(:,1:this%numberOfVibrationalVars+1), &
!          transpose(this%forceConstansProjector%values(:,1:this%numberOfVibrationalVars+1) ) )
    
!   end function MecanicProperties_getTransformationMatrix
  
!   !>
!   !! @brief Remueve los grados de libertad externo del vector de gradiente
!   function MecanicProperties_getGradientProjectedOnExternalDegrees( this, gradientVector ) result( output )
!     implicit none
!     type( MecanicProperties ) :: this
!     type(Vector) :: gradientVector
!     type(Vector) :: output
    
!     integer :: i
!     type(Matrix) :: auxTransformationMatrix
    
!     this%transformationMatrix=MecanicProperties_getTransformationMatrix( this )
!     auxTransformationMatrix=this%transformationMatrix
!     call Matrix_symmetrize( auxTransformationMatrix, "L" )
    
!     do i=1,size(auxTransformationMatrix%values,dim=1)
!        auxTransformationMatrix%values(i,i) = 1.0_8 + auxTransformationMatrix%values(i,i)
!     end do
    
!     call Vector_constructor(output,size(gradientVector%values))
    
!     !!
!     !! Realiza la proyeccion de los grados de libertad externos del
!     !! Vector de gradiente cartesiano
!     !!
!     do i=1, size(gradientVector%values)
!        output%values(i) = dot_product( auxTransformationMatrix%values(i,:),gradientVector%values )
!     end do
    
!     call Matrix_destructor( auxTransformationMatrix )
    
!   end function MecanicProperties_getGradientProjectedOnExternalDegrees
  
!   !>
!   !! @brief Remueve los grados de libertad externos de la matriz Hessiana
!   function MecanicProperties_getHessianeProjectedOnExternalDegrees( this, hessianeMatrix ) result( output )
!     implicit none
!     type( MecanicProperties ) :: this
!     type(Matrix) :: hessianeMatrix
!     type(Matrix) :: output
    
!     integer :: i
!     type(Matrix) :: auxTransformationMatrix
    
!     this%transformationMatrix=MecanicProperties_getTransformationMatrix( this )
!     auxTransformationMatrix=this%transformationMatrix
    
!     do i=1,size(auxTransformationMatrix%values,dim=1)
!        auxTransformationMatrix%values(i,i) = 1.0_8 + auxTransformationMatrix%values(i,i)
!     end do
    
!     output=hessianeMatrix
!     output%values = matmul( hessianeMatrix%values, auxTransformationMatrix%values)
    
!     !!Proyecta las grados de libertad vibracionales sobre los grados de libertad externos
!     output%values=matmul(transpose(auxTransformationMatrix%values),output%values)
    
!     call Matrix_destructor( auxTransformationMatrix )
    
!   end function MecanicProperties_getHessianeProjectedOnExternalDegrees
  
  !>
  !! @brief  Maneja excepciones de la clase
  subroutine MecanicProperties_exception( typeMessage, description, debugDescription)
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
    
  end subroutine MecanicProperties_exception
  
end module MecanicProperties_
