!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://sites.google.com/a/bt.unal.edu.co/andresreyes/home
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!  Prof. G. MERINO's Lab. Universidad de Guanajuato
!!    http://quimera.ugto.mx/qtc/gmerino.html
!!
!!  Authors:
!!    E. F. Posada (efposadac@unal.edu.co)
!!    R. Flores (roberto.floresmoreno.qt@gmail.com)
!!
!!  Contributors:
!!
!!    Todos los derechos reservados, 2011
!!
!!******************************************************************************

module OutputBuilder_
  use Exception_
  use MolecularSystem_
  use CalculateWaveFunction_
  use ParticleManager_
  use BasisSet_
  use Matrix_
  use Vector_
  use String_
  implicit none

  !>
  !! @brief Description
  !!
  !! @author felix
  !!
  !! <b> Creation data : </b> 08-04-11
  !!
  !! <b> History change: </b>
  !!
  !!   - <tt> 08-04-11 </tt>:  felix ( email@server )
  !!        -# description.
  !!   - <tt> 10-31-2014 </tt>:  Mauricio Rodas ( jmrodasr@unal.edu.co )
  !!        -# Adapts this module to Lowdin2 to generate the molden input
  !!   - <tt> 10-31-2014 </tt>:  Jorge Charry ( jacharry@unal.edu.co )
  !!        -# Adapts fully this module to Lowdin2 
  !!   - <tt> 04-20-2015 </tt>:  Jorge Charry ( jacharry@unal.edu.co )
  !!        -# Reorder the coefficients matrix according to the molden format
  !!   - <tt> MM-DD-YYYY </tt>:  authorOfChange ( email@server )
  !!        -# description
  !!
  !<
  type, public :: OutputBuilder
     character(50) :: type
     character(50) :: species
     character(100),allocatable :: fileName(:)
     character(100) :: fileName2
     character(50) :: axisLabel(3)
     character(10) :: wavefunctionType
     integer :: state
     integer :: orbital
     integer :: dimensions
     integer :: outputID
     integer :: auxID
     integer :: pointsPerDim(3)
     real(8) :: cubeSize
     real(8) :: minValue
     real(8) :: maxValue
     type(vector) :: point1
     type(vector) :: point2
     type(vector) :: point3
     type(vector) :: step1
     type(vector) :: step2
  end type OutputBuilder

  type(OutputBuilder), public, allocatable :: outputs_instance(:)

  public :: &
       OutputBuilder_constructor, &
       OutputBuilder_destructor, &
       OutputBuilder_show, &
       OutputBuilder_writeMoldenFile, &
       OutputBuilder_VecGamessFile, &
       OutputBuilder_writeEigenvalues,&
       OutputBuilder_generateAIMFiles, &
       OutputBuilder_generateExtendedWfnFile, &
       OutputBuilder_buildOutput, &
       OutputBuilder_make2DGnuplot, &
       OutputBuilder_make3DGnuplot, &
       OutputBuilder_getPlot, &
       OutputBuilder_getCube, &
       OutputBuilder_casinoFile

  private   

  interface 

     subroutine Molden2AIM(inputFileName,totalEnergy,virial)
       implicit none  
       character(50) :: inputFileName
       real(8) :: totalEnergy, virial
     end subroutine Molden2AIM

  end interface

contains


  !>
  !! @brief Constructor por omision
  !!
  !! @param this
  subroutine OutputBuilder_constructor(this, ID, &
       type, &
       species, &
       plane, &
       axis, &
       state, &
       orbital, &
       dimensions, &
       pointsPerDim, &
       scanStep, &
       cubeSize, &
       minValue, &
       maxValue, &
       offsetX, &
       offsetY, &
       offsetZ, &
       limitX, &
       limitY, &
       limitZ, &
       center, &
       point1, &
       point2, &
       point3)

    type(OutputBuilder) :: this
    integer :: ID
    character(*) :: type
    character(*) :: species
    character(*),optional :: plane
    character(*),optional :: axis
    integer,optional :: state
    integer,optional :: orbital
    integer,optional :: dimensions
    integer,optional :: pointsPerDim
    real(8),optional :: scanStep
    real(8),optional :: cubeSize
    real(8),optional :: minValue
    real(8),optional :: maxValue
    real(8),optional :: offsetX
    real(8),optional :: offsetY
    real(8),optional :: offsetZ
    real(8),optional :: limitX(2)
    real(8),optional :: limitY(2)
    real(8),optional :: limitZ(2)
    real(8),optional :: center(3)
    real(8),optional :: point1(3)
    real(8),optional :: point2(3)
    real(8),optional :: point3(3)

    integer :: i
    real(8) :: auxX, auxY, auxZ, auxReal, auxStep
    real(8) :: auxLimX(2), auxLimY(2), auxLimZ(2)
    character(50) :: auxString
    logical :: existFile

    this%type=type
    ! print *, "this%type", this%type
    this%outputID=ID
    this%species=trim(String_getUppercase(species))
    if( trim(this%species) .eq. "ALL") then
       allocate(this%fileName(MolecularSystem_getNumberOfQuantumSpecies()))
    else
       allocate(this%fileName(1))
    end if

    this%state=1
    this%orbital=0
    this%dimensions=0
    this%cubeSize=10
    this%axisLabel(1:3)=""
    this%minValue=0.0_8
    this%minValue=0.0_8
    call Vector_constructor(this%point1, 3, 0.0_8 )
    call Vector_constructor(this%point2, 3, 0.0_8 )
    call Vector_constructor(this%point3, 3, 0.0_8 )
    call Vector_constructor(this%step1, 3, 0.0_8 )
    call Vector_constructor(this%step2, 3, 0.0_8 )

    if( present(state)) this%state=state
    if( present(orbital)) this%orbital=orbital
    if( present(dimensions)) this%dimensions=dimensions
    if( present(cubeSize)) this%cubeSize=cubeSize
    if( present(pointsPerDim)) this%pointsPerDim(1:3)=pointsPerDim
    if( present(minValue)) this%minValue=minValue
    if( present(maxValue)) this%maxValue=maxValue
    if( present(point1)) this%point1%values=point1
    if( present(point2)) this%point2%values=point2
    if( present(point3)) this%point3%values=point3

    if (this%pointsPerDim(1) .eq. 0) this%pointsPerDim(:)=CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION

    auxString=""
    auxX=0.0
    auxY=0.0
    auxZ=0.0
    auxStep=0.0
    auxReal=0.0

    if( present(offsetX)) auxX=offsetX
    if( present(offsetY)) auxY=offsetY
    if( present(offsetZ)) auxZ=offsetZ
    if( present(limitX)) auxLimX=limitX
    if( present(limitY)) auxLimY=limitY
    if( present(limitZ)) auxLimZ=limitZ

    if(present(plane)) auxString=plane
    if(present(scanStep)) auxStep=scanStep

    if(auxString .ne. "") then
       this%dimensions=3
       if(trim(auxString) .eq. "xy" .or. trim(auxString) .eq. "yx") then
          this%axisLabel(1)="X"
          this%axisLabel(2)="Y"
          this%point1%values(1)=auxLimX(1)
          this%point1%values(2)=auxLimY(1)
          this%point1%values(3)=auxZ

          this%point2%values(1)=auxLimX(2)
          this%point2%values(2)=auxLimY(1)          
          this%point2%values(3)=auxZ

          this%point3%values(1)=auxLimX(1)
          this%point3%values(2)=auxLimY(2)          
          this%point3%values(3)=auxZ
       else if(trim(auxString) .eq. "xz" .or. trim(auxString) .eq. "zx") then
          this%axisLabel(1)="X"
          this%axisLabel(2)="Z"
          this%point1%values(1)=auxLimX(1)
          this%point1%values(2)=auxY
          this%point1%values(3)=auxLimZ(1)

          this%point2%values(1)=auxLimX(2)
          this%point2%values(2)=auxY
          this%point2%values(3)=auxLimZ(1)          

          this%point3%values(1)=auxLimX(1)
          this%point3%values(2)=auxY
          this%point3%values(3)=auxLimZ(2)          
       else if(trim(auxString) .eq. "yz" .or. trim(auxString) .eq. "zy") then
          this%axisLabel(1)="Y"
          this%axisLabel(2)="Z"
          this%point1%values(1)=auxX
          this%point1%values(2)=auxLimY(1)
          this%point1%values(3)=auxLimZ(1)

          this%point2%values(1)=auxX
          this%point2%values(2)=auxLimY(2)
          this%point2%values(3)=auxLimZ(1)          

          this%point3%values(1)=auxX
          this%point3%values(2)=auxLimY(1)
          this%point3%values(3)=auxLimZ(2)          
       else
          call Exception_stopError("Please select a plane (xy,xz or yz) to build the plot", "OutputBuilder_constructor" )
       end if
    end if

    if(present(axis)) auxString=axis
    if(auxString .ne. "") then
       this%dimensions=2
       select case(trim(auxString))
       case ( "x")
          this%axisLabel(1)="X"
          this%point1%values(1)=auxLimX(1)
          this%point2%values(1)=auxLimX(2)          
          this%point1%values(2)=auxY
          this%point2%values(2)=auxY
          this%point1%values(3)=auxZ
          this%point2%values(3)=auxZ
       case ( "y") 
          this%axisLabel(1)="Y"
          this%point1%values(1)=auxX
          this%point2%values(1)=auxX
          this%point1%values(2)=auxLimY(1)
          this%point2%values(2)=auxLimY(2)          
          this%point1%values(3)=auxZ
          this%point2%values(3)=auxZ
       case ( "z") 
          this%axisLabel(1)="Z"
          this%point1%values(1)=auxX
          this%point2%values(1)=auxX
          this%point1%values(2)=auxY
          this%point2%values(2)=auxY
          this%point1%values(3)=auxLimZ(1)
          this%point2%values(3)=auxLimZ(2)          
       case default
          call Exception_stopError( "Please select an axis (x,y or z) to build the plot", "OutputBuilder_constructor" )
       end select

    end if

    if(auxStep .gt. 0.0) then
       if(this%dimensions .eq. 2 .or. this%dimensions .eq. 3) then
          auxReal=sqrt(sum((this%point2%values(:)-this%point1%values(:))*(this%point2%values(:)-this%point1%values(:))))
          this%step1%values(:)=(this%point2%values(:)-this%point1%values(:))/auxReal*auxStep
          this%pointsPerDim(1)=int(auxReal/auxStep)
       end if
       if(this%dimensions .eq. 3) then
          auxReal=sqrt(sum((this%point3%values(:)-this%point1%values(:))*(this%point3%values(:)-this%point1%values(:))))
          this%step2%values(:)=(this%point3%values(:)-this%point1%values(:))/auxReal*auxStep
          this%pointsPerDim(2)=int(auxReal/auxStep)
       end if
    else
       if(this%dimensions .eq. 2 .or. this%dimensions .eq. 3) &
            this%step1%values(:)=(this%point2%values(:)-this%point1%values(:))/this%pointsPerDim(1)
       if(this%dimensions .eq. 3) &
            this%step2%values(:)=(this%point3%values(:)-this%point1%values(:))/this%pointsPerDim(2)
    end if

    if(auxStep .gt. 0.0 .and. this%cubeSize .gt. 0.0) this%pointsPerDim(1:3)=int((this%cubeSize*2.0)/auxStep)

    if( present(center) ) auxReal=sum(center*center)
    if(auxReal .gt. 0.0) this%point1%values=center

    if ( trim(CONTROL_instance%UNITS) == "ANGS") then
       this%point1%values= this%point1%values/ANGSTROM
       this%point2%values= this%point2%values/ANGSTROM
       this%point3%values= this%point3%values/ANGSTROM
       this%cubeSize=this%cubeSize/ANGSTROM
       this%step1%values= this%step1%values/ANGSTROM
       this%step2%values= this%step2%values/ANGSTROM
    end if

    !! By default, we work with HF-KS wavefunctions
    this%wavefunctionType="HF"
    !! Check if there are CI density matrices
    inquire(FILE = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci", EXIST = existFile )
    if(existFile .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0) this%wavefunctionType="CI"

    this%auxID=1
    !!Check for other outputs of the same type
    do i=1, this%outputID-1
       if( trim(outputs_instance(i)%type) .eq. trim(this%type) .and. &
            trim(outputs_instance(i)%species) .eq. trim(this%species) .and. &
            outputs_instance(i)%dimensions .eq. this%dimensions .and. &
            outputs_instance(i)%orbital .eq. this%orbital) this%auxID=this%auxID+1
    end do

  end subroutine OutputBuilder_constructor


  !>
  !! @brief Destructor por omision
  !!
  !! @param this
  !<
  subroutine OutputBuilder_destructor(this)
    implicit none
    type(OutputBuilder) :: this

    call Vector_destructor(this%point1)
    call Vector_destructor(this%point2)
    call Vector_destructor(this%point3)

  end subroutine OutputBuilder_destructor

  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine OutputBuilder_show(this)
    implicit none
    type(OutputBuilder) :: this
    integer :: l

    print *, "--------------------------------------------------------"
    write (*,"(A20,I5,T2,A18)") "Output Number: ", this%outputID, this%type

    ! TODO Fix this line.
    ! if (this%filename2 /= "") print *, "FileName 2: ", this%fileName2
    if (this%species /= "ALL") then
       write (*,"(A20,A)") "for species: ", trim(this%species)
    else
       write (*,"(T20,A)") "for all species "
    end if
    if (this%state /= 1) write (*,"(A20,I10)") "for excited state: ", this%state

    select case(trim(this%type))

    case ( "MOLDENFILE") 

    case ("VECGAMESSFILE")

    case ("CASINOFILE")

    case ("EIGENGAMESSFILE")

    case ("FCHKFILE")

    case ( "WFNFILE") 

    case ( "NBO47FILE") 

    case ( "WFXFILE" ) 

    case ( "EXTENDEDWFNFILE") 

    case ( "DENSITYPLOT") 
       write (*,"(A20,I10)") "dimensions: ", this%dimensions
       write (*,"(A20,F10.5,F10.5,F10.5)") "Point 1 (a.u.): ", this%point1%values(1), this%point1%values(2), this%point1%values(3)
       write (*,"(A20,F10.5,F10.5,F10.5)") "Point 2 (a.u.): ", this%point2%values(1), this%point2%values(2), this%point2%values(3)
       if (this%dimensions .eq. 2) then
          write (*,"(A20,F10.5)") "Step size: ", sqrt(sum(this%step1%values*this%step1%values))
          write (*,"(A20,I10)") "No. steps: ", this%pointsPerDim(1)
       end if
       if (this%dimensions .eq. 3) then
          write (*,"(A20,F10.5,F10.5,F10.5)") "Point 3 (a.u.): ", this%point3%values(1), this%point3%values(2), this%point3%values(3)          
          write (*,"(A20,F10.5,F10.5)") "Step sizes: ", sqrt(sum(this%step1%values*this%step1%values)), sqrt(sum(this%step2%values*this%step2%values))
          write (*,"(A20,2I10)") "No. steps: ", this%pointsPerDim(1), this%pointsPerDim(2)
       end if
    case ( "DENSITYCUBE") 
       write (*,"(A20,F10.5)") "cube size (a.u.): ", this%cubeSize
       write (*,"(A20,3F10.5)") "cube center (a.u.): ", this%point1%values(1:3)
       write (*,"(A20,3I10)") "No. steps: ", this%pointsPerDim(1:3)

    case ( "ORBITALPLOT") 
       if(this%orbital .eq. 0) then
          write (*,"(A40)") "for the highest occupied orbital"
       else
          write (*,"(A20,I10)") "for orbital: ", this%orbital
       end if
       write (*,"(A20,I10)") "dimensions: ", this%dimensions
       write (*,"(A20,F10.5,F10.5,F10.5)") "Point 1 (a.u.): ", this%point1%values(1), this%point1%values(2), this%point1%values(3)
       write (*,"(A20,F10.5,F10.5,F10.5)") "Point 2 (a.u.): ", this%point2%values(1), this%point2%values(2), this%point2%values(3)
       if (this%dimensions .eq. 2) then
          write (*,"(A20,F10.5)") "Step size: ", sqrt(sum(this%step1%values*this%step1%values))
          write (*,"(A20,I10)") "No. steps: ", this%pointsPerDim(1)
       end if
       if (this%dimensions .eq. 3) then
          write (*,"(A20,F10.5,F10.5,F10.5)") "Point 3 (a.u.): ", this%point3%values(1), this%point3%values(2), this%point3%values(3)
          write (*,"(A20,F10.5,F10.5)") "Step sizes: ", sqrt(sum(this%step1%values*this%step1%values)), sqrt(sum(this%step2%values*this%step2%values))
          write (*,"(A20,2I10)") "No. steps: ", this%pointsPerDim(1), this%pointsPerDim(2)
       end if

    case ( "ORBITALCUBE") 
       if(this%orbital .eq. 0) then
          write (*,"(A40)") "for the highest occupied orbital"
       else
          write (*,"(A20,I10)") "for orbital: ", this%orbital
       end if
       write (*,"(A20,F10.5)") "cube size (a.u.): ", this%cubeSize
       write (*,"(A20,3F10.5)") "cube center (a.u.): ", this%point1%values(1:3)
       write (*,"(A20,3I10)") "No. steps: ", this%pointsPerDim(1:3)

    case default

    end select
    do l=1,size(this%fileName)
       write (*,"(A20,A)") "FileName: ", this%fileName(l)
    end do

    print *, "--------------------------------------------------------"
    print *, ""

  end subroutine OutputBuilder_show


  !>
  !! @brief Muestra informacion del objeto
  !!
  !! @param this 
  !<
  subroutine OutputBuilder_buildOutput(this)
    implicit none
    type(OutputBuilder) :: this

    select case( this%type )

    case ( "MOLDENFILE") 
       call OutputBuilder_writeMoldenFile (this)

    case ("VECGAMESSFILE")
       call OutputBuilder_VecGamessFile (this)

    case ("CASINOFILE")
       call OutputBuilder_casinoFile (this)

    case ("EIGENGAMESSFILE")
       call OutputBuilder_writeEigenvalues (this)

    case ("FCHKFILE")
       call OutputBuilder_writeFchkFile (this)

    case ( "WFNFILE") 
       call OutputBuilder_writeMoldenFile (this)
       call OutputBuilder_generateAIMFiles (this)

    case ( "NBO47FILE") 
       call OutputBuilder_writeMoldenFile (this)
       call OutputBuilder_generateAIMFiles (this)

    case ( "WFXFILE" ) 

       call OutputBuilder_writeMoldenFile (this)
       call OutputBuilder_generateAIMFiles (this)

    case ( "EXTENDEDWFNFILE") 
       call OutputBuilder_writeMoldenFile (this)
       call OutputBuilder_generateAIMFiles (this)
       call OutputBuilder_generateExtendedWfnFile (this)

    case ( "DENSITYPLOT") 
       if(this%maxValue .eq. 0.0) this%maxValue=0.5
       call OutputBuilder_getPlot(this)

    case ( "DENSITYCUBE") 
       call OutputBuilder_getCube(this)
       !
    case ( "ORBITALPLOT") 
       if(this%maxValue .eq. 0.0 .and. this%minValue .eq. 0.0) then
          this%maxValue=1.0
          this%minValue=-1.0
       end if
       call OutputBuilder_getPlot(this)
       !
    case ( "ORBITALCUBE") 
       call OutputBuilder_getCube(this)
       !
       !     case ( "fukuiPlot") 
       !        if (this%dimensions == 2) call OutputBuilder_get2DPlot(this)
       !        if (this%dimensions == 3) call OutputBuilder_get3DPlot(this)
       !
       !     case ( "fukuiCube") 
       !        call OutputBuilder_getCube(this)
       !
    case default
       call Exception_stopError("The output type "//this%type//" you requested has not been implemented yet", "OutputBuilder_buildOutput" )

    end select
  end subroutine OutputBuilder_buildOutput

  subroutine OutputBuilder_writeMoldenFile(this)
    implicit none
    type(OutputBuilder) :: this

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: numberOfSpecies
    integer :: state,numberOfStates
    integer :: occupationTotal
    logical :: wasPress
    character(50) :: auxString
    character(10) :: symbol
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Matrix),allocatable :: coefficientsOfcombination(:,:)
    type(Vector),allocatable :: energyOfMolecularOrbital(:,:)
    type(Vector),allocatable :: fractionalOccupations(:,:)
    character(10),allocatable :: labels(:)
    integer :: numberOfContractions
    integer :: totalNumberOfParticles, n

    !     if ( CONTROL_instance%ARE_THERE_DUMMY_ATOMS ) then
    !        auxString=MolecularSystem_getNameOfSpecies( 1 )
    !        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"mol"
    !        open(10,file=this%fileName,status='replace',action='write')
    !        write (10,"(A)") "Hola soy un archivo de molden"
    !        close(10)

    !     else

    localizationOfCenters=ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    auxMatrix=localizationOfCenters
    allocate( labels( size(auxMatrix%values,dim=1) ) )
    allocate( charges( size(auxMatrix%values,dim=1) ) )
    labels=ParticleManager_getLabelsOfCentersOfOptimization()
    charges=ParticleManager_getChargesOfCentersOfOptimization()
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    numberOfStates=1
    if( this%wavefunctionType .eq. "CI") then
       write (*,"(A50)") "We are printing molden files for the CI states!"
       numberOfStates=CONTROL_instance%CI_STATES_TO_PRINT
    else ! (this%wavefunctionType .eq. "HF")
       numberOfStates=1
    end if
    allocate(fractionalOccupations(numberOfSpecies,numberOfStates))
    allocate(coefficientsOfCombination(numberOfSpecies,numberOfStates))
    allocate(energyOfMolecularOrbital(numberOfSpecies,numberOfStates))

    call CalculateWaveFunction_loadCoefficientsMatrices ( numberOfSpecies, numberOfStates, this%wavefunctionType, coefficientsOfCombination, fractionalOccupations, energyOfMolecularOrbital)


    do state=1,numberOfStates
       do l=1,numberOfSpecies

          if (state .eq. 1) then
             auxString=MolecularSystem_getSymbolOfSpecies( l )
          else
             write(auxString, "(I8)")  state
             auxString=trim(MolecularSystem_getSymbolOfSpecies( l ))//"-"//trim( adjustl(auxString))
          end if

          this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".molden"

          totalNumberOfParticles = 0

          open(10,file=this%fileName(l),status='replace',action='write')
          write(10,"(A)") "[Molden Format]"
          if ( CONTROL_instance%UNITS=="ANGS") then
             write(10,"(A)") "[Atoms] Angs"
          else 
             write(10,"(A)") "[Atoms] AU"
          end if

          auxMatrix%values=0.0
          j=0
          do i=1, size(MolecularSystem_instance%species(l)%particles)
             j=j+1
             origin = MolecularSystem_instance%species(l)%particles(j)%origin 
             auxMatrix%values(j,:)=origin
             symbol=MolecularSystem_instance%species(l)%particles(j)%nickname

             if(scan(symbol,"_") /=0) symbol=symbol(1:scan(symbol,"_")-1)
             if(scan(symbol,"[") /=0) symbol=symbol(scan(symbol,"[")+1:scan(symbol,"]")-1)

             if ( CONTROL_instance%UNITS=="ANGS") origin = origin * ANGSTROM

             totalNumberOfParticles = totalNumberOfParticles + 1

             write (10,"(A,I8,I8,F15.8,F15.8,F15.8)") trim(symbol), j,&
                  int(abs(molecularSystem_instance%allParticles( MolecularSystem_instance%species(l)%particles(j)%owner )%particlePtr%charge)) ,&
                  origin(1), origin(2), origin(3)
             ! int(abs(MolecularSystem_instance%species(l)%particles(j)%totalCharge)), &

          end do


          if ( CONTROL_instance%MOLDEN_FILE_FORMAT /= "QUANTUM" ) then
             m=j
             do k=1,size(localizationOfCenters%values,dim=1)

                wasPress=.false.
                do i=1,j
                   if(  abs( auxMatrix%values(i,1) - localizationOfCenters%values(k,1)) < 1.0D-9 .and. &
                        abs( auxMatrix%values(i,2) - localizationOfCenters%values(k,2)) < 1.0D-9 .and. &
                        abs( auxMatrix%values(i,3) - localizationOfCenters%values(k,3)) < 1.0D-9  ) then
                      wasPress=.true.
                   end if
                end do

                if( .not.wasPress) then
                   m=m+1

                   totalNumberOfParticles = totalNumberOfParticles + 1
                   origin=localizationOfCenters%values(k,:)
                   if ( CONTROL_instance%UNITS=="ANGS") origin = origin * ANGSTROM
                   symbol=labels(k)
                   if(scan(symbol,"_") /=0) symbol=symbol(1:scan(symbol,"_")-1)

                   write (10,"(A,I8,I8,F15.8,F15.8,F15.8,I8)") trim(symbol), m,int(abs(charges(k))), origin(1), origin(2), origin(3)

                end if

             end do
          end if
          !          print *, "totalNumberOfParticles ", totalNumberOfParticles
          !         print *, "particles for specie", size(MolecularSystem_instance%species(l)%particles)

          write(10,"(A)") "[GTO]"
          j=0
          do i=1,size(MolecularSystem_instance%species(l)%particles)

             !              if ( trim(MolecularSystem_instance%species(l)%particles(i)%symbol) == trim(auxString) ) then
             j=j+1

             write(10,"(I3,I2)") j,0
             call BasisSet_showInSimpleForm( MolecularSystem_instance%species(l)%particles(i)%basis,&
                  trim(MolecularSystem_instance%species(l)%particles(i)%nickname),10 )
             write(10,*) ""

          end do

          if ( totalNumberOfParticles > size(MolecularSystem_instance%species(l)%particles) ) then
             if ( CONTROL_instance%MOLDEN_FILE_FORMAT == "MIXED" ) then
                do n = 1, ( totalNumberOfParticles - size(MolecularSystem_instance%species(l)%particles) )
                   write(10,"(I3,I2)") j+n,0
                   write(10,"(A,I1,F5.2)") " s  ",1,1.00
                   write(10,"(ES19.10,ES19.10)") 1.00,1.00
                   write(10,*) ""
                end do
             end if
          end if
          !              end if
          write(10,*) ""

          write(10,"(A)") "[MO]"

          numberOfContractions = MolecularSystem_getTotalNumberOfContractions(l)
          occupationTotal=MolecularSystem_getOcupationNumber(l)

          call MolecularSystem_changeOrbitalOrder(coefficientsOfcombination(l,state),l,"LOWDIN","MOLDEN")

          do j=1,size(energyOfMolecularOrbital(l,state)%values)
             write (10,"(A5,ES15.5)") "Ene= ",energyOfMolecularOrbital(l,state)%values(j)

             write (10,"(A11)") "Spin= Alpha"

             write (10,"(A6,F15.10)") "Occup= ",fractionalOccupations(l,state)%values(j)

             i = 0
             do k=1,size(coefficientsOfCombination(l,state)%values,dim=1)
                i = i + 1
                write(10,"(I4,A2,F15.8)") k,"  ", coefficientsOfCombination(l,state)%values(k,j)
             end do

             if ( totalNumberOfParticles > size(MolecularSystem_instance%species(l)%particles) ) then
                if ( CONTROL_instance%MOLDEN_FILE_FORMAT == "MIXED" ) then
                   do n = 1, ( totalNumberOfParticles - size(MolecularSystem_instance%species(l)%particles) )
                      write(10,"(I4,A2,ES15.8)") i+n,"  ", 0.0
                   end do
                end if
             end if

          end do


          close(10)
       end do
    end do

    ! call Matrix_destructor( localizationOfCenters )
    ! call Matrix_destructor( auxMatrix )
    ! deallocate(labels)

    !     end if

  end subroutine OutputBuilder_writeMoldenFile


!!!!!!!!!!!!!!!!!!!!!!!!!GAMESS .VEC FILE LAURA 

  subroutine OutputBuilder_VecGamessFile(this)
    implicit none
    type(OutputBuilder) :: this

    integer :: i
    integer :: j
    integer :: l
    integer :: m
    integer :: specieID
    character(10) :: auxString
    type(Matrix) :: coefficientsOfcombination
    integer :: wfnUnit
    character(100) :: wfnFile
    integer :: numberOfContractions
    character(50) :: arguments(2)

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction                                                                                     
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    do l=1,MolecularSystem_getNumberOfQuantumSpecies()

       auxString=MolecularSystem_getSymbolOfSpecies( l )

       this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".vec"

       open(29,file=this%fileName(l),status='replace',action='write')

       specieID = int( MolecularSystem_getSpecieID(nameOfSpecie = trim(auxString)) )
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
       arguments(2) = MolecularSystem_getNameOfSpecies(specieID)

       arguments(1) = "COEFFICIENTS"
       coefficientsOfcombination = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

       !! Build a vector of labels of contractions
       call MolecularSystem_changeOrbitalOrder(coefficientsOfcombination,l,"LOWDIN","GAMESS")

       do i =1, numberOfContractions
          j =1
          !if (mod(numberOfContractions,2)) then
          if (mod(numberOfContractions,2) == 1 ) then
!!!Se activa cuando el numberOfContractions es impar                                  
             write (29,"(I2,I3)",advance='no') mod(i,100),j
             do m=1,numberOfContractions

                if (mod(m,5)==0) then
                   write (29,"(ES15.8)") coefficientsOfCombination%values(m,i)
                   j=j+1
                   if (m<numberOfContractions) then
                      write (29,"(I2,I3)",advance='no') mod(i,100),j
                   end if
                else
                   write (29,"(ES15.8)",advance='no') coefficientsOfCombination%values(m,i)
                end if
             end do
             !write (29, "(A)", advance='yes')" "
             if (m<numberOfContractions) then
                write (29,"(A)", advance='no')" "
                !write (29,"(A)")" "
             end if

          else
             ! !!!Se activa cuando el numberOfContractions es par                                  
             write (29,"(I2,I3)",advance='no') mod(i,100),j
             do m=1,numberOfContractions

                if (mod(m,5)==0) then
                   write (29,"(ES15.8)") coefficientsOfCombination%values(m,i)
                   j=j+1
                   if (m<numberOfContractions) then
                      write (29,"(I2,I3)",advance='no') mod(i,100),j
                   end if
                else
                   write (29,"(ES15.8)",advance='no') coefficientsOfCombination%values(m,i)
                end if
             end do
             !write (29, "(A)", advance='yes')" "
             if (m<numberOfContractions) then
                write (29,"(A)", advance='no')" "
                !write (29,"(A)")" "
             end if

          end if

          if (.not. mod(m-1,5)==0)write (29,"(A)", advance='yes')" "
       end do

       close(29)
    end do


    !        call Matrix_destructor( localizationOfCenters )
    !        call Matrix_destructor( auxMatrix )
    !        deallocate(labels)

    !     end if

  end subroutine OutputBuilder_VecGamessFile

!!!!!!!!!!END GAMESS .VEC FILE LAURA

  subroutine OutputBuilder_casinoFile(this)
    implicit none
    type(OutputBuilder) :: this

    integer :: i
    integer :: j
    integer :: l
    integer :: g, h, m
    integer :: specieID
    type(Matrix) :: coefficientsOfcombination
    real(8), allocatable :: superMatrix(:,:)
    integer :: wfnUnit
    character(50) :: wfnFile
    integer :: numberOfContractions, superSize
    integer :: numberOfContractionsA, numberOfContractionsB
    integer :: numberOfShellsA, numberOfShellsB, totalShells
    character(50) :: arguments(20)
    integer :: i0, j0, maxl, shellCode
    real(8) :: puntualInteractionEnergy

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction                                                                                     
    open(unit = wfnUnit, file = trim(wfnFile), status = "old", form = "unformatted")


    this%fileName = trim(CONTROL_instance%INPUT_FILE)//"casino"
    open(29,file=this%fileName(1),status='replace',action='write')

    select case ( MolecularSystem_getNumberOfQuantumSpecies() ) 
    case (1) 
       numberOfContractionsA = MolecularSystem_getTotalNumberOfContractions(1)
       numberOfContractionsB = 0
       numberOfShellsA = MolecularSystem_getNumberOfContractions(1)
       numberOfShellsB = 0

    case (2) 
       numberOfContractionsA = MolecularSystem_getTotalNumberOfContractions(1)
       numberOfContractionsB = MolecularSystem_getTotalNumberOfContractions(2)
       numberOfShellsA = MolecularSystem_getNumberOfContractions(1)
       numberOfShellsB = MolecularSystem_getNumberOfContractions(2)
    case default
       call Exception_stopError("The maximum number of quantum species cannot be greater than two", "OutputBuilder_casinoFile" )
    end select

    totalShells = numberOfShellsA + numberOfShellsB

    superSize = 0
    maxl = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       specieID = l 
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
       superSize = superSize + numberOfContractions
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             maxl = max( maxl, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%angularMoment)
          end do
       end do
    end do


    !! Basic info
    write (29,*) " Title"
    write (29,*) ""
    write (29,*) "BASIC_INFO"
    write (29,*) "---------"
    write (29,*) "Generated by:"
    write (29,*) "LOWDIN"
    write (29,*) "Method:"
    write (29,*) CONTROL_instance%METHOD
    write (29,*) "DFT Functional:"
    write (29,*) CONTROL_instance%ELECTRON_EXCHANGE_CORRELATION_FUNCTIONAL 
    write (29,*) "Periodicity:"
    write (29,*) "0"
    write (29,*) "Spin unrestricted:"
    if ( CONTROL_instance%IS_OPEN_SHELL ) write (29,*) ".true."
    if ( .not. CONTROL_instance%IS_OPEN_SHELL ) write (29,*) ".false."
    write (29,*) "nuclear-nuclear repulsion energy (au/atom):"
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=puntualInteractionEnergy, arguments=["PUNTUALINTERACTIONENERGY"])
    write (29,*) puntualInteractionEnergy
    write (29,*) "Number of electrons per primitive cell:" !! ?
    write (29,*) "2"
    write (29,*) ""

    !! Geometry
    write (29,*) "GEOMETRY"
    write (29,*) "---------"
    write (29,*) "Number of atoms:" !! centers?
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       m = m + size(MolecularSystem_instance%species(l)%particles)
    end do
    write (29,"(T4,I4)") m
    write (29,*) "Atomic positions (au):" !! centers?
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          m = m + 1 
          write (29, "(3ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%origin(1:3)
       end do
    end do
    write (29,*) "Atomic numbers for each atom:" 
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          m = m + 1 
          if (mod(m,8)==0) then
             write (29,"(I10)") int(MolecularSystem_instance%species(l)%particles(g)%charge)
          else
             write (29,"(I10)",advance="no") int(MolecularSystem_instance%species(l)%particles(g)%charge)
          end if
       end do
    end do
    if (.not. mod(m,8)==0)  write (29,"(A)", advance='yes') " "
    !write (29,*) "_ii_ _ii_"
    !write (29,"(2I10)") 1,0
    write (29,*) "Valence charges for each atom:" !! what?
    !write (29,*) " 1.0000000000000E+00 0.0000000000000E+00"
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          m = m + 1 
          if (mod(m,4)==0) then
             write (29,"(ES20.13)") MolecularSystem_instance%species(l)%particles(g)%charge
          else
             write (29,"(ES20.13)",advance="no") MolecularSystem_instance%species(l)%particles(g)%charge
          end if
       end do
    end do
    if (.not. mod(m,8)==0)  write (29,"(A)", advance='yes') " "
    write (29,*) ""
    !! Basis set

    write (29,*) "BASIS SET"
    write (29,*) "---------"
    write (29,*) "Number of Gaussian centres"
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       m = m + size(MolecularSystem_instance%species(l)%particles)
    end do
    write (29,"(T4,I4)") m
    write (29,*) "Number of shells per primitive cell" !! total?
    write (29,"(T4,I4)") totalShells
    write (29,*) "Number of basis functions ('AO') per primitive cell"
    write (29,"(T4,I4)") superSize
    write (29,*) "Number of Gaussian primitives per primitive cell"
    write (29,"(T4,I4)") totalShells
    write (29,*) "Highest shell angular momentum (s/p/d/f... 1/2/3/4...)"
    write (29,"(T4,I4)") maxl+1
    write (29,*) "Code for shell types (s/sp/p/d/f... 1/2/3/4/5...) "
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             m = m + 1 
             shellCode = 0
             if ( MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%angularMoment == 0 ) shellCode = 1
             if ( MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%angularMoment > 0 ) shellCode = 2

             if (mod(m,8)==0) then
                write (29,"(I10)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%angularMoment + shellCode
             else
                write (29,"(I10)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%angularMoment + shellCode
             end if
          end do
       end do
    end do
    if (.not. mod(m,8)==0)  write (29,"(A)", advance='yes') " "

    write (29,*) "Number of primitive Gaussians in each shell"
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             m = m + 1 
             if (mod(m,8)==0) then
                write (29,"(I10)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%length
             else
                write (29,"(I10)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%length
             end if
          end do
       end do
    end do
    if (.not. mod(m,8)==0)  write (29,"(A)", advance='yes') " "

    write (29,*) "Sequence number of first shell on each centre"
    write (29,"(3I10)") 1,numberOfShellsA, numberOfShellsA+numberOfShellsB+1
    write (29,*) "Exponents of Gaussian primitives"
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             do i = 1, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%length
                m = m + 1 
                if (mod(m,4)==0) then
                   write (29,"(ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%orbitalExponents(i)
                else
                   write (29,"(ES20.13)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%orbitalExponents(i)
                end if
             end do
          end do
       end do
    end do
    if (.not. mod(m,4)==0)  write (29,"(A)", advance='yes') " "
    write (29,*) "Normalised contraction coefficients" !! check this...
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             do i = 1, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%length
                !            do j = 1, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%numCartesianOrbital
                !              m = m + 1 
                !              if (mod(m,4)==0) then
                !                write (29,"(ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%primNormalization(i,j)
                !              else
                !                write (29,"(ES20.13)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%primNormalization(i,j)
                !              end if
                !            end do

                !            do j = 1, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%numCartesianOrbital
                m = m + 1 
                if (mod(m,4)==0) then
                   write (29,"(ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%primNormalization(i,1)
                else
                   write (29,"(ES20.13)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%primNormalization(i,1)
                end if
                !            end do

             end do
          end do
       end do
    end do
    if (.not. mod(m,4)==0)  write (29,"(A)", advance='yes') " "
    !    m = 0
    !    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
    !      do g = 1,  size(MolecularSystem_instance%species(l)%particles)
    !        do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
    !          do i = 1, MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%length
    !            m = m + 1 
    !            if (mod(m,4)==0) then
    !              write (29,"(ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%contNormalization(i)
    !            else
    !              write (29,"(ES20.13)",advance="no") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%contNormalization(i)
    !            end if
    !          end do
    !        end do
    !      end do 
    !    end do 
    write (29,*) "Position of each shell (au)"
    m = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()
       do g = 1,  size(MolecularSystem_instance%species(l)%particles)
          do h = 1, size(MolecularSystem_instance%species(l)%particles(g)%basis%contraction)
             m = m + 1 
             write (29,"(3ES20.13)") MolecularSystem_instance%species(l)%particles(g)%basis%contraction(h)%origin(1:3)
          end do
       end do
    end do
    write (29,"(A)", advance='yes')" "


    write (29,*) "MULTIDETERMINANT INFORMATION"
    write (29,*) "----------------------------"
    write (29,"(A2)") "GS"
    write (29,*) ""

    !! coefficients
    write (29,*) "EIGENVECTOR COEFFICIENTS"
    write (29,*) "------------------------"



    !! Save the MO coefficients in a supermatrix from for all quantum species (2...)
    if ( allocated (superMatrix) ) deallocate (superMatrix)
    allocate (superMatrix(superSize,superSize)) 
    superMatrix = 0

    i0 = 0
    j0 = 0
    do l = 1,MolecularSystem_getNumberOfQuantumSpecies()

       specieID = l 
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
       arguments(2) = MolecularSystem_getNameOfSpecies(specieID)
       arguments(1) = "COEFFICIENTS"
       coefficientsOfcombination = Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

       do i =1, numberOfContractions
          do j =1, numberOfContractions
             superMatrix(i+i0,j+j0) = coefficientsOfCombination%values(i,j)
          end do
       end do
       !! starting positron for the next species
       i0 = i-1
       j0 = j-1
       print *, "i0 j0", i0, j0

    end do

    do i =1, superSize
       j =1
       !if (mod(numberOfContractions,2)) then
       if (mod(superSize,2) == 1 ) then
!!!Se activa cuando el numberOfContractions es impar                                  
          do m=1,superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m <= superSize) then
             write (29,"(A)", advance='no')" "
             !write (29,"(A)")" "
          end if

       else
          !!Se activa cuando el numberOfContractions es par                                  
          do m=1, superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m <= superSize) then
             write (29,"(A)", advance='no') " "
             !write (29,"(A)")" "
          end if

       end if

       if (.not. mod(m-1,4)==0)  write (29,"(A)", advance='yes') " "
    end do

    !! write it twice... why?

    do i = numberOfContractionsB + 1, superSize
       j =1
       !if (mod(numberOfContractions,2)) then
       if (mod(superSize,2) == 1 ) then
!!!Se activa cuando el numberOfContractions es impar                                  
          do m=1,superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m < superSize) then
             write (29,"(A)", advance='no')" "
             !write (29,"(A)")" "
          end if

       else
          !!Se activa cuando el numberOfContractions es par                                  
          do m=1, superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m < superSize) then
             write (29,"(A)", advance='no') " "
             !write (29,"(A)")" "
          end if

       end if

       if (.not. mod(m-1,4)==0)  write (29,"(A)", advance='yes') " "
    end do

    do i = 1, numberOfContractionsA
       j =1
       !if (mod(numberOfContractions,2)) then
       if (mod(superSize,2) == 1 ) then
!!!Se activa cuando el numberOfContractions es impar                                  
          do m=1,superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m < superSize) then
             write (29,"(A)", advance='no')" "
             !write (29,"(A)")" "
          end if

       else
          !!Se activa cuando el numberOfContractions es par                                  
          do m=1, superSize

             if (mod(m,4)==0) then
                write (29,"(ES20.13)") superMatrix(m,i)
                j = j + 1
             else
                write (29,"(ES20.13)",advance='no') superMatrix(m,i)
             end if
          end do
          !write (29, "(A)", advance='yes')" "
          if (m < superSize) then
             write (29,"(A)", advance='no') " "
             !write (29,"(A)")" "
          end if

       end if

       if (.not. mod(m-1,4)==0)  write (29,"(A)", advance='yes') " "
    end do

    write (29,"(A)") ""
    close(20)
    close(29)

    call Exception_sendWarning("The order of the coefficients only works until P orbitals", "OutputBuilder_casinoFile" )

  end subroutine OutputBuilder_casinoFile

  !!Escribe los valores propios en el archivo eigenvalues.dat para que puedan ser leidos por GAMESS   Laura

  subroutine OutputBuilder_writeEigenvalues(this)
    implicit none
    type(OutputBuilder) :: this

    integer :: j
    integer :: l
    integer :: specieID
    character(10) :: auxString
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Vector) :: energyOfMolecularOrbital
    character(10),allocatable :: labels(:)
    integer :: wfnUnit
    character(100) :: wfnFile
    integer :: numberOfContractions
    character(50) :: arguments(2)
    integer :: totalNumberOfParticles

    ! auxString="speciesName"

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    ! auxString=MolecularSystem_getNameOfSpecies( 1 )
    ! this%fileName=trim(CONTROL_instance%INPUT_FILE)//".eigen"
    ! open(129,file=this%fileName,status='replace',action='write')
    ! close(129)

    localizationOfCenters=ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    auxMatrix=localizationOfCenters
    allocate( labels( size(auxMatrix%values,dim=1) ) )
    allocate( charges( size(auxMatrix%values,dim=1) ) )
    labels=ParticleManager_getLabelsOfCentersOfOptimization()
    charges=ParticleManager_getChargesOfCentersOfOptimization()

    !! Open file for wavefunction                                                                                     
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")


    do l=1,MolecularSystem_getNumberOfQuantumSpecies()

       totalNumberOfParticles = 0

       auxString=MolecularSystem_getSymbolOfSpecies( l )
       specieID = MolecularSystem_getSpecieID(auxString)
       this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".eigen"
       open(129,file=this%fileName(l),status='replace',action='write')

       specieID = int( MolecularSystem_getSpecieID(nameOfSpecie = trim(auxString)) )
       numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
       arguments(2) = MolecularSystem_getNameOfSpecies(specieID)

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = numberOfContractions, &
            unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
            output = energyOfMolecularOrbital )

       do j=1,size(energyOfMolecularOrbital%values)
          write (129,"(F15.12)") energyOfMolecularOrbital%values(j)
       end do
       close(129)
    end do

    call Matrix_destructor( localizationOfCenters )
    call Matrix_destructor( auxMatrix )
    deallocate(labels)


  end  subroutine OutputBuilder_writeEigenvalues


  subroutine OutputBuilder_writeFchkFile(this)
    implicit none
    type(OutputBuilder) :: this

    integer :: h,i,j,k,l
    integer :: numberOfSpecies
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Vector) :: energyOfMolecularOrbital
    type(Matrix) :: coefficientsOfcombination
    character(10),allocatable :: labels(:)

    integer :: wfnUnit
    character(100) :: wfnFile
    integer :: numberOfPrimitives
    integer :: numberOfContractions
    integer :: numberOfShells
    integer :: numberOfAtoms
    integer :: numberOfAtomShells
    integer :: numberOfShellPrimitives
    real(8) :: particlesPerOrbital
    type(matrix) :: densityMatrix
    real(8) :: densityElement
    character(50) :: nameOfSpecies, symbolOfSpecies
    character(40) :: header
    character(50) :: arguments(2)


    localizationOfCenters=ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    auxMatrix=localizationOfCenters
    allocate( labels( size(auxMatrix%values,dim=1) ) )
    allocate( charges( size(auxMatrix%values,dim=1) ) )
    labels=ParticleManager_getLabelsOfCentersOfOptimization()
    charges=ParticleManager_getChargesOfCentersOfOptimization()
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()


    !! Open file for wavefunction                                                                                     
    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    ! do state=1,numberOfStates
    do l=1,numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecies(l)
       symbolOfSpecies=MolecularSystem_getSymbolOfSpecies(l)
       particlesPerOrbital=MolecularSystem_getLambda(l)
       ! if (state .eq. 1) then
       !    auxString=nameOfSpecies
       ! else
       !    write(auxString, "(I8)")  state
       !    auxString=trim(nameOfSpecies)//"-"//trim( adjustl(auxString))
       ! end if

       ! this%fileName=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".fchk"
       this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)//".fchk"

       numberOfAtoms=size(MolecularSystem_instance%species(l)%particles)
       numberOfShells=MolecularSystem_getNumberOfContractions(l)
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(l)

       arguments(2) = MolecularSystem_getNameOfSpecies(l)
       arguments(1) = "COEFFICIENTS"
       coefficientsOfcombination = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
            columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = numberOfContractions, &
            unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
            output = energyOfMolecularOrbital )

       call Matrix_constructor(densityMatrix, int(numberOfContractions,8), int(numberOfContractions,8), 0.0_8)

       open(10,file=this%fileName(l),status='replace',action='write')

       write(10,"(A)") this%fileName(l)

       write(10,"(A10,A30,A30)") "SP", trim(CONTROL_instance%METHOD), trim(MolecularSystem_instance%species(l)%particles(1)%basis%name)

       if(particlesPerOrbital .ne. 2.0) then
          header="Number of electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", MolecularSystem_getOcupationNumber(l)
          header="Number of alpha electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", MolecularSystem_getOcupationNumber(l)
          header="Number of beta electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", 0
       else
          header="Number of electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", MolecularSystem_getOcupationNumber(l)*2
          header="Number of alpha electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", MolecularSystem_getOcupationNumber(l)
          header="Number of beta electrons"
          write(10,"(A40,3X,A1,5X,I12)") header, "I", MolecularSystem_getOcupationNumber(l)
       end if

       header="Number of atoms"
       write(10,"(A40,3X,A1,5X,I12)") header, "I", numberOfAtoms

       header="Number of basis functions"
       write(10,"(A40,3X,A1,5X,I12)") header, "I", numberOfContractions
!!! Falta restar los orbitales removidos 
       header="Number of independent functions"
       write(10,"(A40,3X,A1,5X,I12)") header, "I", numberOfContractions

       header="Atomic numbers"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "I", "N=", numberOfAtoms
       do j=1, numberOfAtoms
          if(mod(j,6).eq.0 .or. j.eq.numberOfAtoms ) then
             ! write(10,"(I12)") MolecularSystem_instance%species(l)%particles(j)%internalSize
             write(10,"(I12)") int(abs(MolecularSystem_instance%species(l)%particles(j)%totalCharge))
          else
             ! write(10,"(I12)", advance="no") MolecularSystem_instance%species(l)%particles(j)%internalSize
             write(10,"(I12)", advance="no") int(abs(MolecularSystem_instance%species(l)%particles(j)%totalCharge))
          end if
       end do

       header="Current cartesian coordinates"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "R", "N=", numberOfAtoms*3
       k=0
       do j=1, numberOfAtoms
          origin = MolecularSystem_instance%species(l)%particles(j)%origin 
          do i=1,3
             k=k+1
             if(mod(k,5).eq.0 .or. (j.eq.numberOfAtoms .and. i.eq.3)) then
                write(10,"(ES16.8)") origin(i)
             else
                write(10,"(ES16.8)", advance="no") origin(i)
             end if
          end do
       end do

       header="Shell types"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "I", "N=", numberOfShells
       k=0
       do j=1, numberOfAtoms
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells 
             k=k+1
             if(mod(k,6).eq.0 .or. (j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells)) then
                write(10,"(6I12)") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%angularMoment
             else
                write(10,"(6I12)", advance="no") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%angularMoment
             end if
          end do
       end do

       header="Shell to atom map"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "I", "N=", numberOfShells
       k=0
       do j=1, numberOfAtoms
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells 
             k=k+1
             if(mod(k,6).eq.0 .or.(j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells)) then
                write(10,"(6I12)") j
             else
                write(10,"(6I12)", advance="no") j
             end if
          end do
       end do

       header="Coordinates of each shell"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "R", "N=", numberOfShells*3
       k=0
       do j=1, numberOfAtoms
          origin = MolecularSystem_instance%species(l)%particles(j)%origin 
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells 
             do h=1, 3
                k=k+1
                if(mod(k,5).eq.0 .or.(j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells .and. h.eq.3)) then
                   write(10,"(ES16.8)") origin(h)
                else
                   write(10,"(ES16.8)", advance="no") origin(h)
                end if
             end do
          end do
       end do

       numberOfPrimitives=0
       header="Number of primitives per shell"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "I", "N=", numberOfShells
       k=0
       do j=1, numberOfAtoms
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells
             numberOfShellPrimitives=MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%length
             numberOfPrimitives=numberOfPrimitives+numberOfShellPrimitives
             k=k+1
             if(mod(k,6).eq.0 .or.(j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells)) then
                write(10,"(6I12)") numberOfShellPrimitives
             else
                write(10,"(6I12)", advance="no") numberOfShellPrimitives
             end if
          end do
       end do

       header="Primitive exponents"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "R", "N=", numberOfPrimitives
       k=0
       do j=1, numberOfAtoms
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells 
             numberOfShellPrimitives=MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%length
             do h=1, numberOfShellPrimitives
                k=k+1
                if(mod(k,5).eq.0 .or.(j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells .and. h.eq. numberOfShellPrimitives)) then
                   write(10,"(ES16.8)") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%orbitalExponents(h)
                else
                   write(10,"(ES16.8)", advance="no") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%orbitalExponents(h)
                end if
             end do
          end do
       end do

       header="Contraction coefficients"
       write(10,"(A40,3X,A1,3X,A2,I12)") header, "R", "N=", numberOfPrimitives
       k=0
       do j=1, numberOfAtoms
          numberOfAtomShells=size(MolecularSystem_instance%species(l)%particles(j)%basis%contraction(:))
          do i=1, numberOfAtomShells 
             numberOfShellPrimitives=MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%length
             do h=1, numberOfShellPrimitives
                k=k+1
                if(mod(k,5).eq.0 .or.(j.eq.numberOfAtoms .and. i.eq.numberOfAtomShells .and. h.eq. numberOfShellPrimitives)) then
                   write(10,"(ES16.8)") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%contractionCoefficients(h)
                else
                   write(10,"(ES16.8)", advance="no") MolecularSystem_instance%species(l)%particles(j)%basis%contraction(i)%contractionCoefficients(h)
                end if
             end do
          end do
       end do

       header="Alpha Orbital Energies"
       write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions
       k=0
       do j=1, numberOfContractions
          k=k+1
          if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions)) then
             write(10,"(ES16.8)") energyOfMolecularOrbital%values(j)
          else
             write(10,"(ES16.8)", advance="no") energyOfMolecularOrbital%values(j)
          end if
       end do

       if(particlesPerOrbital .ne. 2.0) then
          header="Beta Orbital Energies"
          write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions
          k=0
          do j=1, numberOfContractions
             k=k+1
             if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions)) then
                write(10,"(ES16.8)") 0.0
             else
                write(10,"(ES16.8)", advance="no") 0.0
             end if
          end do
       end if


       call MolecularSystem_changeOrbitalOrder( coefficientsOfCombination, l, "LOWDIN", "FCHK" )

       header="Alpha MO coefficients"
       write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions**2
       k=0

       do j=1, numberOfContractions
          do i=1, numberOfContractions
             k=k+1
             if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions .and. i.eq.numberOfContractions)) then
                write(10,"(ES16.8)") coefficientsOfCombination%values(i,j)
             else
                write(10,"(ES16.8)", advance="no") coefficientsOfCombination%values(i,j)
             end if
          end do
       end do

       if(particlesPerOrbital .ne. 2.0) then
          header="Beta MO coefficients"
          write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions**2
          k=0

          do j=1, numberOfContractions
             do i=1, numberOfContractions
                k=k+1
                if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions .and. i.eq.numberOfContractions)) then
                   write(10,"(ES16.8)") 0.0
                else
                   write(10,"(ES16.8)", advance="no") 0.0
                end if
             end do
          end do
       end if

       !Build density matrix with the new order
       do i=1, numberOfContractions
          do j=1, numberOfContractions
             do k=1, MolecularSystem_getOcupationNumber(l)
                densityMatrix%values(i,j) =  &
                     densityMatrix%values( i,j ) + &
                     coefficientsOfCombination%values(i,k)*coefficientsOfCombination%values(j,k)*MolecularSystem_getEta(l)
             end do
          end do
       end do


       header="Total SCF Density"
       write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions*(numberOfContractions+1)/2
       k=0
       do j=1, numberOfContractions
          do i=1, j
             k=k+1
             if(i.eq.j)then
                densityElement=densityMatrix%values(i,i)
             else
                densityElement=densityMatrix%values(i,j)
             end if

             if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions .and. i.eq.j)) then
                write(10,"(ES16.8)") densityElement
             else
                write(10,"(ES16.8)", advance="no") densityElement
             end if
          end do
       end do

       if(particlesPerOrbital .ne. 2.0) then
          header="Spin SCF Density"
          write(10,"(A40,3X,A1,3X,A2,I12)") header , "R", "N=", numberOfContractions*(numberOfContractions+1)/2
          k=0
          do j=1, numberOfContractions
             do i=1, j
                k=k+1
                if(i.eq.j)then
                   densityElement=densityMatrix%values(i,i)
                else
                   densityElement=densityMatrix%values(i,j)
                end if

                if(mod(k,5).eq.0 .or.(j.eq.numberOfContractions .and. i.eq.j)) then
                   write(10,"(ES16.8)") densityElement
                else
                   write(10,"(ES16.8)", advance="no") densityElement
                end if
             end do
          end do
       end if
    end do
    ! end do
  end subroutine OutputBuilder_writeFchkFile



  !**
  ! @brief Call the molden2aim library to generate the wfn, wfx or NBO47 files from a molden file.
  !**

  subroutine OutputBuilder_generateAIMFiles (this)
    implicit none
    type(OutputBuilder) :: this
    character(50) :: auxString
    character(100) :: initialSettingsFile
    character(100) :: moldenFileName
    integer :: l
    character(100) :: wfnFile
    character(2) :: wfnStatus, wfxStatus, nboStatus
    character(10) :: extension
    integer :: wfnUnit
    real(8) :: totalEnergy, virial

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    initialSettingsFile = "m2a.ini"

    select case (this%type) 
    case ( "WFNFILE" )
       wfnStatus="1"
       nboStatus="-1"
       wfxStatus="-1"
       extension=".wfn"
    case ( "NBO47FILE" )
       wfnStatus="-1"
       nboStatus="1"
       wfxStatus="-1"
       extension=".47"
    case ( "WFXFILE" ) 
       wfnStatus="-1"
       nboStatus="-1"
       wfxStatus="1"
       extension=".wfx"
    case ( "EXTENDEDWFNFILE" )
       wfnStatus="1"
       nboStatus="-1"
       wfxStatus="-1"
       extension=".wfn"
    end select

    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=totalEnergy, arguments=["TOTALENERGY"])
    call Vector_getFromFile(unit=wfnUnit, binary=.true., value=virial, arguments=["VIRIAL"])
    close(wfnUnit)
    open (35,file=initialSettingsFile,status='unknown',action='write')
    write(35,"(A)") " ######################################################################## "
    Write(35,"(A)") " #  In the following 6 parameters "
    Write(35,"(A)") " #     >0:  always performs the operation without asking the user "
    Write(35,"(A)") " #     =0:  asks the user whether to perform the operation "
    Write(35,"(A)") " #     <0:  always neglect the operation without asking the user "
    Write(35,"(A)") " molden=1           ! Generating a standard Molden file in Cart. function "
    Write(35,"(A)") " wfn="//wfnStatus//"              ! Generating a WFN file "
    Write(35,"(A)") " wfncheck=-1         ! Checking normalization for WFN "
    Write(35,"(A)") " wfx="//wfxStatus//"             ! Generating a WFX file (not implemented) "
    Write(35,"(A)") " wfxcheck=-1        ! Checking normalization for WFX (not implemented) "
    Write(35,"(A)") " nbo="//nboStatus//"              ! Generating a NBO .47 file "
    Write(35,"(A)") " nbocheck=-1         ! Checking normalization for NBO's .47 "
    Write(35,"(A)") " ######################################################################## "
    Write(35,"(A)") " #  Which quantum chemistry program is used to generate the MOLDEN file? "
    Write(35,"(A)") " #  1: ORCA "
    Write(35,"(A)") " #  5: ACES2 "
    Write(35,"(A)") " #  0: other programs "
    Write(35,"(A)") " # "
    Write(35,"(A)") " #  If non-zero value is given "
    Write(35,"(A)") " # "
    Write(35,"(A)") " program=0 "
    Write(35,"(A)") " ######################################################################## "
    Write(35,"(A)") " #  Which orbirals will be printed in the WFN/WFX file? "
    Write(35,"(A)") " # =0: print only the orbitals with occ. number > 5.0d-8 "
    Write(35,"(A)") " # <0: print only the orbitals with occ. number > 0.1 (debug only) "
    Write(35,"(A)") " # >0: print all the orbitals "
    Write(35,"(A)") " iallmo=1 "
    Write(35,"(A)") " ######################################################################## "
    Write(35,"(A)") " #  Print supporting information or not "
    Write(35,"(A)") " # =0: print "
    Write(35,"(A)") " nosupp=-1 "
    Write(35,"(A)") " ######################################################################## "
    Write(35,"(A)") " #  The following parameters are used only for debugging. "
    Write(35,"(A)") " clear=1            ! delete temporary files (1) or not (0) "
    Write(35,"(A)") " ######################################################################## "
    close(35)

    do l=1,MolecularSystem_getNumberOfQuantumSpecies()
       auxString=MolecularSystem_getSymbolOfSpecies( l )
       moldenFileName=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".molden"
       call Molden2AIM(moldenFileName, totalEnergy, virial)
       !! Just for printing information 
       this%fileName(l) = trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//extension//" and .molden"
    end do


  end subroutine OutputBuilder_generateAIMFiles

  !! For future implementation

  subroutine OutputBuilder_generateExtendedWfnFile (this)
    implicit none
    type(OutputBuilder) :: this
    integer :: l
    character(50) :: auxString

    do l=1,MolecularSystem_getNumberOfQuantumSpecies()
       auxString=MolecularSystem_getNameOfSpecies( l )
    end do

  end subroutine OutputBuilder_generateExtendedWfnFile

  subroutine OutputBuilder_getCube(this )
    implicit none
    type(OutputBuilder) :: this
    character(50) :: outputID, auxID

    integer :: l, i, j, k, n, w, natom
    integer :: speciesID, state
    integer :: numberOfSteps
    real(8) :: step
    real(8) :: lowerLimit(3)
    Type(Vector) :: val
    Type(Matrix) :: coordinate

    integer :: numberOfSpecies, numberOfStates, orbital
    type(matrix), allocatable :: densityMatrices(:,:), coefficientsMatrices(:,:)

    character(100) :: nameOfSpecies, symbolOfSpecies

    !Writes Gaussian Cube 

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    auxID=String_convertIntegerToString(this%auxID)
    numberOfStates=1
    if (this%wavefunctionType .eq. "CI") then
       numberOfStates=CONTROL_instance%CI_STATES_TO_PRINT
       state=this%state
    else
       state=1
    end if
    
    if (this%type=="DENSITYCUBE") then
       allocate(densityMatrices(numberOfSpecies,numberOfStates))
       call CalculateWaveFunction_loadDensityMatrices ( numberOfSpecies, numberOfStates, this%wavefunctionType, densityMatrices )
       
    else if(this%type=="ORBITALCUBE") then
       allocate(coefficientsMatrices(numberOfSpecies,numberOfStates))
       call CalculateWaveFunction_loadCoefficientsMatrices ( numberOfSpecies, numberOfStates, this%wavefunctionType, coefficientsMatrices)
    end if
    
    l=0
    do speciesID=1, numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       symbolOfSpecies=MolecularSystem_getSymbolOfSpecies(speciesID)
       if(trim(this%species) .eq. trim(symbolOfSpecies) .or. trim(this%species) .eq. "ALL" ) then
          l=l+1   
          outputID=String_convertIntegerToString(this%outputID)

          this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)
          if( state .gt. 1) this%fileName(l)=trim(this%fileName(l))//"-"//trim(String_convertIntegerToString(state))
          if( this%auxID .gt. 1) this%fileName(l)=trim(this%fileName(l))//"."//trim(auxID)
          if( this%type=="DENSITYCUBE") this%fileName(l)=trim(this%fileName(l))//".dens.cub"
          if( this%type=="ORBITALCUBE") then
             orbital=this%orbital
             if(orbital.eq.0) orbital=MolecularSystem_getOcupationNumber(speciesID)
             this%fileName(l)=trim(this%fileName(l))//".orb"//trim(String_convertIntegerToString(orbital))//".cub"
          end if

          open(10,file=this%fileName(l),status='replace',action='write')

          lowerLimit(:)=this%point1%values(:)-this%cubeSize/2
          numberOfSteps=this%pointsPerDim(1) !for now, only cubic cubes
          step= this%cubeSize/numberOfSteps

          natom=MolecularSystem_instance%numberOfPointCharges

          write (10,"(A)") "Gaussian Cube generated with Lowdin Software"
          write (10,"(A)") this%fileName(l)
          if(natom .gt. 0) then
             write (10,"(I8,F20.8,F20.8,F20.8,I8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3), 1
          else
             write (10,"(I8,F20.8,F20.8,F20.8,I8)") 1, lowerLimit(1), lowerLimit(2), lowerLimit(3), 1
          end if
          write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, step, 0.0, 0.0
          write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, 0.0, step, 0.0
          write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, 0.0, 0.0, step

          if(natom .gt. 0) then
             do n = 1, MolecularSystem_instance%numberOfPointCharges
                write (10, "(I8,F20.8,F20.8,F20.8,F20.8)") &
                     int(MolecularSystem_instance%pointCharges(n)%charge), 0.0, MolecularSystem_instance%pointCharges(n)%origin(1:3)
             end do
          else
             write (10, "(I8,I8,F20.8,F20.8,F20.8)") &
                  1, 0, this%point1%values
          end if

          do i=1,numberOfSteps
             do j=1, numberOfSteps
                call Matrix_constructor(coordinate,int(numberOfSteps,8),int(3,8),0.0_8)
                do k=1, numberOfSteps
                   coordinate%values(k,1)=lowerLimit(1)+(i-1)*step
                   coordinate%values(k,2)=lowerLimit(2)+(j-1)*step
                   coordinate%values(k,3)=lowerLimit(3)+(k-1)*step
                end do
                if( this%type=="DENSITYCUBE") then
                   call CalculateWaveFunction_getDensityAt( speciesID, coordinate, densityMatrices(speciesID,state), val )
                else if( this%type=="ORBITALCUBE") then
                   call CalculateWaveFunction_getOrbitalValueAt( speciesID, orbital, coordinate, coefficientsMatrices(speciesID,state), val )
                end if
                write(10,*) ( val%values(w) , w=1,numberOfSteps )
                write(10,*) ( "" )
             end do
          end do
          close(10)
       end if
    end do
  end subroutine OutputBuilder_getCube

  subroutine OutputBuilder_getPlot(this)
    type(OutputBuilder) :: this
    character(50) :: outputID, auxID

    integer :: i,j,l,n, speciesID, orbital, state
    integer :: numberOfSteps,numberOfSteps2
    integer :: numberOfSpecies, numberOfStates
    type(matrix), allocatable :: densityMatrices(:,:), coefficientsMatrices(:,:)
    real(8) :: maxValue, minValue
    real(8) :: plotDistance1, plotDistance2, initialValue1, initialValue2
    Type(Vector) :: val
    Type(Matrix) :: coordinate

    character(100) :: nameOfSpecies, symbolOfSpecies
    character(50) :: title, x_title, y_title, z_title

    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
    outputID=String_convertIntegerToString(this%outputID)
    auxID=String_convertIntegerToString(this%auxID)
    numberOfSteps=this%pointsPerDim(1)

    !Define graph display distances, check which axes are changing
    plotDistance1=sqrt(sum((this%point2%values(:)-this%point1%values(:))**2))
    initialValue1=-0.5*plotDistance1
    x_title="distance/a.u."
    if(abs(abs(this%point2%values(1)-this%point1%values(1))-plotDistance1) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
       initialValue1=this%point1%values(1)
       x_title="X/a.u."
    else if(abs(abs(this%point2%values(2)-this%point1%values(2))-plotDistance1) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
       initialValue1=this%point1%values(2)
       x_title="Y/a.u."
    else if(abs(abs(this%point2%values(3)-this%point1%values(3))-plotDistance1) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
       initialValue1=this%point1%values(3)
       x_title="Z/a.u."
    end if

    if(this%dimensions.eq.3) then
       numberOfSteps2=this%pointsPerDim(2)
       plotDistance2=sqrt(sum((this%point3%values(:)-this%point1%values(:))**2))
       initialValue2=-0.5*plotDistance2
       y_title="distance2/a.u."
       if(abs(abs(this%point3%values(1)-this%point1%values(1))-plotDistance2) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
          initialValue2=this%point1%values(1)
          y_title="X/a.u."          
       else if(abs(abs(this%point3%values(2)-this%point1%values(2))-plotDistance2) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
          initialValue2=this%point1%values(2)
          y_title="Y/a.u."          
       else if(abs(abs(this%point3%values(3)-this%point1%values(3))-plotDistance2) .lt. CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
          initialValue2=this%point1%values(3) 
          y_title="Z/a.u."          
       end if
    end if
    if(this%axisLabel(1) .ne. "") x_title=this%axisLabel(1)//"/a.u."
    if(this%axisLabel(2) .ne. "") y_title=this%axisLabel(2)//"/a.u."

    numberOfStates=1
    if (this%wavefunctionType .eq. "CI") then
       numberOfStates=CONTROL_instance%CI_STATES_TO_PRINT
       state=this%state
    else
       state=1
    end if
    if (this%type=="DENSITYPLOT") then
       allocate(densityMatrices(numberOfSpecies,numberOfStates))
       call CalculateWaveFunction_loadDensityMatrices ( numberOfSpecies, numberOfStates, this%wavefunctionType, densityMatrices )
       if(this%dimensions.eq.2) y_title="density/a.u.^{-3}"
       if(this%dimensions.eq.3) z_title=""
       
    else if(this%type=="ORBITALPLOT") then
       allocate(coefficientsMatrices(numberOfSpecies,numberOfStates))
       call CalculateWaveFunction_loadCoefficientsMatrices ( numberOfSpecies, numberOfStates, this%wavefunctionType, coefficientsMatrices)
       if(this%dimensions.eq.2) y_title="orbital/a.u.^{-3/2}"
       if(this%dimensions.eq.3) z_title=""
       
    end if

    l=0
    do speciesID=1, numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecies(speciesID)
       symbolOfSpecies=MolecularSystem_getSymbolOfSpecies(speciesID)
       if(trim(this%species) .eq. trim(symbolOfSpecies) .or. trim(this%species) .eq. "ALL" ) then
          l=l+1   

          !Set up filename
          this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(symbolOfSpecies)
          if( state .gt. 1) this%fileName(l)=trim(this%fileName(l))//"-"//trim(String_convertIntegerToString(state))
          if( this%auxID .gt. 1) this%fileName(l)=trim(this%fileName(l))//"."//trim(auxID)
          if( this%dimensions.eq.2) this%fileName(l)=trim(this%fileName(l))//".2D"
          if( this%dimensions.eq.3) this%fileName(l)=trim(this%fileName(l))//".3D"
          if( this%type=="DENSITYPLOT") this%fileName(l)=trim(this%fileName(l))//".dens"
          if( this%type=="ORBITALPLOT") then
             orbital=this%orbital
             if(orbital.eq.0) orbital=MolecularSystem_getOcupationNumber(speciesID)
             this%fileName(l)=trim(this%fileName(l))//".orb"//trim(String_convertIntegerToString(orbital))
          end if
          open(10,file=this%fileName(l),status='replace',action='write')

          if( this%dimensions.eq.3) then 
             call Matrix_constructor(coordinate,int((numberOfSteps+1)*(numberOfSteps2+1),8),int(3,8),0.0_8)
             n=0
             do i=0,numberOfSteps
                do j=0,numberOfSteps2
                   n=n+1
                   coordinate%values(n,:)=this%point1%values(:)+i*this%step1%values(:)+j*this%step2%values(:)
                end do
             end do
             
          else if( this%dimensions.eq.2) then
             call Matrix_constructor(coordinate,int(numberOfSteps+1,8),int(3,8),0.0_8)
             do i=0,numberOfSteps
                coordinate%values(i+1,:)=this%point1%values(:)+i*this%step1%values(:)
             end do
          end if

          if( this%type=="DENSITYPLOT") call CalculateWaveFunction_getDensityAt( speciesID, coordinate, densityMatrices(speciesID,state), val )
          if( this%type=="ORBITALPLOT") call CalculateWaveFunction_getOrbitalValueAt( speciesID, orbital, coordinate, coefficientsMatrices(speciesID,state), val )

          if( this%dimensions.eq.3) then 
          
             if( this%type=="DENSITYPLOT") write (10,"(A10,A20,A20,A20)") "#","X","Y","Density"
             if( this%type=="ORBITALPLOT") write (10,"(A10,A20,A20,A20)") "#","X","Y","Orbital"
             n=0
             maxValue=0.0
             minValue=0.0
             do i=0,numberOfSteps
                write (10,*) ""
                do j=0,numberOfSteps2
                   n=n+1
                   if(abs(val%values(n))>1.0E-99_8) then
                      write (10,"(T10,F20.8,F20.8,E20.8)") initialValue1+i*Vector_norm(this%step1), initialValue2+j*Vector_norm(this%step2), val%values(n)
                   else
                      write (10,"(T10,F20.8,F20.8,E20.8)") initialValue1+i*Vector_norm(this%step1), initialValue2+j*Vector_norm(this%step2), 0.0
                   end if
                   if (val%values(n) > maxValue) maxValue = val%values(n) 
                   if (val%values(n) < minValue) minValue = val%values(n) 
                end do
             end do
             !!large values lead to bad looking contour plots
             if(maxValue .gt. this%maxValue) maxValue=this%maxValue
             if(minValue .lt. this%minValue) minValue=this%minValue
             
             title=""
             call OutputBuilder_make3DGnuplot( this%fileName(l), title, x_title, y_title, z_title, minValue, maxValue)
          else if( this%dimensions.eq.2) then
             if( this%type=="DENSITYPLOT") write (10,"(A10,A20,A20)") "#","X","Density"
             if( this%type=="ORBITALPLOT") write (10,"(A10,A20,A20)") "#","X","Orbital"
             n=0
             do i=0,numberOfSteps
                n=n+1
                if(abs(val%values(n))>1.0E-99_8) then
                   write (10,"(T10,E20.8,E20.8)") initialValue1+i*Vector_norm(this%step1), val%values(n) 
                else
                   write (10,"(T10,E20.8,E20.8)") initialValue1+i*Vector_norm(this%step1), 0.0
                end if
             end do
             if( this%type=="DENSITYPLOT") title=trim(nameOfSpecies)//" state "//trim(String_convertIntegerToString(state))//" density" 
             if( this%type=="ORBITALPLOT") title=trim(nameOfSpecies)//" state "//trim(String_convertIntegerToString(state))//" orbital"//trim(String_convertIntegerToString(orbital)) 
             call OutputBuilder_make2DGnuplot( this%fileName(l), title, x_title, y_title)
          end if
          close(10)
       end if
    end do

  end subroutine OutputBuilder_getPlot



  subroutine OutputBuilder_make2DGnuplot(fileName, title, x_title, y_title,&
       x_format, y_format, x_range, y_range, numOfGraphs)
    implicit none
    character(*) :: fileName
    character(*) :: title
    character(*) :: x_title
    character(*) :: y_title
    character(*), optional :: x_format
    character(*), optional :: y_format
    character(*), optional :: x_range
    character(*), optional :: y_range
    integer, optional :: numOfGraphs

    integer :: i
    character(20) :: charNumOfGraph
    character(20) :: auxXformat
    character(20) :: auxYformat
    character(20) :: auxXRange
    character(20) :: auxYRange

    integer :: auxNumOfGraphs

    auxXformat="%.1f"
    if(present(x_format)) auxXformat=trim(x_format)

    auxYformat="%3.1E"
    if(present(y_format)) auxYformat=trim(y_format)

    auxXRange=" [] "
    if(present(x_range)) auxXRange=' ['//trim(x_range)//'] '

    auxYRange="[] "
    if(present(y_range)) auxYRange='['//trim(y_range)//'] '

    auxNumOfGraphs=1
    if(present(numOfGraphs)) auxNumOfGraphs=numOfGraphs

    open ( 10,FILE=trim(fileName)//".gnp", STATUS='REPLACE',ACTION='WRITE')
    write (10,"(A)") 'set term post eps enh color dashed rounded dl 4 "Times-Bold" 15'
    write (10,"(A)") 'set output "'//trim(fileName)//'.eps"'
    write (10,"(A)") 'set encoding iso_8859_1'
    write (10,"(A)") 'set title "'//trim(title)//'"'
    write (10,"(A)") 'set xlabel "'//trim(x_title)//'"'
    write (10,"(A)") 'set format x "'//trim(auxXformat)//'"'
    write (10,"(A)") 'set ylabel "'//trim(y_title)//'"'
    write (10,"(A)") 'set format y "'//trim(auxYformat)//'"'
    if( auxNumOfGraphs >1) then
       write (10,"(A$)") 'plot '//trim(auxXRange)//trim(auxYRange)//' "'//trim(fileName)//'" using 1:2 w l title "" '
       do i=2, auxNumOfGraphs
          charNumOfGraph=String_convertIntegerToString(i+1)
          write (10,"(A$)") ', "'//trim(fileName)//'.dat"'//' using 1:'//trim(charNumOfGraph)//' w l  title "" '
       end do
       write (10,"(A)") ""
    else
       write (10,"(A)") 'plot '//trim(auxXRange)//trim(auxYRange)//' "'//trim(fileName)//'" w l title "" '
    end if
    write (10,"(A)") 'set output'
    close(10)

    !     status= system("gnuplot "//trim(fileName)//".gnp")
    call system("gnuplot "//trim(fileName)//".gnp")

  end subroutine OutputBuilder_make2DGnuplot


  subroutine OutputBuilder_make3DGnuplot(fileName, title, x_title, y_title, z_title, minValue, maxValue)

    implicit none
    character(*) :: fileName
    character(*) :: title
    character(*) :: x_title
    character(*) :: y_title
    character(*) :: z_title
    real(8) :: minValue
    real(8) :: maxValue

    open ( 100,FILE=trim(fileName)//'.gnp', STATUS='REPLACE',ACTION='WRITE')
    write (100,"(A)") 'set term post eps enh color "Helvetica" 16 size 7cm,5cm'
    write (100,"(A)") 'set encoding iso_8859_1'
    write (100,"(A)") 'set output "'//trim(fileName)//'.eps"'

    write (100,"(A)") 'set table "'//trim(fileName)//'.table"'
    write (100,"(A)") 'splot "'//trim(fileName)//'" u 1:2:3'
    write (100,"(A)") 'unset table'

    if(minValue.lt.0 .and. maxValue.gt.0) then
       write (100,"(A,I5)") 'levels=', 11
    else
       write (100,"(A,I5)") 'levels=', 10
    end if

    write (100,"(A,E20.8)") 'maxValue=', maxValue
    write (100,"(A,E20.8)") 'minValue=', minValue
    write (100,"(A)") 'step=(maxValue-minValue)/levels'

    write (100,"(A)") 'set contour base'
    write (100,"(A)") 'set cntrparam level incremental minValue, step , maxValue'
    write (100,"(A)") 'unset surface'

    write (100,"(A)") 'set table "'//trim(fileName)//'.cont"'
    write (100,"(A)") 'splot "'//trim(fileName)//'" u 1:2:3'
    write (100,"(A)") 'unset table'

    write (100,"(A)") 'reset'
    write (100,"(A)") 'unset key'

    write (100,"(A)") 'set cbrange [minValue:maxValue]'
    write (100,"(A)") 'set palette maxcolors levels'
    if(minValue.lt.0 .and. maxValue.gt.0) then
       write (100,"(A)") 'set cbtics (minValue, 0.0, maxValue)'
    else
       write (100,"(A)") 'set cbtics step'
    end if
    write (100,"(A)") 'set format cb "%3.1E"'

    if(minValue.lt.0 .and. maxValue.gt.0) then
       write (100,"(A)") 'set palette defined (minValue "blue", 0.0 "white", maxValue "red")'
    else if(minValue.ge.0) then
       write (100,"(A)") 'set palette defined (minValue "white", maxValue "red")'
    else
       write (100,"(A)") 'set palette defined (minValue "blue", maxValue "white")'
    end if

    write (100,"(A)") 'set grid front'

    write (100,"(A)") 'set format x "%.0f"'
    write (100,"(A)") 'set format y "%.0f"'
    write (100,"(A)") 'set xlabel "'//trim(x_title)//'"'
    write (100,"(A)") 'set ylabel "'//trim(y_title)//'"'

    write (100,"(A)")  'plot "'//trim(fileName)//'.table" with image, "'//trim(fileName)//'.cont" w l lt -1 lw 1.5'

    close(100)

    !     status= system("gnuplot "//trim(fileName)//".gnp")
    call system("gnuplot "//trim(fileName)//".gnp")

  end subroutine OutputBuilder_make3DGnuplot

end module OutputBuilder_

