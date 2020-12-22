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
     character(50),allocatable :: fileName(:)
     character(50) :: fileName2
     integer :: state
     integer :: orbital
     integer :: dimensions
     integer :: outputID
     integer :: auxID
     real(8) :: cubeSize
     type(vector) :: point1
     type(vector) :: point2
     type(vector) :: point3
     logical :: isInstanced
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
       OutputBuilder_make2DGraph, &
       OutputBuilder_make3DGraph, &
       OutputBuilder_get2DPlot, &
       OutputBuilder_get3DPlot, &
       OutputBuilder_getDensityPlot, &
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
  subroutine OutputBuilder_constructor(this, ID, type ,species, state, orbital, dimensions, cubeSize, point1, point2, point3  )
    type(OutputBuilder) :: this
    integer :: ID
    character(*) :: type
    character(*) :: species
    integer,optional :: state
    integer,optional :: orbital
    integer,optional :: dimensions
    real(8),optional :: cubeSize
    type(Vector),optional :: point1
    type(Vector),optional :: point2
    type(Vector),optional :: point3
    integer :: i
    
    this%type=type
    ! print *, "this%type", this%type
    this%outputID=ID
    this%species=trim(String_getUppercase(species))
    if( trim(this%species) .eq. "ALL" ) then
       allocate(this%fileName(MolecularSystem_getNumberOfQuantumSpecies()))
    else
       allocate(this%fileName(1))
    end if
    this%state=1
    if( present(state)) this%state=state
    this%orbital=0
    if( present(orbital)) this%orbital=orbital
    this%dimensions=0
    if( present(dimensions)) this%dimensions=dimensions
    this%cubeSize=0
    if( present(cubeSize)) this%cubeSize=cubeSize
    
    call Vector_constructor(this%point1, 3, 0.0_8 )
    call Vector_constructor(this%point2, 3, 0.0_8 )
    call Vector_constructor(this%point3, 3, 0.0_8 )

    if( present(point1)) this%point1%values=point1%values
    if( present(point2)) this%point2%values=point2%values
    if( present(point3)) this%point3%values=point3%values

    if ( trim(CONTROL_instance%UNITS) == "ANGS") then
       this%point1%values= this%point1%values / AMSTRONG
       this%point2%values= this%point2%values / AMSTRONG
       this%point3%values= this%point3%values / AMSTRONG
    end if

    this%auxID=1
    !!Check for other outputs of the same type
    do i=1, this%outputID-1
       if( trim(outputs_instance(i)%type) .eq. trim(this%type)  ) this%auxID=this%auxID+1
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

  !!>
  !! @brief Indica si el objeto ha sido instanciado o no
  !!
  !<
!  function OutputBuilder_isInstanced( this ) result( output )
!    implicit  none
!    type(OutputBuilder), intent(in) :: this
!    logical :: output
!
!    output = this%isInstanced
!
!  end function OutputBuilder_isInstanced

  !>
  !! @brief  Maneja excepciones de la clase
  !<
  subroutine OutputBuilder_exception( typeMessage, description, debugDescription)
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

  end subroutine OutputBuilder_exception


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

    do l=1,size(this%fileName)
       write (*,"(A20,A)") "FileName: ", this%fileName(l)
    end do
    ! TODO Fix this line.
    ! if (this%filename2 /= "") print *, "FileName 2: ", this%fileName2
    if (this%species /= "ALL") write (*,"(A20,A10)") "for species: ", this%species
    if (this%orbital /= 0) write (*,"(A20,I10)") "for orbital: ", this%orbital
    if (this%state /= 1) write (*,"(A20,I10)") "for excited state: ", this%state
    if (this%dimensions /= 0) write (*,"(A20,I2)") "dimensions: ", this%dimensions
    if (this%cubeSize /= 0.0_8) write (*,"(A20,F15.5)") "cube size (a.u.): ", this%cubeSize
    if (this%dimensions >= 1) write (*,"(A20,F10.5,F10.5,F10.5)") "Point 1 (a.u.): ", this%point1%values(1), this%point1%values(2), this%point1%values(3)
    if (this%dimensions >= 2) write (*,"(A20,F10.5,F10.5,F10.5)") "Point 2 (a.u.): ", this%point2%values(1), this%point2%values(2), this%point2%values(3)
    if (this%dimensions >= 3) write (*,"(A20,F10.5,F10.5,F10.5)") "Point 3 (a.u.): ", this%point3%values(1), this%point3%values(2), this%point3%values(3)
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

     case ( "moldenFile") 
        call OutputBuilder_writeMoldenFile (this)

     case ("VecGamessFile")
        call OutputBuilder_VecGamessFile (this)

     case ("casinoFile")
        call OutputBuilder_casinoFile (this)

     case ("EigenGamessFile")
        call OutputBuilder_writeEigenvalues (this)

     case ("fchkFile")
        call OutputBuilder_writeFchkFile (this)
        
     case ( "wfnFile") 
        call OutputBuilder_writeMoldenFile (this)
        call OutputBuilder_generateAIMFiles (this)

    case ( "NBO47File") 
        call OutputBuilder_writeMoldenFile (this)
        call OutputBuilder_generateAIMFiles (this)

    case ( "wfxFile" ) 

        call OutputBuilder_writeMoldenFile (this)
        call OutputBuilder_generateAIMFiles (this)

    case ( "extendedwfnFile") 
        call OutputBuilder_writeMoldenFile (this)
        call OutputBuilder_generateAIMFiles (this)
        call OutputBuilder_generateExtendedWfnFile (this)

    case ( "densityPlot") 
        if (this%dimensions == 2) call OutputBuilder_getDensityPlot(this)
        if (this%dimensions == 3) call OutputBuilder_getDensityPlot(this)

   case ( "densityCube") 
      call OutputBuilder_getDensityCube(this)
!
   case ( "orbitalPlot") 
      if (this%dimensions == 2) call OutputBuilder_get2DPlot(this)
      if (this%dimensions == 3) call OutputBuilder_get3DPlot(this)
!
!     case ( "orbitalCube") 
!        call OutputBuilder_getCube(this)
!
!     case ( "fukuiPlot") 
!        if (this%dimensions == 2) call OutputBuilder_get2DPlot(this)
!        if (this%dimensions == 3) call OutputBuilder_get3DPlot(this)
!
!     case ( "fukuiCube") 
!        call OutputBuilder_getCube(this)
!
     case default
        call OutputBuilder_exception(ERROR, "The output type you requested has not been implemented yet", "OutputBuilder_buildOutput" )

     end select
   end subroutine OutputBuilder_buildOutput


  
  subroutine OutputBuilder_writeMoldenFile(this)
    implicit none
    type(OutputBuilder) :: this
    type(MolecularSystem) :: MolecularSystemInstance

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: numberOfSpecies
    integer :: state,numberOfStates
    real :: occupation
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
    integer :: wfnUnit, occupationsUnit
    character(50) :: wfnFile, occupationsFile, fileName
    integer :: numberOfContractions
    character(50) :: arguments(20)
    character(19) , allocatable :: labelsOfContractions(:)
    integer :: counter, auxcounter
    character(6) :: nickname
    character(4) :: shellCode
    character(2) :: space
    integer :: totalNumberOfParticles, n
    logical :: existFile
    
    !     if ( CONTROL_instance%ARE_THERE_DUMMY_ATOMS ) then
    !        auxString=MolecularSystem_getNameOfSpecie( 1 )
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

    occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    inquire(FILE = occupationsFile, EXIST = existFile )

    !! Check if there are CI fractional occupations or build the occupations vector
    if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE"  .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0 .and. existFile) then

       print *, "              We are printing the molden files for the CI states!"
       
       numberOfStates=CONTROL_instance%CI_STATES_TO_PRINT
       allocate(fractionalOccupations(numberOfSpecies,numberOfStates))
       allocate(energyOfMolecularOrbital(numberOfSpecies,numberOfStates))
       allocate(coefficientsOfCombination(numberOfSpecies,numberOfStates))
       occupationsUnit = 29

       open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")
       do state=1,numberOfStates
          do l=1,numberOfSpecies
             write(auxstring,*) state
             
             arguments(1) = "OCCUPATIONS"//trim(adjustl(auxstring))
             arguments(2) = MolecularSystem_getNameOfSpecie( l )
             call  Vector_getFromFile(elementsNum=MolecularSystem_getTotalNumberOfContractions(l),&
                  unit=occupationsUnit,&
                  arguments=arguments(1:2),&
                  output=fractionalOccupations(l,state))
             
             arguments(1) = "NATURALORBITALS"//trim(adjustl(auxstring)) 
             coefficientsOfCombination(l,state)=Matrix_getFromFile(unit=occupationsUnit,&
                  rows= int(MolecularSystem_getTotalNumberOfContractions(l),4), &
                  columns= int(MolecularSystem_getTotalNumberOfContractions(l),4), &
                  arguments=arguments(1:2))

             call Vector_constructor( energyOfMolecularOrbital(l,state), MolecularSystem_getTotalNumberOfContractions(l) )
             energyOfMolecularOrbital(l,state)%values=0.0
             
          end do
       end do
       close(occupationsUnit)
       
    else
    !! Open file for wavefunction and load results                                                                                     
    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")
       
    numberOfStates=1
    allocate(fractionalOccupations(numberOfSpecies,1))
    allocate(energyOfMolecularOrbital(numberOfSpecies,1))
    allocate(coefficientsOfCombination(numberOfSpecies,1))
    do l=1,numberOfSpecies
       call Vector_constructor( fractionalOccupations(l,1), &
            MolecularSystem_getTotalNumberOfContractions(l) )
       fractionalOccupations(l,1)%values=0.0
       do i=1, MolecularSystem_getOcupationNumber(l)
          fractionalOccupations(l,1)%values(i)=1.0_8 * MolecularSystem_getLambda(l)
       end do
       arguments(2) = MolecularSystem_getNameOfSpecie(l)
       arguments(1) = "COEFFICIENTS"
       coefficientsOfcombination = &
            Matrix_getFromFile(unit=wfnUnit, &
            rows= int(MolecularSystem_getTotalNumberOfContractions(l),4), &
            columns= int(MolecularSystem_getTotalNumberOfContractions(l),4),&
            binary=.true., &
            arguments=arguments(1:2))

       arguments(1) = "ORBITALS"
       call Vector_getFromFile( elementsNum = MolecularSystem_getTotalNumberOfContractions(l), &
            unit = wfnUnit,&
            binary = .true.,&
            arguments = arguments(1:2), &
            output = energyOfMolecularOrbital(l,1) )

    end do
    close(wfnUnit)
       
    end if
    


    do state=1,numberOfStates
       do l=1,numberOfSpecies

          if (state .eq. 1) then
             auxString=MolecularSystem_getNameOfSpecie( l )
          else
             write(auxString, "(I8)")  state
             auxString=trim(MolecularSystem_getNameOfSpecie( l ))//"-"//trim( adjustl(auxString))
          end if
          
          this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".molden"

          totalNumberOfParticles = 0

          open(10,file=this%fileName(l),status='replace',action='write')
          write(10,"(A)") "[Molden Format]"
          ! if ( CONTROL_instance%UNITS=="ANGS") then
          !   write(10,"(A)") "[Atoms] Angs"
          ! else 
          write(10,"(A)") "[Atoms] AU"
          ! end if

          auxMatrix%values=0.0
          j=0
          do i=1, size(MolecularSystem_instance%species(l)%particles)
             j=j+1
             origin = MolecularSystem_instance%species(l)%particles(j)%origin 
             auxMatrix%values(j,:)=origin
             symbol=MolecularSystem_instance%species(l)%particles(j)%nickname

             if(scan(symbol,"_") /=0) symbol=symbol(1:scan(symbol,"_")-1)
             if(scan(symbol,"[") /=0) symbol=symbol(scan(symbol,"[")+1:scan(symbol,"]")-1)

             if ( CONTROL_instance%UNITS=="ANGS") origin = origin * AMSTRONG

             totalNumberOfParticles = totalNumberOfParticles + 1

             write (10,"(A,I8,I8,F15.8,F15.8,F15.8)") trim(symbol), j,&
                  int(abs(MolecularSystem_instance%species(l)%particles(j)%totalCharge)), origin(1), origin(2), origin(3)

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
                  if ( CONTROL_instance%UNITS=="ANGS") origin = origin * AMSTRONG
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

             occupation=fractionalOccupations(l,state)%values(j)

             i = 0
             do k=1,size(coefficientsOfCombination(l,state)%values,dim=1)
                i = i + 1
                write(10,"(I4,A2,F15.10)") k,"  ", coefficientsOfCombination(l,state)%values(k,j)
             end do

              if ( totalNumberOfParticles > size(MolecularSystem_instance%species(l)%particles) ) then
                if ( CONTROL_instance%MOLDEN_FILE_FORMAT == "MIXED" ) then
                  do n = 1, ( totalNumberOfParticles - size(MolecularSystem_instance%species(l)%particles) )
                    write(10,"(I4,A2,F15.10)") i+n,"  ", 0.0_8
                  end do
                end if
              end if

          end do


          close(10)
       end do
    end do

    call Matrix_destructor( localizationOfCenters )
    call Matrix_destructor( auxMatrix )
    deallocate(labels)

    !     end if

  end subroutine OutputBuilder_writeMoldenFile


!!!!!!!!!!!!!!!!!!!!!!!!!GAMESS .VEC FILE LAURA 

  subroutine OutputBuilder_VecGamessFile(this)
    implicit none
    type(OutputBuilder) :: this
    type(MolecularSystem) :: MolecularSystemInstance

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: specieID
    logical :: wasPress
    character(10) :: auxString
    character(10) :: symbol
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Matrix) :: coefficientsOfcombination
    character(10),allocatable :: labels(:)
    integer :: wfnUnit
    character(50) :: wfnFile
    integer :: numberOfContractions
    character(50) :: arguments(20)
    character(19) , allocatable :: labelsOfContractions(:)
    integer :: counter, auxcounter
    character(6) :: nickname
    character(4) :: shellCode
    character(2) :: space
    integer :: totalNumberOfParticles, n

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

!! Open file for wavefunction                                                                                     
        open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

        do l=1,MolecularSystem_getNumberOfQuantumSpecies()

           auxString=MolecularSystem_getNameOfSpecie( l )

       this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".vec"

       open(29,file=this%fileName(l),status='replace',action='write')

           !! Build a vector of labels of contractions
     if(allocated(labelsOfContractions)) deallocate(labelsOfContractions)
           allocate(labelsOfContractions(numberOfContractions))

           labelsOfContractions =  MolecularSystem_getlabelsofcontractions( specieID )

           !! Swap some columns according to the molden format
           do k=1,size(coefficientsOfCombination%values,dim=1)
    !! Take the shellcode
                read (labelsOfContractions(k), "(I5,A2,A6,A2,A4)"), counter, space, nickname, space, shellcode 

    !! Reorder the D functions
                !! counter:  1,  2,  3,  4,  5,  6
                !! Lowdin:  XX, XY, XZ, YY, YZ, ZZ
                !! Molden:  XX, YY, ZZ, XY, XZ, YZ 
                !!  1-1, 2-4, 3-5, 4-2, 5-6, 6-3
                !!  2-4, 3-5, 5-6

    if ( shellcode == "Dxx" ) then 
        auxcounter = counter
        !! Swap XY and YY
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+3)
        !! Swap XZ and ZZ
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+5)
        !! Swap YZ and XZ'
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
                end if

    !! Reorder the F functions
                !! counter:   1,   2,   3,   4,   5,   6,   7,   8    9,  10
                !! Lowdin:  XXX, XXY, XXZ, XYY, XYZ, XZZ, YYY, YYZ, YZZ, ZZZ
                !! Molden:  XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ
                !! Gamess:  XXX, YYY, ZZZ, XXY, XXZ, XYY, YYZ, XZZ, YZZ, XYZ
                
              if ( shellcode == "Fxxx" ) then 
        auxcounter = counter
                    call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+9)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)
                  call Matrix_swapRows(  coefficientsOfCombination, auxcounter+3 , auxcounter+4)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+5)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+8 , auxcounter+6)
              call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)
                end if

      end do


     
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
    type(MolecularSystem) :: MolecularSystemInstance

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: g, h, m
    integer :: specieID
    logical :: wasPress
    character(10) :: auxString
    character(10) :: symbol
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Matrix) :: coefficientsOfcombination
    real(8), allocatable :: superMatrix(:,:)
    character(10),allocatable :: labels(:)
    integer :: wfnUnit
    character(50) :: wfnFile
    integer :: numberOfContractions, superSize
    integer :: numberOfContractionsA, numberOfContractionsB
    integer :: numberOfShellsA, numberOfShellsB, totalShells
    character(50) :: arguments(20)
    character(19) , allocatable :: labelsOfContractions(:)
    integer :: counter, auxcounter
    character(6) :: nickname
    character(2) :: space
    integer :: i0, j0, maxl, shellCode
    integer :: totalNumberOfParticles, n
    real(8) :: puntualInteractionEnergy

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

    !! Open file for wavefunction                                                                                     
    open(unit = wfnUnit, file = trim(wfnFile), status = "old", form = "unformatted")
 

    this%fileName = trim(CONTROL_instance%INPUT_FILE)//"casino"
    open(29,file=this%fileName,status='replace',action='write')

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
        call OutputBuilder_exception(ERROR, "The maximum number of quantum species cannot be greater than two", "OutputBuilder_casinoFile" )
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
      arguments(2) = MolecularSystem_getNameOfSpecie(specieID)
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

    call OutputBuilder_exception(WARNING, "The order of the coefficients only works until P orbitals", "OutputBuilder_casinoFile" )
        
  end subroutine OutputBuilder_casinoFile
  
  !!Escribe los valores propios en el archivo eigenvalues.dat para que puedan ser leidos por GAMESS   Laura

   subroutine OutputBuilder_writeEigenvalues(this)
    implicit none
    type(OutputBuilder) :: this
    type(MolecularSystem) :: MolecularSystemInstance

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: m
    integer :: specieID
    real :: occupation
    integer :: occupationTotal
    logical :: wasPress
    character(10) :: auxString
    character(10) :: symbol
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Vector) :: energyOfMolecularOrbital
    type(Matrix) :: coefficientsOfcombination
    character(10),allocatable :: labels(:)
    integer :: wfnUnit
    character(50) :: wfnFile
    integer :: numberOfContractions
    character(50) :: arguments(20)
    character(19) , allocatable :: labelsOfContractions(:)
    integer :: counter, auxcounter
    character(6) :: nickname
    character(4) :: shellCode
    character(2) :: space
    integer :: totalNumberOfParticles, n

    ! auxString="speciesName"

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

        ! auxString=MolecularSystem_getNameOfSpecie( 1 )
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

           auxString=MolecularSystem_getNameOfSpecie( l )
           specieID = MolecularSystem_getSpecieID(auxString)
           this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".eigen"
           open(129,file=this%fileName(l),status='replace',action='write')

            specieID = int( MolecularSystem_getSpecieID(nameOfSpecie = trim(auxString)) )
            numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
            arguments(2) = MolecularSystem_getNameOfSpecie(specieID)

           arguments(1) = "ORBITALS"
           call Vector_getFromFile( elementsNum = numberOfContractions, &
                unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                output = energyOfMolecularOrbital )

           do j=1,size(energyOfMolecularOrbital%values)
              write (129,"(F15.12)") ,energyOfMolecularOrbital%values(j)
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
    type(MolecularSystem) :: MolecularSystemInstance

    integer :: h,i,j,k,l,n,m
    integer :: numberOfSpecies
    integer :: state,numberOfStates
    character(50) :: auxString
    real(8) :: origin(3)
    real(8), allocatable :: charges(:)
    type(Matrix) :: localizationOfCenters
    type(Matrix) :: auxMatrix
    type(Vector) :: energyOfMolecularOrbital
    type(Matrix) :: coefficientsOfcombination
    type(Matrix),allocatable :: fractionalOccupations(:)
    character(10),allocatable :: labels(:)

    integer :: wfnUnit, occupationsUnit
    character(50) :: wfnFile, occupationsFile, fileName
    integer :: numberOfBetaElectrons
    integer :: numberOfPrimitives
    integer :: numberOfContractions
    integer :: numberOfShells
    integer :: numberOfAtoms
    integer :: numberOfAtomShells
    integer :: numberOfShellPrimitives
    real(8) :: particlesPerOrbital
    type(matrix) :: densityMatrix
    real(8) :: densityElement
    character(50) :: nameOfSpecies
    character(40) :: header
    character(50) :: arguments(20)
    logical :: existFile


    localizationOfCenters=ParticleManager_getCartesianMatrixOfCentersOfOptimization()
    auxMatrix=localizationOfCenters
    allocate( labels( size(auxMatrix%values,dim=1) ) )
    allocate( charges( size(auxMatrix%values,dim=1) ) )
    labels=ParticleManager_getLabelsOfCentersOfOptimization()
    charges=ParticleManager_getChargesOfCentersOfOptimization()
    numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()

       
    !! Check if there are CI fractional occupations or build the occupations vector
    ! allocate(fractionalOccupations(numberOfSpecies))

    ! occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    ! inquire(FILE = occupationsFile, EXIST = existFile )

    ! if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE"  .and. CONTROL_instance%CI_STATES_TO_PRINT .gt. 0 .and. existFile) then

    !    print *, "              We are printing the fchk files for the CI states!"
       
    !    numberOfStates=CONTROL_instance%CI_STATES_TO_PRINT
    !    occupationsUnit = 29

    !    open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")
    !    do l=1,numberOfSpecies
    !       arguments(1) = "OCCUPATIONS"
    !       arguments(2) = MolecularSystem_getNameOfSpecie( l )
    !       fractionalOccupations(l)= Matrix_getFromFile(unit=occupationsUnit,&
    !            rows=int(MolecularSystem_getTotalNumberOfContractions(l),4),&
    !            columns=int(numberOfStates,4),&
    !            arguments=arguments(1:2))
    !    end do
    !    close(occupationsUnit)     
    ! else
    !    numberOfStates=1
    !    do l=1,numberOfSpecies
    !       call Matrix_constructor( fractionalOccupations(l), int(MolecularSystem_getTotalNumberOfContractions(l),8), int(numberOfStates,8), 0.0_8)
    !       do i=1, MolecularSystem_getOcupationNumber(l)
    !          fractionalOccupations(l)%values(i,1)=1.0_8 * MolecularSystem_getLambda(l)
    !       end do
    !    end do
    ! end if
    

    !! Open file for wavefunction                                                                                     
    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

    ! do state=1,numberOfStates
    do l=1,numberOfSpecies
       nameOfSpecies=MolecularSystem_getNameOfSpecie(l)
       particlesPerOrbital=MolecularSystem_getLambda(l)
       ! if (state .eq. 1) then
       !    auxString=nameOfSpecies
       ! else
       !    write(auxString, "(I8)")  state
       !    auxString=trim(nameOfSpecies)//"-"//trim( adjustl(auxString))
       ! end if

       ! this%fileName=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".fchk"
       this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".fchk"

       numberOfAtoms=size(MolecularSystem_instance%species(l)%particles)
       numberOfShells=MolecularSystem_getNumberOfContractions(l)
       numberOfContractions=MolecularSystem_getTotalNumberOfContractions(l)

       arguments(2) = MolecularSystem_getNameOfSpecie(l)
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
             write(10,"(I12)") MolecularSystem_instance%species(l)%particles(j)%internalSize
          else
             write(10,"(I12)", advance="no") MolecularSystem_instance%species(l)%particles(j)%internalSize
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

       
       call MolecularSystem_changeOrbitalOrder( coefficientsOfCombination, l, "LOWDIN", "MOLDEN" )

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
    type(MolecularSystem) :: MolecularSystemInstance
    character(50) :: auxString
    character(50) :: initialSettingsFile
    character(50) :: moldenFileName
    integer :: l
    character(50) :: wfnFile
    character(2) :: wfnStatus, wfxStatus, nboStatus
    character(10) :: extension
    character(50) :: arguments(20)
    integer :: wfnUnit
    real(8) :: totalEnergy, virial

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
    initialSettingsFile = "m2a.ini"

    select case (this%type) 
      case ( "wfnFile" )
    wfnStatus="1"
    nboStatus="-1"
    wfxStatus="-1"
        extension=".wfn"
      case ( "NBO47File" )
    wfnStatus="-1"
    nboStatus="1"
    wfxStatus="-1"
    extension=".47"
      case ( "wfxFile" ) 
    wfnStatus="-1"
    nboStatus="-1"
      wfxStatus="1"
    extension=".wfx"
  case ( "extendedwfnFile" )
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
        auxString=MolecularSystem_getNameOfSpecie( l )
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
    character(50) :: initialWfnFile
    character(50) :: auxString

    do l=1,MolecularSystem_getNumberOfQuantumSpecies()
        auxString=MolecularSystem_getNameOfSpecie( l )
    end do

  end subroutine OutputBuilder_generateExtendedWfnFile

   subroutine OutputBuilder_get3DPlot(this)
     type(OutputBuilder) :: this
     character(50) :: outputID, auxID
     character(50) :: orbitalNum

     integer :: speciesID
     character(50) :: nameOfSpecies
     integer :: i,j
     integer :: numberOfSteps
     type(vector) :: step1
     type(vector) :: step2
     real(8) :: val, val2
     real(8) :: maxValue, maxValue2
     real(8) :: minValue, minValue2
     real(8) :: coordinate(3)

     character(50) :: title, title2
     character(50) :: x_title
     character(50) :: y_title
     character(50) :: z_title
   
     call Vector_Constructor(step1, 3)
     call Vector_Constructor(step2, 3)
     speciesID = MolecularSystem_getSpecieIDFromSymbol( trim(this%species) )
     nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

     this%fileName2=""
     numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
     step1%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
     step2%values(:)=(this%point3%values(:)-this%point1%values(:))/numberOfSteps
     outputID=String_convertIntegerToString(this%outputID)
     auxID=String_convertIntegerToString(this%auxID)

     x_title="x/a.u."
     y_title="y/a.u."
     z_title=""

     select case( this%type )

     case ( "orbitalPlot") 
        orbitalNum=String_convertIntegerToString(this%orbital)
        this%fileName(1)=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".3D.orb"//trim(orbitalNum)
        open(10,file=this%fileName(1),status='replace',action='write')
        write (10,"(A10,A20,A20,A20)") "#", "X","Y","OrbitalValue"
        title=trim(this%species)//" Orbital Number: "//trim(orbitalNum) 

     case ( "fukuiPlot") 
        this%fileName(1)=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".3D.fkpos"
        this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".3D.fkneg"
        open(10,file=this%fileName(1),status='replace',action='write')
        write (10,"(A10,A20,A20,A20)") "#", "X","Y","PositiveFukuiValue"
        title=trim(this%species)//" Positive Fukui"
        open(11,file=this%fileName2,status='replace',action='write')
        write (11,"(A10,A20,A20,A20)") "#", "X","Y","NegativeFukuiValue"
        title2=trim(this%species)//" Negative Fukui"

     case default
        call OutputBuilder_exception(ERROR, "The output plot type you requested has not been implemented yet", "OutputBuilder_get3DPlot" )

     end select

     val=0.0_8     
     val2=0.0_8     
     maxValue=0.0_8
     minValue=0.0_8 
     maxValue2=0.0_8
     minValue2=0.0_8 
     do i=0,numberOfSteps
        write (10,*) ""
        if (this%type .eq. "fukuiPlot") write(11,*) ""
        do j=0,numberOfSteps
           coordinate(:)=i*step1%values(:)+j*step2%values(:)+this%point1%values(:)
           select case( this%type )
           case ( "orbitalPlot") 
              val=CalculateWaveFunction_getOrbitalValueAt(nameOfSpecies, this%orbital, coordinate )  
           case ( "fukuiPlot") 
!!              val=CalculateProperties_getFukuiAt( this%species, "positive", coordinate )  
!!              val2=CalculateProperties_getFukuiAt( this%species, "negative", coordinate )  
           case default
           end select
           write (10,"(T10,F20.8,F20.8,F20.8)") i*Vector_norm(step1),j*Vector_norm(step2),val 
           if (val > maxValue) maxValue = val
           if (val < minValue) minValue = val
           if (this%type .eq. "fukuiPlot" ) then
              write (11,"(T10,F20.8,F20.8,F20.8)") i*Vector_norm(step1),j*Vector_norm(step2),val2 
              if (val2 > maxValue2) maxValue2 = val2
              if (val2 < minValue2) minValue2 = val2
           end if
           ! print *, coordinate, val
        end do
     end do

     call OutputBuilder_make3DGraph( this%fileName(1), title, x_title, y_title, z_title, minValue, maxValue)
     close(10)

     if (this%type .eq. "fukuiPlot" ) then
        call OutputBuilder_make3DGraph( this%fileName2, title2, x_title, y_title, z_title, minValue2, maxValue2)
        close(11)
     end if

     call Vector_Destructor(step1)
     call Vector_Destructor(step2)

   end subroutine OutputBuilder_get3DPlot

   subroutine OutputBuilder_get2DPlot(this)
     implicit none
     type(outputBuilder) :: this
     character(50) :: outputID, auxID
     character(50) :: orbitalNum

     integer :: i, speciesID
     character(50) :: nameOfSpecies

     integer :: numberOfSteps
     type(vector) :: step
     real(8) :: val, val2
     real(8) :: coordinate(3)

     character(50) :: title
     character(50) :: x_title
     character(50) :: y_title

     call Vector_Constructor(step, 3)
     speciesID = MolecularSystem_getSpecieIDFromSymbol( trim(this%species) )
     nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)

     this%fileName2=""
     numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
     step%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
     outputID=String_convertIntegerToString(this%outputID)
     auxID=String_convertIntegerToString(this%auxID)

     x_title="distance/a.u."
     select case( this%type )

     case ( "orbitalPlot") 
        orbitalNum=String_convertIntegerToString(this%orbital)
        this%fileName(1)=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".2D.orb"//trim(orbitalNum)
        open(10,file=this%fileName(1),status='replace',action='write')
        write (10,"(A10,A20,A20)") "#", "X","OrbitalValue"
        title=trim(this%species)//" Orbital Number "//trim(orbitalNum) 
        y_title="orbitalValue/a.u.^{-3/2}"

     case ( "fukuiPlot") 
        this%fileName(1)=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".2D.fkpos"
        this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".2D.fkneg"

        open(10,file=this%fileName(1),status='replace',action='write')
        write (10,"(A10,A20,A20)") "#","X","PositiveFukuiValue"
        title=trim(this%species)//" positive fukui" 
        y_title="density/a.u.^{-3}"

        open(11,file=this%fileName2,status='replace',action='write')
        write (11,"(A10,A20,A20)") "#","X","NegativeFukuiValue"
     
     case default
        call OutputBuilder_exception(ERROR, "The output plot type you requested has not been implemented yet", "OutputBuilder_get3DPlot" )

     end select

     do i=0,numberOfSteps
        coordinate(:)=i*step%values(:)+this%point1%values(:)
        select case( this%type )
        case ( "orbitalPlot") 
           val=CalculateWaveFunction_getOrbitalValueAt(nameOfSpecies, this%orbital, coordinate )  
        case ( "fukuiPlot") 
!!           val=CalculateProperties_getFukuiAt( this%species, "positive", coordinate )  
!!           val2=CalculateProperties_getFukuiAt( this%species, "negative", coordinate )  
        case default
        end select
        write (10,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val 
        if (this%type .eq. "fukuiPlot") write (11,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val2 
     end do

     close(10)

     call OutputBuilder_make2DGraph( this%fileName(1), title, x_title, y_title)
!!     if (this%type .eq. "fukuiPlot") then
!!        close(11)
!!        title=trim(this%species)//" negative fukui" 
!!        call OutputBuilder_make2DGraph( this%fileName2, title, x_title, y_title)
!!     end if
     call Vector_Destructor ( step)

   end subroutine OutputBuilder_get2DPlot


  subroutine OutputBuilder_getDensityCube(this )
    implicit none
    type(OutputBuilder) :: this
    character(50) :: outputID, auxID
    real(8):: cubeSize

    integer :: i, j, k, n, w, natom
    integer :: atomicCharge
    integer :: speciesID
    integer :: numberOfSteps
    real(8) :: step
    real(8) :: lowerLimit(3)
    real(8), allocatable :: val(:)
    real(8) :: coordinate(3)

    integer :: wfnunit, occupationsUnit 
    integer :: numberOfOrbitals
    type(matrix) :: densityMatrix

    character(50) :: arguments(20), wfnFile, occupationsFile, auxstring, nameOfSpecies
    logical :: existFile

    !Writes Gaussian Cube 
    
    speciesID = MolecularSystem_getSpecieIDFromSymbol( trim(this%species) )

    nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
    numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(speciesID)

    outputID=String_convertIntegerToString(this%outputID)
    auxID=String_convertIntegerToString(this%auxID)
  
    ! Check if there are CI density matrices and read those or the HF matrix
    occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
    inquire(FILE = occupationsFile, EXIST = existFile )
    
    if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE"  .and. existFile ) then
       print *, "We are printing a density file for ", trim(nameOfSpecies), " in the CI state No. ", this%state

       occupationsUnit = 29

       open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")


       write(auxstring,*) this%state
       arguments(2) = nameOfSpecies
       arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 

       densityMatrix= Matrix_getFromFile(unit=occupationsUnit, rows= int(numberOfOrbitals,4), &
            columns= int(numberOfOrbitals,4), binary=.false., arguments=arguments(1:2))


       close(occupationsUnit)     
    else

       !! Read density matrix
       !! Open file for wavefunction
       wfnFile = "lowdin.wfn"
       wfnUnit = 20
       open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

       arguments(2) = nameOfSpecies
       arguments(1) = "DENSITY"

       densityMatrix = &
            Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
            columns=int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))

       close (wfnUnit)

    end if

    
    this%fileName(1)=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%species)//".dens.cub"
    open(10,file=this%fileName(1),status='replace',action='write')

    lowerLimit(:)=this%point1%values(:)-this%cubeSize/2
    numberOfSteps=CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
    step= this%cubeSize/numberOfSteps

    allocate (val (numberOfSteps) )

    ! do n=1, size(MolecularSystem_instance%particlesPtr)
    !    if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-" .or. &
    !         trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-ALPHA" .and. &
    !         MolecularSystem_instance%particlesPtr(k)%isQuantum ) then
    !       natom = natom +1
    !    end if
    ! end do
    natom=1
    
    write (10,"(A)") "Gaussian Cube generated with Lowdin Software"
    write (10,"(A)") this%fileName(1)
    write (10,"(I8,F20.8,F20.8,F20.8,I8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3), 1
    write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, step, 0.0, 0.0
    write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, 0.0, step, 0.0
    write (10,"(I8,F20.8,F20.8,F20.8)") numberOfSteps, 0.0, 0.0, step

    write (10, "(I8,I8,F20.8,F20.8,F20.8)") &
         1, 1, this%point1%values
    ! do n=1, size(MolecularSystem_instance%particlesPtr)
    !    if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-" .or. &
    !         trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-ALPHA" .and. &
    !         MolecularSystem_instance%particlesPtr(n)%isQuantum ) then
    !       atomicCharge=-MolecularSystem_instance%particlesPtr(n)%totalCharge
    !       write (10, "(I8,F20.8,F20.8,F20.8,F20.8)") &
    !            atomicCharge, 0.0, MolecularSystem_instance%particlesPtr(n)%origin(1:3)
    !    end if
    ! end do

    do i=1,numberOfSteps
       coordinate(1)=lowerLimit(1)+(i-1)*step
       do j=1, numberOfSteps
          coordinate(2)=lowerLimit(2)+(j-1)*step
          do k=1, numberOfSteps
             coordinate(3)=lowerLimit(3)+(k-1)*step

             val(k)=CalculateWaveFunction_getDensityAt( nameOfSpecies, coordinate, densityMatrix )
          end do
          write(10,*) ( val(w) , w=1,numberOfSteps )
          write(10,*) ( "" )
       end do
    end do

    deallocate (val)
    close(10)
    
  end subroutine OutputBuilder_getDensityCube

  subroutine OutputBuilder_getDensityPlot(this)
     type(OutputBuilder) :: this
     character(50) :: outputID, auxID
     character(50) :: orbitalNum

     integer :: i,j,l, speciesID, wfnunit, occupationsUnit 
     integer :: numberOfSteps, numberOfOrbitals
     integer :: numberOfSpecies
     type(vector) :: step1, step2
     type(matrix) :: densityMatrix, auxMatrix
     real(8) :: val, maxValue, minValue
     real(8) :: coordinate(3), plotDistance1, plotDistance2

     character(50) :: arguments(20), wfnFile, occupationsFile, auxstring, nameOfSpecies
     character(50) :: title, x_title, y_title, z_title
     logical :: existFile

     call Vector_Constructor(step1, 3)
     call Vector_Constructor(step2, 3)

     numberOfSpecies=MolecularSystem_getNumberOfQuantumSpecies()
     outputID=String_convertIntegerToString(this%outputID)
     auxID=String_convertIntegerToString(this%auxID)

     l=0
     do speciesID=1, numberOfSpecies
        nameOfSpecies=MolecularSystem_getNameOfSpecie(speciesID)
        if(trim(this%species) .eq. trim(nameOfSpecies) .or. trim(this%species) .eq. "ALL" ) then
           l=l+1   
           numberOfOrbitals=MolecularSystem_getTotalNumberOfContractions(speciesID)

           occupationsFile = trim(CONTROL_instance%INPUT_FILE)//"Matrices.ci"
           inquire(FILE = occupationsFile, EXIST = existFile )

           ! Check if there are CI density matrices and read those or the HF matrix
           if ( CONTROL_instance%CONFIGURATION_INTERACTION_LEVEL /= "NONE"  .and. existFile) then
              print *, "We are printing a density file for ", trim(nameOfSpecies), " in the CI state No. ", this%state

              occupationsUnit = 29

              open(unit = occupationsUnit, file=trim(occupationsFile), status="old", form="formatted")

              write(auxstring,*) this%state
              arguments(2) = nameOfSpecies
              arguments(1) = "DENSITYMATRIX"//trim(adjustl(auxstring)) 

              densityMatrix= Matrix_getFromFile(unit=occupationsUnit, rows= int(numberOfOrbitals,4), &
                   columns= int(numberOfOrbitals,4), binary=.false., arguments=arguments(1:2))

              close(occupationsUnit)     
           else

              !! Read density matrix
              !! Open file for wavefunction
              wfnFile = "lowdin.wfn"
              wfnUnit = 20
              open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

              arguments(2) = nameOfSpecies
              arguments(1) = "DENSITY"

              densityMatrix = &
                   Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfOrbitals,4), &
                   columns=int(numberOfOrbitals,4), binary=.true., arguments=arguments(1:2))

              close (wfnUnit)

           end if

           ! call Matrix_show(densityMatrix)

           !Define graph parameters
           numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
           plotDistance1=sqrt(sum((this%point2%values(:)-this%point1%values(:))**2))
           plotDistance2=sqrt(sum((this%point3%values(:)-this%point1%values(:))**2))
           step1%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
           step2%values(:)=(this%point3%values(:)-this%point1%values(:))/numberOfSteps
           
           write(auxstring,*) this%state
           title=trim(nameOfSpecies)//" state "//auxstring//" density" 

           val=0.0_8     
           maxValue=0.0_8
           minValue=0.0_8 

           !Write density grids according to the number of dimensions chosen
           if(this%dimensions.eq.3)then
              !!check if there is another density plot with the same same
              if( this%auxID .eq. 1) then
                 this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".3D.dens"
              else
                 this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".3D-"//trim(auxID)//".dens"
              end if
              x_title="x/a.u."
              y_title="y/a.u."
              z_title=""
              open(10,file=this%fileName(l),status='replace',action='write')
              write (10,"(A10,A20,A20,A20)") "#","X","Y","Density"
              do i=0,numberOfSteps
                 write (10,*) ""
                 do j=0,numberOfSteps
                    coordinate(:)=this%point1%values(:)+i*step1%values(:)+j*step2%values(:)
                    val=CalculateWaveFunction_getDensityAt( nameOfSpecies, coordinate, densityMatrix )  

                    write (10,"(T10,F20.8,F20.8,E20.8)") -plotDistance1*0.5+i*Vector_norm(step1), -plotDistance2*0.5+j*Vector_norm(step2),val 
                    if (val > maxValue) maxValue = val
                    if (val < minValue) minValue = val
                    ! print *, coordinate, val
                 end do
              end do

              !!large density values lead to bad looking plots
              if(maxValue .gt. 0.5) maxValue=0.5
              
              call OutputBuilder_make3DGraph( this%fileName(l), title, x_title, y_title, z_title, 0.0_8, maxValue)
              close(10)

           elseif(this%dimensions.eq.2) then
              !!check if there is another density plot with the same same
              if( this%auxID .eq. 1) then
                 this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".2D.dens"
              else
                 this%fileName(l)=trim(CONTROL_instance%INPUT_FILE)//trim(nameOfSpecies)//".2D-"//trim(auxID)//".dens"
              end if

              x_title="distance/a.u."
              y_title="density/a.u.^{-3}"
              open(10,file=this%fileName(l),status='replace',action='write')

              write (10,"(A10,A20,A20)") "#","X","Density"
              do i=0,numberOfSteps
                 coordinate(:)=this%point1%values(:)+i*step1%values(:)
                 val=CalculateWaveFunction_getDensityAt( nameOfSpecies, coordinate, densityMatrix )  

                 write (10,"(T10,F20.8,E20.8)") -plotDistance1*0.5+i*Vector_norm(step1),val 
                 ! print *, coordinate, val
              end do

              call OutputBuilder_make2DGraph( this%fileName(l), title, x_title, y_title)
              close(10)

           end if
        end if
     end do

     call Vector_Destructor(step1)
     call Vector_Destructor(step2)

   end subroutine OutputBuilder_getDensityPlot



   subroutine OutputBuilder_make2DGraph(fileName, title, x_title, y_title,&
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

     integer :: i,status
     character(20) :: charNumOfGraph
     character(20) :: auxXformat
     character(20) :: auxYformat
     character(20) :: auxXRange
     character(20) :: auxYRange

     integer :: auxNumOfGraphs

     auxXformat="%.2f"
     if(present(x_format)) auxXformat=trim(x_format)

     auxYformat="%.2f"
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
        write (10,"(A$)") 'plot '//trim(auxXRange)//trim(auxYRange)//' "'//trim(fileName)//'" using 1:2 w l title "" smooth csplines'
        do i=2, auxNumOfGraphs
           charNumOfGraph=String_convertIntegerToString(i+1)
           write (10,"(A$)") ', "'//trim(fileName)//'.dat"'//' using 1:'//trim(charNumOfGraph)//' w l  title "" smooth csplines'
        end do
        write (10,"(A)") ""
     else
        write (10,"(A)") 'plot '//trim(auxXRange)//trim(auxYRange)//' "'//trim(fileName)//'" w l title "" smooth csplines'
     end if
     write (10,"(A)") 'set output'
     close(10)

!     status= system("gnuplot "//trim(fileName)//".gnp")
     call system("gnuplot "//trim(fileName)//".gnp")

   end subroutine OutputBuilder_make2DGraph


   subroutine OutputBuilder_make3DGraph(fileName, title, x_title, y_title, z_title, minValue, maxValue)

     implicit none
     character(*) :: fileName
     character(*) :: title
     character(*) :: x_title
     character(*) :: y_title
     character(*) :: z_title
     real(8) :: minValue
     real(8) :: maxValue
     real(8) :: maxMinDiff
     
     integer :: status
     integer :: levels

     maxMinDiff=maxValue-minValue
     levels=10
     
     open ( 100,FILE=trim(fileName)//'.gnp', STATUS='REPLACE',ACTION='WRITE')
     write (100,"(A)") 'set term post eps enh color "Helvetica" 16 size 7cm,5cm'
     write (100,"(A)") 'set encoding iso_8859_1'
     write (100,"(A)") 'set output "'//trim(fileName)//'.eps"'
     
     write (100,"(A)") 'set table "'//trim(fileName)//'.table"'
     write (100,"(A)") 'splot "'//trim(fileName)//'" u 1:2:3'
     write (100,"(A)") 'unset table'

     write (100,"(A,I5)") 'levels=', levels
     write (100,"(A,I5)") 'numColors=', 5
     write (100,"(A,E20.8)") 'maxValue=', maxValue
     write (100,"(A,E20.8)") 'minValue=', minValue
     write (100,"(A)") 'step=(maxValue-minValue)/levels'
     write (100,"(A)") 'colorStep=(maxValue-minValue)/numColors'
     
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
     write (100,"(A)") 'set cbtics step'
     write (100,"(A)") 'set format cb "%3.1E"'
     
     write (100,"(A)") 'set palette defined (minValue "white",minValue+colorStep "blue",minValue+colorStep*2 "green",minValue+colorStep*3 "yellow",minValue+colorStep*4 "orange",maxValue "red")'
     write (100,"(A)") 'set grid front'

     
     write (100,"(A)") 'set format x "%.0f"'
     write (100,"(A)") 'set format y "%.0f"'
     write (100,"(A)") 'set xlabel "X (a.u.)"'
     write (100,"(A)") 'set ylabel "Y (a.u.)"'

     write (100,"(A)")  'plot "'//trim(fileName)//'.table" with image, "'//trim(fileName)//'.cont" w l lt -1 lw 1.5'
     
     
     ! write (100,"(A)") 'set output "'//trim(fileName)//'.eps"'
     ! write (100,"(A)") 'set xlabel "'//trim(x_title)//'"'
     ! write (100,"(A)") 'set ylabel "'//trim(y_title)//'"'
     ! write (100,"(A)") 'set zlabel "'//trim(z_title)//'"'
     ! write (100,"(A)") 'set pm3d '

     ! if (minValue < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. maxValue > -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     ! write (100,"(A)") 'set palette model RGB defined ('//String_convertRealToString(minValue)//&
     !      ' "violet", '//String_convertRealToString(minValue/5)//&
     !      ' "blue", '//String_convertRealToString(minValue/25)//&
     !      ' "green", '//String_convertRealToString(0.0_8) // &
     !      ' "white", '//String_convertRealToString(maxValue/25)//&
     !      ' "yellow", '//String_convertRealToString(maxValue/5)//&
     !      ' "orange", '//String_convertRealToString(maxValue)//' "red") '
     ! end if

     ! if (minValue > -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     ! write (100,"(A)") 'set palette model RGB defined ('//String_convertRealToString(maxValue/15625) // &
     !      ' "white", '//String_convertRealToString(maxValue/3125)//&
     !      ' "violet", '//String_convertRealToString(maxValue/625)//&
     !      ' "blue", '//String_convertRealToString(maxValue/125)//&
     !      ' "green", '//String_convertRealToString(maxValue/25)//&
     !      ' "yellow", '//String_convertRealToString(maxValue/5)//&
     !      ' "orange", '//String_convertRealToString(maxValue)//' "red") '
     ! end if

     ! if (maxValue < CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     ! write (100,"(A)") 'set palette model RGB defined ('//String_convertRealToString(minValue/15625) // &
     !      ' "white", '//String_convertRealToString(minValue/3125)//&
     !      ' "violet", '//String_convertRealToString(minValue/625)//&
     !      ' "blue", '//String_convertRealToString(minValue/125)//&
     !      ' "green", '//String_convertRealToString(minValue/25)//&
     !      ' "yellow", '//String_convertRealToString(minValue/5)//&
     !      ' "orange", '//String_convertRealToString(minValue)//' "red") '
     ! end if

     ! write (100,"(A)") 'set format "%.0f"'
     ! write (100,"(A)") 'set format cb "%.1f"'
     ! write (100,"(A)") 'set cbtics '//String_convertRealToString(maxValue/5)
     ! write (100,"(A)") 'set xyplane at 0'
     ! write (100,"(A)") 'set surface'
     ! write (100,"(A)") 'set border 4095'
     ! write (100,"(A)") 'set cntrparam cubicspline'
     ! write (100,"(A)") 'set cntrparam points 20'
     ! write (100,"(A)") 'set cntrparam levels 10'
     ! write (100,"(A)") 'set rmargin -1'
     ! write (100,"(A)") 'set lmargin -1'
     ! write (100,"(A)") 'set tmargin -1'
     ! write (100,"(A)") 'set bmargin -1'
     ! write (100,"(A)") 'unset ztics'
     ! write (100,"(A)") 'set multiplot title "'//trim(title)//'" layout 1,2'

     ! write (100,"(A)") 'set colorbox vertical'
     ! write (100,"(A)") 'set colorbox user origin 0.48,0.25 size 0.04,0.5'
     ! write (100,"(A)") 'set view 50,160'
     ! write (100,"(A)") 'splot "'//trim(fileName)//'" u 1:2:3 notitle w pm3d'

     ! write (100,"(A)") 'unset colorbox'
     ! write (100,"(A)") 'set view 0,0'
     ! write (100,"(A)") 'splot "'//trim(fileName)//'"  u 1:2:3 notitle w pm3d'

     ! write (100,"(A)") 'unset multiplot'

     close(100)

!     status= system("gnuplot "//trim(fileName)//".gnp")
     call system("gnuplot "//trim(fileName)//".gnp")

   end subroutine OutputBuilder_make3DGraph

end module OutputBuilder_

!    subroutine OutputBuilder_get2DDensityPlot(this)
!      implicit none
!      type(outputBuilder) :: this
!      character(50) :: outputID
!      character(50) :: orbitalNum

!      integer :: i
!      integer :: numberOfSteps
!      type(vector) :: step
!      real(8) :: val, val2
!      real(8) :: coordinate(3)

!      character(50) :: title
!      character(50) :: x_title
!      character(50) :: y_title

!      stop "trololo 2D"

! !      call Vector_Constructor(step, 3)

! !      this%fileName2=""
! !      numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
! !      step%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
! !      outputID=String_convertIntegerToString(this%outputID)

! !      select case( this%type )
! !      case ( "densityPlot") 


! !      case ( "orbitalPlot") 
! !         orbitalNum=String_convertIntegerToString(this%orbital)
! !         this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.orb"//trim(orbitalNum)
! !         open(10,file=this%fileName,status='replace',action='write')
! !         write (10,"(A10,A20,A20)") "#", "X","OrbitalValue"
! !         title=trim(this%specie)//" Orbital Number "//trim(orbitalNum) 
! !         y_title="orbitalValue/a.u.^{-3/2}"

! !      case ( "fukuiPlot") 
! !         this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.fkpos"
! !         this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.fkneg"

! !         open(10,file=this%fileName,status='replace',action='write')
! !         write (10,"(A10,A20,A20)") "#","X","PositiveFukuiValue"
! !         title=trim(this%specie)//" positive fukui" 
! !         y_title="density/a.u.^{-3}"

! !         open(11,file=this%fileName,status='replace',action='write')
! !         write (11,"(A10,A20,A20)") "#","X","NegativeFukuiValue"
     
! !      case default
! !         call OutputBuilder_exception(ERROR, "The output plot type you requested has not been implemented yet", "OutputBuilder_get3DPlot" )

! !      end select

! !      do i=0,numberOfSteps
! !         coordinate(:)=i*step%values(:)+this%point1%values(:)
! !         val=CalculateWaveFunction_getDensityAt( this%specie, coordinate )  
! !         select case( this%type )
! !         case ( "densityPlot") 
! !            val=CalculateWaveFunction_getDensityAt( this%specie, coordinate )  
! !         case ( "orbitalPlot") 
! !            val=CalculateWaveFunction_getOrbitalValueAt( this%specie, this%orbital, coordinate )  
! !         case ( "fukuiPlot") 
! ! !!           val=CalculateProperties_getFukuiAt( this%specie, "positive", coordinate )  
! ! !!           val2=CalculateProperties_getFukuiAt( this%specie, "negative", coordinate )  
! !         case default
! !         end select
! !         write (10,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val 
! !         if (this%type .eq. "fukuiPlot") write (11,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val2 
! !      end do

! !      close(10)

! !      call OutputBuilder_make2DGraph( this%fileName, title, x_title, y_title)
! ! !!     if (this%type .eq. "fukuiPlot") then
! ! !!        close(11)
! ! !!        title=trim(this%specie)//" negative fukui" 
! ! !!        call OutputBuilder_make2DGraph( this%fileName2, title, x_title, y_title)
! ! !!     end if
! !      call Vector_Destructor ( step)

!    end subroutine OutputBuilder_get2DDensityPlot

  ! subroutine OutputBuilder_getCube(this )
  !   implicit none
  !   type(output) :: this
  !   character(50) :: outputID
  !   real(8):: cubeSize
  !   character(50) :: orbitalNum

  !   integer :: i, j, k, n, w, natom
  !   integer :: atomicCharge
  !   integer :: specieID
  !   integer :: numberOfSteps
  !   real(8) :: step(3)
  !   real(8) :: lowerLimit(3)
  !   real(8), allocatable :: val(:), val2(:)
  !   real(8) :: coordinate(3)

  !   !Writes Gaussian Cube 
  !   this%fileName=""
  !   ! this%fileName2=""
  !   outputID=String_convertIntegerToString(this%outputID)
  !   specieID= MolecularSystem_getSpecieID( nameOfSpecie=this%specie)

  !   ! if (.not. allocated(CalculateProperties_instance%densityCube) ) call CalculateProperties_buildDensityCubesLimits(CalculateProperties_instance)

  !   ! if  (this%type .eq. "densityCube" .and. .not. CalculateProperties_instance%densityCube(specieID)%areValuesCalculated ) then
  !   !    call CalculateProperties_buildDensityCubes(CalculateProperties_instance)
  !   ! end if

  !   lowerLimit=this%point1
  !   numberOfSteps=CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
  !   step= this%cubeSize/numberOfSteps


  !   allocate (val (int(numberOfSteps(3))) )

  !   select case( this%type )
  !   case ( "densityCube") 
  !      this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".dens.cub"
  !      open(10,file=this%fileName,status='replace',action='write')

  !   ! case ( "orbitalCube") 
  !   !    orbitalNum=String_convertIntegerToString(this%orbital)
  !   !    this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".orb"//trim(orbitalNum)//".cub"
  !   !    open(10,file=this%fileName,status='replace',action='write')

  !   ! case ( "fukuiCube") 
  !   !    this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".fkpos.cub"
  !   !    this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".fkneg.cub"

  !   !    open(10,file=this%fileName,status='replace',action='write')
  !   !    open(11,file=this%fileName2,status='replace',action='write')

  !   case default
  !      call OutputBuilder_exception(ERROR, "The output cube type you requested has not been implemented yet", "OutputBuilder_getCube" )

  !   end select

  !   ! do n=1, size(MolecularSystem_instance%particlesPtr)
  !   !    if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-" .or. &
  !   !         trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-ALPHA" .and. &
  !   !         MolecularSystem_instance%particlesPtr(k)%isQuantum ) then
  !   !       natom = natom +1
  !   !    end if
  !   ! end do

  !   write (10,"(A)") "Gaussian Cube generated with Lowdin Software"
  !   write (10,"(A)") this%fileName
  !   write (10,"(I8,F20.8,F20.8,F20.8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3)
  !   write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(1)), step(1), 0.0, 0.0
  !   write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(2)), 0.0, step(2), 0.0
  !   write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(3)), 0.0, 0.0, step(3)
  !   ! do n=1, size(MolecularSystem_instance%particlesPtr)
  !   !    if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-" .or. &
  !   !         trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "E-ALPHA" .and. &
  !   !         MolecularSystem_instance%particlesPtr(n)%isQuantum ) then
  !   !       atomicCharge=-MolecularSystem_instance%particlesPtr(n)%totalCharge
  !   !       write (10, "(I8,F20.8,F20.8,F20.8,F20.8)") &
  !   !            atomicCharge, 0.0, MolecularSystem_instance%particlesPtr(n)%origin(1:3)
  !   !    end if
  !   ! end do

  !   ! if  (this%type .eq. "fukuiCube") then
  !   !    write (11,"(A)") "Gassian Cube generated with Lowdin Software"
  !   !    write (11,"(A)") this%fileName2
  !   !    write (11,"(I8,F20.8,F20.8,F20.8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3)
  !   !    write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(1)), step(1), 0.0, 0.0
  !   !    write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(2)), 0.0, step(2), 0.0
  !   !    write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(3)), 0.0, 0.0, step(3)
  !   !    do n=1, size(MolecularSystem_instance%particlesPtr)
  !   !       if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-" .or. &
  !   !            trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-ALPHA" .and. &
  !   !            MolecularSystem_instance%particlesPtr(n)%isQuantum ) then
  !   !          atomicCharge=-MolecularSystem_instance%particlesPtr(n)%totalCharge
  !   !          write (11, "(I8,F20.8,F20.8,F20.8,F20.8)") &
  !   !               atomicCharge, 0.0, MolecularSystem_instance%particlesPtr(n)%origin(1:3)
  !   !       end if
  !   !    end do
  !   ! end if
  
  !   do i=1,numberOfSteps
  !      coordinate(1)=lowerLimit(1)+(i-1)*step(1)
  !      do j=1, numberOfSteps
  !         coordinate(2)=lowerLimit(2)+(j-1)*step(2)
  !         do k=1, numberOfSteps
  !            coordinate(3)=lowerLimit(3)+(k-1)*step(3)
  !            select case (this%type)                   
  !            case ( "densityCube") 
  !               val(k)=CalculateProperties_instance%densityCube(specieID)%values(i,j,k)
  !               ! case ( "orbitalCube") 
  !               !     val(k)=MolecularSystem_getOrbitalValueAt( this%specie, this%orbital, coordinate )  
  !               !  case ( "fukuiCube") 
  !               !     val(k)=CalculateProperties_getFukuiAt( this%specie, "positive", coordinate )  
  !               !     val2(k)=CalculateProperties_getFukuiAt( this%specie, "negative", coordinate )  
  !            case default
  !            end select
  !         end do
  !         write(10,*) ( val(w) , w=1,numberOfSteps(3) )
  !         if (this%type .eq. "fukuiCube") write(11,*) ( val2(w) , w=1,numberOfSteps(3) )
  !      end do
  !   end do

  !   deallocate (val, val2)

  !   close(10)
  !   if  (this%type .eq. "fukuiCube" ) close(11)

  ! end subroutine OutputBuilder_getCube

