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
!!		R. Flores (roberto.floresmoreno.qt@gmail.com)
!!
!!	Contributors:
!!
!!		Todos los derechos reservados, 2011
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
     character(50) :: specie
     character(50) :: fileName
     character(50) :: fileName2
     integer :: orbital
     integer :: dimensions
     integer :: outputID
     real(8) :: cubeSize
     type(vector) :: point1
     type(vector) :: point2
     type(vector) :: point3
     logical :: isInstanced
  end type OutputBuilder

  public :: &
       OutputBuilder_constructor, &
       OutputBuilder_destructor, &
       OutputBuilder_show, &
       OutputBuilder_writeMoldenFile, &
       OutputBuilder_VecGamessFile, &
       OutputBuilder_generateAIMFiles, &
       OutputBuilder_generateExtendedWfnFile, &
       OutputBuilder_buildOutput, &
       OutputBuilder_make2DGraph, &
       OutputBuilder_make3DGraph
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
  !<
  subroutine OutputBuilder_constructor(this, ID, type ,specie, orbital, dimensions, cubeSize, point1, point2, point3  )
    character(*) :: type
    integer :: ID
    character(*) :: specie
    integer :: orbital
    integer :: dimensions
    real(8) :: cubeSize
    type(Vector) :: point1
    type(Vector) :: point2
    type(Vector) :: point3
    type(OutputBuilder) :: this

    this%type=type
    this%outputID=ID
    this%specie=specie
    this%orbital=orbital
    this%dimensions=dimensions
    this%cubeSize=cubeSize
    call Vector_copyConstructor(this%point1, point1)
    call Vector_copyConstructor(this%point2, point2)
    call Vector_copyConstructor(this%point3, point3)

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

    print *, "--------------------------------------------------------"
    print *, "Output Number: ", this%outputID
    print *, "FileName: ", this%fileName
    if (this%filename2 /= "") print *, "FileName 2: ", this%fileName2
    print *, this%type
    if (this%specie /= "") write (6,"(A30,A10)") "for specie: ", this%specie
    if (this%orbital /= 0) write (6,"(A30,I10)") "for orbital: ", this%orbital
    if (this%dimensions /= 0) write (6,"(A30,I2)") "number of dimensions: ", this%dimensions
    if (this%cubeSize /= 0.0_8) write (6,"(A30,F15.12)") "cube size in a.u.: ", this%cubeSize
    if (this%dimensions >= 1) write (6,"(A30,F15.12,F15.12,F15.12)") "Point 1: ", this%point1%values(1), this%point1%values(2), this%point1%values(3)
    if (this%dimensions >= 2) write (6,"(A30,F15.12,F15.12,F15.12)") "Point 2: ", this%point2%values(1), this%point2%values(2), this%point2%values(3)
    if (this%dimensions >= 3) write (6,"(A30,F15.12,F15.12,F15.12)") "Point 3: ", this%point3%values(1), this%point3%values(2), this%point3%values(3)
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
        if (this%dimensions == 2) call OutputBuilder_get2DPlot(this)
        if (this%dimensions == 3) call OutputBuilder_get3DPlot(this)

!    case ( "densityCube") 
!       call OutputBuilder_getCube(this)
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

    auxString="speciesName"

    wfnFile = "lowdin.wfn"
    wfnUnit = 20
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

!! Open file for wavefunction                                                                                     
        open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")


        do l=1,MolecularSystem_getNumberOfQuantumSpecies()

	   totalNumberOfParticles = 0

           auxString=MolecularSystem_getNameOfSpecie( l )
           specieID = MolecularSystem_getSpecieID(auxString)
           this%fileName=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".molden"
           open(10,file=this%fileName,status='replace',action='write')
           write(10,"(A)") "[Molden Format]"
           write(10,"(A)") "[Atoms] AU"
           auxMatrix%values=0.0
           j=0
           do i=1, size(MolecularSystem_instance%species(l)%particles)
!              if ( trim(MolecularSystem_instance%species(l)%particles(i)%symbol) == trim(auxString) ) then
                 j=j+1
                 origin = MolecularSystem_instance%species(l)%particles(j)%origin 
                 auxMatrix%values(j,:)=origin
                 symbol=MolecularSystem_instance%species(l)%particles(j)%nickname

                 if(scan(symbol,"_") /=0) symbol=symbol(1:scan(symbol,"_")-1)
                 if(scan(symbol,"[") /=0) symbol=symbol(scan(symbol,"[")+1:scan(symbol,"]")-1)

!                 if ( CONTROL_instance%UNITS=="ANGSTROMS") origin = origin * AMSTRONG

		totalNumberOfParticles = totalNumberOfParticles + 1
#ifdef intel
                 write (10,"(A,I,I,<3>F15.8)") trim(symbol), j,&
                      int(abs(MolecularSystem_instance%species(l)%particles(j)%totalCharge)), origin(1), origin(2), origin(3)
#else

                 write (10,"(A,I8,I8,3F15.8)") trim(symbol), j,&
                      int(abs(MolecularSystem_instance%species(l)%particles(j)%charge)), origin(1), origin(2), origin(3)
#endif

!                end if
              end do

           m=j
           do k=1,size(localizationOfCenters%values,dim=1)

              wasPress=.false.
              do i=1,j
                 if( 	abs( auxMatrix%values(i,1) - localizationOfCenters%values(k,1)) < 1.0D-9 .and. &
                      abs( auxMatrix%values(i,2) - localizationOfCenters%values(k,2)) < 1.0D-9 .and. &
                      abs( auxMatrix%values(i,3) - localizationOfCenters%values(k,3)) < 1.0D-9  ) then
                    wasPress=.true.
                 end if
              end do

              if( .not.wasPress) then
                 m=m+1

		totalNumberOfParticles = totalNumberOfParticles + 1
                 origin=localizationOfCenters%values(k,:)
!                 if ( CONTROL_instance%UNITS=="ANGSTROMS") origin = origin * AMSTRONG
                 symbol=labels(k)
                 if(scan(symbol,"_") /=0) symbol=symbol(1:scan(symbol,"_")-1)
#ifdef intel
                 write (10,"(A,I,I,<3>F15.8)") trim(symbol), m,int(abs(charges(k))), origin(1), origin(2), origin(3)
#else
                 write (10,"(A,I8,I8,3F15.8)") trim(symbol), m,int(abs(charges(k))), origin(1), origin(2), origin(3)
#endif
              end if

           end do

!          print *, "totalNumberOfParticles ", totalNumberOfParticles
!         print *, "particles for specie", size(MolecularSystem_instance%species(l)%particles)

           write(10,"(A)") "[GTO]"
           j=0
           do i=1,size(MolecularSystem_instance%species(l)%particles)

!              if (	trim(MolecularSystem_instance%species(l)%particles(i)%symbol) == trim(auxString) ) then
                 j=j+1

                 write(10,"(I3,I2)") j,0
                 call BasisSet_showInSimpleForm( MolecularSystem_instance%species(l)%particles(i)%basis,&
                      trim(MolecularSystem_instance%species(l)%particles(i)%nickname),10 )
                 write(10,*) ""

		if ( totalNumberOfParticles > size(MolecularSystem_instance%species(l)%particles) ) then

			do n = 1, ( totalNumberOfParticles - size(MolecularSystem_instance%species(l)%particles) )
				write(10,"(I3,I2)") j+n,0
				write(10,"(A,I1,F5.2)") "s  ",1,1.00
				write(10,"(ES19.10,ES19.10)") 0.00,0.00
				write(10,*) ""
			end do 
		end if
!              end if

           end do

           write(10,"(A)") "[MO]"

           specieID = int( MolecularSystem_getSpecieID(nameOfSpecie = trim(auxString)) )
           numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
           arguments(2) = MolecularSystem_getNameOfSpecie(specieID)
           occupationTotal=MolecularSystem_getOcupationNumber( specieID )
           occupation =1.0/MolecularSystem_instance%species(specieID)%particlesFraction


           arguments(1) = "COEFFICIENTS"
           coefficientsOfcombination = &
                Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

           arguments(1) = "ORBITALS"
           call Vector_getFromFile( elementsNum = numberOfContractions, &
                unit = wfnUnit, binary = .true., arguments = arguments(1:2), &
                output = energyOfMolecularOrbital )

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

             	if ( shellcode == "Fxxx" ) then 
		    auxcounter = counter
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+1 , auxcounter+6)
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+2 , auxcounter+9)
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+4 , auxcounter+6)
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+5 , auxcounter+9)
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+6 , auxcounter+9)
	            call Matrix_swapRows(  coefficientsOfCombination, auxcounter+7 , auxcounter+8)

                end if

	    end do

           do j=1,size(energyOfMolecularOrbital%values)
              write (10,"(A5,ES15.5)") "Ene= ",energyOfMolecularOrbital%values(j)

              write (10,"(A11)") "Spin= Alpha"

              if ( j <= occupationTotal) then 
                 write (10,"(A,F7.4)") "Occup= ",occupation
              else
                 write (10,"(A)") "Occup=0.0000"
              end if
              do k=1,size(coefficientsOfCombination%values,dim=1)
                 write(10,"(I4,F15.6)") k,coefficientsOfCombination%values(k,j)
              end do

		if ( totalNumberOfParticles > size(MolecularSystem_instance%species(l)%particles) ) then
			do n = 1, ( totalNumberOfParticles - size(MolecularSystem_instance%species(l)%particles) )
                 		write(10,"(I4,F15.6)") k-1+n,0.0_8
			end do
		end if

           end do

           close(10)
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

    auxString="speciesName"

    wfnFile = "lowdin.wfn"
    wfnUnit = 20

!! Open file for wavefunction                                                                                     
        open(unit=wfnUnit, file=trim(wfnFile), status="old", form="unformatted")

        do l=1,MolecularSystem_getNumberOfQuantumSpecies()

	   auxString=MolecularSystem_getNameOfSpecie( l )

           this%fileName=trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//".vec"

           open(29,file=this%fileName,status='replace',action='write')
           
           specieID = int( MolecularSystem_getSpecieID(nameOfSpecie = trim(auxString)) )
           numberOfContractions = MolecularSystem_getTotalNumberOfContractions(specieID)
           arguments(2) = MolecularSystem_getNameOfSpecie(specieID)
 
           arguments(1) = "COEFFICIENTS"
           coefficientsOfcombination = &
                Matrix_getFromFile(unit=wfnUnit, rows= int(numberOfContractions,4), &
                columns= int(numberOfContractions,4), binary=.true., arguments=arguments(1:2))

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
                write (29,"(I3,I3)",advance='no') i,j
                do m=1,numberOfContractions
                   if (mod(m,5)==0) then
                      write (29,"(ES15.8)") coefficientsOfCombination%values(m,i)
                      j=j+1
                      if (m<numberOfContractions) then
                         write (29,"(I3,I3)",advance='no') i,j
                      end if
                   else
                      write (29,"(ES15.8)",advance='no') coefficientsOfCombination%values(m,i)
                   end if
                end do
                 write (29, "(A)", advance='yes')" "
                if (m<numberOfContractions) then
                   !   write (29,"(A)", advance='no')" "
                    write (29,"(A)")" "
                end if
                
             end do

           close(29)
        end do

!        call Matrix_destructor( localizationOfCenters )
!        call Matrix_destructor( auxMatrix )
!        deallocate(labels)

!     end if

      call OutputBuilder_exception(WARNING, "The order of the coefficients only works until F orbitals", "OutputBuilder_VecGamessFile" )
        
  end subroutine OutputBuilder_VecGamessFile

  !!!!!!!!!!END GAMESS .VEC FILE


  
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
    end do

    !! Just for printing information 
    this%fileName = trim(CONTROL_instance%INPUT_FILE)//trim(auxString)//extension//" and .molden"
 
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

!!   subroutine OutputBuilder_getCube(this )
!!     implicit none
!!     type(output) :: this
!!     character(50) :: outputID
!!     real(8):: cubeSize
!!     character(50) :: orbitalNum
!!
!!     integer :: i, j, k, n, w, natom
!!     integer :: atomicCharge
!!     integer :: specieID
!!     real(8) :: numberOfSteps(3)
!!     real(8) :: step(3)
!!     real(8) :: lowerLimit(3)
!!     real(8), allocatable :: val(:), val2(:)
!!     real(8) :: coordinate(3)
!!
!!     !Writes Gaussian Cube 
!!     this%fileName=""
!!     this%fileName2=""
!!     outputID=String_convertIntegerToString(this%outputID)
!!     specieID= MolecularSystem_getSpecieID( nameOfSpecie=this%specie)
!!
!!     if (.not. allocated(CalculateProperties_instance%densityCube) ) call CalculateProperties_buildDensityCubesLimits(CalculateProperties_instance)
!!
!!     if  (this%type .eq. "densityCube" .and. .not. CalculateProperties_instance%densityCube(specieID)%areValuesCalculated ) then
!!        call CalculateProperties_buildDensityCubes(CalculateProperties_instance)
!!     end if
!!
!!     lowerLimit=CalculateProperties_instance%densityCube(specieID)%lowerLimit%values
!!     numberOfSteps=CalculateProperties_instance%densityCube(specieID)%numberOfPoints%values
!!     step=CalculateProperties_instance%densityCube(specieID)%stepSize%values
!!
!!     allocate (val (int(numberOfSteps(3))) , val2(int(numberOfSteps(3))))
!!
!!     select case( this%type )
!!     case ( "densityCube") 
!!        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".dens.cub"
!!        open(10,file=this%fileName,status='replace',action='write')
!!
!!!!     case ( "orbitalCube") 
!!!!        orbitalNum=String_convertIntegerToString(this%orbital)
!!!!        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".orb"//trim(orbitalNum)//".cub"
!!!!        open(10,file=this%fileName,status='replace',action='write')
!!!!
!!!!     case ( "fukuiCube") 
!!!!        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".fkpos.cub"
!!!!        this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".fkneg.cub"
!!
!!!!        open(10,file=this%fileName,status='replace',action='write')
!!!!        open(11,file=this%fileName2,status='replace',action='write')
!!
!!     case default
!!        call OutputBuilder_exception(ERROR, "The output cube type you requested has not been implemented yet", "OutputBuilder_getCube" )
!!
!!     end select
!!
!!     do n=1, size(MolecularSystem_instance%particlesPtr)
!!        if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-" .or. &
!!             trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-ALPHA" .and. &
!!             MolecularSystem_instance%particlesPtr(k)%isQuantum ) then
!!           natom = natom +1
!!        end if
!!     end do
!!
!!     write (10,"(A)") "Gaussian Cube generated with Lowdin Software"
!!     write (10,"(A)") this%fileName
!!     write (10,"(I8,F20.8,F20.8,F20.8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3)
!!     write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(1)), step(1), 0.0, 0.0
!!     write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(2)), 0.0, step(2), 0.0
!!     write (10,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(3)), 0.0, 0.0, step(3)
!!     do n=1, size(MolecularSystem_instance%particlesPtr)
!!        if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-" .or. &
!!             trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-ALPHA" .and. &
!!             MolecularSystem_instance%particlesPtr(n)%isQuantum ) then
!!           atomicCharge=-MolecularSystem_instance%particlesPtr(n)%totalCharge
!!           write (10, "(I8,F20.8,F20.8,F20.8,F20.8)") &
!!                atomicCharge, 0.0, MolecularSystem_instance%particlesPtr(n)%origin(1:3)
!!        end if
!!     end do
!!
!!     if  (this%type .eq. "fukuiCube") then
!!        write (11,"(A)") "Gassian Cube generated with Lowdin Software"
!!        write (11,"(A)") this%fileName2
!!        write (11,"(I8,F20.8,F20.8,F20.8)") natom, lowerLimit(1), lowerLimit(2), lowerLimit(3)
!!        write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(1)), step(1), 0.0, 0.0
!!        write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(2)), 0.0, step(2), 0.0
!!        write (11,"(I8,F20.8,F20.8,F20.8)") int(numberOfSteps(3)), 0.0, 0.0, step(3)
!!        do n=1, size(MolecularSystem_instance%particlesPtr)
!!           if ( trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-" .or. &
!!                trim(MolecularSystem_instance%particlesPtr(n)%symbol) == "e-ALPHA" .and. &
!!                MolecularSystem_instance%particlesPtr(n)%isQuantum ) then
!!              atomicCharge=-MolecularSystem_instance%particlesPtr(n)%totalCharge
!!              write (11, "(I8,F20.8,F20.8,F20.8,F20.8)") &
!!                   atomicCharge, 0.0, MolecularSystem_instance%particlesPtr(n)%origin(1:3)
!!           end if
!!        end do
!!     end if
!!   
!!     do i=1,numberOfSteps(1)
!!        coordinate(1)=lowerLimit(1)+(i-1)*step(1)
!!        do j=1, numberOfSteps(2)
!!           coordinate(2)=lowerLimit(2)+(j-1)*step(2)
!!           do k=1, numberOfSteps(3)
!!              coordinate(3)=lowerLimit(3)+(k-1)*step(3)
!!              select case (this%type)                   
!!              case ( "densityCube") 
!!                 val(k)=CalculateProperties_instance%densityCube(specieID)%values(i,j,k)
!!!!             case ( "orbitalCube") 
!!!!                 val(k)=MolecularSystem_getOrbitalValueAt( this%specie, this%orbital, coordinate )  
!!!!              case ( "fukuiCube") 
!!!!                 val(k)=CalculateProperties_getFukuiAt( this%specie, "positive", coordinate )  
!!!!                 val2(k)=CalculateProperties_getFukuiAt( this%specie, "negative", coordinate )  
!!              case default
!!              end select
!!           end do
!!           write(10,*) ( val(w) , w=1,numberOfSteps(3) )
!!           if (this%type .eq. "fukuiCube") write(11,*) ( val2(w) , w=1,numberOfSteps(3) )
!!        end do
!!     end do
!!
!!     deallocate (val, val2)
!!
!!     close(10)
!!     if  (this%type .eq. "fukuiCube" ) close(11)
!!
!!   end subroutine OutputBuilder_getCube

   subroutine OutputBuilder_get3DPlot(this)
     type(OutputBuilder) :: this
     character(50) :: outputID
     character(50) :: orbitalNum

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

     this%fileName2=""
     numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
     step1%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
     step2%values(:)=(this%point3%values(:)-this%point1%values(:))/numberOfSteps
     outputID=String_convertIntegerToString(this%outputID)

     x_title="x/a.u."
     y_title="y/a.u."
     z_title=""

     select case( this%type )
     case ( "densityPlot") 
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".3D.dens"
        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20,A20)") "#","X","Y","Density"
        title=trim(this%specie)//" density" 

     case ( "orbitalPlot") 
        orbitalNum=String_convertIntegerToString(this%orbital)
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".3D.orb"//trim(orbitalNum)
        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20,A20)") "#", "X","Y","OrbitalValue"
        title=trim(this%specie)//" Orbital Number: "//trim(orbitalNum) 

     case ( "fukuiPlot") 
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".3D.fkpos"
        this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".3D.fkneg"
        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20,A20)") "#", "X","Y","PositiveFukuiValue"
        title=trim(this%specie)//" Positive Fukui"
        open(11,file=this%fileName2,status='replace',action='write')
        write (11,"(A10,A20,A20,A20)") "#", "X","Y","NegativeFukuiValue"
        title2=trim(this%specie)//" Negative Fukui"

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
           case ( "densityPlot") 
              val=CalculateWaveFunction_getDensityAt( this%specie, coordinate )  
           case ( "orbitalPlot") 
              val=CalculateWaveFunction_getOrbitalValueAt( this%specie, this%orbital, coordinate )  
           case ( "fukuiPlot") 
!!              val=CalculateProperties_getFukuiAt( this%specie, "positive", coordinate )  
!!              val2=CalculateProperties_getFukuiAt( this%specie, "negative", coordinate )  
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

     call OutputBuilder_make3DGraph( this%fileName, title, x_title, y_title, z_title, minValue, maxValue)
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
     character(50) :: outputID
     character(50) :: orbitalNum

     integer :: i
     integer :: numberOfSteps
     type(vector) :: step
     real(8) :: val, val2
     real(8) :: coordinate(3)

     character(50) :: title
     character(50) :: x_title
     character(50) :: y_title

     call Vector_Constructor(step, 3)

     this%fileName2=""
     numberOfSteps= CONTROL_instance%NUMBER_OF_POINTS_PER_DIMENSION
     step%values(:)=(this%point2%values(:)-this%point1%values(:))/numberOfSteps
     outputID=String_convertIntegerToString(this%outputID)

     x_title="distance/a.u."
     select case( this%type )
     case ( "densityPlot") 
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.dens"
        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20)") "#","X","Density"
        title=trim(this%specie)//" density" 
        y_title="density/a.u.^{-3}"

     case ( "orbitalPlot") 
        orbitalNum=String_convertIntegerToString(this%orbital)
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.orb"//trim(orbitalNum)
        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20)") "#", "X","OrbitalValue"
        title=trim(this%specie)//" Orbital Number "//trim(orbitalNum) 
        y_title="orbitalValue/a.u.^{-3/2}"

     case ( "fukuiPlot") 
        this%fileName=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.fkpos"
        this%fileName2=trim(CONTROL_instance%INPUT_FILE)//"out"//trim(outputID)//"."//trim(this%specie)//".2D.fkneg"

        open(10,file=this%fileName,status='replace',action='write')
        write (10,"(A10,A20,A20)") "#","X","PositiveFukuiValue"
        title=trim(this%specie)//" positive fukui" 
        y_title="density/a.u.^{-3}"

        open(11,file=this%fileName,status='replace',action='write')
        write (11,"(A10,A20,A20)") "#","X","NegativeFukuiValue"
     
     case default
        call OutputBuilder_exception(ERROR, "The output plot type you requested has not been implemented yet", "OutputBuilder_get3DPlot" )

     end select

     do i=0,numberOfSteps
        coordinate(:)=i*step%values(:)+this%point1%values(:)
        val=CalculateWaveFunction_getDensityAt( this%specie, coordinate )  
        select case( this%type )
        case ( "densityPlot") 
           val=CalculateWaveFunction_getDensityAt( this%specie, coordinate )  
        case ( "orbitalPlot") 
           val=CalculateWaveFunction_getOrbitalValueAt( this%specie, this%orbital, coordinate )  
        case ( "fukuiPlot") 
!!           val=CalculateProperties_getFukuiAt( this%specie, "positive", coordinate )  
!!           val2=CalculateProperties_getFukuiAt( this%specie, "negative", coordinate )  
        case default
        end select
        write (10,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val 
        if (this%type .eq. "fukuiPlot") write (11,"(T10,F20.8,F20.8)")  i*Vector_norm(step),val2 
     end do

     close(10)

     call OutputBuilder_make2DGraph( this%fileName, title, x_title, y_title)
!!     if (this%type .eq. "fukuiPlot") then
!!        close(11)
!!        title=trim(this%specie)//" negative fukui" 
!!        call OutputBuilder_make2DGraph( this%fileName2, title, x_title, y_title)
!!     end if
     call Vector_Destructor ( step)

   end subroutine OutputBuilder_get2DPlot



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

     integer :: status


     open ( 10,FILE=trim(fileName)//".gnp", STATUS='REPLACE',ACTION='WRITE')
     write (10,"(A)") 'set terminal postscript enhanced eps size 7,3.5 "Helvetica" 25'
     write (10,"(A)") 'set output "'//trim(fileName)//'.eps"'
     write (10,"(A)") 'set encoding iso_8859_1'
     write (10,"(A)") 'set xlabel "'//trim(x_title)//'"'
     write (10,"(A)") 'set ylabel "'//trim(y_title)//'"'
     write (10,"(A)") 'set zlabel "'//trim(z_title)//'"'
     write (10,"(A)") 'set pm3d '

     if (minValue < CONTROL_instance%DOUBLE_ZERO_THRESHOLD .and. maxValue > -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     write (10,"(A)") 'set palette model RGB defined ('//String_convertRealToString(minValue)//&
          ' "violet", '//String_convertRealToString(minValue/5)//&
          ' "blue", '//String_convertRealToString(minValue/25)//&
          ' "green", '//String_convertRealToString(0.0_8) // &
          ' "white", '//String_convertRealToString(maxValue/25)//&
          ' "yellow", '//String_convertRealToString(maxValue/5)//&
          ' "orange", '//String_convertRealToString(maxValue)//' "red") '
     end if

     if (minValue > -CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     write (10,"(A)") 'set palette model RGB defined ('//String_convertRealToString(maxValue/15625) // &
          ' "white", '//String_convertRealToString(maxValue/3125)//&
          ' "violet", '//String_convertRealToString(maxValue/625)//&
          ' "blue", '//String_convertRealToString(maxValue/125)//&
          ' "green", '//String_convertRealToString(maxValue/25)//&
          ' "yellow", '//String_convertRealToString(maxValue/5)//&
          ' "orange", '//String_convertRealToString(maxValue)//' "red") '
     end if

     if (maxValue < CONTROL_instance%DOUBLE_ZERO_THRESHOLD) then
     write (10,"(A)") 'set palette model RGB defined ('//String_convertRealToString(minValue/15625) // &
          ' "white", '//String_convertRealToString(minValue/3125)//&
          ' "violet", '//String_convertRealToString(minValue/625)//&
          ' "blue", '//String_convertRealToString(minValue/125)//&
          ' "green", '//String_convertRealToString(minValue/25)//&
          ' "yellow", '//String_convertRealToString(minValue/5)//&
          ' "orange", '//String_convertRealToString(minValue)//' "red") '
     end if

     write (10,"(A)") 'set format "%.0f"'
     write (10,"(A)") 'set format cb "%.1f"'
     write (10,"(A)") 'set cbtics '//String_convertRealToString(maxValue/5)
     write (10,"(A)") 'set xyplane at 0'
     write (10,"(A)") 'set surface'
     write (10,"(A)") 'set border 4095'
     write (10,"(A)") 'set cntrparam cubicspline'
     write (10,"(A)") 'set cntrparam points 20'
     write (10,"(A)") 'set cntrparam levels 10'
     write (10,"(A)") 'set rmargin -1'
     write (10,"(A)") 'set lmargin -1'
     write (10,"(A)") 'set tmargin -1'
     write (10,"(A)") 'set bmargin -1'
     write (10,"(A)") 'unset ztics'
     write (10,"(A)") 'set multiplot title "'//trim(title)//'" layout 1,2'

     write (10,"(A)") 'set colorbox vertical'
     write (10,"(A)") 'set colorbox user origin 0.48,0.25 size 0.04,0.5'
     write (10,"(A)") 'set view 50,160'
     write (10,"(A)") 'splot "'//trim(fileName)//'" u 1:2:3 notitle w pm3d'

     write (10,"(A)") 'unset colorbox'
     write (10,"(A)") 'set view 0,0'
     write (10,"(A)") 'splot "'//trim(fileName)//'"  u 1:2:3 notitle w pm3d'

     write (10,"(A)") 'unset multiplot'

     close(10)

!     status= system("gnuplot "//trim(fileName)//".gnp")
     call system("gnuplot "//trim(fileName)//".gnp")

   end subroutine OutputBuilder_make3DGraph


! set xtics -3,2,3 nomirror tc lt 0
! set ytics -3,2,3 nomirror tc lt 0
! set mxtics 0.5
! set mytics 0.5
! set xlabel "x /a.u"
! set ylabel "y /a.u"
! set xrange [-3.5:3.5]
! set yrange [-3.5:3.5]

end module OutputBuilder_
