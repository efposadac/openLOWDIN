!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!> @brief This module handles all the information related with the basis set of each particle of each specie.
!! Only handles contracted-gaussian-functions-based basis set.
!! @author E. F. Posada (efposadac@unal.edu.co)
!! @version 1.0
module BasisSet_
  use ContractedGaussian_
  use CONTROL_
  implicit none
  
  !> The basis set is an array of contracted gaussians!
  type :: BasisSet
     type(contractedGaussian), allocatable :: contraction(:) !< Array of contractions for one particle.
     character(30) :: name
     real(8) :: origin(3)
     integer :: length
     integer :: ttype     
     integer :: contractionLength !< Dispuesta por razones de conveniencia
     integer :: numberOfPrimitives
  end type BasisSet
  
  public :: &
       BasisSet_load, &
       BasisSet_showInSimpleForm  

  private :: &
       BasisSet_exception
  
contains

  !> @brief Load basis sets form lowdin.bas file and deMon2K format
  !! @author E. F. Posada, 2013
  !! @version 1.0
  subroutine BasisSet_load(this, formatType, basisName, ParticleName, origin, unit)
    implicit none
    
    type(BasisSet) :: this
    character(*) :: formatType
    character(*), optional :: basisName
    character(*), optional :: particleName
    real(8), optional :: origin(3)
    integer, optional :: unit
    
    integer :: numberOfSpecies
    integer :: numberOfParticles
    integer :: numberOfShells
    integer :: i, j, k, l !< Iterators
    logical :: existFile, found
    character(20) :: token
    character(10) :: symbol
    character(30) :: particleSelected
    integer :: status
    
    !! Setting name
    if(present(basisName)) this%name = trim(basisName)
    if(present(particleName)) particleSelected = trim(particleName)
    if(present(origin)) this%origin = origin
    
    select case(trim(formatType))
       
    case("LOWDIN.BAS")
       
       read(unit,*) this%length
       
       allocate(this%contraction(this%length))
       
       do i = 1, this%length
          
          read(unit,*) this%contraction(i)%angularMoment, this%contraction(i)%length
          
          allocate(this%contraction(i)%orbitalExponents(this%contraction(i)%length))
          allocate(this%contraction(i)%contractionCoefficients(this%contraction(i)%length))
          
          !! Calculates the number of cartesian orbitals, by dimensionality
          select case(CONTROL_instance%DIMENSIONALITY)
             
          case(3)
             this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 )*( this%contraction(i)%angularMoment + 2_8 ) ) / 2_8
          case(2)
             this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 ) )
          case(1)
             this%contraction(i)%numCartesianOrbital = 1 
          case default
             call BasisSet_exception( ERROR, "Class object Basis set in load function",&
                  "This Dimensionality is not avaliable") 
          end select
               
          read(unit,*) this%contraction(i)%origin
          read(unit,*) this%contraction(i)%orbitalExponents
          read(unit,*) this%contraction(i)%contractionCoefficients
          
          !! Normalize
          allocate(this%contraction(i)%contNormalization(this%contraction(i)%numCartesianOrbital))
          allocate(this%contraction(i)%primNormalization(this%contraction(i)%length, &
               this%contraction(i)%length*this%contraction(i)%numCartesianOrbital))
          
          call ContractedGaussian_normalizePrimitive(this%contraction(i))
          call ContractedGaussian_normalizeContraction(this%contraction(i))
          
          !! DEBUG
          !! call ContractedGaussian_showInCompactForm(this%contraction(i))
          
       end do
       
    case("DEMON2K")          
       
       !! Open BasisSet file from library
       inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), exist = existFile)
       
       if(existFile) then
          
          !! Open File
          open(unit=30, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), status="old",form="formatted")
          rewind(30)
          
          found = .false.
          
          !! Open element and Find Element Basis set
          do while(found .eqv. .false.)
             
             read(30,*, iostat=status) token, symbol
             
             !! Some debug information in case of error!
             if (status > 0 ) then
                
                call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")
                
             end if
             
             if (status == -1 ) then
                
                call BasisSet_exception(ERROR, "The basisSet: "//trim(this%name)//" for: "//trim(particleSelected)//" was not found!","BasisSet module at Load function.")
                
             end if
             
             if(trim(token(1:2)) == "O-") then
                
                if(trim(symbol) == trim(particleSelected)) then
                   
                   found = .true.
                   
                end if
                
             end if
             
          end do
          
          !! Neglect any coment
          token = "#"
          do while(trim(token(1:1)) == "#")
             
             read(30,*) token
             
          end do
          
          !! Start reading basis set
          backspace(30)
          
          read(30,*, iostat=status) this%length
          
          !! Some debug information in case of error!
          if (status > 0 ) then
             
             call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")
             
          end if
          
          allocate(this%contraction(this%length))
          
          do i = 1, this%length
             
             read(30,*,iostat=status) this%contraction(i)%id, &
                  this%contraction(i)%angularMoment, &
                  this%contraction(i)%length
             
             !! Some debug information in case of error!
             if (status > 0 ) then
                
                call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")
                   
             end if
             
             allocate(this%contraction(i)%orbitalExponents(this%contraction(i)%length))
             allocate(this%contraction(i)%contractionCoefficients(this%contraction(i)%length))
             
             do j = 1, this%contraction(i)%length
                
                read(30,*,iostat=status) this%contraction(i)%orbitalExponents(j), &
                     this%contraction(i)%contractionCoefficients(j)
                
                !! Some debug information in case of error!
                if (status > 0 ) then
                   
                   call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")
                   
                end if
                
             end do
             
             !! Ajust and normalize contractions
             this%contraction(i)%origin = this%origin
             
             !! Calculates the number of cartesian orbitals, by dimensionality
             select case(CONTROL_instance%DIMENSIONALITY)
                
             case(3)
                this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 )*( this%contraction(i)%angularMoment + 2_8 ) ) / 2_8
             case(2)
                this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 ) )
             case(1)
                this%contraction(i)%numCartesianOrbital = 1 
             case default
                call BasisSet_exception( ERROR, "Class object Basis set in load function",&
                     "This Dimensionality is not avaliable") 
             end select
             
             !! Normalize
             allocate(this%contraction(i)%contNormalization(this%contraction(i)%numCartesianOrbital))
             allocate(this%contraction(i)%primNormalization(this%contraction(i)%length, &
                  this%contraction(i)%length*this%contraction(i)%numCartesianOrbital))
             
             this%contraction(i)%contNormalization = 1.0_8
             this%contraction(i)%primNormalization = 1.0_8

             call ContractedGaussian_normalizePrimitive(this%contraction(i))
             call ContractedGaussian_normalizeContraction(this%contraction(i))
             
             !! DEBUG
             !! call ContractedGaussian_showInCompactForm(this%contraction(i))
             
          end do
          
          close(30)
          
          !!DONE
          
       else
          
          call BasisSet_exception(ERROR, "The basisSet file: "//trim(basisName)//" was not found!","BasisSet module at Load function.")
             
       end if
       
    end select
    
  end subroutine BasisSet_load

  !>
  !! @brief Shows information of the object
  !! @param this: basis set
  subroutine BasisSet_showInCompactForm( this, nameOfOwner )
    implicit none
    
    type(BasisSet) , intent(in) :: this
    character(*) ::nameOfOwner
    integer ::  i, from, to
    
    print *,""
    write(6,"(T5,A10,A11,A15)") trim(nameOfOwner)," BASIS SET: ", trim(this%name)
    
    do i =1, this%length
       call ContractedGaussian_showInCompactForm(this%contraction(i) )
    end do

  end subroutine BasisSet_showInCompactForm

  !<
  !! Define el destructor para clase
  !!
  !! @param thisPtr Funcion base
  !>
  subroutine BasisSet_showInSimpleForm( this, nameOfOwner, unidOfOutput )
  	implicit none
  	type(BasisSet) , intent(in) :: this
  	character(*) ::nameOfOwner
  	integer :: unidOfOutput
  
  	integer ::  i
  	
  	do i =1, this%length
  		call ContractedGaussian_showInSimpleForm( this%contraction(i),unidOfOutput )
  	end do
  
  end subroutine BasisSet_showInSimpleForm


  !>
  !! @brief Saves the basisSet structure to file.
  !! @param this Basis-set object
  !! @author E. F. Posada, 2013
  subroutine BasisSet_saveToFile(this, unit)
    implicit none
    
    type(BasisSet), intent(in) :: this
    integer :: unit
    
    integer :: i
    
    write(unit,*) this%name
    write(unit,*) this%origin
    write(unit,*) this%length
    write(unit,*) this%ttype
    write(unit,*) this%contractionLength
    write(unit,*) this%numberOfPrimitives

    do i = 1, size(this%contraction)
       call ContractedGaussian_saveToFile(this%contraction(i), unit)
    end do
    
  end subroutine BasisSet_saveToFile

  !>
  !! @brief Loads the basisSet structure from file.
  !! @param this Basis-set object
  !! @author E. F. Posada, 2013
  subroutine BasisSet_loadFromFile(this, unit)
    implicit none
    
    type(BasisSet) :: this
    integer :: unit
    
    integer :: i
    
    read(unit,*) this%name
    read(unit,*) this%origin
    read(unit,*) this%length
    read(unit,*) this%ttype
    read(unit,*) this%contractionLength
    read(unit,*) this%numberOfPrimitives

    allocate(this%contraction(this%length))

    do i = 1, size(this%contraction)
       call ContractedGaussian_loadFromFile(this%contraction(i), unit)
    end do
    
  end subroutine BasisSet_loadFromFile
 
  !>
  !! @brief handle exceptions
  subroutine BasisSet_exception( typeMessage, description, debugDescription)
    implicit none
    integer:: typeMessage
    character(*):: description
    character(*):: debugDescription
    
    type(Exception):: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )
    
  end subroutine BasisSet_exception
  
end module BasisSet_
