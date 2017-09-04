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

!> @brief This module handles the information related with the effective core potentials.
!! @author I. Ortiz-Verano (ieortizv@unal.edu.co)
!! @version 1.0

module EffectiveCorePotentials_
  use ContractedEcpGaussian_
  use ContractedGaussian_
  use CONTROL_
!  use InputManager_ 
  implicit none

  type :: EffectiveCorePotentials
     !     type(contractedEcpGaussian),allocatable :: contraction(:)
     type(contractedEcpGaussian),allocatable :: contraction(:)
     !     type(contractedEcpGaussian),allocatable :: nkParameters(:)   !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     !     type(contractedEcpGaussian),allocatable :: zetaParameters(:) !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     !     type(contractedEcpGaussian),allocatable :: dkParameters(:)   !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     character(30) :: name
!     character(30) :: baseName 
     real(8) :: origin(3)
     integer :: ttype   !!!!!!!!!!!?????????????
     integer :: numberOfCoreElectrons
     integer :: length
     integer :: contractionLength
     integer :: numberOfPrimitives
  end type EffectiveCorePotentials

  public :: &
       EffectiveCorePotentials_load, &
       EffectiveCorePotentials_getNumberOfCoreElectrons
  !      EffectiveCorePotentials_loadECP, &
  !      EffectiveCorePotentials_showInSimpleForm  
  
  private :: &
       EffectiveCorePotentials_exception

contains


  !> @brief Load number of core electrons from ECP contined in basis set file in deMon2K format. If an atom hasn't ECP (like lighter elements), no information about this element will appear
  !! @author I. Ortiz-Verano, 2017
  !! @version 1.0
  function EffectiveCorePotentials_getNumberOfCoreElectrons(symbol, basisName) result(output)
    implicit none
    character(*) :: symbol
    character(*) :: basisName
    integer :: output

    character(10) :: token, ssymbol, nelec
    logical :: existFile
    logical :: foundECP, foundElement
    integer :: status
    integer :: numberOfCoreElectrons

    ! Open effective core potentials from basis set file at library
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), exist = existFile)

    if(existFile) then

       !! Open File
       open(unit=30, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), status="old",form="formatted")
       rewind(30)

       foundECP = .false.

       !! Find ECP group in basis set file
       do while(foundECP .eqv. .false.)

          read(30,*, iostat=status) token

          if (status > 0 ) then   

             call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(basisName)//" Please check that this file contains pseudopotentials!","EffectiveCorePotentials module at getNumberOfCoreElectrons function.")

          end if

          if (status == -1 ) then

             call EffectiveCorePotentials_exception(ERROR, "The basisSet: "//trim(basisName)//" for: "//trim(symbol)//" was not found!","BasisSet module at getNumberOfCoreElectrons function.")

          end if

          if(trim(token(1:3)) == "ECP") then

             foundECP = .true.

          end if

       end do

       foundElement = .false.

       backspace(30)

       do while(foundElement .eqv. .false.)

          read(30,*, iostat=status) ssymbol

          !! Some debug information in case of error!
          if (status > 0 ) then

             call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basis set file: "//trim(basisName)//" Please check this file.!","BasisSet module at getNumberOfCoreElectrons function.")

          end if

          if (status == -1 ) then

             call EffectiveCorePotentials_exception(ERROR, "The ECP: "//trim(basisName)//" for: "//trim(symbol)//" was not found!","BasisSet module at getNumberOfCoreElectrons function.")

          end if


          if((trim(ssymbol) == trim(symbol))) then

             backspace(30)
             read(30,*, iostat=status) ssymbol, token, nelec
             foundElement = .true.
             read(nelec, *) output

          end if

       end do

       print*, "NELEC: ", symbol, output

    else

       call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(basisName)//" This file don't exists!","EffectiveCorePotentials module at getNumberOfCoreElectrons function.")

    end if

  end function EffectiveCorePotentials_getNumberOfCoreElectrons


  !> @brief Load ECP parameters form ECP section in basis set file in deMon2K format
  !! @author I. Ortiz-Verano, 2017
  !! @version 1.0
  subroutine EffectiveCorePotentials_load(basisName, symbol, origin, unit)
    implicit none

    type(EffectiveCorePotentials) :: this
!    character(*),optional :: name
    character(*),optional :: basisName
    character(*), optional :: symbol
    real(8), optional :: origin(3)
    integer, optional :: unit    
    !    character(*),optional :: particleName

    integer :: i, j, k, l
    logical :: existFile, foundECP, foundElement, isNewECP
    character(30) :: token, name, ssymbol, nelec
    integer :: status
    integer :: maxAngularMoment
    integer :: numberOfCoreElectrons

    !    Followin variables are not of the ECP kind.
    !    integer :: nkParameters    !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
    !    real(8) :: zetakParameters !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
    !    real(8) :: dkParameters    !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283

    !! Setting name
    if(present(basisName)) this%name = trim(basisName)
    if(present(symbol)) symbol = trim(symbol)
    !    if(present(particleName)) this%particleSelected = trim(particleName)
    if(present(origin)) this%origin = origin

    !! Looking for the basis set file
    inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), exist = existFile)

    if(existFile) then

       !! Open basis set file that could contain ecp information
       open(unit=30, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), status="old", form="formatted")

       rewind(30)

       foundECP = .false.

       do while (foundECP .eqv. .false.)

          read(30,*, iostat=status) token

          if (status > 0 ) then   

             call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(basisName)//" Please check that this file contains pseudopotentials!","EffectiveCorePotentials module at Load function.")

          end if

          if (status == -1 ) then

             call EffectiveCorePotentials_exception(ERROR, "The basisSet: "//trim(basisName)//" for: "//trim(symbol)//" was not found!","BasisSet module at Load function.")

          end if



          if(trim(token(1:3)) == "ECP") then

             foundECP = .true.

             print*, "Now I'll looking for ECP for: ", symbol

          end if

       end do

       foundElement = .false.

       backspace(30)

       do while(foundElement .eqv. .false.)

          read(30,*, iostat=status) ssymbol

          !! Some debug information in case of error!
          if (status > 0 ) then

             call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basis set file: "//trim(basisName)//" Please check this file.!","BasisSet module at Load subroutine.")

          end if

          if (status == -1 ) then

             call EffectiveCorePotentials_exception(ERROR, "The ECP: "//trim(basisName)//" for: "//trim(symbol)//" was not found!","BasisSet module at Load subroutine.")

          end if


          if((trim(ssymbol) == trim(symbol))) then

             !             print*, "Symbol= ",ssymbol

             backspace(30)
             read(30,*, iostat=status) ssymbol, token, nelec
             foundElement = .true.

             read(nelec, *) numberOfCoreElectrons

             !   print*, "Number of Core electrons: ", numberOfCoreElectrons

             read(30,*,iostat=status) ssymbol, token

             !             print*, ssymbol, token

             select case(numberOfCoreElectrons)
             case(10)
                maxAngularMoment = 2
                this%length = 3
             case(18)
                maxAngularMoment = 2
                this%length = 4
             case(28)
                maxAngularMoment = 3
                this%length = 4
             case(36)
                maxAngularMoment = 3
                this%length = 4
             case(46)
                maxAngularMoment = 3
                this%length = 4
             case(60)
                maxAngularMoment = 4
                this%length = 5
             case(68)
                maxAngularMoment = 4
                this%length = 5
             case(78)
                maxAngularMoment = 5
                this%length = 5
             case default
                call EffectiveCorePotentials_exception(ERROR, "Max angular moment can't be stablished for: "//trim(symbol)//" at "//trim(basisName)//" file!"," BasisSet module at Load subroutine.")
                !!error
                
                print*, "error leyendo máximo momento angular"
             end select
             
             allocate(this%contraction(this%length))
             
             
             do i = 1, this%length
                
                
                !! Determine length of each contracted gaussian
                if (i == 1) then
                   !! When i = 1, this "if" find the length of this contraction
                   if((trim(token(1:2)) == "ul")) then
                      isNewECP = .false.
                      k=0
                      do while(isNewECP .eqv. .false.)
                         read(30,*, iostat=status) ssymbol, token
                         k=k+1
                         if((trim(ssymbol) == trim(symbol))) then
                            isNewECP = .true.
                            this%contraction(i)%length=k-1
                            !                            print*, "Contraction: ",i," length of the contraction: ", this%contraction(i)%length
                            exit
                         end if
                      end do
                      
                      allocate(this%contraction(i)%nkParameters(this%contraction(i)%length))
                      allocate(this%contraction(i)%zetakParameters(this%contraction(i)%length))
                      allocate(this%contraction(i)%dkParameters(this%contraction(i)%length))

                      do l=1, this%contraction(i)%length+1
                         backspace(30)
                      end do
                      
                      print*, "Contraction: ", i
                      do j = 1, this%contraction(i)%length
                         
                         read(30,*,iostat=status) this%contraction(i)%nkParameters(j), &
                              this%contraction(i)%zetakParameters(j), &
                              this%contraction(i)%dkParameters(j)
                         
                         print*, "nk, zetak y dk: ",&
                              this%contraction(i)%nkParameters(j), &
                              this%contraction(i)%zetakParameters(j), &
                              this%contraction(i)%dkParameters(j)
                         
                         !! Some debug information in case of error!
                         if (status > 0 ) then
                            
                            call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP at basisSet file: "//trim(this%name)//" Please check that file!","EffectiveCorePotentials module at Load function.")
                            
                         end if
                         
                      end do
                      
                   end if
                   
                   !! When i > 1:
                else if(i .ge. 2) then
                   read(30,*, iostat=status) token
                   if((trim(ssymbol) == trim(symbol))) then
                      isNewECP = .false.
                      k=1
                      do while(isNewECP .eqv. .false.)
                         read(30,*, iostat=status) ssymbol, token
                         k=k+1
                         if(((trim(token) == "S") .OR. &
                              (trim(token) == "P") .OR. &
                              (trim(token) == "D") .OR. &
                              (trim(token) == "F") .OR. &
                              (trim(token) == "G") .OR. &                              
                              (trim(token) == "nelec"))) then
                            isNewECP = .true.
                            this%contraction(i)%length=k-2
                            !                            print*, "Contraction: ",i," length of the contraction: ", this%contraction(i)%length
                            exit
                         end if
                      end do
                      allocate(this%contraction(i)%nkParameters(this%contraction(i)%length))
                      allocate(this%contraction(i)%zetakParameters(this%contraction(i)%length))
                      allocate(this%contraction(i)%dkParameters(this%contraction(i)%length))
                      
                      do l=1, this%contraction(i)%length+1
                         backspace(30)
                      end do

                      print*, "Contraction: ", i
                      do j = 1, this%contraction(i)%length
                         
                         read(30,*,iostat=status) this%contraction(i)%nkParameters(j), &
                              this%contraction(i)%zetakParameters(j), &
                              this%contraction(i)%dkParameters(j)
                         
                         
                         print*, "nk, zetak y dk: ",&
                              this%contraction(i)%nkParameters(j), &
                              this%contraction(i)%zetakParameters(j), &
                              this%contraction(i)%dkParameters(j)
                         
                         !! Some debug information in case of error!
                         if (status > 0 ) then
                            
                            call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP at basisSet file: "//trim(this%name)//" Please check that file!","EffectiveCorePotentials module at Load function.")

                         end if
                         
                      end do
                   end if
                   
                end if
                
                
             end do
             
             print*, "Max angular moment of the core for ", symbol, " : ", maxAngularMoment
             
          end if
          
       end do
       
    end if
    
  end subroutine EffectiveCorePotentials_load
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !>
  !! @brief Saves the EffectiveCorePotentials structure to file.
  !! @param this EffectiveCorePotentials object
  !! @author I. Ortiz-Verano, 2017
  subroutine EffectiveCorePotentials_saveToFile(this, unit)
    implicit none

      type(EffectiveCorePotentials), intent(in) :: this
      integer :: unit

      integer :: i

      write(unit,*) this%name
      write(unit,*) this%origin
      write(unit,*) this%length
      write(unit,*) this%ttype
!      write(unit,*) this%contractionLength
      write(unit,*) this%numberOfPrimitives

     do i = 1, size(this%contraction)
        call ContractedEcpGaussian_saveToFile(this%contraction(i), unit)
     end do

   end subroutine EffectiveCorePotentials_saveToFile


     !>
  !! @brief Shows information of the object
  !! @param this: basis set
  subroutine EffectiveCorePotentials_showInCompactForm( this, nameOfOwner)
    implicit none
    
    type(EffectiveCorePotentials) , intent(in) :: this
    character(*) ::nameOfOwner
    integer ::  i

    print *, nameOfOwner
    print *, this%name

    write(6,"(T5,A10,A11,A15)") trim(nameOfOwner),"    ECP: ", trim(this%name) !!!!!!!!!!!!!!!
    
    do i =1, this%length
       call ContractedEcpGaussian_showInCompactForm(this%contraction(i) )
    end do

  end subroutine EffectiveCorePotentials_showInCompactForm

  !     !>
  !     !! @brief Loads the effective core potentials structure from file.
  !     !! @param this EffectiveCorePotentials object
  !     !! @author I. Ortiz-Verano, 2017
  !     subroutine EffectiveCorePotentials_loadFromFile(this, unit)
  !       implicit none

  !       type(EffectiveCorePotentials) :: this
  !       integer :: unit

  !       integer :: i

  !       read(unit,*) this%name
  !       read(unit,*) this%origin
  !       read(unit,*) this%length
  !       read(unit,*) this%ttype
  ! !      read(unit,*) this%contractionLength
  !       read(unit,*) this%numberOfPrimitives

  !       allocate(this%contraction(this%length))

  !      do i = 1, size(this%contraction)
  !         call ContractedGaussian_loadFromFile(this%contraction(i), unit)
  !      end do

  !     end subroutine EffectiveCorePotentials_loadFromFile

  !>
  !! @brief handle exceptions
  subroutine EffectiveCorePotentials_exception( typeMessage, description, debugDescription)
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

  end subroutine EffectiveCorePotentials_exception

end module EffectiveCorePotentials_

