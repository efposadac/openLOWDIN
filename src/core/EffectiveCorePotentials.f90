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
  ! use InputManager_ 
  implicit none

  type :: EffectiveCorePotentials
     !     type(contractedEcpGaussian),allocatable :: contraction(:)
     type(contractedEcpGaussian),allocatable :: contraction(:)
     !     type(contractedEcpGaussian),allocatable :: nkParameters(:)   !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     !     type(contractedEcpGaussian),allocatable :: zetaParameters(:) !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     !     type(contractedEcpGaussian),allocatable :: dkParameters(:)   !< equation (16) from J. Chem. Phys. 82, 1, 1985, 270-283
     character(30) :: name 
     real(8) :: origin(3)
     integer :: ttype   !!!!!!!!!!!?????????????
     integer :: numberOfCoreElectrons
     integer :: length
     integer :: contractionLenght
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
    character(*),optional :: basisName
    character(*), optional :: symbol
    integer, optional :: unit
    real(8), optional :: origin(3)
    !    character(*),optional :: particleName

    integer :: i, j
    logical :: existFile, foundECP, foundElement
    character(30) :: token, name, ssymbol, nelec
    integer :: status
    integer :: maxAngularMoment
    integer :: numberOfCoreElectrons

    !    Followin variables are note of the ECP kind.
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

!             print*, "2: ", ssymbol, token, nelec

             read(nelec, *) numberOfCoreElectrons

             read(30,*,iostat=status) symbol, token

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

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!                Falta asignar this%contraction(i)%lenght
!!!!!!                IOSTATE 186
!!!!!!                http://www.msg.ucsf.edu/local/programs/IBM_Compilers/Fortran/html/pgs/lr76.htm
!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! do while ()
                !    this%contraction(i)%lenght=3
                ! end do
                
                

                print*, "this:  ",i-1
                print*, token

                if((trim(token) == "ul")) then
                   print*, "Estoy en el if"
                   backspace(30)


                allocate(this%contraction(i)%nkParameters(this%contraction(i)%length))
                allocate(this%contraction(i)%zetakParameters(this%contraction(i)%length))
                allocate(this%contraction(i)%dkParameters(this%contraction(i)%length))
                
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


                 !! Some debug information in case of error!
                if (status > 0 ) then

                   call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP at basisSet file: "//trim(this%name)//" Please check that file!","EffectiveCorePotentials module at Load function.")

                end if

                ! allocate(this%contraction(i)%orbitalExponents(this%contraction(i)%length))
                ! allocate(this%contraction(i)%contractionCoefficients(this%contraction(i)%length))

                ! do j = 1, this%contraction(i)%length

                !    read(30,*,iostat=status) this%contraction(i)%orbitalExponents(j), &
                !         this%contraction(i)%contractionCoefficients(j)

                !    !! Some debug information in case of error!
                !    if (status > 0 ) then

                !       call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")

                !    end if

                ! end do

             end do



             print*, "Max angular moment of the core for ", ssymbol, " : ", maxAngularMoment

          end if

       end do


    end if



  end subroutine EffectiveCorePotentials_load


  !   subroutine EffectiveCorePotentials_load(this, basisName, ParticleName, origin, unit)
  !     implicit none

  !     type(BasisSet) :: this
  ! !    character(*) :: formatType
  !     character(*), optional :: basisName
  !     character(*), optional :: particleName
  !     real(8), optional :: origin(3)
  !     integer, optional :: unit

  !     integer :: i, j !< Iterators
  !     logical :: existFile, found
  !     logical :: foundECP
  !     character(20) :: token
  !     character(10) :: symbol
  !     character(30) :: particleSelected
  !     integer :: status
  !     integer :: nelec, numberOfCoreElectrons

  !     !! Setting name
  !     if(present(basisName)) this%name = trim(basisName)
  !     if(present(particleName)) particleSelected = trim(particleName)
  !     if(present(origin)) this%origin = origin


  !        !! Open effective core potentials from basis set file at library
  !        inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), exist = existFile)

  !        if(existFile) then

  !           !! Open File
  !           open(unit=30, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), status="old",form="formatted")
  !           rewind(30)

  !           foundECP = .false.

  !           !! Find ECP group in basis set file
  !           do while(foundECP .eqv. .false.)

  !              read(30,*, iostat=status) token

  ! !!!!!!!  To do: include debug information

  ! !             if (status > 0 ) then
  ! !                
  ! !                call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(this%name)//" Please check that this file contains pseudopotentials!","EffectiveCorePotentials module at Load function.")
  ! !                
  ! !             end if
  ! !             
  ! !             if (status == -1 ) then
  ! !                
  ! !                call EffectiveCorePotentials_exception(ERROR, "The basisSet: "//trim(this%name)//" for: "//trim(particleSelected)//" was not found!","BasisSet module at Load function.")
  ! !                
  ! !             end if

  !              if(trim(token(1:3)) == "ECP") then

  !                 foundECP = .true.

  ! !!!!!!!!!!!!!!!!!!!! ¿Cómo hago para que solo lea desde aquí la unidad 30?
  !                 do while(foundECP .eqv. .true.)

  !                    backspace(30)

  !                    read(30,*, iostat=status) symbol, token, nelec

  !                    !! Some debug information in case of error!
  !                    if (status > 0 ) then

  !                       call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basis set file: "//trim(this%name)//" Please check this file.!","BasisSet module at Load function.")

  !                    end if

  !                    if (status == -1 ) then

  !                       call BasisSet_exception(ERROR, "The ECP: "//trim(this%name)//" for: "//trim(particleSelected)//" was not found!","BasisSet module at Load function.")

  !                    end if

  !                 found = .false.

  !                 if((trim(symbol) == trim(particleSelected)) .and. (trim(token) == "nelec")) then

  !                    found = .true.

  !                    numberOfCoreElectrons = nelec

  !                 end if

  !              else

  !                 call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(this%name)//" This file don't contains ECP for " //trim(particleSelected)//" !","EffectiveCorePotentials module at Load function.")

  !              end if
  !           end do


  !           !! Neglect any coment
  !           token = "#"
  !           do while(trim(token(1:1)) == "#")

  !              read(30,*) token

  !           end do

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! from here again

  !           !! Start reading ECP from basis set file

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! From Hay1985:
  ! !          \sum_k d_k r^{nk} e^{-\zeta_k r^2}
  ! !          momento angular l se refiere a $r^2 (U_l - U_L)$
  ! !          ul (potenciales $U_L$) se refiere a $r^2 (U_L - N_c/r)$5

  !           backspace(30)

  !           read(30,*, iostat=status) this%length

  !           !! Some debug information in case of error!
  !           if (status > 0 ) then

  !              call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")

  !           end if

  !           allocate(this%contraction(this%length))

  !           do i = 1, this%length

  !              read(30,*,iostat=status) this%contraction(i)%id, &
  !                   this%contraction(i)%angularMoment, &
  !                   this%contraction(i)%length

  !              !! Some debug information in case of error!
  !              if (status > 0 ) then

  !                 call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")

  !              end if

  !              allocate(this%contraction(i)%orbitalExponents(this%contraction(i)%length))
  !              allocate(this%contraction(i)%contractionCoefficients(this%contraction(i)%length))

  !              do j = 1, this%contraction(i)%length

  !                 read(30,*,iostat=status) this%contraction(i)%orbitalExponents(j), &
  !                      this%contraction(i)%contractionCoefficients(j)

  !                 !! Some debug information in case of error!
  !                 if (status > 0 ) then

  !                    call BasisSet_exception(ERROR, "ERROR reading basisSet file: "//trim(this%name)//" Please check that file!","BasisSet module at Load function.")

  !                 end if

  !              end do

  !              !! Ajust and normalize contractions
  !              this%contraction(i)%origin = this%origin

  !              !! Calculates the number of cartesian orbitals, by dimensionality
  !              select case(CONTROL_instance%DIMENSIONALITY)

  !              case(3)
  !                 this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 )*( this%contraction(i)%angularMoment + 2_8 ) ) / 2_8
  !              case(2)
  !                 this%contraction(i)%numCartesianOrbital = ( ( this%contraction(i)%angularMoment + 1_8 ) )
  !              case(1)
  !                 this%contraction(i)%numCartesianOrbital = 1 
  !              case default
  !                 call BasisSet_exception( ERROR, "Class object Basis set in load function",&
  !                      "This Dimensionality is not avaliable") 
  !              end select

  !              !! Normalize
  !              allocate(this%contraction(i)%contNormalization(this%contraction(i)%numCartesianOrbital))
  !              allocate(this%contraction(i)%primNormalization(this%contraction(i)%length, &
  !                   this%contraction(i)%length*this%contraction(i)%numCartesianOrbital))

  !              this%contraction(i)%contNormalization = 1.0_8
  !              this%contraction(i)%primNormalization = 1.0_8

  !              call ContractedGaussian_normalizePrimitive(this%contraction(i))
  !              call ContractedGaussian_normalizeContraction(this%contraction(i))

  !              !! DEBUG
  !              !! call ContractedGaussian_showInCompactForm(this%contraction(i))

  !           end do

  !           close(30)

  !           !!DONE

  !        else

  !           call EffectiveCorePotentials_exception(ERROR, "The basisSet file: "//trim(basisName)//" was not found!","EffectiveCorePotentials module at Load function.")

  !        end if

  ! !    end select

  !   end subroutine EffectiveCorePotentials_load

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   !>
  !   !! @brief Saves the EffectiveCorePotentials structure to file.
  !   !! @param this EffectiveCorePotentials object
  !   !! @author I. Ortiz-Verano, 2017
  !   subroutine EffectiveCorePotentials_saveToFile(this, unit)
  !     implicit none

  !     type(EffectiveCorePotentials), intent(in) :: this
  !     integer :: unit

  !     integer :: i

  !     write(unit,*) this%name
  !     write(unit,*) this%origin
  !     write(unit,*) this%length
  !     write(unit,*) this%ttype
  !     write(unit,*) this%contractionLength
  !     write(unit,*) this%numberOfPrimitives

  ! !    do i = 1, size(this%contraction)
  ! !       call ContractedGaussian_saveToFile(this%contraction(i), unit)
  ! !    end do

  !   end subroutine EffectiveCorePotentials_saveToFile

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

