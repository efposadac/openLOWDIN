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
!  use ContractedGaussian_
  use CONTROL_
  use InputManager_ 
  implicit none
  
  type :: EffectiveCorePotentials
     character(30) :: name !!!!!! cómo se va a reconocer?
     real(8) :: origin(3)
     integer :: ttype   !!!!!!!!!!!?????????????
     integer :: numberOfCoreElectrons
     integer :: nkParameter    !< equation (16) from J. Chem. Phys. 82, 1, 1985
     real(8) :: zetakParameter !< equation (16) from J. Chem. Phys. 82, 1, 1985
     real(8) :: dkParameter    !< equation (16) from J. Chem. Phys. 82, 1, 1985
     character(1) :: maxAngularMoment
     character(1) :: angularMoment
  end type EffectiveCorePotentials
  
  public :: &
       EffectiveCorePotentials_load, &
       EffectiveCorePotentials_loadECP, &
       EffectiveCorePotentials_showInSimpleForm  

  private :: &
       EffectiveCorePotentials_exception
  
contains


  !> @brief Load ECP form basis set file in deMon2K format
  !! @author I. Ortiz-Verano, 2017
  !! @version 1.0
  subroutine EffectiveCorePotentials_load(this, basisName, ParticleName, origin, unit)
    implicit none
    
    type(BasisSet) :: this
!    character(*) :: formatType
    character(*), optional :: basisName
    character(*), optional :: particleName
    real(8), optional :: origin(3)
    integer, optional :: unit
    
    integer :: i, j !< Iterators
    logical :: existFile, found
    logical :: foundECP
    character(20) :: token
    character(10) :: symbol
    character(30) :: particleSelected
    integer :: status
    integer :: nelec, numberOfCoreElectrons
    
    !! Setting name
    if(present(basisName)) this%name = trim(basisName)
    if(present(particleName)) particleSelected = trim(particleName)
    if(present(origin)) this%origin = origin
    
       
       !! Open effective core potentials from basis set file at library
       inquire(file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), exist = existFile)
       
       if(existFile) then
          
          !! Open File
          open(unit=30, file=trim(CONTROL_instance%DATA_DIRECTORY)//trim(CONTROL_instance%BASIS_SET_DATABASE)//trim(basisName), status="old",form="formatted")
          rewind(30)
          
          foundECP = .false.

          !! Find ECP group in basis set file
          do while(foundECP .eqv. .false.)
                          
             read(30,*, iostat=status) token

!!!!!!!  To do: include debug information

!             if (status > 0 ) then
!                
!                call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(this%name)//" Please check that this file contains pseudopotentials!","EffectiveCorePotentials module at Load function.")
!                
!             end if
!             
!             if (status == -1 ) then
!                
!                call EffectiveCorePotentials_exception(ERROR, "The basisSet: "//trim(this%name)//" for: "//trim(particleSelected)//" was not found!","BasisSet module at Load function.")
!                
!             end if
             
             if(trim(token(1:3)) == "ECP") then

                foundECP = .true.
                
!!!!!!!!!!!!!!!!!!!! ¿Cómo hago para que solo lea desde aquí la unidad 30?
                do while(foundECP .eqv. .true.)

                   backspace(30)
                   
                   read(30,*, iostat=status) symbol, token, nelec
                   
                   !! Some debug information in case of error!
                   if (status > 0 ) then

                      call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basis set file: "//trim(this%name)//" Please check this file.!","BasisSet module at Load function.")
                
                   end if

                   if (status == -1 ) then

                      call BasisSet_exception(ERROR, "The ECP: "//trim(this%name)//" for: "//trim(particleSelected)//" was not found!","BasisSet module at Load function.")

                   end if
                
                found = .false.

                if((trim(symbol) == trim(particleSelected)) .and. (trim(token) == "nelec")) then
                   
                   found = .true.

                   numberOfCoreElectrons = nelec
                   
                end if
                
             else

                call EffectiveCorePotentials_exception(ERROR, "ERROR reading ECP from basisSet file: "//trim(this%name)//" This file don't contains ECP for " //trim(particleSelected)//" !","EffectiveCorePotentials module at Load function.")
                
             end if
          end do

         
          !! Neglect any coment
          token = "#"
          do while(trim(token(1:1)) == "#")
             
             read(30,*) token
             
          end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! from here again
          
          !! Start reading ECP from basis set file
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! From Hay1985:
!          \sum_k d_k r^{nk} e^{-\zeta_k r^2}
!          momento angular l se refiere a $r^2 (U_l - U_L)$
!          ul (potenciales $U_L$) se refiere a $r^2 (U_L - N_c/r)$5
          
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
          
          call EffectiveCorePotentials_exception(ERROR, "The basisSet file: "//trim(basisName)//" was not found!","EffectiveCorePotentials module at Load function.")
             
       end if
       
!    end select
    
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
    write(unit,*) this%contractionLength
    write(unit,*) this%numberOfPrimitives

!    do i = 1, size(this%contraction)
!       call ContractedGaussian_saveToFile(this%contraction(i), unit)
!    end do
    
  end subroutine EffectiveCorePotentials_saveToFile

  !>
  !! @brief Loads the effective core potentials structure from file.
  !! @param this EffectiveCorePotentials object
  !! @author I. Ortiz-Verano, 2017
  subroutine EffectiveCorePotentials_loadFromFile(this, unit)
    implicit none
    
    type(EffectiveCorePotentials) :: this
    integer :: unit
    
    integer :: i
    
    read(unit,*) this%name
    read(unit,*) this%origin
    read(unit,*) this%length
    read(unit,*) this%ttype
    read(unit,*) this%contractionLength
    read(unit,*) this%numberOfPrimitives

    allocate(this%contraction(this%length))

!    do i = 1, size(this%contraction)
!       call ContractedGaussian_loadFromFile(this%contraction(i), unit)
!    end do
    
  end subroutine EffectiveCorePotentials_loadFromFile
 
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
