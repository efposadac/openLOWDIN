
module Tensor_
  use MolecularSystem_
  use Vector_
  use ReadTransformedIntegrals_
  use ReadIntegrals_
  implicit none

  type :: Tensor
      type(Vector) :: container
      logical :: isInterSpecies!=.false.
      logical :: isMolecular!=.true.
      integer :: otherSpeciesID!=0
      logical :: IsInstanced

  end type Tensor

  type(Tensor), public, target :: int1
  type(Tensor), public, target :: int2

    !>
  !! @brief Abstract class
  !! @author Carlos Andres Ortiz-Mahecha (CAOM)

  interface Tensor_index
    module procedure Tensor_index2, Tensor_index4Intra, Tensor_index4Inter
  end interface

  private :: &
    Tensor_index2, &
    Tensor_index4Intra, &
    Tensor_index4Inter


contains

  !>
  !! @brief Constructor of the class
  !! @author CAOM
  
  subroutine Tensor_constructor(this, speciesID, otherSpeciesID, isMolecular, isInterspecies)
      implicit none

      type(Tensor), intent(inout), pointer :: this
      ! type(Tensor), intent(inout), pointer, optional :: this2
      logical, optional :: isInterspecies
      logical, optional :: isMolecular
      integer, optional :: otherSpeciesID
      integer :: speciesID
      integer :: nao, sze, onao, osze, tsze

      
      if ( .not. present(isMolecular) ) isMolecular = .true.
      ! if ( .not. present(this2) ) then 
      !   if(allocated(this2%container%values)) deallocate(this2%container%values)
      ! end if
      if ( .not. present(otherSpeciesID) ) otherSpeciesID = 0!this%otherSpeciesID
      if ( .not. present(isInterspecies) ) isInterspecies = .false.!this%isInterspecies



      if (isMolecular) then

        call ReadTransformedIntegrals_readOneSpecies(speciesID, this%container)

        ! this%container => int1%container

          if (present(otherSpeciesID)) then 

            if (otherSpeciesID > 1) then

              call ReadTransformedIntegrals_readTwoSpecies(speciesID, otherSpeciesID, this%container)

              ! this2%container => int2%container

            end if

          end if

      else

          nao = MolecularSystem_getTotalNumberOfContractions(speciesID)

        if(isInterspecies) then

          sze = nao * (nao + 1) / 2
          sze = sze * (sze + 1) / 2

          if(allocated(this%container%values)) deallocate(this%container%values)
          allocate(this%container%values(sze))
          this%container%values = 0.0_8
          
          call ReadIntegrals_intraSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), this%container)

          ! this%container => int1%container
  
        else

          onao = MolecularSystem_getTotalNumberOfContractions(otherspeciesID)
          sze = nao * (nao + 1) / 2
          osze = onao * (onao + 1) / 2
          tsze = sze * osze
    
          if(allocated(this%container%values)) deallocate(this%container%values)
          allocate(this%container%values(tsze))
          this%container%values = 0.0_8

          call ReadIntegrals_interSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), &
              trim(MolecularSystem_getNameOfSpecie(otherSpeciesID)), osze, this%container)

          ! this%container => int1%container

        end if
      end if

    

  end subroutine Tensor_constructor

  subroutine Tensor_destructor()
      implicit none
      
      if(allocated(int1%container%values)) deallocate(int1%container%values)
      if(allocated(int2%container%values)) deallocate(int2%container%values)
    
  end subroutine Tensor_destructor

  ! function Tensor_indexnumber(this, a, b, r, s) result(output)
  !     implicit none
  !     type(Tensor) :: this
  !     integer :: a, b, r, s
  !     real(8) :: output

  !     ! Convert 4 intex to 1
  !     index = 
  !     ! Get value from container
  !     this%container%container(index)

    
  ! end function Tensor_indexnumber
    
  function Tensor_index2(i, j) result(output)
    implicit none
    integer :: i, j
    integer :: output

    if(i > j) then
       output = i * (i + 1) / 2 + j
    else
       output = j * (j + 1) / 2 + i
    end if

  end function Tensor_index2

  function Tensor_index4Intra(i, j, k, l) result(output)
    implicit none
    integer :: i, j, k, l
    integer :: ii, jj, kk, ll
    integer :: output

    integer ij, kl

    ii = i - 1
    jj = j - 1
    kk = k - 1
    ll = l - 1

    ij = Tensor_index2(ii, jj)
    kl = Tensor_index2(kk, ll)

    output = Tensor_index2(ij, kl) + 1

  end function Tensor_index4Intra

  function Tensor_index4Inter(i, j, k, l, w) result(output)
    implicit none
    integer :: i, j, k, l
    integer :: w
    integer :: ii, jj, kk, ll
    integer :: output

    integer ij, kl

    ii = i - 1
    jj = j - 1
    kk = k - 1
    ll = l - 1

    ij = Tensor_index2(ii, jj)
    kl = Tensor_index2(kk, ll)

    output = ij * w + kl + 1

  end function Tensor_index4Inter

end module Tensor_



! TEST
! type(Tensor) :: test
! call Tensor_constructor(test, isInterspecie=.false., isMolecular=.false., speciesID=speciesID, otherSpeciesID=0)
! integral = Tensor_getValue(this, a, b, r, s)
! integral = Tensor_getValue(this, a, b)
