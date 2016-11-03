module Tensor_
  use MolecularSystem_
  use Vector_
  use ReadTransformedIntegrals_
  use ReadIntegrals_
  implicit none

  type :: Tensor
      type(Vector) :: container
      logical :: isInterSpecies=.false.
      logical :: isMolecular=.true.
      integer :: otherSpeciesID=0
      logical :: IsInstanced

  end type Tensor

  type(Tensor), allocatable :: ints(:)


contains

! option 1
!  type(Tensor) :: mytensor
!  call Tensor_constructor(mytensor%container, speciesID, isMolecular, this2, otherSpeciesID, isInterspecies)

! option 2
!  type(Tensor) :: mytensor
!  call Tensor_constructor(mytensor, speciesID, isMolecular, this2, otherSpeciesID, isInterspecies)
!  subroutine Tensor_constructor(this, speciesID, isMolecular, this2, otherSpeciesID, isInterspecies)
!      type(Tensor), intent(inout):: this
!      this%isMolecular = isMolecular 
!      call ReadTransformedIntegrals_readOneSpecies(speciesID, this%container)


  subroutine Tensor_constructor(this, speciesID, isMolecular, this2, otherSpeciesID, isInterspecies)
      implicit none

      type(Vector), intent(inout):: this
      type(Vector), intent(inout), optional :: this2
      logical, optional :: isInterspecies
      logical, optional :: isMolecular
      real(8), allocatable :: intls(:)
      integer, optional :: otherSpeciesID
      integer :: speciesID
      integer :: nao, sze, onao, osze, tsze

!!      if ( .not. present(isMolecular) ) isMolecular = .false.      
!!      if ( .not. present(isMolecular) ) isMolecular = this%isMolecular      


      if (isMolecular) then

        call ReadTransformedIntegrals_readOneSpecies(speciesID, this)

        
          if (present(otherSpeciesID)) then 

            if (otherSpeciesID > 1) then

              call ReadTransformedIntegrals_readTwoSpecies(speciesID, otherSpeciesID, this2)

            end if

          end if

      else

          nao = MolecularSystem_getTotalNumberOfContractions(speciesID)

        if(isInterspecies) then

          sze = nao * (nao + 1) / 2
          sze = sze * (sze + 1) / 2

          if(allocated(intls)) deallocate(intls)
          allocate(intls(sze))
          intls = 0.0_8
          
          call ReadIntegrals_intraSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), this)

          this%values = intls

        else

          onao = MolecularSystem_getTotalNumberOfContractions(otherspeciesID)
          sze = nao * (nao + 1) / 2
          osze = onao * (onao + 1) / 2
          tsze = sze * osze
    
          if(allocated(intls)) deallocate(intls)
          allocate(intls(tsze))
          intls = 0.0_8

          call ReadIntegrals_interSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), &
              trim(MolecularSystem_getNameOfSpecie(otherSpeciesID)), osze, this)

          this2%values = intls

        end if
      end if

    

  end subroutine Tensor_constructor

  subroutine Tensor_destructor()
      implicit none
      
      if(allocated(ints)) deallocate(ints)
    
  end subroutine Tensor_destructor

  ! function Tensor_getValue4(this, a, b, r, s) result(output)
  !     implicit none
  !     type(Tensor) :: this
  !     integer :: a, b, r, s
  !     real(8) :: output

  !     ! Convert 4 intex to 1
  !     index = 
  !     ! Get value from container
  !     this%container(index)

    
  ! end function Tensor_getValue4

  !....
  ! t(a, b, r, s)

    
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
