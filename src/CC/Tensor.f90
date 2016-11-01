module Tensor_
  use Vector_
  use ReadTransformedIntegrals_
  implicit none

  type :: Tensor
      type(Vector) :: container
      logical :: isInterSpecies=.false.
      integer :: otherSpeciesID = 0

  end type Tensor



contains

  subroutine Tensor_constructor(this, speciesID, isInterspecies, isMolecular, otherSpeciesID)
      implicit none

      type(Vector), intent(in):: this
      logical, optional :: isInterspecies
      logical, optional :: isMolecular
      integer, optional :: otherSpeciesID
      integer :: speciesID
      ! call ReadTransformedIntegrals_readOneSpecies(speciesID, )
    

  end subroutine Tensor_constructor

  ! subroutine Tensor_destructor(args)
  !     implicit none
  !     real :: args
    
  ! end subroutine Tensor_destructor

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
