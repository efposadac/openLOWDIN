!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!    http://www.qcc.unal.edu.co/
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief Tensor module
!!        This module calls and initializes all information necessary to use transformedintegrals and return a one index from four indices.
!! @author  Carlos Andres Ortiz Mahecha (CAOM) (caraortizmah@unal.edu.co)
!!
!! <b> Creation date : </b> 2016-11-02
!!
!! <b> History: </b>
!!
!!   - <tt> 2016-10-26 </tt>: (CAOM) ( caraortizmah@unal.edu.co )
!!        -# Development of Tensor module:
!!                This Tensor ... jue jue 
!!   - <tt> data </tt>:  
!!
!!
!! @warning <em>  All characters and events in this module -- even those based on real source code -- are entirely fictional. </br>
!!                All celebrity lines are impersonated.....poorly. </br> 
!!                The following module contains corase language and due to it's cintent should not be viewed by anyone. </em>
!!
!!
!!
module Tensor_
  use MolecularSystem_
  use Vector_
  use ReadTransformedIntegrals_
  use ReadIntegrals_
  implicit none

  type :: Tensor
      type(Vector) :: container
      integer :: speciesID
      integer :: otherSpeciesID
      logical :: isInterSpecies = .false.
      logical :: isMolecular = .false.
      logical :: IsInstanced = .false.
  end type Tensor

  type(Tensor), public, target :: int1
  type(Tensor), public, target :: int2

    !>
  !! @brief Abstract class
  !! @author Carlos Andres Ortiz-Mahecha (CAOM)

  interface Tensor_index
    module procedure Tensor_index4Intra, Tensor_index4Inter
  end interface

  interface Tensor_getValue
    module procedure Tensor_getValue_intra, Tensor_getValue_inter
  end interface

  private :: &
    Tensor_index2, &
    Tensor_index4Intra, &
    Tensor_index4Inter, &
    Tensor_getValue_intra, &
    Tensor_getValue_inter


contains

  !>
  !! @brief Constructor of the class
  !! @author CAOM
  
  subroutine Tensor_constructor(this, speciesID, otherSpeciesID, isMolecular)
      implicit none

      type(Tensor), intent(inout) :: this
      integer, intent(in) :: speciesID
      integer, optional, intent(in) :: otherSpeciesID
      logical, optional, intent(in) :: isMolecular

      integer sze, nao, onao, osze, tsze
      
      this%speciesID = speciesID

      if(present(otherspeciesID)) then
        this%otherspeciesID = otherspeciesID
        this%isInterspecies = .true.
      end if

      if(present(isMolecular)) this%isMolecular = isMolecular

      print*," Begin Tensor_constructor", this%isMolecular,  this%isInterSpecies

      print*, "size of container", size(this%container%values)

      if (this%isMolecular) then

        print*," Begin Tensor_constructor"

          if (this%isInterSpecies) then
            print*, "test two species",this%speciesID, this%otherSpeciesID
            call ReadTransformedIntegrals_readTwoSpecies(this%speciesID, this%otherSpeciesID, this%container)
          else
            print*, "test one species"
            call ReadTransformedIntegrals_readOneSpecies(this%speciesID, this%container)
          end if

      else
  
        nao = MolecularSystem_getTotalNumberOfContractions(speciesID)

        if (this%isInterspecies) then

          sze = nao * (nao + 1) / 2
          sze = sze * (sze + 1) / 2
        
          call ReadIntegrals_intraSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), this%container)
  
        else

          onao = MolecularSystem_getTotalNumberOfContractions(otherspeciesID)
          sze = nao * (nao + 1) / 2
          osze = onao * (onao + 1) / 2
          tsze = sze * osze

          call ReadIntegrals_interSpecies(trim(MolecularSystem_getNameOfSpecie(speciesID)), &
              trim(MolecularSystem_getNameOfSpecie(otherSpeciesID)), osze, this%container)

        end if

      end if

      this%IsInstanced = .true.

  end subroutine Tensor_constructor

  subroutine Tensor_destructor(this)
      implicit none

      type(Tensor), intent(inout) :: this

      call Vector_destructor(this%container)

      this%otherspeciesID = -1
      this%isInterSpecies = .false.
      this%isMolecular = .false.
      this%IsInstanced = .false.

  end subroutine Tensor_destructor

  function Tensor_getValue_intra(this, a, b, r, s, w) result(output)
      implicit none
      type(Tensor), intent(in) :: this
      integer, intent(in) :: a, b, r, s, w

      integer(8) :: index
      real(8) :: output

      ! index = Tensor_index4Intra(a, b, r, s)
      index = IndexMap_tensorR4ToVectorB(a, b, r, s, w)
      output = this%container%values(index)
      
  end function Tensor_getValue_intra

  function Tensor_getValue_inter(this, a, b, r, s, w, w2) result(output)
      implicit none
      type(Tensor), intent(in) :: this
      integer, intent(in) :: a, b, r, s, w, w2

      integer(8) :: index
      real(8) :: output

      ! index = Tensor_index4Inter(a, b, r, s, w)
      index = IndexMap_tensorR4ToVectorB(a, b, r, s, w, w2)
      output = this%container%values(index)
      
  end function Tensor_getValue_inter

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

  function IndexMap_tensorR4ToVectorB( i, j, k, l, basisSizeA, basisSizeB ) result ( output )
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k
    integer, intent(in) :: l
    integer, optional :: basisSizeA
    integer, optional :: basisSizeB

    integer(kind=8) :: output
    integer :: auxSize
    integer :: ij, kl

      if ( .not. present ( basisSizeB ) ) then

        ij = int(IndexMap_tensorR2ToVector( i, j, basisSizeA),4)
        kl = int(IndexMap_tensorR2ToVector( k, l, basisSizeA),4)

        auxSize = ( basisSizeA * ( basisSizeA + 1 ) ) / 2
        output = IndexMap_tensorR2ToVector (ij, kl, auxSize)
      else 

        ij = int(IndexMap_tensorR2ToVector( i, j, basisSizeA),4) - 1
        kl = int(IndexMap_tensorR2ToVector( k, l, basisSizeB),4)

        auxSize = ( basisSizeB * ( basisSizeB + 1 ) ) / 2
        output = int(auxSize,8) * int(ij,8) + int(kl,8)
                       
      end if

  end function IndexMap_TensorR4ToVectorB

  function IndexMap_tensorR2ToVector( i, j, basisSizeA ) result ( output )
      implicit none
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, optional :: basisSizeA

      integer(kind=8) :: output

        if ( i > j ) then
           output = i - j + ( ( ( 2 * basisSizeA * (j -1 )) - ( j * j) + (3*j) ) / 2 )
        else 
           output = j - i + ( ( ( 2 * basisSizeA * (i -1 )) - ( i * i) + (3*i) ) / 2 )
        end if

  end function IndexMap_TensorR2ToVector

  function Tix2(i, j, w) result(output)
    implicit none
    integer :: i, j, w
    integer :: output

    if(i > j) then
       output = (w * j) - (w - i) - ((j * (j + 1)) / 2)
    else
       output = (w * i) - (w - j) - ((i * (i + 1)) / 2)
    end if

  end function Tix2

end module Tensor_