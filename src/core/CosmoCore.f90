module CosmoCore_
	use Units_
  implicit none

  type, public :: surfaceSegment
     real(8), allocatable :: xs(:)
     real(8), allocatable :: ys(:)
     real(8), allocatable :: zs(:)
     real(8), allocatable :: area(:)
     integer :: sizeSurface
  end type surfaceSegment

  ! >Singleton
  type(surfaceSegment), public, target :: surfaceSegment_instance

contains

  !----------------------subroutines------------------------------
  subroutine CosmoCore_lines(surface) 
    !subrutina que cuenta las lineas en el archivo vectors.out

    ! surface segment atributes
    implicit none 

    integer :: n
    character(len=60):: cmd

    type(surfaceSegment), intent(inout) :: surface

    cmd = "cat vectors.vec | grep '[^ ]' | wc -l > nlines.txt"
    call system(cmd) 
    open(1,file='nlines.txt')
    read(1,*) n
    cmd = 'rm nlines.txt'
    call system(cmd)
    surface%sizeSurface=n
    write(*,*) "la superficie tiene", n, "segmentos"
    return

  end subroutine CosmoCore_lines

  !----------------------subroutines------------------------------
  subroutine CosmoCore_Filler(surface)

    ! subrutina que construye la matriz c de acuerdo a Paifeng Su
    ! a partir de el archivo generado por el cálculo usanto gepol y alimenta la
    ! instancia surfaceSegment

    implicit none 

    integer :: i, j

    type(surfaceSegment), intent(inout) :: surface

    real(8), dimension(surface%sizeSurface) :: x !segment x cordinate
    real(8), dimension(surface%sizeSurface) :: y !segment y cordinate
    real(8), dimension(surface%sizeSurface) :: z !segment z cordinate
    real(8), dimension(surface%sizeSurface) :: a	!segment area

    write(*,*)"estamos adentro del filler"

    ! llenado de surface con la información que está en vectors.vec


100 format (2X,F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
    open(55,file='vectors.vec',status='OLD') 
    read(55,100) (x(i),y(i),z(i),a(i),i=1,surface%sizeSurface)

    !asignando espacio en memoria para los parametros

    allocate(surface%xs(surface%sizeSurface))
    allocate(surface%ys(surface%sizeSurface))
    allocate(surface%zs(surface%sizeSurface))
    allocate(surface%area(surface%sizeSurface))

    ! write(*,*)"tipo superficie"
    !! llenando surface con la informacion leida

    do i=1,surface%sizeSurface        
       surface%xs(i)=x(i)/AMSTRONG
       surface%ys(i)=y(i)/AMSTRONG
       surface%zs(i)=z(i)/AMSTRONG
       surface%area(i)=a(i)/((AMSTRONG)**2)
       ! write(*,100)surface%xs(i),surface%ys(i),surface%zs(i),surface%area(i)
    end do

  end subroutine CosmoCore_Filler
end module CosmoCore_
