module CosmoTools_
  use Matrix_
  implicit none 


contains


  !!----------------------subroutines------------------------------

  subroutine CosmoTools_caller()
    !subrutina que llama el ejecutable de gepol
    implicit none
    character(len=60) :: cmd

    cmd = "gepol.x < gepol.inp > gepol.out"
    call system(cmd) 
    write(*,*) "Calculating tesseras"
    !		cmd = "rm gepol.out"
    call system(cmd) 


  end subroutine CosmoTools_caller

  !!----------------------subroutines------------------------------
  subroutine CosmoTools_lines(n) 
    !!subrutina que cuenta las lineas en el archivo vectors.out
    implicit none 
    integer,intent(out) :: n
    character(len=60):: cmd

    cmd = "cat vectors.vec | grep '[^ ]' | wc -l > nlines.txt"
    call system(cmd) 
    open(1,file='nlines.txt')
    read(1,*) n
    !   cmd = 'rm nlines.txt'
    call system(cmd)
    return

  end subroutine CosmoTools_lines

  !!----------------------subroutines------------------------------
  subroutine CosmoTools_Cmatrix(np)

    !subrutina que construye la matriz c de acuerdo a Paifeng Su
    !a partir de el archivo generado por el cálculo usanto gepol
    implicit none 

    integer :: i, j, np
    real(8), dimension(np) :: x !segment x cordinate
    real(8), dimension(np) :: y !segment y cordinate
    real(8), dimension(np) :: z !segment z cordinate
    real(8), dimension(np) :: a	!segment area

    type(Matrix) :: cmat

    !real(8), allocatable :: cmat(:,:)
    !allocate(amat(np,np))

    !! llamado al constructor de matrices

    call Matrix_constructor(cmat, int(np,8), int(np,8))


100 format(F12.8,2X,F12.8,2X,F12.8,2X,F12.8)

    open(55,file='vectors.vec',status='OLD') 
    read(55,100) (x(i),y(i),z(i),a(i),i=1,np)


    do i=1,np
       do j=1,np
          if (i==j) then
             cmat%values(i,j)=3.8*a(i)
          else
             cmat%values(i,j)=((sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)))**-1
          end if
       end do
    end do

    close(55) 

    !! la matriz se imprime en pantalla
	 write(*,*) "Imprimiendo la matriz cmat"
   call Matrix_show(cmat)
		
    !! almacena la matriz
		!! preguntar como funciona y si es necesario?
	 	!	 call Matrix_writeToFile (cmat, unit=56)

    !do i=1,np
    !write(56,120) (amat(i,j), j=1,np)
    !end do


  end subroutine CosmoTools_Cmatrix
	
  !!----------------------subroutines------------------------------
!  subroutine CosmoTools_matrixinversion()
!		call Matrix
!
!  end subroutine CosmoTools_matrixinversion()


end module CosmoTools_
