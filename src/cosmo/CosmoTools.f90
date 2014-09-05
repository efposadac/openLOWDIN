module CosmoTools_
  use Matrix_
	use LapackInterface_
  use Particle_
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
		type(Matrix) :: cmatinv

    !real(8), allocatable :: cmat(:,:)
    !allocate(amat(np,np))

    !! llamado al constructor de matrices

    call Matrix_constructor(cmat, int(np,8), int(np,8))
    call Matrix_constructor(cmatinv, int(np,8), int(np,8))


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
	
	cmatinv=Matrix_inverse(cmat)

  end subroutine CosmoTools_Cmatrix
	
	
  !!------------------------subroutine---------------------

	subroutine CosmoTools_clasicalpot(n,np)
	!!esta subrutina calcula los potenciales clasicos
	!!a partir de las cargas clasicas (z), sus posiciones (pz)y 
	!!las posiciones de los segmentos superficiales (ps)

	implicit none
	integer, intent(in):: n, np
	!!contadores para las cargas clasicas y los segmentos
	integer :: i, j
	real(8) :: V, vij
	real(8) :: xs,ys,zs

	real(8), dimension(n) :: clasical
	real(8), dimension(n,3)	:: clasical_positions
	
	V=0
	vij=0

	!! construyendo las matrices y vectores con 0
	do i=1,n
		clasical(i)=0
	end do

	do i=1,n
		do j=1,3
			clasical_positions(i,j)=0
		end do
	end do
	

	!! llenando los vectores y las matrices
	!! el de las superficies se llena a partir de archivos
	!! uno se puede llenar usando particle_information
	
	if (particle%isQuantum /= ".true.")

	then

		do i=1,n
			clasical(i)=particle%totalCharge
		end do
	
		do i=1,n
			do j=1,3
				clasical_positions(i,j)=particle%origin
			end do
		end do
		
		!!leyendo el archivo y sacando la información 

100 format(F12.8,2X,F12.8,2X,F12.8)
    open(55,file='vectors.vec',status='OLD') 
    read(55,100) (xs(j),ys(j),zs(j),j=1,np)

	
		do i=1,n
			do j=1,np
				vij=clasical(i)/sqrt((clasical_positions(i,1)-xs(j))**2&
				+(clasical_positions(i,2)-ys(j))**2+&
				(clasical_positions(i,3)-zs(j))**2)
				V=V+vij
			end do
		end do

		close(55) 

		write(*,*) "clasical potential" 
	end if

	end subroutine CosmoTools_clasicalpot



end module CosmoTools_
