module CosmoTools_
  use Matrix_
	use LapackInterface_
  use ParticleManager_
	use CONTROL_
	implicit none

  	type, public :: surfaceSegment
    	 real(8), allocatable :: xs(:)
			 real(8), allocatable :: ys(:)
			 real(8), allocatable :: zs(:)
			 real(8), allocatable :: area(:)
			 integer :: size
		end type surfaceSegment

  ! >Singleton
  type(surfaceSegment), public, target :: surfaceSegment_instance

contains


  !----------------------subroutines------------------------------
	subroutine CosmoTools_constructor(surface,cmatinv)
	! subroutine CosmoTools_constructor(surface)
		implicit none
		integer :: n
		
		type(surfaceSegment), intent(inout) :: surface
		type(Matrix), intent(inout) :: cmatinv

		call CosmoTools_caller()
		call CosmoTools_lines(surface)
		call CosmoTools_Filler(surface)
		! call CosmoTools_Cmatrix(surface)
		call CosmoTools_Cmatrix(surface,cmatinv)
		
	end subroutine CosmoTools_constructor
  !----------------------subroutines------------------------------

  subroutine CosmoTools_caller()
    ! subrutina que llama el ejecutable de gepol
    implicit none
    character(len=60) :: cmd

    cmd = "gepol.x < gepol.inp > gepol.out"
    call system(cmd) 
    write(*,*) "Calculating tesseras"
    		cmd = "rm gepol.out"
    call system(cmd) 


  end subroutine CosmoTools_caller

  !----------------------subroutines------------------------------
  subroutine CosmoTools_lines(surface) 
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
		surface%size=n
    write(*,*) "la superficie tiene", n, "segmentos"
    return


  end subroutine CosmoTools_lines

  !----------------------subroutines------------------------------
  subroutine CosmoTools_Filler(surface)

    ! subrutina que construye la matriz c de acuerdo a Paifeng Su
    ! a partir de el archivo generado por el cálculo usanto gepol y alimenta la
		! instancia surfaceSegment

    implicit none 

    integer :: i, j

		type(surfaceSegment), intent(inout) :: surface
    
		real(8), dimension(surface%size) :: x !segment x cordinate
    real(8), dimension(surface%size) :: y !segment y cordinate
    real(8), dimension(surface%size) :: z !segment z cordinate
    real(8), dimension(surface%size) :: a	!segment area

		write(*,*)"estamos adentro del filler"

		! llenado de surface con la información que está en vectors.vec


100 format (2X,F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
    open(55,file='vectors.vec',status='OLD') 
    read(55,100) (x(i),y(i),z(i),a(i),i=1,surface%size)

		!asignando espacio en memoria para los parametros

		allocate(surface%xs(surface%size))
		allocate(surface%ys(surface%size))
		allocate(surface%zs(surface%size))
		allocate(surface%area(surface%size))

		! write(*,*)"tipo superficie"
		!! llenando surface con la informacion leida

		do i=1,surface%size        
			surface%xs(i)=x(i)/AMSTRONG
			surface%ys(i)=y(i)/AMSTRONG
			surface%zs(i)=z(i)/AMSTRONG
			surface%area(i)=a(i)/((AMSTRONG)**2)
		! write(*,*)surface%xs(i),surface%ys(i),surface%zs(i),surface%area(i)
		end do
    
  end subroutine CosmoTools_Filler

  !!------------------------subroutine---------------------

	subroutine CosmoTools_Cmatrix(surface,cmatinv)
	! subroutine CosmoTools_Cmatrix(surface)
    implicit none

    integer :: i, j

    type(Matrix) :: cmat
		type(Matrix), intent(out) :: cmatinv
		type(Matrix) :: unity

		type(surfaceSegment),intent(in) :: surface

    ! llamado al constructor de matrices

    call Matrix_constructor(cmat, int(surface%size,8), int(surface%size,8))
    call Matrix_constructor(cmatinv, int(surface%size,8), int(surface%size,8))
    call Matrix_constructor(unity, int(surface%size,8), int(surface%size,8))

		do i=1,surface%size
       do j=1,surface%size
          if (i==j) then
             cmat%values(i,j)=3.8*surface%area(i)
          else
             cmat%values(i,j)=((sqrt((surface%xs(i)-surface%xs(j))**2+&
						 (surface%ys(i)-surface%ys(j))**2+&
						 (surface%zs(i)-surface%zs(j))**2)))**-1
          end if
       end do
    end do

    close(55) 
	! calculando la matriz inversa
	cmatinv=Matrix_inverse(cmat)
	write(*,*)"esta es la cmatinv incial"
	call Matrix_show(cmatinv)
	! verificando que la inversa esté ok
 	! unity=Matrix_product(cmatinv,cmat)	
	! call Matrix_show(unity)
	

  end subroutine CosmoTools_Cmatrix
	
	
  !!------------------------subroutine---------------------

	subroutine CosmoTools_clasical(surface,np,cmatinv)
	!!esta subrutina calcula las cargas clasicas a partir de
	!!a partir de las cargas clasicas (z), sus posiciones (pz)y 
	!!las posiciones de los segmentos superficiales (ps)

	implicit none
	type(surfaceSegment), intent(in) :: surface
	type(Matrix),intent(inout) :: cmatinv
	type(Matrix) :: clasical_charge
	! type(Matrix) :: clasical_positions
	
	type(Matrix) :: v
	type(Matrix) :: q

	

	integer(8), intent(in):: np 
	!!contador para particulas clasicas
	integer :: i, j, k
	!! parametro lambda segun Su-Li
	real(8) :: lambda


	! real(8), allocatable :: clasical_charge(:)
	real(8), allocatable :: clasical_positions(:,:)
	! !!vector que almacena los potenciales clasicos
	! real(8), allocatable :: v(:)
	! !!vector que almacena las cargas discretizadas de las particulas clasicas
	! real(8), allocatable :: q(:)


	logical:: verifier

	
	!! inicializando

	verifier=.false.
	lambda=0.0

	
	!asignando espacio en memoria para los parametros

	! allocate(clasical_charge(np))
	! allocate(v(surface%size))
	allocate(clasical_positions(np,3))
	! allocate(q(surface%size))
	
	!llamando al constructor de matrices, creando matrices unidimensionales
	call Matrix_constructor(clasical_charge, np, 1)
	call Matrix_constructor(v, int(surface%size,8), 1)
	call Matrix_constructor(q, int(surface%size,8), int(surface%size,8))
	

	lambda=(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC-1)/(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC+0.5)

	write(*,*) "esto es lambda", lambda
	
	
	! write(*,*) "estamos dentro de clasical"
	do i=1,np
	!se alimenta verifier con la informacion del particle manager sobre si es
	!cuantica o clasica. En caso de ser clasica se construye el potencial clasico
    verifier = ParticleManager_instance(i)%particlePtr%isQuantum

		! write(*,*)"testest",verifier, i

		if (verifier == .false.)	then

			! write(*,*)"particula clasica"
		
			clasical_charge%values(i,1)= ParticleManager_instance(i)%particlePtr%totalCharge
		
			! write(*,*)"las cargas",clasical_charge(i)

			clasical_positions(i,:)=ParticleManager_instance(i)%particlePtr%origin(:)
			! write(*,*)"los origenes",clasical_positions(i,:)

			!Do que construye el vector potencial como el valor de la carga clasica
			!sobre la distancia euclidiana para cada una de las cargas clasicas
			
			do j=1,surface%size
		 		v%values(j,1)=clasical_charge%values(i,1)/sqrt((clasical_positions(i,1)-surface%xs(j))**2&
			 	+(clasical_positions(i,2)-surface%ys(j))**2 &
				+(clasical_positions(i,3)-surface%zs(j))**2)
			!!primero se multiplica por el parametro lamda
				do k=1,surface%size
					cmatinv%values(j,k)=lambda*cmatinv%values(j,k)* - 1.0
				end do
			end do
			!! luego se construye q muajaja
			q=Matrix_product(cmatinv,v)
		end if
 	!para determinar las cargas se necesita multiplicar ese potencial total por la
 	!matriz c inversa, creería yo que es similiar al hacer una transformación
 	!unitaria, y la constante dialectrica. Es necesario llamar al inverso
 	!y a la constante dialectrica, la cual es determinada por el usuario.
 	!debe ser un cosmo_instance%dialectric_constant().

	!!calcula el producto lambda parameter /times cmatinv v(j)

	end do
	
	call Matrix_show(q)
			
	end subroutine CosmoTools_clasical

  !!------------------------subroutine---------------------



end module CosmoTools_
