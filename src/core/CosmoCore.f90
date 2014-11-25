module CosmoCore_
  use Units_
  use Matrix_
  use CONTROL_
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
  
	subroutine CosmoCore_constructor(surface,cmatinv)
    implicit none
    integer :: n

    type(surfaceSegment), intent(inout) :: surface
    type(Matrix), intent(inout) :: cmatinv

    call CosmoCore_caller(surface)
    call CosmoCore_lines(surface)
    call CosmoCore_Filler(surface)
    ! call CosmoTools_Cmatrix(surface,cmatinv)

  end subroutine CosmoCore_constructor
  
	!----------------------subroutines------------------------------

  subroutine CosmoCore_caller()
    ! subrutina que llama el ejecutable de gepol
    implicit none
    character(len=60) :: cmd

    cmd = "gepol.x < gepol.inp > gepol.out"
    call system(cmd) 
    write(*,*) "Calculating tesseras"
    cmd = "rm gepol.out"
    call system(cmd) 


  end subroutine CosmoCore_caller

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

    ! write(*,*)"estamos adentro del filler"

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

  !----------------------subroutines------------------------------
  subroutine CosmoCore_cmat(surface,cmat_inv)
    ! subroutine CosmoTools_Cmatrix(surface)
    implicit none

    integer :: i, j

    type(Matrix) :: cmat
    type(Matrix), intent(out) :: cmat_inv

    type(surfaceSegment),intent(in) :: surface

    ! llamado al constructor de matrices

    call Matrix_constructor(cmat, int(surface%sizeSurface,8), int(surface%sizeSurface,8))
    call Matrix_constructor(cmat_inv, int(surface%sizeSurface,8), int(surface%sizeSurface,8))

    do i=1,surface%sizeSurface
       do j=1,surface%sizeSurface
          if (i==j) then
             cmat%values(i,j)=3.8*surface%area(i)
          else
             cmat%values(i,j)=((sqrt((surface%xs(i)-surface%xs(j))**2+&
                  (surface%ys(i)-surface%ys(j))**2+&
                  (surface%zs(i)-surface%zs(j))**2)))**-1
          end if
       end do
    end do

    ! calculando la matriz inversa
    cmat_inv=Matrix_inverse(cmat)

  end subroutine CosmoCore_cmat
  
	!!------------------------subroutine---------------------

  subroutine CosmoCore_clasical(surface,np,cmatinv,q)
    !!esta subrutina calcula las cargas clasicas a partir de
    !!a partir de las cargas clasicas (z), sus posiciones (pz)y 
    !!las posiciones de los segmentos superficiales (ps)

    implicit none
    type(surfaceSegment), intent(in) :: surface

    type(Matrix),intent(inout) :: cmatinv

    !!matrices necesarias para el calculo
    type(Matrix) :: clasical_charge
    type(Matrix) :: aux_surface
    type(Matrix) :: v
    type(Matrix),intent(out) :: q

    !contador particulas 
    integer(8), intent(in):: np 

    !!contadores 
    integer :: i, j, k

		!!entero

		integer(8) ::segments

    !! parametro lambda segun Su-Li
    real(8) :: lambda

    !arreglo para las posiciones clasicas
    real(8), allocatable :: clasical_positions(:,:)
    real(8), allocatable :: q_clasical(:)

    logical:: verifier


    !! inicializando

    verifier=.false.
    lambda=0.0
		
		write(*,*) "surfacesize", int(surface%sizeSurface,8)

		segments=int(surface%sizeSurface,8)

		write(*,*)"segments", segments


    !asignando espacio en memoria para los parametros

    allocate(clasical_positions(np,3))
		allocate(q_clasical(segments))

    !llamando al constructor de matrices, creando matrices unidimensionales
    call Matrix_constructor(clasical_charge, np, 1)
    call Matrix_constructor(v, int(surface%sizeSurface,8), 1)
    call Matrix_constructor(q, int(surface%sizeSurface,8), 1)


    lambda=-(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC-1)/(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC+0.5)

    ! write(*,*) "esto es lambda", lambda


    do i=1,np

       !se alimenta verifier con la informacion del particle manager sobre si es
       !cuantica o clasica. En caso de ser clasica se construye el potencial clasico
       !y en el caso contrario el potencial cuantico

       verifier = ParticleManager_instance(i)%particlePtr%isQuantum

       ! write(*,*)"testest",verifier, i

       if (verifier == .false.)	then
          ! write(*,*)"particula clasica", i

          clasical_charge%values(i,1)= ParticleManager_instance(i)%particlePtr%totalCharge

          ! write(*,*)"las cargas", clasical_charge%values(i,1)

          clasical_positions(i,:)=ParticleManager_instance(i)%particlePtr%origin(:)

          ! write(*,*)"los origenes",clasical_positions(i,:)

          !Do que construye el vector potencial como el valor de la carga clasica
          !sobre la distancia euclidiana para cada una de las cargas clasicas
          !teniendo en cuenta el factor de atenuación lamda

          do j=1,surface%sizeSurface
             v%values(j,1)=clasical_charge%values(i,1)/sqrt((clasical_positions(i,1)-surface%xs(j))**2&
                  +(clasical_positions(i,2)-surface%ys(j))**2 &
                  +(clasical_positions(i,3)-surface%zs(j))**2)
             do k=1,surface%sizeSurface
                cmatinv%values(j,k)=lambda*cmatinv%values(j,k)
             end do
          end do
       end if
    end do

    ! luego se construye q
    q=Matrix_product(cmatinv,v)

    call Matrix_show(q)
    ! write(*,*)"q" 
    ! write(*,"(F15.10)") q%values(:,1)

		do i=1,segments
			q_clasical(i)=q%values(i,1)
			write(*,*)"q_clasical,q%values"
			write(*,*)q_clasical(i)
		end do

  end subroutine CosmoCore_clasical

  !----------------------subroutines------------------------------

  subroutine CosmoCore_q_builder(cmatinv, cosmo_ints, ints, q_charges)
    implicit none
    !! que estructruras se usan?
    !! son tres: una matriz (la de integrales), el inverso de la matriz c y un
    !vector donde almacenar las cargas puntuales, a la vez se necesita que le
    !pase esa información al que calcula las integrales para que funcione la
    !cosa

    real(8), allocatable, intent(inout) ::  cosmo_ints(:)
    real(8), allocatable ::  cmatinvs(:,:)
    real(8), allocatable, intent(inout) ::  q_charges(:)

    type(Matrix), intent(inout) :: cmatinv
    type(Matrix) :: q_charge

    type(Matrix) :: cosmo_pot

		real(8) :: lambda

    integer ,intent(in) :: ints

    integer :: i,j

		! primero se multiplica cmatinv por el lambda y luego por el vector

		if(allocated(q_charges)) deallocate(q_charges)
		allocate(q_charges(ints))
		
		if(allocated(cmatinvs)) deallocate(cmatinvs)
		allocate(cmatinvs(int(ints,8),int(ints,8)))

    call Matrix_constructor(q_charge, int(ints,8), 1)
    call Matrix_constructor(cosmo_pot, int(ints,8), 1)
    
		lambda=-(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC-1)/(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC+0.5)
		
    do i=1,ints
			cosmo_pot%values(i,1)=cosmo_ints(i)*-1
			do j=1,ints
       cmatinv%values(i,j)=cmatinv%values(i,j)*lambda
			 cmatinvs(i,j)=cmatinv%values(i,j)
			end do
    end do

    q_charge=Matrix_product(cmatinv,cosmo_pot)

		! call Matrix_show(q_charge)

		do i=1,ints
			q_charges(i)=q_charge%values(i,1)
		end do

		! write(*,*)q_charges(:)


		
  end subroutine CosmoCore_q_builder
  !----------------------subroutines------------------------------

  subroutine CosmoCore_q_int_builder(cosmo_pot, q_charge, ints, cosmo_pot_int)
    implicit none
    !! esta subrutina se encarga de construir los elementos matriciales de las
    !integrales asociadas a cosmo de acuerdo al paper de su-li y suma todos los
		!potenciales por par de funciones para luego ordenarlas 

    type(Matrix), intent(inout) :: cosmo_pot
    type(Matrix), intent(inout) :: q_charge
    type(Matrix) :: cosmo_pot_ints
		
		real(8) :: cosmo_pot_int
		real(8), allocatable :: cosmo_pot_ints_diito(:)

    integer :: i 
    integer ,intent(in) :: ints
		
		if(allocated(cosmo_pot_ints_diito)) deallocate(cosmo_pot_ints_diito)
		allocate(cosmo_pot_ints_diito(int(ints,8)))

    call Matrix_constructor(cosmo_pot_ints, int(ints,8), 1)

    do i = 1, ints
       cosmo_pot_ints%values(i,1)=cosmo_pot%values(i,1)*q_charge%values(i,1)
			 cosmo_pot_ints_diito(i)=cosmo_pot_ints%values(i,1)
    end do
		
		cosmo_pot_int=SUM(cosmo_pot_ints_diito)


			

  end subroutine CosmoCore_q_int_builder

  !----------------------subroutines------------------------------

  subroutine CosmoCore_q_int_staker(cosmo_pot_ints, cosmo_pot_ints_tot, cartesians, ints)
    implicit none

    type(Matrix), intent(inout) :: cosmo_pot_ints
    type(Matrix), intent(inout) :: cosmo_pot_ints_tot

    integer :: cartesians
    integer :: ints

    call Matrix_constructor(cosmo_pot_ints, int(ints,8), int(cartesians,8))

    !! como meter de nuevo todas las integrales en este arreglo para guardarlas en
    !el mugre archivo???

  end subroutine CosmoCore_q_int_staker



end module CosmoCore_
