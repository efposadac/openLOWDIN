module CosmoTools_
  use CosmoCore_
  use Matrix_
  use IntegralManager_
  use LapackInterface_
  use ParticleManager_
  use CONTROL_
  implicit none

contains


  !----------------------subroutines------------------------------
  subroutine CosmoTools_constructor(surface,cmatinv)
    ! subroutine CosmoTools_constructor(surface)
    implicit none
    integer :: n

    type(surfaceSegment), intent(inout) :: surface
    type(Matrix), intent(inout) :: cmatinv

    call CosmoTools_caller()
    ! call CosmoCore_lines(surface)
    ! call CosmoCore_Filler(surface)
    ! call CosmoTools_Cmatrix(surface,cmatinv)

  end subroutine CosmoTools_constructor


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

    call Matrix_constructor(cmat, int(surface%sizeSurface,8), int(surface%sizeSurface,8))
    call Matrix_constructor(cmatinv, int(surface%sizeSurface,8), int(surface%sizeSurface,8))
    call Matrix_constructor(unity, int(surface%sizeSurface,8), int(surface%sizeSurface,8))

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

    close(55) 
    ! calculando la matriz inversa
    cmatinv=Matrix_inverse(cmat)

    ! write(*,*)"esta es la cmatinv incial"
    ! do j=1,surfaceSegment%sizeSurface
    ! 	write(*,"(5(F15.10))") cmatinv%values(j,:)
    ! end do

    ! call Matrix_show(cmatinv)
    ! verificando que la inversa esté ok
    ! unity=Matrix_product(cmatinv,cmat)
    ! write(*,*)"esta es la unity"

    ! call Matrix_show(unity)
    ! do j=1,surfaceSegment%sizeSurface
    ! 	write(*,"(5(F15.10))") unity%values(j,:)
    ! end do


  end subroutine CosmoTools_Cmatrix


  !!------------------------subroutine---------------------

  subroutine CosmoTools_quantum(surface,np,cmatinv,qe)

    implicit none
    type(surfaceSegment), intent(in) :: surface
    type(Matrix),intent(in) :: cmatinv

    !!matrices necesarias para el calculo
    type(Matrix) :: quantum_charge
    type(Matrix) :: aux_surface
    type(Matrix) :: ve
    type(Matrix) :: qe

    !contador particulas 
    integer(8), intent(in):: np 
    integer :: sizeSurface_aux

    !!contadores 
    integer :: i, j, k


    !arreglo para las posiciones clasicas
    real(8), allocatable :: particle_positions(:,:)

    logical:: verifier
    !! inicializando

    verifier=.false.

    sizeSurface_aux=surface%sizeSurface

    !asignando espacio en memoria para los parametros

    allocate(particle_positions(np,3))

    !llamando al constructor de matrices, creando matrices unidimensionales

    call Matrix_constructor(quantum_charge, np, 1)
    call Matrix_constructor(ve, int(surface%sizeSurface,8), 1)
    call Matrix_constructor(qe, int(surface%sizeSurface,8), 1)
    call Matrix_constructor(aux_surface, int(surface%sizeSurface,8), 4_8)


    do k=1,surface%sizeSurface
       aux_surface%values(k,1)=surface%xs(k)
       aux_surface%values(k,2)=surface%ys(k)
       aux_surface%values(k,3)=surface%zs(k)
       aux_surface%values(k,4)=surface%area(k)
    end do



    !Do que construye el vector potencial como una integral monoelectronica
    !respecto a las posiciones de las particulas cuanticas sobre la distancia
    !euclidiana para cada una de las posiciones en las cuales están centrandas
    !las funciones bases 

    call system(" lowdin-ints.x COSMO ")

    ! Hace el llamado al integral Manager para que compute las integrales de
    ! atracción





    ! call Matrix_show(q)
    ! call Matrix_show(cmatinv)

  end subroutine CosmoTools_quantum





end module CosmoTools_
