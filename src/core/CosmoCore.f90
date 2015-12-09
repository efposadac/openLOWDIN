module CosmoCore_
  use Units_
  use Matrix_
  use CONTROL_
  use MolecularSystem_
  use String_
  use Particle_
  use ParticleManager_


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


    call CosmoCore_caller()
    call CosmoCore_lines(surface)
    call CosmoCore_Filler(surface)
    call CosmoCore_cmat(surface,cmatinv)

  end subroutine CosmoCore_constructor

  !----------------------subroutines------------------------------

  subroutine CosmoCore_caller()
    ! subrutina que llama el ejecutable de gepol
    implicit none
    character(len=60) :: cmd


    cmd = "gepol.x < gepol.inp > gepol.out"
    call system(cmd) 
    ! cmd = "rm gepol.out"
    call system(cmd)

    ! write(*,*)"generada gepol surface"


  end subroutine CosmoCore_caller

  !----------------------subroutines------------------------------
  subroutine CosmoCore_lines(surface) 
    !subrutina que cuenta las lineas en el archivo vectors.out

    ! surface segment atributes
    implicit none 

    integer :: n
    character(len=60):: cmd

    type(surfaceSegment), intent(inout) :: surface


    ! cmd = "cat *.sup | grep '[^ ]' | wc -l > nlines.txt"
    cmd = "cat vectors.vec | grep '[^ ]' | wc -l > nlines.txt"
    call system(cmd) 
    open(1,file='nlines.txt')
    read(1,*) n
    close(1)
    cmd = 'rm nlines.txt'
    call system(cmd)
    surface%sizeSurface=n
    return
    ! write(*,*)"superficie segmentos",n


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

    allocate(surface%xs(surface%sizeSurface))
    allocate(surface%ys(surface%sizeSurface))
    allocate(surface%zs(surface%sizeSurface))
    allocate(surface%area(surface%sizeSurface))

    ! Reading surface from .vec
100 format (2X,F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
    open(55,file='vectors.vec',status='unknown') 
    read(55,100) (x(i),y(i),z(i),a(i),i=1,surface%sizeSurface)
    close(55)


    do i=1,surface%sizeSurface        
       surface%xs(i)=x(i)/AMSTRONG
       surface%ys(i)=y(i)/AMSTRONG
       surface%zs(i)=z(i)/AMSTRONG
       surface%area(i)=a(i)/((AMSTRONG)**2)
    end do

    !Reading surface from .sup
    !
    ! 100 format (F10.7,2X,F10.7,2X,F10.7,2X,F10.7)
    ! 		open(unit=55, file=trim(CONTROL_instance%INPUT_FILE)//"sup", status='old',	action='read') 
    !     read(55,*) (a(i),x(i),y(i),z(i),i=1,surface%sizeSurface)
    !
    ! do i=1,surface%sizeSurface        
    !    surface%xs(i)=x(i)
    !    surface%ys(i)=y(i)
    !    surface%zs(i)=z(i)
    !    surface%area(i)=a(i)
    ! end do
    ! close(55)
    !
    ! write(*,*)"tipo superficie"
    !! llenando surface con la informacion leida
    ! 	 write(*,*)"como lee los numeros"
    !
    !! gepol matrix
    ! write(*,*)"surface%sizeSurface",surface%sizeSurface
    ! write(*,*)"datos leidos"




  end subroutine CosmoCore_Filler

  !----------------------subroutines------------------------------

  subroutine CosmoCore_cmat(surface,cmat_inv)
    implicit none

    integer :: i, j

    type(Matrix) :: cmat
    type(Matrix), intent(out) :: cmat_inv

    type(surfaceSegment),intent(in) :: surface

    ! llamado al constructor de matrices

    call Matrix_constructor(cmat, int(surface%sizeSurface,8), int(surface%sizeSurface,8))
    call Matrix_constructor(cmat_inv, int(surface%sizeSurface,8), int(surface%sizeSurface,8))

    do i=1,surface%sizeSurface
       cmat%values(i,i)=1.07*sqrt(4*3.14159265/surface%area(i))

       do j=i+1,surface%sizeSurface
          cmat%values(i,j)=1/((sqrt((surface%xs(i)-surface%xs(j))**2+&
               (surface%ys(i)-surface%ys(j))**2+&
               (surface%zs(i)-surface%zs(j))**2)))
          cmat%values(j,i)=cmat%values(i,j)
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
    !! q_N

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

    integer(8) :: segments
    real(8) :: q_verifier

    !! parametro lambda segun Su-Li
    real(8) :: lambda

    !arreglo para las posiciones clasicas
    real(8), allocatable :: clasical_positions(:,:)
    real(8), allocatable :: q_clasical(:)
    type(Matrix):: cmatinv_aux

    logical:: verifier


    !! inicializando

    verifier=.false.
    lambda=0.0
    q_verifier=0.0

    ! write(*,*) "surfacesize", int(surface%sizeSurface,8)

    segments=int(surface%sizeSurface,8)

    ! write(*,*)"segments", segments


    !asignando espacio en memoria para los parametros

    allocate(clasical_positions(np,3))
    allocate(q_clasical(segments))

    !llamando al constructor de matrices, creando matrices unidimensionales
    call Matrix_constructor(clasical_charge, np, 1_8)
    call Matrix_constructor(v, int(surface%sizeSurface,8), 1_8)
    call Matrix_constructor(q, int(surface%sizeSurface,8), 1_8)
    call Matrix_constructor(cmatinv_aux, int(segments,8), int(segments,8))

    write(*,*)"Constante dialectrica = ", CONTROL_instance%COSMO_SOLVENT_DIALECTRIC


    lambda=-(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC-1)/(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC)

    write(*,*)"lambda= ",lambda

    !
    open(unit=77, file="cosmo.clasical", status="unknown",form="unformatted")


    do i=1,np

       !se alimenta verifier con la informacion del particle manager sobre si es
       !cuantica o clasica. En caso de ser clasica se construye el potencial clasico
       !y en el caso contrario el potencial cuantico

       verifier = ParticleManager_instance(i)%particlePtr%isQuantum

       ! write(*,*)"testest",verifier, i


       if (verifier .eqv. .false.)	then
          ! write(*,*)"particula clasica", i


          clasical_charge%values(i,1)= ParticleManager_instance(i)%particlePtr%totalCharge

          ! write(*,*)"las cargas", clasical_charge%values(i,1)

          clasical_positions(i,:)=ParticleManager_instance(i)%particlePtr%origin(:)

          ! write(*,*)"los origenes",clasical_positions(i,:)

          !Do que construye el vector potencial como el valor de la carga clasica
          !sobre la distancia euclidiana para cada una de las cargas clasicas
          !teniendo en cuenta el factor de atenuación lamda


          do j=1,segments

             v%values(j,1)=v%values(j,1)+(clasical_charge%values(i,1)/sqrt((clasical_positions(i,1)-surface%xs(j))**2&
                  +(clasical_positions(i,2)-surface%ys(j))**2 &
                  +(clasical_positions(i,3)-surface%zs(j))**2))
          end do


       end if
    end do

    ! write(*,*)"potential"
    ! call Matrix_show(v)


    do k=1,segments
       do j=1,segments
          cmatinv_aux%values(j,k)=cmatinv%values(j,k)
       end do
    end do


    do i = 1, segments
       q_clasical(i) = 0.0
       do j = 1, segments
          q_clasical(i) =q_clasical(i)+cmatinv_aux%values(i,j)*v%values(j,1)*lambda
       end do
    end do


    do j=1,segments
       q_verifier=q_verifier+q_clasical(j)
    end do

    write(77) q_clasical

    close(77)

    write(*,*) "Cosmo Classical Charges value"
    write(*,*) q_verifier

    ! write(*,*) "cmatinv clasical"
    ! call Matrix_show(cmatinv_aux)

  end subroutine CosmoCore_clasical

  !----------------------subroutines------------------------------

  subroutine CosmoCore_q_builder(cmatinv, cosmo_ints, ints, q_charges,specieid)
    implicit none

    !! quantum charges
    !! que estructruras se usan?
    !! son tres: una matriz (la de integrales), el inverso de la matriz c y un
    !vector donde almacenar las cargas puntuales, a la vez se necesita que le
    !pase esa información al que calcula las integrales para que funcione la
    !cosa


    real(8), allocatable, intent(inout) ::  cosmo_ints(:)
    real(8), allocatable, intent(inout) ::  q_charges(:)

    type(Matrix), intent(inout) :: cmatinv
    type(Matrix) :: q_charge
    type(Matrix) ::  cmatinvs

    type(Matrix) :: cosmo_pot

    real(8) :: lambda

    integer :: specieid
    integer :: charge


    character(30) :: specieName

    integer ,intent(in) :: ints
    !!cantidad de segmentos

    integer :: i,j

    ! primero se multiplica cmatinv por el lambda y luego por el vector


    if(allocated(q_charges)) deallocate(q_charges)
    allocate(q_charges(ints))

    call Matrix_constructor(cmatinvs, int(ints,8), int(ints,8))
    call Matrix_constructor(q_charge, int(ints,8), 1_8)
    call Matrix_constructor(cosmo_pot, int(ints,8), 1_8)

    specieName=MolecularSystem_getNameOfSpecie(specieid)

    charge=MolecularSystem_getCharge(MolecularSystem_getSpecieID(specieName))

    lambda=-(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC-1)/(CONTROL_instance%COSMO_SOLVENT_DIALECTRIC)

    do i=1,ints

       cosmo_pot%values(i,1)=cosmo_ints(i)
       ! cosmo_pot%values(i,1)=cosmo_ints(i)*charge
       do j=1,ints
          cmatinvs%values(i,j)=cmatinv%values(i,j)*lambda
       end do
    end do


    q_charge=Matrix_product(cmatinvs,cosmo_pot)


    do i=1,ints
       q_charges(i)=q_charge%values(i,1)
    end do

    ! write(78,*)q_charges(:)

    ! write(*,*)"quantum charges",specieName
    ! write(*,*)q_charges(:)


  end subroutine CosmoCore_q_builder
  !----------------------subroutines------------------------------

  subroutine CosmoCore_q_int_builder(integrals_file,charges_file,surface,charges,integrals,f_aux,g_aux,labels_aux)
    !!calcula los productos punto

    implicit none

    character(100), intent(in):: integrals_file,charges_file
    integer,intent(in) :: f_aux
    integer,intent(in) :: g_aux
    integer, allocatable,intent(in),optional :: labels_aux(:)

    integer :: surface, charges, integrals

    real(8), allocatable :: a_mat(:,:)
    real(8), allocatable :: ints_mat(:,:)
    real(8), allocatable :: ints_mat_aux(:,:)


    real(8), allocatable :: cosmo_int(:)
    integer :: i,j,k,l,m,n
    integer :: ii,g,h,hh,jj,mm

    ! write(*,*)"integrals times charges",(integrals*charges),"integrals ",	integrals_file, "charges ", charges_file

    ! if(allocated(cosmo_int)) deallocate(cosmo_int)
    ! allocate (cosmo_int(integrals*charges))

    if(allocated(a_mat)) deallocate(a_mat)
    allocate (a_mat(surface,charges))

    if(allocated(ints_mat)) deallocate(ints_mat)
    allocate (ints_mat(surface,integrals))

    open(unit=90, file=trim(integrals_file), status='old', form="unformatted") 
    open(unit=100, file=trim(charges_file), status='old', form="unformatted")
    open(unit=115, file='cosmo_int_new', status='new', form="unformatted")


    !!lectura de los archivos


    do n=1,integrals
       read(90)(ints_mat(i,n),i=1,surface)
    end do


    do n=1,charges
       read(100)(a_mat(i,n),i=1,surface)
    end do

    !!calculo del producto punto

    cosmo_int(:)=0.0_8

		!!sumatoria sobre los i segmentos
    m=1
    do n=1,integrals
       do k=1,charges
          write(115) dot_product(ints_mat(:,n),a_mat(:,k))
          m=m+1
       end do
    end do

    close(unit=115)


    close(unit=90)
    close(unit=100)
    open(unit=115, file='cosmo_int_new', status='old', form="unformatted")


    if (f_aux == g_aux) then


       if(trim(charges_file)=="cosmo.clasical") then


          allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(specieID = f_aux), MolecularSystem_getTotalNumberOfContractions(specieID = f_aux)))
          ii = 0
          mm = 1
          do g = 1, size(MolecularSystem_instance%species(f_aux)%particles)
             do h = 1, size(MolecularSystem_instance%species(f_aux)%particles(g)%basis%contraction)

                hh = h
                ii = ii + 1
                jj = ii - 1

                do i = g, size(MolecularSystem_instance%species(f_aux)%particles)
                   do j = hh, size(MolecularSystem_instance%species(f_aux)%particles(i)%basis%contraction)

                      jj = jj + 1


                      !!saving integrals on Matrix

                      do k = labels_aux(ii), labels_aux(ii) + (MolecularSystem_instance%species(f_aux)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                         do l = labels_aux(jj), labels_aux(jj) + (MolecularSystem_instance%species(f_aux)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)

                            ! ints_mat_aux(k, l) = cosmo_int(mm)
                            read(115) ints_mat_aux(k, l) 

                            ints_mat_aux(l, k) = ints_mat_aux(k, l)
                            mm = mm + 1

                         end do
                      end do

                   end do
                   hh = 1
                end do

             end do
          end do



          write(*,"(A,I6,A,A,A)")" Stored ",size(ints_mat_aux,DIM=1)**2," Quantum potential vs clasical charges ",trim(MolecularSystem_instance%species(f_aux)%name),&
               " in file cosmo.opints"

          ! write(*,*)"quantum potential",trim(MolecularSystem_instance%species(f_aux)%name)
          ! write(*,*)ints_mat(:,:)
          !
          ! write(*,*)"clasical charges"
          ! write(*,*)a_mat(:,:)

          write(40) int(size(ints_mat_aux),8)
          write(40) ints_mat_aux

       else


          open(unit=110, file=trim(MolecularSystem_instance%species(f_aux)%name)//"_qq.inn", status='unknown', form="unformatted")
          write(110) m-1
          do n=1,integrals
             do k=1,charges
                write(110) dot_product(ints_mat(:,n),a_mat(:,k))
             end do
          end do
          ! write(110) cosmo_int(:)
          close(unit=110)

          write(*,*)"same specie inner product",MolecularSystem_instance%species(f_aux)%name,m-1

          ! ints_mat(:,:)=1.0_8

          ! m=1
          ! do n=1,integrals
          !    do k=1,charges
          !       cosmo_int(m)=dot_product(ints_mat(:,n),a_mat(:,k))
          !       m=m+1
          !    end do
          ! end do
          !
	  !
          ! open(unit=110, file=trim(MolecularSystem_instance%species(f_aux)%name)//"_qq.chr", status='unknown', form="unformatted")
          ! write(110) m-1
          ! write(110) cosmo_int(:)
          ! close(unit=110)
          ! write(*,*)"charge file ordenado :",trim(MolecularSystem_instance%species(f_aux)%name)//"_qq.chr"
          !
       end if

    else	
       ! write(*,*)trim(MolecularSystem_instance%species(f_aux)%name)//trim(MolecularSystem_instance%species(g_aux)%name)

       ! write(*,*)"(ints_mat(:,n)"
       ! write(*,*)ints_mat(:,:)
       !
       ! write(*,*)"a_mat(:,k))"
       ! write(*,*)a_mat(:,:)
       !
			 write(*,*)"integrals", integrals, " charges ", charges
       open(unit=110, file=trim(MolecularSystem_instance%species(f_aux)%name)//trim(MolecularSystem_instance%species(g_aux)%name)//"_qq.cup", status='unknown', form="unformatted")
       write(110) m-1
       do n=1,integrals
          do k=1,charges
             write(110) dot_product(ints_mat(:,n),a_mat(:,k))
          end do
       end do
       ! write(110) cosmo_int(:)
       close(unit=110)
       write(*,*)"other species inner product :",trim(MolecularSystem_instance%species(f_aux)%name),"charges ,",trim(MolecularSystem_instance%species(g_aux)%name)," potentials",m-1


    end if

    close(unit=115,status='delete')

  end subroutine CosmoCore_q_int_builder

  !----------------------subroutines------------------------------

  subroutine CosmoCore_nucleiPotentialNucleiCharges(surface,output)

    type(surfaceSegment), intent(in) :: surface

    integer(8) :: np 
    integer:: segments,j,i

    type(Matrix) :: clasical_charge

    real(8), allocatable :: q_clasical(:)
    real(8), allocatable :: clasical_positions(:,:)

    real(8), intent(out) :: output

    logical:: verifier

    output=0.0_8
    verifier=.false.
    np=MolecularSystem_instance%numberOfParticles
    segments=int(surface%sizeSurface,8)

    allocate(q_clasical(segments))
    allocate(clasical_positions(np,3))

    open(unit=77, file="cosmo.clasical", status="unknown",form="unformatted")
    read(77)(q_clasical(i),i=1,segments)

    close(77)


    call Matrix_constructor(clasical_charge, np, int(1,8))


    do i=1,np

       verifier = ParticleManager_instance(i)%particlePtr%isQuantum

       if (verifier .eqv. .false.)	then


          clasical_charge%values(i,1)= ParticleManager_instance(i)%particlePtr%totalCharge
          clasical_positions(i,:)=ParticleManager_instance(i)%particlePtr%origin(:)

          do j=1,segments
             output=output+(clasical_charge%values(i,1)*q_clasical(j)/sqrt((clasical_positions(i,1)-surface%xs(j))**2&
                  +(clasical_positions(i,2)-surface%ys(j))**2 &
                  +(clasical_positions(i,3)-surface%zs(j))**2))
          end do


       end if
    end do

    output=0.5_8*output

  end subroutine CosmoCore_nucleiPotentialNucleiCharges

  !----------------------subroutines------------------------------

  subroutine CosmoCore_nucleiPotentialQuantumCharges(surface_aux,charges_file,charges,labels_aux,f_aux)

    type(surfaceSegment), intent(in) :: surface_aux
    integer, intent(in):: charges
    integer, intent(in):: f_aux
    integer, allocatable, intent(in) :: labels_aux(:)

    integer(8) :: np 
    integer:: segments,j,i,k,n
    integer:: ii,g,h,hh,jj,l,m

    type(Matrix) :: clasical_charge

    real(8), allocatable :: a_mat(:,:)
    real(8), allocatable :: clasical_positions(:,:)
    real(8), allocatable :: cosmo_int(:)
    real(8), allocatable :: ints_mat_aux(:,:)

    character(100), intent(in):: charges_file


    np=MolecularSystem_instance%numberOfPointCharges
    segments=int(surface_aux%sizeSurface,8)


    allocate(cosmo_int(charges))
    allocate(a_mat(segments,charges))
    allocate(clasical_positions(np,3))
    allocate(ints_mat_aux(MolecularSystem_getTotalNumberOfContractions(specieID = f_aux), MolecularSystem_getTotalNumberOfContractions(specieID = f_aux)))

    open(unit=100, file=trim(charges_file), status='old', form="unformatted")

    do n=1,charges
       read(100)(a_mat(i,n),i=1,segments)
    end do

    close(unit=100)

    call Matrix_constructor(clasical_charge, np, int(1,8))


    cosmo_int(:)=0.0_8

    do k=1,charges


       do i=1,np

          clasical_charge%values(i,1)= MolecularSystem_instance%pointCharges(i)%charge
          clasical_positions(i,:)=MolecularSystem_instance%pointCharges(i)%origin(:)

          do j=1,segments
             cosmo_int(k)=cosmo_int(k)+(clasical_charge%values(i,1)*a_mat(j,k)/sqrt((clasical_positions(i,1)-surface_aux%xs(j))**2&
                  +(clasical_positions(i,2)-surface_aux%ys(j))**2 &
                  +(clasical_positions(i,3)-surface_aux%zs(j))**2))
          end do
       end do
    end do

    m = 0
    ii = 0
    do g = 1, size(MolecularSystem_instance%species(f_aux)%particles)
       do h = 1, size(MolecularSystem_instance%species(f_aux)%particles(g)%basis%contraction)

          hh = h
          ii = ii + 1
          jj = ii - 1

          do i = g, size(MolecularSystem_instance%species(f_aux)%particles)
             do j = hh, size(MolecularSystem_instance%species(f_aux)%particles(i)%basis%contraction)

                jj = jj + 1


                !!saving integrals on Matrix
                do k = labels_aux(ii), labels_aux(ii) + (MolecularSystem_instance%species(f_aux)%particles(g)%basis%contraction(h)%numCartesianOrbital - 1)
                   do l = labels_aux(jj), labels_aux(jj) + (MolecularSystem_instance%species(f_aux)%particles(i)%basis%contraction(j)%numCartesianOrbital - 1)
                      m = m + 1

                      ints_mat_aux(k, l) = cosmo_int(m)
                      ints_mat_aux(l, k) = ints_mat_aux(k, l)

                   end do
                end do

             end do
             hh = 1
          end do

       end do
    end do
    write(*,"(A,I6,A,A,A)")" Stored ",size(ints_mat_aux,DIM=1)**2," Clasical potential vs quantum charges ",trim(MolecularSystem_instance%species(f_aux)%name),&
         " in file lowdin.opints"
    write(40) int(size(ints_mat_aux),8)
    write(40) ints_mat_aux


  end subroutine CosmoCore_nucleiPotentialQuantumCharges

end module CosmoCore_
