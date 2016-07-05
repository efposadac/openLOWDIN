!!******************************************************************************
!!  This code is part of LOWDIN Quantum chemistry package                 
!!  
!!  this program has been developed under direction of:
!!
!!  Prof. A REYES' Lab. Universidad Nacional de Colombia
!!    http://www.qcc.unal.edu.co
!!  Prof. R. FLORES' Lab. Universidad de Guadalajara
!!    http://www.cucei.udg.mx/~robertof
!!
!!    Todos los derechos reservados, 2013
!!
!!******************************************************************************

module ReadIntegrals_
  use CONTROL_
  use Vector_
  implicit none


contains

  subroutine ReadIntegrals_intraSpecies(nameOfSpecies, integrals)
    implicit none
    character(*) :: nameOfSpecies
    real(8), allocatable, target :: integrals(:)

    integer :: index
    integer :: status, i

    real(8) :: int_value
    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: p(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: q(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid
    integer :: reclen
    logical :: disk = .false.

    if(.not. allocated(integrals)) disk = .true.
    
    if( disk ) then
       inquire(iolength=reclen) int_value
       open(unit=50,FILE=trim(nameOfSpecies)//".dints",ACCESS="direct",FORM="Unformatted",RECL=reclen, STATUS="unknown")
    end if

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, p, q, r, s, integral, i, index)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    if ( trim(nameOfSpecies) == "E-BETA" ) then
       open( UNIT=unitid,FILE=trim(fileid)//trim("E-ALPHA")//".ints", status='old',access='stream', form='Unformatted')
    else 
       open( unit=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//".ints", status='old',access='stream', form='Unformatted')
    end if


    loadintegrals : do

       read(UNIT=unitid, iostat=status) p(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            q(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)


       do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
          if( p(i) == -1 ) exit loadintegrals


          index = ReadIntegrals_index4Intra(int(p(i), 4), int(q(i), 4), int(r(i), 4), int(s(i), 4))

          if (disk) then
            write(50, rec=index) integral(i)
          else
             integrals(index) = integral(i)
          endif

       end do

    end do loadintegrals

    close (unitid)

    !$OMP END PARALLEL
    
    if( disk ) close(50)

  end subroutine ReadIntegrals_intraSpecies


  subroutine ReadIntegrals_interSpecies(nameOfSpecies, nameOfOtherSpecies, w, integrals)

    implicit none
    character(*) :: nameOfSpecies, nameOfOtherSpecies
    integer :: w
    real(8), allocatable, target :: integrals(:)

    integer :: index
    integer :: i

    real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: p(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: q(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
    integer :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)

    !! OpenMP related variables
    character(50) :: fileid
    integer :: nthreads
    integer :: threadid
    integer :: unitid


    !! Read integrals

    !$OMP PARALLEL private(fileid, nthreads, threadid, unitid, p, q, r, s, integral, i, index)
    nthreads = OMP_GET_NUM_THREADS()
    threadid =  OMP_GET_THREAD_NUM()
    unitid = 40 + threadid

    write(fileid,*) threadid
    fileid = trim(adjustl(fileid))

    !! open file for integrals
    open(UNIT=unitid,FILE=trim(fileid)//trim(nameOfSpecies)//"."//trim(nameOfOtherSpecies)//".ints", &
         STATUS='OLD', ACCESS='stream', FORM='Unformatted')

    loadintegrals : do

       read(unitid)   p(1:CONTROL_instance%INTEGRAL_STACK_SIZE), q(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
            integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)

       do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE

          if (p(i) == -1) exit loadintegrals

          index = ReadIntegrals_index4Inter(int(p(i), 4), int(q(i), 4), int(r(i), 4), int(s(i), 4), w)

          integrals(index) = integral(i)

       end do

    end do loadintegrals

    close (unitid)

    !$OMP END PARALLEL

  end subroutine ReadIntegrals_interSpecies



  function ReadIntegrals_index2(i, j) result(output)
    implicit none
    integer :: i, j
    integer :: output

    if(i > j) then
       output = i * (i + 1) / 2 + j
    else
       output = j * (j + 1) / 2 + i
    end if

  end function ReadIntegrals_index2

  function ReadIntegrals_index4Intra(i, j, k, l) result(output)
    implicit none
    integer :: i, j, k, l
    integer :: ii, jj, kk, ll
    integer :: output

    integer ij, kl

    ii = i - 1
    jj = j - 1
    kk = k - 1
    ll = l - 1

    ij = ReadIntegrals_index2(ii, jj)
    kl = ReadIntegrals_index2(kk, ll)

    output = ReadIntegrals_index2(ij, kl) + 1

  end function ReadIntegrals_index4Intra

  function ReadIntegrals_index4Inter(i, j, k, l, w) result(output)
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

    ij = ReadIntegrals_index2(ii, jj)
    kl = ReadIntegrals_index2(kk, ll)

    output = ij * w + kl + 1

  end function ReadIntegrals_index4Inter

end module ReadIntegrals_
