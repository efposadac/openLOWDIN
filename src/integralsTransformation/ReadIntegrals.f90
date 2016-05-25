module ReadIntegrals_
    use CONTROL_
    use Vector_
    implicit none


contains
    
    subroutine ReadIntegrals_intraSpecies(nameOfSpecie, integrals, nproc)
        implicit none
        character(*) :: nameOfSpecie
        real(8), allocatable, target :: integrals(:)
        integer :: nproc

        integer :: ifile, unit, status, i
        character(50) :: sfile

        real(8) :: integral(CONTROL_instance%INTEGRAL_STACK_SIZE)
        integer*2 :: p(CONTROL_instance%INTEGRAL_STACK_SIZE)
        integer*2 :: q(CONTROL_instance%INTEGRAL_STACK_SIZE)
        integer*2 :: r(CONTROL_instance%INTEGRAL_STACK_SIZE)
        integer*2 :: s(CONTROL_instance%INTEGRAL_STACK_SIZE)

        print*, "nproc", nproc

        !! Read integrals
        do ifile = 1, nproc

            write(sfile,*) ifile
            sfile = trim(adjustl(sfile))
            unit = ifile+50

            if ( trim(nameOfSpecie) == "E-BETA" ) then
                open( UNIT=unit,FILE=trim(sfile)//trim("E-ALPHA")//".ints", status='old',access='sequential', form='Unformatted')
            else 
                open( UNIT=unit,FILE=trim(sfile)//trim(nameOfSpecie)//".ints", status='old',access='sequential', form='Unformatted')
            end if

            loadintegrals : do
  
                read(UNIT=unit, iostat=status) &
                    p(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                    q(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                    r(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                    s(1:CONTROL_instance%INTEGRAL_STACK_SIZE), &
                    integral(1:CONTROL_instance%INTEGRAL_STACK_SIZE)
        
        
                do i = 1, CONTROL_instance%INTEGRAL_STACK_SIZE
                    if( p(i) == -1 ) exit loadintegrals
                    integrals(index4(p(i), q(i), r(i), s(i)) + 1) = integral(i)
                end do

            end do loadintegrals
        end do

    end subroutine ReadIntegrals_intraSpecies


    function index2(i, j) result(output)
        implicit none
        integer*2 :: i, j
        integer :: output

        if(i > j) then
            output = i * (i + 1) / 2 + j
        else
            output = j * (j + 1) / 2 + i
        end if

    end function index2

    function index4(i, j, k, l) result(output)
        implicit none
        integer*2 :: i, j, k, l
        integer :: output
        
        integer ij, kl

        i = i - 1_2
        j = j - 1_2
        k = k - 1_2
        l = l - 1_2

        ij = index2(i, j)
        kl = index2(k, l)
        output = index2(int(ij, 2), int(kl, 2))

    end function index4

end module ReadIntegrals_