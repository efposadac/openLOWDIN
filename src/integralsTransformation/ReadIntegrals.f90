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
                    integrals(ReadIntegrals_index4(int(p(i), 4), int(q(i), 4), int(r(i), 4), int(s(i), 4))) = integral(i)
                end do

            end do loadintegrals
        end do

    end subroutine ReadIntegrals_intraSpecies


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

    function ReadIntegrals_index4(i, j, k, l) result(output)
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

    end function ReadIntegrals_index4

end module ReadIntegrals_