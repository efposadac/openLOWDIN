!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadalajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief COSMO and COSMO-APMO program.
!!        This module allows to make calculations in the COSMO-APMO framework
!! @author D. Gonzalez.
!!
!! <b> Creation date : </b> 2014-21-08
!!
!! <b> History: </b>
!!
!!   - <tt> 2008-05-25 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Creacion de modulo y procedimientos  para calculos con solvente implicito
!!
!! @warning This programs only works linked to lowdincore library, and using lowdin-ints.x and lowdin-SCF.x programs, 
!!          all those tools are provided by LOWDIN quantum chemistry package
!!
program Cosmo
  use CONTROL_
  use MolecularSystem_
	use String_
		implicit none 
		integer :: n

  !!Start time
  call Stopwatch_constructor(lowdin_stopwatch)
  call Stopwatch_start(lowdin_stopwatch)

  !!Load CONTROL Parameters
  call MolecularSystem_loadFromFile( "LOWDIN.DAT" )

  !!Load the system in lowdin.sys format
  call MolecularSystem_loadFromFile( "LOWDIN.SYS" )


		call caller()
		write(*,*)"llamado realizado a gepol"
		
		call lines(n)	
		
		write(*,*)"Se crearon ", n, "segmentos"
		
		write(*,*)"construyendo matriz C"
		
		call Cmatrix(n)
		
		write(*,*)"Finalizado"


end program

!!----------------------subroutines------------------------------

subroutine caller()
!subrutina que llama el ejecutable de gepol
implicit none
character(len=60) :: cmd

    cmd = "gepol.x < gepol.inp > gepol.out"
    call system(cmd) 
  	write(*,*) "Calculating tesseras"
!		cmd = "rm gepol.out"
    call system(cmd) 

	
end subroutine

!!----------------------subroutines------------------------------
subroutine lines(n) 
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

end subroutine

!!----------------------subroutines------------------------------
subroutine Cmatrix(np)

!subrutina que construye la matriz c de acuerdo a Paifeng Su
!a partir de el archivo generado por el cálculo usanto gepol
implicit none 

integer :: i, j, np
real(8), dimension(np) :: x !segment x cordinate
real(8), dimension(np) :: y !segment y cordinate
real(8), dimension(np) :: z !segment z cordinate
real(8), dimension(np) :: a	!segment area


real(8), allocatable :: amat(:,:)

allocate(amat(np,np))


100	format(F12.8,2X,F12.8,2X,F12.8,2X,F12.8)
120 format(F12.8, F12.8)

open(55,file='vectors.vec',status='OLD') 

read(55,100) (x(i),y(i),z(i),a(i),i=1,np)


	do i=1,np
		do j=1,np
			if (i==j) then
				amat(i,j)=3.8*a(i)
			else
				amat(i,j)=((sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)))**-1
			end if
		end do
	end do

close(55) 

!! la matriz es almacenada en el archivo amat.mat

open(56,file="amat.mat",status="replace")

	do i=1,np
		write(56,120) (amat(i,j), j=1,np)
	end do
		
close(56)

end subroutine 
