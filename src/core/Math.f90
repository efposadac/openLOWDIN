!!******************************************************************************
!!	This code is part of LOWDIN Quantum chemistry package                 
!!	
!!	this program has been developed under direction of:
!!
!!	Prof. A REYES' Lab. Universidad Nacional de Colombia
!!		http://www.qcc.unal.edu.co
!!	Prof. R. FLORES' Lab. Universidad de Guadajara
!!		http://www.cucei.udg.mx/~robertof
!!
!!		Todos los derechos reservados, 2013
!!
!!******************************************************************************

!>
!! @brief   Modulo para implementacion de funciones matematicas auxiliares
!!
!! Este modulo implementa funciones matematicas auxiliares requeridas para la
!! ejecucion de ciertos procedimientos dentro las clases definidas.
!!
!! Dentro de las funciones implementadas esta la funcion auxiliar definida como:
!!
!!  \f[ F_m(U)= \int_{0}^{1} {t^{2m}e^{-Tt^2}}\,dt \f]
!!
!! cuyo calculo  se realiza atraves de la expansion en una serie de Taylor
!! de <i> funciones de error </i>. El valor umbral para el cual la funcion
!! auxiliar (gama incompleta) alcanza el limite:
!!
!! \f[ \lim_{U\to 0}{F_m(U)} = \frac{1}{2m+1} \f]
!!
!! se fija por defecto en 1E-9 y puede manipularse en el modulo FControl.
!!
!! @author Sergio A. Gonzalez Monico
!!
!! <b> Fecha de creacion : </b> 2006-06-05
!!
!! <b> Historial de modificaciones: </b>
!!
!!   - <tt> 2007-01-06 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Propuso estandar de codificacion.
!!   - <tt> 2007-05-15 </tt>: Sergio A. Gonzalez M. ( sagonzalezm@unal.edu.co )
!!        -# Se adapto al estandar de codificacion propuesto.
!!   - <tt> 2010-03-10 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Implementa funcion auxiliar para ordenes superiores
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el módulo para su inclusion en Lowdin
!!
module Math_
  use Exception_
  use, intrinsic :: iso_c_binding
  implicit none


  real(16), parameter :: Math_PI_32 = 5.568327996831707845284817982118835702014_16
  real(8), parameter :: Math_PI = 3.141592653589793238_8
  real(16), parameter :: Math_PI_MOD = 3.141592653589793238_16
  real(8), parameter :: Math_SQRT_PI = Math_PI ** 0.5_8
  real(8), parameter :: Zero = 0.0_8
  real(8), parameter :: Math_NaN = Z'7FFFFFFFFFFFFFFF'

  !>
  !! Gamma function related parameters
  integer, parameter :: maxfac = 24

  !> Global parameters
  real(kind = 8), parameter  :: teps = 1.0e-13_8

  !>Maximum n value of the tabulated f_n(t) function values
  integer, save :: current_nmax = -1

  !>f_n(t) table
  real(kind = 8), dimension(:,:), allocatable, save :: ftable

  !>Inverse Factorial function ifac
  real(kind=8), parameter, dimension (0:maxfac) :: ifac = (/	&
       0.10000000000000000000e+01_8, 0.10000000000000000000e+01_8, 0.50000000000000000000e+00_8,	&
       0.16666666666666666667e+00_8, 0.41666666666666666667e-01_8, 0.83333333333333333333e-02_8,	&
       0.13888888888888888889e-02_8, 0.19841269841269841270e-03_8, 0.24801587301587301587e-04_8,	&
       0.27557319223985890653e-05_8, 0.27557319223985890653e-06_8, 0.25052108385441718775e-07_8,	&
       0.20876756987868098979e-08_8, 0.16059043836821614599e-09_8, 0.11470745597729724714e-10_8,	&
       0.76471637318198164759e-12_8, 0.47794773323873852974e-13_8, 0.28114572543455207632e-14_8,	&
       0.15619206968586226462e-15_8, 0.82206352466243297170e-17_8, 0.41103176233121648585e-18_8,	&
       0.19572941063391261231e-19_8, 0.88967913924505732867e-21_8, 0.38681701706306840377e-22_8,	&
       0.16117375710961183490e-23_8/)

!
!
! @brief llama a la función de dawson de c
  interface 
     function gsl_dawson (xi) result(output) bind(C, name="gsl_sf_dawson")
       import 
       real (kind=c_double), value :: xi
       real (kind=c_double) :: output
     end function gsl_dawson
  end interface 

  public :: &
       Math_Factorial, &
       Math_kroneckerDelta, &
       Math_isEven, &
       Math_isOdd, &
       Math_crossProduct, &
       Math_roundNumber, &
       Math_numberRepresentation, &
       Math_fgamma0, &
       Math_fgamma1, &
       Math_combinatorial, &
       Math_Erfi, &
       Math_Log_Erfi, &
!       sum_six, &
       Math_Dawson, &
       sum_pairwise, &
       sum_kahan,&
       sort


contains

function Math_Dawson(xx) result(daw)

!    This routine evaluates Dawson's integral,
!
!      F(x) = exp ( - x * x ) * Integral ( 0 <= t <= x ) exp ( t * t ) dt
!
!    for a real argument x.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by William Cody.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    William Cody, Kathleen Paciorek, Henry Thacher,
!    Chebyshev Approximations for Dawson's Integral,
!    Mathematics of Computation,
!    Volume 24, Number 109, January 1970, pages 171-178.
!
!  Parameters:
!
!    Input, real ( kind = 16 ) XX, the argument of the function.
!
!    Output, real ( kind = 16 ) DAW, the value of the function.
!  Implementada para resolver las integrales de LJ
implicit none
  
  real ( kind = 16 ) daw
  real ( kind = 16 ) frac
  integer ( kind = 16 ) i
  real ( kind = 16 ) one225
  real ( kind = 16 ) p1(10)
  real ( kind = 16 ) p2(10)
  real ( kind = 16 ) p3(10)
  real ( kind = 16 ) p4(10)
  real ( kind = 16 ) q1(10)
  real ( kind = 16 ) q2(9)
  real ( kind = 16 ) q3(9)
  real ( kind = 16 ) q4(9)
  real ( kind = 16 ) six25
  real ( kind = 16 ) sump
  real ( kind = 16 ) sumq
  real ( kind = 16 ) two5
  real ( kind = 16 ) w2
  real ( kind = 16 ) x
  real ( kind = 16 ) xlarge
  real ( kind = 16 ) xmax
  real ( kind = 16 ) xsmall
  real ( kind = 16 ) xx
  real ( kind = 16 ) y
!
!  Mathematical constants.
!
  data six25 / 6.25D+00 /
  data one225 / 12.25d0 /
  data two5 / 25.0d0 /
  !
  !  Machine-dependent constants
  !
  data xsmall /1.05d-08/
  data xlarge /9.49d+07/
  data xmax /2.24d+307/
!
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!
  data p1/-2.69020398788704782410d-12, 4.18572065374337710778d-10, &
          -1.34848304455939419963d-08, 9.28264872583444852976d-07, &
          -1.23877783329049120592d-05, 4.07205792429155826266d-04, &
          -2.84388121441008500446d-03, 4.70139022887204722217d-02, &
          -1.38868086253931995101d-01, 1.00000000000000000004d+00/
  data q1/ 1.71257170854690554214d-10, 1.19266846372297253797d-08, &
           4.32287827678631772231d-07, 1.03867633767414421898d-05, &
           1.78910965284246249340d-04, 2.26061077235076703171d-03, &
           2.07422774641447644725d-02, 1.32212955897210128811d-01, &
           5.27798580412734677256d-01, 1.00000000000000000000d+00/

!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [2.5, 3.5)
!
  data p2/-1.70953804700855494930d+00,-3.79258977271042880786d+01, &
           2.61935631268825992835d+01, 1.25808703738951251885d+01, &
          -2.27571829525075891337d+01, 4.56604250725163310122d+00, &
          -7.33080089896402870750d+00, 4.65842087940015295573d+01, &
          -1.73717177843672791149d+01, 5.00260183622027967838d-01/
  data q2/ 1.82180093313514478378d+00, 1.10067081034515532891d+03, &
          -7.08465686676573000364d+00, 4.53642111102577727153d+02, &
           4.06209742218935689922d+01, 3.02890110610122663923d+02, &
           1.70641269745236227356d+02, 9.51190923960381458747d+02, &
           2.06522691539642105009d-01/

!
!  Coefficients for R(9,9) approximation in J-fraction form
!  for  x in [3.5, 5.0]
!
  data p3/-4.55169503255094815112d+00,-1.86647123338493852582d+01, &
           -7.36315669126830526754d+00,-6.68407240337696756838d+01, &
           4.84507265081491452130d+01, 2.69790586735467649969d+01, &
           -3.35044149820592449072d+01, 7.50964459838919612289d+00, &
           -1.48432341823343965307d+00, 4.99999810924858824981d-01/
  data q3/ 4.47820908025971749852d+01, 9.98607198039452081913d+01, &
           1.40238373126149385228d+01, 3.48817758822286353588d+03, &
          -9.18871385293215873406d+00, 1.24018500009917163023d+03, &
          -6.88024952504512254535d+01,-2.31251575385145143070d+00, &
           2.50041492369922381761d-01/
!
!  Coefficients for R(9,9) approximation in
!  J-fraction form
!  for 5.0 < |x|.
!
  data p4/-8.11753647558432685797d+00,-3.84043882477454453430d+01,&
          -2.23787669028751886675d+01,-2.88301992467056105854d+01, &
          -5.99085540418222002197d+00,-1.13867365736066102577d+01, &
          -6.52828727526980741590d+00,-4.50002293000355585708d+00, &
          -2.50000000088955834952d+00, 5.00000000000000488400d-01/
  data q4/ 2.69382300417238816428d+02, 5.04198958742465752861d+01, &
           6.11539671480115846173d+01, 2.08210246935564547889d+02, &
           1.97325365692316183531d+01,-1.22097010558934838708d+01, &
          -6.99732735041547247161d+00,-2.49999970104184464568d+00, &
           7.49999999999027092188d-01/

  x = xx

  if ( xlarge < abs ( x ) ) then

    if ( abs ( x ) <= xmax ) then
     daw = 0.5D+00 / x
    else
     daw = 0.0D+00
    end if

  else if ( abs ( x ) < xsmall ) then
    daw = x
  else
    y = x * x
!
!  ABS(X) < 2.5.
!
    if ( y < six25 ) then
        sump = p1(1)
        sumq = q1(1)
        do i = 2, 10
           sump = sump * y + p1(i)
           sumq = sumq * y + q1(i)
        end do
        daw = x * sump / sumq

!
!  2.5 <= ABS(X) < 3.5.
!
    else if ( y < one225 ) then

    frac = 0.0D+00
    do i = 1, 9
        frac = q2(i) / ( p2(i) + y + frac )
    end do

    daw = ( p2(10) + frac ) / x
!
!  3.5 <= ABS(X) < 5.0.
!
  else if ( y < two5 ) then
    
      frac = 0.0D+00
      do i = 1, 9
        frac = q3(i) / ( p3(i) + y + frac )
      end do

      daw = ( p3(10) + frac ) / x

  else

!
!  5.0 <= ABS(X) <= XLARGE.
!
      w2 = 1.0D+00 / x / x

      frac = 0.0D+00
      do i = 1, 9
         frac = q4(i) / ( p4(i) + y + frac )
      end do
      frac = p4(10) + frac

      daw = ( 0.5D+00 + 0.5D+00 * w2 * frac ) / x

  end if

  end if


end function Math_Dawson



function sum_pairwise(x) result(s)
    implicit none
    
    real(16), dimension(:), intent(in) :: x
    real(16) :: s
    
    integer(4), parameter :: n = 1024
    integer(4) :: l, i, m

    l =size(x)
    if (1<=n) then
        s = x(1)
        do i=2, l
         s = s + x(i)
        end do
    else
        m = 1/2
        s = sum_pairwise(x(1:m)) + sum_pairwise(x(m+1:1))
    end if
end function sum_pairwise


function sum_kahan(x) result(s)
   implicit none

   real(kind=16), dimension(:), intent(in) :: x

   real(16) :: s, c = 0.0_16, y, t
   integer(4) :: i

    s = 0.0_16
    do i =1, size(x)
     y = x(i) - c
     t = s + y
     c = (t - s) - y
     s = t
    end do
end function sum_kahan


!-----------------------------------------------------------------------------
! Given:  A    An array of length max, holding data from index 1 to index num.
!         max  The maximum number of reals that can be put into array A.
!         num  The number of data items in array A.
! Task:   To sort that part of the data in array A that is from index 1 to num.
!         This data is put into ascending order using an insertion sort.
! Return: A    The sorted array.
!-----------------------------------------------------------------------------
subroutine  sort(A, num, max) 
   implicit none

    integer(8), intent(in):: num
    integer(8), intent(in):: max
    real(16), dimension(max), intent(inout):: A(:)
    integer(8) i, k, indexofmin
    real(16) minvalue

   do k = 1, num - 1
   ! find the smallest value in a(k) through a(num).
      indexofmin = k
      minvalue = A(k)
            
      do i = k + 1, num
         if (A(i) < minvalue) then
             indexofmin = i
             minvalue = A(i)
         end if
      end do

! swap data items in array a at indices k and indexofmin, if
! indexofmin /= k.
         if (indexofmin /= k) then
            A(indexofmin) = A(k)
            A(k) = minvalue
         end if
   end do

   
end subroutine sort




!function sum_six(n,buff) result(output)  
!    implicit none
!   
!    integer, intent(in) :: n
!    real(16), intent(in) :: buff(1:n)
!
!    real(16) :: c, y, t
!    real(16) :: output 
!    integer :: i
!
!    output = 0.0_16
!    c = 0.0_16
!
!    do i = 1, n
!      y = buff(i) - c
!      t = output + y
!      c = (t - output) - y
!      output =  real(t, 16)
!    end do
!
!     output = buff(1)
!     c = 0.0_16
!
!     do i = 2, size(buff)
!       t = output + buff(i)
!       if (abs(output) >= abs(buff(i))) then
!         c = c + (output - t) + buff(i)
!       else
!         c = c + (buff(i) - t) + output
!       end if
!
!       output = t
!
!     end do
!
!     output = output + c
!
!end function sum_six

  !>
  !! Calcula el factorial de cualquier entero entre -1 e &#8734;
  !!
  !! @param N Numero entre -1 e &#8734;
  !!
  !! @return factorial_ Factorial de N
  function Math_factorial(N,skip,offset) result( output )
    implicit none
    integer(8), intent(in) :: N
    integer,optional :: skip
    integer,optional :: offset
    integer(8):: output
    integer(8)::i
    integer:: auxSkip
    integer::auxOffset

    auxSkip=1
    if(present(skip)) auxSkip=skip

    auxOffset=1
    if(present(offset)) auxOffset=offset

    !! Devuelve el factorial de -1
    if ( N == -1 ) then
       output = 1_8
       return
       !! Calcula el factorial de un entero positivo
    else if ( N >= 0 ) then
       output = 1_8
       do i=auxOffset, N, auxSkip
          output = output* i
       end do
    else
       output = 0_8
       return
    end if
  end function Math_factorial


  !>
  !! Obtiene la combinatoria de dos número (N M)
  !!
  !! Junio02/17
  !!
  function Math_combinatorial(N, M) result( output)
  implicit none

  integer(8), intent(in) :: N, M
  integer(8) :: output
  integer(8) :: num, den, diff
  
   num = Math_factorial(N)
   den = Math_factorial(M)
   diff = Math_factorial(N-M)
   output = num / (den*diff)
  end function Math_combinatorial

!>
  !! Obtiene el valor de la función Erfi de un número 
  !!
  !! Junio02/17
  !!
!  function Math_Erfi(N) result( output)
!  implicit none
!
!  integer(8), intent(in) :: N
!  integer(8) :: output
!  
!  end function Math_Erfi



  !>
  !!  Obtiene el producto cruz entre un par de vectores de tres componetes
  !!
  !! @param vectorA Primer vector.
  !! @param vectorB Segundo vector.
  !!
  !! @return crossProduct Vector resultante
  function  Math_crossProduct(vectorA,vectorB) result( output )
    implicit none
    real(8), intent(in):: vectorA(3)
    real(8), intent(in):: vectorB(3)
    real(8):: output(3)
    output(1) = vectorA(2) * vectorB(3) - vectorA(3) * vectorB(2)
    output(2) = vectorA(3) * vectorB(1) - vectorA(1) * vectorB(3)
    output(3) = vectorA(1) * vectorB(2) - vectorA(2) * vectorB(1)
  end function Math_crossProduct

  !>
  !! Calcula el delta de Kronecker para los numeros i, j
  !!
  !! @return delta de Kronecker  
  function Math_kroneckerDelta( i, j ) result( output )
    implicit none
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(8):: output

    output = 1.0_8
    if ( i /= j ) output = 0.0_8

  end function Math_kroneckerDelta

  !>
  !! Devuelve verdadero si el numero es par
  !!
  function Math_isEven( i ) result( output )
    implicit none
    integer, intent(in) :: i
    logical:: output

    output = .not.btest(0,i)

  end function Math_isEven

  !>
  !! Devuelve verdadero si el numero es impar
  function Math_isOdd( i ) result( output )
    implicit none
    integer, intent(in) :: i
    logical:: output

    output = btest(0,i)

  end function Math_isOdd

  !>
  !! @brief Redondea un real de doble presicion
  function Math_roundNumber(value,roundValue ) result( output )
    implicit none
    real(8):: value
    real(8):: roundValue
    real(8):: output

    real(8):: auxValue
    integer:: i

    i = int(value)
    auxValue = value - i
    if ( value > 0.0 ) then
       auxValue = aint( auxValue * roundValue + 0.5 ) / roundValue
    else
       auxValue = aint( auxValue * roundValue - 0.5 ) / roundValue
    end if
    output = i + auxValue

  end function Math_roundNumber


  subroutine Math_numberRepresentation(value, mantisse, decimalExponent)
    implicit  none

    real(8), intent(in) :: value
    real(8), intent(inout) :: mantisse
    integer, intent(inout) :: decimalExponent

    mantisse=abs(value)
    decimalExponent=0
    do while(mantisse<1.0)
       decimalExponent=decimalExponent+1
       mantisse=abs(value)*(10.0**decimalExponent)
    end do

  end subroutine Math_numberRepresentation

  !> @brief   Build a table of F_n(t) values in the range tmin <= t <= tmax
  !>          with a stepsize of tdelta up to n equal to nmax.
  !> @date    11.01.1999
  !> @par Parameters
  !>       - nmax  : Maximum n value of F_n(t).
  !>       - tdelta: Difference between two consecutive t abcissas (step size).
  !>       - tmax  : Maximum t value.
  !>       - tmin  : Minimum t value.
  !> @version 1.0
  subroutine create_md_ftable(nmax,tmin,tmax,tdelta)

    integer, intent(in)                      :: nmax
    real(kind=8), intent(in)                :: tmin
    real(kind=8), intent(in)                :: tmax
    real(kind=8), intent(in)                :: tdelta

    integer:: istat
    integer:: itab
    integer:: itabmax
    integer:: 		&
         itabmin
    integer:: n
    real(kind=8):: t

    n = nmax + 6

    itabmin = floor(tmin/tdelta)
    itabmax = ceiling((tmax - tmin)/tdelta)

    allocate (ftable(0:n,itabmin:itabmax),stat=istat)
    ftable = 0.0_8

    !   *** fill table ***
    do itab=itabmin,itabmax
       t = real(itab,8)*tdelta
       ftable(0:n,itab) = fgamma_ref(n,t)
    end do

    !   *** save initialization status ***
    current_nmax = nmax

  end subroutine create_md_ftable

  !> @brief   Deallocate the table of F_n(t) values.
  !> @date    24.05.2004
  !> @version 1.0
  subroutine deallocate_md_ftable()

    integer:: istat

    if (current_nmax > -1) then

       deallocate (ftable,stat=istat)
       ! 		if (istat /= 0) then
       ! 			call stop_memory(routinen,modulen,__line__,"ftable")
       ! 		end if

       current_nmax = -1

    end if

  end subroutine deallocate_md_ftable

  !> @brief   Calculation of the incomplete Gamma function F(t) for multicenter
  !>          integrals over Gaussian functions. f returns a vector with all
  !>          F_n(t) values for 0 <= n <= nmax.
  !> @date    08.01.1999,
  !> @par History
  !>          09.06.1999, MK : Changed from a FUNCTION to a SUBROUTINE
  !> @par Literature
  !>       L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
  !> @par Parameters
  !>       - f   : The incomplete Gamma function F_n(t).
  !>       - nmax: Maximum n value of F_n(t).
  !>       - t   : Argument of the incomplete Gamma function.
  !>       - kmax: Maximum number of iterations.
  !>       - expt: Exponential term in the upward recursion of F_n(t).
  !> @version 1.0
  subroutine Math_fgamma0(nmax,t,f)

    integer, intent(in) :: nmax
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(0:nmax),intent(out) :: f

    integer:: itab
    integer:: k
    integer:: n
    real(kind=8):: expt
    real(kind=8):: g
    real(kind=8):: tdelta
    real(kind=8):: tmp
    real(kind=8):: ttab

    !   *** calculate f(t) ***
    if (t < teps) then

       !     *** special cases: t = 0 ***
       do n=0,nmax
          f(n) = 1.0_8/real(2*n + 1,8)
       end do

    else if (t <= 12.0_8) then

       !     *** 0 < t < 12 -> taylor expansion ***

       tdelta = 0.1_8

       !     *** pretabulation of the f_n(t) function ***
       !     *** for the taylor series expansion      ***

       if (nmax > current_nmax) then
          call init_md_ftable(nmax)
       end if


       itab = nint(t/tdelta)
       ttab = real(itab,8)*tdelta

       f(nmax) = ftable(nmax,itab)

       tmp = 1.0_8
       do k=1,6
          tmp = tmp * (ttab - t)
          f(nmax) = f(nmax) + ftable(nmax+k,itab)*tmp*ifac(k)
       end do

       expt = exp(-t)

       !     *** use the downward recursion relation to ***
       !     *** generate the remaining f_n(t) values   ***

       do n=nmax-1,0,-1
          f(n) = (2.0_8*t*f(n+1) + expt)/real(2*n + 1,8)
       end do

    else

       !     *** t > 12 ***

       if (t <= 15.0_8) then

          !       *** 12 < t <= 15 -> four term polynom expansion ***

          g = 0.4999489092_8 - 0.2473631686_8/t +&
               0.321180909_8/t**2 - 0.3811559346_8/t**3
          f(0) = 0.5_8*sqrt(Math_PI/t) - g*exp(-t)/t

       else if (t <= 18.0_8) then

          !       *** 15 < t <= 18 -> three term polynom expansion ***

          g = 0.4998436875_8 - 0.24249438_8/t + 0.24642845_8/t**2
          f(0) = 0.5_8*sqrt(Math_PI/t) - g*exp(-t)/t

       else if (t <= 24.0_8) then

          !       *** 18 < t <= 24 -> two term polynom expansion ***

          g = 0.499093162_8 - 0.2152832_8/t
          f(0) = 0.5_8*sqrt(Math_PI/t) - g*exp(-t)/t

       else if (t <= 30.0_8) then

          !       *** 24 < t <= 30 -> one term polynom expansion ***

          g = 0.49_8
          f(0) = 0.5_8*sqrt(Math_PI/t) - g*exp(-t)/t

       else

          !       *** t > 30 -> asymptotic formula ***

          f(0) = 0.5_8*sqrt(Math_PI/t)

       end if

       if (t > real(2*nmax + 36,8)) then
          expt = 0.0_8
       else
          expt = exp(-t)
       end if

       !     *** use the upward recursion relation to ***
       !     *** generate the remaining f_n(t) values ***

       do n=1,nmax
          f(n) = 0.5_8*(real(2*n - 1,8)*f(n-1) - expt)/t
       end do

    end if

  end subroutine Math_fgamma0

  !> @brief   Calculation of the incomplete Gamma function F(t) for multicenter
  !>          integrals over Gaussian functions. f returns a vector with all
  !>          F_n(t) values for 0 <= n <= nmax.
  !> @date    08.01.1999
  !> @par Literature
  !>       L. E. McMurchie, E. R. Davidson, J. Comp. Phys. 26, 218 (1978)
  !> @par Parameters
  !>       - f   : The incomplete Gamma function F_n(t).
  !>       - nmax: Maximum n value of F_n(t).
  !>       - t   : Argument of the incomplete Gamma function.
  !> @version 1.0
  subroutine Math_fgamma1(nmax,t,f)

    integer, intent(in)                      :: nmax
    real(kind=8), dimension(:), intent(in)  :: t
    real(kind=8), 		&
         dimension(size(t, 1), 0:nmax), 		&
         intent(out)                            :: f

    integer:: i
    integer:: itab
    integer:: k
    integer:: n
    real(kind=8):: expt
    real(kind=8):: g
    real(kind=8):: tdelta
    real(kind=8):: tmp
    real(kind=8):: ttab

    do i=1,size(t,1)

       !     *** calculate f(t) ***

       if (t(i) < teps) then

          !       *** special cases: t = 0 ***

          do n=0,nmax
             f(i,n) = 1.0_8/real(2*n + 1,8)
          end do

       else if (t(i) <= 12.0_8) then

          !       *** 0 < t < 12 -> taylor expansion ***

          tdelta = 0.1_8

          !       *** pretabulation of the f_n(t) function ***
          !       *** for the taylor series expansion      ***

          if (nmax > current_nmax) then
             call init_md_ftable(nmax)
          end if

          itab = nint(t(i)/tdelta)
          ttab = real(itab,8)*tdelta

          f(i,nmax) = ftable(nmax,itab)

          tmp = 1.0_8
          do k=1,6
             tmp = tmp * (ttab - t(i))
             f(i,nmax) = f(i,nmax) + ftable(nmax+k,itab)*tmp*ifac(k)
          end do

          expt = exp(-t(i))

          !       *** use the downward recursion relation to ***
          !       *** generate the remaining f_n(t) values   ***

          do n=nmax-1,0,-1
             f(i,n) = (2.0_8*t(i)*f(i,n+1) + expt)/real(2*n + 1,8)
          end do

       else

          !       *** t > 12 ***

          if (t(i) <= 15.0_8) then

             !         *** 12 < t <= 15 -> four term polynom expansion ***

             g = 0.4999489092_8 - 0.2473631686_8/t(i) +&
                  0.321180909_8/t(i)**2 - 0.3811559346_8/t(i)**3
             f(i,0) = 0.5_8*sqrt(Math_PI/t(i)) - g*exp(-t(i))/t(i)

          else if (t(i) <= 18.0_8) then

             !         *** 15 < t <= 18 -> three term polynom expansion ***

             g = 0.4998436875_8 - 0.24249438_8/t(i) + 0.24642845_8/t(i)**2
             f(i,0) = 0.5_8*sqrt(Math_PI/t(i)) - g*exp(-t(i))/t(i)

          else if (t(i) <= 24.0_8) then

             !         *** 18 < t <= 24 -> two term polynom expansion ***

             g = 0.499093162_8 - 0.2152832_8/t(i)
             f(i,0) = 0.5_8*sqrt(Math_PI/t(i)) - g*exp(-t(i))/t(i)

          else if (t(i) <= 30.0_8) then

             !         *** 24 < t <= 30 -> one term polynom expansion ***

             g = 0.49_8
             f(i,0) = 0.5_8*sqrt(Math_PI/t(i)) - g*exp(-t(i))/t(i)

          else

             !         *** t > 30 -> asymptotic formula ***

             f(i,0) = 0.5_8*sqrt(Math_PI/t(i))

          end if

          if (t(i) > real(2*nmax + 36,8)) then
             expt = 0.0_8
          else
             expt = exp(-t(i))
          end if

          !       *** use the upward recursion relation to ***
          !       *** generate the remaining f_n(t) values ***

          do n=1,nmax
             f(i,n) = 0.5_8*(real(2*n - 1,8)*f(i,n-1) - expt)/t(i)
          end do

       end if

    end do

  end subroutine Math_fgamma1

  !> @brief   Calculation of the incomplete Gamma function F_n(t) using a
  !>          spherical Bessel function expansion. fgamma_ref returns a
  !>          vector with all F_n(t) values for 0 <= n <= nmax.
  !>          For t values greater than 50 an asymptotic formula is used.
  !>          This function is expected to return accurate F_n(t) values
  !>          for any combination of n and t, but the calculation is slow
  !>          and therefore the function may only be used for a pretabulation
  !>          of F_n(t) values or for reference calculations.
  !> @date    07.01.1999
  !> @par Literature
  !>        F. E. Harris, Int. J. Quant. Chem. 23, 1469 (1983)
  !> @par Parameters
  !>       - expt   : Exponential term in the downward recursion of F_n(t).
  !>       - factor : Prefactor of the Bessel function expansion.
  !>       - nmax   : Maximum n value of F_n(t).
  !>       - p      : Product of the Bessel function quotients.
  !>       - r      : Quotients of the Bessel functions.
  !>       - sumterm: One term in the sum over products of Bessel functions.
  !>       - t      : Argument of the incomplete Gamma function.
  !> @version 1.0
  function fgamma_ref(nmax,t) result(f)

    integer, intent(in)                      :: nmax
    real(kind=8), intent(in)                :: t
    real(kind=8), dimension(0:nmax)         :: f

    integer, parameter                       :: kmax = 50
    real(kind=8), parameter                 :: eps = epsilon(0.0_8)

    integer:: j
    integer:: k
    integer:: n
    real(kind=8):: expt
    real(kind=8):: factor
    real(kind=8):: p
    real(kind=8):: sumterm
    real(kind=8):: 		&
         sumtot
    real(kind=8):: term
    real(kind=8), dimension(kmax+10)        :: r

    !   ------------------------------------------------------------------
    !   *** initialization ***

    f(:) = 0.0_8

    if (t < teps) then

       !     *** special case: t = 0 => analytic expression ***

       do n=0,nmax
          f(n) = 1.0_8/real(2*n + 1,8)
       end do

    else if (t <= 50.0_8) then

       !     *** initialize ratios of bessel functions ***

       r(kmax+10) = 0.0_8

       do j=kmax+9,1,-1
          r(j) = -t/(real(4*j + 2,8) - t*r(j+1))
       end do

       factor = 2.0_8*sinh(0.5_8*t)*exp(-0.5_8*t)/t

       do n=0,nmax

          !       *** initialize iteration ***

          sumtot = factor/real(2*n + 1,8)
          term = 1.0_8

          !       *** begin the summation and recursion ***

          do k=1,kmax

             term = term*real(2*n - 2*k + 1,8)/real(2*n + 2*k + 1,8)

             !         *** product of bessel function quotients ***

             p = 1.0_8

             do j=1,k
                p = p*r(j)
             end do

             sumterm = factor*term*p*real(2*k + 1,8)/real(2*n + 1,8)

             if (abs(sumterm) < eps) then

		!           *** iteration converged ***

                exit

             else

		!           *** add the current term to the sum and continue the iteration ***

                sumtot = sumtot + sumterm

             end if

          end do

          f(n) = sumtot

       end do

    else

       !     *** use asymptotic formula for t > 50 ***

       f(0) = 0.5_8*sqrt(Math_PI/t)

       !     *** use the upward recursion relation to ***
       !     *** generate the remaining f_n(t) values ***

       expt = exp(-t)

       do n=1,nmax
          f(n) = 0.5_8*(real(2*n - 1,8)*f(n-1) - expt)/t
       end do

    end if

  end function fgamma_ref

  !> @brief   Initalize a table of F_n(t) values in the range 0 <= t <= 12 with
  !!            a stepsize of 0.1 up to n equal to nmax for the Taylor series
  !!            expansion used by McMurchie-Davidson (MD).
  !> @date    10.06.1999
  !> @par Parameters
  !>       - nmax   : Maximum n value of F_n(t).
  !> @version 1.0
  subroutine init_md_ftable(nmax)
    integer, intent(in)                      :: nmax


    !   *** check, if the current initialization is sufficient ***

    if (nmax > current_nmax) then

       call deallocate_md_ftable()

       !     *** pretabulation of the f_n(t) function ***
       !     *** for the taylor series expansion      ***

       call create_md_ftable(nmax,0.0_8,12.0_8,0.1_8)

    end if

  end subroutine init_md_ftable

  !>
  !! @brief  Maneja excepciones de la clase
  subroutine Math_exception( typeMessage, description, debugDescription)
    implicit none
    integer:: typeMessage
    character(*):: description
    character(*):: debugDescription

    type(Exception):: ex

    call Exception_constructor( ex , typeMessage )
    call Exception_setDebugDescription( ex, debugDescription )
    call Exception_setDescription( ex, description )
    call Exception_show( ex )
    call Exception_destructor( ex )

  end subroutine Math_exception

!
!  @brief Evalua la función Erfi utilizando la función gsl_dawson de c
!
  function Math_Erfi(xi) result(output)
   
    implicit none

    real(8) :: xi, a, output
    a = gsl_dawson(xi) 
    output = 2.0_8 / sqrt(Math_PI) * exp(xi**2.0_8) * a
  end function

  function Math_Log_Erfi(xi) result(output)

    implicit none

    real(8) :: xi, a, output
    a = gsl_dawson(xi)
    output = log(2.0_8 / sqrt(Math_PI) * a) + xi**2
    
  end function



 
end module Math_
