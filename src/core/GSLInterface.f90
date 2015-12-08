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
!! @brief Clase estatica que contiene la interface para utilizar GSL
!!
!! @author Mauricio Rodas
!!
!! <b> Fecha de creacion : </b> 2015-09-28
!!   - <tt> 2015-09-28 </tt>: Mauricio Rodas ( jmrodasr@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
module GSLInterface_
  implicit none
  
  private
  public ::                  &
       tolow_minimize

  integer, public, parameter ::     &
       MINMETHOD_STEEPEST_DESCENT = 1, &
       MINMETHOD_FR_CG            = 2, &
       MINMETHOD_PR_CG            = 3, &
       MINMETHOD_BFGS             = 4, &
       MINMETHOD_BFGS2            = 5, &
       MINMETHOD_NMSIMPLEX        = 6


  interface tolow_minimize
     function low_minimize(method, dim, x, step, tolgrad, toldr, maxiter, f, write_iter_info, minimum)
       integer :: low_minimize
       integer, intent(in)    :: method
       integer, intent(in)    :: dim
       real(8), intent(inout) :: x
       real(8), intent(in)    :: step
       real(8), intent(in)    :: tolgrad
       real(8), intent(in)    :: toldr
       integer, intent(in)    :: maxiter
       interface
          subroutine f(n, x, val, getgrad, grad)
            integer, intent(in)    :: n
            real(8), intent(in)    :: x(n)
            real(8), intent(inout) :: val
            integer, intent(in)    :: getgrad
            real(8), intent(inout) :: grad(n)
          end subroutine f
          subroutine write_iter_info(iter, n, val, maxdr, maxgrad, x)
            integer, intent(in) :: iter
            integer, intent(in) :: n
            real(8), intent(in) :: val
            real(8), intent(in) :: maxdr
            real(8), intent(in) :: maxgrad
            real(8), intent(in) :: x(n)
          end subroutine write_iter_info
       end interface
       real(8), intent(out)   :: minimum
     end function low_minimize
  end interface tolow_minimize

end module GSLInterface_
