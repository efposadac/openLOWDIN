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
!! @brief Clase estatica que contiene la interface para utilizar lapack
!!
!! @author Nestor Aguirre
!!
!! <b> Fecha de creacion : </b> 2008-08-19
!!   - <tt> 2008-08-19 </tt>: Nestor Aguirre ( nfaguirrec@unal.edu.co )
!!        -# Creacion del archivo y las funciones basicas
!!           tomadas de (http://www.netlib.org/blas/)
!!   - <tt> 2011-02-11 </tt>: Fernando Posada ( efposadac@unal.edu.co )
!!        -# Adapta el modulo para su inclusion en Lowdin
module ArpackInterface_
  implicit none
  
  interface ArpackInterface
     
     subroutine dsaupd ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,IPNTR, WORKD, WORKL, LWORKL, INFO )
       character  bmat*1, which*2
       integer    ido, info, ldv, lworkl, n, ncv, nev
       Double precision  tol
       integer    iparam(*), ipntr(*)
       Double precision  resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
     end subroutine dsaupd

    subroutine DSEUPD(RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, &
                         WORKD, WORKL, LWORKL, INFO )

      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Double precision  sigma, tol
      integer    iparam(*), ipntr(*)
      logical    select(ncv)
      Double precision d(nev), resid(n), v(ldv,ncv), z(ldz, nev), workd(2*n), workl(lworkl)
    end subroutine dseupd
     
  end interface
  
end module ArpackInterface_
