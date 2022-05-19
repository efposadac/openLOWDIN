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
module JadamiluInterface_

  use, intrinsic :: iso_c_binding
  implicit none
  
  interface JadamiluInterface
     
      SUBROUTINE DPJDREVCOM ( N , A , JA , IA , EIGS , RES , X , LX , NEIG  , &
                            SIGMA , ISEARCH , NINIT , MADSPACE , ITER , &
                            TOL , SHIFT , DROPTOL , MEM , ICNTL , &
                            IJOB , NDX1 , NDX2 , IPRINT , INFO , GAP )
        implicit none

        integer(8) N , LX , NEIG , ISEARCH , NINIT , MADSPACE , INFO
        integer(8) ITER , ICNTL (5) , IJOB , NDX1 , NDX2 , IPRINT
        DOUBLE PRECISION SIGMA , TOL , SHIFT , DROPTOL , MEM , GAP
        integer(8) JA (*) , IA (*)
        DOUBLE PRECISION A (*), X(*)
        DOUBLE PRECISION EIGS (* ) , RES (*) 
      END SUBROUTINE DPJDREVCOM 

  end interface
  
end module JadamiluInterface_
