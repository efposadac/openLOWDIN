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
module LapackInterface_
  implicit none
  
  
  interface LapackInterface
     
     !>
     !! @brief Calcula todos los valores propios y opcionalmente los eigenvectores
     !! de una matriz real y simetrica. Los eigenvetores estan ortonormalizados.
     !! Para detalles consultar documentacion correspondiente.
     subroutine dsyev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
       character :: JOBZ, UPLO
       integer :: INFO, LDA, LWORK, N
       real(8) :: A(LDA,*), W(*), WORK(*)
     end subroutine dsyev

     !>
     !! @brief LInear
     subroutine dgesv(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       integer :: N, NRHS, LDA, LDB, INFO
       real(8) :: A(LDA,*), B(*)
       integer :: IPIV(*)
     end subroutine dgesv

     !>
     !! @brief Calcula todos los valores propios y opcionalmente los eigenvectores
     !! por la izquieda o derecha de una matriz cuadrada real y no-simetrica.
     !! Los eigenvectores se normalizan para tener norma euclidiana igual a 1
     !! y la   mas grande componente real. Para detalles consultar documentacion
     !! correspondiente.
     subroutine dgeev(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR,&
          WORK, LWORK, INFO)
       character :: JOBVL, JOBVR
       integer :: INFO, LDA, LDVL, LDVR, LWORK, N
       real(8) :: A(LDA,*), VL(LDVL,*), VR(LDVR,*), WI(*), WORK(*), WR(*)
     end subroutine dgeev
     
     !>
     !! @brief Factoriza una matriz de N x M mediante factorizacion LU empleando pivoteo parcial
     subroutine dgetrf( M, N, A, LDA, IPIV, INFO )
       integer :: INFO, LDA, M, N
       integer :: IPIV( * )
       real(8) :: A(LDA, *)
     end subroutine dgetrf
		
     !>
     !! @brief Calcula la inversa de una matriz previamente factorizada por pivoteo parcial
     subroutine  dgetri( N, A, LDA, IPIV, WORK, LWORK, INFO )
       integer :: INFO, LDA, LWORK, N
       integer :: IPIV(*)
       real(8) :: A(LDA,*), WORK(*)
     end subroutine dgetri
     
     !>
     !! @brief Calcula la descomposicion en valore simples de una matriz real de MxN
     subroutine  dgesvd( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,&
          WORK, LWORK, INFO )
       character :: JOBU, JOBVT
       integer :: INFO, LDA, LDU, LDVT, LWORK, M, N
       real(8) :: A(  LDA,  *  ),  S( * ), U( LDU, * ), VT(LDVT, * ), WORK( * )
     end subroutine dgesvd

     !>                                                                         
     !! @brief Calcula la descomposicion en valore simples de una matriz real de MxN 
     subroutine dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
       character :: TRANSA, TRANSB
       integer :: M, N, K
       integer :: LDA, LDB, LDC
       real(8) :: A(LDA,*), B(LDB,*), C(LDC,*)
       real(8) :: ALPHA, BETA
     end subroutine dgemm
     
  end interface
  
end module LapackInterface_
