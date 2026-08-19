!  ============================== subroutine =============================!
!=========================================================================!
subroutine mt_amux1(n,x,y,a,ja,ia)
implicit none
!
!-------------------------------------------------------------------------!
!                                                                         !
!  Y = A * X                                                              !
!  input:                                                                 !
!     n     = row dimension of A                                          !
!     x     = array of length equal to the column dimension of matrix A   !
!     a, ja, ia = input matrix in compressed sparse row format.           !
!  Output:                                                                !
!     y     = real array of length n, containing the product y=Ax         !
!                                                                         !
!-------------------------------------------------------------------------!
!
integer :: n,i,k1,k2
integer :: ja(*)
integer :: ia(*)
real*8  :: x(*)
real*8  :: a(*)
real*8, dimension(n) :: y
!  ...
!  ...
do i= 1,n
   k1 = ia(i)
   k2 = ia(i+1)-1
   y(i) = dot_product( a(k1:k2),x(ja(k1:k2)) )
enddo
!=====
!=====
return
    End subroutine
    
! - - - - - - - - - - 
subroutine resi_normP(n,nintf,x,b,a,ja,ia,res0)

use omp_lib
!USE MD_OpenMP
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,b
!real*8, dimension(nintf) :: b
real*8 temp,res0
!  ...
res0 = 0.d0

!$omp PARALLEL DO private(k1,k2,temp,j) reduction(+:res0)
do i= 1,nintf
   k1 = ia(i)
   k2 = ia(i+1)-1
   
   temp = b(i)
   do j=k1,k2
   temp = temp-a(j)*x(ja(j))
   enddo
!   temp = dot_product(a(k1:k2),x(ja(k1:k2)))
   res0 = res0 +  temp*temp    !(b(i)-temp)**2.d0
enddo
!$omp end PARALLEL DO
!=====
return
      
    end subroutine
! - - - - - - - - -- 
! - - - - - - - - - - 
subroutine resi_normP_MPI(n,nintf,x,b,a,ja,ia,res0)

!use omp_lib
!USE MD_OpenMP
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,b
!real*8, dimension(nintf) :: b
real*8 temp,res0
!  ...
res0 = 0.d0

!!$omp PARALLEL DO private(k1,k2,temp,j) reduction(+:res0)
do i= 1,nintf
   k1 = ia(i)
   k2 = ia(i+1)-1
   
   temp = b(i)
   do j=k1,k2
   temp = temp-a(j)*x(ja(j))
   enddo
!   temp = dot_product(a(k1:k2),x(ja(k1:k2)))
   res0 = res0 +  temp*temp    !(b(i)-temp)**2.d0
enddo
!!$omp end PARALLEL DO
!=====
return
      
    end subroutine
! - - - - - - - - -- 
    
subroutine resi_P(n,nintf,x,b,r,a,ja,ia)

use omp_lib
!USE MD_OpenMP
      
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,r,b
!real*8, dimension(nintf) :: b
real*8 temp
!  ...

!$omp PARALLEL DO private(k1,k2,temp,j) 
do i= 1,nintf
   k1 = ia(i)
   k2 = ia(i+1)-1
   temp = b(i)
   do j=k1,k2
    temp = temp -a(j)*x(ja(j))   
   enddo
   
!   temp = dot_product(a(k1:k2),x(ja(k1:k2)))
   r(i) = temp  !b(i)-temp
enddo
!$omp end PARALLEL DO
!=====
return
      
    end subroutine

!=========================================================================!
subroutine resi_P_MPI(n,nintf,x,b,r,a,ja,ia)

!use omp_lib
!USE MD_OpenMP
      
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,r,b
!real*8, dimension(nintf) :: b
real*8 temp
!  ...

!!$omp PARALLEL DO private(k1,k2,temp,j) 
do i= 1,nintf
   k1 = ia(i)
   k2 = ia(i+1)-1
   temp = b(i)
   do j=k1,k2
    temp = temp -a(j)*x(ja(j))   
   enddo
   
!   temp = dot_product(a(k1:k2),x(ja(k1:k2)))
   r(i) = temp  !b(i)-temp
enddo
!!$omp end PARALLEL DO
!=====
return
      
end subroutine

!=========================================================================!
subroutine mt_amux2(n,n1,nnz,x,y,a,ja,ia)

use omp_lib
!USE MD_OpenMP

implicit none
!
!-------------------------------------------------------------------------!
!                                                                         !
!  Y = A * X                                                              !
!  input:                                                                 !
!     A(n,n1)                                          !
!     x(n1) = array of length equal to the column dimension of matrix A   !
!     a, ja, ia = input matrix in compressed sparse row format.           !
!  Output:                                                                !
!     y(n)  = real array of length n, containing the product y=Ax         !
!                                                                         !
!-------------------------------------------------------------------------!
!
integer :: n,n1,nnz,i,k1,k2
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8 :: y(*),x(*)
!  ...
!  ...

!$omp PARALLEL DO private(k1,k2)
do i= 1,n
   k1 = ia(i)
   k2 = ia(i+1)-1
   y(i) = dot_product( a(k1:k2),x(ja(k1:k2)) )
enddo
!$omp end PARALLEL DO
!=====
!=====
return
    end subroutine
! - - - - - - - - - - - - -- 
!=========================================================================!
!=========================================================================!
subroutine mt_amux2_MPI(n,n1,nnz,x,y,a,ja,ia)

!use omp_lib
!USE MD_OpenMP

implicit none
!
!-------------------------------------------------------------------------!
!                                                                         !
!  Y = A * X                                                              !
!  input:                                                                 !
!     A(n,n1)                                          !
!     x(n1) = array of length equal to the column dimension of matrix A   !
!     a, ja, ia = input matrix in compressed sparse row format.           !
!  Output:                                                                !
!     y(n)  = real array of length n, containing the product y=Ax         !
!                                                                         !
!-------------------------------------------------------------------------!
!
integer :: n,n1,nnz,i,k1,k2
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8 :: y(*),x(*)
!  ...
!  ...

!!$omp PARALLEL DO private(k1,k2)
do i= 1,n
   k1 = ia(i)
   k2 = ia(i+1)-1
   y(i) = dot_product( a(k1:k2),x(ja(k1:k2)) )
enddo
!!$omp end PARALLEL DO
!=====
!=====
return
end subroutine
! - - - - - - - - - - - - -- 
!=========================================================================!
subroutine mt_amux2P(n,nintf,n1,nintf1,nnz,x,y,a,ja,ia)
implicit none
!
!-------------------------------------------------------------------------!
!                                                                         !
!  Y = A * X                                                              !
!  input:                                                                 !
!     A(n,n1)                                          !
!     x(n1) = array of length equal to the column dimension of matrix A   !
!     a, ja, ia = input matrix in compressed sparse row format.           !
!  Output:                                                                !
!     y(n)  = real array of length n, containing the product y=Ax         !
!                                                                         !
!-------------------------------------------------------------------------!
!
integer :: n,n1,nnz,i,k1,k2,nintf,nintf1,j
integer :: ja(nnz)
integer :: ia(n+1)
real*8  :: a(nnz)
real*8 :: y(n),x(n1)
real*8 temp
!  ...
!  ...
!write(*,*)'1'
!i=484935
!   k1 = ia(i)
!   k2 = ia(i+1)-1
!   write(*,*)'k1,k2',k1,k2
   
!            temp = 0.d0
!      do j=k1,k2
!          write(*,*)'a,ja,x',a(j),ja(j),x(ja(j))
!      temp = temp + a(j)*x(ja(j))
!      enddo
!      write(*,*)'y',y(i)
!      pause

do i= 1,nintf
!    write(*,*)'i=',i
   k1 = ia(i)
   k2 = ia(i+1)-1
!   IF(k1.eq.0.or.k2.eq.0.or.ja(k1).eq.0.or.ja(k2).eq.0) THEN
!       write(*,*)'err',k1,k2,ja(k1),ja(k2)
!   endif
!   if(minval(ja(k1:k2)).eq.0) then
!       write(*,*)'min',minval(ja(k1:k2))
!   endif

   
         temp = 0.d0
      do j=k1,k2
 !         if(ja(j).EQ.0.or.ja(j).GT.n1) then
 !             write(*,*)'err1',ja(j)
 !             cycle
 !         endif
          
      temp = temp + a(j)*x(ja(j))
      enddo
   y(i) = temp
!       write(*,*)'i=/',i
 !  y(i) = dot_product(a(k1:k2),x(ja(k1:k2)))
enddo
!write(*,*)'2'
!=====
!=====
return
end subroutine

! = = = = (P1-⑤에서 06_solver_pcg_ilu.f90 로부터 이식 — BiCGSTAB 본체 SpMV) = = = =
      SUBROUTINE amux0P(nintf,n,nnz,x,y,a,ja,ia)
!-----------------------------------------------------------
!     Y = A * X
!     input:
!       n     = row dimension of A
!       x     = array of length equal to the column dimension of matrix A
!       a, ja, ia = input matrix in compressed sparse row format.
!     output:
!       y     = real array of length n, containing the product y=Ax
!-------------------------------------------------------------------
      IMPLICIT NONE
!      
      INTEGER  nintf,n,nnz
	  INTEGER  ja(nnz),ia(nintf+1)
      REAL*8 a(nnz),x(n)
! 
      REAL*8 y(n)
! tmp
      INTEGER i, k
      REAL*8 tmp
	  
!
      DO i= 1,nintf
        tmp = 0.d0
        DO k=ia(i),ia(i+1)-1
          tmp = tmp + a(k)*x(ja(k))
        ENDDO
        y(i) = tmp
      ENDDO
!	  
      RETURN
      END
