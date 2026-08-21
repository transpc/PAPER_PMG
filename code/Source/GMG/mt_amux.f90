!  ============================== subroutine =============================!
!=========================================================================!
    
! - - - - - - - - - - 
subroutine resi_normP(n,nintf,x,b,a,ja,ia,res0)

use omp_lib
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,b
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
   res0 = res0 +  temp*temp    !(b(i)-temp)**2.d0
enddo
!$omp end PARALLEL DO
!=====
return
      
    end subroutine
! - - - - - - - - -- 
    
subroutine resi_P(n,nintf,x,b,r,a,ja,ia)

use omp_lib
      
implicit none
!
integer :: n,nintf,i,k1,k2,j
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8, dimension(n) :: x,r,b
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
   
   r(i) = temp  !b(i)-temp
enddo
!$omp end PARALLEL DO
!=====
return
      
    end subroutine

!=========================================================================!

!=========================================================================!
subroutine mt_amux2(n,n1,nnz,x,y,a,ja,ia)

use omp_lib

implicit none
!
!-------------------------------------------------------------------------!
!                                                                         !
!  input:                                                                 !
!     A(n,n1)                                          !
!     a, ja, ia = input matrix in compressed sparse row format.           !
!  Output:                                                                !
!                                                                         !
!-------------------------------------------------------------------------!
!
integer :: n,n1,nnz,i,k1,k2
integer :: ja(*)
integer :: ia(*)
real*8  :: a(*)
real*8 :: y(*),x(*)
!  ...

!$omp PARALLEL DO private(k1,k2)
do i= 1,n
   k1 = ia(i)
   k2 = ia(i+1)-1
   y(i) = dot_product( a(k1:k2),x(ja(k1:k2)) )
enddo
!$omp end PARALLEL DO
!=====
return
    end subroutine
! - - - - - - - - - - - - -- 
!=========================================================================!
!=========================================================================!
! - - - - - - - - - - - - -- 
!=========================================================================!

! = = = = (P1-⑤에서 06_solver_pcg_ilu.f90 로부터 이식 — BiCGSTAB 본체 SpMV) = = = =
      SUBROUTINE amux0P(nintf,n,nnz,x,y,a,ja,ia)
!-----------------------------------------------------------
!     input:
!       a, ja, ia = input matrix in compressed sparse row format.
!     output:
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
