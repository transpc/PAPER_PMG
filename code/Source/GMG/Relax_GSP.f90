    
! - - - - - - - - - - - - - - - - - - - -- 
      subroutine Relax_GS(maxit,relax,n,nnz,ia,ja,ju,au,u,b)
      
      use omp_lib
!      USE MD_OpenMP    
! ---
      implicit none
! ---
      integer n,maxit,nnz
      integer i,j,k,iter,j1,j2
      integer ia(*),ja(*),ju(*)
      real*8 temp,temp1,relax
      real*8 u(*),b(*),au(*)    

! ---
!---Gauss-Seidel method -Parallel
  Do iter=1,maxit

!$omp PARALLEL DO private(j1,j2,temp,temp1,j)
   Do i=1,n
      j1 = ia(i)
	  j2 = ia(i+1)-1
      j = ju(i)
      temp1 = au(j)
! ---
     temp = dot_product(au(j1:j2),u(ja(j1:j2)))
	 u(i)=(b(i)-temp)/temp1 + u(i)
   End do 
       
!$omp end PARALLEL DO
   
  End do 
  
      return
    end
    
! 
! - - - - - - - - - - - - - - - - - - - -- 
! - - - - - - - - - - - - - - - - - - - -- 
    
! 
! - - - - - - - - - - - - - - - - - - - -- 
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
! = = = = = = = = = = = = = = = = = = = = 
SUBROUTINE Smooth_GS2(maxit,ista,iend,b,u,au,ia,ja,ju,diagr)
! ---
! note that it is for GS forward 
! maxit = 1
use omp_lib
! 
IMPLICIT NONE
! ---
INTEGER(4) maxit,ista,iend
INTEGER(4) ia(*)
INTEGER(4) ja(*)
INTEGER(4) ju(*)
REAL(8) b(*)
REAL(8) u(*)
REAL(8) au(*)
REAL(8) diagr(*)
! temp
INTEGER(4) i,j,iter,j1,j2,k
REAL(8) temp,temp1   
! ---
!---Gauss-Seidel method
!  DO iter=1,maxit
!$omp PARALLEL DO private(j1,j2,j,temp,k) 
   DO i=ista,iend
      j1 = ia(i)
	  j2 = ia(i+1)-1
!      j = ju(i)
!      temp1 = au(j)
!
      temp = b(i)
      do k=j1,j2
         j = ja(k)
         temp = temp -au(k)*u(j) 
      enddo
      
! ---
!    temp = DOT_PRODUCT( au(j1:j2),u(ja(j1:j2)) )
    u(i)=temp*diagr(i) + u(i) !(b(i)-temp)/temp1 + u(i)
   ENDDO 
     
!$omp end PARALLEL DO
!ENDDO 
  
RETURN
    END
! = = = = = = = = = = 
! = = = = = = = = = = = = = = = = = = = = 
! = = = = = = = = = = 
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
    
      SUBROUTINE smoothing_fine(iter,ndom,nintf,nnode,nnz,ia,ja,ju,           &
                            au,u,b,nnbd,nbdom,spt,rpt,sintf,rintf)
!     fine 레벨 스무딩 — POL(Chebyshev) 전용 (P3 단일화)
      USE MD_MG_index, ONLY: icheb,rcheb

      IMPLICIT NONE

      INTEGER(4) ndom,iter
      INTEGER(4) nintf,nnode,nnz,nnbd
      INTEGER(4) ia(nnode+1),ja(nnz),ju(nnode)
      INTEGER(4) spt(nnbd+1),rpt(nnbd+1),sintf(spt(nnbd+1)-1),rintf(rpt(nnbd+1)-1),nbdom(nnbd)
      REAL(8) au(nnz),u(nnode),b(nnode)

      INTEGER(4) i

      DO i = 1,iter

         call poly_cheb_smooth(ndom,icheb,rcheb,nintf,nnode,nnz,ia,ja,au,b,nnbd,spt,rpt,nbdom,sintf,rintf,u)

         IF(nnbd.NE.0) THEN
         call send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
         ENDIF

      enddo
!
      RETURN
      END
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
! = = = = = = = = = = = 
    
! = = = = = = = = = = = 
    

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
! = = = = = = = = = = = 
        



