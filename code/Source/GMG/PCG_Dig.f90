!**************************************************      - - - - !
! parallel CG using Diagonal Preconditioning - - - - - - - - - - !
! Modify by S.T Ha - - - - - - - - - - - - - - - - - - - - - - - !
! June 2021 - - - - - - - - - - - - - - - - - - - - - - - - - - -!

!***********************************************************************
      SUBROUTINE amux0_PCG(nintf,n,nnz,x,y,a,ja,ia)
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
! - - - - - - - - - 
