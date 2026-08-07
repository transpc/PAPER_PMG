
      subroutine Relax_GSs(maxit,relax,n,nnz,ia,ja,ju,au,u,b,diagr)
      
! ---
      implicit none
! ---
      integer maxit,n,nnz
      real*8 relax
      real*8 u(n),b(n),diagr(n),au(nnz)
!
      integer i,j,jj,iter
      integer ia(n+1),ja(nnz),ju(n)
      real*8 temp
      real*8 r(n)
      REAL(8) :: res, XNORM2
      
! ---
!---Gauss-Seidel method
  call Relax_GS0(maxit,relax,n,nnz,ia,ja,au,u,b,diagr)
! Do iter=1,maxit
!  Do i=1,n
!     temp = b(i)
!     do jj=ia(i),ia(i+1)-1
!        j=ja(jj)
!        temp = temp - au(jj)*u(j)
!     enddo
!     u(i)=relax*temp*diagr(i) + u(i)
!  End do
! End do
   
      IF(maxit > 100) THEN
         r = 0.0d0
         CALL mt_amux2(n,n,nnz,u,r,au,ja,ia)
         r = b-r
         res = XNORM2(n,r)
         PRINT *, 'iter : ',iter, 'res :',res
      ENDIF
  
      return
      end

      subroutine Relax_GS0(maxit,relax,n,nnz,ia,ja,au,u,b,diagr)
      
! ---
      implicit none
! ---
      integer maxit,n,nnz
      integer ia(n+1)
      integer ja(nnz)
      real*8 relax
      real*8 au(nnz)    
      real*8 u(n),b(n),diagr(n)
!
      integer i,j,k,iter
      integer jj,nn
      integer j0,j1,j2,j3,j4,j5,j6,j7,j8
      real*8 temp
      
! ---
!---Gauss-Seidel method
  Do iter=1,maxit
   Do i=1,n
     nn=ia(i+1)-ia(i)
     jj=ia(i)
     if(nn.eq.1) then
         j0=ja(jj  )
         temp= b(i)           &
              -au(jj  )*u(j0)
     elseif(nn.eq.2) then
         j0=ja(jj  )
         j1=ja(jj+1)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1)
     elseif(nn.eq.3) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2)
     elseif(nn.eq.4) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3)
     elseif(nn.eq.5) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4)
     elseif(nn.eq.6) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5)
     elseif(nn.eq.7) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6)
     elseif(nn.eq.8) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7)
     elseif(nn.eq.9) then
         j0=ja(jj  )
         j1=ja(jj+1)
         j2=ja(jj+2)
         j3=ja(jj+3)
         j4=ja(jj+4)
         j5=ja(jj+5)
         j6=ja(jj+6)
         j7=ja(jj+7)
         j8=ja(jj+8)
         temp= b(i)           &
              -au(jj  )*u(j0) &
              -au(jj+1)*u(j1) &
              -au(jj+2)*u(j2) &
              -au(jj+3)*u(j3) &
              -au(jj+4)*u(j4) &
              -au(jj+5)*u(j5) &
              -au(jj+6)*u(j6) &
              -au(jj+7)*u(j7) &
              -au(jj+8)*u(j8)
     else
      temp = b(i)
      do jj=ia(i),ia(i+1)-1
         j=ja(jj)
         temp=temp-au(jj)*u(j)
      enddo
     endif
     u(i)=relax*temp*diagr(i) + u(i)
   End do
  End do
!
      return
      end
