! = = = = = = = = = = = = = = = = = = = = = = = = = = = 

   subroutine poly_cheb_smooth(np,icheb,rcheb,nintf,n,nnz,ia,ja,au,b,nnbd,si,ri,nbdom,sintf,rintf,x)
   
! smoothing using fourth-kind of Chebyshev 
! refs: 
! Optimal polynomial smoothers for multigrid V-cycles -> method 1: icheb(3) = 1
! later: modify the beta parameter (in paper). 
! OPTIMAL POLYNOMIAL SMOOTHERS FOR PARALLEL AMG -> method-2: 
   
! notes that this time only using M^-1 by Diagonal 
! later using ILU for M^-1
   
!   use md_MPI, only: myrank
   
   implicit none
   
 ! inp:
   integer(4) icheb(5)
   real(8) rcheb(5)
   integer(4) n,nnz,nintf
   integer(4) ia(n+1), ja(nnz)
   real(8) au(nnz),b(n)
  INTEGER(4) nnbd,np
  INTEGER(4) si(nnbd+1),ri(nnbd+1),nbdom(nnbd),sintf(si(nnbd+1)-1),rintf(ri(nnbd+1)-1)
! out
   real(8) x(n)
! temp
   integer(4) :: i, j, k, max_iter, k_max, imethod
!   real(8), parameter :: tol = 1.0e-8_dp
   real(8),dimension(:),allocatable:: r
   real(8),dimension(:),allocatable:: z,y
   real(8) :: eig_min, eig_max, theta, sigma, alpha, beta
   real(8) xtemp,rk,a,ro,ro_new,xtmp,xtmp1
! 
   Allocate(r(n))
   Allocate(z(n),y(n))
  
  ! Estimate the extremal eigenvalues of A
!  call estimate_eigenvalues(a, n, eig_min, eig_max)
  
! set parameter:

  eig_max = rcheb(1)
  
  max_iter = icheb(1)
  k_max = icheb(2)
  imethod = icheb(3)  
!
  

! - - - - - 
   
! this is for method 2: 
   
   a = 0.3
!   ro_new = 0.d0
!   alpha = 0.d0
!   beta = 0.d0
      
! ---
   
      if(nnbd >0) then
      CALL send_receive(nnbd,n,si,ri,sintf,rintf,nbdom,x)
      endif
     call amux0_PCG(nintf,n,nnz,x,r,au,ja,ia) !!A*x=r
     
     r(1:nintf) = (b(1:nintf) - r(1:nintf))/eig_max
     z(1:nintf) = 2.d0/(1.d0+a)*r(1:nintf)
     ro = (1.d0-a)/(1.d0+a)
     
! 
  ! Polynomial smoothing iterations
  do k = 1, max_iter

      x(1:nintf) = x(1:nintf) + z(1:nintf)
      
      if(k .LT. max_iter) then
      
      ro_new = 2.d0*(1.d0+a)/(1.d0-a) - ro
      ro_new = 1.d0/ro_new
      
      if(nnbd >0) then
      CALL send_receive(nnbd,n,si,ri,sintf,rintf,nbdom,z)
      endif
     call amux0_PCG(nintf,n,nnz,z,y,au,ja,ia) !!A*z=y
     
      r(1:nintf) = r(1:nintf) -y(1:nintf)/eig_max
      
      alpha = ro*ro_new
      beta = 4.d0*ro_new/(1.d0-a)
      
      z(1:nintf) = alpha*z(1:nintf) + beta*r(1:nintf)   
      endif
      
  enddo   
   
   

! ---
! - - - - - 
   deallocate(r)
   
   deallocate(z,y)
!
   return 
    END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = 
subroutine eig_value
!
!  POL(Chebyshev) 스무딩의 스펙트럼 상계 — Gershgorin 행합 (G3 확정, P3 단일화)
!  참 상계(λ_max 이상 보장)·분할 무관·통신 1회. 구 Lanczos 추정 경로는
!  분할 의존 요동으로 np-붕괴를 유발해 제거됨 (Finding 문서 참조).
!
USE MD_matrix, only: nnz, ia,ja,au
use MD_MG_index, only: rcheb
use MD_MPI, only: nintf,nprcs

IMPLICIT NONE

real(8) eig_max, xtmp
integer (4) i, j

   eig_max = 0.d0
   do i = 1, nintf
      xtmp = 0.d0
      do j = ia(i), ia(i+1)-1
         xtmp = xtmp + DABS(au(j))
      enddo
      eig_max = MAX(eig_max, xtmp)
   enddo
   if (nprcs > 1) then
      call allreduce_max_r1(eig_max, xtmp)
      eig_max = xtmp
   endif
   rcheb(1) = eig_max
   rcheb(2) = eig_max/30.d0

return
    END
