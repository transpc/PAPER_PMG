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
  if(icheb(3) == 2) goto 10 
  
  
! this is for method -1: 
     z = 0.d0
     
      if(nnbd >0) then
      CALL send_receive(nnbd,n,si,ri,sintf,rintf,nbdom,x)
      endif
      call amux0_PCG(nintf,n,nnz,x,r,au,ja,ia) !!A*x=r
     
      r (1:nintf) = b(1:nintf) - r(1:nintf)
! 
  ! Polynomial smoothing iterations
  do k = 1, max_iter
      rk = dble(k)
      
      alpha = (2.d0*rk-3.d0)/(2.d0*rk+1.d0)
      beta = (8.d0*rk-4.d0)/(2.d0*rk+1.d0)
      
      z(1:nintf)  = alpha*z(1:nintf) + beta/eig_max*r(1:nintf) 
      
      x(1:nintf) = x(1:nintf) + z(1:nintf) 
      
      IF(k.LT.max_iter) THEN
          
      if(nnbd >0) then
      CALL send_receive(nnbd,n,si,ri,sintf,rintf,nbdom,z)
      endif
      call amux0_PCG(nintf,n,nnz,z,y,au,ja,ia) !!A*z=y
     
      r(1:nintf) = r(1:nintf) - y(1:nintf)
    
     ENDIF
      
  enddo
	
!
  goto 11
  
! = = = = = = = = = = = = = = = = =  = = = = 
10 continue 
   
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
   
   
11 continue 

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

USE MD_geometry, only: nnode
USE MD_matrix, only: nnz, ia,ja,au
use MD_MG_index, only: icheb, rcheb
use MD_MPI, only: nintf,nprcs
USE MD_MPI_ARP, only: nnbdA,sptA,rptA,sintfA,rintfA,nbdomA
	  
IMPLICIT NONE
      
! temp:
real(8) eig_max, eig_min
integer (4) k

k= icheb(4)

call  lanczos_eig_max(nprcs,k,nintf,nnode,nnz,ia,ja,au,nnbdA,nbdomA,sptA,rptA,sintfA,rintfA,eig_max)

eig_min = eig_max/30.d0

rcheb(1) = eig_max*1.1
rcheb(2) = eig_min

return
    END 
    
! - - - - - - - - - - - - - - - - - - - - - - - 
    
subroutine lanczos_eig_max(np,k,nintf,n,nnz,ia,ja,au,nnbd,nbdom,si,ri,sintf,rintf,eig_max)

! this sub. to find Maximum value of eigenvalue of matrix A(n,n)
! lanczos theory is used; the code is modified from chatGPT

!use md_geometry, only: coord
!use md_MPI, only: myrank

  implicit none
  
! int: 
  integer(4) n,nnz,k,nintf
  integer(4) ia(n+1),ja(nnz)
  real(8) au(nnz)
  INTEGER(4) nnbd,np
  INTEGER(4) si(nnbd+1),ri(nnbd+1),nbdom(nnbd),sintf(si(nnbd+1)-1),rintf(ri(nnbd+1)-1)
  
!out 
  real(8) eig_max
  
 ! temp
!  integer, parameter :: k = 10   ! Number of Lanczos steps
  integer :: i, j, iter
  real(8) :: beta, alpha, tol
  real(8) T(k,k)
  real(8), allocatable :: v(:,:), w(:), v_old(:)    !, y(:)
  real(8) xtmp,xtmp1

  ! Allocate arrays
  allocate(v(n,k))
  allocate(w(n))
  allocate(v_old(n))
!  Allocate(y(n))

  ! Initialize the first vector v(:,1)
!  call random_number(v(:,1))
! notes that in parallel, if set v is difference with that of serial, 
! the results may difference (max_eig)
  
  do i = 1,nintf
      
  v(i,1) = 1.d0    !dsin(xtmp) set to 1 for consisten with parallel
  enddo
  
!
! cal. xtmp = nor2(v,v)
  
      xtmp = 0.d0
      DO i=1,nintf
         xtmp = xtmp + v(i,1)*v(i,1)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreduce_r_s(xtmp,xtmp1)
         xtmp = xtmp1
      ENDIF
    
      xtmp = DSQRT(xtmp)
! ...
!
  v(1:nintf,1) = v(1:nintf,1) / xtmp
  
! - - - - - - - - - - 
  T = 0.d0

  ! Lanczos iteration
  do iter = 1, k
    ! w = A * v(:,iter)
      
      if(nnbd >0) then
      CALL send_receive(nnbd,n,si,ri,sintf,rintf,nbdom,v(1:n,iter))
      endif
     
      call amux0_PCG(nintf,n,nnz,v(:,iter),w,au,ja,ia) !!A*v=w
      
!

    if (iter > 1) then
      w(1:nintf) = w(1:nintf) - beta * v_old(1:nintf)
    end if
!
!
! cal. xtmp = dot_product(v,w):    alpha = dot_product(v(:,iter), w)
  
      xtmp = 0.d0
      DO i=1,nintf
         xtmp = xtmp + v(i,iter)*w(i)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreduce_r_s(xtmp,xtmp1)
         xtmp = xtmp1
      ENDIF
! 
    alpha = xtmp   
    
    w(1:nintf) = w(1:nintf) - alpha * v(1:nintf,iter) 
!
    ! Reorthogonalize
    do j = 1, iter       
!
! ! cal. xtmp = dot_product(v(:,j),w): 
      xtmp = 0.d0
      DO i=1,nintf
         xtmp = xtmp + v(i,j)*w(i)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreduce_r_s(xtmp,xtmp1)
         xtmp = xtmp1
      ENDIF
! - - - - - 
      w(1:nintf) = w(1:nintf) -  xtmp*v(1:nintf,j)
      
    end do

!
!     beta = sqrt(dot_product(w, w))
      xtmp = 0.d0
      DO i=1,nintf
         xtmp = xtmp + w(i)*w(i)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreduce_r_s(xtmp,xtmp1)
         xtmp = xtmp1
      ENDIF

      beta =   dsqrt(xtmp)
!
!    if (beta < tol) exit

    if (iter < k) then
      v_old(1:nintf) = v(1:nintf,iter)
      v(1:nintf,iter+1) = w(1:nintf) / beta
    end if

    T(iter, iter) = alpha
	
    if (iter < k) then
      T(iter, iter+1) = beta
      T(iter+1, iter) = beta
    end if
	
  end do
  
  ! Compute eigenvalues of the tridiagonal matrix T
  call compute_eigenvalues(T, k,eig_max)
  
! test
!  write(*,*)'eig_max = ',eig_max
! 
  deallocate(v,w,v_old)
  
  return
  
  end
  
! = = = = = = = = = = = = = = = = = = = = = = = = 
  subroutine compute_eigenvalues(T, k,eig_max)
  
  use MD_MPI, only: myrank
  implicit none
  
  integer (4) k
    real(8), intent(in) :: T(k,k)

    integer :: info
! out
	real(8) :: eig_max
    
    real(8) d(k),e(k-1),z(k,k),eigvals(k)
  integer (4) i
    
! 
    d= 0.d0
    e = 0.d0
    z = 0.d0
    eigvals = 0.d0
    
    
    do i=1,k
        d(i) = T(i,i)
    enddo
    
    do i=1,k-1
        e(i) = T(i,i+1)
    enddo
    

     call dsteqr('N',k,d,e,z,k, eigvals,  info)
     
!    call dsyev('V', 'U', k, T, k, eigvals, work, lwork, info)
	
 !    if(myrank == 0) then
         
 !   if (info == 0) then
 !     print *, "Maximum eigenvalue: ", maxval(d(1:k))
 !     print *, "Minimum eigenvalue: ", minval(abs(d(1:k)))
 !     print *, "condition: ", maxval(d(1:k))/minval(abs(d(1:k)))
 !   else
 !     print *, "Error in eigenvalue computation, info: ", info
 !   end if
 !    endif
     
	
	eig_max = maxval(d(1:k))    !eigvals(k) -> spectral radius 
    
! test 
!    write(*,*)'test for eigr'
!    eig_max = 1.59367410523645
	
  end subroutine compute_eigenvalues

