!=========================================================================!
!------!-----------!----------------!---------------!------------!--------!
!            CONTAIN SOLOVERS FOR SYSTEM OF ALGEBRAIC EQUATIONS           !
!                                                    Modified: Nov. 2016  !
!=========================================================================!
subroutine mt_amux(n,x,y,a,ja,ia)
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
end subroutine
!=====================================================================!
!=====================================================================!
subroutine mt_lusol(n,y,x,alu,ja,ia,ju)
implicit none
!
!---------------------------------------------------------------------!
!                                                                     !
!  Performing a forward followed by a backward solve                  !
!  for LU matrix as produced by  ILUT                                 !
!  input:                                                             !
!    n     = row dimension of A                                       !
!    y     = the right hand side of the linear system                 !
!    alu, ja, ia = input matrix in compressed sparse row format.      !
!  output:                                                            !
!    x     = the solution                                             !
!                                                                     !
!---------------------------------------------------------------------!
!
integer :: n,i,k1,k2
real*8  :: tmp
integer :: ja(*)
integer :: ju(*)
integer :: ia(*)
real*8  :: alu(*)
real*8, dimension(n) :: x,y      
!
!--forward solver
!
x(1) = y(1)
do i=2,n
   k1 = ia(i)
   k2 = ju(i)-1
   tmp = y(i) - dot_product( alu(k1:k2) , x(ja(k1:k2)) )
   x(i) = tmp
enddo
!
!--backward solve
!
do i=n,1,-1
   k1 = ju(i)+1
   k2 = ia(i+1)-1
   tmp = x(i) - dot_product( alu(k1:k2) , x(ja(k1:k2)) )
   x(i) = tmp*alu(ju(i))
enddo
!=====
!=====
return
end subroutine
!=====================================================================!
!=====================================================================!
subroutine solve_cg(n,ierr,err_crt,maxit,ia,ja,ju,au,alu,rhs,u)
implicit none
!
!---------------------------------------------------------------------!
!   Preconditioned conjugate gradient solver by ILUPC                 !
!   [A]{x} = {b}                                                      !
!       1. r := b-[A]{x}                                              !
!          z := [M]^-1*{r}                                            !
!          p := z                                                     !
!       For k=0,maxit                                                 !
!           2.  alpha  = (r,z)/(A*p,p)                                !
!           3.  x      = x + alpha*p                                  !
!           4.  r      = r - alpha*([A]{p})                           !
!           5.  z      = [M]^-1{r}                                    !
!           6.  beta   = (r,z)/(ro*zo)                                !
!           7.  p      = z + beta*p                                   !
!       End for                                                       !
!---------------------------------------------------------------------!
!
integer :: n,ierr,iter,maxit
integer :: ia(*)
integer :: ju(*)
integer :: ja(*)
real*8 :: au(*)
real*8 :: alu(*)
real*8 :: rhs(*)
real*8 :: u(n)
real*8 :: err_crt,ak0,bk0,ak,bk,err0,err,small,XNORM2
real*8, dimension(:), allocatable :: r,z,p,sol
!=====
!=====
allocate( r(n), z(n), p(n), sol(n))
small=1.d-20
iter=0
!  ...
!  Predictor  By Diag.
!  ...
!sol = rhs(1:n)/au(ju(1:n))
sol = u
call mt_amux(n,sol,r,au,ja,ia)
!  ...
!  -1-

r(1:n)    = rhs(1:n)-r(1:n)
!err0 = norm2(r)
err0 = XNORM2(n,r)

!write(*,*)'ro=',err0
!pause
if(err0.lt.small) goto 100
!write(100,*) 0, log10(1.d0)
!  ...
!  
z = 0.d0
call mt_lusol(n,r,z,alu,ja,ia,ju)
!  ...
p = z
!  ...
!  ...
ak0 = dot_product(r,z)
bk0 = ak0
do iter=1,maxit
!  ...
!  -2-
   call mt_amux(n,p,z,au,ja,ia)
   ak = ak0/dot_product(p,z)
!  ...
!  -3-
   sol = sol+ak*p
!  ...
!  -4-
   r = r-ak*z
!  ...
!  -5-
   call mt_lusol(n,r,z,alu,ja,ia,ju)
!  ...
!  -6-
   ak0 = dot_product(r,z)
   bk  = ak0/bk0
   bk0 = ak0
!  ...
!  -7-
   p   = z+bk*p
!  ...
!  Check convergence
!   err = norm2(r)
   err = XNORM2(n,r)
   
   if(mod(iter,100)==0) then
      write(*,124),iter,err/(err0)
   endif
   if(dabs(err)>1.e20) then
      write(*,125)
      stop
   endif
   if ( err/(err0)<err_crt ) exit
!
!    write(100,*)iter,log10(err/err0)
!
enddo
100   continue
u(1:n) = sol
!write(*,124),iter,err/(err0)
deallocate(r,z,p,sol)
!=====
!=====
123 format('+---------+------------+----------+------------+')
124 format('!   + CG  ::   its = ',i7,'  |   Error = ',E15.6 ,'   !')
125 format('! CG solver :: Blow-up                         !')
return  
end subroutine
!=====================================================================!
!=====================================================================!

!=====================================================================!
subroutine solve_cg2(n,nnz,err_crt,maxit,ia,ja,ju,au,alu,rhs,u)
implicit none
!
!---------------------------------------------------------------------!
!   Preconditioned conjugate gradient solver by ILUPC                 !
!   [A]{x} = {b}                                                      !
!       1. r := b-[A]{x}                                              !
!          z := [M]^-1*{r}                                            !
!          p := z                                                     !
!       For k=0,maxit                                                 !
!           2.  alpha  = (r,z)/(A*p,p)                                !
!           3.  x      = x + alpha*p                                  !
!           4.  r      = r - alpha*([A]{p})                           !
!           5.  z      = [M]^-1{r}                                    !
!           6.  beta   = (r,z)/(ro*zo)                                !
!           7.  p      = z + beta*p                                   !
!       End for                                                       !
!---------------------------------------------------------------------!
!
integer :: n,ierr,iter,maxit,nnz
integer :: ia(n+1)
integer :: ju(n)
integer :: ja(nnz)
real*8 :: au(nnz)
real*8 :: alu(nnz)
real*8 :: rhs(n)
real*8 :: u(n)
real*8 :: err_crt,ak0,bk0,ak,bk,err0,err,small,XNORM2
real*8, dimension(n) :: r,z,p,sol
!=====
!=====
!allocate( r(n), z(n), p(n), sol(n))
small=1.d-20
iter=0
!  ...
!  Predictor  By Diag.
!  ...
!sol = rhs(1:n)/au(ju(1:n))
sol = u
call mt_amux(n,sol,r,au,ja,ia)
!  ...
!  -1-

r(1:n)    = rhs(1:n)-r(1:n)
!err0 = norm2(r)
err0 = XNORM2(n,r)

!write(*,*)'ro=',err0
!pause
if(err0.lt.small) goto 100
!write(100,*) 0, log10(1.d0)
!  ...
!
z = 0.d0
call mt_lusol(n,r,z,alu,ja,ia,ju)
!  ...
p = z
!  ...
!  ...
ak0 = dot_product(r,z)
bk0 = ak0
do iter=1,maxit
!  ...
!  -2-
   call mt_amux(n,p,z,au,ja,ia)
   ak = ak0/dot_product(p,z)
!  ...
!  -3-
   sol = sol+ak*p
!  ...
!  -4-
   r = r-ak*z
!  ...
!  -5-
   call mt_lusol(n,r,z,alu,ja,ia,ju)
!  ...
!  -6-
   ak0 = dot_product(r,z)
   bk  = ak0/bk0
   bk0 = ak0
!  ...
!  -7-
   p   = z+bk*p
!  ...
!  Check convergence
!   err = norm2(r)
   err = XNORM2(n,r)

!   if(mod(iter,100)==0) then
!      write(*,124),iter,err/(err0)
!   endif
   if(dabs(err)>1.e10) then
      write(999,125)
      stop
   endif
   if ( err/(err0)<err_crt ) exit
!
!    write(100,*)iter,log10(err/err0)
!
enddo
100   continue
u(1:n) = sol
!write(*,*),iter,err/(err0)
!deallocate(r,z,p,sol)
!=====
!=====
123 format('+---------+------------+----------+------------+')
124 format('!   + CG  ::   its = ',i7,'  |   Error = ',E15.6 ,'   !')
125 format('! CG solver :: Blow-up                         !')
return
end subroutine
!=====================================================================!
!=====================================================================!

      
