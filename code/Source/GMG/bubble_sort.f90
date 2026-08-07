      subroutine bubble_sort(n,ne)
      implicit none
      integer i,j,n
      integer ne(*)
      integer id
! ---------------------------------------------------------!
            if(n.le.1) return
			
            do i=n,2,-1
               do j=1,i-1
               if(ne(j+1).lt.ne(j)) then
               id=ne(j+1)
               ne(j+1)=ne(j)
               ne(j)=id
               end if
               end do
            end do
!==============check---
!            do i=1,n-1
!               if(ne(i+1).lt.ne(i)) stop" sort- error "
!            end do
! ---
            return
    end
    
! = = = = = = = = = = = = = = 
!  - - - - - - - - - - - - - - - - - - - - 
      subroutine bubble_sort_real(n,ne,dx)
      implicit none
      integer i,j,n
      integer ne(*)
      real*8 dx(*)
      integer id
! ---------------------------------------------------------!
            do i=n,2,-1
               do j=1,i-1
               if(dx(ne(j+1)).lt.dx(ne(j))) then
               id=ne(j+1)
               ne(j+1)=ne(j)
               ne(j)=id
               end if
               end do
            end do
!==============check---
 !           do i=1,n-1
 !              if(dx(ne(i+1)).lt.dx(ne(i))) stop" sort- error "
 !           end do
! ---
            return
    end
    
! = = = = = = = = = = = 
      subroutine bubble_sort_2(n,ne,u)
      implicit none
      integer i,j,n
! out
      REAL(8) u(*)
      integer ne(*)
!temp
      integer id
      REAL(8) xt
      
! ---------------------------------------------------------!
            if(n.le.1) return
			
            do i=n,2,-1
               do j=1,i-1
               if(ne(j+1).lt.ne(j)) then
               id=ne(j+1)
               ne(j+1)=ne(j)
               ne(j)=id
               xt = u(j+1)
               u(j+1) = u(j)
               u(j) = xt
               end if
               end do
            end do
!==============check---
!            do i=1,n-1
!               if(ne(i+1).lt.ne(i)) stop" sort- error "
!            end do
! ---
            return
    end