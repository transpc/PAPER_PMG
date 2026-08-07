!=========================================================================!
!   ILU PRECONDITIONER FOR MATRIX SOLVER                                  !
!                                                    Modified: Nov. 2016  !    
!=========================================================================!
subroutine pc_ilu1(n,ierr,ia,ja,ju,alu)
implicit none
      
integer :: n,ierr
integer :: ia(*)
integer :: ja(*)
integer :: ju(*)      
real*8  :: alu(*)
    
integer, dimension(:), allocatable :: jwk
integer :: k,j,j1,j2,j3,j4,j5,jrow,jj,jw
real*8  :: tl
!=====
!=====
allocate( jwk(n) )
jwk  = 0
ierr = 0
!   
!/ilu operation 
!
do k=1,n
   j1 = ia(k)
   j2 = ia(k+1)-1
   jwk(ja(j1:j2)) = [j1:j2]
   
   do j=j1,j2
      jrow = ja(j)
      if(jrow>=k) exit
      tl = alu(j)*alu(ju(jrow))
      alu(j) = tl
      do jj = ju(jrow)+1,ia(jrow+1)-1
         jw = jwk(ja(jj))
         if(jw/=0) alu(jw)=alu(jw)-tl*alu(jj)
      enddo
   enddo
   
   if(jrow/=k .or. alu(j)==0.0) then
      write(*,*),'ILUPC: jrow=',jrow,'k=',k, 'alu(j)', alu(j)
      ierr   = 1
      alu(j) = 1.d0
      return
   endif
!...   
   alu(j) = 1.d0/alu(j)
   jwk(ja(j1:j2)) = 0
enddo
!=====
!=====
deallocate(jwk)
return
end subroutine
!=========================================================================!
!   ILU PRECONDITIONER FOR MATRIX SOLVER                                  !
!                                                    Modified: Nov. 2016  !
!=========================================================================!
subroutine pc_ilu2(n,ierr,ia,ja,ju,alu)

implicit none
      
integer :: n,ierr
integer :: ia(*)
integer :: ja(*)
integer :: ju(*)
real*8  :: alu(*)
    
integer, dimension(:), allocatable :: jwk
integer :: k,j,j1,j2,j3,j4,j5,jrow,jj,jw
real*8  :: tl
!=====
!=====
allocate(jwk(n))
jwk  = 0
ierr = 0
!  
!/ilu operation
!
do k=1,n
   j1 = ia(k)
   j2 = ia(k+1)-1
   jwk(ja(j1:j2)) = [j1:j2]
   
   do j=j1,j2 
      jrow = ja(j)
      if(jrow>=k) exit
      tl = alu(j)*alu(ju(jrow))
      alu(j) = tl
      do jj=ju(jrow)+1,ia(jrow+1)-1
         jw = jwk(ja(jj))
         if(jw/=0) alu(jw)=alu(jw)-tl*alu(jj)
      enddo
   enddo
   if(jrow/=k .or. alu(j)==0.d0) then
      ierr   = 1
      alu(j) = 1.d0
   endif
      
   alu(j) = 1.d0/alu(j)
   jwk(ja(j1:j2)) = 0
enddo
!=====
!=====
deallocate(jwk)
return
end subroutine
      
