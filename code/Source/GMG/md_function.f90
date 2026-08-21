module md_function
! 
implicit none
!
      contains
!
! --- irow 
! find the row when a-array poisition are known
      FUNCTION irow(i,init,n,ia)
      integer irow
      integer i,n,j,j1,j2,init
      integer ia(n+1)
! ---
      irow = 0
      Do j = init,n
        j1 = ia(j)
        j2 = ia(j+1)
       if(i.ge.j1.and.i.lt.j2) then
        irow = j
        exit
       end if
      End do
      if(irow.eq.0) then
      write(*,*)'function irow is error'
      stop
      end if
! -------------------
      return
    End function
! --- csr ---
! find the poisition in a-array when the row and the colume are known
      FUNCTION icsr(ir,ic,nr,nnz,ia,ja)
      integer icsr
      integer ir,ic,nr,nnz
      integer ia(nr+1),ja(nnz)	  
      integer i,j,i1,i2
! --- 
      icsr = 0
      i1 = ia(ir)
      i2 = ia(ir+1)-1
      do i = i1,i2
         j = ja(i)
         if(j.eq.ic) then
         icsr = i
         exit
         end if
      end do
      if(icsr.eq.0) then
      write(*,*)'function icsr is error'
      stop
      end if
!  ---------------
      return
      End function
      
!
    end module

! 
! a function definds the norm2 of a vector
    
     FUNCTION XNORM2(n,x)

     implicit none
     integer n,i
     real*8 x(n)
     real*8 XNORM2
     
! 
     XNORM2 = 0.d0
     do i=1,n
         XNORM2 = XNORM2 + x(i)*x(i)
     end do
     XNORM2 = dsqrt(XNORM2)
     
     return
    END
    
