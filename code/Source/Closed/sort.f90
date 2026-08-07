      subroutine sortx_i(a,indx,n)
      implicit none
      integer n
      integer a(n),indx(n)
      integer i,j,k,k1
!
      do i=2,n
         k=a(i)
         if(k.lt.a(i-1)) then
            k1=indx(i)
            a(i)=a(i-1)
            indx(i)=indx(i-1)
            do j=i-2,1,-1
               if(k.lt.a(j)) then
                  a(j+1)=a(j)
                  indx(j+1)=indx(j)
               else
                  goto 100
                endif
            enddo
100         continue
            a(j+1)=k
            indx(j+1)=k1
         endif
      enddo
!
      return
      END SUBROUTINE sortx_i
      subroutine sortx_r(a,indx,n)
      implicit none
      integer n
      integer indx(n)
      real*8  a(n),k
      integer i,j,k1
!
      do i=1,n
         indx(i)=i
      enddo
      do i=2,n
         k=a(i)
         if(a(i-1).gt.k) then
            k1=indx(i)
            a(i)=a(i-1)
            indx(i)=indx(i-1)
            do j=i-2,1,-1
               if(a(j).gt.k) then
                  a(j+1)=a(j)
                  indx(j+1)=indx(j)
               else
                  goto 100
                endif
            enddo
100         continue
            a(j+1)=k
            indx(j+1)=k1
         endif
      enddo
!
      return
      END SUBROUTINE sortx_r
!
      SUBROUTINE sortex_r(a,indx,n)
      implicit none
      integer n
      integer indx(n)
      real*8  a(n),k
      integer i,j,k1
!
      do i=1,n
         indx(i)=i
      enddo
      do i=2,n
         k=a(i)
         if(a(i-1).ge.k) then
            k1=indx(i)
            a(i)=a(i-1)
            indx(i)=indx(i-1)
            do j=i-2,1,-1
               if(a(j).ge.k) then
                  a(j+1)=a(j)
                  indx(j+1)=indx(j)
               else
                  goto 100
                endif
            enddo
100         continue
            a(j+1)=k
            indx(j+1)=k1
         endif
      enddo
!
      END SUBROUTINE sortex_r
      subroutine sortbig(a,indx,n)
      implicit none
      integer blk
      parameter(blk=100)
      integer n
      integer a(n),indx(n)
      integer i,j,k,k1,ip,ip1,ip2
!
      do i=2,blk
         k=a(i)
         if(k.lt.a(i-1)) then
            k1=indx(i)
            a(i)=a(i-1)
            indx(i)=indx(i-1)
            do j=i-2,1,-1
               if(k.lt.a(j)) then
                  a(j+1)=a(j)
                  indx(j+1)=indx(j)
               else
                  goto 100
                endif
            enddo
100         continue
            a(j+1)=k
            indx(j+1)=k1
         endif
      enddo
      do i=blk+1,n
         k=a(i)
         if(k.lt.a(i-1)) then
            k1=indx(i)
            a(i)=a(i-1)
            indx(i)=indx(i-1)
            ip1=1
            ip2=i-2
200         continue
            if(ip1.gt.ip2) goto 210
            ip=(ip1+ip2)/2
            if(k.lt.a(ip)) then
               ip2=ip-1
            elseif(k.gt.a(ip)) then
               ip1=ip+1
            else
               ip2=ip
               goto 210
            endif
            goto 200
210         continue
            ip=ip2+1
            do j=i-1,ip,-1
               a(j+1)=a(j)
               indx(j+1)=indx(j)
            enddo
               a(ip)=k
               indx(ip)=k1
         endif
      enddo
!
      return
      END SUBROUTINE sortbig

