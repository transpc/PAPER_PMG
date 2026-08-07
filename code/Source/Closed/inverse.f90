      SUBROUTINE luinverse3m(a,ip,b,n,m)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: n,m 
!.....Output
      INTEGER :: ip(n,2)
      REAL(8) :: a(n,3,3)
      REAL(8) :: b(n,3,m)
!.....Local variables
      INTEGER :: i,l
      INTEGER :: ip1,ip2
      REAL(8) :: pv,t
      REAL(8) :: a11,a12,a13
      REAL(8) :: a21,a22,a23
      REAL(8) :: a31,a32,a33
      REAL(8) :: t1,t2,t3
!
      DO i=1,n
         ip1=1
         pv=a(i,1,1)
         t2=a(i,2,1)
         t3=a(i,3,1)
         IF(abs(pv).lt.abs(t2)) then
            pv=t2
            ip1=2
         ENDIF
         IF(abs(pv).lt.abs(t3)) then
            pv=t3
            ip1=3
         ENDIF
         ip(i,1)=ip1
         pv=1.d0/pv
         IF(ip1.ne.1) then
            a12=a(i,ip1,2)
            a13=a(i,ip1,3)
            a(i,ip1,1)=a(i,1,1)
            a(i,ip1,2)=a(i,1,2)
            a(i,ip1,3)=a(i,1,3)
            a(i,1,2)=a12
            a(i,1,3)=a13
            DO l=1,m
               t=b(i,ip1,l)
               b(i,ip1,l)=b(i,1,l)
               b(i,1,l)=t
            ENDDO
         else
            a12=a(i,1,2)
            a13=a(i,1,3)
         ENDIF
         a11=pv
         a(i,1,1)=a11
         a(i,2,1)=a(i,2,1)*pv
         a(i,3,1)=a(i,3,1)*pv
!
         a(i,2,2)=a(i,2,2)-a(i,2,1)*a12
         a(i,3,2)=a(i,3,2)-a(i,3,1)*a12
!
         ip2=2
         pv=a(i,2,2)
         t3=a(i,3,2)
         IF(abs(pv).lt.abs(t3)) then
            pv=t3
            ip2=3
         ENDIF
         ip(i,2)=ip2
         pv=1.d0/pv
         IF(ip2.ne.2) then
            a21=a(i,ip2,1)
            a23=a(i,ip2,3)
            a(i,ip2,1)=a(i,2,1)
            a(i,ip2,2)=a(i,2,2)
            a(i,ip2,3)=a(i,2,3)
            a(i,2,1)=a21
            a(i,2,3)=a23
            DO l=1,m
               t=b(i,ip2,l)
               b(i,ip2,l)=b(i,2,l)
               b(i,2,l)=t
            ENDDO
         ELSE
            a21=a(i,2,1)
            a23=a(i,2,3)
         ENDIF
         a22=pv
         a(i,2,2)=a22
         a32=a(i,3,2)*pv
         a31=a(i,3,1)
         a23=a23-a21*a13
         a(i,3,2)=a32
         a(i,2,3)=a23
!
         a33=a(i,3,3)-a31*a13-a32*a23
         a33=1.d0/a33
         a(i,3,3)=a33
!
      ENDDO
      DO l=1,m
         DO i=1,n
            t1=b(i,1,l)
            t2=b(i,2,l)
            t3=b(i,3,l)
            t2=t2-a(i,2,1)*t1
            t3=t3-a(i,3,1)*t1-a(i,3,2)*t2
!
            b(i,3,l)=t3*a(i,3,3)
            t2=t2-a(i,2,3)*b(i,3,l)
            b(i,2,l)=t2*a(i,2,2)
            t1=t1-a(i,1,2)*b(i,2,l)-a(i,1,3)*b(i,3,l)
            b(i,1,l)=t1*a(i,1,1)
         ENDDO
      ENDDO
!
      END SUBROUTINE luinverse3m
!
!----------------------------------------------------------------------
!
      SUBROUTINE solve3m(a,ip,b,n,m)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: n,m
      INTEGER :: ip(n,2)
      REAL(8) :: a(n,3,3)
!.....Output
      REAL(8) :: b(n,3,m)
!.....Local variables
      INTEGER :: i,l
      INTEGER :: ip1,ip2
      REAL(8) :: t1,t2,t3
!
      DO l=1,m
!DIR$ NOVECTOR
         DO i=1,n
            ip1=ip(i,1)
            ip2=ip(i,2)
            t1=b(i,ip1,l)
            b(i,ip1,l)=b(i,1,l)
            t2=b(i,ip2,l)
            b(i,ip2,l)=b(i,2,l)
            t3=b(i,3,l)
            t2=t2-a(i,2,1)*t1
            t3=t3-a(i,3,1)*t1-a(i,3,2)*t2
!
            b(i,3,l)=t3*a(i,3,3)
            t2=t2-a(i,2,3)*b(i,3,l)
            b(i,2,l)=t2*a(i,2,2)
            t1=t1-a(i,1,2)*b(i,2,l)-a(i,1,3)*b(i,3,l)
            b(i,1,l)=t1*a(i,1,1)
         ENDDO
      ENDDO
!
      END SUBROUTINE solve3m
!----------------------------------------------------------------------
      SUBROUTINE luinverse31(a,ip,b,n,npb)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: n
      INTEGER :: npb(n)
!.....Output
      INTEGER :: ip(n,2)
      REAL(8) :: a(n,3,3)
      REAL(8) :: b(n,3)
!.....Local variables
      INTEGER :: i
      INTEGER :: ip1,ip2
      REAL(8) :: pv,t
      REAL(8) :: a11,a12,a13
      REAL(8) :: a21,a22,a23
      REAL(8) :: a31,a32,a33
      REAL(8) :: t1,t2,t3
!
      DO i=1,n
         IF(npb(i).ne.0) CYCLE
         ip1=1
         pv=a(i,1,1)
         t2=a(i,2,1)
         t3=a(i,3,1)
         IF(abs(pv).lt.abs(t2)) then
            pv=t2
            ip1=2
         ENDIF
         IF(abs(pv).lt.abs(t3)) then
            pv=t3
            ip1=3
         ENDIF
         ip(i,1)=ip1
         pv=1.d0/pv
         IF(ip1.ne.1) then
            a12=a(i,ip1,2)
            a13=a(i,ip1,3)
            a(i,ip1,1)=a(i,1,1)
            a(i,ip1,2)=a(i,1,2)
            a(i,ip1,3)=a(i,1,3)
            a(i,1,2)=a12
            a(i,1,3)=a13
            t=b(i,ip1)
            b(i,ip1)=b(i,1)
            b(i,1)=t
         else
            a12=a(i,1,2)
            a13=a(i,1,3)
         ENDIF
         a11=pv
         a(i,1,1)=a11
         a(i,2,1)=a(i,2,1)*pv
         a(i,3,1)=a(i,3,1)*pv
!
         a(i,2,2)=a(i,2,2)-a(i,2,1)*a12
         a(i,3,2)=a(i,3,2)-a(i,3,1)*a12
!
         ip2=2
         pv=a(i,2,2)
         t3=a(i,3,2)
         IF(abs(pv).lt.abs(t3)) then
            pv=t3
            ip2=3
         ENDIF
         ip(i,2)=ip2
         pv=1.d0/pv
         IF(ip2.ne.2) then
            a21=a(i,ip2,1)
            a23=a(i,ip2,3)
            a(i,ip2,1)=a(i,2,1)
            a(i,ip2,2)=a(i,2,2)
            a(i,ip2,3)=a(i,2,3)
            a(i,2,1)=a21
            a(i,2,3)=a23
            t=b(i,ip2)
            b(i,ip2)=b(i,2)
            b(i,2)=t
         ELSE
            a21=a(i,2,1)
            a23=a(i,2,3)
         ENDIF
         a22=pv
         a(i,2,2)=a22
         a32=a(i,3,2)*pv
         a31=a(i,3,1)
         a23=a23-a21*a13
         a(i,3,2)=a32
         a(i,2,3)=a23
!
         a33=a(i,3,3)-a31*a13-a32*a23
         a33=1.d0/a33
         a(i,3,3)=a33
      ENDDO
      DO i=1,n
         IF(npb(i).ne.0) CYCLE
         t1=b(i,1)
         t2=b(i,2)
         t3=b(i,3)
         t2=t2-a(i,2,1)*t1
         t3=t3-a(i,3,1)*t1-a(i,3,2)*t2
!
         b(i,3)=t3*a(i,3,3)
         t2=t2-a(i,2,3)*b(i,3)
         b(i,2)=t2*a(i,2,2)
         t1=t1-a(i,1,2)*b(i,2)-a(i,1,3)*b(i,3)
         b(i,1)=t1*a(i,1,1)
      ENDDO
!
      END SUBROUTINE luinverse31
!
!----------------------------------------------------------------------
!
      SUBROUTINE solve31(a,ip,b,n,npb)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: n
      INTEGER :: npb(n)
      INTEGER :: ip(n,2)
      REAL(8) :: a(n,3,3)
!.....Output
      REAL(8) :: b(n,3)
!.....Local variables
      INTEGER :: i
      INTEGER :: ip1,ip2
      REAL(8) :: t1,t2,t3
!
!DIR$ NOVECTOR
      DO i=1,n
         IF(npb(i).ne.0) CYCLE
         ip1=ip(i,1)
         ip2=ip(i,2)
         t1=b(i,ip1)
         b(i,ip1)=b(i,1)
         t2=b(i,ip2)
         b(i,ip2)=b(i,2)
         t3=b(i,3)
         t2=t2-a(i,2,1)*t1
         t3=t3-a(i,3,1)*t1-a(i,3,2)*t2
!
         b(i,3)=t3*a(i,3,3)
         t2=t2-a(i,2,3)*b(i,3)
         b(i,2)=t2*a(i,2,2)
         t1=t1-a(i,1,2)*b(i,2)-a(i,1,3)*b(i,3)
         b(i,1)=t1*a(i,1,1)
      ENDDO
!
      END SUBROUTINE solve31
!      
!----------------------------------------------------------------------
!
      SUBROUTINE luinverse6m(a,b,ldb,m)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: m,ldb
!.....Output
      REAL(8) :: a(6,6)
      REAL(8) :: b(ldb,6)
!.....Local variables
      INTEGER :: l
      REAL(8) :: pv
      REAL(8) :: t1,t2,t3,t4,t5,t6
      INTEGER :: ip1,ip2,ip3,ip4,ip5
!
      pv=a(1,1)
      ip1=1
      IF(abs(pv).lt.abs(a(2,1))) then
         ip1=2
         pv=a(2,1)
      ENDIF
      IF(abs(pv).lt.abs(a(3,1))) then
         ip1=3
         pv=a(3,1)
      ENDIF
      IF(abs(pv).lt.abs(a(4,1))) then
         ip1=4
         pv=a(4,1)
      ENDIF
      IF(abs(pv).lt.abs(a(5,1))) then
         ip1=5
         pv=a(5,1)
      ENDIF
      IF(abs(pv).lt.abs(a(6,1))) then
         ip1=6
         pv=a(6,1)
      ENDIF
      IF(ip1.ne.1) then
         t2=a(ip1,2)
         t3=a(ip1,3)
         t4=a(ip1,4)
         t5=a(ip1,5)
         t6=a(ip1,6)
         a(ip1,1)=a(1,1)
         a(ip1,2)=a(1,2)
         a(ip1,3)=a(1,3)
         a(ip1,4)=a(1,4)
         a(ip1,5)=a(1,5)
         a(ip1,6)=a(1,6)
         a(1,2)=t2
         a(1,3)=t3
         a(1,4)=t4
         a(1,5)=t5
         a(1,6)=t6
      ENDIF
      pv=1.d0/pv
      a(1,1)=pv
      a(2,1)=a(2,1)*pv
      a(3,1)=a(3,1)*pv
      a(4,1)=a(4,1)*pv
      a(5,1)=a(5,1)*pv
      a(6,1)=a(6,1)*pv
!-----
      t1=a(1,2)
      a(2,2)=a(2,2)-a(2,1)*t1
      a(3,2)=a(3,2)-a(3,1)*t1
      a(4,2)=a(4,2)-a(4,1)*t1
      a(5,2)=a(5,2)-a(5,1)*t1
      a(6,2)=a(6,2)-a(6,1)*t1

      ip2=2
      pv=a(2,2)
      IF(abs(pv).lt.abs(a(3,2))) then
         ip2=3
         pv=a(3,2)
      ENDIF
      IF(abs(pv).lt.abs(a(4,2))) then
         ip2=4
         pv=a(4,2)
      ENDIF
      IF(abs(pv).lt.abs(a(5,2))) then
         ip2=5
         pv=a(5,2)
      ENDIF
      IF(abs(pv).lt.abs(a(6,2))) then
         ip2=6
         pv=a(6,2)
      ENDIF
!
      IF(ip2.ne.2) then
         t1=a(ip2,1)
         t3=a(ip2,3)
         t4=a(ip2,4)
         t5=a(ip2,5)
         t6=a(ip2,6)
         a(ip2,1)=a(2,1)
         a(ip2,2)=a(2,2)
         a(ip2,3)=a(2,3)
         a(ip2,4)=a(2,4)
         a(ip2,5)=a(2,5)
         a(ip2,6)=a(2,6)
         a(2,1)=t1
         a(2,3)=t3
         a(2,4)=t4
         a(2,5)=t5
         a(2,6)=t6
      ENDIF
!
      pv=1.d0/pv
      a(2,2)=pv
      a(3,2)=a(3,2)*pv
      a(4,2)=a(4,2)*pv
      a(5,2)=a(5,2)*pv
      a(6,2)=a(6,2)*pv
      t1=a(2,1)
      a(2,3)=a(2,3)-t1*a(1,3)
      a(2,4)=a(2,4)-t1*a(1,4)
      a(2,5)=a(2,5)-t1*a(1,5)
      a(2,6)=a(2,6)-t1*a(1,6)
      t1=a(1,3)
      t2=a(2,3)
      a(3,3)=a(3,3)-a(3,1)*t1-a(3,2)*t2
      a(4,3)=a(4,3)-a(4,1)*t1-a(4,2)*t2
      a(5,3)=a(5,3)-a(5,1)*t1-a(5,2)*t2
      a(6,3)=a(6,3)-a(6,1)*t1-a(6,2)*t2
!
      ip3=3
      pv=a(3,3)
      IF(abs(pv).lt.abs(a(4,3))) then
         ip3=4
         pv=a(4,3)
      ENDIF
      IF(abs(pv).lt.abs(a(5,3))) then
         ip3=5
         pv=a(5,3)
      ENDIF
      IF(abs(pv).lt.abs(a(6,3))) then
         ip3=6
         pv=a(6,3)
      ENDIF
!
      IF(ip3.ne.3) then
         t1=a(ip3,1)
         t2=a(ip3,2)
         t4=a(ip3,4)
         t5=a(ip3,5)
         t6=a(ip3,6)
         a(ip3,1)=a(3,1)
         a(ip3,2)=a(3,2)
         a(ip3,4)=a(3,4)
         a(ip3,3)=a(3,3)
         a(ip3,5)=a(3,5)
         a(ip3,6)=a(3,6)
         a(3,1)=t1
         a(3,2)=t2
         a(3,4)=t4
         a(3,5)=t5
         a(3,6)=t6
      ENDIF
!
      pv=1.d0/pv
      a(3,3)=pv
      a(4,3)=a(4,3)*pv
      a(5,3)=a(5,3)*pv
      a(6,3)=a(6,3)*pv
      t1=a(3,1)
      t2=a(3,2)
      a(3,4)=a(3,4)-t1*a(1,4)-t2*a(2,4)
      a(3,5)=a(3,5)-t1*a(1,5)-t2*a(2,5)
      a(3,6)=a(3,6)-t1*a(1,6)-t2*a(2,6)
      t1=a(1,4)
      t2=a(2,4)
      t3=a(3,4)
      a(4,4)=a(4,4)-a(4,1)*t1-a(4,2)*t2-a(4,3)*t3
      a(5,4)=a(5,4)-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3
      a(6,4)=a(6,4)-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3
!
      ip4=4
      pv=a(4,4)
      IF(abs(pv).lt.abs(a(5,4))) then
         ip4=5
         pv=a(5,4)
      ENDIF
      IF(abs(pv).lt.abs(a(6,4))) then
         ip4=6
         pv=a(6,4)
      ENDIF
!-----
      IF(ip4.ne.4) then
         t1=a(ip4,1)
         t2=a(ip4,2)
         t3=a(ip4,3)
         t5=a(ip4,5)
         t6=a(ip4,6)
         a(ip4,1)=a(4,1)
         a(ip4,2)=a(4,2)
         a(ip4,3)=a(4,3)
         a(ip4,4)=a(4,4)
         a(ip4,5)=a(4,5)
         a(ip4,6)=a(4,6)
         a(4,1)=t1
         a(4,2)=t2
         a(4,3)=t3
         a(4,5)=t5
         a(4,6)=t6
      ENDIF
!-----
      pv=1.d0/pv
      a(4,4)=pv
      a(5,4)=a(5,4)*pv
      a(6,4)=a(6,4)*pv
      t1=a(4,1)
      t2=a(4,2)
      t3=a(4,3)
      a(4,5)=a(4,5)-t1*a(1,5)-t2*a(2,5)-t3*a(3,5)
      a(4,6)=a(4,6)-t1*a(1,6)-t2*a(2,6)-t3*a(3,6)
      t1=a(1,5)
      t2=a(2,5)
      t3=a(3,5)
      t4=a(4,5)
      a(5,5)=a(5,5)-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
      a(6,5)=a(6,5)-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4
!
      ip5=5
      pv=a(5,5)
      IF(abs(pv).lt.abs(a(6,5))) then
         pv=a(6,5)
         ip5=6
      ENDIF
!
      IF(ip5.ne.5) then
         t1=a(ip5,1)
         t2=a(ip5,2)
         t3=a(ip5,3)
         t4=a(ip5,4)
         t6=a(ip5,6)
         a(ip5,1)=a(5,1)
         a(ip5,2)=a(5,2)
         a(ip5,3)=a(5,3)
         a(ip5,4)=a(5,4)
         a(ip5,5)=a(5,5)
         a(ip5,6)=a(5,6)
         a(5,1)=t1
         a(5,2)=t2
         a(5,3)=t3
         a(5,4)=t4
         a(5,6)=t6
      ENDIF
!
!-----
      pv=1.d0/pv
      a(5,5)=pv
      a(6,5)=a(6,5)*pv
      t1=a(1,6)
      t2=a(2,6)
      t3=a(3,6)
      t4=a(4,6)
      t5=a(5,6)
      t6=a(6,6)
      t5=t5-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
      a(5,6)=t5
      t6=t6-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4-a(6,5)*t5
      t6=1.d0/t6
      a(6,6)=t6
!-----
!DIR$ NOVECTOR
      DO l=1,m
         t1=b(l,ip1)
         b(l,ip1)=b(l,1)
         t2=b(l,ip2)
         b(l,ip2)=b(l,2)
         t3=b(l,ip3)
         b(l,ip3)=b(l,3)
         t4=b(l,ip4)
         b(l,ip4)=b(l,4)
         t5=b(l,ip5)
         b(l,ip5)=b(l,5)
         t6=b(l,6)
         t2=t2-a(2,1)*t1
         t3=t3-a(3,1)*t1-a(3,2)*t2
         t4=t4-a(4,1)*t1-a(4,2)*t2-a(4,3)*t3
         t5=t5-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
         t6=t6-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4-a(6,5)*t5
!
         b(l,6)=t6*a(6,6)
         t5=t5-a(5,6)*b(l,6)
         b(l,5)=t5*a(5,5)
         t4=t4-a(4,5)*b(l,5)-a(4,6)*b(l,6)
         b(l,4)=t4*a(4,4)
         t3=t3-a(3,4)*b(l,4)-a(3,5)*b(l,5)-a(3,6)*b(l,6)
         b(l,3)=t3*a(3,3)
         t2=t2-a(2,3)*b(l,3)-a(2,4)*b(l,4)-a(2,5)*b(l,5)-a(2,6)*b(l,6)
         b(l,2)=t2*a(2,2)
         t1=t1-a(1,2)*b(l,2)-a(1,3)*b(l,3)-a(1,4)*b(l,4)-a(1,5)*b(l,5)-a(1,6)*b(l,6)
         b(l,1)=t1*a(1,1)
      ENDDO
!
      END SUBROUTINE luinverse6m 
!
!----------------------------------------------------------------------
!
      SUBROUTINE luinverse6mn(aa,b,nblk,ldb,n,npb,nrhs)
!
      IMPLICIT NONE
!     
!.....Input
      INTEGER :: nblk,n,m,ldb
      INTEGER :: npb(*),nrhs(*)
      REAL(8) :: aa(nblk,6,6)
      REAL(8) :: a(6,6)
      REAL(8) :: b(ldb,6,nblk)
!.....Local variables
      INTEGER :: i,l,ix
      REAL(8) :: pv
      REAL(8) :: a12
      REAL(8) :: a13,a23
      REAL(8) :: a21
      REAL(8) :: t
      REAL(8) :: t1,t2,t3,t4,t5,t6
      INTEGER :: ip1,ip2,ip3,ip4,ip5
!
      DO i=1,n
         IF(npb(i).ne.0) CYCLE
         m=1+3*nrhs(i)
!DIR$ NOVECTOR
         DO ix=1,6
            a(ix,1)=aa(i,ix,1)
            a(ix,2)=aa(i,ix,2)
            a(ix,3)=aa(i,ix,3)
            a(ix,4)=aa(i,ix,4)
            a(ix,5)=aa(i,ix,5)
            a(ix,6)=aa(i,ix,6)
         ENDDO
!
         pv=a(1,1)
         ip1=1
         IF(abs(pv).lt.abs(a(2,1))) then
            ip1=2
            pv=a(2,1)
         ENDIF
         IF(abs(pv).lt.abs(a(3,1))) then
            ip1=3
            pv=a(3,1)
         ENDIF
         IF(abs(pv).lt.abs(a(4,1))) then
            ip1=4
            pv=a(4,1)
         ENDIF
         IF(abs(pv).lt.abs(a(5,1))) then
            ip1=5
            pv=a(5,1)
         ENDIF
         IF(abs(pv).lt.abs(a(6,1))) then
            ip1=6
            pv=a(6,1)
         ENDIF
         IF(ip1.ne.1) then
            t2=a(ip1,2)
            t3=a(ip1,3)
            t4=a(ip1,4)
            t5=a(ip1,5)
            t6=a(ip1,6)
            a(ip1,1)=a(1,1)
            a(ip1,2)=a(1,2)
            a(ip1,3)=a(1,3)
            a(ip1,4)=a(1,4)
            a(ip1,5)=a(1,5)
            a(ip1,6)=a(1,6)
            a(1,2)=t2
            a(1,3)=t3
            a(1,4)=t4
            a(1,5)=t5
            a(1,6)=t6
            DO l=1,m
               t=b(l,ip1,i)
               b(l,ip1,i)=b(l,1,i)
               b(l,1,i)=t
            ENDDO
         ENDIF
         pv=1.d0/pv
         a(1,1)=pv
         a12=a(1,2)
         DO ix=2,6
            a(ix,1)=a(ix,1)*pv
            a(ix,2)=a(ix,2)-a(ix,1)*a12
         ENDDO
!
         ip2=2
         pv=a(2,2)
         IF(abs(pv).lt.abs(a(3,2))) then
            ip2=3
            pv=a(3,2)
         ENDIF
         IF(abs(pv).lt.abs(a(4,2))) then
            ip2=4
            pv=a(4,2)
         ENDIF
         IF(abs(pv).lt.abs(a(5,2))) then
            ip2=5
            pv=a(5,2)
         ENDIF
         IF(abs(pv).lt.abs(a(6,2))) then
            ip2=6
            pv=a(6,2)
         ENDIF
!
         IF(ip2.ne.2) then
            t1=a(ip2,1)
            t3=a(ip2,3)
            t4=a(ip2,4)
            t5=a(ip2,5)
            t6=a(ip2,6)
            a(ip2,1)=a(2,1)
            a(ip2,2)=a(2,2)
            a(ip2,3)=a(2,3)
            a(ip2,4)=a(2,4)
            a(ip2,5)=a(2,5)
            a(ip2,6)=a(2,6)
            a(2,1)=t1
            a(2,3)=t3
            a(2,4)=t4
            a(2,5)=t5
            a(2,6)=t6
            DO l=1,m
               t=b(l,ip2,i)
               b(l,ip2,i)=b(l,2,i)
               b(l,2,i)=t
            ENDDO
         ENDIF
!
         pv=1.d0/pv
         a(2,2)=pv
         a21=a(2,1)
         DO ix=3,6
            a(ix,2)=a(ix,2)*pv
            a(2,ix)=a(2,ix)-a21*a(1,ix)
         ENDDO
         a13=a(1,3)
         a23=a(2,3)
         DO ix=3,6
            a(ix,3)=a(ix,3)-a(ix,1)*a13-a(ix,2)*a23
         ENDDO
!
         ip3=3
         pv=a(3,3)
         IF(abs(pv).lt.abs(a(4,3))) then
            ip3=4
            pv=a(4,3)
         ENDIF
         IF(abs(pv).lt.abs(a(5,3))) then
            ip3=5
            pv=a(5,3)
         ENDIF
         IF(abs(pv).lt.abs(a(6,3))) then
            ip3=6
            pv=a(6,3)
         ENDIF
!
         IF(ip3.ne.3) then
            t1=a(ip3,1)
            t2=a(ip3,2)
            t4=a(ip3,4)
            t5=a(ip3,5)
            t6=a(ip3,6)
            a(ip3,1)=a(3,1)
            a(ip3,2)=a(3,2)
            a(ip3,4)=a(3,4)
            a(ip3,3)=a(3,3)
            a(ip3,5)=a(3,5)
            a(ip3,6)=a(3,6)
            a(3,1)=t1
            a(3,2)=t2
            a(3,4)=t4
            a(3,5)=t5
            a(3,6)=t6
            DO l=1,m
               t=b(l,ip3,i)
               b(l,ip3,i)=b(l,3,i)
               b(l,3,i)=t
            ENDDO
         ENDIF
!
         pv=1.d0/pv
         a(3,3)=pv
         a(4,3)=a(4,3)*pv
         a(5,3)=a(5,3)*pv
         a(6,3)=a(6,3)*pv
         t1=a(3,1)
         t2=a(3,2)
         a(3,4)=a(3,4)-t1*a(1,4)-t2*a(2,4)
         a(3,5)=a(3,5)-t1*a(1,5)-t2*a(2,5)
         a(3,6)=a(3,6)-t1*a(1,6)-t2*a(2,6)
         t1=a(1,4)
         t2=a(2,4)
         t3=a(3,4)
         a(4,4)=a(4,4)-a(4,1)*t1-a(4,2)*t2-a(4,3)*t3
         a(5,4)=a(5,4)-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3
         a(6,4)=a(6,4)-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3
!
         ip4=4
         pv=a(4,4)
         IF(abs(pv).lt.abs(a(5,4))) then
            ip4=5
            pv=a(5,4)
         ENDIF
         IF(abs(pv).lt.abs(a(6,4))) then
            ip4=6
            pv=a(6,4)
         ENDIF
!-----
         IF(ip4.ne.4) then
            t1=a(ip4,1)
            t2=a(ip4,2)
            t3=a(ip4,3)
            t5=a(ip4,5)
            t6=a(ip4,6)
            a(ip4,1)=a(4,1)
            a(ip4,2)=a(4,2)
            a(ip4,3)=a(4,3)
            a(ip4,4)=a(4,4)
            a(ip4,5)=a(4,5)
            a(ip4,6)=a(4,6)
            a(4,1)=t1
            a(4,2)=t2
            a(4,3)=t3
            a(4,5)=t5
            a(4,6)=t6
            DO l=1,m
               t=b(l,ip4,i)
               b(l,ip4,i)=b(l,4,i)
               b(l,4,i)=t
            ENDDO
         ENDIF
!-----
         pv=1.d0/pv
         a(4,4)=pv
         a(5,4)=a(5,4)*pv
         a(6,4)=a(6,4)*pv
         t1=a(4,1)
         t2=a(4,2)
         t3=a(4,3)
         a(4,5)=a(4,5)-t1*a(1,5)-t2*a(2,5)-t3*a(3,5)
         a(4,6)=a(4,6)-t1*a(1,6)-t2*a(2,6)-t3*a(3,6)
         t1=a(1,5)
         t2=a(2,5)
         t3=a(3,5)
         t4=a(4,5)
         a(5,5)=a(5,5)-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
         a(6,5)=a(6,5)-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4
!
         ip5=5
         pv=a(5,5)
         IF(abs(pv).lt.abs(a(6,5))) then
            pv=a(6,5)
            ip5=6
         ENDIF
!
         IF(ip5.ne.5) then
            t1=a(ip5,1)
            t2=a(ip5,2)
            t3=a(ip5,3)
            t4=a(ip5,4)
            t6=a(ip5,6)
            a(ip5,1)=a(5,1)
            a(ip5,2)=a(5,2)
            a(ip5,3)=a(5,3)
            a(ip5,4)=a(5,4)
            a(ip5,5)=a(5,5)
            a(ip5,6)=a(5,6)
            a(5,1)=t1
            a(5,2)=t2
            a(5,3)=t3
            a(5,4)=t4
            a(5,6)=t6
            DO l=1,m
               t=b(l,ip5,i)
               b(l,ip5,i)=b(l,5,i)
               b(l,5,i)=t
            ENDDO
         ENDIF
!
!-----
         pv=1.d0/pv
         a(5,5)=pv
         a(6,5)=a(6,5)*pv
         t1=a(1,6)
         t2=a(2,6)
         t3=a(3,6)
         t4=a(4,6)
         t5=a(5,6)
         t6=a(6,6)
         t5=t5-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
         a(5,6)=t5
         t6=t6-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4-a(6,5)*t5
         t6=1.d0/t6
         a(6,6)=t6
!-----
         DO l=1,m
            t1=b(l,1,i)
            t2=b(l,2,i)
            t3=b(l,3,i)
            t4=b(l,4,i)
            t5=b(l,5,i)
            t6=b(l,6,i)
            t2=t2-a(2,1)*t1
            t3=t3-a(3,1)*t1-a(3,2)*t2
            t4=t4-a(4,1)*t1-a(4,2)*t2-a(4,3)*t3
            t5=t5-a(5,1)*t1-a(5,2)*t2-a(5,3)*t3-a(5,4)*t4
            t6=t6-a(6,1)*t1-a(6,2)*t2-a(6,3)*t3-a(6,4)*t4-a(6,5)*t5
!
            b(l,6,i)=t6*a(6,6)
            t5=t5-a(5,6)*b(l,6,i)
            b(l,5,i)=t5*a(5,5)
            t4=t4-a(4,5)*b(l,5,i)-a(4,6)*b(l,6,i)
            b(l,4,i)=t4*a(4,4)
            t3=t3-a(3,4)*b(l,4,i)-a(3,5)*b(l,5,i)-a(3,6)*b(l,6,i)
            b(l,3,i)=t3*a(3,3)
            t2=t2-a(2,3)*b(l,3,i)-a(2,4)*b(l,4,i)-a(2,5)*b(l,5,i)-a(2,6)*b(l,6,i)
            b(l,2,i)=t2*a(2,2)
            t1=t1-a(1,2)*b(l,2,i)-a(1,3)*b(l,3,i)-a(1,4)*b(l,4,i)-a(1,5)*b(l,5,i)-a(1,6)*b(l,6,i)
            b(l,1,i)=t1*a(1,1)
         ENDDO
      ENDDO
!
      END SUBROUTINE luinverse6mn 
