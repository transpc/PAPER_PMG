      SUBROUTINE lusol0(n,n1,maxmt_lu0,maxmt_lu1,        &
                        diag,                              &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                        y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n1,maxmt_lu0,maxmt_lu1
      INTEGER :: ia0(n+1),ia1(n+1)
      INTEGER :: ja0(maxmt_lu0),ja1(maxmt_lu1)
      REAL(8) :: diag(n)
      REAL(8) :: diag_lu(n),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) :: y(n)
!.....Output
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,j,jj
      INTEGER :: j0,j1,j2,j3
      INTEGER :: l,m
      REAL(8) :: t
!
!.....forward solver
!
      DO i=1,n
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         IF(l.eq.0) THEN
            x(i)=y(i)*diag(i)
         ELSEIF(l.eq.1) THEN
            t=y(i)*diag(i)
            j0=ja0(jj)
            x(i)=t-alu0(jj)*x(j0)
         ELSEIF(l.eq.2) THEN
            t=y(i)*diag(i)
            j0=ja0(jj)
            j1=ja0(jj+1)
            x(i)=t-alu0(jj)*x(j0)-alu0(jj+1)*x(j1)
         ELSEIF(l.eq.3) THEN
            t=y(i)*diag(i)
            j0=ja0(jj)
            j1=ja0(jj+1)
            j2=ja0(jj+2)
            x(i)=t-alu0(jj)*x(j0)-alu0(jj+1)*x(j1)-alu0(jj+2)*x(j2)
         ELSEIF(l.eq.4) THEN
            t=y(i)*diag(i)
            j0=ja0(jj)
            j1=ja0(jj+1)
            j2=ja0(jj+2)
            j3=ja0(jj+3)
            x(i)=t-alu0(jj)*x(j0)-alu0(jj+1)*x(j1)-alu0(jj+2)*x(j2)-alu0(jj+3)*x(j3)
         ELSE
            t=y(i)*diag(i)
            DO m=0,l-1
               j=ja0(jj+m)
               t=t-alu0(jj+m)*x(j)
            ENDDO
            x(i)=t
         ENDIF
      ENDDO
!
!.....backward solver
!
      x(n)=x(n)*diag_lu(n)
      DO i=n-1,1,-1
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         IF(l.eq.0) THEN
            x(i)=x(i)*diag_lu(i)
         ELSEIF(l.eq.1) THEN
            j0=ja1(jj)
            x(i)=(x(i)-alu1(jj)*x(j0))*diag_lu(i)
         ELSEIF(l.eq.2) THEN
            j0=ja1(jj)
            j1=ja1(jj+1)
            x(i)=(x(i)-alu1(jj)*x(j0)-alu1(jj+1)*x(j1))*diag_lu(i)
         ELSEIF(l.eq.3) THEN
            j0=ja1(jj)
            j1=ja1(jj+1)
            j2=ja1(jj+2)
            x(i)=(x(i)-alu1(jj)*x(j0)-alu1(jj+1)*x(j1)-alu1(jj+2)*x(j2))*diag_lu(i)
         ELSE
            t=x(i)
            DO m=0,l-1
               j=ja1(jj+m)
               t=t-alu1(jj+m)*x(j)
            ENDDO
            x(i)=t*diag_lu(i)
         ENDIF
      ENDDO
!
      END SUBROUTINE lusol0
!
!-----------------------------------------------------------------------
!
      SUBROUTINE lusol0r(n,n1,maxmt_lu0,maxmt_lu1,        &
                         diag,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         y,x,perm)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER n,n1,maxmt_lu0,maxmt_lu1
      INTEGER ia0(n+1),ia1(n+1)
      INTEGER ja0(maxmt_lu0),ja1(maxmt_lu1)
      INTEGER perm(n)
      REAL(8) diag(n)
      REAL(8) diag_lu(n),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) y(n)
!.....Output
      REAL(8) x(n1)
!.....Local variables
      INTEGER i,i1,j,jj,ip
      INTEGER l,m
      REAL(8) t
!.....Local arrays
      REAL(8) x1(n)
!
!.....forward solver
!
      DO i=1,n
         ip=perm(i)
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         t=y(ip)*diag(ip)
         DO m=0,l-1
            j=ja0(jj+m)
            t=t-alu0(jj+m)*x1(j)
         ENDDO
         x1(i)=t
      ENDDO
!
!.....backward solver
!
      DO i=1,n
         i1=n-i+1
         ip=perm(i1)
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t=x1(i1)
         DO m=0,l-1
            j=ja1(jj+m)
            t=t-alu1(jj+m)*x1(j)
         ENDDO
         x1(i1)=t*diag_lu(i1)
         x(ip)=x1(i1)
      ENDDO
!
      END SUBROUTINE lusol0r
!-----------------------------------------------------------------------
!
      SUBROUTINE lusol0r1(n,maxmt1_lu0,maxmt1_lu1,y,x,diag_lu,alu0,alu1,ja0,ja1,ia0,ia1,perm)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER n,maxmt1_lu0,maxmt1_lu1
      INTEGER ia0(n+1),ia1(n+1)
      INTEGER ja0(maxmt1_lu0),ja1(maxmt1_lu1)
      INTEGER perm(n)
      REAL(8) alu0(maxmt1_lu0),alu1(maxmt1_lu1)
      REAL(8) diag_lu(n)
      REAL(8) y(n)
!.....Output
      REAL(8) x(n)
!.....Local variables
      INTEGER i,i1,j,jj,ip
      INTEGER l,m
      REAL(8) t,diag
!.....Local arrays
      REAL(8) x1(n)
!
!.....forward solver
!
      DO i=1,n
         ip=perm(i)
         t=y(ip)
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         DO m=0,l-1
            j=ja0(jj+m)
            t=t-alu0(jj+m)*x1(j)
         ENDDO
         x1(i)=t
      ENDDO
!
!.....backward solver
!
      DO i=1,n
         i1=n-i+1
         ip=perm(i1)
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         diag=diag_lu(i1)
         t=x1(i1)
         DO m=0,l-1
            j=ja1(jj+m)
            t=t-alu1(jj+m)*x1(j)
         ENDDO
         x1(i1)=t*diag
         x(ip)=x1(i1)
      ENDDO
!
      END SUBROUTINE lusol0r1
!
      SUBROUTINE lusol00(n,maxmt_lu0,maxmt_lu1,           &
                         diag,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         n1,y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER n,n1,maxmt_lu0,maxmt_lu1
      INTEGER ia0(n+1),ia1(n+1)
      INTEGER ja0(maxmt_lu0),ja1(maxmt_lu1)
      REAL(8) diag(n)
      REAL(8) diag_lu(n),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) y(n1)
!.....Output
      REAL(8) x(n1)
!.....Local variables
      INTEGER i,j,jj
      INTEGER l,m
      REAL(8) t
!
!.....forward solver
!
      DO i=1,n
         t=y(i)*diag(i)
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         DO m=0,l-1
            j=ja0(jj+m)
            t=t-alu0(jj+m)*x(j)
         ENDDO
         x(i)=t
      ENDDO
!
!.....backward solver
!
      x(n)=x(n)*diag_lu(n)
      DO i=n-1,1,-1
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t=x(i)
         DO m=0,l-1
            j=ja1(jj+m)
            t=t-alu1(jj+m)*x(j)
         ENDDO
         x(i)=t*diag_lu(i)
      ENDDO
!
      END SUBROUTINE lusol00
!
      SUBROUTINE lusol00r(n,maxmt_lu0,maxmt_lu1,           &
                          diag,                              &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                          n1,y,x,perm)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER n,n1,maxmt_lu0,maxmt_lu1
      INTEGER ia0(n+1),ia1(n+1)
      INTEGER ja0(maxmt_lu0),ja1(maxmt_lu1)
      INTEGER perm(n)
      REAL(8) diag(n)
      REAL(8) diag_lu(n),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) y(n)
!.....Output
      REAL(8) x(n1)
!.....Local variables
      INTEGER i,j,jj,ip
      INTEGER l,m
      REAL(8) t
!.....Local arrays
      REAL(8) x1(n)
!
!.....forward solver
!
      DO i=1,n
         ip=perm(i)
         t=y(ip)*diag(ip)
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         DO m=0,l-1
            j=ja0(jj+m)
            t=t-alu0(jj+m)*x1(j)
         ENDDO
         x1(i)=t
      ENDDO
!
!.....backward solver
!
      ip=perm(n)
      x1(n)=x1(n)*diag_lu(n)
      x(ip)=x1(n)
      DO i=n-1,1,-1
         ip=perm(i)
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t=x1(i)
         DO m=0,l-1
            j=ja1(jj+m)
            t=t-alu1(jj+m)*x1(j)
         ENDDO
         x1(i)=t*diag_lu(i)
         x(ip)=x1(i)
      ENDDO
!
      END SUBROUTINE lusol00r
!
      SUBROUTINE lusol_tri(n,                 &
                           diag_lu,alu0,alu1, &
                           n1,y,x)
!
      IMPLICIT NONE
!.....Input
      INTEGER n,n1
      REAL(8) diag_lu(n),alu0(n),alu1(n)
      REAL(8) y(n)
!.....Output
      REAL(8) x(n1)
!.....Local variables
      INTEGER i
      REAL(8) t
!
      x(1)=y(1)
      DO i=2,n
         t=alu0(i)
         x(i)=y(i)-t*x(i-1)
      ENDDO
      x(n)=x(n)*diag_lu(n)
      DO i=n-1,1,-1
         x(i)=(x(i)-alu1(i)*x(i+1))*diag_lu(i)
      ENDDO
!
      END SUBROUTINE lusol_tri
