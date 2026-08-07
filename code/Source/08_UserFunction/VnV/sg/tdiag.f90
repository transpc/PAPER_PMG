!
      SUBROUTINE tdiag(a,b,c,d,nn)
!
      IMPLICIT NONE
!
      INTEGER nn,i,k
      REAL(8) a(nn),b(nn),c(nn),d(nn)
      REAL(8) xm
!
      IF(nn.eq.1)THEN
         d(1)=d(1)/b(1)
         RETURN
      ENDIF
!
      DO i=2,nn
         k=i-1
         IF(b(k).eq.0.0d0) THEN
            PRINT *, '### Zero diagonal coefficient in tridiag.'
            STOP
         ENDIF
         xm=a(i)/b(k)
         b(i)=b(i)-xm*c(k)
         d(i)=d(i)-xm*d(k)
      ENDDO
!
      d(nn)=d(nn)/b(nn)
!
      DO i=2,nn
         k=nn+1-i
         d(k)=(d(k)-c(k)*d(k+1))/b(k)
      ENDDO
!
      RETURN
      ENDSUBROUTINE tdiag
!