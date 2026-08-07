!
      SUBROUTINE rv_tdiag(a,b,c,d,nn)
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
      ENDSUBROUTINE rv_tdiag
!
!
      SUBROUTINE rv_tdiag_2d(a,b,c,d,n,nn)
!
      IMPLICIT NONE
!     input
      INTEGER n
      REAL(8) a(n,nn),b(n,nn),c(n,nn)
!     output
      REAL(8) d(n,nn)
!
      INTEGER nn,i,j,k
      REAL(8) xm
!
      IF(nn.eq.1)THEN
         DO j=1,n
            d(j,1)=d(j,1)/b(j,1)
         ENDDO 
         RETURN
      ENDIF
!
      DO j=1,n
         b(j,1)=1.d0/b(j,1)
      ENDDO
      DO i=2,nn
         k=i-1
!DIR$ SIMD
         DO j=1,n
            xm=a(j,i)*b(j,k)
            b(j,i)=1.d0/(b(j,i)-xm*c(j,k))
            d(j,i)=d(j,i)-xm*d(j,k)
         ENDDO
      ENDDO
!
!DIR$ SIMD
      DO j=1,n
         d(j,nn)=d(j,nn)*b(j,nn)
      ENDDO
      DO i=nn-1,1,-1
!DIR$ SIMD
         DO j=1,n
            d(j,i)=(d(j,i)-c(j,i)*d(j,i+1))*b(j,i)
         ENDDO
      ENDDO
!
      RETURN
      ENDSUBROUTINE rv_tdiag_2d
