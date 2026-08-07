!
!-----------------------------------------------------------------------
!
      SUBROUTINE amux0(ncell,n,maxmt,x,y,a,ja,ia,iaa,ngroup,nbgroup)

      IMPLICIT NONE
!
!.....Input
      INTEGER :: ncell,n,maxmt
      INTEGER :: ngroup,nbgroup(3,ngroup)
      INTEGER,DIMENSION(2,ngroup+1) :: ia,iaa
      INTEGER,DIMENSION(maxmt) :: ja
      REAL(8),DIMENSION(maxmt) :: a
      REAL(8),DIMENSION(n) :: x
!.....Output
      REAL(8),DIMENSION(ncell) :: y
!.....Local variables
      INTEGER :: i,i0,i1,i2,nn,nn1,nn2
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11
!
      do i=1,ngroup
         nn =nbgroup(1,i)
         nn1=nbgroup(2,i)
         nn2=nbgroup(3,i)
         i0=ia(2,i)
         i1=iaa(1,i)
         i2=iaa(2,i)
         j0 =i0
         j1 =i0+   nn2
         j2 =i0+ 2*nn2
         j3 =i0+ 3*nn2
         j4 =i0+ 4*nn2
         j5 =i0+ 5*nn2
         j6 =i0+ 6*nn2
         j7 =i0+ 7*nn2
         j8 =i0+ 8*nn2
         j9 =i0+ 9*nn2
         j10=i0+10*nn2
         j11=i0+11*nn2
          IF    (nn1.eq.1) THEN
            CALL amux1(nn,maxmt,x,y(i2),ncell,n, &
            ja(j0),                              &
            a(j0))
          ELSEIF(nn1.eq.2) THEN
            CALL amux2(nn,maxmt,x,y(i2),ncell,n, &
            ja(j0),ja(j1),                       &
            a(j0),a(j1))
          ELSEIF(nn1.eq.3) THEN
            CALL amux3(nn,maxmt,x,y(i2),ncell,n, &
            ja(j0),ja(j1),ja(j2),                &
            a(j0),a(j1),a(j2))
          ELSEIF(nn1.eq.4) THEN
            CALL amux4(nn,maxmt,x,y(i2),ncell,n, &
            ja(j0),ja(j1),ja(j2),ja(j3),         &
            a(j0),a(j1),a(j2),a(j3))
          ELSEIF(nn1.eq.5) THEN
            CALL amux5(nn,maxmt,x,y(i2),ncell,n, &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),  &
            a(j0),a(j1),a(j2),a(j3),a(j4))
          ELSEIF(nn1.eq.6) THEN
            CALL amux6(nn,maxmt,x,y(i2),ncell,n,       &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5), &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5))
          ELSEIF(nn1.eq.7) THEN
            CALL amux7(nn,maxmt,x,y(i2),ncell,n,              &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6), &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6))
          ELSEIF(nn1.eq.8) THEN
            CALL amux8(nn,maxmt,x,y(i2),ncell,n,                     &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7), &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7))
          ELSEIF(nn1.eq.9) THEN
            CALL amux9(nn,maxmt,x,y(i2),ncell,n,                            &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7),ja(j8), &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7),a(j8))
          ELSEIF(nn1.eq.10) THEN
            CALL amux10(nn,maxmt,x,y(i2),ncell,n,                           &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7),ja(j8), &
            ja(j9),                                                         &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7),a(j8),          &
            a(j9))
          ELSEIF(nn1.eq.11) THEN
            CALL amux11(nn,maxmt,x,y(i2),ncell,n,                           &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7),ja(j8), &
            ja(j9),ja(j10),                                                 &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7),a(j8),          &
            a(j9),a(j10))
          ELSEIF(nn1.eq.12) THEN
            CALL amux12(nn,maxmt,x,y(i2),ncell,n,                           &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7),ja(j8), &
            ja(j9),ja(j10),ja(j11),                                         &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7),a(j8),          &
            a(j9),a(j10),a(j11))
          ELSE
            CALL amuxv(nn,nn1,nn2,maxmt,x,y(i2),a(i0),ja(i0),ncell,n)
          ENDIF
      ENDDO
!
      END SUBROUTINE amux0
!!DEC$ ATTRIBUTES INLINE :: amux1
      SUBROUTINE amux1(m,maxmt,x,y,n,nn, &
                       ja0,              &
                       a0)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0
      REAL(8),DIMENSION(maxmt) :: a0
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi
!!DIR$ ASSUME_ALIGNED a0:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         y(i)= a0(i)*x(j0)
      ENDDO
!
      END SUBROUTINE amux1
!!DEC$ ATTRIBUTES INLINE :: amux2
      SUBROUTINE amux2(m,maxmt,x,y,n,nn, &
                       ja0,ja1,          &
                       a0,a1)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1
      REAL(8),DIMENSION(maxmt) :: a0,a1
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)
      ENDDO
!
      END SUBROUTINE amux2
!!DEC$ ATTRIBUTES INLINE :: amux3
      SUBROUTINE amux3(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2, &
                       a0,a1,a2)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,    &
!                     y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)
      ENDDO
!
      END SUBROUTINE amux3
!!DEC$ ATTRIBUTES INLINE :: amux4
      SUBROUTINE amux4(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3, &
                       a0,a1,a2,a3)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,     &
!                     y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)
      ENDDO
!
      END SUBROUTINE amux4
!!DEC$ ATTRIBUTES INLINE :: amux5
      SUBROUTINE amux5(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3,ja4, &
                       a0,a1,a2,a3,a4)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         y(i)=a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4)
      ENDDO
!
      END SUBROUTINE amux5
!!DEC$ ATTRIBUTES INLINE :: amux6
      SUBROUTINE amux6(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3,ja4,ja5, &
                       a0,a1,a2,a3,a4,a5)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)
      ENDDO
!
      END SUBROUTINE amux6
!!DEC$ ATTRIBUTES INLINE :: amux7
      SUBROUTINE amux7(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3,ja4,ja5,ja6, &
                       a0,a1,a2,a3,a4,a5,a6)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)+a6(i)*x(j6)  
      ENDDO
!
      END SUBROUTINE amux7
!!DEC$ ATTRIBUTES INLINE :: amux8
      SUBROUTINE amux8(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7, &
                       a0,a1,a2,a3,a4,a5,a6,a7)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6,a7
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi,ja7:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,a7:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         j7=ja7(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)+a6(i)*x(j6)+a7(i)*x(j7)
      ENDDO
!
      END SUBROUTINE amux8
!!DEC$ ATTRIBUTES INLINE :: amux9
      SUBROUTINE amux9(m,maxmt,x,y,n,nn, &
                       ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8, &
                       a0,a1,a2,a3,a4,a5,a6,a7,a8)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7,j8
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi,ja7:avi,ja8:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,a7:avx,a8:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         j7=ja7(i)
         j8=ja8(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)+a6(i)*x(j6)+a7(i)*x(j7)+a8(i)*x(j8)
      ENDDO
!
      END SUBROUTINE amux9
!!DEC$ ATTRIBUTES INLINE :: amux10
      SUBROUTINE amux10(m,maxmt,x,y,n,nn, &
                        ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9, &
                        a0,a1,a2,a3,a4,a5,a6,a7,a8,a9)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7,j8,j9
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi,ja7:avi,ja8:avi,ja9:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,a7:avx,a8:avx,a9:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         j7=ja7(i)
         j8=ja8(i)
         j9=ja9(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)+a6(i)*x(j6)+a7(i)*x(j7)+a8(i)*x(j8)+a9(i)*x(j9)
      ENDDO
!
      END SUBROUTINE amux10
!!DEC$ ATTRIBUTES INLINE :: amux11
      SUBROUTINE amux11(m,maxmt,x,y,n,nn, &
                        ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9,ja10, &
                        a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9,ja10
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi,ja7:avi,ja8:avi,ja9:avi,ja10:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,a7:avx,a8:avx,a9:avx,a10:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         j7=ja7(i)
         j8=ja8(i)
         j9=ja9(i)
         j10=ja10(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
                +a5(i)*x(j5)+a6(i)*x(j6)+a7(i)*x(j7)+a8(i)*x(j8)+a9(i)*x(j9) &
                +a10(i)*x(j10)
      ENDDO
!
      END SUBROUTINE amux11
!!DEC$ ATTRIBUTES INLINE :: amux12
      SUBROUTINE amux12(m,maxmt,x,y,n,nn, &
                        ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9,ja10,ja11, &
                        a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)
      IMPLICIT NONE
!
      INTEGER :: m,nn,n,maxmt
      INTEGER,DIMENSION(maxmt) :: ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8,ja9,ja10,ja11
      REAL(8),DIMENSION(maxmt) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
      REAL(8),DIMENSION(nn) :: x
      REAL(8),DIMENSION(n) :: y
!
      INTEGER :: i
      INTEGER :: j0,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja0:avi,ja1:avi,ja2:avi,ja3:avi,ja4:avi,ja5:avi,ja6:avi,ja7:avi,ja8:avi,ja9:avi,ja10:avi,ja11:avi
!!DIR$ ASSUME_ALIGNED a0:avx,a1:avx,a2:avx,a3:avx,a4:avx,a5:avx,a6:avx,a7:avx,a8:avx,a9:avx,a10:avx,a11:avx,y:avx
!!DIR$ SIMD
      DO i=1,m
         j0=ja0(i)
         j1=ja1(i)
         j2=ja2(i)
         j3=ja3(i)
         j4=ja4(i)
         j5=ja5(i)
         j6=ja6(i)
         j7=ja7(i)
         j8=ja8(i)
         j9=ja9(i)
         j10=ja10(i)
         j11=ja11(i)
         y(i)= a0(i)*x(j0)+a1(i)*x(j1)+a2(i)*x(j2)+a3(i)*x(j3)+a4(i)*x(j4) &
              +a5(i)*x(j5)+a6(i)*x(j6)+a7(i)*x(j7)+a8(i)*x(j8)+a9(i)*x(j9) &
              +a10(i)*x(j10)+a11(i)*x(j11)
      ENDDO
!
      END SUBROUTINE amux12
!
!!DEC$ ATTRIBUTES INLINE :: amuxv
      SUBROUTINE amuxv(nn,nn1,nn2,maxmt,x,y,a,ja,n0,n1)
      IMPLICIT NONE
!
      INTEGER nn,nn1,nn2,maxmt,n0,n1
      INTEGER ja(maxmt)
      REAL(8) a(maxmt)
      REAL(8) x(n1),y(n0)
!
      INTEGER i,i1,jj
      INTEGER j0
      REAL(8) t0
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED ja:avi
!!DIR$ ASSUME_ALIGNED a:avx,y:avx
!!DIR$ SIMD
      DO i=1,nn
         i1=(i-1)*nn2
         t0=0.d0
         DO jj=1,nn1
            j0=ja(i1)
            t0=t0+a(i1)*x(j0)
!           i1=i1+nn2
         ENDDO
         y(i)=t0
      ENDDO
!
      END SUBROUTINE amuxv
!!DEC$ ATTRIBUTES INLINE :: amuxnn
      SUBROUTINE amuxnn(nn,nn1,nn2,maxmt,x,y,a,ja,jaa,n0,n1)
      IMPLICIT NONE
!
      INTEGER nn,nn1,nn2,maxmt,n0,n1
      INTEGER jaa(n1)
      INTEGER ja(maxmt)
      REAL(8) a(maxmt)
      REAL(8) x(n1),y(n0)
!
      INTEGER i,i1,jj
      INTEGER ip0,ip1,ip2,ip3
      INTEGER j0,j1,j2,j3
      REAL(8) t0,t1,t2,t3
!
      DO i=1,nn-3,4
         i1=i
         ip0=jaa(i)
         ip1=jaa(i+1)
         ip2=jaa(i+2)
         ip3=jaa(i+3)
         t0=0.d0
         t1=0.d0
         t2=0.d0
         t3=0.d0
!!DIR$ NOVECTOR
         DO jj=1,nn1
            j0=ja(i1)
            j1=ja(i1+1)
            j2=ja(i1+2)
            j3=ja(i1+3)
            t0=t0+a(i1  )*x(j0)
            t1=t1+a(i1+1)*x(j1)
            t2=t2+a(i1+2)*x(j2)
            t3=t3+a(i1+3)*x(j3)
            i1=i1+nn2
         ENDDO
         y(ip0)=t0
         y(ip1)=t1
         y(ip2)=t2
         y(ip3)=t3
      ENDDO
      IF(mod(nn,4).eq.1) THEN
         i=nn
         i1=i
         ip0=jaa(i)
         t0=0.d0
!!DIR$ NOVECTOR
         DO jj=1,nn1
            j0=ja(i1)
            t0=t0+a(i1  )*x(j0)
            i1=i1+nn2
         ENDDO
         y(ip0)=t0
      ELSEIF(mod(nn,4).eq.2) THEN
         i=nn-1
         i1=i
         ip0=jaa(i)
         ip1=jaa(i+1)
         t0=0.d0
         t1=0.d0
!!DIR$ NOVECTOR
         DO jj=1,nn1
            j0=ja(i1)
            j1=ja(i1+1)
            t0=t0+a(i1  )*x(j0)
            t1=t1+a(i1+1)*x(j1)
            i1=i1+nn2
         ENDDO
         y(ip0)=t0
         y(ip1)=t1
      ELSEIF(mod(nn,4).eq.3) THEN
         i=nn-2
         i1=i
         ip0=jaa(i)
         ip1=jaa(i+1)
         ip2=jaa(i+2)
         t0=0.d0
         t1=0.d0
         t2=0.d0
!!DIR$ NOVECTOR
         DO jj=1,nn1
            j0=ja(i1)
            j1=ja(i1+1)
            j2=ja(i1+2)
            t0=t0+a(i1  )*x(j0)
            t1=t1+a(i1+1)*x(j1)
            t2=t2+a(i1+2)*x(j2)
            i1=i1+nn2
         ENDDO
         y(ip0)=t0
         y(ip1)=t1
         y(ip2)=t2
      ENDIF
!
      END SUBROUTINE amuxnn
