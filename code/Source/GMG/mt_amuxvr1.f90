       subroutine mt_amuxvr1(n,nnz,x,b,y,a,ja,ia,iaa,ngroup,nbgroup)
       implicit none
!
       integer :: n,nnz
       integer :: ja(nnz)
       integer :: ia(2,ngroup+1)
       integer :: iaa(2,ngroup+1)
       integer :: ngroup,nbgroup(3,ngroup)
       real*8 :: x(n)
       real*8 :: b(iaa(2,ngroup+1)-1)
       real*8 :: y(iaa(2,ngroup+1)-1)
       real*8  :: a(nnz)
!
       integer :: i
       integer :: nn,nn1,nn2,npad
       integer :: i0,i1,i2
       integer :: j0,j1,j2,j3,j4,j5,j6,j7,j8
!  ...
!  ...
     npad=iaa(2,ngroup+1)-1
     do i=1,ngroup
         nn =nbgroup(1,i)
         nn1=nbgroup(2,i)
         nn2=nbgroup(3,i)
         i0=ia(2,i)
         i1=iaa(1,i)
         i2=iaa(2,i)
         j0=ia(2,i)
         j1=j0+  nn2
         j2=j0+2*nn2
         j3=j0+3*nn2
         j4=j0+4*nn2
         j5=j0+5*nn2
         j6=j0+6*nn2
         j7=j0+7*nn2
         j8=j0+8*nn2
!        write(*,101) i,nn,nn1,nn2,i0,i1
101   format(20(i5,2x))
         IF(nn.ge.32) THEN
          IF    (nn1.eq.1) THEN
            call mt_amuxrr10(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0), &
            ja(j0))
          ELSEIF(nn1.eq.2) THEN
            call mt_amuxrr20(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1), &
            ja(j0),ja(j1))
          ELSEIF(nn1.eq.3) THEN
            call mt_amuxrr3(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2), &
            ja(j0),ja(j1),ja(j2))
          ELSEIF(nn1.eq.4) THEN
            call mt_amuxrr4(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3), &
            ja(j0),ja(j1),ja(j2),ja(j3))
          ELSEIF(nn1.eq.5) THEN
            call mt_amuxrr5(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3),a(j4), &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4))
          ELSEIF(nn1.eq.6) THEN
            call mt_amuxrr6(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5), &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5))
          ELSEIF(nn1.eq.7) THEN
            call mt_amuxrr7(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6), &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6))
          ELSEIF(nn1.eq.8) THEN
            call mt_amuxrr8(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7), &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7))
          ELSEIF(nn1.eq.9) THEN
            call mt_amuxrr9(nn,nnz,x,b(i2),y(i2),nn2,n,npad, &
            a(j0),a(j1),a(j2),a(j3),a(j4),a(j5),a(j6),a(j7),a(j8), &
            ja(j0),ja(j1),ja(j2),ja(j3),ja(j4),ja(j5),ja(j6),ja(j7),ja(j8))
          ELSE
            call mt_amuxrr99(nn,nn1,nnz,x,b(i2),y(i2),a(i0),ja(i0),nn2,n,npad)
          ENDIF
         ELSE
          call mt_amuxrrn(nn,nn1,nnz,x,b(i2),y(i2),a(i0),ja(i0),nn2,n,npad)
         ENDIF
     enddo
!    stop 99
     return
     end subroutine
     
!DEC$ ATTRIBUTES INLINE :: amuxrr10
    subroutine mt_amuxrr10(n,nnz,x,b,y,nn,n0,n1, &
    a0, &
    ja0)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz)
    REAL(8) :: a0(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0
!
!DIR$ ASSUME_ALIGNED a0:64, &
                     ja0:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      y(i)= b(i)-a0(i)*x(j0)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr20
    subroutine mt_amuxrr20(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1, &
    ja0,ja1)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz)
    REAL(8) :: a0(nnz),a1(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64, &
                     ja0:64,ja1:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      y(i)= b(i)-a0(i)*x(j0)-a1(i)*x(j1)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr3
    subroutine mt_amuxrr3(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2, &
    ja0,ja1,ja2)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64, &
                     ja0:64,ja1:64,ja2:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      y(i)= b(i)-a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr4
    subroutine mt_amuxrr4(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3, &
    ja0,ja1,ja2,ja3)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD 
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      y(i)= b(i)-a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr5
    subroutine mt_amuxrr5(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3,a4, &
    ja0,ja1,ja2,ja3,ja4)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz),ja4(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz),a4(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3,j4
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64,a4:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64,ja4:64, &
                     y:64,b:64
!DIR$ UNROLL(2)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      j4=ja4(i)
      y(i)= b(i)-a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)-a4(i)*x(j4)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr6
    subroutine mt_amuxrr6(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3,a4,a5, &
    ja0,ja1,ja2,ja3,ja4,ja5)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz),ja4(nnz),ja5(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz),a4(nnz),a5(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3,j4,j5
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64,a4:64,a5:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64,ja4:64,ja5:64, &
                     y:64,b:64
!DIR$ UNROLL(2)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      j4=ja4(i)
      j5=ja5(i)
      y(i)= b(i)                                                        &
           -a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)-a4(i)*x(j4) &
           -a5(i)*x(j5)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr7
    subroutine mt_amuxrr7(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3,a4,a5,a6, &
    ja0,ja1,ja2,ja3,ja4,ja5,ja6)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz),ja4(nnz),ja5(nnz),ja6(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz),a4(nnz),a5(nnz),a6(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3,j4,j5,j6
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64,a4:64,a5:64,a6:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64,ja4:64,ja5:64,ja6:64, &
                     y:64,b:64
!DIR$ UNROLL(2)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      j4=ja4(i)
      j5=ja5(i)
      j6=ja6(i)
      y(i)= b(i)                                                        &
           -a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)-a4(i)*x(j4) &
           -a5(i)*x(j5)-a6(i)*x(j6)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr8
    subroutine mt_amuxrr8(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3,a4,a5,a6,a7, &
    ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz),ja4(nnz),ja5(nnz),ja6(nnz),ja7(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz),a4(nnz),a5(nnz),a6(nnz),a7(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3,j4,j5,j6,j7
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64,a4:64,a5:64,a6:64,a7:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64,ja4:64,ja5:64,ja6:64,ja7:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      j4=ja4(i)
      j5=ja5(i)
      j6=ja6(i)
      j7=ja7(i)
      y(i)= b(i)                                                        &
           -a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)-a4(i)*x(j4) &
           -a5(i)*x(j5)-a6(i)*x(j6)-a7(i)*x(j7)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr9
    subroutine mt_amuxrr9(n,nnz,x,b,y,nn,n0,n1, &
    a0,a1,a2,a3,a4,a5,a6,a7,a8, &
    ja0,ja1,ja2,ja3,ja4,ja5,ja6,ja7,ja8)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja0(nnz),ja1(nnz),ja2(nnz),ja3(nnz),ja4(nnz),ja5(nnz),ja6(nnz),ja7(nnz),ja8(nnz)
    REAL(8) :: a0(nnz),a1(nnz),a2(nnz),a3(nnz),a4(nnz),a5(nnz),a6(nnz),a7(nnz),a8(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i
    integer :: j0,j1,j2,j3,j4,j5,j6,j7,j8
!
!DIR$ ASSUME_ALIGNED a0:64,a1:64,a2:64,a3:64,a4:64,a5:64,a6:64,a7:64,a8:64, &
                     ja0:64,ja1:64,ja2:64,ja3:64,ja4:64,ja5:64,ja6:64,ja7:64,ja8:64, &
                     y:64,b:64
!DIR$ UNROLL(4)
!DIR$ SIMD
    do i=1,n
      j0=ja0(i)
      j1=ja1(i)
      j2=ja2(i)
      j3=ja3(i)
      j4=ja4(i)
      j5=ja5(i)
      j6=ja6(i)
      j7=ja7(i)
      j8=ja8(i)
      y(i)= b(i)                                                        &
           -a0(i)*x(j0)-a1(i)*x(j1)-a2(i)*x(j2)-a3(i)*x(j3)-a4(i)*x(j4) &
           -a5(i)*x(j5)-a6(i)*x(j6)-a7(i)*x(j7)-a8(i)*x(j8)
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrr99
    subroutine mt_amuxrr99(n,nn1,nnz,x,b,y,a,ja,nn,n0,n1)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja(nnz)
    REAL(8) :: a(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i,jj,j,j1
    REAL(8) temp
!
!DIR$ ASSUME_ALIGNED a:64,ja:64,b:64,y:64
!DIR$ SIMD
    do i=1,n
      j1=0
      temp = b(i) 
      do jj=1,nn1
        j=ja(i+j1)
       temp=temp-a(i+j1)*x(j)
       j1=j1+nn
      enddo
      y(i) = temp
    enddo
    return
    end subroutine
!DEC$ ATTRIBUTES INLINE :: amuxrrn
    subroutine mt_amuxrrn(n,nn1,nnz,x,b,y,a,ja,nn,n0,n1)
    implicit none
!
    integer :: n,nn1,nnz,nn,n0,n1
    integer :: ja(nnz)
    REAL(8) :: a(nnz)
    REAL(8) :: x(n0),y(n1),b(n1)
!
    integer :: i,jj,j,j1
    REAL(8) temp
!
    do i=1,n
      j1=0
      temp = b(i) 
!DIR$ NOVECTOR
      do jj=1,nn1
        j=ja(i+j1)
        temp=temp-a(i+j1)*x(j)
        j1=j1+nn
      enddo
      y(i) = temp
    enddo
    return
    end subroutine
