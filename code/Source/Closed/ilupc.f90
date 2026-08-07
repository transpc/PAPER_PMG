!
      SUBROUTINE ilupc10(ncell,maxmt,maxmt_lu0,maxmt_lu1, &
                         ia,ja,ju,au,                      &
                         diag_lu,alu0,alu1,ia0,ia1)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: ncell,maxmt,maxmt_lu0,maxmt_lu1
      INTEGER :: ia(ncell+1),ju(ncell)
      INTEGER :: ja(maxmt)
      INTEGER :: ia0(ncell+1),ia1(ncell+1)
      REAL(8) :: au(maxmt)
!.....Output
      REAL(8) :: diag_lu(ncell)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: jj,kk
      INTEGER :: i1,k0,k1
      REAL(8) :: t
!.....Local arrays
      REAL(8) :: vi(ncell)
!
      DO i=1,ncell
         DO jj=ia(i),ia(i+1)-1
            j=ja(jj)
            vi(j)=au(jj)
         ENDDO
!
         DO jj=ia(i),ju(i)-1
            j=ja(jj)
            i1=ncell-j+1
            k1=ia1(i1)
            t=vi(j)*diag_lu(j)
            vi(j)=t
!DIR$ IVDEP
            DO kk=ju(j)+1,ia(j+1)-1 
               k=ja(kk)
               vi(k)=vi(k)-t*alu1(k1)
               k1=k1+1
            ENDDO
         ENDDO
!
         k0=ia0(i)
         i1=ncell-i+1
         k1=ia1(i1)
         DO jj=ia(i),ju(i)-1
            j=ja(jj)
            alu0(k0)=vi(j)
            vi(j)=0.d0
            k0=k0+1
         ENDDO 
            jj=ju(i)
            j=ja(jj)
            diag_lu(i)=1.d0/vi(j)
            vi(j)=0.d0
         DO jj=ju(i)+1,ia(i+1)-1
            j=ja(jj)
            alu1(k1)=vi(j)
            vi(j)=0.d0
            k1=k1+1
         ENDDO 
      ENDDO
!
      END SUBROUTINE ilupc10
!
      SUBROUTINE ilupc11(ncell,maxmt,maxmt_lu0,maxmt_lu1,  &
                         diag,au,ia,ja,ju,iend,             &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: ncell,maxmt,maxmt_lu0,maxmt_lu1
      INTEGER :: ia(ncell+1),ju(ncell),iend(ncell)
      INTEGER :: ja(maxmt)
      INTEGER :: ia0(ncell+1),ia1(ncell+1)
      REAL(8) :: diag(ncell)
      REAL(8) :: au(maxmt)
!.....Output
      INTEGER :: ja0(maxmt_lu0),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(ncell)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: jj,kk
      INTEGER :: k0,k1
      REAL(8) :: t
!.....Local arrays
      REAL(8) :: vi(ncell)
!
      ia0(1)=1
      ia1(1)=1
      DO i=1,ncell
         DO jj=ia(i),iend(i)
            j=ja(jj)
            vi(j)=au(jj)*diag(i)
         ENDDO
         vi(i)=1.d0
!
         k0=ia0(i)
         DO jj=ia(i),ju(i)-1
            j=ja(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
            vi(j)=0.d0
            k0=k0+1
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO 
         ia0(i+1)=k0
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju(i)+1,iend(i)
            j=ja(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
      END SUBROUTINE ilupc11
!
      SUBROUTINE ilupc11p(n,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                          diag,au,ia_a,ja_a,iend,               &
                          ia_r,ja_r,ju_r,                       &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,    &
                          perm,permi)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n)
      REAL(8) :: diag(n),au(maxmt)
!.....Output
      INTEGER :: ia0(n+1),ja0(maxmt_lu0),ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
!.....Local variables
      INTEGER :: i,j,k,i1,j1
      INTEGER :: jj,kk
      INTEGER :: k0,k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,diagr
!.....Local arrays
      REAL(8) :: vi(n)
!
      ip=1
      ia0(1)=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            vi(j1)=au(jj)*diagr
            ip=ip+1
         ENDDO
         vi(i)=1.d0
!
         k0=ia0(i)
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
            vi(j)=0.d0
            k0=k0+1
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1 
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         ia0(i+1)=k0
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
      END SUBROUTINE ilupc11p
!
      SUBROUTINE factor_solve_tri(n,maxmt,ia,au,    &
                                  diag_lu,alu0,alu1, &
                                  n1,y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n1,maxmt
      INTEGER :: ia(n+1)
      REAL(8) :: au(maxmt)
      REAL(8) :: y(n)
!.....Output
      REAL(8) :: diag_lu(n),alu0(n),alu1(n)
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,jj
      REAL(8) :: t
!
      diag_lu(1)=1.d0/au(1)
      x(1)=y(1)
      DO i=2,n
         jj=ia(i)
         t=au(jj)*diag_lu(i-1)
         diag_lu(i)=1.d0/(au(jj+1)-t*au(jj-1))
         alu0(i)=t
         alu1(i-1)=au(jj-1)
         x(i)=y(i)-t*x(i-1)
      ENDDO
      x(n)=x(n)*diag_lu(n)
      DO i=n-1,1,-1
         x(i)=(x(i)-alu1(i)*x(i+1))*diag_lu(i)
      ENDDO
!
      END SUBROUTINE factor_solve_tri
!
      SUBROUTINE factor_solve_tri0(n,maxmt,ia,au, &
                                   n1,y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n1,maxmt
      INTEGER :: ia(n+1)
      REAL(8) :: au(maxmt)
      REAL(8) :: y(n)
!.....Output
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,jj
      REAL(8) :: t
      REAL(8) :: diag(n)
!
      diag(1)=1.d0/au(1)
      x(1)=y(1)
      DO i=2,n
         jj=ia(i)
         t=au(jj)*diag(i-1)
         diag(i)=1.d0/(au(jj+1)-t*au(jj-1))
         x(i)=y(i)-t*x(i-1)
      ENDDO
      x(n)=x(n)*diag(n)
      DO i=n-1,2,-1
         jj=ia(i)
         x(i)=(x(i)-au(jj+2)*x(i+1))*diag(i)
      ENDDO
      x(1)=(x(1)-au(2)*x(2))*diag(1)
!
      END SUBROUTINE factor_solve_tri0
!
      SUBROUTINE factor_solve_tri02(n,n_pad,maxmt,ia,au, &
                                    n1,y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_pad,n1,maxmt
      INTEGER :: ia(n+1)
      REAL(8) :: au(maxmt)
      REAL(8) :: y(n_pad,2)
!.....Output
      REAL(8) :: x(n1,2)
!.....Local variables
      INTEGER :: i,jj
      REAL(8) :: t
      REAL(8) :: diag(n)
!
      diag(1)=1.d0/au(1)
      x(1,1)=y(1,1)
      x(1,2)=y(1,2)
      DO i=2,n
         jj=ia(i)
         t=au(jj)*diag(i-1)
         diag(i)=1.d0/(au(jj+1)-t*au(jj-1))
         x(i,1)=y(i,1)-t*x(i-1,1)
         x(i,2)=y(i,2)-t*x(i-1,2)
      ENDDO
!
      x(n,1)=x(n,1)*diag(n)
      x(n,2)=x(n,2)*diag(n)
      DO i=n-1,2,-1
         jj=ia(i)
         x(i,1)=(x(i,1)-au(jj+2)*x(i+1,1))*diag(i)
         x(i,2)=(x(i,2)-au(jj+2)*x(i+1,2))*diag(i)
      ENDDO
      x(1,1)=(x(1,1)-au(2)*x(2,1))*diag(1)
      x(1,2)=(x(1,2)-au(2)*x(2,2))*diag(1)
!
      END SUBROUTINE factor_solve_tri02
!
      SUBROUTINE factor_solve_tri03(n,n_pad,maxmt,ia,au, &
                                    n1,y,x)
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,n_pad,n1,maxmt
      INTEGER :: ia(n+1)
      REAL(8) :: au(maxmt)
      REAL(8) :: y(n_pad,3)
!.....Output
      REAL(8) :: x(n1,3)
!.....Local variables
      INTEGER :: i,jj
      REAL(8) :: t
      REAL(8) :: diag(n)
!
      diag(1)=1.d0/au(1)
      x(1,1)=y(1,1)
      x(1,2)=y(1,2)
      x(1,3)=y(1,3)
      DO i=2,n
         jj=ia(i)
         t=au(jj)*diag(i-1)
         diag(i)=1.d0/(au(jj+1)-t*au(jj-1))
         x(i,1)=y(i,1)-t*x(i-1,1)
         x(i,2)=y(i,2)-t*x(i-1,2)
         x(i,3)=y(i,3)-t*x(i-1,3)
      ENDDO
!
      x(n,1)=x(n,1)*diag(n)
      x(n,2)=x(n,2)*diag(n)
      x(n,3)=x(n,3)*diag(n)
      DO i=n-1,2,-1
         jj=ia(i)
         x(i,1)=(x(i,1)-au(jj+2)*x(i+1,1))*diag(i)
         x(i,2)=(x(i,2)-au(jj+2)*x(i+1,2))*diag(i)
         x(i,3)=(x(i,3)-au(jj+2)*x(i+1,3))*diag(i)
      ENDDO
      x(1,1)=(x(1,1)-au(2)*x(2,1))*diag(1)
      x(1,2)=(x(1,2)-au(2)*x(2,2))*diag(1)
      x(1,3)=(x(1,3)-au(2)*x(2,3))*diag(1)
!
      END SUBROUTINE factor_solve_tri03
!
      SUBROUTINE factor(n,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                        diag,au,ia_a,ja_a,iend,               &
                        ia_r,ja_r,ju_r,                       &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      REAL(8) :: diag(n),au(maxmt)
!.....Output
      INTEGER :: ia0(n+1),ja0(maxmt_lu0),ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
!.....Local variables
      INTEGER :: i,j,k,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k0,k1
      INTEGER :: l,m
      REAL(8) :: t,diagr
!.....Local arrays
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
!
      ia0(1)=1
      ia1(1)=1
      DO i=1,n
         diagr=diag(i)
         j1=ia_a(i)
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.iend(i)) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja_a(j1)) THEN
               vi(j)=au(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         k0=ia0(i)
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
            vi(j)=0.d0
            k0=k0+1
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         ia0(i+1)=k0
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO
         ia1(i+1)=k1
      ENDDO
!
      END SUBROUTINE factor
!
      SUBROUTINE factor_solve(n,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                              diag,au,ia_a,ja_a,iend,               &
                              ia_r,ja_r,ju_r,                       &
                              diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,    &
                              n1,y,x)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n1,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n)
!.....Output
      INTEGER :: ia0(n+1),ja0(maxmt_lu0),ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,j,k,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k0,k1
      INTEGER :: l,m
      REAL(8) :: t,diagr
!.....Local arrays
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
!
      ia0(1)=1
      ia1(1)=1
      DO i=1,n
         diagr=diag(i)
         j1=ia_a(i)
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.iend(i)) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja_a(j1)) THEN
               vi(j)=au(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         k0=ia0(i)
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
            vi(j)=0.d0
            k0=k0+1
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO
         ia0(i+1)=k0
!
!.....forward solve
!
         t=y(i)*diagr
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         DO m=0,l-1
            j=ja0(jj+m)
            t=t-alu0(jj+m)*x(j)
         ENDDO
         x(i)=t
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
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
      END SUBROUTINE factor_solve
!
      SUBROUTINE factor_solve0(n,maxmt,maxmt_r,maxmt_lu1, &
                               diag,au,ia_a,ja_a,iend,     &
                               ia_r,ja_r,ju_r,             &
                               n1,y,x)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n1,maxmt,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n)
!.....Output
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,j,k,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1
      INTEGER :: l,m
      REAL(8) :: t,t1,diagr
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
!
      ia1(1)=1
      DO i=1,n
         diagr=diag(i)
         j1=ia_a(i)
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.iend(i)) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja_a(j1)) THEN
               vi(j)=au(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(i)*diagr
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x(j)
            vi(j)=0.d0
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         x(i)=t1
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
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
      END SUBROUTINE factor_solve0
!
      SUBROUTINE factor_solvev2(n,n_pad,maxmt,maxmt_r,maxmt_lu1, &
                                diag,au,ia_a,ja_a,iend,          &
                                ia_r,ja_r,ju_r,                  &
                                n1,y,x)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n_pad,n1,maxmt,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n_pad,2)
!.....Output
      REAL(8) :: x(n1,2)
!.....Local variables
      INTEGER :: i,j,k,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1
      INTEGER :: l,m
      REAL(8) :: t,diagr
      REAL(8) :: t1,t2
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
!
      ia1(1)=1
      DO i=1,n
         diagr=diag(i)
         j1=ia_a(i)
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.iend(i)) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja_a(j1)) THEN
               vi(j)=au(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(i,1)*diagr
         t2=y(i,2)*diagr
!
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x(j,1)
            t2=t2-t*x(j,2)
            vi(j)=0.d0
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         x(i,1)=t1
         x(i,2)=t2
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
!
      x(n,1)=x(n,1)*diag_lu(n)
      x(n,2)=x(n,2)*diag_lu(n)
      DO i=n-1,1,-1
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t1=x(i,1)
         t2=x(i,2)
         DO m=0,l-1
            j=ja1(jj+m)
            t1=t1-alu1(jj+m)*x(j,1)
            t2=t2-alu1(jj+m)*x(j,2)
         ENDDO
         x(i,1)=t1*diag_lu(i)
         x(i,2)=t2*diag_lu(i)
      ENDDO
!
      END SUBROUTINE factor_solvev2
!
      SUBROUTINE factor_solvev3(n,n_pad,maxmt,maxmt_r,maxmt_lu1, &
                                diag,au,ia_a,ja_a,iend,     &
                                ia_r,ja_r,ju_r,             &
                                n1,y,x)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n_pad,n1,maxmt,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n_pad,3)
!.....Output
      REAL(8) :: x(n1,3)
!.....Local variables
      INTEGER :: i,j,k,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1
      INTEGER :: l,m
      REAL(8) :: t,diagr
      REAL(8) :: t1,t2,t3
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
!
      ia1(1)=1
      DO i=1,n
         diagr=diag(i)
         j1=ia_a(i)
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.iend(i)) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja_a(j1)) THEN
               vi(j)=au(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(i,1)*diagr
         t2=y(i,2)*diagr
         t3=y(i,3)*diagr
!
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x(j,1)
            t2=t2-t*x(j,2)
            t3=t3-t*x(j,3)
            vi(j)=0.d0
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         x(i,1)=t1
         x(i,2)=t2
         x(i,3)=t3
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
!
      x(n,1)=x(n,1)*diag_lu(n)
      x(n,2)=x(n,2)*diag_lu(n)
      x(n,3)=x(n,3)*diag_lu(n)
      DO i=n-1,1,-1
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t1=x(i,1)
         t2=x(i,2)
         t3=x(i,3)
         DO m=0,l-1
            j=ja1(jj+m)
            t1=t1-alu1(jj+m)*x(j,1)
            t2=t2-alu1(jj+m)*x(j,2)
            t3=t3-alu1(jj+m)*x(j,3)
         ENDDO
         x(i,1)=t1*diag_lu(i)
         x(i,2)=t2*diag_lu(i)
         x(i,3)=t3*diag_lu(i)
      ENDDO
!
      END SUBROUTINE factor_solvev3
!
      SUBROUTINE factorp(n,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                         diag,au,ia_a,ja_a,iend,                      &
                         ia_r,ja_r,ju_r,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,           &
                         perm,permi,index)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n),index(maxmt2)
      REAL(8) :: diag(n),au(maxmt)
!.....Output
      INTEGER :: ia0(n+1),ja0(maxmt_lu0),ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
!.....Local variables
      INTEGER :: i,j,k,i1,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k0,k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,diagr
!.....Local arrays
      INTEGER :: ja2(n)
      REAL(8) :: vi(n)
      REAL(8) :: a1(n)
!
      ip=1
      ia0(1)=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         l=iend(i1)-ia_a(i1)+1
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja2(j1)) THEN
               vi(j)=a1(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         k0=ia0(i)
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
            vi(j)=0.d0
            k0=k0+1
            l=ia1(j+1)-ia1(j)
            kk=ia1(j)
!DIR$ IVDEP
            DO m=0,l-1 
               k=ja1(kk+m)
               vi(k)=vi(k)-t*alu1(kk+m)
            ENDDO
         ENDDO
         ia0(i+1)=k0
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
      END SUBROUTINE factorp
!
      SUBROUTINE factor_solvep(n,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                               diag,au,ia_a,ja_a,iend,                      &
                               ia_r,ja_r,ju_r,                              &
                               diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,           &
                               n1,y,x,perm,permi,index)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n1,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ia_r(n+1),ju_r(n),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n),index(maxmt2)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n)
!.....Output
      INTEGER :: ia0(n+1),ja0(maxmt_lu0),ia1(n+1),ja1(maxmt_lu1)
      REAL(8) :: diag_lu(n)
      REAL(8) :: alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,j,k,i1,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k0,k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,t1,diagr
!.....Local arrays
      INTEGER :: ja2(n)
      REAL(8) :: vi(n)
      REAL(8) :: x1(n)
      REAL(8) :: a1(n)
!
      ip=1
      ia0(1)=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         l=iend(i1)-ia_a(i1)+1
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja2(j1)) THEN
               vi(j)=a1(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
!        t1=y(ip1)*diagr
         k0=ia0(i)
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            ja0(k0)=j
            alu0(k0)=t
!           t1=t1-t*x1(j)
            vi(j)=0.d0
            k0=k0+1
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1 
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO
         ia0(i+1)=k0
!
!.....forward solve
!
         l=ia0(i+1)-ia0(i)
         jj=ia0(i)
         t1=y(ip1)*diagr
!DIR$ IVDEP
         DO m=0,l-1
            j=ja0(jj+m)
            t1=t1-alu0(jj+m)*x1(j)
         ENDDO
         x1(i)=t1
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
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
      END SUBROUTINE factor_solvep
!
      SUBROUTINE factor_solve0p(n,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                diag,au,ia_a,ja_a,iend,            &
                                ia_r,ja_r,ju_r,                    &
                                n1,y,x,perm,permi,index)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n1,maxmt,maxmt2,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),iend(n),ja_a(maxmt)
      INTEGER :: ju_r(n),ia_r(n+1),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n),index(maxmt2)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n)
!.....Output
      REAL(8) :: x(n1)
!.....Local variables
      INTEGER :: i,j,k,i1,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,t1,diagr
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      INTEGER :: ja2(n)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
      REAL(8) :: x1(n)
      REAL(8) :: a1(n)
!
      ip=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         l=iend(i1)-ia_a(i1)+1
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja2(j1)) THEN
               vi(j)=a1(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(ip1)*diagr
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x1(j)
            vi(j)=0.d0
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1 
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO
         x1(i)=t1
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
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
      END SUBROUTINE factor_solve0p
!
      SUBROUTINE factor_solve0pv2(n,n_pad,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                  diag,au,ia_a,ja_a,iend,                 &
                                  ia_r,ja_r,ju_r,                         &
                                  n1,y,x,perm,permi,index)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n_pad,n1,maxmt,maxmt2,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),ja_a(maxmt),iend(n)
      INTEGER :: ju_r(n),ia_r(n+1),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n),index(maxmt2)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n_pad,2)
!.....Output
      REAL(8) :: x(n1,2)
!.....Local variables
      INTEGER :: i,j,k,i1,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,diagr
      REAL(8) :: t1,t2
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      INTEGER :: ja2(n)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
      REAL(8) :: x1(n,2)
      REAL(8) :: a1(n)
!
      ip=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         l=iend(i1)-ia_a(i1)+1
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja2(j1)) THEN
               vi(j)=a1(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(ip1,1)*diagr
         t2=y(ip1,2)*diagr
!
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x1(j,1)
            t2=t2-t*x1(j,2)
            vi(j)=0.d0
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1 
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO
         x1(i,1)=t1
         x1(i,2)=t2
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
!
      ip=perm(n)
      x1(n,1)=x1(n,1)*diag_lu(n)
      x1(n,2)=x1(n,2)*diag_lu(n)
      x(ip,1)=x1(n,1)
      x(ip,2)=x1(n,2)
      DO i=n-1,1,-1
         ip=perm(i)
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t1=x1(i,1)
         t2=x1(i,2)
         DO m=0,l-1
            j=ja1(jj+m)
            t1=t1-alu1(jj+m)*x1(j,1)
            t2=t2-alu1(jj+m)*x1(j,2)
         ENDDO
         x1(i,1)=t1*diag_lu(i)
         x1(i,2)=t2*diag_lu(i)
         x(ip,1)=x1(i,1)
         x(ip,2)=x1(i,2)
      ENDDO
!
      END SUBROUTINE factor_solve0pv2
!
      SUBROUTINE factor_solve0pv3(n,n_pad,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                  diag,au,ia_a,ja_a,iend,                 &
                                  ia_r,ja_r,ju_r,                         &
                                  n1,y,x,perm,permi,index)
!
      IMPLICIT NONE
!
      REAL(8),PARAMETER :: threshold=1.d-20
!.....Input
      INTEGER :: n,n_pad,n1,maxmt,maxmt2,maxmt_r,maxmt_lu1
      INTEGER :: ia_a(n+1),ja_a(maxmt),iend(n)
      INTEGER :: ju_r(n),ia_r(n+1),ja_r(maxmt_r)
      INTEGER :: perm(n),permi(n),index(maxmt2)
      REAL(8) :: diag(n),au(maxmt)
      REAL(8) :: y(n_pad,3)
!.....Output
      REAL(8) :: x(n1,3)
!.....Local variables
      INTEGER :: i,j,k,i1,j1,j2
      INTEGER :: jj,kk
      INTEGER :: k1,ip,ip1
      INTEGER :: l,m
      REAL(8) :: t,diagr
      REAL(8) :: t1,t2,t3
!.....Local arrays
      INTEGER :: ia1(n+1),ja1(maxmt_lu1)
      INTEGER :: ja2(n)
      REAL(8) :: diag_lu(n),alu1(maxmt_lu1)
!     REAL(8) :: vi(n),vii(n)
      REAL(8) :: vi(n)
      REAL(8) :: x1(n,3)
      REAL(8) :: a1(n)
!
      ip=1
      ia1(1)=1
      DO i=1,n
         ip1=perm(i)
         diagr=diag(ip1)
         i1=perm(i)
         l=iend(i1)-ia_a(i1)+1
         DO jj=ia_a(i1),iend(i1)
            j=ja_a(jj)
            j1=permi(j)
            k=index(ip)
            a1(k)=au(jj)
            ja2(k)=j1
            ip=ip+1
         ENDDO
         j1=1
         DO jj=ia_r(i),ia_r(i+1)-1
            IF(j1.gt.l) THEN
               DO j2=jj,ia_r(i+1)-1
                  j=ja_r(j2)
                  vi(j)=0.d0
               ENDDO
               exit
            ENDIF
            j=ja_r(jj)
            IF(ja_r(jj).eq.ja2(j1)) THEN
               vi(j)=a1(j1)*diagr
               j1=j1+1
            ELSE
               vi(j)=0.d0
            ENDIF
         ENDDO
!
         t1=y(ip1,1)*diagr
         t2=y(ip1,2)*diagr
         t3=y(ip1,3)*diagr
!
         DO jj=ia_r(i),ju_r(i)-1
            j=ja_r(jj)
            t=vi(j)*diag_lu(j)
            IF(abs(t).lt.threshold) THEN
               vi(j)=0.d0
               cycle
            ENDIF
            t1=t1-t*x1(j,1)
            t2=t2-t*x1(j,2)
            t3=t3-t*x1(j,3)
            vi(j)=0.d0
!DIR$ IVDEP
            DO kk=ia1(j),ia1(j+1)-1 
               k=ja1(kk)
               vi(k)=vi(k)-t*alu1(kk)
            ENDDO
         ENDDO
         x1(i,1)=t1
         x1(i,2)=t2
         x1(i,3)=t3
!
         diag_lu(i)=1.d0/vi(i)
!
         k1=ia1(i)
         DO jj=ju_r(i)+1,ia_r(i+1)-1
            j=ja_r(jj)
            IF(abs(vi(j)).ge.threshold) THEN
               ja1(k1)=j
               alu1(k1)=vi(j)
               k1=k1+1
            ENDIF
            vi(j)=0.d0
         ENDDO 
         ia1(i+1)=k1
      ENDDO
!
!.....backward solve
!
      ip=perm(n)
      x1(n,1)=x1(n,1)*diag_lu(n)
      x1(n,2)=x1(n,2)*diag_lu(n)
      x1(n,3)=x1(n,3)*diag_lu(n)
      x(ip,1)=x1(n,1)
      x(ip,2)=x1(n,2)
      x(ip,3)=x1(n,3)
      DO i=n-1,1,-1
         ip=perm(i)
         l=ia1(i+1)-ia1(i)
         jj=ia1(i)
         t1=x1(i,1)
         t2=x1(i,2)
         t3=x1(i,3)
         DO m=0,l-1
            j=ja1(jj+m)
            t1=t1-alu1(jj+m)*x1(j,1)
            t2=t2-alu1(jj+m)*x1(j,2)
            t3=t3-alu1(jj+m)*x1(j,3)
         ENDDO
         x1(i,1)=t1*diag_lu(i)
         x1(i,2)=t2*diag_lu(i)
         x1(i,3)=t3*diag_lu(i)
         x(ip,1)=x1(i,1)
         x(ip,2)=x1(i,2)
         x(ip,3)=x1(i,3)
      ENDDO
!
      END SUBROUTINE factor_solve0pv3
