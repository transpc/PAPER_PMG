      SUBROUTINE pbcg_ilu(eps,maxiter,ncell,ncell_pad,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                          diag,                                                      &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                         &
                          iap,jap,ap,jaar,iaa,ngroup,nbgroup,                        &
                          lev_typet,perm_r,                                          &
                          neq,arhsu,solu,izone,isPSolve)
!
!     Bi-conjugate gradient matrix solver with ILU preconditioning
!
      USE Zinterface
      USE Zcore    , ONLY: np,myrank
      USE Zbicg    , ONLY: pbcgind
      USE Zio_unit , ONLY: unit_log
!      
      IMPLICIT NONE
!
!.....Input
      INTEGER :: maxiter,ncell,ncell_pad,maxmt_pad,maxmt_lu0,maxmt_lu1
      INTEGER :: neq,lev_typet
      INTEGER :: izone
      INTEGER :: ia0(ncell+1),ia1(ncell+1),ja0(maxmt_lu0),ja1(maxmt_lu1)
      INTEGER :: iap(2,ngroup+1),jap(maxmt_pad),iaa(2,ngroup+1),ngroup,nbgroup(3,ngroup)
      INTEGER :: jaar(ncell)
      INTEGER :: perm_r(ncell)
      LOGICAL,OPTIONAL :: isPSolve
      REAL(8) :: eps
      REAL(8) :: diag(ncell)
      REAL*8  :: diag_lu(ncell),alu0(maxmt_lu0),alu1(maxmt_lu1)
      REAL(8) :: ap(maxmt_pad)
      REAL(8) :: arhsu(ncell)
!.....Output
      REAL(8) :: solu(neq)
!.....Local variable
      INTEGER :: i,i1,n
      INTEGER :: its
      REAL(8) :: beta,alpha,alphad,omega,omegan,omegad
      REAL(8) :: rho,rhold
      REAL(8) :: ro,ro0
      REAL(8) :: eps1
      REAL(8) :: omega2(2)
!.....Local arrays
      REAL(8) :: rr(ncell_pad)
      REAL(8) :: r0(ncell),p0(ncell),rb(ncell),s0(ncell),v0(ncell)
      REAL(8) :: y0(neq),z0(neq)
!
      n=neq
!
      its = 0
!
!.....Predictor  By Diag. pc
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED solu:avx,arhsu:avx,diag:avx
!!      DO i=1,ncell
!!         solu(i)=arhsu(i)*diag(i)
!!      ENDDO
      IF(ncell.gt.0) THEN
         IF(lev_typet.eq.0) THEN
            CALL lusol0(ncell,n,maxmt_lu0,maxmt_lu1,       &
                        diag,                              &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                        arhsu,solu)
         ELSE
            CALL lusol0r(ncell,n,maxmt_lu0,maxmt_lu1,       &
                         diag,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         arhsu,solu,perm_r)
         ENDIF
      ENDIF
!
!.....communicating procedure
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
!.....compute initial residual vector-> arhsu-A*sol=arhsu------------
!
      CALL amux0(ncell_pad,n,maxmt_pad,solu,rr,ap,jap,iap,iaa,ngroup,nbgroup)
!
      ro0=0.d0
!!DIR$ ASSUME_ALIGNED jaar:avx,arhsu:avx
      DO i=1,ncell
         i1=jaar(i)
         r0(i)=arhsu(i)-rr(i1)
         ro0=ro0+r0(i)**2
      ENDDO
!
      IF(np.gt.1)CALL allreducei_r1(ro0)
      ro0=SQRT(ro0)
! 
      IF(ro0.lt.1.d-20) RETURN
!
      DO i=1,ncell
         rb(i)=r0(i)
         p0(i)=0.0d0
         v0(i)=0.0d0
      ENDDO
!
!.....p0=r0? and v0=Ap0?
!
      eps1=eps
      rho=1.0d0
      alpha=1.0d0
      omega=1.0d0
!
   10 its=its+1
!      
      rhold=rho
!
      rho=0.0d0
      DO i=1,ncell
         rho=rho+rb(i)*r0(i)
      ENDDO
!
      IF(np.gt.1) CALL allreducei_r1(rho)   
!
      IF(rhold.eq.0.0d0.or.omega.eq.0.0d0) GOTO 990
!
      beta=rho/rhold*alpha/omega
!
      DO i=1,ncell
         p0(i)=r0(i)+beta*(p0(i)-omega*v0(i))
      ENDDO
!
!.....LUy0=p0
!
      IF(ncell.gt.0) THEN
         IF(lev_typet.eq.0) THEN
            CALL lusol0(ncell,n,maxmt_lu0,maxmt_lu1,       &
                        diag,                              &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                        p0,y0)
         ELSE
            CALL lusol0r(ncell,n,maxmt_lu0,maxmt_lu1,       &
                         diag,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         p0,y0,perm_r)
         ENDIF
      ENDIF

      IF(np.gt.1) CALL communicate(y0,izone)
!
      CALL amux0(ncell_pad,n,maxmt_pad,y0,rr,ap,jap,iap,iaa,ngroup,nbgroup)
!
      alphad=0.0d0
!!DIR$ ASSUME_ALIGNED jaar:avx
      DO i=1,ncell
         i1=jaar(i)
         v0(i)=rr(i1)
         alphad=alphad+rb(i)*v0(i)
      ENDDO
!     
      IF(np.gt.1) CALL allreducei_r1(alphad)
!
      IF(alphad.eq.0.0d0) GOTO 990
      alpha=rho/alphad
!
      DO i=1,ncell
         s0(i)=r0(i)-alpha*v0(i)
      ENDDO
!
!.....solve LUz0=s0
!
      IF(ncell.gt.0) THEN
         IF(lev_typet.eq.0) THEN
            CALL lusol0(ncell,n,maxmt_lu0,maxmt_lu1,       &
                        diag,                              &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                        s0,z0)
         ELSE
            CALL lusol0r(ncell,n,maxmt_lu0,maxmt_lu1,       &
                         diag,                              &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         s0,z0,perm_r)
         ENDIF
      ENDIF
!
      IF(np.gt.1) CALL communicate(z0,izone)
!
      CALL amux0(ncell_pad,n,maxmt_pad,z0,rr,ap,jap,iap,iaa,ngroup,nbgroup)
      omegan=0.0d0
      omegad=0.0d0
!!DIR$ ASSUME_ALIGNED jaar:avx
      DO i=1,ncell
         i1=jaar(i)
         r0(i)=rr(i1)
         omegan=omegan+r0(i)*s0(i)
         omegad=omegad+r0(i)**2
      ENDDO
!      
      IF(np.gt.1) THEN
         omega2(1)=omegan
         omega2(2)=omegad
         CALL allreducei_r(omega2,2)
         omegan=omega2(1)
         omegad=omega2(2)
      ENDIF
!
      IF(omegad.eq.0.0d0) GOTO 990
      omega=omegan/omegad
!
      ro=0.d0
!!DIR$ ASSUME_ALIGNED solu:avx
      DO i=1,ncell
         solu(i)=solu(i)+alpha*y0(i)+omega*z0(i)
         r0(i)  =s0(i)  -omega*r0(i)
         ro=ro+r0(i)**2
      ENDDO
!
      IF(np.gt.1) CALL allreducei_r1(ro)
      ro=SQRT(ro)
!
      IF(ro/ro0.le.eps) GOTO 990
      IF(its.ge.maxiter) GOTO 991
      GOTO 10
!
  199 FORMAT('   iters =',i3,' norms=',1pe12.4,' ->',1pe12.4,'eps1=',1pe12.4)
!
  991 IF (myrank.eq.0)THEN
         WRITE(*,*)'          Iteration number for PBCG_ILU exceeds ', maxiter
         WRITE(unit_log,*)'          Iteration number for PBCG_ILU exceeds ', maxiter
      ENDIF
!!DIR$ ASSUME_ALIGNED solu:avx
      DO i=1,n
         solu(i)=0.d0
      ENDDO
      pbcgind=1
      RETURN
!      
  990 CONTINUE
  if(myrank == 0 .and. PRESENT(isPSolve)) then
  if(isPSolve == .true.) then
   write(301,*)its,ro/ro0
  end if
  end if
!
      IF(np.gt.1) CALL communicate(solu,izone)
!      
      END SUBROUTINE pbcg_ilu
!
!-----------------------------------------------------------------------
!
      SUBROUTINE amux0lld(ncell,n,maxmt1,x,y,a,ja,ia)

      IMPLICIT NONE
!
!     Y = A * X
!     input:
!       n     = row DIMENSION of A
!       x     = array of length equal to the column DIMENSION of matrix A
!       a, ja, ia = input matrix in compressed sparse row format.
!     output:
!       y     = REAL array of length n, containing the product y=Ax
!
!     input
      INTEGER ncell,n,maxmt1,ja(maxmt1),ia(ncell+1)
      REAL(8) a(maxmt1),x(n)
!     output
      REAL(8) y(n)
!     local variables
      INTEGER i,k,l,m
      INTEGER j0,j1,j2,j3,j4,j5
      REAL(8) tmp
      REAL(8) t0,t1,t2,t3,t4,t5
!
      DO i= 1,ncell
         l=ia(i+1)-ia(i)
         tmp=0.d0
!
         IF    (l.eq.1) THEN
            k=ia(i)
            j0=ja(k  )
            t0=a(k  )*x(j0)
            tmp=tmp+t0
         ELSEIF(l.eq.2) THEN
            k=ia(i)
            j0=ja(k  )
            j1=ja(k+1)
            t0=a(k  )*x(j0)
            t1=a(k+1)*x(j1)
            tmp=tmp+t0+t1
         ELSEIF(l.eq.3) THEN
            k=ia(i)
            j0=ja(k  )
            j1=ja(k+1)
            j2=ja(k+2)
            t0=a(k  )*x(j0)
            t1=a(k+1)*x(j1)
            t2=a(k+2)*x(j2)
            tmp=tmp+t0+t1+t2
         ELSEIF(l.eq.4) THEN
            k=ia(i)
            j0=ja(k  )
            j1=ja(k+1)
            j2=ja(k+2)
            j3=ja(k+3)
            t0=a(k  )*x(j0)
            t1=a(k+1)*x(j1)
            t2=a(k+2)*x(j2)
            t3=a(k+3)*x(j3)
            tmp=tmp+t0+t1+t2+t3
         ELSEIF(l.eq.5) THEN
            k=ia(i)
            j0=ja(k  )
            j1=ja(k+1)
            j2=ja(k+2)
            j3=ja(k+3)
            j4=ja(k+4)
            t0=a(k  )*x(j0)
            t1=a(k+1)*x(j1)
            t2=a(k+2)*x(j2)
            t3=a(k+3)*x(j3)
            t4=a(k+4)*x(j4)
            tmp=tmp+t0+t1+t2+t3+t4
         ELSEIF(l.eq.6) THEN
            k=ia(i)
            j0=ja(k  )
            j1=ja(k+1)
            j2=ja(k+2)
            j3=ja(k+3)
            j4=ja(k+4)
            j5=ja(k+5)
            t0=a(k  )*x(j0)
            t1=a(k+1)*x(j1)
            t2=a(k+2)*x(j2)
            t3=a(k+3)*x(j3)
            t4=a(k+4)*x(j4)
            t5=a(k+5)*x(j5)
            tmp=tmp+t0+t1+t2+t3+t4+t5
         ELSE
           DO m=0,l-1
              k=ia(i  )+m
              j0=ja(k  )
              t0=a(k  )*x(j0)
              tmp=tmp+t0
           ENDDO
        ENDIF
        y(i)=tmp
        ENDDO
!
      RETURN
      END SUBROUTINE amux0lld
!
!-----------------------------------------------------------------------
!
      SUBROUTINE ilupc(ncell,maxmt1,ia,ja,ju,au,alu,iend)
!
      IMPLICIT NONE
!
!     input
      INTEGER ncell,maxmt1
      INTEGER ia(ncell+1),ja(maxmt1),ju(ncell),iend(ncell)
      REAL(8) au(maxmt1)
!     output
      REAL(8) alu(maxmt1)
!     local variables
      INTEGER i,j,j3
      INTEGER ip1,ip2
      REAL(8) tl
!
!.....ilu operation
!
      DO i=1,ncell
         DO j=ia(i),iEND(i)
            alu(j)=au(j)
         ENDDO
         DO j=ia(i),ju(i)-1
            j3=ja(j)
            ip1=ju(j3)
            tl=au(j)*alu(ip1)
            alu(j)=tl
            ip1=ip1+1
            ip2=j+1
110         CONTINUE
            IF(ip1.gt.iEND(j3)) GOTO 100
            IF(ip2.gt.iEND(i)) GOTO 100
            IF(ja(ip1).eq.ja(ip2)) THEN
               alu(ip2)=alu(ip2)-tl*alu(ip1)
               ip1=ip1+1
               ip2=ip2+1
               GOTO 110 
            ELSEIF(ja(ip1).gt.ja(ip2)) THEN
               ip2=ip2+1
               GOTO 110 
            ELSE
               ip1=ip1+1
               GOTO 110 
            ENDIF
100         CONTINUE
         ENDDO
         j=ju(i)
         alu(j)=1.0/alu(j)
      ENDDO   
!
      RETURN
      END SUBROUTINE ilupc
!
      SUBROUTINE ilupc1(ncell,maxmt1,maxmt_lu0,maxmt_lu1,ia,ja,ju,au,iend, &
                        diag_lu,alu0,alu1,ia0,ia1)
!
      IMPLICIT NONE
!
!     input
      INTEGER ncell,maxmt1,maxmt_lu0,maxmt_lu1
      INTEGER ia(ncell+1),ju(ncell),iend(ncell)
      INTEGER ja(maxmt1)
      INTEGER ia0(ncell+1),ia1(ncell+1)
      REAL(8) au(maxmt1)
!     output
      REAL(8) diag_lu(ncell)
      REAL(8) alu0(maxmt_lu0),alu1(maxmt_lu1)
!     local variables
      INTEGER i,j,i1,i2,i3
      INTEGER :: jj,jl
      INTEGER ip1,ip2,ip11
      REAL(8) tl
!
!.....ilu operation
!
!
      DO i=1,ncell
         i2=ncell-i+1
         jl=ia0(i)
         DO j=ia(i),ju(i)-1
            alu0(jl)=au(j)
            jl=jl+1
         ENDDO
            j=ju(i)
            diag_lu(i)=au(j)
         jj=ia1(i2+1)-1
         DO j=ju(i)+1,iend(i)
            alu1(jj)=au(j)
            jj=jj-1
         ENDDO
!!
         DO j=ia(i),ju(i)-1
            jl=ia0(i)+j-ia(i)
            i1=ja(j)
            i3=ncell-i1+1
            ip1=ju(i1)
            tl=alu0(jl)*diag_lu(i1)
            alu0(jl)=tl
            jl=jl+1
!!!!!!!!!!!!!!!!!
            ip11=ia1(i3+1)-1
            ip1=ip1+1
            ip2=j+1
110         CONTINUE
            IF(ip1.gt.iEND(i1)) GOTO 100
            IF(ip2.ge.ju(i)) GOTO 100
            IF(ja(ip1).eq.ja(ip2)) THEN
               alu0(jl)=alu0(jl)-tl*alu1(ip11)
            jl=jl+1
               ip11=ip11-1
               ip1=ip1+1
               ip2=ip2+1
               GOTO 110 
            ELSEIF(ja(ip1).gt.ja(ip2)) THEN
            jl=jl+1
               ip2=ip2+1
               GOTO 110 
            ELSE
               ip11=ip11-1
               ip1=ip1+1
               GOTO 110 
            ENDIF
100         CONTINUE
!!!!!!!!!!!!!!!!!
! diag
310         CONTINUE
            IF(ip1.gt.iEND(i1)) GOTO 300
            IF(ip2.gt.ju(i)) GOTO 300
            IF(ja(ip1).lt.ja(ip2)) THEN
               ip11=ip11-1
               ip1=ip1+1
               GOTO 310 
            ELSEIF(ja(ip1).eq.ja(ip2)) THEN
               diag_lu(i)=diag_lu(i)-tl*alu1(ip11)
               ip11=ip11-1
               ip1=ip1+1
               ip2=ip2+1
            ELSEIF(ja(ip1).gt.ja(ip2)) THEN
               ip2=ip2+1
            ENDIF
300         CONTINUE
!!!!!!!!!!!!!!!!!
         jj=ia1(i2+1)-1
210         CONTINUE
            IF(ip1.gt.iEND(i1)) GOTO 200
            IF(ip2.gt.iEND(i)) GOTO 200
            IF(ja(ip1).eq.ja(ip2)) THEN
               alu1(jj)=alu1(jj)-tl*alu1(ip11)
               jj=jj-1
               ip11=ip11-1
               ip1=ip1+1
               ip2=ip2+1
               GOTO 210 
            ELSEIF(ja(ip1).gt.ja(ip2)) THEN
               jj=jj-1
               ip2=ip2+1
               GOTO 210 
            ELSE
               ip11=ip11-1
               ip1=ip1+1
               GOTO 210 
            ENDIF
200         CONTINUE
         ENDDO
         diag_lu(i)=1.d0/diag_lu(i)
      ENDDO
!
!
      RETURN
      END SUBROUTINE ilupc1
!
      SUBROUTINE gener_vect_size(n,max_neigh,ia,ngroup)
!
      IMPLICIT NONE
!
      INTEGER :: n,max_neigh
      INTEGER :: ngroup
      INTEGER :: ia(n+1)
!
      INTEGER :: i,nn
      INTEGER :: ic(max_neigh)
!
      do i=1,max_neigh
         ic(i)=0
      enddo
      DO i=1,n
         nn=ia(i+1)-ia(i)
         ic(nn)=ic(nn)+1
      ENDDO
      ngroup=0
      DO i=1,max_neigh
         if(ic(i).ne.0) then
          ngroup=ngroup+1
         endif
      ENDDO
!
      END SUBROUTINE gener_vect_size
!
      SUBROUTINE gener_vect_u(n,nnz_pad,n_pad,n_padv,max_neigh, &
                              ia,iaa,iap,                &
                              jaa,jaar,ngroup,nbgroup)
!
      USE Zvec_param     , ONLY: vl_f
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: n,max_neigh
      INTEGER :: ngroup
      INTEGER :: ia(n+1)
!.....Output
      INTEGER :: nnz_pad,n_pad,n_padv
      INTEGER :: iaa(2,ngroup+1),iap(2,ngroup+1)
      INTEGER :: jaa(n),jaar(n)
      INTEGER :: nbgroup(3,ngroup)
!.....Local variables
      INTEGER :: i,ip,nn,nn1,nnp
!.....Local arrays
      INTEGER :: indx(max_neigh),ic(max_neigh)
!
      do i=1,max_neigh
         ic(i)=0
      enddo
      DO i=1,n
         nn=ia(i+1)-ia(i)
         ic(nn)=ic(nn)+1
      ENDDO
      ip=1
      iaa(1,1)=1
      iaa(2,1)=1
      iap(1,1)=1
      iap(2,1)=1
      DO i=1,max_neigh
        if(ic(i).ne.0) then
         nn=ic(i)
         nn1=i
         nnp=((nn-1)/vl_f+1)*vl_f
         nbgroup(1,ip)=nn
         nbgroup(2,ip)=nn1
         nbgroup(3,ip)=nnp
         iaa(1,ip+1)=iaa(1,ip)+nn
         iaa(2,ip+1)=iaa(2,ip)+nnp
         iap(1,ip+1)=iap(1,ip)+nn*nn1
         iap(2,ip+1)=iap(2,ip)+nnp*nn1
         ip=ip+1
        endif
      ENDDO
      nnz_pad=iap(2,ngroup+1)-1
      n_padv=iaa(2,ngroup+1)-1
      IF(n.eq.0) THEN
         n_pad=0
      ELSE
         n_pad=((n-1)/vl_f+1)*vl_f
      ENDIF
!
      DO ip=1,ngroup
         nn1=nbgroup(2,ip)
         indx(nn1)=iaa(1,ip)
         ic(nn1)=0
      ENDDO
      DO i=1,n
         nn=ia(i+1)-ia(i)
         ip=indx(nn)
         jaa(ip+ic(nn))=i
         ic(nn)=ic(nn)+1
      ENDDO
!
      DO ip=1,ngroup
         nn=nbgroup(2,ip)
         indx(nn)=iaa(2,ip)
         ic(nn)=0
      ENDDO
      DO i=1,n
         nn=ia(i+1)-ia(i)
         ip=indx(nn)
         jaar(i)=ip+ic(nn)
         ic(nn)=ic(nn)+1
      ENDDO
!
      END SUBROUTINE gener_vect_u
!
      SUBROUTINE copy_ja_vector(n,maxmt1,maxmt_pad,ja,ia,jap,iap, &
                                jaa,iaa,ngroup,nbgroup)
!
      IMPLICIT NONE
!
!     input
      INTEGER :: n,maxmt1,maxmt_pad
      INTEGER :: ja(maxmt1),ia(n+1)
      INTEGER :: ngroup,nbgroup(3,ngroup)
      INTEGER :: iap(2,ngroup+1)
      INTEGER :: jaa(n),iaa(2,ngroup+1)
!     ouput
      INTEGER :: jap(maxmt_pad)
!
      INTEGER :: i,j,k,j0,ip,ic1,ic2
      INTEGER :: nn,nn1,nnp
!
      DO i=1,ngroup
           nn =nbgroup(1,i)
           nn1=nbgroup(2,i)
           nnp=nbgroup(3,i)
           ic1=iap(2,i)-1
           j0=iaa(1,i)-1
!DIR$ SIMD
           do k=1,nn
              ip=jaa(j0+k)
              ic2=ic1+k
              do j=ia(ip),ia(ip)+nn1-1
                 jap(ic2)=ja(j)
                 ic2=ic2+nnp
              enddo
           enddo
      ENDDO
!
      END SUBROUTINE copy_ja_vector
!
      SUBROUTINE copy_a_vector(ncell,maxmt1,maxmt_pad,        &
                               ia,au,                         & 
                               iap,ap,jaa,iaa,ngroup,nbgroup)
!
      IMPLICIT NONE
!
!     input
      INTEGER :: ncell,maxmt1,maxmt_pad
      INTEGER :: ngroup,nbgroup(3,ngroup)
      INTEGER :: ia(ncell+1)
      INTEGER :: iap(2,ngroup+1)
      INTEGER :: jaa(ncell),iaa(2,ngroup+1)
      REAL(8) :: au(maxmt1)
!     ouput
      REAL(8) :: ap(maxmt_pad)
!
      INTEGER i,j,k,nn,nn1,ip,ic1,ic2
      INTEGER nnp,j0
! 
      DO i=1,ngroup
           nn =nbgroup(1,i)
           nn1=nbgroup(2,i)
           nnp=nbgroup(3,i)
           ic1=iap(2,i)-1
           j0=iaa(1,i)-1
!DIR$ SIMD
           do k=1,nn
              ip=jaa(j0+k)
              ic2=ic1+k
              do j=ia(ip),ia(ip)+nn1-1
                 ap(ic2)=au(j)
                 ic2=ic2+nnp
              enddo
           enddo
      ENDDO
!
      END SUBROUTINE copy_a_vector
