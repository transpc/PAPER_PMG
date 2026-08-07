!
      SUBROUTINE pbcg_diag(eps,maxiter,ncell,neq,maxmt1,myrank,ia,ja,ju, &
                           au,arhsu,solu,np,                             &
                           izone)
!
!     Bi-conjugate gradient matrix solver with diagonal preconditioning
!
      USE Zinterface
      USE Zbicg        , ONLY: pbcgind
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!     input
      INTEGER neq
      INTEGER maxiter,ncell,maxmt1,myrank,np
      INTEGER ia(ncell+1),ja(maxmt1),ju(ncell)
      INTEGER izone
      REAL(8) eps
      REAL(8) au(maxmt1),arhsu(ncell)
!     output
      REAL(8) solu(neq)
!     local variable
      INTEGER i,n
      INTEGER its
      REAL(8) beta,alpha,alphad,omega,omegan,omegad
      REAL(8) rho,rhold
      REAL(8) ro,ro0
      REAL(8) eps1
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE::r0,z0,p0,rb,s0,y0,v0
!
!.....Array for CSR matrix
!
!
      n=neq
      ALLOCATE(r0(n),z0(n),p0(n),rb(n),s0(n),y0(n),v0(n))      
      r0(:)=0.0d0
      z0(:)=0.0d0
      p0(:)=0.0d0
      rb(:)=0.0d0
      s0(:)=0.0d0
      y0(:)=0.0d0
      v0(:)=0.0d0
!
!.....Predictor  By Diag. pc
!
      DO i=1,ncell
         solu(i)=arhsu(i)/au(ju(i))
      ENDDO
      its = 0
!
!.....communicating procedure
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
      CALL amux0lld(ncell,n,maxmt1,solu,rb,au,ja,ia)
!
      DO i=1,ncell
         arhsu(i)=arhsu(i)-rb(i)
      ENDDO
!
!.....compute initial residual vector-> arhsu-A*sol=arhsu
!
      ro0=0.d0
      DO i=1,ncell
         ro0=ro0+arhsu(i)*arhsu(i)
      ENDDO
!
      IF(np.gt.1) CALL allreducei_r1(ro0)
      ro0=DSQRT(ro0)
!
      IF(ro0.lt.1.d-20)THEN
         DEALLOCATE(r0,z0,p0,rb,s0,y0,v0)      
         RETURN
      ENDIF
!
      DO i=1,ncell
         r0(i)=arhsu(i)
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
!.....Diag PC y0=p0/Diag(*)
!
      DO i=1,ncell
         y0(i)=p0(i)/au(ju(i))
      ENDDO
!
      IF(np.gt.1) CALL communicate(y0,izone)
!      
      CALL amux0lld(ncell,n,maxmt1,y0,v0,au,ja,ia)
!
      alphad=0.0d0
      DO i=1,ncell
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
!......solve LUz0=s0
!      Diag PC z0=s0/Diag(*)
!
      DO i=1,ncell
         z0(i)=s0(i)/au(ju(i))
      ENDDO
!
      IF(np.gt.1) CALL communicate(z0,izone)
!
      CALL amux0lld(ncell,n,maxmt1,z0,r0,au,ja,ia)
!      
      omegan=0.0d0
      omegad=0.0d0
!
      DO i=1,ncell
         omegan=omegan+r0(i)*s0(i)
         omegad=omegad+r0(i)*r0(i)
      ENDDO     
!
      IF(np.gt.1)THEN
         CALL allreducei_r1(omegan)
         CALL allreducei_r1(omegad)
      ENDIF
!
      IF(omegad.eq.0.0d0) GOTO 990
      omega=omegan/omegad
!
      DO i=1,ncell
         solu(i)=solu(i)+alpha*y0(i)+omega*z0(i)
      ENDDO
!
      DO i=1,ncell
         r0(i)=s0(i)-omega*r0(i)
      ENDDO
!
      ro=0.d0
      DO i=1,ncell
         ro=ro+r0(i)*r0(i)
      ENDDO
!
      IF(np.gt.1) CALL allreducei_r1(ro)
      ro=DSQRT(ro)
!
      IF(ro/ro0.le.eps) GOTO 990
      IF(its.ge.maxiter) GOTO 991
      GOTO 10
!
  199 FORMAT('   iters =',i3,' norms=',1pe12.4,' ->',1pe12.4,'eps1=',1pe12.4)
!
  991 IF(myrank.eq.0)THEN
         WRITE(*,*)'          Iteration number for PBCG_DIAG exceeds ', maxiter
         WRITE(unit_log,*)'          Iteration number for PBCG_DIAG exceeds ', maxiter
      ENDIF
      solu(:)=0.0d0
      pbcgind=1  
!     
  990 CONTINUE
!
      IF(np.gt.1) CALL communicate(solu,izone)
!      
      DEALLOCATE(r0,z0,p0,rb,s0,y0,v0)      
!
      RETURN
      END SUBROUTINE pbcg_diag

!c***********************************************************************
!      SUBROUTINE amux0(ncell,n,maxmt1,x,y,a,ja,ia)
!!c-----------------------------------------------------------------------
!!c     Y = A * X
!!c     input:
!!c       n     = row DIMENSION of A
!!c       x     = array of length equal to the column DIMENSION of matrix A
!!c       a, ja, ia = input matrix in compressed sparse row format.
!!c     output:
!!c       y     = REAL array of length n, containing the product y=Ax
!!c-----------------------------------------------------------------------
!      INTEGER i, k
!      INTEGER ncell,n,maxmt1,ja(maxmt1),ia(ncell+1)
!      REAL(8) a(maxmt1),tmp
!      REAL(8) x(n),y(n)
!!c
!      DO i= 1,ncell
!        tmp = 0.d0
!        DO k=ia(i),ia(i+1)-1
!          tmp = tmp + a(k)*x(ja(k))
!        ENDDO
!        y(i) = tmp
!      ENDDO
!      return
!      END
