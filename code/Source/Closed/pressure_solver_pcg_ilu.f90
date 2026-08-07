!
      SUBROUTINE pcg_ilu(eps,maxiter,nintf,neq,maxmt1,myrank,                         &
                         maxmt_lu0,maxmt_lu1,                                         &
                         ia,ja,ju,au,arhsu,solu,np,diag_lu,alu0,alu1,ja0,ja1,ia0,ia1, &
                         izone)
!
!     Conjugate gradient matrix solver with diagonal preconditioning
!
      USE Zinterface
      USE Zio_unit     , ONLY: unit_log
      IMPLICIT NONE
!
!     input
      INTEGER neq
      INTEGER maxiter,maxmt1,nintf,myrank,np
      INTEGER maxmt_lu0,maxmt_lu1
      INTEGER ia(nintf+1),ju(nintf),ja(maxmt1)
      INTEGER ja0(maxmt_lu0),ja1(maxmt_lu1)
      INTEGER ia0(nintf+1),ia1(nintf+1)
      INTEGER izone
      REAL(8) eps
      REAL(8) au(maxmt1),arhsu(nintf)      
      REAL*8  diag_lu(nintf)
      REAL(8) alu0(maxmt_lu0),alu1(maxmt_lu1)
!     output
      REAL(8) solu(neq)
!     local variable
      INTEGER i,n
      INTEGER Allocatestatus,iter
      INTEGER tag
      REAL(8) bknum,bkden,akden,bk,ak
      REAL(8) err0,err
!     local arrays
      REAL(8) diag(nintf)
      REAL(8),DIMENSION(:),ALLOCATABLE::r,z,p0
!
      n=neq
      ALLOCATE(r(n),z(n),p0(n),stat=Allocatestatus)
      IF(AllocateStatus/=0) STOP "**Not enough memory "
!
      r(:)=0.0d0
      z(:)=0.0d0
      p0(:)=0.0d0
!
!.....Predictor  By Diag. pc
!
      r=0.d0 !! in fact, this is not necessary.
      DO i=1,nintf
         diag(i)=1.d0/au(ju(i))
         solu(i)=arhsu(i)*diag(i)
      ENDDO
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
      CALL amux0lld(nintf,n,maxmt1,solu,r,au,ja,ia)
!
      DO i=1,nintf
         r(i)=arhsu(i)-r(i)
      ENDDO
!
!.....ILU PC
!
      IF(nintf.gt.0) THEN
         CALL lusol0(nintf,n,maxmt_lu0,maxmt_lu1,       &
                     diag,                              &
                     diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                     r,z)
      ENDIF
!
      err0 = 0.d0
!
      DO i=1,nintf
         err0=err0+r(i)*r(i)
      ENDDO
!
      IF(np.gt.1) CALL allreducei_r1(err0)
      err0=DSQRT(err0)
!      
      IF(err0.eq.0.0d0) GOTO 502
!
      DO iter=1,maxiter
         tag = iter
!
!.......choi---calculmate coefficient bk and direction vector p
!
         bknum=0.d0
!
         DO i=1,nintf
            bknum=bknum+r(i)*z(i)
         ENDDO
!
         IF(np.gt.1) CALL allreducei_r1(bknum)
!
         IF(iter.eq.1)THEN
            p0=z
         ELSE
            IF(bkden.eq.0.0d0) GOTO 502
            bk=bknum/bkden
            p0=z+bk*p0
         ENDIF
         bkden = bknum
!
!.......choi    calculmate coefficient ak, new itermate x, new residual r
!
!.......Exchange p0 to compute A*p0
!
         IF(np.gt.1) CALL communicate(p0,izone)
!
         CALL amux0lld (nintf,n,maxmt1,p0,z,au,ja,ia) !!A*p0=z
!
         akden=0.d0
! 
         DO i=1,nintf
            akden=akden+p0(i)*z(i)
         ENDDO
!
         IF(np.gt.1) CALL allreducei_r1(akden)
!
         IF(akden.eq.0.0d0) GOTO 502
         ak=bknum/akden
!
         DO i=1,nintf
            solu(i)=solu(i)+ak*p0(i)
            r(i)=r(i)-ak*z(i)
         ENDDO
!
!.......ILU PC
!
         IF(nintf.gt.0) THEN
            CALL lusol0(nintf,n,maxmt_lu0,maxmt_lu1,       &
                        diag,                              &
                        diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                        r,z)
         ENDIF
!
         err=0.d0
!
         DO i=1,nintf
            err=err+r(i)*r(i)
         ENDDO
! 
         IF(np.gt.1) CALL allreducei_r1(err)
         err=DSQRT(err)
!
         IF(DABS(err).gt.1.d20)THEN
            IF(myrank.eq.0)THEN
               WRITE(*,*) 'Solution was not converged, residual=',err
               WRITE(unit_log,*) 'Solution was not converged, residual=',err
            ENDIF  
            STOP
         ENDIF
         IF(err/err0.le.eps.and.iter.ge.10) GOTO 502
!
         IF(iter.eq.maxiter) THEN
            IF (myrank.eq.0)THEN
               WRITE(*,*)'          Iteration number for PCG exceeds ', maxiter
               WRITE(unit_log,*)'          Iteration number for PCG exceeds ', maxiter
            ENDIF
            solu(:)=0.0d0
         ENDIF
!
      ENDDO ! Main -loop
!
  502 CONTINUE
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
      DEALLOCATE(r,z,p0)
!
      END SUBROUTINE pcg_ilu
