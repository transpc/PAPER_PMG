!
      SUBROUTINE bicgstab(eps,maxiter,miniter,svector,source,flag1)
!
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zbicg     , ONLY: pbcgind
      USE Zmpi      , ONLY: au,ju_a,ia_a,ja_a
!
      IMPLICIT NONE
!
      INTEGER iter,i,j,k,maxiter,miniter
      INTEGER flag1
      INTEGER cupid__ekd__,cupid__ekd__release,cupid__ekd__debug      
!
      LOGICAL,SAVE:: INITIAL
!
      REAL(8) poiss_diag(ncell_fluid)
      REAL(8) rms,eps
      REAL(8) resl,alf,bet,beto,gam,om,sum,vres,vv,res0,ureso
      REAL(8) svector(ncell_fp),source(ncell_fluid)
!
      COMMON/Zdist/cupid__ekd__,cupid__ekd__release,cupid__ekd__debug   
!
      DATA INITIAL /.TRUE./
!
      REAL(8),SAVE,ALLOCATABLE:: rvec(:),rtvec(:),pvec(:),uvec(:),vvec(:),zvec(:),dvec(:)
      IF(initial) THEN 
         ALLOCATE(rvec(ncell_fp),rtvec(ncell_fp),pvec(ncell_fp),uvec(ncell_fp),vvec(ncell_fp),zvec(ncell_fp),dvec(ncell_fp))
         initial=.false.
      ENDIF     
!
!.....Get diagonal element
!
      DO i=1,ncell_fluid
         k=ju_a(j)
         poiss_diag(i)=au(k)
      ENDDO
!
!.....Calculate initial residual vector and the norm
!
      res0=0.d0
!
      DO i=1,ncell_fluid
         sum=0.d0
         DO j=ia_a(i),ia_a(i+1)-1
            k=ja_a(j)
            sum=sum+au(j)*svector(k)
         ENDDO
         rvec(i)=source(i)-sum
         res0=res0+dabs(rvec(i))
      ENDDO
!
!.....Diagonal inverse
!
      DO i=1,ncell_fluid
         dvec(i)=1.d0/(poiss_diag(i)+1.0e-20)
      ENDDO
!
      rtvec=rvec
      pvec=0.d0
      uvec=0.d0
      vvec=0.d0
!
      alf=1.d0
      beto=1.d0
      gam=1.d0
!
!.....Inner iteration
!
      DO iter=1,maxiter
!
!........Calculate beta and omega
!
         bet=0.d0
!
         DO i=1,ncell_fluid
            bet=bet+rtvec(i)*rvec(i)
         ENDDO
!
         om=bet*gam/(alf*beto+1.0e-20)
         beto=bet
!
!........Calculate p
!
         DO i=1,ncell_fluid
            pvec(i)=rvec(i)+om*(pvec(i)-alf*uvec(i))
         ENDDO
!
         zvec=0.d0
!
!.......sgs preconditioner
!
         CALL pregs(pvec,zvec)
!
!........Calculate A*zvec -> u
!
         DO i=1,ncell_fluid
            sum=0.d0
            DO j=ia_a(i),ia_a(i+1)-1
               k=ja_a(j)
               sum=sum+au(j)*zvec(k)
            ENDDO
            uvec(i)=sum
         ENDDO
         ureso=0.d0
!
         DO i=1,ncell_fluid
            ureso=ureso+uvec(i)*rtvec(i)
         ENDDO
!
         gam=bet/(ureso+1.0e-20)
!
!........Calculate s and x
!
         DO i=1,ncell_fluid
            svector(i)=svector(i)+gam*zvec(i)
            rvec(i)=rvec(i)-gam*uvec(i)
         ENDDO
!
!........Solve M*zvec = rvec,
!
         zvec=0.d0
!
!........sgs preconditioner
!
         CALL pregs(rvec,zvec)
!
!........calculate A*zvec -> vvec
!
         DO i=1,ncell_fluid
            sum=0.d0
            DO j=ia_a(i),ia_a(i+1)-1
               k=ja_a(j)
               sum=sum+au(j)*zvec(k)
            ENDDO
            vvec(i)=sum
         ENDDO
!
!........Calculate alpha
!
         vres=0.d0
         vv=0.d0
         DO i=1,ncell_fluid
            vres=vres+vvec(i)*rvec(i)
            vv=vv+vvec(i)*vvec(i)
         ENDDO
!
         alf=vres/(vv+1.0e-20)
         resl=0.d0
!
         DO i=1,ncell_fluid
            svector(i)=svector(i)+alf*zvec(i)
            rvec(i)=rvec(i)-alf*vvec(i)
            resl=resl+rvec(i)*rvec(i)
         ENDDO
!
         rms=DSQRT((resl/ncell_fluid))
!
         IF((iter>miniter).and.(rms<eps))EXIT
      ENDDO
!
      IF(iter.le.maxiter)THEN
!        WRITE(*,*) '**iter = ', iter
      ELSE
         IF(cupid__ekd__.eq.0)WRITE(*,*)'          ** maximum iteration of BiCGSTAB **' 
         flag1=100
         svector(:)=0.d0
         pbcgind=1         
      ENDIF
!
      CONTAINS
!
!========================== SGS-Preconditioner =========================
!
!     GS matrix solver for preconditioning of the cgs solver
!
!=======================================================================
!
      SUBROUTINE pregs(bvector,xvector)
!
      IMPLICIT NONE
!
      INTEGER ::i,j,k
!
      REAL(8) bvector(*),xvector(*)
      REAL(8) ::sum
!
!.....forward sweep
!
      DO i=1,ncell_fluid
         sum=0.d0
         DO j=ia_a(i),ju_a(i)-1
            k=ja_a(j)
            sum=sum+au(j)*xvector(k)
         ENDDO
         DO j=ju_a(i)+1,ia_a(i+1)-1
            k=ja_a(j)
            sum=sum+au(j)*xvector(k)
         ENDDO
         xvector(i)=(bvector(i)-sum)*dvec(i)
      ENDDO
!
!.....backward sweep
!
      DO i=ncell_fluid,1,-1
         sum=0.d0
         DO j=ia_a(i),ju_a(i)-1
            k=ja_a(j)
            sum=sum+au(j)*xvector(k)
         ENDDO
         DO j=ju_a(i)+1,ia_a(i+1)-1
            k=ja_a(j)
            sum=sum+au(j)*xvector(k)
         ENDDO
         xvector(i)=(bvector(i)-sum)*dvec(i)
      ENDDO
      END SUBROUTINE pregs
!
      END SUBROUTINE bicgstab
