      SUBROUTINE solve_pbcg_mg(ierr)
!
!     Bi-conjugate gradient matrix solver with PMG preconditioning
!
      USE MD_MPI,       ONLY: nintf,myrank
      USE MD_geometry,  ONLY: nnode
      USE MD_matrix,    ONLY: nnz,ia,ja, au,u,b
      USE MD_parameter, ONLY: maxit, ndom
      USE MD_MG_index,  ONLY: crit_bcg_mg
	  
	  
!      
      IMPLICIT NONE
!
!.....Output: u
      integer :: ierr
!.....Local variable
      INTEGER :: i, n,np
      INTEGER :: its
      INTEGER :: maxiter
      REAL(8) :: beta,alpha,alphad,omega,omegan,omegad
      REAL(8) :: rho,rhold,rho1
      REAL(8) :: ro,ro0,ro01
      REAL(8) :: eps1,eps
	  REAL(8) :: alphad1,omegan1,omegad1,ro1
      REAL(8) :: scale_mg, scale_mg_inv

      REAL(8),DIMENSION(:),ALLOCATABLE:: r0,p0,rb,s0,v0,y0,z0,solu
!

      eps = crit_bcg_mg               ! this is criterion of solver
!	  crit = crit_1              ! criterion for MG pre.
	  maxiter = maxit
      scale_mg = 1.d1
	  
      n = nnode
	  np = ndom
!
      ALLOCATE(r0(nintf),p0(nintf),rb(nintf),s0(nintf),v0(nintf))
      ALLOCATE(y0(n),z0(n))
	  ALLOCATE(solu(n))
! C011-3r1: 첫 외부 반복에서 u = y0 로 읽히는 예조건 상태 벡터가 할당 직후의
! 힙 가비지였음 (C009 "예조건자 시드 실행 문맥 의존"의 근원 — 결과가 할당
! 순서·실행 환경의 힙 레이아웃에 결박). 0 초기화로 결정화, 2회차부터의
! warm start(직전 예조건 출력 재사용)는 그대로 유지
      y0 = 0.d0
      z0 = 0.d0
	  
      its = 0
!	  solu = u
    
!
!.....Predictor  By Diag. pc
!
!#include '../00_Module/c_Solver/avx.h'
!!DIR$ ASSUME_ALIGNED solu:avx,arhsu:avx,diag:avx

! for simple, just pre- by diagonal instead of ILU	 
! 
      IF(nintf.gt.0) THEN
	  
 !     DO i=1,nintf
 !        solu(i) = b(i)/au(ju(i))
 !     ENDDO
       solu = u
! for ILU smoothing: 
 !       alu = 0.d0            
 !     CALL ilupcp_new(nintf,nnode,nnz,ia,ja,ju,au,alu)
      
 !       alu = au           
 !     CALL ilupcp(nintf,nnode,nnz,ia,ja,ju,alu)
      
 !     CALL lusol0P(nintf,n,nnz,b,solu,alu,ja,ia,ju) ! z = M^-1*r
	  
      ENDIF
!
!.....communicating procedure
!
      IF(np.gt.1) CALL communicate_s(solu)
!
!.....compute initial residual vector-> arhsu-A*sol=arhsu------------
!
!      CALL amux0(ncell_pad,n,maxmt_pad,solu,rr,ap,jap,iap,iaa,ngroup,nbgroup)
	  
         CALL amux0P(nintf,n,nnz,solu,r0,au,ja,ia) !!A*solu=r0
!
      ro0=0.d0
!!DIR$ ASSUME_ALIGNED jaar:avx,arhsu:avx

      DO i=1,nintf
!         i1=jaar(i)
         r0(i) = b(i)-r0(i)
         ro0=ro0+r0(i)**2.d0
      ENDDO
!
!      IF(np.gt.1)CALL allreducei_r1(ro0)
!      ro0=SQRT(ro0)
	  
         IF(np.gt.1)THEN
            CALL allreduce_r_s(ro0,ro01)
            ro0 = ro01
         ENDIF
         ro0 = DSQRT(ro0)
		 
! 
      IF(ro0.lt.1.d-20) goto 990
!
      DO i=1,nintf
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
      DO i=1,nintf
         rho=rho+rb(i)*r0(i)
      ENDDO
!
         IF(np.gt.1)THEN
            CALL allreduce_r_s(rho,rho1)
            rho = rho1
         ENDIF  
!
! breakdown 가드 (E, LOOP F): rho/omega 소멸 시 0-나눗셈 -> NaN 침묵 정체 방지
      IF(DABS(rhold).LT.1.d-300 .OR. DABS(omega).LT.1.d-300) GOTO 991
!
      beta=rho/rhold*alpha/omega
!
      DO i=1,nintf
         p0(i)=r0(i)+beta*(p0(i)-omega*v0(i))
      ENDDO
!
!.....LUy0=p0
!
      IF(nintf.gt.0) THEN
	  
! test
!      CALL lusol0P(nintf,n,nnz,p0,y0,alu,ja,ia,ju) ! y0 = M^-1*p0
      
      if(scale_mg .lt. 1.d10) then
      scale_mg = scale_mg*10.d0
      scale_mg_inv = 1.d0/scale_mg
      endif
      
      
      u = y0
      b = scale_mg*p0
      
      CALL SOLVER_NEW(ierr)
!/
      
      y0 = scale_mg_inv*u
	  
      ENDIF

      IF(np.gt.1) CALL communicate_s(y0)
!
         CALL amux0P(nintf,n,nnz,y0,v0,au,ja,ia) !!A*y0=v0
		 
!      CALL amux0(ncell_pad,n,maxmt_pad,y0,rr,ap,jap,iap,iaa,ngroup,nbgroup)
!
      alphad=0.0d0
!!DIR$ ASSUME_ALIGNED jaar:avx
      DO i=1,nintf
!         i1=jaar(i)
!         v0(i)=rr(i1)
         alphad=alphad+rb(i)*v0(i)
      ENDDO
!     
!      IF(np.gt.1) CALL allreducei_r1(alphad)
	  
         IF(np.gt.1)THEN
            CALL allreduce_r_s(alphad,alphad1)
            alphad = alphad1
         ENDIF
!
      IF(DABS(alphad).LT.1.d-300) GOTO 991     ! breakdown 가드 (E)
      alpha=rho/alphad
!
      DO i=1,nintf
         s0(i)=r0(i)-alpha*v0(i)
      ENDDO
!
!.....solve LUz0=s0
!
      IF(nintf.gt.0) THEN

! test
!      CALL lusol0P(nintf,n,nnz,s0,z0,alu,ja,ia,ju) ! z0 = M^-1*s0
      
      
      u = z0
      b = scale_mg*s0
      
      CALL SOLVER_NEW(ierr)
!/
      
      z0 = scale_mg_inv*u
	  
      ENDIF
!
      IF(np.gt.1) CALL communicate_s(z0)
!
         CALL amux0P(nintf,n,nnz,z0,r0,au,ja,ia) !!A*z0=r0
		 
!      CALL amux0(ncell_pad,n,maxmt_pad,z0,rr,ap,jap,iap,iaa,ngroup,nbgroup)
       omegan=0.0d0
       omegad=0.0d0
!!DIR$ ASSUME_ALIGNED jaar:avx
      DO i=1,nintf
!         i1=jaar(i)
!         r0(i)=rr(i1)
         omegan=omegan+r0(i)*s0(i)
         omegad=omegad+r0(i)**2.d0
      ENDDO
!      
      IF(np.gt.1) THEN
!         omega2(1)=omegan
!         omega2(2)=omegad
!         CALL allreducei_r(omega2,2)
!         omegan=omega2(1)
!         omegad=omega2(2)

            CALL allreduce_r_s(omegan,omegan1)
            omegan = omegan1
			
            CALL allreduce_r_s(omegad,omegad1)
            omegad = omegad1
      ENDIF
!
!      IF(omegad.eq.0.0d0) GOTO 990
      omega = omegan/omegad
!
      ro=0.d0
!!DIR$ ASSUME_ALIGNED solu:avx
      DO i=1,nintf
         solu(i)=solu(i)+alpha*y0(i)+omega*z0(i)
         r0(i)  =s0(i)  -omega*r0(i)
         ro=ro+r0(i)**2.d0
      ENDDO
!
         IF(np.gt.1)THEN
            CALL allreduce_r_s(ro,ro1)
            ro = ro1
         ENDIF
         ro = DSQRT(ro)
		 
!      IF(np.gt.1) CALL allreducei_r1(ro)
!      ro=SQRT(ro)
!
      IF(ro.NE.ro) GOTO 991                       ! NaN 가드 (E)
      IF(ro/ro0.le.eps) GOTO 990
      IF(its.ge.maxiter) GOTO 991
      GOTO 10
!
!  199 FORMAT('   iters =',i3,' norms=',1pe12.4,' ->',1pe12.4,'eps1=',1pe12.4)
!
  991 IF (myrank.eq.0)THEN
         WRITE(*,*)'      Iteration number for PBCG_MG exceeds ', maxiter
         WRITE(999,*)'    Iteration number for PBCG_MG exceeds ', maxiter
      ENDIF
	  
	  ierr = 1
	DEALLOCATE(r0,p0,rb,s0,v0,y0,z0,solu)  
	RETURN
!      
  990 CONTINUE
  if(myrank == 0 ) then
   write(501,*)its,ro/ro0
  end if
!
      IF(np.gt.1) CALL communicate_s(solu)
!      
	  ierr = 0
	  u = solu      ! output is u
	  
	DEALLOCATE(r0,p0,rb,s0,v0,y0,z0,solu) 
! - - - - - - - - - - - - - - - - - - - 

    END SUBROUTINE solve_pbcg_mg
!

!-----------------------------------------------------------------------
!
! below are for subs. of communication: 
    
!------------------------------------------------------------------------------
!
      SUBROUTINE communicate_s(ss)
      
!DEC$IF defined (mpi_flag)
!     for fluid
      USE MD_MPI , ONLY: nnbd,nbdom,spt,rpt,sintf,rintf
      USE MD_GEOMETRY,  ONLY: nnode

!
      IMPLICIT NONE
!
      REAL(8) ss(*)
!
         CALL communicatez(nnbd,nbdom,spt,rpt,sintf,rintf,ss,nnode)

!DEC$ENDIF
!
      END SUBROUTINE communicate_s
!
!------------------------------------------------------------------------------
!
      SUBROUTINE communicatez(niut,iut,si,ri,sintf,rintf,ss,n)
      
!DEC$IF defined (mpi_flag)
!
!     This routine is a general "communicate" function for MPI
!
!
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!     
!     input
      INTEGER niut,n
      INTEGER iut(niut)
      INTEGER si(niut+1),ri(niut+1)
      INTEGER rintf(ri(niut+1)-1),sintf(si(niut+1)-1)      
!     output
      REAL(8) ss(n)
!     local variables
      INTEGER i
      INTEGER ierr,tag
!     local arays
      INTEGER status(MPI_STATUS_SIZE)
      INTEGER request(2*niut)
      REAL(8) rvar(ri(niut+1)-1),svar(si(niut+1)-1)
!
      DO i=1,si(niut+1)-1
         svar(i)=ss(sintf(i))
      ENDDO
!
      tag = 1
      DO i=1,niut
         CALL MPI_ISEND(svar(si(i)),si(i+1)-si(i), &
              MPI_DOUBLE_PRECISION,iut(i)-1,tag, &
              MPI_COMM_WORLD,request(i),ierr)
     
         CALL MPI_IRECV(rvar(ri(i)),ri(i+1)-ri(i), &
              MPI_DOUBLE_PRECISION,iut(i)-1,tag, &
              MPI_COMM_WORLD,request(i+niut),ierr)
      ENDDO
!
      DO i=1,niut
         CALL MPI_WAIT(request(i),status,ierr)
         CALL MPI_WAIT(request(i+niut),status,ierr)
      ENDDO
!
      DO i=1,ri(niut+1)-1     !n-ncell
         ss(rintf(i))=rvar(i)
      ENDDO
!
!DEC$ENDIF
!
      END SUBROUTINE communicatez
!
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!
      SUBROUTINE allreduce_r_s(a,a2)
!
!     This routine calls MPI_ALLREDUCE library for REAL single variable
!
!DEC$IF defined (mpi_flag)
      IMPLICIT NONE
!
      INCLUDE 'mpif.h'
!
      INTEGER ierr
!
      REAL(8) a,a1,a2
!
!   call mpi_barrier(mpi_comm_world,ierr)
   
      a1=a
      a2=0.0d0
      CALL MPI_ALLREDUCE(a1,a2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF
!
      RETURN
      END SUBROUTINE allreduce_r_s
!
!------------------------------------------------------------------------------
