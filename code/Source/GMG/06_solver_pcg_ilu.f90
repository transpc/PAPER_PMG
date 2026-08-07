!
      SUBROUTINE pcg_ilu_s(np,myrank,crit,maxiter,nintf,neq,nnz,      &
                       ia,ja,ju,au,alu,arhsu,solu)
!
!

      USE Zinterface
!      use MD_MG_matrix
      use MD_geometry
      
      IMPLICIT NONE
!
!     input
      INTEGER neq,np,myrank,maxiter,nnz,nintf
      INTEGER ia(nintf+1),ju(nintf),ja(nnz)
      REAL(8) crit
      REAL(8) au(nnz),alu(nnz),arhsu(neq)      
!     output
      REAL(8) solu(neq)
!     local variable
      INTEGER i,n
      INTEGER Allocatestatus,iter
      INTEGER :: izone=0
      REAL(8) bknum,bkden,akden,bk,ak,bknum1,akden1
      REAL(8) err0,err1,err,eps
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE::r,z,p0
!
      n=neq
      ALLOCATE(r(n),z(n),p0(n),stat=Allocatestatus)
      IF(AllocateStatus/=0) STOP "**Not enough memory "


	  eps = 1.d-16
!
!      r(:)=0.0d0
!      z(:)=0.0d0
!      p0(:)=0.0d0
!
!.....Predictor  By Diag. pc
!
!      r=0.d0 !! in fact, this is not necessary.
!      DO i=1,nintf
!         solu(i)=arhsu(i)/au(ju(i))
!      ENDDO
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
      CALL amux0P(nintf,n,nnz,solu,r,au,ja,ia)
!
       r(1:nintf)=arhsu(1:nintf)-r(1:nintf)
	   
!      DO i=1,nintf
!         r(i)=arhsu(i)-r(i)
!      ENDDO
!
!.....ILU PC
!
      IF(nintf.gt.0) CALL lusol0P(nintf,n,nnz,r,z,alu,ja,ia,ju)
!
      err0 = 0.d0
!
      DO i=1,nintf
         err0=err0+r(i)*r(i)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreducei_r1(err0)
      ENDIF
      err0=DSQRT(err0)
!      
      IF(err0.LE.eps) GOTO 502
	  
      p0=0.d0
	  bkden = 1.d0
!
      DO iter=1,maxiter
!
!calculmate coefficient bk and direction vector p
!
         bknum=0.d0
!
         DO i=1,nintf
            bknum=bknum+r(i)*z(i)
         ENDDO
!
         IF(np.gt.1)THEN
            CALL allreducei_r1(bknum)
         ENDIF
!
            bk=bknum/bkden
            p0=z+bk*p0

            bkden = bknum
!
!calculmate coefficient ak, new itermate x, new residual r
!
!Exchange p0 to compute A*p0
!
         IF(np.gt.1) CALL communicate(p0,izone)
!
         CALL amux0P(nintf,n,nnz,p0,z,au,ja,ia) !!A*p0=z
!
         akden=0.d0
! 
         DO i=1,nintf
            akden=akden+p0(i)*z(i)
         ENDDO
!
         IF(np.gt.1)THEN
            CALL allreducei_r1(akden)
         ENDIF
!
         ak=bknum/akden
!
         DO i=1,nintf
            solu(i)=solu(i)+ak*p0(i)
            r(i)=r(i)-ak*z(i)
         ENDDO
!
!.......ILU PC
!
          CALL lusol0P(nintf,n,nnz,r,z,alu,ja,ia,ju)
!
         err=0.d0
!
         DO i=1,nintf
            err=err+r(i)*r(i)
         ENDDO
! 
         IF(np.gt.1)THEN
            CALL allreducei_r1(err)
         ENDIF
         err=DSQRT(err)
!
         IF(err.gt.1.d20)THEN
            IF(myrank.eq.0)THEN
               WRITE(*,*) 'blow-up in the solver, residual=',err
 !              WRITE(97,*) 'Solution was not converged, residual=',err
            ENDIF  
            STOP
         ENDIF
         IF(err/err0.le.crit) GOTO 502
!
!
      ENDDO ! Main -loop
	  
      WRITE(*,*) 'Solution was not converged',iter,err	
      pause
      
	  
!
  502 CONTINUE
!
      IF(np.gt.1) CALL communicate(solu,izone)
      IF(myrank.EQ.0) THEN
      WRITE(*,*) 'Solution is converged',iter,err/err0	
      ENDIF

!
      DEALLOCATE(r,z,p0)
!
      RETURN
      END SUBROUTINE pcg_ilu_s
	  
!***********************************************************************
      SUBROUTINE amux0P(nintf,n,nnz,x,y,a,ja,ia)
!-----------------------------------------------------------
!     Y = A * X
!     input:
!       n     = row dimension of A
!       x     = array of length equal to the column dimension of matrix A
!       a, ja, ia = input matrix in compressed sparse row format.
!     output:
!       y     = real array of length n, containing the product y=Ax
!-------------------------------------------------------------------
      IMPLICIT NONE
!      
      INTEGER  nintf,n,nnz
	  INTEGER  ja(nnz),ia(nintf+1)
      REAL*8 a(nnz),x(n)
! 
      REAL*8 y(n)
! tmp
      INTEGER i, k
      REAL*8 tmp
	  
!
      DO i= 1,nintf
        tmp = 0.d0
        DO k=ia(i),ia(i+1)-1
          tmp = tmp + a(k)*x(ja(k))
        ENDDO
        y(i) = tmp
      ENDDO
!	  
      RETURN
      END
! - - - - - - - - - 
!=====================================================================!
subroutine lusol0P(nintf,n,nnz,y,x,alu,ja,ia,ju)
implicit none
!
!---------------------------------------------------------------------!
!                                                                     !
!  Performing a forward followed by a backward solve                  !
!  for LU matrix as produced by  ILUT                                 !
!  input:                                                             !
!    n     = row dimension of A                                       !
!    y     = the right hand side of the linear system                 !
!    alu, ja, ia = input matrix in compressed sparse row format.      !
!  output:                                                            !
!    x     = the solution                                             !
!                                                                     !
!---------------------------------------------------------------------!
!
INTEGER :: n,nintf,nnz
INTEGER :: ja(nnz)
INTEGER :: ju(nintf)
INTEGER :: ia(nintf+1)
REAL*8  :: alu(nnz)
REAL*8  :: x(n)
REAL*8  :: y(n)  
! temp
INTEGER :: i,k1,k2
REAL*8  :: tmp
!
!--forward solver
!
x(1) = y(1)
do i=2,nintf
   k1 = ia(i)
   k2 = ju(i)-1
   tmp = y(i) - dot_product(alu(k1:k2),x(ja(k1:k2)))
   x(i) = tmp
enddo
!
!--backward solve
!
do i=nintf,1,-1
   k1 = ju(i)+1
   k2 = ia(i+1)-1
   tmp = x(i) - dot_product(alu(k1:k2),x(ja(k1:k2)))
   x(i) = tmp*alu(ju(i))
enddo
!=====
!=====
return
end subroutine

!-----------------------------------------------------
!
      SUBROUTINE ilupcp(nintf,nnode,nnz,ia,ja,ju,alu)
!
      IMPLICIT NONE
!
!     input
      INTEGER nintf,nnode,nnz
      INTEGER ia(nintf+1),ja(nnz),ju(nintf)
!     output
      REAL(8) alu(nnz)
!     local variables
      INTEGER i,j,j1,j2,j3
      INTEGER ip1,ip2
      REAL(8) tl
!
!.....ilu operation
!
      DO i=1,nintf
         j1=ia(i)
         j2=ia(i+1)-1
         DO j=j1,ju(i)-1
            j3=ja(j)
            ip1=ju(j3)
            tl=alu(j)*alu(ip1)
            alu(j)=tl
            ip1=ip1+1
            ip2=j+1
110         CONTINUE
            IF(ip1.gt.ia(j3+1)-1) GOTO 100
            IF(ip2.gt.j2) GOTO 100
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
    END SUBROUTINE ilupcp
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
					
      SUBROUTINE pbicg_ilu_s(np,myrank,eps,maxiter,nintf,neq,nnz,      &
                       ia,ja,ju,au,alu,arhsu0,solu,ierror)
!
!     Bi-conjugate gradient matrix solver with ILU preconditioning
!
      USE Zinterface
      use MD_geometry
!      
      IMPLICIT NONE
!
!     input
      INTEGER neq,nnz,ierror
      INTEGER maxiter,nintf,myrank,np
      INTEGER ia(nintf+1),ja(nnz),ju(nintf)
      REAL(8) eps
      REAL(8) au(nnz),alu(nnz),arhsu0(neq)
!     output
      REAL(8) solu(neq)
!     local variable
      INTEGER i,n,j,k
      INTEGER its
      INTEGER :: izone=0
      REAL(8) beta,alpha,alphad,omega,omegan,omegad
      REAL(8) rho,rhold,rho1,alphad1,omegan1,omegad1
      REAL(8) ro,r1,ro0,ro1
      REAL(8) t0
      REAL(8) eps1
!     local arrays
      REAL(8),DIMENSION(:),ALLOCATABLE::r0,z0,p0,rb,s0,y0,v0,arhsu

!
      n=neq
      ALLOCATE(r0(n),z0(n),p0(n),rb(n),s0(n),y0(n),v0(n),arhsu(n))      
      DO i=1,n
         r0(i)=0.0d0
         z0(i)=0.0d0
         p0(i)=0.0d0
         rb(i)=0.0d0
         s0(i)=0.0d0
         y0(i)=0.0d0
         v0(i)=0.0d0
         arhsu(i) = arhsu0(i)
      ENDDO
!
!.....Predictor  By Diag. pc
!
      DO i=1,nintf
         solu(i)=arhsu(i)/au(ju(i))
      ENDDO
!
      its = 0
!
!.....communicating procedure
!
      IF(np.gt.1) CALL communicate(solu,izone)
!
	   CALL amux0P(nintf,n,nnz,solu,rb,au,ja,ia) !!A*p0=z
!
      DO i=1,nintf
         arhsu(i)=arhsu(i)-rb(i)
      ENDDO
!
!.....compute initial residual vector-> arhsu-A*sol=arhsu------------

      ro0=0.d0
      DO i=1,nintf
         ro0=ro0+arhsu(i)*arhsu(i)
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreducei_r1(ro0)
      ENDIF
      ro0=DSQRT(ro0)
! 
      IF(ro0.lt.1.d-20)THEN
         DEALLOCATE(r0,z0,p0,rb,s0,y0,v0)      
         RETURN
      ENDIF
!
      DO i=1,nintf
         r0(i)=arhsu(i)
         rb(i)=r0(i)
         p0(i)=0.0d0
         v0(i)=0.0d0
      ENDDO
!
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
         CALL allreducei_r1(rho)	   
      ENDIF
!
      IF(rhold.eq.0.0d0.or.omega.eq.0.0d0) GOTO 990
!
      beta=rho/rhold*alpha/omega
!
      DO i=1,nintf
         p0(i  )=r0(i  )+beta*(p0(i  )-omega*v0(i  ))
      ENDDO
!
!.....LUy0=p0
!
          CALL lusol0P(nintf,n,nnz,p0,y0,alu,ja,ia,ju)

      IF(np.gt.1) CALL communicate(y0,izone)
	   CALL amux0P(nintf,n,nnz,y0,v0,au,ja,ia) !!A*p0=z
!
!
      alphad=0.0d0
      DO i=1,nintf
         alphad=alphad+rb(i)*v0(i)
      ENDDO      
!     
      IF(np.gt.1)THEN
         CALL allreducei_r1(alphad)
      ENDIF
!
      IF(alphad.eq.0.0d0) GOTO 990
      alpha=rho/alphad
!
      DO i=1,nintf
         s0(i)=r0(i)-alpha*v0(i)
      ENDDO
!
!.....solve LUz0=s0
!
	  
	  CALL lusol0P(nintf,n,nnz,s0,z0,alu,ja,ia,ju)
!
      IF(np.gt.1) CALL communicate(z0,izone)
	  
	   CALL amux0P(nintf,n,nnz,z0,r0,au,ja,ia) !!A*p0=z
!
!
      omegan=0.0d0
      omegad=0.0d0
      DO i=1,nintf
         omegan=omegan+r0(i)*s0(i)
         omegad=omegad+r0(i)*r0(i)
      ENDDO
!      
      IF(np.gt.1) THEN
         CALL allreducei_r1(omegan)
         CALL allreducei_r1(omegad)     
      ENDIF
!
      IF(omegad.eq.0.0d0) GOTO 990
      omega=omegan/omegad
!
      ro=0.d0
      DO i=1,nintf
         solu(i  )=solu(i  )+alpha*y0(  i)+omega*z0(i  )
         r0(i  )  =s0(i  )-omega*r0(i  )
         t0=r0(i  )**2
         ro=ro+t0
      ENDDO
!
      IF(np.gt.1)THEN
         CALL allreducei_r1(ro)
         ro=DSQRT(ro)
      ELSE
         ro=DSQRT(ro)
      ENDIF
!
      IF(ro/ro0.le.eps) GOTO 990
      IF(its.ge.maxiter) GOTO 991
      GOTO 10
!
!
  991 IF (myrank.eq.0)THEN
         WRITE(*,*)'PBCG_ILU exceeds ', maxiter,ro/ro0
 !  
      ENDIF
      solu(:)=0.0d0
      ierror = 1
!      pbcgind=1
!      
  990 CONTINUE
!
      IF(np.gt.1) CALL communicate(solu,izone)
      
  if(myrank == 0 ) then
   write(601,*)its,ro/ro0
  end if
!      
      DEALLOCATE(r0,z0,p0,rb,s0,y0,v0,arhsu)      
!
      RETURN
      END
!
!-----------------------------------------------------------------------
!-----------------------------------------------------
!
      SUBROUTINE ilupcp_new(nintf,nnode,nnz,ia,ja,ju,au,alu)
!
      IMPLICIT NONE
!
!     input
      INTEGER nintf,nnode,nnz
      INTEGER ia(nintf+1),ja(nnz),ju(nintf)
      REAL(8) au(nnz)
!     output
      REAL(8) alu(nnz)
!     local variables
      INTEGER i,j,j1,j2,j3
      INTEGER ip1,ip2
      REAL(8) tl
!
!.....ilu operation
!
      DO i=1,nintf
         j1=ia(i)
         j2=ia(i+1)-1
         DO j=j1,ju(i)-1
            j3=ja(j)
            ip1=ju(j3)
            tl=au(j)*alu(ip1)
            alu(j)=tl
            ip1=ip1+1
            ip2=j+1
110         CONTINUE
            IF(ip1.gt.ia(j3+1)-1) GOTO 100
            IF(ip2.gt.j2) GOTO 100
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
    END SUBROUTINE ilupcp_new
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

!-----------------------------------------------------------------------

SUBROUTINE ilu0_gpt(nintf,n, nnz, ia, ja, ju, alu)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nintf,n       ! Number of rows/columns in the matrix
    INTEGER, INTENT(IN) :: nnz     ! Number of non-zero elements
    INTEGER, INTENT(IN) :: ia(nintf+1) ! Row pointer array for CSR format
    INTEGER, INTENT(IN) :: ja(nnz) ! Column index array for CSR format
    INTEGER, INTENT(IN) :: ju(nintf)   ! Pointer to the diagonal entry of each row
!    REAL(8), INTENT(IN):: au(nnz) ! Original matrix A values in CSR format
    REAL(8)             :: alu(nnz) ! Output L and U in CSR format (initially set to zero)

    INTEGER :: i, j, k, col         ! Loop counters and temporary variables
    INTEGER :: row_start, row_end   ! Indices for the start and end of each row
    REAL(8) :: sum                  ! Temporary variable to store summations

    ! Step 1: Initialize ALU to be the same as AU
!    alu = au

    ! Step 2: Perform ILU(0) factorization row by row
    DO i = 1, nintf
        row_start = ia(i)
        row_end = ju(i) - 1

        ! Step 2.1: Compute the L entries for this row
        DO k = row_start, row_end
            col = ja(k)
            IF (col < i) THEN
                ! Compute L(i, col) entry
                alu(k) = alu(k) / alu(ju(col))  ! Divide by U(col, col)
                
                ! Update other entries in row with this L(i, col)
                DO j = k + 1, row_end
                    IF (ja(j) == col) THEN
                        alu(j) = alu(j) - alu(k) * alu(ju(col) + (col - i))
                    END IF
                END DO
            END IF
        END DO

        ! Step 2.2: Compute the diagonal U(i, i) entry
        k = ju(i)
        sum = alu(k)
        DO j = row_start, k - 1
            col = ja(j)
            IF (col < i) THEN
                sum = sum - alu(j) * alu(ju(col) + (i-col))
            END IF
        END DO
        alu(k) = sum  ! Store U(i, i)

        ! Step 2.3: Compute the U entries for this row (above diagonal)
        row_start = ju(i) + 1
        DO k = row_start, ia(i + 1) - 1
            col = ja(k)
            IF (col > i) THEN
                DO j = ia(i), k - 1
                    alu(k) = alu(k) - alu(j) * alu(ju(ja(j)) + (col - ja(j)))
                END DO
            END IF
        END DO
!
         j=ju(i)
         alu(j)=1.0/alu(j)
         
    END DO
END SUBROUTINE ilu0_gpt