!**************************************************      - - - - !
! parallel CG using Diagonal Preconditioning - - - - - - - - - - !
! Modify by S.T Ha - - - - - - - - - - - - - - - - - - - - - - - !
! June 2021 - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  SUBROUTINE PCG_Dig(maxit,nintf,neq,niut,nnz,myrank,ia,ja,ju,     &
                  si,ri,sintf,rintf,iut,au,arhsu,solu,crit,ierror)

      USE MD_parameter, only: ndom
      IMPLICIT NONE
!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF
!
      INTEGER maxit,myrank,ierror
	  INTEGER nintf,neq,niut,nnz
	  INTEGER ia(nintf+1),ju(nintf),ja(nnz)
	  INTEGER si(niut+1),ri(niut+1),iut(niut),sintf(si(niut+1)-1),rintf(neq-nintf)
	  REAL*8 crit
	  REAL*8 au(nnz),arhsu(neq)
! 
      REAL*8 solu(neq)
! 
!  temp.
      INTEGER i,j,n,Allocatestatus,iter,np,alstatus,tag

!DEC$IF defined (mpi_flag)
      INTEGER status(MPI_STATUS_SIZE),ierr
!DEC$ENDIF

      REAL*8 bknum,bkden,akden,bk,ak,bknum1,akden1
      REAL*8 err0,err1,err,res0,res1
      REAL*8,DIMENSION(:),ALLOCATABLE::r,z,p0,svar,rvar,dig
      INTEGER request(2*niut)

      n=neq
      np = ndom
      Allocate(r(n),z(n),p0(n),stat=Allocatestatus)
      if(AllocateStatus/=0) stop "**Not enough memory "

!-----Diag PC z=r/Diag(*)
      ALLOCATE(Dig(nintf))
      do i=1,nintf
       dig(i)=au(ju(i))
      enddo	  
	  
!----Predictor  By Diag. pc
!      r=0.d0 !! 
!      do i=1,nintf
!        solu(i)=arhsu(i)/au(ju(i))
!      enddo

!-----Exchange solu to compute A*solu
      allocate(svar(si(niut+1)-1),stat=alstatus)
      if (alstatus/=0) stop 'not enough 12memory'
      allocate(rvar(ri(niut+1)-1),stat=alstatus)
      if (alstatus/=0) stop 'not enough 13memory'
      
      do i=1,si(niut+1)-1
        svar(i)=solu(sintf(i))
      enddo
!	  write(*,*)'nsend',si(1+1)-si(1)
!DEC$IF defined (mpi_flag)
      tag=1
      do i=1,niut
        call MPI_ISEND(svar(si(i)),si(i+1)-si(i),    &
             MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(ri(i)),ri(i+1)-ri(i),    &
             MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+niut),ierr)
      enddo
	  
      do i=1,niut
 !       call MPI_WAIT(request(i),status,ierr)          
        call MPI_WAIT(request(i+niut),status,ierr)
      enddo
!DEC$ENDIF
!
      do i=1,n-nintf
        solu(rintf(i))=rvar(i)
      enddo

      call amux0_PCG(nintf,n,nnz,solu,r,au,ja,ia)
!-----------------------------------------------
!      do i=1,nintf
       r(1:nintf)=arhsu(1:nintf)-r(1:nintf)
!      enddo
	  
!-----Diag PC z=r/Diag(*)
!      do i=1,nintf
       z(1:nintf)=r(1:nintf)/Dig(1:nintf)
!      enddo

      err0 = 0.d0
      do i=1,nintf
        err0=err0+r(i)*r(i)
      enddo
!DEC$IF defined (mpi_flag)
      call MPI_ALLREDUCE(err0,err1,1,MPI_DOUBLE_PRECISION,    &
                     MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF
      IF(np.GT.1) THEN
      err0=dsqrt(err1)
      ELSE
      err0 = DSQRT(err0)
      ENDIF 

      if (err0.LT.1.D-15) goto 502

! loop
      DO iter=1,maxit
        tag = iter 
!---calculmate coefficient bk and direction vector p
        bknum=0.d0
        do i=1,nintf
          bknum=bknum+r(i)*z(i)
        enddo
!DEC$IF defined (mpi_flag)
        call MPI_ALLREDUCE(bknum,bknum1,1,MPI_DOUBLE_PRECISION,   &
                    MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF

      IF(np.GT.1) THEN
      bknum=bknum1
      ELSE
      bknum = bknum
      ENDIF 

        if ( iter.eq.1) then
         p0=z
        else
         bk = bknum/bkden
         p0=z+bk*p0
        endif
        bkden = bknum
		
!    calculmate coefficient ak, new itermate x, new residual r

!-----Exchange p0 to compute A*p0
!DEC$IF defined (mpi_flag)
        do i=1,niut
	    call MPI_WAIT(request(i),status,ierr)
	    enddo	
	    do i=1,si(niut+1)-1
	    svar(i)=p0(sintf(i))
	    enddo

        tag=1
	    do i=1,niut
	    call MPI_ISEND(svar(si(i)),si(i+1)-si(i),      &
               MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
               MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(ri(i)),ri(i+1)-ri(i),      &
               MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
               MPI_COMM_WORLD,request(i+niut),ierr)
        enddo
		
	  do i=1,niut
	  call MPI_WAIT(request(i+niut),status,ierr)
      enddo

!DEC$ENDIF
    
	 do i=1,ri(niut+1)-1
	  p0(rintf(i))=rvar(i)
	 enddo

     call amux0_PCG(nintf,n,nnz,p0,z,au,ja,ia) !!A*p0=z
!-----------------------------------------------------
        akden=0.d0
        do i=1,nintf
         akden=akden+p0(i)*z(i)
        enddo
!DEC$IF defined (mpi_flag)
        call MPI_ALLREDUCE(akden,akden1,1,MPI_DOUBLE_PRECISION,  &
                   MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF
        
      IF(np.GT.1) THEN
      akden1=akden1
      ELSE
      akden1 = akden
      ENDIF
      
        ak = bknum/akden1

        do i=1,nintf
          solu(i) = solu(i) + ak*p0(i)
          r(i) = r(i) - ak*z(i)
        enddo
!-------Diag PC z=r/Diag(*)
!        do i=1,nintf
         z(1:nintf)=r(1:nintf)/Dig(1:nintf)
!        enddo

        err=0.d0
        do i=1,nintf
          err=err+r(i)*r(i)
        enddo
!DEC$IF defined (mpi_flag)
        call MPI_ALLREDUCE(err,err1,1,MPI_DOUBLE_PRECISION,     &
                        MPI_SUM,MPI_COMM_WORLD,ierr)
!DEC$ENDIF
      IF(np.GT.1) THEN
      err=dsqrt(err1)
      ELSE
      err = dsqrt(err)
      ENDIF
		
!--------------------------------------------------
       if (abs(err/err0).gt.1.d10) then
           ierror = 1
          if (myrank.eq.0) then
             write(*,*) 'blow-up in piacg solver, res=',err/err0
          endif  
!          stop
       endif
!       if(myrank.eq.0)write(*,*)'Iter=',iter,'Residual=',err/err0
       if (err/err0.LT.crit) go to 502
!----------------------------------------------------
      ENDDO
      
!      ierror = 1
      if (myrank.eq.0) then
	  write(*,*)'DigCG is divergence',iter,err/err0
      endif
      
502   CONTINUE
!DEC$IF defined (mpi_flag)
      do i=1,niut
        call MPI_WAIT(request(i),status,ierr)
      enddo
  
      do i=1,si(niut+1)-1
        svar(i)=solu(sintf(i))
      enddo
	  
      tag=1
      do i=1,niut
        call MPI_ISEND(svar(si(i)),si(i+1)-si(i),   &
            MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
            MPI_COMM_WORLD,request(i),ierr)
			
        call MPI_IRECV(rvar(ri(i)),ri(i+1)-ri(i),   &
            MPI_DOUBLE_PRECISION,iut(i)-1,tag,      &
            MPI_COMM_WORLD,request(i+niut),ierr)
      enddo
	  
      do i=1,niut
        call MPI_WAIT(request(i+niut),status,ierr)
      enddo
	  
      do i=1,n-nintf
        solu(rintf(i))=rvar(i)
      enddo
	
      IF(myrank.eq.0) WRITE(*,*)'Iter=',iter,'Residual=',err/err0
      DO i=1,niut
        CALL MPI_WAIT(request(i),status,ierr)
      ENDDO
!DEC$ENDIF
! 
      DEALLOCATE(r,z,p0,svar,rvar,Dig)

!
      RETURN
      END

!***********************************************************************
      SUBROUTINE amux0_PCG(nintf,n,nnz,x,y,a,ja,ia)
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
