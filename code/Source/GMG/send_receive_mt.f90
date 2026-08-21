   
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
   SUBROUTINE send_receive_mtf(nnbd,nnsend,nnrecv,spt,rpt,sintf,rintf,          &
              sia,ria,nbdom,ia,au)
    
   Implicit NONE

!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF

   INTEGER(4):: nnbd,nnsend,nnrecv
   INTEGER(4):: spt(*),rpt(*),sintf(*),rintf(*),nbdom(*)
   INTEGER(4):: sia(*),ria(*)
   INTEGER(4):: ia(*)
   REAL(8):: au(*)
   INTEGER(4):: tag,ierr,request(2*nnbd)
   
!DEC$IF defined (mpi_flag)          
   INTEGER(4):: status(mpi_status_size)
!DEC$ENDIF

! temp:
   INTEGER(4):: i,i1,i2,j,j1,j2,nnzs,nnzr
   REAL(8),DIMENSION(:),ALLOCATABLE:: svar,rvar
! 
   nnzs = sia(nnsend+1)-1
   nnzr = ria(nnrecv+1)-1
   
   ALLOCATE( svar(nnzs),rvar(nnzr))
   svar = 0.d0
   rvar = 0.d0
   !-------------------------------------------------------------------
   !%MPI send & receive data
   !%send & receive data
   !DEC$IF defined (mpi_flag)
   !-------------------------------------------------------------------
   !%copy data to temporary array
   DO i=1,spt(nnbd+1)-1
      i1 = sia(i)
	  i2 = sia(i+1)-1
	  
	  j = sintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
	  
      svar(i1:i2)=au(j1:j2)
	  
   ENDDO

   !-------------------------------------------------------------------
   tag=1
         do i=1,nnbd
        call MPI_ISEND(svar(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo 
   
!   
   do i=1,rpt(nnbd+1)-1
   
      i1 = ria(i)
	  i2 = ria(i+1)-1
	  
	  j = rintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
   
      au(j1:j2)=rvar(i1:i2)
   enddo
   
   
!DEC$ENDIF
   !-------------------------------------------------------------------
   DEALLOCATE(svar,rvar)
   return
   
    End subroutine    
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    
SUBROUTINE MD_S_R_MT(ilv,ia,au)
    
USE MD_MPI_MG, ONLY: inbdc,ibdomc,isintfc,irintfc,           &
                     nnsend_m,nnrecv_m,isptc,irptc,          &
                     isiac,iriac
USE MD_parameter, ONLY: ndom

IMPLICIT NONE

INTEGER(4)::ilv
INTEGER(4)::ia(*)
REAL(8) au(*)
! temp:
INTEGER(4)::nnbd,nnsend,nnrecv
INTEGER(4),DIMENSION(:),ALLOCATABLE::nbdom,spt,rpt,sintf,rintf,sia,ria
ALLOCATE(nbdom(ndom),spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m))
ALLOCATE(sia(nnsend_m+1),ria(nnrecv_m+1))
!-------------------------------------------------------------------
   nbdom = 0
   spt = 0
   rpt = 0
   sintf = 0
   rintf = 0
   sia = 0
   ria = 0
! 
   nnbd = inbdc(ilv)
   nbdom(1:nnbd) = ibdomc(1:nnbd,ilv)
   spt(1:nnbd+1) = isptc(1:nnbd+1,ilv)
   rpt(1:nnbd+1) = irptc(1:nnbd+1,ilv)
   
   nnsend = spt(nnbd+1) - 1
   nnrecv = rpt(nnbd+1) - 1
   
! 
   IF(nnsend.GT.nnsend_m.OR.nnrecv.GT.nnrecv_m) THEN
       WRITE(*,*)'nnsend_m is small-1',nnsend,nnrecv,nnsend_m,nnrecv_m
       STOP
   ENDIF
   
   sintf(1:nnsend) = isintfc(1:nnsend,ilv)
   rintf(1:nnrecv) = irintfc(1:nnrecv,ilv)
   
   sia(1:nnsend+1) = isiac(1:nnsend+1,ilv)
   ria(1:nnrecv+1) = iriac(1:nnrecv+1,ilv)
   
CALL send_receive_mtc(nnbd,nnsend,nnrecv,spt,rpt,sintf,rintf,                 &
              sia,ria,nbdom,ia,au,ndom,nnsend_m,nnrecv_m)

!-------------------------------------------------------------------
DEALLOCATE(nbdom,spt,rpt,sintf,rintf,sia,ria)

return
   
    End subroutine
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
   SUBROUTINE send_receive_mtc(nnbd,nnsend,nnrecv,spt,rpt,sintf,rintf,          &
              sia,ria,nbdom,ia,au,ndom,nnsend_m,nnrecv_m)
    
   Implicit NONE

!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF

   INTEGER(4):: nnbd,nnsend,nnrecv,ndom,nnsend_m,nnrecv_m
   INTEGER(4):: spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m),nbdom(ndom)
   INTEGER(4):: sia(nnsend_m+1),ria(nnrecv_m+1)
   INTEGER(4):: ia(*)
   REAL(8):: au(*)
   INTEGER(4):: tag,ierr,request(2*nnbd)
   
!DEC$IF defined (mpi_flag)          
   INTEGER(4):: status(mpi_status_size)
!DEC$ENDIF

! temp:
   INTEGER(4):: i,i1,i2,j,j1,j2,nnzs,nnzr
   REAL(8),DIMENSION(:),ALLOCATABLE:: svar,rvar
! 
   nnzs = sia(nnsend+1)-1
   nnzr = ria(nnrecv+1)-1
   
   ALLOCATE( svar(nnzs),rvar(nnzr))
   svar = 0.d0
   rvar = 0.d0
   !-------------------------------------------------------------------
   !%MPI send & receive data
   !DEC$IF defined (mpi_flag)
   !-------------------------------------------------------------------
   !%copy data to temporary array
   DO i=1,spt(nnbd+1)-1
      i1 = sia(i)
	  i2 = sia(i+1)-1
	  
	  j = sintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
	  
      svar(i1:i2)=au(j1:j2)
	  
   ENDDO

   !-------------------------------------------------------------------
   !%send & receive data

   tag=1
         do i=1,nnbd
        call MPI_ISEND(svar(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo 
   
!   
   do i=1,rpt(nnbd+1)-1
   
      i1 = ria(i)
	  i2 = ria(i+1)-1
	  
	  j = rintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
   
      au(j1:j2)=rvar(i1:i2)
   enddo
   
!DEC$ENDIF   
   !-------------------------------------------------------------------
   DEALLOCATE(svar,rvar)
   return
   
    End subroutine    
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
