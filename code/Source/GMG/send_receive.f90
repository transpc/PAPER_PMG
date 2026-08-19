   SUBROUTINE send_receive(nnbd,nnode,spt,rpt,sintf,rintf,nbdom,u)
    
   USE MD_MPI, ONLY: myrank
   IMPLICIT NONE
      
!DEC$IF defined (mpi_flag)
   INCLUDE 'mpif.h'
!DEC$ENDIF

   Integer nnbd,nnode
   Integer spt(nnbd+1),rpt(nnbd+1),sintf(spt(nnbd+1)-1),rintf(rpt(nnbd+1)-1),nbdom(nnbd)
   Real*8 u(nnode)
!
   integer::tag,ierr,request(2*nnbd)
          
!DEC$IF defined (mpi_flag)          
       integer::status(mpi_status_size)
!DEC$ENDIF
! temp:
   Integer i,isend,nsend,psend,irecv,nrecv,precv
 !  REAL*8,DIMENSION(:),ALLOCATABLE:: svar,rvar
   Real*8 svar(spt(nnbd+1)),rvar(rpt(nnbd+1))
!   allocate( svar(spt(nnbd+1)-1),rvar(rpt(nnbd+1)-1))
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !%copy data to temporary array

   do i=1,spt(nnbd+1)-1
      svar(i)=u(sintf(i))
   enddo
!
   svar(spt(nnbd+1)) = 0.d0
!DEC$IF defined (mpi_flag) 
!-------------------------------------------------------------------
!   call mpi_barrier(mpi_comm_world,ierr)
   
         tag=1
      do i=1,nnbd
        call MPI_ISEND(svar(spt(i)),spt(i+1)-spt(i),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(rpt(i)),rpt(i+1)-rpt(i),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo
      
!   call mpi_barrier(mpi_comm_world,ierr)   
 !DEC$ENDIF  
   do i=1,rpt(nnbd+1)-1
      u(rintf(i))=rvar(i)
   enddo

!-------------------------------------------------------------------

  ! DEALLOCATE(svar,rvar)
   return
   
    End subroutine
    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    
!DIR$ ATTRIBUTES NOINLINE :: send_receive_C
   SUBROUTINE send_receive_C(nnbd,ista,nnsend,nnrecv,spt,rpt,sintf,rintf,nbdom,u,ndom,nnsend_m,nnrecv_m)
    
   USE MD_MPI, ONLY: myrank
   IMPLICIT NONE
      
!DEC$IF defined (mpi_flag)
   INCLUDE 'mpif.h'
!DEC$ENDIF

   Integer nnbd,ista,nnsend,nnrecv,ndom,nnsend_m,nnrecv_m
   Integer spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m),nbdom(ndom)
   REAL(8) u(*)
   integer::tag,ierr,request(2*nnbd)
          
!DEC$IF defined (mpi_flag)          
       integer::status(mpi_status_size)
!DEC$ENDIF
! temp:
   Integer i,isend,nsend,psend,irecv,nrecv,precv,itmp
   Real(8) svar(nnsend+1),rvar(nnrecv+1)
   !-------------------------------------------------------------------
   !-------------------------------------------------------------------
   !%copy data to temporary array
   itmp = ista - 1
   do i=1,spt(nnbd+1)-1
      svar(i)=u(sintf(i) + itmp )
   enddo
!
   svar(nnsend+1) = 0.d0
   rvar(nnrecv+1) = 0.d0
!DEC$IF defined (mpi_flag) 
!-------------------------------------------------------------------
!   call mpi_barrier(mpi_comm_world,ierr) 
   
         tag=1
      do i=1,nnbd
        call MPI_ISEND(svar(spt(i)),spt(i+1)-spt(i),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar(rpt(i)),rpt(i+1)-rpt(i),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo

!   call mpi_barrier(mpi_comm_world,ierr) 
   
 !DEC$ENDIF  
   do i=1,rpt(nnbd+1)-1
      u(rintf(i) + itmp )=rvar(i)
   enddo

!-------------------------------------------------------------------

   return
   
    End subroutine
    
! = = = = = = = = = = = 
    SUBROUTINE MD_S_R_NEW(id,ilv,ista,u)
    
USE MD_MPI_ARP, ONLY: inbdcA, ibdomcA, isptcA, irptcA, isintfcA, irintfcA,  &
                      inbdcR, ibdomcR, isptcR, irptcR, isintfcR, irintfcR,  &
                      inbdcP, ibdomcP, isptcP, irptcP, isintfcP, irintfcP,  &
                      nnsend_mA,nnrecv_mA
USE MD_parameter, ONLY: ndom

IMPLICIT NONE

INTEGER(4)::ilv,ista,id
REAL(8) u(*)
! temp:
INTEGER(4)::nnbd,nnsend,nnrecv
!INTEGER(4),DIMENSION(:),ALLOCATABLE::nbdom,spt,rpt,sintf,rintf
INTEGER(4)::nbdom(ndom),spt(ndom),rpt(ndom),sintf(nnsend_mA),rintf(nnrecv_mA)

!ALLOCATE(nbdom(ndom),spt(ndom),rpt(ndom),sintf(nnsend_mA),rintf(nnrecv_mA))
!-------------------------------------------------------------------
!nbdom = 0
!spt = 0
!rpt = 0
!sintf = 0
!rintf = 0

!-------------------------------------------------------------------
   IF(id.EQ.1) THEN
   nnbd = inbdcA(ilv) 
   IF(nnbd.LE.0) GOTO 100
   
   nbdom(1:nnbd) = ibdomcA(1:nnbd,ilv)
   spt(1:nnbd+1) = isptcA(1:nnbd+1,ilv)
   rpt(1:nnbd+1) = irptcA(1:nnbd+1,ilv)
   
   nnsend = spt(nnbd+1) - 1
   nnrecv = rpt(nnbd+1) - 1
   
   sintf(1:nnsend) = isintfcA(1:nnsend,ilv)
   rintf(1:nnrecv) = irintfcA(1:nnrecv,ilv)
   ELSEIF(id.EQ.2) THEN
       
   nnbd = inbdcR(ilv) 
   IF(nnbd.LE.0) GOTO 100
   
   nbdom(1:nnbd) = ibdomcR(1:nnbd,ilv)
   spt(1:nnbd+1) = isptcR(1:nnbd+1,ilv)
   rpt(1:nnbd+1) = irptcR(1:nnbd+1,ilv)
   
   nnsend = spt(nnbd+1) - 1
   nnrecv = rpt(nnbd+1) - 1
   
   sintf(1:nnsend) = isintfcR(1:nnsend,ilv)
   rintf(1:nnrecv) = irintfcR(1:nnrecv,ilv)
   ELSE 
       
   nnbd = inbdcP(ilv) 
   IF(nnbd.LE.0) GOTO 100
   
   nbdom(1:nnbd) = ibdomcP(1:nnbd,ilv)
   spt(1:nnbd+1) = isptcP(1:nnbd+1,ilv)
   rpt(1:nnbd+1) = irptcP(1:nnbd+1,ilv)
   
   nnsend = spt(nnbd+1) - 1
   nnrecv = rpt(nnbd+1) - 1
   
   sintf(1:nnsend) = isintfcP(1:nnsend,ilv)
   rintf(1:nnrecv) = irintfcP(1:nnrecv,ilv)
       ENDIF
       
   
   IF(nnbd.NE.0) THEN
   CALL send_receive_C(nnbd,ista,nnsend,nnrecv,spt,rpt,sintf,rintf,nbdom,u,ndom,nnsend_mA,nnrecv_mA)
   ENDIF
   
100 CONTINUE
    
!-------------------------------------------------------------------
!DEALLOCATE(nbdom,spt,rpt,sintf,rintf)
return
   
End subroutine
   

