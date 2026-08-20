   SUBROUTINE send_receive_csr(ndim,nmax,nnbd,nnode,nnodegl,spt,rpt,sintf,rintf,          &
          nnzs,nnzr,nnsend,nnrecv,sia,ria,nbdom,coord,nnz,ia,ja)
    
    use MD_MPI, only: myrank
   Implicit NONE

!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF

   INTEGER(4):: nnbd,nnode,nnz,ndim,nmax,nnodegl,nnsend,nnrecv,nnzs,nnzr
   INTEGER(4):: spt(*),rpt(*),sintf(*),rintf(*),nbdom(*)
   INTEGER(4):: sia(*),ria(*),ia(*),ja(*)
   REAL(8) coord(ndim,nnodegl)

   INTEGER(4)::tag,ierr,request(2*nnbd)
   
!DEC$IF defined (mpi_flag)          
   INTEGER(4):: status(mpi_status_size)
!DEC$ENDIF

! temp:
   integer i,j,i1,i2,j1,j2,nd
   Integer isend,nsend,psend,irecv,nrecv,precv
   REAL(8),DIMENSION(:),ALLOCATABLE::svar1,svar2,svar3,rvar1,rvar2,rvar3
   REAL(8) xtmp1(ndim,nmax),xtmp2(ndim,nmax)
   INTEGER ni(nmax),ne(nmax)
   ALLOCATE(svar1(nnzs),svar2(nnzs),svar3(nnzs),rvar1(nnzr),rvar2(nnzr),rvar3(nnzr))
   svar1 = 0.d0
   svar2 = 0.d0
   svar3 = 0.d0
   rvar1 = 0.d0
   rvar2 = 0.d0
   rvar3 = 0.d0
   !-------------------------------------------------------------------
   !%MPI send & receive data
   !DEC$IF defined (mpi_flag) 
!   call mpi_barrier(mpi_comm_world,ierr)
   !-------------------------------------------------------------------
   !%copy data to temporary array
   do i=1,spt(nnbd+1)-1
      i1 = sia(i)
	  i2 = sia(i+1)-1
	  
	  j = sintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
	  
      svar1(i1:i2)=coord(1,(ja(j1:j2)))
      svar2(i1:i2)=coord(2,(ja(j1:j2)))

      IF(ndim.EQ.3) THEN 
      svar3(i1:i2)=coord(3,(ja(j1:j2)))
      ENDIF
	  
   enddo

!   call mpi_barrier(mpi_comm_world,ierr)
  !         
         tag=1
      do i=1,nnbd
        call MPI_ISEND(svar1(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        
        call MPI_IRECV(rvar1(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
        
      enddo
      
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo

! y-coord:
   tag=1

      do i=1,nnbd
        call MPI_ISEND(svar2(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar2(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo   
   
! Z-comp
   
   IF(ndim.EQ.3) THEN
   tag=1

       do i=1,nnbd
        call MPI_ISEND(svar3(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar3(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo      
       
   ENDIF
!
   do i=1,rpt(nnbd+1)-1
   
      i1 = ria(i)
	  i2 = ria(i+1)-1
      
      nd = i2-i1+1
      IF(nd.EQ.1) CYCLE
      xtmp2(1,1:nd) = rvar1(i1:i2)
      xtmp2(2,1:nd) = rvar2(i1:i2)
      
      IF(ndim.EQ.3) THEN
      xtmp2(3,1:nd) = rvar3(i1:i2)
      ENDIF
	  
	  j = rintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
      ni(1:nd) = ja(j1:j2)
      ne(1:nd)=[1:nd]
      
      xtmp1(1,1:nd) = coord(1,ni(1:nd))
      xtmp1(2,1:nd) = coord(2,ni(1:nd))
 
      IF(ndim.EQ.3) THEN      
      xtmp1(3,1:nd) = coord(3,ni(1:nd))
      ENDIF

      call ordering_csr(ndim,nd,nmax,ne,xtmp1,xtmp2)
      ja(j1:j2) = ni(ne(1:nd))

   enddo
!   call mpi_barrier(mpi_comm_world,ierr)

!DEC$ENDIF
   !-------------------------------------------------------------------

      DEALLOCATE(svar1,svar2,svar3,rvar1,rvar2,rvar3)

   return
   
    End subroutine
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
      subroutine ordering_csr(ndim,nd,nmax,ne,coord1,coord2)
      implicit none
	  integer ndim,nd,nmax
      integer ne(nmax)
	  REAL*8 coord1(ndim,nmax),coord2(ndim,nmax)
!
      integer id,i,j,imark(nd),j0
	  REAL*8 eps,eps0,xtemp1(ndim),xtemp2(ndim),eps1
! ---------------------------------------------------------!
			
		eps0 = 1.d-10
                eps1 = 1.d10
		imark = 0
		
		DO i = 1,nd
		   xtemp2(1:ndim) = coord2(1:ndim,i)
                   eps1 = 1.d10
		   
		   DO j = 1,nd
		   IF(imark(j).EQ.1) CYCLE
		   
		   xtemp1(1:ndim) = coord1(1:ndim,j)
		   
		   eps = (xtemp1(1)-xtemp2(1))**2.d0+(xtemp1(2)-xtemp2(2))**2.d0
		   IF(ndim.EQ.3) eps = eps +(xtemp1(3)-xtemp2(3))**2.d0
		 
                   IF(eps1.GT.eps) THEN
                   eps1 = eps
                   j0 = j
                   endif
		   IF(eps.LT.eps0) GOTO 10
		   
		   ENDDO
		   
		   WRITE(*,*)'error in ordering_csr',i,j0,eps1
                   write(*,*)xtemp2(1:ndim)
                   write(*,*)coord1(1:ndim,j0)
                   stop
		   
10        CONTINUE

          ne(i) = j
		  imark(j) = 1
		  
		ENDDO
        
! ---
        return
    end
    
! = = = = = = = = 
   SUBROUTINE send_receive_csrc(ndim,nmax,nnbd,nnode,nnodegl,spt,rpt,sintf,rintf,          &
          nnzs,nnzr,nnsend,nnrecv,sia,ria,nbdom,coord,nnz,ia,ja,ndom,nnsend_m,nnrecv_m)
    
    use MD_MPI, only: myrank
   Implicit NONE

!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF

   INTEGER(4):: nnbd,nnode,nnz,ndim,nmax,nnodegl,nnsend,nnrecv,nnzs,nnzr,ndom,nnsend_m,nnrecv_m
   INTEGER(4):: spt(ndom),rpt(ndom),sintf(nnsend_m),rintf(nnrecv_m),nbdom(ndom)
   INTEGER(4):: sia(nnsend_m+1),ria(nnrecv_m+1)
   INTEGER(4):: ia(*),ja(*)
   REAL(8) coord(ndim,nnodegl)

   INTEGER(4)::tag,ierr,request(2*nnbd)
   
   integer::alstatus
   
!DEC$IF defined (mpi_flag)          
   INTEGER(4):: status(mpi_status_size)
!DEC$ENDIF

! temp:
   integer i,j,i1,i2,j1,j2,nd
   Integer isend,nsend,psend,irecv,nrecv,precv
   REAL(8),DIMENSION(:),ALLOCATABLE::svar1,svar2,svar3,rvar1,rvar2,rvar3
   REAL(8) xtmp1(ndim,nmax),xtmp2(ndim,nmax)
   INTEGER ni(nmax),ne(nmax)
   ALLOCATE(svar1(nnzs),svar2(nnzs),svar3(nnzs),rvar1(nnzr),rvar2(nnzr),rvar3(nnzr),stat=alstatus)
   
     IF (alstatus/=0) THEN
         WRITE(*,*)'not allocated, var1',myrank
     ENDIF
     
   svar1 = 0.d0
   svar2 = 0.d0
   svar3 = 0.d0
   rvar1 = 0.d0
   rvar2 = 0.d0
   rvar3 = 0.d0
   
   !-------------------------------------------------------------------
   !%MPI send & receive data
   !DEC$IF defined (mpi_flag) 
!   call mpi_barrier(mpi_comm_world,ierr)
   !-------------------------------------------------------------------
   !%copy data to temporary array
   do i=1,spt(nnbd+1)-1
      i1 = sia(i)
	  i2 = sia(i+1)-1
	  
	  j = sintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
	  
      svar1(i1:i2)=coord(1,(ja(j1:j2)))
      svar2(i1:i2)=coord(2,(ja(j1:j2)))

      IF(ndim.EQ.3) THEN 
      svar3(i1:i2)=coord(3,(ja(j1:j2)))
      ENDIF
	  
   enddo

!   call mpi_barrier(mpi_comm_world,ierr)
  !         
         tag=1
      do i=1,nnbd
          i1 = sia(spt(i))
          i2 = sia(spt(i+1))
          j1 = ria(rpt(i))
          j2 = ria(rpt(i+1))
          
        call MPI_ISEND(svar1(i1),i2-i1,    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        
        call MPI_IRECV(rvar1(j1),j2-j1,    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
        
      enddo
      
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo

! y-coord:
   tag=1

      do i=1,nnbd
        call MPI_ISEND(svar2(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar2(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo   
   
! Z-comp
   
   IF(ndim.EQ.3) THEN
   tag=1

       do i=1,nnbd
        call MPI_ISEND(svar3(sia(spt(i))),sia(spt(i+1))-sia(spt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i),ierr)
        call MPI_IRECV(rvar3(ria(rpt(i))),ria(rpt(i+1))-ria(rpt(i)),    &
             MPI_DOUBLE_PRECISION,nbdom(i)-1,tag,      &
             MPI_COMM_WORLD,request(i+nnbd),ierr)
      enddo
	  
      do i=1,nnbd
        call MPI_WAIT(request(i),status,ierr)
        call MPI_WAIT(request(i+nnbd),status,ierr)
      enddo      
       
   ENDIF
!
   do i=1,rpt(nnbd+1)-1
   
      i1 = ria(i)
	  i2 = ria(i+1)-1
      
      nd = i2-i1+1
      IF(nd.EQ.1) CYCLE
      xtmp2(1,1:nd) = rvar1(i1:i2)
      xtmp2(2,1:nd) = rvar2(i1:i2)
      
      IF(ndim.EQ.3) THEN
      xtmp2(3,1:nd) = rvar3(i1:i2)
      ENDIF
	  
	  j = rintf(i)
	  j1 = ia(j)
	  j2 = ia(j+1)-1
      ni(1:nd) = ja(j1:j2)
      ne(1:nd)=[1:nd]
      
      xtmp1(1,1:nd) = coord(1,ni(1:nd))
      xtmp1(2,1:nd) = coord(2,ni(1:nd))
 
      IF(ndim.EQ.3) THEN      
      xtmp1(3,1:nd) = coord(3,ni(1:nd))
      ENDIF

      call ordering_csr(ndim,nd,nmax,ne,xtmp1,xtmp2)
      ja(j1:j2) = ni(ne(1:nd))

   enddo
!   call mpi_barrier(mpi_comm_world,ierr)

!DEC$ENDIF
   !-------------------------------------------------------------------

      DEALLOCATE(svar1,svar2,svar3,rvar1,rvar2,rvar3)

   return
   
    End subroutine
   

