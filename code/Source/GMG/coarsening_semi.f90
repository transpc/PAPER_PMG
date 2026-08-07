  subroutine coarsening_semi(ndim,nnodf,nmax,nnz,ia,ja,nnodc,imap,imapinv,teta,coord)
!
USE MD_MPI, ONLY: myrank

  implicit none

! inlet
  Integer(4)  nnodf,nnz,ndim,nmax
  integer(4)  ia(nnodf+1),ja(nnz)
  REAL(8) coord(ndim,nnodf)  
  REAL(8) teta
! outlet
  Integer (kind=4) nnodc
  Integer (kind=4) imap(nnodc),imapinv(nnodf)
! temporary
  INTEGER(4) alstatus
  integer (kind = 4) nvpe, nn, ne,i,j,k,m,n,id,k1,k2,j1,j2,nnd,k3,j3,ii
  INTEGER (kind = 4) mi(nmax),ni(nmax)
  REAL(8) dx(nmax),x0,y0,z0,xc,yc,zc,dmin
  integer (kind = 4),dimension(:),allocatable :: i_test
	
! ------------------------------
  allocate(i_test(nnodf),stat=alstatus)
  
     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-coarsening-i_test'
         STOP
     ENDIF
  
   nnd = 0
   i_test = 2

! - - - - - - - - - - - - - - - - - - - - -  
! step 1: coarsening nodes by MIS
! - - - - - - - - - - - - - - - - - - - - -
   
 DO j = 1,nnodf
 
  IF(i_test(j) == 0) CYCLE  
  IF(i_test(j) == 1) CYCLE 
   nnd = nnd + 1
   i_test(j) = 1
   
! remove the neighbour nodes
! based on distant:

   j1 = ia(j)
   j2 = ia(j+1)-1
   
     nn = j2-j1+1
     
     IF(nn.EQ.0) CYCLE
     
	 x0 = coord(1,j)
	 y0 = coord(2,j)
     
     IF(ndim.EQ.3) THEN
	 z0 = coord(3,j)
     ENDIF
     
	 
	  dmin = 1.d10
      
      DO k = 1,nn
          
      id = ja(k+j1-1) 
      ni(k) = id
! 
      xc = coord(1,id)
	  yc = coord(2,id)
      
      IF(ndim.eq.3) THEN
	  zc = coord(3,id)
      ENDIF
      
	  IF(ndim.EQ.2) THEN
	  dx(k) = DSQRT((xc-x0)**2.d0+(yc-y0)**2.d0)
      ELSE
	  dx(k) = DSQRT((xc-x0)**2.d0+(yc-y0)**2.d0+(zc-z0)**2.d0)
      ENDIF
      
	  dmin = min(dmin,dx(k))
      ENDDO
!
      dmin = 1.d0/teta*dmin
	  
      DO k=1,nn
      id = ja(k+j1-1) 
!        id = iworkf(k,j)
        IF(i_test(id).EQ.0) CYCLE
        IF(i_test(id).EQ.1) CYCLE !
		
        IF(dx(k).LE.dmin) THEN
        i_test(id) = 0 
        
!        ELSE
! NEW
!        i_test(id) = 1
        ENDIF
        

      ENDDO
! 
  End do
! 
!  write(*,*) 'initial coarse nodes',nnd

! modifying: Change some nodes from F-> C (no connection)
    
  nnd = 0
 DO j = 1,nnodf
 
  IF(i_test(j) == 1) THEN  
  
   nnd = nnd + 1
   imap(nnd) = j
   imapinv(j) = nnd
   
  ELSE
   j1 = ia(j)
   j2 = ia(j+1)-1
   
!     nn = j2-j1+1  
     
!      nn = iwkf(j)
      DO k=j1,j2      !1,nn
      id = ja(k)
      m=i_test(id)
      IF(m.EQ.1) GOTO 11
      ENDDO

   nnd = nnd + 1
   imap(nnd) = j
   imapinv(j) = nnd
   
   i_test(j) = 1
   
  ENDIF    
      
11 CONTINUE
      
 ENDDO
  
 IF(myrank == 0) THEN
  WRITE(999,*)'number of coarse nodes',nnd     ! notes this is called 2 times (first is serial and second is for parallel)
 ENDIF
 
  
  IF(nnd.GT.nnodc) THEN
      WRITE(999,*)'initial nelem1 is small,nnd=',nnd,nnodc
!      PAUSE
      STOP
  ELSE
      
  nnodc = nnd
  ENDIF
   
	
  Deallocate(i_test)
! 

  Return
  End
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
