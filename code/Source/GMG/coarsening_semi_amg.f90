!  subroutine coarsening_semi_amg(ndim,nnodf,mxnbne,iwkf,iworkf,nnodc,imap,imapinv,teta,coord)
  subroutine coarsening_semi_amg(ndim,nnodf,mxnbne,nnz,ia,ja,nnodc,imap,imapinv,teta,coord)  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! this coarsening is based on amg concept - - - - - - - - - - - - - - - !
! 1: using the seed-point (the first cells on coarse level) - - - - - - !
! 2: coarsening process is based on maximunn connectivity  - - - - - - -!
! the selected cells are from the front-set (updated every cells) - - - !
! 3: some cells from F -> C because of strong-influence definition  - - !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

  implicit none

! inlet
  Integer (kind=4) nnodf,mxnbne,ndim,nnz
!  integer (kind=4) iwkf(nnodf),iworkf(mxnbne,nnodf)
  integer(4)  ia(nnodf+1),ja(nnz)
  REAL(8) coord(ndim,nnodf)  
  REAL(8) teta

  
! NEW
  INTEGER(4) iamg
! iamg = 0-> using all connectivity:

! outlet
  Integer (kind=4) nnodc
  Integer (kind=4) imap(nnodc),imapinv(nnodf)
! temporary
  INTEGER(4) alstatus
  integer (kind = 4) nvpe, nn, ne,i,j,k,m,n,id,k1,k2,j1,j2,nnd,k3,j3,ii
  INTEGER (kind = 4) mi(mxnbne),ni(mxnbne)
  REAL(8) dx(mxnbne),x0,y0,z0,xc,yc,zc,dmin
  integer (kind = 4),dimension(:),allocatable :: i_test,imark
  
  integer (kind=4),dimension(:),allocatable:: iwkf
  integer (kind=4),dimension(:,:),allocatable::  iworkf  
! NEW
  INTEGER(4) lamda_m,isp,numC,numF,nft,nftm,jd,nft0
  INTEGER(4),DIMENSION(:,:),ALLOCATABLE :: iworks
  INTEGER(4),DIMENSION(:),ALLOCATABLE :: ifront,id1
  
	
! ------------------------------
  allocate(i_test(nnodf),stat=alstatus)
  
     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-coarsening-i_test'
         STOP
     ENDIF
  
! iwk---
    allocate(iwkf(nnodf),iworkf(mxnbne,nnodf))
	do i = 1, nnodf
	j1 = ia(i)
	j2 = ia(i+1)
	
	nn = j2-j1
!
    iwkf(i) = nn
    iworkf(1:nn,i) = ja(j1:j2-1)
	ENDDO
	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
! step-0: defind strong-cells: 
    ALLOCATE(iworks(mxnbne,nnodf))
	iworks = 0
  
  DO j = 1,nnodf
     nn = iwkf(j)
     
     IF(nn.EQ.0) THEN
         WRITE(*,*)'no connection'
         STOP
     ENDIF
     
	 x0 = coord(1,j)
	 y0 = coord(2,j)
     
     IF(ndim.EQ.3) THEN
	 z0 = coord(3,j)
     ENDIF
	 
	  dmin = 1.d10
      
      DO k = 1,nn
	  
      id = iworkf(k,j) 
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
      
        IF(dx(k).LE.dmin) iworks(k,j) = 1
		
      ENDDO
! 
  End do
! - - - - - - - - - - - - - - - - - - - - - 
  
! step-1: finding a seed-point:

    lamda_m = maxval(iwkf)
	
   DO i = 1,nnodf
      m = iwkf(i)
	  IF(m.EQ.lamda_m) THEN
	    isp = i
		EXIT
		ENDIF
   ENDDO


! step-2: coloring cells: C-cells and F-cells. 

   nftm = INT((nnodf)**(0.87))
   nftm = 5*nftm                 ! notes
   
   ALLOCATE(imark(nnodf))
   imark = 0
   
   ALLOCATE(ifront(nftm))
   ifront = 0
   
   nnd = 0
   i_test = 2
   nft = 1
   ifront(1) = isp
   
10 CONTINUE

! select seed-point from front:
   isp = ifront(1)
   
   DO i=1,nft-1
   ifront(i) = ifront(i+1)
   ENDDO
   ifront(nft) = 0
   nft = nft -1
!
   nnd = nnd + 1
   i_test(isp) = 1
   
! step 2.1: remove the neighbour nodes (strong cells)

     nn = iwkf(isp)
      
    DO k = 1,nn
	  
      id = iworkf(k,isp) 
	  
	  IF(i_test(id).NE.2) CYCLE
	  
	  IF(iworks(k,isp).EQ.0) THEN     ! weak cell
	  
	    IF(imark(id).EQ.1) CYCLE
	    CALL add2front(id,nft,nftm,ifront,nnodf,iwkf) ! add this weak cell to front
		imark(id) = 1
		
	  ELSE 
		
		i_test(id) = 0 
		
	    DO j = 1,iwkf(id)         ! add neighnor of this cell to front 
		  jd = iworkf(j,id)
		  
		  IF(i_test(jd).NE.2) CYCLE
          
	      IF(imark(jd).EQ.1) CYCLE		  
	      CALL add2front(jd,nft,nftm,ifront,nnodf,iwkf) 	
		  imark(jd) = 1		  
		  
		ENDDO
		
	  ENDIF

    ENDDO

! step-2.2: updating front.
	
    CALL update_front(nft,nftm,ifront,nnodf,i_test)
    
	IF(nft.EQ.0) GOTO 100
    
! extend size if necessary
    i = mxnbne**2
    nft0 = max(nftm-i,i)
    IF(nft.GE.nft0) THEN
	WRITE(*,*)'nftm is small'
!   
    ALLOCATE(id1(nftm))
    id1 = ifront
    DEALLOCATE(ifront)
    i = INT(1.2*nftm)
    ALLOCATE(ifront(i))
    ifront(1:nftm) = id1(1:nftm)
    ifront(nftm+1:i) = 0
    
    DEALLOCATE(id1)
    nftm = i
    
    ENDIF
    
    GOTO 10
    
100  CONTINUE
! 
  write(*,*) 'initial coarse nodes',nnd
  DEALLOCATE(ifront,iworks,imark)
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! modifying: Change some nodes from F-> C (no connection)
    
   nnd = 0
   imap = 0
   imapinv = 0
   numC = 0
   numF = 0
   
 DO j = 1,nnodf
 
  IF(i_test(j) == 1) THEN  
  
   nnd = nnd + 1
   imap(nnd) = j
   imapinv(j) = nnd
   numC = numC+1
   
  ELSE
   numF = numF+1
      nn = iwkf(j)
      DO k=1,nn
      id = iworkf(k,j)
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
  
  WRITE(*,*)'number of coarse nodes',nnd
  
  IF((numF+numC).NE.nnodf) THEN
  WRITE(*,*)'coarsening error,',numF+numC,nnodf
  STOP
  ENDIF
  
  IF(nnd.GT.nnodc) THEN
      WRITE(*,*)'initial nelem1 is small,nnd=',nnd,nnodc
      STOP
  ELSE
      
  nnodc = nnd
  ENDIF
   
	
  Deallocate(i_test)
  Deallocate(iwkf,iworkf)
! 

  Return
  End
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
   SUBROUTINE add2front(id,nft,nftm,ifront,nnode,iwk)
   
! add node id to array ifront:
   IMPLICIT NONE
   
   INTEGER(4) id,nft,nftm,nnode
   INTEGER(4) ifront(nftm),iwk(nnode)
   
! temporary
   INTEGER(4) i,jd,nk,nn
   INTEGER(4),DIMENSION(:),ALLOCATABLE::ifront0
   
   
!
    ALLOCATE(ifront0(nft+1))
	ifront0(1:nft) = ifront(1:nft)
	
      nk = iwk(id)
	  
   DO i = 1,nft
      jd = ifront(i)
	  nn = iwk(jd)
	  
	IF(nk.LT.nn) CYCLE
	IF(nk.EQ.nn.AND.id.GT.jd) CYCLE
! add
    ifront(i) = id
	ifront(i+1:nft+1) = ifront0(i:nft)
	
	GOTO 11
	
   ENDDO
   
   ifront(nft+1) = id
   
11 CONTINUE

   nft = nft + 1
   
   DEALLOCATE(ifront0)
   
   RETURN 
   
   END SUBROUTINE
   
! - - - - - - - - - - - - - - - - - - - - - - - - - !
    SUBROUTINE update_front(nft,nftm,ifront,nnode,imark)
	
	IMPLICIT NONE
	
	INTEGER(4) nft,nftm,nnode
	INTEGER(4) ifront(nftm),imark(nnode)
	
! temporary
   INTEGER(4) i,jd,nnd
   INTEGER(4),DIMENSION(:),ALLOCATABLE::ifront0
   
   ALLOCATE(ifront0(nft+1))
   
   ifront0(1:nft) = ifront(1:nft)
   
   ifront = 0
   
   nnd = 0
   
   DO i = 1,nft
   
   jd = ifront0(i)
   IF(imark(jd).EQ.2) THEN
   nnd = nnd+1
   ifront(nnd) = jd
   ENDIF
   
   ENDDO
   
   nft = nnd
   
   DEALLOCATE(ifront0)
   
   RETURN
   
    END SUBROUTINE
    
! = = = = = = = = = = = = = = = = = = = = = = = = = = = !