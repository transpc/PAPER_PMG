SUBROUTINE Domain_infor_FVM_coarsest(ilv,np,mxnbne,nnode,nnodet,nnei,inei,cnode,                            &
           nbdom,nnbdom,cext,cinter,cintf,iperm,jperm,ri,si,rint,sint,nnodegl,                 &
		   nnode0,cnode0,nnzi0,iai0,jai0)
		   
	USE MD_MPI_ARP, ONLY: inbdomA,nnbdomA,riA,siA,rintA,sintA,                                 &
	                      inbdomP,nnbdomP,riP,siP,rintP,sintP,cext_tmp

implicit none

! inlet: 
INTEGER np,nnode,mxnbne,nnzi,nnodet,ilv
INTEGER nnei(nnode),inei(mxnbne,nnode),cnode(nnode)
INTEGER nnode0,nnzi0
INTEGER cnode0(nnode0),iai0(nnode0+1),jai0(nnzi0)
! out
INTEGER nbdom(np,np),nnbdom(np),cext(np),cinter(np),cintf(np)
INTEGER iperm(np,nnode),jperm(np,nnodet),ri(np,np),si(np,np),rint(np,nnodet),sint(np,nnodet)
INTEGER nnodegl(np)

! temp 

integer i,j,k,idom,nd,ie,ne,nn,proc,prc,cnt,ip,jp,id,jd,neigh,nnd,nk,i1,i2,next_m
integer color,col1,col2,col3,col4,index,sumc
INTEGER(4)::alstatus
integer,dimension(:),allocatable::sort
integer,dimension(:,:),allocatable::index_node,jwk,lcnode3,rnbcnt,snbcnt
integer,dimension(:,:,:),allocatable::nbrecv,nbsend
integer,dimension(:,:),allocatable::node_intf
INTEGER(4),DIMENSION(:,:), ALLOCATABLE:: lcelem
INTEGER(4),DIMENSION(:), ALLOCATABLE:: lnum
INTEGER(4) imark(np,np)

! - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - !

nbdom=0
nnbdom=0 
!----------------------------------------------------------------------
!1-%Find local interface node (inside interface)
allocate(index_node(np,nnode),node_intf(np,nnode),stat=alstatus)
ALLOCATE(lnum(np),lcelem(np,nnodet),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MG-MPI-index-node'
         STOP
     ENDIF
     
index_node = 0
node_intf = 0
lnum = 0
lcelem = 0

DO proc = 1,np
    nd = 0
   DO i = 1,nnode
   IF(cnode(i).NE.proc) CYCLE
    nd = nd+1
    lcelem(proc,nd) = i
!
   nnd = nnei(i)
   DO j = 1,nnd
   id = inei(j,i)
   IF(cnode(id).NE.proc) THEN
   index_node(proc,i) = 1

   EXIT
   ENDIF  
   ENDDO
   ENDDO
   lnum(proc) = nd
   
ENDDO

!----------------------------------------------------------------------
Allocate(jwk(np,nnode),lcnode3(np,nnodet),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MG-MPI-jwk'
         STOP
     ENDIF

!3-%find cinter,cintf,cext& lcnode1,2,3
cinter=0
cintf=0
cext=0
lcnode3=0
jwk=0

do ip=1,np
    DO i=1,lnum(ip)
    jd= lcelem(ip,i)
   IF(index_node(ip,jd).EQ.0) CYCLE

   nnd = nnei(jd)
   DO j = 1,nnd
   id = inei(j,jd)
   IF(cnode(id).NE.ip) THEN
   IF(jwk(ip,id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(ip,id)=1
! new
   jp = cnode(id)
   IF(node_intf(jp,id).EQ.0) THEN
          cintf(jp)=cintf(jp)+1
          node_intf(jp,id) = 1
   ENDIF
   ENDIF
!   
   ENDIF
   
   ENDDO
   
   ENDDO
enddo	  

! - - - - - - - - - - - -----
! neighbor nodes for matrix A
! - - - - - - - - - - - - - -
      nn=0
      DO i=1,np
         if(cext(i).gt.nn) nn=cext(i)
      ENDDO

      nn = ilv*nn
      nn = MAX(nn,20)
      nn = MAX(nn,np)
      
      ALLOCATE(inbdomA(np,np),nnbdomA(np))
      ALLOCATE(riA(np,np),siA(np,np),rintA(np,nn),sintA(np,nn))
	
      CALL Neighbor_node_ARP(np,nnode,nnodet,nn,cnode,cext,lcnode3,inbdomA,nnbdomA,riA,siA,rintA,sintA)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

! add new cext from P(iJ):
   
do ip=1,np
   do jd=1,nnode0
   IF(cnode0(jd).NE.ip) CYCLE

   I = jd
   
   i1 = iai0(I)
   i2 = iai0(I+1)-1
   DO J = i1,i2
   id = jai0(J)
   IF(cnode(id).NE.ip) THEN
   IF(jwk(ip,id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(ip,id)=1
      
!!notes
    jp= cnode(id)
   IF(node_intf(jp,id)==1) CYCLE
      cintf(jp)=cintf(jp)+1
!   
      node_intf(jp,id)=1   
!     

   ENDIF
   ENDIF
   ENDDO
!
   ENDDO
enddo

!-----------------------------------------------
!%mapping: iperm::global->local
!%mapping: jperm::local->global
iperm=0
jperm=0

do ip=1,np
    DO i=1,lnum(ip)
    jd=lcelem(ip,i)
!
    IF(node_intf(ip,jd).EQ.1) CYCLE    ! new
        cinter(ip)=cinter(ip)+1
!       
        iperm(ip,jd)=cinter(ip)
        jperm(ip,cinter(ip))=jd
	ENDDO
ENDDO

!----------------------------------------------------------------------
allocate(sort(np))
sort=cinter

! 4-nbdom

imark = 0

do ip=1,np
   do jp=1,np
       IF(jp==ip) CYCLE
       IF(imark(ip,jp)==1) CYCLE
      index=0
      do i=1,cext(ip)
         if(cnode(lcnode3(ip,i))==jp) THEN
             index=1
             EXIT 
         ENDIF
         
      enddo
      if(index==1)then

         nnbdom(ip)=nnbdom(ip)+1
         nbdom(ip,nnbdom(ip))=jp
         imark(ip,jp) = 1
! new          
          IF(imark(jp,ip)==0) THEN
         nnbdom(jp)=nnbdom(jp)+1
         nbdom(jp,nnbdom(jp))=ip
         imark(jp,ip) = 1   
         ENDIF
              
      endif
   enddo
enddo
!----------------------------------------------------------------------
!%cext=total num of exteria nodes...
!%array for recv&send variables  in "SERIAL"(Global mesh)

      next_m=0
      DO i=1,np
         if(cext(i).gt.next_m) next_m=cext(i)
      ENDDO
      
      next_m = ilv*next_m 
      next_m = MAX(next_m,20)
      next_m = MAX(next_m,np)
      
allocate(rnbcnt(np,np),nbrecv(np,np,next_m),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-nbrecv'
         STOP
     ENDIF
     
rnbcnt=0
nbrecv=0
do prc=1,np
   do ip=1,nnbdom(prc)
      neigh=nbdom(prc,ip)
      do id=1,cext(prc)
         jd=lcnode3(prc,id)
         if(cnode(jd)==neigh) then
            rnbcnt(prc,neigh)=rnbcnt(prc,neigh)+1
            cnt=rnbcnt(prc,neigh)
            nbrecv(prc,neigh,cnt)=jd
         endif
      enddo
   enddo
enddo
!----------------------------------------------------------------------
!%copy recv to send
allocate(snbcnt(np,np),nbsend(np,np,next_m),stat=alstatus)
     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-nbsend'
         STOP
     ENDIF
     
do jp=1,np
   do ip=1,np
      snbcnt(jp,ip)=rnbcnt(ip,jp)
   enddo
enddo
do jp=1,np
   do ip=1,np
      cnt=rnbcnt(ip,jp)
      if(cnt>0)then
         do j=1,cnt
            nbsend(jp,ip,j)=nbrecv(ip,jp,j)
          enddo
      endif
   enddo
enddo
!----------------------------------------------
! ri and si !
do prc=1,np
   ri(prc,1)=1
   si(prc,1)=1
enddo
jwk=0

do prc=1,np
   do jp=1,nnbdom(prc)
      si(prc,jp+1)=si(prc,jp)+snbcnt(prc,nbdom(prc,jp))
      do k=1,snbcnt(prc,nbdom(prc,jp))
         nd=nbsend(prc,nbdom(prc,jp),k)
         if(jwk(prc,nd)==0) then
            sort(prc)=sort(prc)+1
            nn=sort(prc) !!temporary
            iperm(prc,nd)=nn
            jperm(prc,nn)=nd
            jwk(prc,nd)=1
         endif
         sint(prc,si(prc,jp)-1+k)=nd
      enddo
   enddo
enddo
!write(*,*)'nintf=',sort(1)
!---------------------------------------------------------
do prc=1,np
   do jp=1,nnbdom(prc)
      ri(prc,jp+1)=ri(prc,jp)+rnbcnt(prc,nbdom(prc,jp))
      do k=1,rnbcnt(prc,nbdom(prc,jp))
         nd=nbrecv(prc,nbdom(prc,jp),k)
         if(jwk(prc,nd).eq.0) then
            sort(prc)=sort(prc)+1
            nn=sort(prc) !!temporary
            iperm(prc,nd)=nn
            jperm(prc,nn)=nd
            jwk(prc,nd)=1
         endif
         rint(prc,ri(prc,jp)-1+k)=nd
      enddo
   enddo
enddo
!write(*,*)'neq=',sort(1)
!--------------------------------------------
! 
i = MAXVAL(sort)
IF(i.GT.nnodet) THEN
WRITE(*,*)'error,nnodet small',nnodet,i
ENDIF

! cheking: 
  Do ip = 1, np
      i = cinter(ip)+cintf(ip)
      IF(lnum(ip).NE.i) THEN
          write(*,*)'error in lnum',lnum(ip),i
          STOP
      ENDIF
  ENDDO
  
! NEW: added more nodes for Garlekin F
!------------
!1-index nodes:
      jwk = 0
      
      DO proc=1,np
          
          DO i=1,sort(proc)
              id = jperm(proc,i)
              jwk(proc,id)=1
          ENDDO
          
      ENDDO
      
! 3-add: 
      DO proc=1,np
          
         nn=cext(proc)
         DO j=1,nn
            jd =lcnode3(proc,j)
            
            nnd = nnei(jd)
            DO i = 1,nnd
            id = inei(i,jd)
            
                IF(jwk(proc,id).EQ.1) CYCLE
                jwk(proc,id)=1
                sort(proc)=sort(proc)+1
                nk=sort(proc) !!temporary
                iperm(proc,id)=nk
                jperm(proc,nk)=id
            ENDDO
         ENDDO
         
      ENDDO
      
      nnodegl = sort
!     write(*,*)'nnodegl=',sort(1)      
! - - - - - - 
i = MAXVAL(sort)
IF(i.GT.nnodet) THEN
WRITE(*,*)'error,nnodet small',nnodet,i
ENDIF
! - - - - - - - - - - - - - - - - - - 
! NEW for neighbor nodes of matrix P:
! - - - - - - - - - - - - - - - - - -
     DEALLOCATE(lcnode3)
	 ALLOCATE(lcnode3(np,next_m),cext_tmp(np))
	 lcnode3 = 0
	 cext_tmp = 0
!
     CALL Ext_nodes_P(np,next_m,nnode,nnode0,nnzi0,cnode0,cnode,iai0,jai0,jwk,cext_tmp,lcnode3)
     
	 ALLOCATE(inbdomP(np,np),nnbdomP(np),riP(np,np),siP(np,np),rintP(np,next_m),sintP(np,next_m))
     CALL Neighbor_node_ARP(np,nnode,next_m,next_m,cnode,cext_tmp,lcnode3,inbdomP,nnbdomP,riP,siP,rintP,sintP)	
	 
	 DEALLOCATE(cext_tmp)
! - - - - - - - - - - - - - - - 
deallocate(index_node,sort,node_intf)
deallocate(jwk,lcnode3)
deallocate(rnbcnt,nbrecv)
deallocate(snbcnt,nbsend)
DEALLOCATE(lnum,lcelem)

RETURN

ENDSUBROUTINE
