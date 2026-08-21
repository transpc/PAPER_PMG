SUBROUTINE Domain_infor_FVM_coarse(ilv,np,mxnbne,nnode,nnodet,nnode1,nnzi,iar,jar,icoarse,nnei,inei,cnode,          &
           nbdom,nnbdom,cext,cinter,cintf,jperm,ri,si,rint,sint,nnodegl,                 &
		   nnode0,cnode0,nnzi0,iai0,jai0)
		   
	USE MD_MPI_ARP, ONLY: inbdomA,nnbdomA,riA,siA,rintA,sintA,                                 &
	                      inbdomR,nnbdomR,riR,siR,rintR,sintR,                                 &
	                      inbdomP,nnbdomP,riP,siP,rintP,sintP,cext_tmp						  
						  
implicit none

! inlet: 
INTEGER np,nnode,mxnbne,nnodet,ilv
INTEGER nnode1,nnzi
INTEGER iar(nnode1+1),jar(nnzi),icoarse(nnode)
INTEGER nnei(nnode),inei(mxnbne,nnode),cnode(nnode)
INTEGER nnode0,nnzi0
INTEGER cnode0(nnode0),iai0(nnode0+1),jai0(nnzi0)
! out
INTEGER nbdom(np,np),nnbdom(np),cext(np),cinter(np),cintf(np)
INTEGER jperm(np,nnodet),ri(np,np),si(np,np),rint(np,nnodet),sint(np,nnodet)
INTEGER nnodegl(np)

! temp 
integer i,j,k,idom,nd,ie,ne,nn,proc,prc,cnt,ip,jp,id,jd,neigh,nnd,nk,i1,i2,next_m
integer color,col1,col2,col3,col4,index,sumc
INTEGER(4)::alstatus
integer,dimension(:),allocatable::sort
integer,dimension(:),allocatable::index_node,node_intf   ! 소유 랭크 기준 셀별 플래그
integer,dimension(:),allocatable::jwk        ! 랭크별 고스트 마커 — 1D 로 재사용 (구간마다 복원/원복)
integer,dimension(:,:),allocatable::lcnode3,rnbcnt
integer,dimension(:,:,:),allocatable::nbrecv
INTEGER(4),DIMENSION(:,:), ALLOCATABLE:: lcelem
INTEGER(4),DIMENSION(:), ALLOCATABLE:: lnum
INTEGER(4) imark(np,np)

! - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - !

nbdom=0
nnbdom=0
!----------------------------------------------------------------------
!1-%Find local interface node (inside interface)
allocate(index_node(nnode),node_intf(nnode),stat=alstatus)
ALLOCATE(lnum(np),lcelem(np,nnodet),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MG-MPI-index-node'
         STOP
     ENDIF
     
index_node = 0
node_intf = 0
lnum = 0

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
   index_node(i) = 1

   EXIT
   ENDIF  
   ENDDO
   ENDDO
   lnum(proc) = nd
   
ENDDO
!----------------------------------------------------------------------
Allocate(jwk(nnode),lcnode3(np,nnodet),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MG-MPI-jwk'
         STOP
     ENDIF

!3-%find cinter,cintf,cext& lcnode1,2,3
cinter=0
cintf=0
cext=0
jwk=0

do ip=1,np
    DO i=1,lnum(ip)
    jd=lcelem(ip,i)
   IF(index_node(jd).EQ.0) CYCLE

   nnd = nnei(jd)
   DO j = 1,nnd
   id = inei(j,jd)
   IF(cnode(id).NE.ip) THEN
   IF(jwk(id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(id)=1
! new
   jp = cnode(id)
   IF(node_intf(id).EQ.0) THEN
          cintf(jp)=cintf(jp)+1
          node_intf(id) = 1
   ENDIF
   ENDIF
!   
   ENDIF
   
   ENDDO
   
   ENDDO
!  이 ip 가 표시한 고스트 마커만 원복 (다음 ip 재사용)
   DO idom=1,cext(ip)
      jwk(lcnode3(ip,idom))=0
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
! add new cext from R(Ij):
   
do ip=1,np
!  직전 단계에서 이 ip 가 만든 고스트 마커 복원
   DO idom=1,cext(ip)
      jwk(lcnode3(ip,idom))=1
   ENDDO
    DO ie = 1,lnum(ip)
    jd=lcelem(ip,ie)

   I = icoarse(jd)
   IF(I.EQ.0) CYCLE
   i1 = iar(I)
   i2 = iar(I+1)-1
   DO J = i1,i2
   id = jar(J)
   IF(cnode(id).NE.ip) THEN
   IF(jwk(id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(id)=1
!
    jp= cnode(id)

   IF(node_intf(id)==1) CYCLE
      cintf(jp)=cintf(jp)+1  
      node_intf(id)=1   

   ENDIF
   ENDIF
   ENDDO
!
   ENDDO
!  이 ip 가 표시한 고스트 마커만 원복 (다음 ip 재사용)
   DO idom=1,cext(ip)
      jwk(lcnode3(ip,idom))=0
   ENDDO
enddo	

! add new cext from P(iJ):
   
do ip=1,np
!  직전 단계에서 이 ip 가 만든 고스트 마커 복원
   DO idom=1,cext(ip)
      jwk(lcnode3(ip,idom))=1
   ENDDO
   do jd=1,nnode0
   IF(cnode0(jd).NE.ip) CYCLE

   I = jd
   
   i1 = iai0(I)
   i2 = iai0(I+1)-1
   DO J = i1,i2
   id = jai0(J)
   IF(cnode(id).NE.ip) THEN
   IF(jwk(id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(id)=1
      
!!notes
    jp= cnode(id)
   IF(node_intf(id)==1) CYCLE
      cintf(jp)=cintf(jp)+1
!   
      node_intf(id)=1   
!     

   ENDIF
   ENDIF
   ENDDO
!
   ENDDO
!  이 ip 가 표시한 고스트 마커만 원복 (다음 ip 재사용)
   DO idom=1,cext(ip)
      jwk(lcnode3(ip,idom))=0
   ENDDO
enddo  

!-----------------------------------------------
!%mapping: jperm::local->global
! jperm 전면 초기화 없음 — 유효 길이는 cinter/sort 가 한정하고, 판독은 전부 그 범위 안
! (포이즌 값 주입으로 np=1~64 전 게이트·합성 검증). 미기록 페이지는 실체화되지 않음

do ip=1,np
    DO i=1,lnum(ip)
    jd=lcelem(ip,i)
   IF(cnode(jd).NE.ip) CYCLE
!
    IF(node_intf(jd).EQ.1) CYCLE    ! new
        cinter(ip)=cinter(ip)+1
!       
        jperm(ip,cinter(ip))=jd
	ENDDO
ENDDO

!-------------------------------------------------NEW-
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
      
allocate(rnbcnt(np,np),nbrecv(next_m,np,np),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI2-index-nbrecv'
         STOP
     ENDIF
     
rnbcnt=0
do prc=1,np
   do ip=1,nnbdom(prc)
      neigh=nbdom(prc,ip)
      do id=1,cext(prc)
         jd=lcnode3(prc,id)
         if(cnode(jd)==neigh) then
            rnbcnt(prc,neigh)=rnbcnt(prc,neigh)+1
            cnt=rnbcnt(prc,neigh)
            nbrecv(cnt,prc,neigh)=jd
         endif
      enddo
   enddo
enddo
!----------------------------------------------------------------------
! send 목록 = recv 목록의 전치 — 별도 배열 없이 rnbcnt/nbrecv 를 (이웃,자기) 순서로 직접 참조
!----------------------------------------------
! ri and si !
do prc=1,np
   ri(prc,1)=1
   si(prc,1)=1
enddo
! send/recv 두 패스를 prc 단위로 합침 — prc 별 jperm 추가 순서(送 후 受) 불변
do prc=1,np
   nk = sort(prc)                      ! 이 prc 표시 시작 위치 (원복 기준)
   do jp=1,nnbdom(prc)
      si(prc,jp+1)=si(prc,jp)+rnbcnt(nbdom(prc,jp),prc)
      do k=1,rnbcnt(nbdom(prc,jp),prc)
         nd=nbrecv(k,nbdom(prc,jp),prc)
         if(jwk(nd)==0) then
            sort(prc)=sort(prc)+1
            nn=sort(prc) !!temporary
            jperm(prc,nn)=nd
            jwk(nd)=1
         endif
         sint(prc,si(prc,jp)-1+k)=nd
      enddo
   enddo
!---------------------------------------------------------
   do jp=1,nnbdom(prc)
      ri(prc,jp+1)=ri(prc,jp)+rnbcnt(prc,nbdom(prc,jp))
      do k=1,rnbcnt(prc,nbdom(prc,jp))
         nd=nbrecv(k,prc,nbdom(prc,jp))
         if(jwk(nd).eq.0) then
            sort(prc)=sort(prc)+1
            nn=sort(prc) !!temporary
            jperm(prc,nn)=nd
            jwk(nd)=1
         endif
         rint(prc,ri(prc,jp)-1+k)=nd
      enddo
   enddo
!  원복
   do k=nk+1,sort(prc)
      jwk(jperm(prc,k))=0
   enddo
enddo
!--------------------------------------------
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
! 표시 패스와 확장 패스를 proc 단위로 합침 (proc 간 의존 없음)
      DO proc=1,np
          
          DO i=1,sort(proc)
              id = jperm(proc,i)
              jwk(id)=1
          ENDDO
          
! 3-add: 
         nn=cext(proc)
         DO j=1,nn
            jd =lcnode3(proc,j)
            
            nnd = nnei(jd)
            DO i = 1,nnd
            id = inei(i,jd)
            
                IF(jwk(id).EQ.1) CYCLE
                jwk(id)=1
                sort(proc)=sort(proc)+1
                nk=sort(proc) !!temporary
                jperm(proc,nk)=id
            ENDDO
         ENDDO
!  원복
         DO i=1,sort(proc)
             jwk(jperm(proc,i))=0
         ENDDO
         
      ENDDO
      
      nnodegl = sort
! - - - - - - 
i = MAXVAL(sort)
IF(i.GT.nnodet) THEN
WRITE(*,*)'error,nnodet small',nnodet,i
ENDIF
! - - - - - - - - - - - - - - - - - - 
! NEW for neighbor nodes of matrix R:
! - - - - - - - - - - - - - - - - - -
     DEALLOCATE(lcnode3)
	 ALLOCATE(lcnode3(np,next_m),cext_tmp(np))
	 cext_tmp = 0
!
     CALL Ext_nodes_R(np,next_m,nnode,nnode1,nnzi,cnode,icoarse,iar,jar,cext_tmp,lcnode3)
     
	 ALLOCATE(inbdomR(np,np),nnbdomR(np),riR(np,np),siR(np,np),rintR(np,next_m),sintR(np,next_m))
     CALL Neighbor_node_ARP(np,nnode,next_m,next_m,cnode,cext_tmp,lcnode3,inbdomR,nnbdomR,riR,siR,rintR,sintR)	 
	 
! - - - - - - - - - - - - - - - - - - 
! NEW for neighbor nodes of matrix P:
! - - - - - - - - - - - - - - - - - -
	 cext_tmp = 0
!
     CALL Ext_nodes_P(np,next_m,nnode,nnode0,nnzi0,cnode0,cnode,iai0,jai0,cext_tmp,lcnode3)
	 
	 ALLOCATE(inbdomP(np,np),nnbdomP(np),riP(np,np),siP(np,np),rintP(np,next_m),sintP(np,next_m))
     CALL Neighbor_node_ARP(np,nnode,next_m,next_m,cnode,cext_tmp,lcnode3,inbdomP,nnbdomP,riP,siP,rintP,sintP)	
	 
	 DEALLOCATE(cext_tmp)
! - - - - - - 

deallocate(index_node,sort,node_intf)
deallocate(jwk,lcnode3)
deallocate(rnbcnt,nbrecv)
DEALLOCATE(lnum,lcelem)

RETURN

ENDSUBROUTINE
