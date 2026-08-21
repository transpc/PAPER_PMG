subroutine Neighbor_node_ARP(np,nnode,nnodet,nn,cnode,cext,lcnode3,nbdom,nnbdom,ri,si,rint,sint)
implicit none

! inlet: 
INTEGER(4) np,nnode,nn,nnodet
INTEGER(4) cnode(nnode),lcnode3(np,nnodet),cext(np)
! out
INTEGER(4) nbdom(np,np),nnbdom(np)
INTEGER(4) ri(np,np),si(np,np),rint(np,nn),sint(np,nn)

! temp 
integer(4) i,j,k,nd,nk,inb,nnb,prc,cnt,ip,jp,id,jd,neigh,index

integer(4),dimension(:,:),allocatable::rnbcnt,imark
integer(4),dimension(:,:,:),allocatable::nbrecv

! - - - - - - - - - - - - - - - - - - - - - - - - - - -
nbdom=0
nnbdom=0
ri = 0
si = 0
rint = 0
sint = 0

Allocate(imark(np,np))
imark = 0
! 

! test
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
! - - - - - - - - - - 
! test
allocate(rnbcnt(np,np),nbrecv(nn,np,np))

!  
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
!
! test
!---------------------------
! send 목록 = recv 목록의 전치 — 별도 배열 없이 rnbcnt/nbrecv 를 (이웃,자기) 순서로 직접 참조
!
! test
!----------------------------------------------
! ri and si !
do prc=1,np
   ri(prc,1)=1
   si(prc,1)=1
enddo


do prc=1,np
   do jp=1,nnbdom(prc)
   
      inb = nbdom(prc,jp)
	  nnb = rnbcnt(inb,prc)
	  nk = si(prc,jp)
	  
      si(prc,jp+1)=nk + nnb
	  
      do k=1, nnb
         nd=nbrecv(k,inb,prc)
         sint(prc,nk-1+k)=nd
! test
         IF(nk-1+k.GT.nn) write(*,*) 'PMG error-1 (ARP)',nk-1+k,nn
      enddo
   enddo
enddo
! test
!----------------------------

do prc=1,np
   do jp=1,nnbdom(prc)
   
      inb = nbdom(prc,jp)
	  nnb = rnbcnt(prc,inb)
	  nk = ri(prc,jp)
	  
      ri(prc,jp+1)= nk + nnb
	  
      do k=1,nnb                              !rnbcnt(prc,nbdom(prc,jp))
         nd=nbrecv(k,prc,inb)
         rint(prc,nk - 1+ k) = nd
! test
         IF(nk-1+k.GT.nn) write(*,*) 'PMG error-2 (ARP)',nk-1+k,nn
      enddo
   enddo
enddo
! test
! - - - - - - 
     
deallocate(rnbcnt,nbrecv,imark)
!
RETURN

    ENDSUBROUTINE

    
! = = = = = = = = = = = = = = = = = = = = = 

! - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - !
SUBROUTINE Ext_nodes_R(np,nn,nnode,nnode1,nnzr,cnode,icoarse,iar,jar,cext,lcnode3)

!
IMPLICIT NONE

! input
INTEGER (4) np,nnode,nnode1,nnzr,nn
INTEGER (4) cnode(nnode),icoarse(nnode)
INTEGER (4) iar(nnode1+1),jar(nnzr)
! output
INTEGER (4) cext(np),lcnode3(np,nn)
! temp
INTEGER (4) ip,jd,I,J,i1,i2,id
INTEGER (4),ALLOCATABLE:: jwk(:)      ! 랭크별 고스트 마커 — 1D 로 재사용 (ip 마다 원복)

! 
ALLOCATE(jwk(nnode))
jwk = 0
!
do ip=1,np
   do jd=1,nnode
   IF(cnode(jd).NE.ip) CYCLE

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

   ENDIF
   ENDIF
   ENDDO

   Enddo
!  다음 ip 를 위해 이 ip 가 표시한 것만 원복
   DO jd=1,cext(ip)
      jwk(lcnode3(ip,jd))=0
   ENDDO
enddo
! 
DEALLOCATE(jwk)

RETURN
END

! - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - !
SUBROUTINE Ext_nodes_P(np,nn,nnode,nnode0,nnzi0,cnode0,cnode,iai0,jai0,cext,lcnode3)

!
IMPLICIT NONE

! input
INTEGER (4) np,nnode,nnode0,nnzi0,nn
INTEGER (4) cnode0(nnode0),cnode(nnode)
INTEGER (4) iai0(nnode0+1),jai0(nnzi0)
! output
INTEGER (4) cext(np),lcnode3(np,nn)
! temp
INTEGER (4) ip,jd,I,J,i1,i2,id
INTEGER (4),ALLOCATABLE:: jwk(:)      ! 랭크별 고스트 마커 — 1D 로 재사용 (ip 마다 원복)

! 
ALLOCATE(jwk(nnode))
jwk = 0
!
do ip=1,np
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

   ENDIF
   ENDIF
   ENDDO

   Enddo
!  다음 ip 를 위해 이 ip 가 표시한 것만 원복
   DO jd=1,cext(ip)
      jwk(lcnode3(ip,jd))=0
   ENDDO
enddo
! 
DEALLOCATE(jwk)

RETURN
END
