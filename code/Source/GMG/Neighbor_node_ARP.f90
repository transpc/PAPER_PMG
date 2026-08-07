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

integer(4),dimension(:,:),allocatable::rnbcnt,snbcnt,imark
integer(4),dimension(:,:,:),allocatable::nbrecv,nbsend

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
! write(999,*)'FVM-1,nn=',nn
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
! write(999,*)'FVM-2,nn=',nn
allocate(rnbcnt(np,np),nbrecv(np,np,nn))

!  
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
!
! test
! write(999,*)'FVM-3,nn=',nn
!---------------------------
!%copy recv to send
allocate(snbcnt(np,np),nbsend(np,np,nn))
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
!
! test
! write(999,*)'FVM-4,nn=',nn
!----------------------------------------------
! ri and si !
do prc=1,np
   ri(prc,1)=1
   si(prc,1)=1
enddo


do prc=1,np
   do jp=1,nnbdom(prc)
   
      inb = nbdom(prc,jp)
	  nnb = snbcnt(prc,inb)
	  nk = si(prc,jp)
	  
      si(prc,jp+1)=nk + nnb             ! snbcnt(prc,nbdom(prc,jp))
	  
      do k=1, nnb     !snbcnt(prc,nbdom(prc,jp))
         nd=nbsend(prc,inb,k)
         sint(prc,nk-1+k)=nd
! test
         IF(nk-1+k.GT.nn) write(999,*) 'error-1',nk-1+k,nn
      enddo
   enddo
enddo
!write(*,*)'nintf=',sort(1)
! test
! write(999,*)'FVM-5,nnb=',nnb
!----------------------------

do prc=1,np
   do jp=1,nnbdom(prc)
   
      inb = nbdom(prc,jp)
	  nnb = rnbcnt(prc,inb)
	  nk = ri(prc,jp)
	  
      ri(prc,jp+1)= nk + nnb
	  
      do k=1,nnb                              !rnbcnt(prc,nbdom(prc,jp))
         nd=nbrecv(prc,inb,k)
         rint(prc,nk - 1+ k) = nd
! test
         IF(nk-1+k.GT.nn) write(999,*) 'error-2',nk-1+k,nn
      enddo
   enddo
enddo
!write(*,*)'nintf=',sort(1)
! test
! write(999,*)'FVM-6,nnb=',nnb
! - - - - - - 
     
deallocate(rnbcnt,nbrecv,imark)
deallocate(snbcnt,nbsend)
!
RETURN

    ENDSUBROUTINE

    
! = = = = = = = = = = = = = = = = = = = = = 

! - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - !
SUBROUTINE Ext_nodes_R(np,nn,nnode,nnode1,nnzr,cnode,icoarse,iar,jar,jwk,cext,lcnode3)

!
IMPLICIT NONE

! input
INTEGER (4) np,nnode,nnode1,nnzr,nn
INTEGER (4) cnode(nnode),icoarse(nnode)
INTEGER (4) iar(nnode1+1),jar(nnzr),jwk(np,nnode)
! output
INTEGER (4) cext(np),lcnode3(np,nn)
! temp
INTEGER (4) ip,jd,I,J,i1,i2,id

! 
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
   IF(jwk(ip,id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(ip,id)=1

   ENDIF
   ENDIF
   ENDDO

   Enddo
enddo
! 

RETURN
END

! - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - - - - - !
SUBROUTINE Ext_nodes_P(np,nn,nnode,nnode0,nnzi0,cnode0,cnode,iai0,jai0,jwk,cext,lcnode3)

!
IMPLICIT NONE

! input
INTEGER (4) np,nnode,nnode0,nnzi0,nn
INTEGER (4) cnode0(nnode0),cnode(nnode)
INTEGER (4) iai0(nnode0+1),jai0(nnzi0),jwk(np,nnode)
! output
INTEGER (4) cext(np),lcnode3(np,nn)
! temp
INTEGER (4) ip,jd,I,J,i1,i2,id

! 
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
   IF(jwk(ip,id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(ip,id)=1

   ENDIF
   ENDIF
   ENDDO

   Enddo
enddo
! 

RETURN
END
