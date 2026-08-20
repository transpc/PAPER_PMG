subroutine Domain_infor_FVM_fine(np,nf_max,nnode,nnodet,num_neigh,e_neigh,celem,lnum,lcelem,nbdom,nnbdom,    &
           cext,cinter,cintf,iperm,jperm,ri,si,rint,sint,nnodegl,                                   &
           nnode1,nnzi,iar,jar,icoarse)
    
	USE MD_MPI_ARP, ONLY: inbdomA,nnbdomA,riA,siA,rintA,sintA,                                 &
	                      inbdomR,nnbdomR,riR,siR,rintR,sintR, cext_tmp
    
implicit none

! inlet: 
INTEGER np,nf_max,nnode,nnodet
INTEGER e_neigh(nf_max,nnode),celem(nnode),num_neigh(nnode)
INTEGER nnode1,nnzi
INTEGER iar(nnode1+1),jar(nnzi),icoarse(nnode)
! out
INTEGER lnum(np),nbdom(np,np),nnbdom(np),cext(np),cinter(np),cintf(np),lcelem(np,nnodet)
INTEGER iperm(np,nnode),jperm(np,nnodet),ri(np,np),si(np,np),rint(np,nnodet),sint(np,nnodet)
INTEGER nnodegl(np)

! temp 
integer i,j,k,idom,nd,ie,ne,nn,proc,prc,cnt,ip,jp,id,jd,neigh,nk,n,nvpe,i1,i2,next_m
INTEGER(4)::alstatus
integer color,col1,col2,col3,col4,index,sumc,col(nf_max)
integer,dimension(:),allocatable::sort
integer,dimension(:,:),allocatable::index_elem,lcnode3,rnbcnt,snbcnt
integer,dimension(:,:,:),allocatable::nbrecv,nbsend
integer(4),dimension(:,:),allocatable::jwk
INTEGER(4) imark(np,np)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

allocate(index_elem(np,nnodet),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-element'
         STOP
     ENDIF

lnum=0
nbdom=0
nnbdom=0
!----------------------------------------------------------------------
!2-%Find local element 
do proc=1,np
   do ie=1,nnode
!      
       IF(celem(ie).NE.proc) CYCLE
       
       lnum(proc)=lnum(proc)+1
       k = lnum(proc)
       
       lcelem(proc,k)=ie
        
      nvpe = num_neigh(ie)
      col(1:nvpe) = celem(e_neigh(1:nvpe,ie))
      sumc = 0
!      
      do i=1,nvpe
          sumc = sumc +abs(col(i)-proc)
!          
      enddo

      if(sumc.eq.0)then
         index_elem(proc,k)=0
         ELSE
         index_elem(proc,k)=1 !interface element
      endif
   enddo
enddo
!----------------------------------------------------------------------
Allocate(lcnode3(np,nnodet),stat=alstatus)
Allocate(jwk(np,nnode),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-jwk'
         STOP
     ENDIF

!3-%find cinter,cintf,cext& lcnode1,2,3
cinter=0
cintf=0
cext=0
!
lcnode3=0
jwk=0

do ip=1,np
   do ie=1,lnum(ip)
      if(index_elem(ip,ie)==1)then
         cintf(ip)=cintf(ip)+1    
         ne=lcelem(ip,ie)
         jwk(ip,ne)=1    
         
         nvpe = num_neigh(ne)
         do id=1,nvpe
            jd=e_neigh(id,ne)
            
            if(jwk(ip,jd)==0)then 
               if(celem(jd).NE.ip)then
                  cext(ip)=cext(ip)+1
                  lcnode3(ip,cext(ip))=jd
                  jwk(ip,jd)=1
               endif
            endif ! assign old node
         enddo
      endif ! end of interface element
   enddo
enddo
! - - - - - - - - - - - -----
! neighbor nodes for matrix A
! - - - - - - - - - - - - - -
      nn=0
      DO i=1,np
         if(cext(i).gt.nn) nn=cext(i)
      ENDDO

      nn = 2*nn 
      nn = MAX(nn,20)
      nn = MAX(nn,np)
      
      ALLOCATE(inbdomA(np,np),nnbdomA(np))
      ALLOCATE(riA(np,np),siA(np,np),rintA(np,nn),sintA(np,nn))
	
      CALL Neighbor_node_ARP(np,nnode,nnodet,nn,celem,cext,lcnode3,inbdomA,nnbdomA,riA,siA,rintA,sintA)

! - - - - - - - - - - - - - - - - - - - - - - - - - !
! add new cext from R(Ij):for all nodes

do ip=1,np
   do ie=1,lnum(ip)
       jd = lcelem(ip,ie)

   I = icoarse(jd)
   IF(I.EQ.0) CYCLE
   i1 = iar(I)
   i2 = iar(I+1)-1
   DO J = i1,i2
   id = jar(J)
   IF(celem(id).NE.ip) THEN
   IF(jwk(ip,id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(ip,id)=1

!     IF(jwk(jd,ip).EQ.0) THEN
!      cintf(ip)=cintf(ip)+1

!      jwk(jd,ip) = 1

!     ENDIF

!notes-> update for intf
    jp= celem(id)
   IF(jwk(jp,id)==1) CYCLE
      cintf(jp)=cintf(jp)+1
!   
      jwk(jp,id)=1      
!     

   ENDIF
   ENDIF
   ENDDO

   Enddo
enddo

!-----------------------------------------------
!%mapping: iperm::global->local
!%mapping: jperm::local->global
iperm=0
jperm=0
do ip=1,np
   do ie=1,lnum(ip)
      if(index_elem(ip,ie)==0)then
         ne=lcelem(ip,ie)
         IF(jwk(ip,ne).EQ.1) CYCLE     ! notes for R
               cinter(ip)=cinter(ip)+1
!           
               iperm(ip,ne)=cinter(ip)
               jperm(ip,cinter(ip))=ne
!           
      endif ! end of internal element
   enddo
enddo
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
         if(celem(lcnode3(ip,i))==jp) THEN
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
      
      next_m = 2*next_m 
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
         if(celem(jd)==neigh) then
            rnbcnt(prc,neigh)=rnbcnt(prc,neigh)+1
            cnt=rnbcnt(prc,neigh)
            nbrecv(prc,neigh,cnt)=jd
         endif
      enddo
   enddo
enddo
!----------------------------------------------------------------------
!%copy recv to send
allocate(snbcnt(np,np),nbsend(np,np,nn),stat=alstatus)
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
write(999,*)'neq=',sort(1)
! 
! NEW: added more element and nodes for Garlekin F

!1-neiboring Element of each node:     
!------------
!2-index element:
!      index_elem = 0
      jwk = 0
      
      DO proc=1,np
          
          DO i=1,sort(proc)
              id = jperm(proc,i)
              jwk(proc,id)=1
          ENDDO
          
      ENDDO
      
! 3-
      DO proc=1,np
         nn=cext(proc)
         ip=lnum(proc) 
         DO j=1,nn
            ie =lcnode3(proc,j)
! 
                ip = ip + 1
                lcelem(proc,ip)=ie
! 
                nvpe = num_neigh(ie)
                DO i=1,nvpe
                    id=e_neigh(i,ie)
                    IF(jwk(proc,id).EQ.1) CYCLE
                      jwk(proc,id)=1
                      sort(proc)=sort(proc)+1
                      nk=sort(proc) !!temporary
                      iperm(proc,id)=nk
                      jperm(proc,nk)=id
                ENDDO
!           
         ENDDO
         lnum(proc) = ip
         
      ENDDO
      
      nnodegl = sort
     write(999,*)'nnodegl=',sort(1)   
     
i = MAXVAL(sort)
IF(i.GT.nnodet) THEN
WRITE(*,*)'PMG error: nnodet small',nnodet,i
STOP
ENDIF

! NEW for neighbor nodes of matrix R:
! - - - - - - - - - - - - - - - - - -
     DEALLOCATE(lcnode3)
	 ALLOCATE(lcnode3(np,next_m),cext_tmp(np))
	 lcnode3 = 0
	 cext_tmp = 0
!
     CALL Ext_nodes_R(np,next_m,nnode,nnode1,nnzi,celem,icoarse,iar,jar,jwk,cext_tmp,lcnode3)
     
	 ALLOCATE(inbdomR(np,np),nnbdomR(np),riR(np,np),siR(np,np),rintR(np,next_m),sintR(np,next_m))
     CALL Neighbor_node_ARP(np,nnode,next_m,next_m,celem,cext_tmp,lcnode3,inbdomR,nnbdomR,riR,siR,rintR,sintR)	 
	 
	 DEALLOCATE(cext_tmp)
! - - - - - - 
! - - - - - - 
deallocate(index_elem,sort)
deallocate(rnbcnt,nbrecv)
deallocate(snbcnt,nbsend)
deallocate(lcnode3)
deallocate(jwk)

     RETURN

END
