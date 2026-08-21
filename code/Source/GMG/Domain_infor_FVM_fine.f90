subroutine Domain_infor_FVM_fine(np,nf_max,nnode,nnodet,num_neigh,e_neigh,celem,lnum,lcelem,nbdom,nnbdom,    &
           cext,cinter,cintf,jperm,ri,si,rint,sint,nnodegl,                                   &
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
INTEGER jperm(np,nnodet),ri(np,np),si(np,np),rint(np,nnodet),sint(np,nnodet)
INTEGER nnodegl(np)

! temp 
integer i,j,k,idom,nd,ie,ne,nn,proc,prc,cnt,ip,jp,id,jd,neigh,nk,n,nvpe,i1,i2,next_m
INTEGER(4)::alstatus
integer color,col1,col2,col3,col4,index,sumc,col(nf_max)
integer,dimension(:),allocatable::sort
integer,dimension(:),allocatable::index_cell          ! 셀별 인터페이스 플래그 (소유 랭크 기준)
integer,dimension(:,:),allocatable::lcnode3,rnbcnt
integer,dimension(:,:,:),allocatable::nbrecv
integer(4),dimension(:),allocatable::jwk,mark_own   ! jwk: 랭크별 고스트 마커(1D 재사용), mark_own: 소유자 인터페이스 표시
INTEGER(4) imark(np,np)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !

allocate(index_cell(nnode),stat=alstatus)

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
         index_cell(ie)=0
         ELSE
         index_cell(ie)=1 !interface element
      endif
   enddo
enddo
!----------------------------------------------------------------------
Allocate(lcnode3(np,nnodet),stat=alstatus)
Allocate(jwk(nnode),mark_own(nnode),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-jwk'
         STOP
     ENDIF

!3-%find cinter,cintf,cext& lcnode1,2,3
cinter=0
cintf=0
cext=0
!
jwk=0
mark_own=0

do ip=1,np
   do ie=1,lnum(ip)
      ne=lcelem(ip,ie)
      if(index_cell(ne)==1)then
         cintf(ip)=cintf(ip)+1    
         mark_own(ne)=1    
         
         nvpe = num_neigh(ne)
         do id=1,nvpe
            jd=e_neigh(id,ne)
            
            if(celem(jd).NE.ip)then
               if(jwk(jd)==0)then 
                  cext(ip)=cext(ip)+1
                  lcnode3(ip,cext(ip))=jd
                  jwk(jd)=1
               endif
            endif ! assign old node
         enddo
      endif ! end of interface element
   enddo
!  이 ip 가 표시한 고스트 마커만 원복 (다음 ip 재사용)
   do ie=1,cext(ip)
      jwk(lcnode3(ip,ie))=0
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
!  A 단계에서 이 ip 가 만든 고스트 마커 복원
   do ie=1,cext(ip)
      jwk(lcnode3(ip,ie))=1
   enddo
   do ie=1,lnum(ip)
       jd = lcelem(ip,ie)

   I = icoarse(jd)
   IF(I.EQ.0) CYCLE
   i1 = iar(I)
   i2 = iar(I+1)-1
   DO J = i1,i2
   id = jar(J)
   IF(celem(id).NE.ip) THEN
   IF(jwk(id).EQ.0) THEN
      cext(ip)=cext(ip)+1
      lcnode3(ip,cext(ip))=id
      jwk(id)=1

!notes-> update for intf
    jp= celem(id)
   IF(mark_own(id)==1) CYCLE
      cintf(jp)=cintf(jp)+1
!   
      mark_own(id)=1      
!     

   ENDIF
   ENDIF
   ENDDO

   Enddo
!  원복
   do ie=1,cext(ip)
      jwk(lcnode3(ip,ie))=0
   enddo
enddo

!-----------------------------------------------
!%mapping: jperm::local->global
jperm=0
do ip=1,np
   do ie=1,lnum(ip)
      ne=lcelem(ip,ie)
      if(index_cell(ne)==0)then
         IF(mark_own(ne).EQ.1) CYCLE     ! notes for R
               cinter(ip)=cinter(ip)+1
!           
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
      
allocate(rnbcnt(np,np),nbrecv(next_m,np,np),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-pre-MPI1-index-nbrecv'
         STOP
     ENDIF
     
rnbcnt=0
do prc=1,np
   do ip=1,nnbdom(prc)
      neigh=nbdom(prc,ip)
      do id=1,cext(prc)
         jd=lcnode3(prc,id)
         if(celem(jd)==neigh) then
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
write(999,*)'neq=',sort(1)
! 
! NEW: added more element and nodes for Garlekin F

!1-neiboring Element of each node:     
!------------
!2-index element:
! 표시 패스와 확장 패스를 proc 단위로 합침 (proc 간 의존 없음)
      DO proc=1,np
          
          DO i=1,sort(proc)
              id = jperm(proc,i)
              jwk(id)=1
          ENDDO
          
! 3-
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
                    IF(jwk(id).EQ.1) CYCLE
                      jwk(id)=1
                      sort(proc)=sort(proc)+1
                      nk=sort(proc) !!temporary
                      jperm(proc,nk)=id
                ENDDO
!           
         ENDDO
         lnum(proc) = ip
!  원복
         DO i=1,sort(proc)
             jwk(jperm(proc,i))=0
         ENDDO
         
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
	 cext_tmp = 0
!
     CALL Ext_nodes_R(np,next_m,nnode,nnode1,nnzi,celem,icoarse,iar,jar,cext_tmp,lcnode3)
     
	 ALLOCATE(inbdomR(np,np),nnbdomR(np),riR(np,np),siR(np,np),rintR(np,next_m),sintR(np,next_m))
     CALL Neighbor_node_ARP(np,nnode,next_m,next_m,celem,cext_tmp,lcnode3,inbdomR,nnbdomR,riR,siR,rintR,sintR)	 
	 
	 DEALLOCATE(cext_tmp)
! - - - - - - 
deallocate(index_cell,sort)
deallocate(rnbcnt,nbrecv)
deallocate(lcnode3)
deallocate(jwk,mark_own)

     RETURN

END
