SUBROUTINE subdomain_infor_mg
    USE md_geometry, ONLY: nelem,nelem_mg,num_neigh_mg,neigh_mg
    USE MD_MG_coord, ONLY: nmax1,nelem1,nelem2,nnei1,inei1,                  &
                           imapc,icoarse,icoarse1,icoarsef,ialv,inmax
    USE MD_MG_index, ONLY: nlevel,n_GC,nlevel_N,ioplv,n2_min, report_text,   &
                            mxnbne_mg,isend_m,irecv_m,                       &
                            isetup_comm,stg_iintf,stg_inodegl,stg_inbdc,     &
                            stg_ialvP,stg_inmax,stg_nnzc0,stg_nnzi,stg_nnzr, &
                            stg_fibuf,stg_frbuf,stg_ficnt,stg_frcnt,rt_ascii,&
                            stg_mg,stg_mg_init,stg_pushi,stg_pushr
    USE MD_parameter, ONLY: nf_max,ndim,ndom,nvpe
    USE MD_MG_matrix, ONLY: nnz1,nnzi1,iai1,jai1,iar1,jar1,                  &
                            iar2,jar2,Xintp1,Xrest1,ia1,ja1,Xrest2,          &
							iai,jai,iar,jar,iac,jac,Xrest,Xintp
                            
    USE Zmpi, ONLY: celem
    USE Zcoord1, ONLY : xloc_tmp
    
	USE MD_MPI_ARP, ONLY: inbdomA,nnbdomA,riA,siA,rintA,sintA,                                 &   ! NEW
	                      inbdomR,nnbdomR,riR,siR,rintR,sintR,                                 &
	                      inbdomP,nnbdomP,riP,siP,rintP,sintP
    
    IMPLICIT NONE

INTEGER(4):: i,j,k,ie,i1,i2,i3,i4,i5,i6,nd,nnd,j1,j2
INTEGER(4):: prc,np,ilv
INTEGER(4):: iu                  ! NEWUNIT 용 (PMG_infor)
INTEGER(4):: iu_prc(ndom)        ! NEWUNIT 용 — coarse 구간은 np 개 파일을 동시에 열어둠
INTEGER(4):: nintr,nintf,nneib,nelemt
INTEGER(4):: alstatus
CHARACTER(len=64) :: fout

INTEGER(4),DIMENSION(:),ALLOCATABLE::lnum,nnbdom,cinter,cintf,cext,sort,nnodegl_mg
INTEGER(4),DIMENSION(:,:),ALLOCATABLE::iperm,jperm,nbdom
INTEGER(4),DIMENSION(:,:),ALLOCATABLE::si,ri,sint,rint
INTEGER(4),DIMENSION(:,:),ALLOCATABLE::lcelem
! for PMG
INTEGER(4):: nnbdom1(ndom),cinter1(ndom),cintf1(ndom),cext1(ndom)
INTEGER(4):: nbdom1(ndom,ndom)
INTEGER(4),DIMENSION(:),ALLOCATABLE::celem1,celem0,imap
INTEGER(4),DIMENSION(:,:),ALLOCATABLE::iperm1,jperm1
INTEGER(4),DIMENSION(:,:),ALLOCATABLE::si1,ri1,sint1,rint1
REAL(8),DIMENSION(:,:),ALLOCATABLE::coord0,coord1
! for write out
INTEGER(4):: nnodep,nelemp
INTEGER(4):: nintr1,nintf1,nnodep1,nneib1
! new
INTEGER(4) ncolf1,ncolc1,ntmp,ntmpf,nnzt1,nnzr1,nnzr2,nelem0,nnzt
INTEGER(4) ialv_P(nlevel+1,ndom),iintf(nlevel,ndom),               &
           inodegl(nlevel,ndom),inbdc(nlevel,ndom)
INTEGER(4) nnzc0(ndom),nnzi(ndom),nnzr(ndom),nnodep0(ndom),nnodep1gl(ndom)
INTEGER(4) nc_min,ilv_test
INTEGER(4) kci,kcr,icnt_t          ! C011-3: finest 스테이징 pack 커서·카운트
REAL(8) tmp

!-------------------------------------------

! copy data
  nelem = nelem_mg
  np=ndom
!
  CALL read_mesh_FVM
  CALL Prep_fine_FVM
  CALL PREP_GMG

! 
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! FOR THE FINEST GRID
! 1: pre-P on global grid:
! 2: Write out the data
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

! NEW CHANGE. ioplv is only set on global coarsing

 ioplv = 0
 ilv_test = 0
  
500 CONTINUE
ilv_test = 0

  ilv = 1
CALL Predict_nelemt(np,ilv,nelem,nelemt) 

ALLOCATE(lnum(np),nbdom(np,np),nnbdom(np))
ALLOCATE(cext(np),cinter(np),cintf(np),sort(np),nnodegl_mg(np))
ALLOCATE(iperm(np,nelem),jperm(np,nelemt),stat=alstatus)
ALLOCATE(ri(np,np),si(np,np),rint(np,nelemt),sint(np,nelemt),stat=alstatus)
ALLOCATE(lcelem(np,nelemt),stat=alstatus)

IF (alstatus/=0) THEN
    WRITE(*,*)'not enough memory,serial-subdomain-fine'
    CALL STOP_MPI(report_text)
    STOP
ENDIF
  
! local array for iar
      ilv = 2
      nelem1 = ialv(ilv+1)-ialv(ilv)
      
      ncolf1 = ialv(ilv+1) - ialv(1)
      ncolc1 = ncolf1 -nelem

      ntmp = ncolc1 - nelem1
      
      nnzt1 = iar(ntmp+1)-1
      nnzr1 = iar(ncolc1+1)-iar(ntmp+1)
      
      ntmpf = ialv(ilv) - ialv(1) - nelem
      
      ALLOCATE(iar1(nelem1+1),jar1(nnzr1),Xrest1(nnzr1))

       iar1(1:nelem1+1)=iar(ntmp+1:ncolc1+1)- nnzt1

       jar1 (1:nnzr1)=jar(nnzt1+1:nnzt1+nnzr1) - ntmpf
       
       Xrest1(1:nnzr1) = Xrest(nnzt1+1:nnzt1+nnzr1)
       
      ALLOCATE(icoarse(nelem))   
      icoarse(1:nelem) = icoarsef(1:nelem)

CALL Domain_infor_FVM_fine(ndom,nf_max,nelem,nelemt,num_neigh_mg,neigh_mg,celem,lnum,lcelem,   &
         nbdom,nnbdom,cext,cinter,cintf,iperm,jperm,ri,si,rint,sint,nnodegl_mg,        &
         nelem1,nnzr1,iar1,jar1,icoarse) 

! write out data for the finest level
!%Output local domain
!/
allocate(isend_m(ndom),irecv_m(ndom))
isend_m = 1
irecv_m = 1
!/
! C011-3 통신 모드: 카운트 선패스 → 스테이징 버퍼 할당 (reader 가 SCATTERV 로 분배)
IF(isetup_comm.EQ.1) THEN
   IF(ALLOCATED(stg_ficnt)) DEALLOCATE(stg_ficnt,stg_frcnt)
   IF(ALLOCATED(stg_fibuf)) DEALLOCATE(stg_fibuf,stg_frbuf)
   ALLOCATE(stg_ficnt(np),stg_frcnt(np))
   DO prc=1,np
      icnt_t = 8
      DO i=1,cinter(prc)+cintf(prc)+cext(prc)
         icnt_t = icnt_t + 1 + num_neigh_mg(jperm(prc,i))
      ENDDO
      IF(nnbdom(prc).NE.0)  icnt_t = icnt_t + nnbdom(prc) + 2*(nnbdom(prc)+1)    &
                                   + (ri(prc,nnbdom(prc)+1)-1) + (si(prc,nnbdom(prc)+1)-1)
      IF(nnbdomA(prc).NE.0) icnt_t = icnt_t + nnbdomA(prc) + 2*(nnbdomA(prc)+1)  &
                                   + (riA(prc,nnbdomA(prc)+1)-1) + (siA(prc,nnbdomA(prc)+1)-1)
      IF(nnbdomR(prc).NE.0) icnt_t = icnt_t + nnbdomR(prc) + 2*(nnbdomR(prc)+1)  &
                                   + (riR(prc,nnbdomR(prc)+1)-1) + (siR(prc,nnbdomR(prc)+1)-1)
      stg_ficnt(prc) = icnt_t
      stg_frcnt(prc) = nnodegl_mg(prc)*ndim
   ENDDO
   ALLOCATE(stg_fibuf(SUM(stg_ficnt)),stg_frbuf(SUM(stg_frcnt)))
   kci = 0
   kcr = 0
ENDIF
!/
DO prc=1,np
   nnodep=cinter(prc)+cintf(prc)+cext(prc)   !total number of nodes
   nelemp=lnum(prc)                          !number of elements
   nintr=cinter(prc)                         !number of internal nodes
   nintf=cinter(prc)+cintf(prc)              !number of interface nodes
   nneib=nnbdom(prc)                         !number of neighboring domains
   nnd = nnodegl_mg(prc)

  IF(isetup_comm.EQ.0) THEN
   ! I0.3 포맷 — 읽기측 2_read_mesh_MPI.f90 과 동일 유지 (np>999 자동 확장)
   WRITE(fout,'(A,I0.3,A)') 'MG_tmp/part', prc, '.out'
   OPEN(newunit=iu_prc(prc),file=fout,status='replace')

   WRITE(iu_prc(prc),*) nelemp,nintr,nintf,nnodep,nneib,nnd

   WRITE(iu_prc(prc),*) nnbdomA(prc),   nnbdomR(prc)         ! NEW

   DO i=1,nnodep
      ie=jperm(prc,i)
      j = num_neigh_mg(ie)
      WRITE(iu_prc(prc),*) j, (iperm(prc,neigh_mg(k,ie)),k=1,j)
   ENDDO

   DO i=1,nnd
      WRITE(iu_prc(prc),*) (xloc_tmp(jperm(prc,i),j),j=1,ndim)
   ENDDO
  ELSE
!  통신 모드 pack — 파일 레코드와 동일 순서·동일 값 (coord 는 rt_ascii 라운딩)
   stg_fibuf(kci+1) = nelemp;  stg_fibuf(kci+2) = nintr
   stg_fibuf(kci+3) = nintf;   stg_fibuf(kci+4) = nnodep
   stg_fibuf(kci+5) = nneib;   stg_fibuf(kci+6) = nnd
   stg_fibuf(kci+7) = nnbdomA(prc); stg_fibuf(kci+8) = nnbdomR(prc)
   kci = kci + 8
   DO i=1,nnodep
      ie=jperm(prc,i)
      j = num_neigh_mg(ie)
      stg_fibuf(kci+1) = j
      DO k=1,j
         stg_fibuf(kci+1+k) = iperm(prc,neigh_mg(k,ie))
      ENDDO
      kci = kci + 1 + j
   ENDDO
   DO i=1,nnd
      DO j=1,ndim
         stg_frbuf(kcr+j) = rt_ascii(xloc_tmp(jperm(prc,i),j))
      ENDDO
      kcr = kcr + ndim
   ENDDO
  ENDIF
!
  IF(nnbdom(prc).NE.0) THEN
   IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(nbdom(prc,i),i=1,nnbdom(prc))
   WRITE(iu_prc(prc),*)(ri(prc,i),i=1,nnbdom(prc)+1)
   WRITE(iu_prc(prc),*)(si(prc,i),i=1,nnbdom(prc)+1)
   WRITE(iu_prc(prc),*)(iperm(prc,rint(prc,i)),i=1,ri(prc,nnbdom(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm(prc,sint(prc,i)),i=1,si(prc,nnbdom(prc)+1)-1)
   ELSE
   DO i=1,nnbdom(prc)
      stg_fibuf(kci+i) = nbdom(prc,i)
   ENDDO
   kci = kci + nnbdom(prc)
   DO i=1,nnbdom(prc)+1
      stg_fibuf(kci+i) = ri(prc,i)
   ENDDO
   kci = kci + nnbdom(prc)+1
   DO i=1,nnbdom(prc)+1
      stg_fibuf(kci+i) = si(prc,i)
   ENDDO
   kci = kci + nnbdom(prc)+1
   DO i=1,ri(prc,nnbdom(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,rint(prc,i))
   ENDDO
   kci = kci + ri(prc,nnbdom(prc)+1)-1
   DO i=1,si(prc,nnbdom(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,sint(prc,i))
   ENDDO
   kci = kci + si(prc,nnbdom(prc)+1)-1
   ENDIF
!/
   isend_m(prc) = max(isend_m(prc),si(prc,nnbdom(prc)+1)-1)
   irecv_m(prc) = max(irecv_m(prc),ri(prc,nnbdom(prc)+1)-1)
!/
  ENDIF

!
! NEW for SR-for A
  IF(nnbdomA(prc).NE.0) THEN
   IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(inbdomA(prc,i),i=1,nnbdomA(prc))
   WRITE(iu_prc(prc),*)(riA(prc,i),i=1,nnbdomA(prc)+1)
   WRITE(iu_prc(prc),*)(siA(prc,i),i=1,nnbdomA(prc)+1)
   WRITE(iu_prc(prc),*)(iperm(prc,rintA(prc,i)),i=1,riA(prc,nnbdomA(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm(prc,sintA(prc,i)),i=1,siA(prc,nnbdomA(prc)+1)-1)
   ELSE
   DO i=1,nnbdomA(prc)
      stg_fibuf(kci+i) = inbdomA(prc,i)
   ENDDO
   kci = kci + nnbdomA(prc)
   DO i=1,nnbdomA(prc)+1
      stg_fibuf(kci+i) = riA(prc,i)
   ENDDO
   kci = kci + nnbdomA(prc)+1
   DO i=1,nnbdomA(prc)+1
      stg_fibuf(kci+i) = siA(prc,i)
   ENDDO
   kci = kci + nnbdomA(prc)+1
   DO i=1,riA(prc,nnbdomA(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,rintA(prc,i))
   ENDDO
   kci = kci + riA(prc,nnbdomA(prc)+1)-1
   DO i=1,siA(prc,nnbdomA(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,sintA(prc,i))
   ENDDO
   kci = kci + siA(prc,nnbdomA(prc)+1)-1
   ENDIF
  ENDIF
!
! NEW for SR-for R
  IF(nnbdomR(prc).NE.0) THEN
   IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(inbdomR(prc,i),i=1,nnbdomR(prc))
   WRITE(iu_prc(prc),*)(riR(prc,i),i=1,nnbdomR(prc)+1)
   WRITE(iu_prc(prc),*)(siR(prc,i),i=1,nnbdomR(prc)+1)
   WRITE(iu_prc(prc),*)(iperm(prc,rintR(prc,i)),i=1,riR(prc,nnbdomR(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm(prc,sintR(prc,i)),i=1,siR(prc,nnbdomR(prc)+1)-1)
   ELSE
   DO i=1,nnbdomR(prc)
      stg_fibuf(kci+i) = inbdomR(prc,i)
   ENDDO
   kci = kci + nnbdomR(prc)
   DO i=1,nnbdomR(prc)+1
      stg_fibuf(kci+i) = riR(prc,i)
   ENDDO
   kci = kci + nnbdomR(prc)+1
   DO i=1,nnbdomR(prc)+1
      stg_fibuf(kci+i) = siR(prc,i)
   ENDDO
   kci = kci + nnbdomR(prc)+1
   DO i=1,riR(prc,nnbdomR(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,rintR(prc,i))
   ENDDO
   kci = kci + riR(prc,nnbdomR(prc)+1)-1
   DO i=1,siR(prc,nnbdomR(prc)+1)-1
      stg_fibuf(kci+i) = iperm(prc,sintR(prc,i))
   ENDDO
   kci = kci + siR(prc,nnbdomR(prc)+1)-1
   ENDIF
  ENDIF
!

  IF(isetup_comm.EQ.0) CLOSE(iu_prc(prc))

ENDDO

! NEW for ArP - - - - - - - - - - 
! DEALLOCATE
    DEALLOCATE(inbdomA,nnbdomA)
    DEALLOCATE(riA,siA,rintA,sintA)
    
    DEALLOCATE(inbdomR,nnbdomR)
    DEALLOCATE(riR,siR,rintR,sintR)	
! - - - - - - - - - - - - - - - - - - 

   iintf(1,1:np) = cinter(1:np)+cintf(1:np)
   inodegl(1,1:np) = nnodegl_mg(1:np)
   inbdc(1,1:np) =  nnbdom(1:np)
   ialv_P(1,1:np) = 1
   ialv_P(2,1:np) = ialv_P(1,1:np) + lnum(1:np) 
   nnzc0 = 0
   nnzr = 0
   nnzi = 0
   
   nnodep0(1:np) = lnum(1:np)
! 
DEALLOCATE(lnum,nbdom,nnbdom)
DEALLOCATE(cext,cinter,cintf,sort,nnodegl_mg)
DEALLOCATE(ri,si,rint,sint)
DEALLOCATE(lcelem)

DEALLOCATE(icoarse)

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! the coarse level
ALLOCATE(celem0(nelem),coord0(ndim,nelem))
celem0 = celem

DO i = 1,nelem
coord0(1:ndim,i) = xloc_tmp(i,1:ndim)
ENDDO

nelem0 = nelem

IF(isetup_comm.EQ.0) THEN
DO prc=1,np
   WRITE(fout,'(A,I0.3,A)') 'MG_tmp/part_MG', prc, '.out'
   OPEN(newunit=iu_prc(prc),file=fout,status='replace')
ENDDO
ELSE
CALL stg_mg_init(np)          ! C011-4: prc별 성장형 스트림 (재진입 시 리셋)
ENDIF
!/ 
  mxnbne_mg = 0
!/
DO ILV = 2,nlevel_N

      nelem1 = ialv(ilv+1)-ialv(ilv)
      
      ncolf1 = ialv(ilv+1) - ialv(1)
      ncolc1 = ncolf1 -nelem
! AC
      ntmp = ncolc1 - nelem1
      nnzt = iac(ntmp+1)-1
      
      nnz1 = iac(ncolc1+1)-iac(ntmp+1)
      
      ALLOCATE(ia1(nelem1+1),ja1(nnz1))
      
      ia1(1:nelem1+1)=iac(ntmp+1:ncolc1+1) - nnzt
      ja1 (1:nnz1) = jac(nnzt+1:nnzt+nnz1) - ntmp
! P
      ncolf1 = ialv(ilv) - ialv(1)
      ncolc1 = ncolf1 -nelem
      ntmp = ncolf1 - nelem0
      nnzt = iai(ntmp+1)-1
      nnzi1 = iai(ncolf1+1)-iai(ntmp+1)
      ALLOCATE(iai1(nelem0+1),jai1(nnzi1),Xintp1(nnzi1))

      iai1(1:nelem0+1) = iai(ntmp+1:ncolf1+1) - nnzt
      jai1 (1:nnzi1) = jai(nnzt+1:nnzt+nnzi1) - ncolc1
      Xintp1(1:nnzi1) = Xintp(nnzt+1:nnzt+nnzi1)
! R
     IF(ilv.EQ.nlevel_N) GOTO 10
     
      nelem2 = ialv(ilv+2)-ialv(ilv+1)
      ncolf1 = ialv(ilv+2) - ialv(1)
      ncolc1 = ncolf1 -nelem
      ntmp = ncolc1 - nelem2
      nnzt = iar(ntmp+1)-1
      nnzr2 = iar(ncolc1+1)-iar(ntmp+1)
      
      ntmpf = ialv(ilv+1) - ialv(1) - nelem1
      
      ALLOCATE(iar2(nelem2+1),jar2(nnzr2),Xrest2(nnzr2))

       iar2(1:nelem2+1)=iar(ntmp+1:ncolc1+1)- nnzt

       jar2 (1:nnzr2)=jar(nnzt+1:nnzt+nnzr2) - ntmpf
       Xrest2(1:nnzr2) = Xrest(nnzt+1:nnzt+nnzr2) 
       
10     CONTINUE
       
! local array
nmax1 = inmax(ilv)
ALLOCATE(nnei1(nelem1),inei1(nmax1,nelem1))   
DO i = 1,nelem1
    i1 = ia1(i)
    i2 = ia1(i+1)-1
    nd = i2-i1
    nnei1(i) = nd
    
    nnd = 0
     DO j=i1,i2
     k = ja1(j)
     IF(k.EQ.i) CYCLE
     nnd=nnd+1
     inei1(nnd,i) = k
     ENDDO
    
ENDDO
!
ALLOCATE(celem1(nelem1),imap(nelem1),coord1(ndim,nelem1),icoarse1(nelem1))

ncolf1 = ialv(ilv+1) - ialv(1)
ncolc1 = ncolf1 -nelem
ntmp = ncolc1 - nelem1
imap(1:nelem1) = imapc(ntmp+1:ntmp+nelem1)
ntmp = ncolf1 - nelem1
IF(ilv.NE.nlevel_N) THEN
icoarse1(1:nelem1) = icoarsef(ntmp+1:ntmp+nelem1)
ENDIF

DO i=1,nelem1
    j = imap(i)
    celem1(i) = celem0(j)
    coord1(1:ndim,i) = coord0(1:ndim,j)
ENDDO

DEALLOCATE(imap,coord0)

CALL Predict_nelemt(np,ilv,nelem1,nelemt)

ALLOCATE(iperm1(np,nelem1),jperm1(np,nelemt),stat=alstatus)
ALLOCATE(ri1(np,np),si1(np,np))
ALLOCATE(rint1(np,nelemt),sint1(np,nelemt),stat=alstatus)

     IF (alstatus/=0) THEN
         WRITE(*,*)'not enough memory,serial-subdomain-coarse -1'
         CALL STOP_MPI(report_text)
         STOP
     ENDIF
!
    
IF(ilv.NE.nlevel_N) THEN
CALL Domain_infor_FVM_coarse(ilv,ndom,nmax1,nelem1,nelemt,nelem2,nnzr2,iar2,jar2,icoarse1,nnei1,inei1,celem1,      &
         nbdom1,nnbdom1,cext1,cinter1,cintf1,iperm1,jperm1,ri1,si1,rint1,sint1,nnodep1gl,            &
         nelem0,celem0,nnzi1,iai1,jai1)
ELSE
CALL Domain_infor_FVM_coarsest(ilv,ndom,nmax1,nelem1,nelemt,nnei1,inei1,celem1,nbdom1,nnbdom1,                        &
         cext1,cinter1,cintf1,iperm1,jperm1,ri1,si1,rint1,sint1,nnodep1gl,                           &
         nelem0,celem0,nnzi1,iai1,jai1)
ENDIF

DEALLOCATE(celem0,icoarse1,nnei1,inei1)
  
! NEW: finding minimun cells - - - - - - - - - - - - !
  nc_min = 1000
  ilv_test = 0
  
  IF(ioplv.EQ.1) THEN
  DO prc=1,np
     nintf1=cinter1(prc)+cintf1(prc)
     nc_min = MIN(nintf1,nc_min)
  ENDDO
!
    i1 = MAX(INT(n2_min/2),1)
  IF(nc_min.LE.i1)   ilv_test = 1
  IF(nc_min.EQ.0)   ilv_test = 2   
  ENDIF
  
  IF(ilv_test.EQ.2) GOTO 100
! IF: ilv_test = 1-> select this level is the coarsest
! IF: ilv_test = 2-> select previous one (run again)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - !

! wrting out file 
DO prc=1,np

   nnodep1=cinter1(prc)+cintf1(prc)+cext1(prc)
   nintr1=cinter1(prc)
   nintf1=cinter1(prc)+cintf1(prc)
   nneib1=nnbdom1(prc)
   nnd = nnodep1gl(prc)
!
   iintf(ilv,prc) = nintf1
   inodegl(ilv,prc) = nnd
   inbdc(ilv,prc) =  nneib1
   ialv_P(ilv+1,prc) = ialv_P(ilv,prc) + nnodep1
!
  IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)'coarse'
   WRITE(iu_prc(prc),*) nintf1,nnodep1,nneib1,nnd

! NEW
   IF(ilv.NE.nlevel_N) THEN
      WRITE(iu_prc(prc),*)  nnbdomA(prc),nnbdomR(prc),nnbdomP(prc)
   ELSE
      WRITE(iu_prc(prc),*)  nnbdomA(prc),nnbdomP(prc)
   ENDIF
!
   DO i=1,nnd
       j=jperm1(prc,i)     ! to global of coarse 1
       WRITE(iu_prc(prc),*) coord1(1:ndim,j)
   ENDDO
  ELSE
!  통신 모드 pack — 파일 레코드와 동일 순서 (마커 'coarse' 등은 무데이터라 생략)
   CALL stg_pushi(prc,nintf1)
   CALL stg_pushi(prc,nnodep1)
   CALL stg_pushi(prc,nneib1)
   CALL stg_pushi(prc,nnd)
   CALL stg_pushi(prc,nnbdomA(prc))
   IF(ilv.NE.nlevel_N) CALL stg_pushi(prc,nnbdomR(prc))
   CALL stg_pushi(prc,nnbdomP(prc))
   DO i=1,nnd
       j=jperm1(prc,i)     ! to global of coarse 1
       DO j1=1,ndim
          CALL stg_pushr(prc,rt_ascii(coord1(j1,j)))
       ENDDO
   ENDDO
  ENDIF
   
! test
   IF(si1(prc,nnbdom1(prc)+1)-1.GE.nelemt) THEN
   WRITE(*,*)'PMG error: nelemt is small'
   STOP
   ENDIF
!
IF(nnbdom1(prc).NE.0) THEN
  IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(nbdom1(prc,i),i=1,nnbdom1(prc))
   WRITE(iu_prc(prc),*)(ri1(prc,i),i=1,nnbdom1(prc)+1)
   WRITE(iu_prc(prc),*)(si1(prc,i),i=1,nnbdom1(prc)+1)
   WRITE(iu_prc(prc),*)(iperm1(prc,rint1(prc,i)),i=1,ri1(prc,nnbdom1(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm1(prc,sint1(prc,i)),i=1,si1(prc,nnbdom1(prc)+1)-1)
  ELSE
   DO i=1,nnbdom1(prc)
      CALL stg_pushi(prc,nbdom1(prc,i))
   ENDDO
   DO i=1,nnbdom1(prc)+1
      CALL stg_pushi(prc,ri1(prc,i))
   ENDDO
   DO i=1,nnbdom1(prc)+1
      CALL stg_pushi(prc,si1(prc,i))
   ENDDO
   DO i=1,ri1(prc,nnbdom1(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,rint1(prc,i)))
   ENDDO
   DO i=1,si1(prc,nnbdom1(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,sint1(prc,i)))
   ENDDO
  ENDIF
!/
   isend_m(prc) = max(isend_m(prc),si1(prc,nnbdom1(prc)+1)-1)
   irecv_m(prc) = max(irecv_m(prc),ri1(prc,nnbdom1(prc)+1)-1)
!/
ENDIF
!
! NEW for SR for A
IF(nnbdomA(prc).NE.0) THEN
  IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(inbdomA(prc,i),i=1,nnbdomA(prc))
   WRITE(iu_prc(prc),*)(riA(prc,i),i=1,nnbdomA(prc)+1)
   WRITE(iu_prc(prc),*)(siA(prc,i),i=1,nnbdomA(prc)+1)
   WRITE(iu_prc(prc),*)(iperm1(prc,rintA(prc,i)),i=1,riA(prc,nnbdomA(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm1(prc,sintA(prc,i)),i=1,siA(prc,nnbdomA(prc)+1)-1)
  ELSE
   DO i=1,nnbdomA(prc)
      CALL stg_pushi(prc,inbdomA(prc,i))
   ENDDO
   DO i=1,nnbdomA(prc)+1
      CALL stg_pushi(prc,riA(prc,i))
   ENDDO
   DO i=1,nnbdomA(prc)+1
      CALL stg_pushi(prc,siA(prc,i))
   ENDDO
   DO i=1,riA(prc,nnbdomA(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,rintA(prc,i)))
   ENDDO
   DO i=1,siA(prc,nnbdomA(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,sintA(prc,i)))
   ENDDO
  ENDIF
ENDIF
	
! NEW for SR for R
IF(ilv.NE.nlevel_N) THEN
IF(nnbdomR(prc).NE.0) THEN
  IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(inbdomR(prc,i),i=1,nnbdomR(prc))
   WRITE(iu_prc(prc),*)(riR(prc,i),i=1,nnbdomR(prc)+1)
   WRITE(iu_prc(prc),*)(siR(prc,i),i=1,nnbdomR(prc)+1)
   WRITE(iu_prc(prc),*)(iperm1(prc,rintR(prc,i)),i=1,riR(prc,nnbdomR(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm1(prc,sintR(prc,i)),i=1,siR(prc,nnbdomR(prc)+1)-1)
  ELSE
   DO i=1,nnbdomR(prc)
      CALL stg_pushi(prc,inbdomR(prc,i))
   ENDDO
   DO i=1,nnbdomR(prc)+1
      CALL stg_pushi(prc,riR(prc,i))
   ENDDO
   DO i=1,nnbdomR(prc)+1
      CALL stg_pushi(prc,siR(prc,i))
   ENDDO
   DO i=1,riR(prc,nnbdomR(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,rintR(prc,i)))
   ENDDO
   DO i=1,siR(prc,nnbdomR(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,sintR(prc,i)))
   ENDDO
  ENDIF
ENDIF
!
ENDIF

! NEW for SR for P
IF(nnbdomP(prc).NE.0) THEN
  IF(isetup_comm.EQ.0) THEN
   WRITE(iu_prc(prc),*)(inbdomP(prc,i),i=1,nnbdomP(prc))
   WRITE(iu_prc(prc),*)(riP(prc,i),i=1,nnbdomP(prc)+1)
   WRITE(iu_prc(prc),*)(siP(prc,i),i=1,nnbdomP(prc)+1)
   WRITE(iu_prc(prc),*)(iperm1(prc,rintP(prc,i)),i=1,riP(prc,nnbdomP(prc)+1)-1)
   WRITE(iu_prc(prc),*)(iperm1(prc,sintP(prc,i)),i=1,siP(prc,nnbdomP(prc)+1)-1)
  ELSE
   DO i=1,nnbdomP(prc)
      CALL stg_pushi(prc,inbdomP(prc,i))
   ENDDO
   DO i=1,nnbdomP(prc)+1
      CALL stg_pushi(prc,riP(prc,i))
   ENDDO
   DO i=1,nnbdomP(prc)+1
      CALL stg_pushi(prc,siP(prc,i))
   ENDDO
   DO i=1,riP(prc,nnbdomP(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,rintP(prc,i)))
   ENDDO
   DO i=1,siP(prc,nnbdomP(prc)+1)-1
      CALL stg_pushi(prc,iperm1(prc,sintP(prc,i)))
   ENDDO
  ENDIF
ENDIF

   IF(isetup_comm.EQ.0) WRITE(iu_prc(prc),*)'P-1'

   nnodep = nnodep0(prc)
   DO i=1,nnodep
       j=jperm(prc,i)   ! to global
       i1 = iai1(j)
       i2 = iai1(j+1)-1
       nnd = i2-i1+1
      IF(isetup_comm.EQ.0) THEN
       WRITE(iu_prc(prc),*) nnd,(iperm1(prc,jai1(k)),k=i1,i2)

       WRITE(iu_prc(prc),*) (Xintp1(k),k=i1,i2)
      ELSE
       CALL stg_pushi(prc,nnd)
       DO k=i1,i2
          CALL stg_pushi(prc,iperm1(prc,jai1(k)))
       ENDDO
       DO k=i1,i2
          CALL stg_pushr(prc,rt_ascii(Xintp1(k)))
       ENDDO
      ENDIF

       nnzi(prc) = nnzi(prc) + nnd

   ENDDO

   nnodep0(prc) = nnodep1

! R
   IF(isetup_comm.EQ.0) WRITE(iu_prc(prc),*)'R-1'
   DO i=1,nnodep1
       j=jperm1(prc,i)   ! to global
       i1 = iar1(j)
       i2 = iar1(j+1)-1
       nnd = i2-i1+1
      IF(isetup_comm.EQ.0) THEN
       WRITE(iu_prc(prc),*) nnd,(iperm(prc,jar1(k)),k=i1,i2)

       WRITE(iu_prc(prc),*) (Xrest1(k),k=i1,i2)
      ELSE
       CALL stg_pushi(prc,nnd)
       DO k=i1,i2
          CALL stg_pushi(prc,iperm(prc,jar1(k)))
       ENDDO
       DO k=i1,i2
          CALL stg_pushr(prc,rt_ascii(Xrest1(k)))
       ENDDO
      ENDIF

       nnzr(prc) = nnzr(prc) + nnd

   ENDDO

! Ac
   IF(isetup_comm.EQ.0) WRITE(iu_prc(prc),*)'Ac-1'
   DO i=1,nnodep1
       j=jperm1(prc,i)   ! to global
       i1 = ia1(j)
       i2 = ia1(j+1)-1
       nnd = i2-i1+1
      IF(isetup_comm.EQ.0) THEN
       WRITE(iu_prc(prc),*) nnd,(iperm1(prc,ja1(k)),k=i1,i2)
      ELSE
       CALL stg_pushi(prc,nnd)
       DO k=i1,i2
          CALL stg_pushi(prc,iperm1(prc,ja1(k)))
       ENDDO
      ENDIF

       mxnbne_mg = max(nnd,mxnbne_mg)

       nnzc0(prc) = nnzc0(prc) + nnd

   ENDDO

! !    GLOBAL

   IF (n_GC.EQ.0) CYCLE

   IF ((ilv_test.EQ.1).OR.(ilv.EQ.nlevel_N)) THEN

    IF(isetup_comm.EQ.0) THEN
     write(iu_prc(prc),*)'A_GC'
     write(iu_prc(prc),*) nnodep1,nelem1,nnz1
     DO i=1,nnodep1
       j=jperm1(prc,i)   ! to global
       write(iu_prc(prc),*) j
     ENDDO

     DO i=1,nelem1
       i1 = ia1(i)
       i2 = ia1(i+1)-1
       nnd = i2-i1+1
       write(iu_prc(prc),*) nnd,(ja1(k),k=i1,i2),coord1(1:ndim,i)

       mxnbne_mg = max(nnd,mxnbne_mg)

     ENDDO
    ELSE
     CALL stg_pushi(prc,nnodep1)
     CALL stg_pushi(prc,nelem1)
     CALL stg_pushi(prc,nnz1)
     DO i=1,nnodep1
       CALL stg_pushi(prc,jperm1(prc,i))
     ENDDO
     DO i=1,nelem1
       i1 = ia1(i)
       i2 = ia1(i+1)-1
       nnd = i2-i1+1
       CALL stg_pushi(prc,nnd)
       DO k=i1,i2
          CALL stg_pushi(prc,ja1(k))
       ENDDO
       DO k=1,ndim
          CALL stg_pushr(prc,rt_ascii(coord1(k,i)))
       ENDDO

       mxnbne_mg = max(nnd,mxnbne_mg)

     ENDDO
    ENDIF

   ENDIF

ENDDO

! NEW for A, R,P
! allocate:
    DEALLOCATE(inbdomA,nnbdomA)
    DEALLOCATE(riA,siA,rintA,sintA)
!
    IF(ilv.NE.nlevel_N) THEN
    DEALLOCATE(inbdomR,nnbdomR)
    DEALLOCATE(riR,siR,rintR,sintR)	
    ENDIF
!
    DEALLOCATE(inbdomP,nnbdomP)
    DEALLOCATE(riP,siP,rintP,sintP)	
! - - - - - 
IF(ilv.EQ.nlevel_N) GOTO 100
IF(ilv_test.EQ.1) GOTO 100

! updating for next level
  nelem0 = nelem1
  ALLOCATE(celem0(nelem0),coord0(ndim,nelem0))
  celem0 = celem1
  coord0 = coord1
  DEALLOCATE(celem1,coord1)
  
  Xrest1(1:nnzr2) = Xrest2(1:nnzr2)
  iar1(1:nelem2+1) = iar2(1:nelem2+1)
  jar1(1:nnzr2) = jar2(1:nnzr2)
  
  DEALLOCATE(iperm,jperm)
  ALLOCATE(iperm(np,nelem1),jperm(np,nelemt))
  
  iperm = iperm1
  jperm = jperm1
  
  DEALLOCATE(iperm1,jperm1,Xrest2,iar2,jar2)
  
  DEALLOCATE(si1,ri1,sint1,rint1)
  DEALLOCATE(ia1,ja1)

100 CONTINUE
     
  DEALLOCATE(iai1,jai1,Xintp1)
  
IF(ilv_test.EQ.1.OR.ilv_test.EQ.2) GOTO 200   

ENDDO

200 CONTINUE

IF(ilv_test.EQ.1.) nlevel_N = ilv
IF(ilv_test.EQ.2) nlevel_N = ilv-1
    
DEALLOCATE(ia1,ja1)
    
IF(isetup_comm.EQ.0) THEN
DO prc=1,np
   CLOSE(iu_prc(prc))
ENDDO
ENDIF
!
   IF(isetup_comm.EQ.0) THEN
   OPEN(newunit=iu,file='MG_tmp/PMG_infor',status='replace')
   DO prc = 1,ndom
      DO ilv=1,nlevel_N
      WRITE(iu,*) iintf(ilv,prc),inodegl(ilv,prc),inbdc(ilv,prc),inmax(ilv)
      ENDDO

      DO ilv = 1,nlevel_N+1
      WRITE(iu,*) ialv_P(ilv,prc)
      ENDDO

      WRITE(iu,*) nnzc0(prc),nnzi(prc),nnzr(prc)

   ENDDO
   CLOSE(iu)
   ELSE
!  통신 모드: PMG_infor 내용을 스테이징에 적재 (reader 가 BCAST 로 분배).
!  ioplv=1 재진입(GOTO 500)으로 2회 실행될 수 있어 재할당 가드 필요
   IF(ALLOCATED(stg_iintf)) DEALLOCATE(stg_iintf,stg_inodegl,stg_inbdc,     &
                                       stg_ialvP,stg_inmax,stg_nnzc0,       &
                                       stg_nnzi,stg_nnzr)
   ALLOCATE(stg_iintf(nlevel_N,ndom),stg_inodegl(nlevel_N,ndom),            &
            stg_inbdc(nlevel_N,ndom),stg_ialvP(nlevel_N+1,ndom),            &
            stg_inmax(nlevel_N),stg_nnzc0(ndom),stg_nnzi(ndom),stg_nnzr(ndom))
   stg_iintf(1:nlevel_N,:)    = iintf(1:nlevel_N,:)
   stg_inodegl(1:nlevel_N,:)  = inodegl(1:nlevel_N,:)
   stg_inbdc(1:nlevel_N,:)    = inbdc(1:nlevel_N,:)
   stg_ialvP(1:nlevel_N+1,:)  = ialv_P(1:nlevel_N+1,:)
   stg_inmax(1:nlevel_N)      = inmax(1:nlevel_N)
   stg_nnzc0 = nnzc0
   stg_nnzi  = nnzi
   stg_nnzr  = nnzr
   ENDIF

! NEW
   IF(ilv_test.EQ.2) THEN
    IF(ALLOCATED(iar2))      DEALLOCATE(iar2)
    IF(ALLOCATED(jar2)) DEALLOCATE(jar2)
    IF(ALLOCATED(Xrest2)) DEALLOCATE(Xrest2)
    IF(ALLOCATED(celem1)) DEALLOCATE(celem1)
    IF(ALLOCATED(coord1)) DEALLOCATE(coord1) 
    IF(ALLOCATED(iperm)) DEALLOCATE(iperm)    
    IF(ALLOCATED(jperm)) DEALLOCATE(jperm)
    IF(ALLOCATED(iperm1)) DEALLOCATE(iperm1)    
    IF(ALLOCATED(jperm1)) DEALLOCATE(jperm1) 
    IF(ALLOCATED(si1)) DEALLOCATE(si1)  
    IF(ALLOCATED(ri1)) DEALLOCATE(ri1)    
    IF(ALLOCATED(sint1)) DEALLOCATE(sint1) 
    IF(ALLOCATED(rint1)) DEALLOCATE(rint1) 
    IF(ALLOCATED(iar1))      DEALLOCATE(iar1)
    IF(ALLOCATED(jar1)) DEALLOCATE(jar1)
    IF(ALLOCATED(Xrest1)) DEALLOCATE(Xrest1)

   GOTO 500
   ENDIF
   
! 
   cinter1(1:np) = cinter1(1:np) + cintf1(1:np)
   i1 = minval(cinter1(1:np)) 
   i2 = maxval(cinter1(1:np)) 
  
WRITE(999,*)'coarsest level info.'
WRITE(999,*)'minimun cells',i1
WRITE(999,*)'maximun cells',i2
IF(i1.EQ.0) THEN
    WRITE(999,*)'coarsest level zero cells'
ENDIF

! for optimal nlv_go
i1 = SUM(cext1(1:np))
tmp = REAL(i1)/REAL(nnz1)
WRITE(999,*)'optimal nlv_gobal, nnz,ext,beta',nnz1, i1, tmp
IF(tmp.GT.0.3)  WRITE(999,*)'increse "nlv_glo" in mg.in file'
IF(tmp.LT.0.1)  WRITE(999,*)'reduce "nlv_glo" in mg.in file'
!/
!
DEALLOCATE(imapc,icoarsef,ialv,iperm,jperm,num_neigh_mg,neigh_mg)
DEALLOCATE(iar1,jar1,Xrest1,iai,jai,iar,jar,iac,jac)
DEALLOCATE(Xrest,Xintp)
DEALLOCATE(coord1,inmax)

RETURN
END

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

SUBROUTINE Predict_nelemt(np,ilv,nelem,nelemt)

IMPLICIT NONE
!
INTEGER(4)::np,ilv,nelem,nelemt

IF(ilv.LE.2) THEN
    
IF(np.LE.3) THEN
    nelemt = nelem
ELSEIF(np.LE.10) THEN
    nelemt = INT(nelem/np*4)
ELSEIF(np.LE.50) THEN
    nelemt = INT(nelem/np*8)
ELSE
    nelemt = INT(nelem/np*20)    
ENDIF

ELSEIF(ilv.EQ.3) THEN
    
IF(np.LE.5) THEN
    nelemt = nelem
ELSEIF(np.LE.10) THEN
    nelemt = INT(nelem/np*8)
ELSEIF(np.LE.50) THEN
    nelemt = INT(nelem/np*10)
ELSE
    nelemt = INT(nelem/np*50)    
ENDIF

ELSE

    nelemt = nelem
    
ENDIF

IF(nelemt.LE.2000) nelemt = 2000

! 
RETURN
END
    
