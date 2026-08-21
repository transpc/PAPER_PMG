SUBROUTINE read_mesh_MPI
	  
USE MD_parameter, ONLY: ndim, ndom, nf_max
USE MD_geometry, ONLY: nelem,nnode,nnodegl,nelemgl,coord,num_neigh,             &
                       e_neigh
USE MD_MPI, ONLY: nnbd,nintf,nintr,myrank,spt,rpt,sintf,rintf,nbdom
USE MD_MPI_MG, ONLY: nnbd1,nintf1,nnzr1,spt1,rpt1,sintf1,rintf1,nbdom1,         &
                     iintf,inodegl, nintfs, inbdc,                              &
                     ibdomc,isptc,irptc,isintfc,irintfc,isiac,iriac,            &
                     nnsend_m,nnrecv_m   
USE MD_MG_coord, ONLY: nnode1,nelem1,nnode1gl,coord1,ialv,                     &
                       nnods,ncolc,ncolf,coordc,inmax
USE MD_MG_matrix, ONLY: nnzi1,iai1,jai1,iar1,jar1,Xintp1,Xrest1,               &
                        nnz1,ia1,ja1,ju1,au1,diagrc,r,rt,rc,rs,e,et,es,        &
                        auc,aus,Xrest,Xintp,iac,jac,juc,ias,jas,jus,           &
                        iai,jai,iar,jar,nnzc0,nnzi,nnzr,nnzs
USE MD_MG_index, ONLY: mxnbne,nlevel,n_GC,nlevel_N,mxnbne_mg,isend_m,irecv_m,   &
                       isetup_comm,stg_iintf,stg_inodegl,stg_inbdc,             &
                       stg_ialvP,stg_inmax,stg_nnzc0,stg_nnzi,stg_nnzr,         &
                       stg_fibuf,stg_frbuf,stg_ficnt,stg_frcnt,stg_mg
USE MD_MG_Global_C, ONLY: nlv_glo  ,nnodeG,nnzG,nnodeC,imapG,imapGZ,           &
                          iaG,jaG,juG,eG,eG0,rG,rG0,auG,auG0,                  &
                          coordG
! NEW
USE MD_MPI_ARP, ONLY: nnbdA, nbdomA, sptA, rptA, rintfA, sintfA,               &
                      nnbdR, nbdomR, sptR, rptR, rintfR, sintfR,               &
					  inbdcA,inbdcR,inbdcP,                                    &
					  ibdomcA,isptcA,irptcA,isintfcA,irintfcA,                 &
					  ibdomcR,isptcR,irptcR,isintfcR,irintfcR,                 &
					  ibdomcP,isptcP,irptcP,isintfcP,irintfcP,                 &
                      nnsend_mA,nnrecv_mA,nnsend_mR,nnrecv_mR,nnsend_mP,nnrecv_mP	
	  
IMPLICIT NONE

!DEC$IF defined (mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ENDIF

character(len=64)::fout
integer:: ie,ip,i,j,k,nnd,iu
integer::alstatus,i1,i2,i3,i4,ierr
real*8 tmp(10*mxnbne)
INTEGER(4)::ilv,ntmpc,ntmp,ntmpf,ncolf1,ncolc1,ncolc2,nnode0,nnzt,nnzt1
INTEGER(4), DIMENSION(:), ALLOCATABLE::id
! C011-3: finest 통신 모드 수신 버퍼·커서
INTEGER(4)::kci,kcr,icnt_my,rcnt_my
INTEGER(4), DIMENSION(:), ALLOCATABLE::fibuf_my,idis_t,rdis_t
REAL(8), DIMENSION(:), ALLOCATABLE::frbuf_my

!%read grid data
! I0.3: 3자리 zero-padding, np>999 이면 자릿수 자동 확장 — 쓰기측
! 6_subdomain_infor_mg.f90 의 파일명 생성과 반드시 동일 포맷 유지
WRITE(fout,'(A,I0.3,A)') 'MG_tmp/part', myrank+1, '.out'

! update optimized nlevel_N

!DEC$IF defined (mpi_flag)
    CALL MPI_BCAST(nlevel_N,1,mpi_INTEGER,0,mpi_comm_world,ierr)  
    CALL MPI_BCAST(nlv_glo,1,mpi_INTEGER,0,mpi_comm_world,ierr)
    CALL MPI_BCAST(mxnbne_mg,1,mpi_INTEGER,0,mpi_comm_world,ierr)
!/
    IF(myrank /= 0) THEN
        allocate(isend_m(ndom),irecv_m(ndom))
    ENDIF
!
    CALL MPI_BCAST(isend_m,ndom,mpi_INTEGER,0,mpi_comm_world,ierr)
    CALL MPI_BCAST(irecv_m,ndom,mpi_INTEGER,0,mpi_comm_world,ierr)
!DEC$ENDIF
    nlevel = nlevel_N
    
    mxnbne_mg = max(1,mxnbne_mg)    
    ALLOCATE(id(mxnbne_mg))
!/ 
    DO i = 1,ndom
        IF(myrank == i-1) THEN
        nnsend_m = isend_m(i)
        nnrecv_m = irecv_m(i)
        EXIT
        ENDIF
    ENDDO

    DEALLOCATE(isend_m, irecv_m)
    
        nnsend_m = max(nnsend_m,1)
        nnrecv_m = max(nnrecv_m,1)
!/
! 1: for the finest level: - - - - - - - - - - - - - - - 

IF(isetup_comm.EQ.1) THEN
! C011-3 통신 모드: 카운트 SCATTER → 페이로드 SCATTERV → 파일 READ 순서 그대로
! 커서 unpack (part###.out 미생성). REAL(8) 은 pack 시점에 rt_ascii 라운딩 완료
   IF(myrank.NE.0) THEN
      IF(ALLOCATED(stg_ficnt)) DEALLOCATE(stg_ficnt,stg_frcnt)
      ALLOCATE(stg_ficnt(1),stg_frcnt(1))          ! 비루트 더미 (MPI 송신 인자용)
   ENDIF
   IF(ndom.EQ.1) THEN
      icnt_my = stg_ficnt(1)
      rcnt_my = stg_frcnt(1)
   ENDIF
!DEC$IF defined (mpi_flag)
   IF(ndom.GT.1) THEN
   CALL MPI_SCATTER(stg_ficnt,1,mpi_INTEGER,icnt_my,1,mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_SCATTER(stg_frcnt,1,mpi_INTEGER,rcnt_my,1,mpi_INTEGER,0,mpi_comm_world,ierr)
   ENDIF
!DEC$ENDIF
   ALLOCATE(fibuf_my(icnt_my),frbuf_my(MAX(1,rcnt_my)))
   IF(ndom.EQ.1) THEN
      fibuf_my(1:icnt_my) = stg_fibuf(1:icnt_my)
      frbuf_my(1:rcnt_my) = stg_frbuf(1:rcnt_my)
   ENDIF
!DEC$IF defined (mpi_flag)
   IF(ndom.GT.1) THEN
      IF(myrank.EQ.0) THEN
         ALLOCATE(idis_t(ndom),rdis_t(ndom))
         idis_t(1) = 0
         rdis_t(1) = 0
         DO k=2,ndom
            idis_t(k) = idis_t(k-1) + stg_ficnt(k-1)
            rdis_t(k) = rdis_t(k-1) + stg_frcnt(k-1)
         ENDDO
      ELSE
         ALLOCATE(idis_t(1),rdis_t(1))
         IF(.NOT.ALLOCATED(stg_fibuf)) ALLOCATE(stg_fibuf(1),stg_frbuf(1))   ! 비루트 더미
      ENDIF
      CALL MPI_SCATTERV(stg_fibuf,stg_ficnt,idis_t,mpi_INTEGER,             &
                        fibuf_my,icnt_my,mpi_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_SCATTERV(stg_frbuf,stg_frcnt,rdis_t,MPI_DOUBLE_PRECISION,    &
                        frbuf_my,rcnt_my,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
      DEALLOCATE(idis_t,rdis_t)
   ENDIF
!DEC$ENDIF
   DEALLOCATE(stg_fibuf,stg_frbuf,stg_ficnt,stg_frcnt)
   nelem   = fibuf_my(1)
   nintr   = fibuf_my(2)
   nintf   = fibuf_my(3)
   nnode   = fibuf_my(4)
   nnbd    = fibuf_my(5)
   nnodegl = fibuf_my(6)
   nnbdA   = fibuf_my(7)
   nnbdR   = fibuf_my(8)
   kci = 8
   kcr = 0
ELSE
! NEWUNIT: 런타임이 미사용 unit 을 배정 — 고정 unit(999 등)과의 충돌 원천 차단
open(newunit=iu,file=fout,status='old',action='read',iostat=alstatus)
IF(alstatus/=0) THEN
WRITE(*,*)'read_mesh_MPI: cannot open ',TRIM(fout),' rank',myrank
STOP
ENDIF
read(iu,*)nelem,nintr,nintf,nnode,nnbd,nnodegl
! NEW
read(iu,*)nnbdA,nnbdR
ENDIF
!----------------------------------------------------
IF(nelem.NE.nnode) THEN
    WRITE(*,*)'PMG error: nelem=/nnode'
    STOP
ENDIF

! 
nelemgl = nnodegl
!%read connectivity

ALLOCATE(num_neigh(nelem),e_neigh(nf_max,nelem),stat=alstatus)
IF (alstatus/=0) STOP 'not enough inod memory'
IF(isetup_comm.EQ.1) THEN
DO ie=1,nelem
   j = fibuf_my(kci+1)
   DO k=1,j
      e_neigh(k,ie) = fibuf_my(kci+1+k)
   ENDDO
   kci = kci + 1 + j
   num_neigh(ie) = j
ENDDO
ELSE
DO ie=1,nelem
   READ(iu,*) j,e_neigh(1:j,ie)
   num_neigh(ie) = j
ENDDO
ENDIF
!
!%read coordinates
ALLOCATE(coord(ndim,nnodegl),stat=alstatus)
IF(alstatus/=0) STOP 'not enough connect memory'
IF(isetup_comm.EQ.1) THEN
DO ip=1,nnodegl
   coord(1:ndim,ip) = frbuf_my(kcr+1:kcr+ndim)
   kcr = kcr + ndim
ENDDO
ELSE
DO ip=1,nnodegl
   READ(iu,*) coord(1:ndim,ip)
ENDDO
ENDIF
!--------------------------------------------
!%read neighboring data
i = MAX(nnbd,1)
ALLOCATE(rpt(i+1),spt(i+1))
ALLOCATE(nbdom(i))
nbdom = 0
rpt = 0
spt = 0

IF(nnbd.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
nbdom(1:nnbd)  = fibuf_my(kci+1:kci+nnbd)
kci = kci + nnbd
rpt(1:nnbd+1)  = fibuf_my(kci+1:kci+nnbd+1)
kci = kci + nnbd+1
spt(1:nnbd+1)  = fibuf_my(kci+1:kci+nnbd+1)
kci = kci + nnbd+1
ELSE
READ(iu,*) nbdom(1:nnbd)
READ(iu,*) rpt(1:nnbd+1)
READ(iu,*) spt(1:nnbd+1)
ENDIF
ENDIF
!
i = MAX(1,nnode-nintf)
j = MAX(1,spt(nnbd+1)-1)

ALLOCATE(rintf(i))
ALLOCATE(sintf(j))

rintf = 0
sintf = 0
!
IF(nnbd.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
! 파일 모드와 동일 계약: rint 개수 = rpt(nnbd+1)-1 = nnode-nintf (writer 측 ri 와 일치)
rintf(1:(nnode-nintf))   = fibuf_my(kci+1:kci+(nnode-nintf))
kci = kci + (nnode-nintf)
sintf(1:(spt(nnbd+1)-1)) = fibuf_my(kci+1:kci+(spt(nnbd+1)-1))
kci = kci + (spt(nnbd+1)-1)
ELSE
READ(iu,*) rintf(1:(nnode-nintf))
READ(iu,*) sintf(1:(spt(nnbd+1)-1))
ENDIF
ENDIF
! - - - - - - - - - NEW for A  - - - - - - - - - 
i = MAX(nnbdA,1)

ALLOCATE(nbdomA(i))
ALLOCATE(rptA(i+1),sptA(i+1))
nbdomA = 0
rptA = 0
sptA = 0
IF(nnbdA.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
nbdomA(1:nnbdA)  = fibuf_my(kci+1:kci+nnbdA)
kci = kci + nnbdA
rptA(1:nnbdA+1)  = fibuf_my(kci+1:kci+nnbdA+1)
kci = kci + nnbdA+1
sptA(1:nnbdA+1)  = fibuf_my(kci+1:kci+nnbdA+1)
kci = kci + nnbdA+1
ELSE
READ(iu,*) nbdomA(1:nnbdA)
READ(iu,*) rptA(1:nnbdA+1)
READ(iu,*) sptA(1:nnbdA+1)
ENDIF
ENDIF
!
i = MAX(1,rptA(nnbdA+1)-1)
j = MAX(1,sptA(nnbdA+1)-1)

ALLOCATE(rintfA(i))
ALLOCATE(sintfA(j))

rintfA = 0
sintfA = 0
!
IF(nnbdA.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
rintfA(1:(rptA(nnbdA+1)-1)) = fibuf_my(kci+1:kci+(rptA(nnbdA+1)-1))
kci = kci + (rptA(nnbdA+1)-1)
sintfA(1:(sptA(nnbdA+1)-1)) = fibuf_my(kci+1:kci+(sptA(nnbdA+1)-1))
kci = kci + (sptA(nnbdA+1)-1)
ELSE
READ(iu,*) rintfA(1:(rptA(nnbdA+1)-1))
READ(iu,*) sintfA(1:(sptA(nnbdA+1)-1))
ENDIF
ENDIF
! - - - - - - - - - NEW for R  - - - - - - - - - 
i = MAX(nnbdR,1)

ALLOCATE(nbdomR(i))
ALLOCATE(rptR(i+1),sptR(i+1))
nbdomR = 0
rptR = 0
sptR = 0
IF(nnbdR.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
nbdomR(1:nnbdR)  = fibuf_my(kci+1:kci+nnbdR)
kci = kci + nnbdR
rptR(1:nnbdR+1)  = fibuf_my(kci+1:kci+nnbdR+1)
kci = kci + nnbdR+1
sptR(1:nnbdR+1)  = fibuf_my(kci+1:kci+nnbdR+1)
kci = kci + nnbdR+1
ELSE
READ(iu,*) nbdomR(1:nnbdR)
READ(iu,*) rptR(1:nnbdR+1)
READ(iu,*) sptR(1:nnbdR+1)
ENDIF
ENDIF
!
i = MAX(1,rptR(nnbdR+1)-1)
j = MAX(1,sptR(nnbdR+1)-1)

ALLOCATE(rintfR(i))
ALLOCATE(sintfR(j))

rintfR = 0
sintfR = 0
!
IF(nnbdR.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
rintfR(1:(rptR(nnbdR+1)-1)) = fibuf_my(kci+1:kci+(rptR(nnbdR+1)-1))
kci = kci + (rptR(nnbdR+1)-1)
sintfR(1:(sptR(nnbdR+1)-1)) = fibuf_my(kci+1:kci+(sptR(nnbdR+1)-1))
kci = kci + (sptR(nnbdR+1)-1)
ELSE
READ(iu,*) rintfR(1:(rptR(nnbdR+1)-1))
READ(iu,*) sintfR(1:(sptR(nnbdR+1)-1))
ENDIF
ENDIF
!
IF(isetup_comm.EQ.1) THEN
! unpack 정합 검사: 커서가 정확히 수신 길이에서 끝나야 함 (프로토콜 자기 검증)
IF(kci.NE.icnt_my .OR. kcr.NE.rcnt_my) THEN
WRITE(*,*)'read_mesh_MPI: finest unpack mismatch rank',myrank,kci,icnt_my,kcr,rcnt_my
STOP
ENDIF
DEALLOCATE(fibuf_my,frbuf_my)
ELSE
CLOSE(iu)
ENDIF

!/ delete the tmp. file
!/

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! 2: for the coarse levels: - - - - - - - - - - - - - - - 
IF(isetup_comm.EQ.1) THEN
! 통신 모드 (LOOP C011-2): rank0 스테이징(= PMG_infor 파일 내용과 1:1, 전부 정수라
! 라운딩 무관)을 전 rank 로 BCAST 하고 아래에서 자기 열(myrank+1)을 추출.
! 파일 모드의 O(np) 더미 스킵과 공유 파일 동시 OPEN 이 소거된다.
   IF(myrank.NE.0) THEN
      IF(ALLOCATED(stg_iintf)) DEALLOCATE(stg_iintf,stg_inodegl,stg_inbdc,   &
                                          stg_ialvP,stg_inmax,stg_nnzc0,     &
                                          stg_nnzi,stg_nnzr)
      ALLOCATE(stg_iintf(nlevel,ndom),stg_inodegl(nlevel,ndom),              &
               stg_inbdc(nlevel,ndom),stg_ialvP(nlevel+1,ndom),              &
               stg_inmax(nlevel),stg_nnzc0(ndom),stg_nnzi(ndom),stg_nnzr(ndom))
   ENDIF
!DEC$IF defined (mpi_flag)
   IF(ndom.GT.1) THEN
   CALL MPI_BCAST(stg_iintf,  nlevel*ndom,    mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_inodegl,nlevel*ndom,    mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_inbdc,  nlevel*ndom,    mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_ialvP,  (nlevel+1)*ndom,mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_inmax,  nlevel,         mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_nnzc0,  ndom,           mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_nnzi,   ndom,           mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_BCAST(stg_nnzr,   ndom,           mpi_INTEGER,0,mpi_comm_world,ierr)
   ENDIF
!DEC$ENDIF
ELSE
OPEN(newunit=iu,file='MG_tmp/PMG_infor',status='old',action='read',iostat=alstatus)
IF(alstatus/=0) THEN
WRITE(*,*)'read_mesh_MPI: cannot open MG_tmp/PMG_infor rank',myrank
STOP
ENDIF

ntmp = nlevel*2+2
ntmp = ntmp*myrank

DO i = 1,ntmp
READ(iu,*) j
ENDDO
ENDIF

! from MPI-MG
ALLOCATE(iintf(nlevel),inodegl(nlevel),inbdc(nlevel))
ALLOCATE(ibdomc(ndom,nlevel),isptc(ndom,nlevel),irptc(ndom,nlevel))
ALLOCATE(inmax(nlevel))
!
ibdomc = 0
isptc = 0
irptc = 0
inmax = 0
!
IF(isetup_comm.EQ.1) THEN
DO ilv=1,nlevel
iintf(ilv)   = stg_iintf(ilv,myrank+1)
inodegl(ilv) = stg_inodegl(ilv,myrank+1)
inbdc(ilv)   = stg_inbdc(ilv,myrank+1)
inmax(ilv)   = stg_inmax(ilv)
ENDDO
ELSE
DO ilv=1,nlevel
READ(iu,*) iintf(ilv),inodegl(ilv),inbdc(ilv),inmax(ilv)    ! reading *
ENDDO
ENDIF





!

! from MG_coord
ALLOCATE(ialv(nlevel+1))

IF(isetup_comm.EQ.1) THEN
ialv(1:nlevel+1) = stg_ialvP(1:nlevel+1,myrank+1)
ELSE
DO ilv = 1,nlevel+1
    READ(iu,*) ialv(ilv)                               ! reading *
ENDDO
ENDIF

ncolf = ialv(nlevel+1)-1
ncolc = ncolf -nnode
ntmp = SUM(inodegl(2:nlevel))
ALLOCATE(coordc(ndim,ntmp))

! From MG_matrix
IF(isetup_comm.EQ.1) THEN
nnzc0 = stg_nnzc0(myrank+1)
nnzi  = stg_nnzi(myrank+1)
nnzr  = stg_nnzr(myrank+1)
! 스테이징 사용 완료 — 전 rank 해제 (재진입 대비 가드는 할당측에 있음)
DEALLOCATE(stg_iintf,stg_inodegl,stg_inbdc,stg_ialvP,stg_inmax,             &
           stg_nnzc0,stg_nnzi,stg_nnzr)
ELSE
READ(iu,*) nnzc0,nnzi,nnzr                          ! reading *

CLOSE(iu)
ENDIF

DO ilv=2,nlevel
    j = ialv(ilv+1)-ialv(ilv)-iintf(ilv)
    IF(j.GT.nnsend_m) THEN
     nnsend_m = j
    endif
    
     IF(j.GT.nnrecv_m) THEN
     nnsend_m = j
     endif  
ENDDO

    
ALLOCATE(isintfc(nnsend_m,nlevel),irintfc(nnrecv_m,nlevel))
ALLOCATE(isiac(nnsend_m+1,nlevel),iriac(nnrecv_m+1,nlevel))
isintfc = 0
irintfc = 0
isiac = 0
iriac = 0

ALLOCATE(iac(ncolc+1),juc(ncolc),jac(nnzc0),auc(nnzc0))
ALLOCATE(iai(ncolf+1),jai(nnzi),Xintp(nnzi))
ALLOCATE(iar(ncolc+1),jar(nnzr),Xrest(nnzr))
iac(1) = 1
iai(1) = 1
iar(1) = 1
!
ALLOCATE(diagrc(ncolc))
diagrc = 0.d0

! local
ALLOCATE(nbdom1(ndom),rpt1(ndom),spt1(ndom))
i = maxval(inodegl(2:nlevel))
ALLOCATE(coord1(ndim,i))
ALLOCATE(rintf1(nnrecv_m),sintf1(nnsend_m))

! iai1/iar1 은 레벨 루프의 재사용 버퍼 — 행 수가 레벨별 nnode1(=ialv 차분)이므로
! 최대 레벨 폭으로 할당해야 함. nnode(fine)로 잡으면 극소 도메인(대규모 np)에서
! coarse 레벨 폭 > nnode 가 되어 오버런 → 힙 오염 (np=900 비결정 크래시, LOG C010-3)
i = MAXVAL(ialv(2:nlevel+1)-ialv(1:nlevel))
IF(i.LT.nnode) i = nnode
ALLOCATE(iai1(i+1),jai1(nnzi),Xintp1(nnzi))
ALLOCATE(iar1(i+1),jar1(nnzr),Xrest1(nnzr))
!/
!/
ALLOCATE(ia1(nnode+1),ja1(nnzc0),ju1(nnode),au1(nnzc0))

! initial and set for ilv = 1
ntmpc = 0
IF(nnbd.NE.0) THEN
ibdomc(1:nnbd,1) = nbdom(1:nnbd)
isptc(1:nnbd+1,1) = spt(1:nnbd+1)
irptc(1:nnbd+1,1) = rpt(1:nnbd+1)
isintfc(1:(spt(nnbd+1)-1),1) = sintf(1:(spt(nnbd+1)-1))
irintfc(1:(nnode-nintf),1) = rintf(1:(nnode-nintf))
ENDIF
! NEW for A, R P ----------------------
nnsend_mA = nnsend_m
nnrecv_mA = nnrecv_m
nnsend_mR = nnsend_m
nnrecv_mR = nnrecv_m
nnsend_mP = nnsend_m
nnrecv_mP = nnrecv_m
ALLOCATE(inbdcA(nlevel),inbdcR(nlevel),inbdcP(nlevel))
ALLOCATE(ibdomcA(ndom,nlevel),isptcA(ndom,nlevel),irptcA(ndom,nlevel))
ALLOCATE(ibdomcR(ndom,nlevel),isptcR(ndom,nlevel),irptcR(ndom,nlevel))
ALLOCATE(ibdomcP(ndom,nlevel),isptcP(ndom,nlevel),irptcP(ndom,nlevel))
ALLOCATE(isintfcA(nnsend_mA,nlevel),irintfcA(nnrecv_mA,nlevel))
ALLOCATE(isintfcR(nnsend_mR,nlevel),irintfcR(nnrecv_mR,nlevel))
ALLOCATE(isintfcP(nnsend_mP,nlevel),irintfcP(nnrecv_mP,nlevel))
! - - - - - - - - - - - - - - - - - - - 
inbdcA = 0
inbdcR = 0
inbdcP = 0
ibdomcA = 0
isptcA = 0
irptcA = 0
ibdomcR = 0
isptcR = 0
irptcR = 0
isintfcA = 0
irintfcA = 0
isintfcR = 0
irintfcR = 0
isintfcP = 0
irintfcP = 0
! - - - - - - - - - - - - - - 
IF(isetup_comm.EQ.1) THEN
! C011-4 통신 모드: rank0 의 prc별 성장형 스트림(stg_mg)을 연접 → 카운트 SCATTER
! → 페이로드 SCATTERV → 파일 READ 순서 그대로 커서 unpack (part_MG 미생성)
   IF(myrank.EQ.0) THEN
      IF(ALLOCATED(stg_ficnt)) DEALLOCATE(stg_ficnt,stg_frcnt)
      ALLOCATE(stg_ficnt(ndom),stg_frcnt(ndom))
      DO k=1,ndom
         stg_ficnt(k) = stg_mg(k)%ni
         stg_frcnt(k) = stg_mg(k)%nr
      ENDDO
      IF(ALLOCATED(stg_fibuf)) DEALLOCATE(stg_fibuf,stg_frbuf)
      ALLOCATE(stg_fibuf(SUM(stg_ficnt)),stg_frbuf(MAX(1,SUM(stg_frcnt))))
      i1 = 0
      i2 = 0
      DO k=1,ndom
         stg_fibuf(i1+1:i1+stg_ficnt(k)) = stg_mg(k)%ib(1:stg_ficnt(k))
         stg_frbuf(i2+1:i2+stg_frcnt(k)) = stg_mg(k)%rb(1:stg_frcnt(k))
         i1 = i1 + stg_ficnt(k)
         i2 = i2 + stg_frcnt(k)
      ENDDO
      DEALLOCATE(stg_mg)
   ELSE
      IF(ALLOCATED(stg_ficnt)) DEALLOCATE(stg_ficnt,stg_frcnt)
      ALLOCATE(stg_ficnt(1),stg_frcnt(1))
   ENDIF
   IF(ndom.EQ.1) THEN
      icnt_my = stg_ficnt(1)
      rcnt_my = stg_frcnt(1)
   ENDIF
!DEC$IF defined (mpi_flag)
   IF(ndom.GT.1) THEN
   CALL MPI_SCATTER(stg_ficnt,1,mpi_INTEGER,icnt_my,1,mpi_INTEGER,0,mpi_comm_world,ierr)
   CALL MPI_SCATTER(stg_frcnt,1,mpi_INTEGER,rcnt_my,1,mpi_INTEGER,0,mpi_comm_world,ierr)
   ENDIF
!DEC$ENDIF
   ALLOCATE(fibuf_my(icnt_my),frbuf_my(MAX(1,rcnt_my)))
   IF(ndom.EQ.1) THEN
      fibuf_my(1:icnt_my) = stg_fibuf(1:icnt_my)
      frbuf_my(1:rcnt_my) = stg_frbuf(1:rcnt_my)
   ENDIF
!DEC$IF defined (mpi_flag)
   IF(ndom.GT.1) THEN
      IF(myrank.EQ.0) THEN
         ALLOCATE(idis_t(ndom),rdis_t(ndom))
         idis_t(1) = 0
         rdis_t(1) = 0
         DO k=2,ndom
            idis_t(k) = idis_t(k-1) + stg_ficnt(k-1)
            rdis_t(k) = rdis_t(k-1) + stg_frcnt(k-1)
         ENDDO
      ELSE
         ALLOCATE(idis_t(1),rdis_t(1))
         IF(.NOT.ALLOCATED(stg_fibuf)) ALLOCATE(stg_fibuf(1),stg_frbuf(1))   ! 비루트 더미
      ENDIF
      CALL MPI_SCATTERV(stg_fibuf,stg_ficnt,idis_t,mpi_INTEGER,             &
                        fibuf_my,icnt_my,mpi_INTEGER,0,mpi_comm_world,ierr)
      CALL MPI_SCATTERV(stg_frbuf,stg_frcnt,rdis_t,MPI_DOUBLE_PRECISION,    &
                        frbuf_my,rcnt_my,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)
      DEALLOCATE(idis_t,rdis_t)
   ENDIF
!DEC$ENDIF
   DEALLOCATE(stg_fibuf,stg_frbuf,stg_ficnt,stg_frcnt)
   kci = 0
   kcr = 0
ELSE
WRITE(fout,'(A,I0.3,A)') 'MG_tmp/part_MG', myrank+1, '.out'
OPEN(newunit=iu,file=fout,status='old',action='read',iostat=alstatus)
IF(alstatus/=0) THEN
WRITE(*,*)'read_mesh_MPI: cannot open ',TRIM(fout),' rank',myrank
STOP
ENDIF
ENDIF
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
DO ilv = 2,nlevel

IF(isetup_comm.EQ.0) READ(iu,*)                       ! reading ('coarse' 마커)
nintf1 = iintf(ilv)
nnode1 = ialv(ilv+1)-ialv(ilv)
nnbd1 = inbdc(ilv)
nnode1gl = inodegl(ilv)
nelem1 = nnode1
nnode0 = ialv(ilv)-ialv(ilv-1)

! 1:reading local

IF(isetup_comm.EQ.1) THEN
i1 = fibuf_my(kci+1)
i2 = fibuf_my(kci+2)
i3 = fibuf_my(kci+3)
i4 = fibuf_my(kci+4)
kci = kci + 4
ELSE
READ(iu,*) i1,i2,i3,i4
ENDIF
IF((i1.NE.nintf1).OR.(i2.NE.nnode1).OR.(i3.NE.nnbd1).OR.(i4.NE.nnode1gl)) THEN
WRITE(*,*)'error nintf1'
STOP
ENDIF
! NEW for A and R, P - - - - - -
IF(isetup_comm.EQ.1) THEN
IF(ilv.NE.nlevel) THEN
inbdcA(ilv) = fibuf_my(kci+1)
inbdcR(ilv) = fibuf_my(kci+2)
inbdcP(ilv) = fibuf_my(kci+3)
kci = kci + 3
ELSE
inbdcA(ilv) = fibuf_my(kci+1)
inbdcP(ilv) = fibuf_my(kci+2)
kci = kci + 2
ENDIF
ELSE
IF(ilv.NE.nlevel) THEN
READ(iu,*) inbdcA(ilv),inbdcR(ilv),inbdcP(ilv)
ELSE
READ(iu,*) inbdcA(ilv),inbdcP(ilv)
ENDIF
ENDIF
! - - - - - - - - - - - - - - -
IF(isetup_comm.EQ.1) THEN
DO i=1,nnode1gl
   coord1(1:ndim,i) = frbuf_my(kcr+1:kcr+ndim)
   kcr = kcr + ndim
ENDDO
ELSE
DO i=1,nnode1gl
   READ(iu,*) coord1(1:ndim,i)                        ! reading
ENDDO
ENDIF

IF(nnbd1.NE.0) THEN
IF(isetup_comm.EQ.1) THEN
nbdom1(1:nnbd1) = fibuf_my(kci+1:kci+nnbd1)
kci = kci + nnbd1
rpt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
spt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
rintf1(1:(nnode1-nintf1))   = fibuf_my(kci+1:kci+(nnode1-nintf1))
kci = kci + (nnode1-nintf1)
sintf1(1:(spt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(spt1(nnbd1+1)-1))
kci = kci + (spt1(nnbd1+1)-1)
ELSE
READ(iu,*) nbdom1(1:nnbd1)                            ! reading
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(nnode1-nintf1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))
ENDIF

IF((nnode1-nintf1).NE.(rpt1(nnbd1+1)-1)) THEN
WRITE(*,*)'error/nnode1-nintf1.NE.rpt1(nnbd1+1)-1'
STOP
ENDIF

ENDIF

! 2: adding to global
  coordc(1:ndim,ntmpc+1:ntmpc+nnode1gl) = coord1(1:ndim,1:nnode1gl)
  ntmpc = ntmpc + nnode1gl

IF(nnbd1.NE.0) THEN
  ibdomc(1:nnbd1,ilv) = nbdom1(1:nnbd1)
  
  isptc(1:nnbd1+1,ilv) = spt1(1:nnbd1+1)
  irptc(1:nnbd1+1,ilv) = rpt1(1:nnbd1+1)
  
  isintfc(1:(spt1(nnbd1+1)-1),ilv) = sintf1(1:(spt1(nnbd1+1)-1))
  irintfc(1:(nnode1-nintf1),ilv) = rintf1(1:(nnode1-nintf1)) 

ENDIF

! NEW for A - - - - - - - - - - - - - -
 nnbd1 =  inbdcA(ilv)
IF(nnbd1.NE.0) THEN

IF(isetup_comm.EQ.1) THEN
nbdom1(1:nnbd1) = fibuf_my(kci+1:kci+nnbd1)
kci = kci + nnbd1
rpt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
spt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
rintf1(1:(rpt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(rpt1(nnbd1+1)-1))
kci = kci + (rpt1(nnbd1+1)-1)
sintf1(1:(spt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(spt1(nnbd1+1)-1))
kci = kci + (spt1(nnbd1+1)-1)
ELSE
READ(iu,*) nbdom1(1:nnbd1)                            ! reading
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))
ENDIF

ENDIF

! adding to global A
IF(nnbd1.NE.0) THEN
  ibdomcA(1:nnbd1,ilv) = nbdom1(1:nnbd1)
  
  isptcA(1:nnbd1+1,ilv) = spt1(1:nnbd1+1)
  irptcA(1:nnbd1+1,ilv) = rpt1(1:nnbd1+1)
  
  isintfcA(1:(spt1(nnbd1+1)-1),ilv) = sintf1(1:(spt1(nnbd1+1)-1))
  irintfcA(1:(rpt1(nnbd1+1)-1),ilv) = rintf1(1:(rpt1(nnbd1+1)-1))
ENDIF 

! NEW for R - - - - - - - - - - - - - -
IF(ilv.NE.nlevel) THEN
 nnbd1 =  inbdcR(ilv)
IF(nnbd1.NE.0) THEN

IF(isetup_comm.EQ.1) THEN
nbdom1(1:nnbd1) = fibuf_my(kci+1:kci+nnbd1)
kci = kci + nnbd1
rpt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
spt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
rintf1(1:(rpt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(rpt1(nnbd1+1)-1))
kci = kci + (rpt1(nnbd1+1)-1)
sintf1(1:(spt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(spt1(nnbd1+1)-1))
kci = kci + (spt1(nnbd1+1)-1)
ELSE
READ(iu,*) nbdom1(1:nnbd1)                            ! reading
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))
ENDIF

ENDIF

! adding to global R
IF(nnbd1.NE.0) THEN
  ibdomcR(1:nnbd1,ilv) = nbdom1(1:nnbd1)
  
  isptcR(1:nnbd1+1,ilv) = spt1(1:nnbd1+1)
  irptcR(1:nnbd1+1,ilv) = rpt1(1:nnbd1+1)
  
  isintfcR(1:(spt1(nnbd1+1)-1),ilv) = sintf1(1:(spt1(nnbd1+1)-1))
  irintfcR(1:(rpt1(nnbd1+1)-1),ilv) = rintf1(1:(rpt1(nnbd1+1)-1))
ENDIF 
  
ENDIF

! NEW for P - - - - - - - - - - - - - -
 nnbd1 =  inbdcP(ilv)
IF(nnbd1.NE.0) THEN

IF(isetup_comm.EQ.1) THEN
nbdom1(1:nnbd1) = fibuf_my(kci+1:kci+nnbd1)
kci = kci + nnbd1
rpt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
spt1(1:nnbd1+1) = fibuf_my(kci+1:kci+nnbd1+1)
kci = kci + nnbd1+1
rintf1(1:(rpt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(rpt1(nnbd1+1)-1))
kci = kci + (rpt1(nnbd1+1)-1)
sintf1(1:(spt1(nnbd1+1)-1)) = fibuf_my(kci+1:kci+(spt1(nnbd1+1)-1))
kci = kci + (spt1(nnbd1+1)-1)
ELSE
READ(iu,*) nbdom1(1:nnbd1)                            ! reading
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))
ENDIF

ENDIF

! adding to global P
IF(nnbd1.NE.0) THEN
  ibdomcP(1:nnbd1,ilv) = nbdom1(1:nnbd1)
  
  isptcP(1:nnbd1+1,ilv) = spt1(1:nnbd1+1)
  irptcP(1:nnbd1+1,ilv) = rpt1(1:nnbd1+1)
  
  isintfcP(1:(spt1(nnbd1+1)-1),ilv) = sintf1(1:(spt1(nnbd1+1)-1))
  irintfcP(1:(rpt1(nnbd1+1)-1),ilv) = rintf1(1:(rpt1(nnbd1+1)-1))
ENDIF 

! FOR P,R & AC
! - - - -
! 1: reading local:
IF(isetup_comm.EQ.0) read(iu,*)                   ! reading ('P-1' 마커)

iai1(1) = 1
Do i=1,nnode0
IF(isetup_comm.EQ.1) THEN
nnd = fibuf_my(kci+1)
id(1:nnd) = fibuf_my(kci+2:kci+1+nnd)
kci = kci + 1 + nnd
tmp(1:nnd) = frbuf_my(kcr+1:kcr+nnd)
kcr = kcr + nnd
ELSE
read(iu,*) nnd,id(1:nnd)
read(iu,*) tmp(1:nnd)
ENDIF
 k=iai1(i)
 do j=1,nnd
 jai1(k)=id(j)
 Xintp1(k) = tmp(j)
  k=k+1
 enddo
 iai1(i+1)=k
enddo

! R
IF(isetup_comm.EQ.0) READ(iu,*)                   ! reading ('R-1' 마커)
iar1(1) = 1
DO i=1,nnode1
IF(isetup_comm.EQ.1) THEN
nnd = fibuf_my(kci+1)
id(1:nnd) = fibuf_my(kci+2:kci+1+nnd)
kci = kci + 1 + nnd
tmp(1:nnd) = frbuf_my(kcr+1:kcr+nnd)
kcr = kcr + nnd
ELSE
READ(iu,*) nnd,id(1:nnd)
READ(iu,*) tmp(1:nnd)
ENDIF
 k=iar1(i)
 DO j=1,nnd
 jar1(k)=id(j)
 Xrest1(k) = tmp(j)
  k=k+1
 ENDDO
 iar1(i+1)=k
ENDDO

! AC
IF(isetup_comm.EQ.0) READ(iu,*)                  ! reading ('Ac-1' 마커)
ia1(1) = 1
DO i=1,nnode1
IF(isetup_comm.EQ.1) THEN
nnd = fibuf_my(kci+1)
id(1:nnd) = fibuf_my(kci+2:kci+1+nnd)
kci = kci + 1 + nnd
ELSE
READ(iu,*) nnd,id(1:nnd)
ENDIF
k=ia1(i)
 DO j=1,nnd
 ja1(k)=id(j)
 IF(id(j).EQ.i) ju1(i) = k
  k=k+1
 ENDDO
 ia1(i+1)=k
ENDDO

! 2: adding to global
ncolf1 = ialv(ilv+1) - ialv(1)
ncolc1 = ncolf1 -nnode
ntmp = ncolc1 - nnode1
ntmpf = ialv(ilv) - ialv(1)-nnode0
ncolc2 = ialv(ilv) - ialv(1)

!P
nnzt1 = iai(ntmpf+1)-1
nnzi1 = iai1(nnode0+1)-1     !iai(ncolc2+1)-iai(ntmpf+1)

iai(ntmpf+1:ncolc2+1) = iai1(1:nnode0+1) + nnzt1
jai(nnzt1+1:nnzt1+nnzi1) = jai1 (1:nnzi1) + ntmp

Xintp(nnzt1+1:nnzt1+nnzi1) = Xintp1 (1:nnzi1)

!R
nnzt1 = iar(ntmp+1)-1
nnzr1 = iar1(nnode1+1)-1      !iar(ncolc1+1)-iar(ntmp+1)

iar(ntmp+1:ncolc1+1) = iar1(1:nnode1+1) + nnzt1
jar(nnzt1+1:nnzt1+nnzr1) = jar1 (1:nnzr1) + ntmpf

Xrest(nnzt1+1:nnzt1+nnzr1) = Xrest1 (1:nnzr1)

! AC

nnzt = iac(ntmp+1)-1
nnz1 = ia1(nnode1+1)-1     !iac(ncolc1+1)-iac(ntmp+1)
      
iac(ntmp+1:ncolc1+1) = ia1(1:nnode1+1) + nnzt
jac(nnzt+1:nnzt+nnz1) = ja1(1:nnz1) + ntmp
juc (ntmp+1:ncolc1) = ju1 (1:nnode1) + nnzt

ENDDO

! Deallocate local
DEALLOCATE(nbdom1,rpt1,spt1)
DEALLOCATE(coord1)
DEALLOCATE(rintf1,sintf1)

DEALLOCATE(iai1,jai1,Xintp1)
DEALLOCATE(iar1,jar1,Xrest1)
      
! 3:coarsest-global = = = = = = = = = = = = = = = = = = = = = = = 

  IF(n_GC.EQ.1) THEN
  IF(isetup_comm.EQ.1) THEN
  nnodeC = fibuf_my(kci+1)
  nnodeG = fibuf_my(kci+2)
  nnzG   = fibuf_my(kci+3)
  kci = kci + 3
  ELSE
  READ(iu,*)                                         ! reading ('A_GC' 마커)
  READ(iu,*) nnodeC,nnodeG,nnzG
  ENDIF
  
!/
!/
  ALLOCATE(iaG(nnodeG+1),jaG(nnzG),juG(nnodeG),auG(nnzG),auG0(nnzG),stat=alstatus)
     IF (alstatus/=0) THEN
         WRITE(*,*)'not allocated, iaG',myrank
     ENDIF  
!/
     i = nnodeC
     if (i == 0) i = 1
     
  ALLOCATE(imapG(i),eG(nnodeG),rG(nnodeG),rG0(nnodeG),stat=alstatus)
     IF (alstatus/=0) THEN
         WRITE(*,*)'not allocated, rG',myrank
     ENDIF
  ALLOCATE(coordG(ndim,nnodeG),stat=alstatus)
     IF (alstatus/=0) THEN
         WRITE(*,*)'not allocated, coordG',myrank
     ENDIF
!/  
!imapG
  IF(nnodeC .GT. 0) THEN
  IF(isetup_comm.EQ.1) THEN
  imapG(1:nnodeC) = fibuf_my(kci+1:kci+nnodeC)
  kci = kci + nnodeC
  ELSE
  DO i=1,nnodeC
   READ(iu,*) imapG(i)
  ENDDO
  ENDIF
  ENDIF


! AC_G

 iaG(1) = 1
 DO i=1,nnodeG
  IF(isetup_comm.EQ.1) THEN
   nnd = fibuf_my(kci+1)
   id(1:nnd) = fibuf_my(kci+2:kci+1+nnd)
   kci = kci + 1 + nnd
   coordG(1:ndim,i) = frbuf_my(kcr+1:kcr+ndim)
   kcr = kcr + ndim
  ELSE
   READ(iu,*) nnd,id(1:nnd),coordG(1:ndim,i)
  ENDIF
   k=iaG(i)
   DO j=1,nnd
    jaG(k)=id(j)
    IF(id(j).EQ.i) juG(i) = k
    k=k+1
   ENDDO
 
   iaG(i+1)=k
 
  ENDDO

   IF(iaG(nnodeG+1).GT.(nnzG+1))THEN
    WRITE(*,*)'PMG error in read AcG'
   ENDIF

  ENDIF

IF(isetup_comm.EQ.1) THEN
! unpack 정합 검사 (coarse 스트림 — 프로토콜 자기 검증)
IF(kci.NE.icnt_my .OR. kcr.NE.rcnt_my) THEN
WRITE(*,*)'read_mesh_MPI: coarse unpack mismatch rank',myrank,kci,icnt_my,kcr,rcnt_my
STOP
ENDIF
DEALLOCATE(fibuf_my,frbuf_my)
ELSE
CLOSE(iu)
ENDIF
!/
     nnods = ialv(nlevel+1)-ialv(nlevel)
     ncolf = ialv(nlevel+1)-ialv(1)
     ncolc = ncolf - nelem
     nintfs = iintf(nlevel)
     
     nnzs = ia1(nnods+1)-1
           
     ALLOCATE(r(nelem),rt(ncolf),rc(ncolc))
     ALLOCATE(rs(nnods),es(nnods),e(ncolc),et(ncolf))
     
     nnzs = ia1(nnods+1)-1
     ALLOCATE(aus(nnzs))
     ALLOCATE(ias(nnods+1),jas(nnzs),jus(nnods)) 
          ias = ia1
          jas = ja1
          jus = ju1
!/
DEALLOCATE(ia1,ja1,ju1,au1)
DEALLOCATE(id)
!------------------------

RETURN 

END SUBROUTINE



