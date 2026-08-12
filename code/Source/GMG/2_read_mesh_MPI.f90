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
USE MD_MG_index, ONLY: mxnbne,nlevel,n_GC,nlevel_N,mxnbne_mg,isend_m,irecv_m
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

! NEWUNIT: 런타임이 미사용 unit 을 배정 — 고정 unit(999 등)과의 충돌 원천 차단
open(newunit=iu,file=fout,status='old',action='read',iostat=alstatus)
IF(alstatus/=0) THEN
WRITE(*,*)'read_mesh_MPI: cannot open ',TRIM(fout),' rank',myrank
STOP
ENDIF
read(iu,*)nelem,nintr,nintf,nnode,nnbd,nnodegl
! NEW
read(iu,*)nnbdA,nnbdR
!----------------------------------------------------
IF(nelem.NE.nnode) THEN
    WRITE(999,*)'nelem=/nnode'!
    STOP
ENDIF

! 
nelemgl = nnodegl
!%read connectivity

ALLOCATE(num_neigh(nelem),e_neigh(nf_max,nelem),stat=alstatus)
IF (alstatus/=0) STOP 'not enough inod memory'
DO ie=1,nelem
   READ(iu,*) j,e_neigh(1:j,ie)
   num_neigh(ie) = j
ENDDO
!
!----------------------------------------------------------------------
!%read coordinates
ALLOCATE(coord(ndim,nnodegl),stat=alstatus)
IF(alstatus/=0) STOP 'not enough connect memory'
DO ip=1,nnodegl
   READ(iu,*) coord(1:ndim,ip)
ENDDO
!--------------------------------------------
!%read neighboring data
i = MAX(nnbd,1)
ALLOCATE(rpt(i+1),spt(i+1))
ALLOCATE(nbdom(i))
nbdom = 0
rpt = 0
spt = 0

IF(nnbd.NE.0) THEN
READ(iu,*) nbdom(1:nnbd)
READ(iu,*) rpt(1:nnbd+1)
READ(iu,*) spt(1:nnbd+1)
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
    
READ(iu,*) rintf(1:(nnode-nintf))
READ(iu,*) sintf(1:(spt(nnbd+1)-1))
ENDIF
! - - - - - - - - - NEW for A  - - - - - - - - - 
i = MAX(nnbdA,1)

ALLOCATE(nbdomA(i))
ALLOCATE(rptA(i+1),sptA(i+1))
nbdomA = 0
rptA = 0
sptA = 0
IF(nnbdA.NE.0) THEN
READ(iu,*) nbdomA(1:nnbdA)
READ(iu,*) rptA(1:nnbdA+1)
READ(iu,*) sptA(1:nnbdA+1)
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

READ(iu,*) rintfA(1:(rptA(nnbdA+1)-1))
READ(iu,*) sintfA(1:(sptA(nnbdA+1)-1))

ENDIF
! - - - - - - - - - NEW for R  - - - - - - - - - 
i = MAX(nnbdR,1)

ALLOCATE(nbdomR(i))
ALLOCATE(rptR(i+1),sptR(i+1))
nbdomR = 0
rptR = 0
sptR = 0
IF(nnbdR.NE.0) THEN
READ(iu,*) nbdomR(1:nnbdR)
READ(iu,*) rptR(1:nnbdR+1)
READ(iu,*) sptR(1:nnbdR+1)
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

READ(iu,*) rintfR(1:(rptR(nnbdR+1)-1))
READ(iu,*) sintfR(1:(sptR(nnbdR+1)-1))

ENDIF
! 
CLOSE(iu)

!/ delete the tmp. file
!call system('del fout')
!/

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
! 2: for the coarse levels: - - - - - - - - - - - - - - - 
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
DO ilv=1,nlevel
READ(iu,*) iintf(ilv),inodegl(ilv),inbdc(ilv),inmax(ilv)    ! reading *
ENDDO

!IF(nnbd.EQ.0) THEN
!  nnsend_m = 0
!  nnrecv_m = 0
!ELSE
!  nnsend_m = spt(nnbd+1)-1                                    ! using data from finest level.
!  nnrecv_m = rpt(nnbd+1)-1
!ENDIF


!IF((nnode-nintf).NE.nnrecv_m) THEN
!WRITE(*,*)'error/nnode-nintf.NE.nnrecv_m'
!STOP
!ENDIF


!nnsend_m = 2*nnsend_m                                    ! using data from finest level.
!nnrecv_m = 2*nnrecv_m
!
!nnsend_m = MAX(1,nnsend_m)                                    ! using data from finest level.
!nnrecv_m = MAX(1,nnrecv_m)

! from MG_coord
ALLOCATE(ialv(nlevel+1))

DO ilv = 1,nlevel+1
    READ(iu,*) ialv(ilv)                               ! reading *
ENDDO

ncolf = ialv(nlevel+1)-1
ncolc = ncolf -nnode
ntmp = SUM(inodegl(2:nlevel))
ALLOCATE(coordc(ndim,ntmp))

! From MG_matrix
READ(iu,*) nnzc0,nnzi,nnzr                          ! reading *
    
CLOSE(iu)

DO ilv=2,nlevel
    j = ialv(ilv+1)-ialv(ilv)-iintf(ilv)
    IF(j.GT.nnsend_m) THEN
!     write(*,*)'nsend_m',ilv,j,nnsend_m
     nnsend_m = j
    endif
    
     IF(j.GT.nnrecv_m) THEN
!        write(*,*)'nnrecv_m',ilv,j,nnrecv_m
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
!ALLOCATE(coord1(ndim,nnode))
ALLOCATE(rintf1(nnrecv_m),sintf1(nnsend_m))

! iai1/iar1 은 레벨 루프의 재사용 버퍼 — 행 수가 레벨별 nnode1(=ialv 차분)이므로
! 최대 레벨 폭으로 할당해야 함. nnode(fine)로 잡으면 극소 도메인(대규모 np)에서
! coarse 레벨 폭 > nnode 가 되어 오버런 → 힙 오염 (np=900 비결정 크래시, LOG C010-3)
i = MAXVAL(ialv(2:nlevel+1)-ialv(1:nlevel))
IF(i.LT.nnode) i = nnode
ALLOCATE(iai1(i+1),jai1(nnzi),Xintp1(nnzi))
ALLOCATE(iar1(i+1),jar1(nnzr),Xrest1(nnzr))
!/
!  IF(minval([nnode, nnzc0]) == 0) THEN
!      write(myrank+1000,*)'error for nnode, zero point',myrank,nnode
!      write(*,*)'error for nnode, zero point',myrank,nnode
!  ENDIF
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
! ALLOCATE : 
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
! - - - - - - - - - - - - - - 
WRITE(fout,'(A,I0.3,A)') 'MG_tmp/part_MG', myrank+1, '.out'
OPEN(newunit=iu,file=fout,status='old',action='read',iostat=alstatus)
IF(alstatus/=0) THEN
WRITE(*,*)'read_mesh_MPI: cannot open ',TRIM(fout),' rank',myrank
STOP
ENDIF
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !
DO ilv = 2,nlevel
    
READ(iu,*)                                            ! reading
nintf1 = iintf(ilv)
nnode1 = ialv(ilv+1)-ialv(ilv)
nnbd1 = inbdc(ilv)
nnode1gl = inodegl(ilv)
nelem1 = nnode1
nnode0 = ialv(ilv)-ialv(ilv-1)

! 1:reading local

READ(iu,*) i1,i2,i3,i4
IF((i1.NE.nintf1).OR.(i2.NE.nnode1).OR.(i3.NE.nnbd1).OR.(i4.NE.nnode1gl)) THEN
WRITE(*,*)'error nintf1'
STOP
ENDIF
! NEW for A and R, P - - - - - - 
IF(ilv.NE.nlevel) THEN
READ(iu,*) inbdcA(ilv),inbdcR(ilv),inbdcP(ilv)
ELSE 
READ(iu,*) inbdcA(ilv),inbdcP(ilv)
ENDIF
! - - - - - - - - - - - - - - -   
DO i=1,nnode1gl
   READ(iu,*) coord1(1:ndim,i)                        ! reading
ENDDO

IF(nnbd1.NE.0) THEN
    
READ(iu,*) nbdom1(1:nnbd1)                            ! reading 
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(nnode1-nintf1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))

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

READ(iu,*) nbdom1(1:nnbd1)                            ! reading 
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))

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

READ(iu,*) nbdom1(1:nnbd1)                            ! reading 
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))

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

READ(iu,*) nbdom1(1:nnbd1)                            ! reading 
READ(iu,*) rpt1(1:nnbd1+1)
READ(iu,*) spt1(1:nnbd1+1)

READ(iu,*) rintf1(1:(rpt1(nnbd1+1)-1))                  ! reading
READ(iu,*) sintf1(1:(spt1(nnbd1+1)-1))

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
read(iu,*)                                        ! reading

iai1(1) = 1
Do i=1,nnode0
read(iu,*) nnd,id(1:nnd)
read(iu,*) tmp(1:nnd)
 k=iai1(i)
 do j=1,nnd
 jai1(k)=id(j)
 Xintp1(k) = tmp(j)
  k=k+1
 enddo
 iai1(i+1)=k
enddo

! R
READ(iu,*)                                        ! reading
iar1(1) = 1
DO i=1,nnode1
READ(iu,*) nnd,id(1:nnd)
READ(iu,*) tmp(1:nnd)
 k=iar1(i)
 DO j=1,nnd
 jar1(k)=id(j)
 Xrest1(k) = tmp(j)
  k=k+1
 ENDDO
 iar1(i+1)=k
ENDDO

! AC
READ(iu,*)                                       ! reading
ia1(1) = 1
DO i=1,nnode1
READ(iu,*) nnd,id(1:nnd)
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
  READ(iu,*)                                         ! reading
  READ(iu,*) nnodeC,nnodeG,nnzG
  
!/
!  IF(minval([nnodeC, nnodeG, nnzG]) == 0) THEN
!      write(myrank,*)'error for n_GC, zero point',myrank,nnodeC,nnodeG, nnzG
!  ENDIF
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
  DO i=1,nnodeC  
   READ(iu,*) imapG(i)
  ENDDO
  ENDIF
  

! AC_G

 iaG(1) = 1
 DO i=1,nnodeG
   READ(iu,*) nnd,id(1:nnd),coordG(1:ndim,i)
   k=iaG(i)
   DO j=1,nnd
    jaG(k)=id(j)
    IF(id(j).EQ.i) juG(i) = k
    k=k+1
   ENDDO
 
   iaG(i+1)=k
 
  ENDDO

   IF(iaG(nnodeG+1).GT.(nnzG+1))THEN
    WRITE(*,*)'error in read AcG'
    PAUSE
   ENDIF

  ENDIF
  
CLOSE(iu)
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



