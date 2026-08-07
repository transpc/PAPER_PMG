      SUBROUTINE read_tpfh2o
! 
!     This routine reads the steamtable file, tpfh2o
!
      USE STM_TBL_cupid , ONLY: st_tbl,                    &
                                nt,np,ns,ns2,ndxstd,nfluid
      USE Zio_unit      , ONLY: unit_log
!
      IMPLICIT NONE 
!
      INCLUDE 'stcom.h' 
      INCLUDE 'mxnfcd.h'
! 
      COMMON/std2xc/ntd2x,npd2x,nsd2x,ns2d2x,klpd2x,klp2d2,llpd2x,nt5d2x,jpld2x 
      INTEGER ntd2x,npd2x,nsd2x,ns2d2x,klpd2x,klp2d2,llpd2x,nt5d2x,jpld2x
      SAVE/std2xc/ 
!   
      INTEGER ist,isttbl,matmax,nUSE,nwleft
      INTEGER i
      INTEGER istmry,ist_tbl(2,1)
      INTEGER h2oin
!
      REAL stceqv(lstcom) 
!
      EQUIVALENCE(stceqv(1),ttrip) 
!
      CHARACTER(80) record(2) 
! 
      DATA istmry/2500000/                                                   !!! binary 물성치 table 크기에 따라 수정이 필요할 수 있음 (초기값16156)
      DATA tpfnam/'tpfh2o','tpfd2o','tpfco2','tpfhe','tpfh2','tpfo2',              &  !new gas
                  'tpfn2','tpfna','tpfar','tpfdair','tpflbe','tpfh2on','tpfflibe', &
                  'tpfflinabe','tpfh2onew','tpfr12','tpfr134a'/                       !17 fluids 
      DATA fsymbl/'h2o','d2o','co2','he','h2','o2','n2','na','ar',                 &  !new gas
                  'dair','lbe','h2on','flibe','flinabe','h2onew','r12','r134a'/
      DATA wmoles/18.01534d0,20.031d0,44.0098d0,4.0026d0,2.01594d0,   &
                  31.9988d0,28.01348d0,22.98977d0,39.948d0,28.9586d0, &
                  210.0d0,18.01534d0,99.d0,99.d0,18.01534d0,          &
                  72.6d0, 102.03d0/   
      DATA h2oin /0/
! 
!.....Test which thermodynamic property files are needed and bring them into memory. 
!.....Get table storage. 
!
      ALLOCATE(st_tbl(istmry))
      st_tbl(:)=0.d0
!.....This is default set in read_flow nfluid=1
!     If no default you must set it in somaflow.in under
!     &problem_description  namelist
!     nfluid=1
      ndxstd=17                   !!! 17개 이상 물성치 추가시 이 변수를 수정해야 하는지 확인이 반드시 필요함.
! 
!
!.....Find the maximum fluid number required for the problem by checking 
!     the material number of each volume. 
!
      matmax=1
! 
!.....Set base pointer for steam tables DATA block. 
!
      ist=0 
!
!.....Zero out first matmax words of steam tables DATA block. 
!
      st_tbl(2)=0.0d0 
! 
!.....Set pointer word to 1 for fluids that are to be USEd. 
!
      ist_tbl(2,1)=1 
! 
!.....Set pointer to first available word in steam tables DATA block for 
!     storage of steam tables DATA; also calculate number of words 
!     available for steam tables DATA. 
!
      isttbl=2
      nwleft=istmry-matmax 
! 
!.....Load DATA into steam tables DATA block for each required fluid. 
!
      OPEN(unit=16,file=tpfnam(nfluid),status='old',form='unformatted',iostat=h2oin)
!      
      IF(h2oin.ne.0)then
         WRITE(*,*)'Program was terminated:',tpfnam(nfluid),' is missing!'
         WRITE(unit_log,*)'Program was terminated:',tpfnam(nfluid),' is missing!'
         STOP
      ENDIF
! 
      nuse=nwleft-15
      IF(nfluid.eq.15)THEN
         CALL newstread_cupid(16,nuse,record,st_tbl(isttbl+15)) 
      ELSE
         CALL stread_cupid(16,nuse,record,st_tbl(isttbl+15)) 
      ENDIF
!
!.....Save /stcom/ DATA for this fluid. 
!
      DO i=1,lstcom
         st_tbl(isttbl+i-1)=stceqv(i) 
      END DO 
! 
!.....Set pointer relative to beginning of steam tables DATA block to DATA for this fluid.
!
      ist_tbl(2,1)=isttbl 
!
      IF(nfluid.ne.2)THEN
         nt=ntt
         np=npp
         ns=nst
         ns2=nsp
      ELSE
         ntd2x=ntt
         npd2x=npp
         nsd2x=nst 
         ns2d2x=nsp 
         klpd2x=it3bp 
         klp2d2=it4bp 
         llpd2x=it5bp 
         nt5d2x=nprpnt 
         jpld2x=it3p0 
      ENDIF
!
!.....Reformat steam table arrays
!
      CALL copy_tpfh2o(st_tbl(isttbl+15+ntt+npp),        &
                       st_tbl(isttbl+15+ntt+npp+13*nst))
!
!     Close steam tables DATA file for this fluid. 
!
      CLOSE(unit=16) 
! 
      END SUBROUTINE read_tpfh2o
!
      SUBROUTINE copy_tpfh2o(aa3,aa4)
!
      USE STM_TBL_cupid  , ONLY: a31,a3,a41,a4
!
      IMPLICIT NONE
!
      INCLUDE 'stcom.h'
!
!.....Input
      REAL(8) :: aa3(13,nst),aa4(13,nsp)
!.....Local variables
      INTEGER :: i,j
!
      ALLOCATE(a31(nst),a3(6,2,nst))
      ALLOCATE(a41(nsp),a4(6,2,nsp))
      DO i=1,nst
         a31(i)=aa3(1,i)
         DO j=1,6
            a3(j,1,i)=aa3(1+j,i)
            a3(j,2,i)=aa3(7+j,i)
         ENDDO
      ENDDO
      DO i=1,nsp
         a41(i)=aa4(1,i)
         DO j=1,6
            a4(j,1,i)=aa4(1+j,i)
            a4(j,2,i)=aa4(7+j,i)
         ENDDO
      ENDDO
!
      END SUBROUTINE copy_tpfh2o

