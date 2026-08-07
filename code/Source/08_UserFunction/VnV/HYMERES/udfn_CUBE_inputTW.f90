!
      SUBROUTINE udfn_CUBE_inputTW         !need to refer udfn_H2P1_inputTw.
!
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord1      , ONLY: xloc
      USE Zmodel       , ONLY: cube_tw
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf
      USE Ztimecon     , ONLY: time
      USE Zvec_index   , ONLY: left_nf,nbcon_nf
!
      IMPLICIT NONE
!
      INTEGER i,j,ii,k
      INTEGER nf_number,istart,isize,i1,istart2,i2
      INTEGER n_inputTw
      INTEGER,SAVE:: n_inputTw1
      REAL(8) tempTw(7),FACT,tempTwf(2)
      REAL(8),SAVE:: locTw(7)
      REAL(8),SAVE,ALLOCATABLE:: inputTw(:,:)
!
      LOGICAL, SAVE::initialTw
      DATA initialTw/.TRUE./
!
      IF(initialTw)THEN
         initialTw=.FALSE.
!         locTw=(/4.98d0, 5.55d0, 7.15d0, 8.9592d0, 10.9592d0, 12.95936d0, 13.9592d0/)
         locTw=(/2.9592d0, 4.98d0, 6.9592d0, 8.9592d0, 10.9592d0, 12.95936d0, 13.9592d0/)
!
         IF(vv_prob.eq.'ST2-CT-01')THEN
            n_inputTw=8154
         ELSEIF(vv_prob.eq.'ST2-CT-02')THEN
            n_inputTw=8270
         ELSEIF(vv_prob.eq.'ST2-CT-03')THEN
            n_inputTw=11898
         ENDIF
         n_inputTw1=n_inputTw-1
         ALLOCATE(inputTw(10,n_inputTw))
         ALLOCATE(cube_tw(istart_nf(2,6)))
!
         OPEN(875,file='CUBE_TW.dat',status='old')
         READ(875,*) inputTw(:,:)
         CLOSE(875)
      ENDIF
!
      tempTw(:)=0.0d0
      tempTwf(:)=0.0d0
      IF(time.le.inputTw(1,1))THEN
         DO j=1,7
            tempTw(j)=inputTw(j+1,1)
         ENDDO
         tempTwf(1)=inputTw(9,1)
         tempTwf(2)=inputTw(10,1)
      ELSE
         DO i=1,n_inputTw1
            IF(time.gt.inputTw(1,i).and.time.le.inputTw(1,i+1))THEN
               FACT=(time-inputTw(1,i))/(inputTw(1,i+1)-inputTw(1,i))
               DO j=1,7
                  tempTw(j)=inputTw(j+1,i)+FACT*(inputTw(j+1,i+1)-inputTw(j+1,i))
               ENDDO
               tempTwf(1)=inputTw(9,i)+FACT*(inputTw(9,i+1)-inputTw(9,i))
               tempTwf(2)=inputTw(10,i)+FACT*(inputTw(10,i+1)-inputTw(10,i))
            ENDIF
         ENDDO
      ENDIF
!
      cube_tw(:)=0.0d0
      nf_number=6
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=-nbcon_nf(i2)
         IF(k.eq.3)THEN
            IF(xloc(ii,3).le.locTw(1))THEN        !4.98
               cube_tw(i)=tempTw(1)
            ELSEIF(xloc(ii,3).gt.locTw(7))THEN    !13.9592
               cube_tw(i)=tempTw(7)
            ELSE
               DO j=1,6
                  IF(xloc(ii,3).gt.locTw(j).and.xloc(ii,3).le.locTw(j+1))THEN
                     FACT=(xloc(ii,3)-locTw(j))/(locTw(j+1)-locTw(j))
                     cube_tw(i)=tempTw(j)+FACT*(tempTw(j+1)-tempTw(j))
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF(k.eq.4)THEN
            IF(xloc(ii,3).gt.locTw(2))THEN
               cube_tw(i)=tempTwf(1)         ! operating floor
            ELSE
               cube_tw(i)=tempTwf(2)         ! pedestal floor
            ENDIF
         ENDIF
         cube_tw(i)=cube_tw(i)+273.15d0
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_CUBE_inputTW
