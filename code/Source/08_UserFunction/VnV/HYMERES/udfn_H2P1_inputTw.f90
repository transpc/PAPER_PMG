!
      SUBROUTINE udfn_H2P1_inputTw
!
      USE Zconst1      , ONLY: vv_prob
      USE Zconst2      , ONLY: grav
      USE Zcoord1      , ONLY: xloc
      USE Znum_cell    , ONLY: istart_nf
      USE Ztimecon     , ONLY: time
      USE Zmodel       , ONLY: h2p1_tw
      USE Zvec_index   , ONLY: left_nf
!
      IMPLICIT NONE
!
      INTEGER i,j,ii
      INTEGER nf_number,istart,isize,i1
      INTEGER n_inputTw
      INTEGER,SAVE:: n_inputTw1
      REAL(8) tempTw(9),FACT
      REAL(8),SAVE:: locTw(9)
      REAL(8),SAVE,ALLOCATABLE:: inputTw(:,:)
!
      LOGICAL, SAVE::initialTw
      DATA initialTw/.TRUE./
!
      IF(initialTw)THEN
         initialTw=.FALSE.
         locTw=(/0.153d0, 0.312d0, 1.0d0, 1.8d0, 4.45d0, 5.775d0, &
                 7.1d0, 7.8d0, 7.959d0/)
!
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or.vv_prob.eq.'VD_h2p1_0')THEN       !2000s
            n_inputTw=1043
         ELSEIF(vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x')THEN   !2400s
            n_inputTw=1205
         ELSEIF(vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x')THEN   !2400s
            n_inputTw=1213
         ELSEIF(vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x')THEN   !4500s
            n_inputTw=2295
         ELSEIF(vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x')THEN   !5000s
            n_inputTw=2501
         ENDIF
         n_inputTw1=n_inputTw-1
         ALLOCATE(inputTw(10,n_inputTw))
         ALLOCATE(h2p1_tw(istart_nf(2,6)))
!
         IF(vv_prob.eq.'VD_h2p1_0')THEN
            OPEN(875,file='inputTw.dat',status='old',form='unformatted')
            DO i=1,n_inputTw
               READ(875) inputTw(1,i),inputTw(2,i),inputTw(3,i),inputTw(4,i),inputTw(5,i), &
                         inputTw(6,i),inputTw(7,i),inputTw(8,i),inputTw(9,i),inputTw(10,i)
            ENDDO
         ELSE
            OPEN(875,file='inputTw.dat',status='old')
            READ(875,*) inputTw(:,:)
         ENDIF
         CLOSE(875)
      ENDIF
!
      tempTw(:)=0.0d0
      IF(time.le.inputTw(1,1))THEN
         DO j=1,9
            tempTw(j)=inputTw(j+1,1)
         ENDDO
      ELSE
         DO i=1,n_inputTw1
            IF(time.gt.inputTw(1,i).and.time.le.inputTw(1,i+1))THEN
               FACT=(time-inputTw(1,i))/(inputTw(1,i+1)-inputTw(1,i))
               DO j=1,9
                  tempTw(j)=inputTw(j+1,i)+FACT*(inputTw(j+1,i+1)-inputTw(j+1,i))
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      h2p1_tw(:)=0.0d0
      nf_number=6
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         ii=left_nf(i1)
         IF(grav(2).gt.grav(3))THEN               !height:z
            IF(xloc(ii,3).le.locTw(1))THEN        !0.153
               h2p1_tw(i)=tempTw(1)
            ELSEIF(xloc(ii,3).gt.locTw(9))THEN    !7.959
               h2p1_tw(i)=tempTw(9)
            ELSE
               DO j=1,8
                  IF(xloc(ii,3).gt.locTw(j).and.xloc(ii,3).le.locTw(j+1))THEN
                     FACT=(xloc(ii,3)-locTw(j))/(locTw(j+1)-locTw(j))
                     h2p1_tw(i)=tempTw(j)+FACT*(tempTw(j+1)-tempTw(j))
                  ENDIF
               ENDDO
            ENDIF
         ELSEIF(grav(2).lt.grav(3))THEN           !height:y
            IF(xloc(ii,2).le.locTw(1))THEN         !0.153
               h2p1_tw(i)=tempTw(1)
            ELSEIF(xloc(ii,2).gt.locTw(9))THEN    !7.959
               h2p1_tw(i)=tempTw(9)
            ELSE
               DO j=1,8
                  IF(xloc(ii,2).gt.locTw(j).and.xloc(ii,2).le.locTw(j+1))THEN
                     FACT=(xloc(ii,2)-locTw(j))/(locTw(j+1)-locTw(j))
                     h2p1_tw(i)=tempTw(j)+FACT*(tempTw(j+1)-tempTw(j))
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
         h2p1_tw(i)=h2p1_tw(i)+273.16d0
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_H2P1_inputTw