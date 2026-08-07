!
      SUBROUTINE udfn_tw_profile(twall_ctw)
!
!     temperature profile for Twall
!
      USE Zvec_param    , ONLY: nf_ctw
      USE Zvec_index    , ONLY: left_nf,nbcon_nf
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf
      USE Zb_condition  , ONLY: twall
      USE Zconst1       , ONLY: vv_prob
      USE Zconst2       , ONLY: grav
      USE Zcoord1       , ONLY: xloc
      USE Zmodel        , ONLY: h2p1_tw,h2p1_tp1,h2p1_tp2,cube_tw
      USE Zuserdefined  , ONLY: udfl_tw_profile
!
      IMPLICIT NONE
!
!.....Output
      REAL(8),DIMENSION(nf_ctw) :: twall_ctw
!.....Local variables
      INTEGER i,j,k,ii
      INTEGER :: nf_number,len,istart,istart2,i1,i2
!.....Local arrays
      REAL(8),DIMENSION(13),SAVE :: zkuhn
      REAL(8),DIMENSION(12),SAVE :: zTwall
      DATA zkuhn/2.418d0, 2.315d0, 2.181d0, 2.047d0, 1.897d0, 1.709d0,        &
                 1.531d0, 1.313d0, 1.097d0, 0.837d0, 0.569d0, 0.237d0, 0.0d0/
      DATA zTwall/366.15d0, 366.15d0, 365.45d0, 364.25d0, 364.45d0, 363.75d0, &
                  362.55d0, 362.05d0, 361.85d0, 361.85d0, 361.85d0, 361.85d0/
!
      nf_number=6
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len   =istart_nf(2,nf_number)
!
      IF(udfl_tw_profile) THEN
         IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
            vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
            vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
            vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
            vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
            vv_prob.eq.'VD_h2p1_0')THEN
            DO i=1,len
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=-nbcon_nf(i2)
               IF(k.eq.3)THEN
                  twall_ctw(i)=h2p1_tw(i)
               ELSEIF(k.eq.4)THEN
                   IF(grav(2).gt.grav(3))THEN             !height:z
                     IF(xloc(ii,3).lt.5.5d0)THEN
                        twall_ctw(i)=h2p1_tp1
                     ELSE
                        twall_ctw(i)=h2p1_tp2
                     ENDIF
                  ELSEIF(grav(2).lt.grav(3))THEN         !height:y
                     IF(xloc(ii,2).lt.5.5d0)THEN
                        twall_ctw(i)=h2p1_tp1
                     ELSE
                        twall_ctw(i)=h2p1_tp2
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
!         
         ELSEIF(vv_prob.eq.'kuhn_111')THEN
!
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               DO j=1,12
                  IF(xloc(ii,3).le.zkuhn(j).and.xloc(ii,3).gt.zkuhn(j+1)) THEN
                     twall_ctw(i)=zTwall(j)
                     EXIT
                  ENDIF
               ENDDO         
            ENDDO  
!
         ELSEIF(vv_prob.eq.'ST2-CT-01' .or. &
                vv_prob.eq.'ST2-CT-02' .or. &
                vv_prob.eq.'ST2-CT-03'       )THEN
!            IF(k.eq.3)THEN
!               twall(k)=cube_tw(i)
!            ENDIF
            DO i=1,len
               twall_ctw(i)=cube_tw(i)   ! allocated? 
            ENDDO
            
         ENDIF
      ELSE
         DO i=1,len
            i2=istart2+i
            k=-nbcon_nf(i2)
            twall_ctw(i)=twall(k)
         ENDDO
      ENDIF
!
      END SUBROUTINE udfn_tw_profile
