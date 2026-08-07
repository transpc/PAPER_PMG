!
      SUBROUTINE udfn_sij
!
!     Modifies the momentum source terms at the free surface cells
!
      USE Zconst1        , ONLY: vv_prob
      USE Zporous        , ONLY: sgap,s_ij_non_i,s_ij_non_k
      USE Znum_cell      , ONLY: istart_nf
      !OPR1000 rod-scale 
      USE Zcoord2        , ONLY: cell_leng
      USE Zvec_index     , ONLY: left_nf,right_non
      USE Zporous        , ONLY: chn_type
      USE Zvec_geo       , ONLY: xn_nf
!
      IMPLICIT NONE
!            
      INTEGER :: i
      INTEGER :: ii,kk
      INTEGER :: nf_number,istart,len,i1
      REAL(8) :: xn_x,xn_y,xn_z
      !OPR1000 rod-scale 
      REAL(8) rod0,pit0,pit1,gap0,gap1
!
      nf_number=0
      istart=istart_nf(1,nf_number)
      len   =istart_nf(2,nf_number)
          
      DO i=1,len
         s_ij_non_i(i)=0.d0
         s_ij_non_k(i)=0.d0
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         IF(abs(xn_nf(i1,1)).eq.1.0d0)THEN
!====> ii
            s_ij_non_i(i)=sgap(ii,1)
!====> kk
            s_ij_non_k(i)=sgap(kk,1)
         ELSEIF(abs(xn_nf(i1,2)).eq.1.0d0)THEN
!====> ii
            s_ij_non_i(i)=sgap(ii,2)
!====> kk
            s_ij_non_k(i)=sgap(kk,2)
         ENDIF
      ENDDO
!
!.....iSMR Subchannel scale should be deleted.      !PSH
!
      IF(vv_prob.eq.'KSMR'    )then
         s_ij_non_i(:)=0.d0
         s_ij_non_k(:)=0.d0

         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
            
         DO i=1,len
            s_ij_non_i(i)=0.d0
            s_ij_non_k(i)=0.d0
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(chn_type(ii).eq.3.and.chn_type(kk).eq.2) THEN
               s_ij_non_i(i)=0.008222d0
            ELSEIF(chn_type(ii).eq.2.and.chn_type(kk).eq.1) THEN
               s_ij_non_i(i)=0.003372d0
            ELSEIF(chn_type(ii).eq.2.and.chn_type(kk).eq.2) THEN
               s_ij_non_i(i)=0.008222d0
            ELSEIF(chn_type(ii).eq.1.and.chn_type(kk).eq.1) THEN
               s_ij_non_i(i)=0.003372d0
            ENDIF
            !
            IF(chn_type(kk).eq.3.and.chn_type(ii).eq.2) THEN
               s_ij_non_k(i)=0.008222d0
            ELSEIF(chn_type(kk).eq.2.and.chn_type(ii).eq.1) THEN
               s_ij_non_k(i)=0.003372d0
            ELSEIF(chn_type(kk).eq.2.and.chn_type(ii).eq.2) THEN
               s_ij_non_k(i)=0.008222d0
            ELSEIF(chn_type(kk).eq.1.and.chn_type(ii).eq.1) THEN
               s_ij_non_k(i)=0.003372d0
            ENDIF
         ENDDO      
      ENDIF
!
!.....rod-scale (OPR1000/APR1400) TO BE Deleted !!!
!
      IF(vv_prob.eq.'APR1400_fullcore_modmesh01'           .or. &
         vv_prob.eq.'APR1400_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                       )then
          
         s_ij_non_i(:)=0.d0
         s_ij_non_k(:)=0.d0
         ! rod0 (Rod radius)                  : 0.00475
         ! len0 (Subchannel length)           : 0.012852
         ! len1 (Subchannel length at FA wall): 0.00798
         ! len11 (Normal gap)   : 0.012852-0.00475*2=0.003352
         ! len22 (Gap in edge)  : 0.00798-0.00475   =0.00323
         rod0 =0.004750d0
         pit0 =0.012852d0
         pit1 =0.007980d0
         gap0=pit0-rod0*2.d0
         gap1=pit1-rod0
!====> non
         nf_number=0
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            
            xn_x=xn_nf(i1,1)
            xn_y=xn_nf(i1,2)
            xn_z=xn_nf(i1,3)
            IF(abs(xn_z.gt.0.5d0)) CYCLE
!
!====> ii
!
            IF(chn_type(ii).eq.0) CYCLE
!
!...........center channel            
            IF(chn_type(ii).eq.1)THEN
               s_ij_non_i(i)=gap0
!
!...........side channel
            ELSEIF(chn_type(ii).eq.2)THEN
               s_ij_non_i(i)=gap0
              !horizontal
               IF(cell_leng(ii,1).gt.cell_leng(ii,2))then
                  IF(kk.ne.0 .and. chn_type(kk).ne.0)then
                     IF(DABS(xn_x).gt.0.5d0)then
                        s_ij_non_i(i)=gap1
                     ELSEIF(DABS(xn_y).gt.0.5d0)then
                        IF(chn_type(kk).eq.2)THEN
                           s_ij_non_i(i)=pit0
                        ELSEIF(chn_type(kk).eq.1)THEN
                           s_ij_non_i(i)=gap0
                        ELSE
                           STOP 'check sij for horizontal sc_type=2'
                        ENDIF
                     ENDIF
                  ENDIF
              !vertical   
               ELSEIF(cell_leng(ii,1).lt.cell_leng(ii,2))then
                  IF(kk.ne.0 .and. chn_type(kk).ne.0)then
                     IF(DABS(xn_y).gt.0.5d0)then
                        s_ij_non_i(i)=gap1
                     ELSEIF(DABS(xn_x).gt.0.5d0)then
                        IF(chn_type(kk).eq.2)THEN
                           s_ij_non_i(i)=pit0
                        ELSEIF(chn_type(kk).eq.1)THEN
                           s_ij_non_i(i)=gap0
                        ELSE
                           STOP 'check sij for vertical sc_type=2'
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
                  STOP 'check sij for horizontal/vertical sc_type=2'
               ENDIF
!
!...........corner channel
            ELSEIF(chn_type(ii).eq.3)THEN
               s_ij_non_i(i)=pit1
               IF(kk.ne.0)then
                  IF(chn_type(kk).eq.2)THEN
                     s_ij_non_i(i)=gap1
                  ENDIF
               ENDIF
!
!...........guide tube channel (center)
            ELSEIF(chn_type(ii).eq.4)THEN
               s_ij_non_i(i)=gap0*0.02d0 !2% of gap0 (gap0*0.02=0.000067)
!
!...........guide tube channel (perpendicular)
            ELSEIF(chn_type(ii).eq.5)THEN
               s_ij_non_i(i)=0.003398d0
               IF(chn_type(kk).eq.4)THEN
                  s_ij_non_i(i)=gap0*0.02d0 !2% of gap0 (gap0*0.02=0.000067)
               ELSEIF(chn_type(kk).eq.1)THEN
                  s_ij_non_i(i)=gap0
               ENDIF
!
!...........guide tube channel (diagonal)
            ELSEIF(chn_type(ii).eq.6)THEN
               s_ij_non_i(i)=0.003398d0
               IF(chn_type(kk).eq.1)THEN
                  s_ij_non_i(i)=gap0
               ENDIF
            ENDIF
!
!====> kk
!
            IF(chn_type(kk).eq.0) CYCLE
!
!...........center channel            
            IF(chn_type(kk).eq.1)THEN
               s_ij_non_k(i)=gap0
!
!...........side channel
            ELSEIF(chn_type(kk).eq.2)THEN
               s_ij_non_k(i)=gap0
              !horizontal
               IF(cell_leng(kk,1).gt.cell_leng(kk,2))then
                  IF(ii.ne.0 .and. chn_type(ii).ne.0)then
                     IF(DABS(xn_x).gt.0.5d0)then
                        s_ij_non_k(i)=gap1
                     ELSEIF(DABS(xn_y).gt.0.5d0)then
                        IF(chn_type(ii).eq.2)THEN
                           s_ij_non_k(i)=pit0
                        ELSEIF(chn_type(ii).eq.1)THEN
                           s_ij_non_k(i)=gap0
                        ELSE
                           STOP 'check sij of kk for horizontal sc_type=2'
                        ENDIF
                     ENDIF
                  ENDIF
              !vertical   
               ELSEIF(cell_leng(kk,1).lt.cell_leng(kk,2))then
                  IF(ii.ne.0 .and. chn_type(ii).ne.0)then
                     IF(DABS(xn_y).gt.0.5d0)then
                        s_ij_non_k(i)=gap1
                     ELSEIF(DABS(xn_x).gt.0.5d0)then
                        IF(chn_type(ii).eq.2)THEN
                           s_ij_non_k(i)=pit0
                        ELSEIF(chn_type(ii).eq.1)THEN
                           s_ij_non_k(i)=gap0
                        ELSE
                           STOP 'check sij of kk for vertical sc_type=2'
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
                  STOP 'check sij of kk for horizontal/vertical sc_type=2'
               ENDIF
!
!...........corner channel
            ELSEIF(chn_type(kk).eq.3)THEN
               s_ij_non_k(i)=pit1
               IF(ii.ne.0)then
                  IF(chn_type(ii).eq.2)THEN
                     s_ij_non_k(i)=gap1
                  ENDIF
               ENDIF
!
!...........guide tube channel (center)
            ELSEIF(chn_type(kk).eq.4)THEN
               s_ij_non_k(i)=gap0*0.02d0 !2% of gap0 (gap0*0.02=0.000067)
!
!...........guide tube channel (perpendicular)
            ELSEIF(chn_type(kk).eq.5)THEN
               s_ij_non_k(i)=0.003398d0
               IF(chn_type(ii).eq.4)THEN
                  s_ij_non_k(i)=gap0*0.02d0 !2% of gap0 (gap0*0.02=0.000067)
               ELSEIF(chn_type(ii).eq.1)THEN
                  s_ij_non_k(i)=gap0
               ENDIF
!
!...........guide tube channel (diagonal)
            ELSEIF(chn_type(kk).eq.6)THEN
               s_ij_non_k(i)=0.003398d0
               IF(chn_type(ii).eq.1)THEN
                  s_ij_non_k(i)=gap0
               ENDIF
            ENDIF
!
         ENDDO   
      ENDIF
!
      RETURN
      END SUBROUTINE udfn_sij
