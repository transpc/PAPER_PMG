!
      SUBROUTINE udfn_porous_hyd(x,hd)
!
!.....This routine change the cell value of somaGrid.
!
      USE Zparam      , ONLY: ndim
      USE Zconst1     , ONLY: vv_prob
!
      IMPLICIT NONE
!     
      REAL(8) x(ndim),hd
!     
!       
!.....rbht_1196_3d
!
      IF(vv_prob.eq.'rbht_1196_3d'.or.vv_prob.eq.'rbht_1196_3d2') THEN
!        corner
         IF((x(1).lt.0.02676.or.x(1).gt.0.06455).and.(x(2).lt.0.02676.or.x(2).gt.0.06455)) THEN
            hd=0.009993562
!        side (horizontal)
         ELSEIF((x(1).gt.0.02676.and.x(1).lt.0.06455).and.(x(2).lt.0.02676.or.x(2).gt.0.06455)) THEN
            hd=0.010810752
!        side (vertical)
         ELSEIF((x(1).lt.0.02676.or.x(1).gt.0.06455).and.(x(2).gt.0.02676.and.x(2).lt.0.06455)) THEN
            hd=0.010810752
!        center            
         ELSEIF((x(1).gt.0.02676.and.x(1).lt.0.06455).and.(x(2).gt.0.02676.and.x(2).lt.0.06455)) THEN
            hd=0.011773735
         ENDIF
      ENDIF  
!
      IF(vv_prob.eq.'rbht_1196fine2'.or.vv_prob.eq.'rbht_1196_3d_uniform') THEN
         hd=0.0106429d0
      ENDIF                       
!
      RETURN
      ENDSUBROUTINE udfn_porous_hyd
!
!======================================================================
!
      SUBROUTINE udfn_porous_hyd1(ncell,nzone,hd)
!
!.....This routine change the cell value of somaGrid.
!
      USE Znum_cell   , ONLY: i_neigh_tmp,j_neigh_tmp,j_nbcon_tmp
      USE Zconst1     , ONLY: vv_prob
      USE Zporous , ONLY: chn_type_tmp
!
      IMPLICIT NONE
!     
!     input
      INTEGER ncell
      INTEGER nzone(ncell)
      REAL(8) hd(ncell)
!
      INTEGER i,j,k

!
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1'.or. &
         vv_prob.eq.'opr1000_mc_rv'         .or. &
         vv_prob.eq.'opr1000_rv_lbloca'     .or. &
         vv_prob.eq.'opr1000_rv'                  )THEN
         hd(:)=0.012637d0
         DO i=1,ncell
            IF(nzone(i).eq.6 .or. nzone(i).eq.8)THEN
               hd(i)=0.012637d0
            ENDIF
         ENDDO
      ENDIF
      IF(vv_prob.eq.'apr1400_mc_rv' .or. &
         vv_prob.eq.'apr1400_rv')THEN
         hd(:)=0.012637d0
         DO i=1,ncell
            IF(nzone(i).eq.6 .or. nzone(i).eq.8)THEN
               hd(i)=0.012637d0
            ENDIF
         ENDDO
      ENDIF

      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
         vv_prob.eq.'OPR1000_single_assem'                       ) then
         hd(:)=0.012637d0
if(0)then
         DO i=1,ncell
            IF(chn_type_tmp(i).eq.1)THEN
               hd(i)=0.012637d0
            ELSEIF(chn_type_tmp(i).eq.2)THEN
               hd(i)=0.017984d0
              !search lateral neighbors
               k=0
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1 !-2 
                  IF(j_nbcon_tmp(j).eq.-1)then
                     k=k+1
                  ENDIF
               ENDDO
               IF(k.eq.1)then !wall cell
                  hd(i)=0.00967d0
               ELSEIF(k.ge.2)then
                  STOP 'check nbcon=-1 for sc_type=2'
               ENDIF
               !*** shroud neighboring
            ELSEIF(chn_type_tmp(i).eq.3)THEN
               hd(i)=0.024622d0
              !search lateral neighbors
               k=0
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1 !-2 
                  IF(j_nbcon_tmp(j).eq.-1)then
                     k=k+1
                  ENDIF
               ENDDO
               IF(k.eq.2)then
                  hd(i)=0.00785d0
               ELSEIF(k.eq.1)then
                  hd(i)=0.012d0
               ELSEIF(k.ge.3)then
                  STOP 'check nbcon=-1 for sc_type=3'
               ENDIF
            ELSEIF(chn_type_tmp(i).eq.4)THEN
               hd(i)=0.0001263747d0
            ELSEIF(chn_type_tmp(i).eq.5)THEN
               hd(i)=0.0076547307d0
            ELSEIF(chn_type_tmp(i).eq.6)THEN
               hd(i)=0.0136019357d0 !hydro_guide_corner_new
            ENDIF
         ENDDO
else
         DO i=1,ncell
            IF(chn_type_tmp(i).eq.1)THEN
               hd(i)=0.012637d0
            ELSEIF(chn_type_tmp(i).eq.2)THEN
               hd(i)=0.017984d0
              !search lateral neighbors
               DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1 
                  IF(chn_type_tmp(j_neigh_tmp(j)).eq.2)then
                     hd(i)=0.012637d0
                  ENDIF
               ENDDO
               !*** shroud neighboring
            ELSEIF(chn_type_tmp(i).eq.3)THEN
               hd(i)=0.024622d0
            ELSEIF(chn_type_tmp(i).eq.4)THEN
               hd(i)=0.0001263747d0
            ELSEIF(chn_type_tmp(i).eq.5)THEN
               hd(i)=0.0076547307d0
            ELSEIF(chn_type_tmp(i).eq.6)THEN
               hd(i)=0.0136019357d0 !hydro_guide_corner_new
            ENDIF
         ENDDO
endif

      ENDIF      

      RETURN
      ENDSUBROUTINE udfn_porous_hyd1

              
