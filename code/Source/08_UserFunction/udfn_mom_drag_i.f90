!
      SUBROUTINE udfn_mom_drag_i(vfgl_i,n)
!
!     User-defined momentum drag (only when 'udfl_mom_drag_i' is used.)
! 
      USE VOL_DATA     , ONLY: cell
      USE Zconst1      , ONLY: vv_prob
!
      IMPLICIT NONE
!
      INTEGER n
      INTEGER i
!      
      REAL(8) vfgl_i(*)
!
      IF(vv_prob.eq.'ST2-CT-01'.or.vv_prob.eq.'ST2-CT-02'.or. &
         vv_prob.eq.'ST2-CT-03')THEN
         DO i=1,n
            IF(cell%alphag(i).gt.0.99999d0) vFgl_i(i)=dmax1(vFgl_i(i),1000.0d0)
         ENDDO
      ENDIF
!
!.....PAFS-POOL, fluidic_device, SMALL-POOL-3D
!      
      IF(vv_prob.eq.'PAFS-POOL'.or.vv_prob.eq.'fluidic_device')THEN 
         DO i=1,n
            IF(cell%alphag(i).gt.0.99999d0) vFgl_i(i)=dmax1(vFgl_i(i),1000.0d0)
         ENDDO 
!
!.....UPTF
!
      ELSEIF(vv_prob.eq.'UPTF-RV')THEN 
       DO i=1,n
         IF(cell%regime(i).ne.11) vFgl_i(i)=dmin1(100.d0,dmax1(vFgl_i(i),10.0d0))
!         IF(cell%regime(i).ne.11) vFgl_i(i)=dmin1(60.d0,dmax1(vFgl_i(i),10.0d0))
       ENDDO         
!          
!.....dam_break
!
      ELSEIF(vv_prob.eq.'dam_break')THEN 
       DO i=1,n
         IF(cell%regime(i).eq.3)THEN
            vFgl_i(i)=100000.0d0
         ENDIF
       ENDDO
      ENDIF
!
      RETURN
      END SUBROUTINE udfn_mom_drag_i
