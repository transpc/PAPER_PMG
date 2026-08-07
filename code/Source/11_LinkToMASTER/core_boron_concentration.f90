!-----------------------------------------------------------------------------------
      SUBROUTINE averaging_cboron
!
      USE VOL_DATA  ,ONLY: cell
      USE Zzone     ,ONLY: ncell_fluid,nzone
      USE Zcore     ,ONLY: np
      USE Zcoord3   ,ONLY: volp
      USE master4   ,ONLY: ppm_mas
!      
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i
      REAL(8) :: cboron_sum,cboron_avg,volp_sum 
!.....Local arrays
      REAL(8) :: tmp(2)
!      
      cboron_sum=0.d0
      volp_sum=0.d0
      DO i=1,ncell_fluid
         IF(nzone(i).eq.6)THEN
            cboron_sum=cboron_sum+cell%cboron(i)*volp(i)
            volp_sum=volp_sum+volp(i)
         ENDIF 
      ENDDO
      IF(np.gt.1)THEN
        tmp(1)=cboron_sum
        tmp(2)=volp_sum
        CALL allreducei_r(tmp,2)  
        cboron_sum=tmp(1)
        volp_sum  =tmp(2)
      ENDIF  
      cboron_avg=cboron_sum/volp_sum
      ppm_mas=cboron_avg*1.d6
!     
      END SUBROUTINE averaging_cboron
!
!-----------------------------------------------------------------------------------
!
      SUBROUTINE averaging_cboron_rod
!
      USE VOL_DATA   ,ONLY: cell
      USE Zzone      ,ONLY: ncell_fluid
      USE Zcore      ,ONLY: np
      USE Zcoord3    ,ONLY: volp
      USE master4    ,ONLY: ppm_mas
      USE Zporous    ,ONLY: chn_type
!      
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER i
      REAL(8) cboron_sum,cboron_avg,volp_sum 
!.....Local arrays
      REAL(8) :: tmp(2)
!      
      cboron_sum=0.d0
      volp_sum=0.d0
      DO i=1,ncell_fluid
         IF(chn_type(i).ne.0)THEN
            cboron_sum=cboron_sum+cell%cboron(i)*volp(i)
            volp_sum=volp_sum+volp(i)
         ENDIF 
      ENDDO
      IF(np.gt.1)THEN
        tmp(1)=cboron_sum
        tmp(2)=volp_sum
        CALL allreducei_r(tmp,2)  
        cboron_sum=tmp(1)
        volp_sum  =tmp(2)
      ENDIF  
      cboron_avg=cboron_sum/volp_sum
      ppm_mas=cboron_avg*1.d6
!     
      END SUBROUTINE averaging_cboron_rod
