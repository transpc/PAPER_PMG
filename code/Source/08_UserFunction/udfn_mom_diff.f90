!
      SUBROUTINE udf_mom_diff(diff_drp)
!      
!     User-defined momentum diffusion for droplets (only when 'udfl_mom_diff' is used.)
!     
      USE Zparam   , ONLY: ndim
      USE Zmpi     , ONLY: ncell_fp
!
      IMPLICIT NONE
!
      REAL(8) diff_drp(ncell_fp,ndim)
!      
      diff_drp(:,:)=0.0d0
!
      RETURN
      END SUBROUTINE udf_mom_diff
