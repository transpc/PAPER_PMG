      MODULE Zrv_subchan
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,Allocatable:: subchannel_type_tmp(:),subchannel_type(:)
      INTEGER,Allocatable:: mv_fac(:),mv_non(:)
      REAL(8),Allocatable:: q_CHF_Biasi(:),APR_DNBR(:),q_local(:),qvol_liq_DNBR(:,:),qvol_DNBR(:),qvol_DNBR_dummy(:),rod_power(:,:)

      INTEGER n_sg,n_mv
      REAL(8) h_sg,h_mv
      REAL(8),Allocatable:: sg_loc(:),mv_loc(:)
      
      INTEGER nz_nk,nz_th0
      REAL(8),ALLOCATABLE:: dz_nk(:),hz_nk(:)
      REAL(8) qsum00
!
      END MODULE Zrv_subchan