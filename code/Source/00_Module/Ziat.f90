      MODULE Ziat
! 
      IMPLICIT NONE
      SAVE
!
      REAL(8) dbubble_init
      REAL(8),Allocatable::ia_conv(:),ia_old(:),ia(:),ia_b(:),dsm_b(:)
      REAL(8),Allocatable::iat_size(:),iat_coal(:),iat_break(:),iat_nucl(:)
!
      REAL(8) r_db_min
      REAL(8) r_db_max      
      REAL(8) r_ddrop
      CHARACTER(30) s_bubble_diameter
      CHARACTER(30) s_drop_diameter    !drop diamter
      REAL(8) r_dh_hibiki !s_bubble_diameter='hibiki'  
!      
      END MODULE Ziat