      MODULE Wall_DATA
!
      TYPE WAL_DATA;SEQUENCE
         REAL(8),ALLOCATABLE :: twall_partition(:)
         REAL(8),ALLOCATABLE :: wall_fluxl_diff(:),wall_fluxg_diff(:),wall_fluxd_diff(:)
         REAL(8),ALLOCATABLE :: ddepartw(:),ratio_evap(:),bfreq(:)
      ENDTYPE
!
      TYPE(WAL_DATA)::face
!
      ENDMODULE
