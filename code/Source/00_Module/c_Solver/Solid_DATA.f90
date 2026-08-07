      MODULE Solid_DATA
!
      TYPE SOL_DATA;SEQUENCE
         REAL(8),ALLOCATABLE :: tsol(:)
         REAL(8),ALLOCATABLE :: rhos(:)
         REAL(8),ALLOCATABLE :: conds(:),cps(:),rhocps(:)
         REAL(8),ALLOCATABLE :: tsol_o(:)
         REAL(8),ALLOCATABLE :: tsol_max(:),tpellet_surf(:),temp_rod(:,:),hconv_rod_g(:),hconv_rod_l(:)
         INTEGER(4),ALLOCATABLE :: matnum(:),idummy(:)
      ENDTYPE
!
      TYPE(SOL_DATA) :: solid
!
      ENDMODULE
