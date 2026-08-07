      MODULE Ztplot
!      
      IMPLICIT NONE
      SAVE
!
      INTEGER tplot_num,tplot_prop
      INTEGER,ALLOCATABLE :: tplot_cell(:),tplot_cell_loc(:),tplot_cell_rank(:)
      REAL(8) tplot_dt,time2view
      REAL(8),ALLOCATABLE :: tplot_x(:),tplot_y(:),tplot_z(:) 
!
      END MODULE Ztplot
           
