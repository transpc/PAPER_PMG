MODULE Zio_unit
   IMPLICIT NONE
!...IO units
!......Reserved unit numbers for Fortran
   INTEGER, PARAMETER :: unit_keyborad=5,unit_screen=6
!......CUPID inherent unit numbers
   INTEGER :: unit_log, unit_somaflow, unit_grid, unit_rv,unit_chn
   INTEGER :: unit_saveout, unit_restart
   INTEGER :: unit_restart_rv, unit_gap_cond
   INTEGER :: unit_tplotv(200), unit_tplots(200)
!......User defined unit numbers
   
END MODULE Zio_unit