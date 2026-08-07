!
      SUBROUTINE int_htc_interp_temp(i,ag,hil_01,hig_01,hil_09,hig_09)
!
!     Linear interpolation of Hik for temperature when alphag < 0.1 or  alphag > 0.9
!
      USE VOL_DATA                
      USE Zdhda       , ONLY: dHldag,dHgdag,dHldtl,dHgdtg,dHfgdtg
      USE Zmodel      , ONLY: H_il_min,H_ig_min,dtl,dtg
      USE Ztimecon    , ONLY: alpha_min
      USE Zqvol       , ONLY: H_ig,H_il
!
      IMPLICIT NONE
!      
      INTEGER i
!
      REAL(8) ag
      REAL(8) Hila,Higa
      REAL(8) hil_01,hig_01,hil_09,hig_09
!
      dHldtl(i)=0.0d0
      dHgdtg(i)=0.0d0
!
!.....Subcooled liquid only
!
      IF(ag.le.2.0d0*alpha_min)THEN
         IF(cell%tl(i).lt.cell%ts(i)) H_il(i)=0.0d0
      ENDIF
!
!.....Superheated steam only
!
      IF(ag.ge.(1.0d0-2.0d0*alpha_min))THEN
         IF(cell%tg(i).gt.cell%ts(i)) H_ig(i)=0.0d0
      ENDIF
!
      IF(ag.le.0.1d0)THEN
!
!........Linear interpolation for liquid temperature when alphag < 0.1
!
         IF(cell%tl(i).lt.cell%ts(i))THEN
            Hila=0.0d0
         ELSEIF(cell%tl(i).lt.cell%ts(i)+dtl)THEN
            dHldtl(i)=H_il_min/dtl
            Hila=dHldtl(i)*(cell%tl(i)-cell%ts(i))
         ELSE
            Hila=H_il_min
         ENDIF
         dHldag(i)=10.d0*(hil_01-Hila)
         H_il(i)=Hila+dHldag(i)*ag
         Higa=H_ig_min
         dHgdag(i)=10.d0*(hig_01-Higa)
         H_ig(i)=Higa+dHgdag(i)*ag
!
      ELSEIF(ag.ge.0.9d0)THEN
!
!........Linear interpolation for gas temperature when alphag > 0.9
!
         IF(cell%tg(i).ge.cell%ts(i))THEN
            Higa=0.0d0
         ELSEIF(cell%tg(i).gt.cell%ts(i)-dtg)THEN
            dHgdtg(i)=H_ig_min/dtg
            Higa=dHgdtg(i)*(cell%ts(i)-cell%tg(i))
         ELSE
            Higa=H_ig_min
         ENDIF
!  
         dHgdag(i)=10.d0*(Higa-hig_09)
         H_ig(i)=Higa+dHgdag(i)*(ag-1.0d0)
         Hila=H_il_min
         dHldag(i)=10.d0*(Hila-hil_09)
         H_il(i)=Hila+dHldag(i)*(ag-1.0d0)
!
      ENDIF
      dHfgdtg(i)=dHgdtg(i)
!
      RETURN
      END SUBROUTINE int_htc_interp_temp
