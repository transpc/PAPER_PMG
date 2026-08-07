!
      SUBROUTINE udfn_flashing_hif(i,H_il_i,H_ig_i,ag)
!
!     User-defined H_il,H_ig for flashing
!        
      USE VOL_DATA             
      USE Zmodel   , ONLY: h_il_min
      USE Zqvol    , ONLY: h_il
!
      IMPLICIT NONE
!
      INTEGER i
      REAL(8) H_il_i,H_ig_i,coeff,ag
!
      IF(ag.gt.0.999d0)THEN
         IF(cell%tl(i).gt.cell%ts(i)) H_il_i=DMAX1(H_il_i,1.0e6)
      ENDIF
      IF(ag.lt.0.001d0)THEN
         IF(cell%tg(i).lt.cell%ts(i)) H_ig_i=DMAX1(H_ig_i,1.0e6)
      ENDIF
      coeff=0.05d0
      IF(cell%regime(i).eq.11.or.cell%regime(i).eq.12)THEN
         H_il_i=DMAX1(H_il_i,cell%alphal(i)*cell%rhol(i)*cell%cpl(i)/coeff) 
         H_il_i=DMIN1(1.0e8,DMAX1(H_il_i,H_il_min))
         H_il_i=0.99d0*H_il(i)+0.01d0*H_il_i
      ENDIF  
!
      RETURN
      END SUBROUTINE udfn_flashing_hif
