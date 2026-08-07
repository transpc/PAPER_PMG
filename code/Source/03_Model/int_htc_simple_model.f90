!
      SUBROUTINE int_htc_simple_model(ag,al)
!
!     This routine calculates the heat transfer coefficient at the interface when iHTC=0
!
      USE VOL_DATA              
      USE Zmpi      , ONLY: ncell_fp
      USE Zmodel    , ONLY: h_il_coeff_a,h_il_coeff_b,  &
                             h_ig_coeff_a,h_ig_coeff_b, &
                             h_fg_coeff_a,h_fg_coeff_b
      USE Zqvol     , ONLY: h_ig,h_il,h_gf
      USE Ztimecon  , ONLY: alpha_min
      USE Zzone     , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) ag(ncell_fp),al(ncell_fp),agl
!

      DO i=1,ncell_fluid
         agl=ag(i)*al(i)
         H_il(i)=H_il_coeff_a+H_il_coeff_b*DMAX1(1.0e-5,agl)
         H_ig(i)=H_ig_coeff_a+H_ig_coeff_b*DMAX1(1.0e-5,agl)
         H_gf(i)=H_fg_coeff_a+H_fg_coeff_b*DMAX1(1.0e-5,agl)         
! 
         IF(ag(i).le.2.0d0*alpha_min)THEN
!
!...........Make HTC(l->i) zero for the subcooled liquid if void fraction is very low
!         
            IF(cell%tl(i).lt.cell%ts(i))THEN
               H_il(i)=0.0d0
            ENDIF
         ELSE
!
!...........Define the minimum HTC(l->i) for the superheated liquid
!          
            IF(cell%tl(i).gt.cell%ts(i))THEN
               H_il(i)=H_il_coeff_a+H_il_coeff_b*DMAX1(1.0e-5,agl)
               !IF(cell%tl(i)-cell%ts(i).lt.1.d0) THEN
               !   hr=1/10.d0+(1.d0-1/10.d0)/1.d0*(cell%tl(i)-cell%ts(i))
               !   H_il(i)=hr*(H_il_coeff_a+H_il_coeff_b*DMAX1(1.0e-5,agl))
               !ENDIF   
!
!...........suppress the minimum HTC(l->i) for the subcooled liquid
!             
            ELSE
               H_il(i)=1.0d0/10.d0*(H_il_coeff_a+H_il_coeff_b*DMAX1(1.0e-5,agl))
            ENDIF
         ENDIF
!
         IF(ag(i).ge.(1.0d0-2.0d0*alpha_min))THEN
!
!...........Make HTC(g->i) zero for the superheated liquid if void fraction is very high
!   
            IF(cell%tg(i).gt.cell%ts(i))THEN
               H_ig(i)=0.0d0
            ENDIF
         ELSE
         
            IF(cell%tg(i).gt.cell%ts(i))THEN
!
!...........suppress the minimum HTC(g->i) for the superheated vapor
! 
               H_ig(i)=1.0d0/10.d0*(H_ig_coeff_a+H_ig_coeff_b*DMAX1(1.0e-5,agl))
!
!...........Define the minimum HTC(g->i) for the subcooled vapor
!               
            ELSE
               H_ig(i)=H_ig_coeff_a+H_ig_coeff_b*DMAX1(1.0e-5,agl)
               !IF(cell%ts(i)-cell%tg(i).le.1.d0) THEN
               !   hr=1/10.d0+(1.d0-1/10.d0)/1.d0*(cell%ts(i)-cell%tg(i))
               !   H_ig(i)=hr*(H_ig_coeff_a+H_ig_coeff_b*DMAX1(1.0e-5,agl))
               !ENDIF 
            ENDIF
         ENDIF
!
!........Relaxation of Hik
!  
!!!        H_il(i)=(1.0d0-relax_hik)*H_il(i)+relax_hik*hil_o(i)
!!!        H_ig(i)=(1.0d0-relax_hik)*H_ig(i)+relax_hik*hig_o(i)
!!!        hil_o(i)=H_il(i)
!!!      hig_o(i)=H_ig(i)
!         H_gf(i)=H_fg_coeff_a+H_fg_coeff_b*DMAX1(1.0e-5,agl) 
      ENDDO
!
      RETURN
      END SUBROUTINE int_htc_simple_model
