!
      FUNCTION IAT_rc_co_i(vreli,rholi,sigmali,psii,alphavi,ia_oldi,epsiloni,dsmi)
!
!     This routine calculates IAT source by random collision coalescence (2)
!
      USE Zparam,       ONLY:Kc1, kc2, kc3, weber_cr
!  
      IMPLICIT NONE
!      
      REAL(8) vreli,weber,root_weratio
      REAL(8) rholi,sigmali
      REAL(8) psii,alphavi,ia_oldi,epsiloni,dsmi
      REAL(8) IAT_rc_co_i
      REAL(8) coeff
      REAL(8) Gfactor_i
!
      weber=rholi*vreli**2.d0*dsmi/sigmali
      Root_weratio=(weber/weber_cr)**0.5d0
!
      iat_rc_co_i=0.d0
!
      IF(ia_oldi.gt.1.d-4)THEN
         coeff=(alphavi/ia_oldi)**2.d0
      ELSE
!
!........Small void or small area concentration --> bubbly
!
         coeff=dsmi/6.d0
         coeff=0.d0
      ENDIF
!
      IF(ia_oldi.gt.1.0d-4.and.epsiloni.gt.1.d-8)IAT_rc_co_i=-1.d0/3.d0/psii*coeff*Kc1  &
          *epsiloni**(1.d0/3.d0)*alphavi**2.d0/dsmi**(11.d0/3.d0)/(Gfactor_i(alphavi)+Kc2*alphavi*Root_weratio)  &
          *DEXP(-Kc3*Root_weratio) 
!
      RETURN
      ENDFUNCTION IAT_rc_co_i
!      
!
      FUNCTION Gfactor_i(alphavi)
!
!     This function calculates gfactor for subroutine IAT_rc_co_i
!      
      USE Zparam,       ONLY:alpha_max
!
      IMPLICIT NONE
!      
      REAL(8) alphavi
      REAL(8) gfactor_i
!
      Gfactor_i=1.d0-(alphavi/alpha_max)**(1.d0/3.d0)
!
      RETURN
      ENDFUNCTION Gfactor_i
!
!
      FUNCTION IAT_TI_BK_i(vreli,rholi,sigmali,psii,alphavi,ia_oldi,epsiloni,dsmi)

!
!     This routine calculates IAT source by turbulence induced breakup (3)
!
      USE Zparam,       ONLY:weber_cr, kb1, kb2
!
      IMPLICIT NONE    
!        
      REAL(8) vreli,weber,weratio,root_weratio
      REAL(8) rholi,sigmali
      REAL(8) psii,alphavi,ia_oldi,epsiloni,dsmi
      REAL(8) IAT_ti_bk_i
      REAL(8) coeff
!
      weber=rholi*vreli**2.*dsmi/sigmali
      Weratio=Weber/weber_cr
      Root_weratio=weratio**0.5
!
      IAT_ti_bk_i=0.d0
      IF(ia_oldi.gt.1.d-4)THEN
         coeff=(alphavi/ia_oldi)**2.d0
      ELSE
         coeff=dsmi/6.d0
         coeff=0.d0
      ENDIF
      IF(ia_oldi.gt.1.0d-4.and.epsiloni.gt.1.d-8)IAT_ti_bk_i=1.d0/3.d0/psii*coeff*Kb1 &
         *epsiloni**(1.d0/3.d0)*alphavi*(1.d0-alphavi)/dsmi**(11.d0/3.d0) &
         *1.0d0/(1.d0+Kb2*(1.d0-alphavi)*Root_weratio)*DEXP(-Weratio)
!
      RETURN
      ENDFUNCTION IAT_TI_BK_i
!
!
      FUNCTION IAT_nc_ph_i(d_departurei,n_nucleatei,f_departurei,A_heati,Voli)
!
!     This routine calculates IAT source by nucleation induced phase change (4)
!     Currently, not used because heat_partition model give this term.
!
      USE Zparam,       ONLY:pi
!
      IMPLICIT NONE
!      
      REAL(8) d_departurei,n_nucleatei,f_departurei
      REAL(8) A_heati,Voli
      REAL(8) IAT_nc_ph_i
!
      IAT_nc_ph_i=pi*d_departurei*n_nucleatei*f_departurei*A_heati*Voli
!
      RETURN
      ENDFUNCTION IAT_nc_ph_i
!
!
      FUNCTION D_departure_i(sigmai,gravity,rholi,rhogi,cpli,twalli,tsati,hlgi)
!
!     This routine calculates bubble departure diameter for subroutine IAT_nc_ph_i.
!
      IMPLICIT NONE
!            
      REAL(8) d_departure_i
      REAL(8) sigmai,gravity,rholi,rhogi,cpli,tsati,hlgi,twalli
!
      d_departure_i=0.d0
      IF(twalli.gt.tsati)d_departure_i=1.5d0*1.e-4*DSQRT(sigmai/gravity/(rholi-rhogi)) &
         *(rholi*cpli*(twalli-tsati)/rhogi/hlgi)**(5.d0/4.d0)
!
      RETURN
      ENDFUNCTION D_departure_i
!
!
      FUNCTION N_nucleate_i(twalli,tsati)
!
!     This routine calculates nucleate site density for subroutine IAT_nc_ph_i.
!
      IMPLICIT NONE   
!         
      REAL(8) n_nucleate_i
      REAL(8) twalli,tsati

      N_nucleate_i=0.d0
      IF(twalli.gt.tsati)N_nucleate_i=(185.d0*(twalli-tsati))**1.805d0
!
      RETURN
      ENDFUNCTION N_nucleate_i
!
!
      FUNCTION F_departure_i(gravity,rholi,rhogi,d_departurei,twalli,tsati)
!
!     This routine calculates departure frequency for subroutine IAT_nc_ph_i.
!
      IMPLICIT NONE
!      
      REAL(8) f_departure_i
      REAL(8) gravity,rholi,rhogi,d_departurei
      REAL(8) twalli,tsati
!
      F_departure_i=0.d0
      IF(twalli.gt.tsati)F_departure_i=DSQRT(4.d0*gravity*(rholi-rhogi)/3.d0*d_departurei*rholi)
!
      RETURN
      ENDFUNCTION F_departure_i
!
!
      FUNCTION A_heat_i(n_nucleatei,d_departurei)
!
!     This routine calculates two-phase heat transfer area for subroutine IAT_nc_ph_i.
!
      USE Zparam,       ONLY:pi,kfactor_iat
!
      IMPLICIT NONE
!      
      REAL(8) n_nucleatei,d_departurei
      REAL(8) A_heat_i
!
      A_heat_i=n_nucleatei*pi*d_departurei/4.d0*Kfactor_iat
!
      RETURN
      ENDFUNCTION A_heat_i
