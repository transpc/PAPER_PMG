!
      SUBROUTINE udfn_turb_zeq
!
!     User-defined Zero equation. Only for s_turb_zero='user'
!
      USE VOL_DATA     , ONLY: cell
      USE Zconst1      , ONLY: vv_prob
      USE Zturb        , ONLY: yplus
      USE Zvector      , ONLY: ul_o
      USE Zzone        , ONLY: ncell_fluid
!
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) Rel,utau_cell,fric
!      
!.....PAFS-POOL, fluidic_device
!
      IF(vv_prob.eq.'PAFS-POOL'.or.vv_prob.eq.'fluidic_device'.or.vv_prob.eq.'UPTF-RV')THEN
         DO i=1,ncell_fluid
            IF(cell%regime(i).eq.11.or.cell%regime(i).eq.12.or.cell%regime(i).eq.21.or.cell%regime(i).eq.22.or.cell%regime(i).eq.3)THEN
               Rel=DMAX1(1.0d0,cell%rhol(i)*ul_o(i)*(2.0d0*0.11d0)/cell%lviscosl(i))
               CALL mom_wall_kakac(Rel,fric)  
               utau_cell=ul_o(i)*(fric/2.0d0)**0.5d0
               yplus(i)=0.020431d0*utau_cell*cell%rhol(i)/cell%lviscosl(i)
               cell%tviscosl(i)=0.4d0*cell%lviscosl(i)*yplus(i)*(1.0d0-DEXP(-0.0017d0*(yplus(i)**2.0d0)))   
               cell%tviscosl(i)=DMIN1(5000.0d0*cell%lviscosl(i),cell%tviscosl(i))     
            ELSE
               cell%tviscosl(i)=0.0d0
            ENDIF
         ENDDO
      ENDIF
!      
      RETURN
      END SUBROUTINE udfn_turb_zeq
      
      
