!
      SUBROUTINE turb_ke_vis_lg
!
!     This routine calculates liquid and gas turbulent viscosity 
!     for dispersed phase by assumption as follows:
!     Bubble flow - alphag < 0.3
!     Mist flow - alphag > 0.8
!     Bubbly transition - 0.3 < alphag < 0.4
!     Mist transition   - 0.7 < alphag < 0.8
!
      USE VOL_DATA              
      USE Zbc_index , ONLY: icell_type
      USE Zturb     , ONLY: wvis_liq,wvis_gas
      USE Zzone     , ONLY: ncell_fluid
!
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) vis0,dab,dam,sl,rhogl,rholg
      REAL(8) ab1,ab2,am1,am2
!      
      DATA ab1,ab2,am1,am2/0.3d0,0.4d0,0.8d0,0.9d0/
!
      dab=ab2-ab1
      dam=am2-am1
!
      DO i=1,ncell_fluid
!
!........Bubbly flow
!
         IF(cell%alphag(i).le.ab1) THEN
            rhogl=cell%rhog(i)/cell%rhol(i)
!            
            IF(icell_type(i).eq.1) THEN
               wvis_gas(i)=wvis_liq(i)*rhogl
               wvis_gas(i)=DMAX1(cell%lviscosg(i),wvis_gas(i))
            ENDIF         
!            
          cell%tviscosg(i)=cell%tviscosl(i)*rhogl
          cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
!
!........Bubbly transition flow
!
         ELSEIF(cell%alphag(i).le.ab2) THEN
            rhogl=cell%rhog(i)/cell%rhol(i)
            IF(icell_type(i).eq.1) THEN
               vis0=wvis_liq(i)*rhogl
               sl=(wvis_gas(i)-vis0)/dab
               wvis_gas(i)=sl*(cell%alphag(i)-ab1)+vis0
               wvis_gas(i)=DMAX1(cell%lviscosg(i),wvis_gas(i))
            ENDIF         
!            
            vis0=cell%tviscosl(i)*rhogl
            sl=(cell%tviscosg(i)-vis0)/dab
            cell%tviscosg(i)=sl*(cell%alphag(i)-ab1)+vis0
            vis0=cell%eviscosl(i)*rhogl
            sl=(cell%eviscosg(i)-vis0)/dab
            cell%eviscosg(i)=cell%lviscosg(i)+cell%tviscosg(i)
          ENDIF
!
!........Mist flow
!
         IF(cell%alphag(i).ge.am2) THEN
!             
            rholg=cell%rhol(i)/cell%rhog(i)
!            
            IF(icell_type(i).eq.1) THEN
               wvis_liq(i)=wvis_gas(i)*rholg
               wvis_liq(i)=DMAX1(cell%lviscosl(i),wvis_liq(i))
            ENDIF
!            
            cell%tviscosl(i)=cell%tviscosg(i)*rholg
            cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
!
!........Mist transition flow
!
         ELSEIF(cell%alphag(i).ge.am1) THEN
!             
            rholg=cell%rhol(i)/cell%rhog(i)
!            
            IF(icell_type(i).eq.1) THEN
               vis0=wvis_gas(i)*rholg
               sl=(vis0-wvis_liq(i))/dam
               wvis_liq(i)=sl*(cell%alphag(i)-am2)+vis0
               wvis_liq(i)=DMAX1(cell%lviscosl(i),wvis_liq(i))
            ENDIF 
!            
            vis0=cell%tviscosg(i)*rholg
            sl=(vis0-cell%tviscosg(i))/dam
            cell%tviscosl(i)=sl*(cell%alphag(i)-am2)+vis0
            vis0=cell%eviscosg(i)*rholg
            sl=(vis0-cell%eviscosg(i))/dam
            cell%eviscosl(i)=cell%lviscosl(i)+cell%tviscosl(i)
!            
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_vis_lg
