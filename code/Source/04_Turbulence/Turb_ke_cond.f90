!
      SUBROUTINE turb_ke_cond
!
!     This routine calculates wall conductivities and effective 
!     conductivities.
!      
      USE VOL_DATA              
      USE Zparam    , ONLY: prt
      USE Zbc_index , ONLY: icell_type
      USE Zturb     , ONLY: wvis_liq,wvis_gas,wcd_liq,wcd_gas
      USE Zzone     , ONLY: ncell_fluid
!
      INTEGER i
!
!.....Turbulent thermal conductivity at wall
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1) THEN
            wcd_liq(i)=cell%lcondl(i)+wvis_liq(i)*cell%cpl(i)/prt
            wcd_gas(i)=cell%lcondg(i)+wvis_gas(i)*cell%cpg(i)/prt
         ENDIF
      ENDDO
!
!.....Effective thermal conductivity
!
      DO i=1,ncell_fluid
          cell%condg(i)=cell%lcondg(i)+cell%tviscosg(i)*cell%cpg(i)/prt
          cell%condl(i)=cell%lcondl(i)+cell%tviscosl(i)*cell%cpl(i)/prt
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_cond
