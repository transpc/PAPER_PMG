!
      SUBROUTINE udfn_heat_wallHTC_porous
!
!     This SUBROUTINE calculates heat trasnfer coefficient using basic heat 
!     transfer correlation instead of heat_partition_porous in porous cells.
!     (Dittus-Boelter, Colburn, etc.)
!      
      USE VOL_DATA                 
      USE SOLID_DATA    , ONLY:solid
      USE Zparam        , ONLY:pi
      USE Zconst1       , ONLY: vv_prob
      USE Zcoord3       , ONLY:aporous
      USE Znum_cell     , ONLY:n_fluid
      USE Zvector       , ONLY:vl_o
      USE Zqvol         , ONLY:qporous_liq, qporous_gas
      USE Zzone         , ONLY:ncell_cond,nmaterial_c
!
      IMPLICIT NONE
!
      INTEGER i,ii
!      
      REAL(8) h_mac,q1,qb,Re,Pr,twall_porous
      REAL(8) D_tube,Pit_tube,D_h,delvl,hconvl
!      
      IF(vv_prob.eq.'stern'.or.vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro') THEN
!   
!DIR$ SIMD
         DO i=1,ncell_cond
            IF(nmaterial_c(i).lt.0)THEN
!
            ii=n_fluid(i)
!
!........Assign wall temperature
!
            twall_porous=solid%tsol(i)
!
!........Parameters of correlations
!
            IF(vv_prob.eq.'stern') THEN
            D_tube=0.033d0
            Pit_tube=0.0715d0
               D_h=(4.0d0*Pit_tube*Pit_tube-pi*D_tube*D_tube)/(Pit_tube*4.0+pi*D_tube)
            ELSEIF(vv_prob.eq.'atlas_mc_porous')THEN
               D_tube=9.5d-3
               Pit_tube=12.85d-3
               D_h=12.63d-3
            ELSEIF(vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro')THEN
               D_tube=9.5d-3
               Pit_tube=12.85d-3
               D_h=12.63d-3         
            ENDIF
            delvl=dsqrt(DOT_PRODUCT(vl_o(ii,:),vl_o(ii,:)))
            Pr=cell%lviscosl(ii)*cell%cpl(ii)/cell%lcondl(ii)
            Re=cell%rhol(ii)*delvl*D_h/cell%lviscosl(ii)
!        IF(.true.)then
           !Dittus-Boelter
             h_mac=0.023d0*Re**0.8d0*Pr**0.4d0 * (cell%lcondl(ii)/D_h)
!        ELSE  
           !Colburn
!          h_mac=0.33d0*Re**0.6d0*Pr**0.3d0 * (cell%lcondl(ii)/D_h)	
!        ENDIF  
!!!!
!!!!........heat transfer coefficient
!!!!
!!!            hconvl=h_mac
!!!            hconvl=DMAX1(1.0d2,DMIN1(hconvl,1.d8))
!!!            hconvg=0.0d0
!!!!
!!!!........Heat partition
!!!!
!!!            q1=hconvl*Aporous(ii)*(twall_porous-cell%tl(ii)) * 2.0d0
!!!            qb=0.0d0
!!!!
!!!!........Assign partitioned heat to equation variables
!!!!
!!!            qporous_liq(ii)=qporous_liq(ii)+q1
!!!            qporous_gas(ii)=qporous_gas(ii)+qb
!!!   !        hconvl=h_mac 
            q1=h_mac*Aporous(ii)*(twall_porous-cell%tl(ii)) * 2.0d0
            qb=0.0d0          
            qporous_liq(ii)=q1
            qporous_gas(ii)=qb
            hconvl=h_mac    
            ENDIF
         ENDDO
!      
      ELSEIF(vv_prob.eq.'apr1400_lbloca') THEN
!   
         DO i=1,ncell_cond
            IF(nmaterial_c(i).lt.0)THEN
!
            ii=n_fluid(i)
!
!........Assign wall temperature
!
            twall_porous=solid%tsol(i)
!
!........Parameters of correlations
!
            D_tube=0.033d0
            Pit_tube=0.0715d0
            D_h=0.321d0  !dc distance
            delvl=dsqrt(DOT_PRODUCT(vl_o(ii,:),vl_o(ii,:)))
            Pr=cell%lviscosl(ii)*cell%cpl(ii)/cell%lcondl(ii)
            Re=cell%rhol(ii)*delvl*D_h/cell%lviscosl(ii)
!        IF(.true.)then
           !Dittus-Boelter
             h_mac=0.023d0*Re**0.8d0*Pr**0.4d0 * (cell%lcondl(ii)/D_h)
!        ELSE  
!          !Colburn
!          h_mac=0.33d0*Re**0.6d0*Pr**0.3d0 * (cell%lcondl(ii)/D_h)	
!        ENDIF  
!
!........Heat partition
!
            h_mac=cell%condg(ii)*cell%alphag(ii)+cell%condl(ii)*(1.0d0-cell%alphag(ii)) 
            q1=h_mac*Aporous(ii)*(twall_porous-cell%tl(ii))
            qb=0.0d0
   !
!........Assign partitioned heat to equation variables
!
            qporous_liq(ii)=qporous_liq(ii)+q1*(1.0d0-cell%alphag(ii))
            qporous_gas(ii)=qporous_gas(ii)+q1*cell%alphag(ii)
!               
            ENDIF
         ENDDO
!         
      ELSE
!      
         WRITE(*,"(11x,a,a,a)")'udfn_heat_wallHTC_porous is not specified for ',vv_prob,' !'  
!      
      ENDIF
!

      RETURN
      END SUBROUTINE udfn_heat_wallHTC_porous
