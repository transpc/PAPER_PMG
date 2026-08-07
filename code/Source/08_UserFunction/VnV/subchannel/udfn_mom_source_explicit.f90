!
      SUBROUTINE  udfn_subchannel_mom_source_explicit_main
!
!     Modifies the momentum source terms at the free surface cells
!
      IMPLICIT NONE
!                  
!.....Friction model for Porous (Explicit)
!
      CALL udfn_subchannel_mom_source_explicit

      RETURN
      ENDSUBROUTINE  udfn_subchannel_mom_source_explicit_main
!
!======================================================================
!======================================================================
!
      SUBROUTINE  udfn_subchannel_mom_source_explicit
!
!     Modifies the momentum source terms at the free surface cells
!
      USE VOL_DATA             
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: ndim
      USE Zconst1       , ONLY: vv_prob                     !PSH
      USE Zconst2       , ONLY: hydraulicd,sl
      USE Zcoord1       , ONLY: xloc
      USE Zcoord2       , ONLY: cell_leng      
      USE Zcoord3       , ONLY: porosity
      USE Zporous       , ONLY: kloss_grid,kloss_cross,n_sg,h_sg,sg_loc
      USE Zm_src        , ONLY: src_gas,src_liq
      USE Zporous       , ONLY: s_subchannel_fric_cross,s_2p_multiplier,s_subchannel_fric_axial, &
                                l_spacer_grid,chn_type      !PSH
      USE Zvector       , ONLY: vl_o,vg_o,ul_o
!
      IMPLICIT NONE
!            
      INTEGER i,ix,k
!      
      REAL(8) hd,rel,reg,fric !,kloss_grid
      REAL(8) fric1,fric2,KK,bb !nakayama   
      REAL(8) velo,velo_tmp
      REAL(8) velog,velol,velog_tmp,velol_tmp      
      REAL(8) TPM
     !spacer grid
      REAL(8) cell_bot,cell_top
!
!-----------------------------------------------------------------------------
      DO i=1,ncell_fluid
         ! Vessel problem (ONLY core considered)
         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel'.or.vv_prob.eq.'KSMR')then         !PSH
            IF(chn_type(i).eq.0)CYCLE
         ENDIF 
         hd=hydraulicd(i) 
         velo_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
         IF(ndim.eq.3) velo_tmp=velo_tmp+vl_o(i,3)*vl_o(i,3)
         velo=DSQRT(velo_tmp)
!
         velol_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
         IF(ndim.eq.3) velol_tmp=velol_tmp+vl_o(i,3)*vl_o(i,3)
         velol=DSQRT(velol_tmp)
!
         velog_tmp=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)
         IF(ndim.eq.3) velog_tmp=velog_tmp+vg_o(i,3)*vg_o(i,3)
         velog=DSQRT(velog_tmp)         
!                  
         DO ix=1,ndim
!
!...........Friction Model: Cross Flows
!
            IF(ix.ne.ndim) THEN 
               SELECT CASE(s_subchannel_fric_cross)  
!
!              MATRA Form Loss         
               CASE('simple_formloss')                  
                  src_gas(i,ix)=src_gas(i,ix)-(kloss_cross/2.d0)*cell%rhog(i)*velog*vg_o(i,ix)*sl(i,ix)*cell%alphag(i)
                  src_liq(i,ix)=src_liq(i,ix)-(kloss_cross/2.d0)*cell%rhol(i)*velol*vl_o(i,ix)*sl(i,ix)*cell%alphal(i)                  
               END SELECT
            ENDIF   
!
!...........Friction Model: Axial Flow
!
            IF(ix.eq.ndim) THEN
!
!..............Two-Phase Multiplier: none, default, armand model
!
               SELECT CASE(s_2p_multiplier)

               CASE('none')
                  TPM=1.d0 
               CASE('default')
                  TPM=cell%rhol(i)/cell%rhom(i)
               CASE('armand')
                  IF(cell%alphag(i).le.0.6d0) THEN
                     TPM=(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**1.42
                  ELSEIF(cell%alphag(i).le.0.9d0) THEN
                     TPM=0.478d0*(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**2.2
                  ELSE
                     TPM=1.73d0*(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**1.64
                  ENDIF
               END SELECT
!
!..............Friction Model: Axial Flow
!
               SELECT CASE(s_subchannel_fric_axial)
!
!              MATRA               
               CASE('matra')
                  Rel=DMAX1(1.0d0,(cell%rhol(i)*velol*(1.0d0*hd)/cell%lviscosl(i)))
                  CALL mom_wall_matra(Rel,fric)
                  src_liq(i,ix)=src_liq(i,ix)-TPM*fric/(2.0d0*hd)*cell%rhol(i)*velol*vl_o(i,ix)*cell%alphal(i)
!
                  Reg=DMAX1(1.0d0,(cell%rhog(i)*velog*(1.0d0*hd)/cell%lviscosg(i)))
                  CALL mom_wall_matra(Reg,fric)
                  src_gas(i,ix)=src_gas(i,ix)-TPM*fric/(2.0d0*hd)*cell%rhog(i)*velog*vg_o(i,ix)*cell%alphag(i) 
!                  
!              Chandesris Friction
               CASE('chandesris') 
                  Rel=DMAX1(1.0d0,cell%rhol(i)*velo*hd/cell%lviscosl(i))
                  fric=Rel**(-0.25d0)/4.d0                 
                  src_liq(i,ix)=src_liq(i,ix)-TPM*1.d0*fric/(2.0d0*hd)*cell%rhol(i)*velo*vl_o(i,ix)
!              Nakayama               
               CASE('nakayama')
                  KK=porosity(i)*hd*hd/32.d0 
                  bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*(hd*DMAX1(1.e-8,velo)/cell%lviscosl(i)*cell%rhol(i))**0.25d0) 
                  fric1=porosity(i)*cell%lviscosl(i)/KK
                  fric2=porosity(i)*porosity(i)*bb*velo*cell%rhol(i)
                  fric=fric1+fric2
                  src_liq(i,ix)=src_liq(i,ix)-TPM*fric*vl_o(i,ix)
               END SELECT
!               
            ENDIF
         ENDDO
         
      ENDDO  
!
!.....Additional Friction Model: Spacer Grid
!        - ONLY axial flow considered (mostly ndim=3)
!
      IF(l_spacer_grid) THEN
!
         DO i=1,ncell_fluid
            DO k=1,n_sg
               cell_bot=xloc(i,ndim)-0.5d0*cell_leng(i,ndim)
               cell_top=xloc(i,ndim)+0.5d0*cell_leng(i,ndim)
               IF(sg_loc(k).ge.cell_bot .and. sg_loc(k).lt.cell_top)THEN
                  src_liq(i,ndim)=src_liq(i,ndim)-(kloss_grid/(2.0d0*h_sg(k)))*cell%rhol(i)*ul_o(i)*vl_o(i,ndim)
               ENDIF
            ENDDO                  
         ENDDO                  
!           
      ENDIF 
!
!-----------------------------------------------------------------------------
      RETURN
      END SUBROUTINE  udfn_subchannel_mom_source_explicit
