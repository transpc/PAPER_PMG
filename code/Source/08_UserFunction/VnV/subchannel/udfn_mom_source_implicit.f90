!
      SUBROUTINE udfn_subchannel_mom_source
!
!     Define the momentum friction Models
!
      USE Zporous , ONLY: s_subchannel_fric,s_subchannel_mixing
!
      IMPLICIT NONE
!      
!.....Friction model (Form Loss, corss/axial-flow, Spacer Grid, Mixing Vane)
!        - implicit/semi-implicit: fric_model(i,ix) is calculated and then goes to 3x3 matrix
!        - explicit: src_gas,src_liq are calculated and then goes to source term
!
      IF(s_subchannel_fric.eq.'aniso_fric_sem' .or. &
         s_subchannel_fric.eq.'aniso_fric_imp'      )then
!
         CALL udfn_subchannel_mom_source_implicit_main      
!
      ELSEIF(s_subchannel_fric.eq.'aniso_fric_exp')THEN
!
         CALL udfn_subchannel_mom_source_explicit_main      
!
      ENDIF 
!
!.....Turbulent mixing model (EM/EVVD)
!     source term (src_liq, src_gas are calculated)
!
      IF(s_subchannel_mixing=='em' .or. &
         s_subchannel_mixing=='evvd'     ) THEN
!         
         CALL udfn_subchannel_mom_source_turbmixing
!
      ELSE
!         WRITE(97,*) 'NO subchannel turbulent mixing model is selected !!!'
!         WRITE(* ,*) 'NO subchannel turbulent mixing model is selected !!!'
!         STOP
!         
      ENDIF
!
      END SUBROUTINE udfn_subchannel_mom_source
!
!======================================================================
!======================================================================
!
      SUBROUTINE udfn_subchannel_mom_source_implicit_main
!
!     Define the implicit momentum friction Models 
!        - fric_model_liq(i,ix), cell%vfwl_x(i)/y(i)/z(i) are treated implicitly
!
      USE VOL_DATA    , ONLY: cell
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim
      USE Ztimecon    , ONLY: time,itim
      USE Znum_cell   , ONLY: i_neigh,neigh
      USE Zconst1     , ONLY: vv_prob
      USE Zconst2     , ONLY: hydraulicd,sl,stime_vflat
      USE Zcoord1     , ONLY: xloc
      USE Zcoord2     , ONLY: cell_leng      
      USE Zcoord3     , ONLY: porosity,volp
      USE Zporous     , ONLY: fric_model_liq,fric_model_gas
      USE Zporous     , ONLY: kloss_grid
      USE Zporous     , ONLY: s_ij_non_i,s_ij_non_k     
      USE Zporous     , ONLY: s_subchannel_fric_cross,s_subchannel_fric_axial, &
                              s_2p_multiplier,l_spacer_grid,l_mixing_vane
      USE Zvector     , ONLY: vl_o,vg_o,ul_o,ug_o
      USE Zvec_geo    , ONLY: xn_nf
      USE Zporous     , ONLY: sg_loc,mv_loc,chn_type,n_sg,h_sg,n_mv,h_mv
!
      IMPLICIT NONE
!            
!.....External function
      INTEGER get_nf_number_j
!.....Local variables
      INTEGER i,ix,j0,nf_number
      REAL(8) hd
      REAL(8) rel,reg,fric !,kloss_grid
      REAL(8) fric1,fric2,KK,bb !nakayama
      REAL(8) TPM     
      REAL(8) s_ij,xn1,xn2
      ! OPR1000 rod-scale (mixing vane)
      INTEGER j,k,k1,i1
      REAL(8) tr
      REAL(8) F_latconv
      REAL(8) cell_bot,cell_top
!      
!.....Subchannel pressure drop model
!
      DO i=1,ncell_fluid   
         hd=hydraulicd(i) 
 
         ! Vessel problem (ONLY core considered)
         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel'.or.vv_prob.eq.'KSMR')then     !PSH
            IF(chn_type(i).eq.0)CYCLE
         ENDIF     
!
         DO ix=1,ndim
!
!...........Friction Model: Cross Flows
!
            IF(ix.ne.ndim) THEN   
!
               SELECT CASE(s_subchannel_fric_cross)
!
!              MATRA Form Loss         
               CASE('simple_formloss')
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-(kloss_grid/2.d0)*sl(i,ix)*cell%rhol(i)*abs(vl_o(i,ix))*cell%alphal(i)!*ul_o(i)*cell%alphal(i) !*vl_o(i,ix)
                  fric_model_gas(i,ix)=fric_model_gas(i,ix)-(kloss_grid/2.d0)*sl(i,ix)*cell%rhog(i)*abs(vg_o(i,ix))*cell%alphag(i)!*ug_o(i)*cell%alphag(i) !*vl_o(i,ix)
               END SELECT
!            
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
                     TPM=(1.d0-cell%quala(i))*(1.d0-cell%quala(i))/(1.d0-cell%alphag(i))**1.42d0
                  ELSEIF(cell%alphag(i).le.0.9d0) THEN
                     TPM=0.478d0*(1.d0-cell%quala(i))*(1.d0-cell%quala(i))/(1.d0-cell%alphag(i))**2.2d0
                  ELSE
                     TPM=1.73d0*(1.d0-cell%quala(i))*(1.d0-cell%quala(i))/(1.d0-cell%alphag(i))**1.64d0
                  ENDIF
               END SELECT
!
!..............Friction Model: Axial Flow
!
               SELECT CASE(s_subchannel_fric_axial)
!
!              MATRA               
               CASE('matra')
                  Rel=DMAX1(1.d0,(cell%rhol(i)*ul_o(i)*hd/cell%lviscosl(i)))
                  CALL mom_wall_matra(Rel,fric)
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-TPM*fric/(2.d0*hd)*cell%rhol(i)*abs(vl_o(i,ix))*cell%alphal(i)!*ul_o(i)*cell%alphal(i) !*vl_o(i,ix)

                  Reg=DMAX1(1.d0,(cell%rhog(i)*ug_o(i)*hd/cell%lviscosg(i)))
                  CALL mom_wall_matra(Reg,fric)
                  fric_model_gas(i,ix)=fric_model_gas(i,ix)-TPM*fric/(2.d0*hd)*cell%rhog(i)*abs(vg_o(i,ix))*cell%alphag(i)!*ug_o(i)*cell%alphag(i) !*vl_o(i,ix)
!               
!              Chandesris Friction
               CASE('chandesris') 
                  Rel=DMAX1(1.d0,cell%rhol(i)*ul_o(i)*hd/cell%lviscosl(i))
                  fric=Rel**(-0.25d0)/4.d0                 
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-TPM*fric/(2.d0*hd)*cell%rhol(i)*ul_o(i) !*vl_o(i,ix)                  
!
!              Nakayama               
               CASE('nakayama')
                  KK=porosity(i)*hd*hd/32.d0
                  bb=0.3164d0/(2.d0*hd*porosity(i)*porosity(i)*(hd*DMAX1(1.e-8,ul_o(i))/  &
                     cell%lviscosl(i)*cell%rhol(i))**0.25d0)                   
                  fric1=porosity(i)*cell%lviscosl(i)/KK
                  fric2=porosity(i)*porosity(i)*bb*ul_o(i)*cell%rhol(i)
                  fric=fric1+fric2
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-TPM*fric !*vl_o(i,ix)
!
!              Otherwise 
               CASE('takeda')
                  fric1=150.d0*cell%lviscosl(i)*((1.d0-porosity(i))/porosity(i)/hd)**2.d0
                  fric2=1.75d0*cell%rhol(i)*(1.d0-porosity(i))/porosity(i)/hd*ul_o(i)
                  fric=porosity(i)*(fric1+fric2)
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-TPM*fric !*vl_o(i,ix)               
!                    
!              CTF               
               CASE('ctf')
                  Rel=DMAX1(1.d0,(cell%rhol(i)*ul_o(i)*hd/cell%lviscosl(i)))
                  CALL mom_wall_ctf(Rel,fric)
                  fric_model_liq(i,ix)=fric_model_liq(i,ix)-TPM*fric/(2.d0*hd)*cell%rhol(i)*dabs(vl_o(i,ix))
                  Reg=DMAX1(1.d0,(cell%rhog(i)*ug_o(i)*hd/cell%lviscosg(i)))
                  CALL mom_wall_ctf(Reg,fric)
                  fric_model_gas(i,ix)=fric_model_gas(i,ix)-TPM*fric/(2.d0*hd)*cell%rhog(i)*dabs(vg_o(i,ix))
!                    
               END SELECT
            ENDIF 
         ENDDO
      ENDDO
!
!.....Additional Friction Model: Spacer Grid
!        - ONLY axial flow considered (mostly ndim=3)
!
      DO i=1,ncell_fluid
         cell%vfwg_x(i)=0.d0
         cell%vfwg_y(i)=0.d0
         cell%vfwl_x(i)=0.d0
         cell%vfwl_y(i)=0.d0
         cell%vfwg_z(i)=0.d0
         cell%vfwl_z(i)=0.d0
      ENDDO

      IF(l_spacer_grid)then

         DO i=1,ncell_fluid   
            ! Vessel problem (ONLY core considered)
            IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')then
               IF(chn_type(i).eq.0)CYCLE
            ENDIF      
!
!           OPR1000, APR1400_fullcore
            IF(vv_prob.eq.'APR1400_fullcore'                    .or. &
               vv_prob.eq.'APR1400_fullcore_modmesh01'          .or. &
               vv_prob.eq.'APR1400_fullcore_modmesh02_rv'       .or. &
               vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'       .or. &
               vv_prob.eq.'OPR1000_single_assem'                .or. &
               vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel'      )THEN
                   
              !Lamp until stime_vflat (=0.1sec)
               IF(time.le.stime_vflat) THEN 
                  tr=time/stime_vflat
               ELSE
                  tr=1.0
               ENDIF
                   
               kloss_grid=0.57d0
               DO k1=1,n_sg
                  cell_bot=xloc(i,ndim)-0.5d0*cell_leng(i,ndim)
                  cell_top=xloc(i,ndim)+0.5d0*cell_leng(i,ndim)
                  IF((sg_loc(k1)+0.5d0*h_sg(k1)).ge.cell_bot .and. (sg_loc(k1)+0.5d0*h_sg(k1)).lt.cell_top)THEN
                     IF(l_mixing_vane)then
                        cell%vfwl_z(i)=-kloss_grid/(2.d0)*cell%rhol(i)*dabs(vl_o(i,3))*vl_o(i,3)*tr
                     ELSE
                       fric_model_liq(i,ndim)=fric_model_liq(i,ndim)-(kloss_grid/(2.d0*h_sg(k1)))*cell%rhol(i)*ul_o(i)*tr
                     ENDIF
                  ENDIF
               ENDDO                  
            ENDIF !ELSEIF(vv_prob.eq.
         ENDDO !DO i=1,ncell_fluid   
      ENDIF !IF(l_spacer_grid)then
!
!.....Additional Friction Model: Mixing vane
!        - ONLY cross flow considered (mostly ndim.NE.3)
!
      IF(l_mixing_vane) THEN
        !Lamp until stime_vflat (=0.1sec)
         IF(time.le.stime_vflat) THEN 
            tr=time/stime_vflat
         ELSE
            tr=1.0
         ENDIF

         DO i=1,ncell_fluid   
!
!           rod-scale 
            IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel' .or. &
               vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'        .or. &
               vv_prob.eq.'OPR1000_single_assem'                       )THEN

             ! core region ONLY
               IF(chn_type(i).eq.0)CYCLE

               F_latconv=0.27d0 !original

               cell_bot=xloc(i,ndim)-0.5d0*cell_leng(i,ndim)
               cell_top=xloc(i,ndim)+0.5d0*cell_leng(i,ndim)

               IF(chn_type(i).eq.5)THEN     !guide_side
                  F_latconv=F_latconv*0.576d0
               ELSEIF(chn_type(i).eq.6)THEN !guide_cor
                  F_latconv=F_latconv*1.05d0
               ENDIF

               IF(chn_type(i).le.6)THEN
                  DO k1=1,n_mv
                     IF((mv_loc(k1)+0.5d0*h_mv).ge.cell_bot .and. &
                        (mv_loc(k1)+0.5d0*h_mv).le.cell_top        )THEN
                        IF(itim.gt.2)then
                           j0=i_neigh(i)-1
                           DO j=i_neigh(i),i_neigh(i+1)-1
                              nf_number=get_nf_number_j(j)
                              IF(nf_number.eq.0) THEN          
                                 CALL get_vector_disp(j-j0,i,i1)
                                 IF(i1.lt.0) THEN
                                    s_ij=s_ij_non_k(-i1)
                                    xn1=-xn_nf(-i1,1)
                                    xn2=-xn_nf(-i1,2)
                                 ELSE
                                    s_ij=s_ij_non_i(i1)
                                    xn1=xn_nf(i1,1)
                                    xn2=xn_nf(i1,2)
                                 ENDIF
                                 k=neigh(j)
                                 IF(chn_type(i).eq.2)then
                                    IF(chn_type(k).eq.2)then
                                       ! horizontal
                                       IF(cell_leng(i,1).gt.cell_leng(i,2))then
                                          IF(dabs(xn2).gt.0.5d0)then
                                             cycle
                                             cell%vfwl_x(i)=0.d0
                                             cell%vfwl_y(i)=0.d0
                                          ENDIF
                                       ! vertical
                                       ELSE
                                          IF(dabs(xn1).gt.0.5d0)then
                                             cycle
                                             cell%vfwl_x(i)=0.d0
                                             cell%vfwl_y(i)=0.d0
                                          ENDIF
                                       ENDIF
                                    ENDIF 
                                 ENDIF
                                 IF(chn_type(i).eq.3)then
                                    IF(chn_type(k).eq.3)then
                                       cycle
                                       cell%vfwl_x(i)=0.d0
                                       cell%vfwl_y(i)=0.d0
                                    ENDIF
                                 ENDIF
                                 IF(chn_type(k).le.6 .and. chn_type(k).ne.4 .and. chn_type(k).ne.5)THEN
                                    IF(abs(xn1).gt.0.5d0)THEN
                                       cell%vfwl_x(i)=((F_latconv)**2.d0 * abs(vl_o(i,ndim))*vl_o(i,ndim)*cell%rhol(i)*s_ij*h_mv/volp(i)) *tr
                                    ELSEIF(abs(xn2).gt.0.5d0)THEN
                                       cell%vfwl_y(i)=((F_latconv)**2.d0 * abs(vl_o(i,ndim))*vl_o(i,ndim)*cell%rhol(i)*s_ij*h_mv/volp(i)) *tr
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDDO
                        ENDIF !itim.gt.2
                     ENDIF ! IF((mv_loc(k1)+ : Mixing Vane Location
                  ENDDO ! DO k1=1,n_mv
               ENDIF ! subchannel to include mixing vane model
            ENDIF !IF(vv_prob.eq.
!
         ENDDO ! DO i=1,ncell_fluid
      ENDIF
!
      END SUBROUTINE udfn_subchannel_mom_source_implicit_main
!
!======================================================================
!======================================================================
!
!
      SUBROUTINE udfn_subchannel_mom_source_turbmixing
!
!     Calculates turbulent mixing term using EM or EVVD model
!        - results are summed into source term of momentum equation
!        - CALL either EV or EVVD model and then calculates src_liq, src_gas
!
      USE Zinterface
      USE Zmpi        , ONLY: ncell_fp
      USE Zzone       , ONLY: ncell_fluid
      USE Zcore       , ONLY: np
      USE Zparam      , ONLY: ndim
      USE Zvec_param  , ONLY: nf_non
      USE Zcoord2     , ONLY: cell_leng
      USE Zcoord3     , ONLY: volp
      USE Zporous     , ONLY: s_ij_non_i,s_ij_non_k,cell_area, &
                              s_gapij_non_i,s_gapij_non_k,     &
                              s_subchannel_mixing
      USE Zporous     , ONLY: chn_type
      USE Zrv_model   , ONLY: rv_ht_str
!
      IMPLICIT NONE
!            
!.....Local variables
      INTEGER :: i
      LOGICAL,SAVE :: init_sij=.true.
!      
!.....Define s_ij
!
      IF(init_sij) THEN

         SELECT CASE(s_subchannel_mixing)
!
!        EM  
         CASE('em')
            ALLOCATE(s_gapij_non_i(nf_non),s_gapij_non_k(nf_non))     
            CALL udfn_s_gapij
!      
            init_sij=.false.
!               
!        EVVD
         CASE('evvd') 
            ALLOCATE(s_ij_non_i(nf_non),s_ij_non_k(nf_non))     
            CALL udfn_sij
!
!...........cell_area for EVVD model
!     
            ALLOCATE(cell_area(ncell_fp))
            cell_area=0.d0

            DO i=1,ncell_fluid
               IF(rv_ht_str.ne.0)THEN
                  IF(chn_type(i).eq.0)CYCLE
               ENDIF
               cell_area(i)=volp(i)/cell_leng(i,ndim)
            ENDDO  
            IF(np.gt.1) CALL communicate_1d(cell_area)
!      
            init_sij=.false.
!      
         END SELECT
!
      ENDIF !IF(init_sij)
!
!.....EM/EVVD model for turbulence mixing
!
      SELECT CASE(s_subchannel_mixing)
!
!     EM  
      CASE('em')
         CALL udfn_mom_source_em
!               
!     EVVD
      CASE('evvd') 
         CALL udfn_mom_source_evvd
!
      END SELECT
!
      END SUBROUTINE udfn_subchannel_mom_source_turbmixing
