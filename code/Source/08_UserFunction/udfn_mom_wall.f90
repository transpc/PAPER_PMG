!
      SUBROUTINE udfn_mom_wall
!
!     User-defined wall friction (only meaningful when 'udfl_mom_wall' is used.)
!     
      USE Zinterface
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid,nzone
      USE Zcore           , ONLY: np
      USE Zparam          , ONLY: ndim
      USE Znum_cell       , ONLY: istart_nf
      USE Zvec_index      , ONLY: right_non
      USE Znum_cell       , ONLY: istart_nf,istart_nb1,ia_nb,icell_nb
      USE Ztimecon        , ONLY: time
      USE Zndforce        , ONLY: relax_cd
      USE Zb_condition    , ONLY: vb_liq
      USE Zconst1         , ONLY: vv_prob
      USE Zcoord1         , ONLY: xloc
      USE Zcoord3         , ONLY: porosity
      USE Zfluidic_device , ONLY: flow_sp,flow_fd,k_dp,k_sp
      USE Zmodel          , ONLY: vfwl_k,fsar
      USE Zmodel          , ONLY: s_wall_fric
      USE Zvector         , ONLY: ul_o,ug_o,vl_o
      USE Zbc_index       , ONLY: npb
      USE Zvec_geo        , ONLY: sv_nf,svp_nf
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,ii,kk,nb
      INTEGER :: nf_number,istart,len,istart1,i0,i1
      LOGICAL, SAVE :: INITIAL=.TRUE.
      LOGICAL, SAVE :: initial_kdp=.TRUE.
!      
!.....Fluidic device
!
      INTEGER, SAVE :: num_sp,num_fd,num_dc,n_sp(40),n_fd(50),n_dc(200)
! 
      REAL(8) :: vlv_fr
      REAL(8) :: k_sp_high,k_fd,k_dp_high,k_dp_low,flow_sp_sin,flow_sp_cos,tan_theta,theta,fr,k_dp_i
      REAL(8) :: vfwg_i,vfwl_i,HydraulicD
      REAL(8) :: fric,Rel,Reg
!
!.....Stern Moderator
!
      REAL(8) :: D_tube,Pit_tube, Vfs, Vm, Re_tube, Kloss, delvl
      REAL(8) :: DelP_cf, DelP_axi
!
!.....OPR1000 friction
!     
      REAL(8) :: y0,wid_fa,widx0,widy0,widx1,widy1
      INTEGER,ALLOCATABLE,SAVE :: assem_ix(:),assem_iy(:)
!
!.....SMR
!      
      IF(vv_prob.eq.'KSMR'.or.vv_prob.eq.'KSMR-SG-porous'.or.vv_prob.eq.'KSMR-SG-pid'.or.vv_prob.eq.'KSMR-PZR') CALL PressureDrop()
      IF(vv_prob.eq.'Nuscale-RVV'.or.vv_prob.eq.'Nuscale-PZR') CALL PressureDrop_NuScale()      
!
!.....PAFS-POOL, SMALL-POOL, PAFS-3D: pool incuding porous cells 
!
      IF(vv_prob.eq.'PAFS-POOL') THEN
         DO i=1,ncell_fluid         
            vfwg_i=0.d0
            vfwl_i=0.d0
!
!...........Define the hydraulic diameter
!
            HydraulicD=0.11d0  !This value is geometry-dependent. If a model is changed, the value should be changed.
!
!...........Define wall drag coeffs. (vfwl_i,vfwg_i)
!            
            IF(cell%regime(i).eq.11 .or. cell%regime(i).eq.12 .or. cell%regime(i).eq.21 .or.cell%regime(i).eq.22 .or. cell%regime(i).eq.3)THEN
               Rel=MAX(1.d0,(cell%rhol(i)*ul_o(i)*(2.d0*HydraulicD)/cell%lviscosl(i)))
               IF(s_wall_fric.eq.'kakac') CALL mom_wall_kakac(Rel,fric)
               vfwl_i=fric/(2.d0*HydraulicD)*cell%rhol(i)*ul_o(i)
            ELSEIF(cell%regime(i).eq.13)THEN
               Reg=MAX(1.d0,cell%rhog(i)*ug_o(i)*(2.d0*HydraulicD)/cell%lviscosg(i))
               IF(s_wall_fric.eq.'kakac') CALL mom_wall_kakac(Reg,fric)
               vfwg_i=fric/(2.d0*HydraulicD)*cell%rhog(i)*ug_o(i)
               vfwl_i=cell%alphal(i)*fric/(2.d0*HydraulicD)*cell%rhol(i)*ul_o(i)
            ENDIF
!
!...........Define the maximum limit of wall drag coeffs.
!        
            IF(vv_prob.eq.'PAFS-POOL') THEN
               IF(nzone(i).ge.2.d0 .and. nzone(i).le.3.d0)THEN
                  vfwl_i=MIN(200.d0*vfwl_i,vfwl_i+0.5d0/HydraulicD*cell%rhol(i)*ul_o(i)*cell%alphal(i))
                  vfwg_i=MIN(200.d0*vfwg_i,vfwg_i+0.5d0/HydraulicD*cell%rhog(i)*ug_o(i)*cell%alphag(i)) 
               ELSEIF(xloc(i,2).le.1.8d0 .and. xloc(i,1).le.4.d0)THEN
                  vfwl_i=MIN(200.d0*vfwl_i,vfwl_i+0.5d0/HydraulicD*cell%rhol(i)*ul_o(i)*cell%alphal(i))
                  vfwg_i=MIN(200.d0*vfwg_i,vfwg_i+0.5d0/HydraulicD*cell%rhog(i)*ug_o(i)*cell%alphag(i))
               ENDIF 
            ENDIF
!
!...........Modifiy vfwl_i to stabilize the free surface response (Should be changed w.r.t. problems)
!
            IF(cell%regime(i).eq.3)THEN             
               vfwl_i=50.d0*vfwl_i
            ENDIF 
!
            cell%vfwg(i)=relax_cd*cell%vfwg(i)+(1.d0-relax_cd)*vfwg_i
            cell%vfwl(i)=relax_cd*cell%vfwl(i)+(1.d0-relax_cd)*vfwl_i
!            
        ENDDO
      ENDIF  
!
!.....stern
!
      IF(vv_prob.eq.'stern') THEN
         DelP_cf=0.d0
         DelP_axi=0.d0
! bug vfwg_i not computed no idea if zero value is OK???
         vfwg_i=0.d0
         DO i=1,ncell_fluid
            delvl=SQRT(dot_product(vl_o(i,:),vl_o(i,:)))
            IF(porosity(i).eq.1.d0)THEN
               vfwl_i=cell%rhol(i)*delvl*10.d0*(1.d0-porosity(i))
            ELSEIF(porosity(i).ne.1.d0)THEN !Cross flow
               D_tube=0.033d0
               Pit_tube=0.0715d0
               Vm=vb_liq(1,2)
               Vfs=ul_o(i)*porosity(i)
               Re_tube=cell%rhol(i)*Vm*D_tube/cell%lviscosl(i)
               Kloss=4.54d0/Pit_tube*Re_tube**(-0.172d0)
               DelP_cf=Kloss*(cell%rhol(i)*Vfs/2.d0)*porosity(i)
               vfwl_i=(DelP_cf+DelP_axi)
            ENDIF
            cell%vfwg(i)=relax_cd*cell%vfwg(i)+(1.d0-relax_cd)*vfwg_i
            cell%vfwl(i)=relax_cd*cell%vfwl(i)+(1.d0-relax_cd)*vfwl_i          
        ENDDO
      ENDIF  
!
!.....fluidic_device
!
      IF(vv_prob.eq.'fluidic_device')THEN
!
         IF(INITIAL)THEN
            num_sp=0
            num_fd=0
            num_dc=0
            DO i=1,ncell_fluid
!
               IF(xloc(i,1).gt.1.04115.and.xloc(i,1).lt.1.38995.and.  &
                  xloc(i,2).gt.1.04115.and.xloc(i,2).lt.1.38995)THEN
                  IF(xloc(i,ndim).gt.0.32.and.xloc(i,ndim).lt.0.53)THEN
                     num_sp=num_sp+1
                     n_sp(num_sp)=i
                  ENDIF
                  IF(xloc(i,ndim).gt.0.0.and.xloc(i,ndim).lt.0.125)THEN
                     num_dc=num_dc+1
                     n_dc(num_dc)=i
                  ENDIF
               ENDIF
               IF(xloc(i,ndim).gt.0.0.and.xloc(i,ndim).lt.0.25)THEN
                  IF(xloc(i,1).gt.1.12972.and.xloc(i,1).lt.1.30138)THEN
                     IF((xloc(i,2).gt.0.61555.and.xloc(i,2).lt.1.04115).or.   &
                        (xloc(i,2).gt.1.38995.and.xloc(i,2).lt.1.81555))THEN
                        num_fd=num_fd+1
                        n_fd(num_fd)=i
                     ENDIF
                  ENDIF
                  IF(xloc(i,2).gt.1.12972.and.xloc(i,2).lt.1.30138)THEN
                     IF((xloc(i,1).gt.0.61555.and.xloc(i,1).lt.1.04115).or.   &
                        (xloc(i,1).gt.1.38995.and.xloc(i,1).lt.1.81555))THEN
                        num_fd=num_fd+1
                        n_fd(num_fd)=i
                     ENDIF
                  ENDIF
               ENDIF
               IF(xloc(i,ndim).lt.0.d0.and.xloc(i,ndim).gt.-0.176d0)THEN
                  num_dc=num_dc+1
                  n_dc(num_dc)=i
               ENDIF
            ENDDO
            INITIAL=.FALSE.
         ENDIF
!
!........Adjust surface vectors of pressure boundary cells to simulate valve characteristics
!
         vlv_fr=0.1d0*time
         vlv_fr=MIN(1.d0,vlv_fr)
!
!........Coded for general case eventhough nbr_face=6 and 3D pblm
!
         DO nf_number=0,8
            istart =istart_nf(1,nf_number)
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            IF(ndim.eq.2) THEN
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  IF(npb(ii).eq.0) cycle
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     svp_nf(i0,1)=sv_nf(i0,1)*vlv_fr
                     svp_nf(i0,2)=sv_nf(i0,2)*vlv_fr
                  ENDDO
               ENDDO
            ELSE
               DO nb=1,len
                  i1=istart1+nb
                  ii=icell_nb(i1)
                  IF(npb(ii).eq.0) cycle
                  DO i=ia_nb(i1),ia_nb(i1+1)-1
                     i0=istart+i
                     svp_nf(i0,1)=sv_nf(i0,1)*vlv_fr
                     svp_nf(i0,2)=sv_nf(i0,2)*vlv_fr
                     svp_nf(i0,3)=sv_nf(i0,3)*vlv_fr
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
!
         nf_number=0
         istart =istart_nf(1,nf_number)
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         IF(ndim.eq.2) THEN
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               IF(npb(ii).gt.0) cycle
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  kk=right_non(i)
                  IF(npb(kk).eq.0) CYCLE
                  i0=istart+i
                  svp_nf(i0,1)=sv_nf(i0,1)*vlv_fr
                  svp_nf(i0,2)=sv_nf(i0,2)*vlv_fr
               ENDDO
            ENDDO
         ELSE
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
               IF(npb(ii).gt.0) cycle
               DO i=ia_nb(i1),ia_nb(i1+1)-1
                  kk=right_non(i)
                  IF(npb(kk).eq.0) CYCLE
                  i0=istart+i
                  svp_nf(i0,1)=sv_nf(i0,1)*vlv_fr
                  svp_nf(i0,2)=sv_nf(i0,2)*vlv_fr
                  svp_nf(i0,3)=sv_nf(i0,3)*vlv_fr
               ENDDO
            ENDDO
         ENDIF
!
!........Flow resistance model
!
         DO i=1,ncell_fluid
            cell%vfwl(i)=0.d0
            cell%vfwg(i)=0.d0
         ENDDO
         k_sp_high=3.6d4
         k_dp_high=1.19d5
         k_dp_low=8.6d5
         k_fd=6.d3
!
         flow_sp_sin=dsind(20.d0)*flow_sp
         flow_sp_cos=dcosd(20.d0)*flow_sp
         IF (flow_sp_cos.gt.0.9d0*flow_fd) THEN
            k_dp_i=k_dp_high
         ELSEIF (flow_sp_cos.gt.0.1d0*flow_fd) THEN
            tan_theta=flow_sp_sin/(flow_fd-flow_sp_cos)
            theta=DATAnd(tan_theta)
            fr=theta/90.d0
            k_dp_i=fr*k_dp_high+(1.d0-fr)*k_dp_low
         ELSE
            k_dp_i=k_dp_low
         ENDIF
!         
         IF(initial_kdp)THEN
            k_dp=k_dp_i
            initial_kdp=.FALSE.
         ENDIF   
!         
         k_dp=0.9d0*k_dp+0.1d0*k_dp_i
         k_sp=k_sp_high
         DO i=1,num_sp
            cell%vfwl(n_sp(i))=k_sp*ul_o(n_sp(i))
            cell%vfwg(n_sp(i))=k_sp*ug_o(n_sp(i))*cell%rhog(i)/cell%rhol(i)
         ENDDO
         DO i=1,num_fd
            cell%vfwl(n_fd(i))=k_fd*ul_o(n_fd(i))
            cell%vfwg(n_fd(i))=k_fd*ug_o(n_fd(i))*cell%rhog(i)/cell%rhol(i)
         ENDDO
         DO i=1,num_dc
            cell%vfwl(n_dc(i))=k_dp*ul_o(n_dc(i))
            cell%vfwg(n_dc(i))=k_dp*ug_o(n_dc(i))*1.d-3
         ENDDO
!
      ENDIF  
!
!
!.....T_blowdown
!
      IF(vv_prob.eq.'T_blowdown') THEN
         DO i=1,ncell_fluid
            cell%vfwl(i)=1.d5*ul_o(i)
            cell%vfwg(i)=1.d5*ug_o(i)*1.d-3
         ENDDO
      ENDIF  
!
!.....OPR1000 Full vessel      
      IF(vv_prob.eq.'OPR1000_fullvessel_1x1'            .or. &
         vv_prob.eq.'opr1000_mc_rv'                     .or. &
         vv_prob.eq.'opr1000_rv'                        .or. & 
         vv_prob.eq.'opr1000_rv_lbloca'                        .or. &         
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv'     .or. &
         vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel'    ) THEN         
         IF(initial)THEN
            initial=.FALSE.
            ALLOCATE(assem_ix(ncell_fluid),assem_iy(ncell_fluid))
            assem_ix(:)=16
            assem_iy(:)=16
            wid_fa=0.208736
            y0=wid_fa*7.d0+wid_fa*0.5d0
            DO i=1,ncell_fluid
                IF(xloc(i,3).gt.0.d0.and.xloc(i,3).lt.4.53d0) THEN
                  DO i1=1,15
                     widx0=-y0+wid_fa*dble(i1-1)
                     widx1=-y0+wid_fa*dble(i1  )
                     widy0= y0-wid_fa*dble(i1-1)
                     widy1= y0-wid_fa*dble(i1  )
                     IF(xloc(i,1).ge.widx0 .and. xloc(i,1).le.widx1) THEN
                        assem_ix(i)=i1
                     ENDIF
                     IF(xloc(i,2).le.widy0 .and. xloc(i,2).ge.widy1) THEN
                        assem_iy(i)=i1
                     ENDIF
                  ENDDO
                ENDIF   
             ENDDO
             vfwl_k(:)=vfwl_k(:)
         ENDIF
         DO i=1,ncell_fluid
            !Upper head            
            IF(nzone(i).eq.1) THEN
               vfwl_i=vfwl_k(1)*cell%rhol(i)*ul_o(i) 
            !Upper guide structure assembly
            ELSEIF(nzone(i).eq.2) THEN
               vfwl_i=vfwl_k(2)*cell%rhol(i)*ul_o(i) 
           !Upper plenum connected hot leg                
            ELSEIF(nzone(i).eq.3) THEN
               vfwl_i=vfwl_k(3)*cell%rhol(i)*ul_o(i)
            !Hot leg    
            ELSEIF(nzone(i).eq.4) THEN
               vfwl_i=vfwl_k(4)*cell%rhol(i)*ul_o(i) 
            !Upper plenum below hot leg
            ELSEIF(nzone(i).eq.5) THEN
               vfwl_i=vfwl_k(5)*cell%rhol(i)*ul_o(i)
            !Core                  
            ELSEIF(nzone(i).eq.6) THEN
               vfwl_i=vfwl_k(6)*cell%rhol(i)*ul_o(i)
            !Reflector(core)               
            ELSEIF(nzone(i).eq.7) THEN
               vfwl_i=vfwl_k(7)*cell%rhol(i)*ul_o(i)
            !!core inlet               
            ELSEIF(nzone(i).eq.8) THEN
               vfwl_i=vfwl_k(8)/fsar(assem_ix(i),assem_iy(i))**vfwl_k(9)*cell%rhol(i)*ul_o(i)
            !Reflector(Core inlet)
            ELSEIF(nzone(i).eq.9) THEN
               vfwl_i=vfwl_k(10)*cell%rhol(i)*ul_o(i)    
            !Lower support plate              
            ELSEIF(nzone(i).eq.10) THEN
               vfwl_i=vfwl_k(11)*cell%rhol(i)*ul_o(i)  
            !Flow skirt
            ELSEIF(nzone(i).eq.11) THEN
               vfwl_i=vfwl_k(12)*cell%rhol(i)*ul_o(i)  
            !Lower head
            ELSEIF(nzone(i).eq.12) THEN
               vfwl_i=vfwl_k(13)*cell%rhol(i)*ul_o(i)  
            !Downcomer
            ELSEIF(nzone(i).eq.13) THEN
               vfwl_i=vfwl_k(14)*cell%rhol(i)*ul_o(i)  
               vfwl_i=0.d0
            !Downcomer connected to cold leg
            ELSEIF(nzone(i).eq.14) THEN
               vfwl_i=vfwl_k(15)*cell%rhol(i)*ul_o(i)  
            !Cold leg
            ELSEIF(nzone(i).eq.15) THEN
               vfwl_i=vfwl_k(16)*cell%rhol(i)*ul_o(i)  
            !Vessel annulus above cold leg
            ELSEIF(nzone(i).eq.16) THEN
               vfwl_i=vfwl_k(17)*cell%rhol(i)*ul_o(i) 
            ENDIF
            !relaxation
            cell%vfwl(i)=relax_cd*cell%vfwl(i)+(1.d0-relax_cd)*vfwl_i
         ENDDO
         IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel') then
            IF(np.gt.1) CALL communicate_1d(cell%vfwl)
         ENDIF   
      ENDIF
!
!.....APR1400 Full vessel      
      IF(vv_prob.eq.'apr1400_mc_rv' .or. &
         vv_prob.eq.'apr1400_rv') THEN         
         IF(initial)THEN
            initial=.FALSE.
            ALLOCATE(assem_ix(ncell_fluid),assem_iy(ncell_fluid))
            assem_ix(:)=18
            assem_iy(:)=18
            wid_fa=0.208736
            y0=wid_fa*8.d0+wid_fa*0.5d0
            DO i=1,ncell_fluid
                IF(xloc(i,3).gt.0.0d0.and.xloc(i,3).lt.4.53d0) THEN
                  DO i1=1,17
                     widx0=-y0+wid_fa*dble(i1-1)
                     widx1=-y0+wid_fa*dble(i1  )
                     widy0= y0-wid_fa*dble(i1-1)
                     widy1= y0-wid_fa*dble(i1  )
                  IF(xloc(i,1).ge.widx0 .and. xloc(i,1).le.widx1)then
                     assem_ix(i)=i1
                  ENDIF
                  IF(xloc(i,2).le.widy0 .and. xloc(i,2).ge.widy1)then
                     assem_iy(i)=i1
                  ENDIF
                  ENDDO
                ENDIF   
             ENDDO
             vfwl_k(:)=vfwl_k(:)
         ENDIF
         DO i=1,ncell_fluid
            !Upper head            
            IF(nzone(i).eq.1) THEN
               vfwl_i=vfwl_k(1)*cell%rhol(i)*ul_o(i) 
            !Upper guide structure assembly
            ELSEIF(nzone(i).eq.2) THEN
               vfwl_i=vfwl_k(2)*cell%rhol(i)*ul_o(i) 
           !Upper plenum connected hot leg                
            ELSEIF(nzone(i).eq.3) THEN
               vfwl_i=vfwl_k(3)*cell%rhol(i)*ul_o(i)
            !Hot leg    
            ELSEIF(nzone(i).eq.4) THEN
               vfwl_i=vfwl_k(4)*cell%rhol(i)*ul_o(i) 
            !Upper plenum below hot leg
            ELSEIF(nzone(i).eq.5) THEN
               vfwl_i=vfwl_k(5)*cell%rhol(i)*ul_o(i)
            !Core                  
            ELSEIF(nzone(i).eq.6) THEN
               vfwl_i=vfwl_k(6)*cell%rhol(i)*ul_o(i)
            !Reflector(core)               
            ELSEIF(nzone(i).eq.7) THEN
               vfwl_i=vfwl_k(7)*cell%rhol(i)*ul_o(i)
            !!core inlet               
            ELSEIF(nzone(i).eq.8) THEN
               vfwl_i=vfwl_k(8)/fsar(assem_ix(i),assem_iy(i))**vfwl_k(9)*cell%rhol(i)*ul_o(i)
            !Reflector(Core inlet)
            ELSEIF(nzone(i).eq.9) THEN
               vfwl_i=vfwl_k(10)*cell%rhol(i)*ul_o(i)    
            !Lower support plate              
            ELSEIF(nzone(i).eq.10) THEN
               vfwl_i=vfwl_k(11)*cell%rhol(i)*ul_o(i)  
            !Flow skirt
            ELSEIF(nzone(i).eq.11) THEN
               vfwl_i=vfwl_k(12)*cell%rhol(i)*ul_o(i)  
            !Lower head
            ELSEIF(nzone(i).eq.12) THEN
               vfwl_i=vfwl_k(13)*cell%rhol(i)*ul_o(i)  
            !Downcomer
            ELSEIF(nzone(i).eq.13) THEN
               vfwl_i=vfwl_k(14)*cell%rhol(i)*ul_o(i)  
               vfwl_i=0.0d0
            !Downcomer connected to cold leg
            ELSEIF(nzone(i).eq.14) THEN
               vfwl_i=vfwl_k(15)*cell%rhol(i)*ul_o(i)  
            !Cold leg
            ELSEIF(nzone(i).eq.15) THEN
               vfwl_i=vfwl_k(16)*cell%rhol(i)*ul_o(i)  
            !Vessel annulus above cold leg
            ELSEIF(nzone(i).eq.16) THEN
               vfwl_i=vfwl_k(17)*cell%rhol(i)*ul_o(i) 
            ENDIF
            !relaxation
            cell%vfwl(i)=relax_cd*cell%vfwl(i)+(1.d0-relax_cd)*vfwl_i
         ENDDO
      ENDIF
!
!
      END SUBROUTINE udfn_mom_wall
!
!------------------------------------------------------------------------------------------------------------
!
      SUBROUTINE mom_wall_kakac(Re,fric)
!
!     Calcute a coefficient for Wall friction coefficitions. Only meaningful when dfs_wall_fric is used.
!   
      IMPLICIT NONE
!      
!.....Input
      REAL(8) Re
!.....Output
      REAL(8) fric
!.....Local variables
      REAL(8) fric_a,fric_b,fric_m
!      
      IF(Re.le.2100.d0)THEN
         fric_a=0.d0
         fric_b=16.d0
         fric_m=1.d0
      ELSEIF(Re.le.4000.d0)THEN
         fric_a=0.0054d0
         fric_b=2.3d-8
         fric_m=-2.d0/3.d0
      ELSE
         fric_a=1.28d-3
         fric_b=0.1143d0
         fric_m=3.2154d0
      ENDIF
!      
      fric=fric_a+fric_b/((Re)**(1.d0/fric_m))  
!      
      END SUBROUTINE mom_wall_kakac

!      
!------------------------------------------------------------------------------------------------------------
!
      SUBROUTINE mom_wall_matra(Re,fric)
!
!     Calcute a coefficient for Wall friction coefficitions. Only meaningful when dfs_wall_fric is used.
!   
      IMPLICIT NONE
!      
!.....Input
      REAL(8) Re
!.....Output
      REAL(8) fric
!.....Local variables
      REAL(8) fric_a,fric_b,fric_m
!      
      IF(Re.le.2300.d0)THEN
         fric_a=64.d0
         fric_b=-1.d0
         fric_m=0.d0
      ELSEIF(Re.gt.2300.d0.and.Re.le.30000.d0)THEN
         fric_a=0.316d0
         fric_b=-0.25d0
         fric_m=0.d0
      ELSE
         fric_a=0.184d0
         fric_b=-0.20d0
         fric_m=0.d0
      ENDIF
      
!      fric_a=0.184d0
!      fric_b=-0.20d0
!      fric_m=0.d0
       fric=fric_a*(Re**fric_b)+fric_m      
!      
      END SUBROUTINE mom_wall_matra
!
!------------------------------------------------------------------------------------------------------------------
!
      SUBROUTINE mom_wall_ctf(Re,fric)
!
!     Calcute a coefficient for Wall friction coefficitions. Only
!     meaningful when dfs_wall_fric is used.
!   
      IMPLICIT NONE
!      
!.....Input
      REAL(8) Re
!.....Output
      REAL(8) fric
!.....Local variables
      REAL(8) fric_a,fric_b,fric_m
!      
      IF(Re.le.2300.d0)THEN
         fric_a=64.d0
         fric_b=-1.d0
         fric_m=0.d0
      ELSE
         fric_a=0.204d0
         fric_b=-0.2d0
         fric_m=0.d0
      ENDIF

       fric=fric_a*((Re)**(fric_b))+fric_m
!      
      END SUBROUTINE mom_wall_ctf
