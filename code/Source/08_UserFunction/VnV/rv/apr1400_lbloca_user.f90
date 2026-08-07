!
      SUBROUTINE apr1400_lbloca_nbcon(iflag)
!
!     Define a boundary conditoin for the ROCOM problem
!
      USE Zinterface
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid
      USE Zcore           , ONLY: myrank,np  
      USE Zparam          , ONLY: nin_max,nb_max
      USE Zbc_index       , ONLY: nbcon,npb,ngrad,icell_type,iface_wall,npin
      USE Zb_condition    , ONLY: vin_liq,pbnd,t_liq_nd,t_gas_nd,alpha_gas_nd,quala_nd      
      USE Zconst1         , ONLY: restart
      USE Znum_cell       , ONLY: i_neigh
      USE Zcoord3         , ONLY: volp
      USE Zapr1400_lbloca , ONLY: icell_loca,jth_loca,vol_loca,vol_loca_tot,num_loca,         &
                                   p_loca,rhom_loca,q_loca,h_loca,cpg_loca,vg_loca,area_loca, &
                                   vin_liq_loca,m_max,&
                                   ar_liq_loca,ar_gas_loca,p_loca,pps_loca,rhog_loca,eg_loca,     &
                                   el_loca,quala_loca,quals_loca,tl_loca,alphag_loca,alphal_loca, &
                                   cvao_loca,dcva_loca,uao_loca,ra_loca,ul_o_loca,ug_o_loca,      &
                                   vly_o_loca,vgy_o_loca,vol_g_tot,vol_l_tot,                     &
                                   num_ambient,icell_ambient,jth_ambient,pbnd_f
      USE Zvec_major      , ONLY: flux_g_nf,flux_l_nf,flux_d_nf
      USE Zvec_geo        , ONLY: sa_nf
      USE Zncg            , ONLY: tao,                                                  &
                                   cvao_npin,uao_npin,dcva_npin,ra_npin,qn_nvin,qn_npin
      USE Zncg            , ONLY: n_ncg_sp,imp_ncg,ncg_species,qn_cell0,qn_nvin,qn_npin,ncg_diff
                                       
      USE Zb_condition    , ONLY: alpha_gas_nd,alpha_liq_nd,alpha_drp_nd,&
                                  e_gas_nd,e_liq_nd,e_drp_nd,            &
                                  rho_gas_nd,rho_liq_nd,rho_drp_nd,      &
                                  t_gas_nd,t_liq_nd,t_drp_nd,            &
                                  quala_nd,pbnd
      USE Zio_unit        , ONLY: unit_log
!            
      IMPLICIT NONE
!      
!.....Input
      INTEGER iflag
!.....Local variables
      INTEGER i1
      INTEGER i,j,k,j0,m,n,nloca,nvin
      INTEGER nbcon_loca(5)                      !4dvi, 1break of 1 cold legs
      INTEGER ipbc,ipbc_ambient
      INTEGER nbcon_ambient_num,nbcon_ambient(6) ! 2hot legs, 4 cold legs
      REAL(8)  tmp1,tmp2,quala_ambient(6)    
      REAL(8) qn_npin1d(8)
!          
      DATA nbcon_loca/1,7,8,9,10/                !cold leg-1, dvi-1,2,3,4
      DATA ipbc,ipbc_ambient/0,0/ !/0,1/         !break, 5 or 6 legs during refill, 0=neumann(good),1=normal(bad)
      DATA nbcon_ambient_num/6/                  !6=close 1break and open again as a cold leg, 5=maintain 1break for nbcon==1
      DATA nbcon_ambient/5,6,2,3,4,1/  
      DATA quala_ambient/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/ !1.0 is not converged.    
!
      IF(iflag.eq.0)THEN    
!       
!........locate loca, dvi       
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Control B.C. for apr1400_lbloca.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Control B.C. for apr1400_lbloca.'
         ALLOCATE(vol_loca_tot(m_max),num_loca(m_max))
         num_loca(1:m_max)=0
         vol_loca_tot(1:m_max)=0.d0
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               DO m=1,m_max
                  IF(nbcon(j).eq.nbcon_loca(m))THEN
                     num_loca(m)=num_loca(m)+1
                     IF(m.ge.2.and.m.le.m_max)THEN
                        WRITE(*,"(11x,a,1i2,a,3i4,1i8,1i3)")'SI-',m-1,'B.C. for apr1400_lbloca=',nbcon(j),npb(i),j-j0,i,myrank
                        IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i2,a,4i8)")'SI-',m-1,'B.C. for apr1400_lbloca=',nbcon(j),npb(i),j-j0,i
                     ENDIF
                  ENDIF
               ENDDO   
            ENDDO
         ENDDO 
!
!........locate ambient 
         ALLOCATE(num_ambient(nbcon_ambient_num))        
         num_ambient(:)=0
         DO n=1,nbcon_ambient_num
            DO i=1,ncell_fluid
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).eq.nbcon_ambient(n)) num_ambient(n)=num_ambient(n)+1
               ENDDO
            ENDDO      
         ENDDO     
         ALLOCATE(icell_ambient(nbcon_ambient_num,MAXVAL(num_ambient(:))),jth_ambient(nbcon_ambient_num,MAXVAL(num_ambient(:))))
         num_ambient(:)=0
         DO n=1,nbcon_ambient_num
            DO i=1,ncell_fluid
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).eq.nbcon_ambient(n))THEN
                     num_ambient(n)=num_ambient(n)+1
                     icell_ambient(n,num_ambient(n))=i
                     jth_ambient(n,num_ambient(n))=j-j0
                  ENDIF   
               ENDDO
            ENDDO 
         ENDDO
!
!........allocate memory for break, dvi
         nloca=MAXVAL(num_loca(1:m_max))
         IF(np.gt.1) CALL allreducei_max_i1(nloca)
!               
         ALLOCATE(p_loca(m_max),rhom_loca(m_max),q_loca(m_max),h_loca(m_max),cpg_loca(m_max),vg_loca(m_max),area_loca(m_max))
!
         ALLOCATE(ar_liq_loca(m_max),ar_gas_loca(m_max),pps_loca   (m_max),rhog_loca  (m_max))
         ALLOCATE(eg_loca    (m_max),el_loca    (m_max),quala_loca (m_max),quals_loca (m_max),tl_loca    (m_max))
         ALLOCATE(alphag_loca(m_max),alphal_loca(m_max),cvao_loca  (m_max),dcva_loca  (m_max),uao_loca   (m_max))
         ALLOCATE(ra_loca    (m_max))
         ALLOCATE(ul_o_loca  (m_max),ug_o_loca (m_max))   
         ALLOCATE(vly_o_loca (m_max),vgy_o_loca(m_max))
         ALLOCATE(vol_g_tot  (m_max),vol_l_tot (m_max))   
!         
         IF(nloca.gt.0)THEN
            ALLOCATE(icell_loca(nloca,m_max),jth_loca(nloca,m_max),vol_loca(nloca,m_max))
            ALLOCATE(vin_liq_loca(nloca,m_max))
         ENDIF
!               
         DO m=1,m_max
            num_loca(m)=0
            vol_loca_tot(m)=0.d0
            IF(nloca.le.0)CYCLE
            DO n=1,nloca
               icell_loca(n,m)=0
               jth_loca(n,m)=0
               vol_loca(n,m)    =0.d0
               vin_liq_loca(n,m)=0.d0
             ENDDO 
         ENDDO
!         
!........Collect informations of 1 cold leg and 4 DVI and close 4 DVIs
         CALL nbcon_change_start    
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               DO m=1,m_max
                  IF(nbcon(j).eq.nbcon_loca(m))THEN
                     num_loca(m)=num_loca(m)+1
                     icell_loca(num_loca(m),m)=i
                     jth_loca(num_loca(m),m)=j-j0
                     vol_loca(num_loca(m),m)=volp(i)
                     vol_loca_tot(m)=vol_loca_tot(m)+volp(i)
                     IF(m.eq.1)vin_liq_loca(num_loca(m),m)=vin_liq(nbcon(j)) ! vin_norm(nvin)=1
                     IF(m.ge.2.and.m.le.m_max)THEN
                        CALL get_vector_disp(j-j0,i,i1)
                        i1=abs(i1)
                        nbcon(j)=-1
                        npb(i)=0
                        ngrad(i)=1 
                        icell_type(i)=1
                        IF(iface_wall(i).le.0)iface_wall(i)=j-j0
                        flux_g_nf(i1)=0.d0
                        flux_l_nf(i1)=0.d0
                        flux_d_nf(i1)=0.d0
                     ENDIF   
                  ENDIF  
               ENDDO
            ENDDO
         ENDDO
         CALL nbcon_change_end
         IF(np.gt.1) CALL communicate_1d_int(icell_type, &
                                             npb)
!
      ELSEIF(iflag.eq.1)THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'Define break and SI.'
!      
!........close all opening(4 cold legs, 2 hot legs)    
         CALL nbcon_change_start    
         DO i=1,ncell_fluid
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).gt.0.and.nbcon(j).le.nb_max)THEN
                  CALL get_vector_disp(j-j0,i,i1)
                  i1=abs(i1)
                  nbcon(j)=-1
                  npb(i)=0
                  flux_g_nf(i1)=0.d0
                  flux_l_nf(i1)=0.d0
                  flux_d_nf(i1)=0.d0
                  icell_type(i)=1
                  ngrad(i)=1  
                  IF(iface_wall(i).le.0)iface_wall(i)=j
               ENDIF   
            ENDDO
         ENDDO
!         
         CALL apr1400_lbloca_clean_bc_user
!         
!........open 4 DVIs
         nvin=0            
         DO m=2,m_max 
            IF(num_loca(m).le.0)CYCLE
            DO n=1,num_loca(m)
                i=icell_loca(n,m)
                j0=i_neigh(i)-1
                j=jth_loca(n,m)
                nvin=nvin+1
                nbcon(j+j0)=nvin
                npb(i)=0
                icell_type(i)=2
                ngrad(i)=1 
                IF(iface_wall(i).le.0)iface_wall(i)=j
            ENDDO          
         ENDDO
!         
!........open  a cold leg as a break 
         DO m=1,1 
            IF(num_loca(m).le.0)CYCLE
            npin=1
            DO n=1,num_loca(m)
               i=icell_loca(n,m)
               j0=i_neigh(i)-1
               j=jth_loca(n,m)
               nbcon(j+j0)=nin_max+npin
               npb(i)=nbcon(j+j0)-nin_max
               icell_type(i)=3   
               ngrad(i)=2 
               CALL get_vector_disp(j,i,i1)
               flux_g_nf(i1)=0.d0
               flux_l_nf(i1)=0.d0
               flux_d_nf(i1)=0.d0
               IF(iface_wall(i).le.0)iface_wall(i)=j
            ENDDO          
            IF(restart.eq.0)THEN
               pbnd(npin)=cell%p(i)
            ELSE
               pbnd(npin)=pbnd_f
            ENDIF   
            t_liq_nd(npin)=373.15d0      ! apr1400_lbloca_debug
            t_gas_nd(npin)=393.d0        !424(0.5MPa) 393.15(0.2MPa)
            quala_nd(npin)=0.d0          ! 1.d0 or 0.99d0 killed the code.
            IF(ipbc.eq.0)THEN            !good  
               alpha_gas_nd(npin)=100.d0 
               WRITE(*,"(a)")'Neumann P.B. is used at break.'
            ELSE                         !bad  
               alpha_gas_nd(npin)=1.0d0   
               alpha_liq_nd(npin)=0.0d0
               alpha_drp_nd(npin)=0.0d0
               quala_nd(npin)=0.0d0
               WRITE(*,"(a)")'Normal P.B. is used at break.'
            ENDIF            
               n_ncg_sp=1
               imp_ncg=0
               ncg_diff=0
               ncg_species(1)=6          !air 1~8
               qn_cell0(1)=1.d0
               qn_nvin(:,1)=1.d0
               qn_npin(:,1)=1.d0
               CALL ncg_cell
               k=npin
               qn_npin1d(:)=qn_npin(k,:)
               WRITE(*,"(a,1i3,3f7.1)")'npin,pbnd,t_liq_nd,t_gas_nd=',npin,pbnd(k)*1.d-6,t_liq_nd(k),t_gas_nd(k)
               CALL convert_temp2erg(pbnd(k),t_liq_nd(k),t_gas_nd(k),quala_nd(k),e_liq_nd(k),       &
                                      e_gas_nd(k),rho_liq_nd(k),rho_gas_nd(k),tmp1,tmp2,             &
                                      tao,cvao_npin(k),uao_npin(k),dcva_npin(k),ra_npin(k),qn_npin1d)
               t_drp_nd(k)   = t_liq_nd(k)
               rho_drp_nd(k) = rho_liq_nd(k)
               e_drp_nd(k)   = e_liq_nd(k)                   
            ENDDO          
!    
         CALL nbcon_change_end    
         IF(np.gt.1) CALL communicate_1d_int(icell_type, &
                                             npb)
!
      ELSEIF(iflag.eq.2)THEN
!         
!........Calculate physical properties at a break site for blow down and refill phase.
!........In reflood phase, this can be passed due to lack of a break site.
!        
         DO m=1,m_max
            area_loca(m)=0.d0
         ENDDO
         DO m=1,m_max
            IF(num_loca(m).le.0)CYCLE
            DO n=1,num_loca(m)
                i=icell_loca(n,m)
                j0=i_neigh(i)-1
                j=jth_loca(n,m)             
                CALL get_vector_disp(j,i,i1)
                i1=abs(i1)
                area_loca(m)=area_loca(m)+sa_nf(i1)
            ENDDO 
         ENDDO
  
         
      ELSEIF(iflag.eq.3)THEN
!
!........Close break and open cold and hot legs for reflood phase.
! 
         IF(myrank.eq.0)WRITE(*,*)'***************************************'
         IF(myrank.eq.0)WRITE(*,*)'*** Do the apr1400_lbloca_nbcon(3)! ***'
         IF(myrank.eq.0)WRITE(*,*)'***************************************'
         IF(myrank.eq.0)WRITE(unit_log,*)'***************************************'
         IF(myrank.eq.0)WRITE(unit_log,*)'*** Do the apr1400_lbloca_nbcon(3)! ***'
         IF(myrank.eq.0)WRITE(unit_log,*)'***************************************'
!
         CALL nbcon_change_start
!........close break site(1 cold legs); nbcon=1   
         IF(nbcon_ambient_num.eq.6)THEN         
         npin=0
            DO i=1,ncell_fluid
               j0=i_neigh(i)-1
               DO j=i_neigh(i),i_neigh(i+1)-1
                  IF(nbcon(j).eq.nin_max+1)THEN
                     CALL get_vector_disp(j-j0,i,i1)
                     i1=abs(i1)
                     nbcon(j)=-1
                     npb(i)=0
                     flux_g_nf(i1)=0.d0
                     flux_l_nf(i1)=0.d0
                     flux_d_nf(i1)=0.d0
                     icell_type(i)=1
                     ngrad(i)=1  
                     IF(iface_wall(i).le.0)iface_wall(i)=j                     
                  ENDIF   
               ENDDO
            ENDDO
         ELSE
            npin=1 
         ENDIF
!
!........open hot leg                  
         IF(MAXVAL(num_ambient(:)).gt.0)THEN    
            WRITE(*,*)'nbcon_ambient_num,myrank=',nbcon_ambient_num,myrank
            DO n=1,nbcon_ambient_num         !5,6,2,3,4,1
               npin=npin+1                   !npin=1; break site at cold leg         
               DO m=1,num_ambient(n)
                  i=icell_ambient(n,m)
                  j0=i_neigh(i)-1
                  j=jth_ambient(n,m)
                  nbcon(j+j0)=nin_max+npin
                  npb(i)=nbcon(j+j0)-nin_max
                  icell_type(i)=3   
                  ngrad(i)=2 
                  IF(iface_wall(i).le.0)iface_wall(i)=j                  
               ENDDO 
               IF(n.le.2)THEN
                  pbnd(npin)=pbnd_f*1.d0     !5,6-apr1400_lbloca_sensitivity
               ELSE
                  pbnd(npin)=pbnd_f*1.d0     !1.2--> 1.1/1.05 2-apr1400_lbloca_sensitivity
               ENDIF   
               t_liq_nd(npin)=373.15d0 
               t_gas_nd(npin)=424.d0          !0.5MPa
               quala_nd(npin)=quala_ambient(n)    !apr1400_lbloca_debug
               IF(ipbc_ambient.eq.0)THEN                    
                  alpha_gas_nd(npin)=100.d0   !good 
                  WRITE(*,"(a)")'Neumann P.B. is used at ambient opening.'
               ELSE             
                  alpha_gas_nd(npin)=1.0d0    !bad
                  alpha_liq_nd(npin)=0.0d0
                  alpha_drp_nd(npin)=0.0d0                   
                  WRITE(*,"(a)")'Normal P.B. is used at break.'
               ENDIF            
               n_ncg_sp=1
               imp_ncg=0
               ncg_diff=0
               ncg_species(1)=6               !air 1~8
               qn_cell0(1)=1.d0
               qn_nvin(:,1)=1.d0
               qn_npin(:,1)=1.d0
               CALL ncg_cell
               k=npin
               qn_npin1d(:)=qn_npin(k,:)
               WRITE(*,"(a,1i3,3f7.1)")'npin,pbnd,t_liq_nd,t_gas_nd=',npin,pbnd(k)*1.d-6,t_liq_nd(k),t_gas_nd(k)
               CALL convert_temp2erg(pbnd(k),t_liq_nd(k),t_gas_nd(k),quala_nd(k),e_liq_nd(k),       &
                                      e_gas_nd(k),rho_liq_nd(k),rho_gas_nd(k),tmp1,tmp2,             &
                                      tao,cvao_npin(k),uao_npin(k),dcva_npin(k),ra_npin(k),qn_npin1d)
               t_drp_nd(k)   = t_liq_nd(k)
               rho_drp_nd(k) = rho_liq_nd(k)
               e_drp_nd(k)   = e_liq_nd(k)          
            ENDDO 
         ENDIF   
         CALL nbcon_change_end
         IF(np.gt.1) CALL communicate_1d_int(icell_type, &
                                             npb)
!
      ELSE
!
         IF(myrank.eq.0)WRITE(*,*)'iflag of apr1400_lbloca_nbcon should be 0,1,2,3!!!'
         PAUSE
         STOP    
      ENDIF   
!         
      END SUBROUTINE apr1400_lbloca_nbcon     
!
!--------------------------------------------------------------------------------------------------------
!
      SUBROUTINE apr1400_lbloca_user
!
!      This routine declares flow properties at pressure outlet boundaries
!
!
      USE Zconst2        , ONLY: dt
      USE Zcore          , ONLY: myrank
      USE Ztimecon       , ONLY: time
      USE Znum_cell       , ONLY: i_neigh
      USE STM_TBL_cupid  , ONLY: st_tbl,             &
                                 nt,np,ns,ns2,ndxstd
      USE Zb_condition   , ONLY: cb_pl,cb_pg,cb_pd,alphab_liq,alphab_gas,alphab_drp,qualab,rhob_liq,rhob_gas,rhob_drp, &
                                  eb_liq,eb_gas,eb_drp,tb_liq,tb_gas,tb_drp,p_fb,vin_liq,vin_gas,vin_Drp
      USE Zbc_index      , ONLY: nbcon,npb,nvin,vin_norm
      USE Zpress         , ONLY: p
      USE Zpress_coeff   , ONLY: coefp_l,coefp_g,coefp_d
      USE Zapr1400_lbloca, ONLY: icell_loca,jth_loca,num_loca,                   &
                                 area_loca,                                      &
                                 m_max,                                          &
                                 mflux_sit,mflux_sip,mflux_dvi,                  &
                                 sit_mass,sit_mass_spipe,sit_pre,sit_detect,     &
                                 hpsip_pre,hpsip_delay,hpsip_detect,hpsip_avail, &
                                 mflux_sit_int,sit_time,hpsip_time,              &
                                 break_opt,flux_break,                           &
                                 mflux_dvi_int,mflux_sit_int,                    &
                                 nsit,sit_mrate,sit_mrate_time,                  &
                                 mflux_sit_phase1,mflux_sit_phase2,              &
                                 tl_si,tg_si   
      USE Zvec_geo       , ONLY: sa_nf
!
      USE Zuserdefined    , ONLY: user_iary,user_rary
      USE Zconst1         , ONLY: restart
      USE Zncg            , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin, &
                                  qn_nvin,qn_npin      
      
      USE Zncg            , ONLY: n_ncg_sp,imp_ncg,ncg_species,qn_nvin,qn_npin,ncg_diff
      USE Zio_unit        , ONLY: unit_log
      
!                                
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,j,k,n,m,j0
      INTEGER it
!      
      INTEGER i1,i_si
      INTEGER,SAVE :: sit_myrank,hpsip_myrank,iprn,isit
!
      LOGICAL erx            
      LOGICAL :: INITIAL=.TRUE.
      LOGICAL :: sit_print_opt=.TRUE.
      LOGICAL :: hpsip_print_opt=.TRUE.
!      
      REAL(8) sit_time_elapse
      REAL(8) sa_jk
      REAL(8) mflux_dvi_ind,mflux_sit_ind,vin_si
      REAL(8) tmp1,tmp2     
      REAL(8) tmp(2)
      REAL(8) qn_cell0(8)      
!
      REAL betafs,betags,cpfs,cpgs,entfs,entgs,hsubfs,hsubgs,         &
           kapafs,kapags,psats,s(36),tsat,usubfs,usubgs,vsubfs,vsubgs 
      EQUIVALENCE(s( 1),tsat),    &
                 (s(10),psats),   &
                 (s(11),vsubfs),  &
                 (s(12),vsubgs),  &
                 (s(13),usubfs),  &
                 (s(14),usubgs),  &
                 (s(15),hsubfs),  &
                 (s(16),hsubgs),  &
                 (s(17),betafs),  &
                 (s(18),betags),  &
                 (s(19),kapafs),  &
                 (s(20),kapags),  &
                 (s(21),cpfs),    &
                 (s(22),cpgs),    &
                 (s(25),entfs),   &
                 (s(26),entgs)
!         
!........Define APR1400 SIT,SIP...........................................................................         
!
      IF(INITIAL)THEN 
         iprn=0
         break_opt=1 !indicator to start break
         IF(restart.eq.0)THEN
            vin_liq(:)=0.d0
            sit_detect=0
            hpsip_detect=0
            flux_break=0.d0
            sit_time=0.d0
            hpsip_time=0.d0  
            mflux_sit_int=0.d0    
            mflux_dvi_int=0.d0            
         ELSE
            sit_detect=user_iary(1)
            hpsip_detect=user_iary(2)
            flux_break=user_rary(1)
            sit_time=user_rary(2)
            hpsip_time=user_rary(3)
            mflux_sit_int=user_rary(4)                 
            mflux_dvi_int=user_rary(5)            
         ENDIF
      ENDIF
!        
!....time dependent SIT mass flow rate      
      sit_time_elapse=time-sit_time
      IF(initial)THEN  
         isit=1
         mflux_sit_phase1=sit_mrate(isit)
      ELSE
         IF(sit_detect.gt.0.and.(sit_time_elapse.gt.sit_mrate_time(isit)))THEN
             IF(isit.lt.nsit) isit=isit+1
             mflux_sit_phase1=sit_mrate(isit)
         ENDIF         
      ENDIF   
!
!........Henry and Fauske critical velocity model..................................................
!
!
!.....Define Break
!
      CALL apr1400_lbloca_nbcon(2)      !physical properties at a break site
      nvin=0
!
!.....Define SI...................................................................................................
!
      mflux_sit(:)=0.d0
      mflux_sip(:)=0.d0
      mflux_dvi(:)=0.d0
      mflux_sit_ind=0.d0
      mflux_dvi_ind=0.d0      
!            
      DO m=2,m_max                     !5
!      
         i_si=m-1    
         IF(num_loca(m).le.0)CYCLE     !no si
         DO n=1,num_loca(m)
            k=icell_loca(n,m)
            IF(p(k).le.sit_pre.and.sit_detect.eq.0)THEN 
               sit_detect=1
               sit_time=time
               WRITE(*,"(11x,a,1i3)")'==>SIT is actuated!',myrank
               sit_myrank=myrank
            ENDIF
            IF(p(k).le.hpsip_pre.and.hpsip_detect.eq.0)THEN
               hpsip_detect=1
               hpsip_time=time+hpsip_delay
               WRITE(*,"(11x,a,1f10.1,a,1i3)")'==>HPSIP will be actuated at',hpsip_time,'s',myrank
               hpsip_myrank=myrank
            ENDIF
         ENDDO
!            
         DO n=1,num_loca(m)
            k=icell_loca(n,m)
            j=jth_loca(n,m)
            CALL get_vector_disp(j,k,i1)
            i1=abs(i1)
            sa_jk=sa_nf(i1)                
            IF(sit_detect.eq.1.and.time.ge.sit_time)THEN
               IF(mflux_sit_int.lt.sit_mass_spipe)THEN
                  mflux_sit(i_si)=mflux_sit_phase1            !kg/s
               ELSEIF(mflux_sit_int.lt.sit_mass)THEN
                  mflux_sit(i_si)=mflux_sit_phase2            !kg/s 
               ELSE
                  mflux_sit(i_si)=0.d0 
               ENDIF   
            ELSE
               mflux_sit(i_si)=0.d0   
            ENDIF
            IF(hpsip_detect.eq.1.and.time.ge.hpsip_time)THEN
               mflux_sip(i_si)=74.d0*hpsip_avail(i_si)       !kg/s
            ELSE
               mflux_sip(i_si)=0.d0   
            ENDIF
            mflux_dvi(i_si)=mflux_sit(i_si)+mflux_sip(i_si)
!
!...........Set SI b.c.
      
            nvin=nvin+1
            cb_pl(nvin)=coefp_l(k)
            cb_pg(nvin)=coefp_g(k)
            cb_pd(nvin)=coefp_d(k)
		      IF(initial)THEN
	            npb(k)=0
               j0=i_neigh(k)-1
	            nbcon(j+j0)=nvin
	            vin_norm(nvin)=1  
		         alphab_liq(nvin)=1.d0
		         alphab_gas(nvin)=0.d0
		         alphab_drp(nvin)=0.d0
		         tb_liq(nvin)=tl_si                  !383.15 !302.15d0 !apr1400_lbloca_debug_si
		         tb_gas(nvin)=tg_si                  !383.15 !302.15d0
		         tb_drp(nvin)=tb_liq(nvin)
		         p_fb(nvin)=0.2d6                    !p(k)!SPACE-sk34
		         qualab(nvin)=0.d0
		         IF(1)THEN    
			         n_ncg_sp=1
			         imp_ncg=0
			         ncg_diff=0
			         ncg_species(1)=6            !air 1~8
			         qn_cell0(1)=1.d0
			         qn_nvin(:,1)=1.d0
			         qn_npin(:,1)=1.d0
			         CALL ncg_cell
			         i=nvin
                  qn_cell0(:)=qn_nvin(i,:)                 
			         CALL convert_temp2erg(p_fb(i),tb_liq(i),tb_gas(i),qualab(i),eb_liq(i),eb_gas(i),rhob_liq(i),rhob_gas(i),tmp1,tmp2, &
						             tao,cvao_nvin(i),uao_nvin(i),dcva_nvin(i),ra_nvin(i),qn_cell0)
			         rhob_drp(nvin)=rhob_liq(nvin)
			         eb_drp(nvin)  =eb_liq(nvin)  
		         ELSE
!            Initialize s for sth2x3_cupid
!                 s(:)=0.d0
                  s(1)=tb_liq(nvin)
                  s(2)=p_fb(nvin)
                  CALL sth2x3_cupid(s,it,erx,                          &
                                    st_tbl(ndxstd),                    &
                                    st_tbl(ndxstd+nt),                 &
                                    st_tbl(ndxstd+nt+np+13*ns+13*ns2))
                  IF(erx)then
                     print *, '#### ERROR: sth2x3_cupid called from apr1400_lbloca_user'
                     pause
                     stop
                  ENDIF
                  rhob_liq(nvin)=1.d0/vsubfs
                  eb_liq(nvin)=usubfs
                  rhob_drp(nvin)=rhob_liq(nvin)
                  eb_drp(nvin)= eb_liq(nvin)
                  s(2)=p_fb(nvin)
                  CALL sth2x2_cupid(st_tbl(ndxstd),st_tbl(ndxstd+nt),s,erx)
                  rhob_gas(nvin)=1.d0/vsubgs
                  eb_gas(nvin)=usubgs
               ENDIF  
            ENDIF     !initial                      
!                        
!...........Consider direction of SI,             
            vin_si=-mflux_dvi(i_si)/(area_loca(m)*rhob_liq(nvin))     !vin_norm(nvin)=1 
            vin_liq(nvin)=vin_liq(nvin)*0.99d0+vin_si*0.01d0          !apr1400_lbloca_debug_si
            vin_gas(nvin)=0.d0
            vin_drp(nvin)=0.d0
            IF(DABS(vin_si).gt.0.d0)THEN
               mflux_sit_ind=mflux_sit_ind+mflux_sit(i_si)*dt*DABS(vin_liq(nvin)/vin_si)
               mflux_dvi_ind=mflux_dvi_ind+mflux_dvi(i_si)*dt*DABS(vin_liq(nvin)/vin_si)
            ENDIF
            IF(iprn.eq.0)WRITE(*,"(11x,1i3,a,2f10.2,2i3)")i_si,'-SI: vel,p=',vin_liq(nvin),p_fb(nvin)/1.d6,nvin,myrank
         ENDDO
         IF(iprn.eq.0)WRITE(*,"(13x,a/13x,2i3,3e12.3)")'sit_det,hpsip_det,sit_time,hpsip_time,mflux_sit_int=',sit_detect,hpsip_detect,sit_time,hpsip_time,mflux_sit_int
      ENDDO      
      IF(np.gt.1) THEN
         tmp(1)=mflux_sit_ind
         tmp(2)=mflux_dvi_ind
         CALL allreducei_r(tmp,2)
         mflux_sit_ind=tmp(1)
         mflux_sit_ind=tmp(2)
      ENDIF
      mflux_sit_int=mflux_sit_int+mflux_sit_ind
      mflux_dvi_int=mflux_dvi_int+mflux_dvi_ind
!
!.....for the restart
!
      user_iary(1)=sit_detect              
      user_iary(2)=hpsip_detect
      user_rary(1)=flux_break   
      user_rary(2)=sit_time     
      user_rary(3)=hpsip_time   
      user_rary(4)=mflux_sit_int
      user_rary(5)=mflux_dvi_int
      !user_rary(6)=vl_choke     
      !user_rary(7)=vg_choke     
      !user_rary(8)=vl_choke_o   
      !user_rary(9)=vg_choke_o                    
!
!.....Control io etc.
!
      IF(sit_print_opt)THEN
         IF(np.gt.1) THEN
            CALL allreducei_max_i1(sit_myrank)
            CALL allreducei_max_i1(sit_detect)
            CALL allreducei_max_r1(sit_time)
         ENDIF
         IF(sit_detect.eq.1.and.time.ge.sit_time)THEN 
            sit_print_opt=.FALSE.
            IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.3,a,1i3,a)")'****** SIT is actuated at',sit_time,'s at rank',sit_myrank,' *****'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.3,a,1i3,a)")'********** SIT is actuated at',sit_time,'s at rank',sit_myrank,' **********'
         ENDIF
      ENDIF   
      IF(hpsip_print_opt)THEN
         IF(np.gt.1)THEN
            CALL allreducei_max_i1(hpsip_myrank)
            CALL allreducei_max_i1(hpsip_detect)               
            CALL allreducei_max_r1(hpsip_time)
         ENDIF
         IF(hpsip_detect.eq.1.and.time.ge.hpsip_time)THEN
            hpsip_print_opt=.FALSE.
            IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.3,a,1i3)")'***** HPSIP is actuated at',hpsip_time,'s at rank',hpsip_myrank,' *****'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.3,a,1i3)")'********** HPSIP is actuated at',hpsip_time,'s at rank',hpsip_myrank,' **********'
         ENDIF
      ENDIF 
!
      iprn=iprn+1
      IF(iprn.eq.10000)iprn=0               !apr1400_lbloca_debug_print
!
      IF(initial)initial=.FALSE.    
!
      END SUBROUTINE apr1400_lbloca_user
!
!------------------------------------------------------------------------------------------------
!
      SUBROUTINE apr1400_lbloca_clean_bc_user
!
      USE Zbc_index    , ONLY: nvin,npin,vin_norm
      USE Zb_condition , ONLY: e_gas_nd,e_liq_nd,         &
                               eb_gas,eb_liq,             &
                               rho_gas_nd,rho_liq_nd,     &
                               rhob_gas,rhob_liq,         &
                               t_gas_nd,t_liq_nd,         &
                               tb_gas,tb_liq,             &
                               qualab,quala_nd,pbnd,p_fb, &
                               vb_liq,vb_gas,vb_drp,      &
                               vin_liq,vin_gas,vin_Drp
      USE Zncg         , ONLY: tao,cvao_nvin,uao_nvin,dcva_nvin,ra_nvin,             &
                                cvao_npin,uao_npin,dcva_npin,ra_npin,qn_nvin,qn_npin
!
      IMPLICIT NONE                                
!
      tao          =0.d0 
!      
      nvin         =0
      vin_norm(:)  =0.d0
      vin_gas(:)   =0.d0
      vin_liq(:)   =0.d0
      vin_drp(:)   =0.d0
      vb_gas(:,:)  =0.d0
      vb_liq(:,:)  =0.d0
      vb_drp(:,:)  =0.d0
      p_fb(:)      =0.d0
      tb_liq(:)    =0.d0
      tb_gas(:)    =0.d0
      qualab(:)    =0.d0
      eb_liq(:)    =0.d0
      eb_gas(:)    =0.d0
      rhob_liq(:)  =0.d0
      rhob_gas(:)  =0.d0
      cvao_nvin(:) =0.d0
      uao_nvin(:)  =0.d0
      dcva_nvin(:) =0.d0
      ra_nvin(:)   =0.d0
      qn_nvin(:,:) =0.d0
!      
      npin         =0
      pbnd(:)      =0.d0
      t_liq_nd(:)  =0.d0
      t_gas_nd(:)  =0.d0
      quala_nd(:)  =0.d0
      e_liq_nd(:)  =0.d0
      e_gas_nd(:)  =0.d0
      rho_liq_nd(:)=0.d0 
      rho_gas_nd(:)=0.d0 
      cvao_npin(:) =0.d0
      uao_npin(:)  =0.d0
      dcva_npin(:) =0.d0
      ra_npin(:)   =0.d0
      qn_npin(:,:) =0.d0
!      
      END SUBROUTINE apr1400_lbloca_clean_bc_user
!------------------------------------------------------------------------------------------------------------------------
!
      SUBROUTINE apr1400_lbloca_out_user_sub
!
!     This routine writes the calculation results for 2D_loca problem
!
      USE VOL_DATA     , ONLY: cell
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst2      , ONLY: dt
      USE Zcore        , ONLY: np
      USE Zb_condition , ONLY: eb_liq,rhob_liq
      USE Zcoord3      , ONLY: vol
      USE Zqvol        , ONLY: qvol_liq
      USE Zsbloca      , ONLY: break_flow,break_flow_eng_l,break_flow_eng_g,si_flow,si_flow_eng, &
                                q_liq,break_flow_eng,break_flow_int,break_flow_eng_int,          &
                                break_flow_eng_l_int,break_flow_eng_g_int,si_flow_int,           &
                                si_flow_eng_int,q_liq_int
!
      USE Zvec_major   , ONLY: flux_l_nf,flux_g_nf,flux_d_nf
!
      USE Zapr1400_lbloca  , ONLY: icell_loca,jth_loca,num_loca,&
                                   m_max,       &
                                   break_opt
!
      USE Zconst1  , ONLY: restart
      USE Zuserdefined ,ONLY: user_rary    
!
      IMPLICIT NONE
!      
!.....Local variables
      INTEGER :: i,j,k
      INTEGER :: i1
      INTEGER :: nvin_tmp,n,m,i_si 
      LOGICAL :: initial=.true.
      REAL(8) :: delta
      REAL(8) tmp(5)
!      
      IF(initial)THEN
         IF(restart.eq.0)THEN
            break_flow_int      =0.d0
            break_flow_eng_int  =0.d0
            break_flow_eng_g_int=0.d0
            break_flow_eng_l_int=0.d0
            si_flow_int    =0.d0
            si_flow_eng_int=0.d0
            q_liq_int      =0.d0
         ELSE
            break_flow_int      = user_rary(11)
            break_flow_eng_int  = user_rary(12)
            break_flow_eng_g_int= user_rary(13)
            break_flow_eng_l_int= user_rary(14)
            si_flow_int         = user_rary(15)
            si_flow_eng_int     = user_rary(16)
         ENDIF
      ENDIF
!
      break_flow      =0.d0
      break_flow_eng_g=0.d0
      break_flow_eng_l=0.d0
      si_flow    =0.d0
      si_flow_eng=0.d0
      nvin_tmp=0
!
      IF(break_opt.eq.0)GOTO 1
!
      IF(num_loca(1).gt.0)THEN
         m=1
         DO n=1,num_loca(m)
            k=icell_loca(n,m)
            j=jth_loca(n,m)
            CALL get_vector_disp(j,k,i1)
            IF(i1.gt.0)THEN
                delta=1.d0           
            ELSE
                delta=-1.d0
                i1=-i1
            ENDIF
            break_flow      = break_flow+(cell%alphal(k)*cell%rhol(k)*flux_l_nf(i1)       &
                                        +cell%alphad(k)*cell%rhol(k)*flux_d_nf(i1)        &
                                        +cell%alphag(k)*cell%rhog(k)*flux_g_nf(i1))*delta
            break_flow_eng_g= break_flow_eng_g+(cell%alphag(k)*cell%rhog(k)*cell%eg(k)*flux_g_nf(i1))*delta
            break_flow_eng_l= break_flow_eng_l+(cell%alphal(k)*cell%rhol(k)*cell%el(k)*flux_l_nf(i1)       &
                                              +cell%alphad(k)*cell%rhol(k)*cell%el(k)*flux_d_nf(i1))*delta
         ENDDO
      ENDIF
!            
      DO m=2,m_max                 !5
!      
         i_si=m-1    
         IF(num_loca(m).le.0)CYCLE !no si
         DO n=1,num_loca(m)
            nvin_tmp=nvin_tmp+1
            k=icell_loca(n,m)
            j=jth_loca(n,m)            
            CALL get_vector_disp(j,k,i1)
            i1=abs(i1)
            si_flow    =si_flow-rhob_liq(nvin_tmp)*flux_l_nf(i1)
            si_flow_eng=si_flow_eng-rhob_liq(nvin_tmp)*eb_liq(nvin_tmp)*flux_l_nf(i1)
         ENDDO
      ENDDO    
    1 CONTINUE      
!
      q_liq=0.d0
      DO i=1,ncell_fluid
         q_liq=q_liq+qvol_liq(i)*vol(i)
      ENDDO
      IF(np.gt.1) CALL allreducei_r1(q_liq)
!
      IF(np.gt.1)THEN
         tmp(1)=break_flow
         tmp(2)=break_flow_eng_g
         tmp(3)=break_flow_eng_l
         tmp(4)=si_flow
         tmp(5)=si_flow_eng
         CALL allreducei_r(tmp,5)
         break_flow      =tmp(1)
         break_flow_eng_g=tmp(2)
         break_flow_eng_l=tmp(3)
         si_flow         =tmp(4)
         si_flow_eng     =tmp(5)
      ENDIF
!
      break_flow_eng      =break_flow_eng_g    +break_flow_eng_l
      break_flow_int      =break_flow_int      +break_flow      *dt
      break_flow_eng_int  =break_flow_eng_int  +break_flow_eng  *dt
      break_flow_eng_g_int=break_flow_eng_g_int+break_flow_eng_g*dt
      break_flow_eng_l_int=break_flow_eng_l_int+break_flow_eng_l*dt
      si_flow_int    =si_flow_int    +si_flow*dt
      si_flow_eng_int=si_flow_eng_int+si_flow_eng*dt
      q_liq_int      =q_liq_int      +q_liq      *dt
!
      IF(initial)initial=.FALSE.
!
!.....for the restart
!
      user_rary(11)=break_flow_int      
      user_rary(12)=break_flow_eng_int  
      user_rary(13)=break_flow_eng_g_int
      user_rary(14)=break_flow_eng_l_int
      user_rary(15)=si_flow_int         
      user_rary(16)=si_flow_eng_int 
!
      END SUBROUTINE apr1400_lbloca_out_user_sub
!----------------------------------------------------------------------
!
      SUBROUTINE apr1400_lbloca_out_user
!
!     Save output for TECPLOT
!
      USE VOL_DATA , ONLY: cell
      USE Zzone    , ONLY: ncell_fluid                            
      USE Zcore    , ONLY: np,myrank      
      USE Ztimecon , ONLY: time
      USE Zare     , ONLY: are_liq,are_gas,are_drp
      USE Zcoord3  , ONLY: vol
      USE Zsbloca  , ONLY: p_break,t_break,q_break,h_break,a_break,si_flow,si_flow_int,si_flow_eng_int, &
                           rhom_break,break_flow_int,break_flow_eng_int,break_flow_eng_l_int,           &
                           break_flow_eng_g_int,q_liq_int,q_liq,eng_gg,eng_gg_int,ge_err,ge_err_int,    &
                           pw_l,pw_g,pw_g_int,break_flow,break_flow_eng,si_flow_eng
      USE Zapr1400_lbloca  , ONLY:flux_break,pres_break
!
      USE Zconst1  , ONLY: restart
      USE Zuserdefined ,ONLY: user_rary         
!
      IMPLICIT NONE 
!
      CHARACTER*20 varname(35)      
!
!.....Local variables
      INTEGER i
      LOGICAL :: INITIAL=.TRUE.
!      
! bug  not saved variables
      REAL(8),SAVE :: rcs_mass0, rcs_energy0
      REAL(8) rcs_liq_vol, rcs_mass, rcs_energy, rcs_mass_s, rcs_energy_s, err_mass, err_eng
      REAL(8) rcs_energy_g, rcs_energy_l, rcs_energy_s_g, rcs_energy_s_l
      REAL(8) break_flux
      REAL(8) tmp(12)
!
      DATA varname / 'time', 'p_break', 't_break', 'q_break', 'rhom_break', 'h_break', 'a_break',            &
                      'break_flow', 'break_flow_int', 'si_flow', 'si_flow_int', 'rcs_liq_vol', 'rcs_mass',   &
                      'rcs_mass_s', 'q_liq', 'q_liq_int', 'ge_err', 'ge_err_int', 'break_flow_eng',          &
                      'break_flow_eng_int', 'si_flow_eng', 'si_flow_eng_int', 'rcs_energy', 'rcs_energy_s',  &
                      'err_mass', 'err_eng','rcs_energy_g','rcs_energy_s_g','rcs_energy_l','rcs_energy_s_l', &
                      'eng_gg','eng_gg_int','pw_g','pw_l','break_flux'/
!
!.....SBLOCA output
!
      rcs_liq_vol=0.d0
      rcs_mass=0.d0
      rcs_energy=0.d0
      rcs_energy_g=0.d0
      rcs_energy_l=0.d0
      DO i=1,ncell_fluid
         rcs_liq_vol =rcs_liq_vol +vol(i)*cell%alphal(i)
         rcs_mass    =rcs_mass    +(cell%alphal(i)*cell%rhol(i)+cell%alphag(i)*cell%rhog(i)+cell%alphad(i)*cell%rhod(i))*vol(i)
         rcs_energy_g=rcs_energy_g+vol(i)*are_gas(i)
         rcs_energy_l=rcs_energy_l+vol(i)*(are_liq(i)+are_drp(i))
         rcs_energy  =rcs_energy_g+rcs_energy_l 
      END DO
!
      break_flux=flux_break
      p_break=pres_break
!      
      IF(np.gt.1)THEN
         tmp( 1)=rcs_liq_vol
         tmp( 2)=rcs_mass
         tmp( 3)=rcs_energy_g
         tmp( 4)=rcs_energy_l
         tmp( 5)=rcs_energy
         tmp( 6)=p_break
         tmp( 7)=t_break
         tmp( 8)=q_break
         tmp( 9)=rhom_break
         tmp(10)=h_break
         tmp(11)=a_break
         tmp(12)=break_flux
         CALL allreducei_r(tmp,12)
         rcs_liq_vol =tmp( 1)
         rcs_mass    =tmp( 2)
         rcs_energy_g=tmp( 3)
         rcs_energy_l=tmp( 4)
         rcs_energy  =tmp( 5)
         p_break     =tmp( 6)
         t_break     =tmp( 7)
         q_break     =tmp( 8)
         rhom_break  =tmp( 9)
         h_break     =tmp(10)
         a_break     =tmp(11)
         break_flux  =tmp(12)
      ENDIF
!
      IF(INITIAL)THEN
         IF(restart.eq.0)THEN
            rcs_mass0    =rcs_mass
            rcs_energy0  =rcs_energy
            user_rary(17)=rcs_mass0
            user_rary(18)=rcs_energy0
         ELSE
            rcs_mass0  =user_rary(17)
            rcs_energy0=user_rary(18)         
         ENDIF
      ENDIF
      rcs_mass_s    =rcs_mass0  -break_flow_int + si_flow_int
      rcs_energy_s  =rcs_energy0-break_flow_eng_int  +si_flow_eng_int+q_liq_int+pw_g_int
      rcs_energy_s_l=rcs_energy0-break_flow_eng_l_int+si_flow_eng_int+q_liq_int-eng_gg_int
      rcs_energy_s_g=-break_flow_eng_g_int+eng_gg_int+pw_g_int
      err_mass      =(rcs_mass_s-rcs_mass)/rcs_mass
      err_eng       =(rcs_energy_s-rcs_energy)/rcs_energy
!
      IF(initial.and.myrank.eq.0)WRITE(334,5000)(varname(i),i=1,35)
5000  FORMAT(1x,35(A20))
      
      IF(myrank.eq.0)THEN
         WRITE(334,100)time, p_break, t_break, q_break, rhom_break, h_break, a_break,               &
                       break_flow, break_flow_int, si_flow, si_flow_int, rcs_liq_vol, rcs_mass,     &
                       rcs_mass_s, q_liq, q_liq_int, ge_err, ge_err_int, break_flow_eng,            &
                       break_flow_eng_int, si_flow_eng, si_flow_eng_int, rcs_energy, rcs_energy_s,  &
                       err_mass, err_eng,rcs_energy_g,rcs_energy_s_g,rcs_energy_l,rcs_energy_s_l,   &
                       eng_gg,eng_gg_int,pw_g,pw_l,break_flux
      ENDIF
100   FORMAT(35(e20.10,1x))
!
      IF(initial)initial=.FALSE.
!
      END SUBROUTINE apr1400_lbloca_out_user
