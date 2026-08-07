!
      SUBROUTINE apr1400_lbloca_ctrl(iflag)
!
      USE VOL_DATA        , ONLY: cell         
      USE Znum_cell       , ONLY: i_neigh
      USE Zparam          , ONLY: pi,mesh_openfoam
      USE Zrv_choke       , ONLY: p_avg
      USE Zbc_index       , ONLY: nvin,npin,nbcon
      USE Zconst2         , ONLY: gfactor,dt,dtr
      USE Zmpi            , ONLY: jperm     
      USE Zb_condition    , ONLY: pbnd,vin_liq
      USE Zgrad_ls_c3d    , ONLY: lsindex
      USE Zrv_hts_2d      , ONLY: qvol_norm_2d,nqvol,qvol_time,l_ht_str_2d_qcell
      USE Zconst1         , ONLY: restart 
      USE Ztimecon        , ONLY: time
      USE Zwall_HTC       , ONLY: reflood,inline_bundle
      USE Zcoord1         , ONLY: xloc
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank,np      
      USE Zapr1400_lbloca , ONLY: pres_break,pres_ambient,m_max,rpv_status,&
                                  tbreak_s,tbreak_f,tbreak_m,dt_min,pbnd_s,pbnd_f,pbnd_m,time_restart,vin_liq_init,&
                                  nsit,sit_mrate,sit_mrate_time,topenleg_s,topenleg_f
      USE Zapr1400_lbloca , ONLY: tbreak_s,pres_ambient,tl_si,tg_si,& 
                                   hpsip_pre,hpsip_delay,hpsip_avail,&
                                   sit_pre,sit_mass,sit_mass_spipe,mflux_sit_phase1,mflux_sit_phase2,sit_avail
      USE Zuserdefined    , ONLY: udfl_wallHTC_porous,udfl_porous_property,udfl_mat_prop,&
                                  user_iary
      USE Ztplot          , ONLY: tplot_cell_loc
      USE Zio_unit        , ONLY: unit_log
      USE unitManager     , ONLY: createUnit      
      USE Zrv_model       , ONLY: rv_choke
!
      IMPLICIT NONE
!
      INTEGER :: i,j,k,j0  
      INTEGER :: err,iflag,runit    
      INTEGER,SAVE :: iprn          
      LOGICAL :: INITIAL=.TRUE.      
      LOGICAL :: initial_apr1400=.TRUE.      
      LOGICAL :: initial_reflood1=.true.,initial_reflood2=.true.
      REAL(8) :: time_decay
!      
      udfl_wallHTC_porous=.true.   
      udfl_porous_property=.true.  
      udfl_mat_prop=.true.         
!      
     IF(initial)THEN
!            
         initial=.FALSE.
         iprn=0
!
!........reflood period and status         
         rpv_status=0
         topenleg_s=1.0d6
         topenleg_f=1.0d6
!         
!........turn least square method off and set gfactor to 0  at 1D cells         
         IF(mesh_openfoam.ne.1.and.ncell_fluid_all.lt.30000)THEN
           DO i=1,ncell_fluid
              IF(DABS(xloc(i,1)).gt.2.5d0.or.DABS(xloc(i,2)).gt.2.5d0)THEN
                 j0=i_neigh(i)-1
                 k=0
                 DO j=i_neigh(i),i_neigh(i+1)-1
                    IF(nbcon(j).lt.0) k=k+1
                 ENDDO
                 IF(k.gt.2)THEN
                    lsindex(i)=0
                    gfactor(i)=0.0d0 
                    WRITE(*,*)'lsindex=0 at i=',jperm(i)
                 ENDIF  
              ENDIF   
           ENDDO                 
         ENDIF
         IF(mesh_openfoam.ne.1.and.ncell_fluid_all.gt.30000)THEN
           DO i=1,ncell_fluid
              IF(DABS(xloc(i,1)).gt.2.3d0.or.DABS(xloc(i,2)).gt.2.3d0)THEN
                 j0=i_neigh(i)-1
                 k=0
                 DO j=i_neigh(i),i_neigh(i+1)-1
                    IF(nbcon(j).lt.0)k=k+1
                 ENDDO
                 IF(k.gt.2)THEN
                    lsindex(i)=0
                    gfactor(i)=0.0d0
                    WRITE(*,*)'lsindex=0 at i=',jperm(i)                        
                 ENDIF   
              ENDIF   
           ENDDO                 
         ENDIF
!
!........set break time, SIT mass flow rate
         tbreak_s=100.d0
         runit=createUnit('somaAddition.in')
         runit=813
         OPEN(runit,file='somaAddition.in',status='old',iostat=err)
         IF(err.ne.0.and.myrank.eq.0)THEN
            WRITE(*,"(11x,a)")'somaAddition.in is missing!'
            WRITE(unit_log,"(11x,a)")'somaAddition.in is missing!'
            PAUSE
            STOP
         ELSE
            hpsip_avail(:)=0
            sit_avail(:)=0
            READ(runit,*)tbreak_s
            READ(runit,*)tbreak_s
            READ(runit,*)tbreak_s,pres_ambient
            READ(runit,*)tl_si,tg_si 
            READ(runit,*)hpsip_pre,hpsip_delay
            READ(runit,*)hpsip_avail(1),hpsip_avail(2),hpsip_avail(3),hpsip_avail(4)
            READ(runit,*)sit_pre,sit_mass,sit_mass_spipe,mflux_sit_phase1,mflux_sit_phase2
            READ(runit,*)sit_avail(1),sit_avail(2),sit_avail(3),sit_avail(4)
            READ(runit,*)nsit
            ALLOCATE(sit_mrate(nsit),sit_mrate_time(nsit))
            DO i=1,nsit
               READ(runit,*)sit_mrate_time(i),sit_mrate(i)
            ENDDO
            CLOSE(runit)   
         ENDIF
         IF(np.gt.1) CALL broadcast_r1(tbreak_s)
         tbreak_m=tbreak_s+0.1d0
         tbreak_f=tbreak_m+0.1d0
!         
!........change power time table according to the break time 
         IF(.not.l_ht_str_2d_qcell)THEN            
            DO i=1,nqvol 
               IF(qvol_norm_2d(1,i).lt.0.1)EXIT
            ENDDO
            time_decay=qvol_time(i)
            k=i
            DO i=1,nqvol
               IF(i.ge.k)qvol_time(i)=qvol_time(i)-time_decay+tbreak_f !apr1400_lbloca_sensitivity
               IF(myrank.eq.0.and.(i.le.5.or.i.eq.nqvol))THEN
                  WRITE(*,"(11x,a,1i3,1f10.3,1f10.3)")'iqvol,time,qvolf=',i,qvol_time(i),qvol_norm_2d(1,i)
                  WRITE(unit_log,"(11x,a,1i3,1f10.3,1f10.3)")'iqvol,time,qvolf=',i,qvol_time(i),qvol_norm_2d(1,i)
               ENDIF                  
            ENDDO 
         ENDIF  
         IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.2,a)")'Full close time is ',tbreak_m,' s.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.2,a)")'Full close time is ',tbreak_m,' s.'
         IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.2,a)")'Full break time is ',tbreak_f,' s.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.2,a)")'Full bBreak time is ',tbreak_f,' s.'
!
!........set pressure at break opening and ambient hole            
         pbnd_s=pbnd(1)
         pbnd_m=pbnd_s            
         pbnd_f=pres_ambient
         IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.2,a)")'Ambient pressure is ',pres_ambient/1.d6,' MPa.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.2,a)")'Ambient pressure is ',pres_ambient/1.d6,' MPa.'
!
!........store the cells and faces of 1 break, 4 DVIs, 1 ambient hole           
         CALL apr1400_lbloca_nbcon(0)
!         
         ALLOCATE(vin_liq_init(1:nvin))
         vin_liq_init(1:nvin)=vin_liq(1:nvin)
         time_restart=time
!         
     ENDIF !initial  
!     
!.....closing 4 cold legs, opening 1 cold leg as 1 break of pressure boundary           
!     
      IF(time.gt.tbreak_s.and.time.lt.tbreak_m.and.DABS(vin_liq(1)).gt.0.01d0)THEN ! closing
          DO i=1,nvin
            vin_liq(i)=vin_liq_init(i)*(tbreak_m-time)/(tbreak_m-tbreak_s)
          ENDDO
          IF(myrank.eq.0.and.iprn.eq.0)WRITE(*,"(11x,a,1f10.2,1i3)")'--vin_liq(nvin)=',vin_liq(nvin),nvin
          iprn=iprn+1
          IF(iprn.eq.100)iprn=0
      ELSEIF(time.ge.tbreak_m.and.time.lt.tbreak_f)THEN ! breaking
         pbnd(1)=(pbnd_f-pbnd_m)/(tbreak_f-tbreak_m)*(time-tbreak_f)+pbnd_f  
         IF(myrank.eq.0.and.iprn.eq.0)WRITE(*,"(11x,a,1pe12.3,1i3)")'--pbnd(npin) at break=',pbnd(1),npin
         iprn=iprn+1
         IF(iprn.eq.100)iprn=0                       
      ELSEIF(time.ge.tbreak_f)THEN ! broken
         pbnd(1)=pbnd_f
      ENDIF
!         
!.....LBLOCA process      
!
      IF(time.ge.tbreak_m)THEN
!         
!........blow down: close hot and cold legs, define 1 break, 4 DVIs          
         IF(initial_apr1400)THEN
            initial_apr1400=.FALSE.
!           
            CALL apr1400_lbloca_nbcon(1)
            rv_choke=1 !CHOKE-pik, find throat and calculate choke
            GOTO 1
! 
            pres_break=0.0d0
            IF(tplot_cell_loc(1).gt.0)THEN
               pres_break=cell%p(tplot_cell_loc(1))
            ELSE
               pres_break=cell%p(1)
            ENDIF   
            IF(np.gt.1) CALL allreducei_max_r1(pres_break)
         ENDIF  
!         
!........reflood: close 1 break or maintain 1 break, and open hot and cold legs   
         IF(restart.ne.0)rpv_status=user_iary(3)        
         IF(rpv_status.ne.1)THEN
            !IF(pres_break.le.pbnd_f)THEN
            IF(pres_break.le.pbnd_f*1.05d0)THEN               
               IF(myrank.eq.0)WRITE(*,"(11x,a,2f10.1,a)")'pres_break,pbnd_f=',pres_break*1.0d-6,pbnd_f*1.0d-6,'MPa.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a,2f10.1,a)")'pres_break,pbnd_f=',pres_break*1.0d-6,pbnd_f*1.0d-6,'MPa.'
               rpv_status=1
            ENDIF
            IF(np.gt.1) CALL allreducei_max_i1(rpv_status)
            user_iary(3)=rpv_status
         ENDIF
         IF(initial_reflood1)THEN
            IF(rpv_status.eq.1.and.reflood.ne.1)THEN
                initial_reflood1=.FALSE.
               IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.1,a)")'Reflood and inline_bundle is on at time of',time,'s.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.1,a)")'Reflood and inline_bundle is on at time of',time,'s.'
               topenleg_s=time+0.1d0
               topenleg_f=time+0.3d0
            ENDIF   
         ENDIF   
         IF(initial_reflood2)THEN
            IF(time.gt.topenleg_s)THEN
               initial_reflood2=.FALSE.
               reflood=1
               inline_bundle=1
               CALL apr1400_lbloca_nbcon(3)
               rv_choke=0
               IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.1,a)")'A hot leg is opened at',time,'s.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.1,a)")'A hot leg is opened at',time,'s.'
            ENDIF
         ENDIF   
         
!         
!........Critical flow rate model and SI flow                     
         IF(iflag.eq.1)THEN
!           
            CALL apr1400_lbloca_user
!            
            pres_break=0.d0
            pres_break=p_avg
            IF(np.gt.1) CALL allreducei_max_r1(pres_break)   
         ELSE
            IF(myrank.eq.0)WRITE(*,"(11x,a,1i3)")'Finish user_def_input due to iflag=',iflag
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i3)")'Finish user_def_input due to iflag=',iflag
            RETURN 
         ENDIF   
      ELSE
         pres_break=0.d0
         IF(tplot_cell_loc(1).gt.0)THEN
            pres_break=cell%p(tplot_cell_loc(1))
         ELSE
            pres_break=cell%p(1)
         ENDIF 
         IF(np.gt.1) CALL allreducei_max_r1(pres_break)
      ENDIF   
!      
!.....reduce dt at initial status, early break stage, and early reflood stage         
!      
    1 CONTINUE        
      IF(time.gt.tbreak_m-0.01d0.and.time.lt.tbreak_f+0.1d0)THEN !early break
        dt_min=1.d-4
        time=time-dt
        dt=MIN(dt,dt_min)
        time=time+dt       
        dtr=1.d0/dt
      ENDIF   
      IF(time.gt.topenleg_s-0.1d0.and.time.lt.topenleg_f)THEN !early reflood
        dt_min=1.d-4
        time=time-dt
        dt=MIN(dt,dt_min)
        time=time+dt       
        dtr=1.d0/dt
      ENDIF
!      
      RETURN
!      
      ENDSUBROUTINE apr1400_lbloca_ctrl
         
