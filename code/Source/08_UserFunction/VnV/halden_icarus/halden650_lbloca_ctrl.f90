!
      SUBROUTINE halden650_lbloca_ctrl(iflag) !pik-halden
!
      USE VOL_DATA        , ONLY: cell         
      USE Zparam          , ONLY: ndim
      USE Zbc_index       , ONLY: nvin
      USE Zconst2         , ONLY: dt,dtr
      USE Zb_condition    , ONLY: pbnd,vin_liq,vin_gas
      USE Ztimecon        , ONLY: time,smac  
      USE Zimplicit       , ONLY: imp_mom_conv,imp_mom_diff,imp_alpha,imp_scalar_conv,&
                                  imp_scalar_diff,imp_ke_diff,iter_scalar,skip_imp_scalar,iter_mom
      USE Zcore           , ONLY: myrank,np      
      USE Ztplot          , ONLY: tplot_cell_loc
      USE Zio_unit        , ONLY: unit_log
      USE unitManager     , ONLY: createUnit      
      USE Zb_condition    , ONLY: cb_pl,cb_pg,alphab_liq,alphab_gas,rhob_liq,rhob_gas, &
                                   eb_liq,eb_gas,tb_liq,tb_gas,vb_liq
      USE Zpress_coeff    , ONLY: coefp_l,coefp_g                                
      USE Znum_cell       , ONLY: istart_nf,istart_nbcon_nf
      USE Zvec_index      , ONLY: left_nf,nbcon_nf
      USE Zbc_index       , ONLY: vin_norm  
      USE Zzone           , ONLY: ncell_fluid
!
      IMPLICIT NONE
!      
      INTEGER :: i,k
      INTEGER :: err,iflag,runit    
      INTEGER :: nf_number,istart,istart2,isize,i1,i2,ii
!  
      INTEGER,SAVE :: iprn        
      LOGICAL,SAVE :: INITIAL=.TRUE.      
      LOGICAL,SAVE :: initial_halden650=.TRUE.      
!    
      REAL(8),SAVE:: pres_break,pres_ambient,vel_break,dpres_break      
      REAL(8),SAVE:: tbreak_s,tbreak_f,tbreak_m,dt_min,pbnd_s,pbnd_f,pbnd_m 
!     REAL(8),SAVE:: topenleg_s,topenleg_f      
      REAL(8),SAVE:: stime_spray,ftime_spray,dtime_spray,ptime_spray 
      REAL(8),SAVE,ALLOCATABLE :: vin_liq_init(:)      
!
      IF(smac.eq.3)THEN
         imp_mom_conv=1
         imp_mom_diff=1
         imp_alpha=1
         imp_scalar_conv=1
         imp_scalar_diff=1
         imp_ke_diff=1
         iter_scalar=2
         skip_imp_scalar=1
         iter_mom=1
      ELSE
         imp_mom_conv=0
         imp_mom_diff=0
         imp_alpha=0
         imp_scalar_conv=0
         imp_scalar_diff=0
         imp_ke_diff=0
         iter_scalar=0
         skip_imp_scalar=0
         iter_mom=0
      ENDIF
!
     IF(initial)THEN
!
         initial=.FALSE.
         stime_spray=238.0d0
         ftime_spray=518.15d0
         ptime_spray=20.0d0
         dtime_spray=0.5d0
         iprn=0
!........set break time, SIT mass flow rate
         runit=createUnit('somaAddition.in')
         runit=813
         OPEN(runit,file='somaAddition.in',status='old',iostat=err)
         IF(err.ne.0.and.myrank.eq.0)THEN
            WRITE(*,"(11x,a)")'somaAddition.in is missing!'
            WRITE(unit_log,"(11x,a)")'somaAddition.in is missing!'
            PAUSE
            STOP
         ELSE
            READ(runit,*)tbreak_s,pres_ambient
            READ(runit,*)vel_break
            CLOSE(runit)   
         ENDIF
         IF(np.gt.1) CALL broadcast_r1(tbreak_s)
         tbreak_m=tbreak_s+0.1d0
         tbreak_f=tbreak_m+0.1d0
!
!........set pressure at break opening and ambient hole            
         pbnd_s=pbnd(1)
         pbnd_m=pbnd_s            
         pbnd_f=pres_ambient
         IF(myrank.eq.0)WRITE(*,"(11x,a,1f10.2,a)")'Ambient pressure is ',pres_ambient/1.d6,' MPa.'
         IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1f10.2,a)")'Ambient pressure is ',pres_ambient/1.d6,' MPa.'
!
!........save the cells and faces of 1 break, 1 outlet           
         ALLOCATE(vin_liq_init(1:nvin))
         vin_liq_init(1:nvin)=vin_liq(1:nvin)
         vin_liq(2)=0.0d0 !initially no spray
         vb_liq(2,:)=0.0d0         
!         
     ENDIF !initial 
!     
!.....closing 1 outlet ann 1 inlet, opening 1 inlet as 1 break           
!     
      IF(time.gt.tbreak_s.and.time.lt.tbreak_m)THEN 
          vin_liq(1)=-vin_liq_init(1)/(tbreak_m-tbreak_s)*(time-tbreak_m)
          IF(myrank.eq.0.and.iprn.eq.0)WRITE(*,"(11x,a,1e12.5,1i3)")'--1vin_liq(1)=',vin_liq(1)
          iprn=iprn+1
          IF(iprn.eq.10)iprn=0
      ELSEIF(time.ge.tbreak_m.and.time.lt.tbreak_f)THEN
      ELSEIF(time.ge.tbreak_f)THEN
      ENDIF
!         
!.....LBLOCA process      
!
      IF(time.ge.tbreak_m)THEN
!........blow down: close hot 
         pres_break=0.0d0
         IF(tplot_cell_loc(1).gt.0.AND.tplot_cell_loc(1).le.ncell_fluid)THEN
            pres_break=cell%p(tplot_cell_loc(1))
         ELSE
            pres_break=cell%p(1)
         ENDIF 
         CALL allreducei_max_r1(pres_break)
!         
         IF(initial_halden650)THEN
            initial_halden650=.FALSE.
            CALL halden650_lbloca_nbcon(-1)
            dpres_break=pres_break-pres_ambient
            vin_liq_init(1)=vel_break*(pres_break-pres_ambient)/dpres_break
         ENDIF  
!              
         IF(iflag.eq.1)THEN
!            
!...........blow down at inlet1            
            IF(pres_break.le.pres_ambient)THEN !inject
               vin_liq_init(1)=-vel_break*(pres_break-pres_ambient)/dpres_break  
               vin_liq_init(1)=0.0d0 
               vin_liq(1)=vin_liq_init(1)
               vin_gas(1)=vin_liq_init(1) 
            ELSE                               !break
               vin_liq_init(1)=vel_break*(pres_break-pres_ambient)/dpres_break  
               vin_liq(1)=vin_liq_init(1)
               vin_gas(1)=vin_liq_init(1)        
            ENDIF  
!            
!...........spray at inlet 2            
            IF(time.ge.stime_spray.and.time.lt.ftime_spray)THEN !spray
               IF(time.lt.stime_spray+dtime_spray)THEN
                  vin_liq(2)=vin_liq_init(2) !8.034d-4 !0.005/6.223d-3/1000.d0 !rho*vel*A=mass flow
               ELSEIF(time.lt.stime_spray+ptime_spray)THEN
                  vin_liq(2)=0.0d0
               ELSE
                  stime_spray=stime_spray+ptime_spray
                  vin_liq(2)=0.0d0
               ENDIF
            ELSE
               vin_liq(2)=0.0d0
            ENDIF    
!
!...........control inlet1,2 parameter            
            nf_number=2
            istart=istart_nf(1,nf_number)
            istart2=istart_nbcon_nf(nf_number)
            isize =istart_nf(2,nf_number)
            DO i=1,isize
               i1=istart+i
               i2=istart2+i
               ii=left_nf(i1)
               k=nbcon_nf(i2)
               IF(k.ne.1)CYCLE        
               IF(time.le.tbreak_f)THEN
                  IF(vin_norm(k).eq.0)THEN
                     vb_liq(k,1)=-vin_liq_init(k)*(time-tbreak_m)/(tbreak_f-tbreak_m)
                     vb_liq(k,2:ndim)=0.0d0
                  ELSE
                     vin_liq(k)=vin_liq_init(k)*(time-tbreak_m)/(tbreak_f-tbreak_m)
                     vin_gas(k)=vin_liq(k)
                  ENDIF 
               ENDIF
               cb_pl(k)=coefp_l(ii)
               cb_pg(k)=coefp_g(ii)
               alphab_liq(k)=cell%alphal(ii)
               alphab_gas(k)=cell%alphag(ii)
               rhob_liq(k)=cell%rhol(ii)
               rhob_gas(k)=cell%rhog(ii)
               eb_liq(k)=cell%el(ii)
               eb_gas(k)=cell%eg(ii)
               tb_liq(k)=cell%tl(ii)
               tb_gas(k)=cell%tg(ii)                  
            ENDDO   
            IF(myrank.eq.0.and.iprn.eq.0)THEN
               WRITE(*,"(11x,a,2e11.3,2f7.3)")'--vin1,2,p_brk,p_amb=',vin_liq(1),vin_liq(2),pres_break/1.d6,pres_ambient/1.d6
            ENDIF               
!            
         ELSE !iflag.eq.1
            IF(myrank.eq.0)WRITE(*,"(11x,a,1i3)")'Finish user_def_input due to iflag=',iflag
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a,1i3)")'Finish user_def_input due to iflag=',iflag
            RETURN 
         ENDIF !iflag.eq.1  
         iprn=iprn+1
         IF(iprn.eq.1000)iprn=0             
      ELSE  !time.ge.tbreak_m
  
      ENDIF !time.ge.tbreak_m  
!
!.....print spary and control print
!      
      IF(myrank.eq.0)THEN
         IF(vin_liq(2).ne.0.0d0)WRITE(*,"(a,4e12.3)")'SPRAY=',time,stime_spray+dtime_spray,stime_spray+ptime_spray,vin_liq(2)
      ENDIF   
!      
!.....reduce dt at initial status, early break stage, and early reflood stage         
!      
      IF(time.gt.tbreak_m-0.01d0.and.time.lt.tbreak_f+0.1d0)THEN
        dt_min=1.d-5
        time=time-dt
        dt=MIN(dt,dt_min)
        time=time+dt       
        dtr=1.d0/dt
      ENDIF   
!     IF(time.gt.topenleg_s.and.time.lt.topenleg_f)THEN
!       dt_min=1.d-5
!       time=time-dt
!       dt=MIN(dt,dt_min)
!       time=time+dt       
!       dtr=1.d0/dt
!     ENDIF
!      
      RETURN
!      
      ENDSUBROUTINE halden650_lbloca_ctrl
         
