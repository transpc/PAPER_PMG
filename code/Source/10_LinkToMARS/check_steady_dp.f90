!
      SUBROUTINE check_steady_dp
!      
      USE MASTER4      , ONLY: dpres
      USE Zcore        , ONLY: myrank
      USE Zconst1      , ONLY: cplmars
      USE Ztimecon     , ONLY: time 
!      
      IMPLICIT NONE
!      
      INCLUDE 'c3com.h' !n_volleg,p_volleg
!      
      INTEGER :: idpres,k,n_volleg
      INTEGER,SAVE::run_mode3_judge,run_mode2_judge,run_mode2=0,run_mode3=0
      LOGICAL,SAVE::initial,inital_c3RPV3d,inital_c3p,initial_two
      REAL(8) :: dp_mcc(2)
      REAL(8),SAVE::time_c3RPV3d,time_c3p,time_two
      REAL(8),SAVE::dpnew,dpold,dpold_old,dpcri,time_judge,time_old,dtime,pnew,pold     
!
      !!!INTEGER,SAVE::count=0,loop=0
      !!!REAL(8),SAVE::time_avg,dpsum,dpavg,dpind(1000)
!      
      DATA initial,inital_c3RPV3d,inital_c3p,initial_two/.TRUE.,.TRUE.,.TRUE.,.TRUE./
!
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(dpres(8))
      ENDIF
!      
!......pressure drop at 1D RPV      
      dpres(1)=c3p_volleg(2)-c3p_volleg(1)
      dpres(2)=c3p_volleg(3)-c3p_volleg(1)
      dpres(3)=c3p_volleg(5)-c3p_volleg(4)
      dpres(4)=c3p_volleg(6)-c3p_volleg(4)      
!      
!.....pressure drop at 3D RPV     
      dpres(5)=c3p_volleg(8)-c3p_volleg(7)
      dpres(6)=c3p_volleg(9)-c3p_volleg(7)
      dpres(7)=c3p_volleg(11)-c3p_volleg(10)
      dpres(8)=c3p_volleg(12)-c3p_volleg(10)
!
      IF(cplmars.ne.2)RETURN
      IF(i3run_mode.lt.2)RETURN
!
!.....branch i3run_mode into 2 or 3, when it was 4 initially
      IF(i3run_mode.eq.4)THEN
         IF(initial_two)THEN
            initial_two=.FALSE.
            time_two=time
            run_mode2=0
            run_mode3=1
         ENDIF
         IF((time-time_two.gt.c3run_mode3_dur).and.run_mode3_judge.eq.0)THEN
            run_mode2=1
            run_mode3=0
            time_two=time
            inital_c3p=.TRUE.
            inital_c3RPV3d=.TRUE.
            run_mode3_judge=1
         ENDIF    
         IF((time-time_two.gt.c3run_mode2_dur).and.run_mode2_judge.eq.0)THEN
            run_mode2=0
            run_mode3=1
            time_two=time
            inital_c3p=.TRUE.
            inital_c3RPV3d=.TRUE.
            run_mode2_judge=1
         ENDIF             
      ENDIF   
      IF(i3run_mode.eq.2)THEN
         run_mode2=1
         run_mode3=0
      ELSEIF(i3run_mode.eq.3)THEN
         run_mode2=0
         run_mode3=1
      ENDIF   
!      
!.....adjust kfactor of 1D lower plenum referring pressure loss between 3D hot leg and 3D cold leg.
      IF(run_mode2.eq.1)THEN
         IF(inital_c3RPV3d)THEN
            inital_c3RPV3d=.FALSE.
            time_c3RPV3d=time
            pnew=1.d7
            pold=pnew            
            dpnew=1.d7
            dpold=dpnew
            dpold_old=dpold
            dpcri=c3run_mode2_dpcri
            time_judge=time
            run_mode2_judge=1            
            !!!time_avg=time
         ENDIF
!
         !!!IF(time.gt.time_avg)THEN
         !!!   time_avg=time_avg+0.001d0 
         !!!   count=count+1
         !!!   loop=loop+1
         !!!   IF(loop.gt.1000)loop=1
         !!!   dpsum=dpsum+dpres(5) !add
         !!!   IF(count.gt.1000)dpsum=dpsum-dpind(loop)
         !!!   dpind(loop)=dpres(5)
         !!!   dpavg=dpsum/MIN(count,1000)
         !!!ENDIF   
!         
         IF(time-time_c3RPV3d.gt.5.0d0)THEN 
         ELSEIF(time-time_c3RPV3d.gt.0.0d0)THEN 
            idpres=dpres(5)/100
            c3RPV3d_ctl_val=idpres*100
            !!!c3RPV3d_ctl_val=dpavg !no gain
         ELSE
         ENDIF
      ELSE
         c3RPV3d_ctl_val=-1.0d0
      ENDIF
!  
!.....adjust pressure at 3D hot legs referring pressure at 1D hot legs.    
      IF(run_mode3.eq.1)THEN
         IF(inital_c3p)THEN
            inital_c3p=.FALSE.
            time_c3p=time
            pnew=1.d7
            pold=pnew
            dpnew=1.d7
            dpold=dpnew
            dpold_old=dpold
            dpcri=c3run_mode3_dpcri
            time_judge=time
            run_mode3_judge=1
         ENDIF
         IF(time-time_c3p.gt.0.0d0)THEN
            n_volleg=i3n_volleg/2
            dp_mcc(1)=c3p_volleg(1)-c3p_volleg(1+n_volleg)
            dp_mcc(2)=c3p_volleg(4)-c3p_volleg(4+n_volleg)
            DO k=1,i3n_tmdpvol2nd
               c3p_tmdpvol2nd(k)=c3p_tmdpvol2nd(k)+DMIN1(0.002d+7,DMAX1(-0.002d+7,dp_mcc(k)/2.d0))
            ENDDO 
         ENDIF
         !!!IF(i3run_mode.eq.4)THEN !CUPVOL pressure is not equal to 1D
         !!!   idpres=dpres(5)/100
         !!!   c3RPV3d_ctl_val=idpres*100
         !!!ENDIF   
      ELSE
      ENDIF 
!  
!........judge pressure at 3D hot legs using differential
      IF(time.gt.time_judge)THEN
         dtime=time-time_old
         time_judge=time_judge+0.5d0
         time_old=time
         IF(run_mode3.eq.1)THEN  
            pold=pnew
            pnew=c3p_tmdpvol2nd(1)
            dpold_old=dpold
            dpold=dpnew
            dpnew=(pnew-pold)/dtime
            IF(myrank.eq.0)WRITE(*,"(a)")'****************************************************************************'
            IF(myrank.eq.0)WRITE(*,"(a,4(1pe14.5))")'MODE3: dp1,2,3,cri=',dpold_old,dpold,dpnew,dpcri
            IF(myrank.eq.0)WRITE(*,"(a,3(1pe14.5),1i3)")'t,t_ref,t_dur,judge=',time,time_two,c3run_mode3_dur,run_mode3_judge
            IF(myrank.eq.0)WRITE(*,"(a,4(1pe14.5))")'c3p_tmdpvol2nd(1,2)=',c3p_tmdpvol2nd(1),c3p_tmdpvol2nd(2)
            IF(myrank.eq.0)WRITE(*,"(a)")'****************************************************************************'
            IF(DABS(dpold_old).lt.dpcri.and.DABS(dpold).lt.dpcri.and.DABS(dpnew).lt.dpcri)THEN
               run_mode3_judge=0
            ENDIF 
         ELSEIF(run_mode2.eq.1)THEN
            pold=pnew
            pnew=c3RPV3d_ctl_val
            dpold_old=dpold
            dpold=dpnew
            dpnew=(pnew-pold)/dtime            
            IF(myrank.eq.0)WRITE(*,"(a)")'#############################################################################'
            IF(myrank.eq.0)WRITE(*,"(a,4(1pe14.5))")'MODE2: dp1,2,3,cri=',dpold_old,dpold,dpnew,dpcri
            IF(myrank.eq.0)WRITE(*,"(a,3(1pe14.5),1i3)")'t,t_ref,t_dur,judge=',time,time_two,c3run_mode2_dur,run_mode2_judge            
            IF(myrank.eq.0)WRITE(*,"(a,4(1pe14.5))")'c3RPV3d_ctl_val,dpres(5)=',c3RPV3d_ctl_val,dpres(5)
            IF(myrank.eq.0)WRITE(*,"(a)")'#############################################################################'               
            IF(DABS(dpold_old).lt.dpcri.and.DABS(dpold).lt.dpcri.and.DABS(dpnew).lt.dpcri)THEN
               run_mode2_judge=0
            ENDIF 
         ENDIF      
      ENDIF   
!      
      RETURN               
      ENDSUBROUTINE check_steady_dp 
!---------------------------------------------------------------------------
      SUBROUTINE print_ss_status_dp
!
      USE Zcore      ,ONLY:myrank
      USE Zmars      ,only: time_mars
      USE Zconst1    ,ONLY:cplmars
      USE MASTER4    ,ONLY:dpres,dmaster_pass
!      
      IMPLICIT NONE
!      
      INCLUDE 'c3com.h' 
!
      CHARACTER*14 varname(10),varname2(12)      
!
      INTEGER i,j
!
      LOGICAL,SAVE::initial
!
      REAL(8),SAVE::count,c3kfactor_ctl_val_sum
!      
      DATA initial/.true./ 
      DATA varname/'Time(s)','K_factor',&
                   'dpMARS_1a','dpMARS_1b','dpMARS_2a','dpMARS_2b', &
                   'dpCUP_1a','dpCUP_1b','dpCUP_2a','dpCUP_2b'/
      DATA varname2/ 'MF_HL1', 'MF_CL1a', 'MF_CL1b','MF_HL2', 'MF_CL2a', 'MF_CL2b',&
                    'MFcup_HL1', 'MFcup_CL1a', 'MFcup_CL1b','MFcup_HL2', 'MFcup_CL2a', 'MFcup_CL2b'/       
!
      IF(initial)THEN
         initial=.FALSE.
         count=0.0d0
         c3kfactor_ctl_val_sum=0.0d0
         IF(myrank.eq.0)THEN
            OPEN(813,file='ss_status_dp.dat')
            IF(cplmars.le.0)THEN
               WRITE(813,"(12a14)")(varname(i),i=1,6),(varname2(i),i=1,6)
            ELSE
               WRITE(813,"(22a14)")(varname(i),i=1,10),(varname2(i),i=1,12)
            ENDIF
         ENDIF   
      ENDIF
!      
      IF(myrank.eq.0)THEN
!
         IF(dmaster_pass.eq.1)THEN
            c3kfactor_ctl_val_sum=c3kfactor_ctl_val_sum+c3kfactor_ctl_val
            count=count+1
            c3kfactor_ctl_val=c3kfactor_ctl_val_sum/count
         ENDIF  
!         
         IF(cplmars.le.0)THEN
            WRITE(813,"(2e14.6,10e14.6)")time_mars,c3kfactor_ctl_val,(dpres(j),j=1,4),(c3mflow_junleg(j),j=1,6)     
         ELSE
            WRITE(813,"(2e14.6,20e14.6)")time_mars,c3kfactor_ctl_val,(dpres(j),j=1,8),(c3mflow_junleg(j),j=1,12)     
         ENDIF   
      ENDIF   
!      
      RETURN
      END SUBROUTINE print_ss_status_dp 
