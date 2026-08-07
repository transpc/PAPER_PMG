!
      SUBROUTINE rv_ihtc_main(ag)
!
      USE VOL_DATA        , ONLY: cell
      USE Zzone           , ONLY: ncell_fluid    
      USE Zqvol           , ONLY: H_ig,H_il,H_gf,hig_o,hil_o,hgf_o,gamma_wall
      USE Zvector         , ONLY: ul_o,ug_o,vl_o,vg_o
      USE Zrv_htc         , ONLY: HTC_MIN,hig_big,hil_big  
      USE Zrv_weight      , ONLY: liquid,gas,                                    &
                                  bubbly,slug,churn,annular,mpr,                 &
                                  invann,invchn,invslg,mist,mpo,                 &
                                  wf_liquid,wf_gas,                              &
                                  wf_bubbly,wf_slug,wf_churn,wf_annular,wf_MPR,  &
                                  wf_invann,wf_invchn,wf_invslg,wf_mist,wf_MPO,  &
                                  pfnrgj      
      USE Zrv_flowmap     , ONLY: ia_tot,alp_tb
      USE Zconst2         , ONLY: dt,hydraulicd      
      USE Zparam          , ONLY: ndim
      USE Zcoord2         , ONLY: cell_leng     
!      
      IMPLICIT NONE
!
!.....Input
      REAL(8) ag(ncell_fluid)
!.....Local variables
      LOGICAL,SAVE :: limit      
      INTEGER :: i,reg_idx
      INTEGER :: churn_option,interregime_option      
      REAL(8) :: vg,vl,vr,vfg0,we,vfg1,vfg2,dh    
      REAL(8) :: agx,alx      
      REAL(8) :: wchurn,htc_Hil_s,htc_Hig_s,htc_Hil_a,htc_Hig_a 
      REAL(8) :: wfactor,ia_churn1,ia_invchn,ia_preCHF,ia_posDRY    
      REAL(8) :: alpha_bub         
      REAL(8) :: dstarc,hilsave,dtsl,dtsg,trmm1,trmm,xliqh,xvaph,percnt,avelfg
!.....Local arrays
      REAL(8),DIMENSION(17) :: htc_Hig,htc_Hil  !array=flow regimes    
      REAL(8),DIMENSION(ncell_fluid) :: timinv,tcouri,scrchh,yeta, &
                                        hig_save,hil_save
!
      LOGICAL,SAVE :: initial=.true.
      ia_churn1=0.d0      
      ia_invchn=0.d0 
!
      IF(initial) THEN
         ALLOCATE(ia_tot(ncell_fluid),alp_tb(ncell_fluid))
         ia_tot(:)=0.d0         
         alp_tb(:)=0.d0 
!
         timinv(:)  =0.d0
         tcouri(:)  =0.d0
         scrchh(:)  =0.d0
         yeta(:)    =0.d0
         hig_save(:)=0.d0
         hil_save(:)=0.d0 
!                          
         initial=.false.
      ENDIF   
!
      DO i=1,ncell_fluid
!
!........Initialize
!
         vg=ug_o(i)
         vl=ul_o(i)
         htc_Hig(:)=HTC_MIN  !array=flow regimes
         htc_Hil(:)=HTC_MIN
!
         dh=hydraulicd(i)         
!
!........MARS manual for bubble diameter (dbb)
!
         alpha_bub=MAX(cell%alphag(i),1.d-5)
         vr=SQRT(DOT_PRODUCT((vl_o(i,:)-vg_o(i,:)),(vl_o(i,:)-vg_o(i,:))))
         IF(cell%alphag(i).ge.1.d-5)then
            vfg0=vr
         ELSE
            vfg0=vr*cell%alphag(i)*1.d-5
         ENDIF
         we=5.d0
!         vfg1=we*cell%sigma(i)/(cell%rhol(i)*MIN(0.005d0*alpha_bub**(1.d0/3.d0),hydraulicd(i)))
         vfg1=we*cell%sigma(i)/(cell%rhol(i)*MIN(0.005d0*alpha_bub**(1.d0/3.d0),dh))
         vfg2=MAX(vfg0*vfg0,vfg1)
         cell%dbb(i)=we*cell%sigma(i)/MAX(1.d-8,(cell%rhol(i)*vfg2))            
!      
!........Define vertically stratified flow option
!
!         CALL rv_ihtc_vst(2)
!         vst=0
!         IF(cell%wf_vst(i).gt.0.d0.and.cell%wf_vst(i).le.1.d0) vst=1         
!
!........Weight for the flow regimes
!
         CALL rv_ihtc_weight(i)    
!
!........HTC: Liquid (single-phase)
!
         IF(liquid.eq.1) THEN
            reg_idx=1     
            IF(cell%Tl(i)<cell%Ts(i)) THEN           !subcooled liquid
               htc_Hil(reg_idx)=MAX(0.d0,HTC_MIN)
            ELSE                                     !superheated liquid
               htc_Hil(reg_idx)=hil_big
            ENDIF
            htc_Hig(reg_idx)=hig_big                  !big value
!
            IF(cell%vst(i).eq.1) THEN
               agx=ag(i)
               alx=1.d0-agx            
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx)) !vertically stratified
            ENDIF  
         ENDIF
!
!........HTC: Gas (single-phase)
!
         IF(gas.eq.1) THEN
            reg_idx=7         
            IF(cell%Tg(i)>cell%Ts(i)) THEN           ! superheated gas
               htc_Hig(reg_idx)=MAX(0.d0,HTC_MIN)
            ELSE                                     ! subcooled gas
               htc_Hig(reg_idx)=hig_big               ! SPACE uses Dittus-Boelter.
            END IF
            htc_Hil(reg_idx)=hil_big                  ! A BIG BIG BIG value
!
            IF(cell%vst(i).eq.1) THEN
               agx=ag(i)
               alx=1.d0-agx
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx)) !vertically stratified
            ENDIF   
         ENDIF
!
!........HTC: Bubbly
!
         IF(bubbly.eq.1) THEN
            reg_idx=2    
            CALL rv_ihtc_ia_preCHF(i,reg_idx)
!            CALL rv_ihtc_vst(i,vst)                         
!              
            agx=ag(i)
            alx=1.d0-agx
!           CALL rv_iht_bbl(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
            CALL rv_iht_bbl(i,agx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))  !vertically stratified
            ENDIF   
         ENDIF
!
!........HTC: Slug
!
         IF(slug.eq.1) THEN
            reg_idx=3
            CALL rv_ihtc_ia_preCHF(i,reg_idx)   
!            CALL rv_ihtc_vst(i,vst)                                                                         
!       
            agx=ag(i)
            alx=1.d0-agx
!           CALL rv_iht_slg(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
            CALL rv_iht_slg(i,agx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))  !vertically stratified
            ENDIF   
         ENDIF   
!
!........HTC: Slug-annular transition
!
         IF(churn.eq.1) THEN
            reg_idx=4   
            CALL rv_ihtc_ia_preCHF(i,reg_idx)            
!            CALL rv_ihtc_vst(i,vst)                                     
!                 
            agx=MAX(cell%alpha_de(i),MIN(ag(i),cell%alpha_sa(i)))
            churn_option=2
            SELECT CASE(churn_option)
               CASE(1)
!                 linear interpolation
                  wchurn=MAX(0.d0,MIN((agx-cell%alpha_de(i))/MAX(1.0D-6, cell%alpha_sa(i)-cell%alpha_de(i)),1.d0))
               CASE(2)
!                 MARS interpolation
                  wchurn=MAX(0.d0,MIN(20.d0*(agx-cell%alpha_de(i)),1.d0))
            ENDSELECT 
!  
!           slug
            agx=cell%alpha_de(i)
            alx=1.d0-agx
!           CALL rv_iht_slg(i,agx,alx,vg,vl,htc_Hil_s,htc_Hig_s)
            CALL rv_iht_slg(i,agx,vg,vl,htc_Hil_s,htc_Hig_s)
!
!           annular 
            agx=cell%alpha_sa(i)
            alx=1.d0-agx
            CALL rv_iht_anm(i,agx,alx,vg,vl,htc_Hil_a,htc_Hig_a)
!
            htc_Hig(reg_idx)=htc_Hig_s**(1.d0-wchurn)*htc_Hig_a**wchurn
            htc_Hil(reg_idx)=htc_Hil_s**(1.d0-wchurn)*htc_Hil_a**wchurn
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))   !vertically stratified
            ENDIF   
         ENDIF
!
!........HTC: Annular-mist
!
         IF(annular.eq.1) THEN
            reg_idx=5
            CALL rv_ihtc_ia_preCHF(i,reg_idx)
!            CALL rv_ihtc_vst(i,vst)                                     
!
            agx=ag(i)
            alx=1.d0-agx
            CALL rv_iht_anm(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))    !vertically stratified        
            ENDIF   
         ENDIF   
!
!........HTC: MPR
!
         IF(mpr.eq.1) THEN
            reg_idx=6
            CALL rv_ihtc_ia_preCHF(i,reg_idx)                                 
!            CALL rv_ihtc_vst(i,vst)                                     
!
            IF(cell%Tg(i)>cell%Ts(i)) THEN           !superheated gas
               htc_Hig(reg_idx)=MAX(0.d0,HTC_MIN)
            ELSE                                     
               htc_Hig(reg_idx)=hig_big               !subcooled gas (SPACE uses Dittus-Boelter)
            END IF
            htc_Hil(reg_idx)=hil_big                  !big value
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))     !vertically stratified
            ENDIF   
         ENDIF  
!                    
!........HTC: Inverted annular
!
         IF(invann.eq.1) THEN
            reg_idx=8   
            CALL rv_ihtc_ia_posDRY(i,reg_idx)
!            CALL rv_ihtc_vst(i,vst)                                     
!     
            agx=ag(i)            
            alx=1.d0-agx
!           CALL rv_iht_invann(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
            CALL rv_iht_invann(i,agx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))      !vertically stratified
            ENDIF   
         ENDIF
!
!........HTC: Inverted Annular-Inverted Slug Transition
!
         IF(invchn.eq.1) THEN
            reg_idx=9
            CALL rv_ihtc_ia_posDRY(i,reg_idx)            
!            CALL rv_ihtc_vst(i,vst)                                     
!    
            agx=MAX(cell%alpha_bs(i),MIN(ag(i),cell%alpha_cd(i)))
            alx=1.0-agx
!            
            churn_option=2
            SELECT CASE(churn_option)
               CASE(1)
!                 linear interpolation
                  wchurn=MAX(0.d0,MIN((agx-cell%alpha_bs(i))/MAX(1.0D-6,cell%alpha_cd(i)-cell%alpha_bs(i)),1.d0))
               CASE(2)
!                 MARS interpolation
                  wchurn=MAX(0.d0,MIN(20.d0*(agx-cell%alpha_bs(i)),1.d0))
            ENDSELECT 
!  
!           inverted annular
            agx=cell%alpha_bs(i)
            alx=1.d0-agx
!           CALL rv_iht_invann(i,agx,alx,vg,vl,htc_Hil_a,htc_Hig_a)
            CALL rv_iht_invann(i,agx,vg,vl,htc_Hil_a,htc_Hig_a)
!
!           inverted slug
            agx=cell%alpha_cd(i)
            alx=1.d0-agx
!           CALL rv_iht_invslg(i,agx,alx,vg,vl,htc_Hil_s,htc_Hig_s)
            CALL rv_iht_invslg(i,agx,vg,vl,htc_Hil_s,htc_Hig_s)
!
            htc_Hig(reg_idx)=htc_Hig_a**(1.d0-wchurn)*htc_Hig_s**wchurn
            htc_Hil(reg_idx)=htc_Hil_a**(1.d0-wchurn)*htc_Hil_s**wchurn
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
            ENDIF   
         ENDIF          
!
!........HTC: Inverted slug
!
         IF(invslg.eq.1) THEN
            reg_idx=10
            CALL rv_ihtc_ia_posDRY(i,reg_idx)
!            CALL rv_ihtc_vst(i,vst)                                     
!     
            agx=ag(i)
            alx=1.d0-agx
!           CALL rv_iht_invslg(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
            CALL rv_iht_invslg(i,agx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))      !vertically stratified
            ENDIF   
         ENDIF     
!
!........Mist
!
         IF(mist.eq.1) THEN
            reg_idx=11
            CALL rv_ihtc_ia_posDRY(i,reg_idx) 
!            CALL rv_ihtc_vst(i,vst)                                                        
!   
            agx=ag(i)            
            alx=1.d0-agx
            CALL rv_iht_mist(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))      !vertically stratified
            ENDIF   
         ENDIF              
!
!........MPO
!
         IF(mpo.eq.1) THEN
            reg_idx=12
            CALL rv_ihtc_ia_posDRY(i,reg_idx)
!            CALL rv_ihtc_vst(i,vst)                                     
!       
            IF(cell%Tg(i)>cell%Ts(i)) THEN           !superheated gas
               htc_Hig(reg_idx)=MAX(0.d0,HTC_MIN)
            ELSE                                     
               htc_Hig(reg_idx)=hig_big               !subcooled gas(SPACE uses Dittus-Boelter)
            END IF
            htc_Hil(reg_idx)=hil_big                  !big value
!
            IF(cell%vst(i).eq.1) THEN
               CALL rv_iht_vstrat(i,agx,alx,vg,vl,htc_Hil(reg_idx),htc_Hig(reg_idx))      !vertically stratified
            ENDIF   
         ENDIF 
!
!........INTER-REGIME SMOOTHING --> Final Heat Transfer: H_ig, H_il
!
         pfnrgj=cell%wf_dry(i)
         interregime_option=2
         SMOOTHING: SELECT CASE (interregime_option)
            CASE(1)
!              interpolation : linear
               H_ig(i) = (wf_liquid*htc_Hig( 1)+wf_gas    *htc_Hig( 6)+wf_bubbly*htc_Hig( 2)+wf_slug  *htc_Hig( 3)+ &
                          wf_churn *htc_Hig( 4)+wf_annular*htc_Hig( 5)+wf_mpr   *htc_Hig( 7))*(1.d0-pfnrgj)        &
                        +(wf_liquid*htc_Hig( 1)+wf_gas    *htc_Hig( 6)+wf_invann*htc_Hig( 8)+wf_invchn*htc_Hig( 9)+ &
                          wf_invslg*htc_Hig(10)+wf_mist   *htc_Hig(11)+wf_mpo   *htc_Hig(12))*pfnrgj
               H_il(i) = (wf_liquid*htc_Hil( 1)+wf_gas    *htc_Hil( 6)+wf_bubbly*htc_Hil( 2)+wf_slug  *htc_Hil( 3)+ &
                          wf_churn *htc_Hil( 4)+wf_annular*htc_Hil( 5)+wf_mpr   *htc_Hil( 7))*(1.d0-pfnrgj)        &
                        +(wf_liquid*htc_Hil( 1)+wf_gas    *htc_Hil( 6)+wf_invann*htc_Hil( 8)+wf_invchn*htc_Hil( 9)+ &
                          wf_invslg*htc_Hil(10)+wf_mist   *htc_Hil(11)+wf_mpo   *htc_Hil(12))*pfnrgj
            CASE(2)
!              interpolation : power-law
               H_ig(i) = (EXP( wf_liquid*LOG(htc_Hig( 1))+wf_gas   *LOG(htc_Hig( 7))+wf_bubbly *LOG(htc_Hig( 2))+ &
                               wf_slug  *LOG(htc_Hig( 3))+wf_churn *LOG(htc_Hig( 4))+wf_annular*LOG(htc_Hig( 5))+ &
                               wf_mpr   *LOG(htc_Hig( 6)) ))*(1.d0-pfnrgj)                                         &
                        +(EXP( wf_liquid*LOG(htc_Hig( 1))+wf_gas   *LOG(htc_Hig( 7))+wf_invann *LOG(htc_Hig( 8))+ & 
                               wf_invchn*LOG(htc_Hig( 9))+wf_invslg*LOG(htc_Hig(10))+wf_mist   *LOG(htc_Hig(11))+ &
                               wf_mpo   *LOG(htc_Hig(12)) ))*pfnrgj
               H_il(i) = (EXP( wf_liquid*LOG(htc_Hil( 1))+wf_gas   *LOG(htc_Hil( 7))+wf_bubbly *LOG(htc_Hil( 2))+ &
                               wf_slug  *LOG(htc_Hil( 3))+wf_churn *LOG(htc_Hil( 4))+wf_annular*LOG(htc_Hil( 5))+ &
                               wf_mpr   *LOG(htc_Hil( 6)) ))*(1.d0-pfnrgj)                                         &
                        +(EXP( wf_liquid*LOG(htc_Hil( 1))+wf_gas   *LOG(htc_Hil( 7))+wf_invann *LOG(htc_Hil( 8))+ &
                               wf_invchn*LOG(htc_Hil( 9))+wf_invslg*LOG(htc_Hil(10))+wf_mist   *LOG(htc_Hil(11))+ &
                               wf_mpo   *LOG(htc_Hil(12)) ))*pfnrgj
         ENDSELECT SMOOTHING 
!                          
!........Total interfacial area
!
         IF(churn.eq.1)then
             wfactor=(cell%alphag(i)-cell%alpha_de(i))/(cell%alpha_sa(i)-cell%alpha_de(i))
             ia_churn1=(cell%ia_slug_sb(i)+cell%ia_slug_tb(i))*(1.0-wfactor)+ &
                       (cell%ia_annular_ann(i)+cell%ia_annular_drp(i))*wfactor
         ENDIF
         IF(invchn.eq.1)then
             wfactor=(cell%alphag(i)-cell%alpha_bs(i))/(cell%alpha_cd(i)-cell%alpha_bs(i))
             ia_invchn=(cell%ia_invann_sb(i)+cell%ia_invann_ann(i))* &
                       (1.0-wfactor)+(cell%ia_invslg_drp(i)+cell%ia_invslg_ann(i))*wfactor
         ENDIF
         ia_preCHF = wf_bubbly * cell%ia_bubbly(i)                              &
                   + wf_slug   *(cell%ia_slug_sb(i)    +cell%ia_slug_tb(i))     &
                   + wf_churn  * ia_churn1                                      &
                   + wf_annular*(cell%ia_annular_ann(i)+cell%ia_annular_drp(i)) &
                   + wf_MPR    * cell%ia_mpr(i)         
         ia_posDRY = wf_invann *(cell%ia_invann_sb(i)  +cell%ia_invann_ann(i))  &
                   + wf_invchn * ia_invchn                                      &
                   + wf_invslg *(cell%ia_invslg_drp(i) +cell%ia_invslg_ann(i))  &
                   + wf_mist   * cell%ia_mist(i)                                &
                   + wf_MPO    * cell%ia_mpo(i)
         ia_tot(i) = ia_preCHF*(1.d0-pfnrgj)+ia_posDRY*pfnrgj   
!
!........Direct heat transfer
!
         H_gf(i)=H_ig(i)          
!         
         hig_save(i)=h_ig(i)
         hil_save(i)=h_il(i)  
!                  
      ENDDO 
!
      DO i=1,ncell_fluid
         dh=hydraulicd(i)
         dtsl=cell%ts(i)-cell%tl(i)
         dtsg=cell%ts(i)-cell%tg(i)
         dstarc=dh*DSQRT(9.80665d0*MAX(1.d0,cell%rhol(i)-cell%rhog(i))/cell%sigma(i))
         timinv(i)=DSQRT(MAX(9.80665d0,0.516d0*dstarc)/dh)*dt  !   9.80665/19=0.516
         tcouri(i)=MIN(DABS(vg_o(i,ndim)),DABS(vl_o(i,ndim)))*0.7d0*dt/cell_leng(i,ndim)   ! cell_leng(3,i) : z-direction cell length
         scrchh(i)=MAX(0.10536d0,(1.d-7+MIN(DABS(vg_o(i,ndim)),DABS(vl_o(i,ndim))))/MAX(DABS(vg_o(i,ndim)),DABS(vl_o(i,ndim)),1.d-7))
!
!        hil time smoothing
         IF(hil_o(i)*h_il(i).gt.0.d0.and.h_il(i).gt.hil_o(i))THEN
            yeta(i)=EXP(-MIN(0.693d0,MAX(tcouri(i)*MAX(0.01d0,cell%alphal(i)),1.d0-MIN(1.d0,cell%alphal(i)*1.d7),MIN(timinv(i),scrchh(i)))))
            yeta(i)=yeta(i)*(1.d0+MAX(-0.5d0,0.25d0*MIN(0.d0,dtsl)))
            hilsave=h_il(i)
            h_il(i)=h_il(i)*(hil_o(i)/h_il(i))**yeta(i)
            IF(dtsl.lt.0.d0) h_il(i)=0.5d0*(h_il(i)+MAX(h_il(i),hilsave))
         ENDIF
!
!        hig time smoothing
         IF(hig_o(i)*h_ig(i).gt.0.d0.and.h_ig(i).gt.hig_o(i))THEN
            yeta(i)=EXP(-MIN(0.693d0,MAX(tcouri(i)*MAX(0.01d0,cell%alphag(i)),1.d0-MIN(1.d0,cell%alphag(i)*1.d5),MIN(timinv(i),scrchh(i)))))
            yeta(i)=yeta(i)*(1.d0-2.5d0*MAX(0.d0,MIN(0.2d0,dtsg)))
            h_ig(i)=h_ig(i)*(hig_o(i)/h_ig(i))**yeta(i)
         ENDIF
!
!        hgf time smoothing
         IF(hgf_o(i)*h_gf(i).gt.0.d0.and.h_gf(i).gt.hgf_o(i))THEN
            yeta(i)=EXP(-MIN(0.693d0,MAX(tcouri(i)*MAX(0.01d0,cell%alphag(i)),1.d0-MIN(1.d0,cell%alphag(i)*1.d5),MIN(timinv(i),scrchh(i)))))
            yeta(i)=yeta(i)*(1.d0-2.5d0*MAX(0.d0,MIN(0.2d0,dtsg)))
            h_gf(i)=h_gf(i)*(hgf_o(i)/h_gf(i))**yeta(i)
         ENDIF
!
      IF(1)THEN !This makes trouble during single phase water flow
         h_il(i)=MAX(0.d0,h_il(i)) !0.1d0
         h_ig(i)=MAX(0.d0,h_ig(i)) !10.d0
         h_gf(i)=MAX(0.d0,h_gf(i)) !10.d0
      ENDIF
!
!        check on 0.5 vaporization/condensation limit
!
         limit=.false.
         IF(h_il(i)*h_ig(i).gt.0.d0)THEN
            trmm1=-(H_ig(i)*cell%pps(i)/cell%p(i)*(cell%ts(i)-cell%tg(i))+H_il(i)*(cell%ts(i)-cell%tl(i)))
            IF(trmm1.ge.0.d0)THEN
               xliqh=cell%hl(i)
               xvaph=cell%hgsat(i)
            ELSE
               xliqh=cell%hlsat(i)
               xvaph=cell%hg(i)
            ENDIF
            trmm=(gamma_wall(i)+trmm1/MAX(1.d-12,xvaph-xliqh))*dt
            IF(trmm.ne.0.d0)THEN
               IF(trmm.gt.0.d0)THEN
                  trmm1=0.5d0*cell%alphal(i)*cell%rhol(i)
                  IF(trmm.gt.trmm1.and.trmm1.ne.0.d0)THEN
                     avelfg=trmm1/trmm
                     limit=.true.
                  ENDIF
               ELSE
                  trmm1=0.5d0*cell%alphag(i)*cell%rhog(i)
                  IF(-trmm.gt.trmm1.and.trmm1.ne.0.d0)THEN
                     avelfg=-trmm1/trmm
                     limit=.true.
                  ENDIF
               ENDIF
!
            IF(0)THEN
               IF(limit)THEN
                  IF(cell%alphag(i).lt.0.5d0)THEN
                     percnt=MAX(0.01d0,MIN(1.d0,1.d-5/cell%dtsdp(i)))
                     h_il(i)=MAX(0.1d0,percnt*h_il(i)*avelfg)
                     IF(cell%tg(i).lt.cell%ts(i))THEN
                        h_ig(i)=MAX(10.d0,percnt*h_ig(i)*avelfg)
                     ELSE
                        h_ig(i)=MAX(1.d10,h_ig(i),hig_save(i))
                     ENDIF
                  ELSE
                     IF(cell%tg(i).gt.cell%ts(i)) h_ig(i)=MAX(10.d0,h_ig(i)*avelfg)
                  ENDIF
               ENDIF
            ELSE ! Below is not effective.             
               IF(limit)THEN
                  IF(cell%alphag(i).lt.0.5d0)THEN
                     percnt=MAX(0.01d0,MIN(1.d0,1.d-5/cell%dtsdp(i)))
                     h_il(i)=MAX(0.d0,percnt*h_il(i)*avelfg)
                     IF(cell%tg(i).lt.cell%ts(i))THEN
                        h_ig(i)=MAX(0.d0,percnt*h_ig(i)*avelfg)
                     ELSE
                        h_ig(i)=MAX(1.d10,h_ig(i),hig_save(i))
                     ENDIF
                  ELSE
                     IF(cell%tg(i).gt.cell%ts(i)) h_ig(i)=MAX(0.d0,h_ig(i)*avelfg)
                  ENDIF
               ENDIF
            ENDIF               
!                              
            ENDIF
         ENDIF
!         
         hil_o(i)=h_il(i)
         hig_o(i)=h_ig(i)
         hgf_o(i)=h_gf(i)
!
      ENDDO      
!
      END SUBROUTINE rv_ihtc_main
