!
      SUBROUTINE rv_flow_regime
!
!     This routine assign topology and calculates interfacial area according 
!     to the assigned topology.
!
      USE VOL_DATA  
      USE WALL_DATA                       
      USE Zzone      , ONLY: ncell_fluid
      USE Zparam     , ONLY: ndim,pi
      USE Zcoord2    , ONLY: cell_leng
      USE STM_TBL_cupid, ONLY: pcrit
      USE Zconst2    , ONLY: hydraulicd,ggc
      USE Zare       , ONLY: ar_gas,ar_liq
      USE Zvector    , ONLY: vl_o,vg_o
      USE Zwall_HTC  , ONLY: reflood,reflod        
      USE Zqvol      , ONLY: qporous_gas
!
      IMPLICIT NONE

      INTEGER i
!      
!.....pcrit is constant read in stread.f90 why change?
!     REAL(8),PARAMETER :: pcrit=22.06d6 !critical pressure of H2O     
      REAL(8),PARAMETER :: alpha_am=0.9999d0
      REAL(8) alphag,alphal,rhol,rhog,mu
      REAL(8) dh,length
      REAL(8) tgsat, del_rho
      REAL(8) gmg,gml,gm,vm,vtb,qwg
      REAL(8) pratio,twindo
      REAL(8) d_star
      REAL(8) alpha_bs_s,alpha_crit_f,alpha_crit_e
      REAL(8) vg,vl
!     REAL(8) agf0,agftb,agfbb,agfann,agfdrp,agfbub,agf,agfvstst,agfvstsb
      REAL(8) f9,f99
      REAL(8) ia_churn1,ia_invchn,ia_unvst,wfactor,ia_prechf,ia_posdry,ia_sum(ncell_fluid)
!      
!      LOGICAL liquid,gas,bubbly,slug,drift,churn,annular,mpr,invann,invchn,invslg,mist,mpo,vst
      LOGICAL prechf,poschf,posdry
!      
      ggc=9.8d0               !next
      DO i=1,ncell_fluid
         IF(reflood.eq.1) THEN
!         IF(mode(7).eq.1.or.mode(8).eq.1) THEN
            reflod=.TRUE.
         ELSE
            reflod=.FALSE.
         ENDIF
!         qwg=face%wall_fluxg_diff(i)
         qwg=qporous_gas(i)
         prechf=.true.
         poschf=.false.
         posdry=.false.        
!
         cell%ia_bubbly(i)=0.d0
         cell%ia_slug_tb(i)=0.d0
         cell%ia_slug_sb(i)=0.d0
         cell%ia_churn(i)=0.d0
         cell%ia_annular_drp(i)=0.d0
         cell%ia_annular_ann(i)=0.d0 
         cell%ia_mpr(i)=0.d0                                        
         cell%ia_invann_ann(i)=0.d0
         cell%ia_invann_sb(i)=0.d0
         cell%ia_invchn(i)=0.d0
         cell%ia_invslg_drp(i)=0.d0
         cell%ia_invslg_ann(i)=0.d0
         cell%ia_mist(i)=0.d0
         cell%ia_mpo(i)=0.d0
         cell%ia_VST(i)=0.d0
         cell%ia_vst_st(i)=0.d0
         cell%ia_vst_sb(i)=0.d0
!               
         ia_churn1=0.d0
         ia_invchn=0.d0
         ia_unvst=0.d0
         wfactor=0.d0
         ia_prechf=0.d0
         ia_posdry=0.d0
         ia_sum(i)=0.d0
!
!........local variable from cell index
!         
         alphag=cell%alphag(i)
         alphal=1.d0-cell%alphag(i)
         rhog=cell%rhog(i)
         rhol=cell%rhol(i)
         dh=hydraulicd(i)
         cell%length(i)=0.5d0*(cell_leng(i,1)+cell_leng(i,2))
         length=cell%length(i)
         mu=cell%lviscosl(i)
         del_rho=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i))
!
!........Basic parameter
!         
         vg=dsqrt(dot_product(vg_o(i,:),vg_o(i,:)))  
         vl=dsqrt(dot_product(vl_o(i,:),vl_o(i,:)))  
         gmg=ar_gas(i)*vg
         gml=ar_liq(i)*vl
         !gmg=cell%alphag(i)*cell%rhog(i)*vg   !remove later
         !gml=cell%alphal(i)*cell%rhol(i)*vl   !remove later
         gm=gmg+gml
         vm=gm*cell%rhomr(i)
         !vm=gm/(alphag*cell%rhog(i)+(1.d0-alphag)*cell%rhol(i))   !remove later
         vtb=0.35d0*dsqrt(ggc*hydraulicd(i)*del_rho/rhol)
!
!........Flow direction, cell%fdir
!        - 1: upward
!        - 2: downward
!        - 3: count-current
!
         cell%fdir(i)=0
!
!........co-current
         IF(vl_o(i,ndim)*vg_o(i,ndim).gt.0.d0)then
!
!...........upward
            IF(vl_o(i,ndim).gt.0.d0)then
               cell%fdir(i)=1
!
!...........downward
            ELSEIF(vl_o(i,ndim).lt.0.d0)then
               cell%fdir(i)=2
            ENDIF
!
!........count-current
         ELSEIF(vl_o(i,ndim)*vg_o(i,ndim).lt.0.d0)then
            cell%fdir(i)=3
!
!........zero component
         ELSE
            IF(vl_o(i,ndim).gt.0.d0 .or. vg_o(i,ndim).gt.0.d0) cell%fdir(i)=1
            IF(vl_o(i,ndim).lt.0.d0 .or. vg_o(i,ndim).lt.0.d0) cell%fdir(i)=2
            IF(vl_o(i,ndim).eq.0.d0 .or. vg_o(i,ndim).eq.0.d0) then
!               write(*,*) 'Flow direction is not defined because of zero component'
!               write(*,*) 'liquid velocity is',vl_o(i,ndim)
!               write(*,*) 'gas    velocity is',vg_o(i,ndim)
            ENDIF
         ENDIF
!
!........Define post-CHF/post-dryout         
!
         tgsat=cell%tg(i)-cell%ts(i)-1.0d0
!         tgsat=cell%tg(i)-cell%tst(i)-1.0d0         
!         reflod=.false.
!         IF(cell%ireflod(i))reflod=.true.
         IF(reflod) tgsat=tgsat-30.d0
!         IF((cell%qwg(i).gt.0.d0 .or. reflod) .and. tgsat.gt.0.d0)poschf=.true. !next
         IF((qwg.gt.0.d0 .or. reflod) .and. tgsat.gt.0.d0)poschf=.true. !next
!
         cell%wf_dry(i)=0.d0
         IF(poschf) THEN
            pratio=cell%p(i)/pcrit
            IF(pratio.lt.0.025d0) THEN
               twindo=0.06666667d0
            ELSEIF(pratio.gt.0.25d0) THEN
               twindo=0.016666667d0
            ELSE
               twindo=1.d0/(15.d0+200.d0*(pratio-0.025d0))
            ENDIF
            cell%wf_dry(i)=DMAX1(0.d0,DMIN1(1.d0,twindo*tgsat))
         ENDIF
         IF(reflod)then
            IF(poschf)cell%wf_dry(i)=DMAX1(0.d0,DMIN1(1.d0,(1.d0-DEXP(-0.5d0*tgsat))*1.0000454d0))
         ENDIF
!
!........Define regime criteria
!
!........alpha_bs
         d_star=dh*dsqrt(ggc*del_rho/cell%sigma(i))
         d_star=dmax1(22.22d0,dh*dsqrt(ggc*del_rho/cell%sigma(i)))
         alpha_bs_s=DMAX1(0.25d0*DMIN1(1.0d0,(0.045d0*d_star)**8), 1.d-3)
         alpha_bs_s=DMAX1(0.25d0,alpha_bs_s)         
         IF(gm.le.2000.0d0)then
            cell%alpha_bs(i)=alpha_bs_s
         ELSEIF(gm.gt.2000.d0 .and. gm.lt.3000.d0)then
            cell%alpha_bs(i)=alpha_bs_s+(0.5d0-alpha_bs_s)*1.d-3*(gm-2000.d0)
         ELSE
            cell%alpha_bs(i)=0.5d0
         ENDIF
!
!........special update of 'z_trans' for bubbly flow over 2000 kg/(m2.s) (mass flux)         
         IF(gm.gt.2000.d0) THEN
            cell%wf_dry(i)=DMAX1(0.d0,DMIN1(1.d0,cell%wf_dry(i)*(0.4d0-cell%alpha_bs(i))*10.d0))
         ENDIF
         IF(cell%wf_dry(i).eq.0.d0)then
            poschf=.false.
         ELSE
            poschf=.true.
         ENDIF
         IF(cell%wf_dry(i).ge.1.d0)posdry=.true.
!
!........alpha_cd
         cell%alpha_cd(i)=cell%alpha_bs(i)+0.2d0
!
!........alpha_sa
         IF(vg_o(i,ndim).ge.0.d0 .and. vl_o(i,ndim).ge.0.0d0) THEN  !upflow
            alpha_crit_f=DSQRT(ggc*dh*del_rho/rhog)/DMAX1(1.0d-9,vg_o(i,ndim))
         ELSE   !downflow or count-current flow
            alpha_crit_f=0.75d0
         ENDIF
         alpha_crit_e=3.2d0*(ggc*cell%sigma(i)*del_rho/(rhog*rhog))**0.25/DMAX1(1.0d-9,DABS(vg_o(i,ndim)))
         cell%alpha_sa(i)=DMAX1(0.8d0,DMIN1(alpha_crit_f,alpha_crit_e,0.9d0))
!
!........alpha_de
         cell%alpha_de(i)=DMAX1(cell%alpha_bs(i),cell%alpha_sa(i)-0.05d0)
!
!........alpha_gs
         f99=(alphag-cell%alpha_bs(i))/(cell%alpha_sa(i)-cell%alpha_bs(i))
         f9=dexp(-8.d0*f99)
         cell%alpha_gs(i)=cell%alpha_bs(i)*f9
!
!         CALL rv_ihtc_vst(i,vst)    !Added by LSJ for VST and wf_vst calculation     
!
!........check flow regimes and calculate interfacial area
!
!        1) pre-CHF
!           [0] Single-phase liquid
!           [1] bubbly (BBY)
!           [2] slug  (SLG)
!           [3] slg/anm transition (SLG/ANM)
!           [4] annular mist (ANM)
!           [5] pre-mist (MPR)
!        2) post-dryout
!           [6] inverted annular (IAN)
!           [7] inverted annular - inverted slug transition (IAN/ISLG)
!           [8] inverted slug (ISLG)
!           [9] mist (MST)
!           [10]pos-mist (MPO)
!        3) Transition
!           [11]bubbly[1] / inverted annular[6] (BBY-IAN)
!           [12]slug  [2] / inverted annular - inverted slug transition[7] (IAN/ISL-SLG)
!           [13]slug  [2] / inverted slug [8] (SLG/ISL)
!           [14]slg/anm transition[3] / inverted slug[8] (ISL-SLG/ANM)
!           [15]annular mist [4] / mist [9] (ANM/MST)
!
!........pre-CHF
         cell%regime(i)=0
         IF(.not.poschf)then
!...........Single-phse liquid
            IF(alphag.le.1.d-6)then
               cell%regime(i)=0
!...........Bubbly 
            ELSEIF(alphag.le.cell%alpha_bs(i))then
               cell%regime(i)=1
!...........slug 
            ELSEIF(alphag.lt.cell%alpha_de(i)) THEN
               cell%regime(i)=2
!...........slug/annular transition (churn)
            ELSEIF(alphag.lt.cell%alpha_sa(i)) THEN
               cell%regime(i)=3
!...........annular-mist (ANM)
            ELSEIF(alphag.lt.alpha_am) THEN
               cell%regime(i)=4
!...........mist pre-chf (MPR)
            ELSE
               cell%regime(i)=5
            ENDIF
!
!........post-dryout
         ELSEIF(posdry)then
!...........inverse-annular [IAN]
            IF(alphag.le.cell%alpha_bs(i))then
               cell%regime(i)=6
!...........IAN / ISLG transition
            ELSEIF(alphag.lt.cell%alpha_cd(i)) THEN
               cell%regime(i)=7
!...........inverted slug (ISLG)
            ELSEIF(alphag.lt.cell%alpha_sa(i)) THEN
               cell%regime(i)=8
!...........mist (MST)
            ELSEIF(alphag.lt.alpha_am) THEN
               cell%regime(i)=9
!...........mist post-chf (MPO)
            ELSE
               cell%regime(i)=10
            ENDIF
!
!...........pre-post transition
         ELSE
!...........bubbly[1] to inverted annular[6] (BBY-IAN)
            IF(alphag.lt.cell%alpha_bs(i)) THEN
               cell%regime(i)=11 
!...........slug[2] to inverted annular/inverted slug[7] (IAN/ISL-SLG)
            ELSEIF(alphag.lt.cell%alpha_cd(i)) THEN
               cell%regime(i)=12 
!
!...........slug[2] to inverted slug [8] (SLG/ISL)
            ELSEIF(alphag<=cell%alpha_de(i)) THEN
               cell%regime(i)=13 
!
!...........slug/annular-mist[3] to inverted slug[8] (ISL-SLG/ANM)
            ELSEIF(alphag<=cell%alpha_sa(i)) THEN
               cell%regime(i)=14 
!...........annular-mist[4] to mist[9] (ANM/MST)
            ELSE
               cell%regime(i)=15 
            ENDIF           
         ENDIF
!
      ENDDO
! 
!     CALL rv_ihtc_vst(1)       
      CALL rv_ihtc_vst
!
      RETURN
      END SUBROUTINE rv_flow_regime
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      REAL(8) FUNCTION weight1(alphag, alph_left, alph_right,vtran_max)
!
!.....Weighting factor by an independent variable's position and its interval.
!     When alphag is within its left and right boundary, the weighting factor is zero.
!     If alphag goes beyond its interval, the weighting factor linearly fades away to
!     zero within the fading width vt_l and vt_r.
!
      IMPLICIT NONE
!      
      REAL(8), INTENT(IN) :: alphag       ! a variable.
      REAL(8), INTENT(IN) :: alph_left    ! left limit for the variable.
      REAL(8), INTENT(IN) :: alph_right   ! right limit for the variable.
!
      REAL(8) :: vtran_max                         ! transition width
      !REAL(8), PARAMETER :: vtran_max = 1.0D-12    ! transition width
      REAL(8), PARAMETER :: vtran_min = 1.0D-12 ! transition width
      REAL(8) :: wl, wr, vt_l, vt_r
   
      IF(alph_left <= 0.0 .or. 1.0 <= alph_left)THEN
         vt_l = vtran_min
      ELSE
         vt_l = vtran_max
      ENDIF
   
      IF(alph_right <= 0.0 .or. 1.0 <= alph_right)THEN
         vt_r = vtran_min
      ELSE
         vt_r = vtran_max
      ENDIF
   
      wl = 0.5 + 0.5*( alphag - alph_left ) / vt_l
      wl = DMAX1( 0.0, DMIN1( wl , 1.0 ) )
      wr = 0.5 + 0.5*( alphag - alph_right ) / vt_r
      wr = DMAX1( 0.0, DMIN1( wr , 1.0 ) )
      weight1 = wl - wr
!      
      RETURN
      ENDFUNCTION weight1
