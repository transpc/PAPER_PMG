!
!     SUBROUTINE rv_ihtc_vst(iflag)
      SUBROUTINE rv_ihtc_vst
!      
      USE VOL_DATA
      USE Zzone        , ONLY: ncell_fluid     
      USE Zparam       , ONLY: ndim 
      USE Zcoord1      , ONLY: xloc
      USE Zcoord2      , ONLY: cell_leng 
      USE Zconst2      , ONLY: hydraulicd,ggc         
      USE Znum_cell    , ONLY: i_neigh,neigh
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zrv_model    , ONLY: rv_valve
!            
!      
      IMPLICIT NONE      
!      
!     INTEGER iflag
      INTEGER i,j,k,j0
      INTEGER vst,vst_condition1,vst_condition2  
!
      REAL(8) agfvstst
      REAL(8) alphaglower
      REAL(8) alphag,alphal
      REAL(8) vg,vl
      REAL(8) dh
!
!.....Vertically Stratified
!
      DO i=1,ncell_fluid
         ggc=9.81d0 
         dh=hydraulicd(i)     
!         cell%length(i)=0.5d0*(cell_leng(i,1)+cell_leng(i,2)) 
         cell%length(i)=cell_leng(i,ndim) ! yjm mod
!
!        vst condition check      
         vst=0
         vst_condition1=0
         vst_condition2=0
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            k=neigh(j)
            IF(k.le.0.or.k.gt.ncell_fluid)CYCLE
            IF(i.eq.1.or.k.eq.0) CYCLE
!
!            IF(xloc(k,3).lt.xloc(i,3)) THEN
!               IF(cell%alphag(k).le.cell%alpha_bs(k)) THEN
!                  vst_condition1=1
!                  alphaglower=cell%alphag(k)
!                  alphag=cell%alphag(k)
!                  alphal=cell%alphal(k)
!                  vg=DABS(vg_o(ndim,i))
!                  vl=DABS(vl_o(ndim,i)) 
!               ELSE
!                  vst_condition1=0
!               ENDIF  
!!               IF(cell%alphag(k).lt.0.1d-2) THEN 
!!                  vst=1
!!               ENDIF
!            ENDIF               
!
!            IF(xloc(k,3).gt.xloc(i,3)) THEN
!               IF(cell%alphag(k).ge.cell%alpha_sa(k)) THEN
!!               IF(cell%alphag(k).ge.cell%alpha_bs(k)) THEN
!                  vst_condition2=1
!                  alphagupper=cell%alphag(k)
!               ELSE
!                  vst_condition2=0
!               ENDIF 
!            ENDIF 
!
            IF(xloc(k,3).lt.xloc(i,3)) THEN
               IF(cell%alphal(k).gt.1.d0-cell%alpha_bs(k)) THEN
!               IF(cell%alphal(k).gt.cell%alpha_sa(k)) THEN
                  vst_condition1=1
                  alphaglower=cell%alphag(k)
                  alphag=cell%alphag(k)
                  alphal=cell%alphal(k)
                  vg=DABS(vg_o(i,ndim))
                  vl=DABS(vl_o(i,ndim)) 
               ELSE
                  vst_condition1=0
               ENDIF  
            ENDIF 
!                
            IF(xloc(k,3).gt.xloc(i,3)) THEN
!!               IF(cell%alphal(k).gt.1.d0-cell%alpha_bs(k)) THEN
!               IF(cell%alphal(k).gt.cell%alpha_sa(k)) THEN
!!               IF(cell%alphal(k).lt.cell%alpha_sa(k)) THEN
!                  vst_condition2=0
!               ELSE
!                  vst_condition2=1 
!               ENDIF  
               IF(cell%alphal(k).lt.1.d0-cell%alpha_bs(k)) THEN
                  vst_condition2=1
               ELSE
                  vst_condition2=0 
               ENDIF  

            ENDIF 
!            
!            IF(cell%alphal(i).ge.1.d0-cell%alpha_sa(i)) THEN
!               vst_condition2=1
!               alphagupper=cell%alphag(i)
!            ELSE
!               vst_condition2=0
!            ENDIF                 
!            
         ENDDO  
!         
!        inlet cell number: LSJ
         IF(vst_condition1+vst_condition2.eq.2) vst=1
!         
!........valve model
!         
         IF(rv_valve.eq.1) CALL valve_model_vst(i,vst)
!         
         cell%vst(i)=vst
         cell%wf_vst(i)=0.d0        
         cell%ia_vst_st(i)=0.d0
         cell%ia_vst_sb(i)=0.d0
         IF(cell%vst(i).eq.1) THEN
            cell%wf_vst(i)=1.d0 
!
!........stratified surface portion         
            agfvstst=1.d0/cell%length(i)
            cell%ia_vst_st(i)=agfvstst
!!
!!........small bubble portion     
!            del_rho=DMAX1(1.d-5,cell%rhol(i)-cell%rhog(i))
!            tmp1=DMIN1(1.53d0,0.345*DSQRT(dh/DSQRT(cell%sigma(i)/(ggc*del_rho))))
!            tmp2=DSQRT(DSQRT(cell%sigma(i)*ggc*del_rho)/cell%rhol(i))
!            vgjs=tmp1*tmp2
!            jg=alphag*vg
!            jl=alphal*vl
!            alphagulev=jg/(1.2d0*(jg+jl)+vgjs)
!            alphagulev=DMIN1(cell%alphag(i),DMIN1(alphagupper,DMAX1(alphagulev,alphaglower)))
!!
!            fmix=1.d0
!            IF(alphagulev.lt.cell%alpha_bs(i)) THEN
!!
!               dalpha=DMAX1(1.d-5,alphagupper-alphagulev)
!               fmix=DMIN1(1.d0,DMAX1(0.d0,(alphagupper-cell%alphag(i))/dalpha))
!!               
!               alphagfsb=alphagulev*fmix
!               IF(iflag.eq.1) THEN
!                  agfvstsb=1.d0 !3.6d0/cell%dbb(i)*alphagfsb
!               ELSE
!                  agfvstsb=3.6d0/cell%dbb(i)*alphagfsb
!               ENDIF   
!               cell%ia_vst_sb(i)=agfvstsb
!            ENDIF
!          
         ENDIF         
!         
      ENDDO
!
!      DO i=2,ncell_fluid-1               
!         IF(cell%vst(i).eq.1) THEN
!            IF(cell%vst(i).eq.cell%vst(i-1)) THEN
!               cell%vst(i-1)=0
!               cell%ia_vst_st(i-1)=0.d0
!               cell%ia_vst_sb(i-1)=0.d0
!               cell%wf_vst(i-1)=0.d0                 
!               
!            ENDIF   
!            IF(cell%vst(i).eq.cell%vst(i+1)) THEN
!               cell%vst(i+1)=0
!               cell%ia_vst_st(i+1)=0.d0
!               cell%ia_vst_sb(i+1)=0.d0
!               cell%wf_vst(i+1)=0.d0                  
!            ENDIF            
!         ENDIF
!      ENDDO
      DO i=1,ncell_fluid      
         IF(cell%vst(i).eq.1) THEN    
            IF(i.ne.1) THEN        
               IF    (xloc(i,3).gt.xloc(i-1,3)) THEN
                  IF(cell%vst(i-1).eq.1) cell%vst(i-1)=0
               ELSEIF(xloc(i,3).lt.xloc(i-1,3)) THEN
                  IF(cell%vst(i-1).eq.1) cell%vst(i)=0
               ENDIF 
            ENDIF
            IF(i.ne.ncell_fluid) THEN
               IF(xloc(i,3).gt.xloc(i+1,3)) THEN
                  IF(cell%vst(i+1).eq.1) cell%vst(i+1)=0
               ELSEIF(xloc(i,3).lt.xloc(i+1,3)) THEN
                  IF(cell%vst(i+1).eq.1) cell%vst(i)=0
               ENDIF             
            ENDIF   
         ENDIF
      ENDDO
!
!      
      RETURN
      END SUBROUTINE rv_ihtc_vst
    
