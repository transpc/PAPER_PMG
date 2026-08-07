! 
      SUBROUTINE udfn_calc_HTC_int_i(ag,al)
!
!     Calculate heat transfer coefficients
!
      USE VOL_DATA                 
      USE Zparam       , ONLY: ndim
      USE Zconst1      , ONLY: vv_prob
      USE Zdaint_ag    , ONLY: daint1_ag,daint2_ag
      USE Zdhda        , ONLY: dhldag,dhgdag
      USE Zflowregime  , ONLY: alphag_cm,alphag_bc,gamma_1,gamma_2
      USE Zmodel       , ONLY: h_fg_coeff_a,h_fg_coeff_b,h_il_coeff_a,h_il_coeff_b,  &
                                h_ig_coeff_a,h_ig_coeff_b,h_il_min,h_ig_min,h_fg_min, &
                                h_il_min_ag99,h_il_min_user
      USE Zmpi         , ONLY: ncell_fp
      USE Zqvol        , ONLY: h_il,h_ig,h_gf
      USE Ztimecon     , ONLY: alpha_min,smac 
      USE Zmodel       , ONLY: i_weight
      USE Zuserdefined , ONLY: udfl_flashing_hif
      USE Zvector      , ONLY: vrel_o,vl_o,vg_o,ul_o,ug_o
      USE Zvoid        , ONLY: gamma_void,dagdx
      USE Zzone        , ONLY: ncell_fluid
!      
      IMPLICIT NONE
!
      INTEGER i
!      
      REAL(8) H_il_i,H_ig_i,H_gf_i
      REAL(8) ag(ncell_fp),al(ncell_fp)
      REAL(8) Nu,Re,Pr,st,hig1,hig2,hig3,hil1,hil2,hil3,delv
      REAL(8) weight,H_il_1,H_il_2,H_il_3,H_ig_1,H_ig_2,H_ig_3
      REAL(8) dHldag_1,dHldag_2,dHgdag_1,dHgdag_2,dHldag_3,dHgdag_3
      REAL(8) utangent(ndim)  
      REAL(8) quala_coeff
          
!
      IF(vv_prob.eq.'DIVA-NEW'.or.vv_prob.eq.'psbt_sngl'.or.vv_prob.eq.'flashing')THEN
!      
!.........udfl_calc_HTC_int_i: To determine the H_ig
!  
           quala_coeff=1.e5
               
!DIR$ SIMD
         DO i=1,ncell_fluid
            delv=vrel_o(i)
!
!...........Bubble
!
            Re=cell%rhol(i)*delv*cell%D1(i)/cell%lviscosl(i)
            Pr=cell%lviscosl(i)*cell%cpl(i)/cell%lcondl(i)
            Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
            hil1=cell%lcondl(i)*Nu/cell%D1(i)
            hig1=1000.0d0
!
!...........Droplet
!
            Re=cell%rhog(i)*delv*cell%D2(i)/cell%lviscosg(i)
            Pr=cell%lviscosg(i)*cell%cpg(i)/cell%lcondg(i)
            Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
            hil2=cell%lcondl(i)*Nu/cell%D2(i)
            hig2=1000.0d0                  
!
!...........Limit value
!
            hig1=dmax1(hig1,1.0d0)
            hig2=dmax1(hig2,1.0d0)
            hil1=dmax1(hil1,1.0d0)
            hil2=dmax1(hil2,1.0d0)
            hil3=dmax1(hil1,1.0d0)      !PAFS
            hig3=dmax1(hig1,1.0d0)      !PAFS
!
            H_il_1=hil1*cell%aint1(i) 
            H_il_2=hil2*cell%aint2(i) 
            H_ig_1=hig1*cell%aint1(i) 
            H_ig_2=hig2*cell%aint2(i) 
            H_il_3=hil3*cell%aint3(i)   !PAFS
            H_ig_3=hig3*cell%aint3(i)   !PAFS
!
            dHldag_1=hil1*daint1_ag(i)
            dHldag_2=hil2*daint2_ag(i)
            dHgdag_1=hig1*daint1_ag(i)
            dHgdag_2=hig2*daint2_ag(i)
            dHldag_3=0.0d0
            dHgdag_3=0.0d0
!
!...........H_ig,H_il for each regime
!
            IF(cell%regime(i).eq.11)THEN    
               H_il_i=H_il_1
               H_ig_i=H_ig_1
               dhldag(i)=dHldag_1
               dhgdag(i)=dHgdag_1      
            ELSEIF(cell%regime(i).eq.13)THEN
               H_il_i=H_il_2
               H_ig_i=H_ig_2
               dhldag(i)=dHldag_2
               dhgdag(i)=dHgdag_2
            ELSEIF(cell%regime(i).eq.12)THEN
               weight=(alphag_cm-ag(i))/(alphag_cm-alphag_bc)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_2
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_2
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_2
            ELSEIF(cell%regime(i).eq.3)THEN
               H_gf_i=H_fg_coeff_a+H_fg_coeff_b*max(1.0e-2,ag(i)*al(i))
               H_il_i=H_il_3
               H_ig_i=H_ig_3
               dhldag(i)=dHldag_3
               dhgdag(i)=dHgdag_3   
            ELSEIF(cell%regime(i).eq.21)THEN
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_3
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_3     
            ELSEIF(cell%regime(i).eq.23)THEN
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               H_il_i=weight*H_il_2+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_2+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dHldag_2+(1.0d0-weight)*dHldag_3
               dhgdag(i)=weight*dHgdag_2+(1.0d0-weight)*dHgdag_3              
            ELSEIF(cell%regime(i).eq.22)THEN
               weight=(alphag_cm-cell%alphag(i))/(alphag_cm-alphag_bc)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_2
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_2
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_2   
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               H_il_i=weight*H_il_i+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_i+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dhldag(i)+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dhgdag(i)+(1.0d0-weight)*dHgdag_2             
            ELSE
               H_il_i=H_il_coeff_a+H_il_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_ig_i=H_ig_coeff_a+H_ig_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_gf_i=H_fg_coeff_a+H_fg_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_il_i=H_il_1 
               H_ig_i=H_ig_1 
               dhldag(i)=0.0d0
               dhgdag(i)=0.0d0        
            ENDIF
            H_il_i=dmax1(H_il_min,H_il_i)
            H_ig_i=dmax1(H_ig_min,cell%quala(i)*quala_coeff,H_ig_i)  !quala_coeff is determined at UDF.
            H_gf_i=dmax1(H_fg_min,H_ig_i)
!
!...........H_ig,H_il modification for ag>0.999 or ag<0.001
!
            IF(ag(i).lt.0.001d0)THEN
               H_ig_i=dmax1(H_ig_i,1.0e6)
            ENDIF        
!
!...........H_ig,H_il for flashing
!
            IF(udfl_flashing_hif)THEN
               IF(cell%alphag(i).le.0.15d0)THEN
                  H_il_i=dmax1(H_il_i,1.0e7*(cell%tl(i)-cell%ts(i)))
               ELSEIF(cell%alphag(i).le.0.25d0)THEN
                  H_il_i=dmax1(H_il_i,(1.0e7-(1.0e7-1.0e3)/(0.25d0-0.15d0)*(cell%alphag(i)-0.15d0))*(cell%tl(i)-cell%ts(i)))
               ELSE
                  H_il_i=dmax1(H_il_i,1.0e3*(cell%tl(i)-cell%ts(i)))
               ENDIF
            ENDIF 
!               
!..............Subcooled water only
!
            IF(ag(i).le.2.0d0*alpha_min)THEN
               IF(smac.eq.3)THEN
                  H_ig_i=0.0d0
                  H_gf_i=0.0d0
               ENDIF
               IF(cell%tl(i).lt.cell%ts(i)) H_il_i=0.0d0
            ENDIF
!               
!...........Superheated steam only
!
            IF(ag(i).ge.(1.0d0-2.0d0*alpha_min))THEN
               IF(smac.eq.3) H_il_i=0.0d0
               IF(cell%tg(i).ge.cell%ts(i))THEN
                  H_ig_i=0.0d0
                  H_gf_i=0.0d0
               ENDIF
            ENDIF
!               
!...........Update H_ig,H_il,H_gf
!
            H_il(i)=H_il_i
            H_ig(i)=H_ig_i
            H_gf(i)=H_gf_i
         ENDDO
!      
      ELSE 
!      
!DIR$ SIMD
         DO i=1,ncell_fluid
            delv=vrel_o(i)
!
!...........Bubble
!
            Re=cell%rhol(i)*delv*cell%D1(i)/cell%lviscosl(i)
            Pr=cell%lviscosl(i)*cell%cpl(i)/cell%lcondl(i)
            Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
            hil1=cell%lcondl(i)*Nu/cell%D1(i)
            hig1=1000.0d0
!
!...........Bubble
!
            Re=cell%rhog(i)*delv*cell%D2(i)/cell%lviscosg(i)
            Pr=cell%lviscosg(i)*cell%cpg(i)/cell%lcondg(i)
            Nu=2.0d0+0.6d0*Re**0.5d0*Pr**0.3d0
            hil2=cell%lcondl(i)*Nu/cell%D2(i)
            hig2=cell%lcondg(i)*Nu/cell%D2(i)
            hig2=1000.0d0
!
!...........Interface
!
            utangent(:)=(vl_o(i,:)-vg_o(i,:))-dot_product(vl_o(i,:)-vg_o(i,:),dagdx(i,:))  &
                        /dmax1(alpha_min,dsqrt(dot_product(dagdx(i,:),dagdx(i,:))))            &
                        *dagdx(i,:)/dmax1(alpha_min,dsqrt(dot_product(dagdx(i,:),dagdx(i,:))))
            st=0.0045d0*((cell%rhog(i)*ug_o(i)*cell%lviscosl(i))/(cell%rhol(i)*DMAX1(ul_o(i),1d-3)*cell%lviscosg(i)))**(1.0d0/3.0d0)
            !st=0.0045d0*((cell%rhog(i)*ug_o(i)*cell%lviscosl(i))/(cell%rhol(i)*ul_o(i)*cell%lviscosg(i)))**(1.0d0/3.0d0)
            hil3=st*cell%rhol(i)*cell%cpl(i)*dsqrt(dot_product(utangent(:),utangent(:)))
            hig3=st*cell%rhog(i)*cell%cpg(i)*dsqrt(dot_product(utangent(:),utangent(:)))
            hil3=dmax1(hil3,1.0d0) 
            hig3=dmax1(hig3,1.0d0) 
!
!...........Limit value
!
            hig1=dmax1(hig1,1.0d0)
            hig2=dmax1(hig2,1.0d0)
            hil1=dmax1(hil1,1.0d0)
            hil2=dmax1(hil2,1.0d0)
!        
            H_il_1=hil1*cell%aint1(i) 
            H_il_2=hil2*cell%aint2(i) 
            H_ig_1=hig1*cell%aint1(i) 
            H_ig_2=hig2*cell%aint2(i) 
            H_il_3=hil3*cell%aint3(i)   
            H_ig_3=hig3*cell%aint3(i)   
!        
            dHldag_1=hil1*daint1_ag(i)
            dHldag_2=hil2*daint2_ag(i)
            dHgdag_1=hig1*daint1_ag(i)
            dHgdag_2=hig2*daint2_ag(i)
            dHldag_3=0.0d0
            dHgdag_3=0.0d0
!
!...........H_ig,H_il for each regime
!        
            IF(cell%regime(i).eq.11)THEN
               H_il_i=H_il_1
               H_ig_i=H_ig_1
               dhldag(i)=dHldag_1
               dhgdag(i)=dHgdag_1      
            ELSEIF(cell%regime(i).eq.13)THEN
               H_il_i=H_il_2
               H_ig_i=H_ig_2
               dhldag(i)=dHldag_2
               dhgdag(i)=dHgdag_2
            ELSEIF(cell%regime(i).eq.12)THEN
               weight=(alphag_cm-ag(i))/(alphag_cm-alphag_bc)
               IF(i_weight.ge.1) CALL udfn_weight1(weight,ag(i),alphag_cm,alphag_bc)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_2
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_2
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_2
            ELSEIF(cell%regime(i).eq.3)THEN
               H_il_i=H_il_3
               H_ig_i=H_ig_3
               dhldag(i)=dHldag_3
               dhgdag(i)=dHgdag_3   
            ELSEIF(cell%regime(i).eq.21)THEN
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               IF(i_weight.ge.1) CALL udfn_weight1(weight,gamma_void(i),gamma_2,gamma_1)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_3
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_3     
            ELSEIF(cell%regime(i).eq.23)THEN
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               IF(i_weight.ge.1) CALL udfn_weight1(weight,gamma_void(i),gamma_2,gamma_1)
               H_il_i=weight*H_il_2+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_2+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dHldag_2+(1.0d0-weight)*dHldag_3
               dhgdag(i)=weight*dHgdag_2+(1.0d0-weight)*dHgdag_3              
            ELSEIF(cell%regime(i).eq.22)THEN
               weight=(alphag_cm-cell%alphag(i))/(alphag_cm-alphag_bc)
               IF(i_weight.ge.1) CALL udfn_weight1(weight,ag(i),alphag_cm,alphag_bc)
               H_il_i=weight*H_il_1+(1.0d0-weight)*H_il_2
               H_ig_i=weight*H_ig_1+(1.0d0-weight)*H_ig_2
               dhldag(i)=weight*dHldag_1+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dHgdag_1+(1.0d0-weight)*dHgdag_2   
               weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
               IF(i_weight.ge.1) CALL udfn_weight1(weight,gamma_void(i),gamma_2,gamma_1)
               H_il_i=weight*H_il_i+(1.0d0-weight)*H_il_3
               H_ig_i=weight*H_ig_i+(1.0d0-weight)*H_ig_3
               dhldag(i)=weight*dhldag(i)+(1.0d0-weight)*dHldag_2
               dhgdag(i)=weight*dhgdag(i)+(1.0d0-weight)*dHgdag_2             
            ELSE
               H_il_i=H_il_coeff_a+H_il_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_ig_i=H_ig_coeff_a+H_ig_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_gf_i=H_fg_coeff_a+H_fg_coeff_b*dmax1(1.0e-2,ag(i)*al(i))
               H_il_i=H_il_1 
               H_ig_i=H_ig_1 
               dhldag(i)=0.0d0
               dhgdag(i)=0.0d0        
            ENDIF
            H_il_i=dmax1(H_il_min,H_il_i)
            H_ig_i=dmax1(H_ig_min,cell%quala(i)*1.e7,H_ig_i)
            H_gf_i=dmax1(H_fg_min,H_ig_i)
!
!...........H_ig,H_il modification for ag>0.999
!         
            IF(ag(i).gt.0.999d0)THEN
               H_il_i=dmax1(H_il_i,H_il_min_ag99)
            ENDIF
!
!...........H_ig,H_il for flashing
!         
            IF(cell%tl(i).gt.cell%ts(i))THEN
               IF(udfl_flashing_hif)THEN
                  CALL udfn_flashing_hif1(i,H_il_i,H_ig_i,ag(i))
               ELSE
                  H_il_i=dmax1(H_il_i,H_il_min_user) !pfactor-subo-speed
               ENDIF
            ENDIF
!               
!...........Subcooled water only
!         
            IF(ag(i).le.2.0d0*alpha_min)THEN
               IF(smac.eq.3) THEN
                  H_ig_i=0.0d0
                  H_gf_i=0.0d0
               ENDIF
               IF(cell%tl(i).lt.cell%ts(i)) H_il_i=0.0d0
           ENDIF
!               
!..........Superheated steam only
!         
           IF(ag(i).ge.(1.0d0-2.0d0*alpha_min))THEN
              IF(smac.eq.3) H_il_i=0.0d0
              IF(cell%tg(i).ge.cell%ts(i))THEN
                 H_ig_i=0.0d0
                 H_gf_i=0.0d0
              ENDIF
           ENDIF
!               
!..........Update H_il,H_ig,H_gf
!        
           H_il(i)=H_il_i
           H_ig(i)=H_ig_i
           H_gf(i)=H_gf_i
         ENDDO
!
      ENDIF
!
      RETURN
      END SUBROUTINE udfn_calc_HTC_int_i
!
      SUBROUTINE udfn_weight1(weight,x,x2,x1)
!
!     User-defined weight for heat transfer coefficients
!     (only when "i_weight" is used)
!
      IMPLICIT NONE
!
      REAL(8) weight,x,x2,x1
!
      weight=DMIN1(1.0d0,DMAX1(0.0d0,DEXP(-8.0d0*(x-x1)/(x2-x1))))
!      
      RETURN
      END SUBROUTINE udfn_weight1
!
      SUBROUTINE udfn_flashing_hif1(i,H_il_i,H_ig_i,ag)
!
!     User-defined H_il,H_ig for flashing
!        
      USE VOL_DATA             
      USE Zmodel   , ONLY: h_il_min
      USE Zqvol    , ONLY: h_il
!
      IMPLICIT NONE
!
      INTEGER i
      REAL(8) H_il_i,H_ig_i,coeff,ag
!
      IF(ag.gt.0.999d0)THEN
         IF(cell%tl(i).gt.cell%ts(i)) H_il_i=DMAX1(H_il_i,1.0e6)
      ENDIF
      IF(ag.lt.0.001d0)THEN
         IF(cell%tg(i).lt.cell%ts(i)) H_ig_i=DMAX1(H_ig_i,1.0e6)
      ENDIF
      coeff=0.05d0
      IF(cell%regime(i).eq.11.or.cell%regime(i).eq.12)THEN
         H_il_i=DMAX1(H_il_i,cell%alphal(i)*cell%rhol(i)*cell%cpl(i)/coeff) 
         H_il_i=DMIN1(1.0e8,DMAX1(H_il_i,H_il_min))
         H_il_i=0.99d0*H_il(i)+0.01d0*H_il_i
      ENDIF  
!
      RETURN
      END SUBROUTINE udfn_flashing_hif1
