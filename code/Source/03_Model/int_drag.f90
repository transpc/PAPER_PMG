!
      SUBROUTINE int_drag
!
!     This routine calculates interfacial drag coefficients using minimum value,
!     vFgl_1_min=1.0d0,vFgl_2_min=1.0d0,vFgl_3_min=1.0d2 for bubbly,mist,sharp.
!
      USE VOL_DATA                 
      USE Zzone        , ONLY: ncell_fluid
      USE Zconst1      , ONLY: mdrag
      USE Zdaint_ag    , ONLY: daint1_ag_bc,daint2_ag_cm,aint_bc,aint_cm
      USE Zflowregime  , ONLY: alphag_bc,alphag_cm,gamma_1,gamma_2,vFgl_1_min,vFgl_2_min,vFgl_3_min
      USE Zmodel       , ONLY: drag_coeff_a,drag_coeff_b,cb_bubble
      USE Zmpi         , ONLY: ncell_fp      
      USE Ztimecon     , ONLY: alpha_min
      USE Zuserdefined , ONLY: udfl_mom_drag_i
      USE Zmodel       , ONLY: i_weight
      USE Zvector      , ONLY: vrel_o,vg_o,vl_o
      USE Zvoid        , ONLY: dagdx,gamma_void
      
!
      IMPLICIT NONE
!
      INTEGER i,n
!      
      REAL(8) vfgl_i(ncell_fp),vfgd_i(ncell_fp)
      REAL(8) drag_coeff_cb1(ncell_fp)
      REAL(8) Re,Cb,cosine,Ci,Ci_tan,Ci_ort,vFgl_1,vFgl_2,vFgl_3,weight
      REAL(8) vFgl_1c,vFgl_2c,vFgl_1_bc,vFgl_2_cm,dFdag_1_bc,dFdag_2_cm,dFdag
!
      n=ncell_fluid
!
      IF(mdrag.eq.0)THEN
!
!........Use of linear equation using input
!      
       DO i=1,n 
         vFgl_i(i)=Drag_coeff_a+Drag_coeff_b*DMAX1(1.0e-2,cell%alphag(i)*cell%alphal(i))
         drag_coeff_Cb1(i)=0.44d0
         vFgd_i(i)=Drag_coeff_a+Drag_coeff_b*DMAX1(1.0e-5,cell%alphag(i)*cell%alphal(i))      
       ENDDO
      ELSEIF(mdrag.eq.-1)THEN
!
!........Use of constant value using input
!        
      DO i=1,n 
         vFgl_i(i)=Drag_coeff_a
         drag_coeff_Cb1(i)=0.44d0
         vFgd_i(i)=Drag_coeff_a+Drag_coeff_b*DMAX1(1.0e-5,cell%alphag(i)*cell%alphal(i))      
       ENDDO
      ELSE
!
      DO i=1,n 
!
!........Bubbly Flow
!
         Re=cell%rhol(i)*vrel_o(i)*cell%D1(i)/cell%lviscosl(i)
         IF(Re.le.1000.0d0)THEN
            Cb=DMAX1(0.44d0,24.0d0/DMAX1(Re,alpha_min)*(1.0d0+0.15d0*Re**0.687d0))
         ELSE
            Cb=0.44d0
         ENDIF
         vFgl_1c=1.0d0/8.0d0*cell%rhol(i)*Cb*vrel_o(i)
         vFgl_1=vFgl_1c*cell%aint1(i)
!
!........Mist Flow
!
         Re=cell%rhog(i)*vrel_o(i)*cell%D2(i)/cell%lviscosg(i)
         IF(Re.le.1000.0d0)THEN
            Cb=DMAX1(0.44d0,24.0d0/DMAX1(Re,alpha_min)*(1.0d0+0.15d0*Re**0.687d0))
         ELSE
            Cb=0.44d0
         ENDIF
         vFgl_2c=1.0d0/8.0d0*cell%rhog(i)*Cb*vrel_o(i)
         vFgl_2=vFgl_2c*cell%aint2(i)
!
!........Sharp Interface
!
         cosine=DOT_PRODUCT(vl_o(i,:)-vg_o(i,:),dagdx(i,:))/DMAX1(vrel_o(i),alpha_min) &
                /DMAX1(alpha_min,DSQRT(DOT_PRODUCT(dagdx(i,:),dagdx(i,:))))
         cosine=DMIN1(1.0d0,DABS(cosine))
         Ci_tan=0.005d0
         Ci_ort=1.0d0
         Ci=Ci_tan+(Ci_ort-Ci_tan)*cosine
         vFgl_3=1.0d0/2.0d0*cell%aint3(i)*cell%rhog(i)*Ci*vrel_o(i)
!      
         vFgl_1=DMAX1(vFgl_1,vFgl_1_min)
         vFgl_2=DMAX1(vFgl_2,vFgl_2_min)
         vFgl_3=DMAX1(vFgl_3,vFgl_3_min)
!
!........Calculates drag coefficient according to the topology map
!        
         IF(cell%regime(i).eq.11)THEN
            vFgl_i(i)=vFgl_1
         ELSEIF(cell%regime(i).eq.13)THEN
            vFgl_i(i)=vFgl_2
         ELSEIF(cell%regime(i).eq.3)THEN
            vFgl_i(i)=vFgl_3
         ELSEIF(cell%regime(i).eq.12)THEN
            IF(i_weight.ge.1)THEN               
!
!..............Apply the weight factor by udf
!
               CALL udfn_weight(weight,cell%alphag(i),alphag_cm,alphag_bc)
               vFgl_i(i)=weight*vFgl_1+(1.0d0-weight)*vFgl_2
!               
            ELSE
               vFgl_1_bc=vFgl_1c*aint_bc(i)
               vFgl_2_cm=vFgl_2c*aint_cm(i)
               dFdag_1_bc=vFgl_1c*daint1_ag_bc(i)
               dfdag_2_cm=vFgl_2c*daint2_ag_cm(i)
!
!..............Apply cubic interpolation
!               
               CALL cub_interp(cell%alphag(i),alphag_bc,alphag_cm,vFgl_1_bc,vFgl_2_cm,dFdag_1_bc,dfdag_2_cm,vFgl_i(i),dFdag)
            ENDIF
         ELSEIF(cell%regime(i).eq.21)THEN
!
!...........Apply the weight factor by udf
!         
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            IF(i_weight.ge.1) CALL udfn_weight(weight,gamma_void(i),gamma_2,gamma_1)           
            vFgl_i(i)=weight*vFgl_1+(1.0d0-weight)*vFgl_3
!
!            IF(udfl_mom_drag_increase) vFgl_i(i)=DMAX1(1.0d0,vFgl_i(i))*20.d0
         ELSEIF(cell%regime(i).eq.23)THEN
!
!...........Apply the weight factor by udf
!          
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            IF(i_weight.ge.1) CALL udfn_weight(weight,gamma_void(i),gamma_2,gamma_1)
            vFgl_i(i)=weight*vFgl_2+(1.0d0-weight)*vFgl_3
!            if(udfl_mom_drag_increase) vFgl_i(i)=DMAX1(1.0d0,vFgl_i(i))*20.d0
         ELSEIF(cell%regime(i).eq.22)THEN
!
!...........Apply the weight factor by udf
!          
            weight=(alphag_cm-cell%alphag(i))/(alphag_cm-alphag_bc)
            IF(i_weight.ge.1) CALL udfn_weight(weight,cell%alphag(i),alphag_cm,alphag_bc)
            vFgl_i(i)=weight*vFgl_1+(1.0d0-weight)*vFgl_2
            weight=(gamma_2-gamma_void(i))/(gamma_2-gamma_1)
            IF(i_weight.ge.1) CALL udfn_weight(weight,gamma_void(i),gamma_2,gamma_1)
            vFgl_i(i)=weight*vFgl_i(i)+(1.0d0-weight)*vFgl_3
!            
!            IF(udfl_mom_drag_increase) vFgl_i(i)=DMAX1(1.0d0,vFgl_i(i))*20.d0
         ENDIF
         drag_coeff_Cb1(i)=Cb       
         vFgd_i(i)=Drag_coeff_a+Drag_coeff_b*DMAX1(1.0e-5,cell%alphag(i)*cell%alphal(i))      
       ENDDO
      ENDIF
!
!      vFgd_i(:)=Drag_coeff_a+Drag_coeff_b*DMAX1(1.0e-5,cell%alphag(:)*cell%alphal(:))      
!
!.....Linear drag coefficient between gas and droplets. 
!
!
!.....User-defined drag coefficient between gas and liquid. 
!
      IF(udfl_mom_drag_i) CALL udfn_mom_drag_i(vfgl_i,n)
!
      DO i=1,n
         cell%vfgl(i)=vfgl_i(i)
         cell%vfgd(i)=vfgd_i(i)
      ENDDO
      DO i=1,n
         Cb_bubble(i)=DMIN1(2.0d0,Drag_coeff_Cb1(i))
      ENDDO
!
      RETURN
      END SUBROUTINE int_drag
!
