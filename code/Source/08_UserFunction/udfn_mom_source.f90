!
      SUBROUTINE udfn_mom_source
!
!     Modifies the momentum source terms at the free surface cells
!
      USE VOL_DATA    , ONLY: cell         
      USE Zzone       , ONLY: ncell_fluid
      USE Zparam      , ONLY: ndim,pi
      USE Zconst1     , ONLY: vv_prob
      USE Zcoord1     , ONLY: xloc
      USE Zcoord2     , ONLY: cell_leng
      USE Zcoord3     , ONLY: floss,porosity,vol
      USE Zvector     , ONLY: ul_o,ug_o,vl_o,vg_o
      USE Zm_src      , ONLY: src_gas,src_liq
      USE Zporous     , ONLY: fric_model_liq      
      USE KSMR        , ONLY: zone_comp,vol_sg
      USE Zconst2     , ONLY: hydraulicd !icarus2002
!
      IMPLICIT NONE
!            
!.....Local variables
      INTEGER i,ix
!      
      REAL(8) rusquarg,rusquarl,velo,velol,velog,velo_tmp
      REAL(8) velol_tmp,velog_tmp,hd,kloss_grid
      REAL(8) cell_height !icarus2002
!      
!     KSMR-SG-porous      
      INTEGER fric_model,TPM_case,de_model,iter_newton,fric_scheme
      REAL(8) fch(ndim),fric(ndim),cah,cvh,ash,TPM,Re,visc,pv,de,d0,ntube,cleng_avg,c1,c2,c3,c4, &
             ksi,gch,fmch,mum,gg,fma,ff,fd,fcoeff,nn,retr,pvd0   
!      
      IF(vv_prob.eq.'KSMR-SG-porous') THEN
!
!        Explicit friction/Implicit friction (1=explicit anisotropic, 2=implicit anisotropic friction)
         fric_scheme=2                                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input          
!
!        Two-phase multiplier
         TPM_case=1 !0=none, 1=default, 2=armand               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input
         SELECT CASE(TPM_case)
         CASE(1)
            TPM=1.d0 
         CASE(2)
            TPM=cell%rhol(i)/cell%rhom(i)
         CASE(3)
            IF(cell%alphag(i).le.0.6d0) THEN
               TPM=(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**1.42
            ELSEIF(cell%alphag(i).le.0.9d0) THEN
               TPM=0.478d0*(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**2.2
            ELSE
               TPM=1.73d0*(1.d0-cell%quals(i))*(1.d0-cell%quals(i))/(1.d0-cell%alphag(i))**1.64
            ENDIF
         END SELECT         
!
!        Friction model (1=model1, 2=model2, 3=model3)
         fric_model=2                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input 
!
!        Equivalent diameter model (1=triangular pitch, 2=square pitch)
         de_model=2                                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input
         pv=0.024d0  ! radial pitch                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input
         d0=0.015d0  ! tube diameter                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input
         ntube=3952.d0 !total number of horizontal tubes       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! User-input
         IF(de_model.eq.1) THEN !triangular pitch (ICTYPE=4) 
            de=1.10266*pv*pv/d0-d0 
            c1=1.d0
            c2=1.155d0
            c3=4.619d0  !direct cross flow in triangualr area
            c4=pv
            ksi=c1*(1.d0-c2*pi*0.25d0*(d0/pv)**2)/(1.d0-d0/pv)
            gch=c3*ksi*ksi*dsqrt(3.d0)*pv/4.d0/pi/d0
         ELSE !square pitch (ICTYPE=1)
            de=1.27324*pv*pv/d0-d0
            c1=1.d0
            c2=1.d0            
            c3=4.d0     !direct cross flow in square area
            ksi=c1*(1.d0-c2*pi*0.25d0*(d0/pv)**2)/(1.d0-d0/pv)
            gch=c3*ksi*ksi*dsqrt(3.d0)*pv/2.d0/pi/d0
         ENDIF         
!         
         DO i=1,ncell_fluid
            IF(zone_comp(i).ne.10) CYCLE 
!            
!           absolute mixed velocity
            velol=0.d0 
            velog=0.d0 
            DO ix=1,ndim 
               velol=velol+vl_o(i,ix)*vl_o(i,ix)
               velog=velog+vl_o(i,ix)*vl_o(i,ix)
            ENDDO
            velo=cell%alphal(i)*velol+cell%alphag(i)*velog
!
!           mixed viscosity
            visc=cell%alphag(i)*cell%lviscosg(i)+(1.d0-cell%alphag(i))*cell%lviscosl(i)
!
            cah=4.d0*porosity(i)*pv/pi/d0
            cvh=(porosity(i)*pv/(pv-d0))**2
            
            cleng_avg=0.d0
            DO ix=1,ndim 
               cleng_avg=cleng_avg+cell_leng(i,ix)
            ENDDO            
            cleng_avg=cleng_avg/ndim
           
            ash=ntube*pi*d0*cleng_avg           ! total surface area of horizontal tubes in a cell
            ash=vol(i)/vol_sg*4109.6            ! total surface area of horizontal tubes in a cell    
!            
!           friction model            
            DO ix=1,ndim
               IF(fric_model.eq.1) THEN 
                  Re=DMAX1(1.d0,cell%rhom(i)*velo*de/visc)
                  fch(ix)=0.1957d0*Re**(-0.1853)
               ELSE
                  mum=cell%quals(i)/cell%lviscosg(i)+(1.d0-cell%quals(i))/cell%lviscosl(i)
                  gg=cell%rhom(i)*DABS(vl_o(i,ix))
                  RE=DMAX1(1.d0,1.d0*de*gg/mum)                    
                  IF(ix.ne.ndim) THEN 
                     fch(ix)=0.1957d0*Re**(-0.1853)
                  ELSE
                     IF(fric_model.eq.2) THEN 
!                       Newton-Rhapson Method                        
                        fMa=0.5d0 !0.5d0 !0.1d0
                        DO iter_newton=1,100
                           ff=1.d0/fMa+0.868589d0*DLOG((0.1524d-5)/3.7d0/de+2.51d0/RE/fMa)
                           fd=-1.d0/fMa**2-0.868589d0*2.51d0/(RE*(0.1524d-5)/3.7d0/de*fMa**2 + 2.51d0*fMa)
                           fMa = fMa-ff/fd
                           IF(ff/fd.le.1.d-4) THEN
                               EXIT
                           ELSE
                               print*,'error',ff/fd
                           ENDIF 
                        ENDDO
                        fmch=fMa**2
                     ELSEIF(fric_model.eq.3) THEN   
                        pvd0=pv/d0 
!                        
!                       Newton-Rhapson Method                        
                        fMa=0.27d0 !0.4, 0.35
                        DO iter_newton=1,100
                           ff=fMa+6.947812d0*DLOG((0.1524d-5)/3.7d0/de+0.314d0/fMa)
                           fd=1.d0-6.947812d0*0.314d0/((0.1524d-5)/3.7d0/de*fMa**2 + 0.314d0*fMa)
                           fMa = fMa-ff/fd
                           IF(ff/fd.le.1.d-4) THEN
                               EXIT
                           ELSE
                               print*,'error',ff/fd
                           ENDIF    
                        ENDDO
                        retr=fMa**2
                        IF(RE.le.retr) THEN
                           fmch=70.d0/RE*(1.d0/pvd0)**1.6  
                        ELSE
                           IF(de_model.eq.1) THEN !ICTYPE=4 (triangular)
                             fcoeff=0.25+0.1175/(c4/d0-1.d0)**1.08  !f2
                             nn=-0.16 
                           ELSEIF(de_model.eq.2) THEN !ICTYPE=1 (square)
                             fcoeff=0.044+0.08*pvd0/(pvd0-1.d0)**(0.43+1.13/pvd0)  !f1
                             nn=-0.15
                           ENDIF
                           fmch=fcoeff*(pvd0*RE)**nn
                        ENDIF
                     ENDIF
                     fch(ix)=gch*fmch 
                  ENDIF
               ENDIF
               fric(ix)=fch(ix)*cah*cvh*ash*cell%rhol(i)*TPM
               
               IF(fric_scheme.eq.1) THEN  !explicit
                   src_liq(i,ix)=src_liq(i,ix)-fric(ix)*velo*vl_o(i,ix)
!                   src_gas(i,ix)=src_liq(i,ix)-fric(ix)*velo*vg_o(i,ix)
               ELSEIF(fric_scheme.eq.2) THEN !implicit
                   fric_model_liq(i,ix)=-fric(ix)*velo !dABS(vl_o(i,ix)) !velo
!                   fric_model_gas(i,ix)=-fric(ix)*velo
               ENDIF
            ENDDO   
!            
         ENDDO        
      ENDIF        
!           
!........CUPID-RV
!      
         IF(vv_prob.eq.'fs_31203'    .or.&
            vv_prob.eq.'fs_31302'    .or.&            
            vv_prob.eq.'fs_31701'    .or.&            
            vv_prob.eq.'fs_31805'    )THEN 
            DO i=1,ncell_fluid
!               rusquarg=cell%alphag(i)*cell%rhog(i)*ug_o(i)**2
!               rusquarl=cell%alphal(i)*cell%rhol(i)*ul_o(i)**2
               velo_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
               IF(ndim.eq.3) velo_tmp=velo_tmp+vl_o(i,3)*vl_o(i,3)
               velo=DSQRT(velo_tmp)            
               DO ix=1,ndim
!                  src_gas(i,ix)=src_gas(i,ix)-floss(i,ix)*rusquarg
!                  src_liq(i,ix)=src_liq(i,ix)-floss(i,ix)*rusquarl
                  src_liq(i,ix)=src_liq(i,ix)-(floss(i,ix)/0.182828)*cell%rhol(i)*velo*vl_o(i,ix)*cell%alphal(i)
               ENDDO
            ENDDO      
         ENDIF     
!
         IF(vv_prob.eq.'fs_31203_3D'    .or.&
            vv_prob.eq.'fs_31302_3D'    .or.&            
            vv_prob.eq.'fs_31701_3D'    .or.&            
            vv_prob.eq.'fs_31805_3D'    )THEN 
            DO i=1,ncell_fluid
               velo_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
               IF(ndim.eq.3) velo_tmp=velo_tmp+vl_o(i,3)*vl_o(i,3)
               velo=DSQRT(velo_tmp)            
               DO ix=1,ndim
                  src_liq(i,ix)=src_liq(i,ix)-(floss(i,ix)/0.182828)*cell%rhol(i)*velo*vl_o(i,ix)*cell%alphal(i)
               ENDDO
!
               velo_tmp=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)
               IF(ndim.eq.3) velo_tmp=velo_tmp+vg_o(i,3)*vg_o(i,3)
               velo=DSQRT(velo_tmp)            
               DO ix=1,ndim
                  src_gas(i,ix)=src_gas(i,ix)-(floss(i,ix)/0.182828)*cell%rhog(i)*velo*vg_o(i,ix)*cell%alphag(i)
               ENDDO            
            ENDDO      
         ENDIF           
!
         IF(vv_prob.eq.'rbht1196_3d'.or.vv_prob.eq.'rbht1196_1d')THEN 
            DO i=1,ncell_fluid
               velo_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)
               IF(ndim.eq.3) velo_tmp=velo_tmp+vl_o(i,3)*vl_o(i,3)
               velol=DSQRT(velo_tmp)            
!
               velo_tmp=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)
               IF(ndim.eq.3) velo_tmp=velo_tmp+vg_o(i,3)*vg_o(i,3)
               velog=DSQRT(velo_tmp)                        
               DO ix=1,ndim
                  src_liq(i,ix)=src_liq(i,ix)-(floss(i,ix)/cell_leng(i,ndim)/2.d0)*cell%rhol(i)*velol*vl_o(i,ix)*cell%alphal(i)
                  src_gas(i,ix)=src_gas(i,ix)-(floss(i,ix)/cell_leng(i,ndim)/2.d0)*cell%rhog(i)*velog*vg_o(i,ix)*cell%alphag(i)
               ENDDO
            ENDDO      
         ENDIF  
!     
!........rocom, rocom_mc, rocom_2d
!      
         IF(vv_prob.eq.'rocom' .or. vv_prob.eq.'rocom_mc' .or. vv_prob.eq.'rocom_2d')THEN 
            DO i=1,ncell_fluid
               rusquarg=cell%alphag(i)*cell%rhog(i)*ug_o(i)**2
               rusquarl=cell%alphal(i)*cell%rhol(i)*ul_o(i)**2
               DO ix=1,ndim
                  src_gas(i,ix)=src_gas(i,ix)-floss(i,ix)*rusquarg
                  src_liq(i,ix)=src_liq(i,ix)-floss(i,ix)*rusquarl
               ENDDO
            ENDDO      
         ENDIF
!
!........sgp
!
         IF(vv_prob.eq.'sgp_separator')THEN
            CALL udfn_sg_mom_wall
            DO i=1,ncell_fluid
               src_gas(i,1)=src_gas(i,1)-cell%vfwg_x(i)
               src_gas(i,2)=src_gas(i,2)-cell%vfwg_y(i)
               src_gas(i,3)=src_gas(i,3)-cell%vfwg_z(i)
               src_liq(i,1)=src_liq(i,1)-cell%vfwl_x(i)
               src_liq(i,2)=src_liq(i,2)-cell%vfwl_y(i)
               src_liq(i,3)=src_liq(i,3)-cell%vfwl_z(i)
            ENDDO
         ENDIF
!
!........icarus2002
!
         IF(vv_prob.eq.'icarus2002')THEN
            DO i=1,ncell_fluid
!              
!              absolute velocity 
               velol_tmp=vl_o(i,1)*vl_o(i,1)+vl_o(i,2)*vl_o(i,2)+vl_o(i,3)*vl_o(i,3)
               velog_tmp=vg_o(i,1)*vg_o(i,1)+vg_o(i,2)*vg_o(i,2)+vg_o(i,3)*vg_o(i,3)
               velol=DSQRT(velol_tmp)
               velog=DSQRT(velog_tmp)
               velo=cell%alphag(i)*velog+(1.d0-cell%alphag(i))*velol
!
!               MATRA friction model in axial direction
               ix=3 
               hd=hydraulicd(i)
               TPM=5.d0
!
! MARTA friction (rv_fric_w=1을 사용해서 현재의 MATRA 모델은 사용하지 않고 오직 spacer grid 모델만 사용함)              
!
               !Reg=DMAX1(1.0d0,(cell%rhog(i)*velog*(1.0d0*hd)/cell%lviscosg(i)))
               !CALL mom_wall_matra(Reg,fric(ix))
               !src_gas(i,ix)=src_gas(i,ix)-TPM*fric(ix)/(2.0d0*hd)*cell%rhog(i)*velog*vg_o(i,ix)*cell%alphag(i)               
!
!              friction model for spacer grid               
               IF(xloc(i,3).ge.0.3d0.and.xloc(i,3).le.0.325d0)THEN  
                  kloss_grid=7.0d0
                  cell_height=0.025d0
                  src_gas(i,ix)=src_gas(i,ix)-(kloss_grid/(2.0d0*cell_height))*cell%rhog(i)*velo*vg_o(i,ix)
               ENDIF               
!               
            ENDDO
         ENDIF         
!
      RETURN
      END SUBROUTINE udfn_mom_source
