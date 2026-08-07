!
      SUBROUTINE rv_hts_2d
!
      USE Zinterface
      USE Zrv_mpi         , ONLY: ncell_fuel_rod_p
      USE Zcore           , ONLY: np
      USE Zconst2         , ONLY: dt
      USE Zrv_hts_2d      , ONLY: nr_2d,nz0_2d,nmat_2d,t_fuel,vl_2f,vr_2f,                    &
                                  qvol_2f,arl_2f,arr_2f,azl_2f,azr_2f,arw_2f,                 &
                                  hlr_2f,hgr_2f,hstr_2f,hspr_2f,tlr_2f,tgr_2f,tstr_2f,tspr_2f
      USE Zrv_model       , ONLY: rv_gapcond
      USE Zrv_ncell       , ONLY: ncell_fuel_rod,nz_fine,nz_fuel_rod,neigh_fuel_rod,nrod_fuel_rod
      USE Zrv_gap_cond    , ONLY: nr_gapi,cond_gap,h_gap,width_gap,block_gap,rhocp_gap
      
      USE Zrv_ncell       , ONLY: cupid_cell_hts2d,nz_fine !pik-radiation_component
      USE Zrad_comp       , ONLY: qrad_rod !pik-radiation_component
!
      USE Zrv_hts_2d,      ONLY:rod_azr_2f,rod_arr_2f,rod_vr_2f,rod_azl_2f, & 
                                rod_arl_2f,rod_vl_2f,rod_arw_2f,ri_2d_opt
!      
      IMPLICIT NONE
!
!.....Local variables
      INTEGER No_r
      INTEGER :: nz,n1,n2,i,j,k
      INTEGER :: m
      REAL(8) :: al,ar,at,ab
      REAL(8) :: g,q
!.....Local arrays
      REAL(8),DIMENSION(ncell_fuel_rod) :: azr_2f1,arr_2f1,vr_2f1,azl_2f1,arl_2f1,vl_2f1,arw_2f1
      REAL(8) :: cond_tr(ncell_fuel_rod,nr_2d-1),cond_bl(ncell_fuel_rod,nr_2d-1)
      REAL(8) :: a(ncell_fuel_rod,nr_2d),b(ncell_fuel_rod,nr_2d),c(ncell_fuel_rod,nr_2d),d(ncell_fuel_rod,nr_2d)
      REAL(8) :: rcp(ncell_fuel_rod,nr_2d),cond(ncell_fuel_rod_p,nr_2d)
!.....Gap conductance
      LOGICAL,SAVE:: initial_gap=.true.
!
      rcp(:,:)=0.d0
      cond(:,:)=0.d0
!
      CALL rv_mat_prop_2d(nmat_2d,t_fuel,rcp,cond)
!
!.....Gap Conductance Model: LSJ
!
      IF(rv_gapcond) THEN
         IF(initial_gap) THEN
            ALLOCATE(cond_gap(ncell_fuel_rod),h_gap(ncell_fuel_rod),width_gap(ncell_fuel_rod),block_gap(ncell_fuel_rod))
            block_gap=0.d0
            width_gap=100.d0
            initial_gap=.false.
         ENDIF   
!
         CALL rv_gap_conductance
!         
         DO k=1,ncell_fuel_rod
            cond(k,nr_gapi)=cond_gap(k)
            rcp(k,nr_gapi)=rhocp_gap !rv_gapcond.in
         ENDDO
!         
      ENDIF   
!      
!
      IF(ncell_fuel_rod_p.gt.ncell_fuel_rod) THEN
         IF(np.gt.1) CALL communicate_rv_2d(t_fuel, &
                                            cond)
      ENDIF
!
      nz=nz0_2d*nz_fine  
!
!.....Left wall: adiabatic boundary
      i=1
!
!.....geometrical data
      IF(ri_2d_opt.eq.1)THEN
         DO k=1,ncell_fuel_rod
            j=nz_fuel_rod(k)
            No_r=nrod_fuel_rod(k)
            azr_2f1(k)=rod_azr_2f(j,No_r,i)
            arr_2f1(k)=rod_arr_2f(j,No_r,i)
            vr_2f1 (k)=rod_vr_2f (j,No_r,i)
         ENDDO
      ELSE
         DO k=1,ncell_fuel_rod
            j=nz_fuel_rod(k)
            azr_2f1(k)=azr_2f(j,i)
            arr_2f1(k)=arr_2f(j,i)
            vr_2f1 (k)=vr_2f (j,i)
         ENDDO
      ENDIF
!         
      DO k=1,ncell_fuel_rod
         j=nz_fuel_rod(k)
         n1=neigh_fuel_rod(k,1)
         n2=neigh_fuel_rod(k,2)
!
         g=vr_2f1(k)*rcp(k,i)
         q=vr_2f1(k)*qvol_2f(k,i)
         ar=cond(k,i)*arr_2f1(k)
         a(k,i)=0.d0
         c(k,i)=ar*dt
         b(k,i)=-c(k,i)-g
         IF(j.eq.1)THEN
            cond_tr(k,i)=2.d0*cond(k,i)*cond(n2,i)/(cond(n2,i)+cond(k ,i))
            at=cond_tr(k,i)*azr_2f1(k)
            d(k,i)=(at*dt-g)*t_fuel(k,i)-(at*t_fuel(n2,i)+q)*dt
         ELSEIF(j.eq.nz)THEN
            cond_bl(k,i)=2.d0*cond(k,i)*cond(n1,i)/(cond(k ,i)+cond(n1,i))
            ab=cond_bl(k,i)*azr_2f1(k)
            d(k,i)=(ab*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+q)*dt
         ELSE
            cond_tr(k,i)=2.d0*cond(k,i)*cond(n2,i)/(cond(n2,i)+cond(k ,i))
            cond_bl(k,i)=2.d0*cond(k,i)*cond(n1,i)/(cond(k ,i)+cond(n1,i))
            at=cond_tr(k,i)*azr_2f1(k)
            ab=cond_bl(k,i)*azr_2f1(k)
            d(k,i)=((ab+at)*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+at*t_fuel(n2,i)+q)*dt
         ENDIF
      ENDDO !k
!
      DO i=2,nr_2d-1
!
!........geometrical data
         IF(ri_2d_opt.eq.1)THEN
            DO k=1,ncell_fuel_rod
               j=nz_fuel_rod(k)
               No_r=nrod_fuel_rod(k)
               azr_2f1(k)=rod_azr_2f(j,No_r,i)
               arr_2f1(k)=rod_arr_2f(j,No_r,i)
               vr_2f1 (k)=rod_vr_2f (j,No_r,i)
               azl_2f1(k)=rod_azl_2f(j,No_r,i)
               arl_2f1(k)=rod_arl_2f(j,No_r,i)
               vl_2f1 (k)=rod_vl_2f (j,No_r,i)
            ENDDO
         ELSE
            DO k=1,ncell_fuel_rod
               j=nz_fuel_rod(k)
               azr_2f1(k)=azr_2f(j,i)
               arr_2f1(k)=arr_2f(j,i)
               vr_2f1 (k)=vr_2f (j,i)
               azl_2f1(k)=azl_2f(j,i)
               arl_2f1(k)=arl_2f(j,i)
               vl_2f1 (k)=vl_2f (j,i)
            ENDDO
         ENDIF
!         
         DO k=1,ncell_fuel_rod
            j=nz_fuel_rod(k)
            n1=neigh_fuel_rod(k,1)
            n2=neigh_fuel_rod(k,2)
!
            g=vl_2f1(k)*rcp(k,i-1)+vr_2f1(k)*rcp(k,i)
            q=vl_2f1(k)*qvol_2f(k,i-1)+vr_2f1(k)*qvol_2f(k,i)                 !Volumetric heat source
            al=cond(k,i-1)*arl_2f1(k)
            ar=cond(k,i  )*arr_2f1(k)
            a(k,i)=al*dt
            c(k,i)=ar*dt
            b(k,i)=-a(k,i)-c(k,i)-g
            IF(j.eq.1)THEN
               cond_tr(k,i)=2.d0*cond(k,i)*cond(n2,i)/(cond(n2,i)+cond(k ,i))
               at=cond_tr(k,i-1)*azl_2f1(k)+cond_tr(k,i)*azr_2f1(k)
               d(k,i)=(at*dt-g)*t_fuel(k,i)-(at*t_fuel(n2,i)+q)*dt
            ELSEIF(j.eq.nz)THEN
               cond_bl(k,i)=2.d0*cond(k,i)*cond(n1,i)/(cond(k ,i)+cond(n1,i))
               ab=cond_bl(k,i-1)*azl_2f1(k)+cond_bl(k,i)*azr_2f1(k)
               d(k,i)=(ab*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+q)*dt
            ELSE
               cond_tr(k,i)=2.d0*cond(k,i)*cond(n2,i)/(cond(n2,i)+cond(k ,i))
               cond_bl(k,i)=2.d0*cond(k,i)*cond(n1,i)/(cond(k ,i)+cond(n1,i))
               at=cond_tr(k,i-1)*azl_2f1(k)+cond_tr(k,i)*azr_2f1(k)
               ab=cond_bl(k,i-1)*azl_2f1(k)+cond_bl(k,i)*azr_2f1(k)
               d(k,i)=((ab+at)*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+at*t_fuel(n2,i)+q)*dt
            ENDIF
         ENDDO !k
      ENDDO !i
!
!.....Right wall: hf+tf boundary
      i=nr_2d
!
!.....geometrical data
      IF(ri_2d_opt.eq.1)THEN
         DO k=1,ncell_fuel_rod
            j=nz_fuel_rod(k)
            No_r=nrod_fuel_rod(k)
            azl_2f1(k)=rod_azl_2f(j,No_r,i)
            arl_2f1(k)=rod_arl_2f(j,No_r,i)
            vl_2f1 (k)=rod_vl_2f (j,No_r,i)
            arw_2f1(k)=rod_arw_2f(j,No_r)
         ENDDO
      ELSE
         DO k=1,ncell_fuel_rod
            j=nz_fuel_rod(k)
            azl_2f1(k)=azl_2f(j,i)
            arl_2f1(k)=arl_2f(j,i)
            vl_2f1 (k)=vl_2f (j,i)
            arw_2f1(k)=arw_2f(j)
         ENDDO
      ENDIF
!         
      DO k=1,ncell_fuel_rod
         m=cupid_cell_hts2d(k) !pik-radiation_component
         j=nz_fuel_rod(k)
         n1=neigh_fuel_rod(k,1)
         n2=neigh_fuel_rod(k,2)
!
         g=vl_2f1(k)*rcp(k,i-1)
         q=vl_2f1(k)*qvol_2f(k,i-1)
         al=cond(k,i-1)*arl_2f1(k)
         a(k,i)=al*dt
         c(k,i)=0.d0
         b(k,i)=-a(k,i)-g-(hlr_2f(k)+hgr_2f(k)+hstr_2f(k)+hspr_2f(k))*arw_2f1(k)*dt
         IF(j.eq.1)THEN
            at=cond_tr(k,i-1)*azl_2f1(k)
            d(k,i)=(at*dt-g)*t_fuel(k,i)-(at*t_fuel(n2,i)+q+                                         &
                     (hlr_2f(k)*tlr_2f(k)+                                                           &
                     hgr_2f(k)*tgr_2f(k)+hstr_2f(k)*tstr_2f(k)+hspr_2f(k)*tspr_2f(k))*arw_2f1(k))*dt &
                     -qrad_rod(m)/nz_fine*dt !pik-radiation_component
         ELSEIF(j.eq.nz)THEN
            ab=cond_bl(k,i-1)*azl_2f1(k)
            d(k,i)=(ab*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+q+                                                &
                     (hlr_2f(k)*tlr_2f(k)+hgr_2f(k)*tgr_2f(k)+hstr_2f(k)*tstr_2f(k)+hspr_2f(k)*tspr_2f(k))* &
                     arw_2f1(k))*dt                                                                         &
                     -qrad_rod(m)/nz_fine*dt !pik-radiation_component
         ELSE
            at=cond_tr(k,i-1)*azl_2f1(k)
            ab=cond_bl(k,i-1)*azl_2f1(k)
            d(k,i)=((ab+at)*dt-g)*t_fuel(k,i)-(ab*t_fuel(n1,i)+at*t_fuel(n2,i)+q+                           &
                     (hlr_2f(k)*tlr_2f(k)+hgr_2f(k)*tgr_2f(k)+hstr_2f(k)*tstr_2f(k)+hspr_2f(k)*tspr_2f(k))* &
                     arw_2f1(k))*dt                                                                         &
                     -qrad_rod(m)/nz_fine*dt !pik-radiation_component
               
         ENDIF
      ENDDO !k
!
!.....Solve tridiagonal matrix vector mode
!
      CALL rv_tdiag_2d(a,b,c,d,ncell_fuel_rod,nr_2d)
!
      DO i=1,nr_2d
         DO k=1,ncell_fuel_rod
            t_fuel(k,i)=d(k,i)
         ENDDO
      ENDDO
!
      END SUBROUTINE rv_hts_2d
