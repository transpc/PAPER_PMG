!
      SUBROUTINE udfn_sg_primary_ini
!
      USE Zparam,   ONLY:pi
      USE Zsg,      ONLY:n_group,n_1d,max_1d,vn_1d,p_1d,vin_1d,en_1d,eo_1d,do_tube,     &
                          th_tube,cf,rhoin_1d,visin_1d,dp_primary,pr_flow_rated,tube_length,ar_tube,   &
                          di_tube,pin_1d,ein_1d,pout_1d,h_tube,t_1d,nr_tube,dr_tube,dr_tube2,   &
                          p_tube,vol_tube,t_tube,pr_flow,q_rated,eout_1d,hyd_d,pitch,      &
                          ng_tube,ht_area
!
      IMPLICIT NONE 
!
      INTEGER i,j
      REAL(8) re,pr_flow_tot,dx,d_e
!
      CALL udfn_prop_in
!
      di_tube=do_tube-2.0d0*th_tube
      cf=0.184d0
      re=rhoin_1d*di_tube/visin_1d
      re=dp_primary*2.0d0*di_tube/rhoin_1d*re**0.2d0
!
!.....Initialize inlet veocities of each U-tube group
!
      DO j=1,2
         pr_flow_tot=0.0d0
         DO i=1,n_group
            vin_1d(i)=re/cf/tube_length(i)
            vin_1d(i)=vin_1d(i)**(1.0d0/1.8d0)
            pr_flow(i)=rhoin_1d*ar_tube(i)*vin_1d(i)
            pr_flow_tot=pr_flow_tot+pr_flow(i)
         ENDDO
         IF(j.eq.1) cf=cf*(pr_flow_tot/pr_flow_rated)**1.8
      ENDDO
!
      DO i=1,n_group
         vn_1d(i,:)=vin_1d(i)
      ENDDO 
!
!.....Outlet parameters
!
      pout_1d=pin_1d+dp_primary
      eout_1d=ein_1d-q_rated/pr_flow_tot
      d_e=ein_1d-eout_1d
!
!.....Initialize pressure and internal energy
!
      DO i=1,n_group
         dx=0.0d0
         DO j=1,n_1d(i)
            IF(j.eq.1)THEN
               dx=0.5d0*h_tube(i,j)
            ELSE
               dx=dx+0.5d0*(h_tube(i,j-1)+h_tube(i,j))
            ENDIF
            p_1d(i,j)=pin_1d-dp_primary*dx/tube_length(i)
            en_1d(i,j)=ein_1d-d_e*dx/tube_length(i)
         ENDDO
      ENDDO
!
!.....Initialize fluid properties
!
      CALL udfn_prop_liquid
      eo_1d(:,:)=en_1d(:,:)
!
!.....Set parameters for heat conduction calculation
!
      ALLOCATE(p_tube(nr_tube+1),vol_tube(nr_tube),t_tube(nr_tube,n_group,max_1d))
      p_tube(:)=0.0d0
      vol_tube(:)=0.0d0
      t_tube(:,:,:)=0.0d0
!
      dr_tube=th_tube/nr_tube
      dr_tube2=2.0d0*dr_tube
      hyd_d(:)=4.0d0*(pitch*pitch-pi*do_tube*do_tube/4.0d0)/(pi*do_tube)
!
      p_tube(1)=pi*di_tube
      DO i=2,nr_tube+1
         p_tube(i)=pi*(di_tube+(i-1)*dr_tube2)
      ENDDO
!
      DO i=1,nr_tube
         vol_tube(i)=0.5d0*(p_tube(i)+p_tube(i+1))*dr_tube
      ENDDO
!
!....Heat transfer area
!
      DO i=1,n_group
         DO j=1,n_1d(i)
            ht_area(i,j)=p_tube(nr_tube+1)*h_tube(i,j)*ng_tube(i)
         ENDDO
      ENDDO
!
!.....Initialize utube temperature
!
      DO i=1,n_group
         DO j=1,n_1d(i)
            t_tube(:,i,j)=t_1d(i,j)
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE udfn_sg_primary_ini
