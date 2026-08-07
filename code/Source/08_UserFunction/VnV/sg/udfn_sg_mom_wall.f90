!
      SUBROUTINE udfn_sg_mom_wall
!
!.....Wall drag coefficients for gas and liquid phase.
!
      USE VOL_DATA       ,ONLY: cell
      USE Zzone          ,ONLY: ncell_fluid
      USE Zparam         ,ONLY: pi
      USE Zvector        ,ONLY: vl_o,vg_o,ul_o,ug_o
      USE Zcoord1        ,ONLY: xloc
      USE Zcoord3        ,ONLY: porosity,volp
      USE Zsg            ,ONLY: mult_cell,mult_3d_cell1,mult_3d_cell2,                          &
                                f1_mult,f2_mult,ih,iavb,izp,ihp,hyd_d,ht_area,pitch,do_tube,    &
                                igr,j1d,idc,iavb,izp,ihp,h_tube
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: g,m,i,m1,m2,k
      REAL(8) :: vfwg_x_i,vfwg_y_i,vfwg_z_i,vfwl_x_i,vfwl_y_i,vfwl_z_i
      REAL(8):: a,b,c,rlg,denom,velx,vely,velz,p1,p2,fx,fy,fz,TPM,CV,CA
      REAL(8):: qu,rl,rg,rm,ag,al,pr,visl,visg,vp,poro,ul,ug
      REAL(8):: vxl,vxg,vyl,vyg,vzl,vzg,rex,rey,rez
      REAL(8):: vfwg_i,vfwl_i,arl,arg,cs,afad2
!
!.....Local arrays
      INTEGER :: itube(ncell_fluid)
!
      cell%vfwg_x(:)=0.0d0
      cell%vfwg_y(:)=0.0d0
      cell%vfwg_z(:)=0.0d0
      cell%vfwl_x(:)=0.0d0
      cell%vfwl_y(:)=0.0d0
      cell%vfwl_z(:)=0.0d0
!
      p1=pitch/(pitch-do_tube)
      p2=pitch/(pi*do_tube)
!
!.....Input parameters for the tube support plate
!
      cs=0.55d0
      afad2=5.0d0
!
      itube(:)=0
!
      DO i=1,ncell_fluid
         g=igr(i)
         IF(g.eq.0) CYCLE
         m=j1d(i)
!
         vfwg_x_i=0.0d0
         vfwg_y_i=0.0d0
         vfwg_z_i=0.0d0
         vfwl_x_i=0.0d0
         vfwl_y_i=0.0d0
         vfwl_z_i=0.0d0
         vfwg_i=0.0d0
         vfwl_i=0.0d0
!
         IF(mult_cell(g,m).eq.0)THEN
            qu=cell%quals(i)
            rl=cell%rhol(i)
            rg=cell%rhog(i)
            rm=cell%rhom(i)
            ag=cell%alphag(i)
            al=cell%alphal(i)
            pr=cell%p(i)
            vxl=vl_o(i,1)
            vxg=vg_o(i,1)
            vyl=vl_o(i,2)
            vyg=vg_o(i,2)
            vzl=vl_o(i,3)
            vzg=vg_o(i,3)
            visl=cell%lviscosl(i)
            visg=cell%lviscosg(i)
            vp=volp(i)
            poro=porosity(i)
            ul=ul_o(i)
            ug=ug_o(i)
         ELSE
            k=mult_cell(g,m)
            m1=mult_3d_cell1(mult_cell(g,m))
            m2=mult_3d_cell2(mult_cell(g,m))
            qu=f1_mult(k)*cell%quals(m1)+f2_mult(k)*cell%quals(m2)
            rl=f1_mult(k)*cell%rhol(m1)+f2_mult(k)*cell%rhol(m2)
            rg=f1_mult(k)*cell%rhog(m1)+f2_mult(k)*cell%rhog(m2)
            rm=f1_mult(k)*cell%rhom(m1)+f2_mult(k)*cell%rhom(m2)
            ag=f1_mult(k)*cell%alphag(m1)+f2_mult(k)*cell%alphag(m2)
            al=f1_mult(k)*cell%alphal(m1)+f2_mult(k)*cell%alphal(m2)
            vxl=f1_mult(k)*vl_o(m1,1)+f2_mult(k)*vl_o(m2,1)
            vxg=f1_mult(k)*vg_o(m1,1)+f2_mult(k)*vg_o(m2,1)
            vyl=f1_mult(k)*vl_o(m1,2)+f2_mult(k)*vl_o(m2,2)
            vyg=f1_mult(k)*vg_o(m1,2)+f2_mult(k)*vg_o(m2,2)
            vzl=f1_mult(k)*vl_o(m1,3)+f2_mult(k)*vl_o(m2,3)
            vzg=f1_mult(k)*vg_o(m1,3)+f2_mult(k)*vg_o(m2,3)
            visl=f1_mult(k)*cell%lviscosl(m1)+f2_mult(k)*cell%lviscosl(m2)
            visg=f1_mult(k)*cell%lviscosg(m1)+f2_mult(k)*cell%lviscosg(m2)
            vp=volp(m1)+volp(m2)
            poro=f1_mult(k)*porosity(m1)+f2_mult(k)*porosity(m2)
            ul=f1_mult(k)*ul_o(m1)+f2_mult(k)*ul_o(m2)
            ug=f1_mult(k)*ug_o(m1)+f2_mult(k)*ug_o(m2)
         ENDIF
!   
!........Two-phase multiplier(TPM) (Thom Correlation)
!
         rlg=rl/rg
         a=1.0d0
         b=1.0531*rlg+1.05455d0
         IF(qu.le.0.01d0)THEN
            c=-0.0048d0*rlg+1.10755d0
         ELSE
            c=-0.0016d0*rlg+1.0492d0
         ENDIF
         TPM=a+b*qu**c
         TPM=TPM*rm/rl
         IF(pr.le.5.0d6.or.pr.ge.7.6d6) TPM=1.0d0
!
         TPM=TPM/vp
         denom=ag*visg+al*visl
         arl=al*rl
         arg=ag*rg
!
!........Vertical U-tube region
!
         IF(ih(g,m).eq.1)THEN
!
!...........Axial flow pressure drop
!
            velz=DSQRT(ag*vzg*vzg+al*vzl*vzl)
            rez=DMAX1(1.d0,rm*velz*hyd_d(1)/denom)
            fz=0.024d0*rez**(-0.2d0)*ht_area(g,m)
            vfwl_z_i=fz*arl*DABS(vzl)*vzl*TPM
            vfwg_z_i=fz*arg*DABS(vzg)*vzg*TPM
!
!...........Cross flow pressure drop
!
            velx=DSQRT(ag*vxg*vxg+al*vxl*vxl)
            vely=DSQRT(ag*vyg*vyg+al*vyl*vyl)
            rex=DMAX1(1.d0,rm*velx*hyd_d(1)/denom)
            rey=DMAX1(1.d0,rm*vely*hyd_d(1)/denom)
            fx=0.432d0*rex**(-0.205d0)*ht_area(g,m)
            fy=0.432d0*rey**(-0.205d0)*ht_area(g,m)
            CV=poro*p1
            CV=CV*CV
            CA=4.d0*poro*p2
            fx=fx*CA*CV
            fy=fy*CA*CV
            vfwl_x_i=fx*arl*DABS(vxl)*vxl*TPM
            vfwg_x_i=fx*arg*DABS(vxg)*vxg*TPM
            vfwl_y_i=fy*arl*DABS(vyl)*vyl*TPM
            vfwg_y_i=fy*arg*DABS(vyg)*vyg*TPM
!
!........Horizontal U-tube region
!
         ELSEIF(ih(g,m).eq.2)THEN
!
!...........Axial/Cross flow pressure drop
!
            velx=DSQRT(ag*vxg*vxg+al*vxl*vxl)
            vely=DSQRT(ag*vyg*vyg+al*vyl*vyl)
            velz=DSQRT(ag*vzg*vzg+al*vzl*vzl)
            rex=DMAX1(1.d0,rm*velx*hyd_d(1)/denom)
            rey=DMAX1(1.d0,rm*vely*hyd_d(1)/denom)
            rez=DMAX1(1.d0,rm*velz*hyd_d(1)/denom)
            fx=0.432d0*rex**(-0.205d0)*ht_area(g,m)
            fy=0.432d0*rey**(-0.205d0)*ht_area(g,m)
            fz=0.432d0*rez**(-0.205d0)*ht_area(g,m)
            CV=poro*p1
            CV=CV*CV
            CA=4.d0*poro*p2
            fx=fx*CA*CV
            fy=fy*CA*CV
            fz=fz*CA*CV
            vfwg_x_i=fx*arg*DABS(vxg)*vxg*TPM
            vfwl_x_i=fx*arl*DABS(vxl)*vxl*TPM
            vfwg_y_i=fy*arg*DABS(vyg)*vyg*TPM
            vfwl_y_i=fy*arl*DABS(vyl)*vyl*TPM
            vfwg_z_i=fz*arg*DABS(vzg)*vzg*TPM
            vfwl_z_i=fz*arl*DABS(vzl)*vzl*TPM
!
         ENDIF
!
!
!........Anti-vibration bars
!
         IF(iavb(g,m).eq.1)THEN
         ENDIF
!
!
!........U-tube support plate (horizontal)
!
         IF(ihp(g,m).eq.1)THEN
!
            vfwg_i=cs*afad2*h_tube(g,m)*arg*ug*TPM
            vfwl_i=cs*afad2*h_tube(g,m)*arl*ul*TPM
!
         ENDIF
!
!
!........U-tube support plate (vertical)
!
         IF(izp(g,m).eq.1)THEN
!
            vfwg_i=cs*afad2*h_tube(g,m)*arg*ug*TPM
            vfwl_i=cs*afad2*h_tube(g,m)*arl*ul*TPM
!
         ENDIF
!
         IF(mult_cell(g,m).eq.0)THEN
            cell%vfwg_x(i)=vfwg_x_i
            cell%vfwg_y(i)=vfwg_y_i
            cell%vfwg_z(i)=vfwg_z_i
            cell%vfwl_x(i)=vfwl_x_i
            cell%vfwl_y(i)=vfwl_y_i
            cell%vfwl_z(i)=vfwl_z_i
            IF(ihp(g,m).eq.1.or.izp(g,m).eq.1)THEN
               cell%vfwg(i)=vfwg_i
               cell%vfwl(i)=vfwl_i
            ENDIF
         ELSE
            cell%vfwg_x(m1)=vfwg_x_i
            cell%vfwg_y(m1)=vfwg_y_i
            cell%vfwg_z(m1)=vfwg_z_i
            cell%vfwl_x(m1)=vfwl_x_i
            cell%vfwl_y(m1)=vfwl_y_i
            cell%vfwl_z(m1)=vfwl_z_i
            cell%vfwg_x(m2)=vfwg_x_i
            cell%vfwg_y(m2)=vfwg_y_i
            cell%vfwg_z(m2)=vfwg_z_i
            cell%vfwl_x(m2)=vfwl_x_i
            cell%vfwl_y(m2)=vfwl_y_i
            cell%vfwl_z(m2)=vfwl_z_i
            IF(ihp(g,m).eq.1.or.izp(g,m).eq.1)THEN
               cell%vfwg(m1)=vfwg_i
               cell%vfwl(m1)=vfwl_i
               cell%vfwg(m2)=vfwg_i
               cell%vfwl(m2)=vfwl_i 
            ENDIF
         ENDIF
!
         itube(i)=1
!
      ENDDO
!
!.....Downcomer region (artificial values are assigned)
!
      DO i=1,ncell_fluid
         IF(itube(i).eq.0.and.xloc(i,3).lt.12.33)THEN
            IF(idc(i).ge.2)THEN
               cell%vfwl(i)=3.0d0*cell%rhol(i)*ug_o(i)
               cell%vfwg(i)=3.0d0*cell%rhog(i)*ul_o(i)
            ELSE
               cell%vfwl(i)=3.0d0*cell%alphal(i)*cell%rhol(i)*ug_o(i)
               cell%vfwg(i)=3.0d0*cell%alphag(i)*cell%rhog(i)*ul_o(i)
            ENDIF
         ENDIF
      ENDDO
!
      END SUBROUTINE udfn_sg_mom_wall
!
