!
      SUBROUTINE turb_ke_product(ik,pro,tw,strn)
!
!     This routine calculates mean shear generation rate of k
!     at the cv center.
!
      USE VOL_DATA                 
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim,cappa
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas
      USE Zbc_index    , ONLY: icell_type
      USE Zconst1      , ONLY: lowreynolds
      USE Zndforce     , ONLY: d_bfc
      USE Zturb        , ONLY: utau,utaug
      USE Zvector      , ONLY: vl_o,vg_o
!
      IMPLICIT NONE
!
!     input
      INTEGER ik
      REAL(8) tw(ncell_fp),strn(ncell_fp)
!     output
      REAL(8) pro(ncell_fp)
!     local variables
      INTEGER i
!      
!     REAL(8) vii,vij,vxy,vxz,vyz,vti,dvtdn
      REAL(8) vii,vij,vxy,vxz,vyz
!     local arrays
      REAL(8) :: tviscos(ncell_fp)
      REAL(8) :: dvdx(ncell_fp,ndim,ndim)
!
!
!.....Initialize kinetic enerfy production
!
      pro(:)=0.0d0
!
!.....Calculate velocity gradients
!
      IF(ik.eq.1) THEN
         CALL grad_vel(2,vl_o,dvdx,vb_liq,vin_liq)
         tviscos(:)=cell%tviscosl(:)
      ELSEIF(ik.eq.2) THEN
         CALL grad_vel(1,vg_o,dvdx,vb_gas,vin_gas)
         tviscos(:)=cell%tviscosg(:)
      ENDIF
!
!.....Gradient of tangential velocity at wall
!
!      dvtdx(:,:)=0.0d0
!!      
!      DO i=1,ncell_fluid
!!          
!        IF(icell_type(i).eq.1) THEN
!!            
!            DO j=1,num_neigh(i)
!!                
!               IF(nbcon(j,i).eq.0) THEN
!                  k=neigh(j,i)
!                  vti=fac1(j,i)*vt(i)+fac(j,i)*vt(k)
!               ELSEIF(nbcon(j,i).gt.nin_max) THEN
!                  vti=vt(i)
!               ELSEIF(nbcon(j,i).gt.0) THEN
!                  vti=vt(i)
!               ELSEIF(nbcon(j,i).ge.201) THEN !mcc-pik
!                  vti=vt(i)
!               ELSE
!                  vti=0.0d0
!               ENDIF
!!               
!               DO ix=1,ndim
!                  dvtdx(i,ix)=dvtdx(i,ix)+vti*sv(j,i,ix)
!               ENDDO
!!               
!            ENDDO
!!            
!            DO ix=1,ndim
!               dvtdx(i,ix)=dvtdx(i,ix)/vol(i)
!            ENDDO
!!            
!        ENDIF
!!        
!      ENDDO
!
!.....Production of turbulence kinetic energy
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
!         
            vii=2.d0*(dvdx(i,1,1)**2+dvdx(i,2,2)**2)
!         
            vxy=dvdx(i,2,1)+dvdx(i,1,2)
            vij=vxy**2
            strn(i)=SQRT(vii+vij)
!
            pro(i)=vii+vij
            pro(i)=tviscos(i)*pro(i)
!
!........Wall boundary condition
!
            IF(icell_type(i).eq.1.and.lowreynolds.eq.0) THEN
!            j=iface_wall(i)
!            dvtdn=0.0d0
!            DO ix=1,ndim
!               dvtdn=dvtdn+dvtdx(i,ix)*xn(ix,j,i)
!            ENDDO
!            pro(i)=tw(i)*DABS(dvtdn)
               IF(ik.eq.1)THEN
                  pro(i)=tw(i)*utau(i)/cappa/d_bfc(i)
               ELSEIF(ik.eq.2)THEN
                  pro(i)=tw(i)*utaug(i)/cappa/d_bfc(i)           
               ENDIF
            ENDIF
!
         ENDDO
      ELSE
         DO i=1,ncell_fluid
!         
            vii=2.d0*(dvdx(i,1,1)**2+dvdx(i,2,2)**2+dvdx(i,3,3)**2)
!         
            vxy=dvdx(i,2,1)+dvdx(i,1,2)
            vxz=dvdx(i,3,1)+dvdx(i,1,3)
            vyz=dvdx(i,3,2)+dvdx(i,2,3)
            vij=vxy**2+vxz**2+vyz**2
            strn(i)=SQRT(0.5d0*vij)
!
            pro(i)=vii+vij
            pro(i)=tviscos(i)*pro(i)
!
!........Wall boundary condition
!
            IF(icell_type(i).eq.1.and.lowreynolds.eq.0) THEN
!            j=iface_wall(i)
!            dvtdn=0.0d0
!            DO ix=1,ndim
!               dvtdn=dvtdn+dvtdx(i,ix)*xn(ix,j,i)
!            ENDDO
!            pro(i)=tw(i)*DABS(dvtdn)
               IF(ik.eq.1)THEN
                  pro(i)=tw(i)*utau(i)/cappa/d_bfc(i)
               ELSEIF(ik.eq.2)THEN
                  pro(i)=tw(i)*utaug(i)/cappa/d_bfc(i)           
               ENDIF
            ENDIF
!
         ENDDO
      ENDIF
!
      END SUBROUTINE turb_ke_product
