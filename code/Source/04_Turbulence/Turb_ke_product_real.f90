!
      SUBROUTINE turb_ke_product_real(ik,pro,tw,strn,ustar,w_real)
!
!     This routine calculates mean shear generation rate of k
!     at the cv center.
!
      USE VOL_DATA                 
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim,ns,cappa
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas
      USE Zndforce     , ONLY: d_bfc
      USE Zturb        , ONLY: utau,utaug
      USE Zvector      , ONLY: vl_o,vg_o
      USE Zbc_index    , ONLY: icell_type
!
      IMPLICIT NONE
!     input
      REAL(8) :: tw(ncell_fp)
!     output
      REAL(8) :: pro(ncell_fp),strn(ncell_fp),ustar(ncell_fp),w_real(ncell_fp)
!
      INTEGER ik,i,ix
!      
      REAL(8) vii,vij,vxy,vxz,vyz
      REAL(8) vxy2,vxz2,vyz2,vij2 ! for Realizable K-e
      REAL(8) S_bar,SSS,s11,s12,s13,s21,s22,s23,s31,s32,s33  !skchoi        ! for Realizable K-e
!     local arrays
      REAL(8) :: fie(ns),tviscos(ncell_fp)
      REAL(8) :: dvdx(ncell_fp,ndim,ndim)
!
!
!.....Initialize kinetic enerfy production
!
      pro(:)=0.0d0
      strn(:)=0.0d0
      ustar(:)=0.0d0
      w_real(:)=0.0d0
      dvdx(:,:,:)=0.0d0
      fie(:)=0.0d0
      tviscos(:)=0.0d0
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
!.....Production of turbulence kinetic energy
!
      DO i=1,ncell_fluid
         vii=0.0d0
         vij=0.0d0
!         
         DO ix=1,ndim
            vii=vii+2.d0*dvdx(i,ix,ix)*dvdx(i,ix,ix)
         ENDDO
!         
         vxy=dvdx(i,2,1)+dvdx(i,1,2)
         vij=vxy*vxy    
         IF(ndim.eq.3) THEN
            vxz=dvdx(i,ndim,1)+dvdx(i,1,ndim)
            vyz=dvdx(i,ndim,2)+dvdx(i,2,ndim)
            vij=vij+vxz*vxz+vyz*vyz            
         ENDIF
         strn(i)=DSQRT(vii+vij)     
!
!........Realizable k-e
!         
         vxy2=dvdx(i,2,1)-dvdx(i,1,2)
         vij2=vxy2*vxy2   
         IF(ndim.eq.3) THEN
            vxz2=dvdx(i,ndim,1)-dvdx(i,1,ndim)
            vyz2=dvdx(i,ndim,2)-dvdx(i,2,ndim)
            vij2=vij2+vxz2*vxz2+vyz2*vyz2            
         ENDIF         
         ustar(i)=DSQRT(0.5d0*(vii+vij)+0.5d0*vij2)    !skchoi ! U* for Realizable k-e model
!
         s_bar=DSQRT(0.5d0*(vii+vij))
         s11=dvdx(i,1,1)
         s12=0.5d0*(dvdx(i,2,1)+dvdx(i,1,2))
         s22=dvdx(i,2,2)
         s21=0.5d0*(dvdx(i,1,2)+dvdx(i,2,1))
         SSS=s11*(s11**2+3.d0*s12**2)+s22*(s22**2+3.d0*s21**2)               
         IF(ndim.eq.3)THEN
            s13=0.5d0*(dvdx(i,3,1)+dvdx(i,1,3))
            s23=0.5d0*(dvdx(i,3,2)+dvdx(i,2,3))
            s33=dvdx(i,3,3)
            s31=0.5d0*(dvdx(i,1,3)+dvdx(i,3,1))
            s32=0.5d0*(dvdx(i,2,3)+dvdx(i,3,2))            
            SSS=SSS+s11*3.d0*s13**2+s22*3.d0*s23**2 &   
               +s33*(s33**2+3.d0*s31**2+3.d0*s32**2) &
               +6.d0*s12*s23*s31
         ENDIF
         IF(S_bar.eq.0)THEN
            w_real(i)=0.0d0
         ELSE
            w_real(i)=SSS/S_bar**3 !skchoi
         ENDIF
!         
         pro(i)=vii+vij
         pro(i)=tviscos(i)*pro(i)
!
!........Wall boundary condition
!
         IF(icell_type(i).eq.1) THEN
            IF(ik.eq.1)THEN
               pro(i)=tw(i)*utau(i)/cappa/d_bfc(i)
            ELSEIF(ik.eq.2)THEN
               pro(i)=tw(i)*utaug(i)/cappa/d_bfc(i)           
            ENDIF            
         ENDIF
!
      ENDDO
!
      RETURN
      END SUBROUTINE turb_ke_product_real
