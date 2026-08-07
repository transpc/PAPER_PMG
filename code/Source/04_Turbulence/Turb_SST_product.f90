!
      SUBROUTINE turb_SST_product(ik,pro,strn)
!
!     This routine calculates mean shear generation rate of k
!     at the cv center.
!
      USE VOL_DATA                 
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zparam       , ONLY: ndim,cappa
      USE Zb_condition , ONLY: vb_liq,vin_liq,vb_gas,vin_gas
      USE Zvector      , ONLY: vl_o,vg_o
!
      IMPLICIT NONE
!
      INTEGER ik,i
!      
      REAL(8) vii,vij,vxy,vxz,vyz
!     local arrays 
      REAL(8) :: pro(ncell_fp),strn(ncell_fp)
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
!.....Production of turbulence kinetic energy
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
!         
            vii=2.d0*(dvdx(i,1,1)**2+dvdx(i,2,2)**2)
!         
            vxy=dvdx(i,2,1)+dvdx(i,1,2)
            vij=vxy**2
!
            pro(i)=vii+vij
            pro(i)=tviscos(i)*pro(i)
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
            strn(i)=SQRT(vii+vij)
!
            pro(i)=vii+vij
            pro(i)=tviscos(i)*pro(i)
!!
!!........Wall boundary condition - not needed for k-w model
!!
!         IF(icell_type(i).eq.1.and.lowreynolds.ne.1) THEN
!            j=iface_wall(i)
!            dvtdn(i)=0.0d0
!            DO ix=1,ndim
!               dvtdn(i)=dvtdn(i)+dvtdx(ix,i)*xn(ix,j,i)
!            ENDDO
!!            pro(i)=tauw(i)*DABS(dvtdn(i))
!!            pro(i)=tauw(i)*utau(i)/cappa/d_bfc(i)
!
!         ENDIF
!
         ENDDO
      ENDIF
!
      END SUBROUTINE turb_SST_product
