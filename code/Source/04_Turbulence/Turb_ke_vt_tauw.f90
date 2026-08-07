      SUBROUTINE turb_ke_vt_tauw_liq(tw,vt,v,ke,yp)
!
!     This routine calculates tangential velocity of the cell next to wall.
!
      USE VOL_DATA              
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zparam    , ONLY: ndim,cmu,cappa,clog
      USE Zbc_index , ONLY: icell_type,iface_wall0
      USE Zconst1   , ONLY: iturb
      USE Zface     , ONLY: Kepsilon_real
      USE Zturb     , ONLY: cmu_real
      USE Zvec_geo  , ONLY: xn_nf
!
      IMPLICIT NONE
!.....Input
      REAL(8) v(ncell_fp,ndim),ke(ncell_fp),yp(ncell_fp)
!.....Output
      REAL(8) tw(ncell_fp),vt(ncell_fp)
!.....Local variable 
      INTEGER i,i0
!.....Local arrays 
      REAL(8) xn(ndim),vn(ndim)
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1) THEN
            i0=iface_wall0(i)
            IF(ndim.eq.2) THEN
               xn(1)=xn_nf(i0,1)
               xn(2)=xn_nf(i0,2)
               vn(1)=v(i,1)
               vn(2)=v(i,2)
            ELSE
               xn(1)=xn_nf(i0,1)
               xn(2)=xn_nf(i0,2)
               xn(3)=xn_nf(i0,3)
               vn(1)=v(i,1)
               vn(2)=v(i,2)
               vn(3)=v(i,3)
            ENDIF
            CALL v_tangent(vt(i),xn,vn)
!
            IF(iturb.ne.Kepsilon_real)THEN      !std k-e & RNG k-e model
               tw(i)=cell%rhol(i)*cmu**0.25d0*DSQRT(ke(i))*cappa*vt(i)/DLOG(clog*yp(i))
            ELSE
               tw(i)=cell%rhol(i)*cmu_real(i)**0.25d0*DSQRT(ke(i))*cappa*vt(i)/DLOG(clog*yp(i))
            ENDIF                
!
!...........Tangential velocity of neighbor cells to be used for gradient calculation
!
!           DO j=i_neigh(i),i_neigh(i+1)-1
!              IF(j_nbcon(j).eq.0) THEN
!                   
!                 k=neigh(j)
!                  
!                 IF(icell_type(k).ne.1) THEN
!                    CALL v_tangent(vt(k),xn,v(k,:))
!                 ENDIF
!                  
!              ENDIF
!           ENDDO
!
        ENDIF
      ENDDO
!
      END SUBROUTINE turb_ke_vt_tauw_liq
!
      SUBROUTINE turb_ke_vt_tauw_gas(tw,vt,v,ke,yp)
!
!     This routine calculates tangential velocity of the cell next to wall.
!
      USE VOL_DATA              
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zparam    , ONLY: ndim,cmu,cappa,clog
      USE Zbc_index , ONLY: nbcon,icell_type,iface_wall0
      USE Zconst1   , ONLY: iturb
      USE Zface     , ONLY: Kepsilon_real
      USE Znum_cell , ONLY: i_neigh,neigh
      USE Zturb     , ONLY: cmug_real
      USE Zvec_geo  , ONLY: xn_nf
!
      IMPLICIT NONE
!.....Input
      REAL(8) :: v(ncell_fp,ndim),ke(ncell_fp),yp(ncell_fp)
!.....Output
      REAL(8) :: tw(ncell_fp),vt(ncell_fp)
!.....Local variable 
      INTEGER :: i,j,k,i0,j0
!.....Local arrays 
      REAL(8) xn(ndim),vn(ndim)
!
      DO i=1,ncell_fluid
         IF(icell_type(i).eq.1) THEN
!             
            i0=iface_wall0(i)
            IF(ndim.eq.2) THEN
               xn(1)=xn_nf(i0,1)
               xn(2)=xn_nf(i0,2)
               vn(1)=v(i,1)
               vn(2)=v(i,2)
            ELSE
               xn(1)=xn_nf(i0,1)
               xn(2)=xn_nf(i0,2)
               xn(3)=xn_nf(i0,3)
               vn(1)=v(i,1)
               vn(2)=v(i,2)
               vn(3)=v(i,3)
            ENDIF
            CALL v_tangent(vt(i),xn,vn)
!
            IF(iturb.ne.Kepsilon_real)THEN      !std k-e & RNG k-e model
               tw(i)=cell%rhog(i)*cmu**0.25d0*SQRT(ke(i))*cappa*vt(i)/LOG(clog*yp(i))
            ELSE
               tw(i)=cell%rhog(i)*cmug_real(i)**0.25d0*SQRT(ke(i))*cappa*vt(i)/LOG(clog*yp(i))
            ENDIF                
!
!...........Tangential velocity of neighbor cells to be used for gradient calculation
!
            j0=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
               IF(nbcon(j).eq.0) THEN
                  k=neigh(j)
                  IF(icell_type(k).ne.1) THEN
                     IF(ndim.eq.2) THEN
                        vn(1)=v(k,1)
                        vn(2)=v(k,2)
                     ELSE
                        vn(1)=v(k,1)
                        vn(2)=v(k,2)
                        vn(3)=v(k,3)
                     ENDIF
                     CALL v_tangent(vt(k),xn,vn)
                  ENDIF
!                  
               ENDIF
            ENDDO
!
        ENDIF
      ENDDO
!
      END SUBROUTINE turb_ke_vt_tauw_gas
