      SUBROUTINE turb_SST_vt_tauw_liq(vt,v)
!
!     This routine calculates tangential velocity of the cell next to wall.
!
      USE VOL_DATA              
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zparam    , ONLY: ndim,cmu,cappa,clog
      USE Zbc_index , ONLY: icell_type,iface_wall0
      USE Zndforce  , ONLY: d_bfc
      USE Zturb     , ONLY: utau,tauw
      USE Zvec_geo  , ONLY: xn_nf
!
      IMPLICIT NONE
!      
!     input
      REAL(8) v(ncell_fp,ndim)
!     output
      REAL(8) vt(ncell_fp)
!     local variable
      INTEGER i,i0
!     local arrays
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
            tauw(i)=cell%lviscosl(i)*DABS(vt(i))/d_bfc(i)
            utau(i)=DSQRT(DABS(tauw(i))/cell%rhol(i))             
!!
!!...........Tangential velocity of neighbor cells to be used for gradient calculation
!!
!            DO m=1,num_neigh(i)
!               IF(nbcon(m,i).eq.0) THEN
!!                   
!                  k=neigh(m,i)
!!                  
!                  IF(icell_type(k).ne.1) THEN
!                     CALL v_tangent(vt(k),xn(:,j,i),v(:,k))
!                  ENDIF
!!                  
!               ENDIF
!            ENDDO
!
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE turb_SST_vt_tauw_liq
      SUBROUTINE turb_SST_vt_tauw_gas(vt,v)
!
!     This routine calculates tangential velocity of the cell next to wall.
!
      USE VOL_DATA              
      USE Zmpi      , ONLY: ncell_fp
      USE Zzone     , ONLY: ncell_fluid
      USE Zparam    , ONLY: ndim,cmu,cappa,clog
      USE Zbc_index , ONLY: icell_type,iface_wall0
      USE Zndforce  , ONLY: d_bfc
      USE Zturb     , ONLY: utaug,tauwg
      USE Zvec_geo  , ONLY: xn_nf
!
      IMPLICIT NONE
!      
!     input
      REAL(8) v(ncell_fp,ndim)
!     output
      REAL(8) vt(ncell_fp)
!     local variable
      INTEGER i,i0
!     local arrays
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
            tauwg(i)=cell%lviscosg(i)*DABS(vt(i))/d_bfc(i)
            utaug(i)=DSQRT(DABS(tauwg(i))/cell%rhog(i))             
!!
!!...........Tangential velocity of neighbor cells to be used for gradient calculation
!!
!            DO m=1,num_neigh(i)
!               IF(nbcon(m,i).eq.0) THEN
!!                   
!                  k=neigh(m,i)
!!                  
!                  IF(icell_type(k).ne.1) THEN
!                     CALL v_tangent(vt(k),xn(:,j,i),v(:,k))
!                  ENDIF
!!                  
!               ENDIF
!            ENDDO
!
         ENDIF
      ENDDO
!
      RETURN
      END SUBROUTINE turb_SST_vt_tauw_gas
