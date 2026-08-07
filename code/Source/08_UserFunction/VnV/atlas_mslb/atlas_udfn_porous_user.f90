!
      SUBROUTINE atlas_udfn_porous_user(vol_tmp,poro_tmp,nmaterial_tmp,nzone_tmp) 
!.....This routine change the cell value of somaGrid.
!
      USE Zzone       , ONLY: ncell_fluid_all,ncell_cond_all      
      USE Zcore       , ONLY: myrank
      USE Zparam      , ONLY: nn,ndim
      USE Znum_cell   , ONLY: i_neigh_tmp, &
                              perm_tmp1
      USE Zconst1     , ONLY: vv_prob 
      USE Zcoord1     , ONLY: xloc_tmp
!
      IMPLICIT NONE
!           
!     output
      INTEGER nmaterial_tmp(nn),nzone_tmp(nn)
      REAL(8) vol_tmp(nn),poro_tmp(nn)
!     local variables
      INTEGER i,j
      REAL(8) corevol
      REAL(8) rdc,zh,zl
!     local arrays
      REAL(8) x(ndim),porosity,r,porosity_core
!
      IF(vv_prob.eq.'atlas_mc_porous')THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'##atlas_udfn_porous_user is called.' 
          zh=-0.249d0
          zl=-3.046d0
          rdc=0.168d0
          porosity_core=0.58d0      
      ELSEIF(vv_prob.eq.'pwr_mc_poro')THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'##pwr_mc_poro: atlas_udfn_porous_user is called.'   
          zh=-1.1d0    !OLD
          zl=-4.9d0
          rdc=1.994d0         
          porosity_core=0.58d0      
      ELSEIF(vv_prob.eq.'apr1400_mc_poro')THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'##apr1400_mc_poro: atlas_udfn_porous_user is called.'   
          zh=-5.4d0     !NEW 
          zl=-9.22d0
          rdc=1.77d0        
          porosity_core=0.58d0
      ELSEIF(vv_prob.eq.'opr1000_mc_poro')THEN
         IF(myrank.eq.0)WRITE(*,"(11x,a)")'##opr1000_mc_poro: atlas_udfn_porous_user is called.'      
          zh=4.0d0
          zl=0.2d0
          rdc=1.7526d0 
          porosity_core=0.58d0      
      ENDIF    
!     
      IF(myrank.eq.0)THEN
          corevol=0.0d0
          DO i=1,nn
             x(:)=xloc_tmp(i,:)
             r=DSQRT(x(1)*x(1)+x(2)*x(2))
             IF(r.lt.rdc.and.x(ndim).lt.zh.and.x(ndim).gt.zl)THEN
                corevol=corevol+vol_tmp(i)*poro_tmp(i)
             ENDIF   
          ENDDO 
          IF(myrank.eq.0)WRITE(*,"(11x,a,1pe17.5,a)")'total core volume=',corevol,'m3'       
      ENDIF
!
      ncell_fluid_all=0
      ncell_cond_all=0
      DO i=1,nn  
          DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
             perm_tmp1(j)=1.0d0
          ENDDO
          x(:)=xloc_tmp(i,:)
          r=DSQRT(x(1)*x(1)+x(2)*x(2))
          IF(r.lt.rdc.and.x(ndim).lt.zh.and.x(ndim).gt.zl)THEN
             porosity=porosity_core
          ELSE
             porosity=1.0d0   
          ENDIF
          poro_tmp(i)=porosity
          IF(porosity.lt.0.999d0)then
              nmaterial_tmp(i)=-4 !
              nzone_tmp(i)=2      !
              DO j=i_neigh_tmp(i),i_neigh_tmp(i+1)-1
                 perm_tmp1(j)=1.0d0
               ENDDO
              !nmaterial: 0=fluid, -1~-10=conductor, -11~=rod ==> related to writing somagrid.in
              !nmaterial: 1. UO2, 2. Zircaloy, 3. Inconel, 4. Stainless steel, 5. Carbon steel
              !nzone: related to subdomain
          ENDIF
         IF(nmaterial_tmp(i).le.0) ncell_fluid_all=ncell_fluid_all+1
         IF(nmaterial_tmp(i).ne.0) ncell_cond_all=ncell_cond_all+1             
     ENDDO     
!       
      RETURN
      ENDSUBROUTINE atlas_udfn_porous_user
