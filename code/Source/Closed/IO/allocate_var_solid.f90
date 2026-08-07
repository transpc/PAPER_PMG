      SUBROUTINE allocate_var_solid(ncell_cond,ncell_ps)
!
!     This routine allocates variables containing physical information.
!
      USE Zmpi          , ONLY: maxmt_ncond,maxmt_ps
      USE Znum_cell     , ONLY: n_fluid
      USE Zparam        , ONLY: ndim
!     USE Zcoord1       , ONLY: xloc_c,xloc_m_c
      USE Zcoord1       , ONLY: xloc_c
      USE Zcoord2       , ONLY: fac_c,fac1_c
      USE Zcoord3       , ONLY: volp_c
      USE Zcoord4       , ONLY: sap_c,dji_a_c
      USE Zqvol         , ONLY: qvol_ice_solid,qconv_sol
      USE Zzone         , ONLY: nmaterial_c,nzone_c
!
      IMPLICIT NONE
!
!     input
      INTEGER ncell_cond,ncell_ps
!
!.....Zqvol
!
      ALLOCATE(qvol_ice_solid(ncell_ps))
      ALLOCATE(qconv_sol(ncell_cond))
!
      qvol_ice_solid(:)=0.0d0
      qconv_sol(:)=0.0d0
!
      ALLOCATE(fac_c(maxmt_ps),fac1_c(maxmt_ps))
      fac_c(:)=0.0d0
      fac1_c(:)=0.0d0
      ALLOCATE(volp_c(ncell_cond),sap_c(maxmt_ncond),nmaterial_c(ncell_ps))
      ALLOCATE(xloc_c(ncell_ps,ndim),n_fluid(ncell_cond),dji_a_c(maxmt_ncond),nzone_c(ncell_ps))
      sap_c(:)=0.0d0
      nmaterial_c(:)=0
      xloc_c(:,:)=0.0d0
      n_fluid(:)=0
      dji_a_c(:)=0.0d0
      nzone_c(:)=0
!
      END SUBROUTINE allocate_var_solid
