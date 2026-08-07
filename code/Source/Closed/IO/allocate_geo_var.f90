!
      SUBROUTINE allocate_geo_var(ncell_fps,ncell_fluids)
!
!     This routine allocates variables containing geomeric information.
!
      USE Zmpi         , ONLY: maxmt_fp,maxmt_fluid
      USE Zparam       , ONLY: ndim,nb_max
      USE Zbc_index    , ONLY: ngrad,icell_type,iface_wall,iface_wall1,cellvin,cellpin,vin_norm,vin_mfr, &
                               iface_wall0
      USE Zcoord1      , ONLY: xloc_m
      USE Zcoord2      , ONLY: fac,fac1,cell_leng
      USE Zcoord3      , ONLY: svp,volp,aporous,floss,volr,volpr
      USE Zcoord4      , ONLY: sa,saa,sap,dji,dji_a,dji_x,sad,dnj      
      USE Zgrad_ls_c2d , ONLY: a11_2,a12_2,a22_2,a_2,det_2
      USE Zgrad_ls_c3d , ONLY: a11_3,a12_3,a22_3,a13_3,a23_3,a33_3,a_3,det_3,lsindex
      USE Znormal      , ONLY: xn
!
      IMPLICIT NONE
!
      INTEGER ncell_fps,ncell_fluids
      INTEGER n,n1
!
      n=ncell_fps
      n1=ncell_fluids
!
!.....Zcoord1
!
      ALLOCATE(xloc_m(maxmt_fp,ndim))
      xloc_m(:,:)=0.0d0
!
!.....Zcoord2
!
      ALLOCATE(fac(maxmt_fp),fac1(maxmt_fp),cell_leng(n,ndim))
      fac(:)=0.0d0
      fac1(:)=0.0d0
      cell_leng(:,:)=0.0d0
!
!.....Zcoord3
!
      ALLOCATE(svp(maxmt_fluid,ndim),volp(n),aporous(n),floss(n,ndim))                
      ALLOCATE(volr(n),volpr(n))
      svp(:,:)=0.0d0
      volp(:)=0.0d0
      aporous(:)=0.0d0
      floss(:,:)=0.0d0
      volr(:)=0.0d0
      volpr(:)=0.0d0
!
!.....Zcoord4
!
      ALLOCATE(sa(maxmt_fluid),saa(maxmt_fluid),sap(maxmt_fluid),dji(maxmt_fluid),dji_a(maxmt_fluid),dji_x(maxmt_fluid,ndim))
      ALLOCATE(sad(maxmt_fluid),dnj(maxmt_fluid,ndim))
      sa(:)=0.0d0
      saa(:)=0.0d0
      sap(:)=0.0d0
      dji(:)=0.0d0
      dji_a(:)=0.0d0
      dji_x(:,:)=0.0d0
      sad(:)=0.0d0
      dnj(:,:)=0.0d0
!
!.....Znormal
!
      ALLOCATE(xn(maxmt_fp,ndim))
      xn(:,:)=0.0d0
!
!.....Zbc_index
!
      ALLOCATE(ngrad(n),icell_type(n),iface_wall(n),iface_wall0(n),iface_wall1(n1))
      ALLOCATE(cellvin(nb_max),cellpin(nb_max),vin_norm(nb_max),vin_mfr(nb_max))
      ngrad(:)=0
      icell_type(:)=0
      iface_wall(:)=0
      iface_wall0(:)=0
      iface_wall1(:)=0
      cellvin(:)=0
      cellpin(:)=0
      vin_norm(:)=0
      vin_mfr(:)=0
!
!.....Zgrad_ls_c2d
!
      ALLOCATE(a11_2(n),a12_2(n),a22_2(n),a_2(n),det_2(n))
      a11_2(:)=0.0d0
      a12_2(:)=0.0d0
      a22_2(:)=0.0d0
      a_2(:)=0.0d0
      det_2(:)=0.0d0
!
!.....Zgrad_ls_c3d
!
      ALLOCATE(a11_3(n1),a12_3(n1),a13_3(n1),a22_3(n1),a23_3(n1),a33_3(n1),a_3(n),det_3(n),lsindex(n))
      a11_3(:)=0.0d0
      a12_3(:)=0.0d0
      a13_3(:)=0.0d0
      a22_3(:)=0.0d0
      a23_3(:)=0.0d0
      a33_3(:)=0.0d0
      a_3(:)=0.0d0
      det_3(:)=0.0d0
      lsindex(:)=0
!
      RETURN
      END SUBROUTINE allocate_geo_var
