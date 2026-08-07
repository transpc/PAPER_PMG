      SUBROUTINE nbcon_change_start
!     mapping index of _inl,_out,_adw
      USE Zmpi          , ONLY: maxmt_fluid
      USE Zzone         , ONLY: ncell_fluid
      USE Zvec_param    , ONLY: nf_inl,nf_out,nf_adw,nf_flux,                &
                                nf_nbcon_change,nf_nbcon_change_flux,        &
                                nf_inl_old,nf_out_old,nf_adw_old,nf_flux_old
      USE Znum_cell     , ONLY: i_neigh,                                            &
                                istart_nf,istart_nbcon_nf,istart_svp_nf,            &
                                istart_nf_old,istart_nbcon_nf_old,istart_svp_nf_old
      USE Zbc_index     , ONLY: nbcon,nbcon_old
!
      IMPLICIT NONE
!.....Local variables
      INTEGER i,j,j0,nf_number
!
      ALLOCATE(nbcon_old(maxmt_fluid))
!
!.....Copy old nbcon
!
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            nbcon_old(j)=nbcon(j)
         ENDDO
      ENDDO
      nf_nbcon_change=nf_inl+nf_out+nf_adw
      nf_nbcon_change_flux=nf_inl+nf_out
      nf_inl_old=nf_inl
      nf_out_old=nf_out
      nf_adw_old=nf_adw
      nf_flux_old=nf_flux
      DO nf_number=2,4
         istart_nf_old(1,nf_number)=istart_nf(1,nf_number)
         istart_nf_old(2,nf_number)=istart_nf(2,nf_number)
         istart_nbcon_nf_old(nf_number)=istart_nbcon_nf(nf_number)
         istart_svp_nf_old(nf_number)=istart_svp_nf(nf_number)
      ENDDO
!     write(*,*) 'old',nf_inl_old,nf_out_old,nf_adw_old,nf_flux_old
!
      END SUBROUTINE nbcon_change_start
!
      SUBROUTINE nbcon_change_end
!
      USE Zmpi        , ONLY: ncell_fp,maxmt_fluid,maxmt_nfluid,maxmt_cell
      USE Zzone       , ONLY: ncell_fluid,ncell_fluid_all
      USE Zparam      , ONLY: ndim,nin_max,nb_max,nb_mars,nb_sym
      USE Zturb       , ONLY: walln,walln2,wallnr
!     USE Zcoord2     , ONLY: fac,fac1
      USE Znum_cell   , ONLY: i_neigh,i_neigh_tmp,j_nbcon_tmp, &
                              istart_nf,istart_nbcon_nf,       &
                              istart_nf_old
!
      USE Zbc_index   , ONLY: nbcon,nbcon_old,npb,ngrad,icell_type,iface_wall, &
                              iface_wall0
!
      USE Zvec_param  , ONLY: nf_non,nf_mcc,nf_inl,nf_out,nf_flux,nf_flux_old,&
                              nf_nbcon_change,nf_nbcon_change_flux,           &
                              nf_inl_old,nf_flux_old
      USE Zvec_index  , ONLY: left_nf,nbcon_nf, &
                              delta_npb
      USE Zvec_major  ,ONLY: flux_l_nf,flux_g_nf,flux_d_nf,                    &
                             flux_l_nf_o,flux_g_nf_o,flux_d_nf_o,              &
                             liq_conv_nf,vap_conv_nf,drp_conv_nf,              &
                             ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf,                 &
                             al_conv_nf,ad_conv_nf,void_conv_nf,quala_conv_nf, &
                             lbor_conv_nf
      USE Zvec_scalar , ONLY: arli_nf,argi_nf,ardi_nf
      USE Zscalar_coeff , ONLY: sfg_nf,sfl_nf,sfd_nf,   &
                                sfg6_nf,sfl6_nf,sfd6_nf
      USE Zvec_geo    , ONLY: xn_nf,sv_nf,svp_nf,  &
                              xfc_nf,              &
                              dxfc_nf,             &
                              sap_nf,sa_nf,saa_nf, &
                              perm_out,perm_nf,    &
                              dji_nf,djia_nf,      &
                              dji_x_nf,            &        
                              fac_nf,fac1_nf
     USE Zuserdefined , ONLY: vel_bc_profile_inl
     USE Zbc_index    , ONLY: l_horizontal_outlet_init
     USE Zcore        , ONLY: np
!
      IMPLICIT NONE
!.....Local variables
      INTEGER :: i,j,ii,j0,i0,j1
      INTEGER :: itype
      INTEGER :: nf_number,istart,isize,istart2,i1,i2
      REAL*8  a,dr_min,awalln
!.....Local arrays
      REAL*8 x_nf_old(nf_nbcon_change),y_nf_old(nf_nbcon_change),z_nf_old(nf_nbcon_change)
      REAL*8 buff_l(nf_non+nf_mcc),buff_g(nf_non+nf_mcc),buff_d(nf_non+nf_mcc)
!
!.....Signal nbcon_change for i_horizontal_outlet
!
      l_horizontal_outlet_init=.TRUE. 
!     get the index mapping new to old inl,out,adw and new  nbcon_nf,left_nf,jneigh_nf
      CALL vectorize_nbcon_change_index
!===> recompute npb,ngrad
      DO i=1,ncell_fluid
         npb(i)=0
         ngrad(i)=0
      ENDDO
      DO i=1,ncell_fluid
!DIR$ NOVECTOR
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).gt.nin_max.and.nbcon(j).ne.nb_sym.and.nbcon(j).le.nb_max)THEN
!
!..............Pressure boundary
!
               ngrad(i)=2
            ELSEIF(ngrad(i).ne.2.and.nbcon(j).ne.0.and.nbcon(j).ne.nb_sym.and.nbcon(j).lt.nb_mars)THEN
!
!..............Except pressure, cell, symmetric boundary
!  
               ngrad(i)=1
           ENDIF
         ENDDO
      ENDDO
      nf_number=3
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         npb(ii)=nbcon_nf(i2)-nin_max
      ENDDO
!      
      IF(ALLOCATED(delta_npb)) DEALLOCATE(delta_npb)
      ALLOCATE(delta_npb(ncell_fp))
      DO i=1,ncell_fp
         IF(npb(i).ne.0)THEN
            delta_npb(i)=0.0d0
         ELSE
            delta_npb(i)=1.0d0
         ENDIF
      ENDDO
!.....Get global nbcon neeeded in write_fieldview
      IF(np.gt.1) THEN
         CALL gatherv_csr_i(nbcon,maxmt_fluid,j_nbcon_tmp,maxmt_cell,  &
                            ncell_fluid_all,maxmt_nfluid,i_neigh_tmp,0)
      ELSE
         DO i=1,maxmt_fluid
            j_nbcon_tmp(i)=nbcon(i)
         ENDDO
      ENDIF
!.....Reallocate vel_bc_profile_inl, and reinitialize
      IF(nf_inl.ne.nf_inl_old) THEN
         IF(nf_inl_old.gt.0) DEALLOCATE(vel_bc_profile_inl)
         IF(nf_inl.gt.0) ALLOCATE(vel_bc_profile_inl(nf_inl)) 
      ENDIF
      CALL udfn_vel_bc_profile_inl
!.....Reallocate sfl_nf,sfg_nf,sfd_nf
      IF(nf_flux.ne.nf_flux_old) THEN
         DEALLOCATE(sfl_nf,sfg_nf,sfd_nf)
         DEALLOCATE(sfl6_nf,sfg6_nf,sfd6_nf)
         ALLOCATE(sfl_nf(nf_flux,5),sfg_nf(nf_flux,5),sfd_nf(nf_flux,5))
         ALLOCATE(sfl6_nf(nf_flux),sfg6_nf(nf_flux),sfd6_nf(nf_flux))
      ENDIF
      CALL vectorize_nbcon_change_copy2d_nf(dji_nf,djia_nf)
      CALL vectorize_nbcon_change_copy2d_nf(fac_nf,fac1_nf)
      CALL vectorize_nbcon_change_copy_ndim_nf(dji_x_nf)
      CALL vectorize_nbcon_change_copy_ndim_nf(xn_nf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!........Define cell type depending on face condition
!        icell_type(i) = :
!        0 , If the cell has only computing faces
!        1 , If the cell has has wall boundary faces
!        2 , If the cell has has inlet boundary faces
!        3 , If the cell has has outlet boundary faces
!        4 , If the cell has has symmetry boundary faces
!        5 , If the cell has has mcc boundary faces
!
      DO i=1,ncell_fluid
         itype=0
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(nbcon(j).lt.0) THEN
               itype=1
            ELSE
               IF(nbcon(j).gt.0.and.nbcon(j).le.nin_max)THEN
                  IF(itype.eq.0.or.itype.ge.2)itype=2
               ELSEIF(nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)THEN
                  IF(itype.eq.0.or.itype.ge.3)itype=3
               ELSEIF(nbcon(j).eq.101)THEN
                  IF(itype.eq.0.or.itype.ge.4)itype=4
               ELSEIF(nbcon(j).eq.201)THEN
                  IF(itype.eq.0.or.itype.ge.5)itype=5
               ENDIF
            ENDIF
         ENDDO
         icell_type(i)=itype
!            
!...........Save the face index and calculate distance from the wall
!when a cell has wall boundary face
!
         dr_min=huge(0.d0)
         j0=0
         i0=0
         j1=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
             IF((icell_type(i).eq.1 .and.  nbcon(j).lt.0)                               .or. &
                (icell_type(i).eq.2 .and. (nbcon(j).gt.0.and.nbcon(j).le.nin_max))      .or. &
                (icell_type(i).eq.3 .and. (nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)) .or. &
                (icell_type(i).eq.4 .and.  nbcon(j).eq.101))THEN
               CALL get_vector_disp(j-j1,i,i1)
               i1=ABS(i1)
               IF(dji_nf(i1).lt.dr_min) THEN 
                  dr_min=dji_nf(i1)
                  j0=j-j1
                  i0=i1
               ENDIF
             ENDIF
         ENDDO
         iface_wall(i)=j0
         iface_wall0(i)=i0
      ENDDO
      IF(ndim.eq.2) THEN
        DO i=1,ncell_fluid
         j0=iface_wall(i)
         i0=iface_wall0(i)
         IF(j0.ne.0) THEN
            CALL get_vector_disp(j0,i,i1) 
            i1=ABS(i1)
            a=dji_x_nf(i1,1)*xn_nf(i1,1)+dji_x_nf(i1,2)*xn_nf(i1,2)
!           awalln   =fac (j0,i)*a
!           walln2(i)=fac1(j0,i)*a
            awalln   =fac_nf (i0)*a
            walln2(i)=fac1_nf(i0)*a
            j1=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
              IF(j-j1.ne.j0) THEN
               IF(((icell_type(i).eq.1) .and. nbcon(j).lt.0)                               .or. &
                   (icell_type(i).eq.2 .and. (nbcon(j).gt.0.and.nbcon(j).le.nin_max))      .or. &
                   (icell_type(i).eq.3 .and. (nbcon(j).gt.nin_max.and.nbcon(j).le.nb_max)) .or. &
                  ((icell_type(i).eq.4) .and. nbcon(j).eq.101)) THEN
                  awalln=DMIN1(dji_nf(i1),awalln)
               ENDIF
              ENDIF
            ENDDO
            walln(i)=awalln
            wallnr(i)=1.d0/awalln
         ENDIF
        ENDDO
      ELSE
        DO i=1,ncell_fluid
         j0=iface_wall(i)
         i0=iface_wall0(i)
         IF(j0.ne.0) THEN
            CALL get_vector_disp(j0,i,i1) 
            i1=ABS(i1)
            a=dji_x_nf(i1,1)*xn_nf(i1,1)+dji_x_nf(i1,2)*xn_nf(i1,2)+dji_x_nf(i1,3)*xn_nf(i1,3)
!           awalln   =fac (j0,i)*a
!           walln2(i)=fac1(j0,i)*a
            awalln   =fac_nf (i0)*a
            walln2(i)=fac1_nf(i0)*a
            j1=i_neigh(i)-1
            DO j=i_neigh(i),i_neigh(i+1)-1
              IF(j-j1.ne.j0) THEN
               IF(((icell_type(i).eq.1) .and. nbcon(j).lt.0)                                                  .or. &
                  ((icell_type(i).eq.2 .or. icell_type(i).eq.3) .and. (nbcon(j).ge.1.and.nbcon(j).le.nb_max)) .or. &
                  ((icell_type(i).eq.4) .and. nbcon(j).eq.101)) THEN
                  awalln=DMIN1(dji_nf(i1),awalln)
               ENDIF
              ENDIF
            ENDDO
            walln(i)=awalln
            wallnr(i)=1.d0/awalln
         ENDIF
        ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3d_nf(sap_nf,sa_nf,saa_nf)
!===> shuffle flux_l,flux_g,flux_d in,out  no adw
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=flux_l_nf(i1)
         y_nf_old(i)=flux_g_nf(i1)
         z_nf_old(i)=flux_d_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=flux_l_nf(i)
            buff_g(i)=flux_g_nf(i)
            buff_d(i)=flux_d_nf(i)
         ENDDO
         IF(ALLOCATED(flux_l_nf)) DEALLOCATE(flux_l_nf)
         IF(ALLOCATED(flux_g_nf)) DEALLOCATE(flux_g_nf)
         IF(ALLOCATED(flux_d_nf)) DEALLOCATE(flux_d_nf)
         ALLOCATE(flux_l_nf(nf_flux),flux_g_nf(nf_flux),flux_d_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            flux_l_nf(i)=buff_l(i)
            flux_g_nf(i)=buff_g(i)
            flux_d_nf(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(flux_l_nf,flux_g_nf,flux_d_nf, &
                                            x_nf_old,y_nf_old,z_nf_old)
!===> shuffle flux_l_o,flux_g_o,flux_d_o in,out no adw
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=flux_l_nf_o(i1)
         y_nf_old(i)=flux_g_nf_o(i1)
         z_nf_old(i)=flux_d_nf_o(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=flux_l_nf_o(i)
            buff_g(i)=flux_g_nf_o(i)
            buff_d(i)=flux_d_nf_o(i)
         ENDDO
         IF(ALLOCATED(flux_l_nf_o)) DEALLOCATE(flux_l_nf_o)
         IF(ALLOCATED(flux_g_nf_o)) DEALLOCATE(flux_g_nf_o)
         IF(ALLOCATED(flux_d_nf_o)) DEALLOCATE(flux_d_nf_o)
         ALLOCATE(flux_l_nf_o(nf_flux),flux_g_nf_o(nf_flux),flux_d_nf_o(nf_flux))
         DO i=1,nf_non+nf_mcc
            flux_l_nf_o(i)=buff_l(i)
            flux_g_nf_o(i)=buff_g(i)
            flux_d_nf_o(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(flux_l_nf_o,flux_g_nf_o,flux_d_nf_o, &
                                            x_nf_old,y_nf_old,z_nf_old)
!
!===> ecnvc_l,ecnvc_g,ecnvc_d
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=ecnvc_l_nf(i1)
         y_nf_old(i)=ecnvc_g_nf(i1)
         z_nf_old(i)=ecnvc_d_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=ecnvc_l_nf(i)
            buff_g(i)=ecnvc_g_nf(i)
            buff_d(i)=ecnvc_d_nf(i)
         ENDDO
         IF(ALLOCATED(ecnvc_l_nf)) DEALLOCATE(ecnvc_l_nf)
         IF(ALLOCATED(ecnvc_g_nf)) DEALLOCATE(ecnvc_g_nf)
         IF(ALLOCATED(ecnvc_d_nf)) DEALLOCATE(ecnvc_d_nf)
         ALLOCATE(ecnvc_l_nf(nf_flux),ecnvc_g_nf(nf_flux),ecnvc_d_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            ecnvc_l_nf(i)=buff_l(i)
            ecnvc_g_nf(i)=buff_g(i)
            ecnvc_d_nf(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(ecnvc_l_nf,ecnvc_g_nf,ecnvc_d_nf, &
                                            x_nf_old,y_nf_old,z_nf_old)
!===> liq_conv,vap_conv,drp_conv
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=liq_conv_nf(i1)
         y_nf_old(i)=vap_conv_nf(i1)
         z_nf_old(i)=drp_conv_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=liq_conv_nf(i)
            buff_g(i)=vap_conv_nf(i)
            buff_d(i)=drp_conv_nf(i)
         ENDDO
         IF(ALLOCATED(liq_conv_nf)) DEALLOCATE(liq_conv_nf)
         IF(ALLOCATED(vap_conv_nf)) DEALLOCATE(vap_conv_nf)
         IF(ALLOCATED(drp_conv_nf)) DEALLOCATE(drp_conv_nf)
         ALLOCATE(liq_conv_nf(nf_flux),vap_conv_nf(nf_flux),drp_conv_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            liq_conv_nf(i)=buff_l(i)
            vap_conv_nf(i)=buff_g(i)
            drp_conv_nf(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(liq_conv_nf,vap_conv_nf,drp_conv_nf, &
                                            x_nf_old,y_nf_old,z_nf_old)
!===> al_conv,ad_conv,void_conv
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=al_conv_nf(i1)
         y_nf_old(i)=ad_conv_nf(i1)
         z_nf_old(i)=void_conv_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=al_conv_nf(i)
            buff_g(i)=ad_conv_nf(i)
            buff_d(i)=void_conv_nf(i)
         ENDDO
         IF(ALLOCATED(al_conv_nf)) DEALLOCATE(al_conv_nf)
         IF(ALLOCATED(ad_conv_nf)) DEALLOCATE(ad_conv_nf)
         IF(ALLOCATED(void_conv_nf)) DEALLOCATE(void_conv_nf)
         ALLOCATE(al_conv_nf(nf_flux),ad_conv_nf(nf_flux),void_conv_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            al_conv_nf(i)=buff_l(i)
            ad_conv_nf(i)=buff_g(i)
            void_conv_nf(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(al_conv_nf,ad_conv_nf,void_conv_nf, &
                                            x_nf_old,y_nf_old,z_nf_old)
!===> quala_conv_nf
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=quala_conv_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=quala_conv_nf(i)
         ENDDO
         IF(ALLOCATED(quala_conv_nf)) DEALLOCATE(quala_conv_nf)
         ALLOCATE(quala_conv_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            quala_conv_nf(i)=buff_l(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy1v_nf(quala_conv_nf, &
                                            x_nf_old)
!===> lbor,dbor IF allocated in boron_convection
      IF(ALLOCATED(lbor_conv_nf)) THEN
!        copy to old
         istart=istart_nf_old(1,2)
         DO i=1,nf_nbcon_change_flux
            i1=istart+i
            x_nf_old(i)=lbor_conv_nf(i1)
         ENDDO
         IF(nf_flux.ne.nf_flux_old) THEN
            DO i=1,nf_non+nf_mcc
               buff_l(i)=lbor_conv_nf(i)
            ENDDO
            IF(ALLOCATED(lbor_conv_nf)) DEALLOCATE(lbor_conv_nf)
            ALLOCATE(lbor_conv_nf(nf_flux))
            DO i=1,nf_non+nf_mcc
               lbor_conv_nf(i)=buff_l(i)
            ENDDO
         ENDIF
         CALL vectorize_nbcon_change_copy1v_nf(lbor_conv_nf, &
                                               x_nf_old)
      ENDIF
!===> arli,argi,ardi
!     copy to old
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change_flux
         i1=istart+i
         x_nf_old(i)=arli_nf(i1)
         y_nf_old(i)=argi_nf(i1)
         z_nf_old(i)=ardi_nf(i1)
      ENDDO
      IF(nf_flux.ne.nf_flux_old) THEN
         DO i=1,nf_non+nf_mcc
            buff_l(i)=arli_nf(i)
            buff_g(i)=argi_nf(i)
            buff_d(i)=ardi_nf(i)
         ENDDO
         IF(ALLOCATED(arli_nf)) DEALLOCATE(arli_nf)
         IF(ALLOCATED(argi_nf)) DEALLOCATE(argi_nf)
         IF(ALLOCATED(ardi_nf)) DEALLOCATE(ardi_nf)
         ALLOCATE(arli_nf(nf_flux),argi_nf(nf_flux),ardi_nf(nf_flux))
         DO i=1,nf_non+nf_mcc
            arli_nf(i)=buff_l(i)
            argi_nf(i)=buff_g(i)
            ardi_nf(i)=buff_g(i)
         ENDDO
      ENDIF
      CALL vectorize_nbcon_change_copy3v_nf(arli_nf,argi_nf,ardi_nf,    &
                                            x_nf_old,y_nf_old,z_nf_old)
!
      CALL vectorize_nbcon_change_copy_ndim_nf(sv_nf)
      CALL vectorize_nbcon_change_copy_ndim_nf(svp_nf)
      CALL vectorize_nbcon_change_copy_ndim_nf(xfc_nf)
      CALL vectorize_nbcon_change_copy_ndim_nf(dxfc_nf)
      IF(ALLOCATED(perm_out)) DEALLOCATE(perm_out)
      ALLOCATE(perm_out(nf_out)) 
      CALL vectorize_nbcon_change_copy1d_nf3(perm_nf,perm_out)
!
      DEALLOCATE(nbcon_old)
!
      END SUBROUTINE nbcon_change_end
!
      SUBROUTINE vectorize_nbcon_change_index
!
!     mapping index of _adw,_inl,_out
!
      USE Zmpi          , ONLY: ia_nrhs
      USE Zzone         , ONLY: ncell_fluid
      USE Zparam        , ONLY: nin_max,nb_max
!
      USE Zbc_index     , ONLY: nbcon,nbcon_old
!
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_mcc,nf_inl,nf_out,nf_adw,nf_flux,nf_flux1,nf_flux2, &
                                nf_fluxk,nf_fluxk1,nf_fluxk2,                                         &
                                nf_nbcon_change,                                                      &
                                nf_tot,nf_tot_nb1
      USE Znum_cell     , ONLY: nb_nf,i_neigh,istart_nf,istart_nbcon_nf, &
                                istart_nb1,ia_nb,icell_nb,         &
                                istart_nf_old,index_nbcon_nf
      USE Zvec_index    , ONLY: left_nf,jneigh_nf,nbcon_nf
!
      IMPLICIT NONE
!
!.....External function
      INTEGER :: get_nf_number
!.....Local variables
      INTEGER :: nf_number,len,istart,istart1,istart2,i1,i2
      INTEGER :: i,j,j0
      INTEGER :: ii
!.....Local arrays
      INTEGER :: no(2:4),nw(2:4),nf_number_old
      INTEGER :: icount(2:4),icount_nb(2:4)
      INTEGER :: istart_nb1_old(2,-1:nb_nf)
      INTEGER :: icell_nb_old0(nf_tot_nb1)
      INTEGER :: ia_nb_old(nf_tot_nb1)
      INTEGER :: jneigh_nf_old(nf_tot)
!
      DO nf_number=-1,nb_nf
         istart_nb1_old(1,nf_number)=istart_nb1(1,nf_number)
         istart_nb1_old(2,nf_number)=istart_nb1(2,nf_number)
      ENDDO
      DO i=1,nf_tot
         jneigh_nf_old(i)=jneigh_nf(i)
      ENDDO
      DO i=1,nf_tot_nb1
         ia_nb_old(i)=ia_nb(i)
         icell_nb_old0(i)=icell_nb(i)
      ENDDO
!
      DO nf_number=2,4
         istart_nb1(2,nf_number)=0
      ENDDO
      nf_inl=0
      nf_out=0
      nf_adw=0
      DO i=1,ncell_fluid
         DO nf_number=2,4
            icount(nf_number)=0
         ENDDO
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF    (nbcon(j).gt.0 .and. nbcon(j).le.nin_max)THEN
               nf_inl=nf_inl+1
               nf_number=2
               icount(nf_number)=icount(nf_number)+1
            ELSEIF(nbcon(j).gt.nin_max .and. nbcon(j).le.nb_max)THEN
               nf_out=nf_out+1
               nf_number=3
               icount(nf_number)=icount(nf_number)+1
            ELSEIF(nbcon(j).eq.-1)THEN
               nf_adw=nf_adw+1
               nf_number=4
               icount(nf_number)=icount(nf_number)+1
            ENDIF
         ENDDO
         DO nf_number=2,4
            if(icount(nf_number).gt.0) istart_nb1(2,nf_number)=istart_nb1(2,nf_number)+1
         ENDDO
      ENDDO
!        DO nf_number=2,4
!           write(*,*) 'new_nb',nf_number,istart_nb(2,nf_number)
!        ENDDO
      nf_flux =nf_non+nf_mcc+nf_inl+nf_out
      nf_flux1=nf_non+nf_mcc+nf_inl
      nf_flux2=nf_non+nf_mcc
      nf_fluxk =nf_nonk+nf_non+nf_mcc+nf_inl+nf_out
      nf_fluxk1=nf_nonk+nf_non+nf_mcc+nf_inl
      nf_fluxk2=nf_nonk+nf_non+nf_mcc
      istart_nf(2,2)=nf_inl
      istart_nf(2,3)=nf_out
      istart_nf(2,4)=nf_adw
      istart_nf(1,3)=istart_nf(1,2)+nf_inl
      istart_nf(1,4)=istart_nf(1,3)+nf_out
      istart_nbcon_nf(3)=istart_nbcon_nf(2)+nf_inl
      istart_nbcon_nf(4)=istart_nbcon_nf(3)+nf_out
!.....Build csr pointer for scalar matrix to access non,mcc,inl
      DO i=1,ncell_fluid
         i2=0
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            IF(    (nbcon(j).ge.0 .and. nbcon(j).le.nin_max) &
               .or. nbcon(j).ge.201) THEN
               i2=i2+1
            ENDIF
         ENDDO
         ia_nrhs(i)=i2
      ENDDO
!     write(*,*) 'new====>',nf_inl,nf_out,nf_adw
      IF(.not. ALLOCATED(index_nbcon_nf)) ALLOCATE(index_nbcon_nf(nf_nbcon_change))
!     compute total size of nb,nb1 data
      nf_tot_nb1=0
      DO nf_number=-1,nb_nf
         IF(istart_nb1(2,nf_number).ne.0) THEN
            nf_tot_nb1=nf_tot_nb1+istart_nb1(2,nf_number)+1 ! add 1 per entry csr style format
         ENDIF
      ENDDO
!     write(*,*)'nf_tot_nb',nf_tot_nb1
      DEALLOCATE(icell_nb,ia_nb)
      ALLOCATE(icell_nb(nf_tot_nb1))
      ALLOCATE(ia_nb(nf_tot_nb1))
      DO nf_number=2,nb_nf
         IF(istart_nb1(2,nf_number-1).ne.0) THEN
            istart_nb1(1,nf_number)=istart_nb1(1,nf_number-1)+istart_nb1(2,nf_number-1)+1 ! add 1 per entry csr style format
         ELSE
            istart_nb1(1,nf_number)=istart_nb1(1,nf_number-1)
         ENDIF
      ENDDO
!
!     copy untouched nf_number -1,0,1
      DO i=istart_nb1(1,-1)+1,istart_nb1(1,2)
         ia_nb(i)=ia_nb_old(i)
         icell_nb(i)=icell_nb_old0(i)
      ENDDO
!     compute for nbcon 2,3,4
!
      DO nf_number=2,4
         icount_nb(nf_number)=0
         IF(istart_nb1(2,nf_number).ne.0) THEN
            istart=istart_nb1(1,nf_number)
            i1=istart+1
            ia_nb(i1)=1
         ENDIF
      ENDDO
      DO i=1,ncell_fluid
         DO nf_number=2,4
            icount(nf_number)=0
         ENDDO
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            nf_number=get_nf_number(nbcon(j))
            IF(nf_number.ge.2 .and. nf_number.le.4)THEN
               icount(nf_number)=icount(nf_number)+1
            ENDIF
         ENDDO
         DO nf_number=2,4
            IF(icount(nf_number).gt.0) THEN
              istart1=istart_nb1(1,nf_number)
              icount_nb(nf_number)=icount_nb(nf_number)+1
              ii=icount_nb(nf_number)
              i1=istart1+ii
              ia_nb(i1+1)=ia_nb(i1)+icount(nf_number)
              icell_nb(i1)=i
            ENDIF
         ENDDO
      ENDDO
      DO nf_number=2,4
         istart1=istart_nb1(1,nf_number)
         len    =istart_nb1(2,nf_number)
         IF(len.ne.0) THEN
            i1=istart1+len+1
            icell_nb(i1)=0
         ENDIF
      ENDDO
200   format(8i6)
!     copy untouched nf_number 5=>nb_nf
      DO i=istart_nb1(1,5)+1,istart_nb1(1,nb_nf)
         i1=istart_nb1_old(1,5)+i-istart_nb1(1,5)
         ia_nb(i)=ia_nb_old(i1)
         icell_nb(i)=icell_nb_old0(i1)
      ENDDO
!
      nw(2)=0
      nw(3)=0
      nw(4)=0
      no(2)=0
      no(3)=0
      no(4)=0
!     i=9940
!     if(myrank.eq.0) write(*,*) '**',(nbcon_old(j),j=i_neigh(i),i_neigh(i+1)-1)
      DO i=1,ncell_fluid
         j0=i_neigh(i)-1
         DO j=i_neigh(i),i_neigh(i+1)-1
            nf_number=get_nf_number(nbcon(j))
            IF    (nf_number.eq.2)THEN
               nf_number_old=get_nf_number(nbcon_old(j))
               istart=istart_nf(1,nf_number)
               istart2=istart_nf(1,nf_number)-(nf_non+nf_mcc)
               no(nf_number_old)=no(nf_number_old)+1
               nw(2)=nw(2)+1
               i2=istart2+nw(2)
               i1=istart+nw(2)
               nbcon_nf(i2)=nbcon(j)
               istart2=istart_nf_old(1,nf_number_old)-(nf_non+nf_mcc)
               index_nbcon_nf(i2)=istart2+no(nf_number_old)
               left_nf(i1)=i
               jneigh_nf(i1)=j-j0
            ELSEIF(nf_number.eq.3)THEN
               nf_number_old=get_nf_number(nbcon_old(j))
               istart=istart_nf(1,nf_number)
               istart2=istart_nf(1,nf_number)-(nf_non+nf_mcc)
               no(nf_number_old)=no(nf_number_old)+1
               nw(3)=nw(3)+1
               i2=istart2+nw(3)
               i1=istart+nw(3)
               nbcon_nf(i2)=nbcon(j)
               istart2=istart_nf_old(1,nf_number_old)-(nf_non+nf_mcc)
               index_nbcon_nf(i2)=istart2+no(nf_number_old)
               left_nf(i1)=i
               jneigh_nf(i1)=j-j0
            ELSEIF(nf_number.eq.4)THEN
               nf_number_old=get_nf_number(nbcon_old(j))
               istart=istart_nf(1,nf_number)
               istart2=istart_nf(1,nf_number)-(nf_non+nf_mcc)
               no(nf_number_old)=no(nf_number_old)+1
               nw(4)=nw(4)+1
               i2=istart2+nw(4)
               i1=istart+nw(4)
               nbcon_nf(i2)=nbcon(j)
               istart2=istart_nf_old(1,nf_number_old)-(nf_non+nf_mcc)
               index_nbcon_nf(i2)=istart2+no(nf_number_old)
               left_nf(i1)=i
               jneigh_nf(i1)=j-j0
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_index
!
      SUBROUTINE vectorize_nbcon_change_copy1d_nf3(x_nf,x_out)
!
      USE Znum_cell     , ONLY: index_nbcon_nf,          &
                                istart_nf,istart_nbcon_nf
      USE Zvec_param    , ONLY: nf_out,nf_nbcon_change
!
      IMPLICIT NONE
!.....Input
      REAL*8 x_nf(nf_nbcon_change)
!.....Output
      REAL*8 x_out(nf_out)
!.....Local  variables
      INTEGER i,i1,j1
      INTEGER i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2,i2
!.....Local arrays
      REAL*8 x_nf_old(nf_nbcon_change)
!
!.....Copy to old
!
      DO i=1,nf_nbcon_change
         x_nf_old(i)=x_nf(i)
      ENDDO
      j1=istart_nf(1,2)
      DO nf_number=2,4
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            i1_old=index_nbcon_nf(i2)
            x_nf(i1-j1)=x_nf_old(i1_old)
         ENDDO
      ENDDO
      nf_number=3
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i1=istart+i
         x_out(i)=x_nf(i1-j1)
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_copy1d_nf3
!
      SUBROUTINE vectorize_nbcon_change_copy2d_nf(x_nf,y_nf)
!
      USE Znum_cell     , ONLY: index_nbcon_nf,                         &
                                istart_nf,istart_nbcon_nf,istart_nf_old
      USE Zvec_param    , ONLY: nf_nbcon_change,nf_tot
!
      IMPLICIT NONE
!.....Output
      REAL*8 x_nf(nf_tot)
      REAL*8 y_nf(nf_tot)
!.....Local variables
      INTEGER :: i,i1,i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2,i2
!.....Local arrays
      REAL*8 x_nf_old(nf_nbcon_change),y_nf_old(nf_nbcon_change)
!
!.....Copy to old
!
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change
         i1=istart+i
         x_nf_old(i)=x_nf(i1)
         y_nf_old(i)=y_nf(i1)
      ENDDO
      DO nf_number=2,4
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            i1_old=index_nbcon_nf(i2)
            x_nf(i1)=x_nf_old(i1_old)
            y_nf(i1)=y_nf_old(i1_old)
         ENDDO
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_copy2d_nf
!
      SUBROUTINE vectorize_nbcon_change_copy3d_nf(x_nf,y_nf,z_nf)
!
      USE Znum_cell     , ONLY: index_nbcon_nf,                         &
                                istart_nf,istart_nbcon_nf,istart_nf_old
      USE Zvec_param    , ONLY: nf_nbcon_change,nf_tot
!
      IMPLICIT NONE
!.....Output
      REAL*8 x_nf(nf_tot),y_nf(nf_tot),z_nf(nf_tot)
!.....Local variables
      INTEGER :: i,i1,i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2,i2
!.....Local arrays
      REAL*8 x_nf_old(nf_nbcon_change),y_nf_old(nf_nbcon_change),z_nf_old(nf_nbcon_change)
!
!.....Copy to old
!
      istart=istart_nf_old(1,2)
      DO i=1,nf_nbcon_change
         i1=istart+i
         x_nf_old(i)=x_nf(i1)
         y_nf_old(i)=y_nf(i1)
         z_nf_old(i)=z_nf(i1)
      ENDDO
      DO nf_number=2,4
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            i1_old=index_nbcon_nf(i2)
            x_nf(i1)=x_nf_old(i1_old)
            y_nf(i1)=y_nf_old(i1_old)
            z_nf(i1)=z_nf_old(i1_old)
         ENDDO
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_copy3d_nf
!
      SUBROUTINE vectorize_nbcon_change_copy_ndim_nf(x_nf)
!
      USE Zparam      , ONLY: ndim
      USE Znum_cell     , ONLY: index_nbcon_nf,                         &
                                istart_nf,istart_nbcon_nf,istart_nf_old
      USE Zvec_param    , ONLY: nf_nbcon_change,nf_tot
!
      IMPLICIT NONE
!.....Output
      REAL*8 x_nf(nf_tot,ndim)
!.....Local variables
      INTEGER :: i,i1,i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2,i2
!.....Local arrays
      REAL*8 x_nf_old(nf_nbcon_change,ndim)
!
!.....Copy to old
!
      istart=istart_nf_old(1,2)
      IF(ndim.eq.2) THEN
         DO i=1,nf_nbcon_change
            i1=istart+i
            x_nf_old(i,1)=x_nf(i1,1)
            x_nf_old(i,2)=x_nf(i1,2)
         ENDDO
         DO nf_number=2,4
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            istart2=istart_nbcon_nf(nf_number)
            DO i=1,isize
               i1=istart+i
               i2=istart2+i
               i1_old=index_nbcon_nf(i2)
               x_nf(i1,1)=x_nf_old(i1_old,1)
               x_nf(i1,2)=x_nf_old(i1_old,2)
            ENDDO
         ENDDO
      ELSE
         DO i=1,nf_nbcon_change
            i1=istart+i
            x_nf_old(i,1)=x_nf(i1,1)
            x_nf_old(i,2)=x_nf(i1,2)
            x_nf_old(i,3)=x_nf(i1,3)
         ENDDO
         DO nf_number=2,4
            istart=istart_nf(1,nf_number)
            isize =istart_nf(2,nf_number)
            istart2=istart_nbcon_nf(nf_number)
            DO i=1,isize
               i1=istart+i
               i2=istart2+i
               i1_old=index_nbcon_nf(i2)
               x_nf(i1,1)=x_nf_old(i1_old,1)
               x_nf(i1,2)=x_nf_old(i1_old,2)
               x_nf(i1,3)=x_nf_old(i1_old,3)
            ENDDO
         ENDDO
      ENDIF
!
      END SUBROUTINE vectorize_nbcon_change_copy_ndim_nf
!
      SUBROUTINE vectorize_nbcon_change_copy1v_nf(x_nf,     &
                                                  x_nf_old)
!
      USE Znum_cell     , ONLY: index_nbcon_nf,                         &
                                istart_nf,istart_nbcon_nf,istart_nf_old
      USE Zvec_param    , ONLY: nf_non,nf_mcc,nf_nbcon_change,nf_tot
!
      IMPLICIT NONE
!.....Input
      REAL*8 x_nf_old(nf_nbcon_change)
!.....Output
      REAL*8 x_nf(nf_tot)
!.....Local variables
      INTEGER :: i,i1,i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2_old,istart2,i2
!
      istart2_old=istart_nf_old(1,4)-(nf_non+nf_mcc)
      DO nf_number=2,3
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            i1_old=index_nbcon_nf(i2)
            if(i1_old.le.istart2_old) THEN
               x_nf(i1)=x_nf_old(i1_old)
            ELSE
               x_nf(i1)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_copy1v_nf
!
      SUBROUTINE vectorize_nbcon_change_copy3v_nf(x_nf,y_nf,z_nf,             &
                                                  x_nf_old,y_nf_old,z_nf_old)
!
      USE Znum_cell     , ONLY: index_nbcon_nf,                         &
                                istart_nf,istart_nbcon_nf,istart_nf_old
      USE Zvec_index    , ONLY: left_nf,jneigh_nf
      USE Zvec_param    , ONLY: nf_non,nf_mcc,nf_nbcon_change,nf_tot
!
      IMPLICIT NONE
!.....Input
      REAL*8 x_nf_old(nf_nbcon_change),y_nf_old(nf_nbcon_change),z_nf_old(nf_nbcon_change)
!.....Output
      REAL*8 x_nf(nf_tot),y_nf(nf_tot),z_nf(nf_tot)
!.....Local variables
      INTEGER :: i,i1,ii,jj,i1_old
      INTEGER :: nf_number,isize,istart
      INTEGER :: istart2_old,istart2,i2
!
      istart2_old=istart_nf_old(1,4)-(nf_non+nf_mcc)
      DO nf_number=2,3
         istart=istart_nf(1,nf_number)
         isize =istart_nf(2,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         DO i=1,isize
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            jj=jneigh_nf(i1)
            i1_old=index_nbcon_nf(i2)
            if(i1_old.le.istart2_old) THEN
               x_nf(i1)=x_nf_old(i1_old)
               y_nf(i1)=y_nf_old(i1_old)
               z_nf(i1)=z_nf_old(i1_old)
            ELSE
               x_nf(i1)=0.d0
               y_nf(i1)=0.d0
               z_nf(i1)=0.d0
            ENDIF
         ENDDO
      ENDDO
!
      END SUBROUTINE vectorize_nbcon_change_copy3v_nf
