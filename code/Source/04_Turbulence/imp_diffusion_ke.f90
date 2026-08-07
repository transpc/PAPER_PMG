!
      SUBROUTINE imp_diffusion_ke(flag,bm,cm,turb_o,keb,diff_ke,turb_new,ar,flux_nf,ab,rb, &
                                  eps_imp,max_iter,iml)
!
!     This routine calculates implicit momentum diffusion
!
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp,au,ju_a,ia_a
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim,nin_max
      USE Zvec_param   , ONLY: nf_flux,nf_nonk,nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf,                 &
                               nf_number_nb,lens,nf_number_id,istart_nfs, &
                               right_nb_k
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zconst1      , ONLY: lowreynolds,iturb,turb_phase
      USE Zconst2      , ONLY: dt
      USE Zcoord3      , ONLY: volr
      USE Zgradoption  , ONLY: non_orth_turb      
      USE Zimplicit    , ONLY: imp_ke_diff,imp_ke_conv
      USE Zndforce     , ONLY: d_bfc
      USE Z2nd_order   , ONLY: turb_conv_2nd
      USE Zvec_geo     , ONLY: saa_nf,djia_nf,sap_nf,   &
                               fac_non,fac1_non,dnj_non
!
      IMPLICIT NONE
!
!.....Input 
      INTEGER :: flag,max_iter,iml
      REAL(8) :: eps_imp
      REAL(8),DIMENSION(ncell_fluid) :: bm,cm
      REAL(8),DIMENSION(ncell_fp) :: ar,turb_o,diff_ke
      REAL(8),DIMENSION(nf_flux) :: flux_nf
      REAL(8),DIMENSION(nin_max) :: ab,rb,keb
!.....Output
      REAL(8),DIMENSION(ncell_fp) :: turb_new
!.....Local variables
      INTEGER :: i,k
      INTEGER :: ii,jj,kk
      INTEGER :: nv,nf_number,istart0,istart,len,istart2,i0,i1,i2
      REAL(8) :: dkeix,dkeiy,dkeiz,dkej
      REAL(8) :: ar2,arv2
      REAL(8) :: vt,cf_ke
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: diag_ke,src_ke
      REAL(8),DIMENSION(ncell_fp,ndim) :: dkedx
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: off_diag_ke_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_ke_non_k
      REAL(8),DIMENSION(nf_non) :: cv_ke_diag_non
      REAL(8),DIMENSION(nf_inl) :: cv_kekeb_inl
      REAL(8),DIMENSION(nf_non+nf_inl) :: dkej_nf
      REAL(8),DIMENSION(nf_nonk+nf_non+nf_inl) :: cv_ke_diag_nf
      REAL(8),DIMENSION(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw) :: cf_ke_diag_nf
!
!.....Calculate ke/dp gradient at cell center for non-orthogonal grid
!
      IF(non_orth_turb.eq.1) THEN
         CALL grad_scalar(turb_o,dkedx,ncell_fp)
      ELSEIF(non_orth_turb.eq.2) THEN
         CALL grad_scalar(turb_o,dkedx,ncell_fp)  ! will be replaced with frink method
      ENDIF
      IF(np.gt.1) CALL communicate_2d(dkedx)
!
      DO i=1,ncell_fluid
         diag_ke(i)=0.d0
         src_ke(i) =0.d0
      ENDDO
!
!.....Implicit Diffusion
!
      IF(imp_ke_diff.eq.1) THEN            
!
!.....Build summation info for non,inl,adw,fsw,ctw,chw
!
         IF( iturb.eq.1                          .or. &
            (iturb.eq.2 .and. lowreynolds.eq.1)) THEN
            nf_number_nb=5
            nf_number_id(0)=0
            nf_number_id(1)=2
            nf_number_id(2)=4
            nf_number_id(3)=5
            nf_number_id(4)=6
            nf_number_id(5)=7
            istart_nfs(0)=0
            istart_nfs(1)=istart_nfs(0)+nf_non
            istart_nfs(2)=istart_nfs(1)+nf_inl
            istart_nfs(3)=istart_nfs(2)+nf_adw
            istart_nfs(4)=istart_nfs(3)+nf_fsw
            istart_nfs(5)=istart_nfs(4)+nf_ctw
            lens         =istart_nfs(5)+nf_chw
         ELSE
            nf_number_nb=1
            nf_number_id(0)=0
            nf_number_id(1)=2
            istart_nfs(0)=0
            istart_nfs(1)=istart_nfs(0)+nf_non
            lens         =istart_nfs(1)+nf_inl
         ENDIF
!
!........Cells non
!         
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            IF(non_orth_turb.gt.0) THEN
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  cf_ke=(fac1_non(i)*diff_ke(ii) +fac_non(i)*diff_ke(kk))*saa_nf(i1)
                  dkeix= fac1_non(i)*dkedx(ii,1)+fac_non(i)*dkedx(kk,1)
                  dkeiy= fac1_non(i)*dkedx(ii,2)+fac_non(i)*dkedx(kk,2)
                  dkej=dkeix*dnj_non(i,1)+dkeiy*dnj_non(i,2)
                  dkej_nf(i0)=cf_ke*dkej
                  cf_ke_diag_nf(i0)=cf_ke/djia_nf(i1)
                  off_diag_ke_non_i(i)=cf_ke_diag_nf(i0)
               ENDDO 
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  off_diag_ke_non_k(i)=-cf_ke_diag_nf(k)
               ENDDO
            ELSE
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  cf_ke=(fac1_non(i)*diff_ke(ii)+fac_non(i)*diff_ke(kk))*saa_nf(i1)
                  cf_ke_diag_nf(i0)=cf_ke/djia_nf(i1)
                  off_diag_ke_non_i(i)=-cf_ke_diag_nf(i0)*volr(ii)
               ENDDO 
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  ii=right_non(k)
                  off_diag_ke_non_k(i)=-cf_ke_diag_nf(k)*volr(ii)
               ENDDO
            ENDIF 
         ELSE
            IF(non_orth_turb.gt.0) THEN
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  cf_ke=(fac1_non(i)*diff_ke(ii) +fac_non(i)*diff_ke(kk))*saa_nf(i1)
                  dkeix= fac1_non(i)*dkedx(ii,1)+fac_non(i)*dkedx(kk,1)
                  dkeiy= fac1_non(i)*dkedx(ii,2)+fac_non(i)*dkedx(kk,2)
                  dkeiz= fac1_non(i)*dkedx(ii,3)+fac_non(i)*dkedx(kk,3)
                  dkej=dkeix*dnj_non(i,1)+dkeiy*dnj_non(i,2)+dkeiz*dnj_non(i,3)
                  dkej_nf(i0)=cf_ke*dkej
                  cf_ke_diag_nf(i0)=cf_ke/djia_nf(i1)
                  off_diag_ke_non_i(i)=-cf_ke_diag_nf(i0)*volr(ii)
               ENDDO 
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  ii=right_non(k)
                  off_diag_ke_non_k(i)=-cf_ke_diag_nf(k)*volr(ii)
               ENDDO
            ELSE
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  kk=right_non(i)
                  cf_ke=(fac1_non(i)*diff_ke(ii)+fac_non(i)*diff_ke(kk))*saa_nf(i1)
                  cf_ke_diag_nf(i0)=cf_ke/djia_nf(i1)
                  off_diag_ke_non_i(i)=-cf_ke_diag_nf(i0)*volr(ii)
               ENDDO 
               DO i=1,nf_nonk
                  k=right_nb_k(i)
                  ii=right_non(k)
                  off_diag_ke_non_k(i)=-cf_ke_diag_nf(k)*volr(ii)
               ENDDO
            ENDIF 
         ENDIF 
!
!........Cells inl
!
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
          DO i=1,len  
             i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            cf_ke=diff_ke(ii)*sap_nf(i1)
            cf_ke_diag_nf(i0)=cf_ke
            dkej_nf(i0)=cf_ke*keb(k)
         ENDDO
!
!........The rest
!
         IF(iturb.eq.1)THEN
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               IF(iml.eq.1)THEN
                  DO i=1,len  
                     i0=istart0+i
                     i1=istart+i
                     cf_ke_diag_nf(i0)=0.d0
                  ENDDO 
               ELSEIF(iml.eq.2)THEN
                  IF(turb_phase.eq.2)THEN
                     DO i=1,len  
                        i0=istart0+i
                        i1=istart+i
                        ii=left_nf(i1)
                        cf_ke=60.d0*cell%lviscosl(ii)/cell%rhol(ii)/(0.075d0*d_bfc(ii)*d_bfc(ii))
                        cf_ke_diag_nf(i0)=cf_ke*diff_ke(ii)*sap_nf(i1)
                     ENDDO 
                  ELSEIF(turb_phase.eq.1)THEN
                     DO i=1,len  
                        i0=istart0+i
                        i1=istart+i
                        ii=left_nf(i1)
                        cf_ke=60.d0*cell%lviscosg(ii)/cell%rhog(ii)/(0.075d0*d_bfc(ii)*d_bfc(ii))
                        cf_ke_diag_nf(i0)=cf_ke*diff_ke(ii)*sap_nf(i1)
                     ENDDO 
                  ENDIF
               ENDIF
            ENDDO 
         ELSEIF(iturb.eq.2 .and. lowreynolds.eq.1) THEN
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len  
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  cf_ke_diag_nf(i0)=diff_ke(ii)*sap_nf(i1)
               ENDDO 
            ENDDO 
         ENDIF
         CALL sum_nf(1,1,                   &
                     cf_ke_diag_nf,diag_ke)
         IF(non_orth_turb.gt.0) THEN
            nf_number_nb=1
            nf_number_id(1)=0
            istart_nfs(1)=0
            lens=istart_nfs(1)+nf_non
            CALL sum_nf(1,1,            &
                        dkej_nf,src_ke)
         ENDIF
!
         nf_number_nb=0
         nf_number_id(0)=2
         istart_nfs(0)=nf_non
         lens         =istart_nfs(0)+nf_inl
         CALL sum_nf(1,1,            &
                     dkej_nf,src_ke)
!
         IF(      iturb.eq.1 &
            .or. (iturb.eq.2 .and. lowreynolds.eq.1)) THEN
            nf_number_nb=3
            nf_number_id(0)=4
            nf_number_id(1)=5
            nf_number_id(2)=6
            nf_number_id(3)=7
            istart_nfs(0)=nf_non+nf_inl
            istart_nfs(1)=istart_nfs(0)+nf_adw
            istart_nfs(2)=istart_nfs(1)+nf_fsw
            istart_nfs(3)=istart_nfs(2)+nf_ctw
            lens         =istart_nfs(3)+nf_chw
            CALL sum_nf(1,1,                  &
                        cf_ke_diag_nf,src_ke)
         ENDIF
      ENDIF 
!           
!.....Implicit Convection
!
      IF(imp_ke_conv.eq.1) THEN
!
!........Build summation info for non,inl
!
         nf_number_nb=1
         nf_number_id(-1)=-1
         nf_number_id(0)=0
         nf_number_id(1)=2
         istart_nfs(-1)=0
         istart_nfs(0)=istart_nfs(-1)+nf_nonk
         istart_nfs(1)=istart_nfs(0) +nf_non
         lens         =istart_nfs(1) +nf_inl
!
!........Cells non
!
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            ii=left_nf(i1)
            kk=right_non(i)
            IF(flux_nf(i1).lt.0.d0) THEN
               arv2=flux_nf(i1)*ar(kk)
               cv_ke_diag_non(i)=arv2
               cv_ke_diag_nf(i0)=-arv2
               off_diag_ke_non_i(i)=off_diag_ke_non_i(i)+arv2*volr(ii)
            ELSE
               arv2=flux_nf(i1)*ar(ii)
               cv_ke_diag_non(i)=arv2
               cv_ke_diag_nf(i0)=0.d0
            ENDIF
         ENDDO
!
         nv=-1
         nf_number=nf_number_id(nv)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            k=right_nb_k(i)
            ii=right_non(k)
            kk=left_nf(k)
            IF(flux_nf(k).lt.0.d0) THEN
               cv_ke_diag_nf(i)=0.d0
            ELSE
               arv2=flux_nf(k)*ar(kk)
               cv_ke_diag_nf(i)=arv2
               off_diag_ke_non_k(i)=off_diag_ke_non_k(i)-arv2*volr(ii)
            ENDIF
         ENDDO
!
!........Cells inl
!                          
         nv=1
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            ar2=ab(k)*rb(k)
            IF(flux_nf(i1).lt.0.d0) THEN
               arv2=flux_nf(i1)*ar2
               cv_ke_diag_nf(i0)=-arv2
            ELSE
               cv_ke_diag_nf(i0)=0.d0
            ENDIF 
         ENDDO
!
         CALL sum_nf(1,0,                   &
                     cv_ke_diag_nf,diag_ke)
!
!........Build summation info for inl
!
         nf_number_nb=0
         nf_number_id(0)=2
         istart_nfs(0)=0
         lens         =istart_nfs(0)+nf_inl
!
!........Cells inl
!                          
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         istart2=istart_nbcon_nf(nf_number)
         len   =istart_nf(2,nf_number)
         DO i=1,len  
            i0=istart0+i
            i1=istart+i
            i2=istart2+i
            ii=left_nf(i1)
            k=nbcon_nf(i2)
            ar2=ab(k)*rb(k)
            IF(flux_nf(i1).lt.0.d0) THEN
               arv2=flux_nf(i1)*ar2
               cv_kekeb_inl(i0)=-arv2*keb(k)
            ELSE
               cv_kekeb_inl(i0)=0.d0
            ENDIF 
         ENDDO
!
         CALL sum_nf(1,-1,                &
                     cv_kekeb_inl,src_ke)
      ENDIF
!
      vt=1.d0/dt
      DO i=1,ncell_fluid
         diag_ke(i)=ar(i)*vt-cm(i)+diag_ke(i)*volr(i)
         src_ke(i) =ar(i)*vt*turb_o(i)+bm(i)+src_ke(i)*volr(i)
      ENDDO
!
!.....2nd order convection
!
      IF(turb_conv_2nd.gt.0) CALL turb_2nd_conv_imp(turb_o,flux_nf,src_ke,cv_ke_diag_non)
      
!
!.....Build directly solverCSR  array here
!
      IF(flag.eq.2) THEN
         DO i=1,ncell_fluid
            IF(cell%alphag(i).lt.1.d-15) THEN
               diag_ke(i)=1.d0
               src_ke(i)=turb_o(i)
            ENDIF
         ENDDO
      ELSEIF(flag.eq.1) THEN
         DO i=1,ncell_fluid
            IF(cell%alphal(i).lt.1.d-15) THEN
               diag_ke(i)=1.d0
               src_ke(i)=turb_o(i)
            ENDIF
         ENDDO
      ENDIF
      CALL csr_build_a(diag_ke,off_diag_ke_non_i,off_diag_ke_non_k)
      IF(flag.eq.2) THEN
         DO i=1,ncell_fluid
            IF(cell%alphag(i).lt.1.d-15) THEN
               DO jj=ia_a(i),ju_a(i)-1
                  au(jj)=0.d0
               ENDDO
               DO jj=ju_a(i)+1,ia_a(i+1)-1
                  au(jj)=0.d0
               ENDDO
            ENDIF
         ENDDO
      ELSEIF(flag.eq.1) THEN
         DO i=1,ncell_fluid
            IF(cell%alphal(i).lt.1.d-15) THEN
               DO jj=ia_a(i),ju_a(i)-1
                  au(jj)=0.d0
               ENDDO
               DO jj=ju_a(i)+1,ia_a(i+1)-1
                  au(jj)=0.d0
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      CALL csr_cg_solvers_scalar(diag_ke,src_ke,turb_new,eps_imp,max_iter)
!
      END SUBROUTINE imp_diffusion_ke
