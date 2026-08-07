!
      SUBROUTINE turb_ke_diffusion(ke,dp,keb,dpb,diff_ke,diff_dp,diff_k,diff_e)
!
!     This routine calculates diffusion terms of k-epsilon transport equation.
!       
      USE Zinterface
      USE VOL_DATA     , ONLY: cell
      USE Zmpi         , ONLY: ncell_fp
      USE Zzone        , ONLY: ncell_fluid
      USE Zcore        , ONLY: np
      USE Zparam       , ONLY: ndim,nin_max
      USE Zvec_param   , ONLY: nf_non,nf_inl,nf_adw,nf_fsw,nf_ctw,nf_chw
      USE Znum_cell    , ONLY: istart_nf,istart_nbcon_nf, &
                                nf_number_nb,lens,nf_number_id,istart_nfs
      USE Zvec_index   , ONLY: left_nf,right_non,nbcon_nf
      USE Zconst1      , ONLY: lowreynolds,iturb,turb_phase
      USE Zgradoption  , ONLY: non_orth_turb
      USE Zndforce     , ONLY: d_bfc
      USE Zvec_geo     , ONLY: sa_nf,sap_nf,            &
                               fac_non,fac1_non,dnj_non
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: keb(nin_max),dpb(nin_max)
      REAL(8) :: ke(ncell_fp),dp(ncell_fp)
      REAL(8) :: diff_ke(ncell_fp),diff_dp(ncell_fp)
!.....Output
      REAL(8) :: diff_k(ncell_fluid),diff_e(ncell_fluid)
!     local variables
      INTEGER :: i,k
      INTEGER :: ii,kk
      INTEGER :: nv,nf_number,len,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: ke1,ke2,dp1,dp2,diff_kei,diff_dpi
      REAL(8) :: dkej,ddpj
      REAL(8) :: dkeix,dkeiy,dkeiz
      REAL(8) :: ddpix,ddpiy,ddpiz
!.....Local arrays 
      REAL(8) :: dkedx(ncell_fp,ndim),ddpdx(ncell_fp,ndim)
!.....Local vector arrays 
      REAL(8) :: diff_k_nf(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw)
      REAL(8) :: diff_e_nf(nf_non+nf_inl+nf_adw+nf_fsw+nf_ctw+nf_chw)
!
!.....Calculate ke/dp gradient at cell center for non-orthogonal grid
!
      IF(non_orth_turb.eq.1) THEN
         CALL grad_scalar(ke,dkedx,ncell_fp)
         IF(np.gt.1) CALL communicate_2d(dkedx)
         CALL grad_scalar(dp,ddpdx,ncell_fp)
         IF(np.gt.1) CALL communicate_2d(ddpdx)
      ELSEIF(non_orth_turb.eq.2) THEN
         CALL grad_scalar(ke,dkedx,ncell_fp)  ! will be replaced with Frink method
         CALL grad_scalar(dp,ddpdx,ncell_fp)
         IF(np.gt.1) CALL communicate_2d(dkedx, &
                                         ddpdx)
      ENDIF
!
!.....Build summation info for non,inl,fsw,ctw,chw
!
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
!
         ke1=ke(ii)
         dp1=dp(ii)
         ke2=ke(kk)
         dp2=dp(kk)
         diff_kei=fac1_non(i)*diff_ke(ii)+fac_non(i)*diff_ke(kk)
         diff_dpi=fac1_non(i)*diff_dp(ii)+fac_non(i)*diff_dp(kk)
         diff_k_nf(i0)=diff_kei*(ke2-ke1)*sap_nf(i1)
         diff_e_nf(i0)=diff_dpi*(dp2-dp1)*sap_nf(i1)
      ENDDO
!
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart =istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      len    =istart_nf(2,nf_number)
      DO i=1,len
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=nbcon_nf(i2)
         ke1=ke(ii)
         dp1=dp(ii)
         ke2=keb(k)
         dp2=dpb(k)
         diff_kei=diff_ke(ii)
         diff_dpi=diff_dp(ii)
         diff_k_nf(i0)=diff_kei*(ke2-ke1)*sap_nf(i1)
         diff_e_nf(i0)=diff_dpi*(dp2-dp1)*sap_nf(i1)
      ENDDO
!
      IF(iturb.eq.1)THEN
         IF(turb_phase.eq.2)THEN 
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  diff_kei=diff_ke(ii)
                  diff_dpi=diff_dp(ii)          
                  ke1=ke(ii)
                  dp1=dp(ii) 
                  dp2=60.d0*cell%lviscosl(ii)/cell%rhol(ii)/(0.075d0*d_bfc(ii)*d_bfc(ii))
                  diff_k_nf(i0)=-diff_kei*ke1*sap_nf(i1)
                  diff_e_nf(i0)= diff_dpi*(dp2-dp1)*sap_nf(i1)
               ENDDO
            ENDDO
         ELSEIF(turb_phase.eq.1)THEN
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  diff_kei=diff_ke(ii)
                  diff_dpi=diff_dp(ii)          
                  ke1=ke(ii)
                  dp1=dp(ii) 
                  dp2=60.d0*cell%lviscosg(ii)/cell%rhog(ii)/(0.075d0*d_bfc(ii)*d_bfc(ii))
                  diff_k_nf(i0)=-diff_kei*ke1*sap_nf(i1)
                  diff_e_nf(i0)= diff_dpi*(dp2-dp1)*sap_nf(i1)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         IF(lowreynolds.ge.1)THEN
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  ii=left_nf(i1)
                  diff_kei=diff_ke(ii)
                  diff_dpi=diff_dp(ii)
                  ke1=ke(ii)
                  dp1=dp(ii) 
                  diff_k_nf(i0)=-diff_kei*ke1*sap_nf(i1)
                  diff_e_nf(i0)=-diff_dpi*dp1*sap_nf(i1)
               ENDDO
            ENDDO
         ELSE
            DO nv=2,5
               nf_number=nf_number_id(nv)
               istart0=istart_nfs(nv)
               istart=istart_nf(1,nf_number)
               len   =istart_nf(2,nf_number)
               DO i=1,len
                  i0=istart0+i
                  i1=istart+i
                  diff_k_nf(i0)=0.d0
                  diff_e_nf(i0)=0.d0
               ENDDO
            ENDDO
         ENDIF
      ENDIF
!               
      IF(non_orth_turb.gt.0)THEN
         nv=0
         nf_number=nf_number_id(nv)
         istart0=istart_nfs(nv)
         istart=istart_nf(1,nf_number)
         len   =istart_nf(2,nf_number)
         IF(ndim.eq.2) THEN
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
               dkeix=fac1_non(i)*dkedx(ii,1)+fac_non(i)*dkedx(kk,1)
               dkeiy=fac1_non(i)*dkedx(ii,2)+fac_non(i)*dkedx(kk,2)
!
               ddpix=fac1_non(i)*ddpdx(ii,1)+fac_non(i)*ddpdx(kk,1)
               ddpiy=fac1_non(i)*ddpdx(ii,2)+fac_non(i)*ddpdx(kk,2)
!
               dkej=dkeix*dnj_non(i,1)+dkeiy*dnj_non(i,2)
               ddpj=ddpix*dnj_non(i,1)+ddpiy*dnj_non(i,2)
!               
               diff_k_nf(i0)=diff_k_nf(i0)+diff_kei*dkej*sa_nf(i1)
               diff_e_nf(i0)=diff_e_nf(i0)+diff_dpi*ddpj*sa_nf(i1)
            ENDDO
         ELSE
            DO i=1,len
               i0=istart0+i
               i1=istart+i
               ii=left_nf(i1)
               kk=right_non(i)
!            
               dkeix=fac1_non(i)*dkedx(ii,1)+fac_non(i)*dkedx(kk,1)
               dkeiy=fac1_non(i)*dkedx(ii,2)+fac_non(i)*dkedx(kk,2)
               dkeiz=fac1_non(i)*dkedx(ii,3)+fac_non(i)*dkedx(kk,3)
!
               ddpix=fac1_non(i)*ddpdx(ii,1)+fac_non(i)*ddpdx(kk,1)
               ddpiy=fac1_non(i)*ddpdx(ii,2)+fac_non(i)*ddpdx(kk,2)
               ddpiz=fac1_non(i)*ddpdx(ii,3)+fac_non(i)*ddpdx(kk,3)
!
               dkej=dkeix*dnj_non(i,1)+dkeiy*dnj_non(i,2)+dkeiz*dnj_non(i,3)
               ddpj=ddpix*dnj_non(i,1)+ddpiy*dnj_non(i,2)+ddpiz*dnj_non(i,3)
!               
               diff_k_nf(i0)=diff_k_nf(i0)+diff_kei*dkej*sa_nf(i1)
               diff_e_nf(i0)=diff_e_nf(i0)+diff_dpi*ddpj*sa_nf(i1)
            ENDDO
         ENDIF
      ENDIF
!
      CALL sum_nf(0,-1,             &
                  diff_k_nf,diff_k, &
                  diff_e_nf,diff_e)
!
      END  SUBROUTINE turb_ke_diffusion
