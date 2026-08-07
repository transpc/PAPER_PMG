!
      SUBROUTINE radiation_model
!
!     This routine calculates incident radiation and radiation heat flux
!     using P-1 radiation model.
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE SOLID_DATA    , ONLY: solid
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Zvec_param    , ONLY: nf_nonk,nf_non,nf_ctw,nf_fsw
      USE Zvec_index    , ONLY: left_nf,right_non,right_fsw
      USE Znum_cell     , ONLY: istart_nf,istart_nbcon_nf, &
                                nf_number_nb,lens,nf_number_id,istart_nfs, &
                                right_nb_k
      USE Zcoord3       , ONLY: volr
      USE Zturb         , ONLY: wallnr
      USE Zmodel        , ONLY: abs_coeff,wall_emiss,max_iter_rad,eps_imp_rad,rad_source,qrad
      USE Zvec_geo      , ONLY: fac_non,fac1_non,sap_nf,f1,f0,saa_nf
      USE Zvec_major    , ONLY: rad_ir_nbcon
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,k,ii,kk
      INTEGER :: nv,nf_number,isize,istart0,istart,istart2,i0,i1,i2
      REAL(8) :: scat_coeff,las_coeff,SIGMA_SBC,Ew,tw
      REAL(8) :: cf_ir
      REAL(8) :: rgamm,ir1,ir2
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: diag_ir,src_ir,bm,cm
      REAL(8),DIMENSION(ncell_fp) :: rad_gamm,rad_ir
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: off_diag_ir_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_ir_non_k
      REAL(8),DIMENSION(nf_ctw) :: twall_ctw,rad_irctw
      REAL(8),DIMENSION(nf_fsw) :: rad_irfsw
      REAL(8),DIMENSION(nf_non+nf_fsw+nf_ctw) :: cf_ir_nf,rad_diff_nf,qrad_nf
      REAL(8),DIMENSION(nf_fsw+nf_ctw) :: cf_irb
!
      scat_coeff=1.0d0-abs_coeff  ! scattering coefficient (scatter=1-absorbtion)   [1/m]
      las_coeff=0.0d0             ! linear-anisotropic scattering coefficient
      SIGMA_SBC=5.67d-8           ! Stefan-Boltzmann constant [W/m2-K4]
      Ew=wall_emiss/(4.0d0-2.0d0*wall_emiss)
!
      cm(:)=1.0d0
      bm(:)=0.0d0
      rad_gamm(:)=0.0d0
      DO i=1,ncell_fluid  
!         IF(nzone(i).ne.izone)CYCLE    !input izone when a user want to activate 'radiation model' in the specific zone
         rad_gamm(i)=3.0d0*abs_coeff+(3.0d0-las_coeff)*scat_coeff
         rad_gamm(i)=1.0d0/rad_gamm(i)                      ! [m]
         bm(i)=-4.0d0*abs_coeff*SIGMA_SBC*(cell%tg(i)**4)   ! right term (source)
         cm(i)=-abs_coeff                                   ! left term  (diagonal)
      ENDDO
!
      IF(np.gt.1) CALL communicate_1d(rad_gamm)
!
!.....Build summation of diagonals for non, fsw, ctw
!
      nf_number_nb=2
      nf_number_id(0)=0
      nf_number_id(1)=5
      nf_number_id(2)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_fsw
      lens         =istart_nfs(2)+nf_ctw
!
!.....computing cell ; G transport equation
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         
         cf_ir=0.0d0
!         IF(nzone(i).ne.izone)CYCLE    !input izone when a user want to activate 'radiation model' in the specific zone
         kk=right_non(i)
         cf_ir=(fac1_non(i)*rad_gamm(ii)+fac_non(i)*rad_gamm(kk))*sap_nf(i1)
         cf_ir_nf(i0)=cf_ir
         off_diag_ir_non_i(i)=cf_ir*volr(ii)
!        off_diag_ir_non_k(i)=cf_ir*volr(kk)
      ENDDO
!
      nf_number=-1
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         k=right_nb_k(i)
         ii=right_non(k)
         kk=left_nf(k)
         cf_ir=(fac1_non(k)*rad_gamm(kk)+fac_non(k)*rad_gamm(ii))*sap_nf(k)
         off_diag_ir_non_k(i)=cf_ir*volr(ii)
      ENDDO
!
!.....Fluid-Solid interface ; G transport equation
!      
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         cf_ir=rad_gamm(ii)*sap_nf(i1)
         cf_ir_nf(i0)=cf_ir
      ENDDO
!
!.....constant temperature wall ; G transport equation
!      
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         cf_ir=rad_gamm(ii)*sap_nf(i1)
         cf_ir_nf(i0)=cf_ir
      ENDDO
!
      CALL sum_nf(0,1,              &
                  cf_ir_nf,diag_ir)
!
!.....Build summation of source term for fsw, ctw
!
      nf_number_nb=1
      nf_number_id(0)=5
      nf_number_id(1)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_fsw  
      lens         =istart_nfs(1)+nf_ctw
!
!.....Fluid-Solid interface ; G transport equation
!     
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         cf_ir=rad_gamm(ii)*sap_nf(i1)
         rad_irfsw(i)=rad_ir_nbcon(i2)
         cf_irb(i0)=cf_ir*rad_irfsw(i)
      ENDDO
!
!.....constant temperature wall ; G transport equation
!      
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         cf_ir=rad_gamm(ii)*sap_nf(i1)
         rad_irctw(i)=rad_ir_nbcon(i2)
         cf_irb(i0)=cf_ir*rad_irctw(i)
      ENDDO   
      CALL sum_nf(0,-1,          &
                  cf_irb,src_ir)
!
      DO i=1,ncell_fluid
         diag_ir(i)=cm(i)-diag_ir(i)*volr(i)
         src_ir(i) =bm(i)-src_ir(i)*volr(i)
      ENDDO
!
!.....Build directly solverCSR  array here
!
      CALL csr_build_a(diag_ir,off_diag_ir_non_i,off_diag_ir_non_k)
!
      CALL csr_cg_solvers_scalar(diag_ir,src_ir,rad_ir,eps_imp_rad,max_iter_rad)
!
!.....Build summation info for non, fsw, ctw
!
      nf_number_nb=2
      nf_number_id(0)=0
      nf_number_id(1)=5
      nf_number_id(2)=6
      istart_nfs(0)=0
      istart_nfs(1)=istart_nfs(0)+nf_non
      istart_nfs(2)=istart_nfs(1)+nf_fsw
      lens         =istart_nfs(2)+nf_ctw      
!
!.....Computing cell ; radiation heat flux
!
      nv=0
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         ii=left_nf(i1)
         kk=right_non(i)
         rgamm=f1(i)*rad_gamm(ii)+f0(i)*rad_gamm(kk)
         ir1=rad_ir(ii)
         ir2=rad_ir(kk)
         qrad_nf(i0)=0.0d0
         rad_diff_nf(i0)=rgamm*(ir2-ir1)*sap_nf(i1)
      ENDDO
!
!.....Fluid-solid interface ; radiation heat flux
!   
      nv=1
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         k=right_fsw(i)
         tw=solid%tsol_o(k)      
         rad_irfsw(i)=rad_gamm(ii)*wallnr(ii)*rad_ir(ii)
         rad_irfsw(i)=rad_irfsw(i)+Ew*4.0d0*SIGMA_SBC*tw**4
         rad_irfsw(i)=rad_irfsw(i)/(rad_gamm(ii)*wallnr(ii)+Ew)
         qrad_nf(i0)=Ew*(4.0d0*SIGMA_SBC*tw**4-rad_irfsw(i))
         rad_diff_nf(i0)=qrad_nf(i0)*saa_nf(i1)
         rad_ir_nbcon(i2)=rad_irfsw(i)
      ENDDO
!
!.....Constant temperature wall ; radiation heat flux
!
      nv=2
      nf_number=nf_number_id(nv)
      istart0=istart_nfs(nv)
      istart=istart_nf(1,nf_number)
      istart2=istart_nbcon_nf(nf_number)
      isize =istart_nf(2,nf_number)
      IF(isize.gt.0) CALL udfn_tw_profile(twall_ctw)
      DO i=1,isize
         i0=istart0+i
         i1=istart+i
         i2=istart2+i
         ii=left_nf(i1)
         rad_irctw(i)=rad_gamm(ii)*wallnr(ii)*rad_ir(ii)
         rad_irctw(i)=rad_irctw(i)+Ew*4.0d0*SIGMA_SBC*twall_ctw(i)**4
         rad_irctw(i)=rad_irctw(i)/(rad_gamm(ii)*wallnr(ii)+Ew)
         qrad_nf(i0)=Ew*(4.0d0*SIGMA_SBC*twall_ctw(i)**4-rad_irctw(i))
         rad_diff_nf(i0)=qrad_nf(i0)*saa_nf(i1)
         rad_ir_nbcon(i2)=rad_irctw(i)
      ENDDO
!
      CALL sum_nf(0,-1,         &
                  qrad_nf,qrad)
      CALL sum_nf(0,-1,                   &
                  rad_diff_nf,rad_source)
!
   END SUBROUTINE radiation_model
