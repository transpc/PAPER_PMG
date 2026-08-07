!
      SUBROUTINE implicit_momentum(iter)
!
!.....This routine calculates implicit momentum convection and diffusion
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zcore         , ONLY: np
      USE Zparam        , ONLY: ndim
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Zare          , ONLY: ar_gas,ar_liq
      USE Zconst2       , ONLY: dt
      USE Zcoord3       , ONLY: volp
      USE Zimplicit     , ONLY: imp_mom_conv,imp_mom_diff,ag_min_m,al_min_m,iter_mom
      USE Zvector       , ONLY: vg_t,vl_t,vd_t,vg_n,vl_n,vd_n
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: iter
!.....Local variables
      INTEGER :: i
      REAL(8) :: vt
!.....Local arrays
      REAL(8),DIMENSION(ncell_fluid) :: diag_g,diag_l
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src_g,src_l
!.....Local vector arrays
      REAL(8),DIMENSION(nf_non) :: off_diag_g_non_i,off_diag_l_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_g_non_k,off_diag_l_non_k
!
      IF(iter.gt.1.and.iter.eq.iter_mom) THEN
         IF(ndim.eq.2) THEN
            DO i=1,ncell_fp
               vg_n(i,1)=vg_t(i,1)
               vg_n(i,2)=vg_t(i,2)
               vl_n(i,1)=vl_t(i,1)
               vl_n(i,2)=vl_t(i,2)
               vd_n(i,1)=vd_t(i,1)
               vd_n(i,2)=vd_t(i,2)
            ENDDO
         ELSE
            DO i=1,ncell_fp
               vg_n(i,1)=vg_t(i,1)
               vg_n(i,2)=vg_t(i,2)
               vg_n(i,3)=vg_t(i,3)
               vl_n(i,1)=vl_t(i,1)
               vl_n(i,2)=vl_t(i,2)
               vl_n(i,3)=vl_t(i,3)
               vd_n(i,1)=vd_t(i,1)
               vd_n(i,2)=vd_t(i,2)
               vd_n(i,3)=vd_t(i,3)
            ENDDO
         ENDIF
         RETURN
      ENDIF
!
!.....Initial values for the diagonal and source terms
!
      DO i=1,ncell_fluid
         vt=volp(i)/dt
         diag_g(i)=ar_gas(i)*vt
         diag_l(i)=ar_liq(i)*vt
      ENDDO
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fluid
            src_g(i,1)=diag_g(i)*vg_t(i,1)
            src_g(i,2)=diag_g(i)*vg_t(i,2)
            src_l(i,1)=diag_l(i)*vl_t(i,1)
            src_l(i,2)=diag_l(i)*vl_t(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fluid
            src_g(i,1)=diag_g(i)*vg_t(i,1)
            src_g(i,2)=diag_g(i)*vg_t(i,2)
            src_g(i,3)=diag_g(i)*vg_t(i,3)
            src_l(i,1)=diag_l(i)*vl_t(i,1)
            src_l(i,2)=diag_l(i)*vl_t(i,2)
            src_l(i,3)=diag_l(i)*vl_t(i,3)
         ENDDO
      ENDIF
!
      IF(iter.gt.1.and.np.gt.1) CALL communicate_2d(vg_n, &
                                                    vl_n)
!
!.....Write A-Matrix and b-vector (Ax=b) for momentum diffusion and momentum convection
!
      IF(imp_mom_diff.eq.1) THEN
         CALL imp_mom_diffusion(diag_g,diag_l,src_g,src_l,iter,    &
                                off_diag_g_non_i,off_diag_g_non_k, &
                                off_diag_l_non_i,off_diag_l_non_k)
      ELSE
!
!........When imp_mom_diffusion not called, imp_mom_convection does add
!
         DO i=1,nf_non
            off_diag_g_non_i(i)=0.d0
            off_diag_l_non_i(i)=0.d0
         ENDDO
         DO i=1,nf_nonk
            off_diag_g_non_k(i)=0.d0
            off_diag_l_non_k(i)=0.d0
         ENDDO
      ENDIF
      IF(imp_mom_conv.eq.1) CALL imp_mom_convection(diag_g,src_g,diag_l,src_l,iter,    &
                                                    off_diag_g_non_i,off_diag_l_non_i, &
                                                    off_diag_g_non_k,off_diag_l_non_k)
!
!.....Limit value control
!
      DO i=1,ncell_fluid
         IF(cell%alphag_o(i).lt.ag_min_m) diag_g(i)=1.d0
         IF(cell%alphal_o(i).lt.al_min_m) diag_l(i)=1.d0
      ENDDO
!
!.....Apply a CG solver to get the solution (x=A-1b)
!
!.....build directly solver CSR array here
!
      CALL csr_build_a(diag_g,off_diag_g_non_i,off_diag_g_non_k)
!
      CALL csr_cg_solvers_vector(diag_g,src_g,vg_n)
!
!.....build directly solver CSR array here
!
      CALL csr_build_a(diag_l,off_diag_l_non_i,off_diag_l_non_k)
!
      CALL csr_cg_solvers_vector(diag_l,src_l,vl_n)
!
      IF(ndim.eq.2) THEN
         DO i=1,ncell_fp
            vd_n(i,1)=vd_t(i,1)
            vd_n(i,2)=vd_t(i,2)
         ENDDO
      ELSE
         DO i=1,ncell_fp
            vd_n(i,1)=vd_t(i,1)
            vd_n(i,2)=vd_t(i,2)
            vd_n(i,3)=vd_t(i,3)
         ENDDO
      ENDIF
!
      END SUBROUTINE implicit_momentum
