      SUBROUTINE ncg_transport_imp
!
!     This routine calculates ncg transport implicitly
!
      USE Zinterface
      USE VOL_DATA      , ONLY: cell
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Zare          , ONLY: ar_gas
      USE Zconst2       , ONLY: dtr
      USE Zcoord3       , ONLY: volp
      USE Zimplicit     , ONLY: eps_imp_scalar,max_iter_scalar
      USE Zncg          , ONLY: n_ncg_sp,qn_cell,qn_cell_o,ncg_diff
      USE Ztimecon      , ONLY: alpha_min
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER i,nc
      REAL(8) vt
!.....Local array
      REAL(8),DIMENSION(ncell_fluid) ::  diag,src
      REAL(8),DIMENSION(ncell_fp) ::  sol,qn
!.....Local vector array
      REAL(8),DIMENSION(nf_non) :: off_diag_non_i
      REAL(8),DIMENSION(nf_nonk) :: off_diag_non_k
!      
      IF(ncg_diff.gt.0) THEN
         IF(np.gt.1) CALL communicate_1d(ar_gas,     &
                                         cell%quala, &
                                         cell%mdiff)
      ELSE
         IF(np.gt.1) CALL communicate_1d(ar_gas,     &
                                         cell%quala)
      ENDIF
!
      DO nc=1,n_ncg_sp
!
         DO i=1,ncell_fluid
            qn(i)=qn_cell_o(i,nc)
         ENDDO
!
         IF(np.gt.1) CALL communicate_1d(qn)
!
         DO i=1,ncell_fluid
            vt=volp(i)*dtr
            diag(i)=ar_gas(i)*cell%quala(i)*vt
            src(i)=diag(i)*qn(i)
         ENDDO
!
!
!.....Define A-matrix and b-vector for implicit ncg convection (Ax=b)
!
         CALL ncg_convection_imp(diag,src,qn,nc, &
                                 off_diag_non_i,off_diag_non_k)
!
!.....Limit value control
!
         DO i=1,ncell_fluid
            IF(    cell%alphag(i)  .le.alpha_min &
               .or.cell%alphag_o(i).le.alpha_min &
               .or.cell%quala(i)   .le.alpha_min) THEN
               diag(i)=1.0d0
               src(i)=qn(i)
            ENDIF
         ENDDO
!
!........Apply a CG solver to get the solutions
!
!........Build directly solverCSR  array here
!
         CALL csr_build_a(diag,off_diag_non_i,off_diag_non_k)
!
         CALL csr_cg_solvers_scalar(diag,src,sol,eps_imp_scalar,max_iter_scalar)
!
         DO i=1,ncell_fluid
            qn_cell(i,nc)=MIN(MAX(sol(i),0.0d0),1.0d0)
         ENDDO
!
      ENDDO
!
      END SUBROUTINE ncg_transport_imp
