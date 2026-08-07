!
      SUBROUTINE implicit_scalar(eg_1,el_1,xn_1,eg_2,el_2,xn_2,a,iter)
!
!.....This routine calculates implicit scalar convection and diffusion
!
      USE VOL_DATA
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zvec_param    , ONLY: nf_nonk,nf_non
      USE Zconst2       , ONLY: dtr
      USE Zcoord3       , ONLY: volp
      USE Zimplicit     , ONLY: imp_scalar_diff,imp_scalar_conv,eps_imp_scalar,max_iter_scalar
!
      IMPLICIT NONE
!
!.....Input
      INTEGER iter
      REAL(8) vt
      REAL(8) a(ncell_fluid,3)
      REAL(8) eg_2(ncell_fluid),el_2(ncell_fluid),xn_2(ncell_fluid)      
!.....Output
      REAL(8) eg_1(ncell_fluid),el_1(ncell_fluid),xn_1(ncell_fluid)
!.....Local variables
      INTEGER i
!.....Local array
      REAL(8) diag_g(ncell_fluid),src_g(ncell_fluid)
      REAL(8) diag_l(ncell_fluid),src_l(ncell_fluid)
      REAL(8) diag_x(ncell_fluid),src_x(ncell_fluid)
      REAL(8) solg(ncell_fp),soll(ncell_fp),solx(ncell_fp)
!.....Local vector array
      REAL(8) :: off_diag_g_non_i(nf_non),off_diag_l_non_i(nf_non),off_diag_x_non_i(nf_non)
      REAL(8) :: off_diag_g_non_k(nf_nonk),off_diag_l_non_k(nf_nonk),off_diag_x_non_k(nf_nonk)
!
!.....Initial values for the diagonal and source terms
!      
      DO i=1,ncell_fluid
         vt=volp(i)*dtr 
         diag_g(i)=a(i,1)*vt
         diag_l(i)=a(i,2)*vt
         diag_x(i)=a(i,3)*vt
!
         IF(iter.eq.1)THEN
            src_g(i)=0.0d0
            src_l(i)=0.0d0
            src_x(i)=0.0d0
         ELSE
            src_g(i)=diag_g(i)*(eg_2(i)-eg_1(i))
            src_l(i)=diag_l(i)*(el_2(i)-el_1(i))
            src_x(i)=diag_x(i)*(xn_2(i)-xn_1(i))
         ENDIF
!
      ENDDO      
!
!.....Define A-matrix and b-vector for implicit gas energy diffusion (Ax=b)
!
!.....Define A-matrix and b-vector for implicit liquid energy diffusion (Ax=b)
!
      IF(imp_scalar_conv.eq.1)THEN
         CALL imp_scalar_convection(diag_g,diag_l,diag_x,src_g,src_l,src_x,             &
                                    off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i, &
                                    off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k, &
                                    iter)
      ELSE
        DO i=1,nf_non
           off_diag_g_non_i(i)=0.0d0
           off_diag_l_non_i(i)=0.0d0
           off_diag_x_non_i(i)=0.0d0
        ENDDO
        DO i=1,nf_nonk
           off_diag_g_non_k(i)=0.0d0
           off_diag_l_non_k(i)=0.0d0
           off_diag_x_non_k(i)=0.0d0
        ENDDO
      ENDIF  
! 
      IF(imp_scalar_diff.eq.1)THEN
         CALL imp_eng_diffusion(diag_g,diag_l,diag_x,src_g,src_l,src_x,             &
                                off_diag_g_non_i,off_diag_l_non_i,off_diag_x_non_i, &
                                off_diag_g_non_k,off_diag_l_non_k,off_diag_x_non_k, &
                                iter)
      ENDIF    
!
!.....Apply a CG solver to get the solutions
!
      CALL csr_build_a(diag_g,off_diag_g_non_i,off_diag_g_non_k)
      CALL csr_cg_solvers_scalar(diag_g,src_g,solg,eps_imp_scalar,max_iter_scalar)
      CALL csr_build_a(diag_l,off_diag_l_non_i,off_diag_l_non_k)
      CALL csr_cg_solvers_scalar(diag_l,src_l,soll,eps_imp_scalar,max_iter_scalar)
      CALL csr_build_a(diag_x,off_diag_x_non_i,off_diag_x_non_k)
      CALL csr_cg_solvers_scalar(diag_x,src_x,solx,eps_imp_scalar,max_iter_scalar)
!
      DO i=1,ncell_fluid
         eg_1(i)=eg_1(i)+solg(i)
         el_1(i)=el_1(i)+soll(i) 
         xn_1(i)=xn_1(i)+solx(i)
      ENDDO
!
      RETURN
      END SUBROUTINE implicit_scalar
!
