!
      SUBROUTINE csr_cg_solvers_scalar(diag,src,sol,eps,maxiter)
!
!     CG solver based on CSR format for scalars, such as eg_t, el_t
!     Used only for SMAC=3 (or implicit calculation)
!
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid
      USE Zcore         , ONLY: np
      USE Zconst1       , ONLY: parallel
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: maxiter
      REAL(8) :: eps
      REAL(8) :: diag(ncell_fluid),src(ncell_fluid)
!.....Output
      REAL(8) :: sol(ncell_fp)
!.....Local variables
      INTEGER :: i
!
!.....Check zero right hand side for serial run
!
      IF(np.eq.1) THEN
         DO i=1,ncell_fluid
            IF(src(i).ne.0.d0) goto 100
         ENDDO
         DO i=1,ncell_fp
            sol(i)=0.d0
         ENDDO
         RETURN
100      CONTINUE
      ENDIF
!
      IF(parallel.ge.1) THEN
         CALL cupid_solvers_scalar(diag,src,sol,eps,maxiter)
      ELSEIF(parallel.eq.2) THEN
!        CALL mg_solve(src,sol)
      ENDIF
!
      END SUBROUTINE csr_cg_solvers_scalar
