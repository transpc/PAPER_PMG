!
      SUBROUTINE csr_cg_solvers_vector(diag,src,sol)
!  
!     CG solver based on CSR format for vectors, such as vl_n, vg_n
!     Used only for SMAC=3 (or implicit calculation)
!
      USE Zmpi          , ONLY: ncell_fp
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_pad
      USE Zcore         , ONLY: np
      USE Zparam        , ONLY: ndim
      USE Zconst1       , ONLY: parallel
!
      IMPLICIT NONE
!
!.....Input
      REAL(8),DIMENSION(ncell_fluid) :: diag
      REAL(8),DIMENSION(ncell_fluid_pad,ndim) :: src
!.....Output
      REAL(8) :: sol(ncell_fp,ndim)
!.....Local variables
      INTEGER :: i,ix
      INTEGER :: iswitch_factor(ndim),iswitch_factor_global
!
!.....Check zero right hand side for serial run
!
      IF(np.eq.1) THEN
         iswitch_factor_global=0
         DO ix=1,ndim
            iswitch_factor(ix)=0
            DO i=1,ncell_fluid
               IF(src(i,ix).ne.0.d0) goto 100
            ENDDO
            iswitch_factor_global=iswitch_factor_global+1
            iswitch_factor(ix)=1
            DO i=1,ncell_fp
               sol(i,ix)=0.d0
            ENDDO
100         CONTINUE
         ENDDO
         IF(iswitch_factor_global.eq.ndim) RETURN
      ENDIF
!
      IF(parallel.ge.1) THEN
         CALL cupid_solvers_vector(diag,src,sol,iswitch_factor)
      ELSEIF(parallel.eq.2) THEN
         DO ix=1,ndim
!           CALL mg_solve(src(1,ix),sol(1,ix))
         ENDDO
      ENDIF
!
      END SUBROUTINE csr_cg_solvers_vector
