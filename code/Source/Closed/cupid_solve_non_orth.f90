!
      SUBROUTINE cupid_solve_non_orth(epsFactor,diag,src,sol,izone)
!
!     This routine establish pressure matrix which is solved using conjugate gradient solvers
!     For parallel computing, CZR format is applied for pressure matrix.
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                                au,ia_a,ja_a,ju_a,                            &
                                diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,            &
                                ap,jap,iap,jaar,iaa,ngroup,nbgroup,           &
                                levt,lev_typet,                               &
                                perm_r
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_padv
!
      IMPLICIT NONE
!      
!.....Input
      INTEGER :: izone
      REAL(8) :: epsFactor
      REAL(8) :: diag(ncell_fluid)
      REAL(8) :: src(ncell_fluid)
!.....Output
      REAL(8) :: sol(ncell_fp)
!.....Local variables
      LOGICAL :: isPSolve=.true.
!
      IF(levt.eq.-1) THEN
         IF(lev_typet.eq.0) THEN
            CALL lusol00(ncell_fluid,maxmt_lu0,maxmt_lu1,   &
                         diag,                        &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                         ncell_fp,src,sol)
         ELSEIF(lev_typet.gt.0) THEN
            CALL lusol00r(ncell_fluid,maxmt_lu0,maxmt_lu1,   &
                          diag,                        &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1, &
                          ncell_fp,src,sol,perm_r)
         ENDIF
      ELSEIF(levt.eq.-2) THEN
            CALL lusol_tri(ncell_fluid,       &
                           diag_lu,alu0,alu1, &
                           ncell_fp,src,sol)
      ELSE
            CALL csr_cg_solver(epsFactor,ncell_fluid,ncell_fp,ncell_fluid_padv,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                               diag,au,ia_a,ja_a,ju_a,                                                              &
                               diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                                   &
                               ap,iap,jap,jaar,iaa,ngroup,nbgroup,                                                  &
                               lev_typet,perm_r,                                                                    &
                               src,sol,izone,isPSolve) 
      ENDIF
!
      END SUBROUTINE cupid_solve_non_orth
