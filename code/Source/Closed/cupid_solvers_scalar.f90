!
      SUBROUTINE cupid_solvers_scalar(diag,src,sol,eps,maxiter)
!
!     CG solver based on CSR format for scalars, such as eg_t, el_t
!     Used only for SMAC=3 (or implicit calculation)
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1,  &
                                au,ia_a,ja_a,ju_a,iend,                        &
                                diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,             &
                                ap,jap,iap,jaa,jaar,iaa,ngroup,nbgroup,        &
                                maxmt2,maxmt_r,ia_r,ja_r,ju_r,                 &
                                perm_r,permi_r,index_r,                        &
                                levt,lev_typet
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_padv
      USE Zcore         , ONLY: np
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
      INTEGER :: izone=0
!
      IF(levt.eq.0) THEN
         IF(lev_typet.eq.0) THEN
!
!...........CG solver ilu0 with METIS,CUTHILL reordering
!
            CALL ilupc11(ncell_fluid,maxmt,maxmt_lu0,maxmt_lu1, &
                         diag,au,ia_a,ja_a,ju_a,iend,           &
                         diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
         ELSE
!
!...........CG solver ilu0 with METIS,CUTHILL reordering
!
            CALL ilupc11p(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                          diag,au,ia_a,ja_a,iend,                        &
                          ia_r,ja_r,ju_r,                                &
                          diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,             &
                          perm_r,permi_r)
         ENDIF
      ELSEIF(levt.gt.0) THEN
         IF(np.eq.1) THEN
            IF(lev_typet.eq.0) THEN
!
!..............CG solver ilup
!
               CALL factor(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                           diag,au,ia_a,ja_a,iend,                        &
                           ia_r,ja_r,ju_r,                                &
                           diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
            ELSE
!
!..............CG solver ilup with METIS,CUTHILL reordering
!
               CALL factorp(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                            diag,au,ia_a,ja_a,iend,                               &
                            ia_r,ja_r,ju_r,                                       &
                            diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                    &
                            perm_r,permi_r,index_r)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_scalar: ilup solver not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ELSEIF(levt.eq.-1) THEN
         IF(np.eq.1) THEN
!
!...........Direct solver no reordering
!
            IF(lev_typet.eq.0) THEN
               CALL factor_solve0(ncell_fluid,maxmt,maxmt_r,maxmt_lu1, &
                                  diag,au,ia_a,ja_a,iend,              &
                                  ia_r,ja_r,ju_r,                      &
                                  ncell_fp,src,sol)
            ELSE
!
!..............Direct solver with METIS,CUTHILL ordering save LU
!
               CALL factor_solve0p(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                   diag,au,ia_a,ja_a,iend,                     &
                                   ia_r,ja_r,ju_r,                             &
                                  ncell_fp,src,sol,perm_r,permi_r,index_r)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_scalar: Direct solver not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ELSEIF(levt.eq.-2) THEN
         IF(np.eq.1) THEN
            CALL factor_solve_tri0(ncell_fluid,maxmt, &
                                   ia_a,au,             &
                                   ncell_fp,src,sol)
         ELSE
            WRITE(*,*) 'cupid_solvers_scalar: Tridiagonal solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ENDIF
!
      IF(levt.ge.0) THEN
!
!........Copy a array in block format for vector processing
!
         CALL copy_a_vector(ncell_fluid,maxmt,maxmt_pad,   &
                            ia_a,au,                       &
                            iap,ap,jaa,iaa,ngroup,nbgroup)
!
!........Apply ILU pre-conditioned bi-conjugate gradient solver
!
         CALL pbcg_ilu(eps,maxiter,ncell_fluid,ncell_fluid_padv,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                       diag,                                                                   &
                       diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                      &
                       iap,jap,ap,jaar,iaa,ngroup,nbgroup,                                     &
                       lev_typet,perm_r,                                                       &
                       ncell_fp,src,sol,izone)
      ENDIF
!
      END SUBROUTINE cupid_solvers_scalar
