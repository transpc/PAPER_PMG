!
      SUBROUTINE cupid_solvers(epsFactor,diag,non_orth,src,sol,izone)
!
!.....This routine select solver based on somaflow input 
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_fp,maxmt,maxmt_pad,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                                au,ia_a,ja_a,ju_a,iend,                                      &
                                ia_r,ja_r,ju_r,                                              &
                                diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                           &
                                ap,iap,jap,jaa,jaar,iaa,ngroup,nbgroup,                      &
                                levt,lev_typet,perm_r,permi_r,index_r
      USE Zzone         , ONLY: ncell_fluid,ncell_fluid_padv
      USE Zcore         , ONLY: np,myrank
      USE Zbicg         , ONLY: psolve,pbcgind
!
      IMPLICIT NONE
!
!.....Input
      INTEGER :: izone
      INTEGER :: non_orth
      REAL(8) :: epsFactor
      REAL(8),DIMENSION(ncell_fluid) :: diag,src
!.....Output
      REAL(8),DIMENSION(ncell_fp) :: sol
!.....Local variables
      LOGICAL :: isPSolve=.true.
!
      IF(levt.eq.0) THEN
         IF(MOD(psolve,2).eq.0)THEN
            IF(lev_typet.eq.0) THEN
!
!..............CG solver ilu0
!
               CALL ilupc11(ncell_fluid,maxmt,maxmt_lu0,maxmt_lu1,   &
                            diag,au,ia_a,ja_a,ju_a,iend,             &
                            diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
            ELSEIF(lev_typet.gt.0) THEN
!
!..............CG solver ilu0 with METIS,CUTHILL reordering
!
               CALL ilupc11p(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                             diag,au,ia_a,ja_a,iend,                        &
                             ia_r,ja_r,ju_r,                                &
                             diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,             &
                             perm_r,permi_r)
            ENDIF
         ENDIF
      ELSEIF(levt.gt.0) THEN
         IF(MOD(psolve,2).eq.0)THEN
            IF(np.eq.1) THEN
               IF(lev_typet.eq.0) THEN
!
!.................CG solver ilup
!
                  CALL factor(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                              diag,au,ia_a,ja_a,iend,                        &
                              ia_r,ja_r,ju_r,                                &
                              diag_lu,alu0,alu1,ia0,ia1,ja0,ja1)
               ELSE
!
!.................CG solver ilup with METIS,CUTHILL reordering
!
                  CALL factorp(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                               diag,au,ia_a,ja_a,iend,                               &
                               ia_r,ja_r,ju_r,                                       &
                               diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                    &
                               perm_r,permi_r,index_r)
               ENDIF
            ELSE
               WRITE(*,*) 'cupid_solvers: ilup solvet not yet implemented for MPI'
               CALL finalize_mpi
               STOP  
            ENDIF
         ENDIF
      ELSEIF(levt.eq.-1) THEN
         IF(np.eq.1) THEN
            IF(lev_typet.eq.0) THEN
               IF(non_orth.gt.0)THEN
!
!.................Direct solver no reordering save LU
!
                  CALL factor_solve(ncell_fluid,maxmt,maxmt_r,maxmt_lu0,maxmt_lu1, &
                                    diag,au,ia_a,ja_a,iend,                        &
                                    ia_r,ja_r,ju_r,                                &
                                    diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,             &
                                    ncell_fp,src,sol)
               ELSE
!
!.................Direct solver no reordering
!
                  CALL factor_solve0(ncell_fluid,maxmt,maxmt_r,maxmt_lu1, &
                                     diag,au,ia_a,ja_a,iend,              &
                                     ia_r,ja_r,ju_r,                      &
                                     ncell_fp,src,sol)

               ENDIF
            ELSEIF(lev_typet.gt.0) THEN
               IF(non_orth.gt.0)THEN
!
!.................Direct solver with METIS,CUTHILL ordering save LU
!
                  CALL factor_solvep(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu0,maxmt_lu1, &
                                     diag,au,ia_a,ja_a,iend,                               &
                                     ia_r,ja_r,ju_r,                                       &
                                     diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                    &
                                     ncell_fp,src,sol,perm_r,permi_r,index_r)
               ELSE
!
!.................Direct solver with METIS,CUTHILL ordering save LU
!
                  CALL factor_solve0p(ncell_fluid,maxmt,maxmt2,maxmt_r,maxmt_lu1, &
                                      diag,au,ia_a,ja_a,iend,                     &
                                      ia_r,ja_r,ju_r,                             &
                                      ncell_fp,src,sol,perm_r,permi_r,index_r)
                  ENDIF
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers: Direct solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP  
         ENDIF
      ELSEIF(levt.eq.-2) THEN
         IF(np.eq.1) THEN
            IF(non_orth.gt.0)THEN
!
!..............Tridiagonal solver save LU
!
               CALL factor_solve_tri(ncell_fluid,maxmt,        &
                                     ia_a,au,              &
                                     diag_lu,alu0,alu1,  &
                                     ncell_fp,src,sol)
            ELSE
!
!..............Tridiagonal solver
!
               CALL factor_solve_tri0(ncell_fluid,maxmt,  &
                                      ia_a,au,        &
                                      ncell_fp,src,sol)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers: Tridiagonal solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP  
         ENDIF
      ENDIF
!
      IF(levt.ge.0) THEN
!
!........Copy a array in block format for vector processing
!
         CALL copy_a_vector(ncell_fluid,maxmt,maxmt_pad,      &
                            ia_a,au,                       &
                            iap,ap,jaa,iaa,ngroup,nbgroup)
!
         IF(izone.eq.1) THEN
            IF(pbcgind.gt.0)THEN
               pbcgind=0
               IF(myrank.eq.0)WRITE(*,*)'Iteration failed already before solid!!!'
            ENDIF
         ENDIF
!
!........Solve the pressure matrix
!
         CALL csr_cg_solver(epsFactor,ncell_fluid,ncell_fp,ncell_fluid_padv,maxmt,maxmt_pad,maxmt_lu0,maxmt_lu1, &
                            diag,au,ia_a,ja_a,ju_a,                                                              &
                            diag_lu,alu0,alu1,ia0,ia1,ja0,ja1,                                                   &
                            ap,iap,jap,jaar,iaa,ngroup,nbgroup,                                                  &
                            lev_typet,perm_r,                                                                    &
                            src,sol,izone,isPSolve)
      ENDIF
!
      END SUBROUTINE cupid_solvers
