!
      SUBROUTINE cupid_solvers_c(epsFactor,diag,src,sol)
!
!.....This routine select solver based on somaflow input 
!
      USE Zinterface
      USE Zmpi          , ONLY: ncell_ps,maxmt_c,maxmt_pad_c,maxmt_lu0_c,maxmt_lu1_c,   &
                                maxmt2_c,maxmt_r_c,                                     &
                                au_c,ju_a_c,ja_a_c,ia_a_c,iend_c,                       &
                                ia_r_c,ja_r_c,ju_r_c,                                   &
                                diag_lu_c,alu0_c,alu1_c,ja0_c,ja1_c,ia0_c,ia1_c,        &
                                ap_c,iap_c,jap_c,jaa_c,jaar_c,iaa_c,ngroup_c,nbgroup_c, &
                                levt_c,lev_typet_c,perm_r_c,permi_r_c,index_r_c
      USE Zzone         , ONLY: ncell_cond,ncell_cond_padv
      USE Zcore         , ONLY: np,myrank
      USE Zbicg         , ONLY: psolve,pbcgind
!
      IMPLICIT NONE
!
!.....Input
      REAL(8) :: epsFactor
      REAL(8),DIMENSION(ncell_cond) :: diag,src
!.....Output
      REAL(8),DIMENSION(ncell_ps) :: sol
!.....Local variables
      INTEGER :: izone=1
      LOGICAL :: isPSolve=.true.
!
      IF(levt_c.eq.0) THEN
         IF(MOD(psolve,2).eq.0)THEN
            IF(lev_typet_c.eq.0) THEN
!
!.................CG solver ilu0
!
               CALL ilupc11(ncell_cond,maxmt_c,maxmt_lu0_c,maxmt_lu1_c,      &
                            diag,au_c,ia_a_c,ja_a_c,ju_a_c,iend_c,           &
                            diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c)
            ELSEIF(lev_typet_c.gt.0) THEN
!
!.................CG solver ilu0 with METIS,CUTHILL reordering
!
                  CALL ilupc11p(ncell_cond,maxmt_c,maxmt_r_c,maxmt_lu0_c,maxmt_lu1_c, &
                                diag,au_c,ia_a_c,ja_a_c,iend_c,                       &
                                ia_r_c,ja_r_c,ju_r_c,                                 &
                                diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c,      &
                                perm_r_c,permi_r_c)
            ENDIF
         ENDIF
      ELSEIF(levt_c.gt.0) THEN
         IF(MOD(psolve,2).eq.0)THEN
            IF(np.eq.1) THEN
               IF(lev_typet_c.eq.0) THEN

!.................CG solver ilup
!
                  CALL factor(ncell_cond,maxmt_c,maxmt_r_c,maxmt_lu0_c,maxmt_lu1_c, &
                              diag,au_c,ia_a_c,ja_a_c,iend_c,                       &
                              ia_r_c,ja_r_c,ju_r_c,                                 &
                              diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c)
               ELSE
!
!.................CG solver ilup with METIS,CUTHILL reordering
!
                  CALL factorp(ncell_cond,maxmt_c,maxmt2_c,maxmt_r_c,maxmt_lu0_c,maxmt_lu1_c, &
                               diag,au_c,ia_a_c,ja_a_c,iend_c,                                &
                               ia_r_c,ja_r_c,ju_r_c,                                          &
                               diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c,               &
                               perm_r_c,permi_r_c,index_r_c)
               ENDIF
            ELSE
               WRITE(*,*) 'cupid_solvers_c: ilup solvet not yet implemented for MPI'
               CALL finalize_mpi
               STOP
            ENDIF
         ENDIF
      ELSEIF(levt_c.eq.-1) THEN
         IF(np.eq.1) THEN
            IF(lev_typet_c.eq.0) THEN
!
!.................Direct solver no reordering
!
                  CALL factor_solve0(ncell_cond,maxmt_c,maxmt_r_c,maxmt_lu1_c, &
                                     diag,au_c,ia_a_c,ja_a_c,iend_c,           &
                                     ia_r_c,ja_r_c,ju_r_c,                     &
                                     ncell_ps,src,sol)

            ELSEIF(lev_typet_c.gt.0) THEN
!
!.................Direct solver with METIS,CUTHILL ordering 
!
                  CALL factor_solve0p(ncell_cond,maxmt_c,maxmt2_c,maxmt_r_c,maxmt_lu1_c, &
                                      diag,au_c,ia_a_c,ja_a_c,iend_c,                    &
                                      ia_r_c,ja_r_c,ju_r_c,                              &
                                      ncell_ps,src,sol,perm_r_c,permi_r_c,index_r_c)
            ENDIF
         ELSE
            WRITE(*,*) 'cupid_solvers_c: Direct solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ELSEIF(levt_c.eq.-2) THEN
         IF(np.eq.1) THEN
!
!...........Tridiagonal solver
!
            CALL factor_solve_tri0(ncell_cond,maxmt_c, &
                                   ia_a_c,au_c,        &
                                   ncell_ps,src,sol)
         ELSE
            WRITE(*,*) 'cupid_solvers_c: Tridiagonal solvet not yet implemented for MPI'
            CALL finalize_mpi
            STOP
         ENDIF
      ENDIF
!
      IF(levt_c.ge.0) THEN
!
!........Copy a array in block format for vector processing
!
         CALL copy_a_vector(ncell_cond,maxmt_c,maxmt_pad_c,            &
                            ia_a_c,au_c,                               &
                            iap_c,ap_c,jaa_c,iaa_c,ngroup_c,nbgroup_c)
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
         CALL csr_cg_solver(epsFactor,ncell_cond,ncell_ps,ncell_cond_padv,maxmt_c,maxmt_pad_c,maxmt_lu0_c,maxmt_lu1_c, &
                            diag,au_c,ia_a_c,ja_a_c,ju_a_c,                                                            &
                            diag_lu_c,alu0_c,alu1_c,ia0_c,ia1_c,ja0_c,ja1_c,                                           &
                            ap_c,iap_c,jap_c,jaar_c,iaa_c,ngroup_c,nbgroup_c,                                          &
                            lev_typet_c,perm_r_c,                                                                      &
                            src,sol,izone,isPSolve)
      ENDIF
!
      END SUBROUTINE cupid_solvers_c
