      MODULE Zbicg
! 
      IMPLICIT NONE
      SAVE
!
! bicgstab & parallel-bicg
      INTEGER min_bicg,max_bicg, psolve
      INTEGER lev_type,levmpi_type,lev,levmpi
      INTEGER lev_type_c,levmpi_type_c,lev_c,levmpi_c
      INTEGER pbcgind,pbcgsig,pbcgind_max,pbcgsig_max
      REAL(8) eps_bicg,relax_u,relax_p,relax_e,relax_t
!
      END MODULE Zbicg
