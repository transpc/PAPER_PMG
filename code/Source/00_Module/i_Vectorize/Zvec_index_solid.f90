      MODULE Zvec_index_solid
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: left_solid_k,right_solid_k,                        &
                                          left_solid_nf,jneigh_solid_nf,                     &
                                          right_solid_non,right_solid_fsw,                   &
                                          nbcon_solid_ctw,nbcon_solid_chw,nbcon_solid_chtcw, &
                                          f_fluid_fsw
      REAL(8),DIMENSION(:),ALLOCATABLE:: flux_fsw,tliq_fsw,tsol_fsw,qliq_fsw,qgas_fsw, &
                                         dfilm_fsw,vfilm_fsw,dfilm_ctw,vfilm_ctw
!
      END MODULE Zvec_index_solid      
