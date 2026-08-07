!
!!!!!!!
         ALLOCATE(face%twall_partition(ncell_fp))
         ALLOCATE(face%wall_fluxl_diff(ncell_fp),face%wall_fluxg_diff(ncell_fp),face%wall_fluxd_diff(ncell_fp))
         ALLOCATE(face%ddepartw(ncell_fp),face%ratio_evap(ncell_fp),face%bfreq(ncell_fp))