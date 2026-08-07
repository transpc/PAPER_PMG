       IF(ncell_cond.gt.0) THEN
         ALLOCATE(solid%tsol(ncell_cond))
         ALLOCATE(solid%rhos(ncell_cond))
         ALLOCATE(solid%conds(ncell_cond),solid%cps(ncell_cond),solid%rhocps(ncell_cond))
         ALLOCATE(solid%tsol_o(ncell_cond))
         ALLOCATE(solid%tsol_max(ncell_cond),solid%tpellet_surf(ncell_cond))
         ALLOCATE(solid%temp_rod(ncell_cond,10),solid%hconv_rod_g(ncell_cond))
         ALLOCATE(solid%hconv_rod_l(ncell_cond))
         ALLOCATE(solid%matnum(ncell_cond),solid%idummy(ncell_cond))
       ENDIF
