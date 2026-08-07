!
      SUBROUTINE restart_write_rv(nout)
!
!     This routine writes the restart file for RV model
!
      USE Zconst1      , ONLY: lrestart_overwrite
      USE Zconst2      , ONLY: dt
      USE Zcore        , ONLY: myrank
      USE Ztimecon     , ONLY: time,itim      
      USE Zrv_hts_1d   , ONLY: ncell_hts_1d,nr_1d,t_hts_1d
      USE Zrv_hts_2d   , ONLY: nrod_2d,nz0_2d,nr_2d,t_fuel,     &
                                wet_b,wet_t
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nz_fine,p3d_cupid                 
      use unitManager  , only: createUnit      
      USE Zio_unit     , ONLY: unit_restart_rv
!
      IMPLICIT NONE
!
      INTEGER nz0,nz2d,nz1d
!
      INTEGER i,j,k,nout
!      
      CHARACTER*16 myrank_chr
      CHARACTER(50) f_restart      
!      
      LOGICAL, SAVE::INITIAL_r
!
      DATA INITIAL_r /.TRUE./
!     return
      nz1d=ncell_hts_1d
      nz0=nz0_2d
      nz2d=nz0_2d*nz_fine
!
!.....Write restart file (open separate files by myrank number )
!
      IF(INITIAL_r)THEN
         WRITE(myrank_chr,*)myrank
         myrank_chr=adjustl(myrank_chr)
         f_restart='restart_rv_'//trim(myrank_chr)//'.dat'
         !OPEN(800+myrank,file=f_restart,form='unformatted')
         unit_restart_rv=createUnit("restart_rv_write")
         unit_restart_rv=800 !+myrank
         OPEN(unit_restart_rv,file=f_restart)
         INITIAL_r=.FALSE.                  
      ELSE
         IF(lrestart_overwrite)THEN
            CLOSE(unit_restart_rv)
            WRITE(myrank_chr,*)myrank
            myrank_chr=adjustl(myrank_chr)
            f_restart='restart_rv_'//trim(myrank_chr)//'.dat'
            !OPEN(800+myrank,file=f_restart,form='unformatted')
            unit_restart_rv=800 !+myrank
            OPEN(unit_restart_rv,file=f_restart)
            INITIAL_r=.FALSE.      
         ENDIF                    
      ENDIF
!
      WRITE(unit_restart_rv,5001) nout,itim,time,dt
!
      IF(nrod_2d.gt.0)THEN
      DO k=1,nrod_2d
         WRITE(unit_restart_rv,5002) wet_b(k),wet_t(k)
      ENDDO  
      ENDIF
!
      IF(ncell_fuel_rod.gt.0.and.nr_2d.gt.0)THEN
      DO k=1,ncell_fuel_rod
         DO i=1,nr_2d
            WRITE(unit_restart_rv,5002) t_fuel(k,i)
         ENDDO
      ENDDO       
      ENDIF
!
      IF(nz1d.gt.0.and.nr_1d.gt.0)THEN
      DO j=1,nz1d
         DO i=1,nr_1d
            WRITE(unit_restart_rv,5002)t_hts_1d(j,i)
         ENDDO
      ENDDO 
      ENDIF
!      
      IF(ncell_fuel_rod.gt.0.and.nr_2d.gt.0)THEN
          DO k=1,ncell_fuel_rod
              WRITE(unit_restart_rv,5002)p3d_cupid(k)
          ENDDO       
      ENDIF
!         
!      
   5001 FORMAT(2(i10,1x),2e20.10)   
   5002 FORMAT(5e20.10)    
   5003 FORMAT(e20.10,2(i10,1x))  
!
!
      RETURN
      END SUBROUTINE restart_write_rv
