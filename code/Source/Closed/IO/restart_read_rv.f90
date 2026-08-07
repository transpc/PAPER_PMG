!
      SUBROUTINE restart_read_rv
!
!     This routine reads the restart file for RV model
!
      USE Zconst1      , ONLY: restart,cplmaster
      USE Zconst2      , ONLY: dt
      USE Zio_unit     , ONLY: unit_log,unit_restart_rv
      USE Zcore        , ONLY: myrank
      USE Ztimecon     , ONLY: time,itim,itim_last      
      
      USE Zrv_hts_1d   , ONLY: ncell_hts_1d,nr_1d,t_hts_1d
      USE Zrv_hts_2d   , ONLY: nrod_2d,nz0_2d,nr_2d,t_fuel,     &
                                wet_b,wet_t
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nz_fine,p3d_cupid
      USE Zrv_ncell    , ONLY: cupid_cell_hts2d,qvol_mas      
      use unitManager  , only: createUnit
!
      IMPLICIT NONE
!
      INTEGER nz0,nz2d,nz1d
!
      INTEGER i,j,k,nout,restart_start,restart_step
!      
      CHARACTER*16 myrank_chr
      CHARACTER(50) f_restart      
!      
      LOGICAL, SAVE::INITIAL_r
!
      DATA INITIAL_r /.TRUE./
      nz1d=ncell_hts_1d
      nz0=nz0_2d
      nz2d=nz0_2d*nz_fine 
!
!.....READ restart file (open separate files by myrank number )
!
      IF(INITIAL_r)THEN
         WRITE(myrank_chr,*)myrank
         myrank_chr=adjustl(myrank_chr)
         f_restart='restart_rv_'//trim(myrank_chr)//'.dat'
         unit_restart_rv=createUnit("restart_rv_read")
         !OPEN(800+myrank,file=f_restart,status='old',form='unformatted',iostat=restart_step)
         unit_restart_rv=800 !+myrank
         OPEN(unit_restart_rv,file=f_restart,status='old',iostat=restart_step)
         IF(restart_step.ne.0)then
            WRITE(*,*)'          No restart data file at myrank= ',myrank
            WRITE(unit_log,*)'          No restart data file at myrank= ',myrank
            STOP
         ENDIF        
         INITIAL_r=.FALSE. 
      ENDIF
    !
      DO    
         READ(unit_restart_rv,*) nout,restart_start,time,dt
         
    !
         IF(nrod_2d.gt.0)THEN
            DO k=1,nrod_2d
               READ(unit_restart_rv,*) wet_b(k),wet_t(k)
            ENDDO      
         ENDIF
    !
         IF(ncell_fuel_rod.gt.0.and.nr_2d.gt.0)THEN
            DO k=1,ncell_fuel_rod
               DO i=1,nr_2d
                  READ(unit_restart_rv,*) t_fuel(k,i)
               ENDDO
            ENDDO       
         ENDIF
    !
         IF(nz1d.gt.0.and.nr_1d.gt.0)THEN
            DO j=1,nz1d
               DO i=1,nr_1d
                  READ(unit_restart_rv,*)t_hts_1d(j,i)
               ENDDO
            ENDDO 
         ENDIF
!
         IF(ncell_fuel_rod.gt.0.and.nr_2d.gt.0)THEN
            DO k=1,ncell_fuel_rod
               READ(unit_restart_rv,*)p3d_cupid(k)
               IF(cplmaster.gt.0)then
                  i=cupid_cell_hts2d(k)
                  qvol_mas(i)=p3d_cupid(k)
               ENDIF   
            ENDDO       
         ENDIF
     !         
!
!........Stop reading at designated restarting point                                                                             
!         
         IF(restart.eq.2) THEN
            IF(restart_start.eq.itim_last)EXIT
         ELSE
            IF(restart_start.eq.itim)EXIT
         ENDIF
      ENDDO
      CLOSE(unit_restart_rv)          
      
   5001 FORMAT(2(i10,1x),2e20.10)   
   5002 FORMAT(5e20.10)    
   5003 FORMAT(e20.10,2(i10,1x))  
!
      RETURN
      END SUBROUTINE restart_read_rv
