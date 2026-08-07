!
      SUBROUTINE rv_power_transient                  
!      
      USE Zrv_model    , ONLY: rv_ht_str
      USE Zrv_hts_2d   , ONLY: nr_2d,v_2f,qvol_2f0,qvol_2f,qvol_norm_2d,nqvol,qvol_time,power_2d, &
                               l_ht_str_2d_qcell
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nrod_fuel_rod,num_fuel_rod,nz_fuel_rod
      USE Ztimecon     , ONLY: time
      USE Zcore        , ONLY: np,myrank 
      USE zconst1      , ONLY: cplmaster      
      USE Zporous      , ONLY: l_subchannel
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,kk,jj
      INTEGER,SAVE :: iqvol                          
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: q_tot
!
!.....MASTER Code coupling
!
      IF(cplmaster.gt.0)THEN
        CALL rv_power_transient_master
        RETURN
      ENDIF
!
      IF(rv_ht_str.ne.1)RETURN      
!
!.....Text-based subchannel power input (NOT used MASTER coupling)
!
!
      IF(l_ht_str_2d_qcell)THEN
        CALL rv_power_from_input
        RETURN
      ENDIF
!
      IF(l_subchannel)THEN
         CALL rv_power_master_to_cupid_subchan
         RETURN
      ENDIF
!
      IF(initial)THEN
         initial=.FALSE.
         DO iqvol=1,nqvol
            IF(qvol_time(iqvol).ge.time)EXIT 
         ENDDO        
         IF(iqvol.gt.nqvol)THEN
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Check ht_str_2d.in!!!'
            PAUSE
            STOP
         ELSE
            IF(myrank.eq.0)WRITE(*,"(11x,a,2i5,1f8.1,a)")'--nqvol,iqvol,time=',nqvol,iqvol,qvol_time(iqvol),'s.'
         ENDIF
                   
         q_tot=0.0d0  
         DO j=1,nr_2d
            DO i=1,ncell_fuel_rod
               jj=nz_fuel_rod(i)
               kk=nrod_fuel_rod(i)  
               qvol_2f0(i,j)=qvol_2f0(i,j)/v_2f(jj,j)
               qvol_2f(i,j)=qvol_2f0(i,j)*qvol_norm_2d(kk,iqvol)
               q_tot=q_tot+qvol_2f(i,j)*v_2f(jj,j)*num_fuel_rod(kk)               
            ENDDO
         ENDDO
         IF(np.gt.1) CALL allreducei_r1(q_tot)
         power_2d=q_tot
         IF(myrank.eq.0) THEN
            WRITE(*,"(11x,a,1f12.6,a,1f14.6,a)")'--Power=',q_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
            WRITE(unit_log,"(11x,a,1f12.6,a,1f8.1,a)")'--Power=',q_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
         ENDIF
      ENDIF
!
      IF(iqvol.lt.nqvol)THEN
         IF(time.gt.qvol_time(iqvol))THEN
            iqvol=iqvol+1
            q_tot=0.0d0
            DO j=1,nr_2d
               DO i=1,ncell_fuel_rod
                  jj=nz_fuel_rod(i)
                  kk=nrod_fuel_rod(i)
                  qvol_2f(i,j)=qvol_2f0(i,j)*qvol_norm_2d(kk,iqvol)
                  q_tot=q_tot+qvol_2f(i,j)*v_2f(jj,j)*num_fuel_rod(kk)
               ENDDO
            ENDDO
            IF(np.gt.1) CALL allreducei_r1(q_tot)
            power_2d=q_tot 
            IF(myrank.eq.0) THEN
               WRITE(*,"(11x,a,1f12.6,a,1f6.1,a)")'--1Power=',q_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
               WRITE(unit_log,"(11x,a,1f12.6,a,1f6.1,a)")'--1Power=',q_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
            ENDIF
         ENDIF  
      ENDIF  
!            
      END SUBROUTINE rv_power_transient
!----------------------------------------------------------------------------
      SUBROUTINE rv_power_from_input
      USE Zrv_hts_2d   , ONLY: nr_2d,v_2f,qvol_2f0,qvol_2f,nqvol,qvol_time,power_2d,&
                               nz0_2d,qcell_to_qvol_mul
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nrod_fuel_rod,num_fuel_rod,nz_fuel_rod,num_ch,nz_fine
      USE Ztimecon     , ONLY: time
      USE Zcore        , ONLY: np,myrank 
      USE Zrv_mpi      , ONLY: jperm_fuel_rod      
      USE Zmars        , ONLY: time_scram
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,k,kk,jj,m,ii  
      INTEGER,SAVE :: iqvol                        
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: q_tot,qcell_tot,ftmp 
      REAL(8),DIMENSION(:),ALLOCATABLE :: qcell_tmp,qcell,qcell_row
!           
      IF(initial)THEN
         initial=.FALSE.
         OPEN(53,file='ht_str_2d_qcell.in')
         READ(53,*)nqvol
         IF(ALLOCATED(qvol_time))DEALLOCATE(qvol_time)
         ALLOCATE(qvol_time(0:nqvol))
         qvol_time(:)=0.0d0
         DO i=1,nqvol
             READ(53,*)ftmp
             ftmp=ftmp+time_scram
             IF(ftmp.ge.time)THEN
                BACKSPACE(53)
                iqvol=i-1
                EXIT
             ENDIF   
             DO j=1,num_ch 
                READ(53,*)ftmp
             ENDDO
             iqvol=i
         ENDDO
       
         IF(iqvol.ge.nqvol)THEN
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Check ht_str_2d_qcell.in!!!'
            PAUSE
            STOP
         ELSE
            IF(myrank.eq.0)WRITE(*,"(11x,a,2i5,1f7.1,a)")'--nqvol,iqvol,time=',nqvol,iqvol,ftmp,'s.'
         ENDIF
      ENDIF   
!
      IF(iqvol.lt.nqvol)THEN
         IF(time.ge.qvol_time(iqvol))THEN
            iqvol=iqvol+1      
            ALLOCATE(qcell_tmp(num_ch*nz0_2d*nz_fine),qcell(ncell_fuel_rod),qcell_row(nz0_2d))
            READ(53,*)qvol_time(iqvol)
            qvol_time(iqvol)=qvol_time(iqvol)+time_scram
            qcell_tot=0.0d0
            m=0
            DO i=1,num_ch  
                READ(53,*)(qcell_row(j),j=1,nz0_2d)
                DO j=1,nz0_2d
                   DO k=1,nz_fine
                      m=m+1
                      qcell_tmp(m)=qcell_row(j)/nz_fine
                      qcell_tot=qcell_tot+qcell_tmp(m)
                   ENDDO
                ENDDO
            ENDDO 
            DO ii=1,ncell_fuel_rod
               i=jperm_fuel_rod(ii)
               qcell(ii)=qcell_tmp(i)
            ENDDO   
            DO ii=1,ncell_fuel_rod
               i=jperm_fuel_rod(ii)
               qvol_2f0(ii,:)=qcell_to_qvol_mul(ii,:)*qcell(ii)
            ENDDO
            DEALLOCATE(qcell_tmp,qcell,qcell_row)          
   !      
            q_tot=0.0d0
            DO j=1,nr_2d
               DO i=1,ncell_fuel_rod
                  jj=nz_fuel_rod(i)
                  kk=nrod_fuel_rod(i)
                  qvol_2f(i,j)=qvol_2f0(i,j)/v_2f(jj,j) !*qvol_norm_time(k)
                  q_tot=q_tot+qvol_2f(i,j)*v_2f(jj,j)*num_fuel_rod(kk)
               ENDDO
            ENDDO
            IF(np.gt.1) CALL allreducei_r1(q_tot)
            power_2d=q_tot 
            IF(myrank.eq.0) THEN
               WRITE(*,"(11x,a,2f8.1,a,1f6.1,a)")'--Power and qcell=',q_tot/1.0d6,qcell_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
               WRITE(unit_log,"(11x,a,2f8.1,a,1f6.1,a)")'--Power and qcell=',q_tot/1.0d6,qcell_tot/1.0d6,' MW by',qvol_time(iqvol),' s.'
            ENDIF
         ENDIF
      ELSE
         IF(k.gt.nqvol)CLOSE(53)
      ENDIF   
!      
      END SUBROUTINE rv_power_from_input
