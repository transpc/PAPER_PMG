!
!
      SUBROUTINE rv_power_master_to_cupid                  
!      
      USE Zrv_hts_2d   , ONLY: nr_2d,v_2f,qvol_2f0,qvol_2f,power_2d
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nrod_fuel_rod,num_fuel_rod
      USE Zcore        , ONLY: np
      USE zrv_ncell    , ONLY: p3d_cupid,nz_fuel_rod      
!
      IMPLICIT NONE
!
      INTEGER :: i,j,kk,jj
      REAL(8)::  q_tot
!
!.....Give volumetric power   
!     
      q_tot=0.0d0
      DO j=1,nr_2d
         DO i=1,ncell_fuel_rod
            jj=nz_fuel_rod(i)
            kk=nrod_fuel_rod(i)
            qvol_2f(i,j)=p3d_cupid(i)*qvol_2f0(i,j)
            qvol_2f(i,j)=qvol_2f(i,j)/num_fuel_rod(kk)/v_2f(jj,j)
            q_tot=q_tot+qvol_2f(i,j)*num_fuel_rod(kk)*v_2f(jj,j)
         ENDDO
      ENDDO
!
!.....calculate total power
!
      IF(np.gt.1) CALL allreducei_r1(q_tot)
      power_2d=q_tot
!
      END SUBROUTINE rv_power_master_to_cupid 
!
!======================================================================
!======================================================================
!
!
!
      SUBROUTINE rv_power_master_to_cupid_subchan                 
!      
      USE Zrv_hts_2d   , ONLY: nr_2d,v_2f,qvol_2f0,qvol_2f
      USE Zrv_hts_2d   , ONLY: qvol_norm_2d,power_2d,nqvol,qvol_time
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,nrod_fuel_rod,num_fuel_rod
      USE Zcore        , ONLY: np,myrank
      USE zrv_ncell    , ONLY: p3d_cupid,nz_fuel_rod      
      USE Ztimecon     , ONLY: time
      USE Zio_unit     , ONLY: unit_log
!
      IMPLICIT NONE
!
      INTEGER :: i,j,kk,jj
      INTEGER,SAVE :: k                          
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: q_tot
!
!.....Give volumetric power   
!     
      IF(initial)THEN
         initial=.FALSE.
         q_tot=0.0d0  
         k=1
         DO j=1,nr_2d
            DO i=1,ncell_fuel_rod
               jj=nz_fuel_rod(i)
               kk=nrod_fuel_rod(i)
               qvol_2f(i,j)=p3d_cupid(i)*qvol_2f0(i,j)
               qvol_2f(i,j)=qvol_2f(i,j)/num_fuel_rod(kk)/v_2f(jj,j)*qvol_norm_2d(kk,k)
               q_tot=q_tot+qvol_2f(i,j)*num_fuel_rod(kk)*v_2f(jj,j)
            ENDDO
         ENDDO
!
!........calculate total power
!
         IF(np.gt.1) CALL allreducei_r1(q_tot)
         power_2d=q_tot
!
         IF(myrank.eq.0) THEN
            WRITE(*,"(11x,a,e15.5,a)")  '--Total power at Initial= ',q_tot,' J/s'  
            WRITE(unit_log,"(11x,a,e15.5,a)") '--Total power at Initial= ',q_tot,' J/s' 
         ENDIF
      ENDIF
!
!.....Next Fuel Power step, nqvol
!
      IF(time.gt.qvol_time(k))THEN
         q_tot=0.0d0
         IF(k.lt.nqvol) k=k+1
         DO j=1,nr_2d
            DO i=1,ncell_fuel_rod
               jj=nz_fuel_rod(i)
               kk=nrod_fuel_rod(i)
               qvol_2f(i,j)=p3d_cupid(i)*qvol_2f0(i,j)
               qvol_2f(i,j)=qvol_2f(i,j)/num_fuel_rod(kk)/v_2f(jj,j)*qvol_norm_2d(kk,k)
               q_tot=q_tot+qvol_2f(i,j)*num_fuel_rod(kk)*v_2f(jj,j)
            ENDDO
         ENDDO
!
!........calculate total power
!
         IF(np.gt.1) CALL allreducei_r1(q_tot)
         power_2d=q_tot
!
         IF(myrank.eq.0) THEN
            WRITE(*,"(11x,a,e15.5,a)")  '--Total power at qvol_time(k)= ',q_tot,' J/s'  
            WRITE(unit_log,"(11x,a,e15.5,a)") '--Total power at qvol_time(k)= ',q_tot,' J/s' 
         ENDIF
      ENDIF
!
      END SUBROUTINE rv_power_master_to_cupid_subchan
