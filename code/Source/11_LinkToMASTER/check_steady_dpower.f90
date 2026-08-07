!
      SUBROUTINE check_steady_dpower(dpower_max,power_max)
!
      USE Zrv_hts_2d   , ONLY: nr_2d,qvol_2f    
      USE Zrv_ncell    , ONLY: ncell_fuel_rod
      USE Zcore        , ONLY: np
      USE MASTER4      , ONLY: dpower_pass,mas_dpower,count_dpower 
      !Power ratio (THpower/NKpower)-ljr     
      USE Zconst1      , ONLY: cplmaster
      USE Zporous      , ONLY: qsum00
      USE Zrv_hts_2d   , ONLY: power_2d   

!      
      IMPLICIT NONE
!      
      INTEGER i,j
      INTEGER,SAVE :: count=0
      LOGICAL,SAVE :: initial=.TRUE.
!      
      REAL(8) power_max
      REAL(8) dpower,dpower_max
      REAL(8),ALLOCATABLE,SAVE::qvol_2f_old(:,:)
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(qvol_2f_old(ncell_fuel_rod,nr_2d))
         qvol_2f_old(:,:)=qvol_2f(:,:)
      ENDIF
      dpower_max=0.0d0
      power_max=0.0d0
      
      !origin-pik      
      if(cplmaster.eq.1)then
         DO i=1,ncell_fuel_rod
            DO j=1,nr_2d
               dpower=DABS(qvol_2f(i,j)-qvol_2f_old(i,j))
               dpower_max=DMAX1(dpower_max,dpower)
               power_max=DMAX1(power_max,qvol_2f(i,j))
            ENDDO
         ENDDO            
         IF(np.gt.1)THEN
            CALL allreducei_max_r1(dpower_max)
            CALL allreducei_max_r1(power_max)
         ENDIF
      !power ratio - subchannel-TH -ljr      
      else
         dpower=dabs(power_2d-qsum00)/power_2d*100.
         dpower_max=dpower
      endif
!
      IF(dpower_max.lt.mas_dpower)THEN
         count=count+1
      ELSE
         count=0
      ENDIF   
      IF(count.ge.count_dpower)dpower_pass=1        
!   
      qvol_2f_old(:,:)=qvol_2f(:,:)      
!      
      END SUBROUTINE check_steady_dpower  
!-----------------------------------------------------------------
      SUBROUTINE check_steady_dpower_master(dpower_max,power_max)
!
      USE Zrv_ncell    , ONLY: ncell_fuel_rod,p3d_cupid
      USE Zcore        , ONLY: np
      USE MASTER4      , ONLY: dpower_pass,mas_dpower,count_dpower      
!      
      IMPLICIT NONE
!      
      INTEGER i
      INTEGER,SAVE::count
!      
      LOGICAL,SAVE::initial
!      
      REAL(8) power_max
      REAL(8) dpower,dpower_max
      REAL(8),ALLOCATABLE,SAVE::p3d_cupid_old(:)
!      
      DATA count,initial/0,.TRUE./
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(p3d_cupid_old(ncell_fuel_rod))
         p3d_cupid_old(:)=p3d_cupid(:)
      ENDIF
      dpower_max=0.0d0
      power_max=0.0d0
      DO i=1,ncell_fuel_rod
          dpower=DABS(p3d_cupid(i)-p3d_cupid_old(i))
          dpower_max=DMAX1(dpower_max,dpower)
          power_max=DMAX1(power_max,p3d_cupid(i))
      ENDDO            
      IF(np.gt.1)THEN
         CALL allreducei_max_r1(dpower_max)
         CALL allreducei_max_r1(power_max)
      ENDIF 
!
      IF(dpower_max.lt.mas_dpower)THEN
         count=count+1
      ELSE
         count=0
      ENDIF   
      IF(count.ge.count_dpower)dpower_pass=1        
!   
      p3d_cupid_old(:)=p3d_cupid(:)      
!      
      END SUBROUTINE check_steady_dpower_master
