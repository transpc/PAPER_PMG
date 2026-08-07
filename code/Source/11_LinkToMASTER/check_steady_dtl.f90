!
      SUBROUTINE check_steady_dtl(dtemp_max)
!
      USE VOL_DATA    
      USE Zzone      ,ONLY:ncell_fluid,nzone
      USE Zcore      ,ONLY:np
      USE MASTER4    ,ONLY:dtemp_pass,mas_dtemp,count_dtemp
!      
      IMPLICIT NONE
!      
      INTEGER i
      INTEGER,SAVE :: count=0
      LOGICAL,SAVE :: initial=.true.
!
      REAL(8) dtemp_max,dtemp
      REAL(8) dtemp_max_core
      REAL(8),ALLOCATABLE,SAVE :: tlo(:)
!      
      IF(initial)THEN
         initial=.FALSE.
         ALLOCATE(tlo(ncell_fluid))
         tlo(1:ncell_fluid)=cell%tl(1:ncell_fluid)
      ENDIF
!      
      dtemp_max=0.0d0
      dtemp_max_core=0.0d0
      DO i=1,ncell_fluid
         dtemp=DABS(cell%tl(i)-tlo(i))
         dtemp_max=DMAX1(dtemp_max,dtemp)
         IF(nzone(i).eq.6)dtemp_max_core=DMAX1(dtemp_max_core,dtemp)
      ENDDO
      IF(np.gt.1)THEN
         CALL allreducei_max_r1(dtemp_max)
         CALL allreducei_max_r1(dtemp_max_core)
      ENDIF    
!
      IF(dtemp_max.lt.mas_dtemp)THEN
         count=count+1
      ELSE
         count=0
      ENDIF   
      IF(count.ge.count_dtemp)dtemp_pass=1        
!      
      tlo(1:ncell_fluid)=cell%tl(1:ncell_fluid)
!
      END SUBROUTINE check_steady_dtl                         
