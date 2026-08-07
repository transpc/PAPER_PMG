!
      SUBROUTINE print_ss_status(q_tot,q_por,i_opt,dtemp_max,dpower_max,dmpower_max,dmass_max,dmaster_pass,power_max,mpower_max)
!
      USE Zcore      ,ONLY: myrank
      USE Ztimecon   ,ONLY: time
      USE Zconst2    ,ONLY: dt 
      USE MASTER4    ,ONLY: dpres,mas_dmass_opt
!      
      IMPLICIT NONE
      INCLUDE '../10_LinkToMARS/c3com.h' !n_volleg,n_junleg 
!
      CHARACTER*14 varname(20)      
!
      INTEGER i,j
      INTEGER i_opt,dmaster_pass 
!
      LOGICAL,SAVE :: initial=.true.
!      
      REAL(8) q_tot,q_por
      REAL(8) dtemp_max,dpower_max,dmpower_max,dmass_max,power_max,mpower_max
!      
      DATA varname/'Time(s)','dTl(K)','dMass(kg/s)','dPower(W)','dmPower(W)','dt(s)','tPower(MW)','fPower(MW)','MasterMode','AchieveSteady', &
                   'Power(W)','mPower(W)','dpMARS_1a','dpMARS_1b','dpMARS_2a','dpMARS_2b','dpCUP_1a','dpCUP_1b','dpCUP_2a','dpCUP_2b'/       
!
      IF(initial)THEN
         initial=.FALSE.
         IF(myrank.eq.0)THEN
            OPEN(812,file='ss_status.dat')
            IF(mas_dmass_opt.gt.0)THEN
               WRITE(812,"(20a14)")(varname(i),i=1,20)
            ELSE
               WRITE(812,"(20a14)")(varname(i),i=1,12)
            ENDIF   
         ENDIF
      ENDIF
!         
      IF(myrank.eq.0)THEN
         IF(mas_dmass_opt.gt.0)THEN
            WRITE(812,"(6e14.6,2f14.6,2i14,10e14.6)")time,dtemp_max,dmass_max,dpower_max,dmpower_max,dt,q_tot/1.d6,q_por/1.d6,i_opt,dmaster_pass,power_max,mpower_max,(dpres(j),j=1,8)     
         ELSE
            WRITE(812,"(6e14.6,2f14.6,2i14, 2e14.6)")time,dtemp_max,dmass_max,dpower_max,dmpower_max,dt,q_tot/1.d6,q_por/1.d6,i_opt,dmaster_pass,power_max,mpower_max
         ENDIF   
      ENDIF   
!      
      END SUBROUTINE print_ss_status
    
