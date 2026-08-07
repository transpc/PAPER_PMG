      Function linear_power(x1,x2,y1,y2,x)
!      
      IMPLICIT NONE
!      
      REAL(8) x1,x2,y1,y2,x,linear_power
!      
      linear_power=(y2-y1)/(x2-x1)*(x-x1)+y1
!      
      RETURN
!      
      ENDFunction linear_power
!===============================================================================      
      Function heat_factor(t1)
!      
      USE Zconst1,  ONLY:vv_prob
!      
      IMPLICIT NONE
!      
      INTEGER i
!      
      REAL(8) powerfactor
      REAL(8) t1,hf,linear_power,heat_factor      
!
      LOGICAL,SAVE::initial
      DATA initial/.true./           
!
      INTEGER,SAVE :: interval_max
      REAL(8),ALLOCATABLE:: hfao(:),tao(:) !atlas-next      
      REAL(8),SAVE,ALLOCATABLE::hfa(:),ta(:)
!DEC$IF defined(mpi_flag)
      INCLUDE 'mpif.h'
!DEC$ELSE
!DEC$ENDIF         
      !       
      IF(initial)then
         initial=.FALSE.
         IF(vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
            interval_max=38     
            ALLOCATE(tao(interval_max),hfao(interval_max),ta(interval_max),hfa(interval_max)) 
            tao( 1)=0.00000
            tao( 2)=0.50777
            tao( 3)=1.00777
            tao( 4)=2.00777
            tao( 5)=3.00777
            tao( 6)=4.00777
            tao( 7)=5.00777
            tao( 8)=6.00777
            tao( 9)=7.00777
            tao(10)=8.00777
            tao(11)=9.00777
            tao(12)=10.0078
            tao(13)=15.5082
            tao(14)=20.0080
            tao(15)=25.0082
            tao(16)=30.0079
            tao(17)=35.0079
            tao(18)=40.0079
            tao(19)=45.0082
            tao(20)=50.0080
            tao(21)=60.0076
            tao(22)=70.0082
            tao(23)=80.0081
            tao(24)=90.0081
            tao(25)=99.5082
            tao(26)=100.008
            tao(27)=120.008
            tao(28)=140.008
            tao(29)=160.508
            tao(30)=180.008
            tao(31)=200.010
            tao(32)=220.510
            tao(33)=240.008
            tao(34)=260.509
            tao(35)=280.509
            tao(36)=300.001
            tao(37)=900.000
            tao(38)=5000.00
            hfao( 1)=0.0695225   
            hfao( 2)=0.06904875  
            hfao( 3)=0.067082    
            hfao( 4)=0.063692    
            hfao( 5)=0.06106775  
            hfao( 6)=0.05899475  
            hfao( 7)=0.05730075  
            hfao( 8)=0.05587775  
            hfao( 9)=0.054662    
            hfao(10)=0.0536035   
            hfao(11)=0.05266175  
            hfao(12)=0.05180775  
            hfao(13)=0.0482975   
            hfao(14)=  0.046319  
            hfao(15)=0.04460975  
            hfao(16)=0.04324175  
            hfao(17)=0.04209975  
            hfao(18)=0.041118    
            hfao(19)=0.04025625  
            hfao(20)=  0.03949   
            hfao(21)=0.03817475  
            hfao(22)=0.03707775  
            hfao(23)=0.0361435   
            hfao(24)=0.03533425  
            hfao(25)=0.03465725  
            hfao(26)=0.03462375  
            hfao(27)=0.033425    
            hfao(28)=0.03244925  
            hfao(29)=0.03161125  
            hfao(30)=0.0309285   
            hfao(31)=0.030317    
            hfao(32)=0.0297625   
            hfao(33)=0.029292    
            hfao(34)=0.02884425  
            hfao(35)=0.028446    
            hfao(36)=0.02808875  
            hfao(37)=0.02        
            hfao(38)=0.018       
         ELSEIF(vv_prob.eq.'atlas_mc_porous')THEN
            interval_max=19      
            ALLOCATE(tao(interval_max),hfao(interval_max),ta(interval_max),hfa(interval_max)) 
            tao( 1)=0.0
            tao( 2)=12.1
            tao( 3)=14.0
            tao( 4)=21.8
            tao( 5)=31.0
            tao( 6)=38.8
            tao( 7)=55.7
            tao( 8)=90.4
            tao( 9)=189.4
            tao(10)=270.7
            tao(11)=341.4
            tao(12)=443.9
            tao(13)=592.4
            tao(14)=793.9
            tao(15)=991.9
            tao(16)=1494.0
            tao(17)=1992.5
            tao(18)=2993.0
            tao(19)=5001.2       
            hfao( 1)=1.0000 
            hfao( 2)=1.000000
            hfao( 3)=0.953750
            hfao( 4)=0.827500
            hfao( 5)=0.737500
            hfao( 6)=0.690000
            hfao( 7)=0.620000
            hfao( 8)=0.542500
            hfao( 9)=0.440000
            hfao(10)=0.400000
            hfao(11)=0.378750  
            hfao(12)=0.358750
            hfao(13)=0.338750
            hfao(14)=0.318750
            hfao(15)=0.301250
            hfao(16)=0.268750
            hfao(17)=0.243750
            hfao(18)=0.211250
            hfao(19)=0.180000           
         ELSE
            WRITE(*,"(11x,a,a,a)")'hfactor is not assigned in ',vv_prob,'!!!'
            PAUSE
            STOP
         ENDIF
         powerfactor=1.0d0
         !IF(t1.lt.400)powerfactor=1.2d0 !atlas-mslb-sensitivity
         IF(1)THEN !original heat
            DO i=1,interval_max
               ta(i)=tao(i)
               hfa(i)=hfao(i)*powerfactor !run_power_1.0~1.3
            ENDDO   
         ELSE      !modified heat to fit pressure shape
            DO i=1,interval_max
               ta(i)=tao(i) 
               hfa(i)=hfao(i)
               IF(i.le.14)then
                  hfa(i)=hfao(i)*0.90d0 !atlas-mslb
               ELSE   
                  hfa(i)=0.31613*0.80d0
               ENDIF   
            ENDDO 
         ENDIF
         DEALLOCATE(tao,hfao)
      ENDIF
!      
      IF(t1.le.0.0d0)then
         hf=hfa(1) !1.0d0
      ELSEIF(t1.le.ta(1))then 
         hf=hfa(1)  
      ELSEIF(t1.gt.ta(1).and.t1.le.ta(interval_max))then   
          DO i=1,interval_max-1
             IF(t1.gt.ta(i).and.t1.le.ta(i+1))then   
                hf=linear_power(ta(i),ta(i+1),hfa(i),hfa(i+1),t1)
                EXIT
             ENDIF   
          ENDDO
      ELSE
         hf=hfa(interval_max) 
      ENDIF
!      
      heat_factor=hf  
!      
      RETURN    
      ENDFunction heat_factor    
