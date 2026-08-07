      SUBROUTINE rv_cthxpr(ctemp,cdthex) 
!
!     calculates diametral thermal expansion of zircaloy cladding       
!
      IMPLICIT NONE
!                                                                       
!     input                                                             
!          ctemp     -  cladding temperature (k)                        
!     output                                                            
!          cdthex    -  diametral thermal expansion of zircaloy (m/m)   
!                                                                       
!     taken from matpro-11, revision 1                                  
!                                                                       
      INTEGER i1,i2,i
      REAL(8) dthexp(42),cdthex 
      REAL(8) ctemp
!                                                                       
      DATA dthexp/5.1395d-03,1073.15d0,5.2200d-03,1083.15d0,5.2500d-03, &
      1093.15d0,5.2800d-03,1103.15d0,5.2800d-03,1113.15d0,5.2400d-03,   &
      1123.15d0,5.2200d-03,1133.15d0,5.1500d-03,1143.15d0,5.0800d-03,   &
      1153.15d0,4.9000d-03,1163.15d0,4.7000d-03,1173.15d0,4.4500d-03,   &
      1183.15d0,4.1000d-03,1193.15d0,3.5000d-03,1203.15d0,3.1300d-03,   &
      1213.15d0,2.9700d-03,1223.15d0,2.9200d-03,1233.15d0,2.8700d-03,   &
      1243.15d0,2.8600d-03,1253.15d0,2.8800d-03,1263.15d0,2.9000d-03,   &
      1273.15d0/                                                        
!                                                                       
!                                                                       
      IF(ctemp.le.dthexp(2)) GOTO 200 
      IF(ctemp.ge.dthexp(42)) GOTO 300 
!                                                                       
      IF(ctemp.gt.dthexp(14)) GOTO 104 
      i1=4 
      i2=14 
      GOTO 112 
  104 IF(ctemp.gt.dthexp(28)) GOTO 108 
      i1=16 
      i2=28 
      GOTO 112 
  108 i1=30 
      i2=42 
  112 DO 116 i=i1,i2,2 
         IF(ctemp.lt.dthexp(i)) GOTO 120 
  116 END DO 
      cdthex=dthexp(i2-1) 
      GOTO 400 
  120 cdthex=dthexp(i-3)+(ctemp-dthexp(i-2))*(dthexp(i-1)-dthexp(i-3))/(&
      dthexp(i)-dthexp(i-2))                                            
      GOTO 400 
!                                                                       
  200 cdthex=-2.3730d-04+(ctemp-273.15d0) * 6.7210d-06   !MARS code.  USE THIS!
!  200 cdthex=1.5985d-03+(ctemp-273.15d0) * 6.7210d-06   !MARS report  
      GOTO 400 
!                                                                       
  300 cdthex=-6.800d-03+(ctemp-273.15d0) * 9.70d-06      !MARS code.  USE THIS!
!  300 cdthex=-4.150d-03+(ctemp-273.15d0) * 9.70d-06     !MARS report  
!  
  400 RETURN 
      END SUBROUTINE rv_cthxpr