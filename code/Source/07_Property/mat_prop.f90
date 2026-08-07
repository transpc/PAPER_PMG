!
      SUBROUTINE mat_prop(NoMaterial,temp,CpVol,Condu,iOKr,iOKk)
!
!     Calculates the material properties (k & Cp)
!     Input: NoMaterial,temp (K)
!     Output: CpVol(J/m3.K),Condu(W/m.K),iOKr,iOKk
!
!     iOK=-1: temp. is lower   than the minimum temp.
!     iOK=+1: temp. is greater than the maximum temp.
!
      IMPLICIT NONE
!   
      INTEGER  i1,j1,i
!      INTEGER (4)MaxData(2,5)
      INTEGER (4)MaxData(2,9)
      INTEGER (4)NoMaterial,No_R,No_K,iOKr,iOKk   
      REAL(8) temp,CpVol,Condu,ttemp
      LOGICAL, SAVE:: initial
      DATA initial/.true./
!
!     Array for material tables
!     TR: Temp (K) vs. Volumetric heat capacity (J/m3.K)
!     TK: Temp (K) vs. Thermal conductivity (W/m.K)
!
!     1. UO2
!     2. Zircaloy
!     3. Inconel
!     4. Stainless steel
!     5. Carbon steel
!     6. BN
!     7. Monel400: Nichrome
!     8. Inconel 600
!
      REAL(8) TK1(2,32),TR1(2,26),TK2(2,16),TR2(2,15),TK3(2,7),TR3(2,7)
      REAL(8) TK4(2,8),TR4(2,8),TK5(2,8),TR5(2,8)
      REAL(8) TK6(2,10),TR6(2,10)
      REAL(8) TK7(2,10),TR7(2,10)
      REAL(8) TK8(2,10),TR8(2,10)
      REAL(8) TK9(2,2),TR9(2,4)
!
      REAL(8),SAVE:: TTK1(32,9),TTK2(32,9),DTK(31,9)
      REAL(8),SAVE:: TTR1(32,9),TTR2(32,9),DTR(31,9)
!
!      DATA MaxData/32,26,16,15,7,7,8,8,8,8/
      DATA MaxData/32,26,16,15,7,7,8,8,8,8,10,10,10,10,10,10,2,4/      
!
!
!-----------------------------------------------------------------------
!uo-2 (01)
!
!     temp[k] thermal conductivity (W/m.K)
!32
      DATA TK1/273.15d0,7.293d0,373.15d0,7.293d0,473.15d0,6.697d0,          &
      573.15d0,5.815d0,673.14d0,5.136d0,773.15d0,4.603d0,873.15d0,          &
      4.172d0,973.15d0,3.819d0,1073.15d0,3.527d0,1173.15d0,3.281d0,         &
      1273.15d0,3.075d0,1373.15d0,2.903d0,1473.15d0,2.759d0,1573.15d0,      &
      2.641d0,1673.15d0,2.546d0,1773.15d0,2.474d0,1873.15d0,2.423d0,        &
      1973.15d0,2.392d0,2073.15d0,2.382d0,2173.15d0,2.392d0,2273.15d0,      &
      2.422d0,2373.15d0,2.473d0,2473.15d0,2.546d0,2573.15d0,2.643d0,        &
      2673.15d0,2.762d0,2773.15d0,2.907d0,2873.15d0,3.078d0,2973.15d0,      &
      3.471d0,3073.15d0,3.504d0,3173.15d0,3.762d0,3573.15d0,5.133d0,        &
      4873.15d0,14.70d0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!26
      DATA TR1/273.1d0,2.427d6,400.0d0,2.754d6,500.0d0,2.927d6,600.0d0,     &
      3.043d6,700.0d0,3.139d6,800.0d0,3.178d6,900.0d0,3.236d6,              &
      1000.0d0,3.274d6,1100.0d0,3.313d6,1200.0d0,3.351d6,1300.0d0,3.378d6,  &
      1400.0d0,3.428d6,1500.0d0,3.459d6,1600.0d0,3.502d6,                   &
      1700.0d0,3.582d6,1800.0d0,3.660d6,1900.0d0,3.775d6,2000.0d0,3.992d6,  &
      2100.0d0,4.169d6,2200.0d0,4.366d6,2300.0d0,4.622d6,                   &
      2400.0d0,4.897d6,2500.0d0,5.212d6,2600.0d0,5.585d6,3000.0d0,7.395d6,  &
        4873.1d0,16.00d6/
!
!-----------------------------------------------------------------------
!zircaloy (02)
!
!     temp[k] thermal conductivity (W/m.K)
!16
      DATA TK2/273.15d0,13.6d0,373.15d0,14.1d0,473.15d0,14.8d0,573.15d0,    &
      15.8d0,673.15d0,16.9d0,773.15d0,18.1d0,873.15d0,19.5d0,               &
      973.15d0,21.1d0,1073.15d0,22.8d0,1173.15d0,24.6d0,1273.15d0,26.8d0,   &
      1373.15d0,29.2d0,1473.15d0,31.7d0,1573.15d0,34.4d0,                   &
      1673.15d0,37.3d0,1773.15d0,40.4d0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!15
      DATA TR2/273.15d0,1.881d6,573.15d0,2.079d6,773.15d0,2.211d6,903.15d0, &
      2.290d6,923.15d0,2.376d6,1083.15d0,2.376d6,1103.15d0,                 &
      3.630d6,1123.15d0,4.455d6,1143.15d0,4.950d6,1163.15d0,5.115d6,        &
      1183.15d0,4.950d6,1203.15d0,4.455d6,1213.15d0,3.360d6,1243.15d0,      &
      2.376d6,2073.15d0,2.376d6/
!     
!-----------------------------------------------------------------------
!inconel (03)
!
!     temp (k) thermal conductivity (W/m.K)
!7
      DATA TK3/0.2731500d3,0.1318846d2,0.2942611d3,0.1485516d2,0.3664833d3, &
      0.1572060d2,0.4775944d3,0.1745086d2,0.5887055d3,                      &
      0.1918174d2,0.6998167d3,0.2091262d2,0.8387055d3,0.2516817d2/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!
      DATA TR3/0.2731500d3,0.3746931d7,0.2942611d3,0.3746931d7,0.3664833d3, &
      0.3923651d7,0.4775944d3,0.4100438d7,0.5887055d3,                      &
      0.4276822d7,0.6998167d3,0.4453877d7,0.8387055d3,0.4553877d7/
!
!-----------------------------------------------------------------------
!stainless steel (04)
!
!     temp (k)   thermal conductivity (W/m.K)
!8
      DATA TK4/0.2942611d3,0.1488507d2,0.3664833d3,0.1609382d2,0.4775944d3, &
      0.1800040d2,0.5887055d3,0.1955807d2,0.6998167d3,0.2111574d2,          &
      0.8109278d3,0.2284786d2,0.9220389d3,0.2423107d2,0.1088706d4,          &
      0.2648035d2/
!
!     temp (k)   volumdtric heat capacity (J/m3.K)
!
      DATA TR4/0.2942611d3,0.3819698d7,0.3664833d3,0.3998161d7,0.4775944d3, &
      0.4227193d7,0.5887055d3,0.4355491d7,0.6998167d3,0.4446768d7,          &
      0.8109278d3,0.4563263d7,0.9220389d3,0.4625232d7,0.1088706d4,          &
      0.4750512d7/
!
!-----------------------------------------------------------------------
!carbon steel (05)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK5/0.2942611d3,0.3772982d2,0.3664833d3,0.3876847d2,0.4775944d3, &
      0.3859277d2,0.5887055d3,0.3720956d2,0.6998167d3,0.3530920d2,          &
      0.8109278d3,0.3322816d2,0.9220389d3,0.3080443d2,0.1088706d4,          &
      0.2596133d2/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR5/0.2942611d3,0.3480744d7,0.3664833d3,0.3765106d7,0.4775944d3, &
      0.4108486d7,0.5887055d3,0.4436440d7,0.6998167d3,0.4174881d7,          &
      0.8109278d3,0.5277385d7,0.9220389d3,0.6282777d7,0.1088706d4,          &
      0.5322387d7/
!
!-----------------------------------------------------------------------
! BN (06)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK6/0.293e+03,   0.0720e+02, &
               0.373e+03,   0.0749e+02, &
               0.473e+03,   0.0703e+02, &
               0.573e+03,   0.0574e+02, &
               0.673e+03,   0.0515e+02, &
               0.773e+03,   0.0469e+02, &
               0.873e+03,   0.0469e+02, &
               0.973e+03,   0.0452e+02, &
               0.1073e+04,  0.0417e+02, &
               0.1173e+04,  0.0366e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR6/0.293e+03,   1.434e+06, &  
               0.373e+03,   1.959e+06, &  
               0.473e+03,   2.335e+06, &  
               0.573e+03,   2.610e+06, &  
               0.673e+03,   2.805e+06, &  
               0.773e+03,   2.975e+06, &  
               0.873e+03,   3.106e+06, &  
               0.973e+03,   3.217e+06, &  
               0.1073e+04,  3.279e+06, &  
               0.1173e+04,  3.325e+06/
!
!-----------------------------------------------------------------------
! Monel: Nichrome (07)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK7/0.293e+03,   0.134e+02, & 
               0.373e+03,   0.143e+02, & 
               0.473e+03,   0.155e+02, & 
               0.573e+03,   0.176e+02, & 
               0.673e+03,   0.188e+02, & 
               0.773e+03,   0.208e+02, & 
               0.873e+03,   0.242e+02, & 
               0.973e+03,   0.276e+02, & 
               0.1073e+04,  0.263e+02, & 
               0.1173e+04,  0.288e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR7/0.293e+03,   3.687e+06, & 
               0.373e+03,   3.809e+06, & 
               0.473e+03,   3.933e+06, & 
               0.573e+03,   4.023e+06, & 
               0.673e+03,   4.135e+06, & 
               0.773e+03,   4.285e+06, & 
               0.873e+03,   4.678e+06, & 
               0.973e+03,   4.746e+06, & 
               0.1073e+04,  4.790e+06, & 
               0.1173e+04,  4.821e+06/               
!
!-----------------------------------------------------------------------
! Inconel 600 (08)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK8/0.293e+03,   0.143e+02, &  
               0.373e+03,   0.158e+02, &  
               0.473e+03,   0.178e+02, &  
               0.573e+03,   0.190e+02, &  
               0.673e+03,   0.200e+02, &  
               0.773e+03,   0.211e+02, &  
               0.873e+03,   0.261e+02, &  
               0.973e+03,   0.278e+02, &  
               0.1073e+04,  0.268e+02, &  
               0.1173e+04,  0.281e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR8/0.293e+03,   3.718e+06, &   
               0.373e+03,   3.917e+06, &   
               0.473e+03,   4.076e+06, &   
               0.573e+03,   4.190e+06, &   
               0.673e+03,   4.247e+06, &   
               0.773e+03,   4.054e+06, &   
               0.873e+03,   4.759e+06, &   
               0.973e+03,   4.800e+06, &   
               0.1073e+04,  4.846e+06, &   
               0.1173e+04,  4.846e+06/   
!-----------------------------------------------------------------------
!Soil for ATLAS_CUBE (09)
!
!     temp (k)   thermal conductivity (W/m.K)
!8
      DATA TK9/0.2942611d3,0.529d0,0.3664833d3,0.529d0/
!
!     temp (k)   volumdtric heat capacity (J/m3.K)
!
      DATA TR9/0.2942611d3,1.666D6,0.373d3,1.666D6,0.473d3,1.873D6,1.273d3,4.015D6/
!---------------------------------------------------------------------------------
!
      IF(initial)THEN
         initial=.false.
         DO i1=1,9
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)
            IF(i1.eq.1) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR1(1,j1)
                  TTR2(j1,i1)=TR1(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK1(1,j1)
                  TTK2(j1,i1)=TK1(2,j1)
               ENDDO
            ELSEIF(i1.eq.2) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR2(1,j1)
                  TTR2(j1,i1)=TR2(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK2(1,j1)
                  TTK2(j1,i1)=TK2(2,j1)
               ENDDO
            ELSEIF(i1.eq.3) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR3(1,j1)
                  TTR2(j1,i1)=TR3(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK3(1,j1)
                  TTK2(j1,i1)=TK3(2,j1)
               ENDDO
            ELSEIF(i1.eq.4) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR4(1,j1)
                  TTR2(j1,i1)=TR4(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK4(1,j1)
                  TTK2(j1,i1)=TK4(2,j1)
               ENDDO
            ELSEIF(i1.eq.5) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR5(1,j1)
                  TTR2(j1,i1)=TR5(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK5(1,j1)
                  TTK2(j1,i1)=TK5(2,j1)
               ENDDO
            ELSEIF(i1.eq.6) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR6(1,j1)
                  TTR2(j1,i1)=TR6(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK6(1,j1)
                  TTK2(j1,i1)=TK6(2,j1)
               ENDDO
            ELSEIF(i1.eq.7) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR7(1,j1)
                  TTR2(j1,i1)=TR7(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK7(1,j1)
                  TTK2(j1,i1)=TK7(2,j1)
               ENDDO
            ELSEIF(i1.eq.8) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR8(1,j1)
                  TTR2(j1,i1)=TR8(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK8(1,j1)
                  TTK2(j1,i1)=TK8(2,j1)
               ENDDO
            ELSEIF(i1.eq.9) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR9(1,j1)
                  TTR2(j1,i1)=TR9(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK9(1,j1)
                  TTK2(j1,i1)=TK9(2,j1)
               ENDDO               
            ENDIF
         ENDDO
!
         DO i1=1,9
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)
            DO i=2,No_R
               DTR(i-1,i1)=(TTR2(i,i1)-TTR2(i-1,i1))/(TTR1(i,i1)-TTR1(i-1,i1))
            ENDDO 
            DO i=2,No_K
               DTK(i-1,i1)=(TTK2(i,i1)-TTK2(i-1,i1))/(TTK1(i,i1)-TTK1(i-1,i1))
            ENDDO 
         ENDDO 
!
      ENDIF
!-----------------------------------------------------------------------
!
         i1=NoMaterial
         No_K=MaxData(1,i1)
         No_R=MaxData(2,i1)
!
!.....Volumetric heat capacity
!
      ttemp=temp
      iOKr=0
      IF(ttemp.le.TTR1(1,i1))THEN
         ttemp=TTR1(1,i1)
         iOKr=-1
      ELSEIF(ttemp.ge.TTR1(No_R,i1))THEN
         ttemp=TTR1(No_R,i1)
         iOKr=+1
      ENDIF
!
      DO j1=2,No_R
         IF(ttemp.le.TTR1(j1,i1))EXIT
      ENDDO
!
!     CpVol=TTR2(j1-1,i1)+(TTR2(j1,i1)-TTR2(j1-1,i1))/(TTR1(j1,i1)-TTR1(j1-1,i1))*(ttemp-TTR1(j1-1,i1))
      CpVol=TTR2(j1-1,i1)+DTR(j1-1,i1)*(ttemp-TTR1(j1-1,i1))
!
!.....Thermal conductivity
!
      ttemp=temp
      iOKk=0
      IF(ttemp.le.TTK1(1,i1))THEN
         ttemp=TTK1(1,i1)
         iOKk=-1
      ELSEIF(ttemp.ge.TTK1(No_K,i1))THEN
         ttemp=TTK1(No_K,i1)
         iOKk=+1
      ENDIF
!
      DO i=2,No_K
         IF(ttemp.le.TTK1(j1,i1))EXIT
      ENDDO
!
!     Condu=TTK2(j1-1,i1)+(TTK2(j1,i1)-TTK2(j1-1,i1))/(TTK1(j1,i1)-TTK1(j1-1,i1))*(ttemp-TTK1(j1-1,i1))
      Condu=TTK2(j1-1,i1)+DTK(j1-1,i1)*(ttemp-TTK1(j1-1,i1))
!
      RETURN
      END SUBROUTINE mat_prop
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
!
!     SUBROUTINE mat_prop(NoMaterial,temp,CpVol,Condu,iOKr,iOKk)
      SUBROUTINE mat_prop_2d(nmat_1d,ig_hts_1d,t_hts_1d,rcp,cond)
!
      USE Zrv_hts_1d,     ONLY: ncell_hts_1d,nr_1d,ng_hts_1d
!     Calculates the material properties (k & Cp)
!     Input: NoMaterial,temp (K)
!     Output: CpVol(J/m3.K),Condu(W/m.K),iOKr,iOKk
!
!     iOK=-1: temp. is lower   than the minimum temp.
!     iOK=+1: temp. is greater than the maximum temp.
!
      IMPLICIT NONE
!     input
      INTEGER :: nmat_1d(ng_hts_1d)
      INTEGER :: ig_hts_1d(ncell_hts_1d)
      REAL(8) :: t_hts_1d(ncell_hts_1d,nr_1d)
!     output
      REAL(8) :: rcp(ncell_hts_1d,nr_1d),cond(ncell_hts_1d,nr_1d)
!     local variables        
!   
      INTEGER  i1,j1,i
      INTEGER  k,m
!      INTEGER (4)MaxData(2,5)
      INTEGER (4)MaxData(2,8)
      INTEGER (4)NoMaterial,No_R,No_K,iOKr,iOKk   
      REAL(8) temp,CpVol,Condu,ttemp
      LOGICAL, SAVE:: initial
      DATA initial/.true./
!
!     Array for material tables
!     TR: Temp (K) vs. Volumetric heat capacity (J/m3.K)
!     TK: Temp (K) vs. Thermal conductivity (W/m.K)
!
!     1. UO2
!     2. Zircaloy
!     3. Inconel
!     4. Stainless steel
!     5. Carbon steel
!     6. BN
!     7. Monel400: Nichrome
!     8. Inconel 600
!
      REAL(8) TK1(2,32),TR1(2,26),TK2(2,16),TR2(2,15),TK3(2,7),TR3(2,7)
      REAL(8) TK4(2,8),TR4(2,8),TK5(2,8),TR5(2,8)
      REAL(8) TK6(2,10),TR6(2,10)
      REAL(8) TK7(2,10),TR7(2,10)
      REAL(8) TK8(2,10),TR8(2,10)
!
      REAL(8),SAVE:: TTK1(32,8),TTK2(32,8),DTK(31,8)
      REAL(8),SAVE:: TTR1(32,8),TTR2(32,8),DTR(31,8)
!
!      DATA MaxData/32,26,16,15,7,7,8,8,8,8/
      DATA MaxData/32,26,16,15,7,7,8,8,8,8,10,10,10,10,10,10/      
!
!
!-----------------------------------------------------------------------
!uo-2 (01)
!
!     temp[k] thermal conductivity (W/m.K)
!32
      DATA TK1/273.15d0,7.293d0,373.15d0,7.293d0,473.15d0,6.697d0,          &
      573.15d0,5.815d0,673.14d0,5.136d0,773.15d0,4.603d0,873.15d0,          &
      4.172d0,973.15d0,3.819d0,1073.15d0,3.527d0,1173.15d0,3.281d0,         &
      1273.15d0,3.075d0,1373.15d0,2.903d0,1473.15d0,2.759d0,1573.15d0,      &
      2.641d0,1673.15d0,2.546d0,1773.15d0,2.474d0,1873.15d0,2.423d0,        &
      1973.15d0,2.392d0,2073.15d0,2.382d0,2173.15d0,2.392d0,2273.15d0,      &
      2.422d0,2373.15d0,2.473d0,2473.15d0,2.546d0,2573.15d0,2.643d0,        &
      2673.15d0,2.762d0,2773.15d0,2.907d0,2873.15d0,3.078d0,2973.15d0,      &
      3.471d0,3073.15d0,3.504d0,3173.15d0,3.762d0,3573.15d0,5.133d0,        &
      4873.15d0,14.70d0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!26
      DATA TR1/273.1d0,2.427d6,400.0d0,2.754d6,500.0d0,2.927d6,600.0d0,     &
      3.043d6,700.0d0,3.139d6,800.0d0,3.178d6,900.0d0,3.236d6,              &
      1000.0d0,3.274d6,1100.0d0,3.313d6,1200.0d0,3.351d6,1300.0d0,3.378d6,  &
      1400.0d0,3.428d6,1500.0d0,3.459d6,1600.0d0,3.502d6,                   &
      1700.0d0,3.582d6,1800.0d0,3.660d6,1900.0d0,3.775d6,2000.0d0,3.992d6,  &
      2100.0d0,4.169d6,2200.0d0,4.366d6,2300.0d0,4.622d6,                   &
      2400.0d0,4.897d6,2500.0d0,5.212d6,2600.0d0,5.585d6,3000.0d0,7.395d6,  &
        4873.1d0,16.00d6/
!
!-----------------------------------------------------------------------
!zircaloy (02)
!
!     temp[k] thermal conductivity (W/m.K)
!16
      DATA TK2/273.15d0,13.6d0,373.15d0,14.1d0,473.15d0,14.8d0,573.15d0,    &
      15.8d0,673.15d0,16.9d0,773.15d0,18.1d0,873.15d0,19.5d0,               &
      973.15d0,21.1d0,1073.15d0,22.8d0,1173.15d0,24.6d0,1273.15d0,26.8d0,   &
      1373.15d0,29.2d0,1473.15d0,31.7d0,1573.15d0,34.4d0,                   &
      1673.15d0,37.3d0,1773.15d0,40.4d0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!15
      DATA TR2/273.15d0,1.881d6,573.15d0,2.079d6,773.15d0,2.211d6,903.15d0, &
      2.290d6,923.15d0,2.376d6,1083.15d0,2.376d6,1103.15d0,                 &
      3.630d6,1123.15d0,4.455d6,1143.15d0,4.950d6,1163.15d0,5.115d6,        &
      1183.15d0,4.950d6,1203.15d0,4.455d6,1213.15d0,3.360d6,1243.15d0,      &
      2.376d6,2073.15d0,2.376d6/
!     
!-----------------------------------------------------------------------
!inconel (03)
!
!     temp (k) thermal conductivity (W/m.K)
!7
      DATA TK3/0.2731500d3,0.1318846d2,0.2942611d3,0.1485516d2,0.3664833d3, &
      0.1572060d2,0.4775944d3,0.1745086d2,0.5887055d3,                      &
      0.1918174d2,0.6998167d3,0.2091262d2,0.8387055d3,0.2516817d2/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!
      DATA TR3/0.2731500d3,0.3746931d7,0.2942611d3,0.3746931d7,0.3664833d3, &
      0.3923651d7,0.4775944d3,0.4100438d7,0.5887055d3,                      &
      0.4276822d7,0.6998167d3,0.4453877d7,0.8387055d3,0.4553877d7/
!
!-----------------------------------------------------------------------
!stainless steel (04)
!
!     temp (k)   thermal conductivity (W/m.K)
!8
      DATA TK4/0.2942611d3,0.1488507d2,0.3664833d3,0.1609382d2,0.4775944d3, &
      0.1800040d2,0.5887055d3,0.1955807d2,0.6998167d3,0.2111574d2,          &
      0.8109278d3,0.2284786d2,0.9220389d3,0.2423107d2,0.1088706d4,          &
      0.2648035d2/
!
!     temp (k)   volumdtric heat capacity (J/m3.K)
!
      DATA TR4/0.2942611d3,0.3819698d7,0.3664833d3,0.3998161d7,0.4775944d3, &
      0.4227193d7,0.5887055d3,0.4355491d7,0.6998167d3,0.4446768d7,          &
      0.8109278d3,0.4563263d7,0.9220389d3,0.4625232d7,0.1088706d4,          &
      0.4750512d7/
!
!-----------------------------------------------------------------------
!carbon steel (05)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK5/0.2942611d3,0.3772982d2,0.3664833d3,0.3876847d2,0.4775944d3, &
      0.3859277d2,0.5887055d3,0.3720956d2,0.6998167d3,0.3530920d2,          &
      0.8109278d3,0.3322816d2,0.9220389d3,0.3080443d2,0.1088706d4,          &
      0.2596133d2/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR5/0.2942611d3,0.3480744d7,0.3664833d3,0.3765106d7,0.4775944d3, &
      0.4108486d7,0.5887055d3,0.4436440d7,0.6998167d3,0.4174881d7,          &
      0.8109278d3,0.5277385d7,0.9220389d3,0.6282777d7,0.1088706d4,          &
      0.5322387d7/
!
!-----------------------------------------------------------------------
! BN (06)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK6/0.293e+03,   0.0720e+02, &
               0.373e+03,   0.0749e+02, &
               0.473e+03,   0.0703e+02, &
               0.573e+03,   0.0574e+02, &
               0.673e+03,   0.0515e+02, &
               0.773e+03,   0.0469e+02, &
               0.873e+03,   0.0469e+02, &
               0.973e+03,   0.0452e+02, &
               0.1073e+04,  0.0417e+02, &
               0.1173e+04,  0.0366e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR6/0.293e+03,   1.434e+06, &  
               0.373e+03,   1.959e+06, &  
               0.473e+03,   2.335e+06, &  
               0.573e+03,   2.610e+06, &  
               0.673e+03,   2.805e+06, &  
               0.773e+03,   2.975e+06, &  
               0.873e+03,   3.106e+06, &  
               0.973e+03,   3.217e+06, &  
               0.1073e+04,  3.279e+06, &  
               0.1173e+04,  3.325e+06/
!
!-----------------------------------------------------------------------
! Monel: Nichrome (07)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK7/0.293e+03,   0.134e+02, & 
               0.373e+03,   0.143e+02, & 
               0.473e+03,   0.155e+02, & 
               0.573e+03,   0.176e+02, & 
               0.673e+03,   0.188e+02, & 
               0.773e+03,   0.208e+02, & 
               0.873e+03,   0.242e+02, & 
               0.973e+03,   0.276e+02, & 
               0.1073e+04,  0.263e+02, & 
               0.1173e+04,  0.288e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR7/0.293e+03,   3.687e+06, & 
               0.373e+03,   3.809e+06, & 
               0.473e+03,   3.933e+06, & 
               0.573e+03,   4.023e+06, & 
               0.673e+03,   4.135e+06, & 
               0.773e+03,   4.285e+06, & 
               0.873e+03,   4.678e+06, & 
               0.973e+03,   4.746e+06, & 
               0.1073e+04,  4.790e+06, & 
               0.1173e+04,  4.821e+06/               
!
!-----------------------------------------------------------------------
! Inconel 600 (08)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK8/0.293e+03,   0.143e+02, &  
               0.373e+03,   0.158e+02, &  
               0.473e+03,   0.178e+02, &  
               0.573e+03,   0.190e+02, &  
               0.673e+03,   0.200e+02, &  
               0.773e+03,   0.211e+02, &  
               0.873e+03,   0.261e+02, &  
               0.973e+03,   0.278e+02, &  
               0.1073e+04,  0.268e+02, &  
               0.1173e+04,  0.281e+02/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR8/0.293e+03,   3.718e+06, &   
               0.373e+03,   3.917e+06, &   
               0.473e+03,   4.076e+06, &   
               0.573e+03,   4.190e+06, &   
               0.673e+03,   4.247e+06, &   
               0.773e+03,   4.054e+06, &   
               0.873e+03,   4.759e+06, &   
               0.973e+03,   4.800e+06, &   
               0.1073e+04,  4.846e+06, &   
               0.1173e+04,  4.846e+06/   
!
      IF(initial)THEN
         initial=.false.
         DO i1=1,8
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)
            IF(i1.eq.1) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR1(1,j1)
                  TTR2(j1,i1)=TR1(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK1(1,j1)
                  TTK2(j1,i1)=TK1(2,j1)
               ENDDO
            ELSEIF(i1.eq.2) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR2(1,j1)
                  TTR2(j1,i1)=TR2(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK2(1,j1)
                  TTK2(j1,i1)=TK2(2,j1)
               ENDDO
            ELSEIF(i1.eq.3) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR3(1,j1)
                  TTR2(j1,i1)=TR3(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK3(1,j1)
                  TTK2(j1,i1)=TK3(2,j1)
               ENDDO
            ELSEIF(i1.eq.4) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR4(1,j1)
                  TTR2(j1,i1)=TR4(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK4(1,j1)
                  TTK2(j1,i1)=TK4(2,j1)
               ENDDO
            ELSEIF(i1.eq.5) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR5(1,j1)
                  TTR2(j1,i1)=TR5(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK5(1,j1)
                  TTK2(j1,i1)=TK5(2,j1)
               ENDDO
            ELSEIF(i1.eq.6) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR6(1,j1)
                  TTR2(j1,i1)=TR6(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK6(1,j1)
                  TTK2(j1,i1)=TK6(2,j1)
               ENDDO
            ELSEIF(i1.eq.7) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR7(1,j1)
                  TTR2(j1,i1)=TR7(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK7(1,j1)
                  TTK2(j1,i1)=TK7(2,j1)
               ENDDO
            ELSEIF(i1.eq.8) THEN
               DO j1=1,No_R
                  TTR1(j1,i1)=TR8(1,j1)
                  TTR2(j1,i1)=TR8(2,j1)
               ENDDO 
               DO j1=1,No_K
                  TTK1(j1,i1)=TK8(1,j1)
                  TTK2(j1,i1)=TK8(2,j1)
               ENDDO
            ENDIF
         ENDDO
!
         DO i1=1,8
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)
            DO i=2,No_R
               DTR(i-1,i1)=(TTR2(i,i1)-TTR2(i-1,i1))/(TTR1(i,i1)-TTR1(i-1,i1))
            ENDDO 
            DO i=2,No_K
               DTK(i-1,i1)=(TTK2(i,i1)-TTK2(i-1,i1))/(TTK1(i,i1)-TTK1(i-1,i1))
            ENDDO 
         ENDDO 
!
      ENDIF
!-----------------------------------------------------------------------
!
      DO i=1,nr_1d-1
         DO k=1,ncell_hts_1d
            m=ig_hts_1d(k)
            temp=0.5d0*(t_hts_1d(k,i)+t_hts_1d(k,i+1))
         NoMaterial=nmat_1d(m)
         i1=NoMaterial
         No_K=MaxData(1,i1)
         No_R=MaxData(2,i1)
!
!.....Volumetric heat capacity
!
      ttemp=temp
      iOKr=0
      IF(ttemp.le.TTR1(1,i1))THEN
         ttemp=TTR1(1,i1)
         iOKr=-1
      ELSEIF(ttemp.ge.TTR1(No_R,i1))THEN
         ttemp=TTR1(No_R,i1)
         iOKr=+1
      ENDIF
!
      DO j1=2,No_R
         IF(ttemp.le.TTR1(j1,i1))EXIT
      ENDDO
!
!     CpVol=TTR2(j1-1,i1)+(TTR2(j1,i1)-TTR2(j1-1,i1))/(TTR1(j1,i1)-TTR1(j1-1,i1))*(ttemp-TTR1(j1-1,i1))
      CpVol=TTR2(j1-1,i1)+DTR(j1-1,i1)*(ttemp-TTR1(j1-1,i1))
      rcp(k,i)=CpVol
!
!.....Thermal conductivity
!
      ttemp=temp
      iOKk=0
      IF(ttemp.le.TTK1(1,i1))THEN
         ttemp=TTK1(1,i1)
         iOKk=-1
      ELSEIF(ttemp.ge.TTK1(No_K,i1))THEN
         ttemp=TTK1(No_K,i1)
         iOKk=+1
      ENDIF
!
      DO j1=2,No_K
         IF(ttemp.le.TTK1(j1,i1))EXIT
      ENDDO
!
!     Condu=TTK2(j1-1,i1)+(TTK2(j1,i1)-TTK2(j1-1,i1))/(TTK1(j1,i1)-TTK1(j1-1,i1))*(ttemp-TTK1(j1-1,i1))
      Condu=TTK2(j1-1,i1)+DTK(j1-1,i1)*(ttemp-TTK1(j1-1,i1))
      cond(k,i)=Condu
!
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE mat_prop_2d
