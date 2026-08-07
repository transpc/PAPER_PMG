!
!     SUBROUTINE rv_mat_prop_2d(NoMaterial,temp,CpVol,Condu,iOKr,iOKk)
      SUBROUTINE rv_mat_prop_2d(nmat_2d,t_fuel,rcp,cond)
!
      USE Zcore           , ONLY: myrank
      USE Zrv_mpi         , ONLY: ncell_fuel_rod_p
      USE Zrv_ncell       , ONLY: ncell_fuel_rod,nrod_fuel_rod
      USE Zrv_hts_2d      , ONLY: nr_2d,nrod_2d
      use unitManager     , ONLY: createUnit
      USE Zio_unit        , ONLY: unit_log
!
!     Calculates the material properties (k & Cp)
!     Input: NoMaterial,temp (K)
!     Output: CpVol(J/m3.K),Condu(W/m.K),iOKr,iOKk
!
!     iOK=-1: temp. is lower   than the minimum temp.
!     iOK=+1: temp. is greater than the maximum temp.
!
      IMPLICIT NONE
!     input
      INTEGER :: nmat_2d(nrod_2d,nr_2d)
      REAL(8) :: t_fuel(ncell_fuel_rod_p,nr_2d)
!     output
      REAL(8) :: rcp(ncell_fuel_rod,nr_2d),cond(ncell_fuel_rod_p,nr_2d)
!     local variables
      INTEGER  i,i1,j,k,wunit
      INTEGER  err,No_M
      INTEGER(4), SAVE:: MaxData(2,5)
      INTEGER(4) NoMaterial
      INTEGER No_R,No_K,No_rod
      REAL(8) temp
      REAL(8) tempr(ncell_fuel_rod)
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
!
!      REAL(8),SAVE:: TK1(2,32),TR1(2,26),TK2(2,16),TR2(2,15),TK3(2,7),TR3(2,7)
!      REAL(8),SAVE:: TK4(2,8),TR4(2,8),TK5(2,8),TR5(2,8)
      REAL(8),SAVE:: TK1(2,32,5),TR1(2,32,5)
!
      REAL(8),SAVE :: TTK1(32,5),TTK2(32,5),DTK(31,5)
      REAL(8),SAVE :: TTR1(32,5),TTR2(32,5),DTR(31,5)
!
      DATA MaxData/32,26,16,15,7,7,8,8,8,8/
!
!
!-----------------------------------------------------------------------
!uo-2 (01)
!
!     temp[k] thermal conductivity (W/m.K)
!32
      DATA TK1(:,:,1)/300.d0,8.28385d0,400.d0,7.08576d0,500.d0,6.08653d0,                   &
      600.d0,5.31597d0,700.d0,4.72135d0,800.d0,4.25493d0,900.d0,                            &
      3.8823d0,1000.d0,3.57962d0,1100.d0,3.33042d0,1200.d0,3.12335d0,                       &
      1300.d0,2.95055d0,1364.d0, 2.8554d0, 1400.d0,2.79492d0,1500.d0,2.6453d0,1600.d0,      &
      2.52066d0,1700.d0,2.41917d0,1800.d0,2.3397d0,1834.d0, 2.33823d0, 1900.d0,2.32156d0,   &
      2000.d0,2.30826d0,2100.d0,2.31025d0,2200.d0,2.32839d0,2300.d0,                        &
      2.36339d0,2400.d0,2.46197d0,2500.d0,2.57593d0,2600.d0,2.70575d0,                      &
      2700.d0,2.85174d0,2800.d0,3.01411d0,2900.d0,3.19296d0,3000.d0,                        &
      3.3883d0,3100.d0,3.60004d0,3200.d0,3.82805d0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!26
      DATA TR1(:,:,1)/273.1d0,2.427d6,400.0d0,2.754d6,500.0d0,2.927d6,600.0d0,     &
      3.043d6,700.0d0,3.139d6,800.0d0,3.178d6,900.0d0,3.236d6,                     &
      1000.0d0,3.274d6,1100.0d0,3.313d6,1200.0d0,3.351d6,1300.0d0,3.378d6,         &
      1400.0d0,3.428d6,1500.0d0,3.459d6,1600.0d0,3.502d6,                          &
      1700.0d0,3.582d6,1800.0d0,3.660d6,1900.0d0,3.775d6,2000.0d0,3.992d6,         &
      2100.0d0,4.169d6,2200.0d0,4.366d6,2300.0d0,4.622d6,                          &
      2400.0d0,4.897d6,2500.0d0,5.212d6,2600.0d0,5.585d6,3000.0d0,7.395d6,         &
        4873.1d0,16.00d6, &
        0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0/
!
!-----------------------------------------------------------------------
!zircaloy (02)
!
!     temp[k] thermal conductivity (W/m.K)
!16
      DATA TK1(:,:,2)/273.15d0,13.6d0,373.15d0,14.1d0,473.15d0,14.8d0,573.15d0,    &
      15.8d0,673.15d0,16.9d0,773.15d0,18.1d0,873.15d0,19.5d0,                      &
      973.15d0,21.1d0,1073.15d0,22.8d0,1173.15d0,24.6d0,1273.15d0,26.8d0,          &
      1373.15d0,29.2d0,1473.15d0,31.7d0,1573.15d0,34.4d0,                          &
      1673.15d0,37.3d0,1773.15d0,40.4d0, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0/
!
!     temp[k] volumdtric heat capacity (J/m3.K)
!15
      DATA TR1(:,:,2)/273.15d0,1.881d6,573.15d0,2.079d6,773.15d0,2.211d6,903.15d0, &
      2.290d6,923.15d0,2.376d6,1083.15d0,2.376d6,1103.15d0,                        &
      3.630d6,1123.15d0,4.455d6,1143.15d0,4.950d6,1163.15d0,5.115d6,               &
      1183.15d0,4.950d6,1203.15d0,4.455d6,1213.15d0,3.360d6,1243.15d0,             &
      2.376d6,2073.15d0,2.376d6, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, &
      0.0,0.0 /
!     
!-----------------------------------------------------------------------
!inconel (03)
!
!     temp (k) thermal conductivity (W/m.K)
!7
      DATA TK1(:,:,3)/0.2731500d3,0.1318846d2,0.2942611d3,0.1485516d2,0.3664833d3, &
      0.1572060d2,0.4775944d3,0.1745086d2,0.5887055d3,                             &
      0.1918174d2,0.6998167d3,0.2091262d2,0.8387055d3,0.2516817d2,                 &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0/
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!
      DATA TR1(:,:,3)/0.2731500d3,0.3746931d7,0.2942611d3,0.3746931d7,0.3664833d3, &
      0.3923651d7,0.4775944d3,0.4100438d7,0.5887055d3,                             &
      0.4276822d7,0.6998167d3,0.4453877d7,0.8387055d3,0.4553877d7,&
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0/
!
!-----------------------------------------------------------------------
!stainless steel (04)
!
!     temp (k)   thermal conductivity (W/m.K)
!8
      DATA TK1(:,:,4)/0.2942611d3,0.1488507d2,0.3664833d3,0.1609382d2,0.4775944d3, &
      0.1800040d2,0.5887055d3,0.1955807d2,0.6998167d3,0.2111574d2,                 &
      0.8109278d3,0.2284786d2,0.9220389d3,0.2423107d2,0.1088706d4,                 &
      0.2648035d2,                                                                 &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0 /
!      
!     temp (k)   volumdtric heat capacity (J/m3.K)
!
      DATA TR1(:,:,4)/0.2942611d3,0.3819698d7,0.3664833d3,0.3998161d7,0.4775944d3, &
      0.4227193d7,0.5887055d3,0.4355491d7,0.6998167d3,0.4446768d7,                 &
      0.8109278d3,0.4563263d7,0.9220389d3,0.4625232d7,0.1088706d4,                 &
      0.4750512d7,&
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0 / 
!
!-----------------------------------------------------------------------
!carbon steel (05)
!
!     temp (k) thermal conductivity (W/m.K)
!
      DATA TK1(:,:,5)/0.2942611d3,0.3772982d2,0.3664833d3,0.3876847d2,0.4775944d3, &
      0.3859277d2,0.5887055d3,0.3720956d2,0.6998167d3,0.3530920d2,                 &
      0.8109278d3,0.3322816d2,0.9220389d3,0.3080443d2,0.1088706d4,                 &
      0.2596133d2,                                                                 &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0 /
!
!     temp (k) volumdtric heat capacity (J/m3.K)
!8
      DATA TR1(:,:,5)/0.2942611d3,0.3480744d7,0.3664833d3,0.3765106d7,0.4775944d3, &
      0.4108486d7,0.5887055d3,0.4436440d7,0.6998167d3,0.4174881d7,                 &
      0.8109278d3,0.5277385d7,0.9220389d3,0.6282777d7,0.1088706d4,                 &
      0.5322387d7,                                                                 &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0,      &
      0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0, 0.0,0.0 /
!
!.....Read external table
!
      IF(initial)THEN
         initial=.false.
         ! OPEN(5,file='ht_property_12d.in',status='old',iostat=err)
         wunit=createUnit("ht_property_12d(rv_mat_prop_2d)")
         wunit=5
         OPEN(wunit,file='ht_property_12d.in',status='old',iostat=err)
         IF(err.eq.0)THEN
            READ(wunit,*,iostat=err)No_M
            IF(err.ne.0.or.No_M.le.0)THEN
               IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to 0 of No_M.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to 0 of No_M.'
            ELSEIF(No_M.gt.0)THEN
               IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Use external table of ht_property_12d.in.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'--Use external table of ht_property_12d.in.'
               DO i=1,5
                  DO j=1,32
                     TK1(1,j,i)=0.0d0       
                     TK1(2,j,i)=0.0d0       
                     TR1(1,j,i)=0.0d0       
                     TR1(2,j,i)=0.0d0       
                  ENDDO
               ENDDO
               DO i=1,5
                  MaxData(1,i)=0  
                  MaxData(2,i)=0  
               ENDDO
               DO i=1,No_M
                  READ(wunit,*)No_K,No_R
                  MaxData(1,i)=No_K
                  MaxData(2,i)=No_R
                  DO j=1,No_K
                     READ(wunit,*)TK1(1,j,i),TK1(2,j,i)
                  ENDDO
                  DO j=1,No_R
                     READ(wunit,*)TR1(1,j,i),TR1(2,j,i)
                  ENDDO          
               ENDDO
            ELSE
               IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to unknown reason.'
               IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to unknown reason.'
            ENDIF    
         ELSE
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to lack of ht_property_12d.in.'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'--Use intrinsic table of rv_mat_prop_2d due to lack of ht_property_12d.in.'
         ENDIF  
!
         DO i1=1,5
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)
            DO j=1,No_K          
               TTK1(j,i1)=TK1(1,j,i1)         
               TTK2(j,i1)=TK1(2,j,i1)         
            ENDDO
            DO j=1,No_R          
               TTR1(j,i1)=TR1(1,j,i1)         
               TTR2(j,i1)=TR1(2,j,i1)         
            ENDDO
!
            DO j=2,No_R
               DTR(j-1,i1)=(TTR2(j,i1)-TTR2(j-1,i1))/(TTR1(j,i1)-TTR1(j-1,i1))
            ENDDO 
            DO j=2,No_K
               DTK(j-1,i1)=(TTK2(j,i1)-TTK2(j-1,i1))/(TTK1(j,i1)-TTK1(j-1,i1))
            ENDDO 
         ENDDO 
!
      ENDIF 
!
!.....Calculate conductivity and heat capacity
!      
!
      DO i=1,nr_2d-1

         DO k=1,ncell_fuel_rod
            tempr(k)=0.5d0*(t_fuel(k,i)+t_fuel(k,i+1))
         ENDDO
         DO k=1,ncell_fuel_rod
!            
            No_rod=nrod_fuel_rod(k)            
            NoMaterial=nmat_2d(No_rod,i)
            i1=NoMaterial
            No_K=MaxData(1,i1)
            No_R=MaxData(2,i1)            
!            
            temp=tempr(k)
!   
!.....Volumetric heat capacity
!
            IF(temp.le.TTR1(1,i1))THEN
               rcp(k,i)=TTR2(1,i1)
            ELSEIF(temp.ge.TTR1(No_R,i1))THEN
               rcp(k,i)=TTR2(No_R,i1)
            ELSE
!              CALL searchd1(TTR1(1,i1),No_R,temp,j)
!DIR$ NOVECTOR
               DO j=2,No_R
                  IF(temp.le.TTR1(j,i1))EXIT
               ENDDO
               rcp(k,i)=TTR2(j-1,i1)+DTR(j-1,i1)*(temp-TTR1(j-1,i1))
            ENDIF
!
!.....Thermal conductivity
!
            IF(temp.le.TTK1(1,i1))THEN
               cond(k,i)=TTK2(1,i1)
            ELSEIF(temp.ge.TTK1(No_K,i1))THEN
               cond(k,i)=TTK2(No_K,i1)
            ELSE
!              CALL searchd1(TTK1(1,i1),No_K,temp,j)
!DIR$ NOVECTOR
               DO j=2,No_K
                  IF(temp.le.TTK1(j,i1))EXIT
               ENDDO
               cond(k,i)=TTK2(j-1,i1)+DTK(j-1,i1)*(temp-TTK1(j-1,i1))
            ENDIF
         ENDDO
      ENDDO
!
      RETURN
      END SUBROUTINE rv_mat_prop_2d
!
!------------------------------------------------------------------------------
!DEC$ ATTRIBUTES INLINE :: searchd1
      SUBROUTINE searchd1(a,n,x,ip)
      IMPLICIT NONE
!
      INTEGER n,ip
      INTEGER ip1,ip2
      REAL*8  a(n),x
!
      ip1=1
      ip2=n
200   CONTINUE
      ip=(ip1+ip2)/2
      if(ip.eq.ip1) then
         ip=ip2
         goto 100
      endif
      if(x.lt.a(ip)) then
         ip2=ip
      elseif(x.gt.a(ip)) then
         ip1=ip
      else
         goto 100
      endif
      goto 200
100   CONTINUE
!
      RETURN
      END SUBROUTINE searchd1
