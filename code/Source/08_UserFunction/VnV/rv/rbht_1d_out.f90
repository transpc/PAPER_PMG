!!
      SUBROUTINE rbht_1D_out
!
      USE VOL_DATA        , ONLY: cell            
      USE Zmpi            , ONLY: jperm
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank  
      USE Zconst1         , ONLY: vv_prob
      USE Zqvol           , ONLY: qporous_liq,qporous_gas
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fluid_core_all,num_ch,cupid_cell_channel
      USE Ztimecon        , ONLY: time
      USE Zwall_HTC       , ONLY: twall_rv
      USE Zrv_hts_2d      , ONLY: nz0_2d      
!      
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,m,na,nrv1,nrv2
      LOGICAL,SAVE :: initial=.true.
      REAL(8),SAVE :: print_time
      REAL(8) :: twall_rv_avg1, twall_rv_avg2, twall_rv_avg3, twall_rv_avg4
!.....Local arrays
      INTEGER :: cupid_cell_channel_gl(ncell_fluid_core)
!.....Local allocatable arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: cupid_cell_channel_tmp
      REAL(8),DIMENSION(:),ALLOCATABLE :: dat1,dat2,dat3,dat4,dat5
      REAL(8),DIMENSION(:),ALLOCATABLE :: dat_rv1
!
      na=ncell_fluid_all
      nrv1=ncell_fluid_core
      nrv2=ncell_fluid_core_all
!
      IF (vv_prob.eq.'rbht_1196fine'.or.vv_prob.eq.'rbht_1196fine2') THEN
!       
         IF(initial) THEN
            initial=.false.
            print_time=0.0d0
            IF(myrank.eq.0)OPEN(333, file='rbht_tw.dat')
            IF(myrank.eq.0)OPEN(334, file='rbht_ag.dat')
            IF(myrank.eq.0)OPEN(340, file='rbht_qporous_liq.dat')
            IF(myrank.eq.0)OPEN(341, file='rbht_qporous_gas.dat')
            IF(myrank.eq.0)OPEN(343, file='rbht_tg.dat')
            IF(myrank.eq.0)OPEN(350, file='VD12_tw_ref.dat')
            IF(myrank.eq.0)OPEN(351, file='rbht_selected_tg.dat')
            IF(myrank.eq.0)OPEN(361, file='result_calculated1.dat')
            IF(myrank.eq.0)OPEN(362, file='result_calculated2.dat')
            IF(myrank.eq.0)OPEN(363, file='result_calculated3.dat')
            IF(myrank.eq.0)OPEN(364, file='result_calculated4.dat')

         ENDIF
!         
         IF(time.ge.print_time)THEN
            print_time=print_time+1.0d0
            IF(myrank.eq.0) THEN
               ALLOCATE(dat_rv1(nrv2))
               ALLOCATE(cupid_cell_channel_tmp(nrv2))
               ALLOCATE(dat1(na),dat2(na),dat3(na),dat4(na),dat5(na))
            ELSE
               ALLOCATE(dat_rv1(1))
               ALLOCATE(cupid_cell_channel_tmp(1))
               ALLOCATE(dat1(1),dat2(1),dat3(1),dat4(1),dat5(1))
            ENDIF
            DO i=1,nrv1
               cupid_cell_channel_gl(i)=jperm(cupid_cell_channel(i))
            ENDDO            
            CALL gatherv_r(twall_rv(1,1)        ,nrv1,dat_rv1               ,nrv2,2)
            CALL gatherv_i(cupid_cell_channel_gl,nrv1,cupid_cell_channel_tmp,nrv2,2)
!
            CALL gatherv_r(cell%alphag,ncell_fluid,dat1,na,0)
            CALL gatherv_r(cell%tg    ,ncell_fluid,dat2,na,0)
            CALL gatherv_r(cell%tl    ,ncell_fluid,dat3,na,0)
            CALL gatherv_r(qporous_liq,ncell_fluid,dat4,na,0)
            CALL gatherv_r(qporous_gas,ncell_fluid,dat5,na,0)
!
            IF(myrank.eq.0) THEN
               WRITE(333,1100) time,(dat_rv1(m),m=1,49)
               WRITE(334,1100) time,(dat1(cupid_cell_channel_tmp(m)),m=1,49)
               WRITE(343,1100) time,(dat2(cupid_cell_channel_tmp(m)),m=1,49)
               WRITE(340,1100) time,(dat4(cupid_cell_channel_tmp(m)),m=1,49)
               WRITE(341,1100) time,(dat5(cupid_cell_channel_tmp(m)),m=1,49)
               WRITE(350,1100) time,dat_rv1(23),dat_rv1(30),dat_rv1(37),dat_rv1(45)
               WRITE(361,1100) time,dat_rv1(23)
               WRITE(362,1100) time,dat_rv1(30)
               WRITE(363,1100) time,dat_rv1(37)
               WRITE(364,1100) time,dat_rv1(45)
               WRITE(351,1100) time,dat3(23),dat3(30),dat3(37),dat3(45)
            ENDIF
            DEALLOCATE(dat_rv1)
            DEALLOCATE(dat1,dat2,dat3,dat4,dat5)
            DEALLOCATE(cupid_cell_channel_tmp)            
         ENDIF   
!         
1010     FORMAT(1x,1e20.10,100(i2,2x))  
1100     FORMAT(1x,1e20.10,100(e20.10,2x))   
!
      ELSEIF (vv_prob.eq.'rbht1196_1d'.or.vv_prob.eq.'rbht1196_3d') THEN
!       
         IF(initial) THEN
            initial=.false.
            print_time=0.0d0
            IF(myrank.eq.0) THEN
               OPEN(350, file='VD12_tw_ref.dat')
               OPEN(361, file='result_calculated1.dat')            
               OPEN(362, file='result_calculated2.dat')            
               OPEN(363, file='result_calculated3.dat')            
               OPEN(364, file='result_calculated4.dat') 
            ENDIF   
         ENDIF
!          
         IF(time.ge.print_time)THEN
            print_time=print_time+1.0d0
            IF(myrank.eq.0) THEN
               ALLOCATE(dat_rv1(nrv2))
            ELSE
               ALLOCATE(dat_rv1(1))
            ENDIF            
            CALL gatherv_r(twall_rv(1,1),nrv1,dat_rv1,nrv2,2)
!
            IF(myrank.eq.0) THEN
               twall_rv_avg1=0.d0
               twall_rv_avg2=0.d0
               twall_rv_avg3=0.d0
               twall_rv_avg4=0.d0
               DO i=1,num_ch
                  DO j=1,nz0_2d 
                     m=j+nz0_2d*(i-1)        
                     IF(j.eq.23)twall_rv_avg1=twall_rv_avg1+dat_rv1(m)    
                     IF(j.eq.30)twall_rv_avg2=twall_rv_avg2+dat_rv1(m)
                     IF(j.eq.37)twall_rv_avg3=twall_rv_avg3+dat_rv1(m)
                     IF(j.eq.45)twall_rv_avg4=twall_rv_avg4+dat_rv1(m)
                  ENDDO
               ENDDO            
               twall_rv_avg1=twall_rv_avg1/dble(num_ch)
               twall_rv_avg2=twall_rv_avg2/dble(num_ch)
               twall_rv_avg3=twall_rv_avg3/dble(num_ch)
               twall_rv_avg4=twall_rv_avg4/dble(num_ch)
               WRITE(350,1011) time,twall_rv_avg1,twall_rv_avg2,twall_rv_avg3,twall_rv_avg4  
               WRITE(361,1011) time,twall_rv_avg1
               WRITE(362,1011) time,twall_rv_avg2
               WRITE(363,1011) time,twall_rv_avg3
               WRITE(364,1011) time,twall_rv_avg4
            ENDIF   
            DEALLOCATE(dat_rv1)             
         ENDIF              
!         
      ENDIF   
!         
1011  FORMAT(1e14.4,5x,50(f10.4,3x))          
!
      END SUBROUTINE rbht_1D_out   
