!
      SUBROUTINE fs_1D_out
!
      USE VOL_DATA        , ONLY: cell
      USE Zmpi            , ONLY: jperm
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all
      USE Zcore           , ONLY: myrank  
      USE Zparam          , ONLY: ndim 
      USE Ztimecon        , ONLY: time
      USE Zconst1         , ONLY: vv_prob
      USE Zqvol           , ONLY: qporous_liq,qporous_gas
      USE Zrv_ncell       , ONLY: ncell_fluid_core,ncell_fluid_core_all,num_ch,cupid_cell_channel
      USE Zwall_HTC       , ONLY: twall_rv,hmode_rv
      USE Zrv_hts_2d      , ONLY: nz0_2d      
      USE Zcoord1         , ONLY: xloc_tmp
!      
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,na,nrv1,nrv2,ii
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: temp_an
      REAL(8),SAVE :: print_time
!
!.....Local arrays
      INTEGER,DIMENSION(:),ALLOCATABLE :: dat_rv2_ch
      INTEGER,DIMENSION(:),ALLOCATABLE :: cupid_cell_channel_gl
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: cupid_cell_channel_tmp
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: dat_rv2
      REAL(8),DIMENSION(:),ALLOCATABLE :: dat_rv1_ch,dat_rv1_rod
      REAL(8),DIMENSION(:),ALLOCATABLE :: tw4_nz,tw7_nz,tw11_nz,tw14_nz,tw17_nz,tw19_nz
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat,dat_ch
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: dat_rv1,xloc_rod
!
      na=ncell_fluid_all
      nrv1=ncell_fluid_core
      nrv2=ncell_fluid_core_all
!
      IF(initial) THEN
         print_time=0.0d0
         IF(vv_prob.eq.'fs_31203'    .or.&
            vv_prob.eq.'fs_31302'    .or.&            
            vv_prob.eq.'fs_31701'    .or.&            
            vv_prob.eq.'fs_31805'    ) THEN
            IF(myrank.eq.0) THEN
               OPEN(330, file='VD11_tw_result-24.dat') 
               OPEN(331, file='VD11_tw_result-76.dat') 
               OPEN(332, file='VD11_tw_result-120.dat') 
               OPEN(333, file='VD11_tw_ref.dat')  
               OPEN(334, file='fs_31701_vv.dat')  
               OPEN(335, file='VD11_PCT_time.dat')
            ENDIF
            temp_an=time
         ELSEIF(vv_prob.eq.'fs_31203_3D'    .or.&
                vv_prob.eq.'fs_31302_3D'    .or.&            
                vv_prob.eq.'fs_31701_3D'    .or.&            
                vv_prob.eq.'fs_31805_3D'    ) THEN
            IF(myrank.eq.0) THEN
               OPEN(330, file='VD11_tw_result-24.dat') 
               OPEN(331, file='VD11_tw_result-76.dat') 
               OPEN(332, file='VD11_tw_result-120.dat') 
               OPEN(333, file='VD11_tw_ref.dat')  
               OPEN(334, file='fs_31701_3D_vv.dat')  
               OPEN(335, file='VD11_PCT_time.dat') 
               OPEN(336, file='tw061_time.dat') 
               OPEN(337, file='tw122_time.dat') 
               OPEN(338, file='tw182_time.dat') 
               OPEN(339, file='tw244_time.dat') 
               OPEN(340, file='tw305_time.dat') 
               OPEN(341, file='tw335_time.dat') 
            ENDIF
            temp_an=time
         ELSEIF(vv_prob.eq.'fs_31504') THEN
            IF(myrank.eq.0) THEN
               OPEN(333, file='fs_31504_ref.dat')  
               OPEN(334, file='fs_31504_vv.dat') 
            ENDIF
            temp_an=time-71.0d0
         ENDIF
         IF(myrank.eq.0) THEN
            WRITE(333,"(a)")'Time Tw-061 Tw-122 Tw-182 Tw-244 Tw-305 Tw-335             &
                                  ag-061 ag-122 ag-182 ag-244 ag-305 ag-335             &
                                  mode-061 mode-122 mode-182 mode-244 mode-305 mode-335 &
                                  Tl-061 Tl-122 Tl-182 Tl-244 Tl-305 Tl-335             &
                                  Tg-061 Tg-122 Tg-182 Tg-244 Tg-305 Tg-335'
            WRITE(334,"(a)")'Time Tw-1 Tw-2 Tw-3 Tw-4 Tw-5 Tw-6 &
                                  Ql-1 Ql-2 Ql-3 Ql-4 Ql-5 Ql-6 &
                                  Qg-1 Qg-2 Qg-3 Qg-4 Qg-5 Qg-6 &
                                  Gw-1 Gw-2 Gw-3 Gw-4 Gw-5 Gw-6 &
                                  ag-1 ag-2 ag-3 ag-4 ag-5 ag-6'
         ENDIF
      ENDIF
!         
      IF(time.ge.print_time) THEN !4,7,11,14,17,19->upside down=(20-i+1)->17,14,10,7,4,2
!      
         IF(myrank.eq.0) THEN
            ALLOCATE(dat(na,5),dat_ch(nz0_2d,5))
            ALLOCATE(cupid_cell_channel_tmp(nz0_2d,num_ch))
            ALLOCATE(dat_rv1(nz0_2d,num_ch),dat_rv2(nz0_2d,num_ch))
         ELSE
            ALLOCATE(dat(1,5),dat_ch(1,5))
            ALLOCATE(cupid_cell_channel_tmp(1,1))
            ALLOCATE(dat_rv1(1,1),dat_rv2(1,1))
         ENDIF
!
         ALLOCATE(cupid_cell_channel_gl(nrv1))
         DO i=1,nrv1
            cupid_cell_channel_gl(i)=jperm(cupid_cell_channel(i))
         ENDDO
         CALL gatherv_i(cupid_cell_channel_gl,nrv1,cupid_cell_channel_tmp,nrv2,2)
         CALL gatherv_r(twall_rv(1,1)        ,nrv1,dat_rv1               ,nrv2,2)
         CALL gatherv_i(hmode_rv(1,1)        ,nrv1,dat_rv2               ,nrv2,2)
         DEALLOCATE(cupid_cell_channel_gl)
!
         CALL gatherv_r(cell%alphag,ncell_fluid,dat(1,1),na,0)
         CALL gatherv_r(cell%tg    ,ncell_fluid,dat(1,2),na,0)
         CALL gatherv_r(cell%tl    ,ncell_fluid,dat(1,3),na,0)
         CALL gatherv_r(qporous_liq,ncell_fluid,dat(1,4),na,0)
         CALL gatherv_r(qporous_gas,ncell_fluid,dat(1,5),na,0)
!
         IF(myrank.eq.0) THEN
            ALLOCATE(dat_rv1_ch(nz0_2d),dat_rv2_ch(nz0_2d)) 
            ALLOCATE(tw4_nz(num_ch),tw7_nz(num_ch),tw11_nz(num_ch),tw14_nz(num_ch),tw17_nz(num_ch),tw19_nz(num_ch)) !Twall for all channel for some height
            dat_rv1_ch(:)=0.d0
            dat_rv2_ch(:)=0
            dat_ch(:,1)=0.d0
            dat_ch(:,2)=0.d0
            dat_ch(:,3)=0.d0
            dat_ch(:,4)=0.d0
            dat_ch(:,5)=0.d0
            DO i=1,num_ch
               DO j=1,nz0_2d         
                  ii=cupid_cell_channel_tmp(j,i)
                  dat_rv1_ch(j)=dat_rv1_ch(j)+dat_rv1(j,i)
                  dat_rv2_ch(j)=dat_rv2_ch(j)+dat_rv2(j,i)
                  dat_ch(j,1)=dat_ch(j,1)+dat(ii,1)
                  dat_ch(j,2)=dat_ch(j,2)+dat(ii,2)
                  dat_ch(j,3)=dat_ch(j,3)+dat(ii,3)
                  dat_ch(j,4)=dat_ch(j,4)+dat(ii,4)
                  dat_ch(j,5)=dat_ch(j,5)+dat(ii,5)
                  IF(j.eq.4) THEN
                     tw4_nz(i)=dat_rv1(j,i)
                  ELSEIF(j.eq.7) THEN
                     tw7_nz(i)=dat_rv1(j,i)
                  ELSEIF(j.eq.11) THEN
                     tw11_nz(i)=dat_rv1(j,i)
                  ELSEIF(j.eq.14) THEN
                     tw14_nz(i)=dat_rv1(j,i)
                  ELSEIF(j.eq.17) THEN
                     tw17_nz(i)=dat_rv1(j,i)
                  ELSEIF(j.eq.19) THEN
                     tw19_nz(i)=dat_rv1(j,i)
                  ENDIF
               ENDDO
            ENDDO
            DO j=1,nz0_2d         
               dat_rv1_ch(j)=dat_rv1_ch(j)/num_ch
               dat_rv2_ch(j)=dat_rv2_ch(j)/num_ch
               dat_ch(j,1)=dat_ch(j,1)/num_ch
               dat_ch(j,2)=dat_ch(j,2)/num_ch
               dat_ch(j,3)=dat_ch(j,3)/num_ch
               dat_ch(j,4)=dat_ch(j,4)/num_ch
               dat_ch(j,5)=dat_ch(j,5)/num_ch
            ENDDO            
!...........PCT for each channel
            ALLOCATE(dat_rv1_rod(num_ch),xloc_rod(num_ch,ndim)) !PCT
            dat_rv1_rod(:)=0.0d0
            DO i=1,num_ch
               ii=cupid_cell_channel_tmp(1,i)
               xloc_rod(i,:)=xloc_tmp(ii,:)
               dat_rv1_rod(i)=dat_rv1(1,i)
               DO j=2,nz0_2d         
                  dat_rv1_rod(i)=MAX(dat_rv1_rod(i),dat_rv1(j,i))
               ENDDO
            ENDDO   
         ENDIF
!            
         IF(vv_prob.eq.'fs_31203'.or.vv_prob.eq.'fs_31203_3D'.or.&
            vv_prob.eq.'fs_31302'.or.vv_prob.eq.'fs_31302_3D'.or.&
            vv_prob.eq.'fs_31701'.or.vv_prob.eq.'fs_31701_3D'.or.&
            vv_prob.eq.'fs_31805'.or.vv_prob.eq.'fs_31805_3D') THEN
            temp_an=time
         ELSEIF(vv_prob.eq.'fs_31504') THEN
            temp_an=time-71.0d0
         ELSE
            temp_an=time            
         ENDIF
!         
         IF(vv_prob.eq.'fs_31203'.or.vv_prob.eq.'fs_31203_3D'.or.&
            vv_prob.eq.'fs_31302'.or.vv_prob.eq.'fs_31302_3D'.or.&
            vv_prob.eq.'fs_31701'.or.vv_prob.eq.'fs_31701_3D'.or.&
            vv_prob.eq.'fs_31805'.or.vv_prob.eq.'fs_31805_3D') THEN
            print_time=print_time+0.5d0
            IF(myrank.eq.0) THEN
               WRITE(330,1002) temp_an,dat_rv1_ch(4)
               WRITE(331,1002) temp_an,dat_rv1_ch(11)
               WRITE(332,1002) temp_an,dat_rv1_ch(17)
               WRITE(333,1002) temp_an,                                                                                         &
!                                       dat_rv1(4),dat_rv1(7),dat_rv1(11),dat_rv1(14),dat_rv1(17),dat_rv1(19),                  &
                                       dat_rv1_ch(4),dat_rv1_ch(7),dat_rv1_ch(11),dat_rv1_ch(14),dat_rv1_ch(17),dat_rv1_ch(19), &
                                       dat_ch(4,1),dat_ch(7,1),dat_ch(11,1),dat_ch(14,1),dat_ch(17,1),dat_ch(19,1),             &
                                       dat_rv2(4,1),dat_rv2(7,1),dat_rv2(11,1),dat_rv2(14,1),dat_rv2(17,1),dat_rv2(19,1),       &
                                       dat_ch(4,2),dat_ch(7,2),dat_ch(11,2),dat_ch(14,2),dat_ch(17,2),dat_ch(19,2),             &
                                       dat_ch(4,3),dat_ch(7,3),dat_ch(11,3),dat_ch(14,3),dat_ch(17,3),dat_ch(19,3)
               WRITE(334,1002) temp_an,                                                                               &
                                       dat_rv1(1,1),dat_rv1(2,1),dat_rv1(3,1),dat_rv1(4,1),dat_rv1(5,1),dat_rv1(6,1), &
                                       dat_ch(1,1),dat_ch(2,1),dat_ch(3,1),dat_ch(4,1),dat_ch(5,1),dat(6,1),          &
                                       dat_rv2(1,1),dat_rv2(2,1),dat_rv2(3,1),dat_rv2(4,1),dat_rv2(5,1),dat_rv2(6,1), &
                                       dat_ch(1,4),dat_ch(2,4),dat_ch(3,4),dat_ch(4,4),dat_ch(5,4),dat_ch(6,4),       &
                                       dat_ch(1,5),dat_ch(2,5),dat_ch(3,5),dat_ch(4,5),dat_ch(5,5),dat_ch(6,5)
            ENDIF
            IF(initial) THEN
               IF(myrank.eq.0) THEN
                  WRITE(335,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(335,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(335,1003)num_ch,(xloc_rod(i,ndim),i=1,num_ch)
               ENDIF
            ENDIF
!
            IF(myrank.eq.0) THEN
               WRITE(335,1004)temp_an,(dat_rv1_rod(i),i=1,num_ch)
            ENDIF
!            
            IF(initial) THEN
               IF(myrank.eq.0) THEN
                  WRITE(336,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(336,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(337,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(337,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(338,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(338,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(339,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(339,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(340,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(340,1003)num_ch,(xloc_rod(i,2),i=1,num_ch)
                  WRITE(341,1003)num_ch,(xloc_rod(i,1),i=1,num_ch)
                  WRITE(341,1003)num_ch,(xloc_rod(i,2),i=1,num_ch) 
               ENDIF
            ENDIF               
            IF(myrank.eq.0) THEN
               WRITE(336,1004)temp_an,(tw4_nz(i),i=1,num_ch)         
               WRITE(337,1004)temp_an,(tw7_nz(i),i=1,num_ch)  
               WRITE(338,1004)temp_an,(tw11_nz(i),i=1,num_ch)  
               WRITE(339,1004)temp_an,(tw14_nz(i),i=1,num_ch)  
               WRITE(340,1004)temp_an,(tw17_nz(i),i=1,num_ch)  
               WRITE(341,1004)temp_an,(tw19_nz(i),i=1,num_ch) 
            ENDIF
         ENDIF
!
         DEALLOCATE(dat_rv1,dat_rv2)
         DEALLOCATE(cupid_cell_channel_tmp)
         IF(myrank.eq.0) THEN
            DEALLOCATE(dat_rv1_ch,dat_rv2_ch)
            DEALLOCATE(dat,dat_ch)
            DEALLOCATE(tw4_nz,tw7_nz,tw11_nz,tw14_nz,tw17_nz,tw19_nz)
            DEALLOCATE(dat_rv1_rod,xloc_rod)
         ENDIF
!         
      ENDIF
      IF(initial) initial=.false.
!
 1002 FORMAT(1x,1e20.10,2(2x,6e14.5),(2x,6i3),2(2x,6e14.5))          
 1003 FORMAT(1x,1i20,1000(2x,6e14.5))          
 1004 FORMAT(1x,1e20.10,1000(2x,6e14.5))          
!
      END SUBROUTINE fs_1D_out   
