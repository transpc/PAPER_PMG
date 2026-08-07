!                                  
      SUBROUTINE get_mars_index  
!        
      USE Zcore,only:myrank
!
!.....print mars node information in mars_mass.dat
!
      USE TRP_BLK
      USE VOL_DAT
      USE JUN_DAT
      USE CON_VARC
!
      USE Zmars_index
!     
      IMPLICIT NONE
!     
!DEC$IF defined (MCC_DLL)
      !DEC$ ATTRIBUTES DLLIMPORT:: /c3com/,/c3com_dll/
!DEC$ENDIF      
      
      INCLUDE 'c3com.h'
!
      INTEGER(4) i,j,k,idx_,j1,j2,findopt
      INTEGER(4) loop,store_index,change_opt,mars_vol_number(1000),mars_jun_number(1000),store_number
!      
      logical,save:: initial_p
!      
      DATA initial_p /.true./
!
      ALLOCATE(mars_vol_index(1000))
      ALLOCATE(mars_jun_index(1000))
      ALLOCATE(mars_vol_jun_from(1000))
      ALLOCATE(mars_vol_jun_to(1000))
      mars_vol_index(:)=0
      mars_jun_index(:)=0
      mars_vol_jun_from(:)=0
      mars_vol_jun_to(:)=0
!
!...check mars juction and volume for print in cupid
!
      IF(v_hd%nVols(2).gt.1000.or.j_hd%nJuns(2).gt.1000)then
         WRITE(*,*)'Sizes of mars_jun_index and mars_vol_index should be smaller than 10000!'
         PAUSE
         STOP
      ENDIF
!
      DO i=1,v_hd%nVols(2)
         mars_vol_index(i)=i 
         mars_vol_number(i)=v_da(i)%volno(2)
      ENDDO
!
!.....sorting  mars_vol_index and mars_vol_number as the magnitude of mars_vol_number
      DO loop=1,100000
         change_opt=0
!         
         DO i=1,v_hd%nVols(2)-1
            IF(mars_vol_number(i).gt.mars_vol_number(i+1))then
               store_number=mars_vol_number(i)
               store_index=mars_vol_index(i)
               mars_vol_number(i)=mars_vol_number(i+1)
               mars_vol_index(i)=mars_vol_index(i+1)
               mars_vol_number(i+1)=store_number
               mars_vol_index(i+1)=store_index
               change_opt=1
            ENDIF
         ENDDO
!         
         IF(change_opt.eq.0)exit
      ENDDO
!         
      IF(loop.ge.100000)then
         WRITE(*,*)'Fail to sort mars_vol_index!'
         PAUSE
         STOP
      ENDIF
!
      DO i=1,v_hd%nVols(2)
         j=mars_vol_index(i)
         IF(myrank.eq.0)write(662,1001)j,v_da(j)%volno(2)
      ENDDO   
1001 FORMAT('vol_index,volume',3I12)
!
!......................................................................
!
      DO i=1,j_hd%nJuns(2)
         mars_jun_index(i)=i
         mars_jun_number(i)=j_da(i)%volFrom(2)
      ENDDO
!
      DO loop=1,100000
         change_opt=0
!         
            DO i=1,j_hd%nJUns(2)-1
               IF(mars_jun_number(i).gt.mars_jun_number(i+1))then
                  store_number=mars_jun_number(i)
                  store_index=mars_jun_index(i)
                  mars_jun_number(i)=mars_jun_number(i+1)
                  mars_jun_index(i)=mars_jun_index(i+1)
                  mars_jun_number(i+1)=store_number
                  mars_jun_index(i+1)=store_index
                  change_opt=1
               ENDIF
            ENDDO
!            
            IF(change_opt.eq.0)exit
      ENDDO
!      
      IF(loop.ge.100000)then
         WRITE(*,*)'Fail to sort mars_vol_index!'
         PAUSE
         STOP
      ENDIF
!      
      DO i=1,j_hd%nJUns(2)
         j=mars_jun_index(i)
         IF(myrank.eq.0)write(662,1003)j,j_da(j)%volFrom(2),j_da(j)%number(2)
      ENDDO   
!      
1003 FORMAT('jun_index,volFrom,junction',3I12)
!
!.......................................................
!
      DO i=1,j_hd%nJuns(2)
         mars_jun_index(i)=i
         mars_jun_number(i)=j_da(i)%number(2)
      ENDDO
!
      DO loop=1,100000
         change_opt=0
!         
         DO i=1,j_hd%nJUns(2)-1
            IF(mars_jun_number(i).gt.mars_jun_number(i+1))then
               store_number=mars_jun_number(i)
               store_index=mars_jun_index(i)
               mars_jun_number(i)=mars_jun_number(i+1)
               mars_jun_index(i)=mars_jun_index(i+1)
               mars_jun_number(i+1)=store_number
               mars_jun_index(i+1)=store_index
               change_opt=1
            ENDIF
         ENDDO
!         
         IF(change_opt.eq.0)exit
      ENDDO
!      
      IF(loop.ge.100000)then
         WRITE(*,*)'Fail to sort mars_vol_index!'
         PAUSE
         STOP
      ENDIF
!      
      DO i=1,j_hd%nJUns(2)
         j=mars_jun_index(i)
         IF(myrank.eq.0)write(662,1002)j,j_da(j)%volFrom(2),j_da(j)%volto(2),j_da(j)%number(2),j_da(j)%idxfrom(2),j_da(j)%idxto(2)
      ENDDO   
!      
  1002 FORMAT('jun_index,volFrom,volTo,number,idxfrom,idxto',6I12)

!
!.......................................................
!
      DO i=1,v_hd%nVols(2) !not working because inner junction don't have volfrom and volto
         DO k=1,j_hd%nJuns(2)
            j=mars_jun_index(k)
            j1=j_da(j)%idxFrom(2)
            j2=j_da(j)%idxTo(2)
            IF(j1.le.0.or.j2.le.0)CYCLE
            IF(v_da(j1)%volno(2).eq.mars_vol_number(i))mars_vol_jun_from(i)=j
            IF(v_da(j2)%volno(2).eq.mars_vol_number(i))mars_vol_jun_to(i)=j
         ENDDO
      ENDDO
!
!.....Print control variable number
!
      DO i=1,ctl_hd%nCompnt(2) !not working because inner junction don't have volfrom and volto
         IF(myrank.eq.0)write(662,1004)i,ctl_da(i)%CompntNum(2)
  1004 FORMAT('cntrl_index,cntrl_number=',6I12)
         
      ENDDO
!      
!.....Print trip
!      
      DO i=1,trp_hd%nVarTrip(2)+trp_hd%nLogicTrip(2) 
         IF(myrank.eq.0)write(662,1005)i,trp_da(i)%number(2),trp_da(i)%number(1)
    1005 FORMAT('trip_index,trip_number,trip_flag=',6I12)
      ENDDO
!      
!.....See check_trip_cntrl.f90, 'assign mass flow rates and pressures at legs'
!   
!
!.....find legs as junctions 
      DO k=1,i3n_junleg !i3n_volleg=2*i3n_junleg at cupid_mars.in
         findopt=0
         DO i=1,j_hd%nJUns(2)
            !!!!j=mars_jun_index(i)
            IF(number_junleg(k).eq.j_da(i)%number(2))THEN
               index_junleg(k)=i
               findopt=1
            ENDIF   
         ENDDO  
         IF(findopt.eq.0)THEN
            IF(myrank.eq.0)WRITE(*,"(11x,a,1i10)")'Error in finding junctions at legs in cupid_mars.in!',number_junleg(k)
            PAUSE 
            STOP
         ENDIF
      ENDDO    
!            
!.....find legs as volumes  
      DO k=1,i3n_volleg !i3n_volleg=2*i3n_volleg at cupid_mars.in
         findopt=0
         DO i=1,v_hd%nVols(2)
            !!!j=mars_vol_index(i)
            IF(number_volleg(k).eq.v_da(i)%volno(2))THEN
               index_volleg(k)=i
               findopt=1
            ENDIF   
         ENDDO  
         IF(findopt.eq.0)THEN
            WRITE(*,"(11x,a,1i3,1i10)")'Error in finding volumes at legs!',k,number_volleg(k)
            PAUSE
            STOP
         ENDIF               
      ENDDO 
!            
!.....find tmdpvol2nd  
      DO k=1,i3n_tmdpvol2nd 
         findopt=0
         DO i=1,v_hd%nVols(2)
            !!!j=mars_vol_index(i)
            IF(number_tmdpvol2nd(k).eq.v_da(i)%volno(2))THEN
               index_tmdpvol2nd(k)=i
               findopt=1
            ENDIF   
         ENDDO  
         IF(findopt.eq.0)THEN
            WRITE(*,"(11x,a,1i3,1i10)")'Error in finding tmdpvol2nd!',k,number_tmdpvol2nd(k)
            PAUSE
            STOP
         ENDIF               
      ENDDO 
!
!      IF(myrank.eq.0)WRITE(97,"(11x,a)")'##finish set_mars_cell_no !'
!
      RETURN
      END SUBROUTINE get_mars_index    
     
