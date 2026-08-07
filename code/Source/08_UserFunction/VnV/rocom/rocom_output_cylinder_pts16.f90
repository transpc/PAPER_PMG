!      
!.....These SUBROUTINE is very special print routine for secIFic problem
!.....Don't USE for another problem
!
      SUBROUTINE rocom_output_dc_pts16(time,prnopt,stavg,ftavg)
!
!.....Save lateral output for Origin
!   
      USE VOL_DATA                     
      USE Zparam          , ONLY: ndim,pi
      USE Zconst1         , ONLY: vv_prob
      USE Zcore           , ONLY: np,myrank
      USE Zcoord3         , ONLY: vol
      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min,idx_dc_top,idx_dc_bot,size_idx
      USE Zvector         , ONLY: vl_n
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all         
!
      IMPLICIT NONE 
!
      CHARACTER title(30) * 20
      CHARACTER filename(6)*30, sstime*16
      INTEGER prnopt     
      INTEGER i, j, level
      INTEGER Nradius,levfile      
      INTEGER time_avg_opt
      INTEGER iitime   
      INTEGER,SAVE:: sum_num    
      REAL(8) stavg,ftavg
      REAL(8) time
!     REAL(8) radius(2), h(2), hzero
      REAL(8) vol_sum(2),boron_vol_sum(2),boron_max(2),boron_vol_avg(2)
      LOGICAL,SAVE:: initial_p      
!     DATA  Nradius,hzero,h(1),h(2),radius(1),radius(2)/2, 0.0d0, -1.2195d0, -0.3105d0, 0.437d0, 0.498d0/ !rocom:0.376,1.69, -317.5, -1185.5, 0.437, 0.5
      DATA  Nradius/2/
      DATA initial_p / .true./

      REAL(8),SAVE, ALLOCATABLE:: vli_sum(:,:),cboroni_sum(:)
      REAL(8),SAVE, ALLOCATABLE:: vli_ave(:,:),cboroni_ave(:)
      REAL(8), ALLOCATABLE::vlp(:,:)
      REAL(8), ALLOCATABLE::cboron(:)
      REAL(8), ALLOCATABLE::volall(:)
!
      ALLOCATE(vlp(ncell_fluid_all,ndim))
      ALLOCATE(cboron(ncell_fluid_all))
      ALLOCATE(volall(ncell_fluid_all))
!
      IF(np.gt.1)THEN   
         CALL allgatherv_r(vl_n(1,1),vlp(1,1),ncell_fluid,ncell_fluid_all,0)
         CALL allgatherv_r(vl_n(1,2),vlp(1,2),ncell_fluid,ncell_fluid_all,0)
         IF(ndim.eq.3) CALL allgatherv_r(vl_n(1,ndim),vlp(1,ndim),ncell_fluid,ncell_fluid_all,0)
         CALL allgatherv_r(cell%cboron,cboron,ncell_fluid,ncell_fluid_all,0)
         CALL allgatherv_r(vol,volall,ncell_fluid,ncell_fluid_all,0)
      ELSE
         vlp(:,1:ndim)=vl_n(:,1:ndim)
         cboron(:)=cell%cboron(:)
         volall(:)=vol(:)       
      ENDIF   

      IF(myrank.gt.0)goto 100    
      IF(ndim.ne.3)RETURN 
      
      IF (initial_p) THEN
          initial_p = .false.
          ALLOCATE(vli_sum(ncell_fluid_all,ndim),cboroni_sum(ncell_fluid_all))
          ALLOCATE(vli_ave(ncell_fluid_all,ndim),cboroni_ave(ncell_fluid_all))
          sum_num=0
          cboroni_sum(:)=0.0d0
          vli_sum(:,:)=0.0d0                
      ENDIF
!....initialize the variable for time-averaging
      IF(initial_p.eq..false.)THEN
         IF(time.ge.stavg .and. time.le.ftavg)THEN
             time_avg_opt=1
             sum_num=sum_num+1
             DO i=1,ncell_fluid_all                
                vli_sum(i,:)=vli_sum(i,:)+vlp(i,:)
                cboroni_sum(i)=cboroni_sum(i)+cboron(i)
                vli_ave(i,:)=vli_sum(i,:)/sum_num
                cboroni_ave(i)=cboroni_sum(i)/sum_num                
             ENDDO
          ELSE
             time_avg_opt=0   
             sum_num=1         
             DO i=1,ncell_fluid_all
                vli_sum(i,:)=vlp(i,:)
                cboroni_sum(i)=cboron(i)
                vli_ave(i,:)=vli_sum(i,:)/sum_num
                cboroni_ave(i)=cboroni_sum(i)/sum_num                
            ENDDO   
         ENDIF     
      ENDIF    
!
!.....initialize the variables for space-averaging
!
      vol_sum(:)=0.0d0
      boron_vol_sum(:)=0.0d0
      boron_vol_avg(:)=0.0d0
      boron_max(:)=0.0d0
!      
!.....for the files for every PRINT step
!
      title(1)='"A:Radius",'  
      title(2)='"B:Theta",'  
      title(3)='"C:Height",'
      title(4)=''                 !'"D:Tliquid"'
      title(5)='"E:Cboron"'            
      title(6)='"F:Vliquidx",'  
      title(7)='"G:Vliquidy",'  
      title(8)='"H:Vliquidz",'  
      title(9)=''                 !'"I:Pressure"' 
      title(10)=''                !'"J:Tliquid-K"'      
      IF(prnopt)THEN
         iitime=100000+time*100
         IF(myrank.eq.0)WRITE(sstime,*)iitime
         sstime=adjustl(sstime)
         filename(1)='DC_low_Av_'//trim(sstime)//'.dat' 
         filename(2)='DC_up_Av_'//trim(sstime)//'.dat' 
         filename(3)='DC_low_Ins_'//trim(sstime)//'.dat' 
         filename(4)='DC_up_Ins_'//trim(sstime)//'.dat' 
         IF(myrank.eq.0)THEN
            OPEN(61,file=filename(1)) !
            OPEN(62,file=filename(2)) !
            OPEN(63,file=filename(3)) !
            OPEN(64,file=filename(4)) !
         ENDIF
      ENDIF         
      
!
!.....for the 2 heights of DC
!
      DO level = 1, Nradius 
         levfile = 60 + level
         IF(myrank.eq.0)THEN
            IF(prnopt)WRITE(levfile, 8003)(title(i), i=1, 10)
            IF(prnopt)WRITE(levfile+2, 8003)(title(i), i=1, 10)
         ENDIF
  8003   FORMAT('variables=',10A16)
!    
         DO j = 1, size_idx(level)
            IF(level.eq.1)i=idx_dc_bot(j)
            IF(level.eq.2)i=idx_dc_top(j)
!       
            boron_vol_sum(level)=boron_vol_sum(level)+cboron(i)*volall(i)
            IF(boron_max(level).lt.cboron(i))boron_max(level)=cboron(i) 
            vol_sum(level)=vol_sum(level)+volall(i)                       
!
            IF(prnopt)THEN                
               WRITE(levfile,5001)cboroni_ave(i),&                  !tli_ave(i),
                                   vli_ave(i,1),vli_ave(i,2),vli_ave(i,ndim)             !,pi_ave(i),tli_ave(i)+273.15d0                          
               WRITE(levfile+2,5001)cboron(i),&                     !,tlall(i)-273.15d0
                                   vlp(i,1),vlp(i,2),vlp(i,ndim)                         !,pr(i),tlall(i)
            ENDIF    
         ENDDO !i
         IF(prnopt)THEN 
            WRITE(levfile,*)'time=',time
            CLOSE(levfile)   
            WRITE(levfile+2,*)'time=',time
            CLOSE(levfile+2)   
         ENDIF         
      ENDDO !level
 5001 FORMAT(1x,10e16.5)
!
!.....averaging the tl for the 2-d plane
!        
      DO i=1, nradius
        IF(vol_sum(i).gt.0.d0)THEN
           boron_vol_avg(i)=boron_vol_sum(i)/vol_sum(i)
        ELSE
           PRINT *,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
        ENDIF           
      ENDDO    
      tl_space_avg(1)=boron_vol_avg(2)  
      tl_space_avg(2)=boron_vol_avg(1)     
      tl_space_min(1)=boron_max(2)
      tl_space_min(2)=boron_max(1)  
!
      IF(vv_prob.eq.'pts16')THEN
!         tl_space_avg(4)=cboron(29234)           !300,000 grids
!         tl_space_avg(5)=cboron(46254)
!
!         tl_space_avg(4)=cboron(3722028)          !468 NewCL_BL
!         tl_space_avg(5)=cboron(3667428)          !Lower DC
!
         tl_space_avg(4)=cboron(2319900)          !1243 NewCL_BL
         tl_space_avg(5)=cboron(2538015)          !Lower DC       
      ENDIF
!
100   CONTINUE   
   
      DEALLOCATE(vlp)
      DEALLOCATE(volall)
      DEALLOCATE(cboron)  
!    
      RETURN
      END



!      
!.....These SUBROUTINE is very special print routine for secIFic problem
!.....Don't USE for another problem
!
      SUBROUTINE rocom_output_coreinlet_pts16(time,prnopt,stavg,ftavg)
!
!.....Save lateral output for Origin
!   
      USE VOL_DATA                     
      USE Zparam          , ONLY: ndim,pi
      USE Zcore           , ONLY: np,myrank
      USE Zcoord3         , ONLY: vol
      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min,idx_core,size_idx
      USE Zvector         , ONLY: vl_n
      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all         
!
      IMPLICIT NONE 
!
      CHARACTER title(30) * 20
      CHARACTER filename(6)*20, sstime*16
      INTEGER prnopt     
      INTEGER i, j, loop, level
      INTEGER Nradius      
      INTEGER time_avg_opt
      INTEGER iitime   
      INTEGER,SAVE:: sum_num    
      REAL(8) stavg,ftavg
      REAL(8) time
!     REAL(8) radius, h1, h2, hzero
      REAL(8) yy,zz,rr!
      REAL(8) vol_sum(2),boron_vol_sum(2),boron_max(2),boron_vol_avg(2)
      LOGICAL,SAVE:: initial_p      
!     DATA  Nradius,hzero,h1,h2,radius/1, 0.0d0, -1.2515, -1.2515, 0.437/ !rocom:0.376,1.69, -317.5, -1185.5, 0.437, 0.5
      DATA  Nradius/1/
      DATA initial_p / .true./

      REAL(8),SAVE, ALLOCATABLE:: vli_sum(:,:),cboroni_sum(:)
      REAL(8),SAVE, ALLOCATABLE:: vli_ave(:,:),cboroni_ave(:)
      REAL(8), ALLOCATABLE::vlp(:,:)
      REAL(8), ALLOCATABLE::cboron(:)
      REAL(8), ALLOCATABLE::volall(:)
!
      ALLOCATE(vlp(ncell_fluid_all,ndim))
      ALLOCATE(cboron(ncell_fluid_all))
      ALLOCATE(volall(ncell_fluid_all))
!
      IF(np.gt.1)THEN   
         CALL allgatherv_r(vl_n(1,1),vlp(1,1),ncell_fluid,ncell_fluid_all,0)
         CALL allgatherv_r(vl_n(1,2),vlp(1,2),ncell_fluid,ncell_fluid_all,0)
         IF(ndim.eq.3) CALL allgatherv_r(vl_n(1,ndim),vlp(1,ndim),ncell_fluid,ncell_fluid_all,0)
!         
         CALL allgatherv_r(cell%cboron,cboron,ncell_fluid,ncell_fluid_all,0)
         CALL allgatherv_r(vol,volall,ncell_fluid,ncell_fluid_all,0)
      ELSE
         vlp(:,1:ndim)=vl_n(:,1:ndim)
         cboron(:)=cell%cboron(:)
         volall(:)=vol(:)
      ENDIF   
!
      IF(myrank.gt.0)GOTO 100     
      IF(ndim.ne.3)RETURN 
      
      IF (initial_p) THEN
          initial_p = .false.
          ALLOCATE(vli_sum(ncell_fluid_all,ndim),cboroni_sum(ncell_fluid_all))
          ALLOCATE(vli_ave(ncell_fluid_all,ndim),cboroni_ave(ncell_fluid_all))
          sum_num=0
          cboroni_sum(:)=0.0d0
          vli_sum(:,:)=0.0d0         
      ENDIF
!....initialize the variable for time-averaging
      IF(initial_p.eq..false.)THEN
         IF(time.ge.stavg .and. time.le.ftavg)THEN
             time_avg_opt=1
             sum_num=sum_num+1
             DO i=1,ncell_fluid_all
                cboroni_sum(i)=cboroni_sum(i)+cboron(i)
                vli_sum(i,:)=vli_sum(i,:)+vlp(i,:)
                cboroni_ave(i)=cboroni_sum(i)/sum_num
                vli_ave(i,:)=vli_sum(i,:)/sum_num
             ENDDO
          ELSE
             time_avg_opt=0   
             sum_num=1         
             DO i=1,ncell_fluid_all
                cboroni_sum(i)=cboron(i)
                vli_sum(i,:)=vlp(i,:)                
                cboroni_ave(i)=cboroni_sum(i)/sum_num
                vli_ave(i,:)=vli_sum(i,:)/sum_num
            ENDDO   
         ENDIF     
      ENDIF    
!
!.....initialize the variables for space-averaging
!
      vol_sum(:)=0.0d0
      boron_vol_sum(:)=0.0d0
      boron_vol_avg(:)=0.0d0
      boron_max(:)=0.0d0
!      
!.....for the files for every PRINT step
!
      title(1)='"A:Radius",'  
      title(2)='"B:Theta",'  
      title(3)='"C:Height",'
      title(4)=''             !'"D:tl",'
      title(5)='"E:cBoron"'
      title(6)='"F:Vliquidx"'            
      title(7)='"G:Vliquidy",'  
      title(8)='"H:Vliquidz",'  
      title(9)=''             !'"I:Pressure",'  
      title(10)=''            !'"J:Tliquid-K"'    
        
      IF(prnopt)THEN
         iitime=100000+time*100
         IF(myrank.eq.0)WRITE(sstime,*)iitime
         sstime=adjustl(sstime)
         filename(1)='core_Av_'//trim(sstime)//'.dat' 
         filename(2)='core_ins_'//trim(sstime)//'.dat' 
         IF(myrank.eq.0)THEN
            OPEN(65,file=filename(1)) !
            OPEN(66,file=filename(2)) !
         ENDIF
      ENDIF         
!
!.....for the 2 radius
!
      DO level = 1, Nradius 
         loop = 0
         IF(myrank.eq.0)THEN
            IF(prnopt)WRITE(65, 8003)(title(i), i=1, 10)
            IF(prnopt)WRITE(66, 8003)(title(i), i=1, 10)
         ENDIF
  8003   FORMAT('variables=',10A16)

         DO j = 1, size_idx(3)
            i=idx_core(j)

            boron_vol_sum(level)=boron_vol_sum(level)+cboron(i)*volall(i)
            IF(boron_max(level).lt.cboron(i))boron_max(level)=cboron(i) 
            vol_sum(level)=vol_sum(level)+volall(i)                                 
!
            IF(prnopt)THEN 
                WRITE(65,5001)rr,yy,zz,cboroni_ave(i),&                           !,tli_ave(i)
                               vli_ave(i,1),vli_ave(i,2),vli_ave(i,ndim)           !,pi_ave(i),tli_ave(i)+273.15d0                          
                WRITE(66,5001)rr,yy,zz,cboron(i),&                                !,tlall(i)-273.15d0
                               vlp(i,1),vlp(i,2),vlp(i,ndim)                       !,pr(i),tlall(i)
            ENDIF     
         ENDDO !i
!
         IF(prnopt)THEN 
            WRITE(65,*)'time=',time
            CLOSE(65)   
            WRITE(66,*)'time=',time
            CLOSE(66)  
         ENDIF         
      ENDDO !level
 5001 FORMAT(1x,10e16.5)
!
      DO i=1, nradius
        IF(vol_sum(i).gt.0.d0)THEN
           boron_vol_avg(i)=boron_vol_sum(i)/vol_sum(i)
        ELSE
           PRINT *,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
        ENDIF           
      ENDDO    
      tl_space_avg(3)=boron_vol_avg(1)  
      tl_space_min(3)=boron_max(1)       
! 
100   CONTINUE
!   
      DEALLOCATE(vlp)
      DEALLOCATE(volall)
      DEALLOCATE(cboron)  
!    
      RETURN
      END
!













































!!      
!!.....These SUBROUTINE is very special print routine for secIFic problem
!!.....Don't USE for another problem
!!
!      SUBROUTINE rocom_output_dc_pts16(time,prnopt,stavg,ftavg)
!!
!!.....Save lateral output for Origin
!!   
!      USE VOL_DATA                     
!      USE Zparam          , ONLY: ndim,ns,pi
!      USE Zconst1         , ONLY: vv_prob
!      USE Zcore           , ONLY: np,myrank
!      USE Zcoord1         , ONLY: xloc_tmp,xfc_tmp
!      !USE Zcoord2         , ONLY: xfc
!      USE Zcoord3         , ONLY: vol
!      USE Znum_cell       , ONLY: num_neigh,num_neigh_tmp  
!      USE Zpress          , ONLY: p               
!      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min,rad_group
!      USE Zvector         , ONLY: vl_n
!      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all         
!      USE Znode           , ONLY: num_cell_node,cell_node,xnode,nd_max   
!!
!      IMPLICIT NONE 
!!
!      CHARACTER title(30) * 20
!      CHARACTER filename(6)*30, sstime*16
!      INTEGER prnopt     
!      INTEGER i, j, ic, ix, wopt, loop, sortopt, level
!      INTEGER Nradius,levfile      
!      INTEGER time_avg_opt
!      INTEGER iitime   
!      INTEGER,SAVE:: sum_num    
!      REAL(8) stavg,ftavg
!      REAL(8) time, alphag_sum, ia_sum, vg_sum, dsm_sum
!      REAL(8) hdistance
!      REAL(8) radius(2), h(2), theta, radiusn, radiusmax,radiusmin, height, hzero,hzloc,lzloc
!      REAL(8) scelli, rcelli, minyval
!      REAL(8) xx,yy,zz,rr
!      REAL(8) vol_sum(2),tl_vol_sum(2),tl_vol_avg(2),tl_min(2),boron_vol_sum(2),boron_max(2),boron_vol_avg(2)
!      LOGICAL,SAVE:: initial_p      
!      DATA  Nradius,hzero,h(1),h(2),radius(1),radius(2)/2, 0.0d0, -1.2195d0, -0.3105d0, 0.437d0, 0.498d0/ !rocom:0.376,1.69, -317.5, -1185.5, 0.437, 0.5
!      DATA initial_p / .true./
!
!      REAL(8),SAVE, ALLOCATABLE:: vli_sum(:,:),cboroni_sum(:)
!      REAL(8),SAVE, ALLOCATABLE:: vli_ave(:,:),cboroni_ave(:)
!!      REAL(8),SAVE, ALLOCATABLE:: tli_sum(:),tli_ave(:),pi_sum(:),pi_ave(:)
!      REAL(8), ALLOCATABLE::vlp(:,:)
!      REAL(8), ALLOCATABLE::cboron(:)
!      REAL(8), ALLOCATABLE::volall(:)
!      REAL(8), ALLOCATABLE::temp(:),tempall(:)
!!      REAL(8), ALLOCATABLE::tlall(:),pr(:),xfc_all(:,:,:)
!!
!      ALLOCATE(vlp(ndim,ncell_fluid_all))
!      ALLOCATE(cboron(ncell_fluid_all))
!      ALLOCATE(volall(ncell_fluid_all))
!      ALLOCATE(tempall(ncell_fluid_all))
!      ALLOCATE(temp(ncell_fluid))
!!      ALLOCATE(tlall(ncell_fluid_all)      
!!      ALLOCATE(pr(ncell_fluid_all))      
!!      ALLOCATE(xfc_all(ndim,ns,ncell_fluid_all)) 
!!
!      IF(np.gt.1)THEN   
!         temp(1:ncell_fluid)=vl_n(1,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(1,:)=tempall(:)
!         temp(1:ncell_fluid)=vl_n(2,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(2,:)=tempall(:)
!         temp(1:ncell_fluid)=vl_n(ndim,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(ndim,:)=tempall(:)      
!         CALL allgatherv_r(cell%cboron,cboron,ncell_fluid,ncell_fluid_all,0)
!         CALL allgatherv_r(vol,volall,ncell_fluid,ncell_fluid_all,0)
!!        
!!         DO i=1, ndim
!!            DO j=1, ns
!!               temp(1:ncell_fluid)=xfc(i,j,1:ncell_fluid)
!!               CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!!               xfc_all(i,j,:)=tempall(:)
!!            ENDDO         
!!         ENDDO
!!         CALL allgatherv_r(cell%tl,tlall,ncell_fluid,ncell_fluid_all,0)
!!         CALL allgatherv_r(cell%p,pr,ncell_fluid,ncell_fluid_all,0)
!      ELSE
!         vlp(1:ndim,:)=vl_n(1:ndim,:)
!         cboron(:)=cell%cboron(:)
!         volall(:)=vol(:)
!!         tlall(:)=cell%tl(:)
!!         pr(:)=cell%p(:)         
!!         xfc_all(:,:,:)=xfc(:,:,:)         
!      ENDIF   
!
!      IF(myrank.gt.0)goto 100
!!      
!!.....DO just one time     
!      IF(ndim.ne.3)RETURN 
!      
!      IF (initial_p) THEN
!          initial_p = .false.
!          ALLOCATE(vli_sum(ndim,ncell_fluid_all),cboroni_sum(ncell_fluid_all))
!          ALLOCATE(vli_ave(ndim,ncell_fluid_all),cboroni_ave(ncell_fluid_all))
!!          ALLOCATE(tli_sum(ncell_fluid_all),tli_ave(ncell_fluid_all),pi_sum(ncell_fluid_all),pi_ave(ncell_fluid_all))
!          sum_num=0
!          cboroni_sum(:)=0.0d0
!          vli_sum(:,:)=0.0d0  
!!          tli_sum(:)=0.0d0          
!!          pi_sum(:)=0.0d0                 
!      ENDIF
!!....initialize the variable for time-averaging
!      IF(initial_p.eq..false.)THEN
!         IF(time.ge.stavg .and. time.le.ftavg)THEN
!             time_avg_opt=1
!             sum_num=sum_num+1
!             DO i=1,ncell_fluid_all                
!                vli_sum(i,:)=vli_sum(i,:)+vlp(:,i)
!                cboroni_sum(i)=cboroni_sum(i)+cboron(i)
!!                tli_sum(i)=tli_sum(i)+tlall(i)-273.15d0                
!!                pi_sum(i)=pi_sum(i)+pr(i)/1.0d6                
!                
!!                tli_ave(i)=tli_sum(i)/sum_num
!!                pi_ave(i)=pi_sum(i)/sum_num
!                cboroni_ave(i)=cboroni_sum(i)/sum_num
!                vli_ave(i,:)=vli_sum(i,:)/sum_num
!             ENDDO
!          ELSE
!             time_avg_opt=0   
!             sum_num=1         
!             DO i=1,ncell_fluid_all
!                vli_sum(i,:)=vlp(:,i)
!                cboroni_sum(i)=cboron(i)
!!                pi_sum(i)=pr(i)/1.0d6
!!                tli_sum(i)=tlall(i)-273.15d0                
!                
!!                tli_ave(i)=tli_sum(i)/sum_num
!!                pi_ave(i)=pi_sum(i)/sum_num
!                cboroni_ave(i)=cboroni_sum(i)/sum_num
!                vli_ave(i,:)=vli_sum(i,:)/sum_num
!            ENDDO   
!         ENDIF     
!      ENDIF    
!!
!!.....initialize the variables for space-averaging
!!
!      vol_sum(:)=0.0d0
!!      tl_vol_sum(:)=0.d0 
!!      tl_min(:)=10000.0d0
!      boron_vol_sum(:)=0.0d0
!      boron_vol_avg(:)=0.0d0
!      boron_max(:)=0.0d0
!!      
!!.....for the files for every PRINT step
!!
!      title(1)='"A:Radius",'  
!      title(2)='"B:Theta",'  
!      title(3)='"C:Height",'
!      title(4)=''                 !'"D:Tliquid"'
!      title(5)='"E:Cboron"'            
!      title(6)='"F:Vliquidx",'  
!      title(7)='"G:Vliquidy",'  
!      title(8)='"H:Vliquidz",'  
!      title(9)=''                 !'"I:Pressure"' 
!      title(10)=''                !'"J:Tliquid-K"'      
!      IF(prnopt)THEN
!         iitime=100000+time*100
!         IF(myrank.eq.0)WRITE(sstime,*)iitime
!         sstime=adjustl(sstime)
!         filename(1)='DC_low_Av_'//trim(sstime)//'.dat' 
!         filename(2)='DC_up_Av_'//trim(sstime)//'.dat' 
!         filename(3)='DC_low_Ins_'//trim(sstime)//'.dat' 
!         filename(4)='DC_up_Ins_'//trim(sstime)//'.dat' 
!         IF(myrank.eq.0)THEN
!            OPEN(61,file=filename(1)) !
!            OPEN(62,file=filename(2)) !
!            OPEN(63,file=filename(3)) !
!            OPEN(64,file=filename(4)) !
!         ENDIF
!      ENDIF         
!!
!!.....for the 2 heights of DC
!!
!      DO level = 1, Nradius 
!         loop = 0
!         levfile = 60 + level
!         IF(myrank.eq.0)THEN
!            IF(prnopt)WRITE(levfile, 8003)(title(i), i=1, 10)
!            IF(prnopt)WRITE(levfile+2, 8003)(title(i), i=1, 10)
!         ENDIF
!  8003   FORMAT('variables=',10A16)
!
!!
!!........for the i-th cell
!!         
!         DO i = 1, ncell_fluid_all
!!      
!!...........Calculate the height of the i cell.
!!
!            hzloc=xloc_tmp(ndim,i)
!            lzloc=xloc_tmp(ndim,i)
!            DO j=1,num_neigh_tmp(i)
!               hzloc=DMAX1(hzloc,xfc_tmp(ndim,j,i))
!               lzloc=DMIN1(lzloc,xfc_tmp(ndim,j,i))
!            ENDDO     
!!
!!...........Check the height
!!
!            IF(hzloc.le.h(level).or.lzloc.gt.h(level)) CYCLE            
!!
!!...........Check the radius
!!
!            wopt = 0
!            IF(ndim.eq.3) radiusn=DSQRT(xloc_tmp(1,i)**2.d0+xloc_tmp(2,i)**2.0d0)    
!            IF(radiusn.le.radius(2).and.radiusn.ge.radius(1)) wopt=1            
!!
!            IF(wopt.ne.1)CYCLE
!    rad_group(i)=lzloc
!            loop=loop+1
!!...........sum of tl for the 2-d plane          
!!            tl_vol_sum(level)=tl_vol_sum(level)+(tlall(i)-273.15d0)*volall(i)
!!            IF(tl_min(level).gt.(tlall(i)-273.15d0))tl_min(level)=tlall(i)-273.15d0   
!!       
!            boron_vol_sum(level)=boron_vol_sum(level)+cboron(i)*volall(i)
!            IF(boron_max(level).lt.cboron(i))boron_max(level)=cboron(i) 
!            vol_sum(level)=vol_sum(level)+volall(i)                       
!!
!IF(prnopt)THEN                
!!..............calculate theta using x,y               
!               xx=xloc_tmp(1,i)
!               yy=xloc_tmp(2,i)
!               zz=xloc_tmp(3,i)
!               rr=(xx**2.d0+yy**2.d0)**0.5d0
!               IF(rr.ne.0.d0)THEN
!                  theta=xx/rr
!                  theta=acos(theta) !theta=0~180
!               ELSE
!                  PRINT *,'rr is 0 in SUBROUTINE Rocom_output_cylinder(time)!!!'
!                  STOP
!               ENDIF               
!               IF(yy.gt.0.d0)THEN
!                  theta=theta       !0~180 ==> 0-180
!               ELSEIF(yy.le.0.d0)THEN
!                  theta=-theta !2*pi-theta  !180-0 ==> 180-360
!               ELSE
!                  PRINT *,'theta cannot be defined in SUBROUTINE Rocom_output_cylinder(time)!!!'
!                  STOP
!               ENDIF
!            yy=theta*180.d0/pi                             
!!
!            IF(myrank.eq.0)WRITE(levfile,5001)rr,yy,zz,cboroni_ave(i),&                  !tli_ave(i),
!                                   vli_ave(1,i),vli_ave(2,i),vli_ave(ndim,i)             !,pi_ave(i),tli_ave(i)+273.15d0                          
!            IF(myrank.eq.0)WRITE(levfile+2,5001)rr,yy,zz,cboron(i),&                     !,tlall(i)-273.15d0
!                                   vlp(1,i),vlp(2,i),vlp(ndim,i)                         !,pr(i),tlall(i)
!ENDIF !prnopt    
!!        
!         ENDDO !i
!!.....PRINT the DATA 
!      IF(prnopt)THEN 
!         IF(myrank.eq.0)WRITE(65,*)'time=',time
!         IF(myrank.eq.0)CLOSE(65)   
!         IF(myrank.eq.0)WRITE(66,*)'time=',time
!         IF(myrank.eq.0)CLOSE(66)  
!      ENDIF         
!      ENDDO !level
! 5000 FORMAT(1x,'zone n=4, e=1, DATApacking=point, zonetype=FEQUADRILATERAL')
! 5001 FORMAT(1x,10e16.5)
! 5002 FORMAT(1x,3e16.5)
!
!!
!!.....averaging the tl for the 2-d plane
!!        
!      DO i=1, nradius
!        IF(vol_sum(i).gt.0.d0)THEN
!           boron_vol_avg(i)=boron_vol_sum(i)/vol_sum(i)
!        ELSE
!           IF(myrank.eq.0)PRINT *,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
!        ENDIF           
!      ENDDO    
!      tl_space_avg(1)=boron_vol_avg(2)  
!      tl_space_avg(2)=boron_vol_avg(1)     
!      tl_space_min(1)=boron_max(2)
!      tl_space_min(2)=boron_max(1)  
!!
!      IF(vv_prob.eq.'pts16')THEN
!!         tl_space_avg(4)=cboron(3722028)          !468 NewCL_BL
!!         tl_space_avg(5)=cboron(3667428)          !Lower DC
!!         tl_space_min(6)=xloc_tmp(1,3722028)
!!         tl_space_min(7)=xloc_tmp(2,3722028)
!!         tl_space_min(8)=xloc_tmp(3,3722028)
!!         tl_space_avg(6)=xloc_tmp(1,3667428)
!!         tl_space_avg(7)=xloc_tmp(2,3667428)
!!         tl_space_avg(8)=xloc_tmp(3,3667428)
!!
!         tl_space_avg(4)=cboron(2319900)         !1243 NewCL_BL
!         tl_space_avg(5)=cboron(2538015)         !Lower DC
!         tl_space_min(6)=xloc_tmp(1,2319900)
!         tl_space_min(7)=xloc_tmp(2,2319900)
!         tl_space_min(8)=xloc_tmp(3,2319900)
!         tl_space_avg(6)=xloc_tmp(1,2538015)
!         tl_space_avg(7)=xloc_tmp(2,2538015)
!         tl_space_avg(8)=xloc_tmp(3,2538015)         
!      ENDIF
!!      tl_space_avg(4)=cboron(29234)           !300,000 grids
!!      tl_space_avg(5)=cboron(46254)
!!      tl_space_min(6)=xloc_tmp(1,29234)
!!      tl_space_min(7)=xloc_tmp(2,29234)
!!      tl_space_min(8)=xloc_tmp(3,29234)
!!      tl_space_avg(6)=xloc_tmp(1,46254)
!!      tl_space_avg(7)=xloc_tmp(2,46254)
!!      tl_space_avg(8)= xloc_tmp(3,46254)
!!
!100   continue   
!   
!      DEALLOCATE(vlp)
!!      DEALLOCATE(pr)
!!      DEALLOCATE(tlall)
!      DEALLOCATE(volall)
!      DEALLOCATE(cboron)  
!      DEALLOCATE(temp,tempall)   
!!      DEALLOCATE(xfc_all)  
!!    
!      RETURN
!      END
!
!!      
!!.....These SUBROUTINE is very special print routine for secIFic problem
!!.....Don't USE for another problem
!!
!      SUBROUTINE rocom_output_coreinlet_pts16(time,prnopt,stavg,ftavg)
!!
!!.....Save lateral output for Origin
!!   
!      USE VOL_DATA                     
!      USE Zparam          , ONLY: ndim,ns,pi
!      USE Zcore           , ONLY: np,myrank
!      USE Zcoord1         , ONLY: xloc_tmp,xfc_tmp
!      !USE Zcoord2         , ONLY: xfc
!      USE Zcoord3         , ONLY: vol
!      USE Znum_cell       , ONLY: num_neigh,num_neigh_tmp  
!      USE Zpress          , ONLY: p               
!      USE Zrocom_specific , ONLY: tl_space_avg,tl_space_min,rad_group
!      USE Zvector         , ONLY: vl_n
!      USE Zzone           , ONLY: ncell_fluid,ncell_fluid_all         
!      USE Znode           , ONLY: num_cell_node,cell_node,xnode,nd_max   
!!
!      IMPLICIT NONE 
!!
!      CHARACTER title(30) * 20
!      CHARACTER filename(6)*20, sstime*16
!      INTEGER prnopt     
!      INTEGER i, j, ic, ix, wopt, loop, sortopt, level
!      INTEGER Nradius      
!      INTEGER time_avg_opt
!      INTEGER iitime   
!      INTEGER,SAVE:: sum_num    
!      REAL(8) stavg,ftavg
!      REAL(8) time, alphag_sum, ia_sum, vg_sum, dsm_sum
!      REAL(8) hdistance
!      REAL(8) radius, h1, h2, theta, radiusn, radiusmax,radiusmin, height, hzero,hzloc,lzloc
!      REAL(8) scelli, rcelli, minyval
!      REAL(8) xx,yy,zz,rr!
!      REAL(8) vol_sum(2),tl_vol_sum(2),tl_vol_avg(2),tl_min(2),boron_vol_sum(2),boron_max(2),boron_vol_avg(2)
!      LOGICAL,SAVE:: initial_p      
!      DATA  Nradius,hzero,h1,h2,radius/1, 0.0d0, -1.2515, -1.2515, 0.437/ !rocom:0.376,1.69, -317.5, -1185.5, 0.437, 0.5
!      DATA initial_p / .true./
!
!      REAL(8),SAVE, ALLOCATABLE:: vli_sum(:,:),cboroni_sum(:)
!      REAL(8),SAVE, ALLOCATABLE:: vli_ave(:,:),cboroni_ave(:)
!!      REAL(8),SAVE, ALLOCATABLE:: tli_sum(:),tli_ave(:),pi_sum(:),pi_ave(:)
!      REAL(8), ALLOCATABLE::vlp(:,:)
!      REAL(8), ALLOCATABLE::cboron(:)
!      REAL(8), ALLOCATABLE::volall(:)
!      REAL(8), ALLOCATABLE::temp(:),tempall(:)
!!      REAL(8), ALLOCATABLE::tlall(:),pr(:),xfc_all(:,:,:)
!!
!      ALLOCATE(vlp(ndim,ncell_fluid_all))
!      ALLOCATE(cboron(ncell_fluid_all))
!      ALLOCATE(volall(ncell_fluid_all))
!      ALLOCATE(tempall(ncell_fluid_all))
!      ALLOCATE(temp(ncell_fluid))
!!      ALLOCATE(tlall(ncell_fluid_all)      
!!      ALLOCATE(pr(ncell_fluid_all))      
!!      ALLOCATE(xfc_all(ndim,ns,ncell_fluid_all)) 
!!
!      IF(np.gt.1)THEN   
!         temp(1:ncell_fluid)=vl_n(1,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(1,:)=tempall(:)
!         temp(1:ncell_fluid)=vl_n(2,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(2,:)=tempall(:)
!         temp(1:ncell_fluid)=vl_n(ndim,1:ncell_fluid)
!         CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!         vlp(ndim,:)=tempall(:)      
!!        
!!         DO i=1, ndim
!!            DO j=1, ns
!!               temp(1:ncell_fluid)=xfc(i,j,1:ncell_fluid)
!!               CALL allgatherv_r(temp,tempall,ncell_fluid,ncell_fluid_all,0)
!!               xfc_all(i,j,:)=tempall(:)
!!            ENDDO         
!!         ENDDO
!!
!!         CALL allgatherv_r(cell%p,pr,ncell_fluid,ncell_fluid_all,0)
!!         CALL allgatherv_r(cell%tl,tlall,ncell_fluid,ncell_fluid_all,0)
!         CALL allgatherv_r(cell%cboron,cboron,ncell_fluid,ncell_fluid_all,0)
!         CALL allgatherv_r(vol,volall,ncell_fluid,ncell_fluid_all,0)
!      ELSE
!         vlp(1:ndim,:)=vl_n(1:ndim,:)
!         cboron(:)=cell%cboron(:)
!         volall(:)=vol(:)
!!         tlall(:)=cell%tl(:)
!!         pr(:)=cell%p(:)         
!!         xfc_all(:,:,:)=xfc(:,:,:) 
!      ENDIF   
!
!      IF(myrank.gt.0)goto 100
!!      
!!.....DO just one time     
!      IF(ndim.ne.3)RETURN 
!      
!      IF (initial_p) THEN
!          initial_p = .false.
!          ALLOCATE(vli_sum(ndim,ncell_fluid_all),cboroni_sum(ncell_fluid_all))
!          ALLOCATE(vli_ave(ndim,ncell_fluid_all),cboroni_ave(ncell_fluid_all))
!!          ALLOCATE(tli_sum(ncell_fluid_all),tli_ave(ncell_fluid_all),pi_sum(ncell_fluid_all),pi_ave(ncell_fluid_all))
!          sum_num=0
!!          tli_sum(:)=0.0d0
!          cboroni_sum(:)=0.0d0
!!          pi_sum(:)=0.0d0
!          vli_sum(:,:)=0.0d0         
!      ENDIF
!!....initialize the variable for time-averaging
!      IF(initial_p.eq..false.)THEN
!         IF(time.ge.stavg .and. time.le.ftavg)THEN
!             time_avg_opt=1
!             sum_num=sum_num+1
!             DO i=1,ncell_fluid_all
!!                tli_sum(i)=tli_sum(i)+tlall(i)-273.15d0
!!                pi_sum(i)=pi_sum(i)+pr(i)/1.0d6
!                cboroni_sum(i)=cboroni_sum(i)+cboron(i)
!                vli_sum(i,:)=vli_sum(i,:)+vlp(:,i)
!                
!!                tli_ave(i)=tli_sum(i)/sum_num
!!                pi_ave(i)=pi_sum(i)/sum_num
!                cboroni_ave(i)=cboroni_sum(i)/sum_num
!                vli_ave(i,:)=vli_sum(i,:)/sum_num
!             ENDDO
!          ELSE
!             time_avg_opt=0   
!             sum_num=1         
!             DO i=1,ncell_fluid_all
!!                tli_sum(i)=tlall(i)-273.15d0
!!                pi_sum(i)=pr(i)/1.0d6
!                cboroni_sum(i)=cboron(i)
!                vli_sum(i,:)=vlp(:,i)                
!!                
!!                tli_ave(i)=tli_sum(i)/sum_num
!!                pi_ave(i)=pi_sum(i)/sum_num
!                cboroni_ave(i)=cboroni_sum(i)/sum_num
!                vli_ave(i,:)=vli_sum(i,:)/sum_num
!            ENDDO   
!         ENDIF     
!      ENDIF    
!!
!!.....initialize the variables for space-averaging
!!
!      vol_sum(:)=0.0d0
!!      tl_vol_sum(:)=0.d0 
!!      tl_min(:)=10000.0d0
!      boron_vol_sum(:)=0.0d0
!      boron_vol_avg(:)=0.0d0
!      boron_max(:)=0.0d0
!!      
!!.....for the files for every PRINT step
!!
!      title(1)='"A:Radius",'  
!      title(2)='"B:Theta",'  
!      title(3)='"C:Height",'
!      title(4)=''             !'"D:tl",'
!      title(5)='"E:cBoron"'
!      title(6)='"F:Vliquidx"'            
!      title(7)='"G:Vliquidy",'  
!      title(8)='"H:Vliquidz",'  
!      title(9)=''             !'"I:Pressure",'  
!      title(10)=''            !'"J:Tliquid-K"'    
!        
!      IF(prnopt)THEN
!         iitime=100000+time*100
!         IF(myrank.eq.0)WRITE(sstime,*)iitime
!         sstime=adjustl(sstime)
!         filename(1)='core_Av_'//trim(sstime)//'.dat' 
!         filename(2)='core_ins_'//trim(sstime)//'.dat' 
!         IF(myrank.eq.0)THEN
!            OPEN(65,file=filename(1)) !
!            OPEN(66,file=filename(2)) !
!         ENDIF
!      ENDIF         
!!
!!.....for the 2 radius
!!
!      DO level = 1, Nradius 
!         loop = 0
!         IF(myrank.eq.0)THEN
!            IF(prnopt)WRITE(65, 8003)(title(i), i=1, 10)
!            IF(prnopt)WRITE(66, 8003)(title(i), i=1, 10)
!         ENDIF
!  8003   FORMAT('variables=',10A16)
!
!!
!!........for the i-th cell
!!         
!         DO i = 1, ncell_fluid_all
!!      
!!...........Calculate the height of the i cell.
!!
!            hzloc=xloc_tmp(ndim,i)
!            lzloc=xloc_tmp(ndim,i)
!            DO j=1,num_neigh_tmp(i)
!               hzloc=DMAX1(hzloc,xfc_tmp(ndim,j,i))
!               lzloc=DMIN1(lzloc,xfc_tmp(ndim,j,i))
!            ENDDO     
!!
!!...........Check the height
!!
!            IF(hzloc.le.h1.or.lzloc.gt.h1) CYCLE            
!!
!!...........Check the radius
!!
!            wopt = 0
!            IF(ndim.eq.3) radiusn=DSQRT(xloc_tmp(1,i)**2.d0+xloc_tmp(2,i)**2.0d0)    
!            IF(radiusn.le.radius) wopt=1            
!!
!            IF(wopt.ne.1)CYCLE
!    rad_group(i)=-2.0d0
!            loop=loop+1
!!...........sum of tl for the 2-d plane                 
!!            tl_vol_sum(level)=tl_vol_sum(level)+(tlall(i)-273.15d0)*volall(i)
!!            IF(tl_min(level).gt.(tlall(i)-273.15d0))tl_min(level)=tlall(i)-273.15d0    
!!
!            boron_vol_sum(level)=boron_vol_sum(level)+cboron(i)*volall(i)
!            IF(boron_max(level).lt.cboron(i))boron_max(level)=cboron(i) 
!            vol_sum(level)=vol_sum(level)+volall(i)                                 
!!
!IF(prnopt)THEN 
!
!!..............calculate theta using x,y               
!               xx=xloc_tmp(1,i)
!               yy=xloc_tmp(2,i)
!               zz=xloc_tmp(3,i)
!               rr=(xx**2.d0+yy**2.d0)**0.5d0
!               IF(rr.ne.0.d0)THEN
!                  theta=xx/rr
!                  theta=acos(theta) !theta=0~180
!               ELSE
!                  PRINT *,'rr is 0 in SUBROUTINE Rocom_output_cylinder(time)!!!'
!                  STOP
!               ENDIF               
!               IF(yy.gt.0.d0)THEN
!                  theta=theta       !0~180 ==> 0-180
!               ELSEIF(yy.le.0.d0)THEN
!                  theta=-theta !2*pi-theta  !180-0 ==> 180-360
!               ELSE
!                  PRINT *,'theta cannot be defined in SUBROUTINE Rocom_output_cylinder(time)!!!'
!                  STOP
!               ENDIF
!            yy=theta*180.d0/pi  
!!
!            IF(myrank.eq.0)WRITE(65,5001)rr,yy,zz,cboroni_ave(i),&                           !,tli_ave(i)
!                                   vli_ave(1,i),vli_ave(2,i),vli_ave(ndim,i)                 !,pi_ave(i),tli_ave(i)+273.15d0                          
!            IF(myrank.eq.0)WRITE(66,5001)rr,yy,zz,cboron(i),&                                !,tlall(i)-273.15d0
!                                   vlp(1,i),vlp(2,i),vlp(ndim,i)                             !,pr(i),tlall(i)
!ENDIF !prnopt    
!!        
!         ENDDO !i
!!.....PRINT the DATA 
!      IF(prnopt)THEN 
!         IF(myrank.eq.0)WRITE(65,*)'time=',time
!         IF(myrank.eq.0)CLOSE(65)   
!         IF(myrank.eq.0)WRITE(66,*)'time=',time
!         IF(myrank.eq.0)CLOSE(66)  
!      ENDIF         
!      ENDDO !level
! 5000 FORMAT(1x,'zone n=4, e=1, DATApacking=point, zonetype=FEQUADRILATERAL')
! 5001 FORMAT(1x,10e16.5)
! 5002 FORMAT(1x,3e16.5)
!
!!
!!.....averaging the tl for the 2-d plane
!!        
!!      DO i=1, nradius
!!        IF(vol_sum(i).gt.0.d0)THEN
!!           tl_vol_avg(i)=tl_vol_sum(i)/vol_sum(i)
!!        ELSE
!!           IF(myrank.eq.0)PRINT *,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
!!        ENDIF           
!!      ENDDO    
!!      tl_space_avg(1)=tl_vol_avg(1)     
!!      tl_space_min(1)=tl_min(1)      
!!
!      DO i=1, nradius
!        IF(vol_sum(i).gt.0.d0)THEN
!           boron_vol_avg(i)=boron_vol_sum(i)/vol_sum(i)
!        ELSE
!           IF(myrank.eq.0)PRINT *,'vol_sum(i) is zero for',vol_sum(i),'level in SUBROUTINE rocom_output_cylindeer'           
!        ENDIF           
!      ENDDO    
!      tl_space_avg(3)=boron_vol_avg(1)  
!      tl_space_min(3)=boron_max(1)       
! 
!100   continue   
!   
!      DEALLOCATE(vlp)
!!      DEALLOCATE(pr)
!!      DEALLOCATE(tlall)
!      DEALLOCATE(volall)
!      DEALLOCATE(cboron)  
!      DEALLOCATE(temp,tempall)   
!!      DEALLOCATE(xfc_all)  
!
!!      DEALLOCATE(tli_sum,vli_sum,pi_sum,cboroni_sum)
!!      DEALLOCATE(tli_ave,vli_ave,pi_ave,cboroni_ave)
!    
!      RETURN
!      END
