!
      SUBROUTINE radiation_component_init_zone(nn,nzone_tmp,xloc_tmp) !pik-radiation_component
!   
      USE Zparam           ,ONLY: ndim
      USE Zcore            ,ONLY: myrank,np
      USE Zrad_comp  
!      
      IMPLICIT NONE
!      
      INTEGER:: nn !the geometrical number of cells , solid cell index
      INTEGER:: nzone_tmp(nn)
      REAL(8):: xloc_tmp(nn,ndim)
!      
      INTEGER:: i,j,k
      INTEGER:: pair,loop,ix,k1,k2,iloc
      INTEGER:: sort_opt
      REAL(8):: x(ndim),x1(ndim),d1,d2,dm1,dm2
!
!.....2d rod input
      INTEGER:: num_ch,nz0_2d,ch_opt,idummy,io     
      INTEGER,ALLOCATABLE::nf_input(:,:),rod_in_fluid(:),cell_idx_old(:,:,:)
!
!.....screening input     
      INTEGER:: major_dir,minor_dir,idistance
      INTEGER,ALLOCATABLE::count_gabage(:,:),minor_rad_zone(:),major_rad_zone(:)
      INTEGER,ALLOCATABLE:: rad_zone(:,:),set_idx(:),pair_idx(:),count_rad_zone(:,:)
      REAL(8),ALLOCATABLE::  xbottom(:),xtop(:),xleft(:),xright(:)
      INTEGER,ALLOCATABLE:: nlayer(:,:),mylayer(:,:)
      REAL(8),ALLOCATABLE:: d_idx_tmp(:,:,:),dm_idx_tmp(:,:,:)
      
!
!.....Read radiation information
!     
      OPEN(333,file='radiation_component.in',status='old')
      READ(333,*)input_opt
      READ(333,*)nset,major_dir,minor_dir
      ALLOCATE(nproperty(2,nset))
      ALLOCATE(nlayer(2,nset))
      ALLOCATE(mylayer(2,nset))
      ALLOCATE(epsil(2,nset))
      ALLOCATE(viewf(2,2,nset))      
      ALLOCATE(area(2,nset))
      ALLOCATE(rad_zone(2,nset))
!      
      ALLOCATE(nsize(2,nset))      
!      
      nproperty=0
      nsize=0
      epsil=0.0d0
      area=0.0d0
      DO j=1,nset
         READ(333,*)nproperty(1,j),nproperty(2,j) !,nsize(j)              !fluid solid size 
         READ(333,*)nlayer(1,j),nlayer(2,j),mylayer(1,j),mylayer(2,j)    !minor_dir
         READ(333,*)epsil(1,j),epsil(2,j),area(1,j),area(2,j)           !epilson, radiation, area
         READ(333,*)viewf(1,1,j),viewf(1,2,j),viewf(2,1,j),viewf(2,2,j) !view factor
         READ(333,*)rad_zone(1,j)
         READ(333,*)rad_zone(2,j)
      ENDDO
!
!.....For special output      
      READ(333,*)nloc_out
      ncell_out=nloc_out         
!         
      ALLOCATE(xloc_out_tmp(nloc_out,ndim)) !input
      ALLOCATE(set_idx(nloc_out))
      ALLOCATE(pair_idx(nloc_out))
!         
      ALLOCATE(cell_out_tmp(nloc_out,2))    !output
      ALLOCATE(frac_out_tmp(nloc_out,2))
      ALLOCATE(nproperty_out_tmp(nloc_out))
      xloc_out_tmp=0.0d0
      set_idx=0
      pair_idx=0
      cell_out_tmp=0
      frac_out_tmp=0.0d0
      nproperty_out_tmp=0
!      
      DO i=1,nloc_out
         READ(333,*)(xloc_out_tmp(i,ix),ix=1,3),set_idx(i),pair_idx(i)
      ENDDO  
      CLOSE(333)
!         
      IF(myrank.eq.0)THEN         
!
!........read ht_str_2d.in to check 2d_rod
!         
         OPEN(52,file='ht_str_2d.in',status='old',iostat=io)
         IF(io.gt.0) RETURN
         READ(52,*) num_ch,nz0_2d,ch_opt
         IF(nz0_2d.eq.0)RETURN
         IF(ch_opt.eq.1)THEN
            DO i=1,num_ch
               READ(52,*) idummy
            ENDDO
         ELSE
            READ(52,*) idummy
         ENDIF     
         ALLOCATE(nf_input(num_ch,nz0_2d))
         DO i=1,num_ch
            READ(52,*) (nf_input(i,j),j=1,nz0_2d)
         ENDDO
         ALLOCATE(rod_in_fluid(nn))
         rod_in_fluid(:)=0
         DO i=1,num_ch
            DO j=1,nz0_2d
               rod_in_fluid(nf_input(i,j))=1
            ENDDO
         ENDDO   
         DEALLOCATE(nf_input)
         CLOSE(52)
!
!........pick up cell_idx_tmp
!
!
!........obtain the number of cells, count_rad_zone(pair,j)         
         ALLOCATE(count_rad_zone(2,nset))
         count_rad_zone(:,:)=0.0d0
         DO k=1,nn
            DO pair=1,2
               DO j=1,nset
                  IF(nzone_tmp(k).eq.rad_zone(pair,j))THEN
                     IF(nproperty(pair,j).eq.1.and.rod_in_fluid(k).eq.0)CYCLE
                     count_rad_zone(pair,j)=count_rad_zone(pair,j)+1
                     i=count_rad_zone(pair,j)
                  ENDIF
               ENDDO  
            ENDDO   
         ENDDO 
!
!........obtain the cell index
         nsize_max=0
         DO pair=1,2
            DO j=1,nset
               nsize_max=MAX(nsize_max,count_rad_zone(pair,j))
            ENDDO
         ENDDO    
         ALLOCATE(cell_idx_tmp(nsize_max,2,nset))
         ALLOCATE(d_idx_tmp(nsize_max,2,nset))
!         
         cell_idx_tmp(:,:,:)=0
         d_idx_tmp(:,:,:)=100.0d0
         count_rad_zone(:,:)=0.0d0
         DO k=1,nn
            DO pair=1,2
               DO j=1,nset
                  IF(nzone_tmp(k).eq.rad_zone(pair,j))THEN
                     IF(nproperty(pair,j).eq.1.and.rod_in_fluid(k).eq.0)CYCLE
                     count_rad_zone(pair,j)=count_rad_zone(pair,j)+1
                     i=count_rad_zone(pair,j)
                     cell_idx_tmp(i,pair,j)=k
                  ENDIF
               ENDDO  
            ENDDO   
         ENDDO          
         DEALLOCATE(rod_in_fluid)
!
!........set major rad_zone of a pair solid         
         ALLOCATE(major_rad_zone(nset))
         ALLOCATE(minor_rad_zone(nset))
         DO j=1,nset
            nsize(1,j)=MIN(count_rad_zone(1,j),count_rad_zone(2,j))
            nsize(2,j)=nsize(1,j)
            IF(count_rad_zone(1,j).le.count_rad_zone(2,j))THEN
               major_rad_zone(j)=1
               minor_rad_zone(j)=2
            ELSE
               major_rad_zone(j)=2
               minor_rad_zone(j)=1
            ENDIF   
         ENDDO
!
!........screen 2nd pair according to the bottom and top of major zone along major_dir
!         
         ALLOCATE(xbottom(nset))
         ALLOCATE(xtop(nset))
         ALLOCATE(xleft(nset))
         ALLOCATE(xright(nset))         
         ALLOCATE(count_gabage(2,nset))
         xbottom(:)=+100.0d0
         xtop(:)=-100.0d0
         xleft(:)=+100.0d0
         xright(:)=-100.0d0         
         count_gabage(:,:)=0
!    
!........get bottom and top of major zone         
         DO j=1,nset
            pair=major_rad_zone(j)
            DO i=1,count_rad_zone(pair,j)
               k1=cell_idx_tmp(i,pair,j)             
               xbottom(j)=MIN(xbottom(j),xloc_tmp(k1,major_dir))
               xtop(j)=MAX(xtop(j),xloc_tmp(k1,major_dir))
               xleft(j)=MIN(xleft(j),xloc_tmp(k1,minor_dir))
               xright(j)=MAX(xright(j),xloc_tmp(k1,minor_dir))
            ENDDO
         ENDDO 
         xbottom(:)=xbottom(:)-1.0d-3
         xtop(:)=xtop(:)+1.0d-3
         xleft(:)=xleft(:)-1.0d-3
         xright(:)=xright(:)+1.0d-3         
!
!........color the cells beyond the bottom and top         
         DO j=1,nset
            pair=minor_rad_zone(j)
            DO i=1,count_rad_zone(pair,j)
               k1=cell_idx_tmp(i,pair,j)             
               IF(xloc_tmp(k1,major_dir).lt.xbottom(j).or.&
                  xloc_tmp(k1,major_dir).gt.xtop(j))THEN
                  cell_idx_tmp(i,pair,j)=-1
                  count_gabage(pair,j)=count_gabage(pair,j)+1
               ENDIF   
            ENDDO
         ENDDO 
!
!........fill the colored cells with the next cells         
         DO j=1,nset
            pair=minor_rad_zone(j)
            DO loop=1,10000
               sort_opt=0
               DO i=1,count_rad_zone(pair,j)
                  IF(cell_idx_tmp(i,pair,j).eq.-1)THEN
                     sort_opt=1
                     DO k=i,count_rad_zone(pair,j)-1
                        cell_idx_tmp(k,pair,j)=cell_idx_tmp(k+1,pair,j)
                        cell_idx_tmp(k+1,pair,j)=0
                     ENDDO   
                  ENDIF   
               ENDDO
               IF(sort_opt.eq.0)EXIT
            ENDDO   
            IF(sort_opt.eq.1)THEN
               WRITE(*,*)'1error during sorting cell_idx_tmp pair!!!'
               PAUSE
               STOP
            ELSE   
               count_rad_zone(pair,j)=count_rad_zone(pair,j)-count_gabage(pair,j)
            ENDIF
         ENDDO
!         
         DO j=1,nset
            IF(count_rad_zone(1,j).ne.count_rad_zone(2,j))THEN
               WRITE(*,*)'2error in the number of cells!!!'
               !PAUSE
               !STOP
            ENDIF
         ENDDO   
!
!........re-calculate nsize_max because count_rad_zone is changed        
         nsize_max=0
         DO pair=1,2
            DO j=1,nset
               nsize_max=MAX(nsize_max,count_rad_zone(pair,j))
            ENDDO
         ENDDO  
!
!........sort according to the distance from 0 along the major_dir
!         
         !x1(major_dir)=0.0d0
         x1(major_dir)=minval(xbottom(1:nset))
         DO pair=1,2
            DO j=1,nset
               DO i=1,count_rad_zone(pair,j)
                  k1=cell_idx_tmp(i,pair,j)             
                  d1=0.0d0
                  d1=xloc_tmp(k1,major_dir)-x1(major_dir)
                  idistance=d1*1000                  
                  d_idx_tmp(i,pair,j)=d1/1000.0d0
               ENDDO
            ENDDO 
         ENDDO
!          
         DO pair=1,2
            DO j=1,nset
               DO loop=1,10000
                  sort_opt=0
                  DO i=1,count_rad_zone(pair,j)-1
                     k1=cell_idx_tmp(i,pair,j)
                     k2=cell_idx_tmp(i+1,pair,j)
                     d1=d_idx_tmp(i,pair,j)
                     d2=d_idx_tmp(i+1,pair,j)
                     IF(d1.gt.d2)THEN
                        cell_idx_tmp(i,pair,j)=k2
                        cell_idx_tmp(i+1,pair,j)=k1
                        d_idx_tmp(i,pair,j)=d2
                        d_idx_tmp(i+1,pair,j)=d1
                        sort_opt=1
                     ENDIF   
                  ENDDO   
                  IF(sort_opt.eq.0)EXIT
               ENDDO   
               IF(sort_opt.eq.1)THEN
                  WRITE(*,*)'3error during sorting cell_idx_tmp!!!'
                  PAUSE
                  STOP
               ENDIF   
            ENDDO
         ENDDO
!=========================================================================
!
!........sort according to the distance from 0 along the minor_dir
!     
      IF(minor_dir.ne.0)THEN         
         ALLOCATE(dm_idx_tmp(nsize_max,2,nset))
         x1(minor_dir)=minval(xleft(1:nset))
           
         DO pair=1,2
            DO j=1,nset
               DO i=1,count_rad_zone(pair,j)
                  k1=cell_idx_tmp(i,pair,j)             
                  dm1=0.0d0
                  dm1=xloc_tmp(k1,minor_dir)-x1(minor_dir)
                  idistance=dm1*1000
                  dm_idx_tmp(i,pair,j)=idistance/1000.0d0
               ENDDO
            ENDDO 
         ENDDO
!          
         DO pair=1,2
            DO j=1,nset
               DO loop=1,10000
                  sort_opt=0
                  DO i=1,count_rad_zone(pair,j)-1
                     k1=cell_idx_tmp(i,pair,j)
                     k2=cell_idx_tmp(i+1,pair,j)
                     d1=d_idx_tmp(i,pair,j)
                     d2=d_idx_tmp(i+1,pair,j)                     
                     dm1=dm_idx_tmp(i,pair,j)
                     dm2=dm_idx_tmp(i+1,pair,j)
                     IF(DABS(d1-d2).le.1.d-5)THEN !same level
                        IF(dm1.gt.dm2)THEN
                           cell_idx_tmp(i,pair,j)=k2
                           cell_idx_tmp(i+1,pair,j)=k1
                           dm_idx_tmp(i,pair,j)=dm2
                           dm_idx_tmp(i+1,pair,j)=dm1
                           sort_opt=1
                        ENDIF   
                     ENDIF   
                  ENDDO   
                  IF(sort_opt.eq.0)EXIT
               ENDDO   
               IF(sort_opt.eq.1)THEN
                  WRITE(*,*)'4error during sorting cell_idx_tmp!!!'
                  PAUSE
                  STOP
               ENDIF   
            ENDDO
         ENDDO
!
!........color the cells beyond the bottom and top along the minor_dir
         DO j=1,nset
            DO pair=1,2
               mylayer(pair,j)=mod(mylayer(pair,j),nlayer(pair,j))
            ENDDO
         ENDDO   
         count_gabage(:,:)=0
         DO j=1,nset
            DO pair=1,2
               IF(nlayer(pair,j).eq.1)CYCLE
               DO i=1,count_rad_zone(pair,j)
                  IF(mod(i,nlayer(pair,j)).ne.mylayer(pair,j))THEN            
                     cell_idx_tmp(i,pair,j)=-1
                     count_gabage(pair,j)=count_gabage(pair,j)+1
                  ENDIF   
               ENDDO
            ENDDO
         ENDDO 
!
!........fill the colored cells with the next cells         
         DO j=1,nset
            DO pair=1,2
               DO loop=1,10000
                  sort_opt=0
                  DO i=1,count_rad_zone(pair,j)
                     IF(cell_idx_tmp(i,pair,j).eq.-1)THEN
                        sort_opt=1
                        DO k=i,count_rad_zone(pair,j)-1
                           cell_idx_tmp(k,pair,j)=cell_idx_tmp(k+1,pair,j)
                           cell_idx_tmp(k+1,pair,j)=0
                        ENDDO   
                     ENDIF   
                  ENDDO
                  IF(sort_opt.eq.0)EXIT
               ENDDO   
               IF(sort_opt.eq.1)THEN
                  WRITE(*,*)'5error during sorting cell_idx_tmp pair!!!'
                  PAUSE
                  STOP
               ELSE   
                  count_rad_zone(pair,j)=count_rad_zone(pair,j)-count_gabage(pair,j)
               ENDIF
            ENDDO
         ENDDO
!         
         DO j=1,nset
            IF(count_rad_zone(1,j).ne.count_rad_zone(2,j))THEN
               WRITE(*,*)'6error in the number of cells!!!'
               PAUSE
               STOP
            ENDIF
         ENDDO   
!
!........re-calculate nsize_max because count_rad_zone is changed        
         nsize_max=0
         DO pair=1,2
            DO j=1,nset
               nsize_max=MAX(nsize_max,count_rad_zone(pair,j))
               nsize(pair,j)=count_rad_zone(pair,j)
            ENDDO
         ENDDO  
         
      ENDIF
!=========================================================================
         DEALLOCATE(minor_rad_zone)
         DEALLOCATE(major_rad_zone)
         DEALLOCATE(xtop)
         DEALLOCATE(count_gabage)      
!
!........For special output      
!   
         ALLOCATE(cell_idx_old(nsize_max,2,nset))
         cell_idx_old=cell_idx_tmp
!
!........obtain the cell index and distance factor
!         
         DO iloc=1,nloc_out
            pair=pair_idx(iloc)
            j=set_idx(iloc)
!
!...........obtain d_idx_tmp(distance from xloc_out_tmp)        
            x1(:)=xloc_out_tmp(iloc,:)
            DO i=1,nsize(1,j) !count_rad_zone(pair,j)
               k1=cell_idx_tmp(i,pair,j)             
               d1=0.0d0
               x(:)=xloc_tmp(k1,:)-x1(:)
               DO ix=1,3
                  d1=d1+x(ix)*x(ix)
               ENDDO   
               d_idx_tmp(i,pair,j)=d1
            ENDDO
!
!...........sort according to the distance        
            DO loop=1,10000
               sort_opt=0
               DO i=1,nsize(1,j)-1 !count_rad_zone(pair,j)-1
                  k1=cell_idx_tmp(i,pair,j)
                  k2=cell_idx_tmp(i+1,pair,j)
                  d1=d_idx_tmp(i,pair,j)
                  d2=d_idx_tmp(i+1,pair,j)
                  IF(d1.gt.d2)THEN
                     cell_idx_tmp(i,pair,j)=k2
                     cell_idx_tmp(i+1,pair,j)=k1
                     d_idx_tmp(i,pair,j)=d2
                     d_idx_tmp(i+1,pair,j)=d1
                     sort_opt=1
                  ENDIF   
               ENDDO   
               IF(sort_opt.eq.0)EXIT
            ENDDO   
            IF(sort_opt.eq.1)THEN
               WRITE(*,*)'7error during sorting cell_idx_tmp!!!'
               PAUSE
               STOP
            ENDIF    
!            
!...........obtain the cell index and distance factor            
            cell_out_tmp(iloc,1)=cell_idx_tmp(1,pair,j)  
            cell_out_tmp(iloc,2)=cell_idx_tmp(2,pair,j) 
            nproperty_out_tmp(iloc)=nproperty(pair,j)
            d1=d_idx_tmp(1,pair,j)
            d2=d_idx_tmp(2,pair,j)         
            frac_out_tmp(iloc,1)=d2/(d1+d2)
            frac_out_tmp(iloc,2)=1.0d0-frac_out_tmp(iloc,1)
         ENDDO
!         
         cell_idx_tmp=cell_idx_old
         DEALLOCATE(cell_idx_old)
         DEALLOCATE(d_idx_tmp)
         DEALLOCATE(xloc_out_tmp)
         DEALLOCATE(count_rad_zone)
         DEALLOCATE(set_idx)
         DEALLOCATE(pair_idx)
!         
         DEALLOCATE(xbottom)
      ENDIF !myrank.eq.0         
!   
!.....communicate nsize_max, nsize, cell_idx_tmp
!      
      IF(np.gt.1)THEN
         CALL broadcast_i1(nsize_max)
         DO j=1,nset
            CALL broadcast_i(nsize(:,j),2)
         ENDDO   
      ENDIF   
      IF(myrank.ne.0)THEN
            ALLOCATE(cell_idx_tmp(nsize_max,2,nset))
            cell_idx_tmp(:,:,:)=0
      ENDIF
!      
      END SUBROUTINE radiation_component_init_zone
!
