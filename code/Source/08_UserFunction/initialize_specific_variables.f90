!
      SUBROUTINE initialize_specific_variables
!
!     Define initial user-defined flow properties (only when "udfl_init_variables" is used.)
!     
      USE VOL_DATA     , ONLY: cell
      USE SOLID_DATA   , ONLY: solid
      USE Zmpi         , ONLY: jperm
      USE Zzone        , ONLY: ncell_fluid,ncell_cond,nmaterial_c,nzone,num_max_zone
      USE Zrv_ncell    , ONLY: ncell_fluid_core,n_channel_fluid,nz_fluid,cupid_cell_channel,num_ch
      USE Zcore        , ONLY: myrank,np
      USE Zparam       , ONLY: ndim
      USE Znum_cell    , ONLY: istart_nf,istart_nb1,icell_nb
      USE Zvec_index    , ONLY: left_nf
      USE Zb_condition , ONLY: v_wall
      USE Zbc_index    , ONLY: npb,ngrad
      USE Zconst1      , ONLY: vv_prob
      USE Zcoord1      , ONLY: xloc,xloc_c
      USE Zpress       , ONLY: p
      USE Ztimecon     , ONLY: alpha_min
      USE Zcoord3      , ONLY: floss
      USE Zmodel       , ONLY: wVertical
      USE Zvector      , ONLY: vg_n,vl_n,vd_n
      
      !OPR1000 rod-scale
      USE MASTER4      , ONLY: nchn,nxyf,nz_th
      USE Zrv_mpi      , ONLY: jperm_fuel_rod
      USE Zzone        , ONLY: num_fzone
      USE Zrv_hts_2d   , ONLY: nz0_2d,nrod_2d,nr_2d,t_fuel
      USE Zrv_ncell    , ONLY: p3d_cupid,qvol_mas,                                 &
                               ncell_fuel_rod,ncell_fuel_rod_all,cupid_cell_hts2d, &
                               asm_ni,asm_nz,chn_nx,chn_ny,                        &
                               nrod_fuel_rod,nz_fuel_rod
      USE Zporous  , ONLY: chn_type,chn_type_tmp
      USE Zcoord2      , ONLY: cell_leng
      USE Zvec_geo     , ONLY: xn_nf
      USE Zgrad_ls_c3d , ONLY: lsindex
      USE MASTER4      , ONLY: npinx
      USE Zio_unit     , ONLY: unit_log
      USE Zgrad_ls_c3d , ONLY: lsindex 
!
      IMPLICIT NONE
!
!.....Local variables
      INTEGER :: i,j,nb
      INTEGER :: ii,k,m
      INTEGER :: nf_number,istart,len,istart1,i1
      INTEGER :: err,ch_opt,nz0_opt
      INTEGER :: iOKr,iOKk
      LOGICAL,SAVE :: initial=.true.
      REAL(8) :: xxn,xzn
      REAL(8) :: CpVol,Condu
!.....Local arrays
      REAL(8) :: v0(num_max_zone,ndim)
      REAL(8) :: zF(7),tF(7),zS(5),tS(5),FACT     ! ST2-CT-01
!.....Local allocatable arrays
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: pi,tli,tgi,agi,qualai,flxi,flyi,flzi

      !OPR1000 rod-scale
      INTEGER :: nz,ji,iz,ia,cx,cy,ni
      INTEGER :: dumi
      INTEGER :: itchn,itrod,iqrod
      REAL(8) :: tave0(num_fzone),pave0(num_fzone)
      REAL(8) :: tchn(nchn,nchn,nxyf,nz_th),rchn(nchn,nchn,nxyf,nz_th)
      REAL(8) :: t_fuel0(ncell_fuel_rod_all,nr_2d)
      REAL(8) :: p3d_cupid0(ncell_fuel_rod_all)
      
      !single_assembly
      INTEGER :: ci,x1,y1,xy
      INTEGER :: cj,sc
      INTEGER :: nrod_fuel_rod0(ncell_fuel_rod_all)
      INTEGER :: guide0(nrod_2d)
      
      !Nuscale
      REAL(8) :: pzr_height,pool_height       
      

!
!.....siphon, initialize volume fraction because tank has water leveln
! 
      IF(vv_prob.eq.'siphon')CALL siphon_initialization     
!
!.....UTPF-RV
!      
      IF(vv_prob.eq.'UTPF-RV')THEN
         ngrad(:)=0
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=xn_nf(i1,3)
               IF(ABS(xxn).gt.0.5d0) ngrad(ii)=1 
!                 IF(xloc(i,2).gt.0.5d0) ngrad(i)=0
            ENDDO               
         ENDDO               
      ENDIF
!
      IF(vv_prob.eq.'Nuscale-03Pool')THEN
         ngrad(:)=0
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=xn_nf(i1,2)  !!2d
               IF(ABS(xxn).gt.0.5d0) ngrad(ii)=1 
            ENDDO               
         ENDDO  
      ENDIF       
!
!.....CUPID-RV
!     
      IF(vv_prob.eq.'fs_31203'.or.vv_prob.eq.'fs_31203_3D'.or.&
         vv_prob.eq.'fs_31302'.or.vv_prob.eq.'fs_31302_3D'.or.&
         vv_prob.eq.'fs_31701'.or.vv_prob.eq.'fs_31701_3D'.or.&
         vv_prob.eq.'fs_31805'.or.vv_prob.eq.'fs_31805_3D'.or.&
         vv_prob.eq.'rbht1196_1d'.or.vv_prob.eq.'rbht1196_3d'.or.&
         vv_prob.eq.'icarus2002') THEN    
         floss(:,:)=0.d0
!
         OPEN(50,file='ht_fluid.in',status='old',iostat=err)
         IF(err.gt.0) STOP '### ht_fluid.in is missing'
!
         i=num_ch
         j=nz0_2d
         ALLOCATE(pi(j,i),tli(j,i),tgi(j,i),agi(j,i),qualai(j,i),flxi(j,i),flyi(j,i),flzi(j,i))
!
         READ(50,*)i,j,ch_opt,nz0_opt
         IF(num_ch.ne.i.or.nz0_2d.ne.j)THEN
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'Error in ht_fluid.in!'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Error in ht_fluid.in!'
         ENDIF
         IF(nz0_opt.eq.1)THEN
            IF(myrank.eq.0)WRITE(*,"(11x,a)")'Reading ht_fluid.in...'
            IF(myrank.eq.0)WRITE(unit_log,"(11x,a)")'Reading ht_fluid.in...'
            IF(ch_opt.eq.1)THEN
               DO i=1,num_ch
                  DO j=1,nz0_2d
                     READ(50,*) pi(j,i),tli(j,i),tgi(j,i),agi(j,i),qualai(j,i),flxi(j,i),flyi(j,i),flzi(j,i)
                  ENDDO
               ENDDO
            ELSE
                DO j=1,nz0_2d
                   READ(50,*) pi(j,1),tli(j,1),tgi(j,1),agi(j,1),qualai(j,1),flxi(j,1),flyi(j,1),flzi(j,1)
                ENDDO
                DO i=2,num_ch
                   DO j=1,nz0_2d
                      pi(j,i)     =pi(j,1)
                      tli(j,i)    =tli(j,1)
                      tgi(j,i)    =tgi(j,1)
                      agi(j,i)    =agi(j,1)
                      qualai(j,i) =qualai(j,1)
                      flxi(j,i)   =flxi(j,1)
                      flyi(j,i)   =flyi(j,1)
                      flzi(j,i)   =flzi(j,1)
                   ENDDO
                ENDDO
            ENDIF
!
            DO i=1,ncell_fluid_core
               j=n_channel_fluid(i)
               k=nz_fluid(i)
               m=cupid_cell_channel(i)
               cell%p(m)     =pi(k,j)
               cell%tl(m)    =tli(k,j)
               cell%tg(m)    =tgi(k,j)
               p(m)=cell%p(m)
               cell%alphag(m)=agi(k,j)
               cell%alphal(m)=1.d0-agi(k,j)
               cell%alphad(m)=0.d0
               cell%quala(m) =qualai(k,j)
               floss(m,1)    =flxi(k,j)
               floss(m,2)    =flyi(k,j)
               IF(ndim.eq.3) floss(m,3)=flzi(k,j)
            ENDDO
         ENDIF
!
         DEALLOCATE(pi,tli,tgi,agi,qualai,flxi,flyi,flzi)
!
      ENDIF               
!   
!.....PAFS-POOL
!      
      IF(vv_prob.eq.'PAFS-POOL')THEN
         DO i=1,ncell_fluid
            IF(xloc(i,2).gt.9.84d0)THEN
               cell%alphal(i)=0.001d0
               cell%alphag(i)=0.999d0-alpha_min
               cell%alphad(i)=alpha_min
               cell%quala(i)=0.99d0
               cell%tl(i)=313.15d0
               cell%tg(i)=313.15d0
               IF(npb(i).eq.1)THEN
                  cell%alphal(i)=alpha_min
                  cell%alphag(i)=1.d0-alpha_min-alpha_min
                  cell%alphad(i)=alpha_min
                  cell%quala(i)=0.99d0
                  cell%tl(i)=313.15d0
                  cell%tg(i)=313.15d0
               ENDIF
            ELSE
               cell%alphal(i)=0.999d0-alpha_min
               cell%alphag(i)=0.001d0
               cell%alphad(i)=alpha_min
               cell%quala(i)=0.99d0   
               cell%tl(i)=313.15d0
               cell%tg(i)=313.15d0 
            ENDIF
         ENDDO  
         ngrad(:)=0
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=xn_nf(i1,2)
               IF(ABS(xxn).gt.0.5d0) ngrad(ii)=1 
!                 IF(xloc(i,2).gt.0.5d0) ngrad(i)=0
            ENDDO               
         ENDDO               
      ENDIF
      
!     
!......Nuscale
!      
      IF(vv_prob.eq.'Nuscale-RRV')THEN
         DO i=1,ncell_fluid
            IF(nzone(i).ne.2) CYCLE
            IF(xloc(i,3).gt.8.0d0)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0 
               p(i)=1.d5+996.66*9.8*(8.00-xloc(i,ndim))
               cell%p(i)=p(i)
            ENDIF
         ENDDO  

         DO i=1,ncell_fluid
            IF(nzone(i).ne.3) CYCLE 
            IF(xloc(i,3).gt.18.0d0)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0 
               p(i)=1.d5+996.66*9.8*(18.00-xloc(i,ndim))
            ENDIF
         ENDDO           
!         
!        ls off for vertial wall
         lsindex(:)=0  !all cells are basically LS off.
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=DABS(xn_nf(i1,ndim))
               IF(xzn.gt.0.1d0.or. xxn.gt.0.1d0) lsindex(ii)=1   !ls is on for whole wall face
               IF(xzn.lt.0.1d0.and.xxn.gt.0.1d0) lsindex(ii)=0  !ls is off for vertical wall
            ENDDO               
         ENDDO  
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=DABS(xn_nf(i1,ndim))
               IF(xzn.gt.0.1d0) lsindex(ii)=1  !ls is on for corner cell.
            ENDDO               
         ENDDO           
      ENDIF
!      
      IF(vv_prob.eq.'Nuscale-RVV')THEN
         DO i=1,ncell_fluid
            IF(nzone(i).ne.1) CYCLE 
            pzr_height=14.8d0  !13.9d0 
            IF(xloc(i,3).gt.pzr_height)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.0d0
               cell%tl(i)=620.00d0
               cell%tg(i)=620.00d0
            ELSE
               cell%alphal(i)=0.999d0
               cell%alphag(i)=0.001d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0 
               cell%tl(i)=568.d0 
               cell%tg(i)=568.d0 
               p(i)=150.d5+735.9*9.8*(pzr_height-xloc(i,ndim))  
            ENDIF
         ENDDO  
!
         DO i=1,ncell_fluid
            IF(nzone(i).ne.3) CYCLE 
            pool_height=20.1d0 !17.2d0  
            IF(xloc(i,3).gt.pool_height)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=1.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0 
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0 
               p(i)=1.d5+996.66*9.8*(pool_height-xloc(i,ndim))
            ENDIF
         ENDDO           
!         
!        ls off for vertial wall
         lsindex(:)=0  !all cells are basically LS off.
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xzn=DABS(xn_nf(i1,ndim))
               IF(xzn.gt.0.1d0) lsindex(ii)=1  !ls is on for corner cell.
            ENDDO               
         ENDDO   
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=xn_nf(i1,ndim) 
               IF(DABS(xzn).gt.0.1d0.and.xxn.gt.0.1d0) lsindex(ii)=0  !ls is off for vertical wall
            ENDDO               
         ENDDO          
      ENDIF   
!      
      IF(vv_prob.eq.'Nuscale-PZR')THEN
         DO i=1,ncell_fluid
            IF(nzone(i).ne.1) CYCLE 
            pzr_height=13.9d0  !14.8d0 !14.d0 
            IF(xloc(i,3).gt.pzr_height)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.0d0
               cell%tl(i)=620.00d0
               cell%tg(i)=620.00d0
            ELSEIF(xloc(i,3).gt.13.8.and.xloc(i,3).gt.13.9) THEN
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=610.d0
               cell%tg(i)=610.d0
               p(i)=150.d5+735.9*9.8*(pzr_height-xloc(i,ndim))                  
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=568.d0
               cell%tg(i)=568.d0
               p(i)=150.d5+735.9*9.8*(pzr_height-xloc(i,ndim))  
            ENDIF
         ENDDO  
!
         DO i=1,ncell_fluid
            IF(nzone(i).ne.3) CYCLE 
            pool_height=20.1d0 !17.2d0  
            IF(xloc(i,3).gt.pool_height)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=1.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0
               cell%tl(i)=300.00d0
               cell%tg(i)=300.00d0 
               p(i)=1.d5+996.66*9.8*(pool_height-xloc(i,ndim))
            ENDIF
         ENDDO           
!               
      ENDIF
!      
      IF(vv_prob.eq.'KSMR-PZR')THEN
         DO i=1,ncell_fluid
            pzr_height=16.0d0 
            IF(xloc(i,3).gt.pzr_height)THEN
               cell%alphal(i)=0.d0
               cell%alphag(i)=1.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.0d0
               cell%tl(i)=615.50d0
               cell%tg(i)=615.50d0
            ELSEIF(xloc(i,3).gt.15.8.and.xloc(i,3).le.pzr_height) THEN
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0 
               cell%tl(i)=600.0d0 
               cell%tg(i)=600.0d0 
               p(i)=150.d5+886.15*9.81*(pzr_height-xloc(i,ndim))                     
            ELSE
               cell%alphal(i)=1.d0
               cell%alphag(i)=0.d0
               cell%alphad(i)=0.d0
               cell%quala(i)=0.d0 
               cell%tl(i)=568.d0 
               cell%tg(i)=568.d0 
               p(i)=150.d5+736.24*9.81*(pzr_height-xloc(i,ndim))  
            ENDIF
         ENDDO  
!         
!        ls off for vertial wall
         lsindex(:)=0  !all cells are basically LS off.
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xzn=DABS(xn_nf(i1,ndim))
               IF(xzn.gt.0.1d0) lsindex(ii)=1  !ls is on for corner cell.
            ENDDO               
         ENDDO   
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=xn_nf(i1,ndim) 
               IF(DABS(xzn).gt.0.1d0.and.xxn.gt.0.1d0) lsindex(ii)=0  !ls is off for vertical wall
            ENDDO               
         ENDDO   
!         
      ENDIF          
!   
!.....manometric
! 	  
      IF (vv_prob.eq.'manometric')THEN
         DO i=1,ncell_fluid
            cell%alphag(i)=0.99d0
            cell%alphad(i)=0.d0
            cell%quala(i)=1.d0
            IF(xloc(i,2).le.0.3d0)THEN
               cell%alphag(i)=0.01d0
            ELSEIF(xloc(i,2).le.0.7d0.and.xloc(i,1).lt.0.d0)THEN
               cell%alphag(i)=0.01d0
            ENDIF
            cell%alphal(i)=1.d0-cell%alphag(i)
         ENDDO
      ENDIF
!   
!.....moving_wall
! 
      IF (vv_prob.eq.'moving_wall')THEN
         v_wall(1)=0.d0
         v_wall(2)=-0.1d0
      ENDIF
!   
!.....manoair, manoair_mcc
! 
      IF(vv_prob.eq.'manoair'.or.vv_prob.eq.'manoair_mcc')THEN
         DO i=1,ncell_fluid
            cell%alphag(i)=0.99d0
            cell%alphad(i)=0.d0
            cell%quala(i)=1.d0
            IF(xloc(i,2).le.0.3d0)THEN
               cell%alphag(i)=0.01d0
            ELSEIF (xloc(i,2).le.0.7d0.and.xloc(i,1).lt.0.d0)THEN
               cell%alphag(i)=0.01d0
            ENDIF
            cell%alphal(i)=1.d0-cell%alphag(i)
         ENDDO
!   
!.....manosteam, manosteam_mcc
!          
      ELSEIF(vv_prob.eq.'manosteam'.or.vv_prob.eq.'manosteam_mcc')THEN
         DO i=1,ncell_fluid
            cell%alphag(i)=0.99d0
            cell%alphad(i)=0.d0
            IF(xloc(i,2).le.0.3d0)THEN
               cell%alphag(i)=0.01d0
            ELSEIF(xloc(i,2).le.0.7d0.and.xloc(i,1).lt.0.d0)THEN
               cell%alphag(i)=0.01d0
            ENDIF
            cell%alphal(i)=1.d0-cell%alphag(i)
         ENDDO
      ENDIF                
!   
!.....fluidic_device
!       
      IF(vv_prob.eq.'fluidic_device')THEN
!
!.........Low pressure case:   water volume=51.94, water level=51.94/5.91=8.788
!.........High pressure case1: measured water level=8.91 --> water volume=49.979, water level=49.979/5.91=8.45669 
!
          DO i=1,ncell_fluid
             IF(xloc(i,ndim).gt.8.49d0)THEN
                cell%alphal(i)=0.00001d0
                cell%alphag(i)=0.99999d0
                cell%alphad(i)=0.d0
                p(i)=2.16d6
             ELSEIF(xloc(i,ndim).lt.8.26d0)THEN
                cell%alphal(i)=0.99999d0
                cell%alphag(i)=0.00001d0
                cell%alphad(i)=0.d0
                p(i)=2.16d6+998.96*9.8*(8.26-xloc(i,ndim))
             ELSE
                cell%alphal(i)=(49.979-5.91*8.26)/(8.49-8.26)/5.91
                cell%alphag(i)=1.d0-cell%alphal(i)
                cell%alphad(i)=alpha_min
                p(i)=2.16d6
             ENDIF
!
             cell%tl(i)=293.15d0
             cell%tg(i)=293.15d0
             cell%quala(i)=1.d0
         ENDDO
!
      ENDIF
!
!.....Adward_pipe
!
      IF(vv_prob.eq.'T_blowdown')THEN
!
          DO i=1,ncell_fluid
             IF(xloc(i,1).lt.4.d0)THEN
                cell%alphal(i)=1.d0
                p(i)=0.2d6
                cell%tl(i)=390.d0
                cell%tg(i)=393.d0   ! T_sat at 0.2MPa
!
             ELSE
                cell%alphal(i)=0.d0
                p(i)=0.1d6
                cell%tl(i)=373.15d0
                cell%tg(i)=373.15d0   ! T_sat at 0.1MPa
             ENDIF
             cell%alphag(i)=1.d0-cell%alphal(i)
             cell%alphad(i)=0.d0
!
             cell%quala(i)=0.d0
         ENDDO
!
      ENDIF
!      
      DO i=1,ncell_cond
         IF(ABS(nmaterial_c(i)).lt.50)THEN
            CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iOKr,iOKk)
            solid%rhocps(i)=CpVol
            solid%conds(i)=Condu
         ENDIF
         solid%matnum(i)=ABS(nmaterial_c(i))         
      ENDDO
!
!.....rocom, rocom_mc
!
      IF(vv_prob.eq.'rocom'.or.vv_prob.eq.'rocom_mc') CALL rocom_porous_ring_user
!
!
      IF(vv_prob.eq.'duct')THEN
!         IF(myrank.eq.0)WRITE(97,*)'Initial variables were assigned!!!'
          DO i=1,ncell_fluid
             IF(xloc(i,1).ge.0.d0)THEN
                cell%tl(i)=400.d0
                cell%tg(i)=400.d0
             ELSE
                cell%tl(i)=300.d0
                cell%tg(i)=300.d0
             ENDIF
         ENDDO
      ENDIF   
!
!.....sg_separator
!
      IF(vv_prob.eq.'sgp_separator')THEN
          DO i=1,ncell_fluid
             IF(xloc(i,3).gt.12.5442d0)THEN
                cell%alphal(i)=0.5d0
                cell%tl(i)=550.94
                cell%tg(i)=558.94
!
             ELSE
!                cell%alphal(i)=1.d0
                cell%alphal(i)=0.5d0
                cell%tl(i)=550.0
                cell%tg(i)=558.94
             ENDIF
             cell%alphag(i)=1.d0-cell%alphal(i)
             cell%alphad(i)=0.d0
             cell%quala(i)=0.d0
         ENDDO
      ENDIF
!               
      IF(vv_prob.eq.'kuhn_111'.or.vv_prob.eq.'Nuscale-02'.or.vv_prob.eq.'Nuscale-03Pool')THEN
         ALLOCATE(wVertical(ncell_fluid))
         wVertical(:)=0
         DO nf_number=5,6
            istart1=istart_nb1(1,nf_number)
            len    =istart_nb1(2,nf_number)
            DO nb=1,len
               i1=istart1+nb
               ii=icell_nb(i1)
                wVertical(ii)=1
            ENDDO               
         ENDDO
      ENDIF
!               
      IF(vv_prob.eq.'atlas_mc_porous'.or.vv_prob.eq.'pwr_mc_poro'.or.vv_prob.eq.'apr1400_mc_poro'.or.vv_prob.eq.'opr1000_mc_poro')THEN
         CALL udfn_mat_prop
      ENDIF        
!           
      IF(vv_prob.eq.'vessel_test')THEN
      
         DO i=1,ncell_fluid
            IF(nzone(i).eq.4) THEN
                IF(xloc(i,1).gt.0.d0) THEN
                   vl_n(i,1)=v0(nzone(i),1)
                   vd_n(i,1)=v0(nzone(i),1)
                   vg_n(i,1)=v0(nzone(i),1)
                ELSE
                   vl_n(i,1)=-v0(nzone(i),1)
                   vd_n(i,1)=-v0(nzone(i),1)
                   vg_n(i,1)=-v0(nzone(i),1)
                ENDIF
            ELSEIF(nzone(i).eq.15)THEN
                IF(xloc(i,1).gt.0.d0) THEN
                   vl_n(i,1)=-v0(nzone(i),1)
                   vd_n(i,1)=-v0(nzone(i),1)
                   vg_n(i,1)=-v0(nzone(i),1)
                ELSE
                   vl_n(i,1)= v0(nzone(i),1)
                   vd_n(i,1)= v0(nzone(i),1)
                   vg_n(i,1)= v0(nzone(i),1)
                ENDIF
                IF(xloc(i,2).gt.0.d0) THEN
                   vl_n(i,2)=-v0(nzone(i),2)
                   vd_n(i,2)=-v0(nzone(i),2)
                   vg_n(i,2)=-v0(nzone(i),2)
                ELSE
                   vl_n(i,2)= v0(nzone(i),2)
                   vd_n(i,2)= v0(nzone(i),2)
                   vg_n(i,2)= v0(nzone(i),2)
                ENDIF                
                   vl_n(i,3)= v0(nzone(i),3)
                   vd_n(i,3)= v0(nzone(i),3)
                   vg_n(i,3)= v0(nzone(i),3)
            ENDIF
         ENDDO
      ENDIF
!
      IF(vv_prob.eq.'h2p1_0'.or.vv_prob.eq.'h2p1_0x'.or. &
         vv_prob.eq.'h2p1_1'.or.vv_prob.eq.'h2p1_1x'.or. &
         vv_prob.eq.'h2p1_2'.or.vv_prob.eq.'h2p1_2x'.or. &
         vv_prob.eq.'h2p1_3'.or.vv_prob.eq.'h2p1_3x'.or. &
         vv_prob.eq.'h2p1_4'.or.vv_prob.eq.'h2p1_4x'.or. &
         vv_prob.eq.'VD_h2p1_0')THEN
         CALL udfn_H2P1_initialHeTg
      ENDIF
!
      IF(vv_prob.eq.'ST2-CT-01'.or.vv_prob.eq.'ST2-CT-02'.or. &
         vv_prob.eq.'ST2-CT-03')THEN
!........Initial water level
         DO i=1,ncell_fluid
            IF(nzone(i).eq.2)THEN      !HVT(z0=2.287m / LT=1.54m???)
               IF(xloc(i,3).le.2.8d0)THEN      !!!
                  cell%alphal(i)=0.999d0
                  cell%alphag(i)=1.0d0-cell%alphal(i)
                  cell%quala(i)=0.999 !1.0d0
               ENDIF
            ELSEIF(nzone(i).eq.3)THEN  !IRWST(z0=2.40752m / LT=1.50m???)
               IF(xloc(i,3).le.2.8d0)THEN     !!!
                  cell%alphal(i)=0.999d0
                  cell%alphag(i)=1.0d0-cell%alphal(i)
                  cell%quala(i)=0.999 !1.0d0
               ENDIF
            ENDIF
         ENDDO
!
         zF=(/4.98d0, 5.55d0, 7.15d0, 8.9592d0, 10.9592d0, 12.95936d0, 13.9592d0/)
         zS=(/4.7d0, 4.98d0, 5.55d0, 7.15d0, 7.45d0/)
         IF(vv_prob.eq.'ST2-CT-01')THEN
            tF=(/28.022493d0, 24.88155d0, 39.702728d0, 48.519222d0, &
                51.710503d0, 52.144196d0, 55.87598d0/)
!            tS=(/22.374044d0, 22.315708d0, 22.439453d0, 25.349928d0, 27.315607d0/)
            tS=(/22.9003246d0, 24.7473791d0, 23.3261695d0, 26.5761515d0, 27.7142368d0/)
         ELSEIF(vv_prob.eq.'ST2-CT-02')THEN
            tF=(/28.864323d0, 27.806581d0, 38.165958d0, 43.070412d0, &
                44.13459d0, 44.223705d0, 43.486626d0/)
!            tS=(/28.21413d0, 24.237345d0, 25.639997d0, 24.803741d0, 25.105078d0/)
            tS=(/27.1354862d0, 25.2684085d0, 26.0969843d0, 26.2237241d0, 25.6603546d0/)
         ELSEIF(vv_prob.eq.'ST2-CT-03')THEN
            tF=(/27.168901d0, 26.450195d0, 36.910694d0, 43.3424d0, &
                44.759022d0, 44.613178d0, 43.964611d0/)
!            tS=(/28.427378d0, 22.373602d0, 23.655453d0, 22.338736d0, 22.823328d0/)
            tS=(/26.3002254d0, 23.4418441d0, 24.1281055d0, 24.1313514d0, 23.4716748d0/)
         ENDIF
!........Initial gas temperature : Lv07 of I / Lv06 to Lv01 of A
         DO i=1,ncell_fluid
            IF(xloc(i,3).le.zF(1))THEN
               cell%tg(i)=tF(1)
            ELSEIF(xloc(i,3).gt.zF(7))THEN
               cell%tg(i)=tF(7)
            ELSE
               DO j=1,6
                  IF(xloc(i,3).gt.zF(j).and.xloc(i,3).le.zF(j+1))THEN
                     FACT=(xloc(i,3)-zF(j))/(zF(j+1)-zF(j))
                     cell%tg(i)=tF(j)+FACT*(tF(j+1)-tF(j))
                  ENDIF
               ENDDO
            ENDIF
            cell%tg(i)=cell%tg(i)+273.15d0
         ENDDO
!........Initial solid temperature : PSW-B / SSW01-B / RFP01-B / SG01-B / PZR-B 
         DO i=1,ncell_cond
         IF(xloc_c(i,3).le.zS(1))THEN
            solid%tsol(i)=tS(1)
         ELSEIF(xloc_c(i,3).gt.zS(5))THEN
            solid%tsol(i)=tS(5)
         ELSE
            DO j=1,4
               IF(xloc_c(i,3).gt.zS(j).and.xloc_c(i,3).le.zS(j+1))THEN
                  FACT=(xloc_c(i,3)-zS(j))/(zS(j+1)-zS(j))
                  solid%tsol(i)=tS(j)+FACT*(tS(j+1)-tS(j))
               ENDIF
            ENDDO
         ENDIF
         solid%tsol(i)=solid%tsol(i)+273.15d0
         IF(ABS(nmaterial_c(i)).lt.50) THEN
            CALL mat_prop(ABS(nmaterial_c(i)),solid%tsol(i),CpVol,Condu,iokr,iokk)
            solid%rhocps(i)=CpVol
            solid%conds(i)=Condu
         ENDIF
         solid%matnum(i)=ABS(nmaterial_c(i))         
         ENDDO
         ENDIF
!
!.....OPR1000 rod-scale FA ONLY
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv')THEN

         IF(initial)then
            IF(myrank.eq.0)then
               write(*,*)'          READ Initial Condition for fuel rod'
               write(*,*)'          trod,qrod,tchn'
            ENDIF
            initial=.false.
         ENDIF

         !Core region, channel-wise temperature
         OPEN(11001,file='tchn_core.in',status='unknown',iostat=itchn)
         IF(itchn.ne.0)then
            IF(myrank.eq.0) then
               WRITE(* ,*)'Program was terminated due to lack of <tchn_core.in>!'
               WRITE(unit_log,*)'Program was terminated due to lack of <tchn_core.in>!'
               CALL finalize_mpi
               STOP
            ENDIF
         ENDIF
         tchn=0.d0
         rchn=0.d0
         DO ia=1,nxyf
            DO cy=1,nchn
               DO cx=1,nchn
                  READ(11001,48) dumi,dumi,dumi, &
                                (tchn(cx,cy,ia,k),k=1,nz_th), &
                                (rchn(cx,cy,ia,k),k=1,nz_th)
               ENDDO
            ENDDO
         ENDDO
48   format(3(i4,1x),100(f15.8,1x))

         DO i=1,ncell_fluid
            !IF(chn_type(i).ne.0)then
               ji=jperm(i)
               iz=asm_nz(ji)
               ia=asm_ni(ji)
               cy=chn_ny(ji)
               cx=chn_nx(ji)
               IF(iz.lt.1 .and. iz.gt.nz_th)STOP 'Error asm_nz for tchn'
               IF(ia.lt.1 .and. ia.gt.nxyf) STOP 'Error asm_na for tchn'
               IF(cx.lt.1 .and. cx.gt.nchn) STOP 'Error chn_nx for tchn'
               IF(cy.lt.1 .and. cy.gt.nchn) STOP 'Error chn_ny for tchn'
               cell%tl(i)  =tchn(cx,cy,ia,iz)
               cell%rhol(i)=rchn(cx,cy,ia,iz)
            !ENDIF
         ENDDO

         !Fuel temperature
         OPEN(11002,file='trod_core.in',status='unknown',iostat=itrod)
         IF(itrod.ne.0)then
            IF(myrank.eq.0) then
               WRITE(* ,*)'Program was terminated due to lack of <trod_core.in>!'
               WRITE(unit_log,*)'Program was terminated due to lack of <trod_core.in>!'
               CALL finalize_mpi
               STOP
            ENDIF
         ENDIF

         t_fuel0=0.d0
         DO i=1,ncell_fuel_rod_all
            READ(11002,49)(t_fuel0(i,k),k=1,nr_2d)
         ENDDO
49   format(30(f15.8,1x))

         DO i=1,ncell_fuel_rod
            ji=jperm_fuel_rod(i)
            t_fuel(i,:)=t_fuel0(ji,:)
         ENDDO

         !Fuel power
         OPEN(11003,file='qrod_core.in',status='unknown',iostat=iqrod)
         IF(iqrod.ne.0)then
            IF(myrank.eq.0) then
               WRITE(* ,*)'Program was terminated due to lack of <qrod_core.in>!'
               WRITE(unit_log,*)'Program was terminated due to lack of <qrod_core.in>!'
               CALL finalize_mpi
               STOP
            ENDIF
         ENDIF

         p3d_cupid0=0.d0
         DO i=1,ncell_fuel_rod_all
            READ(11003,*)p3d_cupid0(i)
         ENDDO

         DO i=1,ncell_fuel_rod
            ji=jperm_fuel_rod(i)
            p3d_cupid(i)=p3d_cupid0(ji)
            k=cupid_cell_hts2d(i)
            qvol_mas(k)=p3d_cupid(i)/cell_leng(k,3)
         ENDDO
         
      ENDIF
!
!.....OPR1000 rod-scale
!
      IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')THEN

!........Global region, nzone(i)
         OPEN(11000,file='tave_vssl.in',status='unknown')
         DO i=1,num_fzone
            READ(11000,*)k,m,tave0(i),pave0(i)
         ENDDO

         DO i=1,ncell_fluid
            nz=nzone(i)
            cell%tl(i)=tave0(nz)
         ENDDO

!........Core region, channel-wise temperature
         OPEN(11001,file='tchn_core.in',status='unknown')
         tchn=0.d0
         rchn=0.d0
         DO ia=1,nxyf
            DO cy=1,nchn
               DO cx=1,nchn
                  READ(11001,48) dumi,dumi,dumi, &
                                (tchn(cx,cy,ia,k),k=1,nz_th), &
                                (rchn(cx,cy,ia,k),k=1,nz_th)
               ENDDO
            ENDDO
         ENDDO

         DO i=1,ncell_fluid
            IF(chn_type(i).ne.0)then
               ji=jperm(i)
               iz=asm_nz(ji)
               ia=asm_ni(ji)
               cy=chn_ny(ji)
               cx=chn_nx(ji)
               IF(iz.lt.1 .and. iz.gt.nz_th)STOP 'Error asm_nz for tchn'
               IF(ia.lt.1 .and. ia.gt.nxyf) STOP 'Error asm_na for tchn'
               IF(cx.lt.1 .and. cx.gt.nchn) STOP 'Error chn_nx for tchn'
               IF(cy.lt.1 .and. cy.gt.nchn) STOP 'Error chn_ny for tchn'
               cell%tl(i)  =tchn(cx,cy,ia,iz)
               cell%rhol(i)=rchn(cx,cy,ia,iz)
            ENDIF
         ENDDO

!........Fuel temperature
         OPEN(11002,file='trod_core.in',status='unknown')
         t_fuel0=0.d0
         DO i=1,ncell_fuel_rod_all
            READ(11002,49)(t_fuel0(i,k),k=1,nr_2d)
         ENDDO

         DO i=1,ncell_fuel_rod
            ji=jperm_fuel_rod(i)
            t_fuel(i,:)=t_fuel0(ji,:)
         ENDDO
      ENDIF
!
!.....OPR1000 rod-scale Single assembly
!
      IF(vv_prob.eq.'OPR1000_single_assem')THEN
!
!........Zero power setting at guide tube region
!
         p3d_cupid =1.d0
            
         CALL allgatherv_i(nrod_fuel_rod,nrod_fuel_rod0,ncell_fuel_rod,ncell_fuel_rod_all,3)
         
         guide0=0
         DO i=1,ncell_fuel_rod
            if(nz_fuel_rod(i).eq.1)then
               ji=jperm_fuel_rod(i)
               ni=nrod_fuel_rod0(ji)
               ci=cupid_cell_hts2d(i)
               cj=jperm(ci)
               sc=chn_type_tmp(cj)
               IF(sc.eq.4)then
                  p3d_cupid0(ji)=0.d0
                  guide0(ni)=1
                  x1=ni-1
                  guide0(x1)=1
                  y1=ni-npinx
                  guide0(y1)=1
                  xy=y1-1
                  guide0(xy)=1
               ENDIF
            endif
         ENDDO
         
         IF(np.gt.1) CALL allreducei_i(guide0,nrod_2d)
          
         DO i=1,ncell_fuel_rod
            ni=nrod_fuel_rod(i)
            IF(guide0(ni).ne.0) p3d_cupid(i)=0.d0
            nz=nz_fuel_rod(i)
            IF(nz.eq.nz0_2d) p3d_cupid(i)=0.d0
         ENDDO
         
         DO i=1,ncell_fuel_rod
            ci=cupid_cell_hts2d(i)
            qvol_mas(ci)=p3d_cupid(i)
         ENDDO

      ENDIF

      IF(vv_prob.eq.'ismr_rad')THEN
         lsindex(:)=0  !all cells are basically LS off.
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               !xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=DABS(xn_nf(i1,ndim))
               IF(xzn.gt.0.1d0) lsindex(ii)=1   !ls is on for whole wall face
            ENDDO               
         ENDDO           
         
         !DO i=1,ncell_fluid
         !   IF(nzone(i).ne.1)CYCLE
         !   IF(xloc(i,3).gt.40d0) THEN
         !      cell%alphal(i)=0.0d0 !0.001d0
         !      cell%alphag(i)=1.d0 !0.999d0-1.d-8
         !      cell%alphad(i)=0.d0 !1.d-8
         !      cell%quala(i)=0.d0
         !      cell%tg(i)=373.15
         !      cell%tl(i)=300.0
         !   ELSE
         !      cell%alphal(i)=1.d0 !0.999d0-1.d-8
         !      cell%alphag(i)=0.d0 !0.001d0
         !      cell%alphad(i)=0.d0 !1.d-8
         !      cell%quala(i)=0.d0
         !      cell%tg(i)=373.15
         !      cell%tl(i)=300.0
         !   ENDIF
         !ENDDO    
         
         DO i=1,ncell_fluid
            IF(xloc(i,3).le.40d0.and.nzone(i).eq.1) THEN
               p(i)=101300.d0+998.96*9.8*(40-xloc(i,ndim))     
            ENDIF
         ENDDO

      ENDIF          

      IF(vv_prob.eq.'ismr_2d')THEN
         lsindex(:)=0  !all cells are basically LS off.
         DO nf_number=1,8
            istart=istart_nf(1,nf_number)
            len   =istart_nf(2,nf_number)
            DO i=1,len
               i1=istart+i
               ii=left_nf(i1)
               !xxn=DSQRT(xn_nf(i1,1)*xn_nf(i1,1)+xn_nf(i1,2)*xn_nf(i1,2))
               xzn=DABS(xn_nf(i1,2))
               IF(xzn.gt.0.1d0) lsindex(ii)=1   !ls is on for whole wall face
            ENDDO               
         ENDDO           
         
         DO i=1,ncell_fluid
            IF(nzone(i).ne.2)CYCLE
            IF(xloc(i,2).gt.20d0) THEN
               cell%alphal(i)=0.0d0 !0.001d0
               cell%alphag(i)=1.d0 !0.999d0-1.d-8
               cell%alphad(i)=0.d0 !1.d-8
               cell%quala(i)=1.0d0
               cell%tg(i)=300.0
               cell%tl(i)=300.0
            ELSE
               cell%alphal(i)=1.d0 !0.999d0-1.d-8
               cell%alphag(i)=0.d0 !0.001d0
               cell%alphad(i)=0.d0 !1.d-8
               cell%quala(i)=1.d0
               cell%tg(i)=300.0
               cell%tl(i)=300.0
            ENDIF
         ENDDO    
         
         DO i=1,ncell_fluid
            IF(xloc(i,2).le.20d0.and.nzone(i).eq.2) THEN
               p(i)=101300.d0+998.96*9.8*(20.d0-xloc(i,2))     
            ENDIF
         ENDDO

      ENDIF  
      
!
      END SUBROUTINE initialize_specific_variables
!
      

