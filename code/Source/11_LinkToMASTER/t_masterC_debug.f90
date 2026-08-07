!
      SUBROUTINE t_masterC_assembly_debug(i_opt)
!
!     This routine calculates MASTER dll for Assembly-scale
!
      USE Vol_DATA
      USE Solid_DATA
      USE Zzone      , ONLY: ncell_fluid
      USE Zconst2    , ONLY: dt
      USE Zqvol      , ONLY: qvol_liq
      USE Ztimecon   , ONLY: time
      !Add for CUPID-RV/MASTER - jrlee
      USE Zrv_ncell, ONLY:ncell_fuel_rod,cupid_cell_hts2d, &
                          master_to_assem1, p3d_cupid, qvol_mas,master_to_assem1_rod, master_to_assem1_cell
                          
                          
      USE Zrv_hts_2d, ONLY:t_fuel,nr_2d
      USE zcore     , ONLY:myrank,np
      USE zmpi      , ONLY:jperm
      USE Zcoord3   , ONLY:volp      
      !Add linkmaster
      USE MASTER4   , ONLY: i_flag,TTIME,DTTR,TFC,TFS,TCOO,DCOO,BCOO,          &
                            NZ_TH,ppm,ppm_mas,ncb,ncb_mas,cbnam,cbnam_mas,     &
                            zcb,zcb_mas,PWTH,P3D_TH,power_master,VOL_TH,nxy_th,&
                            rv_model_master,pin3d_th
!DEC$IF defined (master_flag)
      USE MASTER4   , ONLY: ppct_master,RSTFN_INP,iok
      USE Zio_unit  , ONLY: unit_log      
!DEC$ENDIF
      !new mapping - yhy
      USE Zrv_ncell , ONLY: asm_ni,asm_nz 
!
      IMPLICIT NONE
      INCLUDE 'master_c.h'
      INCLUDE 'master_dll.h'
!
      INTEGER i,ic,j,k
      INTEGER i_opt
      INTEGER iz,im
      INTEGER ji,m
!
      LOGICAL,SAVE:: old_mapping,init_oldmap   
      DATA old_mapping/.TRUE./
      DATA init_oldmap/.TRUE./
!
      REAL(8) temp_mix
      REAL(8) p3d_cupid_all
!      
!.....time & dt
      i_flag=i_opt
      TTIME=time
      DTTR=dt
!
      TFC(:,:)=0.0d0
      TFS(:,:)=0.0d0
      TCOO(:,:)=0.0d0
      DCOO(:,:)=0.0d0
      BCOO(:,:)=0.0d0
!
!.....Local cell connectivity (ht_str, cell)
!
      IF(old_mapping .and. init_oldmap)then
         init_oldmap=.false.
         ALLOCATE(master_to_assem1_cell(NXY_TH,NZ_TH))
         ALLOCATE(master_to_assem1_rod(NXY_TH,NZ_TH))
         master_to_assem1_cell=0
         master_to_assem1_rod=0
         DO iz=1,NZ_TH
            DO im=1,NXY_TH
               ic=master_to_assem1(im,iz)
               do i=1,ncell_fuel_rod
                  j=cupid_cell_hts2d(i)
                  if(jperm(j).eq.ic) then
                     master_to_assem1_rod(im,iz)=i
                     master_to_assem1_cell(im,iz)=j
                  endif
               enddo
            ENDDO
         ENDDO
      ENDIF !rv_model_master
!      
      IF(rv_model_master.eq.1)then
         IF(old_mapping)THEN      
            DO iz=1,NZ_TH
               DO im=1,NXY_TH
                  i=master_to_assem1_cell(im,iz)
                  j=master_to_assem1_rod(im,iz)
                  IF(j.eq.0.and.i.ne.0) WRITE(*,*)'irod,icel=',j,i
                  IF(j.ne.0.and.i.eq.0) WRITE(*,*)'irod,icel=',j,i
                  IF(j.eq.0.or. i.eq.0) CYCLE
                  temp_mix=(1.0d0-cell%alphag(i))*cell%tl(i)+cell%alphag(i)*cell%tg(i)
                  TCOO(im,iz)=temp_mix-273.15d0
                  DCOO(im,iz)=cell%rhom(i) 
                  TFC(im,iz)=t_fuel(j,    1)-273.15d0
                  TFS(im,iz)=t_fuel(j,nr_2d)-273.15d0
               ENDDO
            ENDDO 
         ELSE         
            DO m=1,ncell_fuel_rod
               i=cupid_cell_hts2d(m)
               ji=jperm(i)
               im=asm_ni(ji)                            
               IF(im.eq.0) CYCLE 
               iz=asm_nz(ji)                     
               TCOO(im,iz)=cell%tl(i)     -273.15d0
               DCOO(im,iz)=cell%rhol(i)
               BCOO(im,iz)=cell%cboron(i) 
               TFC(im,iz)=t_fuel(m,    1)-273.15d0
               TFS(im,iz)=t_fuel(m,nr_2d)-273.15d0
            ENDDO    
         ENDIF  
         IF(np.gt.1)THEN      
            DO i=1,NZ_TH
               CALL allreducei_r(TCOO(1,i),NXY_TH)
               CALL allreducei_r(DCOO(1,i),NXY_TH)
               CALL allreducei_r(BCOO(1,i),NXY_TH)
               CALL allreducei_r(TFC(1,i),NXY_TH)
               CALL allreducei_r(TFS(1,i),NXY_TH)
            ENDDO
         ENDIF
      ELSEIF(rv_model_master.eq.0)then
         DO iz=1,NZ_TH
            DO im=1,NXY_TH
               k=master_to_assem1_cell(im,iz)
               IF(k.ne.0)then
                  temp_mix=(1.0d0-cell%alphag(k))*cell%tl(k) &
                                 +cell%alphag(k) *cell%tg(k)
                  TCOO(im,iz)=temp_mix-273.15d0
                  DCOO(im,iz)=cell%rhol(k) 
                  BCOO(im,iz)=cell%cboron(k)                   
                  TFC(im,iz)= TCOO(im,iz)+300.d0
                  TFS(im,iz)= TCOO(im,iz)+  5.d0
               ENDIF
            ENDDO
         ENDDO
         IF(np.gt.1)THEN      
            DO i=1,NZ_TH
               CALL allreducei_r(TCOO(1,i),NXY_TH)
               CALL allreducei_r(DCOO(1,i),NXY_TH)
               CALL allreducei_r(BCOO(1,i),NXY_TH)
               CALL allreducei_r(TFC(1,i),NXY_TH)
               CALL allreducei_r(TFS(1,i),NXY_TH)
            ENDDO
         ENDIF
      ENDIF !rv_model_master
!      
!.....Reflector zone - Initial Condition is OK
!
      DO iz=1,NZ_TH
         DO im=1,NXY_TH
            IF(TCOO(im,iz).eq.0.d0)THEN
               TCOO(im,iz)=296.d0
               DCOO(im,iz)=735.18d0
               TFC(im,iz) = TCOO(im,iz)+300.d0
               TFS(im,iz) = TCOO(im,iz)+  5.d0
            ENDIF   
         ENDDO
      ENDDO
!
!.....Reactivity, next 
!   
       ppm=ppm_mas
       ncb=ncb_mas !bug
       DO i=1,ncb_mas
          cbnam(i)=cbnam_mas(i)
          zcb(i)=zcb_mas(i)
       ENDDO
!          
      PWTH=power_master
      P3D_TH=0.0d0
      VOL_TH=0.0d0
      PIN3D_TH=0.0d0
!      
      IF(myrank.eq.0)then
!DEC$IF defined (master_flag) 
!
!........Steady state calculation
         IF(i_flag.eq.1)then
            DO WHILE(1)
            CALL MASTER_DLL(I_FLAG, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_INP)
            IF(iok.eq.2)EXIT
            ENDDO
!
!........Transient calculation 
         ELSEIF(i_flag.eq.2)then
            CALL MASTER_DLL(I_FLAG, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_INP)
!
!........WRITE restart file
         !ELSEIF(i_flag.eq.5)then
         !   RSTFN_INP="MAS_RST.STD.BIN"
         !   WRITE(*,*)'Ready to write MAS_RST.STD.BIN'
         !   WRITE(unit_log,*)'Ready to write MAS_RST.STD.BIN'               
         !   CALL MASTER_DLL(5, TTIME, DTTR,                        &
         !                      TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
         !                      PWTH,                                       &
         !                      NXY_TH, NZ_TH,                              &
         !                      PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
         !                      IOK, RSTFN_INP)
         !   WRITE(*,*)'Finish writing MAS_RST.STD.BIN'
         !   WRITE(unit_log,*)'Finish writing MAS_RST.STD.BIN'
         !   RETURN
         ELSE
            WRITE(*,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC, but i_opt= !!!',i_opt
            WRITE(unit_log,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC !!!',i_opt              
            PAUSE
            STOP
         ENDIF
!DEC$ENDIF 
      ENDIF !myrank.eq.0
!
      IF(np.gt.1)then
         DO i=1,NZ_TH
            CALL allreducei_r(P3D_TH(1,i),NXY_TH)
            CALL allreducei_r(vol_TH(1,i),NXY_TH)
         ENDDO
      ENDIF
!
      IF(old_mapping)THEN
         p3d_cupid=0.d0
         p3d_cupid_all=0.0d0
         qvol_mas=0.0d0
         DO iz=1,NZ_TH
            DO im=1,NXY_TH !2x2 subchannel
                i=master_to_assem1_cell(im,iz)
                j=master_to_assem1_rod(im,iz)
                IF(j.eq.0.or.i.eq.0)CYCLE
                p3d_cupid(j)=p3d_cupid(j)+P3D_TH(im,iz)*1.e6*vol_th(im,iz)/1.e6
                p3d_cupid_all=p3d_cupid_all+P3D_TH(im,iz)*1.e6*vol_th(im,iz)/1.e6
                IF(rv_model_master.eq.1) qvol_mas(i)=qvol_mas(i)+P3D_TH(im,iz)*1.e6
                IF(rv_model_master.eq.0) qvol_mas(i)=qvol_mas(i)+P3D_TH(im,iz)*VOL_TH(im,iz)
            ENDDO
         ENDDO
      ELSE
         p3d_cupid=0.0d0
         p3d_cupid_all=0.0d0
         qvol_mas=0.0d0
         DO m=1,ncell_fuel_rod
            i=cupid_cell_hts2d(m)
            ji=jperm(i)
            im=asm_ni(ji)
            IF(im.eq.0) CYCLE 
            iz=asm_nz(ji)
            p3d_cupid(m)=p3d_cupid(m)+P3D_TH(im,iz)*1.e6*vol_th(im,iz)/1.e6
            p3d_cupid_all=p3d_cupid_all+P3D_TH(im,iz)*1.e6*vol_th(im,iz)/1.e6
            qvol_mas(i)=p3d_cupid(m)
         ENDDO  
      ENDIF      
      IF(np.gt.1) CALL allreducei_r1(p3d_cupid_all)
      IF(myrank.eq.0) WRITE(*,"(a,1f10.2,1x,a)")'p3d_cupid_all at t_masterC_assembly_debug=',p3d_cupid_all/1.e6,'MW'     
!
!.....Qvol (W/m3)
!
      IF(rv_model_master.eq.0)then
         qvol_liq=0.0d0
         DO i=1,ncell_fluid
            IF(asm_ni(i).ne.0)THEN
               qvol_liq(i)=qvol_mas(i)/volp(i)
            ENDIF
         ENDDO
      ENDIF
!
      END SUBROUTINE t_masterC_assembly_debug
