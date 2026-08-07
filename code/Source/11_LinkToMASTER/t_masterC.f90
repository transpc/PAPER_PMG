!
      SUBROUTINE t_masterC(i_opt)
!
      USE Zconst1          , ONLY: cplmaster
      USE Zrv_model        , ONLY: rv_model 
      USE MASTER4          , ONLY: rv_model_master
      USE Zuserdefined     , ONLY: user_iary
      USE Zconst1          , ONLY: restart   
      USE Zcore            , ONLY: myrank
!        
      IMPLICIT NONE
!
      INTEGER i_opt
!
      IF(restart.ne.0.and.user_iary(32).gt.0)i_opt=user_iary(32)
      rv_model_master=rv_model
!  
!.....assembly scale    
!
      IF(cplmaster.eq.1)then 
         CALL averaging_cboron
         CALL CM_Control_Reactivity        
         IF(i_opt.ne.0)then
            CALL t_masterC_assembly(i_opt)
         ELSE
            IF(myrank.eq.0)WRITE(*,"(a)")'i_opt is ZERO in t_masterC in t_masterC!!!' 
            PAUSE
            STOP 
         ENDIF
!
!.....subchannel-scale
!
      ELSEIF(cplmaster.eq.2)then
         CALL averaging_cboron_rod
         CALL CM_Control_Reactivity        
         IF(i_opt.ne.0)then
            CALL t_masterC_rod(i_opt)
         ELSE
            STOP 'i_opt is ZERO in t_masterC'
         ENDIF
      ELSE
        WRITE(*,"(11x,a)")'cplmaster should be 1 or 2!!!'
        PAUSE
        STOP   
      ENDIF
!
      END SUBROUTINE t_masterC
!---------------------------------------------------------------------------------------
      SUBROUTINE t_masterC_assembly(i_opt)
!
!     This routine calculates MASTER dll for Assembly-scale
!
      USE Vol_DATA   , ONLY: cell
      USE Zzone      , ONLY: ncell_fluid,nzone
      USE Zconst2    , ONLY: dt
      USE Zqvol      , ONLY: qvol_liq
      USE Ztimecon   , ONLY: time
      !Add for CUPID-RV/MASTER - jrlee
      USE Zrv_ncell  , ONLY:ncell_fuel_rod,cupid_cell_hts2d,        &
                            master_to_assem1_cell, &
                            p3d_cupid, qvol_mas
      USE Zrv_hts_2d, ONLY: t_fuel,nr_2d
      USE zcore     , ONLY: myrank,np
      USE zmpi      , ONLY: jperm
      USE Zcoord2   , ONLY: cell_leng      
      USE Zcoord3   , ONLY: volp      
      !Add linkmaster
      USE MASTER4   , ONLY: i_flag,TTIME,DTTR,TFC,TFS,TCOO,DCOO,BCOO,          &
                            NZ_TH,ppm,ppm_mas,ncb,ncb_mas,cbnam,cbnam_mas,     &
                            zcb,zcb_mas,PWTH,P3D_TH,power_master,VOL_TH,nxy_th,&
                            rv_model_master,pin3d_th
!DEC$IF defined (master_flag) 
      USE MASTER4   , ONLY: ppct_master,RSTFN_INP,RSTFN_REV,iok
      USE Zio_unit  , ONLY: unit_log      
!DEC$ENDIF 
      !new mapping - yhy
      USE Zrv_ncell , ONLY: asm_ni,asm_ni2,asm_nz,chn_nx,chn_ny
      USE MASTER4   , ONLY: nxyf,npinx,nchn,num_asm,mst_to_asmi
!
      IMPLICIT NONE
      INCLUDE 'master_c.h'
      INCLUDE 'master_dll.h'
!
      INTEGER i
      INTEGER i_opt
!
      REAL(8) temp_mix

!.....CUPID-RV-MASETER
      INTEGER k
      INTEGER iz,im
      REAL(8) p3d_cupid_all
      REAL(8) vol_tot
!.....new mapping - yhy
      REAL(8) tm_cpd_pin(npinx,npinx,nxyf,nz_th)
      REAL(8) rm_cpd_pin(npinx,npinx,nxyf,nz_th)
      REAL(8) cb_cpd_pin(npinx,npinx,nxyf,nz_th)
      REAL(8) tc_cpd_pin(npinx,npinx,nxyf,nz_th)
      REAL(8) ts_cpd_pin(npinx,npinx,nxyf,nz_th)
      REAL(8) tm_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) rm_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) cb_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) tc_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) ts_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) pw_mst(npinx,npinx,nxyf,nz_th)
      REAL(8) tmp0(npinx,npinx,nxyf,nz_th)
!
      INTEGER num_pin,num_chn,npoc
      INTEGER m,ji,ia,ia2,ixa,ixb,iya,iyb,ix,iy
      INTEGER npoc2
      INTEGER ixc,iyc
!
      LOGICAL,SAVE :: master_debug=.false.
      LOGICAL,SAVE :: init_mst=.true.
!      
!.....time & dt
!
      i_flag=i_opt
      TTIME=time
      DTTR=dt
!
!.....Initialize MASTER DLL
!
      IF(init_mst)then
         init_mst=.false.
         IF(myrank.eq.0)THEN
!DEC$IF defined (master_flag) 
            CALL MASTER_DLL(0, TTIME, DTTR,                             &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_INP)
            WRITE(*,*) 'i_masterC: Loading 3D Kinetics DLL (master.dll)'
!DEC$ENDIF 
         ENDIF
      ENDIF
!
!.....Read MASTER restart file
!
      IF(myrank.eq.0)THEN
!DEC$IF defined (master_flag) 
         IF(i_flag.eq.3)THEN
            RSTFN_REV="MAS_RST.STD.BIN"
            WRITE(* ,*)'t_masterC: Ready to read MAS_RST.STD.BIN'
            WRITE(unit_log,*)'t_masterC: Ready to read MAS_RST.STD.BIN'
            CALL MASTER_DLL(I_FLAG, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_REV)
            WRITE(* ,*)'t_masterC: Finish reading MAS_RST.STD.BIN'
            WRITE(unit_log,*)'t_masterC: Finish reading MAS_RST.STD.BIN'
            i_flag=2
         ELSEIF(i_flag.eq.5)then
            RSTFN_INP="MAS_RST.STD.BIN"
            WRITE(* ,*)'Ready to write MAS_RST.STD.BIN'
            WRITE(unit_log,*)'Ready to write MAS_RST.STD.BIN'               
            CALL MASTER_DLL(5, TTIME, DTTR,                        &
                               TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                               PWTH,                                       &
                               NXY_TH, NZ_TH,                              &
                               PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                               IOK, RSTFN_INP)
            WRITE(* ,*)'Finish writing MAS_RST.STD.BIN'
            WRITE(unit_log,*)'Finish writing MAS_RST.STD.BIN'
         ENDIF
!DEC$ENDIF 
      ENDIF
      IF(i_flag.eq.5)RETURN
!
!.....Debug mode (using old mapping in order to reduce memory usage)
!
      IF(master_debug)THEN
         CALL t_masterC_assembly_debug(i_opt)       
         RETURN
      ENDIF   
!    
!.....START t_master
!  
      IF(rv_model_master.eq.1)then
         num_pin=npinx
         num_chn=nchn
         npoc=num_pin/num_chn
!
!
!........Step 1. CUPID to Generic1 variables
!                cell%tl     => tm_cpd_chn
!                cell%rhol   => rm_cpd_chn
!                cell%cboron => cb_cpd_chn
!                t_fuel      => tc_cpd_pin/ts_cpd_pin
!
         tm_cpd_pin=0.0d0
         rm_cpd_pin=0.0d0
         cb_cpd_pin=0.0d0
         tc_cpd_pin=0.0d0
         ts_cpd_pin=0.0d0

         DO m=1,ncell_fuel_rod
            i=cupid_cell_hts2d(m)
            ji=jperm(i)
            ia=asm_ni(ji)                            
            IF(ia.eq.0) CYCLE 
            iz=asm_nz(ji)                     
            ixa=npoc*(chn_nx(ji)-1)+1                
            ixb=ixa+npoc-1                         
            iya=npoc*(chn_ny(ji)-1)+1                
            iyb=iya+npoc-1                         
            DO ix=ixa,ixb                          
               DO iy=iya,iyb
                  tm_cpd_pin(ix,iy,ia,iz)=cell%tl(i)     -273.15d0
                  rm_cpd_pin(ix,iy,ia,iz)=cell%rhol(i)
                  cb_cpd_pin(ix,iy,ia,iz)=cell%cboron(i)   
                  tc_cpd_pin(ix,iy,ia,iz)=t_fuel(m,    1)-273.15d0
                  ts_cpd_pin(ix,iy,ia,iz)=t_fuel(m,nr_2d)-273.15d0
               ENDDO                               
            ENDDO                                  
         ENDDO       
      
         IF(np.gt.1)then
            DO i=1,nz_th
               call allreducei_r(tm_cpd_pin(1,1,1,i),npinx*npinx*nxyf)
               call allreducei_r(rm_cpd_pin(1,1,1,i),npinx*npinx*nxyf)
               call allreducei_r(cb_cpd_pin(1,1,1,i),npinx*npinx*nxyf)
               call allreducei_r(tc_cpd_pin(1,1,1,i),npinx*npinx*nxyf)
               call allreducei_r(ts_cpd_pin(1,1,1,i),npinx*npinx*nxyf)
            ENDDO
         ENDIF
!
!........Step 2. Generic1 to Generic2 variables
!                tm_cpd_chn => tm_mst
!                rm_cpd_chn => rm_mst
!                cb_cpd_chn => cb_mst
!                tc_cpd_chn => tc_mst
!                ts_cpd_chn => ts_mst
!
         tm_mst=0.0d0
         rm_mst=0.0d0
         cb_mst=0.0d0
         tc_mst=0.0d0
         ts_mst=0.0d0
         
         DO iz=1,nz_th
            DO ia=1,num_asm
               DO ix=1,num_pin
                  IF(MOD(ix,npoc).eq.0)THEN
                     ixc=ix/npoc
                  ELSE
                     ixc=ix/npoc+1
                  ENDIF
                  DO iy=1,num_pin
                     IF(MOD(iy,npoc).eq.0)THEN
                        iyc=iy/npoc
                     ELSE
                        iyc=iy/npoc+1
                     ENDIF
                     tm_mst(ixc,iyc,ia,iz)=tm_mst(ixc,iyc,ia,iz)+tm_cpd_pin(ix,iy,ia,iz)
                     rm_mst(ixc,iyc,ia,iz)=rm_mst(ixc,iyc,ia,iz)+rm_cpd_pin(ix,iy,ia,iz)
                     cb_mst(ixc,iyc,ia,iz)=cb_mst(ixc,iyc,ia,iz)+cb_cpd_pin(ix,iy,ia,iz)
                     tc_mst(ixc,iyc,ia,iz)=tc_mst(ixc,iyc,ia,iz)+tc_cpd_pin(ix,iy,ia,iz)
                     ts_mst(ixc,iyc,ia,iz)=ts_mst(ixc,iyc,ia,iz)+ts_cpd_pin(ix,iy,ia,iz)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   
         !
         npoc2=npoc*npoc
         DO iz=1,nz_th
            DO ia=1,num_asm
               DO iy=1,num_pin
                  DO ix=1,num_pin
                     tm_mst(ix,iy,ia,iz)=tm_mst(ix,iy,ia,iz)/npoc2
                     rm_mst(ix,iy,ia,iz)=rm_mst(ix,iy,ia,iz)/npoc2
                     cb_mst(ix,iy,ia,iz)=cb_mst(ix,iy,ia,iz)/npoc2
                     tc_mst(ix,iy,ia,iz)=tc_mst(ix,iy,ia,iz)/npoc2
                     ts_mst(ix,iy,ia,iz)=ts_mst(ix,iy,ia,iz)/npoc2
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
!
!........Step 3. Generic2 to MASTER variables
!                tm_mst => TCOO
!                rm_mst => DCOO
!                cb_mst => BCOO
!                tc_mst => TFC
!                ts_mst => TFS
!
         DO iz=1,NZ_TH
            DO im=1,NXY_TH
               TFC(im,iz)=0.0d0
               TFS(im,iz)=0.0d0
               TCOO(im,iz)=0.0d0
               DCOO(im,iz)=0.0d0
               BCOO(im,iz)=0.0d0
            ENDDO
         ENDDO      

!........1x1 (nchn1d=1)
         IF(nchn.eq.1)then
            DO im=1,NXY_TH
               ia=mst_to_asmi(im)
               IF(ia.eq.0)CYCLE
               DO iz=1,NZ_TH
                  TCOO(im,iz)=tm_mst(1,1,ia,iz)
                  DCOO(im,iz)=rm_mst(1,1,ia,iz)
                  BCOO(im,iz)=cb_mst(1,1,ia,iz)
                  TFC(im,iz) =tc_mst(1,1,ia,iz)
                  TFS(im,iz) =ts_mst(1,1,ia,iz)
               ENDDO
            ENDDO
!........Other
!        2x2 (nchn=2)
!        4x4 (nchn=4)
!        8x8 (nchn=8)
         ELSE
            DO m=1,ncell_fuel_rod
               i=cupid_cell_hts2d(m)
               ji=jperm(i)
               ia=asm_ni(ji)
               ia2=asm_ni2(ji)
               iz=asm_nz(ji)
               ix=chn_nx(ji)
               iy=chn_ny(ji)
               TCOO(ia2,iz)=tm_mst(ix,iy,ia,iz)
               DCOO(ia2,iz)=rm_mst(ix,iy,ia,iz)
               BCOO(ia2,iz)=cb_mst(ix,iy,ia,iz)
               TFC( ia2,iz)=tc_mst(ix,iy,ia,iz)
               TFS( ia2,iz)=ts_mst(ix,iy,ia,iz)
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
         ENDIF
!      
      ELSEIF(rv_model_master.eq.0)then
         DO iz=1,NZ_TH
            DO im=1,NXY_TH
               k=master_to_assem1_cell(im,iz)
               IF(k.ne.0)then
                  temp_mix=(1.0d0-cell%alphag(k))*cell%tl(k) &
                                 +cell%alphag(k) *cell%tg(k)
                  TCOO(im,iz)=temp_mix-273.15d0
                  DCOO(im,iz)=cell%rhol(k) 
                  TFC(im,iz)= TCOO(im,iz)+300.d0
                  TFS(im,iz)= TCOO(im,iz)+  5.d0
               ENDIF
            ENDDO
         ENDDO
         IF(np.gt.1)THEN      
            DO i=1,NZ_TH
               CALL allreducei_r(TCOO(1,i),NXY_TH)
               CALL allreducei_r(DCOO(1,i),NXY_TH)
               CALL allreducei_r(TFC(1,i),NXY_TH)
               CALL allreducei_r(TFS(1,i),NXY_TH)
            ENDDO
         ENDIF
         
      ENDIF
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
         !   WRITE(* ,*)'Ready to write MAS_RST.STD.BIN'
         !   WRITE(unit_log,*)'Ready to write MAS_RST.STD.BIN'               
         !   CALL MASTER_DLL(5, TTIME, DTTR,                        &
         !                      TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
         !                      PWTH,                                       &
         !                      NXY_TH, NZ_TH,                              &
         !                      PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
         !                      IOK, RSTFN_INP)
         !   WRITE(* ,*)'Finish writing MAS_RST.STD.BIN'
         !   WRITE(unit_log,*)'Finish writing MAS_RST.STD.BIN'
         !   RETURN
         ELSE
            WRITE(* ,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC, but i_opt= !!!',i_opt
            WRITE(unit_log,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC !!!',i_opt              
            PAUSE
            STOP
         ENDIF
!DEC$ENDIF 
      ENDIF !myrank.eq.0
!
      tmp0=0.0d0
      IF(np.gt.1)then
         DO i=1,NZ_TH
            call allreducei_r(PIN3D_TH(1,1,1,i),NPINX*NPINX*NXYF)
         ENDDO
      ENDIF
!
!.....Step 4. MASTER to Generic3 variables
!             PIN3D_TH => pw_mst
!
      pw_mst=0.0d0
      DO iz=1,nz_th
         DO ia=1,num_asm
            DO ix=1,num_pin
               IF(MOD(ix,npoc).eq.0)THEN
                  ixc=ix/npoc
               ELSE
                  ixc=ix/npoc+1
               ENDIF
               DO iy=1,num_pin
                  IF(MOD(iy,npoc).eq.0)THEN
                     iyc=iy/npoc
                  ELSE
                     iyc=iy/npoc+1
                  ENDIF
                  pw_mst(ixc,iyc,ia,iz)=pw_mst(ixc,iyc,ia,iz)+PIN3D_TH(ix,iy,ia,iz)
               ENDDO
            ENDDO
         ENDDO
      ENDDO   

      p3d_cupid=0.0d0
      p3d_cupid_all=0.0d0
      qvol_mas=0.0d0
!
!.....Step 5. Generic3 to CUPID-RV variables
!             pw_mst => p3d_cupid (W)   => rv_hts
!                       qvol_mas  (W/m) => paraview
!      
      DO m=1,ncell_fuel_rod
         i=cupid_cell_hts2d(m)
         ji=jperm(i)
         ia=asm_ni(ji)
         iz=asm_nz(ji)
         ix=chn_nx(ji)
         iy=chn_ny(ji)
         p3d_cupid(m)=pw_mst(ix,iy,ia,iz)
         p3d_cupid_all=p3d_cupid_all+p3d_cupid(m)
         qvol_mas(i)=p3d_cupid(m)/cell_leng(i,3)
      ENDDO   
!
      p3d_cupid_all=0.0d0
      vol_tot=0.0d0
      IF(rv_model_master.eq.0)then
         DO i=1,ncell_fluid
            IF(nzone(i).eq.6)THEN
               p3d_cupid_all=p3d_cupid_all+qvol_mas(i)
               vol_tot=vol_tot+volp(i)
            ENDIF
         ENDDO
      ENDIF
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
      RETURN 
      END SUBROUTINE t_masterC_assembly
