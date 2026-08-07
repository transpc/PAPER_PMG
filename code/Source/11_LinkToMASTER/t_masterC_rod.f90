!
      SUBROUTINE t_masterC_rod(i_opt)
!
!     This routine calculates MASTER dll
!
      USE Vol_DATA      , ONLY: cell
      USE zmpi          , ONLY: jperm
      USE Zzone         , ONLY: ncell_fluid
      USE zcore         , ONLY: myrank,np
      USE Zconst1       , ONLY: vv_prob
      USE Zconst2       , ONLY: dt
      USE Ztimecon      , ONLY: time
      USE Zrv_ncell     , ONLY: ncell_fuel_rod,   &
                                cupid_cell_hts2d,p3d_cupid, &
                                qvol_mas, &
                                asm_ni,asm_ni2,asm_nz,chn_nx,chn_ny
      USE Zrv_hts_2d    , ONLY: t_fuel,nr_2d
      USE Zrv_model     , ONLY: rv_model
      USE Zcoord2       , ONLY: cell_leng
      USE Zporous       , ONLY: chn_type
!
!.....from linkmaster
!
      USE MASTER4    , ONLY: NCB,NXY_TH,NZ_TH,CBNAM,            &
                             NXYF,NPINX,                        &
                             TCOO,DCOO,BCOO,TFC,TFS,            &
                             ZCB,                               &
                             PIN3D_TH,                          &
                             TTIME,DTTR,PWTH,PPM,PPCT_MASTER,   &
                             power_master,                      &
                             ncb_mas,cbnam_mas,zcb_mas,ppm_mas, &
                             nchn
!DEC$IF defined (master_flag) 
      USE MASTER4    , ONLY: iok,P3D_TH,VOL_TH,                 & 
                             RSTFN_INP,RSTFN_REV
      USE Zio_unit      , ONLY: unit_log      
!DEC$ENDIF 
!
!.....Irregular axial mapping
      USE Zparam        , ONLY: ndim
      USE Zcoord1       , ONLY: xloc
      USE Zporous   , ONLY: nz_nk,hz_nk,nz_th0
!
      IMPLICIT none
      INCLUDE 'master_c.h'
      INCLUDE 'master_dll.h'
!
      INTEGER i,j
      INTEGER i_opt
!      
!.....CUPID-RV-MASETER
      INTEGER NZ
      REAL(8),ALLOCATABLE::PIN3D_TH0(:,:,:,:)
      REAL(8) p3d_cupid0
!      
!.....NEW mapping
      INTEGER m,ji,ia,ia2,iz,cx,cy
      INTEGER ni
      REAL(8) tm_cpd_chn(nchn,nchn,nxyf,nz_th0)
      REAL(8) rm_cpd_chn(nchn,nchn,nxyf,nz_th0)
      REAL(8) cb_cpd_chn(nchn,nchn,nxyf,nz_th0)
      REAL(8) tm_mst_pin(ncell_fuel_rod)
      REAL(8) rm_mst_pin(ncell_fuel_rod)
      REAL(8) cb_mst_pin(ncell_fuel_rod)

!DEC$IF defined (master_flag) 
      LOGICAL,SAVE :: init_opt1=.true.
!DEC$ENDIF 
      LOGICAL,SAVE :: init_mst=.true.
!
!.....Irregular axial mapping
      INTEGER iz0
      REAL(8) btm,top,h_nk
      
!
!======================================================================
!.....MASTER4 Calculation
!======================================================================
!
!.....time & dt
      TTIME=time-TTIME_0
      DTTR=dt
!
      IF(init_mst) then
         init_mst=.false.
         IF(myrank.eq.0)THEN
!DEC$IF defined (master_flag) 
!write(*,*)'i_opt be4 0 is ',i_opt
            CALL MASTER_DLL(0, TTIME, DTTR,                             &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_INP)
           !WRITE(unit_log,*) 'i_masterC: Loading 3D Kinetics DLL (master.dll)'
            WRITE(97,*) 'i_masterC: Loading 3D Kinetics DLL (master.dll)'
            WRITE(* ,*) 'i_masterC: Loading 3D Kinetics DLL (master.dll)'
!write(*,*)'i_opt aft 0 is ',i_opt
!DEC$ENDIF 
         ENDIF
      ENDIF
!
!.....Read MASTER restart file
!
      IF(myrank.eq.0)THEN
!DEC$IF defined (master_flag) 
         IF(i_opt.eq.3)THEN
            RSTFN_REV="MAS_RST.STD.BIN"
           !WRITE(unit_log,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            WRITE(* ,*)'i_masterC: Ready to read MAS_RST.STD.BIN'
            CALL MASTER_DLL(i_opt, TTIME, DTTR,                        &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_REV)
           !WRITE(unit_log,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            WRITE(97,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            WRITE(* ,*)'i_masterC: Finish reading MAS_RST.STD.BIN'
            i_opt=2
         ENDIF
!DEC$ENDIF 
      ENDIF
!
      IF(rv_model.eq.1)then
!
!........Step 1. CUPID to Generic1 variables
!                cell%tl     => tm_cpd_chn
!                cell%rhol   => rm_cpd_chn
!                cell%cboron => cb_cpd_chn
!                         
         tm_cpd_chn=0.0d0
         rm_cpd_chn=0.0d0
         cb_cpd_chn=0.0d0
         DO i=1,ncell_fluid
            ji=jperm(i)
            ia=asm_ni(ji)
            IF(ia.eq.0)CYCLE
            iz=asm_nz(ji)
            cx=chn_nx(ji)
            cy=chn_ny(ji)
            tm_cpd_chn(cx,cy,ia,iz)=cell%tl(i)-273.15d0
            rm_cpd_chn(cx,cy,ia,iz)=cell%rhol(i)
            cb_cpd_chn(cx,cy,ia,iz)=cell%cboron(i)
         ENDDO
         
         IF(np.gt.1)then
            DO i=1,nz_th0
               call allreducei_r(tm_cpd_chn(1,1,1,i),nchn*nchn*nxyf)
               call allreducei_r(rm_cpd_chn(1,1,1,i),nchn*nchn*nxyf)
               call allreducei_r(cb_cpd_chn(1,1,1,i),nchn*nchn*nxyf)
            ENDDO
         ENDIF
!
!........Step 2. Generic1 to Generic2 variables
!                tm_cpd_chn => tm_mst_pin
!                rm_cpd_chn => rm_mst_pin
!                cb_cpd_chn => cb_mst_pin
!                         
!........Step 3. Generic2 to MASTER variables
!                tm_mst_pin => TCOO
!                rm_mst_pin => DCOO
!                cb_mst_pin => BCOO
!
         TCOO=0.0d0
         DCOO=0.0d0
         BCOO=0.0d0
         TFC=0.0d0
         TFS=0.0d0
         tm_mst_pin=0.0d0
         rm_mst_pin=0.0d0
         cb_mst_pin=0.0d0
         DO m=1,ncell_fuel_rod
            i=cupid_cell_hts2d(m)
            IF(vv_prob.eq.'OPR1000_fullcore_modmesh02_rv_vessel')then
               IF(chn_type(i).eq.0)stop 'check cupid_cell_channel in averaging TFC/TFS'
            ENDIF
            ji=jperm(i)
            ia=asm_ni(ji)
            ia2=asm_ni2(ji)
            IF(ia2.eq.0) stop 'ia2 check' !CYCLE
            iz=asm_nz(ji)
            cx=chn_nx(ji)
            cy=chn_ny(ji)
            IF(iz.eq.0 .or. iz.gt.nz_th0)CYCLE
            tm_mst_pin(m)=0.25d0*(tm_cpd_chn(cx,cy,ia,iz)+tm_cpd_chn(cx+1,cy,ia,iz)+tm_cpd_chn(cx,cy+1,ia,iz)+tm_cpd_chn(cx+1,cy+1,ia,iz))
            rm_mst_pin(m)=0.25d0*(rm_cpd_chn(cx,cy,ia,iz)+rm_cpd_chn(cx+1,cy,ia,iz)+rm_cpd_chn(cx,cy+1,ia,iz)+rm_cpd_chn(cx+1,cy+1,ia,iz))
            cb_mst_pin(m)=0.25d0*(cb_cpd_chn(cx,cy,ia,iz)+cb_cpd_chn(cx+1,cy,ia,iz)+cb_cpd_chn(cx,cy+1,ia,iz)+cb_cpd_chn(cx+1,cy+1,ia,iz))
!
!...........1:1 axially mapping
!
if(0)then            
            TCOO(ia2,iz)=TCOO(ia2,iz)+tm_mst_pin(m)
            DCOO(ia2,iz)=DCOO(ia2,iz)+rm_mst_pin(m)
            BCOO(ia2,iz)=BCOO(ia2,iz)+cb_mst_pin(m)
            
            TFC(ia2,iz)=TFC(ia2,iz)+t_fuel(m,    1)-273.15d0
            TFS(ia2,iz)=TFS(ia2,iz)+t_fuel(m,nr_2d)-273.15d0
else
!
!...........axially IRREGULAR mapping
!
            btm=xloc(i,ndim)-0.5d0*cell_leng(i,ndim)
            top=xloc(i,ndim)+0.5d0*cell_leng(i,ndim)
            DO iz0=1,nz_nk
               h_nk=hz_nk(iz0+1)-hz_nk(iz0)
               IF(btm .ge. hz_nk(iz0) .and. btm .lt. hz_nk(iz0+1) .and. top .ge. hz_nk(iz0+1))THEN
                  TCOO(ia2,iz0)=TCOO(ia2,iz0) + tm_mst_pin(m)   * (hz_nk(iz0+1)-btm)/h_nk
                  DCOO(ia2,iz0)=DCOO(ia2,iz0) + rm_mst_pin(m)   * (hz_nk(iz0+1)-btm)/h_nk
                  BCOO(ia2,iz0)=BCOO(ia2,iz0) + cb_mst_pin(m)   * (hz_nk(iz0+1)-btm)/h_nk
                  TFC(ia2,iz0) =TFC(ia2,iz0)  +(t_fuel(m,    1)-273.15d0) * (hz_nk(iz0+1)-btm)/h_nk
                  TFS(ia2,iz0) =TFS(ia2,iz0)  +(t_fuel(m,nr_2d)-273.15d0) * (hz_nk(iz0+1)-btm)/h_nk
               ELSEIF(btm .lt. hz_nk(iz0) .and. top .ge. hz_nk(iz0+1))THEN
                  TCOO(ia2,iz0)=TCOO(ia2,iz0) + tm_mst_pin(m)
                  DCOO(ia2,iz0)=DCOO(ia2,iz0) + rm_mst_pin(m)
                  BCOO(ia2,iz0)=BCOO(ia2,iz0) + cb_mst_pin(m)
                  TFC(ia2,iz0) =TFC(ia2,iz0)  +(t_fuel(m,    1)-273.15d0)
                  TFS(ia2,iz0) =TFS(ia2,iz0)  +(t_fuel(m,nr_2d)-273.15d0)
               ELSEIF(btm .ge. hz_nk(iz0)  .and.  top .lt. hz_nk(iz0+1))THEN
                  TCOO(ia2,iz0)=TCOO(ia2,iz0) + tm_mst_pin(m)   * (top-btm)/h_nk
                  DCOO(ia2,iz0)=DCOO(ia2,iz0) + rm_mst_pin(m)   * (top-btm)/h_nk
                  BCOO(ia2,iz0)=BCOO(ia2,iz0) + cb_mst_pin(m)   * (top-btm)/h_nk
                  TFC(ia2,iz0) =TFC(ia2,iz0)  +(t_fuel(m,    1)-273.15d0) * (top-btm)/h_nk
                  TFS(ia2,iz0) =TFS(ia2,iz0)  +(t_fuel(m,nr_2d)-273.15d0) * (top-btm)/h_nk
               ELSEIF(btm .lt. hz_nk(iz0)  .and.  top .gt. hz_nk(iz0)  .and.  top .le. hz_nk(iz0+1))THEN
                  TCOO(ia2,iz0)=TCOO(ia2,iz0) + tm_mst_pin(m)   * (top-hz_nk(iz0))/h_nk
                  DCOO(ia2,iz0)=DCOO(ia2,iz0) + rm_mst_pin(m)   * (top-hz_nk(iz0))/h_nk
                  BCOO(ia2,iz0)=BCOO(ia2,iz0) + cb_mst_pin(m)   * (top-hz_nk(iz0))/h_nk
                  TFC(ia2,iz0) =TFC(ia2,iz0)  +(t_fuel(m,    1)-273.15d0) * (top-hz_nk(iz0))/h_nk
                  TFS(ia2,iz0) =TFS(ia2,iz0)  +(t_fuel(m,nr_2d)-273.15d0) * (top-hz_nk(iz0))/h_nk
               ENDIF
            ENDDO
endif            
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
!RV OFF
      ELSE
         DO j=1,NZ_TH
            DO i=1,NXY_TH
               TFC(i,j) =TCOO(i,j)+300.d0
               TFS(i,j) =TCOO(i,j)+  5.d0
            ENDDO
         ENDDO
      ENDIF
!
!.....Reflector region for MASTER Variables
!
      IF(rv_model.eq.1)then
         DO j=1,NZ_TH
            DO i=1,NXY_TH
               IF(TCOO(i,j).eq.0.0d0)then
                  TCOO(i,j)=296.21d0
                  DCOO(i,j)=735.18d0
                  TFC(i,j) =TCOO(i,j)+300.d0
                  TFS(i,j) =TCOO(i,j)+  5.d0
               ELSE
                  TCOO(i,j)=TCOO(i,j)/64.d0
                  DCOO(i,j)=DCOO(i,j)/64.d0
                  TFC(i,j)=TFC(i,j)/64.d0
                  TFS(i,j)=TFS(i,j)/64.d0
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO j=1,NZ_TH
            DO i=1,NXY_TH
               IF(TCOO(i,j).eq.0.0d0)then
                  TCOO(i,j)=296.21d0
                  DCOO(i,j)=735.18d0
                  TFC(i,j) =TCOO(i,j)+300.d0
                  TFS(i,j) =TCOO(i,j)+  5.d0
               ELSE
                  TCOO(i,j)=TCOO(i,j)/64.d0
                  DCOO(i,j)=DCOO(i,j)/64.d0
                  TFC(i,j) =TCOO(i,j)+300.d0
                  TFS(i,j) =TCOO(i,j)+  5.d0
               ENDIF
            ENDDO
         ENDDO
      ENDIF           
!
!.....Control rod
!
      ppm=ppm_mas
      ncb=ncb_mas !bug
      DO i=1,ncb_mas
         cbnam(i)=cbnam_mas(i)
         zcb(i)=zcb_mas(i)
      ENDDO
!
      PIN3D_TH=0.d0
      PWTH=power_master

      IF(myrank.eq.0)then
!DEC$IF defined (master_flag) 

         IF(i_opt.eq.1)then
            DO WHILE(1)
               CALL MASTER_DLL(i_opt, TTIME, DTTR,                        &
                               TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                               PWTH,                                       &
                               NXY_TH, NZ_TH,                              &
                               PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                               IOK, RSTFN_INP)
               IF(iok.eq.2)EXIT
            ENDDO
         ELSEIF(i_opt.eq.2)then
            IF(init_opt1)then
               init_opt1=.false.
               DO WHILE(1)
                  CALL MASTER_DLL(1, TTIME, DTTR,                             &
                                  TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                                  PWTH,                                       &
                                  NXY_TH, NZ_TH,                              &
                                  PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                                  IOK, RSTFN_INP)
                  IF(iok.eq.2)EXIT
               ENDDO         
            ENDIF
            CALL MASTER_DLL(i_opt, TTIME, DTTR,                         &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_INP)
         ELSEIF(i_opt.eq.5)then
            RSTFN_REV="MAS_RST.STD.BIN"
           !WRITE(unit_log,*)'Ready to write MAS_RST.STD.BIN'
            WRITE(97,*)'Ready to write MAS_RST.STD.BIN'
            WRITE(* ,*)'Ready to write MAS_RST.STD.BIN'
            CALL MASTER_DLL(i_opt, TTIME, DTTR,                         &
                            TFC, TFS, TCOO, DCOO, PPM, NCB, CBNAM, ZCB, &
                            PWTH,                                       &
                            NXY_TH, NZ_TH,                              &
                            PPCT_MASTER, P3D_TH, VOL_TH, PIN3D_TH,      &
                            IOK, RSTFN_REV)
           !WRITE(unit_log,*)'Finish writing MAS_RST.STD.BIN'
            WRITE(97,*)'Finish writing MAS_RST.STD.BIN'
            WRITE(* ,*)'Finish writing MAS_RST.STD.BIN'

         ELSE
            WRITE(unit_log,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC !!!',i_opt
            WRITE(97,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC, but i_opt=!!!',i_opt
            WRITE(* ,"(a,1i3)")'i_opt should be 1,2,5 in t_masterC, but i_opt=!!!',i_opt
            STOP
         ENDIF
!DEC$ENDIF 
      ENDIF !myrank.eq.0
      IF(i_opt.eq.5)RETURN
!     
      CALL barrier_mpi
      IF(np.gt.1) CALL broadcast_r1(PPCT_MASTER)
!     
      ALLOCATE(PIN3D_TH0(NPINX,NPINX,NXYF,NZ_TH))
      PIN3D_TH0=0.0d0
      IF(np.gt.1)then
         DO i=1,NZ_TH
            call allreducei_r(PIN3D_TH(1,1,1,i),NPINX*NPINX*NXYF)
         ENDDO
      ENDIF
      DEALLOCATE(PIN3D_TH0)      
!
!........Step 4. MASTER to CUPID-RV variables
!                PIN3D_TH => p3d_cupid

      p3d_cupid=0.0d0
      p3d_cupid0=0.0d0
      qvol_mas=0.0d0

      ! RV ON
      IF(rv_model.eq.1)then
if(0)then      
         DO m=1,ncell_fuel_rod
            i=cupid_cell_hts2d(m)
            ji=jperm(i)
            ni=asm_ni(ji)
            nz=asm_nz(ji)
            cx=chn_nx(ji)
            cy=chn_ny(ji)
            p3d_cupid(m)=PIN3D_TH(cx,cy,ni,nz)
            p3d_cupid0=p3d_cupid0+p3d_cupid(m)
            qvol_mas(i)=p3d_cupid(m)/cell_leng(i,3)
         ENDDO   
!
!NOT perfectly matched along axial direction (nz)
else
         DO m=1,ncell_fuel_rod
            i=cupid_cell_hts2d(m)
            ji=jperm(i)
            ni=asm_ni(ji)
            nz=asm_nz(ji)
            cx=chn_nx(ji)
            cy=chn_ny(ji)
            btm=xloc(i,ndim)-0.5d0*cell_leng(i,ndim)
            top=xloc(i,ndim)+0.5d0*cell_leng(i,ndim)
            DO iz=1,nz_nk
               h_nk=hz_nk(iz+1)-hz_nk(iz)
               IF(btm .ge. hz_nk(iz) .and. btm .lt. hz_nk(iz+1) .and. top .ge. hz_nk(iz+1))THEN
                  p3d_cupid(m)=p3d_cupid(m) + PIN3D_TH(cx,cy,ni,iz) * (hz_nk(iz+1)-btm)/h_nk
               ELSEIF(btm .lt. hz_nk(iz) .and. top .ge. hz_nk(iz+1))THEN
                  p3d_cupid(m)=p3d_cupid(m) + PIN3D_TH(cx,cy,ni,iz)
               ELSEIF(btm .ge. hz_nk(iz)  .and.  top .lt. hz_nk(iz+1))THEN
                  p3d_cupid(m)=p3d_cupid(m) + PIN3D_TH(cx,cy,ni,iz) * (top-btm)/h_nk
               ELSEIF(btm .lt. hz_nk(iz)  .and.  top .gt. hz_nk(iz)  .and.  top .le. hz_nk(iz+1))THEN
                  p3d_cupid(m)=p3d_cupid(m) + PIN3D_TH(cx,cy,ni,iz) * (top-hz_nk(iz))/h_nk
               ENDIF
            ENDDO

            p3d_cupid0=p3d_cupid0+p3d_cupid(m)
            qvol_mas(i)=p3d_cupid(m)/cell_leng(i,3)
         ENDDO   
endif         
      ! RV OFF
      ELSE
!         CALL p3d_to_qvol
!         DO i=1,ncell_fluid
!            IF(nzone(i).eq.2)then
!               p3d_cupid0=p3d_cupid0+qvol_liq(i)*volp(i)
!            ENDIF
!         ENDDO
      ENDIF

      IF(np.gt.1) CALL allreducei_r1(p3d_cupid0)
      IF(myrank.eq.0)then
        !write(unit_log,*)'total power in t_master_rod is',p3d_cupid0,time
         write(97,*)'total power in t_master_rod is',p3d_cupid0,time
         write(* ,*)'total power in t_master_rod is',p3d_cupid0,time
      ENDIF   
!
      CALL barrier_mpi
!
      RETURN
      END SUBROUTINE t_masterC_rod
!
!======================================================================
!======================================================================
!
      SUBROUTINE write_trod
!
!     This routine calculates MASTER dll
!
      USE Vol_DATA  , ONLY: cell
      USE Zmpi      , ONLY: jperm
      USE Zzone     , ONLY: ncell_fluid
      USE Zcore     , ONLY: myrank,np
      USE Zrv_ncell , ONLY: ncell_fuel_rod,nrod_fuel_rod,nz_fuel_rod
      USE Zrv_ncell , ONLY: num_ch,ncell_fluid_core,n_channel_fluid,nz_fluid,cupid_cell_channel
      USE Zrv_ncell , ONLY: ncell_fuel_rod_all
      USE Zrv_ncell , ONLY: p3d_cupid
      USE Zrv_hts_2d, ONLY: t_fuel,nrod_2d,nz0_2d,nr_2d
      USE MASTER4   , ONLY: nchn,nxyf
      USE Zrv_mpi   , ONLY: ncell_fuel_rod_p
      USE Zrv_ncell , ONLY: asm_ni,asm_nz,chn_nx,chn_ny
      USE Zporous   ,ONLY: chn_type
      USE Zporous   ,ONLY: nz_th0

      IMPLICIT NONE

      INTEGER i,j,k
      INTEGER ni,nz
      INTEGER ic
      INTEGER ji,iz,ia,cx,cy
      REAL(8) tcoo0(num_ch,nz0_2d)
      REAL(8) dcoo0(num_ch,nz0_2d)
      REAL(8) t_2d_input0(nrod_2d,nz0_2d,nr_2d)
      REAL(8) tchn(nchn,nchn,nxyf,nz_th0),rchn(nchn,nchn,nxyf,nz_th0)
!.....Local allocatable arrays
      REAL(8),DIMENSION(:),ALLOCATABLE :: p3d_cupid0
      REAL(8),DIMENSION(:,:),ALLOCATABLE :: t_fuel0

      RETURN ! Need to check in detail (near future)

      !channel fluid temperature/density - NOT necessary
      tcoo0=0.0d0
      dcoo0=0.0d0
      DO i=1,ncell_fluid_core
         ni=n_channel_fluid(i)
         nz=nz_fluid(i)
         ic=cupid_cell_channel(i)
         tcoo0(ni,nz)=cell%tl(ic)
         dcoo0(ni,nz)=cell%rhol(ic)
      ENDDO

      DO k=1,nz0_2d
         CALL allreducei_r(tcoo0(1,k),num_ch)
         CALL allreducei_r(dcoo0(1,k),num_ch)
      ENDDO

      !fuel rod temperature (t_2d_input format) - NOT necessary 
      t_2d_input0=0.0d0
      DO i=1,ncell_fuel_rod
         ni=nrod_fuel_rod(i) 
         nz=nz_fuel_rod(i)
         DO k=1,nr_2d
            t_2d_input0(ni,nz,k)=t_fuel(i,k)
         ENDDO
      ENDDO

      DO k=1,nr_2d
         CALL allreducei_r(t_2d_input0(1,1,k),nrod_2d*nz0_2d)
      ENDDO
  
if(0)then
      IF(myrank.eq.0)then
         OPEN(unit=1570,name='ht_str_2d_temp0.in')
         DO i=1,nrod_2d
            DO j=1,nz0_2d
               write(1570,45)(t_2d_input0(i,j,k),k=1,nr_2d)
            ENDDO
         ENDDO
         DO i=1,num_ch
            write(1570,45)(tcoo0(i,j),j=1,nz0_2d)
         ENDDO
         DO i=1,num_ch
            write(1570,45)(dcoo0(i,j),j=1,nz0_2d)
         ENDDO
         CLOSE(1570)
      ENDIF
45 format(50(f15.5,1x))
endif

      !fuel rod temperature - t_fuel format
      IF(np.gt.1) THEN
         ALLOCATE(t_fuel0(ncell_fuel_rod_all,nr_2d))
      ELSE
         ALLOCATE(t_fuel0(1,nr_2d))
      ENDIF
      CALL gatherv_r_2d(t_fuel,ncell_fuel_rod_p,t_fuel0,ncell_fuel_rod,ncell_fuel_rod_all,3)
!
      IF(myrank.eq.0)then
         OPEN(unit=1580,name='trod_core.in')
         DO i=1,ncell_fuel_rod_all
            write(1580,45) (t_fuel0(i,k),k=1,nr_2d)
         ENDDO
         CLOSE(1580)
      ENDIF
      DEALLOCATE(t_fuel0)


      !fuel rod power (W) - p3d_cupid
      IF(np.gt.1) THEN
         ALLOCATE(p3d_cupid0(ncell_fuel_rod_all))
      ELSE
         ALLOCATE(p3d_cupid0(1))
      ENDIF
      CALL gatherv_r(p3d_cupid,ncell_fuel_rod,p3d_cupid0,ncell_fuel_rod_all,3)
!
      IF(myrank.eq.0)then
         OPEN(unit=1581,name='qrod_core.in')
         DO i=1,ncell_fuel_rod_all
            write(1581,*) p3d_cupid0(i)
         ENDDO
         CLOSE(1581)
      ENDIF

      !channel-wise temperature/density
      tchn=0.0d0
      rchn=0.0d0
      DO i=1,ncell_fluid
         IF(chn_type(i).ne.0)then
            ji=jperm(i)
            iz=asm_nz(ji)
            ia=asm_ni(ji)
            cy=chn_ny(ji)
            cx=chn_nx(ji)
            IF(iz.lt.1 .and. iz.gt.nz_th0)STOP 'Error asm_nz for tchn'
            IF(ia.lt.1 .and. ia.gt.nxyf)  STOP 'Error asm_na for tchn'
            IF(cx.lt.1 .and. cx.gt.nchn)  STOP 'Error chn_nx for tchn'
            IF(cy.lt.1 .and. cy.gt.nchn)  STOP 'Error chn_ny for tchn'
            tchn(cx,cy,ia,iz)=cell%tl(i)
            rchn(cx,cy,ia,iz)=cell%rhol(i)
         ENDIF
      ENDDO
      IF(np.gt.1)then
         DO i=1,nz_th0
            call allreducei_r(tchn(1,1,1,i),nchn*nchn*nxyf)
            call allreducei_r(rchn(1,1,1,i),nchn*nchn*nxyf)
         ENDDO
      ENDIF

      IF(myrank.eq.0)then
         OPEN(1590,file='tchn_core.in',status='unknown')
         DO ia=1,nxyf
            DO cy=1,nchn
               DO cx=1,nchn
                  write(1590,46) ia,cx,cy,(tchn(cx,cy,ia,k),k=1,nz_th0),(rchn(cx,cy,ia,k),k=1,nz_th0)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
46  format(3(i4,1x),100(f15.8,1x))
!
      END SUBROUTINE write_trod
